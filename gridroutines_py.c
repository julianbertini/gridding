/*  Python specific libraries - mandatory! */
#include </home/julianbertini/miniconda3/pkgs/python-3.7.7-hcff3b4d_5/include/python3.7m/Python.h>

/* other libraries */
#include <math.h>
#include <stdio.h>



static double *
conv_PyObject_to_array(pyobj)
/* 
    Converts a PyObject that holds a sequence (list) into a C array

    Allocates memory for seq using malloc().
    CALLER RESPONSIBLE FOR FREEING MEMORY
*/

PyObject *pyobj; /* PyObject holding the Python list */

{
    PyObject *pyseq; /* PySequence object containing Python list */
    PyObject *item; /* Temp. variable holding current item in list */
    PyObject *fitem; /* Temp. variable holding current float item in list */
    int seqlen; /* Length of the pyseq list */
    double *seq; /* Pointer to C array */
    double *seqptr; /* auxiliary variable for looping*/

    pyseq = PySequence_Fast(pyobj, "argument must be iterable"); // new ref.
    if (!pyseq) return NULL;

    /* Prepare data as array of doubles */
    seqlen = PySequence_Fast_GET_SIZE(pyseq); 
    seq = malloc(seqlen*sizeof(double));
    seqptr = seq;

    for (int i = 0; i < seqlen; i++)
    {
        item = PySequence_Fast_GET_ITEM(pyseq, i);// borrowed ref.
        if (!item)
        {
            printf("%s", "problem with item\n");
            Py_DECREF(pyseq);
            free(seq);
            return NULL;
        }
        fitem = PyNumber_Float(item); // new ref.
        if (!fitem)
        {
            printf("%s", "problem with fitem\n");
            Py_DECREF(pyseq);
            free(seq);
            PyErr_SetString(PyExc_TypeError, "all items must be numbers");
            return NULL;
        }
        *seqptr = PyFloat_AS_DOUBLE(fitem); // set array value
        seqptr++; // advance pointer
        Py_DECREF(fitem);
    }
    Py_DECREF(pyseq);

    return seq;
}

static int 
conv_array_to_PyObject(seq, pyobj)
/* 
    Converts a C array to a PyObject that holds a sequence (list) 

    Does NOT handle malloc() or free() calls. They must be handled by the caller.

    pyobj is expected to already have PySequence support

    Returns 1 (true) on success, 0 (false) on failure.
*/

double *seq; /* Pointer to C array */
PyObject *pyobj; /* PyObject holding the Python list */

{
    PyObject *item; /* Holds item of PyObject sequence */
    int seqlen; /* Length of the pyseq list */
    double *seqptr = seq; /* auxiliary variable for looping */

    // check that pyobj provides sequence protocol
    if (!PySequence_Check(pyobj)) 
    {
        PyErr_SetString(PyExc_TypeError, "PyObject must provide sequence protocol.");
        return 0;
    }

    /* Prepare data as array of doubles */
    seqlen = PySequence_Size(pyobj); 

    for (int i = 0; i < seqlen; i++)
    {
        item = Py_BuildValue("d", *seqptr); // new ref
        if (PySequence_SetItem(pyobj, i, item) < 0) return 0;
        seqptr++; // advance pointer
    }

    return 1;
}

static PyObject *
calcdcflut(PyObject *self, PyObject *args)

/* Calculation of density compensation factors using 
			lookup-table convolution kernel */

/*------------------------------------------------------------------
	NOTES:

    I think: uses the (1 iter.) method from  Jackson et al. - Julian
    See paper for reference. 

    TODO: need to extend this to use Pipe et al. iterative method
    test this first, though. I think there is a bug below. 

	Assume gridding uses the following formula, which 
	describes the contribution of each data sample 
	to the value at each grid point:

		grid-point value += data value / dcf * kernel(dk)

	where:
		data value is the complex sampled data value.
		dcf is the density compensation factor for the sample point.
		kernel is the convolution kernel function.
		dk is the k-space separation between sample point and
			grid point.

	"grid units"  are integers from 0 to gridsize-1, where
	the grid represents a k-space of -.5 to .5.

	The convolution kernel is passed as a series of values 
	that correspond to "kernel indices" 0 to nkernelpoints-1.
  ------------------------------------------------------------------ */

{
    /* INPUT/OUTPUT PARAMETERS */
    PyObject *py_kx;		/* Array of kx locations of samples. */
    PyObject *py_ky;		/* Array of ky locations of samples. */
    PyObject *py_dcf;		/* Output: Density compensation factors. */
    PyObject *py_kerneltable;	/* 1D array of convolution kernel values, starting
                            at 0, and going to convwidth. */
    int gridsize;		/* Number of points in kx and ky in grid. */
    int nsamples;		/* Number of k-space samples, total. */
    int nkernelpts;		/* Number of points in kernel lookup-table */
    double convwidth;	/* Kernel width, in grid points. */

    /* OTHER VARS */
    int kcount1;		/* Counts through k-space sample locations */
    int kcount2;		/* Counts through k-space sample locations */

    double *dcf;
    double *kx;
    double *ky;
    double *kerneltable;

    double *dcfptr;			/* Aux. pointer, for looping. */
    double *kxptr;			/* Aux. pointer. */
    double *kyptr;			/* Aux. pointer. */

    double kwidth;			/* Conv kernel width, in k-space units */
    double dkx,dky,dk;		/* Delta in x, y and abs(k) for kernel calc.*/
    int kernelind;			/* Index of kernel value, for lookup. */
    double fracdk;			/* Fractional part of lookup index. */
    double fracind;			/* Fractional lookup index. */
    double kern;			/* Kernel value, avoid duplicate calculation.*/

    /* parse input argument  */
    if (!PyArg_ParseTuple(args, "OOOOiiid", &py_kx, &py_ky, &py_dcf, &py_kerneltable, 
    &gridsize, &nsamples, &nkernelpts, &convwidth)) return NULL;

    kwidth = convwidth/(double)(gridsize);	/* Width of kernel, in k-space units. */

    /* convert PyObject lists to C pointers */
    kx = conv_PyObject_to_array(py_kx); // allocates memory
    ky = conv_PyObject_to_array(py_ky); // allocates memory
    dcf = conv_PyObject_to_array(py_dcf); // allocates memory
    kerneltable = conv_PyObject_to_array(py_kerneltable); // allocates memory

    /* ========= Set all DCFs to 1. ========== */

    /* 	DCF = 1/(sum( ksamps convolved with kernel )  */

    dcfptr = dcf; 
    for (kcount1 = 0; kcount1 < nsamples; kcount1++)
        *(dcfptr++) = 1.0;
     
    /* ========= Loop Through k-space Samples ========= */
                    
    dcfptr = dcf;
    kxptr = kx;
    kyptr = ky;
     
    for (kcount1 = 0; kcount1 < nsamples; kcount1++)
        {
        /* printf("Current k-space location = (%f,%f) \n",*kxptr,*kyptr); */

        for (kcount2 = kcount1+1; kcount2 < nsamples; kcount2++)
            {
            dkx = *kxptr-kx[kcount2];
            dky = *kyptr-ky[kcount2];
            dk = sqrt(dkx*dkx+dky*dky);	/* k-space separation*/
            /*
            printf("   Comparing with k=(%f,%f),  sep=%f \n",
                    kx[kcount2],ky[kcount2],dk);	
            */
        
            if (dk < kwidth)	/* sample affects this
                           grid point */
                {
                /* Find index in kernel lookup table */
                fracind = dk/kwidth*(double)(nkernelpts-1);
                kernelind = (int)fracind;
                fracdk = fracind-(double)kernelind;

                /* Linearly interpolate in kernel lut */
                kern = kerneltable[(int)kernelind]*(1-fracdk)+
                        kerneltable[(int)kernelind+1]*fracdk;

                dcf[kcount2] += kern;
                *dcfptr += kern;
                }
            }
        dcfptr++;
        kxptr++;
        kyptr++;
        }

    // this resets dcfptr to the beg. of array (pointed to by dcf)
    dcfptr = dcf; 
    for (kcount1 = 0; kcount1 < nsamples; kcount1++)
    {
        *dcfptr = 1/(*dcfptr);
        dcfptr++;
    }

    if (!conv_array_to_PyObject(dcf, py_dcf)) return NULL;
    
    free(dcf);
    free(kx);
    free(ky);
    free(kerneltable);

    // needed when writing a "void" Python C routine
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *
gridlut(PyObject *self, PyObject *args)

/* Gridding function that uses a lookup table for a circularly
	symmetric convolution kernel, with linear interpolation.  
	See Notes below */

/* INPUT/OUTPUT */


/*------------------------------------------------------------------
	NOTES:

	This uses the following formula, which describes the contribution
	of each data sample to the value at each grid point:

		grid-point value += data value / dcf * kernel(dk)

	where:
		data value is the complex sampled data value.
		dcf is the density compensation factor for the sample point.
		kernel is the convolution kernel function.
		dk is the k-space separation between sample point and
			grid point.

	"grid units"  are integers from 0 to gridsize-1, where
	the grid represents a k-space of -.5 to .5.

	The convolution kernel is passed as a series of values 
	that correspond to "kernel indices" 0 to nkernelpoints-1.
  ------------------------------------------------------------------ */

{
     /* INPUT/OUTPUT PARAMETERS */
    PyObject *py_kx;		/* Array of kx locations of samples. */
    PyObject *py_ky;		/* Array of ky locations of samples. */
    PyObject *py_s_real;		/* Sampled data, real part. */
    PyObject *py_s_imag; 	/* Sampled data, real part. */
    PyObject *py_dcf;		/* Density compensation factors. */
    PyObject *py_sg_real;	/* OUTPUT array, real parts of data. */
    PyObject *py_sg_imag;	/* OUTPUT array, imag parts of data. */	 
    PyObject *py_kerneltable;	/* 1D array of convolution kernel values, starting
                               at 0, and going to convwidth. */

    int nsamples;		/* Number of k-space samples, total. */
    int gridsize;		/* Number of points in kx and ky in grid. */
    int nkernelpts;		/* Number of points in kernel lookup-table */
    double convwidth;	/* Kernel width, in grid points.	*/

    /* OTHER VARS */
    double *kx;		/* Array of kx locations of samples. */
    double *ky;		/* Array of ky locations of samples. */
    double *s_real;		/* Sampled data, real part. */
    double *s_imag; 	/* Sampled data, real part. */
    double *dcf;		/* Density compensation factors. */
    double *sg_real;	/* OUTPUT array, real parts of data. */
    double *sg_imag;	/* OUTPUT array, imag parts of data. */	 
    double *kerneltable;	/* 1D array of convolution kernel values, starting
                               at 0, and going to convwidth. */

    int kcount;		/* Counts through k-space sample locations */
    int gcount1, gcount2;	/* Counters for loops */
    //int col;		/* Grid Columns, for faster lookup. */

    double kwidth;			/* Conv kernel width, in k-space units */
    double dkx,dky,dk;		/* Delta in x, y and abs(k) for kernel calc.*/
    int ixmin,ixmax,iymin,iymax;	/* min/max indices that current k may affect*/
    int kernelind;			/* Index of kernel value, for lookup. */
    double fracdk;			/* Fractional part of lookup index. */
    double fracind;			/* Fractional lookup index. */
    double kern;			/* Kernel value, avoid duplicate calculation.*/

    double *dcfptr;			/* Aux. pointer, for looping. */
    double *kxptr;			/* Aux. pointer. */
    double *kyptr;			/* Aux. pointer. */
    double *srptr;          /* Aux. pointer. */
    double *siptr;          /* Aux. pointer. */
    double *sgrptr;			/* Aux. pointer, for loop. */
    double *sgiptr;			/* Aux. pointer, for loop. */

    int gridsizesq;			/* Square of gridsize */
    int gridcenter;			/* Index in output of kx,ky=0 points. */
    int gptr_cinc;			/* Increment for grid pointer. */

    /* parse input argument  */
    if (!PyArg_ParseTuple(args, "OOOOOOOOiiid", &py_kx, &py_ky, &py_s_real, &py_s_imag, 
    &py_dcf, &py_sg_real, &py_sg_imag, &py_kerneltable, &nsamples, &gridsize, &nkernelpts, &convwidth)) return NULL;

    gridcenter = gridsize/2;	/* Index of center of grid. */
    kwidth = convwidth/(double)(gridsize);	/* Width of kernel, in k-space units. */

    /* convert PyObjects to C pointers 
       all allocate memory; must free once done. */
    kx = conv_PyObject_to_array(py_kx);
    ky = conv_PyObject_to_array(py_ky);
    s_real = conv_PyObject_to_array(py_s_real);
    s_imag = conv_PyObject_to_array(py_s_imag);
    dcf = conv_PyObject_to_array(py_dcf);
    sg_real = conv_PyObject_to_array(py_sg_real);
    sg_imag = conv_PyObject_to_array(py_sg_imag);
    kerneltable = conv_PyObject_to_array(py_kerneltable);


    /* ========= Zero Output Points ========== */

    dcfptr = dcf;
    kxptr = kx;
    kyptr = ky;
    srptr = s_real;
    siptr = s_imag;
    sgrptr = sg_real;
    sgiptr = sg_imag;
    gridsizesq = gridsize*gridsize;
     
    // setting each element in array to zero using the pointers
    for (gcount1 = 0; gcount1 < gridsizesq; gcount1++)
        {
        *sgrptr++ = 0;
        *sgiptr++ = 0;
        }

     
    /* ========= Loop Through k-space Samples ========= */
                    
    for (kcount = 0; kcount < nsamples; kcount++)
        {

        /* ----- Find limit indices of grid points that current
             sample may affect (and check they are within grid) ----- */

        /* This is where the (kx, ky) positions come into play.

           They tell us the bounds of the conv. kernel region, which
           we then use to calculate the data values at the locations
           within the kernel region below. */

        ixmin = (int) ((*kxptr-kwidth)*gridsize + gridcenter);
        if (ixmin < 0) ixmin=0;
        ixmax = (int) ((*kxptr+kwidth)*gridsize + gridcenter)+1;
        if (ixmax >= gridsize) ixmax=gridsize-1;
        iymin = (int) ((*kyptr-kwidth)*gridsize + gridcenter);
        if (iymin < 0) iymin=0;
        iymax = (int) ((*kyptr+kwidth)*gridsize + gridcenter)+1;
        if (iymax >= gridsize) iymax=gridsize-1;

            
        /* Increment for grid pointer at end of column to top of next col.
           
           Recall how in C we save a 2D array as a long 1D array with indexing
           defined by A[row*NUM_COLS + col] 

           The y-dim goes along the cols; the x-dim goes along the rows */

        gptr_cinc = gridsize-(iymax-iymin)-1;	/* 1 bc at least 1 sgrptr++ */
            
        /* start the pointer at the min positions where the kernel has an effect */
        sgrptr = sg_real + (ixmin*gridsize+iymin);
        sgiptr = sg_imag + (ixmin*gridsize+iymin);
                            
        for (gcount1 = ixmin; gcount1 <= ixmax; gcount1++)
            {
            dkx = (double)(gcount1-gridcenter) / (double)gridsize  - *kxptr;
            for (gcount2 = iymin; gcount2 <= iymax; gcount2++)
                {
                dky = (double)(gcount2-gridcenter) / 
                (double)gridsize - *kyptr;

                dk = sqrt(dkx*dkx+dky*dky);	/* k-space separation*/

                if (dk < kwidth)	/* sample affects this
                               grid point */
                    {
                    /* Find index in kernel lookup table */
                    fracind = dk/kwidth*(double)(nkernelpts-1);
                    kernelind = (int)fracind;
                    fracdk = fracind-(double)kernelind;

                    /* Linearly interpolate in kernel lut */
                    kern = kerneltable[(int)kernelind]*(1-fracdk)+
                            kerneltable[(int)kernelind+1]*fracdk;

                    /* This is the conv. step. since the s_real/s_imag contain
                       discrete points, this is like conv. the kernel with a 
                       bunch of delta functions. So we can just mult. the
                       kernel value at the given disp. by the data point value.
                       No need to add a bunch of values together, since data
                       point is like a delta func. */
                    *sgrptr += kern * *srptr * *dcfptr;
                    *sgiptr += kern * *siptr * *dcfptr;
                    }
                sgrptr++;
                sgiptr++;
                }
            sgrptr+= gptr_cinc;
            sgiptr+= gptr_cinc;
            }
        kxptr++;		/* Advance kx pointer */
        kyptr++;   	/* Advance ky pointer */
        dcfptr++;		/* Advance dcf pointer */		
        srptr++;	/* Advance real-sample pointer */
        siptr++;	/* Advance imag-sample pointer */
        }

    if (!conv_array_to_PyObject(sg_real, py_sg_real)) return NULL;
    if (!conv_array_to_PyObject(sg_imag, py_sg_imag)) return NULL;

    free(dcf);
    free(kx);
    free(ky);
    free(s_real);
    free(s_imag);
    free(sg_real);
    free(sg_imag);
    free(kerneltable);

    // needed when writing a "void" Python C routine
    Py_INCREF(Py_None);
    return Py_None;
}

/*** PYTHON MODULE CREATION METHODS  ***/

static PyMethodDef gridmethods[] = {
/* */
    {"calcdcflut",  calcdcflut, METH_VARARGS,
     "Calculates the DCFs using a kernel look-up table"},
    {"gridlut",  gridlut, METH_VARARGS,
     "Grids data using DCFs and a kernel look-up table"},
    {NULL, NULL, 0, NULL} /* Sentinel */
};

PyDoc_STRVAR(gridroutines_doc, "Various gridding routines.");
static struct PyModuleDef gridmodule = {
    PyModuleDef_HEAD_INIT,
    "grid",   /* name of module */
    gridroutines_doc, /* module documentation, may be NULL */
    -1,       /* size of per-interpreter state of the module,
                 or -1 if the module keeps state in global variables. */
    gridmethods
};

PyMODINIT_FUNC
PyInit_grid(void)
{
    return PyModule_Create(&gridmodule);
}


