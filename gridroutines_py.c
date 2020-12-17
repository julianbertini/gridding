/*  Python specific libraries - mandatory! */
#include </home/julianbertini/miniconda3/pkgs/python-3.7.7-hcff3b4d_5/include/python3.7m/Python.h>

/* other libraries */
#include <math.h>
#include <stdio.h>



static double *conv_PyObject_to_array(pyobj)
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

static int conv_array_to_PyObject(seq, pyobj)
/* 
    Converts a C array to a PyObject that holds a sequence (list) 

    Does NOT handle malloc() or free() calls. They must be handled by the caller.
    I.e., the memory must already be assigned for seq.

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
PyObject *kx;		/* Array of kx locations of samples. */
PyObject *ky;		/* Array of ky locations of samples. */
PyObject *dcf;		/* Output: Density compensation factors. */
PyObject *kerneltable;	/* 1D array of convolution kernel values, starting
				        at 0, and going to convwidth. */
int gridsize;		/* Number of points in kx and ky in grid. */
int nsamples;		/* Number of k-space samples, total. */
double convwidth;	/* Kernel width, in grid points.	*/

/* OTHER VARS */
int nkernelpts;		/* Number of points in kernel lookup-table */
int kcount1;		/* Counts through k-space sample locations */
int kcount2;		/* Counts through k-space sample locations */

double *dcfptr;			/* Aux. pointer, for looping.	*/
double *kxptr;			/* Aux. pointer. */
double *kyptr;			/* Aux. pointer. */

double kwidth;			/* Conv kernel width, in k-space units */
double dkx,dky,dk;		/* Delta in x, y and abs(k) for kernel calc.*/
int kernelind;			/* Index of kernel value, for lookup. */
double fracdk;			/* Fractional part of lookup index. */
double fracind;			/* Fractional lookup index. */
double kern;			/* Kernel value, avoid duplicate calculation.*/

kwidth = convwidth/(double)(gridsize);	/* Width of kernel, in k-space units. */

/* parse input argument  */
if (!PyArg_ParseTuple(args, "OOOOiid", &kx, &ky, &dcf, &kerneltable, &gridsize, &nsamples, &convwidth)) return NULL;

/* Here, need to convert Python lists n such into C variables */
// Going to write a helper for this. 

// My hope is that I can just modify the *dcf pointer directly and not have to
// return it to Python. Can easily test this first.

// Commenting everything else out to test initial array/PyObject conversions

// allocate memory
double  *result = conv_PyObject_to_array(kx);
ky = PySequence_Fast(ky, "argument must be iterable"); // new ref
conv_array_to_PyObject(result, ky);
if (result) free(result);


/* ========= Set all DCFs to 1. ========== */

/* 	DCF = 1/(sum( ksamps convolved with kernel ) 	*/

//dcfptr = dcf; 
//for (kcount1 = 0; kcount1 < nsamples; kcount1++)
//        *(dcfptr++) = 1.0;
// 
///* ========= Loop Through k-space Samples ========= */
//                
//dcfptr = dcf;
//kxptr = kx;
//kyptr = ky;
// 
//for (kcount1 = 0; kcount1 < nsamples; kcount1++)
//	{
//	/* printf("Current k-space location = (%f,%f) \n",*kxptr,*kyptr); */
//
//	for (kcount2 = kcount1+1; kcount2 < nsamples; kcount2++)
//		{
//		dkx = *kxptr-kx[kcount2];
//		dky = *kyptr-ky[kcount2];
//		dk = sqrt(dkx*dkx+dky*dky);	/* k-space separation*/
//		/*
//		printf("   Comparing with k=(%f,%f),  sep=%f \n",
//				kx[kcount2],ky[kcount2],dk);	
//		*/
//	
//		if (dk < kwidth)	/* sample affects this
//					   grid point		*/
//		    {
//			/* Find index in kernel lookup table */
//		    fracind = dk/kwidth*(double)(nkernelpts-1);
//		    kernelind = (int)fracind;
//		    fracdk = fracind-(double)kernelind;
//
//			/* Linearly interpolate in kernel lut */
//		    kern = kerneltable[(int)kernelind]*(1-fracdk)+
//		    		kerneltable[(int)kernelind+1]*fracdk;
//
//		    dcf[kcount2] += kern;
//		    *dcfptr += kern;
//		    }
//		}
//	dcfptr++;
//	kxptr++;
//	kyptr++;
//	}
//
//// this resets dcfptr to the beg. of array (pointed to by dcf)
//dcfptr = dcf; 
//for (kcount1 = 0; kcount1 < nsamples; kcount1++)
//    {
//    *dcfptr = 1/(*dcfptr);
//    dcfptr++;
//    }

// needed when writing a "void" routine
Py_INCREF(Py_None);
return Py_None;
}

static PyMethodDef gridmethods[] = {
/* */
    {"calcdcflut",  calcdcflut, METH_VARARGS,
     "Calculates the DCFs using a kernel look-up table"},
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


