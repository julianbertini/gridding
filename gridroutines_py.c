/*  Python specific libraries - mandatory! */
#include </home/julianbertini/miniconda3/pkgs/python-3.7.7-hcff3b4d_5/include/python3.7m/Python.h>

/* other libraries */
#include <math.h>
#include <stdio.h>


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
double *kx;		/* Array of kx locations of samples. */
double *ky;		/* Array of ky locations of samples. */
int nsamples;		/* Number of k-space samples, total. */
double *dcf;		/* Output: Density compensation factors. */
int gridsize;		/* Number of points in kx and ky in grid. */
double convwidth;	/* Kernel width, in grid points.	*/
double *kerneltable;	/* 1D array of convolution kernel values, starting
				at 0, and going to convwidth. */

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
if (!PyArg_ParseTuple(args, "s", &command))
    return NULL;

/* Here, need to convert Python lists n such into C variables */
// Going to write a helper for this. 

// Need helper that will take a Python list object and convert it to C array
// Need helper that will take a C array and convert it to Python object 
// So one for going each way. I will need to do this item by item.

// Then, will need to return the Python list object, which will store the 
// dcf information



/* ========= Set all DCFs to 1. ========== */

/* 	DCF = 1/(sum( ksamps convolved with kernel ) 	*/

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
					   grid point		*/
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

// needed when writing a "void" routine
Py_INCREF(Py_None);
return Py_None;
}


static PyObject *
spam_system(PyObject *self, PyObject *args)
{
    const char *command;
    int sts;

    if (!PyArg_ParseTuple(args, "s", &command))
        return NULL;
    sts = system(command);
    return PyLong_FromLong(sts);
}

static PyMethodDef SpamMethods[] = {
/* */
    {"system",  spam_system, METH_VARARGS,
     "Execute a shell command."},
    {NULL, NULL, 0, NULL}        /* Sentinel */
};

PyDoc_STRVAR(spam_doc, "Nom nom nom nom.");
static struct PyModuleDef spammodule = {
    PyModuleDef_HEAD_INIT,
    "spam",   /* name of module */
    spam_doc, /* module documentation, may be NULL */
    -1,       /* size of per-interpreter state of the module,
                 or -1 if the module keeps state in global variables. */
    SpamMethods
};

PyMODINIT_FUNC
PyInit_spam(void)
{
    return PyModule_Create(&spammodule);
}
