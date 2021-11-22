/*
 * The file to be included by 'new-style' integrators that support both
 * Matlab and Python.
 */
#ifndef ATELEM_C
#define ATELEM_C

#include "atcommon.h"

/*----------------------------------------------------*/
/*            For the integrator code                 */
/*----------------------------------------------------*/

#define atIsFinite isfinite
#define atIsNaN isnan
#define atGetNaN() (NAN)
#define atGetInf() (INFINITY)
#define atMalloc malloc
#define atCalloc calloc
#define atFree free


/*----------------------------------------------------*/
/*            For the Python interface                */
/*----------------------------------------------------*/

#if defined(PYAT)

typedef PyObject atElem;
#define check_error() if (PyErr_Occurred()) return NULL

static int array_imported = 0;

static NUMPY_IMPORT_ARRAY_TYPE init_numpy(void)
{
  import_array();
  return NUMPY_IMPORT_ARRAY_RETVAL;
}

static long atGetLong(const PyObject *element, const char *name)
{
  const PyObject *attr = PyObject_GetAttrString((PyObject *)element, name);
  if (!attr) return 0L;
  Py_DECREF(attr);
  return PyLong_AsLong((PyObject *)attr);
}

static double atGetDouble(const PyObject *element, const char *name)
{
  const PyObject *attr = PyObject_GetAttrString((PyObject *)element, name);
  if (!attr) return 0.0;
  Py_DECREF(attr);
  return PyFloat_AsDouble((PyObject *)attr);
}

static long
atGetOptionalLong(const PyObject *element, const char *name, long default_value)
{
  long l = atGetLong(element, name);
  if (PyErr_Occurred()) {
    PyErr_Clear();
    l = default_value;
  }
  return l;
}

static double
atGetOptionalDouble(const PyObject *element, const char *name,
		    double default_value)
{
  double d = atGetDouble(element, name);
  if (PyErr_Occurred()) {
    PyErr_Clear();
    d = default_value;
  }
  return d;
}

static double*
atGetArrayData(PyArrayObject *array, char *name, int atype, int *msz,
	       int *nsz)
{
  char errmessage[60];
  int ndims;
  npy_intp *dims;
  if (!array_imported) {
    init_numpy();
    array_imported = 1;
  }
  Py_DECREF(array);
  if (!PyArray_Check(array)) {
    snprintf(errmessage, 60, "The attribute %s is not an array.", name);
    PyErr_SetString(PyExc_RuntimeError, errmessage);
    return NULL;
  }
  if (PyArray_TYPE(array) != atype) {
    snprintf(errmessage, 60, "The attribute %s is not a double array.",
	     name);
    PyErr_SetString(PyExc_RuntimeError, errmessage);
    return NULL;
  }
  if ((PyArray_FLAGS(array) & NPY_ARRAY_FARRAY_RO) != NPY_ARRAY_FARRAY_RO) {
    snprintf(errmessage, 60, "The attribute %s is not Fortran-aligned.",
	     name);
    PyErr_SetString(PyExc_RuntimeError, errmessage);
    return NULL;
  }
  ndims = PyArray_NDIM(array);
  dims = PyArray_SHAPE(array);
  *nsz = (ndims >= 2) ? dims[1] : 0;
  *msz = (ndims >= 1) ? dims[0] : 0;
  return (double *) PyArray_DATA(array);
}

static double*
atGetDoubleArraySz(const PyObject *element, char *name, int *msz, int *nsz)
{
  PyArrayObject *array =
    (PyArrayObject *) PyObject_GetAttrString((PyObject *)element, name);
  if (array == NULL) {
    return NULL;
  }
  return (double*) atGetArrayData(array, name, NPY_DOUBLE, msz, nsz);
}

static double* atGetDoubleArray(const PyObject *element, char *name)
{
  int msz, nsz;
  return atGetDoubleArraySz(element, name, &msz, &nsz);
}

static double*
atGetOptionalDoubleArraySz(const PyObject *element, char *name, int *msz,
			   int *nsz)
{
  PyArrayObject *array =
    (PyArrayObject *) PyObject_GetAttrString((PyObject *)element, name);
  if (array == NULL) {
    PyErr_Clear();
    return NULL;
  }
  return (double*) atGetArrayData(array, name, NPY_DOUBLE, msz, nsz);
}

static double* atGetOptionalDoubleArray(const PyObject *element, char *name)
{
  int msz, nsz;
  return atGetOptionalDoubleArraySz(element, name, &msz, &nsz);
}

#endif /* defined(PYAT) */

#if defined(PYAT) || defined(MATLAB_MEX_FILE)
#include "attypes.h"

#ifdef __cplusplus
#define C_LINK extern "C"
#else
#define C_LINK
#endif

C_LINK ExportMode struct elem_type*
trackFunction(const atElem *ElemData, struct elem_type *Elem, double *r_in,
	      int num_particles, struct parameters *Param);

#endif /* defined(PYAT) || defined(MATLAB_MEX_FILE) */

#endif /*ATELEM_C*/
