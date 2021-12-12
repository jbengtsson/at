#define  PY_SSIZE_T_CLEAN
#include <Python.h>
#include <stdarg.h>
#ifdef _OPENMP
  #include <string.h>
  #include <omp.h>
#endif /*_OPENMP*/
#include <stdbool.h> 
#include <math.h>
#include <float.h>

#define  NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/ndarrayobject.h>

#define NUMPY_IMPORT_ARRAY_RETVAL NULL
#define NUMPY_IMPORT_ARRAY_TYPE void *

#define ATPY_PASS "trackFunction"

#if defined(PCWIN) || defined(PCWIN64) || defined(_WIN32)
  #include <windows.h>
  #define LIBRARYHANDLETYPE        HINSTANCE
  #define FREELIBFCN(libfilename)  FreeLibrary((libfilename))
  #define LOADLIBFCN(libfilename)  LoadLibrary((libfilename))
  #define GETTRACKFCN(libfilename) GetProcAddress((libfilename), ATPY_PASS)
  #define SEPARATOR "\\"
  #define OBJECTEXT ".pyd"
#else
  #include <dlfcn.h>
  #define LIBRARYHANDLETYPE void *
  #define FREELIBFCN(libfilename)  dlclose(libfilename)
  #define LOADLIBFCN(libfilename)  dlopen((libfilename), RTLD_LAZY)
  #define GETTRACKFCN(libfilename) dlsym((libfilename), ATPY_PASS)
  #define SEPARATOR "/"
  #define OBJECTEXT ".so"
#endif

#define LIMIT_AMPLITUDE 1

  
/* define the general signature of a pass function */
typedef union elem*
(*track_function)(const PyObject *element, union elem *elemptr,
		  double *r_in, int num_particles, struct parameters *param);

static npy_uint32     num_elements         = 0;
static union elem     **elemdata_list      = NULL;
static PyObject       **element_list       = NULL;
static track_function *integrator_list     = NULL;
static PyObject       **pyintegrator_list  = NULL;
static PyObject       **kwargs_list        = NULL;
static char           integrator_path[300];

/* For AT defined in atpass.c. */
static struct LibraryListElement {
  const char                *MethodName;
  LIBRARYHANDLETYPE         LibraryHandle;
  track_function            FunctionHandle;
  PyObject                  *PyFunctionHandle;
  struct LibraryListElement *Next;
} *LibraryList = NULL;

struct lat_type {
  bool
    *bxlost;
  int
    *ixnturn,
    *ixnelem;
  double
    *dxlostcoord;
  npy_uint32
    num_particles;
  npy_intp
    pdims[1],
    lxdims[2];
  PyObject
    *xnturn,
    *xnelem,
    *xlost,
    *xlostcoord;
};
