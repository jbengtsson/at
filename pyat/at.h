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
#include <dlfcn.h>

#define  NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/ndarrayobject.h>

#define NUMPY_IMPORT_ARRAY_RETVAL NULL
#define NUMPY_IMPORT_ARRAY_TYPE void *

#if 0

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

#endif

#define LIMIT_AMPLITUDE 1

  
/* define the general signature of a pass function */
typedef union elem*
(*track_function)(const PyObject *element, union elem *elemptr,
		  double *r_in, int num_particles, struct parameters *param);

struct parameters {
  int    nturn;
  double RingLength;
  double T0;
};

static struct LibraryListElement {
  const char                *MethodName;
  void                      *LibraryHandle;
  track_function            FunctionHandle;
  struct LibraryListElement *Next;
} *LibraryList = NULL;

struct lat_type {
  bool
    *bxlost;
  unsigned int
    num_refpts;
  int
    num_turns,
    *ixnturn,
    *ixnelem;
  double
    *dxlostcoord;
  npy_uint32
    num_particles,
    np6,
    *refpts,
    keep_lattice,
    losses,
    omp_num_threads;
  npy_intp
    pdims[1],
    lxdims[2];
  PyObject
    *xnturn,
    *xnelem,
    *xlost,
    *xlostcoord;
  struct parameters
    param;
};
