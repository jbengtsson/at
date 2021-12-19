/* Python interface for Tracy-2; compatible with Python 3 only.
   It provides a module 'atpass' containing one method 'atpass'.              */

#include <string>
#include <sstream>
#include <iostream>

#include "at_types.h"
#include "at.h"


static npy_uint32
  num_elements         = 0;
static union elem
  **elemdata_list      = NULL;
static PyObject
  **element_list       = NULL,
  **kwargs_list        = NULL;
static track_function
  *integrator_list     = NULL;
static std::string
  integrator_path,
  integrator_lib;


//------------------------------------------------------------------------------

static PyObject *print_error(int elem_number, PyObject *rout)
{
  printf("Error in tracking element %d.\n", elem_number);
  Py_XDECREF(rout);
  return NULL;
}

static PyObject *set_error(PyObject *errtype, const char *fmt, ...)
{
  char    buffer[132];
  va_list ap;

  va_start(ap, fmt);
  vsprintf(buffer, fmt, ap);
  PyErr_SetString(errtype, buffer);
  va_end(ap);
  return NULL;
}

/* Query Python for the full extension given to shared objects.
   This is useful for Python 3, where the extension may not be trivial.
   If none is defined, return NULL.                                           */
static PyObject* get_ext_suffix(void)
{
  PyObject *sysconfig_module, *get_config_var_fn, *ext_suffix;

  sysconfig_module = PyImport_ImportModule("distutils.sysconfig");
  if (sysconfig_module == NULL) return NULL;
  get_config_var_fn =
    PyObject_GetAttrString(sysconfig_module, "get_config_var");
  Py_DECREF(sysconfig_module);
  if (get_config_var_fn == NULL) return NULL;
  ext_suffix = PyObject_CallFunction(get_config_var_fn, "s", "EXT_SUFFIX");
  Py_DECREF(get_config_var_fn);
  return ext_suffix;
}

static PyObject *isopenmp(PyObject *self)
{
#ifdef _OPENMP
  Py_RETURN_TRUE;
#else
  Py_RETURN_FALSE;
#endif /*_OPENMP)*/
}

//------------------------------------------------------------------------------


#include "integrator_helper.cc"


bool get_at_lat(PyObject *&args, PyArrayObject *&rin, PyObject *&element)
{
  if (!PyArg_ParseTuple(args, "OO!", &element, &PyArray_Type, &rin)) {
    return false;
  }
  if (PyArray_NDIM(rin) != 2) {
    set_error(PyExc_ValueError, "rin is not of two dimensions");
    return false;
  }
  if (PyArray_DIM(rin, 0) != PS_DIM) {
    set_error(PyExc_ValueError, "rin is not 6D");
    return false;
  }
  if (PyArray_TYPE(rin) != NPY_DOUBLE) {
    set_error(PyExc_ValueError, "rin is not a double array");
    return false;
  }
  if ((PyArray_FLAGS(rin) & NPY_ARRAY_FARRAY_RO) != NPY_ARRAY_FARRAY_RO) {
    set_error(PyExc_ValueError, "rin is not Fortran-aligned");
    return false;
  }
  return true;
}

static PyObject *at_elempass(PyObject *self, PyObject *args)
{
  PyObject                  *element, *PyPassMethod;
  PyArrayObject             *rin;
  npy_uint32                num_particles;
  track_function            integrator;
  double                    *drin;
  struct parameters         param;
  struct LibraryListElement *LibraryListPtr;
  const char *  str = NULL;
  if (!get_at_lat(args, rin, element)) return NULL;

  num_particles = (PyArray_SIZE(rin)/PS_DIM);
  drin = static_cast<double*>(PyArray_DATA(rin));

  param.RingLength = 0e0;
  param.T0         = 0e0;
  param.nturn      = 0;

  PyPassMethod = PyObject_GetAttrString(element, "PassMethod");
  if (!PyPassMethod) return NULL;

  std::string  fn_name = PyUnicode_AsUTF8(PyPassMethod);
  integrator = lookup_function(fn_name);
  if(!integrator){
    return NULL;
  }
  std::cerr << "Calling integrator method " << fn_name
	    << " using function " << integrator
	    << std::endl;
  std::cerr.flush();

  union elem *elem_data =
    integrator(element, NULL, drin, num_particles, &param);
  std::cerr << " done" << std::endl;
  std::cerr.flush();

  if (!elem_data) return NULL;
  free(elem_data);
  Py_RETURN_NONE;
}


bool get_at_lat(PyObject *args, PyObject *kwargs,
		PyObject *&lattice, PyArrayObject *&rin, int &num_turns,
		PyArrayObject *&refs, npy_uint32 &keep_lattice,
		npy_uint32 &omp_num_threads, npy_uint32 &losses)
{
  static char
    *kwlist[] =
    {(char*)"line", (char*)"rin", (char*)"nturns", (char*)"refpts",
     (char*)"reuse", (char*)"omp_num_threads", (char*)"losses",
     (char*)NULL};

  if (!PyArg_ParseTupleAndKeywords
      (args, kwargs, "O!O!i|O!III", kwlist, &PyList_Type, &lattice,
       &PyArray_Type, &rin, &num_turns, &PyArray_Type, &refs, &keep_lattice,
       &omp_num_threads, &losses)) {
    return false;
  }
  if (PyArray_DIM(rin, 0) != PS_DIM) {
    PyErr_SetString(PyExc_ValueError, "Numpy array is not 6D");
    return false;
  }
  if (PyArray_TYPE(rin) != NPY_DOUBLE) {
    PyErr_SetString(PyExc_ValueError, "rin is not a double array");
    return false;
  }
  if ((PyArray_FLAGS(rin) & NPY_ARRAY_FARRAY_RO) != NPY_ARRAY_FARRAY_RO) {
    PyErr_SetString(PyExc_ValueError, "rin is not Fortran-aligned");
    return false;
  }
  return true;
}

bool get_lat(PyObject *lattice, double &lattice_length, PyObject *rout)
{
  npy_uint32     elem_index;
  PyObject       **element;
  track_function *integrator, t_integ;
  PyObject       *el;
  PyObject       *PyPassMethod;
  double         length;

  /* Release the stored elements */
  for (elem_index = 0; elem_index < num_elements; elem_index++) {
    free(elemdata_list[elem_index]);
    /* Release the stored elements, may be NULL if */
    Py_XDECREF(element_list[elem_index]);
  }                          /* a previous call was interrupted by an error */
  num_elements = PyList_Size(lattice);

  /* Pointer to Element structures used by the tracking function */
  free(elemdata_list);
  elemdata_list =
    (union elem **)calloc(num_elements, sizeof(union elem *));

  /* Pointer to Element list, make sure all pointers are initially NULL */
  free(element_list);
  element_list = (PyObject **)calloc(num_elements, sizeof(PyObject *));

  /* pointer to the list of C integrators */
  integrator_list =
    (track_function*)realloc(integrator_list,
			      num_elements*sizeof(track_function));

  free(kwargs_list);
  kwargs_list = (PyObject **)calloc(num_elements, sizeof(PyObject *));

  lattice_length = 0e0;
  element = element_list;
  integrator = integrator_list;
  for (elem_index = 0; elem_index < num_elements; elem_index++) {
    el           = PyList_GET_ITEM(lattice, elem_index);
    PyPassMethod = PyObject_GetAttrString(el, "PassMethod");

    if (!PyPassMethod) {
      /* No PassMethod */
      print_error(elem_index, rout);
      return false;
    }

    std::string fh_name = PyUnicode_AsUTF8(PyPassMethod);
    t_integ = lookup_function(fh_name);
    if (!t_integ) {
      /* No trackFunction for the given PassMethod */
      print_error(elem_index, rout);
      return false;
    }

    length = PyFloat_AsDouble(PyObject_GetAttrString(el, "Length"));
    if (PyErr_Occurred())
      PyErr_Clear();
    else
      lattice_length += length;

    //std::cerr << "Element "<< index << "using  integrator method " << fh_name
    //	      << " using function " << (void *) t_integ	 << std::endl;
    //std::cerr.flush();

    *integrator++ = t_integ;
    *element++ = el;

    /* Keep a reference to each element in case of reuse */
    Py_INCREF(el);
    Py_DECREF(PyPassMethod);
  }
  return true;
}

void init_losses(lat_type &lat)
{
  unsigned int  i;
  static double r0[PS_DIM];

  lat.pdims[0] = lat.num_particles;
  lat.lxdims[0] = PS_DIM;
  lat.lxdims[1] = lat.num_particles;
  lat.xnturn = PyArray_EMPTY(1, lat.pdims, NPY_UINT32, 1);
  lat.xnelem = PyArray_EMPTY(1, lat.pdims, NPY_UINT32, 1);
  lat.xlost = PyArray_EMPTY(1, lat.pdims, NPY_BOOL, 1);
  lat.xlostcoord = PyArray_EMPTY(2, lat.lxdims, NPY_DOUBLE, 1);
  lat.ixnturn = static_cast<int*>(PyArray_DATA((PyArrayObject*)lat.xnturn));
  lat.ixnelem = static_cast<int*>(PyArray_DATA((PyArrayObject*)lat.xnelem));
  lat.bxlost = static_cast<bool*>(PyArray_DATA((PyArrayObject*)lat.xlost));
  lat.dxlostcoord =
    static_cast<double*>(PyArray_DATA((PyArrayObject*)lat.xlostcoord));

  for(i = 0; i < lat.num_particles; i++) {
    lat.bxlost[i]  = 0;
    lat.ixnturn[i] = 0;
    lat.ixnelem[i] = 0;
    memcpy(lat.dxlostcoord+PS_DIM*i, r0, PS_DIM*sizeof(double));
  }
}

static void
check_if_lost(lat_type &lat, double *drin, const int num_elem)
{
  unsigned int n, c;
  double       *r6;

  for (c = 0; c < lat.num_particles; c++) {/* Loop over particles */
    if (!lat.bxlost[c]) {  /* No change if already marked */
      r6 = drin+c*PS_DIM;
      for (n = 0; n < PS_DIM; n++) {
	if (!isfinite(r6[n]) || ((fabs(r6[n])>LIMIT_AMPLITUDE)&&n<5)) {
	  lat.bxlost[c] = true;
	  lat.ixnturn[c] = lat.param.nturn;
	  lat.ixnelem[c] = num_elem;
	  memcpy(lat.dxlostcoord+PS_DIM*c, r6, PS_DIM*sizeof(double));
	  r6[0] = NAN;
	  r6[1] = 0;
	  r6[2] = 0;
	  r6[3] = 0;
	  r6[4] = 0;
	  r6[5] = 0;
	  break;
	}
      }
    }
  }
}

static void set_lost(double *drin, npy_uint32 np)
{
  unsigned int n, c;
  double       *r6;

  for (c = 0; c < np; c++) {/* Loop over particles */
    r6 = drin+c*PS_DIM;
    if (isfinite(r6[0])) {  /* No change if already marked */
      for (n = 0; n < PS_DIM; n++) {
	if (!isfinite(r6[n]) || ((fabs(r6[n]) > LIMIT_AMPLITUDE) && n < 5)) {
	  r6[0] = NAN;
	  r6[1] = 0;
	  r6[2] = 0;
	  r6[3] = 0;
	  r6[4] = 0;
	  r6[5] = 0;
	  break;
	}
      }
    }
  }
}

bool track(lat_type &lat, double *drin, double *&drout, PyObject *rout)
{
  unsigned int
    nextrefindex;
  union elem
    **elemdata   = elemdata_list;
  npy_uint32
    elem_index,
    nextref;
  PyObject
    **element    = element_list;
  track_function
    *integrator  = integrator_list;
  PyObject
    **kwargs     = kwargs_list;

  nextrefindex = 0;
  nextref = (nextrefindex < lat.num_refpts)?
    lat.refpts[nextrefindex++] : INT_MAX;
  for (elem_index = 0; elem_index < num_elements; elem_index++) {
    if (elem_index == nextref) {
      memcpy(drout, drin, lat.np6*sizeof(double));
      /*  shift the location to write to in the output array */
      drout += lat.np6;
      nextref = (nextrefindex < lat.num_refpts)?
	lat.refpts[nextrefindex++] : INT_MAX;
    }

    /* the actual integrator call */
    {
      track_function t_integrator = * integrator;
      /*
      std::cerr << "Calling integrator " << (void *) t_integrator
		<< " Elem " << *elemdata
		<< " elemdata " << *element
		<< " prams " << &lat.param
		<< std::endl;
      */
    *elemdata =
      t_integrator(*element, *elemdata, drin, lat.num_particles, &lat.param);
    // std::cerr << "done. returned " << *elemdata << std::endl;
    }
    /* trackFunction failed */
    if (!*elemdata) {
      print_error(elem_index, rout);
      return false;
    }

    if (lat.losses)
      check_if_lost(lat, drin, elem_index);
    else
      set_lost(drin, lat.num_particles);

    element++;
    integrator++;
    elemdata++;
    kwargs++;
  }
  /* the last element in the ring */
  if (num_elements == nextref) {
    memcpy(drout, drin, lat.np6*sizeof(double));
    /*  shift the location to write to in the output array */
    drout += lat.np6;
  }
  return true;
}

PyObject* get_losses(lat_type &lat, PyObject *rout)
{
  PyObject
    *tout = NULL,
    *dict =  NULL;

  tout = PyTuple_New(2);
  if(!tout){
    return NULL;
  }

  dict = PyDict_New();
  if(!dict){
    return NULL;
  }
  PyDict_SetItemString(dict, (char *)"islost", (PyObject *)lat.xlost);
  PyDict_SetItemString(dict, (char *)"turn",   (PyObject *)lat.xnturn);
  PyDict_SetItemString(dict, (char *)"elem",   (PyObject *)lat.xnelem);
  PyDict_SetItemString(dict, (char *)"coord",  (PyObject *)lat.xlostcoord);
  PyTuple_SetItem(tout,  0,  rout);
  PyTuple_SetItem(tout,  1,  dict);
  // Py_DECREF(lat.xlost);
  // Py_DECREF(lat.xnturn);
  // Py_DECREF(lat.xnelem);
  // Py_DECREF(lat.xlostcoord);

  return tout;
}

/* Parse the arguments to atpass, set things up, and execute.
   Arguments:
     line:   sequence of elements
     rin:    numpy 6-vector of initial conditions
     nturns: int number of turns to simulate
     refpts: numpy uint32 array denoting elements at which to return state
     reuse:  whether to reuse the cached state of the ring                    */

static PyObject* at_atpass(PyObject *self, PyObject *args, PyObject *kwargs) {
  static bool
    valid = false;
  static double
    lattice_length = 0e0;

  int
    turn;
  double
    *drin,
    *drout;
  npy_intp
    outdims[4];
  PyObject
    *lattice;
  PyArrayObject
    *rin,
    *refs = NULL;
  PyObject
    *rout;
  struct LibraryListElement
    *LibraryListPtr;

  struct lat_type lat;

#ifdef _OPENMP
  int maxthreads;
#endif /*_OPENMP*/

  lat.keep_lattice    = 0;
  lat.losses          = 0;
  lat.omp_num_threads = 0;

  if (!get_at_lat
      (args, kwargs, lattice, rin, lat.num_turns, refs, lat.keep_lattice,
       lat.omp_num_threads, lat.losses))
    return NULL;

  lat.num_particles = (PyArray_SIZE(rin)/PS_DIM);
  lat.np6           = lat.num_particles*PS_DIM;
  drin              = static_cast<double*>(PyArray_DATA(rin));

  if (refs) {
    if (PyArray_TYPE(refs) != NPY_UINT32) {
      PyErr_SetString(PyExc_ValueError, "refpts is not a uint32 array");
      return NULL;
    }
    lat.refpts = static_cast<unsigned int*>(PyArray_DATA(refs));
    lat.num_refpts = PyArray_SIZE(refs);
  } else {
    lat.refpts = NULL;
    lat.num_refpts = 0;
  }

  outdims[0] = PS_DIM;
  outdims[1] = lat.num_particles;
  outdims[2] = lat.num_refpts;
  outdims[3] = lat.num_turns;

  rout  = PyArray_EMPTY(4, outdims, NPY_DOUBLE, 1);
  drout = static_cast<double*>(PyArray_DATA((PyArrayObject*)rout));

  if (lat.losses) init_losses(lat);

#ifdef _OPENMP
  if ((lat.omp_num_threads > 0) && (num_particles > OMP_PARTICLE_THRESHOLD)) {
    unsigned int nthreads = omp_get_num_procs();
    maxthreads = omp_get_max_threads();
    if (lat.omp_num_threads < nthreads) nthreads = lat.omp_num_threads;
    if (num_particles < nthreads) nthreads = num_particles;
    omp_set_num_threads(nthreads);
  }
#endif /*_OPENMP*/

  if (!(lat.keep_lattice && valid)) {
    if (!get_lat(lattice, lattice_length, rout)) return NULL;
    valid = false;
  }

  lat.param.RingLength = lattice_length;
  lat.param.T0         = lattice_length/C0;

  for (turn = 0; turn < lat.num_turns; turn++) {
    lat.param.nturn = turn;
    if (!track(lat, drin, drout, rout)) return NULL;
  }

  /* Tracking successful: the lattice can be reused */
  valid = true;

#ifdef _OPENMP
  if ((lat.omp_num_threads > 0) && (num_particles > OMP_PARTICLE_THRESHOLD)) {
    omp_set_num_threads(maxthreads);
  }
#endif /*_OPENMP*/

  if (lat.losses) {
    return get_losses(lat, rout);
  } else {
    return rout;
  }
}

//------------------------------------------------------------------------------

static PyMethodDef AtMethods[] =
  {
   {"elempass", (PyCFunction)at_elempass, METH_VARARGS,
    PyDoc_STR
    ("elempass(element, rin)\n\n"
     "Track input particles rin through a single element.\n\n"
     "element: AT element\n"
     "rin:     6 x n_particles Fortran-ordered numpy array.\n"
     "         On return, rin contains the final coordinates of the particles\n"
     )},

   {"atpass", (PyCFunction)at_atpass, METH_VARARGS | METH_KEYWORDS,
    PyDoc_STR
    ("rout = atpass(line, rin, n_turns, refpts=[], reuse=False"
     ", omp_num_threads=0)\n\n"
     "Track input particles rin along line for nturns turns.\n"
     "Record 6D phase space at elements corresponding to refpts for each turn."
     "\n\n"
     "line:    list of elements\n"
     "rin:     6 x n_particles Fortran-ordered numpy array.\n"
     "         On return, rin contains the final coordinates of the particles\n"
     "n_turns: number of turns to be tracked\n"
     "refpts:  numpy array of indices of elements where output is desired\n"
     "         0 means entrance of the first element\n"
     "         len(line) means end of the last element\n"
     "reuse:   if True, use previously cached description of the lattice.\n\n"
     "rout:    6 x n_particles x n_refpts x n_turns Fortran-ordered numpy array"
     "\n"
     "         of particle coordinates\n"
     )},

   {"isopenmp", (PyCFunction)isopenmp, METH_NOARGS,
    PyDoc_STR("isopenmp()\n\n"
	      "Return whether OpenMP is active.\n"
	      )},
   {NULL, NULL, 0, NULL}        /* Sentinel */
  };


/* Use Python calls to establish the location of the at integrators package.  */
static PyObject* get_integrators(void) {
  PyObject *at_module, *os_module, *fileobj, *dirname_function, *dirobj;

  at_module = PyImport_ImportModule(home_dir.c_str());
  if (at_module == NULL) return NULL;
  fileobj = PyObject_GetAttrString(at_module, "__file__");
  Py_DECREF(at_module);
  if (fileobj == NULL) return NULL;
  os_module = PyImport_ImportModule("os.path");
  if (os_module == NULL) return NULL;
  dirname_function = PyObject_GetAttrString(os_module, "dirname");
  Py_DECREF(os_module);
  if (dirname_function == NULL) return NULL;
  dirobj = PyObject_CallFunctionObjArgs(dirname_function, fileobj, NULL);
  Py_DECREF(fileobj);
  Py_DECREF(dirname_function);
  return dirobj;
}

PyMODINIT_FUNC PyInit_atpass(void)
{
  PyObject   *integ_path_obj, *ext_suffix_obj;
  const char *ext_suffix, *integ_path;

  static struct PyModuleDef moduledef =
    {
     PyModuleDef_HEAD_INIT,
     "at",                                                /* m_name */
     PyDoc_STR("Clone of atpass in Accelerator Toolbox"), /* m_doc */
     -1,                                                  /* m_size */
     AtMethods,                                           /* m_methods */
     NULL,                                                /* m_reload */
     NULL,                                                /* m_traverse */
     NULL,                                                /* m_clear */
     NULL,                                                /* m_free */
    };
  PyObject *m = PyModule_Create(&moduledef);

  if (m == NULL) return NULL;
  import_array();

  integ_path_obj = get_integrators();
  if (integ_path_obj == NULL) return NULL;
  ext_suffix_obj = get_ext_suffix();
  if (ext_suffix_obj == NULL) return NULL;
  ext_suffix =
    (ext_suffix_obj == Py_None)? ".so" : PyUnicode_AsUTF8(ext_suffix_obj);
  integ_path = PyUnicode_AsUTF8(integ_path_obj);
  integrator_path = integ_path;
  integrator_lib = ext_suffix;
  Py_DECREF(integ_path_obj);
  Py_DECREF(ext_suffix_obj);

  return m;
}

//------------------------------------------------------------------------------
