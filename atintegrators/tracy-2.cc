#include <math.h>

#include <armadillo>

#include "at_common.h"
#include "at_types.h"

#include "tracy-2.h"

#include "at_lalib.cc"


//------------------------------------------------------------------------------

/*----------------------------------------------------*/
/*            For the Python interface                */
/*----------------------------------------------------*/

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

static long atGetOptionalLong(const PyObject *element, const char *name,
			      long default_value)
{
  long l = atGetLong(element, name);
  if (PyErr_Occurred()) {
    PyErr_Clear();
    l = default_value;
  }
  return l;
}

static double atGetOptionalDouble(const PyObject *element, const char *name,
				  double default_value)
{
  double d = atGetDouble(element, name);
  if (PyErr_Occurred()) {
    PyErr_Clear();
    d = default_value;
  }
  return d;
}

static double* atGetArrayData(PyArrayObject *array, char *name, int atype,
			      int *msz, int *nsz)
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

static double* atGetDoubleArraySz(const PyObject *element, char *name, int *msz,
				  int *nsz)
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

static double* atGetOptionalDoubleArraySz(const PyObject *element, char *name,
					  int *msz, int *nsz)
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


//------------------------------------------------------------------------------

// Interface functions for: Armidillo vectors, STL vectors, and C arrays.

inline std::vector<double> vectostl(const arma::vec &vec)
{ return {vec(x_), vec(px_), vec(y_), vec(py_), vec(delta_), vec(ct_), 1e0}; }

inline arma::vec stltovec(const std::vector<double> &a)
{ return {a[x_], a[px_], a[y_], a[py_], a[delta_], a[ct_], 1e0}; }

arma::vec arrtovec(const double a[])
{ return {a[x_], a[px_], a[y_], a[py_], a[delta_], a[ct_], 1e0}; }

void vectoarr(const arma::vec &vec, double a[])
{
  for (int k = 0; k < PS_DIM; k++)
    a[k] = vec(k);
}

inline std::vector<double> arrtostl(const double a[])
{
  return {a[x_], a[px_], a[y_], a[py_], a[delta_], a[ct_], 1e0};
}

void stltoarr(const std::vector<double> &vec, double a[])
{
  for (int k = 0; k < PS_DIM; k++)
    a[k] = vec[k];
}

arma::mat arrtomat(const double a[])
{
  arma::mat A(PS_DIM, PS_DIM);

  for (int j = 0; j < PS_DIM; j++)
    for (int k = 0; k < PS_DIM; k++)
      A(j, k) = a[j*PS_DIM+k];
  return A;
}

void mattoarr(const arma::mat &A, double a[])
{
  int j, k;

  for (j = 0; j < PS_DIM; j++)
    for (k = 0; k < PS_DIM; k++)
      a[j*PS_DIM+k] = A(j, k);
}

//------------------------------------------------------------------------------
/*
 * consider removing elem von argument list
 */
struct elem_type* init_elem(const PyObject *ElemData, struct elem_type *Elem,
			    const bool len, const bool aper)
{
  double Length, *R1, *R2, *T1, *T2, *EApertures, *RApertures;

  Elem = (struct elem_type*)calloc(1, sizeof(struct elem_type));
  if(!Elem){
    return NULL;
  }
  // std::cerr << __FUNCTION__ << " allocated Elem " << Elem << std::endl;
  if (len) {
    Length = atGetDouble(ElemData, "Length");
    check_error();
  } else
    Length = 0e0;
  R1 = atGetOptionalDoubleArray(ElemData, (char*)"R1");
  check_error();
  R2 = atGetOptionalDoubleArray(ElemData, (char*)"R2");
  check_error();
  T1 = atGetOptionalDoubleArray(ElemData, (char*)"T1");
  check_error();
  T2 = atGetOptionalDoubleArray(ElemData, (char*)"T2");
  check_error();

  Elem->Length     = Length;
  Elem->R1         = R1;
  Elem->R2         = R2;
  Elem->T1         = T1;
  Elem->T2         = T2;

  if (aper) {
    EApertures = atGetOptionalDoubleArray(ElemData, (char*)"EApertures");
    check_error();
    RApertures = atGetOptionalDoubleArray(ElemData, (char*)"RApertures");
    check_error();

    Elem->EApertures = EApertures;
    Elem->RApertures = RApertures;
  }

  return Elem;
}

struct elem_type* init_id(const PyObject *ElemData, struct elem_type *Elem)
{
  Elem = init_elem(ElemData, Elem, false, true);
  if (Elem) {
    Elem->id_ptr = (struct elem_id*)malloc(sizeof(struct elem_id));
    return Elem;
  } else
    return NULL;
}

struct elem_type* init_ap(const PyObject *ElemData, struct elem_type *Elem)
{
  double *limits;

  Elem = (struct elem_type*)calloc(1, sizeof(struct elem_type));
  Elem->Length = 0e0;
  if (Elem) {
    Elem->ap_ptr = (struct elem_ap*)calloc(1, sizeof(struct elem_ap));

    limits = atGetDoubleArray(ElemData, (char*)"Limits");
    check_error();

    Elem->ap_ptr->limits = limits;
    return Elem;
  } else
    return NULL;
}

struct elem_type* init_drift(const PyObject *ElemData, struct elem_type *Elem)
{
  Elem = init_elem(ElemData, Elem, true, true);
  if (Elem) {
    Elem->drift_ptr = (struct elem_drift*)calloc(1, sizeof(struct elem_drift));
    return Elem;
  } else
    return NULL;
}

struct elem_type*
init_mpole(const PyObject *ElemData, struct elem_type *Elem, const bool bend,
	   const bool cbend, const bool incl_E0, const bool incl_E2)
{
  int
    MaxOrder, NumIntSteps,
    FringeQuadEntrance, FringeQuadExit;
  double
    *PolynomA, *PolynomB, *fringeIntM0, *fringeIntP0, *KickAngle;


  double
    h1 = 0.0,            ///< todo: good default?
    h2  =0.0,            ///< todo: good default?
    RefDZ = 0.0,         ///< todo: good default?
    Energy = 1.9,        ///< todo: good default?
    ByError = 0.0,       ///< todo: good default?
    X0ref = 0.0,         ///< todo: good default?
    FringeInt1 = 0.0,    ///< todo: good default?
    FringeInt2 = 0.0,    ///< todo: good default?
    FullGap = 0.0,       ///< todo: good default?
    EntranceAngle = 0.0, ///< todo: good default?
    ExitAngle = 0.0,     ///< todo: good default?
    BendingAngle = 0.0,   ///< todo: good default?
    FringeBendEntrance = 0.0,  ///< todo: good default?
    FringeBendExit = 0.0 ///< todo: good default?
    ;
  elem_mpole
    *mpole;

  Elem = init_elem(ElemData, Elem, true, !cbend);
  if (Elem) {
    Elem->mpole_ptr = (struct elem_mpole*)calloc(1, sizeof(struct elem_mpole));
    mpole           = Elem->mpole_ptr;

    if (incl_E0) {
      Energy = atGetDouble(ElemData,                (char*)"Energy");
      check_error();
    }

    PolynomA = atGetDoubleArray(ElemData,           (char*)"PolynomA");
    check_error();
    PolynomB = atGetDoubleArray(ElemData,           (char*)"PolynomB");
    check_error();
    MaxOrder = atGetLong(ElemData,                  (char*)"MaxOrder");
    check_error();

    NumIntSteps = atGetLong(ElemData,               (char*)"NumIntSteps");
    check_error();

    if (bend) {
      BendingAngle = atGetDouble(ElemData,          (char*)"BendingAngle");
      check_error();
      EntranceAngle = atGetDouble(ElemData,         (char*)"EntranceAngle");
      check_error();
      ExitAngle = atGetDouble(ElemData,             (char*)"ExitAngle");
      check_error();

      FringeBendEntrance =
	atGetOptionalLong(ElemData,                 (char*)"FringeBendEntrance",
			  1);
      check_error();
      FringeBendExit = atGetOptionalLong(ElemData,  (char*)"FringeBendExit", 1);
      check_error();

      FullGap = atGetOptionalDouble(ElemData,       (char*)"FullGap", 0);
      check_error();
      FringeInt1 = atGetOptionalDouble(ElemData,    (char*)"FringeInt1", 0);
      check_error();
      FringeInt2 = atGetOptionalDouble(ElemData,    (char*)"FringeInt2", 0);
      check_error();
    }

    if (cbend) {
      X0ref = atGetOptionalDouble(ElemData,         (char*)"X0ref", 0);
      check_error();
      ByError = atGetOptionalDouble(ElemData,       (char*)"ByError", 0);
      check_error();
      RefDZ = atGetOptionalDouble(ElemData,         (char*)"RefDZ", 0);
      check_error();
    }

    if (incl_E2) {
      h1 = atGetOptionalDouble(ElemData,            (char*)"H1", 0);
      check_error();
      h2 = atGetOptionalDouble(ElemData,            (char*)"H2", 0);
      check_error();
    }

    FringeQuadEntrance =
      atGetOptionalLong(ElemData,                   (char*)"FringeQuadEntrance",
			0);
    check_error();
    FringeQuadExit = atGetOptionalLong(ElemData,    (char*)"FringeQuadExit", 0);
    check_error();
    fringeIntM0 =
      atGetOptionalDoubleArray(ElemData,            (char*)"fringeIntM0");
    check_error();
    fringeIntP0 =
      atGetOptionalDoubleArray(ElemData,            (char*)"fringeIntP0");
    check_error();

    KickAngle = atGetOptionalDoubleArray(ElemData,   (char*)"KickAngle");
    check_error();

    if (incl_E0)
      mpole->Energy             = Energy;

    mpole->PolynomA             = PolynomA;
    mpole->PolynomB             = PolynomB;
    mpole->MaxOrder             = MaxOrder;
    mpole->NumIntSteps          = NumIntSteps;

    if (bend) {
      mpole->BendingAngle       = BendingAngle;
      mpole->EntranceAngle      = EntranceAngle;
      mpole->ExitAngle          = ExitAngle;
      mpole->FringeBendEntrance = FringeBendEntrance;
      mpole->FringeBendExit     = FringeBendExit;
      mpole->FringeInt1         = FringeInt1;
      mpole->FringeInt2         = FringeInt2;
      mpole->FullGap            = FullGap;
    } else {
      mpole->BendingAngle       = 0e0;
      mpole->EntranceAngle      = 0e0;
      mpole->ExitAngle          = 0e0;
      mpole->FringeBendEntrance = 0;
      mpole->FringeBendExit     = 0;
      mpole->FringeInt1         = 0e0;
      mpole->FringeInt2         = 0e0;
      mpole->FullGap            = 0e0;
    }

    if (cbend) {
      mpole->X0ref   = X0ref;
      mpole->ByError = ByError;
      mpole->RefDZ   = RefDZ;
    } else {
      mpole->X0ref   = 0e0;
      mpole->ByError = 0e0;
      mpole->RefDZ   = 0e0;
    }

    if (incl_E2) {
      mpole->H1      = h1;
      mpole->H2      = h2;
    }

    mpole->irho = (Elem->Length != 0e0)? mpole->BendingAngle/Elem->Length : 0e0;

    mpole->FringeQuadEntrance = FringeQuadEntrance;
    mpole->FringeQuadExit     = FringeQuadExit;
    mpole->fringeIntM0        = fringeIntM0;
    mpole->fringeIntP0        = fringeIntP0;
    mpole->KickAngle          = KickAngle;

    return Elem;
  } else
    return NULL;
}

struct elem_type* init_cav(const PyObject *ElemData, struct elem_type *Elem)
{
  double   Length, Voltage, Energy, Frequency, TimeLag;
  elem_cav *cav;

  Elem = (struct elem_type*)malloc(sizeof(struct elem_type));
  if (Elem) {
    Elem          = (struct elem_type*)malloc(sizeof(struct elem_type));
    Elem->cav_ptr = (struct elem_cav*)malloc(sizeof(struct elem_cav));
    cav           = Elem->cav_ptr;

    Length    = atGetDouble(ElemData, (char*)"Length");
    check_error();

    Voltage   = atGetDouble(ElemData, (char*)"Voltage");
    check_error();
    Energy    = atGetDouble(ElemData, (char*)"Energy");
    check_error();
    Frequency = atGetDouble(ElemData, (char*)"Frequency");
    check_error();
    TimeLag   = atGetOptionalDouble(ElemData, (char*)"TimeLag", 0);
    check_error();

    Elem->Length   = Length;

    cav->Voltage   = Voltage;
    cav->Energy    = Energy;
    cav->Frequency = Frequency;
    cav->TimeLag   = TimeLag;

    return Elem;
  } else
    return NULL;
}

struct elem_type* init_wig(const PyObject *ElemData, struct elem_type *Elem,
			   const bool rad)
{
  int
    i, Nstep, Nmeth, NHharm, NVharm;
  double
    *tmppr, kw, *By, *Bx, Lw, Bmax, Energy;
  elem_wig
    *wig;

  Elem = init_elem(ElemData, Elem, true, false);
  if (Elem) {
    Elem->wig_ptr = (struct elem_wig*)malloc(sizeof(struct elem_wig));
    wig           = Elem->wig_ptr;

    Nmeth  = atGetLong(ElemData,        (char*)"Nmeth");
    check_error();
    Nstep  = atGetLong(ElemData,        (char*)"Nstep");
    check_error();
    NHharm = atGetLong(ElemData,        (char*)"NHharm");
    check_error();
    NVharm = atGetLong(ElemData,        (char*)"NVharm");
    check_error();

    Energy = atGetDouble(ElemData,      (char*)"Energy");
    check_error();
    Lw     = atGetDouble(ElemData,      (char*)"Lw");
    check_error();
    Bmax   = atGetDouble(ElemData,      (char*)"Bmax");
    check_error();
    By     = atGetDoubleArray(ElemData, (char*)"By");
    check_error();
    Bx     = atGetDoubleArray(ElemData, (char*)"Bx");
    check_error();

    wig->Pmethod = Nmeth;
    wig->PN      = Nstep;
    wig->NHharm  = NHharm;
    wig->NVharm  = NVharm;

    wig->E0      = Energy/1e9;
    wig->Lw      = Lw;
    wig->PB0     = Bmax;

    wig->Nw      = (int)(Elem->Length/Lw);

    kw           = 2e0*PI/(wig->Lw);
    wig->Zw      = 0e0;
    wig->Aw      = 0e0;

    if (rad) {
      wig->Po = wig->E0/XMC2;
      wig->srCoef =
	(q_e*q_e)*((wig->Po)*(wig->Po)*(wig->Po))
	/(6*PI*epsilon_o*m_e*(clight*clight));
      wig->HSplitPole = 0;
      wig->VSplitPole = 0;
      wig->zStartH = 0e0;
      wig->zEndH = Elem->Length;
      wig->zStartV = 0e0;
      wig->zEndV = Elem->Length;
    }

    tmppr = By;
    for (i = 0; i < NHharm; i++) {
      tmppr++;
      wig->HCw[i] = 0e0;
      wig->HCw_raw[i] = *tmppr;
      tmppr++;
      wig->Hkx[i] = (*tmppr) * kw;
      tmppr++;
      wig->Hky[i] = (*tmppr) * kw;
      tmppr++;
      wig->Hkz[i] = (*tmppr) * kw;
      tmppr++;
      wig->Htz[i] = *tmppr;
      tmppr++;
    }

    tmppr = Bx;
    for (i = 0; i < NVharm; i++) {
      tmppr++;
      wig->VCw[i] = 0e0;
      wig->VCw_raw[i] = *tmppr;
      tmppr++;
      wig->Vkx[i] = (*tmppr) * kw;
      tmppr++;
      wig->Vky[i] = (*tmppr) * kw;
      tmppr++;
      wig->Vkz[i] = (*tmppr) * kw;
      tmppr++;
      wig->Vtz[i] = *tmppr;
      tmppr++;
    }

    for (i = NHharm; i < WHmax; i++) {
      wig->HCw[i]     = 0e0;
      wig->HCw_raw[i] = 0e0;
      wig->Hkx[i]     = 0e0;
      wig->Hky[i]     = 0e0;
      wig->Hkz[i]     = 0e0;
      wig->Htz[i]     = 0e0;
    }
    for (i = NVharm; i < WHmax; i++) {
      wig->VCw[i]     = 0e0;
      wig->VCw_raw[i] = 0e0;
      wig->Vkx[i]     = 0e0;
      wig->Vky[i]     = 0e0;
      wig->Vkz[i]     = 0e0;
      wig->Vtz[i]     = 0e0;
    }

    return Elem;
  } else
    return NULL;
}

struct elem_type* init_M66(const PyObject *ElemData, struct elem_type *Elem)
{
  double *M66;

  Elem = init_elem(ElemData, Elem, false, false);
  if (Elem) {
    Elem->M66_ptr = (struct elem_M66*)malloc(sizeof(struct elem_M66));

    M66 = atGetDoubleArray(ElemData, (char*)"M66");
    check_error();

    Elem->M66_ptr->M66 = M66;

    return Elem;
  } else
    return NULL;
}

struct elem_type* init_corr(const PyObject *ElemData, struct elem_type *Elem)
{
  double Length, *KickAngle;

  Elem = (struct elem_type*)malloc(sizeof(struct elem_type));
  if (Elem) {
    Elem->corr_ptr = (struct elem_corr*)malloc(sizeof(struct elem_corr));

    Length = atGetDouble(ElemData, "Length");
    check_error();
    KickAngle = atGetDoubleArray(ElemData, (char*)"KickAngle");
    check_error();

    Elem->Length              = Length;
    Elem->corr_ptr->KickAngle = KickAngle;

    return Elem;
  } else
    return NULL;
}

struct elem_type* init_H(const PyObject *ElemData, struct elem_type *Elem)
{
  long    max_order, num_int_steps, type, multipole_fringe;
  double  *polynom_a, *polynom_b, phi, gK;
  elem_H  *H;

  Elem = init_elem(ElemData, Elem, true, false);
  if (Elem) {
    Elem->H_ptr = (struct elem_H*)malloc(sizeof(struct elem_H));
    H           = Elem->H_ptr;

    type = atGetLong(ElemData,                     (char*)"Type");
    check_error();
    num_int_steps = atGetLong(ElemData,            (char*)"NumIntSteps");
    check_error();
    max_order = atGetLong(ElemData,                (char*)"MaxOrder");
    check_error();
    multipole_fringe = atGetOptionalLong(ElemData, (char*)"MultipoleFringe", 0);
    check_error();

    polynom_a = atGetDoubleArray(ElemData,         (char*)"PolynomA");
    check_error();
    polynom_b = atGetDoubleArray(ElemData,         (char*)"PolynomB");
    check_error();
    phi = atGetOptionalDouble(ElemData,            (char*)"BendingAngle", 0.0);
    check_error();
    gK = atGetOptionalDouble(ElemData,             (char*)"gK", 0.0);
    check_error();

    H->Type            = type;
    H->NumIntSteps     = num_int_steps;
    H->MaxOrder        = max_order;
    H->MultipoleFringe = multipole_fringe;

    H->PolynomA        = polynom_a;
    H->PolynomB        = polynom_b;
    H->BendingAngle    = phi;
    H->gK              = gK;

    return Elem;
  } else
    return NULL;
}

//------------------------------------------------------------------------------

// Tracy-2/Thor_scsi Interface -- Beam Dynamics Functions.

inline double get_p_s(const std::vector<double> &ps)
{
  double p_s, p_s2;

  if (true)
    // Small angle axproximation.
    p_s = 1e0 + ps[delta_];
  else {
    p_s2 = sqr(1e0+ps[delta_]) - sqr(ps[px_]) - sqr(ps[py_]);
    if (p_s2 >= 0e0)
      p_s = sqrt(p_s2);
    else {
      //      printf("get_p_s: *** Speed of light exceeded!\n");
      p_s = NAN;
    }
  }
  return(p_s);
}

static double B2perp(double bx, double by, double irho,
		     double x, double xpr, double y, double ypr)
/* Calculates sqr(|e x B|) , where e is a unit vector in the direction of
   velocity                                                                   */

{
  double v_norm2 = 1/(sqr(1+x*irho)+ sqr(xpr) + sqr(ypr));

  /* components of the  velocity vector
   * double ex, ey, ez;
   * ex = xpr;
   * ey = ypr;
   * ez = (1+x*irho);
   */

  return ((sqr(by*(1+x*irho)) + sqr(bx*(1+x*irho))
	   + sqr(bx*ypr - by*xpr) )*v_norm2) ;
}

void get_B2(const double h_ref, const double B[], const std::vector<double> &xp,
	    double &B2_perp, double &B2_par)
{
  // compute B_perp^2 and B_par^2
  double xn, e[3];

  xn = 1e0/sqrt(sqr(1e0+xp[x_]*h_ref)+sqr(xp[px_])+sqr(xp[py_]));
  e[X_] = xp[px_]*xn; e[Y_] = xp[py_]*xn; e[Z_] = (1e0+xp[x_]*h_ref)*xn;

  // left-handed coordinate system
  B2_perp =
    sqr(B[Y_]*e[Z_]-B[Z_]*e[Y_]) + sqr(B[X_]*e[Y_]-B[Y_]*e[X_])
    + sqr(B[Z_]*e[X_]-B[X_]*e[Z_]);

//  B2_par = sqr(B[X_]*e[X_]+B[Y_]*e[Y_]+B[Z_]*e[Z_]);
}

void radiate(std::vector<double> &ps, const double L, const double h_ref,
	     const double B[], const double E0)
{
  // M. Sands "The Physics of Electron Storage Rings" SLAC-121, p. 98.
  // ddelta/d(ds) = -C_gamma*E_0^3*(1+delta)^2*(B_perp/(Brho))^2/(2*pi)
  double              p_s0, p_s1, ds, B2_perp = 0e0, B2_par = 0e0;
  std::vector<double> cs;

  const double cl_rad = C_gamma*cube(E0)/(2e0*M_PI*1e27);

  // Large ring: x' and y' unchanged.
  p_s0 = get_p_s(ps); cs = ps; cs[px_] /= p_s0; cs[py_] /= p_s0;

  // H = -p_s => ds = H*L.
  ds = (1e0+cs[x_]*h_ref+(sqr(cs[px_])+sqr(cs[py_]))/2e0)*L;
  get_B2(h_ref, B, cs, B2_perp, B2_par);

  ps[delta_] -= cl_rad*sqr(p_s0)*B2_perp*ds;
  p_s1 = get_p_s(ps); ps[px_] = cs[px_]*p_s1; ps[py_] = cs[py_]*p_s1;

  // if (conf.emittance) is_tps<T>::emittance(conf, B2_perp, ds, p_s0, cs);
}

void Drift(double L, std::vector<double> &ps, const bool exact)
{
  double u;

  if (!exact) {
    // Small angle axproximation.
    u = L/(1e0+ps[delta_]);
    ps[x_]  += u*ps[px_]; ps[y_] += u*ps[py_];
    ps[ct_] += u*(sqr(ps[px_])+sqr(ps[py_]))/(2e0*(1e0+ps[delta_]));
  } else {
    u = L/get_p_s(ps);
    ps[x_]  += u*ps[px_]; ps[y_] += u*ps[py_];
    ps[ct_] += u*(1e0+ps[delta_]) - L;
  }
  if (false) ps[ct_] += L;
}

static double get_psi(double irho, double phi, double gap)
{
  /* Correction for magnet gap (longitudinal fringe field)

     irho h = 1/rho [1/m]
     phi  edge angle
     gap  full gap between poles

     2
     K1*gap*h*(1 + sin phi)
     psi = ----------------------- * (1 - K2*g*gap*tan phi)
     cos phi

     K1 is usually 1/2
     K2 is zero here                                                  */

  double psi;

  const double k1 = 0.5e0, k2 = 0e0;

  if (phi == 0e0)
    psi = 0e0;
  else
    psi = k1*gap*irho*(1e0+sqr(sin(phi*M_PI/180e0)))/cos(phi*M_PI/180e0)
      *(1e0 - k2*gap*irho*tan(phi*M_PI/180e0));

  return psi;
}

void EdgeFocus(const double irho, const double phi, const double gap,
	       std::vector<double> &ps)
{
  ps[px_] += irho*tan(phi*M_PI/180e0)*ps[x_];
  if (false) {
    // warning: => diverging Taylor map (see SSC-141)
    // ps[py_] -=
    //   irho*tan(phi*M_PI/180e0-get_psi(irho, phi, gap))*ps[y_]
    //   /(1e0+ps[delta_]);
    // leading order correction.
    ps[py_] -=
      irho*tan(phi*M_PI/180e0-get_psi(irho, phi, gap))*ps[y_]*(1e0-ps[delta_]);
  } else
    ps[py_] -= irho*tan(phi*M_PI/180e0-get_psi(irho, phi, gap))*ps[y_];
}

void thin_kick(const int Order, const double MB[], const double L,
	       const double h_bend, const double h_ref, std::vector<double> &ps,
	       const double E0, const bool rad)
{
  // The vector potential for the combined-function sector bend is from:
  // C. Iselin "Lie Transformations and Transport Equations for Combined-
  // Function Dipoles" Part. Accel. 17, 143-155 (1985).
  int                 j;
  double              BxoBrho, ByoBrho, ByoBrho1, B[3], u, p_s;
  std::vector<double> ps0;

  if ((h_bend != 0e0) || ((1 <= Order) && (Order <= HOMmax))) {
    ps0 = ps;
    // Compute magnetic field with Horner's rule.
    ByoBrho = MB[Order+HOMmax]; BxoBrho = MB[HOMmax-Order];
    for (j = Order-1; j >= 1; j--) {
      ByoBrho1 = ps0[x_]*ByoBrho - ps0[y_]*BxoBrho + MB[j+HOMmax];
      BxoBrho  = ps0[y_]*ByoBrho + ps0[x_]*BxoBrho + MB[HOMmax-j];
      ByoBrho  = ByoBrho1;
    }

    if (rad) {
      B[X_] = BxoBrho; B[Y_] = ByoBrho + h_bend; B[Z_] = 0e0;
      radiate(ps, L, h_ref, B, E0);
    }

    if (h_ref != 0e0) {
      // Sector bend.
      if (true) {
	ps[px_] -= L*(ByoBrho+(h_bend-h_ref)/2e0+h_ref*h_bend*ps0[x_]
		      -h_ref*ps0[delta_]);
	ps[ct_] += L*h_ref*ps0[x_];
      } else {
	// The Hamiltonian is split into: H_d + H_k; with [H_d, H_d] = 0.
	p_s = get_p_s(ps0); u = L*h_ref*ps0[x_]/p_s;
	ps[x_]  += u*ps0[px_];
	ps[y_]  += u*ps0[py_];
	ps[ct_] += u*(1e0+ps0[delta_]);
	// ps[px_] -= L*(h_bend*(1e0+h_ref*ps0[x_])-h_ref*p_s);
	// Field expansion up to sextupole like terms.
	ByoBrho += h_bend - MB[Quad+HOMmax]*h_ref*sqr(ps0[y_])/2e0;
	ps[px_] -= L*((1e0+h_ref*ps0[x_])*ByoBrho-h_ref*p_s);
	ps[py_] += L*(1e0+h_ref*ps0[x_])*BxoBrho;
      }
    } else
      // Cartesian bend.
      ps[px_] -= L*(h_bend+ByoBrho);
    ps[py_] += L*BxoBrho;
  }
}

void Cav_Pass(const double L, const double f_RF, const double V_RFoE0,
	      const double phi, std::vector<double> &ps)
{
  double delta;

  Drift(L/2e0, ps, false);
  if (V_RFoE0 != 0e0) {
    delta = -V_RFoE0*sin(2e0*M_PI*f_RF/C0*ps[ct_]+phi);
    ps[delta_] += delta;

    // if (globval.radiation) globval.dE -= is_double<T>::cst(delta);

    // if (globval.pathlength) ps[ct_] -= C->Ph/C->Pfreq*c0;
  }
  Drift(L/2e0, ps, false);
}

//------------------------------------------------------------------------------

inline void Drift(double ps[], const double L, const bool exact)
{
  std::vector<double> ps_stl = arrtostl(ps);

  Drift(L, ps_stl, false);
  stltoarr(ps_stl, ps);
}

void thin_kick(double ps[], const double a[], const double b[],
	       const double L, const double irho, const int n_max,
	       const double E0, const bool rad)
{
  double              bn[2*HOMmax+1];
  std::vector<double> ps_stl = arrtostl(ps);

  for (int k = n_max+1; k > 0; k--) {
    bn[HOMmax+k] = b[k-1];
    bn[HOMmax-k] = a[k-1];
  }
  thin_kick(n_max+1, bn, L, irho, irho, ps_stl, E0, rad);
  stltoarr(ps_stl, ps);
}

//------------------------------------------------------------------------------

void IdentityPass(double ps[], const int num_particles,
		  const struct elem_type *Elem)
{
  int    k;
  double *ps_vec;

  for (k = 0; k < num_particles; k++) {	/*Loop over particles  */
    ps_vec = ps+k*PS_DIM;
    if (!atIsNaN(ps_vec[0])) {
      /*  misalignment at entrance  */
      if (Elem->T1) ATaddvv(ps_vec, Elem->T1);
      if (Elem->R1) ATmultmv(ps_vec, Elem->R1);
      /* Check physical apertures */
      if (Elem->RApertures) checkiflostRectangularAp(ps_vec, Elem->RApertures);
      if (Elem->EApertures) checkiflostEllipticalAp(ps_vec, Elem->EApertures);
      /* Misalignment at exit */
      if (Elem->R2) ATmultmv(ps_vec, Elem->R2);
      if (Elem->T2) ATaddvv(ps_vec, Elem->T2);
    }
  }
}

//------------------------------------------------------------------------------

void AperturePass(double ps[], int num_particles, const struct elem_type *Elem)
{
  /* Checks X and Y of each input 6-vector and marks the corresponding element
     in lossflag array with 0 if X,Y are exceed the limits given by limitsptr
     array limitsptr has 4 elements: (MinX, MaxX, MinY, MaxY)                 */
  int    k;
  double *ps_vec;

  for (k = 0; k < num_particles; k++) {
    ps_vec = ps+k*PS_DIM;
    if (!atIsNaN(ps_vec[0])) {
	/*  check if this particle is already marked as lost */
	checkiflostRectangularAp(ps_vec, Elem->ap_ptr->limits);
    }
  }
}

//------------------------------------------------------------------------------

void DriftPass(double *ps, const int num_particles,
	       const struct elem_type *Elem)
/* le - physical length
   ps - 6-by-N matrix of initial conditions reshaped into
   1-d array of 6*N elements                                                  */
{
  int    k;
  double *ps_vec;

#pragma omp parallel for if (num_particles > OMP_PARTICLE_THRESHOLD*10)    \
  default(shared) shared(ps, num_particles) private(k, ps_vec)

  for (k = 0; k < num_particles; k++) { /*Loop over particles  */
    ps_vec = ps+k*PS_DIM;
    if(!atIsNaN(ps_vec[0])) {
      /*  misalignment at entrance  */
      if (Elem->T1) ATaddvv(ps_vec, Elem->T1);
      if (Elem->R1) ATmultmv(ps_vec, Elem->R1);
      /* Check physical apertures at the entrance of the magnet */
      if (Elem->RApertures) checkiflostRectangularAp(ps_vec, Elem->RApertures);
      if (Elem->EApertures) checkiflostEllipticalAp(ps_vec, Elem->EApertures);
      Drift(ps_vec, Elem->Length, false);
      /* Check physical apertures at the exit of the magnet */
      if (Elem->RApertures) checkiflostRectangularAp(ps_vec, Elem->RApertures);
      if (Elem->EApertures) checkiflostEllipticalAp(ps_vec, Elem->EApertures);
      /* Misalignment at exit */
      if (Elem->R2) ATmultmv(ps_vec, Elem->R2);
      if (Elem->T2) ATaddvv(ps_vec, Elem->T2);
    }
  }
}

//------------------------------------------------------------------------------

void CorrectorPass(double ps[], const int num_particles,
		   const struct elem_type *Elem)
/* xkick, ykick - horizontal and vertical kicks in radiand
   r - 6-by-N matrix of initial conditions reshaped into
   1-d array of 6*N elements                                                  */
{
  int    k;
  double  xkick, ykick, *ps_vec = nullptr;

  if (Elem->Length == 0e0)
    for(k = 0; k < num_particles; k++) {
      ps_vec = ps+k*PS_DIM;
      if(!atIsNaN(ps[0])) {
	ps[1] += Elem->corr_ptr->KickAngle[0];
	ps[3] += Elem->corr_ptr->KickAngle[1];
      }
    }
  else {

#pragma omp parallel for if (num_particles > OMP_PARTICLE_THRESHOLD)    \
  default(none) shared(ps, num_particles, len, xkick, ykick) private(c)

    xkick = Elem->corr_ptr->KickAngle[0];
    ykick = Elem->corr_ptr->KickAngle[1];
    for(k = 0; k < num_particles; k++) {
      ps_vec = ps+k*PS_DIM;
      if(!atIsNaN(ps[0])) {
	double
	  p_norm = 1/(1+ps[4]),
	  NormL  = Elem->Length*p_norm;

	ps[5] +=
	  NormL*p_norm*(xkick*xkick/3 + ykick*ykick/3 + ps[1]*ps[1]
			+ ps[3]*ps[3] + ps[1]*xkick + ps[3]*ykick)/2;

	ps[0] += NormL*(ps[1]+xkick/2);
	ps[1] += xkick;
	ps[2] += NormL*(ps[3]+ykick/2);
	ps[3] += ykick;
      }
    }
  }
}

//------------------------------------------------------------------------------

void edge_fringe(double ps[], const double inv_rho,
		 const double edge_angle, const double fint,
		 const double gap, const int method, const bool hor)
{
  std::vector<double> ps_stl = arrtostl(ps);

  EdgeFocus(inv_rho, edge_angle*180e0/M_PI, gap, ps_stl);
  stltoarr(ps_stl, ps);
}

static void edge_fringe_entrance(double ps[], double inv_rho, double edge_angle,
				 double fint, double gap, int method)
{ edge_fringe(ps, inv_rho, edge_angle, fint, gap, method, true); }

static void edge_fringe_exit(double ps[], double inv_rho, double edge_angle,
			     double fint, double gap, int method)
{ edge_fringe(ps, inv_rho, edge_angle, fint, gap, method, false); }

static void QuadFringePassP(double ps[], const double b2)
{
  /* x=ps[0],px=ps[1],y=ps[2],py=ps[3],delta=ps[4],ct=ps[5]
     Lee-Whiting's thin lens limit formula as given in p. 390 of "Beam
     Dynamics..."by E. Forest                                                 */
  double u     = b2/(12.0*(1.0+ps[4]));
  double x2    = ps[0]*ps[0];
  double z2    = ps[2]*ps[2];
  double xz    = ps[0]*ps[2];
  double gx    = u * (x2+3*z2) * ps[0];
  double gz    = u * (z2+3*x2) * ps[2];
  double r1tmp = 0;
  double r3tmp = 0;

  ps[0] += gx;
  r1tmp = 3*u*(2*xz*ps[3]-(x2+z2)*ps[1]);

  ps[2] -= gz;

  r3tmp = 3*u*(2*xz*ps[1]-(x2+z2)*ps[3]);
  ps[5] -= (gz*ps[3] - gx*ps[1])/(1+ps[4]);

  ps[1] += r1tmp;
  ps[3] -= r3tmp;
}

static void QuadFringePassN(double ps[], const double b2)
{
  /* x=ps[0],px=ps[1],y=ps[2],py=ps[3],delta=ps[4],ct=ps[5]
     Lee-Whiting's thin lens limit formula as given in p. 390 of "Beam
     Dynamics..."by E. Forest                                                 */
  double u     = b2/(12.0*(1.0+ps[4]));
  double x2    = ps[0]*ps[0];
  double z2    = ps[2]*ps[2];
  double xz    = ps[0]*ps[2];
  double gx    = u * (x2+3*z2) * ps[0];
  double gz    = u * (z2+3*x2) * ps[2];
  double r1tmp = 0;
  double r3tmp = 0;

  ps[0] -= gx;
  r1tmp = 3*u*(2*xz*ps[3]-(x2+z2)*ps[1]);

  ps[2] += gz;

  r3tmp = 3*u*(2*xz*ps[1]-(x2+z2)*ps[3]);
  ps[5] += (gz*ps[3] - gx*ps[1])/(1+ps[4]);

  ps[1] -= r1tmp;
  ps[3] += r3tmp;
}

/* from elegant code */
static void quadPartialFringeMatrix(double R[6][6], double K1, double inFringe,
				    double *fringeInt, int part)
{
  double J1x, J2x, J3x, J1y, J2y, J3y;
  double K1sqr, expJ1x, expJ1y;

  R[4][4] = R[5][5] = 1;

  K1sqr = K1*K1;

  if (part==1) {
    J1x = inFringe*(K1*fringeInt[1] - 2*K1sqr*fringeInt[3]/3.);
    J2x = inFringe*(K1*fringeInt[2]);
    J3x = inFringe*(K1sqr*(fringeInt[2] + fringeInt[4]));

    K1  = -K1;
    J1y = inFringe*(K1*fringeInt[1] - 2*K1sqr*fringeInt[3]/3.);
    J2y = -J2x;
    J3y = J3x;
  } else {
    J1x = inFringe*(K1*fringeInt[1] + K1sqr*fringeInt[0]*fringeInt[2]/2);
    J2x = inFringe*(K1*fringeInt[2]);
    J3x = inFringe*(K1sqr*(fringeInt[4]-fringeInt[0]*fringeInt[1]));

    K1  = -K1;
    J1y = inFringe*(K1*fringeInt[1] + K1sqr*fringeInt[0]*fringeInt[2]);
    J2y = -J2x;
    J3y = J3x;
  }

  expJ1x  = R[0][0] = exp(J1x);
  R[0][1] = J2x/expJ1x;
  R[1][0] = expJ1x*J3x;
  R[1][1] = (1 + J2x*J3x)/expJ1x;

  expJ1y  = R[2][2] = exp(J1y);
  R[2][3] = J2y/expJ1y;
  R[3][2] = expJ1y*J3y;
  R[3][3] = (1 + J2y*J3y)/expJ1y;

  return;
}

static void linearQuadFringeElegantEntrance
(double* r6, double b2, double *fringeIntM0, double *fringeIntP0)
{
  double R[6][6];
  double *fringeIntM, *fringeIntP;
  double delta, inFringe;
  /* quadrupole linear fringe field, from elegant code */
  inFringe = -1.0;
  fringeIntM = fringeIntP0;
  fringeIntP = fringeIntM0;
  delta = r6[4];
  /* determine first linear matrix for this delta */
  quadPartialFringeMatrix(R, b2/(1+delta), inFringe, fringeIntM, 1);
  r6[0] = R[0][0]*r6[0] + R[0][1]*r6[1];
  r6[1] = R[1][0]*r6[0] + R[1][1]*r6[1];
  r6[2] = R[2][2]*r6[2] + R[2][3]*r6[3];
  r6[3] = R[3][2]*r6[2] + R[3][3]*r6[3];
  /* nonlinear fringe field */
  QuadFringePassP(r6,b2);   /*This is original AT code*/
  /*Linear fringe fields from elegant*/
  inFringe=-1.0;
  /* determine and apply second linear matrix, from elegant code */
  quadPartialFringeMatrix(R, b2/(1+delta), inFringe, fringeIntP, 2);
  r6[0] = R[0][0]*r6[0] + R[0][1]*r6[1];
  r6[1] = R[1][0]*r6[0] + R[1][1]*r6[1];
  r6[2] = R[2][2]*r6[2] + R[2][3]*r6[3];
  r6[3] = R[3][2]*r6[2] + R[3][3]*r6[3];
}

static void linearQuadFringeElegantExit
(double* r6, double b2, double *fringeIntM0, double *fringeIntP0)
{
  double R[6][6];
  double *fringeIntM, *fringeIntP;
  double delta, inFringe;
  /* quadrupole linear fringe field, from elegant code */
  inFringe=1.0;
  fringeIntM = fringeIntM0;
  fringeIntP = fringeIntP0;
  delta = r6[4];
  /* determine first linear matrix for this delta */
  quadPartialFringeMatrix(R, b2/(1+delta), inFringe, fringeIntM, 1);
  r6[0] = R[0][0]*r6[0] + R[0][1]*r6[1];
  r6[1] = R[1][0]*r6[0] + R[1][1]*r6[1];
  r6[2] = R[2][2]*r6[2] + R[2][3]*r6[3];
  r6[3] = R[3][2]*r6[2] + R[3][3]*r6[3];
  /* nonlinear fringe field */
  QuadFringePassN(r6, b2);   /*This is original AT code*/
  /*Linear fringe fields from elegant*/
  inFringe=1.0;
  /* determine and apply second linear matrix, from elegant code */
  quadPartialFringeMatrix(R, b2/(1+delta), inFringe, fringeIntP, 2);
  r6[0] = R[0][0]*r6[0] + R[0][1]*r6[1];
  r6[1] = R[1][0]*r6[0] + R[1][1]*r6[1];
  r6[2] = R[2][2]*r6[2] + R[2][3]*r6[3];
  r6[3] = R[3][2]*r6[2] + R[3][3]*r6[3];
}

void MpolePass(double ps[], const int num_particles,
	       const struct elem_type *Elem, const bool rad)
{
  int    k, m;
  double *ps_vec;

  const elem_mpole *mpole = Elem->mpole_ptr;

  const bool
    useLinFrEleEntrance =
    (mpole->fringeIntM0 != NULL && mpole->fringeIntP0 != NULL
     && mpole->FringeQuadEntrance == 2),
    useLinFrEleExit =
    (mpole->fringeIntM0 != NULL && mpole->fringeIntP0 != NULL
     && mpole->FringeQuadExit == 2);

  const double
    SL = Elem->Length/mpole->NumIntSteps,
    L1 = SL*DRIFT1,
    L2 = SL*DRIFT2,
    K1 = SL*KICK1,
    K2 = SL*KICK2;

  if (mpole->KickAngle) {
    /* Convert corrector component to polynomial coefficients */
    mpole->PolynomB[0] -= sin(mpole->KickAngle[0])/Elem->Length;
    mpole->PolynomA[0] += sin(mpole->KickAngle[1])/Elem->Length;
  }

#pragma omp parallel for if (num_particles > OMP_PARTICLE_THRESHOLD)	\
  default(none)								\
  shared(r, num_particles, Elem)					\
  private(k)

  for (k = 0; k < num_particles; k++) {	/* Loop over particles  */
    ps_vec = ps+k*PS_DIM;
    if (!atIsNaN(ps_vec[0])) {
      /*  misalignment at entrance  */
      if (Elem->T1) ATaddvv(ps_vec, Elem->T1);
      if (Elem->R1) ATmultmv(ps_vec, Elem->R1);

      /* Check physical apertures at the entrance of the magnet */
      if (Elem->RApertures) checkiflostRectangularAp(ps_vec, Elem->RApertures);
      if (Elem->EApertures) checkiflostEllipticalAp(ps_vec, Elem->EApertures);

      if (mpole->irho != 0e0) {
	/* edge focus */
	edge_fringe_entrance
	  (ps_vec, mpole->irho, mpole->EntranceAngle, mpole->FringeInt1,
	   mpole->FullGap, mpole->FringeBendEntrance);
      }

      /* quadrupole gradient fringe entrance*/
      if (mpole->FringeQuadEntrance && mpole->PolynomB[1]!=0) {
	if (useLinFrEleEntrance) /*Linear fringe fields from elegant*/
	  linearQuadFringeElegantEntrance
	    (ps_vec, mpole->PolynomB[1], mpole->fringeIntM0,
	     mpole->fringeIntP0);
	else
	  QuadFringePassP(ps_vec, mpole->PolynomB[1]);
      }

      /* integrator */
      for (m = 0; m < mpole->NumIntSteps; m++) { /* Loop over slices*/
	Drift(ps_vec, L1, false);
	thin_kick(ps_vec, mpole->PolynomA, mpole->PolynomB, K1, mpole->irho,
		  mpole->MaxOrder, mpole->Energy, rad);
	Drift(ps_vec, L2, false);
	thin_kick(ps_vec, mpole->PolynomA, mpole->PolynomB, K2, mpole->irho,
		    mpole->MaxOrder, mpole->Energy, rad);
	Drift(ps_vec, L2, false);
	thin_kick(ps_vec, mpole->PolynomA, mpole->PolynomB, K1, mpole->irho,
		    mpole->MaxOrder, mpole->Energy, rad);
	Drift(ps_vec, L1, false);
      }

      /* quadrupole gradient fringe */
      if (mpole->FringeQuadExit && mpole->PolynomB[1]!=0) {
	if (useLinFrEleExit) /*Linear fringe fields from elegant*/
	  linearQuadFringeElegantExit
	    (ps_vec, mpole->PolynomB[1], mpole->fringeIntM0,
	     mpole->fringeIntP0);
	else
	  QuadFringePassN(ps_vec, mpole->PolynomB[1]);
      }

      if (mpole->irho != 0e0) {
	/* edge focus */
	edge_fringe_exit
	  (ps_vec, mpole->irho, mpole->ExitAngle, mpole->FringeInt2,
	   mpole->FullGap, mpole->FringeBendExit);
      }

      /* Check physical apertures at the exit of the magnet */
      if (Elem->RApertures) checkiflostRectangularAp(ps_vec, Elem->RApertures);
      if (Elem->EApertures) checkiflostEllipticalAp(ps_vec, Elem->EApertures);

      /* Misalignment at exit */
      if (Elem->R2) ATmultmv(ps_vec, Elem->R2);
      if (Elem->T2) ATaddvv(ps_vec, Elem->T2);
    }
  }
  if (mpole->KickAngle) {
    /* Remove corrector component in polynomial coefficients */
    mpole->PolynomB[0] += sin(mpole->KickAngle[0])/Elem->Length;
    mpole->PolynomA[0] -= sin(mpole->KickAngle[1])/Elem->Length;
  }
}

//------------------------------------------------------------------------------

void cav_pass(double ps[], const double L, const double V_RFoE0,
	      const double f_RF, const double lag)
{
  std::vector<double> ps_stl = arrtostl(ps);

  const double phi = -lag*2e0*M_PI*f_RF/C0;

  Cav_Pass(L, f_RF, V_RFoE0, phi, ps_stl);
  stltoarr(ps_stl, ps);
}

void CavityPass(double ps[], const int num_particles,
		const struct elem_type *Elem)
{
  int    k;
  double *ps_vec;

  const elem_cav *cav = Elem->cav_ptr;

  for(k = 0; k < num_particles; k++) {
    ps_vec = ps+k*PS_DIM;
    if(!atIsNaN(ps_vec[0]))
      cav_pass(ps_vec,  Elem->Length, cav->Voltage/cav->Energy, cav->Frequency,
	       cav->TimeLag);
  }
}

//------------------------------------------------------------------------------

void Matrix66Pass(double ps[], const int num_particles,
		  const struct elem_type *Elem)
{
  int    k;
  double *ps_vec;

  for (k = 0; k < num_particles; k++) {	/*Loop over particles  */
    ps_vec = ps+k*PS_DIM;
    if (!atIsNaN(ps_vec[0])) {
      if (Elem->T1 != NULL) ATaddvv(ps_vec, Elem->T1);
      if (Elem->R1 != NULL) ATmultmv(ps_vec, Elem->R1);
      ATmultmv(ps_vec, Elem->M66_ptr->M66);
      if (Elem->R2 != NULL) ATmultmv(ps_vec, Elem->R2);
      if (Elem->T2 != NULL) ATaddvv(ps_vec, Elem->T2);
    }
  }
}

//------------------------------------------------------------------------------

/* track.cc
   tracking routines for exact Hamiltonian from Forest / PTC / Tracy-3
   James Rowland 2010

   Exact integrator for different element types

   This method will work for a drift, a quadrupole, a sextupole or a
   bending magnet. It distinguishes between these using the Class field
   on the element.

   The 'ExactHamiltonianPass' method uses the square root hamiltonian in
   cartesian co-ordinates (see other notes for derivation).
   This is equivalent to setting exact=true in MADX-PTC.
   Multipole fringe fields are also enabled for quadrupoles
   (fringe = true option in MADX-PTC).

   Note that the PolynomB array in the exact cartesian rectangular bend
   refers to the normalized straight multipole components of the vector
   potential, so PolynomB(1) should be set to 1/rho (B_bend / Brho).
   The conjugate momenta in the curvilinear co-ordinate system are not
   the same as in the cartesian system so PolynomB(1) must be set back
   to zero when using a curvilinear symplectic integrator method such
   as the 'BndMPoleSymplectic4E2Pass'. See Forest p362 for a detailed
   explanation of the vector potential in curvilinear co-ordinates.           */

#undef DEBUG_MODE

#ifdef DEBUG_MODE
#define Log(x) printf x
#else
#define Log(x)
#endif

/* Forest-Ruth 4th order coefficients
   could also use 6th order Yoshida */

#define INT_ORDER 4

/* Generated by:

   #!/usr/bin/env python
   c1 = c4 = 1.0/(2.0*(2.0-2.0**(1.0/3.0)))
   c2 = c3 = (1-2**(1.0/3.0))/(2.0*(2.0-2.0**(1.0/3.0)))
   d1 = d3 = 1.0/(2.0-2.0**(1.0/3.0))
   d2 = -(2**(1.0/3.0)/(2.0-2.0**(1.0/3.0)))
   d4 = 0
   print "double c[] = {% .17f, % .17f, % .17f, % .17f};" % (c1, c2, c3, c4)
   print "double d[] = {% .17f, % .17f, % .17f, % .17f};" % (d1, d2, d3, d4)  */

double
c[] =
  { 0.67560359597982889, -0.17560359597982883, -0.17560359597982883,
    0.67560359597982889},
  d[] =
  { 1.35120719195965777, -1.70241438391931532,  1.35120719195965777,
    0.00000000000000000};


template<typename T>
void mpole_fringe(T *x, const int type, const double LR, const double *F,
		  const int nF, const int edge)
{
  // PTC mpole_fringer
  // Forest 13.29
  // not re-derived and checked
  // note this is the sum over n of Forest 13.29
  // one for each multipole component

  T
    I, U, V, DU, DV, DUX, DVX, DUY, DVY,
    FX, FY, FX_X, FX_Y, FY_X, FY_Y,
    RX, IX, DRX, DIX;

  if (edge == 0)
    I = 1;
  else
    I = -1;

  FX   = 0;
  FY   = 0;
  FX_X = 0;
  FX_Y = 0;
  FY_X = 0;
  FY_Y = 0;

  RX   = 1.0;
  IX   = 0.0;

  // invariant is (j is the index, i is the complex unit)
  // RX+IXi = (x + iy)^j
  for(int n = 0; n < nF; n++) {
    int j = n + 1;
    double
      B = F[2*n],
      A = F[2*n+1];

    DRX = RX;
    DIX = IX;

    // complex muls
    RX = DRX * x[x_] - DIX * x[y_];
    IX = DRX * x[y_] + DIX * x[x_];

    if(j == 1 && type == dipole) {
      U  =         - A * IX;
      V  =         + A * RX;
      DU =         - A * DIX;
      DV =         + A * DRX;
    } else {
      U  = B * RX  - A * IX;
      V  = B * IX  + A * RX;
      DU = B * DRX - A * DIX;
      DV = B * DIX + A * DRX;
    }

    T f1 = -I / 4.0 / (j + 1);

    U  = U  * f1;
    V  = V  * f1;
    DU = DU * f1;
    DV = DV * f1;

    DUX =  j * DU;
    DVX =  j * DV;
    DUY = -j * DV;
    DVY =  j * DU;

    double nf = 1.0 * (j + 2) / j;

    FX += U * x[x_] + nf * V * x[y_];
    FY += U * x[y_] - nf * V * x[x_];

    FX_X += DUX * x[x_] + U      + nf * x[y_] * DVX;
    FX_Y += DUY * x[x_] + nf * V + nf * x[y_] * DVY;

    FY_X += DUX * x[y_] - nf * V - nf * x[x_] * DVX;
    FY_Y += DUY * x[y_] + U      - nf * x[x_] * DVY;
  }

  T DEL = 1.0 / (1 + x[delta_]);

  // solve 2x2 matrix equation

  T A = 1 -FX_X * DEL;
  T B =   -FY_X * DEL;
  T D = 1 -FY_Y * DEL;
  T C =   -FX_Y * DEL;

  x[x_] = x[x_] - FX * DEL;
  x[y_] = x[y_] - FY * DEL;

  T pxf = (D * x[px_] - B * x[py_]) / (A * D - B * C);
  T pyf = (A * x[py_] - C * x[px_]) / (A * D - B * C);
  x[py_]      = pyf;
  x[px_]      = pxf;
  x[ct_]      = x[ct_] - (x[px_] * FX + x[py_] * FY) * DEL * DEL;
}

template<typename T>
T pow2(T x) { return x * x; }

/* this is the z momentum */
template<typename T>
T get_pz(T * x)
{ return sqrt(pow2(1 + x[delta_]) - pow2(x[px_]) - pow2(x[py_])); }

/* Forest 10.26, layout rotation
   phi: angle [rad]              */
template <typename T> void Yrot(double phi, T * x)
{
  T c, s;

  c = cos(phi);
  s = sin(phi);
  T x1[6] = {x[0], x[1], x[2], x[3], x[4], x[5]};
  T ps = get_pz(x);
  T p = c*ps - s*x1[px_];
  x[x_] = x1[x_]*ps/p;
  x[px_] = s*ps + c*x1[px_];
  x[y_] += x1[x_]*x1[py_]*s/p;
  x[ct_] += (1.0+x1[delta_])*x1[x_]*s/p;
}

/* Forest 10.23, exact drift
   L: length [m]              */
template<typename T>
void exact_drift(T * x, double L)
{
  T u = L / get_pz(x);

  x[x_] += x[px_] * u;
  x[y_] += x[py_] * u;
  x[ct_] += u * (1.0 + x[delta_]);
}

/* Forest-Ruth 4th order integrator
   x      : phase space (inout)
   L      : length
   F      : multipole coefficients
   nf     : length of F
   slices : number of integration steps  */

template<typename T>
void fr4(T *x, const double L, const double *F, const int nF, const int slices)
{
  int    max_order = nF - 1, s, n, i;
  double ds        = L / slices;

  for(s = 0; s < slices; s++) {
    for(n = 0; n < INT_ORDER; n++) {
      exact_drift(x, c[n] * ds);

      /* multipole summation with horner's rule
	 scaled field = sum_n (b_n+ia_n) (x+iy)^n */

      /* C99 complex numbers don't work in C++
	 forget the C++ complex class
	 complex double f = F[max_order];
	 complex double z = x[x_] + x[y_] * I;    */

      T fr = F[2*max_order];
      T fi = F[2*max_order+1];

      for(i = max_order - 1; i >= 0; i--) {
	/* complex multiplication
	   f = f * z + F[i];       */
	T temp1 = fr * x[x_] - fi * x[y_];
	T temp2 = fr * x[y_] + fi * x[x_];
	fr = temp1 + F[2*i];
	fi = temp2 + F[2*i+1];
      }

      x[px_] -= d[n] * ds *  fr;
      x[py_] -= d[n] * ds * -fi;
    }
  }
}

/* bend fringe */

template <typename T>
T Sec(T x) { return 1.0 / cos(x); }

template<typename T>
void bend_fringe(T *x, double irho, double gK)
{
  T dpx, dpy, dd, b0, px, py, pz, g, K, d, phi, xp, yp, yf, xf, lf, pyf;

  b0 = irho;

  /* gK always multiplied together so put everything in g and set K to one */

  K = 1.0;
  g = gK;

  pz = get_pz(x);
  px = x[px_];
  py = x[py_];
  d  = x[delta_];
  xp = px / pz;
  yp = py / pz;

  phi =
    -b0 * tan( b0 * g * K * (1 + pow2(xp)*(2 + pow2(yp)))*pz
	       - atan(xp / (1 + pow2(yp))));

  /* these are the partial derivatives of phi with respect to px, py and delta
     total horror from Mathematica. This could benefit from some mini-TPSA */

  dpx =
    -((b0*(pow(px,2)*pow(pz,4)
	   *(pow(py,2) - pow(pz,2))
	   - pow(pz,6)*(pow(py,2) + pow(pz,2)) +
	   b0*g*K*px*(pow(pz,2)*pow(pow(py,2) + pow(pz,2),2)
		      *(2*pow(py,2) + 3*pow(pz,2))
		      + pow(px,4)*(3*pow(py,2)*pow(pz,2) + 2*pow(pz,4)) +
		      pow(px,2)*(3*pow(py,6)
				   + 8*pow(py,4)*pow(pz,2)
				   + 9*pow(py,2)*pow(pz,4)
				   + 5*pow(pz,6))))*pow(Sec((b0*g*K*(pow(pz,4)
			  + pow(px,2)*(pow(py,2) + 2*pow(pz,2))))/pow(pz,3)
		 - atan((px*pz)/(pow(py,2) + pow(pz,2)))),2))/
      (pow(pz,5)*(pow(py,4) + pow(px,2)*pow(pz,2) + 2*pow(py,2)*pow(pz,2)
		  + pow(pz,4))));

  dpy =
    -((b0*py*(px*pow(pz,4)*(pow(py,2) + pow(pz,2))
	      + b0*g*K*(-(pow(pz,4)*pow(pow(py,2) + pow(pz,2),2))
			+ pow(px,4)*(3*pow(py,2)*pow(pz,2) + 4*pow(pz,4))
			+ pow(px,2)*(3*pow(py,6)
				       + 10*pow(py,4)*pow(pz,2)
				       + 11*pow(py,2)*pow(pz,4)
				       + 3*pow(pz,6))))*
       pow(Sec((b0*g*K*(pow(pz,4) + pow(px,2) *(pow(py,2) + 2*pow(pz,2))))
	       /pow(pz,3) - atan((px*pz)/(pow(py,2) + pow(pz,2)))),2))/
      (pow(pz,5)*(pow(py,4) + pow(px,2)*pow(pz,2) + 2*pow(py,2)*pow(pz,2)
		  + pow(pz,4))));

  dd =
    (b0*(1 + d)*(px*pow(pz,4)*(pow(py,2) - pow(pz,2)) + b0*g*K*
		 (-(pow(pz,4)*pow(pow(py,2) + pow(pz,2),2))
		  + pow(px,4)*(3*pow(py,2)*pow(pz,2) + 2*pow(pz,4))
		  + pow(px,2)*(3*pow(py,6) + 8*pow(py,4)*pow(pz,2)
				 + 7*pow(py,2)*pow(pz,4) + pow(pz,6))))
     *pow(Sec((b0*g*K*(pow(pz,4)
			 + pow(px,2)*(pow(py,2) + 2*pow(pz,2))))/pow(pz,3)
		- atan((px*pz)/(pow(py,2) + pow(pz,2)))),2))/
    (pow(pz,5)*(pow(py,4) + pow(px,2)*pow(pz,2)
		  + 2*pow(py,2)*pow(pz,2) + pow(pz,4)));

  /* solve quadratic equation in yf (Forest fringe_part_I.pdf) */

  yf = (2 * x[y_]) / (1 + sqrt(1 - 2 * dpy * x[y_]));
  xf = x[x_] + 0.5 * dpx * pow2(yf);
  lf = x[ct_] - 0.5 * dd * pow2(yf);
  pyf = py - phi * yf;

  x[y_]  = yf;
  x[x_]  = xf;
  x[py_] = pyf;
  x[ct_] = lf;
}

template <typename T>
void bend(T *x, const int type, const double L, const double phi,
	  const double gK, const double *F, const int nF, const int slices,
	  const bool mp_fringe)
{
  double
    irho = phi/L,
    /* convert arc length to rectangular length */
    LR = 2e0/irho*sin(phi/2e0);

  Yrot(phi/2e0, x);
  bend_fringe(x, F[0], gK);
  if(mp_fringe) mpole_fringe(x, type, LR, F, nF, 0);
  fr4(x, LR, F, nF, slices);
  if(mp_fringe) mpole_fringe(x, type, LR, F, nF, 1);
  bend_fringe(x, -F[0], gK);
  Yrot(phi/2e0, x);
}

template<typename T>
void track_element(T *x, const elem_type *Elem)
{
  const elem_H *H = Elem->H_ptr;

  Log(("track element\n"));
  switch(H->Type) {
  case element_type(drift):
    Log(("drift %f\n", Elem->Length));
    exact_drift(x, Elem->Length);
    x[ct_] -= Elem->Length;
    break;
  case element_type(dipole):
    Log(("bend %f %f %f\n", Elem->Length, H->BendingAngle,
	 creal(Elem->F[0])));
    bend(x, H->Type, Elem->Length, H->BendingAngle, H->gK, H->F, H->MaxOrder,
	 H->NumIntSteps, H->MultipoleFringe);
    x[ct_] -= Elem->Length;
    break;
  case element_type(multipole):
    Log(("multipole %f %f\n", Elem->Length, creal(H->PolynomB[1])));
    if(H->MultipoleFringe)
      mpole_fringe(x, H->Type, Elem->Length, H->F, H->MaxOrder, 0);
    fr4(x, Elem->Length, H->F, H->MaxOrder, H->NumIntSteps);
    if(H->MultipoleFringe)
      mpole_fringe(x, H->Type, Elem->Length, H->F, H->MaxOrder, 1);
    x[ct_] -= Elem->Length;
    break;
  case element_type(marker):
    Log(("marker\n"));
    break;
  default:
    Log(("unknown element\n"));
    exit(1);
  }
}

void track_map(double * x, lattice * lat, double * map1)
{
  fprintf(stderr,
	  "track_map is not available, rebuild with #define TPSA_MODE\n");
  exit(1);
}

void HamPass(double ps[], const int num_particles, const struct elem_type *Elem)
{
  int    k;
  double *ps_vec;

  for(k = 0; k < Elem->H_ptr->MaxOrder; k++) {
    Elem->H_ptr->F[2*k]   = Elem->H_ptr->PolynomB[k];
    Elem->H_ptr->F[2*k+1] = Elem->H_ptr->PolynomA[k];
  }

  for(k = 0; k < num_particles; k++) {
    ps_vec = ps+k*PS_DIM;
    if(!atIsNaN(ps_vec[0])) {
      /* misalignment at entrance */
      if (Elem->T1) ATaddvv(ps_vec, Elem->T1);
      if (Elem->R1) ATmultmv(ps_vec, Elem->R1);

      track_element(ps_vec, Elem);

      /* misalignment at exit */
      if (Elem->R2) ATmultmv(ps_vec, Elem->R2);
      if (Elem->T2) ATaddvv(ps_vec, Elem->T2);
    }
  }
}

//------------------------------------------------------------------------------

void E1rotation(double ps[],double X0ref, double E1)
/* At Entrance Edge:
   move particles to the field edge and convert coordinates to x, dx/dz, y,
   dy/dz, then convert to x, px, y, py as integration is done with px, py     */
{
  double x0, dxdz0, dydz0, psi, fac;

  dxdz0 = ps[1]/sqrt(sqr(1+ps[4])-sqr(ps[1])-sqr(ps[3]));
  dydz0 = ps[3]/sqrt(sqr(1+ps[4])-sqr(ps[1])-sqr(ps[3]));
  x0 = ps[0];

  psi = atan(dxdz0);
  ps[0] = ps[0]*cos(psi)/cos(E1+psi)+X0ref;
  ps[1] = tan(E1+psi);
  ps[3] = dydz0/(cos(E1)-dxdz0*sin(E1));
  ps[2] += x0*sin(E1)*ps[3];
  ps[5] += x0*tan(E1)/(1-dxdz0*tan(E1))*sqrt(1+sqr(dxdz0)+sqr(dydz0));
  /* convert to px, py */
  fac = sqrt(1+sqr(ps[1])+sqr(ps[3]));
  ps[1] = ps[1]*(1+ps[4])/fac;
  ps[3] = ps[3]*(1+ps[4])/fac;
}

void E2rotation(double ps[], double X0ref, double E2)
/* At Exit Edge:
   move particles to arc edge and convert coordinates to x, px, y, py         */
{
  double x0, dxdz0, dydz0, psi, fac;

  dxdz0 = ps[1]/sqrt(sqr(1+ps[4])-sqr(ps[1])-sqr(ps[3]));
  dydz0 = ps[3]/sqrt(sqr(1+ps[4])-sqr(ps[1])-sqr(ps[3]));
  x0 = ps[0];

  psi = atan(dxdz0);
  fac = sqrt(1+sqr(dxdz0)+sqr(dydz0));
  ps[0] = (ps[0]-X0ref)*cos(psi)/cos(E2+psi);
  ps[1] = tan(E2+psi);
  ps[3] = dydz0/(cos(E2)-dxdz0*sin(E2));
  ps[2] += ps[3]*(x0-X0ref)*sin(E2);
  ps[5] += (x0-X0ref)*tan(E2)/(1-dxdz0*tan(E2))*fac;
  /* convert to px, py */
  fac = sqrt(1+sqr(ps[1])+sqr(ps[3]));
  ps[1] = ps[1]*(1+ps[4])/fac;
  ps[3] = ps[3]*(1+ps[4])/fac;
}

void edgey(double ps[], double inv_rho, double edge_angle)
/* Edge focusing in dipoles with hard-edge field for vertical only */
{
  double psi = inv_rho*tan(edge_angle);

  /*ps[1]+=ps[0]*psi;*/
  ps[3]-=ps[2]*psi;
}

void edgey_fringe(double ps[], double inv_rho, double edge_angle, double fint,
		  double gap)
/* Edge focusing in dipoles with fringe field, for vertical only */
{
  double
    // fx = inv_rho*tan(edge_angle),
    psi_bar =
    edge_angle-inv_rho*gap*fint*(1+sin(edge_angle)*sin(edge_angle))
    /cos(edge_angle)/(1+ps[4]),
    fy = inv_rho*tan(psi_bar);

  /*ps[1]+=ps[0]*fx;*/
  ps[3]-=ps[2]*fy;
}

void CBendPass(double ps[], const int num_particles, elem_type *Elem)
{
  int    k, m;
  double *ps_vec, SL, L1, L2, K1, K2;
  bool   useT1, useT2, useR1, useR2, useFringe1, useFringe2;

  const elem_mpole *mpole = Elem->mpole_ptr;

  SL = Elem->Length/Elem->mpole_ptr->NumIntSteps;
  L1 = SL*DRIFT1;
  L2 = SL*DRIFT2;
  K1 = SL*KICK1;
  K2 = SL*KICK2;

  /* mexPrintf("E0ref=%f\n",X0ref); */
  useT1 = (Elem->T1 != NULL);
  useT2 = (Elem->T2 != NULL);
  useR1 = (Elem->R1 != NULL);
  useR2 = (Elem->R2 != NULL);
  /* if either is 0 - do not calculate fringe effects */
  useFringe1 = (mpole->FringeInt1 != 0 && mpole->FullGap != 0);
  useFringe2 = (mpole->FringeInt2 != 0 && mpole->FullGap != 0);

  for(k = 0; k < num_particles; k++)	{   /* Loop over particles  */
    ps_vec = ps+k*PS_DIM;
    if(!atIsNaN(ps_vec[0])) {
      /*  misalignment at entrance  */
      if(useT1)	ATaddvv(ps_vec, Elem->T1);
      if(useR1)	ATmultmv(ps_vec, Elem->R1);
      /* edge focus */
      if(useFringe1)
	edgey_fringe(ps_vec, mpole->irho+mpole->PolynomB[1]*mpole->X0ref,
		     mpole->EntranceAngle, mpole->FringeInt1, mpole->FullGap);
      else
	edgey(ps_vec, mpole->irho+mpole->PolynomB[1]*mpole->X0ref,
	      mpole->EntranceAngle);
      /* Rotate and translate to straight Cartesian coordinate */
      E1rotation(ps_vec, mpole->X0ref, mpole->EntranceAngle);
      /* integrator */
      for(m = 0; m < mpole->NumIntSteps; m++) {
	/* Loop over slices */
	ps_vec = ps+k*PS_DIM;
	Drift(ps_vec, L1, true);
	// Irho = 0.
	thin_kick(ps_vec, mpole->PolynomA, mpole->PolynomB, K1, mpole->irho,
		  mpole->MaxOrder, mpole->Energy, false);

	Drift(ps_vec, L2, true);
	thin_kick(ps_vec, mpole->PolynomA, mpole->PolynomB, K2, mpole->irho,
		  mpole->MaxOrder, mpole->Energy, false);
	Drift(ps_vec, L2, true);
	thin_kick(ps_vec, mpole->PolynomA, mpole->PolynomB, K1, mpole->irho,
		  mpole->MaxOrder, mpole->Energy, false);
	Drift(ps_vec, L1, true);
      }
      /* Rotate and translate back to curvilinear coordinate */
      E2rotation(ps_vec, mpole->X0ref, mpole->ExitAngle);
      ps_vec[5] -= mpole->RefDZ;
      if(useFringe2)
	edgey_fringe(ps_vec, mpole->irho
		     +mpole->PolynomB[1]*mpole->X0ref,
		     mpole->ExitAngle, mpole->FringeInt2,
		     mpole->FullGap);
      else    /* edge focus */
	edgey(ps_vec,
	      mpole->irho
	      +mpole->PolynomB[1]*mpole->X0ref,
	      mpole->ExitAngle);
      /* Misalignment at exit */
      if(useR2)	ATmultmv(ps_vec, Elem->R2);
      if(useT2)	ATaddvv(ps_vec, Elem->T2);
    }
  }
}

//------------------------------------------------------------------------------

/* This code was modified from the original BndMPoleSymplectic4Pass.c of AT
   to correctly integrate the Hamiltonian in  the curvilinear coordinate
   system of the dipole and to include the second order Transport map of the
   fringe field. New version created by Xiaobiao Huang in March 2009, in
   final verified version in August 2009.                                     */

static void ATbendhxdrift6(double* r, double L,double h)
/* the pseudo-drift element described by Hamiltonian
   H1 = (1+hx) (px^2+py^2)/2(1+delta),                                        */
{
  double hs = h*L;
  double i1pd = 1.0/(1+r[4]);
  double x=r[0],px=r[1],py=r[3];

  r[0] += (1+h*x)*px*i1pd*L+1/4.*hs*L*(sqr(px)-sqr(py))*i1pd*i1pd;
  /* (1.0/h+x)*((1.0+hs*px*i1pd/2.)*(1.0+hs*px*i1pd/2.)-(hs*py*i1pd/2.)
     *(hs*py*i1pd/2.))-1./h;*/
  r[1] -= hs*(sqr(px)+sqr(py))*i1pd/2.0;

  r[2]+= (1.0+h*x)*i1pd*py*L*(1.+px*hs/2.0);
  r[5]+= (1.0+h*x)*i1pd*i1pd*L/2.0*(sqr(px)+sqr(py));
}

static void edge_fringe2A(double* r, double inv_rho, double edge_angle,
			  double fint, double gap,double h1,double K1)
{
  /* Entrance Fringe field transport map to second order in dipoles with
     fringe field                                                             */
  double fx = inv_rho*tan(edge_angle);
  double dpsi =
    inv_rho*gap*fint*(1+sin(edge_angle)*sin(edge_angle))/cos(edge_angle);
  /* /(1+r[4]); */
  double psi_bar = edge_angle-dpsi;
  double fy = inv_rho*tan(psi_bar);
  double h = inv_rho;
  double tpsi=tan(edge_angle), tpsib=tan(psi_bar);
  double spsi=1.0/cos(edge_angle);
  /* spsib=1.0/cos(psi_bar) */
  double T111,T234,T414,   T212,T313,  T133,T423,T211,T233,T413;
  double r0 = r[0], r2 = r[2], r1 = r[1];

  T111 = -0.5*h*tpsi*tpsi;
  /*  T234=  -0.5*h*tpsi*tpsib;    */
  T234 = -0.5*h*tpsi*tpsi;
  T414 = T234;
  T212 = -T111;
  T313 = -T234;
  T133 = 0.5*h*spsi*spsi;  T423=-T133;
  T211 = 0.5*h*h1*spsi*spsi*spsi + K1*tpsi;
  T233 = -0.5*h*h1*spsi*spsi*spsi -K1*tpsi+0.5*h*h*tpsi*(tpsib*tpsib+spsi*spsi);
  T413 = -0.5*h*h1*spsi*spsi*spsi -K1*tpsi;
  /*-0.5*h*h*tpsi*(spsi*spsi+tpsib*tpsib);*/

  r[0] += T111*r[0]*r[0]+T133*r[2]*r[2];
  r[1] += r0*fx + 2*T212*r0*r[1]+2*T234*r[2]*r[3]+T211*r0*r0+T233*r[2]*r[2] ;
  r[2] += 2*T313*r0*r[2];
  r[3] += -r2*fy + 2*T414*r0*r[3]+2*T413*r0*r2+2*T423*r1*r2 ;

}

static void edge_fringe2B(double* r, double inv_rho, double edge_angle,
			  double fint, double gap,double h2,double K1)
{
  /* Exit Fringe field transport map to second order in dipoles with fringe
     field                                                                    */
    double
      fx = inv_rho*tan(edge_angle),
      dpsi =
      inv_rho*gap*fint*(1+sin(edge_angle)*sin(edge_angle))/cos(edge_angle),
      /* /(1+r[4]),  */
      psi_bar = edge_angle-dpsi,
      fy      = inv_rho*tan(psi_bar),
      h       = inv_rho,
      tpsi    = tan(edge_angle), tpsib=tan(psi_bar),
      spsi    = 1.0/cos(edge_angle), /* spsib=1.0/cos(psi_bar) */
      r0      = r[0], r2 = r[2], r1 = r[1],
      T111,
      T234,
      T414,
      T212,
      T313,
      T133,
      T423,
      T211,
      T233,
      T413;

    T111 = 0.5*h*tpsi*tpsi;
    /*  T234=  0.5*h*tpsi*tpsib;    */
    T234 = 0.5*h*tpsi*tpsi;
    T414 = T234;
    T212 = -T111;
    T313 = -T234;
    T133 = -0.5*h*spsi*spsi;  T423=-T133;
    T211 = 0.5*h*h2*spsi*spsi*spsi +K1*tpsi-0.5*h*h*tpsi*tpsi*tpsi;
    T233 = -0.5*h*h2*spsi*spsi*spsi -K1*tpsi-0.5*h*h*tpsi*tpsib*tpsib;
    T413 = -0.5*h*h2*spsi*spsi*spsi -K1*tpsi+0.5*h*h*tpsi*(spsi*spsi);

    r[0] += T111*r[0]*r[0]+T133*r[2]*r[2];
    r[1] += r0*fx + 2*T212*r0*r[1]+2*T234*r[2]*r[3]+T211*r0*r0+T233*r[2]*r[2] ;
    r[2] += 2*T313*r0*r[2];
    r[3] += -r2*fy + 2*T414*r0*r[3]+2*T413*r0*r2+2*T423*r1*r2 ;

}

void MpoleE2Pass(double ps[], const int num_particles, elem_type *Elem)
{
  bool   useT1, useT2, useR1, useR2, useFringe1, useFringe2;
  int    k, m;
  double *ps_vec;

  const elem_mpole *mpole = Elem->mpole_ptr;

  const double
    SL = Elem->Length/mpole->NumIntSteps,
    L1 = SL*DRIFT1,
    L2 = SL*DRIFT2,
    K1 = SL*KICK1,
    K2 = SL*KICK2;

  useT1 = (Elem->T1 != NULL);
  useT2 = (Elem->T2 != NULL);
  useR1 = (Elem->R1 != NULL);
  useR2 = (Elem->R2 != NULL);

  /* if either is 0 - do not calculate fringe effects */
  useFringe1 = (mpole->FringeInt1 != 0 && mpole->FullGap != 0);
  useFringe2 = (mpole->FringeInt2 != 0 && mpole->FullGap != 0);

#pragma omp parallel for if (num_particles > OMP_PARTICLE_THRESHOLD)    \
  default(shared) shared(r, num_particles) private(c, ps_vec, m)

  for(k = 0; k < num_particles; k++) { /* Loop over particles */
    ps_vec = ps+k*PS_DIM;
    if(!atIsNaN(ps_vec[0])) {
      /*  misalignment at entrance  */
      if(useT1) ATaddvv(ps_vec, Elem->T1);
      if(useR1) ATmultmv(ps_vec, Elem->R1);
      /* Check physical apertures at the entrance of the magnet */
      if (Elem->RApertures) checkiflostRectangularAp(ps_vec, Elem->RApertures);
      if (Elem->EApertures) checkiflostEllipticalAp(ps_vec, Elem->EApertures);
      /* edge focus */
      if(useFringe1) {
	edge_fringe2A(ps_vec, mpole->irho, mpole->EntranceAngle,
		      mpole->FringeInt1, mpole->FullGap, mpole->H1,
		      mpole->PolynomB[1]);
      } else {
	edge_fringe2A(ps_vec, mpole->irho, mpole->EntranceAngle, 0, 0,
		      mpole->H1, mpole->PolynomB[1]);
      }
      /* integrator */
      for(m = 0; m < mpole->NumIntSteps; m++) { /* Loop over slices*/
	ps_vec = ps+k*PS_DIM;
	ATbendhxdrift6(ps_vec, L1, mpole->irho);
	thin_kick(ps_vec, mpole->PolynomA, mpole->PolynomB, K1, mpole->irho,
		  mpole->MaxOrder, mpole->Energy, false);
	ATbendhxdrift6(ps_vec,L2, mpole->irho);
	thin_kick(ps_vec, mpole->PolynomA, mpole->PolynomB, K2, mpole->irho,
		  mpole->MaxOrder, mpole->Energy, false);
	ATbendhxdrift6(ps_vec,L2, mpole->irho);
	thin_kick(ps_vec, mpole->PolynomA, mpole->PolynomB, K1, mpole->irho,
		  mpole->MaxOrder, mpole->Energy, false);
	ATbendhxdrift6(ps_vec,L1, mpole->irho);
      }
      /* edge focus */
      if(useFringe2) {
	edge_fringe2B(ps_vec, mpole->irho, mpole->ExitAngle, mpole->FringeInt2,
		      mpole->FullGap, mpole->H2, mpole->PolynomB[1]);
      } else {
	edge_fringe2B(ps_vec, mpole->irho, mpole->ExitAngle, 0, 0, mpole->H2,
		      mpole->PolynomB[1]);
      }
      /* Check physical apertures at the exit of the magnet */
      if (Elem->RApertures) checkiflostRectangularAp(ps_vec, Elem->RApertures);
      if (Elem->EApertures) checkiflostEllipticalAp(ps_vec, Elem->EApertures);
      /* Misalignment at exit */
      if(useR2) ATmultmv(ps_vec, Elem->R2);
      if(useT2) ATaddvv(ps_vec, Elem->T2);
    }
  }
}

//------------------------------------------------------------------------------

/*----------------------------------------------------------------------------
  Modification Log:
  -----------------
  .04  2003-04-29      YK Wu, Duke University, wu@fel.duke.edu
 		       using scientific notation for constants.
                       Checked with TRACY pascal code.
                       Computing differential pathlength only.

  .03  2003-04-28      YK Wu, Duke University, wu@fel.duke.edu
 		       Convert to C code and cross-checked with the pascal
		       version;

  .02  2001-12-xx      Y. K. Wu, Duke University, wu@fel.duke.edu
                       Implementing DA version of the wiggler integrator for
		       Pascal.
                       Gauge is disabled !!! (Dec. 4, 2001)

  .01  2001-02-12      Y. K. Wu, LBNL
                       Implementing a generic wiggler integrator
                       for paraxial-ray Hamiltonian approximation.

 *----------------------------------------------------------------------------
   Accelerator Physics Group, Duke FEL Lab, www.fel.duke.edu

                                                                              */
double sinc(double x)
{
  double x2, result;
  /* Expand sinc(x) = sin(x)/x to x^8 */
  x2 = x*x;
  result = 1e0 - x2/6e0*(1e0 - x2/20e0 *(1e0 - x2/42e0*(1e0-x2/72e0) ) );
  return result;
}

void GWigAx(struct elem_wig *pWig, double *Xvec, double *pax, double *paxpy)
{
  int    i;
  double x, y, z;
  double kx, ky, kz, tz, kw;
  double cx, sxkx, chx, shx;
  double cy, sy, chy, shy, sz;
  double gamma0, beta0;
  double ax, axpy;

  x = Xvec[0];
  y = Xvec[2];
  z = pWig->Zw;

  kw   = 2e0*PI/(pWig->Lw);
  ax   = 0e0;
  axpy = 0e0;

  gamma0   = pWig->E0/XMC2;
  beta0    = sqrt(1e0 - 1e0/(gamma0*gamma0));
  pWig->Aw = (q_e/m_e/clight)/(2e0*PI) * (pWig->Lw) * (pWig->PB0);

  /* Horizontal Wiggler: note that one potentially could have: kx=0 */
  for (i = 0; i < pWig->NHharm; i++) {
    pWig->HCw[i] = pWig->HCw_raw[i]*(pWig->Aw)/(gamma0*beta0);
    kx = pWig->Hkx[i];
    ky = pWig->Hky[i];
    kz = pWig->Hkz[i];
    tz = pWig->Htz[i];

    cx  = cos(kx*x);
    chy = cosh(ky * y);
    sz  = sin(kz * z + tz);
    ax  = ax + (pWig->HCw[i])*(kw/kz)*cx*chy*sz;

    shy = sinh(ky * y);
    if (fabs(kx/kw) > GWIG_EPS) {
      sxkx = sin(kx * x)/kx;
    } else {
      sxkx = x*sinc(kx*x);
    }

    axpy = axpy + pWig->HCw[i]*(kw/kz)*ky*sxkx*shy*sz;
  }

  /* Vertical Wiggler: note that one potentially could have: ky=0 */
  for (i = 0; i < pWig->NVharm; i++ ) {
    pWig->VCw[i] = pWig->VCw_raw[i]*(pWig->Aw)/(gamma0*beta0);
    kx = pWig->Vkx[i];
    ky = pWig->Vky[i];
    kz = pWig->Vkz[i];
    tz = pWig->Vtz[i];

    shx = sinh(kx * x);
    sy  = sin(ky * y);
    sz  = sin(kz * z + tz);
    ax  = ax + pWig->VCw[i]*(kw/kz)*(ky/kx)*shx*sy*sz;

    chx = cosh(kx * x);
    cy  = cos(ky * y);
    axpy = axpy + pWig->VCw[i]*(kw/kz)* pow(ky/kx,2) *chx*cy*sz;
  }

  *pax   = ax;
  *paxpy = axpy;
}

void GWigAy(struct elem_wig *pWig, double *Xvec, double *pay, double *paypx)
{
  int    i;
  double x, y, z;
  double kx, ky, kz, tz, kw;
  double cx, sx, chx, shx;
  double cy, syky, chy, shy, sz;
  double gamma0, beta0;
  double ay, aypx;

  x = Xvec[0];
  y = Xvec[2];
  z = pWig->Zw;

  kw   = 2e0*PI/(pWig->Lw);
  ay   = 0e0;
  aypx = 0e0;

  gamma0  = pWig->E0/XMC2;
  beta0   = sqrt(1e0 - 1e0/(gamma0*gamma0));
  pWig->Aw = (q_e/m_e/clight)/(2e0*PI) * (pWig->Lw) * (pWig->PB0);

  /* Horizontal Wiggler: note that one potentially could have: kx=0 */
  for ( i = 0; i < pWig->NHharm; i++ ){
    pWig->HCw[i] = (pWig->HCw_raw[i])*(pWig->Aw)/(gamma0*beta0);
    kx = pWig->Hkx[i];
    ky = pWig->Hky[i];
    kz = pWig->Hkz[i];
    tz = pWig->Htz[i];

    sx = sin(kx * x);
    shy = sinh(ky * y);
    sz  = sin(kz * z + tz);
    ay  = ay + (pWig->HCw[i])*(kw/kz)*(kx/ky)*sx*shy*sz;

    cx  = cos(kx * x);
    chy = cosh(ky * y);

    aypx = aypx + (pWig->HCw[i])*(kw/kz)*pow(kx/ky,2) * cx*chy*sz;
  }

  /* Vertical Wiggler: note that one potentially could have: ky=0 */
  for (i = 0; i < pWig->NVharm; i++) {
    pWig->VCw[i] = (pWig->VCw_raw[i])*(pWig->Aw)/(gamma0*beta0);
    kx = pWig->Vkx[i];
    ky = pWig->Vky[i];
    kz = pWig->Vkz[i];
    tz = pWig->Vtz[i];

    chx = cosh(kx * x);
    cy  = cos(ky * y);
    sz  = sin(kz * z + tz);
    ay  = ay + (pWig->VCw[i])*(kw/kz)*chx*cy*sz;

    shx = sinh(kx * x);
    if (fabs(ky/kw) > GWIG_EPS) {
      syky  = sin(ky * y)/ky;
    } else {
      syky = y * sinc(ky * y);
    }
    aypx = aypx + (pWig->VCw[i])*(kw/kz)* kx*shx*syky*sz;
  }

  *pay = ay;
  *paypx = aypx;
}

/* This function appears to be unused. */
void GWigGauge(struct elem_wig *pWig, double *X, int flag)
{
  double ax, ay, axpy, aypx;

  GWigAx(pWig, X, &ax, &axpy);
  GWigAy(pWig, X, &ay, &aypx);

  if (flag == Elem_Entrance) {
    /* At the entrance of the wiggler */
    X[1] = X[1] + ax;
    X[3] = X[3] + ay;
  } else if (flag == Elem_Exit) {
    /* At the exit of the wiggler */
    X[1] = X[1] - ax;
    X[3] = X[3] - ay;
  } else {
    printf("  GWigGauge: Unknown flag = %i\n", flag);
  }
}

void GWigMap_2nd(struct elem_wig *pWig, double *X, double dl)
{
  double dld, dl2, dl2d, ax, ay, axpy, aypx;

  dld  = dl/(1.0e0 + X[4]);
  dl2  = 0.5e0 * dl;
  dl2d = dl2/(1.0e0 + X[4]);

  /* Step1: increase a half step in z */
  pWig->Zw = pWig->Zw + dl2;

  /* Step2: a half drift in y */
  GWigAy(pWig, X, &ay, &aypx);
  X[1] = X[1] - aypx;
  X[3] = X[3] - ay;

  X[2] = X[2] + dl2d*X[3];
  X[5] = X[5] + 0.5e0*dl2d*(X[3]*X[3])/(1.0e0+X[4]);

  GWigAy(pWig, X, &ay, &aypx);
  X[1] = X[1] + aypx;
  X[3] = X[3] + ay;

  /* Step3: a full drift in x */
  GWigAx(pWig, X, &ax, &axpy);
  X[1] = X[1] - ax;
  X[3] = X[3] - axpy;

  X[0] = X[0] + dld*X[1];
  /* Full path length
     X[5] = X[5] + dl + 0.5e0*dld*(X[1]*X[1])/(1.0e0+X[4]);
  */
  /* Differential path length only */
  X[5] = X[5] + 0.5e0*dld*(X[1]*X[1])/(1.0e0+X[4]);

  GWigAx(pWig, X, &ax, &axpy);
  X[1] = X[1] + ax;
  X[3] = X[3] + axpy;

  /* Step4: a half drift in y */
  GWigAy(pWig, X, &ay, &aypx);
  X[1] = X[1] - aypx;
  X[3] = X[3] - ay;

  X[2] = X[2] + dl2d*X[3];
  X[5] = X[5] + 0.5e0*dl2d*(X[3]*X[3])/(1.0e0+X[4]);

  GWigAy(pWig, X, &ay, &aypx);
  X[1] = X[1] + aypx;
  X[3] = X[3] + ay;

  /* Step5: increase a half step in z */
  pWig->Zw = pWig->Zw + dl2;

}

void GWigPass_2nd(struct elem_type *Elem, double X[])
{
  int      i, Nstep;
  double   dl;
  elem_wig *pWig = Elem->wig_ptr;

  Nstep = pWig->PN*(pWig->Nw);
  dl    = pWig->Lw/(pWig->PN);

  for (i = 1; i <= Nstep; i++) {
    GWigMap_2nd(pWig, X, dl);
  }
}

void GWigB(struct elem_wig *pWig, double *Xvec, double *B)
/* Compute magnetic field at particle location.
 * Added by M. Borland, August 2007.
 */
{
  int    i;
  double x, y, z;
  double kx, ky, kz, tz, kw;
  double cx, sx, chx, shx;
  double cy, sy, chy, shy;
  double cz;
  /* B0 is a reserved symbol on MacOS, defined in termios.h */
  double _B0;

  x = Xvec[0];
  y = Xvec[2];
  z = pWig->Zw;

  kw   = 2e0*PI/(pWig->Lw);

  B[0] = 0;
  B[1] = 0;

  if (pWig->NHharm && z>=pWig->zStartH && z<=pWig->zEndH) {
    _B0 = pWig->PB0;
    if (!pWig->HSplitPole) {
      /* Normal Horizontal Wiggler: note that one potentially could have:
	 kx=0                                                             */
      for (i = 0; i < pWig->NHharm; i++) {
        kx = pWig->Hkx[i];
        ky = pWig->Hky[i];
        kz = pWig->Hkz[i];
        tz = pWig->Htz[i];

        sx  = sin(kx*x);
        cx  = cos(kx*x);
        chy = cosh(ky * y);
        shy = sinh(ky * y);
        cz = cos(kz*z+tz);

        /* Accumulate field values in user-supplied array (Bx, By) */
        B[0] += _B0*pWig->HCw_raw[i]*kx/ky*sx*shy*cz;
        B[1] -= _B0*pWig->HCw_raw[i]*cx*chy*cz;
      }
    } else {
      /* Split-pole Horizontal Wiggler: note that one potentially could have:
	 ky=0 (caught in main routine)                                        */
      for (i = 0; i < pWig->NHharm; i++) {
        kx = pWig->Hkx[i];
        ky = pWig->Hky[i];
        kz = pWig->Hkz[i];
        tz = pWig->Htz[i];

        shx = sinh(kx*x);
        chx = cosh(kx*x);
        cy  = cos(ky * y);
        sy  = sin(ky * y);
        cz  = cos(kz*z+tz);

        B[0] -= _B0*pWig->HCw_raw[i]*kx/ky*shx*sy*cz;
        B[1] -= _B0*pWig->HCw_raw[i]*chx*cy*cz;
      }
    }
  }

  if (pWig->NVharm && z>=pWig->zStartV && z<=pWig->zEndV) {
    _B0 = pWig->PB0;
    if (!pWig->VSplitPole) {
      /* Normal Vertical Wiggler: note that one potentially could have: ky=0 */
      for (i = 0; i < pWig->NVharm; i++ ) {
        kx = pWig->Vkx[i];
        ky = pWig->Vky[i];
        kz = pWig->Vkz[i];
        tz = pWig->Vtz[i];

        shx = sinh(kx * x);
        chx = cosh(kx * x);
        sy  = sin(ky * y);
        cy  = cos(ky * y);
        cz  = cos(kz*z + tz);

        /* Accumulate field values in user-supplied array (Bx, By) */
        B[0] += _B0*pWig->VCw_raw[i]*chx*cy*cz;
        B[1] -= _B0*pWig->VCw_raw[i]*ky/kx*shx*sy*cz;
      }
    } else {
      /* Split-pole Vertical Wiggler: note that one potentially could have:
	 kx=0 (caught in main routine)                                        */
      for (i = 0; i < pWig->NVharm; i++ ) {
        kx = pWig->Vkx[i];
        ky = pWig->Vky[i];
        kz = pWig->Vkz[i];
        tz = pWig->Vtz[i];

        sx  = sin(kx * x);
        cx  = cos(kx * x);
        shy = sinh(ky * y);
        chy = cosh(ky * y);
        cz  = cos(kz*z + tz);

        /* Accumulate field values in user-supplied array (Bx, By) */
        B[0] += _B0*pWig->VCw_raw[i]*cx*chy*cz;
        B[1] += _B0*pWig->VCw_raw[i]*ky/kx*sx*shy*cz;
      }
    }
  }
}

void GWigRadiationKicks(struct elem_wig *pWig, double *X, double *Bxy,
			double dl)
/* Apply kicks for synchrotron radiation.
   Added by M. Borland, August 2007.                                          */
{
  double irho2, H, dFactor;
  double B2;
  double dDelta;

  /* B^2 in T^2 */
  B2 = (Bxy[0]*Bxy[0]) + (Bxy[1]*Bxy[1]);
  if (B2==0)
    return;

  /* Beam rigidity in T*m */
  H = (pWig->Po)/586.679074042074490;

  /* 1/rho^2 */
  irho2 = B2/(H*H);

  /* (1+delta)^2 */
  dFactor = ((1+X[4])*(1+X[4]));

  /* Classical radiation loss */
  dDelta = -(pWig->srCoef)*dFactor*irho2*dl;
  X[4] += dDelta;
  X[1] *= (1+dDelta);
  X[3] *= (1+dDelta);
}

void GWigPass_2nd(struct elem_type *Elem, double X[], const bool rad)
{
  int      i, Nstep;
  double   dl, B[2], ax, ay, axpy, aypx;
  elem_wig *pWig = Elem->wig_ptr;

  Nstep = pWig->PN*(pWig->Nw);
  dl    = pWig->Lw/(pWig->PN);

  if (rad) {
    GWigAx(pWig, X, &ax, &axpy);
    GWigAy(pWig, X, &ay, &aypx);
    GWigB(pWig, X, B);
    X[1] -= ax;
    X[3] -= ay;
    GWigRadiationKicks(pWig, X, B, dl);
    X[1] += ax;
    X[3] += ay;
  }

  for (i = 1; i <= Nstep; i++) {
    GWigMap_2nd(pWig, X, dl);
    if (rad) {
      GWigAx(pWig, X, &ax, &axpy);
      GWigAy(pWig, X, &ay, &aypx);
      GWigB(pWig, X, B);
      X[1] -= ax;
      X[3] -= ay;
      GWigRadiationKicks(pWig, X, B, dl);
      X[1] += ax;
      X[3] += ay;
    }
  }
}

void GWigPass_4th(struct elem_type *Elem, double X[], const bool rad)
{
  int      i, Nstep;
  double   dl, dl1, dl0, B[2], ax, ay, axpy, aypx;
  elem_wig *pWig = Elem->wig_ptr;

  const double
    x1 =  1.3512071919596576340476878089715e0,
    x0 = -1.7024143839193152680953756179429e0;

  Nstep = pWig->PN*(pWig->Nw);
  dl = pWig->Lw/(pWig->PN);

  dl1 = x1*dl;
  dl0 = x0*dl;

  if (rad) {
    GWigAx(pWig, X, &ax, &axpy);
    GWigAy(pWig, X, &ay, &aypx);
    GWigB(pWig, X, B);
    X[1] -= ax;
    X[3] -= ay;
    GWigRadiationKicks(pWig, X, B, dl);
    X[1] += ax;
    X[3] += ay;
  }

  for (i = 1; i <= Nstep; i++ ) {
    GWigMap_2nd(pWig, X, dl1);
    GWigMap_2nd(pWig, X, dl0);
    GWigMap_2nd(pWig, X, dl1);

    if (rad) {
      GWigAx(pWig, X, &ax, &axpy);
      GWigAy(pWig, X, &ay, &aypx);
      GWigB(pWig, X, B);
      X[1] -= ax;
      X[3] -= ay;
      GWigRadiationKicks(pWig, X, B, dl);
      X[1] += ax;
      X[3] += ay;
    }
  }
}

void WigPass(double ps[], const int num_particles, struct elem_type *Elem)
{
  int    k;
  double *ps_vec;

  for (k = 0; k < num_particles; k++) {
    ps_vec = ps+k*PS_DIM;
    if (!atIsNaN(ps_vec[0])) {
      switch (Elem->wig_ptr->Pmethod) {
      case second:
	GWigPass_2nd(Elem, ps_vec, false);
	break;
      case fourth:
	GWigPass_4th(Elem, ps_vec, false);
	break;
      default:
	printf("Invalid wiggler integration method %d.\n",
	       Elem->wig_ptr->Pmethod);
	break;
      }
    }
  }
}

void WigRadPass(double ps[], const int num_particles, struct elem_type *Elem)
{
  int    k;
  double *ps_vec;

  for(k = 0; k < num_particles; k++) {
    ps_vec = ps+k*PS_DIM;
    if(!atIsNaN(ps_vec[0])) {
      switch (Elem->wig_ptr->Pmethod) {
      case second:
	GWigPass_2nd(Elem, ps_vec, true);
	break;
      case fourth:
	GWigPass_4th(Elem, ps_vec, true);
	break;
      default:
	printf("Invalid wiggler integration method %d.\n",
	       Elem->wig_ptr->Pmethod);
	break;
      }
    }
  }
}

//------------------------------------------------------------------------------
