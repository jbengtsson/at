#include <math.h>

#include <armadillo>

#include "tracy-2.h"

#include "gwig.cc"
#include "track.cc"


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

struct elem_type* init_elem(const PyObject *ElemData, struct elem_type *Elem,
			    const bool len, const bool aper)
{
  double Length, *R1, *R2, *T1, *T2, *EApertures, *RApertures;

  Elem = (struct elem_type*)malloc(sizeof(struct elem_type));

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

  Elem = (struct elem_type*)malloc(sizeof(struct elem_type));
  Elem->Length = 0e0;
  if (Elem) {
    Elem->ap_ptr = (struct elem_ap*)malloc(sizeof(struct elem_ap));

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
    Elem->drift_ptr = (struct elem_drift*)malloc(sizeof(struct elem_drift));
    return Elem;
  } else
    return NULL;
}

struct elem_type*
init_mpole(const PyObject *ElemData, struct elem_type *Elem, const bool bend,
	   const bool cbend, const bool incl_E0, const bool incl_E2)
{
  int
    MaxOrder, NumIntSteps,  FringeBendEntrance, FringeBendExit,
    FringeQuadEntrance, FringeQuadExit;
  double
    BendingAngle, EntranceAngle, ExitAngle, FullGap, FringeInt1, FringeInt2,
    *PolynomA, *PolynomB, *fringeIntM0, *fringeIntP0, *KickAngle, X0ref,
    ByError, RefDZ, Energy, h1, h2;
  elem_mpole
    *mpole;

  Elem = init_elem(ElemData, Elem, true, !cbend);
  if (Elem) {
    Elem->mpole_ptr = (struct elem_mpole*)malloc(sizeof(struct elem_mpole));
    mpole           = Elem->mpole_ptr;

    if (incl_E0) {
      Energy = atGetDouble(ElemData,                 (char*)"Energy");
      check_error();
    }

    PolynomA = atGetDoubleArray(ElemData,            (char*)"PolynomA");
    check_error();
    PolynomB = atGetDoubleArray(ElemData,            (char*)"PolynomB");
    check_error();
    MaxOrder = atGetLong(ElemData,                   (char*)"MaxOrder");
    check_error();

    NumIntSteps = atGetLong(ElemData,                (char*)"NumIntSteps");
    check_error();

    if (bend) {
      BendingAngle = atGetDouble(ElemData,           (char*)"BendingAngle");
      check_error();
      EntranceAngle = atGetDouble(ElemData,          (char*)"EntranceAngle");
      check_error();
      ExitAngle = atGetDouble(ElemData,              (char*)"ExitAngle");
      check_error();

      FringeBendEntrance =
	atGetOptionalLong(ElemData,                  (char*)"FringeBendEntrance",
			  1);
      check_error();
      FringeBendExit = atGetOptionalLong(ElemData,   (char*)"FringeBendExit", 1);
      check_error();

      FullGap = atGetOptionalDouble(ElemData,        (char*)"FullGap", 0);
      check_error();
      FringeInt1 = atGetOptionalDouble(ElemData,     (char*)"FringeInt1", 0);
      check_error();
      FringeInt2 = atGetOptionalDouble(ElemData,     (char*)"FringeInt2", 0);
      check_error();
    }

    if (cbend) {
      X0ref = atGetOptionalDouble(ElemData,          (char*)"X0ref", 0);
      check_error();
      ByError = atGetOptionalDouble(ElemData,        (char*)"ByError", 0);
      check_error();
      RefDZ = atGetOptionalDouble(ElemData,          (char*)"RefDZ", 0);
      check_error();
    }

    if (incl_E2) {
      h1 = atGetOptionalDouble(ElemData,              (char*)"H1", 0);
      check_error();
      h2 = atGetOptionalDouble(ElemData,              (char*)"H2", 0);
      check_error();
    }

    FringeQuadEntrance =
      atGetOptionalLong(ElemData,                    (char*)"FringeQuadEntrance",
			0);
    check_error();
    FringeQuadExit = atGetOptionalLong(ElemData,     (char*)"FringeQuadExit",
				       0);
    check_error();
    fringeIntM0 = atGetOptionalDoubleArray(ElemData, (char*)"fringeIntM0");
    check_error();
    fringeIntP0 = atGetOptionalDoubleArray(ElemData, (char*)"fringeIntP0");
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

#if 1

//------------------------------------------------------------------------------

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

void thin_kick(const int Order, const double MB[], const double L,
	       const double h_bend, const double h_ref, std::vector<double> &ps)
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

    // if (false) {
    //   B[X_] = BxoBrho; B[Y_] = ByoBrho + h_bend; B[Z_] = 0e0;
    //   radiate(ps, L, h_ref, B);
    // }

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

void Drift(double L, std::vector<double> &ps)
{
  double u;

  if (true) {
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

#include <iomanip>

void Cav_Pass(const double L, const double f_RF, const double V_RFoE0,
	      const double phi, std::vector<double> &ps)
{
  double delta;

  Drift(L/2e0, ps);
  if (V_RFoE0 != 0e0) {
    delta = -V_RFoE0*sin(2e0*M_PI*f_RF/C0*ps[ct_]+phi);
    ps[delta_] += delta;

    // if (globval.radiation) globval.dE -= is_double<T>::cst(delta);

    // if (globval.pathlength) ps[ct_] -= C->Ph/C->Pfreq*c0;
  }
  Drift(L/2e0, ps);
  for (int k = 0; k < 6; k++)
    std::cout << std::scientific << std::setprecision(3) << std::setw(11)
	      << ps[k];
  std::cout << "\n";
}

//------------------------------------------------------------------------------

inline void atdrift(double ps[], const double L)
{
  std::vector<double> ps_stl(PS_DIM, 0e0);

  ps_stl = arrtostl(ps);
  Drift(L, ps_stl);
  stltoarr(ps_stl, ps);
}

inline void fastdrift(double ps[], const double L)
{
  std::vector<double> ps_stl = arrtostl(ps);

  Drift(L*(1e0+ps[delta_]), ps_stl);
  stltoarr(ps_stl, ps);
}

void edge_fringe(double ps[], const double inv_rho,
		 const double edge_angle, const double fint,
		 const double gap, const int method, const bool hor)
{
  std::vector<double> ps_stl = arrtostl(ps);

  EdgeFocus(inv_rho, edge_angle*180e0/M_PI, gap, ps_stl);
  stltoarr(ps_stl, ps);
}

void thin_kick(double ps[], const double a[], const double b[],
		 const double L, const double irho, const int n_max)
{
  double              bn[2*HOMmax+1];
  std::vector<double> ps_stl = arrtostl(ps);

  for (int k = n_max+1; k > 0; k--) {
    bn[HOMmax+k] = b[k-1];
    bn[HOMmax-k] = a[k-1];
  }
  thin_kick(n_max+1, bn, L, irho, irho, ps_stl);
  stltoarr(ps_stl, ps);
}

void cav_pass(double ps[], const double L, const double V_RFoE0,
	      const double f_RF, const double lag)
{
  std::vector<double> ps_stl = arrtostl(ps);

#if 0
  const double phi = -lag*2e0*M_PI*f_RF/C0;
#else
  const double phi = -lag*TWOPI*f_RF/C0;
#endif

  Cav_Pass(L, f_RF, V_RFoE0, phi, ps_stl);
  stltoarr(ps_stl, ps);
}

#else

#include "at_pass_obsolete.cc"

#endif

//------------------------------------------------------------------------------

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

//------------------------------------------------------------------------------

void E1rotation(double ps[],double X0ref, double E1)
/* At Entrance Edge:
   move particles to the field edge and convert coordinates to x, dx/dz, y,
   dy/dz, then convert to x, px, y, py as integration is done with px, py     */
{
  double x0, dxdz0, dydz0, psi, fac;

  dxdz0 = ps[1]/sqrt(SQR(1+ps[4])-SQR(ps[1])-SQR(ps[3]));
  dydz0 = ps[3]/sqrt(SQR(1+ps[4])-SQR(ps[1])-SQR(ps[3]));
  x0 = ps[0];

  psi = atan(dxdz0);
  ps[0] = ps[0]*cos(psi)/cos(E1+psi)+X0ref;
  ps[1] = tan(E1+psi);
  ps[3] = dydz0/(cos(E1)-dxdz0*sin(E1));
  ps[2] += x0*sin(E1)*ps[3];
  ps[5] += x0*tan(E1)/(1-dxdz0*tan(E1))*sqrt(1+SQR(dxdz0)+SQR(dydz0));
  /* convert to px, py */
  fac = sqrt(1+SQR(ps[1])+SQR(ps[3]));
  ps[1] = ps[1]*(1+ps[4])/fac;
  ps[3] = ps[3]*(1+ps[4])/fac;
}

void E2rotation(double ps[], double X0ref, double E2)
/* At Exit Edge:
   move particles to arc edge and convert coordinates to x, px, y, py         */
{
  double x0, dxdz0, dydz0, psi, fac;

  dxdz0 = ps[1]/sqrt(SQR(1+ps[4])-SQR(ps[1])-SQR(ps[3]));
  dydz0 = ps[3]/sqrt(SQR(1+ps[4])-SQR(ps[1])-SQR(ps[3]));
  x0 = ps[0];

  psi = atan(dxdz0);
  fac = sqrt(1+SQR(dxdz0)+SQR(dydz0));
  ps[0] = (ps[0]-X0ref)*cos(psi)/cos(E2+psi);
  ps[1] = tan(E2+psi);
  ps[3] = dydz0/(cos(E2)-dxdz0*sin(E2));
  ps[2] += ps[3]*(x0-X0ref)*sin(E2);
  ps[5] += (x0-X0ref)*tan(E2)/(1-dxdz0*tan(E2))*fac;
  /* convert to px, py */
  fac = sqrt(1+SQR(ps[1])+SQR(ps[3]));
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
    fx = inv_rho*tan(edge_angle),
    psi_bar =
    edge_angle-inv_rho*gap*fint*(1+sin(edge_angle)*sin(edge_angle))
    /cos(edge_angle)/(1+ps[4]),
    fy = inv_rho*tan(psi_bar);

  /*ps[1]+=ps[0]*fx;*/
  ps[3]-=ps[2]*fy;
}

void ladrift6(double ps[], double L)
/* large angle drift, X. Huang, 7/31/2018
   Input parameter L is the physical length
   1/(1+delta) normalization is done internally
   Hamiltonian H = (1+\delta)-sqrt{(1+\delta)^2-p_x^2-p_y^2}, change sign for
   $\Delta z$ in AT                                                           */
{
  double
    p_norm = 1./sqrt(SQR(1+ps[4])-SQR(ps[1])-SQR(ps[3])),
    NormL = L*p_norm;

  ps[0]+= NormL*ps[1];
  ps[2]+= NormL*ps[3];
  ps[5]+= L*(p_norm*(1+ps[4])-1.);
}

void bndstrthinkick(double ps[], double* A, double* B, double L, double irho,
		    int max_order)
/*****************************************************************************
  Calculate multipole kick in a straight bending magnet, This is not the usual
  Bends!
  created by X. Huang, 7/31/2018
  The reference coordinate system  is straight in s.
  The B vector does not contain b0, we assume b0=irho

  Note: in the US convention the transverse multipole field is written as:

                         max_order+1
                           ----
                           \                       n-1
	   (B + iB  )/ B rho  =  >   (ia  + b ) (x + iy)
         y    x            /       n    n
	                       ----
                          n=1
	is a polynomial in (x,y) with the highest order = MaxOrder


	Using different index notation

                         max_order
                           ----
                           \                       n
	   (B + iB  )/ B rho  =  >   (iA  + B ) (x + iy)
         y    x            /       n    n
	                       ----
                          n=0

	A,B: i=0 ... max_order
   [0] - dipole, [1] - quadrupole, [2] - sextupole ...
   units for A,B[i] = 1/[m]^(i+1)
	Coeficients are stroed in the PolynomA, PolynomB field of the element
	structure in MATLAB

	A[i] (C++,C) =  PolynomA(i+1) (MATLAB)
	B[i] (C++,C) =  PolynomB(i+1) (MATLAB)
	i = 0 .. MaxOrder
******************************************************************************/
{
  int    i;
  double
    ReSum = B[max_order],
    ImSum = A[max_order],
    ReSumTemp;

  /* recursively calculate the local transvrese magnetic field
     Bx = ReSum, By = ImSum */
  B[0] = irho;
  for(i = max_order-1; i >= 0; i--) {
    ReSumTemp = ReSum*ps[0] - ImSum*ps[2] + B[i];
    ImSum = ImSum*ps[0] +  ReSum*ps[2] + A[i];
    ReSum = ReSumTemp;
  }
  ps[1] -=  L*(ReSum);
  ps[3] +=  L*ImSum;
  ps[5] +=  0; /* pathlength */
}

//------------------------------------------------------------------------------

static void drift6(double* r, double L)
/*   Input parameter L is the physical length
     1/(1+delta) normalization is done internally
*/
{	double p_norm = 1/(1+r[4]); 
	double NormL  = L*p_norm;   
	r[0]+= NormL*r[1]; 
	r[2]+= NormL*r[3];
	r[5]+= NormL*p_norm*(r[1]*r[1]+r[3]*r[3])/2;
}

static double B2perp(double bx, double by, double irho,
		     double x, double xpr, double y, double ypr)
/* Calculates sqr(|e x B|) , where e is a unit vector in the direction of
   velocity                                                                   */
        
{
  double v_norm2 = 1/(SQR(1+x*irho)+ SQR(xpr) + SQR(ypr));
    
  /* components of the  velocity vector
   * double ex, ey, ez;
   * ex = xpr;
   * ey = ypr;
   * ez = (1+x*irho);
   */
    
  return ((SQR(by*(1+x*irho)) + SQR(bx*(1+x*irho))
	   + SQR(bx*ypr - by*xpr) )*v_norm2) ;
}

static void bndthinkickrad(double* r, double* A, double* B, double L,
			   double irho, double E0, int max_order)

/*****************************************************************************
Calculate multipole kick in a curved elemrnt (bending magnet)
The reference coordinate system  has the curvature given by the inverse
(design) radius irho.
IMPORTANT !!!
The magnetic field Bo that provides this curvature MUST NOT be included in the dipole term
PolynomB[1](MATLAB notation)(C: B[0] in this function) of the By field expansion
HOWEVER!!! to calculate the effect of classical radiation the full field must be
used in the square of the |v x B|.
When calling B2perp(Bx, By, ...), use the By = RESum + irho, where ImSum is the sum of
the polynomial terms in PolynomB.

 The kick is given by
 
             e L      L delta      L x
  theta  = - --- B  + -------  -  -----  ,
       x      p   y     rho           2
               0                   rho
 
           e L
  theta  = --- B
       y    p   x
             0
 
 ******************************************************************************/
{
   int i;
   double ImSum = A[max_order];
   double ReSum = B[max_order];
   double ReSumTemp;
   double x ,xpr, y, ypr, p_norm,dp_0, B2P;
   double CRAD = CGAMMA*E0*E0*E0/(TWOPI*1e27);	/* [m]/[GeV^3] M.Sands (4.1) */
   
   /* recursively calculate the local transvrese magnetic field
     Bx = ReSum, By = ImSum                                                   */
   for (i=max_order-1; i>=0; i--) {
   	ReSumTemp = ReSum*r[0] - ImSum*r[2] + B[i];
        ImSum = ImSum*r[0] +  ReSum*r[2] + A[i];
        ReSum = ReSumTemp;
   }
   
   /* calculate angles from momentums 	*/
   p_norm = 1/(1+r[4]);
   x   = r[0];
   xpr = r[1]*p_norm;
   y   = r[2];
   ypr = r[3]*p_norm;
   
   B2P = B2perp(ImSum, ReSum +irho, irho, x , xpr, y ,ypr);
   
   dp_0 = r[4];
   r[4] = r[4] - CRAD*SQR(1+r[4])*B2P*(1 + x*irho + (SQR(xpr)+SQR(ypr))/2 )*L;
   
   /* recalculate momentums from angles after losing energy for radiation */
   p_norm = 1/(1+r[4]);
   r[1] = xpr/p_norm;
   r[3] = ypr/p_norm;
   
   r[1] -=  L*(ReSum-(dp_0-r[0]*irho)*irho);
   r[3] +=  L*ImSum;
   r[5] +=  L*irho*r[0]; /* pathlength */
}

double StrB2perp(double bx, double by, double x, double xpr, double y,
		 double ypr)
/* Calculates sqr(|B x e|) , where e is a unit vector in the direction of
   velocity                                                                   */
{
  double v_norm2;

  v_norm2 = 1/(1 + SQR(xpr) + SQR(ypr));

  /* components of the normalized velocity vector
     double ex, ey, ez;
     ex = xpr; 
     ey = ypr; 
     ez = 1;
  */
  	
  return((SQR(by) + SQR(bx) + SQR(bx*ypr - by*xpr) )*v_norm2) ;
} 

static void strthinkickrad(double ps[], const double A[], const double B[],
			   const double L, const double E0, const int max_order)
/*****************************************************************************
 Calculate and apply a multipole kick to a 6-dimentional
 phase space vector in a straight element ( quadrupole)
 
 IMPORTANT !!!
 he reference coordinate system is straight but the field expansion may still
 ontain dipole terms: PolynomA(1), PolynomB(1) - in MATLAB notation,
 [0], B[0] - C,C++ notation
 
******************************************************************************/
{
  int
    i;
  double
    ReSum = B[max_order],
    ImSum = A[max_order],
    ReSumTemp,
    irho = 0, /*straight elements no curvature.*/
    x, xpr, y, ypr, p_norm,dp_0, B2P,
    CRAD = CGAMMA*E0*E0*E0/(TWOPI*1e27);	/* [m]/[GeV^3] M.Sands (4.1) */
   
  for (i = max_order-1; i >= 0; i--) {
    ReSumTemp = ReSum*ps[0] - ImSum*ps[2] + B[i];
    ImSum = ImSum*ps[0] +  ReSum*ps[2] + A[i];
    ReSum = ReSumTemp;
  }
   
  /* calculate angles from momentums 	*/
  p_norm = 1/(1+ps[4]);
  x   = ps[0];
  xpr = ps[1]*p_norm;
  y   = ps[2];
  ypr = ps[3]*p_norm;
   
  /*B2P = B2perp(ImSum, ReSum +irho, irho, x , xpr, y ,ypr);*/
  B2P = StrB2perp(ImSum, ReSum, x , xpr, y, ypr);
   
  dp_0 = ps[4];
  ps[4] = ps[4] - CRAD*SQR(1+ps[4])*B2P*(1 + x*irho + (SQR(xpr)+SQR(ypr))/2 )*L;
   
  /* recalculate momentums from angles after losing energy for radiation */
  p_norm = 1/(1+ps[4]);
  ps[1] = xpr/p_norm;
  ps[3] = ypr/p_norm;
   
  ps[1] -=  L*(ReSum-(dp_0-ps[0]*irho)*irho);
  ps[3] +=  L*ImSum;
  ps[5] +=  L*irho*ps[0]; /* pathlength */
}

//------------------------------------------------------------------------------

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
    
  double r0=r[0],r2=r[2],r1=r[1];
  T111 = -0.5*h*tpsi*tpsi;
  /*  T234=  -0.5*h*tpsi*tpsib;    */
  T234=  -0.5*h*tpsi*tpsi;
  T414=T234;
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

static void edge_fringe2B(double* r, double inv_rho, double edge_angle, double fint, double gap,double h2,double K1)
{   /* Exit Fringe field transport map to second order in dipoles with fringe field */
    double fx = inv_rho*tan(edge_angle);
    double dpsi = inv_rho*gap*fint*(1+sin(edge_angle)*sin(edge_angle))/cos(edge_angle); /* /(1+r[4]);  */
    double psi_bar = edge_angle-dpsi;
    double fy = inv_rho*tan(psi_bar);
    double h = inv_rho;
    double tpsi=tan(edge_angle), tpsib=tan(psi_bar);
    double spsi=1.0/cos(edge_angle); /* spsib=1.0/cos(psi_bar) */
    double T111,T234,T414,   T212,T313,  T133,T423,T211,T233,T413;
    
    double r0=r[0],r2=r[2],r1=r[1];
    T111 = 0.5*h*tpsi*tpsi;
    /*  T234=  0.5*h*tpsi*tpsib;    */
    T234=  0.5*h*tpsi*tpsi;
    T414=T234;
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

/* This code was modified from the original BndMPoleSymplectic4Pass.c of AT
 to correctly integrate the Hamiltonian in  the curvilinear coordinate 
system of the dipole and to include the second order Transport map of the
 fringe field. New version created by Xiaobiao Huang in March 2009, in final
 verified version in August 2009.                                             */

void edge_fringe2A(double* r, double inv_rho, double edge_angle, double fint,
		   double gap,double h1,double K1);
void edge_fringe2B(double* r, double inv_rho, double edge_angle, double fint,
		   double gap,double h2,double K1);
void ATmultmv(double *r, const double* A);
void ATaddvv(double *r, const double *dr);

/*original kick function by Andrei Terebilo*/
static void bndthinkick0(double* r, double* A, double* B, double L,
			 double irho, int max_order)
{
  int    i;
  double ReSum = B[max_order];
  double ImSum = A[max_order];
  double ReSumTemp;
  /* recursively calculate the local transvrese magnetic field
     Bx = ReSum, By = ImSum
  */
  for(i = max_order-1; i >= 0; i--) {
    ReSumTemp = ReSum*r[0] - ImSum*r[2] + B[i];
    ImSum = ImSum*r[0] +  ReSum*r[2] + A[i];
    ReSum = ReSumTemp;
  }
  r[1] -=  L*(ReSum-(r[4]-r[0]*irho)*irho);
  r[3] +=  L*ImSum;
  r[5] +=  L*irho*r[0]; /* pathlength */
}

static void bndthinkick(double* r, double* A, double* B, double L, double h,
			int max_order)
/*****************************************************************************
(1) PolynomA is neglected.
(2) The vector potential is expanded up to 4th order of x and y. 
(3) Coefficients in PolynomB higher than 4th order is treated as if they are on straight geometry.
(4) The Hamiltonian is H2 = - h x delta - (1+h x)As/Brho-B0 x/Brho            */
{
  int    i;
  double ReSum = 0; /*B[max_order];*/
  double ImSum = 0; /*A[max_order];*/
    
  double ReSumTemp;
  double K1,K2;
    
  K1 = B[1];
  K2 = (max_order>=2) ? B[2] : 0;
    
  ReSum = B[max_order];
  for (i = max_order-1; i >= 0; i--) {
    ReSumTemp = ReSum*r[0] - ImSum*r[2] + B[i];
    ImSum = ImSum*r[0] +  ReSum*r[2] ;
    ReSum = ReSumTemp;
  }
    
  r[1] -=
    L*(-h*r[4] + ReSum + h*(h*r[0]+K1*(r[0]*r[0]-0.5*r[2]*r[2])
			    +K2*(r[0]*r[0]*r[0]-4.0/3.0*r[0]*r[2]*r[2])));
  r[3] +=
    L*(ImSum+h*(K1*r[0]*r[2]+4.0/3.0*K2*r[0]*r[0]*r[2]+(h/6.0*K1-K2/3.0)*r[2]
		*r[2]*r[2]));
  r[5] +=  L*h*r[0]; /* pathlength */
    
}

/* the pseudo-drift element described by Hamiltonian
   H1 = (1+hx) (px^2+py^2)/2(1+delta),                                        */
static void ATbendhxdrift6(double* r, double L,double h)
{
  double hs = h*L;
  double i1pd = 1.0/(1+r[4]);
  double x=r[0],px=r[1],py=r[3];

  r[0] += (1+h*x)*px*i1pd*L+1/4.*hs*L*(px*px-py*py)*i1pd*i1pd;
  /* (1.0/h+x)*((1.0+hs*px*i1pd/2.)*(1.0+hs*px*i1pd/2.)-(hs*py*i1pd/2.)
     *(hs*py*i1pd/2.))-1./h;*/
  r[1] -= hs*(px*px+py*py)*i1pd/2.0;
	
  r[2]+= (1.0+h*x)*i1pd*py*L*(1.+px*hs/2.0);
  r[5]+= (1.0+h*x)*i1pd*i1pd*L/2.0*(px*px+py*py);
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
      atdrift(ps_vec, Elem->Length);
      /* Check physical apertures at the exit of the magnet */
      if (Elem->RApertures) checkiflostRectangularAp(ps_vec, Elem->RApertures);
      if (Elem->EApertures) checkiflostEllipticalAp(ps_vec, Elem->EApertures);
      /* Misalignment at exit */
      if (Elem->R2) ATmultmv(ps_vec, Elem->R2);
      if (Elem->T2) ATaddvv(ps_vec, Elem->T2);
    }
  }
}

void MpolePass(double ps[], const int num_particles,
	       const struct elem_type *Elem)
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
      double norm = 1.0/(1.0+ps_vec[4]), NormL1 = L1*norm, NormL2 = L2*norm;

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
	fastdrift(ps_vec, NormL1);
	thin_kick(ps_vec, mpole->PolynomA, mpole->PolynomB, K1, mpole->irho,
		    mpole->MaxOrder);
	fastdrift(ps_vec, NormL2);
	thin_kick(ps_vec, mpole->PolynomA, mpole->PolynomB, K2, mpole->irho,
		    mpole->MaxOrder);
	fastdrift(ps_vec, NormL2);
	thin_kick(ps_vec, mpole->PolynomA, mpole->PolynomB, K1, mpole->irho,
		    mpole->MaxOrder);
	fastdrift(ps_vec, NormL1);
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

void MpoleRadPass(double ps[], const int num_particles,
		  const struct elem_type *Elem)
{
  int    k, m;
  double *ps_vect, SL, L1, L2, K1, K2;

  const elem_mpole *mpole = Elem->mpole_ptr;

  const bool
    useLinFrEleEntrance =
    (mpole->fringeIntM0 != NULL && mpole->fringeIntP0 != NULL
     && mpole->FringeQuadEntrance == 2),
    useLinFrEleExit =
    (mpole->fringeIntM0 != NULL && mpole->fringeIntP0 != NULL
     && mpole->FringeQuadExit == 2);

  SL = Elem->Length/mpole->NumIntSteps;
  L1 = SL*DRIFT1;
  L2 = SL*DRIFT2;
  K1 = SL*KICK1;
  K2 = SL*KICK2;
    
  if (mpole->KickAngle) {
    /* Convert corrector component to polynomial coefficients */
    mpole->PolynomB[0] -= sin(mpole->KickAngle[0])/Elem->Length;
    mpole->PolynomA[0] += sin(mpole->KickAngle[1])/Elem->Length;
  }

#pragma omp parallel for if (num_particles > OMP_PARTICLE_THRESHOLD)   \
  default(shared) shared(r, num_particles) private(c, ps_vect, m)

  for (k = 0; k < num_particles; k++) { /* Loop over particles  */
    ps_vect = ps+k*PS_DIM;
    if(!atIsNaN(ps_vect[0])) {
      /*  misalignment at entrance  */
      if (Elem->T1) ATaddvv(ps_vect, Elem->T1);
      if (Elem->R1) ATmultmv(ps_vect, Elem->R1);
      /* Check physical apertures at the entrance of the magnet */
      if (Elem->RApertures)
	checkiflostRectangularAp(ps_vect, Elem->RApertures);
      if (Elem->EApertures)
	checkiflostEllipticalAp(ps_vect, Elem->EApertures);
      if (mpole->FringeQuadEntrance && mpole->PolynomB[1] != 0e0) {
	if (useLinFrEleEntrance) /*Linear fringe fields from elegant*/
	  linearQuadFringeElegantEntrance
	    (ps_vect, mpole->PolynomB[1], mpole->fringeIntM0,
	     mpole->fringeIntP0);
	else
	  QuadFringePassP(ps_vect, mpole->PolynomB[1]);
      }
      /* integrator */
      for (m = 0; m < mpole->NumIntSteps; m++) { /* Loop over slices */
	ps_vect = ps+k*PS_DIM;
	ATdrift6(ps_vect, L1);
	strthinkickrad(ps_vect, mpole->PolynomA, mpole->PolynomB, K1,
		       mpole->Energy, mpole->MaxOrder);
	ATdrift6(ps_vect, L2);
	strthinkickrad(ps_vect, mpole->PolynomA, mpole->PolynomB, K2,
		       mpole->Energy, mpole->MaxOrder);
	ATdrift6(ps_vect, L2);
	strthinkickrad(ps_vect, mpole->PolynomA, mpole->PolynomB, K1,
		       mpole->Energy, mpole->MaxOrder);
	ATdrift6(ps_vect,L1);
      }
      if (mpole->FringeQuadExit && mpole->PolynomB[1]!=0) {
	if (useLinFrEleExit) /*Linear fringe fields from elegant*/
	  linearQuadFringeElegantExit(ps_vect, mpole->PolynomB[1],
				      mpole->fringeIntM0, mpole->fringeIntP0);
	else
	  QuadFringePassN(ps_vect, mpole->PolynomB[1]);
      }
      /* Check physical apertures at the exit of the magnet */
      if (Elem->RApertures) checkiflostRectangularAp(ps_vect, Elem->RApertures);
      if (Elem->EApertures) checkiflostEllipticalAp(ps_vect, Elem->EApertures);
      /* Misalignment at exit */
      if (Elem->R2) ATmultmv(ps_vect, Elem->R2);
      if (Elem->T2) ATaddvv(ps_vect, Elem->T2); 
    }
  }
  if (mpole->KickAngle) {
    /* Remove corrector component in polynomial coefficients */
    mpole->PolynomB[0] += sin(mpole->KickAngle[0])/Elem->Length; 
    mpole->PolynomA[0] -= sin(mpole->KickAngle[1])/Elem->Length;
  }
}

void BendRadPass(double ps[], const int num_particles,
		 const struct elem_type *Elem)
        
{	
  int    k, m;
  double *ps_vect, SL, L1, L2, K1, K2;

  const elem_mpole *mpole = Elem->mpole_ptr;

  const bool
    useLinFrEleEntrance =
    (mpole->fringeIntM0 != NULL && mpole->fringeIntP0 != NULL
     && mpole->FringeQuadEntrance == 2),
    useLinFrEleExit =
    (mpole->fringeIntM0 != NULL && mpole->fringeIntP0 != NULL
     && mpole->FringeQuadExit == 2);

  SL = Elem->Length/mpole->NumIntSteps;
  L1 = SL*DRIFT1;
  L2 = SL*DRIFT2;
  K1 = SL*KICK1;
  K2 = SL*KICK2;
    
  if (mpole->KickAngle) {
    /* Convert corrector component to polynomial coefficients */
    mpole->PolynomB[0] -= sin(mpole->KickAngle[0])/Elem->Length;
    mpole->PolynomA[0] += sin(mpole->KickAngle[1])/Elem->Length;
  }
#pragma omp parallel for if (num_particles > OMP_PARTICLE_THRESHOLD)	\
  default(shared) shared(r, num_particles) private(c, ps_vect, m)

  for(k = 0; k < num_particles; k++) { /* Loop over particles */
    ps_vect = ps+k*PS_DIM;
    if(!atIsNaN(ps_vect[0])) {
      /*  misalignment at entrance  */
      if (Elem->T1) ATaddvv(ps_vect, Elem->T1);
      if (Elem->R1) ATmultmv(ps_vect, Elem->R1);
      /* Check physical apertures at the entrance of the magnet */
      if (Elem->RApertures) checkiflostRectangularAp(ps_vect, Elem->RApertures);
      if (Elem->EApertures) checkiflostEllipticalAp(ps_vect, Elem->EApertures);
      /* edge focus */
      edge_fringe_entrance(ps_vect, mpole->irho, mpole->EntranceAngle,
			   mpole->FringeInt1, mpole->FullGap,
			   mpole->FringeBendEntrance);
      /* quadrupole gradient fringe */
      if (mpole->FringeQuadEntrance && mpole->PolynomB[1]!=0) {
	if (useLinFrEleEntrance) /*Linear fringe fields from elegant*/
	  linearQuadFringeElegantEntrance(ps_vect, mpole->PolynomB[1],
					  mpole->fringeIntM0,
					  mpole->fringeIntP0);
	else
	  QuadFringePassP(ps_vect, mpole->PolynomB[1]);
      }
      /* integrator  */
      for(m = 0; m < mpole->NumIntSteps; m++) { /* Loop over slices */
	drift6(ps_vect, L1);
	bndthinkickrad(ps_vect, mpole->PolynomA, mpole->PolynomB, K1,
		       mpole->irho, mpole->Energy, mpole->MaxOrder);
	drift6(ps_vect, L2);
	bndthinkickrad(ps_vect, mpole->PolynomA, mpole->PolynomB, K2,
		       mpole->irho, mpole->Energy, mpole->MaxOrder);
	drift6(ps_vect, L2);
	bndthinkickrad(ps_vect, mpole->PolynomA, mpole->PolynomB, K1,
		       mpole->irho, mpole->Energy, mpole->MaxOrder);
	drift6(ps_vect, L1);
      }
      /* quadrupole gradient fringe */
      if (mpole->FringeQuadExit && mpole->PolynomB[1]!=0) {
	if (useLinFrEleExit) /*Linear fringe fields from elegant*/
	  linearQuadFringeElegantExit(ps_vect, mpole->PolynomB[1],
				      mpole->fringeIntM0, mpole->fringeIntP0);
	else
	  QuadFringePassN(ps_vect, mpole->PolynomB[1]);
      }
      /* edge focus */
      edge_fringe_exit(ps_vect, mpole->irho, mpole->ExitAngle,
		       mpole->FringeInt2, mpole->FullGap,
		       mpole->FringeBendExit);
      /* Check physical apertures at the exit of the magnet */
      if (Elem->RApertures) checkiflostRectangularAp(ps_vect, Elem->RApertures);
      if (Elem->EApertures) checkiflostEllipticalAp(ps_vect, Elem->EApertures);
      /* Misalignment at exit */
      if (Elem->R2) ATmultmv(ps_vect, Elem->R2);
      if (Elem->T2) ATaddvv(ps_vect, Elem->T2);
    }
  }
  if (mpole->KickAngle) {
    /* Remove corrector component in polynomial coefficients */
    mpole->PolynomB[0] += sin(mpole->KickAngle[0])/Elem->Length;
    mpole->PolynomA[0] -= sin(mpole->KickAngle[1])/Elem->Length;
  }
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

void CorrectorPass(double ps[], const int num_particles,
		   const struct elem_type *Elem)
/* xkick, ykick - horizontal and vertical kicks in radiand 
   r - 6-by-N matrix of initial conditions reshaped into 
   1-d array of 6*N elements                                                  */
{
  int    k;
  double *ps_vec, xkick, ykick;

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

void WigPass(double ps[], const int num_particles, struct elem_type *Elem)
{
  int    k;
  double *ps_vec;

  for (k = 0; k < num_particles; k++) {
    ps_vec = ps+k*PS_DIM;
    if (!atIsNaN(ps_vec[0])) {
      switch (Elem->wig_ptr->Pmethod) {
      case second:
	GWigPass_2nd(Elem, ps_vec);
	break;
      case fourth:
	GWigPass_4th(Elem, ps_vec);
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
	GWigRadPass_2nd(Elem, ps_vec);
	break;
      case fourth:
	GWigRadPass_4th(Elem, ps_vec);
	break;
      default:
	printf("Invalid wiggler integration method %d.\n",
	       Elem->wig_ptr->Pmethod);
	break;
      }
    }
  }
}

void HamPass(double ps[], const int num_particles, const struct elem_type *Elem)
{
  int    k;
  double *ps_vect;

  for(k = 0; k < Elem->H_ptr->MaxOrder; k++) {
    Elem->H_ptr->F[2*k]   = Elem->H_ptr->PolynomB[k];
    Elem->H_ptr->F[2*k+1] = Elem->H_ptr->PolynomA[k];
  }

  for(k = 0; k < num_particles; k++) {
    ps_vect = ps+k*PS_DIM;
    if(!atIsNaN(ps_vect[0])) {
      /* misalignment at entrance */
      if (Elem->T1) ATaddvv(ps_vect, Elem->T1);
      if (Elem->R1) ATmultmv(ps_vect, Elem->R1);

      track_element(ps_vect, Elem);

      /* misalignment at exit */
      if (Elem->R2) ATmultmv(ps_vect, Elem->R2);
      if (Elem->T2) ATaddvv(ps_vect, Elem->T2);
    }
  }
}

void CBendPass(double ps[], const int num_particles, elem_type *Elem)
{
  int    k, m;
  double *ps_vect, SL, L1, L2, K1, K2;
  bool   useT1, useT2, useR1, useR2, useFringe1, useFringe2;

  SL = Elem->Length/Elem->mpole_ptr->NumIntSteps;
  L1 = SL*DRIFT1;
  L2 = SL*DRIFT2;
  K1 = SL*KICK1;
  K2 = SL*KICK2;

  /* mexPrintf("E0ref=%f\n",X0ref); */
  if(Elem->T1 == NULL)
    useT1 = false;
  else
    useT1 = true;
  if(Elem->T2 == NULL)
    useT2 = false;
  else
    useT2 = true;
  if(Elem->R1 == NULL)
    useR1 = false;
  else
    useR1 = true;
  if(Elem->R2 == NULL)
    useR2 = false;
  else
    useR2 = true;
  /* if either is 0 - do not calculate fringe effects */
  if(Elem->mpole_ptr->FringeInt1 == 0 || Elem->mpole_ptr->FullGap == 0)
    useFringe1 = false;
  else
    useFringe1 = true;
  if(Elem->mpole_ptr->FringeInt2 == 0 || Elem->mpole_ptr->FullGap == 0)
    useFringe2 = false;
  else
    useFringe2 = true;

  for(k = 0; k < num_particles; k++)	{   /* Loop over particles  */
    ps_vect = ps+k*PS_DIM;
    if(!atIsNaN(ps_vect[0])) {
      /*  misalignment at entrance  */
      if(useT1)	ATaddvv(ps_vect, Elem->T1);
      if(useR1)	ATmultmv(ps_vect, Elem->R1);
      /* edge focus */
      if(useFringe1)
	edgey_fringe(ps_vect,
		     Elem->mpole_ptr->irho+Elem->mpole_ptr->PolynomB[1]
		     *Elem->mpole_ptr->X0ref, Elem->mpole_ptr->EntranceAngle,
		     Elem->mpole_ptr->FringeInt1, Elem->mpole_ptr->FullGap);
      else
	edgey(ps_vect,
	      Elem->mpole_ptr->irho
	      +Elem->mpole_ptr->PolynomB[1]*Elem->mpole_ptr->X0ref,
	      Elem->mpole_ptr->EntranceAngle);
      /* Rotate and translate to straight Cartesian coordinate */
      E1rotation(ps_vect, Elem->mpole_ptr->X0ref,
		 Elem->mpole_ptr->EntranceAngle);
      /* integrator */
      for(m = 0; m < Elem->mpole_ptr->NumIntSteps; m++) {
	/* Loop over slices */
	ps_vect = ps+k*PS_DIM;
	ladrift6(ps_vect, L1);
	bndstrthinkick(ps_vect, Elem->mpole_ptr->PolynomA,
		       Elem->mpole_ptr->PolynomB, K1, Elem->mpole_ptr->irho,
		       Elem->mpole_ptr->MaxOrder);
	ladrift6(ps_vect, L2);
	bndstrthinkick(ps_vect, Elem->mpole_ptr->PolynomA,
		       Elem->mpole_ptr->PolynomB, K2, Elem->mpole_ptr->irho,
		       Elem->mpole_ptr->MaxOrder);
	ladrift6(ps_vect, L2);
	bndstrthinkick(ps_vect, Elem->mpole_ptr->PolynomA,
		       Elem->mpole_ptr->PolynomB, K1, Elem->mpole_ptr->irho,
		       Elem->mpole_ptr->MaxOrder);
	ladrift6(ps_vect, L1);
      }
      /* Rotate and translate back to curvilinear coordinate */
      E2rotation(ps_vect, Elem->mpole_ptr->X0ref, Elem->mpole_ptr->ExitAngle);
      ps_vect[5] -= Elem->mpole_ptr->RefDZ;
      if(useFringe2)
	edgey_fringe(ps_vect, Elem->mpole_ptr->irho
		     +Elem->mpole_ptr->PolynomB[1]*Elem->mpole_ptr->X0ref,
		     Elem->mpole_ptr->ExitAngle, Elem->mpole_ptr->FringeInt2,
		     Elem->mpole_ptr->FullGap);
      else    /* edge focus */
	edgey(ps_vect,
	      Elem->mpole_ptr->irho
	      +Elem->mpole_ptr->PolynomB[1]*Elem->mpole_ptr->X0ref,
	      Elem->mpole_ptr->ExitAngle);
      /* Misalignment at exit */
      if(useR2)	ATmultmv(ps_vect, Elem->R2);
      if(useT2)	ATaddvv(ps_vect, Elem->T2);
    }
  }
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
  if(mpole->FringeInt1 == 0 || mpole->FullGap == 0)
    useFringe1 = false;
  else
    useFringe1=true;
  if(mpole->FringeInt2 == 0 || mpole->FullGap == 0)
    useFringe2 = false;
  else
    useFringe2=true;
    
#pragma omp parallel for if (num_particles > OMP_PARTICLE_THRESHOLD)    \
  default(shared) shared(r, num_particles) private(c, ps_vec, m)

  for(k = 0; k < num_particles; k++) { /* Loop over particles  */
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
	bndthinkick(ps_vec, mpole->PolynomA, mpole->PolynomB, K1, mpole->irho,
		    mpole->MaxOrder);
	ATbendhxdrift6(ps_vec,L2, mpole->irho);
	bndthinkick(ps_vec, mpole->PolynomA, mpole->PolynomB, K2, mpole->irho,
		    mpole->MaxOrder);
	ATbendhxdrift6(ps_vec,L2, mpole->irho);
	bndthinkick(ps_vec, mpole->PolynomA, mpole->PolynomB,  K1, mpole->irho,
		    mpole->MaxOrder);
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
