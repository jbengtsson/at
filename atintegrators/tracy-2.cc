#include "tracy-2.h"

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

struct elem_type* init_mpole(const atElem *ElemData, struct elem_type *Elem,
			     const bool bend)
{
  int
    MaxOrder, NumIntSteps,  FringeBendEntrance, FringeBendExit,
    FringeQuadEntrance, FringeQuadExit;
  double
    Length, BendingAngle, EntranceAngle, ExitAngle, FullGap, FringeInt1,
    FringeInt2, *PolynomA, *PolynomB, *R1, *R2, *T1, *T2, *EApertures,
    *RApertures, *fringeIntM0, *fringeIntP0, *KickAngle;
  elem_mpole *mpole;

  Elem            = (struct elem_type*)malloc(sizeof(struct elem_type));
  Elem->mpole_ptr = (struct elem_mpole*)malloc(sizeof(struct elem_mpole));
  mpole           = Elem->mpole_ptr;

  Length = atGetDouble(ElemData, "Length");
  check_error();
  PolynomA = atGetDoubleArray(ElemData, (char*)"PolynomA");
  check_error();
  PolynomB = atGetDoubleArray(ElemData, (char*)"PolynomB");
  check_error();
  MaxOrder = atGetLong(ElemData, (char*)"MaxOrder");
  check_error();
  NumIntSteps = atGetLong(ElemData, (char*)"NumIntSteps");
  check_error();

  if (bend) {
    BendingAngle = atGetDouble(ElemData, (char*)"BendingAngle");
    check_error();
    EntranceAngle = atGetDouble(ElemData, (char*)"EntranceAngle");
    check_error();
    ExitAngle = atGetDouble(ElemData, (char*)(char*)"ExitAngle");
    check_error();
    FringeBendEntrance =
      atGetOptionalLong(ElemData, (char*)"FringeBendEntrance", 1);
    check_error();
    FringeBendExit = atGetOptionalLong(ElemData, (char*)"FringeBendExit", 1);
    check_error();
    FullGap = atGetOptionalDouble(ElemData, (char*)"FullGap", 0);
    check_error();
    FringeInt1 = atGetOptionalDouble(ElemData, (char*)"FringeInt1", 0);
    check_error();
    FringeInt2 = atGetOptionalDouble(ElemData, (char*)"FringeInt2", 0);
    check_error();
  }

  FringeQuadEntrance =
    atGetOptionalLong(ElemData, (char*)"FringeQuadEntrance", 0);
  check_error();
  FringeQuadExit = atGetOptionalLong(ElemData, (char*)"FringeQuadExit", 0);
  check_error();
  fringeIntM0 = atGetOptionalDoubleArray(ElemData, (char*)"fringeIntM0");
  check_error();
  fringeIntP0 = atGetOptionalDoubleArray(ElemData, (char*)"fringeIntP0");
  check_error();
  
  R1 = atGetOptionalDoubleArray(ElemData, (char*)"R1");
  check_error();
  R2 = atGetOptionalDoubleArray(ElemData, (char*)"R2");
  check_error();
  T1 = atGetOptionalDoubleArray(ElemData, (char*)"T1");
  check_error();
  T2 = atGetOptionalDoubleArray(ElemData, (char*)"T2");
  check_error();

  EApertures = atGetOptionalDoubleArray(ElemData, (char*)"EApertures");
  check_error();
  RApertures = atGetOptionalDoubleArray(ElemData, (char*)"RApertures");
  check_error();
  KickAngle = atGetOptionalDoubleArray(ElemData, (char*)"KickAngle");
  check_error();
        
  Elem->Length       = Length;
  Elem->R1           = R1;
  Elem->R2           = R2;
  Elem->T1           = T1;
  Elem->T2           = T2;
  Elem->EApertures   = EApertures;
  Elem->RApertures   = RApertures;

  mpole->PolynomA    = PolynomA;
  mpole->PolynomB    = PolynomB;
  mpole->MaxOrder    = MaxOrder;
  mpole->NumIntSteps = NumIntSteps;

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

  mpole->irho = (Elem->Length != 0e0)? mpole->BendingAngle/Elem->Length : 0e0;

  mpole->FringeQuadEntrance = FringeQuadEntrance;
  mpole->FringeQuadExit     = FringeQuadExit;
  mpole->fringeIntM0        = fringeIntM0;
  mpole->fringeIntP0        = fringeIntP0;
  mpole->KickAngle          = KickAngle;

  return Elem;
}


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

void Cav_Pass(const double L, const double f_RF, const double VoE0,
	      const double phi, std::vector<double> &ps)
{
  double delta;

  Drift(L/2e0, ps);
  if (VoE0 != 0e0) {
    delta = -VoE0*sin(2e0*M_PI*f_RF/C0*ps[ct_]+phi);
    ps[delta_] += delta;

    // if (globval.radiation) globval.dE -= is_double<T>::cst(delta);

    // if (globval.pathlength) ps[ct_] -= C->Ph/C->Pfreq*c0;
  }
  Drift(L/2e0, ps);
}

#if 1

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

void strthinkick(double ps[], const double a[], const double b[],
		 const double L, const int n_max)
{
  double              bn[2*HOMmax+1];
  std::vector<double> ps_stl = arrtostl(ps);

  for (int k = n_max+1; k > 0; k--) {
    bn[HOMmax+k] = b[k-1];
    bn[HOMmax-k] = a[k-1];
  }
  thin_kick(n_max+1, bn, L, 0e0, 0e0, ps_stl);
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

void bndthinkick(double ps[], const double a[], const double b[],
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

void cav_pass(double ps[], const double L, const double VoE0,
	      const double f_RF, const double lag)
{
  std::vector<double> ps_stl = arrtostl(ps);

  const double phi = -lag*2e0*M_PI*f_RF/C0;

  Cav_Pass(L, f_RF, VoE0, phi, ps_stl);
  stltoarr(ps_stl, ps);
}

#else

static void atdrift(double* r, double L)
/* Input parameter L is the physical length
   1/(1+delta) normalization is done internally                               */
{
  double p_norm = 1/(1+r[4]); 
  double NormL  = L*p_norm;   
  r[0]+= NormL*r[1]; 
  r[2]+= NormL*r[3];
  r[5]+= NormL*p_norm*(r[1]*r[1]+r[3]*r[3])/2;
}

static void fastdrift(double* r, double NormL)

/* NormL=(Physical Length)/(1+delta)  is computed externally to speed up
   calculations  in the loop if momentum deviation (delta) does not change
   such as in 4-th order symplectic integrator w/o radiation                  */

{
  r[0] += NormL*r[1];
  r[2] += NormL*r[3];
  r[5] += NormL*(r[1]*r[1]+r[3]*r[3])/(2*(1+r[4]));
}

static void strthinkick(double* r, const double* A, const double* B, double L,
			int max_order)
/***************************************************************************** 
 Calculate and apply a multipole kick to a 6-dimentional
 phase space vector in a straight element (quadrupole)
 
 IMPORTANT !!!
 The reference coordinate system is straight but the field expansion may still
 contain dipole terms: PolynomA(1), PolynomB(1) - in MATLAB notation,
 A[0], B[0] - C,C++ notation
 
******************************************************************************/
{
  int i;
  double ReSum = B[max_order];
  double ImSum = A[max_order];
  double ReSumTemp;
  for (i = max_order-1; i >= 0; i--) {
    ReSumTemp = ReSum*r[0] - ImSum*r[2] + B[i];
    ImSum = ImSum*r[0] +  ReSum*r[2] + A[i];
    ReSum = ReSumTemp;
  }
  r[1] -= L*ReSum;
  r[3] += L*ImSum;
}

static void edge_fringe(double r[], const double inv_rho,
			const double edge_angle, const double fint,
			const double gap, const int method, const bool hor)
{
  /*     method 0 no fringe field
   *     method 1 legacy version Brown First Order
   *     method 2 SOLEIL close to second order of Brown
   *     method 3 THOMX
   */
  double fringecorr, fx, fy;
  /* Fringe field correction */
  if ((fint==0.0) || (gap==0.0) || (method==0))
    fringecorr = 0.0;
  else {
    double sedge = sin(edge_angle);
    double cedge = cos(edge_angle);
    fringecorr = inv_rho*gap*fint*(1+sedge*sedge)/cedge;
  }
    
  /* Edge angle focusing */
  fx = inv_rho*tan(edge_angle);
  if (method==1)
    fy = inv_rho*tan(edge_angle-fringecorr/(1+r[4]));
  else if (method==2)
    fy = inv_rho*tan(edge_angle-fringecorr/(1+r[4]))/(1+r[4]);
  else if (method==3)
    if (hor)
      fy = inv_rho*tan(edge_angle-fringecorr+r[1]/(1+r[4]));
    else
      fy = inv_rho*tan(edge_angle-fringecorr-r[1]/(1+r[4]));
  else    /* fall back to legacy version */
    fy = inv_rho*tan(edge_angle-fringecorr/(1+r[4]));

  r[1] += r[0]*fx;
  r[3] -= r[2]*fy;
}

static void bndthinkick(double* r, double* A, double* B, double L, double irho,
			int max_order)
/***************************************************************************** 
Calculate multipole kick in a curved elemrnt (bending magnet)
The reference coordinate system  has the curvature given by the inverse 
(design) radius irho.
IMPORTANT !!!
The magnetic field Bo that provides this curvature MUST NOT be included in the dipole term
PolynomB[1](MATLAB notation)(C: B[0] in this function) of the By field expansion

The kick is given by

           e L      L delta      L x
theta  = - --- B  + -------  -  -----  , 
     x     p    y     rho           2
            0                    rho

         e L
theta  = --- B
     y    p   x
           0

*************************************************************************/
{
  int    i;
  double ReSum = B[max_order];
  double ImSum = A[max_order];
  double ReSumTemp;
  /* recursively calculate the local transverse magnetic field
   * Bx = ReSum, By = ImSum
   */
  for (i = max_order-1; i >= 0; i--) {
    ReSumTemp = ReSum*r[0] - ImSum*r[2] + B[i];
    ImSum = ImSum*r[0] +  ReSum*r[2] + A[i];
    ReSum = ReSumTemp;
  }
  r[1] -= L*(ReSum-(r[4]-r[0]*irho)*irho);
  r[3] += L*ImSum;
  r[5] += L*irho*r[0]; /* pathlength */
}

void cav_pass(double r_in[], const double le, const double nv,
	      const double freq, const double lag)
{
  double p_norm, NormL;

  if (le == 0)
    r_in[4] += -nv*sin(TWOPI*freq*(r_in[5]-lag)/C0);
  else {
    p_norm = 1/(1+r_in[4]);
    NormL  = le/2*p_norm;
    /* Prropagate through a drift equal to half cavity length */
    r_in[0] += NormL*r_in[1];
    r_in[2] += NormL*r_in[3];
    r_in[5] += NormL*p_norm*(r_in[1]*r_in[1]+r_in[3]*r_in[3])/2;
    /* Longitudinal momentum kick */
    r_in[4] += -nv*sin(TWOPI*freq*(r_in[5]-lag)/C0);
    p_norm = 1/(1+r_in[4]);
    NormL  = le/2*p_norm;
    /* Prropagate through a drift equal to half cavity length */
    r_in[0] += NormL*r_in[1];
    r_in[2] += NormL*r_in[3];
    r_in[5] += NormL*p_norm*(r_in[1]*r_in[1]+r_in[3]*r_in[3])/2;
  }
}

#endif

static void edge_fringe_entrance(double* r, double inv_rho, double edge_angle,
				 double fint, double gap, int method)
{ edge_fringe(r, inv_rho, edge_angle, fint, gap, method, true); }

static void edge_fringe_exit(double* r, double inv_rho, double edge_angle,
			     double fint, double gap, int method)
{ edge_fringe(r, inv_rho, edge_angle, fint, gap, method, false); }


static void QuadFringePassP(double* r, const double b2)
{
  /* x=r[0],px=r[1],y=r[2],py=r[3],delta=r[4],ct=r[5]
     Lee-Whiting's thin lens limit formula as given in p. 390 of "Beam
     Dynamics..."by E. Forest                                                 */
  double u     = b2/(12.0*(1.0+r[4]));
  double x2    = r[0]*r[0];
  double z2    = r[2]*r[2];
  double xz    = r[0]*r[2];
  double gx    = u * (x2+3*z2) * r[0];
  double gz    = u * (z2+3*x2) * r[2];
  double r1tmp = 0;
  double r3tmp = 0;

  r[0] += gx;
  r1tmp = 3*u*(2*xz*r[3]-(x2+z2)*r[1]);
   
  r[2] -= gz;
   
  r3tmp = 3*u*(2*xz*r[1]-(x2+z2)*r[3]);
  r[5] -= (gz*r[3] - gx*r[1])/(1+r[4]);
   
  r[1] += r1tmp;
  r[3] -= r3tmp;
}

static void QuadFringePassN(double* r, const double b2)
{
  /* x=r[0],px=r[1],y=r[2],py=r[3],delta=r[4],ct=r[5] 
     Lee-Whiting's thin lens limit formula as given in p. 390 of "Beam
     Dynamics..."by E. Forest                                                 */
  double u     = b2/(12.0*(1.0+r[4]));
  double x2    = r[0]*r[0];
  double z2    = r[2]*r[2];
  double xz    = r[0]*r[2];
  double gx    = u * (x2+3*z2) * r[0];
  double gz    = u * (z2+3*x2) * r[2];
  double r1tmp = 0;
  double r3tmp = 0;

  r[0] -= gx;
  r1tmp = 3*u*(2*xz*r[3]-(x2+z2)*r[1]);
   
  r[2] += gz;
   
  r3tmp = 3*u*(2*xz*r[1]-(x2+z2)*r[3]);
  r[5] += (gz*r[3] - gx*r[1])/(1+r[4]);
   
  r[1] -= r1tmp;
  r[3] += r3tmp;
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
	       const struct elem_type *Elem)
{
  int    k;
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

  for(k = 0; k < num_particles; k++) {	/* Loop over particles  */
    ps_vec = ps+k*6;
    if (!atIsNaN(ps_vec[0])) {
      int    m;
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
      for(m = 0; m < mpole->NumIntSteps; m++) { /* Loop over slices*/
	fastdrift(ps_vec, NormL1);
	bndthinkick(ps_vec, mpole->PolynomA, mpole->PolynomB, K1, mpole->irho,
		    mpole->MaxOrder);
	fastdrift(ps_vec, NormL2);
	bndthinkick(ps_vec, mpole->PolynomA, mpole->PolynomB, K2, mpole->irho,
		    mpole->MaxOrder);
	fastdrift(ps_vec, NormL2);
	bndthinkick(ps_vec, mpole->PolynomA, mpole->PolynomB, K1, mpole->irho,
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
