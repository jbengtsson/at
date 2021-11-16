/*   File: atphyslib.c
 *   Common physics functions for Accelerator Toolbox
 *   A.Terebilo   10/28/04
 *
 *   functions edge_fringe2A and edge_fringe2B were added by Xiaobiao Huang, August 2009
 *
 *   Two additional methods for bending magnet fringe fields added, February 2017
 *   method 1 legacy version Brown First Order
 *   Version 2 SOLEIL close to second order of Brown
 *   Version 3 THOMX
 */

#include <math.h>


#include <armadillo>


#define PS_DIM 6


#define sqr(x)  ((x)*(x))
#define cube(x) ((x)*(x)*(x))

#define HOMmax  21

enum spatial_ind
  { X_ = 0,
    Y_ = 1,
    Z_ = 2 };

enum phase_space_ind
  { x_     = 0,
    px_    = 1,
    y_     = 2,
    py_    = 3,
    delta_ = 4,
    ct_    = 5 };

enum MpoleKind
  { All    = 0,
    Dip    = 1,
    Quad   = 2,
    Sext   = 3,
    Oct    = 4,
    Dec    = 5,
    Dodec  = 6 };

double C_u, C_gamma, C_q, cl_rad, q_fluct;


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

#if 1

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

inline void atdrift(double ps[], const double L)
{
  std::vector<double> ps_stl(PS_DIM, 0e0);

  ps_stl = arrtostl(ps);
  Drift(L, ps_stl);
  stltoarr(ps_stl, ps);
}

inline void fastdrift(double ps[], const double L)
{
  std::vector<double> ps_stl(PS_DIM, 0e0);

  ps_stl = arrtostl(ps);
  Drift(L, ps_stl);
  stltoarr(ps_stl, ps);
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
    //   irho*tan(phi*M_PI/180e0-get_psi(irho, phi, gap))*ps[y_]/(1e0+ps[delta_]);
    // Leading order correction.
    ps[py_] -=
      irho*tan(phi*M_PI/180e0-get_psi(irho, phi, gap))*ps[y_]*(1e0-ps[delta_]);
  } else
    ps[py_] -= irho*tan(phi*M_PI/180e0-get_psi(irho, phi, gap))*ps[y_];
}


template<typename T>
void p_rot(double phi, std::vector<double> &ps)
{
  double              c, s, t, pz, p, val;
  std::vector<double> ps1;

  c = cos(phi*M_PI/180e0); s = sin(phi*M_PI/180e0); t = tan(phi*M_PI/180e0);
  pz = get_p_s(ps);

  if (true) {
    ps[px_] = s*pz + c*ps[px_];
  } else {
    // ps1 = ps; p = c*pz - s*ps1[px_];
    // px[x_]   = ps1[x_]*pz/p; px[px_] = s*pz + c*ps1[px_];
    // px[y_]  += ps1[x_]*ps1[py_]*s/p;
    // px[ct_] += (1e0+ps1[delta_])*ps1[x_]*s/p;

    ps1 = ps; val = 1e0 - ps1[px_]*t/pz;
    ps[x_]  = ps1[x_]/(c*val);
    ps[px_] = ps1[px_]*c + s*pz;
    ps[y_]  = ps1[y_] + t*ps1[x_]*ps1[py_]/(pz*val);
    ps[ct_] = ps1[ct_] + ps1[x_]*(1e0+ps1[delta_])*t/(pz*val);
  }
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

void strthinkick(double ps[], const double a[], const double b[],
		 const double L, const int n_max)
{
  double              bn[2*HOMmax+1];
  std::vector<double> ps_stl(PS_DIM, 0e0);

  for (int k = n_max+1; k > 0; k--) {
    bn[HOMmax+k] = b[k-1];
    bn[HOMmax-k] = a[k-1];
  }
  ps_stl = arrtostl(ps);
  thin_kick(n_max+1, bn, L, 0e0, 0e0, ps_stl);
  stltoarr(ps_stl, ps);
}

void bndthinkick(double ps[], const double a[], const double b[],
		 const double L, const double irho, const int n_max)
{
  double              bn[2*HOMmax+1];
  std::vector<double> ps_stl(PS_DIM, 0e0);

  for (int k = n_max+1; k > 0; k--) {
    bn[HOMmax+k] = b[k-1];
    bn[HOMmax-k] = a[k-1];
  }
  ps_stl = arrtostl(ps);
  thin_kick(n_max+1, bn, L, irho, irho, ps_stl);
  stltoarr(ps_stl, ps);
}

#else

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
  for (i=max_order-1; i>=0; i--) {
    ReSumTemp = ReSum*r[0] - ImSum*r[2] + B[i];
    ImSum = ImSum*r[0] +  ReSum*r[2] + A[i];
    ReSum = ReSumTemp;
  }
  r[1] -=  L*ReSum;
  r[3] +=  L*ImSum;
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
  int i;
  double ReSum = B[max_order];
  double ImSum = A[max_order];
  double ReSumTemp;
  /* recursively calculate the local transverse magnetic field
   * Bx = ReSum, By = ImSum
   */
  for (i=max_order-1; i>=0; i--) {
    ReSumTemp = ReSum*r[0] - ImSum*r[2] + B[i];
    ImSum = ImSum*r[0] +  ReSum*r[2] + A[i];
    ReSum = ReSumTemp;
  }
  r[1] -=  L*(ReSum-(r[4]-r[0]*irho)*irho);
  r[3] +=  L*ImSum;
  r[5] +=  L*irho*r[0]; /* pathlength */
}

#endif

static void edge_fringe_entrance(double* r, double inv_rho, double edge_angle,
        double fint, double gap, int method)
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
        fy = inv_rho*tan(edge_angle-fringecorr+r[1]/(1+r[4]));
    else    /* fall back to legacy version */
        fy = inv_rho*tan(edge_angle-fringecorr/(1+r[4]));

    r[1]+=r[0]*fx;
    r[3]-=r[2]*fy;
}

static void edge_fringe_exit(double* r, double inv_rho, double edge_angle,
        double fint, double gap, int method)
{
    /*     method 0 no fringe field
     *     method 1 legacy version Brown First Order
     *     method 2 SOLEIL close to second order of Brown
     *     method 3 THOMX
     */
    /* Fringe field correction */
    double fringecorr, fx, fy;
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
        fy = inv_rho*tan(edge_angle-fringecorr-r[1]/(1+r[4]));
    else    /* fall back to legacy version */
        fy = inv_rho*tan(edge_angle-fringecorr/(1+r[4]));

    r[1]+=r[0]*fx;
    r[3]-=r[2]*fy;
}
