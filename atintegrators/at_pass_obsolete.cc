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

static void thin_kick(double* r, double* A, double* B, double L, double irho,
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
