
#include "atelem.cc"
#include "atlalib.cc"


struct elem {
  double
  Length,
    *R1,
    *R2,
    *T1,
    *T2,
    *EApertures,
    *RApertures;
};

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

void DriftPass(double ps_n[], double le, const double t1[], const double t2[],
	       const double r1[], const double r2[], double RApertures[],
	       double EApertures[], int num_particles)
/* le - physical length
   ps_in - 6-by-N matrix of initial conditions reshaped into 
   1-d array of 6*N elements 
*/
{
  int                 j, k;
  arma::vec           ps_vec(PS_DIM);
  std::vector<double> ps(PS_DIM, 0e0);

#pragma omp parallel for if (num_particles > OMP_PARTICLE_THRESHOLD*10) \
  default(shared) shared(ps_in,num_particles) private(c,r6)
  for (j = 0; j < num_particles; j++) { /*Loop over particles  */
    // r6 = ps_in+c*6;
    for (k = 0; k < PS_DIM; k++)
      ps_vec(k) = ps_n[j*PS_DIM+k];
    if (!isnan(ps_vec(0))) {
      /*  misalignment at entrance  */
      if (t1) ps_vec = ps_vec + arrtovec(t1);
      if (r1) ps_vec = arrtomat(r1)*ps_vec;
#if 0
      /* Check physical apertures at the entrance of the magnet */
      if (RApertures) checkiflostRectangularAp(r6, RApertures);
      if (EApertures) checkiflostEllipticalAp(r6, EApertures);
#endif
      // ATdrift6(r6, le);
      ps = vectostl(ps_vec);
      Drift(le, ps);
      ps_vec = stltovec(ps);
#if 0
      /* Check physical apertures at the exit of the magnet */
      if (RApertures) checkiflostRectangularAp(r6, RApertures);
      if (EApertures) checkiflostEllipticalAp(r6, EApertures);
#endif
      /* Misalignment at exit */
      if (r2) ps_vec = arrtomat(r2)*ps_vec;
      if (t2) ps_vec = ps_vec + arrtovec(t2);
      for (k = 0; k < PS_DIM; k++)
	ps_n[j*PS_DIM+k] = ps_vec(k);
    }
  }
}

struct elem*
trackFunction(const atElem *ElemData,struct elem *Elem, double ps_in[],
	      int num_particles, struct parameters *Param)
{
  /*  if (ElemData) {*/
  if (!Elem) {
    double Length;
    double *R1, *R2, *T1, *T2, *EApertures, *RApertures;
    Length=atGetDouble(ElemData,"Length"); check_error();
    R1=atGetOptionalDoubleArray(ElemData, (char*)"R1"); check_error();
    R2=atGetOptionalDoubleArray(ElemData, (char*)"R2"); check_error();
    T1=atGetOptionalDoubleArray(ElemData, (char*)"T1"); check_error();
    T2=atGetOptionalDoubleArray(ElemData, (char*)"T2"); check_error();
    EApertures=atGetOptionalDoubleArray(ElemData, (char*)"EApertures");
    check_error();
    RApertures=atGetOptionalDoubleArray(ElemData, (char*)"RApertures");
    check_error();
    Elem = (struct elem*)atMalloc(sizeof(struct elem));
    Elem->Length=Length;
    Elem->R1=R1;
    Elem->R2=R2;
    Elem->T1=T1;
    Elem->T2=T2;
    Elem->EApertures=EApertures;
    Elem->RApertures=RApertures;
  }
  DriftPass(ps_in, Elem->Length, Elem->T1, Elem->T2, Elem->R1, Elem->R2,
	    Elem->RApertures, Elem->EApertures, num_particles);
  /*  }
      else {
      atFree(Elem->T1);
      atFree(Elem->T2);
      atFree(Elem->R1);
      atFree(Elem->R2);
      atFree(Elem->EApertures);
      atFree(Elem->RApertures);
      }*/
  return Elem;
}
