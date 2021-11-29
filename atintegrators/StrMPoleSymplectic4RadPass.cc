#include "elem.cc"
#include "tracy-2.cc"

struct elem
{
  double Length;
  double *PolynomA;
  double *PolynomB;
  int MaxOrder;
  int NumIntSteps;
  double Energy;
  /* Optional fields */
  int FringeQuadEntrance;
  int FringeQuadExit;
  double *fringeIntM0;
  double *fringeIntP0;
  double *R1;
  double *R2;
  double *T1;
  double *T2;
  double *RApertures;
  double *EApertures;
  double *KickAngle;
};


struct elem_type*
trackFunction(const PyObject *ElemData, struct elem_type *Elem, double ps[],
	      const int num_particles, const struct parameters *Param)
{
  if (!Elem) Elem = init_mpole(ElemData, Elem, false, false, true);
  if (Elem) {
    MpoleRadPass(ps, num_particles, Elem);
    return Elem;
  } else
    return NULL;
}
