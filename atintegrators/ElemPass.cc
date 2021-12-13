#include "tracy-2.cc"

extern "C" struct elem_type*
corr_pass(const PyObject *ElemData, struct elem_type *Elem, double ps[],
	  const int num_particles, const struct parameters *Param)
{
  if (!Elem) Elem = init_corr(ElemData, Elem);
  if (Elem) {
    CorrectorPass(ps, num_particles, Elem);
    return Elem;
  } else
    return NULL;
}

extern "C" struct elem_type*
id_pass(const PyObject *ElemData, struct elem_type *Elem, double ps[],
	const int num_particles, const struct parameters *Param)
{
  if (!Elem) Elem = init_id(ElemData, Elem);
  if (Elem) {
    IdentityPass(ps, num_particles, Elem);
    return Elem;
  } else
    return NULL;
}

extern "C" struct elem_type*
aper_pass(const PyObject *ElemData, struct elem_type *Elem, double ps[],
	  const int num_particles, const struct parameters *Param)
{
  if (!Elem) Elem = init_ap(ElemData, Elem);
  if (Elem) {
    AperturePass(ps, num_particles, Elem);
    return Elem;
  } else
    return NULL;
}

extern "C" struct elem_type*
drift_pass(const PyObject *ElemData, struct elem_type *Elem, double ps[],
	   const int num_particles, const struct parameters *Param)
{
  if (!Elem) Elem = init_drift(ElemData, Elem);
  if (Elem) {
    DriftPass(ps, num_particles, Elem);
    return Elem;
  } else
    return NULL;
}

extern "C" struct elem_type*
mpole_pass(const PyObject *ElemData, struct elem_type *Elem, double ps[],
	   const int num_particles, const struct parameters *Param)
{
  if (!Elem) Elem = init_mpole(ElemData, Elem, false, false, false, false);
  if (Elem) {
    MpolePass(ps, num_particles, Elem, false);
    return Elem;
  } else
    return NULL;
}

extern "C" struct elem_type*
mpole_rad_pass(const PyObject *ElemData, struct elem_type *Elem, double ps[],
	       const int num_particles, const struct parameters *Param)
{
  if (!Elem) Elem = init_mpole(ElemData, Elem, false, false, true, false);
  if (Elem) {
    MpolePass(ps, num_particles, Elem, true);
    return Elem;
  } else
    return NULL;
}

extern "C" struct elem_type*
bend_pass(const PyObject *ElemData, struct elem_type *Elem, double *ps,
	  const int num_particles, const struct parameters *Param)
{
  if (!Elem) Elem = init_mpole(ElemData, Elem, true, false, false, false);
  if (Elem) {
    MpolePass(ps, num_particles, Elem, false);
    return Elem;
  } else
    return NULL;
}

extern "C" struct elem_type*
bend_rad_pass(const PyObject *ElemData, struct elem_type *Elem, double ps[],
	      const int num_particles, const struct parameters *Param)
{
  if (!Elem) Elem = init_mpole(ElemData, Elem, true, false, true, false);
  if (Elem) {
    MpolePass(ps, num_particles, Elem, true);
    return Elem;
  } else
    return NULL;
}

extern "C" struct elem_type*
bend_exact_pass(const PyObject *ElemData, struct elem_type *Elem,double ps[],
		const int num_particles, const struct parameters *Param)
{
  if (!Elem) Elem = init_mpole(ElemData, Elem, true, false, false, true);
  if (Elem) {
    MpoleE2Pass(ps, num_particles, Elem);
    return Elem;
  } else
    return NULL;
}

extern "C" struct elem_type*
cbend_pass(const PyObject *ElemData, struct elem_type *Elem, double ps[],
	   const int num_particles, const struct parameters *Param)
{
  if (!Elem) Elem = init_mpole(ElemData, Elem, true, true, false, false);
  if (Elem) {
    CBendPass(ps, num_particles, Elem);
    return Elem;
  } else
    return NULL;
}

extern "C" struct elem_type*
cav_pass(const PyObject *ElemData, struct elem_type *Elem, double *ps,
	 const int num_particles, const struct parameters *Param)
{
  if (!Elem) Elem = init_cav(ElemData, Elem);
  if (Elem) {
    CavityPass(ps, num_particles, Elem);
    return Elem;
  } else
    return NULL;
}

extern "C" struct elem_type*
wig_pass(const PyObject *ElemData, struct elem_type *Elem, double ps[],
	 const int num_particles, const struct parameters *Param)
{
  if (!Elem) Elem = init_wig(ElemData, Elem, false);
  if (Elem) {
    WigPass(ps, num_particles, Elem);
    return Elem;
  } else
    return NULL;
}

extern "C" struct elem_type*
wig_rad_pass(const PyObject *ElemData, struct elem_type *Elem, double ps[],
	     const int num_particles, const struct parameters *Param)

{
  if (!Elem) Elem = init_wig(ElemData, Elem, true);
  if (Elem) {
    WigRadPass(ps, num_particles, Elem);
    return Elem;
  } else
    return NULL;
}

extern "C" struct elem_type*
M66_pass(const PyObject *ElemData,	struct elem_type *Elem, double ps[],
	 const int num_particles, const struct parameters *Param)
{
  if (!Elem) Elem = init_M66(ElemData, Elem);
  if (Elem) {
    Matrix66Pass(ps, num_particles, Elem);
    return Elem;
  } else
    return NULL;
}

extern "C" struct elem_type*
H_exact_pass(const PyObject *ElemData, struct elem_type *Elem, double ps[],
	     const int num_particles, const struct parameters *Param)
{
  if (!Elem) Elem = init_H(ElemData, Elem);
  if (Elem) {
    HamPass(ps, num_particles, Elem);
    return Elem;
  } else
    return NULL;
}
