/* CavityPass.c
 * Accelerator Toolbox
 * Revision 3/10/04
 * A.Terebilo terebilo@ssrl.slac.stanford.edu
 */

#include "atelem.cc"
#include "atphyslib.cc"


void CavityPass(double r_in[], const double le, const double nv,
		const double freq, const double lag, const int num_particles)
/* le - physical length
 * nv - peak voltage (V) normalized to the design enegy (eV)
 * r is a 6-by-N matrix of initial conditions reshaped into
 * 1-d array of 6*N elements
 */
{
  int c, c6;

  for(c = 0; c < num_particles; c++) {
    c6 = c*6;
    if(!atIsNaN(r_in[c6]))
      cav_pass(&r_in[c6], le, nv, freq, lag);
  }
}

struct elem_type* init_cav(const atElem *ElemData, struct elem_type *Elem)
{
  double   Length, Voltage, Energy, Frequency, TimeLag;
  elem_cav *cav;

  Elem          = (struct elem_type*)malloc(sizeof(struct elem_type));
  Elem->cav_ptr = (struct elem_cav*)malloc(sizeof(struct elem_cav));
  cav           = Elem->cav_ptr;

  Length    = atGetDouble(ElemData,"Length"); check_error();
  Voltage   = atGetDouble(ElemData,"Voltage"); check_error();
  Energy    = atGetDouble(ElemData,"Energy"); check_error();
  Frequency = atGetDouble(ElemData,"Frequency"); check_error();
  TimeLag   = atGetOptionalDouble(ElemData,"TimeLag",0); check_error();

  Elem->Length    = Length;
  cav->Voltage    = Voltage;
  cav->Energy     = Energy;
  cav->Frequency  = Frequency;
  cav->TimeLag    = TimeLag;

  return Elem;
}

struct elem_type*
trackFunction(const atElem *ElemData, struct elem_type *Elem,
	      double *r_in, int num_particles, struct parameters *Param)
{
  elem_cav *cav;

  if (!Elem) Elem = init_cav(ElemData, Elem);

  cav = Elem->cav_ptr;

  CavityPass(r_in, Elem->Length, cav->Voltage/cav->Energy, cav->Frequency,
	     cav->TimeLag, num_particles);

  return Elem;
}
