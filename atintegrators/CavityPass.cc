/* CavityPass.c
 * Accelerator Toolbox
 * Revision 3/10/04
 * A.Terebilo terebilo@ssrl.slac.stanford.edu
 */

#include "atelem.cc"
#include "atphyslib.cc"


struct elem
{
  double Length;
  double Voltage;
  double Energy;
  double Frequency;
  /* optional fields */
  double TimeLag;
};


#if defined(MATLAB_MEX_FILE) || defined(PYAT)
ExportMode struct elem*
trackFunction(const atElem *ElemData,struct elem *Elem,
	      double *r_in, int num_particles, struct parameters *Param)
{
  if (!Elem) {
    double Length, Voltage, Energy, Frequency, TimeLag;
    Length=atGetDouble(ElemData,"Length"); check_error();
    Voltage=atGetDouble(ElemData,"Voltage"); check_error();
    Energy=atGetDouble(ElemData,"Energy"); check_error();
    Frequency=atGetDouble(ElemData,"Frequency"); check_error();
    TimeLag=atGetOptionalDouble(ElemData,"TimeLag",0); check_error();
    Elem = (struct elem*)atMalloc(sizeof(struct elem));
    Elem->Length=Length;
    Elem->Voltage=Voltage;
    Elem->Energy=Energy;
    Elem->Frequency=Frequency;
    Elem->TimeLag=TimeLag;
  }
  CavityPass(r_in,Elem->Length,Elem->Voltage/Elem->Energy,Elem->Frequency,
	     Elem->TimeLag,num_particles);
  return Elem;
}

MODULE_DEF(CavityPass)        /* Dummy module initialisation */

#endif /*defined(MATLAB_MEX_FILE) || defined(PYAT)*/

#if defined(MATLAB_MEX_FILE)

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  if (nrhs == 2) {
    double *r_in;
    const mxArray *ElemData = prhs[0];
    int num_particles = mxGetN(prhs[1]);
    double Length=atGetDouble(ElemData,"Length");
    double Voltage=atGetDouble(ElemData,"Voltage");
    double Energy=atGetDouble(ElemData,"Energy");
    double Frequency=atGetDouble(ElemData,"Frequency");
    double TimeLag=atGetOptionalDouble(ElemData,"TimeLag",0);
    if (mxGetM(prhs[1]) != 6)
      mexErrMsgIdAndTxt("AT:WrongArg","Second argument must be a 6 x N matrix");
    /* ALLOCATE memory for the output array of the same size as the input  */
    plhs[0] = mxDuplicateArray(prhs[1]);
    r_in = mxGetDoubles(plhs[0]);
    CavityPass(r_in,Length,Voltage/Energy,Frequency,TimeLag,num_particles);
  }
  else if (nrhs == 0) {   /* return list of required fields */
    plhs[0] = mxCreateCellMatrix(4,1);
    mxSetCell(plhs[0],0,mxCreateString("Length"));
    mxSetCell(plhs[0],1,mxCreateString("Voltage"));
    mxSetCell(plhs[0],2,mxCreateString("Energy"));
    mxSetCell(plhs[0],3,mxCreateString("Frequency"));
    if (nlhs>1) /* optional fields */
      {   plhs[1] = mxCreateCellMatrix(1,1);
	mxSetCell(plhs[1],0,mxCreateString("TimeLag"));
      }
  }
  else {
    mexErrMsgIdAndTxt("AT:WrongArg","Needs 0 or 2 arguments");
  }
}
#endif /* MATLAB_MEX_FILE */
