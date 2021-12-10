/* findmpoleraddifmatrix.c

   mex-function to calculate radiation diffusion matrix B defined in [2] 
   for multipole elements in MATLAB Accelerator Toolbox
   A.Terebilo 8/14/00

   References
   [1] M.Sands 'The Physics of Electron Storage Rings
   [2] Ohmi, Kirata, Oide 'From the beam-envelope matrix to synchrotron
   radiation integrals', Phys.Rev.E  Vol.49 p.751 (1994)
*/

#include "tracy-2.cc"

/* Fourth order-symplectic integrator constants */

#define DRIFT1    0.6756035959798286638
#define DRIFT2   -0.1756035959798286639
#define KICK1     1.351207191959657328
#define KICK2    -1.702414383919314656

/* Physical constants used in the calculations */

#define TWOPI	  6.28318530717959
#define CGAMMA 	  8.846056192e-05 	/* [m]/[GeV^3] Ref[1] (4.1)      */
#define M0C2      5.10999060e5		/* Electron rest mass [eV]       */
#define LAMBDABAR 3.86159323e-13	/* Compton wavelength/2pi [m]    */
#define CER   	  2.81794092e-15	/* Classical electron radius [m] */
#define CU        1.323094366892892	/* 55/(24*sqrt(3)) factor        */


static void edgefringeB(double* r, double *B, double inv_rho, double edge_angle,
			double fint, double gap)
{
  int    m;
  double fx, fy, psi;

  if (inv_rho<=0.0) return; /* Skip if not a bending element*/
    
  fx = inv_rho*tan(edge_angle);
  psi = inv_rho*gap*fint*(1+pow(sin(edge_angle),2))/cos(edge_angle);
  if(fint >0 && gap >0)
    fy = inv_rho*tan(edge_angle-psi/(1+r[4]));
  else
    fy = fx;
        
  /*  Propagate B */

  for(m=0;m<6;m++)
    {
      B[1+6*m] += fx*B[6*m];
      B[3+6*m] -= fy*B[2+6*m];
    }
  if(fint >0 && gap >0)
    for(m=0;m<6;m++)    
      B[3+6*m] -= B[4+6*m]*r[2]*
	(inv_rho*inv_rho+fy*fy)*psi/pow((1+r[4]),2)/inv_rho;
	

  for(m=0;m<6;m++)
    {
      B[m+6*1] += fx*B[m+6*0];
      B[m+6*3] -= fy*B[m+6*2];
    }
  if(fint >0 && gap >0)
    for(m=0;m<6;m++)
      B[m+6*3] -= B[m+6*4]*r[2]*
	(inv_rho*inv_rho+fy*fy)*psi/pow((1+r[4]),2)/inv_rho;

  /* Propagate particle */
  r[1]+=r[0]*fx;
  r[3]-=r[2]*fy;

}

static void thinkickM(double* orbit_in, double* A, double* B, double L,
		      double irho, int max_order, double *M66)
/* Calculate the symplectic (no radiation) transfer matrix of a 
   thin multipole kick near the entrance point orbit_in
   For elements with straight coordinate system irho = 0
   For curved elements the B polynomial (PolynomB in MATLAB) 
   MUST NOT include the guide field  By0 = irho * E0 /(c*e)

   M is a (*double) pointer to a preallocated 1-dimentional array 
   of 36 elements of matrix M arranged column-by-column 
*/
{
  int m,n;
   
  double ReSumNTemp;
  double ImSumN = max_order*A[max_order];
  double ReSumN = max_order*B[max_order];	

  /* Recursively calculate the derivatives
     ReSumN = (irho/B0)*Re(d(By + iBx)/dx)
     ImSumN = (irho/B0)*Im(d(By + iBx)/dy)
  */
  for(n=max_order-1;n>0;n--)
    {	ReSumNTemp = (ReSumN*orbit_in[0] - ImSumN*orbit_in[2]) + n*B[n];
      ImSumN = ImSumN*orbit_in[0] +  ReSumN*orbit_in[2] + n*A[n];
      ReSumN = ReSumNTemp;
    }

  /* Initialize M66 to a 6-by-6 identity matrix */
  for(m=0;m<36;m++)
    M66[m]= 0;
  /* Set diagonal elements to 1 */
  for(m=0;m<6;m++)
    M66[m*7] = 1;
	
  /* The relationship between indexes when a 6-by-6 matrix is 
     represented in MATLAB as one-dimentional array containing
     36 elements arranged column-by-column is
     [i][j] <---> [i+6*j] 
  */

  M66[1]   = -L*ReSumN;				/* [1][0] */
  M66[13]  =  L*ImSumN;				/* [1][2] */
  M66[3]   =  L*ImSumN;				/* [3][0] */
  M66[15]  =  L*ReSumN;				/* [3][2] */
  M66[25]  =  L*irho;				/* [1][4] */
  M66[1]  += -L*irho*irho;			/* [1][0] */
  M66[5]   =  L*irho;				/* [5][0] */

}



static void thinkickB(double* orbit_in, double* A, double* B, double L,
		      double irho, int max_order, double E0, double *B66)

/* Calculate Ohmi's diffusion matrix of a thin multipole  element 
   For elements with straight coordinate system irho = 0
   For curved elements the B polynomial (PolynomB in MATLAB) 
   MUST NOT include the guide field  By0 = irho * E0 /(c*e)
   The result is stored in a preallocated 1-dimentional array B66
   of 36 elements of matrix B arranged column-by-column
*/

{
  double BB,B2P,B3P;
  int    i;
  double p_norm = 1/(1+orbit_in[4]);
  double p_norm2 = sqr(p_norm);
  double ImSum = A[max_order];
  double ReSum = B[max_order];	
  double ReSumTemp;
	
  /* recursively calculate the local transvrese magnetic field
     ReSum = irho*By/B0
     ImSum = irho*Bx/B0
  */

  for(i=max_order-1;i>=0;i--)
    {	ReSumTemp = ReSum*orbit_in[0] - ImSum*orbit_in[2] + B[i];
      ImSum = ImSum*orbit_in[0] +  ReSum*orbit_in[2] + A[i];
      ReSum = ReSumTemp;
    }
	
	
  /* calculate |B x n|^3 - the third power of the B field component 
     orthogonal to the normalized velocity vector n
  */
  B2P = B2perp(ImSum, ReSum +irho, irho, orbit_in[0] , orbit_in[1]*p_norm , 
	       orbit_in[2] , orbit_in[3]*p_norm );
  B3P = B2P*sqrt(B2P);

  BB =
    CU * CER * LAMBDABAR *  pow(E0/M0C2,5) * L * B3P * sqr(sqr(1+orbit_in[4]))*
    (1+orbit_in[0]*irho + (sqr(orbit_in[1])+sqr(orbit_in[3]))*p_norm2/2);

	
  /* When a 6-by-6 matrix is represented in MATLAB as one-dimentional 
     array containing 36 elements arranged column-by-column,
     the relationship between indexes  is
     [i][j] <---> [i+6*j] 

  */
	
  /* initialize B66 to 0 */
  for(i=0;i<36;i++)
    B66[i] = 0;
	
  /* Populate B66 */
  B66[7]      = BB*sqr(orbit_in[1])*p_norm2;
  B66[19]     = BB*orbit_in[1]*orbit_in[3]*p_norm2;
  B66[9]      = B66[19];
  B66[21]     = BB*sqr(orbit_in[3])*p_norm2;
  B66[10]     = BB*orbit_in[1]*p_norm;
  B66[25]     = B66[10];
  B66[22]     = BB*orbit_in[3]*p_norm;
  B66[27]     = B66[22];
  B66[28]     = BB;
}


static void drift_propagateB(double *orb_in, double L,  double *B)
{	/* Propagate cumulative Ohmi's diffusion matrix B through a drift
	   B is a (*double) pointer to 1-dimentional array 
	   containing 36 elements of matrix elements arranged column-by-column
	   as in MATLAB representation 

	   The relationship between indexes when a 6-by-6 matrix is 
	   represented in MATLAB as one-dimentional array containing
	   36 elements arranged column-by-column is
	   [i][j] <---> [i+6*j] 
	*/
		
  int m;
		
  double DRIFTMAT[36];
  for (m=0;m<36;m++) DRIFTMAT[m] = 0.0;
  /* Set diagonal elements to 1	*/
  for (m=0;m<6;m++) DRIFTMAT[m*7] = 1.0;

  DRIFTMAT[6]  =  L/(1+orb_in[4]);
  DRIFTMAT[20] =  DRIFTMAT[6];
  DRIFTMAT[24] = -L*orb_in[1]/sqr(1+orb_in[4]);
  DRIFTMAT[26] = -L*orb_in[3]/sqr(1+orb_in[4]);
  DRIFTMAT[11] =  L*orb_in[1]/sqr(1+orb_in[4]);
  DRIFTMAT[23] =  L*orb_in[3]/sqr(1+orb_in[4]);	
  DRIFTMAT[29] =
    -L*(sqr(orb_in[1])+sqr(orb_in[3]))/((1+orb_in[4])*sqr(1+orb_in[4]));

  ATsandwichmmt(DRIFTMAT,B);
}


static void FindElemB(double *orbit_in, double le, double irho, double *A,
		      double *B, double *T1, double* T2,double *R1, double *R2,
		      double entrance_angle, double exit_angle,
		      double fringe_int1, double fringe_int2, double full_gap,
		      int max_order, int num_int_steps,
		      double E0, double *BDIFF)

{	/* Find Ohmi's diffusion matrix BDIFF of a thick multipole
	   BDIFF - cumulative Ohmi's diffusion is initialized to 0
	   BDIFF is preallocated 1 dimensional array to store 6-by-6 matrix 
	   columnwise
	*/
	
  int    m;	
  double MKICK[36], BKICK[36];

  /* 4-th order symplectic integrator constants */
  double SL, L1, L2, K1, K2;

  SL = le/num_int_steps;
  L1 = SL*DRIFT1;
  L2 = SL*DRIFT2;
  K1 = SL*KICK1;
  K2 = SL*KICK2;
	
	
  /* Allocate memory for thin kick matrix MKICK
     and a diffusion matrix BKICK
  */
  for(m=0; m < 6; m++)
    {
      MKICK[m] = 0;
      BKICK[m] = 0;
    }
	
  /* Transform orbit to a local coordinate system of an element
     BDIFF stays zero	*/
  if (T1) ATaddvv(orbit_in, T1);	
  if (R1) ATmultmv(orbit_in, R1);
    
  /* Propagate orbit_in and BDIFF through the entrance edge */
  edgefringeB(orbit_in, BDIFF, irho, entrance_angle, fringe_int1, full_gap);

  /* Propagate orbit_in and BDIFF through a 4-th orderintegrator */

  for(m=0; m < num_int_steps; m++) /* Loop over slices	*/			
    {
      drift_propagateB(orbit_in,L1, BDIFF);
      Drift(orbit_in, L1, false);
				
      thinkickM(orbit_in, A, B, K1, irho, max_order, MKICK);
      thinkickB(orbit_in, A, B, K1, irho, max_order, E0, BKICK);
      ATsandwichmmt(MKICK, BDIFF);
      ATaddmm(BKICK, BDIFF);
      thin_kick(orbit_in, A, B, K1, irho, max_order, E0, true);
		
      drift_propagateB(orbit_in, L2, BDIFF);
      Drift(orbit_in, L2, false);
				
      thinkickM(orbit_in, A, B, K2, irho, max_order, MKICK);
      thinkickB(orbit_in, A, B, K2, irho, max_order, E0, BKICK);
      ATsandwichmmt(MKICK, BDIFF);
      ATaddmm(BKICK, BDIFF);
      thin_kick(orbit_in, A, B, K2, irho, max_order, E0, true);
	
      drift_propagateB(orbit_in,L2, BDIFF);
      Drift(orbit_in, L2, false);
				
      thinkickM(orbit_in, A, B, K1, irho, max_order, MKICK);
      thinkickB(orbit_in, A, B, K1, irho, max_order, E0, BKICK);
      ATsandwichmmt(MKICK, BDIFF);
      ATaddmm(BKICK, BDIFF);
      thin_kick(orbit_in, A, B, K1, irho, max_order, E0, true);

      drift_propagateB(orbit_in, L1, BDIFF);
      Drift(orbit_in, L1, false);
    }  
		
  edgefringeB(orbit_in, BDIFF, irho, exit_angle, fringe_int2, full_gap);

  if (R2) ATmultmv(orbit_in,R2);	
  if (T2) ATaddvv(orbit_in,T2);	
}


static double *diffmatrix(const PyObject *ElemData, double *orb, double energy,
			  double *bdiff)
{
  double
    Length, ba, irho, FringeInt1,  FringeInt2, FullGap, *PolynomA,
    *PolynomB;  

  int    MaxOrder, NumIntSteps;
  double EntranceAngle, ExitAngle;
  double *R1, *R2, *T1, *T2;

  /* Required fields */
  Length=atGetOptionalDouble(ElemData,"Length", 0.0);
  /* If ELEMENT has a zero length, return zeros matrix end exit */
  if (Length == 0.0) return bdiff;
  PolynomA=atGetOptionalDoubleArray(ElemData, (char*)"PolynomA");
  PolynomB=atGetOptionalDoubleArray(ElemData, (char*)"PolynomB");
  /* If the ELEMENT does not have PolynomA and PolynomB return zero matrix and
     exit */
  if (PolynomA == NULL ||  PolynomB == NULL) return bdiff;

  MaxOrder=atGetLong(ElemData,"MaxOrder"); check_error();
  NumIntSteps=atGetLong(ElemData,"NumIntSteps"); check_error();
  if (atIsNaN(energy)) {
    energy=atGetDouble(ElemData,"Energy"); check_error();
  }

  /* Optional fields */
  ba=atGetOptionalDouble(ElemData,"BendingAngle",0.0);
  irho = ba/Length;
  EntranceAngle=atGetOptionalDouble(ElemData,"EntranceAngle",0.0);
  ExitAngle=atGetOptionalDouble(ElemData,"ExitAngle",0.0);
  FullGap=atGetOptionalDouble(ElemData,"FullGap",0); check_error();
  FringeInt1=atGetOptionalDouble(ElemData,"FringeInt1",0); check_error();
  FringeInt2=atGetOptionalDouble(ElemData,"FringeInt2",0); check_error();
  R1=atGetOptionalDoubleArray(ElemData, (char*)"R1"); check_error();
  R2=atGetOptionalDoubleArray(ElemData, (char*)"R2"); check_error();
  T1=atGetOptionalDoubleArray(ElemData, (char*)"T1"); check_error();
  T2=atGetOptionalDoubleArray(ElemData, (char*)"T2"); check_error();

  FindElemB(orb, Length, irho, PolynomA, PolynomB, 
            T1, T2, R1, R2,
            EntranceAngle, 	ExitAngle,
            FringeInt1, FringeInt2, FullGap,
            MaxOrder, NumIntSteps, energy, bdiff);
  return bdiff;
}


#if defined(PYAT)

#define MODULE_NAME diffmatrix
#define MODULE_DESCR "Computation of the radiation diffusion matrix"

static PyObject *compute_diffmatrix(PyObject *self, PyObject *args) {
  PyObject      *pyElem, *pyMatrix;
  PyArrayObject *pyOrbit;
  double        *orb0, *bdiff, *retval;
  double        orb[6];
  double        energy;
  npy_intp      outdims[2] = {6, 6};
  int           i;

  if (!PyArg_ParseTuple
      (args, "OO!d", &pyElem, &PyArray_Type, &pyOrbit, &energy)) {
    return NULL;
  }
  if (PyArray_DIM(pyOrbit,0) != 6) {
    PyErr_SetString(PyExc_ValueError, "Orbit is not a (6,) array");
    return NULL;
  }
  if (PyArray_TYPE(pyOrbit) != NPY_DOUBLE) {
    PyErr_SetString(PyExc_ValueError, "Orbit is not a double array");
    return NULL;
  }
  if ((PyArray_FLAGS(pyOrbit) & NPY_ARRAY_FARRAY_RO) != NPY_ARRAY_FARRAY_RO) {
    PyErr_SetString(PyExc_ValueError, "Orbit is not Fortran-aligned");
    return NULL;
  }

  // C -> C++.
  orb0 = static_cast<double*>(PyArray_DATA(pyOrbit));
  /* make local copy of the input closed orbit vector */
  for (i=0;i<6;i++) orb[i] = orb0[i];

  /* ALLOCATE memory for the output array */
  pyMatrix = PyArray_ZEROS(2, outdims, NPY_DOUBLE, 1);
  bdiff = static_cast<double*>(PyArray_DATA((PyArrayObject *)pyMatrix));

  retval = diffmatrix(pyElem, orb, energy, bdiff);
  if (retval == NULL) {
    Py_DECREF(pyMatrix);
    return NULL;
  }
  return pyMatrix;
}

static PyMethodDef AtMethods[] =
  { {"find_mpole_raddiff_matrix",
     (PyCFunction)compute_diffmatrix, METH_VARARGS,
     PyDoc_STR
     ("diffmatrix=find_mpole_raddiff_matrix(element, orbit, energy)\n\n"
      "element:   element structure with field names consistent with\n"
      "           a multipole transverse field model.\n"
      "orbit:     (6,) vector of the closed orbit at the entrance\n"
      "energy:    ring energy\n\n"
      
      "calculate radiation diffusion matrix B defined in [2]\n"
      "for multipole elements in MATLAB Accelerator Toolbox\n"
      "A.Terebilo 8/14/00\n\n"
      
      "References\n"
      "[1] M.Sands 'The Physics of Electron Storage Rings\n"
      "[2] Ohmi, Kirata, Oide 'From the beam-envelope matrix to synchrotron\n"
      "radiation integrals', Phys.Rev.E  Vol.49 p.751 (1994)\n"
      )},
    {NULL, NULL, 0, NULL}        /* Sentinel */
};

PyMODINIT_FUNC MOD_INIT(MODULE_NAME)
{

#if PY_MAJOR_VERSION >= 3
  static struct PyModuleDef moduledef =
    { PyModuleDef_HEAD_INIT,
      STR(MODULE_NAME),             /* m_name */
      PyDoc_STR(MODULE_DESCR),      /* m_doc */
      -1,                           /* m_size */
      AtMethods,                    /* m_methods */
      NULL,                         /* m_reload */
      NULL,                         /* m_traverse */
      NULL,                         /* m_clear */
      NULL,                         /* m_free */
    };

  PyObject *m = PyModule_Create(&moduledef);
#else
  PyObject *m = Py_InitModule3(STR(MODULE_NAME), AtMethods,
			       MODULE_DESCR);
#endif
  if (m == NULL) return MOD_ERROR_VAL;
  import_array();
  return MOD_SUCCESS_VAL(m);
}

#endif /*PYAT*/
