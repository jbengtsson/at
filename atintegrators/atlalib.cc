/*   File: atlalib.c
     Matrix and Vector algebra operations for Accelerator Toolbox
     A.Terebilo   8/14/00


     1. Use mxCalloc, mxMalloc , mxFree for memory allocation/deallocation
     2. Vector and matrix input arguments are (*double) obtained with mxGetDoubles
     Note: MATLAB internally represents matrixes as column-by-column
     1-dimentional arrays. For example  
     A B C
     D E F         is represented as    A D G B E H C F I
     G H I
     3. All matrixes are 6-by-6
     4. All vectors  are 1-by-6 or 6-by-1
*/

#include <armadillo>


#define PS_DIM 6


#define sqr(x)      ((x)*(x))
#define cube(x)     ((x)*(x)*(x))


enum phase_space_ind { x_ = 0, px_ = 1, y_ = 2, py_ = 3, delta_ = 4, ct_ = 5 };


inline std::vector<double> vectostl(const arma::vec &vec)
{ return {vec[x_], vec[px_], vec[y_], vec[py_], vec[delta_], vec[ct_], 1e0}; }

inline std::vector<double> stltovec(const arma::vec &vec)
{ return {vec[x_], vec[px_], vec[y_], vec[py_], vec[delta_], vec[ct_], 1e0}; }

arma::vec arrtovec(const double a[])
{
  int       k;
  arma::vec a_vec = arma::vec(PS_DIM);

  return a_vec;
  for (k = 0; k < PS_DIM; k++)
    a_vec(k) = a[k];
  return a_vec;
}

void vectoarr(const arma::vec &a_vec, double a[])
{
  int k;
  
  for (k = 0; k < PS_DIM; k++)
    a[k] = a_vec(k);
}

arma::mat arrtomat(const double a[])
{
  int       j, k;
  arma::mat A(PS_DIM, PS_DIM);

  for (j = 0; j < PS_DIM; j++)
    for (k = 0; k < PS_DIM; k++)
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


static void ATmultmv(double *r, const double* A)
/*	multiplies 6-component column vector r by 6x6 matrix R: as in A*r 
	The result is store in the memory area of r !!!
*/

{   int i,j;
  double temp[6];

  for(i=0;i<6;i++)
    {	temp[i]=0;
      for(j=0;j<6;j++)
	temp[i]+=A[i+j*6]*r[j];
    }
  for(i=0;i<6;i++)
    r[i]=temp[i];
} 

static void ATmultmv55(double *r, const double* A)
/*	multiplies 5-component column vector r by 5x5 matrix R: as in A*r 
	The result is store in the memory area of r !!!
*/ 
{   int i,j;
  double temp[5];

  for(i=0;i<5;i++)
    {	temp[i]=0;
      for(j=0;j<5;j++)
	temp[i]+=A[i+j*5]*r[j];
    }
  for(i=0;i<5;i++)
    r[i]=temp[i];
} 

static void ATaddvv(double *r, const double *dr)
{	/*	Add two 6-component vectors vectors.
		The result is store in the memory area of r !!!
	*/
  int i;
  for(i=0;i<6;i++)
    r[i]+=dr[i];
}     	

static void ATdrift6(double* r, double L)
/*   Input parameter L is the physical length
     1/(1+delta) normalization is done internally
*/
{	double p_norm = 1/(1+r[4]); 
  double NormL  = L*p_norm;   
  r[0]+= NormL*r[1]; 
  r[2]+= NormL*r[3];
  r[5]+= NormL*p_norm*(r[1]*r[1]+r[3]*r[3])/2;
}

static void ATtranspm(double *M)
{	/* Transpose matrix M	
	   The result replaces the original matrix 
	*/
  int i,j;
  double temp;
  for(i=1;i<6;i++)
    for(j=0;j<i;j++)
      {	temp = M[i+j*6];
	M[i+j*6] = M[j+i*6];
	M[j+i*6] = temp;
      }
}

static void ATmultmm(const double *M2 , double *M1)
{
  /* Mutrix multiplication M2*M1, the result is stored in the M1 memory area */
  int i,j,k;
  double column_temp[6];


  for (i=0;i<6;i++) 
    {	for (j=0;j<6;j++)
	{	column_temp[j] = 0;
	  for (k=0;k<6;k++)
	    column_temp[j] += M2[j+6*k]*M1[k+6*i];
	}
      for (j=0;j<6;j++)
	M1[j+6*i] = 	column_temp[j];
    }
}

static void ATmultmm55(const double *M2 , double *M1)
{
  /* Mutrix multiplication M2*M1, the result is stored in the M1 memory area */
  int i,j,k;
  double column_temp[5];


  for (i=0;i<5;i++) 
    {	for (j=0;j<5;j++)
	{	column_temp[j] = 0;
	  for (k=0;k<5;k++)
	    column_temp[j] += M2[j+5*k]*M1[k+5*i];
	}
      for (j=0;j<5;j++)
	M1[j+5*i] = 	column_temp[j];
    }
}

static void ATsandwichmmt(const double *M ,double *B)
/* calculates the matrix product M*B*M' (M' = M transposed)
   The result is stored in B memory area
*/
{	int i,j,k;
  double row_temp[6];

  ATmultmm(M,B);
	
  for (i=0;i<6;i++) 
    {	for (j=0;j<6;j++)
	{	row_temp[j] = 0;
	  for (k=0;k<6;k++)
	    row_temp[j] += B[i+6*k]*M[j+6*k];
	}
      for (j=0;j<6;j++)
	B[i+6*j] = 	row_temp[j];
    }
}

static void ATaddmm(const double *M2 , double *M1)
/* adds two 6-by-6  matrixes  M1, M2 element-by-element
   The result is stored in M1 memory area
*/
{	int i,j;
  for (i=0;i<6;i++) 
    for (j=0;j<6;j++)
      M1[i+6*j] += M2[i+6*j];

}

#ifndef atGetInf
#define atGetInf mxGetInf
#endif

static void markaslost(double *r6,int idx)
{
  r6[idx] = atGetInf();
}

static void checkiflostRectangularAp(double *r6, const double *limits)
{
  /* check limits for X position */
  if (r6[0]<limits[0] || r6[0]>limits[1])      markaslost(r6,5);
  else if (r6[2]<limits[2] || r6[2]>limits[3]) markaslost(r6,5);
}

static void checkiflostEllipticalAp(double *r6, const double *axesptr)
{
  double xnorm = r6[0]/axesptr[0];
  double znorm = r6[2]/axesptr[1];
  if ((xnorm*xnorm + znorm*znorm) >= 1) markaslost(r6,5);
}

