#include <math.h>

#include <armadillo>

#define PS_DIM 6

#define TWOPI  6.28318530717959
#define C0     2.99792458e8


#define sqr(x)  ((x)*(x))
#define cube(x) ((x)*(x)*(x))

#define HOMmax  21

#define DRIFT1    0.6756035959798286638
#define DRIFT2   -0.1756035959798286639
#define KICK1     1.351207191959657328
#define KICK2    -1.702414383919314656


double C_u, C_gamma, C_q, cl_rad, q_fluct;


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

struct elem_id { };

struct elem_drift { };

/* struct elem_mpole */
/* { */
/*   int */
/*     MaxOrder, */
/*     NumIntSteps, */
/*     FringeQuadEntrance, */
/*     FringeQuadExit; */
/*   double */
/*     *PolynomA, */
/*     *PolynomB, */
/*     *fringeIntM0, */
/*     *fringeIntP0, */
/*     *KickAngle; */
/* }; */

struct elem_mpole {
  int
    MaxOrder,
    NumIntSteps,
    // Dipole.
    FringeBendEntrance,
    FringeBendExit,
    // Quadrupole or gradient dipole.
    FringeQuadEntrance,
    FringeQuadExit;
  double
    irho,               // 1/rho; non zero for dipoles.
    // Dipole.
    BendingAngle,       // [rad].
    EntranceAngle,      // [rad].
    ExitAngle,          // [rad].
    FringeInt1,
    FringeInt2,
    FullGap,
    // General multipole.
    *PolynomA,
    *PolynomB,
    *fringeIntM0,
    *fringeIntP0,
      *KickAngle;
};

struct elem_cav {
  double
    Voltage,            // V0/E0.
    Energy,
    Frequency,
    TimeLag;
};

struct elem_type {
  double
    Length,
    *R1,
    *R2,
    *T1,
    *T2,
    *EApertures,
    *RApertures;
  union {
    elem_id    *id_ptr;
    elem_drift *drift_ptr;
    elem_mpole *mpole_ptr;
    /* elem_bend  *bend_ptr; */
    elem_cav   *cav_ptr;
  };
};
