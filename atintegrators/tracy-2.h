
#define TWOPI     6.28318530717959

#define HOMmax    21

#define DRIFT1      0.6756035959798286638
#define DRIFT2     -0.1756035959798286639
#define KICK1       1.351207191959657328
#define KICK2      -1.702414383919314656

#define WHmax      200
#define GWIG_EPS   (1e-6)
#define second     2
#define fourth     4

#define atIsFinite isfinite
#define atIsNaN    isnan
#define atGetNaN() (NAN)
#define atGetInf() (INFINITY)
#define atMalloc   malloc
#define atCalloc   calloc
#define atFree     free

#define sqr(x)  ((x)*(x))
#define cube(x) ((x)*(x)*(x))


static const int
  Elem_Entrance =  1,
  Elem_Exit     = -1;

static const double
  q_e       = 1.602176462e-19, /* electron charge, [C] */
  m_e       = 9.10938188e-31,  /* electron mass, [kg] */
  clight    = 2.99792458e8,    /* speed of light [m/s] */
  r_e       = 2.817940285e-15, /* electron classic radius, [m] */
  XMC2      = 0.510998902e-03, /* mc^2 in GeV */
  PI        = 3.141592653589793238462643383279502884197e0,
  epsilon_o = 8.854187817e-12, /* Vacuum permittivity */
  C_gamma   = 8.846056192e-05; /* [m]/[GeV^3] Ref[1] (4.1) */

double C_u, C_q, cl_rad, q_fluct;


enum spatial_ind
  {
   X_ = 0,
   Y_,
   Z_
  };

enum phase_space_ind
  {
   x_ = 0,
    px_,
    y_,
    py_,
    delta_,
    ct_
  };

// Used by track_element.
enum element_type
  {
   drift = 0,
   dipole,
   multipole,
   marker
  };

enum elem_kind
  {
   corr_ = 0,
   ident_,
   aper_,
   drift_,
   mpole_, mpole_rad_,
   bend_, bend_rad_,
   bend_exact,
   cbend_,
   cav_,
   wig_, wig_rad,
   M66_,
   H_exact
  };

enum MpoleKind
  {
   All = 0,
   Dip,
   Quad,
   Sext,
   Oct,
   Dec5,
   Dodec
  };

struct elem_id { };

struct elem_ap {
  double *limits;
};

struct elem_drift { };

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
    Energy,
    irho,               // 1/rho; non zero for dipoles.
    // Dipole.
    BendingAngle,       // [rad].
    EntranceAngle,      // [rad].
    ExitAngle,          // [rad].
    FringeInt1,
    FringeInt2,
    FullGap,
    *PolynomA,          // General multipole.
    *PolynomB,          //        ''
    *fringeIntM0,
    *fringeIntP0,
    *KickAngle,
    X0ref,              // Cartesian bend.
    ByError,            //      ''
    RefDZ,              //      ''
    H1,
    H2;                 // E2 mod.
};

struct elem_cav {
  double
    Voltage,            // V0/E0.
    Energy,
    Frequency,          // [Hz].
    TimeLag;
};

struct elem_wig {
  int
    Pmethod,        /* Integration Method */
    PN,             /* Number of integration steps */
    Nw,             /* Number of periods */
    NHharm,         /* No. of horizontal harmonics */
    NVharm,         /* No. of vertical harmonics */

    // For radiation calculations.
    HSplitPole,
    VSplitPole;

  double
    E0,             /* Energy of ring, [GeV] */
    PB0,            /* B0 in [Tesla] */
    Lw,             /* Wiggler Period [m] */

    Zw,             /* Longitudinal variable [m] */
    Aw,             /* Wiggler parameter */
    HCw[WHmax],
    VCw[WHmax],
    HCw_raw[WHmax],
    VCw_raw[WHmax],
    Hkx[WHmax],
    Hky[WHmax],
    Hkz[WHmax],
    Htz[WHmax],
    Vkx[WHmax],
    Vky[WHmax],
    Vkz[WHmax],
    Vtz[WHmax],

    // For radiation calculations.
    zStartH,
    zStartV,  /* Start and end z coordinates of the wiggler field, which
		 are computed                                                 */
    zEndH,
    zEndV,    /* based on the phase of the first harmonic to get matched
		 dispersion.                                                  */
    srCoef,
    Po;       /* beta*gamma for reference particle */
};

struct elem_M66 {
  double *M66;
};

struct elem_corr {
  double *KickAngle;
};


#define FMAX 32

struct elem_H {
  int
    Type,            /* type is defined in track.h:
			- 0: drift
			- 1: dipole
			- 2: multipole
			- 3: marker                   */
    NumIntSteps,
    MaxOrder,
    MultipoleFringe; /* bool, whether to calculate multipole fringe */
  double
    *PolynomA,
    *PolynomB,
    F[FMAX],         // Temporary storage.
    BendingAngle,    /* required for bend */
    gK;              /* g * K, required for bend */
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
  elem_kind
    kind;
  union {
    elem_id    *id_ptr;
    elem_ap    *ap_ptr;
    elem_drift *drift_ptr;
    elem_corr  *corr_ptr;
    elem_mpole *mpole_ptr;
    elem_cav   *cav_ptr;
    elem_M66   *M66_ptr;
    elem_H     *H_ptr;
    elem_wig   *wig_ptr;
  };
};

// For HamPass.
struct lattice {
  elem_type *next;
  int       N;
};

