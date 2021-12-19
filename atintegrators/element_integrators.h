/* no guard macros on purpose */
/*
 * can the marco implementations simplified by a more consistent naming
 * of the internal functions ?
 */

ELEM_PASS(corr,        CorrectorPass,   "")
ELEM_PASS(id,          IdentityPass,    "")
ELEM_PASS(drift,       DriftPass,       "")
ELEM_PASS(cav,         CavityPass,      "")
ELEM_PASS(M66,         Matrix66Pass,    "")

/* these below require a bit of helper init functions */
ELEM_PASS(aper,        AperturePass,    "")
#if 1
ELEM_PASS(bend,        MpolePassNoRad, "BndMPoleSymplectic4Pass")
ELEM_PASS(bend_rad,    MpolePassRad,   "BndMPoleSymplectic4RadPass")
ELEM_PASS(bend_exact,  MpoleE2Pass,    "BndMPoleSymplectic4E2Pass")
ELEM_PASS(mpole,       MpolePass,      "StrMPoleSymplectic4Pass")
ELEM_PASS(mpole_rad,   MpolePassRad,   "StrMPoleSymplectic4RadPass")
ELEM_PASS(cbend,       CBendPass,      "BndStrMPoleSymplectic4Pass")
ELEM_PASS(H_exact,     HamPass,        "ExactHamiltonianPass")
/* these use */
ELEM_PASS(wig,         WigPass,        "GWigSymplecticPass")
ELEM_PASS(wig_rad,     WigPass,        "GWigSymplecticRadPass")
#endif
/* no guard macros on purpose */
