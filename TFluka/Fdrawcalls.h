#include "Rtypes.h"
//
// Structure saving parameters of Fluka procedures
// mgdraw(Int_t& icode, Int_t& mreg)
// bxdraw(Int_t& icode, Int_t& mreg, Int_t& newreg,
//        Double_t& xsco, Double_t& ysco, Double_t& zsco)
// eedraw(Int_t& icode)
// endraw(Int_t& icode, Int_t& mreg,
//        Double_t& rull, Double_t& xsco, Double_t& ysco, Double_t& zsco)
// sogdraw
// usdraw(Int_t& icode, Int_t& mreg,
//        Double_t& xsco, Double_t& ysco, Double_t& zsco)
//
//*---------------------- icode values ----------------------------------*
//*                                                                      *
//*    Icode = 1:  call from Kaskad                                      *
//*    Icode = 2:  call from Emfsco                                      *
//*    Icode = 3:  call from Kasneu                                      *
//*    Icode = 4:  call from Kashea                                      *
//*    Icode = 5:  call from Kasoph                                      *
//*                                                                      *
//*    Boundary- (X) crossing                                            *
//*    ======================                                            *
//*    Icode = 1x: call from Kaskad                                      *
//*    Icode = 19: boundary crossing                                     *
//*    Icode = 2x: call from Emfsco                                      *
//*    Icode = 29: boundary crossing                                     *
//*    Icode = 3x: call from Kasneu                                      *
//*    Icode = 39: boundary crossing                                     *
//*    Icode = 4x: call from Kashea                                      *
//*    Icode = 49: boundary crossing                                     *
//*    Icode = 59: call from Kasoph                                      *
//*    Icode = 59: boundary crossing                                     *
//*                                                                      *
//*    ENergy deposition DRAWing:                                        *
//*    =========================                                         *
//*                                                                      *
//*     Icode = 1x: call from Kaskad                                     *
//*             10: elastic interaction recoil                           *
//*             11: inelastic interaction recoil                         *
//*             12: stopping particle                                    *
//*             13: pseudo-neutron deposition                            *
//*             14: escape                                               *
//*             15: time kill                                            *
//*     Icode = 2x: call from Emfsco                                     *
//*             20: local energy deposition (i.e. photoelectric)         *
//*             21: below threshold, iarg=1                              *
//*             22: below threshold, iarg=2                              *
//*             23: escape                                               *
//*             24: time kill                                            *
//*     Icode = 3x: call from Kasneu                                     *
//*             30: target recoil                                        *
//*             31: below threshold                                      *
//*             32: escape                                               *
//*             33: time kill                                            *
//*     Icode = 4x: call from Kashea                                     *
//*             40: escape                                               *
//*             41: time kill                                            *
//*     Icode = 5x: call from Kasoph                                     *
//*             50: optical photon absorption                            *
//*             51: escape                                               *
//*             52: time kill                                            *
//*                                                                      *
//*     USer dependent DRAWing:                                          *
//*     ======================                                           *
//*                                                                      *
//*     Icode = 10x: call from Kaskad                                    *
//*             100: elastic   interaction secondaries                   *
//*             101: inelastic interaction secondaries                   *
//*             102: particle decay  secondaries                         *
//*             103: delta ray  generation secondaries                   *
//*             104: pair production secondaries                         *
//*             105: bremsstrahlung  secondaries                         *
//*     Icode = 20x: call from Emfsco                                    *
//*             208: bremsstrahlung secondaries                          *
//*             210: Moller secondaries                                  *
//*             212: Bhabha secondaries                                  *
//*             214: in-flight annihilation secondaries                  *
//*             215: annihilation at rest   secondaries                  *
//*             217: pair production        secondaries                  *
//*             219: Compton scattering     secondaries                  *
//*             221: photoelectric          secondaries                  *
//*             225: Rayleigh scattering    secondaries                  *
//*     Icode = 30x: call from Kasneu                                    *
//*             300: interaction secondaries                             *
//*     Icode = 40x: call from Kashea                                    *
//*             400: delta ray  generation secondaries                   *
//*     For all interactions secondaries are put on FINUC common (kp=1,np)
//*     but for KASHEA delta ray generation where only the secondary     *
//*     electron is present and stacked on STACK common for kp=lstack    *
//*                                                                      *
//*---------------------- end of icode values ---------------------------*
//

typedef struct {
   Int_t    icode;
   Int_t    mreg;
   Int_t    newreg;
   Double_t rull;
   Double_t xsco;
   Double_t ysco;
   Double_t zsco;
} Fdrawcalls;

