#include "cfortran.h"
#include "Rtypes.h"

#include "Fdimpar.h"

extern "C" {
//*$ create trackr.add
//*copy trackr
//*                                                                      *
//*=== trackr ===========================================================*
//*                                                                      *
//*----------------------------------------------------------------------*
//*                                                                      *
//*     tracks recording       by  alfredo ferrari, infn - milan         *
//*                                                                      *
//*     last change    31 january 2001    by   alfredo ferrari           *
//*                                                                      *
//*            included in :                                             *
//*                          electr                                      *
//*                          emfsco                                      *
//*                          kaskad (new version)                        *
//*                          kashea                                      *
//*                          kasneu                                      *
//*                          geoden (new version)                        *
//*                          mageas                                      *
//*                          magmov                                      *
//*                          magnew                                      *
//*                          move                                        *
//*                          photon                                      *
//*                          usrsco                                      *
//*                                                                      *
//*          ntrack = number of track segments                           *
//*          mtrack = number of energy deposition events along the track *
//*   0 < i < ntrack                                                     *
//*          xtrack = end x-point of the ith track segment               *
//*          ytrack = end y-point of the ith track segment               *
//*          ztrack = end z-point of the ith track segment               *
//*   1 < i < ntrack                                                     *
//*          ttrack = length of the ith track segment                    *
//*   1 < j < mtrack                                                     *
//*          dtrack = energy deposition of the jth deposition event      *
//*          dptrck = momentum loss of the jth deposition event          *
//*                                                                      *
//*          jtrack = identity number of the particle                    *
//*          etrack = total energy of the particle                       *
//*          ptrack = momentum of the particle (not always defined, if   *
//*                 < 0 must be obtained from etrack)                    *
//*      cx,y,ztrck = direction cosines of the current particle          *
//*      cx,y,ztrpl = polarization cosines of the current particle       *
//*          wtrack = weight of the particle                             *
//*          wscrng = scoring weight: it can differ from wtrack if some  *
//*                   biasing techniques are used (for example inelastic *
//*                   interaction length biasing)                        *
//*          ctrack = total curved path                                  *
//*          cmtrck = cumulative curved path since particle birth        *
//*          zfftrk = <z_eff> of the particle                            *
//*          zfrttk = actual z_eff of the particle                       *
//*          atrack = age of the particle                                *
//*          akshrt = kshrt amplitude for k0/k0bar                       *
//*          aklong = klong amplitude for k0/k0bar                       *
//*          wninou = neutron algebraic balance of interactions (both    *
//*                   for "high" energy particles and "low" energy       *
//*                   neutrons)                                          *
//*          spausr = user defined spare variables for the current       *
//*                   particle                                           *
//*          sttrck = macroscopic total cross section for low energy     *
//*                   neutron collisions                                 *
//*          satrck = macroscopic absorption cross section for low energy*
//*                   neutron collisions (it can be negative for pnab>1) *
//*          ktrack = if > 0 neutron group of the particle (neutron)     *
//*                                                                      *
//*          ntrack > 0, mtrack > 0 : energy loss distributed along the  *
//*                                   track                              *
//*          ntrack > 0, mtrack = 0 : no energy loss along the track     *
//*          ntrack = 0, mtrack = 0 : local energy deposition (the       *
//*                                   value and the point are not re-    *
//*                                   corded in trackr)                  *
//*          mmtrck = flag recording the material index for low energy   *
//*                   neutron collisions                                 *
//*          lt1trk = initial lattice cell of the current track          *
//*                  (or lattice cell for a point energy deposition)     *
//*          lt2trk = final   lattice cell of the current track          *
//*          ihspnt = current geometry history pointer (not set if -1)   *
//*          ltrack = flag recording the generation number               *
//*          llouse = user defined flag for the current particle         *
//*          ispusr = user defined spare flags for the current particle  *
//*          lfsssc = logical flag for inelastic interactions ending with*
//*                   fission (used also for low energy neutrons)        *
//*                                                                      *
//*----------------------------------------------------------------------*
//

//
// TFluka specific:
// ispusr[mkbmx2 - 1] : track index in vmcstack
// ispusr[mkbmx2 - 2] : flag for "interrupted" track
//
    
const Int_t mxtrck = 2500;

typedef struct {
   Double_t xtrack[mxtrck+1];
   Double_t ytrack[mxtrck+1];
   Double_t ztrack[mxtrck+1];
   Double_t ttrack[mxtrck];
   Double_t dtrack[mxtrck];
   Double_t dptrck[mxtrck][3];
   Double_t etrack;
   Double_t ptrack;
   Double_t cxtrck;
   Double_t cytrck;
   Double_t cztrck;
   Double_t wtrack;
   Double_t cxtrpl;
   Double_t cytrpl;
   Double_t cztrpl;
   Double_t zfftrk;
   Double_t zfrttk;
   Double_t atrack;
   Double_t ctrack;
   Double_t cmtrck;
   Double_t akshrt;
   Double_t aklong;
   Double_t wscrng;
   Double_t wninou;
   Double_t spausr[mkbmx1];
   Double_t sttrck;
   Double_t satrck;
   Int_t    ntrack;
   Int_t    mtrack;
   Int_t    jtrack;
   Int_t    ktrack;
   Int_t    mmtrck;
   Int_t    lt1trk;
   Int_t    lt2trk;
   Int_t    ihspnt;
   Int_t    ltrack;
   Int_t    llouse;
   Int_t    ispusr[mkbmx2];
   Int_t    lfsssc;
   Int_t    lpkill;
} trackrCommon;
#define TRACKR COMMON_BLOCK(TRACKR,trackr)
COMMON_BLOCK_DEF(trackrCommon,TRACKR);
//static union { Double_t spause; Double_t spausr[0];};
//static union { Int_t    ispuse; Int_t    ispusr[0];};
}
