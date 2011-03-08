#include <TParticle.h>

#include "AliLog.h"
#include "AliStack.h"
#include "AliESDtrack.h"
#include "AliTOFT0maker.h"

#include "AliRsnMonitorTrack.h"

ClassImp(AliRsnMonitorTrack)

AliRsnMonitorTrack::AliRsnMonitorTrack() :
   fCutsPassed(kFALSE),
   fPrim(kFALSE),
   fPDG(0),
   fPDGM(0),
   fMother(-1),
   fStatus(0),
   fLength(0.0),
   fCharge(0),
   fITSsa(kFALSE),
   fTOFok(kFALSE),
   fITSchi2(1E10),
   fITSsignal(0.0),
   fTPCcount(0),
   fTPCchi2(1E10),
   fTPCsignal(1E10),
   fTOFsignal(0.0)
{
//
// Unique constructor
//
   
   Int_t k;

   for (k = 0; k < 2; k++) fDCA[k] = 0.0;
   for (k = 0; k < 6; k++) fITSmap[k] = 0;
   for (k = 0; k < AliPID::kSPECIES; k++) {
      fITSnsigma[k] = 0.0;
      fTPCnsigma[k] = 0.0;
      fTOFref   [k] = 0.0;
      fTOFsigma [k] = 0.0;
   }
   for (k = 0; k < 3; k++) {
      fPsim[k] = 0.0;
      fPrec[k] = 0.0;
      fPtpc[k] = 0.0;
   }
}

AliRsnMonitorTrack::AliRsnMonitorTrack(const AliRsnMonitorTrack& copy) :
   TObject(copy),
   fCutsPassed(copy.fCutsPassed),
   fPrim(copy.fPrim),
   fPDG(copy.fPDG),
   fPDGM(copy.fPDGM),
   fMother(copy.fMother),
   fStatus(copy.fStatus),
   fLength(copy.fLength),
   fCharge(copy.fCharge),
   fITSsa(copy.fITSsa),
   fTOFok(copy.fTOFok),
   fITSchi2(copy.fITSchi2),
   fITSsignal(copy.fITSsignal),
   fTPCcount(copy.fTPCcount),
   fTPCchi2(copy.fTPCchi2),
   fTPCsignal(copy.fTPCsignal),
   fTOFsignal(copy.fTOFsignal)
{
//
// Copy constructor
//

   Int_t k;

   for (k = 0; k < 2; k++) fDCA[k] = copy.fDCA[k];
   for (k = 0; k < 6; k++) fITSmap[k] = copy.fITSmap[k];
   for (k = 0; k < AliPID::kSPECIES; k++) {
      fITSnsigma[k] = copy.fITSnsigma[k];
      fTPCnsigma[k] = copy.fTPCnsigma[k];
      fTOFref   [k] = copy.fTOFref   [k];
      fTOFsigma [k] = copy.fTOFsigma [k];
   }
   for (k = 0; k < 3; k++) {
      fPsim[k] = copy.fPsim[k];
      fPrec[k] = copy.fPrec[k];
      fPtpc[k] = copy.fPtpc[k];
   }
}

//_____________________________________________________________________________________________
void AliRsnMonitorTrack::Reset()
{
//
// Generic reset method, to set all fields to meaningless values
//

   Int_t k;

   fCutsPassed = kFALSE;
   fITSsa      = kFALSE;
   fTOFok      = kFALSE;
   fPrim       = kFALSE;
   fPDG        = 0;
   fMother     = -1;
   fPDGM       = 0;
   fStatus     = 0;
   fLength     = 0.0;
   fCharge     = 0;

   for (k = 0; k < 2; k++) fDCA[k] = 1E10;

   for (k = 0; k < 6; k++) fITSmap[k] = kFALSE;
   fITSchi2   = 1E10;
   fITSsignal = 0.0;

   fTPCchi2   = 1E10;
   fTPCsignal = 1E10;
   fTPCcount  =    0;

   fTOFsignal = 1E10;

   for (k = 0; k < AliPID::kSPECIES; k++) {
      fITSnsigma[k] = 1E10;
      fTPCnsigma[k] = 1E10;
      fTOFref   [k] = 1E10;
      fTOFsigma [k] = 1E10;
   }

   for (k = 0; k < 3; k++) fPsim[k] = fPrec[k] = fPtpc[k] = 0.0;
}

//__________________________________________________________________________________________________
Bool_t AliRsnMonitorTrack::AdoptMC(Int_t label, AliStack *stack)
{
//
// Get info from MC for a given track in the stack
//

   if (!stack) return  kFALSE;

   Int_t nPart = stack->GetNtrack();
   if (label < 0 || label > nPart) {
      AliError(Form("Label = %d -- MAX = %d", label, nPart));
      return kFALSE;
   }

   TParticle *mc = stack->Particle(label);
   if (!mc) return kFALSE;

   // 'direct' data
   fPDG     = (Int_t)mc->GetPdgCode();
   fMother  = (Int_t)mc->GetFirstMother();
   fPrim    = (Bool_t)stack->IsPhysicalPrimary(label);
   fPDGM    = 0;
   fPsim[0] = mc->Px();
   fPsim[1] = mc->Py();
   fPsim[2] = mc->Pz();

   // assign mother (if any)
   if (fMother >= 0 && fMother < nPart) {
      TParticle *m = stack->Particle(fMother);
      if (m) fPDGM = (Int_t)TMath::Abs(m->GetPdgCode());
   }

   return kTRUE;
}
