#include <Riostream.h>

#include "AliVTrack.h"
#include "AliESDtrack.h"
#include "AliAODTrack.h"

#include "AliPIDResponse.h"
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"


#include "AliRsnValuePID.h"

ClassImp(AliRsnValuePID)

//__________________________________________________________________________________________________
AliRsnValuePID::AliRsnValuePID() :
   AliRsnValue(),
   fSpecies(AliPID::kUnknown),
   fValuePID(kValues),
   fPID(0x0)
{
//
// Dummy constructor
//

   Int_t i;
   for (i = 0; i < AliPID::kSPECIES; i++) fTOFtimes[i] = fTOFsigma[i] = 0.0;
   
   SetTargetType(AliRsnTarget::kDaughter);
}

//__________________________________________________________________________________________________
AliRsnValuePID::AliRsnValuePID
(const char *name, EValuePID type, AliPID::EParticleType species, Int_t nbins, Double_t min, Double_t max) :
   AliRsnValue(name, nbins, min, max),
   fSpecies(species),
   fValuePID(type),
   fPID(0x0)
{
//
// Constructor 1 (fixed bins with number of bins, or no bins)
//

   Int_t i;
   for (i = 0; i < AliPID::kSPECIES; i++) fTOFtimes[i] = fTOFsigma[i] = 0.0;
   
   SetTargetType(AliRsnTarget::kDaughter);
}

//__________________________________________________________________________________________________
AliRsnValuePID::AliRsnValuePID
(const char *name, EValuePID type, AliPID::EParticleType species, Double_t min, Double_t max, Double_t step) :
   AliRsnValue(name, min, max, step),
   fSpecies(species),
   fValuePID(type),
   fPID(0x0)
{
//
// Constructor 2 (fixed bins with step)
//

   Int_t i;
   for (i = 0; i < AliPID::kSPECIES; i++) fTOFtimes[i] = fTOFsigma[i] = 0.0;
   
   SetTargetType(AliRsnTarget::kDaughter);
}

//__________________________________________________________________________________________________
AliRsnValuePID::AliRsnValuePID
(const char *name, EValuePID type, AliPID::EParticleType species, Int_t nbins, Double_t *array) :
   AliRsnValue(name, nbins, array),
   fSpecies(species),
   fValuePID(type),
   fPID(0x0)
{
//
// Constructor 3 (variable bins)
//
   
   Int_t i;
   for (i = 0; i < AliPID::kSPECIES; i++) fTOFtimes[i] = fTOFsigma[i] = 0.0;
   
   SetTargetType(AliRsnTarget::kDaughter);
}

//__________________________________________________________________________________________________
AliRsnValuePID::AliRsnValuePID(const AliRsnValuePID& copy) :
   AliRsnValue(copy),
   fSpecies(copy.fSpecies),
   fValuePID(copy.fValuePID),
   fPID(copy.fPID)
{
//
// Copy constructor
//

   Int_t i;
   for (i = 0; i < AliPID::kSPECIES; i++) fTOFtimes[i] = fTOFsigma[i] = 0.0;
}

//__________________________________________________________________________________________________
AliRsnValuePID& AliRsnValuePID::operator=(const AliRsnValuePID& copy)
{
//
// Assignment operator
//

   AliRsnValue::operator=(copy);
   
   fValuePID = copy.fValuePID;
   fSpecies = copy.fSpecies;
   fPID = copy.fPID;
   
   Int_t i;
   for (i = 0; i < AliPID::kSPECIES; i++) fTOFtimes[i] = fTOFsigma[i] = 0.0;
   
   return (*this);
}

//__________________________________________________________________________________________________
Bool_t AliRsnValuePID::Eval(TObject *object, Bool_t)
{
//
// Evaluation function
//

   if (fValuePID != kITSsignal && fValuePID != kTPCsignal && !fPID) InitializePID();
   if (!fPID) {
      AliError("PID not correctly initialized");
      fComputedValue = 0.0;
      return kFALSE;
   }
   if (!TargetOK(object)) return kFALSE;
   
   AliVTrack *vtrack = fDaughter->GetRefVtrack();
   
   switch (fValuePID) {
      case kITSsignal:
         fComputedValue = vtrack->GetITSsignal();
         return kTRUE;
      case kITSnsigma:
         fComputedValue = fPID->NumberOfSigmasITS(vtrack, fSpecies);
         return kTRUE;
      case kTPCsignal:
         fComputedValue = vtrack->GetTPCsignal();
         return kTRUE;
      case kTPCnsigma:
         fComputedValue = fPID->NumberOfSigmasTPC(vtrack, fSpecies);
         return kTRUE;
      case kTOFsignal:
         if (!TOFComputations(vtrack)) return kFALSE;
         fComputedValue = (vtrack->GetTOFsignal() - fPID->GetTOFResponse().GetStartTime(vtrack->P()));
         return kTRUE;
      case kTOFnsigma:
         if (!TOFComputations(vtrack)) return kFALSE;
         fComputedValue = fPID->NumberOfSigmasTOF(vtrack, fSpecies);
         return kTRUE;
      case kTOFtime:
         if (!TOFComputations(vtrack)) return kFALSE;
         fComputedValue = fTOFtimes[(Int_t)fSpecies];
         return kTRUE;
      case kTOFsigma:
         if (!TOFComputations(vtrack)) return kFALSE;
         fComputedValue = fTOFsigma[(Int_t)fSpecies];
         return kTRUE;
      default:
         AliError("Unrecognized option");
         fComputedValue = 0.0;
         return kFALSE;
   }
}

//__________________________________________________________________________________________________
void AliRsnValuePID::Print(Option_t *) const
{
//
// Printout
//

   AliRsnValue::Print();
}

//__________________________________________________________________________________________________
void AliRsnValuePID::InitializePID()
{
//
// Initialize PID object
//

   AliAnalysisManager   *man = AliAnalysisManager::GetAnalysisManager();
   AliInputEventHandler *inh = (AliInputEventHandler*)man->GetInputEventHandler();
   
   fPID = inh->GetPIDResponse();
}

//__________________________________________________________________________________________________
Bool_t AliRsnValuePID::TOFComputations(AliVTrack *vtrack)
{
//
// Make TOF computations
//

   if (vtrack->InheritsFrom(AliESDtrack::Class())) {
      AliESDtrack *track = (AliESDtrack*)vtrack;
      track->GetIntegratedTimes(fTOFtimes);
      Int_t i;
      for (i = 0; i < AliPID::kSPECIES; i++) {
         fTOFsigma[i] = fPID->GetTOFResponse().GetExpectedSigma(track->GetP(), fTOFtimes[i], AliPID::ParticleMass(i));
      }
      return kTRUE;
   } else if (vtrack->InheritsFrom(AliAODTrack::Class())) {
      AliAODTrack *track = (AliAODTrack*)vtrack;
      AliAODPid *pidObj = track->GetDetPid();
      if (!pidObj) return kFALSE;
      pidObj->GetIntegratedTimes(fTOFtimes);
      pidObj->GetTOFpidResolution(fTOFsigma);
      return kTRUE;
   } else {
      return kFALSE;
   }
}

