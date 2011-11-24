//
// Class AliRsnCutPID
//
// General implementation of a single cut strategy, which can be:
// - a value contained in a given interval  [--> IsBetween()   ]
// - a value equal to a given reference     [--> MatchesValue()]
//
// In all cases, the reference value(s) is (are) given as data members
// and each kind of cut requires a given value type (Int, UInt, Double),
// but the cut check procedure is then automatized and chosen thanks to
// an enumeration of the implemented cut types.
// At the end, the user (or any other point which uses this object) has
// to use the method IsSelected() to check if this cut has been passed.
//
// authors: Martin Vala (martin.vala@cern.ch)
//          Alberto Pulvirenti (alberto.pulvirenti@ct.infn.it)
//
#include "TMath.h"

#include "AliExternalTrackParam.h"

#include "AliRsnDaughter.h"
#include "AliRsnCutPID.h"

ClassImp(AliRsnCutPID)

//_________________________________________________________________________________________________
AliRsnCutPID::AliRsnCutPID() :
   AliRsnCut(),
   fPerfect(kFALSE),
   fUseDefault(kTRUE)
{
//
// Default constructor.
// Sets the cut to realistic PID with default weights,
// and defines the 'fMinI' value of the base class as the PID
// to which we want to compare this object.
//

   Int_t i;

   for (i = 0; i < kDetectors; i++) {
      fUseDetector[i] = kFALSE;
      fPtThreshold[i] = 0.0;
      fGoAboveThreshold[i] = kTRUE;
   }

   for (i = 0; i < AliPID::kSPECIES; i++) {
      fWeight[i] = 0.0;
      fPrior[i] = 1.0;
   }

   SetTargetType(AliRsnTarget::kDaughter);
}

//_________________________________________________________________________________________________
AliRsnCutPID::AliRsnCutPID(const char *name, AliPID::EParticleType pid, Double_t probMin, Bool_t perfectPID) :
   AliRsnCut(name, AliRsnCut::kDaughter, (Int_t)pid),
   fPerfect(perfectPID),
   fUseDefault(kTRUE)
{
//
// Default constructor.
// Sets the cut to realistic PID with default weights,
// and defines the 'fMinI' value of the base class as the PID
// to which we want to compare this object.
//

   Int_t i;

   for (i = 0; i < kDetectors; i++) {
      fUseDetector[i] = kFALSE;
      fPtThreshold[i] = 0.0;
      fGoAboveThreshold[i] = kTRUE;
   }

   for (i = 0; i < AliPID::kSPECIES; i++) {
      fWeight[i] = 0.0;
      fPrior[i] = 1.0;
   }

   fMinD = probMin;
   fMaxD = 1.000001;

   SetTargetType(AliRsnTarget::kDaughter);
}

//_____________________________________________________________________________
Bool_t AliRsnCutPID::CheckThreshold(EDetector det, Double_t value) const
{
//
// Checks if the passed value (track pT) stays in the
// interval where the detector should be accepted
//

   if (!CheckBounds(det)) return kFALSE;

   if (fGoAboveThreshold[det]) return (value >= fPtThreshold[det]);
   else return (value <= fPtThreshold[det]);
}

//_________________________________________________________________________________________________
Bool_t AliRsnCutPID::ComputeWeights(AliRsnDaughter *daughter)
{
//
// Compute the PID weights using the class settings.
// If the argument is an ESD track, customized weights can be computed
// It the argument is a track (ESD or AOD), at least default weights
// can be computed, otherwise, no weights can be defined.
// The return value tells if the operation was successful.
//

   Int_t  i, j;
   Bool_t useDefault = fUseDefault;
   Bool_t perfectPID = fPerfect;
   if (perfectPID && !daughter->GetRefMC()) return kFALSE;
   if (!daughter->Ref2ESDtrack()) useDefault = kTRUE;
   if (!daughter->Ref2ESDtrack() && !daughter->Ref2AODtrack()) return kFALSE;

   // if perfect PID ise required,
   // compare the PDG code of the type stored in 'fMinI' of the cut
   // and that of the particle which is checked, looking at its MC
   if (perfectPID) {
      i = TMath::Abs(AliPID::ParticleCode(fMinI));
      j = daughter->GetPDG();
      return (i == j);
   }

   // if default weights are (or need to be) used,
   // they are taken directly and function exits
   if (useDefault) {
      if (daughter->Ref2ESDtrack())
         daughter->Ref2ESDtrack()->GetESDpid(fWeight);
      else {
         for (i = 0; i < AliPID::kSPECIES; i++)
            fWeight[i] = daughter->Ref2AODtrack()->PID()[i];
      }
      return kTRUE;
   }

   // if we arrive here, this means that we have an ESD track
   // and we want to customize the PID
   AliESDtrack *track = daughter->Ref2ESDtrack();
   Double_t     w[kDetectors][AliPID::kSPECIES];
   track->GetITSpid(w[kITS]);
   track->GetTPCpid(w[kTPC]);
   track->GetTRDpid(w[kTRD]);
   track->GetTOFpid(w[kTOF]);
   track->GetHMPIDpid(w[kHMPID]);

   // if a detector is excluded or the track has not the right pT
   // all related weights are set to 1 in order not to contribute
   // to final multiplication
   for (i = 0; i < kDetectors; i++) {
      if (!fUseDetector[i] || !CheckThreshold((EDetector)i, track->Pt())) {
         for (j = 0; j < AliPID::kSPECIES; j++) {
            w[i][j] = 1.0;
         }
      }
   }

   // multiplicate all weights to compute final one
   for (i = 0; i < AliPID::kSPECIES; i++) {
      fWeight[i] = w[kITS][i] * w[kTPC][i] * w[kTRD][i] * w[kTOF][i] * w[kHMPID][i];
   }

   return kTRUE;
}

//_________________________________________________________________________________________________
Int_t AliRsnCutPID::RealisticPID(AliRsnDaughter * const daughter, Double_t &prob)
{
//
// Combines the weights collected from chosen detectors with the priors
// and gives the corresponding particle with the largest probability,
// in terms of the AliPID particle type enumeration.
// The argument, passed by reference, gives the corresponding probability,
// since the cut could require that it is larger than a threshold.
//

   // try to compute the weights
   if (!ComputeWeights(daughter)) {
      prob = -1.0;
      return AliPID::kUnknown;
   }

   // combine with priors and normalize
   Int_t    i;
   Double_t sum = 0.0, w[AliPID::kSPECIES];
   for (i = 0; i < AliPID::kSPECIES; i++) {
      w[i] = fWeight[i] * fPrior[i];
      sum += w[i];
   }
   if (sum <= 0.0) {
      AliError("Sum = 0");
      prob = -1.0;
      return AliPID::kUnknown;
   }
   for (i = 0; i < AliPID::kSPECIES; i++) w[i] /= sum;

   // find the largest probability and related PID
   Int_t ibest = 0;
   prob = w[0];
   for (i = 1; i < AliPID::kSPECIES; i++) {
      if (w[i] > prob) {
         prob = w[i];
         ibest = i;
      }
   }

   // return the value, while the probability
   // will be consequentially set
   return ibest;
}

//_________________________________________________________________________________________________
Int_t AliRsnCutPID::PerfectPID(AliRsnDaughter * const daughter)
{
//
// If MC is present, retrieve the particle corresponding to the used track
// (using the fRefMC data member) and then return the true particle type
// taken from the PDG code of the reference particle itself, converted
// into the enumeration from AliPID object.
//

   // works only if the MC is present
   if (!daughter->GetRefMC()) return AliPID::kUnknown;

   // get the PDG code of the particle
   Int_t pdg = daughter->GetPDG();

   // loop over all species listed in AliPID to find the match
   Int_t i;
   for (i = 0; i < AliPID::kSPECIES; i++) {
      if (AliPID::ParticleCode(i) == TMath::Abs(pdg)) return i;
   }

   return AliPID::kUnknown;
}

//_________________________________________________________________________________________________
Bool_t AliRsnCutPID::IsSelected(TObject *object)
{
//
// Cut checker.
//

   // convert the object into the unique correct type

   if (!TargetOK(object)) {
      AliError(Form("[%s]: this cut works only with AliRsnDaughter objects", GetName()));
      return kTRUE;
   }

   // if target is OK, do a dynamic cast
   AliRsnDaughter *daughter = fDaughter;

   // depending on the PID type, recalls the appropriate method:
   // in case of perfect PID, checks only if the PID type is
   // corresponding to the request in the cut (fMinI)
   // while in case of realistic PID checks also the probability
   // to be within the required interval
   if (fPerfect && daughter) {
      fCutValueI = PerfectPID(daughter);
      return OkValueI();
   } else if (daughter) {
      fCutValueI = RealisticPID(daughter, fCutValueD);
      return OkValueI() && OkRangeD();
   } else
      return kFALSE;
}

//__________________________________________________________________________________________________
void AliRsnCutPID::IncludeDetector(EDetector det, Double_t threshold, Bool_t goAbove)
{
//
// Include a detector for a customized weight computing
// and specify also its eventual threshold and if the detector
// must be used above or below the threshold.
// By default the threshold is zero and detector is always used.
//

   if (!CheckBounds(det)) return;

   fUseDetector[det] = kTRUE;
   fPtThreshold[det] = threshold;
   fGoAboveThreshold[det] = goAbove;
}
