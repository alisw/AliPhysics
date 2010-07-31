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
  AliRsnCut(AliRsnCut::kDaughter),
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
  
  for (i = 0; i < kDetectors; i++) 
  {
    fUseDetector[i] = kFALSE;
    fPtThreshold[i] = 0.0;
    fGoAboveThreshold[i] = kTRUE;
  }
  
  for (i = 0; i < AliPID::kSPECIES; i++)
  {
    fWeight[i] = 0.0;
    fPrior[i] = 1.0;
  }
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
  
  for (i = 0; i < kDetectors; i++) 
  {
    fUseDetector[i] = kFALSE;
    fPtThreshold[i] = 0.0;
    fGoAboveThreshold[i] = kTRUE;
  }
  
  for (i = 0; i < AliPID::kSPECIES; i++)
  {
    fWeight[i] = 0.0;
    fPrior[i] = 1.0;
  }
  
  fMinD = probMin;
  fMaxD = 1.000001;
}

//_____________________________________________________________________________
Bool_t AliRsnCutPID::CheckThreshold(EDetector det, Double_t value)
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
  if (perfectPID && !daughter->GetRefMC()) perfectPID = kFALSE;
  if (!daughter->GetRefESDtrack()) useDefault = kTRUE;
  if (!daughter->GetRefESDtrack() && !daughter->GetRefAODtrack()) return kFALSE;
  
  // if perfect PID ise required, this overcomes all
  // in this case the weight of the correct species is set to 1
  // and the others to 0
  // of course this happens only if there is a reference MC
  if (perfectPID)
  {
    j = TMath::Abs(daughter->GetRefMC()->Particle()->GetPdgCode());
    for (i = 0; i < AliPID::kSPECIES; i++)
    {
      if (AliPID::ParticleCode((AliPID::EParticleType)i) == j) fWeight[i] = 1.0;
      else fWeight[i] = 0.0;
    }
    return kTRUE;
  }
  
  // if default weights are (or need to be) used,
  // they are taken directly and function exits
  if (useDefault)
  {
    if (daughter->GetRefESDtrack())
      daughter->GetRefESDtrack()->GetESDpid(fWeight);
    else
    {
      for (i = 0; i < AliPID::kSPECIES; i++)
        fWeight[i] = daughter->GetRefAODtrack()->PID()[i];
    }
    return kTRUE;
  }
  
  // if we arrive here, this means that we have an ESD track
  // and we want to customize the PID
  AliESDtrack *track = daughter->GetRefESDtrack();
  Double_t     w[kDetectors][AliPID::kSPECIES];
  track->GetITSpid(w[kITS]);
  track->GetTPCpid(w[kTPC]);
  track->GetTRDpid(w[kTRD]);
  track->GetTOFpid(w[kTOF]);
  track->GetHMPIDpid(w[kHMPID]);

  // if a detector is excluded or the track has not the right pT
  // all related weights are set to 1 in order not to contribute
  // to final multiplication
  for (i = 0; i < kDetectors; i++) 
  {
    if (!fUseDetector[i] || !CheckThreshold((EDetector)i, track->Pt())) 
    {
      for (j = 0; j < AliPID::kSPECIES; j++) {
        w[i][j] = 1.0;
      }
    }
  }

  // multiplicate all weights to compute final one
  for (i = 0; i < AliPID::kSPECIES; i++) 
  {
    fWeight[i] = w[kITS][i] * w[kTPC][i] * w[kTRD][i] * w[kTOF][i] * w[kHMPID][i];
  }
  
  return kTRUE;
}

//_________________________________________________________________________________________________
Bool_t AliRsnCutPID::IsSelected(TObject *obj1, TObject* /*obj2*/)
{
//
// Cut checker.
//

  // coherence check
  if (!AliRsnCut::TargetOK(obj1))
  {
    AliError(Form("Wrong target. Skipping cut", GetName()));
    return kTRUE;
  }
  
  // convert the object into the unique correct type
  AliRsnDaughter *daughter = dynamic_cast<AliRsnDaughter*>(obj1);
  
  // try to compute the weights
  if (!ComputeWeights(daughter)) return kFALSE;
  
  // combine with priors and get the majority
  Int_t    i;
  Double_t sum = 0.0, w[AliPID::kSPECIES];
  for (i = 0; i < AliPID::kSPECIES; i++)
  {
    w[i] = fWeight[i] * fPrior[i];
    sum += w[i];
  }
  if (sum <= 0.0)
  {
    AliError("Sum = 0");
    return kFALSE;
  }
  for (i = 0; i < AliPID::kSPECIES; i++) w[i] /= sum;
  
  // find the largest probability and related PID
  // and assign them to the mother class members which
  // are checked for the cut
  fCutValueI = 0;
  fCutValueD = w[0];
  for (i = 1; i < AliPID::kSPECIES; i++)
  {
    if (w[i] > fCutValueD) 
    {
      fCutValueD = w[i];
      fCutValueI = i;
    }
  }
  
  // if the best probability is too small, the cut is failed anyway
  if (!OkRangeD()) return kFALSE;
  
  // if the best probability is OK, the cut is passed
  // if it correspond to the right particle
  return OkValue();
}

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
