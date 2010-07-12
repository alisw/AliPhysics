//
// Class AliRsnCutBetheBloch
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
#include "AliRsnCutBetheBloch.h"

ClassImp(AliRsnCutBetheBloch)

//_________________________________________________________________________________________________
AliRsnCutBetheBloch::AliRsnCutBetheBloch() :
    AliRsnCut(),
    fCorrect(kTRUE),
    fMIP(50.0),
    fType(AliPID::kUnknown)
{
//
// Default constructor.
//

  fConst[0] = fConst[1] = fConst[2] = fConst[3] = fConst[4] = 0.0;
}

//_________________________________________________________________________________________________
AliRsnCutBetheBloch::AliRsnCutBetheBloch
(const char *name, Double_t fractionRange, AliPID::EParticleType type, Double_t mip, Bool_t correct) :
    AliRsnCut(name, 0.0, fractionRange),
    fCorrect(correct),
    fMIP(mip),
    fType(type)
{
//
// Main constructor.
// the cut range is the relative fraction of the value:
// BB*(1-fraction) < TPC < BB*(1+fraction)
// which means:
// -fraction < (TPC - BB)/BB < fraction
//

  fConst[0] = fConst[1] = fConst[2] = fConst[3] = fConst[4] = 0.0;
}

//_____________________________________________________________________________
Double_t AliRsnCutBetheBloch::BetheBloch(AliRsnDaughter * const trackRef)
{
//
// Computes the theoretical dE/dx according to
// a given mass hypothesis, from which betaGamma is computed
//
// This is the empirical ALEPH parameterization of the Bethe-Bloch formula.
// It is normalized to 1 at the minimum.
//
// The default values for the kp* parameters are for ALICE TPC.
// The value is computed in MIP units, multiplied by 50 to have it in energy.
//

  AliPID pid;
  Double_t mass = pid.ParticleMass(fType);

  // get the track momentum at the inner wall of TPC: if absent cut is not passed
  AliExternalTrackParam track(*trackRef->GetRefESD()->GetInnerParam());

  Double_t betaGamma = track.P() / mass;
  Double_t beta = betaGamma / TMath::Sqrt(1.0 + betaGamma * betaGamma);
  Double_t aa = TMath::Power(beta, fConst[3]);
  Double_t bb = TMath::Power(1.0 / betaGamma, fConst[4]);

  bb = TMath::Log(fConst[2] + bb);

  Double_t out = (fConst[1] - aa - bb) * fConst[0] / aa;

  if (fCorrect) {
    Double_t kMeanCorr = 0.1;
    Double_t meanCorr = (1 + (out - 1) * kMeanCorr);
    out *= meanCorr;
  }

  return out * fMIP;
}

//_____________________________________________________________________________
Double_t AliRsnCutBetheBloch::RelDiff(AliRsnDaughter *track)
{
//
// Relative difference between BB value and TPC signal
//

  if (!track->GetRefESD()) return -99999.9;

  // compute Bethe-Bloch with the given mass hypothesis
  Double_t bb = BetheBloch(track);
  return TMath::Abs((track->GetRefESD()->GetTPCsignal() - bb) / bb);
}

//_________________________________________________________________________________________________
Bool_t AliRsnCutBetheBloch::IsSelected(ETarget tgt, AliRsnDaughter *track)
{
//
// Cut checker.
//

  // coherence check
  if (tgt != AliRsnCut::kParticle) {
    AliError(Form("[%s] Wrong target. Skipping cut", GetName()));
    return kTRUE;
  }

  // if the required PID of the track is not the same as the
  // reference of the cut, the cut is automatically skipped
  if (track->RequiredPID() != fType) return kTRUE;

  // retrieve the TPC signal
  AliESDtrack *esd = track->GetRefESD();
  if (!esd) {
    AliError("ESD information unavailable");
    return kTRUE;
  }
  if (!track->GetRefESD()->GetInnerParam()) {
    AliDebug(AliLog::kDebug+2, "Rejecting a track with no info at the TPC inner wall");
    return kFALSE;
  }

  // the cut range is the relative fraction of the value:
  // BB*(1-fraction) < TPC < BB*(1+fraction)
  // which means:
  // -fraction < (TPC - BB)/BB < fraction
  // so we must compute the cut value accordingly
  fCutValueD = RelDiff(track);

  // then, this cut is checked inside the range
  return OkRange();
}

//_________________________________________________________________________________________________
Bool_t AliRsnCutBetheBloch::IsSelected(ETarget /*tgt*/, AliRsnPairParticle* /*pair*/)
{
//
// Cut checker
//

  AliWarning("Cannot apply this cut to pairs");
  return kTRUE;
}

//_________________________________________________________________________________________________
Bool_t AliRsnCutBetheBloch::IsSelected(ETarget /*tgt*/, AliRsnEvent* /*event*/)
{
//
// Cut checker
//

  AliWarning("Cannot apply this cut to events");
  return kTRUE;
}

//_________________________________________________________________________________________________
Bool_t AliRsnCutBetheBloch::IsSelected(ETarget /*tgt*/, AliRsnEvent* /*ev1*/, AliRsnEvent* /*ev2*/)
{
//
// Cut checker
//

  AliWarning("Cannot apply this cut to event mixing");
  return kTRUE;
}
