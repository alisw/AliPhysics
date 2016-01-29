///
/// \file AliFemtoV0TrackCutNSigmaFilter.cxx
///


#include "AliFemtoV0TrackCutNSigmaFilter.h"
#include "AliESDtrack.h"

#include <TH1D.h>

#include <cstdio>

#ifdef __ROOT__
  /// \cond CLASSIMP
  ClassImp(AliFemtoV0TrackCutNSigmaFilter);
  /// \endcond
#endif


AliFemtoV0TrackCutNSigmaFilter::AliFemtoV0TrackCutNSigmaFilter():
  AliFemtoV0TrackCut()
  , fUseCustomPionNSigmaFilter(false)
  , fUseCustomKaonNSigmaFilter(false)
  , fUseCustomProtonNSigmaFilter(false)

  , fPionNSigmaFilter(NULL)
  , fKaonNSigmaFilter(NULL)
  , fProtonNSigmaFilter(NULL)

{
  /* no-op */
}

AliFemtoV0TrackCutNSigmaFilter::~AliFemtoV0TrackCutNSigmaFilter()
{
  delete fPionNSigmaFilter;
  delete fKaonNSigmaFilter;
  delete fProtonNSigmaFilter;
  delete fMinvPurityAidHistoV0;
}

bool AliFemtoV0TrackCutNSigmaFilter::Pass(const AliFemtoV0* aV0)
{
  // test the particle and return
  // true if it meets all the criteria
  // false if it doesn't meet at least one of the criteria

  const Float_t pt = aV0->PtV0(),
               eta = aV0->EtaV0();

  //kinematic cuts
  if (TMath::Abs(eta) > fEta) return false;  //put in kinematic cuts by hand
  if (pt < fPtMin || fPtMax < pt) return false;
  if (TMath::Abs(aV0->EtaPos()) > fMaxEtaDaughters) return false;
  if (TMath::Abs(aV0->EtaNeg()) > fMaxEtaDaughters) return false;

  const Float_t pt_pos = aV0->PtPos(),
                pt_neg = aV0->PtNeg();

  if (pt_pos < fPtMinPosDaughter || fPtMaxPosDaughter < pt_pos) return false;
  if (pt_neg < fPtMaxNegDaughter || fPtMaxNegDaughter < pt_neg) return false;

  const Float_t mass_lambda = aV0->MassLambda();

  //V0 from kinematics information
  switch (fParticleType) {
  case kLambdaMC:
    if (!(mass_lambda > fInvMassLambdaMin && mass_lambda < fInvMassLambdaMax) || !(aV0->PosNSigmaTPCP() == 0)) {
      return false;
    } else {
      return true;
    }
  case kAntiLambdaMC:
    if (!(mass_lambda > fInvMassLambdaMin && mass_lambda < fInvMassLambdaMax) || !(aV0->NegNSigmaTPCP() == 0)) {
      return false;
    } else {
      return true;
    }
  case kAll:
    if (!(aV0->MassK0Short() > fInvMassK0sMin && aV0->MassK0Short() < fInvMassK0sMax) || !(aV0->NegNSigmaTPCP() == 0)) {
      return false;
    } else {
      return true;
    }
  default:
    break;
  }

  //quality cuts
  if (aV0->OnFlyStatusV0() != fOnFlyStatus) return false;
  if (aV0->StatusNeg() == 999 || aV0->StatusPos() == 999) return false;
  if (aV0->TPCNclsPos() < fTPCNclsDaughters) return false;
  if (aV0->TPCNclsNeg() < fTPCNclsDaughters) return false;
  if (aV0->NdofPos() > fNdofDaughters) return false;
  if (aV0->NdofNeg() > fNdofDaughters) return false;
  if (!(aV0->StatusNeg() & fStatusDaughters)) return false;
  if (!(aV0->StatusPos() & fStatusDaughters)) return false;


  //DCA between daughter particles
  if (TMath::Abs(aV0->DcaV0Daughters()) > fMaxDcaV0Daughters)
    return false;

  //DCA of daughters to primary vertex
  if (TMath::Abs(aV0->DcaPosToPrimVertex()) < fMinDcaDaughterPosToVert || TMath::Abs(aV0->DcaNegToPrimVertex()) < fMinDcaDaughterNegToVert)
    return false;

  //DCA V0 to prim vertex
  if (TMath::Abs(aV0->DcaV0ToPrimVertex()) > fMaxDcaV0 || TMath::Abs(aV0->DcaV0ToPrimVertex()) < fMinDcaV0)
    return false;

  //cos pointing angle
  if (aV0->CosPointingAngle() < fMinCosPointingAngle)
    return false;

  //decay length
  if (aV0->DecayLengthV0() > fMaxDecayLength)
    return false;

  bool pid_check = false;
  // Looking for lambdas = proton + pim
  if (fParticleType == kLambda) {
    if (IsProtonNSigma(aV0->PtPos(), aV0->PosNSigmaTPCP(), aV0->PosNSigmaTOFP())) //proton
      if (IsPionNSigma(aV0->PtNeg(), aV0->NegNSigmaTPCPi(), aV0->NegNSigmaTOFPi())) { //pion-
        pid_check = true;
        if(fBuildPurityAidV0) fMinvPurityAidHistoV0->Fill(aV0->MassLambda());
        //invariant mass lambda
        if (aV0->MassLambda() < fInvMassLambdaMin || aV0->MassLambda() > fInvMassLambdaMax)
          return false;
      }

  }//Looking for antilambdas =  antiproton + pip
  else if (fParticleType == kAntiLambda) {
    if (IsProtonNSigma(aV0->PtNeg(), aV0->NegNSigmaTPCP(), aV0->NegNSigmaTOFP())) //antiproton
      if (IsPionNSigma(aV0->PtPos(), aV0->PosNSigmaTPCPi(), aV0->PosNSigmaTOFPi())) { //pion+
        pid_check = true;
        if(fBuildPurityAidV0) fMinvPurityAidHistoV0->Fill(aV0->MassAntiLambda());
        //invariant mass antilambda
        if (aV0->MassAntiLambda() < fInvMassLambdaMin || aV0->MassAntiLambda() > fInvMassLambdaMax)
          return false;
      }
  }//Looking for K0s = pip + pim
  else if (fParticleType == kK0s) {
    if (IsPionNSigma(aV0->PtNeg(), aV0->NegNSigmaTPCPi(), aV0->NegNSigmaTOFPi())) //pion-
      if (IsPionNSigma(aV0->PtPos(), aV0->PosNSigmaTPCPi(), aV0->PosNSigmaTOFPi())) { //pion+
        pid_check = true;
        if(fBuildPurityAidV0) fMinvPurityAidHistoV0->Fill(aV0->MassK0Short());
        //invariant mass K0s
        if (aV0->MassK0Short() < fInvMassK0sMin || aV0->MassK0Short() > fInvMassK0sMax)
          return false;
      }
  }

  if (!pid_check) {
    return false;
  }

  return true;
}


AliFemtoString AliFemtoV0TrackCutNSigmaFilter::Report()
{
  // Prepare report from the execution
  TString report;
  report += TString::Format("Usings custom Pion NSigma Filter:\t%i\n", fUseCustomPionNSigmaFilter)
          + TString::Format("Usings custom Kaon NSigma Filter:\t%i\n", fUseCustomKaonNSigmaFilter)
          + TString::Format("Usings custom Proton NSigma Filter:\t%i\n", fUseCustomProtonNSigmaFilter)
          + AliFemtoV0TrackCut::Report();
  return AliFemtoString(report);
}


TList* AliFemtoV0TrackCutNSigmaFilter::AppendSettings(TList *settings, const TString &prefix) const
{
  settings->AddVector(
    new TObjString(prefix + TString::Format("AliFemtoV0TrackCutNSigmaFilter.InvMassLambdaMin=%lf", fInvMassLambdaMin)),

  NULL);

  return settings;
}

bool AliFemtoV0TrackCutNSigmaFilter::IsKaonNSigma(float mom, float nsigmaTPCK, float nsigmaTOFK)
{
  if (fUseCustomKaonNSigmaFilter) {
    return fKaonNSigmaFilter->Pass(mom, nsigmaTPCK, nsigmaTOFK);
  } else {
    return AliFemtoV0TrackCut::IsKaonNSigma(mom, nsigmaTPCK, nsigmaTOFK);
  }
}


bool AliFemtoV0TrackCutNSigmaFilter::IsPionNSigma(float mom, float nsigmaTPCPi, float nsigmaTOFPi)
{
  if (fUseCustomPionNSigmaFilter) {
    return fPionNSigmaFilter->Pass(mom, nsigmaTPCPi, nsigmaTOFPi);
  } else {
    return AliFemtoV0TrackCut::IsPionNSigma(mom, nsigmaTPCPi, nsigmaTOFPi);
  }
}

bool AliFemtoV0TrackCutNSigmaFilter::IsProtonNSigma(float mom, float nsigmaTPCP, float nsigmaTOFP)
{
  if (fUseCustomProtonNSigmaFilter) {
    return fProtonNSigmaFilter->Pass(mom, nsigmaTPCP, nsigmaTOFP);
  } else {
    return AliFemtoV0TrackCut::IsProtonNSigma(mom,nsigmaTPCP,nsigmaTOFP);
  }
}


void AliFemtoV0TrackCutNSigmaFilter::CreateCustomNSigmaFilter(DaughterParticleType aDaughterType)
{
  switch (aDaughterType) {
  case kPion:
    fUseCustomPionNSigmaFilter = true;
    fPionNSigmaFilter = new AliFemtoNSigmaFilter();
    break;

  case kKaon:
    fUseCustomKaonNSigmaFilter = true;
    fKaonNSigmaFilter = new AliFemtoNSigmaFilter();
    break;

  case kProton:
    fUseCustomProtonNSigmaFilter = true;
    fProtonNSigmaFilter = new AliFemtoNSigmaFilter();
    break;

  default:
    cerr << "E-AliFemtoV0TrackCutNSigmaFilter::CreateCustomNSigmaFilter: Invalid DaughterParticleType"
            "selection '" << aDaughterType << "'.  No custom filter will be initialized!!!!!" << endl;
  }
}

void AliFemtoV0TrackCutNSigmaFilter::AddTPCAndTOFNSigmaCut(DaughterParticleType aDaughterType, double aMomMin, double aMomMax, double aNSigmaValueTPC, double aNSigmaValueTOF)
{
  if(aDaughterType == kPion) {fPionNSigmaFilter->AddTPCAndTOFCut(aMomMin,aMomMax,aNSigmaValueTPC,aNSigmaValueTOF);}
  else if(aDaughterType == kKaon) {fKaonNSigmaFilter->AddTPCAndTOFCut(aMomMin,aMomMax,aNSigmaValueTPC,aNSigmaValueTOF);}
  else if(aDaughterType == kProton) {fProtonNSigmaFilter->AddTPCAndTOFCut(aMomMin,aMomMax,aNSigmaValueTPC,aNSigmaValueTOF);}
  else {/*cerr << "Invalid DaughterParticleType selection in AddTPCAndTOFNSigmaCut.  No cut will be added to the filter!!!!!" << endl;*/}
}


void AliFemtoV0TrackCutNSigmaFilter::AddTPCNSigmaCut(DaughterParticleType aDaughterType, double aMomMin, double aMomMax, double aNSigmaValueTPC)
{
  if(aDaughterType == kPion) {fPionNSigmaFilter->AddTPCCut(aMomMin,aMomMax,aNSigmaValueTPC);}
  else if(aDaughterType == kKaon) {fKaonNSigmaFilter->AddTPCCut(aMomMin,aMomMax,aNSigmaValueTPC);}
  else if(aDaughterType == kProton) {fProtonNSigmaFilter->AddTPCCut(aMomMin,aMomMax,aNSigmaValueTPC);}
  else {/*cerr << "Invalid DaughterParticleType selection in AddTPCNSigmaCut.  No cut will be added to the filter!!!!!" << endl;*/}
}


void AliFemtoV0TrackCutNSigmaFilter::AddTOFNSigmaCut(DaughterParticleType aDaughterType, double aMomMin, double aMomMax, double aNSigmaValueTOF)
{
  if(aDaughterType == kPion) {fPionNSigmaFilter->AddTOFCut(aMomMin,aMomMax,aNSigmaValueTOF);}
  else if(aDaughterType == kKaon) {fKaonNSigmaFilter->AddTOFCut(aMomMin,aMomMax,aNSigmaValueTOF);}
  else if(aDaughterType == kProton) {fProtonNSigmaFilter->AddTOFCut(aMomMin,aMomMax,aNSigmaValueTOF);}
  else {/*cerr << "Invalid DaughterParticleType selection in AddTOFNSigmaCut.  No cut will be added to the filter!!!!!" << endl;*/}
}
