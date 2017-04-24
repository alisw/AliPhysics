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

  , fUseCustomK0sRejectionFilters(false)
  , fUseCustomLambdaRejectionFilters(false)
  , fUseCustomAntiLambdaRejectionFilters(false)

  , fK0sRejectionFilters(0)
  , fLambdaRejectionFilters(0)
  , fAntiLambdaRejectionFilters(0)

  , fBuildCTau(false)
  , fCTau(0)


{
  /* no-op */
}

AliFemtoV0TrackCutNSigmaFilter::AliFemtoV0TrackCutNSigmaFilter(const AliFemtoV0TrackCutNSigmaFilter& aCut) :
  AliFemtoV0TrackCut(aCut),

  fUseCustomPionNSigmaFilter(aCut.fUseCustomPionNSigmaFilter),
  fUseCustomKaonNSigmaFilter(aCut.fUseCustomKaonNSigmaFilter),
  fUseCustomProtonNSigmaFilter(aCut.fUseCustomProtonNSigmaFilter),

  fUseCustomK0sRejectionFilters(aCut.fUseCustomK0sRejectionFilters),
  fUseCustomLambdaRejectionFilters(aCut.fUseCustomLambdaRejectionFilters),
  fUseCustomAntiLambdaRejectionFilters(aCut.fUseCustomAntiLambdaRejectionFilters),

  fK0sRejectionFilters(aCut.fK0sRejectionFilters),
  fLambdaRejectionFilters(aCut.fLambdaRejectionFilters),
  fAntiLambdaRejectionFilters(aCut.fAntiLambdaRejectionFilters),
  fBuildCTau(aCut.fBuildCTau)
{
  if(aCut.fPionNSigmaFilter) fPionNSigmaFilter = new AliFemtoNSigmaFilter(*aCut.fPionNSigmaFilter);
  if(aCut.fKaonNSigmaFilter) fKaonNSigmaFilter = new AliFemtoNSigmaFilter(*aCut.fKaonNSigmaFilter);
  if(aCut.fProtonNSigmaFilter) fProtonNSigmaFilter = new AliFemtoNSigmaFilter(*aCut.fProtonNSigmaFilter);

  if(aCut.fCTau) fCTau = new TH1D(*aCut.fCTau);
  else fCTau = 0;
}


AliFemtoV0TrackCutNSigmaFilter& AliFemtoV0TrackCutNSigmaFilter::operator=(const AliFemtoV0TrackCutNSigmaFilter& aCut)
{
  if(this == &aCut) {return *this;}

  AliFemtoV0TrackCut::operator=(aCut);

  fUseCustomPionNSigmaFilter = aCut.fUseCustomPionNSigmaFilter;
  fUseCustomKaonNSigmaFilter = aCut.fUseCustomKaonNSigmaFilter;
  fUseCustomProtonNSigmaFilter = aCut.fUseCustomProtonNSigmaFilter;

  fUseCustomK0sRejectionFilters = aCut.fUseCustomK0sRejectionFilters;
  fUseCustomLambdaRejectionFilters = aCut.fUseCustomLambdaRejectionFilters;
  fUseCustomAntiLambdaRejectionFilters = aCut.fUseCustomAntiLambdaRejectionFilters;

  fK0sRejectionFilters = aCut.fK0sRejectionFilters;
  fLambdaRejectionFilters = aCut.fLambdaRejectionFilters;
  fAntiLambdaRejectionFilters = aCut.fAntiLambdaRejectionFilters;

  fBuildCTau = aCut.fBuildCTau;

  if(aCut.fPionNSigmaFilter) fPionNSigmaFilter = new AliFemtoNSigmaFilter(*aCut.fPionNSigmaFilter);
  if(aCut.fKaonNSigmaFilter) fKaonNSigmaFilter = new AliFemtoNSigmaFilter(*aCut.fKaonNSigmaFilter);
  if(aCut.fProtonNSigmaFilter) fProtonNSigmaFilter = new AliFemtoNSigmaFilter(*aCut.fProtonNSigmaFilter);

  if(aCut.fCTau) fCTau = new TH1D(*aCut.fCTau);
  else fCTau = 0;

  return *this;
}


AliFemtoV0TrackCutNSigmaFilter* AliFemtoV0TrackCutNSigmaFilter::Clone()
{
  return(new AliFemtoV0TrackCutNSigmaFilter(*this));
}

AliFemtoV0TrackCutNSigmaFilter::~AliFemtoV0TrackCutNSigmaFilter()
{
  delete fPionNSigmaFilter;
  delete fKaonNSigmaFilter;
  delete fProtonNSigmaFilter;
  delete fMinvPurityAidHistoV0;

//TODO            TODO           TODO               TODO 
  if(fUseCustomK0sRejectionFilters)
  {
    fK0sRejectionFilters.clear();
    fK0sRejectionFilters.shrink_to_fit();
  } 

  if(fUseCustomLambdaRejectionFilters)
  {
    fLambdaRejectionFilters.clear();
    fLambdaRejectionFilters.shrink_to_fit();
  } 

  if(fUseCustomAntiLambdaRejectionFilters)
  {
    fAntiLambdaRejectionFilters.clear();
    fAntiLambdaRejectionFilters.shrink_to_fit();
  } 
}

bool AliFemtoV0TrackCutNSigmaFilter::Pass(const AliFemtoV0* aV0)
{
  // test the particle and return
  // true if it meets all the criteria
  // false if it doesn't meet at least one of the criteria

  Float_t pt = aV0->PtV0();
  Float_t eta = aV0->EtaV0();

  //kinematic cuts
  if (TMath::Abs(eta) > fEta) return false;  //put in kinematic cuts by hand
  if (pt < fPtMin || fPtMax < pt) return false;
  if (TMath::Abs(aV0->EtaPos()) > fMaxEtaDaughters) return false;
  if (TMath::Abs(aV0->EtaNeg()) > fMaxEtaDaughters) return false;

  if (aV0->PtPos() < fPtMinPosDaughter) return false;
  if (aV0->PtNeg() < fPtMinNegDaughter) return false;
  if (aV0->PtPos() > fPtMaxPosDaughter) return false;
  if (aV0->PtNeg() > fPtMaxNegDaughter) return false;

  //V0 from kinematics information
  if (fParticleType == kLambdaMC ) {
    if (!(aV0->MassLambda() > fInvMassLambdaMin && aV0->MassLambda() < fInvMassLambdaMax) || !(aV0->PosNSigmaTPCP() == 0)) {
      return false;
    } else {
      return true;
    }
  } else if (fParticleType == kAntiLambdaMC) {
    if (!(aV0->MassLambda() > fInvMassLambdaMin && aV0->MassLambda() < fInvMassLambdaMax) || !(aV0->NegNSigmaTPCP() == 0)) {
      return false;
    } else {
      return true;
    }
  } else if (fParticleType == kK0sMC) {
    if (!(aV0->MassK0Short() > fInvMassK0sMin && aV0->MassK0Short() < fInvMassK0sMax) || !(aV0->NegNSigmaTPCP() == 0))
      return false;
    else {
      return true;
    }
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

  //fiducial volume radius
  if(aV0->RadiusV0()<fRadiusV0Min || aV0->RadiusV0()>fRadiusV0Max)
    return false;

  //DCA between daughter particles
  if (TMath::Abs(aV0->DcaV0Daughters()) > fMaxDcaV0Daughters)
    return false;

  //DCA of daughters to primary vertex
  if (TMath::Abs(aV0->DcaPosToPrimVertex()) < fMinDcaDaughterPosToVert || TMath::Abs(aV0->DcaNegToPrimVertex()) < fMinDcaDaughterNegToVert)
    return false;

  //DCA V0 to prim vertex
  if (TMath::Abs(aV0->DcaV0ToPrimVertex()) > fMaxDcaV0 || TMath::Abs(aV0->DcaV0ToPrimVertex()) < fMinDcaV0)
    return false;

  //becomes obsolete - wrong name of the data memeber and the corresnponding methods (by default is set to fMaxCosPointingAngle = 0.0)
  //cos pointing angle
  if (aV0->CosPointingAngle() < fMaxCosPointingAngle)
    return false;

  //this is the correct name of the data member and the corresponding methods (we are accepting cos(pointing angle bigger than certain minimum)
  //cos pointing angle
  if (aV0->CosPointingAngle() < fMinCosPointingAngle)
    return false;

  //decay length
  if (aV0->DecayLengthV0() > fMaxDecayLength)
    return false;

  //transverse decay vertex
  //TO BE IMPLEMENTED (maybe in future)

  if (fParticleType == kAll)
    return true;


  bool pid_check = false;
  // Looking for lambdas = proton + pim
  if(fParticleType == kLambda) 
  {
    if(IsProtonNSigma(aV0->PtPos(), aV0->PosNSigmaTPCP(), aV0->PosNSigmaTOFP(),fNsigmaPosDaughterTPC,fNsigmaPosDaughterTOF,fRequireTOFProton)) //proton
    {
      if (IsPionNSigma(aV0->PtNeg(), aV0->NegNSigmaTPCPi(), aV0->NegNSigmaTOFPi(),fNsigmaNegDaughterTPC,fNsigmaNegDaughterTOF,fRequireTOFPion)) //pion-
      {
        pid_check = true;
        if(fUseLooseInvMassCut)
        {
          if(aV0->MassLambda()<fLooseInvMassMin || aV0->MassLambda()>fLooseInvMassMax) return false;
        }
        if(fRemoveMisidentified)
        {
          if(IsMisIDK0s(aV0))
          {
            if(fBuildMisIDHistograms)
            {
              fLambdaMassOfMisIDV0->Fill(aV0->MassLambda());
              fK0sMassOfMisIDV0->Fill(aV0->MassK0Short());
            }
            return false;
          }
        }
        if(fBuildPurityAidV0) fMinvPurityAidHistoV0->Fill(aV0->MassLambda());
        //invariant mass lambda
        if(aV0->MassLambda() < fInvMassLambdaMin || aV0->MassLambda() > fInvMassLambdaMax) return false;
      }
    }
  }

  //Looking for antilambdas =  antiproton + pip
  else if(fParticleType == kAntiLambda)
  {
    if(IsProtonNSigma(aV0->PtNeg(), aV0->NegNSigmaTPCP(), aV0->NegNSigmaTOFP(),fNsigmaNegDaughterTPC,fNsigmaNegDaughterTOF,fRequireTOFProton)) //antiproton
    {
      if(IsPionNSigma(aV0->PtPos(), aV0->PosNSigmaTPCPi(), aV0->PosNSigmaTOFPi(),fNsigmaPosDaughterTPC,fNsigmaPosDaughterTOF,fRequireTOFPion)) //pion+
      {
        pid_check = true;
        if(fUseLooseInvMassCut)
        {
          if(aV0->MassAntiLambda()<fLooseInvMassMin || aV0->MassAntiLambda()>fLooseInvMassMax) return false;
        }
        if(fRemoveMisidentified)
        {
          if(IsMisIDK0s(aV0))
          {
            if(fBuildMisIDHistograms)
            {
              fAntiLambdaMassOfMisIDV0->Fill(aV0->MassAntiLambda());
              fK0sMassOfMisIDV0->Fill(aV0->MassK0Short());
            }
            return false;
          }
        }
        if(fBuildPurityAidV0) fMinvPurityAidHistoV0->Fill(aV0->MassAntiLambda());
        //invariant mass antilambda
        if (aV0->MassAntiLambda() < fInvMassLambdaMin || aV0->MassAntiLambda() > fInvMassLambdaMax) return false;
      }
    }
  }

  //Looking for K0s = pip + pim
  else if (fParticleType == kK0s)
  {
    if(IsPionNSigma(aV0->PtNeg(), aV0->NegNSigmaTPCPi(), aV0->NegNSigmaTOFPi(),fNsigmaNegDaughterTPC,fNsigmaNegDaughterTOF,fRequireTOFPion)) //pion-
    {
      if(IsPionNSigma(aV0->PtPos(), aV0->PosNSigmaTPCPi(), aV0->PosNSigmaTOFPi(),fNsigmaPosDaughterTPC,fNsigmaPosDaughterTOF,fRequireTOFPion)) //pion+
      {
        pid_check = true;
        if(fUseLooseInvMassCut)
        {
          if(aV0->MassK0Short()<fLooseInvMassMin || aV0->MassK0Short()>fLooseInvMassMax) return false;
        }
        if(fRemoveMisidentified)
        {
          if(IsMisIDLambda(aV0))
          {
            if(fBuildMisIDHistograms)
            {
              fK0sMassOfMisIDV0->Fill(aV0->MassK0Short());
              fLambdaMassOfMisIDV0->Fill(aV0->MassLambda());
            }
            return false;
          }
          else if(IsMisIDAntiLambda(aV0))
          {
            if(fBuildMisIDHistograms)
            {
              fK0sMassOfMisIDV0->Fill(aV0->MassK0Short());
              fAntiLambdaMassOfMisIDV0->Fill(aV0->MassAntiLambda());
            }
            return false;
          }
        }
        if(fBuildPurityAidV0) fMinvPurityAidHistoV0->Fill(aV0->MassK0Short());
        //invariant mass K0s
        if (aV0->MassK0Short() < fInvMassK0sMin || aV0->MassK0Short() > fInvMassK0sMax) return false;
      }
    }
  }

  if (!pid_check) return false;
  if(fBuildCTau) fCTau->Fill(GetCTau(aV0));
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


bool AliFemtoV0TrackCutNSigmaFilter::IsPionNSigma(float mom, float nsigmaTPCPi, float nsigmaTOFPi, double nsigmacutTPC, double nsigmacutTOF, bool requireTOF)
{
  if (fUseCustomPionNSigmaFilter) {
    return fPionNSigmaFilter->Pass(mom, nsigmaTPCPi, nsigmaTOFPi);
  } else {
    return AliFemtoV0TrackCut::IsPionNSigma(mom,nsigmaTPCPi,nsigmaTOFPi,nsigmacutTPC,nsigmacutTOF,requireTOF);
  }
}

bool AliFemtoV0TrackCutNSigmaFilter::IsProtonNSigma(float mom, float nsigmaTPCP, float nsigmaTOFP, double nsigmacutTPC, double nsigmacutTOF, bool requireTOF)
{
  if (fUseCustomProtonNSigmaFilter) {
    return fProtonNSigmaFilter->Pass(mom, nsigmaTPCP, nsigmaTOFP);
  } else {
    return AliFemtoV0TrackCut::IsProtonNSigma(mom,nsigmaTPCP,nsigmaTOFP,nsigmacutTPC,nsigmacutTOF,requireTOF);
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
    break;
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


void AliFemtoV0TrackCutNSigmaFilter::CreateCustomV0Rejection(AliFemtoV0Type aV0TypeToReject)
{
  switch (aV0TypeToReject) {
  case kK0s:
    fUseCustomK0sRejectionFilters = true;
    fK0sRejectionFilters.emplace_back(); //add AliFemtoNSigmaFilter object for (positive) daughter 1
    fK0sRejectionFilters.emplace_back(); //add AliFemtoNSigmaFilter object for (negative) daughter 2
    break;

  case kLambda:
    fUseCustomLambdaRejectionFilters = true;
    fLambdaRejectionFilters.emplace_back(); //add AliFemtoNSigmaFilter object for (positive) daughter 1
    fLambdaRejectionFilters.emplace_back(); //add AliFemtoNSigmaFilter object for (negative) daughter 2
    break;

  case kAntiLambda:
    fUseCustomAntiLambdaRejectionFilters = true;
    fAntiLambdaRejectionFilters.emplace_back(); //add AliFemtoNSigmaFilter object for (positive) daughter 1
    fAntiLambdaRejectionFilters.emplace_back(); //add AliFemtoNSigmaFilter object for (negative) daughter 2
    break;
      
  default:
    cerr << "E-AliFemtoV0TrackCutNSigmaFilter::CreateCustomV0Rejection: Invalid V0Type"
            "selection '" << aV0TypeToReject << "'.  No rejection filter will be initialized!!!!!" << endl;
    break;
  }
}


void AliFemtoV0TrackCutNSigmaFilter::AddTPCAndTOFNSigmaCutToV0Rejection(AliFemtoV0Type aV0Type, int aDaughterCharge, double aMomMin, double aMomMax, double aNSigmaValueTPC, double aNSigmaValueTOF)
{
  //TODO
  unsigned int iDaughter = -1;
  if(aDaughterCharge>0) iDaughter=0;
  else if(aDaughterCharge<0) iDaughter=1;
  else cerr << "E-AliFemtoV0TrackCutNSigmaFilter::AddTPCAndTOFNSigmaCutToV0Rejection: Invalid DaughterCharge" << endl;

  switch (aV0Type) {
  case kK0s:
    fK0sRejectionFilters[iDaughter].AddTPCAndTOFCut(aMomMin,aMomMax,aNSigmaValueTPC,aNSigmaValueTOF);
    break;

  case kLambda:
    fLambdaRejectionFilters[iDaughter].AddTPCAndTOFCut(aMomMin,aMomMax,aNSigmaValueTPC,aNSigmaValueTOF);
    break;

  case kAntiLambda:
    fAntiLambdaRejectionFilters[iDaughter].AddTPCAndTOFCut(aMomMin,aMomMax,aNSigmaValueTPC,aNSigmaValueTOF);
    break;
      
  default:
    cerr << "E-AliFemtoV0TrackCutNSigmaFilter::AddTPCAndTOFNSigmaCutToV0Rejection: Invalid V0Type"
            "selection '" << aV0Type << "'.  No cut will be initialized!!!!!" << endl;
    break;
  }
}

void AliFemtoV0TrackCutNSigmaFilter::AddTPCAndTOFNSigmaCutToV0Rejection(AliFemtoV0Type aV0Type, double aMomMinPos, double aMomMaxPos, double aNSigmaValueTPCPos, double aNSigmaValueTOFPos,
                                                                                                double aMomMinNeg, double aMomMaxNeg, double aNSigmaValueTPCNeg, double aNSigmaValueTOFNeg)
{
  //TODO
  switch (aV0Type) {
  case kK0s:
    fK0sRejectionFilters[0].AddTPCAndTOFCut(aMomMinPos,aMomMaxPos,aNSigmaValueTPCPos,aNSigmaValueTOFPos);
    fK0sRejectionFilters[1].AddTPCAndTOFCut(aMomMinNeg,aMomMaxNeg,aNSigmaValueTPCNeg,aNSigmaValueTOFNeg);
    break;

  case kLambda:
    fLambdaRejectionFilters[0].AddTPCAndTOFCut(aMomMinPos,aMomMaxPos,aNSigmaValueTPCPos,aNSigmaValueTOFPos);
    fLambdaRejectionFilters[1].AddTPCAndTOFCut(aMomMinNeg,aMomMaxNeg,aNSigmaValueTPCNeg,aNSigmaValueTOFNeg);
    break;

  case kAntiLambda:
    fAntiLambdaRejectionFilters[0].AddTPCAndTOFCut(aMomMinPos,aMomMaxPos,aNSigmaValueTPCPos,aNSigmaValueTOFPos);
    fAntiLambdaRejectionFilters[1].AddTPCAndTOFCut(aMomMinNeg,aMomMaxNeg,aNSigmaValueTPCNeg,aNSigmaValueTOFNeg);
    break;
      
  default:
    cerr << "E-AliFemtoV0TrackCutNSigmaFilter::AddTPCAndTOFNSigmaCutToV0Rejection: Invalid V0Type"
            "selection '" << aV0Type << "'.  No cut will be initialized!!!!!" << endl;
    break;
  }
}


void AliFemtoV0TrackCutNSigmaFilter::AddTPCNSigmaCutToV0Rejection(AliFemtoV0Type aV0Type, int aDaughterCharge, double aMomMin, double aMomMax, double aNSigmaValueTPC)
{
  //TODO
  unsigned int iDaughter = -1;
  if(aDaughterCharge>0) iDaughter=0;
  else if(aDaughterCharge<0) iDaughter=1;
  else cerr << "E-AliFemtoV0TrackCutNSigmaFilter::AddTPCNSigmaCutToV0Rejection: Invalid DaughterCharge" << endl;

  switch (aV0Type) {
  case kK0s:
    fK0sRejectionFilters[iDaughter].AddTPCCut(aMomMin,aMomMax,aNSigmaValueTPC);
    break;

  case kLambda:
    fLambdaRejectionFilters[iDaughter].AddTPCCut(aMomMin,aMomMax,aNSigmaValueTPC);
    break;

  case kAntiLambda:
    fAntiLambdaRejectionFilters[iDaughter].AddTPCCut(aMomMin,aMomMax,aNSigmaValueTPC);
    break;
      
  default:
    cerr << "E-AliFemtoV0TrackCutNSigmaFilter::AddTPCNSigmaCutToV0Rejection: Invalid V0Type"
            "selection '" << aV0Type << "'.  No cut will be initialized!!!!!" << endl;
    break;
  }
}

void AliFemtoV0TrackCutNSigmaFilter::AddTPCNSigmaCutToV0Rejection(AliFemtoV0Type aV0Type, double aMomMinPos, double aMomMaxPos, double aNSigmaValueTPCPos,
                                                                                          double aMomMinNeg, double aMomMaxNeg, double aNSigmaValueTPCNeg)
{
  switch (aV0Type) {
  case kK0s:
    fK0sRejectionFilters[0].AddTPCCut(aMomMinPos,aMomMaxPos,aNSigmaValueTPCPos);
    fK0sRejectionFilters[1].AddTPCCut(aMomMinNeg,aMomMaxNeg,aNSigmaValueTPCNeg);
    break;

  case kLambda:
    fLambdaRejectionFilters[0].AddTPCCut(aMomMinPos,aMomMaxPos,aNSigmaValueTPCPos);
    fLambdaRejectionFilters[1].AddTPCCut(aMomMinNeg,aMomMaxNeg,aNSigmaValueTPCNeg);
    break;

  case kAntiLambda:
    fAntiLambdaRejectionFilters[0].AddTPCCut(aMomMinPos,aMomMaxPos,aNSigmaValueTPCPos);
    fAntiLambdaRejectionFilters[1].AddTPCCut(aMomMinNeg,aMomMaxNeg,aNSigmaValueTPCNeg);
    break;
      
  default:
    cerr << "E-AliFemtoV0TrackCutNSigmaFilter::AddTPCNSigmaCutToV0Rejection: Invalid V0Type"
            "selection '" << aV0Type << "'.  No cut will be initialized!!!!!" << endl;
    break;
  }
}


void AliFemtoV0TrackCutNSigmaFilter::AddTOFNSigmaCutToV0Rejection(AliFemtoV0Type aV0Type, int aDaughterCharge, double aMomMin, double aMomMax, double aNSigmaValueTOF)
{
  //TODO
  unsigned int iDaughter = -1;
  if(aDaughterCharge>0) iDaughter=0;
  else if(aDaughterCharge<0) iDaughter=1;
  else cerr << "E-AliFemtoV0TrackCutNSigmaFilter::AddTOFNSigmaCutToV0Rejection: Invalid DaughterCharge" << endl;

  switch (aV0Type) {
  case kK0s:
    fK0sRejectionFilters[iDaughter].AddTOFCut(aMomMin,aMomMax,aNSigmaValueTOF);
    break;

  case kLambda:
    fLambdaRejectionFilters[iDaughter].AddTOFCut(aMomMin,aMomMax,aNSigmaValueTOF);
    break;

  case kAntiLambda:
    fAntiLambdaRejectionFilters[iDaughter].AddTOFCut(aMomMin,aMomMax,aNSigmaValueTOF);
    break;
      
  default:
    cerr << "E-AliFemtoV0TrackCutNSigmaFilter::AddTOFNSigmaCutToV0Rejection: Invalid V0Type"
            "selection '" << aV0Type << "'.  No cut will be initialized!!!!!" << endl;
    break;
  }
}

void AliFemtoV0TrackCutNSigmaFilter::AddTOFNSigmaCutToV0Rejection(AliFemtoV0Type aV0Type, double aMomMinPos, double aMomMaxPos, double aNSigmaValueTOFPos,
                                                                                          double aMomMinNeg, double aMomMaxNeg, double aNSigmaValueTOFNeg)
{
  switch (aV0Type) {
  case kK0s:
    fK0sRejectionFilters[0].AddTOFCut(aMomMinPos,aMomMaxPos,aNSigmaValueTOFPos);
    fK0sRejectionFilters[1].AddTOFCut(aMomMinNeg,aMomMaxNeg,aNSigmaValueTOFNeg);
    break;

  case kLambda:
    fLambdaRejectionFilters[0].AddTOFCut(aMomMinPos,aMomMaxPos,aNSigmaValueTOFPos);
    fLambdaRejectionFilters[1].AddTOFCut(aMomMinNeg,aMomMaxNeg,aNSigmaValueTOFNeg);
    break;

  case kAntiLambda:
    fAntiLambdaRejectionFilters[0].AddTOFCut(aMomMinPos,aMomMaxPos,aNSigmaValueTOFPos);
    fAntiLambdaRejectionFilters[1].AddTOFCut(aMomMinNeg,aMomMaxNeg,aNSigmaValueTOFNeg);
    break;
      
  default:
    cerr << "E-AliFemtoV0TrackCutNSigmaFilter::AddTOFNSigmaCutToV0Rejection: Invalid V0Type"
            "selection '" << aV0Type << "'.  No cut will be initialized!!!!!" << endl;
    break;
  }
}

/*
bool AliFemtoV0TrackCutNSigmaFilter::IsMisIDK0s(const AliFemtoV0* aV0)
{
  //TODO
  if(fUseCustomK0sRejectionFilters)
  {
    if(aV0->MassK0Short()>fInvMassRejectK0sMin && aV0->MassK0Short()<fInvMassRejectK0sMax)
    {
      if(fK0sRejectionFilters[0].Pass(aV0->PtPos(), aV0->PosNSigmaTPCPi(), aV0->PosNSigmaTOFPi()))
      {
        if(fK0sRejectionFilters[1].Pass(aV0->PtNeg(), aV0->NegNSigmaTPCPi(), aV0->NegNSigmaTOFPi()))
        {
          return true;
        }
      }
    }
    return false;
  }
  else
  {
    return AliFemtoV0TrackCut::IsMisIDK0s(aV0);
  }
}
*/
bool AliFemtoV0TrackCutNSigmaFilter::IsMisIDK0s(const AliFemtoV0* aV0)
{
  //TODO test vs above commented out
  if(fUseCustomK0sRejectionFilters)
  {
    if(aV0->MassK0Short()<fInvMassRejectK0sMin || aV0->MassK0Short()>fInvMassRejectK0sMax) return false;
    if(!fK0sRejectionFilters[0].Pass(aV0->PtPos(), aV0->PosNSigmaTPCPi(), aV0->PosNSigmaTOFPi())) return false;
    if(!fK0sRejectionFilters[1].Pass(aV0->PtNeg(), aV0->NegNSigmaTPCPi(), aV0->NegNSigmaTOFPi())) return false;

    //TODO: Currently testing the addition of the following final cut
    // At this point, the V0 passes as both a (Anti)Lambda (particle of interest) and a K0Short (MisID to be removed)
    // For now, use mass to decide winner
    double tK0ShortMass = 0.497614, tLambdaMass = 1.115683;
    if((fParticleType==kLambda) && (TMath::Abs(aV0->MassLambda()-tLambdaMass) < TMath::Abs(aV0->MassK0Short()-tK0ShortMass))) return false;
    else if((fParticleType==kAntiLambda) && (TMath::Abs(aV0->MassAntiLambda()-tLambdaMass) < TMath::Abs(aV0->MassK0Short()-tK0ShortMass))) return false;
    else {}

    return true;
  }
  else return AliFemtoV0TrackCut::IsMisIDK0s(aV0);
}

bool AliFemtoV0TrackCutNSigmaFilter::IsMisIDLambda(const AliFemtoV0* aV0)
{
  if(fUseCustomLambdaRejectionFilters)
  {
    if(aV0->MassLambda()<fInvMassRejectLambdaMin || aV0->MassLambda()>fInvMassRejectLambdaMax) return false;
    if(!fLambdaRejectionFilters[0].Pass(aV0->PtPos(), aV0->PosNSigmaTPCP(), aV0->PosNSigmaTOFP())) return false;
    if(!fLambdaRejectionFilters[1].Pass(aV0->PtNeg(), aV0->NegNSigmaTPCPi(), aV0->NegNSigmaTOFPi())) return false;

    //TODO: Currently testing the addition of the following final cut
    // At this point, the V0 passes as both a K0Short (particle of interest) and a Lambda (MisID to be removed)
    // For now, use mass to decide winner
    double tK0ShortMass = 0.497614, tLambdaMass = 1.115683;
    if(TMath::Abs(aV0->MassK0Short()-tK0ShortMass) < TMath::Abs(aV0->MassLambda()-tLambdaMass)) return false;

    return true;
  }
  else return AliFemtoV0TrackCut::IsMisIDLambda(aV0);
}
bool AliFemtoV0TrackCutNSigmaFilter::IsMisIDAntiLambda(const AliFemtoV0* aV0)
{
  if(fUseCustomAntiLambdaRejectionFilters)
  {
    if(aV0->MassAntiLambda()<fInvMassRejectAntiLambdaMin || aV0->MassAntiLambda()>fInvMassRejectAntiLambdaMax) return false;
    if(!fAntiLambdaRejectionFilters[0].Pass(aV0->PtPos(), aV0->PosNSigmaTPCPi(), aV0->PosNSigmaTOFPi())) return false;
    if(!fAntiLambdaRejectionFilters[1].Pass(aV0->PtNeg(), aV0->NegNSigmaTPCP(), aV0->NegNSigmaTOFP())) return false;

    //TODO: Currently testing the addition of the following final cut
    // At this point, the V0 passes as both a K0Short (particle of interest) and a AntiLambda (MisID to be removed)
    // For now, use mass to decide winner
    double tK0ShortMass = 0.497614, tLambdaMass = 1.115683;
    if(TMath::Abs(aV0->MassK0Short()-tK0ShortMass) < TMath::Abs(aV0->MassAntiLambda()-tLambdaMass)) return false;

    return true;
  }
  else return AliFemtoV0TrackCut::IsMisIDAntiLambda(aV0);

}

void AliFemtoV0TrackCutNSigmaFilter::SetCTauHistoV0(int aNbins, double aMin, double aMax)
{
  fBuildCTau = true;
  TString tName = "fCTau";
  if(fParticleType==kLambda || fParticleType==kLambdaMC) tName += TString("Lambda");
  else if(fParticleType==kAntiLambda || fParticleType==kAntiLambdaMC) tName += TString("AntiLambda");
  else if(fParticleType==kK0s || fParticleType==kK0sMC) tName += TString("K0Short");
  else {}

  fCTau = new TH1D(tName, "CTau of V0", aNbins, aMin, aMax);
  fCTau->Sumw2();
}

double AliFemtoV0TrackCutNSigmaFilter::GetCTau(const AliFemtoV0* aV0)
{
  double tPtotV0 = aV0->PtotV0();
  double tMassV0;

  //TODO should I use PDG mass values, or aV0->MassLambda(), MassAntiLambda(), MassK0Short()?
  if(fParticleType == kLambda || fParticleType == kAntiLambda) tMassV0 = 1.115683;
  else if(fParticleType == kK0s) tMassV0 = 0.497614;
  else tMassV0 = 0.;

  double tCTau = aV0->DecayLengthV0()*tMassV0/tPtotV0;

  return tCTau;
}

TList *AliFemtoV0TrackCutNSigmaFilter::GetOutputList()
{
  TList* tOutputList = AliFemtoV0TrackCut::GetOutputList();  //add all of typical objects
  if(fBuildCTau) tOutputList->Add(fCTau);

  return tOutputList;
}

