///
/// \file AliFemtoXiTrackCutNSigmaFilter.cxx
///


#include "AliFemtoXiTrackCutNSigmaFilter.h"

#ifdef __ROOT__
  /// \cond CLASSIMP
  ClassImp(AliFemtoXiTrackCutNSigmaFilter);
  /// \endcond
#endif


AliFemtoXiTrackCutNSigmaFilter::AliFemtoXiTrackCutNSigmaFilter() : 
  AliFemtoXiTrackCut(),

  fUseCustomPionNSigmaFilter(false),
  fUseCustomKaonNSigmaFilter(false),
  fUseCustomProtonNSigmaFilter(false),
  fUseCustomBacPionNSigmaFilter(false),

  fPionNSigmaFilter(NULL),
  fKaonNSigmaFilter(NULL),
  fProtonNSigmaFilter(NULL),
  fBacPionNSigmaFilter(NULL),

  fUseCustomK0sRejectionFilters(false),
  fUseCustomLambdaRejectionFilters(false),
  fUseCustomAntiLambdaRejectionFilters(false),

  fK0sRejectionFilters(0),
  fLambdaRejectionFilters(0),
  fAntiLambdaRejectionFilters(0)
{
  // Default constructor
}


AliFemtoXiTrackCutNSigmaFilter::AliFemtoXiTrackCutNSigmaFilter(const AliFemtoXiTrackCutNSigmaFilter& aCut) :
  AliFemtoXiTrackCut(aCut),

  fUseCustomPionNSigmaFilter(aCut.fUseCustomPionNSigmaFilter),
  fUseCustomKaonNSigmaFilter(aCut.fUseCustomKaonNSigmaFilter),
  fUseCustomProtonNSigmaFilter(aCut.fUseCustomProtonNSigmaFilter),
  fUseCustomBacPionNSigmaFilter(aCut.fUseCustomBacPionNSigmaFilter),

  fUseCustomK0sRejectionFilters(aCut.fUseCustomK0sRejectionFilters),
  fUseCustomLambdaRejectionFilters(aCut.fUseCustomLambdaRejectionFilters),
  fUseCustomAntiLambdaRejectionFilters(aCut.fUseCustomAntiLambdaRejectionFilters),

  fK0sRejectionFilters(aCut.fK0sRejectionFilters),
  fLambdaRejectionFilters(aCut.fLambdaRejectionFilters),
  fAntiLambdaRejectionFilters(aCut.fAntiLambdaRejectionFilters)
{
  if(aCut.fPionNSigmaFilter) fPionNSigmaFilter = new AliFemtoNSigmaFilter(*aCut.fPionNSigmaFilter);
  if(aCut.fKaonNSigmaFilter) fKaonNSigmaFilter = new AliFemtoNSigmaFilter(*aCut.fKaonNSigmaFilter);
  if(aCut.fProtonNSigmaFilter) fProtonNSigmaFilter = new AliFemtoNSigmaFilter(*aCut.fProtonNSigmaFilter);
  if(aCut.fBacPionNSigmaFilter) fBacPionNSigmaFilter = new AliFemtoNSigmaFilter(*aCut.fBacPionNSigmaFilter);

  if(aCut.fMinvPurityAidHistoXi) fMinvPurityAidHistoXi = new TH1D(*aCut.fMinvPurityAidHistoXi);
  else fMinvPurityAidHistoXi = 0;
}

AliFemtoXiTrackCutNSigmaFilter& AliFemtoXiTrackCutNSigmaFilter::operator=(const AliFemtoXiTrackCutNSigmaFilter& aCut)
{
  //assignment operator
  if (this == &aCut) return *this;

  AliFemtoXiTrackCut::operator=(aCut);

  fUseCustomPionNSigmaFilter = aCut.fUseCustomPionNSigmaFilter;
  fUseCustomKaonNSigmaFilter = aCut.fUseCustomKaonNSigmaFilter;
  fUseCustomProtonNSigmaFilter = aCut.fUseCustomProtonNSigmaFilter;
  fUseCustomBacPionNSigmaFilter = aCut.fUseCustomBacPionNSigmaFilter;

  fUseCustomK0sRejectionFilters = aCut.fUseCustomK0sRejectionFilters;
  fUseCustomLambdaRejectionFilters = aCut.fUseCustomLambdaRejectionFilters;
  fUseCustomAntiLambdaRejectionFilters = aCut.fUseCustomAntiLambdaRejectionFilters;

  fK0sRejectionFilters = aCut.fK0sRejectionFilters;
  fLambdaRejectionFilters = aCut.fLambdaRejectionFilters;
  fAntiLambdaRejectionFilters = aCut.fAntiLambdaRejectionFilters;

  if(aCut.fPionNSigmaFilter) fPionNSigmaFilter = new AliFemtoNSigmaFilter(*aCut.fPionNSigmaFilter);
  if(aCut.fKaonNSigmaFilter) fKaonNSigmaFilter = new AliFemtoNSigmaFilter(*aCut.fKaonNSigmaFilter);
  if(aCut.fProtonNSigmaFilter) fProtonNSigmaFilter = new AliFemtoNSigmaFilter(*aCut.fProtonNSigmaFilter);
  if(aCut.fBacPionNSigmaFilter) fBacPionNSigmaFilter = new AliFemtoNSigmaFilter(*aCut.fBacPionNSigmaFilter);

  if(aCut.fMinvPurityAidHistoXi) fMinvPurityAidHistoXi = new TH1D(*aCut.fMinvPurityAidHistoXi);
  else fMinvPurityAidHistoXi = 0;

  return *this;
}

AliFemtoXiTrackCutNSigmaFilter* AliFemtoXiTrackCutNSigmaFilter::Clone()
{
  return(new AliFemtoXiTrackCutNSigmaFilter(*this));
}

AliFemtoXiTrackCutNSigmaFilter::~AliFemtoXiTrackCutNSigmaFilter()
{
  delete fPionNSigmaFilter;
  delete fKaonNSigmaFilter;
  delete fProtonNSigmaFilter;
  delete fBacPionNSigmaFilter;

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


bool AliFemtoXiTrackCutNSigmaFilter::PassV0(const AliFemtoXi* aXi)
{
  // test the particle and return
  // true if it meets all the criteria
  // false if it doesn't meet at least one of the criteria

  Float_t pt = aXi->PtV0();
  Float_t eta = aXi->EtaV0();

  //kinematic cuts
  if (TMath::Abs(eta) > fEta) return false;  //put in kinematic cuts by hand
  if (pt < fPtMin || fPtMax < pt) return false;
  if (TMath::Abs(aXi->EtaPos()) > fMaxEtaDaughters) return false;
  if (TMath::Abs(aXi->EtaNeg()) > fMaxEtaDaughters) return false;

  if (aXi->PtPos() < fPtMinPosDaughter) return false;
  if (aXi->PtNeg() < fPtMinNegDaughter) return false;
  if (aXi->PtPos() > fPtMaxPosDaughter) return false;
  if (aXi->PtNeg() > fPtMaxNegDaughter) return false;

  //V0 from kinematics information
  if (fParticleType == kLambdaMC ) {
    if (!(aXi->MassLambda() > fInvMassLambdaMin && aXi->MassLambda() < fInvMassLambdaMax) || !(aXi->PosNSigmaTPCP() == 0)) {
      return false;
    } else {
      return true;
    }
  } else if (fParticleType == kAntiLambdaMC) {
    if (!(aXi->MassLambda() > fInvMassLambdaMin && aXi->MassLambda() < fInvMassLambdaMax) || !(aXi->NegNSigmaTPCP() == 0)) {
      return false;
    } else {
      return true;
    }
  } else if (fParticleType == kK0sMC) {
    if (!(aXi->MassK0Short() > fInvMassK0sMin && aXi->MassK0Short() < fInvMassK0sMax) || !(aXi->NegNSigmaTPCP() == 0))
      return false;
    else {
      return true;
    }
  }



  //quality cuts
  if (aXi->OnFlyStatusV0() != fOnFlyStatus) return false;
  if (aXi->StatusNeg() == 999 || aXi->StatusPos() == 999) return false;
  if (aXi->TPCNclsPos() < fTPCNclsDaughters) return false;
  if (aXi->TPCNclsNeg() < fTPCNclsDaughters) return false;
  if (aXi->NdofPos() > fNdofDaughters) return false;
  if (aXi->NdofNeg() > fNdofDaughters) return false;
  if (!(aXi->StatusNeg() & fStatusDaughters)) return false;
  if (!(aXi->StatusPos() & fStatusDaughters)) return false;

  //fiducial volume radius
  if(aXi->RadiusV0()<fRadiusV0Min || aXi->RadiusV0()>fRadiusV0Max)
    return false;

  //DCA between daughter particles
  if (TMath::Abs(aXi->DcaV0Daughters()) > fMaxDcaV0Daughters)
    return false;

  //DCA of daughters to primary vertex
  if (TMath::Abs(aXi->DcaPosToPrimVertex()) < fMinDcaDaughterPosToVert || TMath::Abs(aXi->DcaNegToPrimVertex()) < fMinDcaDaughterNegToVert)
    return false;

  //DCA V0 to prim vertex
  if (TMath::Abs(aXi->DcaV0ToPrimVertex()) > fMaxDcaV0 || TMath::Abs(aXi->DcaV0ToPrimVertex()) < fMinDcaV0)
    return false;

  //becomes obsolete - wrong name of the data memeber and the corresnponding methods (by default is set to fMaxCosPointingAngle = 0.0)
  //cos pointing angle
  if (aXi->CosPointingAngle() < fMaxCosPointingAngle)
    return false;

  //this is the correct name of the data member and the corresponding methods (we are accepting cos(pointing angle bigger than certain minimum)
  //cos pointing angle
  if (aXi->CosPointingAngle() < fMinCosPointingAngle)
    return false;

  //decay length
  if (aXi->DecayLengthV0() > fMaxDecayLength)
    return false;

  //transverse decay vertex
  //TO BE IMPLEMENTED (maybe in future)

  if (fParticleType == kAll)
    return true;


  bool pid_check = false;
  // Looking for lambdas = proton + pim
  if(fParticleType == kLambda) 
  {
    if(IsProtonNSigma(aXi->PtPos(), aXi->PosNSigmaTPCP(), aXi->PosNSigmaTOFP(),fNsigmaPosDaughterTPC,fNsigmaPosDaughterTOF,fRequireTOFProton)) //proton
    {
      if (IsPionNSigma(aXi->PtNeg(), aXi->NegNSigmaTPCPi(), aXi->NegNSigmaTOFPi(),fNsigmaNegDaughterTPC,fNsigmaNegDaughterTOF,fRequireTOFPion)) //pion-
      {
        pid_check = true;
        if(fUseLooseInvMassCut)
        {
          if(aXi->MassLambda()<fLooseInvMassMin || aXi->MassLambda()>fLooseInvMassMax) return false;
        }
        if(fRemoveMisidentified)
        {
          if(IsMisIDK0s(aXi))
          {
            if(fBuildMisIDHistograms)
            {
              fLambdaMassOfMisIDV0->Fill(aXi->MassLambda());
              fK0sMassOfMisIDV0->Fill(aXi->MassK0Short());
            }
            return false;
          }
        }
        if(fBuildPurityAidV0) fMinvPurityAidHistoV0->Fill(aXi->MassLambda());
        //invariant mass lambda
        if(aXi->MassLambda() < fInvMassLambdaMin || aXi->MassLambda() > fInvMassLambdaMax) return false;
      }
    }
  }

  //Looking for antilambdas =  antiproton + pip
  else if(fParticleType == kAntiLambda)
  {
    if(IsProtonNSigma(aXi->PtNeg(), aXi->NegNSigmaTPCP(), aXi->NegNSigmaTOFP(),fNsigmaNegDaughterTPC,fNsigmaNegDaughterTOF,fRequireTOFProton)) //antiproton
    {
      if(IsPionNSigma(aXi->PtPos(), aXi->PosNSigmaTPCPi(), aXi->PosNSigmaTOFPi(),fNsigmaPosDaughterTPC,fNsigmaPosDaughterTOF,fRequireTOFPion)) //pion+
      {
        pid_check = true;
        if(fUseLooseInvMassCut)
        {
          if(aXi->MassAntiLambda()<fLooseInvMassMin || aXi->MassAntiLambda()>fLooseInvMassMax) return false;
        }
        if(fRemoveMisidentified)
        {
          if(IsMisIDK0s(aXi))
          {
            if(fBuildMisIDHistograms)
            {
              fAntiLambdaMassOfMisIDV0->Fill(aXi->MassAntiLambda());
              fK0sMassOfMisIDV0->Fill(aXi->MassK0Short());
            }
            return false;
          }
        }
        if(fBuildPurityAidV0) fMinvPurityAidHistoV0->Fill(aXi->MassAntiLambda());
        //invariant mass antilambda
        if (aXi->MassAntiLambda() < fInvMassLambdaMin || aXi->MassAntiLambda() > fInvMassLambdaMax) return false;
      }
    }
  }

  //Looking for K0s = pip + pim
  else if (fParticleType == kK0s)
  {
    if(IsPionNSigma(aXi->PtNeg(), aXi->NegNSigmaTPCPi(), aXi->NegNSigmaTOFPi(),fNsigmaNegDaughterTPC,fNsigmaNegDaughterTOF,fRequireTOFPion)) //pion-
    {
      if(IsPionNSigma(aXi->PtPos(), aXi->PosNSigmaTPCPi(), aXi->PosNSigmaTOFPi(),fNsigmaPosDaughterTPC,fNsigmaPosDaughterTOF,fRequireTOFPion)) //pion+
      {
        pid_check = true;
        if(fUseLooseInvMassCut)
        {
          if(aXi->MassK0Short()<fLooseInvMassMin || aXi->MassK0Short()>fLooseInvMassMax) return false;
        }
        if(fRemoveMisidentified)
        {
          if(IsMisIDLambda(aXi))
          {
            if(fBuildMisIDHistograms)
            {
              fK0sMassOfMisIDV0->Fill(aXi->MassK0Short());
              fLambdaMassOfMisIDV0->Fill(aXi->MassLambda());
            }
            return false;
          }
          else if(IsMisIDAntiLambda(aXi))
          {
            if(fBuildMisIDHistograms)
            {
              fK0sMassOfMisIDV0->Fill(aXi->MassK0Short());
              fAntiLambdaMassOfMisIDV0->Fill(aXi->MassAntiLambda());
            }
            return false;
          }
        }
        if(fBuildPurityAidV0) fMinvPurityAidHistoV0->Fill(aXi->MassK0Short());
        //invariant mass K0s
        if (aXi->MassK0Short() < fInvMassK0sMin || aXi->MassK0Short() > fInvMassK0sMax) return false;
      }
    }
  }

  if (!pid_check) return false;
  return true;
}


bool AliFemtoXiTrackCutNSigmaFilter::Pass(const AliFemtoXi* aXi)
{
  // test the particle and return 
  // true if it meets all the criteria
  // false if it doesn't meet at least one of the criteria
   
  Float_t pt = aXi->PtXi();
  Float_t eta = aXi->EtaXi();
  
  if(aXi->ChargeXi()==0)
    return false;


  //ParticleType selection
  //If fParticleTypeXi=kAll, any charged candidate will pass 
  if(fParticleTypeXi == AliFemtoXiTrackCut::kXiPlus && aXi->ChargeXi() == -1) 
    return false;

  if(fParticleTypeXi == AliFemtoXiTrackCut::kXiMinus && aXi->ChargeXi() == 1) 
    return false;

  //kinematic cuts
  if(TMath::Abs(eta) > fMaxEtaXi) return false;  //put in kinematic cuts by hand
  if(pt < fMinPtXi) return false;
  if(pt > fMaxPtXi) return false;
  if(TMath::Abs(aXi->EtaBac()) > fMaxEtaBac) return false;
  if(aXi->PtBac()< fMinPtBac) return false; 
  if(aXi->PtBac()> fMaxPtBac) return false; 

  

    //Xi from kinematics information
    if (fParticleTypeXi == AliFemtoXiTrackCut::kXiMinusMC || fParticleTypeXi == AliFemtoXiTrackCut::kXiPlusMC) {
      if(!(aXi->MassXi()>fInvMassXiMin && aXi->MassXi()<fInvMassXiMax) || !(aXi->BacNSigmaTPCPi()==0))
	return false; 
      else
	{
	  return true;  
	} 
    }


    //quality cuts
    if(aXi->StatusBac() == 999) return false;
    if(aXi->TPCNclsBac()<fTPCNclsBac) return false;
    if(aXi->NdofBac()>fNdofBac) return false;
    if(!(aXi->StatusBac()&fStatusBac)) return false;


    //DCA Xi to prim vertex
    if(TMath::Abs(aXi->DcaXiToPrimVertex())>fMaxDcaXi)
      return false;

    //DCA Xi bachelor to prim vertex
    if(TMath::Abs(aXi->DcaBacToPrimVertex())<fMinDcaXiBac)
      return false;

    //DCA Xi daughters
    if(TMath::Abs(aXi->DcaXiDaughters())>fMaxDcaXiDaughters)
      return false;
    
    //cos pointing angle
    if(aXi->CosPointingAngleXi()<fMinCosPointingAngleXi)
      return false; 
    
    //decay length
    if(aXi->DecayLengthXi()>fMaxDecayLengthXi)
      return false;
    
 
  if(fParticleTypeXi == kAll)
    return true;



  bool pid_check=false;
  // Looking for Xi
  if (fParticleTypeXi == AliFemtoXiTrackCut::kXiMinus || fParticleTypeXi == AliFemtoXiTrackCut::kXiPlus) {
    if (IsBacPionNSigma(aXi->PtBac(), aXi->BacNSigmaTPCPi(), aXi->BacNSigmaTOFPi())) //bachelor pion
	{
	  pid_check=true;
	}

  }

  if (!pid_check) return false;
  
  if(!PassV0(aXi))
    return false;

  if(fBuildPurityAidXi) {fMinvPurityAidHistoXi->Fill(aXi->MassXi());}

   //invariant mass Xi
  if(aXi->MassXi()<fInvMassXiMin || aXi->MassXi()>fInvMassXiMax)
   {
     return false;
   }  
  
  return true;
}


AliFemtoString AliFemtoXiTrackCutNSigmaFilter::Report()
{
  TString report;
  report += TString::Format("Usings custom Pion NSigma Filter:\t%i\n", fUseCustomPionNSigmaFilter)
          + TString::Format("Usings custom Kaon NSigma Filter:\t%i\n", fUseCustomKaonNSigmaFilter)
          + TString::Format("Usings custom Proton NSigma Filter:\t%i\n", fUseCustomProtonNSigmaFilter)
          + TString::Format("Usings custom BacPion NSigma Filter:\t%i\n", fUseCustomBacPionNSigmaFilter)
          + AliFemtoXiTrackCut::Report()
          + AliFemtoV0TrackCut::Report();
  return AliFemtoString(report);
}

TList *AliFemtoXiTrackCutNSigmaFilter::ListSettings()
{
  // return a list of settings in a writable form
  TList *tListSetttings = new TList();
 
  return tListSetttings;
}

TList *AliFemtoXiTrackCutNSigmaFilter::GetOutputList()
{
  TList *tOutputList = AliFemtoCutMonitorHandler::GetOutputList();  //add all of the typical objects

  if(fBuildPurityAidXi) tOutputList->Add(fMinvPurityAidHistoXi);
  if(fBuildPurityAidV0) tOutputList->Add(fMinvPurityAidHistoV0);

  return tOutputList;
}

bool AliFemtoXiTrackCutNSigmaFilter::IsKaonNSigma(float mom, float nsigmaTPCK, float nsigmaTOFK)
{
  if (fUseCustomKaonNSigmaFilter) {
    return fKaonNSigmaFilter->Pass(mom, nsigmaTPCK, nsigmaTOFK);
  } else {
    return AliFemtoV0TrackCut::IsKaonNSigma(mom, nsigmaTPCK, nsigmaTOFK);
  }
}


bool AliFemtoXiTrackCutNSigmaFilter::IsPionNSigma(float mom, float nsigmaTPCPi, float nsigmaTOFPi, double nsigmacutTPC, double nsigmacutTOF, bool requireTOF)
{
  if (fUseCustomPionNSigmaFilter) {
    return fPionNSigmaFilter->Pass(mom, nsigmaTPCPi, nsigmaTOFPi);
  } else {
    return AliFemtoV0TrackCut::IsPionNSigma(mom,nsigmaTPCPi,nsigmaTOFPi,nsigmacutTPC,nsigmacutTOF,requireTOF);
  }
}

bool AliFemtoXiTrackCutNSigmaFilter::IsProtonNSigma(float mom, float nsigmaTPCP, float nsigmaTOFP, double nsigmacutTPC, double nsigmacutTOF, bool requireTOF)
{
  if (fUseCustomProtonNSigmaFilter) {
    return fProtonNSigmaFilter->Pass(mom, nsigmaTPCP, nsigmaTOFP);
  } else {
    return AliFemtoV0TrackCut::IsProtonNSigma(mom,nsigmaTPCP,nsigmaTOFP,nsigmacutTPC,nsigmacutTOF,requireTOF);
  }
}

bool AliFemtoXiTrackCutNSigmaFilter::IsBacPionNSigma(float mom, float nsigmaTPCPi, float nsigmaTOFPi, double nsigmacutTPC, double nsigmacutTOF, bool requireTOF)
{
  if (fUseCustomBacPionNSigmaFilter) {
    return fBacPionNSigmaFilter->Pass(mom, nsigmaTPCPi, nsigmaTOFPi);
  } else {
    return AliFemtoV0TrackCut::IsPionNSigma(mom,nsigmaTPCPi,nsigmaTOFPi,nsigmacutTPC,nsigmacutTOF,requireTOF);
  }
}

void AliFemtoXiTrackCutNSigmaFilter::CreateCustomNSigmaFilter(XiDaughterParticleType aDaughterType)
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

  case kBacPion:
    fUseCustomBacPionNSigmaFilter = true;
    fBacPionNSigmaFilter = new AliFemtoNSigmaFilter();
    break;

  default:
    cerr << "E-AliFemtoXiTrackCutNSigmaFilter::CreateCustomNSigmaFilter: Invalid DaughterParticleType"
            "selection '" << aDaughterType << "'.  No custom filter will be initialized!!!!!" << endl;
    break;
  }
}

void AliFemtoXiTrackCutNSigmaFilter::AddTPCAndTOFNSigmaCut(XiDaughterParticleType aDaughterType, double aMomMin, double aMomMax, double aNSigmaValueTPC, double aNSigmaValueTOF)
{
  if(aDaughterType == kPion) {fPionNSigmaFilter->AddTPCAndTOFCut(aMomMin,aMomMax,aNSigmaValueTPC,aNSigmaValueTOF);}
  else if(aDaughterType == kKaon) {fKaonNSigmaFilter->AddTPCAndTOFCut(aMomMin,aMomMax,aNSigmaValueTPC,aNSigmaValueTOF);}
  else if(aDaughterType == kProton) {fProtonNSigmaFilter->AddTPCAndTOFCut(aMomMin,aMomMax,aNSigmaValueTPC,aNSigmaValueTOF);}
  else if(aDaughterType==kBacPion) fBacPionNSigmaFilter->AddTPCAndTOFCut(aMomMin,aMomMax,aNSigmaValueTPC,aNSigmaValueTOF);
  else {/*cerr << "Invalid DaughterParticleType selection in AddTPCAndTOFNSigmaCut.  No cut will be added to the filter!!!!!" << endl;*/}
}

void AliFemtoXiTrackCutNSigmaFilter::AddTPCNSigmaCut(XiDaughterParticleType aDaughterType, double aMomMin, double aMomMax, double aNSigmaValueTPC)
{
  if(aDaughterType == kPion) {fPionNSigmaFilter->AddTPCCut(aMomMin,aMomMax,aNSigmaValueTPC);}
  else if(aDaughterType == kKaon) {fKaonNSigmaFilter->AddTPCCut(aMomMin,aMomMax,aNSigmaValueTPC);}
  else if(aDaughterType == kProton) {fProtonNSigmaFilter->AddTPCCut(aMomMin,aMomMax,aNSigmaValueTPC);}
  else if(aDaughterType==kBacPion) fBacPionNSigmaFilter->AddTPCCut(aMomMin,aMomMax,aNSigmaValueTPC);
  else {/*cerr << "Invalid DaughterParticleType selection in AddTPCNSigmaCut.  No cut will be added to the filter!!!!!" << endl;*/}
}

void AliFemtoXiTrackCutNSigmaFilter::AddTOFNSigmaCut(XiDaughterParticleType aDaughterType, double aMomMin, double aMomMax, double aNSigmaValueTOF)
{
  if(aDaughterType == kPion) {fPionNSigmaFilter->AddTOFCut(aMomMin,aMomMax,aNSigmaValueTOF);}
  else if(aDaughterType == kKaon) {fKaonNSigmaFilter->AddTOFCut(aMomMin,aMomMax,aNSigmaValueTOF);}
  else if(aDaughterType == kProton) {fProtonNSigmaFilter->AddTOFCut(aMomMin,aMomMax,aNSigmaValueTOF);}
  else if(aDaughterType==kBacPion) fBacPionNSigmaFilter->AddTOFCut(aMomMin,aMomMax,aNSigmaValueTOF);
  else {/*cerr << "Invalid DaughterParticleType selection in AddTOFNSigmaCut.  No cut will be added to the filter!!!!!" << endl;*/}
}



void AliFemtoXiTrackCutNSigmaFilter::CreateCustomV0Rejection(AliFemtoV0Type aV0TypeToReject)
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
    cerr << "E-AliFemtoXiTrackCutNSigmaFilter::CreateCustomV0Rejection: Invalid V0Type"
            "selection '" << aV0TypeToReject << "'.  No rejection filter will be initialized!!!!!" << endl;
    break;
  }
}


void AliFemtoXiTrackCutNSigmaFilter::AddTPCAndTOFNSigmaCutToV0Rejection(AliFemtoV0Type aV0Type, int aDaughterCharge, double aMomMin, double aMomMax, double aNSigmaValueTPC, double aNSigmaValueTOF)
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
    cerr << "E-AliFemtoXiTrackCutNSigmaFilter::AddTPCAndTOFNSigmaCutToV0Rejection: Invalid V0Type"
            "selection '" << aV0Type << "'.  No cut will be initialized!!!!!" << endl;
    break;
  }
}

void AliFemtoXiTrackCutNSigmaFilter::AddTPCAndTOFNSigmaCutToV0Rejection(AliFemtoV0Type aV0Type, double aMomMinPos, double aMomMaxPos, double aNSigmaValueTPCPos, double aNSigmaValueTOFPos,
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
    cerr << "E-AliFemtoXiTrackCutNSigmaFilter::AddTPCAndTOFNSigmaCutToV0Rejection: Invalid V0Type"
            "selection '" << aV0Type << "'.  No cut will be initialized!!!!!" << endl;
    break;
  }
}


void AliFemtoXiTrackCutNSigmaFilter::AddTPCNSigmaCutToV0Rejection(AliFemtoV0Type aV0Type, int aDaughterCharge, double aMomMin, double aMomMax, double aNSigmaValueTPC)
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
    cerr << "E-AliFemtoXiTrackCutNSigmaFilter::AddTPCNSigmaCutToV0Rejection: Invalid V0Type"
            "selection '" << aV0Type << "'.  No cut will be initialized!!!!!" << endl;
    break;
  }
}

void AliFemtoXiTrackCutNSigmaFilter::AddTPCNSigmaCutToV0Rejection(AliFemtoV0Type aV0Type, double aMomMinPos, double aMomMaxPos, double aNSigmaValueTPCPos,
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
    cerr << "E-AliFemtoXiTrackCutNSigmaFilter::AddTPCNSigmaCutToV0Rejection: Invalid V0Type"
            "selection '" << aV0Type << "'.  No cut will be initialized!!!!!" << endl;
    break;
  }
}


void AliFemtoXiTrackCutNSigmaFilter::AddTOFNSigmaCutToV0Rejection(AliFemtoV0Type aV0Type, int aDaughterCharge, double aMomMin, double aMomMax, double aNSigmaValueTOF)
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
    cerr << "E-AliFemtoXiTrackCutNSigmaFilter::AddTOFNSigmaCutToV0Rejection: Invalid V0Type"
            "selection '" << aV0Type << "'.  No cut will be initialized!!!!!" << endl;
    break;
  }
}

void AliFemtoXiTrackCutNSigmaFilter::AddTOFNSigmaCutToV0Rejection(AliFemtoV0Type aV0Type, double aMomMinPos, double aMomMaxPos, double aNSigmaValueTOFPos,
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
    cerr << "E-AliFemtoXiTrackCutNSigmaFilter::AddTOFNSigmaCutToV0Rejection: Invalid V0Type"
            "selection '" << aV0Type << "'.  No cut will be initialized!!!!!" << endl;
    break;
  }
}

/*
bool AliFemtoXiTrackCutNSigmaFilter::IsMisIDK0s(const AliFemtoV0* aV0)
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
bool AliFemtoXiTrackCutNSigmaFilter::IsMisIDK0s(const AliFemtoXi* aXi)
{
  //TODO test vs above commented out
  if(fUseCustomK0sRejectionFilters)
  {
    if(aXi->MassK0Short()<fInvMassRejectK0sMin || aXi->MassK0Short()>fInvMassRejectK0sMax) return false;
    if(!fK0sRejectionFilters[0].Pass(aXi->PtPos(), aXi->PosNSigmaTPCPi(), aXi->PosNSigmaTOFPi())) return false;
    if(!fK0sRejectionFilters[1].Pass(aXi->PtNeg(), aXi->NegNSigmaTPCPi(), aXi->NegNSigmaTOFPi())) return false;

    //TODO: Currently testing the addition of the following final cut
    // At this point, the V0 passes as both a (Anti)Lambda (particle of interest) and a K0Short (MisID to be removed)
    // For now, use mass to decide winner
    double tK0ShortMass = 0.497614, tLambdaMass = 1.115683;
    if((fParticleType==kLambda) && (TMath::Abs(aXi->MassLambda()-tLambdaMass) < TMath::Abs(aXi->MassK0Short()-tK0ShortMass))) return false;
    else if((fParticleType==kAntiLambda) && (TMath::Abs(aXi->MassAntiLambda()-tLambdaMass) < TMath::Abs(aXi->MassK0Short()-tK0ShortMass))) return false;
    else {}

    return true;
  }
  else return AliFemtoV0TrackCut::IsMisIDK0s(aXi);
}

bool AliFemtoXiTrackCutNSigmaFilter::IsMisIDLambda(const AliFemtoXi* aXi)
{
  if(fUseCustomLambdaRejectionFilters)
  {
    if(aXi->MassLambda()<fInvMassRejectLambdaMin || aXi->MassLambda()>fInvMassRejectLambdaMax) return false;
    if(!fLambdaRejectionFilters[0].Pass(aXi->PtPos(), aXi->PosNSigmaTPCP(), aXi->PosNSigmaTOFP())) return false;
    if(!fLambdaRejectionFilters[1].Pass(aXi->PtNeg(), aXi->NegNSigmaTPCPi(), aXi->NegNSigmaTOFPi())) return false;

    //TODO: Currently testing the addition of the following final cut
    // At this point, the V0 passes as both a K0Short (particle of interest) and a Lambda (MisID to be removed)
    // For now, use mass to decide winner
    double tK0ShortMass = 0.497614, tLambdaMass = 1.115683;
    if(TMath::Abs(aXi->MassK0Short()-tK0ShortMass) < TMath::Abs(aXi->MassLambda()-tLambdaMass)) return false;

    return true;
  }
  else return AliFemtoV0TrackCut::IsMisIDLambda(aXi);
}
bool AliFemtoXiTrackCutNSigmaFilter::IsMisIDAntiLambda(const AliFemtoXi* aXi)
{
  if(fUseCustomAntiLambdaRejectionFilters)
  {
    if(aXi->MassAntiLambda()<fInvMassRejectAntiLambdaMin || aXi->MassAntiLambda()>fInvMassRejectAntiLambdaMax) return false;
    if(!fAntiLambdaRejectionFilters[0].Pass(aXi->PtPos(), aXi->PosNSigmaTPCPi(), aXi->PosNSigmaTOFPi())) return false;
    if(!fAntiLambdaRejectionFilters[1].Pass(aXi->PtNeg(), aXi->NegNSigmaTPCP(), aXi->NegNSigmaTOFP())) return false;

    //TODO: Currently testing the addition of the following final cut
    // At this point, the V0 passes as both a K0Short (particle of interest) and a AntiLambda (MisID to be removed)
    // For now, use mass to decide winner
    double tK0ShortMass = 0.497614, tLambdaMass = 1.115683;
    if(TMath::Abs(aXi->MassK0Short()-tK0ShortMass) < TMath::Abs(aXi->MassAntiLambda()-tLambdaMass)) return false;

    return true;
  }
  else return AliFemtoV0TrackCut::IsMisIDAntiLambda(aXi);

}

