///
/// \file AliFemtoXiTrackCutNSigmaFilter.cxx
///


#include "AliFemtoXiTrackCutNSigmaFilter.h"
#include <TObjArray.h>

#ifdef __ROOT__
  /// \cond CLASSIMP
  ClassImp(AliFemtoXiTrackCutNSigmaFilter);
  /// \endcond
#endif


AliFemtoXiTrackCutNSigmaFilter::AliFemtoXiTrackCutNSigmaFilter() : 
  AliFemtoXiTrackCut(),

  fIsDaughterV0FilterUpdated(false),
  fIsDaughterBacPionFilterUpdated(false),

  fDaughterV0Filter(NULL),
  fDaughterBacPionFilter(NULL),


  fUseCustomBacPionNSigmaFilter(false),
  fBacPionNSigmaFilter(NULL),

  fCTauXi(0)
{
  // Default constructor
  fDaughterV0Filter = new AliFemtoV0TrackCutNSigmaFilter();
  fCTauXi = new TH1D("fCTauXi", "CTau of Xi", 500, 0.0, 50.0);
  fCTauXi->Sumw2();
}


AliFemtoXiTrackCutNSigmaFilter::AliFemtoXiTrackCutNSigmaFilter(const AliFemtoXiTrackCutNSigmaFilter& aCut) :
  AliFemtoXiTrackCut(aCut),

  fIsDaughterV0FilterUpdated(aCut.fIsDaughterV0FilterUpdated),
  fDaughterBacPionFilter(aCut.fDaughterBacPionFilter),
  fUseCustomBacPionNSigmaFilter(aCut.fUseCustomBacPionNSigmaFilter)
{
  if(aCut.fDaughterV0Filter) fDaughterV0Filter = new AliFemtoV0TrackCutNSigmaFilter(*aCut.fDaughterV0Filter);
  if(aCut.fDaughterBacPionFilter) fDaughterBacPionFilter = new AliFemtoESDTrackCutNSigmaFilter(*aCut.fDaughterBacPionFilter);

  if(aCut.fBacPionNSigmaFilter) fBacPionNSigmaFilter = new AliFemtoNSigmaFilter(*aCut.fBacPionNSigmaFilter);

  if(aCut.fMinvPurityAidHistoXi) fMinvPurityAidHistoXi = new TH1D(*aCut.fMinvPurityAidHistoXi);
  else fMinvPurityAidHistoXi = 0;

  if(aCut.fCTauXi) fCTauXi = new TH1D(*aCut.fCTauXi);
  else fCTauXi = 0;
}

AliFemtoXiTrackCutNSigmaFilter& AliFemtoXiTrackCutNSigmaFilter::operator=(const AliFemtoXiTrackCutNSigmaFilter& aCut)
{
  //assignment operator
  if (this == &aCut) return *this;

  AliFemtoXiTrackCut::operator=(aCut);

  fIsDaughterV0FilterUpdated = aCut.fIsDaughterV0FilterUpdated;
  fDaughterBacPionFilter = aCut.fDaughterBacPionFilter;
  fUseCustomBacPionNSigmaFilter = aCut.fUseCustomBacPionNSigmaFilter;

  if(aCut.fDaughterV0Filter) fDaughterV0Filter = new AliFemtoV0TrackCutNSigmaFilter(*aCut.fDaughterV0Filter);
  if(aCut.fDaughterBacPionFilter) fDaughterBacPionFilter = new AliFemtoESDTrackCutNSigmaFilter(*aCut.fDaughterBacPionFilter);

  if(aCut.fBacPionNSigmaFilter) fBacPionNSigmaFilter = new AliFemtoNSigmaFilter(*aCut.fBacPionNSigmaFilter);

  if(aCut.fMinvPurityAidHistoXi) fMinvPurityAidHistoXi = new TH1D(*aCut.fMinvPurityAidHistoXi);
  else fMinvPurityAidHistoXi = 0;

  if(aCut.fCTauXi) fCTauXi = new TH1D(*aCut.fCTauXi);
  else fCTauXi = 0;

  return *this;
}

AliFemtoXiTrackCutNSigmaFilter* AliFemtoXiTrackCutNSigmaFilter::Clone()
{
  return(new AliFemtoXiTrackCutNSigmaFilter(*this));
}

AliFemtoXiTrackCutNSigmaFilter::~AliFemtoXiTrackCutNSigmaFilter()
{
  delete fDaughterV0Filter;
  delete fDaughterBacPionFilter;

  delete fBacPionNSigmaFilter;
}


bool AliFemtoXiTrackCutNSigmaFilter::PassV0(const AliFemtoXi* aXi)
{
  if(!fIsDaughterV0FilterUpdated) UpdateDaughterV0Filter();
  return fDaughterV0Filter->Pass(aXi);
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
  
  fCTauXi->Fill(GetCTauXi(aXi));
  return true;
}


AliFemtoString AliFemtoXiTrackCutNSigmaFilter::Report()
{
  TString report;
  report += TString::Format("Usings custom BacPion NSigma Filter:\t%i\n", fUseCustomBacPionNSigmaFilter)
          + AliFemtoXiTrackCut::Report()
          + fDaughterV0Filter->Report();
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
  if(fBuildPurityAidV0) tOutputList->Add(fDaughterV0Filter->GetMinvPurityAidHistoV0());

  if(fDaughterV0Filter->GetCTauHist()!=NULL) tOutputList->Add(fDaughterV0Filter->GetCTauHist());
  tOutputList->Add(fCTauXi);

  if(fDaughterV0Filter->GetBuildMisIDHistograms()) tOutputList->Add(fDaughterV0Filter->GetMisIDHistos());

  return tOutputList;
}

void AliFemtoXiTrackCutNSigmaFilter::UpdateDaughterV0Filter()
{
  //TODO find more robust way of doing this, so it can evolve is AliFemtoV0Track cut 
  //  gains more attributes
  fDaughterV0Filter->SetInvariantMassLambda(fInvMassLambdaMin,fInvMassLambdaMax);
  fDaughterV0Filter->SetInvariantMassK0s(fInvMassK0sMin,fInvMassK0sMax);

  fDaughterV0Filter->SetMinDaughtersToPrimVertex(fMinDcaDaughterPosToVert,fMinDcaDaughterNegToVert);
  fDaughterV0Filter->SetMaxDcaV0Daughters(fMaxDcaV0Daughters);
  fDaughterV0Filter->SetMaxDcaV0(fMaxDcaV0);
  fDaughterV0Filter->SetMinDcaV0(fMinDcaV0);
  fDaughterV0Filter->SetMaxCosPointingAngle(fMaxCosPointingAngle);
  fDaughterV0Filter->SetMinCosPointingAngle(fMinCosPointingAngle);
  fDaughterV0Filter->SetMaxV0DecayLength(fMaxDecayLength);
  fDaughterV0Filter->SetMinTransverseDistancePrimSecVtx(fMinTransverseDistancePrimSecVtx);
  fDaughterV0Filter->SetParticleType(fParticleType);
  fDaughterV0Filter->SetEta(fEta);
  fDaughterV0Filter->SetPt(fPtMin,fPtMax);
  fDaughterV0Filter->SetEtaDaughters(fMaxEtaDaughters);
  fDaughterV0Filter->SetTPCnclsDaughters(fTPCNclsDaughters);
  fDaughterV0Filter->SetNdofDaughters(fNdofDaughters);
  fDaughterV0Filter->SetStatusDaughters(fStatusDaughters);
  fDaughterV0Filter->SetPtPosDaughter(fPtMinPosDaughter,fPtMaxPosDaughter);
  fDaughterV0Filter->SetPtNegDaughter(fPtMinNegDaughter,fPtMaxNegDaughter);
  fDaughterV0Filter->SetOnFlyStatus(fOnFlyStatus);
  fDaughterV0Filter->SetIgnoreOnFlyStatus(fIgnoreOnFlyStatus);
  fDaughterV0Filter->SetMinAvgSeparation(fMinAvgSepDaughters);

  fDaughterV0Filter->SetRadiusV0Min(fRadiusV0Min);
  fDaughterV0Filter->SetRadiusV0Max(fRadiusV0Max);

  fDaughterV0Filter->SetRequireTOFPion(fRequireTOFPion);
  fDaughterV0Filter->SetRequireTOFProton(fRequireTOFProton);
  
  fDaughterV0Filter->SetNsigmaPosDaughter(fNsigmaPosDaughterTPC);
  fDaughterV0Filter->SetNsigmaNegDaughter(fNsigmaNegDaughterTPC);
  fDaughterV0Filter->SetNsigmaPosDaughter(fNsigmaPosDaughterTPC,fNsigmaPosDaughterTOF);
  fDaughterV0Filter->SetNsigmaNegDaughter(fNsigmaNegDaughterTPC,fNsigmaNegDaughterTOF);


  fIsDaughterV0FilterUpdated = true;
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
  case kKaon:
  case kProton:
    fDaughterV0Filter->CreateCustomNSigmaFilter((AliFemtoV0TrackCutNSigmaFilter::DaughterParticleType)aDaughterType);
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
  if(aDaughterType==kPion || aDaughterType==kKaon || aDaughterType==kProton) fDaughterV0Filter->AddTPCAndTOFNSigmaCut((AliFemtoV0TrackCutNSigmaFilter::DaughterParticleType)aDaughterType, aMomMin, aMomMax, aNSigmaValueTPC, aNSigmaValueTOF);
  else if(aDaughterType==kBacPion) fBacPionNSigmaFilter->AddTPCAndTOFCut(aMomMin,aMomMax,aNSigmaValueTPC,aNSigmaValueTOF);
  else {/*cerr << "Invalid DaughterParticleType selection in AddTPCAndTOFNSigmaCut.  No cut will be added to the filter!!!!!" << endl;*/}
}

void AliFemtoXiTrackCutNSigmaFilter::AddTPCNSigmaCut(XiDaughterParticleType aDaughterType, double aMomMin, double aMomMax, double aNSigmaValueTPC)
{
  if(aDaughterType==kPion || aDaughterType==kKaon || aDaughterType==kProton) fDaughterV0Filter->AddTPCNSigmaCut((AliFemtoV0TrackCutNSigmaFilter::DaughterParticleType)aDaughterType, aMomMin, aMomMax, aNSigmaValueTPC);
  else if(aDaughterType==kBacPion) fBacPionNSigmaFilter->AddTPCCut(aMomMin,aMomMax,aNSigmaValueTPC);
  else {/*cerr << "Invalid DaughterParticleType selection in AddTPCNSigmaCut.  No cut will be added to the filter!!!!!" << endl;*/}
}

void AliFemtoXiTrackCutNSigmaFilter::AddTOFNSigmaCut(XiDaughterParticleType aDaughterType, double aMomMin, double aMomMax, double aNSigmaValueTOF)
{
  if(aDaughterType==kPion || aDaughterType==kKaon || aDaughterType==kProton) fDaughterV0Filter->AddTOFNSigmaCut((AliFemtoV0TrackCutNSigmaFilter::DaughterParticleType)aDaughterType, aMomMin, aMomMax, aNSigmaValueTOF);
  else if(aDaughterType==kBacPion) fBacPionNSigmaFilter->AddTOFCut(aMomMin,aMomMax,aNSigmaValueTOF);
  else {/*cerr << "Invalid DaughterParticleType selection in AddTOFNSigmaCut.  No cut will be added to the filter!!!!!" << endl;*/}
}



void AliFemtoXiTrackCutNSigmaFilter::CreateCustomV0Rejection(AliFemtoV0Type aV0TypeToReject)
{
  UpdateDaughterV0Filter();
  fIsDaughterV0FilterUpdated = false;  //In case any further changes are made before running,
                                       //but fDaughterV0Filter needs to know its type to know which MisID histos to build
  if(!fDaughterV0Filter->GetBuildMisIDHistograms()) fDaughterV0Filter->SetBuildMisIDHistograms(true);
  fDaughterV0Filter->CreateCustomV0Rejection(aV0TypeToReject);
}


void AliFemtoXiTrackCutNSigmaFilter::AddTPCAndTOFNSigmaCutToV0Rejection(AliFemtoV0Type aV0Type, int aDaughterCharge, double aMomMin, double aMomMax, double aNSigmaValueTPC, double aNSigmaValueTOF)
{
  fDaughterV0Filter->AddTPCAndTOFNSigmaCutToV0Rejection(aV0Type,aDaughterCharge,aMomMin,aMomMax,aNSigmaValueTPC,aNSigmaValueTOF);
}

void AliFemtoXiTrackCutNSigmaFilter::AddTPCAndTOFNSigmaCutToV0Rejection(AliFemtoV0Type aV0Type, double aMomMinPos, double aMomMaxPos, double aNSigmaValueTPCPos, double aNSigmaValueTOFPos,
                                                                                                double aMomMinNeg, double aMomMaxNeg, double aNSigmaValueTPCNeg, double aNSigmaValueTOFNeg)
{
  fDaughterV0Filter->AddTPCAndTOFNSigmaCutToV0Rejection(aV0Type, aMomMinPos, aMomMaxPos, aNSigmaValueTPCPos, aNSigmaValueTOFPos,
                                                                 aMomMinNeg, aMomMaxNeg, aNSigmaValueTPCNeg, aNSigmaValueTOFNeg);
}


void AliFemtoXiTrackCutNSigmaFilter::AddTPCNSigmaCutToV0Rejection(AliFemtoV0Type aV0Type, int aDaughterCharge, double aMomMin, double aMomMax, double aNSigmaValueTPC)
{
  fDaughterV0Filter->AddTPCNSigmaCutToV0Rejection(aV0Type,aDaughterCharge,aMomMin,aMomMax,aNSigmaValueTPC);
}

void AliFemtoXiTrackCutNSigmaFilter::AddTPCNSigmaCutToV0Rejection(AliFemtoV0Type aV0Type, double aMomMinPos, double aMomMaxPos, double aNSigmaValueTPCPos,
                                                                                          double aMomMinNeg, double aMomMaxNeg, double aNSigmaValueTPCNeg)
{
  fDaughterV0Filter->AddTPCNSigmaCutToV0Rejection(aV0Type, aMomMinPos, aMomMaxPos, aNSigmaValueTPCPos,
                                                           aMomMinNeg, aMomMaxNeg, aNSigmaValueTPCNeg);
}


void AliFemtoXiTrackCutNSigmaFilter::AddTOFNSigmaCutToV0Rejection(AliFemtoV0Type aV0Type, int aDaughterCharge, double aMomMin, double aMomMax, double aNSigmaValueTOF)
{
  fDaughterV0Filter->AddTOFNSigmaCutToV0Rejection(aV0Type, aDaughterCharge, aMomMin, aMomMax, aNSigmaValueTOF);
}

void AliFemtoXiTrackCutNSigmaFilter::AddTOFNSigmaCutToV0Rejection(AliFemtoV0Type aV0Type, double aMomMinPos, double aMomMaxPos, double aNSigmaValueTOFPos,
                                                                                          double aMomMinNeg, double aMomMaxNeg, double aNSigmaValueTOFNeg)
{
  fDaughterV0Filter->AddTOFNSigmaCutToV0Rejection(aV0Type, aMomMinPos, aMomMaxPos, aNSigmaValueTOFPos,
                                                           aMomMinNeg, aMomMaxNeg, aNSigmaValueTOFNeg);
}



double AliFemtoXiTrackCutNSigmaFilter::GetCTauXi(const AliFemtoXi* aXi)
{
  double tPtotXi = aXi->MomXi().Mag();
  //TODO should I use PDG mass values, or aXi->MassXi()?
  double tMassXi = 1.32171;
  double tCTauXi = aXi->DecayLengthXi()*tMassXi/tPtotXi;

  return tCTauXi;
}

void AliFemtoXiTrackCutNSigmaFilter::SetMinvPurityAidHistoV0(const char* name, const char* title, const int& nbins, const float& aInvMassMin, const float& aInvMassMax)
{
  fBuildPurityAidV0 = true;
  fDaughterV0Filter->SetMinvPurityAidHistoV0(name,title,nbins,aInvMassMin,aInvMassMax);
}

void AliFemtoXiTrackCutNSigmaFilter::SetCTauHistoV0(int aNbins, double aMin, double aMax)
{
  fDaughterV0Filter->SetCTauHistoV0(aNbins,aMin,aMax);
}



