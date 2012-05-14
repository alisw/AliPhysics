/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *               
 *                                                                        *               
 * Author: The ALICE Off-line Project.                                    *               
 * Contributors are mentioned in the code where appropriate.              *               
 *                                                                        *               
 * Permission to use, copy, modify and distribute this software and its   *               
 * documentation strictly for non-commercial purposes is hereby granted   *               
 * without fee, provided that the above copyright notice appears in all   *               
 * copies and that both the copyright notice and this permission notice   *               
 * appear in the supporting documentation. The authors make no claims     *               
 * about the suitability of this software for any purpose. It is          *               
 * provided "as is" without express or implied warranty.                  *               
 **************************************************************************/

// Class for heavy-flavour electron v2 with EMCal triggered events
// Author: Denise Godoy


#include "TChain.h"
#include "TTree.h"
#include "TH2F.h"
#include "TMath.h"
#include "TCanvas.h"
#include "THnSparse.h"
#include "TLorentzVector.h"
#include "TString.h"
#include "TFile.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"

#include "AliESDEvent.h"
#include "AliESDHandler.h"
#include "AliAODEvent.h"
#include "AliAODHandler.h"

#include "AliAnalysisTaskElecV2.h"
#include "TGeoGlobalMagField.h"
#include "AliLog.h"
#include "AliAnalysisTaskSE.h"
#include "TRefArray.h"
#include "TVector.h"
#include "AliESDInputHandler.h"
#include "AliESDpid.h"
#include "AliESDtrackCuts.h"
#include "AliPhysicsSelection.h"
#include "AliESDCaloCluster.h"
#include "AliAODCaloCluster.h"
#include "AliEMCALRecoUtils.h"
#include "AliEMCALGeometry.h"
#include "AliGeomManager.h"
#include "stdio.h"
#include "TGeoManager.h"
#include "iostream"
#include "fstream"

#include "AliEMCALTrack.h"
#include "AliMagF.h"

#include "AliKFParticle.h"
#include "AliKFVertex.h"

#include "AliPID.h"
#include "AliPIDResponse.h"
#include "AliHFEcontainer.h"
#include "AliHFEcuts.h"
#include "AliHFEpid.h"
#include "AliHFEpidBase.h"
#include "AliHFEpidQAmanager.h"
#include "AliHFEtools.h"
#include "AliCFContainer.h"
#include "AliCFManager.h"

#include "AliEventplane.h"
#include "AliCentrality.h"

ClassImp(AliAnalysisTaskElecV2)
//________________________________________________________________________
AliAnalysisTaskElecV2::AliAnalysisTaskElecV2(const char *name) 
  : AliAnalysisTaskSE(name)
  ,fESD(0)
  ,fOutputList(0)
  ,fTrackCuts(0)
  ,fCuts(0)
  ,fIdentifiedAsOutInz(kFALSE)
  ,fPassTheEventCut(kFALSE)
  ,fRejectKinkMother(kFALSE)
  ,fVz(0.0)
  ,fCFM(0)	
  ,fPID(0)
  ,fPIDqa(0)	       
  ,fOpeningAngleCut(0.1)
  ,fInvmassCut(0.01)	
  ,fNoEvents(0)
  ,fTrkpt(0)
  ,fTrkEovPBef(0)	 
  ,fTrkEovPAft(0)	
  ,fdEdxBef(0)	 
  ,fdEdxAft(0)	 
  ,fInvmassLS(0)		
  ,fInvmassULS(0)		
  ,fOpeningAngleLS(0)	
  ,fOpeningAngleULS(0)	
  ,fPhotoElecPt(0)
  ,fSemiInclElecPt(0)
  ,fTrackPtBefTrkCuts(0)	 
  ,fTrackPtAftTrkCuts(0)
  ,fTPCnsigma(0)
  ,fCent(0)
  ,fevPlaneV0A(0)
  ,fevPlaneV0C(0)
  ,fevPlaneTPC(0)
  ,fTPCsubEPres(0)
  ,fEPres(0)
  ,fCorr(0)
  ,feTPCV2(0)
  ,feV2(0)
  ,fphoteV2(0)
  ,fChargPartV2(0)
{
  //Named constructor
  
  fPID = new AliHFEpid("hfePid");
  fTrackCuts = new AliESDtrackCuts();
  
  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  // Output slot #0 id reserved by the base class for AOD
  // Output slot #1 writes into a TH1 container
  // DefineOutput(1, TH1I::Class());
  DefineOutput(1, TList::Class());
  //  DefineOutput(3, TTree::Class());
}

//________________________________________________________________________
AliAnalysisTaskElecV2::AliAnalysisTaskElecV2() 
  : AliAnalysisTaskSE("DefaultAnalysis_AliAnalysisElecHadCorrel")
  ,fESD(0)
  ,fOutputList(0)
  ,fTrackCuts(0)
  ,fCuts(0)
  ,fIdentifiedAsOutInz(kFALSE)
  ,fPassTheEventCut(kFALSE)
  ,fRejectKinkMother(kFALSE)
  ,fVz(0.0)
  ,fCFM(0)	
  ,fPID(0)       
  ,fPIDqa(0)	       
  ,fOpeningAngleCut(0.1)
  ,fInvmassCut(0.01)	
  ,fNoEvents(0)
  ,fTrkpt(0)
  ,fTrkEovPBef(0)	 
  ,fTrkEovPAft(0)	 
  ,fdEdxBef(0)	 
  ,fdEdxAft(0)	 
  ,fInvmassLS(0)		
  ,fInvmassULS(0)		
  ,fOpeningAngleLS(0)	
  ,fOpeningAngleULS(0)	
  ,fPhotoElecPt(0)
  ,fSemiInclElecPt(0)
  ,fTrackPtBefTrkCuts(0)	 
  ,fTrackPtAftTrkCuts(0)	 	  
  ,fTPCnsigma(0)
  ,fCent(0)
  ,fevPlaneV0A(0)
  ,fevPlaneV0C(0)
  ,fevPlaneTPC(0)
  ,fTPCsubEPres(0)
  ,fEPres(0)
  ,fCorr(0)
  ,feTPCV2(0)
  ,feV2(0)
  ,fphoteV2(0)
  ,fChargPartV2(0)
{
	//Default constructor
	fPID = new AliHFEpid("hfePid");

	fTrackCuts = new AliESDtrackCuts();
	
	// Constructor
	// Define input and output slots here
	// Input slot #0 works with a TChain
	DefineInput(0, TChain::Class());
	// Output slot #0 id reserved by the base class for AOD
	// Output slot #1 writes into a TH1 container
	// DefineOutput(1, TH1I::Class());
	DefineOutput(1, TList::Class());
	//DefineOutput(3, TTree::Class());
}
//_________________________________________

AliAnalysisTaskElecV2::~AliAnalysisTaskElecV2()
{
  //Destructor 
  
  delete fOutputList;
  delete fPID;
  delete fCFM;
  delete fPIDqa;
  delete fTrackCuts;
}
//_________________________________________

void AliAnalysisTaskElecV2::UserExec(Option_t*)
{
  //Main loop
  //Called for each event
  
  // create pointer to event
  fESD = dynamic_cast<AliESDEvent*>(InputEvent());
  if (!fESD) {
    printf("ERROR: fESD not available\n");
    return;
  }
  
  if(!fCuts){
    AliError("HFE cuts not available");
    return;
  }
  
  if(!fPID->IsInitialized()){ 
    // Initialize PID with the given run number
    AliWarning("PID not initialised, get from Run no");
    fPID->InitializePID(fESD->GetRunNumber());
  }
 
 
  Int_t fNOtrks =  fESD->GetNumberOfTracks();
  const AliESDVertex *pVtx = fESD->GetPrimaryVertex();
  
  Double_t pVtxZ = -999;
  pVtxZ = pVtx->GetZ();
  
  if(TMath::Abs(pVtxZ)>10) return;
  fNoEvents->Fill(0);
  
  if(fNOtrks<2) return;
  
  AliPIDResponse *pidResponse = fInputHandler->GetPIDResponse();
  if(!pidResponse){
    AliDebug(1, "Using default PID Response");
    pidResponse = AliHFEtools::GetDefaultPID(kFALSE, fInputEvent->IsA() == AliAODEvent::Class()); 
  }
  
  fPID->SetPIDResponse(pidResponse);
  
  fCFM->SetRecEventInfo(fESD);
  
  Float_t cent = -1.;
  AliCentrality *centrality = fESD->GetCentrality(); 
  cent = centrality->GetCentralityPercentile("V0M");
  fCent->Fill(cent);
  
  if(cent>90.) return;
	
  //Event planes
  
  Double_t evPlaneV0A = TVector2::Phi_0_2pi(fESD->GetEventplane()->GetEventplane("V0A",fESD,2));
  if(evPlaneV0A > TMath::Pi()) evPlaneV0A = evPlaneV0A - TMath::Pi();
  fevPlaneV0A->Fill(evPlaneV0A);
  
  Double_t evPlaneV0C = TVector2::Phi_0_2pi(fESD->GetEventplane()->GetEventplane("V0C",fESD,2));
  if(evPlaneV0C > TMath::Pi()) evPlaneV0C = evPlaneV0C - TMath::Pi();
  fevPlaneV0C->Fill(evPlaneV0C);
  
  AliEventplane* esdTPCep = fESD->GetEventplane();
  TVector2 *standardQ = 0x0;
  Double_t qx = -999., qy = -999.;
  
  standardQ = esdTPCep->GetQVector(); 
  if(standardQ)
  {
	  qx = standardQ->X();
	  qy = standardQ->Y();
  }
  TVector2 qVectorfortrack;
  qVectorfortrack.Set(qx,qy);
  Float_t evPlaneTPC = TVector2::Phi_0_2pi(qVectorfortrack.Phi())/2.;
  fevPlaneTPC->Fill(evPlaneTPC);
  
  TVector2 *qsub1a = esdTPCep->GetQsub1();
  TVector2 *qsub2a = esdTPCep->GetQsub2();
  Double_t evPlaneResTPC = -999.;
  if(qsub1a && qsub2a)
  {
	  evPlaneResTPC = TMath::Cos(2.*TVector2::Phi_0_2pi(qsub1a->Phi()/2.- qsub2a->Phi()/2.));
  }
    
  fTPCsubEPres->Fill(evPlaneResTPC,cent);
  
  Double_t evPlaneRes[4]={GetCos2DeltaPhi(evPlaneV0A,evPlaneV0C),GetCos2DeltaPhi(evPlaneV0A,evPlaneTPC),GetCos2DeltaPhi(evPlaneV0C,evPlaneTPC),cent};
  fEPres->Fill(evPlaneRes);
  
  // Track loop 
  for (Int_t iTracks = 0; iTracks < fESD->GetNumberOfTracks(); iTracks++) {
    AliESDtrack* track = fESD->GetTrack(iTracks);
    if (!track) {
      printf("ERROR: Could not receive track %d\n", iTracks);
      continue;
    }
    
    if(TMath::Abs(track->Eta())>0.7) continue;
    
    fTrackPtBefTrkCuts->Fill(track->Pt());		
    // RecKine: ITSTPC cuts  
    if(!ProcessCutStep(AliHFEcuts::kStepRecKineITSTPC, track)) continue;
    
    //RecKink
    if(fRejectKinkMother) { // Quick and dirty fix to reject both kink mothers and daughters
      if(track->GetKinkIndex(0) != 0) continue;
    } 
    
    // RecPrim
    if(!ProcessCutStep(AliHFEcuts::kStepRecPrim, track)) continue;
    
    // HFEcuts: ITS layers cuts
    if(!ProcessCutStep(AliHFEcuts::kStepHFEcutsITS, track)) continue;
    
    // HFE cuts: TPC PID cleanup
    if(!ProcessCutStep(AliHFEcuts::kStepHFEcutsTPC, track)) continue;
    
    fTrackPtAftTrkCuts->Fill(track->Pt());		
    
    Double_t clsE = -999., p = -999., EovP=-999., pt = -999., dEdx=-999., fTPCnSigma=0, phi=-999.;
    
    // Track extrapolation
    
    pt = track->Pt();
    if(pt<2) continue;
    fTrkpt->Fill(pt);
    
    Int_t clsId = track->GetEMCALcluster();
    if (clsId>0){
      AliESDCaloCluster *cluster = fESD->GetCaloCluster(clsId);
      if(cluster && cluster->IsEMCAL()){
	clsE = cluster->E();
      }
    }
        
    p = track->P();
    phi = track->Phi();
    dEdx = track->GetTPCsignal();
    EovP = clsE/p;
    fTPCnSigma = fPID->GetPIDResponse() ? fPID->GetPIDResponse()->NumberOfSigmasTPC(track, AliPID::kElectron) : 1000;
    fdEdxBef->Fill(p,dEdx);
    fTPCnsigma->Fill(p,fTPCnSigma);
       
    Double_t corr[8]={phi,fTPCnSigma,cent,pt,EovP,GetDeltaPhi(phi,evPlaneTPC),GetDeltaPhi(phi,evPlaneV0A),GetDeltaPhi(phi,evPlaneV0C)};
    fCorr->Fill(corr);
       
    if(fTPCnSigma >= 1.5 && fTPCnSigma <= 3)fTrkEovPBef->Fill(pt,EovP);
    Int_t pidpassed = 1;
    
    //--- track accepted
    AliHFEpidObject hfetrack;
    hfetrack.SetAnalysisType(AliHFEpidObject::kESDanalysis);
    hfetrack.SetRecTrack(track);
    hfetrack.SetPbPb();
    if(!fPID->IsSelected(&hfetrack, NULL, "", fPIDqa)) pidpassed = 0;
    
    Double_t corrV2[7]={phi,cent,pt,EovP,GetCos2DeltaPhi(phi,evPlaneTPC),GetCos2DeltaPhi(phi,evPlaneV0A),GetCos2DeltaPhi(phi,evPlaneV0C)};
    fChargPartV2->Fill(corrV2); 

    if(fTPCnSigma >= -0.5){
      	Double_t qX = standardQ->X() - esdTPCep->GetQContributionX(track); 
	Double_t qY = standardQ->Y() - esdTPCep->GetQContributionY(track); 
	TVector2 newQVectorfortrack;
	newQVectorfortrack.Set(qX,qY);
	Double_t corrV2TPC = -999.;
	corrV2TPC = TVector2::Phi_0_2pi(newQVectorfortrack.Phi())/2; 

	Double_t correctedV2[5]={cent,pt,GetCos2DeltaPhi(phi,corrV2TPC),GetCos2DeltaPhi(phi,evPlaneV0A),GetCos2DeltaPhi(phi,evPlaneV0C)};

	feTPCV2->Fill(correctedV2);
    }
    
    if(pidpassed==0) continue;
    
    Double_t qX = standardQ->X() - esdTPCep->GetQContributionX(track); 
    Double_t qY = standardQ->Y() - esdTPCep->GetQContributionY(track); 
    TVector2 newQVectorfortrack;
    newQVectorfortrack.Set(qX,qY);
    Double_t corrV2TPC = -999.;
    corrV2TPC = TVector2::Phi_0_2pi(newQVectorfortrack.Phi())/2; 
        
    Double_t correctedV2[5]={cent,pt,GetCos2DeltaPhi(phi,corrV2TPC),GetCos2DeltaPhi(phi,evPlaneV0A),GetCos2DeltaPhi(phi,evPlaneV0C)};
    
    feV2->Fill(correctedV2);
    
    fTrkEovPAft->Fill(pt,EovP);
    fdEdxAft->Fill(p,dEdx);
    
    Bool_t fFlagPhotonicElec = kFALSE;
    SelectPhotonicElectron(iTracks,track,fFlagPhotonicElec);
    
    if(fFlagPhotonicElec){
      fphoteV2->Fill(correctedV2);
      fPhotoElecPt->Fill(pt);
    }
    
    if(!fFlagPhotonicElec) fSemiInclElecPt->Fill(pt);
 }
 PostData(1, fOutputList);
}
//_________________________________________
void AliAnalysisTaskElecV2::UserCreateOutputObjects()
{
  //--------Initialize PID
  fPID->SetHasMCData(kFALSE);
  if(!fPID->GetNumberOfPIDdetectors()) 
    {
      fPID->AddDetector("TPC", 0);
      fPID->AddDetector("EMCAL", 1);
    }
  
  fPID->SortDetectors(); 
  fPIDqa = new AliHFEpidQAmanager();
  fPIDqa->Initialize(fPID);
  
  //--------Initialize correction Framework and Cuts
  fCFM = new AliCFManager;
  const Int_t kNcutSteps = AliHFEcuts::kNcutStepsMCTrack + AliHFEcuts::kNcutStepsRecTrack + AliHFEcuts::kNcutStepsDETrack;
  fCFM->SetNStepParticle(kNcutSteps);
  for(Int_t istep = 0; istep < kNcutSteps; istep++)
    fCFM->SetParticleCutsList(istep, NULL);
  
  if(!fCuts){
    AliWarning("Cuts not available. Default cuts will be used");
    fCuts = new AliHFEcuts;
    fCuts->CreateStandardCuts();
  }
  fCuts->Initialize(fCFM);
  
  //---------Output Tlist
  fOutputList = new TList();
  fOutputList->SetOwner();
  fOutputList->Add(fPIDqa->MakeList("PIDQA"));
  
  fNoEvents = new TH1F("fNoEvents","",1,0,1) ;
  fOutputList->Add(fNoEvents);
  
  fTrkpt = new TH1F("fTrkpt","track pt",100,0,50);
  fOutputList->Add(fTrkpt);
  
  fTrackPtBefTrkCuts = new TH1F("fTrackPtBefTrkCuts","track pt before track cuts",100,0,50);
  fOutputList->Add(fTrackPtBefTrkCuts);
  
  fTrackPtAftTrkCuts = new TH1F("fTrackPtAftTrkCuts","track pt after track cuts",100,0,50);
  fOutputList->Add(fTrackPtAftTrkCuts);
  
  fTPCnsigma = new TH2F("fTPCnsigma", "TPC - n sigma",100,0,50,200,-10,10);
  fOutputList->Add(fTPCnsigma);
  
  fTrkEovPBef = new TH2F("fTrkEovPBef","track E/p before HFE pid",100,0,50,100,0,2);
  fOutputList->Add(fTrkEovPBef);
  
  fTrkEovPAft = new TH2F("fTrkEovPAft","track E/p after HFE pid",100,0,50,100,0,2);
  fOutputList->Add(fTrkEovPAft);
  
  fdEdxBef = new TH2F("fdEdxBef","track dEdx vs p before HFE pid",100,0,50,150,0,150);
  fOutputList->Add(fdEdxBef);
  
  fdEdxAft = new TH2F("fdEdxAft","track dEdx vs p after HFE pid",100,0,50,150,0,150);
  fOutputList->Add(fdEdxAft);
  
  fInvmassLS = new TH1F("fInvmassLS", "Inv mass of LS (e,e); mass(GeV/c^2); counts;", 500,0,0.5);
  fOutputList->Add(fInvmassLS);
  
  fInvmassULS = new TH1F("fInvmassULS", "Inv mass of ULS (e,e); mass(GeV/c^2); counts;", 500,0,0.5);
  fOutputList->Add(fInvmassULS);
  
  fOpeningAngleLS = new TH1F("fOpeningAngleLS","Opening angle for LS pairs",100,0,1);
  fOutputList->Add(fOpeningAngleLS);
  
  fOpeningAngleULS = new TH1F("fOpeningAngleULS","Opening angle for ULS pairs",100,0,1);
  fOutputList->Add(fOpeningAngleULS);
  
  fPhotoElecPt = new TH1F("fPhotoElecPt", "photonic electron pt",100,0,50);
  fOutputList->Add(fPhotoElecPt);
  
  fSemiInclElecPt = new TH1F("fSemiInclElecPt", "Semi-inclusive electron pt",100,0,50);
  fOutputList->Add(fSemiInclElecPt);
  
  fCent = new TH1F("fCent","Centrality",100,0,100) ;
  fOutputList->Add(fCent);
  
  fevPlaneV0A = new TH1F("fevPlaneV0A","V0A EP",100,0,TMath::Pi()) ;
  fOutputList->Add(fevPlaneV0A);
  
  fevPlaneV0C = new TH1F("fevPlaneV0C","V0C EP",100,0,TMath::Pi()) ;
  fOutputList->Add(fevPlaneV0C);
  
  fevPlaneTPC = new TH1F("fevPlaneTPC","TPC EP",100,0,TMath::Pi()) ;
  fOutputList->Add(fevPlaneTPC);
    
  fTPCsubEPres = new TH2F("fTPCsubEPres","TPC subevent plane resolution",100,-1,1,90,0,90);
  fOutputList->Add(fTPCsubEPres);
  
  Int_t binsv1[4]={100,100,100,90}; // V0A-V0C, V0A-TPC, V0C-TPC, cent
  Double_t xminv1[4]={-1,-1,-1,0};
  Double_t xmaxv1[4]={1,1,1,90}; 
  fEPres = new THnSparseD ("fEPres","EP resolution",4,binsv1,xminv1,xmaxv1);
  fOutputList->Add(fEPres);
	
  Int_t binsv2[8]={100,100,90,100,100,100,100,100}; // phi,fTPCnSigma,cent, pt, EovP, TPCdeltaPhi, V0AdeltaPhi, V0CdeltaPhi
  Double_t xminv2[8]={0,-3.5,0,0,0,0,0,0};
  Double_t xmaxv2[8]={2*TMath::Pi(),3.5,90,50,3,TMath::Pi(),TMath::Pi(),TMath::Pi()}; 
  fCorr = new THnSparseD ("fCorr","Correlations",8,binsv2,xminv2,xmaxv2);
  fOutputList->Add(fCorr);
    
  Int_t binsv3[5]={90,100,100,100,100}; // cent, pt, TPCcos2DeltaPhi, V0Acos2DeltaPhi, V0Ccos2DeltaPhi
  Double_t xminv3[5]={0,0,-1,-1,-1};
  Double_t xmaxv3[5]={90,50,1,1,1}; 
  feV2 = new THnSparseD ("feV2","inclusive electron v2",5,binsv3,xminv3,xmaxv3);
  fOutputList->Add(feV2);
  
  Int_t binsv4[5]={90,100,100,100,100}; // cent, pt, TPCdeltaPhi, V0AdeltaPhi, V0CdeltaPhi
  Double_t xminv4[5]={0,0,-1,-1,-1};
  Double_t xmaxv4[5]={90,50,1,1,1}; 
  fphoteV2 = new THnSparseD ("fphoteV2","photonic electron v2",5,binsv4,xminv4,xmaxv4);
  fOutputList->Add(fphoteV2);
  
  Int_t binsv5[7]={100,90,100,100,100,100,100}; // phi, cent, pt, EovP, TPCdeltaPhi, V0AdeltaPhi, V0CdeltaPhi
  Double_t xminv5[7]={0,0,0,0,-1,-1,-1};
  Double_t xmaxv5[7]={2*TMath::Pi(),90,50,3,1,1,1}; 
  fChargPartV2 = new THnSparseD ("fChargPartV2","Charged particle v2",7,binsv5,xminv5,xmaxv5);
  fOutputList->Add(fChargPartV2);
  
  Int_t binsv6[5]={90,100,100,100,100}; // cent, pt, TPCdeltaPhi, V0AdeltaPhi, V0CdeltaPhi
  Double_t xminv6[5]={0,0,-1,-1,-1};
  Double_t xmaxv6[5]={90,50,1,1,1}; 
  feTPCV2 = new THnSparseD ("feTPCV2","inclusive electron v2 (TPC)",5,binsv6,xminv6,xmaxv6);
  fOutputList->Add(feTPCV2);
  
  PostData(1,fOutputList);
}

//________________________________________________________________________
void AliAnalysisTaskElecV2::Terminate(Option_t *)
{
  // Info("Terminate");
	AliAnalysisTaskSE::Terminate();
}

//________________________________________________________________________
Bool_t AliAnalysisTaskElecV2::ProcessCutStep(Int_t cutStep, AliVParticle *track)
{
  // Check single track cuts for a given cut step
  const Int_t kMCOffset = AliHFEcuts::kNcutStepsMCTrack;
  if(!fCFM->CheckParticleCuts(cutStep + kMCOffset, track)) return kFALSE;
  return kTRUE;
}
//_________________________________________
void AliAnalysisTaskElecV2::SelectPhotonicElectron(Int_t itrack, AliESDtrack *track, Bool_t &fFlagPhotonicElec)
{
  //Identify non-heavy flavour electrons using Invariant mass method
  
  fTrackCuts->SetAcceptKinkDaughters(kFALSE);
  fTrackCuts->SetRequireTPCRefit(kTRUE);
  fTrackCuts->SetEtaRange(-0.7,0.7);
  fTrackCuts->SetRequireSigmaToVertex(kTRUE);
  fTrackCuts->SetMaxChi2PerClusterTPC(3.5);
  fTrackCuts->SetMinNClustersTPC(100);
  
  const AliESDVertex *pVtx = fESD->GetPrimaryVertex();
  
  Bool_t flagPhotonicElec = kFALSE;
  
  for(Int_t jTracks = itrack+1; jTracks<fESD->GetNumberOfTracks(); jTracks++){
    AliESDtrack* trackAsso = fESD->GetTrack(jTracks);
    if (!trackAsso) {
      printf("ERROR: Could not receive track %d\n", jTracks);
      continue;
    }
    
    Double_t dEdxAsso = -999., ptAsso=-999., openingAngle = -999.;
    Double_t mass=999., width = -999;
    Bool_t fFlagLS=kFALSE, fFlagULS=kFALSE;
    
    dEdxAsso = trackAsso->GetTPCsignal();
    ptAsso = trackAsso->Pt();
    Int_t chargeAsso = trackAsso->Charge();
    Int_t charge = track->Charge();
    
    if(ptAsso <0.3) continue;
    if(!fTrackCuts->AcceptTrack(trackAsso)) continue;
    if(dEdxAsso <70 || dEdxAsso>100) continue; //11a pass1
    
    Int_t fPDGe1 = 11; Int_t fPDGe2 = 11;
    if(charge>0) fPDGe1 = -11;
    if(chargeAsso>0) fPDGe2 = -11;
    
    if(charge == chargeAsso) fFlagLS = kTRUE;
    if(charge != chargeAsso) fFlagULS = kTRUE;
    
    AliKFParticle ge1(*track, fPDGe1);
    AliKFParticle ge2(*trackAsso, fPDGe2);
    AliKFParticle recg(ge1, ge2);
    
    if(recg.GetNDF()<1) continue;
    Double_t chi2recg = recg.GetChi2()/recg.GetNDF();
    if(TMath::Sqrt(TMath::Abs(chi2recg))>3.) continue;
    
    AliKFVertex primV(*pVtx);
    primV += recg;
    recg.SetProductionVertex(primV);
    
    recg.SetMassConstraint(0,0.0001);
    
    openingAngle = ge1.GetAngle(ge2);
    if(fFlagLS) fOpeningAngleLS->Fill(openingAngle);
    if(fFlagULS) fOpeningAngleULS->Fill(openingAngle);
    
    if(openingAngle > fOpeningAngleCut) continue;
    
    recg.GetMass(mass,width);
    
    if(fFlagLS) fInvmassLS->Fill(mass);
    if(fFlagULS) fInvmassULS->Fill(mass);
          
    if(mass<fInvmassCut && fFlagULS && !flagPhotonicElec){
      flagPhotonicElec = kTRUE;
    }
    
  }
  fFlagPhotonicElec = flagPhotonicElec;
  
}
//_________________________________________
Double_t AliAnalysisTaskElecV2::GetCos2DeltaPhi(Double_t phiA,Double_t phiB) const
{
  //Get cos[2(phi-psi_EP)] or cos[2(psi_subEP1 - psi_subEP2)]
  Double_t dPhi = TVector2::Phi_0_2pi(phiA - phiB); 
  if(dPhi > TMath::Pi()) dPhi = dPhi - TMath::Pi();
  Double_t cos2DeltaPhi = TMath::Cos(2*dPhi);
  
  return cos2DeltaPhi;
}

//_________________________________________
Double_t AliAnalysisTaskElecV2::GetDeltaPhi(Double_t phiA,Double_t phiB) const
{
  //Get phi-psi_EP
  Double_t dPhi = TVector2::Phi_0_2pi(phiA - phiB); 
  if(dPhi > TMath::Pi()) dPhi = dPhi - TMath::Pi();
  
  return dPhi;
}
