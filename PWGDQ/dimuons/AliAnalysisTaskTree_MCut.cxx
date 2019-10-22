/**************************************************************************
 * Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
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

/* $Id: AliAnalysisTaskTree_MCut.cxx $ */

//-----------------------------------------------------------------------------
// Analysis task to compute muon/dimuon kinematic distributions
// The output is a list of histograms.
// The macro class can run on AOD or in the train with the ESD filter.
// R. Arnaldi
//
//-----------------------------------------------------------------------------

//#ifndef AliAnalysisTaskTree_MCut_CXX
//#define AliAnalysisTaskTree_MCut_CXX

// ROOT includes
#include "TROOT.h"
#include "TH1.h"
#include "TH2.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TMath.h"
#include "TObjString.h"
#include "TObjArray.h"
#include "TF1.h"
#include "TList.h"
#include "TStyle.h"
#include "TLorentzVector.h"
#include "TString.h"
#include "TTree.h"
#include "TFile.h"

#include "AliInputEventHandler.h"
#include "AliAODHeader.h"
#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskSE.h"
#include "AliMuonTrackCuts.h"   

#include "AliAnalysisTaskTree_MCut.h"

Double_t CostHE(AliAODTrack*, AliAODTrack*);
Double_t CostCS(AliAODTrack*, AliAODTrack*);
Double_t PhiHE(AliAODTrack*, AliAODTrack*);
Double_t PhiCS(AliAODTrack*, AliAODTrack*);

ClassImp(AliAnalysisTaskTree_MCut)
//__________________________________________________________________________
AliAnalysisTaskTree_MCut::AliAnalysisTaskTree_MCut() :
  AliAnalysisTaskSE(),
  fOutputTree(0x0),
  fNevt(0x0),
  fBeamEnergy(0.),
  fkAnalysisType(0x0),
  fPeriod(0x0),
  fCountTotEv(0x0),
  fCountTrigger(0x0),
  fCountCINT7(0x0),
  fCountCMUL7(0x0),
  fCountCMLL7(0x0),
  fCountCMSL7(0x0),
  fCountCMSH7(0x0),
  fNMuons(0x0),
  fNTracklets(0x0),
  fNContributors(0x0),
  fNDimu(0x0),
  fIsPhysSelected(0x0),
  fAODEvent(0x0),
//  fTrigClass(0x0),
  finpmask(0)
{
  //
  //Default ctor
  //    
  fMuonTrackCuts = new AliMuonTrackCuts("StandardMuonTrackCuts", "StandardMuonTrackCuts");
  fMuonTrackCuts->SetFilterMask(AliMuonTrackCuts::kMuPdca);  
  fMuonTrackCuts->SetAllowDefaultParams(kTRUE);

  fVertex[0]=999.; fVertex[1]=999.; fVertex[2]=999.;
  for(Int_t i=0; i<1500;i++){
    fPt[i]=999.;
    fE[i]=999.;
    fPx[i]=999; 
    fPy[i]=999; 
    fPz[i]=999; 
    fY[i]=999.; 
    fEta[i]=999.; 
    fMatchTrig[i]=999.; 
    fTrackChi2[i]=999.; 
    fMatchTrigChi2[i]=999.;
    fDCA[i]=999.;
    fCharge[i]=999;
    fRAtAbsEnd[i]=999;
    fpDCA[i] = 999;   
  }
  for(Int_t i=0; i<400;i++){  
    fDimuPt[i]=999.; 
    fDimuPx[i]=999.; 
    fDimuPy[i]=999.; 
    fDimuPz[i]=999.; 
    fDimuY[i]=999.; 
    fDimuMass[i]=999.;
    fDimuCharge[i]=999;
    fDimuMatch[i]=999;
    fDimuCostHE[i] = 999;
    fDimuPhiHE[i] = 999;
    fDimuCostCS[i] = 999;
    fDimuPhiCS[i] = 999;
    for(Int_t k=0;k<2;k++) fDimuMu[i][k]=999;
  }   
}

//__________________________________________________________________________
AliAnalysisTaskTree_MCut::AliAnalysisTaskTree_MCut(const char *name) :
  AliAnalysisTaskSE(name),
  fOutputTree(0x0),
  fNevt(0x0),
  fBeamEnergy(0.),
  fkAnalysisType(0x0),
  fPeriod(0x0),
  fCountTotEv(0x0),
  fCountTrigger(0x0),
  fCountCINT7(0x0),
  fCountCMUL7(0x0),
  fCountCMLL7(0x0),
  fCountCMSL7(0x0),
  fCountCMSH7(0x0),  
  fNMuons(0x0),
  fNTracklets(0x0),
  fNContributors(0x0),
  fNDimu(0x0),
  fIsPhysSelected(0x0),
  fAODEvent(0x0),
//  fTrigClass(0x0),
  finpmask(0)
{
  //
  // Constructor. Initialization of Inputs and Outputs
  //
  Info("AliAnalysisTaskTree_MCut","Calling Constructor");
  
  fMuonTrackCuts = new AliMuonTrackCuts("StandardMuonTrackCuts", "TestStandardMuonTrackCuts");
  fMuonTrackCuts->SetFilterMask(AliMuonTrackCuts::kMuPdca);  
  fMuonTrackCuts->SetAllowDefaultParams(kTRUE);

  fVertex[0]=999.; fVertex[1]=999.; fVertex[2]=999.;
  for(Int_t i=0; i<1500;i++){
    fPt[i]=999.;
    fE[i]=999.;
    fPx[i]=999; 
    fPy[i]=999; 
    fPz[i]=999; 
    fY[i]=999.; 
    fEta[i]=999.; 
    fMatchTrig[i]=999.; 
    fTrackChi2[i]=999.; 
    fMatchTrigChi2[i]=999.;
    fDCA[i]=999.;
    fCharge[i]=999;
    fRAtAbsEnd[i]=999;
    fpDCA[i] = 999;   
  }
  for(Int_t i=0; i<400;i++){  
    fDimuPt[i]=999.; 
    fDimuPx[i]=999.; 
    fDimuPy[i]=999.; 
    fDimuPz[i]=999.; 
    fDimuY[i]=999.; 
    fDimuMass[i]=999.;
    fDimuCharge[i]=999;
    fDimuMatch[i]=999;
    fDimuCostHE[i] = 999;
    fDimuPhiHE[i] = 999;
    fDimuCostCS[i] = 999;
    fDimuPhiCS[i] = 999;
    for(Int_t k=0;k<2;k++) fDimuMu[i][k]=999;
  } 
  
  DefineOutput(1,TTree::Class());
}

//___________________________________________________________________________
AliAnalysisTaskTree_MCut& AliAnalysisTaskTree_MCut::operator=(const AliAnalysisTaskTree_MCut& c) 
{
  //
  // Assignment operator
  //
  if (this!=&c) {
    AliAnalysisTaskSE::operator=(c) ;
  }
  return *this;
}

//___________________________________________________________________________
AliAnalysisTaskTree_MCut::AliAnalysisTaskTree_MCut(const AliAnalysisTaskTree_MCut& c) :
  AliAnalysisTaskSE(c),  
  fOutputTree(c.fOutputTree),
  fNevt(c.fNevt),
  fBeamEnergy(c.fBeamEnergy),
  fkAnalysisType(c.fkAnalysisType),
  fPeriod(c.fPeriod),
  fCountTotEv(c.fCountTotEv),
  fCountTrigger(c.fCountTrigger),
  fCountCINT7(c.fCountCINT7),
  fCountCMUL7(c.fCountCMUL7),
  fCountCMLL7(c.fCountCMLL7),
  fCountCMSL7(c.fCountCMSL7),
  fCountCMSH7(c.fCountCMSH7),  
  fNMuons(c.fNMuons),
  fNTracklets(c.fNTracklets),
  fNContributors(c.fNContributors),
  fNDimu(c.fNDimu),
  fIsPhysSelected(c.fIsPhysSelected),
  fAODEvent(c.fAODEvent),
//  fTrigClass(c.fTrigClass),
  finpmask(c.finpmask),
  fMuonTrackCuts(c.fMuonTrackCuts) 

 {
  //
  // Copy Constructor									
  //
}

//___________________________________________________________________________
AliAnalysisTaskTree_MCut::~AliAnalysisTaskTree_MCut() {
  //
  //destructor
  //
  Info("~AliAnalysisTaskTree_MCut","Calling Destructor");
  if (AliAnalysisManager::GetAnalysisManager()->GetAnalysisType() != AliAnalysisManager::kProofAnalysis) delete fOutputTree;
}

//___________________________________________________________________________
void AliAnalysisTaskTree_MCut::NotifyRun()
{
  fMuonTrackCuts->SetAllowDefaultParams(kTRUE);
  fMuonTrackCuts->SetRun((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()));
}

//___________________________________________________________________________
void AliAnalysisTaskTree_MCut::UserCreateOutputObjects(){

  if (fOutputTree) return; 
   
  OpenFile(1,"RECREATE");
  fOutputTree = new TTree("ppTree","Data Tree");

  fOutputTree->Branch("FiredTriggerClasses",fTrigClass,"FiredTriggerClasses/C");
  fOutputTree->Branch("inpmask",&finpmask,"inpmask/i"); 

  fOutputTree->Branch("NMuons",&fNMuons,"NMuons/I");
  fOutputTree->Branch("Vertex",fVertex,"Vertex[3]/D");

  fOutputTree->Branch("Pt",fPt,"Pt[NMuons]/D");
  fOutputTree->Branch("E",fE,"E[NMuons]/D");
  fOutputTree->Branch("Px",fPx,"Px[NMuons]/D");
  fOutputTree->Branch("Py",fPy,"Py[NMuons]/D");
  fOutputTree->Branch("Pz",fPz,"Pz[NMuons]/D");
  fOutputTree->Branch("Y",fY,"Y[NMuons]/D");
  fOutputTree->Branch("Eta",fEta,"Eta[NMuons]/D");
  fOutputTree->Branch("MatchTrig",fMatchTrig,"MatchTrig[NMuons]/I");
  fOutputTree->Branch("TrackChi2",fTrackChi2,"TrackChi2[NMuons]/D");
  fOutputTree->Branch("MatchTrigChi2",fMatchTrigChi2,"MatchTrigChi2[NMuons]/D");
  fOutputTree->Branch("Charge",fCharge,"Charge[NMuons]/I");
  fOutputTree->Branch("RAtAbsEnd",fRAtAbsEnd,"RAtAbsEnd[NMuons]/D");
  fOutputTree->Branch("pDCA",fpDCA,"pDCA[NMuons]/I");  

  fOutputTree->Branch("NDimu",&fNDimu,"NDimu/I");
  fOutputTree->Branch("DimuMu",fDimuMu,"DimuMu[NDimu][2]/I");
  fOutputTree->Branch("DimuPt",fDimuPt,"DimuPt[NDimu]/D");
  fOutputTree->Branch("DimuPx",fDimuPx,"DimuPx[NDimu]/D");
  fOutputTree->Branch("DimuPy",fDimuPy,"DimuPy[NDimu]/D");
  fOutputTree->Branch("DimuPz",fDimuPz,"DimuPz[NDimu]/D");
  fOutputTree->Branch("DimuY",fDimuY,"DimuY[NDimu]/D");
  fOutputTree->Branch("DimuMass",fDimuMass,"DimuMass[NDimu]/D");
  fOutputTree->Branch("DimuCharge",fDimuCharge,"DimuCharge[NDimu]/I");
  fOutputTree->Branch("DimuMatch",fDimuMatch,"DimuMatch[NDimu]/I");
  fOutputTree->Branch("DimuCostHE",fDimuCostHE,"DimuCostHE[NDimu]/D");
  fOutputTree->Branch("DimuPhiHE",fDimuPhiHE,"DimuPhiHE[NDimu]/D");
  fOutputTree->Branch("DimuCostCS",fDimuCostCS,"DimuCostCS[NDimu]/D");
  fOutputTree->Branch("DimuPhiCS",fDimuPhiCS,"DimuPhiCS[NDimu]/D");
  fOutputTree->Branch("IsPhysSelected",&fIsPhysSelected,"IsPhysSelected/O");

  fOutputTree->ls(); 

 PostData(1,fOutputTree); 
 
} 

//_________________________________________________
void AliAnalysisTaskTree_MCut::UserExec(Option_t *)
{

  fNMuons=0; 
  fNTracklets=-1;
  fNContributors=-1;
  fNDimu=0;
  fVertex[0]=999.; fVertex[1]=999.; fVertex[2]=999.;
  for(Int_t i=0; i<1500;i++){
    fPt[i]=999.;
    fE[i]=999.;
    fPx[i]=999; 
    fPy[i]=999; 
    fPz[i]=999; 
    fY[i]=999.; 
    fEta[i]=999.; 
    fMatchTrig[i]=999.; 
    fTrackChi2[i]=999.; 
    fMatchTrigChi2[i]=999.;
    fDCA[i]=999.;
    fCharge[i]=999;
    fRAtAbsEnd[i]=999;
    fpDCA[i] = 999.;  
  }
  for(Int_t i=0; i<400;i++){  
    fDimuPt[i]=999.; 
    fDimuPx[i]=999.; 
    fDimuPy[i]=999.; 
    fDimuPz[i]=999.; 
    fDimuY[i]=999.; 
    fDimuMass[i]=999.;
    fDimuCharge[i]=999.;
    fDimuMatch[i]=0;
    fDimuCostHE[i] = 999.;
    fDimuPhiHE[i] = 999.;
    fDimuCostCS[i] = 999.;
    fDimuPhiCS[i] = 999.;
    for(Int_t k=0;k<2;k++) fDimuMu[i][k]=999;
  } 
 
//
// Execute analysis for current event
//
  fAODEvent = dynamic_cast<AliAODEvent*> (InputEvent());
  if ( ! fAODEvent ) {
    AliError ("AOD event not found. Nothing done!");
    return;
  }
  
   AliAODHeader *aodheader=dynamic_cast<AliAODHeader*>(fAODEvent->GetHeader());
    TString firedtrigger = aodheader->GetFiredTriggerClasses();
    sprintf(fTrigClass,"%s",firedtrigger.Data());
    
   finpmask = aodheader->GetL0TriggerInputs(); 
    
  //   to apply physics selection
    UInt_t fSelectMask = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
    fIsPhysSelected = fSelectMask & (AliVEvent::kMuonUnlikeLowPt7 | AliVEvent::kMuonLikeLowPt7 | AliVEvent::kMuonSingleLowPt7 |   AliVEvent::kMuonSingleHighPt7  | AliVEvent::kINT7inMUON  | AliVEvent::kINT7); 	      

  AliAODVertex *PrimVertex =  fAODEvent->GetPrimaryVertex();
  fVertex[0]=PrimVertex->GetX();
  fVertex[1]=PrimVertex->GetY();
  fVertex[2]=PrimVertex->GetZ();
   
  //------------------------------------------------------
  // build dimuon object starting from single muons
  // Keep only dimuons with 
  // -- mass > mcut
  // -- 2.5<y<4
  // -- matchtrigger = 2
  //------------------------------------------------------
  
   Double_t MCut = 4;
  
   Int_t numdimu = 0;
   Int_t nummu = 0;  
   Int_t ntracks = fAODEvent->GetNumberOfTracks(); 
   if(ntracks!=0) {

   Int_t LabelOld1[1500];
   Int_t LabelOld2[1500];   
   Bool_t GoodMuon[1500]={kFALSE};
	      
   for (Int_t i=0;i<ntracks;i++){
	
     AliAODTrack *mu0=(AliAODTrack*)fAODEvent->GetTrack(i);
     if (!mu0->IsMuonTrack()) continue; 

     for(Int_t j=i+1;j<ntracks;j++){	     
	  
       AliAODTrack *mu1=(AliAODTrack*) fAODEvent->GetTrack(j);
       if (!mu1->IsMuonTrack()) continue;
	  
       AliAODDimuon *dimu=new AliAODDimuon(mu0,mu1);
	      
       if(dimu->Mass()>MCut && (mu0->GetMatchTrigger()>1 && mu1->GetMatchTrigger()>1) && dimu->Y()>-4 && dimu->Y()<-2.5){ 

	 fDimuMass[numdimu] = dimu->Mass();
	 fDimuPt[numdimu] = dimu->Pt();
	 fDimuPx[numdimu] = dimu->Px();
	 fDimuPy[numdimu] = dimu->Py();
	 fDimuPz[numdimu] = dimu->Pz();
	 fDimuY[numdimu] = dimu->Y();
	 fDimuCharge[numdimu]= dimu->Charge();
	 if(mu0->GetMatchTrigger()>1 && mu1->GetMatchTrigger()>1) fDimuMatch[numdimu]=2; 
	 else fDimuMatch[numdimu]=-1; 

	 fDimuCostHE[numdimu] = CostHE(mu0,mu1);
	 fDimuPhiHE[numdimu] = PhiHE(mu0,mu1);
	 fDimuCostCS[numdimu] = CostCS(mu0,mu1);
	 fDimuPhiCS[numdimu] = PhiCS(mu0,mu1);

	 LabelOld1[numdimu]=i;
	 LabelOld2[numdimu]=j;
	 GoodMuon[i]=kTRUE;
	 GoodMuon[j]=kTRUE;
	 numdimu++;      
	} // close loop on dimuon cuts
	delete dimu;
      } // close loop on second dimuon tracks
  } // close loop on first dimuon track    
 
 // loop on single muons to keep only muons belonging to a dimuon surviving cuts
  for(Int_t i=0;i<ntracks;i++){
     if(GoodMuon[i]) {
	AliAODTrack *mu0=(AliAODTrack*)fAODEvent->GetTrack(i);
	fCharge[nummu] = mu0->Charge();
	fPt[nummu] = mu0->Pt();
	fPx[nummu] = mu0->Px();
	fPy[nummu] = mu0->Py();
	fPz[nummu] = mu0->Pz();
	fY[nummu]  = mu0->Y();
	fEta[nummu]= mu0->Eta();
	fE[nummu] = mu0->E();
	fMatchTrig[nummu]   = mu0->GetMatchTrigger();
	fMatchTrigChi2[nummu]= mu0->GetChi2MatchTrigger();
	fRAtAbsEnd[nummu]=mu0->GetRAtAbsorberEnd();
     if(fMuonTrackCuts -> IsSelected(mu0)) fpDCA[nummu] = 1;
     nummu++; 
     }  
   }  
 
 // reassign labels to muons belonging to a dimuon surviving cuts  [labels need to be reassigned since only selected tracks are saved and not all of them]
  for(Int_t i = 0;i<numdimu;i++){
    Int_t LabelNew1 = 0;
    Int_t LabelNew2 = 0;
    for(int j=0;j<LabelOld1[i];j++){
      if(GoodMuon[j]) LabelNew1++;    
    }  
    for(int j=0;j<LabelOld2[i];j++){
      if(GoodMuon[j]) LabelNew2++;    
    }  
    fDimuMu[i][0]=LabelNew1;  
    fDimuMu[i][1]=LabelNew2;
   }
  
  fNMuons =nummu;
  fNDimu=numdimu;     
  } // end loop on ntracks !=0   
  fOutputTree->Fill();
  PostData(1,fOutputTree);
  
}
//______________________________________________________________________________
Double_t CostHE(AliAODTrack* Mu0, AliAODTrack* Mu1){
  Double_t EBeam = 6500;
  Double_t mp = 0.93827231;
  Double_t pbeam = TMath::Sqrt(EBeam*EBeam - mp*mp);
  Double_t pla10 = Mu0 -> Px();
  Double_t pla11 = Mu0 -> Py();
  Double_t pla12 = Mu0 -> Pz();
  Double_t e1 = Mu0 -> E();
  Double_t mu1Charge = Mu0 -> Charge();
  Double_t pla20 = Mu1 -> Px();
  Double_t pla21 = Mu1 -> Py();
  Double_t pla22 = Mu1 -> Pz();
  Double_t e2 = Mu1 -> E();
  if(pla10==0 && pla11==0 && pla12==0 && e1==0 && mu1Charge==0 && pla20==0 && pla21==0 && pla22==0 && e2==0.){return -666.;}

  // Fill the Lorentz vector for projectile and target
  // For the moment we consider no crossing angle
  // Projectile runs towards the MUON arm
  TLorentzVector pProjLab(0.,0.,-pbeam,EBeam); // projectile
  TLorentzVector pTargLab(0.,0., pbeam,EBeam); // target
  //
  // --- Get the muons parameters in the LAB frame
  //
  TLorentzVector pMu1Lab(pla10,pla11,pla12,e1);
  TLorentzVector pMu2Lab(pla20,pla21,pla22,e2);
  //
  // --- Obtain the dimuon parameters in the LAB frame
  //
  TLorentzVector pDimuLab = pMu1Lab + pMu2Lab;
  //
  // --- Translate the dimuon parameters in the dimuon rest frame
  //
  TVector3 beta = (-1./pDimuLab.E())*pDimuLab.Vect();
  TLorentzVector pMu1Dimu = pMu1Lab;
  TLorentzVector pMu2Dimu = pMu2Lab;
  pMu1Dimu.Boost(beta);
  pMu2Dimu.Boost(beta);

  TVector3 zaxis;
  zaxis=(pDimuLab.Vect()).Unit();
  //
  // --- Calculation of the polarization angle (Helicity)
  // (angle between mu+ and the z axis defined above)
  //
  Double_t cost;
  if(mu1Charge > 0) {
    cost = zaxis.Dot((pMu1Dimu.Vect()).Unit());
  } else {
    cost = zaxis.Dot((pMu2Dimu.Vect()).Unit());
  }
  return cost;
}
//______________________________________________________________________________
//______________________________________________________________________________
Double_t PhiHE(AliAODTrack* Mu0, AliAODTrack* Mu1){
  // Calculation the Helicity aimuthal angle (adapted from code by R. Arnaldi)
  Double_t EBeam = 6500.;
  if(EBeam <= 0){
    printf("Can not compute phiHE with EBeam=%f\n",EBeam);
    return -999999999;
  }
  Double_t mp = 0.93827231;
  Double_t pbeam = TMath::Sqrt(EBeam*EBeam - mp*mp);
  Double_t pla10 = Mu0 -> Px();
  Double_t pla11 = Mu0 -> Py();
  Double_t pla12 = Mu0 -> Pz();
  Double_t e1 = Mu0 -> E();
  Double_t mu1Charge = Mu0 -> Charge();
  Double_t pla20 = Mu1 -> Px();
  Double_t pla21 = Mu1 -> Py();
  Double_t pla22 = Mu1 -> Pz();
  Double_t e2 = Mu1 -> E();
  //Double_t mu2Charge=Mu1->Charge();
  if(pla10==0 && pla11==0 && pla12==0 && e1==0 && mu1Charge==0 && pla20==0 && pla21==0 && pla22==0 && e2==0.){return -666.;}

  // Fill the Lorentz vector for projectile and target
  // For the moment we consider no crossing angle
  // Projectile runs towards the MUON arm
  TLorentzVector pProjCM(0.,0.,-pbeam,EBeam); // projectile
  TLorentzVector pTargCM(0.,0., pbeam,EBeam); // target
  //
  // --- Get the muons parameters in the CM frame
  //
  TLorentzVector pMu1CM(pla10,pla11,pla12,e1);
  TLorentzVector pMu2CM(pla20,pla21,pla22,e2);
  //
  // --- Obtain the dimuon parameters in the CM frame
  //
  TLorentzVector pDimuCM = pMu1CM + pMu2CM;
  //
  // --- Translate the muon parameters in the dimuon rest frame
  //
  TVector3 zaxis = (pDimuCM.Vect()).Unit();
  //
  // --- Translate the dimuon parameters in the dimuon rest frame
  //
  TVector3 beta = (-1./pDimuCM.E())*pDimuCM.Vect();
  TLorentzVector pMu1Dimu = pMu1CM;
  TLorentzVector pMu2Dimu = pMu2CM;
  pMu1Dimu.Boost(beta);
  pMu2Dimu.Boost(beta);

  TLorentzVector pProjDimu = pProjCM;
  TLorentzVector pTargDimu = pTargCM;
  pProjDimu.Boost(beta);
  pTargDimu.Boost(beta);

  TVector3 yaxis = ((pProjDimu.Vect()).Cross(pTargDimu.Vect())).Unit();
  TVector3 xaxis = (yaxis.Cross(zaxis)).Unit();
  //
  // --- Calculation of the azimuthal angle (Helicity)
  //
   Double_t phi;
   if(mu1Charge>0) phi = TMath::ATan2((pMu1Dimu.Vect()).Dot(yaxis),(pMu1Dimu.Vect()).Dot(xaxis));
   else phi = TMath::ATan2((pMu2Dimu.Vect()).Dot(yaxis),(pMu2Dimu.Vect()).Dot(xaxis));

   return phi;
}
//______________________________________________________________________________
Double_t CostCS(AliAODTrack* Mu0, AliAODTrack* Mu1){
  Double_t EBeam = 6500.;
  Double_t mp = 0.93827231;
  Double_t pbeam = TMath::Sqrt(EBeam*EBeam - mp*mp);
  Double_t pla10 = Mu0 -> Px();
  Double_t pla11 = Mu0 -> Py();
  Double_t pla12 = Mu0 -> Pz();
  Double_t e1 = Mu0 -> E();
  Double_t mu1Charge = Mu0 -> Charge();
  Double_t pla20 = Mu1 -> Px();
  Double_t pla21 = Mu1 -> Py();
  Double_t pla22 = Mu1 -> Pz();
  Double_t e2 = Mu1 -> E();
  Double_t mu2Charge = Mu1 -> Charge();
  if(pla10==0 && pla11==0 && pla12==0 && e1==0 && mu1Charge==0 && pla20==0 && pla21==0 && pla22==0 && e2==0.){return -666.;}

  // Fill the Lorentz vector for projectile and target
  // For the moment we do not consider the crossing angle
  // Projectile runs towards the MUON arm
  TLorentzVector pProjCM(0.,0.,-pbeam,EBeam); // projectile
  TLorentzVector pTargCM(0.,0., pbeam,EBeam); // target
  //
  // --- Get the muons parameters in the CM frame
  //
  TLorentzVector pMu1CM(pla10,pla11,pla12,e1);
  TLorentzVector pMu2CM(pla20,pla21,pla22,e2);
  //
  // --- Obtain the dimuon parameters in the CM frame
  //
  TLorentzVector pDimuCM = pMu1CM + pMu2CM;
  //
  // --- Translate the dimuon parameters in the dimuon rest frame
  //
  TVector3 beta = (-1./pDimuCM.E())*pDimuCM.Vect();
  TLorentzVector pMu1Dimu = pMu1CM;
  TLorentzVector pMu2Dimu = pMu2CM;
  TLorentzVector pProjDimu = pProjCM;
  TLorentzVector pTargDimu = pTargCM;
  pMu1Dimu.Boost(beta);
  pMu2Dimu.Boost(beta);
  pProjDimu.Boost(beta);
  pTargDimu.Boost(beta);
  //
  // --- Determine the z axis for the CS angle
  //
  TVector3 zaxisCS = (((pProjDimu.Vect()).Unit())-((pTargDimu.Vect()).Unit())).Unit();
  //
  // --- Determine the CS angle (angle between mu+ and the z axis defined above)
  //
  Double_t cost;
  if(mu1Charge > 0) {
    cost = zaxisCS.Dot((pMu1Dimu.Vect()).Unit());
    // Theta CS is not properly defined for Like-Sign muons
    if(mu2Charge > 0 && cost<0) cost = -cost;
  } else {
    // Theta CS is not properly defined for Like-Sign muons
    cost = zaxisCS.Dot((pMu2Dimu.Vect()).Unit());
    if(mu2Charge < 0 && cost<0) cost = -cost;
  }
  return cost;
}
//______________________________________________________________________________
Double_t PhiCS(AliAODTrack* Mu0, AliAODTrack* Mu1){
  // Cosinus of the Collins-Soper polar decay angle
  Double_t EBeam = 6500.;
  if(EBeam <= 0){
    printf("Can not compute phiCS with EBeam=%f\n",EBeam);
    return -999999999;
  }
  Double_t mp = 0.93827231;
  Double_t pbeam = TMath::Sqrt(EBeam*EBeam - mp*mp);
  Double_t pla10 = Mu0 -> Px();
  Double_t pla11 = Mu0->Py();
  Double_t pla12 = Mu0->Pz();
  Double_t e1 = Mu0->E();
  Double_t mu1Charge = Mu0 -> Charge();
  Double_t pla20 = Mu1 -> Px();
  Double_t pla21 = Mu1 -> Py();
  Double_t pla22 = Mu1 -> Pz();
  Double_t e2 = Mu1 -> E();
  //Double_t mu2Charge=Mu1->Charge();
  if(pla10==0 && pla11==0 && pla12==0 && e1==0 && mu1Charge==0 && pla20==0 && pla21==0 && pla22==0 && e2==0.){return -666.;}

  // Fill the Lorentz vector for projectile and target
  // For the moment we do not consider the crossing angle
  // Projectile runs towards the MUON arm
  TLorentzVector pProjCM(0.,0.,-pbeam,EBeam); // projectile
  TLorentzVector pTargCM(0.,0., pbeam,EBeam); // target
  //
  // --- Get the muons parameters in the CM frame
  //
  TLorentzVector pMu1CM(pla10,pla11,pla12,e1);
  TLorentzVector pMu2CM(pla20,pla21,pla22,e2);
  //
  // --- Obtain the dimuon parameters in the CM frame
  //
  TLorentzVector pDimuCM = pMu1CM + pMu2CM;
  //
  // --- Translate the dimuon parameters in the dimuon rest frame
  //
  TVector3 beta = (-1./pDimuCM.E())*pDimuCM.Vect();
  TLorentzVector pMu1Dimu = pMu1CM;
  TLorentzVector pMu2Dimu = pMu2CM;
  TLorentzVector pProjDimu = pProjCM;
  TLorentzVector pTargDimu = pTargCM;
  pMu1Dimu.Boost(beta);
  pMu2Dimu.Boost(beta);
  pProjDimu.Boost(beta);
  pTargDimu.Boost(beta);
  //
  // --- Determine the z axis for the CS angle
  //
  TVector3 zaxisCS = (((pProjDimu.Vect()).Unit()) - ((pTargDimu.Vect()).Unit())).Unit();
  //
  // --- Determine the CS angle (angle between mu+ and the z axis defined above)
  //
   TVector3 yaxisCS = (((pProjDimu.Vect()).Unit()).Cross((pTargDimu.Vect()).Unit())).Unit();
   TVector3 xaxisCS = (yaxisCS.Cross(zaxisCS)).Unit();

   Double_t phi;
   if(mu1Charge>0) phi = TMath::ATan2((pMu1Dimu.Vect()).Dot(yaxisCS),((pMu1Dimu.Vect()).Dot(xaxisCS)));
   else phi = TMath::ATan2((pMu2Dimu.Vect()).Dot(yaxisCS),((pMu2Dimu.Vect()).Dot(xaxisCS)));

   return phi;
}

//________________________________________________________________________
void AliAnalysisTaskTree_MCut::Terminate(Option_t *) 
{

 }


