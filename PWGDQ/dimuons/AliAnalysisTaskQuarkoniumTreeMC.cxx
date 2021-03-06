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

//-----------------------------------------------------------------------------
// Analysis task to compute muon/dimuon kinematic distributions
// The output is a tree.
// R. Arnaldi
//
//-----------------------------------------------------------------------------

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
#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskSE.h"
#include "AliMultSelection.h" 
#include "AliAODMCParticle.h"
#include "AliAnalysisTaskQuarkoniumTreeMC.h"

Double_t CostHE_MC(AliAODMCParticle*, AliAODMCParticle*);
Double_t CostCS_MC(AliAODMCParticle*, AliAODMCParticle*);
Double_t PhiHE_MC(AliAODMCParticle*, AliAODMCParticle*);
Double_t PhiCS_MC(AliAODMCParticle*, AliAODMCParticle*);
Double_t CostHE_rec(AliAODTrack*, AliAODTrack*);
Double_t CostCS_rec(AliAODTrack*, AliAODTrack*);
Double_t PhiHE_rec(AliAODTrack*, AliAODTrack*);
Double_t PhiCS_rec(AliAODTrack*, AliAODTrack*);

ClassImp(AliAnalysisTaskQuarkoniumTreeMC)
//__________________________________________________________________________
AliAnalysisTaskQuarkoniumTreeMC::AliAnalysisTaskQuarkoniumTreeMC() :
  AliAnalysisTaskSE(),
  fOutputTree(0x0),
  fBeamEnergy(0.),
  fkAnalysisType(0x0),
  fPeriod(0x0),
  fResonance(0x0),
  fNMuons_gen(0x0),
  fNDimu_gen(0x0),
  fNMuons_rec(0x0),
  fNDimu_rec(0x0),
  fAODEvent(0x0)
{
  //
  //Default ctor
  //  
  fVertex[0]=999.; fVertex[1]=999.; fVertex[2]=999.;
  for(Int_t i=0; i<100;i++){
    fPt_rec[i]=999.;
    fE_rec[i]=999.;
    fPx_rec[i]=999; 
    fPy_rec[i]=999; 
    fPz_rec[i]=999; 
    fY_rec[i]=999.; 
    fEta_rec[i]=999.; 
    fMatchTrig_rec[i]=999.; 
    fTrackChi2_rec[i]=999.; 
    fMatchTrigChi2_rec[i]=999.;
    fCharge_rec[i]=999;
    fRAtAbsEnd_rec[i]=999;

  }
  for(Int_t i=0; i<1000;i++){  
    fDimuPt_gen[i]=999.; 
    fDimuPx_gen[i]=999.; 
    fDimuPy_gen[i]=999.; 
    fDimuPz_gen[i]=999.; 
    fDimuY_gen[i]=999.; 
    fDimuMass_gen[i]=999.;
    fDimuCharge_gen[i]=999;
    fDimuMatch_gen[i]=999;
    fDimuCostHE_gen[i] = 999;
    fDimuPhiHE_gen[i] = 999;
    fDimuCostCS_gen[i] = 999;
    fDimuPhiCS_gen[i] = 999;
    
    fDimuPt_rec[i]=999.; 
    fDimuPx_rec[i]=999.; 
    fDimuPy_rec[i]=999.; 
    fDimuPz_rec[i]=999.; 
    fDimuY_rec[i]=999.; 
    fDimuMass_rec[i]=999.;
    fDimuCharge_rec[i]=999;
    fDimuMatch_rec[i]=999;
    fDimuCostHE_rec[i] = 999;
    fDimuPhiHE_rec[i] = 999;
    fDimuCostCS_rec[i] = 999;
    fDimuPhiCS_rec[i] = 999;
    for(Int_t k=0;k<2;k++) fDimuMu_rec[i][k]=999;

  } 
  
}

//__________________________________________________________________________
AliAnalysisTaskQuarkoniumTreeMC::AliAnalysisTaskQuarkoniumTreeMC(const char *name) :
  AliAnalysisTaskSE(name),
  fOutputTree(0x0),
  fBeamEnergy(0.),
  fkAnalysisType(0x0),
  fPeriod(0x0),
  fResonance(0x0),
  fNMuons_gen(0x0),
  fNDimu_gen(0x0),
  fNMuons_rec(0x0),
  fNDimu_rec(0x0),
  fAODEvent(0x0)
{
  //
  // Constructor. Initialization of Inputs and Outputs
  //
  Info("AliAnalysisTaskQuarkoniumTreeMC","Calling Constructor");
  fVertex[0]=999.; fVertex[1]=999.; fVertex[2]=999.;
  for(Int_t i=0; i<100;i++){

    fPt_rec[i]=999.;
    fE_rec[i]=999.;
    fPx_rec[i]=999; 
    fPy_rec[i]=999; 
    fPz_rec[i]=999; 
    fY_rec[i]=999.; 
    fEta_rec[i]=999.; 
    fMatchTrig_rec[i]=999.; 
    fTrackChi2_rec[i]=999.; 
    fMatchTrigChi2_rec[i]=999.;
    fCharge_rec[i]=999;
    fRAtAbsEnd_rec[i]=999;

  }
  for(Int_t i=0; i<1000;i++){  
    fDimuPt_gen[i]=999.; 
    fDimuPx_gen[i]=999.; 
    fDimuPy_gen[i]=999.; 
    fDimuPz_gen[i]=999.; 
    fDimuY_gen[i]=999.; 
    fDimuMass_gen[i]=999.;
    fDimuCharge_gen[i]=999;
    fDimuMatch_gen[i]=999;
    fDimuCostHE_gen[i] = 999;
    fDimuPhiHE_gen[i] = 999;
    fDimuCostCS_gen[i] = 999;
    fDimuPhiCS_gen[i] = 999;
    
    fDimuPt_rec[i]=999.; 
    fDimuPx_rec[i]=999.; 
    fDimuPy_rec[i]=999.; 
    fDimuPz_rec[i]=999.; 
    fDimuY_rec[i]=999.; 
    fDimuMass_rec[i]=999.;
    fDimuCharge_rec[i]=999;
    fDimuMatch_rec[i]=999;
    fDimuCostHE_rec[i] = 999;
    fDimuPhiHE_rec[i] = 999;
    fDimuCostCS_rec[i] = 999;
    fDimuPhiCS_rec[i] = 999;
    for(Int_t k=0;k<2;k++) fDimuMu_rec[i][k]=999;

  } 

  DefineOutput(1,TTree::Class());
}

//___________________________________________________________________________
AliAnalysisTaskQuarkoniumTreeMC& AliAnalysisTaskQuarkoniumTreeMC::operator=(const AliAnalysisTaskQuarkoniumTreeMC& c) 
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
AliAnalysisTaskQuarkoniumTreeMC::AliAnalysisTaskQuarkoniumTreeMC(const AliAnalysisTaskQuarkoniumTreeMC& c) :
  AliAnalysisTaskSE(c),  
  fOutputTree(c.fOutputTree),
  fBeamEnergy(c.fBeamEnergy),
  fkAnalysisType(c.fkAnalysisType),
  fPeriod(c.fPeriod),
  fResonance(c.fResonance),
  fNMuons_gen(c.fNMuons_gen),
  fNDimu_gen(c.fNDimu_gen),
  fNMuons_rec(c.fNMuons_rec),
  fNDimu_rec(c.fNDimu_rec),
  fAODEvent(c.fAODEvent)
 {
  //
  // Copy Constructor									
  //
}

//___________________________________________________________________________
AliAnalysisTaskQuarkoniumTreeMC::~AliAnalysisTaskQuarkoniumTreeMC() {
  //
  //destructor
  //
  Info("~AliAnalysisTaskQuarkoniumTreeMC","Calling Destructor");
  if (AliAnalysisManager::GetAnalysisManager()->GetAnalysisType() != AliAnalysisManager::kProofAnalysis) delete fOutputTree;
}


//___________________________________________________________________________
void AliAnalysisTaskQuarkoniumTreeMC::UserCreateOutputObjects(){
    
  if (fOutputTree) return; //already initialised ADDED
   
  OpenFile(1,"RECREATE");
  fOutputTree = new TTree("MCTree","Data Tree");

  fOutputTree->Branch("NMuons_gen",&fNMuons_gen,"NMuons_gen/I");
  fOutputTree->Branch("NDimu_gen",&fNDimu_gen,"NDimu_gen/I");
  fOutputTree->Branch("DimuPt_gen",fDimuPt_gen,"DimuPt_gen[NDimu_gen]/D");
  fOutputTree->Branch("DimuPx_gen",fDimuPx_gen,"DimuPx_gen[NDimu_gen]/D");
  fOutputTree->Branch("DimuPy_gen",fDimuPy_gen,"DimuPy_gen[NDimu_gen]/D");
  fOutputTree->Branch("DimuPz_gen",fDimuPz_gen,"DimuPz_gen[NDimu_gen]/D");
  fOutputTree->Branch("DimuY_gen",fDimuY_gen,"DimuY_gen[NDimu_gen]/D");
  fOutputTree->Branch("DimuMass_gen",fDimuMass_gen,"DimuMass_gen[NDimu_gen]/D");
  fOutputTree->Branch("DimuCharge_gen",fDimuCharge_gen,"DimuCharge_gen[NDimu_gen]/I");
  fOutputTree->Branch("DimuMatch_gen",fDimuMatch_gen,"DimuMatch_gen[NDimu_gen]/I");
  fOutputTree->Branch("DimuCostHE_gen",fDimuCostHE_gen,"DimuCostHE_gen[NDimu_gen]/D");
  fOutputTree->Branch("DimuPhiHE_gen",fDimuPhiHE_gen,"DimuPhiHE_gen[NDimu_gen]/D");
  fOutputTree->Branch("DimuCostCS_gen",fDimuCostCS_gen,"DimuCostCS_gen[NDimu_gen]/D");
  fOutputTree->Branch("DimuPhiCS_gen",fDimuPhiCS_gen,"DimuPhiCS_gen[NDimu_gen]/D");

  fOutputTree->Branch("NMuons_rec",&fNMuons_rec,"NMuons_rec/I");
  fOutputTree->Branch("Pt_rec",fPt_rec,"Pt_rec[NMuons_rec]/D");
  fOutputTree->Branch("E_rec",fE_rec,"E_rec[NMuons_rec]/D");
  fOutputTree->Branch("Px_rec",fPx_rec,"Px_rec[NMuons_rec]/D");
  fOutputTree->Branch("Py_rec",fPy_rec,"Py_rec[NMuons_rec]/D");
  fOutputTree->Branch("Pz_rec",fPz_rec,"Pz_rec[NMuons_rec]/D");
  fOutputTree->Branch("Y_rec",fY_rec,"Y_rec[NMuons_rec]/D");
  fOutputTree->Branch("Eta_rec",fEta_rec,"Eta_rec[NMuons_rec]/D");
  fOutputTree->Branch("MatchTrig_rec",fMatchTrig_rec,"MatchTrig_rec[NMuons_rec]/I");
  fOutputTree->Branch("TrackChi2_rec",fTrackChi2_rec,"TrackChi2_rec[NMuons_rec]/D");
  fOutputTree->Branch("MatchTrigChi2_rec",fMatchTrigChi2_rec,"MatchTrigChi2_rec[NMuons_rec]/D");
  fOutputTree->Branch("Charge_rec",fCharge_rec,"Charge_rec[NMuons_rec]/I");
  fOutputTree->Branch("RAtAbsEnd_rec",fRAtAbsEnd_rec,"RAtAbsEnd_rec[NMuons_rec]/D");
 
  fOutputTree->Branch("NDimu_rec",&fNDimu_rec,"NDimu_rec/I");
  fOutputTree->Branch("DimuMu_rec",fDimuMu_rec,"DimuMu_rec[NDimu_rec][2]/I");
  fOutputTree->Branch("DimuPt_rec",fDimuPt_rec,"DimuPt_rec[NDimu_rec]/D");
  fOutputTree->Branch("DimuPx_rec",fDimuPx_rec,"DimuPx_rec[NDimu_rec]/D");
  fOutputTree->Branch("DimuPy_rec",fDimuPy_rec,"DimuPy_rec[NDimu_rec]/D");
  fOutputTree->Branch("DimuPz_rec",fDimuPz_rec,"DimuPz_rec[NDimu_rec]/D");
  fOutputTree->Branch("DimuY_rec",fDimuY_rec,"DimuY_rec[NDimu_rec]/D");
  fOutputTree->Branch("DimuMass_rec",fDimuMass_rec,"DimuMass_rec[NDimu_rec]/D");
  fOutputTree->Branch("DimuCharge_rec",fDimuCharge_rec,"DimuCharge_rec[NDimu_rec]/I");
  fOutputTree->Branch("DimuMatch_rec",fDimuMatch_rec,"DimuMatch_rec[NDimu_rec]/I");
  fOutputTree->Branch("DimuCostHE_rec",fDimuCostHE_rec,"DimuCostHE_rec[NDimu_rec]/D");
  fOutputTree->Branch("DimuPhiHE_rec",fDimuPhiHE_rec,"DimuPhiHE_rec[NDimu_rec]/D");
  fOutputTree->Branch("DimuCostCS_rec",fDimuCostCS_rec,"DimuCostCS_rec[NDimu_rec]/D");
  fOutputTree->Branch("DimuPhiCS_rec",fDimuPhiCS_rec,"DimuPhiCS_rec[NDimu_rec]/D");


  fOutputTree->ls(); 

  PostData(1,fOutputTree); 
 
} 

//_________________________________________________
void AliAnalysisTaskQuarkoniumTreeMC::UserExec(Option_t *)
{
  fNMuons_gen=0; 
  fNDimu_gen=0;
  fNMuons_rec=0; 
  fNDimu_rec=0;
  fVertex[0]=999.; fVertex[1]=999.; fVertex[2]=999.;
  for(Int_t i=0; i<100;i++){

    fPt_rec[i]=999.;
    fE_rec[i]=999.;
    fPx_rec[i]=999; 
    fPy_rec[i]=999; 
    fPz_rec[i]=999; 
    fY_rec[i]=999.; 
    fEta_rec[i]=999.; 
    fMatchTrig_rec[i]=999.; 
    fTrackChi2_rec[i]=999.; 
    fMatchTrigChi2_rec[i]=999.;
    fCharge_rec[i]=999;
    fRAtAbsEnd_rec[i]=999;
  }
  for(Int_t i=0; i<1000;i++){  
    fDimuPt_gen[i]=999.; 
    fDimuPx_gen[i]=999.; 
    fDimuPy_gen[i]=999.; 
    fDimuPz_gen[i]=999.; 
    fDimuY_gen[i]=999.; 
    fDimuMass_gen[i]=999.;
    fDimuCharge_gen[i]=999.;
    fDimuMatch_gen[i]=0;
    fDimuCostHE_gen[i] = 999;
    fDimuPhiHE_gen[i] = 999;
    fDimuCostCS_gen[i] = 999;
    fDimuPhiCS_gen[i] = 999;

    fDimuPt_rec[i]=999.; 
    fDimuPx_rec[i]=999.; 
    fDimuPy_rec[i]=999.; 
    fDimuPz_rec[i]=999.; 
    fDimuY_rec[i]=999.; 
    fDimuMass_rec[i]=999.;
    fDimuCharge_rec[i]=999.;
    fDimuCostHE_rec[i] = 999;
    fDimuPhiHE_rec[i] = 999;
    fDimuCostCS_rec[i] = 999;
    fDimuPhiCS_rec[i] = 999;
    for(Int_t k=0;k<2;k++) fDimuMu_rec[i][k]=999;

  } 
//
// Execute analysis for current event
//
  fAODEvent = dynamic_cast<AliAODEvent*> (InputEvent());
  if ( ! fAODEvent ) {
    AliError ("AOD event not found. Nothing done!");
    return;
  }
  
  TClonesArray *mcarray = dynamic_cast<TClonesArray*>(fAODEvent->FindListObject(AliAODMCParticle::StdBranchName()));
  
  AliAODVertex *PrimVertex =  fAODEvent->GetPrimaryVertex();
  fVertex[0]=PrimVertex->GetX();
  fVertex[1]=PrimVertex->GetY();
  fVertex[2]=PrimVertex->GetZ();
  
  Int_t PDGCode;
  if(fResonance == "JPsi") PDGCode = 443;
  if(fResonance == "Psi2S") PDGCode = 100443;
   
  // generated events
  Int_t numdimu_gen = 0;
  for (Int_t i=0;i<mcarray->GetEntries();i++){
     
       AliAODMCParticle *mcp = (AliAODMCParticle *)mcarray->At(i);
       if(mcp->GetPdgCode()==PDGCode){
	   
         fDimuPt_gen[i]=mcp->Pt(); 
         fDimuPx_gen[i]=mcp->Px(); 
         fDimuPy_gen[i]=mcp->Py(); 
         fDimuPz_gen[i]=mcp->Pz(); 
         fDimuY_gen[i]=mcp->Y(); 
         fDimuMass_gen[i]=mcp->M();
         fDimuCharge_gen[i]=mcp->Charge();
         fDimuMatch_gen[i]=0;
	 	 
	 Int_t d0 = mcp->GetDaughterFirst();
	 Int_t d1 = mcp->GetDaughterLast(); //in muon AOD we have only the 2 muons, in AOD there might be 3 tracks
	 Int_t daught[3];
	 
	 AliAODMCParticle *mcp_daughter;
	 Int_t j=0;
	 for(int i=d0;i<d1;i++){
	   mcp_daughter = (AliAODMCParticle *) mcarray->At(i);
	   if(mcp_daughter->GetPdgCode()!=13 && mcp_daughter->GetPdgCode()!=-13) continue; //skip gamma due to radiative production
	   daught[j]=i;
	   j++;
	 }
  	 AliAODMCParticle *mcp_daught0 = (AliAODMCParticle *)mcarray->At(daught[0]);
	 AliAODMCParticle *mcp_daught1 = (AliAODMCParticle *)mcarray->At(daught[1]);
         fDimuCostHE_gen[i] = CostHE_MC(mcp_daught0,mcp_daught1);
         fDimuPhiHE_gen[i] = PhiHE_MC(mcp_daught0,mcp_daught1);
         fDimuCostCS_gen[i] = CostCS_MC(mcp_daught0,mcp_daught1);
         fDimuPhiCS_gen[i] = PhiCS_MC(mcp_daught0,mcp_daught1);

	 numdimu_gen++;
       }
     }

  fNDimu_gen=numdimu_gen;     

  Int_t numdimu = 0;
  Int_t nummu = 0;

//----------------------------------------------------------
// read muon tracks from AliAOD.root (not AliAOD.Muons.root)
//----------------------------------------------------------  
 
  Int_t ntracks = fAODEvent->GetNumberOfTracks();   
  for (Int_t i=0;i<ntracks;i++){
       
    AliAODTrack *mu0=(AliAODTrack*)fAODEvent->GetTrack(i);
    if(mu0->GetLabel()== -1) {
      printf("negative label\n");
      continue;
    }	
    AliAODMCParticle *mctrack0 = (AliAODMCParticle*) mcarray->At(mu0->GetLabel());     
    Int_t mum_num0 = mctrack0->GetMother();
    if(mu0->GetLabel() != mctrack0->GetLabel())continue;
    if(mctrack0->GetMother()<0) continue;
    AliAODMCParticle *mcpart0 = (AliAODMCParticle *)mcarray->At(mum_num0);
    if(mcpart0->GetPdgCode()==PDGCode){
      fCharge_rec[i] = mu0->Charge();
      fPt_rec[i] = mu0->Pt();
      fPx_rec[i] = mu0->Px();
      fPy_rec[i] = mu0->Py();
      fPz_rec[i] = mu0->Pz();
      fY_rec[i]  = mu0->Y();
      fEta_rec[i]= mu0->Eta();
      fE_rec[i] = mu0->E();
      fMatchTrig_rec[i]   = mu0->GetMatchTrigger();
      fMatchTrigChi2_rec[i]= mu0->GetChi2MatchTrigger();
      fRAtAbsEnd_rec[i]=mu0->GetRAtAbsorberEnd();
     
      if (!mu0->IsMuonTrack()) continue;
      for(Int_t j=i+1;j<ntracks;j++){
	 AliAODTrack *mu1=(AliAODTrack*)fAODEvent->GetTrack(j);
	 if (!mu1->IsMuonTrack()) continue;
	 AliAODDimuon *dimu=new AliAODDimuon(mu0,mu1);
	 if(mu1->GetLabel()== -1) {
           printf("negative label\n");
           continue;
          }	

	  AliAODMCParticle *mctrack1 = (AliAODMCParticle*) mcarray->At(mu1->GetLabel());     
 	  Int_t mum_num1 = mctrack1->GetMother();

	  if(mu1->GetLabel() != mctrack1->GetLabel())continue;
	    if(mctrack1->GetMother()<0) continue;
	      AliAODMCParticle *mcpart1 = (AliAODMCParticle *)mcarray->At(mum_num1);
 	      if(mcpart1->GetPdgCode()==PDGCode){
      	 
               fDimuMass_rec[numdimu] = dimu->Mass();
               fDimuPt_rec[numdimu] = dimu->Pt();
               fDimuPx_rec[numdimu] = dimu->Px();
               fDimuPy_rec[numdimu] = dimu->Py();
               fDimuPz_rec[numdimu] = dimu->Pz();
               fDimuY_rec[numdimu] = dimu->Y();
               fDimuCharge_rec[numdimu]= dimu->Charge();
               fDimuMu_rec[numdimu][0]=i;  fDimuMu_rec[numdimu][1]=j;
	
	       fDimuCostHE_rec[numdimu] = CostHE_rec(mu0,mu1);
	       fDimuPhiHE_rec[numdimu] = PhiHE_rec(mu0,mu1);
	       fDimuCostCS_rec[numdimu] = CostCS_rec(mu0,mu1);
	       fDimuPhiCS_rec[numdimu] = PhiCS_rec(mu0,mu1);

               if(mu0->GetMatchTrigger()>1 || mu1->GetMatchTrigger()>1) fDimuMatch_rec[numdimu]=1;
               if(mu0->GetMatchTrigger()>1 && mu1->GetMatchTrigger()>1) fDimuMatch_rec[numdimu]=2;
	
               numdimu++;      
            }
            delete dimu;
          }
          nummu++;
        }
       }
       fNMuons_rec =nummu;
       fNDimu_rec=numdimu;     
       fOutputTree->Fill();
       PostData(1,fOutputTree);
  
}
//______________________________________________________________________________
Double_t CostHE_MC(AliAODMCParticle* Mu0, AliAODMCParticle* Mu1){
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
Double_t PhiHE_MC(AliAODMCParticle* Mu0, AliAODMCParticle* Mu1){
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
Double_t CostCS_MC(AliAODMCParticle* Mu0, AliAODMCParticle* Mu1){
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
Double_t PhiCS_MC(AliAODMCParticle* Mu0, AliAODMCParticle* Mu1){
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


//______________________________________________________________________________
Double_t CostHE_rec(AliAODTrack* Mu0, AliAODTrack* Mu1){
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
Double_t PhiHE_rec(AliAODTrack* Mu0, AliAODTrack* Mu1){
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
Double_t CostCS_rec(AliAODTrack* Mu0, AliAODTrack* Mu1){
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
Double_t PhiCS_rec(AliAODTrack* Mu0, AliAODTrack* Mu1){
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
void AliAnalysisTaskQuarkoniumTreeMC::Terminate(Option_t *) 
{

 }


