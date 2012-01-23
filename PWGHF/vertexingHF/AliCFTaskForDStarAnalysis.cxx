/**************************************************************************
 * Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
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

/* $Id$ */

//-----------------------------------------------------------------------
//   Class for DStar corrections: 
//   
//   The D0 cutting varibles position in the container and container 
//   binning method from a C.Zampolli example
//   In this way a simple comparison between D0 and D0 from D* is possible. 
//
//-----------------------------------------------------------------------
// Author : A.Grelli, Utrecht University
//
//         a.grelli@uu.nl
//-----------------------------------------------------------------------
 
#include <TCanvas.h>
#include <TParticle.h>
#include <TDatabasePDG.h>
#include <TH1I.h>
#include <TStyle.h>
#include <TFile.h>

#include "AliCFTaskForDStarAnalysis.h"
#include "AliStack.h"
#include "AliMCEvent.h"
#include "AliCFManager.h"
#include "AliCFContainer.h"
#include "AliLog.h"
#include "AliAnalysisManager.h"
#include "AliAODHandler.h"
#include "AliAODEvent.h"
#include "AliAODRecoDecay.h"
#include "AliAODRecoDecayHF.h"
#include "AliAODRecoCascadeHF.h"
#include "AliAODRecoDecayHF2Prong.h"
#include "AliAODMCParticle.h"
#include "AliAODMCHeader.h"
#include "AliESDtrack.h"
#include "TChain.h"
#include "THnSparse.h"
#include "TH2D.h"

//__________________________________________________________________________
AliCFTaskForDStarAnalysis::AliCFTaskForDStarAnalysis() :
  AliAnalysisTaskSE(),
  fCFManager(0x0),
  fHistEventsProcessed(0x0),
  fCorrelation(0x0),
  fCountRecoDStarSel(0),
  fEvents(0),
  fMinITSClusters(5),
  fMinITSClustersSoft(4),
  fAcceptanceUnf(kTRUE)
{
  //
  //Default ctor
  //
}
//___________________________________________________________________________
AliCFTaskForDStarAnalysis::AliCFTaskForDStarAnalysis(const Char_t* name) :
  AliAnalysisTaskSE(name),
  fCFManager(0x0),
  fHistEventsProcessed(0x0),
  fCorrelation(0x0),
  fCountRecoDStarSel(0),
  fEvents(0),
  fMinITSClusters(5),
  fMinITSClustersSoft(4),
  fAcceptanceUnf(kTRUE)
{
  //
  // Constructor. Initialization of Inputs and Outputs
  //
  Info("AliCFTaskForDStarAnalysis","Calling Constructor");
  
  DefineOutput(1,TH1I::Class());
  DefineOutput(2,AliCFContainer::Class());
  DefineOutput(3,THnSparseD::Class());
}

//___________________________________________________________________________
AliCFTaskForDStarAnalysis& AliCFTaskForDStarAnalysis::operator=(const AliCFTaskForDStarAnalysis& c) 
{
  //
  // Assignment operator
  //
  if (this!=&c) {
    AliAnalysisTaskSE::operator=(c) ;
    fCFManager  = c.fCFManager;
    fHistEventsProcessed = c.fHistEventsProcessed;
  }
  return *this;
}

//___________________________________________________________________________
AliCFTaskForDStarAnalysis::AliCFTaskForDStarAnalysis(const AliCFTaskForDStarAnalysis& c) :
  AliAnalysisTaskSE(c),
  fCFManager(c.fCFManager),
  fHistEventsProcessed(c.fHistEventsProcessed),
  fCorrelation(c.fCorrelation),
  fCountRecoDStarSel(c.fCountRecoDStarSel),
  fEvents(c.fEvents),
  fMinITSClusters(c.fMinITSClusters),
  fMinITSClustersSoft(c.fMinITSClustersSoft),
  fAcceptanceUnf(c.fAcceptanceUnf)
{
  //
  // Copy Constructor
  //
}

//___________________________________________________________________________
AliCFTaskForDStarAnalysis::~AliCFTaskForDStarAnalysis() {
  //
  //destructor
  //
  if (fCFManager)           delete fCFManager ;
  if (fHistEventsProcessed) delete fHistEventsProcessed ;
  if (fCorrelation)         delete fCorrelation ;
}

//_________________________________________________
void AliCFTaskForDStarAnalysis::UserExec(Option_t *)
{
  //
  // Main loop function
  //
  
  if (!fInputEvent) {
    Error("UserExec","NO EVENT FOUND!");
    return;
  }
  
  AliAODEvent* aodEvent = dynamic_cast<AliAODEvent*>(fInputEvent);
  
  TClonesArray *arrayDStartoD0pi=0;
  
  if(!aodEvent && AODEvent() && IsStandardAOD()) {
    // In case there is an AOD handler writing a standard AOD, use the AOD 
    // event in memory rather than the input (ESD) event.    
    aodEvent = dynamic_cast<AliAODEvent*> (AODEvent());
    // in this case the braches in the deltaAOD (AliAOD.VertexingHF.root)
    // have to taken from the AOD event hold by the AliAODExtension
    AliAODHandler* aodHandler = (AliAODHandler*) 
      ((AliAnalysisManager::GetAnalysisManager())->GetOutputEventHandler());
    if(aodHandler->GetExtensions()) {
      AliAODExtension *ext = (AliAODExtension*)aodHandler->GetExtensions()->FindObject("AliAOD.VertexingHF.root");
      AliAODEvent *aodFromExt = ext->GetAOD();
      arrayDStartoD0pi=(TClonesArray*)aodFromExt->GetList()->FindObject("Dstar");
    }
  } else {
    arrayDStartoD0pi=(TClonesArray*)aodEvent->GetList()->FindObject("Dstar");
  }
  
  if (!arrayDStartoD0pi) {
    AliError("Could not find array of HF vertices");
    return;
  }
  
  // fix for temporary bug in ESDfilter 
  // the AODs with null vertex pointer didn't pass the PhysSel
  if(!aodEvent->GetPrimaryVertex() || TMath::Abs(aodEvent->GetMagneticField())<0.001) return;

  fEvents++;
  if (fEvents%10000 ==0) AliDebug(2,Form("Event %d",fEvents));

  fCFManager->SetRecEventInfo(aodEvent);
  fCFManager->SetMCEventInfo(aodEvent);
  
  // event selection
  Double_t containerInput[15] ;
  Double_t containerInputMC[15] ;
  
  //loop on the MC event
  
  TClonesArray* mcArray = dynamic_cast<TClonesArray*>(aodEvent->FindListObject(AliAODMCParticle::StdBranchName()));
  if (!mcArray) {
    AliError("Could not find Monte-Carlo in AOD");
    return;
  }
  
  AliAODMCHeader *mcHeader = dynamic_cast<AliAODMCHeader*>(aodEvent->GetList()->FindObject(AliAODMCHeader::StdBranchName()));
  if (!mcHeader) {
    AliError("Could not find MC Header in AOD");
    return;
  }
  
  // AOD primary vertex
  AliAODVertex *vtx1 = (AliAODVertex*)aodEvent->GetPrimaryVertex();
  Double_t zPrimVertex = vtx1->GetZ();
  Double_t zMCVertex = mcHeader->GetVtxZ();
  Bool_t vtxFlag = kTRUE;
  TString title=vtx1->GetTitle();
  if(!title.Contains("VertexerTracks")) vtxFlag=kFALSE;
  
  for (Int_t iPart=0; iPart<mcArray->GetEntriesFast(); iPart++) { 
    AliAODMCParticle* mcPart = dynamic_cast<AliAODMCParticle*>(mcArray->At(iPart));
    if (!mcPart) {
      AliWarning("MC Particle not found in tree, skipping"); 
      continue;
    } 

    // check the MC-level cuts
    if (!fCFManager->CheckParticleCuts(0, mcPart)) continue;  // 0 stands for MC level
    
    // fill the container for Gen-level selection
    Double_t vectorMC[9] = {9999.,9999.,9999.,9999.,9999.,9999.,9999.,9999.,9999.};
    
    if (GetDStarMCParticle(mcPart, mcArray, vectorMC)){
      
      containerInputMC[0] = vectorMC[0];
      containerInputMC[1] = vectorMC[1] ;
      containerInputMC[2] = vectorMC[2] ;
      containerInputMC[3] = vectorMC[3] ;
      containerInputMC[4] = vectorMC[4] ;
      containerInputMC[5] = vectorMC[5] ; 
      containerInputMC[6] = 0.;      
      containerInputMC[7] = 0.;          
      containerInputMC[8] = 0.;
      containerInputMC[9] = -100000.; 
      containerInputMC[10] = 1.01; 
      containerInputMC[11] = vectorMC[6];  
      containerInputMC[12] = zMCVertex;     // z vertex
      containerInputMC[13] = vectorMC[7];   // pt D0 pion
      containerInputMC[14] = vectorMC[8];   // pt D0 kaon  
      
      fCFManager->GetParticleContainer()->Fill(containerInputMC,kStepGenerated);
      
      // check the MC-Acceptance level cuts
      // since standard CF functions are not applicable, using Kine Cuts on daughters
      
      Int_t daughter0 = mcPart->GetDaughter(0);
      Int_t daughter1 = mcPart->GetDaughter(1);
      
      AliDebug(2, Form("daughter0 = %d and daughter1 = %d",daughter0,daughter1));
      
      //D0
      AliAODMCParticle* mcPartDaughter0 = dynamic_cast<AliAODMCParticle*>(mcArray->At(daughter0));
      if (!mcPartDaughter0) {
	AliError("Could not find Monte-Carlo in AOD");
	return;
      }
  
      // Soft Pion
      AliAODMCParticle* mcPartDaughter1 = dynamic_cast<AliAODMCParticle*>(mcArray->At(daughter1));
      if (!mcPartDaughter1) {
	AliError("Could not find Monte-Carlo in AOD");
	return;
      }
      
      // Acceptance variables for the soft pion
      Double_t eta1 = mcPartDaughter1->Eta();
      Double_t pt1  = mcPartDaughter1->Pt(); 
      
      Int_t daughD00 = 0;
      Int_t daughD01 = 0;
      
      // Just to be sure to take the right particles    
      if(TMath::Abs(mcPartDaughter0->GetPdgCode())==421){
	daughD00 = mcPartDaughter0->GetDaughter(0);
	daughD01 = mcPartDaughter0->GetDaughter(1);
      }else{
	daughD00 = mcPartDaughter1->GetDaughter(0);
	daughD01 = mcPartDaughter1->GetDaughter(1);
      }
      
      // the D0 daughters
      AliAODMCParticle* mcD0PartDaughter0 = dynamic_cast<AliAODMCParticle*>(mcArray->At(daughD00));
      AliAODMCParticle* mcD0PartDaughter1 = dynamic_cast<AliAODMCParticle*>(mcArray->At(daughD01));
      
      if (!mcD0PartDaughter0 || !mcD0PartDaughter1) {
	AliWarning("At least one D0 Daughter Particle not found in tree, but it should be, this check was already done...");
	return; 
      }
      
      // D0 daughters - needed for acceptance
      Double_t theD0pt0 =  mcD0PartDaughter0->Pt();
      Double_t theD0pt1 =  mcD0PartDaughter1->Pt();
      Double_t theD0eta0 = mcD0PartDaughter0->Eta();
      Double_t theD0eta1 = mcD0PartDaughter1->Eta();
      
      // ACCEPTANCE REQUESTS ---------
      
      // soft pion 
      Bool_t daught1inAcceptance = (TMath::Abs(eta1) <= 0.9 && pt1 >= 0.05);
      // Do daughters
      Bool_t D0daught0inAcceptance = (TMath::Abs(theD0eta0) <= 0.9 && theD0pt0 >= 0.1); 
      Bool_t D0daught1inAcceptance = (TMath::Abs(theD0eta1) <= 0.9 && theD0pt1 >= 0.1);
      
      if (daught1inAcceptance && D0daught0inAcceptance && D0daught1inAcceptance) {
	
	AliDebug(2, "D* Daughter particles in acceptance");
	
	fCFManager->GetParticleContainer()->Fill(containerInputMC,kStepAcceptance);
	
	// check on the vertex
	if (vtxFlag){
	  // filling the container if the vertex is ok
	  fCFManager->GetParticleContainer()->Fill(containerInputMC,kStepVertex) ;
	  
	  Bool_t refitFlag = kTRUE;
	  for (Int_t iaod =0; iaod<aodEvent->GetNumberOfTracks(); iaod++){
	    AliAODTrack *track = (AliAODTrack*)aodEvent->GetTrack(iaod);
	    
            // refit only for D0 daughters
	    if ((track->GetLabel() == daughD00) || (track->GetLabel() == daughD01)) {
	      if(!(track->GetStatus()&AliESDtrack::kTPCrefit) || !(track->GetStatus()&AliESDtrack::kITSrefit)) {
		refitFlag = kFALSE;
	      }
	    }
	  } 
	  if (refitFlag){ // refit
      
	    fCFManager->GetParticleContainer()->Fill(containerInputMC,kStepRefit);

	  } // end of refit	  
	} // end of vertex
      }	//end of acceptance		
    } // end of MC D*
    else {
      AliDebug(3,"Problems in filling the container");
      continue;
    }
  } // end of MC loop
  
  //rec
  AliDebug(2, Form("Found %d D* candidates",arrayDStartoD0pi->GetEntriesFast()));
  
  //D* and D0 prongs needed to MatchToMC method
  Int_t pdgDgDStartoD0pi[2]={421,211};
  Int_t pdgDgD0toKpi[2]={321,211};
  
  for (Int_t iDStartoD0pi = 0; iDStartoD0pi<arrayDStartoD0pi->GetEntriesFast(); iDStartoD0pi++) {
    
    // D* candidates
    AliAODRecoCascadeHF* dstarD0pi = (AliAODRecoCascadeHF*)arrayDStartoD0pi->At(iDStartoD0pi);
    
    // D0 from the reco cascade
    AliAODRecoDecayHF2Prong* D0particle = (AliAODRecoDecayHF2Prong*)dstarD0pi->Get2Prong();
    Bool_t unsetvtx=kFALSE;
    
    // needed for pointing angle
    if(!D0particle->GetOwnPrimaryVtx()) {
      D0particle->SetOwnPrimaryVtx(vtx1);
      unsetvtx=kTRUE;
    }
    
    // find associated MC particle for D* ->D0toKpi
    Int_t mcLabel = dstarD0pi->MatchToMC(413,421,pdgDgDStartoD0pi,pdgDgD0toKpi,mcArray); 
    
    // find D0->Kpi ... needed in the following
    Int_t mcD0Label = D0particle->MatchToMC(421,mcArray,2,pdgDgD0toKpi); 
    
    if (mcLabel == -1 || mcD0Label ==-1) 
      {
	AliDebug(2,"No MC particle found");
      }
    else {
      
      // the D* and the D0 in MC
      AliAODMCParticle* mcVtxHFDStar = (AliAODMCParticle*)mcArray->At(mcLabel);
      AliAODMCParticle* mcVtxHFD0Kpi = (AliAODMCParticle*)mcArray->At(mcD0Label);
      
      if (!mcVtxHFD0Kpi || !mcVtxHFDStar) {
	AliWarning("Could not find associated MC D0 and/or D* in AOD MC tree");
	continue;
      }
      
      // soft pion
      AliAODTrack *track2 = (AliAODTrack*)dstarD0pi->GetBachelor(); 

      //D0tokpi
      AliAODTrack *track0 = (AliAODTrack*)D0particle->GetDaughter(0);
      AliAODTrack *track1 = (AliAODTrack*)D0particle->GetDaughter(1);
      
      // check if associated MC v0 passes the cuts
      if (!fCFManager->CheckParticleCuts(0 ,mcVtxHFDStar)) { 
	AliDebug(2, "Skipping the particles due to cuts");
	continue; 
      }
      
      // fill the container...
      Double_t pt = TMath::Sqrt(TMath::Power(D0particle->Pt(),2)+TMath::Power(track2->Pt(),2));
      Double_t rapidity = dstarD0pi->YDstar();
      Double_t cosThetaStar = 9999.;
      Double_t pTpi  = 0.;
      Double_t pTK   = 0.;
      Double_t dca   = (D0particle->GetDCA())*1E4;
      Double_t d0pi  = 0.;
      Double_t d0K   = 0.;
      Double_t d0xd0 = (D0particle->Prodd0d0())*1E8;
      Double_t cosPointingAngle = D0particle->CosPointingAngle();
      Double_t phi = D0particle->Phi();
           
      // Select D0 cutting variables
      Int_t pdgCode = mcVtxHFD0Kpi->GetPdgCode();

      // D0 related variables
      if (pdgCode > 0){

	cosThetaStar = D0particle->CosThetaStarD0();
	pTpi    = D0particle->PtProng(0);
	pTK     = D0particle->PtProng(1);
	d0pi    = (D0particle->Getd0Prong(0))*1E4;
	d0K     = (D0particle->Getd0Prong(1))*1E4;
   
      }
      else {

	cosThetaStar = D0particle->CosThetaStarD0bar();
	pTpi    = D0particle->PtProng(1);
	pTK     = D0particle->PtProng(0);
	d0pi    = (D0particle->Getd0Prong(1))*1E4;
	d0K     = (D0particle->Getd0Prong(0))*1E4;
   
      }
      
      // ct of the D0 from D*
      Double_t cT = D0particle->CtD0();
      
      containerInput[0]  = pt;
      containerInput[1]  = rapidity;
      containerInput[2]  = cosThetaStar;
      containerInput[3]  = track2->Pt();
      containerInput[4]  = D0particle->Pt();
      containerInput[5]  = cT*1.E4;  // in micron
      containerInput[6]  = dca;  // in micron
      containerInput[7]  = d0pi;  // in micron
      containerInput[8]  = d0K;  // in micron
      containerInput[9]  = d0xd0;  // in micron^2
      containerInput[10] = cosPointingAngle;  // in micron
      containerInput[11] = phi;  
      containerInput[12] = zPrimVertex;    // z of reconstructed of primary vertex
      containerInput[13] = pTpi;  // D0 pion
      containerInput[14] = pTK;   // D0 kaon  

      fCFManager->GetParticleContainer()->Fill(containerInput,kStepReconstructed) ;   
      
      // refit in ITS and TPC for D0 daughters
      if((!(track0->GetStatus()&AliESDtrack::kTPCrefit)) || (!(track1->GetStatus()&AliESDtrack::kTPCrefit)) || (!(track0->GetStatus()&AliESDtrack::kITSrefit)) || (!(track1->GetStatus()&AliESDtrack::kITSrefit))) {
	continue;
      }

      // reft in ITS for soft pion
      if((!(track2->GetStatus()&AliESDtrack::kITSrefit))) {
	continue;
      }

      // cut in acceptance for the soft pion and for the D0 daughters      
      Bool_t acceptanceProng0 = (TMath::Abs(D0particle->EtaProng(0))<= 0.9 && D0particle->PtProng(0) >= 0.1);
      Bool_t acceptanceProng1 = (TMath::Abs(D0particle->EtaProng(1))<= 0.9 && D0particle->PtProng(1) >= 0.1);
      
      // soft pion acceptance ... is it fine 0.9?????
      Bool_t acceptanceProng2 = (TMath::Abs(track2->Eta())<= 0.9 && track2->Pt() >= 0.05);
      
      if (acceptanceProng0 && acceptanceProng1 && acceptanceProng2) {
	AliDebug(2,"D* reco daughters in acceptance");
	fCFManager->GetParticleContainer()->Fill(containerInput,kStepRecoAcceptance) ;
	
	if(fAcceptanceUnf){	  
	  Double_t fill[4]; //fill response matrix
	  
	  // dimensions 0&1 : pt,eta (Rec)		
	  fill[0] = pt ;
	  fill[1] = rapidity;		
	  // dimensions 2&3 : pt,eta (MC)					
	  fill[2] =  mcVtxHFDStar->Pt();
	  fill[3] =  mcVtxHFDStar->Y();	  
	  fCorrelation->Fill(fill);		
	}  
	
	// cut on the min n. of clusters in ITS for the D0 and soft pion
	Int_t ncls0=0,ncls1=0,ncls2=0;
	for(Int_t l=0;l<6;l++) {
	  if(TESTBIT(track0->GetITSClusterMap(),l)) ncls0++;
	  if(TESTBIT(track1->GetITSClusterMap(),l)) ncls1++;
	  if(TESTBIT(track2->GetITSClusterMap(),l)) ncls2++;
	}	
	// see AddTask for soft pion ITS clusters request
	AliDebug(2, Form("n clusters = %d", ncls0));

	if (ncls0 >= fMinITSClusters && ncls1 >= fMinITSClusters && ncls2>= fMinITSClustersSoft) {
	  fCFManager->GetParticleContainer()->Fill(containerInput,kStepRecoITSClusters) ;
	  
	  // D0 cuts optimized for D* analysis
	  Double_t cuts[7] = {9999999., 1.1, 0., 9999999., 9999999., 0.,0.027}; 

          // needed for cuts	  
          Double_t theD0pt = D0particle->Pt();

	  if (theD0pt <= 1){ // first bin not optimized
	    cuts[0] = 400;
	    cuts[1] = 0.8;
	    cuts[2] = 0.21;
	    cuts[3] = 500;
	    cuts[4] = 500;
	    cuts[5] = -20000;
	    cuts[6] = 0.6;  
	  }
	  else if (theD0pt > 1 && theD0pt <= 2){
	    cuts[0] = 200; 
	    cuts[1] = 0.7; 
	    cuts[2] = 0.8; 
	    cuts[3] = 210;
	    cuts[4] = 210;
	    cuts[5] = -20000;
	    cuts[6] = 0.9;  
	  }
	  else if (theD0pt > 2 && theD0pt <= 3){
	    cuts[0] = 400;
	    cuts[1] = 0.8; 
	    cuts[2] = 0.8;
	    cuts[3] = 420;
	    cuts[4] = 350; 
	    cuts[5] = -8500;
	    cuts[6] = 0.9;   
	  }
	  else if (theD0pt > 3 && theD0pt <= 5){
	    cuts[0] = 160;  
	    cuts[1] = 1.0; 
	    cuts[2] = 1.2;  
	    cuts[3] = 560; 
	    cuts[4] = 420; 
	    cuts[5] = -8500;
	    cuts[6] = 0.9;  
	  }
	  else if (theD0pt > 5){
	    cuts[0] = 800;
	    cuts[1] = 1.0;
	    cuts[2] = 1.2; 
	    cuts[3] = 700;
	    cuts[4] = 700; 
	    cuts[5] = 10000;
	    cuts[6] = 0.9;  
	  }
	  if (dca < cuts[0] 
	      && TMath::Abs(cosThetaStar) < cuts[1]  
	      && pTpi > cuts[2] 
	      && pTK > cuts[2]  
	      && TMath::Abs(d0pi) < cuts[3] 
	      && TMath::Abs(d0K) < cuts[4]  
	      && d0xd0 < cuts[5] 
	      && cosPointingAngle > cuts[6]
	    ){
	    
	    AliDebug(2,"Particle passed D* selection cuts");
	    fCFManager->GetParticleContainer()->Fill(containerInput,kStepRecoCuts) ;   
	    
	    if(!fAcceptanceUnf){ // unfolding
	      
	      Double_t fill[4]; //fill response matrix
	      
	      // dimensions 0&1 : pt,eta (Rec)		    
	      fill[0] = pt ;
	      fill[1] = rapidity;		    
	      // dimensions 2&3 : pt,eta (MC)		    
	      fill[2] =  mcVtxHFDStar->Pt();
	      fill[3] =  mcVtxHFDStar->Y();
	      
	      fCorrelation->Fill(fill);		    
	    }
	  }		
	}
      }
    }
    if(unsetvtx) D0particle->UnsetOwnPrimaryVtx();
  } // end loop on D*->Kpipi
  
  fHistEventsProcessed->Fill(0);
  
  PostData(1,fHistEventsProcessed) ;
  PostData(2,fCFManager->GetParticleContainer()) ;
  PostData(3,fCorrelation) ;
}


//___________________________________________________________________________
void AliCFTaskForDStarAnalysis::Terminate(Option_t*)
{
  // The Terminate()	
  AliAnalysisTaskSE::Terminate();
  
  // draw correlation matrix
  AliCFContainer *cont= dynamic_cast<AliCFContainer*> (GetOutputData(2));
  if(!cont) {
    printf("CONTAINER NOT FOUND\n");
    return;
  } 
}

//___________________________________________________________________________
void AliCFTaskForDStarAnalysis::UserCreateOutputObjects() {

  //
  // useroutput
  //
  
  Info("UserCreateOutputObjects","CreateOutputObjects of task %s\n", GetName());
  
  //slot #1
  OpenFile(1);
  fHistEventsProcessed = new TH1I("CFDSchist0","",1,0,1) ;
}

//________________________________________________________________________________
Bool_t AliCFTaskForDStarAnalysis::GetDStarMCParticle(AliAODMCParticle* mcPart, TClonesArray* mcArray, Double_t* vectorMC)const {
  
  // 
  // fill the D* and D0 MC container
  //
  
  Bool_t isDStar = kFALSE;
  
  if(TMath::Abs(mcPart->GetPdgCode())!=413) return isDStar;
  
  // getting the daughters
  Int_t daughter0 = mcPart->GetDaughter(0);
  Int_t daughter1 = mcPart->GetDaughter(1);
  
  AliDebug(2, Form("daughter0 = %d and daughter1 = %d",daughter0,daughter1));
  if (daughter0 == 0 || daughter1 == 0) {
    AliDebug(2, "Error! the D* MC doesn't have correct daughters!!");
    return isDStar;  
  }
  
  if (TMath::Abs(daughter1 - daughter0) != 1) { // should be everytime true - see PDGdatabooklet
    AliDebug(2, "The D* MC doesn't come from a 2-prong decay, skipping!!");
    return isDStar;  
  }
  
  AliAODMCParticle* mcPartDaughter0 = dynamic_cast<AliAODMCParticle*>(mcArray->At(daughter0));
  AliAODMCParticle* mcPartDaughter1 = dynamic_cast<AliAODMCParticle*>(mcArray->At(daughter1));
  if (!mcPartDaughter0 || !mcPartDaughter1) {
    AliWarning("D*: At least one Daughter Particle not found in tree, skipping"); 
    return isDStar;  
  }
  
  if (!(TMath::Abs(mcPartDaughter0->GetPdgCode())==421 &&
	TMath::Abs(mcPartDaughter1->GetPdgCode())==211) && 
      !(TMath::Abs(mcPartDaughter0->GetPdgCode())==211 &&
	TMath::Abs(mcPartDaughter1->GetPdgCode())==421)) {
    AliDebug(2, "The D* MC doesn't come from a Kpi decay, skipping!!");
    return isDStar;  
  }
  
  Double_t vtx2daughter0[3] = {0,0,0};   // secondary vertex from daughter 0
  Double_t vtx2daughter1[3] = {0,0,0};   // secondary vertex from daughter 1

  // getting vertex from daughters
  mcPartDaughter0->XvYvZv(vtx2daughter0);  
  mcPartDaughter1->XvYvZv(vtx2daughter1);  
  
  // check if the secondary vertex is the same for both
  if (vtx2daughter0[0] != vtx2daughter1[0] && vtx2daughter0[1] != vtx2daughter1[1] && vtx2daughter0[2] != vtx2daughter1[2]) {
    AliError("The D* daughters have different secondary vertex, skipping the track");
    return isDStar;
  }
  
  AliAODMCParticle* neutralDaugh = mcPartDaughter0;
  
  Double_t VectorD0[2] ={0.,0.};
  
  if (!EvaluateIfD0toKpi(neutralDaugh,mcArray,VectorD0)) {
    AliDebug(2, "Error! the D0 MC doesn't have correct daughters!!");
    return isDStar;  
  }
  // get the pT of the daughters
  
  Double_t pTpi = 0.;
  Double_t pTD0 = 0.;
  
  if (TMath::Abs(mcPartDaughter0->GetPdgCode()) == 211) {
    pTpi = mcPartDaughter0->Pt();
    pTD0 = mcPartDaughter1->Pt();
  }
  else {
    pTpi = mcPartDaughter1->Pt();
    pTD0 = mcPartDaughter0->Pt();
  }
  
  vectorMC[0] = mcPart->Pt();
  vectorMC[1] = mcPart->Y() ;
  vectorMC[2] = 0;
  vectorMC[3] = pTpi ;
  vectorMC[4] = pTD0 ;
  vectorMC[5] = 0;
  vectorMC[6] = mcPart->Phi() ;
  vectorMC[7] = VectorD0[0] ;
  vectorMC[8] = VectorD0[1] ;
  
  isDStar = kTRUE;
  
  return isDStar;
}
//________________________________________________________________________________________________

Bool_t AliCFTaskForDStarAnalysis::EvaluateIfD0toKpi(AliAODMCParticle* neutralDaugh, TClonesArray* mcArray, Double_t* VectorD0)const{
  
  //
  // chack wether D0 is decaing into kpi
  //
  
  Bool_t isHadronic = kFALSE;
  
  Int_t daughterD00 = neutralDaugh->GetDaughter(0);
  Int_t daughterD01 = neutralDaugh->GetDaughter(1);
  
  AliDebug(2, Form("daughter0 = %d and daughter1 = %d",daughterD00,daughterD01));
  if (daughterD00 == 0 || daughterD01 == 0) {
    AliDebug(2, "Error! the D0 MC doesn't have correct daughters!!");
    return isHadronic;  
  }
  
  if (TMath::Abs(daughterD01 - daughterD00) != 1) { // should be everytime true - see PDGdatabooklet
    AliDebug(2, "The D0 MC doesn't come from a 2-prong decay, skipping!!");
    return isHadronic;  
  }
  
  AliAODMCParticle* mcPartDaughterD00 = dynamic_cast<AliAODMCParticle*>(mcArray->At(daughterD00));
  AliAODMCParticle* mcPartDaughterD01 = dynamic_cast<AliAODMCParticle*>(mcArray->At(daughterD01));
  if (!mcPartDaughterD00 || !mcPartDaughterD01) {
    AliWarning("D0 MC analysis: At least one Daughter Particle not found in tree, skipping"); 
    return isHadronic;  
  }
  
  if (!(TMath::Abs(mcPartDaughterD00->GetPdgCode())==321 &&
	TMath::Abs(mcPartDaughterD01->GetPdgCode())==211) && 
      !(TMath::Abs(mcPartDaughterD00->GetPdgCode())==211 &&
	TMath::Abs(mcPartDaughterD01->GetPdgCode())==321)) {
    AliDebug(2, "The D0 MC doesn't come from a Kpi decay, skipping!!");
    return isHadronic;  
  }

  Double_t pTD0pi = 0;
  Double_t pTD0K = 0;

 
  if (TMath::Abs(mcPartDaughterD00->GetPdgCode()) == 211) {
    pTD0pi = mcPartDaughterD00->Pt();
    pTD0K = mcPartDaughterD01->Pt();
  }
  else {
    pTD0pi = mcPartDaughterD01->Pt();
    pTD0K  = mcPartDaughterD00->Pt();
  }
 
  isHadronic = kTRUE;

  VectorD0[0] = pTD0pi;
  VectorD0[1] = pTD0K;

  return isHadronic;

}
