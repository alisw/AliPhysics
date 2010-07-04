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

/* $Id:$ */


#include <TChain.h>
#include <TFile.h>
#include <TList.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TProfile.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TTree.h>
#include <TDatabasePDG.h> 

#include "AliAnalysisTaskJetChem.h"
#include "AliAnalysisManager.h"
#include "AliMCEventHandler.h"
#include "AliAODInputHandler.h"
#include "AliAODHandler.h"
#include "AliAODJet.h"
#include "AliAODTrack.h"
#include "AliAODMCParticle.h"

#include "AliGenPythiaEventHeader.h"
#include "AliAnalysisHelperJetTasks.h"
#include "AliESDtrack.h"
#include "AliExternalTrackParam.h"

//#include <iostream> // TEST!!!


//
// Analysis class for jet chemistry studies
// based on AliAnalysisTaskUE by Arian Abrahantes Quintana and Enesto Lopez
// contact: Oliver Busch, o.busch@gsi.de

// This class needs as input AOD with track and Jets
// the output is a list of histograms
//
// AOD can be either connected to the InputEventHandler  
// for a chain of AOD files 
// or 
// to the OutputEventHandler
// for a chain of ESD files, so this case class should be 
// in the train after the Jet finder
//




ClassImp(AliAnalysisTaskJetChem)

////////////////////////////////////////////////////////////////////////


//____________________________________________________________________
AliAnalysisTaskJetChem:: AliAnalysisTaskJetChem(const char* name): AliAnalysisTaskSE(name),
fDebug(kFALSE),
fDeltaAOD(kFALSE),
fDeltaAODBranch(""),
fAODBranch("jets"),
fDeltaAODBranchMC(""),
fAODBranchMC("jetsMC"),
fArrayJetsAOD(0x0), 
fArrayJetsMC(0x0),
fAOD(0x0),            
fAODjets(0x0),
fListOfHistos(0x0),
fJetsOnFly(kFALSE),
fUseLOConeJets(kFALSE),
fUseLOConeMCJets(kFALSE),
fUsePythiaJets(kFALSE),
fConeRadius(0.),
fTrackPtCutJF(0),
fFilterBitJF(0x01),
fRequireITSRefitJF(kFALSE),
fRejectK0TracksJF(kFALSE),
fJetPtCut(5.0),
fJetEtaCut(0.5),
fFilterBit(0x10),
fTrackPtCut(0.),
fTrackEtaCut(0.9),
fUseOnFlyV0s(kFALSE),
fCutnSigdEdx(0.), 
fUseAODMCTracksForUE(kFALSE),
fAreaReg(1.0), 
fAvgTrials(1),
fhPrimVertexNCont(0x0),
fhPrimVertexRho(0x0),
fhPrimVertexZ(0x0),
fhNJets(0x0),
fhNJetsMC(0x0),
fhLeadingEta(0x0),
fhLeadingNTracksVsEta(0x0),
fhLeadingPtVsEta(0x0),
fhLeadingPhi(0x0), 
fhLeadingPt(0x0),
fhLeadingPtDiffr(0x0),
fhLeadingEtaMC(0x0), 
fhLeadingPhiMC(0x0), 
fhLeadingPtMC(0x0),
fhLeadingPtMCDiffr(0x0),
fhPhiEtaTracksNoCut(0x0),
fhPtTracksNoCut(0x0),
fhPhiEtaTracks(0x0),
fhPtTracks(0x0),
fhTrackMult(0x0),
fhEtaMCTracks(0x0),            
fhPhiMCTracks(0x0),            
fhPtMCTracks(0x0),
fhnTracksVsPtLeading(0x0),
fhdNdEtaPhiDist(0x0),        
fhRegionSumPtMaxVsEt(0x0),
fhRegionMultMaxVsEt(0x0),     
fhRegionSumPtMinVsEt(0x0),
fhRegionMultMinVsEt(0x0),     
fhNV0s(0x0),
fhV0onFly(0x0),
fhV0DCADaughters(0x0),
fhV0Radius(0x0),
fhV0DCAToVertex(0x0),
fhV0DCAToVertexK0(0x0),
fhV0InvMassK0(0x0),
fhV0InvMassK0JetEvt(0x0),
fhV0InvMassLambda(0x0),
fhV0InvMassAntiLambda(0x0),
fhV0InvMassLambdaJetEvt(0x0),
fhV0InvMassAntiLambdaJetEvt(0x0),
fhdROpanK0VsPt(0x0),
fhdPhiJetV0(0x0),
fhdPhiJetK0(0x0),
fhdRJetK0(0x0),
fhdNdptV0(0x0),
fhdNdptK0(0x0),
fhPtVsEtaK0(0x0),
fhV0InvMassK0DCA(0x0),
fhV0InvMassK0DCAdEdx(0x0),
fhV0InvMassK0DCAPID(0x0),
fhV0InvMassLambdaDCAdEdx(0x0), 
fhV0InvMassAntiLambdaDCAdEdx(0x0),
fhdNdptK0DCA(0x0),
fhdNdptK0DCAdEdx(0x0),
fhV0InvMassK0Min(0x0),
fhV0InvMassLambdaMin(0x0),
fhV0InvMassAntiLambdaMin(0x0),
fhV0InvMassK0Max(0x0),
fhV0InvMassLambdaMax(0x0),
fhV0InvMassAntiLambdaMax(0x0),
fhV0InvMassK0Jet(0x0),
fhV0InvMassLambdaJet(0x0),
fhV0InvMassAntiLambdaJet(0x0),
fhV0InvMassK0Lambda(0x0), 
fhdNdptK0JetEvt(0x0),
fhdNdzK0(0x0),  
fhdNdzK05to10(0x0),  
fhdNdzK010to20(0x0), 
fhdNdzK020to30(0x0), 
fhdNdzK030to40(0x0), 
fhdNdzK040to60(0x0), 
fhdNdxiK0(0x0),
fhdNdzLambda(0x0),
fhdNdzAntiLambda(0x0),
fhdNdzK0Max(0x0),     
fhdNdxiK0Max(0x0),    
fhdNdzLambdaMax(0x0), 
fhdNdxiLambdaMax(0x0),
fhdNdptK0Max(0x0), 
fhdNdptLambdaMax(0x0),
fhdNdzK0Min(0x0),     
fhdNdxiK0Min(0x0),    
fhdNdzLambdaMin(0x0), 
fhdNdxiLambdaMin(0x0),
fhdNdptK0Min(0x0),
fhdNdptLambdaMin(0x0),
fhdNdzK0Jet(0x0),     
fhdNdxiK0Jet(0x0),    
fhdNdzLambdaJet(0x0), 
fhdNdxiLambdaJet(0x0),
fhdNdptK0Jet(0x0),
fhdNdptLambdaJet(0x0),
fhdEdxVsMomV0(0x0),
fhdEdxVsMomV0pidEdx(0x0),
fhdEdxVsMomV0piPID(0x0),
fhdPhiJetK0MC(0x0),
fhdRJetK0MC(0x0),
fhdRV0MC(0x0),
fhdNdptchPiMCMax(0x0),    
fhdNdptK0MCMax(0x0),      
fhdNdptchKMCMax(0x0),     
fhdNdptpMCMax(0x0),       
fhdNdptpBarMCMax(0x0),    
fhdNdptLambdaMCMax(0x0),  
fhdNdptLambdaBarMCMax(0x0),
fhdNdptchPiMCMin(0x0),    
fhdNdptK0MCMin(0x0),      
fhdNdptchKMCMin(0x0),     
fhdNdptpMCMin(0x0),       
fhdNdptpBarMCMin(0x0),    
fhdNdptLambdaMCMin(0x0),  
fhdNdptLambdaBarMCMin(0x0),
fhdNdptOmegaMCMin(0x0),   
fhdNdptOmegaBarMCMin(0x0),
fhdNdptchPiMCJet(0x0),    
fhdNdptK0MCJet(0x0),      
fhdNdptchKMCJet(0x0),     
fhdNdptpMCJet(0x0),       
fhdNdptpBarMCJet(0x0),    
fhdNdptLambdaMCJet(0x0),  
fhdNdptLambdaBarMCJet(0x0),
fhPIDMC(0x0),
fhPIDMC_quarkEv(0x0),
fhPIDMC_gluonEv(0x0),
fhPIDMCAll(0x0),
fhPIDMCMin(0x0),
fhPIDMCJet(0x0),
fhPIDMCMotherK0(0x0),
fhPIDMCGrandMotherK0(0x0),
fhPIDMCMotherChK(0x0),
fhPIDMCMotherK0Trans(0x0),
fhPIDMCGrandMotherK0Trans(0x0),
fhPIDMCMotherChKTrans(0x0),
fhdNdptgammaMC(0x0),
fhdNdptchPiMC(0x0),
fhdNdptpi0MC(0x0),
fhdNdptK0MC(0x0),
fhdNdptchKMC(0x0),
fhdNdptpMC(0x0),
fhdNdptpBarMC(0x0),
fhdNdptLambdaMC(0x0),
fhdNdptLambdaBarMC(0x0),
fhdNdptOmegaMC(0x0),
fhdNdptOmegaBarMC(0x0),
fhdNdxiMC(0x0),
fhdNdxiK0MC(0x0),
fhdNdxiK0MCJet(0x0),
fhdNdzK0MC(0x0), 
fhdNdzK0MCJet(0x0), 
fhdNdptK0MCJetEvt(0x0),
fhnJetsAODvsMC(0x0),
fhLeadingPtAODvsMC(0x0),
fhLeadingEtaAODvsMC(0x0),
fhLeadingPhiAODvsMC(0x0),
fhnTracksLeadingAODvsMC(0x0),
fhLeadingdRAODMC(0x0),
fhLeadingPtAODvsMCdRcut(0x0),
fhdnTracksVsdPtLeadingAODMC(0x0),
fhnTracksJetVsPtAOD(0x0),
fhnTracksJetVsPtAODquarkEv(0x0),
fhRadiusJetVsPtAOD(0x0), 
fhnTracksJetVsPtMC(0x0),
fhnTracksJetVsPtMCquarkEv(0x0),
fhRadiusJetVsPtMC(0x0),
fhnTracksJetVsPtMCK0(0x0),
fhnTracksJetVsPtMCK0quarkEv(0x0),
fhRadiusJetVsPtMCK0(0x0),
fhnTracksJetVsPtAODK0(0x0),
fhnTracksJetVsPtAODK0quarkEv(0x0),
fhRadiusJetVsPtAODK0(0x0),
fhnTracksJetVsPtAODpKch(0x0),
fhRadiusJetVsPtAODpKch(0x0),
fhPythiaProcess(0x0), 
fhPythiaProcessK0(0x0),
fhPythiaProcessKch(0x0),
fhPythiaProcessp(0x0),
fhPythiaProcesspbar(0x0),
fhdNdzJets5to10(0x0),   
fhdNdzJets10to20(0x0),  
fhdNdzJets20to30(0x0),  
fhdNdzJets30to40(0x0),  
fhdNdzJets40to60(0x0),  
fhdNdxiJets5to10(0x0),  
fhdNdxiJets10to20(0x0),
fhdNdxiJets20to30(0x0), 
fhdNdxiJets30to40(0x0), 
fhdNdxiJets40to60(0x0), 
fhdNdptTracksJetPt5to10(0x0),
fhdNdptTracksJetPt10to20(0x0),
fhdNdptTracksJetPt20to30(0x0),
fhdNdptTracksJetPt30to40(0x0),
fhdNdptTracksJetPt40to60(0x0),
fh1Xsec(0x0),
fh1Trials(0x0),
fpdgdb(0x0){
  // Default constructor 

  fAreaReg = 2*TMath::Pi()/6.0 * 2*fTrackEtaCut;
  fpdgdb = TDatabasePDG::Instance(); 

  // Output slot #1 writes into a TList container, 0 reserved for std AOD output
  DefineOutput(1, TList::Class());
}

//______________________________________________________________
Bool_t AliAnalysisTaskJetChem::UserNotify()
{
  //
  // read the cross sections
  // and number of trials from pyxsec.root
  // 

  fAvgTrials = 1;
  TTree *tree = AliAnalysisManager::GetAnalysisManager()->GetTree();
  Float_t xsection = 0;
  Float_t ftrials  = 1;
  if(tree){
    TFile *curfile = tree->GetCurrentFile();
    if (!curfile) {
      Error("Notify","No current file");
      return kFALSE;
    }
    if(!fh1Xsec||!fh1Trials){
      Printf("%s%d No Histogram fh1Xsec",(char*)__FILE__,__LINE__);
      return kFALSE;
    }
    AliAnalysisHelperJetTasks::PythiaInfoFromFile(curfile->GetName(),xsection,ftrials);
    fh1Xsec->Fill("<#sigma>",xsection);
    // construct average trials 
    Float_t nEntries = (Float_t)tree->GetTree()->GetEntries();
    if(ftrials>=nEntries)fAvgTrials = ftrials/nEntries; 
  }  
  return kTRUE;
}

//____________________________________________________________________
void  AliAnalysisTaskJetChem::UserCreateOutputObjects()
{
  // Create the output container
  //
  AliInfo("UserCreateOutPutObjects()");
  //
  //  Histograms

  OpenFile(1);
  CreateHistos();
  PostData(1, fListOfHistos); // PostData at least once

}

//____________________________________________________________________
void  AliAnalysisTaskJetChem::UserExec(Option_t */*option*/)
{
  // Execute analysis for current event
  //
  if (fDebug > 3) AliInfo( " Processing event..." );

  TObject* handler = AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler();
  
  if( handler && handler->InheritsFrom("AliAODInputHandler") ) {
    fAOD  =  ((AliAODInputHandler*)handler)->GetEvent();
    if (fDebug > 1) AliInfo("  ==== Tracks from AliAODInputHandler");
    // Case when jets are reconstructed on the fly from AOD tracks
    // (the Jet Finder is using the AliJetAODReader) of InputEventHandler
    // and put in the OutputEventHandler AOD. Useful whe you want to reconstruct jets with
    // different parameters to default ones stored in the AOD or to use a different algorithm
    if( fJetsOnFly ) {
      handler = AliAnalysisManager::GetAnalysisManager()->GetOutputEventHandler();
      if( handler && handler->InheritsFrom("AliAODHandler") ) {
        fAODjets = ((AliAODHandler*)handler)->GetAOD();
        if (fDebug > 1) AliInfo("  ==== Jets from AliAODHandler");
      }
    } else {
      fAODjets = fAOD;
      if (fDebug > 1) AliInfo("  ==== Jets from AliAODInputHandler");
    }
  } else {
    handler = AliAnalysisManager::GetAnalysisManager()->GetOutputEventHandler();
    if( handler && handler->InheritsFrom("AliAODHandler") ) {
      fAOD  = ((AliAODHandler*)handler)->GetAOD();
      fAODjets = fAOD;
      if (fDebug > 1) AliInfo("  ==== Tracks and Jets from AliAODHandler");
    } else {
      AliFatal("I can't get any AOD Event Handler");
      return;
    }
  }

  // -------------

  // fetch the pythia header info and get the trials
  AliMCEventHandler* mcHandler = dynamic_cast<AliMCEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
  Float_t nTrials = 1;
  if (mcHandler) {  
    AliMCEvent* mcEvent = mcHandler->MCEvent();
    if (mcEvent) {
      AliGenPythiaEventHeader*  pythiaGenHeader = AliAnalysisHelperJetTasks::GetPythiaEventHeader(mcEvent);
      if(pythiaGenHeader){
	nTrials = pythiaGenHeader->Trials();
      }
    }
  }

  fh1Trials->Fill("#sum{ntrials}",fAvgTrials); 
    
  AnalyseEvent();
  
  // Post the data
  PostData(1, fListOfHistos);
}

//____________________________________________________________________
void  AliAnalysisTaskJetChem::AnalyseEvent()
{
  // check Trigger 

  // TString firedTriggerClasses = fAOD->GetFiredTriggerClasses(); 
  // AliInfo(Form("firedTriggerClasses %s",firedTriggerClasses.Data()));
  // if(firedTriggerClasses.Length() > 0 && !firedTriggerClasses.Contains("CINT1B")) return;  

  AliAODVertex *primVertex = fAOD->GetPrimaryVertex();  // is this from tracks or from SPD ? SPD should only be used if no track vertex
  if(!primVertex){
    AliInfo("no prim Vertex found - skip event");
    fhPrimVertexNCont->Fill(-1);
    return; 
  }
    
  TString primVtxName(primVertex->GetName());
  
  if(primVtxName.CompareTo("TPCVertex",TString::kIgnoreCase) == 1){

    AliInfo("found TPC prim Vertex  - skip event");
    fhPrimVertexNCont->Fill(-1);
    return; 
  }
    
  Int_t nVertContributors = primVertex->GetNContributors();
  fhPrimVertexNCont->Fill(nVertContributors);
  
  //cout<<" prim vertex name "<<primVertex->GetName()<<" nCont "<<nVertContributors<<endl;
  
  if(nVertContributors<1){ // eventually check if not SPD vertex ??? 
    AliInfo("prim vertex no contributors - skip event");
    return;
  }
    
  Double_t vertX = primVertex->GetX();
  Double_t vertY = primVertex->GetY();
  Double_t vertZ = primVertex->GetZ();
  
  Double_t vertRho = TMath::Sqrt(vertX*vertX+vertY*vertY);
  
  fhPrimVertexRho->Fill(vertRho);
  fhPrimVertexZ->Fill(vertZ);
  
  if(TMath::Abs(vertZ)>10){
    AliInfo(Form("prim vertex z=%f - skip event",vertZ));
    return; 
  }
  
  
  Int_t pythiaPID = GetPythiaProcessID();
  
  // ------------------------------------------------
  // Find Leading Jets 1,2,3 
  // (could be skipped if Jets are sort by Pt...)

  //Int_t    index1 = -1;
  //Int_t nTracksLeading = 0;

  Int_t nJetsAOD = 0;
  AliAODJet* leadingJetAOD = NULL; // non-zero if any jet in acc and leading jet survives pt cut 
  Int_t indexLeadingAOD   = -1;
  Double_t ptLeadingAOD   = 0.; 
  Int_t indexMaxRegionAOD = 0; // initialize with '0', '+-1' transverse regions 

  Int_t nJetsMC = 0;
  AliAODJet* leadingJetMC = NULL; // non-zero if leading jet survives acc, pt cut 
  Int_t indexLeadingMC   = -1;
  Double_t ptLeadingMC   = 0.; 
  Int_t indexMaxRegionMC = 0; // initialize with '0', '+-1' transverse regions 

  Double_t   ptLeadingAODAllEta  = 0; 
  AliAODJet* leadingJetAODAllEta = NULL; 
  
  if(fUseLOConeJets)   fArrayJetsAOD  = FindChargedParticleJets(); 
  else{ // Use jets on AOD/DeltaAOD

    if(fDeltaAOD){
      if (fDebug > 1) AliInfo(" ==== Jets From  Delta-AODs !");
      if (fDebug > 1) AliInfo(Form(" ====  Reading Branch: %s  ",fDeltaAODBranch.Data()));
      fArrayJetsAOD = (TClonesArray*)fAODjets->GetList()->FindObject(fDeltaAODBranch.Data());
      if (!fArrayJetsAOD){
	AliFatal(" No jet-array! ");
	return;
      }
    }
    else{
      if (fDebug > 1) AliInfo(" ==== AOD jets: Read Standard-AODs  !");
      if (fDebug > 1) AliInfo(Form(" ====  Reading Branch: %s  ",fAODBranch.Data()));
      
      fArrayJetsAOD = (TClonesArray*)fAODjets->FindListObject(fAODBranch.Data());
    }
  }
  

  if(fUseLOConeMCJets)    fArrayJetsMC = FindChargedParticleJetsMC(); 
  else if(fUsePythiaJets) fArrayJetsMC = GetPythiaJets();
  else{

    if(fDeltaAOD){
      if (fDebug > 1) AliInfo(" ==== MC Jets From  Delta-AODs !");
      if (fDebug > 1) AliInfo(Form(" ====  Reading Branch: %s  ",fDeltaAODBranchMC.Data()));
      fArrayJetsMC = (TClonesArray*)fAODjets->GetList()->FindObject(fDeltaAODBranchMC.Data());
      if (!fArrayJetsMC){
	AliFatal(" No jet-array! ");
	return;
      }
    }
    else{
      if (fDebug > 1) AliInfo(" ==== MC jets: Read Standard-AODs  !");
      if (fDebug > 1) AliInfo(Form(" ====  Reading Branch: %s  ",fAODBranchMC.Data()));
      
      fArrayJetsMC = (TClonesArray*)fAODjets->FindListObject(fAODBranchMC.Data());
    }
  }

  if(fArrayJetsAOD) nJetsAOD = fArrayJetsAOD->GetEntries();
  if(fArrayJetsMC)  nJetsMC  = fArrayJetsMC->GetEntries();

  fhNJets->Fill(nJetsAOD);
  fhNJetsMC->Fill(nJetsMC);
  fhnJetsAODvsMC->Fill(nJetsMC,nJetsAOD);

  if(fDebug>1) AliInfo(Form("AOD %d jets",nJetsAOD)); 


  // for xcheck: find leading jet in large eta range
  for(Int_t i=0; i<nJetsAOD; i++){

    AliAODJet* jet  = (AliAODJet*) fArrayJetsAOD->At(i);
    Double_t jetPt  = jet->Pt();
    
    if(jetPt > ptLeadingAODAllEta){ 
      ptLeadingAODAllEta  = jetPt;
      leadingJetAODAllEta = jet;
    }
  }

  if(leadingJetAODAllEta){
    //cout<<" trackRefs entries "<<leadingJetAODAllEta->GetRefTracks()->GetEntriesFast()<<endl;
    fhLeadingNTracksVsEta->Fill(leadingJetAODAllEta->Eta(),leadingJetAODAllEta->GetRefTracks()->GetEntriesFast());
    fhLeadingPtVsEta->Fill(leadingJetAODAllEta->Eta(),leadingJetAODAllEta->Pt()); 
    
  }

  // find leading jet AOD
  for(Int_t i=0; i<nJetsAOD; ++i){

    AliAODJet* jet  = (AliAODJet*) fArrayJetsAOD->At(i); 
    Double_t jetPt  = jet->Pt();
    Double_t jetEta = jet->Eta();

    if((jetPt > ptLeadingAOD) && (TMath::Abs(jetEta)<fJetEtaCut)){ 
      ptLeadingAOD    = jetPt; 
      indexLeadingAOD = i;
    }
  }

  // find leading jet MC
  for(Int_t i=0; i<nJetsMC; ++i){

    AliAODJet* jet  = (AliAODJet*) fArrayJetsMC->At(i); 
    Double_t jetPt  = jet->Pt();
    Double_t jetEta = jet->Eta();

    if((jetPt > ptLeadingMC) && (TMath::Abs(jetEta)<fJetEtaCut)){
      ptLeadingMC    = jetPt; 
      indexLeadingMC = i;
    }
  }

 
  //cout<<" ptLeadingAOD "<<ptLeadingAOD<<" MC "<<ptLeadingMC<<endl;

  if(indexLeadingAOD>=0){ // event with jet  

    Double_t etaLeadingAOD  = ((AliAODJet*) fArrayJetsAOD->At(indexLeadingAOD))->Eta();
    Double_t phiLeadingAOD  = ((AliAODJet*) fArrayJetsAOD->At(indexLeadingAOD))->Phi();
    Int_t nTracksLeadingAOD = ((AliAODJet*) fArrayJetsAOD->At(indexLeadingAOD))->GetRefTracks()->GetEntriesFast();
    
    if(fDebug>1) AliInfo(Form("\n Pt Leading AOD Jet = %6.1f eta=%5.3f, nTracks %d",ptLeadingAOD,etaLeadingAOD,nTracksLeadingAOD));

    fhLeadingEta->Fill(etaLeadingAOD);
 
    if(TMath::Abs(etaLeadingAOD)<fJetEtaCut){ // leading jet eta cut

      fhnTracksVsPtLeading->Fill(ptLeadingAOD,nTracksLeadingAOD);
      fhLeadingPt->Fill(ptLeadingAOD);
      if(IsDiffractiveEvent(pythiaPID))  fhLeadingPtDiffr->Fill(ptLeadingAOD);

      if(ptLeadingAOD>fJetPtCut){ // leading jet pt cut
	
	leadingJetAOD = (AliAODJet*) fArrayJetsAOD->At(indexLeadingAOD);
	
	fhLeadingPhi->Fill(phiLeadingAOD);
	
	// ----------------------------------------------
	// Find max and min regions
	Double_t sumPtRegionPosit  = 0.;
	Double_t sumPtRegionNegat  = 0.;
	Int_t    nTrackRegionPosit = 0;
	Int_t    nTrackRegionNegat = 0;
      
	Int_t nTracks = fAOD->GetNTracks();
	  
	for (Int_t ipart=0; ipart<nTracks; ++ipart) {
	    
	  AliAODTrack* part = fAOD->GetTrack(ipart);
	  if ( !part->TestFilterBit(fFilterBit) ) continue; // track cut selection
	  if (!part->IsPrimaryCandidate()) continue; // reject whatever is not linked to collision point
	  // PID Selection: Reject everything but hadrons
	  Bool_t isHadron = part->GetMostProbablePID()==AliAODTrack::kPion || 
	    part->GetMostProbablePID()==AliAODTrack::kKaon || 
	    part->GetMostProbablePID()==AliAODTrack::kProton;
	  if (!isHadron ) continue;
	  if (!part->Charge() ) continue; //Only charged
	  if (part->Pt() < fTrackPtCut ) continue;
	  if(TMath::Abs(part->Eta()) > fTrackEtaCut ) continue;
	  
	  TVector3 partVect(part->Px(), part->Py(), part->Pz());
	    
	  Int_t region = IsTrackInsideRegion(leadingJetAOD,&partVect );  

	  if (region > 0) {
	    sumPtRegionPosit += part->Pt();
	    nTrackRegionPosit++;
	  }
	  if (region < 0) {
	    sumPtRegionNegat += part->Pt();
	    nTrackRegionNegat++;
	  }
	} // tracks loop 

	// fill sumPt and mult density for transverse regions
	
	if( sumPtRegionPosit > sumPtRegionNegat ) {
	  FillSumPtRegion( ptLeadingAOD, sumPtRegionPosit/fAreaReg, sumPtRegionNegat/fAreaReg );
	} 
	else {
	  FillSumPtRegion( ptLeadingAOD, sumPtRegionNegat/fAreaReg, sumPtRegionPosit/fAreaReg );
	}
	if (nTrackRegionPosit > nTrackRegionNegat ) {
	  FillMultRegion( ptLeadingAOD, nTrackRegionPosit/fAreaReg, nTrackRegionNegat/fAreaReg);
	} 
	else {
	  FillMultRegion( ptLeadingAOD, nTrackRegionNegat/fAreaReg, nTrackRegionPosit/fAreaReg);
	}
        
	indexMaxRegionAOD = (sumPtRegionPosit > sumPtRegionNegat) ? 1 : -1;
	
      } // leading jet pt cut
    } // leading jet eta cut  
  } // jet event
  

  if(indexLeadingMC>=0){ // event with MC jet  

    Double_t etaLeadingMC  = ((AliAODJet*) fArrayJetsMC->At(indexLeadingMC))->Eta();
    Double_t phiLeadingMC  = ((AliAODJet*) fArrayJetsMC->At(indexLeadingMC))->Phi();
    Int_t nTracksLeadingMC = ((AliAODJet*) fArrayJetsMC->At(indexLeadingMC))->GetRefTracks()->GetEntriesFast();
    
    if(fDebug>1) AliInfo(Form("\n Pt Leading MC Jet = %6.1f eta=%5.3f, nTracks %d",ptLeadingMC,etaLeadingMC,nTracksLeadingMC));

    fhLeadingEtaMC->Fill(etaLeadingMC);
  
    if(TMath::Abs(etaLeadingMC)<fJetEtaCut){ // leading jet eta cut

      fhLeadingPtMC->Fill(ptLeadingMC);
      if(IsDiffractiveEvent(pythiaPID)) fhLeadingPtMCDiffr->Fill(ptLeadingMC);

      if(ptLeadingMC>fJetPtCut){ // leading jet pt cut
	
	leadingJetMC = (AliAODJet*) fArrayJetsMC->At(indexLeadingMC);
   
	fhLeadingPhiMC->Fill(phiLeadingMC); // -pi to pi

	// ----------------------------------------------
	// Find max and min regions
	Double_t sumPtRegionPosit  = 0;
	Double_t sumPtRegionNegat  = 0;
	Int_t    nTrackRegionPosit = 0;
	Int_t    nTrackRegionNegat = 0;

	TClonesArray* farray = (TClonesArray*)fAOD->FindListObject("mcparticles");
	
	Int_t ntrks = farray->GetEntries();
	if (fDebug>1) AliInfo(Form("In UE MC analysis tracks %d \n",ntrks));
	  
	for(Int_t i =0 ; i < ntrks; i++){   
	    
	  AliAODMCParticle* mctrk = (AliAODMCParticle*)farray->At(i);
	  //Cuts
	  if (!(mctrk->IsPhysicalPrimary())) continue;
	  //if (!(mctrk->IsPrimary())) continue;
	  
	  if (mctrk->Charge() == 0 || mctrk->Charge()==-99) continue;
	    
	  if (mctrk->Pt() < fTrackPtCut ) continue;
	  if( TMath::Abs(mctrk->Eta()) > fTrackEtaCut ) continue;
	  
	  Bool_t isHadron = TMath::Abs(mctrk->GetPdgCode())==211 ||
	    TMath::Abs(mctrk->GetPdgCode())==2212 ||
	    TMath::Abs(mctrk->GetPdgCode())==321;
	  
	  if (!isHadron) continue;
	    
	  TVector3 partVect(mctrk->Px(), mctrk->Py(), mctrk->Pz());
	  
	  Int_t region = IsTrackInsideRegion(leadingJetMC,&partVect );  
	  
	  if (region > 0) {
	    sumPtRegionPosit += mctrk->Pt();
	    nTrackRegionPosit++;
	  }
	  if (region < 0) {
	    sumPtRegionNegat += mctrk->Pt();
	    nTrackRegionNegat++;
	  }
	} // AliAODMCParticle loop
      	
	indexMaxRegionMC = (sumPtRegionPosit > sumPtRegionNegat) ? 1 : -1;
	
      } // leading jet pt
    } // leading jet eta 
  } // found jet

 
  Bool_t foundK0AOD   = kFALSE; // event with K0 ? 
  Bool_t foundK0MC    = kFALSE; // event with K0 ? 

  CheckV0s(leadingJetAOD,indexMaxRegionAOD,foundK0AOD); // here leadingJetAOD/MC nonzero if jet passes eta & pt cut
  CheckMCParticles(leadingJetMC,indexMaxRegionMC,foundK0MC);
  CompLeadingJets(leadingJetAOD,leadingJetMC,pythiaPID,foundK0AOD,foundK0MC);
  FillReferencePlotsTracks();
  FillReferenceFF(leadingJetAOD);
  

  if(fUseLOConeJets && fArrayJetsAOD){
    fArrayJetsAOD->Delete(); // no 'Clear': AliAODjet contains TMomentum and TRefArray
    delete fArrayJetsAOD;
  }
  if(fUseLOConeMCJets && fArrayJetsMC){
    fArrayJetsMC->Delete(); // no 'Clear': AliAODjet contains TMomentum and TRefArray
    delete fArrayJetsMC;
  }
}

// __________________________________________________________________


Double_t AliAnalysisTaskJetChem::GetJetRadius(const AliAODJet* jet, const Double_t energyFrac){

  // calc jet radius containing fraction energyFrac of full jet pt  
  
  const Int_t kArraySize = 1000;
  const Int_t kInitVal   = -999;

  Int_t nTracks = jet->GetRefTracks()->GetEntriesFast();

  if(nTracks>kArraySize){
    AliError(Form("nTracks in jet %d exceeds max array size",nTracks));
    return -1;
  }


  Double_t deltaR[kArraySize];
  Double_t pt[kArraySize];
  Int_t index[kArraySize];
  for(int i=0; i<kArraySize; i++) index[i] = kInitVal;
  

  Double_t ptTot = 0;

  for(int i=0; i<nTracks; i++){

    AliAODTrack* track = (AliAODTrack*) jet->GetRefTracks()->At(i);

    TLorentzVector *mom4Jet   = jet->MomentumVector();
    TVector3 mom3Jet = mom4Jet->Vect();

    Double_t trackMom[3];
    track->PxPyPz(trackMom);
    TVector3 mom3Track(trackMom);

    Double_t dR = mom3Jet.DeltaR(mom3Track);

    deltaR[i] = dR;
    pt[i]     = track->Pt();

    ptTot += pt[i];
  }

  //cout<<" ptTot "<<ptTot<<" jetPt "<<jet->Pt()<<endl; // Xcheck 

  TMath::Sort(nTracks,deltaR,index,kFALSE); // sort in decreasing order 

  Double_t ptSum = 0;

  for(int i=0; i<nTracks; i++){

    Int_t ind = index[i];
    
    Double_t ptTrack = pt[ind];
    Double_t dR      = deltaR[ind];

    ptSum += ptTrack;

    if(ptSum >= ptTot*energyFrac) return dR;
  }

  return -1;

}

//____________________________________________________________________
void AliAnalysisTaskJetChem::FillSumPtRegion( Double_t leadingE, Double_t ptMax, Double_t ptMin  )
{
  // Fill sumPt of control regions
  
  fhRegionSumPtMaxVsEt->Fill( leadingE, ptMax );
  fhRegionSumPtMinVsEt->Fill( leadingE, ptMin );
}

//____________________________________________________________________
void AliAnalysisTaskJetChem::FillMultRegion(Double_t leadingE, Double_t nTrackPtmax, Double_t nTrackPtmin)
{
  // Fill Nch multiplicity of control regions
  
  fhRegionMultMaxVsEt->Fill( leadingE, nTrackPtmax );
  fhRegionMultMinVsEt->Fill( leadingE, nTrackPtmin );

}

//____________________________________________________________________
Int_t AliAnalysisTaskJetChem::IsTrackInsideRegion(const AliAODJet* aodjetVect,const TVector3 *partVect) 
{  
  // return de region in delta phi
  // -1 negative delta phi 
  //  1 positive delta phi
  //  0 outside region

  TLorentzVector* jetVectLorentz = aodjetVect->MomentumVector(); 
  TVector3 jetVect = jetVectLorentz->Vect();

  static const Double_t k60rad  = 60.*TMath::Pi()/180.;
  static const Double_t k120rad = 120.*TMath::Pi()/180.;
  
  Int_t region = 0;
  if( TMath::Abs(partVect->Eta()) > fTrackEtaCut ) return 0;
  // transverse regions
  if (jetVect.DeltaPhi(*partVect) < -k60rad && jetVect.DeltaPhi(*partVect) > -k120rad ) region = -1;
  if (jetVect.DeltaPhi(*partVect) > k60rad && jetVect.DeltaPhi(*partVect) < k120rad ) region = 1;
    
  return region;
}

//____________________________________________________________________
TClonesArray* AliAnalysisTaskJetChem::FindChargedParticleJetsMC()
{
  // CDF jet finder: 
  // loop over pt-ordered list of tracks, combine tracks within jet cone, recalc cone axis after each step
  // based on implementation by Arian Abrahantes Quintana and Enesto Lopez
 
  TClonesArray* farray = (TClonesArray*)fAOD->FindListObject("mcparticles");
  if(!farray){
    AliInfo("no mcparticles branch"); 
    return 0;
  }

  Int_t nTracks = farray->GetEntries();

  if( !nTracks ) return 0;
  TObjArray tracks(nTracks);

  for (Int_t ipart=0; ipart<nTracks; ++ipart) {

    AliAODMCParticle* part = (AliAODMCParticle*)farray->At(ipart);

    if(!part->IsPhysicalPrimary()) continue;

    // exclude neutrinos
    Int_t   pdg  = TMath::Abs(part->GetPdgCode());
    if((pdg == 12 || pdg == 14 || pdg == 16)) return kFALSE;

    if( !part->Charge() ) continue; // comment / uncomment here
    fhEtaMCTracks->Fill(part->Eta());

    if(TMath::Abs(part->Eta()) > fTrackEtaCut) continue; 
    fhPtMCTracks->Fill(part->Pt());

    if( part->Pt() < fTrackPtCutJF ) continue;
    fhPhiMCTracks->Fill(part->Phi());

    tracks.AddLast(part);
  
  }

  QSortTracks( tracks, 0, tracks.GetEntriesFast() );
  
  nTracks = tracks.GetEntriesFast();

  if( !nTracks ) return 0;
  TObjArray *jets = new TObjArray(nTracks);
  TIter itrack(&tracks);
  while( nTracks ) {
    // 2- Start with the highest pT particle ...
    Float_t px,py,pz,pt; 

    AliAODMCParticle* track = (AliAODMCParticle*)itrack.Next();

    if( !track ) continue;
    px = track->Px();
    py = track->Py();
    pz = track->Pz();
    pt = track->Pt(); 
    jets->AddLast( new AliAODJet(px,py,pz,pt) ); // Use the energy member to store Pt
    AliAODJet* jet = (AliAODJet*)jets->Last();
    jet->AddTrack(track);
    tracks.Remove(track);
   

    TVector2 jetVect2(jet->Px(),jet->Py());

    // 3- Go to the next highest pT particle not already included...
    AliAODMCParticle* track1;
    while ( (track1  = (AliAODMCParticle*)(itrack.Next())) ) {
      TVector2 trackVect2(track1->Px(),track1->Py());
      Double_t dphi = trackVect2.DeltaPhi(jetVect2);

      Double_t r = TMath::Sqrt( (jet->Eta()-track1->Eta())*(jet->Eta()-track1->Eta()) +
                               dphi*dphi );

      if( r < fConeRadius ) {
        Double_t fPt   = jet->E()+track1->Pt();  // Scalar sum of Pt
        // recalculating the centroid
        Double_t eta = jet->Eta()*jet->E()/fPt + track1->Eta()*track1->Pt()/fPt;
 
	jetVect2.SetMagPhi(jet->E()/fPt,jetVect2.Phi());
	trackVect2.SetMagPhi(track1->Pt()/fPt,trackVect2.Phi());
	
	TVector2 sumVect2 = jetVect2+trackVect2;
	Double_t phi = sumVect2.Phi();

        //jet->SetPtEtaPhiE( 1., eta, phi, fPt );
	((TLorentzVector*) jet->MomentumVector())->SetPtEtaPhiE(1,eta,phi,fPt);

	jet->AddTrack(track1);
        tracks.Remove(track1);
      }
    }
    
    tracks.Compress();

    nTracks = tracks.GetEntries();
    //   4- Continue until all particles are in a jet.
    itrack.Reset();
  } // end while nTracks
  
  // Convert to AODjets....
  Int_t njets = jets->GetEntriesFast();

  TClonesArray* aodjets = new TClonesArray("AliAODJet",njets);
  aodjets->SetOwner(kTRUE);

  Int_t count = 0;
  for(Int_t ijet=0; ijet<njets; ++ijet) {
    AliAODJet* jet = (AliAODJet*)jets->At(ijet);

    if (jet->E() < fJetPtCut) continue;
    Float_t px, py,pz,en; // convert to 4-vector
    px = jet->E() * TMath::Cos(jet->Phi());  // Pt * cos(phi)
    py = jet->E() * TMath::Sin(jet->Phi());  // Pt * sin(phi)
    pz = jet->E() / TMath::Tan(2.0 * TMath::ATan(TMath::Exp(-jet->Eta())));
    en = TMath::Sqrt(px * px + py * py + pz * pz);
    jet->SetPxPyPzE(px,py,pz,en);


    TClonesArray &tmpaodjets = *aodjets;
    new(tmpaodjets[count++]) AliAODJet(*jet);
    //aodjets->AddLast( new AliAODJet(*jet));    
    //aodjets->AddLast( new AliAODJet(px, py, pz, en) );
  }
  // Order jets according to their pT .
  QSortTracks( *aodjets, 0, aodjets->GetEntriesFast() );
  
  // debug
  //if (fDebug>3) AliInfo(Form(" %d Charged jets found\n",njets));
  
  jets->Delete(); // OB - should I cleanup or leave it to garbage collection ? 
  delete jets;

  return aodjets;
}
 
// ___________________________________________________________________

Bool_t AliAnalysisTaskJetChem::IsTrackFromK0(const Int_t indexTrack){

  // check wether track with index is from K0 (= V0 with proper inv mass)

  TClonesArray* tracks = fAOD->GetTracks();

  for(int i=0; i<fAOD->GetNumberOfV0s(); i++){ // loop over V0s
    
    AliAODv0* v0 = fAOD->GetV0(i);
    
    Bool_t isOnFly = v0->GetOnFlyStatus();

    if( (fUseOnFlyV0s && !isOnFly) || (!fUseOnFlyV0s && isOnFly) ) continue;

    Double_t massK0 = v0->MassK0Short();
    
    // FIXME: here should go cuts (at least V0 quality)

    if(IsK0InvMass(massK0)){
      
      AliAODTrack *trackPos = (AliAODTrack *) (v0->GetSecondaryVtx()->GetDaughter(0));
      AliAODTrack *trackNeg = (AliAODTrack *) (v0->GetSecondaryVtx()->GetDaughter(1));   

      Int_t indexPos = tracks->IndexOf(trackPos);
      Int_t indexNeg = tracks->IndexOf(trackNeg);
     

      if(indexPos == indexTrack){ return kTRUE;}
      if(indexNeg == indexTrack){ return kTRUE;}

    }
  }

  return kFALSE;

}

//____________________________________________________________________
TClonesArray*  AliAnalysisTaskJetChem::FindChargedParticleJets()
{

  // CDF jet finder: 
  // loop over pt-ordered list of tracks, combine tracks within jet cone, recalc cone axis after each step 
  // based on implementation by Arian Abrahantes Quintana and Enesto Lopez
  // ref: PHYSICAL REVIEW D 65 092002, CDF Collaboration
 

  //  1 - Order all charged particles according to their pT .
  Int_t nTracks = fAOD->GetNTracks();

  if( !nTracks ) return 0;
  TObjArray tracks(nTracks);

  // standard cuts + ITS refit = stdrd PWG4 cut (compose them for productions with old ESDFilter task) 

  for (Int_t ipart=0; ipart<nTracks; ++ipart) {
    AliAODTrack* part = fAOD->GetTrack( ipart );

    UInt_t status = part->GetStatus();

    if( !part->TestFilterBit(fFilterBitJF) ) continue; // track cut selection
    if(fRequireITSRefitJF &&  ((status&AliESDtrack::kITSrefit)==0)) continue;  

    if(!part->Charge() ) continue;

    if(TMath::Abs(part->Eta()) > fTrackEtaCut) continue; 
    
    if(fRejectK0TracksJF && IsTrackFromK0(ipart)) continue; 

    if( part->Pt() < fTrackPtCutJF ) continue;

    tracks.AddLast(part);
  }

  QSortTracks( tracks, 0, tracks.GetEntriesFast() );
  
  nTracks = tracks.GetEntriesFast();
  
  if( !nTracks ) return 0;
  TObjArray *jets = new TObjArray(nTracks);
  TIter itrack(&tracks);
  while( nTracks ) {
    // 2- Start with the highest pT particle ...
    Float_t px,py,pz,pt; 
    AliAODTrack* track = (AliAODTrack*)itrack.Next();
    if( !track ) continue;
    px = track->Px();
    py = track->Py();
    pz = track->Pz();
    pt = track->Pt(); 
    //jets->AddLast( new TLorentzVector(px, py, pz, pt) ); // Use the energy member to store Pt
    jets->AddLast( new AliAODJet(px,py,pz,pt) ); // Use the energy member to store Pt
   //TLorentzVector* jet = (TLorentzVector*)jets->Last();
    AliAODJet* jet = (AliAODJet*)jets->Last();
    jet->AddTrack(track);
    tracks.Remove( track );

    TVector2 jetVect2(jet->Px(),jet->Py()); // OB

    // 3- Go to the next highest pT particle not already included...
    AliAODTrack* track1;
    while ( (track1  = (AliAODTrack*)(itrack.Next())) ) {
      //Double_t dphi = TVector2::Phi_mpi_pi(jet->Phi()-track1->Phi()); // OB remove
      TVector2 trackVect2(track1->Px(),track1->Py());
      Double_t dphi = trackVect2.DeltaPhi(jetVect2);

      Double_t r = TMath::Sqrt( (jet->Eta()-track1->Eta())*(jet->Eta()-track1->Eta()) +
                               dphi*dphi );

      if( r < fConeRadius ) {
        Double_t fPt   = jet->E()+track1->Pt();  // Scalar sum of Pt
        // recalculating the centroid
        Double_t eta = jet->Eta()*jet->E()/fPt + track1->Eta()*track1->Pt()/fPt;

        //Double_t phi = jet->Phi()*jet->E()/fPt + track1->Phi()*track1->Pt()/fPt; // OB - remove
	// OB - recalc phi via weighted 2vectors
	jetVect2.SetMagPhi(jet->E()/fPt,jetVect2.Phi());
	trackVect2.SetMagPhi(track1->Pt()/fPt,trackVect2.Phi());
	
	TVector2 sumVect2 = jetVect2+trackVect2;
	Double_t phi = sumVect2.Phi();

        //jet->SetPtEtaPhiE( 1., eta, phi, fPt );
	((TLorentzVector*) jet->MomentumVector())->SetPtEtaPhiE(1,eta,phi,fPt);

	jet->AddTrack(track1);
        tracks.Remove(track1);
      }
    }
    
    tracks.Compress();

    nTracks = tracks.GetEntries();
    //   4- Continue until all particles are in a jet.
    itrack.Reset();
  } // end while nTracks
  
  // Convert to AODjets....
  Int_t njets = jets->GetEntriesFast();

  TClonesArray* aodjets = new TClonesArray("AliAODJet",njets);
  aodjets->SetOwner(kTRUE);
  
  Int_t count = 0;

  for(Int_t ijet=0; ijet<njets; ++ijet) {
    //TLorentzVector* jet = (TLorentzVector*)jets->At(ijet);
    AliAODJet* jet = (AliAODJet*)jets->At(ijet);

    if (jet->E() < fJetPtCut) continue;
    Float_t px, py,pz,en; // convert to 4-vector
    px = jet->E() * TMath::Cos(jet->Phi());  // Pt * cos(phi)
    py = jet->E() * TMath::Sin(jet->Phi());  // Pt * sin(phi)
    pz = jet->E() / TMath::Tan(2.0 * TMath::ATan(TMath::Exp(-jet->Eta())));
    en = TMath::Sqrt(px * px + py * py + pz * pz);
    jet->SetPxPyPzE(px,py,pz,en);

    TClonesArray &tmpaodjets = *aodjets;
    new(tmpaodjets[count++]) AliAODJet(*jet);
    //aodjets->AddLast( new AliAODJet(*jet));    
    //aodjets->AddLast( new AliAODJet(px, py, pz, en) );
  }
  // Order jets according to their pT .
  //QSortTracks( *aodjets, 0, aodjets->GetEntriesFast() ); // not for TClonesArray
  
  if (fDebug>3) AliInfo(Form(" %d Charged jets found - after cuts %d \n",njets,count));
  
  jets->Delete(); // OB: cleanup
  delete jets;


  return aodjets;
}

//____________________________________________________________________

TClonesArray*  AliAnalysisTaskJetChem::GetPythiaJets(){

  // return Pythia jets in jet acc

  // note: present verion of AliMCEventHandler expects "AliAOD.root" (not "AliAODs.root"), otherwise galice.root will not be found in proper dir
  AliMCEventHandler* mcHandler = dynamic_cast<AliMCEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
  if (!mcHandler) {
    Printf("ERROR: Could not retrieve MC event handler");
    return 0;
  }
  
  AliMCEvent* mcEvent = mcHandler->MCEvent();
  if (!mcEvent) {
    Printf("ERROR: Could not retrieve MC event");
    return 0;
  }
  
  
  AliGenPythiaEventHeader*  pythiaGenHeader = AliAnalysisHelperJetTasks::GetPythiaEventHeader(mcEvent);
  if(!pythiaGenHeader){
    Printf("ERROR: Could not retrieve pythiaGenHeader");
    return 0;
  }
  
  //Get Jets from MC header
  Int_t nPythiaGenJets = pythiaGenHeader->NTriggerJets();
  
  TClonesArray* aodjets = new TClonesArray("AliAODJet",nPythiaGenJets);
  aodjets->SetOwner(kTRUE);

  int count = 0;
  for(int ip = 0; ip<nPythiaGenJets; ++ip){
    Float_t p[4];
    pythiaGenHeader->TriggerJet(ip,p);
    TVector3 tempVect(p[0],p[1],p[2]);
    if ( TMath::Abs(tempVect.Eta())>fJetEtaCut ) continue;

    TClonesArray &tmpaodjets = *aodjets;
    new(tmpaodjets[count++]) AliAODJet(p[0],p[1],p[2],p[3]);
    //aodjets->AddLast(new AliAODJet(p[0],p[1],p[2],p[3]));
   }
  
  return aodjets;
}


//____________________________________________________________________
void  AliAnalysisTaskJetChem::QSortTracks(TObjArray &a, Int_t first, Int_t last)
{
  // Sort array of TObjArray of tracks by Pt using a quicksort algorithm.
  
  static TObject *tmp;
  static int i;           // "static" to save stack space
  int j;
  
  while (last - first > 1) {
    i = first;
    j = last;
    for (;;) {
      while (++i < last && ((AliVParticle*)a[i])->Pt() > ((AliVParticle*)a[first])->Pt() )
        ;
      while (--j > first && ((AliVParticle*)a[j])->Pt() < ((AliVParticle*)a[first])->Pt() )
        ;
      if (i >= j)
        break;
      
      tmp  = a[i];
      a[i] = a[j];
      a[j] = tmp;
    }
    if (j == first) {
      ++first;
      continue;
    }
    tmp = a[first];
    a[first] = a[j];
    a[j] = tmp;
    if (j - first < last - (j + 1)) {
      QSortTracks(a, first, j);
      first = j + 1;   // QSortTracks(j + 1, last);
    } else {
      QSortTracks(a, j + 1, last);
      last = j;        // QSortTracks(first, j);
    }
  }
}
   
//------------------------------------------------------------------

TH1F* AliAnalysisTaskJetChem::CreatePIDhisto(const char* name){
  
  // create histogram
  
  TH1F* result = new TH1F(name,"",60,0,60);
  result->SetOption("E");

  // bin equal Geant ID

  result->GetXaxis()->SetBinLabel(result->GetXaxis()->FindBin(1),"photon");
  result->GetXaxis()->SetBinLabel(result->GetXaxis()->FindBin(2),"e+");
  result->GetXaxis()->SetBinLabel(result->GetXaxis()->FindBin(3),"e-");
  result->GetXaxis()->SetBinLabel(result->GetXaxis()->FindBin(4),"e-neutrino"); 
  result->GetXaxis()->SetBinLabel(result->GetXaxis()->FindBin(5),"mu+");
  result->GetXaxis()->SetBinLabel(result->GetXaxis()->FindBin(6),"mu-");
  result->GetXaxis()->SetBinLabel(result->GetXaxis()->FindBin(7),"pi0");
  result->GetXaxis()->SetBinLabel(result->GetXaxis()->FindBin(8),"pi+");
  result->GetXaxis()->SetBinLabel(result->GetXaxis()->FindBin(9),"pi-");
  result->GetXaxis()->SetBinLabel(result->GetXaxis()->FindBin(10),"K long");
  result->GetXaxis()->SetBinLabel(result->GetXaxis()->FindBin(11),"K+");
  result->GetXaxis()->SetBinLabel(result->GetXaxis()->FindBin(12),"K-");
  result->GetXaxis()->SetBinLabel(result->GetXaxis()->FindBin(13),"n");
  result->GetXaxis()->SetBinLabel(result->GetXaxis()->FindBin(14),"p");
  result->GetXaxis()->SetBinLabel(result->GetXaxis()->FindBin(15),"anti-proton");
  result->GetXaxis()->SetBinLabel(result->GetXaxis()->FindBin(16),"K short");
  result->GetXaxis()->SetBinLabel(result->GetXaxis()->FindBin(17),"eta");
  result->GetXaxis()->SetBinLabel(result->GetXaxis()->FindBin(18),"Lambda");
  result->GetXaxis()->SetBinLabel(result->GetXaxis()->FindBin(19),"Sigma+");
  result->GetXaxis()->SetBinLabel(result->GetXaxis()->FindBin(20),"Sigma0");
  result->GetXaxis()->SetBinLabel(result->GetXaxis()->FindBin(21),"Sigma-");
  result->GetXaxis()->SetBinLabel(result->GetXaxis()->FindBin(22),"Xi0");
  result->GetXaxis()->SetBinLabel(result->GetXaxis()->FindBin(23),"Xi-");
  result->GetXaxis()->SetBinLabel(result->GetXaxis()->FindBin(24),"Omega-");
  result->GetXaxis()->SetBinLabel(result->GetXaxis()->FindBin(25),"anti-neutron");
  result->GetXaxis()->SetBinLabel(result->GetXaxis()->FindBin(26),"anti-Lambda");
  result->GetXaxis()->SetBinLabel(result->GetXaxis()->FindBin(27),"anti-Sigma-");
  result->GetXaxis()->SetBinLabel(result->GetXaxis()->FindBin(28),"anti-Sigma0");
  result->GetXaxis()->SetBinLabel(result->GetXaxis()->FindBin(29),"anti-Sigma+"); 
  result->GetXaxis()->SetBinLabel(result->GetXaxis()->FindBin(30),"anti-Xi0");
  result->GetXaxis()->SetBinLabel(result->GetXaxis()->FindBin(31),"anti-Xi+");
  result->GetXaxis()->SetBinLabel(result->GetXaxis()->FindBin(32),"anti-Omega+");
  result->GetXaxis()->SetBinLabel(result->GetXaxis()->FindBin(33),"tau+");
  result->GetXaxis()->SetBinLabel(result->GetXaxis()->FindBin(34),"tau-");
  result->GetXaxis()->SetBinLabel(result->GetXaxis()->FindBin(35),"D+");
  result->GetXaxis()->SetBinLabel(result->GetXaxis()->FindBin(36),"D-");
  result->GetXaxis()->SetBinLabel(result->GetXaxis()->FindBin(37),"D0");
  result->GetXaxis()->SetBinLabel(result->GetXaxis()->FindBin(38),"anti-D0");
  result->GetXaxis()->SetBinLabel(result->GetXaxis()->FindBin(39),"Ds+");
  result->GetXaxis()->SetBinLabel(result->GetXaxis()->FindBin(40),"anti Ds-");
  result->GetXaxis()->SetBinLabel(result->GetXaxis()->FindBin(41),"Lamba_c+");
  result->GetXaxis()->SetBinLabel(result->GetXaxis()->FindBin(42),"W+");
  result->GetXaxis()->SetBinLabel(result->GetXaxis()->FindBin(43),"W-");
  result->GetXaxis()->SetBinLabel(result->GetXaxis()->FindBin(44),"Z0");
  result->GetXaxis()->SetBinLabel(result->GetXaxis()->FindBin(45),"d");
  result->GetXaxis()->SetBinLabel(result->GetXaxis()->FindBin(46),"t");
  result->GetXaxis()->SetBinLabel(result->GetXaxis()->FindBin(47),"alpha");
  result->GetXaxis()->SetBinLabel(result->GetXaxis()->FindBin(48),"G_nu");
  result->GetXaxis()->SetBinLabel(result->GetXaxis()->FindBin(kPDGpm311Bin),"K0/#bar{K0}");
  result->GetXaxis()->SetBinLabel(result->GetXaxis()->FindBin(kPDG333Bin),"phi");
  result->GetXaxis()->SetBinLabel(result->GetXaxis()->FindBin(kPDGpm313Bin),"K*(892)0");
  result->GetXaxis()->SetBinLabel(result->GetXaxis()->FindBin(kPDGp323Bin),"K*(892)+");
  result->GetXaxis()->SetBinLabel(result->GetXaxis()->FindBin(kPDGm323Bin),"K*(892)-");
  result->GetXaxis()->SetBinLabel(result->GetXaxis()->FindBin(kPDGNeutrinoBin),"nu");
  result->GetXaxis()->SetBinLabel(result->GetXaxis()->FindBin(kPDGCharmedBaryonBin),"charmed baryon");
  result->GetXaxis()->SetBinLabel(result->GetXaxis()->FindBin(kPDGQuarkBin),"q/#bar{q}");
  result->GetXaxis()->SetBinLabel(result->GetXaxis()->FindBin(kPDGDiQuarkBin),"q #bar{q}");
  result->GetXaxis()->LabelsOption("v"); // "u" ?
  
  return result;
}

//------------------------------------------------------------------

TH1F*  AliAnalysisTaskJetChem::CreatePythiaIDhisto(const char* name){
  
  // create histogram
  
  TH1F* result = new TH1F(name,"",22,0,22);
  result->SetOption("E");

  result->GetXaxis()->SetBinLabel(kPythiaPIDP11Bin,"qq #rightarrow qq");             // ISUB = 11 
  result->GetXaxis()->SetBinLabel(kPythiaPIDP12Bin,"q#bar{q} #rightarrow q#bar{q}"); // ISUB = 12
  result->GetXaxis()->SetBinLabel(kPythiaPIDP13Bin,"q#bar{q} #rightarrow gg");       // ISUB = 13
  result->GetXaxis()->SetBinLabel(kPythiaPIDP28Bin,"qg #rightarrow qg ");            // ISUB = 28 
  result->GetXaxis()->SetBinLabel(kPythiaPIDP53Bin,"gg #rightarrow q#bar{q}");       // ISUB = 53
  result->GetXaxis()->SetBinLabel(kPythiaPIDP68Bin,"gg #rightarrow gg");             // ISUB = 68
  result->GetXaxis()->SetBinLabel(kPythiaPIDP92Bin,"SD");                            // ISUB = 92
  result->GetXaxis()->SetBinLabel(kPythiaPIDP93Bin,"SD");                            // ISUB = 93
  result->GetXaxis()->SetBinLabel(kPythiaPIDP94Bin,"DD");                            // ISUB = 94
  result->GetXaxis()->SetBinLabel(kPythiaPIDP95Bin,"low pt (MPI)");                  // ISUB = 95
  result->GetXaxis()->SetBinLabel(kPythiaPIDPOtherBin,"other");                      // ISUB = XX

  result->GetXaxis()->LabelsOption("v"); // "u" ?
  
  return result;
}


//------------------------------------------------------------------

void AliAnalysisTaskJetChem::FillPythiaIDhisto(TH1F* h, const Int_t PID){

  // fill Pyhtia PID histogram

  Int_t bin = -1; 
 
  if(PID == 11)      bin = kPythiaPIDP11Bin;
  else if(PID == 12) bin = kPythiaPIDP12Bin;
  else if(PID == 13) bin = kPythiaPIDP13Bin;
  else if(PID == 28) bin = kPythiaPIDP28Bin;
  else if(PID == 53) bin = kPythiaPIDP53Bin;
  else if(PID == 68) bin = kPythiaPIDP68Bin;
  else if(PID == 92) bin = kPythiaPIDP92Bin;
  else if(PID == 93) bin = kPythiaPIDP93Bin;
  else if(PID == 94) bin = kPythiaPIDP94Bin;
  else if(PID == 95) bin = kPythiaPIDP95Bin;
  else{ 
    if(PID != -1) AliInfo(Form("unknown PID %d",PID));
    return;
  }
  
  h->Fill(h->GetBinCenter(bin),1);
}


//____________________________________________________________________


void  AliAnalysisTaskJetChem::FillPIDhisto(TH1F* hist, Int_t pdg, Float_t weight){
  
  // convert pdg code to Geant ID and fill corresponding bin 

  Int_t fGID = fpdgdb->ConvertPdgToGeant3(pdg);

  //cout<<" pdg "<<pdg<<" fGID "<<fGID<<endl;

  if(TMath::Abs(pdg) == 311) fGID = kPDGpm311Bin; 
  if(pdg == 333)             fGID = kPDG333Bin; 
  if(TMath::Abs(pdg) == 313) fGID = kPDGpm313Bin; 
  if(pdg == 323)             fGID = kPDGp323Bin; 
  if(pdg == -323)            fGID = kPDGm323Bin; 

  if(TMath::Abs(pdg)==12 || TMath::Abs(pdg)==14 || TMath::Abs(pdg)==16) fGID = kPDGNeutrinoBin;

  if(TMath::Abs(pdg)==4122) fGID = kPDGCharmedBaryonBin;

  if(1<=TMath::Abs(pdg) && TMath::Abs(pdg)<=6) fGID = kPDGQuarkBin; 

  if(TMath::Abs(pdg)==1103 || TMath::Abs(pdg)==2101 || TMath::Abs(pdg)==2103 || TMath::Abs(pdg)==2203 || 
     TMath::Abs(pdg)==3101 || TMath::Abs(pdg)==3103 || TMath::Abs(pdg)==3201 || TMath::Abs(pdg)==3203 || 
     TMath::Abs(pdg)==3303 || TMath::Abs(pdg)==4101 || TMath::Abs(pdg)==4103 || TMath::Abs(pdg)==4201 || 
     TMath::Abs(pdg)==4203 || TMath::Abs(pdg)==4301 || TMath::Abs(pdg)==4303 || TMath::Abs(pdg)==4403 || 
     TMath::Abs(pdg)==5101 || TMath::Abs(pdg)==5103 || TMath::Abs(pdg)==5201 || TMath::Abs(pdg)==5203 || 
     TMath::Abs(pdg)==5301 || TMath::Abs(pdg)==5303 || TMath::Abs(pdg)==5401 || TMath::Abs(pdg)==5403 || 
     TMath::Abs(pdg)==5503)  fGID = kPDGDiQuarkBin;
    

  hist->Fill(fGID,weight);

  if(fGID>hist->GetBinLowEdge(hist->GetNbinsX()+1) || 
     fGID<hist->GetBinLowEdge(1)){
    
    AliError(Form("fGID %d for pdg %d exceeding histo limits ",fGID,pdg));
  }

  if(fGID == 0){
    AliInfo(Form("fGID 0 for pdg %d ",pdg));
  }

}

//____________________________________________________________________
void  AliAnalysisTaskJetChem::CreateHistos(){

  // make histos
  
  fListOfHistos = new TList();
  fListOfHistos->SetOwner(kTRUE);  

  fhPrimVertexNCont = new TH1F("hPrimVertexNCont","",52,-2,50);
  fhPrimVertexNCont->SetXTitle("");
  fhPrimVertexNCont->Sumw2();
  fListOfHistos->Add( fhPrimVertexNCont ); 

  fhPrimVertexRho = new TH1F("hPrimVertexRho","",100,0,1);
  fhPrimVertexRho->SetXTitle("");
  fhPrimVertexRho->Sumw2();
  fListOfHistos->Add( fhPrimVertexRho ); 

  fhPrimVertexZ = new TH1F("hPrimVertexZ","",40,-20,20);
  fhPrimVertexZ->SetXTitle("");
  fhPrimVertexZ->Sumw2();
  fListOfHistos->Add( fhPrimVertexZ ); 

  fhNJets = new TH1F("hNJets","",10, 0, 10);
  fhNJets->SetXTitle("# of jets");
  fhNJets->Sumw2();
  fListOfHistos->Add( fhNJets ); 

  fhNJetsMC = new TH1F("hNJetsMC","",10,0,10);
  fhNJetsMC->SetXTitle("# of jets");
  fhNJetsMC->Sumw2();
  fListOfHistos->Add( fhNJetsMC ); 

  fhLeadingEta = new TH1F("hLeadingEta","Leading Jet eta",12,-0.6,0.6);
  fhLeadingEta->SetXTitle("eta");
  fhLeadingEta->SetYTitle("dN/deta");
  fhLeadingEta->Sumw2();
  fListOfHistos->Add(fhLeadingEta); 

  fhLeadingNTracksVsEta = new TH2F("hLeadingNTracksVsEta","",20,-1.0,1.0,20,0,20);
  fhLeadingNTracksVsEta->SetXTitle("eta");
  fhLeadingNTracksVsEta->SetYTitle("# of tracks");
  fhLeadingNTracksVsEta->Sumw2();
  fListOfHistos->Add(fhLeadingNTracksVsEta); 

  fhLeadingPtVsEta = new TH2F("hLeadingPtVsEta","",20,-1.0,1.0,50,0,50);
  fhLeadingPtVsEta->SetXTitle("eta");
  fhLeadingPtVsEta->SetYTitle("# of tracks");
  fhLeadingPtVsEta->Sumw2();
  fListOfHistos->Add(fhLeadingPtVsEta); 


  fhLeadingPhi = new TH1F("hLeadingPhi","Leading Jet phi",63,0,6.3);
  fhLeadingPhi->SetXTitle("phi");
  fhLeadingPhi->SetYTitle("dN/dphi");
  fhLeadingPhi->Sumw2();
  fListOfHistos->Add(fhLeadingPhi);


  fhLeadingPt  = new TH1F("hLeadingPt","leading Jet p_{T}",50,0,50);
  fhLeadingPt->SetXTitle("p_{T} (GeV/c)");
  fhLeadingPt->SetYTitle("dN/dp_{T} (1/GeV)");
  fhLeadingPt->Sumw2();
  fListOfHistos->Add( fhLeadingPt ); 

  fhLeadingPtDiffr  = new TH1F("hLeadingPtDiffr","leading Jet p_{T}",50,0,50);
  fhLeadingPtDiffr->SetXTitle("P_{T} (GeV/c)");
  fhLeadingPtDiffr->SetYTitle("dN/dp_{T} (1/GeV)");
  fhLeadingPtDiffr->Sumw2();
  fListOfHistos->Add( fhLeadingPtDiffr ); 

  fhLeadingEtaMC = new TH1F("hLeadingEtaMC","Leading Jet eta",12,-0.6,0.6);
  fhLeadingEtaMC->SetXTitle("eta");
  fhLeadingEtaMC->SetYTitle("dN/deta");
  fhLeadingEtaMC->Sumw2();
  fListOfHistos->Add(fhLeadingEtaMC); 

  fhLeadingPhiMC = new TH1F("hLeadingPhiMC","Leading Jet phi",63,0,6.3);
  fhLeadingPhiMC->SetXTitle("phi");
  fhLeadingPhiMC->SetYTitle("dN/dphi");
  fhLeadingPhiMC->Sumw2();
  fListOfHistos->Add(fhLeadingPhiMC);

  fhLeadingPtMC  = new TH1F("hLeadingPtMC","Leading Jet p_{T}",50,0,50);
  fhLeadingPtMC->SetXTitle("p_{T} (GeV/c)");
  fhLeadingPtMC->SetYTitle("dN/dp_{T} (1/GeV)");
  fhLeadingPtMC->Sumw2();
  fListOfHistos->Add( fhLeadingPtMC ); 

  fhLeadingPtMCDiffr  = new TH1F("hLeadingPtMCDiffr","Leading Jet p_{T}",50,0,50);
  fhLeadingPtMCDiffr->SetXTitle("p_{T} (GeV/c)");
  fhLeadingPtMCDiffr->SetYTitle("dN/dp_{T} (1/GeV)");
  fhLeadingPtMCDiffr->Sumw2();
  fListOfHistos->Add( fhLeadingPtMCDiffr ); 

  fhPhiEtaTracksNoCut = new TH2F("hPhiEtaTracksNoCut","phi vs eta tracks",20,-1.0,1.0,63,0,6.3);
  fhPhiEtaTracksNoCut->SetXTitle("eta");
  fhPhiEtaTracksNoCut->SetYTitle("phi");
  fhPhiEtaTracksNoCut->Sumw2();
  fListOfHistos->Add(fhPhiEtaTracksNoCut);

  fhPtTracksNoCut = new TH1F("hPtTracksNoCut","p_{T} tracks",150,0,150);
  fhPtTracksNoCut->SetXTitle("p_{T} (GeV)");
  fhPtTracksNoCut->SetYTitle("dN/dp_{T} (1/GeV)");
  fhPtTracksNoCut->Sumw2();
  fListOfHistos->Add(fhPtTracksNoCut);

  fhPhiEtaTracks = new TH2F("hPhiEtaTracks","phi vs eta tracks",20,-1.0,1.0,63,0,6.3);
  fhPhiEtaTracks->SetXTitle("eta");
  fhPhiEtaTracks->SetYTitle("phi");
  fhPhiEtaTracks->Sumw2();
  fListOfHistos->Add(fhPhiEtaTracks);

  fhPtTracks = new TH1F("hPtTracks","p_{T} tracks",150,0,150);
  fhPtTracks->SetXTitle("P{T} (GeV)");
  fhPtTracks->SetYTitle("dN/dp_{T} (1/GeV)");
  fhPtTracks->Sumw2();
  fListOfHistos->Add(fhPtTracks);

  fhTrackMult = new TH1F("hTrackMult","",150,0,150);
  fhTrackMult->SetXTitle("n Tracks");
  fhTrackMult->SetYTitle("counts");
  fhTrackMult->Sumw2();
  fListOfHistos->Add(fhTrackMult);

  fhEtaMCTracks = new TH1F("hEtaMCTracks","eta tracks",30,-1.5,1.5);
  fhEtaMCTracks->SetXTitle("eta");
  fhEtaMCTracks->SetYTitle("dN/deta");
  fhEtaMCTracks->Sumw2();
  fListOfHistos->Add(fhEtaMCTracks);

  fhPhiMCTracks = new TH1F("hPhiMCTracks","phi tracks",63,0,6.3);
  fhPhiMCTracks->SetXTitle("phi");
  fhPhiMCTracks->SetYTitle("dN/dphi");
  fhPhiMCTracks->Sumw2();
  fListOfHistos->Add(fhPhiMCTracks);

  fhPtMCTracks = new TH1F("hPtMCTracks","p_{T} tracks",50,0,50);
  fhPtMCTracks->SetXTitle("p_{T} (GeV)");
  fhPtMCTracks->SetYTitle("dN/dp_{T} (1/GeV)");
  fhPtMCTracks->Sumw2();
  fListOfHistos->Add(fhPtMCTracks);

  fhnTracksVsPtLeading  = new TH2F("hnTracksVsPtLeading","",50,0.,50.,20,-0.5,19.5);
  fhnTracksVsPtLeading->SetXTitle("p_{T} (GeV/c)");
  fhnTracksVsPtLeading->SetYTitle("n tracks");
  fhnTracksVsPtLeading->Sumw2();
  fListOfHistos->Add( fhnTracksVsPtLeading );

  fhdNdEtaPhiDist  = new TH1F("hdNdEtaPhiDist","Charge particle density |#eta|< 1 vs #Delta#phi",  120, 0.,   2.*TMath::Pi());
  fhdNdEtaPhiDist->SetXTitle("#Delta#phi");
  fhdNdEtaPhiDist->SetYTitle("dN{ch}/d#etad#phi");
  fhdNdEtaPhiDist->Sumw2();
  fListOfHistos->Add( fhdNdEtaPhiDist );

  fhRegionSumPtMaxVsEt = new TH1F("hRegionSumPtMaxVsEt","P{T}^{90, max} vs Leading Jet P{T}",50,0.,50.);
  fhRegionSumPtMaxVsEt->SetXTitle("P{T} (GeV/c)");
  fhRegionSumPtMaxVsEt->Sumw2();
  fListOfHistos->Add( fhRegionSumPtMaxVsEt ); 

  fhRegionMultMaxVsEt = new TH1F("hRegionMultMaxVsEt","N{ch}^{90, max} vs Leading Jet P{T}",50,0.,50.);
  fhRegionMultMaxVsEt->SetXTitle("E (GeV hRegionAveSumPtVsEt/c)");
  fhRegionMultMaxVsEt->Sumw2();
  fListOfHistos->Add( fhRegionMultMaxVsEt );

  fhRegionSumPtMinVsEt = new TH1F("hRegionSumPtMinVsEt","P{T}^{90, min} vs Leading Jet P{T}",50,0.,50.);
  fhRegionSumPtMinVsEt->SetXTitle("P{T} (GeV/c)");
  fhRegionSumPtMinVsEt->Sumw2();
  fListOfHistos->Add( fhRegionSumPtMinVsEt );
    
  fhRegionMultMinVsEt = new TH1F("hRegionMultMinVsEt","N{ch}^{90, min} vs Leading Jet P{T}",50,0.,50.);
  fhRegionMultMinVsEt->SetXTitle("E (GeV/c)");
  fhRegionMultMinVsEt->Sumw2();
  fListOfHistos->Add( fhRegionMultMinVsEt );
  
  // V0s 

  fhNV0s = new TH1F("hNV0s","n V0s",50,0,50);
  fhNV0s->SetXTitle("n V0s");
  fhNV0s->Sumw2();
  fListOfHistos->Add(fhNV0s);

  fhV0onFly = new TH1F("hV0onFly","on-the-fly V0",5,0,5);
  fhV0onFly->SetXTitle("is on-the-fly V0");
  fhV0onFly->Sumw2();
  fListOfHistos->Add(fhV0onFly);

  fhV0DCADaughters = new TH1F("hV0DCADaughters","V0 DCA daughters",200,0,2.0);
  fhV0DCADaughters->SetXTitle("V0 DCA daughters");
  fhV0DCADaughters->Sumw2();
  fListOfHistos->Add(fhV0DCADaughters); 

  fhV0Radius = new TH1F("hV0Radius","V0 radius",2500,0,250);
  fhV0Radius->SetXTitle("V0 radius");
  fhV0Radius->Sumw2();
  fListOfHistos->Add(fhV0Radius); 

  fhV0DCAToVertex = new TH1F("hV0DCAToVertex","",100,0,10);
  fhV0DCAToVertex->SetXTitle("V0 DCA (cm)");
  fhV0DCAToVertex->Sumw2();
  fListOfHistos->Add(fhV0DCAToVertex);

  fhV0DCAToVertexK0 = new TH1F("hV0DCAToVertexK0","",100,0,10); 
  fhV0DCAToVertexK0->SetXTitle("V0 DCA (cm)");
  fhV0DCAToVertexK0->Sumw2();
  fListOfHistos->Add(fhV0DCAToVertexK0);

  fhV0InvMassK0 = new TH1F("hV0InvMassK0","",2000,0,2);
  fhV0InvMassK0->SetXTitle("inv mass (GeV)");
  fhV0InvMassK0->Sumw2();
  fListOfHistos->Add(fhV0InvMassK0); 

  fhV0InvMassK0JetEvt = new TH1F("hV0InvMassK0JetEvt","",2000,0,2);
  fhV0InvMassK0JetEvt->SetXTitle("inv mass (GeV)");
  fhV0InvMassK0JetEvt->Sumw2();
  fListOfHistos->Add(fhV0InvMassK0JetEvt); 

  fhV0InvMassLambda = new TH1F("hV0InvMassLambda","",2000,0,2);
  fhV0InvMassLambda->SetXTitle("inv mass (GeV)");
  fhV0InvMassLambda->Sumw2();
  fListOfHistos->Add(fhV0InvMassLambda);

  fhV0InvMassAntiLambda = new TH1F("hV0InvMassAntiLambda","",2000,0,2);
  fhV0InvMassAntiLambda->SetXTitle("inv mass (GeV)");
  fhV0InvMassAntiLambda->Sumw2();
  fListOfHistos->Add(fhV0InvMassAntiLambda);

  fhV0InvMassLambdaJetEvt = new TH1F("hV0InvMassLambdaJetEvt","",2000,0,2);
  fhV0InvMassLambdaJetEvt->SetXTitle("inv mass (GeV)");
  fhV0InvMassLambdaJetEvt->Sumw2();
  fListOfHistos->Add(fhV0InvMassLambdaJetEvt);

  fhV0InvMassAntiLambdaJetEvt = new TH1F("hV0InvMassAntiLambdaJetEvt","",2000,0,2);
  fhV0InvMassAntiLambdaJetEvt->SetXTitle("inv mass (GeV)");
  fhV0InvMassAntiLambdaJetEvt->Sumw2();
  fListOfHistos->Add(fhV0InvMassAntiLambdaJetEvt);

  fhdROpanK0VsPt = new TH2F("hdROpanK0VsPt","V0 dR vs pt",100,0,10,100,0,10);
  fhdROpanK0VsPt->SetXTitle("opening angle R (rad)");
  fhdROpanK0VsPt->Sumw2();
  fListOfHistos->Add(fhdROpanK0VsPt);

  fhdPhiJetV0 = new TH1F("hdPhiJetV0","",640,-3.2,3.2);
  fhdPhiJetV0->SetXTitle("#Delta #phi V0 - jet");
  fhdPhiJetV0->Sumw2();
  fListOfHistos->Add(fhdPhiJetV0); 

  fhdPhiJetK0 = new TH1F("hdPhiJetK0","",640,-3.2,3.2);
  fhdPhiJetK0->SetXTitle("#Delta #phi V0 - jet");
  fhdPhiJetK0->Sumw2();
  fListOfHistos->Add(fhdPhiJetK0); 

  fhdRJetK0 = new TH1F("hdRJetK0","dN/dR K0-jet",500,0,5);
  fhdRJetK0->SetXTitle("#Delta R K0 - jet");
  fhdRJetK0->SetYTitle("1/N_{jet} dN/dR");
  fhdRJetK0->Sumw2();
  fListOfHistos->Add(fhdRJetK0); 

  fhdNdptV0 = new TH1F("hdNdptV0","dN/dpt V0",100,0,10);
  fhdNdptV0->SetXTitle("p_{T} (GeV/c)");
  fhdNdptV0->SetYTitle("1/N_{jet} dN/dpt (1/GeV)");
  fhdNdptV0->Sumw2();
  fListOfHistos->Add(fhdNdptV0);

  fhdNdptK0 = new TH1F("hdNdptK0","",100,0,10);
  fhdNdptK0->SetXTitle("p_{T} (GeV/c)");
  fhdNdptK0->SetYTitle("1/N_{jet} dN/dpt (1/GeV)");
  fhdNdptK0->Sumw2();
  fListOfHistos->Add(fhdNdptK0); 

  fhPtVsEtaK0 = new TH2F("hPtVsEtaK0","",20,-1,1,100,0,10);
  fhPtVsEtaK0->SetXTitle("#eta");  
  fhPtVsEtaK0->SetYTitle("p_{T} (GeV/c)");
  fhPtVsEtaK0->Sumw2();
  fListOfHistos->Add(fhPtVsEtaK0); 

  fhV0InvMassK0DCA = new TH1F("hV0InvMassK0DCA","",2000,0,2);
  fhV0InvMassK0DCA->SetXTitle("inv mass (GeV)");
  fhV0InvMassK0DCA->Sumw2();
  fListOfHistos->Add(fhV0InvMassK0DCA); 

  fhV0InvMassK0DCAdEdx = new TH1F("hV0InvMassK0DCAdEdx","",2000,0,2);
  fhV0InvMassK0DCAdEdx->SetXTitle("inv mass (GeV)");
  fhV0InvMassK0DCAdEdx->Sumw2();
  fListOfHistos->Add(fhV0InvMassK0DCAdEdx); 

  fhV0InvMassK0DCAPID = new TH1F("hV0InvMassK0DCAPID","",2000,0,2);
  fhV0InvMassK0DCAPID->SetXTitle("inv mass (GeV)");
  fhV0InvMassK0DCAPID->Sumw2();
  fListOfHistos->Add(fhV0InvMassK0DCAPID); 

  fhV0InvMassLambdaDCAdEdx = new TH1F("hV0InvMassLambdaDCAdEdx","",2000,0,2);
  fhV0InvMassLambdaDCAdEdx->SetXTitle("inv mass (GeV)");
  fhV0InvMassLambdaDCAdEdx->Sumw2();
  fListOfHistos->Add(fhV0InvMassLambdaDCAdEdx); 

  fhV0InvMassAntiLambdaDCAdEdx = new TH1F("hV0InvMassAntiLambdaDCAdEdx","",2000,0,2);
  fhV0InvMassAntiLambdaDCAdEdx->SetXTitle("inv mass (GeV)");
  fhV0InvMassAntiLambdaDCAdEdx->Sumw2();
  fListOfHistos->Add(fhV0InvMassAntiLambdaDCAdEdx); 

  fhdNdptK0DCA = new TH1F("hdNdptK0DCA","",100,0,10);
  fhdNdptK0DCA->SetXTitle("p_{T} (GeV/c)");
  fhdNdptK0DCA->SetYTitle("1/N_{jet} dN/dpt (1/GeV)");
  fhdNdptK0DCA->Sumw2();
  fListOfHistos->Add(fhdNdptK0DCA); 

  fhdNdptK0DCAdEdx = new TH1F("hdNdptK0DCAdEdx","",100,0,10);
  fhdNdptK0DCAdEdx->SetXTitle("p_{T} (GeV/c)");
  fhdNdptK0DCAdEdx->SetYTitle("1/N_{jet} dN/dpt (1/GeV)");
  fhdNdptK0DCAdEdx->Sumw2();
  fListOfHistos->Add(fhdNdptK0DCAdEdx); 

  fhV0InvMassK0Min = new TH1F("hV0InvMassK0Min","",2000,0,2);
  fhV0InvMassK0Min->SetXTitle("inv mass (GeV)");
  fhV0InvMassK0Min->Sumw2();
  fListOfHistos->Add(fhV0InvMassK0Min); 

  fhV0InvMassLambdaMin = new TH1F("hV0InvMassLambdaMin","",2000,0,2);
  fhV0InvMassLambdaMin->SetXTitle("inv mass (GeV)");
  fhV0InvMassLambdaMin->Sumw2();
  fListOfHistos->Add(fhV0InvMassLambdaMin); 

  fhV0InvMassAntiLambdaMin = new TH1F("hV0InvMassAntiLambdaMin","",2000,0,2);
  fhV0InvMassAntiLambdaMin->SetXTitle("inv mass (GeV)");
  fhV0InvMassAntiLambdaMin->Sumw2();
  fListOfHistos->Add(fhV0InvMassAntiLambdaMin);

  fhV0InvMassK0Max = new TH1F("hV0InvMassK0Max","",2000,0,2);
  fhV0InvMassK0Max->SetXTitle("inv mass (GeV)");
  fhV0InvMassK0Max->Sumw2();
  fListOfHistos->Add(fhV0InvMassK0Max); 

  fhV0InvMassLambdaMax = new TH1F("hV0InvMassLambdaMax","",2000,0,2);
  fhV0InvMassLambdaMax->SetXTitle("inv mass (GeV)");
  fhV0InvMassLambdaMax->Sumw2();
  fListOfHistos->Add(fhV0InvMassLambdaMax); 

  fhV0InvMassAntiLambdaMax = new TH1F("hV0InvMassAntiLambdaMax","",2000,0,2);
  fhV0InvMassAntiLambdaMax->SetXTitle("inv mass (GeV)");
  fhV0InvMassAntiLambdaMax->Sumw2();
  fListOfHistos->Add(fhV0InvMassAntiLambdaMax); 

  fhV0InvMassK0Jet = new TH1F("hV0InvMassK0Jet","",2000,0,2);
  fhV0InvMassK0Jet->SetXTitle("inv mass (GeV)");
  fhV0InvMassK0Jet->Sumw2();
  fListOfHistos->Add(fhV0InvMassK0Jet); 

  fhV0InvMassLambdaJet = new TH1F("hV0InvMassLambdaJet","",2000,0,2);
  fhV0InvMassLambdaJet->SetXTitle("inv mass (GeV)");
  fhV0InvMassLambdaJet->Sumw2();
  fListOfHistos->Add(fhV0InvMassLambdaJet); 

  fhV0InvMassAntiLambdaJet = new TH1F("hV0InvMassAntiLambdaJet","",2000,0,2);
  fhV0InvMassAntiLambdaJet->SetXTitle("inv mass (GeV)");
  fhV0InvMassAntiLambdaJet->Sumw2();
  fListOfHistos->Add(fhV0InvMassAntiLambdaJet); 

  fhV0InvMassK0Lambda = new TH1F("hV0InvMassK0Lambda","",2000,0,2);
  fhV0InvMassK0Lambda->SetXTitle("inv mass (GeV)");
  fhV0InvMassK0Lambda->Sumw2();
  fListOfHistos->Add(fhV0InvMassK0Lambda); 

  fhdNdptK0JetEvt = new TH1F("hdNdptK0JetEvt","",100,0,10);
  fhdNdptK0JetEvt->SetXTitle("p_{T} (GeV/c)");
  fhdNdptK0JetEvt->SetYTitle("1/N_{jet} dN/dpt (1/GeV)");
  fhdNdptK0JetEvt->Sumw2();
  fListOfHistos->Add(fhdNdptK0JetEvt); 
  
  fhdNdzK0 = new TH1F("hdNdzK0","",150,0,1.5);
  fhdNdzK0->SetXTitle("z");
  fhdNdzK0->SetYTitle("1/N_{jet} dN/dz");
  fhdNdzK0->Sumw2();
  fListOfHistos->Add(fhdNdzK0);

  fhdNdzK05to10 = new TH1F("hdNdzK05to10","",150,0,1.5);
  fhdNdzK05to10->SetXTitle("z");
  fhdNdzK05to10->SetYTitle("1/N_{jet} dN/dz");
  fhdNdzK05to10->Sumw2();
  fListOfHistos->Add(fhdNdzK05to10);

  fhdNdzK010to20 = new TH1F("hdNdzK010to20","",150,0,1.5);
  fhdNdzK010to20->SetXTitle("z");
  fhdNdzK010to20->SetYTitle("1/N_{jet} dN/dz");
  fhdNdzK010to20->Sumw2();
  fListOfHistos->Add(fhdNdzK010to20);

  fhdNdzK020to30 = new TH1F("hdNdzK020to30","",150,0,1.5);
  fhdNdzK020to30->SetXTitle("z");
  fhdNdzK020to30->SetYTitle("1/N_{jet} dN/dz");
  fhdNdzK020to30->Sumw2();
  fListOfHistos->Add(fhdNdzK020to30);

  fhdNdzK030to40 = new TH1F("hdNdzK030to40","",150,0,1.5);
  fhdNdzK030to40->SetXTitle("z");
  fhdNdzK030to40->SetYTitle("1/N_{jet} dN/dz");
  fhdNdzK030to40->Sumw2();
  fListOfHistos->Add(fhdNdzK030to40);

  fhdNdzK040to60 = new TH1F("hdNdzK040to60","",150,0,1.5);
  fhdNdzK040to60->SetXTitle("z");
  fhdNdzK040to60->SetYTitle("1/N_{jet} dN/dz");
  fhdNdzK040to60->Sumw2();
  fListOfHistos->Add(fhdNdzK040to60);

  fhdNdxiK0 = new TH1F("hdNdxiK0","",100,0,10);
  fhdNdxiK0->SetXTitle("xi");
  fhdNdxiK0->SetYTitle("1/N_{jet} dN/dxi");
  fhdNdxiK0->Sumw2();
  fListOfHistos->Add(fhdNdxiK0);
  
  fhdNdzLambda = new TH1F("hdNdzLambda","",150,0,1.5);
  fhdNdzLambda->SetXTitle("z");
  fhdNdzLambda->SetYTitle("1/N_{jet} dN/dz");
  fhdNdzLambda->Sumw2();
  fListOfHistos->Add(fhdNdzLambda);

  fhdNdzAntiLambda = new TH1F("hdNdzAntiLambda","",150,0,1.5);
  fhdNdzAntiLambda->SetXTitle("z");
  fhdNdzAntiLambda->SetYTitle("1/N_{jet} dN/dz");
  fhdNdzAntiLambda->Sumw2();
  fListOfHistos->Add(fhdNdzAntiLambda);

  fhdNdzK0Max = new TH1F("hdNdzK0Max","",150,0,1.5);
  fhdNdzK0Max->SetXTitle("z");
  fhdNdzK0Max->SetYTitle("1/N_{jet} dN/dz");
  fhdNdzK0Max->Sumw2();
  fListOfHistos->Add(fhdNdzK0Max);

  fhdNdxiK0Max = new TH1F("hdNdxiK0Max","",100,0,10);
  fhdNdxiK0Max->SetXTitle("xi");
  fhdNdxiK0Max->SetYTitle("1/N_{jet} dN/dxi");
  fhdNdxiK0Max->Sumw2();
  fListOfHistos->Add(fhdNdxiK0Max);

  fhdNdzLambdaMax = new TH1F("hdNdzLambdaMax","",150,0,1.5);
  fhdNdzLambdaMax->SetXTitle("z");
  fhdNdzLambdaMax->SetYTitle("1/N_{jet} dN/dz");
  fhdNdzLambdaMax->Sumw2();
  fListOfHistos->Add(fhdNdzLambdaMax); 

  fhdNdxiLambdaMax = new TH1F("hdNdxiLambdaMax","",700,0,7);
  fhdNdxiLambdaMax->SetXTitle("xi");
  fhdNdxiLambdaMax->SetYTitle("1/N_{jet} dN/dxi");
  fhdNdxiLambdaMax->Sumw2();
  fListOfHistos->Add(fhdNdxiLambdaMax); 

  fhdNdptK0Max = new TH1F("hdNdptK0Max","",100,0,10);
  fhdNdptK0Max->SetXTitle("p_{T} (GeV/c)");
  fhdNdptK0Max->SetYTitle("1/N_{jet} dN/dpt (1/GeV)");
  fhdNdptK0Max->Sumw2();
  fListOfHistos->Add(fhdNdptK0Max); 

  fhdNdptLambdaMax = new TH1F("hdNdptLambdaMax","",100,0,10);
  fhdNdptLambdaMax->SetXTitle("p_{T} (GeV/c)");
  fhdNdptLambdaMax->SetYTitle("1/N_{jet} dN/dpt (1/GeV)");
  fhdNdptLambdaMax->Sumw2();
  fListOfHistos->Add(fhdNdptLambdaMax); 

  fhdNdzK0Min = new TH1F("hdNdzK0Min","",150,0,1.5);
  fhdNdzK0Min->SetXTitle("z");
  fhdNdzK0Min->SetYTitle("1/N_{jet} dN/dz");
  fhdNdzK0Min->Sumw2();
  fListOfHistos->Add(fhdNdzK0Min); 

  fhdNdxiK0Min = new TH1F("hdNdxiK0Min","",100,0,10);
  fhdNdxiK0Min->SetXTitle("xi");
  fhdNdxiK0Min->SetYTitle("1/N_{jet} dN/dxi");
  fhdNdxiK0Min->Sumw2();
  fListOfHistos->Add(fhdNdxiK0Min); 

  fhdNdzLambdaMin = new TH1F("hdNdzLambdaMin","",150,0,1.5);
  fhdNdzLambdaMin->SetXTitle("z");
  fhdNdzLambdaMin->SetYTitle("1/N_{jet} dN/dz");
  fhdNdzLambdaMin->Sumw2();
  fListOfHistos->Add(fhdNdzLambdaMin); 

  fhdNdxiLambdaMin = new TH1F("hdNdxiLambdaMin","",700,0,7);
  fhdNdxiLambdaMin->SetXTitle("xi");
  fhdNdxiLambdaMin->SetYTitle("1/N_{jet} dN/dxi");
  fhdNdxiLambdaMin->Sumw2();
  fListOfHistos->Add(fhdNdxiLambdaMin); 

  fhdNdptK0Min = new TH1F("hdNdptK0Min","",100,0,10);
  fhdNdptK0Min->SetXTitle("p_{T} (GeV/c)");
  fhdNdptK0Min->SetYTitle("1/N_{jet} dN/dpt (1/GeV)");
  fhdNdptK0Min->Sumw2();
  fListOfHistos->Add(fhdNdptK0Min); 

  fhdNdptLambdaMin = new TH1F("hdNdptLambdaMin","",100,0,10);
  fhdNdptLambdaMin->SetXTitle("p_{T} (GeV/c)");
  fhdNdptLambdaMin->SetYTitle("1/N_{jet} dN/dpt (1/GeV)");
  fhdNdptLambdaMin->Sumw2();
  fListOfHistos->Add(fhdNdptLambdaMin); 

  fhdNdzK0Jet = new TH1F("hdNdzK0Jet","",150,0,1.5);
  fhdNdzK0Jet->SetXTitle("z");
  fhdNdzK0Jet->SetYTitle("1/N_{jet} dN/dz");
  fhdNdzK0Jet->Sumw2();
  fListOfHistos->Add(fhdNdzK0Jet);

  fhdNdxiK0Jet = new TH1F("hdNdxiK0Jet","",100,0,10);
  fhdNdxiK0Jet->SetXTitle("xi");
  fhdNdxiK0Jet->SetYTitle("1/N_{jet} dN/dxi");
  fhdNdxiK0Jet->Sumw2();
  fListOfHistos->Add(fhdNdxiK0Jet); 

  fhdNdzLambdaJet = new TH1F("hdNdzLambdaJet","",150,0,1.5);
  fhdNdzLambdaJet->SetXTitle("z");
  fhdNdzLambdaJet->SetYTitle("1/N_{jet} dN/dz");
  fhdNdzLambdaJet->Sumw2();
  fListOfHistos->Add(fhdNdzLambdaJet); 

  fhdNdxiLambdaJet = new TH1F("hdNdxiLambdaJet","",700,0,7);
  fhdNdxiLambdaJet->SetXTitle("xi");
  fhdNdxiLambdaJet->SetYTitle("1/N_{jet} dN/dxi");
  fhdNdxiLambdaJet->Sumw2();
  fListOfHistos->Add(fhdNdxiLambdaJet); 

  fhdNdptK0Jet = new TH1F("hdNdptK0Jet","",100,0,10);
  fhdNdptK0Jet->SetXTitle("p_{T} (GeV/c)");
  fhdNdptK0Jet->SetYTitle("1/N_{jet} dN/dpt (1/GeV)");
  fhdNdptK0Jet->Sumw2();
  fListOfHistos->Add(fhdNdptK0Jet); 

  fhdNdptLambdaJet = new TH1F("hdNdptLambdaJet","",100,0,10);
  fhdNdptLambdaJet->SetXTitle("p_{T} (GeV/c)");
  fhdNdptLambdaJet->SetYTitle("1/N_{jet} dN/dpt (1/GeV)");
  fhdNdptLambdaJet->Sumw2();
  fListOfHistos->Add(fhdNdptLambdaJet); 

  fhdEdxVsMomV0 = new TH2F("hdEdxVsMomV0","",200,-10,10,200,0,2000);
  fhdEdxVsMomV0->SetXTitle("mom (GeV/c)");
  fhdEdxVsMomV0->SetYTitle("dE/dx");
  fhdEdxVsMomV0->Sumw2();
  fListOfHistos->Add(fhdEdxVsMomV0); 

  fhdEdxVsMomV0pidEdx = new TH2F("hdEdxVsMomV0pidEdx","",200,-10,10,200,0,2000);
  fhdEdxVsMomV0pidEdx->SetXTitle("mom (GeV/c)");
  fhdEdxVsMomV0pidEdx->SetYTitle("dE/dx");
  fhdEdxVsMomV0pidEdx->Sumw2();
  fListOfHistos->Add(fhdEdxVsMomV0pidEdx); 

  fhdEdxVsMomV0piPID = new TH2F("hdEdxVsMomV0piPID","",200,-10,10,200,0,2000);
  fhdEdxVsMomV0piPID->SetXTitle("mom (GeV/c)");
  fhdEdxVsMomV0piPID->SetYTitle("dE/dx");
  fhdEdxVsMomV0piPID->Sumw2();
  fListOfHistos->Add(fhdEdxVsMomV0piPID); 

  fhdPhiJetK0MC = new TH1F("hdPhiJetK0MC","",640,-3.2,3.2);
  fhdPhiJetK0MC->SetXTitle("#Delta #phi K0 - jet");
  fhdPhiJetK0MC->SetYTitle("1/N_{jet} dN/dphi");
  fhdPhiJetK0MC->Sumw2();
  fListOfHistos->Add(fhdPhiJetK0MC); 

  fhdRJetK0MC   = new TH1F("hdRJetK0MC","dN/R K0-jet",500,0,5);
  fhdRJetK0MC->SetXTitle("#Delta R K0 - jet");
  fhdRJetK0MC->SetYTitle("1/N_{jet} dN/dR");
  fhdRJetK0MC->Sumw2();
  fListOfHistos->Add(fhdRJetK0MC); 

  fhdRV0MC =  new TH1F("hdRV0MC","",500,0.,1.);
  fhdRV0MC->SetXTitle("#Delta R");
  fhdRV0MC->SetYTitle("");
  fhdRV0MC->Sumw2();
  fListOfHistos->Add(fhdRV0MC); 

  fhdNdptchPiMCMax = new TH1F("hdNdptchPiMCMax","",100,0,10);
  fhdNdptchPiMCMax->SetXTitle("p_{T} (GeV/c)");
  fhdNdptchPiMCMax->SetYTitle("1/N_{jet} dN/dpt (1/GeV)");
  fhdNdptchPiMCMax->Sumw2();
  fListOfHistos->Add(fhdNdptchPiMCMax); 

  fhdNdptK0MCMax = new TH1F("hdNdptK0MCMax","",100,0,10);
  fhdNdptK0MCMax->SetXTitle("p_{T} (GeV/c)");
  fhdNdptK0MCMax->SetYTitle("1/N_{jet} dN/dpt (1/GeV)");
  fhdNdptK0MCMax->Sumw2();
  fListOfHistos->Add(fhdNdptK0MCMax); 

  fhdNdptchKMCMax = new TH1F("hdNdptchKMCMax","",100,0,10);
  fhdNdptchKMCMax->SetXTitle("p_{T} (GeV/c)");
  fhdNdptchKMCMax->SetYTitle("1/N_{jet} dN/dpt (1/GeV)");
  fhdNdptchKMCMax->Sumw2();
  fListOfHistos->Add(fhdNdptchKMCMax); 

  fhdNdptpMCMax = new TH1F("hdNdptpMCMax","",100,0,10);
  fhdNdptpMCMax->SetXTitle("p_{T} (GeV/c)");
  fhdNdptpMCMax->SetYTitle("1/N_{jet} dN/dpt (1/GeV)");
  fhdNdptpMCMax->Sumw2();
  fListOfHistos->Add(fhdNdptpMCMax); 

  fhdNdptpBarMCMax = new TH1F("hdNdptpBarMCMax","",100,0,10);
  fhdNdptpBarMCMax->SetXTitle("p_{T} (GeV/c)");
  fhdNdptpBarMCMax->SetYTitle("1/N_{jet} dN/dpt (1/GeV)");
  fhdNdptpBarMCMax->Sumw2();
  fListOfHistos->Add(fhdNdptpBarMCMax); 

  fhdNdptLambdaMCMax = new TH1F("hdNdptLambdaMCMax","",100,0,10);
  fhdNdptLambdaMCMax->SetXTitle("p_{T} (GeV/c)");
  fhdNdptLambdaMCMax->SetYTitle("1/N_{jet} dN/dpt (1/GeV)");
  fhdNdptLambdaMCMax->Sumw2();
  fListOfHistos->Add(fhdNdptLambdaMCMax); 

  fhdNdptLambdaBarMCMax = new TH1F("hdNdptLambdaBarMCMax","",100,0,10);
  fhdNdptLambdaBarMCMax->SetXTitle("p_{T} (GeV/c)");
  fhdNdptLambdaBarMCMax->SetYTitle("1/N_{jet} dN/dpt (1/GeV)");
  fhdNdptLambdaBarMCMax->Sumw2();
  fListOfHistos->Add(fhdNdptLambdaBarMCMax); 


  fhdNdptchPiMCMin = new TH1F("hdNdptchPiMCMin","",100,0,10);
  fhdNdptchPiMCMin->SetXTitle("p_{T} (GeV/c)");
  fhdNdptchPiMCMin->SetYTitle("1/N_{jet} dN/dpt (1/GeV)");
  fhdNdptchPiMCMin->Sumw2();
  fListOfHistos->Add(fhdNdptchPiMCMin); 
 
  fhdNdptK0MCMin = new TH1F("hdNdptK0MCMin","",100,0,10);
  fhdNdptK0MCMin->SetXTitle("p_{T} (GeV/c)");
  fhdNdptK0MCMin->SetYTitle("1/N_{jet} dN/dpt (1/GeV)");
  fhdNdptK0MCMin->Sumw2();
  fListOfHistos->Add(fhdNdptK0MCMin); 

  fhdNdptchKMCMin = new TH1F("hdNdptchKMCMin","",100,0,10);
  fhdNdptchKMCMin->SetXTitle("p_{T} (GeV/c)");
  fhdNdptchKMCMin->SetYTitle("1/N_{jet} dN/dpt (1/GeV)");
  fhdNdptchKMCMin->Sumw2();
  fListOfHistos->Add(fhdNdptchKMCMin); 

  fhdNdptpMCMin = new TH1F("hdNdptpMCMin","",100,0,10);
  fhdNdptpMCMin->SetXTitle("p_{T} (GeV/c)");
  fhdNdptpMCMin->SetYTitle("1/N_{jet} dN/dpt (1/GeV)");
  fhdNdptpMCMin->Sumw2();
  fListOfHistos->Add(fhdNdptpMCMin); 

  fhdNdptpBarMCMin = new TH1F("hdNdptpBarMCMin","",100,0,10);
  fhdNdptpBarMCMin->SetXTitle("p_{T} (GeV/c)");
  fhdNdptpBarMCMin->SetYTitle("1/N_{jet} dN/dpt (1/GeV)");
  fhdNdptpBarMCMin->Sumw2();
  fListOfHistos->Add(fhdNdptpBarMCMin); 

  fhdNdptLambdaMCMin = new TH1F("hdNdptLambdaMCMin","",100,0,10);
  fhdNdptLambdaMCMin->SetXTitle("p_{T} (GeV/c)");
  fhdNdptLambdaMCMin->SetYTitle("1/N_{jet} dN/dpt (1/GeV)");
  fhdNdptLambdaMCMin->Sumw2();
  fListOfHistos->Add(fhdNdptLambdaMCMin); 

  fhdNdptLambdaBarMCMin = new TH1F("hdNdptLambdaBarMCMin","",100,0,10);
  fhdNdptLambdaBarMCMin->SetXTitle("p_{T} (GeV/c)");
  fhdNdptLambdaBarMCMin->SetYTitle("1/N_{jet} dN/dpt (1/GeV)");
  fhdNdptLambdaBarMCMin->Sumw2();
  fListOfHistos->Add(fhdNdptLambdaBarMCMin); 

  fhdNdptOmegaMCMin = new TH1F("hdNdptOmegaMCMin","",100,0,10);;
  fhdNdptOmegaMCMin->SetXTitle("p_{T} (GeV/c)");
  fhdNdptOmegaMCMin->SetYTitle("1/N_{event} dN/dpt (1/GeV)");
  fhdNdptOmegaMCMin->Sumw2();
  fListOfHistos->Add(fhdNdptOmegaMCMin); 

  fhdNdptOmegaBarMCMin = new TH1F("hdNdptOmegaBarMCMin","",100,0,10);;
  fhdNdptOmegaBarMCMin->SetXTitle("p_{T} (GeV/c)");
  fhdNdptOmegaBarMCMin->SetYTitle("1/N_{event} dN/dpt (1/GeV)");
  fhdNdptOmegaBarMCMin->Sumw2();
  fListOfHistos->Add(fhdNdptOmegaBarMCMin); 

  fhdNdptchPiMCJet = new TH1F("hdNdptchPiMCJet","",100,0,10);
  fhdNdptchPiMCJet->SetXTitle("p_{T} (GeV/c)");
  fhdNdptchPiMCJet->SetYTitle("1/N_{jet} dN/dpt (1/GeV)");
  fhdNdptchPiMCJet->Sumw2();
  fListOfHistos->Add(fhdNdptchPiMCJet); 

  fhdNdptK0MCJet = new TH1F("hdNdptK0MCJet","",100,0,10);
  fhdNdptK0MCJet->SetXTitle("p_{T} (GeV/c)");
  fhdNdptK0MCJet->SetYTitle("1/N_{jet} dN/dpt (1/GeV)");
  fhdNdptK0MCJet->Sumw2();
  fListOfHistos->Add(fhdNdptK0MCJet); 

  fhdNdptchKMCJet = new TH1F("hdNdptchKMCJet","",100,0,10);
  fhdNdptchKMCJet->SetXTitle("p_{T} (GeV/c)");
  fhdNdptchKMCJet->SetYTitle("1/N_{jet} dN/dpt (1/GeV)");
  fhdNdptchKMCJet->Sumw2();
  fListOfHistos->Add(fhdNdptchKMCJet);

  fhdNdptpMCJet = new TH1F("hdNdptpMCJet","",100,0,10);
  fhdNdptpMCJet->SetXTitle("p_{T} (GeV/c)");
  fhdNdptpMCJet->SetYTitle("1/N_{jet} dN/dpt (1/GeV)");
  fhdNdptpMCJet->Sumw2();
  fListOfHistos->Add(fhdNdptpMCJet); 

  fhdNdptpBarMCJet = new TH1F("hdNdptpBarMCJet","",100,0,10);
  fhdNdptpBarMCJet->SetXTitle("p_{T} (GeV/c)");
  fhdNdptpBarMCJet->SetYTitle("1/N_{jet} dN/dpt (1/GeV)");
  fhdNdptpBarMCJet->Sumw2();
  fListOfHistos->Add(fhdNdptpBarMCJet);

  fhdNdptLambdaMCJet = new TH1F("hdNdptLambdaMCJet","",100,0,10);
  fhdNdptLambdaMCJet->SetXTitle("p_{T} (GeV/c)");
  fhdNdptLambdaMCJet->SetYTitle("1/N_{jet} dN/dpt (1/GeV)");
  fhdNdptLambdaMCJet->Sumw2();
  fListOfHistos->Add(fhdNdptLambdaMCJet); 

  fhdNdptLambdaBarMCJet = new TH1F("hdNdptLambdaBarMCJet","",100,0,10);
  fhdNdptLambdaBarMCJet->SetXTitle("p_{T} (GeV/c)");
  fhdNdptLambdaBarMCJet->SetYTitle("1/N_{jet} dN/dpt (1/GeV)");
  fhdNdptLambdaBarMCJet->Sumw2();
  fListOfHistos->Add(fhdNdptLambdaBarMCJet); 

  fhPIDMC = CreatePIDhisto("hPIDMC");
  fhPIDMC->Sumw2();
  fListOfHistos->Add(fhPIDMC); 

  fhPIDMC_quarkEv = CreatePIDhisto("hPIDMC_quarkEv");
  fhPIDMC_quarkEv->Sumw2();
  fListOfHistos->Add(fhPIDMC_quarkEv); 

  fhPIDMC_gluonEv = CreatePIDhisto("hPIDMC_gluonEv");
  fhPIDMC_gluonEv->Sumw2();
  fListOfHistos->Add(fhPIDMC_gluonEv); 

  fhPIDMCAll = CreatePIDhisto("hPIDMCAll");
  fhPIDMCAll->Sumw2();
  fListOfHistos->Add(fhPIDMCAll); 

  fhPIDMCMin = CreatePIDhisto("hPIDMCMin");
  fhPIDMCMin->Sumw2();
  fListOfHistos->Add(fhPIDMCMin); 

  fhPIDMCJet = CreatePIDhisto("hPIDMCJet");
  fhPIDMCJet->Sumw2();
  fListOfHistos->Add(fhPIDMCJet); 

  fhPIDMCMotherK0 = CreatePIDhisto("hPIDMCMotherK0"); 
  fhPIDMCMotherK0->Sumw2();
  fListOfHistos->Add(fhPIDMCMotherK0); 

  fhPIDMCGrandMotherK0 = CreatePIDhisto("hPIDMCGrandMotherK0"); 
  fhPIDMCGrandMotherK0->Sumw2();
  fListOfHistos->Add(fhPIDMCGrandMotherK0); 

  fhPIDMCMotherChK = CreatePIDhisto("hPIDMCMotherChK"); 
  fhPIDMCMotherChK->Sumw2();
  fListOfHistos->Add(fhPIDMCMotherChK); 

  fhPIDMCMotherK0Trans = CreatePIDhisto("hPIDMCMotherK0Trans"); 
  fhPIDMCMotherK0Trans->Sumw2();
  fListOfHistos->Add(fhPIDMCMotherK0Trans); 

  fhPIDMCGrandMotherK0Trans = CreatePIDhisto("hPIDMCGrandMotherK0Trans"); 
  fhPIDMCGrandMotherK0Trans->Sumw2();
  fListOfHistos->Add(fhPIDMCGrandMotherK0Trans); 

  fhPIDMCMotherChKTrans = CreatePIDhisto("hPIDMCMotherChKTrans"); 
  fhPIDMCMotherChKTrans->Sumw2();
  fListOfHistos->Add(fhPIDMCMotherChKTrans); 

  fhdNdptgammaMC = new TH1F("hdNdptgammaMC","",100,0,10);
  fhdNdptgammaMC->SetXTitle("p_{T} (GeV/c)");
  fhdNdptgammaMC->SetYTitle("1/N_{event} dN/dpt (1/GeV)");
  fhdNdptgammaMC->Sumw2();
  fListOfHistos->Add(fhdNdptgammaMC); 

  fhdNdptchPiMC  = new TH1F("hdNdptchPiMC","",100,0,10);;
  fhdNdptchPiMC->SetXTitle("p_{T} (GeV/c)");
  fhdNdptchPiMC->SetYTitle("1/N_{event} dN/dpt (1/GeV)");
  fhdNdptchPiMC->Sumw2();
  fListOfHistos->Add(fhdNdptchPiMC); 

  fhdNdptpi0MC    = new TH1F("hdNdptpi0MC","",100,0,10);;
  fhdNdptpi0MC->SetXTitle("p_{T} (GeV/c)");
  fhdNdptpi0MC->SetYTitle("1/N_{event} dN/dpt (1/GeV)");
  fhdNdptpi0MC->Sumw2();
  fListOfHistos->Add(fhdNdptpi0MC); 

  fhdNdptK0MC = new TH1F("hdNdptK0MC","",100,0,10);;
  fhdNdptK0MC->SetXTitle("p_{T} (GeV/c)");
  fhdNdptK0MC->SetYTitle("1/N_{event} dN/dpt (1/GeV)");
  fhdNdptK0MC->Sumw2();
  fListOfHistos->Add(fhdNdptK0MC); 

  fhdNdptchKMC    = new TH1F("hdNdptchKMC","",100,0,10);;
  fhdNdptchKMC->SetXTitle("p_{T} (GeV/c)");
  fhdNdptchKMC->SetYTitle("1/N_{event} dN/dpt (1/GeV)");
  fhdNdptchKMC->Sumw2();
  fListOfHistos->Add(fhdNdptchKMC); 

  fhdNdptpMC    = new TH1F("hdNdptpMC","",100,0,10);;
  fhdNdptpMC->SetXTitle("p_{T} (GeV/c)");
  fhdNdptpMC->SetYTitle("1/N_{event} dN/dpt (1/GeV)");
  fhdNdptpMC->Sumw2();
  fListOfHistos->Add(fhdNdptpMC); 

  fhdNdptpBarMC    = new TH1F("hdNdptpBarMC","",100,0,10);;
  fhdNdptpBarMC->SetXTitle("p_{T} (GeV/c)");
  fhdNdptpBarMC->SetYTitle("1/N_{event} dN/dpt (1/GeV)");
  fhdNdptpBarMC->Sumw2();
  fListOfHistos->Add(fhdNdptpBarMC); 

  fhdNdptLambdaMC    = new TH1F("hdNdptLambdaMC","",100,0,10);;
  fhdNdptLambdaMC->SetXTitle("p_{T} (GeV/c)");
  fhdNdptLambdaMC->SetYTitle("1/N_{event} dN/dpt (1/GeV)");
  fhdNdptLambdaMC->Sumw2();
  fListOfHistos->Add(fhdNdptLambdaMC); 

  fhdNdptLambdaBarMC    = new TH1F("hdNdptLambdaBarMC","",100,0,10);;
  fhdNdptLambdaBarMC->SetXTitle("p_{T} (GeV/c)");
  fhdNdptLambdaBarMC->SetYTitle("1/N_{event} dN/dpt (1/GeV)");
  fhdNdptLambdaBarMC->Sumw2();
  fListOfHistos->Add(fhdNdptLambdaBarMC); 

  fhdNdptOmegaMC = new TH1F("hdNdptOmegaMC","",100,0,10);;
  fhdNdptOmegaMC->SetXTitle("p_{T} (GeV/c)");
  fhdNdptOmegaMC->SetYTitle("1/N_{event} dN/dpt (1/GeV)");
  fhdNdptOmegaMC->Sumw2();
  fListOfHistos->Add(fhdNdptOmegaMC); 

  fhdNdptOmegaBarMC = new TH1F("hdNdptOmegaBarMC","",100,0,10);;
  fhdNdptOmegaBarMC->SetXTitle("p_{T} (GeV/c)");
  fhdNdptOmegaBarMC->SetYTitle("1/N_{event} dN/dpt (1/GeV)");
  fhdNdptOmegaBarMC->Sumw2();
  fListOfHistos->Add(fhdNdptOmegaBarMC); 

  fhdNdxiMC = new TH1F("hdNdxiMC","",100,0,10); 
  fhdNdxiMC->SetXTitle("#xi");
  fhdNdxiMC->Sumw2();
  fListOfHistos->Add(fhdNdxiMC);            

  fhdNdxiK0MC = new TH1F("hdNdxiK0MC","",100,0,10); 
  fhdNdxiK0MC->SetXTitle("#xi");
  fhdNdxiK0MC->Sumw2();
  fListOfHistos->Add(fhdNdxiK0MC);            

  fhdNdxiK0MCJet = new TH1F("hdNdxiK0MCJet","",100,0,10); 
  fhdNdxiK0MCJet->SetXTitle("#xi");
  fhdNdxiK0MCJet->Sumw2();
  fListOfHistos->Add(fhdNdxiK0MCJet);

  fhdNdzK0MC = new TH1F("hdNdzK0MC","",150,0,1.5);
  fhdNdzK0MC->SetXTitle("z");
  fhdNdzK0MC->SetYTitle("1/N_{jet} dN/dz");
  fhdNdzK0MC->Sumw2();
  fListOfHistos->Add(fhdNdzK0MC);

  fhdNdzK0MCJet = new TH1F("hdNdzK0MCJet","",150,0,1.5);
  fhdNdzK0MCJet->SetXTitle("z");
  fhdNdzK0MCJet->SetYTitle("1/N_{jet} dN/dz");
  fhdNdzK0MCJet->Sumw2();
  fListOfHistos->Add(fhdNdzK0MCJet);

  fhdNdptK0MCJetEvt = new TH1F("hdNdptK0MCJetEvt","",100,0,10);;
  fhdNdptK0MCJetEvt->SetXTitle("p_{T} (GeV/c)");
  fhdNdptK0MCJetEvt->SetYTitle("1/N_{event} dN/dpt (1/GeV)");
  fhdNdptK0MCJetEvt->Sumw2();
  fListOfHistos->Add(fhdNdptK0MCJetEvt);

  fhnJetsAODvsMC = new TH2F("hnJetsAODvsMC","",20,0,20,20,0.,20.); 
  fhnJetsAODvsMC->SetXTitle("n jets MC");
  fhnJetsAODvsMC->SetXTitle("n jets AOD");
  fhnJetsAODvsMC->Sumw2();
  fListOfHistos->Add(fhnJetsAODvsMC);

  fhLeadingPtAODvsMC = new TH2F("hLeadingPtAODvsMC","",20,0,20,20,0.,20.); 
  fhLeadingPtAODvsMC->SetXTitle("p_{T} MC (GeV/c)");
  fhLeadingPtAODvsMC->SetYTitle("p_{T} AOD (GeV/c)");
  fhLeadingPtAODvsMC->Sumw2();
  fListOfHistos->Add(fhLeadingPtAODvsMC);

  fhLeadingEtaAODvsMC = new TH2F("hLeadingEtaAODvsMC","",20,-1,1,20,-1.,1.); 
  fhLeadingEtaAODvsMC->SetXTitle("#eta MC");
  fhLeadingEtaAODvsMC->SetYTitle("#eta AOD");
  fhLeadingEtaAODvsMC->Sumw2();
  fListOfHistos->Add(fhLeadingEtaAODvsMC);

  fhLeadingPhiAODvsMC = new TH2F("hLeadingPhiAODvsMC","",63,0,6.3,63,0.,6.3); 
  fhLeadingPhiAODvsMC->SetXTitle("#phi MC");
  fhLeadingPhiAODvsMC->SetYTitle("#phi AOD");
  fhLeadingPhiAODvsMC->Sumw2();
  fListOfHistos->Add(fhLeadingPhiAODvsMC);


  fhnTracksLeadingAODvsMC = new TH2F("hnTracksLeadingAODvsMC","",20,0.,20,20,-0.5,19.5);
  fhnTracksLeadingAODvsMC->SetXTitle("nTracks MC");
  fhnTracksLeadingAODvsMC->SetYTitle("nTracks AOD");
  fhnTracksLeadingAODvsMC->Sumw2();
  fListOfHistos->Add(fhnTracksLeadingAODvsMC);

  fhLeadingdRAODMC = new TH1F("hLeadingdRAODMC","",40,0.,4);
  fhLeadingdRAODMC->SetXTitle("#Delta R");
  fhLeadingdRAODMC->Sumw2();
  fListOfHistos->Add(fhLeadingdRAODMC);


  fhLeadingPtAODvsMCdRcut = new TH2F("hLeadingPtAODvsMCdRcut","",40,0,20,40,0.,20.); 
  fhLeadingPtAODvsMCdRcut->SetXTitle("p_{T} MC (GeV/c)");
  fhLeadingPtAODvsMCdRcut->SetYTitle("p_{T} AOD (GeV/c)");
  fhLeadingPtAODvsMCdRcut->Sumw2();
  fListOfHistos->Add(fhLeadingPtAODvsMCdRcut);


  fhdnTracksVsdPtLeadingAODMC = new TH2F("hdnTracksVsdPtLeadingAODMC","",40,-10.,10,20,-10.5,9.5);
  fhdnTracksVsdPtLeadingAODMC->SetXTitle("#Delta Pt AOD-MC (GeV/c)");
  fhdnTracksVsdPtLeadingAODMC->SetYTitle("#Delta N AOD-MC");
  fhdnTracksVsdPtLeadingAODMC->Sumw2();
  fListOfHistos->Add(fhdnTracksVsdPtLeadingAODMC);

  fhnTracksJetVsPtAOD = new TH2F("hnTracksJetVsPtAOD","",50,0,50,20,-0.5,19.5);
  fhnTracksJetVsPtAOD->SetXTitle("p_{T} (GeV/c)");
  fhnTracksJetVsPtAOD->SetYTitle("nTracks");
  fhnTracksJetVsPtAOD->Sumw2();
  fListOfHistos->Add(fhnTracksJetVsPtAOD);

  fhnTracksJetVsPtAODquarkEv = new TH2F("hnTracksJetVsPtAODquarkEv","",50,0,50,20,-0.5,19.5);
  fhnTracksJetVsPtAODquarkEv->SetXTitle("p_{T} (GeV/c)");
  fhnTracksJetVsPtAODquarkEv->SetYTitle("nTracks");
  fhnTracksJetVsPtAODquarkEv->Sumw2();
  fListOfHistos->Add(fhnTracksJetVsPtAODquarkEv);

  fhRadiusJetVsPtAOD = new TH2F("hRadiusJetVsPtAOD","",50,0,50,10,0.,1);
  fhRadiusJetVsPtAOD->SetXTitle("p_{T} (GeV/c)");
  fhRadiusJetVsPtAOD->SetYTitle("radius");
  fhRadiusJetVsPtAOD->Sumw2();
  fListOfHistos->Add(fhRadiusJetVsPtAOD);

  fhnTracksJetVsPtMC = new TH2F("hnTracksJetVsPtMC","",50,0,50,20,-0.5,19.5);
  fhnTracksJetVsPtMC->SetXTitle("p_{T} (GeV/c)");
  fhnTracksJetVsPtMC->SetYTitle("nTracks");
  fhnTracksJetVsPtMC->Sumw2();
  fListOfHistos->Add(fhnTracksJetVsPtMC);

  fhnTracksJetVsPtMCquarkEv = new TH2F("hnTracksJetVsPtMCquarkEv","",50,0,50,20,-0.5,19.5);
  fhnTracksJetVsPtMCquarkEv->SetXTitle("p_{T} (GeV/c)");
  fhnTracksJetVsPtMCquarkEv->SetYTitle("nTracks");
  fhnTracksJetVsPtMCquarkEv->Sumw2();
  fListOfHistos->Add(fhnTracksJetVsPtMCquarkEv);

  fhRadiusJetVsPtMC = new TH2F("hRadiusJetVsPtMC","",50,0,50,10,0.,1);
  fhRadiusJetVsPtMC->SetXTitle("p_{T} (GeV/c)");
  fhRadiusJetVsPtMC->SetYTitle("radius");
  fhRadiusJetVsPtMC->Sumw2();
  fListOfHistos->Add(fhRadiusJetVsPtMC);

  fhnTracksJetVsPtMCK0 = new TH2F("hnTracksJetVsPtMCK0","",50,0,50,20,-0.5,19.5);
  fhnTracksJetVsPtMCK0->SetXTitle("p_{T} (GeV/c)");
  fhnTracksJetVsPtMCK0->SetYTitle("nTracks");
  fhnTracksJetVsPtMCK0->Sumw2();
  fListOfHistos->Add(fhnTracksJetVsPtMCK0);


  fhnTracksJetVsPtMCK0quarkEv = new TH2F("hnTracksJetVsPtMCK0quarkEv","",50,0,50,20,-0.5,19.5);
  fhnTracksJetVsPtMCK0quarkEv->SetXTitle("p_{T} (GeV/c)");
  fhnTracksJetVsPtMCK0quarkEv->SetYTitle("nTracks");
  fhnTracksJetVsPtMCK0quarkEv->Sumw2();
  fListOfHistos->Add(fhnTracksJetVsPtMCK0quarkEv);

  fhRadiusJetVsPtMCK0 = new TH2F("hRadiusJetVsPtMCK0","",50,0,50,10,0.,1);
  fhRadiusJetVsPtMCK0->SetXTitle("p_{T} (GeV/c)");
  fhRadiusJetVsPtMCK0->SetYTitle("radius");
  fhRadiusJetVsPtMCK0->Sumw2();
  fListOfHistos->Add(fhRadiusJetVsPtMCK0);

  fhnTracksJetVsPtAODK0 = new TH2F("hnTracksJetVsPtAODK0","",50,0,50,20,-0.5,19.5);
  fhnTracksJetVsPtAODK0->SetXTitle("p_{T} (GeV/c)");
  fhnTracksJetVsPtAODK0->SetYTitle("nTracks AODK0");
  fhnTracksJetVsPtAODK0->Sumw2();
  fListOfHistos->Add(fhnTracksJetVsPtAODK0);

  fhnTracksJetVsPtAODK0quarkEv = new TH2F("hnTracksJetVsPtAODK0quarkEv","",50,0,50,20,-0.5,19.5);
  fhnTracksJetVsPtAODK0quarkEv->SetXTitle("p_{T} (GeV/c)");
  fhnTracksJetVsPtAODK0quarkEv->SetYTitle("nTracks AODK0quarkEv");
  fhnTracksJetVsPtAODK0quarkEv->Sumw2();
  fListOfHistos->Add(fhnTracksJetVsPtAODK0quarkEv);

  fhRadiusJetVsPtAODK0 = new TH2F("hRadiusJetVsPtAODK0","",50,0,50,10,0.,1);
  fhRadiusJetVsPtAODK0->SetXTitle("p_{T} (GeV/c)");
  fhRadiusJetVsPtAODK0->SetYTitle("radius");
  fhRadiusJetVsPtAODK0->Sumw2();
  fListOfHistos->Add(fhRadiusJetVsPtAODK0);

  fhnTracksJetVsPtAODpKch = new TH2F("hnTracksJetVsPtAODpKch","",50,0,50,20,-0.5,19.5);
  fhnTracksJetVsPtAODpKch->SetXTitle("p_{T} (GeV/c)");
  fhnTracksJetVsPtAODpKch->SetYTitle("nTracks AODpKch");
  fhnTracksJetVsPtAODpKch->Sumw2();
  fListOfHistos->Add(fhnTracksJetVsPtAODpKch);

  fhRadiusJetVsPtAODpKch = new TH2F("hRadiusJetVsPtAODpKch","",50,0,50,20,-0.5,19.5);
  fhRadiusJetVsPtAODpKch->SetXTitle("p_{T} (GeV/c)");
  fhRadiusJetVsPtAODpKch->SetYTitle("Radius AODpKch");
  fhRadiusJetVsPtAODpKch->Sumw2();
  fListOfHistos->Add(fhRadiusJetVsPtAODpKch);

  fhPythiaProcess      = CreatePythiaIDhisto("hPythiaProcess");
  fListOfHistos->Add(fhPythiaProcess);  

  fhPythiaProcessK0   = CreatePythiaIDhisto("hPythiaProcessK0");
  fListOfHistos->Add(fhPythiaProcessK0);  

  fhPythiaProcessKch   = CreatePythiaIDhisto("hPythiaProcessKch");
  fListOfHistos->Add(fhPythiaProcessKch);  

  fhPythiaProcessp    = CreatePythiaIDhisto("hPythiaProcessp");
  fListOfHistos->Add(fhPythiaProcessp);  
  
  fhPythiaProcesspbar = CreatePythiaIDhisto("hPythiaProcesspbar");
  fListOfHistos->Add(fhPythiaProcesspbar);  
 
  fhdNdzJets5to10 = new TH1F("hdNdzJets5to10","",25,0,1.25);
  fhdNdzJets5to10->SetXTitle("z");
  fhdNdzJets5to10->SetYTitle("dN/dz");
  fhdNdzJets5to10->Sumw2();
  fListOfHistos->Add(fhdNdzJets5to10);

  fhdNdzJets10to20 = new TH1F("hdNdzJets10to20","",25,0,1.25);
  fhdNdzJets10to20->SetXTitle("z");
  fhdNdzJets10to20->SetYTitle("dN/dz");
  fhdNdzJets10to20->Sumw2();
  fListOfHistos->Add(fhdNdzJets10to20);

  fhdNdzJets20to30 = new TH1F("hdNdzJets20to30","",25,0,1.25);
  fhdNdzJets20to30->SetXTitle("z");
  fhdNdzJets20to30->SetYTitle("dN/dz");
  fhdNdzJets20to30->Sumw2();
  fListOfHistos->Add(fhdNdzJets20to30);

  fhdNdzJets30to40 = new TH1F("hdNdzJets30to40","",25,0,1.25);
  fhdNdzJets30to40->SetXTitle("z");
  fhdNdzJets30to40->SetYTitle("dN/dz");
  fhdNdzJets30to40->Sumw2();
  fListOfHistos->Add(fhdNdzJets30to40);

  fhdNdzJets40to60 = new TH1F("hdNdzJets40to60","",25,0,1.25);
  fhdNdzJets40to60->SetXTitle("z");
  fhdNdzJets40to60->SetYTitle("dN/dz");
  fhdNdzJets40to60->Sumw2();
  fListOfHistos->Add(fhdNdzJets40to60);
  

  fhdNdxiJets5to10 = new TH1F("hdNdxiJets5to10","",70,0,7);
  fhdNdxiJets5to10->SetXTitle("z");
  fhdNdxiJets5to10->SetYTitle("dN/dz");
  fhdNdxiJets5to10->Sumw2();
  fListOfHistos->Add(fhdNdxiJets5to10);

  fhdNdxiJets10to20 = new TH1F("hdNdxiJets10to20","",70,0,7);
  fhdNdxiJets10to20->SetXTitle("z");
  fhdNdxiJets10to20->SetYTitle("dN/dz");
  fhdNdxiJets10to20->Sumw2();
  fListOfHistos->Add(fhdNdxiJets10to20);

  fhdNdxiJets20to30 = new TH1F("hdNdxiJets20to30","",70,0,7);
  fhdNdxiJets20to30->SetXTitle("z");
  fhdNdxiJets20to30->SetYTitle("dN/dz");
  fhdNdxiJets20to30->Sumw2();
  fListOfHistos->Add(fhdNdxiJets20to30);

  fhdNdxiJets30to40 = new TH1F("hdNdxiJets30to40","",70,0,7);
  fhdNdxiJets30to40->SetXTitle("z");
  fhdNdxiJets30to40->SetYTitle("dN/dz");
  fhdNdxiJets30to40->Sumw2();
  fListOfHistos->Add(fhdNdxiJets30to40);

  fhdNdxiJets40to60 = new TH1F("hdNdxiJets40to60","",70,0,7);
  fhdNdxiJets40to60->SetXTitle("z");
  fhdNdxiJets40to60->SetYTitle("dN/dz");
  fhdNdxiJets40to60->Sumw2();
  fListOfHistos->Add(fhdNdxiJets40to60);

  fhdNdptTracksJetPt5to10 = new TH1F("hdNdptTracksJetPt5to10","",250,0,25);
  fhdNdptTracksJetPt5to10->SetXTitle("p_{T} (GeV)");
  fhdNdptTracksJetPt5to10->SetYTitle("dN/dp_{T} 1/GeV");
  fhdNdptTracksJetPt5to10->Sumw2();
  fListOfHistos->Add(fhdNdptTracksJetPt5to10);

  fhdNdptTracksJetPt10to20 = new TH1F("hdNdptTracksJetPt10to20","",25,0,25);
  fhdNdptTracksJetPt10to20->SetXTitle("p_{T} (GeV)");
  fhdNdptTracksJetPt10to20->SetYTitle("dN/dp_{T} 1/GeV");
  fhdNdptTracksJetPt10to20->Sumw2();
  fListOfHistos->Add(fhdNdptTracksJetPt10to20);

  fhdNdptTracksJetPt20to30 = new TH1F("hdNdptTracksJetPt20to30","",25,0,25);
  fhdNdptTracksJetPt20to30->SetXTitle("p_{T} (GeV)");
  fhdNdptTracksJetPt20to30->SetYTitle("dN/dp_{T} 1/GeV");
  fhdNdptTracksJetPt20to30->Sumw2();
  fListOfHistos->Add(fhdNdptTracksJetPt20to30);

  fhdNdptTracksJetPt30to40 = new TH1F("hdNdptTracksJetPt30to40","",25,0,25);
  fhdNdptTracksJetPt30to40->SetXTitle("p_{T} (GeV)");
  fhdNdptTracksJetPt30to40->SetYTitle("dN/dp_{T} 1/GeV");
  fhdNdptTracksJetPt30to40->Sumw2();
  fListOfHistos->Add(fhdNdptTracksJetPt30to40);

  fhdNdptTracksJetPt40to60 = new TH1F("hdNdptTracksJetPt40to60","",25,0,25);
  fhdNdptTracksJetPt40to60->SetXTitle("p_{T} (GeV)");
  fhdNdptTracksJetPt40to60->SetYTitle("dN/dp_{T} 1/GeV");
  fhdNdptTracksJetPt40to60->Sumw2();
  fListOfHistos->Add(fhdNdptTracksJetPt40to60);


  fh1Xsec = new TProfile("h1Xsec","xsec from pyxsec.root",1,0,1); 
  fh1Xsec->SetXTitle("<#sigma>");
  fh1Xsec->Sumw2();
  fListOfHistos->Add( fh1Xsec ); 
  
  fh1Trials = new TH1F("h1Trials","trials from pyxsec.root",1,0,1);
  fh1Trials->SetXTitle("#sum{ntrials}");
  fh1Trials->Sumw2();
  fListOfHistos->Add( fh1Trials ); 
 
//   fSettingsTree   = new TTree("JetChemAnalysisSettings","Analysis Settings");
//   fSettingsTree->Branch("fUseLOConeJets",&fUseLOConeJets,"UseLOConeJets/O");
//   fSettingsTree->Branch("fUseLOConeMCJets",&fUseLOConeMCJets,"UseLOConeMCJets/O");
//   fSettingsTree->Branch("fUsePythiaJets",&fUsePythiaJets,"UsePythiaJets/O");
//   fSettingsTree->Branch("fConeRadius", &fConeRadius,"Rad/D");
//   fSettingsTree->Branch("fTrackPtCutJF", &fTrackPtCutJF,"TrackPt/D");
//   fSettingsTree->Branch("fFilterBitJF",&fFilterBitJF,"FilterBitJF/i");
//   fSettingsTree->Branch("fRequireITSRefitJF",&fRequireITSRefitJF,"RequireITSRefitJF/O");
//   fSettingsTree->Branch("fRejectK0TracksJF",&fRejectK0TracksJF,"RejectK0TracksJF/O");
//   fSettingsTree->Branch("fJetPtCut",&fJetPtCut,"JetPtCut/D");
//   fSettingsTree->Branch("fJetEtaCut",&fJetEtaCut,"JetEtaCut/D");
//   fSettingsTree->Branch("fFilterBit",&fFilterBit,"FilterBit/i");
//   fSettingsTree->Branch("fTrackPtCut",&fTrackPtCut,"TrackPtCut/D");
//   fSettingsTree->Branch("fTrackEtaCut",&fTrackEtaCut,"TrackEtaCut/D");
//   fSettingsTree->Branch("fUseOnFlyV0s",&fUseOnFlyV0s,"UseOnFlyV0s/O");
//   fSettingsTree->Branch("fCutnSigdEdx",&fCutnSigdEdx,"CutnSigdEdx/D");
//   fSettingsTree->Branch("fUseAODMCTracksForUE",&fUseAODMCTracksForUE,"UseAODMCTracksForUE/O");
//   fSettingsTree->Branch("fAreaReg",&fAreaReg,"AreaReg/D");
//   fSettingsTree->Branch("fAvgTrials",&fAvgTrials,"AvgTrials/D");
  // fListOfHistos->Add(fSettingsTree);
  
}

//____________________________________________________________________
void  AliAnalysisTaskJetChem::Terminate(Option_t */*option*/){

  // Terminate analysis
  //
  
  WriteSettings();
  
  // Update pointers reading them from the output slot
  fListOfHistos = dynamic_cast<TList*> (GetOutputData(1));
  if( !fListOfHistos ) {
    AliError("Histogram List is not available");
    return;
  }

  fhLeadingPt          = (TH1F*)fListOfHistos->FindObject("hLeadingPt");
  fhRegionSumPtMaxVsEt = (TH1F*)fListOfHistos->FindObject("hRegionSumPtMaxVsEt");
  fhRegionSumPtMinVsEt = (TH1F*)fListOfHistos->FindObject("hRegionSumPtMinVsEt");
  fhRegionMultMaxVsEt  = (TH1F*)fListOfHistos->FindObject("hRegionMultMaxVsEt");
  fhRegionMultMinVsEt  = (TH1F*)fListOfHistos->FindObject("hRegionMultMinVsEt");
    
   
  fhRegionSumPtMaxVsEt->Divide(fhLeadingPt);
  fhRegionSumPtMinVsEt->Divide(fhLeadingPt);
  fhRegionMultMaxVsEt->Divide(fhLeadingPt);
  fhRegionMultMinVsEt->Divide(fhLeadingPt);

    
  //Get Normalization
  fh1Xsec                = (TProfile*) fListOfHistos->FindObject("fh1Xsec");
  fh1Trials              = (TH1F*)fListOfHistos->FindObject("fh1Trials");
  
  //std::cout<<" fh1Xsec "<<fh1Xsec<<"  fh1Trials "<<fh1Trials<<std::endl;

  if(fh1Xsec && fh1Trials){

    Double_t xsec = fh1Xsec->GetBinContent(1);
    Double_t ntrials = fh1Trials->GetBinContent(1);
    Double_t normFactor = xsec/ntrials;
    Printf("xSec %f nTrials %f Norm %f \n",xsec,ntrials,normFactor);
    fhLeadingPt->Scale(normFactor);
  }

  
  if(fDebug > 1) AliInfo("End analysis");
}

// ----------------------------------------------------------------------------

void  AliAnalysisTaskJetChem::WriteSettings(){ 

  // fill settings tree
  // not on GRID !

  //fSettingsTree->Fill(); 

}

// ---------------------------------------------------------------------------

Bool_t AliAnalysisTaskJetChem::IsK0InvMass(const Double_t mass) const {

  // K0 mass ?  Use STAR params for bin counting & mass fit 

  const Double_t massK0    = 0.497;        // from fits
  const Double_t sigmaK0   = 0.0046;
  const Double_t nSigmaSignal = 3.5; // STAR parameters for bin counting

  if((massK0-nSigmaSignal*sigmaK0)<=mass &&  
     mass<=(massK0 + nSigmaSignal*sigmaK0)) return kTRUE;

  return kFALSE;
}

// ---------------------------------------------------------------------------

Bool_t AliAnalysisTaskJetChem::IsLambdaInvMass(const Double_t mass) const{

  // Lambda mass ?  


  if(1.1<mass && mass<1.13) return kTRUE; // FIXME - adjust range from fit
  return kFALSE;
}

// ---------------------------------------------------------------------------

Bool_t AliAnalysisTaskJetChem::IsAcceptedDCAK0(/*const Double_t dca*/) const{

  // DCA cut 

  return kTRUE; // FIXME - adjust cut
}


// ---------------------------------------------------------------------------

Bool_t AliAnalysisTaskJetChem::IsAcceptedDCALambda(/*const Double_t dca*/) const {

  // DCA cut

  return kTRUE; // FIXME - adjust cut
}

// ---------------------------------------------------------------------------

Bool_t AliAnalysisTaskJetChem::IsAccepteddEdx(const Double_t mom,const Double_t signal, AliPID::EParticleType n, const Double_t cutnSig) const{

  // apply TPC dE/dx cut similar as in AliTPCpidESD 
  // note: AliTPCidESD uses inner track param for momentum - not avaiable on AOD,  
  //       so we use global track momentum 
  // should use separate parametrisation for MC and data, but probably ALEPH param & 7% resolution used here anyway not the last word 
 
 
  const Double_t kBBMIP(50.);
  const Double_t kBBRes(0.07);
  //const Double_t kBBRange(5.);
  const Double_t kBBp1(0.76176e-1);
  const Double_t kBBp2(10.632);
  const Double_t kBBp3(0.13279e-4);
  const Double_t kBBp4(1.8631);
  const Double_t kBBp5(1.9479);

  Double_t mass=AliPID::ParticleMass(n); 
  Double_t betaGamma = mom/mass;

  const Float_t kmeanCorrection =0.1;
  Double_t bb = AliExternalTrackParam::BetheBlochAleph(betaGamma,kBBp1,kBBp2,kBBp3,kBBp4,kBBp5);
  Double_t meanCorrection =(1+(bb-1)*kmeanCorrection);
  Double_t bethe = bb * meanCorrection; // expected
  Double_t sigma = bethe * kBBRes;
        

  Double_t dedx = signal/kBBMIP; // measured

  Double_t nSig = (TMath::Abs(dedx - bethe))/sigma;
  
  if(nSig > cutnSig) return kFALSE; 

  return kTRUE;
}

// ----------------------------------------------------------------------------

void AliAnalysisTaskJetChem::CheckV0s(AliAODJet* jetVect, Int_t maxPtRegionIndex,Bool_t& foundK0){

  // loop over AOD V0s, fill masses etc

  Int_t nV0 = fAOD->GetNumberOfV0s();
  fhNV0s->Fill(nV0);  

  for(int i=0; i<fAOD->GetNumberOfV0s(); i++){ // loop over V0s
    
    AliAODv0* v0 = fAOD->GetV0(i);

    Bool_t isOnFly = v0->GetOnFlyStatus();
    isOnFly ? fhV0onFly->Fill(1) : fhV0onFly->Fill(0);
    
    if((fUseOnFlyV0s && isOnFly)  || (!fUseOnFlyV0s && !isOnFly) ){ // 'offline' V0s :  using vertex tracks, on-fly: take decay vertex into account during mom fit 
      
      Double_t massK0         = v0->MassK0Short();
      Double_t massLambda     = v0->MassLambda();
      Double_t massAntiLambda = v0->MassAntiLambda();
      Double_t radiusV0       = v0->RadiusV0();

      Double_t etaV0          = v0->Eta();
      Double_t fDCAV0ToVertex = v0->DcaV0ToPrimVertex();
      Double_t fDCADaughters  = v0->DcaV0Daughters();


      Double_t v0Mom[3];
      v0->PxPyPz(v0Mom);
      TVector3 v0MomVect(v0Mom);
      
      AliAODTrack *trackPos = (AliAODTrack *) (v0->GetSecondaryVtx()->GetDaughter(0));
      AliAODTrack *trackNeg = (AliAODTrack *) (v0->GetSecondaryVtx()->GetDaughter(1));   

      Double_t posMom[3], negMom[3];
      trackPos->PxPyPz(posMom);
      trackNeg->PxPyPz(negMom);
      
      TVector3 posMomVec(posMom);
      TVector3 negMomVec(negMom);
      
      Double_t dROpanV0 = posMomVec.DeltaR(negMomVec);
      
      AliAODPid*  aodPidPos = trackPos->GetDetPid();
      AliAODPid*  aodPidNeg = trackNeg->GetDetPid();

      Double_t  dEdxPos = aodPidPos->GetTPCsignal();
      Double_t  dEdxNeg = aodPidNeg->GetTPCsignal();

      Double_t momPos  = trackPos->P();
      Double_t momNeg  = trackNeg->P();
           
      fhdEdxVsMomV0->Fill(momPos,dEdxPos);
      fhdEdxVsMomV0->Fill(momNeg,dEdxNeg);
 
      Double_t pV0       = TMath::Sqrt(v0->Ptot2V0());
      Double_t ptV0      = TMath::Sqrt(v0->Pt2V0());
     
      fhV0InvMassK0->Fill(massK0);
      fhV0InvMassLambda->Fill(massLambda);
      fhV0InvMassAntiLambda->Fill(massAntiLambda);
    
      fhV0DCADaughters->Fill(fDCADaughters);
      fhV0Radius->Fill(radiusV0);
      fhV0DCAToVertex->Fill(fDCAV0ToVertex);
      fhdNdptV0->Fill(ptV0);


      // K0 signal before cuts 

      if(IsK0InvMass(massK0)){
	
	fhdNdptK0->Fill(ptV0);
	
	fhPtVsEtaK0->Fill(etaV0,ptV0);

	Double_t dRV0MC = AssociateV0MC(&v0MomVect,310); // K0
	if(dRV0MC < 0) dRV0MC = 0.99;
	fhdRV0MC->Fill(dRV0MC);
      }
 
      if(IsAcceptedDCAK0(/*fDCAV0ToVertex*/)){
	
	fhV0InvMassK0DCA->Fill(massK0);
	if(IsK0InvMass(massK0)) fhdNdptK0DCA->Fill(ptV0);

	
	if(IsAccepteddEdx(momPos,dEdxPos,AliPID::kPion,fCutnSigdEdx)) fhdEdxVsMomV0pidEdx->Fill(momPos,dEdxPos);
	if(IsAccepteddEdx(momNeg,dEdxNeg,AliPID::kPion,fCutnSigdEdx)) fhdEdxVsMomV0pidEdx->Fill(-1*momNeg,dEdxNeg);

	if(IsAccepteddEdx(momPos,dEdxPos,AliPID::kPion,fCutnSigdEdx) && 
	   IsAccepteddEdx(momNeg,dEdxNeg,AliPID::kPion,fCutnSigdEdx)){ 
	  
	  fhV0InvMassK0DCAdEdx->Fill(massK0);

	  if(IsK0InvMass(massK0)){
	    fhdNdptK0DCAdEdx->Fill(ptV0);
	    fhV0DCAToVertexK0->Fill(fDCAV0ToVertex);
	    fhdROpanK0VsPt->Fill(ptV0,dROpanV0);
	    foundK0 = kTRUE;
	  }
	}
	

	// AOD pid - cuts strongly into signal

	AliAODTrack::AODTrkPID_t mpPIDNeg = trackNeg->GetMostProbablePID();
	AliAODTrack::AODTrkPID_t mpPIDPos = trackPos->GetMostProbablePID();
	
	if(mpPIDNeg == AliAODTrack::kPion) fhdEdxVsMomV0piPID->Fill(momPos,dEdxPos);
	if(mpPIDPos == AliAODTrack::kPion) fhdEdxVsMomV0piPID->Fill(-1*momNeg,dEdxNeg);


	if( (mpPIDNeg == AliAODTrack::kPion) && (mpPIDPos == AliAODTrack::kPion) ){

	  fhV0InvMassK0DCAPID->Fill(massK0);
	}	
      }
      

      if(IsLambdaInvMass(massLambda) || (IsLambdaInvMass(massAntiLambda))){
	
	fhV0InvMassK0Lambda->Fill(massK0);      
      }


      // Lambda 

      if(IsAcceptedDCALambda(/*fDCAV0ToVertex*/) && 
	 IsAccepteddEdx(momPos,dEdxPos,AliPID::kProton,fCutnSigdEdx) && 
	 IsAccepteddEdx(momNeg,dEdxNeg,AliPID::kPion,fCutnSigdEdx)){
    
	fhV0InvMassLambdaDCAdEdx->Fill(massLambda);
      }
      
      if(IsAcceptedDCALambda(/*fDCAV0ToVertex*/) && 
	 IsAccepteddEdx(momNeg,dEdxNeg,AliPID::kProton,fCutnSigdEdx) && 
	 IsAccepteddEdx(momPos,dEdxPos,AliPID::kPion,fCutnSigdEdx)){
    
	fhV0InvMassAntiLambdaDCAdEdx->Fill(massAntiLambda);
      }


      // jet events 

      if(jetVect){
 
	Int_t regionV0vect = IsTrackInsideRegion(jetVect,&v0MomVect);

	// calc xi 
	Double_t jetpt     = jetVect->Pt();
	Double_t dPhiJetV0 = (jetVect->MomentumVector()->Vect()).DeltaPhi(v0MomVect);
	Double_t dRJetV0   = (jetVect->MomentumVector()->Vect()).DeltaR(v0MomVect);		

	Double_t z         = pV0/jetpt;
	Double_t xi        = TMath::Log(1/z);

	fhdPhiJetV0->Fill(dPhiJetV0);

	// K0

	if(IsAcceptedDCAK0(/*fDCAV0ToVertex*/) && 
	   IsAccepteddEdx(momPos,dEdxPos,AliPID::kPion,fCutnSigdEdx) && 
	   IsAccepteddEdx(momNeg,dEdxNeg,AliPID::kPion,fCutnSigdEdx)){
	  
	  fhV0InvMassK0JetEvt->Fill(massK0);

	  if(IsK0InvMass(massK0)){
	    fhdNdptK0JetEvt->Fill(ptV0);
	    fhdNdzK0->Fill(z); 
	    fhdNdxiK0->Fill(xi);

	    fhdPhiJetK0->Fill(dPhiJetV0);
	    fhdRJetK0->Fill(dRJetV0);

	    if(5<jetpt && jetpt<=10)       fhdNdzK05to10->Fill(z);
 	    else if(10<jetpt && jetpt<=20) fhdNdzK010to20->Fill(z);
 	    else if(20<jetpt && jetpt<=30) fhdNdzK020to30->Fill(z);
 	    else if(30<jetpt && jetpt<=40) fhdNdzK030to40->Fill(z);
 	    else if(40<jetpt && jetpt<=60) fhdNdzK040to60->Fill(z);
	  }
	}
	

	// Lambda
	
	if(IsAcceptedDCALambda(/*fDCAV0ToVertex*/) && 
	   IsAccepteddEdx(momPos,dEdxPos,AliPID::kProton,fCutnSigdEdx) && 
	   IsAccepteddEdx(momNeg,dEdxNeg,AliPID::kPion,fCutnSigdEdx)){
	  
	  fhV0InvMassLambdaJetEvt->Fill(massLambda);

	  if(IsLambdaInvMass(massLambda)){
	    fhdNdzLambda->Fill(z); 
	  }
	}
	
	if(IsAcceptedDCALambda(/*fDCAV0ToVertex*/) && 
	   IsAccepteddEdx(momNeg,dEdxNeg,AliPID::kProton,fCutnSigdEdx) && 
	   IsAccepteddEdx(momPos,dEdxPos,AliPID::kPion,fCutnSigdEdx)){
	  
	  fhV0InvMassAntiLambdaJetEvt->Fill(massAntiLambda);

	  if(IsLambdaInvMass(massAntiLambda)){
	    fhdNdzAntiLambda->Fill(z); 
	  }
	}
	
	// fill histos max region 
	
	if(regionV0vect != 0 && regionV0vect == maxPtRegionIndex){ // max region
	
	  // K0

	  if(IsAcceptedDCAK0(/*fDCAV0ToVertex*/) && 
	     IsAccepteddEdx(momPos,dEdxPos,AliPID::kPion,fCutnSigdEdx) && 
	     IsAccepteddEdx(momNeg,dEdxNeg,AliPID::kPion,fCutnSigdEdx)){ // K0 cuts
	    
	    fhV0InvMassK0Max->Fill(massK0);
	    
	    if(IsK0InvMass(massK0)){
	      fhdNdzK0Max->Fill(z);
	      fhdNdxiK0Max->Fill(xi);
	      fhdNdptK0Max->Fill(ptV0);
	    }
	  }

	  // Lambda
	  
	  if(IsAcceptedDCALambda(/*fDCAV0ToVertex*/) && 
	     IsAccepteddEdx(momPos,dEdxPos,AliPID::kProton,fCutnSigdEdx) && 
	     IsAccepteddEdx(momNeg,dEdxNeg,AliPID::kPion,fCutnSigdEdx)){
	    
	    fhV0InvMassLambdaMax->Fill(massLambda);
	  }

	  if(IsAcceptedDCALambda(/*fDCAV0ToVertex*/) && 
	     IsAccepteddEdx(momNeg,dEdxNeg,AliPID::kProton,fCutnSigdEdx) && 
	     IsAccepteddEdx(momPos,dEdxPos,AliPID::kPion,fCutnSigdEdx)){
	    
	    fhV0InvMassAntiLambdaMax->Fill(massAntiLambda);
	  }
	  
	  if(IsAcceptedDCALambda(/*fDCAV0ToVertex*/) && 
	     ((IsAccepteddEdx(momPos,dEdxPos,AliPID::kProton,fCutnSigdEdx) && 
	       IsAccepteddEdx(momNeg,dEdxNeg,AliPID::kPion,fCutnSigdEdx)) ||
	      (IsAccepteddEdx(momNeg,dEdxNeg,AliPID::kProton,fCutnSigdEdx) && 
	       IsAccepteddEdx(momPos,dEdxPos,AliPID::kPion,fCutnSigdEdx)))){
	       
	    if(IsLambdaInvMass(massLambda) || (IsLambdaInvMass(massAntiLambda))){
	      fhdNdzLambdaMax->Fill(z);
	      fhdNdxiLambdaMax->Fill(xi);
	      fhdNdptLambdaMax->Fill(ptV0);
	    }
	  }
	}
	
	// fill histos min region

	if(regionV0vect != 0 && regionV0vect != maxPtRegionIndex){ // min region 
	
	  // K0
  
	  if(IsAcceptedDCAK0(/*fDCAV0ToVertex*/) && 
	     IsAccepteddEdx(momPos,dEdxPos,AliPID::kPion,fCutnSigdEdx) && 
	     IsAccepteddEdx(momNeg,dEdxNeg,AliPID::kPion,fCutnSigdEdx)){
	    
	    fhV0InvMassK0Min->Fill(massK0);
	    
	    if(IsK0InvMass(massK0)){
	      fhdNdzK0Min->Fill(z);
	      fhdNdxiK0Min->Fill(xi);
	      fhdNdptK0Min->Fill(ptV0);
	    }
	  }
       

	  // Lambda
	  
	  
	  if(IsAcceptedDCALambda(/*fDCAV0ToVertex*/) && 
	     IsAccepteddEdx(momPos,dEdxPos,AliPID::kProton,fCutnSigdEdx) && 
	     IsAccepteddEdx(momNeg,dEdxNeg,AliPID::kPion,fCutnSigdEdx)){
	    
	    fhV0InvMassLambdaMin->Fill(massLambda);
	  }

	  if(IsAcceptedDCALambda(/*fDCAV0ToVertex*/) && 
	     IsAccepteddEdx(momNeg,dEdxNeg,AliPID::kProton,fCutnSigdEdx) && 
	     IsAccepteddEdx(momPos,dEdxPos,AliPID::kPion,fCutnSigdEdx)){
	    
	    fhV0InvMassAntiLambdaMin->Fill(massAntiLambda);
	  }
	  
	  if(IsAcceptedDCALambda(/*fDCAV0ToVertex*/) && 
	     ((IsAccepteddEdx(momPos,dEdxPos,AliPID::kProton,fCutnSigdEdx) && 
	       IsAccepteddEdx(momNeg,dEdxNeg,AliPID::kPion,fCutnSigdEdx)) ||
	      (IsAccepteddEdx(momNeg,dEdxNeg,AliPID::kProton,fCutnSigdEdx) && 
	       IsAccepteddEdx(momPos,dEdxPos,AliPID::kPion,fCutnSigdEdx)))){
	    
	    if(IsLambdaInvMass(massLambda) || (IsLambdaInvMass(massAntiLambda))){
	      fhdNdzLambdaMin->Fill(z);
	      fhdNdxiLambdaMin->Fill(xi);
	      fhdNdptLambdaMin->Fill(ptV0);
	    }
	  }
	}
	

	// jet region 

	if(regionV0vect == 0){ // jet region
	
	  //Double_t dRJetV0 = jetVect->DeltaR(v0MomVect); 
	  
	  if(dRJetV0 <= fConeRadius){ 
	  
	    // fill histos jet region  

	    // K0

	    if(IsAcceptedDCAK0(/*fDCAV0ToVertex*/) && 
	       IsAccepteddEdx(momPos,dEdxPos,AliPID::kPion,fCutnSigdEdx) && 
	       IsAccepteddEdx(momNeg,dEdxNeg,AliPID::kPion,fCutnSigdEdx)){ // K0 cuts
	      
	      fhV0InvMassK0Jet->Fill(massK0);
 	    
	      if(IsK0InvMass(massK0)){
		fhdNdzK0Jet->Fill(z);
		fhdNdxiK0Jet->Fill(xi);
		fhdNdptK0Jet->Fill(ptV0);
	      }
	    }


	    // Lambda 
 	  
	    if(IsAcceptedDCALambda(/*fDCAV0ToVertex*/) && 
	       IsAccepteddEdx(momPos,dEdxPos,AliPID::kProton,fCutnSigdEdx) && 
	       IsAccepteddEdx(momNeg,dEdxNeg,AliPID::kPion,fCutnSigdEdx)){
	    
	      fhV0InvMassLambdaJet->Fill(massLambda);
	    }	  
	    if(IsAcceptedDCALambda(/*fDCAV0ToVertex*/) && 
	       IsAccepteddEdx(momNeg,dEdxNeg,AliPID::kProton,fCutnSigdEdx) && 
	       IsAccepteddEdx(momPos,dEdxPos,AliPID::kPion,fCutnSigdEdx)){
	      
	      fhV0InvMassAntiLambdaJet->Fill(massAntiLambda);
	    }
	    
	    if(IsAcceptedDCALambda(/*fDCAV0ToVertex*/) && 
	       ((IsAccepteddEdx(momPos,dEdxPos,AliPID::kProton,fCutnSigdEdx) && 
		 IsAccepteddEdx(momNeg,dEdxNeg,AliPID::kPion,fCutnSigdEdx)) ||
		(IsAccepteddEdx(momNeg,dEdxNeg,AliPID::kProton,fCutnSigdEdx) && 
		 IsAccepteddEdx(momPos,dEdxPos,AliPID::kPion,fCutnSigdEdx)))){
	      
	      if(IsLambdaInvMass(massLambda) || (IsLambdaInvMass(massAntiLambda))){
		fhdNdzLambdaJet->Fill(z);
		fhdNdxiLambdaJet->Fill(xi);
		fhdNdptLambdaJet->Fill(ptV0);
	      }
	    }
	  }
	}
      }
    }
  }
}

// ----------------------------------------------------------------------------

void AliAnalysisTaskJetChem::CheckMCParticles(AliAODJet* jetVect,Int_t maxPtRegionIndex,Bool_t& isK0Event){

  // histos for generated MC  

  TClonesArray* farray = (TClonesArray*)fAOD->FindListObject("mcparticles");
       
  if(!farray){
    AliInfo("no mcparticles branch"); 
    return;
  }

  Int_t ntrks = farray->GetEntries();
  if (fDebug>1) AliInfo(Form("check MC particles, tracks %d \n",ntrks));

  Int_t pythiaPID = GetPythiaProcessID();

  isK0Event          = kFALSE;
  Bool_t ispEvent    = kFALSE;
  Bool_t ispBarEvent = kFALSE;
  Bool_t isKchEvent  = kFALSE;

  Bool_t isQuarkHardScatteringEvent = IsQuarkHardScatteringEvent(pythiaPID);
  Bool_t isGluonHardScatteringEvent = IsGluonHardScatteringEvent(pythiaPID);

  for(Int_t i =0 ; i<ntrks; i++){ // mc tracks loop

    AliAODMCParticle* mctrk = (AliAODMCParticle*)farray->At(i);
    
    Double_t trackPt = mctrk->Pt(); 

    //Cuts
    if (!(mctrk->IsPhysicalPrimary())) continue;
 
    if ((trackPt < fTrackPtCut) || (TMath::Abs(mctrk->Eta()) > fTrackEtaCut)) continue; 

    // OB PID histo
    Int_t pdg = mctrk->GetPdgCode();
    

    FillPIDhisto(fhPIDMC,pdg);
    if(isQuarkHardScatteringEvent) FillPIDhisto(fhPIDMC_quarkEv,pdg);
    if(isGluonHardScatteringEvent) FillPIDhisto(fhPIDMC_gluonEv,pdg);

    if(pdg == 111)              fhdNdptpi0MC->Fill(trackPt);
    if(pdg == 22)               fhdNdptgammaMC->Fill(trackPt);
    if(TMath::Abs(pdg) == 211)  fhdNdptchPiMC->Fill(trackPt);
    if(pdg == 310)              fhdNdptK0MC->Fill(trackPt);
    if(TMath::Abs(pdg) == 321)  fhdNdptchKMC->Fill(trackPt);
    if(pdg == 2212)             fhdNdptpMC->Fill(trackPt);
    if(pdg == -2212)            fhdNdptpBarMC->Fill(trackPt);
    if(pdg == 3122)             fhdNdptLambdaMC->Fill(trackPt);
    if(pdg == -3122)            fhdNdptLambdaBarMC->Fill(trackPt);
    if(pdg == 3332)             fhdNdptOmegaMC->Fill(trackPt);
    if(pdg == -3332)            fhdNdptOmegaBarMC->Fill(trackPt);
    
    if(pdg == 310)             isK0Event   = kTRUE;
    if(TMath::Abs(pdg) == 321) isKchEvent  = kTRUE;
    if(pdg == 2212)            ispEvent    = kTRUE;
    if(pdg == -2212)           ispBarEvent = kTRUE;
    
    
    Int_t pdgMotherChK = -1;
    Int_t pdgMotherK0  = -1;
    Int_t pdgGrandMotherK0 = -1;

    if(TMath::Abs(pdg) == 321){ // chK
      Int_t labelMother = mctrk->GetMother();
      AliAODMCParticle* mctrkMother    = NULL;
      
      for(Int_t k=0 ; k<ntrks; k++){
      
	AliAODMCParticle* mctrk2 = (AliAODMCParticle*)farray->At(k);
	if(mctrk2->GetLabel() == labelMother) mctrkMother = mctrk2;
      }
      
      if(mctrkMother) pdgMotherChK = mctrkMother->GetPdgCode();
      FillPIDhisto(fhPIDMCMotherChK,pdgMotherChK);

      //printf("pdgMotherChK %d \n",pdgMotherChK);
    }
    
    if(pdg == 310){ // K0
      Int_t labelMother = mctrk->GetMother();
      AliAODMCParticle* mctrkMother = NULL;
      
      for(Int_t k=0 ; k<ntrks; k++){
	AliAODMCParticle* mctrk2 = (AliAODMCParticle*)farray->At(k);
	if(mctrk2->GetLabel() == labelMother) mctrkMother = mctrk2;
      }
      
      if(mctrkMother) pdgMotherK0 = mctrkMother->GetPdgCode();
      FillPIDhisto(fhPIDMCMotherK0,pdgMotherK0);

      Int_t labelGrandMother = -1; 
      if(mctrkMother) mctrkMother->GetMother();
      AliAODMCParticle* mctrkGrandMother = NULL;
      
      for(Int_t k=0 ; k<ntrks; k++){
	AliAODMCParticle* mctrk2 = (AliAODMCParticle*)farray->At(k);
	if(mctrk2->GetLabel() == labelGrandMother) mctrkGrandMother = mctrk2;
      }
      
      if(mctrkGrandMother) pdgGrandMotherK0 = mctrkGrandMother->GetPdgCode();
      FillPIDhisto(fhPIDMCGrandMotherK0,pdgGrandMotherK0);
    }
    
    if(jetVect){ // jet event 

      FillPIDhisto(fhPIDMCAll,pdg);
    
      TVector3 partVect(mctrk->Px(), mctrk->Py(), mctrk->Pz());
      
      Int_t region = IsTrackInsideRegion(jetVect, &partVect );  
      
      //Double_t deltaPhi = jetVect->DeltaPhi(partVect); //+k270rad;
      Double_t deltaPhi = (jetVect->MomentumVector()->Vect()).DeltaPhi(partVect);
      if( deltaPhi > 2.*TMath::Pi() )  deltaPhi-= 2.*TMath::Pi();
      //Double_t deltaR = jetVect->DeltaR(partVect); //+k270rad;
      Double_t deltaR = (jetVect->MomentumVector()->Vect()).DeltaR(partVect);

      // calc xi 
      Double_t jetpt = jetVect->Pt(); 
      //Double_t pV0   = partVect.Mag(); 
      //Double_t ptV0  = partVect.Pt();
   

      Double_t z  = trackPt / jetpt;
      Double_t xi = TMath::Log(1/z);

      if(!(mctrk->Charge() == 0 || mctrk->Charge()==-99)) fhdNdxiMC->Fill(xi);

      if(pdg == 310){ // K0
	fhdPhiJetK0MC->Fill(deltaPhi);
	fhdRJetK0MC->Fill(deltaR);
	fhdNdxiK0MC->Fill(xi);
	fhdNdzK0MC->Fill(z);
	fhdNdptK0MCJetEvt->Fill(trackPt);
      }

      
      if(region != 0 && region == maxPtRegionIndex){ // max region
	
	if(TMath::Abs(pdg) == 211)  fhdNdptchPiMCMax->Fill(trackPt);
	if(pdg == 310)              fhdNdptK0MCMax->Fill(trackPt);
	if(TMath::Abs(pdg) == 321)  fhdNdptchKMCMax->Fill(trackPt);
	if(pdg == 2212)             fhdNdptpMCMax->Fill(trackPt);
	if(pdg == -2212)            fhdNdptpBarMCMax->Fill(trackPt);
	if(pdg == 3122)             fhdNdptLambdaMCMax->Fill(trackPt);
	if(pdg == -3122)            fhdNdptLambdaBarMCMax->Fill(trackPt);
      }
      
      if(region != 0 && region != maxPtRegionIndex){ // min region
	
	FillPIDhisto(fhPIDMCMin,pdg);
	
	if(TMath::Abs(pdg) == 211)  fhdNdptchPiMCMin->Fill(trackPt);
	if(pdg == 310)              fhdNdptK0MCMin->Fill(trackPt);
	if(TMath::Abs(pdg) == 321)  fhdNdptchKMCMin->Fill(trackPt);
	if(pdg == 2212)             fhdNdptpMCMin->Fill(trackPt);
	if(pdg == -2212)            fhdNdptpBarMCMin->Fill(trackPt);
	if(pdg == 3122)             fhdNdptLambdaMCMin->Fill(trackPt);
	if(pdg == -3122)            fhdNdptLambdaBarMCMin->Fill(trackPt);

	if(pdg == 3332)  fhdNdptOmegaMCMin->Fill(trackPt);
	if(pdg == -3332) fhdNdptOmegaBarMCMin->Fill(trackPt);
      }

      // trans region
      if(region != 0){
	if(TMath::Abs(pdg) == 321) FillPIDhisto(fhPIDMCMotherChKTrans,pdgMotherChK);
	if(pdg == 310){
	  FillPIDhisto(fhPIDMCMotherK0Trans,pdgMotherK0);
	  FillPIDhisto(fhPIDMCGrandMotherK0Trans,pdgGrandMotherK0);
	}
      }
   
      if(region == 0){ //  jet region ?
	
	FillPIDhisto(fhPIDMCJet,pdg);
	
	//Double_t dRJetV0 = jetVect->DeltaR(partVect); 
	Double_t dRJetV0 =  (jetVect->MomentumVector()->Vect()).DeltaR(partVect);
	
	if(dRJetV0 <= fConeRadius){ 
	  
	  if(pdg == 310){ // K0
	    
	    fhdNdptK0MCJet->Fill(trackPt);
	    fhdNdxiK0MCJet->Fill(xi);
	    fhdNdzK0MCJet->Fill(z);
	  }
	
	  if(TMath::Abs(pdg) == 211)  fhdNdptchPiMCJet->Fill(trackPt);
	  if(TMath::Abs(pdg) == 321)  fhdNdptchKMCJet->Fill(trackPt);
	  if(pdg == 2212)             fhdNdptpMCJet->Fill(trackPt);
	  if(pdg == -2212)            fhdNdptpBarMCJet->Fill(trackPt);
	  if(pdg == 3122)             fhdNdptLambdaMCJet->Fill(trackPt);
	  if(pdg == -3122)            fhdNdptLambdaBarMCJet->Fill(trackPt);  
	}
      }
    }
  }


  FillPythiaIDhisto(fhPythiaProcess,pythiaPID);
  if(isK0Event)   FillPythiaIDhisto(fhPythiaProcessK0,pythiaPID);
  if(isKchEvent)  FillPythiaIDhisto(fhPythiaProcessKch,pythiaPID);
  if(ispEvent)    FillPythiaIDhisto(fhPythiaProcessp,pythiaPID);
  if(ispBarEvent) FillPythiaIDhisto(fhPythiaProcesspbar,pythiaPID);

}
	
// ----------------------------------------------------------------------------

Double_t AliAnalysisTaskJetChem::AssociateV0MC(const TVector3* V0Mom,const Int_t pdg){

  // find closest MC gen. particle for V0 vector

  TClonesArray* farray = (TClonesArray*)fAOD->FindListObject("mcparticles");
       
  if(!farray){
    AliInfo("no mcparticles branch"); 
    return -1;
  }

  Double_t dRmin = -1;

  Int_t ntrks = farray->GetEntries();

  for(Int_t i =0 ; i<ntrks; i++){

    AliAODMCParticle* mctrk = (AliAODMCParticle*)farray->At(i);
    
    //Cuts
    if (!(mctrk->IsPhysicalPrimary())) continue;

    Int_t pdgtrk = mctrk->GetPdgCode();
    
    if(pdgtrk != pdg) continue;

    TVector3 partVect(mctrk->Px(), mctrk->Py(), mctrk->Pz());

    Double_t dR = V0Mom->DeltaR(partVect);

    if(dRmin<0) dRmin = dR; // initialize
    
    if(dR < dRmin) dRmin = dR;
      
  }
  
  return dRmin;
}

// ----------------------------------------------------------------------------

void AliAnalysisTaskJetChem::CompLeadingJets(AliAODJet* jetLeadingAOD,AliAODJet* jetLeadingMC,const Int_t pythiaPID,
					     const Bool_t foundK0AOD,const Bool_t foundK0MC){

  // leading jet properties

  Double_t ptLeadingAOD      = -1;
  Double_t etaLeadingAOD     = -1;
  Double_t phiLeadingAOD     = -1;
  Int_t    nTracksLeadingAOD = -1;

  Double_t ptLeadingMC       = -1; 
  Double_t etaLeadingMC      = -1;
  Double_t phiLeadingMC      = -1;
  Int_t    nTracksLeadingMC  = -1;

  if(jetLeadingAOD){
    if(jetLeadingAOD->Pt()>fJetPtCut && TMath::Abs(jetLeadingAOD->Eta())<fJetEtaCut){

      ptLeadingAOD      = jetLeadingAOD->Pt();
      etaLeadingAOD     = jetLeadingAOD->Eta();
      phiLeadingAOD     = jetLeadingAOD->Phi();
      nTracksLeadingAOD = jetLeadingAOD->GetRefTracks()->GetEntriesFast();

      Double_t radiusAOD = GetJetRadius(jetLeadingAOD,0.8);
      fhnTracksJetVsPtAOD->Fill(ptLeadingAOD,nTracksLeadingAOD);
      fhRadiusJetVsPtAOD->Fill(ptLeadingAOD,radiusAOD);
      if(IsQuarkHardScatteringEvent(pythiaPID)) fhnTracksJetVsPtAODquarkEv->Fill(ptLeadingAOD,nTracksLeadingAOD);
      
      if(foundK0AOD){
	
	fhnTracksJetVsPtAODK0->Fill(ptLeadingAOD,nTracksLeadingAOD);
	if(IsQuarkHardScatteringEvent(pythiaPID)) fhnTracksJetVsPtAODK0quarkEv->Fill(ptLeadingAOD,nTracksLeadingAOD);
	fhRadiusJetVsPtAODK0->Fill(ptLeadingAOD,radiusAOD);
      }
      

      // check if p/Kch in jet 

      Bool_t foundpKch = kFALSE; 
      Int_t nTracksJet = jetLeadingAOD->GetRefTracks()->GetEntriesFast();

      for(int i=0; i<nTracksJet; i++){

	AliAODTrack* track = (AliAODTrack*) jetLeadingAOD->GetRefTracks()->At(i);

	Double_t mom  = track->P();
	
	AliAODPid* aodPid = track->GetDetPid();
	Double_t   dEdx   = aodPid->GetTPCsignal();

	if(IsAccepteddEdx(mom,dEdx,AliPID::kKaon,fCutnSigdEdx) ||  
	   IsAccepteddEdx(mom,dEdx,AliPID::kProton,fCutnSigdEdx)){

	  foundpKch = kTRUE; 
	}
      } // track loop
      
      if(foundpKch){
	fhnTracksJetVsPtAODpKch->Fill(ptLeadingAOD,nTracksLeadingAOD);
	fhRadiusJetVsPtAODpKch->Fill(ptLeadingAOD,radiusAOD);
      }
    }
  }


  if(jetLeadingMC){
    if(jetLeadingMC->Pt()>fJetPtCut && TMath::Abs(jetLeadingMC->Eta())<fJetEtaCut){

      ptLeadingMC      = jetLeadingMC->Pt();
      etaLeadingMC     = jetLeadingMC->Eta();
      phiLeadingMC     = jetLeadingMC->Phi();
      nTracksLeadingMC = jetLeadingMC->GetRefTracks()->GetEntriesFast();    
      
      Double_t radiusMC  = GetJetRadius(jetLeadingMC,0.8);
      fhnTracksJetVsPtMC->Fill(ptLeadingMC,nTracksLeadingMC);
      fhRadiusJetVsPtMC->Fill(ptLeadingMC,radiusMC);
      if(IsQuarkHardScatteringEvent(pythiaPID)) fhnTracksJetVsPtMCquarkEv->Fill(ptLeadingMC,nTracksLeadingMC);
      
      if(foundK0MC){
	
	fhnTracksJetVsPtMCK0->Fill(ptLeadingMC,nTracksLeadingMC);
	if(IsQuarkHardScatteringEvent(pythiaPID)) fhnTracksJetVsPtMCK0quarkEv->Fill(ptLeadingMC,nTracksLeadingMC);
	fhRadiusJetVsPtMCK0->Fill(ptLeadingMC,radiusMC);
      }
    }
  }
  
  if(jetLeadingAOD && jetLeadingMC){
    
    //std::cout<<" comp: leading jetPt AOD "<<ptLeadingAOD<<" MC "<<ptLeadingMC<<std::endl;
    //if(jetLeadingAOD && jetLeadingMC) 
    //std::cout<<" leading jet eta AOD "<<jetLeadingAOD->Eta()<<" MC "<<jetLeadingMC->Eta()<<std::endl;

    
    if(jetLeadingMC->Pt()>fJetPtCut && jetLeadingAOD->Pt()>fJetPtCut && 
       TMath::Abs(jetLeadingMC->Eta())<fJetEtaCut && TMath::Abs(jetLeadingAOD->Eta())<fJetEtaCut){
      
      fhLeadingPtAODvsMC->Fill(ptLeadingMC,ptLeadingAOD);
      fhLeadingEtaAODvsMC->Fill(etaLeadingMC,etaLeadingAOD);
      fhLeadingPhiAODvsMC->Fill(phiLeadingMC,phiLeadingAOD);
      fhnTracksLeadingAODvsMC->Fill(nTracksLeadingMC,nTracksLeadingAOD);
      
      TLorentzVector *mom4MC  = jetLeadingMC->MomentumVector();
      TLorentzVector *mom4AOD = jetLeadingAOD->MomentumVector();
      
      Double_t dR = mom4MC->DeltaR(*mom4AOD);
      fhLeadingdRAODMC->Fill(dR);
      
      if(dR<0.4) fhLeadingPtAODvsMCdRcut->Fill(ptLeadingMC,ptLeadingAOD);
      
      Double_t dPt = ptLeadingAOD - ptLeadingMC;
      Double_t dnTracks = nTracksLeadingAOD - nTracksLeadingMC;
      
      fhdnTracksVsdPtLeadingAODMC->Fill(dPt,dnTracks);  
    }
  }
  
}

// ---------------------------------------------------------------------------

// void AliAnalysisTaskJetChem::CheckK0(){

//   TClonesArray* farray = (TClonesArray*)fAOD->FindListObject("mcparticles");
  
//   if(!farray){
//     AliInfo("no mcparticles branch"); 
//     return;
//   }
  
  
//   Int_t ntrks = farray->GetEntries();
    
//   for(Int_t i=0; i<ntrks; i++){ // trk loop
    
//     AliAODMCParticle* mctrk = (AliAODMCParticle*)farray->At(i);
    
//     Int_t pdg = mctrk->GetPdgCode();
    
//     if(pdg != 310) continue; // K0

//     cout<<" MC K0: physPrim "<<mctrk->IsPhysicalPrimary()<<endl;

//     Int_t labelMother = mctrk->GetMother();

//     cout<<" found K00, label mother "<<labelMother<<endl;

//     AliAODMCParticle* mctrkMother = NULL;
//     Int_t pdgMother = -1;
    
//     for(Int_t k=0 ; k<ntrks; k++){
      
//       mctrkMother = (AliAODMCParticle*)farray->At(k);
//       if(mctrkMother->GetLabel() == labelMother) break;
//     }
    
//     pdgMother = mctrkMother->GetPdgCode();
//     cout<<" K0 mother pdg "<<pdgMother<<" GID "<<fpdgdb->ConvertPdgToGeant3(pdgMother)<<" isPrimary "<<mctrkMother->IsPrimary()<<endl;
//     //cout<<" mother name "<<mctrkMother->GetName()<<endl;
    
//   }
// }

// ---------------------------------------------------------------------------

// void AliAnalysisTaskJetChem::CheckK0Stack(){
  
//   AliMCEventHandler* mcHandler = dynamic_cast<AliMCEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
//   if (!mcHandler) {
//     Printf("ERROR: Could not retrieve MC event handler");
//     return;
//   }
  
//   AliMCEvent* mcEvent = mcHandler->MCEvent();
//   if (!mcEvent) {
//     Printf("ERROR: Could not retrieve MC event");
//     return;
//   }
  
//    AliStack* mcStack = mcEvent->Stack();//Load Stack
//    AliGenPythiaEventHeader*  pythiaGenHeader = AliAnalysisHelperJetTasks::GetPythiaEventHeader(mcEvent);

//    Int_t nTracksMC = mcStack->GetNtrack();

//    Bool_t foundK0 = kFALSE;

//    for (Int_t iTracks = 0; iTracks < nTracksMC; iTracks++) {
//      //Cuts
//      if(!(mcStack->IsPhysicalPrimary(iTracks))) continue;
     
//      TParticle* mctrk = mcStack->Particle(iTracks);
        
//      Int_t pdg = mctrk->GetPdgCode();

//     if ((mctrk->Pt() < 2*fTrackPtCut) || (TMath::Abs(mctrk->Eta()) > fTrackEtaCut )) continue; 

//      if(pdg == 310) foundK0 = kTRUE;

//      if(pdg != 310) continue; // K0 short

//      cout<<" check K0 "<<endl;

//      Int_t indexMother = -999;
//      TParticle* mctrkMother = mctrk;
//      TParticle* mctrkFirstGenMother = NULL; // last mother which is not the primary proton

//      Int_t nAncestors = 0;

//      while(indexMother != -1){

//        indexMother = mctrkMother->GetFirstMother();
     
//        if(indexMother != -1){
// 	 mctrkFirstGenMother = mctrkMother;
// 	 mctrkMother  =  mcStack->Particle(indexMother);
// 	 nAncestors++;
//        }

//        cout<<" nAncestors "<<nAncestors<<" pdg mother "<<mctrkMother->GetPdgCode()<<" name "<<mctrkMother->GetName()<<endl;
//      }
//      cout<<" pdg firstGenMother "<<mctrkFirstGenMother->GetPdgCode()<<" name "<<mctrkFirstGenMother->GetName()<<endl;



//      cout<<" pythiaGenHeader "<<pythiaGenHeader<<endl;
//      cout<<" Pythia Process type "<<pythiaGenHeader->ProcessType()<<" ptHard "<<pythiaGenHeader->GetPtHard()<<endl;


//      fhPythiaProcess->Fill(pythiaGenHeader->ProcessType());
//      if(foundK0) fhPythiaProcess_K0->Fill(pythiaGenHeader->ProcessType());

//      //Int_t indexGrandMother      = mctrkMother->GetFirstMother();
//      //cout<<" indexGrandMother "<<indexGrandMother<<endl;
     
//      //if(indexGrandMother>-1){
//      //  TParticle* mctrkGrandMother = mcStack->Particle(indexGrandMother);
//      //  cout<<" pdg grandMother "<<mctrkGrandMother->GetPdgCode()<<" name "<<mctrkGrandMother->GetName()<<endl;
//      // }

//    }
// }

// ---------------------------------------------------------------------------------

Bool_t  AliAnalysisTaskJetChem::IsQuarkHardScatteringEvent(const Int_t PID){

  // Pythia Manual sec. 8.2.1 : 
  // if Pythia PID = 92,93,94   event is diffractive 
  // if Pythia PID = 13, 28, 68 hard scattering products are gg, qg 
  // if Pythia PID = 11, 12, 53 hard scattering products are qq,  q\bar{q}  

  if(PID == 92 || PID == 93 || PID == 94 || PID==95) return kFALSE;
  else if(PID == 13 || PID == 28 || PID == 68) return kFALSE;
  else if(PID == 11 || PID == 12 || PID == 53) return kTRUE;
  else{
    AliInfo(Form("unknown Pythia PID %d",PID));
  }

  return kFALSE;
}


// ---------------------------------------------------------------------------------

Bool_t  AliAnalysisTaskJetChem::IsGluonHardScatteringEvent(const Int_t PID){

  // Pythia Manual sec. 8.2.1 : 
  // if Pythia PID = 92,93,94   event is diffractive 
  // if Pythia PID = 95: low pt event (MPI) 
  // if Pythia PID = 13, 28, 68 hard scattering products are gg, qg 
  // if Pythia PID = 11, 12, 53 hard scattering products are qq,  q\bar{q}  


  if(PID == 92 || PID == 93 || PID == 94 || PID == 95) return kFALSE;
  else if(PID == 13 || PID == 68) return kTRUE;
  else if(PID == 28) return kFALSE; // mixed gq final state
  else if(PID == 11 || PID == 12 || PID == 53) return kFALSE;
  else{
    AliInfo(Form("unknown Pythia PID %d",PID));
  }

  return kFALSE;
}

// ---------------------------------------------------------------------------------

Bool_t  AliAnalysisTaskJetChem::IsDiffractiveEvent(const Int_t PID){

  // Pythia Manual sec. 8.2.1 : 
  // if Pythia PID = 92,93,94   event is diffractive 
  // if Pythia PID = 13, 28, 68 hard scattering products are gg, qg 
  // if Pythia PID = 11, 12, 53 hard scattering products are qq,  q\bar{q}  

  if(PID == -1) return kFALSE;

  if(PID == 13 || PID == 28 || PID == 68) return kFALSE;
  else if(PID == 11 || PID == 12 || PID == 53) return kFALSE;
  else if(PID == 92 || PID == 93 || PID == 94 || PID==95) return kTRUE;
  else{
    AliInfo(Form("unknown Pythia PID %d",PID));
  }

  return kFALSE;

}


// ----------------------------------------------------------------------------------

Int_t AliAnalysisTaskJetChem::GetPythiaProcessID(){

  // Pythia PID for this event 
  
  AliMCEventHandler* mcHandler = dynamic_cast<AliMCEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
  if (!mcHandler) {
    //Printf("ERROR: Could not retrieve MC event handler");
    return -1;
  }
  
  AliMCEvent* mcEvent = mcHandler->MCEvent();
  if (!mcEvent) {
    AliInfo("could not retrieve MC event");
    return -1;
  }
  
  AliGenPythiaEventHeader*  pythiaGenHeader = AliAnalysisHelperJetTasks::GetPythiaEventHeader(mcEvent);
  
  if (!pythiaGenHeader) {
    AliInfo("Could not retrieve pythiaEventHeader");
    return -1;
  }

  Int_t  pid = pythiaGenHeader->ProcessType(); 

  return pid;

}

// ----------------------------------------------------------------------------------

void AliAnalysisTaskJetChem::GetJetTracksResum(TList* list, AliAODJet* jet, const Double_t radius){
  

  if(!jet) return; // no jet in acc in event

  // list of AOD tracks in jet cone, using cone axis and distance axis-track (and not trackrefs)

  Int_t nTracks = fAOD->GetNTracks();

  if(!nTracks) return;

  Double_t jetMom[3];
  jet->PxPyPz(jetMom);
  TVector3 jet3mom(jetMom);
    
  
  // standard cuts + ITS refit = stdrd PWG4 cut (compose them for productions with old ESDFilter task) 

  for (Int_t itrack=0; itrack<nTracks; itrack++) {

    AliAODTrack* track = fAOD->GetTrack(itrack);
    
    UInt_t status = track->GetStatus();

    if(!track->TestFilterBit(fFilterBitJF)) continue; // track cut selection
    if(fRequireITSRefitJF &&  ((status&AliESDtrack::kITSrefit)==0)) continue;  

    if(TMath::Abs(track->Eta()) > fTrackEtaCut) continue; 
    if( track->Pt() < fTrackPtCutJF ) continue;

    Double_t trackMom[3];
    track->PxPyPz(trackMom);
    TVector3 track3mom(trackMom);

    Double_t dR = jet3mom.DeltaR(track3mom);

    if(dR<radius){

      list->Add(track);
      
    }
  }
}

// ----------------------------------------------------------------------------------

void AliAnalysisTaskJetChem::GetJetTracksTrackrefs(TList* list, AliAODJet* jet){
  
  if(!jet) return; // no jet in acc in event

  // list of AOD tracks in jet cone, using trackrefs
  
  Int_t nTracks = jet->GetRefTracks()->GetEntriesFast();

  if(!nTracks) return;
  
  // standard cuts + ITS refit = stdrd PWG4 cut (compose them for productions with old ESDFilter task) 

  for (Int_t itrack=0; itrack<nTracks; itrack++) {

    AliAODTrack* track = (AliAODTrack*) jet->GetRefTracks()->At(itrack);

    UInt_t status = track->GetStatus();

    if(!track->TestFilterBit(fFilterBitJF)) continue; // track cut selection
    if(fRequireITSRefitJF &&  ((status&AliESDtrack::kITSrefit)==0)) continue;  

    if(TMath::Abs(track->Eta()) > fTrackEtaCut) continue; 
    if( track->Pt() < fTrackPtCutJF ) continue;
        
    list->Add(track);
  }

  //cout<<" trackrefs Size "<<nTracks<<" acc track list size "<<list->GetEntries()<<endl;

}


// ----------------------------------------------------------------------------------

void AliAnalysisTaskJetChem::FillReferenceFF(AliAODJet* jet){
  

  if(!jet) return;

  TList* jetTracks = new TList(); // FIXME - avoid new/delete
  //GetJetTracksResum(jetTracks,jet,0.7); 
  GetJetTracksTrackrefs(jetTracks,jet); 

  Double_t jetpt = jet->Pt();

  TIter next(jetTracks);
  while(AliAODTrack* track = static_cast<AliAODTrack*>(next())){
    
    Double_t trackpt = track->Pt();
    Double_t z       = trackpt/jetpt;
    Double_t xi      = TMath::Log(1/z);
    
    //cout<<" trackpt "<<trackpt<<" jetpt "<<jetpt<<" z "<<z<<" xi "<<xi<<endl;

    if(5<jetpt && jetpt<=10){
      fhdNdptTracksJetPt5to10->Fill(trackpt);
      fhdNdzJets5to10->Fill(z);
      fhdNdxiJets5to10->Fill(xi);
    }
    else if(10<jetpt && jetpt<=20){
      fhdNdptTracksJetPt10to20->Fill(trackpt);
      fhdNdzJets10to20->Fill(z);
      fhdNdxiJets10to20->Fill(xi);
    }
    else if(20<jetpt && jetpt<=30){
      fhdNdptTracksJetPt20to30->Fill(trackpt);
      fhdNdzJets20to30->Fill(z);
      fhdNdxiJets20to30->Fill(xi);
    }
    else if(30<jetpt && jetpt<=40){
      fhdNdptTracksJetPt30to40->Fill(trackpt);
      fhdNdzJets30to40->Fill(z);
      fhdNdxiJets30to40->Fill(xi);
    }
    else if(40<jetpt && jetpt<=60){
      fhdNdptTracksJetPt40to60->Fill(trackpt);
      fhdNdzJets40to60->Fill(z);
      fhdNdxiJets40to60->Fill(xi);
    }
  }

  
  delete jetTracks;
  
}

// -------------------------------------------------------------

void AliAnalysisTaskJetChem::FillReferencePlotsTracks(){

  // eta/phi & pt tracks before/after cuts
  // track multiplicity / evt

  Int_t nTracks     = fAOD->GetNTracks();
  Int_t countTracks = 0;

  // standard cuts + ITS refit = stdrd PWG4 cut (compose them for productions with old ESDFilter task) 

  for(Int_t itrack=0; itrack<nTracks; itrack++){
    
    AliAODTrack* track = fAOD->GetTrack(itrack);
    
    Double_t trackPt  = track->Pt();
    Double_t trackPhi = track->Phi();
    Double_t trackEta = track->Eta();

    fhPhiEtaTracksNoCut->Fill(trackEta,trackPhi);
    fhPtTracksNoCut->Fill(trackPt);
    
    UInt_t status = track->GetStatus();
    
    if(!track->TestFilterBit(fFilterBitJF)) continue; // track cut selection
    if(fRequireITSRefitJF &&  ((status&AliESDtrack::kITSrefit)==0)) continue;  

    fhPhiEtaTracks->Fill(trackEta,trackPhi);
    fhPtTracks->Fill(trackPt);

    if(TMath::Abs(track->Eta()) > fTrackEtaCut) continue; 
    if( track->Pt() < fTrackPtCutJF ) continue;

    countTracks++;
    
  }

  fhTrackMult->Fill(countTracks);
}

