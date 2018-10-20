#include "TChain.h"
#include "TH2D.h"
#include "TH3D.h"
#include "Riostream.h"
#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TList.h"
#include "AliAnalysisTaskMLTreeMaker.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "AliStack.h"
#include "AliESDtrackCuts.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliAODEvent.h"
#include "AliMCEvent.h"
#include "AliAODTrack.h"
#include "AliLog.h"
#include "AliMCEvent.h"
#include "AliAODMCParticle.h"
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliPIDResponse.h"
#include <AliAnalysisManager.h>
#include <AliAODHandler.h>
#include <AliESDInputHandler.h>
#include <AliAODInputHandler.h>
#include <AliMCEventHandler.h>
#include <AliMCEvent.h>
#include <AliMCParticle.h>
#include <AliAODMCParticle.h>
#include <AliAODMCHeader.h>
#include <AliStack.h>
#include <AliESDEvent.h>
#include <AliAODEvent.h>
#include <AliESDtrack.h>
#include <AliAODTrack.h>
#include <AliLog.h>
#include <TClonesArray.h>
#include <TParticle.h>
#include <TMCProcess.h>
#include <vector>
#include "AliPIDResponse.h"
#include "AliTrackReference.h"
#include "AliHeader.h"
#include "AliGenEventHeader.h"
#include "AliGenHijingEventHeader.h"
#include "AliGenCocktailEventHeader.h"
#include "AliCentrality.h"
#include "AliMultSelection.h"



// Authors: Sebastian Lehner (SMI Vienna) - selehner@cern.ch


ClassImp(AliAnalysisTaskMLTreeMaker)

Bool_t cutonTPCsignalN=  kFALSE;       
Int_t num= 0;
Int_t ev=0;


AliAnalysisTaskMLTreeMaker::AliAnalysisTaskMLTreeMaker():
  AliAnalysisTaskSE(),
  eventCuts(0),
  eventplaneCuts(0),
  evfilter(0),
  trcuts(0),
  trfilter(),
  pidcuts(0),      
  cuts(0),
  filter(0), 
  varManager(0), 
  fPIDResponse(0),        
  eta(0),
  phi(0),
  pt(0),        
  charge(0.),   
  NCrossedRowsTPC(0), 
  NClustersTPC(0),        
  HasSPDfirstHit(0),        
  RatioCrossedRowsFindableClusters(0), 
  NTPCSignal(0),
  fGeneratorHashs(0x0),      
  loCuts(kTRUE),        
  runn(0),      
  n(0),
  cent(0),
  fList(0x0),        
  fCentralityPercentileMin(0),
  fCentralityPercentileMax(100),         
  fPtMin(0),
  fPtMax(1000),
  fEtaMin(-10),
  fEtaMax(10),  
  fESigITSMin(-100.),
  fESigITSMax(3.),
  fESigTPCMin(-3.),
  fESigTPCMax(3.),
  fESigTOFMin(-3),
  fESigTOFMax(3),
  fPSigTPCMin(-100.),
  fPSigTPCMax(4.),
  fUsePionPIDTPC(kFALSE),
  fPionSigmas(kFALSE),
  fKaonSigmas(kFALSE),
  fFilterBit(96),
  gMultiplicity(-999),
  mcTrackIndex(0),
  fMcArray(0x0),
  mcEvent(0x0),      
  EsigTPC(0),
  EsigTOF(0),
  EsigITS(0),
  PsigTPC(0),
  PsigTOF(0),
  PsigITS(0),
  KsigTPC(0),
  KsigTOF(0),
  KsigITS(0),
  hasMC(kFALSE),       
  Rej(kFALSE),  
  MCpt(0),
  MCeta(0),
  MCphi(0),        
  MCvertx(0),
  MCverty(0),
  MCvertz(0),
  glabel(0),
  gLabelFirstMother(0),
  gLabelMinFirstMother(0),
  gLabelMaxFirstMother(0),
  iGenIndex(0),
  iPdgFirstMother(0),        
  dcar(),
  dcaz(), 
  vertx(0),
  verty(0),
  vertz(0),
  nITS(0),        
  nITSshared(0),        
  ITS1S(0),
  ITS2S(0),
  ITS3S(0),
  ITS4S(0),
  ITS5S(0),
  ITS6S(0),  
  chi2ITS(0),        
  chi2GlobalPerNDF(0),
  chi2GlobalvsTPC(0),
  fCutMaxChi2TPCConstrainedVsGlobalVertexType(0),
  pdg(0),
  pdgmother(0),        
  fmean(0),
  fwidth(0), 
  fuseCorr(kFALSE),      
  fTree(0),
  fQAHist(0)
{

}

AliAnalysisTaskMLTreeMaker::AliAnalysisTaskMLTreeMaker(const char *name) :
  AliAnalysisTaskSE(name),
  eventCuts(0),
  eventplaneCuts(0),
  evfilter(0),
  trcuts(0),
  trfilter(),
  pidcuts(0),      
  cuts(0),
  filter(0), 
  varManager(0), 
  fPIDResponse(0),        
  eta(0),
  phi(0),
  pt(0),        
  charge(0.),   
  NCrossedRowsTPC(0), 
  NClustersTPC(0),        
  HasSPDfirstHit(0),        
  RatioCrossedRowsFindableClusters(0), 
  NTPCSignal(0),
  fGeneratorHashs(0x0),        
  loCuts(kTRUE),        
  runn(0),      
  n(0),
  cent(0),
  fList(0x0), 
  fCentralityPercentileMin(0),
  fCentralityPercentileMax(100),         
  fPtMin(0),
  fPtMax(1000),
  fEtaMin(-10),
  fEtaMax(10),  
  fESigITSMin(-100.),
  fESigITSMax(3.),
  fESigTPCMin(-3.),
  fESigTPCMax(3.),
  fESigTOFMin(-3),
  fESigTOFMax(3),
  fPSigTPCMin(-100.),
  fPSigTPCMax(4.),
  fUsePionPIDTPC(kFALSE),
  fPionSigmas(kFALSE),
  fKaonSigmas(kFALSE),
  fFilterBit(96),
  gMultiplicity(-999),
  mcTrackIndex(0),
  fMcArray(0x0),
  mcEvent(0x0),      
  EsigTPC(0),
  EsigTOF(0),
  EsigITS(0),
  PsigTPC(0),
  PsigTOF(0),
  PsigITS(0),
  KsigTPC(0),
  KsigTOF(0),
  KsigITS(0),
  hasMC(kFALSE),       
  Rej(kFALSE),  
  MCpt(0),
  MCeta(0),
  MCphi(0),        
  MCvertx(0),
  MCverty(0),
  MCvertz(0), 
  glabel(0),
  gLabelFirstMother(0),
  gLabelMinFirstMother(0),
  gLabelMaxFirstMother(0),
  iGenIndex(0),
  iPdgFirstMother(0),         
  dcar(),
  dcaz(), 
  vertx(0),
  verty(0),
  vertz(0),
  nITS(0),        
  nITSshared(0),        
  ITS1S(0),
  ITS2S(0),
  ITS3S(0),
  ITS4S(0),
  ITS5S(0),
  ITS6S(0),  
  chi2ITS(0),        
  chi2GlobalPerNDF(0),
  chi2GlobalvsTPC(0),
  fCutMaxChi2TPCConstrainedVsGlobalVertexType(0),
  pdg(0),
  pdgmother(0),        
  hasmother(0),
  motherlabel(0),
  fmean(0),
  fwidth(0), 
  fuseCorr(kFALSE),              
  fTree(0),
  fQAHist(0)
{
// SetupTrackCuts(); 
// SetupEventCuts(); 
// AliInfo("Track & Event cuts were set"); 
   
 DefineOutput(1, TList::Class());
}

AliAnalysisTaskMLTreeMaker::~AliAnalysisTaskMLTreeMaker(){
  delete eventCuts;
  delete eventplaneCuts;
  delete evfilter;
  
  delete trcuts;
  delete trfilter;
  delete pidcuts;
  delete cuts;
  delete filter; 

  delete fList;
  delete fQAHist;
  delete fTree;


}

void AliAnalysisTaskMLTreeMaker::UserCreateOutputObjects() {
    
   AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
   AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
   inputHandler->SetNeedField();
   
    
  fList = new TList();
  fList->SetName("output_Tlist");
  fList->SetOwner();
   
  AliInfo("Try to get PIDResponse");
  fPIDResponse = inputHandler->GetPIDResponse();
  
     if (!fPIDResponse){
	   AliError("Failed to get PIDResponse - return");
	   return;}
  
  if(hasMC) std::cout <<"Running on MC!"<< std::endl;
  else std::cout <<"Running on RD!"<< std::endl;

  Bool_t oldStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);

  fTree = new TTree("Track_Tree","Tracks");
  fList->Add(fTree);
  
  fQAHist = new TH1F("h1", "h1 title", 4, 0, 1);
  fList->Add(fQAHist);
  
  fTree->Branch("centrality", &cent);
  fTree->Branch("#tracks", &n);
  fTree->Branch("pt", &pt);
  fTree->Branch("eta", &eta);
  fTree->Branch("phi", &phi);
  fTree->Branch("charge", &charge);
  fTree->Branch("RunNumber", &runn);
  fTree->Branch("EsigTPC", &EsigTPC);
  fTree->Branch("EsigITS", &EsigITS);
  fTree->Branch("EsigTOF", &EsigTOF);
  
  fTree->Branch("PsigTPC", &PsigTPC);
  
  fTree->Branch("NCrossedRowsTPC", &NCrossedRowsTPC);
  fTree->Branch("NClustersTPC", &NClustersTPC);
  fTree->Branch("RatioCrossedRowsFindableClusters", &RatioCrossedRowsFindableClusters);
  fTree->Branch("HasSPDfirstHit", &HasSPDfirstHit);   
  fTree->Branch("NTPCSignal", &NTPCSignal); 

  fTree->Branch("DCAxy", &dcar);
  fTree->Branch("DCAz", &dcaz);
  
  fTree->Branch("vertx", &vertx);
  fTree->Branch("verty", &verty);
  fTree->Branch("vertz", &vertz);
  
  fTree->Branch("nITS", &nITS);
  fTree->Branch("ITS1Shared", &ITS1S);
  fTree->Branch("ITS2Shared", &ITS2S); 
  fTree->Branch("ITS3Shared", &ITS3S); 
  fTree->Branch("ITS4Shared", &ITS4S); 
  fTree->Branch("ITS5Shared", &ITS5S); 
  fTree->Branch("ITS6Shared", &ITS6S);   
  fTree->Branch("nITSshared_frac", &nITSshared);
  fTree->Branch("chi2ITS", &chi2ITS);
//  fTree->Branch("chi2TPC", &chi2TPC);
  fTree->Branch("chi2GlobalvsTPC", &chi2GlobalvsTPC);
  fTree->Branch("chi2GlobalPerNDF", &chi2GlobalPerNDF);
  
  if(hasMC) {
      
    fTree->Branch("Pdg_Mother", &pdgmother);
    fTree->Branch("Mother_label", &motherlabel);
    fTree->Branch("Has_Mother", &hasmother);
    fTree->Branch("Pdg_Mother", &pdgmother);
  
    fTree->Branch("MCpt", &MCpt);
    fTree->Branch("MCeta", &MCeta);
    fTree->Branch("MCphi", &MCphi);
    
    fTree->Branch("MCTrack_vertx", &MCvertx);
    fTree->Branch("MCTrack_verty", &MCverty);
    fTree->Branch("MCTrack_vertz", &MCvertz);
    
    fTree->Branch("Pdg", &pdg);
    fTree->Branch("Label", &glabel);
    fTree->Branch("LabelFirstMother", &gLabelFirstMother);
    fTree->Branch("LabelMinFirstMother", &gLabelMinFirstMother);
    fTree->Branch("LabelMaxFirstMother", &gLabelMaxFirstMother);
    fTree->Branch("GenIndex", &iGenIndex);
    fTree->Branch("PdgFirstMother", &iPdgFirstMother);   
    
  }
  
  PostData(1, fList);
  
  AliInfo("Finished setting up the Output");
  TH1::AddDirectory(oldStatus);
  
  if(hasMC){
  TString generatorName = "Hijing_0;pizero_1;eta_2;etaprime_3;rho_4;omega_5;phi_6;jpsi_7;Pythia CC_8;Pythia BB_8;Pythia B_8";
 

  TObjArray arr = *(generatorName.Tokenize(";"));
  std::cout << "Used Generators: " << std::endl;
  for (int i = 0; i < arr.GetEntries(); ++i){
    TString temp = arr.At(i)->GetName();
    std::cout << "--- " << temp << std::endl;
    fGeneratorHashs.push_back(temp.Hash());
    }
   }
  }
//________________________________________________________________________

void AliAnalysisTaskMLTreeMaker::UserExec(Option_t *) {
  // Called for each event
  
  fQAHist->Fill("Events_all",1);
  
  AliVEvent* event = dynamic_cast<AliVEvent*>(InputEvent()); 
  
  if(!event) {
    AliError("event not available");
    fQAHist->Fill("Events_not_available",1);
    return;
  }

  UInt_t selectedMask=(1<<evfilter->GetCuts()->GetEntries())-1;
  varManager->SetEvent(event);
  if(selectedMask!=(evfilter->IsSelected(event))){
    fQAHist->Fill("Events_not_selected_filter",1);
    return;
  }
 
  
 AliInputEventHandler *eventHandler = nullptr;
 AliInputEventHandler *eventHandlerMC = nullptr;
  
  if ((AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler())->IsA() == AliAODInputHandler::Class()){
    eventHandler = dynamic_cast<AliAODInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
    eventHandlerMC = eventHandler;
  }
     
  if(hasMC){
    fMcArray = eventHandlerMC->MCEvent();
  }   
  
  if(fuseCorr){   
    AliDielectronPID::SetCentroidCorrFunction( (TH1*) fmean->Clone());
    AliDielectronPID::SetWidthCorrFunction( (TH1*) fwidth->Clone());
    ::Info("AliAnalysisTaskMLTreeMaker::UserExec","Setting Correction Histos");
  }

  
  AliMultSelection *MultSelection = 0x0; 
  MultSelection = (AliMultSelection * ) event->FindListObject("MultSelection");
  
  if( !MultSelection) {
   //If you get this warning (and lPercentiles 300) please check that the AliMultSelectionTask actually ran (before your task)
   AliWarning("AliMultSelection object not found!");
  }
  else cent = MultSelection->GetMultiplicityPercentile("V0M",kFALSE);
  
  fQAHist->Fill("Events before cent",1);
  
  if(cent<fCentralityPercentileMin || cent>fCentralityPercentileMax) return;
  
  fQAHist->Fill("Events after cent",1);
  
  Double_t lMultiplicityVar = -1;
  Int_t acceptedTracks = GetAcceptedTracks(event,lMultiplicityVar);
  
  runn = event->GetRunNumber();

  n= acceptedTracks;
  if(acceptedTracks){
  fTree->Fill();
  fQAHist->Fill("Events_track_and_cent_selected",1);
  }

  PostData(1, fList);
}

//~ //________________________________________________________________________

void  AliAnalysisTaskMLTreeMaker::FinishTaskOutput(){
  // Finish task output

  // not implemented ...


}
//~ 

//~ //________________________________________________________________________

void AliAnalysisTaskMLTreeMaker::Terminate(Option_t *) {

}
//~ 


//________________________________________________________________________

Double_t AliAnalysisTaskMLTreeMaker::IsEventAccepted(AliVEvent *event){

  if(event->GetPrimaryVertex()){
    if (TMath::Abs(event->GetPrimaryVertex()->GetZ()) < 10){
      if (event->GetPrimaryVertexSPD()->GetNContributors() > 0)
	return 1;
    }
  }
  return 0;
}


//________________________________________________________________________

Int_t AliAnalysisTaskMLTreeMaker::GetAcceptedTracks(AliVEvent *event, Double_t gCentrality){
  
  ev++;
  Int_t acceptedTracks = 0;
  Bool_t isAOD         = kFALSE;
  Int_t mpdg=0;
  eta.clear();
  phi.clear();
  pt.clear(); 
//  NClustersITS.clear();
  NCrossedRowsTPC.clear();
  NClustersTPC.clear();
  HasSPDfirstHit.clear(); 
  RatioCrossedRowsFindableClusters.clear(); 
  NTPCSignal.clear(); 
  EsigTPC.clear();
  EsigTOF.clear();
  EsigITS.clear();
  PsigTPC.clear();
  PsigTOF.clear();
  PsigITS.clear();
  KsigTPC.clear();
  KsigTOF.clear();
  KsigITS.clear();
  MCpt.clear();
  MCeta.clear();
  MCphi.clear(); 
  dcar.clear();
  dcaz.clear();
  nITS.clear();
  nITSshared.clear();
  chi2ITS.clear();
//  chi2TPC.clear();
//  chi2Global.clear();
  chi2GlobalvsTPC.clear();
  chi2GlobalPerNDF.clear();
  pdg.clear();
  pdgmother.clear();
  hasmother.clear();
  motherlabel.clear();
  label.clear();  
  charge.clear();
  MCvertx.clear();
  MCverty.clear();
  MCvertz.clear();
  glabel.clear();
  gLabelFirstMother.clear();
  gLabelMinFirstMother.clear();
  gLabelMaxFirstMother.clear();
  iGenIndex.clear();
  iPdgFirstMother.clear();
  ITS1S.clear();
  ITS2S.clear();
  ITS3S.clear();
  ITS4S.clear();
  ITS5S.clear();
  ITS6S.clear(); 
   
  
  
  // Loop over tracks in event
  AliMCEvent *mcEvent=0;
//  Int_t temppdg;
//  Int_t tempmpdg;
  AliAODMCParticle* mcMTrack;

  
  // need this to use PID in dielectron framework
  varManager->SetPIDResponse(fPIDResponse);


  for (Int_t iTracks = 0; iTracks < event->GetNumberOfTracks(); iTracks++) {
    
    
      fQAHist->Fill("All tracks",1); 
      AliVTrack* track = dynamic_cast<AliVTrack *>(event->GetTrack(iTracks));
      if (!track) {
	      AliError(Form("Could not receive track %d", iTracks));
	      continue;
      }

      // check for the first track if AOD or ESD track
      if (iTracks==0){
	if( ((TString)track->ClassName()).Contains("AliAODTrack") ){
	  isAOD = kTRUE;
	}
	else{
	  isAOD = kFALSE;
	}
      }
      if(!isAOD) fQAHist->Fill("Not AOD track",1);
      else       fQAHist->Fill("Is AOD track",1);
      
      fQAHist->Fill("After ESD check, bef. MC",1); 
      
      if(hasMC){ 
        mcEvent = MCEvent(); 
        if (!mcEvent) {
          AliError(Form("Could not receive MC -> hasMC set to kFALSE!!"));
          hasMC=kFALSE;
          continue;
        }
        else{
          fQAHist->Fill("After MC check",1); 
//          Rej=kFALSE;
//          if(!(CheckGenerator(TMath::Abs(track->GetLabel())))) Rej=kTRUE;
          }

        }

      fQAHist->Fill("Tracks aft MC, bef tr cuts",1); 
      
    UInt_t selectedMask = (1 << filter->GetCuts()->GetEntries()) - 1;
    if (selectedMask != (filter->IsSelected((AliVParticle*) track))) {
      fQAHist->Fill("Tracks not selected filter",1);
          continue;
      }
      
      Double_t pttemp = track->Pt();
      Double_t etatemp = track->Eta();
       
      //Get PID response for tree - this is w/o postcalibration - the PID response after postcalibration has to be taken from dielectron task
      Double_t tempEsigTPC=fPIDResponse->NumberOfSigmasTPC(track, (AliPID::EParticleType) 0);

        if (fuseCorr){
            tempEsigTPC-=AliDielectronPID::GetCntrdCorr(track);
            tempEsigTPC/=AliDielectronPID::GetWdthCorr(track);
        }

      Double_t tempEsigITS=fPIDResponse->NumberOfSigmasITS(track, (AliPID::EParticleType) 0);
      Double_t tempEsigTOF=fPIDResponse->NumberOfSigmasTOF(track, (AliPID::EParticleType) 0);
      
      fQAHist->Fill("Selected tracks",1); 

      if(hasMC){ 
        AliAODMCParticle* mcTrack = dynamic_cast<AliAODMCParticle *>(mcEvent->GetTrack(TMath::Abs(track->GetLabel())));

        pdg.push_back( mcTrack->PdgCode());
        
        MCpt.push_back(mcTrack->Pt());
        MCeta.push_back(mcTrack->Eta());
        MCphi.push_back(mcTrack->Phi());
        

        //Get production vertex for MC tracks

          Double_t MCvert[3] = {0};
          mcTrack->XvYvZv(MCvert);
          
          label.push_back(track->GetLabel());
          
          MCvertx.push_back(MCvert[0]);
          MCverty.push_back(MCvert[1]);
          MCvertz.push_back(MCvert[2]);

        if(!(mcTrack->GetMother() < 0)) {  
          hasmother.push_back(1);
          AliAODMCParticle* mcmother = dynamic_cast<AliAODMCParticle *>(fMCEvent->GetTrack(mcTrack->GetMother()));
                mpdg= mcmother->PdgCode();
	        pdgmother.push_back(mpdg);

          motherlabel.push_back(abs(mcmother->GetLabel()));
        }
        else{
          hasmother.push_back(0);  
          pdgmother.push_back( -9999);
          motherlabel.push_back(-9999);
        }

      if( abs(mcTrack->PdgCode())==11 && mpdg!=22) fQAHist->Fill("Selected electrons, non-conversion",1);
          
     // infos of first mother:
      Int_t gMotherIndex = mcTrack->GetMother();
      Int_t tempFirstMotherIndex    = 666666666;
      Int_t tempLabelFirstMother=-1;
      Int_t tempPdgFirstMother=-99;
      Int_t tempLabelMinFirstMother=-1;
      Int_t tempLabelMaxFirstMother=-1;
      Int_t nParticles = mcEvent->GetNumberOfTracks();

      AliMCParticle* firstMotherTrack = NULL;
      
      if(gMotherIndex != -1) {
	
  	AliMCParticle* motherTrack = (AliMCParticle*)(mcEvent->GetTrack(gMotherIndex));
  	Int_t temppdgmother = motherTrack->PdgCode();

  	// find first mother
  	tempFirstMotherIndex = motherTrack->GetMother();

  	while(tempFirstMotherIndex>0){
  	  tempLabelFirstMother = tempFirstMotherIndex;
  	  firstMotherTrack = (AliMCParticle*)(mcEvent->GetTrack(tempLabelFirstMother));
  	  tempFirstMotherIndex = firstMotherTrack->GetMother();
  	}

  	if(tempLabelFirstMother != -1) { 	  // if grandmother not primary!
  	  tempPdgFirstMother = firstMotherTrack->PdgCode();
  	}
  	else{     // if grandmother already primary!
  	  tempLabelFirstMother = gMotherIndex; // set mother to first mother
  	  tempPdgFirstMother = temppdgmother;
  	}

  	// find range of -1 - minimum
  	tempLabelMinFirstMother = tempLabelFirstMother;

  	while(tempFirstMotherIndex<0){	
  	  tempLabelMinFirstMother--;
  	  if(tempLabelMinFirstMother<0){
  	    tempFirstMotherIndex = 0;
  	  }
  	  else{
  	    firstMotherTrack = (AliMCParticle*)(mcEvent->GetTrack(tempLabelMinFirstMother));
  	    tempFirstMotherIndex = firstMotherTrack->GetMother();
  	  }
  	}
  	tempLabelMinFirstMother ++; // set back by one
  	tempFirstMotherIndex = -1; // set back to -1

  	// find range of -1 - maximum
  	tempLabelMaxFirstMother = tempLabelFirstMother;
  	while(tempFirstMotherIndex<0){
  	  tempLabelMaxFirstMother++;
  	  if(tempLabelMaxFirstMother > nParticles){
  	    tempFirstMotherIndex = 0;
  	  }
  	  else{
  	    firstMotherTrack = (AliMCParticle*)(mcEvent->GetTrack(tempLabelMaxFirstMother));
  	    tempFirstMotherIndex = firstMotherTrack->GetMother();
  	  }

  	}
  	tempLabelMaxFirstMother --; // set back by one     
          
      }
      
      glabel.push_back(mcTrack->GetLabel());
      gLabelFirstMother.push_back(tempLabelFirstMother);
      gLabelMinFirstMother.push_back(tempLabelMinFirstMother);
      gLabelMaxFirstMother.push_back(tempLabelMaxFirstMother);
      iGenIndex.push_back(CheckGenerator(TMath::Abs(track->GetLabel())));  
      iPdgFirstMother.push_back(tempPdgFirstMother);
               
      } //End if hasMC 
      
                        
       if(!acceptedTracks){   //Get vertex only for first track in event
          Double_t vert[3] = {0};  
          event->GetPrimaryVertex()->GetXYZ(vert); 
          vertx= vert[0];
          verty= vert[1];
          vertz= vert[2];
        }


      //Fill Tree with non MC data
      EsigTPC.push_back(tempEsigTPC);
      EsigITS.push_back(tempEsigITS);
      EsigTOF.push_back(tempEsigTOF);
      
      Double_t tempPsigTPC=fPIDResponse->NumberOfSigmasTPC(track, (AliPID::EParticleType) 2);
      PsigTPC.push_back(tempPsigTPC);
      
//      if(fPionSigmas){
//        Double_t tempPsigTPC=fPIDResponse->NumberOfSigmasTPC(track, (AliPID::EParticleType) 2);
//        Double_t tempPsigITS=fPIDResponse->NumberOfSigmasITS(track, (AliPID::EParticleType) 2);
//        Double_t tempPsigTOF=fPIDResponse->NumberOfSigmasTOF(track, (AliPID::EParticleType) 2);
//        PsigTPC.push_back(tempPsigTPC);
//        PsigITS.push_back(tempPsigITS);
//        PsigTOF.push_back(tempPsigTOF);
//      }
//      if(fKaonSigmas){
//        Double_t tempKsigTPC=fPIDResponse->NumberOfSigmasTPC(track, (AliPID::EParticleType) 3);
//        Double_t tempKsigITS=fPIDResponse->NumberOfSigmasITS(track, (AliPID::EParticleType) 3);
//        Double_t tempKsigTOF=fPIDResponse->NumberOfSigmasTOF(track, (AliPID::EParticleType) 3);
//        KsigTPC.push_back(tempKsigTPC);
//        KsigITS.push_back(tempKsigITS);
//        KsigTOF.push_back(tempKsigTOF);
//      }
      eta.push_back(etatemp);
      phi.push_back(track->Phi());
      pt.push_back(pttemp);
      charge.push_back(track->Charge());   

      NCrossedRowsTPC.push_back(track->GetTPCCrossedRows());
      NClustersTPC.push_back(track->GetTPCNcls());//->GetNumberOfTPCClusters());
      HasSPDfirstHit.push_back(track->HasPointOnITSLayer(0)); 
      RatioCrossedRowsFindableClusters.push_back((Double_t) track->GetTPCCrossedRows()/ (Double_t) track->GetTPCNclsF());       
      NTPCSignal.push_back(track->GetTPCsignalN());
      

      if( track->HasSharedPointOnITSLayer(0) ) {ITS1S.push_back(1);}
      else {ITS1S.push_back(0);}
      if( track->HasSharedPointOnITSLayer(1) ) {ITS2S.push_back(1);}
      else {ITS2S.push_back(0);}
      if( track->HasSharedPointOnITSLayer(2) ) {ITS3S.push_back(1);}
      else {ITS3S.push_back(0);}
      if( track->HasSharedPointOnITSLayer(3) ) {ITS4S.push_back(1);}
      else {ITS4S.push_back(0);}
      if( track->HasSharedPointOnITSLayer(4) ) {ITS5S.push_back(1);}
      else {ITS5S.push_back(0);}
      if( track->HasSharedPointOnITSLayer(5) ) {ITS6S.push_back(1);}
      else {ITS6S.push_back(0);}
      
       //Get DCA position
      if(isAOD){
	Double_t tempdcaD[2] = {0.,0.};
      	GetDCA(event,(AliAODTrack*)track,tempdcaD,0);
	dcar.push_back(tempdcaD[0]);
	dcaz.push_back(tempdcaD[1]);
      }
      else{
	Float_t tempdca[2] = {0.,0.};
      	track->GetImpactParameters( &tempdca[0], &tempdca[1]); //GetImpactParameter is also used in AliESDtrackCuts.cxx to cut on DCA to vertex
	dcar.push_back(tempdca[0]);
	dcaz.push_back(tempdca[1]);
      }

      Int_t tempnits = track->GetNcls(0);    // 0 = ITS 
      nITS.push_back(tempnits);        
      Double_t nitssharedtemp = 0.;
 
      if(tempnits){
        for(int d = 0; d<6;d++){
          nitssharedtemp+= (Double_t) track->HasSharedPointOnITSLayer(d);
        }
        nitssharedtemp/=tempnits;
      }

      nITSshared.push_back(nitssharedtemp);
      
      chi2ITS.push_back(track->GetITSchi2());
//      chi2TPC.push_back(track->GetTPCchi2());//this variable will be always 0 for AODs (not yet in)
      
      if(isAOD){ chi2GlobalPerNDF.push_back(((AliAODTrack*)track)->Chi2perNDF());
                 chi2GlobalvsTPC.push_back(((AliAODTrack*)track)->GetChi2TPCConstrainedVsGlobal());  
      }

      acceptedTracks++;
  }
    
  num++;
  return acceptedTracks;  
}

//______________________________________________________________________________
Bool_t AliAnalysisTaskMLTreeMaker::GetDCA(const AliVEvent* event, const AliAODTrack* track, Double_t* d0z0, Double_t* covd0z0)
// this is a copy of the AliDielectronVarManager
{
  if(track->TestBit(AliAODTrack::kIsDCA)){
    d0z0[0]=track->DCA();
    d0z0[1]=track->ZAtDCA();
    // the covariance matrix is not stored in case of AliAODTrack::kIsDCA
    return kTRUE;
  }

  Bool_t ok=kFALSE;
  if(event) {
    AliExternalTrackParam etp; etp.CopyFromVTrack(track);

    Float_t xstart = etp.GetX();
    if(xstart>3.) {
      d0z0[0]=-999.;
      d0z0[1]=-999.;
      return kFALSE;
    }

    AliAODVertex *vtx =(AliAODVertex*)(event->GetPrimaryVertex());
    Double_t fBzkG = event->GetMagneticField(); // z componenent of field in kG
    ok = etp.PropagateToDCA(vtx,fBzkG,kVeryBig,d0z0,covd0z0);
  }
  if(!ok){
    d0z0[0]=-999.;
    d0z0[1]=-999.;
  }
  return ok;
}


void AliAnalysisTaskMLTreeMaker::SetupTrackCuts(AliDielectronCutGroup* f)
{
filter   = new AliAnalysisFilter("filter","filter");  
filter->AddCuts(f);
}


void AliAnalysisTaskMLTreeMaker::SetupEventCuts(AliDielectronEventCuts* f)
{
  evfilter   = new AliAnalysisFilter("evfilter","evfilter");  
  evfilter->AddCuts(f);
}



int AliAnalysisTaskMLTreeMaker::CheckGenerator(Int_t trackID){     //check if the generator is on the list of generators
  if (fGeneratorHashs.size() == 0) return -1;
  TString genname;
  Bool_t hasGenerator = fMcArray->GetCocktailGenerator(TMath::Abs(trackID), genname); // fMC is AliMCEvent
  
    if(!hasGenerator) {

    Printf("no cocktail header list was found for this track");
    return -2;
  }
  else{

    for ( int i = 0; i < fGeneratorHashs.size(); ++i){
      // std::cout << genname.Hash() << " " << fGeneratorHashs[i] << std::endl;
      if (genname.Hash() == fGeneratorHashs[i]) return i;

    }
    std::cout << genname << std::endl;
    return -3;
  }
  return -4; // should not happen
}
