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
#include "AliAnalysisTaskMLTreeMakerwoCutlib.h"
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


ClassImp(AliAnalysisTaskMLTreeMakerwoCutlib)

Bool_t cutonTPCsignalN=  kFALSE;       
Int_t num= 0;
Int_t ev=0;


AliAnalysisTaskMLTreeMakerwoCutlib::AliAnalysisTaskMLTreeMakerwoCutlib():
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
  enh(0),
  NCrossedRowsTPC(0), 
  NClustersTPC(0),        
  HasSPDfirstHit(0),        
  RatioCrossedRowsFindableClusters(0), 
  NTPCSignal(0),
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
  label(0),      
  motherlabel(0),
  fTree(0),
  fQAHist(0)
{

}

AliAnalysisTaskMLTreeMakerwoCutlib::AliAnalysisTaskMLTreeMakerwoCutlib(const char *name) :
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
  enh(0),
  NCrossedRowsTPC(0), 
  NClustersTPC(0),        
  HasSPDfirstHit(0),        
  RatioCrossedRowsFindableClusters(0), 
  NTPCSignal(0),
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
  label(0),      
  motherlabel(0),
  fTree(0),
  fQAHist(0)
{
 SetupTrackCuts(); 
 SetupEventCuts(); 
 AliInfo("Track & Event cuts were set"); 
   
 DefineOutput(1, TList::Class());
}

AliAnalysisTaskMLTreeMakerwoCutlib::~AliAnalysisTaskMLTreeMakerwoCutlib(){
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

void AliAnalysisTaskMLTreeMakerwoCutlib::UserCreateOutputObjects() {
    
   AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
   AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
   inputHandler->SetNeedField();
   
    
  fList = new TList();
  fList->SetName("output_Tlist");
  fList->SetOwner();
   
  fPIDResponse = inputHandler->GetPIDResponse();
     if (!fPIDResponse){
	   
	   return;}
  
  if (man->GetMCtruthEventHandler()!=0x0) hasMC=kTRUE;
  else hasMC = kFALSE; 


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
      
    fTree->Branch("Pdg", &pdg);
    fTree->Branch("Pdg_Mother", &pdgmother);
    fTree->Branch("Mother_label", &motherlabel);
    fTree->Branch("Label", &label);      
    fTree->Branch("Has_Mother", &hasmother);
    fTree->Branch("IsEnh", &enh);
  
    fTree->Branch("MCpt", &MCpt);
    fTree->Branch("MCeta", &MCeta);
    fTree->Branch("MCphi", &MCphi);
    
    fTree->Branch("MCTrack_vertx", &MCvertx);
    fTree->Branch("MCTrack_verty", &MCverty);
    fTree->Branch("MCTrack_vertz", &MCvertz);
  }
  
  PostData(1, fList);
  
  AliInfo("Finished setting up the Output");
  TH1::AddDirectory(oldStatus);
}

//________________________________________________________________________

void AliAnalysisTaskMLTreeMakerwoCutlib::UserExec(Option_t *) {
  // Main loop

  // Called for each event
  AliVEvent* event = dynamic_cast<AliVEvent*>(InputEvent()); 
  
  fQAHist->Fill("Events_all",1);
  
  if(!event) {
    AliError("event not available");
    return;
  }

  UInt_t selectedMask=(1<<evfilter->GetCuts()->GetEntries())-1;
  varManager->SetEvent(event);
  if(selectedMask!=(evfilter->IsSelected(event))){
    return;
  }
  
  fQAHist->Fill("Events_accepted",1);
  
  if(hasMC){
    AliMCEventHandler* mchandler = dynamic_cast<AliMCEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
    fMcArray = mchandler->MCEvent();
    // get the accepted tracks in main event
  }
  Double_t lMultiplicityVar = -1;
  Int_t acceptedTracks = GetAcceptedTracks(event,lMultiplicityVar);

  AliMultSelection *MultSelection = 0x0; 
  MultSelection = (AliMultSelection * ) event->FindListObject("MultSelection");
  
  if( !MultSelection) {
   //If you get this warning (and lPercentiles 300) please check that the AliMultSelectionTask actually ran (before your task)
   AliWarning("AliMultSelection object not found!");
  }
  
  else cent = MultSelection->GetMultiplicityPercentile("V0M");
  
  
 runn = event->GetRunNumber();

  
  n= acceptedTracks;
  if(acceptedTracks){
    fTree->Fill();
    fQAHist->Fill("Events_track_selected",1);
  }

  PostData(1, fList);
}

//~ //________________________________________________________________________

void  AliAnalysisTaskMLTreeMakerwoCutlib::FinishTaskOutput(){
  // Finish task output

  // not implemented ...


}
//~ 

//~ //________________________________________________________________________

void AliAnalysisTaskMLTreeMakerwoCutlib::Terminate(Option_t *) {

}
//~ 


//________________________________________________________________________

Double_t AliAnalysisTaskMLTreeMakerwoCutlib::IsEventAccepted(AliVEvent *event){

  if(event->GetPrimaryVertex()){
    if (TMath::Abs(event->GetPrimaryVertex()->GetZ()) < 10){
      if (event->GetPrimaryVertexSPD()->GetNContributors() > 0)
	return 1;
    }
  }
  return 0;
}


//________________________________________________________________________

Int_t AliAnalysisTaskMLTreeMakerwoCutlib::GetAcceptedTracks(AliVEvent *event, Double_t gCentrality){
  
  ev++;
  Int_t acceptedTracks = 0;
  Bool_t isAOD         = kFALSE;
  
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
  enh.clear();
  MCvertx.clear();
  MCverty.clear();
  MCvertz.clear();
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
      AliVTrack* track = dynamic_cast<AliVTrack *>(event->GetTrack(iTracks));
      if (!track) {
	      AliError(Form("Could not receive ESD track %d", iTracks));
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
      
      fQAHist->Fill("After ESD check, bef. MC",1); 
      
      if(hasMC){ 
        mcEvent = MCEvent(); 
        if (!mcEvent) {
          AliError(Form("Could not receive MC -> hasMC set to kFALSE!!"));
          hasMC=kFALSE;
          continue;
        }
        else{
          fQAHist->Fill("After MC check, bef. Hij",1); 

          Rej=kFALSE;

                mcMTrack = dynamic_cast<AliAODMCParticle *>(mcEvent->GetTrack(TMath::Abs(track->GetLabel())));
//                temppdg = mcMTrack->PdgCode(); 
                
                if(!(mcMTrack->GetMother() < 0)){       //get direct mother
                    mcTrackIndex = mcMTrack->GetMother(); 
                    mcMTrack = dynamic_cast<AliAODMCParticle *>(mcEvent->GetTrack(mcMTrack->GetMother()));
//                    tempmpdg= mcMTrack->PdgCode(); 
                }
//                else tempmpdg=-9999;
                    
                while(!(mcMTrack->GetMother() < 0)){        //get first mother in chain
                    mcTrackIndex = mcMTrack->GetMother(); 
                    mcMTrack = dynamic_cast<AliAODMCParticle *>(mcEvent->GetTrack(mcMTrack->GetMother()));
                }

                if(!(mcEvent->IsFromBGEvent(abs(mcTrackIndex)))) Rej=kTRUE;

          }


        }

      fQAHist->Fill("Tracks aft MC&Hij, bef tr cuts",1); 

      UInt_t selectedMask=(1<<filter->GetCuts()->GetEntries())-1;
      if(selectedMask!=(filter->IsSelected((AliVParticle*)track))){
          continue;
      }
      
//      // Kinematic cuts
      Double_t pttemp = track->Pt();
      Double_t etatemp = track->Eta();
//      
//      if( pttemp > fPtMax || pttemp < fPtMin ) continue;
//      if( etatemp > fEtaMax || etatemp < fEtaMin ) continue;
// 
      Double_t tempEsigTPC=fPIDResponse->NumberOfSigmasTPC(track, (AliPID::EParticleType) 0);
      Double_t tempEsigITS=fPIDResponse->NumberOfSigmasITS(track, (AliPID::EParticleType) 0);
      Double_t tempEsigTOF=fPIDResponse->NumberOfSigmasTOF(track, (AliPID::EParticleType) 0);
//      
//      if(fUsePionPIDTPC){
//        if (fPIDResponse->NumberOfSigmasTPC(track, (AliPID::EParticleType) 2) > fPSigTPCMin &&  fPIDResponse->NumberOfSigmasTPC(track, (AliPID::EParticleType) 2)  < fPSigTPCMax){ continue;} //exclude pions in TPC
//      }
//
//      if (fPIDResponse->CheckPIDStatus(AliPIDResponse::kTOF,track)==AliPIDResponse::kDetPidOk && (tempEsigTOF < fESigTOFMin || tempEsigTOF > fESigTOFMax)) continue;  
//      if (tempEsigITS < fESigITSMin || tempEsigITS > fESigITSMax) continue;  
//      if (tempEsigTPC < fESigTPCMin || tempEsigTPC > fESigTPCMax) continue;
      
      
      fQAHist->Fill("Selected tracks",1); 

//      printf("Found %d with Mother %d originated from %d generated by %s  \n",temppdg,tempmpdg,mcMTrack->PdgCode(),(mcEvent->GetGenerator(mcTrackIndex)).Data()); 
      //Fill Tree with MC data
      if(hasMC){ 
        AliAODMCParticle* mcTrack = dynamic_cast<AliAODMCParticle *>(mcEvent->GetTrack(TMath::Abs(track->GetLabel())));

        pdg.push_back( mcTrack->PdgCode());
        if(Rej) enh.push_back(1);
        else enh.push_back(0);
        
        
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
	        pdgmother.push_back( mcmother->PdgCode());

          motherlabel.push_back(abs(mcmother->GetLabel()));
        }
        else{
          hasmother.push_back(0);  
          pdgmother.push_back( -9999);
          motherlabel.push_back(-9999);
        }
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
Bool_t AliAnalysisTaskMLTreeMakerwoCutlib::GetDCA(const AliVEvent* event, const AliAODTrack* track, Double_t* d0z0, Double_t* covd0z0)
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


void AliAnalysisTaskMLTreeMakerwoCutlib::SetupTrackCuts()
{
  
filter   = new AliAnalysisFilter("filter","filter");
varManager = new AliDielectronVarManager;  

cuts     = new AliDielectronCutGroup("cuts","cuts",AliDielectronCutGroup::kCompAND);    

AliDielectronVarCuts* trackCutsAOD =new AliDielectronVarCuts("trackCutsAOD","trackCutsAOD");
AliDielectronPID *PIDcut_3 = new AliDielectronPID("PIDcut_3","PIDcut_3");
AliDielectronTrackCuts *trackCutsDiel = new AliDielectronTrackCuts("trackCutsDiel","trackCutsDiel");


//Cuts used for nano AOD filtering (taken from ConfigLMEE_nano_PbPb2015.C on 22082018) - redundant for nano AODs but necessary for MC

//  AliDielectronTrackCuts *trkFilter = new AliDielectronTrackCuts("TrkFilter","TrkFilter");    
//  trkFilter->SetAODFilterBit(AliDielectronTrackCuts::kTPCqualSPDany); // I think we loose the possibility to use prefilter?

  AliDielectronVarCuts *varCuts   = new AliDielectronVarCuts("VarCuts","VarCuts");
  AliDielectronTrackCuts *trkCuts = new AliDielectronTrackCuts("TrkCuts","TrkCuts");
  // specific cuts
  trkCuts->SetITSclusterCut(AliDielectronTrackCuts::kOneOf, 3); // SPD any
  trkCuts->SetRequireITSRefit(kTRUE);
  trkCuts->SetRequireTPCRefit(kTRUE); // not useful when using prefilter

  // standard cuts
  varCuts->AddCut(AliDielectronVarManager::kNclsTPC,      80.0, 160.0);
  varCuts->AddCut(AliDielectronVarManager::kNclsITS,      3.0, 100.0);
  varCuts->AddCut(AliDielectronVarManager::kITSchi2Cl,    0.0,   15.0);
//  varCuts->AddCut(AliDielectronVarManager::kNclsSITS,     0.0,   3.1); // means 0 and 1 shared Cluster    // did not work on ESD when filtering nanoAODs
  varCuts->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0,   8.0);
  varCuts->AddCut(AliDielectronVarManager::kNFclsTPCr,    80.0, 160.0);

  varCuts->AddCut(AliDielectronVarManager::kPt,           0.2, 8.);
  varCuts->AddCut(AliDielectronVarManager::kEta,         -0.8,   0.8);
  varCuts->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
  varCuts->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);
  varCuts->AddCut(AliDielectronVarManager::kKinkIndex0,   0.);


  AliDielectronPID *pidCuts        = new AliDielectronPID("PIDCuts","PIDCuts");
  // TOF
  // pidCuts->AddCut(AliDielectronPID::kTOF,AliPID::kElectron,-4.,4.,0.,0.,kFALSE, AliDielectronPID::kIfAvailable);   //tighter TOF if cut applied below
  // TPC
  pidCuts->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,-4.,4.);
  pidCuts->AddCut(AliDielectronPID::kTPC,AliPID::kPion,-100.,3.5,0.,0.,kTRUE);
  // ITS
  pidCuts->AddCut(AliDielectronPID::kITS,AliPID::kElectron,-4.,4.);


//Define Carsten's Cut set 3 (as specified in e mail from Carsten from 10.01.2018) - if cuts from NanoAOD were tighter they were implemented here (and are commented above)

PIDcut_3->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -1.5, 3.0 , 0. ,100., kFALSE);
PIDcut_3->AddCut(AliDielectronPID::kTPC,AliPID::kPion, -100, 4.0 , 0. ,100., kTRUE);
PIDcut_3->AddCut(AliDielectronPID::kITS,AliPID::kElectron, -3.0, 1.0 , 0. ,100., kFALSE);
PIDcut_3->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3.0 , 3.0 , 0. ,100., kFALSE, AliDielectronPID::kIfAvailable);

trackCutsAOD->AddCut(AliDielectronVarManager::kPt,           0.2, 8.0);
trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);
trackCutsAOD->AddCut(AliDielectronVarManager::kNclsITS,      5.0, 100.0);
trackCutsAOD->AddCut(AliDielectronVarManager::kITSchi2Cl,    0.0,   5.0);
trackCutsAOD->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0,   4.0);
trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCr,    120.0, 160.0);
trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCfCross,     0.95, 1.05);

AliDielectronCutGroup* SharedClusterCut = new AliDielectronCutGroup("SharedClusterCut","SharedClusterCut",AliDielectronCutGroup::kCompOR);
double delta = 0.00001;
AliDielectronVarCuts* trackCutsSharedCluster0 = new AliDielectronVarCuts("trackCutsSharedCluster0", "trackCutsSharedCluster0");
trackCutsSharedCluster0->AddCut(AliDielectronVarManager::kNclsSMapITS, 0-delta, 0+delta);
AliDielectronVarCuts* trackCutsSharedCluster2 = new AliDielectronVarCuts("trackCutsSharedCluster2", "trackCutsSharedCluster2");
trackCutsSharedCluster2->AddCut(AliDielectronVarManager::kNclsSMapITS, 2-delta, 2+delta);
AliDielectronVarCuts* trackCutsSharedCluster4 = new AliDielectronVarCuts("trackCutsSharedCluster4", "trackCutsSharedCluster4");
trackCutsSharedCluster4->AddCut(AliDielectronVarManager::kNclsSMapITS, 4-delta, 4+delta);
AliDielectronVarCuts* trackCutsSharedCluster8 = new AliDielectronVarCuts("trackCutsSharedCluster8", "trackCutsSharedCluster8");
trackCutsSharedCluster8->AddCut(AliDielectronVarManager::kNclsSMapITS, 8-delta, 8+delta);
AliDielectronVarCuts* trackCutsSharedCluster16 = new AliDielectronVarCuts("trackCutsSharedCluster16", "trackCutsSharedCluster16");
trackCutsSharedCluster16->AddCut(AliDielectronVarManager::kNclsSMapITS, 16-delta, 16+delta);
AliDielectronVarCuts* trackCutsSharedCluster32 = new AliDielectronVarCuts("trackCutsSharedCluster32", "trackCutsSharedCluster32");
trackCutsSharedCluster32->AddCut(AliDielectronVarManager::kNclsSMapITS, 32-delta, 32+delta);
SharedClusterCut->AddCut(trackCutsSharedCluster0);
SharedClusterCut->AddCut(trackCutsSharedCluster2);
SharedClusterCut->AddCut(trackCutsSharedCluster4);
SharedClusterCut->AddCut(trackCutsSharedCluster8);
SharedClusterCut->AddCut(trackCutsSharedCluster16);
SharedClusterCut->AddCut(trackCutsSharedCluster32);


trackCutsDiel->SetAODFilterBit(AliDielectronTrackCuts::kGlobalNoDCA);//(1<<4) -> error
trackCutsDiel->SetClusterRequirementITS(AliDielectronTrackCuts::Detector(0),AliDielectronTrackCuts::ITSClusterRequirement(3));//(AliESDtrackCuts::kSPD,AliESDtrackCuts::kFirst) -> error

//Add desired cuts to cutgroup//

//NanoAOD filter cuts

//cuts->AddCut(trkFilter);    //not used on ESDs in ConfigLMEE_nano_PbPb2015.C

cuts->AddCut(varCuts);
cuts->AddCut(trkCuts);
cuts->AddCut(pidCuts);


//Carsten cut set 3
cuts->AddCut(PIDcut_3);
cuts->AddCut(trackCutsDiel);
cuts->AddCut(trackCutsAOD);
cuts->AddCut(SharedClusterCut);

cuts->Print();

filter->AddCuts(cuts);
}


void AliAnalysisTaskMLTreeMakerwoCutlib::SetupEventCuts()
{
  evfilter   = new AliAnalysisFilter("evfilter","evfilter");  
   
  //Event cuts as in NanoAOD filtering (as from ConfigLMEE_nano_PbPb2015.C on 22082018)
  
  AliDielectronEventCuts *eventCuts=new AliDielectronEventCuts("eventCuts","eventCuts");
  eventCuts->SetVertexType(AliDielectronEventCuts::kVtxAny);
  eventCuts->SetRequireVertex();
  eventCuts->SetMinVtxContributors(1);
  eventCuts->SetVertexZ(-10.,+10.);
  // eventCuts->SetCentralityRange(10., 50., kTRUE);
  eventCuts->Print();
  evfilter->AddCuts(eventCuts);

  AliDielectronVarCuts *eventplaneCuts = new AliDielectronVarCuts("eventplaneCuts","eventplaneCuts");
  eventplaneCuts->AddCut(AliDielectronVarManager::kQnTPCrpH2,-999.,kTRUE); // makes sure that the event has an eventplane
  eventplaneCuts->Print();

  evfilter->AddCuts(eventplaneCuts);

}
