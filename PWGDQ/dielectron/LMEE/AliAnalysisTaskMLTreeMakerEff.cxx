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
#include "AliAnalysisTaskMLTreeMakerEff.h"
#include "AliTrackReference.h"
#include "AliHeader.h"
#include "AliGenEventHeader.h"
#include "AliGenHijingEventHeader.h"
#include "AliGenCocktailEventHeader.h"
#include "AliCentrality.h"
#include "AliMultSelection.h"

// Authors: Sebastian Lehner (SMI Vienna) - selehner@cern.ch


ClassImp(AliAnalysisTaskMLTreeMakerEff)

   
  Int_t num1= 0;
  Int_t ev1=0;


AliAnalysisTaskMLTreeMakerEff::AliAnalysisTaskMLTreeMakerEff():
  AliAnalysisTaskSE(),
  fPIDResponse(0),
  NCrossedRowsTPC(0),
  NClustersTPC(0),
  HasSPDfirstHit(0), 
  RatioCrossedRowsFindableClusters(0), 
  NTPCSignal(0),
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
  fESDTrackCuts(0),
  EsigTPC(0),
  EsigTOF(0),
  EsigITS(0),
  PsigTPC(0),
  PsigTOF(0),
  PsigITS(0),
  KsigTPC(0),
  KsigTOF(0),
  KsigITS(0),
  nITS(0),
  nITSshared(0),
  chi2ITS(0),
  chi2TPC(0),
  chi2Global(0),
  chi2GlobalvsTPC(0),
  fCutMaxChi2TPCConstrainedVsGlobalVertexType(0),
  ProdVx(0),
  ProdVy(0),
  ProdVz(0),
  dcar(),
  dcaz(),
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
  mcTrackIndex(0),
  fMcArray(0x0),
  MCpt(0),
  MCeta(0),
  MCphi(0),
  pt(0),      
  eta(0),
  phi(0),
  vertx(0),
  verty(0),
  vertz(0),
  pdg(0),
  pdgmother(0),
  hasmother(0),
  motherlabel(0),
  enh(0),
  charge(0.),      
  IsRec(0),
  pass(0),
  IsPrim(0),
  Rej(kFALSE),
  fTree(0),
  fQAHist(0)
{

}

AliAnalysisTaskMLTreeMakerEff::AliAnalysisTaskMLTreeMakerEff(const char *name) :
  AliAnalysisTaskSE(name),
  fPIDResponse(0),      
  NCrossedRowsTPC(0),
  NClustersTPC(0),
  HasSPDfirstHit(0), 
  RatioCrossedRowsFindableClusters(0), 
  NTPCSignal(0),
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
  fESDTrackCuts(0),
  EsigTPC(0),
  EsigTOF(0),
  EsigITS(0),
  PsigTPC(0),
  PsigTOF(0),
  PsigITS(0),
  KsigTPC(0),
  KsigTOF(0),
  KsigITS(0),
  nITS(0),
  nITSshared(0),
  chi2ITS(0),
  chi2TPC(0),
  chi2Global(0),
  chi2GlobalvsTPC(0),
  fCutMaxChi2TPCConstrainedVsGlobalVertexType(0),
  ProdVx(0),
  ProdVy(0),
  ProdVz(0),
  dcar(),
  dcaz(),
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
  mcTrackIndex(0),
  fMcArray(0x0),
  MCpt(0),
  MCeta(0),
  MCphi(0),
  pt(0),      
  eta(0),
  phi(0),
  vertx(0),
  verty(0),
  vertz(0),
  pdg(0),
  pdgmother(0),
  hasmother(0),
  motherlabel(0),
  enh(0),
  charge(0.),      
  IsRec(0),
  pass(0),
  IsPrim(0),
  Rej(kFALSE),
  fTree(0),
  fQAHist(0)
{  

  DefineOutput(1, TList::Class());
  fESDTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kFALSE,0);
}


void AliAnalysisTaskMLTreeMakerEff::UserCreateOutputObjects() {
    
   AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
   AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
   inputHandler->SetNeedField();
   
  fPIDResponse = inputHandler->GetPIDResponse();
  if (!fPIDResponse)  return;

  fList = new TList();
  fList->SetName("output_Tlist");
  fList->SetOwner();
   
  
  if (man->GetMCtruthEventHandler()==0x0){
      AliError("No MC Truth available");
      return;
  }



  Bool_t oldStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);

  fTree = new TTree("Track_Tree_Gen","Tracks_Gen");
  fList->Add(fTree);
  
  fQAHist = new TH1F("h1", "h1 title", 4, 0, 1);
  fList->Add(fQAHist);
  
  fTree->Branch("centrality", &cent);
  fTree->Branch("#tracks", &n);
  fTree->Branch("charge", &charge);
  fTree->Branch("RunNumber", &runn);
  fTree->Branch("EsigTPC", &EsigTPC);
  fTree->Branch("EsigITS", &EsigITS);
  fTree->Branch("EsigTOF", &EsigTOF);
  
  fTree->Branch("NCrossedRowsTPC", &NCrossedRowsTPC);
  fTree->Branch("NClustersTPC", &NClustersTPC);
  fTree->Branch("RatioCrossedRowsFindableClusters", &RatioCrossedRowsFindableClusters);
  fTree->Branch("HasSPDfirstHit", &HasSPDfirstHit);   
  fTree->Branch("NTPCSignal", &NTPCSignal); 
  
  
   if(fPionSigmas){
    fTree->Branch("PsigTPC", &PsigTPC);
    fTree->Branch("PsigITS", &PsigITS);
    fTree->Branch("PsigTOF", &PsigTOF);
  }
  if(fKaonSigmas){
    fTree->Branch("KsigTPC", &KsigTPC);
    fTree->Branch("KsigITS", &KsigITS);
    fTree->Branch("KsigTOF", &KsigTOF);
  }
  
  fTree->Branch("nITS", &nITS);
  fTree->Branch("nITSshared_frac", &nITSshared);
  fTree->Branch("chi2ITS", &chi2ITS);
  fTree->Branch("chi2TPC", &chi2TPC);
  fTree->Branch("chi2GlobalvsTPC", &chi2GlobalvsTPC);
  
  fTree->Branch("DCAxy", &dcar);
  fTree->Branch("DCAz", &dcaz);
  
  fTree->Branch("vertx", &vertx);
  fTree->Branch("verty", &verty);
  fTree->Branch("vertz", &vertz);
  
  fTree->Branch("ProdVx", &ProdVx);
  fTree->Branch("ProdVy", &ProdVy);
  fTree->Branch("ProdVz", &ProdVz);
      
  fTree->Branch("Pdg", &pdg);
  fTree->Branch("Pdg_Mother", &pdgmother);
  fTree->Branch("Mother_label", &motherlabel);
  fTree->Branch("Has_Mother", &hasmother);
  fTree->Branch("IsEnh", &enh);
  fTree->Branch("IsRec", &IsRec);
  fTree->Branch("PassCuts", &pass);
  fTree->Branch("IsPrim", &IsPrim);
  
  fTree->Branch("MCpt", &MCpt);
  fTree->Branch("MCeta", &MCeta);
  fTree->Branch("MCphi", &MCphi);
  
  fTree->Branch("pt", &pt);
  fTree->Branch("eta", &eta);
  fTree->Branch("phi", &phi);
  
  
  PostData(1, fList);
  
  AliInfo("Finished setting up the Output");
  TH1::AddDirectory(oldStatus);
  
}

//________________________________________________________________________

void AliAnalysisTaskMLTreeMakerEff::UserExec(Option_t *) {
  // Main loop


  // Called for each event
  AliESDEvent* esdevent = dynamic_cast<AliESDEvent*>(InputEvent());
  AliVEvent* event = dynamic_cast<AliVEvent*>(InputEvent()); 
  
  fQAHist->Fill("Events_all",1);
  
  if(!event) {
    AliError("event not available");
    return;
  }

  // check event cuts
  if( IsEventAccepted(esdevent) == 0){ 
    return;
  }

  fQAHist->Fill("Events_accepted",1);

  AliMCEventHandler* mchandler = dynamic_cast<AliMCEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
  fMcArray = mchandler->MCEvent();

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
}


void  AliAnalysisTaskMLTreeMakerEff::FinishTaskOutput(){


}


void AliAnalysisTaskMLTreeMakerEff::Terminate(Option_t *) {


    
}

Double_t AliAnalysisTaskMLTreeMakerEff::IsEventAccepted(AliESDEvent *event){

  
    if (TMath::Abs(event->GetVertex()->GetZ()) < 10  &&  event->GetPrimaryVertexSPD() ){
      if (event->GetPrimaryVertexSPD()->GetNContributors() >0) return 1;
      else return 0;}
    return 0;
}




Int_t AliAnalysisTaskMLTreeMakerEff::GetAcceptedTracks(AliVEvent *event, Double_t gCentrality){

    
  ev1++;
  Int_t acceptedTracks = 0;
  Int_t temprec=0;
  Int_t passtrcuts;
  
  std::vector<std::vector<int>> vec;
  std::vector<int> row(3,0);        // 1 is el, 2 is rec, 3 pos (=label) in esd event of rec
  std::vector<int> ellabels;
  
  MCpt.clear();
  MCeta.clear();
  MCphi.clear(); 
  
  pt.clear();
  eta.clear();
  phi.clear(); 
  
  pdg.clear();
  pdgmother.clear();
  hasmother.clear();
  motherlabel.clear();
  
  charge.clear();
  enh.clear();
  IsRec.clear();
  pass.clear();
  IsPrim.clear();
  
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
  dcar.clear();
  dcaz.clear();
  nITS.clear();
  nITSshared.clear();
  chi2ITS.clear();
  chi2TPC.clear();
  chi2Global.clear();
  chi2GlobalvsTPC.clear();
  
  ProdVx.clear();
  ProdVy.clear();
  ProdVz.clear();
  
  
  // Loop over tracks in event
  AliGenCocktailEventHeader* coHeader;
  AliMCEvent *mcEvent;
  AliESDtrack* esdTrack;
  Int_t temppdg;
  Int_t tempmpdg;
  AliMCParticle* mcMTrack;
  AliMCParticle* mcTrack;
  
  mcEvent = MCEvent(); 
        if (!mcEvent) {
          AliError(Form("Could not receive MC Event!"));
          return -1;
  }
  
  //1st loop to get info for each MC particle if it is an electron and if primary in the 2d vec and to store where they are in the ellabels 
  for (Int_t iTracks = 0; iTracks < mcEvent->GetNumberOfTracks(); iTracks++) {
      mcTrack = dynamic_cast<AliMCParticle *>(mcEvent->GetTrack(iTracks));
      if (!mcTrack) {
	      AliWarning(Form("Could not receive MC particle %d", iTracks));
	      continue;}
      
     row[0]=0; 

     if( abs(mcTrack->PdgCode())==11) {row[0]=1;
              ellabels.push_back(iTracks);
     }
     
     vec.push_back(row); 
  
  }
  //2nd loop over ESD tracks to store for each ESD track that there is a rec track in the 2d Vec 
  for (Int_t iTracks = 0; iTracks < event->GetNumberOfTracks(); iTracks++) {
      AliESDtrack* esdTrack1 = dynamic_cast<AliESDtrack *>(event->GetTrack(iTracks));
      if (!esdTrack1) {
	      AliError(Form("Could not receive ESD track %d", iTracks));
	      continue;
      } 
      if(vec[TMath::Abs(esdTrack1->GetLabel())][0] == 0) continue;      //if not electron go on
      vec[TMath::Abs(esdTrack1->GetLabel())][1] = 1;
      vec[TMath::Abs(esdTrack1->GetLabel())][2] = iTracks;
  }
  
  
// loop over electrons
  for (Int_t i = 0; i < ellabels.size(); i++) {  
      Int_t iTracks = ellabels[i];
      mcTrack = dynamic_cast<AliMCParticle *>(mcEvent->GetTrack(iTracks));
      if (!mcTrack) {
	      AliWarning(Form("Could not receive MC particle %d", iTracks));
	      continue;
      }
        
      if( abs(mcTrack->PdgCode())!=11 ) AliError("This is not an electron!");
      
      
      Rej=kFALSE;
 
      mcMTrack = dynamic_cast<AliMCParticle *>(mcEvent->GetTrack(TMath::Abs(mcTrack->GetLabel())));
      
      // Kinematic MC cuts
      Double_t pttemp = mcTrack->Pt();
      Double_t etatemp = mcTrack->Eta();
      
      fQAHist->Fill("Before MC track cuts",1); 
      
      if( pttemp > fPtMax || pttemp < fPtMin ) continue;
      if( etatemp > fEtaMax || etatemp < fEtaMin ) continue;
      if( abs(mcTrack->PdgCode())!=11)  continue;
      
      fQAHist->Fill("After MC track cuts",1); 

      
      //check if particle is from enh signal
      temppdg = mcMTrack->PdgCode();  
      if(!(mcMTrack->GetMother() < 0)){       //get direct mother
          mcTrackIndex = mcMTrack->GetMother(); 
          mcMTrack = dynamic_cast<AliMCParticle *>(mcEvent->GetTrack(mcMTrack->GetMother()));
          tempmpdg= mcMTrack->PdgCode(); 
      }
      else tempmpdg=-9999;
 
      while(!(mcMTrack->GetMother() < 0)){        //get first mother in chain
          mcTrackIndex = mcMTrack->GetMother(); 
          mcMTrack = dynamic_cast<AliMCParticle *>(mcEvent->GetTrack(mcMTrack->GetMother()));  
      }
 
//      if(!(mcEvent->IsFromBGEvent(abs(mcTrackIndex)))) Rej=kTRUE;
      if(!(IsFromBGEventAOD(mcEvent,abs(mcTrackIndex)))) Rej=kTRUE;
      passtrcuts = 0;          
          
      //check if MC particle was reconstructed      
      temprec=vec[iTracks][1];
      
      if(temprec) {
          
      esdTrack = dynamic_cast<AliESDtrack *>(event->GetTrack(vec[iTracks][2]));
      if (!esdTrack) {
	      AliError(Form("Could not receive ESD track %d", vec[iTracks][2]));
      }   
      
      passtrcuts=1;
      
      if(!fESDTrackCuts->AcceptTrack(esdTrack))   passtrcuts = 0;
      
      fQAHist->Fill("Is Rec bef ESD cut",1); 
      
      // Kinematic cuts
      Double_t pttemprec = esdTrack->Pt();
      Double_t etatemprec = esdTrack->Eta();
      
      if( pttemprec > fPtMax || pttemprec < fPtMin ) passtrcuts = 0;
      if( etatemprec > fEtaMax || etatemprec < fEtaMin ) passtrcuts = 0;
 
      Double_t tempEsigTPC=fPIDResponse->NumberOfSigmasTPC(esdTrack, (AliPID::EParticleType) 0);
      Double_t tempEsigITS=fPIDResponse->NumberOfSigmasITS(esdTrack, (AliPID::EParticleType) 0);
      Double_t tempEsigTOF=fPIDResponse->NumberOfSigmasTOF(esdTrack, (AliPID::EParticleType) 0);
      
      if(fUsePionPIDTPC){
        if (fPIDResponse->NumberOfSigmasTPC(esdTrack, (AliPID::EParticleType) 2) > fPSigTPCMin &&  fPIDResponse->NumberOfSigmasTPC(esdTrack, (AliPID::EParticleType) 2)  < fPSigTPCMax){ passtrcuts = 0;} //exclude pions in TPC
      }

      if (fPIDResponse->CheckPIDStatus(AliPIDResponse::kTOF,esdTrack)==AliPIDResponse::kDetPidOk && (tempEsigTOF < fESigTOFMin || tempEsigTOF > fESigTOFMax)) passtrcuts = 0; 
      if (tempEsigITS < fESigITSMin || tempEsigITS > fESigITSMax) passtrcuts = 0; 
      if (tempEsigTPC < fESigTPCMin || tempEsigTPC > fESigTPCMax) passtrcuts = 0;
      
      
     fQAHist->Fill("Is Rec after ESD cut",1);
      
      
      //Fill Tree with rec data
      EsigTPC.push_back(tempEsigTPC);
      EsigITS.push_back(tempEsigITS);
      EsigTOF.push_back(tempEsigTOF);
      if(fPionSigmas){
        Double_t tempPsigTPC=fPIDResponse->NumberOfSigmasTPC(esdTrack, (AliPID::EParticleType) 2);
        Double_t tempPsigITS=fPIDResponse->NumberOfSigmasITS(esdTrack, (AliPID::EParticleType) 2);
        Double_t tempPsigTOF=fPIDResponse->NumberOfSigmasTOF(esdTrack, (AliPID::EParticleType) 2);
        PsigTPC.push_back(tempPsigTPC);
        PsigITS.push_back(tempPsigITS);
        PsigTOF.push_back(tempPsigTOF);
      }
      if(fKaonSigmas){
        Double_t tempKsigTPC=fPIDResponse->NumberOfSigmasTPC(esdTrack, (AliPID::EParticleType) 3);
        Double_t tempKsigITS=fPIDResponse->NumberOfSigmasITS(esdTrack, (AliPID::EParticleType) 3);
        Double_t tempKsigTOF=fPIDResponse->NumberOfSigmasTOF(esdTrack, (AliPID::EParticleType) 3);
        KsigTPC.push_back(tempKsigTPC);
        KsigITS.push_back(tempKsigITS);
        KsigTOF.push_back(tempKsigTOF);
      }
      eta.push_back(etatemprec);
      phi.push_back(esdTrack->Phi());
      pt.push_back(pttemprec);
      charge.push_back(esdTrack->Charge());   
  
      NCrossedRowsTPC.push_back(esdTrack->GetTPCCrossedRows());
      NClustersTPC.push_back(esdTrack->GetNumberOfTPCClusters());
      HasSPDfirstHit.push_back(esdTrack->HasPointOnITSLayer(0)); 
      RatioCrossedRowsFindableClusters.push_back((Double_t) esdTrack->GetTPCCrossedRows()/ (Double_t) esdTrack->GetTPCNclsF());       
      NTPCSignal.push_back(esdTrack->GetTPCsignalN());

      Float_t tempdca[2] = {0};
      esdTrack->GetImpactParameters( &tempdca[0], &tempdca[1]); //GetImpactParameter is also used in AliESDtrackCuts.cxx to cut on DCA to verte

      dcar.push_back(tempdca[0]);
      dcaz.push_back(tempdca[1]);

      Int_t tempnits = esdTrack->GetNcls(0);    // 0 = ITS 
      nITS.push_back(tempnits);        
      Double_t nitssharedtemp = 0.;
 
      if(tempnits){
        for(int d = 0; d<6;d++){
          nitssharedtemp+= (Double_t) esdTrack->HasSharedPointOnITSLayer(d);
        }
        nitssharedtemp/=tempnits;
      }

      nITSshared.push_back(nitssharedtemp);
      
      chi2ITS.push_back(esdTrack->GetITSchi2());
      chi2TPC.push_back(esdTrack->GetTPCchi2());
      
      fCutMaxChi2TPCConstrainedVsGlobalVertexType = fESDTrackCuts->kVertexTracks | fESDTrackCuts->kVertexSPD;

      const AliESDVertex* vertex = 0;
      if (fCutMaxChi2TPCConstrainedVsGlobalVertexType & fESDTrackCuts->kVertexTracks){
        vertex = esdTrack->GetESDEvent()->GetPrimaryVertexTracks();}
      
      if ((!vertex || !vertex->GetStatus()) && fCutMaxChi2TPCConstrainedVsGlobalVertexType & fESDTrackCuts->kVertexSPD){
	      vertex = esdTrack->GetESDEvent()->GetPrimaryVertexSPD();}
	
      if ((!vertex || !vertex->GetStatus()) && fCutMaxChi2TPCConstrainedVsGlobalVertexType & fESDTrackCuts->kVertexTPC){
	      vertex = esdTrack->GetESDEvent()->GetPrimaryVertexTPC();}

      if (vertex->GetStatus()){
        chi2GlobalvsTPC.push_back(esdTrack->GetChi2TPCConstrainedVsGlobal(vertex));}
      
      
         
      }
      
      
      else{   
      //Fill Tree with non MC dummy data for those not rec
      EsigTPC.push_back(-9999);
      EsigITS.push_back(-9999);
      EsigTOF.push_back(-9999);
      if(fPionSigmas){
        PsigTPC.push_back(-9999);
        PsigITS.push_back(-9999);
        PsigTOF.push_back(-9999);
      }
      if(fKaonSigmas){
        KsigTPC.push_back(-9999);
        KsigITS.push_back(-9999);
        KsigTOF.push_back(-9999);
      }
      eta.push_back(-9999);
      phi.push_back(-9999);
      pt.push_back(-9999);
      charge.push_back(-9999);   
      NCrossedRowsTPC.push_back(-9999);
      NClustersTPC.push_back(-9999);
      HasSPDfirstHit.push_back(-9999); 
      RatioCrossedRowsFindableClusters.push_back(-9999);       
      NTPCSignal.push_back(-9999);
      dcar.push_back(-9999);
      dcaz.push_back(-9999);
      nITS.push_back(-9999);        
      nITSshared.push_back(-9999);
      chi2ITS.push_back(-9999);
      chi2TPC.push_back(-9999);
      chi2GlobalvsTPC.push_back(-9999);
      }
      
      
      //Fill MC Data
      fQAHist->Fill("Tracks aft MC&Hij, bef tr cuts",1); 
      
      ProdVx.push_back(mcTrack->Xv());
      ProdVy.push_back(mcTrack->Yv());
      ProdVz.push_back(mcTrack->Zv());

      pdg.push_back( mcTrack->PdgCode());
      
      if(Rej) enh.push_back(1);
      else enh.push_back(0);
      
      pass.push_back(passtrcuts);
      IsRec.push_back(temprec);
      
      if(mcEvent->IsPhysicalPrimary(iTracks)) IsPrim.push_back(1);
      else IsPrim.push_back(0);
      
      MCpt.push_back(pttemp);
      MCeta.push_back(etatemp);
      MCphi.push_back(mcTrack->Phi());
      
      charge.push_back(mcTrack->Charge());  

      //Get vertex only for first track in event
      if(!acceptedTracks){     
        Double_t vert[3] = {0};
 
        mcTrack->XvYvZv(vert);
        
        vertx= vert[0];
        verty= vert[1];
        vertz= vert[2];
      }
      
      if(!(mcTrack->GetMother() < 0)) {  
        hasmother.push_back(1);
        AliMCParticle* mcmother = dynamic_cast<AliMCParticle *>(fMCEvent->GetTrack(mcTrack->GetMother()));
	         pdgmother.push_back( mcmother->PdgCode());
 
        motherlabel.push_back(abs(mcmother->GetLabel()));
      }
      else{
        hasmother.push_back(0);  
        pdgmother.push_back( -9999);
        motherlabel.push_back(-9999);
      }

    fQAHist->Fill("Is written",1);
      // count tracks
      acceptedTracks++;
  }
    
  num1++;
  return acceptedTracks;  
}


Bool_t AliAnalysisTaskMLTreeMakerEff::IsFromBGEventAOD(AliMCEvent* fAOD, Int_t Index)
{
    //Check if the particle is from Hijing or Enhanced event
    AliAODMCHeader *mcHeader;
    Int_t fNBG =-1;
    
    mcHeader = dynamic_cast<AliAODMCHeader*>(fAOD->GetList()->FindObject(AliAODMCHeader::StdBranchName()));
    if (!mcHeader) {
        AliError("Could not find MC Header in AOD");
        return (0);
    }
    
    TList *List = mcHeader->GetCocktailHeaders();
    AliGenHijingEventHeader* hijingH = dynamic_cast<AliGenHijingEventHeader*>(List->FindObject("Hijing"));
    if (!hijingH){
        AliError("no GenHijing header");
        return (0);
    }
    fNBG = hijingH->NProduced();
    std::cout<<"hijingH->NProduced() = "<<fNBG<<std::endl;
    return (Index < fNBG);
}
