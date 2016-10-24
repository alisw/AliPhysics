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
#include "AliAnalysisTaskMLTreeMaker.h"
#include "AliTrackReference.h"
#include "AliHeader.h"
#include "AliGenEventHeader.h"
#include "AliGenHijingEventHeader.h"
#include "AliGenCocktailEventHeader.h"
#include "AliCentrality.h"

// Authors: Sebastian Lehner (SMI Vienna) - selehner@cern.ch


ClassImp(AliAnalysisTaskMLTreeMaker)

Bool_t cutonTPCsignalN=  kFALSE;       
Int_t num= 0;
Int_t ev=0;

AliAnalysisTaskMLTreeMaker::AliAnalysisTaskMLTreeMaker():
  AliAnalysisTaskSE(),
  fList(0x0),
  fCentralityPercentileMin(0),
  fCentralityPercentileMax(100), 
  fPtMin(0),
  fPtMax(1000),
  fEtaMin(-10),
  fEtaMax(10),
  fFilterBit(96),
  fMcArray(0x0),
  fTree(0),
  fQAHist(0),  
  eta(0),
  phi(0),
  pt(0),
  EsigTPC(0),
  EsigTOF(0),
  EsigITS(0),
  PsigTPC(0),
  PsigTOF(0),
  PsigITS(0),
  KsigTPC(0),
  KsigTOF(0),
  KsigITS(0),
  fESigITSMin(-100.),
  fESigITSMax(1.),
  fESigTPCMin(-3.),
  fESigTPCMax(3.),
  fESigTOFMin(-3),
  fESigTOFMax(3),
  fPSigTPCMin(-100.),
  fPSigTPCMax(4.),
  fPionSigmas(kFALSE),
  fKaonSigmas(kFALSE),
  fUsePionPIDTPC(kFALSE),
  hasMC(kFALSE),
  MCpt(0),
  MCeta(0),
  MCphi(0),
  pdg(0),
  pdgmother(0),
  hasmother(0),      
  dcar(),
  dcaz(),
  nITS(0),
  fESDTrackCuts(0),
  gMultiplicity(-999),
  chi2ITS(0),
  chi2TPC(0),
  chi2Global(0),
  nITSshared(0),
  chi2GlobalvsTPC(0),
  fCutMaxChi2TPCConstrainedVsGlobalVertexType(0),
  motherlabel(0),
  charge(0.),      
  runn(0),      
  Rej(kFALSE),
  fPIDResponse(0),
  n(0),
  cent(0),
  vertx(0),
  verty(0),
  vertz(0),
  mcTrackIndex(0),
  NCrossedRowsTPC(0),
  NClustersTPC(0),
  HasSPDfirstHit(0), 
  RatioCrossedRowsFindableClusters(0), 
  NTPCSignal(0),
  loCuts(kTRUE)      
{

}

AliAnalysisTaskMLTreeMaker::AliAnalysisTaskMLTreeMaker(const char *name) :
  AliAnalysisTaskSE(name),
  fList(0x0),
  fCentralityPercentileMin(0),
  fCentralityPercentileMax(100), 
  fPtMin(0),
  fPtMax(1000),
  fEtaMin(-10),
  fEtaMax(10),
  fFilterBit(96),
  fMcArray(0x0),
  fTree(0),
  fQAHist(0),  
  eta(0),
  phi(0),
  pt(0),
  EsigTPC(0),
  EsigTOF(0),
  EsigITS(0),
  PsigTPC(0),
  PsigTOF(0),
  PsigITS(0),
  KsigTPC(0),
  KsigTOF(0),
  KsigITS(0),
  fESigITSMin(-100.),
  fESigITSMax(1.),
  fESigTPCMin(-3.),
  fESigTPCMax(3.),
  fESigTOFMin(-3),
  fESigTOFMax(3),
  fPSigTPCMin(-100.),
  fPSigTPCMax(4.),
  fPionSigmas(kFALSE),
  fKaonSigmas(kFALSE),
  fUsePionPIDTPC(kFALSE),
  hasMC(kFALSE),
  MCpt(0),
  MCeta(0),
  MCphi(0),
  pdg(0),
  pdgmother(0),
  hasmother(0),      
  dcar(),
  dcaz(),
  nITS(0),
  fESDTrackCuts(0),
  gMultiplicity(-999),
  chi2ITS(0),
  chi2TPC(0),
  chi2Global(0),
  nITSshared(0),
  chi2GlobalvsTPC(0),
  fCutMaxChi2TPCConstrainedVsGlobalVertexType(0),
  motherlabel(0),
  charge(0.),      
  runn(0),      
  Rej(kFALSE),
  fPIDResponse(0),
  n(0),
  cent(0),
  vertx(0),
  verty(0),
  vertz(0),
  mcTrackIndex(0),
  NCrossedRowsTPC(0),
  NClustersTPC(0),
  HasSPDfirstHit(0), 
  RatioCrossedRowsFindableClusters(0), 
  NTPCSignal(0),
  loCuts(kTRUE)          
{

  if(loCuts){
  AliInfo(Form("Loose cuts!!"));
  fESDTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kFALSE,0);
  cutonTPCsignalN = kFALSE; 
  }
    else{  
  //Alberto Style ESD track cuts - according to analysis note v.6   
  AliInfo(Form("Alberto cuts!!"));
  fESDTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kTRUE,1);
  fESDTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kFirst);
  fESDTrackCuts->SetMaxDCAToVertexXYPtDep("0.00515869+0.0101668/pt^1.34489");
  fESDTrackCuts->SetMaxDCAToVertexZ(0.1);
  fESDTrackCuts->SetMinNClustersITS(4);
  fESDTrackCuts->SetMinNCrossedRowsTPC(100);
  fESDTrackCuts->SetMinNClustersTPC(70);
  fESDTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.6);
  cutonTPCsignalN = kTRUE;
  }
  
  DefineOutput(1, TList::Class());
}

//________________________________________________________________________

//~ AliAnalysisTaskMLTreeMaker::~AliAnalysisTaskMLTreeMaker() {

  //~ // Destructor

  //~ // ... not implemented

//~ }


//________________________________________________________________________

void AliAnalysisTaskMLTreeMaker::UserCreateOutputObjects() {
    
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
  
//  fTree->Branch("NClustersITS", &NClustersITS);
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
  
  fTree->Branch("DCAxy", &dcar);
  fTree->Branch("DCAz", &dcaz);
  
  fTree->Branch("vertx", &vertx);
  fTree->Branch("verty", &verty);
  fTree->Branch("vertz", &vertz);
  
  fTree->Branch("nITS", &nITS);
  fTree->Branch("nITSshared_frac", &nITSshared);
  fTree->Branch("chi2ITS", &chi2ITS);
  fTree->Branch("chi2TPC", &chi2TPC);
  fTree->Branch("chi2GlobalvsTPC", &chi2GlobalvsTPC);
  
  if(hasMC) {
      
    fTree->Branch("Pdg", &pdg);
    fTree->Branch("Pdg_Mother", &pdgmother);
    fTree->Branch("Mother_label", &motherlabel);
    fTree->Branch("Has_Mother", &hasmother); 
  
    fTree->Branch("MCpt", &MCpt);
    fTree->Branch("MCeta", &MCeta);
    fTree->Branch("MCphi", &MCphi);
  }
  
  PostData(1, fList);

//  PostData(2, fQAHist);
//  PostData(1, fTree);
  
  AliInfo("Finished setting up the Output");
  TH1::AddDirectory(oldStatus);
}

//________________________________________________________________________

void AliAnalysisTaskMLTreeMaker::UserExec(Option_t *) {
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
  
  if(hasMC){
    AliMCEventHandler* mchandler = dynamic_cast<AliMCEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
    fMcArray = mchandler->MCEvent();
  // get the accepted tracks in main event

  }
  Double_t lMultiplicityVar = -1;
  Int_t acceptedTracks = GetAcceptedTracks(event,lMultiplicityVar);
  //fHistTrackStats->Fill(acceptedTracks,lMultiplicityVar);

  AliCentrality *centrality = esdevent->GetCentrality();
  if (!centrality) AliError(Form("Could not receive Centrality"));  
            
  cent = centrality->GetCentralityPercentile("V0M");
  n= acceptedTracks;
  if(acceptedTracks){
    fTree->Fill();
    fQAHist->Fill("Events_track_selected",1);
  }
}

//~ //________________________________________________________________________

void  AliAnalysisTaskMLTreeMaker::FinishTaskOutput(){
  // Finish task output

  // not implemented ...


}
//~ 

//~ //________________________________________________________________________

void AliAnalysisTaskMLTreeMaker::Terminate(Option_t *) {
  // Draw result to the screen

  // Called once at the end of the query

  // not implemented ...


}
//~ 


//________________________________________________________________________

Double_t AliAnalysisTaskMLTreeMaker::IsEventAccepted(AliESDEvent *event){

  
    if (TMath::Abs(event->GetVertex()->GetZ()) < 10  &&  event->GetPrimaryVertexSPD() ){
      if (event->GetPrimaryVertexSPD()->GetNContributors() >0) return 1;
      else return 0;}
    return 0;
}


//________________________________________________________________________

Int_t AliAnalysisTaskMLTreeMaker::GetAcceptedTracks(AliVEvent *event, Double_t gCentrality){
  ev++;
  Int_t acceptedTracks = 0;

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
  chi2TPC.clear();
  chi2Global.clear();
  chi2GlobalvsTPC.clear();
  pdg.clear();
  pdgmother.clear();
  hasmother.clear();
  motherlabel.clear();
  charge.clear();
  
  // Loop over tracks in event
  AliGenCocktailEventHeader* coHeader;
  AliMCEvent *mcEvent;
  
  for (Int_t iTracks = 0; iTracks < event->GetNumberOfTracks(); iTracks++) {
      AliESDtrack* esdTrack = dynamic_cast<AliESDtrack *>(event->GetTrack(iTracks));
      if (!esdTrack) {
	      AliError(Form("Could not receive ESD track %d", iTracks));
	      continue;
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
      
          if(!iTracks){                             //check for first track if Hijing was the event gen
            Rej=kFALSE;
            coHeader = dynamic_cast<AliGenCocktailEventHeader*> (mcEvent->GenEventHeader());
            if (!coHeader){
              AliError(Form("Could not receive coHeader -> Rej set to kFALSE!! no rejection of enhanced sources!!"));
              //continue;
            }
            else{
                TList* list = coHeader->GetHeaders();
//                if(list->FindObject("Hijing")) Rej = kTRUE;
                
                AliMCParticle* mcTrack = dynamic_cast<AliMCParticle *>(mcEvent->GetTrack(TMath::Abs(esdTrack->GetLabel())));
                
                do{
                    mcTrackIndex = mcTrack->GetMother();                    
                    mcTrack = dynamic_cast<AliMCParticle *>(mcEvent->GetTrack(mcTrack->GetMother()));

                }
                while(!(mcTrack->GetMother() < 0));
                
                if(!(mcEvent->IsFromBGEvent(abs(mcTrackIndex)))) Rej=kTRUE;
            }
          }
          //Reject non-Hijing BG tracks    
          if (Rej){
            AliError(Form("Reject non-Hijing BG tracks!!"));
            continue;
          }
        }
      }
      fQAHist->Fill("Tracks aft MC&Hij, bef tr cuts",1); 
      
      
      if(!fESDTrackCuts->AcceptTrack(esdTrack))   continue;

      
//      Alberto Cut on TPC signal N (number of TPC clusters used for dE/dx)
      if(cutonTPCsignalN && esdTrack->GetTPCsignalN()<50) continue; 
      
      
      // Kinematic cuts
      Double_t pttemp = esdTrack->Pt();
      Double_t etatemp = esdTrack->Eta();
      
      if( pttemp > fPtMax || pttemp < fPtMin ) continue;
      if( etatemp > fEtaMax || etatemp < fEtaMin ) continue;
 
      Double_t tempEsigTPC=fPIDResponse->NumberOfSigmasTPC(esdTrack, (AliPID::EParticleType) 0);
      Double_t tempEsigITS=fPIDResponse->NumberOfSigmasITS(esdTrack, (AliPID::EParticleType) 0);
      Double_t tempEsigTOF=fPIDResponse->NumberOfSigmasTOF(esdTrack, (AliPID::EParticleType) 0);
      
      if(fUsePionPIDTPC){
        if (fPIDResponse->NumberOfSigmasTPC(esdTrack, (AliPID::EParticleType) 2) > fPSigTPCMin &&  fPIDResponse->NumberOfSigmasTPC(esdTrack, (AliPID::EParticleType) 2)  < fPSigTPCMax){ continue;} //exclude pions in TPC
      }

      if (fPIDResponse->CheckPIDStatus(AliPIDResponse::kTOF,esdTrack)==AliPIDResponse::kDetPidOk && (tempEsigTOF < fESigTOFMin || tempEsigTOF > fESigTOFMax)) continue;  
      if (tempEsigITS < fESigITSMin || tempEsigITS > fESigITSMax) continue;  
      if (tempEsigTPC < fESigTPCMin || tempEsigTPC > fESigTPCMax) continue;
      
      fQAHist->Fill("Selected tracks",1);  
      
      //Fill Tree with MC data
      if(hasMC){ 
        AliMCParticle* mcTrack = dynamic_cast<AliMCParticle *>(mcEvent->GetTrack(TMath::Abs(esdTrack->GetLabel())));

        pdg.push_back( mcTrack->PdgCode());

        MCpt.push_back(mcTrack->Pt());
        MCeta.push_back(mcTrack->Eta());
        MCphi.push_back(mcTrack->Phi());
        
        runn.push_back(event->GetRunNumber());
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
      } //End if hasMC 
      
      //Fill Tree with non MC data
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
      eta.push_back(etatemp);
      phi.push_back(esdTrack->Phi());
      pt.push_back(pttemp);
      charge.push_back(esdTrack->Charge());   

      NCrossedRowsTPC.push_back(esdTrack->GetTPCCrossedRows());
      NClustersTPC.push_back(esdTrack->GetNumberOfTPCClusters());
      HasSPDfirstHit.push_back(esdTrack->HasPointOnITSLayer(0)); 
      RatioCrossedRowsFindableClusters.push_back((Double_t) esdTrack->GetTPCCrossedRows()/ (Double_t) esdTrack->GetTPCNclsF());       
      NTPCSignal.push_back(esdTrack->GetTPCsignalN());
      
       //Get DCA position
      //cout<<num<<"  track:  "<<iTracks<<"  "<<" pt "<<pttemp;
      Float_t tempdca[2] = {0};
       esdTrack->GetImpactParameters( &tempdca[0], &tempdca[1]); //GetImpactParameter is also used in AliESDtrackCuts.cxx to cut on DCA to verte

//       cout<<"  dcaxy: "<<tempdca[0]<<endl;  

      dcar.push_back(tempdca[0]);
      dcaz.push_back(tempdca[1]);

      Int_t tempnits = esdTrack->GetNcls(0);    // 0 = ITS 
      nITS.push_back(tempnits);        
      Double_t nitssharedtemp = 0.;
 
      if(tempnits){
        for(int d = 0; d<6;d++){
          nitssharedtemp+= (Double_t) esdTrack->HasSharedPointOnITSLayer(d);
        }
//              if(nitssharedtemp) cout<<"frac: "<<nitssharedtemp<<endl;
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
 
      // count tracks
      acceptedTracks++;
  }
    
  num++;
  return acceptedTracks;  
}

