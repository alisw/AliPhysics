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

using std::vector;

// Authors: Sebastian Lehner (SMI Vienna) - selehner@cern.ch


ClassImp(AliAnalysisTaskMLTreeMaker)

Int_t num= 0;
Int_t ev=0;


AliAnalysisTaskMLTreeMaker::AliAnalysisTaskMLTreeMaker():
AliAnalysisTaskSE(),
  fList(0x0),
  fCentralityPercentileMin(0.),
  fCentralityPercentileMax(100.), 
  fPtMin(0),
  fPtMax(1000),
  fEtaMin(-10),
  fEtaMax(10),
  fFilterBit(96),
  fMcArray(0x0),
  fTree(0),
  fQAHist(0),  
  eta(),
  phi(),
  pt(),
  EsigTPC(),
  EsigTOF(),
  EsigITS(),
  hasMC(kFALSE),
  MCpt(),
  MCeta(),
  MCphi(),
  pdg(),
  pdgmother(0),
  dcar(),
  dcaz(),
  nITS(),
  fESDTrackCuts(0),
  gMultiplicity(-999),
  chi2ITS(),
  chi2TPC(),
  chi2Global(),
  nITSshared(),
  chi2GlobalvsTPC(),
  fCutMaxChi2TPCConstrainedVsGlobalVertexType(),
  tempdca(),
  motherlabel(-999.),
  charge(),      
  IsBG(),
  runn()      

{

}

AliAnalysisTaskMLTreeMaker::AliAnalysisTaskMLTreeMaker(const char *name) :
  AliAnalysisTaskSE(name),
  fList(0x0),
  fCentralityPercentileMin(0.),
  fCentralityPercentileMax(100.), 
  fPtMin(0),
  fPtMax(1000),
  fEtaMin(-10),
  fEtaMax(10),
  fFilterBit(96),
  fMcArray(0x0),
  fTree(0),
  fQAHist(0),      
  eta(),
  phi(),
  pt(),
  EsigTPC(),
  EsigTOF(),
  EsigITS(),
  hasMC(kFALSE),
  MCpt(),
  MCeta(),
  MCphi(),
  pdg(0),
  pdgmother(0),
  dcar(),
  dcaz(),
  nITS(),
  fESDTrackCuts(0),
  gMultiplicity(-999),
  chi2ITS(),
  chi2TPC(),
  chi2Global(),
  nITSshared(),
  chi2GlobalvsTPC(),
  fCutMaxChi2TPCConstrainedVsGlobalVertexType(),
  tempdca(),
  motherlabel(-999.),
  charge(),      
  IsBG(),
  runn()      
{
  //DefineInput(0, TChainvPt::Class());

  //~ // Output slot #0 writes into a TH1 container

  fESDTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kFALSE,1);
  DefineOutput(1, TTree::Class());
  DefineOutput(2, TH1::Class());
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
  
  fTree->Branch("#tracks", &n);
  fTree->Branch("pt", &pt);
  fTree->Branch("eta", &eta);
  fTree->Branch("phi", &phi);
  fTree->Branch("charge", &charge);
  fTree->Branch("ISBG", &IsBG);   
  fTree->Branch("RunNumber", &runn);
  fTree->Branch("EsigTPC", &EsigTPC);
  fTree->Branch("EsigITS", &EsigITS);
  fTree->Branch("EsigTOF", &EsigTOF);
  
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
      AliMCEvent* mcEvent = mchandler->MCEvent();

    fMcArray = mcEvent;
  // get the accepted tracks in main event

  }
  Double_t lMultiplicityVar = -1;
            Int_t acceptedTracks = GetAcceptedTracks(event,lMultiplicityVar);
            //fHistTrackStats->Fill(acceptedTracks,lMultiplicityVar);

//            n=event->GetNumberOfTracks();

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
EsigTPC.clear();
EsigTOF.clear();
EsigITS.clear();
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
IsBG.clear();
    // Loop over tracks in event


    for (Int_t iTracks = 0; iTracks < event->GetNumberOfTracks(); iTracks++) {
//        if(ev!=10) continue;

      AliESDtrack* esdTrack = dynamic_cast<AliESDtrack *>(event->GetTrack(iTracks));
//   AliESDtrack* fixesdTrack = dynamic_cast<AliESDtrack *>(event->GetTrack(14));

     
      if (!esdTrack) {
	AliError(Form("Could not receive ESD track %d", iTracks));
	continue;
      }

    fQAHist->Fill("Tracks_selected_before_Track_Cuts",1);   
    if(!fESDTrackCuts->AcceptTrack(esdTrack))   continue;
      

            // kinematic cuts

      pttemp = esdTrack->Pt();
      etatemp = esdTrack->Eta();
      
      if( pttemp > fPtMax || pttemp < fPtMin ) continue;
      if( etatemp > fEtaMax || etatemp < fEtaMin ) continue;
 
      
      tempEsigTPC=fPIDResponse->NumberOfSigmasTPC(esdTrack, (AliPID::EParticleType) 0);
      tempEsigITS=fPIDResponse->NumberOfSigmasITS(esdTrack, (AliPID::EParticleType) 0);
      tempEsigTOF=fPIDResponse->NumberOfSigmasTOF(esdTrack, (AliPID::EParticleType) 0);



      if (fPIDResponse->NumberOfSigmasTPC(esdTrack, (AliPID::EParticleType) 2) > -100. &&  fPIDResponse->NumberOfSigmasTPC(esdTrack, (AliPID::EParticleType) 2)  < 4.) continue; //exclude pions in TPC

      if (fPIDResponse->CheckPIDStatus(AliPIDResponse::kTOF,esdTrack)==AliPIDResponse::kDetPidOk && (tempEsigTOF<-3. || tempEsigTOF > 3.)) continue;  
      if (tempEsigITS<-4. || tempEsigITS > 1.) continue;  
      if (tempEsigTPC<-1.5 || tempEsigTPC > 3.) continue;
      
      
      if(hasMC){ 
      AliMCEvent *mcEvent = MCEvent(); 
      
      
      if (!mcEvent) {
        AliError(Form("Could not receive MC -> hasMC set to kFALSE!!"));
        hasMC=kFALSE;
        continue;
        }
      
      else{
          
      //reject non Hijing tracks    

//      if (!(mcEvent->IsFromBGEvent(esdTrack->GetLabel()))) IsBG.push_back(1);

//      else IsBG.push_back(0); 

          
      if (!(mcEvent->IsFromBGEvent(esdTrack->GetLabel()))) continue;
          
          
      //Fill Tree with MC data

          
      AliMCParticle* mcTrack = dynamic_cast<AliMCParticle *>(mcEvent->GetTrack(TMath::Abs(esdTrack->GetLabel())));

      pdg.push_back( mcTrack->PdgCode());

      MCpt.push_back(mcTrack->Pt());
      MCeta.push_back(mcTrack->Eta());
      MCphi.push_back(mcTrack->Phi());
      
      runn.push_back(event->GetRunNumber());
      
      if(!acceptedTracks){      //get Vertex only for first track in event

      mcTrack->XvYvZv(vert);
      
      vertx= vert[0];
      verty= vert[1];
      vertz= vert[2];

      }
      
      if(!(mcTrack->GetMother() < 0)) {  
        hasmother.push_back(1);
        AliMCParticle* mcmother = dynamic_cast<AliMCParticle *>(fMCEvent->GetTrack(mcTrack->GetMother()));
	pdgmother.push_back( mcmother->PdgCode());
        motherlabel.push_back(mcmother->GetLabel());
   
                }
      else{
              hasmother.push_back(0);  
              pdgmother.push_back( -9999);
              motherlabel.push_back(-9999);
            }

        }
      }
      
    //Fill Tree with non MC data

      
      EsigTPC.push_back(tempEsigTPC);
      EsigITS.push_back(tempEsigITS);
      EsigTOF.push_back(tempEsigTOF);
      

      eta.push_back(etatemp);
      phi.push_back(esdTrack->Phi());
      pt.push_back(pttemp);
      charge.push_back(esdTrack->Charge());         

       //Get DCA position

//       cout<<num<<"  track:  "<<iTracks<<"  "<<" pt "<<pttemp;

       esdTrack->GetImpactParameters(& tempdca[0],& tempdca[1]);       //GetImpactParameter is also used in AliESDtrackCuts.cxx to cut on DCA to verte

//       cout<<"  dcaxy: "<<tempdca[0]<<endl;  

       dcar.push_back(tempdca[0]);
       dcaz.push_back(tempdca[1]);
      
       
      tempnits=esdTrack->GetNcls(0);    // 0 = ITS 

      nITS.push_back(tempnits);        
      nitssharedtemp=0;
 
      if(tempnits){
        for(int d = 0; d<6;d++){
          nitssharedtemp+= (Int_t) esdTrack->HasSharedPointOnITSLayer(d);
         }
//              if(nitssharedtemp) cout<<"frac: "<<nitssharedtemp<<endl;

      nitssharedtemp/=tempnits;
      }

      nITSshared.push_back(nitssharedtemp);
      
      chi2ITS.push_back(esdTrack->GetITSchi2());
      chi2TPC.push_back(esdTrack->GetTPCchi2());
      
      fCutMaxChi2TPCConstrainedVsGlobalVertexType = fESDTrackCuts->kVertexTracks | fESDTrackCuts->kVertexSPD;

      const AliESDVertex* vertex = 0;
      if (fCutMaxChi2TPCConstrainedVsGlobalVertexType & fESDTrackCuts->kVertexTracks)
                 
	vertex = esdTrack->GetESDEvent()->GetPrimaryVertexTracks();
      
      if ((!vertex || !vertex->GetStatus()) && fCutMaxChi2TPCConstrainedVsGlobalVertexType & fESDTrackCuts->kVertexSPD)
	vertex = esdTrack->GetESDEvent()->GetPrimaryVertexSPD();
	
      if ((!vertex || !vertex->GetStatus()) && fCutMaxChi2TPCConstrainedVsGlobalVertexType & fESDTrackCuts->kVertexTPC)
	vertex = esdTrack->GetESDEvent()->GetPrimaryVertexTPC();
      if (vertex->GetStatus())

	chi2GlobalvsTPC.push_back(esdTrack->GetChi2TPCConstrainedVsGlobal(vertex));
                    

      
 
       
 
      // count tracks

      acceptedTracks++;
      fQAHist->Fill("Tracks_selected",1);  
    }

    
  num++;
  return acceptedTracks;  
}

