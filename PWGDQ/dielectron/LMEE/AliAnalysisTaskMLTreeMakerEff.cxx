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
#include "AliAnalysisTaskMLTreeMakerEff.h"
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
#include "TSystem.h"


// Authors: Sebastian Lehner (SMI Vienna) - selehner@cern.ch


ClassImp(AliAnalysisTaskMLTreeMakerEff)




AliAnalysisTaskMLTreeMakerEff::AliAnalysisTaskMLTreeMakerEff():
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
  TMVAReader(0),
  useTMVA(kFALSE),        
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
  ZDCepA(0),
  ZDCepC(0),
  TPCep(0),
  TPCepA(0),
  TPCepC(0),
  fQnList(0),        
  man(0),        
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
  nITSTMVA(0),
  ITS1SharedTMVA(0),
  ITS2SharedTMVA(0),
  ITS3SharedTMVA(0),
  ITS4SharedTMVA(0),
  ITS5SharedTMVA(0),
  ITS6SharedTMVA(0),
  nITSshared_fracTMVA(0),
  NCrossedRowsTPCTMVA(0),
  NClustersTPCTMVA(0),
  NTPCSignalTMVA(0),
  logDCAxyTMVA (0),
  logDCAzTMVA (0),   
  chi2GlobalPerNDFTMVA(0),
  chi2ITSTMVA(0),
  etaTMVA(0),
  phiTMVA(0),
  ptTMVA(0),   
  centTMVA(0),      
  MVAout(0),
  TMVAWeightFileName(0),      
  fwidthTPC(0), 
  fmeanTPC(0), 
  fwidthITS(0), 
  fmeanITS(0), 
  fwidthTOF(0), 
  fmeanTOF(0), 
  fIsTMVAInit(kFALSE),      
  fuseCorr(kFALSE),      
  fTree(0),
  fQAHist(0)
{

}

AliAnalysisTaskMLTreeMakerEff::AliAnalysisTaskMLTreeMakerEff(const char *name) :
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
  TMVAReader(0),
  useTMVA(kFALSE),        
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
  ZDCepA(0),        
  ZDCepC(0),
  TPCep(0),
  TPCepA(0),
  TPCepC(0),     
  fQnList(0),         
  man(0),         
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
  nITSTMVA(0),
  ITS1SharedTMVA(0),
  ITS2SharedTMVA(0),
  ITS3SharedTMVA(0),
  ITS4SharedTMVA(0),
  ITS5SharedTMVA(0),
  ITS6SharedTMVA(0),
  nITSshared_fracTMVA(0),
  NCrossedRowsTPCTMVA(0),
  NClustersTPCTMVA(0),
  NTPCSignalTMVA(0),
  logDCAxyTMVA (0),
  logDCAzTMVA (0),   
  chi2GlobalPerNDFTMVA(0),
  chi2ITSTMVA(0),
  etaTMVA(0),
  phiTMVA(0),
  ptTMVA(0),   
  centTMVA(0),      
  MVAout(0),
  TMVAWeightFileName(0),      
  fwidthTPC(0), 
  fmeanTPC(0), 
  fwidthITS(0), 
  fmeanITS(0), 
  fwidthTOF(0), 
  fmeanTOF(0), 
  fIsTMVAInit(kFALSE),           
  fuseCorr(kFALSE),      
  fTree(0),
  fQAHist(0)
{
// SetupTrackCuts(); 
// SetupEventCuts(); 
// AliInfo("Track & Event cuts were set"); 
//  if(useTMVA) SetupTMVAReader("TMVAClassification_BDTG.weights_094.xml");
  DefineOutput(1, TList::Class());
}

AliAnalysisTaskMLTreeMakerEff::~AliAnalysisTaskMLTreeMakerEff(){
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

void AliAnalysisTaskMLTreeMakerEff::UserCreateOutputObjects() {

    if (useTMVA) SetupTMVAReader(TMVAWeightFileName);  
  
    
   man=AliAnalysisManager::GetAnalysisManager();
   AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
   inputHandler->SetNeedField();

   
  fList = new TList();
  fList->SetName("output_Tlist");
  fList->SetOwner();
   
  std::cout <<"Running on MC!"<< std::endl;

  Bool_t oldStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);

  fTree = new TTree("Track_Tree","Tracks");
  fList->Add(fTree);
  
  fQAHist = new TH1F("h1", "h1 title", 4, 0, 1);
  fList->Add(fQAHist);
  
 
  fTree->Branch("#tracks", &n);      
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
    
  
  
  PostData(1, fList);
  
  AliInfo("Finished setting up the Output");
  TH1::AddDirectory(oldStatus);
  

  TString generatorName = "Hijing_0;pizero_1;eta_2;etaprime_3;rho_4;omega_5;phi_6;jpsi_7;Pythia CC_8;Pythia BB_8;Pythia B_8;Starlight_0;Hijing_1;";
  TObjArray arr = *(generatorName.Tokenize(";"));
  std::cout << "Used Generators: " << std::endl;
  for (int i = 0; i < arr.GetEntries(); ++i){
    TString temp = arr.At(i)->GetName();
    std::cout << "--- " << temp << std::endl;
    fGeneratorHashs.push_back(temp.Hash());
    }
}
//________________________________________________________________________

void AliAnalysisTaskMLTreeMakerEff::UserExec(Option_t *) {
  // Called for each event
  
  fQAHist->Fill("Events_all",1);
  
  AliVEvent* event = dynamic_cast<AliVEvent*>(InputEvent()); 
  
  if(!event) {
    AliError("event not available");
    fQAHist->Fill("Events_not_available",1);
    return;
  }


  
  AliMultSelection *MultSelection = 0x0; 
  MultSelection = (AliMultSelection * ) event->FindListObject("MultSelection");

  if( !MultSelection) {
   //If you get this warning (and lPercentiles 300) please check that the AliMultSelectionTask actually ran (before your task)
   AliWarning("AliMultSelection object not found!");
  }
  else cent = MultSelection->GetMultiplicityPercentile("V0M",kFALSE);

  
 AliInputEventHandler *eventHandler = nullptr;
 AliInputEventHandler *eventHandlerMC = nullptr;
  
  if ((AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler())->IsA() == AliAODInputHandler::Class()){
    eventHandler = dynamic_cast<AliAODInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
    eventHandlerMC = eventHandler;
  }
     

  fMcArray = eventHandlerMC->MCEvent();


  fQAHist->Fill("Events",1);

  
  Double_t lMultiplicityVar = -1;
  Int_t acceptedTracks = GetAcceptedTracks(event,lMultiplicityVar);
  
  runn = event->GetRunNumber();

  n= acceptedTracks;
     
 
//  AliAnalysisTaskFlowVectorCorrections *flowQnVectorTask = dynamic_cast<AliAnalysisTaskFlowVectorCorrections*> (man->GetTask("FlowQnVectorCorrections"));  
//  AliQnCorrectionsManager *flowQnVectorMgr = flowQnVectorTask->GetAliQnCorrectionsManager();
//  TList *qnlist = flowQnVectorMgr->GetQnVectorList();
//  if(qnlist != NULL)  AliAnalysisTaskMLTreeMakerEff::FillQnEventplanes(qnlist);
  
  fTree->Fill();
  fQAHist->Fill("Events_track_and_cent_selected",1);


  PostData(1, fList);
}

//~ //________________________________________________________________________

void  AliAnalysisTaskMLTreeMakerEff::FinishTaskOutput(){
  // Finish task output

  // not implemented ...


}
//~ 

//~ //________________________________________________________________________

void AliAnalysisTaskMLTreeMakerEff::Terminate(Option_t *) {

}
//~ 


//________________________________________________________________________

Double_t AliAnalysisTaskMLTreeMakerEff::IsEventAccepted(AliVEvent *event){

  if(event->GetPrimaryVertex()){
    if (TMath::Abs(event->GetPrimaryVertex()->GetZ()) < 10){
      if (event->GetPrimaryVertexSPD()->GetNContributors() > 0)
	return 1;
    }
  }
  return 0;
}


//________________________________________________________________________

Int_t AliAnalysisTaskMLTreeMakerEff::GetAcceptedTracks(AliVEvent *event, Double_t gCentrality){
  
  ev++;
  Int_t acceptedTracks = 0;
  Bool_t isAOD         = kFALSE;
  Int_t mpdg=0;

  MCpt.clear();
  MCeta.clear();
  MCphi.clear(); 
  
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
  
   
  
  
  // Loop over tracks in event
  AliMCEvent *mcEvent=0;
  AliAODMCParticle* mcMTrack;
  mcEvent = MCEvent(); 
  if (!mcEvent) {
    AliError(Form("Could not receive MC -> hasMC set to kFALSE!!"));
  }

  for (Int_t iTracks = 0; iTracks < mcEvent->GetNumberOfTracks(); iTracks++) {
    
  
        AliAODMCParticle* mcTrack = dynamic_cast<AliAODMCParticle *>(mcEvent->GetTrack(iTracks));
        
        Int_t genIndexTemp=CheckGenerator(iTracks);
        if(genIndexTemp==11)        iGenIndex.push_back(genIndexTemp); 
        else continue;
        if(!mcTrack->IsPhysicalPrimary()) continue;
        pdg.push_back( mcTrack->PdgCode());
        
        MCpt.push_back(mcTrack->Pt());
        if(mcTrack->Eta()>50 )        MCeta.push_back(50);
        else MCeta.push_back(mcTrack->Eta());
        MCphi.push_back(mcTrack->Phi());
        

        //Get production vertex for MC tracks

          Double_t MCvert[3] = {0};
          mcTrack->XvYvZv(MCvert);
          
          label.push_back(mcTrack->GetLabel());
          
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
      iPdgFirstMother.push_back(tempPdgFirstMother);
               
   


      acceptedTracks++;
  }

  return acceptedTracks;  
}

//______________________________________________________________________________
Bool_t AliAnalysisTaskMLTreeMakerEff::GetDCA(const AliVEvent* event, const AliAODTrack* track, Double_t* d0z0, Double_t* covd0z0)
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

void AliAnalysisTaskMLTreeMakerEff::FillQnEventplanes(TList* qnlist){
  TString qnListDetector;
  TString fgQnVectorNorm="";
  // TPC Eventplane q-Vector
  qnListDetector = "TPC" + fgQnVectorNorm;
  const AliQnCorrectionsQnVector *qVecQnFrameworkTPC = AliDielectronQnEPcorrection::GetQnVectorFromList(qnlist,qnListDetector.Data(),"latest","latest");
  TVector2 *qVectorTPC = new TVector2(-200.,-200.);
  if(qVecQnFrameworkTPC != NULL){
    qVectorTPC->Set(qVecQnFrameworkTPC->Qx(2),qVecQnFrameworkTPC->Qy(2));
    TPCep = TVector2::Phi_mpi_pi(qVectorTPC->Phi())/2;
  }
  delete qVectorTPC;
  
          
  // TPC A-Side/Neg. Eta Eventplane q-Vector
//  qnListDetector = "TPCNegEta" + fgQnVectorNorm;
//  const AliQnCorrectionsQnVector *qVecQnFrameworkTPCaSide = AliDielectronQnEPcorrection::GetQnVectorFromList(qnlist,qnListDetector.Data(),"latest","latest");
//  TVector2 *qVectorTPCaSide = new TVector2(-200.,-200.);
//  cout<<"TPC A-Side test"<<endl;
//  if(qVecQnFrameworkTPCaSide != NULL){
//    qVectorTPCaSide->Set(qVecQnFrameworkTPCaSide->Qx(2),qVecQnFrameworkTPCaSide->Qy(2));
//    TPCepA = TVector2::Phi_mpi_pi(qVectorTPCaSide->Phi())/2;
//    cout<<"qVecQnFrameworkTPCaSide->Qx(2):"<<qVecQnFrameworkTPCaSide->Qx(2)<<endl;
//  }
//  else cout<<"qVecQnFrameworkTPCaSide == NULL!"<<endl;
//  delete qVectorTPCaSide;

  // TPC C-Side/Pos. Eta Eventplane q-Vector
//  qnListDetector = "TPCPosEta" + fgQnVectorNorm;
//  const AliQnCorrectionsQnVector *qVecQnFrameworkTPCcSide = AliDielectronQnEPcorrection::GetQnVectorFromList(qnlist,qnListDetector.Data(),"latest","latest");
//  TVector2 *qVectorTPCcSide = new TVector2(-200.,-200.);
//  if(qVecQnFrameworkTPCcSide != NULL){
//    qVectorTPCcSide->Set(qVecQnFrameworkTPCcSide->Qx(2),qVecQnFrameworkTPCcSide->Qy(2));
//    TPCepC = TVector2::Phi_mpi_pi(qVectorTPCcSide->Phi())/2;
//  }
//  delete qVectorTPCcSide;        
}


Bool_t AliAnalysisTaskMLTreeMakerEff::FillZDCEventPlane(Double_t* ZDCevArr){
  // this is a copy of the AliDielectronVarManager

  AliFlowVector vQarray[2];
  
  AliAnalysisTaskZDCEP *fZDCEPTask = dynamic_cast<AliAnalysisTaskZDCEP*>(AliAnalysisManager::GetAnalysisManager()->GetTask("AnalysisTaskZDCEP"));

  if (fZDCEPTask != NULL) {

    // get ZDC Q-vectors
    TObjArray* dataContainers               = (AliAnalysisManager::GetAnalysisManager())->GetContainers();
    AliAnalysisDataContainer* dataContainer = dynamic_cast<AliAnalysisDataContainer*>(dataContainers->FindObject("ZDCEPExchangeContainer"));
    AliFlowEvent* anEvent                   = dynamic_cast<AliFlowEvent*>(dataContainer->GetData());
    if(anEvent) {
      // Get Q vectors for the subevents
      anEvent->GetZDC2Qsub(vQarray);
     } else { 
      Printf("AliAnalysisTaskMLTreeMakerEff::FillZDCEventPlane: Flowevent not found. Aborting!\n");
      return kFALSE;
    }
  } 
  else {
    Printf("This task needs AliAnalysisTaskZDCEP and it is not present. Aborting!!!");
    return kFALSE;
  }

  // ZDCC = vQarray[0], ZDCA = vQarray[1], see AliFlowEventSimple
  ZDCevArr[0] = TVector2::Phi_mpi_pi(vQarray[0].Phi());
  ZDCevArr[1] = TVector2::Phi_mpi_pi(vQarray[1].Phi());

  return kTRUE;
}

void AliAnalysisTaskMLTreeMakerEff::SetupTrackCuts(AliDielectronCutGroup* f)
{
filter   = new AliAnalysisFilter("filter","filter");  
filter->AddCuts(f);
}


void AliAnalysisTaskMLTreeMakerEff::SetupEventCuts(AliDielectronEventCuts* f)
{
  evfilter   = new AliAnalysisFilter("evfilter","evfilter");  
  evfilter->AddCuts(f);
}


int AliAnalysisTaskMLTreeMakerEff::CheckGenerator(Int_t trackID){     //check if the generator is on the list of generators
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
      if (genname.Hash() == fGeneratorHashs[i]){ 
//        std::cout <<"GenName acc: "<< genname << std::endl;
        return i;
      }
    }
//    std::cout <<"GenName rejected: "<< genname << std::endl;
    return -3;
  }
  return -4; // should not happen
}

  
void AliAnalysisTaskMLTreeMakerEff::SetupTMVAReader(TString weightFile){
  
  TMVAReader = new TMVA::Reader( "!Color:!Silent" );

  TMVAReader->AddVariable( "nITS", &nITSTMVA);
  TMVAReader->AddVariable( "ITS1Shared", &ITS1SharedTMVA);
  TMVAReader->AddVariable( "ITS2Shared", &ITS2SharedTMVA);
  TMVAReader->AddVariable( "ITS3Shared", &ITS3SharedTMVA);
  TMVAReader->AddVariable( "ITS4Shared", &ITS4SharedTMVA);
  TMVAReader->AddVariable( "ITS5Shared", &ITS5SharedTMVA);
  TMVAReader->AddVariable( "ITS6Shared", &ITS6SharedTMVA);
  TMVAReader->AddVariable( "nITSshared_frac", &nITSshared_fracTMVA);
//  TMVAReader->AddVariable( "NCrossedRowsTPC", &NCrossedRowsTPCTMVA);
//  TMVAReader->AddVariable( "NClustersTPC", &NClustersTPCTMVA);
//  TMVAReader->AddVariable( "NTPCSignal", &NTPCSignalTMVA);
  TMVAReader->AddVariable( "log(abs(DCAxy))", &logDCAxyTMVA );
  TMVAReader->AddVariable( "log(abs(DCAz))", &logDCAzTMVA );   
  TMVAReader->AddVariable( "chi2GlobalPerNDF", &chi2GlobalPerNDFTMVA);
  TMVAReader->AddVariable( "chi2ITS", &chi2ITSTMVA);
  TMVAReader->AddVariable( "eta", &etaTMVA);
//  TMVAReader->AddVariable( "phi", &phiTMVA);
  TMVAReader->AddVariable( "pt", &ptTMVA);   
  TMVAReader->AddVariable( "centrality", &centTMVA);

  gSystem->Exec(Form("alien_cp alien:///alice/cern.ch/user/s/selehner/TMVAweights/%s .",weightFile.Data()));
  std::cout<<"Setting weights file: "<<weightFile.Data()<<std::endl;

  TMVAReader->BookMVA( "BDTG method", weightFile.Data() );

}

