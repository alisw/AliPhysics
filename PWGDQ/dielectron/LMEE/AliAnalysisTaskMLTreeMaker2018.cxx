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
#include "AliAnalysisTaskMLTreeMaker2018.h"
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


// Authors: Jerome Jung (IKF Frankfurt) - jerome.jung@cern.ch


ClassImp(AliAnalysisTaskMLTreeMaker2018)


AliAnalysisTaskMLTreeMaker2018::AliAnalysisTaskMLTreeMaker2018():
  AliAnalysisTaskSE(),
  eventCuts(0),
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
  NClustersITS(0),
  NCrossedRowsTPC(0),
  NClustersTPC(0),
  HasSPDfirstHit(0),
  RatioCrossedRowsFindableClusters(0),
  NTPCSignal(0),
  fGeneratorHashs(0x0),
  fGeneratorName(""),
  fGeneratorMCSignalName(""),
  fGeneratorULSSignalName(""),
  fGeneratorMCSignalHashs(0x0),
  fGeneratorULSSignalHashs(0x0),
  runNumber(0),
  nTracks(0),
  nPairs({0,0,0}),
  cent(0),
  man(0),
  fList(0x0),
  fCentralityPercentileMin(0),
  fCentralityPercentileMax(100),
  fPtMin(0),
  fPtMax(1000),
  fEtaMin(-10),
  fEtaMax(10),
  mcEvent(0x0),
  fEvent(0x0),

  isAOD(kTRUE),
  hasMC(kFALSE),
  filterQuality(kFALSE),
  doPairing(kTRUE),
  doULS(kTRUE),
  doLS(kTRUE),
  doSingleLegMCSignal(kFALSE),
  fSingleLegMCSignal(0),
  fPairMCSignal(0),
  fNegPart(0),
  fPosPart(0),
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
  dcaXY(),
  dcaZ(),
  dcaXY_res(),
  dcaZ_res(),
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
  pdg(0),
  pdgmother(0),
  hasmother(0),
  label(0),
  motherlabel(0),
  fTreeTracks(0),
  fTreeULS(0),
  fTreeLSmm(0),
  fTreeLSpp(0),
  fTreePairSignal(0),
  fTreePairs({}),
  fTreePairSignals(0),
  fQAHistEvents(0),
  fQAHistTracks(0),
  eta_tracks1({}),
  phi_tracks1({}),
  pt_tracks1({}),
  charge_tracks1({}),
  MCpt_tracks1({}),
  MCeta_tracks1({}),
  MCphi_tracks1({}),
  MCvertx_tracks1({}),
  MCverty_tracks1({}),
  MCvertz_tracks1({}),
  dcaXY_tracks1({}),
  dcaZ_tracks1({}),
  dcaXY_res_tracks1({}),
  dcaZ_res_tracks1({}),
  verticesX_tracks1({}),
  verticesY_tracks1({}),
  verticesZ_tracks1({}),
  eta_tracks2({}),
  phi_tracks2({}),
  pt_tracks2({}),
  charge_tracks2({}),
  MCpt_tracks2({}),
  MCeta_tracks2({}),
  MCphi_tracks2({}),
  MCvertx_tracks2({}),
  MCverty_tracks2({}),
  MCvertz_tracks2({}),
  dcaXY_tracks2({}),
  dcaZ_tracks2({}),
  dcaXY_res_tracks2({}),    
  dcaZ_res_tracks2({}),
  verticesX_tracks2({}),
  verticesY_tracks2({}),
  verticesZ_tracks2({}),
  primVerticesX({}),
  primVerticesY({}),
  primVerticesZ({}),
  pairMasses({}),
  pairPtees({}),
  pairOpAngs({}),
  pairDCAee({}),
 
  //pairPhiV({}),
  pairCosPointAngs({}),
  pairDecayLengths({}),
  pairRs({}),
  



  fFillPairSignalOffset(0)

{

}

AliAnalysisTaskMLTreeMaker2018::AliAnalysisTaskMLTreeMaker2018(const char *name):
  AliAnalysisTaskSE(name),
  eventCuts(0),
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
  NClustersITS(0),
  NClustersTPC(0),
  HasSPDfirstHit(0),
  RatioCrossedRowsFindableClusters(0),
  NTPCSignal(0),
  fGeneratorHashs(0x0),
  fGeneratorName(""),
  fGeneratorMCSignalName(""),
  fGeneratorULSSignalName(""),
  fGeneratorMCSignalHashs(0x0),
  fGeneratorULSSignalHashs(0x0),
  runNumber(0),
  nTracks(0),
  nPairs({0,0,0}),
  //nPairs(0),
  cent(0),
  man(0),
  fList(0x0),
  fCentralityPercentileMin(0),
  fCentralityPercentileMax(100),
  fPtMin(0),
  fPtMax(1000),
  fEtaMin(-10),
  fEtaMax(10),
  mcEvent(0x0),
  fEvent(0x0),

  isAOD(kTRUE),
  hasMC(kFALSE),
  filterQuality(kFALSE),
  doPairing(kTRUE),
  doULS(kTRUE),
  doLS(kTRUE),
  doSingleLegMCSignal(kFALSE),
  fSingleLegMCSignal(0),
  fPairMCSignal(0),
  fNegPart(0),
  fPosPart(0),
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
  dcaXY(),
  dcaZ(),
  dcaXY_res(),
  dcaZ_res(),
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
  pdg(0),
  pdgmother(0),
  hasmother(0),
  label(0),
  motherlabel(0),
  fTreeTracks(0),
  fTreeULS(0),
  fTreeLSmm(0),
  fTreeLSpp(0),
  fTreePairSignal(0),
  fTreePairs({}),
  fTreePairSignals(0),
  fQAHistEvents(0),
  fQAHistTracks(0),

  eta_tracks1({}),
  phi_tracks1({}),
  pt_tracks1({}),
  charge_tracks1({}),
  MCpt_tracks1({}),
  MCeta_tracks1({}),
  MCphi_tracks1({}),
  MCvertx_tracks1({}),
  MCverty_tracks1({}),
  MCvertz_tracks1({}),
  dcaXY_tracks1({}),
  dcaZ_tracks1({}),
  dcaXY_res_tracks1({}),
  dcaZ_res_tracks1({}),
  verticesX_tracks1({}),
  verticesY_tracks1({}),
  verticesZ_tracks1({}),
  eta_tracks2({}),
  phi_tracks2({}),
  pt_tracks2({}),
  charge_tracks2({}),
  MCpt_tracks2({}),
  MCeta_tracks2({}),
  MCphi_tracks2({}),
  MCvertx_tracks2({}),
  MCverty_tracks2({}),
  MCvertz_tracks2({}),
  dcaXY_tracks2({}),
  dcaZ_tracks2({}),
  dcaXY_res_tracks2({}),    
  dcaZ_res_tracks2({}),
  verticesX_tracks2({}),
  verticesY_tracks2({}),
  verticesZ_tracks2({}),
  primVerticesX({}),
  primVerticesY({}),
  primVerticesZ({}),
  pairMasses({}),
  pairPtees({}),
  pairOpAngs({}),
  pairDCAee({}),

  pairCosPointAngs({}),
  pairDecayLengths({}),
  pairRs({}),

  fFillPairSignalOffset(0)
{

  DefineOutput(1, TList::Class());
}

AliAnalysisTaskMLTreeMaker2018::~AliAnalysisTaskMLTreeMaker2018(){
//  delete eventCuts;
//  delete evfilter;
//  
//  delete trcuts;
//  delete trfilter;
//  delete pidcuts;
//  delete cuts;
//  delete filter; 
//
//  delete fList;
//  delete fQAHist;
//  delete fTree;


}

void AliAnalysisTaskMLTreeMaker2018::UserCreateOutputObjects() {


  man=AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
  inputHandler->SetNeedField();


  fList = new TList();
  fList->SetName("output_Tlist");
  fList->SetOwner();

  AliInfo("Try to get PIDResponse");
  fPIDResponse = inputHandler->GetPIDResponse();

  if (!fPIDResponse){
    AliError("Failed to get PIDResponse - return");
    return;
  }

  if(hasMC) std::cout <<"Running on MC!"<< std::endl;
  else std::cout <<"Running on RD!"<< std::endl;

  if(filterQuality) std::cout <<"Writing track quality variables!"<< std::endl;

  Bool_t oldStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);

  fTreeTracks = new TTree("Track_Tree","Tracks");
  if(doTracks) fList->Add(fTreeTracks);

  fQAHistEvents = new TH1F("QAevents", "QA for events", 4, 0, 1);
  fList->Add(fQAHistEvents);
  fQAHistTracks = new TH1F("QAtracks", "QA for tracks", 4, 0, 1);
  fList->Add(fQAHistTracks);

  fTreeTracks->Branch("RunNumber", &runNumber);
  fTreeTracks->Branch("centrality", &cent);
  fTreeTracks->Branch("#tracks", &nTracks);

  fTreeTracks->Branch("pt", &pt);
  fTreeTracks->Branch("eta", &eta);
  fTreeTracks->Branch("phi", &phi);
  fTreeTracks->Branch("charge", &charge);

  fTreeTracks->Branch("DCAxy", &dcaXY);
  fTreeTracks->Branch("DCAz", &dcaZ);

  fTreeTracks->Branch("DCAxy_res", &dcaXY_res);
  fTreeTracks->Branch("DCAz_res", &dcaZ_res);

  fTreeTracks->Branch("vertx", &vertx);
  fTreeTracks->Branch("verty", &verty);
  fTreeTracks->Branch("vertz", &vertz);

  if(filterQuality){
    fTreeTracks->Branch("chi2GlobalvsTPC", &chi2GlobalvsTPC);
    fTreeTracks->Branch("chi2GlobalPerNDF", &chi2GlobalPerNDF);


    fTreeTracks->Branch("nITS", &nITS);
    fTreeTracks->Branch("ITS1Shared", &ITS1S);
    fTreeTracks->Branch("ITS2Shared", &ITS2S);
    fTreeTracks->Branch("ITS3Shared", &ITS3S);
    fTreeTracks->Branch("ITS4Shared", &ITS4S);
    fTreeTracks->Branch("ITS5Shared", &ITS5S);
    fTreeTracks->Branch("ITS6Shared", &ITS6S);
    fTreeTracks->Branch("nITSshared_frac", &nITSshared);

    fTreeTracks->Branch("chi2ITS", &chi2ITS);
    fTreeTracks->Branch("chi2TPC", &chi2TPC);


    //probably not needed 
    fTreeTracks->Branch("NCrossedRowsTPC", &NCrossedRowsTPC);
    fTreeTracks->Branch("NClustersTPC", &NClustersTPC);
    fTreeTracks->Branch("RatioCrossedRowsFindableClusters", &RatioCrossedRowsFindableClusters);
    fTreeTracks->Branch("HasSPDfirstHit", &HasSPDfirstHit);
  }





  fTreeTracks->SetAutoSave(100000000);
  if(doPairing){

    Int_t nPairTrees = 0;
    if(doULS) nPairTrees++;
    if(doLS) nPairTrees += 2;

    fFillPairSignalOffset = nPairTrees;

    if(fPairMCSignal.size()>0) nPairTrees += fPairMCSignal.size();

    for(int pairtree = 0; pairtree<nPairTrees; pairtree++){
      
      primVerticesX.push_back({0});
      primVerticesY.push_back({0});
      primVerticesZ.push_back({0});

      pt_tracks1.push_back({0});
      eta_tracks1.push_back({0});
      phi_tracks1.push_back({0});
      charge_tracks1.push_back({0});
      dcaXY_tracks1.push_back({0});
      dcaZ_tracks1.push_back({0});
      dcaXY_res_tracks1.push_back({0});
      dcaZ_res_tracks1.push_back({0});
      verticesX_tracks1.push_back({0});
      verticesY_tracks1.push_back({0});
      verticesZ_tracks1.push_back({0});


      pt_tracks2.push_back({0});
      eta_tracks2.push_back({0});
      phi_tracks2.push_back({0});
      charge_tracks2.push_back({0});
      dcaXY_tracks2.push_back({0});
      dcaZ_tracks2.push_back({0});
      dcaXY_res_tracks2.push_back({0});
      dcaZ_res_tracks2.push_back({0});
      verticesX_tracks2.push_back({0});
      verticesY_tracks2.push_back({0});
      verticesZ_tracks2.push_back({0});
  
      nPairs.push_back(0);
      pairMasses.push_back({0});
      pairPtees.push_back({0});
      pairOpAngs.push_back({0});
      pairDCAee.push_back({0});

      pairCosPointAngs.push_back({0});
      pairDecayLengths.push_back({0});
      pairRs.push_back({0});


    }


    if(doULS) {
      fTreeULS = new TTree("Pair_Tree_ULS","Pairs_ULS");
      fList->Add(fTreeULS);
      fTreePairs.push_back(fTreeULS);
    }

    if(doLS) {
      fTreeLSmm = new TTree("Pair_Tree_LSmm","Pairs_LSmm");
      fList->Add(fTreeLSmm);
      fTreePairs.push_back(fTreeLSmm);

      fTreeLSpp = new TTree("Pair_Tree_LSpp","Pairs_LSpp");
      fList->Add(fTreeLSpp);
      fTreePairs.push_back(fTreeLSpp);
    }

    for(int signal = 0; signal<fPairMCSignal.size(); signal++){

      fTreePairSignal = new TTree(Form("PairSignal_Tree_%s", fPairMCSignal.at(signal).GetName()),Form("Pairs_%s", fPairMCSignal.at(signal).GetName()));
      fList->Add(fTreePairSignal);
      fTreePairSignals.push_back(fTreePairSignal);
      fTreePairs.push_back(fTreePairSignal);
    }


    for (int index=0; index<fTreePairs.size(); index++){
      fTreePairs.at(index)->Branch("#Pairs", &nPairs.at(index));    

      fTreePairs.at(index)->Branch("pt_track1", &pt_tracks1.at(index));
      fTreePairs.at(index)->Branch("eta_track1", &eta_tracks1.at(index));
      fTreePairs.at(index)->Branch("phi_track1", &phi_tracks1.at(index));
      fTreePairs.at(index)->Branch("charge_track1", &charge_tracks1.at(index));
      fTreePairs.at(index)->Branch("DCAxy_track1", &dcaXY_tracks1.at(index));
      fTreePairs.at(index)->Branch("DCAz_track1", &dcaZ_tracks1.at(index));
      fTreePairs.at(index)->Branch("DCAxy_res_track1", &dcaXY_res_tracks1.at(index));
      fTreePairs.at(index)->Branch("DCAz_res_track1", &dcaZ_res_tracks1.at(index));
      fTreePairs.at(index)->Branch("vertx_track1", &verticesX_tracks1.at(index));
      fTreePairs.at(index)->Branch("verty_track1", &verticesY_tracks1.at(index));
      fTreePairs.at(index)->Branch("vertz_track1", &verticesZ_tracks1.at(index));

      fTreePairs.at(index)->Branch("pt_track2", &pt_tracks2.at(index));
      fTreePairs.at(index)->Branch("eta_track2", &eta_tracks2.at(index));
      fTreePairs.at(index)->Branch("phi_track2", &phi_tracks2.at(index));
      fTreePairs.at(index)->Branch("charge_track2", &charge_tracks2.at(index));
      fTreePairs.at(index)->Branch("DCAxy_track2", &dcaXY_tracks2.at(index));
      fTreePairs.at(index)->Branch("DCAz_track2", &dcaZ_tracks2.at(index));
      fTreePairs.at(index)->Branch("DCAxy_res_track2", &dcaXY_res_tracks2.at(index));
      fTreePairs.at(index)->Branch("DCAz_res_track2", &dcaZ_res_tracks2.at(index));
      fTreePairs.at(index)->Branch("vertx_track2", &verticesX_tracks2.at(index));
      fTreePairs.at(index)->Branch("verty_track2", &verticesY_tracks2.at(index));
      fTreePairs.at(index)->Branch("vertz_track2", &verticesZ_tracks2.at(index));

      fTreePairs.at(index)->Branch("pairMass", &pairMasses.at(index));
      fTreePairs.at(index)->Branch("pairPtee", &pairPtees.at(index));
      fTreePairs.at(index)->Branch("pairOpAng", &pairOpAngs.at(index));
      fTreePairs.at(index)->Branch("pairDCAee", &pairDCAee.at(index));


      fTreePairs.at(index)->Branch("pairCosPointAng", &pairCosPointAngs.at(index));
      fTreePairs.at(index)->Branch("pairDecayLength", &pairDecayLengths.at(index));
      fTreePairs.at(index)->Branch("pairR", &pairRs.at(index));

      fTreePairs.at(index)->Branch("primVertx", &primVerticesX.at(index));
      fTreePairs.at(index)->Branch("primVerty", &primVerticesY.at(index));
      fTreePairs.at(index)->Branch("primVertz", &primVerticesZ.at(index));

      fTreePairs.at(index)->SetAutoSave(100000000);
    }
  }
  else{
    
    if(hasMC) {

      fTreeTracks->Branch("Pdg", &pdg);
      fTreeTracks->Branch("Label", &glabel);

      fTreeTracks->Branch("Pdg_Mother", &pdgmother);
      fTreeTracks->Branch("Mother_label", &motherlabel);
      fTreeTracks->Branch("Has_Mother", &hasmother);

      fTreeTracks->Branch("LabelFirstMother", &gLabelFirstMother);
      fTreeTracks->Branch("LabelMinFirstMother", &gLabelMinFirstMother);
      fTreeTracks->Branch("LabelMaxFirstMother", &gLabelMaxFirstMother);
      fTreeTracks->Branch("PdgFirstMother", &iPdgFirstMother);


      fTreeTracks->Branch("MCpt", &MCpt);
      fTreeTracks->Branch("MCeta", &MCeta);
      fTreeTracks->Branch("MCphi", &MCphi);

      fTreeTracks->Branch("MCTrack_vertx", &MCvertx);
      fTreeTracks->Branch("MCTrack_verty", &MCverty);
      fTreeTracks->Branch("MCTrack_vertz", &MCvertz);

      fTreeTracks->Branch("GenIndex", &iGenIndex);

    }
  }



  PostData(1, fList);

  AliInfo("Finished setting up the Output");
  TH1::AddDirectory(oldStatus);


  if(hasMC){

    TObjArray arr = *(fGeneratorName.Tokenize(";"));
    //std::cout << "Used Generators: " << std::endl;
    for (int i = 0; i < arr.GetEntries(); ++i){
      TString temp = arr.At(i)->GetName();
      fGeneratorHashs.push_back(temp.Hash());
    }
  
    arr = *(fGeneratorMCSignalName.Tokenize(";"));
    //std::cout << "Used Generators for MCSignals: " << std::endl;
    for (int i = 0; i < arr.GetEntries(); ++i){
      TString temp = arr.At(i)->GetName();
      fGeneratorMCSignalHashs.push_back(temp.Hash());
    }
  
    arr = *(fGeneratorULSSignalName.Tokenize(";"));
    //std::cout << "Used Generators for ULSSignals: " << std::endl;
    for (int i = 0; i < arr.GetEntries(); ++i){
      TString temp = arr.At(i)->GetName();
      fGeneratorULSSignalHashs.push_back(temp.Hash());
    }
  }
}
//________________________________________________________________________

void AliAnalysisTaskMLTreeMaker2018::UserExec(Option_t * option) {
  // Called for each event

  // ##########################################################
  // Set MC event
  if(!AliDielectronMC::Instance()->ConnectMCEvent()) return;

  // ##########################################################

  fQAHistEvents->Fill("Events_all",1);

  AliVEvent* event = dynamic_cast<AliVEvent*>(InputEvent());
  if(!event) {
    AliError("event not available");
    fQAHistEvents->Fill("Events_not_available",1);
    return;
  }
  runNumber = event->GetRunNumber();


  AliInputEventHandler *eventHandler = nullptr;
  AliInputEventHandler *eventHandlerMC = nullptr;
  if ((AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler())->IsA() == AliAODInputHandler::Class()){
    eventHandler = dynamic_cast<AliAODInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
    eventHandlerMC = eventHandler;
  }
  else   if ((AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler())->IsA() == AliESDInputHandler::Class()){
    eventHandler = dynamic_cast<AliESDInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
    eventHandlerMC = eventHandler;
  }
  if (isAOD) fEvent = static_cast<AliAODEvent*>(eventHandler->GetEvent());
  else       fEvent = static_cast<AliESDEvent*>(eventHandler->GetEvent());


  AliMultSelection *MultSelection = 0x0;
  MultSelection = (AliMultSelection * ) event->FindListObject("MultSelection");


  if( !MultSelection) {
   //If you get this warning (and lPercentiles 300) please check that the AliMultSelectionTask actually ran (before your task)
   AliWarning("AliMultSelection object not found!");
   std::cout << "+++ No Multiplicity task +++" << std::endl;

  }
  else cent = MultSelection->GetMultiplicityPercentile("V0M",kFALSE);


  fQAHistEvents->Fill("Events before cent",1);
  if(fCentralityPercentileMax!=0 && (cent<fCentralityPercentileMin || cent>fCentralityPercentileMax)) return;
  fQAHistEvents->Fill("Events after cent",1);


  UInt_t selectedMask=(1<<evfilter->GetCuts()->GetEntries())-1;
  TBits* fUsedVars = new TBits(AliDielectronVarManager::kNMaxValues);
  fUsedVars->SetBitNumber(AliDielectronVarManager::kP, kTRUE);
  varManager->SetFillMap(fUsedVars);
  varManager->SetEvent(event);
  if(selectedMask!=(evfilter->IsSelected(event))){
    fQAHistEvents->Fill("Events_not_selected_filter",1);
    return;
  }

  fQAHistEvents->Fill("Events before get tracks",1);
  Double_t lMultiplicityVar = -1;
  Int_t acceptedTracks = GetAcceptedTracks(event,lMultiplicityVar);
  nTracks= acceptedTracks;
  
  std::tuple<int, int, int> acceptedPairs = GetAcceptedPairs(event,lMultiplicityVar);
  
  fQAHistEvents->Fill("Events after get tracks",1);
  if(acceptedTracks){
    if(doTracks){
      fTreeTracks->Fill();
      
    }
    fQAHistEvents->Fill("Events with tracks",1);
    if(nPairs[0]>0 && doULS){
      fTreeULS->Fill();
    }
    if(nPairs[1]>0 && doLS){
      fTreeLSmm->Fill();
    }
    if(nPairs[2]>0 && doLS){
      fTreeLSpp->Fill();
    }
    for(int i=0; i<fTreePairSignals.size(); i++){
      if(nPairs[i+fFillPairSignalOffset]>0){   
        fTreePairSignals[i]->Fill();
      }
    }
  }
  else fQAHistEvents->Fill("Events without tracks",1);

  PostData(1, fList);
}

//~ //________________________________________________________________________

void  AliAnalysisTaskMLTreeMaker2018::FinishTaskOutput(){
  // Finish task output

  // not implemented ...


}
//~

//~ //________________________________________________________________________

void AliAnalysisTaskMLTreeMaker2018::Terminate(Option_t *) {

}
//~


//________________________________________________________________________

//Double_t AliAnalysisTaskMLTreeMaker2018::IsEventAccepted(AliVEvent *event){
//
////  if(event->GetPrimaryVertex()){
////    if (TMath::Abs(event->GetPrimaryVertex()->GetZ()) < 10){
////      if (event->GetPrimaryVertexSPD()->GetNContributors() > 0)
////    return 1;
////    }
////  }
//  return 0;
//}

//________________________________________________________________________
Int_t AliAnalysisTaskMLTreeMaker2018::FillSimplePairs(std::vector<Particle> parts1, std::vector<Particle> parts2, AliVEvent *event, Int_t index){



  Double_t vert[3] = {0};
  event->GetPrimaryVertex()->GetXYZ(vert);


  std::vector<Float_t> eta_track1;
  std::vector<Float_t> phi_track1;
  std::vector<Float_t> pt_track1;
  std::vector<Int_t> charge_track1;

  std::vector<Float_t> MCpt_track1;
  std::vector<Float_t> MCeta_track1;
  std::vector<Float_t> MCphi_track1;

  std::vector<Float_t> MCvertx_track1;
  std::vector<Float_t> MCverty_track1;
  std::vector<Float_t> MCvertz_track1;

  std::vector<Float_t> dcaXY_track1;    //DCA
  std::vector<Float_t> dcaZ_track1;

  std::vector<Float_t> dcaXY_res_track1;    //DCA
  std::vector<Float_t> dcaZ_res_track1;

  std::vector<Float_t> vertx_track1;
  std::vector<Float_t> verty_track1;
  std::vector<Float_t> vertz_track1;


  std::vector<Float_t> eta_track2;
  std::vector<Float_t> phi_track2;
  std::vector<Float_t> pt_track2;
  std::vector<Int_t> charge_track2;

  std::vector<Float_t> MCpt_track2;
  std::vector<Float_t> MCeta_track2;
  std::vector<Float_t> MCphi_track2;

  std::vector<Float_t> MCvertx_track2;
  std::vector<Float_t> MCverty_track2;
  std::vector<Float_t> MCvertz_track2;

  std::vector<Float_t> dcaXY_track2;    //DCA
  std::vector<Float_t> dcaZ_track2;

  std::vector<Float_t> dcaXY_res_track2;    //DCA
  std::vector<Float_t> dcaZ_res_track2;

  std::vector<Float_t> vertx_track2;
  std::vector<Float_t> verty_track2;
  std::vector<Float_t> vertz_track2;

  std::vector<Float_t> primVertx;
  std::vector<Float_t> primVerty;
  std::vector<Float_t> primVertz;

  std::vector<Float_t> pairMass;
  std::vector<Float_t> pairPtee;
  std::vector<Float_t> pairOpAng;
  std::vector<Float_t> pairDCA;

  std::vector<Float_t> pairCosPointAng;
  std::vector<Float_t> pairDecayLength;
  std::vector<Float_t> pairR;


  Int_t filledPairs = 0;
  for (unsigned int pos_i = 0; pos_i < parts1.size(); ++pos_i){
    for (unsigned int pos_j = 0; pos_j < parts2.size(); ++pos_j){

      
      //if(pos_i == pos_j && parts1[pos_i].fCharge == parts2[pos_j].fCharge ) continue;
      if( pos_j==0 && parts1[pos_i].fCharge == parts2[pos_j].fCharge ){
        pos_j = pos_i+1;
        if(pos_j >= parts2.size()) break;
      } 

      AliDielectronPair candidate;
      candidate.SetKFUsage(kTRUE);
      candidate.SetTracks(static_cast<AliVTrack*>(event->GetTrack(parts1[pos_i].GetTrackID())), -11,
                          static_cast<AliVTrack*>(event->GetTrack(parts2[pos_j].GetTrackID())),  11);

      // Fill AliDielectronPair specific information
      const AliKFParticle &kfPair = candidate.GetKFParticle();
      float DecayLength = kfPair.GetDecayLengthXY();
      float R = kfPair.GetR();
      //kfPair.GetChi2()/kfPair.GetNDF();
      float CosPointAng = event ? candidate.GetCosPointingAngle(event->GetPrimaryVertex()) : -1;

      //fgEvent ? pair->PsiPair(fgEvent->GetMagneticField()) : -5;
      //fgEvent ? pair->PhivPair(fgEvent->GetMagneticField()) : -5;



      TLorentzVector Lvec1;
      TLorentzVector Lvec2;
      Lvec1.SetPtEtaPhiM(parts1[pos_i].fPt, parts1[pos_i].fEta, parts1[pos_i].fPhi, AliPID::ParticleMass(AliPID::kElectron));
      Lvec2.SetPtEtaPhiM(parts2[pos_j].fPt, parts2[pos_j].fEta, parts2[pos_j].fPhi, AliPID::ParticleMass(AliPID::kElectron));
      TLorentzVector LvecM = Lvec1 + Lvec2;
      float mass = LvecM.M();
      float pairpt = LvecM.Pt();
      float opangle = Lvec1.Angle(Lvec2.Vect());

      float dcaee = sqrt( (pow(parts1[pos_i].fDCAxy/sqrt(parts1[pos_i].fDCAxy_res),2) + pow(parts2[pos_j].fDCAxy/sqrt(parts2[pos_j].fDCAxy_res),2)) / 2  );

      eta_track1.push_back(parts1[pos_i].fEta);
      phi_track1.push_back(parts1[pos_i].fPhi);
      pt_track1.push_back(parts1[pos_i].fPt);
      charge_track1.push_back(parts1[pos_i].fCharge);

      dcaXY_track1.push_back(parts1[pos_i].fDCAxy);    //DCA
      dcaZ_track1.push_back(parts1[pos_i].fDCAz);

      dcaXY_res_track1.push_back(parts1[pos_i].fDCAxy_res);    //DCA
      dcaZ_res_track1.push_back(parts1[pos_i].fDCAz_res);

      vertx_track1.push_back(parts1[pos_i].fXv);
      verty_track1.push_back(parts1[pos_i].fYv);
      vertz_track1.push_back(parts1[pos_i].fZv);


      eta_track2.push_back(parts2[pos_j].fEta);
      phi_track2.push_back(parts2[pos_j].fPhi);
      pt_track2.push_back(parts2[pos_j].fPt);
      charge_track2.push_back(parts2[pos_j].fCharge);

      dcaXY_track2.push_back(parts2[pos_j].fDCAxy);    //DCA
      dcaZ_track2.push_back(parts2[pos_j].fDCAz);

      dcaXY_res_track2.push_back(parts2[pos_j].fDCAxy_res);    //DCA
      dcaZ_res_track2.push_back(parts2[pos_j].fDCAz_res);

      vertx_track2.push_back(parts2[pos_j].fXv);
      verty_track2.push_back(parts2[pos_j].fYv);
      vertz_track2.push_back(parts2[pos_j].fZv);

      primVertx.push_back((Float_t)vert[0]);
      primVerty.push_back((Float_t)vert[1]);
      primVertz.push_back((Float_t)vert[2]);

      pairMass.push_back(mass);
      pairPtee.push_back(pairpt);
      pairOpAng.push_back(opangle);
      pairDCA.push_back(dcaee);
      
      pairCosPointAng.push_back(CosPointAng);
      pairDecayLength.push_back(DecayLength);
      pairR.push_back(R);



      filledPairs++;
    }
  }



  if(filledPairs>0){

    nPairs[index] = filledPairs;

    eta_tracks1[index] = eta_track1;
    phi_tracks1[index] = phi_track1;
    pt_tracks1[index] = pt_track1;
    charge_tracks1[index] = charge_track1;

    dcaXY_tracks1[index] = dcaXY_track1;    //DCA
    dcaZ_tracks1[index] = dcaZ_track1;

    dcaXY_res_tracks1[index] = dcaXY_res_track1;    //DCA
    dcaZ_res_tracks1[index] = dcaZ_res_track1;

    verticesX_tracks1[index] = vertx_track1;
    verticesY_tracks1[index] = verty_track1;
    verticesZ_tracks1[index] = vertz_track1;


    eta_tracks2[index] = eta_track2;
    phi_tracks2[index] = phi_track2;
    pt_tracks2[index] = pt_track2;
    charge_tracks2[index] = charge_track2;

    dcaXY_tracks2[index] = dcaXY_track2;    //DCA
    dcaZ_tracks2[index] = dcaZ_track2;

    dcaXY_res_tracks2[index] = dcaXY_res_track2;    //DCA
    dcaZ_res_tracks2[index] = dcaZ_res_track2;

    verticesX_tracks2[index] = vertx_track2;
    verticesY_tracks2[index] = verty_track2;
    verticesZ_tracks2[index] = vertz_track2;

    primVerticesX[index] = primVertx;
    primVerticesY[index] = primVerty;
    primVerticesZ[index] = primVertz;

    pairMasses[index] = pairMass;
    pairPtees[index] = pairPtee;
    pairOpAngs[index] = pairOpAng;
    pairDCAee[index] = pairDCA;

    pairCosPointAngs[index] = pairCosPointAng;
    pairDecayLengths[index] = pairDecayLength;
    pairRs[index] = pairR;

  }
 
  return filledPairs;
}

//________________________________________________________________________
Int_t AliAnalysisTaskMLTreeMaker2018::FillSignalPairs(AliVEvent *event){

  Double_t vert[3] = {0};
  event->GetPrimaryVertex()->GetXYZ(vert);

  //event variables
  std::vector<std::vector<Float_t>> primVertx;
  std::vector<std::vector<Float_t>> primVerty;
  std::vector<std::vector<Float_t>> primVertz;


  //std::vector<std::vector<Float_t>> MCvertx_track1;
  //std::vector<std::vector<Float_t>> MCverty_track1;
  //std::vector<std::vector<Float_t>> MCvertz_track1;
  //std::vector<std::vector<Float_t>> MCvertx_track2;
  //std::vector<std::vector<Float_t>> MCverty_track2;
  //std::vector<std::vector<Float_t>> MCvertz_track2;



  //leg1 variables
  std::vector<std::vector<Float_t>> eta_track1;
  std::vector<std::vector<Float_t>> phi_track1;
  std::vector<std::vector<Float_t>> pt_track1;
  std::vector<std::vector<Int_t>> charge_track1;
  //std::vector<std::vector<Float_t>> MCpt_track1;
  //std::vector<std::vector<Float_t>> MCeta_track1;
  //std::vector<std::vector<Float_t>> MCphi_track1;
  std::vector<std::vector<Float_t>> dcaXY_track1;    //DCA
  std::vector<std::vector<Float_t>> dcaZ_track1;
  std::vector<std::vector<Float_t>> dcaXY_res_track1;    //DCA
  std::vector<std::vector<Float_t>> dcaZ_res_track1;
  std::vector<std::vector<Float_t>> vertx_track1;
  std::vector<std::vector<Float_t>> verty_track1;
  std::vector<std::vector<Float_t>> vertz_track1;



  //leg2 variables
  std::vector<std::vector<Float_t>> eta_track2;
  std::vector<std::vector<Float_t>> phi_track2;
  std::vector<std::vector<Float_t>> pt_track2;
  std::vector<std::vector<Int_t>> charge_track2;
  //std::vector<std::vector<Float_t>> MCpt_track2;
  //std::vector<std::vector<Float_t>> MCeta_track2;
  //std::vector<std::vector<Float_t>> MCphi_track2;
  std::vector<std::vector<Float_t>> dcaXY_track2;    //DCA
  std::vector<std::vector<Float_t>> dcaZ_track2;
  std::vector<std::vector<Float_t>> dcaXY_res_track2;    //DCA
  std::vector<std::vector<Float_t>> dcaZ_res_track2;
  std::vector<std::vector<Float_t>> vertx_track2;
  std::vector<std::vector<Float_t>> verty_track2;
  std::vector<std::vector<Float_t>> vertz_track2;


  //pair variables
  std::vector<Int_t> filledPairMCSignal;
  std::vector<std::vector<Float_t>> pairMass;
  std::vector<std::vector<Float_t>> pairPtee;
  std::vector<std::vector<Float_t>> pairOpAng;
  std::vector<std::vector<Float_t>> pairDCA;

  std::vector<std::vector<Float_t>> pairCosPointAng;
  std::vector<std::vector<Float_t>> pairDecayLength;
  std::vector<std::vector<Float_t>> pairR;


  for(int signal=0; signal < fPairMCSignal.size(); signal++){

    primVertx.push_back({});
    primVerty.push_back({});
    primVertz.push_back({});
 
    eta_track1.push_back({});
    phi_track1.push_back({});
    pt_track1.push_back({});
    charge_track1.push_back({});
    dcaXY_track1.push_back({});    
    dcaZ_track1.push_back({});   
    dcaXY_res_track1.push_back({});
    dcaZ_res_track1.push_back({});
    vertx_track1.push_back({});
    verty_track1.push_back({});
    vertz_track1.push_back({});
   
    eta_track2.push_back({});
    phi_track2.push_back({});
    pt_track2.push_back({});
    charge_track2.push_back({});
    dcaXY_track2.push_back({});    
    dcaZ_track2.push_back({});                  
    dcaXY_res_track2.push_back({});    
    dcaZ_res_track2.push_back({});
    vertx_track2.push_back({});
    verty_track2.push_back({});
    vertz_track2.push_back({});
   
                              
    filledPairMCSignal.push_back(0);                          
    pairMass.push_back({});
    pairPtee.push_back({});
    pairOpAng.push_back({});
    pairDCA.push_back({});

    pairCosPointAng.push_back({});
    pairDecayLength.push_back({});
    pairR.push_back({});

  }

  for (unsigned int neg_i = 0; neg_i < fNegPart.size(); ++neg_i){
    for (unsigned int pos_i = 0; pos_i < fPosPart.size(); ++pos_i){
      AliVParticle* mcPart1 = fMCEvent->GetTrack(fNegPart[neg_i].GetTrackLabel());
      AliVParticle* mcPart2 = fMCEvent->GetTrack(fPosPart[pos_i].GetTrackLabel());

      AliDielectronPair candidate;
      candidate.SetKFUsage(kFALSE);
      candidate.SetTracks(static_cast<AliVTrack*>(event->GetTrack(fNegPart[neg_i].GetTrackID())), -11,
                          static_cast<AliVTrack*>(event->GetTrack(fPosPart[pos_i].GetTrackID())),  11);

      // Fill AliDielectronPair specific information
      const AliKFParticle &kfPair = candidate.GetKFParticle();
      float DecayLength = kfPair.GetDecayLength();
      float R = kfPair.GetR();
      //kfPair.GetChi2()/kfPair.GetNDF();
      float CosPointAng = event ? candidate.GetCosPointingAngle(event->GetPrimaryVertex()) : -1;

      //fgEvent ? pair->PsiPair(fgEvent->GetMagneticField()) : -5;
      //fgEvent ? pair->PhivPair(fgEvent->GetMagneticField()) : -5;


      // Check if electrons are from MCSignal Generator
      if (!fPosPart[pos_i].GetMCSignalPair() || !fNegPart[neg_i].GetMCSignalPair()){
	continue;
      }
      // Apply MC signals
      std::vector<Bool_t> mcSignal_acc(fPairMCSignal.size(), kTRUE); // vector which stores if track is accepted by [i]-th mcsignal

      // Check if it according to mcsignals
      for (unsigned int i = 0; i < fPairMCSignal.size(); ++i){
        mcSignal_acc[i] = AliDielectronMC::Instance()->IsMCTruth(mcPart1, mcPart2, &(fPairMCSignal[i]));
      }
      // check if at least one mc signal is true
      if (CheckIfOneIsTrue(mcSignal_acc) == kFALSE) continue;
      // Construct pair variables from LorentzVectors
      TLorentzVector Lvec1;
      TLorentzVector Lvec2;
      Lvec1.SetPtEtaPhiM(fNegPart[neg_i].fPt, fNegPart[neg_i].fEta, fNegPart[neg_i].fPhi, AliPID::ParticleMass(AliPID::kElectron));
      Lvec2.SetPtEtaPhiM(fPosPart[pos_i].fPt, fPosPart[pos_i].fEta, fPosPart[pos_i].fPhi, AliPID::ParticleMass(AliPID::kElectron));
      //Lvec1.SetPtEtaPhiM(mcPart1->Pt(), mcPart1->Eta(), mcPart1->Phi(), AliPID::ParticleMass(AliPID::kElectron));
      //Lvec2.SetPtEtaPhiM(mcPart2->Pt(), mcPart2->Eta(), mcPart2->Phi(), AliPID::ParticleMass(AliPID::kElectron));
      TLorentzVector LvecM = Lvec1 + Lvec2;
      double mass = LvecM.M();
      double pairpt = LvecM.Pt();
      double opangle = Lvec1.Angle(Lvec2.Vect());
      double weight = 1.;

      float dcaee = sqrt( (pow(fNegPart[neg_i].fDCAxy/sqrt(fNegPart[neg_i].fDCAxy_res),2) + pow(fPosPart[pos_i].fDCAxy/sqrt(fPosPart[pos_i].fDCAxy_res),2)) / 2  );

      //std::cout <<  " &&&&&& new pair: " << mass << "  " << candidate.M()  << std::endl;
      //std::cout <<  " ------ pointing ang: " << candidate.GetCosPointingAngle(event->GetPrimaryVertex())  << std::endl;
      //std::cout <<  " ------ radius: " << candidate.GetKFParticle().GetR()  << std::endl;
      //std::cout <<  " ------ decay length: " << candidate.GetKFParticle().GetDecayLength()  << std::endl;



      for (unsigned int i = 0; i < mcSignal_acc.size(); ++i){
        if (mcSignal_acc[i] == kTRUE){
          filledPairMCSignal[i] = filledPairMCSignal[i] + 1;	
          //pt_track1[i].push_back(mcPart1->Pt());

          eta_track1[i].push_back(fNegPart[neg_i].fEta);
          phi_track1[i].push_back(fNegPart[neg_i].fPhi);
          pt_track1[i].push_back(fNegPart[neg_i].fPt);
          charge_track1[i].push_back(fNegPart[neg_i].fCharge);

          dcaXY_track1[i].push_back(fNegPart[neg_i].fDCAxy);    //DCA
          dcaZ_track1[i].push_back(fNegPart[neg_i].fDCAz);

          dcaXY_res_track1[i].push_back(fNegPart[neg_i].fDCAxy_res);    //DCA
          dcaZ_res_track1[i].push_back(fNegPart[neg_i].fDCAz_res);

	  vertx_track1[i].push_back(fNegPart[neg_i].fXv);
          verty_track1[i].push_back(fNegPart[neg_i].fYv);
          vertz_track1[i].push_back(fNegPart[neg_i].fZv);
    

          eta_track2[i].push_back(fPosPart[pos_i].fEta);
          phi_track2[i].push_back(fPosPart[pos_i].fPhi);
          pt_track2[i].push_back(fPosPart[pos_i].fPt);
          charge_track2[i].push_back(fPosPart[pos_i].fCharge);

          dcaXY_track2[i].push_back(fPosPart[pos_i].fDCAxy);    //DCA
          dcaZ_track2[i].push_back(fPosPart[pos_i].fDCAz);

          dcaXY_res_track2[i].push_back(fPosPart[pos_i].fDCAxy_res);    //DCA
          dcaZ_res_track2[i].push_back(fPosPart[pos_i].fDCAz_res);

          vertx_track2[i].push_back(fPosPart[pos_i].fXv);
          verty_track2[i].push_back(fPosPart[pos_i].fYv);
          vertz_track2[i].push_back(fPosPart[pos_i].fZv);
    
          primVertx[i].push_back((Float_t)vert[0]);
          primVerty[i].push_back((Float_t)vert[1]);
          primVertz[i].push_back((Float_t)vert[2]);



          pairMass[i].push_back(mass);
          pairPtee[i].push_back(pairpt);
          pairOpAng[i].push_back(opangle);
	  pairDCA[i].push_back(dcaee);

          pairCosPointAng[i].push_back(CosPointAng);
          pairDecayLength[i].push_back(DecayLength);
          pairR[i].push_back(R);

        }
      } // end of loop over all MCsignals
    } // end of loop over all positive particles
  } // end of loop over all negative particles

  for(int signal=fFillPairSignalOffset; signal < fPairMCSignal.size() + fFillPairSignalOffset; signal++){ 

    //event variables
    primVerticesX[signal] = primVertx[signal-fFillPairSignalOffset];
    primVerticesY[signal] = primVerty[signal-fFillPairSignalOffset];
    primVerticesZ[signal] = primVertz[signal-fFillPairSignalOffset];

    //leg1 variables
    eta_tracks1[signal] = eta_track1[signal-fFillPairSignalOffset];
    phi_tracks1[signal] = phi_track1[signal-fFillPairSignalOffset];
    pt_tracks1[signal] = pt_track1[signal-fFillPairSignalOffset];
    charge_tracks1[signal] = charge_track1[signal-fFillPairSignalOffset];
    dcaXY_tracks1[signal] = dcaXY_track1[signal-fFillPairSignalOffset];    //DCA
    dcaZ_tracks1[signal] = dcaZ_track1[signal-fFillPairSignalOffset];
    dcaXY_res_tracks1[signal] = dcaXY_res_track1[signal-fFillPairSignalOffset];    //DCA
    dcaZ_res_tracks1[signal] = dcaZ_res_track1[signal-fFillPairSignalOffset];
    verticesX_tracks1[signal] = vertx_track1[signal-fFillPairSignalOffset];
    verticesY_tracks1[signal] = verty_track1[signal-fFillPairSignalOffset];
    verticesZ_tracks1[signal] = vertz_track1[signal-fFillPairSignalOffset];


    //leg2 variables
    eta_tracks2[signal] = eta_track2[signal-fFillPairSignalOffset];
    phi_tracks2[signal] = phi_track2[signal-fFillPairSignalOffset];
    pt_tracks2[signal] = pt_track2[signal-fFillPairSignalOffset];
    charge_tracks2[signal] = charge_track2[signal-fFillPairSignalOffset];
    dcaXY_tracks2[signal] = dcaXY_track2[signal-fFillPairSignalOffset];    //DCA
    dcaZ_tracks2[signal] = dcaZ_track2[signal-fFillPairSignalOffset];
    dcaXY_res_tracks2[signal] = dcaXY_res_track2[signal-fFillPairSignalOffset];    //DCA
    dcaZ_res_tracks2[signal] = dcaZ_res_track2[signal-fFillPairSignalOffset];
    verticesX_tracks2[signal] = vertx_track2[signal-fFillPairSignalOffset];
    verticesY_tracks2[signal] = verty_track2[signal-fFillPairSignalOffset];
    verticesZ_tracks2[signal] = vertz_track2[signal-fFillPairSignalOffset];

    //pair variables
    nPairs[signal] = filledPairMCSignal[signal-fFillPairSignalOffset];
    pairMasses[signal] = pairMass[signal-fFillPairSignalOffset];
    pairPtees[signal] = pairPtee[signal-fFillPairSignalOffset];
    pairOpAngs[signal] = pairOpAng[signal-fFillPairSignalOffset];
    pairDCAee[signal] = pairDCA[signal-fFillPairSignalOffset];

    pairCosPointAngs[signal] = pairCosPointAng[signal-fFillPairSignalOffset];
    pairDecayLengths[signal] = pairDecayLength[signal-fFillPairSignalOffset];
    pairRs[signal] = pairR[signal-fFillPairSignalOffset];


  }


  return 0;
}

//________________________________________________________________________

std::tuple<int, int, int> AliAnalysisTaskMLTreeMaker2018::GetAcceptedPairs(AliVEvent *event, Double_t gCentrality){

  // ##########################################################
  // DO PAIRING
  // ##########################################################

  std::fill(nPairs.begin(), nPairs.end(), 0);

  for (int i = 0; i < pt_tracks1.size(); i++) {

    primVerticesX[i].clear();
    primVerticesY[i].clear();
    primVerticesZ[i].clear();

    //MCvertx_tracks1[i].clear();
    //MCverty_tracks1[i].clear();
    //MCvertz_tracks1[i].clear();
    //MCvertx_tracks2[i].clear();
    //MCverty_tracks2[i].clear();
    //MCvertz_tracks2[i].clear();
 
    eta_tracks1[i].clear();
    phi_tracks1[i].clear();
    pt_tracks1[i].clear();
    charge_tracks1[i].clear();
    //MCpt_tracks1[i].clear();
    //MCeta_tracks1[i].clear();
    //MCphi_tracks1[i].clear();
    dcaXY_tracks1[i].clear();
    dcaZ_tracks1[i].clear();
    dcaXY_res_tracks1[i].clear();
    dcaZ_res_tracks1[i].clear();
    verticesX_tracks1[i].clear();
    verticesY_tracks1[i].clear();
    verticesZ_tracks1[i].clear();


    eta_tracks2[i].clear();
    phi_tracks2[i].clear();
    pt_tracks2[i].clear();
    charge_tracks2[i].clear();
    //MCpt_tracks2[i].clear();
    //MCeta_tracks2[i].clear();
    //MCphi_tracks2[i].clear();
    dcaXY_tracks2[i].clear();
    dcaZ_tracks2[i].clear();
    dcaXY_res_tracks2[i].clear();
    dcaZ_res_tracks2[i].clear();
    verticesX_tracks2[i].clear();
    verticesY_tracks2[i].clear();
    verticesZ_tracks2[i].clear();
 

    pairMasses[i].clear();
    pairPtees[i].clear();
    pairOpAngs[i].clear();

    pairCosPointAngs[i].clear();
    pairDecayLengths[i].clear();
    pairRs[i].clear();

  }

  //std::cout << "+++ Pairs +++" << std::endl;

  Int_t numberPairsULS = 0;
  Int_t numberPairsLSmm = 0;
  Int_t numberPairsLSpp = 0;
  if(doULS){
    //std::cout << "Filling ULS pairs" << std::endl;
    numberPairsULS = FillSimplePairs(fNegPart, fPosPart, event, 0);
  }
  if(doLS){
    //std::cout << "Filling LSmm pairs" << std::endl;
    numberPairsLSmm = FillSimplePairs(fNegPart, fNegPart, event, 1);
    
    //std::cout << "Filling LSmm pairs" << std::endl;
    numberPairsLSpp = FillSimplePairs(fPosPart, fPosPart, event, 2);
  }
  if(fPairMCSignal.size() > 0){
    //std::cout << "Filling MCsignal pairs" << std::endl;
    Int_t check =  FillSignalPairs(event);
  }

  return {numberPairsULS, numberPairsLSmm, numberPairsLSpp };
}

// ############################################################################
// ############################################################################
AliAnalysisTaskMLTreeMaker2018::Particle AliAnalysisTaskMLTreeMaker2018::CreateParticle(AliVParticle* mcPart1){
  double  pt1      = mcPart1->Pt();
  double  eta1     = mcPart1->Eta();
  double  phi1     = mcPart1->Phi();
  short   charge1  = mcPart1->Charge();

  float x1 = mcPart1->Xv();
  float y1 = mcPart1->Yv();
  float z1 = mcPart1->Zv();

  Particle part(pt1, eta1, phi1, charge1, x1, y1, z1);
  part.DielectronPairFromSameMother.resize(fDielectronPairNotFromSameMother.size(), false);

  return part;
}

void AliAnalysisTaskMLTreeMaker2018::CheckIfFromMotherWithDielectronAsDaughter(Particle& part){
  std::cout << "Do nothing" << std::endl;
//  if (isAOD && doULS && doLS){
//
//    for (unsigned int k = 0; k < fDielectronPairNotFromSameMother.size(); ++k){
//      if (fDielectronPairNotFromSameMother[k] == true){
//        AliAODMCParticle* mother = dynamic_cast<AliAODMCParticle*> (fMCEvent->GetTrack(part.GetMotherID()));
//        // int number_of_daugthers = mother->GetNDaughters() ;
//        int LabelFirstDaughter = mother->GetDaughterFirst();
//        int LabelLastDaughter = mother->GetDaughterLast();
//        // std::cout << "number_of_daughters = " << number_of_daugthers << "  first_daugther = " << LabelFirstDaughter << "  last_daugther = " << LabelLastDaughter << std::endl;
//
//        bool ele_from_same_mother = false;
//        bool pos_from_same_mother = false;
//        for (int daughter_i = LabelFirstDaughter; daughter_i <= LabelLastDaughter; daughter_i++){
//          int pdgCode = fMCEvent->GetTrack(daughter_i)->PdgCode();
//          if      (pdgCode == 11)  ele_from_same_mother = true;
//          else if (pdgCode == -11) pos_from_same_mother = true;
//          // std::cout << "daugther[" << daughter_i << "] with pdgcode: " << pdgCode << std::endl;
//        }
//        if (ele_from_same_mother == true && pos_from_same_mother == true) {
//          part.DielectronPairFromSameMother[k] = true;
//          // std::cout << "dielectron pair from same mother" << std::endl;
//        }
//        else{
//          part.DielectronPairFromSameMother[k] = false;
//        }
//      }
//      else{
//        part.DielectronPairFromSameMother[k] = false;
//      }
//    }
//  }
//  if (part.DielectronPairFromSameMother.size() != fDielectronPairNotFromSameMother.size()){
//    std::cout << "ERROR IN SOME PART" << std::endl;
//    // vec = std::vector<bool>(fDielectronPairNotFromSameMother.size(), false);
//  }
//  //part.DielectronPairFromSameMother = vec;

}


// ############################################################################
// ############################################################################
Bool_t AliAnalysisTaskMLTreeMaker2018::CheckIfOneIsTrue(std::vector<Bool_t>& vec){
  bool min_one_is_true = kFALSE;
  unsigned int size = vec.size();
  for (unsigned int i = 0; i < size; ++i){
    if (vec[i] == kTRUE) {min_one_is_true = kTRUE; break;}
  }
  return min_one_is_true;
}



//________________________________________________________________________

Int_t AliAnalysisTaskMLTreeMaker2018::GetAcceptedTracks(AliVEvent *event, Double_t gCentrality){

  Int_t acceptedTracks = 0;
  Int_t mpdg=0;
  eta.clear();
  phi.clear();
  pt.clear();
  NClustersITS.clear();
  NCrossedRowsTPC.clear();
  NClustersTPC.clear();
  HasSPDfirstHit.clear();
  RatioCrossedRowsFindableClusters.clear();
  NTPCSignal.clear();
  MCpt.clear();
  MCeta.clear();
  MCphi.clear();
  dcaXY.clear();
  dcaZ.clear();
  dcaXY_res.clear();
  dcaZ_res.clear();

  nITS.clear();
  nITSshared.clear();
  chi2ITS.clear();
  chi2TPC.clear();
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


  fNegPart.clear();
  fPosPart.clear();


  // Loop over tracks in event
  AliMCEvent *mcEvent=0;
  AliAODMCParticle* mcMTrack;


  // need this to use PID in dielectron framework
  varManager->SetPIDResponse(fPIDResponse);

  for (Int_t iTracks = 0; iTracks < event->GetNumberOfTracks(); iTracks++) {

    fQAHistTracks->Fill("All tracks",1);
    AliVTrack* track = dynamic_cast<AliVTrack *>(event->GetTrack(iTracks));
    if (!track) {
      AliError(Form("Could not receive track %d", iTracks));
      continue;
    }

    //std::cout << "Track issue: " << track->GetLabel() << "  " << iTracks << "  " <<  event->GetNumberOfTracks() << std::endl;

    if(!isAOD) fQAHistTracks->Fill("Not AOD track",1);
    else       fQAHistTracks->Fill("Is AOD track",1);

    fQAHistTracks->Fill("After AOD check, bef. MC",1);

    if(hasMC){

      mcEvent = MCEvent();
      if (!mcEvent) {
        AliError(Form("Could not receive MC -> hasMC set to kFALSE!!"));
	std::cout << "Could not receive MC -> hasMC set to kFALSE!!" << std::endl;

        hasMC=kFALSE;
        continue;
      }
      else{
        fQAHistTracks->Fill("After MC check",1);

//      if(CheckGenerator(TMath::Abs(track->GetLabel()))<0) continue;
      }

    }

    fQAHistTracks->Fill("Tracks aft MC Gen, bef tr cuts",1);

    UInt_t selectedMask = (1 << filter->GetCuts()->GetEntries()) - 1;
    if (selectedMask != (filter->IsSelected((AliVParticle*) track))) {
      fQAHistTracks->Fill("Tracks not selected filter",1);
      continue;
    }

    Double_t pttemp = track->Pt();
    Double_t etatemp = track->Eta();
    Int_t pdgtemp = -99999;


    if (etatemp <= fEtaMin + PRECISION)
      continue;
    if (etatemp >= fEtaMax - PRECISION)
      continue;
    if (pttemp <= fPtMin + PRECISION)
      continue;
    if (pttemp >= fPtMax - PRECISION)
      continue;

    fQAHistTracks->Fill("After track cuts",1);

    if(hasMC){

      AliAODMCParticle* mcTrack = dynamic_cast<AliAODMCParticle *>(mcEvent->GetTrack(TMath::Abs(track->GetLabel())));

      pdgtemp = mcTrack->PdgCode();
      if (fabs(pdgtemp) != 11)
        continue;
      pdg.push_back( pdgtemp );


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

      if( abs(mcTrack->PdgCode())==11 && mpdg!=22) fQAHistTracks->Fill("Selected electrons, non-conversion",1);

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

        if(tempLabelFirstMother != -1) {          // if grandmother not primary!
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
      iGenIndex.push_back(CheckGenerator(TMath::Abs(track->GetLabel()), fGeneratorHashs,kTRUE));
      iPdgFirstMother.push_back(tempPdgFirstMother);

      } //End if hasMC


      if(!acceptedTracks){   //Get vertex only for first track in event
        Double_t vert[3] = {0};
        event->GetPrimaryVertex()->GetXYZ(vert);
        vertx= vert[0];
        verty= vert[1];
        vertz= vert[2];
      }


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


      Double_t d0z0[2]={-999.0,-999.0};
      Double_t dcaRes[3] = {-999.,-999.,-999.};
      //Get DCA position
      if(isAOD){
        
        GetDCA(event,(AliAODTrack*)track,d0z0,dcaRes);
        dcaXY.push_back((Float_t)d0z0[0]);
        dcaZ.push_back((Float_t)d0z0[1]);

        dcaXY_res.push_back((Float_t)dcaRes[0]);
        dcaZ_res.push_back((Float_t)dcaRes[2]);
      }
      else{
	Float_t d0z0_tmp[2]={-999.0,-999.0};
        track->GetImpactParameters( &d0z0_tmp[0], &d0z0_tmp[1]); //GetImpactParameter is also used in AliESDtrackCuts.cxx to cut on DCA to vertex

        d0z0[0] = (Double_t)d0z0_tmp[0];
	d0z0[1] = (Double_t)d0z0_tmp[1];

        dcaXY.push_back((Float_t)d0z0[0]);
        dcaZ.push_back((Float_t)d0z0[1]);
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
      chi2TPC.push_back(track->GetTPCchi2());//this variable will be always 0 for AODs (not yet in)

      if(isAOD){ chi2GlobalPerNDF.push_back(((AliAODTrack*)track)->Chi2perNDF());
                 chi2GlobalvsTPC.push_back(((AliAODTrack*)track)->GetChi2TPCConstrainedVsGlobal());
      }


       // ##########################################################

      int abslabel = TMath::Abs(track->GetLabel());

      // ##########################################################
      // Apply MC signals
      std::vector<Bool_t> mcSignal_acc(fSingleLegMCSignal.size(), kFALSE); // vector which stores if track is accepted by [i]-th mcsignal
      CheckSingleLegMCsignals(mcSignal_acc, abslabel);  // needs to be tested JEROME

      // ##########################################################
      // check if at least one mc signal is true otherwise skip this particle
      if (doSingleLegMCSignal && CheckIfOneIsTrue(mcSignal_acc) == kFALSE) continue;

      // ##########################################################
      // check if correct generator used
      bool generatorForMCSignal  = CheckGenerator(iTracks, fGeneratorMCSignalHashs, kTRUE);
      bool generatorForULSSignal = CheckGenerator(iTracks, fGeneratorULSSignalHashs, kTRUE);
 
      

      if (!generatorForMCSignal && !generatorForULSSignal) continue;



      // Create summary particle from track info 
      int motherID = TMath::Abs(fMCEvent->GetTrack(abslabel)->GetMother());
      Particle part  = CreateParticle(track);
      part.isMCSignal = mcSignal_acc;
//      part.SetTrackID(fabs(track->GetLabel())); //iTracks
      part.SetTrackID(iTracks); //iTracks
      part.SetTrackLabel(fabs(track->GetLabel())); //iTracks

      part.SetMotherID(motherID);
      part.SetULSSignalPair(generatorForULSSignal);
      part.SetMCSignalPair(generatorForMCSignal);


      part.SetDCA(d0z0[0], d0z0[1]);
      part.SetDCAres(dcaRes[0], dcaRes[2]);
      // ##########################################################
      // check if electron comes from a mother with ele+pos as daughters
//      CheckIfFromMotherWithDielectronAsDaughter(part);


      if(doPairing == true && part.fCharge < 0) fNegPart.push_back(part);
      if(doPairing == true && part.fCharge > 0) fPosPart.push_back(part);

      acceptedTracks++;
  }

  return acceptedTracks;
}

//______________________________________________________________________________
Bool_t AliAnalysisTaskMLTreeMaker2018::GetDCA(const AliVEvent* event, const AliAODTrack* track, Double_t* d0z0, Double_t* covd0z0)
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


void AliAnalysisTaskMLTreeMaker2018::SetupTrackCuts(AliDielectronCutGroup* f)
{
  filter   = new AliAnalysisFilter("filter","filter");
  filter->AddCuts(f);
}


void AliAnalysisTaskMLTreeMaker2018::SetupEventCuts(AliDielectronEventCuts* f)
{
  evfilter   = new AliAnalysisFilter("evfilter","evfilter");
  evfilter->AddCuts(f);
}



bool AliAnalysisTaskMLTreeMaker2018::CheckGenerator(int trackID, std::vector<unsigned int> vecHashes, Bool_t isGen){
  //std::cout << "Entering Generator" << std::endl;
  //if(fRejectParticleFromOOB){
  //  if(isAOD){//for AOD
  //    AliAODMCHeader* mcHeader = (AliAODMCHeader*)fEvent->GetList()->FindObject(AliAODMCHeader::StdBranchName());
  //    if (!mcHeader) {
  //        AliError("Could not find MC headr in AOD");
  //        return false;
  //    }
  //    TClonesArray* mcArray = dynamic_cast<TClonesArray*>(fEvent->FindListObject(AliAODMCParticle::StdBranchName()));
  //    if (!mcArray) {
  //        AliError("Could not find MC array in AOD");
  //        return false;
  //    }
  //    if(isGen && AliAnalysisUtils::IsParticleFromOutOfBunchPileupCollision(TMath::Abs(trackID), mcHeader, mcArray)) return false;//particles from pileup collision should NOT be used.
  //  }
  //  else{//for ESD
  //      if(isGen && AliAnalysisUtils::IsParticleFromOutOfBunchPileupCollision(TMath::Abs(trackID), fMC)) return false;//particles from pileup collision should NOT be used.
  //  }
  //}

  if (vecHashes.size() == 0){
    return true;
  }

  TString genname="";
  Bool_t hasGenerator = fMCEvent->GetCocktailGenerator(TMath::Abs(trackID), genname);
  //AliMCParticle* p = (AliMCParticle*)fMC->GetTrack(TMath::Abs(trackID));
  //Int_t genID = p->GetGeneratorIndex();
  //AliInfo(Form("genID = %d , generator name = %s",genID,genname.Data()));

  if(!hasGenerator) {
    Printf("no cocktail header list was found for this track");
    return false;
  }
  else{
    for (unsigned int i = 0; i < vecHashes.size(); ++i){
      if (genname.Hash() == vecHashes[i]) return true;
    }//end of vecHashes loop
    return false;
  }

  return false; // should not happen
}




// ############################################################################
// ############################################################################
void AliAnalysisTaskMLTreeMaker2018::CheckSingleLegMCsignals(std::vector<Bool_t>& vec, const int tracklabel){
  for (unsigned int i = 0; i < fSingleLegMCSignal.size(); ++i){
    vec.at(i) = AliDielectronMC::Instance()->IsMCTruth(tracklabel, &(fSingleLegMCSignal[i]), 1);
  }
}

