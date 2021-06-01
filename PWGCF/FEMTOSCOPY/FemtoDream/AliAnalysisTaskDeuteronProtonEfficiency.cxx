#include "TChain.h" 
#include "TTree.h"
#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "THnSparse.h"
#include "TCanvas.h"
#include "AliESDtrackCuts.h"
#include "AliESDpid.h"
#include "AliESDVertex.h"
#include "AliVertexerTracks.h"
#include "TFile.h"
#include "AliAnalysisUtils.h"
#include "AliEventCuts.h"

#include "AliMCParticle.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliCentrality.h"
#include "TDatabasePDG.h"
#include "AliMultSelection.h"

#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliESDv0.h"
#include <TRandom3.h>

#include "AliKFParticle.h"
#include "AliKFVertex.h"
#include "AliKFParticleBase.h"
#include "Riostream.h"
#include <iostream>
#include <fstream>
#include <string>

#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliStack.h"
#include "AliGenHijingEventHeader.h"
#include "AliPID.h"
#include "TPDGCode.h"
#include "TDatabasePDG.h"
#include "AliPIDResponse.h"




#include "AliAnalysisTaskDeuteronProtonEfficiency.h"


using namespace std;
ClassImp(AliAnalysisTaskDeuteronProtonEfficiency)



AliAnalysisTaskDeuteronProtonEfficiency::AliAnalysisTaskDeuteronProtonEfficiency() : AliAnalysisTaskSE(),
  fESD(0),
  fESDHandler(0),
  fEvent(0x0),
  mcEvent(0),
  fMCtrue(0),
  fESDtrackCutsProton(0),
  fESDtrackCutsDeuteron(0),
  fHistList(0),
  fHistPtProtonGen(0),
  fHistEtaProtonGen(0),
  fHistPtAntiProtonGen(0),
  fHistEtaAntiProtonGen(0),
  fHistPtDeuteronGen(0),
  fHistEtaDeuteronGen(0),
  fHistPtAntiDeuteronGen(0),
  fHistEtaAntiDeuteronGen(0),
  fHistPtHelium3Gen(0),
  fHistEtaHelium3Gen(0),
  fHistPtAntiHelium3Gen(0),
  fHistEtaAntiHelium3Gen(0),
  fHistSEDPairGen(0),
  fHistSEDAntiPairGen(0),
  fHistPtParticlesGen(0),
  fHistPtAntiParticlesGen(0),
  fHistPtProtonRec(0),
  fHistEtaProtonRec(0),
  fHistPtAntiProtonRec(0),
  fHistEtaAntiProtonRec(0),
  fHistPtDeuteronRec(0),
  fHistEtaDeuteronRec(0),
  fHistPtAntiDeuteronRec(0),
  fHistEtaAntiDeuteronRec(0),
  fHistPtHelium3Rec(0),
  fHistEtaHelium3Rec(0),
  fHistPtAntiHelium3Rec(0),
  fHistEtaAntiHelium3Rec(0),
  fHistSEDPairRec(0),
  fHistSEDAntiPairRec(0),
  fHistPtParticlesRec(0),
  fHistPtAntiParticlesRec(0),
  fHistEventCounter(0),
  fPIDResponse(0)
{


}



AliAnalysisTaskDeuteronProtonEfficiency::AliAnalysisTaskDeuteronProtonEfficiency(const char *name) : AliAnalysisTaskSE(name),
  fESD(0),
  fESDHandler(0),
  fEvent(0x0),
  mcEvent(0),
  fMCtrue(0),
  fESDtrackCutsProton(0),
  fESDtrackCutsDeuteron(0),
  fHistList(0),
  fHistPtProtonGen(0),
  fHistEtaProtonGen(0),
  fHistPtAntiProtonGen(0),
  fHistEtaAntiProtonGen(0),
  fHistPtDeuteronGen(0),
  fHistEtaDeuteronGen(0),
  fHistPtAntiDeuteronGen(0),
  fHistEtaAntiDeuteronGen(0),
  fHistPtHelium3Gen(0),
  fHistEtaHelium3Gen(0),
  fHistPtAntiHelium3Gen(0),
  fHistEtaAntiHelium3Gen(0),
  fHistSEDPairGen(0),
  fHistSEDAntiPairGen(0),
  fHistPtParticlesGen(0),
  fHistPtAntiParticlesGen(0),
  fHistPtProtonRec(0),
  fHistEtaProtonRec(0),
  fHistPtAntiProtonRec(0),
  fHistEtaAntiProtonRec(0),
  fHistPtDeuteronRec(0),
  fHistEtaDeuteronRec(0),
  fHistPtAntiDeuteronRec(0),
  fHistEtaAntiDeuteronRec(0),
  fHistPtHelium3Rec(0),
  fHistEtaHelium3Rec(0),
  fHistPtAntiHelium3Rec(0),
  fHistEtaAntiHelium3Rec(0),
  fHistSEDPairRec(0),
  fHistSEDAntiPairRec(0),
  fHistPtParticlesRec(0),
  fHistPtAntiParticlesRec(0),
  fHistEventCounter(0),
  fPIDResponse(0)
{

  DefineInput(0,TChain::Class());
  DefineOutput(1,TList::Class());


  // Proton cuts
  fESDtrackCutsProton = new AliESDtrackCuts("AliESDtrackCutsProton","AliESDtrackCutsProton");
  fESDtrackCutsProton = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts(); // FilterBit 128
  fESDtrackCutsProton->SetPtRange(0.4,4.0);
  fESDtrackCutsProton->SetEtaRange(-0.8,0.8);
  fESDtrackCutsProton->SetMinNClustersTPC(80);
  fESDtrackCutsProton->SetMinNCrossedRowsTPC(70);
  fESDtrackCutsProton->SetMinRatioCrossedRowsOverFindableClustersTPC(0.83);
  fESDtrackCutsProton->SetAcceptSharedTPCClusters(false);
  fESDtrackCutsProton->SetMaxDCAToVertexXY(0.1);
  fESDtrackCutsProton->SetMaxDCAToVertexZ(0.2);

  // Deuteron cuts
  fESDtrackCutsDeuteron = new AliESDtrackCuts("AliESDtrackCutsDeuteron","AliESDtrackCutsDeuteron");
  fESDtrackCutsDeuteron->SetEtaRange(-0.8,0.8);
  fESDtrackCutsDeuteron->SetPtRange(0.4,4.0);
  fESDtrackCutsDeuteron->SetMinNClustersTPC(80);
  fESDtrackCutsDeuteron->SetMinNCrossedRowsTPC(70);
  fESDtrackCutsDeuteron->SetMinRatioCrossedRowsOverFindableClustersTPC(0.83);
  fESDtrackCutsDeuteron->SetAcceptSharedTPCClusters(false);
  fESDtrackCutsDeuteron->SetMaxDCAToVertexXY(0.1);
  fESDtrackCutsDeuteron->SetMaxDCAToVertexZ(0.2);
  fESDtrackCutsDeuteron->SetRequireITSPid(true); // FilterBit 256
  fESDtrackCutsDeuteron->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);	// FilterBit 256


}



// destructor
AliAnalysisTaskDeuteronProtonEfficiency::~AliAnalysisTaskDeuteronProtonEfficiency(){

  if(fHistList)
    {

      fHistList->Clear();
      delete fHistList;

    }

}






void AliAnalysisTaskDeuteronProtonEfficiency::UserCreateOutputObjects()
{

  fHistList = new TList();
  fHistList->SetOwner();

  // histograms for generated particles
  fHistPtProtonGen = new TH1F("fHistPtProtonGen","#it{p}_{T} distribution of generated Protons",240,0.0,6.0);
  fHistPtProtonGen->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  fHistPtProtonGen->GetYaxis()->SetTitle("entries");

  fHistEtaProtonGen = new TH1F("fHistEtaProtonGen","#eta distribution of generated Protons",200,-1.5,1.5);
  fHistEtaProtonGen->GetXaxis()->SetTitle("#eta");
  fHistEtaProtonGen->GetYaxis()->SetTitle("entries");

  fHistPtAntiProtonGen = new TH1F("fHistPtAntiProtonGen","#it{p}_{T} distribution of generated Antiprotons",240,0.0,6.0);
  fHistPtAntiProtonGen->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  fHistPtAntiProtonGen->GetYaxis()->SetTitle("entries");

  fHistEtaAntiProtonGen = new TH1F("fHistEtaAntiProtonGen","#eta distribution of generated Antiprotons",200,-1.5,1.5);
  fHistEtaAntiProtonGen->GetXaxis()->SetTitle("#eta");
  fHistEtaAntiProtonGen->GetYaxis()->SetTitle("entries");

  fHistPtDeuteronGen = new TH1F("fHistPtDeuteronGen","#it{p}_{T} distribution of generated Deuterons",240,0.0,6.0);
  fHistPtDeuteronGen->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  fHistPtDeuteronGen->GetYaxis()->SetTitle("entries");

  fHistEtaDeuteronGen = new TH1F("fHistEtaDeuteronGen","#eta distribution of generated Deuterons",200,-1.5,1.5);
  fHistEtaDeuteronGen->GetXaxis()->SetTitle("#eta");
  fHistEtaDeuteronGen->GetYaxis()->SetTitle("entries");

  fHistPtAntiDeuteronGen = new TH1F("fHistPtAntiDeuteronGen","#it{p}_{T} distribution of generated Antideuterons",240,0.0,6.0);
  fHistPtAntiDeuteronGen->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  fHistPtAntiDeuteronGen->GetYaxis()->SetTitle("entries");

  fHistEtaAntiDeuteronGen = new TH1F("fHistEtaAntiDeuteronGen","#eta distribution of generated Antideuterons",200,-1.5,1.5);
  fHistEtaAntiDeuteronGen->GetXaxis()->SetTitle("#eta");
  fHistEtaAntiDeuteronGen->GetYaxis()->SetTitle("entries");

  fHistPtHelium3Gen = new TH1F("fHistPtHelium3Gen","#it{p}_{T} distribution of the generated Pair",240,0.0,6.0);
  fHistPtHelium3Gen->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  fHistPtHelium3Gen->GetYaxis()->SetTitle("entries");

  fHistEtaHelium3Gen = new TH1F("fHistEtaHelium3Gen","#eta distribution of the generated Pair",200,-1.5,1.5);
  fHistEtaHelium3Gen->GetXaxis()->SetTitle("#eta");
  fHistEtaHelium3Gen->GetYaxis()->SetTitle("entries");

  fHistPtAntiHelium3Gen = new TH1F("fHistPtAntiHelium3Gen","#it{p}_{T} distribution of the generated Antipair",240,0.0,6.0);
  fHistPtAntiHelium3Gen->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  fHistPtAntiHelium3Gen->GetYaxis()->SetTitle("entries");

  fHistEtaAntiHelium3Gen = new TH1F("fHistEtaAntiHelium3Gen","#eta distribution of the generated Antipair",200,-1.5,1.5);
  fHistEtaAntiHelium3Gen->GetXaxis()->SetTitle("#eta");
  fHistEtaAntiHelium3Gen->GetYaxis()->SetTitle("entries");

  fHistSEDPairGen = new TH1F("fHistSEDPairGen","same-event distribution of generated d-p pairs",750,0.0,3.0);
  fHistSEDPairGen->GetXaxis()->SetTitle("#it{k*} (GeV/#it{c})");
  fHistSEDPairGen->GetYaxis()->SetTitle("entries");

  fHistSEDAntiPairGen = new TH1F("fHistSEDAntiPairGen","same-event distribution of generated #bar{d}-#bar{p} pairs",750,0.0,3.0);
  fHistSEDAntiPairGen->GetXaxis()->SetTitle("#it{k*} (GeV/#it{c})");
  fHistSEDAntiPairGen->GetYaxis()->SetTitle("entries");

  fHistPtParticlesGen = new TH2F("fHistPtParticlesGen","#it{p}_{T} distribution of generated protons and deuterons",240,0.0,6.0,240,0.0,6.0);
  fHistPtParticlesGen->GetXaxis()->SetTitle("#it{p}_{T} of deuterons (GeV/#it{c})");
  fHistPtParticlesGen->GetYaxis()->SetTitle("#it{p}_{T} of protons (GeV/#it{c})");

  fHistPtAntiParticlesGen = new TH2F("fHistPtAntiParticlesGen","#it{p}_{T} distribution of generated antiprotons and antideuterons",240,0.0,6.0,240,0.0,6.0);
  fHistPtAntiParticlesGen->GetXaxis()->SetTitle("#it{p}_{T} of antideuterons (GeV/#it{c})");
  fHistPtAntiParticlesGen->GetYaxis()->SetTitle("#it{p}_{T} of antiprotons (GeV/#it{c})");

  //  histograms for reconstructed particles
  fHistPtProtonRec = new TH1F("fHistPtProtonRec","#it{p}_{T} distribution of reconstructed Protons",240,0.0,6.0);
  fHistPtProtonRec->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  fHistPtProtonRec->GetYaxis()->SetTitle("entries");

  fHistEtaProtonRec = new TH1F("fHistEtaProtonRec","#eta distribution of reconstructed Protons",200,-1.5,1.5);
  fHistEtaProtonRec->GetXaxis()->SetTitle("#eta");
  fHistEtaProtonRec->GetYaxis()->SetTitle("entries");

  fHistPtAntiProtonRec = new TH1F("fHistPtAntiProtonRec","#it{p}_{T} distribution of reconstructed Antiprotons",240,0.0,6.0);
  fHistPtAntiProtonRec->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  fHistPtAntiProtonRec->GetYaxis()->SetTitle("entries");

  fHistEtaAntiProtonRec = new TH1F("fHistEtaAntiProtonRec","#eta distribution of reconstruced Antiprotons",200,-1.5,1.5);
  fHistEtaAntiProtonRec->GetXaxis()->SetTitle("#eta");
  fHistEtaAntiProtonRec->GetYaxis()->SetTitle("entries");

  fHistPtDeuteronRec = new TH1F("fHistPtDeuteronRec","#it{p}_{T} distribution of reconstructed Deuterons",240,0.0,6.0);
  fHistPtDeuteronRec->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  fHistPtDeuteronRec->GetYaxis()->SetTitle("entries");

  fHistEtaDeuteronRec = new TH1F("fHistEtaDeuteronRec","#eta distribution of reconstructed Deuterons",200,-1.5,1.5);
  fHistEtaDeuteronRec->GetXaxis()->SetTitle("#eta");
  fHistEtaDeuteronRec->GetYaxis()->SetTitle("entries");

  fHistPtAntiDeuteronRec = new TH1F("fHistPtAntiDeuteronRec","#it{p}_{T} distribution of reconstructed Antideuterons",240,0.0,6.0);
  fHistPtAntiDeuteronRec->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  fHistPtAntiDeuteronRec->GetYaxis()->SetTitle("entries");

  fHistEtaAntiDeuteronRec = new TH1F("fHistEtaAntiDeuteronRec","#eta distribution of reconstructed Antideuterons",200,-1.5,1.5);
  fHistEtaAntiDeuteronRec->GetXaxis()->SetTitle("#eta");
  fHistEtaAntiDeuteronRec->GetYaxis()->SetTitle("entries");

  fHistPtHelium3Rec = new TH1F("fHistPtHelium3Rec","#it{p}_{T} distribution of the reconstructed Pair",240,0.0,6.0);
  fHistPtHelium3Rec->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  fHistPtHelium3Rec->GetYaxis()->SetTitle("entries");

  fHistEtaHelium3Rec = new TH1F("fHistEtaHelium3Rec","#eta distribution of the reconstructed Pair",200,-1.5,1.5);
  fHistEtaHelium3Rec->GetXaxis()->SetTitle("#eta");
  fHistEtaHelium3Rec->GetYaxis()->SetTitle("entries");

  fHistPtAntiHelium3Rec = new TH1F("fHistPtAntiHelium3Rec","#it{p}_{T} distribution of the reconstructed Antipair",240,0.0,6.0);
  fHistPtAntiHelium3Rec->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  fHistPtAntiHelium3Rec->GetYaxis()->SetTitle("entries");

  fHistEtaAntiHelium3Rec = new TH1F("fHistEtaAntiHelium3Rec","#eta distribution of the reconstructed Antipair",200,-1.5,1.5);
  fHistEtaAntiHelium3Rec->GetXaxis()->SetTitle("#eta");
  fHistEtaAntiHelium3Rec->GetYaxis()->SetTitle("entries");

  fHistSEDPairRec = new TH1F("fHistSEDPairRec","same-event distribution of reconstructed d-p pairs",750,0.0,3.0);
  fHistSEDPairRec->GetXaxis()->SetTitle("#it{k*} (GeV/#it{c})");
  fHistSEDPairRec->GetYaxis()->SetTitle("entries");

  fHistSEDAntiPairRec = new TH1F("fHistSEDAntiPairRec","same-event distribution of reconstructed #bar{d}-#bar{p} pairs",750,0.0,3.0);
  fHistSEDAntiPairRec->GetXaxis()->SetTitle("#it{k*} (GeV/#it{c})");
  fHistSEDAntiPairRec->GetYaxis()->SetTitle("entries");

  fHistPtParticlesRec = new TH2F("fHistPtParticlesRec","#it{p}_{T} distribution of reconstructed protons and deuterons",240,0.0,6.0,240,0.0,6.0);
  fHistPtParticlesRec->GetXaxis()->SetTitle("#it{p}_{T} of deuterons (GeV/#it{c})");
  fHistPtParticlesRec->GetYaxis()->SetTitle("#it{p}_{T} of protons (GeV/#it{c})");

  fHistPtAntiParticlesRec = new TH2F("fHistPtAntiParticlesRec","#it{p}_{T} distribution of reconstructed antiprotons and antideuterons",240,0.0,6.0,240,0.0,6.0);
  fHistPtAntiParticlesRec->GetXaxis()->SetTitle("#it{p}_{T} of antideuterons (GeV/#it{c})");
  fHistPtAntiParticlesRec->GetYaxis()->SetTitle("#it{p}_{T} of antiprotons (GeV/#it{c})");

  // histogram to count the number of events
  fHistEventCounter = new TH1F("fHistEventCounter","simple event counter",1,0,1);
  fHistEventCounter->GetXaxis()->SetBinLabel(1,"number of events");
  fHistEventCounter->GetYaxis()->SetTitle("entries");

  fHistList->Add(fHistPtProtonGen);
  fHistList->Add(fHistEtaProtonGen);
  fHistList->Add(fHistPtAntiProtonGen);
  fHistList->Add(fHistEtaAntiProtonGen);
  fHistList->Add(fHistPtDeuteronGen);
  fHistList->Add(fHistEtaDeuteronGen);
  fHistList->Add(fHistPtAntiDeuteronGen);
  fHistList->Add(fHistEtaAntiDeuteronGen);
  fHistList->Add(fHistPtHelium3Gen);
  fHistList->Add(fHistEtaHelium3Gen);
  fHistList->Add(fHistPtAntiHelium3Gen);
  fHistList->Add(fHistEtaAntiHelium3Gen);
  fHistList->Add(fHistSEDPairGen);
  fHistList->Add(fHistSEDAntiPairGen);
  fHistList->Add(fHistPtParticlesGen);
  fHistList->Add(fHistPtAntiParticlesGen);

  fHistList->Add(fHistPtProtonRec);
  fHistList->Add(fHistEtaProtonRec);
  fHistList->Add(fHistPtAntiProtonRec);
  fHistList->Add(fHistEtaAntiProtonRec);
  fHistList->Add(fHistPtDeuteronRec);
  fHistList->Add(fHistEtaDeuteronRec);
  fHistList->Add(fHistPtAntiDeuteronRec);
  fHistList->Add(fHistEtaAntiDeuteronRec);
  fHistList->Add(fHistPtHelium3Rec);
  fHistList->Add(fHistEtaHelium3Rec);
  fHistList->Add(fHistPtAntiHelium3Rec);
  fHistList->Add(fHistEtaAntiHelium3Rec);
  fHistList->Add(fHistSEDPairRec);
  fHistList->Add(fHistSEDAntiPairRec);
  fHistList->Add(fHistPtParticlesRec);
  fHistList->Add(fHistPtAntiParticlesRec);

  fHistList->Add(fHistEventCounter);

//  AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
//  AliInputEventHandler *inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());

  PostData(1,fHistList);

} // end of UserCreateOutputObjects 



void AliAnalysisTaskDeuteronProtonEfficiency::UserExec(Option_t *)
{

  fMCtrue = kTRUE;
  AliMCEventHandler *eventHandler = dynamic_cast<AliMCEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());

  if(!eventHandler)
    {
    
      fMCtrue = kFALSE;
      Error("AliAnalysisTaskDeuteronProtonEfficiency.cxx","No MC Handler found.");

    }


  mcEvent = 0x0;
  AliStack* stack = 0x0;

  if (eventHandler) mcEvent = eventHandler->MCEvent();
  if (!mcEvent) if (fMCtrue) return;

  if (fMCtrue)
    {
      stack = mcEvent->Stack();
      if (!stack) return;
    }


  fESDHandler= dynamic_cast<AliESDInputHandler*>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
  if (!fESDHandler) {
    AliError("Could not get ESD InputHandler");
    return;
  } 
  fESD=dynamic_cast<AliESDEvent*>(InputEvent());
  if (!fESD) return;


  fPIDResponse = dynamic_cast<AliPIDResponse*>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()->GetPIDResponse());

  // Check the PID response
  if(!fPIDResponse) {
    AliError("Cannot get pid response");
    return;
  }


  // Masses (GeV/c2)
  double ProtonMass   = 0.93827;
  double DeuteronMass = 1.87561;

  // PDG codes
  int ProtonPDG	      =	       2212;
  int DeuteronPDG     =  1000010020;
/*  int PionPDG	      =		211;
  int KaonPDG	      =		321;
  int TritonPDG	      =	 1000010030;
  int Helium3PDG      =	 1000020030;
  int HypertritonPDG  =	 1010010030;
  int ElectronPDG     =		 11;
  int MyonPDG	      =		 13;
*/

  // Single-particle and pair cuts
  double EtaLimit1 = -0.8;
  double EtaLimit2 = +0.8;
  double PairPtLimit1 = 1.0;
  double PairPtLimit2 = 2.0;
  double ProtonPtLimit1 = 0.0;
  double DeuteronPtLimit1 = 0.0;
  double ProtonPtThreshold = 1.4;
  double DeuteronPtThreshold = 0.7;

  // Variables for generated particles
  double MomentumProtonGen[3]    = {0,0,0};
  double MomentumDeuteronGen[3]  = {0,0,0};
  double EtaProtonGen    = 0.0;
  double EtaDeuteronGen  = 0.0;
  double EtaHelium3Gen   = 0.0;
  double PtProtonGen	 = 0.0;
  double PtDeuteronGen	 = 0.0;
  double PtHelium3Gen	 = 0.0;

  // Variables for reconstructed particles
  double MomentumProtonRec[3]    = {0,0,0};
  double MomentumDeuteronRec[3]  = {0,0,0};
  double EtaProtonRec    = 0.0;
  double EtaDeuteronRec  = 0.0;
  double EtaHelium3Rec   = 0.0;
  double PtProtonRec	 = 0.0;
  double PtDeuteronRec	 = 0.0;
  double PtHelium3Rec	 = 0.0;

  double RelativeMomentum = 0.0;

  int nStackLoop1 = 0;
  int nStackLoop2 = 0;

  // Variables for PID
  double nSigmaTPCTOFproton = 0.0;
  double nSigmaTPCTOFdeuteron = 0.0;
  double nSigmaTPCproton = 0.0;
  double nSigmaTPCdeuteron = 0.0;
  double nSigmaTOFproton = 0.0;
  double nSigmaTOFdeuteron = 0.0;

  bool TPC1_OK = false;
  bool TPC2_OK = false;
  bool TOF1_OK = false;
  bool TOF2_OK = false;


  if(fMCtrue){

    // loop for the generated particles
    for(nStackLoop1 = 0; nStackLoop1 < stack->GetNtrack(); nStackLoop1++)
      {

	AliMCParticle *Particle1 = (AliMCParticle*)(mcEvent->GetTrack(nStackLoop1));
	if(!Particle1->IsPhysicalPrimary()) continue;

	if(Particle1->PdgCode() == ProtonPDG)
	  {

	    for(nStackLoop2 = 0; nStackLoop2 < stack->GetNtrack(); nStackLoop2++)
	      {

		AliMCParticle *Particle2 = (AliMCParticle*)(mcEvent->GetTrack(nStackLoop2));
		if(!Particle2->IsPhysicalPrimary()) continue;

		if(Particle2->PdgCode() == DeuteronPDG)
		  {

		    MomentumProtonGen[0] = Particle1->Px();
		    MomentumProtonGen[1] = Particle1->Py();
		    MomentumProtonGen[2] = Particle1->Pz();

		    MomentumDeuteronGen[0] = Particle2->Px();
		    MomentumDeuteronGen[1] = Particle2->Py();
		    MomentumDeuteronGen[2] = Particle2->Pz();
      
		    TLorentzVector LorentzVectorProtonGen;
		    TLorentzVector LorentzVectorDeuteronGen;
		    TLorentzVector LorentzVectorHelium3Gen;

		    LorentzVectorProtonGen.SetXYZM(MomentumProtonGen[0],MomentumProtonGen[1],MomentumProtonGen[2],ProtonMass);
		    LorentzVectorDeuteronGen.SetXYZM(MomentumDeuteronGen[0],MomentumDeuteronGen[1],MomentumDeuteronGen[2],DeuteronMass);
		    LorentzVectorHelium3Gen = LorentzVectorProtonGen + LorentzVectorDeuteronGen;

		    PtProtonGen = LorentzVectorProtonGen.Pt();
		    PtDeuteronGen = LorentzVectorDeuteronGen.Pt();
		    PtHelium3Gen = LorentzVectorHelium3Gen.Pt();

		    EtaProtonGen = LorentzVectorProtonGen.Eta();
		    EtaDeuteronGen = LorentzVectorDeuteronGen.Eta();
		    EtaHelium3Gen = LorentzVectorHelium3Gen.Eta();

		    // single particle cuts and pair cuts
		    if(PtDeuteronGen  < DeuteronPtLimit1)				    continue;
		    if(PtProtonGen < ProtonPtLimit1)					    continue;
		    if((EtaProtonGen < EtaLimit1)	  || (EtaProtonGen > EtaLimit2))    continue;
		    if((EtaDeuteronGen < EtaLimit1)	  || (EtaDeuteronGen > EtaLimit2))  continue;
		    if((PtHelium3Gen < PairPtLimit1)	  || (PtHelium3Gen > PairPtLimit2)) continue;

		    RelativeMomentum = CalculateRelativeMomentum(LorentzVectorHelium3Gen,LorentzVectorProtonGen,LorentzVectorDeuteronGen);

		    fHistPtProtonGen->Fill(PtProtonGen);
		    fHistEtaProtonGen->Fill(EtaProtonGen);
		    fHistPtDeuteronGen->Fill(PtDeuteronGen);
		    fHistEtaDeuteronGen->Fill(EtaDeuteronGen);
		    fHistPtHelium3Gen->Fill(PtHelium3Gen);
		    fHistEtaHelium3Gen->Fill(EtaHelium3Gen);
		    fHistSEDPairGen->Fill(RelativeMomentum);
		    fHistPtParticlesGen->Fill(PtDeuteronGen,PtProtonGen);

		  } // end of "if deuteron PDG"

	      } // end of stack loop 2

	  }// end of "if proton PDG"

      } // end of stack loop 1



    // loop for the generated antiparticles
    for(nStackLoop1 = 0; nStackLoop1 < stack->GetNtrack(); nStackLoop1++)
      {

	AliMCParticle *Particle1 = (AliMCParticle*)(mcEvent->GetTrack(nStackLoop1));
	if(!Particle1->IsPhysicalPrimary()) continue;
      
	if(Particle1->PdgCode() == -ProtonPDG)
	  {

	    for(nStackLoop2 = 0; nStackLoop2 < stack->GetNtrack(); nStackLoop2++)
	      {

		AliMCParticle *Particle2 = (AliMCParticle*)(mcEvent->GetTrack(nStackLoop2));
		if(!Particle2->IsPhysicalPrimary()) continue;

		if(Particle2->PdgCode() == -DeuteronPDG)
		  {

		    MomentumProtonGen[0] = Particle1->Px();
		    MomentumProtonGen[1] = Particle1->Py();
		    MomentumProtonGen[2] = Particle1->Pz();

		    MomentumDeuteronGen[0] = Particle2->Px();
		    MomentumDeuteronGen[1] = Particle2->Py();
		    MomentumDeuteronGen[2] = Particle2->Pz();

		    TLorentzVector LorentzVectorProtonGen;
		    TLorentzVector LorentzVectorDeuteronGen;
		    TLorentzVector LorentzVectorHelium3Gen;

		    LorentzVectorProtonGen.SetXYZM(MomentumProtonGen[0],MomentumProtonGen[1],MomentumProtonGen[2],ProtonMass);
		    LorentzVectorDeuteronGen.SetXYZM(MomentumDeuteronGen[0],MomentumDeuteronGen[1],MomentumDeuteronGen[2],DeuteronMass);

		    LorentzVectorHelium3Gen = LorentzVectorProtonGen + LorentzVectorDeuteronGen;

		    EtaProtonGen = LorentzVectorProtonGen.Eta();
		    EtaDeuteronGen = LorentzVectorDeuteronGen.Eta();
		    EtaHelium3Gen = LorentzVectorHelium3Gen.Eta();

		    PtProtonGen = LorentzVectorProtonGen.Pt();
		    PtDeuteronGen = LorentzVectorDeuteronGen.Pt();
		    PtHelium3Gen = LorentzVectorHelium3Gen.Pt();

		    // single anti-particle cuts and anti-pair cuts
		    if(PtDeuteronGen  < DeuteronPtLimit1)				    continue;
		    if(PtProtonGen < ProtonPtLimit1)					    continue;
		    if((EtaProtonGen < EtaLimit1)	  || (EtaProtonGen > EtaLimit2))    continue;
		    if((EtaDeuteronGen < EtaLimit1)	  || (EtaDeuteronGen > EtaLimit2))  continue;
		    if((PtHelium3Gen < PairPtLimit1)	  || (PtHelium3Gen > PairPtLimit2)) continue;

		    RelativeMomentum = CalculateRelativeMomentum(LorentzVectorHelium3Gen,LorentzVectorProtonGen,LorentzVectorDeuteronGen);

		    fHistPtAntiProtonGen->Fill(PtProtonGen);
		    fHistEtaAntiProtonGen->Fill(EtaProtonGen);
		    fHistPtAntiDeuteronGen->Fill(PtDeuteronGen);
		    fHistEtaAntiDeuteronGen->Fill(EtaDeuteronGen);
		    fHistPtAntiHelium3Gen->Fill(PtHelium3Gen);
		    fHistEtaAntiHelium3Gen->Fill(EtaHelium3Gen);
		    fHistSEDAntiPairGen->Fill(RelativeMomentum);
		    fHistPtAntiParticlesGen->Fill(PtDeuteronGen,PtProtonGen);

		  } // end of "if antideuteron PDG"

	      } // end of stack loop 2

	  } // end of "if antiproton PDG"

      } // end of stack loop 1

    // loop for reconstructed particles
    for(int nTrack1 = 0; nTrack1 < fESD->GetNumberOfTracks(); nTrack1++)
      {

	AliESDtrack *Track1 = dynamic_cast<AliESDtrack*>(fESD->GetTrack(nTrack1));

	for(int nTrack2 = nTrack1 + 1; nTrack2 < fESD->GetNumberOfTracks(); nTrack2++)
	  {

	    AliESDtrack *Track2 = dynamic_cast<AliESDtrack*>(fESD->GetTrack(nTrack2));

	    if(!fESDtrackCutsProton->AcceptTrack(Track1)) continue;
	    if(!fESDtrackCutsDeuteron->AcceptTrack(Track2)) continue;

	    if(Track1->GetSign() < 1) continue;
	    if(Track2->GetSign() < 1) continue;

	    AliMCParticle *Particle1 = (AliMCParticle*)(mcEvent->GetTrack(TMath::Abs(Track1->GetLabel())));
	    AliMCParticle *Particle2 = (AliMCParticle*)(mcEvent->GetTrack(TMath::Abs(Track2->GetLabel())));

	    int Particle1PDG = Particle1->PdgCode();
	    int Particle2PDG = Particle2->PdgCode();
/*
	    if((Particle1PDG == TritonPDG)	|| (Particle1PDG == -TritonPDG))      continue;
	    if((Particle1PDG == Helium3PDG)	|| (Particle1PDG == -Helium3PDG))     continue;
	    if((Particle1PDG == HypertritonPDG) || (Particle1PDG == -HypertritonPDG)) continue;
	    if((Particle1PDG == PionPDG)	|| (Particle1PDG == -PionPDG))	      continue;
	    if((Particle1PDG == KaonPDG)	|| (Particle1PDG == -KaonPDG))	      continue;
	    if((Particle1PDG == ElectronPDG)	|| (Particle1PDG == -ElectronPDG))    continue;
	    if((Particle1PDG == MyonPDG)	|| (Particle1PDG == -MyonPDG))	      continue;

	    if((Particle2PDG == TritonPDG)	|| (Particle2PDG == -TritonPDG))      continue;
	    if((Particle2PDG == Helium3PDG)	|| (Particle2PDG == -Helium3PDG))     continue;
	    if((Particle2PDG == HypertritonPDG) || (Particle2PDG == -HypertritonPDG)) continue;
	    if((Particle2PDG == PionPDG)	|| (Particle2PDG == -PionPDG))	      continue;
	    if((Particle2PDG == KaonPDG)	|| (Particle2PDG == -KaonPDG))	      continue;
	    if((Particle2PDG == ElectronPDG)	|| (Particle2PDG == -ElectronPDG))    continue;
	    if((Particle2PDG == MyonPDG)	|| (Particle2PDG == -MyonPDG))	      continue;
*/
	    if(!(Particle1PDG == ProtonPDG))	continue;
	    if(!(Particle2PDG == DeuteronPDG))	continue;

	    AliPIDResponse::EDetPidStatus status1TPC = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTPC,Track1);
	    AliPIDResponse::EDetPidStatus status2TPC = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTPC,Track2);

	    AliPIDResponse::EDetPidStatus status1TOF = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTOF,Track1);
	    AliPIDResponse::EDetPidStatus status2TOF = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTOF,Track2);

	    if(status1TPC == AliPIDResponse::kDetPidOk) TPC1_OK = true;
	    if(status2TPC == AliPIDResponse::kDetPidOk) TPC2_OK = true;

	    if(status1TOF == AliPIDResponse::kDetPidOk) TOF1_OK = true;
	    if(status2TOF == AliPIDResponse::kDetPidOk) TOF2_OK = true;

	    if(!(TPC1_OK)) continue;
	    if(!(TPC2_OK)) continue;
	    if(!(TOF1_OK)) continue;
	    if(!(TOF2_OK)) continue;

	    nSigmaTPCproton = fPIDResponse->NumberOfSigmasTPC(Track1,AliPID::kProton);
	    nSigmaTPCdeuteron = fPIDResponse->NumberOfSigmasTPC(Track2,AliPID::kDeuteron);

	    nSigmaTOFproton = fPIDResponse->NumberOfSigmasTOF(Track1,AliPID::kProton);
	    nSigmaTOFdeuteron = fPIDResponse->NumberOfSigmasTOF(Track2,AliPID::kDeuteron);

	    nSigmaTPCTOFproton = TMath::Sqrt(pow(nSigmaTPCproton,2.) + pow(nSigmaTOFproton,2.));
	    nSigmaTPCTOFdeuteron = TMath::Sqrt(pow(nSigmaTPCdeuteron,2.) + pow(nSigmaTOFdeuteron,2.));

	    MomentumProtonRec[0] = Track1->Px();
	    MomentumProtonRec[1] = Track1->Py();
	    MomentumProtonRec[2] = Track1->Pz();

	    MomentumDeuteronRec[0] = Track2->Px();
	    MomentumDeuteronRec[1] = Track2->Py();
	    MomentumDeuteronRec[2] = Track2->Pz();

	    TLorentzVector LorentzVectorProtonRec;
	    TLorentzVector LorentzVectorDeuteronRec;
	    TLorentzVector LorentzVectorHelium3Rec;

	    LorentzVectorProtonRec.SetXYZM(MomentumProtonRec[0],MomentumProtonRec[1],MomentumProtonRec[2],ProtonMass);
	    LorentzVectorDeuteronRec.SetXYZM(MomentumDeuteronRec[0],MomentumDeuteronRec[1],MomentumDeuteronRec[2],DeuteronMass);

	    LorentzVectorHelium3Rec = LorentzVectorProtonRec + LorentzVectorDeuteronRec;

	    EtaProtonRec = LorentzVectorProtonRec.Eta();
	    EtaDeuteronRec = LorentzVectorDeuteronRec.Eta();
	    EtaHelium3Rec = LorentzVectorHelium3Rec.Eta();

	    PtProtonRec = LorentzVectorProtonRec.Pt();
	    PtDeuteronRec = LorentzVectorDeuteronRec.Pt();
	    PtHelium3Rec = LorentzVectorHelium3Rec.Pt();

	    if(PtProtonRec <= ProtonPtThreshold){

	      if(!(nSigmaTPCproton < 3.0)) continue;

	    }else{
	    
	      if(!(nSigmaTPCTOFproton < 3.0)) continue;

	    }


	    if(PtDeuteronRec <= DeuteronPtThreshold){

	      if(!(nSigmaTPCdeuteron < 3.0)) continue;

	    }else{
	    
	      if(!(nSigmaTPCTOFdeuteron < 3.0)) continue;

	    }

	    if((PtHelium3Rec < PairPtLimit1)	  || (PtHelium3Rec > PairPtLimit2)) continue;
	    if((!Particle1->IsPhysicalPrimary())  || (!Particle2->IsPhysicalPrimary())) continue;


	    RelativeMomentum = CalculateRelativeMomentum(LorentzVectorHelium3Rec,LorentzVectorProtonRec,LorentzVectorDeuteronRec);

	    fHistPtProtonRec->Fill(PtProtonRec);
	    fHistEtaProtonRec->Fill(EtaProtonRec);
	    fHistPtDeuteronRec->Fill(PtDeuteronRec);
	    fHistEtaDeuteronRec->Fill(EtaDeuteronRec);
	    fHistPtHelium3Rec->Fill(PtHelium3Rec);
	    fHistEtaHelium3Rec->Fill(EtaHelium3Rec);
	    fHistSEDPairRec->Fill(RelativeMomentum);
	    fHistPtParticlesRec->Fill(PtDeuteronRec,PtProtonRec);

	  } // end of track2 loop

      } // end of track1 loop


    // loop for reconstructed antiparticles
    for(int nTrack1 = 0; nTrack1 < fESD->GetNumberOfTracks(); nTrack1++)
      {

	AliESDtrack *Track1 = dynamic_cast<AliESDtrack*>(fESD->GetTrack(nTrack1));

	for(int nTrack2 = nTrack1 + 1; nTrack2 < fESD->GetNumberOfTracks(); nTrack2++)
	  {

	    AliESDtrack *Track2 = dynamic_cast<AliESDtrack*>(fESD->GetTrack(nTrack2));

	    if(!fESDtrackCutsProton->AcceptTrack(Track1)) continue;
	    if(!fESDtrackCutsDeuteron->AcceptTrack(Track2)) continue;

	    if(Track1->GetSign() > -1) continue;
	    if(Track2->GetSign() > -1) continue;

	    AliMCParticle *Particle1 = (AliMCParticle*)(mcEvent->GetTrack(TMath::Abs(Track1->GetLabel())));
	    AliMCParticle *Particle2 = (AliMCParticle*)(mcEvent->GetTrack(TMath::Abs(Track2->GetLabel())));

	    int Particle1PDG = Particle1->PdgCode();
	    int Particle2PDG = Particle2->PdgCode();
/*
	    if((Particle1PDG == TritonPDG)	|| (Particle1PDG == -TritonPDG))      continue;
	    if((Particle1PDG == Helium3PDG)	|| (Particle1PDG == -Helium3PDG))     continue;
	    if((Particle1PDG == HypertritonPDG) || (Particle1PDG == -HypertritonPDG)) continue;
	    if((Particle1PDG == PionPDG)	|| (Particle1PDG == -PionPDG))	      continue;
	    if((Particle1PDG == KaonPDG)	|| (Particle1PDG == -KaonPDG))	      continue;
	    if((Particle1PDG == ElectronPDG)	|| (Particle1PDG == -ElectronPDG))    continue;
	    if((Particle1PDG == MyonPDG)	|| (Particle1PDG == -MyonPDG))	      continue;

	    if((Particle2PDG == TritonPDG)	|| (Particle2PDG == -TritonPDG))      continue;
	    if((Particle2PDG == Helium3PDG)	|| (Particle2PDG == -Helium3PDG))     continue;
	    if((Particle2PDG == HypertritonPDG) || (Particle2PDG == -HypertritonPDG)) continue;
	    if((Particle2PDG == PionPDG)	|| (Particle2PDG == -PionPDG))	      continue;
	    if((Particle2PDG == KaonPDG)	|| (Particle2PDG == -KaonPDG))	      continue;
	    if((Particle2PDG == ElectronPDG)	|| (Particle2PDG == -ElectronPDG))    continue;
	    if((Particle2PDG == MyonPDG)	|| (Particle2PDG == -MyonPDG))	      continue;
*/
	    if(!(Particle1PDG == -ProtonPDG))	continue;
	    if(!(Particle2PDG == -DeuteronPDG))	continue;

	    AliPIDResponse::EDetPidStatus status1TPC = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTPC,Track1);
	    AliPIDResponse::EDetPidStatus status2TPC = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTPC,Track2);

	    AliPIDResponse::EDetPidStatus status1TOF = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTOF,Track1);
	    AliPIDResponse::EDetPidStatus status2TOF = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTOF,Track2);

	    if(status1TPC == AliPIDResponse::kDetPidOk) TPC1_OK = true;
	    if(status2TPC == AliPIDResponse::kDetPidOk) TPC2_OK = true;

	    if(status1TOF == AliPIDResponse::kDetPidOk) TOF1_OK = true;
	    if(status2TOF == AliPIDResponse::kDetPidOk) TOF2_OK = true;

	    if(!(TPC1_OK)) continue;
	    if(!(TPC2_OK)) continue;
	    if(!(TOF1_OK)) continue;
	    if(!(TOF2_OK)) continue;

	    nSigmaTPCproton = fPIDResponse->NumberOfSigmasTPC(Track1,AliPID::kProton);
	    nSigmaTPCdeuteron = fPIDResponse->NumberOfSigmasTPC(Track2,AliPID::kDeuteron);

	    nSigmaTOFproton = fPIDResponse->NumberOfSigmasTOF(Track1,AliPID::kProton);
	    nSigmaTOFdeuteron = fPIDResponse->NumberOfSigmasTOF(Track2,AliPID::kDeuteron);

	    nSigmaTPCTOFproton = TMath::Sqrt(pow(nSigmaTPCproton,2.) + pow(nSigmaTOFproton,2.));
	    nSigmaTPCTOFdeuteron = TMath::Sqrt(pow(nSigmaTPCdeuteron,2.) + pow(nSigmaTOFdeuteron,2.));

	    MomentumProtonRec[0] = Track1->Px();
	    MomentumProtonRec[1] = Track1->Py();
	    MomentumProtonRec[2] = Track1->Pz();

	    MomentumDeuteronRec[0] = Track2->Px();
	    MomentumDeuteronRec[1] = Track2->Py();
	    MomentumDeuteronRec[2] = Track2->Pz();

	    TLorentzVector LorentzVectorProtonRec;
	    TLorentzVector LorentzVectorDeuteronRec;
	    TLorentzVector LorentzVectorHelium3Rec;

	    LorentzVectorProtonRec.SetXYZM(MomentumProtonRec[0],MomentumProtonRec[1],MomentumProtonRec[2],ProtonMass);
	    LorentzVectorDeuteronRec.SetXYZM(MomentumDeuteronRec[0],MomentumDeuteronRec[1],MomentumDeuteronRec[2],DeuteronMass);

	    LorentzVectorHelium3Rec = LorentzVectorProtonRec + LorentzVectorDeuteronRec;

	    EtaProtonRec = LorentzVectorProtonRec.Eta();
	    EtaDeuteronRec = LorentzVectorDeuteronRec.Eta();
	    EtaHelium3Rec = LorentzVectorHelium3Rec.Eta();

	    PtProtonRec = LorentzVectorProtonRec.Pt();
	    PtDeuteronRec = LorentzVectorDeuteronRec.Pt();
	    PtHelium3Rec = LorentzVectorHelium3Rec.Pt();

	    if(PtProtonRec <= ProtonPtThreshold){

	      if(!(nSigmaTPCproton < 3.0)) continue;

	    }else{
	    
	      if(!(nSigmaTPCTOFproton < 3.0)) continue;

	    }


	    if(PtDeuteronRec <= DeuteronPtThreshold){

	      if(!(nSigmaTPCdeuteron < 3.0)) continue;

	    }else{
	    
	      if(!(nSigmaTPCTOFdeuteron < 3.0)) continue;

	    }

	    if((PtHelium3Rec < PairPtLimit1)	  || (PtHelium3Rec > PairPtLimit2)) continue;
	    if((!Particle1->IsPhysicalPrimary())  || (!Particle2->IsPhysicalPrimary())) continue;


	    RelativeMomentum = CalculateRelativeMomentum(LorentzVectorHelium3Rec,LorentzVectorProtonRec,LorentzVectorDeuteronRec);

	    fHistPtAntiProtonRec->Fill(PtProtonRec);
	    fHistEtaAntiProtonRec->Fill(EtaProtonRec);
	    fHistPtAntiDeuteronRec->Fill(PtDeuteronRec);
	    fHistEtaAntiDeuteronRec->Fill(EtaDeuteronRec);
	    fHistPtAntiHelium3Rec->Fill(PtHelium3Rec);
	    fHistEtaAntiHelium3Rec->Fill(EtaHelium3Rec);
	    fHistSEDAntiPairRec->Fill(RelativeMomentum);
	    fHistPtAntiParticlesRec->Fill(PtDeuteronRec,PtProtonRec);

	  } // end of track2 loop

      } // end of track1 loop


  fHistEventCounter->Fill(0.5);

  } // end of "if is MC"

  PostData(1,fHistList);


} // end of UserExec




void AliAnalysisTaskDeuteronProtonEfficiency::Terminate(Option_t *)
{



} // end of Terminate



double AliAnalysisTaskDeuteronProtonEfficiency::CalculateRelativeMomentum(TLorentzVector &Pair, TLorentzVector &Part1, TLorentzVector &Part2)
{

  double beta = Pair.Beta();
  double betaX = beta * cos(Pair.Phi()) * sin(Pair.Theta());
  double betaY = beta * sin(Pair.Phi()) * sin(Pair.Theta());
  double betaZ = beta * cos(Pair.Theta());

  TLorentzVector Part1CMS = Part1;
  TLorentzVector Part2CMS = Part2;

  Part1CMS.Boost(-betaX,-betaY,-betaZ);
  Part2CMS.Boost(-betaX,-betaY,-betaZ);

  TLorentzVector RelativeMomentum = Part1CMS - Part2CMS;

  return 0.5 * RelativeMomentum.P();

} // end of CalculateRelativeMomentum


