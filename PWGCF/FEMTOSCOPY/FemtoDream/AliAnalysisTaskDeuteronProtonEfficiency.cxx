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
  fHistPtProtonGenPDG(0),
  fHistPtProtonGenPrimary(0),
  fHistPtProtonGenEtaCut(0),
  fHistPtProtonGenPairPtCut(0),
  fHistEtaProtonGenPDG(0),
  fHistEtaProtonGenPrimary(0),
  fHistEtaProtonGenEtaCut(0),
  fHistEtaProtonGenPairPtCut(0),
  fHistPtAntiProtonGenPDG(0),
  fHistPtAntiProtonGenPrimary(0),
  fHistPtAntiProtonGenEtaCut(0),
  fHistPtAntiProtonGenPairPtCut(0),
  fHistEtaAntiProtonGenPDG(0),
  fHistEtaAntiProtonGenPrimary(0),
  fHistEtaAntiProtonGenEtaCut(0),
  fHistEtaAntiProtonGenPairPtCut(0),
  fHistPtDeuteronGenPDG(0),
  fHistPtDeuteronGenPrimary(0),
  fHistPtDeuteronGenEtaCut(0),
  fHistPtDeuteronGenPairPtCut(0),
  fHistEtaDeuteronGenPDG(0),
  fHistEtaDeuteronGenPrimary(0),
  fHistEtaDeuteronGenEtaCut(0),
  fHistEtaDeuteronGenPairPtCut(0),
  fHistPtAntiDeuteronGenPDG(0),
  fHistPtAntiDeuteronGenPrimary(0),
  fHistPtAntiDeuteronGenEtaCut(0),
  fHistPtAntiDeuteronGenPairPtCut(0),
  fHistEtaAntiDeuteronGenPDG(0),
  fHistEtaAntiDeuteronGenPrimary(0),
  fHistEtaAntiDeuteronGenEtaCut(0),
  fHistEtaAntiDeuteronGenPairPtCut(0),
  fHistPtHelium3GenPairPtCut(0),
  fHistEtaHelium3GenPairPtCut(0),
  fHistPtParticlesGen(0),
  fHistSEDPairGen(0),
  fHistPtAntiHelium3GenPairPtCut(0),
  fHistEtaAntiHelium3GenPairPtCut(0),
  fHistPtAntiParticlesGen(0),
  fHistSEDAntiPairGen(0),
  fHistPtProtonRecPDG(0),
  fHistPtProtonRecPrimary(0),
  fHistPtProtonRecPairPtCut(0),
  fHistPtProtonRecTrackCuts(0),
  fHistEtaProtonRecPDG(0),
  fHistEtaProtonRecPrimary(0),
  fHistEtaProtonRecTrackCuts(0),
  fHistEtaProtonRecPairPtCut(0),
  fHistPtAntiProtonRecPDG(0),
  fHistPtAntiProtonRecPrimary(0),
  fHistPtAntiProtonRecTrackCuts(0),
  fHistPtAntiProtonRecPairPtCut(0),
  fHistEtaAntiProtonRecPDG(0),
  fHistEtaAntiProtonRecPrimary(0),
  fHistEtaAntiProtonRecTrackCuts(0),
  fHistEtaAntiProtonRecPairPtCut(0),
  fHistPtDeuteronRecPDG(0),
  fHistPtDeuteronRecPrimary(0),
  fHistPtDeuteronRecTrackCuts(0),
  fHistPtDeuteronRecPairPtCut(0),
  fHistEtaDeuteronRecPDG(0),
  fHistEtaDeuteronRecPrimary(0),
  fHistEtaDeuteronRecTrackCuts(0),
  fHistEtaDeuteronRecPairPtCut(0),
  fHistPtAntiDeuteronRecPDG(0),
  fHistPtAntiDeuteronRecPrimary(0),
  fHistPtAntiDeuteronRecTrackCuts(0),
  fHistPtAntiDeuteronRecPairPtCut(0),
  fHistEtaAntiDeuteronRecPDG(0),
  fHistEtaAntiDeuteronRecPrimary(0),
  fHistEtaAntiDeuteronRecTrackCuts(0),
  fHistEtaAntiDeuteronRecPairPtCut(0),
  fHistPtHelium3RecPairPtCut(0),
  fHistEtaHelium3RecPairPtCut(0),
  fHistPtParticlesRec(0),
  fHistSEDPairRec(0),
  fHistPtAntiHelium3RecPairPtCut(0),
  fHistEtaAntiHelium3RecPairPtCut(0),
  fHistPtAntiParticlesRec(0),
  fHistSEDAntiPairRec(0),
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
  fHistPtProtonGenPDG(0),
  fHistPtProtonGenPrimary(0),
  fHistPtProtonGenEtaCut(0),
  fHistPtProtonGenPairPtCut(0),
  fHistEtaProtonGenPDG(0),
  fHistEtaProtonGenPrimary(0),
  fHistEtaProtonGenEtaCut(0),
  fHistEtaProtonGenPairPtCut(0),
  fHistPtAntiProtonGenPDG(0),
  fHistPtAntiProtonGenPrimary(0),
  fHistPtAntiProtonGenEtaCut(0),
  fHistPtAntiProtonGenPairPtCut(0),
  fHistEtaAntiProtonGenPDG(0),
  fHistEtaAntiProtonGenPrimary(0),
  fHistEtaAntiProtonGenEtaCut(0),
  fHistEtaAntiProtonGenPairPtCut(0),
  fHistPtDeuteronGenPDG(0),
  fHistPtDeuteronGenPrimary(0),
  fHistPtDeuteronGenEtaCut(0),
  fHistPtDeuteronGenPairPtCut(0),
  fHistEtaDeuteronGenPDG(0),
  fHistEtaDeuteronGenPrimary(0),
  fHistEtaDeuteronGenEtaCut(0),
  fHistEtaDeuteronGenPairPtCut(0),
  fHistPtAntiDeuteronGenPDG(0),
  fHistPtAntiDeuteronGenPrimary(0),
  fHistPtAntiDeuteronGenEtaCut(0),
  fHistPtAntiDeuteronGenPairPtCut(0),
  fHistEtaAntiDeuteronGenPDG(0),
  fHistEtaAntiDeuteronGenPrimary(0),
  fHistEtaAntiDeuteronGenEtaCut(0),
  fHistEtaAntiDeuteronGenPairPtCut(0),
  fHistPtHelium3GenPairPtCut(0),
  fHistEtaHelium3GenPairPtCut(0),
  fHistPtParticlesGen(0),
  fHistSEDPairGen(0),
  fHistPtAntiHelium3GenPairPtCut(0),
  fHistEtaAntiHelium3GenPairPtCut(0),
  fHistPtAntiParticlesGen(0),
  fHistSEDAntiPairGen(0),
  fHistPtProtonRecPDG(0),
  fHistPtProtonRecPrimary(0),
  fHistPtProtonRecTrackCuts(0),
  fHistPtProtonRecPairPtCut(0),
  fHistEtaProtonRecPDG(0),
  fHistEtaProtonRecPrimary(0),
  fHistEtaProtonRecPairPtCut(0),
  fHistEtaProtonRecTrackCuts(0),
  fHistPtAntiProtonRecPDG(0),
  fHistPtAntiProtonRecPrimary(0),
  fHistPtAntiProtonRecTrackCuts(0),
  fHistPtAntiProtonRecPairPtCut(0),
  fHistEtaAntiProtonRecPDG(0),
  fHistEtaAntiProtonRecPrimary(0),
  fHistEtaAntiProtonRecTrackCuts(0),
  fHistEtaAntiProtonRecPairPtCut(0),
  fHistPtDeuteronRecPDG(0),
  fHistPtDeuteronRecPrimary(0),
  fHistPtDeuteronRecTrackCuts(0),
  fHistPtDeuteronRecPairPtCut(0),
  fHistEtaDeuteronRecPDG(0),
  fHistEtaDeuteronRecPrimary(0),
  fHistEtaDeuteronRecTrackCuts(0),
  fHistEtaDeuteronRecPairPtCut(0),
  fHistPtAntiDeuteronRecPDG(0),
  fHistPtAntiDeuteronRecPrimary(0),
  fHistPtAntiDeuteronRecTrackCuts(0),
  fHistPtAntiDeuteronRecPairPtCut(0),
  fHistEtaAntiDeuteronRecPDG(0),
  fHistEtaAntiDeuteronRecPrimary(0),
  fHistEtaAntiDeuteronRecTrackCuts(0),
  fHistEtaAntiDeuteronRecPairPtCut(0),
  fHistPtHelium3RecPairPtCut(0),
  fHistEtaHelium3RecPairPtCut(0),
  fHistPtParticlesRec(0),
  fHistSEDPairRec(0),
  fHistPtAntiHelium3RecPairPtCut(0),
  fHistEtaAntiHelium3RecPairPtCut(0),
  fHistPtAntiParticlesRec(0),
  fHistSEDAntiPairRec(0),
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
//  fESDtrackCutsProton->SetMaxDCAToVertexXY(0.1);
//  fESDtrackCutsProton->SetMaxDCAToVertexZ(0.2);

  // Deuteron cuts
  fESDtrackCutsDeuteron = new AliESDtrackCuts("AliESDtrackCutsDeuteron","AliESDtrackCutsDeuteron");
  fESDtrackCutsDeuteron->SetEtaRange(-0.8,0.8);
  fESDtrackCutsDeuteron->SetPtRange(0.4,4.0);
  fESDtrackCutsDeuteron->SetMinNClustersTPC(80);
  fESDtrackCutsDeuteron->SetMinNCrossedRowsTPC(70);
  fESDtrackCutsDeuteron->SetMinRatioCrossedRowsOverFindableClustersTPC(0.83);
  fESDtrackCutsDeuteron->SetAcceptSharedTPCClusters(false);
//  fESDtrackCutsDeuteron->SetMaxDCAToVertexXY(0.1);
//  fESDtrackCutsDeuteron->SetMaxDCAToVertexZ(0.2);
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
  
  // generated protons
  fHistPtProtonGenPDG = new TH1F("fHistPtProtonGenPDG","generated p (PDG)",240,0.0,6.0);
  fHistPtProtonGenPDG->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  fHistPtProtonGenPDG->GetYaxis()->SetTitle("entries");

  fHistPtProtonGenPrimary = new TH1F("fHistPtProtonGenPrimary","generated p (PDG + primary)",240,0.0,6.0);
  fHistPtProtonGenPrimary->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  fHistPtProtonGenPrimary->GetYaxis()->SetTitle("entries");

  fHistPtProtonGenEtaCut = new TH1F("fHistPtProtonGenEtaCut","generated p (PDG + primary + #it{#eta} cut)",240,0.0,6.0);
  fHistPtProtonGenEtaCut->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  fHistPtProtonGenEtaCut->GetYaxis()->SetTitle("entries");

  fHistPtProtonGenPairPtCut = new TH1F("fHistPtProtonGenPairPtCut","generated p (PDG + primary + #it{#eta} cut + pair #it{p}_{T} cut)",240,0.0,6.0);
  fHistPtProtonGenPairPtCut->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  fHistPtProtonGenPairPtCut->GetYaxis()->SetTitle("entries");

  fHistEtaProtonGenPDG = new TH1F("fHistEtaProtonGenPDG","generated p (PDG)",200,-1.5,1.5);
  fHistEtaProtonGenPDG->GetXaxis()->SetTitle("#it{#eta}");
  fHistEtaProtonGenPDG->GetYaxis()->SetTitle("entries");

  fHistEtaProtonGenPrimary = new TH1F("fHistEtaProtonGenPrimary","generated p (PDG + primary)",200,-1.5,1.5);
  fHistEtaProtonGenPrimary->GetXaxis()->SetTitle("#it{#eta}");
  fHistEtaProtonGenPrimary->GetYaxis()->SetTitle("entries");

  fHistEtaProtonGenEtaCut = new TH1F("fHistEtaProtonGenEtaCut","generated p (PDG + primary + #it{#eta} cut)",200,-1.5,1.5);
  fHistEtaProtonGenEtaCut->GetXaxis()->SetTitle("#it{#eta}");
  fHistEtaProtonGenEtaCut->GetYaxis()->SetTitle("entries");

  fHistEtaProtonGenPairPtCut = new TH1F("fHistEtaProtonGenPairPtCut","generated p (PDG + primary + #it{#eta} cut + pair #it{p}_{T} cut)",200,-1.5,1.5);
  fHistEtaProtonGenPairPtCut->GetXaxis()->SetTitle("#it{#eta}");
  fHistEtaProtonGenPairPtCut->GetYaxis()->SetTitle("entries");
  
  // generated antiprotons
  fHistPtAntiProtonGenPDG = new TH1F("fHistPtAntiProtonGenPDG","generated #bar{p} (PDG)",240,0.0,6.0);
  fHistPtAntiProtonGenPDG->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  fHistPtAntiProtonGenPDG->GetYaxis()->SetTitle("entries");

  fHistPtAntiProtonGenPrimary = new TH1F("fHistPtAntiProtonGenPrimary","generated #bar{p} (PDG + primary)",240,0.0,6.0);
  fHistPtAntiProtonGenPrimary->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  fHistPtAntiProtonGenPrimary->GetYaxis()->SetTitle("entries");

  fHistPtAntiProtonGenEtaCut = new TH1F("fHistPtAntiProtonGenEtaCut","generated #bar{p} (PDG + primary + #it{#eta} cut)",240,0.0,6.0);
  fHistPtAntiProtonGenEtaCut->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  fHistPtAntiProtonGenEtaCut->GetYaxis()->SetTitle("entries");

  fHistPtAntiProtonGenPairPtCut = new TH1F("fHistPtAntiProtonGenPairPtCut","generated #bar{p} (PDG + primary + #it{#eta} cut + pair #it{p}_{T} cut)",240,0.0,6.0);
  fHistPtAntiProtonGenPairPtCut->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  fHistPtAntiProtonGenPairPtCut->GetYaxis()->SetTitle("entries");

  fHistEtaAntiProtonGenPDG = new TH1F("fHistEtaAntiProtonGenPDG","generated #bar{p} (PDG)",200,-1.5,1.5);
  fHistEtaAntiProtonGenPDG->GetXaxis()->SetTitle("#it{#eta}");
  fHistEtaAntiProtonGenPDG->GetYaxis()->SetTitle("entries");

  fHistEtaAntiProtonGenPrimary = new TH1F("fHistEtaAntiProtonGenPrimary","generated #bar{p} (PDG + primary)",200,-1.5,1.5);
  fHistEtaAntiProtonGenPrimary->GetXaxis()->SetTitle("#it{#eta}");
  fHistEtaAntiProtonGenPrimary->GetYaxis()->SetTitle("entries");

  fHistEtaAntiProtonGenEtaCut = new TH1F("fHistEtaAntiProtonGenEtaCut","generated #bar{p} (PDG + primary + #it{#eta} cut)",200,-1.5,1.5);
  fHistEtaAntiProtonGenEtaCut->GetXaxis()->SetTitle("#it{#eta}");
  fHistEtaAntiProtonGenEtaCut->GetYaxis()->SetTitle("entries");

  fHistEtaAntiProtonGenPairPtCut = new TH1F("fHistEtaAntiProtonGenPairPtCut","generated #bar{p} (PDG + primary + #it{#eta} cut + pair #it{p}_{T} cut)",200,-1.5,1.5);
  fHistEtaAntiProtonGenPairPtCut->GetXaxis()->SetTitle("#it{#eta}");
  fHistEtaAntiProtonGenPairPtCut->GetYaxis()->SetTitle("entries");

  // generated deuterons
  fHistPtDeuteronGenPDG = new TH1F("fHistPtDeuteronGenPDG","generated d (PDG)",240,0.0,6.0);
  fHistPtDeuteronGenPDG->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  fHistPtDeuteronGenPDG->GetYaxis()->SetTitle("entries");

  fHistPtDeuteronGenPrimary = new TH1F("fHistPtDeuteronGenPrimary","generated d (PDG + primary)",240,0.0,6.0);
  fHistPtDeuteronGenPrimary->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  fHistPtDeuteronGenPrimary->GetYaxis()->SetTitle("entries");

  fHistPtDeuteronGenEtaCut = new TH1F("fHistPtDeuteronGenEtaCut","generated d (PDG + primary + #it{#eta} cut)",240,0.0,6.0);
  fHistPtDeuteronGenEtaCut->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  fHistPtDeuteronGenEtaCut->GetYaxis()->SetTitle("entries");

  fHistPtDeuteronGenPairPtCut = new TH1F("fHistPtDeuteronGenPairPtCut","generated d (PDG + primary + #it{#eta} cut + pair #it{p}_{T} cut)",240,0.0,6.0);
  fHistPtDeuteronGenPairPtCut->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  fHistPtDeuteronGenPairPtCut->GetYaxis()->SetTitle("entries");

  fHistEtaDeuteronGenPDG = new TH1F("fHistEtaDeuteronGenPDG","generated d (PDG)",200,-1.5,1.5);
  fHistEtaDeuteronGenPDG->GetXaxis()->SetTitle("#it{#eta}");
  fHistEtaDeuteronGenPDG->GetYaxis()->SetTitle("entries");

  fHistEtaDeuteronGenPrimary = new TH1F("fHistEtaDeuteronGenPrimary","generated d (PDG + primary)",200,-1.5,1.5);
  fHistEtaDeuteronGenPrimary->GetXaxis()->SetTitle("#it{#eta}");
  fHistEtaDeuteronGenPrimary->GetYaxis()->SetTitle("entries");

  fHistEtaDeuteronGenEtaCut = new TH1F("fHistEtaDeuteronGenEtaCut","generated d (PDG + primary + #it{#eta} cut)",200,-1.5,1.5);
  fHistEtaDeuteronGenEtaCut->GetXaxis()->SetTitle("#it{#eta}");
  fHistEtaDeuteronGenEtaCut->GetYaxis()->SetTitle("entries");

  fHistEtaDeuteronGenPairPtCut = new TH1F("fHistEtaDeuteronGenPairPtCut","generated d (PDG + primary + #it{#eta} cut + pair #it{p}_{T} cut)",200,-1.5,1.5);
  fHistEtaDeuteronGenPairPtCut->GetXaxis()->SetTitle("#it{#eta}");
  fHistEtaDeuteronGenPairPtCut->GetYaxis()->SetTitle("entries");

  // generated antideuterons
  fHistPtAntiDeuteronGenPDG = new TH1F("fHistPtAntiDeuteronGenPDG","generated #bar{d} (PDG)",240,0.0,6.0);
  fHistPtAntiDeuteronGenPDG->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  fHistPtAntiDeuteronGenPDG->GetYaxis()->SetTitle("entries");

  fHistPtAntiDeuteronGenPrimary = new TH1F("fHistPtAntiDeuteronGenPrimary","generated #bar{d} (PDG + primary)",240,0.0,6.0);
  fHistPtAntiDeuteronGenPrimary->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  fHistPtAntiDeuteronGenPrimary->GetYaxis()->SetTitle("entries");

  fHistPtAntiDeuteronGenEtaCut = new TH1F("fHistPtAntiDeuteronGenEtaCut","generated #bar{d} (PDG + primary + #it{#eta} cut)",240,0.0,6.0);
  fHistPtAntiDeuteronGenEtaCut->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  fHistPtAntiDeuteronGenEtaCut->GetYaxis()->SetTitle("entries");

  fHistPtAntiDeuteronGenPairPtCut = new TH1F("fHistPtAntiDeuteronGenPairPtCut","generated #bar{d} (PDG + primary + #it{#eta} cut + pair #it{p}_{T} cut)",240,0.0,6.0);
  fHistPtAntiDeuteronGenPairPtCut->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  fHistPtAntiDeuteronGenPairPtCut->GetYaxis()->SetTitle("entries");

  fHistEtaAntiDeuteronGenPDG = new TH1F("fHistEtaAntiDeuteronGenPDG","generated #bar{d} (PDG)",200,-1.5,1.5);
  fHistEtaAntiDeuteronGenPDG->GetXaxis()->SetTitle("#it{#eta}");
  fHistEtaAntiDeuteronGenPDG->GetYaxis()->SetTitle("entries");

  fHistEtaAntiDeuteronGenPrimary = new TH1F("fHistEtaAntiDeuteronGenPrimary","generated #bar{d} (PDG + primary)",200,-1.5,1.5);
  fHistEtaAntiDeuteronGenPrimary->GetXaxis()->SetTitle("#it{#eta}");
  fHistEtaAntiDeuteronGenPrimary->GetYaxis()->SetTitle("entries");

  fHistEtaAntiDeuteronGenEtaCut = new TH1F("fHistEtaAntiDeuteronGenEtaCut","generated #bar{d} (PDG + primary + #it{#eta} cut)",200,-1.5,1.5);
  fHistEtaAntiDeuteronGenEtaCut->GetXaxis()->SetTitle("#it{#eta}");
  fHistEtaAntiDeuteronGenEtaCut->GetYaxis()->SetTitle("entries");

  fHistEtaAntiDeuteronGenPairPtCut = new TH1F("fHistEtaAntiDeuteronGenPairPtCut","generated #bar{d} (PDG + primary + #it{#eta} cut + pair #it{p}_{T} cut)",200,-1.5,1.5);
  fHistEtaAntiDeuteronGenPairPtCut->GetXaxis()->SetTitle("#it{#eta}");
  fHistEtaAntiDeuteronGenPairPtCut->GetYaxis()->SetTitle("entries");

  // generated pairs
  fHistPtHelium3GenPairPtCut = new TH1F("fHistPtHelium3GenPairPtCut","generated pairs (pair #it{p}_{T} cut)",240,0.0,6.0);
  fHistPtHelium3GenPairPtCut->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  fHistPtHelium3GenPairPtCut->GetYaxis()->SetTitle("entries");

  fHistEtaHelium3GenPairPtCut = new TH1F("fHistEtaHelium3GenPairPtCut","generated pairs (pair #it{p}_{T} cut)",200,-1.5,1.5);
  fHistEtaHelium3GenPairPtCut->GetXaxis()->SetTitle("#it{#eta}");
  fHistEtaHelium3GenPairPtCut->GetYaxis()->SetTitle("entries");

  fHistSEDPairGen = new TH1F("fHistSEDPairGen","same-event distribution of generated d-p pairs",750,0.0,3.0);
  fHistSEDPairGen->GetXaxis()->SetTitle("#it{k*} (GeV/#it{c})");
  fHistSEDPairGen->GetYaxis()->SetTitle("entries");

  fHistPtParticlesGen = new TH2F("fHistPtParticlesGen","#it{p}_{T} distribution of generated p and d",240,0.0,6.0,240,0.0,6.0);
  fHistPtParticlesGen->GetXaxis()->SetTitle("#it{p}_{T} of deuterons (GeV/#it{c})");
  fHistPtParticlesGen->GetYaxis()->SetTitle("#it{p}_{T} of protons (GeV/#it{c})");

  // generated antipairs
  fHistPtAntiHelium3GenPairPtCut = new TH1F("fHistPtAntiHelium3GenPairPtCut","generated antipairs (pair #it{p}_{T} cut)",240,0.0,6.0);
  fHistPtAntiHelium3GenPairPtCut->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  fHistPtAntiHelium3GenPairPtCut->GetYaxis()->SetTitle("entries");

  fHistEtaAntiHelium3GenPairPtCut = new TH1F("fHistEtaAntiHelium3GenPairPtCut","generated antipairs (#it{p}_{T} cut)",200,-1.5,1.5);
  fHistEtaAntiHelium3GenPairPtCut->GetXaxis()->SetTitle("#it{#eta}");
  fHistEtaAntiHelium3GenPairPtCut->GetYaxis()->SetTitle("entries");

  fHistSEDAntiPairGen = new TH1F("fHistSEDAntiPairGen","same-event distribution of generated #bar{d}-#bar{p} pairs",750,0.0,3.0);
  fHistSEDAntiPairGen->GetXaxis()->SetTitle("#it{k*} (GeV/#it{c})");
  fHistSEDAntiPairGen->GetYaxis()->SetTitle("entries");

  fHistPtAntiParticlesGen = new TH2F("fHistPtAntiParticlesGen","#it{p}_{T} distribution of generated #bar{p} and #bar{d}",240,0.0,6.0,240,0.0,6.0);
  fHistPtAntiParticlesGen->GetXaxis()->SetTitle("#it{p}_{T} of antideuterons (GeV/#it{c})");
  fHistPtAntiParticlesGen->GetYaxis()->SetTitle("#it{p}_{T} of antiprotons (GeV/#it{c})");

  //  histograms for reconstructed particles
  
  // reconstructed protons
  fHistPtProtonRecPDG = new TH1F("fHistPtProtonRecPDG","reconstructed p (PDG)",240,0.0,6.0);
  fHistPtProtonRecPDG->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  fHistPtProtonRecPDG->GetYaxis()->SetTitle("entries");

  fHistPtProtonRecPrimary = new TH1F("fHistPtProtonRecPrimary","reconstructed p (PDG + primary)",240,0.0,6.0);
  fHistPtProtonRecPrimary->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  fHistPtProtonRecPrimary->GetYaxis()->SetTitle("entries");

  fHistPtProtonRecTrackCuts = new TH1F("fHistPtProtonRecTrackCuts","reconstructed p (PDG + Primary + track cuts)",240,0.0,6.0);
  fHistPtProtonRecTrackCuts->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  fHistPtProtonRecTrackCuts->GetYaxis()->SetTitle("entries");

  fHistPtProtonRecPairPtCut = new TH1F("fHistPtProtonRecPairPtCut","reconstructed p (PDG + primary + track cuts + pair #it{p}_{T})",240,0.0,6.0);
  fHistPtProtonRecPairPtCut->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  fHistPtProtonRecPairPtCut->GetYaxis()->SetTitle("entries");

  fHistEtaProtonRecPDG = new TH1F("fHistEtaProtonRecPDG","reconstructed p (PDG)",200,-1.5,1.5);
  fHistEtaProtonRecPDG->GetXaxis()->SetTitle("#it{#eta}");
  fHistEtaProtonRecPDG->GetYaxis()->SetTitle("entries");

  fHistEtaProtonRecPrimary = new TH1F("fHistEtaProtonRecPrimary","reconstructed p (PDG + primary)",200,-1.5,1.5);
  fHistEtaProtonRecPrimary->GetXaxis()->SetTitle("#it{#eta}");
  fHistEtaProtonRecPrimary->GetYaxis()->SetTitle("entries");

  fHistEtaProtonRecTrackCuts = new TH1F("fHistEtaProtonRecTrackCuts","reconstructed p (PDG + primary + track cuts)",200,-1.5,1.5);
  fHistEtaProtonRecTrackCuts->GetXaxis()->SetTitle("#it{#eta}");
  fHistEtaProtonRecTrackCuts->GetYaxis()->SetTitle("entries");

  fHistEtaProtonRecPairPtCut = new TH1F("fHistEtaProtonRecPairPtCut","reconstructed p (PDG + primary + track cuts + pair #it{p}_{T})",200,-1.5,1.5);
  fHistEtaProtonRecPairPtCut->GetXaxis()->SetTitle("#it{#eta}");
  fHistEtaProtonRecPairPtCut->GetYaxis()->SetTitle("entries");

  // reconstructed antiprotons
  fHistPtAntiProtonRecPDG = new TH1F("fHistPtAntiProtonRecPDG","reconstructed #bar{p} (PDG)",240,0.0,6.0);
  fHistPtAntiProtonRecPDG->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  fHistPtAntiProtonRecPDG->GetYaxis()->SetTitle("entries");

  fHistPtAntiProtonRecPrimary = new TH1F("fHistPtAntiProtonRecPrimary","reconstructed #bar{p} (PDG + primary)",240,0.0,6.0);
  fHistPtAntiProtonRecPrimary->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  fHistPtAntiProtonRecPrimary->GetYaxis()->SetTitle("entries");

  fHistPtAntiProtonRecTrackCuts = new TH1F("fHistPtAntiProtonRecTrackCuts","reconstructed #bar{p} (PDG + primary + track cuts)",240,0.0,6.0);
  fHistPtAntiProtonRecTrackCuts->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  fHistPtAntiProtonRecTrackCuts->GetYaxis()->SetTitle("entries");

  fHistPtAntiProtonRecPairPtCut = new TH1F("fHistPtAntiProtonRecPairPtCut","reconstructed #bar{p} (PDG + primary + track cuts  + pair #it{p}_{T} cut)",240,0.0,6.0);
  fHistPtAntiProtonRecPairPtCut->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  fHistPtAntiProtonRecPairPtCut->GetYaxis()->SetTitle("entries");

  fHistEtaAntiProtonRecPDG = new TH1F("fHistEtaAntiProtonRecPDG","reconstruced #bar{p} (PDG)",200,-1.5,1.5);
  fHistEtaAntiProtonRecPDG->GetXaxis()->SetTitle("#it{#eta}");
  fHistEtaAntiProtonRecPDG->GetYaxis()->SetTitle("entries");

  fHistEtaAntiProtonRecPrimary = new TH1F("fHistEtaAntiProtonRecPrimary","reconstruced #bar{p} (PDG + primary)",200,-1.5,1.5);
  fHistEtaAntiProtonRecPrimary->GetXaxis()->SetTitle("#it{#eta}");
  fHistEtaAntiProtonRecPrimary->GetYaxis()->SetTitle("entries");

  fHistEtaAntiProtonRecTrackCuts = new TH1F("fHistEtaAntiProtonRecTrackCuts","reconstruced #bar{p} (PDG + primary + track cuts)",200,-1.5,1.5);
  fHistEtaAntiProtonRecTrackCuts->GetXaxis()->SetTitle("#it{#eta}");
  fHistEtaAntiProtonRecTrackCuts->GetYaxis()->SetTitle("entries");

  fHistEtaAntiProtonRecPairPtCut = new TH1F("fHistEtaAntiProtonRecPairPtCut","reconstruced #bar{p} (PDG + primary + track cuts + pair #it{p}_{T} cut)",200,-1.5,1.5);
  fHistEtaAntiProtonRecPairPtCut->GetXaxis()->SetTitle("#it{#eta}");
  fHistEtaAntiProtonRecPairPtCut->GetYaxis()->SetTitle("entries");

  // reconstructed deuterons
  fHistPtDeuteronRecPDG = new TH1F("fHistPtDeuteronRecPDG","reconstructed d (PDG)",240,0.0,6.0);
  fHistPtDeuteronRecPDG->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  fHistPtDeuteronRecPDG->GetYaxis()->SetTitle("entries");

  fHistPtDeuteronRecPrimary = new TH1F("fHistPtDeuteronRecPrimary","reconstructed d (PDG + primary)",240,0.0,6.0);
  fHistPtDeuteronRecPrimary->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  fHistPtDeuteronRecPrimary->GetYaxis()->SetTitle("entries");

  fHistPtDeuteronRecTrackCuts = new TH1F("fHistPtDeuteronRecTrackCuts","reconstructed d (PDG + primary + track cuts)",240,0.0,6.0);
  fHistPtDeuteronRecTrackCuts->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  fHistPtDeuteronRecTrackCuts->GetYaxis()->SetTitle("entries");

  fHistPtDeuteronRecPairPtCut = new TH1F("fHistPtDeuteronRecPairPtCut","reconstructed d (PDG + primary + track cuts + pair #it{p}_{T} cut)",240,0.0,6.0);
  fHistPtDeuteronRecPairPtCut->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  fHistPtDeuteronRecPairPtCut->GetYaxis()->SetTitle("entries");

  fHistEtaDeuteronRecPDG = new TH1F("fHistEtaDeuteronRecPDG","reconstructed d (PDG)",200,-1.5,1.5);
  fHistEtaDeuteronRecPDG->GetXaxis()->SetTitle("#it{#eta}");
  fHistEtaDeuteronRecPDG->GetYaxis()->SetTitle("entries");

  fHistEtaDeuteronRecPrimary = new TH1F("fHistEtaDeuteronRecPrimary","reconstructed d (PDG + primary)",200,-1.5,1.5);
  fHistEtaDeuteronRecPrimary->GetXaxis()->SetTitle("#it{#eta}");
  fHistEtaDeuteronRecPrimary->GetYaxis()->SetTitle("entries");

  fHistEtaDeuteronRecTrackCuts = new TH1F("fHistEtaDeuteronRecTrackCuts","reconstructed d (PDG + primary + track cuts)",200,-1.5,1.5);
  fHistEtaDeuteronRecTrackCuts->GetXaxis()->SetTitle("#it{#eta}");
  fHistEtaDeuteronRecTrackCuts->GetYaxis()->SetTitle("entries");

  fHistEtaDeuteronRecPairPtCut = new TH1F("fHistEtaDeuteronRecPairPtCut","reconstructed d (PDG + primary + track cuts + pair #it{p}_{T} cut)",200,-1.5,1.5);
  fHistEtaDeuteronRecPairPtCut->GetXaxis()->SetTitle("#it{#eta}");
  fHistEtaDeuteronRecPairPtCut->GetYaxis()->SetTitle("entries");

  // reconstructed antideuterons
  fHistPtAntiDeuteronRecPDG = new TH1F("fHistPtAntiDeuteronRecPDG","reconstructed #bar{d} (PDG)",240,0.0,6.0);
  fHistPtAntiDeuteronRecPDG->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  fHistPtAntiDeuteronRecPDG->GetYaxis()->SetTitle("entries");

  fHistPtAntiDeuteronRecPrimary = new TH1F("fHistPtAntiDeuteronRecPrimary","reconstructed #bar{d} (PDG + primary)",240,0.0,6.0);
  fHistPtAntiDeuteronRecPrimary->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  fHistPtAntiDeuteronRecPrimary->GetYaxis()->SetTitle("entries");

  fHistPtAntiDeuteronRecTrackCuts = new TH1F("fHistPtAntiDeuteronRecTrackCuts","reconstructed #bar{d} (PDG + primary + track cuts)",240,0.0,6.0);
  fHistPtAntiDeuteronRecTrackCuts->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  fHistPtAntiDeuteronRecTrackCuts->GetYaxis()->SetTitle("entries");

  fHistPtAntiDeuteronRecPairPtCut = new TH1F("fHistPtAntiDeuteronRecPairPtCut","reconstructed #bar{d} (PDG + primary + track cuts  + pair #it{p}_{T} cut)",240,0.0,6.0);
  fHistPtAntiDeuteronRecPairPtCut->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  fHistPtAntiDeuteronRecPairPtCut->GetYaxis()->SetTitle("entries");

  fHistEtaAntiDeuteronRecPDG = new TH1F("fHistEtaAntiDeuteronRecPDG","reconstruced #bar{d} (PDG)",200,-1.5,1.5);
  fHistEtaAntiDeuteronRecPDG->GetXaxis()->SetTitle("#it{#eta}");
  fHistEtaAntiDeuteronRecPDG->GetYaxis()->SetTitle("entries");

  fHistEtaAntiDeuteronRecPrimary = new TH1F("fHistEtaAntiDeuteronRecPrimary","reconstruced #bar{d} (PDG + primary)",200,-1.5,1.5);
  fHistEtaAntiDeuteronRecPrimary->GetXaxis()->SetTitle("#it{#eta}");
  fHistEtaAntiDeuteronRecPrimary->GetYaxis()->SetTitle("entries");

  fHistEtaAntiDeuteronRecTrackCuts = new TH1F("fHistEtaAntiDeuteronRecTrackCuts","reconstruced #bar{d} (PDG + primary + track cuts)",200,-1.5,1.5);
  fHistEtaAntiDeuteronRecTrackCuts->GetXaxis()->SetTitle("#it{#eta}");
  fHistEtaAntiDeuteronRecTrackCuts->GetYaxis()->SetTitle("entries");

  fHistEtaAntiDeuteronRecPairPtCut = new TH1F("fHistEtaAntiDeuteronRecPairPtCut","reconstruced #bar{d} (PDG + primary + track cuts + pair #it{p}_{T} cut)",200,-1.5,1.5);
  fHistEtaAntiDeuteronRecPairPtCut->GetXaxis()->SetTitle("#it{#eta}");
  fHistEtaAntiDeuteronRecPairPtCut->GetYaxis()->SetTitle("entries");

  // reconstructed pairs
  fHistPtHelium3RecPairPtCut = new TH1F("fHistPtHelium3RecPairPtCut","reconstructed pairs (pair #it{p}_{T} cut)",240,0.0,6.0);
  fHistPtHelium3RecPairPtCut->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  fHistPtHelium3RecPairPtCut->GetYaxis()->SetTitle("entries");

  fHistEtaHelium3RecPairPtCut = new TH1F("fHistEtaHelium3RecPairPtCut","reconstructed pairs (pair #it{p}_{T} cut)",200,-1.5,1.5);
  fHistEtaHelium3RecPairPtCut->GetXaxis()->SetTitle("#it{#eta}");
  fHistEtaHelium3RecPairPtCut->GetYaxis()->SetTitle("entries");

  fHistPtParticlesRec = new TH2F("fHistPtParticlesRec","#it{p}_{T} distribution of reconstructed p and d",240,0.0,6.0,240,0.0,6.0);
  fHistPtParticlesRec->GetXaxis()->SetTitle("#it{p}_{T} of deuterons (GeV/#it{c})");
  fHistPtParticlesRec->GetYaxis()->SetTitle("#it{p}_{T} of protons (GeV/#it{c})");

  fHistSEDPairRec = new TH1F("fHistSEDPairRec","same-event distribution of reconstructed d-p pairs",750,0.0,3.0);
  fHistSEDPairRec->GetXaxis()->SetTitle("#it{k*} (GeV/#it{c})");
  fHistSEDPairRec->GetYaxis()->SetTitle("entries");

  // reconstructed antipairs
  fHistPtAntiHelium3RecPairPtCut = new TH1F("fHistPtAntiHelium3RecPairPtCut","reconstructed antipairs (pair #it{p}_{T} cut)",240,0.0,6.0);
  fHistPtAntiHelium3RecPairPtCut->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  fHistPtAntiHelium3RecPairPtCut->GetYaxis()->SetTitle("entries");

  fHistEtaAntiHelium3RecPairPtCut = new TH1F("fHistEtaAntiHelium3RecPairPtCut","reconstructed antipairs (pair #it{p}_{T} cut)",200,-1.5,1.5);
  fHistEtaAntiHelium3RecPairPtCut->GetXaxis()->SetTitle("#it{#eta}");
  fHistEtaAntiHelium3RecPairPtCut->GetYaxis()->SetTitle("entries");

  fHistPtAntiParticlesRec = new TH2F("fHistPtAntiParticlesRec","#it{p}_{T} distribution of reconstructed antiprotons and antideuterons",240,0.0,6.0,240,0.0,6.0);
  fHistPtAntiParticlesRec->GetXaxis()->SetTitle("#it{p}_{T} of antideuterons (GeV/#it{c})");
  fHistPtAntiParticlesRec->GetYaxis()->SetTitle("#it{p}_{T} of antiprotons (GeV/#it{c})");

  fHistSEDAntiPairRec = new TH1F("fHistSEDAntiPairRec","same-event distribution of reconstructed #bar{d}-#bar{p} pairs",750,0.0,3.0);
  fHistSEDAntiPairRec->GetXaxis()->SetTitle("#it{k*} (GeV/#it{c})");
  fHistSEDAntiPairRec->GetYaxis()->SetTitle("entries");

  // histogram to count the number of events
  fHistEventCounter = new TH1F("fHistEventCounter","simple event counter",1,0,1);
  fHistEventCounter->GetXaxis()->SetBinLabel(1,"number of events");
  fHistEventCounter->GetYaxis()->SetTitle("entries");

  // generated protons
  fHistList->Add(fHistPtProtonGenPDG);
  fHistList->Add(fHistPtProtonGenPrimary);
  fHistList->Add(fHistPtProtonGenEtaCut);
  fHistList->Add(fHistPtProtonGenPairPtCut);
  fHistList->Add(fHistEtaProtonGenPDG);
  fHistList->Add(fHistEtaProtonGenPrimary);
  fHistList->Add(fHistEtaProtonGenEtaCut);
  fHistList->Add(fHistEtaProtonGenPairPtCut);

  // generated antiprotons
  fHistList->Add(fHistPtAntiProtonGenPDG);
  fHistList->Add(fHistPtAntiProtonGenPrimary);
  fHistList->Add(fHistPtAntiProtonGenEtaCut);
  fHistList->Add(fHistPtAntiProtonGenPairPtCut);
  fHistList->Add(fHistEtaAntiProtonGenPDG);
  fHistList->Add(fHistEtaAntiProtonGenPrimary);
  fHistList->Add(fHistEtaAntiProtonGenEtaCut);
  fHistList->Add(fHistEtaAntiProtonGenPairPtCut);

  // generated deuterons
  fHistList->Add(fHistPtDeuteronGenPDG);
  fHistList->Add(fHistPtDeuteronGenPrimary);
  fHistList->Add(fHistPtDeuteronGenEtaCut);
  fHistList->Add(fHistPtDeuteronGenPairPtCut);
  fHistList->Add(fHistEtaDeuteronGenPDG);
  fHistList->Add(fHistEtaDeuteronGenPrimary);
  fHistList->Add(fHistEtaDeuteronGenEtaCut);
  fHistList->Add(fHistEtaDeuteronGenPairPtCut);

  // generated antideuterons
  fHistList->Add(fHistPtAntiDeuteronGenPDG);
  fHistList->Add(fHistPtAntiDeuteronGenPrimary);
  fHistList->Add(fHistPtAntiDeuteronGenEtaCut);
  fHistList->Add(fHistPtAntiDeuteronGenPairPtCut);
  fHistList->Add(fHistEtaAntiDeuteronGenPDG);
  fHistList->Add(fHistEtaAntiDeuteronGenPrimary);
  fHistList->Add(fHistEtaAntiDeuteronGenEtaCut);
  fHistList->Add(fHistEtaAntiDeuteronGenPairPtCut);

  // generated pairs
  fHistList->Add(fHistPtHelium3GenPairPtCut);
  fHistList->Add(fHistEtaHelium3GenPairPtCut);
  fHistList->Add(fHistPtParticlesGen);
  fHistList->Add(fHistSEDPairGen);

  // generated antipairs
  fHistList->Add(fHistPtAntiHelium3GenPairPtCut);
  fHistList->Add(fHistEtaAntiHelium3GenPairPtCut);
  fHistList->Add(fHistPtAntiParticlesGen);
  fHistList->Add(fHistSEDAntiPairGen);

  // reconstructed protons
  fHistList->Add(fHistPtProtonRecPDG);
  fHistList->Add(fHistPtProtonRecPrimary);
  fHistList->Add(fHistPtProtonRecTrackCuts);
  fHistList->Add(fHistPtProtonRecPairPtCut);
  fHistList->Add(fHistEtaProtonRecPDG);
  fHistList->Add(fHistEtaProtonRecPrimary);
  fHistList->Add(fHistEtaProtonRecTrackCuts);
  fHistList->Add(fHistEtaProtonRecPairPtCut);

  // reconstructed antiprotons
  fHistList->Add(fHistPtAntiProtonRecPDG);
  fHistList->Add(fHistPtAntiProtonRecPrimary);
  fHistList->Add(fHistPtAntiProtonRecTrackCuts);
  fHistList->Add(fHistPtAntiProtonRecPairPtCut);
  fHistList->Add(fHistEtaAntiProtonRecPDG);
  fHistList->Add(fHistEtaAntiProtonRecPrimary);
  fHistList->Add(fHistEtaAntiProtonRecTrackCuts);
  fHistList->Add(fHistEtaAntiProtonRecPairPtCut);

  // reconstructed deuterons
  fHistList->Add(fHistPtDeuteronRecPDG);
  fHistList->Add(fHistPtDeuteronRecPrimary);
  fHistList->Add(fHistPtDeuteronRecTrackCuts);
  fHistList->Add(fHistPtDeuteronRecPairPtCut);
  fHistList->Add(fHistEtaDeuteronRecPDG);
  fHistList->Add(fHistEtaDeuteronRecPrimary);
  fHistList->Add(fHistEtaDeuteronRecTrackCuts);
  fHistList->Add(fHistEtaDeuteronRecPairPtCut);

  // reconstructed antideuterons
  fHistList->Add(fHistPtAntiDeuteronRecPDG);
  fHistList->Add(fHistPtAntiDeuteronRecPrimary);
  fHistList->Add(fHistPtAntiDeuteronRecTrackCuts);
  fHistList->Add(fHistPtAntiDeuteronRecPairPtCut);
  fHistList->Add(fHistEtaAntiDeuteronRecPDG);
  fHistList->Add(fHistEtaAntiDeuteronRecPrimary);
  fHistList->Add(fHistEtaAntiDeuteronRecTrackCuts);
  fHistList->Add(fHistEtaAntiDeuteronRecPairPtCut);

  // reconstructed pairs
  fHistList->Add(fHistPtHelium3RecPairPtCut);
  fHistList->Add(fHistEtaHelium3RecPairPtCut);
  fHistList->Add(fHistPtParticlesRec);
  fHistList->Add(fHistSEDPairRec);

  // reconstructed antipairs
  fHistList->Add(fHistPtAntiHelium3RecPairPtCut);
  fHistList->Add(fHistEtaAntiHelium3RecPairPtCut);
  fHistList->Add(fHistPtAntiParticlesRec);
  fHistList->Add(fHistSEDAntiPairRec);

  // event counter
  fHistList->Add(fHistEventCounter);


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
  double ProtonpTPCThreshold = 0.7;
  double DeuteronpTPCThreshold = 1.4;

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
    for(nStackLoop1 = 0; nStackLoop1 < stack->GetNtrack(); nStackLoop1++) // proton loop
      {

	AliMCParticle *Particle1 = (AliMCParticle*)(mcEvent->GetTrack(nStackLoop1));
	// check if particle is a proton
	if(!(Particle1->PdgCode() == ProtonPDG)) continue;

	MomentumProtonGen[0] = Particle1->Px();
	MomentumProtonGen[1] = Particle1->Py();
	MomentumProtonGen[2] = Particle1->Pz();
	  
	TLorentzVector LorentzVectorProtonGen;
	LorentzVectorProtonGen.SetXYZM(MomentumProtonGen[0],MomentumProtonGen[1],MomentumProtonGen[2],ProtonMass);
	PtProtonGen = LorentzVectorProtonGen.Pt();
	EtaProtonGen = LorentzVectorProtonGen.Eta();

	fHistPtProtonGenPDG->Fill(PtProtonGen);
	fHistEtaProtonGenPDG->Fill(EtaProtonGen);

	// check if proton is a primary particle
	if(!Particle1->IsPhysicalPrimary()) continue;

	fHistPtProtonGenPrimary->Fill(PtProtonGen);
	fHistEtaProtonGenPrimary->Fill(EtaProtonGen);

	// apply a cut in pseudo-rapidity
	if((EtaProtonGen < EtaLimit1) || (EtaProtonGen > EtaLimit2)) continue;

	fHistPtProtonGenEtaCut->Fill(PtProtonGen);
	fHistEtaProtonGenEtaCut->Fill(EtaProtonGen);


	    for(nStackLoop2 = 0; nStackLoop2 < stack->GetNtrack(); nStackLoop2++) // deuteron loop
	      {

		AliMCParticle *Particle2 = (AliMCParticle*)(mcEvent->GetTrack(nStackLoop2));
		// check if particle is a deuteron
		if(!(Particle2->PdgCode() == DeuteronPDG)) continue;

		MomentumDeuteronGen[0] = Particle2->Px();
		MomentumDeuteronGen[1] = Particle2->Py();
		MomentumDeuteronGen[2] = Particle2->Pz();

		TLorentzVector LorentzVectorDeuteronGen;
		LorentzVectorDeuteronGen.SetXYZM(MomentumDeuteronGen[0],MomentumDeuteronGen[1],MomentumDeuteronGen[2],DeuteronMass);
		PtDeuteronGen = LorentzVectorDeuteronGen.Pt();
		EtaDeuteronGen = LorentzVectorDeuteronGen.Eta();

		fHistPtDeuteronGenPDG->Fill(PtDeuteronGen);
		fHistEtaDeuteronGenPDG->Fill(EtaDeuteronGen);

		// check if deuteron is a primary particle
		if(!Particle2->IsPhysicalPrimary()) continue;

		fHistPtDeuteronGenPrimary->Fill(PtDeuteronGen);
		fHistEtaDeuteronGenPrimary->Fill(EtaDeuteronGen);

		// apply a cut in pseudo-rapidity
		if((EtaDeuteronGen < EtaLimit1)	|| (EtaDeuteronGen > EtaLimit2))  continue;

		fHistPtDeuteronGenEtaCut->Fill(PtDeuteronGen);
		fHistEtaDeuteronGenEtaCut->Fill(EtaDeuteronGen);

      
		TLorentzVector LorentzVectorHelium3Gen;
		LorentzVectorHelium3Gen = LorentzVectorProtonGen + LorentzVectorDeuteronGen;

		PtHelium3Gen = LorentzVectorHelium3Gen.Pt();
		EtaHelium3Gen = LorentzVectorHelium3Gen.Eta();

		// pair cut
		if((PtHelium3Gen < PairPtLimit1)  || (PtHelium3Gen > PairPtLimit2)) continue;

		RelativeMomentum = CalculateRelativeMomentum(LorentzVectorHelium3Gen,LorentzVectorProtonGen,LorentzVectorDeuteronGen);

		fHistPtProtonGenPairPtCut->Fill(PtProtonGen);
		fHistEtaProtonGenPairPtCut->Fill(EtaProtonGen);
		fHistPtDeuteronGenPairPtCut->Fill(PtDeuteronGen);
		fHistEtaDeuteronGenPairPtCut->Fill(EtaDeuteronGen);
		fHistPtHelium3GenPairPtCut->Fill(PtHelium3Gen);
		fHistEtaHelium3GenPairPtCut->Fill(EtaHelium3Gen);
		fHistSEDPairGen->Fill(RelativeMomentum);
		fHistPtParticlesGen->Fill(PtDeuteronGen,PtProtonGen);

	      } // end of stack loop 2 (deuteron loop)

      } // end of stack loop 1 (proton loop)


    // loop for the generated antiparticles
    for(nStackLoop1 = 0; nStackLoop1 < stack->GetNtrack(); nStackLoop1++) // antiproton loop
      {

	AliMCParticle *Particle1 = (AliMCParticle*)(mcEvent->GetTrack(nStackLoop1));
	// check if particle is an antiproton
	if(!(Particle1->PdgCode() == -ProtonPDG)) continue;

	MomentumProtonGen[0] = Particle1->Px();
	MomentumProtonGen[1] = Particle1->Py();
	MomentumProtonGen[2] = Particle1->Pz();
	  
	TLorentzVector LorentzVectorProtonGen;
	LorentzVectorProtonGen.SetXYZM(MomentumProtonGen[0],MomentumProtonGen[1],MomentumProtonGen[2],ProtonMass);
	PtProtonGen = LorentzVectorProtonGen.Pt();
	EtaProtonGen = LorentzVectorProtonGen.Eta();

	fHistPtAntiProtonGenPDG->Fill(PtProtonGen);
	fHistEtaAntiProtonGenPDG->Fill(EtaProtonGen);

	// check if antiproton is a primary particle
	if(!Particle1->IsPhysicalPrimary()) continue;

	fHistPtAntiProtonGenPrimary->Fill(PtProtonGen);
	fHistEtaAntiProtonGenPrimary->Fill(EtaProtonGen);

	// apply a cut in pseudo-rapidity
	if((EtaProtonGen < EtaLimit1) || (EtaProtonGen > EtaLimit2)) continue;

	fHistPtAntiProtonGenEtaCut->Fill(PtProtonGen);
	fHistEtaAntiProtonGenEtaCut->Fill(EtaProtonGen);


	    for(nStackLoop2 = 0; nStackLoop2 < stack->GetNtrack(); nStackLoop2++) // antideuteron loop
	      {

		AliMCParticle *Particle2 = (AliMCParticle*)(mcEvent->GetTrack(nStackLoop2));
		// check if particle is an anitdeuteron
		if(!(Particle2->PdgCode() == -DeuteronPDG)) continue;

		MomentumDeuteronGen[0] = Particle2->Px();
		MomentumDeuteronGen[1] = Particle2->Py();
		MomentumDeuteronGen[2] = Particle2->Pz();

		TLorentzVector LorentzVectorDeuteronGen;
		LorentzVectorDeuteronGen.SetXYZM(MomentumDeuteronGen[0],MomentumDeuteronGen[1],MomentumDeuteronGen[2],DeuteronMass);
		PtDeuteronGen = LorentzVectorDeuteronGen.Pt();
		EtaDeuteronGen = LorentzVectorDeuteronGen.Eta();

		fHistPtAntiDeuteronGenPDG->Fill(PtDeuteronGen);
		fHistEtaAntiDeuteronGenPDG->Fill(EtaDeuteronGen);

		// check if antideuteron is a primary particle
		if(!Particle2->IsPhysicalPrimary()) continue;

		fHistPtAntiDeuteronGenPrimary->Fill(PtDeuteronGen);
		fHistEtaAntiDeuteronGenPrimary->Fill(EtaDeuteronGen);

		// apply a cut in pseudo-rapidity
		if((EtaDeuteronGen < EtaLimit1)	|| (EtaDeuteronGen > EtaLimit2))  continue;

		fHistPtAntiDeuteronGenEtaCut->Fill(PtDeuteronGen);
		fHistEtaAntiDeuteronGenEtaCut->Fill(EtaDeuteronGen);

      
		TLorentzVector LorentzVectorHelium3Gen;
		LorentzVectorHelium3Gen = LorentzVectorProtonGen + LorentzVectorDeuteronGen;

		PtHelium3Gen = LorentzVectorHelium3Gen.Pt();
		EtaHelium3Gen = LorentzVectorHelium3Gen.Eta();

		// antipair cut
		if((PtHelium3Gen < PairPtLimit1)  || (PtHelium3Gen > PairPtLimit2)) continue;

		RelativeMomentum = CalculateRelativeMomentum(LorentzVectorHelium3Gen,LorentzVectorProtonGen,LorentzVectorDeuteronGen);

		fHistPtAntiProtonGenPairPtCut->Fill(PtProtonGen);
		fHistEtaAntiProtonGenPairPtCut->Fill(EtaProtonGen);
		fHistPtAntiDeuteronGenPairPtCut->Fill(PtDeuteronGen);
		fHistEtaAntiDeuteronGenPairPtCut->Fill(EtaDeuteronGen);
		fHistPtAntiHelium3GenPairPtCut->Fill(PtHelium3Gen);
	  	fHistEtaAntiHelium3GenPairPtCut->Fill(EtaHelium3Gen);
		fHistSEDAntiPairGen->Fill(RelativeMomentum);
		fHistPtAntiParticlesGen->Fill(PtDeuteronGen,PtProtonGen);

	      } // end of stack loop 2 (antideuteron loop)

      } // end of stack loop 1 (antiproton loop)



    // loop for reconstructed particles
    for(int nTrack1 = 0; nTrack1 < fESD->GetNumberOfTracks(); nTrack1++) // proton loop
      {

	TPC1_OK = false;
	TOF1_OK = false;
	TPC2_OK = false;
	TOF2_OK = false;

	AliESDtrack *Track1 = dynamic_cast<AliESDtrack*>(fESD->GetTrack(nTrack1));

	AliMCParticle *Particle1 = (AliMCParticle*)(mcEvent->GetTrack(TMath::Abs(Track1->GetLabel())));
	int Particle1PDG = Particle1->PdgCode();

	MomentumProtonRec[0] = Track1->Px();
	MomentumProtonRec[1] = Track1->Py();
	MomentumProtonRec[2] = Track1->Pz();

	TLorentzVector LorentzVectorProtonRec;
	LorentzVectorProtonRec.SetXYZM(MomentumProtonRec[0],MomentumProtonRec[1],MomentumProtonRec[2],ProtonMass);

	PtProtonRec = LorentzVectorProtonRec.Pt();
	EtaProtonRec = LorentzVectorProtonRec.Eta();

	// check if particle is a proton
	if(!(Particle1PDG == ProtonPDG))	continue;

	fHistPtProtonRecPDG->Fill(PtProtonRec);
	fHistEtaProtonRecPDG->Fill(EtaProtonRec);

	// check if proton is a primary particle
	if(!Particle1->IsPhysicalPrimary()) continue;

	fHistPtProtonRecPrimary->Fill(PtProtonRec);
	fHistEtaProtonRecPrimary->Fill(EtaProtonRec);

	// apply track cuts
	if(!fESDtrackCutsProton->AcceptTrack(Track1)) continue;
	if(Track1->GetSign() < 1) continue;

	AliPIDResponse::EDetPidStatus status1TPC = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTPC,Track1);
	AliPIDResponse::EDetPidStatus status1TOF = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTOF,Track1);

	if(status1TPC == AliPIDResponse::kDetPidOk) TPC1_OK = true;
	if(status1TOF == AliPIDResponse::kDetPidOk) TOF1_OK = true;
	if(!(TPC1_OK)) continue;
	if(!(TOF1_OK)) continue;

	nSigmaTPCproton = fPIDResponse->NumberOfSigmasTPC(Track1,AliPID::kProton);
	nSigmaTOFproton = fPIDResponse->NumberOfSigmasTOF(Track1,AliPID::kProton);
	nSigmaTPCTOFproton = TMath::Sqrt(pow(nSigmaTPCproton,2.) + pow(nSigmaTOFproton,2.));

	Double_t ProtonpTPC = Track1->GetInnerParam()->GetP();

	if(ProtonpTPC <= ProtonpTPCThreshold){

	  if(!(nSigmaTPCproton < 3.0)) continue;

	}else{
	    
	  if(!(nSigmaTPCTOFproton < 3.0)) continue;

        }

	fHistPtProtonRecTrackCuts->Fill(PtProtonRec);
	fHistEtaProtonRecTrackCuts->Fill(EtaProtonRec);



	for(int nTrack2 = nTrack1 + 1; nTrack2 < fESD->GetNumberOfTracks(); nTrack2++)	// deuteron loop
	  {

	    AliESDtrack *Track2 = dynamic_cast<AliESDtrack*>(fESD->GetTrack(nTrack2));

	    AliMCParticle *Particle2 = (AliMCParticle*)(mcEvent->GetTrack(TMath::Abs(Track2->GetLabel())));
	    int Particle2PDG = Particle2->PdgCode();

	    MomentumDeuteronRec[0] = Track2->Px();
	    MomentumDeuteronRec[1] = Track2->Py();
	    MomentumDeuteronRec[2] = Track2->Pz();

	    TLorentzVector LorentzVectorDeuteronRec;
	    LorentzVectorDeuteronRec.SetXYZM(MomentumDeuteronRec[0],MomentumDeuteronRec[1],MomentumDeuteronRec[2],DeuteronMass);

	    PtDeuteronRec = LorentzVectorDeuteronRec.Pt();
	    EtaDeuteronRec = LorentzVectorDeuteronRec.Eta();

	    // check if particle is a deuteron
	    if(!(Particle2PDG == DeuteronPDG))	continue;

	    fHistPtDeuteronRecPDG->Fill(PtDeuteronRec);
	    fHistEtaDeuteronRecPDG->Fill(EtaDeuteronRec);

	    // check if deuteron is a primary particle
	    if((!Particle2->IsPhysicalPrimary())) continue;

	    fHistPtDeuteronRecPrimary->Fill(PtDeuteronRec);
	    fHistEtaDeuteronRecPrimary->Fill(EtaDeuteronRec);

	    // apply track cuts
	    if(!fESDtrackCutsDeuteron->AcceptTrack(Track2)) continue;
	    if(Track2->GetSign() < 1) continue;

	    AliPIDResponse::EDetPidStatus status2TPC = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTPC,Track2);
	    AliPIDResponse::EDetPidStatus status2TOF = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTOF,Track2);

	    if(status2TPC == AliPIDResponse::kDetPidOk) TPC2_OK = true;
	    if(status2TOF == AliPIDResponse::kDetPidOk) TOF2_OK = true;

	    if(!(TPC2_OK)) continue;
	    if(!(TOF2_OK)) continue;

	    nSigmaTPCdeuteron = fPIDResponse->NumberOfSigmasTPC(Track2,AliPID::kDeuteron);
	    nSigmaTOFdeuteron = fPIDResponse->NumberOfSigmasTOF(Track2,AliPID::kDeuteron);
	    nSigmaTPCTOFdeuteron = TMath::Sqrt(pow(nSigmaTPCdeuteron,2.) + pow(nSigmaTOFdeuteron,2.));

	    Double_t DeuteronpTPC = Track2->GetInnerParam()->GetP();

	    if(DeuteronpTPC <= DeuteronpTPCThreshold){

	      if(!(nSigmaTPCdeuteron < 3.0)) continue;

	    }else{
	    
	      if(!(nSigmaTPCTOFdeuteron < 3.0)) continue;

	    }

	    fHistPtDeuteronRecTrackCuts->Fill(PtDeuteronRec);
	    fHistEtaDeuteronRecTrackCuts->Fill(EtaDeuteronRec);

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



	    TLorentzVector LorentzVectorHelium3Rec;
	    LorentzVectorHelium3Rec = LorentzVectorProtonRec + LorentzVectorDeuteronRec;

	    EtaHelium3Rec = LorentzVectorHelium3Rec.Eta();
	    PtHelium3Rec = LorentzVectorHelium3Rec.Pt();

	    // apply pair cut
	    if((PtHelium3Rec < PairPtLimit1)	  || (PtHelium3Rec > PairPtLimit2)) continue;

	    RelativeMomentum = CalculateRelativeMomentum(LorentzVectorHelium3Rec,LorentzVectorProtonRec,LorentzVectorDeuteronRec);

	    fHistPtProtonRecPairPtCut->Fill(PtProtonRec);
	    fHistEtaProtonRecPairPtCut->Fill(EtaProtonRec);
	    fHistPtDeuteronRecPairPtCut->Fill(PtDeuteronRec);
	    fHistEtaDeuteronRecPairPtCut->Fill(EtaDeuteronRec);
	    fHistPtHelium3RecPairPtCut->Fill(PtHelium3Rec);
	    fHistEtaHelium3RecPairPtCut->Fill(EtaHelium3Rec);
	    fHistSEDPairRec->Fill(RelativeMomentum);
	    fHistPtParticlesRec->Fill(PtDeuteronRec,PtProtonRec);

	  } // end of track2 loop (deuteron loop)

      } // end of track1 loop (proton loop)




    // loop for reconstructed antiparticles
    for(int nTrack1 = 0; nTrack1 < fESD->GetNumberOfTracks(); nTrack1++) // antiproton loop
      {

	TPC1_OK = false;
	TOF1_OK = false;
	TPC2_OK = false;
	TOF2_OK = false;

	AliESDtrack *Track1 = dynamic_cast<AliESDtrack*>(fESD->GetTrack(nTrack1));

	AliMCParticle *Particle1 = (AliMCParticle*)(mcEvent->GetTrack(TMath::Abs(Track1->GetLabel())));
	int Particle1PDG = Particle1->PdgCode();

	MomentumProtonRec[0] = Track1->Px();
	MomentumProtonRec[1] = Track1->Py();
	MomentumProtonRec[2] = Track1->Pz();

	TLorentzVector LorentzVectorProtonRec;
	LorentzVectorProtonRec.SetXYZM(MomentumProtonRec[0],MomentumProtonRec[1],MomentumProtonRec[2],ProtonMass);

	PtProtonRec = LorentzVectorProtonRec.Pt();
	EtaProtonRec = LorentzVectorProtonRec.Eta();

	// check if particle is an antiproton
	if(!(Particle1PDG == -ProtonPDG))	continue;

	fHistPtAntiProtonRecPDG->Fill(PtProtonRec);
	fHistEtaAntiProtonRecPDG->Fill(EtaProtonRec);

	// check if antiproton is a primary particle
	if(!Particle1->IsPhysicalPrimary()) continue;

	fHistPtAntiProtonRecPrimary->Fill(PtProtonRec);
	fHistEtaAntiProtonRecPrimary->Fill(EtaProtonRec);

	// apply track cuts
	if(!fESDtrackCutsProton->AcceptTrack(Track1)) continue;
	if(Track1->GetSign() > -1) continue;

	AliPIDResponse::EDetPidStatus status1TPC = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTPC,Track1);
	AliPIDResponse::EDetPidStatus status1TOF = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTOF,Track1);

	if(status1TPC == AliPIDResponse::kDetPidOk) TPC1_OK = true;
	if(status1TOF == AliPIDResponse::kDetPidOk) TOF1_OK = true;
	if(!(TPC1_OK)) continue;
	if(!(TOF1_OK)) continue;

	nSigmaTPCproton = fPIDResponse->NumberOfSigmasTPC(Track1,AliPID::kProton);
	nSigmaTOFproton = fPIDResponse->NumberOfSigmasTOF(Track1,AliPID::kProton);
	nSigmaTPCTOFproton = TMath::Sqrt(pow(nSigmaTPCproton,2.) + pow(nSigmaTOFproton,2.));

	Double_t ProtonpTPC = Track1->GetInnerParam()->GetP();

	if(ProtonpTPC <= ProtonpTPCThreshold){

	  if(!(nSigmaTPCproton < 3.0)) continue;

	}else{
	    
	  if(!(nSigmaTPCTOFproton < 3.0)) continue;

        }

	fHistPtAntiProtonRecTrackCuts->Fill(PtProtonRec);
	fHistEtaAntiProtonRecTrackCuts->Fill(EtaProtonRec);



	for(int nTrack2 = nTrack1 + 1; nTrack2 < fESD->GetNumberOfTracks(); nTrack2++)	// antideuteron loop
	  {

	    AliESDtrack *Track2 = dynamic_cast<AliESDtrack*>(fESD->GetTrack(nTrack2));

	    AliMCParticle *Particle2 = (AliMCParticle*)(mcEvent->GetTrack(TMath::Abs(Track2->GetLabel())));
	    int Particle2PDG = Particle2->PdgCode();

	    MomentumDeuteronRec[0] = Track2->Px();
	    MomentumDeuteronRec[1] = Track2->Py();
	    MomentumDeuteronRec[2] = Track2->Pz();

	    TLorentzVector LorentzVectorDeuteronRec;
	    LorentzVectorDeuteronRec.SetXYZM(MomentumDeuteronRec[0],MomentumDeuteronRec[1],MomentumDeuteronRec[2],DeuteronMass);

	    PtDeuteronRec = LorentzVectorDeuteronRec.Pt();
	    EtaDeuteronRec = LorentzVectorDeuteronRec.Eta();

	    // check if particle is an antideuteron
	    if(!(Particle2PDG == -DeuteronPDG))	continue;

	    fHistPtAntiDeuteronRecPDG->Fill(PtDeuteronRec);
	    fHistEtaAntiDeuteronRecPDG->Fill(EtaDeuteronRec);

	    // check if antideuteron is a primary particle
	    if((!Particle2->IsPhysicalPrimary())) continue;

	    fHistPtAntiDeuteronRecPrimary->Fill(PtDeuteronRec);
	    fHistEtaAntiDeuteronRecPrimary->Fill(EtaDeuteronRec);

	    // apply track cuts
	    if(!fESDtrackCutsDeuteron->AcceptTrack(Track2)) continue;
	    if(Track2->GetSign() > -1) continue;

	    AliPIDResponse::EDetPidStatus status2TPC = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTPC,Track2);
	    AliPIDResponse::EDetPidStatus status2TOF = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTOF,Track2);

	    if(status2TPC == AliPIDResponse::kDetPidOk) TPC2_OK = true;
	    if(status2TOF == AliPIDResponse::kDetPidOk) TOF2_OK = true;

	    if(!(TPC2_OK)) continue;
	    if(!(TOF2_OK)) continue;

	    nSigmaTPCdeuteron = fPIDResponse->NumberOfSigmasTPC(Track2,AliPID::kDeuteron);
	    nSigmaTOFdeuteron = fPIDResponse->NumberOfSigmasTOF(Track2,AliPID::kDeuteron);
	    nSigmaTPCTOFdeuteron = TMath::Sqrt(pow(nSigmaTPCdeuteron,2.) + pow(nSigmaTOFdeuteron,2.));

	    Double_t DeuteronpTPC = Track2->GetInnerParam()->GetP();

	    if(DeuteronpTPC <= DeuteronpTPCThreshold){

	      if(!(nSigmaTPCdeuteron < 3.0)) continue;

	    }else{
	    
	      if(!(nSigmaTPCTOFdeuteron < 3.0)) continue;

	    }

	    fHistPtAntiDeuteronRecTrackCuts->Fill(PtDeuteronRec);
	    fHistEtaAntiDeuteronRecTrackCuts->Fill(EtaDeuteronRec);

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



	    TLorentzVector LorentzVectorHelium3Rec;
	    LorentzVectorHelium3Rec = LorentzVectorProtonRec + LorentzVectorDeuteronRec;

	    EtaHelium3Rec = LorentzVectorHelium3Rec.Eta();
	    PtHelium3Rec = LorentzVectorHelium3Rec.Pt();

	    // apply antipair cut
	    if((PtHelium3Rec < PairPtLimit1)	  || (PtHelium3Rec > PairPtLimit2)) continue;

	    RelativeMomentum = CalculateRelativeMomentum(LorentzVectorHelium3Rec,LorentzVectorProtonRec,LorentzVectorDeuteronRec);

	    fHistPtAntiProtonRecPairPtCut->Fill(PtProtonRec);
	    fHistEtaAntiProtonRecPairPtCut->Fill(EtaProtonRec);
	    fHistPtAntiDeuteronRecPairPtCut->Fill(PtDeuteronRec);
	    fHistEtaAntiDeuteronRecPairPtCut->Fill(EtaDeuteronRec);
	    fHistPtAntiHelium3RecPairPtCut->Fill(PtHelium3Rec);
	    fHistEtaAntiHelium3RecPairPtCut->Fill(EtaHelium3Rec);
	    fHistSEDAntiPairRec->Fill(RelativeMomentum);
	    fHistPtAntiParticlesRec->Fill(PtDeuteronRec,PtProtonRec);

	  } // end of track2 loop (antideuteron loop)

      } // end of track1 loop (antiproton loop)






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


