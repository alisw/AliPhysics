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
  fHistPtProtonRecTrackCuts(0),
  fHistPtProtonRecPairPtCut(0),
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
  fHistdEdx_LHC18a2a(0),
  fHistdEdx_LHC18a2b(0),
  fHistdEdx_LHC18a2b4(0),
  fHistdEdx_LHC20l7a(0),
  GeneratedProtonArray(0),
  GeneratedDeuteronArray(0),
  GeneratedAntiProtonArray(0),
  GeneratedAntiDeuteronArray(0),
  ReconstructedProtonArray(0),
  ReconstructedDeuteronArray(0),
  ReconstructedAntiProtonArray(0),
  ReconstructedAntiDeuteronArray(0),
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
  fHistdEdx_LHC18a2a(0),
  fHistdEdx_LHC18a2b(0),
  fHistdEdx_LHC18a2b4(0),
  fHistdEdx_LHC20l7a(0),
  GeneratedProtonArray(0),
  GeneratedDeuteronArray(0),
  GeneratedAntiProtonArray(0),
  GeneratedAntiDeuteronArray(0),
  ReconstructedProtonArray(0),
  ReconstructedDeuteronArray(0),
  ReconstructedAntiProtonArray(0),
  ReconstructedAntiDeuteronArray(0),
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

  // Deuteron cuts
  fESDtrackCutsDeuteron = new AliESDtrackCuts("AliESDtrackCutsDeuteron","AliESDtrackCutsDeuteron");
//  fESDtrackCutsDeuteron = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts(); // FilterBit 128
  fESDtrackCutsDeuteron->SetEtaRange(-0.8,0.8);
  fESDtrackCutsDeuteron->SetPtRange(0.4,4.0);
  fESDtrackCutsDeuteron->SetMinNClustersTPC(80);
  fESDtrackCutsDeuteron->SetMinNCrossedRowsTPC(70);
  fESDtrackCutsDeuteron->SetMinRatioCrossedRowsOverFindableClustersTPC(0.83);
  fESDtrackCutsDeuteron->SetAcceptSharedTPCClusters(false);
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

  // dE/dx for reconstructed particles and antiparticles
  fHistdEdx_LHC18a2a = new TH2F("fHistdEdx_LHC18a2a","reconstructed particles and antiparticles in LHC18a2a",1000,-5.0,5.0,1000,0.0,1500.0);
  fHistdEdx_LHC18a2a->GetXaxis()->SetTitle("#frac{#it{p}}{#it{z}} (GeV/#it{c})");
  fHistdEdx_LHC18a2a->GetYaxis()->SetTitle("TPC signal (a.u.)");

  fHistdEdx_LHC18a2b = new TH2F("fHistdEdx_LHC18a2b","reconstructed particles and antiparticles in LHC18a2b",1000,-5.0,5.0,1000,0.0,1500.0);
  fHistdEdx_LHC18a2b->GetXaxis()->SetTitle("#frac{#it{p}}{#it{z}} (GeV/#it{c})");
  fHistdEdx_LHC18a2b->GetYaxis()->SetTitle("TPC signal (a.u.)");

  fHistdEdx_LHC18a2b4 = new TH2F("fHistdEdx_LHC18a2b4","reconstructed particles and antiparticles in LHC18a2b4",1000,-5.0,5.0,1000,0.0,1500.0);
  fHistdEdx_LHC18a2b4->GetXaxis()->SetTitle("#frac{#it{p}}{#it{z}} (GeV/#it{c})");
  fHistdEdx_LHC18a2b4->GetYaxis()->SetTitle("TPC signal (a.u.)");

  fHistdEdx_LHC20l7a = new TH2F("fHistdEdx_LHC20l7a","reconstructed particles and antiparticles in LHC20l7a",1000,-5.0,5.0,1000,0.0,1500.0);
  fHistdEdx_LHC20l7a->GetXaxis()->SetTitle("#frac{#it{p}}{#it{z}} (GeV/#it{c})");
  fHistdEdx_LHC20l7a->GetYaxis()->SetTitle("TPC signal (a.u.)");

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

  // dE/dx for particles and antiparticles
  fHistList->Add(fHistdEdx_LHC18a2a);
  fHistList->Add(fHistdEdx_LHC18a2b);
  fHistList->Add(fHistdEdx_LHC18a2b4);
  fHistList->Add(fHistdEdx_LHC20l7a);

  // event counter
  fHistList->Add(fHistEventCounter);

  GeneratedProtonArray = new std::vector<int>;
  GeneratedDeuteronArray = new std::vector<int>;
  GeneratedAntiProtonArray = new std::vector<int>;
  GeneratedAntiDeuteronArray = new std::vector<int>;

  ReconstructedProtonArray = new std::vector<int>;
  ReconstructedDeuteronArray = new std::vector<int>;
  ReconstructedAntiProtonArray = new std::vector<int>;
  ReconstructedAntiDeuteronArray = new std::vector<int>;

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
  double MassProton   = 0.93827;
  double MassDeuteron = 1.87561;

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


  // Particle labels
  int LabelProton	= 0;
  int LabelDeuteron	= 0;
  int LabelAntiProton	= 0;
  int LabelAntiDeuteron	= 0;


  // Single-particle and pair cuts
  double EtaLimit1 = -0.8;
  double EtaLimit2 = +0.8;
  double PairPtLimit1 = 0.0;
  double PairPtLimit2 = 999.0;
  double ProtonpTPCThreshold = 0.7;
  double DeuteronpTPCThreshold = 1.4;

  // Variables for particles
  double MomentumProton[3]  = {0,0,0};
  double EtaProton	    = 0.0;
  double PtProton	    = 0.0;

  double MomentumDeuteron[3]  = {0,0,0};
  double EtaDeuteron	      = 0.0;
  double PtDeuteron	      = 0.0;

  double MomentumAntiProton[3]  = {0,0,0};
  double EtaAntiProton		= 0.0;
  double PtAntiProton		= 0.0;

  double MomentumAntiDeuteron[3]  = {0,0,0};
  double EtaAntiDeuteron	  = 0.0;
  double PtAntiDeuteron		  = 0.0;

  double PtPair		  = 0.0;
  double EtaPair	  = 0.0;
  



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

  double fBetheParametersDeuteron[6];

  // DCA calculation
  Float_t xv_Proton[2];
  Float_t yv_Proton[3];
  Float_t xv_Deuteron[2];
  Float_t yv_Deuteron[3];
  
  Float_t ProtonDCAxy;
  Float_t ProtonDCAz;
  Float_t DeuteronDCAxy;
  Float_t DeuteronDCAz;

  double ProtonDCAxyMax = 0.1;
  double ProtonDCAzMax	= 0.2;
  double DeuteronDCAxyMax = 0.1;
  double DeuteronDCAzMax  = 0.2;

  double RelativeMomentum = 0.0;

  int nStackLoop = 0;
  int nStackLoop1 = 0;
  int nStackLoop2 = 0;



  // Variables for PID
  double nSigmaTPCTOFproton = 0.0;
  double nSigmaTPCTOFdeuteron = 0.0;
  double nSigmaTPCproton = 0.0;
  double nSigmaTPCdeuteron = 0.0;
  double nSigmaTOFproton = 0.0;
  double nSigmaTOFdeuteron = 0.0;

  bool TPC_OK = false;
  bool TOF_OK = false;


  int RunNumber = fESD->GetRunNumber();
  int ParticleLabel = 0;

  bool LHC18a2a = false;
  bool LHC18a2b = false;
  bool LHC18a2b4 = false;
  bool LHC20l7a = false;


  if(fMCtrue){


    // set the appropriate Bethe-Bloch parameters for deuterons for the individual productions

    if(RunNumber >= 270581 && RunNumber <= 282704){
    // LHC18a2a (anchored to 2017 data)

	LHC18a2a = true;
	fBetheParametersDeuteron[0] = 2.63655; 
	fBetheParametersDeuteron[1] = 16.8483;
	fBetheParametersDeuteron[2] = 0.0854068; 
	fBetheParametersDeuteron[3] = 2.03083;
	fBetheParametersDeuteron[4] = -1.90682;
	fBetheParametersDeuteron[5] = 0.06;

    }else if(RunNumber >= 256504 && RunNumber <= 259888){
    // LHC18a2b4 (anchored to 2016 k/l pass2)

	LHC18a2b4 = true;
	fBetheParametersDeuteron[0] = 1.64162; 
	fBetheParametersDeuteron[1] = 24.2932;
	fBetheParametersDeuteron[2] = 0.464873; 
	fBetheParametersDeuteron[3] = 2.0147;
	fBetheParametersDeuteron[4] = -0.773004;
	fBetheParametersDeuteron[5] = 0.06;

    }else if(RunNumber >= 252235 && RunNumber <= 264347){
    // LHC18a2b (anchored to 2016 data)

	LHC18a2b = true;
	fBetheParametersDeuteron[0] = 3.74652; 
	fBetheParametersDeuteron[1] = 12.159;
	fBetheParametersDeuteron[2] = 0.0332405; 
	fBetheParametersDeuteron[3] = 1.77673;
	fBetheParametersDeuteron[4] = -2.31179;
	fBetheParametersDeuteron[5] = 0.06;

    }else if(RunNumber >= 285009 && RunNumber <= 294925){
    // LHC20l7a (anchored to 2018 data)

	LHC20l7a = true;
	// parameters taken from LHC18a2a
	fBetheParametersDeuteron[0] = 2.63655; 
	fBetheParametersDeuteron[1] = 16.8483;
	fBetheParametersDeuteron[2] = 0.0854068; 
	fBetheParametersDeuteron[3] = 2.03083;
	fBetheParametersDeuteron[4] = -1.90682;
	fBetheParametersDeuteron[5] = 0.06;

    }

    GeneratedProtonArray->clear();
    GeneratedDeuteronArray->clear();
    GeneratedAntiProtonArray->clear();
    GeneratedAntiDeuteronArray->clear();

    ReconstructedProtonArray->clear();
    ReconstructedDeuteronArray->clear();
    ReconstructedAntiProtonArray->clear();
    ReconstructedAntiDeuteronArray->clear();




// +++++++++++++++++++++++++++++++++++++++++++++++++++++
// +++++++++++ Loops for generated particles +++++++++++
// +++++++++++++++++++++++++++++++++++++++++++++++++++++

    // loop for generated protons
    for(int track = 0; track < stack->GetNtrack(); track++)
      {

	AliMCParticle *Proton = (AliMCParticle*)(mcEvent->GetTrack(track));
	LabelProton = TMath::Abs(Proton->GetLabel());

	MomentumProton[0] = Proton->Px();
	MomentumProton[1] = Proton->Py();
	MomentumProton[2] = Proton->Pz();
	  
	TLorentzVector LorentzVectorProton;
	LorentzVectorProton.SetXYZM(MomentumProton[0],MomentumProton[1],MomentumProton[2],MassProton);
	PtProton = LorentzVectorProton.Pt();
	EtaProton = LorentzVectorProton.Eta();


	// check if particle is a proton
	if(!(Proton->PdgCode() == ProtonPDG)) continue;

	    fHistPtProtonGenPDG->Fill(PtProton);
	    fHistEtaProtonGenPDG->Fill(EtaProton);


	// check if proton is a primary particle
	if(!Proton->IsPhysicalPrimary()) continue;

	    fHistPtProtonGenPrimary->Fill(PtProton);
	    fHistEtaProtonGenPrimary->Fill(EtaProton);


	// apply a cut in pseudo-rapidity
	if((EtaProton < EtaLimit1) || (EtaProton > EtaLimit2)) continue;

	  fHistPtProtonGenEtaCut->Fill(PtProton);
	  fHistEtaProtonGenEtaCut->Fill(EtaProton);

	GeneratedProtonArray->push_back(LabelProton);

      } // end of loop for generated protons



    // loop for generated deuterons
    for(int track = 0; track < stack->GetNtrack(); track++)
      {

	AliMCParticle *Deuteron = (AliMCParticle*)(mcEvent->GetTrack(track));
	LabelDeuteron = TMath::Abs(Deuteron->GetLabel());

	MomentumDeuteron[0] = Deuteron->Px();
	MomentumDeuteron[1] = Deuteron->Py();
	MomentumDeuteron[2] = Deuteron->Pz();

	TLorentzVector LorentzVectorDeuteron;
	LorentzVectorDeuteron.SetXYZM(MomentumDeuteron[0],MomentumDeuteron[1],MomentumDeuteron[2],MassDeuteron);
	PtDeuteron = LorentzVectorDeuteron.Pt();
	EtaDeuteron = LorentzVectorDeuteron.Eta();


	// check if particle is a deuteron
	if(!(Deuteron->PdgCode() == DeuteronPDG)) continue;

	    fHistPtDeuteronGenPDG->Fill(PtDeuteron);
	    fHistEtaDeuteronGenPDG->Fill(EtaDeuteron);


	// check if deuteron is a primary particle
	if(!Deuteron->IsPhysicalPrimary()) continue;

	    fHistPtDeuteronGenPrimary->Fill(PtDeuteron);
	    fHistEtaDeuteronGenPrimary->Fill(EtaDeuteron);


	// apply a cut in pseudo-rapidity
	if((EtaDeuteron < EtaLimit1)	|| (EtaDeuteron > EtaLimit2))  continue;

	    fHistPtDeuteronGenEtaCut->Fill(PtDeuteron);
	    fHistEtaDeuteronGenEtaCut->Fill(EtaDeuteron);

	GeneratedDeuteronArray->push_back(LabelDeuteron);

      }	// end of loop for generated deuterons


    // loop for generated particle pairs
    for(int track1 = 0; track1 < GeneratedProtonArray->size(); track1++) // proton loop
      {

	AliMCParticle *Proton = (AliMCParticle*)(mcEvent->GetTrack(GeneratedProtonArray->at(track1)));

	MomentumProton[0] = Proton->Px();
	MomentumProton[1] = Proton->Py();
	MomentumProton[2] = Proton->Pz();

	TLorentzVector LorentzVectorProton;
	LorentzVectorProton.SetXYZM(MomentumProton[0],MomentumProton[1],MomentumProton[2],MassProton);
	PtProton = LorentzVectorProton.Pt();
	EtaProton = LorentzVectorProton.Eta();

	for(int track2 = 0; track2 < GeneratedDeuteronArray->size(); track2++) // deuteron loop
	  {

	    AliMCParticle *Deuteron = (AliMCParticle*)(mcEvent->GetTrack(GeneratedDeuteronArray->at(track2)));

	    MomentumDeuteron[0] = Deuteron->Px();
	    MomentumDeuteron[1] = Deuteron->Py();
	    MomentumDeuteron[2] = Deuteron->Pz();

	    TLorentzVector LorentzVectorDeuteron;
	    LorentzVectorDeuteron.SetXYZM(MomentumDeuteron[0],MomentumDeuteron[1],MomentumDeuteron[2],MassDeuteron);
	    PtDeuteron = LorentzVectorDeuteron.Pt();
	    EtaDeuteron = LorentzVectorDeuteron.Eta();
      
	    TLorentzVector LorentzVectorPair;
	    LorentzVectorPair = LorentzVectorProton + LorentzVectorDeuteron;

	    PtPair = LorentzVectorPair.Pt();
	    EtaPair = LorentzVectorPair.Eta();

	    // pair cut
	    if((PtPair < PairPtLimit1)  || (PtPair > PairPtLimit2)) continue;

	    RelativeMomentum = CalculateRelativeMomentum(LorentzVectorPair,LorentzVectorProton,LorentzVectorDeuteron);

		fHistPtProtonGenPairPtCut->Fill(PtProton);
		fHistEtaProtonGenPairPtCut->Fill(EtaProton);
		fHistPtDeuteronGenPairPtCut->Fill(PtDeuteron);
		fHistEtaDeuteronGenPairPtCut->Fill(EtaDeuteron);
		fHistPtHelium3GenPairPtCut->Fill(PtPair);
		fHistEtaHelium3GenPairPtCut->Fill(EtaPair);
		fHistSEDPairGen->Fill(RelativeMomentum);
		fHistPtParticlesGen->Fill(PtDeuteron,PtProton);


	  } // end of loop for generated particle pairs (deuteron loop)

      } // end of loop for generated particle pairs (proton loop)






// +++++++++++++++++++++++++++++++++++++++++++++++++++++
// +++++++++++ Loops for generated antiparticles +++++++
// +++++++++++++++++++++++++++++++++++++++++++++++++++++

    // loop for generated antiprotons
    for(int track = 0; track < stack->GetNtrack(); track++)
      {

	AliMCParticle *AntiProton = (AliMCParticle*)(mcEvent->GetTrack(track));
	LabelAntiProton = TMath::Abs(AntiProton->GetLabel());

	MomentumAntiProton[0] = AntiProton->Px();
	MomentumAntiProton[1] = AntiProton->Py();
	MomentumAntiProton[2] = AntiProton->Pz();
	  
	TLorentzVector LorentzVectorAntiProton;
	LorentzVectorAntiProton.SetXYZM(MomentumAntiProton[0],MomentumAntiProton[1],MomentumAntiProton[2],MassProton);
	PtAntiProton = LorentzVectorAntiProton.Pt();
	EtaAntiProton = LorentzVectorAntiProton.Eta();


	// check if particle is an antiproton
	if(!(AntiProton->PdgCode() == -ProtonPDG)) continue;

	    fHistPtAntiProtonGenPDG->Fill(PtAntiProton);
	    fHistEtaAntiProtonGenPDG->Fill(EtaAntiProton);


	// check if antiproton is a primary particle
	if(!AntiProton->IsPhysicalPrimary()) continue;

	    fHistPtAntiProtonGenPrimary->Fill(PtAntiProton);
	    fHistEtaAntiProtonGenPrimary->Fill(EtaAntiProton);


	// apply a cut in pseudo-rapidity
	if((EtaAntiProton < EtaLimit1) || (EtaAntiProton > EtaLimit2)) continue;

	  fHistPtAntiProtonGenEtaCut->Fill(PtAntiProton);
	  fHistEtaAntiProtonGenEtaCut->Fill(EtaAntiProton);

	GeneratedAntiProtonArray->push_back(LabelAntiProton);

      } // end of loop for generated antiprotons



    // loop for generated antideuterons
    for(int track = 0; track < stack->GetNtrack(); track++)
      {

	AliMCParticle *AntiDeuteron = (AliMCParticle*)(mcEvent->GetTrack(track));
	LabelAntiDeuteron = TMath::Abs(AntiDeuteron->GetLabel());

	MomentumAntiDeuteron[0] = AntiDeuteron->Px();
	MomentumAntiDeuteron[1] = AntiDeuteron->Py();
	MomentumAntiDeuteron[2] = AntiDeuteron->Pz();

	TLorentzVector LorentzVectorAntiDeuteron;
	LorentzVectorAntiDeuteron.SetXYZM(MomentumAntiDeuteron[0],MomentumAntiDeuteron[1],MomentumAntiDeuteron[2],MassDeuteron);
	PtAntiDeuteron = LorentzVectorAntiDeuteron.Pt();
	EtaAntiDeuteron = LorentzVectorAntiDeuteron.Eta();


	// check if particle is an antideuteron
	if(!(AntiDeuteron->PdgCode() == -DeuteronPDG)) continue;

	    fHistPtAntiDeuteronGenPDG->Fill(PtAntiDeuteron);
	    fHistEtaAntiDeuteronGenPDG->Fill(EtaAntiDeuteron);


	// check if antideuteron is a primary particle
	if(!AntiDeuteron->IsPhysicalPrimary()) continue;

	    fHistPtAntiDeuteronGenPrimary->Fill(PtAntiDeuteron);
	    fHistEtaAntiDeuteronGenPrimary->Fill(EtaAntiDeuteron);


	// apply a cut in pseudo-rapidity
	if((EtaAntiDeuteron < EtaLimit1)	|| (EtaAntiDeuteron > EtaLimit2))  continue;

	    fHistPtAntiDeuteronGenEtaCut->Fill(PtAntiDeuteron);
	    fHistEtaAntiDeuteronGenEtaCut->Fill(EtaAntiDeuteron);

	GeneratedAntiDeuteronArray->push_back(LabelAntiDeuteron);

      }	// end of loop for generated antideuterons


    // loop for generated antiparticle pairs
    for(int track1 = 0; track1 < GeneratedAntiProtonArray->size(); track1++)
      {

	AliMCParticle *AntiProton = (AliMCParticle*)(mcEvent->GetTrack(GeneratedAntiProtonArray->at(track1)));

	MomentumAntiProton[0] = AntiProton->Px();
	MomentumAntiProton[1] = AntiProton->Py();
	MomentumAntiProton[2] = AntiProton->Pz();

	TLorentzVector LorentzVectorAntiProton;
	LorentzVectorAntiProton.SetXYZM(MomentumAntiProton[0],MomentumAntiProton[1],MomentumAntiProton[2],MassProton);
	PtAntiProton = LorentzVectorAntiProton.Pt();
	EtaAntiProton = LorentzVectorAntiProton.Eta();

	for(int track2 = 0; track2 < GeneratedAntiDeuteronArray->size(); track2++)
	  {

	    AliMCParticle *AntiDeuteron = (AliMCParticle*)(mcEvent->GetTrack(GeneratedAntiDeuteronArray->at(track2)));

	    MomentumAntiDeuteron[0] = AntiDeuteron->Px();
	    MomentumAntiDeuteron[1] = AntiDeuteron->Py();
	    MomentumAntiDeuteron[2] = AntiDeuteron->Pz();

	    TLorentzVector LorentzVectorAntiDeuteron;
	    LorentzVectorAntiDeuteron.SetXYZM(MomentumAntiDeuteron[0],MomentumAntiDeuteron[1],MomentumAntiDeuteron[2],MassDeuteron);
	    PtAntiDeuteron = LorentzVectorAntiDeuteron.Pt();
	    EtaAntiDeuteron = LorentzVectorAntiDeuteron.Eta();
      
	    TLorentzVector LorentzVectorPair;
	    LorentzVectorPair = LorentzVectorAntiProton + LorentzVectorAntiDeuteron;

	    PtPair = LorentzVectorPair.Pt();
	    EtaPair = LorentzVectorPair.Eta();

	    // pair cut
	    if((PtPair < PairPtLimit1)  || (PtPair > PairPtLimit2)) continue;

	    RelativeMomentum = CalculateRelativeMomentum(LorentzVectorPair,LorentzVectorAntiProton,LorentzVectorAntiDeuteron);

		fHistPtAntiProtonGenPairPtCut->Fill(PtAntiProton);
		fHistEtaAntiProtonGenPairPtCut->Fill(EtaAntiProton);
		fHistPtAntiDeuteronGenPairPtCut->Fill(PtAntiDeuteron);
		fHistEtaAntiDeuteronGenPairPtCut->Fill(EtaAntiDeuteron);
		fHistPtAntiHelium3GenPairPtCut->Fill(PtPair);
		fHistEtaAntiHelium3GenPairPtCut->Fill(EtaPair);
		fHistSEDAntiPairGen->Fill(RelativeMomentum);
		fHistPtAntiParticlesGen->Fill(PtAntiDeuteron,PtAntiProton);


	  } // end of loop for generated antiparticle pairs (antideuteron loop)

      } // end of loop for generated antiparticle pairs (antiproton loop)



// +++++++++++++++++++++++++++++++++++++++++++++++++++++
// +++++++++++ Loops for reconstructed particles +++++++
// +++++++++++++++++++++++++++++++++++++++++++++++++++++



    // loop for reconstructed protons
    for(int track = 0; track < fESD->GetNumberOfTracks(); track++)
      {

	TPC_OK = false;
	TOF_OK = false;

	AliESDtrack *ProtonTrack = dynamic_cast<AliESDtrack*>(fESD->GetTrack(track));
	LabelProton = TMath::Abs(ProtonTrack->GetLabel());
	AliMCParticle *Proton = (AliMCParticle*)(mcEvent->GetTrack(LabelProton));

	MomentumProton[0] = ProtonTrack->Px();
	MomentumProton[1] = ProtonTrack->Py();
	MomentumProton[2] = ProtonTrack->Pz();

	TLorentzVector LorentzVectorProton;
	LorentzVectorProton.SetXYZM(MomentumProton[0],MomentumProton[1],MomentumProton[2],MassProton);

	PtProton = LorentzVectorProton.Pt();
	EtaProton = LorentzVectorProton.Eta();

	// check if particle is a proton
	if(!(Proton->PdgCode() == ProtonPDG))	continue;

	fHistPtProtonRecPDG->Fill(PtProton);
	fHistEtaProtonRecPDG->Fill(EtaProton);

	// check if proton is a primary particle
	if(!Proton->IsPhysicalPrimary()) continue;

	fHistPtProtonRecPrimary->Fill(PtProton);
	fHistEtaProtonRecPrimary->Fill(EtaProton);

	// apply track cuts
	if(!fESDtrackCutsProton->AcceptTrack(ProtonTrack)) continue;
	if(ProtonTrack->GetSign() < 1) continue;

	// nsigma cuts
	AliPIDResponse::EDetPidStatus statusTPC = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTPC,ProtonTrack);
	AliPIDResponse::EDetPidStatus statusTOF = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTOF,ProtonTrack);

	if(statusTPC == AliPIDResponse::kDetPidOk) TPC_OK = true;
	if(statusTOF == AliPIDResponse::kDetPidOk) TOF_OK = true;
	if(!(TPC_OK)) continue;
	if(!(TOF_OK)) continue;

	nSigmaTPCproton = fPIDResponse->NumberOfSigmasTPC(ProtonTrack,AliPID::kProton);
	nSigmaTOFproton = fPIDResponse->NumberOfSigmasTOF(ProtonTrack,AliPID::kProton);
	nSigmaTPCTOFproton = TMath::Sqrt(pow(nSigmaTPCproton,2.) + pow(nSigmaTOFproton,2.));

	Double_t ProtonpTPC = ProtonTrack->GetInnerParam()->GetP();

	if(ProtonpTPC <= ProtonpTPCThreshold){

	  if(!(nSigmaTPCproton < 3.0)) continue;

	}else{
	    
	  if(!(nSigmaTPCTOFproton < 3.0)) continue;

        }

	// DCA cuts
	ProtonTrack->GetImpactParameters(xv_Proton,yv_Proton);
	ProtonDCAxy = xv_Proton[0];
	ProtonDCAz  = xv_Proton[1];
	if((ProtonDCAxy > ProtonDCAxyMax) || (ProtonDCAz > ProtonDCAzMax)) continue;


	fHistPtProtonRecTrackCuts->Fill(PtProton);
	fHistEtaProtonRecTrackCuts->Fill(EtaProton);

	ReconstructedProtonArray->push_back(LabelProton);


      }	// end of loop for reconstructed protons




    // loop for reconstructed deuterons
    for(int track = 0; track < fESD->GetNumberOfTracks(); track++)
      {

	TPC_OK = false;
	TOF_OK = false;

	AliESDtrack *DeuteronTrack = dynamic_cast<AliESDtrack*>(fESD->GetTrack(track));
	LabelDeuteron = TMath::Abs(DeuteronTrack->GetLabel());
	AliMCParticle *Deuteron = (AliMCParticle*)(mcEvent->GetTrack(LabelDeuteron));

	MomentumDeuteron[0] = DeuteronTrack->Px();
	MomentumDeuteron[1] = DeuteronTrack->Py();
	MomentumDeuteron[2] = DeuteronTrack->Pz();

	TLorentzVector LorentzVectorDeuteron;
	LorentzVectorDeuteron.SetXYZM(MomentumDeuteron[0],MomentumDeuteron[1],MomentumDeuteron[2],MassDeuteron);

	PtDeuteron = LorentzVectorDeuteron.Pt();
	EtaDeuteron = LorentzVectorDeuteron.Eta();


	// check if particle is a deuteron
	if(!(Deuteron->PdgCode() == DeuteronPDG)) continue;

	    fHistPtDeuteronRecPDG->Fill(PtDeuteron);
	    fHistEtaDeuteronRecPDG->Fill(EtaDeuteron);


	// check if deuteron is a primary particle
	if((!Deuteron->IsPhysicalPrimary())) continue;

	    fHistPtDeuteronRecPrimary->Fill(PtDeuteron);
	    fHistEtaDeuteronRecPrimary->Fill(EtaDeuteron);


	// apply track cuts
	if(!fESDtrackCutsDeuteron->AcceptTrack(DeuteronTrack)) continue;
	if(DeuteronTrack->GetSign() < 1) continue;

	// fill dE/dx with particles and antiparticles according to the production period
	if(LHC18a2a == true){

	  fHistdEdx_LHC18a2a->Fill(DeuteronTrack->GetInnerParam()->GetP() * DeuteronTrack->GetSign(),DeuteronTrack->GetTPCsignal());

	  }else if(LHC18a2b == true){

	    fHistdEdx_LHC18a2b->Fill(DeuteronTrack->GetInnerParam()->GetP() * DeuteronTrack->GetSign(),DeuteronTrack->GetTPCsignal());

	    }else if(LHC18a2b4 == true){

	      fHistdEdx_LHC18a2b4->Fill(DeuteronTrack->GetInnerParam()->GetP() * DeuteronTrack->GetSign(),DeuteronTrack->GetTPCsignal());

	      }else if(LHC20l7a == true){

		fHistdEdx_LHC20l7a->Fill(DeuteronTrack->GetInnerParam()->GetP() * DeuteronTrack->GetSign(),DeuteronTrack->GetTPCsignal());

		}

	AliPIDResponse::EDetPidStatus statusTPC = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTPC,DeuteronTrack);
	AliPIDResponse::EDetPidStatus statusTOF = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTOF,DeuteronTrack);

	if(statusTPC == AliPIDResponse::kDetPidOk) TPC_OK = true;
	if(statusTOF == AliPIDResponse::kDetPidOk) TOF_OK = true;

	if(!(TPC_OK)) continue;
	if(!(TOF_OK)) continue;

	nSigmaTPCdeuteron = Bethe(*DeuteronTrack,AliPID::ParticleMass(AliPID::kDeuteron),1,fBetheParametersDeuteron);
	nSigmaTOFdeuteron = fPIDResponse->NumberOfSigmasTOF(DeuteronTrack,AliPID::kDeuteron);
	nSigmaTPCTOFdeuteron = TMath::Sqrt(pow(nSigmaTPCdeuteron,2.) + pow(nSigmaTOFdeuteron,2.));

	Double_t DeuteronpTPC = DeuteronTrack->GetInnerParam()->GetP();

	if(DeuteronpTPC <= DeuteronpTPCThreshold){

	  if(!(nSigmaTPCdeuteron < 3.0)) continue;

	}else{
	    
	  if(!(nSigmaTPCTOFdeuteron < 3.0)) continue;

	}

	// DCA cuts
	DeuteronTrack->GetImpactParameters(xv_Deuteron,yv_Deuteron);
	DeuteronDCAxy = xv_Deuteron[0];
	DeuteronDCAz  = xv_Deuteron[1];
	if((DeuteronDCAxy > DeuteronDCAxyMax) || (DeuteronDCAz > DeuteronDCAzMax)) continue;

	fHistPtDeuteronRecTrackCuts->Fill(PtDeuteron);
	fHistEtaDeuteronRecTrackCuts->Fill(EtaDeuteron);

	ReconstructedDeuteronArray->push_back(LabelDeuteron);

      }	// end of loop for reconstructed deuterons




    // loop for reconstructed particle pairs
    for(int track1 = 0; track1 < ReconstructedProtonArray->size(); track1++)
      {

	AliESDtrack *ProtonTrack = (AliESDtrack*)(mcEvent->GetTrack(ReconstructedProtonArray->at(track1)));

	MomentumProton[0] = ProtonTrack->Px();
	MomentumProton[1] = ProtonTrack->Py();
	MomentumProton[2] = ProtonTrack->Pz();

	TLorentzVector LorentzVectorProton;
	LorentzVectorProton.SetXYZM(MomentumProton[0],MomentumProton[1],MomentumProton[2],MassProton);
	PtProton = LorentzVectorProton.Pt();
	EtaProton = LorentzVectorProton.Eta();

	for(int track2 = 0; track2 < ReconstructedDeuteronArray->size(); track2++)
	  {

	    AliESDtrack *DeuteronTrack = (AliESDtrack*)(mcEvent->GetTrack(ReconstructedDeuteronArray->at(track2)));

	    MomentumDeuteron[0] = DeuteronTrack->Px();
	    MomentumDeuteron[1] = DeuteronTrack->Py();
	    MomentumDeuteron[2] = DeuteronTrack->Pz();

	    TLorentzVector LorentzVectorDeuteron;
	    LorentzVectorDeuteron.SetXYZM(MomentumDeuteron[0],MomentumDeuteron[1],MomentumDeuteron[2],MassDeuteron);
	    PtDeuteron = LorentzVectorDeuteron.Pt();
	    EtaDeuteron = LorentzVectorDeuteron.Eta();


	    TLorentzVector LorentzVectorPair;
	    LorentzVectorPair = LorentzVectorProton + LorentzVectorDeuteron;

	    EtaPair = LorentzVectorPair.Eta();
	    PtPair = LorentzVectorPair.Pt();

	    // apply pair cut
	    if((PtPair < PairPtLimit1)	  || (PtPair > PairPtLimit2)) continue;

	    RelativeMomentum = CalculateRelativeMomentum(LorentzVectorPair,LorentzVectorProton,LorentzVectorDeuteron);

	    fHistPtProtonRecPairPtCut->Fill(PtProton);
	    fHistEtaProtonRecPairPtCut->Fill(EtaProton);
	    fHistPtDeuteronRecPairPtCut->Fill(PtDeuteron);
	    fHistEtaDeuteronRecPairPtCut->Fill(EtaDeuteron);
	    fHistPtHelium3RecPairPtCut->Fill(PtPair);
	    fHistEtaHelium3RecPairPtCut->Fill(EtaPair);
	    fHistSEDPairRec->Fill(RelativeMomentum);
	    fHistPtParticlesRec->Fill(PtDeuteron,PtProton);

	  } // end of loop for reconstructed particle pairs (deuteron loop)

      } // end of loop for reconstructed particle pairs (proton loop)



// +++++++++++++++++++++++++++++++++++++++++++++++++++++
// +++++++++++ Loops for reconstructed antiparticles +++
// +++++++++++++++++++++++++++++++++++++++++++++++++++++


    // loop for reconstructed antiprotons
    for(int track = 0; track < fESD->GetNumberOfTracks(); track++)
      {

	TPC_OK = false;
	TOF_OK = false;

	AliESDtrack *AntiProtonTrack = dynamic_cast<AliESDtrack*>(fESD->GetTrack(track));

	LabelAntiProton = TMath::Abs(AntiProtonTrack->GetLabel());
	AliMCParticle *AntiProton = (AliMCParticle*)(mcEvent->GetTrack(LabelAntiProton));

	MomentumAntiProton[0] = AntiProtonTrack->Px();
	MomentumAntiProton[1] = AntiProtonTrack->Py();
	MomentumAntiProton[2] = AntiProtonTrack->Pz();

	TLorentzVector LorentzVectorAntiProton;
	LorentzVectorAntiProton.SetXYZM(MomentumAntiProton[0],MomentumAntiProton[1],MomentumAntiProton[2],MassProton);

	PtAntiProton = LorentzVectorAntiProton.Pt();
	EtaAntiProton = LorentzVectorAntiProton.Eta();

	// check if particle is an antiproton
	if(!(AntiProton->PdgCode() == -ProtonPDG))  continue;

	fHistPtAntiProtonRecPDG->Fill(PtAntiProton);
	fHistEtaAntiProtonRecPDG->Fill(EtaAntiProton);

	// check if antiproton is a primary particle
	if(!AntiProton->IsPhysicalPrimary()) continue;

	fHistPtAntiProtonRecPrimary->Fill(PtAntiProton);
	fHistEtaAntiProtonRecPrimary->Fill(EtaAntiProton);

	// apply track cuts
	if(!fESDtrackCutsProton->AcceptTrack(AntiProtonTrack)) continue;
	if(AntiProtonTrack->GetSign() > -1) continue;

	AliPIDResponse::EDetPidStatus statusTPC = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTPC,AntiProtonTrack);
	AliPIDResponse::EDetPidStatus statusTOF = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTOF,AntiProtonTrack);

	if(statusTPC == AliPIDResponse::kDetPidOk) TPC_OK = true;
	if(statusTOF == AliPIDResponse::kDetPidOk) TOF_OK = true;
	if(!(TPC_OK)) continue;
	if(!(TOF_OK)) continue;

	nSigmaTPCproton = fPIDResponse->NumberOfSigmasTPC(AntiProtonTrack,AliPID::kProton);
	nSigmaTOFproton = fPIDResponse->NumberOfSigmasTOF(AntiProtonTrack,AliPID::kProton);
	nSigmaTPCTOFproton = TMath::Sqrt(pow(nSigmaTPCproton,2.) + pow(nSigmaTOFproton,2.));

	Double_t ProtonpTPC = AntiProtonTrack->GetInnerParam()->GetP();

	if(ProtonpTPC <= ProtonpTPCThreshold){

	  if(!(nSigmaTPCproton < 3.0)) continue;

	}else{
	    
	  if(!(nSigmaTPCTOFproton < 3.0)) continue;

        }

	// DCA cuts
	AntiProtonTrack->GetImpactParameters(xv_Proton,yv_Proton);
	ProtonDCAxy = xv_Proton[0];
	ProtonDCAz  = xv_Proton[1];
	if((ProtonDCAxy > ProtonDCAxyMax) || (ProtonDCAz > ProtonDCAzMax)) continue;


	fHistPtAntiProtonRecTrackCuts->Fill(PtAntiProton);
	fHistEtaAntiProtonRecTrackCuts->Fill(EtaAntiProton);

	ReconstructedAntiProtonArray->push_back(LabelAntiProton);

      } // end of loop for reconstructed antiprotons


      // loop for reconstructed antideuterons
      for(int track = 0; track < fESD->GetNumberOfTracks(); track++)
	{

	  TPC_OK = false;
	  TOF_OK = false;

	  AliESDtrack *AntiDeuteronTrack = dynamic_cast<AliESDtrack*>(fESD->GetTrack(track));
	  LabelAntiDeuteron = TMath::Abs(AntiDeuteronTrack->GetLabel());

	  AliMCParticle *AntiDeuteron = (AliMCParticle*)(mcEvent->GetTrack(LabelAntiDeuteron));

	  MomentumAntiDeuteron[0] = AntiDeuteronTrack->Px();
	  MomentumAntiDeuteron[1] = AntiDeuteronTrack->Py();
	  MomentumAntiDeuteron[2] = AntiDeuteronTrack->Pz();

	  TLorentzVector LorentzVectorAntiDeuteron;
	  LorentzVectorAntiDeuteron.SetXYZM(MomentumAntiDeuteron[0],MomentumAntiDeuteron[1],MomentumAntiDeuteron[2],MassDeuteron);

	  PtAntiDeuteron = LorentzVectorAntiDeuteron.Pt();
	  EtaAntiDeuteron = LorentzVectorAntiDeuteron.Eta();

	  // check if particle is an antideuteron
	  if(!(AntiDeuteron->PdgCode() == -DeuteronPDG)) continue;

	  fHistPtAntiDeuteronRecPDG->Fill(PtAntiDeuteron);
	  fHistEtaAntiDeuteronRecPDG->Fill(EtaAntiDeuteron);

	  // check if antideuteron is a primary particle
	  if((!AntiDeuteron->IsPhysicalPrimary())) continue;

	  fHistPtAntiDeuteronRecPrimary->Fill(PtAntiDeuteron);
	  fHistEtaAntiDeuteronRecPrimary->Fill(EtaAntiDeuteron);

	  // apply track cuts
	  if(!fESDtrackCutsDeuteron->AcceptTrack(AntiDeuteronTrack)) continue;
	  if(AntiDeuteronTrack->GetSign() > -1) continue;

	  AliPIDResponse::EDetPidStatus statusTPC = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTPC,AntiDeuteronTrack);
	  AliPIDResponse::EDetPidStatus statusTOF = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTOF,AntiDeuteronTrack);

	  if(statusTPC == AliPIDResponse::kDetPidOk) TPC_OK = true;
	  if(statusTOF == AliPIDResponse::kDetPidOk) TOF_OK = true;

	  if(!(TPC_OK)) continue;
	  if(!(TOF_OK)) continue;

	  //nSigmaTPCdeuteron = fPIDResponse->NumberOfSigmasTPC(Track2,AliPID::kDeuteron);
	  nSigmaTPCdeuteron = Bethe(*AntiDeuteronTrack,AliPID::ParticleMass(AliPID::kDeuteron),1,fBetheParametersDeuteron);
	  nSigmaTOFdeuteron = fPIDResponse->NumberOfSigmasTOF(AntiDeuteronTrack,AliPID::kDeuteron);
	  nSigmaTPCTOFdeuteron = TMath::Sqrt(pow(nSigmaTPCdeuteron,2.) + pow(nSigmaTOFdeuteron,2.));

	  Double_t DeuteronpTPC = AntiDeuteronTrack->GetInnerParam()->GetP();

	    if(DeuteronpTPC <= DeuteronpTPCThreshold){

	      if(!(nSigmaTPCdeuteron < 3.0)) continue;

	    }else{
	    
	      if(!(nSigmaTPCTOFdeuteron < 3.0)) continue;

	    }

	    // DCA cuts
	  AntiDeuteronTrack->GetImpactParameters(xv_Deuteron,yv_Deuteron);
	  DeuteronDCAxy = xv_Deuteron[0];
	  DeuteronDCAz  = xv_Deuteron[1];
	  if((DeuteronDCAxy > DeuteronDCAxyMax) || (DeuteronDCAz > DeuteronDCAzMax)) continue;

	  fHistPtAntiDeuteronRecTrackCuts->Fill(PtAntiDeuteron);
	  fHistEtaAntiDeuteronRecTrackCuts->Fill(EtaAntiDeuteron);

	  ReconstructedAntiDeuteronArray->push_back(LabelAntiDeuteron);

	} // end of loop for reconstructed antideuterons



    // loop for reconstructed antiparticle pairs
    for(int track1 = 0; track1 < ReconstructedAntiProtonArray->size(); track1++)
      {

	AliESDtrack *AntiProtonTrack = (AliESDtrack*)(mcEvent->GetTrack(ReconstructedAntiProtonArray->at(track1)));

	MomentumAntiProton[0] = AntiProtonTrack->Px();
	MomentumAntiProton[1] = AntiProtonTrack->Py();
	MomentumAntiProton[2] = AntiProtonTrack->Pz();

	TLorentzVector LorentzVectorAntiProton;
	LorentzVectorAntiProton.SetXYZM(MomentumAntiProton[0],MomentumAntiProton[1],MomentumAntiProton[2],MassProton);
	PtAntiProton = LorentzVectorAntiProton.Pt();
	EtaAntiProton = LorentzVectorAntiProton.Eta();

	for(int track2 = 0; track2 < ReconstructedAntiDeuteronArray->size(); track2++)
	  {

	    AliESDtrack *AntiDeuteronTrack = (AliESDtrack*)(mcEvent->GetTrack(ReconstructedAntiDeuteronArray->at(track2)));

	    MomentumAntiDeuteron[0] = AntiDeuteronTrack->Px();
	    MomentumAntiDeuteron[1] = AntiDeuteronTrack->Py();
	    MomentumAntiDeuteron[2] = AntiDeuteronTrack->Pz();

	    TLorentzVector LorentzVectorAntiDeuteron;
	    LorentzVectorAntiDeuteron.SetXYZM(MomentumAntiDeuteron[0],MomentumAntiDeuteron[1],MomentumAntiDeuteron[2],MassDeuteron);
	    PtAntiDeuteron = LorentzVectorAntiDeuteron.Pt();
	    EtaAntiDeuteron = LorentzVectorAntiDeuteron.Eta();


	    TLorentzVector LorentzVectorPair;
	    LorentzVectorPair = LorentzVectorAntiProton + LorentzVectorAntiDeuteron;

	    EtaPair = LorentzVectorPair.Eta();
	    PtPair = LorentzVectorPair.Pt();

	    // apply pair cut
	    if((PtPair < PairPtLimit1)	  || (PtPair > PairPtLimit2)) continue;

	    RelativeMomentum = CalculateRelativeMomentum(LorentzVectorPair,LorentzVectorAntiProton,LorentzVectorAntiDeuteron);

	    fHistPtAntiProtonRecPairPtCut->Fill(PtAntiProton);
	    fHistEtaAntiProtonRecPairPtCut->Fill(EtaAntiProton);
	    fHistPtAntiDeuteronRecPairPtCut->Fill(PtAntiDeuteron);
	    fHistEtaAntiDeuteronRecPairPtCut->Fill(EtaAntiDeuteron);
	    fHistPtAntiHelium3RecPairPtCut->Fill(PtPair);
	    fHistEtaAntiHelium3RecPairPtCut->Fill(EtaPair);
	    fHistSEDAntiPairRec->Fill(RelativeMomentum);
	    fHistPtAntiParticlesRec->Fill(PtAntiDeuteron,PtAntiProton);

	  } // end of loop for reconstructed antiparticle pairs (antideuteron loop)

      } // end of loop for reconstructed antiparticle pairs (antiproton loop)





  fHistEventCounter->Fill(0.5);

  } // end of "if is MC"

  PostData(1,fHistList);


} // end of UserExec




void AliAnalysisTaskDeuteronProtonEfficiency::Terminate(Option_t *)
{



} // end of Terminate



// Calculates number of sigma deviation from expected dE/dx in TPC
double AliAnalysisTaskDeuteronProtonEfficiency::Bethe(const AliESDtrack &track, double mass, int charge, double *params){

  double expected = charge*charge*AliExternalTrackParam::BetheBlochAleph(charge*track.GetInnerParam()->GetP()/mass,params[0],params[1],params[2],params[3],params[4]);
  double sigma = expected * params[5];

  if(TMath::IsNaN(expected)) return -999;
  return (track.GetTPCsignal() - expected) / sigma;


} // end of Bethe



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


