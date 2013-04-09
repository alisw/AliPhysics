#ifndef ALIANALYSISTASKSE_H
#include <Riostream.h>
#include <TROOT.h>
#include <TFile.h>
#include <TCint.h>
#include <TChain.h>
#include <TTree.h>
#include <TKey.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TList.h>
#include <TClonesArray.h>
#include <TObject.h>
#include <TMath.h>
#include <TSystem.h>
#include <TInterpreter.h>
#include <TH1.h>
#include "AliAnalysisTask.h"
#include "AliCentrality.h"
#include "AliStack.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliAODEvent.h"
#include "AliAODHandler.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskSE.h"
#endif

#include <TRandom3.h>
#include "AliGenPythiaEventHeader.h"
#include "AliAODMCHeader.h"
#include "AliMCEvent.h"
#include "AliLog.h"
#include <AliEmcalJet.h>
#include <AliPicoTrack.h>
#include "AliVEventHandler.h"
#include "AliVParticle.h"
#include "AliAnalysisUtils.h"


#include "AliAnalysisTaskChargedJetsPA.h"

//TODO: Not accessing the particles when using MC
//TODO: FillHistogram can be done better with virtual TH1(?)
ClassImp(AliAnalysisTaskChargedJetsPA)

// ######################################################################################## DEFINE HISTOGRAMS
void AliAnalysisTaskChargedJetsPA::Init()
{
  #ifdef DEBUGMODE
    AliInfo("Creating histograms.");
  #endif

  AddHistogram1D<TH1D>("hNumberEvents", "Number of events (0 = before, 1 = after vertex cuts)", "", 2, 0, 2, "#Delta z(cm)","N^{Events}/cut");

  // NOTE: Jet histograms
  if (fAnalyzeJets)
  {
    // ######## Jet spectra
    AddHistogram2D<TH2D>("hJetPt", "Jets p_{T} distribution", "", 500, -50., 200., fNumberOfCentralityBins, 0, 100, "p_{T} (GeV/c)","Centrality","dN^{Jets}/dp_{T}");
    AddHistogram2D<TH2D>("hJetPtBgrdSubtractedKTImprovedCMS", "Jets p_{T} distribution, KT background (Improved CMS) subtracted", "", 500, -50., 200., fNumberOfCentralityBins, 0, 100, "p_{T} (GeV/c)","Centrality","dN^{Jets}/dp_{T}");    

    if(fAnalyzeDeprecatedBackgrounds)
    {
      AddHistogram2D<TH2D>("hJetPtBgrdSubtractedTR", "Jets p_{T} distribution, TR background (Cone R=0.6 around jets excluded) subtracted", "", 500, -50., 200., fNumberOfCentralityBins, 0, 100, "p_{T} (GeV/c)","Centrality","dN^{Jets}/dp_{T}");
      AddHistogram2D<TH2D>("hJetPtBgrdSubtractedRC", "Jets p_{T} distribution, RC background subtracted", "", 500, -50., 200.,fNumberOfCentralityBins, 0, 100, "p_{T} (GeV/c)","Centrality","dN^{Jets}/dp_{T}");
      AddHistogram2D<TH2D>("hJetPtBgrdSubtractedKTPbPb", "Jets p_{T} distribution, KT background (PbPb w/o ghosts) subtracted", "", 500, -50., 200., fNumberOfCentralityBins, 0, 100, "p_{T} (GeV/c)","Centrality","dN^{Jets}/dp_{T}");
      AddHistogram2D<TH2D>("hJetPtBgrdSubtractedKTPbPbWithGhosts", "Jets p_{T} distribution, KT background (PbPb w/ ghosts) subtracted", "", 500, -50., 200., fNumberOfCentralityBins, 0, 100, "p_{T} (GeV/c)","Centrality","dN^{Jets}/dp_{T}");
      AddHistogram2D<TH2D>("hJetPtBgrdSubtractedKTCMS", "Jets p_{T} distribution, KT background (CMS) subtracted", "", 500, -50., 200., fNumberOfCentralityBins, 0, 100, "p_{T} (GeV/c)","Centrality","dN^{Jets}/dp_{T}");    
      AddHistogram2D<TH2D>("hJetPtBgrdSubtractedKTMean", "Jets p_{T} distribution, KT background (Mean) subtracted", "", 500, -50., 200., fNumberOfCentralityBins, 0, 100, "p_{T} (GeV/c)","Centrality","dN^{Jets}/dp_{T}");    
      AddHistogram2D<TH2D>("hJetPtBgrdSubtractedKTTrackLike", "Jets p_{T} distribution, KT background (track-like) subtracted", "", 500, -50., 200., fNumberOfCentralityBins, 0, 100, "p_{T} (GeV/c)","Centrality","dN^{Jets}/dp_{T}");
    }

    // ######## Jet stuff
    AddHistogram1D<TH1D>("hJetCountAll", "Number of Jets", "", 200, 0., 200., "N jets","dN^{Events}/dN^{Jets}");
    AddHistogram1D<TH1D>("hJetCountAccepted", "Number of accepted Jets", "", 200, 0., 200., "N jets","dN^{Events}/dN^{Jets}");
    AddHistogram1D<TH1D>("hLeadingJetPt", "Leading jet p_{T}", "", 500,  0, 100, "p_{T} (GeV/c)","dN^{Jets}/dp_{T}");
    AddHistogram1D<TH1D>("hSecondLeadingJetPt", "Second Leading jet p_{T}", "", 500,  0, 100, "p_{T} (GeV/c)","dN^{Jets}/dp_{T}");
    AddHistogram1D<TH1D>("hJetDeltaPhi", "Jets combinatorial #Delta #phi", "", 250, 0., TMath::Pi(), "#Delta #phi","dN^{Jets}/d(#Delta #phi)");
    AddHistogram1D<TH1D>("hLeadingJetDeltaPhi", "1st and 2nd leading jet #Delta #phi", "", 250, 0., TMath::Pi(), "#Delta #phi","dN^{Jets}/d(#Delta #phi)");

    // ########## Dijet stuff
    AddHistogram1D<TH1D>("hDijetLeadingJetPt", "Dijet leading jet p_{T} distribution", "", 500, 0., 100., "p_{T} (GeV/c)","dN^{Jets}/dp_{T}");
    AddHistogram1D<TH1D>("hDijetConstituentsPt", "Dijet constituents p_{T} distribution", "", 500, 0., 100., "p_{T} (GeV/c)","dN^{Jets}/dp_{T}");
    AddHistogram2D<TH2D>("hDijetPtCorrelation", "Dijet constituents p_{T} correlation", "COLZ", 500, 5., 100., 500, 5., 100., "1st leading jet p_{T} (GeV/c)","2nd leading jet p_{T} (GeV/c)","dN^{Dijets}/d^{2}p_{T}");
  }

  // NOTE: Jet background histograms
  if (fAnalyzeBackground)
  {
    // ########## Default background estimates
    AddHistogram2D<TH2D>("hKTBackgroundImprovedCMS", "KT background density (Improved CMS approach)", "LEGO2", 400, 0., 40., fNumberOfCentralityBins, 0, 100, "#rho (GeV/c)","Centrality", "dN^{Events}/d#rho");
    AddHistogram2D<TH2D>("hDeltaPtKTImprovedCMS", "Background fluctuations #delta p_{T} (KT, Improved CMS-like)", "", 1201, -40.0, 40.0, fNumberOfCentralityBins, 0, 100, "#delta p_{T} (GeV/c)","Centrality","dN^{Jets}/d#delta p_{T}");
    AddHistogram2D<TH2D>("hDeltaPtNoBackground", "Background fluctuations #delta p_{T} (No background)", "", 1201, -40.0, 40.0, fNumberOfCentralityBins, 0, 100, "#delta p_{T} (GeV/c)","Centrality","dN^{Jets}/d#delta p_{T}");
    AddHistogram2D<TH2D>("hDeltaPtNoBackgroundNoEmptyCones", "Background fluctuations #delta p_{T} (No background, no empty cones)", "", 1201, -40.0, 40.0, fNumberOfCentralityBins, 0, 100, "#delta p_{T} (GeV/c)","Centrality","dN^{Jets}/d#delta p_{T}");

    AddHistogram2D<TProfile2D>("hJetPtSubtractedRhoKTImprovedCMS", "Mean subtracted KT (CMS w/o signal) background from jets", "COLZ", 600, 0, 150, fNumberOfCentralityBins, 0, 100, "Jet p_{T}", "Centrality", "#rho mean");
    AddHistogram1D<TProfile>("hKTMeanBackgroundImprovedCMS", "KT background mean (Improved CMS approach)", "", 100, 0, 100, "Centrality", "#rho mean");

    AddHistogram2D<TH2D>("hDijetBackground", "Background density (dijets excluded)", "", 200, 0., 20., fNumberOfCentralityBins, 0, 100, "#rho (GeV/c)","Centrality", "dN^{Events}/d#rho");
    AddHistogram2D<TH2D>("hDijetBackgroundPerpendicular", "Background density (dijets excluded)", "", 200, 0., 20., fNumberOfCentralityBins, 0, 100, "#rho (GeV/c)","Centrality", "dN^{Events}/d#rho");

    if(fAnalyzeDeprecatedBackgrounds)
    {
      // ########## Different background estimates
      AddHistogram2D<TH2D>("hRCBackground", "RC background density (Signal jets excluded)", "LEGO2", 400, 0., 40., fNumberOfCentralityBins, 0, 100, "#rho (GeV/c)","Centrality", "dN^{Events}/d#rho");
      AddHistogram2D<TH2D>("hKTBackgroundPbPb", "KT background density (PbPb approach, no ghosts)", "LEGO2", 400, 0., 40., fNumberOfCentralityBins, 0, 100, "#rho (GeV/c)","Centrality", "dN^{Events}/d#rho");
      AddHistogram2D<TH2D>("hKTBackgroundPbPbWithGhosts", "KT background density (PbPb approach w/ ghosts)", "LEGO2", 400, 0., 40., fNumberOfCentralityBins, 0, 100, "#rho (GeV/c)","Centrality", "dN^{Events}/d#rho");
      AddHistogram2D<TH2D>("hKTBackgroundCMS", "KT background density (CMS approach)", "LEGO2", 400, 0., 40., fNumberOfCentralityBins, 0, 100, "#rho (GeV/c)","Centrality", "dN^{Events}/d#rho");
      AddHistogram2D<TH2D>("hKTBackgroundMean", "KT background density (Mean approach)", "LEGO2", 400, 0., 40., fNumberOfCentralityBins, 0, 100, "#rho (GeV/c)","Centrality", "dN^{Events}/d#rho");
      AddHistogram2D<TH2D>("hKTBackgroundTrackLike", "KT background density (Track-like approach)", "LEGO2", 400, 0., 40., fNumberOfCentralityBins, 0, 100, "#rho (GeV/c)","Centrality", "dN^{Events}/d#rho");

      AddHistogram2D<TH2D>("hTRBackgroundNoExcl", "TR background density (No signal excluded)", "LEGO2", 400, 0., 40., fNumberOfCentralityBins, 0, 100, "#rho (GeV/c)","Centrality", "dN^{Events}/d#rho");
      AddHistogram2D<TH2D>("hTRBackgroundCone02", "TR background density (Cones R=0.2 around signal jets excluded)", "LEGO2", 400, 0., 40., fNumberOfCentralityBins, 0, 100, "#rho (GeV/c)","Centrality", "dN^{Events}/d#rho");
      AddHistogram2D<TH2D>("hTRBackgroundCone04", "TR background density (Cones R=0.4 around signal jets excluded)", "LEGO2", 400, 0., 40., fNumberOfCentralityBins, 0, 100, "#rho (GeV/c)","Centrality", "dN^{Events}/d#rho");
      AddHistogram2D<TH2D>("hTRBackgroundCone06", "TR background density (Cones R=0.6 around signal jets excluded)", "LEGO2", 400, 0., 40., fNumberOfCentralityBins, 0, 100, "#rho (GeV/c)","Centrality", "dN^{Events}/d#rho");
      AddHistogram2D<TH2D>("hTRBackgroundCone08", "TR background density (Cones R=0.8 around signal jets excluded)", "LEGO2", 400, 0., 40., fNumberOfCentralityBins, 0, 100, "#rho (GeV/c)","Centrality", "dN^{Events}/d#rho");
      AddHistogram2D<TH2D>("hTRBackgroundExact",  "TR background density (signal jets exactly excluded)", "LEGO2", 400, 0., 40., fNumberOfCentralityBins, 0, 100, "#rho (GeV/c)","Centrality", "dN^{Events}/d#rho");

      // ########## Delta Pt
      AddHistogram2D<TH2D>("hDeltaPtKTPbPb", "Background fluctuations #delta p_{T} (KT, PbPb w/o ghosts)", "", 1201, -40.0, 40.0, fNumberOfCentralityBins, 0, 100, "#delta p_{T} (GeV/c)","Centrality","dN^{Jets}/d#delta p_{T}");
      AddHistogram2D<TH2D>("hDeltaPtKTPbPbWithGhosts", "Background fluctuations #delta p_{T} (KT, PbPb w/ ghosts)", "", 1201, -40.0, 40.0, fNumberOfCentralityBins, 0, 100, "#delta p_{T} (GeV/c)","Centrality","dN^{Jets}/d#delta p_{T}");
      AddHistogram2D<TH2D>("hDeltaPtKTCMS", "Background fluctuations #delta p_{T} (KT, CMS-like)", "", 1201, -40.0, 40.0, fNumberOfCentralityBins, 0, 100, "#delta p_{T} (GeV/c)","Centrality","dN^{Jets}/d#delta p_{T}");
      AddHistogram2D<TH2D>("hDeltaPtKTMean", "Background fluctuations #delta p_{T} (KT, Mean)", "", 1201, -40.0, 40.0, fNumberOfCentralityBins, 0, 100, "#delta p_{T} (GeV/c)","Centrality","dN^{Jets}/d#delta p_{T}");
      AddHistogram2D<TH2D>("hDeltaPtKTTrackLike", "Background fluctuations #delta p_{T} (KT, track-like)", "", 1201, -40.0, 40.0, fNumberOfCentralityBins, 0, 100, "#delta p_{T} (GeV/c)","Centrality","dN^{Jets}/d#delta p_{T}");
      AddHistogram2D<TH2D>("hDeltaPtTR", "Background fluctuations #delta p_{T} (TR, cone R=0.6)", "", 1201, -40.0, 40.0, fNumberOfCentralityBins, 0, 100,  "#delta p_{T} (GeV/c)","Centrality","dN^{Jets}/d#delta p_{T}");
      AddHistogram2D<TH2D>("hDeltaPtRC", "Background fluctuations #delta p_{T} (RC)", "", 1201, -40.0, 40.0, fNumberOfCentralityBins, 0, 100,  "#delta p_{T} (GeV/c)","Centrality","dN^{Jets}/d#delta p_{T}");


      // ########## Profiles for background means vs. centrality
      AddHistogram1D<TProfile>("hKTMeanBackgroundPbPb", "KT background mean (PbPb approach w/o ghosts)", "", 100, 0, 100, "Centrality", "#rho mean");
      AddHistogram1D<TProfile>("hKTMeanBackgroundPbPbWithGhosts", "KT background mean (PbPb approach)", "", 100, 0, 100, "Centrality", "#rho mean");
      AddHistogram1D<TProfile>("hKTMeanBackgroundCMS", "KT background mean (CMS approach)", "", 100, 0, 100, "Centrality", "#rho mean");
      AddHistogram1D<TProfile>("hKTMeanBackgroundMean", "KT background mean (Mean approach)", "",  100, 0, 100, "Centrality", "#rho mean");
      AddHistogram1D<TProfile>("hKTMeanBackgroundTPC", "KT background mean (Track-like approach)", "", 100, 0, 100, "Centrality", "#rho mean");
      AddHistogram1D<TProfile>("hTRMeanBackground", "TR background mean", "", 100, 0, 100, "Centrality", "#rho mean");
      AddHistogram1D<TProfile>("hRCMeanBackground", "RC background mean", "", 100, 0, 100, "Centrality", "#rho mean");
    }
  }


  // NOTE: Track & Cluster & QA histograms
  if (fAnalyzeQA)
  {
    AddHistogram1D<TH1D>("hNumberEvents", "Number of events (0 = before, 1 = after vertex cuts)", "", 2, 0, 2, "#Delta z(cm)","N^{Events}/cut");
    AddHistogram1D<TH1D>("hVertexX", "X distribution of the vertex", "", 2000, -1., 1., "#Delta x(cm)","dN^{Events}/dx");
    AddHistogram1D<TH1D>("hVertexY", "Y distribution of the vertex", "", 2000, -1., 1., "#Delta y(cm)","dN^{Events}/dy");
    AddHistogram2D<TH2D>("hVertexXY", "XY distribution of the vertex", "COLZ", 500, -1., 1., 500, -1., 1.,"#Delta x(cm)", "#Delta y(cm)","dN^{Events}/dxdy");
    AddHistogram1D<TH1D>("hVertexZ", "Z distribution of the vertex", "", 100, -10., 10., "#Delta z(cm)","dN^{Events}/dz");
    AddHistogram1D<TH1D>("hVertexR", "R distribution of the vertex", "", 100, 0., 1., "#Delta r(cm)","dN^{Events}/dr");
    AddHistogram1D<TH1D>("hCentralityV0M", "Centrality distribution V0M", "", 100, 0., 100., "Centrality","dN^{Events}");
    AddHistogram1D<TH1D>("hCentralityV0A", "Centrality distribution V0A", "", 100, 0., 100., "Centrality","dN^{Events}");
    AddHistogram1D<TH1D>("hCentralityV0C", "Centrality distribution V0C", "", 100, 0., 100., "Centrality","dN^{Events}");

    AddHistogram2D<TH2D>("hTrackCountAcc", "Number of tracks in acceptance vs. centrality", "LEGO2", 750, 0., 750., 100, 0, 100, "N tracks","Centrality", "dN^{Events}/dN^{Tracks}");
    AddHistogram1D<TH1D>("hTrackPt", "Tracks p_{T} distribution", "", 1000, 0., 250., "p_{T} (GeV/c)","dN^{Tracks}/dp_{T}");
    AddHistogram1D<TH1D>("hTrackPtNegEta", "Tracks p_{T} distribution (negative #eta)", "", 1000, 0., 250., "p_{T} (GeV/c)","dN^{Tracks}/dp_{T}");
    AddHistogram1D<TH1D>("hTrackPtPosEta", "Tracks p_{T} distribution (positive #eta)", "", 1000, 0., 250., "p_{T} (GeV/c)","dN^{Tracks}/dp_{T}");      
    AddHistogram1D<TH1D>("hTrackCharge", "Charge", "", 11, -5, 5, "Charge (e)","dN^{Tracks}/dq");
    AddHistogram1D<TH1D>("hTrackPhi", "Track #phi distribution", "", 360, 0, TMath::TwoPi(), "#phi","dN^{Tracks}/d#phi");
    AddHistogram2D<TH2D>("hTrackPhiEta", "Track angular distribution", "LEGO2", 100, 0., 2*TMath::Pi(),100, -2.5, 2.5, "#phi","#eta","dN^{Tracks}/(d#phi d#eta)");

    AddHistogram2D<TH2D>("hTrackPhiPtCut", "Track #phi distribution for different pT cuts", "LEGO2", 360, 0, TMath::TwoPi(), 20, 0, 20, "#phi", "p_{T} lower cut", "dN^{Tracks}/d#phi dp_{T}");      
    AddHistogram2D<TH2D>("hTrackPhiLabel", "Track #phi distribution in different labels", "LEGO2", 360, 0, TMath::TwoPi(), 3, 0, 3, "#phi", "Label", "dN^{Tracks}/d#phi");
    AddHistogram1D<TH1D>("hTrackEta", "Track #eta distribution", "", 180, -fTrackEtaWindow, +fTrackEtaWindow, "#eta","dN^{Tracks}/d#eta");
    if (fAnalyzeJets)
    {
      // ######## Jet QA
      AddHistogram1D<TH1D>("hJetArea", "Jets area distribution", "", 200, 0., 2., "Area","dN^{Jets}/dA");
      AddHistogram2D<TH2D>("hJetAreaVsPt", "Jets area vs. p_{T} distribution", "COLZ", 100, 0., 2., 200, 0., 200., "Area", "p_{T}", "dN^{Jets}/dA dp_{T}");
      AddHistogram2D<TH2D>("hJetPhiEta", "Jets angular distribution", "LEGO2", 360, 0., 2*TMath::Pi(),100, -0.6, 0.6, "#phi","#eta","dN^{Jets}/(d#phi d#eta)");
      AddHistogram2D<TH2D>("hJetPtVsConstituentCount", "Jets number of constituents vs. jet p_{T}", "COLZ", 400, 0., 200., 100, 0., 100., "p_{T}","N^{Tracks}","dN^{Jets}/(dp_{T} dN^{tracks})");
    }
  }


  // NOTE: Pythia histograms
  if (fAnalyzePythia)
  {
    AddHistogram1D<TH1D>("hPythiaPtHard", "Pythia p_{T} hard distribution", "", 2000, 0, 400, "p_{T} hard","dN^{Events}/dp_{T,hard}");
    AddHistogram1D<TProfile>("hPythiaXSection", "Pythia cross section distribution", "", fNumPtHardBins+2, -1, fNumPtHardBins+1, "p_{T} hard bin","dN^{Events}/dp_{T,hard}");
    AddHistogram1D<TH1D>("hPythiaNTrials", "Pythia trials (no correction for manual cuts)", "", fNumPtHardBins+2, -1, fNumPtHardBins+1, "p_{T} hard bin", "Trials");
  }

  // register Histograms
  for (Int_t i = 0; i < fHistCount; i++)
  {
    fOutputList->Add(fHistList->At(i));
  }
  
  PostData(1,fOutputList); // important for merging

}

//________________________________________________________________________
AliAnalysisTaskChargedJetsPA::AliAnalysisTaskChargedJetsPA(const char *name, const char* trackArrayName, const char* jetArrayName, const char* backgroundJetArrayName) : AliAnalysisTaskSE(name), fOutputList(0), fAnalyzeJets(1), fAnalyzeQA(1), fAnalyzeBackground(1), fAnalyzeDeprecatedBackgrounds(1), fAnalyzePythia(0), fHasTracks(0), fHasJets(0), fHasBackgroundJets(0), fIsMC(0), fJetArray(0), fTrackArray(0), fBackgroundJetArray(0), fJetArrayName(0), fTrackArrayName(0), fBackgroundJetArrayName(0), fNumPtHardBins(11), fUsePtHardBin(-1), fRandConeRadius(0.4), fSignalJetRadius(0.4), fBackgroundJetRadius(0.4), fTRBackgroundConeRadius(0.6), fNumberRandCones(8), fNumberExcludedJets(-1), fDijetMaxAngleDeviation(10.0), fPhysicalJetRadius(0.6), fSignalJetEtaWindow(0.5), fBackgroundJetEtaWindow(0.5), fTrackEtaWindow(0.9), fVertexWindow(10.0), fVertexMaxR(1.0), fMinTrackPt(0.150), fMinJetPt(1.0), fMinJetArea(0.5), fMinBackgroundJetPt(0.0), fMinDijetLeadingPt(10.0), fNumberOfCentralityBins(100), fCentralityType("V0A"), fFirstLeadingJet(0), fSecondLeadingJet(0), fNumberSignalJets(0), fCrossSection(0.0), fTrials(0.0),  fRandom(0), fHelperClass(0), fInitialized(0), fTaskInstanceCounter(0), fHistList(0), fHistCount(0), fIsDEBUG(0)
{
  #ifdef DEBUGMODE
    AliInfo("Calling constructor.");
  #endif

  // Every instance of this task gets his own number
  static Int_t instance = 0;
  fTaskInstanceCounter = instance;
  instance++;

  fTrackArrayName = new TString(trackArrayName);
  if (fTrackArrayName->Contains("MCParticles") || fTrackArrayName->Contains("mcparticles")) //TODO: Not working for now
    fIsMC = kTRUE;

  fJetArrayName = new TString(jetArrayName);
  if (strcmp(fJetArrayName->Data(),"") == 0)
    fAnalyzeJets = kFALSE;
  else
    fAnalyzeJets = kTRUE;
    
  fBackgroundJetArrayName = new TString(backgroundJetArrayName);
  if (strcmp(fBackgroundJetArrayName->Data(),"") == 0)
    fAnalyzeBackground = kFALSE;
  else
    fAnalyzeBackground = kTRUE;

  DefineOutput(1, TList::Class());
 
  fHistList = new TList();

  for(Int_t i=0;i<1024;i++)
    fSignalJets[i] = NULL;

  #ifdef DEBUGMODE
    AliInfo("Constructor done.");
  #endif
  
}

//________________________________________________________________________
inline Double_t AliAnalysisTaskChargedJetsPA::GetConePt(Double_t eta, Double_t phi, Double_t radius)
{
  Double_t tmpConePt = 0.0;

  for (Int_t i = 0; i < fTrackArray->GetEntries(); i++)
  {
    AliVTrack* tmpTrack = static_cast<AliVTrack*>(fTrackArray->At(i));
    if (IsTrackInAcceptance(tmpTrack))
      if(IsTrackInCone(tmpTrack, eta, phi, radius))
        tmpConePt = tmpConePt + tmpTrack->Pt();
  }
  return tmpConePt;
}


//________________________________________________________________________
inline Double_t AliAnalysisTaskChargedJetsPA::GetPtHard()
{
  Double_t tmpPtHard = -1.0;

  if (!MCEvent())
    AliError("MCEvent not accessible although demanded!");
  else
  {
    AliGenPythiaEventHeader* pythiaHeader = dynamic_cast<AliGenPythiaEventHeader*>(MCEvent()->GenEventHeader());
    if (!pythiaHeader)
    {
      // Check if AOD
      AliAODMCHeader* aodMCH = dynamic_cast<AliAODMCHeader*>(InputEvent()->FindListObject(AliAODMCHeader::StdBranchName()));

      for(UInt_t i = 0;i<aodMCH->GetNCocktailHeaders();i++)
      {
        pythiaHeader = dynamic_cast<AliGenPythiaEventHeader*>(aodMCH->GetCocktailHeader(i));
        if(pythiaHeader) break;
      }
      if(!pythiaHeader)
        AliError("Pythia Header not accessible!");
      else
        tmpPtHard = pythiaHeader->GetPtHard();
    }
    else
      tmpPtHard = pythiaHeader->GetPtHard();
  }
  return tmpPtHard;
}


//________________________________________________________________________
inline Int_t AliAnalysisTaskChargedJetsPA::GetPtHardBin()
{
  // ########## PT HARD BIN EDGES
  const Int_t kPtHardLowerEdges[] =  { 0, 5,11,21,36,57, 84,117,152,191,234};
  const Int_t kPtHardHigherEdges[] = { 5,11,21,36,57,84,117,152,191,234,1000000};

  Int_t tmpPtHardBin = 0;
  Double_t tmpPtHard = GetPtHard();
 
  for (tmpPtHardBin = 0; tmpPtHardBin <= fNumPtHardBins; tmpPtHardBin++)
    if (tmpPtHard >= kPtHardLowerEdges[tmpPtHardBin] && tmpPtHard < kPtHardHigherEdges[tmpPtHardBin])
      break;

  return tmpPtHardBin;
}


//________________________________________________________________________
inline Bool_t AliAnalysisTaskChargedJetsPA::IsTrackInCone(AliVTrack* track, Double_t eta, Double_t phi, Double_t radius)
{
  // This is to use a full cone in phi even at the edges of phi (2pi -> 0) (0 -> 2pi)
  Double_t trackPhi = 0.0;
  if (track->Phi() > (TMath::TwoPi() - (radius-phi)))
    trackPhi = track->Phi() - TMath::TwoPi();
  else if (track->Phi() < (phi+radius - TMath::TwoPi()))
    trackPhi = track->Phi() + TMath::TwoPi();
  else
    trackPhi = track->Phi();
  
  if ( TMath::Abs(trackPhi-phi)*TMath::Abs(trackPhi-phi) + TMath::Abs(track->Eta()-eta)*TMath::Abs(track->Eta()-eta) <= radius*radius)
    return kTRUE;
  
  return kFALSE;
}

//________________________________________________________________________
inline Bool_t AliAnalysisTaskChargedJetsPA::IsTrackInAcceptance(AliVParticle* track)
{
  if (track != 0)
    if (TMath::Abs(track->Eta()) <= fTrackEtaWindow)
      if (track->Pt() >= fMinTrackPt)
        return kTRUE;

  return kFALSE;
}

//________________________________________________________________________
inline Bool_t AliAnalysisTaskChargedJetsPA::IsTrackInJet(AliEmcalJet* jet, Int_t trackIndex)
{
  for (Int_t i = 0; i < jet->GetNumberOfTracks(); ++i)
  {
    Int_t jetTrack = jet->TrackAt(i);
    if (jetTrack == trackIndex)
      return kTRUE;
  }
  return kFALSE;
}

//________________________________________________________________________
inline Bool_t AliAnalysisTaskChargedJetsPA::IsJetOverlapping(AliEmcalJet* jet1, AliEmcalJet* jet2)
{
  for (Int_t i = 0; i < jet1->GetNumberOfTracks(); ++i)
  {
    Int_t jet1Track = jet1->TrackAt(i);
    for (Int_t j = 0; j < jet2->GetNumberOfTracks(); ++j)
    {
      Int_t jet2Track = jet2->TrackAt(j);
      if (jet1Track == jet2Track)
        return kTRUE;
    }
  }
  return kFALSE;
}


//________________________________________________________________________
inline Bool_t AliAnalysisTaskChargedJetsPA::IsBackgroundJetInAcceptance(AliEmcalJet *jet)
{   
  if (jet != 0)
    if (TMath::Abs(jet->Eta()) <= fBackgroundJetEtaWindow)
      if (jet->Pt() >= fMinBackgroundJetPt)
        return kTRUE;

  return kFALSE;
}

//________________________________________________________________________
inline Bool_t AliAnalysisTaskChargedJetsPA::IsSignalJetInAcceptance(AliEmcalJet *jet)
{   
  if (jet != 0)
    if (TMath::Abs(jet->Eta()) <= fSignalJetEtaWindow)
      if (jet->Pt() >= fMinJetPt)
        if (jet->Area() >= fMinJetArea)
          return kTRUE;
  return kFALSE;
}

//________________________________________________________________________
inline Bool_t AliAnalysisTaskChargedJetsPA::IsDijet(AliEmcalJet *jet1, AliEmcalJet *jet2)
{   
  // Output from GetDeltaPhi is < pi in any case
  if ((jet1 != 0) && (jet2 != 0))
    if((TMath::Pi() - GetDeltaPhi(jet1->Phi(),jet2->Phi())) < fDijetMaxAngleDeviation)
      if ((jet1->Pt() > fMinDijetLeadingPt) && (jet2->Pt() > fMinDijetLeadingPt)) //TODO: Introduce recoil jet?
        return kTRUE;

  return kFALSE;
}

//________________________________________________________________________
void AliAnalysisTaskChargedJetsPA::ExecOnce()
{
  #ifdef DEBUGMODE
    AliInfo("Starting ExecOnce.");
  #endif
  fInitialized = kTRUE;

  // Check for track array
  if (strcmp(fTrackArrayName->Data(), "") != 0)
  {
    fTrackArray = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(fTrackArrayName->Data()));
    fHasTracks = kTRUE;
    if (!fTrackArray) 
    {
      AliWarning(Form("%s: Could not retrieve tracks %s! This is OK, if tracks are not demanded.", GetName(), fTrackArrayName->Data())); 
      fHasTracks = kFALSE;
    } 
    else
    {
      TClass *cl = fTrackArray->GetClass();
      if (!cl->GetBaseClass("AliVParticle"))
      {
      	AliError(Form("%s: Collection %s does not contain AliVParticle objects!", GetName(), fTrackArrayName->Data())); 
      	fTrackArray = 0;
        fHasTracks = kFALSE;
      }
    }
  }

  // Check for jet array
  if (strcmp(fJetArrayName->Data(), "") != 0)
  {
    fJetArray = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(fJetArrayName->Data()));
    fHasJets = kTRUE;

    if (!fJetArray) 
    {
      AliWarning(Form("%s: Could not retrieve jets %s! This is OK, if jets are not demanded.", GetName(), fJetArrayName->Data())); 
      fHasJets = kFALSE;
    } 
    else
    {
      if (!fJetArray->GetClass()->GetBaseClass("AliEmcalJet")) 
      {
        AliError(Form("%s: Collection %s does not contain AliEmcalJet objects!", GetName(), fJetArrayName->Data())); 
        fJetArray = 0;
        fHasJets = kFALSE;
      }
    }
  }

  // Check for background object
  if (strcmp(fBackgroundJetArrayName->Data(), "") != 0)
  {
    fBackgroundJetArray = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(fBackgroundJetArrayName->Data()));
    fHasBackgroundJets = kTRUE;
    if (!fBackgroundJetArray)
    {
      AliInfo(Form("%s: Could not retrieve background jets %s! This is OK, if background is not demanded.", GetName(), fBackgroundJetArrayName->Data())); 
      fHasBackgroundJets = kFALSE;
    }
  }

  // Look, if initialization is OK
  if (!fHasTracks && fAnalyzeBackground)
  {
    AliError(Form("%s: Tracks NOT successfully casted although demanded! Deactivating background analysis",GetName()));
    fAnalyzeBackground = kFALSE;
  }
  if ((!fHasJets && fAnalyzeJets) || (!fHasJets && fAnalyzeBackground))
  {
    AliError(Form("%s: Jets NOT successfully casted although demanded!  Deactivating jet- and background analysis",GetName()));
    fAnalyzeJets = kFALSE;
    fAnalyzeBackground = kFALSE;
  }
  if (!fHasBackgroundJets && fAnalyzeBackground)
  {
    AliError(Form("%s: Background NOT successfully casted although demanded!  Deactivating background analysis",GetName()));
    fAnalyzeBackground = kFALSE;
  }

  // Initialize helper class (for vertex selection)
  fHelperClass = new AliAnalysisUtils();

  // Histogram init
  Init();

  #ifdef DEBUGMODE
    AliInfo("ExecOnce done.");
  #endif

}

//________________________________________________________________________
void AliAnalysisTaskChargedJetsPA::GetSignalJets()
{
  // Reset vars
  fFirstLeadingJet = NULL;
  fSecondLeadingJet = NULL;
  fNumberSignalJets = 0;

  TList tmpJets;
  for (Int_t i = 0; i < fJetArray->GetEntries(); i++)
  {
    AliEmcalJet* jet = static_cast<AliEmcalJet*>(fJetArray->At(i));
    if (!jet)
    {
      AliError(Form("%s: Could not receive jet %d", GetName(), i));
      continue;
    }
    if (!IsSignalJetInAcceptance(jet))
      continue;

    for (Int_t j = 0; j <= tmpJets.GetEntries(); j++)
    {
      if (j>tmpJets.GetEntries()-1) // When passed last item add the jet at the end
      {
        tmpJets.Add(jet);
        break;
      }

      AliEmcalJet* listJet = static_cast<AliEmcalJet*>(tmpJets.At(j));
     
      if(jet->Pt() < listJet->Pt()) // Insert jet before that one in list if pt smaller
      {
        tmpJets.AddAt(jet, j);
        break;
      }
    }
  }

  for (Int_t i = 0; i < tmpJets.GetEntries(); i++)
  {
    AliEmcalJet* jet = static_cast<AliEmcalJet*>(tmpJets.At(i));
    fSignalJets[fNumberSignalJets] = jet;
    fNumberSignalJets++;
  }
  
  if (fNumberSignalJets > 0)
    fFirstLeadingJet  = static_cast<AliEmcalJet*>(tmpJets.At(0));
  if (fNumberSignalJets > 1)
    fSecondLeadingJet = static_cast<AliEmcalJet*>(tmpJets.At(1));

}

//________________________________________________________________________
Int_t AliAnalysisTaskChargedJetsPA::GetLeadingJets(TClonesArray* jetArray, Int_t* jetIDArray, Bool_t isSignalJets)
{
// Writes first two leading jets into already registered array jetIDArray

  if (!jetArray)
  {
    AliError("Could not get the jet array to get leading jets from it!");
    return 0;
  }

  Float_t maxJetPts[] = {0, 0};
  jetIDArray[0] = -1;
  jetIDArray[1] = -1;

  Int_t jetCount = jetArray->GetEntries();
  Int_t jetCountAccepted = 0;

  for (Int_t i = 0; i < jetCount; i++)
  {
    AliEmcalJet* jet = static_cast<AliEmcalJet*>(jetArray->At(i));
    if (!jet) 
    {
      AliError(Form("%s: Could not receive jet %d", GetName(), i));
      continue;
    }

    if(isSignalJets)
    {
      if (!IsSignalJetInAcceptance(jet)) continue;
    }
    else
    {
      if (!IsBackgroundJetInAcceptance(jet)) continue;
    }    

    if (jet->Pt() > maxJetPts[0]) 
    {
      maxJetPts[1] = maxJetPts[0];
      jetIDArray[1] = jetIDArray[0];
      maxJetPts[0] = jet->Pt();
      jetIDArray[0] = i;
    }
    else if (jet->Pt() > maxJetPts[1]) 
    {
      maxJetPts[1] = jet->Pt();
      jetIDArray[1] = i;
    }
    jetCountAccepted++;
  }
  return jetCountAccepted;
}


//________________________________________________________________________
Double_t AliAnalysisTaskChargedJetsPA::GetCorrectedJetPt(AliEmcalJet* jet, Double_t background)
{
  #ifdef DEBUGMODE
    AliInfo("Getting corrected jet spectra.");
  #endif

  if(!jet)
  {
    AliError("Jet pointer passed to GetCorrectedJet() not valid!");
    return -1.0;
  }

  Double_t correctedPt = -1.0;

  // if the passed background is not valid, do not subtract it
  if(background < 0)
    background = 0;

  // Subtract background
  correctedPt = jet->Pt() - background * jet->Area();

  #ifdef DEBUGMODE
    AliInfo("Got corrected jet spectra.");
  #endif 

  return correctedPt;
}



//________________________________________________________________________
void AliAnalysisTaskChargedJetsPA::GetDeltaPt(Double_t& deltaPt, Double_t rho, Bool_t leadingJetExclusion)
{
  #ifdef DEBUGMODE
    AliInfo("Getting Delta Pt.");
  #endif

  // Define an invalid delta pt
  deltaPt = -10000.0;

  // Define eta range
  Double_t etaMin, etaMax;
  etaMin = -(fTrackEtaWindow-fRandConeRadius);
  etaMax = +(fTrackEtaWindow-fRandConeRadius);

  // Define random cone
  Bool_t coneValid = kTRUE;
  Double_t tmpRandConeEta = etaMin + fRandom->Rndm()*(etaMax-etaMin);
  Double_t tmpRandConePhi = fRandom->Rndm()*TMath::TwoPi();

  // if there is a jet, check for overlap if demanded
  if(leadingJetExclusion)
  {
    for (Int_t i = 0; i<fNumberSignalJets; i++)
    {
      AliEmcalJet* tmpJet = fSignalJets[i];

      Double_t excludedJetPhi = tmpJet->Phi();
      Double_t excludedJetEta = tmpJet->Eta();
      Double_t tmpDeltaPhi = GetDeltaPhi(tmpRandConePhi, excludedJetPhi);

      // Check, if cone has overlap with jet
      if ( tmpDeltaPhi*tmpDeltaPhi + TMath::Abs(tmpRandConeEta-excludedJetEta)*TMath::Abs(tmpRandConeEta-excludedJetEta) <= fRandConeRadius*fRandConeRadius)
      {
        // Define probability to exclude the RC
        Double_t probability = 1 - (fNumberSignalJets-1)/fNumberSignalJets;

        // Only exclude cone with a given probability
        if (fRandom->Rndm()<=probability)
        {
          coneValid = kFALSE;
          break;
        }
      }
    }
  }

  // Get the cones' pt and calculate delta pt
  if (coneValid)
    deltaPt = GetConePt(tmpRandConeEta,tmpRandConePhi,fRandConeRadius) - (rho*fRandConeRadius*fRandConeRadius*TMath::Pi());

  #ifdef DEBUGMODE
    AliInfo("Got Delta Pt.");
  #endif
}

//________________________________________________________________________
void AliAnalysisTaskChargedJetsPA::GetKTBackgroundDensity(Int_t numberExcludeLeadingJets, Double_t& rhoPbPb, Double_t& rhoPbPbWithGhosts, Double_t& rhoCMS, Double_t& rhoImprovedCMS, Double_t& rhoMean, Double_t& rhoTrackLike)
{
  #ifdef DEBUGMODE
    AliInfo("Getting KT background density.");
  #endif

  static Double_t tmpRhoPbPb[1024];
  static Double_t tmpRhoPbPbWithGhosts[1024];
  static Double_t tmpRhoMean[1024];
  static Double_t tmpRhoCMS[1024];
  static Double_t tmpRhoImprovedCMS[1024];
  Double_t tmpCoveredArea = 0.0;
  Double_t tmpSummedArea = 0.0;
  Double_t tmpPtTrackLike = 0.0;
  Double_t tmpAreaTrackLike = 0.0;

  // Setting invalid values
  rhoPbPb = 0.0;
  rhoPbPbWithGhosts = 0.0;
  rhoCMS = 0.0;
  rhoImprovedCMS = 0.0;
  rhoMean = 0.0;
  rhoTrackLike = 0.0;

  Int_t rhoPbPbJetCount = 0;
  Int_t rhoPbPbWithGhostsJetCount = 0;
  Int_t rhoCMSJetCount = 0;
  Int_t rhoImprovedCMSJetCount = 0;
  Int_t rhoMeanJetCount = 0;


  // Find 2 leading KT jets for the original PbPb approach
  Int_t leadingKTJets[]   = {-1, -1};
  GetLeadingJets(fBackgroundJetArray, &leadingKTJets[0], kFALSE);

  // Exclude UP TO numberExcludeLeadingJets
  if(numberExcludeLeadingJets==-1)
    numberExcludeLeadingJets = fNumberSignalJets;
  if (fNumberSignalJets < numberExcludeLeadingJets)
    numberExcludeLeadingJets = fNumberSignalJets;

  for (Int_t i = 0; i < fBackgroundJetArray->GetEntries(); i++)
  {
    AliEmcalJet* backgroundJet = static_cast<AliEmcalJet*>(fBackgroundJetArray->At(i));

    if (!backgroundJet)
    {
      AliError(Form("%s: Could not receive jet %d", GetName(), i));
      continue;
    } 

    // Search for overlap with signal jets
    Bool_t isOverlapping = kFALSE;
    for(Int_t j=0;j<numberExcludeLeadingJets;j++)
    {
      AliEmcalJet* signalJet = fSignalJets[j];
     
      if(IsJetOverlapping(signalJet, backgroundJet))
      {
        isOverlapping = kTRUE;
        break;
      }
    }

    tmpSummedArea += backgroundJet->Area();
    if(backgroundJet->Pt() > 0.150)
      tmpCoveredArea += backgroundJet->Area();

    if (!IsBackgroundJetInAcceptance(backgroundJet))
      continue;

    // PbPb approach (take ghosts into account)
    if ((i != leadingKTJets[0]) && (i != leadingKTJets[1]))
    {
      tmpRhoPbPbWithGhosts[rhoPbPbWithGhostsJetCount] = backgroundJet->Pt() / backgroundJet->Area();
      rhoPbPbWithGhostsJetCount++;
    }

    if(backgroundJet->Pt() > 0.150)
    {
      // CMS approach: don't take ghosts into acount
      tmpRhoCMS[rhoCMSJetCount] = backgroundJet->Pt() / backgroundJet->Area();
      rhoCMSJetCount++;

      // Improved CMS approach: like CMS but excluding signal
      if(!isOverlapping)
      {
        tmpRhoImprovedCMS[rhoImprovedCMSJetCount] = backgroundJet->Pt() / backgroundJet->Area();
        rhoImprovedCMSJetCount++;
      }

      // PbPb w/o ghosts approach (just neglect ghosts)
      if ((i != leadingKTJets[0]) && (i != leadingKTJets[1]))
      {  
        tmpRhoPbPb[rhoPbPbJetCount] = backgroundJet->Pt() / backgroundJet->Area();
        rhoPbPbJetCount++;
      }
    }

    // (no overlap with signal jets)
    if(!isOverlapping)
    {
      // Mean approach
      tmpRhoMean[rhoMeanJetCount] = backgroundJet->Pt() / backgroundJet->Area();
      rhoMeanJetCount++;
      
      // Track like approach approach
      tmpPtTrackLike += backgroundJet->Pt();
      tmpAreaTrackLike += backgroundJet->Area();
    }

  }

  if (tmpAreaTrackLike > 0)
    rhoTrackLike = tmpPtTrackLike/tmpAreaTrackLike;
  if (rhoPbPbJetCount > 0)
    rhoPbPb = TMath::Median(rhoPbPbJetCount, tmpRhoPbPb);
  if (rhoPbPbWithGhostsJetCount > 0)
    rhoPbPbWithGhosts = TMath::Median(rhoPbPbWithGhostsJetCount, tmpRhoPbPbWithGhosts);
  if (rhoCMSJetCount > 0)
    rhoCMS = TMath::Median(rhoCMSJetCount, tmpRhoCMS) * tmpCoveredArea/tmpSummedArea;
  if (rhoImprovedCMSJetCount > 0)
    rhoImprovedCMS = TMath::Median(rhoImprovedCMSJetCount, tmpRhoImprovedCMS) * tmpCoveredArea/tmpSummedArea;
  if (rhoMeanJetCount > 0)
    rhoMean = TMath::Mean(rhoMeanJetCount, tmpRhoMean);

  #ifdef DEBUGMODE
    AliInfo("Got KT background density.");
  #endif
}



//________________________________________________________________________
Int_t AliAnalysisTaskChargedJetsPA::GetRCBackgroundDensity(Int_t numberExcludeLeadingJets, Double_t& rhoMean, Double_t& rhoMedian, Double_t etaMin, Double_t etaMax, Int_t numberRandCones)
{
  #ifdef DEBUGMODE
    AliInfo("Getting RC background density.");
  #endif

  if(numberRandCones == 0)
    numberRandCones = fNumberRandCones;

  std::vector<AliEmcalJet> tmpCones(numberRandCones);

  // Setting invalid values
  rhoMean = 0.0;
  rhoMedian = 0.0;

  // Exclude UP TO numberExcludeLeadingJets
  if(numberExcludeLeadingJets==-1)
    numberExcludeLeadingJets = fNumberSignalJets;
  if (fNumberSignalJets < numberExcludeLeadingJets)
    numberExcludeLeadingJets = fNumberSignalJets;

  // Search given amount of RCs
  Int_t numAcceptedRCs = 0;
  for(Int_t i=0;i<numberRandCones;i++)
  {
    Double_t tmpRandConeEta = 0.0;
    Double_t tmpRandConePhi = 0.0;

    // Search random cone in acceptance with no overlap with already excluded jets (leading jets and random cones)

    // Check if etaMin/etaMax is given correctly
    if(etaMin < -fSignalJetEtaWindow)
      etaMin = -fSignalJetEtaWindow;
    if(etaMax > fSignalJetEtaWindow)
      etaMax = fSignalJetEtaWindow;

    // Set the random cone position
    if ((etaMin == 0) && (etaMax == 0))
      tmpRandConeEta = (fTrackEtaWindow-fRandConeRadius)*(2.0*fRandom->Rndm()-1.0); // full RC is in acceptance
    else
      tmpRandConeEta = etaMin + fRandom->Rndm()*(etaMax-etaMin);

    tmpRandConePhi = fRandom->Rndm()*TMath::TwoPi();

    // Exclude signal jets
    Bool_t coneValid = kFALSE;
    for(Int_t j=0;j<numberExcludeLeadingJets;j++)
    {
      AliEmcalJet* signalJet = fSignalJets[j];

      Double_t tmpDeltaPhi = GetDeltaPhi(tmpRandConePhi, signalJet->Phi());
      
      if ( tmpDeltaPhi*tmpDeltaPhi + TMath::Abs(signalJet->Eta()-tmpRandConeEta)*TMath::Abs(signalJet->Eta()-tmpRandConeEta) <= (fRandConeRadius+fPhysicalJetRadius)*(fRandConeRadius+fPhysicalJetRadius))
      {
        coneValid = kFALSE;
        break;
      }
    }

    // RC is accepted, so save it
    if(coneValid)
    {
      AliEmcalJet tmpJet(GetConePt(tmpRandConeEta, tmpRandConePhi, fRandConeRadius), tmpRandConeEta, tmpRandConePhi, 0.0);
      tmpCones[numAcceptedRCs] = tmpJet;
      numAcceptedRCs++;
    }
  }

  // Calculate Rho and the mean from the RCs (no excluded jets are considered!)
  if(numAcceptedRCs > 0)
  {
    std::vector<Double_t> tmpRho(numAcceptedRCs);
    for (Int_t i=0; i<numAcceptedRCs;i++)
      tmpRho[i] = tmpCones[i].Pt()/(fRandConeRadius*fRandConeRadius*TMath::Pi());

    rhoMean = TMath::Mean(tmpRho.begin(), tmpRho.end());
    rhoMedian = 0.0; // NOT IMPLEMENTED because TMath::Median is not working with iterators
  }
    
  #ifdef DEBUGMODE
    AliInfo("Got RC background density.");
  #endif
  return numAcceptedRCs;
}

//________________________________________________________________________
void AliAnalysisTaskChargedJetsPA::GetTRBackgroundDensity(Int_t numberExcludeLeadingJets, Double_t& rhoNoExclusion, Double_t& rhoConeExclusion02, Double_t& rhoConeExclusion04, Double_t& rhoConeExclusion06, Double_t& rhoConeExclusion08, Double_t& rhoExactExclusion)
{
  #ifdef DEBUGMODE
    AliInfo("Getting TR background density.");
  #endif

  Double_t summedTracksPtCone02 = 0.0;
  Double_t summedTracksPtCone04 = 0.0;
  Double_t summedTracksPtCone06 = 0.0;
  Double_t summedTracksPtCone08 = 0.0;
  Double_t summedTracksPtWithinJets = 0.0;
  Double_t summedTracksPt = 0.0;
  
  // Setting invalid values
  rhoNoExclusion = 0.0;
  rhoConeExclusion02 = 0.0;
  rhoConeExclusion04 = 0.0;
  rhoConeExclusion06 = 0.0;
  rhoConeExclusion08 = 0.0; 
  rhoExactExclusion  = 0.0;

  // Exclude UP TO numberExcludeLeadingJets
  if(numberExcludeLeadingJets==-1)
    numberExcludeLeadingJets = fNumberSignalJets;
  if (fNumberSignalJets < numberExcludeLeadingJets)
    numberExcludeLeadingJets = fNumberSignalJets;

  for (Int_t i = 0; i < fTrackArray->GetEntries(); i++)
  {
    AliVTrack* tmpTrack = static_cast<AliVTrack*>(fTrackArray->At(i));
    Bool_t trackWithinJet = kFALSE; Bool_t trackWithin02Cone = kFALSE; Bool_t trackWithin04Cone = kFALSE; Bool_t trackWithin06Cone = kFALSE; Bool_t trackWithin08Cone = kFALSE;

    if (IsTrackInAcceptance(tmpTrack))
    {
      // Check if tracks overlaps with jet
      for(Int_t j=0;j<numberExcludeLeadingJets;j++)
      {
        AliEmcalJet* signalJet = fSignalJets[j];

        // Exact jet exclusion
        if (IsTrackInJet(signalJet, i))
          trackWithinJet = kTRUE;

        // Cone exclusions
        if (IsTrackInCone(tmpTrack, signalJet->Eta(), signalJet->Phi(), 0.2))
        {
          trackWithin02Cone = kTRUE;
          trackWithin04Cone = kTRUE;
          trackWithin06Cone = kTRUE;
          trackWithin08Cone = kTRUE;
          break;
        }
        else if (IsTrackInCone(tmpTrack, signalJet->Eta(), signalJet->Phi(), 0.4))
        {
          trackWithin04Cone = kTRUE;
          trackWithin06Cone = kTRUE;
          trackWithin08Cone = kTRUE;
        }
        else if (IsTrackInCone(tmpTrack, signalJet->Eta(), signalJet->Phi(), 0.6))
        {
          trackWithin06Cone = kTRUE;
          trackWithin08Cone = kTRUE;
        }
        else if (IsTrackInCone(tmpTrack, signalJet->Eta(), signalJet->Phi(), 0.8))
        {
          trackWithin08Cone = kTRUE;
        }
      }

      if(!trackWithin08Cone)
      {
        summedTracksPtCone08 += tmpTrack->Pt();
      }
      if(!trackWithin06Cone)
      {
        summedTracksPtCone06 += tmpTrack->Pt();
      }
      if(!trackWithin04Cone)
      {
        summedTracksPtCone04 += tmpTrack->Pt();
      }
      if(!trackWithin02Cone)
      {
        summedTracksPtCone02 += tmpTrack->Pt();
      }
      if(!trackWithinJet)
      {
        summedTracksPtWithinJets += tmpTrack->Pt();
      }
      summedTracksPt += tmpTrack->Pt();

    }
  }

  // Calculate the correct area where the tracks were taking from

  Double_t tmpFullTPCArea = (2.0*fTrackEtaWindow) * TMath::TwoPi();
  Double_t tmpAreaCone02     = tmpFullTPCArea;
  Double_t tmpAreaCone04     = tmpFullTPCArea;
  Double_t tmpAreaCone06     = tmpFullTPCArea;
  Double_t tmpAreaCone08     = tmpFullTPCArea;
  Double_t tmpAreaWithinJets = tmpFullTPCArea;
  std::vector<Double_t> tmpEtas(numberExcludeLeadingJets);
  std::vector<Double_t> tmpPhis(numberExcludeLeadingJets);

  for(Int_t i=0;i<numberExcludeLeadingJets;i++)
  {
    AliEmcalJet* tmpJet = fSignalJets[i];
    tmpEtas[i] = tmpJet->Eta();
    tmpPhis[i] = tmpJet->Phi();
    tmpAreaWithinJets -= tmpJet->Area();
  }

  tmpAreaCone02 -= tmpFullTPCArea * MCGetOverlapMultipleCirclesRectancle(numberExcludeLeadingJets, tmpEtas, tmpPhis, 0.2, -fTrackEtaWindow, +fTrackEtaWindow, 0., TMath::TwoPi());
  tmpAreaCone04 -= tmpFullTPCArea * MCGetOverlapMultipleCirclesRectancle(numberExcludeLeadingJets, tmpEtas, tmpPhis, 0.4, -fTrackEtaWindow, +fTrackEtaWindow, 0., TMath::TwoPi());
  tmpAreaCone06 -= tmpFullTPCArea * MCGetOverlapMultipleCirclesRectancle(numberExcludeLeadingJets, tmpEtas, tmpPhis, 0.6, -fTrackEtaWindow, +fTrackEtaWindow, 0., TMath::TwoPi());
  tmpAreaCone08 -= tmpFullTPCArea * MCGetOverlapMultipleCirclesRectancle(numberExcludeLeadingJets, tmpEtas, tmpPhis, 0.8, -fTrackEtaWindow, +fTrackEtaWindow, 0., TMath::TwoPi());
 
  rhoConeExclusion02 = summedTracksPtCone02/tmpAreaCone02;
  rhoConeExclusion04 = summedTracksPtCone04/tmpAreaCone04;
  rhoConeExclusion06 = summedTracksPtCone06/tmpAreaCone06;
  rhoConeExclusion08 = summedTracksPtCone08/tmpAreaCone08;
  rhoExactExclusion  = summedTracksPtWithinJets/tmpAreaWithinJets;
  rhoNoExclusion     = summedTracksPt/tmpFullTPCArea;


  #ifdef DEBUGMODE
    AliInfo("Got TR background density.");
  #endif
}

//________________________________________________________________________
void AliAnalysisTaskChargedJetsPA::GetTRBackgroundDensity(Int_t numberExcludeLeadingJets, Double_t& rhoMean, Double_t& area, AliEmcalJet* excludeJet1, AliEmcalJet* excludeJet2, Bool_t doSearchPerpendicular)
{
  #ifdef DEBUGMODE
    AliInfo("Getting TR background density.");
  #endif

  // Setting invalid values
  Double_t summedTracksPt = 0.0;
  rhoMean = 0.0;
  area = -1.0;

  Double_t tmpRadius = 0.0;
  if (doSearchPerpendicular)
    tmpRadius = 0.4*TMath::Pi(); // exclude 90 degrees around jets
  else
    tmpRadius = 0.8;
    
  numberExcludeLeadingJets = 2; // dijet is excluded here in any case



  if (!fTrackArray || !fJetArray)
  {
    AliError("Could not get the track/jet array to calculate track rho!");
    return;
  }

  Int_t   trackCount = fTrackArray->GetEntries();
  Int_t   trackCountAccepted = 0;
  for (Int_t i = 0; i < trackCount; i++)
  {
    AliVTrack* tmpTrack = static_cast<AliVTrack*>(fTrackArray->At(i));
    if (IsTrackInAcceptance(tmpTrack))
    {
      if (IsTrackInCone(tmpTrack, excludeJet1->Eta(), excludeJet1->Phi(), tmpRadius))
        continue;

      if (numberExcludeLeadingJets > 1)
        if (IsTrackInCone(tmpTrack, excludeJet2->Eta(), excludeJet2->Phi(), tmpRadius))
          continue;

        // Add track pt to array
        summedTracksPt = summedTracksPt + tmpTrack->Pt();
        trackCountAccepted++;
    }
  }

  if (trackCountAccepted > 0)
  {
    Double_t tmpArea = 2.0*fTrackEtaWindow*TMath::TwoPi()  - 2*(tmpRadius*tmpRadius * TMath::Pi()); //TPC area - excluding jet area
    rhoMean = summedTracksPt/tmpArea;
    area = tmpArea;
  }

  #ifdef DEBUGMODE
    AliInfo("Got TR background density.");
  #endif
}

//________________________________________________________________________
void AliAnalysisTaskChargedJetsPA::Calculate(AliVEvent* event)
{
  #ifdef DEBUGMODE
    AliInfo("Starting Calculate().");
  #endif
  ////////////////////// NOTE: initialization & casting

  // Check, if analysis should be done in pt hard bins
  if(fUsePtHardBin != -1)
    if(GetPtHardBin() != fUsePtHardBin)
      return;


  // Additional cuts
  FillHistogram("hNumberEvents", 0.5); // number of events before manual cuts

  if(!fHelperClass->IsVertexSelected2013pA(event))
    return;
 
  FillHistogram("hNumberEvents", 1.5); // number of events after manual cuts

  #ifdef DEBUGMODE
    AliInfo("Calculate()::Init done.");
  #endif

  ////////////////////// NOTE: Get Centrality, (Leading)Signal jets and Background

  // Get centrality
  AliCentrality* tmpCentrality = NULL;
  tmpCentrality = event->GetCentrality();
  Double_t centralityPercentile = 0.0;
  Double_t centralityPercentileV0A = 0.0;
  Double_t centralityPercentileV0C = 0.0;
  Double_t centralityPercentileV0M = 0.0;
  if (tmpCentrality != NULL)
  {
    centralityPercentile = tmpCentrality->GetCentralityPercentile(fCentralityType.Data());
    centralityPercentileV0A = tmpCentrality->GetCentralityPercentile("V0A");
    centralityPercentileV0C = tmpCentrality->GetCentralityPercentile("V0C");
    centralityPercentileV0M = tmpCentrality->GetCentralityPercentile("V0M");
  }

  if((centralityPercentile < 0.0) || (centralityPercentile > 100.0))
  {
    AliWarning(Form("Centrality value not valid (c=%E), setting to failsafe c=1.0.",centralityPercentile)); 
    centralityPercentile = 1.0;
  }
  // Get jets
  if (fAnalyzeBackground || fAnalyzeJets)
    GetSignalJets();

  // Get background estimates
  Double_t              backgroundKTImprovedCMS = -1.0;
  Double_t              backgroundDijet = -1.0;
  Double_t              backgroundDijetPerpendicular = -1.0;

  Double_t              backgroundKTPbPb = -1.0;
  Double_t              backgroundKTPbPbWithGhosts = -1.0;
  Double_t              backgroundKTCMS = -1.0;
  Double_t              backgroundKTMean = -1.0;
  Double_t              backgroundKTTrackLike = -1.0;
  Double_t              backgroundTRNoExcl = -1.0;
  Double_t              backgroundTRCone02 = -1.0;
  Double_t              backgroundTRCone04 = -1.0;
  Double_t              backgroundTRCone06 = -1.0;
  Double_t              backgroundTRCone08 = -1.0;
  Double_t              backgroundTRExact  = -1.0;
  Double_t              backgroundRC = -1.0;

  // Calculate background for different jet exclusions

  if (fAnalyzeBackground)
  {
    Double_t dummy = 0.0;

    GetKTBackgroundDensity    (fNumberExcludedJets, backgroundKTPbPb, backgroundKTPbPbWithGhosts, backgroundKTCMS, backgroundKTImprovedCMS, backgroundKTMean, backgroundKTTrackLike);
//    cout << "My task brings rho= " << backgroundKTImprovedCMS << endl; // DEBUG
    GetRCBackgroundDensity    (fNumberExcludedJets, backgroundRC, dummy);
    GetTRBackgroundDensity    (fNumberExcludedJets, backgroundTRNoExcl, backgroundTRCone02, backgroundTRCone04, backgroundTRCone06, backgroundTRCone08, backgroundTRExact);
  }

  #ifdef DEBUGMODE
    AliInfo("Calculate()::Centrality&SignalJets&Background-Calculation done.");
  #endif

  if (fAnalyzeQA)
  {
    FillHistogram("hVertexX",event->GetPrimaryVertex()->GetX());
    FillHistogram("hVertexY",event->GetPrimaryVertex()->GetY());
    FillHistogram("hVertexXY",event->GetPrimaryVertex()->GetX(), event->GetPrimaryVertex()->GetY());
    FillHistogram("hVertexZ",event->GetPrimaryVertex()->GetZ());
    FillHistogram("hVertexR",TMath::Sqrt(event->GetPrimaryVertex()->GetX()*event->GetPrimaryVertex()->GetX() + event->GetPrimaryVertex()->GetY()*event->GetPrimaryVertex()->GetY()));
    FillHistogram("hCentralityV0M",centralityPercentileV0M);
    FillHistogram("hCentralityV0A",centralityPercentileV0A);
    FillHistogram("hCentralityV0C",centralityPercentileV0C);

    Int_t trackCountAcc = 0;
    Int_t nTracks = fTrackArray->GetEntries();
    for (Int_t i = 0; i < nTracks; i++)
    {
      AliVTrack* track = static_cast<AliVTrack*>(fTrackArray->At(i));
      if (IsTrackInAcceptance(track))
      {
        FillHistogram("hTrackPhiEta", track->Phi(),track->Eta(), 1);
        FillHistogram("hTrackPt", track->Pt());
        if(track->Eta() >= 0)
          FillHistogram("hTrackPtPosEta", track->Pt());
        else
          FillHistogram("hTrackPtNegEta", track->Pt());        
                
        FillHistogram("hTrackEta", track->Eta());
        FillHistogram("hTrackPhi", track->Phi());
        FillHistogram("hTrackPhiLabel", track->Phi(), (static_cast<AliPicoTrack*>(track))->GetTrackType());
        for(Int_t j=0;j<20;j++)
          if(track->Pt() > j)
            FillHistogram("hTrackPhiPtCut", track->Phi(), track->Pt());

        FillHistogram("hTrackCharge", track->Charge());
        trackCountAcc++;
      }
    }
    FillHistogram("hTrackCountAcc", trackCountAcc, centralityPercentileV0M);

  }
  #ifdef DEBUGMODE
    AliInfo("Calculate()::QA done.");
  #endif

  ////////////////////// NOTE: Jet analysis and calculations

  if (fAnalyzeJets)
  {
    // ### SIGNAL JET ANALYSIS
    for (Int_t i = 0; i<fNumberSignalJets; i++)
    {
      AliEmcalJet* tmpJet = fSignalJets[i];

      // Jet spectra
      FillHistogram("hJetPt", tmpJet->Pt(), centralityPercentile);
      FillHistogram("hJetPtBgrdSubtractedKTImprovedCMS", GetCorrectedJetPt(tmpJet, backgroundKTImprovedCMS), centralityPercentile);
      FillHistogram("hJetPtSubtractedRhoKTImprovedCMS", tmpJet->Pt(), centralityPercentile, backgroundKTImprovedCMS);
      
      if(fAnalyzeDeprecatedBackgrounds)
      {
        FillHistogram("hJetPtBgrdSubtractedTR", GetCorrectedJetPt(tmpJet, backgroundTRCone06), centralityPercentile);
        FillHistogram("hJetPtBgrdSubtractedRC", GetCorrectedJetPt(tmpJet, backgroundRC), centralityPercentile);
        FillHistogram("hJetPtBgrdSubtractedKTPbPb", GetCorrectedJetPt(tmpJet, backgroundKTPbPb), centralityPercentile);
        FillHistogram("hJetPtBgrdSubtractedKTPbPbWithGhosts", GetCorrectedJetPt(tmpJet, backgroundKTPbPbWithGhosts), centralityPercentile);
        FillHistogram("hJetPtBgrdSubtractedKTCMS", GetCorrectedJetPt(tmpJet, backgroundKTCMS), centralityPercentile);
        FillHistogram("hJetPtBgrdSubtractedKTMean", GetCorrectedJetPt(tmpJet, backgroundKTMean), centralityPercentile);
        FillHistogram("hJetPtBgrdSubtractedKTTrackLike", GetCorrectedJetPt(tmpJet, backgroundKTTrackLike), centralityPercentile);
      }

      if(fAnalyzeQA)
      {
        FillHistogram("hJetArea", tmpJet->Area());
        FillHistogram("hJetAreaVsPt", tmpJet->Area(), tmpJet->Pt());
        FillHistogram("hJetPtVsConstituentCount", tmpJet->Pt(),tmpJet->GetNumberOfTracks());
        FillHistogram("hJetPhiEta", tmpJet->Phi(),tmpJet->Eta());
      }
      // Signal jet vs. signal jet - "Combinatorial"
      for (Int_t j = i+1; j<fNumberSignalJets; j++)
        FillHistogram("hJetDeltaPhi", GetDeltaPhi(tmpJet->Phi(), fSignalJets[j]->Phi()));

    }

    // ### DIJETS
    if(fNumberSignalJets >= 2)
    {
      FillHistogram("hLeadingJetDeltaPhi", GetDeltaPhi(fFirstLeadingJet->Phi(), fSecondLeadingJet->Phi()));

      if (IsDijet(fFirstLeadingJet, fSecondLeadingJet))
      {
        FillHistogram("hDijetConstituentsPt", fFirstLeadingJet->Pt());
        FillHistogram("hDijetConstituentsPt", fSecondLeadingJet->Pt());

        FillHistogram("hDijetLeadingJetPt", fFirstLeadingJet->Pt());
        FillHistogram("hDijetPtCorrelation", fFirstLeadingJet->Pt(), fSecondLeadingJet->Pt());
        Double_t dummyArea = 0;
        GetTRBackgroundDensity (2, backgroundDijet, dummyArea, fFirstLeadingJet, fSecondLeadingJet, kFALSE);
        GetTRBackgroundDensity (2, backgroundDijetPerpendicular, dummyArea, fFirstLeadingJet, fSecondLeadingJet, kTRUE);
      }
    }

    // ### SOME JET PLOTS
    FillHistogram("hJetCountAll", fJetArray->GetEntries());
    FillHistogram("hJetCountAccepted", fNumberSignalJets);
    if (fFirstLeadingJet)
      FillHistogram("hLeadingJetPt", fFirstLeadingJet->Pt());
    if (fSecondLeadingJet)
      FillHistogram("hSecondLeadingJetPt", fSecondLeadingJet->Pt());

  } //endif AnalyzeJets

  #ifdef DEBUGMODE
    AliInfo("Calculate()::Jets done.");
  #endif
  ////////////////////// NOTE: Background analysis

  if (fAnalyzeBackground)
  {
    // Calculate background in centrality classes
    FillHistogram("hKTBackgroundImprovedCMS", backgroundKTImprovedCMS, centralityPercentile);
    FillHistogram("hKTMeanBackgroundImprovedCMS", centralityPercentile, backgroundKTImprovedCMS);

    // In case of dijets -> look at the background
    if (backgroundDijet >= 0)
      FillHistogram("hDijetBackground", backgroundDijet, centralityPercentile); 
    if (backgroundDijetPerpendicular >= 0)
      FillHistogram("hDijetBackgroundPerpendicular", backgroundDijetPerpendicular, centralityPercentile); 
    
    if(fAnalyzeDeprecatedBackgrounds)
    {
      FillHistogram("hKTBackgroundPbPb", backgroundKTPbPb, centralityPercentile);
      FillHistogram("hKTBackgroundPbPbWithGhosts", backgroundKTPbPbWithGhosts, centralityPercentile);
      FillHistogram("hKTBackgroundCMS", backgroundKTCMS, centralityPercentile);
      FillHistogram("hKTBackgroundMean", backgroundKTMean, centralityPercentile);
      FillHistogram("hKTBackgroundTrackLike", backgroundKTTrackLike, centralityPercentile);

      FillHistogram("hTRBackgroundNoExcl", backgroundTRNoExcl, centralityPercentile);
      FillHistogram("hTRBackgroundCone02", backgroundTRCone02, centralityPercentile);
      FillHistogram("hTRBackgroundCone04", backgroundTRCone04, centralityPercentile);
      FillHistogram("hTRBackgroundCone06", backgroundTRCone06, centralityPercentile);
      FillHistogram("hTRBackgroundCone08", backgroundTRCone08, centralityPercentile);
      FillHistogram("hTRBackgroundExact", backgroundTRExact, centralityPercentile);

      FillHistogram("hRCBackground", backgroundRC, centralityPercentile);

      // Calculate background profiles in terms of centrality
      FillHistogram("hKTMeanBackgroundPbPb", centralityPercentile,  backgroundKTPbPb);
      FillHistogram("hKTMeanBackgroundPbPbWithGhosts", centralityPercentile,  backgroundKTPbPbWithGhosts);
      FillHistogram("hKTMeanBackgroundCMS", centralityPercentile, backgroundKTCMS);
      FillHistogram("hKTMeanBackgroundMean", centralityPercentile, backgroundKTMean);
      FillHistogram("hKTMeanBackgroundTPC", centralityPercentile, backgroundKTTrackLike);
      FillHistogram("hTRMeanBackground", centralityPercentile,  backgroundTRCone06);
    }


    // Calculate the delta pt

    Double_t tmpDeltaPtNoBackground = 0.0;
    Double_t tmpDeltaPtKTImprovedCMS = 0.0;

    Double_t tmpDeltaPtKTPbPb = 0.0;
    Double_t tmpDeltaPtKTPbPbWithGhosts = 0.0;
    Double_t tmpDeltaPtKTCMS = 0.0;
    Double_t tmpDeltaPtKTMean = 0.0;
    Double_t tmpDeltaPtKTTrackLike = 0.0;
    Double_t tmpDeltaPtRC = 0.0;
    Double_t tmpDeltaPtTR = 0.0;

    GetDeltaPt(tmpDeltaPtNoBackground, 0.0);
    GetDeltaPt(tmpDeltaPtKTImprovedCMS, backgroundKTImprovedCMS);

    GetDeltaPt(tmpDeltaPtKTPbPb, backgroundKTPbPb);
    GetDeltaPt(tmpDeltaPtKTPbPbWithGhosts, backgroundKTPbPbWithGhosts);
    GetDeltaPt(tmpDeltaPtKTCMS, backgroundKTCMS);
    GetDeltaPt(tmpDeltaPtKTMean, backgroundKTMean);
    GetDeltaPt(tmpDeltaPtKTTrackLike, backgroundKTTrackLike);
    GetDeltaPt(tmpDeltaPtRC, backgroundRC);
    GetDeltaPt(tmpDeltaPtTR, backgroundTRCone06);


    // If valid, fill the delta pt histograms

    if(tmpDeltaPtKTImprovedCMS > -10000.0)
      FillHistogram("hDeltaPtKTImprovedCMS", tmpDeltaPtKTImprovedCMS, centralityPercentile);
    if(tmpDeltaPtNoBackground > -10000.0)
      FillHistogram("hDeltaPtNoBackground", tmpDeltaPtNoBackground, centralityPercentile);
    if(tmpDeltaPtNoBackground > 0.000001)
      FillHistogram("hDeltaPtNoBackgroundNoEmptyCones", tmpDeltaPtNoBackground, centralityPercentile);

    if(fAnalyzeDeprecatedBackgrounds)
    {
      if(tmpDeltaPtKTPbPb > -10000.0)
        FillHistogram("hDeltaPtKTPbPb", tmpDeltaPtKTPbPb, centralityPercentile);
      if(tmpDeltaPtKTPbPbWithGhosts > -10000.0)
        FillHistogram("hDeltaPtKTPbPbWithGhosts", tmpDeltaPtKTPbPbWithGhosts, centralityPercentile);
      if(tmpDeltaPtKTCMS > -10000.0)
        FillHistogram("hDeltaPtKTCMS", tmpDeltaPtKTCMS, centralityPercentile);
      if(tmpDeltaPtKTMean > -10000.0)
        FillHistogram("hDeltaPtKTMean", tmpDeltaPtKTMean, centralityPercentile);
      if(tmpDeltaPtKTTrackLike > -10000.0)
        FillHistogram("hDeltaPtKTTrackLike", tmpDeltaPtKTTrackLike, centralityPercentile);

      if(tmpDeltaPtRC > -10000.0)
        FillHistogram("hDeltaPtRC", tmpDeltaPtRC, centralityPercentile);
      if(tmpDeltaPtTR > -10000.0)
        FillHistogram("hDeltaPtTR", tmpDeltaPtTR, centralityPercentile);
    }
  }
  
  #ifdef DEBUGMODE
    AliInfo("Calculate()::Background done.");
  #endif
  
  ////////////////////// NOTE: Pythia histograms
  if(fAnalyzePythia)
  {
    FillHistogram("hPythiaPtHard", GetPtHard());
    FillHistogram("hPythiaNTrials", GetPtHardBin()-0.1, fTrials);
    FillHistogram("hPythiaXSection", GetPtHardBin()-0.1, fCrossSection);

    #ifdef DEBUGMODE
      AliInfo("Calculate()::Pythia done.");
    #endif
  }
  #ifdef DEBUGMODE
    AliInfo("Calculate() done.");
  #endif
}

//________________________________________________________________________
Bool_t AliAnalysisTaskChargedJetsPA::Notify()
{
  // Implemented Notify() to read the cross sections
  // and number of trials from pyxsec.root
  // 
  #ifdef DEBUGMODE
    AliInfo("Notify started.");
  #endif

  if(fAnalyzePythia)
  {
    TTree *tree = AliAnalysisManager::GetAnalysisManager()->GetTree();
    TFile *currFile = tree->GetCurrentFile();

    TString file(currFile->GetName());

    if(file.Contains("root_archive.zip#")){
      Ssiz_t pos1 = file.Index("root_archive",12,TString::kExact);
      Ssiz_t pos = file.Index("#",1,pos1,TString::kExact);
      file.Replace(pos+1,20,"");
    }
    else {
      // not an archive take the basename....
      file.ReplaceAll(gSystem->BaseName(file.Data()),"");
    }
   
    TFile *fxsec = TFile::Open(Form("%s%s",file.Data(),"pyxsec.root")); // problem that we cannot really test the existance of a file in a archive so we have to lvie with open error message from root
    if(!fxsec){
      // next trial fetch the histgram file
      fxsec = TFile::Open(Form("%s%s",file.Data(),"pyxsec_hists.root"));
      if(!fxsec){
          // not a severe condition but inciate that we have no information
        return kFALSE;
      }
      else{
        // find the tlist we want to be independtent of the name so use the Tkey
        TKey* key = (TKey*)fxsec->GetListOfKeys()->At(0); 
        if(!key){
          fxsec->Close();
          return kFALSE;
        }
        TList *list = dynamic_cast<TList*>(key->ReadObj());
        if(!list){
          fxsec->Close();
          return kFALSE;
        }
        fCrossSection = ((TProfile*)list->FindObject("h1Xsec"))->GetBinContent(1);
        fTrials  = ((TH1F*)list->FindObject("h1Trials"))->GetBinContent(1);
        fxsec->Close();
      }
    } // no tree pyxsec.root
    else {
      TTree *xtree = (TTree*)fxsec->Get("Xsection");
      if(!xtree){
        fxsec->Close();
        return kFALSE;
      }
      UInt_t   ntrials  = 0;
      Double_t  xsection  = 0;
      xtree->SetBranchAddress("xsection",&xsection);
      xtree->SetBranchAddress("ntrials",&ntrials);
      xtree->GetEntry(0);
      fTrials = ntrials;
      fCrossSection = xsection;
      fxsec->Close();
    }
    #ifdef DEBUGMODE
      AliInfo("Notify ended.");
    #endif
  }
  return kTRUE;
}

//________________________________________________________________________
inline Double_t AliAnalysisTaskChargedJetsPA::EtaToTheta(Double_t arg)
  {return 2.*atan(exp(-arg));} 
//________________________________________________________________________
inline Double_t AliAnalysisTaskChargedJetsPA::ThetaToEta(Double_t arg)
{
  if ((arg > TMath::Pi()) || (arg < 0.0))
  {
    AliError(Form("ThetaToEta got wrong input! (%f)", arg));
    return 0.0;
  }
  return -log(tan(arg/2.));
}
//________________________________________________________________________
inline Double_t AliAnalysisTaskChargedJetsPA::GetDeltaPhi(Double_t phi1, Double_t phi2)
  {return min(TMath::Abs(phi1-phi2),TMath::TwoPi() - TMath::Abs(phi1-phi2));}

//________________________________________________________________________
Double_t AliAnalysisTaskChargedJetsPA::MCGetOverlapCircleRectancle(Double_t cPosX, Double_t cPosY, Double_t cRadius, Double_t rPosXmin, Double_t rPosXmax, Double_t rPosYmin, Double_t rPosYmax)
{
  const Int_t kTests = 1000;
  Int_t hits = 0;
  TRandom3 randomGen(0);
 
  // Loop over kTests-many tests
  for (Int_t i=0; i<kTests; i++)
  {
    //Choose random position in rectangle for the tester
    Double_t tmpTestX = randomGen.Uniform(rPosXmin, rPosXmax);
    Double_t tmpTestY = randomGen.Uniform(rPosYmin, rPosYmax);

    //Check, if tester is in circle. If yes, increment circle counter.
    Double_t tmpDistance = TMath::Sqrt( (tmpTestX - cPosX)*(tmpTestX - cPosX) + (tmpTestY - cPosY)*(tmpTestY - cPosY) );
    if(tmpDistance < cRadius)
      hits++;
  }

  // return ratio
  return (static_cast<Double_t>(hits)/static_cast<Double_t>(kTests));
}

//________________________________________________________________________
Double_t AliAnalysisTaskChargedJetsPA::MCGetOverlapMultipleCirclesRectancle(Int_t numCircles, std::vector<Double_t> cPosX, std::vector<Double_t> cPosY, Double_t cRadius, Double_t rPosXmin, Double_t rPosXmax, Double_t rPosYmin, Double_t rPosYmax)
{

  const Int_t kTests = 1000;
  Int_t hits = 0;
  TRandom3 randomGen(0);
 
  // Loop over kTests-many tests
  for (Int_t i=0; i<kTests; i++)
  {
    //Choose random position in rectangle for the tester
    Double_t tmpTestX = randomGen.Uniform(rPosXmin, rPosXmax);
    Double_t tmpTestY = randomGen.Uniform(rPosYmin, rPosYmax);

    //Check, if tester is in one of the circles. If yes, increment circle counter.
    for(Int_t j=0; j<numCircles; j++)
    {
      Double_t tmpDistance = TMath::Sqrt( (tmpTestX - cPosX[j])*(tmpTestX - cPosX[j]) + (tmpTestY - cPosY[j])*(tmpTestY - cPosY[j]) );
      if(tmpDistance < cRadius)
      {
        hits++;
        break;
      }
    }
  }

  // return ratio
  return (static_cast<Double_t>(hits)/static_cast<Double_t>(kTests));

}

//________________________________________________________________________
inline void AliAnalysisTaskChargedJetsPA::FillHistogram(const char * key, Double_t x)
{
  TH1* tmpHist = static_cast<TH1*>(fOutputList->FindObject(GetHistoName(key)));
  if(!tmpHist)
  {
    AliWarning(Form("Cannot find histogram <%s> ",key)) ;
    return;
  }

  tmpHist->Fill(x);
}

//________________________________________________________________________
inline void AliAnalysisTaskChargedJetsPA::FillHistogram(const char * key, Double_t x, Double_t y)
{
  TH1* tmpHist = static_cast<TH1*>(fOutputList->FindObject(GetHistoName(key)));
  if(!tmpHist)
  {
    AliWarning(Form("Cannot find histogram <%s> ",key));
    return;
  }

  if (tmpHist->IsA()->GetBaseClass("TH1"))
    static_cast<TH1*>(tmpHist)->Fill(x,y); // Fill x with y
  else if (tmpHist->IsA()->GetBaseClass("TH2"))
    static_cast<TH2*>(tmpHist)->Fill(x,y); // Fill x,y with 1
}

//________________________________________________________________________
inline void AliAnalysisTaskChargedJetsPA::FillHistogram(const char * key, Double_t x, Double_t y, Double_t add)
{
  TH2* tmpHist = static_cast<TH2*>(fOutputList->FindObject(GetHistoName(key)));
  if(!tmpHist)
  {
    AliWarning(Form("Cannot find histogram <%s> ",key));
    return;
  }
  
  tmpHist->Fill(x,y,add);
}
//________________________________________________________________________
template <class T> T* AliAnalysisTaskChargedJetsPA::AddHistogram1D(const char* name, const char* title, const char* options, Int_t xBins, Double_t xMin, Double_t xMax, const char* xTitle, const char* yTitle)
{
  T* tmpHist = new T(GetHistoName(name), GetHistoName(title), xBins, xMin, xMax);

  tmpHist->GetXaxis()->SetTitle(xTitle);
  tmpHist->GetYaxis()->SetTitle(yTitle);
  tmpHist->SetOption(options);
  tmpHist->SetMarkerStyle(kFullCircle);
  tmpHist->Sumw2();

  fHistList->Add(tmpHist);
  fHistCount++;
  
  return tmpHist;
}

//________________________________________________________________________
template <class T> T* AliAnalysisTaskChargedJetsPA::AddHistogram2D(const char* name, const char* title, const char* options, Int_t xBins, Double_t xMin, Double_t xMax, Int_t yBins, Double_t yMin, Double_t yMax, const char* xTitle, const char* yTitle, const char* zTitle)
{
  T* tmpHist = new T(GetHistoName(name), GetHistoName(title), xBins, xMin, xMax, yBins, yMin, yMax);
  tmpHist->GetXaxis()->SetTitle(xTitle);
  tmpHist->GetYaxis()->SetTitle(yTitle);
  tmpHist->GetZaxis()->SetTitle(zTitle);
  tmpHist->SetOption(options);
  tmpHist->SetMarkerStyle(kFullCircle);
  tmpHist->Sumw2();

  fHistList->Add(tmpHist);
  fHistCount++;

  return tmpHist;
}

//________________________________________________________________________
void AliAnalysisTaskChargedJetsPA::Terminate(Option_t *)
{
  PostData(1, fOutputList);

  // Mandatory
  fOutputList = dynamic_cast<TList*> (GetOutputData(1)); // '1' refers to the output slot
  if (!fOutputList) {
    printf("ERROR: Output list not available\n");
    return;
  }
}

//________________________________________________________________________
AliAnalysisTaskChargedJetsPA::~AliAnalysisTaskChargedJetsPA()
{
  // Destructor. Clean-up the output list, but not the histograms that are put inside
  // (the list is owner and will clean-up these histograms). Protect in PROOF case.
  if (fOutputList && !AliAnalysisManager::GetAnalysisManager()->IsProofMode()) {
    delete fOutputList;
  }
}

//________________________________________________________________________
void AliAnalysisTaskChargedJetsPA::UserCreateOutputObjects()
{
  // called once to create user defined output objects like histograms, plots etc. 
  // and to put it on the output list.
  // Note: Saving to file with e.g. OpenFile(0) is must be before creating other objects.

  fRandom = new TRandom3(0);
  
  fOutputList = new TList();
  fOutputList->SetOwner(); // otherwise it produces leaks in merging

  PostData(1, fOutputList);
}

//________________________________________________________________________
void AliAnalysisTaskChargedJetsPA::UserExec(Option_t *) 
{
  #ifdef DEBUGMODE
    AliInfo("UserExec() started.");
  #endif

  if (!InputEvent())
  {
    AliError("??? Event pointer == 0 ???");
    return;
  }

  if (!fInitialized)
    ExecOnce(); // Get tracks, jets, background from arrays if not already given + Init Histos
  
  Calculate(InputEvent());
        
  PostData(1, fOutputList);
}
