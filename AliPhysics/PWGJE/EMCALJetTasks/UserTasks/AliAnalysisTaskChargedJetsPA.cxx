#ifndef ALIANALYSISTASKSE_H
#include <Riostream.h>
#include <TROOT.h>
#include <TFile.h>
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

#include "AliAnalysisDataContainer.h"
#include <THn.h>
#include "TFormula.h"
#include "AliESDtrackCuts.h"
#include <time.h>
#include <TRandom3.h>
#include "AliGenPythiaEventHeader.h"
#include "AliAODMCHeader.h"
#include "AliMCEvent.h"
#include "AliLog.h"
#include <AliEmcalJet.h>
#include <AliPicoTrack.h>
#include "AliVEventHandler.h"
#include "AliVParticle.h"
#include "AliAODMCParticle.h"
#include "AliAnalysisUtils.h"
#include "AliRhoParameter.h"
#include "TVector3.h"

#include "AliAnalysisTaskChargedJetsPA.h"
using std::min;
using std::cout;
using std::endl;

const double AliAnalysisTaskChargedJetsPA::kMaxChi2 = 200;

//TODO: FillHistogram can be done better with virtual TH1(?)
ClassImp(AliAnalysisTaskChargedJetsPA)

// ######################################################################################## DEFINE HISTOGRAMS
void AliAnalysisTaskChargedJetsPA::Init()
{
  #ifdef DEBUGMODE
    AliInfo("Creating histograms.");
  #endif

  SetCurrentOutputList(0);

  TH1* tmpHisto = AddHistogram1D<TH1D>("hVertexAcceptance", "Accepted vertices for different conditions", "", 4, 0, 4, "stage","N^{Events}/cut");
  tmpHisto->GetXaxis()->SetBinLabel(1, "Triggered all");
  tmpHisto->GetXaxis()->SetBinLabel(2, "Triggered w/ vertex");
  tmpHisto->GetXaxis()->SetBinLabel(3, "Pile-up corrected all");
  tmpHisto->GetXaxis()->SetBinLabel(4, "Pile-up corrected w vertex");

  TH2* tmpHisto2D = AddHistogram2D<TH2D>("hCentrality", Form("Accepted events in centrality (%s)", fCentralityType.Data()), "COLZ", 102, 0., 102., 4, 0,4,"Centrality","Cut stage","dN^{Events}");
  tmpHisto2D->GetYaxis()->SetBinLabel(1, "Before cuts");
  tmpHisto2D->GetYaxis()->SetBinLabel(2, "After pile up");
  tmpHisto2D->GetYaxis()->SetBinLabel(3, "After vertex demand");
  tmpHisto2D->GetYaxis()->SetBinLabel(4, "After vertex cuts");

  tmpHisto = AddHistogram1D<TH1D>("hTrackAcceptance", "Accepted tracks (0 = before cuts, 1 = after eta, 2 = after pT)", "", 3, 0, 3, "stage","N^{Tracks}/cut");
  tmpHisto->GetXaxis()->SetBinLabel(1, "Before cuts");
  tmpHisto->GetXaxis()->SetBinLabel(2, "After eta");
  tmpHisto->GetXaxis()->SetBinLabel(3, "After p_{T}");

  tmpHisto = AddHistogram1D<TH1D>("hJetAcceptance", "Accepted jets (0 = before cuts, 1 = after eta, 2 = after pT, 3 = after area)", "", 4, 0, 4, "stage","N^{Jets}/cut");
  tmpHisto->GetXaxis()->SetBinLabel(1, "Before cuts");
  tmpHisto->GetXaxis()->SetBinLabel(2, "After eta");
  tmpHisto->GetXaxis()->SetBinLabel(3, "After p_{T}");
  tmpHisto->GetXaxis()->SetBinLabel(4, "After area");

  tmpHisto2D = AddHistogram2D<TH2D>("hJetPtCutStages", "Jets p_{T} distribution", "", 500, -50., 200., 4, 0, 4, "p_{T} (GeV/c)","Cut stage","dN^{Jets}/dp_{T}");
  tmpHisto2D->GetYaxis()->SetBinLabel(1, "Before cuts");
  tmpHisto2D->GetYaxis()->SetBinLabel(2, "After eta");
  tmpHisto2D->GetYaxis()->SetBinLabel(3, "After p_{T}");
  tmpHisto2D->GetYaxis()->SetBinLabel(4, "After area");

  AddHistogram1D<TH1D>("hVertexX", "X distribution of the vertex", "", 2000, -1., 1., "#Delta x(cm)","dN^{Events}/dx");
  AddHistogram1D<TH1D>("hVertexY", "Y distribution of the vertex", "", 2000, -1., 1., "#Delta y(cm)","dN^{Events}/dy");
  AddHistogram2D<TH2D>("hVertexXY", "XY distribution of the vertex", "COLZ", 500, -1., 1., 500, -1., 1.,"#Delta x(cm)", "#Delta y(cm)","dN^{Events}/dxdy");
  AddHistogram1D<TH1D>("hVertexZ", "Z distribution of the vertex (after std. vertex cut)", "", 200, -20., 20., "#Delta z(cm)","dN^{Events}/dz");
  AddHistogram1D<TH1D>("hVertexR", "R distribution of the vertex", "", 100, 0., 1., "#Delta r(cm)","dN^{Events}/dr");
  AddHistogram1D<TH1D>("hCentralityV0M", "Centrality distribution V0M", "", fNumberOfCentralityBins, 0., 100., "Centrality","dN^{Events}");
  AddHistogram1D<TH1D>("hCentralityCL1", "Centrality distribution CL1", "", fNumberOfCentralityBins, 0., 100., "Centrality","dN^{Events}");
  AddHistogram1D<TH1D>("hCentralityV0A", "Centrality distribution V0A", "", fNumberOfCentralityBins, 0., 100., "Centrality","dN^{Events}");
  AddHistogram1D<TH1D>("hCentralityV0C", "Centrality distribution V0C", "", fNumberOfCentralityBins, 0., 100., "Centrality","dN^{Events}");
  AddHistogram1D<TH1D>("hCentralityZNA", "Centrality distribution ZNA", "", fNumberOfCentralityBins, 0., 100., "Centrality","dN^{Events}");

  if(fDoJetAnalysis)
  {
    // Background corrected jet spectra
    AddHistogram2D<TH2D>("hJetPtNoBgrdSubtracted", "Jets p_{T} distribution, no bgrd. subtracted", "", 500, -50., 200., fNumberOfCentralityBins, 0, 100, "p_{T} (GeV/c)","Centrality","dN^{Jets}/dp_{T}");    
    AddHistogram2D<TH2D>("hJetPtBgrdSubtractedExternal", "Jets p_{T} distribution, external bgrd. subtracted", "", 500, -50., 200., fNumberOfCentralityBins, 0, 100, "p_{T} (GeV/c)","Centrality","dN^{Jets}/dp_{T}");    
    AddHistogram2D<TH2D>("hJetPtBgrdSubtractedPP", "Jets p_{T} distribution, pp background subtracted", "", 500, -50., 200., fNumberOfCentralityBins, 0, 100, "p_{T} (GeV/c)","Centrality","dN^{Jets}/dp_{T}");
    AddHistogram2D<TH2D>("hJetPtBgrdSubtractedExternal_Phi1", "Jets p_{T} distribution, external background (Improved CMS) subtracted (1st part of azimuth)", "", 500, -50., 200., fNumberOfCentralityBins, 0, 100, "p_{T} (GeV/c)","Centrality","dN^{Jets}/dp_{T}");    
    AddHistogram2D<TH2D>("hJetPtBgrdSubtractedExternal_Phi2", "Jets p_{T} distribution, external background (Improved CMS) subtracted (2nd part of azimuth)", "", 500, -50., 200., fNumberOfCentralityBins, 0, 100, "p_{T} (GeV/c)","Centrality","dN^{Jets}/dp_{T}");    
    AddHistogram2D<TH2D>("hJetPtBgrdSubtractedKTImprovedCMS", "Jets p_{T} distribution, KT background (Improved CMS) subtracted", "", 500, -50., 200., fNumberOfCentralityBins, 0, 100, "p_{T} (GeV/c)","Centrality","dN^{Jets}/dp_{T}");
    AddHistogram2D<TH2D>("hJetPtBgrdSubtractedKTImprovedCMS_Biased_10GeV", "Jets p_{T} distribution, KT background (Improved CMS) subtracted, leading track bias 10 GeV", "", 500, -50., 200., fNumberOfCentralityBins, 0, 100, "p_{T} (GeV/c)","Centrality","dN^{Jets}/dp_{T}");
    AddHistogram2D<TH2D>("hJetPtBgrdSubtractedKTImprovedCMS_Biased_5GeV", "Jets p_{T} distribution, KT background (Improved CMS) subtracted, leading track bias 5 GeV", "", 500, -50., 200., fNumberOfCentralityBins, 0, 100, "p_{T} (GeV/c)","Centrality","dN^{Jets}/dp_{T}");
    AddHistogram2D<TH2D>("hJetPtBgrdSubtractedKTImprovedCMS_Biased_2GeV", "Jets p_{T} distribution, KT background (Improved CMS) subtracted, leading track bias 2 GeV", "", 500, -50., 200., fNumberOfCentralityBins, 0, 100, "p_{T} (GeV/c)","Centrality","dN^{Jets}/dp_{T}");
    AddHistogram2D<TH2D>("hJetPtBgrdSubtractedTR", "Jets p_{T} distribution, TR background (Cone R=0.6 around jets excluded) subtracted", "", 500, -50., 200., fNumberOfCentralityBins, 0, 100, "p_{T} (GeV/c)","Centrality","dN^{Jets}/dp_{T}");
    AddHistogram2D<TH2D>("hJetPtBgrdSubtractedKTPbPb", "Jets p_{T} distribution, KT background (PbPb w/o ghosts) subtracted", "", 500, -50., 200., fNumberOfCentralityBins, 0, 100, "p_{T} (GeV/c)","Centrality","dN^{Jets}/dp_{T}");
    AddHistogram2D<TH2D>("hJetPtBgrdSubtractedKTPbPbWithGhosts", "Jets p_{T} distribution, KT background (PbPb w/ ghosts) subtracted", "", 500, -50., 200., fNumberOfCentralityBins, 0, 100, "p_{T} (GeV/c)","Centrality","dN^{Jets}/dp_{T}");
    AddHistogram2D<TH2D>("hJetPtBgrdSubtractedKTCMS", "Jets p_{T} distribution, KT background (CMS) subtracted", "", 500, -50., 200., fNumberOfCentralityBins, 0, 100, "p_{T} (GeV/c)","Centrality","dN^{Jets}/dp_{T}");    
    AddHistogram2D<TH2D>("hJetPtBgrdSubtractedKTMean", "Jets p_{T} distribution, KT background (Mean) subtracted", "", 500, -50., 200., fNumberOfCentralityBins, 0, 100, "p_{T} (GeV/c)","Centrality","dN^{Jets}/dp_{T}");    
    AddHistogram2D<TH2D>("hJetPtBgrdSubtractedKTTrackLike", "Jets p_{T} distribution, KT background (track-like) subtracted", "", 500, -50., 200., fNumberOfCentralityBins, 0, 100, "p_{T} (GeV/c)","Centrality","dN^{Jets}/dp_{T}");

    AddHistogram2D<TProfile2D>("hJetPtSubtractedRhoExternal", "Mean subtracted KT (External) background from jets", "COLZ", 600, 0, 150, fNumberOfCentralityBins, 0, 100, "Jet p_{T}", "Centrality", "#rho mean");
    AddHistogram2D<TProfile2D>("hJetPtSubtractedRhoKTImprovedCMS", "Mean subtracted KT (CMS w/o signal) background from jets", "COLZ", 600, 0, 150, fNumberOfCentralityBins, 0, 100, "Jet p_{T}", "Centrality", "#rho mean");
    AddHistogram2D<TProfile2D>("hJetPtSubtractedRhoPP", "Mean subtracted KT (pp from Michal) background from jets", "COLZ", 600, 0, 150, fNumberOfCentralityBins, 0, 100, "Jet p_{T}", "Centrality", "#rho mean");

    // Jet QA plots
    AddHistogram2D<TH2D>("hJetConstituentPt0GeV", "Jet constituents p_{T} distribution (p_{T,jet} > 0 GeV)", "", 500, -50., 200., fNumberOfCentralityBins, 0, 100, "p_{T} (GeV/c)","Centrality","dN^{Tracks}/dp_{T}");
    AddHistogram2D<TH2D>("hJetConstituentPt1GeV", "Jet constituents p_{T} distribution (p_{T,jet} > 1 GeV)", "", 500, -50., 200., fNumberOfCentralityBins, 0, 100, "p_{T} (GeV/c)","Centrality","dN^{Tracks}/dp_{T}");
    AddHistogram2D<TH2D>("hJetConstituentPt2GeV", "Jet constituents p_{T} distribution (p_{T,jet} > 2 GeV)", "", 500, -50., 200., fNumberOfCentralityBins, 0, 100, "p_{T} (GeV/c)","Centrality","dN^{Tracks}/dp_{T}");
    AddHistogram2D<TH2D>("hJetConstituentPt3GeV", "Jet constituents p_{T} distribution (p_{T,jet} > 3 GeV)", "", 500, -50., 200., fNumberOfCentralityBins, 0, 100, "p_{T} (GeV/c)","Centrality","dN^{Tracks}/dp_{T}");
    AddHistogram2D<TH2D>("hJetConstituentPt4GeV", "Jet constituents p_{T} distribution (p_{T,jet} > 4 GeV)", "", 500, -50., 200., fNumberOfCentralityBins, 0, 100, "p_{T} (GeV/c)","Centrality","dN^{Tracks}/dp_{T}");
    AddHistogram2D<TH2D>("hJetConstituentPt5GeV", "Jet constituents p_{T} distribution (p_{T,jet} > 5 GeV)", "", 500, -50., 200., fNumberOfCentralityBins, 0, 100, "p_{T} (GeV/c)","Centrality","dN^{Tracks}/dp_{T}");
    AddHistogram2D<TH2D>("hJetConstituentPt7GeV", "Jet constituents p_{T} distribution (p_{T,jet} > 7 GeV)", "", 500, -50., 200., fNumberOfCentralityBins, 0, 100, "p_{T} (GeV/c)","Centrality","dN^{Tracks}/dp_{T}");
    AddHistogram2D<TH2D>("hJetConstituentPt10GeV", "Jet constituents p_{T} distribution (p_{T,jet} > 10 GeV)", "", 500, -50., 200., fNumberOfCentralityBins, 0, 100, "p_{T} (GeV/c)","Centrality","dN^{Tracks}/dp_{T}");
    AddHistogram2D<TH2D>("hJetConstituentPtVsJetPt", "Jet constituents p_{T} distribution", "", 500, -50., 200., 200, 0, 200, "#it{p}_{T} (GeV/c)","#it{p}_{T}^{jet} (GeV/c)","dN^{Tracks}/dp_{T}");
    AddHistogram1D<TH1D>("hJetCountAll", "Number of Jets", "", 200, 0., 200., "N jets","dN^{Events}/dN^{Jets}");
    AddHistogram1D<TH1D>("hJetCountAccepted", "Number of accepted Jets", "", 200, 0., 200., "N jets","dN^{Events}/dN^{Jets}");
    AddHistogram2D<TH2D>("hJetCount", "Correlation jets/accepted jets", "", 200, 0., 200., 200, 0., 200., "N jets","N jets accepted", "d^{2}N^{Events}/dN^{Jets dN^{Jets, acc}}");
    AddHistogram1D<TH1D>("hLeadingJetPt", "Leading jet p_{T}", "", 500, -50., 200., "p_{T} (GeV/c)","dN^{Jets}/dp_{T}");
    AddHistogram1D<TH1D>("hSecondLeadingJetPt", "Second leading jet p_{T}", "", 500, -50., 200., "p_{T} (GeV/c)","dN^{Jets}/dp_{T}");
    AddHistogram1D<TH1D>("hCorrectedLeadingJetPt", "Corrected leading jet p_{T}", "", 500, -50., 200., "p_{T} (GeV/c)","dN^{Jets}/dp_{T}");
    AddHistogram1D<TH1D>("hCorrectedSecondLeadingJetPt", "Corrected second leading jet p_{T}", "", 500, -50., 200., "p_{T} (GeV/c)","dN^{Jets}/dp_{T}");
    AddHistogram1D<TH1D>("hJetDeltaPhi", "Jets combinatorial #Delta #phi", "", 250, 0., TMath::Pi(), "#Delta #phi","dN^{Jets}/d(#Delta #phi)");
    AddHistogram1D<TH1D>("hLeadingJetDeltaPhi", "1st and 2nd leading jet #Delta #phi", "", 250, 0., TMath::Pi(), "#Delta #phi","dN^{Jets}/d(#Delta #phi)");

    // Background distributions

    AddHistogram2D<TH2D>("hKTBackgroundExternal", "KT background density (External task)", "LEGO2", 400, 0., 40., fNumberOfCentralityBins, 0, 100, "#rho (GeV/c)","Centrality", "dN^{Events}/d#rho");
    AddHistogram2D<TH2D>("hKTBackgroundExternalVsPt", "KT background density (External task)", "LEGO2", 400, 0., 40., 200, 0, 200, "#rho (GeV/c)","Raw jet p_{T}", "dN^{Events}/d#rho");
    AddHistogram2D<TH2D>("hKTBackgroundExternal20GeV", "KT background density (External task, jet p_{T} > 20 GeV)", "LEGO2", 400, 0., 40., fNumberOfCentralityBins, 0, 100, "#rho (GeV/c)","Centrality", "dN^{Events}/d#rho");
    AddHistogram2D<TH2D>("hKTBackgroundImprovedCMS", "KT background density (Improved CMS approach)", "LEGO2", 400, 0., 40., fNumberOfCentralityBins, 0, 100, "#rho (GeV/c)","Centrality", "dN^{Events}/d#rho");
    AddHistogram2D<TH2D>("hPPBackground", "PP background density (Michals approach)", "LEGO2", 400, 0., 40., fNumberOfCentralityBins, 0, 100, "#rho (GeV/c)","Centrality", "dN^{Events}/d#rho");
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

    // Delta pt distributions
    AddHistogram2D<TH2D>("hDeltaPtPP", "Background fluctuations #delta p_{T} (PP approach)", "", 1801, -40.0, 80.0, fNumberOfCentralityBins, 0, 100, "#delta p_{T} (GeV/c)","Centrality","dN^{Jets}/d#delta p_{T}");
    AddHistogram2D<TH2D>("hDeltaPtExternalBgrd", "Background fluctuations #delta p_{T} (KT, External)", "", 1801, -40.0, 80.0, fNumberOfCentralityBins, 0, 100, "#delta p_{T} (GeV/c)","Centrality","dN^{Jets}/d#delta p_{T}");
    AddHistogram2D<TH2D>("hDeltaPtExternalBgrdVsPt", "Background fluctuations #delta p_{T} (KT, External, in p_{T} bins)", "", 1801, -40.0, 80.0, 200, 0, 200, "#delta p_{T} (GeV/c)","Raw jet p_{T}","dN^{Jets}/d#delta p_{T}");
    AddHistogram2D<TH2D>("hDeltaPtExternalBgrdPartialExclusion", "Background fluctuations #delta p_{T} (KT, External, partial jet exclusion)", "", 1801, -40.0, 80.0, fNumberOfCentralityBins, 0, 100, "#delta p_{T} (GeV/c)","Centrality","dN^{Jets}/d#delta p_{T}");
    AddHistogram2D<TH2D>("hDeltaPtKTImprovedCMS", "Background fluctuations #delta p_{T} (KT, Improved CMS-like)", "", 1801, -40.0, 80.0, fNumberOfCentralityBins, 0, 100, "#delta p_{T} (GeV/c)","Centrality","dN^{Jets}/d#delta p_{T}");
    AddHistogram2D<TH2D>("hDeltaPtNoBackground", "Background fluctuations #delta p_{T} (No background)", "", 1801, -40.0, 80.0, fNumberOfCentralityBins, 0, 100, "#delta p_{T} (GeV/c)","Centrality","dN^{Jets}/d#delta p_{T}");
    AddHistogram2D<TH2D>("hDeltaPtKTPbPb", "Background fluctuations #delta p_{T} (KT, PbPb w/o ghosts)", "", 1801, -40.0, 80.0, fNumberOfCentralityBins, 0, 100, "#delta p_{T} (GeV/c)","Centrality","dN^{Jets}/d#delta p_{T}");
    AddHistogram2D<TH2D>("hDeltaPtKTPbPbWithGhosts", "Background fluctuations #delta p_{T} (KT, PbPb w/ ghosts)", "", 1801, -40.0, 80.0, fNumberOfCentralityBins, 0, 100, "#delta p_{T} (GeV/c)","Centrality","dN^{Jets}/d#delta p_{T}");
    AddHistogram2D<TH2D>("hDeltaPtKTCMS", "Background fluctuations #delta p_{T} (KT, CMS-like)", "", 1801, -40.0, 80.0, fNumberOfCentralityBins, 0, 100, "#delta p_{T} (GeV/c)","Centrality","dN^{Jets}/d#delta p_{T}");
    AddHistogram2D<TH2D>("hDeltaPtKTMean", "Background fluctuations #delta p_{T} (KT, Mean)", "", 1801, -40.0, 80.0, fNumberOfCentralityBins, 0, 100, "#delta p_{T} (GeV/c)","Centrality","dN^{Jets}/d#delta p_{T}");
    AddHistogram2D<TH2D>("hDeltaPtKTTrackLike", "Background fluctuations #delta p_{T} (KT, track-like)", "", 1801, -40.0, 80.0, fNumberOfCentralityBins, 0, 100, "#delta p_{T} (GeV/c)","Centrality","dN^{Jets}/d#delta p_{T}");
    AddHistogram2D<TH2D>("hDeltaPtTR", "Background fluctuations #delta p_{T} (TR, cone R=0.6)", "", 1801, -40.0, 80.0, fNumberOfCentralityBins, 0, 100,  "#delta p_{T} (GeV/c)","Centrality","dN^{Jets}/d#delta p_{T}");

    // Track QA plots
    AddHistogram2D<TH2D>("hTrackCountAcc", "Number of tracks in acceptance vs. centrality", "LEGO2", 750, 0., 750., fNumberOfCentralityBins, 0, 100, "N tracks","Centrality", "dN^{Events}/dN^{Tracks}");
    AddHistogram2D<TH2D>("hTrackPt", "Tracks p_{T} distribution", "", 1000, 0., 250., fNumberOfCentralityBins, 0, 100, "p_{T} (GeV/c)", "Centrality", "dN^{Tracks}/dp_{T}");
    AddHistogram2D<TH2D>("hTrackPtNegEta", "Tracks p_{T} distribution (negative #eta)", "", 1000, 0., 250., fNumberOfCentralityBins, 0, 100, "p_{T} (GeV/c)","Centrality","dN^{Tracks}/dp_{T}");
    AddHistogram2D<TH2D>("hTrackPtPosEta", "Tracks p_{T} distribution (positive #eta)", "", 1000, 0., 250., fNumberOfCentralityBins, 0, 100, "p_{T} (GeV/c)","Centrality","dN^{Tracks}/dp_{T}");
    AddHistogram1D<TH1D>("hTrackCharge", "Charge", "", 11, -5, 5, "Charge (e)","dN^{Tracks}/dq");
    AddHistogram1D<TH1D>("hTrackPhi", "Track #phi distribution", "", 360, 0, TMath::TwoPi(), "#phi","dN^{Tracks}/d#phi");
    AddHistogram2D<TH2D>("hTrackPhiEta", "Track angular distribution", "LEGO2", 100, 0., 2*TMath::Pi(),100, -2.5, 2.5, "#phi","#eta","dN^{Tracks}/(d#phi d#eta)");
    AddHistogram2D<TH2D>("hTrackPtPhiEta", "Track p_{T} angular distribution", "LEGO2", 100, 0., 2*TMath::Pi(),100, -2.5, 2.5, "#phi","#eta","dp_{T}^{Tracks}/(d#phi d#eta)");
    AddHistogram2D<TH2D>("hTrackPhiPtCut", "Track #phi distribution for different pT cuts", "LEGO2", 360, 0, TMath::TwoPi(), 20, 0, 20, "#phi", "p_{T} lower cut", "dN^{Tracks}/d#phi dp_{T}");
    AddHistogram2D<TH2D>("hTrackPhiTrackType", "Track #phi distribution for different track types", "LEGO2", 360, 0, TMath::TwoPi(), 3, 0, 3, "#phi", "Label", "dN^{Tracks}/d#phi");
    AddHistogram2D<TH2D>("hTrackPtTrackType", "Track p_{T} distribution for different track types", "LEGO2", 1000, 0., 250., 3, 0, 3, "p_{T} (GeV/c)", "Label", "dN^{Tracks}/dp_{T}");
    AddHistogram2D<TH2D>("hTrackEta", "Track #eta distribution", "COLZ", 180, fMinEta, fMaxEta, fNumberOfCentralityBins, 0., 100., "#eta", "Centrality", "dN^{Tracks}/d#eta");

    // Jet QA plots
    AddHistogram1D<TH1D>("hRawJetArea", "Jets area distribution w/o area cut", "", 200, 0., 2., "Area","dN^{Jets}/dA");
    AddHistogram2D<TH2D>("hJetArea", "Jets area distribution", "COLZ", 200, 0., 2.,  500, -50., 200, "Area","Jet p_{T}","dN^{Jets}/dA");
    AddHistogram2D<TH2D>("hRawJetPhiEta", "Raw Jets angular distribution w/o #eta cut", "LEGO2", 360, 0., 2*TMath::Pi(),100, -1.0, 1.0, "#phi","#eta","dN^{Jets}/(d#phi d#eta)");
    AddHistogram2D<TH2D>("hJetEta", "Jets #eta distribution", "COLZ", 180, fMinEta, fMaxEta, fNumberOfCentralityBins, 0., 100., "#eta", "Centrality", "dN^{Jets}/d#eta");
    AddHistogram2D<TH2D>("hJetEta2GeVTracks", "Jets #eta distribution, track p_{T} > 2 GeV", "COLZ", 180, fMinEta, fMaxEta, fNumberOfCentralityBins, 0., 100., "#eta", "Centrality", "dN^{Jets}/d#eta");
    AddHistogram2D<TH2D>("hJetEta4GeVTracks", "Jets #eta distribution, track p_{T} > 4 GeV", "COLZ", 180, fMinEta, fMaxEta, fNumberOfCentralityBins, 0., 100., "#eta", "Centrality", "dN^{Jets}/d#eta");
    AddHistogram2D<TH2D>("hJetPhiEta", "Jets angular distribution", "LEGO2", 360, 0., 2*TMath::Pi(),100, -1.0, 1.0, "#phi","#eta","dN^{Jets}/(d#phi d#eta)");
    AddHistogram2D<TH2D>("hJetPtPhiEta", "Jets p_{T} angular distribution", "LEGO2", 360, 0., 2*TMath::Pi(),100, -1.0, 1.0, "#phi","#eta","dp_{T}^{Jets}/(d#phi d#eta)");
    AddHistogram2D<TH2D>("hJetPtVsConstituentCount", "Jets number of constituents vs. jet p_{T}", "COLZ", 400, 0., 200., 100, 0., 100., "p_{T}","N^{Tracks}","dN^{Jets}/(dp_{T} dN^{tracks})");

    // ######## Jet constituent analysis

    if(fAnalyzeJetConstituents)
    {
      {
        //                        jet pt,  const pT,  const count,  RC const count,  PC const count
        Int_t    bins [5]     = { 30,         50,           30,              30,              30};
        Double_t minEdges[5]  = { 0,              0.1,            0,               0,             0};
        Double_t maxEdges[5]  = { 150,          150,           30,              30,              30};
        TString axisName[5]  = {"jet p_{T}","Constituent p_{T}", "Constituent count","RC constituent count","PC constituent count"};
        TString axisTitle[5]  = {"jet p_{T}","Constituent p_{T}", "Constituent count","RC constituent count","PC constituent count"};
        THnF * histJetConstituents = new THnF("hJetConstituents", "Jet constituent count/p_{T} in jet, RC, and PC", 5, bins, minEdges, maxEdges);
        BinLogAxis(histJetConstituents,1);
        for (Int_t iaxis=0; iaxis<5;iaxis++){
          histJetConstituents->GetAxis(iaxis)->SetName(axisName[iaxis]);
          histJetConstituents->GetAxis(iaxis)->SetTitle(axisTitle[iaxis]);
        }
        fCurrentOutputList->Add(histJetConstituents);
      }

      {
        //                        jet pt,  const pt,   const count      distance
        Int_t    bins [4]     = { 30,         50,        30,     50};
        Double_t minEdges[4]  = { 0,           0.1,       0,       0};
        Double_t maxEdges[4]  = { 150,          150,     30,      0.5};
        TString axisName[4]  = {"jet p_{T}","Constituent p_{T}","Constituent count","Distance from jet axis"};
        TString axisTitle[4]  = {"jet p_{T}","Constituent p_{T}","Constituent count","Distance from jet axis"};
        THnF * histJetConstituentDistance = new THnF("hJetConstituentDistance", "Jet constituent distance vs. jet and constituent p_{T}", 4, bins, minEdges, maxEdges);
        BinLogAxis(histJetConstituentDistance,1);
        for (Int_t iaxis=0; iaxis<4;iaxis++){
          histJetConstituentDistance->GetAxis(iaxis)->SetName(axisName[iaxis]);
          histJetConstituentDistance->GetAxis(iaxis)->SetTitle(axisTitle[iaxis]);
        }
        fCurrentOutputList->Add(histJetConstituentDistance);
      }
    }

    // ######## Jet profiles
    if(fAnalyzeJetProfile)
    {
      SetCurrentOutputList(1);
      AddHistogram2D<TH2D>("hJetProfile10GeV", "Jet profile, cone p_{T}/jet p_{T} vs. jet radius, jet p_{T} > 10 GeV", "", 12, 0, 0.6,200, 0., 2., "Cone radius","dN^{Jets}/dR", "Ratio");
      AddHistogram2D<TH2D>("hJetProfile20GeV", "Jet profile, cone p_{T}/jet p_{T} vs. jet radius, jet p_{T} > 20 GeV", "", 12, 0, 0.6,200, 0., 2., "Cone radius","dN^{Jets}/dR", "Ratio");
      AddHistogram2D<TH2D>("hJetProfile30GeV", "Jet profile, cone p_{T}/jet p_{T} vs. jet radius, jet p_{T} > 30 GeV", "", 12, 0, 0.6,200, 0., 2., "Cone radius","dN^{Jets}/dR", "Ratio");
      AddHistogram2D<TH2D>("hJetProfile40GeV", "Jet profile, cone p_{T}/jet p_{T} vs. jet radius, jet p_{T} > 40 GeV", "", 12, 0, 0.6,200, 0., 2., "Cone radius","dN^{Jets}/dR", "Ratio");
      AddHistogram2D<TH2D>("hJetProfile50GeV", "Jet profile, cone p_{T}/jet p_{T} vs. jet radius, jet p_{T} > 50 GeV", "", 12, 0, 0.6,200, 0., 2., "Cone radius","dN^{Jets}/dR", "Ratio");
      AddHistogram2D<TH2D>("hJetProfile60GeV", "Jet profile, cone p_{T}/jet p_{T} vs. jet radius, jet p_{T} > 60 GeV", "", 12, 0, 0.6,200, 0., 2., "Cone radius","dN^{Jets}/dR", "Ratio");
      AddHistogram2D<TH2D>("hJetProfile70GeV", "Jet profile, cone p_{T}/jet p_{T} vs. jet radius, jet p_{T} > 70 GeV", "", 12, 0, 0.6,200, 0., 2., "Cone radius","dN^{Jets}/dR", "Ratio");
      SetCurrentOutputList(0);
    }
  }
  // ######## Jet track cuts
  if(fAnalyzeTrackcuts)
  {
    SetCurrentOutputList(2);

    AddCutHistogram("hCutsNumberClusters", "Trackcut histogram: Number of clusters", "Number of clusters", 40, 20, 160);
    AddCutHistogram("hCutsChi2TPC", "Trackcut histogram: #chi^{2} per TPC cluster", "#chi^{2}", 40, 0, 8);
    AddCutHistogram("hCutsChi2ITS", "Trackcut histogram: #chi^{2} per ITS cluster", "#chi^{2}", 25, 0., 50);
    AddCutHistogram("hCutsChi2Constrained", "Trackcut histogram: #chi^{2} for global constrained tracks", "#chi^{2}", 60, 0, 60);
    AddCutHistogram("hCutsDCAXY", "Trackcut histogram: Max. DCA xy for prim. vertex", "DCA xy", 20, 0, 4);
    AddCutHistogram("hCutsDCAZ", "Trackcut histogram: Max. DCA z for prim. vertex", "DCA z", 20, 0, 4);
    AddCutHistogram("hCutsSPDHit", "Trackcut histogram: Hit in SPD layer", "Hit or not", 2, -0.5, 1.5);
    AddCutHistogram("hCutsNumberCrossedRows", "Trackcut histogram: Number of crossed rows", "Number of crossed rows", 40, 20, 160);
    AddCutHistogram("hCutsNumberCrossedRowsOverFindableClusters", "Trackcut histogram: Number of crossed rows over findable clusters", "Number of crossed rows over findable clusters", 26, 0.4, 1.8);
    AddCutHistogram("hCutsSharedTPC", "Trackcut histogram: Shared TPC clusters", "Shared fraction", 40, 0, 1);
    AddCutHistogram("hCutsTPCRefit", "Trackcut histogram: TPC refit", "Has TPC refit", 2, -0.5, 1.5);
    AddCutHistogram("hCutsAcceptKinks", "Trackcut histogram: Kink in track", "Kink in track", 2, -0.5, 1.5);
    AddCutHistogram("hCutsTPCLength", "Trackcut histogram: TPC length", "TPC length", 40, 0, 170);
    AddCutHistogram("hCutsTrackConstrained", "Trackcut histogram: Tracks constrained to vertex", "Track is constrained", 2, -0.5, 1.5);
    AddCutHistogram("hCutsTPCITSMatching", "Trackcut histogram: TPC-ITS matching", "Track is matched", 2, -0.5, 1.5);
    AddCutHistogram("hCutsClustersPtDependence", "Trackcut histogram: pT dependence for number of clusters/crossed rows cut.", "Value at 20 GeV: 90, 100, 110, or 120", 4, -0.5, 3.5);

    const int nbPt=100;
    const double ptMax=50;
    AddHistogram2D<TH2D>("hCutsITSTPC_NMatch", "Number matches", "", nbPt,0,ptMax,kMaxMatch+1,-0.5,kMaxMatch+0.5, "p_{T}","N matches");
    AddHistogram2D<TH2D>("hCutsITSTPC_BestMatch", "Best match chi2", "", nbPt,0,ptMax,2*int(TMath::Max(1.1,kMaxChi2)),0,kMaxChi2, "p_{T}","chi2");
    AddHistogram2D<TH2D>("hCutsITSTPC_BestMatch_cuts", "Best match chi2", "", nbPt,0,ptMax,2*int(TMath::Max(1.1,kMaxChi2)),0,kMaxChi2, "p_{T}","chi2");
    AddHistogram2D<TH2D>("hCutsITSTPC_AllMatch", "All matches chi2", "", nbPt,0,ptMax,2*int(TMath::Max(1.1,kMaxChi2)),0,kMaxChi2, "p_{T}","chi2");
    AddHistogram2D<TH2D>("hCutsITSTPC_AllMatchGlo", "All matches chi2", "", nbPt,0,ptMax,2*int(TMath::Max(1.1,kMaxChi2)),0,kMaxChi2, "p_{T}","chi2");
    AddHistogram2D<TH2D>("hCutsITSTPC_PtCorr_ITSTPC", "PtCorr", "", nbPt,0,ptMax,nbPt,0,ptMax, "p_{T}","p_{T}");
    AddHistogram2D<TH2D>("hCutsITSTPC_dPtRel_ITSTPC", "dPt/pt", "", nbPt,0,ptMax,2*nbPt+1,-0.4*ptMax,0.4*ptMax, "p_{T}","1/pt");
    AddHistogram2D<TH2D>("hCutsITSTPC_dInvPtRel_ITSTPC", "pt*dPt^{-1}", "", nbPt,0,ptMax,2*nbPt+1,-0.4*ptMax,0.4*ptMax, "p_{T}","1/pt");

    AddHistogram2D<TH2D>("hCutsITSTPC_NMatchBg", "Number matches", "", nbPt,0,ptMax,kMaxMatch+1,-0.5,kMaxMatch+0.5, "p_{T}","N matches");
    AddHistogram2D<TH2D>("hCutsITSTPC_BestMatchBg", "Best match chi2", "", nbPt,0,ptMax,2*int(TMath::Max(1.1,kMaxChi2)),0,kMaxChi2, "p_{T}","chi2");
    AddHistogram2D<TH2D>("hCutsITSTPC_BestMatchBg_cuts", "Best match chi2", "", nbPt,0,ptMax,2*int(TMath::Max(1.1,kMaxChi2)),0,kMaxChi2, "p_{T}","chi2");
    AddHistogram2D<TH2D>("hCutsITSTPC_AllMatchBg", "All matches chi2", "", nbPt,0,ptMax,2*int(TMath::Max(1.1,kMaxChi2)),0,kMaxChi2, "p_{T}","chi2");
    AddHistogram2D<TH2D>("hCutsITSTPC_AllMatchGloBg", "All matches chi2", "", nbPt,0,ptMax,2*int(TMath::Max(1.1,kMaxChi2)),0,kMaxChi2, "p_{T}","chi2");
    AddHistogram2D<TH2D>("hCutsITSTPC_PtCorrBg_ITSTPC", "PtCorr", "", nbPt,0,ptMax,nbPt,0,ptMax, "p_{T}","p_{T}");
    AddHistogram2D<TH2D>("hCutsITSTPC_dPtRelBg_ITSTPC", "dPt/pt", "", nbPt,0,ptMax,2*nbPt+1,-0.4*ptMax,0.4*ptMax, "p_{T}","1/pt");
    AddHistogram2D<TH2D>("hCutsITSTPC_dInvPtRelBg_ITSTPC", "pt*dPt^{-1}", "", nbPt,0,ptMax,2*nbPt+1,-0.4*ptMax,0.4*ptMax, "p_{T}","1/pt");

    SetCurrentOutputList(0);
  }

  PostData(1, fOutputLists[0]);
  if(fAnalyzeJetProfile)
    PostData(2, fOutputLists[1]);
  if(fAnalyzeTrackcuts)
  {
    if(fAnalyzeJetProfile)
      PostData(3, fOutputLists[2]);
    else
      PostData(2, fOutputLists[1]);
  }

}

//________________________________________________________________________
AliAnalysisTaskChargedJetsPA::AliAnalysisTaskChargedJetsPA(const char *name, const char* trackArrayName, const char* jetArrayName, const char* backgroundJetArrayName, Bool_t analyzeJetProfile, Bool_t analyzeTrackcuts) : AliAnalysisTaskSE(name), fOutputLists(), fCurrentOutputList(0), fDoJetAnalysis(1), fAnalyzeJetProfile(0), fAnalyzeTrackcuts(0), fAnalyzeJetConstituents(1), fParticleLevel(0), fUseDefaultVertexCut(1), fUsePileUpCut(1), fSetCentralityToOne(0), fNoExternalBackground(0), fBackgroundForJetProfile(0), fPartialAnalysisNParts(1), fPartialAnalysisIndex(0), fJetArray(0), fTrackArray(0), fBackgroundJetArray(0), fJetArrayName(), fTrackArrayName(), fBackgroundJetArrayName(), fRhoTaskName(), fRandConeRadius(0.4), fRandConeNumber(10), fSignalJetRadius(0.4), fBackgroundJetRadius(0.4), fNumberExcludedJets(-1), fMinEta(-0.9), fMaxEta(0.9), fMinJetEta(-0.5), fMaxJetEta(0.5), fMinTrackPt(0.150), fMinJetPt(5.0), fMinJetArea(0.5), fMinBackgroundJetPt(0.0), fMinNCrossedRows(70), fUsePtDepCrossedRowsCut(0), fNumberOfCentralityBins(20), fCentralityType("V0A"), fMatchTr(), fMatchChi(), fPrimaryVertex(0), fFirstLeadingJet(0), fSecondLeadingJet(0), fFirstLeadingKTJet(0), fSecondLeadingKTJet(0), fNumberSignalJets(0), fNumberSignalJetsAbove5GeV(0), fRandom(0), fHelperClass(0), fInitialized(0), fTaskInstanceCounter(0), fIsDEBUG(0), fIsPA(1), fNoTerminate(1), fEventCounter(0), fTempExcludedRCs(0), fTempAllRCs(1), fTempOverlapCounter(0), fTempMeanExclusionProbability(0), fHybridESDtrackCuts(0), fHybridESDtrackCuts_variedPtDep(0), fHybridESDtrackCuts_variedPtDep2(0)
{
  #ifdef DEBUGMODE
    AliInfo("Calling constructor.");
  #endif

  // Every instance of this task gets his own number
  static Int_t instance = 0;
  fTaskInstanceCounter = instance;
  instance++;

  fAnalyzeJetProfile = analyzeJetProfile;
  fAnalyzeTrackcuts  = analyzeTrackcuts; 

  // Save the observables array names
  fTrackArrayName  = trackArrayName;
  fJetArrayName = jetArrayName;
  fBackgroundJetArrayName = backgroundJetArrayName;

  if (fTrackArrayName.Contains("MCParticles") || fTrackArrayName.Contains("mcparticles"))
    fParticleLevel = kTRUE;

  DefineOutput(1, TList::Class());
  if(fAnalyzeJetProfile)
    DefineOutput(2, TList::Class());
  if(fAnalyzeTrackcuts)
  {
    if(fAnalyzeJetProfile)
      DefineOutput(3, TList::Class());
    else
      DefineOutput(2, TList::Class());
  }

  #ifdef DEBUGMODE
    AliInfo("Constructor done.");
  #endif
}

//________________________________________________________________________
void AliAnalysisTaskChargedJetsPA::InitializeTrackcuts()
{
  AliESDtrackCuts* commonTrackCuts = new AliESDtrackCuts;
  commonTrackCuts->SetMaxChi2PerClusterTPC(4);
  commonTrackCuts->SetMaxChi2PerClusterITS(36);
  commonTrackCuts->SetAcceptKinkDaughters(kFALSE);
  commonTrackCuts->SetRequireTPCRefit(kTRUE);
  commonTrackCuts->SetRequireITSRefit(kTRUE);
  commonTrackCuts->SetRequireSigmaToVertex(kFALSE);
  commonTrackCuts->SetMaxDCAToVertexXY(2.4);
  commonTrackCuts->SetMaxDCAToVertexZ(3.2);
  commonTrackCuts->SetDCAToVertex2D(kTRUE);
  commonTrackCuts->SetMaxFractionSharedTPCClusters(0.4);
  commonTrackCuts->SetMaxChi2TPCConstrainedGlobal(36);
  commonTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kAny);

  AliESDtrackCuts*    fTrackCutsPA_global = NULL;
  AliESDtrackCuts*    fTrackCutsPA_complementary = NULL;
  AliESDtrackCuts*    fTrackCutsPP_global = NULL;
  AliESDtrackCuts*    fTrackCutsPP_complementary = NULL;
  AliESDtrackCuts*    fTrackCutsPP_global_variedPtDep = NULL;
  AliESDtrackCuts*    fTrackCutsPP_complementary_variedPtDep = NULL;
  AliESDtrackCuts*    fTrackCutsPP_global_variedPtDep2 = NULL;
  AliESDtrackCuts*    fTrackCutsPP_complementary_variedPtDep2 = NULL;

  //pPb
  fTrackCutsPA_global = static_cast<AliESDtrackCuts*>(commonTrackCuts->Clone("fTrackCutsPA_global"));
  fTrackCutsPA_global->SetMinNCrossedRowsTPC(fMinNCrossedRows);
  fTrackCutsPA_global->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
  fTrackCutsPA_complementary = static_cast<AliESDtrackCuts*>(fTrackCutsPA_global->Clone("fTrackCutsPA_complementary"));
  fTrackCutsPA_complementary->SetRequireITSRefit(kFALSE);
  fTrackCutsPA_complementary->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kOff);

  //pp 
  fTrackCutsPP_global = static_cast<AliESDtrackCuts*>(commonTrackCuts->Clone("fTrackCutsPP_global"));
  TFormula *f1NClustersTPCLinearPtDep = new TFormula("f1NClustersTPCLinearPtDep","70.+30./20.*x");
  fTrackCutsPP_global->SetMinNClustersTPCPtDep(f1NClustersTPCLinearPtDep,20.);
  fTrackCutsPP_global->SetMinNClustersTPC(70);
  fTrackCutsPP_global->SetRequireTPCStandAlone(kTRUE); //cut on NClustersTPC and chi2TPC Iter1
  fTrackCutsPP_global->SetEtaRange(-0.9,0.9);
  fTrackCutsPP_global->SetPtRange(0.15, 1e15);
  fTrackCutsPP_complementary = static_cast<AliESDtrackCuts*>(fTrackCutsPP_global->Clone("fTrackCutsPP_complementary"));
  fTrackCutsPP_complementary->SetRequireITSRefit(kFALSE);
  fTrackCutsPP_complementary->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kOff);

  //pp, different pT dependence of number clusters cut, No. I

  fTrackCutsPP_global_variedPtDep = static_cast<AliESDtrackCuts*>(commonTrackCuts->Clone("fTrackCutsPP_global_variedPtDep"));
  TFormula *f1NClustersTPCLinearPtDep2 = new TFormula("f1NClustersTPCLinearPtDep2","70.+15./20.*x");
  fTrackCutsPP_global_variedPtDep->SetMinNClustersTPCPtDep(f1NClustersTPCLinearPtDep2,20.);
  fTrackCutsPP_global_variedPtDep->SetMinNClustersTPC(70);
  fTrackCutsPP_global_variedPtDep->SetRequireTPCStandAlone(kTRUE); //cut on NClustersTPC and chi2TPC Iter1
  fTrackCutsPP_global_variedPtDep->SetEtaRange(-0.9,0.9);
  fTrackCutsPP_global_variedPtDep->SetPtRange(0.15, 1e15);
  fTrackCutsPP_complementary_variedPtDep = static_cast<AliESDtrackCuts*>(fTrackCutsPP_global_variedPtDep->Clone("fTrackCutsPP_complementary_variedPtDep"));
  fTrackCutsPP_complementary_variedPtDep->SetRequireITSRefit(kFALSE);
  fTrackCutsPP_complementary_variedPtDep->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kOff);

  //pp, different pT dependence of number clusters cut, No. II

  fTrackCutsPP_global_variedPtDep2 = static_cast<AliESDtrackCuts*>(commonTrackCuts->Clone("fTrackCutsPP_global_variedPtDep2"));
  TFormula *f1NClustersTPCLinearPtDep3 = new TFormula("f1NClustersTPCLinearPtDep3","70.+45./20.*x");
  fTrackCutsPP_global_variedPtDep2->SetMinNClustersTPCPtDep(f1NClustersTPCLinearPtDep3,20.);
  fTrackCutsPP_global_variedPtDep2->SetMinNClustersTPC(70);
  fTrackCutsPP_global_variedPtDep2->SetRequireTPCStandAlone(kTRUE); //cut on NClustersTPC and chi2TPC Iter1
  fTrackCutsPP_global_variedPtDep2->SetEtaRange(-0.9,0.9);
  fTrackCutsPP_global_variedPtDep2->SetPtRange(0.15, 1e15);
  fTrackCutsPP_complementary_variedPtDep2 = static_cast<AliESDtrackCuts*>(fTrackCutsPP_global_variedPtDep2->Clone("fTrackCutsPP_complementary_variedPtDep2"));
  fTrackCutsPP_complementary_variedPtDep2->SetRequireITSRefit(kFALSE);
  fTrackCutsPP_complementary_variedPtDep2->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kOff);

  fHybridESDtrackCuts = new AliESDHybridTrackcuts();
  if(fIsPA)
  {
    fHybridESDtrackCuts->SetMainCuts(fTrackCutsPA_global);
    fHybridESDtrackCuts->SetAdditionalCuts(fTrackCutsPA_complementary);
  }
  else
  {
    fHybridESDtrackCuts_variedPtDep = new AliESDHybridTrackcuts();
    fHybridESDtrackCuts_variedPtDep2 = new AliESDHybridTrackcuts();

    fHybridESDtrackCuts->SetMainCuts(fTrackCutsPP_global);
    fHybridESDtrackCuts->SetAdditionalCuts(fTrackCutsPP_complementary);
    fHybridESDtrackCuts_variedPtDep->SetMainCuts(fTrackCutsPP_global_variedPtDep);
    fHybridESDtrackCuts_variedPtDep->SetAdditionalCuts(fTrackCutsPP_complementary_variedPtDep);
    fHybridESDtrackCuts_variedPtDep2->SetMainCuts(fTrackCutsPP_global_variedPtDep2);
    fHybridESDtrackCuts_variedPtDep2->SetAdditionalCuts(fTrackCutsPP_complementary_variedPtDep2);
  }

  delete commonTrackCuts;
}

//________________________________________________________________________
void AliAnalysisTaskChargedJetsPA::CreateCutHistograms()
{

  AliESDEvent* fESD = dynamic_cast<AliESDEvent*>( InputEvent() );
  if (!fESD)
  {
    AliError("For cut analysis, ESDs must be processed!");
    return;
  }

  SetCurrentOutputList(2);

  Float_t dca[2], cov[3]; // dca_xy, dca_z, sigma_xy, sigma_xy_z, sigma_z for the vertex cut
  for (Int_t i=0;i < fESD->GetNumberOfTracks(); i++)
  {
    AliESDtrack* track = fESD->GetTrack(i);

    // Basics kinematic variables
    Double_t pT                  = track->Pt();
    Double_t eta                 = track->Eta();
    Double_t phi                 = track->Phi();

    // Number of clusters
    Double_t nclsITS             = track->GetITSclusters(0);

    // Crossed rows
    Double_t ncrTPC              = track->GetTPCCrossedRows();
    Double_t nCRoverFC           = 0;
    if(track->GetTPCNclsF())
      nCRoverFC = track->GetTPCCrossedRows()/track->GetTPCNclsF();

    // Chi2 of tracks
    Double_t chi2ITS             = 999.; 
    if (nclsITS)
      chi2ITS = track->GetITSchi2()/nclsITS;
    Double_t chi2TPC             = 999.;
    Double_t chi2TPCConstrained  = track->GetChi2TPCConstrainedVsGlobal(static_cast<const AliESDVertex*>(fPrimaryVertex));

    // Misc
    Double_t SharedTPCClusters = 999.;
    Double_t nClustersTPC = 0;

    if(fHybridESDtrackCuts->GetMainCuts()->GetRequireTPCStandAlone())
    {
      nClustersTPC = track->GetTPCNclsIter1();
      if(nClustersTPC)
        chi2TPC = track->GetTPCchi2Iter1()/nClustersTPC;
    }
    else 
    {
      nClustersTPC = track->GetTPCclusters(0);
      if(nClustersTPC)
        chi2TPC = track->GetTPCchi2()/nClustersTPC;
    }

    if(nClustersTPC)
      SharedTPCClusters = static_cast<Double_t>(track->GetTPCnclsS())/static_cast<Double_t>(nClustersTPC);


    Double_t tpcLength   = 0.;
    if (track->GetInnerParam() && track->GetESDEvent()) {
      tpcLength = track->GetLengthInActiveZone(1, 1.8, 220, track->GetESDEvent()->GetMagneticField());
    }
    track->GetImpactParameters(dca, cov);

    // Basic kinematic cuts
    if((pT<0.15) || (TMath::Abs(eta)>0.9))
      continue;

    Int_t trackType = 0;

    // ################################################################
    // ################################################################

    if(fIsPA)
    {
      trackType = fHybridESDtrackCuts->AcceptTrack(track);
      Double_t tmpThreshold90  = 70. + 20./20. * pT;
      Double_t tmpThreshold100 = 70. + 30./20. * pT;
      Double_t tmpThreshold110 = 70. + 40./20. * pT;
      Double_t tmpThreshold120 = 70. + 50./20. * pT;

      if(pT>20.)
      {
        tmpThreshold90 = 70. + 20.;
        tmpThreshold100 = 70. + 30.;
        tmpThreshold110 = 70. + 40.;
        tmpThreshold120 = 70. + 50.;
      }

      if (trackType)
      {
        if(ncrTPC>=tmpThreshold90)
          FillCutHistogram("hCutsClustersPtDependence", 0, pT, eta, phi, trackType-1);
        if(ncrTPC>=tmpThreshold100)
          FillCutHistogram("hCutsClustersPtDependence", 1, pT, eta, phi, trackType-1);
        if(ncrTPC>=tmpThreshold110)
          FillCutHistogram("hCutsClustersPtDependence", 2, pT, eta, phi, trackType-1);
        if(ncrTPC>=tmpThreshold120)
          FillCutHistogram("hCutsClustersPtDependence", 3, pT, eta, phi, trackType-1);
      }

      if(fUsePtDepCrossedRowsCut && (ncrTPC<tmpThreshold100)) // pT dep crossed rows cut is not fulfilled
        continue; // next track
    }
    else
    {
      trackType = fHybridESDtrackCuts_variedPtDep->AcceptTrack(track);
      if (trackType)
        FillCutHistogram("hCutsClustersPtDependence", 0, pT, eta, phi, trackType-1);

      trackType = fHybridESDtrackCuts->AcceptTrack(track);
      if (trackType)
        FillCutHistogram("hCutsClustersPtDependence", 1, pT, eta, phi, trackType-1);

      trackType = fHybridESDtrackCuts_variedPtDep2->AcceptTrack(track);
      if (trackType)
        FillCutHistogram("hCutsClustersPtDependence", 2, pT, eta, phi, trackType-1);
    }

    // ################################################################
    // ################################################################
    Int_t minNclsTPC = fHybridESDtrackCuts->GetMainCuts()->GetMinNClusterTPC();
    Int_t minNclsTPC_Additional = fHybridESDtrackCuts->GetAdditionalCuts()->GetMinNClusterTPC();
    fHybridESDtrackCuts->GetMainCuts()->SetMinNClustersTPC(0);
    fHybridESDtrackCuts->GetAdditionalCuts()->SetMinNClustersTPC(0);

    trackType = fHybridESDtrackCuts->AcceptTrack(track);
    if (trackType)
      FillCutHistogram("hCutsNumberClusters", nClustersTPC, pT, eta, phi, trackType-1);

    fHybridESDtrackCuts->GetMainCuts()->SetMinNClustersTPC(minNclsTPC);
    fHybridESDtrackCuts->GetAdditionalCuts()->SetMinNClustersTPC(minNclsTPC_Additional);
    // ################################################################
    // ################################################################
    Float_t maxChi2 = fHybridESDtrackCuts->GetMainCuts()->GetMaxChi2PerClusterTPC();
    Float_t maxChi2_Additional = fHybridESDtrackCuts->GetAdditionalCuts()->GetMaxChi2PerClusterTPC();
    fHybridESDtrackCuts->GetMainCuts()->SetMaxChi2PerClusterTPC(999.);
    fHybridESDtrackCuts->GetAdditionalCuts()->SetMaxChi2PerClusterTPC(999.);

    trackType = fHybridESDtrackCuts->AcceptTrack(track);
    if (trackType)
      FillCutHistogram("hCutsChi2TPC", chi2TPC, pT, eta, phi, trackType-1);

    fHybridESDtrackCuts->GetMainCuts()->SetMaxChi2PerClusterTPC(maxChi2);
    fHybridESDtrackCuts->GetAdditionalCuts()->SetMaxChi2PerClusterTPC(maxChi2_Additional);

    // ################################################################
    // ################################################################
    Float_t maxChi2TPCConstrained = fHybridESDtrackCuts->GetMainCuts()->GetMaxChi2TPCConstrainedGlobal();
    Float_t maxChi2TPCConstrained_Additional = fHybridESDtrackCuts->GetAdditionalCuts()->GetMaxChi2TPCConstrainedGlobal();
    fHybridESDtrackCuts->GetMainCuts()->SetMaxChi2TPCConstrainedGlobal(999.);
    fHybridESDtrackCuts->GetAdditionalCuts()->SetMaxChi2TPCConstrainedGlobal(999.);

    trackType = fHybridESDtrackCuts->AcceptTrack(track);
    if (trackType)
      FillCutHistogram("hCutsChi2Constrained", chi2TPCConstrained, pT, eta, phi, trackType-1);

    fHybridESDtrackCuts->GetMainCuts()->SetMaxChi2TPCConstrainedGlobal(maxChi2TPCConstrained);
    fHybridESDtrackCuts->GetAdditionalCuts()->SetMaxChi2TPCConstrainedGlobal(maxChi2TPCConstrained_Additional);

    // ################################################################
    // ################################################################
    Float_t maxDcaZ = fHybridESDtrackCuts->GetMainCuts()->GetMaxDCAToVertexZ();
    Float_t maxDcaZ_Additional = fHybridESDtrackCuts->GetAdditionalCuts()->GetMaxDCAToVertexZ();
    fHybridESDtrackCuts->GetMainCuts()->SetMaxDCAToVertexZ(999.);
    fHybridESDtrackCuts->GetAdditionalCuts()->SetMaxDCAToVertexZ(999.);

    trackType = fHybridESDtrackCuts->AcceptTrack(track);
    if (trackType)
      FillCutHistogram("hCutsDCAZ", TMath::Abs(dca[1]), pT, eta, phi, trackType-1);

    fHybridESDtrackCuts->GetMainCuts()->SetMaxDCAToVertexZ(maxDcaZ);
    fHybridESDtrackCuts->GetAdditionalCuts()->SetMaxDCAToVertexZ(maxDcaZ_Additional);

    // ################################################################
    // ################################################################
    Float_t maxDcaXY = fHybridESDtrackCuts->GetMainCuts()->GetMaxDCAToVertexXY();
    Float_t maxDcaXY_Additional = fHybridESDtrackCuts->GetAdditionalCuts()->GetMaxDCAToVertexXY();
    fHybridESDtrackCuts->GetMainCuts()->SetMaxDCAToVertexXY(999.);
    fHybridESDtrackCuts->GetAdditionalCuts()->SetMaxDCAToVertexXY(999.);

    trackType = fHybridESDtrackCuts->AcceptTrack(track);
    if (trackType)
      FillCutHistogram("hCutsDCAXY", TMath::Abs(dca[0]), pT, eta, phi, trackType-1);

    fHybridESDtrackCuts->GetMainCuts()->SetMaxDCAToVertexXY(maxDcaXY);
    fHybridESDtrackCuts->GetAdditionalCuts()->SetMaxDCAToVertexXY(maxDcaXY_Additional);

    // ################################################################
    // ################################################################
    AliESDtrackCuts::ITSClusterRequirement clusterReq = fHybridESDtrackCuts->GetMainCuts()->GetClusterRequirementITS(AliESDtrackCuts::kSPD);
    AliESDtrackCuts::ITSClusterRequirement clusterReq_Additional = fHybridESDtrackCuts->GetAdditionalCuts()->GetClusterRequirementITS(AliESDtrackCuts::kSPD);
    fHybridESDtrackCuts->GetMainCuts()->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kOff);
    fHybridESDtrackCuts->GetAdditionalCuts()->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kOff);

    Int_t hasPoint = 0;
    if (track->HasPointOnITSLayer(0) || track->HasPointOnITSLayer(1)) hasPoint = 1;
    trackType = fHybridESDtrackCuts->AcceptTrack(track);
    if (trackType)
      FillCutHistogram("hCutsSPDHit", hasPoint, pT, eta, phi, trackType-1);

    fHybridESDtrackCuts->GetMainCuts()->SetClusterRequirementITS(AliESDtrackCuts::kSPD, clusterReq);
    fHybridESDtrackCuts->GetAdditionalCuts()->SetClusterRequirementITS(AliESDtrackCuts::kSPD, clusterReq_Additional);

    // ################################################################
    // ################################################################
    Float_t minNcrTPC = fHybridESDtrackCuts->GetMainCuts()->GetMinNCrossedRowsTPC();
    Float_t minNcrTPC_Additional = fHybridESDtrackCuts->GetAdditionalCuts()->GetMinNCrossedRowsTPC();
    fHybridESDtrackCuts->GetMainCuts()->SetMinNCrossedRowsTPC(0);
    fHybridESDtrackCuts->GetAdditionalCuts()->SetMinNCrossedRowsTPC(0);

    trackType = fHybridESDtrackCuts->AcceptTrack(track);
    if (trackType)
      FillCutHistogram("hCutsNumberCrossedRows", ncrTPC, pT, eta, phi, trackType-1);

    fHybridESDtrackCuts->GetMainCuts()->SetMinNCrossedRowsTPC(minNcrTPC);
    fHybridESDtrackCuts->GetAdditionalCuts()->SetMinNCrossedRowsTPC(minNcrTPC_Additional);

    // ################################################################
    // ################################################################
    Float_t minCRoverFC = fHybridESDtrackCuts->GetMainCuts()->GetMinRatioCrossedRowsOverFindableClustersTPC();
    Float_t minCRoverFC_Additional = fHybridESDtrackCuts->GetAdditionalCuts()->GetMinRatioCrossedRowsOverFindableClustersTPC();
    fHybridESDtrackCuts->GetMainCuts()->SetMinRatioCrossedRowsOverFindableClustersTPC(0.);
    fHybridESDtrackCuts->GetAdditionalCuts()->SetMinRatioCrossedRowsOverFindableClustersTPC(0.);

    trackType = fHybridESDtrackCuts->AcceptTrack(track);
    if (trackType)
      FillCutHistogram("hCutsNumberCrossedRowsOverFindableClusters", nCRoverFC, pT, eta, phi, trackType-1);

    fHybridESDtrackCuts->GetMainCuts()->SetMinRatioCrossedRowsOverFindableClustersTPC(minCRoverFC);
    fHybridESDtrackCuts->GetAdditionalCuts()->SetMinRatioCrossedRowsOverFindableClustersTPC(minCRoverFC_Additional);

    // ################################################################
    // ################################################################
    Float_t maxSharedTPC = fHybridESDtrackCuts->GetMainCuts()->GetMaxFractionSharedTPCClusters();
    Float_t maxSharedTPC_Additional = fHybridESDtrackCuts->GetAdditionalCuts()->GetMaxFractionSharedTPCClusters();
    fHybridESDtrackCuts->GetMainCuts()->SetMaxFractionSharedTPCClusters(999.);
    fHybridESDtrackCuts->GetAdditionalCuts()->SetMaxFractionSharedTPCClusters(999.);

    trackType = fHybridESDtrackCuts->AcceptTrack(track);
    if (trackType)
      FillCutHistogram("hCutsSharedTPC", SharedTPCClusters, pT, eta, phi, trackType-1);

    fHybridESDtrackCuts->GetMainCuts()->SetMaxFractionSharedTPCClusters(maxSharedTPC);
    fHybridESDtrackCuts->GetAdditionalCuts()->SetMaxFractionSharedTPCClusters(maxSharedTPC_Additional);

    // ################################################################
    // ################################################################
    Bool_t reqTPCRefit = fHybridESDtrackCuts->GetMainCuts()->GetRequireTPCRefit();
    Bool_t reqTPCRefit_Additional = fHybridESDtrackCuts->GetAdditionalCuts()->GetRequireTPCRefit();

    fHybridESDtrackCuts->GetMainCuts()->SetRequireTPCRefit(1);
    fHybridESDtrackCuts->GetAdditionalCuts()->SetRequireTPCRefit(1);
    trackType = fHybridESDtrackCuts->AcceptTrack(track);
    if (trackType)
      FillCutHistogram("hCutsTPCRefit", 1, pT, eta, phi, trackType-1);
    else // track is not accepted as global hybrid with TPC refit requirement
    {
      fHybridESDtrackCuts->GetMainCuts()->SetRequireTPCRefit(0);
      fHybridESDtrackCuts->GetAdditionalCuts()->SetRequireTPCRefit(0);
      trackType = fHybridESDtrackCuts->AcceptTrack(track);
      if (trackType)
        FillCutHistogram("hCutsTPCRefit", 0, pT, eta, phi, trackType-1);
    }
    fHybridESDtrackCuts->GetMainCuts()->SetRequireTPCRefit(reqTPCRefit);
    fHybridESDtrackCuts->GetAdditionalCuts()->SetRequireTPCRefit(reqTPCRefit_Additional);

    // ################################################################
    // ################################################################
    Bool_t accKinks = fHybridESDtrackCuts->GetMainCuts()->GetAcceptKinkDaughters();
    Bool_t accKinks_Additional = fHybridESDtrackCuts->GetAdditionalCuts()->GetAcceptKinkDaughters();

    fHybridESDtrackCuts->GetMainCuts()->SetAcceptKinkDaughters(0);
    fHybridESDtrackCuts->GetAdditionalCuts()->SetAcceptKinkDaughters(0);
    trackType = fHybridESDtrackCuts->AcceptTrack(track);
    if (trackType) // A passing track has no kinks
      FillCutHistogram("hCutsAcceptKinks", 0, pT, eta, phi, trackType-1);
    else
    {
      fHybridESDtrackCuts->GetMainCuts()->SetAcceptKinkDaughters(1);
      fHybridESDtrackCuts->GetAdditionalCuts()->SetAcceptKinkDaughters(1);
      trackType = fHybridESDtrackCuts->AcceptTrack(track);
      if (trackType) // A passing track has kinks
        FillCutHistogram("hCutsAcceptKinks", 1, pT, eta, phi, trackType-1);
    }
    fHybridESDtrackCuts->GetMainCuts()->SetAcceptKinkDaughters(accKinks);
    fHybridESDtrackCuts->GetAdditionalCuts()->SetAcceptKinkDaughters(accKinks_Additional);

    // ################################################################
    // ################################################################
    Float_t maxChi2ITS = fHybridESDtrackCuts->GetMainCuts()->GetMaxChi2PerClusterITS();
    Float_t maxChi2ITS_Additional = fHybridESDtrackCuts->GetAdditionalCuts()->GetMaxChi2PerClusterITS();
    fHybridESDtrackCuts->GetMainCuts()->SetMaxChi2PerClusterITS(999.);
    fHybridESDtrackCuts->GetAdditionalCuts()->SetMaxChi2PerClusterITS(999.);

    trackType = fHybridESDtrackCuts->AcceptTrack(track);
    if (trackType)
      FillCutHistogram("hCutsChi2ITS", chi2ITS, pT, eta, phi, trackType-1);

    fHybridESDtrackCuts->GetMainCuts()->SetMaxChi2PerClusterITS(maxChi2ITS);
    fHybridESDtrackCuts->GetAdditionalCuts()->SetMaxChi2PerClusterITS(maxChi2ITS_Additional);

    // ################################################################
    // ################################################################
    Float_t minTpcLength = fHybridESDtrackCuts->GetMainCuts()->GetMinLengthActiveVolumeTPC(); // Active length TPC
    Float_t minTpcLength_Additional = fHybridESDtrackCuts->GetAdditionalCuts()->GetMinLengthActiveVolumeTPC();
    fHybridESDtrackCuts->GetMainCuts()->SetMinLengthActiveVolumeTPC(0);
    fHybridESDtrackCuts->GetAdditionalCuts()->SetMinLengthActiveVolumeTPC(0);

    trackType = fHybridESDtrackCuts->AcceptTrack(track);
    if (trackType)
      FillCutHistogram("hCutsTPCLength", tpcLength, pT, eta, phi, trackType-1);

    fHybridESDtrackCuts->GetMainCuts()->SetMinLengthActiveVolumeTPC(minTpcLength);
    fHybridESDtrackCuts->GetAdditionalCuts()->SetMinLengthActiveVolumeTPC(minTpcLength_Additional);

    // ################################################################
    // ################################################################
    Bool_t isMatched = kFALSE;
    Float_t chi2tpc = fHybridESDtrackCuts->GetMainCuts()->GetMaxChi2TPCConstrainedGlobal();
    Float_t chi2its = fHybridESDtrackCuts->GetMainCuts()->GetMaxChi2PerClusterITS();
    Float_t chi2tpc_Additional = fHybridESDtrackCuts->GetAdditionalCuts()->GetMaxChi2TPCConstrainedGlobal();
    Float_t chi2its_Additional = fHybridESDtrackCuts->GetAdditionalCuts()->GetMaxChi2PerClusterITS();
    
    fHybridESDtrackCuts->GetMainCuts()->SetMaxChi2TPCConstrainedGlobal(99999.);
    fHybridESDtrackCuts->GetMainCuts()->SetMaxChi2PerClusterITS(999999.);
    fHybridESDtrackCuts->GetAdditionalCuts()->SetMaxChi2TPCConstrainedGlobal(99999.);
    fHybridESDtrackCuts->GetAdditionalCuts()->SetMaxChi2PerClusterITS(999999.);
    
    trackType = fHybridESDtrackCuts->AcceptTrack(track);
    if (trackType)
      FillCutHistogram("hCutsTPCITSMatching", isMatched, pT, eta, phi, trackType-1);

    fHybridESDtrackCuts->GetMainCuts()->SetMaxChi2TPCConstrainedGlobal(chi2tpc);
    fHybridESDtrackCuts->GetMainCuts()->SetMaxChi2PerClusterITS(chi2its);
    fHybridESDtrackCuts->GetAdditionalCuts()->SetMaxChi2TPCConstrainedGlobal(chi2tpc_Additional);
    fHybridESDtrackCuts->GetAdditionalCuts()->SetMaxChi2PerClusterITS(chi2its_Additional);

    isMatched=kTRUE;
    trackType = fHybridESDtrackCuts->AcceptTrack(track);
    if (trackType)
      FillCutHistogram("hCutsTPCITSMatching", isMatched, pT, eta, phi, trackType-1);

    // ################################################################
    // ################################################################
    if((fHybridESDtrackCuts->GetMainCuts()->GetClusterRequirementITS(AliESDtrackCuts::kSPD) == AliESDtrackCuts::kOff)
    || (fHybridESDtrackCuts->GetAdditionalCuts() && (fHybridESDtrackCuts->GetAdditionalCuts()->GetClusterRequirementITS(AliESDtrackCuts::kSPD) == AliESDtrackCuts::kOff))) 
    {
      Bool_t isConstrainedWithITSRefit = static_cast<Bool_t>(track->GetConstrainedParam()) && ((track->GetStatus())&AliESDtrack::kITSrefit);
      if (trackType)
        FillCutHistogram("hCutsTrackConstrained", isConstrainedWithITSRefit, pT, eta, phi, trackType-1);
    }

  }

  CreateITSTPCMatchingHistograms();
  SetCurrentOutputList(0);
}

//________________________________________________________________________
void AliAnalysisTaskChargedJetsPA::CreateITSTPCMatchingHistograms()
{
  //
  // check how many its-sa tracks get matched to TPC
  //
  Bool_t fExcludeMomFromChi2ITSTPC = kFALSE; // ITS->TPC : exclude momentum from matching chi2 calculation

  AliESDEvent* fESD = dynamic_cast<AliESDEvent*>( InputEvent() );
  if (!fESD)
  {
    AliError("For cut analysis, ESDs must be processed!");
    return;
  }

  int ntr = fESD->GetNumberOfTracks();
  //
  // initialize histograms
  //
  TH2D * hNMatch         = (TH2D*) fCurrentOutputList->FindObject("hCutsITSTPC_NMatch");
  TH2D * hBestMatch      = (TH2D*) fCurrentOutputList->FindObject("hCutsITSTPC_BestMatch");
  TH2D * hBestMatch_cuts = (TH2D*) fCurrentOutputList->FindObject("hCutsITSTPC_BestMatch_cuts");
  TH2D * hAllMatch       = (TH2D*) fCurrentOutputList->FindObject("hCutsITSTPC_AllMatch");
  TH2D * hAllMatchGlo    = (TH2D*) fCurrentOutputList->FindObject("hCutsITSTPC_AllMatchGlo");  
  TH2D * hPtCorr_ITSTPC  = (TH2D*) fCurrentOutputList->FindObject("hCutsITSTPC_PtCorr_ITSTPC");
  TH2D * hdPtRel_ITSTPC  = (TH2D*) fCurrentOutputList->FindObject("hCutsITSTPC_dPtRel_ITSTPC");
  TH2D * hdInvPtRel_ITSTPC = (TH2D*) fCurrentOutputList->FindObject("hCutsITSTPC_dInvPtRel_ITSTPC");

  //
  TH2D * hNMatchBg          = (TH2D*) fCurrentOutputList->FindObject("hCutsITSTPC_NMatchBg");
  TH2D * hBestMatchBg       = (TH2D*) fCurrentOutputList->FindObject("hCutsITSTPC_BestMatchBg");
  TH2D * hBestMatchBg_cuts  = (TH2D*) fCurrentOutputList->FindObject("hCutsITSTPC_BestMatchBg_cuts");
  TH2D * hAllMatchBg        = (TH2D*) fCurrentOutputList->FindObject("hCutsITSTPC_AllMatchBg");
  TH2D * hAllMatchGloBg     = (TH2D*) fCurrentOutputList->FindObject("hCutsITSTPC_AllMatchGloBg");    
  TH2D * hdPtRelBg_ITSTPC    = (TH2D*) fCurrentOutputList->FindObject("hCutsITSTPC_dPtRelBg_ITSTPC");
  TH2D * hdInvPtRelBg_ITSTPC = (TH2D*) fCurrentOutputList->FindObject("hCutsITSTPC_dInvPtRelBg_ITSTPC");

  if(!(hNMatch && hBestMatch && hBestMatch_cuts && hAllMatch && hAllMatchGlo && hPtCorr_ITSTPC && hdPtRel_ITSTPC && hdInvPtRel_ITSTPC && hNMatchBg && hBestMatchBg && hBestMatchBg_cuts && hAllMatchBg && hAllMatchGloBg && hdPtRelBg_ITSTPC && hdInvPtRelBg_ITSTPC))
  {
    cout << " === ERROR: At least one of the ITSTPC histograms not found! ===\n";
    cout << Form(" === Details: %p-%p-%p-%p-%p-%p-%p-%p-%p-%p-%p-%p-%p-%p-%p", hNMatch, hBestMatch, hBestMatch_cuts, hAllMatch, hAllMatchGlo, hPtCorr_ITSTPC, hdPtRel_ITSTPC, hdInvPtRel_ITSTPC, hNMatchBg, hBestMatchBg, hBestMatchBg_cuts, hAllMatchBg, hAllMatchGloBg, hdPtRelBg_ITSTPC, hdInvPtRelBg_ITSTPC) << endl;
    fCurrentOutputList->Print();
    return;    
  }
  //
  for (int it=0;it<ntr;it++) {
    AliESDtrack* trSA = fESD->GetTrack(it);
    if (!trSA->IsOn(AliESDtrack::kITSpureSA) || !trSA->IsOn(AliESDtrack::kITSrefit)) continue;
    double pt = trSA->Pt();

    // OB - fiducial eta and pt cuts
    Double_t etaSA = trSA->Eta();

    if(TMath::Abs(etaSA)>0.8) continue;

    //
    Int_t nmatch = 0;
    for (int i=kMaxMatch;i--;) {fMatchChi[i]=0; fMatchTr[i]=0;}
    for (int it1=0;it1<ntr;it1++){
      if (it1==it) continue;

      AliESDtrack* trESD = fESD->GetTrack(it1);
      if (!trESD->IsOn(AliESDtrack::kTPCrefit)) continue;

      Match(trSA,trESD, nmatch, fExcludeMomFromChi2ITSTPC);
    }
    //
    
    hNMatch->Fill(pt,nmatch);

    if (nmatch>0){
      hBestMatch->Fill(pt,fMatchChi[0]);
      hPtCorr_ITSTPC->Fill(pt,fMatchTr[0]->Pt()); 
      hdPtRel_ITSTPC->Fill(pt,(pt-fMatchTr[0]->Pt())/pt); 
      hdInvPtRel_ITSTPC->Fill(pt,pt*( 1/pt - (1/fMatchTr[0]->Pt()) )); 
    }
    
    if (nmatch>0 && fHybridESDtrackCuts){
      
      if(fHybridESDtrackCuts->AcceptTrack(fMatchTr[0])){
        hBestMatch_cuts->Fill(pt,fMatchChi[0]);
      }
    }
    
    //
    for (int imt=nmatch;imt--;) {
      hAllMatch->Fill(pt,fMatchChi[imt]);
      if (fMatchTr[imt]->IsOn(AliESDtrack::kITSrefit)) hAllMatchGlo->Fill(pt,fMatchChi[imt]);
    }
    //
    nmatch = 0;
    for (int i=kMaxMatch;i--;) {fMatchChi[i]=0; fMatchTr[i]=0;}
    for (int it1=0;it1<ntr;it1++) {
      if (it1==it) continue;
      AliESDtrack* trESD = fESD->GetTrack(it1);
      if (!trESD->IsOn(AliESDtrack::kTPCrefit)) continue;

      Match(trSA,trESD, nmatch, fExcludeMomFromChi2ITSTPC, TMath::Pi());
    }
    //

    hNMatchBg->Fill(pt,nmatch);

    if (nmatch>0){
      hBestMatchBg->Fill(pt,fMatchChi[0]);
      hdPtRelBg_ITSTPC->Fill(pt,(pt-fMatchTr[0]->Pt())/pt); 
      hdInvPtRelBg_ITSTPC->Fill(pt,pt*( 1/pt - (1/fMatchTr[0]->Pt()) )); 
    }

    if (nmatch>0 && fHybridESDtrackCuts){
      if(fHybridESDtrackCuts->AcceptTrack(fMatchTr[0])){
        hBestMatchBg_cuts->Fill(pt,fMatchChi[0]);
      }
    }

    for (int imt=nmatch;imt--;) {
      hAllMatchBg->Fill(pt,fMatchChi[imt]);
      if (fMatchTr[imt]->IsOn(AliESDtrack::kITSrefit)) hAllMatchGloBg->Fill(pt,fMatchChi[imt]);
    }
    //
  }
}

//________________________________________________________________________
void AliAnalysisTaskChargedJetsPA::Match(AliESDtrack* tr0, AliESDtrack* tr1, Int_t& nmatch, Bool_t excludeMom, Double_t rotate)
{
  //
  // check if two tracks are matching, possible rotation for combinatoric backgr.
  // 
  AliESDEvent* fESD = dynamic_cast<AliESDEvent*>( InputEvent() );
  if (!fESD)
  {
    AliError("For cut analysis, ESDs must be processed!");
    return;
  }

  Float_t bField = fESD->GetMagneticField();
  //
  const AliExternalTrackParam* trtpc0 = tr1->GetInnerParam();
  if (!trtpc0) return;
  AliExternalTrackParam trtpc(*trtpc0);
  //
  if (TMath::Abs(rotate)>1e-5) {
    const double *par = trtpc.GetParameter();
    const double *cov = trtpc.GetCovariance();
    double alp = trtpc.GetAlpha() + rotate;
    trtpc.Set(trtpc.GetX(),alp,par,cov);
  }
  //
  if (!trtpc.Rotate(tr0->GetAlpha())) return;
  if (!trtpc.PropagateTo(tr0->GetX(),bField)) return;
  double chi2 = tr0->GetPredictedChi2(&trtpc);

  //std::cout<<" in Match, nmatch "<<nmatch<<" par[4] before "<<trtpc.GetParameter()[4]<<" chi2 "<<chi2<<endl;

  // OB chi2 excluding pt 
  if(excludeMom){
    ((double*)trtpc.GetParameter())[4] = tr0->GetParameter()[4]; // set ITS mom equal TPC mom
    chi2 = tr0->GetPredictedChi2(&trtpc);

    //std::cout<<" in Match, nmatch "<<nmatch<<" par[4] after "<<trtpc.GetParameter()[4]<<" tr0 mom "<<tr0->GetParameter()[4]
    //         <<" chi2 "<<chi2<<std::endl;
  }


  if (chi2>kMaxChi2) return;

  // std::cout<<" found good match, tr1 "<<tr1<<" chi2 "<<chi2<<std::endl;
  // std::cout<<" before: fMatchChi[0]  "<<fMatchChi[0]<<" [1] "<<fMatchChi[1]
  //          <<" [2]  "<<fMatchChi[2]<<" [3] "<<fMatchChi[3]
  //          <<" [4]  "<<fMatchChi[4]<<std::endl; 

  // std::cout<<" before: fMatchTr[0]  "<<fMatchTr[0]<<" [1] "<<fMatchTr[1]
  //          <<" [2]  "<<fMatchTr[2]<<" [3] "<<fMatchTr[3]
  //          <<" [4]  "<<fMatchTr[4]<<std::endl; 

  //
  int ins;
  for (ins=0;ins<nmatch;ins++) if (chi2<fMatchChi[ins]) break;
  if (ins>=kMaxMatch) return;
  
  for (int imv=nmatch;imv>ins;imv--) {
    if (imv>=kMaxMatch) continue;
    fMatchTr[imv]  = fMatchTr[imv-1];
    fMatchChi[imv] = fMatchChi[imv-1];
  }
  fMatchTr[ins] = tr1;
  fMatchChi[ins] = chi2;
  nmatch++;
  if (nmatch>=kMaxMatch) nmatch = kMaxMatch;
  //
}

//________________________________________________________________________
Double_t AliAnalysisTaskChargedJetsPA::GetExternalRho()
{
  // Get rho from event.
  AliRhoParameter *rho = 0;
  if (!fRhoTaskName.IsNull()) {
    rho = dynamic_cast<AliRhoParameter*>(InputEvent()->FindListObject(fRhoTaskName.Data()));
    if (!rho) {
      AliWarning(Form("%s: Could not retrieve rho with name %s!", GetName(), fRhoTaskName.Data())); 
      return 0;
    }
  }
  else
    return 0;

  return (rho->GetVal());
}

//________________________________________________________________________
inline Bool_t AliAnalysisTaskChargedJetsPA::IsEventInAcceptance(AliVEvent* event)
{
  if (!event)
    return kFALSE;

  fPrimaryVertex = event->GetPrimaryVertex();

  // ### Create plot for vertex acceptance
  fHelperClass->SetMaxVtxZ(10000.);

  // Check vertex existance
  Bool_t hasVertex = kFALSE;
  if(fUseDefaultVertexCut)
  {
    if(fHelperClass->IsVertexSelected2013pA(event))
      hasVertex = kTRUE;
  }
  else
  {
    if(!fPrimaryVertex || (fPrimaryVertex->GetNContributors()<2) || (TMath::Sqrt(fPrimaryVertex->GetX()*fPrimaryVertex->GetX() + fPrimaryVertex->GetY()*fPrimaryVertex->GetY()) > 1.0)) 
      hasVertex = kTRUE;
  }

  // All triggered events
  FillHistogram("hVertexAcceptance", 0.5); // all triggered events all
  if(hasVertex)
    FillHistogram("hVertexAcceptance", 1.5); // all triggered events w/ vertex

  // Pile-up corrected events
  if(!fHelperClass->IsPileUpEvent(event))
  {
    FillHistogram("hVertexAcceptance", 2.5); // pile-up corr. events all
    if(hasVertex)
      FillHistogram("hVertexAcceptance", 3.5); // pile-up corr. events w/ vertex
  }

  fHelperClass->SetMaxVtxZ(10.);

  // ### Get centrality values
  AliCentrality* tmpCentrality = event->GetCentrality();
  Double_t centralityPercentile = -1.0;
  if (tmpCentrality != NULL)
    centralityPercentile = tmpCentrality->GetCentralityPercentile(fCentralityType.Data());
  FillHistogram("hCentrality",centralityPercentile);

  if(fSetCentralityToOne)
    centralityPercentile = 1.0;

  if((centralityPercentile < 0.0) || (centralityPercentile > 101.0))
    AliWarning(Form("Centrality value not valid (c=%E)",centralityPercentile)); 

  FillHistogram("hCentrality", centralityPercentile, 0.5); // before any cuts

  // ### CUT STAGE 1: Pile-up events (only holds for pPb)
  if(fUsePileUpCut)
    if(fHelperClass->IsPileUpEvent(event))
      return kFALSE;

  FillHistogram("hCentrality", centralityPercentile, 1.5); // after pileup cut

  // ### CUT STAGE 2: Existence of primary vertex

  if(!hasVertex)
    return kFALSE; 

  FillHistogram("hVertexZ",fPrimaryVertex->GetZ());
  FillHistogram("hCentrality", centralityPercentile, 2.5); // after vertex existance cut

  // ### CUT STAGE 3: Position of primary vertex
  if((TMath::Abs(fPrimaryVertex->GetZ()) > 10.0))
    return kFALSE;

  FillHistogram("hCentrality", centralityPercentile, 3.5); // after vertex position cut

  return kTRUE;
}

//________________________________________________________________________
inline Bool_t AliAnalysisTaskChargedJetsPA::IsTrackInAcceptance(AliVParticle* track)
{
  FillHistogram("hTrackAcceptance", 0.5);
  if (track != 0)
  {
    if ((track->Eta() < fMaxEta) && (track->Eta() >= fMinEta))
    {
      FillHistogram("hTrackAcceptance", 1.5);
      if (track->Pt() >= fMinTrackPt)
      {
        FillHistogram("hTrackAcceptance", 2.5);
        return kTRUE;
      }
    }
  }
  return kFALSE;
}

//________________________________________________________________________
inline Bool_t AliAnalysisTaskChargedJetsPA::IsBackgroundJetInAcceptance(AliEmcalJet *jet)
{   
  if (jet != 0)
    if ((jet->Eta() >= fMinJetEta) && (jet->Eta() < fMaxJetEta))
      if (jet->Pt() >= fMinBackgroundJetPt)
        return kTRUE;

  return kFALSE;
}

//________________________________________________________________________
inline Bool_t AliAnalysisTaskChargedJetsPA::IsSignalJetInAcceptance(AliEmcalJet *jet, Bool_t usePtCut)
{
  Bool_t acceptedWithPtCut = kFALSE;
  Bool_t acceptedWithoutPtCut = kFALSE;

  FillHistogram("hJetAcceptance", 0.5);
  if (jet != 0)
    if ((jet->Eta() >= fMinJetEta) && (jet->Eta() < fMaxJetEta))
    {
      FillHistogram("hJetAcceptance", 1.5);
      if (jet->Pt() >= fMinJetPt) // jet fulfills pt cut
      {
        FillHistogram("hJetAcceptance", 2.5);
        if (jet->Area() >= fMinJetArea)
        {
          FillHistogram("hJetAcceptance", 3.5);
          acceptedWithPtCut = kTRUE;
        }
      }
      else if(!usePtCut) // jet does not fulfill pt cut
      {
        if (jet->Area() >= fMinJetArea)
          acceptedWithoutPtCut = kTRUE;
      }
    }

  if(usePtCut)
    return (acceptedWithPtCut);
  else
    return (acceptedWithPtCut || acceptedWithoutPtCut);
}

//________________________________________________________________________
void AliAnalysisTaskChargedJetsPA::ExecOnce()
{
  #ifdef DEBUGMODE
    AliInfo("Starting ExecOnce.");
  #endif
  fInitialized = kTRUE;

  // Check for track array
  if (strcmp(fTrackArrayName.Data(), "") != 0)
  {
    fTrackArray = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(fTrackArrayName.Data()));
    if (!fTrackArray) 
      AliWarning(Form("%s: Could not retrieve tracks %s!", GetName(), fTrackArrayName.Data())); 
    else
    {
      TClass *cl = fTrackArray->GetClass();
      if (!cl->GetBaseClass("AliVParticle"))
      {
        AliError(Form("%s: Collection %s does not contain AliVParticle objects!", GetName(), fTrackArrayName.Data())); 
        fTrackArray = 0;
      }
    }
  }

  // Check for jet array
  if (strcmp(fJetArrayName.Data(), "") != 0)
  {
    fJetArray = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(fJetArrayName.Data()));
    if (!fJetArray) 
      AliWarning(Form("%s: Could not retrieve jets %s!", GetName(), fJetArrayName.Data())); 
    else
    {
      if (!fJetArray->GetClass()->GetBaseClass("AliEmcalJet")) 
      {
        AliError(Form("%s: Collection %s does not contain AliEmcalJet objects!", GetName(), fJetArrayName.Data())); 
        fJetArray = 0;
      }
    }
  }

  // Check for background object
  if (strcmp(fBackgroundJetArrayName.Data(), "") != 0)
  {
    fBackgroundJetArray = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(fBackgroundJetArrayName.Data()));
    if (!fBackgroundJetArray)
      AliInfo(Form("%s: Could not retrieve background jets %s!", GetName(), fBackgroundJetArrayName.Data())); 
  }

  // Initialize helper class (for vertex selection & pile up correction)
  fHelperClass = new AliAnalysisUtils();
  fHelperClass->SetCutOnZVertexSPD(kFALSE);
  // Histogram init
  Init();
  // Trackcut initialization
  InitializeTrackcuts();

  #ifdef DEBUGMODE
    AliInfo("ExecOnce done.");
  #endif

}

//________________________________________________________________________
void AliAnalysisTaskChargedJetsPA::GetLeadingJets()
{
  // Reset vars
  fFirstLeadingJet = NULL;
  fSecondLeadingJet = NULL;
  fFirstLeadingKTJet = NULL;
  fSecondLeadingKTJet = NULL;

  fNumberSignalJets = 0;
  fNumberSignalJetsAbove5GeV = 0;

  Int_t jetIDArray[]   = {-1, -1};
  Float_t maxJetPts[] = {0, 0};
  jetIDArray[0] = -1;
  jetIDArray[1] = -1;

  Int_t jetIDArrayKT[]   = {-1, -1};
  Float_t maxJetPtsKT[] = {0, 0};
  jetIDArrayKT[0] = -1;
  jetIDArrayKT[1] = -1;

  // Find leading signal jets
  for (Int_t i = 0; i < fJetArray->GetEntries(); i++)
  {
    AliEmcalJet* jet = static_cast<AliEmcalJet*>(fJetArray->At(i));
    if (!jet) 
    {
      AliError(Form("%s: Could not receive jet %d", GetName(), i));
      continue;
    }

    if (!IsSignalJetInAcceptance(jet)) continue;

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
    fNumberSignalJets++;
    if(jet->Pt() >= 5.)
      fNumberSignalJetsAbove5GeV++;
  }

  // Find leading background jets
  for (Int_t i = 0; i < fBackgroundJetArray->GetEntries(); i++)
  {
    AliEmcalJet* jet = static_cast<AliEmcalJet*>(fBackgroundJetArray->At(i));
    if (!jet) 
    {
      AliError(Form("%s: Could not receive jet %d", GetName(), i));
      continue;
    }

    if (!IsBackgroundJetInAcceptance(jet)) continue;

    if (jet->Pt() > maxJetPtsKT[0]) 
    {
      maxJetPtsKT[1] = maxJetPtsKT[0];
      jetIDArrayKT[1] = jetIDArrayKT[0];
      maxJetPtsKT[0] = jet->Pt();
      jetIDArrayKT[0] = i;
    }
    else if (jet->Pt() > maxJetPtsKT[1]) 
    {
      maxJetPtsKT[1] = jet->Pt();
      jetIDArrayKT[1] = i;
    }
  }

  if (jetIDArray[0] > -1)
    fFirstLeadingJet  = static_cast<AliEmcalJet*>(fJetArray->At(jetIDArray[0]));
  if (jetIDArray[1] > -1)
    fSecondLeadingJet = static_cast<AliEmcalJet*>(fJetArray->At(jetIDArray[1]));
  if (jetIDArrayKT[0] > -1)
    fFirstLeadingKTJet  = static_cast<AliEmcalJet*>(fBackgroundJetArray->At(jetIDArrayKT[0]));
  if (jetIDArrayKT[1] > -1)
    fSecondLeadingKTJet = static_cast<AliEmcalJet*>(fBackgroundJetArray->At(jetIDArrayKT[1]));
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
inline Double_t AliAnalysisTaskChargedJetsPA::GetCorrectedConePt(Double_t eta, Double_t phi, Double_t radius, Double_t background)
{
  Double_t tmpConePt = 0.0;

  for (Int_t i = 0; i < fTrackArray->GetEntries(); i++)
  {
    AliVTrack* tmpTrack = static_cast<AliVTrack*>(fTrackArray->At(i));
    if (IsTrackInAcceptance(tmpTrack))
      if(IsTrackInCone(tmpTrack, eta, phi, radius))
        tmpConePt = tmpConePt + tmpTrack->Pt();
  }
  Double_t realConeArea = (1.0*(fMaxEta-fMinEta)) * TMath::TwoPi() * MCGetOverlapCircleRectancle(eta, phi, radius, fMinEta, fMaxEta, 0., TMath::TwoPi());
  tmpConePt -= background * realConeArea; // subtract background

  return tmpConePt;
}

//________________________________________________________________________
inline Int_t AliAnalysisTaskChargedJetsPA::GetConeConstituentCount(Double_t eta, Double_t phi, Double_t radius)
{
  Int_t tmpConeCount = 0.0;

  for (Int_t i = 0; i < fTrackArray->GetEntries(); i++)
  {
    AliVTrack* tmpTrack = static_cast<AliVTrack*>(fTrackArray->At(i));
    if (IsTrackInAcceptance(tmpTrack))
      if(IsTrackInCone(tmpTrack, eta, phi, radius))
        tmpConeCount++;
  }
 
  return tmpConeCount;
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
Double_t AliAnalysisTaskChargedJetsPA::GetCorrectedJetPt(AliEmcalJet* jet, Double_t background)
{
  #ifdef DEBUGMODE
    AliInfo("Getting corrected jet spectra.");
  #endif

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
Double_t AliAnalysisTaskChargedJetsPA::GetDeltaPt(Double_t rho, Double_t overlappingJetExclusionProbability)
{
  #ifdef DEBUGMODE
    AliInfo("Getting Delta Pt.");
  #endif

  // Define an invalid delta pt
  Double_t deltaPt = -10000.0;

  // Define the ratio of excluded RCs over all RCs
  Double_t ratioExcludedRCs = static_cast<Double_t>(fTempExcludedRCs)/fTempAllRCs;

  // Define eta range
  Double_t etaMin, etaMax;
  etaMin = fMinEta+fRandConeRadius;
  etaMax = fMaxEta-fRandConeRadius;

  // Define random cone
  Bool_t coneValid = kTRUE;
  Double_t tmpRandConeEta = etaMin + fRandom->Rndm()*(etaMax-etaMin);
  Double_t tmpRandConePhi = fRandom->Rndm()*TMath::TwoPi();

  // Check if a signal jet is overlapping with random cone
  if(overlappingJetExclusionProbability)
  {
    // Calculate the mean exclusion probability
    fTempOverlapCounter++;
    fTempMeanExclusionProbability += overlappingJetExclusionProbability;
    // For all jets, check overlap
    for (Int_t i = 0; i<fJetArray->GetEntries(); i++)
    {
      AliEmcalJet* tmpJet = static_cast<AliEmcalJet*>(fJetArray->At(i));
      if (!IsSignalJetInAcceptance(tmpJet, kTRUE)) continue;

      // Check overlap
      Double_t tmpDeltaPhi = GetDeltaPhi(tmpRandConePhi, tmpJet->Phi());
      if ( tmpDeltaPhi*tmpDeltaPhi + (tmpRandConeEta-tmpJet->Eta())*(tmpRandConeEta-tmpJet->Eta()) <= fRandConeRadius*fRandConeRadius )
      {
        fTempAllRCs++;
        // If an overlap is given, discard or accept it according to the exclusion prob. 
        if(ratioExcludedRCs < fTempMeanExclusionProbability/fTempOverlapCounter) // to few RCs excluded -> exclude this one
        {
          coneValid = kFALSE;
          fTempExcludedRCs++;
        }
        else  // to many RCs excluded -> take this one
          coneValid = kTRUE;
      }
    }
  }


  // Get the cones' pt and calculate delta pt
  if (coneValid)
    deltaPt = GetConePt(tmpRandConeEta,tmpRandConePhi,fRandConeRadius) - (rho*fRandConeRadius*fRandConeRadius*TMath::Pi());

  #ifdef DEBUGMODE
    AliInfo("Got Delta Pt.");
  #endif
  return deltaPt;
}

//________________________________________________________________________
void AliAnalysisTaskChargedJetsPA::GetKTBackgroundDensityAll(Int_t numberExcludeLeadingJets, Double_t& rhoPbPb, Double_t& rhoPbPbWithGhosts, Double_t& rhoCMS, Double_t& rhoImprovedCMS, Double_t& rhoMean, Double_t& rhoTrackLike)
{
  #ifdef DEBUGMODE
    AliInfo("Getting ALL KT background density.");
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

    tmpSummedArea += backgroundJet->Area();
    if(backgroundJet->Pt() > 0.150)
      tmpCoveredArea += backgroundJet->Area();

    if (!IsBackgroundJetInAcceptance(backgroundJet))
      continue;

    // Search for overlap with signal jets
    Bool_t isOverlapping = kFALSE;
    for(Int_t j=0;j<numberExcludeLeadingJets;j++)
    {
      AliEmcalJet* signalJet = fFirstLeadingJet;
      if(j==1)
        signalJet = fSecondLeadingJet;

      if(signalJet->Pt() < 5.0)
          continue;

      if(IsJetOverlapping(signalJet, backgroundJet))
      {
        isOverlapping = kTRUE;
        break;
      }
    }

    Double_t tmpRho = 0.0;
    if(backgroundJet->Area())
      tmpRho = backgroundJet->Pt() / backgroundJet->Area();

    // PbPb approach (take ghosts into account)
    if((backgroundJet != fFirstLeadingKTJet) || (backgroundJet != fSecondLeadingKTJet))
    {
      tmpRhoPbPbWithGhosts[rhoPbPbWithGhostsJetCount] = tmpRho;
      rhoPbPbWithGhostsJetCount++;
    }

    if(backgroundJet->Pt() > 0.150)
    {
      // CMS approach: don't take ghosts into acount
      tmpRhoCMS[rhoCMSJetCount] = tmpRho;
      rhoCMSJetCount++;

      // Improved CMS approach: like CMS but excluding signal
      if((backgroundJet != fFirstLeadingKTJet) || (backgroundJet != fSecondLeadingKTJet))
      {
        tmpRhoImprovedCMS[rhoImprovedCMSJetCount] = tmpRho;
        rhoImprovedCMSJetCount++;
      }

      // PbPb w/o ghosts approach (just neglect ghosts)
      if((backgroundJet != fFirstLeadingKTJet) || (backgroundJet != fSecondLeadingKTJet))
      {  
        tmpRhoPbPb[rhoPbPbJetCount] = tmpRho;
        rhoPbPbJetCount++;
      }
    }

    // (no overlap with signal jets)
    if(!isOverlapping)
    {
      // Mean approach
      tmpRhoMean[rhoMeanJetCount] = tmpRho;
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
  {
    rhoImprovedCMS = TMath::Median(rhoImprovedCMSJetCount, tmpRhoImprovedCMS) * tmpCoveredArea/tmpSummedArea;
  }
  if (rhoMeanJetCount > 0)
    rhoMean = TMath::Mean(rhoMeanJetCount, tmpRhoMean);

  #ifdef DEBUGMODE
    AliInfo("Got ALL KT background density.");
  #endif
}

//________________________________________________________________________
void AliAnalysisTaskChargedJetsPA::GetTRBackgroundDensity(Int_t numberExcludeLeadingJets, Double_t& rhoNoExclusion, Double_t& rhoConeExclusion02, Double_t& rhoConeExclusion04, Double_t& rhoConeExclusion06, Double_t& rhoConeExclusion08, Double_t& rhoExactExclusion)
{
  #ifdef DEBUGMODE
    AliInfo("Getting TR background density.");
  #endif

  Double_t summedTracksPtCone04 = 0.0;
  Double_t summedTracksPtCone02 = 0.0;
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
    numberExcludeLeadingJets = fNumberSignalJetsAbove5GeV;
  if (fNumberSignalJets < numberExcludeLeadingJets)
    numberExcludeLeadingJets = fNumberSignalJetsAbove5GeV;
  if(numberExcludeLeadingJets>2)
  {
    AliWarning(Form("Warning: GetTRBackgroundDensity() can only exclude up to 2 leading jets! Demanded %i", numberExcludeLeadingJets) );
    numberExcludeLeadingJets = 2;
  }

  for (Int_t i = 0; i < fTrackArray->GetEntries(); i++)
  {
    AliVTrack* tmpTrack = static_cast<AliVTrack*>(fTrackArray->At(i));
    Bool_t trackWithinJet = kFALSE; Bool_t trackWithin02Cone = kFALSE; Bool_t trackWithin04Cone = kFALSE; Bool_t trackWithin06Cone = kFALSE; Bool_t trackWithin08Cone = kFALSE;

    if (IsTrackInAcceptance(tmpTrack))
    {
      // Check if tracks overlaps with jet
      for(Int_t j=0;j<numberExcludeLeadingJets;j++)
      {
        AliEmcalJet* signalJet = fFirstLeadingJet;
        if(j==1)
          signalJet = fSecondLeadingJet;

        if(signalJet->Pt() < 5.0)
          continue;

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

  Double_t tmpFullTPCArea = (1.0*(fMaxEta-fMinEta)) * TMath::TwoPi();
  Double_t tmpAreaCone02     = tmpFullTPCArea;
  Double_t tmpAreaCone04     = tmpFullTPCArea;
  Double_t tmpAreaCone06     = tmpFullTPCArea;
  Double_t tmpAreaCone08     = tmpFullTPCArea;
  Double_t tmpAreaWithinJets = tmpFullTPCArea;
  std::vector<Double_t> tmpEtas(fNumberSignalJetsAbove5GeV);
  std::vector<Double_t> tmpPhis(fNumberSignalJetsAbove5GeV);

  Int_t iSignal = 0;
  for(Int_t i=0;i<numberExcludeLeadingJets;i++)
  {
    AliEmcalJet* signalJet = fFirstLeadingJet;
    if(i==1)
      signalJet = fSecondLeadingJet;

    if(signalJet->Pt() < 5.0)
      continue;

    tmpEtas[iSignal] = signalJet->Eta();
    tmpPhis[iSignal] = signalJet->Phi();
    tmpAreaWithinJets -= signalJet->Area();

    iSignal++;
  }

  tmpAreaCone02 -= tmpFullTPCArea * MCGetOverlapMultipleCirclesRectancle(fNumberSignalJetsAbove5GeV, tmpEtas, tmpPhis, 0.2, fMinEta, fMaxEta, 0., TMath::TwoPi());
  tmpAreaCone04 -= tmpFullTPCArea * MCGetOverlapMultipleCirclesRectancle(fNumberSignalJetsAbove5GeV, tmpEtas, tmpPhis, 0.4, fMinEta, fMaxEta, 0., TMath::TwoPi());
  tmpAreaCone06 -= tmpFullTPCArea * MCGetOverlapMultipleCirclesRectancle(fNumberSignalJetsAbove5GeV, tmpEtas, tmpPhis, 0.6, fMinEta, fMaxEta, 0., TMath::TwoPi());
  tmpAreaCone08 -= tmpFullTPCArea * MCGetOverlapMultipleCirclesRectancle(fNumberSignalJetsAbove5GeV, tmpEtas, tmpPhis, 0.8, fMinEta, fMaxEta, 0., TMath::TwoPi());
 
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
void AliAnalysisTaskChargedJetsPA::GetPPBackgroundDensity(Double_t& background)
{
  // This is the background that was used for the pp 7 TeV ALICE paper
  // The background is estimated using the leading jet

  background = 0;

  AliEmcalJet* jet = NULL;
  if(fFirstLeadingJet)
    jet = fFirstLeadingJet;
  else
    return;

  Double_t jetMom[3] = { jet->Px(), jet->Py(), jet->Pz() };
  TVector3 jet3mom1(jetMom);
  TVector3 jet3mom2(jetMom);

  jet3mom1.RotateZ(TMath::Pi());
  jet3mom2.RotateZ(-TMath::Pi());

  for (int i = 0; i < fTrackArray->GetEntries(); i++)
  {
    AliVTrack* track = static_cast<AliVTrack*>(fTrackArray->At(i));
    if (!IsTrackInAcceptance(track))
      continue;

    Double_t trackMom[3] = { track->Px(), track->Py(), track->Pz() };
    TVector3 track3mom(trackMom);

    Double_t dR1 = jet3mom1.DeltaR(track3mom);
    Double_t dR2 = jet3mom2.DeltaR(track3mom);

    if (dR1 <= fSignalJetRadius || dR2 <= fSignalJetRadius)
      background += track3mom.Pt();
  }

  background /= (2 * TMath::Pi() * fSignalJetRadius * fSignalJetRadius);
}

//________________________________________________________________________
void AliAnalysisTaskChargedJetsPA::Calculate(AliVEvent* event)
{
  #ifdef DEBUGMODE
    AliInfo("Starting Calculate().");
  #endif
  ////////////////////// NOTE: initialization & casting

  fEventCounter++;

  // This is to take only every Nth event
  if((fEventCounter+fPartialAnalysisIndex) % fPartialAnalysisNParts != 0)
    return;

  if(!IsEventInAcceptance(event))
    return;

  #ifdef DEBUGMODE
    AliInfo("Calculate()::Init done.");
  #endif

  ////////////////////// NOTE: Create cut histograms

  if(fAnalyzeTrackcuts)
    CreateCutHistograms();

  ////////////////////// NOTE: Get Centrality, (Leading)Signal jets and Background

  // Get centrality
  AliCentrality* tmpCentrality = event->GetCentrality();
  Double_t centralityPercentile = -1.0;
  Double_t centralityPercentileCL1 = 0.0;
  Double_t centralityPercentileV0A = 0.0;
  Double_t centralityPercentileV0C = 0.0;
  Double_t centralityPercentileV0M = 0.0;
  Double_t centralityPercentileZNA = 0.0;
  if (tmpCentrality != NULL)
  {
    centralityPercentile    = tmpCentrality->GetCentralityPercentile(fCentralityType.Data());
    centralityPercentileCL1 = tmpCentrality->GetCentralityPercentile("CL1");
    centralityPercentileV0A = tmpCentrality->GetCentralityPercentile("V0A");
    centralityPercentileV0C = tmpCentrality->GetCentralityPercentile("V0C");
    centralityPercentileV0M = tmpCentrality->GetCentralityPercentile("V0M");
    centralityPercentileZNA = tmpCentrality->GetCentralityPercentile("ZNA");
  }

  if(fSetCentralityToOne)
    centralityPercentile = 1.0;

  ////////////////////// NOTE: Get event QA histograms

  FillHistogram("hVertexX",fPrimaryVertex->GetX());
  FillHistogram("hVertexY",fPrimaryVertex->GetY());
  FillHistogram("hVertexXY",fPrimaryVertex->GetX(), fPrimaryVertex->GetY());
  FillHistogram("hVertexR",TMath::Sqrt(fPrimaryVertex->GetX()*fPrimaryVertex->GetX() + fPrimaryVertex->GetY()*fPrimaryVertex->GetY()));
  FillHistogram("hCentralityCL1",centralityPercentileCL1);
  FillHistogram("hCentralityV0M",centralityPercentileV0M);
  FillHistogram("hCentralityV0A",centralityPercentileV0A);
  FillHistogram("hCentralityV0C",centralityPercentileV0C);
  FillHistogram("hCentralityZNA",centralityPercentileZNA);

  if(!fDoJetAnalysis)
    return;

  GetLeadingJets();


/*
  //DEBUG
  if(fFirstLeadingJet->Pt()>=80.)
  {
    const char* fname = CurrentFileName();
    TObjString* tmpStr = new TObjString(Form("jet pT=%3.3f, fname=%s, entry=%i", fFirstLeadingJet->Pt(), fname, AliAnalysisManager::GetAnalysisManager()->GetCurrentEntry()));
    fCurrentOutputList->Add(tmpStr);
  }
  //DEBUG
*/

  // ##################### Calculate background densities
  Double_t              backgroundKTImprovedCMS = -1.0;
  Double_t              backgroundExternal = -1.0;
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
  Double_t              backgroundPP       = -1.0;
  Double_t              backgroundJetProfile = -1.0;

  // Get background estimates
  GetKTBackgroundDensityAll (fNumberExcludedJets, backgroundKTPbPb, backgroundKTPbPbWithGhosts, backgroundKTCMS, backgroundKTImprovedCMS, backgroundKTMean, backgroundKTTrackLike);
  GetTRBackgroundDensity    (fNumberExcludedJets, backgroundTRNoExcl, backgroundTRCone02, backgroundTRCone04, backgroundTRCone06, backgroundTRCone08, backgroundTRExact);
  GetPPBackgroundDensity(backgroundPP);

  backgroundExternal = GetExternalRho();
  if(fNoExternalBackground)
    backgroundExternal = 0;

  if(fBackgroundForJetProfile==0)
    backgroundJetProfile = backgroundExternal;
  else if(fBackgroundForJetProfile==1)
    backgroundJetProfile = backgroundKTImprovedCMS;
  else if(fBackgroundForJetProfile==2)
    backgroundJetProfile = backgroundKTCMS;
  else if(fBackgroundForJetProfile==3)
    backgroundJetProfile = backgroundPP;
  else if(fBackgroundForJetProfile==4)
    backgroundJetProfile = backgroundTRCone06;
  else if(fBackgroundForJetProfile==5)
    backgroundJetProfile = 0;

  #ifdef DEBUGMODE
    AliInfo("Calculate()::Centrality&SignalJets&Background-Calculation done.");
  #endif

  // ##################### Fill event QA histograms

  Int_t trackCountAcc = 0;
  Int_t nTracks = fTrackArray->GetEntries();
  for (Int_t i = 0; i < nTracks; i++)
  {
    AliVTrack* track = static_cast<AliVTrack*>(fTrackArray->At(i));

    if (track != 0)
      if (track->Pt() >= fMinTrackPt)
      {
        FillHistogram("hTrackPhiEta", track->Phi(),track->Eta(), 1);
        FillHistogram("hTrackPtPhiEta", track->Phi(),track->Eta(), track->Pt());
      }

    if (IsTrackInAcceptance(track))
    {
      FillHistogram("hTrackPt", track->Pt(), centralityPercentile);

      if(track->Eta() >= 0)
        FillHistogram("hTrackPtPosEta", track->Pt(), centralityPercentile);
      else
        FillHistogram("hTrackPtNegEta", track->Pt(), centralityPercentile);
              
      FillHistogram("hTrackEta", track->Eta(), centralityPercentile);
      FillHistogram("hTrackPhi", track->Phi());
      
      if(static_cast<AliPicoTrack*>(track))
      {
        FillHistogram("hTrackPhiTrackType", track->Phi(), (static_cast<AliPicoTrack*>(track))->GetTrackType());
        FillHistogram("hTrackPtTrackType", track->Pt(), (static_cast<AliPicoTrack*>(track))->GetTrackType());
      }

      for(Int_t j=0;j<20;j++)
        if(track->Pt() > j)
          FillHistogram("hTrackPhiPtCut", track->Phi(), track->Pt());

      FillHistogram("hTrackCharge", track->Charge());
      trackCountAcc++;
    }
  }
  FillHistogram("hTrackCountAcc", trackCountAcc, centralityPercentile);

  #ifdef DEBUGMODE
    AliInfo("Calculate()::QA done.");
  #endif

  // ##################### Fill jet histograms

  FillHistogram("hJetCountAll", fJetArray->GetEntries());
  FillHistogram("hJetCountAccepted", fNumberSignalJets);
  FillHistogram("hJetCount", fJetArray->GetEntries(), fNumberSignalJets);
  if (fFirstLeadingJet)
  {
    FillHistogram("hLeadingJetPt", fFirstLeadingJet->Pt());
    FillHistogram("hCorrectedLeadingJetPt", GetCorrectedJetPt(fFirstLeadingJet,backgroundExternal));
  }
  if (fSecondLeadingJet)
  {
    FillHistogram("hSecondLeadingJetPt", fSecondLeadingJet->Pt());
    FillHistogram("hCorrectedSecondLeadingJetPt", GetCorrectedJetPt(fSecondLeadingJet,backgroundExternal));
  }

  for (Int_t i = 0; i<fJetArray->GetEntries(); i++)
  {
    AliEmcalJet* tmpJet = static_cast<AliEmcalJet*>(fJetArray->At(i));
    if (!tmpJet)
      continue;

    // ### JETS BEFORE ANY CUTS
    if (tmpJet->Area() >= fMinJetArea)
      FillHistogram("hRawJetPhiEta", tmpJet->Phi(), tmpJet->Eta());
    if ((tmpJet->Eta() >= fMinJetEta) && (tmpJet->Eta() < fMaxJetEta))
      FillHistogram("hRawJetArea", tmpJet->Area());

    FillHistogram("hJetPtCutStages", tmpJet->Pt(), 0.5);
    if ((tmpJet->Eta() >= fMinJetEta) && (tmpJet->Eta() < fMaxJetEta))
    {
      FillHistogram("hJetPtCutStages", tmpJet->Pt(), 1.5);
      if (tmpJet->Pt() >= fMinJetPt)
      {
        FillHistogram("hJetPtCutStages", tmpJet->Pt(), 2.5);
        if (tmpJet->Area() >= fMinJetArea)
        {
          FillHistogram("hJetPtCutStages", tmpJet->Pt(), 3.5);
        }
      }
    }

    // ### JETS AFTER CUTS
    if(IsSignalJetInAcceptance(tmpJet))
    {
      // Background corrected jet spectra
      FillHistogram("hJetPtNoBgrdSubtracted", GetCorrectedJetPt(tmpJet, 0.0), centralityPercentile);
      FillHistogram("hJetPtBgrdSubtractedExternal", GetCorrectedJetPt(tmpJet, backgroundExternal), centralityPercentile);
      FillHistogram("hJetPtBgrdSubtractedKTImprovedCMS", GetCorrectedJetPt(tmpJet, backgroundKTImprovedCMS), centralityPercentile);
      FillHistogram("hJetPtBgrdSubtractedPP", GetCorrectedJetPt(tmpJet, backgroundPP), centralityPercentile);
      if(tmpJet->Phi() >= TMath::Pi())
        FillHistogram("hJetPtBgrdSubtractedExternal_Phi2", GetCorrectedJetPt(tmpJet, backgroundExternal), centralityPercentile);
      else          
        FillHistogram("hJetPtBgrdSubtractedExternal_Phi1", GetCorrectedJetPt(tmpJet, backgroundExternal), centralityPercentile);
      FillHistogram("hJetPtBgrdSubtractedTR", GetCorrectedJetPt(tmpJet, backgroundTRCone06), centralityPercentile);
      FillHistogram("hJetPtBgrdSubtractedKTPbPb", GetCorrectedJetPt(tmpJet, backgroundKTPbPb), centralityPercentile);
      FillHistogram("hJetPtBgrdSubtractedKTPbPbWithGhosts", GetCorrectedJetPt(tmpJet, backgroundKTPbPbWithGhosts), centralityPercentile);
      FillHistogram("hJetPtBgrdSubtractedKTCMS", GetCorrectedJetPt(tmpJet, backgroundKTCMS), centralityPercentile);
      FillHistogram("hJetPtBgrdSubtractedKTMean", GetCorrectedJetPt(tmpJet, backgroundKTMean), centralityPercentile);
      FillHistogram("hJetPtBgrdSubtractedKTTrackLike", GetCorrectedJetPt(tmpJet, backgroundKTTrackLike), centralityPercentile);

      FillHistogram("hJetPtSubtractedRhoExternal", tmpJet->Pt(), centralityPercentile, tmpJet->Pt() - GetCorrectedJetPt(tmpJet, backgroundExternal));
      FillHistogram("hJetPtSubtractedRhoKTImprovedCMS", tmpJet->Pt(), centralityPercentile, tmpJet->Pt() - GetCorrectedJetPt(tmpJet, backgroundKTImprovedCMS));
      FillHistogram("hJetPtSubtractedRhoPP", tmpJet->Pt(), centralityPercentile, tmpJet->Pt() - GetCorrectedJetPt(tmpJet, backgroundPP));
      for(Int_t j=0; j<fRandConeNumber; j++)
        FillHistogram("hDeltaPtExternalBgrdVsPt", GetDeltaPt(backgroundExternal), GetCorrectedJetPt(tmpJet, backgroundExternal));
      FillHistogram("hKTBackgroundExternalVsPt", backgroundExternal, GetCorrectedJetPt(tmpJet, backgroundExternal));

      // ###### CONSTITUENT ANALYSIS

      if(fAnalyzeJetConstituents)
      {
        THnF* tmpConstituentHist = static_cast<THnF*>(fCurrentOutputList->FindObject("hJetConstituents"));
        THnF* tmpConstituentDistanceHist = static_cast<THnF*>(fCurrentOutputList->FindObject("hJetConstituentDistance"));

        for(Int_t j=0; j<tmpJet->GetNumberOfTracks(); j++)
        {
          AliVParticle* tmpTrack = tmpJet->TrackAt(j, fTrackArray);
          // Define random cone  
          Double_t tmpRandConeEta = fMinJetEta + fRandom->Rndm()*TMath::Abs(fMaxJetEta-fMinJetEta);
          Double_t tmpRandConePhi = fRandom->Rndm()*TMath::TwoPi();

          Double_t tmpPConeEta = tmpJet->Eta();
          Double_t tmpPConePhi = tmpJet->Phi() + TMath::Pi();

          if(tmpPConePhi>=TMath::TwoPi())
            tmpPConePhi = tmpPConePhi - TMath::TwoPi();

          Double_t tmpRCcount  = GetConeConstituentCount(tmpRandConeEta, tmpRandConePhi, fSignalJetRadius);
          Double_t tmpPCcount  = GetConeConstituentCount(tmpPConeEta, tmpPConePhi, fSignalJetRadius);

          Double_t tmpDistance = TMath::Sqrt( (tmpJet->Eta()-tmpTrack->Eta())*(tmpJet->Eta()-tmpTrack->Eta()) 
                                            + (tmpJet->Phi()-tmpTrack->Phi())*(tmpJet->Phi()-tmpTrack->Phi()) ); // distance between jet axis and track

          Double_t tmpVec1[5] = {tmpJet->Pt(), tmpTrack->Pt(), static_cast<Double_t>(tmpJet->GetNumberOfTracks()), tmpRCcount, tmpPCcount};
          Double_t tmpVec2[4] = {tmpJet->Pt(), tmpTrack->Pt(), static_cast<Double_t>(tmpJet->GetNumberOfTracks()), tmpDistance};


          tmpConstituentHist->Fill(tmpVec1);
          tmpConstituentDistanceHist->Fill(tmpVec2);

          FillHistogram("hJetConstituentPtVsJetPt", tmpTrack->Pt(), tmpJet->Pt());
        }
      }
      
      FillHistogram("hJetPtVsConstituentCount", tmpJet->Pt(),tmpJet->GetNumberOfTracks());

      // Leading track biased jets
      Double_t leadingTrackPt = 0.0;
      for(Int_t j=0; j<tmpJet->GetNumberOfTracks(); j++)
      {
        if(tmpJet->TrackAt(j, fTrackArray)->Pt() > leadingTrackPt)
          leadingTrackPt = tmpJet->TrackAt(j, fTrackArray)->Pt();
      }

      if(leadingTrackPt >= 10)
        FillHistogram("hJetPtBgrdSubtractedKTImprovedCMS_Biased_10GeV", GetCorrectedJetPt(tmpJet, backgroundKTImprovedCMS), centralityPercentile);
      else if(leadingTrackPt >= 5)
        FillHistogram("hJetPtBgrdSubtractedKTImprovedCMS_Biased_5GeV", GetCorrectedJetPt(tmpJet, backgroundKTImprovedCMS), centralityPercentile);
      else if(leadingTrackPt >= 2)
        FillHistogram("hJetPtBgrdSubtractedKTImprovedCMS_Biased_2GeV", GetCorrectedJetPt(tmpJet, backgroundKTImprovedCMS), centralityPercentile);


      // Fill jet constituent histograms
      for(Int_t j=0; j<tmpJet->GetNumberOfTracks(); j++)
      {
        FillHistogram("hJetConstituentPt0GeV", tmpJet->TrackAt(j, fTrackArray)->Pt(), centralityPercentile);
        if(tmpJet->Pt() >= 1.0)
          FillHistogram("hJetConstituentPt1GeV", tmpJet->TrackAt(j, fTrackArray)->Pt(), centralityPercentile);
        if(tmpJet->Pt() >= 2.0)
          FillHistogram("hJetConstituentPt2GeV", tmpJet->TrackAt(j, fTrackArray)->Pt(), centralityPercentile);
        if(tmpJet->Pt() >= 3.0)
          FillHistogram("hJetConstituentPt3GeV", tmpJet->TrackAt(j, fTrackArray)->Pt(), centralityPercentile);
        if(tmpJet->Pt() >= 4.0)
          FillHistogram("hJetConstituentPt4GeV", tmpJet->TrackAt(j, fTrackArray)->Pt(), centralityPercentile);
        if(tmpJet->Pt() >= 5.0)
          FillHistogram("hJetConstituentPt5GeV", tmpJet->TrackAt(j, fTrackArray)->Pt(), centralityPercentile);
        if(tmpJet->Pt() >= 7.0)
          FillHistogram("hJetConstituentPt7GeV", tmpJet->TrackAt(j, fTrackArray)->Pt(), centralityPercentile);
        if(tmpJet->Pt() >= 10.0)
          FillHistogram("hJetConstituentPt10GeV", tmpJet->TrackAt(j, fTrackArray)->Pt(), centralityPercentile);
      }

      if(tmpJet->Pt() >= 5.0)
      {
        Double_t lowestTrackPt = 1e99;
        Double_t highestTrackPt = 0.0;
        for(Int_t j=0; j<tmpJet->GetNumberOfTracks(); j++)
        {
//          FillHistogram("hJetConstituentPt", tmpJet->TrackAt(j, fTrackArray)->Pt(), centralityPercentile);
          // Find the lowest pT of a track in the jet
          if (tmpJet->TrackAt(j, fTrackArray)->Pt() < lowestTrackPt)
            lowestTrackPt = tmpJet->TrackAt(j, fTrackArray)->Pt();
          if (tmpJet->TrackAt(j, fTrackArray)->Pt() > highestTrackPt)
            highestTrackPt = tmpJet->TrackAt(j, fTrackArray)->Pt();
        }
        FillHistogram("hJetArea", tmpJet->Area(), tmpJet->Pt());
        // Signal jet vs. signal jet - "Combinatorial"
        for (Int_t j = 0; j<fJetArray->GetEntries(); j++)
        {
          AliEmcalJet* tmpJet2 = static_cast<AliEmcalJet*>(fJetArray->At(j));
          if (!tmpJet2)
            continue;
          if(tmpJet2->Pt() >= 5.0)
            FillHistogram("hJetDeltaPhi", GetDeltaPhi(tmpJet->Phi(), tmpJet2->Phi()));
        }

        FillHistogram("hJetPhiEta", tmpJet->Phi(),tmpJet->Eta());
        FillHistogram("hJetPtPhiEta", tmpJet->Phi(),tmpJet->Eta(),tmpJet->Pt());
        FillHistogram("hJetEta", tmpJet->Eta(), centralityPercentile);

        if(lowestTrackPt>=2.0)
          FillHistogram("hJetEta2GeVTracks", tmpJet->Eta(), centralityPercentile);
        if(lowestTrackPt>=4.0)
          FillHistogram("hJetEta4GeVTracks", tmpJet->Eta(), centralityPercentile);
      }
    }
  } // end of jet loop

  if(fAnalyzeJetProfile)
    CreateJetProfilePlots(backgroundJetProfile);

  #ifdef DEBUGMODE
    AliInfo("Calculate()::Jets done.");
  #endif

  // ##################### Fill background plots

  FillHistogram("hKTBackgroundExternal", backgroundExternal, centralityPercentile);
  if(fFirstLeadingJet && (fFirstLeadingJet->Pt()>=20.))
    FillHistogram("hKTBackgroundExternal20GeV", backgroundExternal, centralityPercentile);

  FillHistogram("hKTBackgroundImprovedCMS", backgroundKTImprovedCMS, centralityPercentile);
  FillHistogram("hPPBackground", backgroundPP, centralityPercentile);
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

  // Calculate the delta pt
  for(Int_t i=0; i<fRandConeNumber; i++)
  {
    Double_t tmpRatio =1./10.;
    if(fNumberSignalJets)
      tmpRatio =1./fNumberSignalJets;

    Double_t tmpDeltaPtNoBackground = GetDeltaPt(0.0);
    Double_t tmpDeltaPtExternalBgrd = GetDeltaPt(backgroundExternal);
    Double_t tmpDeltaPtExternalBgrdPartialExclusion = GetDeltaPt(backgroundExternal, tmpRatio);
    Double_t tmpDeltaPtPP = GetDeltaPt(backgroundPP);
    Double_t tmpDeltaPtKTImprovedCMS = GetDeltaPt(backgroundKTImprovedCMS);

    Double_t tmpDeltaPtKTPbPb = 0;
    Double_t tmpDeltaPtKTPbPbWithGhosts = 0;
    Double_t tmpDeltaPtKTCMS = 0;
    Double_t tmpDeltaPtKTMean = 0;
    Double_t tmpDeltaPtKTTrackLike = 0;
    Double_t tmpDeltaPtTR = 0;

    tmpDeltaPtKTPbPb = GetDeltaPt(backgroundKTPbPb);
    tmpDeltaPtKTPbPbWithGhosts = GetDeltaPt(backgroundKTPbPbWithGhosts);
    tmpDeltaPtKTCMS = GetDeltaPt(backgroundKTCMS);
    tmpDeltaPtKTMean = GetDeltaPt(backgroundKTMean);
    tmpDeltaPtKTTrackLike = GetDeltaPt(backgroundKTTrackLike);
    tmpDeltaPtTR = GetDeltaPt(backgroundTRCone06);


    // If valid, fill the delta pt histograms

    if(tmpDeltaPtExternalBgrd > -10000.0)
      FillHistogram("hDeltaPtExternalBgrd", tmpDeltaPtExternalBgrd, centralityPercentile);
    if(tmpDeltaPtKTImprovedCMS > -10000.0)
      FillHistogram("hDeltaPtKTImprovedCMS", tmpDeltaPtKTImprovedCMS, centralityPercentile);
    if(tmpDeltaPtExternalBgrdPartialExclusion > -10000.0)
      FillHistogram("hDeltaPtExternalBgrdPartialExclusion", tmpDeltaPtExternalBgrdPartialExclusion, centralityPercentile);
    if(tmpDeltaPtPP > -10000.0)
      FillHistogram("hDeltaPtPP", tmpDeltaPtPP, centralityPercentile);

    if(tmpDeltaPtNoBackground > -10000.0)
      FillHistogram("hDeltaPtNoBackground", tmpDeltaPtNoBackground, centralityPercentile);

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
    if(tmpDeltaPtTR > -10000.0)
      FillHistogram("hDeltaPtTR", tmpDeltaPtTR, centralityPercentile);
  }

  #ifdef DEBUGMODE
    AliInfo("Calculate()::Background done.");
  #endif
  
  #ifdef DEBUGMODE
    AliInfo("Calculate() done.");
  #endif
}

//________________________________________________________________________
Bool_t AliAnalysisTaskChargedJetsPA::UserNotify()
{
  return kTRUE;
}

//________________________________________________________________________
void AliAnalysisTaskChargedJetsPA::CreateJetProfilePlots(Double_t bgrd)
{
  for (Int_t i = 0; i<fJetArray->GetEntries(); i++)
  {
    AliEmcalJet* tmpJet = static_cast<AliEmcalJet*>(fJetArray->At(i));
    if (!tmpJet)
      continue;
    if(!IsSignalJetInAcceptance(tmpJet))
      continue;

    SetCurrentOutputList(1);
    // Jet profile analysis
    if(TMath::Abs(tmpJet->Eta()) <= 0.3)
    {
      if(tmpJet->Pt()>=70.0)
      {
        FillHistogram("hJetProfile70GeV", 0.05-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.05, bgrd))/GetCorrectedJetPt(tmpJet, bgrd));
        FillHistogram("hJetProfile70GeV", 0.10-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.10, bgrd))/GetCorrectedJetPt(tmpJet, bgrd));
        FillHistogram("hJetProfile70GeV", 0.15-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.15, bgrd))/GetCorrectedJetPt(tmpJet, bgrd));
        FillHistogram("hJetProfile70GeV", 0.20-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.20, bgrd))/GetCorrectedJetPt(tmpJet, bgrd));
        FillHistogram("hJetProfile70GeV", 0.25-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.25, bgrd))/GetCorrectedJetPt(tmpJet, bgrd));
        FillHistogram("hJetProfile70GeV", 0.30-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.30, bgrd))/GetCorrectedJetPt(tmpJet, bgrd));
        FillHistogram("hJetProfile70GeV", 0.35-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.35, bgrd))/GetCorrectedJetPt(tmpJet, bgrd));
        FillHistogram("hJetProfile70GeV", 0.40-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.40, bgrd))/GetCorrectedJetPt(tmpJet, bgrd));
        FillHistogram("hJetProfile70GeV", 0.45-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.45, bgrd))/GetCorrectedJetPt(tmpJet, bgrd));
        FillHistogram("hJetProfile70GeV", 0.50-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.50, bgrd))/GetCorrectedJetPt(tmpJet, bgrd));
        FillHistogram("hJetProfile70GeV", 0.55-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.55, bgrd))/GetCorrectedJetPt(tmpJet, bgrd));
        FillHistogram("hJetProfile70GeV", 0.60-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.60, bgrd))/GetCorrectedJetPt(tmpJet, bgrd));
      }
      else if(GetCorrectedJetPt(tmpJet, bgrd)>=60.0)
      {
        FillHistogram("hJetProfile60GeV", 0.05-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.05, bgrd))/GetCorrectedJetPt(tmpJet, bgrd));
        FillHistogram("hJetProfile60GeV", 0.10-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.10, bgrd))/GetCorrectedJetPt(tmpJet, bgrd));
        FillHistogram("hJetProfile60GeV", 0.15-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.15, bgrd))/GetCorrectedJetPt(tmpJet, bgrd));
        FillHistogram("hJetProfile60GeV", 0.20-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.20, bgrd))/GetCorrectedJetPt(tmpJet, bgrd));
        FillHistogram("hJetProfile60GeV", 0.25-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.25, bgrd))/GetCorrectedJetPt(tmpJet, bgrd));
        FillHistogram("hJetProfile60GeV", 0.30-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.30, bgrd))/GetCorrectedJetPt(tmpJet, bgrd));
        FillHistogram("hJetProfile60GeV", 0.35-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.35, bgrd))/GetCorrectedJetPt(tmpJet, bgrd));
        FillHistogram("hJetProfile60GeV", 0.40-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.40, bgrd))/GetCorrectedJetPt(tmpJet, bgrd));
        FillHistogram("hJetProfile60GeV", 0.45-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.45, bgrd))/GetCorrectedJetPt(tmpJet, bgrd));
        FillHistogram("hJetProfile60GeV", 0.50-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.50, bgrd))/GetCorrectedJetPt(tmpJet, bgrd));
        FillHistogram("hJetProfile60GeV", 0.55-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.55, bgrd))/GetCorrectedJetPt(tmpJet, bgrd));
        FillHistogram("hJetProfile60GeV", 0.60-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.60, bgrd))/GetCorrectedJetPt(tmpJet, bgrd));
      }
      else if(GetCorrectedJetPt(tmpJet, bgrd)>=50.0)
      {
        FillHistogram("hJetProfile50GeV", 0.05-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.05, bgrd))/GetCorrectedJetPt(tmpJet, bgrd));
        FillHistogram("hJetProfile50GeV", 0.10-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.10, bgrd))/GetCorrectedJetPt(tmpJet, bgrd));
        FillHistogram("hJetProfile50GeV", 0.15-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.15, bgrd))/GetCorrectedJetPt(tmpJet, bgrd));
        FillHistogram("hJetProfile50GeV", 0.20-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.20, bgrd))/GetCorrectedJetPt(tmpJet, bgrd));
        FillHistogram("hJetProfile50GeV", 0.25-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.25, bgrd))/GetCorrectedJetPt(tmpJet, bgrd));
        FillHistogram("hJetProfile50GeV", 0.30-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.30, bgrd))/GetCorrectedJetPt(tmpJet, bgrd));
        FillHistogram("hJetProfile50GeV", 0.35-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.35, bgrd))/GetCorrectedJetPt(tmpJet, bgrd));
        FillHistogram("hJetProfile50GeV", 0.40-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.40, bgrd))/GetCorrectedJetPt(tmpJet, bgrd));
        FillHistogram("hJetProfile50GeV", 0.45-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.45, bgrd))/GetCorrectedJetPt(tmpJet, bgrd));
        FillHistogram("hJetProfile50GeV", 0.50-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.50, bgrd))/GetCorrectedJetPt(tmpJet, bgrd));
        FillHistogram("hJetProfile50GeV", 0.55-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.55, bgrd))/GetCorrectedJetPt(tmpJet, bgrd));
        FillHistogram("hJetProfile50GeV", 0.60-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.60, bgrd))/GetCorrectedJetPt(tmpJet, bgrd));
      }
      else if(GetCorrectedJetPt(tmpJet, bgrd)>=40.0)
      {
        FillHistogram("hJetProfile40GeV", 0.05-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.05, bgrd))/GetCorrectedJetPt(tmpJet, bgrd));
        FillHistogram("hJetProfile40GeV", 0.10-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.10, bgrd))/GetCorrectedJetPt(tmpJet, bgrd));
        FillHistogram("hJetProfile40GeV", 0.15-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.15, bgrd))/GetCorrectedJetPt(tmpJet, bgrd));
        FillHistogram("hJetProfile40GeV", 0.20-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.20, bgrd))/GetCorrectedJetPt(tmpJet, bgrd));
        FillHistogram("hJetProfile40GeV", 0.25-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.25, bgrd))/GetCorrectedJetPt(tmpJet, bgrd));
        FillHistogram("hJetProfile40GeV", 0.30-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.30, bgrd))/GetCorrectedJetPt(tmpJet, bgrd));
        FillHistogram("hJetProfile40GeV", 0.35-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.35, bgrd))/GetCorrectedJetPt(tmpJet, bgrd));
        FillHistogram("hJetProfile40GeV", 0.40-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.40, bgrd))/GetCorrectedJetPt(tmpJet, bgrd));
        FillHistogram("hJetProfile40GeV", 0.45-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.45, bgrd))/GetCorrectedJetPt(tmpJet, bgrd));
        FillHistogram("hJetProfile40GeV", 0.50-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.50, bgrd))/GetCorrectedJetPt(tmpJet, bgrd));
        FillHistogram("hJetProfile40GeV", 0.55-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.55, bgrd))/GetCorrectedJetPt(tmpJet, bgrd));
        FillHistogram("hJetProfile40GeV", 0.60-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.60, bgrd))/GetCorrectedJetPt(tmpJet, bgrd));
      }
      else if(GetCorrectedJetPt(tmpJet, bgrd)>=30.0)
      {
        FillHistogram("hJetProfile30GeV", 0.05-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.05, bgrd))/GetCorrectedJetPt(tmpJet, bgrd));
        FillHistogram("hJetProfile30GeV", 0.10-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.10, bgrd))/GetCorrectedJetPt(tmpJet, bgrd));
        FillHistogram("hJetProfile30GeV", 0.15-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.15, bgrd))/GetCorrectedJetPt(tmpJet, bgrd));
        FillHistogram("hJetProfile30GeV", 0.20-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.20, bgrd))/GetCorrectedJetPt(tmpJet, bgrd));
        FillHistogram("hJetProfile30GeV", 0.25-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.25, bgrd))/GetCorrectedJetPt(tmpJet, bgrd));
        FillHistogram("hJetProfile30GeV", 0.30-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.30, bgrd))/GetCorrectedJetPt(tmpJet, bgrd));
        FillHistogram("hJetProfile30GeV", 0.35-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.35, bgrd))/GetCorrectedJetPt(tmpJet, bgrd));
        FillHistogram("hJetProfile30GeV", 0.40-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.40, bgrd))/GetCorrectedJetPt(tmpJet, bgrd));
        FillHistogram("hJetProfile30GeV", 0.45-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.45, bgrd))/GetCorrectedJetPt(tmpJet, bgrd));
        FillHistogram("hJetProfile30GeV", 0.50-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.50, bgrd))/GetCorrectedJetPt(tmpJet, bgrd));
        FillHistogram("hJetProfile30GeV", 0.55-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.55, bgrd))/GetCorrectedJetPt(tmpJet, bgrd));
        FillHistogram("hJetProfile30GeV", 0.60-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.60, bgrd))/GetCorrectedJetPt(tmpJet, bgrd));
      }
      else if(GetCorrectedJetPt(tmpJet, bgrd)>=20.0)
      {
        FillHistogram("hJetProfile20GeV", 0.05-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.05, bgrd))/GetCorrectedJetPt(tmpJet, bgrd));
        FillHistogram("hJetProfile20GeV", 0.10-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.10, bgrd))/GetCorrectedJetPt(tmpJet, bgrd));
        FillHistogram("hJetProfile20GeV", 0.15-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.15, bgrd))/GetCorrectedJetPt(tmpJet, bgrd));
        FillHistogram("hJetProfile20GeV", 0.20-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.20, bgrd))/GetCorrectedJetPt(tmpJet, bgrd));
        FillHistogram("hJetProfile20GeV", 0.25-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.25, bgrd))/GetCorrectedJetPt(tmpJet, bgrd));
        FillHistogram("hJetProfile20GeV", 0.30-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.30, bgrd))/GetCorrectedJetPt(tmpJet, bgrd));
        FillHistogram("hJetProfile20GeV", 0.35-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.35, bgrd))/GetCorrectedJetPt(tmpJet, bgrd));
        FillHistogram("hJetProfile20GeV", 0.40-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.40, bgrd))/GetCorrectedJetPt(tmpJet, bgrd));
        FillHistogram("hJetProfile20GeV", 0.45-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.45, bgrd))/GetCorrectedJetPt(tmpJet, bgrd));
        FillHistogram("hJetProfile20GeV", 0.50-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.50, bgrd))/GetCorrectedJetPt(tmpJet, bgrd));
        FillHistogram("hJetProfile20GeV", 0.55-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.55, bgrd))/GetCorrectedJetPt(tmpJet, bgrd));
        FillHistogram("hJetProfile20GeV", 0.60-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.60, bgrd))/GetCorrectedJetPt(tmpJet, bgrd));
      }
      else if(GetCorrectedJetPt(tmpJet, bgrd)>=10.0)
      {
        FillHistogram("hJetProfile10GeV", 0.05-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.05, bgrd))/GetCorrectedJetPt(tmpJet, bgrd));
        FillHistogram("hJetProfile10GeV", 0.10-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.10, bgrd))/GetCorrectedJetPt(tmpJet, bgrd));
        FillHistogram("hJetProfile10GeV", 0.15-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.15, bgrd))/GetCorrectedJetPt(tmpJet, bgrd));
        FillHistogram("hJetProfile10GeV", 0.20-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.20, bgrd))/GetCorrectedJetPt(tmpJet, bgrd));
        FillHistogram("hJetProfile10GeV", 0.25-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.25, bgrd))/GetCorrectedJetPt(tmpJet, bgrd));
        FillHistogram("hJetProfile10GeV", 0.30-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.30, bgrd))/GetCorrectedJetPt(tmpJet, bgrd));
        FillHistogram("hJetProfile10GeV", 0.35-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.35, bgrd))/GetCorrectedJetPt(tmpJet, bgrd));
        FillHistogram("hJetProfile10GeV", 0.40-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.40, bgrd))/GetCorrectedJetPt(tmpJet, bgrd));
        FillHistogram("hJetProfile10GeV", 0.45-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.45, bgrd))/GetCorrectedJetPt(tmpJet, bgrd));
        FillHistogram("hJetProfile10GeV", 0.50-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.50, bgrd))/GetCorrectedJetPt(tmpJet, bgrd));
        FillHistogram("hJetProfile10GeV", 0.55-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.55, bgrd))/GetCorrectedJetPt(tmpJet, bgrd));
        FillHistogram("hJetProfile10GeV", 0.60-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.60, bgrd))/GetCorrectedJetPt(tmpJet, bgrd));
      }
    }
    SetCurrentOutputList(0);
  }
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
  TH1* tmpHist = static_cast<TH1*>(fCurrentOutputList->FindObject(GetHistoName(key)));
  if(!tmpHist)
  {
    AliError(Form("Cannot find histogram <%s> ",key)) ;
    return;
  }

  tmpHist->Fill(x);
}

//________________________________________________________________________
inline void AliAnalysisTaskChargedJetsPA::FillHistogram(const char * key, Double_t x, Double_t y)
{
  TH1* tmpHist = static_cast<TH1*>(fCurrentOutputList->FindObject(GetHistoName(key)));
  if(!tmpHist)
  {
    AliError(Form("Cannot find histogram <%s> ",key));
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
  TH2* tmpHist = static_cast<TH2*>(fCurrentOutputList->FindObject(GetHistoName(key)));
  if(!tmpHist)
  {
    AliError(Form("Cannot find histogram <%s> ",key));
    return;
  }
  
  tmpHist->Fill(x,y,add);
}

//________________________________________________________________________
inline void AliAnalysisTaskChargedJetsPA::FillCutHistogram(const char * key, Double_t cut, Double_t pT, Double_t eta, Double_t phi, Int_t isAdditionalTrack)
{
  THnF* tmpHist = static_cast<THnF*>(fCurrentOutputList->FindObject(GetHistoName(key)));
  if(!tmpHist)
  {
    AliError(Form("Cannot find histogram <%s> ",key));
    return;
  }

  Double_t tmpVec[5] = {cut, pT, eta, phi, static_cast<Double_t>(isAdditionalTrack)};
  tmpHist->Fill(tmpVec);
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

  fCurrentOutputList->Add(tmpHist);

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

  fCurrentOutputList->Add(tmpHist);

  return tmpHist;
}

//________________________________________________________________________
THnF* AliAnalysisTaskChargedJetsPA::AddCutHistogram(const char* name, const char* title, const char* cutName, Int_t nBins, Double_t xMin, Double_t xMax)
{
  //                        Cut,      pT,  eta,           phi,  type 
  Int_t    bins [5]     = { nBins,   100,   20,            18,     2};
  Double_t minEdges[5]  = { xMin,    0.1,   -1,             0,  -0.5};
  Double_t maxEdges[5]  = { xMax,     40,   +1, 2*TMath::Pi(),   1.5};

  TString axisName[5]  = {cutName,"#it{p}_{T}","#eta","#phi","Track type"};
  TString axisTitle[5] = {cutName,"#it{p}_{T}","#eta","#phi","Track type"};

  THnF * histo = new THnF(name, title, 5, bins, minEdges, maxEdges);
  BinLogAxis(histo, 1);

  for (Int_t iaxis=0; iaxis<5;iaxis++){
    histo->GetAxis(iaxis)->SetName(axisName[iaxis]);
    histo->GetAxis(iaxis)->SetTitle(axisTitle[iaxis]);
  }

  fCurrentOutputList->Add(histo);
  return histo;
}

//________________________________________________________________________
void AliAnalysisTaskChargedJetsPA::BinLogAxis(const THn *h, Int_t axisNumber)
{
  // Method for the correct logarithmic binning of histograms
  TAxis *axis = h->GetAxis(axisNumber);
  int bins = axis->GetNbins();

  Double_t from = axis->GetXmin();
  Double_t to = axis->GetXmax();
  Double_t *newBins = new Double_t[bins + 1];
   
  newBins[0] = from;
  Double_t factor = pow(to/from, 1./bins);
  
  for (int i = 1; i <= bins; i++) {
   newBins[i] = factor * newBins[i-1];
  }
  axis->Set(bins, newBins);
  delete [] newBins;
}

//________________________________________________________________________
void AliAnalysisTaskChargedJetsPA::Terminate(Option_t *)
{
  if(fNoTerminate)
    return;

  fOutputLists[0] = dynamic_cast<TList*> (GetOutputData(1)); // >1 refers to output slots
  PostData(1, fOutputLists[0]);

  if(fAnalyzeJetProfile)
  {
    fOutputLists[1] = dynamic_cast<TList*> (GetOutputData(2)); // >1 refers to output slots
    PostData(2, fOutputLists[1]);
  }
  if(fAnalyzeTrackcuts)
  {
    if(fAnalyzeJetProfile)
    {
      fOutputLists[2] = dynamic_cast<TList*> (GetOutputData(3)); // >1 refers to output slots
      PostData(3, fOutputLists[2]);
    }
    else
    {
      fOutputLists[1] = dynamic_cast<TList*> (GetOutputData(2)); // >1 refers to output slots
      PostData(2, fOutputLists[1]);
    }
  }
}

//________________________________________________________________________
AliAnalysisTaskChargedJetsPA::~AliAnalysisTaskChargedJetsPA()
{
  // Destructor. Clean-up the output list, but not the histograms that are put inside
  // (the list is owner and will clean-up these histograms). Protect in PROOF case.

  if(fNoTerminate)
    return;

  delete fHybridESDtrackCuts;
  delete fHybridESDtrackCuts_variedPtDep;

  for(Int_t i=0; i<static_cast<Int_t>(fOutputLists.size()); i++)
    if (fOutputLists[i] && !AliAnalysisManager::GetAnalysisManager()->IsProofMode())
      delete fOutputLists[i];

}

//________________________________________________________________________
void AliAnalysisTaskChargedJetsPA::UserCreateOutputObjects()
{
  // called once to create user defined output objects like histograms, plots etc. 
  // and to put it on the output list.
  // Note: Saving to file with e.g. OpenFile(0) is must be before creating other objects.

  fRandom = new TRandom3(0);


  Int_t tmpListCount = 1;
  if(fAnalyzeJetProfile)
    tmpListCount++;
  if(fAnalyzeTrackcuts)
    tmpListCount++;

  fOutputLists.resize(tmpListCount);
  for(Int_t i=0; i<tmpListCount; i++)
  {
    fOutputLists[i] = new TList();
    fOutputLists[i]->SetOwner(); // otherwise it produces leaks in merging
    PostData(i+1, fOutputLists[i]);
  }
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

  PostData(1, fOutputLists[0]);
  if(fAnalyzeJetProfile)
    PostData(2, fOutputLists[1]);
  if(fAnalyzeTrackcuts)
  {
    if(fAnalyzeJetProfile)
      PostData(3, fOutputLists[2]);
    else
      PostData(2, fOutputLists[1]);
  }

}
