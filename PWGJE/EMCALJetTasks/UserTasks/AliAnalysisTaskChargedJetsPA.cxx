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

//TODO: Not accessing the particles when using MC
//TODO: FillHistogram can be done better with virtual TH1(?)
ClassImp(AliAnalysisTaskChargedJetsPA)

// ######################################################################################## DEFINE HISTOGRAMS
void AliAnalysisTaskChargedJetsPA::Init()
{
  #ifdef DEBUGMODE
    AliInfo("Creating histograms.");
  #endif

  TH1D* tmpHisto = AddHistogram1D<TH1D>("hNumberEvents", "Number of events (0 = before cuts, 1 = after cuts)", "", 2, 0, 2, "stage","N^{Events}/cut");
  tmpHisto->GetXaxis()->SetBinLabel(1, "Before cuts");
  tmpHisto->GetXaxis()->SetBinLabel(2, "After cuts");

  tmpHisto = AddHistogram1D<TH1D>("hEventAcceptance", "Accepted events (0 = before cuts, 1 = after pile up, 2 = after vertex)", "", 3, 0, 3, "stage","N^{Events}/cut");
  tmpHisto->GetXaxis()->SetBinLabel(1, "Before cuts");
  tmpHisto->GetXaxis()->SetBinLabel(2, "After pile up");
  tmpHisto->GetXaxis()->SetBinLabel(3, "After vertex");

  tmpHisto = AddHistogram1D<TH1D>("hTrackAcceptance", "Accepted tracks (0 = before cuts, 1 = after eta, 2 = after pT)", "", 3, 0, 3, "stage","N^{Tracks}/cut");
  tmpHisto->GetXaxis()->SetBinLabel(1, "Before cuts");
  tmpHisto->GetXaxis()->SetBinLabel(2, "After eta");
  tmpHisto->GetXaxis()->SetBinLabel(3, "After p_{T}");

  tmpHisto = AddHistogram1D<TH1D>("hJetAcceptance", "Accepted jets (0 = before cuts, 1 = after eta, 2 = after pT, 3 = after area)", "", 4, 0, 4, "stage","N^{Jets}/cut");
  tmpHisto->GetXaxis()->SetBinLabel(1, "Before cuts");
  tmpHisto->GetXaxis()->SetBinLabel(2, "After eta");
  tmpHisto->GetXaxis()->SetBinLabel(3, "After p_{T}");
  tmpHisto->GetXaxis()->SetBinLabel(4, "After area");

  // NOTE: Jet histograms
  if (fAnalyzeJets)
  {
    // ######## Jet spectra
    AddHistogram1D<TH1D>("hRawJetPt", "Raw jets p_{T} distribution (before cuts)", "", 500, 0., 250., "p_{T} (GeV/c)", "dN^{Jets}/dp_{T}");
    AddHistogram2D<TH2D>("hJetPt", "Jets p_{T} distribution", "", 500, -50., 200., fNumberOfCentralityBins, 0, 100, "p_{T} (GeV/c)","Centrality","dN^{Jets}/dp_{T}");
    AddHistogram2D<TH2D>("hJetPtBgrdSubtractedKTImprovedCMS", "Jets p_{T} distribution, KT background (Improved CMS) subtracted", "", 500, -50., 200., fNumberOfCentralityBins, 0, 100, "p_{T} (GeV/c)","Centrality","dN^{Jets}/dp_{T}");
    AddHistogram2D<TH2D>("hJetPtBgrdSubtractedKTImprovedCMS_Phi1", "Jets p_{T} distribution, KT background (Improved CMS) subtracted (1st part of azimuth)", "", 500, -50., 200., fNumberOfCentralityBins, 0, 100, "p_{T} (GeV/c)","Centrality","dN^{Jets}/dp_{T}");    
    AddHistogram2D<TH2D>("hJetPtBgrdSubtractedKTImprovedCMS_Phi2", "Jets p_{T} distribution, KT background (Improved CMS) subtracted (2nd part of azimuth)", "", 500, -50., 200., fNumberOfCentralityBins, 0, 100, "p_{T} (GeV/c)","Centrality","dN^{Jets}/dp_{T}");    

    AddHistogram2D<TH2D>("hJetPtBgrdSubtractedPP", "Jets p_{T} distribution, pp background subtracted", "", 500, -50., 200., fNumberOfCentralityBins, 0, 100, "p_{T} (GeV/c)","Centrality","dN^{Jets}/dp_{T}");
    AddHistogram2D<TH2D>("hJetPtBgrdSubtractedExternal", "Jets p_{T} distribution, external bgrd. subtracted", "", 500, -50., 200., fNumberOfCentralityBins, 0, 100, "p_{T} (GeV/c)","Centrality","dN^{Jets}/dp_{T}");    



    AddHistogram2D<TProfile2D>("hJetPtSubtractedRhoKTImprovedCMS", "Mean subtracted KT (CMS w/o signal) background from jets", "COLZ", 600, 0, 150, fNumberOfCentralityBins, 0, 100, "Jet p_{T}", "Centrality", "#rho mean");
    AddHistogram2D<TProfile2D>("hJetPtSubtractedRhoExternal", "Mean subtracted KT (External) background from jets", "COLZ", 600, 0, 150, fNumberOfCentralityBins, 0, 100, "Jet p_{T}", "Centrality", "#rho mean");
    AddHistogram2D<TProfile2D>("hJetPtSubtractedRhoPP", "Mean subtracted KT (pp from Michal) background from jets", "COLZ", 600, 0, 150, fNumberOfCentralityBins, 0, 100, "Jet p_{T}", "Centrality", "#rho mean");

    if(fAnalyzeDeprecatedBackgrounds)
    {
      AddHistogram2D<TH2D>("hJetPtBgrdSubtractedTR", "Jets p_{T} distribution, TR background (Cone R=0.6 around jets excluded) subtracted", "", 500, -50., 200., fNumberOfCentralityBins, 0, 100, "p_{T} (GeV/c)","Centrality","dN^{Jets}/dp_{T}");
      AddHistogram2D<TH2D>("hJetPtBgrdSubtractedKTPbPb", "Jets p_{T} distribution, KT background (PbPb w/o ghosts) subtracted", "", 500, -50., 200., fNumberOfCentralityBins, 0, 100, "p_{T} (GeV/c)","Centrality","dN^{Jets}/dp_{T}");
      AddHistogram2D<TH2D>("hJetPtBgrdSubtractedKTPbPbWithGhosts", "Jets p_{T} distribution, KT background (PbPb w/ ghosts) subtracted", "", 500, -50., 200., fNumberOfCentralityBins, 0, 100, "p_{T} (GeV/c)","Centrality","dN^{Jets}/dp_{T}");
      AddHistogram2D<TH2D>("hJetPtBgrdSubtractedKTCMS", "Jets p_{T} distribution, KT background (CMS) subtracted", "", 500, -50., 200., fNumberOfCentralityBins, 0, 100, "p_{T} (GeV/c)","Centrality","dN^{Jets}/dp_{T}");    
      AddHistogram2D<TH2D>("hJetPtBgrdSubtractedKTMean", "Jets p_{T} distribution, KT background (Mean) subtracted", "", 500, -50., 200., fNumberOfCentralityBins, 0, 100, "p_{T} (GeV/c)","Centrality","dN^{Jets}/dp_{T}");    
      AddHistogram2D<TH2D>("hJetPtBgrdSubtractedKTTrackLike", "Jets p_{T} distribution, KT background (track-like) subtracted", "", 500, -50., 200., fNumberOfCentralityBins, 0, 100, "p_{T} (GeV/c)","Centrality","dN^{Jets}/dp_{T}");
    }

    // ######## Jet profiles
    if(fAnalyzeJetProfile)
    {
      AddHistogram2D<TH2D>("hJetProfile10GeV", "Jet profile, cone p_{T}/jet p_{T} vs. jet radius, jet p_{T} > 10 GeV", "", 12, 0, 0.6,200, 0., 2., "Cone radius","dN^{Jets}/dR", "Ratio");
      AddHistogram2D<TH2D>("hJetProfile20GeV", "Jet profile, cone p_{T}/jet p_{T} vs. jet radius, jet p_{T} > 20 GeV", "", 12, 0, 0.6,200, 0., 2., "Cone radius","dN^{Jets}/dR", "Ratio");
      AddHistogram2D<TH2D>("hJetProfile30GeV", "Jet profile, cone p_{T}/jet p_{T} vs. jet radius, jet p_{T} > 30 GeV", "", 12, 0, 0.6,200, 0., 2., "Cone radius","dN^{Jets}/dR", "Ratio");
      AddHistogram2D<TH2D>("hJetProfile40GeV", "Jet profile, cone p_{T}/jet p_{T} vs. jet radius, jet p_{T} > 40 GeV", "", 12, 0, 0.6,200, 0., 2., "Cone radius","dN^{Jets}/dR", "Ratio");
      AddHistogram2D<TH2D>("hJetProfile50GeV", "Jet profile, cone p_{T}/jet p_{T} vs. jet radius, jet p_{T} > 50 GeV", "", 12, 0, 0.6,200, 0., 2., "Cone radius","dN^{Jets}/dR", "Ratio");
      AddHistogram2D<TH2D>("hJetProfile60GeV", "Jet profile, cone p_{T}/jet p_{T} vs. jet radius, jet p_{T} > 60 GeV", "", 12, 0, 0.6,200, 0., 2., "Cone radius","dN^{Jets}/dR", "Ratio");
      AddHistogram2D<TH2D>("hJetProfile70GeV", "Jet profile, cone p_{T}/jet p_{T} vs. jet radius, jet p_{T} > 70 GeV", "", 12, 0, 0.6,200, 0., 2., "Cone radius","dN^{Jets}/dR", "Ratio");
    }
    // ######## Jet stuff
    AddHistogram2D<TH2D>("hJetConstituentPt", "Jet constituents p_{T} distribution", "", 500, -50., 200., fNumberOfCentralityBins, 0, 100, "p_{T} (GeV/c)","Centrality","dN^{Tracks}/dp_{T}");
    AddHistogram1D<TH1D>("hJetCountAll", "Number of Jets", "", 200, 0., 200., "N jets","dN^{Events}/dN^{Jets}");
    AddHistogram1D<TH1D>("hJetCountAccepted", "Number of accepted Jets", "", 200, 0., 200., "N jets","dN^{Events}/dN^{Jets}");
    AddHistogram2D<TH2D>("hJetCount", "Correlation jets/accepted jets", "", 200, 0., 200., 200, 0., 200., "N jets","N jets accepted", "d^{2}N^{Events}/dN^{Jets dN^{Jets, acc}}");
    AddHistogram1D<TH1D>("hLeadingJetPt", "Leading jet p_{T}", "", 500, -50., 200., "p_{T} (GeV/c)","dN^{Jets}/dp_{T}");
    AddHistogram1D<TH1D>("hSecondLeadingJetPt", "Second leading jet p_{T}", "", 500, -50., 200., "p_{T} (GeV/c)","dN^{Jets}/dp_{T}");
    AddHistogram1D<TH1D>("hCorrectedLeadingJetPt", "Corrected leading jet p_{T}", "", 500, -50., 200., "p_{T} (GeV/c)","dN^{Jets}/dp_{T}");
    AddHistogram1D<TH1D>("hCorrectedSecondLeadingJetPt", "Corrected second leading jet p_{T}", "", 500, -50., 200., "p_{T} (GeV/c)","dN^{Jets}/dp_{T}");
    AddHistogram1D<TH1D>("hJetDeltaPhi", "Jets combinatorial #Delta #phi", "", 250, 0., TMath::Pi(), "#Delta #phi","dN^{Jets}/d(#Delta #phi)");
    AddHistogram1D<TH1D>("hLeadingJetDeltaPhi", "1st and 2nd leading jet #Delta #phi", "", 250, 0., TMath::Pi(), "#Delta #phi","dN^{Jets}/d(#Delta #phi)");
  }

  // NOTE: Jet background histograms
  if (fAnalyzeBackground)
  {
    // ########## Default background estimates
    AddHistogram2D<TH2D>("hKTBackgroundImprovedCMS", "KT background density (Improved CMS approach)", "LEGO2", 400, 0., 40., fNumberOfCentralityBins, 0, 100, "#rho (GeV/c)","Centrality", "dN^{Events}/d#rho");
    AddHistogram2D<TH2D>("hKTBackgroundImprovedCMSExternal", "KT background density (Improved CMS approach from external task)", "LEGO2", 400, 0., 40., fNumberOfCentralityBins, 0, 100, "#rho (GeV/c)","Centrality", "dN^{Events}/d#rho");
    AddHistogram2D<TH2D>("hPPBackground", "PP background density (Michals approach)", "LEGO2", 400, 0., 40., fNumberOfCentralityBins, 0, 100, "#rho (GeV/c)","Centrality", "dN^{Events}/d#rho");

    AddHistogram2D<TH2D>("hDeltaPtExternalBgrd", "Background fluctuations #delta p_{T} (KT, External)", "", 1201, -40.0, 40.0, fNumberOfCentralityBins, 0, 100, "#delta p_{T} (GeV/c)","Centrality","dN^{Jets}/d#delta p_{T}");
    AddHistogram2D<TH2D>("hDeltaPtKTImprovedCMS", "Background fluctuations #delta p_{T} (KT, Improved CMS-like)", "", 1201, -40.0, 40.0, fNumberOfCentralityBins, 0, 100, "#delta p_{T} (GeV/c)","Centrality","dN^{Jets}/d#delta p_{T}");
    AddHistogram2D<TH2D>("hDeltaPtKTImprovedCMSPartialExclusion", "Background fluctuations #delta p_{T} (KT, Improved CMS-like, partial jet exclusion)", "", 1201, -40.0, 40.0, fNumberOfCentralityBins, 0, 100, "#delta p_{T} (GeV/c)","Centrality","dN^{Jets}/d#delta p_{T}");
    AddHistogram2D<TH2D>("hDeltaPtKTImprovedCMSPartialExclusion_Signal", "Background fluctuations #delta p_{T} (KT, Improved CMS-like, partial jet exclusion w/ 1/N_{sig} probability)", "", 1201, -40.0, 40.0, fNumberOfCentralityBins, 0, 100, "#delta p_{T} (GeV/c)","Centrality","dN^{Jets}/d#delta p_{T}");
    AddHistogram2D<TH2D>("hDeltaPtKTImprovedCMSFullExclusion", "Background fluctuations #delta p_{T} (KT, Improved CMS-like, full leading jet exclusion)", "", 1201, -40.0, 40.0, fNumberOfCentralityBins, 0, 100, "#delta p_{T} (GeV/c)","Centrality","dN^{Jets}/d#delta p_{T}");
    AddHistogram2D<TH2D>("hDeltaPtNoBackground", "Background fluctuations #delta p_{T} (No background)", "", 1201, -40.0, 40.0, fNumberOfCentralityBins, 0, 100, "#delta p_{T} (GeV/c)","Centrality","dN^{Jets}/d#delta p_{T}");
    AddHistogram2D<TH2D>("hDeltaPtNoBackgroundNoEmptyCones", "Background fluctuations #delta p_{T} (No background, no empty cones)", "", 1201, -40.0, 40.0, fNumberOfCentralityBins, 0, 100, "#delta p_{T} (GeV/c)","Centrality","dN^{Jets}/d#delta p_{T}");

    AddHistogram1D<TProfile>("hKTMeanBackgroundImprovedCMS", "KT background mean (Improved CMS approach)", "", 100, 0, 100, "Centrality", "#rho mean");

    if(fAnalyzeDeprecatedBackgrounds)
    {
      // ########## Different background estimates
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

      // ########## Profiles for background means vs. centrality
      AddHistogram1D<TProfile>("hKTMeanBackgroundPbPb", "KT background mean (PbPb approach w/o ghosts)", "", fNumberOfCentralityBins, 0, 100, "Centrality", "#rho mean");
      AddHistogram1D<TProfile>("hKTMeanBackgroundPbPbWithGhosts", "KT background mean (PbPb approach)", "", fNumberOfCentralityBins, 0, 100, "Centrality", "#rho mean");
      AddHistogram1D<TProfile>("hKTMeanBackgroundCMS", "KT background mean (CMS approach)", "", fNumberOfCentralityBins, 0, 100, "Centrality", "#rho mean");
      AddHistogram1D<TProfile>("hKTMeanBackgroundMean", "KT background mean (Mean approach)", "",  fNumberOfCentralityBins, 0, 100, "Centrality", "#rho mean");
      AddHistogram1D<TProfile>("hKTMeanBackgroundTPC", "KT background mean (Track-like approach)", "", fNumberOfCentralityBins, 0, 100, "Centrality", "#rho mean");
      AddHistogram1D<TProfile>("hTRMeanBackground", "TR background mean", "", fNumberOfCentralityBins, 0, 100, "Centrality", "#rho mean");
    }
  }
  // NOTE: Jet constituent correlations
  if(fAnalyzeMassCorrelation)
  {
    AddHistogram1D<TH1D>("hJetMassFromConstituents", "Jet mass by mass of charged constituents", "", 200, 0., 200., "N jets","dN^{Jets}/dN^{mass}");
    AddHistogram1D<TH1D>("hJetMass", "Jet mass from fastjet", "", 200, 0., 200., "N jets","dN^{Jets}/dN^{mass}");

    AddHistogram2D<TH2D>("hJetPtVsMass", "Correlation jet pt/summed constituent jet", "", 400, 0., 100., 400, 0., 100., "p_{T}","Mass", "d^{2}N}/dN^{p_{T}dM");
    AddHistogram2D<TH2D>("hJetPtVsJetMass", "Correlation jet pt/summed jet", "", 400, 0., 100., 400, 0., 100., "p_{T}","Mass", "d^{2}N}/dN^{p_{T}dM");

    AddHistogram2D<TH2D>("hJetPtVsProtonCount", "Correlation jet pt/amount of proton in jet", "", 400, 0., 100., 20, 0., 1., "p_{T}","Proton percentage", "d^{2}N}/dN^{p_{T}dPercentage");
    AddHistogram2D<TH2D>("hJetPtVsPionCount", "Correlation jet pt/amount of pion in jet", "", 400, 0., 100., 20, 0., 1., "p_{T}","Pion percentage", "d^{2}N}/dN^{p_{T}dPercentage");
    AddHistogram2D<TH2D>("hJetPtVsKaonCount", "Correlation jet pt/amount of kaon in jet", "", 400, 0., 100., 20, 0., 1., "p_{T}","Kaon percentage", "d^{2}N}/dN^{p_{T}dPercentage");
    AddHistogram2D<TH2D>("hJetPtVsElectronCount", "Correlation jet pt/amount of proton in jet", "", 400, 0., 100., 20, 0., 1., "p_{T}","Electron percentage", "d^{2}N}/dN^{p_{T}dPercentage");
    AddHistogram2D<TH2D>("hJetPtVsOthersCount", "Correlation jet pt/amount of other particles in jet", "", 400, 0., 100., 20, 0., 1., "p_{T}","Others percentage", "d^{2}N}/dN^{p_{T}dPercentage");
    AddHistogram1D<TH1D>("hJetConstituentProtonCount", "Proton count in jets", "", 100, 0., 100., "N protons","dN^{Jets}/dN^{protons}");
    AddHistogram1D<TH1D>("hJetConstituentPionCount", "Pion count in jets", "", 100, 0., 100., "N pions","dN^{Jets}/dN^{pions}");
    AddHistogram1D<TH1D>("hJetConstituentKaonCount", "Kaon count in jets", "", 100, 0., 100., "N kaons","dN^{Jets}/dN^{kaons}");
    AddHistogram1D<TH1D>("hJetConstituentElectronCount", "Electron count in jets", "", 100, 0., 100., "N electrons","dN^{Jets}/dN^{electrons}");
    AddHistogram1D<TH1D>("hJetConstituentOthersCount", "Others count in jets", "", 100, 0., 100., "N others","dN^{Jets}/dN^{others}");

    AddHistogram2D<TH2D>("hJetPtVsMass_6_14", "Correlation jet pt/summed constituent jet", "", 400, 0., 100., 400, 0., 100., "p_{T}","Mass", "d^{2}N}/dN^{p_{T}dM");
    AddHistogram2D<TH2D>("hJetPtVsJetMass_6_14", "Correlation jet pt/summed jet", "", 400, 0., 100., 400, 0., 100., "p_{T}","Mass", "d^{2}N}/dN^{p_{T}dM");
    AddHistogram2D<TH2D>("hJetPtVsProtonCount_6_14", "Correlation jet pt/amount of proton in jet", "", 400, 0., 100., 8, 6., 14., "p_{T}","Proton count", "d^{2}N}/dN^{p_{T}dN");
    AddHistogram2D<TH2D>("hJetPtVsPionCount_6_14", "Correlation jet pt/amount of pion in jet", "", 400, 0., 100., 8, 6., 14., "p_{T}","Pion count", "d^{2}N}/dN^{p_{T}dN");
    AddHistogram2D<TH2D>("hJetPtVsKaonCount_6_14", "Correlation jet pt/amount of kaon in jet", "", 400, 0., 100., 8, 6., 14., "p_{T}","Kaon count", "d^{2}N}/dN^{p_{T}dN");
    AddHistogram2D<TH2D>("hJetPtVsElectronCount_6_14", "Correlation jet pt/amount of proton in jet", "", 400, 0., 100., 8, 6., 14., "p_{T}","Electron count", "d^{2}N}/dN^{p_{T}dN");
    AddHistogram2D<TH2D>("hJetPtVsOthersCount_6_14", "Correlation jet pt/amount of other particles in jet", "", 400, 0., 100., 8, 6., 14., "p_{T}","Others count", "d^{2}N}/dN^{p_{T}dN");
    AddHistogram1D<TH1D>("hJetConstituentProtonCount_6_14", "Proton count in jets", "", 100, 0., 100., "N protons","dN^{Jets}/dN^{protons}");
    AddHistogram1D<TH1D>("hJetConstituentPionCount_6_14", "Pion count in jets", "", 100, 0., 100., "N pions","dN^{Jets}/dN^{pions}");
    AddHistogram1D<TH1D>("hJetConstituentKaonCount_6_14", "Kaon count in jets", "", 100, 0., 100., "N kaons","dN^{Jets}/dN^{kaons}");
    AddHistogram1D<TH1D>("hJetConstituentElectronCount_6_14", "Electron count in jets", "", 100, 0., 100., "N electrons","dN^{Jets}/dN^{electrons}");
    AddHistogram1D<TH1D>("hJetConstituentOthersCount_6_14", "Others count in jets", "", 100, 0., 100., "N others","dN^{Jets}/dN^{others}");

    AddHistogram2D<TH2D>("hJetPtVsMass_2_X", "Correlation jet pt/summed constituent jet", "", 400, 0., 100., 400, 0., 100., "p_{T}","Mass", "d^{2}N}/dN^{p_{T}dM");
    AddHistogram2D<TH2D>("hJetPtVsJetMass_2_X", "Correlation jet pt/summed jet", "", 400, 0., 100., 400, 0., 100., "p_{T}","Mass", "d^{2}N}/dN^{p_{T}dM");
    AddHistogram2D<TH2D>("hJetPtVsProtonCount_2_X", "Correlation jet pt/amount of proton in jet", "", 400, 0., 100., 68, 2., 70., "p_{T}","Proton count", "d^{2}N}/dN^{p_{T}dN");
    AddHistogram2D<TH2D>("hJetPtVsPionCount_2_X", "Correlation jet pt/amount of pion in jet", "", 400, 0., 100., 68, 2., 70., "p_{T}","Pion count", "d^{2}N}/dN^{p_{T}dN");
    AddHistogram2D<TH2D>("hJetPtVsKaonCount_2_X", "Correlation jet pt/amount of kaon in jet", "", 400, 0., 100., 68, 2., 70., "p_{T}","Kaon count", "d^{2}N}/dN^{p_{T}dN");
    AddHistogram2D<TH2D>("hJetPtVsElectronCount_2_X", "Correlation jet pt/amount of proton in jet", "", 400, 0., 100., 68, 2., 70., "p_{T}","Electron count", "d^{2}N}/dN^{p_{T}dN");
    AddHistogram2D<TH2D>("hJetPtVsOthersCount_2_X", "Correlation jet pt/amount of other particles in jet", "", 400, 0., 100., 68, 2., 70., "p_{T}","Others count", "d^{2}N}/dN^{p_{T}dN");
    AddHistogram1D<TH1D>("hJetConstituentProtonCount_2_X", "Proton count in jets", "", 100, 0., 100., "N protons","dN^{Jets}/dN^{protons}");
    AddHistogram1D<TH1D>("hJetConstituentPionCount_2_X", "Pion count in jets", "", 100, 0., 100., "N pions","dN^{Jets}/dN^{pions}");
    AddHistogram1D<TH1D>("hJetConstituentKaonCount_2_X", "Kaon count in jets", "", 100, 0., 100., "N kaons","dN^{Jets}/dN^{kaons}");
    AddHistogram1D<TH1D>("hJetConstituentElectronCount_2_X", "Electron count in jets", "", 100, 0., 100., "N electrons","dN^{Jets}/dN^{electrons}");
    AddHistogram1D<TH1D>("hJetConstituentOthersCount_2_X", "Others count in jets", "", 100, 0., 100., "N others","dN^{Jets}/dN^{others}");

    AddHistogram2D<TH2D>("hJetPtVsMass_2_7", "Correlation jet pt/summed constituent jet", "", 400, 0., 100., 400, 0., 100., "p_{T}","Mass", "d^{2}N}/dN^{p_{T}dM");
    AddHistogram2D<TH2D>("hJetPtVsJetMass_2_7", "Correlation jet pt/summed jet", "", 400, 0., 100., 400, 0., 100., "p_{T}","Mass", "d^{2}N}/dN^{p_{T}dM");
    AddHistogram2D<TH2D>("hJetPtVsProtonCount_2_7", "Correlation jet pt/amount of proton in jet", "", 400, 0., 100., 6, 2., 8., "p_{T}","Proton count", "d^{2}N}/dN^{p_{T}dN");
    AddHistogram2D<TH2D>("hJetPtVsPionCount_2_7", "Correlation jet pt/amount of pion in jet", "", 400, 0., 100., 6, 2., 8., "p_{T}","Pion count", "d^{2}N}/dN^{p_{T}dN");
    AddHistogram2D<TH2D>("hJetPtVsKaonCount_2_7", "Correlation jet pt/amount of kaon in jet", "", 400, 0., 100., 6, 2., 8., "p_{T}","Kaon count", "d^{2}N}/dN^{p_{T}dN");
    AddHistogram2D<TH2D>("hJetPtVsElectronCount_2_7", "Correlation jet pt/amount of proton in jet", "", 400, 0., 100., 6, 2., 8., "p_{T}","Electron count", "d^{2}N}/dN^{p_{T}dN");
    AddHistogram2D<TH2D>("hJetPtVsOthersCount_2_7", "Correlation jet pt/amount of other particles in jet", "", 400, 0., 100., 6, 2., 8., "p_{T}","Others count", "d^{2}N}/dN^{p_{T}dN");
    AddHistogram1D<TH1D>("hJetConstituentProtonCount_2_7", "Proton count in jets", "", 100, 0., 100., "N protons","dN^{Jets}/dN^{protons}");
    AddHistogram1D<TH1D>("hJetConstituentPionCount_2_7", "Pion count in jets", "", 100, 0., 100., "N pions","dN^{Jets}/dN^{pions}");
    AddHistogram1D<TH1D>("hJetConstituentKaonCount_2_7", "Kaon count in jets", "", 100, 0., 100., "N kaons","dN^{Jets}/dN^{kaons}");
    AddHistogram1D<TH1D>("hJetConstituentElectronCount_2_7", "Electron count in jets", "", 100, 0., 100., "N electrons","dN^{Jets}/dN^{electrons}");
    AddHistogram1D<TH1D>("hJetConstituentOthersCount_2_7", "Others count in jets", "", 100, 0., 100., "N others","dN^{Jets}/dN^{others}");

  }

  // NOTE: Track & Cluster & QA histograms
  if (fAnalyzeQA)
  {
    AddHistogram1D<TH1D>("hVertexX", "X distribution of the vertex", "", 2000, -1., 1., "#Delta x(cm)","dN^{Events}/dx");
    AddHistogram1D<TH1D>("hVertexY", "Y distribution of the vertex", "", 2000, -1., 1., "#Delta y(cm)","dN^{Events}/dy");
    AddHistogram2D<TH2D>("hVertexXY", "XY distribution of the vertex", "COLZ", 500, -1., 1., 500, -1., 1.,"#Delta x(cm)", "#Delta y(cm)","dN^{Events}/dxdy");
    AddHistogram1D<TH1D>("hVertexZBeforeVertexCut", "Z distribution of the vertex (before std. vertex cut)", "", 200, -20., 20., "#Delta z(cm)","dN^{Events}/dz");
    AddHistogram1D<TH1D>("hVertexZAfterVertexCut", "Z distribution of the vertex (after std. vertex cut)", "", 200, -20., 20., "#Delta z(cm)","dN^{Events}/dz");
    AddHistogram1D<TH1D>("hVertexR", "R distribution of the vertex", "", 100, 0., 1., "#Delta r(cm)","dN^{Events}/dr");
    AddHistogram1D<TH1D>("hCentralityV0M", "Centrality distribution V0M", "", fNumberOfCentralityBins, 0., 100., "Centrality","dN^{Events}");
    AddHistogram1D<TH1D>("hCentralityV0A", "Centrality distribution V0A", "", fNumberOfCentralityBins, 0., 100., "Centrality","dN^{Events}");
    AddHistogram1D<TH1D>("hCentralityV0C", "Centrality distribution V0C", "", fNumberOfCentralityBins, 0., 100., "Centrality","dN^{Events}");
    AddHistogram1D<TH1D>("hCentralityZNA", "Centrality distribution ZNA", "", fNumberOfCentralityBins, 0., 100., "Centrality","dN^{Events}");
    AddHistogram1D<TH1D>("hCentrality", Form("Centrality distribution %s", fCentralityType.Data()), "", fNumberOfCentralityBins, 0., 100., "Centrality","dN^{Events}");


    AddHistogram2D<TH2D>("hTrackCountAcc", "Number of tracks in acceptance vs. centrality", "LEGO2", 750, 0., 750., fNumberOfCentralityBins, 0, 100, "N tracks","Centrality", "dN^{Events}/dN^{Tracks}");
    AddHistogram2D<TH2D>("hTrackPt", "Tracks p_{T} distribution", "", 1000, 0., 250., fNumberOfCentralityBins, 0, 100, "p_{T} (GeV/c)", "Centrality", "dN^{Tracks}/dp_{T}");
    AddHistogram2D<TH2D>("hTrackPtNegEta", "Tracks p_{T} distribution (negative #eta)", "", 1000, 0., 250., fNumberOfCentralityBins, 0, 100, "p_{T} (GeV/c)","Centrality","dN^{Tracks}/dp_{T}");
    AddHistogram2D<TH2D>("hTrackPtPosEta", "Tracks p_{T} distribution (positive #eta)", "", 1000, 0., 250., fNumberOfCentralityBins, 0, 100, "p_{T} (GeV/c)","Centrality","dN^{Tracks}/dp_{T}");
    AddHistogram1D<TH1D>("hTrackCharge", "Charge", "", 11, -5, 5, "Charge (e)","dN^{Tracks}/dq");
    AddHistogram1D<TH1D>("hTrackPhi", "Track #phi distribution", "", 360, 0, TMath::TwoPi(), "#phi","dN^{Tracks}/d#phi");
    AddHistogram2D<TH2D>("hTrackPhiEta", "Track angular distribution", "LEGO2", 100, 0., 2*TMath::Pi(),100, -2.5, 2.5, "#phi","#eta","dN^{Tracks}/(d#phi d#eta)");

    AddHistogram2D<TH2D>("hTrackPhiPtCut", "Track #phi distribution for different pT cuts", "LEGO2", 360, 0, TMath::TwoPi(), 20, 0, 20, "#phi", "p_{T} lower cut", "dN^{Tracks}/d#phi dp_{T}");
    AddHistogram2D<TH2D>("hTrackPhiTrackType", "Track #phi distribution for different track types", "LEGO2", 360, 0, TMath::TwoPi(), 3, 0, 3, "#phi", "Label", "dN^{Tracks}/d#phi");
    AddHistogram2D<TH2D>("hTrackEta", "Track #eta distribution", "COLZ", 180, -fTrackEtaWindow, +fTrackEtaWindow, fNumberOfCentralityBins, 0., 100., "#eta", "Centrality", "dN^{Tracks}/d#eta");
    if (fAnalyzeJets)
    {
      // ######## Jet QA
      AddHistogram1D<TH1D>("hRawJetArea", "Jets area distribution w/o area cut", "", 200, 0., 2., "Area","dN^{Jets}/dA");
      AddHistogram2D<TH2D>("hJetArea", "Jets area distribution", "COLZ", 200, 0., 2., 150, 0.,150., "Area","Jet p_{T}","dN^{Jets}/dA");
      AddHistogram2D<TH2D>("hRawJetPhiEta", "Raw Jets angular distribution w/o #eta cut", "LEGO2", 360, 0., 2*TMath::Pi(),100, -1.0, 1.0, "#phi","#eta","dN^{Jets}/(d#phi d#eta)");
      AddHistogram2D<TH2D>("hJetEta", "Jets #eta distribution", "COLZ", 180, -fTrackEtaWindow, +fTrackEtaWindow, fNumberOfCentralityBins, 0., 100., "#eta", "Centrality", "dN^{Jets}/d#eta");
      AddHistogram2D<TH2D>("hJetEta2GeVTracks", "Jets #eta distribution, track p_{T} > 2 GeV", "COLZ", 180, -fTrackEtaWindow, +fTrackEtaWindow, fNumberOfCentralityBins, 0., 100., "#eta", "Centrality", "dN^{Jets}/d#eta");
      AddHistogram2D<TH2D>("hJetEta4GeVTracks", "Jets #eta distribution, track p_{T} > 4 GeV", "COLZ", 180, -fTrackEtaWindow, +fTrackEtaWindow, fNumberOfCentralityBins, 0., 100., "#eta", "Centrality", "dN^{Jets}/d#eta");
      AddHistogram2D<TH2D>("hJetPhiEta", "Jets angular distribution", "LEGO2", 360, 0., 2*TMath::Pi(),100, -1.0, 1.0, "#phi","#eta","dN^{Jets}/(d#phi d#eta)");
      AddHistogram2D<TH2D>("hJetPtVsConstituentCount", "Jets number of constituents vs. jet p_{T}", "COLZ", 400, 0., 200., 100, 0., 100., "p_{T}","N^{Tracks}","dN^{Jets}/(dp_{T} dN^{tracks})");
    }
  }

  // register Histograms
  for (Int_t i = 0; i < fHistCount; i++)
  {
    fOutputList->Add(fHistList->At(i));
  }
  
  PostData(1,fOutputList); // important for merging

}

//________________________________________________________________________
AliAnalysisTaskChargedJetsPA::AliAnalysisTaskChargedJetsPA(const char *name, const char* trackArrayName, const char* jetArrayName, const char* backgroundJetArrayName) : AliAnalysisTaskSE(name), fOutputList(0), fAnalyzeJets(1), fAnalyzeJetProfile(1), fAnalyzeQA(1), fAnalyzeBackground(1), fAnalyzeDeprecatedBackgrounds(1), fAnalyzePythia(0), fAnalyzeMassCorrelation(0), fHasTracks(0), fHasJets(0), fHasBackgroundJets(0), fIsKinematics(0), fUseDefaultVertexCut(1), fUsePileUpCut(1), fSetCentralityToOne(0), fNoExternalBackground(0), fPartialAnalysisNParts(1), fPartialAnalysisIndex(0), fJetArray(0), fTrackArray(0), fBackgroundJetArray(0), fJetArrayName(0), fTrackArrayName(0), fBackgroundJetArrayName(0), fNumPtHardBins(11), fUsePtHardBin(-1), fRhoTaskName(), fNcoll(6.88348), fRandConeRadius(0.4), fSignalJetRadius(0.4), fBackgroundJetRadius(0.4), fTRBackgroundConeRadius(0.6), fNumberRandCones(8), fNumberExcludedJets(-1), fDijetMaxAngleDeviation(10.0), fPhysicalJetRadius(0.6), fSignalJetEtaWindow(0.5), fBackgroundJetEtaWindow(0.5), fTrackEtaWindow(0.9), fMinTrackPt(0.150), fMinJetPt(0.15), fMinJetArea(0.5), fMinBackgroundJetPt(0.0), fMinDijetLeadingPt(10.0), fNumberOfCentralityBins(20), fCentralityType("V0A"), fFirstLeadingJet(0), fSecondLeadingJet(0), fNumberSignalJets(0), fCrossSection(0.0), fTrials(0.0), fRandom(0), fHelperClass(0), fInitialized(0), fTaskInstanceCounter(0), fHistList(0), fHistCount(0), fIsDEBUG(0), fEventCounter(0)
{
  #ifdef DEBUGMODE
    AliInfo("Calling constructor.");
  #endif

  // Every instance of this task gets his own number
  static Int_t instance = 0;
  fTaskInstanceCounter = instance;
  instance++;

  fTrackArrayName = new TString(trackArrayName);
  if (fTrackArrayName->Contains("MCParticles") || fTrackArrayName->Contains("mcparticles"))
    fIsKinematics = kTRUE;

  fJetArrayName = new TString(jetArrayName);
  if (strcmp(fJetArrayName->Data(),"") == 0)
  {
    fAnalyzeJets = kFALSE;
    fAnalyzeJetProfile = kFALSE;
  }
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
  Double_t realConeArea = (2.0*fTrackEtaWindow) * TMath::TwoPi() * MCGetOverlapCircleRectancle(eta, phi, radius, -fTrackEtaWindow, +fTrackEtaWindow, 0., TMath::TwoPi());
  tmpConePt -= background * realConeArea; // subtract background

  return tmpConePt;
}


//________________________________________________________________________
inline Double_t AliAnalysisTaskChargedJetsPA::GetPtHard()
{
  #ifdef DEBUGMODE
    AliInfo("Starting GetPtHard.");
  #endif
  AliGenPythiaEventHeader* pythiaHeader = dynamic_cast<AliGenPythiaEventHeader*>(MCEvent()->GenEventHeader());
  if (MCEvent()) 
    if (!pythiaHeader)
    {
      // Check if AOD
      AliAODMCHeader* aodMCH = dynamic_cast<AliAODMCHeader*>(InputEvent()->FindListObject(AliAODMCHeader::StdBranchName()));

      if (aodMCH)
      {
        for(UInt_t i = 0;i<aodMCH->GetNCocktailHeaders();i++)
        {
          pythiaHeader = dynamic_cast<AliGenPythiaEventHeader*>(aodMCH->GetCocktailHeader(i));
          if (pythiaHeader) break;
        }
      }
    }

  #ifdef DEBUGMODE
    AliInfo("Ending GetPtHard.");
  #endif
  if (pythiaHeader)
    return pythiaHeader->GetPtHard();

  AliWarning(Form("In task %s: GetPtHard() failed!", GetName()));
  return -1.0;
}


//________________________________________________________________________
inline Double_t AliAnalysisTaskChargedJetsPA::GetPythiaTrials()
{
  #ifdef DEBUGMODE
    AliInfo("Starting GetPythiaTrials.");
  #endif
  AliGenPythiaEventHeader* pythiaHeader = dynamic_cast<AliGenPythiaEventHeader*>(MCEvent()->GenEventHeader());
  if (MCEvent()) 
    if (!pythiaHeader)
    {
      // Check if AOD
      AliAODMCHeader* aodMCH = dynamic_cast<AliAODMCHeader*>(InputEvent()->FindListObject(AliAODMCHeader::StdBranchName()));

      if (aodMCH)
      {
        for(UInt_t i = 0;i<aodMCH->GetNCocktailHeaders();i++)
        {
          pythiaHeader = dynamic_cast<AliGenPythiaEventHeader*>(aodMCH->GetCocktailHeader(i));
          if (pythiaHeader) break;
        }
      }
    }

  #ifdef DEBUGMODE
    AliInfo("Ending GetPythiaTrials.");
  #endif
  if (pythiaHeader)
    return pythiaHeader->Trials();

  AliWarning(Form("In task %s: GetPythiaTrials() failed!", GetName()));
  return -1.0;
}



//________________________________________________________________________
inline Int_t AliAnalysisTaskChargedJetsPA::GetPtHardBin()
{
  #ifdef DEBUGMODE
    AliInfo("Starting GetPtHardBin.");
  #endif
  // ########## PT HARD BIN EDGES
  const Int_t kPtHardLowerEdges[] =  { 0, 5,11,21,36,57, 84,117,152,191,234};
  const Int_t kPtHardHigherEdges[] = { 5,11,21,36,57,84,117,152,191,234,1000000};

  Int_t tmpPtHardBin = 0;
  Double_t tmpPtHard = GetPtHard();
 
  for (tmpPtHardBin = 0; tmpPtHardBin <= fNumPtHardBins; tmpPtHardBin++)
    if (tmpPtHard >= kPtHardLowerEdges[tmpPtHardBin] && tmpPtHard < kPtHardHigherEdges[tmpPtHardBin])
      break;

  #ifdef DEBUGMODE
    AliInfo("Ending GetPtHardBin.");
  #endif
  return tmpPtHardBin;
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
inline Bool_t AliAnalysisTaskChargedJetsPA::IsEventInAcceptance(AliVEvent* event)
{
  if (!event)
    return kFALSE;

  FillHistogram("hEventAcceptance", 0.5); // number of events before manual cuts
  if(fUsePileUpCut)
    if(fHelperClass->IsPileUpEvent(event))
      return kFALSE;

  FillHistogram("hEventAcceptance", 1.5); // number of events after pileup cuts

  if(fAnalyzeQA)
    FillHistogram("hVertexZBeforeVertexCut",event->GetPrimaryVertex()->GetZ());

  if(fUseDefaultVertexCut)
  {
    if(!fHelperClass->IsVertexSelected2013pA(event))
      return kFALSE;
  }
  else // Failsafe vertex cut
  {
    if(TMath::Abs(event->GetPrimaryVertex()->GetZ()) > 10.0)
      return kFALSE;
  }

  if(fAnalyzeQA)
    FillHistogram("hVertexZAfterVertexCut",event->GetPrimaryVertex()->GetZ());

  FillHistogram("hEventAcceptance", 2.5); // number of events after vertex cut

  return kTRUE;
}

//________________________________________________________________________
inline Bool_t AliAnalysisTaskChargedJetsPA::IsTrackInAcceptance(AliVParticle* track)
{
  FillHistogram("hTrackAcceptance", 0.5);
  if (track != 0)
  {
    if(fIsKinematics)
    {
      // TODO: Only working for AOD MC
      if((!track->Charge()) || (!(static_cast<AliAODMCParticle*>(track))->IsPhysicalPrimary()) )
        return kFALSE;
    }
    if (TMath::Abs(track->Eta()) <= fTrackEtaWindow)
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
    if (TMath::Abs(jet->Eta()) <= fBackgroundJetEtaWindow)
      if (jet->Pt() >= fMinBackgroundJetPt)
        return kTRUE;

  return kFALSE;
}

//________________________________________________________________________
inline Bool_t AliAnalysisTaskChargedJetsPA::IsSignalJetInAcceptance(AliEmcalJet *jet)
{   
  FillHistogram("hJetAcceptance", 0.5);
  if (jet != 0)
    if (TMath::Abs(jet->Eta()) <= fSignalJetEtaWindow)
    {
      FillHistogram("hJetAcceptance", 1.5);
      if (jet->Pt() >= fMinJetPt)
      {
        FillHistogram("hJetAcceptance", 2.5);
        if (jet->Area() >= fMinJetArea)
        {
          FillHistogram("hJetAcceptance", 3.5);
          return kTRUE;
        }
      }
    }
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

  // Initialize helper class (for vertex selection & pile up correction)
  fHelperClass = new AliAnalysisUtils();
  fHelperClass->SetCutOnZVertexSPD(kFALSE);
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
    AliError("Jet pointer passed to GetCorrectedJetPt() not valid!");
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
Double_t AliAnalysisTaskChargedJetsPA::GetDeltaPt(Double_t rho, Double_t leadingJetExclusionProbability)
{
  #ifdef DEBUGMODE
    AliInfo("Getting Delta Pt.");
  #endif

  // Define an invalid delta pt
  Double_t deltaPt = -10000.0;

  // Define eta range
  Double_t etaMin, etaMax;
  etaMin = -(fTrackEtaWindow-fRandConeRadius);
  etaMax = +(fTrackEtaWindow-fRandConeRadius);

  // Define random cone
  Bool_t coneValid = kTRUE;
  Double_t tmpRandConeEta = etaMin + fRandom->Rndm()*(etaMax-etaMin);
  Double_t tmpRandConePhi = fRandom->Rndm()*TMath::TwoPi();

  // if there is a jet, check for overlap if demanded
  if(leadingJetExclusionProbability)
  {
    AliEmcalJet* tmpLeading = dynamic_cast<AliEmcalJet*>(fJetArray->At(0));
    // Get leading jet (regardless of pT)
    for (Int_t i = 1; i<fJetArray->GetEntries(); i++)
    {
      AliEmcalJet* tmpJet = static_cast<AliEmcalJet*>(fJetArray->At(i));
      // if jet is in acceptance and higher, take as new leading
      if (tmpJet)
        if ((TMath::Abs(tmpJet->Eta()) <= fSignalJetEtaWindow) && (tmpJet->Area() >= fMinJetArea))
          if((!tmpLeading) || (tmpJet->Pt() > tmpLeading->Pt()))
            tmpLeading = tmpJet;
    }
    if(tmpLeading)
    {
      Double_t excludedJetPhi = tmpLeading->Phi();
      Double_t excludedJetEta = tmpLeading->Eta();
      Double_t tmpDeltaPhi = GetDeltaPhi(tmpRandConePhi, excludedJetPhi);

      // Check, if cone has overlap with jet
      if ( tmpDeltaPhi*tmpDeltaPhi + TMath::Abs(tmpRandConeEta-excludedJetEta)*TMath::Abs(tmpRandConeEta-excludedJetEta) <= fRandConeRadius*fRandConeRadius)
      {
        // Define probability to exclude the RC
        Double_t probability = leadingJetExclusionProbability;

        // Only exclude cone with a given probability
        if (fRandom->Rndm()<=probability)
          coneValid = kFALSE;
      }
    }
  }


  // Get the cones' pt and calculate delta pt
  if (coneValid)
    deltaPt = GetConePt(tmpRandConeEta,tmpRandConePhi,fRandConeRadius) - (rho*fRandConeRadius*fRandConeRadius*TMath::Pi());

  return deltaPt;
  #ifdef DEBUGMODE
    AliInfo("Got Delta Pt.");
  #endif
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

    tmpSummedArea += backgroundJet->Area();
    if(backgroundJet->Pt() > 0.150)
      tmpCoveredArea += backgroundJet->Area();

    if (!IsBackgroundJetInAcceptance(backgroundJet))
      continue;

    // Search for overlap with signal jets
    Bool_t isOverlapping = kFALSE;
    for(Int_t j=0;j<numberExcludeLeadingJets;j++)
    {
      AliEmcalJet* signalJet = fSignalJets[j];
     
      if(signalJet->Pt() >= 5.0)
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
    if ((i != leadingKTJets[0]) && (i != leadingKTJets[1]))
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
      if(!isOverlapping)
      {
        tmpRhoImprovedCMS[rhoImprovedCMSJetCount] = tmpRho;
        rhoImprovedCMSJetCount++;
      }

      // PbPb w/o ghosts approach (just neglect ghosts)
      if ((i != leadingKTJets[0]) && (i != leadingKTJets[1]))
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
void AliAnalysisTaskChargedJetsPA::GetKTBackgroundDensity(Int_t numberExcludeLeadingJets, Double_t& rhoImprovedCMS)
{
  #ifdef DEBUGMODE
    AliInfo("Getting KT background density.");
  #endif

  static Double_t tmpRhoImprovedCMS[1024];
  Double_t tmpCoveredArea = 0.0;
  Double_t tmpSummedArea = 0.0;

  // Setting invalid values
  rhoImprovedCMS = 0.0;

  Int_t rhoImprovedCMSJetCount = 0;

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
      if(signalJet->Pt() >= 5.0)     
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

    Double_t tmpRho = backgroundJet->Pt() / backgroundJet->Area();

    if(backgroundJet->Pt() > 0.150)
      if(!isOverlapping)
      {
        tmpRhoImprovedCMS[rhoImprovedCMSJetCount] = tmpRho;
        rhoImprovedCMSJetCount++;
      }
  }

  if (rhoImprovedCMSJetCount > 0)
  {
    rhoImprovedCMS = TMath::Median(rhoImprovedCMSJetCount, tmpRhoImprovedCMS) * tmpCoveredArea/tmpSummedArea;
  }
  #ifdef DEBUGMODE
    AliInfo("Got KT background density.");
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
    numberExcludeLeadingJets = fNumberSignalJets;
  if (fNumberSignalJets < numberExcludeLeadingJets)
    numberExcludeLeadingJets = fNumberSignalJets;

  Int_t fSignalJetCount5GeV = 0;
  for(Int_t j=0;j<numberExcludeLeadingJets;j++)
  {
    AliEmcalJet* signalJet = fSignalJets[j];
    if(signalJet->Pt() < 5.0)
      continue;
    fSignalJetCount5GeV++;
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
        AliEmcalJet* signalJet = fSignalJets[j];

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

  Double_t tmpFullTPCArea = (2.0*fTrackEtaWindow) * TMath::TwoPi();
  Double_t tmpAreaCone02     = tmpFullTPCArea;
  Double_t tmpAreaCone04     = tmpFullTPCArea;
  Double_t tmpAreaCone06     = tmpFullTPCArea;
  Double_t tmpAreaCone08     = tmpFullTPCArea;
  Double_t tmpAreaWithinJets = tmpFullTPCArea;
  std::vector<Double_t> tmpEtas(fSignalJetCount5GeV);
  std::vector<Double_t> tmpPhis(fSignalJetCount5GeV);

  Int_t iSignal = 0;
  for(Int_t i=0;i<numberExcludeLeadingJets;i++)
  {
    AliEmcalJet* tmpJet = fSignalJets[i];

    if(tmpJet->Pt() < 5.0)
      continue;

    tmpEtas[iSignal] = tmpJet->Eta();
    tmpPhis[iSignal] = tmpJet->Phi();
    tmpAreaWithinJets -= tmpJet->Area();

    iSignal++;
  }

  tmpAreaCone02 -= tmpFullTPCArea * MCGetOverlapMultipleCirclesRectancle(fSignalJetCount5GeV, tmpEtas, tmpPhis, 0.2, -fTrackEtaWindow, +fTrackEtaWindow, 0., TMath::TwoPi());
  tmpAreaCone04 -= tmpFullTPCArea * MCGetOverlapMultipleCirclesRectancle(fSignalJetCount5GeV, tmpEtas, tmpPhis, 0.4, -fTrackEtaWindow, +fTrackEtaWindow, 0., TMath::TwoPi());
  tmpAreaCone06 -= tmpFullTPCArea * MCGetOverlapMultipleCirclesRectancle(fSignalJetCount5GeV, tmpEtas, tmpPhis, 0.6, -fTrackEtaWindow, +fTrackEtaWindow, 0., TMath::TwoPi());
  tmpAreaCone08 -= tmpFullTPCArea * MCGetOverlapMultipleCirclesRectancle(fSignalJetCount5GeV, tmpEtas, tmpPhis, 0.8, -fTrackEtaWindow, +fTrackEtaWindow, 0., TMath::TwoPi());
 
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
void AliAnalysisTaskChargedJetsPA::GetPPBackgroundDensity(Double_t& background, AliEmcalJet* jet)
{
  // This is the background that was used for the pp 7 TeV ALICE paper
  Double_t jetMom[3] = { jet->Px(), jet->Py(), jet->Pz() };
  TVector3 jet3mom1(jetMom);
  TVector3 jet3mom2(jetMom);
  background = 0;

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

  // Check, if analysis should be done in pt hard bins
  if(fUsePtHardBin != -1)
    if(GetPtHardBin() != fUsePtHardBin)
      return;

  // This is to take only every Nth event
  if((fEventCounter+fPartialAnalysisIndex) % fPartialAnalysisNParts != 0)
    return;

  FillHistogram("hNumberEvents",0.5);

  if(!IsEventInAcceptance(event))
    return;

  FillHistogram("hNumberEvents",1.5);

  #ifdef DEBUGMODE
    AliInfo("Calculate()::Init done.");
  #endif

  ////////////////////// NOTE: Get Centrality, (Leading)Signal jets and Background

  // Get centrality
  AliCentrality* tmpCentrality = NULL;
  tmpCentrality = event->GetCentrality();
  Double_t centralityPercentile = -1.0;
  Double_t centralityPercentileV0A = 0.0;
  Double_t centralityPercentileV0C = 0.0;
  Double_t centralityPercentileV0M = 0.0;
  Double_t centralityPercentileZNA = 0.0;
  if (tmpCentrality != NULL)
  {
    centralityPercentile = tmpCentrality->GetCentralityPercentile(fCentralityType.Data());
    centralityPercentileV0A = tmpCentrality->GetCentralityPercentile("V0A");
    centralityPercentileV0C = tmpCentrality->GetCentralityPercentile("V0C");
    centralityPercentileV0M = tmpCentrality->GetCentralityPercentile("V0M");
    centralityPercentileZNA = tmpCentrality->GetCentralityPercentile("ZNA");
  }

  if((centralityPercentile < 0.0) || (centralityPercentile > 100.0))
    AliWarning(Form("Centrality value not valid (c=%E)",centralityPercentile)); 

  if(fSetCentralityToOne)
    centralityPercentile = 1.0;

  // Get jets
  if (fAnalyzeBackground || fAnalyzeJets)
    GetSignalJets();

  // Get background estimates
  Double_t              backgroundKTImprovedCMS = -1.0;
  Double_t              backgroundKTImprovedCMSExternal = -1.0;
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
  // Calculate background for different jet exclusions

  if (fAnalyzeBackground)
  {

    if(fAnalyzeDeprecatedBackgrounds)
      GetKTBackgroundDensityAll    (fNumberExcludedJets, backgroundKTPbPb, backgroundKTPbPbWithGhosts, backgroundKTCMS, backgroundKTImprovedCMS, backgroundKTMean, backgroundKTTrackLike);
    else
      GetKTBackgroundDensity       (fNumberExcludedJets, backgroundKTImprovedCMS);

    if(fAnalyzeDeprecatedBackgrounds)
      GetTRBackgroundDensity    (fNumberExcludedJets, backgroundTRNoExcl, backgroundTRCone02, backgroundTRCone04, backgroundTRCone06, backgroundTRCone08, backgroundTRExact);

    backgroundKTImprovedCMSExternal = GetExternalRho();
    if(fNoExternalBackground)
      backgroundKTImprovedCMSExternal = 0;

  }

  #ifdef DEBUGMODE
    AliInfo("Calculate()::Centrality&SignalJets&Background-Calculation done.");
  #endif

  if (fAnalyzeQA)
  {
    FillHistogram("hVertexX",event->GetPrimaryVertex()->GetX());
    FillHistogram("hVertexY",event->GetPrimaryVertex()->GetY());
    FillHistogram("hVertexXY",event->GetPrimaryVertex()->GetX(), event->GetPrimaryVertex()->GetY());
    FillHistogram("hVertexR",TMath::Sqrt(event->GetPrimaryVertex()->GetX()*event->GetPrimaryVertex()->GetX() + event->GetPrimaryVertex()->GetY()*event->GetPrimaryVertex()->GetY()));
    FillHistogram("hCentralityV0M",centralityPercentileV0M);
    FillHistogram("hCentralityV0A",centralityPercentileV0A);
    FillHistogram("hCentralityV0C",centralityPercentileV0C);
    FillHistogram("hCentralityZNA",centralityPercentileZNA);
    FillHistogram("hCentrality",centralityPercentile);

    Int_t trackCountAcc = 0;
    Int_t nTracks = fTrackArray->GetEntries();
    for (Int_t i = 0; i < nTracks; i++)
    {
      AliVTrack* track = static_cast<AliVTrack*>(fTrackArray->At(i));

      if (track != 0)
        if (track->Pt() >= fMinTrackPt)
          FillHistogram("hTrackPhiEta", track->Phi(),track->Eta(), 1);

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
          FillHistogram("hTrackPhiTrackType", track->Phi(), (static_cast<AliPicoTrack*>(track))->GetTrackType());

        for(Int_t j=0;j<20;j++)
          if(track->Pt() > j)
            FillHistogram("hTrackPhiPtCut", track->Phi(), track->Pt());

        FillHistogram("hTrackCharge", track->Charge());
        trackCountAcc++;
      }
    }
    FillHistogram("hTrackCountAcc", trackCountAcc, centralityPercentile);


/*
    // This is code that is run locally to get some special events
    TFile* fileOutput = new TFile("SpecialEvents.root", "UPDATE");
    if(fSecondLeadingJet&&(fSecondLeadingJet->Pt()>10.))
    {
      cout << "Event found\n";
      TH2D* tmpEvent = new TH2D(Form("Event%lu", fEventCounter), "", 40, -0.9, 0.9, 40, 0., TMath::TwoPi());
      tmpEvent->GetXaxis()->SetTitle("#eta");
      tmpEvent->GetYaxis()->SetTitle("#phi");
      tmpEvent->SetOption("LEGO2");
      tmpEvent->Sumw2();

      for (Int_t i = 0; i < nTracks; i++)
      {
        AliVTrack* track = static_cast<AliVTrack*>(fTrackArray->At(i));

        if (IsTrackInAcceptance(track))
        {
          tmpEvent->Fill(track->Eta(), track->Phi(), track->Pt());
        }
      }
      tmpEvent->Write(0, kOverwrite);
    }
    fileOutput->Close();
*/
  }
  #ifdef DEBUGMODE
    AliInfo("Calculate()::QA done.");
  #endif

  ////////////////////// NOTE: Jet analysis and calculations

  if (fAnalyzeJets)
  {
    for (Int_t i = 0; i<fJetArray->GetEntries(); i++)
    {
      AliEmcalJet* tmpJet = static_cast<AliEmcalJet*>(fJetArray->At(i));
      if (!tmpJet)
        continue;

      FillHistogram("hRawJetPt", tmpJet->Pt());
      if (tmpJet->Pt() >= 5.0)
      {
      // ### RAW JET ANALYSIS
        if (tmpJet->Area() >= fMinJetArea)
          FillHistogram("hRawJetPhiEta", tmpJet->Phi(), tmpJet->Eta());
        if (TMath::Abs(tmpJet->Eta()) <= fSignalJetEtaWindow)
          FillHistogram("hRawJetArea", tmpJet->Area());
      }

      if(IsSignalJetInAcceptance(tmpJet))
      {
      // ### SIGNAL JET ANALYSIS
        // Jet spectra
        FillHistogram("hJetPt", tmpJet->Pt(), centralityPercentile);
        FillHistogram("hJetPtBgrdSubtractedKTImprovedCMS", GetCorrectedJetPt(tmpJet, backgroundKTImprovedCMS), centralityPercentile);
        if(tmpJet->Phi() >= TMath::Pi())
          FillHistogram("hJetPtBgrdSubtractedKTImprovedCMS_Phi2", GetCorrectedJetPt(tmpJet, backgroundKTImprovedCMS), centralityPercentile);
        else          
          FillHistogram("hJetPtBgrdSubtractedKTImprovedCMS_Phi1", GetCorrectedJetPt(tmpJet, backgroundKTImprovedCMS), centralityPercentile);

        // pp background
        GetPPBackgroundDensity(backgroundPP, tmpJet);
        FillHistogram("hJetPtBgrdSubtractedPP", GetCorrectedJetPt(tmpJet, backgroundPP), centralityPercentile);

        FillHistogram("hJetPtBgrdSubtractedExternal", GetCorrectedJetPt(tmpJet, backgroundKTImprovedCMSExternal), centralityPercentile);
        FillHistogram("hJetPtSubtractedRhoKTImprovedCMS", tmpJet->Pt(), centralityPercentile, tmpJet->Pt() - GetCorrectedJetPt(tmpJet, backgroundKTImprovedCMS));
        FillHistogram("hJetPtSubtractedRhoExternal", tmpJet->Pt(), centralityPercentile, tmpJet->Pt() - GetCorrectedJetPt(tmpJet, backgroundKTImprovedCMSExternal));
        FillHistogram("hJetPtSubtractedRhoPP", tmpJet->Pt(), centralityPercentile, tmpJet->Pt() - GetCorrectedJetPt(tmpJet, backgroundPP));

        if(fAnalyzeDeprecatedBackgrounds)
        {
          FillHistogram("hJetPtBgrdSubtractedTR", GetCorrectedJetPt(tmpJet, backgroundTRCone06), centralityPercentile);
          FillHistogram("hJetPtBgrdSubtractedKTPbPb", GetCorrectedJetPt(tmpJet, backgroundKTPbPb), centralityPercentile);
          FillHistogram("hJetPtBgrdSubtractedKTPbPbWithGhosts", GetCorrectedJetPt(tmpJet, backgroundKTPbPbWithGhosts), centralityPercentile);
          FillHistogram("hJetPtBgrdSubtractedKTCMS", GetCorrectedJetPt(tmpJet, backgroundKTCMS), centralityPercentile);
          FillHistogram("hJetPtBgrdSubtractedKTMean", GetCorrectedJetPt(tmpJet, backgroundKTMean), centralityPercentile);
          FillHistogram("hJetPtBgrdSubtractedKTTrackLike", GetCorrectedJetPt(tmpJet, backgroundKTTrackLike), centralityPercentile);
        }

        // Jet profile analysis
        if(TMath::Abs(tmpJet->Eta()) <= 0.3)
        {
          if(tmpJet->Pt()>=70.0)
          {
            FillHistogram("hJetProfile70GeV", 0.05-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.05, backgroundKTImprovedCMSExternal))/GetCorrectedJetPt(tmpJet, backgroundKTImprovedCMSExternal));
            FillHistogram("hJetProfile70GeV", 0.10-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.10, backgroundKTImprovedCMSExternal))/GetCorrectedJetPt(tmpJet, backgroundKTImprovedCMSExternal));
            FillHistogram("hJetProfile70GeV", 0.15-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.15, backgroundKTImprovedCMSExternal))/GetCorrectedJetPt(tmpJet, backgroundKTImprovedCMSExternal));
            FillHistogram("hJetProfile70GeV", 0.20-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.20, backgroundKTImprovedCMSExternal))/GetCorrectedJetPt(tmpJet, backgroundKTImprovedCMSExternal));
            FillHistogram("hJetProfile70GeV", 0.25-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.25, backgroundKTImprovedCMSExternal))/GetCorrectedJetPt(tmpJet, backgroundKTImprovedCMSExternal));
            FillHistogram("hJetProfile70GeV", 0.30-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.30, backgroundKTImprovedCMSExternal))/GetCorrectedJetPt(tmpJet, backgroundKTImprovedCMSExternal));
            FillHistogram("hJetProfile70GeV", 0.35-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.35, backgroundKTImprovedCMSExternal))/GetCorrectedJetPt(tmpJet, backgroundKTImprovedCMSExternal));
            FillHistogram("hJetProfile70GeV", 0.40-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.40, backgroundKTImprovedCMSExternal))/GetCorrectedJetPt(tmpJet, backgroundKTImprovedCMSExternal));
            FillHistogram("hJetProfile70GeV", 0.45-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.45, backgroundKTImprovedCMSExternal))/GetCorrectedJetPt(tmpJet, backgroundKTImprovedCMSExternal));
            FillHistogram("hJetProfile70GeV", 0.50-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.50, backgroundKTImprovedCMSExternal))/GetCorrectedJetPt(tmpJet, backgroundKTImprovedCMSExternal));
            FillHistogram("hJetProfile70GeV", 0.55-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.55, backgroundKTImprovedCMSExternal))/GetCorrectedJetPt(tmpJet, backgroundKTImprovedCMSExternal));
            FillHistogram("hJetProfile70GeV", 0.60-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.60, backgroundKTImprovedCMSExternal))/GetCorrectedJetPt(tmpJet, backgroundKTImprovedCMSExternal));
          }
          else if(GetCorrectedJetPt(tmpJet, backgroundKTImprovedCMSExternal)>=60.0)
          {
            FillHistogram("hJetProfile60GeV", 0.05-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.05, backgroundKTImprovedCMSExternal))/GetCorrectedJetPt(tmpJet, backgroundKTImprovedCMSExternal));
            FillHistogram("hJetProfile60GeV", 0.10-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.10, backgroundKTImprovedCMSExternal))/GetCorrectedJetPt(tmpJet, backgroundKTImprovedCMSExternal));
            FillHistogram("hJetProfile60GeV", 0.15-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.15, backgroundKTImprovedCMSExternal))/GetCorrectedJetPt(tmpJet, backgroundKTImprovedCMSExternal));
            FillHistogram("hJetProfile60GeV", 0.20-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.20, backgroundKTImprovedCMSExternal))/GetCorrectedJetPt(tmpJet, backgroundKTImprovedCMSExternal));
            FillHistogram("hJetProfile60GeV", 0.25-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.25, backgroundKTImprovedCMSExternal))/GetCorrectedJetPt(tmpJet, backgroundKTImprovedCMSExternal));
            FillHistogram("hJetProfile60GeV", 0.30-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.30, backgroundKTImprovedCMSExternal))/GetCorrectedJetPt(tmpJet, backgroundKTImprovedCMSExternal));
            FillHistogram("hJetProfile60GeV", 0.35-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.35, backgroundKTImprovedCMSExternal))/GetCorrectedJetPt(tmpJet, backgroundKTImprovedCMSExternal));
            FillHistogram("hJetProfile60GeV", 0.40-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.40, backgroundKTImprovedCMSExternal))/GetCorrectedJetPt(tmpJet, backgroundKTImprovedCMSExternal));
            FillHistogram("hJetProfile60GeV", 0.45-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.45, backgroundKTImprovedCMSExternal))/GetCorrectedJetPt(tmpJet, backgroundKTImprovedCMSExternal));
            FillHistogram("hJetProfile60GeV", 0.50-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.50, backgroundKTImprovedCMSExternal))/GetCorrectedJetPt(tmpJet, backgroundKTImprovedCMSExternal));
            FillHistogram("hJetProfile60GeV", 0.55-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.55, backgroundKTImprovedCMSExternal))/GetCorrectedJetPt(tmpJet, backgroundKTImprovedCMSExternal));
            FillHistogram("hJetProfile60GeV", 0.60-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.60, backgroundKTImprovedCMSExternal))/GetCorrectedJetPt(tmpJet, backgroundKTImprovedCMSExternal));
          }
          else if(GetCorrectedJetPt(tmpJet, backgroundKTImprovedCMSExternal)>=50.0)
          {
            FillHistogram("hJetProfile50GeV", 0.05-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.05, backgroundKTImprovedCMSExternal))/GetCorrectedJetPt(tmpJet, backgroundKTImprovedCMSExternal));
            FillHistogram("hJetProfile50GeV", 0.10-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.10, backgroundKTImprovedCMSExternal))/GetCorrectedJetPt(tmpJet, backgroundKTImprovedCMSExternal));
            FillHistogram("hJetProfile50GeV", 0.15-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.15, backgroundKTImprovedCMSExternal))/GetCorrectedJetPt(tmpJet, backgroundKTImprovedCMSExternal));
            FillHistogram("hJetProfile50GeV", 0.20-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.20, backgroundKTImprovedCMSExternal))/GetCorrectedJetPt(tmpJet, backgroundKTImprovedCMSExternal));
            FillHistogram("hJetProfile50GeV", 0.25-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.25, backgroundKTImprovedCMSExternal))/GetCorrectedJetPt(tmpJet, backgroundKTImprovedCMSExternal));
            FillHistogram("hJetProfile50GeV", 0.30-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.30, backgroundKTImprovedCMSExternal))/GetCorrectedJetPt(tmpJet, backgroundKTImprovedCMSExternal));
            FillHistogram("hJetProfile50GeV", 0.35-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.35, backgroundKTImprovedCMSExternal))/GetCorrectedJetPt(tmpJet, backgroundKTImprovedCMSExternal));
            FillHistogram("hJetProfile50GeV", 0.40-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.40, backgroundKTImprovedCMSExternal))/GetCorrectedJetPt(tmpJet, backgroundKTImprovedCMSExternal));
            FillHistogram("hJetProfile50GeV", 0.45-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.45, backgroundKTImprovedCMSExternal))/GetCorrectedJetPt(tmpJet, backgroundKTImprovedCMSExternal));
            FillHistogram("hJetProfile50GeV", 0.50-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.50, backgroundKTImprovedCMSExternal))/GetCorrectedJetPt(tmpJet, backgroundKTImprovedCMSExternal));
            FillHistogram("hJetProfile50GeV", 0.55-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.55, backgroundKTImprovedCMSExternal))/GetCorrectedJetPt(tmpJet, backgroundKTImprovedCMSExternal));
            FillHistogram("hJetProfile50GeV", 0.60-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.60, backgroundKTImprovedCMSExternal))/GetCorrectedJetPt(tmpJet, backgroundKTImprovedCMSExternal));
          }
          else if(GetCorrectedJetPt(tmpJet, backgroundKTImprovedCMSExternal)>=40.0)
          {
            FillHistogram("hJetProfile40GeV", 0.05-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.05, backgroundKTImprovedCMSExternal))/GetCorrectedJetPt(tmpJet, backgroundKTImprovedCMSExternal));
            FillHistogram("hJetProfile40GeV", 0.10-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.10, backgroundKTImprovedCMSExternal))/GetCorrectedJetPt(tmpJet, backgroundKTImprovedCMSExternal));
            FillHistogram("hJetProfile40GeV", 0.15-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.15, backgroundKTImprovedCMSExternal))/GetCorrectedJetPt(tmpJet, backgroundKTImprovedCMSExternal));
            FillHistogram("hJetProfile40GeV", 0.20-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.20, backgroundKTImprovedCMSExternal))/GetCorrectedJetPt(tmpJet, backgroundKTImprovedCMSExternal));
            FillHistogram("hJetProfile40GeV", 0.25-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.25, backgroundKTImprovedCMSExternal))/GetCorrectedJetPt(tmpJet, backgroundKTImprovedCMSExternal));
            FillHistogram("hJetProfile40GeV", 0.30-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.30, backgroundKTImprovedCMSExternal))/GetCorrectedJetPt(tmpJet, backgroundKTImprovedCMSExternal));
            FillHistogram("hJetProfile40GeV", 0.35-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.35, backgroundKTImprovedCMSExternal))/GetCorrectedJetPt(tmpJet, backgroundKTImprovedCMSExternal));
            FillHistogram("hJetProfile40GeV", 0.40-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.40, backgroundKTImprovedCMSExternal))/GetCorrectedJetPt(tmpJet, backgroundKTImprovedCMSExternal));
            FillHistogram("hJetProfile40GeV", 0.45-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.45, backgroundKTImprovedCMSExternal))/GetCorrectedJetPt(tmpJet, backgroundKTImprovedCMSExternal));
            FillHistogram("hJetProfile40GeV", 0.50-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.50, backgroundKTImprovedCMSExternal))/GetCorrectedJetPt(tmpJet, backgroundKTImprovedCMSExternal));
            FillHistogram("hJetProfile40GeV", 0.55-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.55, backgroundKTImprovedCMSExternal))/GetCorrectedJetPt(tmpJet, backgroundKTImprovedCMSExternal));
            FillHistogram("hJetProfile40GeV", 0.60-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.60, backgroundKTImprovedCMSExternal))/GetCorrectedJetPt(tmpJet, backgroundKTImprovedCMSExternal));
          }
          else if(GetCorrectedJetPt(tmpJet, backgroundKTImprovedCMSExternal)>=30.0)
          {
            FillHistogram("hJetProfile30GeV", 0.05-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.05, backgroundKTImprovedCMSExternal))/GetCorrectedJetPt(tmpJet, backgroundKTImprovedCMSExternal));
            FillHistogram("hJetProfile30GeV", 0.10-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.10, backgroundKTImprovedCMSExternal))/GetCorrectedJetPt(tmpJet, backgroundKTImprovedCMSExternal));
            FillHistogram("hJetProfile30GeV", 0.15-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.15, backgroundKTImprovedCMSExternal))/GetCorrectedJetPt(tmpJet, backgroundKTImprovedCMSExternal));
            FillHistogram("hJetProfile30GeV", 0.20-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.20, backgroundKTImprovedCMSExternal))/GetCorrectedJetPt(tmpJet, backgroundKTImprovedCMSExternal));
            FillHistogram("hJetProfile30GeV", 0.25-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.25, backgroundKTImprovedCMSExternal))/GetCorrectedJetPt(tmpJet, backgroundKTImprovedCMSExternal));
            FillHistogram("hJetProfile30GeV", 0.30-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.30, backgroundKTImprovedCMSExternal))/GetCorrectedJetPt(tmpJet, backgroundKTImprovedCMSExternal));
            FillHistogram("hJetProfile30GeV", 0.35-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.35, backgroundKTImprovedCMSExternal))/GetCorrectedJetPt(tmpJet, backgroundKTImprovedCMSExternal));
            FillHistogram("hJetProfile30GeV", 0.40-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.40, backgroundKTImprovedCMSExternal))/GetCorrectedJetPt(tmpJet, backgroundKTImprovedCMSExternal));
            FillHistogram("hJetProfile30GeV", 0.45-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.45, backgroundKTImprovedCMSExternal))/GetCorrectedJetPt(tmpJet, backgroundKTImprovedCMSExternal));
            FillHistogram("hJetProfile30GeV", 0.50-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.50, backgroundKTImprovedCMSExternal))/GetCorrectedJetPt(tmpJet, backgroundKTImprovedCMSExternal));
            FillHistogram("hJetProfile30GeV", 0.55-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.55, backgroundKTImprovedCMSExternal))/GetCorrectedJetPt(tmpJet, backgroundKTImprovedCMSExternal));
            FillHistogram("hJetProfile30GeV", 0.60-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.60, backgroundKTImprovedCMSExternal))/GetCorrectedJetPt(tmpJet, backgroundKTImprovedCMSExternal));
          }
          else if(GetCorrectedJetPt(tmpJet, backgroundKTImprovedCMSExternal)>=20.0)
          {
            FillHistogram("hJetProfile20GeV", 0.05-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.05, backgroundKTImprovedCMSExternal))/GetCorrectedJetPt(tmpJet, backgroundKTImprovedCMSExternal));
            FillHistogram("hJetProfile20GeV", 0.10-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.10, backgroundKTImprovedCMSExternal))/GetCorrectedJetPt(tmpJet, backgroundKTImprovedCMSExternal));
            FillHistogram("hJetProfile20GeV", 0.15-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.15, backgroundKTImprovedCMSExternal))/GetCorrectedJetPt(tmpJet, backgroundKTImprovedCMSExternal));
            FillHistogram("hJetProfile20GeV", 0.20-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.20, backgroundKTImprovedCMSExternal))/GetCorrectedJetPt(tmpJet, backgroundKTImprovedCMSExternal));
            FillHistogram("hJetProfile20GeV", 0.25-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.25, backgroundKTImprovedCMSExternal))/GetCorrectedJetPt(tmpJet, backgroundKTImprovedCMSExternal));
            FillHistogram("hJetProfile20GeV", 0.30-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.30, backgroundKTImprovedCMSExternal))/GetCorrectedJetPt(tmpJet, backgroundKTImprovedCMSExternal));
            FillHistogram("hJetProfile20GeV", 0.35-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.35, backgroundKTImprovedCMSExternal))/GetCorrectedJetPt(tmpJet, backgroundKTImprovedCMSExternal));
            FillHistogram("hJetProfile20GeV", 0.40-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.40, backgroundKTImprovedCMSExternal))/GetCorrectedJetPt(tmpJet, backgroundKTImprovedCMSExternal));
            FillHistogram("hJetProfile20GeV", 0.45-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.45, backgroundKTImprovedCMSExternal))/GetCorrectedJetPt(tmpJet, backgroundKTImprovedCMSExternal));
            FillHistogram("hJetProfile20GeV", 0.50-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.50, backgroundKTImprovedCMSExternal))/GetCorrectedJetPt(tmpJet, backgroundKTImprovedCMSExternal));
            FillHistogram("hJetProfile20GeV", 0.55-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.55, backgroundKTImprovedCMSExternal))/GetCorrectedJetPt(tmpJet, backgroundKTImprovedCMSExternal));
            FillHistogram("hJetProfile20GeV", 0.60-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.60, backgroundKTImprovedCMSExternal))/GetCorrectedJetPt(tmpJet, backgroundKTImprovedCMSExternal));
          }
          else if(GetCorrectedJetPt(tmpJet, backgroundKTImprovedCMSExternal)>=10.0)
          {
            FillHistogram("hJetProfile10GeV", 0.05-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.05, backgroundKTImprovedCMSExternal))/GetCorrectedJetPt(tmpJet, backgroundKTImprovedCMSExternal));
            FillHistogram("hJetProfile10GeV", 0.10-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.10, backgroundKTImprovedCMSExternal))/GetCorrectedJetPt(tmpJet, backgroundKTImprovedCMSExternal));
            FillHistogram("hJetProfile10GeV", 0.15-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.15, backgroundKTImprovedCMSExternal))/GetCorrectedJetPt(tmpJet, backgroundKTImprovedCMSExternal));
            FillHistogram("hJetProfile10GeV", 0.20-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.20, backgroundKTImprovedCMSExternal))/GetCorrectedJetPt(tmpJet, backgroundKTImprovedCMSExternal));
            FillHistogram("hJetProfile10GeV", 0.25-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.25, backgroundKTImprovedCMSExternal))/GetCorrectedJetPt(tmpJet, backgroundKTImprovedCMSExternal));
            FillHistogram("hJetProfile10GeV", 0.30-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.30, backgroundKTImprovedCMSExternal))/GetCorrectedJetPt(tmpJet, backgroundKTImprovedCMSExternal));
            FillHistogram("hJetProfile10GeV", 0.35-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.35, backgroundKTImprovedCMSExternal))/GetCorrectedJetPt(tmpJet, backgroundKTImprovedCMSExternal));
            FillHistogram("hJetProfile10GeV", 0.40-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.40, backgroundKTImprovedCMSExternal))/GetCorrectedJetPt(tmpJet, backgroundKTImprovedCMSExternal));
            FillHistogram("hJetProfile10GeV", 0.45-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.45, backgroundKTImprovedCMSExternal))/GetCorrectedJetPt(tmpJet, backgroundKTImprovedCMSExternal));
            FillHistogram("hJetProfile10GeV", 0.50-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.50, backgroundKTImprovedCMSExternal))/GetCorrectedJetPt(tmpJet, backgroundKTImprovedCMSExternal));
            FillHistogram("hJetProfile10GeV", 0.55-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.55, backgroundKTImprovedCMSExternal))/GetCorrectedJetPt(tmpJet, backgroundKTImprovedCMSExternal));
            FillHistogram("hJetProfile10GeV", 0.60-0.05/2, (GetCorrectedConePt(tmpJet->Eta(), tmpJet->Phi(), 0.60, backgroundKTImprovedCMSExternal))/GetCorrectedJetPt(tmpJet, backgroundKTImprovedCMSExternal));
          }
        }
        FillHistogram("hJetPtVsConstituentCount", tmpJet->Pt(),tmpJet->GetNumberOfTracks());

        if((fAnalyzeQA) && (tmpJet->Pt() >= 5.0))
        {
          Double_t lowestTrackPt = 1e99;
          Double_t highestTrackPt = 0.0;
          for(Int_t j=0; j<tmpJet->GetNumberOfTracks(); j++)
          {
            FillHistogram("hJetConstituentPt", tmpJet->TrackAt(j, fTrackArray)->Pt(), centralityPercentile);
            // Find the lowest pT of a track in the jet
            if (tmpJet->TrackAt(j, fTrackArray)->Pt() < lowestTrackPt)
              lowestTrackPt = tmpJet->TrackAt(j, fTrackArray)->Pt();
            if (tmpJet->TrackAt(j, fTrackArray)->Pt() > highestTrackPt)
              highestTrackPt = tmpJet->TrackAt(j, fTrackArray)->Pt();
          }
          FillHistogram("hJetArea", tmpJet->Area(), tmpJet->Pt());
          // Signal jet vs. signal jet - "Combinatorial"
          for (Int_t j = i+1; j<fNumberSignalJets; j++)
            FillHistogram("hJetDeltaPhi", GetDeltaPhi(tmpJet->Phi(), fSignalJets[j]->Phi()));

          FillHistogram("hJetPhiEta", tmpJet->Phi(),tmpJet->Eta());
          FillHistogram("hJetEta", tmpJet->Eta(), centralityPercentile);

          if(lowestTrackPt>=2.0)
            FillHistogram("hJetEta2GeVTracks", tmpJet->Eta(), centralityPercentile);
          if(lowestTrackPt>=4.0)
            FillHistogram("hJetEta4GeVTracks", tmpJet->Eta(), centralityPercentile);
        }


        // Jet constituent mass analyses (for PYTHIA)
        if (fAnalyzeMassCorrelation)
        {
          Double_t jetMass = 0.0;
          Int_t constituentCount = tmpJet->GetNumberOfTracks();
          Int_t electronCount = 0;
          Int_t pionCount = 0;
          Int_t protonCount = 0;
          Int_t kaonCount = 0;
          Int_t othersCount = 0;

          for(Int_t j=0; j<tmpJet->GetNumberOfTracks(); j++)
          {
            AliVParticle* tmpParticle = tmpJet->TrackAt(j, fTrackArray);

            jetMass += tmpParticle->M();
            if(TMath::Abs(tmpParticle->PdgCode()) == 11) // electron, positron
              electronCount++;
            else if(TMath::Abs(tmpParticle->PdgCode()) == 211) // p-, p+
              pionCount++;
            else if(TMath::Abs(tmpParticle->PdgCode()) == 2212) // p, pbar
              protonCount++;
            else if(TMath::Abs(tmpParticle->PdgCode()) == 321) // kaon+,kaon-
              kaonCount++;
            else
              othersCount++;
          }

          FillHistogram("hJetMassFromConstituents", jetMass);
          FillHistogram("hJetMass", tmpJet->M());

          FillHistogram("hJetPtVsMass", tmpJet->Pt(), jetMass);
          FillHistogram("hJetPtVsJetMass", tmpJet->Pt(), tmpJet->M());

          FillHistogram("hJetPtVsProtonCount", tmpJet->Pt(), protonCount/(static_cast<Double_t>(constituentCount)));
          FillHistogram("hJetPtVsPionCount", tmpJet->Pt(), pionCount/(static_cast<Double_t>(constituentCount)));
          FillHistogram("hJetPtVsKaonCount", tmpJet->Pt(), kaonCount/(static_cast<Double_t>(constituentCount)));
          FillHistogram("hJetPtVsElectronCount", tmpJet->Pt(), electronCount/(static_cast<Double_t>(constituentCount)));
          FillHistogram("hJetPtVsOthersCount", tmpJet->Pt(), othersCount/(static_cast<Double_t>(constituentCount)));
          FillHistogram("hJetConstituentProtonCount", protonCount);
          FillHistogram("hJetConstituentPionCount", pionCount);
          FillHistogram("hJetConstituentKaonCount", kaonCount);
          FillHistogram("hJetConstituentElectronCount", electronCount);
          FillHistogram("hJetConstituentOthersCount", othersCount);

          // Results for the "normal" jet (avoiding single particle jets)
          if((constituentCount>=6) && (constituentCount<=14))
          {
            FillHistogram("hJetPtVsMass_6_14", tmpJet->Pt(), jetMass);
            FillHistogram("hJetPtVsJetMass_6_14", tmpJet->Pt(), tmpJet->M());

            FillHistogram("hJetPtVsProtonCount_6_14", tmpJet->Pt(), protonCount);
            FillHistogram("hJetPtVsPionCount_6_14", tmpJet->Pt(), pionCount);
            FillHistogram("hJetPtVsKaonCount_6_14", tmpJet->Pt(), kaonCount);
            FillHistogram("hJetPtVsElectronCount_6_14", tmpJet->Pt(), electronCount);
            FillHistogram("hJetPtVsOthersCount_6_14", tmpJet->Pt(), othersCount);
            FillHistogram("hJetConstituentProtonCount_6_14", protonCount);
            FillHistogram("hJetConstituentPionCount_6_14", pionCount);
            FillHistogram("hJetConstituentKaonCount_6_14", kaonCount);
            FillHistogram("hJetConstituentElectronCount_6_14", electronCount);
            FillHistogram("hJetConstituentOthersCount_6_14", othersCount);
          }
          if((constituentCount>=2))
          {
            FillHistogram("hJetPtVsMass_2_X", tmpJet->Pt(), jetMass);
            FillHistogram("hJetPtVsJetMass_2_X", tmpJet->Pt(), tmpJet->M());

            FillHistogram("hJetPtVsProtonCount_2_X", tmpJet->Pt(), protonCount);
            FillHistogram("hJetPtVsPionCount_2_X", tmpJet->Pt(), pionCount);
            FillHistogram("hJetPtVsKaonCount_2_X", tmpJet->Pt(), kaonCount);
            FillHistogram("hJetPtVsElectronCount_2_X", tmpJet->Pt(), electronCount);
            FillHistogram("hJetPtVsOthersCount_2_X", tmpJet->Pt(), othersCount);
            FillHistogram("hJetConstituentProtonCount_2_X", protonCount);
            FillHistogram("hJetConstituentPionCount_2_X", pionCount);
            FillHistogram("hJetConstituentKaonCount_2_X", kaonCount);
            FillHistogram("hJetConstituentElectronCount_2_X", electronCount);
            FillHistogram("hJetConstituentOthersCount_2_X", othersCount);
          }
          if((constituentCount>=2) && (constituentCount<=7))
          {
            FillHistogram("hJetPtVsMass_2_7", tmpJet->Pt(), jetMass);
            FillHistogram("hJetPtVsJetMass_2_7", tmpJet->Pt(), tmpJet->M());

            FillHistogram("hJetPtVsProtonCount_2_7", tmpJet->Pt(), protonCount);
            FillHistogram("hJetPtVsPionCount_2_7", tmpJet->Pt(), pionCount);
            FillHistogram("hJetPtVsKaonCount_2_7", tmpJet->Pt(), kaonCount);
            FillHistogram("hJetPtVsElectronCount_2_7", tmpJet->Pt(), electronCount);
            FillHistogram("hJetPtVsOthersCount_2_7", tmpJet->Pt(), othersCount);
            FillHistogram("hJetConstituentProtonCount_2_7", protonCount);
            FillHistogram("hJetConstituentPionCount_2_7", pionCount);
            FillHistogram("hJetConstituentKaonCount_2_7", kaonCount);
            FillHistogram("hJetConstituentElectronCount_2_7", electronCount);
            FillHistogram("hJetConstituentOthersCount_2_7", othersCount);
          }
        }
      }
    }

    // ### SOME JET PLOTS
    FillHistogram("hJetCountAll", fJetArray->GetEntries());
    FillHistogram("hJetCountAccepted", fNumberSignalJets);
    FillHistogram("hJetCount", fJetArray->GetEntries(), fNumberSignalJets);
    if (fFirstLeadingJet)
    {
      FillHistogram("hLeadingJetPt", fFirstLeadingJet->Pt());
      FillHistogram("hCorrectedLeadingJetPt", GetCorrectedJetPt(fFirstLeadingJet,backgroundKTImprovedCMS));
    }
    if (fSecondLeadingJet)
    {
      FillHistogram("hSecondLeadingJetPt", fSecondLeadingJet->Pt());
      FillHistogram("hCorrectedSecondLeadingJetPt", GetCorrectedJetPt(fSecondLeadingJet,backgroundKTImprovedCMS));
    }
  } //endif AnalyzeJets

  #ifdef DEBUGMODE
    AliInfo("Calculate()::Jets done.");
  #endif
  ////////////////////// NOTE: Background analysis

  if (fAnalyzeBackground)
  {
    // Calculate background in centrality classes
    FillHistogram("hKTBackgroundImprovedCMS", backgroundKTImprovedCMS, centralityPercentile);
    FillHistogram("hKTBackgroundImprovedCMSExternal", backgroundKTImprovedCMSExternal, centralityPercentile);
    FillHistogram("hPPBackground", backgroundPP, centralityPercentile);
    FillHistogram("hKTMeanBackgroundImprovedCMS", centralityPercentile, backgroundKTImprovedCMS);

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

      // Calculate background profiles in terms of centrality
      FillHistogram("hKTMeanBackgroundPbPb", centralityPercentile,  backgroundKTPbPb);
      FillHistogram("hKTMeanBackgroundPbPbWithGhosts", centralityPercentile,  backgroundKTPbPbWithGhosts);
      FillHistogram("hKTMeanBackgroundCMS", centralityPercentile, backgroundKTCMS);
      FillHistogram("hKTMeanBackgroundMean", centralityPercentile, backgroundKTMean);
      FillHistogram("hKTMeanBackgroundTPC", centralityPercentile, backgroundKTTrackLike);
      FillHistogram("hTRMeanBackground", centralityPercentile,  backgroundTRCone06);
    }


    // Calculate the delta pt
    Double_t tmpDeltaPtNoBackground = GetDeltaPt(0.0);
    Double_t tmpDeltaPtKTImprovedCMS = GetDeltaPt(backgroundKTImprovedCMS);
    Double_t tmpDeltaPtExternalBgrd = GetDeltaPt(backgroundKTImprovedCMSExternal);


    Double_t tmpDeltaPtKTImprovedCMSPartialExclusion = 0.0;
    if(fNcoll)
      tmpDeltaPtKTImprovedCMSPartialExclusion = GetDeltaPt(backgroundKTImprovedCMS, 1.0/fNcoll);
    else
      tmpDeltaPtKTImprovedCMSPartialExclusion = GetDeltaPt(backgroundKTImprovedCMS, 1.0);

    Double_t tmpDeltaPtKTImprovedCMSPartialExclusion_Signal = 0.0;
    if(fNumberSignalJets)
      tmpDeltaPtKTImprovedCMSPartialExclusion_Signal = GetDeltaPt(backgroundKTImprovedCMS, 1.0/fNumberSignalJets);
    else
      tmpDeltaPtKTImprovedCMSPartialExclusion_Signal = GetDeltaPt(backgroundKTImprovedCMS, 1.0);
 
    Double_t tmpDeltaPtKTImprovedCMSFullExclusion = GetDeltaPt(backgroundKTImprovedCMS, 1.0);

    Double_t tmpDeltaPtKTPbPb = 0;
    Double_t tmpDeltaPtKTPbPbWithGhosts = 0;
    Double_t tmpDeltaPtKTCMS = 0;
    Double_t tmpDeltaPtKTMean = 0;
    Double_t tmpDeltaPtKTTrackLike = 0;
    Double_t tmpDeltaPtTR = 0;

    if(fAnalyzeDeprecatedBackgrounds)
    {
      tmpDeltaPtKTPbPb = GetDeltaPt(backgroundKTPbPb);
      tmpDeltaPtKTPbPbWithGhosts = GetDeltaPt(backgroundKTPbPbWithGhosts);
      tmpDeltaPtKTCMS = GetDeltaPt(backgroundKTCMS);
      tmpDeltaPtKTMean = GetDeltaPt(backgroundKTMean);
      tmpDeltaPtKTTrackLike = GetDeltaPt(backgroundKTTrackLike);
      tmpDeltaPtTR = GetDeltaPt(backgroundTRCone06);
    }

    // If valid, fill the delta pt histograms

    if(tmpDeltaPtExternalBgrd > -10000.0)
      FillHistogram("hDeltaPtExternalBgrd", tmpDeltaPtExternalBgrd, centralityPercentile);
    if(tmpDeltaPtKTImprovedCMS > -10000.0)
      FillHistogram("hDeltaPtKTImprovedCMS", tmpDeltaPtKTImprovedCMS, centralityPercentile);
    if(tmpDeltaPtKTImprovedCMSPartialExclusion > -10000.0)
      FillHistogram("hDeltaPtKTImprovedCMSPartialExclusion", tmpDeltaPtKTImprovedCMSPartialExclusion, centralityPercentile);
    if(tmpDeltaPtKTImprovedCMSPartialExclusion_Signal > -10000.0)
      FillHistogram("hDeltaPtKTImprovedCMSPartialExclusion_Signal", tmpDeltaPtKTImprovedCMSPartialExclusion_Signal, centralityPercentile);
    if(tmpDeltaPtKTImprovedCMSFullExclusion > -10000.0)
      FillHistogram("hDeltaPtKTImprovedCMSFullExclusion", tmpDeltaPtKTImprovedCMSFullExclusion, centralityPercentile);

    if(tmpDeltaPtNoBackground > 0.000001)
      FillHistogram("hDeltaPtNoBackgroundNoEmptyCones", tmpDeltaPtNoBackground, centralityPercentile);
    else if(tmpDeltaPtNoBackground > -10000.0)
      FillHistogram("hDeltaPtNoBackground", tmpDeltaPtNoBackground, centralityPercentile);


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
      if(tmpDeltaPtTR > -10000.0)
        FillHistogram("hDeltaPtTR", tmpDeltaPtTR, centralityPercentile);
    }
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
  // Implemented Notify() to read the cross sections
  // and number of trials from pyxsec.root
  // 
  #ifdef DEBUGMODE
    AliInfo("UserNotify started.");
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
  }
  #ifdef DEBUGMODE
    AliInfo("UserNotify ended.");
  #endif
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
    AliError(Form("Cannot find histogram <%s> ",key)) ;
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
  TH2* tmpHist = static_cast<TH2*>(fOutputList->FindObject(GetHistoName(key)));
  if(!tmpHist)
  {
    AliError(Form("Cannot find histogram <%s> ",key));
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
