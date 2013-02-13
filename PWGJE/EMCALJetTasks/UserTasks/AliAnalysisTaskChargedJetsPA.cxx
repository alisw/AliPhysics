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
  // NOTE: Track & Cluster & QA histograms
  if (fAnalyzeQA)
  {
  
    AddHistogram1D<TH1D>("hNumberEvents", "Number of events (0 = before, 1 = after vertex cuts)", "", 2, 0, 2, "#Delta z(cm)","N^{Events}/cut");
      
    AddHistogram1D<TH1D>("hAppliedEtaCorrectionFactor", "Applied #eta correction factor for the k_{T} background", "", 500, 0.5, 1.5, "Correction factor","dN^{Jets}/df");
    AddHistogram1D<TH1D>("hAppliedEtaCorrectionFactor2", "Applied #eta correction factor for the k_{T} background 2", "", 500, 0.5, 1.5, "Correction factor","dN^{Jets}/df");
    AddHistogram1D<TH1D>("hVertexZ", "Z distribution of the vertex", "", 400, -40., 40., "#Delta z(cm)","dN^{Events}/dz");

    AddHistogram1D<TH1D>("hVertexR", "R distribution of the vertex", "", 100, 0., 1., "#Delta r(cm)","dN^{Events}/dr");
    AddHistogram1D<TH1D>("hCentrality", "Centrality distribution", "", 5, 0., 100., "Centrality (classes)","dN^{Events}");

    AddHistogram2D<TH2D>("hTrackCountAcc", "Number of tracks in acceptance vs. centrality", "LEGO2", 750, 0., 750., 5, 0, 100, "N tracks","Centrality", "dN^{Events}/dN^{Tracks}");
    AddHistogram2D<TH2D>("hTrackPhiEta", "Track angular distribution", "LEGO2", 100, 0., 2*TMath::Pi(),100, -2.5, 2.5, "#phi","#eta","dN^{Tracks}/(d#phi d#eta)");
    AddHistogram1D<TH1D>("hTrackPt", "Tracks p_{T} distribution", "", 20000, 0., 200., "p_{T} (GeV/c)","dN^{Tracks}/dp_{T}");
    AddHistogram1D<TH1D>("hTrackCharge", "Charge", "", 11, -5, 5, "Charge (e)","dN^{Tracks}/dq");
    AddHistogram1D<TH1D>("hTrackEta", "Track #eta distribution", "", 180, -fTrackEtaWindow, +fTrackEtaWindow, "#eta","dN^{Tracks}/d#eta");

    AddHistogram2D<TH2D>("hClusterCountAcc", "Number of clusters in acceptance vs. centrality", "LEGO2", 750, 0., 750., 5, 0, 100, "N clusters","Centrality", "dN^{Events}/dN^{Clusters}");
    AddHistogram1D<TH1D>("hClusterE", "Clusters energy distribution", "", 20000, 0., 200., "p_{T} (GeV/c)","dN^{Cluster}/dp_{T}");
  }

  // NOTE: Pythia histograms
  if (fAnalyzePythia)
  {
    AddHistogram1D<TH1D>("hPythiaPtHard", "Pythia p_{T} hard distribution", "", 2000, 0, 400, "p_{T} hard","dN^{Events}/dp_{T,hard}");
    AddHistogram1D<TProfile>("hPythiaXSection", "Pythia cross section distribution", "", fNumPtHardBins+2, -1, fNumPtHardBins+1, "p_{T} hard bin","dN^{Events}/dp_{T,hard}");
    AddHistogram1D<TH1D>("hPythiaNTrials", "Pythia trials (no correction for manual cuts)", "", fNumPtHardBins+2, -1, fNumPtHardBins+1, "p_{T} hard bin", "Trials");
  }

  // NOTE: Jet histograms
  if (fAnalyzeJets)
  {
    // ######## Jet spectra
    AddHistogram1D<TH1D>("hJetPt", "Jets p_{T} distribution", "", 1000, 0., 200., "p_{T} (GeV/c)","dN^{Jets}/dp_{T}");
    AddHistogram1D<TH1D>("hJetPtBgrdSubtractedRC", "Jets p_{T} distribution, RC background subtracted", "", 500, -50., 200., "p_{T} (GeV/c)","dN^{Jets}/dp_{T}");
    AddHistogram1D<TH1D>("hJetPtBgrdSubtractedKT", "Jets p_{T} distribution, KT background subtracted, corrected for eta dependence)", "", 500, -50., 200., "p_{T} (GeV/c)","dN^{Jets}/dp_{T}");
    AddHistogram1D<TH1D>("hJetPtBgrdSubtractedKTNoEtaCorr", "Jets p_{T} distribution, KT background subtracted", "", 500, -50., 200., "p_{T} (GeV/c)","dN^{Jets}/dp_{T}");
    AddHistogram1D<TH1D>("hJetPtBgrdSubtractedKT2", "Jets p_{T} distribution, KT background 2 subtracted, corrected for eta dependence)", "", 500, -50., 200., "p_{T} (GeV/c)","dN^{Jets}/dp_{T}");
    AddHistogram1D<TH1D>("hJetPtBgrdSubtractedKT2NoEtaCorr", "Jets p_{T} distribution, KT background 2 subtracted", "", 500, -50., 200., "p_{T} (GeV/c)","dN^{Jets}/dp_{T}");

    AddHistogram1D<TH1D>("hJetPtBgrdSubtractedTR", "Jets p_{T} distribution, Track background subtracted", "", 500, -50., 200., "p_{T} (GeV/c)","dN^{Jets}/dp_{T}");


    AddHistogram1D<TH1D>("hJetArea", "Jets area distribution", "", 200, 0., 2., "Area","dN^{Jets}/dA");
    AddHistogram2D<TH2D>("hJetPtArea", "Jets p_{T} distribution", "LEGO2", 1000, 0., 200.,100, 0., 1., "p_{T} (GeV/c)","Area","dN^{Jets}/(dp_{T}dA)");
    AddHistogram1D<TH1D>("hJetDeltaPhi", "Jets combinatorial #Delta #phi", "", 250, 0., TMath::Pi(), "#Delta #phi","dN^{Jets}/d(#Delta #phi)");
    AddHistogram2D<TH2D>("hJetDeltaPhiPt", "Jets combinatorial #Delta #phi vs. p_{T}", "LEGO2", 250, 0., TMath::Pi(), 20, 0.,100., "#Delta #phi","max(p_{T,1},p_{T,2}) (GeV/c)","dN^{Jets}/d(#Delta #phi)dp_{T}");
    AddHistogram1D<TH1D>("hLeadingJetDeltaPhi", "1st and 2nd leading jet #Delta #phi", "", 250, 0., TMath::Pi(), "#Delta #phi","dN^{Jets}/d(#Delta #phi)");
    AddHistogram2D<TH2D>("hLeadingJetDeltaPhiPt", "1st and 2nd leading jet #Delta #phi vs. p_{T}", "LEGO2", 250, 0., TMath::Pi(),20, 0.,100., "#Delta #phi","1st leading p_{T} (GeV/c)","dN^{Jets}/d(#Delta #phi)dp_{T}");
    AddHistogram2D<TH2D>("hJetPtEta", "Jets p_{T} distribution", "LEGO2", 1000, 0., 200.,100, -0.6, 0.6, "p_{T} (GeV/c)","#eta","dN^{Jets}/(dp_{T}d#eta)");
    AddHistogram2D<TH2D>("hJetPtPhi", "Jets p_{T} #phi distribution", "LEGO2", 1000, 0., 200.,100, 0.0, TMath::TwoPi(), "p_{T} (GeV/c)","#phi","dN^{Jets}/(dp_{T}d#phi)");
    AddHistogram2D<TH2D>("hJetPtCentrality", "Jets p_{T} distribution", "LEGO2", 1000, 0., 200.,5, 0, 100, "p_{T} (GeV/c)","Centrality","dN^{Jets}/dp_{T}");
    AddHistogram2D<TH2D>("hJetPhiEta", "Jets angular distribution", "LEGO2", 100, 0., 2*TMath::Pi(),100, -0.6, 0.6, "#phi","#eta","dN^{Jets}/(d#phi d#eta)");
    AddHistogram1D<TH1D>("hJetCountAll", "Number of Jets", "", 200, 0., 200., "N jets","dN^{Events}/dN^{Jets}");
    AddHistogram1D<TH1D>("hJetCountAccepted", "Number of accepted Jets", "", 200, 0., 200., "N jets","dN^{Events}/dN^{Jets}");

    AddHistogram1D<TH1D>("hLeadingJetPt", "Leading jet p_{T}", "", 500,  0, 100, "p_{T} (GeV/c)","dN^{Jets}/dp_{T}");
    AddHistogram1D<TH1D>("hSecondLeadingJetPt", "Second Leading jet p_{T}", "", 500,  0, 100, "p_{T} (GeV/c)","dN^{Jets}/dp_{T}");

    AddHistogram1D<TH1D>("hDijetConstituentsPt", "Dijet constituents p_{T} distribution", "", 500, 0., 100., "p_{T} (GeV/c)","dN^{Jets}/dp_{T}");
    AddHistogram1D<TH1D>("hDijetLeadingJetPt", "Dijet leading jet p_{T} distribution", "", 500, 0., 100., "p_{T} (GeV/c)","dN^{Jets}/dp_{T}");
    AddHistogram2D<TH2D>("hDijetPtCorrelation", "Dijet constituents p_{T} correlation", "COLZ", 500, 5., 100., 500, 5., 100., "1st leading jet p_{T} (GeV/c)","2nd leading jet p_{T} (GeV/c)","dN^{Dijets}/d^{2}p_{T}");

    AddHistogram1D<TH1D>("hDijet2ConstituentsPt", "Dijet2 constituents p_{T} distribution", "", 500, 0., 100., "p_{T} (GeV/c)","dN^{Jets}/dp_{T}");
    AddHistogram1D<TH1D>("hDijet2LeadingJetPt", "Dijet2 leading jet p_{T} distribution", "", 500, 0., 100., "p_{T} (GeV/c)","dN^{Jets}/dp_{T}");
    AddHistogram2D<TH2D>("hDijet2PtCorrelation", "Dijet2 constituents p_{T} correlation", "COLZ", 500, 5., 100., 500, 5., 100., "1st leading jet p_{T} (GeV/c)","2nd leading jet p_{T} (GeV/c)","dN^{Dijets}/d^{2}p_{T}");

  }
  // NOTE: Jet background histograms
  if (fAnalyzeBackground)
  {
    // ########## Delta Pt
    AddHistogram1D<TH1D>("hDeltaPtKT", "Background fluctuations #delta p_{T} (KT, 0 jets excluded)", "", 500, -20., 80., "#delta p_{T} (GeV/c)","dN^{Jets}/d#delta p_{T}");
    AddHistogram1D<TH1D>("hDeltaPtKT1Excl", "Background fluctuations #delta p_{T} (KT, 1 jets excluded)", "", 500, -20., 80., "#delta p_{T} (GeV/c)","dN^{Jets}/d#delta p_{T}");
    AddHistogram1D<TH1D>("hDeltaPtKT2Excl", "Background fluctuations #delta p_{T} (KT, 2 jets excluded)", "", 500, -20., 80., "#delta p_{T} (GeV/c)","dN^{Jets}/d#delta p_{T}");

    Double_t dptEtaMin = -(fTrackEtaWindow-fRandConeRadius) + 2*(fTrackEtaWindow-fRandConeRadius)/fBackgroundEtaBins *  fKTDeltaPtEtaBin;
    Double_t dptEtaMax = -(fTrackEtaWindow-fRandConeRadius) + 2*(fTrackEtaWindow-fRandConeRadius)/fBackgroundEtaBins * (fKTDeltaPtEtaBin+1);

    AddHistogram1D<TH1D>("hDeltaPtKTEta", Form("Background fluctuations #delta p_{T} (KT, 0 jets excluded, #eta=%1.3f to %1.3f)", dptEtaMin,dptEtaMax), "", 500, -20., 80., "#delta p_{T} (GeV/c)","dN^{Jets}/d#delta p_{T}");
    AddHistogram1D<TH1D>("hDeltaPtKTEta1Excl", Form("Background fluctuations #delta p_{T} (KT, 1 jets excluded, #eta=%1.3f to %1.3f)", dptEtaMin,dptEtaMax), "", 500, -20., 80., "#delta p_{T} (GeV/c)","dN^{Jets}/d#delta p_{T}");
    AddHistogram1D<TH1D>("hDeltaPtKTEta2Excl", Form("Background fluctuations #delta p_{T} (KT, 2 jets excluded, #eta=%1.3f to %1.3f)", dptEtaMin,dptEtaMax), "", 500, -20., 80., "#delta p_{T} (GeV/c)","dN^{Jets}/d#delta p_{T}");

    AddHistogram1D<TH1D>("hDeltaPtKT2Eta2Excl", Form("Background fluctuations #delta p_{T} (KT 2, 2 jets excluded, #eta=%1.3f to %1.3f)", dptEtaMin,dptEtaMax), "", 500, -20., 80., "#delta p_{T} (GeV/c)","dN^{Jets}/d#delta p_{T}");

    AddHistogram1D<TH1D>("hDeltaPtRC", "Background fluctuations #delta p_{T} (RC, 0 jets excluded)", "", 500, -20., 80., "#delta p_{T} (GeV/c)","dN^{Jets}/d#delta p_{T}");
    AddHistogram1D<TH1D>("hDeltaPtRC1Excl", "Background fluctuations #delta p_{T} (RC, 1 jets excluded)", "", 500, -20., 80., "#delta p_{T} (GeV/c)","dN^{Jets}/d#delta p_{T}");
    AddHistogram1D<TH1D>("hDeltaPtRC2Excl", "Background fluctuations #delta p_{T} (RC, 2 jets excluded)", "", 500, -20., 80., "#delta p_{T} (GeV/c)","dN^{Jets}/d#delta p_{T}");

    AddHistogram1D<TH1D>("hDeltaPtTR", "Background fluctuations #delta p_{T} (TR, 0 jets excluded)", "", 500, -20., 80., "#delta p_{T} (GeV/c)","dN^{Jets}/d#delta p_{T}");
    AddHistogram1D<TH1D>("hDeltaPtTR1Excl", "Background fluctuations #delta p_{T} (TR, 1 jets excluded)", "", 500, -20., 80., "#delta p_{T} (GeV/c)","dN^{Jets}/d#delta p_{T}");
    AddHistogram1D<TH1D>("hDeltaPtTR2Excl", "Background fluctuations #delta p_{T} (TR, 2 jets excluded)", "", 500, -20., 80., "#delta p_{T} (GeV/c)","dN^{Jets}/d#delta p_{T}");



    AddHistogram2D<TH2D>("hKTJetPhiEta", "KT Jets angular distribution", "LEGO2", 100, 0., 2*TMath::Pi(),100, -0.6, 0.6, "#phi","#eta","dN^{Jets}/(d#phi d#eta)");
    AddHistogram2D<TH2D>("hKTLeadingJetPhiEta", "KT Leading jets angular distribution", "LEGO2", 100, 0., 2*TMath::Pi(),100, -0.6, 0.6, "#phi","#eta","dN^{Jets}/(d#phi d#eta)");

    AddHistogram1D<TH1D>("hDijetBackground", "Background density (dijets excluded)", "", 400, 0., 40., "#rho (GeV/c)","dN^{Events}/d#rho");
    AddHistogram1D<TH1D>("hDijetBackgroundMostCentral", "Background density (0-20 centrality, dijets excluded)", "", 400, 0., 40., "#rho (GeV/c)","dN^{Events}/d#rho");
    AddHistogram2D<TH2D>("hDijetBackgroundVsCentrality", "Background density vs. centrality (dijets excluded)", "", 200, 0., 20., 5, 0, 100, "#rho (GeV/c)","Centrality", "dN^{Events}/d#rho");

    AddHistogram1D<TH1D>("hDijetBackgroundPerpendicular", "Background density (dijets excluded)", "", 400, 0., 40., "#rho (GeV/c)","dN^{Events}/d#rho");
    AddHistogram1D<TH1D>("hDijetBackgroundPerpendicularMostCentral", "Background density (0-20 centrality, dijets excluded)", "", 400, 0., 40., "#rho (GeV/c)","dN^{Events}/d#rho");
    AddHistogram2D<TH2D>("hDijetBackgroundPerpendicularVsCentrality", "Background density vs. centrality (dijets excluded)", "", 200, 0., 20., 5, 0, 100, "#rho (GeV/c)","Centrality", "dN^{Events}/d#rho");
  
    AddHistogram2D<TH2D>("hRCBackground", "RC background density (2 leading jets excluded, mean(8 RCs))", "LEGO2", fBackgroundEtaBins, -(fTrackEtaWindow-fRandConeRadius), +(fTrackEtaWindow-fRandConeRadius), 400, 0., 40., "#eta", "#rho (GeV/c)","dN^{Events}/d#rho");

    AddHistogram1D<TH1D>("hAccConesInRCBackground", Form("Number of cones used for RC background (|#eta| < %1.1f)", fSignalJetEtaWindow), "", 8, 0, 8, "Used cones", "dN^{Events}/dN^{Cones}");

    AddHistogram2D<TH2D>("hRCBackgroundMostCentral", "RC background density (0-20 centrality, 2 leading jets excluded)", "LEGO2", fBackgroundEtaBins, -(fTrackEtaWindow-fRandConeRadius), +(fTrackEtaWindow-fRandConeRadius), 400, 0., 40., "#eta",  "#rho (GeV/c)","dN^{Events}/d#rho");
    AddHistogram2D<TH2D>("hRCBackgroundMostPeripheral", "RC background density (80-100 centrality, 2 leading jets excluded)", "LEGO2", fBackgroundEtaBins, -(fTrackEtaWindow-fRandConeRadius), +(fTrackEtaWindow-fRandConeRadius),  400, 0., 40., "#eta", "#rho (GeV/c)","dN^{Events}/d#rho");
    AddHistogram2D<TH2D>("hRCBackgroundVsCentrality", "RC background density vs centrality (2 leading jets excluded)", "LEGO2", 200, 0., 20., 5, 0, 100, "#rho (GeV/c)","Centrality", "dN^{Events}/d#rho");

    
    AddHistogram2D<TH2D>("hKTBackground", "KT background density (2 leading jets excluded)", "LEGO2", fBackgroundEtaBins, -(fTrackEtaWindow-fRandConeRadius), +(fTrackEtaWindow-fRandConeRadius), 400, 0., 40., "#eta", "#rho (GeV/c)","dN^{Events}/d#rho");
    AddHistogram2D<TH2D>("hKTBackgroundMostCentral", "KT background density (0-20 centrality, 2 leading jets excluded)", "LEGO2", fBackgroundEtaBins, -(fTrackEtaWindow-fRandConeRadius), +(fTrackEtaWindow-fRandConeRadius), 400, 0., 40., "#eta",  "#rho (GeV/c)","dN^{Events}/d#rho");
    AddHistogram2D<TH2D>("hKTBackgroundMostPeripheral", "KT background density (80-100 centrality, 2 leading jets excluded)", "LEGO2", fBackgroundEtaBins, -(fTrackEtaWindow-fRandConeRadius), +(fTrackEtaWindow-fRandConeRadius),  400, 0., 40., "#eta", "#rho (GeV/c)","dN^{Events}/d#rho");
    AddHistogram2D<TH2D>("hKTBackgroundVsCentrality", "KT background density vs centrality (2 leading jets excluded)", "LEGO2", 200, 0., 20., 5, 0, 100, "#rho (GeV/c)","Centrality", "dN^{Events}/d#rho");
  
    AddHistogram2D<TH2D>("hKTBackground2", "KT background 2 density (2 leading jets excluded)", "LEGO2", fBackgroundEtaBins, -(fTrackEtaWindow-fRandConeRadius), +(fTrackEtaWindow-fRandConeRadius), 400, 0., 40., "#eta", "#rho (GeV/c)","dN^{Events}/d#rho");
    AddHistogram2D<TH2D>("hKTBackground2MostCentral", "KT background 2 density (0-20 centrality, 2 leading jets excluded)", "LEGO2", fBackgroundEtaBins, -(fTrackEtaWindow-fRandConeRadius), +(fTrackEtaWindow-fRandConeRadius), 400, 0., 40., "#eta",  "#rho (GeV/c)","dN^{Events}/d#rho");
    AddHistogram2D<TH2D>("hKTBackground2MostPeripheral", "KT background 2 density (80-100 centrality, 2 leading jets excluded)", "LEGO2", fBackgroundEtaBins, -(fTrackEtaWindow-fRandConeRadius), +(fTrackEtaWindow-fRandConeRadius),  400, 0., 40., "#eta", "#rho (GeV/c)","dN^{Events}/d#rho");
    AddHistogram2D<TH2D>("hKTBackground2VsCentrality", "KT background 2 density vs centrality (2 leading jets excluded)", "LEGO2", 200, 0., 20., 5, 0, 100, "#rho (GeV/c)","Centrality", "dN^{Events}/d#rho");


    AddHistogram2D<TH2D>("hTrackBackground", "Track background density (2 leading jets excluded)", "LEGO2", fBackgroundEtaBins, -(fTrackEtaWindow-fRandConeRadius), +(fTrackEtaWindow-fRandConeRadius), 400, 0., 40., "#eta", "#rho (GeV/c)","dN^{Events}/d#rho");
    AddHistogram2D<TH2D>("hTrackBackgroundMostCentral", "Track background density (0-20 centrality, 2 leading jets excluded)", "LEGO2", fBackgroundEtaBins,  -(fTrackEtaWindow-fRandConeRadius), +(fTrackEtaWindow-fRandConeRadius), 400, 0., 40., "#eta", "#rho (GeV/c)","dN^{Events}/d#rho");
    AddHistogram2D<TH2D>("hTrackBackgroundMostPeripheral", "Track background density (80-100 centrality, 2 leading jets excluded)", "LEGO2", fBackgroundEtaBins, -(fTrackEtaWindow-fRandConeRadius), +(fTrackEtaWindow-fRandConeRadius), 400, 0., 40., "#eta", "#rho (GeV/c)","dN^{Events}/d#rho");
    AddHistogram2D<TH2D>("hTrackBackgroundVsCentrality", "Track background density vs centrality (2 leading jets excluded)", "LEGO2", 200, 0., 20., 5, 0, 100, "#rho (GeV/c)","Centrality", "dN^{Events}/d#rho");

  }



  // register Histograms
  for (Int_t i = 0; i < fHistCount; i++)
  {
    fOutputList->Add(fHistList->At(i));
  }
  
  PostData(1,fOutputList); // important for merging

}

//________________________________________________________________________
AliAnalysisTaskChargedJetsPA::AliAnalysisTaskChargedJetsPA(const char *name, const char* trackArrayName, const char* clusterArrayName, const char* jetArrayName, const char* backgroundJetArrayName) : AliAnalysisTaskSE(name), fOutputList(0), fAnalyzeQA(1), fAnalyzeJets(1), fAnalyzeBackground(1), fAnalyzePythia(0), fHasTracks(0), fHasClusters(0), fHasJets(0), fHasBackgroundJets(0), fIsMC(0), fJetArray(0), fTrackArray(0), fClusterArray(0), fBackgroundJetArray(0), fJetArrayName(0), fTrackArrayName(0), fClusterArrayName(0), fBackgroundJetArrayName(0), fNumPtHardBins(11), fRandConeRadius(0.4), fSignalJetRadius(0.4), fBackgroundJetRadius(0.4),  fKTDeltaPtEtaBin(3), fTrackBackgroundConeRadius(0.4), fNumberRandCones(8), fNumberExcludedJets(2), fDijetMaxAngleDeviation(10.0), fBackgroundEtaBins(5), fJetBgrdCorrectionFactors(0), fSignalJetEtaWindow(0.5), fBackgroundJetEtaWindow(0.5), fTrackEtaWindow(0.9), fClusterEtaWindow(0.7), fVertexWindow(10.0), fVertexMaxR(1.0), fMinTrackPt(0.150), fMinClusterPt(0.300), fMinJetPt(1.0), fMinJetArea(0.4), fMinBackgroundJetPt(0.15), fMinDijetLeadingPt(10.0), fFirstLeadingJet(0), fSecondLeadingJet(0), fNumberSignalJets(0), fCrossSection(0.0), fTrials(0.0),  fRandom(0), fHelperClass(0), fInitialized(0), fTaskInstanceCounter(0), fHistList(0), fHistCount(0)
{
  #ifdef DEBUGMODE
    AliInfo("Calling constructor.");
  #endif

  // Constructor
  // Define input and output slots here (never in the dummy constructor)
  // Input slot #0 works with a TChain - it is connected to the default input container
  // Output slot #1 writes into a TH1 container
  // Constructor

  // Every instance of this task gets his own number
  static Int_t instance = 0;
  fTaskInstanceCounter = instance;
  instance++;

  fTrackArrayName = new TString(trackArrayName);
  fClusterArrayName = new TString(clusterArrayName);
  if (strcmp(fTrackArrayName->Data(),"") == 0)
    fAnalyzeQA = kFALSE;
  else
  {
    fAnalyzeQA = kTRUE;
    if (fTrackArrayName->Contains("MCParticles")) //TODO: Hardcoded for now
      fIsMC = kTRUE;
  }

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
      AliError("Pythia Header not accessible!");
    else
      tmpPtHard = pythiaHeader->GetPtHard();
  }
  return tmpPtHard;
}

//________________________________________________________________________
inline Int_t AliAnalysisTaskChargedJetsPA::GetPtHardBin()
{
  // ########## PT HARD BIN EDGES
  const Int_t localkPtHardLowerEdges[] = { 0, 5,11,21,36,57, 84,117,152,191,234};
  const Int_t localkPtHardHigherEdges[] = { 5,11,21,36,57,84,117,152,191,234,1000000};

  Int_t tmpPtHardBin = 0;
  Double_t tmpPtHard = GetPtHard();
 
  for (tmpPtHardBin = 0; tmpPtHardBin <= fNumPtHardBins; tmpPtHardBin++)
    if (tmpPtHard >= localkPtHardLowerEdges[tmpPtHardBin] && tmpPtHard < localkPtHardHigherEdges[tmpPtHardBin])
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
inline Bool_t AliAnalysisTaskChargedJetsPA::IsClusterInAcceptance(AliVCluster* cluster)
{
  if (cluster != 0)
 //   if (TMath::Abs(cluster->Eta()) <= fClusterEtaWindow)
//      if (cluster->Phi() <= 187.0/360.0 * TMath::TwoPi());
//        if (cluster->Phi() >= 80.0/360.0 * TMath::TwoPi());
          if (cluster->E() >= fMinClusterPt)
            return kTRUE;

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
      AliInfo(Form("%s: Could not retrieve tracks %s! This is OK, if tracks are not demanded.", GetName(), fTrackArrayName->Data())); 
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
  // Check for cluster array
  if (strcmp(fClusterArrayName->Data(), "") != 0)
  {
    fClusterArray = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(fClusterArrayName->Data()));
    fHasClusters = kTRUE;
    if (!fClusterArray) 
    {
      AliInfo(Form("%s: Could not retrieve clusters %s! This is OK, if clusters are not demanded.", GetName(), fClusterArrayName->Data())); 
      fHasClusters = kFALSE;
    } 
    else
    {
      TClass *cl = fClusterArray->GetClass();
      if (!cl->GetBaseClass("AliVCluster"))
      {
      	AliError(Form("%s: Collection %s does not contain AliVCluster objects!", GetName(), fClusterArrayName->Data())); 
      	fClusterArray = 0;
        fHasClusters = kFALSE;
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
      AliInfo(Form("%s: Could not retrieve jets %s! This is OK, if jets are not demanded.", GetName(), fJetArrayName->Data())); 
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
  if ((!fHasTracks && fAnalyzeQA) || (!fHasTracks && fAnalyzeBackground))
  {
    AliError(Form("%s: Tracks NOT successfully casted although demanded! Deactivating QA and background analysis",GetName()));
    fAnalyzeQA = kFALSE;
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

  Float_t maxJetPts[] = {0, 0};
  Int_t jetIDArray[]   = {-1, -1};
  Int_t jetCount = fJetArray->GetEntries();

  // Go through all jets and save signal jets and the two leading ones
  for (Int_t i = 0; i < jetCount; i++)
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
    fSignalJets[fNumberSignalJets] = jet;
    fNumberSignalJets++;
  }
  
  if (fNumberSignalJets > 0)
    fFirstLeadingJet  = static_cast<AliEmcalJet*>(fJetArray->At(jetIDArray[0]));
  if (fNumberSignalJets > 1)
    fSecondLeadingJet = static_cast<AliEmcalJet*>(fJetArray->At(jetIDArray[1]));

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
Double_t AliAnalysisTaskChargedJetsPA::GetJetBackgroundCorrFactor(Double_t eta, Double_t background)
{
  Double_t tmpCorrFactor = 1.0;

  if(fJetBgrdCorrectionFactors)
    tmpCorrFactor = fJetBgrdCorrectionFactors->GetBinContent
                          (
                            fJetBgrdCorrectionFactors->GetXaxis()->FindBin(eta),
                            fJetBgrdCorrectionFactors->GetYaxis()->FindBin(background)
                          );

  return tmpCorrFactor;
}
//________________________________________________________________________
Double_t AliAnalysisTaskChargedJetsPA::GetCorrectedJetPt(AliEmcalJet* jet, Double_t background, Bool_t useEtaCorrection)
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

  // Get correction factor from saved histo or similar in dependence of jet eta and background density
  Double_t corrfactor = 1.0;
  if(useEtaCorrection)
  {
    corrfactor = GetJetBackgroundCorrFactor(jet->Eta(), background);
  }

  // Get Eta corrected background
  Double_t tmpCorrectedBackground = background * corrfactor;

  // Subtract background
  correctedPt = jet->Pt() - tmpCorrectedBackground * jet->Area();

  #ifdef DEBUGMODE
    AliInfo("Got corrected jet spectra.");
  #endif 

  return correctedPt;
}


//________________________________________________________________________
void AliAnalysisTaskChargedJetsPA::GetDeltaPt(Double_t& deltaPt, Double_t rho, Int_t numberExcludeLeadingJets, Int_t usedEtaBin, Bool_t useEtaCorrection)
{
  #ifdef DEBUGMODE
    AliInfo("Getting Delta Pt.");
  #endif

  // Define the tmp delta pt
  deltaPt = -10000.0;

  // Exclude UP TO numberExcludeLeadingJets
  if (fNumberSignalJets < 2)
    numberExcludeLeadingJets = fNumberSignalJets;
  
  Double_t etaMin, etaMax;
  if (usedEtaBin==-1)
  {
    etaMin = -(fTrackEtaWindow-fRandConeRadius);
    etaMax = +(fTrackEtaWindow-fRandConeRadius);
  }
  else
  {
    etaMin = -(fTrackEtaWindow-fRandConeRadius) + 2*(fTrackEtaWindow-fRandConeRadius)/fBackgroundEtaBins *  usedEtaBin;
    etaMax = -(fTrackEtaWindow-fRandConeRadius) + 2*(fTrackEtaWindow-fRandConeRadius)/fBackgroundEtaBins * (usedEtaBin+1);
  }  


  Double_t tmpRandConeEta = 0.0;
  Double_t tmpRandConePhi = 0.0;

  Bool_t coneValid = kTRUE;


  tmpRandConeEta = etaMin + fRandom->Rndm()*(etaMax-etaMin);
  tmpRandConePhi = fRandom->Rndm()*TMath::TwoPi();

  // Apply eta correction on demand
  if(useEtaCorrection)
    rho = GetJetBackgroundCorrFactor(tmpRandConeEta, rho)*rho;

  // Go through all excluded leading jets and check if there's an overlap
  for(Int_t j=0;j<numberExcludeLeadingJets;j++)
  {
    AliEmcalJet* tmpJet = NULL;

    if (j==0)
      tmpJet = fFirstLeadingJet;
    else if (j==1)
      tmpJet = fSecondLeadingJet;
    else
      AliFatal("Trying to exclude more than 2 jets for delta pt -- not implemented.");

    Double_t excludedJetPhi = tmpJet->Phi();
    Double_t excludedJetEta = tmpJet->Eta();
    Double_t tmpDeltaPhi = GetDeltaPhi(tmpRandConePhi, excludedJetPhi);

    // Check, if cone has overlap with jet
    if ( tmpDeltaPhi*tmpDeltaPhi + TMath::Abs(tmpRandConeEta-excludedJetEta)*TMath::Abs(tmpRandConeEta-excludedJetEta) <= fRandConeRadius*fRandConeRadius)
    {
      coneValid = kFALSE;
      break;
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
void AliAnalysisTaskChargedJetsPA::GetKTBackgroundDensity(Int_t numberExcludeLeadingJets, Double_t& rhoMedian, Double_t& areaMean, Double_t etaMin, Double_t etaMax)
{
  #ifdef DEBUGMODE
    AliInfo("Getting KT background density.");
  #endif

  // static declaration. Advantage: more speed. Disadvantage: Problematic for events with more than 1024 jets :)
  static Double_t tmpRhos[1024];
  static Double_t tmpAreas[1024];
  Int_t maxJetIds[]   = {-1, -1}; // Indices for excludes jets (up to two)

  // Setting invalid values
  rhoMedian = -1.0;
  areaMean= -1.0;

  // Exclude UP TO numberExcludeLeadingJets
  Int_t numberBgrdJets = GetLeadingJets(fBackgroundJetArray, &maxJetIds[0], kFALSE);
  if (numberBgrdJets < numberExcludeLeadingJets)
    numberExcludeLeadingJets = numberBgrdJets;
  if ((etaMin == 0) && (etaMax == 0))
  {
    etaMin = -fBackgroundJetEtaWindow;
    etaMax = +fBackgroundJetEtaWindow;
  }

  Int_t jetCountAccepted = 0;
  Int_t jetCount = fBackgroundJetArray->GetEntries();

  for (Int_t i = 0; i < jetCount; i++)
  {
    AliEmcalJet* jet = static_cast<AliEmcalJet*>(fBackgroundJetArray->At(i));
    if (!jet)
    {
      AliError(Form("%s: Could not receive jet %d", GetName(), i));
      continue;
    } 

    // exclude leading jets
    if (numberExcludeLeadingJets > 0)
      if (i == maxJetIds[0])
        continue;
    if (numberExcludeLeadingJets > 1)
      if (i == maxJetIds[1])
        continue;
      


    if (!IsBackgroundJetInAcceptance(jet))
      continue;
    if (!((jet->Eta() >= etaMin) && (jet->Eta() < etaMax)))
      continue;

    
    tmpRhos[jetCountAccepted] = jet->Pt() / jet->Area();
    tmpAreas[jetCountAccepted] = jet->Area();
    jetCountAccepted++;
  }

  if (jetCountAccepted > 0)
  {
    rhoMedian = TMath::Median(jetCountAccepted, tmpRhos);
    areaMean   = TMath::Mean(jetCountAccepted, tmpAreas);
  }
  #ifdef DEBUGMODE
    AliInfo("Got KT background density.");
  #endif
}

//________________________________________________________________________
void AliAnalysisTaskChargedJetsPA::GetKTBackground2Density(Int_t numberExcludeLeadingJets, Double_t& rhoMedian, Double_t& areaMean, Double_t etaMin, Double_t etaMax)
{
  #ifdef DEBUGMODE
    AliInfo("Getting KT background 2 density.");
  #endif

  // static declaration. Advantage: more speed. Disadvantage: Problematic for events with more than 1024 jets :)
  static Double_t tmpRhos[1024];
  static Double_t tmpAreas[1024];

  // Setting invalid values
  rhoMedian = -1.0;
  areaMean= -1.0;

  if ((etaMin == 0) && (etaMax == 0))
  {
    etaMin = -fBackgroundJetEtaWindow;
    etaMax = +fBackgroundJetEtaWindow;
  }

  Int_t jetCountAccepted = 0;
  Int_t jetCount = fBackgroundJetArray->GetEntries();

  for (Int_t i = 0; i < jetCount; i++)
  {
    Bool_t jetValid = kTRUE;
    AliEmcalJet* jet = static_cast<AliEmcalJet*>(fBackgroundJetArray->At(i));
    if (!jet)
    {
      AliError(Form("%s: Could not receive jet %d", GetName(), i));
      continue;
    } 

    if (!((jet->Eta() >= etaMin) && (jet->Eta() < etaMax)))
      continue;
    if (!IsBackgroundJetInAcceptance(jet))
      continue;

    // Look, if theres an overlap of leading jets/ kT jet. If yes, exclude this jet
    for(Int_t j=0;j<numberExcludeLeadingJets;j++)
    {
      AliEmcalJet* tmpLeadingJet = NULL;

      if (j==0)
        tmpLeadingJet = fFirstLeadingJet;
      else if (j==1)
        tmpLeadingJet = fSecondLeadingJet;
      else
        AliFatal("Trying to exclude more than 2 jets in KT background 2 -- not implemented.");

      if (tmpLeadingJet)
      {
        Double_t tmpDeltaPhi = GetDeltaPhi(jet->Phi(), tmpLeadingJet->Phi());
        if ( tmpDeltaPhi*tmpDeltaPhi + TMath::Abs(jet->Eta()-tmpLeadingJet->Eta())*TMath::Abs(jet->Eta()-tmpLeadingJet->Eta()) <= fBackgroundJetRadius*fBackgroundJetRadius)
        {
          jetValid = kFALSE;
          break;
        }
      }
    }

    if(!jetValid)
      continue;
   
    tmpRhos[jetCountAccepted] = jet->Pt() / jet->Area();
    tmpAreas[jetCountAccepted] = jet->Area();
    jetCountAccepted++;
  }

  if (jetCountAccepted > 0)
  {
    rhoMedian = TMath::Median(jetCountAccepted, tmpRhos);
    areaMean   = TMath::Mean(jetCountAccepted, tmpAreas);
  }
  #ifdef DEBUGMODE
    AliInfo("Got KT background 2 density.");
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
  rhoMean = -1.0;
  rhoMedian = -1.0;

  // Exclude UP TO numberExcludeLeadingJets
  if (fNumberSignalJets < 2)
    numberExcludeLeadingJets = fNumberSignalJets;

  // Search given amount of RCs
  Int_t numAcceptedRCs = 0;
  for(Int_t i=0;i<numberRandCones;i++)
  {
    Double_t tmpRandConeEta = 0.0;
    Double_t tmpRandConePhi = 0.0;
    Double_t excludedJetEta = 0.0;
    Double_t excludedJetPhi = 0.0;

    // Search random cone in acceptance with no overlap with already excluded jets (leading jets and random cones)
    Bool_t coneValid = kTRUE;

    // Set the random cone position
    if ((etaMin == 0) && (etaMax == 0))
      tmpRandConeEta = (fTrackEtaWindow-fRandConeRadius)*(2.0*fRandom->Rndm()-1.0); // full RC is in acceptance
    else
      tmpRandConeEta = etaMin + fRandom->Rndm()*(etaMax-etaMin);

    tmpRandConePhi = fRandom->Rndm()*TMath::TwoPi();

    // Go through all excluded leading jets and check if there's an overlap
     
    for(Int_t j=0;j<numberExcludeLeadingJets;j++)
    {
      AliEmcalJet* tmpJet = NULL;

      if (j==0)
        tmpJet = fFirstLeadingJet;
      else if (j==1)
        tmpJet = fSecondLeadingJet;
      else
        AliFatal("Trying to exclude more than 2 jets in RC background -- not implemented.");

      excludedJetPhi = tmpJet->Phi();
      excludedJetEta = tmpJet->Eta();
      Double_t tmpDeltaPhi = GetDeltaPhi(tmpRandConePhi, excludedJetPhi);
      
      if ( tmpDeltaPhi*tmpDeltaPhi + TMath::Abs(tmpRandConeEta-excludedJetEta)*TMath::Abs(tmpRandConeEta-excludedJetEta) <= fRandConeRadius*fRandConeRadius)
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
void AliAnalysisTaskChargedJetsPA::GetTrackBackgroundDensity(Int_t numberExcludeLeadingJets, Double_t& rhoMean, Double_t& area, Double_t etaMin, Double_t etaMax)
{
  #ifdef DEBUGMODE
    AliInfo("Getting track background density.");
  #endif

  Double_t summedTracksPt = 0.0;

  if ((etaMin == 0) && (etaMax == 0))
  {
    etaMin = -fTrackEtaWindow;
    etaMax = +fTrackEtaWindow;
  }

  // Setting invalid values
  rhoMean = -1.0;
  area = -1.0;
  // Exclude UP TO numberExcludeLeadingJets
  if (fNumberSignalJets < 2)
    numberExcludeLeadingJets = fNumberSignalJets;


  Int_t   trackCount = fTrackArray->GetEntries();
  Int_t   trackCountAccepted = 0;
  for (Int_t i = 0; i < trackCount; i++)
  {
    Bool_t  trackValid = kTRUE;
    AliVTrack* tmpTrack = static_cast<AliVTrack*>(fTrackArray->At(i));
    if (IsTrackInAcceptance(tmpTrack))
      if ((tmpTrack->Eta() >= etaMin) && (tmpTrack->Eta() < etaMax))
      {
        for (Int_t j = 0; j < numberExcludeLeadingJets; j++)
        {
          AliEmcalJet* tmpJet = NULL;
          if (j==0)
            tmpJet = fFirstLeadingJet;
          else if (j==1)
            tmpJet = fSecondLeadingJet;
          else
            AliFatal("Trying to exclude more than 2 jets in track background -- not implemented.");

          if (IsTrackInCone(tmpTrack, tmpJet->Eta(), tmpJet->Phi(), fTrackBackgroundConeRadius))
          {
            trackValid = kFALSE;
            break;
          }
        }
        if (trackValid)
        {
          // Add track pt to array
          summedTracksPt = summedTracksPt + tmpTrack->Pt();
          trackCountAccepted++;
        }
      }
  }

  if (trackCountAccepted > 0)
  {
    Double_t tmpArea = 0.0;

    tmpArea = (2.0*fTrackEtaWindow) * TMath::TwoPi() * (etaMax-etaMin)/(2.0*fTrackEtaWindow); // area of the used eta strip
    
    // Now: exclude the part of the leading jet that is in the strip
    if (numberExcludeLeadingJets == 2)
      tmpArea = tmpArea*(1.0-MCGetOverlapCircleRectancle(fFirstLeadingJet->Eta(), fFirstLeadingJet->Phi(), fTrackBackgroundConeRadius, etaMin, etaMax, 0., TMath::TwoPi()) -MCGetOverlapCircleRectancle(fSecondLeadingJet->Eta(), fSecondLeadingJet->Phi(), fTrackBackgroundConeRadius, etaMin, etaMax, 0., TMath::TwoPi()));
    else if (numberExcludeLeadingJets == 1)
      tmpArea = tmpArea*(1.0-MCGetOverlapCircleRectancle(fFirstLeadingJet->Eta(), fFirstLeadingJet->Phi(), fTrackBackgroundConeRadius, etaMin, etaMax, 0., TMath::TwoPi()));
   
    rhoMean = summedTracksPt/tmpArea;
    area  = tmpArea;
  }

  #ifdef DEBUGMODE
    AliInfo("Got track background density.");
  #endif
}

//________________________________________________________________________
void AliAnalysisTaskChargedJetsPA::GetTrackBackgroundDensity(Int_t numberExcludeLeadingJets, Double_t& rhoMean, Double_t& area, AliEmcalJet* excludeJet1, AliEmcalJet* excludeJet2, Bool_t doSearchPerpendicular)
{
  #ifdef DEBUGMODE
    AliInfo("Getting track background density.");
  #endif

  // Setting invalid values
  Double_t summedTracksPt = 0.0;
  rhoMean = -1.0;
  area = -1.0;

  Double_t tmpRadius = 0.0;
  if (doSearchPerpendicular)
    tmpRadius = 0.5*TMath::Pi(); // exclude 90 degrees around jets
  else
    tmpRadius = fSignalJetRadius;
    
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
    AliInfo("Got track background density.");
  #endif
}

//________________________________________________________________________
void AliAnalysisTaskChargedJetsPA::Calculate(AliVEvent* event)
{
  #ifdef DEBUGMODE
    AliInfo("Starting Calculate().");
  #endif
  ////////////////////// NOTE: initialization & casting

  if (!event) {
    AliError("??? Event pointer == 0 ???");
    return;
  }
 
  if (!fInitialized)
    ExecOnce(); // Get tracks, jets, background from arrays if not already given + Init Histos
  
  // Additional cuts
  FillHistogram("hNumberEvents", 0.5); // number of events before manual cuts

  if(!fHelperClass->IsVertexSelected2013pA(event))
    return;

  FillHistogram("hNumberEvents", 1.5); // number of events after manual cuts

  #ifdef DEBUGMODE
    AliInfo("Calculate()::Init done.");
  #endif

  ////////////////////// NOTE: Get Centrality, (Leading)Signal jets and Background

  // Get centrality (V0A)
  AliCentrality* tmpCentrality = NULL;
  tmpCentrality = event->GetCentrality();
  Double_t centralityPercentile = 0.0;
  if (tmpCentrality != NULL)
    centralityPercentile = tmpCentrality->GetCentralityPercentile("V0A");

  // Get jets
  if (fAnalyzeBackground || fAnalyzeJets)
    GetSignalJets();

  // Get background

  // Background with N excluded leading jets
  std::vector<Double_t>  ktBackgroundRhoMedian(fBackgroundEtaBins+1);
  std::vector<Double_t>  ktBackgroundAreaMean(fBackgroundEtaBins+1);
  std::vector<Double_t>  ktBackground2RhoMedian(fBackgroundEtaBins+1);
  std::vector<Double_t>  ktBackground2AreaMean(fBackgroundEtaBins+1);
  std::vector<Double_t>  rcBackgroundRhoMean(fBackgroundEtaBins+1);
  std::vector<Double_t>  rcBackgroundRhoMedian(fBackgroundEtaBins+1);
  std::vector<Double_t>  trackBackgroundRhoMean(fBackgroundEtaBins+1);
  std::vector<Double_t>  trackBackgroundArea(fBackgroundEtaBins+1);
  Double_t  dijetBackground = -1.0; // calculation only done in events with dijets I!
  Double_t  dijetBackgroundPerpendicular = -1.0; // calculation only done in events with dijets I!
  if (fAnalyzeBackground)
  {

    // Get backgrounds in bins of eta
    for(Int_t i = 0; i<fBackgroundEtaBins; i++)
    {
      // scheme: etaMin = RangeMin + l*binN; etaMax = RangeMin + l*(binN+1)

      Double_t etaMin = -(fTrackEtaWindow-fRandConeRadius) + 2*(fTrackEtaWindow-fRandConeRadius)/fBackgroundEtaBins *  i;
      Double_t etaMax = -(fTrackEtaWindow-fRandConeRadius) + 2*(fTrackEtaWindow-fRandConeRadius)/fBackgroundEtaBins * (i+1);
      GetRCBackgroundDensity (fNumberExcludedJets, rcBackgroundRhoMean[i],rcBackgroundRhoMedian[i],  etaMin, etaMax);
      GetTrackBackgroundDensity (fNumberExcludedJets, trackBackgroundRhoMean[i], trackBackgroundArea[i], etaMin, etaMax);
      GetKTBackgroundDensity (fNumberExcludedJets, ktBackgroundRhoMedian[i], ktBackgroundAreaMean[i], etaMin, etaMax);
      GetKTBackground2Density (fNumberExcludedJets, ktBackground2RhoMedian[i], ktBackground2AreaMean[i], etaMin, etaMax);

    }
    Int_t tmpNRCs = 0;

    // All eta in one bin
    tmpNRCs = GetRCBackgroundDensity (fNumberExcludedJets, rcBackgroundRhoMean[fBackgroundEtaBins], rcBackgroundRhoMedian[fBackgroundEtaBins]);
    FillHistogram("hAccConesInRCBackground", tmpNRCs);
    GetTrackBackgroundDensity (fNumberExcludedJets, trackBackgroundRhoMean[fBackgroundEtaBins], trackBackgroundArea[fBackgroundEtaBins]);
    GetKTBackgroundDensity (fNumberExcludedJets, ktBackgroundRhoMedian[fBackgroundEtaBins], ktBackgroundAreaMean[fBackgroundEtaBins]);
    GetKTBackground2Density (fNumberExcludedJets, ktBackground2RhoMedian[fBackgroundEtaBins], ktBackground2AreaMean[fBackgroundEtaBins]);

  }

  #ifdef DEBUGMODE
    AliInfo("Calculate()::Centrality&SignalJets&Background-Calculation done.");
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

  ////////////////////// NOTE: Track & QA histograms

  if (fAnalyzeQA)
  {
    FillHistogram("hVertexZ",event->GetPrimaryVertex()->GetZ());
    FillHistogram("hVertexR",TMath::Sqrt(event->GetPrimaryVertex()->GetX()*event->GetPrimaryVertex()->GetX() + event->GetPrimaryVertex()->GetY()*event->GetPrimaryVertex()->GetY()));
    FillHistogram("hCentrality",centralityPercentile);

    Int_t trackCountAcc = 0;
    Int_t nTracks = fTrackArray->GetEntries();
    for (Int_t i = 0; i < nTracks; i++)
    {
      AliVTrack* track = static_cast<AliVTrack*>(fTrackArray->At(i));
      if (IsTrackInAcceptance(track))
      {
        FillHistogram("hTrackPhiEta", track->Phi(),track->Eta(), 1);
        FillHistogram("hTrackPt", track->Pt());
        FillHistogram("hTrackEta", track->Eta());
        FillHistogram("hTrackCharge", track->Charge());
        trackCountAcc++;
      }
    }
    FillHistogram("hTrackCountAcc", trackCountAcc, centralityPercentile);

    if (fHasClusters)
    {
      Int_t clusterCountAcc = 0;
      Int_t nClusters = fClusterArray->GetEntries();
      for (Int_t i = 0; i < nClusters; i++)
      {
        AliVCluster* cluster = static_cast<AliVCluster*>(fClusterArray->At(i));
        if (IsClusterInAcceptance(cluster))
        {
          FillHistogram("hClusterE", cluster->E());
          clusterCountAcc++;
        }
      }
      FillHistogram("hClusterCountAcc", clusterCountAcc, centralityPercentile);
    }
  }
  #ifdef DEBUGMODE
    AliInfo("Calculate()::QA done.");
  #endif

  ////////////////////// NOTE: Jet analysis and calculations

  if (fAnalyzeJets)
  {
    FillHistogram("hJetCountAll", fJetArray->GetEntries());
    FillHistogram("hJetCountAccepted", fNumberSignalJets);
    if (fFirstLeadingJet)
      FillHistogram("hLeadingJetPt", fFirstLeadingJet->Pt());
    if (fSecondLeadingJet)
      FillHistogram("hSecondLeadingJetPt", fSecondLeadingJet->Pt());

    // ### Dijets I ###
    if(fNumberSignalJets >= 2)
    {
      FillHistogram("hLeadingJetDeltaPhi", GetDeltaPhi(fFirstLeadingJet->Phi(), fSecondLeadingJet->Phi()));
      FillHistogram("hLeadingJetDeltaPhiPt", GetDeltaPhi(fFirstLeadingJet->Phi(), fSecondLeadingJet->Phi()), fFirstLeadingJet->Pt());

      if (IsDijet(fFirstLeadingJet, fSecondLeadingJet)) // Gettin' the money
      {
        FillHistogram("hDijetConstituentsPt", fFirstLeadingJet->Pt()); FillHistogram("hDijetConstituentsPt", fSecondLeadingJet->Pt());
        FillHistogram("hDijetLeadingJetPt", fFirstLeadingJet->Pt());
        FillHistogram("hDijetPtCorrelation", fFirstLeadingJet->Pt(), fSecondLeadingJet->Pt());
        Double_t dummyArea = 0;
        GetTrackBackgroundDensity (2, dijetBackground, dummyArea, fFirstLeadingJet, fSecondLeadingJet, kFALSE);
        GetTrackBackgroundDensity (2, dijetBackgroundPerpendicular, dummyArea, fFirstLeadingJet, fSecondLeadingJet, kTRUE);
      }
    }

    // SIGNAL JET ANALYSIS
    for (Int_t i = 0; i<fNumberSignalJets; i++)
    {
      AliEmcalJet* tmpJet = fSignalJets[i];

      FillHistogram("hJetPtArea", tmpJet->Pt(), tmpJet->Area());
      FillHistogram("hJetPtEta", tmpJet->Pt(), tmpJet->Eta());
      FillHistogram("hJetPtPhi", tmpJet->Pt(), tmpJet->Phi());
      FillHistogram("hJetPtCentrality", tmpJet->Pt(), centralityPercentile);
      FillHistogram("hJetArea", tmpJet->Area());
      FillHistogram("hJetPt", tmpJet->Pt());
      FillHistogram("hJetPhiEta", tmpJet->Phi(),tmpJet->Eta());
      
      // Background subtracted spectra
      
      FillHistogram("hJetPtBgrdSubtractedRC", GetCorrectedJetPt(tmpJet, rcBackgroundRhoMean[fBackgroundEtaBins]));
      FillHistogram("hJetPtBgrdSubtractedKT", GetCorrectedJetPt(tmpJet, ktBackgroundRhoMedian[fBackgroundEtaBins], kTRUE));
      FillHistogram("hJetPtBgrdSubtractedKT2", GetCorrectedJetPt(tmpJet, ktBackground2RhoMedian[fBackgroundEtaBins], kTRUE));
      FillHistogram("hJetPtBgrdSubtractedKTNoEtaCorr", GetCorrectedJetPt(tmpJet, ktBackgroundRhoMedian[fBackgroundEtaBins]));
      FillHistogram("hJetPtBgrdSubtractedKT2NoEtaCorr", GetCorrectedJetPt(tmpJet, ktBackground2RhoMedian[fBackgroundEtaBins]));
      FillHistogram("hJetPtBgrdSubtractedTR", GetCorrectedJetPt(tmpJet, trackBackgroundRhoMean[fBackgroundEtaBins]));


      Double_t tmpCorrFactor = GetJetBackgroundCorrFactor(tmpJet->Eta(), ktBackgroundRhoMedian[fBackgroundEtaBins]);
      FillHistogram("hAppliedEtaCorrectionFactor", tmpCorrFactor);
      tmpCorrFactor = GetJetBackgroundCorrFactor(tmpJet->Eta(), ktBackground2RhoMedian[fBackgroundEtaBins]);
      FillHistogram("hAppliedEtaCorrectionFactor2", tmpCorrFactor);

      // Signal jet vs. signal jet
      for (Int_t j = i+1; j<fNumberSignalJets; j++)
      {
        AliEmcalJet* tmpJet2 = fSignalJets[j];
        FillHistogram("hJetDeltaPhi", GetDeltaPhi(tmpJet->Phi(), tmpJet2->Phi()));
        FillHistogram("hJetDeltaPhiPt", GetDeltaPhi(tmpJet->Phi(), tmpJet2->Phi()), max(tmpJet->Pt(), tmpJet2->Pt()));

        // ### Dijets II ###
        if (IsDijet(tmpJet, tmpJet2)) // Gettin' the money
        {
          FillHistogram("hDijet2ConstituentsPt", tmpJet->Pt()); FillHistogram("hDijet2ConstituentsPt", tmpJet2->Pt());
          FillHistogram("hDijet2LeadingJetPt", fFirstLeadingJet->Pt());
          FillHistogram("hDijet2PtCorrelation", tmpJet->Pt(), tmpJet2->Pt());
        }
      }
    }
  } //endif AnalyzeJets

  #ifdef DEBUGMODE
    AliInfo("Calculate()::Jets done.");
  #endif
  ////////////////////// NOTE: Background analysis

  if (fAnalyzeBackground)
  {

    Int_t leadingJetIds[] = {-1, -1};
    GetLeadingJets(fBackgroundJetArray, &leadingJetIds[0], kFALSE);

    for (Int_t i = 0; i < fBackgroundJetArray->GetEntries(); i++)
    {
      AliEmcalJet* jet = static_cast<AliEmcalJet*>(fBackgroundJetArray->At(i));
      if (!jet)
      {
        AliError(Form("%s: Could not receive kt jet %d", GetName(), i));
        continue;
      }
      if (!IsBackgroundJetInAcceptance(jet))
        continue;
      if (!((jet->Eta() >= -fBackgroundJetEtaWindow) && (jet->Eta() < fBackgroundJetEtaWindow)))
        continue;
      
      FillHistogram("hKTJetPhiEta", jet->Phi(),jet->Eta());
      if(i==leadingJetIds[0])
        FillHistogram("hKTLeadingJetPhiEta", jet->Phi(),jet->Eta());
            
    }

    // ############# RC, Track, and KT background calculations
    Double_t etaMin = 0;
    for (Int_t i=0;i<fBackgroundEtaBins;i++)
    {
      etaMin = -(fTrackEtaWindow-fRandConeRadius) + 2*(fTrackEtaWindow-fRandConeRadius)/fBackgroundEtaBins *  (i+0.5);
      FillHistogram("hRCBackground", etaMin, rcBackgroundRhoMean[i]);
      FillHistogram("hTrackBackground", etaMin, trackBackgroundRhoMean[i]);
      FillHistogram("hKTBackground", etaMin, ktBackgroundRhoMedian[i]);
      FillHistogram("hKTBackground2", etaMin, ktBackground2RhoMedian[i]);
      if(centralityPercentile <= 20.)
      {
        FillHistogram("hRCBackgroundMostCentral", etaMin, rcBackgroundRhoMean[i]);
        FillHistogram("hTrackBackgroundMostCentral", etaMin, trackBackgroundRhoMean[i]);
        FillHistogram("hKTBackgroundMostCentral", etaMin, ktBackgroundRhoMedian[i]);
        FillHistogram("hKTBackground2MostCentral", etaMin, ktBackground2RhoMedian[i]);
      }
      else if(centralityPercentile >= 80.)
      {
        FillHistogram("hRCBackgroundMostPeripheral", etaMin, rcBackgroundRhoMean[i]);
        FillHistogram("hTrackBackgroundMostPeripheral", etaMin, trackBackgroundRhoMean[i]);
        FillHistogram("hKTBackgroundMostPeripheral", etaMin, ktBackgroundRhoMedian[i]);
        FillHistogram("hKTBackground2MostPeripheral", etaMin, ktBackground2RhoMedian[i]);
      }
    }

    FillHistogram("hRCBackgroundVsCentrality", rcBackgroundRhoMean[fBackgroundEtaBins], centralityPercentile);
    FillHistogram("hTrackBackgroundVsCentrality", trackBackgroundRhoMean[fBackgroundEtaBins], centralityPercentile);
    FillHistogram("hKTBackgroundVsCentrality", ktBackgroundRhoMedian[fBackgroundEtaBins], centralityPercentile);
    FillHistogram("hKTBackground2VsCentrality", ktBackground2RhoMedian[fBackgroundEtaBins], centralityPercentile);

    if (dijetBackground >= 0)
    {
      // Background in Dijet events
      FillHistogram("hDijetBackground", dijetBackground); 
      if(centralityPercentile <= 20.)
        FillHistogram("hDijetBackgroundMostCentral", dijetBackground); 
      FillHistogram("hDijetBackgroundVsCentrality", dijetBackground, centralityPercentile);
    }
    if (dijetBackgroundPerpendicular >= 0)
    {
      // Background in Dijet events
      FillHistogram("hDijetBackgroundPerpendicular", dijetBackgroundPerpendicular); 
      if(centralityPercentile <= 20.)
        FillHistogram("hDijetBackgroundPerpendicularMostCentral", dijetBackgroundPerpendicular); 
      FillHistogram("hDijetBackgroundPerpendicularVsCentrality", dijetBackgroundPerpendicular, centralityPercentile);
    }

    // ########## Delta pT calculations (most central, kt is eta corrected)
    if (centralityPercentile <= 20.)
    {
      Double_t tmpDeltaPtKT, tmpDeltaPtKT2Excl, tmpDeltaPtKT1Excl;
      Double_t tmpDeltaPtKTEta, tmpDeltaPtKTEta2Excl, tmpDeltaPtKTEta1Excl, tmpDeltaPtKT2Eta2Excl;
      Double_t tmpDeltaPtRC, tmpDeltaPtRC2Excl, tmpDeltaPtRC1Excl;
      Double_t tmpDeltaPtTR, tmpDeltaPtTR2Excl, tmpDeltaPtTR1Excl;

      GetDeltaPt(tmpDeltaPtKT, ktBackgroundRhoMedian[fBackgroundEtaBins], 0, -1, kTRUE);
      GetDeltaPt(tmpDeltaPtKTEta, ktBackgroundRhoMedian[fKTDeltaPtEtaBin], 0, fKTDeltaPtEtaBin);
      GetDeltaPt(tmpDeltaPtRC, rcBackgroundRhoMean[fBackgroundEtaBins], 0);
      GetDeltaPt(tmpDeltaPtTR, trackBackgroundRhoMean[fBackgroundEtaBins], 0);

      GetDeltaPt(tmpDeltaPtKT1Excl, ktBackgroundRhoMedian[fBackgroundEtaBins], 1, -1, kTRUE);
      GetDeltaPt(tmpDeltaPtKTEta1Excl, ktBackgroundRhoMedian[fKTDeltaPtEtaBin], 1, fKTDeltaPtEtaBin);
      GetDeltaPt(tmpDeltaPtRC1Excl, rcBackgroundRhoMean[fBackgroundEtaBins], 1);
      GetDeltaPt(tmpDeltaPtTR1Excl, trackBackgroundRhoMean[fBackgroundEtaBins], 1);

      GetDeltaPt(tmpDeltaPtKT2Excl, ktBackgroundRhoMedian[fBackgroundEtaBins], 2, -1, kTRUE);
      GetDeltaPt(tmpDeltaPtKTEta2Excl, ktBackgroundRhoMedian[fKTDeltaPtEtaBin], 2, fKTDeltaPtEtaBin);
      GetDeltaPt(tmpDeltaPtRC2Excl, rcBackgroundRhoMean[fBackgroundEtaBins], 2);
      GetDeltaPt(tmpDeltaPtTR2Excl, trackBackgroundRhoMean[fBackgroundEtaBins], 2);

      GetDeltaPt(tmpDeltaPtKT2Eta2Excl, ktBackground2RhoMedian[fKTDeltaPtEtaBin], 2, fKTDeltaPtEtaBin);

      // kT Background
      if(tmpDeltaPtKT > -10000.0)
        FillHistogram("hDeltaPtKT", tmpDeltaPtKT);
      if(tmpDeltaPtKT1Excl > -10000.0)
        FillHistogram("hDeltaPtKT1Excl", tmpDeltaPtKT1Excl);
      if(tmpDeltaPtKT2Excl > -10000.0)
        FillHistogram("hDeltaPtKT2Excl", tmpDeltaPtKT2Excl);

      if(tmpDeltaPtKT > -10000.0)
        FillHistogram("hDeltaPtKTEta", tmpDeltaPtKTEta);
      if(tmpDeltaPtKTEta1Excl > -10000.0)
        FillHistogram("hDeltaPtKTEta1Excl", tmpDeltaPtKTEta1Excl);
      if(tmpDeltaPtKTEta2Excl > -10000.0)
        FillHistogram("hDeltaPtKTEta2Excl", tmpDeltaPtKTEta2Excl);
      if(tmpDeltaPtKT2Eta2Excl > -10000.0)
        FillHistogram("hDeltaPtKT2Eta2Excl", tmpDeltaPtKT2Eta2Excl);

      // RC Background
      if(tmpDeltaPtRC > -10000.0)
        FillHistogram("hDeltaPtRC", tmpDeltaPtRC);
      if(tmpDeltaPtRC1Excl > -10000.0)
        FillHistogram("hDeltaPtRC1Excl", tmpDeltaPtRC1Excl);
      if(tmpDeltaPtRC2Excl > -10000.0)
        FillHistogram("hDeltaPtRC2Excl", tmpDeltaPtRC2Excl);
      // TR Background
      if(tmpDeltaPtTR > -10000.0)
        FillHistogram("hDeltaPtTR", tmpDeltaPtTR);
      if(tmpDeltaPtTR1Excl > -10000.0)
        FillHistogram("hDeltaPtTR1Excl", tmpDeltaPtTR1Excl);
      if(tmpDeltaPtTR2Excl > -10000.0)
        FillHistogram("hDeltaPtTR2Excl", tmpDeltaPtTR2Excl);
    }
  }
  
  #ifdef DEBUGMODE
    AliInfo("Calculate()::Background done.");
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
inline void AliAnalysisTaskChargedJetsPA::FillHistogram(const char * key, Double_t x)
{
  TH1* tmpHist = static_cast<TH1*>(fOutputList->FindObject(GetHistoName(key)));
  if(!tmpHist)
  {
    AliInfo(Form("Cannot find histogram <%s> ",key)) ;
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
    AliInfo(Form("Cannot find histogram <%s> ",key));
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
    AliInfo(Form("Cannot find histogram <%s> ",key));
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

  Calculate(InputEvent());
        
  PostData(1, fOutputList);
}
