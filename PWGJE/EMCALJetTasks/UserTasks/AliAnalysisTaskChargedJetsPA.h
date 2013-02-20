#ifndef ALIANALYSISTASKCHARGEDJETSPA_H
#define ALIANALYSISTASKCHARGEDJETSPA_H

//#define DEBUGMODE

class TH1F;
class TH2F;
class TList;
class TClonesArray;
class TString;
class AliEmcalJet;
class AliRhoParameter;
class AliVParticle;
class AliLog;
class AliAnalysisUtils;

class AliAnalysisTaskChargedJetsPA : public AliAnalysisTaskSE {
 public:
  // ######### CONTRUCTORS/DESTRUCTORS AND STD FUNCTIONS
  AliAnalysisTaskChargedJetsPA(const char *name, const char* trackArrayName, const char* jetArrayName, const char* backgroundJetArrayName);
  AliAnalysisTaskChargedJetsPA();
  virtual ~AliAnalysisTaskChargedJetsPA();
  virtual void     UserCreateOutputObjects();
  virtual void     UserExec(Option_t *option);
  virtual void     Terminate(Option_t *);

  // ######### SETTERS/GETTERS
  void        SetAnalyzeJets(Bool_t val) {fAnalyzeJets = val;}
  void        SetAnalyzeBackground(Bool_t val) {fAnalyzeBackground = val;}
  void        SetAnalyzePythia(Bool_t val) {fAnalyzePythia = val;}

  void        SetTrackMinPt(Double_t minPt) {fMinJetPt = minPt;}
  void        SetSignalJetMinPt(Double_t minPt) {fMinJetPt = minPt;}
  void        SetSignalJetMinArea(Double_t minArea) {fMinJetArea = minArea;}
  void        SetBackgroundJetMinPt(Double_t minPt) {fMinBackgroundJetPt = minPt;}
  void        SetDijetLeadingMinPt(Double_t minPt) {fMinDijetLeadingPt = minPt;}
  void        SetNumberOfPtHardBins(Int_t count) {fNumPtHardBins = count;}
  void        SetNumberOfRandConesPerEvent(Int_t count) {fNumberRandCones = count;}
  void        SetRandConeRadius(Double_t radius) {fRandConeRadius = radius;}
  void        SetSignalJetRadius(Double_t radius) {fSignalJetRadius = radius;}
  void        SetBackgroundJetRadius(Double_t radius) {fBackgroundJetRadius = radius;}
  void        SetTRBackgroundConeRadius(Double_t radius) {fTRBackgroundConeRadius = radius;}
  void        SetCentralityType(const char* type) {fCentralityType = type;}

  void        SetDijetMaxAngleDeviation(Double_t degrees) {fDijetMaxAngleDeviation = degrees/360.0 * TMath::TwoPi();} // degrees are more comfortable
  void        SetAcceptanceWindows(Double_t trackEta, Double_t vertexZ, Double_t vertexMaxR, Double_t signalJetRadius, Double_t bgrdJetRadius){fVertexWindow = vertexZ; fVertexMaxR = vertexMaxR; fTrackEtaWindow = trackEta; fSignalJetRadius = signalJetRadius; fBackgroundJetRadius = bgrdJetRadius; fSignalJetEtaWindow = fTrackEtaWindow-fSignalJetRadius; fBackgroundJetEtaWindow = fTrackEtaWindow-fBackgroundJetRadius;}
  void        SetKTEtaCorrectionFactors(TH2D* histo);
  void        SetRCEtaCorrectionFactors(TH2D* histo);
  void        SetTREtaCorrectionFactors(TH2D* histo);
  Int_t       GetInstanceCounter() {return fTaskInstanceCounter;}

 private:
  enum EtaCorrectionMode {kNoEtaCorrection, kKTEtaCorrection, kRCEtaCorrection, kTREtaCorrection};

  // ######### MAIN CALCULATION FUNCTIONS
  void        GetSignalJets();
  Int_t       GetLeadingJets(TClonesArray* jetArray, Int_t* jetIDArray, Bool_t isSignalJets);
  Double_t    GetBackgroundEtaCorrFactor(EtaCorrectionMode mode, Double_t eta, Double_t background);
  Double_t    GetBackgroundEtaBinCorrFactor(EtaCorrectionMode mode, Int_t eta, Double_t background);
  Double_t    GetCorrectedJetPt(AliEmcalJet* jet, Double_t background, EtaCorrectionMode mode);
  void        GetDeltaPt(Double_t& deltaPt, Double_t rho, EtaCorrectionMode mode, Bool_t leadingJetExclusion = kTRUE);

  void        GetKTBackgroundDensity(Int_t numberExcludeLeadingJets, Double_t& rhoMedian, Double_t& areaMean, Double_t etaMin = 0, Double_t etaMax = 0);
  Int_t       GetRCBackgroundDensity(Int_t numberExcludeLeadingJets, Double_t& rhoMean, Double_t& rhoMedian, Double_t etaMin = 0, Double_t etaMax = 0, Int_t numberRandCones = 0);
  void        GetTRBackgroundDensity(Int_t numberExcludeLeadingJets, Double_t& rhoMean, Double_t& area, Double_t etaMin = 0, Double_t etaMax = 0);
  void        GetTRBackgroundDensity(Int_t numberExcludeLeadingJets, Double_t& rhoMean, Double_t& area, AliEmcalJet* excludeJet1, AliEmcalJet* excludeJet2, Bool_t doSearchPerpendicular);
  Double_t    GetConePt(Double_t eta, Double_t phi, Double_t radius);
  Double_t    GetPtHard();
  Int_t       GetPtHardBin();
  void        GetPerpendicularCone(Double_t vecPhi, Double_t vecTheta, Double_t& conePt);

  // ######### CHECK FUNCTIONS
  Bool_t      IsTrackInAcceptance(AliVParticle* track);
  Bool_t      IsTrackInCone(AliVTrack* track, Double_t eta, Double_t phi, Double_t radius);

  Bool_t      IsBackgroundJetInAcceptance(AliEmcalJet* jet);
  Bool_t      IsSignalJetInAcceptance(AliEmcalJet* jet);
  Bool_t      IsDijet(AliEmcalJet* jet1, AliEmcalJet* jet2);
  
  // ######### HELPER FUNCTIONS
  Double_t    EtaToTheta(Double_t arg);
  Double_t    ThetaToEta(Double_t arg);
  Double_t    GetDeltaPhi(Double_t phi1, Double_t phi2);
  Double_t    MCGetOverlapCircleRectancle(Double_t cPosX, Double_t cPosY, Double_t cRadius, Double_t rPosXmin, Double_t rPosXmax, Double_t rPosYmin, Double_t rPosYmax);

  // ######### HISTOGRAM FUNCTIONS
  void        FillHistogram(const char * key, Double_t x);
  void        FillHistogram(const char * key, Double_t x, Double_t y);
  void        FillHistogram(const char * key, Double_t x, Double_t y, Double_t add);
  const char* GetHistoName(const char* name)
  {
    if (fIsMC)    
      return Form("H%d_%s_MC", fTaskInstanceCounter, name);
    return Form("H%d_%s", fTaskInstanceCounter, name);
  }
  template <class T> T* AddHistogram1D(const char* name = "CustomHistogram", const char* title = "NO_TITLE", const char* options = "", Int_t xBins = 100, Double_t xMin = 0.0, Double_t xMax = 20.0, const char* xTitle = "x axis", const char* yTitle = "y axis");
  template <class T> T* AddHistogram2D(const char* name = "CustomHistogram", const char* title = "NO_TITLE", const char* options = "", Int_t xBins = 100, Double_t xMin = 0.0, Double_t xMax = 20.0, Int_t yBins = 100, Double_t yMin = 0.0, Double_t yMax = 20.0, const char* xTitle = "x axis", const char* yTitle = "y axis", const char* zTitle = "z axis");

  // ######### STANDARD FUNCTIONS
  Bool_t    Notify();
  void      Calculate(AliVEvent* event);
  void      ExecOnce();
  void      Init ();

  TList*              fOutputList;            //! Output list
  // ########## USAGE TRIGGERS 
  Bool_t              fAnalyzeJets;           // trigger if jets should be processed
  Bool_t              fAnalyzeBackground;     // trigger if background should be processed
  Bool_t              fAnalyzePythia;         // trigger if pythia properties should be processed
  Bool_t              fHasTracks;             // trigger if tracks are actually valid
  Bool_t              fHasJets;               // trigger if jets are actually valid
  Bool_t              fHasBackgroundJets;     // trigger if background is actually valid
  Bool_t              fIsMC;                  // trigger if data is MC (for naming reasons)

  // ########## SOURCE INFORMATION
  TClonesArray*       fJetArray;              //! object containing the jets
  TClonesArray*       fTrackArray;            //! object containing the tracks
  TClonesArray*       fBackgroundJetArray;    //! object containing background jets
  TString*            fJetArrayName;          // name of object containing the jets
  TString*            fTrackArrayName;        // name of object containing the tracks
  TString*            fBackgroundJetArrayName;// name of object containing event wise bckgrds
  Int_t               fNumPtHardBins;         // Number of used pt hard bins

  // ########## JET/DIJET/RC PROPERTIES
  Double_t            fRandConeRadius;        // Radius for the random cones
  Double_t            fSignalJetRadius;       // Radius for the signal jets
  Double_t            fBackgroundJetRadius;   // Radius for the KT background jets
  Double_t            fTRBackgroundConeRadius;// Radius for the jets excluded in track background
  Int_t               fNumberRandCones;       // Number of random cones to be put into one event
  Int_t               fNumberExcludedJets;    // Number of jets to be excluded from backgrounds
  Double_t            fDijetMaxAngleDeviation;// Max angle deviation from pi between two jets to be accept. as dijet
  TH2D*               fJetKTEtaCorrection;    // Correction factors in bins of rho and eta to correct the eta dependence of the jet background
  TH2D*               fJetRCEtaCorrection;    // Correction factors in bins of rho and eta to correct the eta dependence of the jet background
  TH2D*               fJetTREtaCorrection;    // Correction factors in bins of rho and eta to correct the eta dependence of the jet background

  // ########## CUTS 
  Double_t            fSignalJetEtaWindow;    // +- window in eta for signal jets
  Double_t            fBackgroundJetEtaWindow;// +- window in eta for background jets
  Double_t            fTrackEtaWindow;        // +- window in eta for tracks
  Double_t            fVertexWindow;          // +- window in Z for the vertex
  Double_t            fVertexMaxR;            // +- window in R for the vertex (distance in xy-plane)
  Double_t            fMinTrackPt;            // Min track pt to be accepted
  Double_t            fMinJetPt;              // Min jet pt to be accepted
  Double_t            fMinJetArea;            // Min jet area to be accepted
  Double_t            fMinBackgroundJetPt;    // Min jet pt to be accepted as background jet
  Double_t            fMinDijetLeadingPt;     // Min jet pt to be accepted as constituent of dijet 
  TString             fCentralityType;        // Used centrality estimate (V0A, V0C, V0M, ...)

  // ########## EVENT PROPERTIES
  AliEmcalJet*        fFirstLeadingJet;       //! leading jet in event
  AliEmcalJet*        fSecondLeadingJet;      //! next to leading jet in event
  Int_t               fNumberSignalJets;      //! Number of signal jets in event
  AliEmcalJet*        fSignalJets[1024];      //! memory for signal jet pointers
  Double_t            fCrossSection;          //! value is filled, if pythia header is accessible
  Double_t            fTrials;                //! value is filled, if pythia header is accessible

  // ########## GENERAL VARS
  TRandom3*           fRandom;                //! A random number
  AliAnalysisUtils*   fHelperClass;           //! Vertex selection helper
  Bool_t              fInitialized;           //! trigger if tracks/jets are loaded
  Int_t               fTaskInstanceCounter;   // for naming reasons
  TList*              fHistList;              // Histogram list
  Int_t               fHistCount;             // Histogram count

  AliAnalysisTaskChargedJetsPA(const AliAnalysisTaskChargedJetsPA&);
  AliAnalysisTaskChargedJetsPA& operator=(const AliAnalysisTaskChargedJetsPA&);

  ClassDef(AliAnalysisTaskChargedJetsPA, 2); // Charged jet analysis for pA

};
#endif
