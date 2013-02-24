#ifndef ALIANALYSISTASKQUALITYASSURANCEPA_H
#define ALIANALYSISTASKQUALITYASSURANCEPA_H

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

class AliAnalysisTaskQualityAssurancePA : public AliAnalysisTaskSE 
{
 public:

  AliAnalysisTaskQualityAssurancePA() : AliAnalysisTaskSE(), fOutputList(0), fAnalyzeQA(1), fAnalyzeJets(1), fAnalyzePythia(0), fHasTracks(0), fHasClusters(0), fHasJets(0), fIsMC(0), fJetArray(0), fTrackArray(0), fClusterArray(0), fJetArrayName(0), fTrackArrayName(0), fClusterArrayName(0), fRunNumbers(0), fNumPtHardBins(11), fSignalJetRadius(0.4), fNumberExcludedJets(2), fSignalJetEtaWindow(0.5), fTrackEtaWindow(0.9), fClusterEtaWindow(0.7), fVertexWindow(10.0), fVertexMaxR(1.0), fMinTrackPt(0.150), fMinClusterPt(0.300), fMinJetPt(1.0), fMinJetArea(0.4), fFirstLeadingJet(0), fSecondLeadingJet(0), fNumberSignalJets(0), fCrossSection(0.0), fTrials(0.0), fRandom(0), fHelperClass(0), fInitialized(0), fTaskInstanceCounter(0), fHistList(0), fHistCount(0) {}

  AliAnalysisTaskQualityAssurancePA(const char *name, const char* trackArrayName, const char* clusterArrayName, const char* jetArrayName);

  // Standard  functions
  virtual           ~AliAnalysisTaskQualityAssurancePA();
  virtual void      UserCreateOutputObjects();
  virtual void      UserExec(Option_t *option);
  virtual void      Terminate(Option_t *);
  
  // Setters
  void              SetAnalyzeTracks(Bool_t val) {fAnalyzeQA = val;}
  void              SetAnalyzeJets(Bool_t val) {fAnalyzeJets = val;}
  void              SetAnalyzePythia(Bool_t val) {fAnalyzePythia = val;}

  void              SetTrackMinPt(Double_t minPt) {fMinJetPt = minPt;}
  void              SetSignalJetMinPt(Double_t minPt) {fMinJetPt = minPt;}
  void              SetSignalJetMinArea(Double_t minArea) {fMinJetArea = minArea;}
  void              SetRunNumbers(const char* runNumbers) {*fRunNumbers = runNumbers;}
  void              SetNumberOfPtHardBins(Int_t count) {fNumPtHardBins = count;}

  void              SetSignalJetRadius(Double_t radius) {fSignalJetRadius = radius;}
  void              SetAcceptanceWindows(Double_t trackEta, Double_t vertexZ, Double_t vertexMaxR, Double_t signalJetRadius) {fVertexWindow = vertexZ; fVertexMaxR = vertexMaxR; fTrackEtaWindow = trackEta; fSignalJetRadius = signalJetRadius; fSignalJetEtaWindow = fTrackEtaWindow-fSignalJetRadius;}

  // Getters
  Int_t             GetInstanceCounter() {return fTaskInstanceCounter;}

 private:
  // ######### MAIN CALCULATION FUNCTIONS
  void              GetSignalJets();
  Int_t             GetLeadingJets(TClonesArray* jetArray, Int_t* jetIDArray);
  Double_t          GetPtHard();
  Int_t             GetPtHardBin();

  // ######### CHECK FUNCTIONS
  Bool_t            IsTrackInAcceptance(AliVParticle* track);
  Bool_t            IsClusterInAcceptance(AliVCluster* cluster);
  Bool_t            IsSignalJetInAcceptance(AliEmcalJet* jet);

  // ######### HISTOGRAM FUNCTIONS
  void              FillHistogram(const char* runNumber, const char * key, Double_t x);
  void              FillHistogram(const char* runNumber, const char * key, Double_t x, Double_t y);
  void              FillHistogram(const char* runNumber, const char * key, Double_t x, Double_t y, Double_t add);
  const char*       GetHistoName(const char* runNumber, const char* name)
  {
    if (fIsMC)    
      return Form("H%d_%s_%s_MC", fTaskInstanceCounter, name, runNumber);
    return Form("H%d_%s_%s", fTaskInstanceCounter, name, runNumber);
  }
  template <class T> T* AddHistogram1D(const char* runNumber, const char* name = "CustomHistogram", const char* title = "NO_TITLE", const char* options = "", Int_t xBins = 100, Double_t xMin = 0.0, Double_t xMax = 20.0, const char* xTitle = "x axis", const char* yTitle = "y axis");
  template <class T> T* AddHistogram2D(const char* runNumber, const char* name = "CustomHistogram", const char* title = "NO_TITLE", const char* options = "", Int_t xBins = 100, Double_t xMin = 0.0, Double_t xMax = 20.0, Int_t yBins = 100, Double_t yMin = 0.0, Double_t yMax = 20.0, const char* xTitle = "x axis", const char* yTitle = "y axis", const char* zTitle = "z axis");


  // ######### STANDARD FUNCTIONS
  Bool_t            Notify();
  void              Calculate(AliVEvent* event);
  void              ExecOnce();
  void              Init ();

  TList*              fOutputList;            //! Output list

  // ########## USAGE TRIGGERS 
  Bool_t              fAnalyzeQA;             // trigger if tracks should be analyzed
  Bool_t              fAnalyzeJets;           // trigger if jets should be processed
  Bool_t              fAnalyzePythia;         // trigger if pythia properties should be processed
  Bool_t              fHasTracks;             // trigger if tracks are actually valid
  Bool_t              fHasClusters;           // trigger if clusters are actually valid
  Bool_t              fHasJets;               // trigger if jets are actually valid
  Bool_t              fIsMC;                  // trigger if data is MC (for naming reasons)

  // ########## SOURCE INFORMATION
  TClonesArray*       fJetArray;              //! object containing the jets
  TClonesArray*       fTrackArray;            //! object containing the tracks
  TClonesArray*       fClusterArray;          //! object containing the clusters
  TString*            fJetArrayName;          // name of object containing the jets
  TString*            fTrackArrayName;        // name of object containing the tracks
  TString*            fClusterArrayName;      // name of object containing the tracks
  TString*            fRunNumbers;            // analyzed run numbers
  Int_t               fNumPtHardBins;         // Number of used pt hard bins

  // ########## JET/DIJET/RC PROPERTIES
  Double_t            fSignalJetRadius;       // Radius for the signal jets
  Int_t               fNumberExcludedJets;    // Number of jets to be excluded

  // ########## CUTS 
  Double_t            fSignalJetEtaWindow;    // +- window in eta for signal jets
  Double_t            fTrackEtaWindow;        // +- window in eta for tracks
  Double_t            fClusterEtaWindow;      // +- window in eta for clusters
  Double_t            fVertexWindow;          // +- window in Z for the vertex
  Double_t            fVertexMaxR;            // +- window in R for the vertex (distance in xy-plane)
  Double_t            fMinTrackPt;            // Min track pt to be accepted
  Double_t            fMinClusterPt;          // Min track pt to be accepted
  Double_t            fMinJetPt;              // Min jet pt to be accepted
  Double_t            fMinJetArea;            // Min jet area to be accepted

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

  AliAnalysisTaskQualityAssurancePA(const AliAnalysisTaskQualityAssurancePA&); // not implemented
  AliAnalysisTaskQualityAssurancePA& operator=(const AliAnalysisTaskQualityAssurancePA&); // not implemented

  ClassDef(AliAnalysisTaskQualityAssurancePA, 2); // QA helper task for pA

};
#endif
