#ifndef ALIANALYSISTASKCHARGEDJETSPA_H
#define ALIANALYSISTASKCHARGEDJETSPA_H

//#define DEBUGMODE

class TList;
class TClonesArray;
class TString;
class AliEmcalJet;
class AliRhoParameter;
class AliVParticle;
class AliLog;
class AliAnalysisUtils;
class TRandom3;
class AliESDtrack;

#include "AliESDtrackCuts.h"
#include "THn.h"

// Class for the hybrid track cuts
class AliESDHybridTrackcuts {
 public:
  AliESDHybridTrackcuts() : fMainCuts(0), fAdditionalCuts(0) {};
  AliESDHybridTrackcuts(const AliESDHybridTrackcuts& pd) : fMainCuts(0), fAdditionalCuts(0) {((AliESDHybridTrackcuts &) pd).Copy(*this);};
  virtual ~AliESDHybridTrackcuts() {delete fAdditionalCuts; delete fMainCuts;};

  virtual void Copy(AliESDHybridTrackcuts &c) const
  {
    AliESDHybridTrackcuts& target = (AliESDHybridTrackcuts &) c;
    target.fAdditionalCuts = fAdditionalCuts;
    target.fMainCuts = fMainCuts;
  };

  // Meta function to accept hybrid tracks
  Int_t AcceptTrack (const AliESDtrack* esdTrack) {
    if (AcceptNormalTrack(esdTrack))
      return 1;
    else if (AcceptAdditionalTrack(esdTrack))
      return 2;
    else
      return 0;
  };

  Bool_t  AcceptNormalTrack (const AliESDtrack* esdTrack) {if(fMainCuts) return fMainCuts->AcceptTrack(esdTrack); else return kFALSE;};
  Bool_t  AcceptAdditionalTrack (const AliESDtrack* esdTrack) {if(fAdditionalCuts) return fAdditionalCuts->AcceptTrack(esdTrack); else return kFALSE;};
  void    SetMainCuts (AliESDtrackCuts* trackCuts) {fMainCuts = trackCuts;};
  void    SetAdditionalCuts (AliESDtrackCuts* trackCuts) {fAdditionalCuts = trackCuts;};
  AliESDtrackCuts* GetAdditionalCuts () {return fAdditionalCuts;};
  AliESDtrackCuts* GetMainCuts () {return fMainCuts;};

  AliESDHybridTrackcuts &operator=(const AliESDHybridTrackcuts &c) {if (this != &c) ((AliESDHybridTrackcuts &) c).Copy(*this); return *this;};
 protected:
  AliESDtrackCuts*  fMainCuts;          // trackcut object for the global tracks
  AliESDtrackCuts*  fAdditionalCuts;    // trackcut object for the complementary tracks

// ClassDef(AliESDHybridTrackcuts, 2); 
};

class AliAnalysisTaskChargedJetsPA : public AliAnalysisTaskSE {
 public:
  enum{
    kMaxMatch=5
  };
  static const double kMaxChi2;
  // ######### CONTRUCTORS/DESTRUCTORS AND STD FUNCTIONS
  AliAnalysisTaskChargedJetsPA() : AliAnalysisTaskSE(), fOutputLists(), fCurrentOutputList(0), fDoJetAnalysis(1), fAnalyzeJetProfile(0), fAnalyzeTrackcuts(0), fAnalyzeJetConstituents(1), fParticleLevel(0), fUseDefaultVertexCut(1), fUsePileUpCut(1), fSetCentralityToOne(0), fNoExternalBackground(0), fBackgroundForJetProfile(0), fPartialAnalysisNParts(1), fPartialAnalysisIndex(0), fJetArray(0), fTrackArray(0), fBackgroundJetArray(0), fJetArrayName(), fTrackArrayName(), fBackgroundJetArrayName(), fRhoTaskName(), fRandConeRadius(0.4), fRandConeNumber(10), fSignalJetRadius(0.4), fBackgroundJetRadius(0.4), fNumberExcludedJets(-1), fMinEta(-0.9), fMaxEta(0.9), fMinJetEta(-0.5), fMaxJetEta(0.5), fMinTrackPt(0.150), fMinJetPt(5.0), fMinJetArea(0.5), fMinBackgroundJetPt(0.0), fMinNCrossedRows(70), fUsePtDepCrossedRowsCut(0), fNumberOfCentralityBins(20), fCentralityType("V0A"), fMatchTr(), fMatchChi(), fPrimaryVertex(0), fFirstLeadingJet(0), fSecondLeadingJet(0), fFirstLeadingKTJet(0), fSecondLeadingKTJet(0), fNumberSignalJets(0), fNumberSignalJetsAbove5GeV(0), fRandom(0), fHelperClass(0), fInitialized(0), fTaskInstanceCounter(0), fIsDEBUG(0), fIsPA(1), fNoTerminate(1), fEventCounter(0), fTempExcludedRCs(0), fTempAllRCs(1), fTempOverlapCounter(0), fTempMeanExclusionProbability(0), fHybridESDtrackCuts(0), fHybridESDtrackCuts_variedPtDep(0), fHybridESDtrackCuts_variedPtDep2(0)
  {
  // dummy
  }

  AliAnalysisTaskChargedJetsPA(const char *name, const char* trackArrayName, const char* jetArrayName, const char* backgroundJetArrayName, Bool_t analyzeJetProfile, Bool_t analyzeTrackcuts);
  virtual ~AliAnalysisTaskChargedJetsPA();
  virtual void     UserCreateOutputObjects();
  virtual void     UserExec(Option_t *option);
  virtual Bool_t   UserNotify();
  virtual void     Terminate(Option_t *);

  // ######### SETTERS/GETTERS
  void        SetDoJetAnalysis(Bool_t val) {fDoJetAnalysis = val;}
  void        SetAnalyzeJetProfile(Bool_t val) {fAnalyzeJetProfile = val;}
  void        SetAnalyzeTrackcuts(Bool_t val) {fAnalyzeTrackcuts = val;}
  void        SetAnalyzeJetConstituents(Bool_t val) {fAnalyzeJetConstituents = val;}
  void        SetAnalyzePartialEvents(Int_t nParts, Int_t index) {fPartialAnalysisNParts = nParts; fPartialAnalysisIndex = index;}
  void        SetUseDefaultVertexCut (Bool_t val) {fUseDefaultVertexCut = val;}
  void        SetUsePileUpCut (Bool_t val) {fUsePileUpCut = val;}
  void        SetNoTerminate (Bool_t val) {fNoTerminate = val;}
  void        SetIsPA (Bool_t val) {fIsPA = val;}
  void        SetCentralityToOne (Bool_t val) {fSetCentralityToOne = val;}
  void        SetNoExternalBackground (Bool_t val) {fNoExternalBackground = val;}
  void        SetBackgroundForJetProfile (Bool_t val) {fBackgroundForJetProfile = val;}
  void        SetMinNCrossedRows(Int_t val) {fMinNCrossedRows = val;}
  void        SetUsePtDepCrossedRowsCut(Bool_t val) {fUsePtDepCrossedRowsCut = val;}

  void        SetNumberOfCentralityBins(Int_t val) {fNumberOfCentralityBins = val;} 
  void        SetTrackMinPt(Double_t minPt) {fMinJetPt = minPt;}
  void        SetSignalJetMinPt(Double_t minPt) {fMinJetPt = minPt;}
  void        SetSignalJetMinArea(Double_t minArea) {fMinJetArea = minArea;}
  void        SetBackgroundJetMinPt(Double_t minPt) {fMinBackgroundJetPt = minPt;}
  void        SetRandConeRadius(Double_t radius) {fRandConeRadius = radius;}
  void        SetRandConeNumber(Int_t number) {fRandConeNumber = number;}
  void        SetSignalJetRadius(Double_t radius) {fSignalJetRadius = radius;}
  void        SetBackgroundJetRadius(Double_t radius) {fBackgroundJetRadius = radius;}
  void        SetCentralityType(const char* type) {fCentralityType = type;}
  void        SetExternalRhoTaskName(const char* name) {fRhoTaskName = name;}
  void        SetAcceptanceEta(Double_t minEta, Double_t maxEta) {fMinEta = minEta; fMaxEta = maxEta;}
  void        SetAcceptanceJetEta(Double_t minEta, Double_t maxEta) {fMinJetEta = minEta; fMaxJetEta = maxEta;}
  Int_t       GetInstanceCounter() {return fTaskInstanceCounter;}
  void        SetCurrentOutputList(Int_t i)
  {
    if(i==0)
    {
      fCurrentOutputList = fOutputLists[0];
    }
    else if(i==1)
    {
      if(fAnalyzeJetProfile)
        fCurrentOutputList = fOutputLists[1];
      else
        AliError("Non-existing output list demanded!");
    }
    else if(i==2)
    {
      if(!fAnalyzeJetProfile && fAnalyzeTrackcuts)
        fCurrentOutputList = fOutputLists[1];
      else if(fAnalyzeJetProfile && fAnalyzeTrackcuts)
        fCurrentOutputList = fOutputLists[2];
      else
        AliError("Non-existing output list demanded!");
    }
  }

 private:

  // ######### MAIN CALCULATION FUNCTIONS
  void        InitializeTrackcuts();
  void        GetLeadingJets();
  Double_t    GetCorrectedJetPt(AliEmcalJet* jet, Double_t background);
  Double_t    GetDeltaPt(Double_t rho, Double_t overlappingJetExclusionProbability = 0);

  void        GetKTBackgroundDensityAll(Int_t numberExcludeLeadingJets, Double_t& rhoPbPb, Double_t& rhoPbPbWithGhosts, Double_t& rhoCMS, Double_t& rhoImprovedCMS, Double_t& rhoMean, Double_t& rhoTrackLike);
  void        GetTRBackgroundDensity(Int_t numberExcludeLeadingJets, Double_t& rhoNoExclusion, Double_t& rhoConeExclusion02, Double_t& rhoConeExclusion04, Double_t& rhoConeExclusion06, Double_t& rhoConeExclusion08, Double_t& rhoExactExclusion);
  void        GetPPBackgroundDensity(Double_t& background);

  Double_t    GetConePt(Double_t eta, Double_t phi, Double_t radius);
  Double_t    GetCorrectedConePt(Double_t eta, Double_t phi, Double_t radius, Double_t background);
  Int_t       GetConeConstituentCount(Double_t eta, Double_t phi, Double_t radius);
  Double_t    GetExternalRho();
  void        CreateJetProfilePlots(Double_t bgrd);
  void        CreateCutHistograms();
  void        CreateITSTPCMatchingHistograms();
  void        GetPerpendicularCone(Double_t vecPhi, Double_t vecTheta, Double_t& conePt);

  // ######### CHECK FUNCTIONS
  Bool_t      IsTrackInAcceptance(AliVParticle* track);
  Bool_t      IsTrackInCone(AliVTrack* track, Double_t eta, Double_t phi, Double_t radius);
  Bool_t      IsTrackInJet(AliEmcalJet* jet, Int_t trackIndex);
  Bool_t      IsJetOverlapping(AliEmcalJet* jet1, AliEmcalJet* jet2);

  Bool_t      IsEventInAcceptance(AliVEvent* event);
  Bool_t      IsBackgroundJetInAcceptance(AliEmcalJet* jet);
  Bool_t      IsSignalJetInAcceptance(AliEmcalJet* jet, Bool_t usePtCut = kFALSE);
  
  // ######### HELPER FUNCTIONS
  Double_t    EtaToTheta(Double_t arg);
  Double_t    ThetaToEta(Double_t arg);
  Double_t    GetDeltaPhi(Double_t phi1, Double_t phi2);
  Double_t    MCGetOverlapCircleRectancle(Double_t cPosX, Double_t cPosY, Double_t cRadius, Double_t rPosXmin, Double_t rPosXmax, Double_t rPosYmin, Double_t rPosYmax);
  Double_t    MCGetOverlapMultipleCirclesRectancle(Int_t numCircles, std::vector<Double_t> cPosX, std::vector<Double_t> cPosY, Double_t cRadius, Double_t rPosXmin, Double_t rPosXmax, Double_t rPosYmin, Double_t rPosYmax);
  void        Match(AliESDtrack* tr0, AliESDtrack* tr1, Int_t& nmatch, Bool_t excludeMom = kFALSE, Double_t rotate=0);

  // ######### HISTOGRAM FUNCTIONS
  void        FillHistogram(const char * key, Double_t x);
  void        FillHistogram(const char * key, Double_t x, Double_t y);
  void        FillHistogram(const char * key, Double_t x, Double_t y, Double_t add);
  void        FillCutHistogram(const char * key, Double_t cut, Double_t pT, Double_t eta, Double_t phi, Int_t isAdditionalTrack);

  const char* GetHistoName(const char* name)
  {
    if (fParticleLevel)    
      return Form("%s_MC", name);
    return Form("%s", name);
  }
  template <class T> T* AddHistogram1D(const char* name = "CustomHistogram", const char* title = "NO_TITLE", const char* options = "", Int_t xBins = 100, Double_t xMin = 0.0, Double_t xMax = 20.0, const char* xTitle = "x axis", const char* yTitle = "y axis");
  template <class T> T* AddHistogram2D(const char* name = "CustomHistogram", const char* title = "NO_TITLE", const char* options = "", Int_t xBins = 100, Double_t xMin = 0.0, Double_t xMax = 20.0, Int_t yBins = 100, Double_t yMin = 0.0, Double_t yMax = 20.0, const char* xTitle = "x axis", const char* yTitle = "y axis", const char* zTitle = "z axis");
  THnF* AddCutHistogram(const char* name, const char* title, const char* cutName, Int_t nBins, Double_t xMin, Double_t xMax);
  void  BinLogAxis(const THn *h, Int_t axisNumber);

  // ######### STANDARD FUNCTIONS
  void      Calculate(AliVEvent* event);
  void      ExecOnce();
  void      Init ();

  std::vector<TList*> fOutputLists;           //! Output lists
  TList*              fCurrentOutputList;     //! Currently selected list where the histograms will be saved to
  // ########## USAGE TRIGGERS 
  Bool_t              fDoJetAnalysis;         // trigger if jets/tracks etc. should be analyzed
  Bool_t              fAnalyzeJetProfile;     // trigger if jet profile should be analyzed
  Bool_t              fAnalyzeTrackcuts;      // trigger if trackcuts should be analyzed
  Bool_t              fAnalyzeJetConstituents;// trigger if constituents should be analyzed
  Bool_t              fParticleLevel;         // trigger if data is kinematics only (for naming reasons)
  Bool_t              fUseDefaultVertexCut;   // trigger if automatic vertex cut from helper class should be done
  Bool_t              fUsePileUpCut;          // trigger if pileup cut should be done
  Bool_t              fSetCentralityToOne;    // trigger if centrality val. should be set to one for every event (failsafe)
  Bool_t              fNoExternalBackground;  // External background is set to 0 (e.g. for PP)
  Int_t               fBackgroundForJetProfile; // Which background will be subtracted for the profile calculation (0=External,1=Improved CMS,2=CMS,3=PP,4=TR,5=None)
  Int_t               fPartialAnalysisNParts; // take only every Nth event
  Int_t               fPartialAnalysisIndex;  // using e.g. only every 5th event, this specifies which one
  

  // ########## SOURCE INFORMATION
  TClonesArray*       fJetArray;              //! object containing the jets
  TClonesArray*       fTrackArray;            //! object containing the tracks
  TClonesArray*       fBackgroundJetArray;    //! object containing background jets
  TString             fJetArrayName;          // name of object containing the jets
  TString             fTrackArrayName;        // name of object containing the tracks
  TString             fBackgroundJetArrayName;// name of object containing event wise bckgrds
  TString             fRhoTaskName;           // name of rho task for this analysis
  // ########## JET/DIJET/RC PROPERTIES
  Double_t            fRandConeRadius;        // Radius for the random cones
  Int_t               fRandConeNumber;        // Number of random cones thrown per event
  Double_t            fSignalJetRadius;       // Radius for the signal jets
  Double_t            fBackgroundJetRadius;   // Radius for the KT background jets
  Int_t               fNumberExcludedJets;    // Number of jets to be excluded from backgrounds
  // ########## CUTS 
  Double_t            fMinEta;                // min eta of tracks
  Double_t            fMaxEta;                // max eta of tracks
  Double_t            fMinJetEta;             // min eta of jets
  Double_t            fMaxJetEta;             // max eta of jets
  Double_t            fMinTrackPt;            // Min track pt to be accepted
  Double_t            fMinJetPt;              // Min jet pt to be accepted
  Double_t            fMinJetArea;            // Min jet area to be accepted
  Double_t            fMinBackgroundJetPt;    // Min jet pt to be accepted as background jet
  Int_t               fMinNCrossedRows;       // Min number of crossed TPC rows for trackcut analysis
  Bool_t              fUsePtDepCrossedRowsCut;// Trigger if linear pT dep. for crossed rows cut should be applied
  Int_t               fNumberOfCentralityBins;// Number of centrality bins used for histograms
  TString             fCentralityType;        // Used centrality estimate (V0A, V0C, V0M, ...)

  AliESDtrack*        fMatchTr[kMaxMatch];    //! Helper variables track matching
  Double_t            fMatchChi[kMaxMatch];   //! Helper variables track matching


  // ########## EVENT PROPERTIES
  const AliVVertex*   fPrimaryVertex;         //! Vertex found per event
  AliEmcalJet*        fFirstLeadingJet;       //! leading jet in event
  AliEmcalJet*        fSecondLeadingJet;      //! next to leading jet in event
  AliEmcalJet*        fFirstLeadingKTJet;     //! leading kT jet in event
  AliEmcalJet*        fSecondLeadingKTJet;    //! next to leading kT jet in event
  Int_t               fNumberSignalJets;      // Number of signal jets in event
  Int_t               fNumberSignalJetsAbove5GeV; // Number of signal jets in event > 5GeV
  // ########## GENERAL VARS
  TRandom3*           fRandom;                //! A random number
  AliAnalysisUtils*   fHelperClass;           //! Vertex selection helper
  Bool_t              fInitialized;           // trigger if tracks/jets are loaded
  Int_t               fTaskInstanceCounter;   // for naming reasons
  Bool_t              fIsDEBUG;               // Debug trigger
  Bool_t              fIsPA;                  // pPb trigger
  Bool_t              fNoTerminate;           // don't use terminate routines
  ULong_t             fEventCounter;          // Internal event counter
  // ########### HELPER VARS
  Int_t               fTempExcludedRCs;       // used for delta pt signal exclusion
  Int_t               fTempAllRCs;            // used for delta pt signal exclusion
  Int_t               fTempOverlapCounter;    // used for delta pt signal exclusion
  Double_t            fTempMeanExclusionProbability;   // used for delta pt signal exclusion

  AliESDHybridTrackcuts* fHybridESDtrackCuts; //! these trackcuts are applied
  AliESDHybridTrackcuts* fHybridESDtrackCuts_variedPtDep; //! these trackcuts are applied
  AliESDHybridTrackcuts* fHybridESDtrackCuts_variedPtDep2; //! these trackcuts are applied

  AliAnalysisTaskChargedJetsPA(const AliAnalysisTaskChargedJetsPA&);
  AliAnalysisTaskChargedJetsPA& operator=(const AliAnalysisTaskChargedJetsPA&);

  ClassDef(AliAnalysisTaskChargedJetsPA, 3); // Charged jet analysis for pA

};
#endif
