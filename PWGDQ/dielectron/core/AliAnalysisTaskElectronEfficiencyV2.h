# ifndef AliAnalysisTaskElectronEfficiencyV2_H
# define AliAnalysisTaskElectronEfficiencyV2_H

//###########################################################
//#                                                         #
//#             Single Electron Efficiency Task             #
//#             Pair Efficiency Task                        #
//#             Track Resolution                            #
//#                                                         #
//#  Authors:                                               #
//#   Carsten Klein, Uni Frankfurt / Carsten.Klein@cern.ch  #
//#                                                         #
//###########################################################

#include "AliAnalysisTaskSE.h"
#include "AliDielectronSignalMC.h"
#include "TString.h"
#include "THnSparse.h"
class TH1F;
class TH2F;
class TH3D;


class AliTriggerAnalysis;
class AliAnalysisFilter;
class AliPIDResponse;

class TFile;
class TTree;
class TList;

class TObjArray;


class AliAnalysisTaskElectronEfficiencyV2 : public AliAnalysisTaskSE {
public:
  // two class constructors
  AliAnalysisTaskElectronEfficiencyV2() ;

  AliAnalysisTaskElectronEfficiencyV2(const char* name);
  // class destructor
   virtual ~AliAnalysisTaskElectronEfficiencyV2();
   // called once at beginning or runtime
   virtual void UserCreateOutputObjects();
   // called for each event
   virtual void UserExec(Option_t* option);
   // called at end of analysis
   virtual void Terminate(Option_t* option);

   enum Detector {kITS, kTPC, kTOF};
   Bool_t               GetEnablePhysicsSelection() const   {return fSelectPhysics; }
   Int_t                GetTriggerMask() const              {return fTriggerMask; }
   AliAnalysisCuts*     GetEventFilter()                    {return fEventFilter;}

   std::vector<double>  GetPtBins() const                   {return fPtBins;}
   std::vector<double>  GetEtaBins() const                  {return fEtaBins;}
   std::vector<double>  GetPhiBins() const                  {return fPhiBins;}

   // MC Signal setter
   void   AddSingleLegMCSignal(AliDielectronSignalMC signal1)         {fSingleLegMCSignal.push_back(signal1);}
   void   AddPairMCSignal(AliDielectronSignalMC signal1)              {fPairMCSignal.push_back(signal1);}
   void   AddMCSignalsWhereDielectronPairNotFromSameMother(std::vector<bool> vec) {fDielectronPairNotFromSameMother = vec;}

   // PID correction functions
   void   SetCentroidCorrFunction(Detector det, TObject *fun, UInt_t varx, UInt_t vary=0, UInt_t varz=0);
   void   SetWidthCorrFunction   (Detector det, TObject *fun, UInt_t varx, UInt_t vary=0, UInt_t varz=0);


   // Generator
   void   SetGeneratorName         (TString generatorName) { fGeneratorName = generatorName;}
   void   SetGeneratorMCSignalName (TString generatorName) { fGeneratorMCSignalName  = generatorName;}
   void   SetGeneratorULSSignalName(TString generatorName) { fGeneratorULSSignalName = generatorName;}

   void   SetCheckGenID(Bool_t flag) { fCheckGenID = flag;}
   void   SetGeneratorIndex         (std::vector<UInt_t> generatorIndex) { fGeneratorIndex = generatorIndex;}
   void   SetGeneratorMCSignalIndex (std::vector<UInt_t> generatorIndex) { fGeneratorMCSignalIndex  = generatorIndex;}
   void   SetGeneratorULSSignalIndex(std::vector<UInt_t> generatorIndex) { fGeneratorULSSignalIndex = generatorIndex;}


   // Event setter
   void   SetEnablePhysicsSelection(Bool_t selectPhysics)   {fSelectPhysics = selectPhysics;}
   void   SetTriggerMask(Int_t triggermask)                 {fTriggerMask = triggermask;}
   void   SetEventFilter(AliAnalysisCuts * const filter)    {fEventFilter = filter;}
   void   SetCentrality(double cent_min, double cent_max)   {fMinCentrality = cent_min; fMaxCentrality = cent_max;} // To ignore centrality use SetCentrality(-1, -1)
   void   SetCentralityEstimator(TString estimator)         {fCentralityEst = estimator;}

   void   SetCentralityFile(std::string filename) {fCentralityFilename = filename; }

  void   SetCentralityFileFromAlien(std::string filename) {fCentralityFilenameFromAlien = filename; }

   // Support Histos
   void   SetSupportHistoMCSignalAndCutsetting(int nMCSignal, int nCutsetting) {fSupportMCSignal = nMCSignal; fSupportCutsetting = nCutsetting;}

   // Resolution setter
   void   SetResolutionFile(std::string filename) {fResoFilename = filename; }
   void   SetResolutionFileFromAlien(std::string filename) {fResoFilenameFromAlien = filename; }
   void   SetSmearGenerated(bool setSmearingGen) { fDoGenSmearing = setSmearingGen; }
   void   SetResolutionDeltaPtBinsLinear (const double min, const double max, const unsigned int steps){SetBinsLinear("ptDelta_reso", min, max, steps);}
   void   SetResolutionRelPtBinsLinear   (const double min, const double max, const unsigned int steps){SetBinsLinear("ptRel_reso", min, max, steps);}
   void   SetResolutionEtaBinsLinear  (const double min, const double max, const unsigned int steps){SetBinsLinear("eta_reso", min, max, steps);}
   void   SetResolutionPhiBinsLinear  (const double min, const double max, const unsigned int steps){SetBinsLinear("phi_reso", min, max, steps);}
   void   SetResolutionThetaBinsLinear(const double min, const double max, const unsigned int steps){SetBinsLinear("theta_reso", min, max, steps);}

   // single electron binning setter
   void   SetPtBins(std::vector<double> ptBins)             {fPtBins = ptBins;}
   void   SetPtBinsLinear(const double min, const double max, const unsigned int steps) {SetBinsLinear("pt", min, max, steps);}
   void   SetEtaBins(std::vector<double> etaBins)           {fEtaBins = etaBins;}
   void   SetEtaBinsLinear(const double min, const double max, const unsigned int steps){SetBinsLinear("eta", min, max, steps);}
   void   SetPhiBins(std::vector<double> phiBins)           {fEtaBins = phiBins;}
   void   SetPhiBinsLinear(const double min, const double max, const unsigned int steps){SetBinsLinear("phi", min, max, steps);}
   void   SetThetaBinsLinear(const double min, const double max, const unsigned int steps){SetBinsLinear("theta", min, max, steps);}

   // pair binning setter
   void   SetMassBins(std::vector<double> massBins){fMassBins=massBins;}
   void   SetMassBinsLinear(const double min, const double max, const unsigned int steps){SetBinsLinear("mass", min, max, steps);}
   void   SetPairPtBins(std::vector<double> pairptBins){ fPairPtBins = pairptBins;}
   void   SetPairPtBinsLinear(const double min, const double max, const unsigned int steps){SetBinsLinear("pairpt", min, max, steps);}
   void   SetPhiVBins(std::vector<double> phivBins){fPhiVBins=phivBins;}
   void   SetPhiVBinsLinear(const double min, const double max, const unsigned int steps){SetBinsLinear("phiv", min, max, steps);}

   // Pair related setter
   void   SetDoPairing(Bool_t doPairing) {fDoPairing = doPairing;}
   void   SetULSandLS(Bool_t doULSandLS) {fDoULSandLS = doULSandLS;}
   void   SetDeactivateLS(Bool_t deactivateLS) {fDeactivateLS = deactivateLS;}
   void   SetKinematicCuts(double ptMin, double ptMax, double etaMin, double etaMax) {fPtMin = ptMin; fPtMax = ptMax; fEtaMin = etaMin; fEtaMax = etaMax;}
   void   SetFillPhiV(Bool_t doPhiV) {fDoFillPhiV = doPhiV;}
   void   SetPhiVCut(Bool_t apply, Double_t maxMee, Double_t minphiv){fApplyPhivCut = apply; fMaxMee = maxMee; fMinPhiV = minphiv;}

   // Single leg from Pair related setter
   void   SetWriteLegsFromPair(bool enable){fWriteLegsFromPair = enable;}
   void   SetPtMinLegsFromPair(const double ptMin){fPtMinLegsFromPair = ptMin;}
   void   SetPtMaxLegsFromPair(const double ptMax){fPtMaxLegsFromPair = ptMax;}
   void   SetEtaMinLegsFromPair(const double etaMin){fEtaMinLegsFromPair = etaMin;}
   void   SetEtaMaxLegsFromPair(const double etaMax){fEtaMaxLegsFromPair = etaMax;}
   void   SetPhiMinLegsFromPair(const double phiMin){fPhiMinLegsFromPair = phiMin;}
   void   SetPhiMaxLegsFromPair(const double phiMax){fPhiMaxLegsFromPair = phiMax;}
   void   SetOpAngleMinLegsFromPair(const double opAngleMin){fOpAngleMinLegsFromPair = opAngleMin;}
   void   SetOpAngleMaxLegsFromPair(const double opAngleMax){fOpAngleMaxLegsFromPair = opAngleMax;}
   void   SetPtNBinsLegsFromPair(const int ptNBins){fPtNBinsLegsFromPair = ptNBins;}
   void   SetEtaNBinsLegsFromPair(const int etaNBins){fEtaNBinsLegsFromPair = etaNBins;}
   void   SetPhiNBinsLegsFromPair(const int phiNBins){fPhiNBinsLegsFromPair = phiNBins;}
   void   SetOpAngleNBinsLegsFromPair(const int opAngleNBins){fOpAngleNBinsLegsFromPair = opAngleNBins;}

   // Set Cocktail waiting
   void SetDoCocktailWeighting(bool doCocktailWeight) { fDoCocktailWeighting = doCocktailWeight; }
   void SetCocktailWeighting(std::string CocktailFilename) { fCocktailFilename = CocktailFilename; }
   void SetCocktailWeightingFromAlien(std::string CocktailFilenameFromAlien) { fCocktailFilenameFromAlien = CocktailFilenameFromAlien; }

   // Generator related setter
   void   SetMinPtGen(double ptMin)   {fPtMinGen = ptMin;}; // Look only at particles which are above a threshold. (reduces computing time/less tracks when looking at secondaries)
   void   SetMaxPtGen(double ptMax)   {fPtMaxGen = ptMax;};
   void   SetMinEtaGen(double etaMin) {fEtaMinGen = etaMin;}; // Look only at particles which are above a threshold. (reduces computing time/less tracks when looking at secondaries)
   void   SetMaxEtaGen(double etaMax) {fEtaMaxGen = etaMax;};

   // Track cuts setter
   void   AddTrackCuts(AliAnalysisFilter* filter) {fTrackCuts.push_back(filter);}

  class Particle{
  public:
    Particle() :
      fPt(-99), fEta(-99), fPhi(-99), fCharge(-99), fPt_smeared(0.), fEta_smeared(0.), fPhi_smeared(0.), fTrackID(0), fMotherID(0), fMCSignalPair(false), fULSSignalPair(false), isMCSignal(), isReconstructed(), DielectronPairFromSameMother() {}
    Particle(double pt, double eta, double phi, short charge) :
      fPt(pt), fEta(eta), fPhi(phi), fCharge(charge), fPt_smeared(0.), fEta_smeared(0.), fPhi_smeared(0.), fTrackID(0), fMotherID(0), fMCSignalPair(false), fULSSignalPair(false), isMCSignal(), isReconstructed(), DielectronPairFromSameMother() {}

    void SetTrackID(int id) {fTrackID = id;}
    void SetMotherID(int id) {fMotherID = id;}
    void SetMCSignalPair (bool value) {fMCSignalPair = value;}
    void SetULSSignalPair(bool value) {fULSSignalPair = value;}
    void SetDielectronPairFromSameMother(std::vector<Bool_t> vec){DielectronPairFromSameMother = vec;}

    int  GetTrackID() {return fTrackID;}
    int  GetMotherID() {return fMotherID;}
    bool GetMCSignalPair() {return fMCSignalPair;}
    bool GetULSSignalPair() {return fULSSignalPair;}

    double  fPt;
    double  fEta;
    double  fPhi;
    short   fCharge;
    double  fPt_smeared;
    double  fEta_smeared;
    double  fPhi_smeared;
    int     fTrackID;
    int     fMotherID;
    bool    fMCSignalPair;
    bool    fULSSignalPair;
    std::vector<Bool_t> isMCSignal;
    std::vector<Bool_t> isReconstructed;
    std::vector<Bool_t> DielectronPairFromSameMother;
  };

private:
  enum {kAllEvents=0, kPhysicsSelectionEvents, kFilteredEvents , kCentralityEvents, kLastBin};

  void    SetBinsLinear(const std::string variable, const double min, const double max, const unsigned int steps);

  void    SetPIDResponse(AliPIDResponse *fPIDRespIn)        {fPIDResponse = fPIDRespIn;}
  void    CheckSingleLegMCsignals(std::vector<Bool_t>& vec, const int track);
  void    CheckPairMCsignals(std::vector<Bool_t>& vec, AliVParticle* part1, AliVParticle* part2);
  bool    CheckGenerator(int trackID, std::vector<unsigned int> vecHashes);
  bool    CheckGeneratorIndex(int trackID, std::vector<unsigned int> vecGenIDs);
  void    CheckIfFromMotherWithDielectronAsDaughter(Particle& part);
  Bool_t  CheckIfOneIsTrue(std::vector<Bool_t>& vec);

  Particle    CreateParticle(AliVParticle* part);

  void    CreateSupportHistos();

  void    FillTrackHistograms(AliVParticle* track, AliVParticle* mcTrack);

  TLorentzVector ApplyResolution(double pt, double eta, double phi, short ch);
  Double_t GetSmearing(TObjArray *arr, Double_t x);

  double GetWeight(Particle part1, Particle part2, double motherpt);
  double PhivPair(Double_t MagField, Int_t charge1, Int_t charge2, TVector3 dau1, TVector3 dau2);
  Double_t CalculateNbins();

  AliAnalysisCuts*  fEventFilter; // event filter

  TFile* fResoFile;
  std::string fResoFilename;
  std::string fResoFilenameFromAlien;
  TObjArray* fArrResoPt;
  TObjArray* fArrResoEta;
  TObjArray* fArrResoPhi_Pos;
  TObjArray* fArrResoPhi_Neg;

  TList* fOutputList;
  TList* fSingleElectronList;
  TList* fPairList;
  TList* fResolutionList;

  TH2D* fPGen_DeltaP;
  TH2D* fPGen_PrecOverPGen;
  TH2D* fPtGen_DeltaPt;
  TH2D* fPtGen_DeltaPtOverPtGen;
  TH2D* fPtGen_PtRecOverPtGen;
  TH2D* fPtGen_DeltaPt_wGenSmeared;
  TH2D* fPtGen_DeltaPtOverPtGen_wGenSmeared;
  TH2D* fPtGen_PtRecOverPtGen_wGenSmeared;
  TH2D* fPGen_DeltaEta;
  TH2D* fPtGen_DeltaEta;
  TH2D* fPGen_DeltaTheta;
  TH2D* fPGen_DeltaPhi_Ele;
  TH2D* fPGen_DeltaPhi_Pos;
  TH2D* fPtGen_DeltaPhi_Ele;
  TH2D* fPtGen_DeltaPhi_Pos;
  TH2D* fThetaGen_DeltaTheta;
  TH2D* fPhiGen_DeltaPhi;


  std::vector<double> fPtBins;
  std::vector<double> fEtaBins;
  std::vector<double> fPhiBins;
  std::vector<double> fThetaBins;
  std::vector<double> fResolutionDeltaPtBins;
  std::vector<double> fResolutionRelPtBins;
  std::vector<double> fResolutionEtaBins;
  std::vector<double> fResolutionPhiBins;
  std::vector<double> fResolutionThetaBins;
  std::vector<double> fMassBins;
  std::vector<double> fPairPtBins;
  std::vector<double> fPhiVBins;
  bool fDoGenSmearing;

  double  fPtMin; // Kinematic cut for pairing
  double  fPtMax; // Kinematic cut for pairing
  double  fEtaMin; // Kinematic cut for pairing
  double  fEtaMax; // Kinematic cut for pairing

  double  fPtMinGen;
  double  fPtMaxGen;
  double  fEtaMinGen;
  double  fEtaMaxGen;

  std::vector<AliDielectronSignalMC> fSingleLegMCSignal;
  std::vector<AliDielectronSignalMC> fPairMCSignal;
  std::vector<bool> fDielectronPairNotFromSameMother; // this is used to get electrons from charmed mesons in a environment where GEANT is doing the decay of D mesons, like in LHC18b5a

  TString fGeneratorName;
  TString fGeneratorMCSignalName;
  TString fGeneratorULSSignalName;
  std::vector<unsigned int> fGeneratorHashs;
  std::vector<unsigned int> fGeneratorMCSignalHashs;
  std::vector<unsigned int> fGeneratorULSSignalHashs;

  Bool_t fCheckGenID;
  std::vector<UInt_t> fGeneratorIndex;
  std::vector<UInt_t> fGeneratorMCSignalIndex;
  std::vector<UInt_t> fGeneratorULSSignalIndex;

  AliPIDResponse* fPIDResponse;
  AliVEvent*      fEvent;
  AliMCEvent*     fMC;
  AliVTrack*      fTrack;
  Bool_t          isAOD;

  Bool_t fSelectPhysics;
  Int_t  fTriggerMask;
  std::vector<AliAnalysisFilter*> fTrackCuts;
  TBits*  fUsedVars;                // used variables by AliDielectronVarManager

  int fSupportMCSignal; // Setting for which the support histograms are filled
  int fSupportCutsetting; // Setting for which the support histograms are filled

  TH1F* fHistEvents;
  TH1F* fHistEventStat;
  TH1F* fHistCentralityRaw;
  TH1F* fHistCentrality;
  TH1F* fHistVertex;
  TH1F* fHistVertexContibutors;
  TH1F* fHistNTracks;
  Double_t fMinCentrality;
  Double_t fMaxCentrality;
  TString fCentralityEst; // Which centrality estimator to use.

  TFile* fCentralityFile;
  std::string fCentralityFilename;
  std::string fCentralityFilenameFromAlien;
  TH1F* fHistCentralityCorrection;
  Double_t fNBinsCentralityCorr;
  Double_t fEntriesCentralityCorr;
  TList* fOutputListSupportHistos;

  std::vector<TH3D*> fHistGenPosPart;
  std::vector<TH3D*> fHistGenNegPart;
  std::vector<TH3D*> fHistGenSmearedPosPart;
  std::vector<TH3D*> fHistGenSmearedNegPart;
  std::vector<TH3D*> fHistRecPosPart;
  std::vector<TH3D*> fHistRecNegPart;

  std::vector<TH2D*> fHistGenPair;
  std::vector<TH2D*> fHistGenSmearedPair;
  std::vector<TObject*> fHistRecPair;
  std::vector<TH2D*> fHistGenPair_ULSandLS;
  std::vector<TH2D*> fHistGenSmearedPair_ULSandLS;
  std::vector<TH2D*> fHistRecPair_ULSandLS;

  bool fWriteLegsFromPair;
  double fPtMinLegsFromPair;
  double fPtMaxLegsFromPair;
  double fEtaMinLegsFromPair;
  double fEtaMaxLegsFromPair;
  double fPhiMinLegsFromPair;
  double fPhiMaxLegsFromPair;
  double fOpAngleMinLegsFromPair;
  double fOpAngleMaxLegsFromPair;
  int fPtNBinsLegsFromPair;
  int fEtaNBinsLegsFromPair;
  int fPhiNBinsLegsFromPair;
  int fOpAngleNBinsLegsFromPair;
  std::vector<THnSparseF*> fTHnSparseGenSmearedLegsFromPair;
  std::vector<THnSparseF*> fTHnSparseRecLegsFromPair;

  Bool_t fDoFillPhiV;
  Bool_t fApplyPhivCut;
  Double_t fMaxMee;
  Double_t fMinPhiV;

  Bool_t fDoPairing;
  Bool_t fDoULSandLS;
  Bool_t fDeactivateLS;
  std::vector<Particle> fGenNegPart;
  std::vector<Particle> fGenPosPart;
  std::vector<Particle> fRecNegPart;
  std::vector<Particle> fRecPosPart;

  bool fDoCocktailWeighting;
  std::string fCocktailFilename;
  std::string fCocktailFilenameFromAlien;
  TFile* fCocktailFile;
  TH1F* fPtPion;
  TH1F* fPtEta;
  TH1F* fPtEtaPrime;
  TH1F* fPtRho;
  TH1F* fPtOmega;
  TH1F* fPtPhi;
  TH1F* fPtJPsi;

  TH1* fPostPIDCntrdCorrTPC;     // post pid correction object for centroids in TPC
  TH1* fPostPIDWdthCorrTPC;      // post pid correction object for widths in TPC
  TH1* fPostPIDCntrdCorrITS;     // post pid correction object for centroids in ITS
  TH1* fPostPIDWdthCorrITS;      // post pid correction object for widths in ITS
  TH1* fPostPIDCntrdCorrTOF;     // post pid correction object for centroids in TOF
  TH1* fPostPIDWdthCorrTOF;      // post pid correction object for widths in TOF



  AliAnalysisTaskElectronEfficiencyV2(const AliAnalysisTaskElectronEfficiencyV2&); // not implemented
  AliAnalysisTaskElectronEfficiencyV2& operator=(const AliAnalysisTaskElectronEfficiencyV2&); // not implemented

  ClassDef(AliAnalysisTaskElectronEfficiencyV2, 7);
};


# endif

