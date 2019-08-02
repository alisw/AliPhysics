# ifndef AliAnalysisTaskEtaReconstruction_H
# define AliAnalysisTaskEtaReconstruction_H

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
class TH1F;
class TH2F;
class TH3D;

class TTree;

class AliTriggerAnalysis;
class AliAnalysisFilter;
class AliPIDResponse;

class TFile;
class TTree;
class TList;

class TObjArray;


class AliAnalysisTaskEtaReconstruction : public AliAnalysisTaskSE {
public:
  // two class constructors
  AliAnalysisTaskEtaReconstruction() ;

  AliAnalysisTaskEtaReconstruction(const char* name);
  // class destructor
   virtual ~AliAnalysisTaskEtaReconstruction();
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
   void   AddFourPairULSMCSignal(AliDielectronSignalMC signal1)          {fFourPairULSMCSignal.push_back(signal1);}
   void   AddFourPairLSMCSignal (AliDielectronSignalMC signal1)           {fFourPairLSMCSignal.push_back(signal1);}
   void   AddMCSignalsWhereDielectronPairNotFromSameMother(std::vector<bool> vec) {fDielectronPairNotFromSameMother = vec;}

   // PID correction functions
   // void   SetCentroidCorrFunction(Detector det, TObject *fun, UInt_t varx, UInt_t vary=0, UInt_t varz=0);
   // void   SetWidthCorrFunction   (Detector det, TObject *fun, UInt_t varx, UInt_t vary=0, UInt_t varz=0);


   // Generator
   void   SetGeneratorName         (TString generatorName) { fGeneratorName = generatorName;}
   void   SetGeneratorMCSignalName (TString generatorName) { fGeneratorMCSignalName  = generatorName;}

   // Event setter
   void   SetEnablePhysicsSelection(Bool_t selectPhysics)   {fSelectPhysics = selectPhysics;}
   void   SetTriggerMask(Int_t triggermask)                 {fTriggerMask = triggermask;}
   void   SetEventFilter(AliAnalysisCuts * const filter)    {fEventFilter = filter;}
   void   SetCentrality(double cent_min, double cent_max)             {fMinCentrality = cent_min; fMaxCentrality = cent_max;} //For pp and pPb analysis use SetCentrality(-1, -1)

   void   SetCentralityFile(std::string filename) {fCentralityFilename = filename; }

   // Support Histos
   void   SetSupportHistoMCSignalAndCutsetting(int nMCSignal, int nCutsetting) {fSupportMCSignal = nMCSignal; fSupportCutsetting = nCutsetting;}

   // Resolution setter
   void   SetResolutionFile(std::string filename) {fResoFilename = filename; }
   void   SetResolutionFileFromAlien(std::string filename) {fResoFilenameFromAlien = filename; }
   void   SetSmearGenerated(bool setSmearingGen) { fDoGenSmearing = setSmearingGen; }

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

   // Pair related setter
   void   SetDoPairing(Bool_t doPairing) {fDoPairing = doPairing;}
   void   SetDoFourPairing(Bool_t doFourPairing) {fDoFourPairing = doFourPairing;}
   void   SetKinematicCuts(double ptMin, double ptMax, double etaMin, double etaMax) {fPtMin = ptMin; fPtMax = ptMax; fEtaMin = etaMin; fEtaMax = etaMax;}

   // Single leg from Pair related setter
   // void   SetWriteTreeLegFromPair(bool enable){fWriteLegFromPair = enable;}

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
      fPt(-99), fEta(-99), fPhi(-99), fCharge(-99), fPt_smeared(0.), fEta_smeared(0.), fPhi_smeared(0.), fTrackID(0), fMotherID(0), fGrandMotherID(0), fMCSignalPair(false), isMCSignal(), isReconstructed(), DielectronPairFromSameMother() {}
    Particle(double pt, double eta, double phi, short charge) :
      fPt(pt), fEta(eta), fPhi(phi), fCharge(charge), fPt_smeared(0.), fEta_smeared(0.), fPhi_smeared(0.), fTrackID(0), fMotherID(0), fGrandMotherID(0), fMCSignalPair(false), isMCSignal(), isReconstructed(), DielectronPairFromSameMother() {}

    void SetTrackID(int id) {fTrackID = id;}
    void SetMotherID(int id) {fMotherID = id;}
    void SetGrandMotherID(int id) {fGrandMotherID = id;}
    void SetMCSignalPair (bool value) {fMCSignalPair = value;}
    void SetDielectronPairFromSameMother(std::vector<Bool_t> vec){DielectronPairFromSameMother = vec;}

    int  GetTrackID() {return fTrackID;}
    int  GetMotherID() {return fMotherID;}
    int  GetGrandMotherID() {return fGrandMotherID;}
    bool GetMCSignalPair() {return fMCSignalPair;}

    double  fPt;
    double  fEta;
    double  fPhi;
    short   fCharge;
    double  fPt_smeared;
    double  fEta_smeared;
    double  fPhi_smeared;
    int     fTrackID;
    int     fMotherID;
    int     fGrandMotherID;
    bool    fMCSignalPair;
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
  void    CheckIfFromMotherWithDielectronAsDaughter(Particle& part);
  Bool_t  CheckIfOneIsTrue(std::vector<Bool_t>& vec);

  Particle    CreateParticle(AliVParticle* part);

  void    CreateSupportHistos();

  void    FillTrackHistograms(AliVParticle* track, AliVParticle* mcTrack);

  TLorentzVector ApplyResolution(double pt, double eta, double phi, short ch);
  Double_t GetSmearing(TObjArray *arr, Double_t x);

  double GetWeight(Particle part1, Particle part2, double motherpt);

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
  TList* fFourPairList;

  std::vector<double> fPtBins;
  std::vector<double> fEtaBins;
  std::vector<double> fPhiBins;
  std::vector<double> fThetaBins;
  std::vector<double> fMassBins;
  std::vector<double> fPairPtBins;
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
  std::vector<AliDielectronSignalMC> fFourPairULSMCSignal;
  std::vector<AliDielectronSignalMC> fFourPairLSMCSignal;
  std::vector<bool> fDielectronPairNotFromSameMother; // this is used to get electrons from charmed mesons in a environment where GEANT is doing the decay of D mesons, like in LHC18b5a

  TString fGeneratorName;
  TString fGeneratorMCSignalName;
  std::vector<unsigned int> fGeneratorHashs;
  std::vector<unsigned int> fGeneratorMCSignalHashs;

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
  TH1F* fHistCentrality;
  TH1F* fHistVertex;
  TH1F* fHistVertexContibutors;
  TH1F* fHistNTracks;
  Double_t fMinCentrality;
  Double_t fMaxCentrality;

  TFile* fCentralityFile;
  std::string fCentralityFilename;
  TH1F* fHistCentralityCorrection;
  TList* fOutputListSupportHistos;

  std::vector<TH3D*> fHistGenPosPart;
  std::vector<TH3D*> fHistGenNegPart;
  std::vector<TH3D*> fHistGenNeuPart;
  std::vector<TH3D*> fHistGenSmearedPosPart;
  std::vector<TH3D*> fHistGenSmearedNegPart;
  std::vector<TH3D*> fHistRecPosPart;
  std::vector<TH3D*> fHistRecNegPart;

  std::vector<TH2D*> fHistGenPair;
  std::vector<TH2D*> fHistGenULSFourPair;
  std::vector<TH2D*> fHistGenLSFourPair;
  std::vector<TH2D*> fHistGenSmearedPair;
  std::vector<TH2D*> fHistRecPair;


  Bool_t fDoPairing;
  Bool_t fDoFourPairing;
  std::vector<Particle> fGenNegPart;
  std::vector<Particle> fGenPosPart;
  std::vector<Particle> fGenNeuPart;
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

  // TH1* fPostPIDCntrdCorrTPC;     // post pid correction object for centroids in TPC
  // TH1* fPostPIDWdthCorrTPC;      // post pid correction object for widths in TPC
  // TH1* fPostPIDCntrdCorrITS;     // post pid correction object for centroids in ITS
  // TH1* fPostPIDWdthCorrITS;      // post pid correction object for widths in ITS
  // TH1* fPostPIDCntrdCorrTOF;     // post pid correction object for centroids in TOF
  // TH1* fPostPIDWdthCorrTOF;      // post pid correction object for widths in TOF


  AliAnalysisTaskEtaReconstruction(const AliAnalysisTaskEtaReconstruction&); // not implemented
  AliAnalysisTaskEtaReconstruction& operator=(const AliAnalysisTaskEtaReconstruction&); // not implemented

  ClassDef(AliAnalysisTaskEtaReconstruction, 1);
};


# endif
