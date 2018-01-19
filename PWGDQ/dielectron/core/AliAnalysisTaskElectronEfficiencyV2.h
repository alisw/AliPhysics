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

   Bool_t               GetEnablePhysicsSelection() const   {return fSelectPhysics; }
   Int_t                GetTriggerMask() const              {return fTriggerMask; }
   AliAnalysisCuts*     GetEventFilter()                    {return fEventFilter;}

   std::vector<double>  GetPtBins() const                   {return fPtBins;}
   std::vector<double>  GetEtaBins() const                  {return fEtaBins;}
   std::vector<double>  GetPhiBins() const                  {return fPhiBins;}

   // MC Signal setter
   void   AddSingleLegMCSignal(AliDielectronSignalMC signal1)         {fSingleLegMCSignal.push_back(signal1);}
   void   AddPairMCSignal(AliDielectronSignalMC signal1)              {fPairMCSignal.push_back(signal1);}

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
   void   SetMassBinsLinear(const double min, const double max, const unsigned int steps){SetBinsLinear("mass", min, max, steps);}
   void   SetPairPtBinsLinear(const double min, const double max, const unsigned int steps){SetBinsLinear("pairpt", min, max, steps);}

   // Pair related setter
   void   SetDoPairing(Bool_t doPairing) {fDoPairing = doPairing;}
   void   SetULSandLS(Bool_t doULSandLS) {fDoULSandLS = doULSandLS;}
   void   SetKinematicCuts(double ptMin, double ptMax, double etaMin, double etaMax) {fPtMin = ptMin; fPtMax = ptMax; fEtaMin = etaMin; fEtaMax = etaMax;}

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
      fPt(-99), fEta(-99), fPhi(-99), fCharge(-99), fPt_smeared(0.), fEta_smeared(0.), fPhi_smeared(0.), fTrackID(0), isMCSignal(), isReconstructed() {}
    Particle(double pt, double eta, double phi, short charge) :
      fPt(pt), fEta(eta), fPhi(phi), fCharge(charge), fPt_smeared(0.), fEta_smeared(0.), fPhi_smeared(0.), fTrackID(0), isMCSignal(), isReconstructed() {}
    void SetTrackID(int id) {fTrackID = id;}
    int  GetTrackID() {return fTrackID;}

    double  fPt;
    double  fEta;
    double  fPhi;
    short   fCharge;
    double  fPt_smeared;
    double  fEta_smeared;
    double  fPhi_smeared;
    int     fTrackID;
    std::vector<Bool_t> isMCSignal;
    std::vector<Bool_t> isReconstructed;
  };

private:
  enum {kAllEvents=0, kPhysicsSelectionEvents, kFilteredEvents , kCentralityEvents, kLastBin};

  void    SetBinsLinear(const std::string variable, const double min, const double max, const unsigned int steps);

  void    SetPIDResponse(AliPIDResponse *fPIDRespIn)        {fPIDResponse = fPIDRespIn;}
  void    CheckSingleLegMCsignals(std::vector<Bool_t>& vec, const int track);
  void    CheckPairMCsignals(std::vector<Bool_t>& vec, AliVParticle* part1, AliVParticle* part2);

  Bool_t  CheckIfOneIsTrue(std::vector<Bool_t>& vec);

  Particle    CreateParticle(AliVParticle* part);

  void    CreateSupportHistos();

  void    FillTrackHistograms(AliVParticle* track);

  TLorentzVector ApplyResolution(double pt, double eta, double phi, short ch);
  Double_t GetSmearing(TObjArray *arr, Double_t x);


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
  TH2D* fPtGen_DeltaPt;
  TH2D* fPtGen_DeltaPtOverPtGen;
  TH2D* fPGen_PrecOverPGen;
  TH2D* fPtGen_PtRecOverPtGen;
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
  std::vector<TH3D*> fHistGenSmearedPosPart;
  std::vector<TH3D*> fHistGenSmearedNegPart;
  std::vector<TH3D*> fHistRecPosPart;
  std::vector<TH3D*> fHistRecNegPart;

  std::vector<TH2D*> fHistGenPair;
  std::vector<TH2D*> fHistGenSmearedPair;
  std::vector<TH2D*> fHistRecPair;
  std::vector<TH2D*> fHistGenPair_ULSandLS;
  std::vector<TH2D*> fHistGenSmearedPair_ULSandLS;
  std::vector<TH2D*> fHistRecPair_ULSandLS;

  Bool_t fDoPairing;
  Bool_t fDoULSandLS;
  std::vector<Particle> fGenNegPart;
  std::vector<Particle> fGenPosPart;
  std::vector<Particle> fRecNegPart;
  std::vector<Particle> fRecPosPart;

  AliAnalysisTaskElectronEfficiencyV2(const AliAnalysisTaskElectronEfficiencyV2&); // not implemented
  AliAnalysisTaskElectronEfficiencyV2& operator=(const AliAnalysisTaskElectronEfficiencyV2&); // not implemented

  ClassDef(AliAnalysisTaskElectronEfficiencyV2, 1);
};


# endif
