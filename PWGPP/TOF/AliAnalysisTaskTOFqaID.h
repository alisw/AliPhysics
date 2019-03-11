#ifndef ALIANALYSISTASKTOFQAID_h
#define ALIANALYSISTASKTOFQAID_h

class TString;
class THashList;
class AliAnalysisFilter;
class AliCDBManager;
class AliTOFcalib;
class AliTOFT0maker;
class AliTOFT0v1;
class AliTOFHeader;
class AliTOFChannelOnlineStatusArray;
class AliTOFcalib;
class AliESDpid;

#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskTOFqaID : public AliAnalysisTaskSE {
  public:
  enum ETrackCutSetTOFqa_t { kRun1Cuts = 0,
    kStd2010,
    kStd2010crossedRows,
    kStd2011,
    kStd2011crossedRows,
    kStd2015crossedRows,
    kNCutSetTOFqa };

  AliAnalysisTaskTOFqaID();
  AliAnalysisTaskTOFqaID(const char* name);
  AliAnalysisTaskTOFqaID(const AliAnalysisTaskTOFqaID& copy);
  AliAnalysisTaskTOFqaID& operator=(const AliAnalysisTaskTOFqaID& copy);
  virtual ~AliAnalysisTaskTOFqaID();

  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t* option);
  virtual void Terminate(Option_t*);

  Int_t GetStripIndex(const Int_t* in);
  void SetTrackFilter(AliAnalysisFilter* filter) { fTrackFilter = filter; };
  void SetMinPtCut(Float_t minpt) { fMatchingMomCut = minpt; }
  void SetMaxEtaCut(Float_t maxeta) { fMatchingEtaCut = maxeta; }
  void EnableAdvancedCheck(Bool_t enable) { fEnableAdvancedCheck = enable; };
  void EnableChargeSplit(Bool_t enable) { fEnableChargeSplit = enable; };
  void SetExpTimeHistoRange(Float_t min, Float_t max)
  {
    fExpTimeRangeMin = min;
    fExpTimeRangeMax = max;
  };
  void SetExpTimeHistoSmallRange(Float_t min, Float_t max)
  {
    fExpTimeSmallRangeMin = min;
    fExpTimeSmallRangeMax = max;
  };
  void SetExpTimeBinWidth(Float_t width) { fExpTimeBinWidth = width; };
  void SetResolutionMinP(Double_t min) { fResolutionMinP = min; };
  void SetResolutionMaxP(Double_t max) { fResolutionMaxP = max; };
  Bool_t SetSelectMCspecies(Bool_t enableMC, Int_t absPdgCode)
  {
    fIsMC = enableMC;
    fSelectedPdg = absPdgCode;
    return kTRUE;
  };
  TString GetSpeciesName(Int_t absPdgCode);
  void HistogramMakeUp(TH1* hist, Color_t color = -1, Int_t markerStyle = -1);
  Double_t GetPhiAtTPCouterRadius(AliESDtrack* track);
  void SetOCDBInfo(const char* cdbLocation, UInt_t runN)
  {
    fOCDBLocation = cdbLocation;
    fRunNumber = runN;
  }
  void SetVerboseMode(Bool_t verbose = kTRUE) { fVerbose = verbose; };                   // Setter for verbose flag
  void SetUseTOFT0CalibMode(Bool_t usecalib = kTRUE) { fUseTOFT0CalibMode = usecalib; }; // Setter for TOFT0V1 calib. mode

  // Binning arrays
  TArrayD fVariableBinsPt;   // array of bins for pt and p
  TArrayD fVariableBinsMult; // array of bins for multiplicity (e.g. TOF hits)

  protected:
  void AddTofBaseHisto(THashList* list, Int_t charge, TString suffix);
  void AddMatchingEffHisto(THashList* list, Int_t charge, TString suffix);
  void AddPidHisto(THashList* list, Int_t charge, TString suffix);
  void AddStartTimeHisto(THashList* list, TString suffix);
  void AddTrdHisto();
  void AddTofTrgHisto(TString suffix);

  void FillTofBaseHisto(AliESDtrack* track, Int_t charge, TString suffix);
  void FillMatchedTrkHisto(Int_t charge, TString suffix);
  void FillPrimaryTrkHisto(Int_t charge, TString suffix);
  void FillPidHisto(AliESDtrack* track, Int_t charge, TString suffix);
  void FillStartTimeHisto(TString suffix);
  void FillTrdHisto(AliESDtrack* track, Int_t charge);
  void FillStartTimeMaskHisto(TString suffix);
  void FillTofTrgHisto(TString suffix);

  Bool_t ComputeTimeZeroByTOF1GeV();
  Bool_t SelectMCspecies(AliMCEvent* ev, AliESDtrack* track);
  Bool_t ComputeMatchingEfficiency(THashList* list, TString variable);
  Bool_t IsTPCTOFMatched(AliESDtrack* track, Bool_t checkMatchLabel);
  Bool_t IsInTRD(AliESDtrack* track);
  Bool_t IsEventSelected(AliESDEvent* event);
  void LoadChannelMapsFromOCDB();
  Bool_t IsChannelGood(Int_t channel = -1);

  private:
  Int_t fRunNumber;                //run number
  AliESDEvent* fESD;               //ESD object
  AliMCEvent* fMCevent;            //MC event
  AliAnalysisFilter* fTrackFilter; //track filter object
  AliESDVertex* fVertex;           //pointer to the vertex object
  AliESDpid* fESDpid;              //pointer to the PID object
  //AliTOFT0v1 *        fTOFT0v1; // TOF-T0 v1
  AliTOFHeader* fTOFHeader; // TOF header needed for trigger info
  Int_t fNTOFtracks[3];     //number of tracks matching with TOF
  //Int_t fNPrimaryTracks; //number of primary tracks
  Float_t fT0[3];              //event time
  Float_t fSigmaSpecie[5];     //number of TOF PID sigmas, ie.fSigmaPion, fSigmaKaon, fSigmaProton;
  Double_t fTrkExpTimes[5];    //expected times from tracking for 5 mass hypothesis
  Double_t fThExpTimes[5];     //theoretical expected times for 5 mass hypothesis
  Bool_t fEnableAdvancedCheck; //flag to enable advanced checks
  Bool_t fEnableChargeSplit;   //flag to enable split for sign of charge

  Float_t fExpTimeBinWidth;                             //bin width for t-texp histos
  Float_t fExpTimeRangeMin, fExpTimeRangeMax;           //range of t-texp histogram
  Float_t fExpTimeSmallRangeMin, fExpTimeSmallRangeMax; //reduced range of t-texp histogram
  Int_t fnExpTimeBins;
  Int_t fnExpTimeSmallBins;

  TArrayD fMyTimeZeroTOF;                        //timeZero by TOF recomputed one from not calib., the other two for calib. mode
  TArrayD fMyTimeZeroTOFsigma;                   //timeZero sigma by TOF recomputed
  TArrayD fMyTimeZeroTOFtracks;                  // number of tracks used to recompute TOF_T0
  Bool_t fMyTimeZeroTOFstatus;                   // Status of the computed TOF_T0 (kTRUE -> OK, kFALSE -> not OK)
  Bool_t fIsMC;                                  //flag for MC
  Bool_t fVerbose;                               //Flag for verbose mode
  Bool_t fUseTOFT0CalibMode;                     // Flag to use AliTOFT0v1 in calibration mode and user mode
  Int_t fSelectedPdg;                            //pdg code of the selected specie (for MC only)
  Double_t fP;                                   //momentum
  Double_t fPt;                                  //transverse momentum
  Double_t fEta;                                 //psedorapidity
  Double_t fPhi;                                 //phi at vertex
  Double_t fTPCOuterPhi;                         //phi at outer tpc radius
  Double_t fL;                                   //integrated track lenght
  Double_t fMatchingMomCut;                      //minimum pT cut for matching eff vs eta, phi
  Double_t fMatchingEtaCut;                      //simmetric eta cut for matching eff vs pt, eta, phi
  Double_t fTof;                                 // Signal measured by TOF in ps
  Short_t fMCTOFMatch;                           //status of matching in MC
  TString fOCDBLocation;                         // ocdb path
  AliTOFChannelOnlineStatusArray* fChannelArray; //array of channel status
  AliTOFcalib* fCalib;                           //TOF calibration object
  Double_t fResolutionMinP;                      //Minimum momentum to compute the TOF resolution
  Double_t fResolutionMaxP;                      //Maximum momentum to compute the TOF resolution
  //output objects
  THashList* fHlist;         //list of general histos
  THashList* fHlistTimeZero; //list of timeZero related histos
  THashList* fHlistPID;      //list of PID-related histos
  THashList* fHlistTRD;      //list of general histos for positive tracks
  THashList* fHlistTrigger;  //list of general histos for TOF trg infos

  static const Int_t fnBinsEta = 200; // binning for eta
  static const Int_t fnBinsPhi = 72;  // binning for phi and phi_TPCouter
  static const Int_t fnBinsT0 = 140;  // binning for T0
  static const Double_t fBinsEta[2];  // binning for eta - max and min
  static const Double_t fBinsPhi[2];  // binning for phi and phi_TPCouter - max and min
  static const Double_t fBinsT0[2];   // binning for T0

  void SetVariableBinning(); // sets the array with variable binning

  ClassDef(AliAnalysisTaskTOFqaID, 9); // Analysis for the TOF QA
};

#endif
