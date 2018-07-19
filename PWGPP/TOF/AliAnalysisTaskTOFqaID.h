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
  AliAnalysisTaskTOFqaID(const char *name);
  AliAnalysisTaskTOFqaID(const AliAnalysisTaskTOFqaID& copy);
  AliAnalysisTaskTOFqaID& operator= (const AliAnalysisTaskTOFqaID& copy);
  virtual ~AliAnalysisTaskTOFqaID();

  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);

  Int_t   GetStripIndex(const Int_t * in);
  void    SetTrackFilter(AliAnalysisFilter *filter) { fTrackFilter = filter; return; };
  void    SetMinPtCut(Float_t minpt) { fMatchingMomCut = minpt; return; }
  void    SetMaxEtaCut(Float_t maxeta) { fMatchingEtaCut = maxeta; return; }
  void    EnableAdvancedCheck(Bool_t enable) { fEnableAdvancedCheck = enable; return; };
  void    EnableChargeSplit(Bool_t enable) { fEnableChargeSplit = enable; return; };
  void    SetExpTimeHistoRange(Float_t min, Float_t max) { fExpTimeRangeMin = min; fExpTimeRangeMax = max; return;};
  void    SetExpTimeHistoSmallRange(Float_t min, Float_t max) { fExpTimeSmallRangeMin = min; fExpTimeSmallRangeMax = max; return;};
  void    SetExpTimeBinWidth(Float_t width) { fExpTimeBinWidth = width; return;};
  Bool_t  SetSelectMCspecies(Bool_t enableMC, Int_t absPdgCode) {fIsMC = enableMC; fSelectedPdg = absPdgCode; return kTRUE;};
  TString GetSpeciesName(Int_t absPdgCode);
  void    HistogramMakeUp(TH1* hist, Color_t color =-1, Int_t markerStyle = -1,  TString drawOpt = "", TString newName = "", TString newTitle = "", TString xTitle = "", TString yTitle = "");
  Double_t GetPhiAtTPCouterRadius(AliESDtrack * track);
  void    SetOCDBInfo(const char *cdbLocation, UInt_t runN) {fOCDBLocation=cdbLocation; fRunNumber=runN;}

 protected:
  void    AddTofBaseHisto(THashList *list, Int_t charge, TString suffix);
  void    AddMatchingEffHisto(THashList *list, Int_t charge, TString suffix);
  void    AddPidHisto(THashList *list, Int_t charge, TString suffix);
  void    AddStartTimeHisto(THashList *list, TString suffix);
  void    AddTrdHisto();
  void    AddTofTrgHisto(TString suffix);

  void    FillTofBaseHisto(AliESDtrack * track, Int_t charge, TString suffix);
  void    FillMatchedTrkHisto(Int_t charge, TString suffix);
  void    FillPrimaryTrkHisto(Int_t charge, TString suffix);
  void    FillPidHisto(AliESDtrack * track, Int_t charge, TString suffix);
  void    FillStartTimeHisto(TString suffix);
  void    FillTrdHisto(AliESDtrack * track, Int_t charge);
  void    FillStartTimeMaskHisto(TString suffix);
  void    FillTofTrgHisto(TString suffix);

  Bool_t  ComputeTimeZeroByTOF1GeV();
  Bool_t  SelectMCspecies(AliMCEvent * ev, AliESDtrack * track);
  Bool_t  ComputeMatchingEfficiency(THashList* list, TString variable);
  Bool_t  IsTPCTOFMatched(AliESDtrack * track, Bool_t checkMatchLabel);
  Bool_t  IsInTRD(AliESDtrack * track);
  Bool_t  IsEventSelected(AliESDEvent * event);
  void    LoadChannelMapsFromOCDB();
  Bool_t  IsChannelGood(Int_t channel = -1);

 private:
  Int_t               fRunNumber; //run number
  AliESDEvent *       fESD;    //ESD object
  AliMCEvent *        fMCevent;    //MC event
  AliAnalysisFilter * fTrackFilter; //track filter object
  AliESDVertex *      fVertex; //pointer to the vertex object
  AliESDpid *         fESDpid; //pointer to the PID object
  //AliTOFT0v1 *        fTOFT0v1; // TOF-T0 v1
  AliTOFHeader *      fTOFHeader; // TOF header needed for trigger info
  Int_t               fNTOFtracks[3]; //number of tracks matching with TOF
  //Int_t fNPrimaryTracks; //number of primary tracks
  Float_t             fT0[3]; //event time
  Float_t             fSigmaSpecie[5]; //number of TOF PID sigmas, ie.fSigmaPion, fSigmaKaon, fSigmaProton;
  Double_t            fTrkExpTimes[5]; //expected times from tracking for 5 mass hypothesis
  Double_t            fThExpTimes[5]; //theoretical expected times for 5 mass hypothesis
  Bool_t              fEnableAdvancedCheck; //flag to enable advanced checks
  Bool_t              fEnableChargeSplit; //flag to enable split for sign of charge

  Float_t             fExpTimeBinWidth;//bin width for t-texp histos
  Float_t             fExpTimeRangeMin, fExpTimeRangeMax; //range of t-texp histogram
  Float_t             fExpTimeSmallRangeMin, fExpTimeSmallRangeMax; //reduced range of t-texp histogram
  Int_t               fnExpTimeBins;
  Int_t               fnExpTimeSmallBins;

  Double_t            fMyTimeZeroTOF, fMyTimeZeroTOFsigma; //timeZero by TOF recomputed
  Int_t               fMyTimeZeroTOFtracks; // number of tracks used to recompute TOF_T0
  Bool_t              fIsMC; //flag for MC
  Int_t               fSelectedPdg; //pdg code of the selected specie (for MC only)
  Double_t            fP; //momentum
  Double_t            fPt; //transverse momentum
  Double_t            fEta; //psedorapidity
  Double_t            fPhi; //phi at vertex
  Double_t            fTPCOuterPhi; //phi at outer tpc radius
  Double_t            fL; //integrated track lenght
  Double_t            fMatchingMomCut;//minimum pT cut for matching eff vs eta, phi
  Double_t            fMatchingEtaCut;//simmetric eta cut for matching eff vs pt, eta, phi
  Double_t            fTof;
  Short_t             fMCTOFMatch; //status of matching in MC
  TString             fOCDBLocation;       // ocdb path
  AliTOFChannelOnlineStatusArray * fChannelArray; //array of channel status
  AliTOFcalib *       fCalib; //TOF calibration object
  //output objects
  THashList *             fHlist;  //list of general histos
  THashList *             fHlistTimeZero; //list of timeZero related histos
  THashList *             fHlistPID; //list of PID-related histos
  THashList *             fHlistTRD;  //list of general histos for positive tracks
  THashList *             fHlistTrigger;  //list of general histos for TOF trg infos

  static const Int_t fnBinsPt = 300; // binning for pt and p 
  static const Int_t fnBinsEta = 200; // binning for eta
  static const Int_t fnBinsPhi = 72; // binning for phi and phi_TPCouter
  static const Double_t fBinsPt[2]; // binning for pt and p - max and min
  static const Double_t fBinsEta[2]; // binning for eta - max and min
  static const Double_t fBinsPhi[2]; // binning for phi and phi_TPCouter - max and min
  Double_t fVariableBinsPt[fnBinsPt + 1]; // array of bins for pt and p

  void SetPtVariableBinning(); // sets the array with variable binning in p and pT
  
  ClassDef(AliAnalysisTaskTOFqaID, 6); // example of analysis
};

#endif
