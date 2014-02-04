#ifndef ALIANALYSISTASKTOFQAID_h
#define ALIANALYSISTASKTOFQAID_h

class TString;
class TList;
class AliAnalysisFilter;
class AliCDBManager;
class AliTOFcalib;
class AliTOFT0maker;
class AliTOFT0v1;
class AliTOFHeader;

#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskTOFqaID : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskTOFqaID();
  AliAnalysisTaskTOFqaID(const char *name);
  AliAnalysisTaskTOFqaID(const AliAnalysisTaskTOFqaID& copy);
  AliAnalysisTaskTOFqaID& operator= (const AliAnalysisTaskTOFqaID& copy);
  virtual ~AliAnalysisTaskTOFqaID();
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);
    
  Int_t   GetStripIndex(const Int_t * in);
  void    SetTrackFilter(AliAnalysisFilter *filter) {fTrackFilter = filter;};
  void    EnableAdvancedCheck(Bool_t enable){fEnableAdvancedCheck=enable;};
  void    SetExpTimeHistoRange(Float_t min, Float_t max){fExpTimeRangeMin=min; fExpTimeRangeMax=max;return;};
  void    SetExpTimeHistoSmallRange(Float_t min, Float_t max){fExpTimeSmallRangeMin=min; fExpTimeSmallRangeMax=max;return;};
  void    SetExpTimeBinWidth(Float_t width){fExpTimeBinWidth=width;return;};
  Bool_t  SetSelectMCspecies(Bool_t enableMC, Int_t absPdgCode){fIsMC=enableMC; fSelectedPdg=absPdgCode; return kTRUE;};
  TString GetSpeciesName(Int_t absPdgCode);
  void    HistogramMakeUp(TH1* hist, Color_t color, Int_t markerStyle,  TString drawOpt, TString newName, TString newTitle, TString xTitle, TString yTitle);
  Double_t GetPhiAtTPCouterRadius(AliESDtrack * track);

 protected:
  void    AddTofBaseHisto(TList *list, Int_t charge, TString suffix);
  void    AddMatchingEffHisto(TList *list, Int_t charge, TString suffix);
  void    AddPidHisto(TList *list, Int_t charge, TString suffix);
  void    AddStartTimeHisto(TList *list, TString suffix);
  void    AddTrdHisto();
  void    AddTofTrgHisto(TString suffix);

  void    FillTofBaseHisto(AliESDtrack * track, Int_t charge, TString suffix);
  void    FillMatchedTrkHisto(Int_t charge, TString suffix);
  void    FillPrimaryTrkHisto(Int_t charge, TString suffix);
  void    FillPidHisto(AliESDtrack * track,Int_t charge, TString suffix);
  void    FillStartTimeHisto(TString suffix);
  void    FillTrdHisto(AliESDtrack * track, Int_t charge);
  void    FillStartTimeMaskHisto(TString suffix);
  void    FillTofTrgHisto(TString suffix);

  Bool_t  ComputeTimeZeroByTOF1GeV();
  Bool_t  SelectMCspecies(AliMCEvent * ev, AliESDtrack * track);
  Bool_t  ComputeMatchingEfficiency(TList* list, TString variable);
  Bool_t  IsTPCTOFMatched(AliESDtrack * track);
  Bool_t  IsInTRD(AliESDtrack * track);
  Bool_t  IsEventSelected(AliESDEvent * event);

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
  Double_t            fMatchingMomCut;//pT cut for eta, phi matching eff 
  Double_t            fTof;
  //output objects
  TList *             fHlist;  //list of general histos
  TList *             fHlistTimeZero; //list of timeZero related histos
  TList *             fHlistPID; //list of PID-related histos
  TList *             fHlistTRD;  //list of general histos for positive tracks
  TList *             fHlistTrigger;  //list of general histos for TOF trg infos

  ClassDef(AliAnalysisTaskTOFqaID, 1); // example of analysis
};

#endif



