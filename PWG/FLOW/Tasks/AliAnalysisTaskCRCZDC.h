/* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/**********************************
 * analysis task for CRC with ZDC *
 *                                *
 * author: Jacopo Margutti        *
 *         (margutti@nikhef.nl)   *
 **********************************/

#ifndef AliAnalysisTaskCRCZDC_H
#define AliAnalysisTaskCRCZDC_H

#include "AliAnalysisTaskSE.h"
#include "AliFlowTrackSimple.h"
#include "TGrid.h"

class AliCFManager;
class AliFlowEventCuts;
class AliFlowTrackCuts;
class AliFlowEventSimpleMaker;
class AliFlowEvent;
class AliFlowVector;
class TList;
class TF1;
class TRandom3;
class AliAnalysisTaskSE;
class TString;
class AliESDpid;
class AliGenEventHeader;
class AliGenPythiaEventHeader;
class AliGenHijingEventHeader;
class AliFlowTrack;
class AliAnalysisUtils;
class AliMultSelection;
class AliCentrality;
class AliStack;
class TROOT;
class TSystem;
class TFile;
class TH1F;
class TH2F;
class TProfile;
class TProfile2D;
class TProfile3D;
class TH3D;
class TH3F;

class AliAnalysisTaskCRCZDC : public AliAnalysisTaskSE {

public:

  enum kAnalysisInput{kESD=1, kAOD=2};
  AliAnalysisTaskCRCZDC();
  AliAnalysisTaskCRCZDC(const char *name, TString RPtype = "", Bool_t QAon = kFALSE, UInt_t seed=666, Bool_t bCandidates=kFALSE, Int_t StepZDCRecenter=-1);
  virtual ~AliAnalysisTaskCRCZDC();

  enum DataSet { k2010,
    k2011,
    k2015,
    k2015v6,
    k2015pidfix,
    kAny,
    k2018q,
    k2018r //@Shi
  };

  enum CentrEstimator {
    kV0M,
    kCL0,
    kCL1,
    kTRK
  };

  enum AnalysisType {
    kMCkine,
    kMCAOD,
    kAUTOMATIC,
    kMCESD,
    kTrackQA,
    kTracklets
  };

  virtual void InitializeRunArrays();

  // Implementation of interface methods
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *option);
  virtual void NotifyRun();

  void    SetAnalysisType(AnalysisType type) { this->fAnalysisType = type; }

  void    SetRPType(TString rptype) { this->fRPType = rptype; }
  TString GetRPType() const         { return this->fRPType; }

  void    SetMinMult(Int_t multmin)    {this->fMinMult = multmin; }
  Int_t   GetMinMult() const           {return this->fMinMult; }
  void    SetMaxMult(Int_t multmax)    {this->fMaxMult = multmax; }
  Int_t   GetMaxMult() const           {return this->fMaxMult; }

  void SetSubeventEtaRange(Double_t minA, Double_t maxA, Double_t minB, Double_t maxB)
  {this->fMinA = minA; this->fMaxA = maxA; this->fMinB = minB; this->fMaxB = maxB; }
  Double_t GetMinA() const {return this->fMinA;}
  Double_t GetMaxA() const {return this->fMaxA;}
  Double_t GetMinB() const {return this->fMinB;}
  Double_t GetMaxB() const {return this->fMaxB;}

  void DefineDeadZone( Double_t etaMin, Double_t etaMax, Double_t phiMin, Double_t phiMax )
  {this->fExcludedEtaMin = etaMin; this->fExcludedEtaMax = etaMax;
    this->fExcludedPhiMin = phiMin; this->fExcludedPhiMax = phiMax; }

  void          SetCutsEvent(AliFlowEventCuts* cutsEvent) {fCutsEvent=cutsEvent;}
  AliFlowEventCuts* GetCutsEvent() const {return fCutsEvent;}
  void          SetCutsRP(AliFlowTrackCuts* cutsRP);
  AliFlowTrackCuts* GetCutsRP() const {return fCutsRP;} //to be reimplemented
  void          SetCutsPOI(AliFlowTrackCuts* cutsPOI);
  AliFlowTrackCuts* GetCutsPOI() const {return fCutsPOI;} //to be reimplemented

  void          SetCFManager1(AliCFManager* cfmgr) {this->fCFManager1 = cfmgr; }
  AliCFManager* GetCFManager1() const {return this->fCFManager1; }
  void          SetCFManager2(AliCFManager* cfmgr) {this->fCFManager2 = cfmgr; }
  AliCFManager* GetCFManager2() const       {return this->fCFManager2; }
  TList*        GetQAList()      const      {return fQAList; }
  void          SetQAOn(Bool_t kt)        {fQAon = kt; }
  Bool_t        GetQAOn()   const         {return fQAon; }
  Bool_t        SelectPileup(AliAODEvent* aod);
  void          SetupPileUpRemovalFunctions18qPass3();
  void          SetupPileUpRemovalFunctions18rPass3();


  void          SetShuffleTracks(Bool_t b)  {fShuffleTracks=b;}

  // setters for common constants
  void SetNbinsMult( Int_t i ) { fNbinsMult = i; }
  void SetNbinsPt( Int_t i )   { fNbinsPt = i; }
  void SetNbinsPhi( Int_t i )  { fNbinsPhi = i; }
  void SetNbinsEta( Int_t i )  { fNbinsEta = i; }
  void SetNbinsQ( Int_t i )    { fNbinsQ = i; }
  void SetNbinsMass( Int_t i ) { fNbinsMass = i; }

  void SetMultMin( Double_t i ) { fMultMin = i; }
  void SetMultMax( Double_t i ) { fMultMax = i; }
  void SetPtMin( Double_t i )   { fPtMin = i; }
  void SetPtMax( Double_t i )   { fPtMax = i; }
  void SetPhiMin( Double_t i )  { fPhiMin = i; }
  void SetPhiMax( Double_t i )  { fPhiMax = i; }
  void SetEtaMin( Double_t i )  { fEtaMin = i; }
  void SetEtaMax( Double_t i )  { fEtaMax = i; }
  void SetQMin( Double_t i )    { fQMin = i; }
  void SetQMax( Double_t i )    { fQMax = i; }
  void SetMassMin( Double_t i ) { fMassMin = i; }
  void SetMassMax( Double_t i ) { fMassMax = i; }
  void SetHistWeightvsPhiMin( Double_t i ) {fHistWeightvsPhiMin=i;}
  void SetHistWeightvsPhiMax( Double_t i ) {fHistWeightvsPhiMax=i;}
  void SetCutTPC(Bool_t cut) {fCutTPC = cut;}
  // end setters common constants

  // setters for adding by hand flow values (afterburner)
  void SetAfterburnerOn(Bool_t b=kTRUE) {fAfterburnerOn=b;}
  void SetNonFlowNumberOfTrackClones(Int_t n) {fNonFlowNumberOfTrackClones=n;}
  void SetPtDifferentialV2( TF1 *gPtV2) {
    fDifferentialV2 = gPtV2;}
  void SetFlow( Double_t v1, Double_t v2, Double_t v3=0.0, Double_t v4=0.0, Double_t v5=0.0)
  {fV1=v1;fV2=v2;fV3=v3;fV4=v4;fV5=v5;}

  virtual void  SetDebugLevel(Int_t level) {fDebug = level;}
  void SetInput(int input) {fAnalysisInput = input;}
  void SetMCInput() {fIsMCInput = kTRUE;}
  void SetUseMCCen( Bool_t kB ) { fUseMCCen = kB; }
  void SetRejectPileUp( Bool_t kB ) { fRejectPileUp = kB; }
  void SetRejectPileUpTight( Bool_t kB ) { fRejectPileUpTight = kB; }
  void SetResetNegativeZDC( Bool_t kB ) { fResetNegativeZDC = kB; }
  void SetCorrectPhiTracklets( Bool_t kB ) { fCorrectPhiTracklets = kB; }
  void SetCentralityRange(Float_t centrlow=0., Float_t centrup=100.) {fCentrLowLim=centrlow;
    fCentrUpLim=centrup;}
  void SetCentralityEstimator(CentrEstimator centrest) {fCentrEstimator=centrest;}
  void SetDataSet(DataSet cDataSet) {fDataSet = cDataSet;}
  void SetZDCGainAlpha( Float_t a ) { fZDCGainAlpha = a; }
  void SetTowerEqList(TList* const kList) {this->fTowerEqList = (TList*)kList->Clone(); fUseTowerEq=kTRUE;};
  TList* GetTowerEqList() const {return this->fTowerEqList;};
  void SetBadTowerCalibList(TList* const kList) {this->fBadTowerCalibList = (TList*)kList->Clone(); fUseBadTowerCalib=kTRUE;};
  TList* GetBadTowerCalibList() const {return this->fBadTowerCalibList;};
  void SetVZEROGainEqList(TList* const kList) {this->fVZEROGainEqList = (TList*)kList->Clone();};
  TList* GetVZEROGainEqList() const {return this->fVZEROGainEqList;};
  void SetVZEROQVecRecList(TList* const kList) {this->fVZEROQVecRecList = (TList*)kList->Clone();};
  TList* GetVZEROQVecRecList() const {return this->fVZEROQVecRecList;};
  void SetZDCSpectraCorrList(TList* const kList) {this->fZDCSpectraCorrList = (TList*)kList->Clone(); fUseZDCSpectraCorr=kTRUE;};
  TList* GetZDCSpectraCorrList() const {return this->fZDCSpectraCorrList;};
  void SetFillZNCenDisRbR(Bool_t bFillZNCenDisRbR) {fFillZNCenDisRbR = bFillZNCenDisRbR;}; //@shi
  
  virtual Int_t GetCenBin(Double_t Centrality);
  Double_t GetWDist(const AliVVertex* v0, const AliVVertex* v1);
  Bool_t plpMV(const AliAODEvent* aod);
  Double_t GetBadTowerResp(Double_t Et, TH2D* BadTowerCalibHist);
  void SetWhichVZERORings(int minVZC, int maxVZC, int minVZA, int maxVZA) {fMinRingVZC = minVZC; fMaxRingVZC = maxVZC; fMinRingVZA = minVZA; fMaxRingVZA = maxVZA;}

  //@shi setters for the step of recentering for ZDC
  void SetStepZDCRecenter( Int_t a ) { fStepZDCRecenter = a; }
  void SetStoreCalibZDCRecenter( Bool_t b ) { fStoreCalibZDCRecenter = b; }
  void SetZDCCalibList(TList* const kList) {this->fZDCCalibList = (TList*)kList->Clone();};
  TList* GetZDCCalibList() const {return this->fZDCCalibList;};

  //@shi only for last check of recentering as the whole calib file is too large
  void SetZDCCalibListStep3CommonPart(TList* const kList) {this->fZDCCalibListStep3CommonPart = (TList*)kList->Clone();}; // run number independent calib hists
  TList* GetZDCCalibListStep3CommonPart() const {return this->fZDCCalibListStep3CommonPart;};
  void SetZDCCalibListStep3RunByRun(TList* const kList) {this->fZDCCalibListStep3RunByRun = (TList*)kList->Clone();}; // run number dependent
  TList* GetZDCCalibListStep3RunByRun() const {return this->fZDCCalibListStep3RunByRun;};
  
private:
  AliAnalysisTaskCRCZDC(const AliAnalysisTaskCRCZDC& dud);
  AliAnalysisTaskCRCZDC& operator=(const AliAnalysisTaskCRCZDC& dud);

  AnalysisType  fAnalysisType;      // can be MC, ESD or AOD
  TString       fRPType;            // can be Global or Tracklet or FMD
  AliCFManager* fCFManager1;        // correction framework manager
  AliCFManager* fCFManager2;        // correction framework manager
  AliFlowEventCuts* fCutsEvent;     //event cuts
  AliFlowTrackCuts* fCutsRP;        //cuts for RPs
  AliFlowTrackCuts* fCutsPOI;       //cuts for POIs
  TList*            fCutContainer;  //contains the cut objects
  TList*            fQAList;        // QA histogram list
  AliAnalysisUtils* fAnalysisUtil;  ///< Event selection
  Int_t         fMinMult;           // Minimum multiplicity from tracks selected using CORRFW
  Int_t         fMaxMult;           // Maximum multiplicity from tracks selected using CORRFW
  Double_t      fMinA;              // Minimum of eta range for subevent A
  Double_t      fMaxA;              // Maximum of eta range for subevent A
  Double_t      fMinB;              // Minimum of eta range for subevent B
  Double_t      fMaxB;              // Maximum of eta range for subevent B

  // mc event handlers
  AliGenEventHeader*        fGenHeader;       //!
  AliGenPythiaEventHeader*  fPythiaGenHeader; //!
  AliGenHijingEventHeader*  fHijingGenHeader; //!
  AliFlowTrack*             fFlowTrack;       //!

  Bool_t fQAon;                     // flag to set the filling of the QA hostograms
  Bool_t fLoadCandidates;           // true if reciving candidates collection

  // setters for common constants
  //histogram sizes
  Int_t  fNbinsMult; // histogram size
  Int_t  fNbinsPt;   // histogram size
  Int_t  fNbinsPhi;  // histogram size
  Int_t  fNbinsEta;  // histogram size
  Int_t  fNbinsQ;    // histogram size
  Int_t  fNbinsMass; // histogram size

  // Histograms limits
  Double_t  fMultMin;  // histogram limit
  Double_t  fMultMax;  // histogram limit
  Double_t  fPtMin;    // histogram limit
  Double_t  fPtMax;    // histogram limit
  Double_t  fPhiMin;   // histogram limit
  Double_t  fPhiMax;   // histogram limit
  Double_t  fEtaMin;   // histogram limit
  Double_t  fEtaMax;   // histogram limit
  Double_t  fQMin;     // histogram limit
  Double_t  fQMax;     // histogram limit
  Double_t  fMassMin;  // histogram limit
  Double_t  fMassMax;  // histogram limit
  Double_t fHistWeightvsPhiMin; //histogram limit
  Double_t fHistWeightvsPhiMax; //histogram limit
  // end common constants

  // Excluding a range
  Double_t  fExcludedEtaMin;  // excluded region limit
  Double_t  fExcludedEtaMax;  // excluded region limit
  Double_t  fExcludedPhiMin;  // excluded region limit
  Double_t  fExcludedPhiMax;  // excluded region limit
  // End of excluding a range

  // values afterburner
  Bool_t    fAfterburnerOn;              // do we afterburn?
  Int_t     fNonFlowNumberOfTrackClones; // number of times to clone the particles (nonflow)
  Double_t  fV1;        // Add Flow. Must be in range [0,0.5].
  Double_t  fV2;        // Add Flow. Must be in range [0,0.5].
  Double_t  fV3;        // Add Flow. Must be in range [0,0.5].
  Double_t  fV4;        // Add Flow. Must be in range [0,0.5].
  Double_t  fV5;        // Add Flow. Must be in range [0,0.5].
  TF1 *fDifferentialV2; // pt-differential v2

  AliFlowEvent* fFlowEvent; //flowevent
  Bool_t fShuffleTracks;    //serve the tracks shuffled

  TRandom3* fMyTRandom3;     // TRandom3 generator

  //******************************************************************************************************

  Int_t    fAnalysisInput;      // analysis input
  Bool_t   fIsMCInput;          // true when input is MC
  Bool_t   fUseMCCen;           // use GetZNCentroidInPbPb, with correction from MC
  Float_t  fCentrLowLim;	// centrality lower limit
  Float_t  fCentrUpLim;		// centrality upper limit
  CentrEstimator  fCentrEstimator;     // string for the centrality estimator
  Bool_t   fRejectPileUp;
  Bool_t   fRejectPileUpTight;
  Bool_t   fResetNegativeZDC;
  //
  TList       *fOutput;	   	//! list send on output slot 0
  //@Shi
  TList       *fOutputRecenter1; //! list send to output slot 1 to store 0-29th run numbers
  TList       *fOutputRecenter2; //! list send to output slot 2 to store 30-59th run numbers (if all stored together, max limit of a tlist is exceed)
  TList       *fOutputRecenter3; //! list send to output slot 3 to store the rest of run numbers (if all stored together, max limit of a tlist is exceed)

  //
  TH1F *fhZNCPM[5];		//! ZNC PM high res.
  TH1F *fhZNAPM[5];		//! ZNA PM high res.
  //
  TH1F *fhZNCPMQiPMC[4];	//! PMQi/PMC for ZNC
  TH1F *fhZNAPMQiPMC[4];	//! PMQi/PMC for ZNA
  //
  TH1F *fhZPCPM[5];		//! ZPC PM high res.
  TH1F *fhZPAPM[5];		//! ZPA PM high res.
  //
  TH1F *fhZPCPMQiPMC[4];	//! PMQi/PMC for ZPC
  TH1F *fhZPAPMQiPMC[4];	//! PMQi/PMC for ZPA
  //
  TH2F *fhZNCvsZNA;		//! ZNC vs ZNA;
  TH2F *fhZDCCvsZDCCA;		//! ZDCC vs ZDCCA
  TH2F *fhZNCvsZPC;		//! ZNC vs ZPC;
  TH2F *fhZNAvsZPA;		//! ZNA vs ZPA;
  TH2F *fhZNvsZP;		//! ZNC+ZNA vs ZPC+ZPA;
  TH2F *fhZNvsVZERO;		//! ZN vs VZERO;
  TH2F *fhZDCvsVZERO;		//! ZDC vs VZERO;
  TH3D *fhZNCenDis[2];		//! ZN centroid vs centrality

  //
  TH1F *fhAsymm;		//! ZN asymmetry
  TH2F *fhZNAvsAsymm;		//! ZNA vs asymmetry
  TH2F *fhZNCvsAsymm;		//! ZNC vs asymmetry
  //
  TH2F *fhZNCvscentrality;	//! ZNC vs. centrality
  TH2F *fhZNAvscentrality;	//! ZNA vs. centrality
  TH2F *fhZPCvscentrality;	//! ZNC vs. centrality
  TH2F *fhZPAvscentrality;	//! ZNA vs. centrality

  //@Shi define ZN and ZP corelation hists (begin)
  TH3D *fZPAvsZNASignal; //! ZPA vs. ZNA
  TH3D *fZPCvsZNCSignal; //! ZPC vs. ZNC
  TProfile *fZNenergyBeforeCalibration; //! test histogram for gain equalization
  TProfile *fZNenergyAfterCalibration; //! test histogram for gain equalization
  //@Shi define ZN and ZP corelation hists (end)
  
  TH3D *fhZNSpectra;   	//! ZNA vs. centrality
  TH3D *fhZNSpectraCor;   	//! ZNA vs. centrality
  TH3D *fhZNSpectraPow;   	//! ZNA vs. centrality
  TH3D *fhZNBCCorr;   	//! ZNA vs. centrality
  
  //@Shi add fhZPSpectra
  TH3D *fhZPSpectra;   	//! ZP vs. centrality

  TH3D *fQATrackTPCNcls; //!
  TH3D *fQATrackITSNcls; //!
  TH3D *fQATrackTPCchi2; //!
  TH3D *fQATrackITSchi2; //!
  TH3D *fQATrackTPCScls; //!
  TH3D *fQATrackITSScls; //!

  const static Int_t fCRCMaxnRun = 211;
  
  //@Shi add run by run ZN centroid vs centrality
  TH3D *fhZNCenDisRbR[fCRCMaxnRun][2];      //! ZN centroid vs centrality run by run

  //@Shi Add run by run recentering histograms for ZDC (begin)
  Int_t fStepZDCRecenter = -1;
  Bool_t fStoreCalibZDCRecenter = kFALSE;
  TList *fZDCCalibList; //
  TList *fZDCCalibListStep3CommonPart; //
  TList *fZDCCalibListStep3RunByRun; //
  TArrayD fAvVtxPosX;    // Run list
  TArrayD fAvVtxPosY;    // Run list
  TArrayD fAvVtxPosZ;    // Run list
  Double_t fVtxPosCor[3];
  TProfile *fAve_VtxX; //!
  TProfile *fAve_VtxY; //!
  TProfile *fAve_VtxZ; //!
  TProfile *fRun_VtxXQPreCalib[fCRCMaxnRun][4];  //!
  TProfile *fRun_VtxYQPreCalib[fCRCMaxnRun][4];  //!
  TProfile *fRun_VtxZQPreCalib[fCRCMaxnRun][4];  //!
  TProfile *fRun_VtxXQCalibStep1[fCRCMaxnRun][4];  //!
  TProfile *fRun_VtxYQCalibStep1[fCRCMaxnRun][4];  //!
  TProfile *fRun_VtxZQCalibStep1[fCRCMaxnRun][4];  //!
  TProfile *fRun_VtxXQCalibStep2[fCRCMaxnRun][4];  //!
  TProfile *fRun_VtxYQCalibStep2[fCRCMaxnRun][4];  //!
  TProfile *fRun_VtxZQCalibStep2[fCRCMaxnRun][4];  //!
  TProfile *fRun_CentQCalib[fCRCMaxnRun][4];  //!
  TProfile *fRun_VtxXQCalib[fCRCMaxnRun][4];  //!
  TProfile *fRun_VtxYQCalib[fCRCMaxnRun][4];  //!
  TProfile *fRun_VtxZQCalib[fCRCMaxnRun][4];  //!
  TProfile *fRun_CentQCalib2[fCRCMaxnRun][4];  //!
  
  TProfile *fCorrQAReCRe;  //!
  TProfile *fCorrQAReCIm;  //!
  TProfile *fCorrQAImCRe;  //!
  TProfile *fCorrQAImCIm;  //!
  TProfile *fRun_CentQ[fCRCMaxnRun][4];  //!
  TProfile3D *fRun_VtxXYZQ[fCRCMaxnRun][4];  //!
  const static Int_t fnCentBinForRecentering = 20; // this means that a wider centrality bin is used {0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100}
  TProfile3D *fCent_VtxXYZQ[fnCentBinForRecentering][4];  //!
  
  TProfile* fAvr_Run_CentQ[4]; //!
  TProfile3D *fAvr_Run_VtxXYZQ[4]; //!
  TProfile3D *fAvr_Cent_VtxXYZQ[4]; //!
  //@Shi Add run by run recentering histograms for ZDC (end)
  
//  TH3D *fhZNSpectraRbR[fCRCMaxnRun]; //! ZNA vs. centrality

  const static Int_t fCRCnTow = 5;
  const static Int_t fnCen = 10;
  Int_t fCRCnRun;
  Float_t fZDCGainAlpha;
  DataSet fDataSet;
  Int_t fRunList[fCRCMaxnRun];                   //! Run list
//  TProfile *fhnTowerGain[fCRCnTow]; //! towers gain
  TList *fCRCQVecListRun[fCRCMaxnRun]; //! Q Vectors list per run
  //@Shi add fRecenter1ListRunbyRun and fRecenter2ListRunbyRun and fRecenter3ListRunbyRun
  TList *fRecenter1ListRunbyRun[fCRCMaxnRun]; //! Recenter list run by run 1
  TList *fRecenter2ListRunbyRun[fCRCMaxnRun]; //! Recenter list run by run 2
  TList *fRecenter3ListRunbyRun[fCRCMaxnRun]; //! Recenter list run by run 3
  TProfile *fZNCTower[fCRCMaxnRun][fCRCnTow];		//! ZNC tower spectra
  TProfile *fZNATower[fCRCMaxnRun][fCRCnTow];		//! ZNA tower spectra
  //@Shi add fZPCTower and fZPATower
  TProfile *fZPCTower[fCRCMaxnRun][fCRCnTow];		//! ZPC tower spectra
  TProfile *fZPATower[fCRCMaxnRun][fCRCnTow];		//! ZPA tower spectra
  TClonesArray* fStack; //!
  TList *fSpectraMCList;   //! list with pt spectra
  TH1F *fPtSpecGen[2][10];		//! PtSpecGen
  TH1F *fPtSpecFB32[2][10];		//! PtSpecRec FB32
  TH1F *fPtSpecFB96[2][10];		//! PtSpecRec FB96
  TH1F *fPtSpecFB128[2][10];  //! PtSpecRec FB128
  TH1F *fPtSpecFB768[2][10];  //! PtSpecRec FB768
  Bool_t fCutTPC;
  TH1F *fCenDis; //! centrality distribution
  TH1F *fPileUpCount; //! centrality distribution
  TH1F *fPileUpMultSelCount; //! centrality distribution
  TF1 *fMultTOFLowCut; //!
  TF1 *fMultTOFHighCut; //!
  TProfile2D *fVZEROMult; //!
  TH1F *fRecordNegativeEZNA; //! record # of event with negative EZNA value which causes Power(EZNA, fZDCGainAlpha) problem (Shi)
  TH1F *fRecordNegativeEZNC; //! record # of event with negative EZNC value which causes Power(EZNC, fZDCGainAlpha) problem (Shi)

  Bool_t fFillZNCenDisRbR; //@Shi flag for run by run ZN centroid distribution. Do not turn on for large dataset when running on grid. It takes too much memory 
  AliMultSelection* fMultSelection; //! MultSelection (RUN2 centrality estimator)
  Bool_t fUseTowerEq; //
  TList *fTowerEqList;   // list with weights
  TH1D *fTowerGainEq[2][5]; //!
  TH1D *fHZDCCparameters; //! ZDC-C gain eq and recentering parameters using least square fit for 2018 dataset
  TH1D *fHZDCAparameters; //! ZDC-A gain eq and recentering parameters using least square fit for 2018 dataset
  TList *fBadTowerStuffList; //! list for storing calib files
  Bool_t fUseBadTowerCalib; //
  TList *fBadTowerCalibList; // list with original calib files
  TList *fVZEROGainEqList; //
  TList *fVZEROQVecRecList; //
  TList *fVZEROStuffList; //!
  Bool_t fUseZDCSpectraCorr; //
  TList *fZDCSpectraCorrList; //
  TH1D *SpecCorMu1[8]; //!
  TH1D *SpecCorMu2[8]; //!
  TH1D *SpecCorAv[8]; //!
  TH1D *SpecCorSi[8]; //!
  TH2D *fBadTowerCalibHist[100]; //!
  TH2D *fVZEROGainEqHist; //!
  TH2F *fHCorrectV0ChWeghts; //! V0 gain eq weight file for 2018 dataset  
  const static Int_t fkVZEROnHar = 4;
//  TProfile3D *fVZEROQVectorRecQx[fkVZEROnHar]; //!
//  TProfile3D *fVZEROQVectorRecQy[fkVZEROnHar]; //!
  TProfile3D *fVZEROQVectorRecQxStored[fkVZEROnHar]; //!
  TProfile3D *fVZEROQVectorRecQyStored[fkVZEROnHar]; //!
  const static Int_t fkVZEROnQAplots = 8;
  TProfile2D *fVZEROQVectorRecFinal[fkVZEROnHar][fkVZEROnQAplots]; //!
  Int_t fCachedRunNum;   //
  Int_t fMinRingVZC; //
  Int_t fMaxRingVZC; //
  Int_t fMinRingVZA; //
  Int_t fMaxRingVZA; //
  Bool_t fCorrectPhiTracklets; //

  // TrackQA
  TList *fTrackQAList; //!
  const static Int_t fKNFBs = 4;
  TH3D* fTrackQADCAxy[fKNFBs][4]; //!
  TH3D* fTrackQADCAz[fKNFBs][4]; //!
  TH2D* fTrackQApT[fKNFBs][4]; //!
  TProfile2D* fTrackQADphi[fKNFBs][4]; //!
  TH2D* fEbEQRe[fKNFBs][4]; //!
  TH2D* fEbEQIm[fKNFBs][4]; //!
  TH2D* fEbEQMu[fKNFBs][4]; //!

  // Functions for Pile Up Event Removal 2018 period:
  TF1                   *fV0CutPU;      //
  TF1                   *fSPDCutPU;     //
  TF1                   *fMultCutPU;    //
  TF1                   *fCenCutLowPU;  //
  TF1                   *fCenCutHighPU; //
  
  ClassDef(AliAnalysisTaskCRCZDC,12);

};

#endif
