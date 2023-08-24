#ifndef AliAnalysisTaskCMWESE_cxx
#define AliAnalysisTaskCMWESE_cxx
#include "AliAnalysisTaskSE.h"
#include "AliEventCuts.h"
#include "AliAODTrack.h"
#include "TProfile3D.h"
#include "TProfile2D.h"
#include "TProfile.h"
#include "TComplex.h"
#include "TList.h"
#include "TFile.h"
#include "TSpline.h"
#include "TProfile.h"
#include "TH1.h"
#include "TH2.h"

//class AliAnalysisUtils;

#include "AliAnalysisTaskSE.h"
#include "AliEventCuts.h"

class AliAnalysisTaskCMWESE : public AliAnalysisTaskSE
{
public:

  AliAnalysisTaskCMWESE();
  AliAnalysisTaskCMWESE(const char *name);
  AliAnalysisTaskCMWESE(const char *name, TString PR, bool NUE, bool NUA, bool V0Calib);

  virtual ~AliAnalysisTaskCMWESE();

  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *);

  int     GetDebug(){return fDebug;}
  void SetDebug(int x){fDebug = x;}

  int     GerHarmonic(){return fHarmonic;}
  void  SetHarmonic(double x)  {fHarmonic = x;}

  TString GetTrigger(){return fTrigger;}
  void SetTrigger(TString x){fTrigger = x;}

  TString GetCentEst(){return fCentEst;}
  void SetCentEst(TString x){fCentEst = x;}

  int    GetFilterBit(){return fFltbit;}
  void SetFilterBit(int x){fFltbit = x;}

  int    GetNclsCut(){return fNclsCut;}
  void SetNclsCut(int x){fNclsCut = x;}

  float GetChi2High(){return fChi2Hg;}
  void SetChi2High(float x){fChi2Hg = x;}

  float GetChi2Low(){return fChi2Lo;}
  void SetChi2Low(float x){fChi2Lo = x;}

  float GetDCAcutZ(){return fDcaCutz;}
  void SetDCAcutZ(float x){fDcaCutz = x;}

  float GetDCAcutXY(){return fDcaCutxy;}
  void SetDCAcutXY(float x){fDcaCutxy = x;}

  float GetPtMin(){return fPtMin;}
  void SetPtMin(float x){fPtMin = x;}
  
  float GetPtMax(){return fPtMax;}
  void SetPtMax(float x){fPtMax = x;}

  int    GetCentBinLow(){return fCbinLo;}
  void SetCentBinLow(int x){fCbinLo = x;}
  
  int    GetCentBinHigh(){return fCbinHg;}
  void SetCentBinHigh(int x){fCbinHg = x;}

  TString GetPeriod(){return fPeriod;}
  void   SetPeriod(TString x) { fPeriod = x; }

  TString GetMultComp(){return fMultComp;}
  void   SetMultComp(TString x) { fMultComp = x; }

  float   GetEtaGap(){return fEtaGap;}
  void   SetEtaGap(float x) { fEtaGap = x; }

  bool  GetV0CalibOn(){return fV0CalibOn;}
  void SetV0CalibOn(bool x){fV0CalibOn = x;}

  bool   GetTPCCalibOn(){return fTPCCalibOn;}
  void SetTPCCalibOn(bool x){fTPCCalibOn = x;}

  bool GetV0QAOn(){return fQAV0;}
  void SetV0QAOn(bool x){fQAV0 = x;}

  bool GetTPCQAOn(){return fQATPC;}
  void SetTPCQAOn(bool x){fQATPC = x;}

  bool GetNUEOn(){return fDoNUE;}
  void SetNUEOn(bool x){fDoNUE = x;}

  bool GetNUAOn(){return fDoNUA;}
  void SetNUAOn(bool x){fDoNUA = x;}
  
  float GetCentCut(){return fCentCut;}
  void SetCentCut(float x){fCentCut = x;}

  float GetVzCut(){return fZvtxCut;}
  void SetVzCut(float x){fZvtxCut = x;}

  int    GetPUSyst(){return fPUSyst;}
  void SetPUSyst(int x){fPUSyst = x;}

private:

  static const int NCENTBINS = 10;
  static const int NPOIBINS  = 2;
  static const int NQNBINS  = 10;

  void          ResetHists();
  bool          DoCumulants();
  void          CalcRefFlow();
  void          CalcIntCov();
  double      GetNUECor(int charge, double pt);
  double      GetNUACor(int charge, double phi, double eta, double vz);
  int             GetRunNumBin(int runNum);
  // pile-up
  bool          RejectEvtMultComp      (AliAODEvent* fAOD);
  bool          RejectEvtTFFit(AliAODEvent* fAOD);
  bool          RejectEvtTPCITSfb32TOF       (AliAODEvent* fAOD);
  bool          AODPileupCheck (AliAODEvent* fAOD);
  bool          PileUpMultiVertex (AliAODEvent* fAOD);
  double      GetWDist(const AliVVertex* v0, const AliVVertex* v1);
  bool          RemovalForLHC18           (AliAODEvent* fAOD);
  bool          RemovalForRun1 (AliAODEvent* fAOD, AliAnalysisUtils* fUtils);
  bool          RemovalForpPb          (AliAODEvent* fAOD);
  bool          AcceptAODTrack(AliAODEvent* fAOD, AliAODTrack *track, AliAODVertex* fVtx);
  bool          AnalyzeAOD(AliAODEvent* fAOD, AliAODVertex* fVtx);
  bool          CalcQnVectorV0(AliAODEvent* fAOD, AliAODVertex* fVtx, double mAch, double Mult);
  bool          CalcQnVectorTPC(AliAODEvent* fAOD);
  int             GetQnPercV0(double qn_thtsEvt, double mAch, double Mult);
  bool          GetV0CalibHisto(AliAODEvent* fAOD, AliAODVertex* fVtx);
  int             GetQnPercTPC(AliAODEvent* fAOD);
  double      GetEventPlane(double qx, double qy);
  int             GetPercCode(double perc);
  void          OpenInfoCalbration(Int_t run);
  void          OpenInfoCalbration18(Int_t run);

    // Cuts and options
    int                     fDebug; // debug level controls amount of output statements
    double              fHarmonic; // value of harmonic
    TString             fTrigger; // flag of trigger; 0 = kINT7; 1 = kMB; 2 = kMB+kCentral+kSemiCentral
    TString             fCentEst;
    int                     fFltbit; // AOD filter bit selection
    int                     fNclsCut; // ncls cut for all tracks 
    float                  fChi2Hg; // upper limmit for chi2
    float                  fChi2Lo; // lower limmit for chi2
    float                  fDcaCutz; // dcaz cut for all tracks
    float                  fDcaCutxy; // dcaxy cut for all tracks
    float                  fPtMin; // minimum pt for Q-vector components
    float                  fPtMax; // maximum pt for Q-vector components
    int                     fCbinHg; // higher centrality bin for histogram array
    int                     fCbinLo; // lower centrality bin for histogram array
    TString             fPeriod; // period
    TString             fMultComp; // Multiplicity Comparison pile-up
    double              fEtaGap; // value for the Eta Gap in c2 calculation
    bool                  fV0CalibOn; // switch for v0 qn calib
    bool                  fTPCCalibOn; // switch for tpc qn calib
    bool                  fQAV0; // flag for V0  qn QA
    bool                  fQATPC; // flag for TPC qn QA
    bool                  fDoNUE; // switch for NUE
    bool                  fDoNUA; // switch for NUA
    float                  fCentCut; // centrality restriction for V0M and TRK
    // Global Variables Unchanged in an Evt
    int                     fRunNum; // runnumber
    int                     fRunNumBin; // runnumer bin; 10:139510...; 11:170387...; 15HIR:246994...
    int                     fVzBin; // vertex z bin
    int                     fCentBin; // centrality bin: 0-10
    double              fCent; // value of centrality 
    int                     fQnBin; // qn bin: 0-10
    const float        fEtaCut; // eta cut
    const float        fDedxCut; //dedx cut
    float                  fZvtxCut; // z-vertex selection for collision  
    int                     fPUSyst; // 0: default; 1: Tight PU cut; 2: TOF PU cut; 3: Fill PU effect QA; 4: with new Para       
    // Weight List   
    TList*                fListNUE; // read list for NUE
    TList*                fListNUA; // read lists for NUA
    TList*                fListVZEROCALIB; // read list fpr V0 Calib
    TList*                fListwCent;
    TH1D*               hwCent;
    // Q QStar event-wise
    TComplex        fNegEtaQ;
    TComplex        fNegEtaQStar;
    TComplex        fPosEtaQ;
    TComplex        fPosEtaQStar;
    double              fNegEtaMQ;
    double              fPosEtaMQ;

    // Read Files for NUA/NUE/VoCalib
    TH1D*             hNUEweightPlus;
    TH1D*             hNUEweightMinus;
    TH2D*             hNUAweightPlus;
    TH2D*             hNUAweightMinus;
    TH3F*              hCorrectNUAPos; // Protty
    TH3F*              hCorrectNUANeg; // Protty
    TH2D*             hMultV0Read;
    TH2F*             fHCorrectV0ChWeghts;
    TH2D*             hQnPercentile;
    TH1D*             hQnPercentile_centThisEvt;
    TSpline3*        sp; 
    TF1*                fSPDCutPU;
    TF1*                fV0CutPU;
    TF1*                fCenCutLowPU;
    TF1*                fCenCutHighPU;
    TF1*                fMultCutPU;
    // Output QA
    TList*               fOutputList;          
    TH1D*             hEvtCount;   
    TH1I*               hRunNumBin;  
    TH1D*             hPt;
    TH2D*             hPDedx;

    // Update Evt-by-Evt, will not be saved
    TH2D*             hReQ_thisEvt;
    TH2D*             hImQ_thisEvt;
    TH2D*             hReQ2_thisEvt;
    TH2D*             hImQ2_thisEvt; 
    TH2D*             hMQ_thisEvt;
    TH2D*             hMQ_weight_thisEvt;
    TH2D*             hReQPos_thisEvt;
    TH2D*             hImQPos_thisEvt;
    TH2D*             hMQPos_thisEvt;
    TH2D*             hMQPos_weight_thisEvt;
    TH2D*             hReQNeg_thisEvt;
    TH2D*             hImQNeg_thisEvt;
    TH2D*             hMQNeg_thisEvt;
    TH2D*             hMQNeg_weight_thisEvt;
    TProfile*          pRefFlow_thisEvt;
    TProfile*          pIntd2_thisEvt;
    TProfile*          pbufferPt_thisEvt;
    TProfile*          pbufferEta_thisEvt;

    // Read Files for V0Calib
    TProfile3D*     pV0XMeanRead[3]; 
    TProfile3D*     pV0YMeanRead[3];
    // Run2 A.Dorbin
    TH1D*             hMultV0; //Dobrin
    TH1D*             hQxnmV0[2];
    TH1D*             hQynmV0[2];
    double             fMultV0Ch[64];
    double             fV0XMean[3];
    double             fV0YMean[3];
    TSpline3*        splQ2c[90]; //A.Dobrin

    // Output QA
    TH1D*             hCent[2];
    TH1D*             hCentRBR[138];
    TH1D*             hVz[2];
    TH2D*             hCentQA[8];
    TH2D*             hMultCentQA[2];
    TH2D*             hMultMultQA[6];

    // track-wise QA
    TH1D*             hEta[2];
    TH1D*             hPhi[2];
    TH2D*             hEtaPhi[2];
    TH1D*             hDcaXy[2];
    TH1D*             hDcaZ[2];
    TH1D*             hNhits[2];

    // Qn & psi QA
    TH2D*             hQxCentRecenter[3];
    TH2D*             hQxVtxRecenter[3];
    TH2D*             hQyCentRecenter[3];
    TH2D*             hQyVtxRecenter[3];
    TH2D*             hQnCentRecenter[3];
    TH1D*             hPsiV0Recenter[NCENTBINS][3];

    // physics 
    TProfile*         pRefFlow[NCENTBINS];     
    TProfile*         pIntd2[NCENTBINS];
    TProfile*         pIntd2Ach[NCENTBINS];
    TProfile*         pAch[NCENTBINS];
    TH1D*            hMultMb[NQNBINS];
    TH1D*            hAchMb[NQNBINS];

    // pt,eta vs Ach
    TProfile* pPosPtAch[NCENTBINS];
    TProfile* pNegPtAch[NCENTBINS];
    TProfile* pPosEtaAch[NCENTBINS];
    TProfile* pNegEtaAch[NCENTBINS];

    AliAnalysisTaskCMWESE(const AliAnalysisTaskCMWESE&);
    AliAnalysisTaskCMWESE& operator=(const AliAnalysisTaskCMWESE&);

    ClassDef(AliAnalysisTaskCMWESE, 1);
  };

#endif