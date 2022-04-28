#ifndef AliAnalysisTaskCMWESETrkSyst_cxx
#define AliAnalysisTaskCMWESETrkSyst_cxx
#include "AliAnalysisTaskSE.h"
#include "AliEventCuts.h"
#include "AliGFWCuts.h"
#include "AliGFWWeights.h"
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

class AliAnalysisTaskCMWESETrkSyst : public AliAnalysisTaskSE
{
public:

  AliAnalysisTaskCMWESETrkSyst();
  AliAnalysisTaskCMWESETrkSyst(const char *name);
  AliAnalysisTaskCMWESETrkSyst(const char *name, TString PR, bool NUE, bool NUA, bool V0Calib, bool MCSyst, bool FBSyst);

  virtual ~AliAnalysisTaskCMWESETrkSyst();

  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *);

  int     GetDebug(){return fDebug;}
  void SetDebug(int x){fDebug = x;}

  int     GerHarmonic(){return fHarmonic;}
  void  SetHarmonic(double x)  {fHarmonic = x;}

  TString GetTrigger(){return fTrigger;}
  void SetTrigger(TString x){fTrigger = x;}

  int    GetFilterBit(){return fFltbit;}
  void SetFilterBit(int x){fFltbit = x;}

  int    GetFilterBitSyst(){return fFltbitSyst;}
  void SetFilterBitSyst(int x){fFltbitSyst = x;}

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

  bool GetV0CalibOn(){return fV0CalibOn;}
  void SetV0CalibOn(bool x){fV0CalibOn = x;}

  bool GetV0QAOn(){return fQAV0;}
  void SetV0QAOn(bool x){fQAV0 = x;}

  bool GetNUEOn(){return fDoNUE;}
  void SetNUEOn(bool x){fDoNUE = x;}

  bool GetNUAOn(){return fDoNUA;}
  void SetNUAOn(bool x){fDoNUA = x;}

  bool GetMCSystOn(){return fMCSyst;}
  void SetMCSystOn(bool x){fMCSyst = x;}

  bool GetFBSystOn(){return fFBSyst;}
  void SetFBSystOn(bool x){fFBSyst = x;}
  
private:

  static const int NCENTBINS = 10;
  static const int NPOIBINS  = 2;
  static const int NQNBINS  = 10;
  static const int NRUNNUM = 138;
  void          ResetHists();
  bool          DoCumulants();
  void          CalcRefFlow();
  void          CalcIntCov();
  double      GetNUECor(int charge, double pt, int SystSwitch);
  double      GetNUACor(int charge, double phi, double eta, double vz, int SystSwitch);
  bool          LoadEmilNUAWeights(const Int_t &lRunNo);
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
  bool          AcceptAODTrack(AliAODEvent* fAOD, AliAODTrack *track, AliAODVertex* fVtx, int SystSwitch);
  bool          AnalyzeAOD(AliAODEvent* fAOD, AliAODVertex* fVtx, int SystSwitch);
  bool          CalcQnVectorV0(AliAODEvent* fAOD, AliAODVertex* fVtx, int SystSwitch);
  bool          CalcQnVectorTPC(AliAODEvent* fAOD);
  int             GetQnPercV0(double qn_thtsEvt, int V0lb);
  bool          GetV0CalibHisto(AliAODEvent* fAOD, AliAODVertex* fVtx);
  int             GetQnPercTPC(AliAODEvent* fAOD);
  double      GetEventPlane(double qx, double qy);
  int             GetPercCode(double perc);
  void          OpenInfoCalbration(Int_t run);
    // Cuts and options
    AliGFWCuts*    fGFWSelection;    
    int                     fDebug; // debug level controls amount of output statements
    double              fHarmonic; // value of harmonic
    TString             fTrigger; // flag of trigger
    int                     fFltbit; // AOD filter bit selection
    int                     fFltbitSyst; // FB for systematic
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
    bool                  fQAV0; // flag for V0  qn QA
    bool                  fDoNUE; // switch for NUE
    bool                  fDoNUA; // switch for NUA
    bool                  fMCSyst;
    bool                  fFBSyst;
    // Global Variables Unchanged in an Evt
    int                     fRunNum; // runnumber
    int                     fRunNumBin; // runnumer bin; 10:139510...; 11:170387...; 15HIR:246994...
    int                     fVzBin; // vertex z bin
    int                     fCentBinV0M; // centrality bin: 0-10
    int                     fCentBinSPD1; // centrality bin: 0-10
    double              fCentV0M; // value of centrality 
    double              fCentSPD1; // value of centrality 
    int                     fQnBinV0C; // qn bin: 0-10
    int                     fQnBinV0A; // qn bin: 0-10
    const float        fEtaCut; // eta cut
    const float        fZvtxCut; // z-vertex selection for collision  
    // Weight List   
    TList*                fListNUE; // read list for NUE
    TList*                fListNUA; // read lists for NUA
    TList*                fListVZEROCALIB; // read list fpr V0 Calib
    TList*                fListNUESyst;
    TList*                fListNUESystFB;
    TList*                fListNUASystFB;
    AliGFWWeights *fWeights;//! This should be stored in TList

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
    TH2D*             hQnV0CPercentile;
    TH2D*             hQnV0APercentile;
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
    TH1D*             hPtMCSyst;
    TH2D*             hPDedxMCSyst;
    TH1D*             hPtFBSyst;
    TH2D*             hPDedxFBSyst;

    // Update Evt-by-Evt, will not be saved
    TH2D*             hReQ_thisEvt;
    TH2D*             hImQ_thisEvt;
    TH2D*             hMQ_thisEvt;
    TH2D*             hMQ_weight_thisEvt;
    TProfile*          pRefFlow_thisEvt;
    TProfile*          pIntd2q2c_thisEvt;
    TProfile*          pIntd2q2a_thisEvt;

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
    TSpline3*        splQ2c[80]; //A.Dobrin
    TSpline3*        splQ2a[80]; //A.Dobrin
    // Output QA
    TH1D*             hCentV0M[2];
    TH1D*             hCentSPD1[2];
    TH1D*             hVz[2];
    TH2D*             hCentQA[8];
    TH2D*             hMultCentQA[2];
    TH2D*             hMultMultQA[6];

    // track-wise QA
    TH1D*             hEta[2];
    TH1D*             hPhi[2];
    TH1D*             hDcaXy[2];
    TH1D*             hDcaZ[2];
    TH1D*             hNhits[2];

    TH1D*             hEtaMCSyst[2];
    TH1D*             hPhiMCSyst[2];
    TH1D*             hDcaXyMCSyst[2];
    TH1D*             hDcaZMCSyst[2];
    TH1D*             hNhitsMCSyst[2];

    TH1D*             hEtaFBSyst[2];
    TH1D*             hPhiFBSyst[2];
    TH1D*             hDcaXyFBSyst[2];
    TH1D*             hDcaZFBSyst[2];
    TH1D*             hNhitsFBSyst[2];

    // Qn & psi QA
    TH2D*             hQxCentRecenter[3];
    TH2D*             hQxVtxRecenter[3];
    TH2D*             hQyCentRecenter[3];
    TH2D*             hQyVtxRecenter[3];
    TH2D*             hQnCentRecenter[3];
    TH1D*             hPsiV0Recenter[NCENTBINS][3];

    // Syst 1 (MC)
    TProfile*         pRefFlowDef[NCENTBINS];     
    TProfile*         pIntd2Def[NCENTBINS];
    TProfile*         pIntd2AchDef[NCENTBINS];
    TProfile*         pAchDef[NCENTBINS];

    // Syst 1 (MC)
    TProfile*         pRefFlowMCSyst[NCENTBINS];     
    TProfile*         pIntd2MCSyst[NCENTBINS];
    TProfile*         pIntd2AchMCSyst[NCENTBINS];
    TProfile*         pAchMCSyst[NCENTBINS];

    // Syst 2 (FB)
    TProfile*         pRefFlowFBSyst[NCENTBINS];     
    TProfile*         pIntd2FBSyst[NCENTBINS];
    TProfile*         pIntd2AchFBSyst[NCENTBINS];
    TProfile*         pAchFBSyst[NCENTBINS];

    AliAnalysisTaskCMWESETrkSyst(const AliAnalysisTaskCMWESETrkSyst&);
    AliAnalysisTaskCMWESETrkSyst& operator=(const AliAnalysisTaskCMWESETrkSyst&);

    ClassDef(AliAnalysisTaskCMWESETrkSyst, 1);
  };

#endif