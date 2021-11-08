#ifndef AliAnalysisTaskEPCalib_cxx
#define AliAnalysisTaskEPCalib_cxx

//class TList;
//class TH1F;
//class TH2F;
//class TProfile;
//class AliAnalysisUtils;

//#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskEPCalib : public AliAnalysisTaskSE
{
 public:

  AliAnalysisTaskEPCalib();
  AliAnalysisTaskEPCalib(const char *name);
  virtual ~AliAnalysisTaskEPCalib();

  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *);

  bool    GetTPCEstOn(){return fTPCEstOn;}
  void    SetTPCEstOn(bool x){fTPCEstOn = x;}

  bool    GetTPCNUAWeight(){return fTPCNUAWeight;}
  void    SetTPCNUAWeight(bool x){fTPCNUAWeight = x;}

  bool    GetFillTPCQMean(){return fFillTPCQMean;}
  void    SetFillTPCQMean(bool x){fFillTPCQMean = x;}

  bool    GetFillTPCShift(){return fFillTPCShift;}
  void    SetFillTPCShift(bool x){fFillTPCShift = x;}

  bool    GetTPCCalib(){return fTPCCalib;}
  void    SetTPCCalib(bool x){fTPCCalib = x;}

  bool    GetVZEROEstOn(){return fVZEROEstOn;}
  void    SetVZEROEstOn(bool x){fVZEROEstOn = x;}

  bool    GetVZEROGainEq(){return fVZEROGainEq;}
  void    SetVZEROGainEq(bool x){fVZEROGainEq = x;}

  bool    GetFillVZEROQMean(){return fFillVZEROQMean;}
  void    SetFillVZEROQMean(bool x){fFillVZEROQMean = x;}

  bool    GetVZEROCalib(){return fVZEROCalib;}
  void    SetVZEROCalib(bool x){fVZEROCalib = x;}

  bool    GetfQAV0(){return fQAV0;}
  void    SetfQAV0(bool x){fQAV0 = x;}

  int       GetDebug(){return fDebug;}
  void     SetDebug(int x){fDebug = x;}
  
  double  GerHarmonic(){return fHarmonic;}
  void      SetHarmonic(double x)  {fHarmonic = x;}

  int       GetTrigger(){return fTrigger;}
  void    SetTrigger(int x){fTrigger = x;}

  int       GetFilterBit(){return fFltbit;}
  void    SetFilterBit(int x){fFltbit = x;}
    
  int       GetNclsCut(){return fNclsCut;}
  void    SetNclsCut(int x){fNclsCut = x;}
    
  float    GetChi2High(){return fChi2Hg;}
  void    SetChi2High(float x){fChi2Hg = x;}
    
  float    GetChi2Low(){return fChi2Lo;}
  void    SetChi2Low(float x){fChi2Lo = x;}
    
  float    GetDCAcutZ(){return fDcaCutz;}
  void    SetDCAcutZ(float x){fDcaCutz = x;}
    
  float    GetDCAcutXY(){return fDcaCutxy;}
  void    SetDCAcutXY(float x){fDcaCutxy = x;}
    
  float    GetPtMin(){return fPtMin;}
  void    SetPtMin(float x){fPtMin = x;}
      
  float    GetPtMax(){return fPtMax;}
  void    SetPtMax(float x){fPtMax = x;}
    
  int       GetCentBinLow(){return fCbinLo;}
  void    SetCentBinLow(int x){fCbinLo = x;}
      
  int       GetCentBinHigh(){return fCbinHg;}
  void    SetCentBinHigh(int x){fCbinHg = x;}

  TString GetPeriod(){return fPeriod;}
  void      SetPeriod(TString x) { fPeriod = x; }

  TString GetMultComp(){return fMultComp;}
  void      SetMultComp(TString x) { fMultComp = x; }
  
  float      GetCentCut(){return fCentCut;}
  void      SetCentCut(float x){fCentCut = x;}

 private:

  static const int NCENTBINS = 10;
  static const int NRUNNUM=68;

  int                     GetRunNumBin(int runNum);
  // pile-up        
  bool                  RejectEvtTFFit(AliAODEvent* fAOD);
  bool                  RemovalForRun1 (AliAODEvent* fAOD, AliAnalysisUtils* fUtils);
  bool                  AcceptAODTrack(AliAODEvent* fAOD, AliAODTrack *track, AliAODVertex* fVtx);
  void                  TPCPlane(AliAODEvent* fAOD);
  void                  V0Plane(AliAODEvent* fAOD);
  double              GetEventPlane(double qx, double qy);
  double              GetNUACor(int charge, double phi, double eta, double vz);

  bool                  fTPCEstOn;
  bool                  fTPCNUAWeight;
  bool                  fFillTPCQMean;
  bool                  fFillTPCShift;
  bool                  fTPCCalib ;
  bool                  fVZEROEstOn ;
  bool                  fVZEROGainEq;
  bool                  fFillVZEROQMean;
  bool                  fVZEROCalib;
  bool                  fQAV0;

  // Cuts and options
  int                     fDebug; // debug level controls amount of output statements
  double              fHarmonic;
  int                     fTrigger; // flag of trigger; 0 = kINT7; 1 = kMB; 2 = kMB+kCentral+kSemiCentral
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
  float                  fCentCut; // centrality restriction for V0M and TRK
  // Global Variables Unchanged in an Evt
  int                     fRunNum; // runnumber
  int                     fRunNumBin; // runnumer bin; 10:139510...; 11:170387...; 15HIR:246994...
  int                     fVzBin; // vertex z bin
  int                     fCentBin; // centrality bin: 0-10
  double              fCent; // value of centrality 
  double              fCentSPD;
  const float        fEtaCut; // eta cut
  const float        fDedxCut; //dedx cut
  const float        fZvtxCut; // z-vertex selection for collision  
  TList*               fListNUA1;
  TList*               fListNUA2;
  TList*               fListNUA3;
  TH2D*             hNUAweightPlus;
  TH2D*             hNUAweightMinus;
  TH3F*              hCorrectNUAPos;
  TH3F*              hCorrectNUANeg;
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

  TString            fRunNumList[NRUNNUM];
  // Output QA
  TH1D*             hCent[2];
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

  //V0 Calib
  TH1D*             hPsiVZERODirectGet[NCENTBINS][3];
  TH1D*             hPsiV0Cor[NCENTBINS][3];
  TH1D*             hMultV0Fill[NRUNNUM];
  TH1D*             hMultV0GE[NRUNNUM];
  TH1D*             hMultV0Read[NRUNNUM];
  TProfile*          pV0XMeanFill[NRUNNUM][3];
  TProfile*          pV0YMeanFill[NRUNNUM][3];
  TH1D*             hQxnmV0[NRUNNUM][3];
  TH1D*             hQynmV0[NRUNNUM][3];
  TH2D*             hQnCentCor[3];
  TH2D*             hQxCentCor[3];
  TH2D*             hQyCentCor[3];
  TH2D*             hQxVtxCor[3];
  TH2D*             hQyVtxCor[3];

  // TPC Calib
  TH1D*             hPsiTPC[NCENTBINS][12];
  TH2D*             hQxTPCCent[9];
  TH2D*             hQyTPCCent[9];
  TH2D*             hQnTPCCent[9];
  TProfile*          pTPCCosMeanFill[NRUNNUM][4];
  TProfile*          pTPCSinMeanFill[NRUNNUM][4];
  TH1D*             hTPCCosMeanRead[NRUNNUM][4];
  TH1D*             hTPCSinMeanRead[NRUNNUM][4];
  TProfile2D*     pTPCShiftFillCoeffCos[3];
  TProfile2D*     pTPCShiftFillCoeffSin[3];
  TH2D*             hTPCShiftReadCoeffCos[3];
  TH2D*             hTPCShiftReadCoeffSin[3];
  AliAnalysisTaskEPCalib(const AliAnalysisTaskEPCalib&);
  AliAnalysisTaskEPCalib& operator=(const AliAnalysisTaskEPCalib&);

  ClassDef(AliAnalysisTaskEPCalib, 1);
};
#endif
