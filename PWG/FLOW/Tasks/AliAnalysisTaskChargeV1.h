/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. */
/* See cxx source for full Copyright notice */
/* $Id$ */

#ifndef AliAnalysisTaskChargeV1_H
#define AliAnalysisTaskChargeV1_H

#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskChargeV1 : public AliAnalysisTaskSE
{
public:
  AliAnalysisTaskChargeV1();
  AliAnalysisTaskChargeV1(const char *name);
  virtual ~AliAnalysisTaskChargeV1();

  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *option);
  // Calib
  void SetListForNUE(TList *flist) { this->fListNUE = (TList *)flist->Clone(); }
  void SetListForNUA(TList *flist) { this->fListNUA = (TList *)flist->Clone(); }
  void SetNUEOn(bool doNUE) { this->IsDoNUE = doNUE; }
  void SetNUAOn(bool doNUA) { this->IsDoNUA = doNUA; }
  void SetListForZDCCalib(TList *flist) { this->fListZDCCalib = (TList *)flist->Clone(); }

private:
  bool GetVZEROPlane();
  double GetNUACor(int charge, double phi, double eta, double vz);
  double GetNUECor(int charge, double pt);
  bool GetZDCPlaneLsFit();
  inline double GetEventPlane(double qx, double qy, double harmonic);
  int GetRunNumBin(int runNum);
  bool LoadCalibHistForThisRun();
  bool RejectEvtTFFit();
  void ResetHists();
  bool QC2SE();
  bool CheckPIDofParticle(AliAODTrack *ftrack, int pidToCheck);

  AliAODEvent *fAOD;  // input event
  TList *fOutputList; // output list
  TList *fQAList;
  AliPIDResponse *fPID; // PID Handler
  double mHarmonic;
  unsigned int fFilterBit;
  double fPtMin;
  double fPtMax;
  double fEtaMax;
  double fNhitsMin;
  double fChi2Max;
  double fDeDxMin;
  double fNSigmaTPCCut;
  double fNSigmaTOFCut;
  TH3D **hYield;
  TProfile **pC2;
  TProfile **pC2QRe;
  TProfile **pC2QIm;
  TProfile3D **pD2pQStar;
  TProfile3D **pD2pRe;
  TProfile3D **pD2pIm;
  TProfile3D **v1p_qc;
  TProfile3D **v1t_qc;
  int runNum;
  int oldRunNum;
  int runNumBin;
  double vtx[3];
  int vzBin;
  double cent;
  double centSPD1;
  int centBin;
  // Event-wise
  TH1I *hEvtCount;
  TH1I *hRunNumBin;
  TH1D *hCent;
  TH2D *hCentCorr[3];
  TH2D *hVxy[2];
  TH1D *hVz[2];
  TH2D *fHist2DMultCentQA[2];
  // Track-wise
  TH1D *hPt[2];
  TH1D *hEta[2];
  TH2D *hBeforePhi[8];
  TH2D *hAfterPhi[8];
  TH1D *hNhits[2];
  TH2D *hPDedx;
  // pile up function
  TF1 *fSPDCutPU;
  TF1 *fV0CutPU;
  TF1 *fCenCutLowPU;
  TF1 *fCenCutHighPU;
  TF1 *fMultCutPU;

  ////////////////////////
  // NUE
  ////////////////////////
  bool IsDoNUE;
  TList *fListNUE; // read list for NUE
  TH1D *hNUEweightPlus;
  TH1D *hNUEweightMinus;

  ////////////////////////
  // NUA
  ////////////////////////
  bool IsDoNUA;
  TList *fListNUA;      // read lists for NUA
                        // 18q
  TH3F *hCorrectNUAPos; // Protty
  TH3F *hCorrectNUANeg; // Protty
  // TPC plane
  TH1D *pos1Plane;
  TH1D *neg1Plane;
  TProfile *Res1Square;
  TProfile2D **TPCcos_t;
  TProfile2D **TPCcos_p;

  // qc
  int nCentrality;
  TH2D *hMQ_thisEvt;
  TH2D *hReQ_thisEvt;
  TH2D *hImQ_thisEvt;
  TProfile *pC2_thisEvt;
  TProfile *pC2QRe_thisEvt;
  TProfile *pC2QIm_thisEvt;
  TProfile3D *pD2pQStar_thisEvt;
  TProfile3D *pD2pRe_thisEvt;
  TProfile3D *pD2pIm_thisEvt;
  TProfile3D *v1p_thisEvt;
  TProfile3D *v1t_thisEvt;
  TH3D *hMp_thisEvt;
  TH3D *hRep_thisEvt;
  TH3D *hImp_thisEvt;
  double PTGAP;
  double ETAGAP;
  int NPTBINS;
  double PTBINLOW;
  double PTBINUP;
  int NETABINS;
  double ETABINLOW;
  double ETABINUP;
  int NPOIBINS;

  // v1tp
  TProfile *ptEta;
  TProfile2D **px_P;
  TProfile2D **px_T;
  TProfile2D **v1_t;
  TProfile2D **v1_p;
  TProfile *ResQ;
  TH1D *Psi_P;
  TH1D *Psi_T;
  TH1D *Psi_PT;

  ////////////////////////
  // ZDC
  ////////////////////////
  TList *fListZDCCalib;
  // 18q 18r 15o
  TH1D *fHZDCCparameters;
  TH1D *fHZDCAparameters;
  TProfile *fProfileZDCPsi1Correlation; // ZNC-ZNA 1st
  TProfile *fProfileZDCPsi2Correlation; // ZNC-ZNA 2nd
  TH1D **fHist2Psi1ZNCCent;
  TH1D **fHist2Psi1ZNACent;

  //  QA
  //  ZDC
  //  ZNC

  TProfile *fProfileZNCTowerMeanEnegry[2]; // ![0]Raw [1]GE
  TProfile *fProfileZNCQxCent[2];          // ![0]GE  [1]RC
  TProfile *fProfileZNCQyCent[2];          // ![0]GE  [1]RC
  TH1D *fHist2CalibPsi1ZNCCent[3];         // ![0]GE  [1]RC [2]SF
  // ZNA
  TProfile *fProfileZNATowerMeanEnegry[2];
  TProfile *fProfileZNAQxCent[2];
  TProfile *fProfileZNAQyCent[2];
  TH1D *fHist2CalibPsi1ZNACent[3];
  // ZNC-ZNA
  TProfile *fProfileZDCQxAQxCCent[2];
  TProfile *fProfileZDCQxAQyCCent[2];
  TProfile *fProfileZDCQyAQxCCent[2];
  TProfile *fProfileZDCQyAQyCCent[2];


  double fPsi1ZNC;
  double fPsi1ZNA;

  //PID QA
  TH2D* fHistPIDPt;
  TH2D* fHistPIDEta;
  TH2D* fHistPIDPhi;
  TH2D* fHist2ProtonSigTPC;
  TH2D* fHist2ProtonSigTOF;
  TH2D* fHist2PionSigTPC;
  TH2D* fHist2PionSigTOF;
  TH2D* fHist2KionSigTPC;
  TH2D* fHist2KionSigTOF;


  // ZDC v1
  double Qtx;
  double Qty;
  double Qpx;
  double Qpy;
  TProfile2D **ZDCpx_P;
  TProfile2D **ZDCpx_T;
  TProfile2D **ZDCv1_t;
  TProfile2D **ZDCv1_p;
  TProfile   *ZDCResQ;
  TProfile2D **ZDCcos_t;
  TProfile2D **ZDCcos_p;

  AliAnalysisTaskChargeV1(const AliAnalysisTaskChargeV1 &);            // not implemented
  AliAnalysisTaskChargeV1 &operator=(const AliAnalysisTaskChargeV1 &); // not implemented

  ClassDef(AliAnalysisTaskChargeV1, 1);
};

#endif
