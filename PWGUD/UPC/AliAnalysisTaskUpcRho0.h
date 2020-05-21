/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice */
/* $Id$ */

#ifndef ALIANALYSISTASKUPCRHO0_H
#define ALIANALYSISTASKUPCRHO0_H

class TClonesArray;
class TFile;
class TTree;
class TList;
class TH1;
class TList;
class AliPIDResponse;
class AliESDEvent;
class TBits;
class AliTOFTriggerMask;

#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskUpcRho0 : public AliAnalysisTaskSE {
  public:
  	AliAnalysisTaskUpcRho0();
  	AliAnalysisTaskUpcRho0(const char *name, Bool_t);
	virtual ~AliAnalysisTaskUpcRho0();

	virtual void Init();
	virtual void UserCreateOutputObjects();
	virtual void UserExec(Option_t *option);
	virtual void Terminate(Option_t *){};

	void SetIsMC(Bool_t _isMC){ isMC = _isMC; }
	void SetDebugMode(Bool_t _debugMode){ debugMode = _debugMode; }
	void SetEfficiencyFileName(TString _fEfficiencyFileName){ fEfficiencyFileName = _fEfficiencyFileName; isUsingEffi = kTRUE; }
	void SetTOFFileName(TString _name) {fTOFFileName = _name; isUsingTOFeff = kTRUE;}
 	void SetTrigger(TString _fTriggerName){ fTriggerName = _fTriggerName; }
 	void SetTPCNcls(Int_t _fTPCNcls) {fTPCNcls = _fTPCNcls;}
 	void SetOption(TString _fOption){fOption = _fOption;}

  private:
  	Bool_t Is0STPfired(Int_t *, Int_t *);
  	Bool_t IsTriggered(AliESDEvent *);

  	Bool_t isMC;
  	Bool_t debugMode;
  	TString fTriggerName;
  	Int_t fTPCNcls; // number of TPC clusters
  	TString fOption;

  	// tree
  	TTree *fRhoTree;
	// tree variables and branches
	Int_t RunNum_T;
	UShort_t BunchCrossNumber_T;
	UInt_t OrbitNumber_T;
	UInt_t PeriodNumber_T;
	Bool_t LikeSign_T;
	Float_t Mass_T;
	Float_t Pt_T;
	Float_t Rapidity_T;
	Int_t V0Adecision_T;
	Int_t V0Cdecision_T;
	Int_t ADAdecision_T;
	Int_t ADCdecision_T;
	Bool_t UBAfired_T;
	Bool_t UBCfired_T;
	Bool_t VBAfired_T;
	Bool_t VBCfired_T;
	Float_t ZNAenergy_T;
	Float_t ZNCenergy_T;
	Float_t ZPAenergy_T;
	Float_t ZPCenergy_T;
	Float_t ZDCAtime_T[4];
	Float_t ZDCCtime_T[4];
	Float_t PIDTPCPion_T[2];
	Float_t PIDTPCElectron_T[2];
	Int_t TPCsignal_T[2];
	Float_t TrackP_T[2];
	Float_t Vertex_T[3];
	Int_t VtxContrib_T;
	Float_t VtxChi2_T,VtxNDF_T;
	Float_t SpdVertex_T[3];
	Int_t SpdVtxContrib_T;
	Int_t Ntracklets_T;
	Float_t Phi_T;
	Float_t TrackEta_T[2];
	Float_t TrackPhi_T[2];
	Float_t TrackPx_T[2];
	Float_t TrackPy_T[2];
	Float_t TrackPz_T[2];
	Bool_t ChipCut_T;
	Int_t ITSModuleInner_T[2];
	Int_t ITSModuleOuter_T[2];

	// MC tree
	TTree *fMCTree;
	Int_t RunNum_MC_T;
	Float_t Mass_MC_T;
	Float_t Pt_MC_T;
	Float_t Rapidity_MC_T;
	Float_t Phi_MC_T;

	AliPIDResponse *fPIDResponse;
	TClonesArray *GenPart_T;

	TList *fListHist;
	TH1I *fHistTriggersPerRun;
	TH1I *fITSmodule;
	TH1I *fFOchip;
	TH1I *fFOcount;
	TH1F *TPCclustersP;
	TH1F *TPCclustersN;
	TH1F *fDeltaPhiRho;
	TH1F *fDeltaPhiEe;
	TH2F *dEdx;
	TH2F *EtaPhiP;
	TH2F *EtaPhiN;

	// dEdx histograms
	TH2F *fHistdEdxVsP1;
	TH2F *fHistdEdxVsP2;
	TH2F *fHistdEdxVsP3;
	TH2F *fHistdEdxVsP4;
	TH2F *fHistdEdxVsP5;
	TH2F *fHistdEdxVsP6;
	TH2F *fHistdEdxVsP7;
	TH2F *fHistdEdxVsP8;
	TH2F *fHistdEdxVsP9;

	TH2F *fFOcorr;
	TH1F *fGoodTracks;
	TH1F *fTrackChi2;

	// SPD effi
	Bool_t isUsingEffi;
	TString fEfficiencyFileName;
	TFile *fSPDfile;
	TH1D *hBCmod4;
	TH2D *hSPDeff;

	// TOF effi
	Bool_t isUsingTOFeff;
	TFile *fTOFfile;
	TString fTOFFileName;
	TH2F *hTOFeff;
	AliTOFTriggerMask *fTOFmask;

	AliAnalysisTaskUpcRho0(const AliAnalysisTaskUpcRho0&); //not implemented
	AliAnalysisTaskUpcRho0& operator =(const AliAnalysisTaskUpcRho0&); //not implemented
  
 	ClassDef(AliAnalysisTaskUpcRho0, 7); 

};

#endif