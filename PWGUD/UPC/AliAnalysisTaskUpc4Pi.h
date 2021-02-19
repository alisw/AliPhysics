/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice */
/* $Id$ */

#ifndef ALIANALYSISTASKUPC4PI_H
#define ALIANALYSISTASKUPC4PI_H

class TObjString;
class TClonesArray;
class TFile;
class TTree;
class TList;
class TH1;
class TH2;
class TList;
class AliPIDResponse;
class AliESDEvent;
class TBits;
class AliTOFTriggerMask;

#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskUpc4Pi : public AliAnalysisTaskSE {
  public:
  	AliAnalysisTaskUpc4Pi();
  	AliAnalysisTaskUpc4Pi(const char *name, Bool_t);
	virtual ~AliAnalysisTaskUpc4Pi();

	virtual void Init();
	virtual void UserCreateOutputObjects();
	virtual void UserExec(Option_t *option);
	virtual void Terminate(Option_t *);
	
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
	static const int Maxtrk=8;
	Int_t ntrk;
  	TTree *f4PiTree;	//! Tree with all candidates
  	TTree *f4PiTree1;	//! Tree with 0 sum charge and only 4 tracks
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
	Float_t PIDTPCPion_T[8];
	Float_t PIDTPCElectron_T[8];
	Int_t TPCsignal_T[8];
	Float_t TrackP_T[8];
	Float_t Vertex_T[3];
	Int_t VtxContrib_T;
	Float_t VtxChi2_T,VtxNDF_T;
	Float_t SpdVertex_T[3];
	Int_t SpdVtxContrib_T;
	Int_t Ntracklets_T;
	Float_t Phi_T;
	Int_t NTracks_T;
	Float_t TrackEta_T[8];
	Float_t TrackPhi_T[8];
	Float_t TrackPx_T[8];
	Float_t TrackPy_T[8];
	Float_t TrackPz_T[8];
	Bool_t ChipCut_T;
	Bool_t TriggerTOF_T;
	Bool_t TriggerSPD_T;
	Int_t ITSModuleInner_T[8];
	Int_t ITSModuleOuter_T[8];

	// MC tree
	TTree *fMCTree;		//! Tree for MC (not used now)
	Int_t RunNum_MC_T;
	Float_t Mass_MC_T;
	Float_t Pt_MC_T;
	Float_t Rapidity_MC_T;
	Float_t Phi_MC_T;

	AliPIDResponse *fPIDResponse;
	TClonesArray *GenPart_T;

	TList *fListHist;	//! List of output Histograms
	TH1I *fHistTriggersPerRun;	//! Triggers hist
	TH1I *fITSmodule;		//! ITS hist
	TH1I *fFOchip;			//! FO chip hist
	TH1I *fFOcount;			//! FO count hist
	TH1F *TPCclustersP;		//! positive TPC clusters hist
	TH1F *TPCclustersN;		//! negative TPC clusters hist
	TH1F *fDeltaPhiRho;		//! Delta Phi Rho hist
	TH1F *fDeltaPhiEe;		//! Delta Phi Ee hist
	TH2F *dEdx;			//! dEdx hist
	TH2F *EtaPhiP;			//! EtaPhiP hist
	TH2F *EtaPhiN;			//! EtaPhiN hist

	// dEdx histograms
	TH2F *fHistdEdxVsP1;		//! dEdx hist1
	TH2F *fHistdEdxVsP2;		//! dEdx hist2
	TH2F *fHistdEdxVsP3;		//! dEdx hist3
	TH2F *fHistdEdxVsP4;		//! dEdx hist4
	TH2F *fHistdEdxVsP5;		//! dEdx hist5
	TH2F *fHistdEdxVsP6;		//! dEdx hist6
	TH2F *fHistdEdxVsP7;		//! dEdx hist7
	TH2F *fHistdEdxVsP8;		//! dEdx hist8
	TH2F *fHistdEdxVsP9;		//! dEdx hist9

	TH2F *fFOcorr;			//! FO corr hist
	TH1F *fGoodTracks;		//! Good Tracks number hist
	TH1F *fTrackChi2;		//! Track Chi2 hist

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

	AliAnalysisTaskUpc4Pi(const AliAnalysisTaskUpc4Pi&); //not implemented
	AliAnalysisTaskUpc4Pi& operator =(const AliAnalysisTaskUpc4Pi&); //not implemented
  
 	ClassDef(AliAnalysisTaskUpc4Pi, 3); 

};

#endif
