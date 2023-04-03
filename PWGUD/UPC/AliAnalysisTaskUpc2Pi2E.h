/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice */
/* $Id$ */

#ifndef ALIANALYSISTASKUPC2PI2E_H
#define ALIANALYSISTASKUPC2PI2E_H

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
class TBits;
class AliESDtrackCuts;

#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskUpc2Pi2E : public AliAnalysisTaskSE {
  public:
  	AliAnalysisTaskUpc2Pi2E();
  	AliAnalysisTaskUpc2Pi2E(const char *name, Bool_t);
	virtual ~AliAnalysisTaskUpc2Pi2E();

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
        void    SetCrossed(Int_t spd[4], TBits &crossed);
        Int_t   GetChipId(Int_t index, Int_t &chipId2, Bool_t debug = 0);
        Bool_t  IsSTGFired(TBits bits, Int_t dphiMin = 4, Int_t dphiMax = 10, Bool_t tolerance = 1);
  private:
  	Bool_t Is0STPfired(Int_t *, Int_t *);
  	Bool_t IsTriggered(AliESDEvent *);
	AliESDtrackCuts *fTrackCutsBit4;

  	Bool_t isMC;
  	Bool_t debugMode;
  	TString fTriggerName;
  	Int_t fTPCNcls; // number of TPC clusters
  	TString fOption;

  	// tree
	static const int Maxtrk=20;
	Int_t ntrk;
  	TTree *f2Pi2ETree;	//! Tree with all candidates
  	TTree *f2Pi2ETree1;	//! Tree with 0 sum charge and only 4 tracks
	// tree variables and branches
	std::vector<int> FORChip;
	UShort_t BunchCrossNumber_T;
	UInt_t OrbitNumber_T;
	UInt_t PeriodNumber_T;
	Int_t RunNum_T;
	Bool_t LikeSign_T;
	Float_t Mass_T;
	Float_t Pt_T;
	Int_t Trigger_T;
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
	Float_t PIDTPCPion_T[20];
	Float_t PIDTPCMuon_T[20];
	Float_t PIDTPCKaon_T[20];
	Float_t PIDTPCProton_T[20];
	Float_t PIDTPCElectron_T[20];
	Int_t TPCsignal_T[20];
	Bool_t TPCrefit_T[20];
	Bool_t ITSrefit_T[20];
	Bool_t ITSO_T[20];
	Bool_t ITSI_T[20];
	Bool_t ITSSA_T[20];
	Bool_t TrackCuts_T[20];
	Float_t TrackDCAxy_T[20];
	Float_t TrackDCAz_T[20];
	Bool_t fMatchingSPD_T[20];
	Int_t TPCNcls_T[20];
	Int_t ITSNcls_T[20];
	Float_t TPCchi2_T[20];
	Float_t TrackP_T[20];
	Int_t TrackC_T[20];
	Float_t Vertex_T[3];
	Int_t VtxContrib_T;
	Float_t VtxChi2_T,VtxNDF_T;
	Float_t SpdVertex_T[3];
	Int_t SpdVtxContrib_T;
	Int_t Ntracklets_T;
	Float_t Phi_T;
	Int_t NTracks_T;
	Float_t TrackEta_T[20];
	Float_t TrackPhi_T[20];
	Float_t TrackPx_T[20];
	Float_t TrackPy_T[20];
	Float_t TrackPz_T[20];
	Bool_t ChipCut_T;
	Bool_t TriggerTOF_T;
	Bool_t TriggerSPD_T;
	Int_t ITSModuleInner_T[20];
	Int_t ITSModuleOuter_T[20];
	Int_t fFOmodules_T[240];

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
        TBits fFOCrossFiredChips;

	AliAnalysisTaskUpc2Pi2E(const AliAnalysisTaskUpc2Pi2E&); //not implemented
	AliAnalysisTaskUpc2Pi2E& operator =(const AliAnalysisTaskUpc2Pi2E&); //not implemented
  
 	ClassDef(AliAnalysisTaskUpc2Pi2E, 6); 

};

#endif
