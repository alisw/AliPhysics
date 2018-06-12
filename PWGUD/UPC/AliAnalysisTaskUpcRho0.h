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
// class TH2;
class TList;
class AliPIDResponse;
// class AliAODEvent;
class AliESDEvent;
// class AliTOFTriggerMask;
class TBits;

#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskUpcRho0 : public AliAnalysisTaskSE {
  public:
  	AliAnalysisTaskUpcRho0();
  	AliAnalysisTaskUpcRho0(const char *name, Bool_t);
	virtual ~AliAnalysisTaskUpcRho0();

	virtual void Init();
	virtual void UserCreateOutputObjects();
	virtual void UserExec(Option_t *option);
	// virtual void RunAODtrig();
	// virtual void RunAODhist();
	// virtual void RunAODtree();
	// virtual void RunAODMC(AliAODEvent *aod);
	// virtual void RunAODsystematics(AliAODEvent *aod);
	// virtual void RunESDtrig();
	// virtual void RunESDhist();
	// virtual void RunESDtree();
	// virtual void RunESDMC(AliESDEvent *esd);
	virtual void Terminate(Option_t *){};

	void SetIsMC(Bool_t _isMC){ isMC = _isMC; }
	void SetEfficiencyFile(TFile *_fSPDfile){ fSPDfile = _fSPDfile; };
 
  private:
  	Bool_t isMC;
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
	Float_t SpdVertex_T[3];
	Float_t DeltaPhi_T;
	Int_t Ntracklets_T;
	Float_t Phi_T;
	Float_t TrackEta_T[2];
	Float_t TrackPhi_T[2];
	Float_t TrackPx_T[2];
	Float_t TrackPy_T[2];
	Float_t TrackPz_T[2];
	Bool_t ChipCut_T;
	Int_t ITSModule_T;

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
	TH2F *dEdx;
	TH2F *EtaPhiP;
	TH2F *EtaPhiN;

	TH2F *fFOcorr;

	// SPD effi
	TFile *fSPDfile;
	TH1D *hBCmod4;
	TH2D *hSPDeff;

	AliAnalysisTaskUpcRho0(const AliAnalysisTaskUpcRho0&); //not implemented
	AliAnalysisTaskUpcRho0& operator =(const AliAnalysisTaskUpcRho0&); //not implemented
  
 	ClassDef(AliAnalysisTaskUpcRho0, 7); 

};

#endif