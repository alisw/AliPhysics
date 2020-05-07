/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. */
/* See cxx source for full Copyright notice */
/* Add a description of your MPI analysis */

#ifndef AliAnalysisTaskMcKno_H
#define AliAnalysisTaskMcKno_H

class AliESDtrackCuts;
class AliESDEvent;
class TList;
class TH1D;
class TH2D;
class TH3D;
class TH1I;
class TProfile;
class THnSparse;

#include "AliAnalysisTaskSE.h"


#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliMCParticle.h"
#include "AliStack.h"
#include "TParticle.h"
#include "AliGenEventHeader.h"



class AliAnalysisTaskMcKno : public AliAnalysisTaskSE
{
public:
	AliAnalysisTaskMcKno();
	AliAnalysisTaskMcKno(const char *name);
	virtual                 ~AliAnalysisTaskMcKno();

	virtual void            UserCreateOutputObjects();
	virtual void            UserExec(Option_t* option);
	virtual void            Terminate(Option_t* option);
	void       GetLeadingObject(Bool_t isMC);
	void       GetDetectorResponse();
	void       GetBinByBinCorrections();
	void       GetMultiplicityDistributions();
	void       GetMultiplicityDistributionsData();
	void       GetMB();
	void       SetPtMin(Double_t val)              {fPtMin = val;}   // Set pT cut for associated particles
	void       SetLeadingPtMin(Double_t PtLmin)    {fLeadPtCutMin = PtLmin;}   // use differnet ptcuts
	void       SetLeadingPtMax(Double_t PtLmax)    {fLeadPtCutMax = PtLmax;}   // use differnet ptcuts
	void       SetV0Mmin(Double_t V0Mmin)          {fV0Mmin = V0Mmin;}   // Set V0M min value
	void       SetV0Mmax(Double_t V0Mmax)          {fV0Mmax = V0Mmax;}   // Set V0M max value
	void       SetUseMC(Bool_t mc = kFALSE)        {fUseMC = mc;}   // use to analyse MC data
	void       SetMCclosureTest(Bool_t mcc = kFALSE)    {fIsMCclosure = mcc;}
	void       SetIspPb(Bool_t pPb = kFALSE)    {fIspPb = pPb;}
	void       SetNchTScut(Bool_t TPConly = kTRUE)    {fIsTPConly = TPConly;}
	bool       HasRecVertex();
	//Systematic ============================
	void       SetTPCclustersVar1(Bool_t TPCclustersVar1 = kFALSE) {fTPCclustersVar1 = TPCclustersVar1;}
        void       SetTPCclustersVar2(Bool_t TPCclustersVar2 = kFALSE) {fTPCclustersVar2 = TPCclustersVar2;}
        void       SetNcrVar1(Bool_t NcrVar1 = kFALSE) {fNcrVar1 = NcrVar1;}
        void       SetNcrVar2(Bool_t NcrVar2 = kFALSE) {fNcrVar2 = NcrVar2;}
        void       SetChisqTPCVar1(Bool_t ChisqTPCVar1 = kFALSE) {fChisqTPCVar1 = ChisqTPCVar1;}
        void       SetChisqTPCVar2(Bool_t ChisqTPCVar2 = kFALSE) {fChisqTPCVar2 = ChisqTPCVar2;}
        void       SetChisqITSVar1(Bool_t ChisqITSVar1 = kFALSE) {fChisqITSVar1 = ChisqITSVar1;}
        void       SetChisqITSVar2(Bool_t ChisqITSVar2 = kFALSE) {fChisqITSVar2 = ChisqITSVar2;}
        void       SetChisqITSmTPCVar1(Bool_t ChisqITSmTPCVar1 = kFALSE) {fChisqITSmTPCVar1 = ChisqITSmTPCVar1;}
        void       SetChisqITSmTPCVar2(Bool_t ChisqITSmTPCVar2 = kFALSE) {fChisqITSmTPCVar2 = ChisqITSmTPCVar2;}
        void       SetDcazVar1(Bool_t DcazVar1 = kFALSE) {fDcazVar1 = DcazVar1;}
        void       SetDcazVar2(Bool_t DcazVar2 = kFALSE) {fDcazVar2 = DcazVar2;}
        void       SetGeoTPCVar1(Bool_t GeoTPCVar1 = kFALSE) {fGeoTPCVar1 = GeoTPCVar1;}
        void       SetGeoTPCVar2(Bool_t GeoTPCVar2 = kFALSE) {fGeoTPCVar2 = GeoTPCVar2;}
        void       SetGeoTPCVar3(Bool_t GeoTPCVar3 = kFALSE) {fGeoTPCVar3 = GeoTPCVar3;}
        void       SetGeoTPCVar4(Bool_t GeoTPCVar4 = kFALSE) {fGeoTPCVar4 = GeoTPCVar4;}
	void       SetSPDreqVar1(Bool_t SPDreqVar1 = kFALSE) {fSPDreqVar1 = SPDreqVar1;}
        //Systematic ============================
	virtual    Double_t DeltaPhi(Double_t phia, Double_t phib,
			Double_t rangeMin = -TMath::Pi()/2, Double_t rangeMax = 3*TMath::Pi()/2 );
protected:



private:
	AliESDEvent*            fESD;                                        //! input ESD event
	AliEventCuts        fEventCuts;
	AliStack*    fMCStack;                                                 //! MC stack
	AliMCEvent*  fMC;                                               //! MC Event
	Bool_t       fUseMC;                // analyze MC events
	Bool_t       fIsMCclosure;
	Bool_t       fIspPb;
	Bool_t       fIsTPConly;
	
	// Systematic------------------------------------
	Bool_t       fTPCclustersVar1;
	Bool_t       fTPCclustersVar2;
	Bool_t       fNcrVar1;
	Bool_t       fNcrVar2;
	Bool_t       fGeoTPCVar1;
	Bool_t       fGeoTPCVar2;
	Bool_t       fGeoTPCVar3;
	Bool_t       fGeoTPCVar4;
	Bool_t       fChisqTPCVar1;
	Bool_t       fChisqTPCVar2;
	Bool_t       fChisqITSVar1;
	Bool_t       fChisqITSVar2;
	Bool_t       fChisqITSmTPCVar1;
	Bool_t       fChisqITSmTPCVar2;
	Bool_t       fDcazVar1;
	Bool_t       fDcazVar2;
	Bool_t       fSPDreqVar1;
	// Systematic------------------------------------
	
	AliAnalysisFilter*  fLeadingTrackFilter;
	AliAnalysisFilter*  fTrackFilter;
	AliAnalysisFilter*  fTrackFilterwoDCA;
	TList*                  fOutputList;                                      //! output list in the root file
	

	Double_t fEtaCut;
	Double_t fPtMin;	
	Double_t fLeadPtCutMin;
	Double_t fLeadPtCutMax;
	Double_t fV0Mmin;
	Double_t fV0Mmax;
	Double_t fGenLeadPhi; 
	Double_t fGenLeadPt;
	Int_t    fGenLeadIn;
	Double_t fRecLeadPhi; 
	Double_t fRecLeadPt;
	Int_t    fRecLeadIn;
	Double_t ftrackmult08;
	Double_t fv0mpercentile;
	Double_t fv0mpercentilebefvtx;
	Float_t fdcaxy;
	Float_t fdcaz;	
	AliMultSelection *fMultSelection;
	AliMultSelection *fMultSelectionbefvtx;
       
	// KNO
	TH1D * hNchTSGen;
	TH1D * hNchTSGenTest;
	TH1D * hNchGen;
	TH1D * hNchGenTest;
	TH1D * hNchTSRec;
	TH1D * hNchTSRecTest;	
        TH1D * hNchData;
        TH1D * hNchTSData;
	TH2D * hNchResponse;
	TH1D * hNchRec;
	TH1D * hNchRecTest;
	TH1D * hPtInPrim;
	TH1D * hPtInPrim_pion;
	TH1D * hPtInPrim_kaon;
	TH1D * hPtInPrim_proton;
	TH1D * hPtInPrim_sigmap;
	TH1D * hPtInPrim_sigmam;
	TH1D * hPtInPrim_omega;
	TH1D * hPtInPrim_xi;
	TH1D * hPtInPrim_rest;
	TH1D * hPtOut;
	TH1D * hPtOutPrim;
	TH1D * hPtOutPrim_pion;
	TH1D * hPtOutPrim_kaon;
	TH1D * hPtOutPrim_proton;
	TH1D * hPtOutPrim_sigmap;
	TH1D * hPtOutPrim_sigmam;
	TH1D * hPtOutPrim_omega;
	TH1D * hPtOutPrim_xi;
	TH1D * hPtOutPrim_rest;
	TH1D * hPtOutSec; 
	TH1D * hCounter;
	TH1D * hRefMult08;
	TH1D * hV0Mmult;
	TH1D * hV0Mmultbefvtx;

	TH2D * hRefMultvsV0Mmult;
	TH2D * hV0MmultvsUE;
	TH2D * hRefmultvsUE;
	TH2D * hITSclustersvsUE;
	TH2D * hITSclustersvsNch;
	
	TH2D * hPtVsUEGenTest[3];//UE->NchTS
	TH2D * hPtVsUERecTest[3];//UE->NchTS
	TH2D * hPtVsUEData[3];//UE->NchTS

	TH2D * hPtVsNchGenTest[3];
	TH2D * hPtVsNchRecTest[3];
	TH2D * hPtVsNchData[3];
	TH1D * hPhiGen[3];
	TH1D * hPhiRec[3];

	TH2D * hPtVsV0MData;//V0M
	TH2D * hDphiVsUEGenTest; //UE->NchTS
	TH2D * hDphiVsUERecTest;//UE->NchTS
	TH2D * hDphiVsUEData;//UE->NchTS

	TH2D * hDphiVsNchGenTest;
	TH2D * hDphiVsNchRecTest;
	TH2D * hDphiVsNchData;

	//multiplicity percentile

	TH3D * hPtVsUEvsNchData_V0M[3];//UE->NchTS

	TH3D * hDphiVsUEvsNchData_V0M;//UE->NchTS

	//TH3D * hV0MVsUEvsRef;//UE->NchTS

	TH2D * hPTVsDCAData;
	TH2D * hPTVsDCAcentData;

	TH2F * hptvsdcaPrim;
	TH2F * hptvsdcaDecs;
	TH2F * hptvsdcaMatl;
	TH2F * hptvsdcacentralPrim;
	TH2F * hptvsdcacentralDecs;
	TH2F * hptvsdcacentralMatl;
	TH2F * hptvsdcaAll;
	TH2F * hptvsdcacentralAll;

	AliAnalysisTaskMcKno(const AliAnalysisTaskMcKno&);                  // not implemented
	AliAnalysisTaskMcKno& operator=(const AliAnalysisTaskMcKno&);       // not implemented

	ClassDef(AliAnalysisTaskMcKno, 3);
};

#endif
