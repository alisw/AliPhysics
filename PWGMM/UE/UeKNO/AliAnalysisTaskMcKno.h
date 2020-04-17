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
	Float_t fdcaxy;
	Float_t fdcaz;	
	AliMultSelection *fMultSelection;

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

	TH2D * hRefMultvsV0Mmult;
	TH2D * hV0MmultvsUE;
	TH2D * hRefmultvsUE;

	TH2D * hPtVsUEGenTest[3];//UE->NchTS
	TH2D * hPtVsUERecTest[3];//UE->NchTS
	TH2D * hPtVsUEData[3];//UE->NchTS

	TH2D * hPtVsNchGenTest[3];
	TH2D * hPtVsNchRecTest[3];
	TH2D * hPtVsNchData[3];
	TH1D * hPhiGen[3];
	TH1D * hPhiRec[3];

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
