/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. */
/* See cxx source for full Copyright notice */
/* Add a description of your MPI analysis */

#ifndef AliAnalysisTaskKnoUeChecks_H
#define AliAnalysisTaskKnoUeChecks_H

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



class AliAnalysisTaskKnoUeChecks : public AliAnalysisTaskSE
{
    public:
	AliAnalysisTaskKnoUeChecks();
	AliAnalysisTaskKnoUeChecks(const char *name);

	virtual                 ~AliAnalysisTaskKnoUeChecks();

	virtual void            UserCreateOutputObjects();
	virtual void            UserExec(Option_t* option);
	virtual void            Terminate(Option_t* option);

	//Getter
	void       GetLeadingObject(Bool_t isMC);
	void       GetTrackingEfficiencyTPConly();
	void       GetDetectorResponseDataDriven();
	void       GetMultiplicityDistributionsDataDriven();
	void       GetDetectorResponse();
	void       GetMultiplicityDistributions();
	void       GetMultiplicityDistributionsData();
	bool       HasRecVertex();
	virtual    Double_t DeltaPhi(Double_t phia, Double_t phib,Double_t rangeMin = -TMath::Pi()/2, Double_t rangeMax = 3*TMath::Pi()/2 );

	//Setter
	void       SetUseMC(Bool_t mc = kFALSE)             {fUseMC = mc;}              // use to analyse MC data
	void       SetMCclosureTest(Bool_t mcc = kFALSE)    {fIsMCclosure = mcc;}
	void       SetPtMin(Double_t val)                   {fPtMin = val;}             // use differnet ptcuts
	void       SetLeadingPtMin(Double_t PtLmin)         {fLeadPtCutMin = PtLmin;}   // use differnet ptcuts
	void       SetLeadingPtMax(Double_t PtLmax)         {fLeadPtCutMax = PtLmax;}   // use differnet ptcuts




    protected:


    private:
	//in header file the order of declaration shoukd be same as the initialization in .cxx file 
	AliESDEvent         *fESD;                  //! input ESD event
	AliMCEvent          *fMC;                   //! MC Event
	AliStack            *fMCStack;              //! MC stack
	Bool_t              fUseMC;                 // analyze MC events
	Bool_t              fIsMCclosure;
	AliEventCuts        fEventCuts;
	AliAnalysisFilter   *fLeadingTrackFilter;
	AliAnalysisFilter   *fTrackFilter;
	TList               *fOutputList;           //! output list in the root file

	Double_t  fEtaCut;
	Double_t  fPtMin;
	Double_t  fLeadPtCutMin;
	Double_t  fLeadPtCutMax;
	Double_t  fGenLeadPhi;
	Double_t  fGenLeadPt;
	Int_t     fGenLeadIn;
	Double_t  fRecLeadPhi;
	Double_t  fRecLeadPt;
	Int_t     fRecLeadIn;



	//histograms
	TH1I *hCounter;
	//vertex Z position
	TH1D *hZvtxAllMeasured;
	TH1D *hZvtxTrigMeasured;
	TH1D *hZvtxGoodVtxMeasured;
	TH1D *hZvtxCutAccMeasured;
	TH1D *hZvtxAllGen;
	TH1D *hZvtxCutGen;
	TH1D *hZvtxTrigGen;
	TH1D *hZvtxGoodVtxGen;
	TH1D *hZvtxCutAccGen;
	// KNO
	TH1D *hNchTSData;
	TH1D *hNchTSminData;
	TH1D *hNchTSmaxData;
	TH1D *hPhiData_TS1;
	TH1D *hPhiData_TS2;

	TH1D *hPhiGen[3];
	TH1D *hPhiRec[3];
	TH1D *hNchTSGen;
	TH1D *hNchTSRec;
	TH2D *hNchTSResponse;
	TH1D *hNchTSGenTest;
	TH1D *hNchTSRecTest;

	TH1D *hPhiGen_TS1;
	TH1D *hPhiGen_TS2;
	TH1D *hPhiRec_TS1;
	TH1D *hPhiRec_TS2;
	TH1D *hPhiGenTest_TS1;
	TH1D *hPhiGenTest_TS2;
	TH1D *hPhiRecTest_TS1;
	TH1D *hPhiRecTest_TS2;

	TH1D *hNchTSminGen;
	TH1D *hNchTSminRec;
	TH2D *hNchTSminResponse;
	TH1D *hNchTSminGenTest;
	TH1D *hNchTSminRecTest;

	TH1D *hNchTSmaxGen;
	TH1D *hNchTSmaxRec;
	TH2D *hNchTSmaxResponse;
	TH1D *hNchTSmaxGenTest;
	TH1D *hNchTSmaxRecTest;

	//data driven test
	TH1D *hPtPrimGen;
	TH1D *hPtPrimRec;
	TH1D *hNchTSDataTrue;
	TH1D *hNchTSDataMeasured;
	TH2D *hNchTSDataResponse;
	TH1D *hNchTSDataTrueTest;
	TH1D *hNchTSDataMeasuredTest;
	
	TH1D *hNchTSminDataTrue;
	TH1D *hNchTSminDataMeasured;
	TH2D *hNchTSminDataResponse;
	TH1D *hNchTSminDataTrueTest;
	TH1D *hNchTSminDataMeasuredTest;
	
	TH1D *hNchTSmaxDataTrue;
	TH1D *hNchTSmaxDataMeasured;
	TH2D *hNchTSmaxDataResponse;
	TH1D *hNchTSmaxDataTrueTest;
	TH1D *hNchTSmaxDataMeasuredTest;



	AliAnalysisTaskKnoUeChecks(const AliAnalysisTaskKnoUeChecks&);                  // not implemented
	AliAnalysisTaskKnoUeChecks& operator=(const AliAnalysisTaskKnoUeChecks&);       // not implemented

	ClassDef(AliAnalysisTaskKnoUeChecks, 3);
};

#endif
