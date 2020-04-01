/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. */
/* See cxx source for full Copyright notice */
/* Add a description of your MPI analysis */

#ifndef AliAnalysisTaskMcKnoUe_H
#define AliAnalysisTaskMcKnoUe_H

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



class AliAnalysisTaskMcKnoUe : public AliAnalysisTaskSE
{
public:
	AliAnalysisTaskMcKnoUe();
	AliAnalysisTaskMcKnoUe(const char *name);
	virtual                 ~AliAnalysisTaskMcKnoUe();

	virtual void            UserCreateOutputObjects();
	virtual void            UserExec(Option_t* option);
	virtual void            Terminate(Option_t* option);
	void       GetLeadingObject(Bool_t isMC);
	void       GetDetectorResponse();
	void       GetBinByBinCorrections();
	void       GetUEObservables();
        void       GetUEObservablesData();
	void       GetPtLeadingMisRecCorrection();
	void       GetMeanUEObservables(std::vector<Double_t> &gen, std::vector<Double_t> &rec);
	void       GetMultiplicityDistributions();
	void       GetMultiplicityDistributionsData();
	void       SetPtMin(Double_t val)              {fPtMin = val;}   // use differnet ptcuts
	void       SetUseMC(Bool_t mc = kFALSE)              {fUseMC = mc;}   // use to analyse MC data
	void       SetMCclosureTest(Bool_t mcc = kFALSE)    {fIsMCclosure = mcc;}
	void       SetParametrizationEfficiency(Bool_t ispy = kTRUE)  {fIsPythia = ispy;}
	bool       HasRecVertex();
	virtual    Double_t DeltaPhi(Double_t phia, Double_t phib,
			Double_t rangeMin = -TMath::Pi()/2, Double_t rangeMax = 3*TMath::Pi()/2 );
protected:



private:
	AliESDEvent*            fESD;                                        //! input ESD event
	Bool_t       fIsPythia; 
	AliEventCuts        fEventCuts;
	AliStack*    fMCStack;                                                 //! MC stack
	AliMCEvent*  fMC;                                               //! MC Event
	Bool_t       fUseMC;                // analyze MC events
	Bool_t       fIsMCclosure;         
	AliAnalysisFilter*  fLeadingTrackFilter;
	AliAnalysisFilter*  fTrackFilter;
	TList*                  fOutputList;                                      //! output list in the root file

	Double_t fEtaCut;
	Double_t fPtMin;
	Double_t fLeadPtCutMin;
	Double_t fLeadPtCutMax;
	Double_t fGenLeadPhi; 
	Double_t fGenLeadPt;
	Int_t    fGenLeadIn;
	Double_t fRecLeadPhi; 
	Double_t fRecLeadPt;
	Int_t    fRecLeadIn;

	// KNO
	TH1D * hPhiGen[3];
	TH1D * hNchTSGen;
	TH1D * hNchTSGenTest;
	TH1D * hPhiRec[3];
	TH1D * hNchTSRec;
	TH1D * hNchTSRecTest;
	TH1D * hNchTSData;
	TH2D * hNchResponse;

	// UE 
	TH1D * hPtInPrim;
	TH1D * hPtOut;
	TH1D * hPtOutPrim; 
	TH1D * hPtOutSec; 
	TH1D * hCounter;
	TH2D * hNumDenMC[3];
	TH2D * hSumPtMC[3];
	TH2D * hNumDenMCMatch[3];
	TH2D * hSumPtMCMatch[3];
	TH2D * hNumDenMCDd[3];
	TH2D * hSumPtMCDd[3];
	TH2D * hNumDenMCMatchDd[3];
	TH2D * hSumPtMCMatchDd[3];

	TH1D * hPtLeadingTrue;
	TH1D * hPtLeadingMeasured;
	TH1D * hPtLeadingData;
	TH2D * hPtVsPtLeadingMeasured[3];
	TH2D * hPtVsPtLeadingData[3];
	TH2D * hPtVsPtLeadingTrue[3];
	TProfile * pNumDenMeasured[3];
	TProfile * pNumDenData[3];
	TProfile * pNumDenTrue[3];
	TProfile * pSumPtMeasured[3];
	TProfile * pSumPtData[3];
	TProfile * pSumPtTrue[3];

	TProfile * pNumDenMeasuredAll[3];
	TProfile * pNumDenTrueAll[3];
	TProfile * pSumPtMeasuredAll[3];
	TProfile * pSumPtTrueAll[3];

	TProfile * pNumDenMeasuredPS[3];
	TProfile * pNumDenTruePS[3];
	TProfile * pSumPtMeasuredPS[3];
	TProfile * pSumPtTruePS[3];

	TProfile * pNumDenMeasuredPSV[3];
	TProfile * pNumDenTruePSV[3];
	TProfile * pSumPtMeasuredPSV[3];
	TProfile * pSumPtTruePSV[3];

	TProfile * pNumDenMeasuredGood[3];
	TProfile * pNumDenTrueGood[3];
	TProfile * pSumPtMeasuredGood[3];
	TProfile * pSumPtTrueGood[3];

	TH2D * hPtVsUEGenTest[3];
	TH2D * hPtVsUERecTest[3];
	TH2D * hPtVsUEData[3];

	TH1D * hPtLeadingRecPS;
	TH1D * hPtLeadingRecPSV;
	TH1D * hPtLeadingRecGood;
	TH1D * hPtLeadingGenPS;
	TH1D * hPtLeadingGenPSV;
	TH1D * hPtLeadingGenGood;
	TH1D * hPtLeadingRecAll;
	TH1D * hPtLeadingGenAll;

	AliAnalysisTaskMcKnoUe(const AliAnalysisTaskMcKnoUe&);                  // not implemented
	AliAnalysisTaskMcKnoUe& operator=(const AliAnalysisTaskMcKnoUe&);       // not implemented

	ClassDef(AliAnalysisTaskMcKnoUe, 3);
};

#endif
