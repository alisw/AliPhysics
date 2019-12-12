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
	void       GetPtLeadingMisRecCorrection();
	void       GetMultiplicityDistributions();
        void       SetPtMin(Double_t val)              {fPtMin = val;}   // use differnet ptcuts
	void       SetUseMC(Bool_t mc = kFALSE)              {fUseMC = mc;}   // use to analyse MC data
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

	TH1D * hPtLeadingTrue;
	TH1D * hPtLeadingMeasured;
	TH2D * hPtVsPtLeadingMeasured[3];
	TProfile * pNumDenMeasured[3];
	TProfile * pNumDenTrue[3];
	TProfile * pSumPtMeasured[3];
	TProfile * pSumPtTrue[3];

	TH2D * hPtVsUEGenTest[3];
	TH2D * hPtVsUERecTest[3];

	TH1D * hPtLeadingRecPS;
	TH1D * hPtLeadingRecPSV;
	TH1D * hPtLeadingGenPS;
	TH1D * hPtLeadingGenPSV;

	AliAnalysisTaskMcKnoUe(const AliAnalysisTaskMcKnoUe&);                  // not implemented
	AliAnalysisTaskMcKnoUe& operator=(const AliAnalysisTaskMcKnoUe&);       // not implemented

	ClassDef(AliAnalysisTaskMcKnoUe, 3);
};

#endif
