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
	void       GetMultiplicityDistributions();
	void       SetUseMC(Bool_t mc = kFALSE)              {fUseMC = mc;}   // use to analyse MC data
	virtual    Double_t DeltaPhi(Double_t phia, Double_t phib,
			Double_t rangeMin = -TMath::Pi()/2, Double_t rangeMax = 3*TMath::Pi()/2 );
protected:



private:
	AliESDEvent*            fESD;                                        //! input ESD event
	AliEventCuts        fEventCuts;
	AliStack*    fStack;                                                 //! MC stack
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

	TH1D * hPhiGen[3];
	TH1D * hNchTSGen;
	TH1D * hNchTSGenTest;
	TH1D * hPhiRec[3];
	TH1D * hNchTSRec;
        TH1D * hNchTSRecTest;
	TH2D * hNchResponse;

	AliAnalysisTaskMcKnoUe(const AliAnalysisTaskMcKnoUe&);                  // not implemented
	AliAnalysisTaskMcKnoUe& operator=(const AliAnalysisTaskMcKnoUe&);       // not implemented

	ClassDef(AliAnalysisTaskMcKnoUe, 3);
};

#endif
