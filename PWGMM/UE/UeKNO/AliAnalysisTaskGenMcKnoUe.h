/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. */
/* See cxx source for full Copyright notice */
/* Add a description of your MPI analysis */

#ifndef AliAnalysisTaskGenMcKnoUe_H
#define AliAnalysisTaskGenMcKnoUe_H

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



class AliAnalysisTaskGenMcKnoUe : public AliAnalysisTaskSE
{
	public:
		AliAnalysisTaskGenMcKnoUe();
		AliAnalysisTaskGenMcKnoUe(const char *name);

		virtual                 ~AliAnalysisTaskGenMcKnoUe();

		virtual void            UserCreateOutputObjects();
		virtual void            UserExec(Option_t* option);
		virtual void            Terminate(Option_t* option);

		Bool_t IsMCEventSelected(TObject* obj);

		void       GetGenLeadingObject();
		void       GetGenUEObservables();
		void       MakeALICE3Analysis();
		virtual    Double_t DeltaPhi(Double_t phia, Double_t phib,Double_t rangeMin = -TMath::Pi()/2.0, Double_t rangeMax = 3.0*TMath::Pi()/2.0 );
		void       SetPtMin(Double_t val){fPtMin = val;}
		void       SetIsPP(Bool_t val){fIsPP = val;}

	protected:

	private:

		AliMCEvent*  fMC;                     //! MC Event
		AliInputEventHandler*    fMcHandler;  //!<!
		AliStack*    fMCStack;                //! MC stack
		Double_t fEtaCut;
		Bool_t fIsPP;
		Double_t fPtMin;
		TList*  fOutputList;  //! output list in the root file

		Double_t fGenLeadPhi; 
		Double_t fGenLeadPt;
		Int_t    fGenLeadIn;
		// UE 
		TH1D * hPtLeadingTrue;
		TH2D * hPtLVsV0A;
		TH2D * hetaphi;
		TH2D * hnchmpi;
		TH3D * hnchmpirho;
		TH2D * hnchrho;
		TH2D * hmpirho;

		TH3D * hPtVsPtLeadingTrue[3];

		AliAnalysisTaskGenMcKnoUe(const AliAnalysisTaskGenMcKnoUe&);
		AliAnalysisTaskGenMcKnoUe& operator=(const AliAnalysisTaskGenMcKnoUe&);

		ClassDef(AliAnalysisTaskGenMcKnoUe, 3);
};

#endif
