#ifndef AliAnalysisTaskDStarCorrelations_H
#define AliAnalysisTaskDStarCorrelations_H

/**************************************************************************
 * Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id$ */

//-----------------------------------------------------------------------
//          
//
//						   Author S.Bjelogrlic
//                         Utrecht University 
//                      sandro.bjelogrlic@cern.ch
//
//-----------------------------------------------------------------------

#include <TH2F.h>
#include <TH1D.h>
#include <TH3D.h>
#include "AliAnalysisTaskSE.h"
#include "AliAODEvent.h"
#include "AliRDHFCutsDStartoKpipi.h"
#include "AliEventPoolManager.h"
#include "AliAODRecoCascadeHF.h"
#include "AliHFAssociatedTrackCuts.h"
#include "AliNormalizationCounter.h"

class TParticle ;
class TClonesArray ;
class AliAODMCParticle;
class AliAODEvent;
class AliVParticle;
class TObjArray;
class AliEventPoolManager;
class AliESDEvent;



class AliAnalysisTaskDStarCorrelations : public AliAnalysisTaskSE
{

	public :
	AliAnalysisTaskDStarCorrelations();
	AliAnalysisTaskDStarCorrelations(const Char_t* name,AliRDHFCutsDStartoKpipi* cuts, AliHFAssociatedTrackCuts *AsscCuts);
	virtual ~AliAnalysisTaskDStarCorrelations();
	
	
	
	// Implementation of interface methods  
	virtual void UserCreateOutputObjects();
	virtual void Init();
	virtual void LocalInit() {Init();}
	virtual void UserExec(Option_t *option);
	virtual void Terminate(Option_t *option);
	
	
	// methods to get the tracks to correlate
	TObjArray*  AcceptAndReduceTracks(AliAODEvent* inputEvent);
	TObjArray*  AcceptAndReduceKZero(AliAODEvent* inputEvent, Int_t loopindex, Int_t plotassociation);

	void DefineHistoForAnalysis();
	
	// correlators
	void FillCorrelations(Double_t ptTrig, Double_t phiTrig, Double_t etaTrig, Double_t phiTrack, Double_t etaTrack);
	void FillSideBandCorrelations(Double_t ptTrig, Double_t phiTrig, Double_t etaTrig, Double_t phiTrack, Double_t etaTrack);
	void FillMCTagCorrelations(Double_t ptTrig, Double_t phiTrig,  Double_t etaTrig, Double_t ptTrack, Double_t phiTrack, Double_t etaTrack, Int_t mcSource);
	
	// setters
	void SetMonteCarlo(Bool_t k) {fmontecarlo = k;}
	void SetUseMixing (Bool_t j) {fmixing = j;}
	void SetCorrelator(Int_t l) {fselect = l;}
	


	private:

	AliAnalysisTaskDStarCorrelations(const AliAnalysisTaskDStarCorrelations &source);
	AliAnalysisTaskDStarCorrelations& operator=(const AliAnalysisTaskDStarCorrelations& source); 

	TObject* fhandler; //! Analysis Handler
	AliEventPoolManager*     fPoolMgr;         //! event pool manager
	TClonesArray* fmcArray; //mcarray
	AliNormalizationCounter *fCounter; // counter
	
	
	Int_t fselect; // select what to correlate with a D* 1-chargedtracks,2-chargedkaons,3-k0s
	Bool_t fmontecarlo;//switch for MC
	Bool_t fmixing;// switch for event mixing
	Int_t fEvents; //! number of event
	Int_t fDebug; //! debug level
	
	
	TList *fOutput;                  //! user output data
	AliRDHFCutsDStartoKpipi *fCuts;  // Cuts D*
	AliHFAssociatedTrackCuts *fCuts2; // cuts for associated 
					  
	ClassDef(AliAnalysisTaskDStarCorrelations,1); // class for D meson correlations				  
};

#endif
