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
#include "AliHFCorrelator.h"

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
	

	void DefineHistoForAnalysis();
	
	// correlators
	void FillMCTagCorrelations(Double_t ptTrig, Double_t DelPhi,  Double_t DelEta, Double_t ptTrack, Bool_t *mcSource);	
	void FillMCTagLeadingCorrelations(Double_t ptTrig, Double_t DelPhi,  Double_t DelEta, Bool_t *mcSource);
	// checker for event mixing
	void EventMixingChecks(AliAODEvent * AOD); 
	// setters
	void SetMonteCarlo(Bool_t k) {fmontecarlo = k;}
	void SetUseMixing (Bool_t j) {fmixing = j;}
	void SetCorrelator(Int_t l) {fselect = l;} // select 1 for hadrons, 2 for Kaons, 3 for Kzeros
	void SetUseDisplacement(Int_t m) {fDisplacement=m;} // select 0 for no displ, 1 for abs displ, 2 for d0/sigma_d0
	void SetRunPbPb(Bool_t system){fSystem=system;} // select between pp (kFALSE) or PbPb (kTRUE)
    void SetLevelOfDebug(Int_t debug){fDebugLevel=debug;} // set debug level
	void SetUseReconstruction(Bool_t reco){fReco = reco;}

	


	private:

	AliAnalysisTaskDStarCorrelations(const AliAnalysisTaskDStarCorrelations &source);
	AliAnalysisTaskDStarCorrelations& operator=(const AliAnalysisTaskDStarCorrelations& source); 

	TObject* fhandler; //! Analysis Handler
	TClonesArray* fmcArray; //mcarray
	AliNormalizationCounter *fCounter; // counter
    AliHFCorrelator * fCorrelator; // object for correlations

	
	
	Int_t fselect; // select what to correlate with a D* 1-chargedtracks,2-chargedkaons,3-k0s
	Bool_t fmontecarlo;//switch for MC
	Bool_t fmixing;// switch for event mixing
	Bool_t fSystem; // pp or PbPb
	Bool_t fReco; // use reconstruction or MC truth
	
	Int_t fEvents; //! number of event
	Int_t fDebugLevel; //! debug level
	Int_t fDisplacement; // set 0 for no displacement cut, 1 for absolute d0, 2 for d0/sigma_d0
	
	
	TList *fOutput;                  //! user output data
	TList *fOutputMC;                //! outpu for MC
	AliRDHFCutsDStartoKpipi *fCuts;  // Cuts D*
	AliHFAssociatedTrackCuts *fCuts2; // cuts for associated 
					  
	ClassDef(AliAnalysisTaskDStarCorrelations,3); // class for D meson correlations				  
};

#endif
