/*
 * AliAnalysisTaskLeuteronNanoAOD.h
 *
 *  Created on:	06 November 2019
 *	Author:	Michael Jung
 */

#ifndef PWGCF_FEMTOSCOPY_FEMTODREAM_ALIANALYSISTASKLEUTERONNANOAOD_H_
#define PWGCF_FEMTOSCOPY_FEMTODREAM_ALIANALYSISTASKLEUTERONNANOAOD_H_

#include "AliAnalysisTaskSE.h"
#include "AliFemtoDreamEvent.h"
#include "AliFemtoDreamEventCuts.h"
#include "AliVTrack.h"
#include "AliFemtoDreamTrack.h"
#include "AliFemtoDreamTrackCuts.h"
#include "AliFemtoDreamv0.h"
#include "AliFemtoDreamv0Cuts.h"
#include "AliFemtoDreamCollConfig.h"
#include "AliFemtoDreamPairCleaner.h"
#include "AliFemtoDreamPartCollection.h"

class AliAnalysisTaskLeuteronNanoAOD : public AliAnalysisTaskSE {

  public:

    AliAnalysisTaskLeuteronNanoAOD();									      // class constructor without parameters
    AliAnalysisTaskLeuteronNanoAOD(const char* name,bool isMC,bool isHighMultV0, bool BruteForceDebugging);   // class constructor with parameters
    AliAnalysisTaskLeuteronNanoAOD& operator = (const AliAnalysisTaskLeuteronNanoAOD &task);		      // copy assignment operator
    AliAnalysisTaskLeuteronNanoAOD(const AliAnalysisTaskLeuteronNanoAOD &task);				      // copy constructor
    virtual ~AliAnalysisTaskLeuteronNanoAOD();								      // class destructor

    virtual void UserCreateOutputObjects();			// is called only once -> define output objects within this function
    virtual void UserExec(Option_t *option);			// is called in every event -> define what to search for in the events
    Float_t CalculateMassSqTOF(AliFemtoDreamTrack *track);	// calculate the mass^2 of the particle using TOF
    virtual void Terminate(Option_t *option){};			// is called only once -> terminates the analysis

    void SetEventCuts(AliFemtoDreamEventCuts *evtCuts){
      fEventCuts = evtCuts;
    };
    void SetTrackCutsPart1(AliFemtoDreamTrackCuts *trkCuts){	// particle 1 will be Protons
      fTrackCutsPart1 = trkCuts;
    };
    void SetTrackCutsPart2(AliFemtoDreamTrackCuts *trkCuts){	// particle 2 will be Antiprotons
      fTrackCutsPart2 = trkCuts;
    };
    void SetTrackCutsPart3(AliFemtoDreamTrackCuts *trkCuts){	// particle 3 will be Deuterons
      fTrackCutsPart3 = trkCuts;
    };
    void SetTrackCutsPart4(AliFemtoDreamTrackCuts *trkCuts){	// particle 4 will be Antideuterons
      fTrackCutsPart4 = trkCuts;
    };
    void Setv0CutsPart5(AliFemtoDreamv0Cuts *v0Cuts){		// particle 5 will be Lambdas
      fv0CutsPart5 = v0Cuts;
    };
    void Setv0CutsPart6(AliFemtoDreamv0Cuts *v0Cuts){		// particle 6 will be Antilambdas
      fv0CutsPart6 = v0Cuts;
    };
    void SetCollectionConfig(AliFemtoDreamCollConfig *config){
      fConfig = config;
    };

  private:

    void ResetGlobalTrackReference();				// is called in every event -> resets the track information
    void StoreGlobalTrackReference(AliVTrack *track);		// is called for every track -> stores the track information
    bool fIsMC;							// run over data "fIsMC(false)" or over Monte Carlo data "fIsMC(true)"
    bool fIsHighMultV0;
    bool fBruteForceDebugging;
    int fTrackBufferSize;						

    TList			    *fEventList;		// list for the event cuts
    TList			    *fProtonList;		// list for the proton cuts
    TList			    *fAntiprotonList;		// list for the antiproton cuts
    TList			    *fDeuteronList;		// list for the deuteron cuts
    TH2F			    *fDeuteronMassSqTOF;	// TH2F for calculation of deuteron mass2
    TList			    *fAntideuteronList;		// list for the antideuteron cuts
    TH2F			    *fAntideuteronMassSqTOF;	// TH2F for calculation of antideuteron mass2
    TList			    *fLambdaList;		// list for the lambda cuts
    TList			    *fAntilambdaList;		// list for the antilambda cuts
    TList			    *fPairCleanerList;		// list for the pair cleaner
    TList			    *fResultsList;		// list for the results
    TList			    *fResultsQAList;		// list for the QA of the results
    
    AliFemtoDreamEvent		    *fEvent;
    AliFemtoDreamTrack		    *fTrack;
    AliFemtoDreamv0		    *fFemtov0;

    AliFemtoDreamEventCuts	    *fEventCuts;		// cuts for the events
    AliFemtoDreamTrackCuts	    *fTrackCutsPart1;		// cuts for the tracks of particle 1 (Protons)
    AliFemtoDreamTrackCuts	    *fTrackCutsPart2;		// cuts for the tracks of particle 2 (Antiprotons)
    AliFemtoDreamTrackCuts	    *fTrackCutsPart3;		// cuts for the tracks of particle 3 (Deuterons)
    AliFemtoDreamTrackCuts	    *fTrackCutsPart4;		// cuts for the tracks of particle 4 (Antideuterons)
    AliFemtoDreamv0Cuts		    *fv0CutsPart5;		// cuts for the tracks of particle 5 (Lambdas)
    AliFemtoDreamv0Cuts		    *fv0CutsPart6;		// cuts for the tracks of particle 6 (Antilambdas)

    AliFemtoDreamCollConfig	    *fConfig;			// store the configurations needed for the calculation of the correlation function
    AliFemtoDreamPairCleaner	    *fPairCleaner;
    AliFemtoDreamPartCollection	    *fPartColl;
    AliVTrack			    **fGTI;			// global track information (GTI)

    ClassDef(AliAnalysisTaskLeuteronNanoAOD,1);
	
};

#endif /* PWGCF_FEMTOSCOPY_FEMTODREAM_ALIANALYSISTASKLEUTERONNANOAOD_H_ */

