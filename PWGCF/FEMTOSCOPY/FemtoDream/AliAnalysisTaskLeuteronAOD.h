/*
 * AliAnalysisTaskLeuteronAOD.h
 *
 *  Created on:	19 December 2019
 *	Author:	Michael Jung
 */

#ifndef PWGCF_FEMTOSCOPY_FEMTODREAM_ALIANALYSISTASKLEUTERONAOD_H_
#define PWGCF_FEMTOSCOPY_FEMTODREAM_ALIANALYSISTASKLEUTERONAOD_H_

#include "Rtypes.h"

#include "AliAnalysisTaskSE.h"
#include "AliFemtoDreamEvent.h"
#include "AliFemtoDreamEventCuts.h"
#include "AliAODTrack.h"
#include "AliFemtoDreamTrack.h"
#include "AliFemtoDreamTrackCuts.h"
#include "AliFemtoDreamv0.h"
#include "AliFemtoDreamv0Cuts.h"
#include "AliFemtoDreamCollConfig.h"
#include "AliFemtoDreamPairCleaner.h"
#include "AliFemtoDreamPartCollection.h"
#include "TF1.h"

class AliAnalysisTaskLeuteronAOD : public AliAnalysisTaskSE {

  public:

    AliAnalysisTaskLeuteronAOD();									// class constructor without parameters
    AliAnalysisTaskLeuteronAOD(const char* name,bool isMC,bool isHighMultV0,bool BruteForceDebugging,bool isSidebandSignal, bool isUpperSideband, bool isLowerSideband,bool isPbPb,bool doEventQAPlots, bool doResultsQAPlots,bool isCentral);	// class constructor with parameters
    AliAnalysisTaskLeuteronAOD& operator = (const AliAnalysisTaskLeuteronAOD &task);			// copy assignment operator
    AliAnalysisTaskLeuteronAOD(const AliAnalysisTaskLeuteronAOD &task);					// copy constructor
    virtual ~AliAnalysisTaskLeuteronAOD();								// class destructor

    virtual void UserCreateOutputObjects();	      // is called only once -> define output objects within this function
    virtual void UserExec(Option_t *option);	      // is called in every event -> define what to search for in the events 
    Float_t CalculateMassSqTOF(AliFemtoDreamTrack *track);      // calculate the mass^2 of the particle using TOF
    bool CheckTPCDeuteronPID(AliFemtoDreamTrack *track, double nSigma);
    bool CheckDeuteronMassSquarePID(AliFemtoDreamTrack *track,double nSigma);
    bool CheckAntiDeuteronMassSquarePID(AliFemtoDreamTrack *track,double nSigma);
    Double_t GetDeuteronMass2Mean_pp(float pT);
    Double_t GetLimit(float pT, double mean, double sign,double offset,double lastpar);
    Double_t GetAntideuteronMass2Mean_pp(float pT);
    virtual void Terminate(Option_t *option){};	      // is called only once -> terminates the analysis

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
    void SetTrackCutsPart3Mass(AliFemtoDreamTrackCuts *trkCuts){	// particle 3 will be DeuteronsMass
      fTrackCutsPart3Mass = trkCuts;
    };
    void SetTrackCutsPart3Sigma(AliFemtoDreamTrackCuts *trkCuts){	// particle 3 will be DeuteronsSigma
      fTrackCutsPart3Sigma = trkCuts;
    };
    void SetTrackCutsPart4(AliFemtoDreamTrackCuts *trkCuts){	// particle 4 will be Antideuterons
      fTrackCutsPart4 = trkCuts;
    };
    void SetTrackCutsPart4Mass(AliFemtoDreamTrackCuts *trkCuts){	// particle 4 will be AntideuteronsMass
      fTrackCutsPart4Mass = trkCuts;
    };
    void SetTrackCutsPart4Sigma(AliFemtoDreamTrackCuts *trkCuts){	// particle 4 will be AntideuteronsSigma
      fTrackCutsPart4Sigma = trkCuts;
    };
    void Setv0CutsPart5(AliFemtoDreamv0Cuts *v0Cuts){		// particle 3 will be Lambdas
      fv0CutsPart5 = v0Cuts;
    };
    void Setv0CutsPart6(AliFemtoDreamv0Cuts *v0Cuts){		// particle 4 will be Antilambdas
      fv0CutsPart6 = v0Cuts;
    };
    void SetCollectionConfig(AliFemtoDreamCollConfig *config){
      fConfig = config;
    };

  private:

    void ResetGlobalTrackReference();				// is called in every event -> resets the track information
    void StoreGlobalTrackReference(AliAODTrack *track);		// is called for every track -> stores the track information
    bool fIsMC;							// run over data "fIsMC(false)" or over Monte Carlo data "fIsMC(true)"
    bool fIsHighMultV0;
    bool fBruteForceDebugging;
    bool fisSidebandSignal;
    bool fisUpperSideband;
    bool fisLowerSideband;
    bool fisPbPb;
    bool fisCentral;
    int fTrackBufferSize;			  

    TList			    *fEventList;		// list for the event cuts
    TList			    *fProtonList;		// list for the proton cuts
    TList			    *fAntiprotonList;		// list for the antiproton cuts
    TList			    *fDeuteronList;		// list for the deuteron cuts
    TH2F                            *fDeuteronMassSqTOF;        // TH2F for calculation of deuteron mass2
    TH2F                            *fDeuteronMassSqTOFFullPt;  // TH2F for calculation of deuteron mass2 full pt
    TH2F                            *fDeuteronTPCnSigma;	// TH2F for calculation of nsigma TPC for purity
    TList			    *fAntideuteronList;		// list for the antideuteron cuts
    TH2F                            *fAntideuteronMassSqTOF;	// TH2F for calculation of antideuteron mass2
    TH2F                            *fAntideuteronMassSqTOFFullPt;  // TH2F for calculation of antideuteron mass2 full pt
    TH2F                            *fAntideuteronTPCnSigma;	// TH2F for calculation of nsigma TPC for purity
    TList			    *fLambdaList;		// list for the lambda cuts
    TList			    *fAntilambdaList;		// list for the antilambda cuts
    TList			    *fPairCleanerList;		// list for the pair cleaner
    TList			    *fResultsList;		// list for the results
    TList			    *fResultsQAList;		// list for the QA of the results
    TH1F			    *fSimpleEventCounter;	// count the number of events
    TH1F			    *fEventCentrality;		// centrality of PbPb events

    AliFemtoDreamEvent		    *fEvent;
    AliFemtoDreamTrack		    *fTrack;
    AliFemtoDreamv0		    *fFemtov0;

    AliFemtoDreamEventCuts	    *fEventCuts;		// cuts for the events
    AliFemtoDreamTrackCuts	    *fTrackCutsPart1;		// cuts for the tracks of particle 1 (Protons)
    AliFemtoDreamTrackCuts	    *fTrackCutsPart2;		// cuts for the tracks of particle 2 (Antiprotons)
    AliFemtoDreamTrackCuts	    *fTrackCutsPart3;		// cuts for the tracks of particle 3 (Deuterons)
    AliFemtoDreamTrackCuts	    *fTrackCutsPart3Mass;	// cuts for the tracks of particle 3 (DeuteronsMass)
    AliFemtoDreamTrackCuts	    *fTrackCutsPart3Sigma;	// cuts for the tracks of particle 3 (DeuteronsSigma)
    AliFemtoDreamTrackCuts	    *fTrackCutsPart4;		// cuts for the tracks of particle 4 (Antideuterons)
    AliFemtoDreamTrackCuts	    *fTrackCutsPart4Mass;	// cuts for the tracks of particle 4 (AntideuteronsMass)
    AliFemtoDreamTrackCuts	    *fTrackCutsPart4Sigma;	// cuts for the tracks of particle 4 (AntideuteronsSigma)
    AliFemtoDreamv0Cuts		    *fv0CutsPart5;		// cuts for the tracks of particle 5 (Lambdas)
    AliFemtoDreamv0Cuts		    *fv0CutsPart6;		// cuts for the tracks of particle 6 (Antilambdas)

    AliFemtoDreamCollConfig	    *fConfig;			// store the configurations needed for the calculation of the correlation function
    bool fEnableEventQAPlots;
    bool fEnableResultsQAPlots;
    AliFemtoDreamPairCleaner	    *fPairCleaner;
    AliFemtoDreamPartCollection	    *fPartColl;
    AliAODTrack			    **fGTI;			// global track information (GTI)

    ClassDef(AliAnalysisTaskLeuteronAOD,1);
	
};

#endif /* PWGCF_FEMTOSCOPY_FEMTODREAM_ALIANALYSISTASKLEUTERONAOD_H_ */

