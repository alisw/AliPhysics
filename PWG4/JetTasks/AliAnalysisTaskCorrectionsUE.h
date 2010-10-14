#ifndef ALIANALYSISTASKCORRECTIONSUE_H
#define ALIANALYSISTASKCORRECTIONSUE_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////////////////////////////////////////////////////////////////////////
//
// Analysis class to Correct Underlying Event studies
//
// This class needs as input ESDs.
// The output is an analysis-specific container.
//
// The class is used to get the contamination from secondaries
// from tracks DCA distribution 
// as function of track pT and pseudo-rapidity.
// It provides additional information for the corrections 
// that can not be retrieved by the AliAnalysisTaskLeadingTackUE
// task, which is running on AODs.
// 
////////////////////////////////////////////////////////////////////////

#include "AliAnalysisTask.h"

class AliAnalyseLeadingTrackUE;
class AliESDtrackCuts;
class AliInputEventHandler;
class AliAODEvent;
class AliCFContainer;
class AliESDEvent;
class AliESDtrack;
class AliAODInputHandler;
class AliESDInputHandler;
class AliMCEventHandler;
class AliMCEvent;
class TH1F;
class TH2F;
class TH3F;
class TH1I;
class TProfile;
class TTree;
class TVector3;

class  AliAnalysisTaskCorrectionsUE : public AliAnalysisTask
  {
  public:
    // track cuts steps
    enum CFSteps {
    	kCFStepAll     = 0,
	kCFStepCuts    = 1, // standard ITS+TPC 2009 cuts w.o. DCA cut and SPD cluster requirement
	kCFStepSPD     = 2, // add SPD cluster requirement
	kCFStepDCA     = 3  // add pT dependent DCA cut
    	};

    AliAnalysisTaskCorrectionsUE(const char* name="AliAnalysisTaskCorrectionsUE");
    virtual           ~AliAnalysisTaskCorrectionsUE() {if ( fListOfHistos ) delete fListOfHistos; }
    AliAnalysisTaskCorrectionsUE(const  AliAnalysisTaskCorrectionsUE &det);
    AliAnalysisTaskCorrectionsUE&   operator=(const  AliAnalysisTaskCorrectionsUE &det);
       
    // return instance of the singleton
    static  AliAnalysisTaskCorrectionsUE* Instance();
      
    // Implementation of interace methods
    virtual     Bool_t Notify();
    virtual     void   ConnectInputData(Option_t *);
    virtual     void   CreateOutputObjects();
    virtual     void   Exec(Option_t *option);
    virtual     void   Terminate(Option_t *);

    //  Setters/Getters
    virtual     void   SetDebugLevel( Int_t level )  { fDebug = level; }
    virtual     void   SetMode(Int_t mode)           { fMode  = mode;  }

    //Event QA
    void  SetZVertex( Double_t val )          { fZVertex = val; }
    void  SetTracksInVertex( Int_t val ){ fnTracksVertex = val; }

    // Track selection cuts
    void  SetTrackEtaCut( Double_t val ) { fTrackEtaCut = val; }
    void  SetTrackPtCut( Double_t val )  { fTrackPtCut = val; }
	
  protected:
    static     AliAnalysisTaskCorrectionsUE*     fgTaskCorrectionsUE;        // Pointer to single instance
  private:
    void       AddSettingsTree();                                  // add list of settings to output list
    // Analysis methods
    void       AnalyseCorrectionMode();                            // main algorithm to get correction maps
    void       CreateContainer();  			           // create the output CF container
    void       FillContainer(AliESDtrack *track, Int_t step,Bool_t mcvertex, Double_t matchLeading); // fill container 
    AliAnalyseLeadingTrackUE*       fAnalyseUE;         //! points to AliAnalyseLeadingTrackUE class
    Int_t                           fDebug;             //  Debug flag
    AliESDEvent*                    fESDEvent;          //! ESD Event
    AliESDInputHandler*             fESDHandler;        //! ESD Input Handler
    AliInputEventHandler*           fInputHandler;      //  Input event handler
    TList*                          fListOfHistos;      //  Output list of histograms
    AliMCEvent*                     fMcEvent;           //  pointer to MC event
    AliMCEventHandler*              fMcHandler;         //  pointer to MC handler
    Int_t 	                    fMode;              //  fMode = 0: data-like analysis 
    				                        //  fMode = 1: corrections analysis	
    AliCFContainer*		    fOutCFcont;         //  output CF container   
    TH1F*			    fhEntries;          //  count events	 
    TH1F* 			    fhFakes;           	//  counts the amount of fake tracks 
    TH1F*                           fhPtMCAll;          //  pT distribution of all accepted MC tracks 
    TH1F*                           fhPtMCPrim;         //  pT distribution MC primaries
    TH1F*                           fhPtMCSec;          //  pT distribution MC secondaries
    TH1F*                           fhPtMCPrimFake;     //  pT distribution MC fake primaries
    TH1F*                           fhPtMCSecFake;      //  pT distribution MC fake secondaries

    // Cuts Events QA
    Int_t          fnTracksVertex;        // QA tracks pointing to principal vertex (= 3 default) 
    Double_t       fZVertex;              // Position of Vertex in Z direction
    TH1F*          fhVertexContributors;  // Plot number of contributors in vertex  
    TH3F*          fhVertexReso;          //  vertex resolution in XY and Z vs. number of contributors
    // Track cuts
    Double_t       fTrackEtaCut;          // Eta cut on tracks in the regions (fRegionType=1)
    Double_t       fTrackPtCut;           // Pt cut of tracks in the regions
    AliESDtrackCuts* fEsdTrackCuts;       // ITS+TPC 2009 cuts (no SPD requirement, no DCA cut) 
    AliESDtrackCuts* fEsdTrackCutsSPD;     // Require 1 cluser in SPD
    AliESDtrackCuts* fEsdTrackCutsSDD;     // Require 1 cluser in 1st layer SDD
    AliESDtrackCuts* fEsdTrackCutsDCA;     // Add pT dependent DCA cut
    ClassDef( AliAnalysisTaskCorrectionsUE, 5); // Analysis task to correct Underlying Event analysis
  };

#endif

    
