#ifndef ALIANALYSISTASKLEADINGTRACKUE_H
#define ALIANALYSISTASKLEADINGTRACKUE_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////////////////////////////////////////////////////////////////////////
//
// Analysis class for Underlying Event studies w.r.t. leading track
//
// Look for correlations on the tranverse regions w.r.t
// the leading track in the event
//
// This class needs input AODs.
// The output is a list of analysis-specific containers.
//
// The AOD can be either connected to the InputEventHandler  
// for a chain of AOD files 
// or 
// to the OutputEventHandler
// for a chain of ESD files,
// in this case the class should be in the train after the jet-finder
//
//    Authors:
//    Arian Abrahantes Quintana 
//    Jan Fiete Grosse-Oetringhaus
//    Ernesto Lopez Torres
//    Sara Vallero
// 
////////////////////////////////////////////////////////////////////////

#include "AliAnalysisTask.h"
#include "AliUEHist.h"

class AliAODEvent;
class AliAnalyseLeadingTrackUE;
class AliInputEventHandler;
class AliMCEvent;
class AliMCEventHandler;
class AliUEHistograms;
class AliVParticle;
class TH1D;
class TObjArray;

class  AliAnalysisTaskLeadingTrackUE : public AliAnalysisTask
  {
  public:
    AliAnalysisTaskLeadingTrackUE(const char* name="AliAnalysisTaskLeadingTrackUE");
    virtual           ~AliAnalysisTaskLeadingTrackUE();
       
      
    // Implementation of interace methods
    virtual     Bool_t Notify();
    virtual     void   ConnectInputData(Option_t *);
    virtual     void   CreateOutputObjects();
    virtual     void   Exec(Option_t *option);
    virtual     void   Terminate(Option_t *);

    void FillReducedEfficiency(Int_t eventId, AliUEHist::CFStep step,const TObjArray* ltRECO, Bool_t twoStep);

    // Setters/Getters
    // general configuration
    virtual     void    SetDebugLevel( Int_t level )  { fDebug = level; }
    virtual     void    SetMode(Int_t mode)           { fMode  = mode;  }
    virtual     void    SetReduceMemoryFootprint(Bool_t flag) { fReduceMemoryFootprint = flag; }
    
    // histogram settings
    void   		SetPtRangeInHist( Int_t bin, Double_t min, Double_t max ) {
      				fBinsPtInHist = bin; 
      				fMinJetPtInHist = min; 
      				fMaxJetPtInHist = max; 
    				}
    void SetTrackingEfficiency( const TH1D* hist) { fkTrackingEfficiency = hist; }

    // for event QA
    void   SetTracksInVertex( Int_t val ){ fnTracksVertex = val; }
    void   SetZVertex( Double_t val )    { fZVertex = val; }
    
    // track cuts
    void   SetTrackEtaCut(Double_t val)       { fTrackEtaCut = val; }
    void   SetLeadingTrackEtaCut( Double_t val )    { fLeadingTrackEtaCut = val; }
    void   SetFilterBit( UInt_t val )        { fFilterBit = val;  }
    void   SetEventSelectionBit( UInt_t val )        { fSelectBit = val;  }
    void   SetUseChargeHadrons( Bool_t val ) { fUseChargeHadrons = val; }
    
  protected:
  static AliAnalysisTaskLeadingTrackUE*     fgTaskLeadingTrackUE;       // Pointer to single instance

  private:
    AliAnalysisTaskLeadingTrackUE(const  AliAnalysisTaskLeadingTrackUE &det);
    AliAnalysisTaskLeadingTrackUE&   operator=(const  AliAnalysisTaskLeadingTrackUE &det);
    void            AddSettingsTree();                                  // add list of settings to output list
    // Analysis methods
    void            AnalyseCorrectionMode();                            // main algorithm to get correction maps
    void            AnalyseDataMode();                                  // main algorithm to get raw distributions
    void            Initialize(); 			                // initialize some common pointer



    // General configuration
    Int_t               fDebug;           //  Debug flag
    Int_t 	        fMode;            //  fMode = 0: data-like analysis 
    				          //  fMode = 1: corrections analysis	
    Bool_t              fReduceMemoryFootprint; // reduce memory consumption by writing less debug histograms
    
    // Pointers to external UE classes
    AliAnalyseLeadingTrackUE*     fAnalyseUE;      //! points to class containing common analysis algorithms
    AliUEHistograms*  fHistosUE;       //! points to class to handle histograms/containers  
    
    const TH1D* fkTrackingEfficiency;       // used for study of bias by tracking 

    // Handlers and events
    AliAODEvent*             fAOD;             //! AOD Event 
    TClonesArray*            fArrayMC;         //! Array of MC particles 
    AliInputEventHandler*    fInputHandler;    //! Generic InputEventHandler 
    AliMCEvent*              fMcEvent;         //! MC event
    AliMCEventHandler*       fMcHandler;       //! MCEventHandler 
    
    // Histogram settings
    TList*              fListOfHistos;    //  Output list of containers 
    Int_t          	fBinsPtInHist;     //  # bins for pT histos
    Double_t       	fMinJetPtInHist;   //  min Jet Pt value for histo range
    Double_t       	fMaxJetPtInHist;   //  max Jet Pt value for histo range
    
    
    // Event QA cuts
    Int_t          	fnTracksVertex;        // QA tracks pointing to principal vertex (= 3 default) 
    Double_t       	fZVertex;              // Position of Vertex in Z direction
    
    // Track cuts
    Double_t            fTrackEtaCut;          // Eta cut on inclusive tracks
    Double_t      	fLeadingTrackEtaCut;   // Eta cut on leading track
    UInt_t         	fFilterBit;            // Select tracks from an specific track cut (default 0xFF all track selected)
    UInt_t         	fSelectBit;            // Select events according to AliAnalysisTaskJetServices bit maps 
    Bool_t         	fUseChargeHadrons;     // Only use charge hadrons
    
    // MC cross-section 
    Double_t      	fAvgTrials;        // average trials used to fill the fh1Trials histogram in case we do not have trials on a event by event basis

    ClassDef( AliAnalysisTaskLeadingTrackUE, 1); // Analysis task for Underlying Event analysis w.r.t. leading track
  };

#endif

    
