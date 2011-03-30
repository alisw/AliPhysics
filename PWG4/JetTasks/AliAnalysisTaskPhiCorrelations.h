#ifndef AliAnalysisTaskPhiCorrelations_H
#define AliAnalysisTaskPhiCorrelations_H

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
#include "TString.h"

class AliAODEvent;
class AliAnalyseLeadingTrackUE;
class AliInputEventHandler;
class AliMCEvent;
class AliMCEventHandler;
class AliUEHistograms;
class AliVParticle;
class TH1D;
class TObjArray;
class AliEventPoolManager;

class  AliAnalysisTaskPhiCorrelations : public AliAnalysisTask
  {
  public:
    AliAnalysisTaskPhiCorrelations(const char* name="AliAnalysisTaskPhiCorrelations");
    virtual           ~AliAnalysisTaskPhiCorrelations();
       
      
    // Implementation of interace methods
    virtual     void   ConnectInputData(Option_t *);
    virtual     void   CreateOutputObjects();
    virtual     void   Exec(Option_t *option);

    // Setters/Getters
    // general configuration
    virtual     void    SetDebugLevel( Int_t level )  { fDebug = level; }
    virtual     void    SetMode(Int_t mode)           { fMode  = mode;  }
    virtual     void    SetReduceMemoryFootprint(Bool_t flag) { fReduceMemoryFootprint = flag; }
    
    // histogram settings
    void SetTrackingEfficiency( const TH1D* hist) { fkTrackingEfficiency = hist; }

    // for event QA
    void   SetTracksInVertex( Int_t val ){ fnTracksVertex = val; }
    void   SetZVertex( Double_t val )    { fZVertex = val; }
    
    // track cuts
    void   SetTrackEtaCut( Double_t val )    { fTrackEtaCut = val; }
    void   SetPtMin(Double_t val)            { fPtMin = val; }
    void   SetFilterBit( UInt_t val )        { fFilterBit = val;  }
    
    void   SetEventSelectionBit( UInt_t val )        { fSelectBit = val;  }
    void   SetUseChargeHadrons( Bool_t val ) { fUseChargeHadrons = val; }
    void   SetSelectCharge(Int_t selectCharge) { fSelectCharge = selectCharge; }
    void   SetCentralityMethod(const char* method) { fCentralityMethod = method; }

    
  private:
    AliAnalysisTaskPhiCorrelations(const  AliAnalysisTaskPhiCorrelations &det);
    AliAnalysisTaskPhiCorrelations&   operator=(const  AliAnalysisTaskPhiCorrelations &det);
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
    AliUEHistograms*  fHistos;       //! points to class to handle histograms/containers  
    AliUEHistograms*  fHistosMixed;       //! points to class to handle mixed histograms/containers  
    
    const TH1D* fkTrackingEfficiency;       // used for study of bias by tracking 

    // Handlers and events
    AliAODEvent*             fAOD;             //! AOD Event 
    TClonesArray*            fArrayMC;         //! Array of MC particles 
    AliInputEventHandler*    fInputHandler;    //! Generic InputEventHandler 
    AliMCEvent*              fMcEvent;         //! MC event
    AliMCEventHandler*       fMcHandler;       //! MCEventHandler 
    AliEventPoolManager*     fPoolMgr;         //! event pool manager
    
    // Histogram settings
    TList*              fListOfHistos;    //  Output list of containers 
    
    // Event QA cuts
    Int_t          	fnTracksVertex;        // QA tracks pointing to principal vertex (= 3 default) 
    Double_t       	fZVertex;              // Position of Vertex in Z direction
    TString             fCentralityMethod;     // Method to determine centrality
    
    // Track cuts
    Double_t      	fTrackEtaCut;          // Eta cut on particles
    Double_t            fPtMin;                // Min pT to start correlations
    UInt_t         	fFilterBit;            // Select tracks from an specific track cut (default 0xFF all track selected)
    UInt_t         	fSelectBit;            // Select events according to AliAnalysisTaskJetServices bit maps 
    Bool_t         	fUseChargeHadrons;     // Only use charge hadrons
    
    Int_t fSelectCharge;           // (un)like sign selection when building correlations: 0: no selection; 1: unlike sign; 2: like sign
    
    ClassDef( AliAnalysisTaskPhiCorrelations, 1); // Analysis task for Underlying Event analysis w.r.t. leading track
  };

#endif

    
