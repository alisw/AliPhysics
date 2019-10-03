#ifndef AliAnalysisTaskCorrelationsDev_H
#define AliAnalysisTaskCorrelationsDev_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////////////////////////////////////////////////////////////////////////
//
// Analysis class for development of correlation analysis code
//
//    Authors:
//    Jan Fiete Grosse-Oetringhaus
// 
////////////////////////////////////////////////////////////////////////

#include "AliAnalysisTaskSE.h"
#include "TString.h"

class TNtuple;

class AliAnalyseLeadingTrackUE;

class  AliAnalysisTaskCorrelationsDev : public AliAnalysisTaskSE
{
public:
  AliAnalysisTaskCorrelationsDev(const char* name="AliAnalysisTaskCorrelationsDev");
  virtual           ~AliAnalysisTaskCorrelationsDev();

  // Implementation of interace methods
  virtual     void   UserCreateOutputObjects();
  virtual     void   UserExec(Option_t *option);

  // Setters/Getters
  // general configuration
  virtual     void    SetDebugLevel( Int_t level )  { fDebug = level; }

  void   SetZVertex( Double_t val )    { fZVertex = val; }

  // track cuts
  void   SetTrackEtaCut( Double_t val )    { fTrackEtaCut = val; }
  void   SetTrackEtaCutMin( Double_t val ) { fTrackEtaCutMin = val; }
  void   SetPtMin(Double_t val)            { fPtMin = val; }
  void   SetFilterBit( UInt_t val )        { fFilterBit = val;  }
  void   SetTrackStatus(UInt_t status)     { fTrackStatus = status; }

  void   SetEventCuts(AliAnalysisCuts* eventCuts) { fEventCuts = eventCuts; }
  void   SetEventSelectionBit( UInt_t val )        { fSelectBit = val;  }
  void   SetCentralityMethod(const char* method, Bool_t unchecked = kFALSE) { fCentralityMethod = method; fUseUncheckedCentrality = kFALSE; }
  void   SetUseNewCentralityFramework(Bool_t flag) { fUseNewCentralityFramework = flag; }
  
  void   SetCustomTree(TString eventConfig, Int_t maxTracks, TString trackConfig) { fTreeEventConfig = eventConfig.Tokenize(";"); fTreeMaxTracks = maxTracks; fTreeTrackConfig = trackConfig.Tokenize(";"); }

private:
  AliAnalysisTaskCorrelationsDev(const  AliAnalysisTaskCorrelationsDev &det);
  AliAnalysisTaskCorrelationsDev&   operator=(const  AliAnalysisTaskCorrelationsDev &det);
  void            AddSettingsTree();                                  // add list of settings to output list
  // Analysis methods
  void            Initialize(); 			                // initialize some common pointer
  Double_t        GetCentrality(AliVEvent* inputEvent, TObject* mc);
  Float_t         GetValueInterpreted(TObject* source, const char* command);

  AliAnalyseLeadingTrackUE*     fAnalyseUE;      //! points to class containing common analysis algorithms
  
  // General configuration
  Int_t               fDebug;           //  Debug flag

  // Histogram settings
  TList*              fListOfHistos;    //  Output list of containers

  // Event
  AliAnalysisCuts* fEventCuts;          // additional cuts which have to be passed
  UInt_t     fSelectBit;                // Select events according to AliAnalysisTaskJetServices bit maps
  Double_t   fZVertex;                  // Position of Vertex in Z direction
  TString    fCentralityMethod;         // Method to determine centrality
  Bool_t     fUseUncheckedCentrality;   // use unchecked centrality; default: kFALSE
  Bool_t     fUseNewCentralityFramework;// use the AliMultSelection framework

  // Track cuts
  Double_t   fTrackEtaCut;              // Maximum Eta cut on particles
  Double_t   fTrackEtaCutMin;           // Minimum Eta cut on particles
  Double_t   fPtMin;                    // Min pT to start correlations
  UInt_t     fFilterBit;                // Select tracks from an specific track cut
  UInt_t     fTrackStatus;              // if non-0, the bits set in this variable are required for each track
  
  // custom tree
  TNtuple*   fOutputTree;               //! custom tree output
  Float_t*   fOutputContainer;          //! to fill the tuple 
  TObjArray* fTreeEventConfig;          // list of fields from header kept in tree
  Int_t      fTreeMaxTracks;            // max number of tracks in tree
  TObjArray* fTreeTrackConfig;          // list of fields from track kept in tree

  ClassDef(AliAnalysisTaskCorrelationsDev, 3); // Analysis task for correlation development
};

#endif


