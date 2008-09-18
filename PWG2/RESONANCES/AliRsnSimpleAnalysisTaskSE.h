/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               *
 **************************************************************************/

//----------------------------------------------------------------------------------
//  Class AliRsnSimpleAnalysisTaskSE
// ------------------------
// Reader for conversion of ESD output into the internal format
// used for resonance study.
// ---
// original author: A. Pulvirenti             (email: alberto.pulvirenti@ct.infn.it)
// ---
// adapted for Analysis Framework
// by    : R. Vernet                          (email: renaud.vernet@cern.ch)
//----------------------------------------------------------------------------------

#ifndef AliRsnSimpleAnalysisTaskSE_H
#define AliRsnSimpleAnalysisTaskSE_H

#include "AliAnalysisTaskSE.h"

class TList;
class AliVEvent;
class AliMCEvent;
class AliRsnReader;
class AliRsnEvent;
class AliRsnSimpleAnalyzer;

class AliRsnSimpleAnalysisTaskSE : public AliAnalysisTaskSE
{
  public:

    AliRsnSimpleAnalysisTaskSE();
    AliRsnSimpleAnalysisTaskSE(const char *name);
    virtual ~AliRsnSimpleAnalysisTaskSE() { }

    // Implementation of interface methods
    virtual void UserCreateOutputObjects();
    virtual void UserExec(Option_t *option);

    // setters
    Bool_t Configure(const char *configFile = "RsnConfig.C");
    void   SetReader(AliRsnReader *reader) {fReader = reader;}
    void   SetPID(AliRsnPID *pid) {fPID = pid;}
    void   SetAnalyzer(AliRsnSimpleAnalyzer *analyzer) {fAnalyzer = analyzer;}
    void   PrintSettings();

  private:

    AliRsnSimpleAnalysisTaskSE(const AliRsnSimpleAnalysisTaskSE&) :
        AliAnalysisTaskSE(),fReader(0x0),fPID(0x0),fAnalyzer(0x0),fRsnEvent(0x0),fHistograms(0x0)
    { /*nothing*/ }
    AliRsnSimpleAnalysisTaskSE& operator=(const AliRsnSimpleAnalysisTaskSE&)
    { /*nothing*/ return (*this); }

    AliRsnReader*         fReader;     // read manager
    AliRsnPID*            fPID;        // PID manager
    AliRsnSimpleAnalyzer* fAnalyzer;   // analyzer
    AliRsnEvent*          fRsnEvent;   // current event pointer -> for moving among methods
    TList*                fHistograms; // list of output histograms

    ClassDef(AliRsnSimpleAnalysisTaskSE, 1); // implementation of RsnReader as ReaderTaskSE
};

#endif
