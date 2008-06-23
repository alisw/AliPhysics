/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               *
 **************************************************************************/

//----------------------------------------------------------------------------------
//  Class AliRsnAnalysisSimpleTask
// ------------------------
// Reader for conversion of ESD output into the internal format
// used for resonance study.
// ---
// original author: A. Pulvirenti             (email: alberto.pulvirenti@ct.infn.it)
// ---
// adapted for Analysis Framework
// by    : R. Vernet                          (email: renaud.vernet@cern.ch)
//----------------------------------------------------------------------------------

#ifndef AliRsnAnalysisSimpleTask_H
#define AliRsnAnalysisSimpleTask_H

#include "AliAnalysisTask.h"

class TList;
class AliVEvent;
class AliMCEvent;
class AliRsnPID;
class AliRsnReader;
class AliRsnEvent;
class AliRsnAnalyzerSimple;

class AliRsnAnalysisSimpleTask : public AliAnalysisTask
{
public:

    AliRsnAnalysisSimpleTask();
    AliRsnAnalysisSimpleTask(const char *name);
    virtual ~AliRsnAnalysisSimpleTask() { }
    
    // Implementation of interface methods
    virtual void ConnectInputData (Option_t *);
    virtual void CreateOutputObjects();
    virtual void Exec (Option_t *option);
    virtual void Terminate(Option_t *option);
    
    // setters
    void SetReader(AliRsnReader *reader) {fReader = reader;}
    void SetPID(AliRsnPID *pid) {fPID = pid;}
    void SetAnalyzer(AliRsnAnalyzerSimple *analyzer) {fAnalyzer = analyzer;}
    
private:

    AliRsnAnalysisSimpleTask(const AliRsnAnalysisSimpleTask&);
    AliRsnAnalysisSimpleTask& operator=(const AliRsnAnalysisSimpleTask&);

    AliVEvent*            fEvent;      // input event
	AliMCEvent*           fMC;         // corresponding MC event
    AliRsnReader*         fReader;     // read manager
    AliRsnPID*            fPID;        // particle identification manager
    AliRsnAnalyzerSimple* fAnalyzer;   // analyzer
    AliRsnEvent*          fCurrEvent;  // current event pointer -> for moving among methods
    TList*                fHistograms; // list of output histograms

    ClassDef(AliRsnAnalysisSimpleTask, 1); // implementation of RsnReader as ReaderTaskSE
};

#endif
