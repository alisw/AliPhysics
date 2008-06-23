/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               *
 **************************************************************************/

//-------------------------------------------------------------------------
//                      Class AliRsnAnalysisSimple
//             Reconstruction and analysis of K* Rsn
// ........................................
// ........................................
// ........................................
// ........................................
//
// author: A. Pulvirenti             (email: alberto.pulvirenti@ct.infn.it)
//-------------------------------------------------------------------------

#ifndef AliRsnAnalysisSimple_H
#define AliRsnAnalysisSimple_H

#include "AliRsnPID.h"

class TTree;
class TArrayI;
class TObjArray;
class AliRsnPairSimple;
class AliRsnPID;
class AliRsnEventBuffer;
class AliRsnAnalyzerSimple;

class AliRsnAnalysisSimple : public TObject
{
public:

    AliRsnAnalysisSimple(AliRsnAnalyzerSimple *ana = 0x0, AliRsnPID *pid = 0x0);
    virtual ~AliRsnAnalysisSimple() {Clear();}
    virtual void Clear(Option_t *option = "C");

    /* setters */
    void    SetAnalyzer(AliRsnAnalyzerSimple *analyzer) {fAnalyzer = analyzer;}
    void    SetEventsTree(TTree *tree);
    void    SetFileName(char *fname) {strcpy(fFileName, fname);}
    void    SetStep(Int_t step) {fStep = step;}
    void    SetPID(AliRsnPID *pid) {fPID = pid;}

    /* working routines */
    Bool_t  Initialize();
    Stat_t  Process();
    void    SaveOutput() const;

private:

    AliRsnAnalysisSimple(const AliRsnAnalysisSimple &copy) :
      TObject(copy),fInitialized(kFALSE),fStep(1000),fTree(0x0),fPID(0x0),fAnalyzer(0x0) { }
    AliRsnAnalysisSimple& operator=(const AliRsnAnalysisSimple & /*copy*/) { return (*this); }

    Bool_t                fInitialized;     // flag to check initialization
    Int_t                 fStep;            // progress step
    Char_t                fFileName[250];   // output file name
	TTree                *fTree;            //! TTree of events
    AliRsnPID            *fPID;             //! PID object
	AliRsnAnalyzerSimple *fAnalyzer;        //! analyzer

	// Rsn analysis implementation
	ClassDef(AliRsnAnalysisSimple,1)
};

#endif
