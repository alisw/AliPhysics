/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
/* AliAnalysisTaskExampleV.h
 *
 * Template task producing a P_t spectrum and pseudorapidity distribution.
 * Includes explanations of physics and primary track selections
 *
 * Based on tutorial example from offline pages
 * Edited by Arvinder Palaha
 */
#ifndef ALIANALYSISTASKEX01_H
#define ALIANALYSISTASKEX01_H

class TH1F;
class TList;
class AliVEvent;

#ifndef ALIANALYSISTASKSE_H
#include "AliAnalysisTask.h"
#endif

class AliAnalysisTaskExampleV : public AliAnalysisTask {
 public:
    AliAnalysisTaskExampleV();
    AliAnalysisTaskExampleV(const char *name);
    virtual ~AliAnalysisTaskExampleV();
    
    virtual void     CreateOutputObjects();
    virtual void     Exec(Option_t *option);
    virtual void     Terminate(Option_t *);
    virtual void     ConnectInputData(Option_t*);
    
 private:
    TList           *fOutput;        //! Output list
    TH1F            *fHistPt;        //! Pt spectrum
    TH1F            *fHistEta;       //! pseudorapidity spectrum
    Bool_t          fSkipExec;       // example config (don't exclude from streaming!)
    // NEW HISTO to be declared here
    
    AliVEvent *fV; //!
    
    AliAnalysisTaskExampleV(const AliAnalysisTaskExampleV&); // not implemented
    AliAnalysisTaskExampleV& operator=(const AliAnalysisTaskExampleV&); // not implemented
    
    ClassDef(AliAnalysisTaskExampleV, 1); // example of analysis
};

#endif

