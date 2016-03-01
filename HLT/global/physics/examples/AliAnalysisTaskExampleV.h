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

#ifndef ALIANALYSISTASKSE_H
#include "AliAnalysisTaskSE.h"
#endif

class AliAnalysisTaskExampleV : public AliAnalysisTaskSE {
 public:
    AliAnalysisTaskExampleV();
    AliAnalysisTaskExampleV(const char *name);
    virtual ~AliAnalysisTaskExampleV();
    
    virtual void     UserCreateOutputObjects();
    virtual void     UserExec(Option_t *option);
    virtual void     Terminate(Option_t *);
    
 private:
    TList           *fOutput;        // Output list
    TH1F            *fHistPt;        // Pt spectrum
    TH1F            *fHistEta;       // pseudorapidity spectrum
    // NEW HISTO to be declared here
    
    AliAnalysisTaskExampleV(const AliAnalysisTaskExampleV&); // not implemented
    AliAnalysisTaskExampleV& operator=(const AliAnalysisTaskExampleV&); // not implemented
    
    ClassDef(AliAnalysisTaskExampleV, 1); // example of analysis
};

#endif

