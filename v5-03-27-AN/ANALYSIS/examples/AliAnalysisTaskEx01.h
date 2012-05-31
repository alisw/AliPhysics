/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
/* AliAnalysisTaskEx01.h
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
class AliESDtrackCuts;

#ifndef ALIANALYSISTASKSE_H
#include "AliAnalysisTaskSE.h"
#endif

class AliAnalysisTaskEx01 : public AliAnalysisTaskSE {
 public:
    AliAnalysisTaskEx01();
    AliAnalysisTaskEx01(const char *name);
    virtual ~AliAnalysisTaskEx01();
    
    virtual void     UserCreateOutputObjects();
    virtual void     UserExec(Option_t *option);
    virtual void     Terminate(Option_t *);
    
 private:
    TList           *fOutput;        // Output list
    AliESDtrackCuts *fTrackCuts;     // Track cuts
    TH1F            *fHistPt;        // Pt spectrum
    TH1F            *fHistEta;       // pseudorapidity spectrum
    // NEW HISTO to be declared here
    
    AliAnalysisTaskEx01(const AliAnalysisTaskEx01&); // not implemented
    AliAnalysisTaskEx01& operator=(const AliAnalysisTaskEx01&); // not implemented
    
    ClassDef(AliAnalysisTaskEx01, 1); // example of analysis
};

#endif

