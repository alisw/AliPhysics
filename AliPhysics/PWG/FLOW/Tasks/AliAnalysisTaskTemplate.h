/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice */
/* $Id: $ */

#ifndef AliAnalysisTaskTemplate_H
#define AliAnalysisTaskTemplate_H

// author: Redmer Alexander Bertens, Utrecht University, 2013
// rbertens@cern.ch, rbertens@nihef.nl, r.a.bertens@uu.nl
class AliFlowEventSimple;
class AliFlowAnalysisTemplate;
class TList;

#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskTemplate : public AliAnalysisTaskSE {
    public:
        AliAnalysisTaskTemplate();
        AliAnalysisTaskTemplate(const char *name, Bool_t usePhiWeights);
        virtual         ~AliAnalysisTaskTemplate();
        virtual void    UserCreateOutputObjects();
        virtual void    UserExec(Option_t *option);
        virtual void    Terminate(Option_t *);
        // setters
        void            SetApplyCorrectionForNUA(Bool_t const a)        {fApplyCorrectionForNUA         = a;}
        void            SetHarmonic(Int_t const h)                      {fHarmonic = h;}
        // getters
        Bool_t          GetApplyCorrectionForNUA() const                {return fApplyCorrectionForNUA;}
        Int_t GetHarmonic() const                                       {return fHarmonic;};   
    private:
        AliFlowAnalysisTemplate*        fFlowTask;      // analysis object
        TList*                          fOutputList;    // collection of output
        Bool_t                          fUsePhiWeights; // use phi weights
        TList*                          fListWeights;   // list with weights
        Bool_t                          fApplyCorrectionForNUA;         // apply automatic correction for non-uniform acceptance 
        Int_t                           fHarmonic;      // harmonic

        AliAnalysisTaskTemplate(const AliAnalysisTaskTemplate& aAnalysisTask);                  // not implemented
        AliAnalysisTaskTemplate& operator=(const AliAnalysisTaskTemplate& aAnalysisTask);       // not implemented

        ClassDef(AliAnalysisTaskTemplate, 1); // example of analysis
};
#endif
