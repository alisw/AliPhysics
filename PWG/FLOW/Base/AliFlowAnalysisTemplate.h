/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */


#ifndef AliFlowAnalysisTemplate_H
#define AliFlowAnalysisTemplate_H

// aliroot includes
#include  "TList.h"
// forward declarations ROOT
class TH1D;
class TH1F;
class TH2D;
class TProfile;
class TDirectoryFile;
// forward declarations ALIROOT
class AliFlowTrackSimple;
class AliFlowEventSimple;
class AliFlowCommonHist;
class AliFlowCommonHistResults;

// Description: Template maker to serve as a startign point for a flow analysis.
// Author:      Redmer Alexander Bertens, Utrecht University, 2013
//              rbertens@cern.ch, rbertens@nikhef.nl, r.a.bertens@uu.nl
 
class AliFlowAnalysisTemplate {
    public:
        AliFlowAnalysisTemplate();
        virtual  ~AliFlowAnalysisTemplate();
        void Init();
        void Make(AliFlowEventSimple* anEvent);
        void GetOutputHistograms(TList *outputListHistos);
        void Finish();
        void WriteHistograms(TDirectoryFile *outputFileName) const;
        // setters
        void SetHarmonic(Int_t h)               { fHarmonic = h; }
        void SetApplyCorrectionForNUA(Bool_t a) { fApplyCorrectionForNUA = a; }
        void SetDebug(Bool_t d)                 { fDebug = d; }
        void SetUsePhiWeights(Bool_t p)         { fUsePhiWeights = p; }
        void SetWeightsList(TList* const w)     { fWeightsList = (TList*)w->Clone(); }
        // getters
        TList*    GetHistList() const           { return fHistList; }
        AliFlowCommonHist*        GetCommonHists()    const { return fCommonHists; }
        AliFlowCommonHistResults* GetCommonHistsRes() const { return fCommonHistsRes; }
    private:
        Int_t fDebug;                   // debug flag 
        Int_t fUsePhiWeights;           // use phi weights
        Int_t fApplyCorrectionForNUA;   // apply correction for non-uniform acceptance
        Int_t fHarmonic;                // harmonic 
        TList*     fWeightsList;        // list holding input histograms with phi weights
        TList*     fHistList;           // list to hold all output histograms  
        AliFlowCommonHist*              fCommonHists;           // control histograms
        AliFlowCommonHistResults*       fCommonHistsRes;        // results histograms

        AliFlowAnalysisTemplate(const AliFlowAnalysisTemplate& a);              // not implemented
        AliFlowAnalysisTemplate& operator=(const AliFlowAnalysisTemplate& a);   // not implemented
        ClassDef(AliFlowAnalysisTemplate, 0)
};
 

#endif
