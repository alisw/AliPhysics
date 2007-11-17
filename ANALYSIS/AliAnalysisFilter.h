#ifndef ALIANALYSISFILTER_H
#define ALIANALYSISFILTER_H
 
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//
// Manager class for filter decisions based on cuts
// Author: Andreas Morsch
// andreas.morsch@cern.ch

#include <TNamed.h>

class AliAnalysisCuts;

class AliAnalysisFilter : public TNamed
{
 public:
    AliAnalysisFilter();
    AliAnalysisFilter(const char* name, const char* title = "AnalysisFilter");
    AliAnalysisFilter(const AliAnalysisFilter& obj);  
    virtual ~AliAnalysisFilter() {;}
    virtual UInt_t IsSelected(TObject* obj);
    virtual void AddCuts(AliAnalysisCuts* cuts);
    virtual void Init();
 private:
    TList* fCuts;    // List of cuts
    ClassDef(AliAnalysisFilter, 2); // Manager class for filter decisions
};
 
#endif
