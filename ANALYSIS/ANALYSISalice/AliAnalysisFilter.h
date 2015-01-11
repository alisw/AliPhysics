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
class TList;

class AliAnalysisFilter : public TNamed
{
 public:
    AliAnalysisFilter();
    AliAnalysisFilter(const char* name, const char* title = "AnalysisFilter");
    AliAnalysisFilter(const AliAnalysisFilter& obj);
    AliAnalysisFilter& operator=(const AliAnalysisFilter& other);
    virtual ~AliAnalysisFilter();
    virtual UInt_t IsSelected(TObject* obj);
    virtual UInt_t IsSelected(TList* obj);
    virtual Bool_t IsSelected(char* name);
    virtual void AddCuts(AliAnalysisCuts* cuts);
    virtual void Init();
    TList*  GetCuts() const {return fCuts;}
	    
 private:
    TList* fCuts;    // List of cuts
    ClassDef(AliAnalysisFilter, 2); // Manager class for filter decisions
};
 
#endif
