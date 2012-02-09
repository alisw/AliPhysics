#ifndef ALIANALYSISCUTS_H
#define ALIANALYSISCUTS_H
 
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

// Base class for analysis cuts
// Author Andreas Morsch
// andreas.morsch@cern.ch

#include "AliVCuts.h"

class TList;
class TCollection;

class AliAnalysisCuts : public AliVCuts
{
 public:
    AliAnalysisCuts();
    AliAnalysisCuts(const char* name, const char* title);
    AliAnalysisCuts(const AliAnalysisCuts& obj);
    AliAnalysisCuts& operator=(const AliAnalysisCuts& obj);
    virtual ~AliAnalysisCuts() {;}
    virtual Bool_t IsSelected(TObject* /* obj  */ )  {return kFALSE;}
    virtual Bool_t IsSelected(TList*   /* list */ ) = 0;
    virtual void   Init() {;}
    virtual void   SetFilterMask(UInt_t mask) {fFilterMask = mask;}
    virtual UInt_t GetFilterMask()   const    {return fFilterMask;}
    virtual void   SetSelected(Bool_t dec)    {fSelected = dec;}
    virtual UInt_t Selected()        const    {return fSelected;}	    
    virtual Long64_t Merge(TCollection* /* list */)      { return 0; }
    virtual TList* GetOutput()                { return 0; }
    virtual TObject *GetStatistics(Option_t *) const {return 0;}
 private:
    UInt_t fFilterMask; // Mask to use one of the previous decisions inside a filter
    Bool_t fSelected;   // Final decision on selction
    ClassDef(AliAnalysisCuts, 5); // Base class for filter decisions on ESD objects
};
 
#endif
