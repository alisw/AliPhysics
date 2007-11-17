#ifndef ALIANALYSISCUTS_H
#define ALIANALYSISCUTS_H
 
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

// Base class for analysis cuts
// Author Andreas Morsch
// andreas.morsch@cern.ch

#include <TNamed.h>

class AliAnalysisCuts : public TNamed
{
 public:
    AliAnalysisCuts();
    AliAnalysisCuts(const char* name, const char* title);
    AliAnalysisCuts(const AliAnalysisCuts& obj);  
    virtual ~AliAnalysisCuts() {;}
    virtual Bool_t IsSelected(TObject* obj) = 0;
    virtual void   Init() {;}
 private:
    ClassDef(AliAnalysisCuts, 2); // Base class for filter decisions on ESD objects
};
 
#endif
