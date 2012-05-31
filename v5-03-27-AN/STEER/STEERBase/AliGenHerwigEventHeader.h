#ifndef ALIGENHERWIGEVENTHEADER_H
#define ALIGENHERWIGEVENTHEADER_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include "AliGenEventHeader.h"


class AliGenHerwigEventHeader : public AliGenEventHeader
{
 public:
    AliGenHerwigEventHeader();
    AliGenHerwigEventHeader(const char* name);
    virtual ~AliGenHerwigEventHeader() {}
    Int_t    ProcessType() const  {return fProcessType;}
    void     SetProcessType(Int_t type)  {fProcessType = type;}
    Int_t    Trials() const {return fTrials;}
    void     SetTrials(Int_t trials) {fTrials = trials;}
    Float_t  Weight() const {return fWeight;}
    void     SetWeight(Float_t weight) {fWeight = weight;}
protected:
    Int_t    fProcessType;               // HERWIG process id for this event 
    Int_t    fTrials;                    // Number of trials to fulfill trigger condition
    Float_t  fWeight;                    // Event weight (= cross section in nb for unweighted events)
    ClassDef(AliGenHerwigEventHeader, 2)  // Event header for Herwig event
};
	
	

#endif
