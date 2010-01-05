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
    Int_t    ProcessType()  {return fProcessType;}
    void     SetProcessType(Int_t type)  {fProcessType = type;}
    Int_t    Trials() {return fTrials;}
    void     SetTrials(Int_t trials) {fTrials = trials;}
protected:
    Int_t    fProcessType;               // HERWIG process id for this event 
    Int_t    fTrials;                    // Number of trials to fulfill trigger condition
    ClassDef(AliGenHerwigEventHeader, 1)  // Event header for Herwig event
};
	
	

#endif
