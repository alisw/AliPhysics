#ifndef ALIVEVENTPOOL_H
#define ALIVEVENTPOOL_H
 
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliAnalysisCuts.h 25775 2008-05-15 08:08:39Z morsch $ */

// Base class for event pool.
// This class is needed by the AnalysisManager to steer a mixing analysis.
// Author Andreas Morsch
// andreas.morsch@cern.ch

#include <TNamed.h>
class TChain;

class AliVEventPool : public TNamed
{
 public:
    AliVEventPool();
    AliVEventPool(const char* name, const char* title);

    virtual ~AliVEventPool() {;}
    // Interface
    virtual TChain* GetNextChain()                 = 0;
    virtual void   GetCurrentBin(Float_t* /*bin*/) = 0;
    virtual Int_t  GetDimension()                  = 0;
    virtual void   Init()                          = 0;
    virtual Int_t  BinNumber() const               = 0;
    // Basic functionality
    virtual void SetChain(TChain* chain) {fChain = chain;}
    AliVEventPool(const AliVEventPool& obj);
    AliVEventPool& operator=(const AliVEventPool& other);
 protected:
    TChain* fChain;   // Chain of trees
    ClassDef(AliVEventPool, 1); 
};
 
#endif
