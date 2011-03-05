#ifndef ALIEVENTPOOLLOOP_H
#define ALIEVENTPOOLLOOP_H
 
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

// Realisation of an AliVEventPool which allows the user to
// run the analysis in a loop, i.e. passing several times over 
// the same event chain.
// Author Andreas Morsch
// andreas.morsch@cern.ch

#include <AliVEventPool.h>
class AliRunTagCuts;
class AliLHCTagCuts;
class AliDetectorTagCuts;
class AliEventTagCuts;
class AliTagAnalysis;
class TChain;

class AliEventPoolLoop : public AliVEventPool
{
 public:
    AliEventPoolLoop();
    AliEventPoolLoop(Int_t nit);
    AliEventPoolLoop(const char* name, const char* title);

    virtual ~AliEventPoolLoop() {;}
    // Interface
    virtual TChain* GetNextChain();
    virtual void  GetCurrentBin(Float_t* /*bin*/);
    virtual Int_t GetDimension();
    virtual void  Init();
    virtual Int_t BinNumber() const {return fNIteration;}
	    
 private:
    AliEventPoolLoop(const AliEventPoolLoop& obj);
    AliEventPoolLoop& operator=(const AliEventPoolLoop& other);
 protected:
    Int_t fMaxIterations; // Maximum number of iterations 
    Int_t fNIteration;    // Number of iterations
    TChain* fChainClone; // Clone of the original 
    ClassDef(AliEventPoolLoop, 0); 
};
 
#endif
