/*************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               *
 **************************************************************************/

//-------------------------------------------------------------------------
//                      Class AliRsnAnalysis
//             Reconstruction and analysis of K* Rsn
// ........................................
// ........................................
// ........................................
// ........................................
//
// author: A. Pulvirenti             (email: alberto.pulvirenti@ct.infn.it)
//-------------------------------------------------------------------------

#ifndef AliRsnPairSimple_H
#define AliRsnPairSimple_H

#include <TNamed.h>
#include "AliRsnPID.h"
#include "AliRsnPairParticle.h"

class TH1D;
class AliRsnEvent;
class AliRsnPairDef;
class AliRsnCutMgr;
class AliRsnPID;

class AliRsnPairSimple : public TNamed
{

public:

    AliRsnPairSimple(AliRsnPairDef *pd = 0x0, const char *name  = "", const char *title = "");
    virtual ~AliRsnPairSimple() {Clear();}
    virtual void Clear(Option_t *option = "");

    // getters
    TH1D*               GetHistogram() {return fHistogram;}
    TH1D*               GetHistogramMC() {return fHistogramMC;}
    AliRsnPairDef*      GetPairDef() {return fPairDef;}
    AliRsnPairParticle* GetPair() {return &fPair;}
    Bool_t              StoreOnlyTruePairs() const {return fStoreOnlyTrue;}
    Bool_t              IsForMixing() const {return fForMixing;}

    // setters
    void   SetPIDMethod(AliRsnPID::EMethod method) {fPIDMethod = method;}
    void   SetPairDef(AliRsnPairDef *def) {fPairDef = def;}
    void   SetStoreOnlyTrue(Bool_t doit = kTRUE) {fStoreOnlyTrue = doit;}
    void   SetForMixing(Bool_t doit = kTRUE) {fForMixing = doit;}
    void   SetCutManager(AliRsnCutMgr *cutMgr) {fCuts = cutMgr;}

    // working parameters
    void   InitHistograms();
    Stat_t Process(AliRsnEvent *event1, AliRsnEvent *event2 = 0);
    Stat_t Process(AliRsnDaughter *t1, AliRsnDaughter *t2);

private:

    // private functions
    AliRsnPairSimple(const AliRsnPairSimple &copy);
    const AliRsnPairSimple& operator=(const AliRsnPairSimple &copy);
    const char* GetHistName();
    const char* GetHistTitle();

    // flags
    AliRsnPID::EMethod   fPIDMethod;       // flag to know PID method used
    Bool_t               fForMixing;       // flag is true for objects created for event mixing
    Bool_t               fStoreOnlyTrue;   // output = only spectra of true pairs

    // cut manager
    AliRsnCutMgr        *fCuts;            // cut manager for single particle

    // objects
    AliRsnPairParticle   fPair;            // utility class for pair
    AliRsnPairDef       *fPairDef;         // definitions for pair
    TH1D                *fHistogram;       // invariant mass distribution
    TH1D                *fHistogramMC;     // invariant mass distribution (MC)

    // ROOT dictionary
    ClassDef(AliRsnPairSimple, 1)
};

#endif
