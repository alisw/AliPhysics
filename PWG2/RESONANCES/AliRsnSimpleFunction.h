//
// Class AliRsnSimpleFcn
//
// This class defines a base classe to implement a typical computation
// which uses the internal RSN package event format (AliRsnEvent).
// It contains some default flags which turn out to be useful:
//  - a flag to select only the "true" pairs (tracks from same resonance)
//  - a flag to know if the computation is done over two events (mixing)
//
// Any kind of analysis object should be implemented as inheriting from this
// because the AliRsnAnalyzer which executes the analysis will accept a collection
// of such objects, in order to have a unique format of processing method
//
// The user who implements a kind of computation type should inherit from 
// this class and override the virtual functions defined in it, which 
// initialize the final output histogram and define how to process data.
//
//
// author: A. Pulvirenti             (email: alberto.pulvirenti@ct.infn.it)
//

#ifndef ALIRSNSIMPLEFCN_H
#define ALIRSNSIMPLEFCN_H

#include <TNamed.h>

#include "AliRsnCutMgr.h"
#include "AliRsnPairDef.h"
#include "AliRsnHistoDef.h"
#include "AliRsnPairParticle.h"

class TH1D;
class TH2D;
class AliRsnEvent;

class AliRsnSimpleFunction : public TNamed
{

public:

    AliRsnSimpleFunction
      (const char *name="", 
       AliRsnPairDef *pd=0, AliRsnHistoDef *hd=0, AliRsnCutMgr *cuts=0, 
       Option_t *option="");
    AliRsnSimpleFunction(const AliRsnSimpleFunction &copy);
    const AliRsnSimpleFunction& operator=(const AliRsnSimpleFunction &copy);
    virtual ~AliRsnSimpleFunction() {Clear();}
    virtual void Clear(Option_t *option = "");

    Bool_t              TrueFlag() {return fTrueFlag;}
    Bool_t              MixFlag() {return fMixFlag;}
    AliRsnPairParticle* GetPair() {return &fPair;}
    AliRsnPairDef*      GetPairDef() {return fPairDef;}
    AliRsnHistoDef*     GetHistoDef() {return fHistoDef;}
    AliRsnCutMgr*       GetCutMgr() {return fCuts;}
    TH1D*               GetHistogram1D() {return fHisto1D;}
    TH2D*               GetHistogram2D() {return fHisto2D;}

    void   SetTrueFlag(Bool_t value = kTRUE) {fTrueFlag = value;}
    void   SetMixFlag(Bool_t value = kTRUE) {fMixFlag = value;}
    void   SetPairDef(AliRsnPairDef *def) {fPairDef = def;}
    void   SetHistoDef(AliRsnHistoDef *def) {fHistoDef = def;}
    void   SetCutMgr(AliRsnCutMgr *cutMgr) {fCuts = cutMgr;}

    // virtual working routines
    virtual Bool_t Init();
    virtual Bool_t ProcessOne(AliRsnEvent *event);
    virtual Bool_t ProcessTwo(AliRsnEvent *event1, AliRsnEvent *event2);

protected:

    Bool_t CutPass(AliRsnDaughter *d);
    Bool_t CutPass(AliRsnPairParticle *p);
    Bool_t CutPass(AliRsnEvent *e);

    Bool_t               fTrueFlag;      // flag to store only true pairs
    Bool_t               fMixFlag;       // flag to reserve this object for mixing

    AliRsnPairParticle   fPair;          // utility class for pair
    AliRsnPairDef       *fPairDef;       // definitions for pair
    AliRsnHistoDef      *fHistoDef;      // definitions for histogram
    AliRsnCutMgr        *fCuts;          // cut manager
    
    TH1D                *fHisto1D;       // output histogram (1D)
    TH2D                *fHisto2D;       // output histogram (2D)

    // ROOT dictionary
    ClassDef(AliRsnSimpleFunction, 1)
};

#endif
