//
// Class AliRsnAnalysisPhiKK
//
// Virtual Class derivated from AliRsnVAnalysisTaskSE which will be base class
// for all RSN SE tasks
//
// authors: Martin Vala (martin.vala@cern.ch)
//          Alberto Pulvirenti (alberto.pulvirenti@ct.infn.it)
//
#ifndef ALIRSNANALYSISPHIKK_H
#define ALIRSNANALYSISPHIKK_H

#include "TClonesArray.h"

#include "AliRsnVAnalysisTaskSE.h"

#include "AliRsnDaughter.h"
#include "AliRsnMother.h"
#include "AliRsnPairDef.h"
#include "AliRsnEvent.h"
#include "AliRsnCutSet.h"

class AliRsnPIDDefESD;

class AliRsnAnalysisPhiKK : public AliRsnVAnalysisTaskSE
{
  public:
  
    AliRsnAnalysisPhiKK(const char *name = "RSNphi", Bool_t useKine = kFALSE);
    AliRsnAnalysisPhiKK(const AliRsnAnalysisPhiKK& copy);
    AliRsnAnalysisPhiKK& operator=(const AliRsnAnalysisPhiKK& copy);
    virtual ~AliRsnAnalysisPhiKK() {}

    virtual void            RsnUserCreateOutputObjects();
    virtual void            RsnUserExec(Option_t*);
    virtual void            RsnTerminate(Option_t*);
    virtual Bool_t          EventProcess();
    
    AliRsnCutSet*           GetEventCuts()          {return &fCutEvent;}
    AliRsnCutSet*           GetCommonDaughterCuts() {return &fCutTrackCommon;}
    AliRsnCutSet*           GetPosDaughterCuts()    {return &fCutTrackPos;}
    AliRsnCutSet*           GetNegDaughterCuts()    {return &fCutTrackNeg;}
    AliRsnCutSet*           GetMotherCuts()         {return &fCutPair;}
    AliRsnPairDef*          GetPairDef()            {return &fPairDef;}
    
    static Bool_t           IsTruePair(AliRsnDaughter *d1, AliRsnDaughter *d2);
    void                    AddFunction(AliRsnFunction* const fcn);
    TList*                  GenerateHistograms();

  private:
  
    AliRsnPairDef           fPairDef;           // resonance decay tree (fixed)
    
    AliRsnCutSet            fCutEvent;          // cut set for events
    AliRsnCutSet            fCutTrackCommon;    // cut set for tracks (common)
    AliRsnCutSet            fCutTrackPos;       // cut set for tracks (only pos)
    AliRsnCutSet            fCutTrackNeg;       // cut set for tracks (only neg)
    AliRsnCutSet            fCutPair;           // cut set for pairs

    TClonesArray            fFuncPM;            // collection of functions for unlike-sign
    TClonesArray            fFuncPP;            // collection of functions for like-sign ++
    TClonesArray            fFuncMM;            // collection of functions for like-sign --
    TClonesArray            fFuncTrue;          // collection of functions for unlike-sign true pairs
    
    TList                  *fOutList;           // list of output events

    ClassDef(AliRsnAnalysisPhiKK, 1)
};

#endif
