//
// Class AliRsnAnalysisKStarKpi
//
// Virtual Class derivated from AliRsnVAnalysisTaskSE which will be base class
// for all RSN SE tasks
//
// authors: Martin Vala (martin.vala@cern.ch)
//          Alberto Pulvirenti (alberto.pulvirenti@ct.infn.it)
//
#ifndef AliRsnAnalysisKStarKpi_H
#define AliRsnAnalysisKStarKpi_H

#include "TClonesArray.h"

#include "AliRsnVAnalysisTaskSE.h"

#include "AliRsnDaughter.h"
#include "AliRsnMother.h"
#include "AliRsnPairDef.h"
#include "AliRsnEvent.h"
#include "AliRsnCutSet.h"

class AliRsnPIDDefESD;

class AliRsnAnalysisKStarKpi : public AliRsnVAnalysisTaskSE
{
  public:
  
    AliRsnAnalysisKStarKpi(const char *name = "RSNphi", Bool_t useKine = kFALSE);
    AliRsnAnalysisKStarKpi(const AliRsnAnalysisKStarKpi& copy);
    AliRsnAnalysisKStarKpi& operator=(const AliRsnAnalysisKStarKpi& copy);
    virtual ~AliRsnAnalysisKStarKpi() {}

    virtual void            RsnUserCreateOutputObjects();
    virtual void            RsnUserExec(Option_t*);
    virtual void            RsnTerminate(Option_t*);
    virtual Bool_t          EventProcess();
    
    AliRsnCutSet*           GetEventCuts()          {return &fCutEvent;}
    AliRsnCutSet*           GetCommonDaughterCuts() {return &fCutTrackCommon;}
    AliRsnCutSet*           GetKaonDaughterCuts()   {return &fCutTrackKaon;}
    AliRsnCutSet*           GetPionDaughterCuts()   {return &fCutTrackPion;}
    AliRsnCutSet*           GetMotherCuts()         {return &fCutPair;}
    
    void                    AddFunction(AliRsnFunction* const fcn);
    TList*                  GenerateHistograms();

  private:
  
    TArrayI                 fGoodK;             // indexes of good tracks (kaons)
    TArrayI                 fGoodPi;            // indexes of good tracks (pions)
    AliRsnDaughter          fKaon;              // daughter data-member (pion)
    AliRsnDaughter          fPion;              // daughter data-member (pion)
    AliRsnMother            fMother;            // mother data-member (to save memory)
    AliRsnPairDef           fPairDef;           // resonance decay tree (fixed)
    
    AliRsnCutSet            fCutEvent;          // cut set for events
    AliRsnCutSet            fCutTrackCommon;    // cut set for tracks (common)
    AliRsnCutSet            fCutTrackKaon;      // cut set for tracks (only pos)
    AliRsnCutSet            fCutTrackPion;      // cut set for tracks (only neg)
    AliRsnCutSet            fCutPair;           // cut set for pairs

    TClonesArray            fFuncPM;            // collection of functions for unlike-sign (K+ pi-)
    TClonesArray            fFuncMP;            // collection of functions for unlike-sign (K- pi+)
    TClonesArray            fFuncPP;            // collection of functions for like-sign ++
    TClonesArray            fFuncMM;            // collection of functions for like-sign --
    TClonesArray            fFuncTruePM;        // collection of functions for unlike-sign true pairs (K+ pi-)
    TClonesArray            fFuncTrueMP;        // collection of functions for unlike-sign true pairs (K- pi+)
    
    TList                  *fOutList;           // list of output events

    ClassDef(AliRsnAnalysisKStarKpi, 1)
};

#endif
