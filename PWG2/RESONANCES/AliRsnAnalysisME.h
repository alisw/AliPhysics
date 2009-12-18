//
// Class AliRsnAnalysisME
//
// Virtual Class derivated from AliRsnVAnalysisTaskME which will be base class
// for all RSN SE tasks
//
// authors: Martin Vala (martin.vala@cern.ch)
//          Alberto Pulvirenti (alberto.pulvirenti@ct.infn.it)
//

#ifndef ALIRSNANALYSISME_H
#define ALIRSNANALYSISME_H

#include "AliPID.h"
#include "AliRsnVAnalysisTaskME.h"
#include "AliRsnAnalysisManager.h"
#include "AliRsnEvent.h"
#include "AliRsnPIDIndex.h"

class AliPID;
class AliESDtrackCuts;
class AliRsnAnalysisME : public AliRsnVAnalysisTaskME
{

  public:
    AliRsnAnalysisME(const char *name = "AliRsnAnalysisME", Int_t numOfOutputs = 1);
    AliRsnAnalysisME(const AliRsnAnalysisME& copy);
    virtual ~AliRsnAnalysisME() { ; };

    // Implement this
    virtual void    RsnUserCreateOutputObjects();
    virtual void    RsnUserExec(Option_t*);
    virtual void    RsnTerminate(Option_t*);

    AliRsnAnalysisManager *GetAnalysisManager(Int_t index = 0, TString name = "");
    void SetAnalysisManagerName(const char *name, Int_t index = 0) { fRsnAnalysisManager[index].SetName(name); };

    // Prior probs
    void            SetPriorProbability(AliPID::EParticleType type, Double_t p);
    void            DumpPriors();
    void            GetPriorProbability(Double_t *out)const;

  private:

    AliRsnAnalysisME& operator=(const AliRsnAnalysisME& /*copy*/) { return *this; }

    AliRsnAnalysisManager fRsnAnalysisManager[kMaxNumberOfOutputs];      // analysis main engine
    AliRsnPIDIndex        fPIDIndex;                // utility --> PID sorter
    AliRsnPIDIndex        fPIDIndexMix;             // utility --> PID sorter (mixed event)
    AliRsnEvent           fEvent;                   // utility --> event interface
    AliRsnEvent           fEventMix;                // utility --> event interface (mixed event)

    Double_t              fPrior[AliPID::kSPECIES]; // prior probabilities

    void                  DoMixing(AliVEvent *ev);
    void                  DoAODMixing(AliAODEvent* aod1, AliAODEvent* aod2);
    void                  DoESDMixing(AliESDEvent* esd1, AliESDEvent* esd2);

    ClassDef(AliRsnAnalysisME, 1)
};

#endif
