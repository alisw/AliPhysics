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
    AliRsnAnalysisME(const char *name = "AliRsnAnalysisME");
    AliRsnAnalysisME(const AliRsnAnalysisME& copy);
    virtual ~AliRsnAnalysisME() {;};

    // Implement this
    virtual void    RsnUserCreateOutputObjects();
    virtual void    RsnUserExec(Option_t*);
    virtual void    RsnTerminate(Option_t*);

    AliRsnAnalysisManager *GetAnalysisManager(TString name="");
    void SetAnalysisManagerName(const char*name) { fRsnAnalysisManager.SetName(name);};

    // Prior probs
    void            SetPriorProbability(AliPID::EParticleType type, Double_t p);
    void            DumpPriors();
    void            GetPriorProbability(Double_t *out)const;

  private:

    AliRsnAnalysisME& operator=(const AliRsnAnalysisME& /*copy*/) {return *this;}

    AliRsnAnalysisManager fRsnAnalysisManager;      // analysis main engine
    AliRsnPIDIndex        fPIDIndex;                // utility --> PID sorter
    AliRsnPIDIndex        fPIDIndexMix;             // utility --> PID sorter (mixed event)
    AliRsnEvent           fEvent;                   // utility --> event interface
    AliRsnEvent           fEventMix;                // utility --> event interface (mixed event)

    Double_t              fPrior[AliPID::kSPECIES]; // prior probabilities

    ClassDef(AliRsnAnalysisME, 1)
};

#endif
