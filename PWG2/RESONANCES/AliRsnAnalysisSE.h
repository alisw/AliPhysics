//
// Class AliRsnAnalysisSE
//
// Virtual Class derivated from AliRsnVAnalysisTaskSE which will be base class
// for all RSN SE tasks
//
// authors: Martin Vala (martin.vala@cern.ch)
//          Alberto Pulvirenti (alberto.pulvirenti@ct.infn.it)
//
#ifndef AliRsnAnalysisSE_H
#define AliRsnAnalysisSE_H

#include "AliPID.h"
#include "AliRsnVAnalysisTaskSE.h"
#include "AliRsnAnalysisManager.h"
#include "AliRsnPIDIndex.h"

class AliRsnAnalysisSE : public AliRsnVAnalysisTaskSE
{

public:
    AliRsnAnalysisSE(const char *name = "AliRsnAnalysisSE");
    AliRsnAnalysisSE(const AliRsnAnalysisSE& copy);

    // Implement this
    virtual void    RsnUserCreateOutputObjects();
    virtual void    RsnUserExec(Option_t*);
    virtual void    RsnTerminate(Option_t*);

    AliRsnAnalysisManager *GetAnalysisManager(TString name = "");
    void                   SetAnalysisManagerName(const char*name) {fRsnAnalysisManager.SetName(name);};

    // Prior probs
    void            SetPriorProbability(AliPID::EParticleType type, Double_t p);
    void            DumpPriors();
    void            GetPriorProbability(Double_t *out);

    // ESD cuts
    void            SetESDtrackCuts(AliESDtrackCuts *cuts) {fESDCuts = cuts;}

private:

    AliRsnAnalysisSE& operator=(const AliRsnAnalysisSE& /*copy*/) {return *this;}

    AliRsnAnalysisManager fRsnAnalysisManager;      // analysis main engine
    AliRsnPIDIndex        fPIDIndex;                // utility --> PID sorter
    AliRsnEvent           fEvent;                   // utility --> event interface

    AliESDtrackCuts      *fESDCuts;                 // ESD track cuts
    Double_t              fPrior[AliPID::kSPECIES]; // prior probabilities

    ClassDef(AliRsnAnalysisSE, 1)
};

#endif
