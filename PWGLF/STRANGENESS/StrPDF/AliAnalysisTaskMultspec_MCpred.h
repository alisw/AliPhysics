#ifndef AliAnalysisTaskMultspec_MCpred_H
#define AliAnalysisTaskMultspec_MCpred_H

#include "TString.h"
#include "AliAnalysisTaskSE.h"
#include "THistManager.h"


class AliAnalysisTaskMultspec_MCpred : public AliAnalysisTaskSE {

  public:
    AliAnalysisTaskMultspec_MCpred();
    AliAnalysisTaskMultspec_MCpred(const char *name, TString lExtraOptions = "");
    virtual ~AliAnalysisTaskMultspec_MCpred();

    virtual void UserCreateOutputObjects();
    virtual void UserExec(Option_t *option);
    virtual void Terminate(Option_t *);

  private:

    //outputs
    THistManager *fHistos_misc;                               //!<! Output histos

    AliAnalysisTaskMultspec_MCpred(const AliAnalysisTaskMultspec_MCpred&);            // not implemented
    AliAnalysisTaskMultspec_MCpred& operator=(const AliAnalysisTaskMultspec_MCpred&); // not implemented

    ClassDef(AliAnalysisTaskMultspec_MCpred, 1);
    //version 1
};

#endif
