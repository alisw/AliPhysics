//
// Class AliRsnEventTaskSE
//
// TODO
//
// authors: Martin Vala (martin.vala@cern.ch)
//          Alberto Pulvirenti (alberto.pulvirenti@ct.infn.it)
//
#ifndef ALIRSNEVENTTASK_H
#define ALIRSNEVENTTASK_H

#include <TClonesArray.h>
#include "AliRsnAnalysisTaskSEBase.h"

class TList;
class AliRsnEvent;
class AliRsnEventFunction;


class AliRsnEventTaskSE : public AliRsnAnalysisTaskSEBase
{
  public:
    AliRsnEventTaskSE(const char *name = "AliRsnEventTaskSE");
    AliRsnEventTaskSE(const AliRsnEventTaskSE& copy): AliRsnAnalysisTaskSEBase(copy),
        fEventFunctions("AliRsnEventFunction", 0),fOutList(0x0) {}
    AliRsnEventTaskSE& operator= (const AliRsnEventTaskSE&) {return *this;}
    ~AliRsnEventTaskSE();

    virtual void InitIOVars();
    virtual void UserCreateOutputObjects();
    virtual void UserExec(Option_t *option);
    virtual void Terminate(Option_t *);

    void AddEventFunction(AliRsnEventFunction *fcn);

  private:

    TClonesArray  fEventFunctions;
    TList        *fOutList;              // List of output

    void          ProcessEventAnalysis(AliRsnEvent *curEvent);
    void          PostEventProcess(const Short_t &index=0);

    ClassDef(AliRsnEventTaskSE, 1)
};

#endif
