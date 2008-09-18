//
// Class AliRsnReaderTask
//
// TODO
//
// authors: Martin Vala (martin.vala@cern.ch)
//          Alberto Pulvirenti (alberto.pulvirenti@ct.infn.it)
//
#ifndef ALIRSNREADERTASK_H
#define ALIRSNREADERTASK_H

#include "AliRsnBaseAT.h"

class AliRsnPID;
class AliRsnReader;

class AliRsnReaderTask : public AliRsnBaseAT
{
  public:

    AliRsnReaderTask();
    AliRsnReaderTask(const char *name);
    virtual ~AliRsnReaderTask() {}

    // Implementation of interface methods
    virtual void    InitIOVars();
    virtual void    LocalInit();
    virtual Bool_t  Notify();
    virtual void    CreateOutputObjects();
    virtual void    Exec(Option_t *option);
    virtual void    Terminate(Option_t *);

  private:

    AliRsnReaderTask(const AliRsnReaderTask&);
    AliRsnReaderTask& operator= (const AliRsnReaderTask&);

    TTree*        fRsnTree;     // output tree

    ClassDef(AliRsnReaderTask, 0)    // implementation of RsnReader as AnalysisTask
};

#endif
