#ifndef ALIFMDANALYSISTASKSE_H
#define ALIFMDANALYSISTASKSE_H

#include "AliAnalysisTaskSE.h"
#include "AliFMDAnalysisTaskSharing.h"
#include "AliFMDAnalysisTaskDensity.h"
#include "AliFMDAnalysisTaskBackgroundCorrection.h"
#include "AliFMDAnalysisTaskDndeta.h"

class AliFMDAnalysisTaskSE : public AliAnalysisTaskSE
{
 public:
    AliFMDAnalysisTaskSE();
    AliFMDAnalysisTaskSE(const char* name);
    virtual ~AliFMDAnalysisTaskSE() {;}
 AliFMDAnalysisTaskSE(const AliFMDAnalysisTaskSE& o) : AliAnalysisTaskSE(),
      fListOfHistos(o.fListOfHistos),
      fSharing(o.fSharing),
      fDensity(o.fDensity),
      fBackground(o.fBackground),
      fDndeta(o.fDndeta){}
    AliFMDAnalysisTaskSE& operator=(const AliFMDAnalysisTaskSE&) { return *this; }
    
    virtual void Init();
    virtual void UserCreateOutputObjects();
    virtual void UserExec(Option_t* /*option*/);
    void  Terminate(Option_t */*option*/);
 private:
    
    TList* fListOfHistos;
    AliFMDAnalysisTaskSharing              fSharing;
    AliFMDAnalysisTaskDensity              fDensity;
    AliFMDAnalysisTaskBackgroundCorrection fBackground;
    AliFMDAnalysisTaskDndeta               fDndeta;
    ClassDef(AliFMDAnalysisTaskSE, 1);

};
#endif
