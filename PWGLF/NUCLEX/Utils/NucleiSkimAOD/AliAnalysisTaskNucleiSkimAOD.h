#ifndef AliAnalysisTaskNucleiSkimAOD_H
#define AliAnalysisTaskNucleiSkimAOD_H

#include <Rtypes.h>
#include <TAxis.h>

#include "AliAnalysisTaskSE.h"
#include "AliEventCuts.h"

class AliPIDResponse;
class TList;
class TParticle;

class AliAnalysisTaskNucleiSkimAOD : public AliAnalysisTaskSE {
  public:
    AliAnalysisTaskNucleiSkimAOD(const char *name = "AliAnalysisTaskNucleiSkimAOD");
    virtual ~AliAnalysisTaskNucleiSkimAOD();

    virtual void   UserCreateOutputObjects();
    virtual void   UserExec(Option_t *option);
    virtual void   Terminate(const Option_t*) {}

    AliEventCuts    mEventCuts;

  private:
    // Private methods
    AliAnalysisTaskNucleiSkimAOD(const AliAnalysisTaskNucleiSkimAOD&);            //! Not implemented
    AliAnalysisTaskNucleiSkimAOD& operator=(const AliAnalysisTaskNucleiSkimAOD&); //! Not implemented

    AliPIDResponse* mPIDresponse; //!<! PID response
    TList*          mOutputList;  //!<! List with the output histograms (event counters)


    ClassDef(AliAnalysisTaskNucleiSkimAOD,1)
};

#endif
