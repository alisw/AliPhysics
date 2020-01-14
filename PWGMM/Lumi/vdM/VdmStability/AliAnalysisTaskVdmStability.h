#ifndef AliAnalysisTaskVdmStability_H
#define AliAnalysisTaskVdmStability_H

#include "AliAnalysisTaskSE.h"
#include "AliAODEvent.h"

class AliAnalysisTaskVdmStability : public AliAnalysisTaskSE
{
public:
    //Class constructors
    AliAnalysisTaskVdmStability();
    AliAnalysisTaskVdmStability(const char *taskname);
    //class destructor
    virtual         ~AliAnalysisTaskVdmStability();
    //define output - called once at beginning of runtime
    virtual void    UserCreateOutputObjects();
    //event loop - called for each event
    virtual void    UserExec(Option_t *);
    //terminate - called at the end of analysis
    virtual void    Terminate(Option_t *);

private:
    AliAODEvent*    fAOD;           //!<! Input AOD event
    TList*          fOutputList;     //!<! output list
    TH1F*           fHist;           //!<! dummy histogram



    /// \cond CLASSDEF
    ClassDef(AliAnalysisTaskVdmStability,1);
    /// \endcond};
};
#endif
