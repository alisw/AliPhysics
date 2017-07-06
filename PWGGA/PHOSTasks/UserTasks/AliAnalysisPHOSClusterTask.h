#ifndef AliAnalysisPHOSClusterTask_H
#define AliAnalysisPHOSClusterTask_H

#include <Rtypes.h>
#include <TAxis.h>
#include <TRefArray.h>

#include <vector>

#include "AliAnalysisTaskSE.h"
#include "AliESDtrackCuts.h"
#include "AliEventCuts.h"


class TList;
class TTree;
class TParticle;
class TH1F;
class TH2F;

using std::vector;

class AliAnalysisPHOSClusterTask : public AliAnalysisTaskSE {
  public:
    AliAnalysisPHOSClusterTask();
    AliAnalysisPHOSClusterTask(const char *name);
    virtual ~AliAnalysisPHOSClusterTask();

    virtual void   UserCreateOutputObjects();
    virtual void   UserExec(Option_t *option);
    virtual void   Terminate(const Option_t*) {}


//    void          SetZVertexCut(Double_t Z)         {fZVertex = Z;} // Set z vertex cut parameter


///  private:
    protected:
    // Private methods

    TList*       fOutput; //!
    TH1F*        fTEST;   //!




    private:
    AliAnalysisPHOSClusterTask(const AliAnalysisPHOSClusterTask&);            //! Not implemented
    AliAnalysisPHOSClusterTask& operator=(const AliAnalysisPHOSClusterTask&); //! Not implemented
    ClassDef(AliAnalysisPHOSClusterTask,2);
};

#endif






