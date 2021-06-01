/*!
   \file AliMCSpectraWeightsAnalysisTask.cxx
   \brief A minimal analysis task for AliMCSpectraWeights
   \author Patrick Huhn
   \date 25/10/2019
*/
#ifndef ALIMCSPECTRAWEIGHTSANALYSISTASK_H
#define ALIMCSPECTRAWEIGHTSANALYSISTASK_H
#include "AliAnalysisTaskSE.h"
#include "TH3F.h"
#include "THn.h"
class TArrayD;
class TString;
class AliMCSpectraWeights;
class AliStack;
class AliMCEvent;
class AliVEvent;
class TList;
class AliEventCuts;

/**
 * \class AliMCSpectraWeightsAnalysisTask
 * \brief Minimal analysis task for AliMCSpectraWeights
 */
class AliMCSpectraWeightsAnalysisTask : public AliAnalysisTaskSE {
  public:
    AliMCSpectraWeightsAnalysisTask();
    AliMCSpectraWeightsAnalysisTask(const char* name);
    virtual ~AliMCSpectraWeightsAnalysisTask();

    virtual void UserCreateOutputObjects();
    virtual void UserExec(Option_t* option);
    virtual void Terminate(Option_t* option) {}

    //------ Setter ---------
    void SetDebugLevel(Int_t level) { fDebugLevel = level; }
    void SetMCSpectraWeightObject(AliMCSpectraWeights* obj) {
        fMCSpectraWeights = obj;
    }
    void SetTriggerMask(UInt_t triggermask) { fTriggerMask = triggermask; }
    /// Set the flag for the use of MC.
    void SetUseMC(Bool_t useMC = kTRUE) { fIsMC = useMC; }
    
    //------ Gettter --------
    UInt_t GetTriggerMask() { return fTriggerMask; }
    
   //AddTask
    static AliMCSpectraWeightsAnalysisTask* AddTaskWeights(const char* collisionSystem,
    const char* previousTrain, const char* name = "MCWeightsAnalysisTask", const char* outfile = 0);
    
  private:
    Int_t fDebugLevel;       ///!< Debug level
    TList* fOutputList;      //!<! Output list
    AliVEvent* fEvent;       //!<! Event object (AliVEvent)
    AliMCEvent* fMCEvent;    //!<! MC event
    AliStack* fMCStack;      //!<! MC stack
    AliEventCuts fEventCuts; //!<! event cuts
    Bool_t fIsMC;            //
    UInt_t fTriggerMask;     // trigger mask

    // Particle composition
    AliMCSpectraWeights* fMCSpectraWeights; //->

    AliMCSpectraWeightsAnalysisTask(const AliMCSpectraWeightsAnalysisTask&);
    AliMCSpectraWeightsAnalysisTask&
    operator=(const AliMCSpectraWeightsAnalysisTask&);

    ClassDef(AliMCSpectraWeightsAnalysisTask, 2);
};

#endif /* ALIMCSPECTRAWEIGHTSANALYSISTASK_H */
