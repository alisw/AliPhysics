//
//  AliMCWeightsTask.hpp
//  ParticleComposition
//
//  Created by Patrick Huhn on 05.10.20.
//  Copyright © 2020 Patrick Huhn. All rights reserved.
//

#ifndef AliMCWeightsTask_hpp
#define AliMCWeightsTask_hpp

#include "AliAnalysisTaskSE.h"
#include "TList.h"
#include "AliVEvent.h"
#include "AliMCEvent.h"
#include "AliMCSpectraWeights.h"

enum MCGeneratorType {
    PP_PYTHIA=0,
    PPB_EPOS,
    PBPB_HIJING,
};

class AliMCWeightsTask : public AliAnalysisTaskSE {
public:
    AliMCWeightsTask();
    AliMCWeightsTask(const char* name);
    virtual ~AliMCWeightsTask();


    virtual void UserCreateOutputObjects();
    virtual void UserExec(Option_t* option);
    virtual void Terminate(Option_t* option) {}

    void SetMCSpectraWeightObject (AliMCSpectraWeights* weights) {fMCSpectraWeights=weights;}

    static AliMCWeightsTask* AddTaskAliMCWeightsTask ();

private:
    TList* fOutputList;      //!<! Output list
    AliVEvent* fEvent;       //!<! Event object (AliVEvent)
    AliMCEvent* fMCEvent;    //!<! MC event
    AliMCSpectraWeights*    fMCSpectraWeights; //-> object to determine efficiency scaling

    AliMCWeightsTask(const AliMCWeightsTask&); // not implemented
    AliMCWeightsTask& operator=(const AliMCWeightsTask&); // not implemented

    ClassDef(AliMCWeightsTask, 2);
};


#endif /* AliMCWeightsTask_hpp */
