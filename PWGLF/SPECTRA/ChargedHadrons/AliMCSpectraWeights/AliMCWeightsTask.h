//
//  AliMCWeightsTask.hpp
//  ParticleComposition
//
//  Created by Patrick Huhn on 05.10.20.
//  Copyright © 2020 Patrick Huhn. All rights reserved.
//

#ifndef AliMCWeightsTask_hpp
#define AliMCWeightsTask_hpp

#include <stdio.h>
#include <string>
#include "AliAnalysisTaskSE.h"
#include "TList.h"
#include "AliVEvent.h"
#include "AliMCEvent.h"
#include "AliMCSpectraWeights.h"

//#define __AliMCWeightsTask_DebugPCC__
//#define __AliMCWeightsTask_DebugTiming__

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
    AliMCWeightsTask(const AliMCWeightsTask&) = delete; // not implemented
    AliMCWeightsTask& operator=(const AliMCWeightsTask&) = delete; // not implemented

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



    /// \cond CLASSIMP
    ClassDef(AliMCWeightsTask, 2);
    /// \endcond
};


#endif /* AliMCWeightsTask_hpp */
