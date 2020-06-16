/****************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved.   *
 *                                                                          *
 * Author: Andre Barreiros                                                  *
 * Version 1.0                                                              *
 *                                                                          *
 *                                                                          *
 * Permission to use, copy, modify and distribute this software and its     *
 * documentation strictly for non-commercial purposes is hereby granted     *
 * without fee, provided that the above copyright notice appears in all     *
 * copies and that both the copyright notice and this permission notice     *
 * appear in the supporting documentation. The authors make no claims       *
 * about the suitability of this software for any purpose. It is            *
 * provided "as is" without express or implied warranty.                    *
 ****************************************************************************/

//#include "AliAnalysisTRDEfficiency.h"


//***************************************************************************************
//This AddTask is supposed to set up the main task
//($ALIPHYSICS/PWGGA/GammaConv/AliAnalysisTRDEfficiency.cxx) 
//***************************************************************************************


AliAnalysisTRDEfficiency* AddTask_TRDEfficiency(
                                TString name = "name"
                          ) {


    
    TString fileName = AliAnalysisManager::GetCommonFileName();
    
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
        return 0x0;
    }

    // ================== GetInputEventHandler =============================
    AliVEventHandler *inputHandler=mgr->GetInputEventHandler();
    
    // get the input event handler, again via a static method. 
    // this handler is part of the managing system and feeds events
    // to your task
    if (!mgr->GetInputEventHandler()) {
        return 0x0;
    }
    
    AliAnalysisTRDEfficiency* task = new AliAnalysisTRDEfficiency(name.Data());   
    if(!task) return 0x0;
    task->SelectCollisionCandidates(AliVEvent::kINT7);  // was kAnyINT
    
    
    // add your task to the manager
    mgr->AddTask(task);
    // your task needs input: here we connect the manager to your task
    mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());
    // same for the output
    mgr->ConnectOutput(task,1,mgr->CreateContainer("MyOutputContainer", TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
    
    return task;
}
