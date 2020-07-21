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
                                TString name = "name",
                                TString   photonCutNumberV0Reader       = "00200009227300008250400000"       
                          ) {


    
    TString fileName = AliAnalysisManager::GetCommonFileName();
    
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
        return 0x0;
    }
    
    // ================== GetInputEventHandler =============================
    AliVEventHandler *inputHandler=mgr->GetInputEventHandler();

    //=========  Set Cutnumber for V0Reader ================================
    TString cutnumberPhoton = photonCutNumberV0Reader.Data();
    TString cutnumberEvent = "00000003";
    AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();  // probably dont need this

    //========= Add V0 Reader to  ANALYSIS manager if not yet existent =====
    TString V0ReaderName        = Form("V0ReaderV1_%s_%s",cutnumberEvent.Data(),cutnumberPhoton.Data());
    //TString V0ReaderName("V0ReaderV1");    
    AliV0ReaderV1 *fV0ReaderV1  =  NULL;
    cout << "V0ReaderName " << V0ReaderName << endl;
    if( !(AliV0ReaderV1*)mgr->GetTask(V0ReaderName.Data()) ){
        cout << "V0Reader: " << V0ReaderName.Data() << " not found!!"<< endl;
        return 0x0;
    } else {
        cout << "V0Reader: " << V0ReaderName.Data() << " found!!"<< endl;
    }

    AliAnalysisTRDEfficiency* task = new AliAnalysisTRDEfficiency(name.Data());   
    if(!task) return 0x0;
    task->SelectCollisionCandidates(AliVEvent::kINT7);  // was kAnyINT

    task->SetV0ReaderName(V0ReaderName);

    // add your task to the manager
    mgr->AddTask(task);
    // your task needs input: here we connect the manager to your task
    mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());
    // same for the output
    mgr->ConnectOutput(task,1,mgr->CreateContainer("MyOutputContainer", TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));

    return task;
}
