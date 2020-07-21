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
    //cout << "$$$$$$$$$ part of the task $$$$$$$$" << endl;
    //=========  Set Cutnumber for V0Reader ================================
    TString cutnumberPhoton = "00200009227300008250400000";
    TString cutnumberEvent = "00010113";
    AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
    //========= Add V0 Reader to  ANALYSIS manager if not yet existent =====
    TString V0ReaderName = Form("V0ReaderV1_%s_%s",cutnumberEvent.Data(),cutnumberPhoton.Data());
    if( !(AliV0ReaderV1*)mgr->GetTask(V0ReaderName.Data()) ){
        AliV0ReaderV1 *fV0ReaderV1 = new AliV0ReaderV1(V0ReaderName.Data());
        //cout << "@@@@@@@@@  in the movement  @@@@@@@@@" << endl;
        fV0ReaderV1->SetUseOwnXYZCalculation(kTRUE);
        fV0ReaderV1->SetCreateAODs(kFALSE);// AOD Output
        fV0ReaderV1->SetUseAODConversionPhoton(kTRUE);
        fV0ReaderV1->SetUseMassToZero(kFALSE);
        if (!mgr) {
            Error("AddTask_V0ReaderV1", "No analysis manager found.");
            return 0x0;
        }

        AliConvEventCuts *fEventCuts=NULL;
        if(cutnumberEvent!=""){
            fEventCuts= new AliConvEventCuts(cutnumberEvent.Data(),cutnumberEvent.Data());
            fEventCuts->SetPreSelectionCutFlag(kTRUE);
            fEventCuts->SetV0ReaderName(V0ReaderName);
            if(fEventCuts->InitializeCutsFromCutString(cutnumberEvent.Data())){
                fV0ReaderV1->SetEventCuts(fEventCuts);
                fEventCuts->SetFillCutHistograms("",kTRUE);
            }
        }

        // Set AnalysisCut Number
        AliConversionPhotonCuts *fCuts=NULL;
        if(cutnumberPhoton!=""){
            fCuts= new AliConversionPhotonCuts(cutnumberPhoton.Data(),cutnumberPhoton.Data());
            fCuts->SetPreSelectionCutFlag(kTRUE);
            fCuts->SetIsHeavyIon(0);
            fCuts->SetV0ReaderName(V0ReaderName);
            if(fCuts->InitializeCutsFromCutString(cutnumberPhoton.Data())){
                fV0ReaderV1->SetConversionCuts(fCuts);
                fCuts->SetFillCutHistograms("",kTRUE);
            }
        }
        //TString cutnumberAODBranch = Form("GammaConv_%s_%s_gamma", cutnumberEvent.Data(), cutnumberPhoton.Data());
        //if(inputHandler->IsA()==AliAODInputHandler::Class()){
            // AOD mode
        //    fV0ReaderV1->AliV0ReaderV1::SetDeltaAODBranchName(Form("GammaConv_%s_gamma",cutnumberAODBranch.Data()));
        //}
        fV0ReaderV1->Init();

        //AliLog::SetGlobalLogLevel(AliLog::kInfo);

        //connect input V0Reader
        mgr->AddTask(fV0ReaderV1);
        mgr->ConnectInput(fV0ReaderV1,0,cinput);
    }
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
