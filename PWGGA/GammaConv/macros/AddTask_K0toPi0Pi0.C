/**************************************************************************************
 * Copyright (C) 2017, Copyright Holders of the ALICE Collaboration                   *
 * All rights reserved.                                                               *
 *                                                                                    *
 * Redistribution and use in source and binary forms, with or without                 *
 * modification, are permitted provided that the following conditions are met:        *
 *     * Redistributions of source code must retain the above copyright               *
 *       notice, this list of conditions and the following disclaimer.                *
 *     * Redistributions in binary form must reproduce the above copyright            *
 *       notice, this list of conditions and the following disclaimer in the          *
 *       documentation and/or other materials provided with the distribution.         *
 *     * Neither the name of the <organization> nor the                               *
 *       names of its contributors may be used to endorse or promote products         *
 *       derived from this software without specific prior written permission.        *
 *                                                                                    *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND    *
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED      *
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE             *
 * DISCLAIMED. IN NO EVENT SHALL ALICE COLLABORATION BE LIABLE FOR ANY                *
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES         *
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;       *
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND        *
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT         *
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS      *
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                       *
 **************************************************************************************/

void AddTask_K0toPi0Pi0(Bool_t runLightOutput = kFALSE, 
						TString periodName = "",
						TString periodNameV0Reader = "",
						TString selectConfig = "0") {

  Int_t trainConfig = selectConfig.Atoi();
  // ================== GetAnalysisManager ===============================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTask_K0toPi0Pi0", "No analysis manager found.");
    return ;
  }

  // ================== GetInputEventHandler =============================
  AliVEventHandler *inputHandler=mgr->GetInputEventHandler();
  
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
  
  
  //================================================
  //========= Add task to the ANALYSIS manager =====
  //================================================
  //            find input container
  AliAnalysisTaskK0toPi0Pi0 *task=NULL;
  task = new AliAnalysisTaskK0toPi0Pi0("TaskK0toPi0Pi0" + selectConfig);
  

   //=========  Set Cutnumber for V0Reader ================================
  //TString cutnumberPhoton = "00200008400000002200000000"; 
  TString cutnumberPhoton = "00000000000000000000000000"; 
  TString cutnumberEvent = "00000003"; 
 
//========= Add V0 Reader to  ANALYSIS manager if not yet existent =====
  TString V0ReaderName = Form("V0ReaderV1_%s_%s",cutnumberEvent.Data(),cutnumberPhoton.Data());
  if( !(AliV0ReaderV1*)mgr->GetTask(V0ReaderName.Data()) ){
    AliV0ReaderV1 *fV0ReaderV1 = new AliV0ReaderV1(V0ReaderName.Data());
    if (periodNameV0Reader.CompareTo("") != 0) fV0ReaderV1->SetPeriodName(periodNameV0Reader);
    fV0ReaderV1->SetUseOwnXYZCalculation(kTRUE);
    fV0ReaderV1->SetCreateAODs(kFALSE);// AOD Output
    fV0ReaderV1->SetUseAODConversionPhoton(kTRUE);

    if (!mgr) {
      Error("AddTask_V0ReaderV1", "No analysis manager found.");
      return;
    }

    AliConvEventCuts *fEventCuts=NULL;
    if(cutnumberEvent!=""){
      fEventCuts= new AliConvEventCuts(cutnumberEvent.Data(),cutnumberEvent.Data());
      fEventCuts->SetPreSelectionCutFlag(kTRUE);
      fEventCuts->SetV0ReaderName(V0ReaderName);
      if (periodNameV0Reader.CompareTo("") != 0) fEventCuts->SetPeriodEnum(periodNameV0Reader);
      fEventCuts->SetLightOutput(runLightOutput);
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
      //fCuts->SetIsHeavyIon(isHeavyIon);
      fCuts->SetV0ReaderName(V0ReaderName);
      fCuts->SetLightOutput(runLightOutput);
      if(fCuts->InitializeCutsFromCutString(cutnumberPhoton.Data())){
        fV0ReaderV1->SetConversionCuts(fCuts);
        fCuts->SetFillCutHistograms("",kTRUE);
      }
    }
    
   
    fV0ReaderV1->Init();

    AliLog::SetGlobalLogLevel(AliLog::kInfo);
    //connect input V0Reader
    mgr->AddTask(fV0ReaderV1);
    mgr->ConnectInput(fV0ReaderV1,0,cinput);

  }
  
  task->SetNameV0Reader(V0ReaderName); 
  
  TString defaultConvPhotonCut = "00200009327000008250400000";
  TString defaultCaloPhotonCut = "1111113067002230000"; // default cut string 1111111067032230000, currently energy cut set to 0
  TString defaultPi0Cut = "0163103100000010"; // conv/conv
  TString Pi0Cut = "0163103100000060"; // calo/calo
  TString mixedPi0Cut = "0163103100000060"; // conv/calo
  TString defaultK0Cut = "0163103100000000";
  
  
  
  if (trainConfig == 0){// ================================================= trainConfig = 0 -> 8 TeV MB
  	TString defaultEventCut = "00010113";
  	TString defaultConvPhotonCut = "00200009327000008250400000";
  	TString defaultCaloPhotonCut = "1111113067002230000"; // default cut string 1111111067032230000, currently energy cut set to 0
  	TString defaultPi0Cut = "0163103700000010"; // conv conv
  	TString Pi0Cut = "0163103700000060"; // calo/calo
  	TString mixedPi0Cut = "0163103700000060"; // conv/calo
  	TString defaultK0Cut = "0163103000000000";
  } else if (trainConfig == 1){ // ========================================= traingConfig = 1 -> 8 TeV EMC7 
  	TString defaultEventCut = "00052113";
  	TString defaultConvPhotonCut = "00200009327000008250400000";
  	TString defaultCaloPhotonCut = "1111113067002230000"; // default cut string 1111111067032230000, currently energy cut set to 0
  	TString defaultPi0Cut = "0163103700000010"; //conv conv
  	TString Pi0Cut = "0163103700000060"; // calo/calo
  	TString mixedPi0Cut = "0163103700000060"; // conv/calo
  	TString defaultK0Cut = "0163103000000000";
  } else if (trainConfig == 2){ // ========================================= traingConfig = 2 -> 8 TeV EGA 
  	TString defaultEventCut = "00081113";
  	TString defaultConvPhotonCut = "00200009327000008250400000";
  	TString defaultCaloPhotonCut = "1111113067002230000"; // default cut string 1111111067032230000, currently energy cut set to 0
  	TString defaultPi0Cut = "0163103700000010"; // conv conv
  	TString Pi0Cut = "0163103700000060"; // calo/calo
  	TString mixedPi0Cut = "0163103700000060"; // conv/calo
  	TString defaultK0Cut = "0163103000000000";
  } else if (trainConfig == 3){ // ========================================= traingConfig = 3 -> 7 TeV MB LHC 10 (V0OR)
  	TString defaultEventCut = "00000113";
  	TString defaultConvPhotonCut = "00200009327000008250400000";
  	TString defaultCaloPhotonCut = "1111113067002230000"; // default cut string 1111111067032230000, currently energy cut set to 0
  	TString defaultPi0Cut = "0163103700000010"; // conv conv
  	TString Pi0Cut = "0163103700000060"; // calo/calo
  	TString mixedPi0Cut = "0163103700000060"; // conv/calo
  	TString defaultK0Cut = "0163103000000000";
  }
  
  //create AliCaloTrackMatcher instance, if there is none present
  TString caloCutPos = defaultCaloPhotonCut;
  caloCutPos.Resize(1);
  TString TrackMatcherName = Form("CaloTrackMatcher_%s",caloCutPos.Data());
  if( !(AliCaloTrackMatcher*)mgr->GetTask(TrackMatcherName.Data()) ){
    AliCaloTrackMatcher* fTrackMatcher = new AliCaloTrackMatcher(TrackMatcherName.Data(),caloCutPos.Atoi());
    fTrackMatcher->SetV0ReaderName(V0ReaderName);
    mgr->AddTask(fTrackMatcher);
    mgr->ConnectInput(fTrackMatcher,0,cinput);
  }
  
  // ============= Set the different cut types to their default values ==============
  AliConvEventCuts *analysisEventCuts = new AliConvEventCuts(); 
  analysisEventCuts->SetPeriodEnum(periodName);
  analysisEventCuts->SetV0ReaderName(V0ReaderName);
  analysisEventCuts->InitializeCutsFromCutString(defaultEventCut.Data());
  analysisEventCuts->SetFillCutHistograms("", kTRUE); 
  task->SetEventCuts(analysisEventCuts); 

  AliConversionPhotonCuts *analysisConvoCuts = new AliConversionPhotonCuts();  
  analysisConvoCuts->SetV0ReaderName(V0ReaderName);
  analysisConvoCuts->InitializeCutsFromCutString(defaultConvPhotonCut.Data());
  analysisConvoCuts->SetFillCutHistograms("",kTRUE);
  task->SetConversionPhotonCuts(analysisConvoCuts); 
  
  AliCaloPhotonCuts *analysisCaloCuts = new AliCaloPhotonCuts(); 
  analysisCaloCuts->SetV0ReaderName(V0ReaderName);
  analysisCaloCuts->InitializeCutsFromCutString(defaultCaloPhotonCut.Data());
  analysisCaloCuts->SetFillCutHistograms("");
  task->SetCaloPhotonCuts(analysisCaloCuts); 
  
  
  AliConversionMesonCuts *analysisPi0CutsConvConv = new AliConversionMesonCuts();
  analysisPi0CutsConvConv->InitializeCutsFromCutString(defaultPi0Cut.Data()); 
  analysisPi0CutsConvConv->SetFillCutHistograms(""); 
  task->SetPi0CutsConvConv(analysisPi0CutsConvConv); 
  
  
  AliConversionMesonCuts *analysisPi0CutsCaloCalo = new AliConversionMesonCuts();
  analysisPi0CutsCaloCalo->InitializeCutsFromCutString(Pi0Cut.Data()); 
  analysisPi0CutsCaloCalo->SetFillCutHistograms(""); 
  task->SetPi0CutsCaloCalo(analysisPi0CutsCaloCalo);
  
  AliConversionMesonCuts *analysisPi0CutsConvCalo = new AliConversionMesonCuts();
  analysisPi0CutsConvCalo->InitializeCutsFromCutString(mixedPi0Cut.Data()); 
  analysisPi0CutsConvCalo->SetFillCutHistograms(""); 
  task->SetPi0CutsConvCalo(analysisPi0CutsConvCalo);  
  
  
  AliConversionMesonCuts *analysisK0Cuts = new AliConversionMesonCuts(); 
  analysisK0Cuts->InitializeCutsFromCutString(defaultK0Cut.Data()); 
  analysisK0Cuts->SetFillCutHistograms(""); 
  task->SetK0Cuts(analysisK0Cuts); 
  
  
  //connect containers
  AliAnalysisDataContainer *coutput =
  mgr->CreateContainer("K0toPi0Pi0" + selectConfig, TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:K0toPi0Pi0%s",AliAnalysisManager::GetCommonFileName(),selectConfig.Data()));
  mgr->AddTask(task);
  mgr->ConnectInput(task,0,cinput);
  mgr->ConnectOutput(task,1,coutput);
  
}
