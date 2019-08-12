/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: Friederike Bock, Lucia Leardini , A. Marin                     *
 * Version 1.0                                                            *
 *                                                                        *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

//***************************************************************************************
//CutHandler contains all cuts for a certain analysis and trainconfig,
//it automatically checks length of cutStrings and takes care of the number of added cuts,
//no specification of the variable 'numberOfCuts' needed anymore.
//***************************************************************************************

class CutHandlerConvMaterial{
  public:
    CutHandlerConvMaterial(Int_t nMax=10){
      nCuts=0; nMaxCuts=nMax; validCuts = true;
      eventCutArray = new TString[nMaxCuts]; photonCutArray = new TString[nMaxCuts]; mesonCutArray = new TString[nMaxCuts]; clusterCutArray = new TString[nMaxCuts];
      for(Int_t i=0; i<nMaxCuts; i++) {eventCutArray[i] = ""; photonCutArray[i] = ""; mesonCutArray[i] = ""; clusterCutArray[i] = "";}
    }
    void AddCut(TString eventCut, TString photonCut){
      if(nCuts>=nMaxCuts) {cout << "ERROR in CutHandlerConvMaterial: Exceeded maximum number of cuts!" << endl; validCuts = false; return;}
      if( eventCut.Length()!=8 || photonCut.Length()!=26  ) {cout << "ERROR in CutHandlerConvMaterial: Incorrect length of cut string!" << endl; validCuts = false; return;}
      eventCutArray[nCuts]=eventCut; photonCutArray[nCuts]=photonCut; mesonCutArray[nCuts]=""; clusterCutArray[nCuts]="";
      nCuts++;
      return;
    }
    void AddCut(TString eventCut, TString photonCut, TString mesonCut){
      if(nCuts>=nMaxCuts) {cout << "ERROR in CutHandlerConvMaterial: Exceeded maximum number of cuts!" << endl; validCuts = false; return;}
      if( eventCut.Length()!=8 || photonCut.Length()!=26 || mesonCut.Length()!=16 ) {cout << "ERROR in CutHandlerConvMaterial: Incorrect length of cut string!" << endl; validCuts = false; return;}
      eventCutArray[nCuts]=eventCut; photonCutArray[nCuts]=photonCut; mesonCutArray[nCuts]=mesonCut; clusterCutArray[nCuts]="";
      nCuts++;
      return;
    }
    void AddCut(TString eventCut, TString photonCut, TString mesonCut, TString clusterCut){
      if(nCuts>=nMaxCuts) {cout << "ERROR in CutHandlerConvMaterial: Exceeded maximum number of cuts!" << endl; validCuts = false; return;}
      if( eventCut.Length()!=8 || photonCut.Length()!=26 || mesonCut.Length()!=16 || clusterCut.Length()!=19 ) {cout << "ERROR in CutHandlerConvMaterial: Incorrect length of cut string!" << endl; validCuts = false; return;}
      eventCutArray[nCuts]=eventCut; photonCutArray[nCuts]=photonCut; mesonCutArray[nCuts]=mesonCut; clusterCutArray[nCuts]=clusterCut;
      nCuts++;
      return;
    }

    Bool_t AreValid(){return validCuts;}
    Int_t GetNCuts(){if(validCuts) return nCuts; else return 0;}
    TString GetEventCut(Int_t i){if(validCuts&&i<nMaxCuts&&i>=0) return eventCutArray[i]; else{cout << "ERROR in CutHandlerConvMaterial: GetEventCut wrong index i" << endl;return "";}}
    TString GetPhotonCut(Int_t i){if(validCuts&&i<nMaxCuts&&i>=0) return photonCutArray[i]; else {cout << "ERROR in CutHandlerConvMaterial: GetPhotonCut wrong index i" << endl;return "";}}
    TString GetMesonCut(Int_t i){if(validCuts&&i<nMaxCuts&&i>=0) return mesonCutArray[i]; else {cout << "ERROR in CutHandlerConvMaterial: GetMesonCut wrong index i" << endl;return "";}}
    TString GetClusterCut(Int_t i){if(validCuts&&i<nMaxCuts&&i>=0) return clusterCutArray[i]; else {cout << "ERROR in CutHandlerConvMaterial: GetClusterCut wrong index i" << endl;return "";}}
  private:
    Bool_t validCuts;
    Int_t nCuts; Int_t nMaxCuts;
    TString* eventCutArray;
    TString* photonCutArray;
    TString* mesonCutArray;
    TString* clusterCutArray;
};




//***************************************************************************************
//main function
//***************************************************************************************
void AddTask_MaterialHistos_PbPb(	Int_t   trainConfig             = 1,        // change different set of cuts
                                    Int_t 	isMC 				    = 0,
                                    TString periodname              = "LHC10b", // period name
                                    TString periodNameV0Reader      = "",
                                    TString periodNameAnchor        = "",
                                    Bool_t 	enableV0EffiStudies 	= kFALSE,    // enables V0finding efficiency histograms
                                    TString additionalTrainConfig   = "0"       // additional counter for trainconfig, this has to be always the last parameter
                                ){

    Int_t IsHeavyIon = 1;
    if (additionalTrainConfig.Atoi() > 0){
        trainConfig = trainConfig + additionalTrainConfig.Atoi();
    }

    // ================== GetAnalysisManager ===============================
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
        Error(Form("AddTask_GammaConvV1_%i",trainConfig), "No analysis manager found.");
        return ;
    }

    // ================== GetInputEventHandler =============================
    AliVEventHandler *inputHandler=mgr->GetInputEventHandler();
    Bool_t isMCForOtherSettings = 0;
    if (isMC > 0) isMCForOtherSettings = 1;
    //========= Add PID Reponse to ANALYSIS manager ====
    if(!(AliPIDResponse*)mgr->GetTask("PIDResponseTask")){
        gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
        AddTaskPIDResponse(isMCForOtherSettings);
    }

    //=========  Set Cutnumber for V0Reader ================================
    TString cutnumberEvent = "10000003";
    TString cutnumberPhoton = "00000000000000000500000000";
    if (trainConfig>100 && trainConfig<200 ){
        cutnumberPhoton = "10000000000000000500000000";
    }


    Bool_t enableV0findingEffi = kFALSE;
    if(enableV0EffiStudies > 0){
        enableV0findingEffi = kTRUE;
        cutnumberPhoton = "00000070000000000500004000";
        if(enableV0EffiStudies == 1) cutnumberEvent = "50100013";
        else if(enableV0EffiStudies == 2) cutnumberEvent = "52500013";
    }


    Bool_t  enableConstructGamma = kFALSE;

    AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
    //========= Add V0 Reader to  ANALYSIS manager if not yet existent =====
    TString V0ReaderName = Form("V0ReaderV1_%s_%s",cutnumberEvent.Data(),cutnumberPhoton.Data());
    if( !(AliV0ReaderV1*)mgr->GetTask(V0ReaderName.Data()) ){
        AliV0ReaderV1 *fV0ReaderV1 = new AliV0ReaderV1(V0ReaderName.Data());
        if (periodNameV0Reader.CompareTo("") != 0) fV0ReaderV1->SetPeriodName(periodNameV0Reader);
        fV0ReaderV1->SetUseOwnXYZCalculation(kTRUE);
        fV0ReaderV1->SetUseConstructGamma(enableConstructGamma);
        fV0ReaderV1->SetCreateAODs(kFALSE);// AOD Output
        fV0ReaderV1->SetUseAODConversionPhoton(kTRUE);
        fV0ReaderV1->SetProduceV0FindingEfficiency(enableV0findingEffi);
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
            fCuts->SetIsHeavyIon(IsHeavyIon);
            fCuts->SetV0ReaderName(V0ReaderName);
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

    } else {
        Error("AddTask_V0ReaderV1", "Cannot execute AddTask, V0ReaderV1 already exists.");
    }

  //================================================
  //========= Add task to the ANALYSIS manager =====
  //            find input container

  AliAnalysisTaskMaterialHistos *fMaterialHistos=NULL;
  fMaterialHistos= new AliAnalysisTaskMaterialHistos(Form("MaterialHistos_%i",trainConfig));
  fMaterialHistos->SetIsMC(isMC);
  fMaterialHistos->SetIsHeavyIon(IsHeavyIon);
  fMaterialHistos->SetV0ReaderName(V0ReaderName);

  CutHandlerConvMaterial cuts;
  if( trainConfig == 1){ // kMB, V0M, std cuts
    cuts.AddCut("50000013", "00000009247602008250404000");
  } else  if( trainConfig == 2){ // kMB, V0M, open cuts
    cuts.AddCut("50000013", "00000070000000000500004000");
  } else if ( trainConfig == 3){ // selection for 11h meson cut
    cuts.AddCut("50000013", "00000009247602008850404000");
  } else if ( trainConfig == 4){ // selection for 11h photon cut
    cuts.AddCut("50000013", "00000009247002008850404000");
  } else if ( trainConfig == 5){ // selection for 11h meson cut
    cuts.AddCut("60100013", "00000009247602008850404000");
    cuts.AddCut("61200013", "00000009247602008850404000");
    cuts.AddCut("50100013", "00000009247602008850404000");
    cuts.AddCut("51200013", "00000009247602008850404000");
    cuts.AddCut("52500013", "00000009247602008850404000");
    cuts.AddCut("55800013", "00000009247602008850404000");
  } else if ( trainConfig == 6){ // selection for 11h photon cut
    cuts.AddCut("60100013", "00000009247002008850404000");
    cuts.AddCut("61200013", "00000009247002008850404000");
    cuts.AddCut("50100013", "00000009247002008850404000");
    cuts.AddCut("51200013", "00000009247002008850404000");
    cuts.AddCut("52500013", "00000009247002008850404000");
    cuts.AddCut("55800013", "00000009247002008850404000");
    // Offline V0Finder is used
  } else  if(trainConfig == 101){ //INT7, CL1
    cuts.AddCut("20010103", "10000009247602008250404000");
  } else  if( trainConfig == 102){ // kMB, V0M, open cuts
    cuts.AddCut("50000013", "10000070000000000500004000");
  } else if ( trainConfig == 103){ // selection for 11h meson cut
    cuts.AddCut("50000013", "10000009247602008850404000");
  } else if ( trainConfig == 104){ // selection for 11h photon cut
    cuts.AddCut("50000013", "10000009247002008850404000");
  } else if ( trainConfig == 105){ // selection for 11h meson cut
    cuts.AddCut("60100013", "10000009247602008850404000");
    cuts.AddCut("61200013", "10000009247602008850404000");
    cuts.AddCut("50100013", "10000009247602008850404000");
    cuts.AddCut("51200013", "10000009247602008850404000");
    cuts.AddCut("52500013", "10000009247602008850404000");
    cuts.AddCut("55800013", "10000009247602008850404000");
  } else if ( trainConfig == 106){ // selection for 11h photon cut
    cuts.AddCut("60100013", "10000009247002008850404000");
    cuts.AddCut("61200013", "10000009247002008850404000");
    cuts.AddCut("50100013", "10000009247002008850404000");
    cuts.AddCut("51200013", "10000009247002008850404000");
    cuts.AddCut("52500013", "10000009247002008850404000");
    cuts.AddCut("55800013", "10000009247002008850404000");

  } else {
    Error(Form("GammaConvV1_%i",trainConfig), "wrong trainConfig variable no cuts have been specified for the configuration");
    return;
  }

	if(!cuts.AreValid()){
    cout << "\n\n****************************************************" << endl;
    cout << "ERROR: No valid cuts stored in CutHandlerConvMaterial! Returning..." << endl;
    cout << "****************************************************\n\n" << endl;
    return;
  }

  Int_t numberOfCuts = cuts.GetNCuts();

  TList *EventCutList = new TList();
  TList *ConvCutList = new TList();


  EventCutList->SetOwner(kTRUE);
  AliConvEventCuts **analysisEventCuts        = new AliConvEventCuts*[numberOfCuts];
  ConvCutList->SetOwner(kTRUE);
  AliConversionPhotonCuts **analysisCuts      = new AliConversionPhotonCuts*[numberOfCuts];
  for(Int_t i = 0; i<numberOfCuts; i++){
    analysisEventCuts[i]          = new AliConvEventCuts();
    analysisEventCuts[i]->SetV0ReaderName(V0ReaderName);
    analysisEventCuts[i]->InitializeCutsFromCutString((cuts.GetEventCut(i)).Data());
    EventCutList->Add(analysisEventCuts[i]);
    analysisEventCuts[i]->SetFillCutHistograms("",kTRUE);

    analysisCuts[i]               = new AliConversionPhotonCuts();
    analysisCuts[i]->SetIsHeavyIon(IsHeavyIon);
    analysisCuts[i]->SetV0ReaderName(V0ReaderName);
    analysisCuts[i]->InitializeCutsFromCutString((cuts.GetPhotonCut(i)).Data());
    ConvCutList->Add(analysisCuts[i]);
    analysisCuts[i]->SetFillCutHistograms("",kTRUE);

   }




  fMaterialHistos->SetEventCutList(numberOfCuts,EventCutList);
  fMaterialHistos->SetConversionCutList(numberOfCuts,ConvCutList);

  AliAnalysisDataContainer *coutput =
    mgr->CreateContainer(Form("GammaConvMaterial_%i",trainConfig), TList::Class(),
			 AliAnalysisManager::kOutputContainer,Form("GammaConv_Material_%i.root",trainConfig ));

  mgr->AddTask(fMaterialHistos);
  mgr->ConnectInput(fMaterialHistos,  0, cinput );
  mgr->ConnectOutput(fMaterialHistos,  1, coutput);
  //connect containers
	return;
}
