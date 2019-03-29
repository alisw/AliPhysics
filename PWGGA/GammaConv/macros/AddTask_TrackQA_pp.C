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
//This AddTask is supposed to set up the main task
//($ALIPHYSICS/PWGGA/GammaConv/AliAnalysisTaskTrackQA.cxx) for
//pp together with all supporting classes
//***************************************************************************************



void AddTask_TrackQA_pp( Int_t   trainConfig             = 1,                  // change different set of cuts
                                Int_t   isMC                    = 0,    
				TString   photonCutNumberV0Reader = "",       // 00000008400000000100000000 nom. B, 00000088400000000100000000 low B
				TString periodName              = "",
                                TString periodNameAnchor        = "",
				Bool_t 	enableV0findingEffi     = kFALSE,    // enables V0finding efficiency histograms
				TString additionalTrainConfig   = "0"       // additional counter for trainconfig, this has to be always the last parameter
                              ){

				//	TString fileNameInputForMultWeighing    = "/alice/cern.ch/user/a/amarin/multWeightFile/histosForMultWeighting_pp.root",            
  				//      TString fileNameElecDeDxPostCalibration = "dEdxCorrectionMap_Period_Pass.root",

  Int_t IsHeavyIon = 0;
  AliCutHandlerPCM cuts;
  //  TString fileNameGammaPtWeights      = cuts.GetSpecialFileNameFromString (fileNameExternalInputs, "FGPW:");

  TString addTaskName                 = "AddTask_TrackQA_pp";
  TString sAdditionalTrainConfig      = cuts.GetSpecialSettingFromAddConfig(additionalTrainConfig, "", "", addTaskName);


  if (sAdditionalTrainConfig.Atoi() > 0){
    trainConfig = trainConfig + additionalTrainConfig.Atoi();
    cout << "INFO: " << addTaskName.Data() << " running additionalTrainConfig '" << sAdditionalTrainConfig.Atoi() << "', train config: '" << trainConfig << "'" << endl;

  }else {
    cout << "INFO: " << " NO ADDITIONAL TRAIN CONFIG"<< endl;
  }

  // ================== GetAnalysisManager ===============================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
     Error(Form("%s_%i", addTaskName.Data(),  trainConfig), "No analysis manager found.");
     return ;
  }

  // ================== GetInputEventHandler =============================
  AliVEventHandler *inputHandler=mgr->GetInputEventHandler();

  //========= Add PID Reponse to ANALYSIS manager ====
  if(!(AliPIDResponse*)mgr->GetTask("PIDResponseTask")){
    Error(Form("AddTask_TrackQA_%i",trainConfig), "No PID response has been initialized aborting.");
    return;
  }

  Bool_t isMCForOtherSettings = 0;
  if (isMC > 0) isMCForOtherSettings = 1;


  // // enableConstructGamma is not fully developped . drop this option
  // // Bool_t  enableConstructGamma=kFALSE;
  // // if (trainConfig>200){
  // //   enableConstructGamma=kTRUE;
  // // }
  // if (trainConfig>100 && trainConfig<200 ){
  //   cutnumberPhoton = "10000000000000000500004000";
  // }


 //=========  Set Cutnumber for V0Reader ================================
  TString cutnumberPhoton;
  // if ( trainConfig<100 ){
  cutnumberPhoton = photonCutNumberV0Reader.Data();
  // }  else if (trainConfig>100 && trainConfig<200 ){
  //   cutnumberPhoton = photonCutNumberOfflineV0Reader.Data();
  // }
  TString cutnumberEvent      = "00000003";
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();

  //========= Check V0 Reader in  ANALYSIS manager  =====
  TString V0ReaderName        = Form("V0ReaderV1_%s_%s",cutnumberEvent.Data(),cutnumberPhoton.Data());
  AliV0ReaderV1 *fV0ReaderV1  =  NULL;
  fV0ReaderV1 = (AliV0ReaderV1*)mgr->GetTask(V0ReaderName.Data());

  if( !(fV0ReaderV1) ){
    cout << "V0Reader: " << V0ReaderName.Data() << " not found!!"<< endl;
    return;
  } else {
    cout << "V0Reader: " << V0ReaderName.Data() << " found!!"<< endl;
  }
  fV0ReaderV1->AliV0ReaderV1::SetProduceV0FindingEfficiency(enableV0findingEffi);

  cout<< "  "<< endl;
  cout<< " Done with the V0Reader "<< endl;
  cout<< "  "<< endl;

 //========= Add Pion Selector ====================
  TString PionCuts          = "000000200";            //Pion Cuts
  //  TString PionSelectorName = Form("PionSelector_%s", PionCuts.Data());
  TString PionSelectorName = Form("PionSelector");
  if( !(AliPrimaryPionSelector*)mgr->GetTask(PionSelectorName.Data()) ){
    AliPrimaryPionSelector *fPionSelector = new AliPrimaryPionSelector(PionSelectorName.Data());
    AliPrimaryPionCuts *fPionCuts=0;
    if( PionCuts!=""){
      fPionCuts= new AliPrimaryPionCuts(PionCuts.Data(),PionCuts.Data());
      //      if(runLightOutput>0) fPionCuts->SetLightOutput(kTRUE);
      if(fPionCuts->InitializeCutsFromCutString(PionCuts.Data())){
        fPionSelector->SetPrimaryPionCuts(fPionCuts);
        fPionCuts->SetFillCutHistograms("Pion",kTRUE);
      }
    }

    fPionSelector->Init();
    mgr->AddTask(fPionSelector);

    AliAnalysisDataContainer *cinput1  = mgr->GetCommonInputContainer();
    //connect input V0Reader
    mgr->ConnectInput (fPionSelector,0,cinput1);

  } 


  cout<< "  "<< endl;
  cout<< " Done with general pion selector "<< endl;
  cout<< "  "<< endl;

//========= Add Kaon Selector ====================
  TString KaonCuts          = "000000200";            //Kaon Cuts
  TString KaonSelectorName = Form("KaonSelector");
  if( !(AliPrimaryKaonSelector*)mgr->GetTask(KaonSelectorName.Data()) ){
    AliPrimaryKaonSelector *fKaonSelector = new AliPrimaryKaonSelector(KaonSelectorName.Data());
    AliPrimaryKaonCuts *fKaonCuts=0;
    if( KaonCuts!=""){
      fKaonCuts= new AliPrimaryKaonCuts(KaonCuts.Data(),KaonCuts.Data());
      //      if(runLightOutput>0) fPionCuts->SetLightOutput(kTRUE);
      if(fKaonCuts->InitializeCutsFromCutString(KaonCuts.Data())){
        fKaonSelector->SetPrimaryKaonCuts(fKaonCuts);
        fKaonCuts->SetFillCutHistograms("Kaon",kTRUE);
      }
    }

    fKaonSelector->Init();
    mgr->AddTask(fKaonSelector);

    AliAnalysisDataContainer *cinput1  = mgr->GetCommonInputContainer();
    //connect input V0Reader
    mgr->ConnectInput (fKaonSelector,0,cinput1);

  }

  cout<< "  "<< endl;
  cout<< "  "<< endl;
  cout<< "  "<< endl;

 //========= Add Proton Selector ====================
  TString ProtonCuts          = "000000200";            //Proton Cuts
  TString ProtonSelectorName = Form("ProtonSelector");
  if( !(AliPrimaryProtonSelector*)mgr->GetTask(ProtonSelectorName.Data()) ){
    AliPrimaryProtonSelector *fProtonSelector = new AliPrimaryProtonSelector(ProtonSelectorName.Data());
    AliPrimaryProtonCuts *fProtonCuts=0;
    if( ProtonCuts!=""){
      fProtonCuts= new AliPrimaryProtonCuts(ProtonCuts.Data(),ProtonCuts.Data());
      //     if(runLightOutput>0) fPionCuts->SetLightOutput(kTRUE);
      if(fProtonCuts->InitializeCutsFromCutString(ProtonCuts.Data())){
        fProtonSelector->SetPrimaryProtonCuts(fProtonCuts);
        fProtonCuts->SetFillCutHistograms("Proton",kTRUE);
      }
    }

    fProtonSelector->Init();
    mgr->AddTask(fProtonSelector);

    AliAnalysisDataContainer *cinput1  = mgr->GetCommonInputContainer();
    //connect input V0Reader
    mgr->ConnectInput (fProtonSelector,0,cinput1);
  }

  cout<< "  "<< endl;
  cout<< "  "<< endl;
  cout<< "  "<< endl;
 //========= Add Deuteron Selector ====================
  TString DeuteronCuts          = "000000200";            //Deuteron Cuts
  TString DeuteronSelectorName = Form("DeuteronSelector");
  if( !(AliPrimaryDeuteronSelector*)mgr->GetTask(DeuteronSelectorName.Data()) ){
    AliPrimaryDeuteronSelector *fDeuteronSelector = new AliPrimaryDeuteronSelector(DeuteronSelectorName.Data());
    AliPrimaryDeuteronCuts *fDeuteronCuts=0;
    if( DeuteronCuts!=""){
      fDeuteronCuts= new AliPrimaryDeuteronCuts(DeuteronCuts.Data(),DeuteronCuts.Data());
      //     if(runLightOutput>0) fDeuteronCuts->SetLightOutput(kTRUE);
      if(fDeuteronCuts->InitializeCutsFromCutString(DeuteronCuts.Data())){
        fDeuteronSelector->SetPrimaryDeuteronCuts(fDeuteronCuts);
        fDeuteronCuts->SetFillCutHistograms("Deuteron",kTRUE);
      }
    }

    fDeuteronSelector->Init();
    mgr->AddTask(fDeuteronSelector);

    AliAnalysisDataContainer *cinput1  = mgr->GetCommonInputContainer();
    //connect input V0Reader
    mgr->ConnectInput (fDeuteronSelector,0,cinput1);

  }

  //================================================
  //========= Add task to the ANALYSIS manager =====
  //            find input container

  cout<< "Starting with my own TrakQA Task"<< endl;

  AliAnalysisTaskTrackQA *fTrackQA=NULL;
  fTrackQA= new AliAnalysisTaskTrackQA(Form("TrackQA_%i",trainConfig));
  cout<< " " << endl;
  cout<< " TrackQA-1"<< endl;

  fTrackQA->SetIsMC(isMC);
  fTrackQA->SetIsHeavyIon(IsHeavyIon);
  fTrackQA->SetV0ReaderName(V0ReaderName);
  fTrackQA->SetPionSelectorName(PionSelectorName);
  fTrackQA->SetKaonSelectorName(KaonSelectorName);
  fTrackQA->SetProtonSelectorName(ProtonSelectorName);
  fTrackQA->SetDeuteronSelectorName(DeuteronSelectorName);

  cout<< "  " << endl;
  cout<< " TrackQA-2"<< endl;


  // fTrackQA->SetDoDeDxMaps(doDeDxMaps);
  // if (doMultiplicityWeighting>0 && isMC==0){
  //   doMultiplicityWeighting = 0;
  // }

  //  fTrackQA->SetDoMultWeights(doMultiplicityWeighting);
  // CutHandlerConvMaterial cuts;
  if(trainConfig == 1){
    cuts.AddCutTrackQA("00000103", "429400700", "429400700", "429400700", "429400700");
    //cuts.AddCutTrackQAPion("00000103", "429400700");
  } else if (trainConfig == 2) {
    cuts.AddCutTrackQA("00000103", "329400700", "329400700", "329400700", "329400700");
    //cuts.AddCutTrackQAPion("00000103", "329400700");
  } else if (trainConfig == 3) {
    cuts.AddCutTrackQA("00000103", "329000500", "329000500", "329000500", "329000500");
    //cuts.AddCutTrackQAPion("00000103", "329400700");

  }else {
    cout<< " should not be here" << endl;
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
  // TList *ConvCutList = new TList();
  TList *PionCutList  = new TList();
  TList *KaonCutList  = new TList();
  TList *ProtonCutList  = new TList();
  TList *DeuteronCutList  = new TList();

  EventCutList->SetOwner(kTRUE);
  AliConvEventCuts **analysisEventCuts        = new AliConvEventCuts*[numberOfCuts];

  //ConvCutList->SetOwner(kTRUE);
  //  AliConversionPhotonCuts **analysisCuts      = new AliConversionPhotonCuts*[numberOfCuts];


  PionCutList->SetOwner(kTRUE);
  AliPrimaryPionCuts **analysisPionCuts     = new AliPrimaryPionCuts*[numberOfCuts];

  KaonCutList->SetOwner(kTRUE);
  AliPrimaryKaonCuts **analysisKaonCuts     = new AliPrimaryKaonCuts*[numberOfCuts];

  ProtonCutList->SetOwner(kTRUE);
  AliPrimaryProtonCuts **analysisProtonCuts     = new AliPrimaryProtonCuts*[numberOfCuts];

  DeuteronCutList->SetOwner(kTRUE);
  AliPrimaryDeuteronCuts **analysisDeuteronCuts     = new AliPrimaryDeuteronCuts*[numberOfCuts];




  cout<<"names"<< periodNameAnchor.Data()<< " " << periodName.Data()<< endl;
  for(Int_t i = 0; i<numberOfCuts; i++){
    TString cutName( Form("%s_%s_%s_%s_%s",(cuts.GetEventCut(i)).Data(), 
    			  (cuts.GetPionCut(i)).Data(), 
    			  (cuts.GetKaonCut(i)).Data(),
    			  (cuts.GetProtonCut(i)).Data(),
    			  (cuts.GetDeuteronCut(i)).Data()  ) );
 
  // TString cutName( Form("%s_%s",(cuts.GetEventCut(i)).Data(), 
  // 			  (cuts.GetPionCut(i)).Data() ) );
    
    analysisEventCuts[i]          = new AliConvEventCuts();
    analysisEventCuts[i]->SetV0ReaderName(V0ReaderName);
    analysisEventCuts[i]->InitializeCutsFromCutString((cuts.GetEventCut(i)).Data());
    EventCutList->Add(analysisEventCuts[i]);
    analysisEventCuts[i]->SetFillCutHistograms("",kTRUE);


    analysisPionCuts[i] = new AliPrimaryPionCuts();
    //analysisPionCuts[i] ->SetPionSelectorName(PionSelectorName);
    //    if(runLightOutput>0) analysisPionCuts[i]->SetLightOutput(kTRUE);
    if( !analysisPionCuts[i]->InitializeCutsFromCutString((cuts.GetPionCut(i)).Data())) {
      cout<< "ERROR:  analysisPionCuts [ " <<i<<" ] "<<endl;
      return 0;
    } else {
      PionCutList->Add(analysisPionCuts[i]);
      analysisPionCuts[i]->SetFillCutHistograms("Pion",kFALSE,cutName);
    }

    analysisKaonCuts[i] = new AliPrimaryKaonCuts();
    //    if(runLightOutput>0) analysisKaonCuts[i]->SetLightOutput(kTRUE);
    if( !analysisKaonCuts[i]->InitializeCutsFromCutString((cuts.GetKaonCut(i)).Data())) {
      cout<< "ERROR:  analysisKaonCuts [ " <<i<<" ] "<<endl;
      return 0;
    } else {
      KaonCutList->Add(analysisKaonCuts[i]);
      analysisKaonCuts[i]->SetFillCutHistograms("Kaon",kFALSE,cutName);
    }

    analysisProtonCuts[i] = new AliPrimaryProtonCuts();
    //    if(runLightOutput>0) analysisProtonCuts[i]->SetLightOutput(kTRUE);
    if( !analysisProtonCuts[i]->InitializeCutsFromCutString((cuts.GetProtonCut(i)).Data())) {
      cout<< "ERROR:  analysisProtonCuts [ " <<i<<" ] "<<endl;
      return 0;
    } else {
      ProtonCutList->Add(analysisProtonCuts[i]);
      analysisProtonCuts[i]->SetFillCutHistograms("Proton",kFALSE,cutName);
    }

    analysisDeuteronCuts[i] = new AliPrimaryDeuteronCuts();
    //    if(runLightOutput>0) analysisDeuteronCuts[i]->SetLightOutput(kTRUE);
    if( !analysisDeuteronCuts[i]->InitializeCutsFromCutString((cuts.GetDeuteronCut(i)).Data())) {
      cout<< "ERROR:  analysisDeuteronCuts [ " <<i<<" ] "<<endl;
      return 0;
    } else {
      DeuteronCutList->Add(analysisDeuteronCuts[i]);
      analysisDeuteronCuts[i]->SetFillCutHistograms("Deuteron",kFALSE,cutName);
    }
 
  }

  fTrackQA->SetEventCutList(numberOfCuts,EventCutList);
  fTrackQA->SetPionCutList(PionCutList);
  fTrackQA->SetKaonCutList(KaonCutList);
  fTrackQA->SetProtonCutList(ProtonCutList);
  fTrackQA->SetDeuteronCutList(DeuteronCutList);

  //  fTrackQA->SetConversionCutList(numberOfCuts,ConvCutList);
  mgr->AddTask(fTrackQA);


  AliAnalysisDataContainer *coutput =
    mgr->CreateContainer(Form("TrackQA_%i",trainConfig), TList::Class(),
			 AliAnalysisManager::kOutputContainer,Form("TrackQA_%i.root",trainConfig ));

  mgr->ConnectInput(fTrackQA,  0, cinput );
  mgr->ConnectOutput(fTrackQA,  1, coutput);
  //connect containers
	return;
}
