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
//($ALIPHYSICS/PWGGA/GammaConv/AliAnalysisTaskMaterialHistos.cxx) for
//pp together with all supporting classes
//***************************************************************************************



void AddTask_MaterialHistos_pp( Int_t   trainConfig             = 1,                  // change different set of cuts
                                Int_t   isMC                    = 0,    
				TString   photonCutNumberV0Reader = "",       // 00000008400000000100000000 nom. B, 00000088400000000100000000 low B
				TString   photonCutNumberOfflineV0Reader = "",       // 10000008400000000100000000 nom. B, 10000088400000000100000000 low B
				TString periodName              = "",
                                TString periodNameAnchor        = "",
				Bool_t 	enableV0findingEffi     = kFALSE,    // enables V0finding efficiency histograms
				// settings for weights
				// FPTW:fileNamePtWeights, FMUW:fileNameMultWeights, FMAW:fileNameMatBudWeights, FEPC:fileNamedEdxPostCalib, FGPW:fileNameGammaPtWeights separate with ;
				TString fileNameExternalInputs        = "",
				Int_t   doDeDxMaps              =  0,
 				Bool_t  enableElecDeDxPostCalibration = kFALSE,
				Int_t   doMultiplicityWeighting = 0,
				Int_t   doWeightingGamma        = 0,        // enable Weighting
				Int_t   enableMatBudWeightsPi0  = 0,        // 1 = three radial bins, 2 = 10 radial bins
				TString additionalTrainConfig   = "0"       // additional counter for trainconfig, this has to be always the last parameter
                              ){

				//	TString fileNameInputForMultWeighing    = "/alice/cern.ch/user/a/amarin/multWeightFile/histosForMultWeighting_pp.root",            
  				//      TString fileNameElecDeDxPostCalibration = "dEdxCorrectionMap_Period_Pass.root",

  Int_t IsHeavyIon = 0;
  AliCutHandlerPCM cuts;
  TString fileNamePtWeights           = cuts.GetSpecialFileNameFromString (fileNameExternalInputs, "FPTW:");
  TString fileNameMultWeights         = cuts.GetSpecialFileNameFromString (fileNameExternalInputs, "FMUW:");
  TString fileNameMatBudWeights       = cuts.GetSpecialFileNameFromString (fileNameExternalInputs, "FMAW:");
  TString fileNamedEdxPostCalib       = cuts.GetSpecialFileNameFromString (fileNameExternalInputs, "FEPC:");
  TString fileNameGammaPtWeights      = cuts.GetSpecialFileNameFromString (fileNameExternalInputs, "FGPW:");

  TString addTaskName                 = "AddTask_MaterialHistos_pp";
  TString sAdditionalTrainConfig      = cuts.GetSpecialSettingFromAddConfig(additionalTrainConfig, "", "", addTaskName);


  if (sAdditionalTrainConfig.Atoi() > 0){
    trainConfig = trainConfig + additionalTrainConfig.Atoi();
    cout << "INFO: " << addTaskName.Data() << " running additionalTrainConfig '" << sAdditionalTrainConfig.Atoi() << "', train config: '" << trainConfig << "'" << endl;

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
    Error(Form("AddTask_MaterialHistos_%i",trainConfig), "No PID response has been initialized aborting.");
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
  if ( trainConfig<100 ){
    cutnumberPhoton = photonCutNumberV0Reader.Data();
  }  else if (trainConfig>100 && trainConfig<200 ){
    cutnumberPhoton = photonCutNumberOfflineV0Reader.Data();
  }
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

  //================================================
  //========= Add task to the ANALYSIS manager =====
  //            find input container


  AliAnalysisTaskMaterialHistos *fMaterialHistos=NULL;
  fMaterialHistos= new AliAnalysisTaskMaterialHistos(Form("MaterialHistos_%i",trainConfig));
  fMaterialHistos->SetIsMC(isMC);
  fMaterialHistos->SetIsHeavyIon(IsHeavyIon);
  fMaterialHistos->SetV0ReaderName(V0ReaderName);
  fMaterialHistos->SetDoDeDxMaps(doDeDxMaps);
  if (doMultiplicityWeighting>0 && isMC==0){
    doMultiplicityWeighting = 0;
  }

  fMaterialHistos->SetDoMultWeights(doMultiplicityWeighting);
  // CutHandlerConvMaterial cuts;
  if(trainConfig == 1){
    cuts.AddCutPCMMaterial("00000103", "00000009266302004204400000");
    //     cuts.AddCutPCMMaterial("00000103", "00000009286372004204400000");
  } else if (trainConfig == 2) {
    cuts.AddCutPCMMaterial("00000103", "00000009286342004204400000");
    cuts.AddCutPCMMaterial("00000103", "00000009286342001204400000");
    cuts.AddCutPCMMaterial("00000103", "00000009286342007204400000");
 } else if (trainConfig == 3) {
    cuts.AddCutPCMMaterial("00000103", "00000009286342004254400000");
    cuts.AddCutPCMMaterial("00000103", "00000009286372004254400000");
    cuts.AddCutPCMMaterial("00000103", "00000009286372004204400000");
 } else if (trainConfig == 4) {
    cuts.AddCutPCMMaterial("00000103", "00000009217300005800000000");
    cuts.AddCutPCMMaterial("00000103", "00000009217320005800000000");
    cuts.AddCutPCMMaterial("00000103", "00000009395074005200000000");
    cuts.AddCutPCMMaterial("00000103", "00000009395074005284400000");
 } else if (trainConfig == 5) {
    cuts.AddCutPCMMaterial("00000103", "00000009266320005800000000");
    cuts.AddCutPCMMaterial("00000103", "00000009266320005854400000");
    cuts.AddCutPCMMaterial("00000103", "00000009366320005804400000");
    cuts.AddCutPCMMaterial("00000103", "00000009366300005804400000");
    cuts.AddCutPCMMaterial("00000103", "00000009366300005854400000");
 } else if (trainConfig == 6) {
    cuts.AddCutPCMMaterial("00010103", "00000009247602008250404000"); // INT7
 } else if (trainConfig == 7) {
    cuts.AddCutPCMMaterial("00000103", "00000009266320005800004000");
    cuts.AddCutPCMMaterial("00000103", "00000009266320005854404000");
    cuts.AddCutPCMMaterial("00000103", "00000009366320005804404000");
    cuts.AddCutPCMMaterial("00000103", "00000009366300005804404000");
    cuts.AddCutPCMMaterial("00000103", "00000009366300005854404000");
 } else if (trainConfig == 8  || trainConfig ==  1008 || trainConfig ==  2008 || trainConfig ==  3008 ) {
    cuts.AddCutPCMMaterial("00000103", "00000009266300008804004000");
    cuts.AddCutPCMMaterial("00000103", "00000009266300008800404000");
    cuts.AddCutPCMMaterial("00000103", "00000009266370008804004000");
 } else if (trainConfig == 9) {
    cuts.AddCutPCMMaterial("00010103", "00000009266300008804004000");
    cuts.AddCutPCMMaterial("00010103", "00000009266300008800404000");
 } else if (trainConfig == 10) {
    cuts.AddCutPCMMaterial("00010103", "00000009266300008804004000");
    cuts.AddCutPCMMaterial("00010103", "00000009266300008800404000");
    cuts.AddCutPCMMaterial("00010103", "00000009266370008804004000");
 } else if (trainConfig == 21) {
    //    cuts.AddCutPCMMaterial("00010103", "00000009266300008884404000"); //- AM 28.03.18
    cuts.AddCutPCMMaterial("00010103", "00000009266300008854404000");
    cuts.AddCutPCMMaterial("00010103", "00000009266300008850404000");
    cuts.AddCutPCMMaterial("00010103", "00000009266300008254404000");
    cuts.AddCutPCMMaterial("00010103", "00000009227300008854404000");
    //    cuts.AddCutPCMMaterial("00010103", "00000009266370008854404000");

 } else if (trainConfig == 22) {
    cuts.AddCutPCMMaterial("00000103", "00000009266300008884404000");
    cuts.AddCutPCMMaterial("00000103", "00000009266300008854404000");
    cuts.AddCutPCMMaterial("00000103", "00000009266370008854404000");
 } else if (trainConfig == 23) {
    //    cuts.AddCutPCMMaterial("00010103", "00000009266300008884404000"); //- AM 28.03.18
    cuts.AddCutPCMMaterial("00010103", "0c000009266300008850404000");
    cuts.AddCutPCMMaterial("00010103", "0d000009266300008850404000");
 } else if (trainConfig == 24) {
    cuts.AddCutPCMMaterial("00010103", "0d000009266300008858404000");
    cuts.AddCutPCMMaterial("00010103", "0d0000d9266300008850404000");   // increased pT to 60 MeV for e+e-
    cuts.AddCutPCMMaterial("00010103", "00000009266300008750404000");
    cuts.AddCutPCMMaterial("00010103", "00000009266300008650404000");
 } else if (trainConfig == 25) {
    cuts.AddCutPCMMaterial("00010103", "00000009a27300008250a04120");
    cuts.AddCutPCMMaterial("00010103", "0c000009a27300008250a04120");
    cuts.AddCutPCMMaterial("00010103", "0d000009a27300008250a04120");
  } else if (trainConfig == 26) {   // new Cuts on Psi pair, chi2
    cuts.AddCutPCMMaterial("00010103", "0d00000929730000dgd0404000");
    cuts.AddCutPCMMaterial("00010103", "0d00000926630000dgd0404000");
    cuts.AddCutPCMMaterial("00010103", "0d00000926630000dkd0404000");
  } else if (trainConfig == 27) {  // to calculate MBW for PHOS pi0 region
    cuts.AddCutPCMMaterial("00010103", "0700bb09266300008884404000"); 
    cuts.AddCutPCMMaterial("00010103", "0800bb09266300008850404000");

  } else if (trainConfig == 28) {   // new Cuts on Psi pair, chi2, to calculate MBW for PHOS pi0 region
    cuts.AddCutPCMMaterial("00010103", "0800bb0929730000dgd0404000");
    cuts.AddCutPCMMaterial("00010103", "0800bb0926630000dgd0404000");
    cuts.AddCutPCMMaterial("00010103", "0800bb0926630000dkd0404000");



 } else if (trainConfig == 76) {   // +50 To be used with MBW
    cuts.AddCutPCMMaterial("00010103", "0d00000929730000dgd0404000");
    cuts.AddCutPCMMaterial("00010103", "0d00000926630000dgd0404000");
    cuts.AddCutPCMMaterial("00010103", "0d00000926630000dkd0404000");




   // Offline V0Finder is used

  } else  if(trainConfig == 101){
    cuts.AddCutPCMMaterial("00000103", "10000009266302004204400000");
    //     cuts.AddCutPCMMaterial("00000103", "00000009286372004204400000");
  } else if (trainConfig == 102) {
    cuts.AddCutPCMMaterial("00000103", "10000009286342004204400000");
    cuts.AddCutPCMMaterial("00000103", "10000009286342001204400000");
    cuts.AddCutPCMMaterial("00000103", "10000009286342007204400000");
  } else if (trainConfig == 103) {
    cuts.AddCutPCMMaterial("00000103", "10000009286342004254400000");
    cuts.AddCutPCMMaterial("00000103", "10000009286372004254400000");
    cuts.AddCutPCMMaterial("00000103", "10000009286372004204400000");
  } else if (trainConfig == 104) {
    cuts.AddCutPCMMaterial("00000103", "10000009217300005800000000");
    cuts.AddCutPCMMaterial("00000103", "10000009217320005800000000");
    cuts.AddCutPCMMaterial("00000103", "10000009395074005200000000");
    cuts.AddCutPCMMaterial("00000103", "10000009395074005284400000");
 } else if (trainConfig == 105) {
    cuts.AddCutPCMMaterial("00000103", "10000009266320005800000000");
    cuts.AddCutPCMMaterial("00000103", "10000009266320005854400000");
    cuts.AddCutPCMMaterial("00000103", "10000009366320005804400000");
    cuts.AddCutPCMMaterial("00000103", "10000009366300005804400000");
    cuts.AddCutPCMMaterial("00000103", "10000009366300005854400000");
 } else if (trainConfig == 107) {
    cuts.AddCutPCMMaterial("00000103", "10000009266320005800004000");
    cuts.AddCutPCMMaterial("00000103", "10000009266320005854404000");
    cuts.AddCutPCMMaterial("00000103", "10000009366320005804404000");
    cuts.AddCutPCMMaterial("00000103", "10000009366300005804404000");
    cuts.AddCutPCMMaterial("00000103", "10000009366300005854404000");
 } else if (trainConfig == 108  || trainConfig ==  1108 || trainConfig ==  2108 || trainConfig ==  3108 ) {
    cuts.AddCutPCMMaterial("00000103", "10000009266300008804004000");
    cuts.AddCutPCMMaterial("00000103", "10000009266300008800404000");
    cuts.AddCutPCMMaterial("00000103", "10000009266370008804004000");
 } else if (trainConfig == 109) {
    cuts.AddCutPCMMaterial("00010103", "10000009266300008804004000");
    cuts.AddCutPCMMaterial("00010103", "10000009266300008800404000");
 } else if (trainConfig == 110) {
    cuts.AddCutPCMMaterial("00010103", "10000009266300008804004000");
    cuts.AddCutPCMMaterial("00010103", "10000009266300008800404000");
    cuts.AddCutPCMMaterial("00010103", "10000009266370008804004000");
 } else if (trainConfig == 121) {
    //    cuts.AddCutPCMMaterial("00010103", "10000009266300008884404000"); //-AM 28.03.18
    //    cuts.AddCutPCMMaterial("00010103", "10000009266370008854404000");

    cuts.AddCutPCMMaterial("00010103", "10000009266300008854404000");
    cuts.AddCutPCMMaterial("00010103", "10000009266300008850404000");
    cuts.AddCutPCMMaterial("00010103", "10000009266300008254404000");
    cuts.AddCutPCMMaterial("00010103", "10000009227300008854404000");

 } else if (trainConfig == 122) {
    cuts.AddCutPCMMaterial("00000103", "10000009266300008884404000");
    cuts.AddCutPCMMaterial("00000103", "10000009266300008854404000");
    cuts.AddCutPCMMaterial("00000103", "10000009266370008854404000");
 } else if (trainConfig == 123) {
    cuts.AddCutPCMMaterial("00010103", "1c000009266300008850404000");
    cuts.AddCutPCMMaterial("00010103", "1d000009266300008850404000");
 } else if (trainConfig == 124) {
    cuts.AddCutPCMMaterial("00010103", "1d000009266300008858404000");
    cuts.AddCutPCMMaterial("00010103", "1d0000d9266300008850404000");   // increased pT to 60 MeV for e+e-
    cuts.AddCutPCMMaterial("00010103", "10000009266300008750404000");
    cuts.AddCutPCMMaterial("00010103", "10000009266300008650404000");
 } else if (trainConfig == 125) {
    cuts.AddCutPCMMaterial("00010103", "10000009a27300008250a04120");
    cuts.AddCutPCMMaterial("00010103", "1c000009a27300008250a04120");
    cuts.AddCutPCMMaterial("00010103", "1d000009a27300008250a04120");

 } else if (trainConfig == 126) {   // new Cuts on Psi pair, chi2
    cuts.AddCutPCMMaterial("00010103", "1d00000929730000dgd0404000");
    cuts.AddCutPCMMaterial("00010103", "1d00000926630000dgd0404000");
    cuts.AddCutPCMMaterial("00010103", "1d00000926630000dkd0404000");
 } else if (trainConfig == 127) {  // to calculate MBW for PHOS pi0 region
    cuts.AddCutPCMMaterial("00010103", "1700bb09266300008884404000"); 
    cuts.AddCutPCMMaterial("00010103", "1800bb09266300008850404000");

  } else if (trainConfig == 128) {   // new Cuts on Psi pair, chi2, to calculate MBW for PHOS pi0 region
    cuts.AddCutPCMMaterial("00010103", "1800bb0929730000dgd0404000");
    cuts.AddCutPCMMaterial("00010103", "1800bb0926630000dgd0404000");
    cuts.AddCutPCMMaterial("00010103", "1800bb0926630000dkd0404000");



 } else if (trainConfig == 176) {   // +50 To be used with MBW
    cuts.AddCutPCMMaterial("00010103", "1d00000929730000dgd0404000");
    cuts.AddCutPCMMaterial("00010103", "1d00000926630000dgd0404000");
    cuts.AddCutPCMMaterial("00010103", "1d00000926630000dkd0404000");
    
  } else  if(trainConfig == 111){
    cuts.AddCutPCMMaterial("00000003", "10000070000000000500004000");
  // 7 TeV testconfig for pileup checks
  } else if (trainConfig == 150) {
    cuts.AddCutPCMMaterial("00000103", "00000009266300008884004000");
    cuts.AddCutPCMMaterial("00000103", "00000009266300008804004000");
    cuts.AddCutPCMMaterial("00000103", "10000009266300008884004000");
    cuts.AddCutPCMMaterial("00000103", "10000009266300008804004000");
    // 8 TeV testconfig for pileup checks
  } else if (trainConfig == 151) {
    cuts.AddCutPCMMaterial("00010103", "00000009266300008884004000");
    cuts.AddCutPCMMaterial("00010103", "00000009266300008804004000");
    cuts.AddCutPCMMaterial("00010103", "10000009266300008884004000");
    cuts.AddCutPCMMaterial("00010103", "10000009266300008804004000");

    // ConstructGamma is used

  // } else  if(trainConfig == 201){
  //   cuts.AddCutPCMMaterial("00000103", "00000009266302004204400000");
  //   //     cuts.AddCutPCMMaterial("00000103", "00000009286372004204400000");
  // } else if (trainConfig == 202) {
  //   cuts.AddCutPCMMaterial("00000103", "00000009286342004204400000");
  //   cuts.AddCutPCMMaterial("00000103", "00000009286342001204400000");
  //   cuts.AddCutPCMMaterial("00000103", "00000009286342007204400000");
  // } else if (trainConfig == 203) {
  //   cuts.AddCutPCMMaterial("00000103", "00000009286342004254400000");
  //   cuts.AddCutPCMMaterial("00000103", "00000009286372004254400000");
  //   cuts.AddCutPCMMaterial("00000103", "00000009286372004204400000");



  }else {
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
  cout<<"names"<< periodNameAnchor.Data()<< " " << periodName.Data()<< endl;
  Bool_t initializedMatBudWeigths_existing    = kFALSE;


  for(Int_t i = 0; i<numberOfCuts; i++){
    analysisEventCuts[i]          = new AliConvEventCuts();

    TString dataInputMultHisto  = "";
    TString mcInputMultHisto    = "";
    if (doMultiplicityWeighting>0){
      cout << "INFO enableling mult weighting" << endl;
      if( periodNameAnchor.CompareTo("LHC15n")==0  ||
	  periodNameAnchor.CompareTo("LHC16d")==0  ||
	  periodNameAnchor.CompareTo("LHC17p")==0  ||  
	  periodNameAnchor.CompareTo("LHC17q")==0  ){
	TString cutNumber = cuts.GetEventCut(i);
	TString centCut = cutNumber(0,3);  // first three digits of event cut
	if(doMultiplicityWeighting == 1){
	  dataInputMultHisto = Form("%s_%s", periodNameAnchor.Data(), centCut.Data());
	  mcInputMultHisto   = Form("%s_%s", periodName.Data(), centCut.Data());
	  cout<< "Histogram names data/MC:: "<< dataInputMultHisto.Data()<< " " << mcInputMultHisto.Data()<< endl;
	}else if(doMultiplicityWeighting == 2){
	  dataInputMultHisto = Form("V0M_%s_%s", periodNameAnchor.Data(), centCut.Data());
	  mcInputMultHisto   = Form("V0M_%s_%s", periodName.Data(), centCut.Data());
	  cout<< "Histogram names data/MC:: "<< dataInputMultHisto.Data()<< " " << mcInputMultHisto.Data()<< endl;
	}else{
	  dataInputMultHisto = Form("%s_%s", periodNameAnchor.Data(), centCut.Data());
	  mcInputMultHisto   = Form("%s_%s", periodName.Data(), centCut.Data());
	  cout<< "Histogram names data/MC:: "<< dataInputMultHisto.Data()<< " " << mcInputMultHisto.Data()<< endl;
	}
       }
      analysisEventCuts[i]->SetUseWeightMultiplicityFromFile(kTRUE, fileNameMultWeights, dataInputMultHisto, mcInputMultHisto );
    }

    TString dataInputPtHisto  = "";
    TString mcInputPtHisto    = "";

    if (doWeightingGamma) {
      cout << "INFO enableling gamma pT weighting" << endl;
      if( periodNameAnchor.CompareTo("LHC15n")==0  ||
	  periodNameAnchor.CompareTo("LHC16d")==0  ||
	  periodNameAnchor.CompareTo("LHC17p")==0  ||  
	  periodNameAnchor.CompareTo("LHC17q")==0  ){
	TString cutNumber = cuts.GetEventCut(i);
	TString centCut = cutNumber(0,3);  // first three digits of event cut
	dataInputPtHisto = Form("pT_%s_%s", periodNameAnchor.Data(), centCut.Data());
	mcInputPtHisto   = Form("pT_%s_%s", periodName.Data(), centCut.Data());
	cout<< "Histogram names data/MC:: "<< dataInputPtHisto.Data()<< " " << mcInputPtHisto.Data()<< endl;
       }
      analysisEventCuts[i]->SetUseGammaPtReweightingWithHistogramFromFile(kTRUE, fileNameGammaPtWeights, mcInputPtHisto, dataInputPtHisto);
    }
    analysisEventCuts[i]->SetV0ReaderName(V0ReaderName);
    analysisEventCuts[i]->InitializeCutsFromCutString((cuts.GetEventCut(i)).Data());
    EventCutList->Add(analysisEventCuts[i]);
    analysisEventCuts[i]->SetFillCutHistograms("",kTRUE);

    analysisCuts[i]               = new AliConversionPhotonCuts();

    if (enableMatBudWeightsPi0 > 0){
      cout<< "Material budget weigthing enabled"<< endl;
        if (isMC > 0){
            if (analysisCuts[i]->InitializeMaterialBudgetWeights(enableMatBudWeightsPi0,fileNameMatBudWeights)){
                initializedMatBudWeigths_existing = kTRUE;
		cout<< "Material budget weigthing enabled, went well"<< endl;
	    }
            else {cout << "ERROR The initialization of the materialBudgetWeights did not work out." << endl;}
        }
        else {cout << "ERROR 'enableMatBudWeightsPi0'-flag was set > 0 even though this is not a MC task. It was automatically reset to 0." << endl;}
    }

    if (enableElecDeDxPostCalibration){
      if (isMC == 0){
        if(fileNamedEdxPostCalib.CompareTo("") != 0){
          analysisCuts[i]->SetElecDeDxPostCalibrationCustomFile(fileNamedEdxPostCalib);
          cout << "Setting custom dEdx recalibration file: " << fileNamedEdxPostCalib.Data() << endl;
        }
        analysisCuts[i]->SetDoElecDeDxPostCalibration(enableElecDeDxPostCalibration);
        cout << "Enabled TPC dEdx recalibration." << endl;
      } else{
        cout << "ERROR enableElecDeDxPostCalibration set to True even if MC file. Automatically reset to 0"<< endl;
        enableElecDeDxPostCalibration=kFALSE;
        analysisCuts[i]->SetDoElecDeDxPostCalibration(kFALSE);
      }
    }

    analysisCuts[i]->SetV0ReaderName(V0ReaderName);
    analysisCuts[i]->InitializeCutsFromCutString((cuts.GetPhotonCut(i)).Data());
    ConvCutList->Add(analysisCuts[i]);
    analysisCuts[i]->SetFillCutHistograms("",kTRUE);
  }

  fMaterialHistos->SetEventCutList(numberOfCuts,EventCutList);
  fMaterialHistos->SetConversionCutList(numberOfCuts,ConvCutList);
  if (initializedMatBudWeigths_existing) {
      fMaterialHistos->SetDoMaterialBudgetWeightingOfGammasForTrueMesons(kTRUE);
  }

  mgr->AddTask(fMaterialHistos);


  AliAnalysisDataContainer *coutput =
    mgr->CreateContainer(Form("GammaConvMaterial_%i",trainConfig), TList::Class(),
			 AliAnalysisManager::kOutputContainer,Form("GammaConv_Material_%i.root",trainConfig ));

  mgr->ConnectInput(fMaterialHistos,  0, cinput );
  mgr->ConnectOutput(fMaterialHistos,  1, coutput);
  //connect containers
	return;
}
