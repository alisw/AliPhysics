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

void AddTask_MaterialHistos_pPb(
    Int_t trainConfig = 1, // change different set of cuts
    Int_t isMC = 0,
    TString periodName = "", // period name
    TString periodNameV0Reader = "",
    TString periodNameAnchor = "",
    Int_t doDeDxMaps = 0,
    Bool_t enableV0findingEffi = kFALSE, // enables V0finding efficiency histograms
    Bool_t enableElecDeDxPostCalibration = kFALSE,
    TString fileNameElecDeDxPostCalibration = "dEdxCorrectionMap_Period_Pass.root",
    Int_t doMultiplicityWeighting = 0,
    TString fileNameInputForMultWeighing = "/alice/cern.ch/user/a/amarin/multWeightFile/histosForMultWeighting_pp.root", //
    TString additionalTrainConfig = "0"                                                                                  // additional counter for trainconfig, this has to be always the last parameter
)
{

  Int_t IsHeavyIon = 0;
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

  //=========  Set Cutnumber for V0Reader ================================
  TString cutnumberPhoton = "00000009266300008804004000";
  //"00200008400000002200000000";

  TString cutnumberEvent = "00000103";
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();


  Bool_t  enableConstructGamma=kFALSE;
  if (trainConfig>200){
    enableConstructGamma=kTRUE;
  }

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
  fMaterialHistos->SetDoDeDxMaps(doDeDxMaps);
  if (doMultiplicityWeighting>0 && isMC==0){
    doMultiplicityWeighting = 0;
  }

  fMaterialHistos->SetDoMultWeights(doMultiplicityWeighting);
  CutHandlerConvMaterial cuts;
  if(trainConfig == 1){
    cuts.AddCut("80010103", "00000009266300008804004000");
  } else if (trainConfig == 2) {
    cuts.AddCut("80210103", "00000009266300008804004000");
    cuts.AddCut("82410103", "00000009266300008804004000");
    cuts.AddCut("84610103", "00000009266300008804004000");
    cuts.AddCut("86010103", "00000009266300008804004000");
  } else if (trainConfig == 3) {
    cuts.AddCut("80110103", "00000009266300008804004000");
    cuts.AddCut("81210103", "00000009266300008804004000");
    cuts.AddCut("82310103", "00000009266300008804004000");
    cuts.AddCut("83410103", "00000009266300008804004000");
  } else if (trainConfig == 4) {
    cuts.AddCut("a0110103", "00000009266300008804004000");
    cuts.AddCut("a1210103", "00000009266300008804004000");
    cuts.AddCut("c0110103", "00000009266300008804004000");
    cuts.AddCut("c0210103", "00000009266300008804004000");

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
  cout<<"names"<< periodNameAnchor.Data()<< " " << periodName.Data()<< " " << periodNameV0Reader.Data()<< endl;
  for (Int_t i = 0; i < numberOfCuts; i++)
  {
    analysisEventCuts[i] = new AliConvEventCuts();

    TString dataInputMultHisto = "";
    TString mcInputMultHisto = "";
    if (doMultiplicityWeighting > 0)
    {
      cout << "INFO enableling mult weighting" << endl;
      if (periodNameAnchor.CompareTo("LHC16d") == 0 ||
          periodNameAnchor.CompareTo("LHC17p") == 0 ||
          periodNameAnchor.CompareTo("LHC17q") == 0)
      {
        TString cutNumber = cuts.GetEventCut(i);
        TString centCut = cutNumber(0, 3); // first three digits of event cut
        dataInputMultHisto = Form("%s_%s", periodNameAnchor.Data(), centCut.Data());
        mcInputMultHisto = Form("%s_%s", periodName.Data(), centCut.Data());
        cout << "Histogram names data/MC:: " << dataInputMultHisto.Data() << " " << mcInputMultHisto.Data() << endl;
      }
      analysisEventCuts[i]->SetUseWeightMultiplicityFromFile(kTRUE, fileNameInputForMultWeighing, dataInputMultHisto, mcInputMultHisto);
    }

    analysisEventCuts[i]->SetV0ReaderName(V0ReaderName);
    analysisEventCuts[i]->InitializeCutsFromCutString((cuts.GetEventCut(i)).Data());
    EventCutList->Add(analysisEventCuts[i]);
    analysisEventCuts[i]->SetFillCutHistograms("", kTRUE);

    analysisCuts[i] = new AliConversionPhotonCuts();
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
    analysisCuts[i]->SetFillCutHistograms("", kTRUE);
  }

  fMaterialHistos->SetEventCutList(numberOfCuts, EventCutList);
  fMaterialHistos->SetConversionCutList(numberOfCuts, ConvCutList);
  mgr->AddTask(fMaterialHistos);

  AliAnalysisDataContainer *coutput =
      mgr->CreateContainer(Form("GammaConvMaterial_%i", trainConfig), TList::Class(),
                           AliAnalysisManager::kOutputContainer, Form("GammaConv_Material_%i.root", trainConfig));

  mgr->ConnectInput(fMaterialHistos, 0, cinput);
  mgr->ConnectOutput(fMaterialHistos, 1, coutput);
  //connect containers
  return;
}
