/////////////////////////////////////////////////////////////////////////////////////////////
//
// AddTaskGammaConvFlow_PbPb2 macro
// Author: Andrea Dubla, Redmer A. Bertens, Friederike Bock
//         Commented where necessary
/////////////////////////////////////////////////////////////////////////////////////////////

class AliAnalysisDataContainer;
class AliFlowTrackCuts;
class AliFlowTrackSimpleCuts;
class AliFlowEventCuts;
class AliFlowEventSimpleCuts;

class CutHandlerConvFlow{
  public:
    CutHandlerConvFlow(Int_t nMax=10){
      nCuts=0; nMaxCuts=nMax; validCuts = true;
      eventCutArray = new TString[nMaxCuts]; photonCutArray = new TString[nMaxCuts];
      for(Int_t i=0; i<nMaxCuts; i++) {eventCutArray[i] = ""; photonCutArray[i] = "";}
    }

    void AddCut(TString eventCut, TString photonCut){
      if(nCuts>=nMaxCuts) {cout << "ERROR in CutHandlerConvFlow: Exceeded maximum number of cuts!" << endl; validCuts = false; return;}
      if( eventCut.Length()!=8 || photonCut.Length()!=26) {cout << "ERROR in CutHandlerConvFlow: Incorrect length of cut string!" << endl; validCuts = false; return;}
      eventCutArray[nCuts]=eventCut; photonCutArray[nCuts]=photonCut;
      nCuts++;
      return;
    }
    Bool_t AreValid(){return validCuts;}
    Int_t GetNCuts(){if(validCuts) return nCuts; else return 0;}
    TString GetEventCut(Int_t i){if(validCuts&&i<nMaxCuts&&i>=0) return eventCutArray[i]; else{cout << "ERROR in CutHandlerConvFlow: GetEventCut wrong index i" << endl;return "";}}
    TString GetPhotonCut(Int_t i){if(validCuts&&i<nMaxCuts&&i>=0) return photonCutArray[i]; else {cout << "ERROR in CutHandlerConvFlow: GetPhotonCut wrong index i" << endl;return "";}}
  private:
    Bool_t validCuts;
    Int_t nCuts; Int_t nMaxCuts;
    TString* eventCutArray;
    TString* photonCutArray;
};

void AddTask_GammaConvFlow_PbPb2(
                                  TString uniqueID              = "",
                                  Int_t harmonic                = 2,
                                  Int_t trainConfig             = 1,                            //change different set of cuts
                                  Int_t enableQAMesonTask       = 0,                            //enable QA in AddTask_GammaConvFlow_PbPb2
                                  Int_t enableQAPhotonTask      = 0,                            // enable additional QA task
                                  //TString fileNameInputForWeighting = "MCSpectraInput.root",  // path to file for weigting input
                                  Bool_t doWeighting            = kFALSE,                       //enable Weighting
                                  TString cutnumberAODBranch    = "100000006008400000001500000",
                                  Bool_t BasicHistoSP           = kFALSE,
                                  Bool_t debug                  = kFALSE,
                                  Bool_t UseMassSel             = kFALSE,
                                  Float_t MinMass               = 0,
                                  Float_t MaxMass               = 0.2,
                                  Bool_t UseKappaSel            = kFALSE,
                                  Float_t MinKappa              = 0,
                                  Float_t MaxKappa              = 10,
                                  Int_t FilterVariable          = 1,                              // 1 = Mass, 2 = Kappa, 3 = MCee mass, 4 = MCee kappa, 5 = !MCee mass, 6 = !MCee kappa
                                  const Int_t NFilterBins       = 1,
                                  Double_t MinFilter            = 0.0,
                                  Double_t MaxFilter            = 0.2,
                                  Bool_t isMC                   = kFALSE,
                                  Int_t ApplydPhidRCut         = 0,
                                  Bool_t PerformExtraStudies    = kFALSE,                         // with kTRUE it performs the LTM study and dRdPhi study
                                  TString additionalTrainConfig = "0"                             // additional counter for trainconfig, always has to be last parameter
                               ) {

  Int_t isHeavyIon = 1;

  // make use of train subwagon feature
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

  //=========  Set Cutnumber for V0Reader ================================
  TString cutnumberPhoton = "30000008400100001500000000";
  TString cutnumberEvent = "10000003";
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
  //========= Add V0 Reader to  ANALYSIS manager if not yet existent =====
  TString V0ReaderName = Form("V0ReaderV1_%s_%s",cutnumberEvent.Data(),cutnumberPhoton.Data());
  if( !(AliV0ReaderV1*)mgr->GetTask(V0ReaderName.Data()) ){
      AliV0ReaderV1 *fV0ReaderV1 = new AliV0ReaderV1(V0ReaderName.Data());

      fV0ReaderV1->SetUseOwnXYZCalculation(kTRUE);
      fV0ReaderV1->SetCreateAODs(kFALSE);// AOD Output
      fV0ReaderV1->SetUseAODConversionPhoton(kTRUE);
      fV0ReaderV1->SetUseMassToZero(kFALSE);
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
          fCuts->SetIsHeavyIon(isHeavyIon);
          fCuts->SetV0ReaderName(V0ReaderName);
          if(fCuts->InitializeCutsFromCutString(cutnumberPhoton.Data())){
              fV0ReaderV1->SetConversionCuts(fCuts);
              fCuts->SetFillCutHistograms("",kTRUE);
          }
      }

      if(inputHandler->IsA()==AliAODInputHandler::Class()){
          // AOD mode
        fV0ReaderV1->AliV0ReaderV1::SetDeltaAODBranchName(Form("GammaConv_%s_gamma",cutnumberAODBranch.Data()));
      }
      fV0ReaderV1->Init();

      AliLog::SetGlobalLogLevel(AliLog::kInfo);

      //connect input V0Reader
      mgr->AddTask(fV0ReaderV1);
      mgr->ConnectInput(fV0ReaderV1,0,cinput);
  }


  CutHandlerConvFlow cuts;

  if (trainConfig == 1){
    cuts.AddCut("60100013", "04200009297002003220000000");
  } else if (trainConfig == 2) {
    cuts.AddCut("61200013", "04200009297002003220000000");
  } else if (trainConfig == 3) {
    cuts.AddCut("50100013", "04200009297002003220000000");
  } else if (trainConfig == 4) {
    cuts.AddCut("50200013", "04200009297002003220000000");
  } else if (trainConfig == 5) {
    cuts.AddCut("51200013", "04200009297002003220000000");
  } else if (trainConfig == 6) {
    cuts.AddCut("52400013", "04200009297002003220000000");
  } else if (trainConfig == 7) {
    cuts.AddCut("54600013", "04200009297002003220000000");
  } else if (trainConfig == 8) {
    cuts.AddCut("54800013", "04200009297002003220000000");
  } else if (trainConfig == 9) {
    cuts.AddCut("54500013", "04200009297002003220000000");
  } else if (trainConfig == 10) {
    cuts.AddCut("55600013", "04200009297002003220000000");
  } else if (trainConfig == 11) {
    cuts.AddCut("56800013", "04200009297002003220000000");
  } else if (trainConfig == 12) {
    cuts.AddCut("56700013", "04200009297002003220000000");
  } else if (trainConfig == 13) {
    cuts.AddCut("57800013", "04200009297002003220000000");
  } else if (trainConfig == 14) {
    cuts.AddCut("46900013", "04200009297002003220000000");
  } else if (trainConfig == 15) {
    cuts.AddCut("58900013", "04200009297002003220000000");
  } else  if (trainConfig == 16){
    cuts.AddCut("60100013", "00200009297002008250400000");
  } else if (trainConfig == 17) {
    cuts.AddCut("61200013", "00200009297002008250400000");
  } else if (trainConfig == 18) {
    cuts.AddCut("50100013", "00200009297002008250400000");
  } else if (trainConfig == 19) {
    cuts.AddCut("50200013", "00200009297002008250400000");
  } else if (trainConfig == 20) {
    cuts.AddCut("51200013", "00200009297002008250400000");
  } else if (trainConfig == 21) {
    cuts.AddCut("52400013", "00200009297002008250400000");
  } else if (trainConfig == 22) {
    cuts.AddCut("54600013", "00200009297002008250400000");
  } else if (trainConfig == 23) {
    cuts.AddCut("54800013", "00200009297002008250400000");
  } else if (trainConfig == 24) {
    cuts.AddCut("54500013", "00200009297002008250400000");
  } else if (trainConfig == 25) {
    cuts.AddCut("55600013", "00200009297002008250400000");
  } else if (trainConfig == 26) {
    cuts.AddCut("56800013", "00200009297002008250400000");
  } else if (trainConfig == 27) {
    cuts.AddCut("56700013", "00200009297002008250400000");
  } else if (trainConfig == 28) {
    cuts.AddCut("57800013", "00200009297002008250400000");
  } else if (trainConfig == 29) {
    cuts.AddCut("46900013", "00200009297002008250400000");
  } else if (trainConfig == 30) {
    cuts.AddCut("58900013", "00200009297002008250400000");
  } else if (trainConfig == 31) {
    cuts.AddCut("50800013", "00200009297002008250400000");
  } else if (trainConfig == 32) {
    cuts.AddCut("52500013", "00200009297002008250400000");
  } else if (trainConfig == 33) {
    cuts.AddCut("53500013", "00200009297002008250400000");
  } else if (trainConfig == 34) {
    cuts.AddCut("54500013", "00200009297002008250400000");
  } else if (trainConfig == 35) {
    cuts.AddCut("53400013", "00200009297002008250400000");
  } else if (trainConfig == 36) {
    cuts.AddCut("52300013", "00200009297002008250400000");
  } else if (trainConfig == 37){
    cuts.AddCut("60100013", "00200009297002208250400000");
  } else if (trainConfig == 38) {
    cuts.AddCut("61200013", "00200009297002208250400000");
  } else if (trainConfig == 39) {
    cuts.AddCut("51200013", "00200009297002208250400000");
  } else if (trainConfig == 40) {
    cuts.AddCut("54600013", "00200009297002208250400000");
  } else if (trainConfig == 41) {
    cuts.AddCut("56800013", "00200009297002208250400000");
  } else if (trainConfig == 42) {
    cuts.AddCut("53400013", "00200009297002208250400000");
  } else if (trainConfig == 43) {
    cuts.AddCut("52300013", "00200009297002208250400000");
  } else if (trainConfig == 44) {
    cuts.AddCut("50200013", "00200009297002008250400000");
  } else if (trainConfig == 45) {
    cuts.AddCut("50400013", "00200009297002008250400000");
  } else if (trainConfig == 46) {
    cuts.AddCut("61200013", "00200009697004000500000000");
  } else if (trainConfig == 47) {
    cuts.AddCut("61200013", "00200009697005000500000000");
  } else if (trainConfig == 48) {
    cuts.AddCut("61200013", "00200009797004000500000000");
  } else if (trainConfig == 49) {
    cuts.AddCut("61200013", "00200009797005000500000000");
  } else if (trainConfig == 50) {
    cuts.AddCut("52400013", "00200009007000008250400000");
  } else if (trainConfig == 51) {
    cuts.AddCut("50200013", "00200009007000008250400000");
  } else if (trainConfig == 52) {
    cuts.AddCut("52400013", "00200009007000008250400000");
  } else if (trainConfig == 53) {
    cuts.AddCut("54800013", "00200009007000008250400000");
  //Full standard cut output closed TPC cuts
  } else if (trainConfig == 54) {
    cuts.AddCut("50200013", "00200009297002008250400000");
    cuts.AddCut("52400013", "00200009297002008250400000");
    cuts.AddCut("54800013", "00200009297002008250400000");
  //Full standard cut output open TPC cuts
  } else if (trainConfig == 55) {
    cuts.AddCut("50200013", "00200009007000008250400000");
    cuts.AddCut("52400013", "00200009007000008250400000");
    cuts.AddCut("54800013", "00200009007000008250400000");
  //Full standard cut output open TPC cuts same sign pairing
  } else if (trainConfig == 56) {
    cuts.AddCut("50200013", "20200009007000008250400000");
    cuts.AddCut("52400013", "20200009007000008250400000");
    cuts.AddCut("54800013", "20200009007000008250400000");
  //Full standard cut output open TPC cuts same and unline sign pairing
  } else if (trainConfig == 57) {
    cuts.AddCut("50200013", "30200009007000008250400000");
    cuts.AddCut("52400013", "30200009007000008250400000");
    cuts.AddCut("54800013", "30200009007000008250400000");
  //Full standard cut output open TPC cuts - var #1 for std cuts - TOF and chi2
  } else if (trainConfig == 58) {
    cuts.AddCut("50200013", "00200009007002008750400000");
    cuts.AddCut("52400013", "00200009007002008750400000");
    cuts.AddCut("54800013", "00200009007002008750400000");
  //Full standard cut output open TPC cuts - var #2 for std cuts - tightening all
  } else if (trainConfig == 59) {
    cuts.AddCut("50200013", "00200009007002009760600000");
    cuts.AddCut("52400013", "00200009007002009760600000");
    cuts.AddCut("54800013", "00200009007002009760600000");
  //systematics 0-20%
  } else if (trainConfig == 60) {
    cuts.AddCut("50200013", "04200009007000008250400000"); //eta cut: |eta| <0.75
    cuts.AddCut("50200013", "00500009007000008250400000"); //minR= 10 maxR = 180
    cuts.AddCut("50200013", "00200049007000008250400000"); //singleptcut = 75MeV
    cuts.AddCut("50200013", "00200019007000008250400000"); //singleptcut = 100MeV
    cuts.AddCut("50200013", "00200029007000008250400000"); //singleptcut = 150MeV
  } else if (trainConfig == 61) {
    cuts.AddCut("50200013", "00200009007000009250400000"); // Qtmax = 0.03
    cuts.AddCut("50200013", "00200009007000002250400000"); // Qtmax = 0.06
    cuts.AddCut("50200013", "00200009007000006250400000"); // Qtmax = 0.02
    cuts.AddCut("50200013", "00200009007000008750400000"); //chi2cut = 10
    cuts.AddCut("50200013", "00200009007000008150400000"); //chi2cut = 50
  } else if (trainConfig == 62) {
    cuts.AddCut("50200013", "00200009007000008260400000"); //psipaircut = 0.05
    cuts.AddCut("50200013", "00200009007000008280400000"); //psipaircut = 0.2
    cuts.AddCut("50200013", "00200009007000008270400000"); //psipaircut = 0.07
    cuts.AddCut("50200013", "00200009007000008250600000"); //cos p angle cut = 0.9
    cuts.AddCut("50200013", "00200009007000008250300000"); //cos p angle cut = 0.75
  //systematics 20-40%
  } else if (trainConfig == 63) {
    cuts.AddCut("52400013", "04200009007000008250400000"); //eta cut: |eta| <0.75
    cuts.AddCut("52400013", "00500009007000008250400000"); //minR= 10 maxR = 180
    cuts.AddCut("52400013", "00200049007000008250400000"); //singleptcut = 75MeV
    cuts.AddCut("52400013", "00200019007000008250400000"); //singleptcut = 100MeV
    cuts.AddCut("52400013", "00200029007000008250400000"); //singleptcut = 150MeV
  } else if (trainConfig == 64) {
    cuts.AddCut("52400013", "00200009007000009250400000"); // Qtmax = 0.03
    cuts.AddCut("52400013", "00200009007000002250400000"); // Qtmax = 0.06
    cuts.AddCut("52400013", "00200009007000006250400000"); // Qtmax = 0.02
    cuts.AddCut("52400013", "00200009007000008750400000"); //chi2cut = 10
    cuts.AddCut("52400013", "00200009007000008150400000"); //chi2cut = 50
  } else if (trainConfig == 65) {
    cuts.AddCut("52400013", "00200009007000008260400000"); //psipaircut = 0.05
    cuts.AddCut("52400013", "00200009007000008280400000"); //psipaircut = 0.2
    cuts.AddCut("52400013", "00200009007000008270400000"); //psipaircut = 0.07
    cuts.AddCut("52400013", "00200009007000008250600000"); //cos p angle cut = 0.9
    cuts.AddCut("52400013", "00200009007000008250300000"); //cos p angle cut = 0.75
  //systematics 40-80%
  } else if (trainConfig == 66) {
    cuts.AddCut("54800013", "04200009007000008250400000"); //eta cut: |eta| <0.75
    cuts.AddCut("54800013", "00500009007000008250400000"); //minR= 10 maxR = 180
    cuts.AddCut("54800013", "00200049007000008250400000"); //singleptcut = 75MeV
    cuts.AddCut("54800013", "00200019007000008250400000"); //singleptcut = 100MeV
    cuts.AddCut("54800013", "00200029007000008250400000"); //singleptcut = 150MeV
  } else if (trainConfig == 67) {
    cuts.AddCut("54800013", "00200009007000009250400000"); // Qtmax = 0.03
    cuts.AddCut("54800013", "00200009007000002250400000"); // Qtmax = 0.06
    cuts.AddCut("54800013", "00200009007000006250400000"); // Qtmax = 0.02
    cuts.AddCut("54800013", "00200009007000008750400000"); //chi2cut = 10
    cuts.AddCut("54800013", "00200009007000008150400000"); //chi2cut = 50
  } else if (trainConfig == 68) {
    cuts.AddCut("54800013", "00200009007000008260400000"); //psipaircut = 0.05
    cuts.AddCut("54800013", "00200009007000008280400000"); //psipaircut = 0.2
    cuts.AddCut("54800013", "00200009007000008270400000"); //psipaircut = 0.07
    cuts.AddCut("54800013", "00200009007000008250600000"); //cos p angle cut = 0.9
    cuts.AddCut("54800013", "00200009007000008250300000"); //cos p angle cut = 0.75

  } else if (trainConfig == 69) { //for dir. photon purity studies
    cuts.AddCut("50100013", "00200009007000008250400000"); //chi2cut = 10
  } else if (trainConfig == 70) { //for dir. photon purity studies
    cuts.AddCut("50100013", "00200009007000008750400000"); //chi2cut = 10
  } else if (trainConfig == 71) { //for dir. photon purity studies
    cuts.AddCut("52400013", "00200009007000008750400000"); //chi2cut = 10

  //Additional systematics 0-20% & 20-40%
  } else if (trainConfig == 81) {
    cuts.AddCut("50200013", "00300009007000008250400000"); // 10 < R < 70
    cuts.AddCut("50200013", "00400009007000008250400000"); // 5 < R < 70
    cuts.AddCut("52400013", "00300009007000008250400000"); // 10 < R < 70
    cuts.AddCut("52400013", "00400009007000008250400000"); // 5 < R < 70
  } else if (trainConfig == 82) {
    cuts.AddCut("50200013", "05200009007000008250400000"); // |eta| < 0.5
    cuts.AddCut("50200013", "08200009007000008250400000"); // |eta| < 0.4
    cuts.AddCut("52400013", "05200009007000008250400000"); // |eta| < 0.5
    cuts.AddCut("52400013", "08200009007000008250400000"); // |eta| < 0.4
  //Additional checks purity with semi-closed TPC cuts
  } else if (trainConfig == 83) {
    cuts.AddCut("50200013", "00200009207000008250400000"); // TPCe -5,5
    cuts.AddCut("52400013", "00200009207000008250400000"); // TPCe -5,5
    cuts.AddCut("50200013", "00200009407000008250400000"); // TPCe -6,7
    cuts.AddCut("52400013", "00200009407000008250400000"); // TPCe -6,7
  } else if (trainConfig == 84) {
    cuts.AddCut("50200013", "00200009017000008250400000"); // TPCpi 0,-10
    cuts.AddCut("52400013", "00200009017000008250400000"); // TPCpi 0,-10
    cuts.AddCut("50200013", "00200009097000008250400000"); // TPCpi 1,0.5
    cuts.AddCut("52400013", "00200009097000008250400000"); // TPCpi 1,0.5
  } else if (trainConfig == 85) {
    cuts.AddCut("50200013", "00200009037000008250400000"); // TPCpi 2,5,-10
    cuts.AddCut("52400013", "00200009037000008250400000"); // TPCpi 2,5,-10
    cuts.AddCut("50200013", "00200009047000008250400000"); // TPCpi 3,1
    cuts.AddCut("52400013", "00200009047000008250400000"); // TPCpi 3,1
  } else if (trainConfig == 86) {
    cuts.AddCut("50200013", "00200009067000008250400000"); // TPCpi 2, 0,5
    cuts.AddCut("52400013", "00200009067000008250400000"); // TPCpi 2, 0,5
    cuts.AddCut("50200013", "00200009087000008250400000"); // TPCpi 2, 1
    cuts.AddCut("52400013", "00200009087000008250400000"); // TPCpi 2, 1
  } else if (trainConfig == 87) {
    cuts.AddCut("50200013", "00200009417000008250400000"); // TPCpi 0,-10, TPCe -6,7
    cuts.AddCut("52400013", "00200009417000008250400000"); // TPCpi 0,-10, TPCe -6,7
    cuts.AddCut("50200013", "00200009497000008250400000"); // TPCpi 1,0.5, TPCe -6,7
    cuts.AddCut("52400013", "00200009497000008250400000"); // TPCpi 1,0.5, TPCe -6,7
  } else if (trainConfig == 88) {
    cuts.AddCut("50200013", "00200009437000008250400000"); // TPCpi 2,5,-10, TPCe -6,7
    cuts.AddCut("52400013", "00200009437000008250400000"); // TPCpi 2,5,-10, TPCe -6,7
    cuts.AddCut("50200013", "00200009447000008250400000"); // TPCpi 3, 1     TPCe -6,7
    cuts.AddCut("52400013", "00200009447000008250400000"); // TPCpi 3, 1     TPCe -6,7
  } else if (trainConfig == 89) {
    cuts.AddCut("50200013", "00200009467000008250400000"); // TPCpi 2, 0,5, TPCe -6,7
    cuts.AddCut("52400013", "00200009467000008250400000"); // TPCpi 2, 0,5, TPCe -6,7
    cuts.AddCut("50200013", "00200009487000008250400000"); // TPCpi 2,1, TPCe -6,7
    cuts.AddCut("52400013", "00200009487000008250400000"); // TPCpi 2,1, TPCe -6,7
  } else {
      Error(Form("GammaConvV1_%i",trainConfig), "wrong trainConfig variable no cuts have been specified for the configuration");
      return;
  }

  if(!cuts.AreValid()){
    cout << "\n\n****************************************************" << endl;
    cout << "ERROR: No valid cuts stored in CutHandlerConv! Returning..." << endl;
    cout << "****************************************************\n\n" << endl;
    return;
  }

  const int numberOfCuts = cuts.GetNCuts();

    //================================================
  //========= Add task to the ANALYSIS manager =====
  //================================================
  AliAnalysisTaskGammaConvFlow *task=NULL;
  task= new AliAnalysisTaskGammaConvFlow(Form("GammaConvV1_%i_v%d",trainConfig,harmonic),numberOfCuts);
  task->SetIsHeavyIon(isHeavyIon);
  task->AliAnalysisTaskGammaConvFlow::SetIsMC(isMC);
  task->SetV0ReaderName(V0ReaderName);

  cutsRP = new AliFlowTrackCuts(Form("RFPcuts%s",uniqueID));
  if(!cutsRP) {
      if(debug) cout << " Fatal error: no RP cuts found, could be a library problem! " << endl;
      return 0x0;
  }
  cutsRP = cutsRP->GetStandardVZEROOnlyTrackCuts(); // select vzero tracks
  cutsRP->SetVZEROgainEqualizationPerRing(kFALSE);
  cutsRP->SetApplyRecentering(kTRUE);

  task->SetRPCuts(cutsRP);

  if(UseMassSel==kTRUE)  task->SetMassWindow(MinMass,MaxMass);
  if(UseKappaSel==kTRUE) task->SetKappaWindow(MinKappa,MaxKappa);

  task->SetPerformExtraStudies(PerformExtraStudies);
  task->SetApplydPhidRCut(ApplydPhidRCut);

  AliAnalysisDataContainer *flowEvent[numberOfCuts];
  //======================================================================
  for(Int_t i = 0; i<numberOfCuts; i++){
    if(debug) cout << " === RECEIVED REQUEST FOR FLOW ANALYSIS === " << endl;
    TString totalCutString = Form("%s_%s",(cuts.GetEventCut(i)).Data(),(cuts.GetPhotonCut(i)).Data());
    flowEvent[i] = mgr->CreateContainer(Form("FlowContainer_%s_%s",uniqueID.Data(),totalCutString.Data()), AliFlowEventSimple::Class(), AliAnalysisManager::kExchangeContainer);
    mgr->ConnectOutput(task, 2+i, flowEvent[i]);

    Double_t NFilterBinValues[20];
    if(FilterVariable==7 && NFilterBins == 4){
      task->SetFilterVariable(2,MinFilter,MaxFilter);
      NFilterBinValues[0] = -20;
      NFilterBinValues[1] = -13;
      NFilterBinValues[2] = -11;
      NFilterBinValues[3] = -6;
      NFilterBinValues[4] = -3;
      NFilterBinValues[5] = 5;
      NFilterBinValues[6] = 11;
      NFilterBinValues[7] = 20;
    } else if(FilterVariable==8 && NFilterBins == 4){
      task->SetFilterVariable(2,MinFilter,MaxFilter);
      NFilterBinValues[0] = -20;
      NFilterBinValues[1] = -13;
      NFilterBinValues[2] = -11;
      NFilterBinValues[3] = -6;
      NFilterBinValues[4] = 0;
      NFilterBinValues[5] = 5;
      NFilterBinValues[6] = 11;
      NFilterBinValues[7] = 20;
    }else{
      task->SetFilterVariable(FilterVariable,MinFilter,MaxFilter);
      if(NFilterBins > 1){
        for(Int_t k=0;k<NFilterBins+1;k++){
          NFilterBinValues[k] = MinFilter + k*(MaxFilter-MinFilter)/NFilterBins;
        }
      }else{
        NFilterBinValues[0] = MinFilter;
        NFilterBinValues[1] = MaxFilter;
      }
    }

    for(Int_t j=0;j<NFilterBins;j++){
      AliFlowTrackSimpleCuts *POIfilterVZERO = new AliFlowTrackSimpleCuts();
      if(FilterVariable==7 || FilterVariable==8){
        POIfilterVZERO->SetMassMin(NFilterBinValues[2*j]); POIfilterVZERO->SetMassMax(NFilterBinValues[2*j+1]);
      }else{
        POIfilterVZERO->SetMassMin(NFilterBinValues[j]); POIfilterVZERO->SetMassMax(NFilterBinValues[j+1]);
      }

      if(debug) cout << "    --> Created IO containers " << flowEvent[i] << endl;
      AddSPmethod(Form("SPVZEROQa_in_%s_%i", uniqueID.Data(), j), "Qa", harmonic, flowEvent[i], debug,uniqueID, POIfilterVZERO, trainConfig,BasicHistoSP,totalCutString);
      if(debug) cout << "    --> Hanging SP Qa task ... succes!" << endl;
      AddSPmethod(Form("SPVZEROQb_in_%s_%i", uniqueID.Data(), j), "Qb", harmonic, flowEvent[i], debug,uniqueID, POIfilterVZERO, trainConfig,BasicHistoSP,totalCutString);
      if(debug) cout << "    --> Hanging SP Qb task ... succes!"<< endl;
    }
  }

  //======================================================================
  TList *EventCutList = new TList();
  TList *ConvCutList = new TList();

  TList *HeaderList = new TList();
  TObjString *Header1 = new TObjString("BOX");
  HeaderList->Add(Header1);
  //    TObjString *Header3 = new TObjString("eta_2");
  //    HeaderList->Add(Header3);

  EventCutList->SetOwner(kTRUE);
  AliConvEventCuts **analysisEventCuts = new AliConvEventCuts*[numberOfCuts];
  ConvCutList->SetOwner(kTRUE);
  AliConversionPhotonCuts **analysisCuts = new AliConversionPhotonCuts*[numberOfCuts];

  for(Int_t i = 0; i<numberOfCuts; i++){
    analysisEventCuts[i] = new AliConvEventCuts();
    analysisEventCuts[i]->InitializeCutsFromCutString((cuts.GetEventCut(i)).Data());
    analysisEventCuts[i]->SetV0ReaderName(V0ReaderName);
    EventCutList->Add(analysisEventCuts[i]);
    analysisEventCuts[i]->SetFillCutHistograms("",kFALSE);

    analysisCuts[i] = new AliConversionPhotonCuts();
    analysisCuts[i]->InitializeCutsFromCutString((cuts.GetPhotonCut(i)).Data());
    analysisCuts[i]->SetV0ReaderName(V0ReaderName);
    ConvCutList->Add(analysisCuts[i]);
    analysisCuts[i]->SetFillCutHistograms("",kFALSE);

    analysisEventCuts[i]->SetAcceptedHeader(HeaderList);
  }

  task->SetEventCutList(numberOfCuts,EventCutList);
  task->SetConversionCutList(numberOfCuts,ConvCutList);
//    task->SetMesonCutList(numberOfCuts,MesonCutList);
//    task->SetMoveParticleAccordingToVertex(kTRUE);
  task->SetDoMesonAnalysis(kFALSE);
  task->SetDoMesonQA(enableQAMesonTask); //Attention new switch for Pi0 QA
  task->SetDoPhotonQA(enableQAPhotonTask);  //Attention new switch small for Photon QA

  //connect containers
  AliAnalysisDataContainer *coutput =
  mgr->CreateContainer(Form("GammaConvV1_%i_v%d",trainConfig,harmonic), TList::Class(),
                      AliAnalysisManager::kOutputContainer,Form("GammaConvFlow_%i.root",trainConfig));

  mgr->AddTask(task);
  mgr->ConnectInput(task,0,cinput);
  mgr->ConnectOutput(task,1,coutput);

  return;

}

//_____________________________________________________________________________
void AddSPmethod(char *name, char *Qvector, int harmonic, AliAnalysisDataContainer *flowEvent, bool debug, TString uniqueID,  AliFlowTrackSimpleCuts* POIfilter,Int_t trainConfig, bool BasicHistoSP = kTRUE, TString CutNumberString)
{
  // add sp task and invm filter tasks
  if(debug)  cout << " ******* Switching to SP task ******* " << endl;
  TString fileName = Form("GammaConvFlow_%i.root:SP_V0_%s",trainConfig, CutNumberString.Data());
  if(debug) cout << "    --> fileName " << fileName << endl;
  TString myFolder = fileName;
  if(debug) cout << "    --> myFolder " << myFolder << endl;
  TString myNameSP;
  myNameSP = Form("%sSPv%d%s_%s", name, harmonic, Qvector,CutNumberString.Data());
  if(debug) cout << " myNameSP " << myNameSP << endl;
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  AliAnalysisDataContainer *flowEventOut = mgr->CreateContainer(Form("Filter_%s",myNameSP.Data()),AliFlowEventSimple::Class(),AliAnalysisManager::kExchangeContainer);
  AliAnalysisTaskFilterFE *tskFilter = new AliAnalysisTaskFilterFE(Form("TaskFilter_%s", myNameSP.Data()), NULL, POIfilter);
  // tskFilter->SetSubeventEtaRange(minEtaA, maxEtaA, minEtaB, maxEtaB);
  tskFilter->SetSubeventEtaRange(-10, -1, 1, 10);
  mgr->AddTask(tskFilter);
  mgr->ConnectInput(tskFilter, 0, flowEvent);
  mgr->ConnectOutput(tskFilter, 1, flowEventOut);
  AliAnalysisDataContainer *outSP = mgr->CreateContainer(myNameSP.Data(), TList::Class(), AliAnalysisManager::kOutputContainer, fileName);
  AliAnalysisTaskScalarProduct *tskSP = new AliAnalysisTaskScalarProduct(Form("TaskScalarProduct_%s", myNameSP.Data()), kFALSE);
  tskSP->SetApplyCorrectionForNUA(kTRUE);
  tskSP->SetHarmonic(harmonic);
  tskSP->SetTotalQvector(Qvector);
  //if (bEP) tskSP->SetBehaveAsEP();
  tskSP->SetBookOnlyBasicCCH(BasicHistoSP);
  mgr->AddTask(tskSP);
  mgr->ConnectInput(tskSP, 0, flowEventOut);
  mgr->ConnectOutput(tskSP, 1, outSP);
}
//_____________________________________________________________________________





