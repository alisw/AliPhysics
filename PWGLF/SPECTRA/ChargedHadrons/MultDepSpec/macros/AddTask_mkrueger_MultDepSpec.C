#if defined(__CLING__)
#include "AliMultDepSpecAnalysisTask.h"
#endif

AliMultDepSpecAnalysisTask* AddTask_mkrueger_MultDepSpec(TString controlstring, Int_t cutModeLow = 100, Int_t cutModeHigh = 121)
{

  // settings:
  UInt_t offlineTriggerMask = AliVEvent::kINT7;
  Bool_t isPbPbAnalysis = kFALSE;

  Int_t maxMultiplicity = 100;
  Float_t etaCut = 0.8;

  Float_t lowerPtCut = 0.15;
  Float_t upperPtCut = 50.;

  string colsys = "pp";

  const Int_t nMultSteps = 3;
  Int_t multSteps[] = {100, 0, 0};
  Int_t multBinWidth[] = {1,1,1};

  Int_t nBinsCent = 1;
  Double_t centBinEdgesDummy[2] = {0., 100.};
  Double_t centBinEdgesUsed[9] = {0., 5., 10., 20., 40., 60., 80., 90., 100.};
  Double_t* centBinEdges = centBinEdgesDummy;

  if(controlstring.Contains("pp")){
    if(controlstring.Contains("5TeV"));
    if(controlstring.Contains("7TeV")) offlineTriggerMask = AliVEvent::kMB;
    if(controlstring.Contains("2TeV")) offlineTriggerMask = AliVEvent::kMB;
    if(controlstring.Contains("09TeV")) offlineTriggerMask = AliVEvent::kMB;
  }
  if(controlstring.Contains("XeXe"))  {
    colsys = "XeXe";
    multSteps[0] = 3500;  multBinWidth[0] = 1;
    multSteps[1] = 0;    multBinWidth[1] = 1;
    multSteps[2] = 0;    multBinWidth[2] = 1;
  }
  if(controlstring.Contains("pPb"))  {
    colsys = "pPb";
    multSteps[0] = 200;   multBinWidth[0] = 1;
    multSteps[1] = 0;     multBinWidth[1] = 1;
    multSteps[2] = 0;     multBinWidth[2] = 1;
  }
  if(controlstring.Contains("PbPb")) {
    isPbPbAnalysis = kTRUE;
    colsys = "PbPb";
    multSteps[0] = 4500;   multBinWidth[0] = 1;
    multSteps[1] = 0;    multBinWidth[1] = 1;
    multSteps[2] = 0;    multBinWidth[2] = 1;
    if(controlstring.Contains("2TeV")){offlineTriggerMask = AliVEvent::kMB;}
  }


  // Binning in Multiplicity
  const Int_t nBinsMult = multSteps[0]+ multSteps[1] + multSteps[2] + 1;
  Double_t multBinEdges[nBinsMult+1];
  multBinEdges[0] = -0.5;
  multBinEdges[1] = 0.5;
  Int_t startBin = 1;
  Int_t endBin = 1;
  for(Int_t multStep = 0; multStep < nMultSteps; multStep++){
    endBin += multSteps[multStep];
    for(Int_t multBin = startBin; multBin < endBin; multBin++)  {
      multBinEdges[multBin+1] = multBinEdges[multBin] + multBinWidth[multStep];
    }
    startBin = endBin;
  }


  AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();

  if (!mgr) {
    Error("AddTask_mkrueger_MultDepSpec", "No analysis manager found.");
    return 0;
  }

  // Switch off all AliInfo (too much output!!!)
  AliLog::SetGlobalLogLevel(AliLog::kError);
  mgr->SetDebugLevel(0);

  TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
  Bool_t hasMC = (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler() != 0x0);

  if(controlstring.Contains("PbPb") && hasMC) offlineTriggerMask = AliVEvent::kMB;
  if(controlstring.Contains("highMult")){offlineTriggerMask = AliVEvent::kHighMult;}

  AliMultDepSpecAnalysisTask* mainTask = NULL;

  char taskName[100] = "";

  for(Int_t cutMode = cutModeLow; cutMode <= cutModeHigh; cutMode++){
    sprintf(taskName, "mkrueger_%s_eta_%.2f_cutMode_%d", colsys.c_str(), etaCut, cutMode);

    AliMultDepSpecAnalysisTask* task = new AliMultDepSpecAnalysisTask(taskName);

    task->SetCutMode(cutMode);
    // general cuts
    task->SelectCollisionCandidates(offlineTriggerMask);
    task->SetTriggerMask(offlineTriggerMask);

    task->SetUseMC(hasMC);
    if(type.Contains("ESD")) task->SetUseESD();
    else task->SetUseAOD();
    task->SetBinsMult(nBinsMult, multBinEdges);
    task->SetBinsCent(nBinsCent, centBinEdges);

    // nominal cut-setting:
    task->SetMinEta(-etaCut);
    task->SetMaxEta(etaCut);
    task->SetMinPt(lowerPtCut);
    task->SetMaxPt(upperPtCut);
    // todo add vertex cut


    if(cutMode == cutModeLow) mainTask = task; // has no particular use... just to return a task

  // hang task in train

    mgr->AddTask(task);

    AliAnalysisDataContainer* cinput = mgr->GetCommonInputContainer();
    AliAnalysisDataContainer* coutput = mgr->CreateContainer(taskName,
							    TList::Class(),
							    AliAnalysisManager::kOutputContainer,
							    "AnalysisResults.root");



    mgr->ConnectInput(task, 0, cinput);
    mgr->ConnectOutput(task, 1, coutput);


  }


  return mainTask;

}
