#include "AliAnalysisDataContainer.h"
#include "AliGFWWeights.h"
#include "GFWFlags.h"
class TNamed;
using namespace GFWFlags;
TList *constructList(TFile *infile, TString subfx) {
  if(!infile) {printf("Input file is not open!\n"); return 0; };
  if(infile->IsZombie()) {printf("Input file is a zombie!\n"); return 0; };
  TList *retList = new TList();
  retList->SetOwner(kTRUE);
  TList *LoK = infile->GetListOfKeys();
  UInt_t lEvFlg, lTrFlg;
  AliAnalysisTaskGFWFlow::SetupFlagsByIndex(subfx.Atoi(),lEvFlg,lTrFlg);
  TString syspf = GetSystPF(BitIndex(lEvFlg),BitIndex(lTrFlg));
  for(Int_t i=0; i<LoK->GetEntries(); i++) {
    TString nm = LoK->At(i)->GetName();
    TString cmpsyspf=syspf+";1";
    if(!nm.Contains(syspf.Data())) continue;
    AliGFWWeights *w = (AliGFWWeights*)infile->Get(LoK->At(i)->GetName());
    retList->Add(w);
  };
  return retList;
}
AliAnalysisTaskGFWFlow* AddTaskGFWFlow(TString name = "name", Bool_t ProduceWeights=kFALSE, Bool_t IsMC=kFALSE, Bool_t IsTrain=kTRUE, TString weightpath="", TString centMap="", TString subfx="")
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) return 0x0;
  if (!mgr->GetInputEventHandler())	return 0x0;
  TString fileName = AliAnalysisManager::GetCommonFileName();
  AliAnalysisTaskGFWFlow* task = new AliAnalysisTaskGFWFlow(Form("%s%s",name.Data(),subfx.Data()), ProduceWeights, IsMC, IsTrain);
  if(!task)
    return 0x0;
  //My settings:
  mgr->AddTask(task); // add your task to the manager
  //Connect weights to a container
  if(!ProduceWeights) {
    if(IsTrain) {
      if(centMap.IsNull()) return 0; //AliFatal("Centrality map not specified!\n");
      if(centMap.Contains("alien:")) TGrid::Connect("alien:");
      TFile *tf = TFile::Open(centMap.Data());
      TH1D *cmap = (TH1D*)tf->Get("AMPT_Cent_Map")->Clone("AMPT_Cent_Map");
      cmap->SetDirectory(0);
      if(!cmap) return 0;//AliFatal("Could not find AMPT_Cent_Map in file specified!\n");
      AliAnalysisDataContainer *cInMap = mgr->CreateContainer("CentralityMap",TH1D::Class(),AliAnalysisManager::kInputContainer);
      cInMap->SetData(cmap);
      mgr->ConnectInput(task,1,cInMap);
      printf("Centrality map set!\n");
    } else {
      TObjArray *AllContainers = mgr->GetContainers();
      TString iwName = "InputWeights"+subfx;
      if(!AllContainers->FindObject(iwName.Data())) {
        printf("InputWeights not loaded yet, loading now!\n");
        if(weightpath.EqualTo("")) { printf("Weight path for containers not set!\n"); return NULL; };
        if(weightpath.Contains("alien:")) TGrid::Connect("alien:");
        TFile *tf = TFile::Open(weightpath.Data());
        TList *tl = constructList(tf,subfx);
        if(!tl) return 0;
        AliAnalysisDataContainer *cInWeights = mgr->CreateContainer(iwName.Data(),TList::Class(), AliAnalysisManager::kInputContainer);
        cInWeights->SetData(tl);
        mgr->ConnectInput(task,1,cInWeights);
      } else {
        mgr->ConnectInput(task,1,(AliAnalysisDataContainer*)AllContainers->FindObject(iwName.Data()));
      };
    };
  };

  AliAnalysisDataContainer* cInput0 = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer* cOutput1;
  if(ProduceWeights)
    cOutput1 = mgr->CreateContainer(Form("OutputList%s",subfx.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, mgr->GetCommonFileName());
  else cOutput1 = mgr->CreateContainer(Form("OutCont%s",subfx.Data()), AliGFWFlowContainer::Class(), AliAnalysisManager::kOutputContainer, mgr->GetCommonFileName());
  // Connecting containers to task
  mgr->ConnectInput(task,0,cInput0); // your task needs input: here we connect the manager to your task
  mgr->ConnectOutput(task,1,cOutput1);
  AliAnalysisDataContainer *multidist = mgr->CreateContainer(Form("MultiDist%s",subfx.Data()), TH1D::Class(), AliAnalysisManager::kOutputContainer, mgr->GetCommonFileName());
  mgr->ConnectOutput(task,2,multidist);
  return task;
}
