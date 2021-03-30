AliAnalysisTaskSEPbPbCorrelationsJetV2 *AddTaskPbPbCorrelationsJetV2(TString centMethod = "V0M", Int_t FilterBit = 768, Int_t nTPCcls = 70, Double_t Vz_max = 10, TString sNameList = "default") 
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    printf("Error in adding AnalysisTaskMuonFlow: no Analysis Manager found!\n");
    return NULL;
  }

  //AliAnalysisTaskSEPbPbCorrelationsJetV2 *task = new AliAnalysisTaskSEPbPbCorrelationsJetV2("PbPbCorrelations_JetV2");
  AliAnalysisTaskSEPbPbCorrelationsJetV2 *task = new AliAnalysisTaskSEPbPbCorrelationsJetV2(Form("PbPbCorrelations_JetV2_%s",centMethod.Data()));

  // Set analysis cuts
  task->SetCentMethod(centMethod);
  task->SetRemovePileup(kTRUE);
  task->SetRemovePileup2(kTRUE);
  task->SetRemovePileup3(kTRUE);
  task->SetUseRes(kTRUE);
  task->SetAverageRes(kFALSE); // this is really needed
  task->SetPtOrder(kTRUE); 
 
  Double_t centLimits[] = {20.,30.,60.,90.};
  const Int_t nBinCent = sizeof(centLimits) / sizeof(Double_t) - 1;
  task->SetCentBinning(nBinCent, centLimits);

  printf("Centrality: ");
  for(Int_t ibin = 0; ibin <= nBinCent; ++ibin) printf("%.1f ",(Double_t)centLimits[ibin]);
  printf("%s\n",centMethod.Data());


  const Int_t nZvtxBins  = 1;
  Double_t vertexLimits[nZvtxBins+1] = {-1*Vz_max,Vz_max};
  task->SetZvtxBinning(nZvtxBins, vertexLimits);

  Double_t ptLimits[] = {0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0, 4.5, 5.0, 5.5, 6.0, 7.0, 8.0, 9.0, 10.0, 12.0, 14.0, 17.0, 20.0, 25.0, 30.0, 40.0, 50.0};
  const Int_t nBinPt = sizeof(ptLimits) / sizeof(Double_t) - 1;
  task->SetPtBinning(nBinPt, ptLimits);

  //Double_t trigPtLimits[] = {0.2, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 14.0, 20.0, 50.0};
  Double_t trigPtLimits[] = {0.2, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 12.0, 14.0, 20.0, 50.0};
  const Int_t nBinTrigPt = sizeof(trigPtLimits) / sizeof(Double_t) - 1;
  task->SetTrigPtBinning(nBinTrigPt, trigPtLimits);

  Double_t assocPtLimits[] = {0.5, 1., 1.5, 2., 3., 5.,200.};
  const Int_t nBinAssocPt = sizeof(assocPtLimits) / sizeof(Double_t) - 1;
  task->SetAssocPtBinning(nBinAssocPt, assocPtLimits);

  const Int_t nBinEta = 5;
  Double_t etaLimits[nBinEta+1] = {-4.5, -4., -3.5, -3.0, -2.5, -2.0};
  task->SetEtaBinning(nBinEta, etaLimits);

  task->SetFilterBit(FilterBit);
  task->SetNTPCcls(nTPCcls);

// Input File

  TGrid::Connect("alien://");
  TFile *foadb = TFile::Open("alien:///alice/cern.ch/user/s/sitang/EP_Cali_run2/calibV0HIR.root");
  TFile *fFileRes = TFile::Open("alien:///alice/cern.ch/user/s/sitang/EP_Cali_run2/qVectResMuons_V0A.root");

//  TFile *foadb    = TFile::Open("./Convert_EP_Calibration/calibV0HIR.root");
//  TFile *fFileRes = TFile::Open("./Convert_EP_Calibration/qVectResMuons_V0A.root");

  if(!foadb || !fFileRes)
  {
   printf("ERROR: Required files are not found!\n");
  }

  TList *list_contQ   = (TList*)foadb->Get("list_contQ");  list_contQ->SetName(Form("list_contQ_%s",sNameList.Data()));  
  TList *list_Res     = (TList*)fFileRes->Get("list_Res"); list_Res->SetName(Form("list_Res_%s",sNameList.Data()));

  if(!list_contQ || !list_Res)
  {
   printf("ERROR: No lists are not found!\n");
  }
  
  task->SetListContQ(list_contQ);
  task->SetListRes(list_Res);  

  mgr->AddTask(task);

 // create input container
  AliAnalysisDataContainer *cinput1 = mgr->CreateContainer(Form("contQ_%s",sNameList.Data()),
                                    TList::Class(),
                                    AliAnalysisManager::kInputContainer);
  cinput1->SetData(list_contQ);
  
  AliAnalysisDataContainer *cinput2 = mgr->CreateContainer(Form("Res_%s",sNameList.Data()),
                                    TList::Class(),
                                    AliAnalysisManager::kInputContainer);
  cinput2->SetData(list_Res);

 
  // create output container
  TString outputFileName = AliAnalysisManager::GetCommonFileName();
  AliAnalysisDataContainer *output = mgr->CreateContainer(Form("FlowIncHistos_%s_%s", centMethod.Data(),sNameList.Data()), TList::Class(), AliAnalysisManager::kOutputContainer,
							  //Form("%s:FlowInc_%s", AliAnalysisManager::GetCommonFileName(), centMethod));
							  outputFileName);

  AliAnalysisDataContainer *output1 = mgr->CreateContainer(Form("FlowJetHistos_%s_%s",centMethod.Data(),sNameList.Data()), TList::Class(), AliAnalysisManager::kOutputContainer,
							  //Form("%s:FlowJet_%s", AliAnalysisManager::GetCommonFileName(), centMethod));
							  outputFileName);



  // finaly connect input and output
  mgr->ConnectInput(task, 0,  mgr->GetCommonInputContainer());
  mgr->ConnectInput(task, 1, cinput1);
  mgr->ConnectInput(task, 2, cinput2);
  mgr->ConnectOutput(task, 1, output);
  mgr->ConnectOutput(task, 2, output1);
  
  return task;
}

