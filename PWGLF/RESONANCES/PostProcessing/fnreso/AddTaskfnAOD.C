
void AddTaskfnAOD(const char *suffix = "CMWchrg", Double_t fFilterBit = 32, Double_t kkshmasscut = 1.04, Double_t nsigtpcpion=2.0, Double_t nsigtofpion=3.0, Double_t nsigtpckaon=2.0, Double_t nsigtofkaon=3.0, Double_t dcaxypos=0.06, Double_t dcaxyneg=0.06, Double_t dcav0daugh=1.0, Double_t dcav0pv=0.3, Double_t cospa=0.97, Double_t lowrad=0.5, Double_t lifetime=15, Double_t pidpion=4, Float_t nCRcut=70.0, Float_t ratiocrfccut=0.8, Double_t chi2globalcut=36.0, Double_t chi2cut=36.0)

{
  // standard with task
  printf("===================================================================================\n");
  printf("\n                PID: Initialising AliAnalysisTaskCMWPU                             \n");
  printf("===================================================================================\n");

  TGrid::Connect("alien://");

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  TString     outfileName = AliAnalysisManager::GetCommonFileName();
  AliAnalysisDataContainer  *cinput = mgr->GetCommonInputContainer();  // AOD event


  TString list1OutName = outfileName;        // common outfile filename
  list1OutName        += ":Results";         // This directory contains result histograms

  Int_t gCentMin = 0;
  Int_t gCentMax = 100;

  TString TaskCMWPID;
  TaskCMWPID.Form("gTaskCMWCent%d_%d_%s", gCentMin, gCentMax, suffix);

  AliAnalysisTaskfnAOD *task_CMW = new AliAnalysisTaskfnAOD(TaskCMWPID);


  task_CMW->SetFilterBit(fFilterBit);
  task_CMW->Setkkshmasscut(kkshmasscut);
  task_CMW->SetPIDnsigtpcpion(nsigtpcpion);
  task_CMW->SetPIDnsigtofpion(nsigtofpion);
  task_CMW->SetPIDnsigtpckaon(nsigtpckaon);
  task_CMW->SetPIDnsigtofkaon(nsigtofkaon);
  task_CMW->Setdcaxyposneg(dcaxypos, dcaxyneg);
  task_CMW->Setdcav0daugh(dcav0daugh);
  task_CMW->Setdcav0pv(dcav0pv);
  task_CMW->SetCosPA(cospa);
  task_CMW->SetLowradius(lowrad);
  task_CMW->SetLT(lifetime);
  task_CMW->SetPIDpion(pidpion);
  task_CMW->SetPTC(nCRcut, ratiocrfccut, chi2globalcut, chi2cut);




  ///---> Now Pass data and containers to Analysis Object ----
  
  mgr->AddTask(task_CMW);                        // connect the task to the analysis manager
  mgr->ConnectInput(task_CMW, 0, cinput);        // give AOD event to my Task..!!


  AliAnalysisDataContainer  *cOutPut1;
  TString                  sMyOutName;
  sMyOutName.Form("SimpleTask_%s",suffix);
  
  cOutPut1 = (AliAnalysisDataContainer *) mgr->CreateContainer(sMyOutName,TList::Class(),AliAnalysisManager::kOutputContainer,list1OutName.Data());
  mgr->ConnectOutput(task_CMW, 1, cOutPut1);
  
 
  printf("\n\n ================> AddTaskCMW() Configured properly <==================\n\n");

  
}
