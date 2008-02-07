#include <fstream>
#include <TFile.h>
#include <TSystem.h>

/*
.L $ALICE_ROOT/TPC/macros/testTPC/AliTPCjobs.cxx+
AliTPCJobs jobs;
jobs.fJobFile="job.list"
*/
class AliTPCJobs : public TNamed{
public:
  AliTPCJobs();
  void ProcessAllJobs();
  Bool_t GetNextJob();
  void ProcessJob(TString jobID, TString inputData, TString outputDir, TString   action);
  
  void    SetLock(TString jobID);
  void    SetDone(TString jobID);
  Bool_t  IsLocked(TString jobID);

  TString  fJobFile;    
  TString  fWorkDir;
};


AliTPCJobs::AliTPCJobs(){
  // 
  //
  //
  gSystem->Load("libXrdClient.so");
  gSystem->Load("libNetx.so");
}

void AliTPCJobs::ProcessAllJobs(){
  //
  //
  //
  Int_t counter=0;
  while (GetNextJob()){
    //
    printf("PROCESSING JOB\n",counter);
    counter++;
    if (!GetNextJob()) break;   
  }
}


Bool_t AliTPCJobs::GetNextJob(){
  ifstream in;
  in.open(fJobFile);
  TString id;
  TString inputData;
  TString outputDir;
  TString action;
  Bool_t hasJob=kFALSE;
  while(in.good()) {
    in>>id;
    in>>inputData;
    in>>outputDir;
    in>>action;
    if (!IsLocked(id)){
      hasJob=kTRUE;
      break;
    }
  }
  printf("Process %s\n",id.Data());
  if (hasJob) ProcessJob(id,inputData,outputDir, action);
  return hasJob;
}


void    AliTPCJobs::SetLock(TString jobID){
  printf("touch out/%s.lock\n",jobID.Data());
  gSystem->Exec(Form("touch out/%s.lock",jobID.Data()));
}

void    AliTPCJobs::SetDone(TString jobID){ 
  printf("touch out/%s.done\n",jobID.Data());
  gSystem->Exec(Form("touch out/%s.done", jobID.Data()));
}


Bool_t    AliTPCJobs::IsLocked(TString jobID){
  TString path = "out/";
  path+=jobID;
  path+=".lock";
  Long_t pid; Long_t psize; Long_t pflags; Long_t pmodtime;
  Int_t status = gSystem->GetPathInfo(path,&pid,&psize,&pflags,&pmodtime);
  return (status==0);
}


void AliTPCJobs::ProcessJob(TString jobID, TString inputData, TString outputDir, TString   action){
  //
  //
  // 1. Create lock file
  // 2. Get Input data
  // 3. Process data
  // 4. Create Done file
  SetLock(jobID);
  if (action.Contains("COPY")){
    char command[10000];
      sprintf(command,"xrdcp -d 1 -DIFirstConnectMaxCnt 2 -DIConnectTimeout 2 -DIRequestTimeout 2 -DIMaxRedirectcount 2 -DIRedirCntTimeout 2 %s\t%s\n",inputData.Data(), outputDir.Data());
    printf("Exec\t%s\n", command);
    gSystem->Exec(command);
    //TFile::Cp(inputData.Data(), outputDir.Data());
  }
  SetDone(jobID);
}
