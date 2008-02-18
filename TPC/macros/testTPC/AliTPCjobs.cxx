//#include <fstream>
//#include <TFile.h>
//#include <TSystem.h>


/*

  Process all jobs specified in the job.list file
  Concurent agents acces the list
  Once starting to procces given job the jonb is locked 

*/



/*


.L $ALICE_ROOT/TPC/macros/testTPC/AliTPCjobs.cxx+
 

AliTPCJobs jobs;
jobs.fJobFile="job.list"
jobs.ProcessAllJobs();
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
  ClassDef(AliTPCJobs,0)
};

 ClassImp(AliTPCJobs)

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
    printf("PROCESSING JOB\t%d\n",counter);
    counter++;
    if (!GetNextJob()) break;   
  }
}


Bool_t AliTPCJobs::GetNextJob(){
  //
  // GetNextJob  - get job from the list which is not locked
  //
  ifstream ins;
  ins.open(fJobFile);
  TString id;
  TString inputData;
  TString outputDir;
  TString action;
  Bool_t hasJob=kFALSE;
  while(ins.good()) {
    ins>>id;
    ins>>inputData;
    ins>>outputDir;
    ins>>action;
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
      sprintf(command,"xrdcp  -DIFirstConnectMaxCnt 4 -DIConnectTimeout 4 -DIRequestTimeout 4 -DIMaxRedirectcount 4 -DIRedirCntTimeout 4 %s\t%s\n",inputData.Data(), outputDir.Data());
    printf("Exec\t%s\n", command);
    gSystem->Exec(command);
    //TFile::Cp(inputData.Data(), outputDir.Data());
  }else{
    char command[10000];
    sprintf(command,"$ALICE_ROOT/TPC/macros/testTPC/action.sh %s %s %s %s", jobID.Data(), inputData.Data(), outputDir.Data(), action.Data());
    printf("%s\n\n",command);
    gSystem->Exec(command);
    printf("\n\n");

  }
  
  SetDone(jobID);
}



void AliTPCjobs(){
  //
  //
  //
  AliTPCJobs jobs;
  jobs.fJobFile="job.list";
  jobs.ProcessAllJobs();
}
