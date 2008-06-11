#include <fstream>
#include <TFile.h>
#include <TSystem.h>
#include <TMap.h>

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
  void GetListOfJobs();

  void ProcessAllJobs();
  Bool_t GetNextJob();
  Bool_t GetNextJobOld();
  void ProcessJob(TString jobID, TString inputData, TString outputDir, TString   action);
  
  void    SetLock(TString jobID);
  void    SetDone(TString jobID);
  void    SetFail(TString jobID);

  void    Stage(TString inputData,TString jobID);

  Bool_t  IsLocked(TString jobID);
  Bool_t  IsFail(TString jobID);
  Bool_t  IsStaged(TString inputData);
  Bool_t  IsCastorFile(TString filename);
  //
  //
  TString  fJobFile;
  TString  fWorkDir;
  TObjArray fJobInfo;  // array of job infos

  ClassDef(AliTPCJobs,0)
};

 ClassImp(AliTPCJobs)

AliTPCJobs::AliTPCJobs(){
  // 
  // load neccessary libraries 
  //
  gSystem->Load("libXrdClient.so");
  gSystem->Load("libNetx.so");
}

void AliTPCJobs::GetListOfJobs(){
  //
  // char * fJobFile ="job.list"
  //
  //  Update the list of jobs - using the "local file"
  //
  ifstream ins;
  ins.open(fJobFile);
  char line[10000];
  TObjArray *arrayList=0;
  while(ins.good()) {
    ins.getline(line,10000);
    //printf("%s\n",line);     
    arrayList = (TString(line)).Tokenize(" ");
    if (arrayList->GetEntries()>3){
      TObjString * jobID = (TObjString*)arrayList->At(0);
      TMap *map =  (TMap*)fJobInfo.FindObject(jobID->GetName());
      if (!map) {
	map = new TMap();
	map->SetName(jobID->GetName());
	fJobInfo.AddLast(map);
	map->Add(new TObjString("Input"),arrayList->At(1));
	map->Add(new TObjString("Output"),arrayList->At(2));
	map->Add(new TObjString("Action"),arrayList->At(3));
	map->Add(new TObjString("Status"),new TObjString("0"));	
      }
    }else{
      delete arrayList;
    }
  }
}


void AliTPCJobs::ProcessAllJobs(){
  //
  //
  //
  GetListOfJobs();
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
  TString id;
  TString inputData;
  TString outputDir;
  TString action;
  TString status;  
	
  Bool_t hasJob=kFALSE;
  Int_t entries  = fJobInfo.GetEntries();	
  TMap *map=0;
  for (Int_t i=0; i<entries; i++){
    map = (TMap*)fJobInfo.At(i);
    status=(*map)("Status")->GetName();
    id=map->GetName();
    inputData=(*map)("Input")->GetName();
    outputDir=(*map)("Output")->GetName();
    action=(*map)("Action")->GetName();
    if (status.Contains("1")) continue;
//    if (IsFail(id)) continue;
    if (!IsLocked(id)){
      hasJob=kTRUE;
      ((TObjString*)(*map)("Status"))->SetString("1");
      if (!IsStaged(inputData)){
        //Stage(inputData,id);
       continue;
      }
      break;
    }
  }

  printf("Process %s\n",id.Data());
  if (hasJob) {
    ProcessJob(id,inputData,outputDir, action);    
  }
  return hasJob;
}



Bool_t AliTPCJobs::GetNextJobOld(){
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
    //if (!inputData.Contains(".root")) SetFail(id);
    if (IsFail(id)) continue;
    if (!IsLocked(id)){
      hasJob=kTRUE;
      break;
    }
    if (!IsStaged(inputData)){
      //Stage(inputData,id);
      continue;
    }
  }
  printf("Process %s\n",id.Data());
  if (hasJob) ProcessJob(id,inputData,outputDir, action);
  return hasJob;
}


void    AliTPCJobs::SetLock(TString jobID){
  //
  // 
  //
  printf("touch out/%s.lock\n",jobID.Data());
  gSystem->Exec(Form("touch out/%s.lock",jobID.Data()));
}

void    AliTPCJobs::SetDone(TString jobID){ 
  //
  //
  //
  printf("touch out/%s.done\n",jobID.Data());
  gSystem->Exec(Form("touch out/%s.done", jobID.Data()));
}

void    AliTPCJobs::SetFail(TString jobID){
  //
  //
  // 
  printf("touch out/%s.fail\n",jobID.Data());
  gSystem->Exec(Form("touch out/%s.fail", jobID.Data()));
}

void  AliTPCJobs::Stage(TString inputData, TString jobID){
  //
  // stage file
  //
  inputData.ReplaceAll("root://voalice04.cern.ch:1094/","");
  if ( !IsCastorFile(inputData) ) return;
  char command[1000];
  sprintf(command,"stager_get -M %s \| grep SUBREQUEST_FAILED ",inputData.Data());
  FILE *pipe = gSystem->OpenPipe(command,"r");
  TString result;
  result.Gets(pipe);
  gSystem->ClosePipe(pipe);
  if ( result.Contains("SUBREQUEST_FAILED") ){
      SetFail(jobID);
  }
}

Bool_t    AliTPCJobs::IsLocked(TString jobID){
  TString path = "out/";
  path+=jobID;
  path+=".lock";
  Long_t pid; Long_t psize; Long_t pflags; Long_t pmodtime;
  Int_t status = gSystem->GetPathInfo(path,&pid,&psize,&pflags,&pmodtime);
  return (status==0);
}

Bool_t    AliTPCJobs::IsFail(TString jobID){
  //
  //
  //
  TString path = "out/";
  path+=jobID;
  path+=".fail";
  Long_t pid; Long_t psize; Long_t pflags; Long_t pmodtime;
  Int_t status = gSystem->GetPathInfo(path,&pid,&psize,&pflags,&pmodtime);
  return (status==0);
}

Bool_t  AliTPCJobs::IsStaged(TString inputData){
  //
  // check if file bname is staged
  //
  inputData.ReplaceAll("root://voalice04.cern.ch:1094/","");
  if ( !IsCastorFile(inputData) ) return kTRUE;
  char command[1000];
  sprintf(command,"stager_qry -M %s \| grep /castor \| gawk  \'{ print $3;}\'",inputData.Data());
  FILE *pipe = gSystem->OpenPipe(command,"r");
  TString result;
  result.Gets(pipe);
  gSystem->ClosePipe(pipe);
  if ( result.Contains("STAGED") ) return kTRUE;

  return kFALSE;
}

Bool_t  AliTPCJobs::IsCastorFile(TString filename){
  //
  // check if filename begins with 'castor'
  //
  return filename.BeginsWith("/castor/");
}


void AliTPCJobs::ProcessJob(TString jobID, TString inputData, TString outputDir, TString   action){
  //
  //
  // 1. Create lock file
  // 2. Get Input data
  // 3. Process data
  // 4. Create Done file
  SetLock(jobID);
  if (action.Contains("TEST")){
    char command[10000];
    sprintf(command,"$ALICE_ROOT/TPC/macros/testTPC/action.sh %s %s %s %s", jobID.Data(), inputData.Data(), outputDir.Data(), action.Data());    
    printf("%s\n\n",command); 
    SetDone(jobID);
    return;
  }
  if (action.Contains("COPY")){
    
    char command[10000];
    sprintf(command,"mkdirhier %s",gSystem->DirName(outputDir.Data()));
    printf("Exec\t%s\n", command);
    gSystem->Exec(command);
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

    // create temp work dir
//    TString processDir=fWorkDir+jobID+action;
//    Int_t res = gSystem->mkdir(processDir,kTRUE);
//    if ( res==-1 ){
//	AliWarning(Form("Cannot create dir/already exists: '%s'",processDir.Data()));
//	return;
//    }
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
