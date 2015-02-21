#include <TDirectory.h>
#include <TTreeStream.h>

#include <AliLog.h>

#include "AliMESbaseTask.h"
#include "AliMESeventInfo.h"

ClassImp(AliMESbaseTask)

Int_t AliMESbaseTask::fgDebugUsers(0); 
TTreeSRedirector* AliMESbaseTask::fgDebugStream(NULL); 

//_____________________________________________________________________
AliMESbaseTask::AliMESbaseTask()
  : AliAnalysisTaskSE()
  ,fEvInfo(NULL)
  ,fTracks(NULL)
  ,fMCevInfo(NULL)
  ,fMCtracks(NULL)
{
  // 
  // Constructor
  //
}

//_____________________________________________________________________
AliMESbaseTask::AliMESbaseTask(const char *name)
  : AliAnalysisTaskSE(name)
  ,fEvInfo(NULL)
  ,fTracks(NULL)
  ,fMCevInfo(NULL)
  ,fMCtracks(NULL)
{
  // 
  // Constructor
  //
  DefineOutput(kQA,       TList::Class());
  DefineInput(kEventInfo, AliMESeventInfo::Class());
  DefineInput(kTracks,    TObjArray::Class());
}

//________________________________________________________________________
void AliMESbaseTask::SetMCdata(Bool_t mc)
{ 
  // prepare task for MC processing
  SetBit(kMCdata, mc);
  if(mc){
    DefineInput(kMCeventInfo, AliMESeventInfo::Class());
    DefineInput(kMCtracks,    TObjArray::Class());
  }
}

//_____________________________________________________________________
AliMESbaseTask::~AliMESbaseTask()
{
  // 
  // Destructor
  //
  if(DebugLevel()) CloseDebugStream(); 
}

//________________________________________________________________________
void AliMESbaseTask::UserCreateOutputObjects()
{
  // Build user objects
  BuildQAHistos(); PostData(kQA, fHistosQA);
}

//________________________________________________________________________
void AliMESbaseTask::UserExec(Option_t *)
{
  // link data containers for user processing
/*
  AliInfo("");
  for(Int_t in(0); in<GetNinputs(); in++){
    TObject *o(GetInputData(in));
    if(!o) continue;
    printf("-> slot[%d]=%s\n", in, o->IsA()->GetName());
  }

  for(Int_t is(0); is<GetNoutputs(); is++){
    TObject *o(GetOutputData(is));
    if(!o) continue;
    printf("<- slot[%d]=%s\n", is, o->IsA()->GetName());
  }
  InputEvent();
*/

  fEvInfo   = dynamic_cast<AliMESeventInfo*>(GetInputData(kEventInfo));
  fTracks   = dynamic_cast<TObjArray*>(GetInputData(kTracks));
  if(HasMCdata()){
    fMCevInfo   = dynamic_cast<AliMESeventInfo*>(GetInputData(kMCeventInfo));
    fMCtracks   = dynamic_cast<TObjArray*>(GetInputData(kMCtracks));
  }
  AliDebug(2, Form("Tracks REC[%d] MC[%d]", fTracks->GetEntries(), fMCtracks?fMCtracks->GetEntries():0));


//   PostData(AliMESbaseTask::kQA, fHistosQA);
//   PostData(AliMESbaseTask::kEventInfo+1, fEvInfo);
//   PostData(AliMESbaseTask::kTracks+1, fTracks);
//   PostData(AliMESbaseTask::kMCeventInfo+1, fMCevInfo);
//   PostData(AliMESbaseTask::kMCtracks+1, fMCtracks);

}

//_____________________________________________________________________
void AliMESbaseTask::SetDebugLevel(Int_t level)
{
  //
  // Init debug stream and register user task
  //
  AliAnalysisTaskSE::SetDebugLevel(level);
  if(level>=1 && !fgDebugStream) OpenDebugStream();
  if(level>=1) AddDebugUser(GetName());
}

//_____________________________________________________________________
void AliMESbaseTask::OpenDebugStream()
{
  // open debug stream
  Printf("AliMESbaseTask::OpenDebugStream : Open debug stream to \"MES.Debug.root\"");
  TDirectory *savedir = gDirectory;
  fgDebugStream = new TTreeSRedirector("MES.Debug.root", "RECREATE");
  savedir->cd();
} 

//_____________________________________________________________________
void AliMESbaseTask::CloseDebugStream()
{
  // un-register user and eventually close debug stream

  fgDebugUsers--;
  Printf("AliMESbaseTask::CloseDebugStream : user(s)[%d]", fgDebugUsers);

  if(fgDebugUsers<0) Printf(" - E - AliMESbaseTask::CloseDebugStream : Wrong no. of users[%d] for the debug stream", fgDebugUsers);
  if(fgDebugStream && fgDebugUsers<=0){ 
    delete fgDebugStream;
    fgDebugStream = NULL;
  }
}

//_____________________________________________________________________
void AliMESbaseTask::AddDebugUser(const char *name)
{
  // open debug stream
  fgDebugUsers++;
  Printf("AliMESbaseTask::AddDebugUser : Add task \"%s\" user(s)[%d]", name, fgDebugUsers);
} 

//_____________________________________________________________________
Bool_t AliMESbaseTask::BuildQAHistos()
{
  AliWarning("Function to be implemented in the derived classes");
  return kFALSE;
}
