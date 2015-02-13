#include "AliLog.h"
#include "AliMESchgTask.h"
#include "AliMESeventInfo.h"
#include "AliMEStrackInfo.h"

ClassImp(AliMESchgTask)

//________________________________________________________________________
AliMESchgTask::AliMESchgTask()
  : AliMESbaseTask()
{
  //
  // Constructor
  //
}

//________________________________________________________________________
AliMESchgTask::AliMESchgTask(const char *name)
  : AliMESbaseTask(name)
{
  //
  // Constructor
  //
}

//________________________________________________________________________
AliMESchgTask::~AliMESchgTask()
{
  //
  // Destructor
  //
}

//________________________________________________________________________
void AliMESchgTask::UserCreateOutputObjects()
{
  //define user data containers
  AliMESbaseTask::UserCreateOutputObjects();  

  //define extra user containers
}

//________________________________________________________________________
void AliMESchgTask::UserExec(Option_t *opt)
{
  // Run user analysis. The following objects are allocated after calling AliMESbaseTask::UserExec(opt)
  // fEvInfo  -  reconstructed event information (class AliMESeventInfo)
  // fTracks  -  reconstructed array of tracks (class TObjArray of AliMEStrackInfo)
  // fMCevInfo-  MC event information (class AliMESeventInfo)
  // fMCtracks-  MC array of tracks (class TObjArray of AliMEStrackInfo)
  AliMESbaseTask::UserExec(opt);
}

//________________________________________________________________________
Bool_t AliMESchgTask::PostProcess()
{
  return kTRUE;
}

//________________________________________________________
Bool_t AliMESchgTask::BuildQAHistos()
{
  // Make QA sparse histos for
  fHistosQA = new TList(); fHistosQA->SetOwner(kTRUE);

  return kTRUE;
}
