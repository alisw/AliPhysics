#include <THnSparse.h>

#include <AliLog.h>
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
  
  // ************
  // test
  Double_t vec_hMultEst[2]; // vector used to fill hMultEst
  THnSparseD *hMultEst = (THnSparseD*)fHistosQA->At(0);

  Double_t mult_comb08 = fEvInfo->GetMultiplicity(AliMESeventInfo::kComb);  	// combined multiplicity with |eta| < 0.8
  Double_t mult_V0M = fEvInfo->GetMultiplicity(AliMESeventInfo::kV0M);
  
  vec_hMultEst[0] = mult_comb08;
  vec_hMultEst[1] = mult_V0M;
  
  hMultEst->Fill(vec_hMultEst);
  // ************
  
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
  
  
  // ************
  // test
  const Int_t ndim(2);
  const Int_t cldNbins[ndim]   = {105, 100};
  const Double_t cldMin[ndim]  = {-5, 0.,}, cldMax[ndim]  = {100., 100.,};
  THnSparseD *hMultEst = new THnSparseD("hMultEst","hMultEst;combined 0.8;V0M;",ndim, cldNbins, cldMin, cldMax);
  fHistosQA->AddAt(hMultEst, 0);
  // ************
  

  return kTRUE;
}
