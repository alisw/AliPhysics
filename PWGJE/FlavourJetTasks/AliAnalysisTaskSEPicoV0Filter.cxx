#include <TString.h>
#include <TClonesArray.h>

#include "AliAnalysisManager.h"

#include "AliAODHandler.h"
#include "AliAODEvent.h"

#include "AliPicoHeaderCJ.h"
#include "AliPicoV0RD.h"
#include "AliPicoV0MC.h"

#include "AliAnalysisTaskSEPicoV0Filter.h"

ClassImp(AliAnalysisTaskSEPicoV0Filter)

//_____________________________________________________________________________
AliAnalysisTaskSEPicoV0Filter::AliAnalysisTaskSEPicoV0Filter() :
AliAnalysisTaskSE(),
fIsAnaInfoMC(kFALSE),
fV0s(0),
fPicoHeaderCJ(0),
fPicoV0sClArr(0),
fListUserOutputs(0)
{
//
//  AliAnalysisTaskSEPicoV0Filter::AliAnalysisTaskSEPicoV0Filter
//
}

//_____________________________________________________________________________
AliAnalysisTaskSEPicoV0Filter::AliAnalysisTaskSEPicoV0Filter(const char *name) :
AliAnalysisTaskSE(name),
fIsAnaInfoMC(kFALSE),
fV0s(0),
fPicoHeaderCJ(0),
fPicoV0sClArr(0),
fListUserOutputs(0)
{
//
//  AliAnalysisTaskSEPicoV0Filter::AliAnalysisTaskSEPicoV0Filter
//

//DefineOutput(2, TList::Class());
}

//_____________________________________________________________________________
AliAnalysisTaskSEPicoV0Filter::~AliAnalysisTaskSEPicoV0Filter()
{
//
//  AliAnalysisTaskSEPicoV0Filter::~AliAnalysisTaskSEPicoV0Filter
//

  if (fV0s)             { delete fV0s;             fV0s             = 0; }
  if (fPicoHeaderCJ)    { delete fPicoHeaderCJ;    fPicoHeaderCJ    = 0; }
  if (fPicoV0sClArr)    { delete fPicoV0sClArr;    fPicoV0sClArr    = 0; }
  if (fListUserOutputs) { delete fListUserOutputs; fListUserOutputs = 0; }
}

//_____________________________________________________________________________
void AliAnalysisTaskSEPicoV0Filter::Init()
{
//
//  AliAnalysisTaskSEPicoV0Filter::Init
//

  return;
}

//_____________________________________________________________________________
void AliAnalysisTaskSEPicoV0Filter::UserCreateOutputObjects()
{
//
//  AliAnalysisTaskSEPicoV0Filter::UserCreateOutputObjects
//

  fPicoHeaderCJ = new  AliPicoHeaderCJ();
  fPicoHeaderCJ->SetName("PicoHeaderCJ");
  AddAODBranch("AliPicoHeaderCJ", &fPicoHeaderCJ);

  if (fIsAnaInfoMC) {
    fPicoV0sClArr = new TClonesArray("AliPicoV0MC");
    fPicoV0sClArr->SetName("PicoV0sMC");
  } else {
    fPicoV0sClArr = new TClonesArray("AliPicoV0RD");
    fPicoV0sClArr->SetName("PicoV0sRD");
  }

  AddAODBranch("TClonesArray",    &fPicoV0sClArr);

/*fListUserOutputs = new TList();
  fListUserOutputs->SetOwner();
  CreateUserOutputHistograms();
  PostData(1, fListUserOutputs);*/
  return;
}

//_____________________________________________________________________________
void AliAnalysisTaskSEPicoV0Filter::Terminate(Option_t */*opt*/)
{
//
//  AliAnalysisTaskSEPicoV0Filter::Terminate
//

  return;
}

//_____________________________________________________________________________
void AliAnalysisTaskSEPicoV0Filter::UserExec(Option_t */*opt*/)
{
//
//  AliAnalysisTaskSEPicoV0Filter::Run
//

  Int_t ncs = 0;
  fPicoHeaderCJ->Reset();
  fPicoV0sClArr->Delete();
//=============================================================================

  AliAnalysisManager::GetAnalysisManager()->GetOutputEventHandler()->SetFillAOD(kTRUE);
//=============================================================================

  fPicoHeaderCJ->SetEventInfo(fInputHandler);
  fV0s = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject("PicoV0s")); if (!fV0s) return;
//=============================================================================

  AliPicoV0RD *pV0RD = 0;
  AliPicoV0MC *pV0MC = 0;
  ncs = fPicoV0sClArr->GetEntriesFast();
  for (Int_t i=0; i<fV0s->GetEntriesFast(); i++) {
    if (fIsAnaInfoMC) {
      pV0MC = static_cast<AliPicoV0MC*>(fV0s->At(i)); if (!pV0MC) continue;
      new ((*fPicoV0sClArr)[ncs++]) AliPicoV0MC(*pV0MC);
      pV0MC = 0;
    } else {
      pV0RD = static_cast<AliPicoV0RD*>(fV0s->At(i)); if (!pV0RD) continue;
      new ((*fPicoV0sClArr)[ncs++]) AliPicoV0RD(*pV0RD);
      pV0RD = 0;
    }
  }

  return;
}

//_____________________________________________________________________________
void AliAnalysisTaskSEPicoV0Filter::CreateUserOutputHistograms()
{
//
//  AliAnalysisTaskSEPicoV0Filter::CreateUserOutputHistograms
//

  if (!fListUserOutputs) return;

  Bool_t bStatusTmpH = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);

  TH1::AddDirectory(bStatusTmpH);
  return;
}
