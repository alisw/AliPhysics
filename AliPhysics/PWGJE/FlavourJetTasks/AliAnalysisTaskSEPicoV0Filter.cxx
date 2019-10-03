#include <TClonesArray.h>

#include "AliAnalysisManager.h"

#include "AliAODEvent.h"
#include "AliAODHandler.h"

#include "AliPicoV0RD.h"
#include "AliPicoV0MC.h"
#include "AliPicoHeaderV0.h"

#include "AliAnalysisTaskSEPicoV0Filter.h"

ClassImp(AliAnalysisTaskSEPicoV0Filter)

//_____________________________________________________________________________
AliAnalysisTaskSEPicoV0Filter::AliAnalysisTaskSEPicoV0Filter() :
AliAnalysisTaskSE(),
fIsMC(kFALSE),
fMult(""),
fMultEstDef(""),
fCutMinMult(0.),
fCutMaxMult(0.),
fV0s(nullptr),
fPicoHeader(nullptr),
fPicoV0sClArr(nullptr),
fListUserOutputs(nullptr)
{
//
//  AliAnalysisTaskSEPicoV0Filter::AliAnalysisTaskSEPicoV0Filter
//
}

//_____________________________________________________________________________
AliAnalysisTaskSEPicoV0Filter::AliAnalysisTaskSEPicoV0Filter(const char *name) :
AliAnalysisTaskSE(name),
fIsMC(kFALSE),
fMult(""),
fMultEstDef(""),
fCutMinMult(-99999.),
fCutMaxMult(999999.),
fV0s(nullptr),
fPicoHeader(nullptr),
fPicoV0sClArr(nullptr),
fListUserOutputs(nullptr)
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

  if (fV0s)             { delete fV0s;             fV0s             = nullptr; }
  if (fPicoHeader)      { delete fPicoHeader;      fPicoHeader      = nullptr; }
  if (fPicoV0sClArr)    { delete fPicoV0sClArr;    fPicoV0sClArr    = nullptr; }
  if (fListUserOutputs) { delete fListUserOutputs; fListUserOutputs = nullptr; }
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

  if (fPicoHeader) {
    delete fPicoHeader;
    fPicoHeader = nullptr;
  }

  fPicoHeader = new  AliPicoHeaderV0();
  fPicoHeader->SetName("PicoHeaderV0");

  if (!fMult.IsNull()) {
    const auto aMult(fMult.Tokenize(":"));

    TIter next(aMult);
    TObjString *ps(nullptr);
    while ((ps = static_cast<TObjString*>(next()))) {
      fPicoHeader->AddMultEstimator(ps->String());
    }
  }

  AddAODBranch("AliPicoHeaderV0", &fPicoHeader);
//=============================================================================

  if (fPicoV0sClArr) {
    delete fPicoV0sClArr;
    fPicoV0sClArr = nullptr;
  }

  if (fIsMC) {
    fPicoV0sClArr = new TClonesArray("AliPicoV0MC");
    fPicoV0sClArr->SetName("PicoV0sMC");
  } else {
    fPicoV0sClArr = new TClonesArray("AliPicoV0RD");
    fPicoV0sClArr->SetName("PicoV0sRD");
  }

  AddAODBranch("TClonesArray", &fPicoV0sClArr);
//=============================================================================

/*if (fListUserOutputs) {
    delete fListUserOutputs;
    fListUserOutputs = nullptr;
  }

  fListUserOutputs = new TList();
  fListUserOutputs->SetOwner();
  const auto b(TH1::AddDirectoryStatus());
  TH1::AddDirectory(kFALSE);
//TODO
  TH1::AddDirectory(b);
  PostData(1, fListUserOutputs);*/
//=============================================================================

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

  fPicoHeader->Reset();
  fPicoV0sClArr->Delete();
//=============================================================================

  if (!fV0s) {
    fV0s = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject("PicoV0s"));
    if (!fV0s) return;
  }

  const auto nV0s(fV0s->GetEntriesFast());
  if (nV0s<=0) return;
//=============================================================================

  fPicoHeader->SetEventInfo(fInputHandler);

  if (!fMultEstDef.IsNull()) {
    const auto dMult(fPicoHeader->MultiplicityPercentile(fMultEstDef));
    if ((dMult<fCutMinMult) || (dMult>=fCutMaxMult)) {
      fPicoHeader->Reset();
      return;
    }
  }
//=============================================================================

  auto l(fPicoV0sClArr->GetEntriesFast());

  for (auto i=0; i<nV0s; ++i) {
    if (fIsMC) {
      const auto pV0(static_cast<AliPicoV0MC*>(fV0s->At(i))); if (!pV0) continue;
      new ((*fPicoV0sClArr)[l++]) AliPicoV0MC(*pV0);
    } else {
      const auto pV0(static_cast<AliPicoV0RD*>(fV0s->At(i))); if (!pV0) continue;
      new ((*fPicoV0sClArr)[l++]) AliPicoV0RD(*pV0);
    }
  }
//=============================================================================

  AliAnalysisManager::GetAnalysisManager()->GetOutputEventHandler()->SetFillAOD(kTRUE);
//=============================================================================

  return;
}
