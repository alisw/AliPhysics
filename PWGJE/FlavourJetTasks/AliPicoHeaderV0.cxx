#include "AliVEvent.h"
#include "AliVVertex.h"

#include "AliCentrality.h"
#include "AliEventplane.h"
#include "AliMultSelection.h"
#include "AliInputEventHandler.h"

#include "AliPicoHeaderV0.h"

ClassImp(AliPicoHeaderV0)

//_____________________________________________________________________________
AliPicoHeaderV0::AliPicoHeaderV0(const TString s) :
TNamed(s.Data(), ""),
fPS(0),
fTrg(""),
fEP(-999.),
fAct(nullptr)
{
//
// AliPicoHeaderV0::AliPicoHeaderV0
//

  for (auto &d : fVtx) d = -999.;
}

//_____________________________________________________________________________
AliPicoHeaderV0::AliPicoHeaderV0(const AliPicoHeaderV0 &src) :
TNamed(src),
fPS(src.fPS),
fTrg(src.fTrg),
fEP(src.fEP)
{
//
// AliPicoHeaderV0::AliPicoHeaderV0
//

  if (src.fAct) {
    TObjString *ps(nullptr);
    if (fAct) { delete fAct; fAct = nullptr; }

    fAct = new TMap();
    const auto next(src.fAct->MakeIterator());
    while ((ps = static_cast<TObjString*>((*next)()))) {
      const auto s(ps->String());
      const auto p(static_cast<TParameter<Float_t>*>((*src.fAct)(s.Data())));

      if (p) fAct->Add(new TObjString(s.Data()),
                       new TParameter<Float_t>(s.Data(),p->GetVal()));
    }
  }

  auto l(0);
  for (auto &d : src.fVtx) fVtx[l++] = d;
}

//_____________________________________________________________________________
AliPicoHeaderV0& AliPicoHeaderV0::operator=(const AliPicoHeaderV0 &src)
{
//
// AliPicoHeaderV0::operator=
//

  if (&src==this) return *this;

  TNamed::operator=(src);

  fPS  = src.fPS;
  fTrg = src.fTrg;
  fEP  = src.fEP;

  if (src.fAct) {
    TObjString *ps(nullptr);
    if (fAct) { delete fAct; fAct = nullptr; }

    fAct = new TMap();
    const auto next(src.fAct->MakeIterator());
    while ((ps = static_cast<TObjString*>((*next)()))) {
      const auto s(ps->String());
      const auto p(static_cast<TParameter<Float_t>*>((*src.fAct)(s.Data())));

      if (p) fAct->Add(new TObjString(s.Data()),
                       new TParameter<Float_t>(s.Data(),p->GetVal()));
    }
  }

  auto l(0);
  for (auto &d : src.fVtx) fVtx[l++] = d;

  return *this;
}

//_____________________________________________________________________________
AliPicoHeaderV0::~AliPicoHeaderV0()
{
//
// AliPicoHeaderV0::~AliPicoHeaderV0
//
  if (fAct) { delete fAct; fAct = nullptr; }
}

//_____________________________________________________________________________
void AliPicoHeaderV0::SetEventInfo(AliVEventHandler* const pH)
{
//
// AliPicoHeaderV0::SetEventInfo
//

  const auto pV(pH->GetEvent()); if (!pV) return;
//============================================================================-

  fPS  = pH->IsEventSelected();
  fTrg = pV->GetFiredTriggerClasses();
//============================================================================-

  const auto pVtx(pV->GetPrimaryVertex());
  SetTitle(pVtx->GetTitle());
  pVtx->GetXYZ(fVtx);
//============================================================================-

  if (fAct) {
    TObjString *ps(nullptr);
    const auto next(fAct->MakeIterator());

    if (TString(GetName()).Contains("old")) {
      const auto pm(pV->GetCentrality());
      if (pm) while ((ps = static_cast<TObjString*>((*next)()))) {
        const auto s(ps->String());
        const auto p(static_cast<TParameter<Float_t>*>((*fAct)(s.Data())));
        if (p) p->SetVal(pm->GetCentralityPercentile(s.Data()));
      }
    } else {
      const auto pm(static_cast<AliMultSelection*>(pV->FindListObject("MultSelection")));
      if (pm) while ((ps = static_cast<TObjString*>((*next)()))) {
        const auto s(ps->String());
        const auto p(static_cast<TParameter<Float_t>*>((*fAct)(s.Data())));
        if (p) p->SetVal(pm->GetMultiplicityPercentile(s.Data()));
      }
    }
  }
//============================================================================-

  const auto pEP(pV->GetEventplane());
  if (pEP) fEP = pEP->GetEventplane("Q");
//============================================================================-

  return;
}

//_____________________________________________________________________________
Float_t AliPicoHeaderV0::MultiplicityPercentile(const TString &s)
{
//
//  AliPicoHeaderV0::CentralityPercentile(const TString s) const
//

  if (!fAct) return -999.;
  if (s.IsNull()) return -999.;
//=============================================================================

  const auto p(static_cast<TParameter<Float_t>*>((*fAct)(s.Data())));
  if (!p) return -999.;
//=============================================================================

  return p->GetVal();
}

//_____________________________________________________________________________
void AliPicoHeaderV0::Reset()
{
//
//  AliPicoHeaderV0::Reset
//

  fPS  = 0;
  fTrg = "";

  fEP = -999.;
  for (auto &d : fVtx) d = -999.;
//=============================================================================

  if (fAct) {
    TObjString *ps(nullptr);
    const auto next(fAct->MakeIterator());
    while ((ps = static_cast<TObjString*>((*next)()))) {
      const auto s(ps->String());
      const auto p(static_cast<TParameter<Float_t>*>((*fAct)(s.Data())));
      if (p) p->SetVal(-999.);
    }
  }
//=============================================================================

  return;
}
