/* $Id: AliOfflineTrigger.cxx 35782 2009-10-22 11:54:31Z jgrosseo $ */

#include <AliOfflineTrigger.h>

#include <AliLog.h>

#include <AliESDEvent.h>

#include <AliMultiplicity.h>
#include <AliESDVZERO.h>

ClassImp(AliOfflineTrigger)

AliOfflineTrigger::AliOfflineTrigger() :
  fSPDGFOThreshold(1),
  fV0AThreshold(1),
  fV0CThreshold(1)
{
}

Bool_t AliOfflineTrigger::IsEventTriggered(const AliESDEvent* aEsd, AliPWG0Helper::Trigger trigger) const
{
  // checks if an event has been triggered "offline"

  UInt_t triggerNoFlags = (UInt_t) trigger % (UInt_t) AliPWG0Helper::kStartOfFlags;
  
  switch (triggerNoFlags)
  {
    case AliPWG0Helper::kAcceptAll:
    {
      return kTRUE;
      break;
    }
    case AliPWG0Helper::kMB1:
    {
      if (SPDGFOTrigger(aEsd) || V0Trigger(aEsd, kASide) || V0Trigger(aEsd, kCSide))
        return kTRUE;
      break;
    }
    case AliPWG0Helper::kMB2:
    {
      if (SPDGFOTrigger(aEsd) && (V0Trigger(aEsd, kASide) || V0Trigger(aEsd, kCSide)))
        return kTRUE;
      break;
    }
    case AliPWG0Helper::kMB3:
    {
      if (SPDGFOTrigger(aEsd) && V0Trigger(aEsd, kASide) && V0Trigger(aEsd, kCSide))
        return kTRUE;
      break;
    }
    case AliPWG0Helper::kSPDGFO:
    {
      if (SPDGFOTrigger(aEsd))
        return kTRUE;
      break;
    }
    case AliPWG0Helper::kV0A:
    {
      if (V0Trigger(aEsd, kASide))
        return kTRUE;
      break;
    }
    case AliPWG0Helper::kV0C:
    {
      if (V0Trigger(aEsd, kCSide))
        return kTRUE;
      break;
    }
    case AliPWG0Helper::kZDCA:
    {
      if (ZDCTrigger(aEsd, kASide))
        return kTRUE;
      break;
    }
    case AliPWG0Helper::kZDCC:
    {
      if (ZDCTrigger(aEsd, kCSide))
        return kTRUE;
      break;
    }
    case AliPWG0Helper::kFMD:
    {
      if (FMDTrigger(aEsd))
        return kTRUE;
      break;
    }
    default:
    {
      AliFatal(Form("Trigger type %d not implemented", triggerNoFlags));
    }
  }
  
  return kFALSE;
}

Bool_t AliOfflineTrigger::SPDGFOTrigger(const AliESDEvent* aEsd) const
{
  // Returns if the SPD gave a global Fast OR trigger
  
  Int_t firedChips = 0;
  const AliMultiplicity* mult = aEsd->GetMultiplicity();
  if (!mult)
  {
    AliError("AliMultiplicity not available");
    return kFALSE;
  }
  firedChips = mult->GetNumberOfFiredChips(0) + mult->GetNumberOfFiredChips(1);
  
  if (firedChips >= fSPDGFOThreshold)
    return kTRUE;
  return kFALSE;
}

Bool_t AliOfflineTrigger::V0Trigger(const AliESDEvent* aEsd, AliceSide side) const
{
  // Returns if the V0 triggered
  
  AliESDVZERO* v0Data = aEsd->GetVZEROData();
  if (!v0Data)
  {
    AliError("AliESDVZERO not available");
    return kFALSE;
  }
  
  Int_t aCount = 0;
  Int_t cCount = 0;
  for (Int_t i=0; i<32; i++)
  {
    if (v0Data->BBTriggerV0A(i))
      aCount++;
    if (v0Data->BBTriggerV0C(i))
      cCount++;
  }
  
  if (side == kASide && aCount >= fV0AThreshold)
    return kTRUE;
  if (side == kCSide && cCount >= fV0CThreshold)
    return kTRUE;
  return kFALSE;  
}

Bool_t AliOfflineTrigger::ZDCTrigger(const AliESDEvent* /* aEsd */, AliceSide /* side */) const
{
  // Returns if ZDC triggered
  
  AliFatal("Not implemented");
  
  return kFALSE;
}

Bool_t AliOfflineTrigger::FMDTrigger(const AliESDEvent* /* aEsd */) const
{
  // Returns if the FMD triggered

  AliFatal("Not implemented");
  
  return kFALSE;
}
