/* $Id: AliOfflineTrigger.cxx 35782 2009-10-22 11:54:31Z jgrosseo $ */

/**************************************************************************
 * Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

//-------------------------------------------------------------------------
//                      Implementation of   Class AliOfflineTrigger
//   This class provides offline triggers from data in the ESD
//   Origin: Jan Fiete Grosse-Oetringhaus, CERN
//-------------------------------------------------------------------------

#include <AliOfflineTrigger.h>

#include <AliLog.h>

#include <AliESDEvent.h>

#include <AliMultiplicity.h>
#include <AliESDVZERO.h>
#include <AliESDZDC.h>
#include <AliESDFMD.h>

ClassImp(AliOfflineTrigger)

AliOfflineTrigger::AliOfflineTrigger() :
  fSPDGFOThreshold(1),
  fV0AThreshold(1),
  fV0CThreshold(1),
  fFMDLowCut(0.2),
  fFMDHitCut(0.5)
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
    case AliPWG0Helper::kZDC:
    {
      if (ZDCTrigger(aEsd, kASide) || ZDCTrigger(aEsd, kCentralBarrel) || ZDCTrigger(aEsd, kCSide))
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
    case AliPWG0Helper::kFMDA:
    {
      if (FMDTrigger(aEsd, kASide))
        return kTRUE;
      break;
    }
    case AliPWG0Helper::kFMDC:
    {
      if (FMDTrigger(aEsd, kCSide))
        return kTRUE;
      break;
    }
    case AliPWG0Helper::kFPANY:
    {
      if (SPDGFOTrigger(aEsd) || V0Trigger(aEsd, kASide) || V0Trigger(aEsd, kCSide) || ZDCTrigger(aEsd, kASide) || ZDCTrigger(aEsd, kCentralBarrel) || ZDCTrigger(aEsd, kCSide) || FMDTrigger(aEsd, kASide) || FMDTrigger(aEsd, kCSide))
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

Bool_t AliOfflineTrigger::ZDCTrigger(const AliESDEvent* aEsd, AliceSide side) const
{
  // Returns if ZDC triggered
  
  AliESDZDC* zdcData = aEsd->GetESDZDC();
  if (!zdcData)
  {
    AliError("AliESDZDC not available");
    return kFALSE;
  }
  
  UInt_t quality = zdcData->GetESDQuality();
  
  // from Nora's presentation, general first physics meeting 16.10.09
  static UInt_t zpc  = 0x20;
  static UInt_t znc  = 0x10;
  static UInt_t zem1 = 0x08;
  static UInt_t zem2 = 0x04;
  static UInt_t zpa  = 0x02;
  static UInt_t zna  = 0x01;
  
  if (side == kASide && ((quality & zpa) || (quality & zna)))
    return kTRUE;
  if (side == kCentralBarrel && ((quality & zem1) || (quality & zem2)))
    return kTRUE;
  if (side == kCSide && ((quality & zpc) || (quality & znc)))
    return kTRUE;
  
  return kFALSE;
}

Bool_t AliOfflineTrigger::FMDTrigger(const AliESDEvent* aEsd, AliceSide side) const
{
  // Returns if the FMD triggered
  //
  // Authors: FMD team, Hans Dalsgaard (code merged from FMD/AliFMDOfflineTrigger)

  // Workaround for AliESDEvent::GetFMDData is not const!
  const AliESDFMD* fmdData = (const_cast<AliESDEvent*>(aEsd))->GetFMDData();
  if (!fmdData)
  {
    AliError("AliESDFMD not available");
    return kFALSE;
  }

  Int_t detFrom = (side == kASide) ? 1 : 3;
  Int_t detTo   = (side == kASide) ? 2 : 3;

  Float_t totalMult = 0;
  for (UShort_t det=detFrom;det<=detTo;det++) {
    Int_t nRings = (det == 1 ? 1 : 2);
    for (UShort_t ir = 0; ir < nRings; ir++) {
      Char_t   ring = (ir == 0 ? 'I' : 'O');
      UShort_t nsec = (ir == 0 ? 20  : 40);
      UShort_t nstr = (ir == 0 ? 512 : 256);
      for (UShort_t sec =0; sec < nsec;  sec++) {
        for (UShort_t strip = 0; strip < nstr; strip++) {
          Float_t mult = fmdData->Multiplicity(det,ring,sec,strip);
          if (mult == AliESDFMD::kInvalidMult) continue;
          
          if (mult > fFMDLowCut)
            totalMult = totalMult + mult;
          else
            {
              if( totalMult > fFMDHitCut) {
          return kTRUE;
              }
              else totalMult = 0 ;
            }
        }
      }
    }
  }
  return kFALSE;
}
