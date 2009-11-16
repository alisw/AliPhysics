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

#include <TH1F.h>
#include <TList.h>
#include <TIterator.h>

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
  fFMDLowCut(0.7),
  fFMDHitCut(1.2),
  fHistSPD(0),
  fHistV0A(0),       
  fHistV0C(0),    
  fHistZDC(0),    
  fHistFMDA(0),    
  fHistFMDC(0),   
  fHistFMDSingle(0),
  fHistFMDSum(0)
{
}

void AliOfflineTrigger::EnableHistograms()
{
  // creates the monitoring histograms
  
  fHistSPD = new TH1F("fHistSPD", "SPD GFO;number of fired chips;events", 1202, -1.5, 1200.5);
  fHistV0A = new TH1F("fHistV0A", "V0A;number of BB triggers;events", 34, -1.5, 32.5);
  fHistV0C = new TH1F("fHistV0C", "V0C;number of BB triggers;events", 34, -1.5, 32.5);
  fHistZDC = new TH1F("fHistZDC", "ZDC;trigger bits;events", 8, -1.5, 6.5);
  
  // TODO check limits
  fHistFMDA = new TH1F("fHistFMDA", "FMDA;combinations above threshold;events", 102, -1.5, 100.5);
  fHistFMDC = new TH1F("fHistFMDC", "FMDC;combinations above threshold;events", 102, -1.5, 100.5);
  fHistFMDSingle = new TH1F("fHistFMDSingle", "FMD single;multiplicity value;counts", 1000, 0, 10);
  fHistFMDSum = new TH1F("fHistFMDSum", "FMD sum;multiplicity value;counts", 1000, 0, 10);
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

void AliOfflineTrigger::FillHistograms(const AliESDEvent* aEsd)
{
  // fills the histograms with the info from the ESD
  
  fHistSPD->Fill(SPDFiredChips(aEsd));
  
  fHistV0A->Fill(V0BBTriggers(aEsd, kASide));  
  fHistV0C->Fill(V0BBTriggers(aEsd, kCSide));
  
  AliESDZDC* zdcData = aEsd->GetESDZDC();
  if (zdcData)
  {
    UInt_t quality = zdcData->GetESDQuality();
    
    // from Nora's presentation, general first physics meeting 16.10.09
    static UInt_t zpc  = 0x20;
    static UInt_t znc  = 0x10;
    static UInt_t zem1 = 0x08;
    static UInt_t zem2 = 0x04;
    static UInt_t zpa  = 0x02;
    static UInt_t zna  = 0x01;
   
    fHistZDC->Fill(1, quality & zna);
    fHistZDC->Fill(2, quality & zpa);
    fHistZDC->Fill(3, quality & zem2);
    fHistZDC->Fill(4, quality & zem1);
    fHistZDC->Fill(5, quality & znc);
    fHistZDC->Fill(6, quality & zpc);
  }
  else
  {
    fHistZDC->Fill(-1);
    AliError("AliESDZDC not available");
  }
  
  fHistFMDA->Fill(FMDHitCombinations(aEsd, kASide, kTRUE));
  fHistFMDC->Fill(FMDHitCombinations(aEsd, kCSide, kTRUE));
}

Int_t AliOfflineTrigger::SPDFiredChips(const AliESDEvent* aEsd) const
{
  // returns the number of fired chips in the SPD
  
  const AliMultiplicity* mult = aEsd->GetMultiplicity();
  if (!mult)
  {
    AliError("AliMultiplicity not available");
    return -1;
  }
  return mult->GetNumberOfFiredChips(0) + mult->GetNumberOfFiredChips(1);
}

Bool_t AliOfflineTrigger::SPDGFOTrigger(const AliESDEvent* aEsd) const
{
  // Returns if the SPD gave a global Fast OR trigger
  
  Int_t firedChips = SPDFiredChips(aEsd);
  
  if (firedChips >= fSPDGFOThreshold)
    return kTRUE;
  return kFALSE;
}

Int_t AliOfflineTrigger::V0BBTriggers(const AliESDEvent* aEsd, AliceSide side) const
{
  // returns the number of BB triggers in V0A | V0C
  
  AliESDVZERO* v0Data = aEsd->GetVZEROData();
  if (!v0Data)
  {
    AliError("AliESDVZERO not available");
    return -1;
  }
  
  Int_t count = 0;
  for (Int_t i=0; i<32; i++)
  {
    if (side == kASide && v0Data->BBTriggerV0A(i))
      count++;
    if (side == kCSide && v0Data->BBTriggerV0C(i))
      count++;
  }
  
  return count;
}

Bool_t AliOfflineTrigger::V0Trigger(const AliESDEvent* aEsd, AliceSide side) const
{
  // Returns if the V0 triggered
  
  Int_t count = V0BBTriggers(aEsd, side);
  
  if (side == kASide && count >= fV0AThreshold)
    return kTRUE;
  if (side == kCSide && count >= fV0CThreshold)
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

Int_t AliOfflineTrigger::FMDHitCombinations(const AliESDEvent* aEsd, AliceSide side, Bool_t fillHistograms) const
{
  // returns number of hit combinations agove threshold
  //
  // Authors: FMD team, Hans Dalsgaard (code merged from FMD/AliFMDOfflineTrigger)

  // Workaround for AliESDEvent::GetFMDData is not const!
  const AliESDFMD* fmdData = (const_cast<AliESDEvent*>(aEsd))->GetFMDData();
  if (!fmdData)
  {
    AliError("AliESDFMD not available");
    return -1;
  }

  Int_t detFrom = (side == kASide) ? 1 : 3;
  Int_t detTo   = (side == kASide) ? 2 : 3;

  Int_t triggers = 0;
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
          
          if (fillHistograms)
            fHistFMDSingle->Fill(mult);
          
          if (mult > fFMDLowCut)
            totalMult = totalMult + mult;
          else
          {
            if (totalMult > fFMDHitCut)
              triggers++;
              
            if (fillHistograms)
              fHistFMDSum->Fill(totalMult);
              
            totalMult = 0;
          }
        }
      }
    }
  }
  
  return triggers;
}

Bool_t AliOfflineTrigger::FMDTrigger(const AliESDEvent* aEsd, AliceSide side) const
{
  // Returns if the FMD triggered
  //
  // Authors: FMD team, Hans Dalsgaard (code merged from FMD/AliFMDOfflineTrigger)

  Int_t triggers = FMDHitCombinations(aEsd, side, kFALSE);
  
  if (triggers > 0)
    return kTRUE;
    
  return kFALSE;
}

Long64_t AliOfflineTrigger::Merge(TCollection* list)
{
  // Merge a list of AliMultiplicityCorrection objects with this (needed for
  // PROOF).
  // Returns the number of merged objects (including this).

  if (!list)
    return 0;

  if (list->IsEmpty())
    return 1;

  TIterator* iter = list->MakeIterator();
  TObject* obj;

  // collections of all histograms
  const Int_t nHists = 8;
  TList collections[nHists];

  Int_t count = 0;
  while ((obj = iter->Next())) {

    AliOfflineTrigger* entry = dynamic_cast<AliOfflineTrigger*> (obj);
    if (entry == 0) 
      continue;

    collections[0].Add(entry->fHistSPD);
    collections[1].Add(entry->fHistV0A);
    collections[2].Add(entry->fHistV0C);
    collections[3].Add(entry->fHistZDC);
    collections[4].Add(entry->fHistFMDA);
    collections[5].Add(entry->fHistFMDC);
    collections[6].Add(entry->fHistFMDSingle);
    collections[7].Add(entry->fHistFMDSum);

    count++;
  }

  fHistSPD->Merge(&collections[0]);
  fHistV0A->Merge(&collections[1]);
  fHistV0C->Merge(&collections[2]);
  fHistZDC->Merge(&collections[3]);
  fHistFMDA->Merge(&collections[4]);
  fHistFMDC->Merge(&collections[5]);
  fHistFMDSingle->Merge(&collections[6]);
  fHistFMDSum->Merge(&collections[7]);

  delete iter;

  return count+1;
}

void AliOfflineTrigger::WriteHistograms() const
{
  // write histograms to current directory
  
  if (!fHistSPD)
    return;
    
  fHistSPD->Write();
  fHistV0A->Write();
  fHistV0C->Write();
  fHistZDC->Write();
  fHistFMDA->Write();
  fHistFMDC->Write();
  fHistFMDSingle->Write();
  fHistFMDSum->Write();
}
