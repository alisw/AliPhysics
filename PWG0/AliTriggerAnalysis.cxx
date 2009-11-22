/* $Id: AliTriggerAnalysis.cxx 35782 2009-10-22 11:54:31Z jgrosseo $ */

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
//                      Implementation of   Class AliTriggerAnalysis
//   This class provides function to check if events have been triggered based on the data in the ESD
//   The trigger bits, trigger class inputs and only the data (offline trigger) can be used
//   Origin: Jan Fiete Grosse-Oetringhaus, CERN
//-------------------------------------------------------------------------

#include <Riostream.h>
#include <TH1F.h>
#include <TList.h>
#include <TIterator.h>

#include <AliTriggerAnalysis.h>

#include <AliLog.h>

#include <AliESDEvent.h>

#include <AliMultiplicity.h>
#include <AliESDVZERO.h>
#include <AliESDZDC.h>
#include <AliESDFMD.h>

ClassImp(AliTriggerAnalysis)

AliTriggerAnalysis::AliTriggerAnalysis() :
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

void AliTriggerAnalysis::EnableHistograms()
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

//____________________________________________________________________
const char* AliTriggerAnalysis::GetTriggerName(Trigger trigger) 
{
  // returns the name of the requested trigger
  // the returned string will only be valid until the next call to this function [not thread-safe]
  
  static TString str;
  
  UInt_t triggerNoFlags = (UInt_t) trigger % (UInt_t) kStartOfFlags;
  
  switch (triggerNoFlags)
  {
    case kAcceptAll : str = "ACCEPT ALL (bypass!)"; break;
    case kMB1 : str = "MB1"; break;
    case kMB2 : str = "MB2"; break;
    case kMB3 : str = "MB3"; break;
    case kSPDGFO : str = "SPD GFO"; break;
    case kV0A : str = "V0 A"; break;
    case kV0C : str = "V0 C"; break;
    case kZDC : str = "ZDC"; break;
    case kZDCA : str = "ZDC A"; break;
    case kZDCC : str = "ZDC C"; break;
    case kFMDA : str = "FMD A"; break;
    case kFMDC : str = "FMD C"; break;
    case kFPANY : str = "SPD GFO | V0 | ZDC | FMD"; break;
    default: str = ""; break;
  }
   
  if (trigger & kOfflineFlag)
    str += " OFFLINE";  
  
  return str;
}

Bool_t AliTriggerAnalysis::IsTriggerFired(const AliESDEvent* aEsd, Trigger trigger) const
{
  // checks if an event has been triggered

  if (trigger & kOfflineFlag)
    return IsOfflineTriggerFired(aEsd, trigger);
    
  return IsTriggerBitFired(aEsd, trigger);
}

Bool_t AliTriggerAnalysis::IsTriggerBitFired(const AliESDEvent* aEsd, Trigger trigger) const
{
  // checks if an event is fired using the trigger bits

  return IsTriggerBitFired(aEsd->GetTriggerMask(), trigger);
}

Bool_t AliTriggerAnalysis::IsTriggerBitFired(ULong64_t triggerMask, Trigger trigger) const
{
  // checks if an event is fired using the trigger bits
  //
  // this function needs the branch TriggerMask in the ESD
  
  // definitions from p-p.cfg
  ULong64_t spdFO = (1 << 14);
  ULong64_t v0left = (1 << 10);
  ULong64_t v0right = (1 << 11);

  switch (trigger)
  {
    case kAcceptAll:
    {
      return kTRUE;
      break;
    }
    case kMB1:
    {
      if (triggerMask & spdFO || ((triggerMask & v0left) || (triggerMask & v0right)))
        return kTRUE;
      break;
    }
    case kMB2:
    {
      if (triggerMask & spdFO && ((triggerMask & v0left) || (triggerMask & v0right)))
        return kTRUE;
      break;
    }
    case kMB3:
    {
      if (triggerMask & spdFO && (triggerMask & v0left) && (triggerMask & v0right))
        return kTRUE;
      break;
    }
    case kSPDGFO:
    {
      if (triggerMask & spdFO)
        return kTRUE;
      break;
    }
    default:
      Printf("IsEventTriggered: ERROR: Trigger type %d not implemented in this method", (Int_t) trigger);
      break;
  }

  return kFALSE;
}

Bool_t AliTriggerAnalysis::IsTriggerBitFired(const AliESDEvent* aEsd, ULong64_t tclass) const
{
  // Checks if corresponding bit in mask is on
  
  ULong64_t trigmask = aEsd->GetTriggerMask();
  return (trigmask & (1ull << (tclass-1)));
}

Bool_t AliTriggerAnalysis::IsOfflineTriggerFired(const AliESDEvent* aEsd, Trigger trigger) const
{
  // checks if an event has been triggered "offline"

  UInt_t triggerNoFlags = (UInt_t) trigger % (UInt_t) kStartOfFlags;
  
  switch (triggerNoFlags)
  {
    case kAcceptAll:
    {
      return kTRUE;
      break;
    }
    case kMB1:
    {
      if (SPDGFOTrigger(aEsd) || V0Trigger(aEsd, kASide) || V0Trigger(aEsd, kCSide))
        return kTRUE;
      break;
    }
    case kMB2:
    {
      if (SPDGFOTrigger(aEsd) && (V0Trigger(aEsd, kASide) || V0Trigger(aEsd, kCSide)))
        return kTRUE;
      break;
    }
    case kMB3:
    {
      if (SPDGFOTrigger(aEsd) && V0Trigger(aEsd, kASide) && V0Trigger(aEsd, kCSide))
        return kTRUE;
      break;
    }
    case kSPDGFO:
    {
      if (SPDGFOTrigger(aEsd))
        return kTRUE;
      break;
    }
    case kV0A:
    {
      if (V0Trigger(aEsd, kASide))
        return kTRUE;
      break;
    }
    case kV0C:
    {
      if (V0Trigger(aEsd, kCSide))
        return kTRUE;
      break;
    }
    case kZDC:
    {
      if (ZDCTrigger(aEsd, kASide) || ZDCTrigger(aEsd, kCentralBarrel) || ZDCTrigger(aEsd, kCSide))
        return kTRUE;
      break;
    }
    case kZDCA:
    {
      if (ZDCTrigger(aEsd, kASide))
        return kTRUE;
      break;
    }
    case kZDCC:
    {
      if (ZDCTrigger(aEsd, kCSide))
        return kTRUE;
      break;
    }
    case kFMDA:
    {
      if (FMDTrigger(aEsd, kASide))
        return kTRUE;
      break;
    }
    case kFMDC:
    {
      if (FMDTrigger(aEsd, kCSide))
        return kTRUE;
      break;
    }
    case kFPANY:
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


Bool_t AliTriggerAnalysis::IsTriggerClassFired(const AliESDEvent* aEsd, const Char_t* tclass) const 
{
  // tclass is logical function of inputs, e.g. 01 && 02 || 03 && 11 && 21
  // = L0 inp 1 && L0 inp 2 || L0 inp 3 && L1 inp 1 && L2 inp 1
  // NO brackets in logical function !
  // Spaces between operators and inputs.
  // Not all logical functions are available in CTP= 
  // =any function of first 4 inputs; 'AND' of other inputs, check not done
  // This method will be replaced/complemened by similar one
  // which works withh class and inputs names as in CTP cfg file
  
  TString TClass(tclass);
  TObjArray* tcltokens = TClass.Tokenize(" ");
  Char_t level=((TObjString*)tcltokens->At(0))->String()[0];
  UInt_t input=atoi((((TObjString*)tcltokens->At(0))->String()).Remove(0));
  Bool_t tcl = IsInputFired(aEsd,level,input);
 
  for (Int_t i=1;i<tcltokens->GetEntriesFast();i=i+2) {
    level=((TObjString*)tcltokens->At(i+1))->String()[0];
    input=atoi((((TObjString*)tcltokens->At(i+1))->String()).Remove(0));
    Bool_t inpnext = IsInputFired(aEsd,level,input);
    Char_t op =((TObjString*)tcltokens->At(i))->String()[0];
    if (op == '&') tcl=tcl && inpnext;
    else if (op == '|') tcl =tcl || inpnext;
    else {
       AliError(Form("Syntax error in %s", tclass));
       tcltokens->Delete();
       return kFALSE;
    }
  }
  tcltokens->Delete();
  return tcl;
}

Bool_t AliTriggerAnalysis::IsInputFired(const AliESDEvent* aEsd, Char_t level, UInt_t input) const
{
  // Checks trigger input of any level
  
  switch (level)
  {
    case '0': return IsL0InputFired(aEsd,input);
    case '1': return IsL1InputFired(aEsd,input);
    case '2': return IsL2InputFired(aEsd,input);
    default:
      AliError(Form("Wrong level %i",level));
      return kFALSE;
  }
}

Bool_t AliTriggerAnalysis::IsL0InputFired(const AliESDEvent* aEsd, UInt_t input) const 
{
  // Checks if corresponding bit in mask is on
  
  UInt_t inpmask = aEsd->GetHeader()->GetL0TriggerInputs();
  return (inpmask & (1<<(input-1)));
}

Bool_t AliTriggerAnalysis::IsL1InputFired(const AliESDEvent* aEsd, UInt_t input) const
{
  // Checks if corresponding bit in mask is on
  
  UInt_t inpmask = aEsd->GetHeader()->GetL1TriggerInputs();
  return (inpmask & (1<<(input-1)));
}

Bool_t AliTriggerAnalysis::IsL2InputFired(const AliESDEvent* aEsd, UInt_t input) const 
{
  // Checks if corresponding bit in mask is on
  
  UInt_t inpmask = aEsd->GetHeader()->GetL2TriggerInputs();
  return (inpmask & (1<<(input-1)));
}

void AliTriggerAnalysis::FillHistograms(const AliESDEvent* aEsd) 
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

Int_t AliTriggerAnalysis::SPDFiredChips(const AliESDEvent* aEsd) const
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

Bool_t AliTriggerAnalysis::SPDGFOTrigger(const AliESDEvent* aEsd) const
{
  // Returns if the SPD gave a global Fast OR trigger
  
  Int_t firedChips = SPDFiredChips(aEsd);
  
  if (firedChips >= fSPDGFOThreshold)
    return kTRUE;
  return kFALSE;
}

Int_t AliTriggerAnalysis::V0BBTriggers(const AliESDEvent* aEsd, AliceSide side) const
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

Bool_t AliTriggerAnalysis::V0Trigger(const AliESDEvent* aEsd, AliceSide side) const
{
  // Returns if the V0 triggered
  
  Int_t count = V0BBTriggers(aEsd, side);
  
  if (side == kASide && count >= fV0AThreshold)
    return kTRUE;
  if (side == kCSide && count >= fV0CThreshold)
    return kTRUE;
  return kFALSE;
}

Bool_t AliTriggerAnalysis::ZDCTrigger(const AliESDEvent* aEsd, AliceSide side) const
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

Int_t AliTriggerAnalysis::FMDHitCombinations(const AliESDEvent* aEsd, AliceSide side, Bool_t fillHistograms) const
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

Bool_t AliTriggerAnalysis::FMDTrigger(const AliESDEvent* aEsd, AliceSide side) const
{
  // Returns if the FMD triggered
  //
  // Authors: FMD team, Hans Dalsgaard (code merged from FMD/AliFMDOfflineTrigger)

  Int_t triggers = FMDHitCombinations(aEsd, side, kFALSE);
  
  if (triggers > 0)
    return kTRUE;
    
  return kFALSE;
}

Long64_t AliTriggerAnalysis::Merge(TCollection* list)
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

    AliTriggerAnalysis* entry = dynamic_cast<AliTriggerAnalysis*> (obj);
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

void AliTriggerAnalysis::WriteHistograms() const
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
