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
#include <TH2F.h>
#include <TList.h>
#include <TIterator.h>
#include "TParameter.h"
#include <TMap.h>

#include <AliTriggerAnalysis.h>

#include <AliLog.h>

#include <AliESDEvent.h>

#include <AliMultiplicity.h>
#include <AliESDVZERO.h>
#include <AliESDZDC.h>
#include <AliESDFMD.h>

ClassImp(AliTriggerAnalysis)

AliTriggerAnalysis::AliTriggerAnalysis() :
  fSPDGFOThreshold(2),
  fV0TimeOffset(0),
  fFMDLowCut(0.2),
  fFMDHitCut(0.5),
  fHistBitsSPD(0),
  fHistFiredBitsSPD(0),
  fHistV0A(0),       
  fHistV0C(0),
  fHistZDC(0),    
  fHistFMDA(0),    
  fHistFMDC(0),   
  fHistFMDSingle(0),
  fHistFMDSum(0),
  fTriggerClasses(0)
{
  // constructor
}

AliTriggerAnalysis::~AliTriggerAnalysis()
{
  // destructor
  
  if (fHistBitsSPD)
  {
    delete fHistBitsSPD;
    fHistBitsSPD = 0;
  }

  if (fHistFiredBitsSPD)
  {
    delete fHistFiredBitsSPD;
    fHistFiredBitsSPD = 0;
  }

  if (fHistV0A)
  {
    delete fHistV0A;
    fHistV0A = 0;
  }

  if (fHistV0C)
  {
    delete fHistV0C;
    fHistV0C = 0;
  }

  if (fHistZDC)
  {
    delete fHistZDC;
    fHistZDC = 0;
  }

  if (fHistFMDA)
  {
    delete fHistFMDA;
    fHistFMDA = 0;
  }

  if (fHistFMDC)
  {
    delete fHistFMDC;
    fHistFMDC = 0;
  }

  if (fHistFMDSingle)
  {
    delete fHistFMDSingle;
    fHistFMDSingle = 0;
  }

  if (fHistFMDSum)
  {
    delete fHistFMDSum;
    fHistFMDSum = 0;
  }

  if (fTriggerClasses)
  {
    delete fTriggerClasses;
    fTriggerClasses = 0;
  }
}

void AliTriggerAnalysis::EnableHistograms()
{
  // creates the monitoring histograms
  
  fHistBitsSPD = new TH2F("fHistBitsSPD", "SPD GFO;number of fired chips (offline);number of fired chips (hardware)", 1202, -1.5, 1200.5, 1202, -1.5, 1200.5);
  fHistFiredBitsSPD = new TH1F("fHistFiredBitsSPD", "SPD GFO Hardware;chip number;events", 1200, -0.5, 1199.5);
  fHistV0A = new TH1F("fHistV0A", "V0A;leading time (ns);events", 200, 0, 100);
  fHistV0C = new TH1F("fHistV0C", "V0C;leading time (ns);events", 200, 0, 100);
  fHistZDC = new TH1F("fHistZDC", "ZDC;trigger bits;events", 8, -1.5, 6.5);
  
  // TODO check limits
  fHistFMDA = new TH1F("fHistFMDA", "FMDA;combinations above threshold;events", 102, -1.5, 100.5);
  fHistFMDC = new TH1F("fHistFMDC", "FMDC;combinations above threshold;events", 102, -1.5, 100.5);
  fHistFMDSingle = new TH1F("fHistFMDSingle", "FMD single;multiplicity value;counts", 1000, 0, 10);
  fHistFMDSum = new TH1F("fHistFMDSum", "FMD sum;multiplicity value;counts", 1000, 0, 10);
  
  fTriggerClasses = new TMap;
  fTriggerClasses->SetOwner();
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
    case kSPDGFOBits : str = "SPD GFO Bits"; break;
    case kV0A : str = "V0 A BB"; break;
    case kV0C : str = "V0 C BB"; break;
    case kV0ABG : str = "V0 A BG"; break;
    case kV0CBG : str = "V0 C BG"; break;
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

Bool_t AliTriggerAnalysis::IsTriggerFired(const AliESDEvent* aEsd, Trigger trigger)
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

Bool_t AliTriggerAnalysis::IsOfflineTriggerFired(const AliESDEvent* aEsd, Trigger trigger)
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
      if (SPDGFOTrigger(aEsd, 0) || V0Trigger(aEsd, kASide) == kV0BB || V0Trigger(aEsd, kCSide) == kV0BB)
        return kTRUE;
      break;
    }
    case kMB2:
    {
      if (SPDGFOTrigger(aEsd, 0) && (V0Trigger(aEsd, kASide) == kV0BB || V0Trigger(aEsd, kCSide) == kV0BB))
        return kTRUE;
      break;
    }
    case kMB3:
    {
      if (SPDGFOTrigger(aEsd, 0) && V0Trigger(aEsd, kASide) == kV0BB && V0Trigger(aEsd, kCSide) == kV0BB)
        return kTRUE;
      break;
    }
    case kSPDGFO:
    {
      if (SPDGFOTrigger(aEsd, 0))
        return kTRUE;
      break;
    }
    case kSPDGFOBits:
    {
      if (SPDGFOTrigger(aEsd, 1))
        return kTRUE;
      break;
    }
    case kV0A:
    {
      if (V0Trigger(aEsd, kASide) == kV0BB)
        return kTRUE;
      break;
    }
    case kV0C:
    {
      if (V0Trigger(aEsd, kCSide) == kV0BB)
        return kTRUE;
      break;
    }
    case kV0ABG:
    {
      if (V0Trigger(aEsd, kASide) == kV0BG)
        return kTRUE;
      break;
    }
    case kV0CBG:
    {
      if (V0Trigger(aEsd, kCSide) == kV0BG)
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
      if (SPDGFOTrigger(aEsd, 0) || V0Trigger(aEsd, kASide) == kV0BB || V0Trigger(aEsd, kCSide) == kV0BB || ZDCTrigger(aEsd, kASide) || ZDCTrigger(aEsd, kCentralBarrel) || ZDCTrigger(aEsd, kCSide) || FMDTrigger(aEsd, kASide) || FMDTrigger(aEsd, kCSide))
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
  
  fHistBitsSPD->Fill(SPDFiredChips(aEsd, 0), SPDFiredChips(aEsd, 1, kTRUE));
  
  V0Trigger(aEsd, kASide, kTRUE);
  V0Trigger(aEsd, kCSide, kTRUE);
  
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
  
void AliTriggerAnalysis::FillTriggerClasses(const AliESDEvent* aEsd)
{
  // fills trigger classes map
  
  TParameter<Long64_t>* count = dynamic_cast<TParameter<Long64_t>*> (fTriggerClasses->GetValue(aEsd->GetFiredTriggerClasses().Data()));
  if (!count)
  {
    count = new TParameter<Long64_t>(aEsd->GetFiredTriggerClasses(), 0);
    fTriggerClasses->Add(new TObjString(aEsd->GetFiredTriggerClasses().Data()), count);
  }
  count->SetVal(count->GetVal() + 1);
  
  // TODO add first and last orbit number here
}

Int_t AliTriggerAnalysis::SPDFiredChips(const AliESDEvent* aEsd, Int_t origin, Bool_t fillHists)
{
  // returns the number of fired chips in the SPD
  //
  // origin = 0 --> aEsd->GetMultiplicity()->GetNumberOfFiredChips() (filled from clusters)
  // origin = 1 --> aEsd->GetMultiplicity()->TestFastOrFiredChips() (from hardware bits)
  
  const AliMultiplicity* mult = aEsd->GetMultiplicity();
  if (!mult)
  {
    AliError("AliMultiplicity not available");
    return -1;
  }
  
  if (origin == 0)
    return mult->GetNumberOfFiredChips(0) + mult->GetNumberOfFiredChips(1);
    
  if (origin == 1)
  {
    Int_t nChips = 0;
    for (Int_t i=0; i<1200; i++)
      if (mult->TestFastOrFiredChips(i) == kTRUE)
      {
        nChips++;
        if (fillHists)
          fHistFiredBitsSPD->Fill(i);
      }
    return nChips;
  }
  
  return -1;
}

Bool_t AliTriggerAnalysis::SPDGFOTrigger(const AliESDEvent* aEsd, Int_t origin)
{
  // Returns if the SPD gave a global Fast OR trigger
  
  Int_t firedChips = SPDFiredChips(aEsd, origin);
  
  if (firedChips >= fSPDGFOThreshold)
    return kTRUE;
  return kFALSE;
}

AliTriggerAnalysis::V0Decision AliTriggerAnalysis::V0Trigger(const AliESDEvent* aEsd, AliceSide side, Bool_t fillHists)
{
  // Returns the V0 trigger decision in V0A | V0C
  //
  // Based on algorithm by Cvetan Cheshkov
  
  AliESDVZERO* esdV0 = aEsd->GetVZEROData();
  if (!esdV0)
  {
    AliError("AliESDVZERO not available");
    return kV0Invalid;
  }
  
  Int_t begin = -1;
  Int_t end = -1;
  
  if (side == kASide)
  {
    begin = 32;
    end = 64;
  } 
  else if (side == kCSide)
  {
    begin = 0;
    end = 32;
  }
  else
    return kV0Invalid;
  
  Float_t time = 0;
  Int_t ntime = 0;
  for (Int_t i = begin; i < end; ++i) {
    if (esdV0->GetTime(i) > 1e-6 && esdV0->GetAdc(i) > 6.0) {
      Float_t correctedTime = V0CorrectLeadingTime(i, esdV0->GetTime(i), esdV0->GetAdc(i));
    
      time += correctedTime;
      time += fV0TimeOffset;
      
      ntime++;
    }
  }

  if (ntime > 0) 
    time /= ntime;

  if (fillHists)
  {
    if (side == kASide && fHistV0A)
      fHistV0A->Fill(time);
    if (side == kCSide && fHistV0C)
      fHistV0C->Fill(time);
  }
  
  if (side == kASide)
  {
    if (time > 48 && time < 62)
      return kV0BB;
    if (time > 26 && time < 33) 
      return kV0BG;
  }
  
  if (side == kCSide)
  {
    if (time > 49 && time < 60)
      return kV0BB;
    if (time > 43 && time < 48.5)
      return kV0BG; 
  }
  
  return kV0Empty;
}

Float_t AliTriggerAnalysis::V0CorrectLeadingTime(Int_t i, Float_t time, Float_t adc) const
{
  // Correct for slewing and align the channels
  //
  // Author: Cvetan Cheshkov

  if (time == 0) return 0;

  // Time offsets between channels
  Float_t timeShift[64] = {30.2914 , 30.0019 , 30.7429 , 30.1997 , 30.1511 , 29.6437 , 30.0609 , 29.5452 , 30.1437 , 30.745 , 30.7537 , 30.446 , 30.2771 , 30.838 , 30.3748 , 30.0635 , 30.1786 , 30.282 , 31.0992 , 30.7491 , 30.624 , 30.9268 , 30.6585 , 30.4895 , 31.5815 , 31.3871 , 31.2032 , 31.5778 , 31.0838 , 31.2259 , 31.2122 , 31.5989 , 28.3792 , 28.8325 , 27.8719 , 28.3475 , 26.9925 , 27.9300 , 28.4223 , 28.4996 , 28.2934 , 28.1281 , 27.209 , 28.5327 , 28.1181 , 28.0888 , 29.5111 , 28.6601 , 29.7705 , 29.6531 , 30.3373 , 30.2345 , 30.5935 , 29.8164 , 30.2235 , 29.6505 , 30.1225 , 31.2045 , 30.8399 , 30.6789 , 30.2784 , 31.7028 , 31.4239 , 30.1814};
  time -= timeShift[i];

  // Slewing correction
  if (adc == 0) return time;

  time += 30.;
  if (adc > 300.) adc = 300.;
  if (adc > 70.) {
    return (time -
	    2.93028e+01 +
	    adc*1.25188e-02 -
	    adc*adc*2.68348e-05);
  }
  else {
    return (time -
	    3.52314e+01 +
	    adc*5.99289e-01 -
	    adc*adc*2.74668e-02 +
	    adc*adc*adc*6.61224e-04 -
	    adc*adc*adc*adc*7.77105e-06 +
	    adc*adc*adc*adc*adc*3.51229e-08);
  }
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

Int_t AliTriggerAnalysis::FMDHitCombinations(const AliESDEvent* aEsd, AliceSide side, Bool_t fillHists)
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
          
          if (fillHists)
            fHistFMDSingle->Fill(mult);
          
          if (mult > fFMDLowCut)
            totalMult = totalMult + mult;
          else
          {
            if (totalMult > fFMDHitCut)
              triggers++;
              
            if (fillHists)
              fHistFMDSum->Fill(totalMult);
              
            totalMult = 0;
          }
        }
      }
    }
  }
  
  return triggers;
}

Bool_t AliTriggerAnalysis::FMDTrigger(const AliESDEvent* aEsd, AliceSide side)
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
  const Int_t nHists = 9;
  TList collections[nHists];

  Int_t count = 0;
  while ((obj = iter->Next())) {

    AliTriggerAnalysis* entry = dynamic_cast<AliTriggerAnalysis*> (obj);
    if (entry == 0) 
      continue;

    Int_t n = 0;
    collections[n++].Add(entry->fHistV0A);
    collections[n++].Add(entry->fHistV0C);
    collections[n++].Add(entry->fHistZDC);
    collections[n++].Add(entry->fHistFMDA);
    collections[n++].Add(entry->fHistFMDC);
    collections[n++].Add(entry->fHistFMDSingle);
    collections[n++].Add(entry->fHistFMDSum);
    collections[n++].Add(entry->fHistBitsSPD);
    collections[n++].Add(entry->fHistFiredBitsSPD);

    // merge fTriggerClasses
    TIterator* iter2 = entry->fTriggerClasses->MakeIterator();
    TObjString* obj = 0;
    while ((obj = dynamic_cast<TObjString*> (iter2->Next())))
    {
      TParameter<Long64_t>* param2 = dynamic_cast<TParameter<Long64_t>*> (entry->fTriggerClasses->GetValue(obj));
      
      TParameter<Long64_t>* param1 = dynamic_cast<TParameter<Long64_t>*> (fTriggerClasses->GetValue(obj));
      if (param1)
      {
        param1->SetVal(param1->GetVal() + param2->GetVal());
      }
      else
      {
        param1 = dynamic_cast<TParameter<Long64_t>*> (param2->Clone());
        fTriggerClasses->Add(new TObjString(obj->String()), param1);
      }
    }
    
    delete iter2;
  
    count++;
  }

  Int_t n = 0;
  fHistV0A->Merge(&collections[n++]);
  fHistV0C->Merge(&collections[n++]);
  fHistZDC->Merge(&collections[n++]);
  fHistFMDA->Merge(&collections[n++]);
  fHistFMDC->Merge(&collections[n++]);
  fHistFMDSingle->Merge(&collections[n++]);
  fHistFMDSum->Merge(&collections[n++]);
  fHistBitsSPD->Merge(&collections[n++]);
  fHistFiredBitsSPD->Merge(&collections[n++]);
  
  delete iter;

  return count+1;
}

void AliTriggerAnalysis::SaveHistograms() const
{
  // write histograms to current directory
  
  if (!fHistBitsSPD)
    return;
    
  fHistBitsSPD->Write();
  fHistBitsSPD->ProjectionX();
  fHistBitsSPD->ProjectionY();
  fHistFiredBitsSPD->Write();
  fHistV0A->Write();
  fHistV0C->Write();
  fHistZDC->Write();
  fHistFMDA->Write();
  fHistFMDC->Write();
  fHistFMDSingle->Write();
  fHistFMDSum->Write();
  
  fTriggerClasses->Write("fTriggerClasses", TObject::kSingleKey);
}

void AliTriggerAnalysis::PrintTriggerClasses() const
{
  // print trigger classes
  
  Printf("Trigger Classes:");
  
  TIterator* iter = fTriggerClasses->MakeIterator();
  TObjString* obj = 0;
  while ((obj = dynamic_cast<TObjString*> (iter->Next())))
  {
    TParameter<Long64_t>* param = dynamic_cast<TParameter<Long64_t>*> (fTriggerClasses->GetValue(obj));
    
    Printf("%s: %ld triggers", obj->String().Data(), param->GetVal());
  }
}
