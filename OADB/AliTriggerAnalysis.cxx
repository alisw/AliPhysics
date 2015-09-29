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

#include "TH1F.h"
#include "TH2F.h"
#include "TList.h"
#include "TIterator.h"
#include "TParameter.h"
#include "TMap.h"
#include "TRandom.h"
#include "AliTriggerAnalysis.h"
#include "AliLog.h"
#include "AliVEvent.h"
#include "AliESDEvent.h"
#include "AliMultiplicity.h"
#include "AliESDAD.h"
#include "AliESDVZERO.h"
#include "AliESDZDC.h"
#include "AliESDFMD.h"
#include "AliESDVertex.h"
#include "AliESDtrackCuts.h"
#include "AliDAQ.h"
#include "AliESDTrdTrack.h"
#include "AliVCaloTrigger.h"

ClassImp(AliTriggerAnalysis)

AliTriggerAnalysis::AliTriggerAnalysis() :
fSPDGFOThreshold(2),
fSPDGFOEfficiency(0),
fV0TimeOffset(0),
fV0AdcThr(0),
fV0HwAdcThr(2.5),
fV0HwWinLow(61.5),
fV0HwWinHigh(86.5),
fZDCCutRefSum(-568.5),
fZDCCutRefDelta(-2.1),
fZDCCutSigmaSum(3.25),
fZDCCutSigmaDelta(2.25),
fZDCCutRefSumCorr(-65.5),
fZDCCutRefDeltaCorr(-2.1),
fZDCCutSigmaSumCorr(6.0),
fZDCCutSigmaDeltaCorr(1.2),
fZDCCutZNATimeCorrMin(0.0),
fZDCCutZNATimeCorrMax(2.0),
fZDCCutZNCTimeCorrMin(0.0),
fZDCCutZNCTimeCorrMax(5.0),
fASPDCvsTCut(65),
fBSPDCvsTCut(4),
fTRDptHSE(3.),
fTRDpidHSE(144),
fTRDptHQU(2.),
fTRDpidHQU(164.),
fTRDptHEE(3.),
fTRDpidHEE(144),
fTRDminSectorHEE(6),
fTRDmaxSectorHEE(8),
fTRDptHJT(3.),
fTRDnHJT(3),
fDoFMD(kTRUE),
fFMDLowCut(0.2),
fFMDHitCut(0.5),
fHistBitsSPD(0),
fHistFiredBitsSPD(0),
fHistSPDClsVsTrk(0),
fHistAD(0),
fHistADA(0),
fHistADC(0),
fHistV0A(0),
fHistV0C(0),
fHistZDC(0),
fHistTDCZDC(0),
fHistTimeZDC(0),
fHistTimeCorrZDC(0),
fHistFMDA(0),
fHistFMDC(0),
fHistFMDSingle(0),
fHistFMDSum(0),
fHistT0(0),
fTriggerClasses(0),
fMC(kFALSE),
fEsdTrackCuts(0),
fTPCOnly(kFALSE)
{
  // constructor
}

//-------------------------------------------------------------------------------------------------
AliTriggerAnalysis::~AliTriggerAnalysis(){
  // destructor
  if (fHistBitsSPD)      { delete fHistBitsSPD;      fHistBitsSPD = 0;      }
  if (fHistFiredBitsSPD) { delete fHistFiredBitsSPD; fHistFiredBitsSPD = 0; }
  if (fHistSPDClsVsTrk)  { delete fHistSPDClsVsTrk;  fHistSPDClsVsTrk = 0;  }
  if (fHistAD)           { delete fHistAD;           fHistAD  = 0;          }
  if (fHistADA)          { delete fHistADA;          fHistADA = 0;          }
  if (fHistADC)          { delete fHistADC;          fHistADC = 0;          }
  if (fHistV0A)          { delete fHistV0A;          fHistV0A = 0;          }
  if (fHistV0C)          { delete fHistV0C;          fHistV0C = 0;          }
  if (fHistZDC)          { delete fHistZDC;          fHistZDC = 0;          }
  if (fHistTDCZDC)       { delete fHistTDCZDC;       fHistTDCZDC = 0;       }
  if (fHistTimeZDC)      { delete fHistTimeZDC;      fHistTimeZDC = 0;      }
  if (fHistTimeCorrZDC)  { delete fHistTimeCorrZDC;  fHistTimeCorrZDC = 0;  }
  if (fHistFMDA)         { delete fHistFMDA;         fHistFMDA = 0;         }
  if (fHistFMDC)         { delete fHistFMDC;         fHistFMDC = 0;         }
  if (fHistFMDSingle)    { delete fHistFMDSingle;    fHistFMDSingle = 0;    }
  if (fHistFMDSum)       { delete fHistFMDSum;       fHistFMDSum = 0;       }
  if (fHistT0)           { delete fHistT0;           fHistT0 = 0;           }
  if (fEsdTrackCuts)     { delete fEsdTrackCuts;     fEsdTrackCuts =0;      }
  if (fTriggerClasses)   { fTriggerClasses->DeleteAll(); delete fTriggerClasses; fTriggerClasses = 0; }
}


//-------------------------------------------------------------------------------------------------
void AliTriggerAnalysis::EnableHistograms(Bool_t isLowFlux){
  // creates the monitoring histograms 
  // dynamical range of histograms can be adapted for pp and pPb via isLowFlux flag)
  // TODO check limits for FMD
  
  Int_t nBins  = isLowFlux ?  600 : 1202;
  Int_t nBinsX = isLowFlux ?  100 :  300;
  Int_t nBinsY = isLowFlux ?  500 : 1000;
  Float_t xMax = isLowFlux ?  400 : 2999.5;
  Float_t yMax = isLowFlux ? 4000 : 9999.5;
  
  // do not add these hists to the directory
  Bool_t oldStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE); 
  fHistBitsSPD      = new TH2F("fHistBitsSPD", "SPD GFO;number of fired chips (offline);number of fired chips (hardware)", nBins, -1.5, -1.5 + nBins, nBins, -1.5, -1.5+nBins);
  fHistFiredBitsSPD = new TH1F("fHistFiredBitsSPD", "SPD GFO Hardware;chip number;events", 1200, -0.5, 1199.5);
  fHistSPDClsVsTrk  = new TH2F("fHistSPDClsVsTrk", "SPD Clusters vs Tracklets; n tracklets; n clusters", nBinsX, -0.5, xMax, nBinsY, -0.5, yMax);
  fHistAD           = new TH2F("fHistAD", "ADC+ADA vs ADC+ADA;ADC-ADA time (ns);ADC+ADA time (ns)", 300, -150, 150, 300, -50, 250);
  fHistADA          = new TH1F("fHistADA", "ADA;mean time (ns);events", 2000, -100, 100);
  fHistADC          = new TH1F("fHistADC", "ADC;mean time (ns);events", 2000, -100, 100);
  fHistV0A          = new TH1F("fHistV0A", "V0A;mean time (ns);events", 400, -100, 100);
  fHistV0C          = new TH1F("fHistV0C", "V0C;mean time (ns);events", 400, -100, 100);
  fHistZDC          = new TH1F("fHistZDC", "ZDC;trigger bits;events", 8, -1.5, 6.5);
  fHistTDCZDC       = new TH1F("fHistTDCZDC", "ZDC;TDC bits;events", 32, -0.5, 32-0.5);
  fHistTimeZDC      = new TH2F("fHistTimeZDC", "ZDC;TDC timing C-A;TDC timing C+A", 120,-30,30,120,-600,-540);
  fHistTimeCorrZDC  = new TH2F("fHistTimeCorrZDC", "ZDC;Corrected TDC timing C-A; Corrected TDC timing C+A", 120,-30,30,260,-100,30);
  fHistFMDA         = new TH1F("fHistFMDA", "FMDA;combinations above threshold;events", 102, -1.5, 100.5);
  fHistFMDC         = new TH1F("fHistFMDC", "FMDC;combinations above threshold;events", 102, -1.5, 100.5);
  fHistFMDSingle    = new TH1F("fHistFMDSingle", "FMD single;multiplicity value;counts", 1000, 0, 10);
  fHistFMDSum       = new TH1F("fHistFMDSum", "FMD sum;multiplicity value;counts", 1000, 0, 10);
  fHistT0           = new TH1F("fHistT0", "T0;time (ns);events", 100, -25, 25);
  fTriggerClasses = new TMap;
  fTriggerClasses->SetOwner();
  TH1::AddDirectory(oldStatus);
}


//-------------------------------------------------------------------------------------------------
const char* AliTriggerAnalysis::GetTriggerName(Trigger trigger){
  // returns the name of the requested trigger
  // the returned string will only be valid until the next call to this function [not thread-safe]
  
  static TString str;
  
  UInt_t triggerNoFlags = (UInt_t) trigger % (UInt_t) kStartOfFlags;
  
  switch (triggerNoFlags)  {
    case kAcceptAll :      str = "ACCEPT ALL (bypass!)";      break;
    case kMB1 :            str = "MB1";                       break;
    case kMB2 :            str = "MB2";                       break;
    case kMB3 :            str = "MB3";                       break;
    case kSPDGFO :         str = "SPD GFO";                   break;
    case kSPDGFOBits :     str = "SPD GFO Bits";              break;
    case kSPDGFOL0 :       str = "SPD GFO L0 (first layer)";  break;
    case kSPDGFOL1 :       str = "SPD GFO L1 (second layer)"; break;
    case kSPDClsVsTrkBG :  str = "Cluster vs Tracklets";      break;
    case kADA :            str = "AD A BB";                   break;
    case kADC :            str = "AD C BB";                   break;
    case kADABG :          str = "AD A BG";                   break;
    case kADCBG :          str = "AD C BG";                   break;
    case kV0A :            str = "V0 A BB";                   break;
    case kV0C :            str = "V0 C BB";                   break;
    case kV0OR :           str = "V0 OR BB";                  break;
    case kV0AND :          str = "V0 AND BB";                 break;
    case kV0ABG :          str = "V0 A BG";                   break;
    case kV0CBG :          str = "V0 C BG";                   break;
    case kZDC :            str = "ZDC";                       break;
    case kZDCA :           str = "ZDC A";                     break;
    case kZDCC :           str = "ZDC C";                     break;
    case kZNA :            str = "ZN A";                      break;
    case kZNC :            str = "ZN C";                      break;
    case kZNABG :          str = "ZN A BG";                   break;
    case kZNCBG :          str = "ZN C BG";                   break;
    case kFMDA :           str = "FMD A";                     break;
    case kFMDC :           str = "FMD C";                     break;
    case kFPANY :          str = "SPD GFO | V0 | ZDC | FMD";  break;
    case kNSD1 :           str = "NSD1";                      break;
    case kMB1Prime:        str = "MB1prime";                  break;
    case kZDCTDCA :        str = "ZDC TDC A";                 break;
    case kZDCTDCC :        str = "ZDC TDC C";                 break;
    case kZDCTime :        str = "ZDC Time Cut";              break;
    case kCentral :        str = "V0 Central";                break;
    case kSemiCentral :    str = "V0 Semi-central";           break;
    case kEmcalL0 :        str = "EMCAL";                     break;
    case kTRDHCO :         str = "TRD cosmics";               break;
    case kTRDHJT :         str = "TRD Jet";                   break;
    case kTRDHSE :         str = "TRD Single Electron";       break;
    case kTRDHQU :         str = "TRD Quarkonia";             break;
    case kTRDHEE :         str = "TRD Dielectron";            break;
    default:               str = "";                          break;
  }
  if (trigger & kOfflineFlag) str += " OFFLINE";  
  if (trigger & kOneParticle) str += " OneParticle";  
  if (trigger & kOneTrack)    str += " OneTrack";  
  return str;
}


//-------------------------------------------------------------------------------------------------
Bool_t AliTriggerAnalysis::IsTriggerFired(const AliESDEvent* aEsd, Trigger trigger){
  // checks if an event has been triggered
  if (trigger & kOfflineFlag) return IsOfflineTriggerFired(aEsd, trigger);
  return IsTriggerBitFired(aEsd, trigger);
}


//-------------------------------------------------------------------------------------------------
Bool_t AliTriggerAnalysis::IsTriggerBitFired(const AliESDEvent* /*aEsd*/, Trigger /*trigger*/) const { 
  AliFatal("This IsTriggerBitFired function is obsolete.\n"); 
  return 0; 
}


//-------------------------------------------------------------------------------------------------
Bool_t AliTriggerAnalysis::IsTriggerBitFired(const AliESDEvent* aEsd, ULong64_t tclass) const {
  // Checks if corresponding bit in mask is on
  ULong64_t trigmask = aEsd->GetTriggerMask();
  return (trigmask & (1ull << (tclass-1)));
}


//-------------------------------------------------------------------------------------------------
Int_t AliTriggerAnalysis::EvaluateTrigger(const AliESDEvent* aEsd, Trigger trigger){
  // evaluates a given trigger
  // trigger combinations are not supported, for that see IsOfflineTriggerFired

  UInt_t triggerNoFlags = (UInt_t) trigger % (UInt_t) kStartOfFlags;
  Bool_t offline = trigger & kOfflineFlag;
  
  if (!offline) {
    if ( triggerNoFlags==kT0BG
      || triggerNoFlags==kT0Pileup
      || triggerNoFlags==kSPDClsVsTrkBG
      || triggerNoFlags==kZDCA
      || triggerNoFlags==kZDCC
      || triggerNoFlags==kZDCTDCA
      || triggerNoFlags==kZDCTDCC
      || triggerNoFlags==kZDCTime
      || triggerNoFlags==kZNA
      || triggerNoFlags==kZNC
      || triggerNoFlags==kZNABG
      || triggerNoFlags==kZNCBG
      || triggerNoFlags==kFMDA
      || triggerNoFlags==kFMDC
      || triggerNoFlags==kTPCLaserWarmUp
      || triggerNoFlags==kTPCHVdip
      || triggerNoFlags==kIncompleteEvent
      || triggerNoFlags==kEMCAL
      || triggerNoFlags==kEmcalL0
      || triggerNoFlags==kEmcalL1GammaHigh
      || triggerNoFlags==kEmcalL1GammaLow
      || triggerNoFlags==kEmcalL1JetHigh
      || triggerNoFlags==kEmcalL1JetLow
      || triggerNoFlags==kTRDHCO
      || triggerNoFlags==kTRDHJT
      || triggerNoFlags==kTRDHSE
      || triggerNoFlags==kTRDHQU
      || triggerNoFlags==kTRDHEE
      ) AliFatal(Form("Online trigger not available for trigger %d", triggerNoFlags));
  } else {
    if (  triggerNoFlags==kCTPV0A 
        ||triggerNoFlags==kCTPV0C
        ||triggerNoFlags==kCentral
        ||triggerNoFlags==kSemiCentral
      ) AliFatal(Form("Offline trigger not available for trigger %d", triggerNoFlags));
  }
  
  switch (triggerNoFlags) {
    case kCTPV0A:          return aEsd->GetHeader()->IsTriggerInputFired("V0A");
    case kCTPV0C:          return aEsd->GetHeader()->IsTriggerInputFired("V0C");
    case kSPDGFO:          return SPDFiredChips(aEsd, !offline, kFALSE, 0); 
    case kSPDGFOL0:        return SPDFiredChips(aEsd, !offline, kFALSE, 1);
    case kSPDGFOL1:        return SPDFiredChips(aEsd, !offline, kFALSE, 2);
    case kSPDClsVsTrkBG:   return IsSPDClusterVsTrackletBG(aEsd);
    case kADA:             return ADTrigger(aEsd, kASide, !offline) == kADBB; 
    case kADC:             return ADTrigger(aEsd, kCSide, !offline) == kADBB;
    case kADABG:           return ADTrigger(aEsd, kASide, !offline) == kADBG;
    case kADCBG:           return ADTrigger(aEsd, kCSide, !offline) == kADBG;
    case kV0A:             return V0Trigger(aEsd, kASide, !offline) == kV0BB; 
    case kV0C:             return V0Trigger(aEsd, kCSide, !offline) == kV0BB;
    case kV0ABG:           return V0Trigger(aEsd, kASide, !offline) == kV0BG;
    case kV0CBG:           return V0Trigger(aEsd, kCSide, !offline) == kV0BG;
    case kT0:              return T0Trigger(aEsd, !offline) == kT0BB;
    case kT0BG:            return T0Trigger(aEsd, !offline) == kT0DecBG;
    case kT0Pileup:        return T0Trigger(aEsd, !offline) == kT0DecPileup;
    case kZDCA:            return ZDCTrigger(aEsd, kASide);
    case kZDCC:            return ZDCTrigger(aEsd, kCSide);
    case kZDCTDCA:         return ZDCTDCTrigger(aEsd, kASide);
    case kZDCTDCC:         return ZDCTDCTrigger(aEsd, kCSide);
    case kZDCTime:         return ZDCTimeTrigger(aEsd);
    case kZNA:             return ZDCTDCTrigger(aEsd,kASide,kTRUE,kFALSE,kFALSE);
    case kZNC:             return ZDCTDCTrigger(aEsd,kCSide,kTRUE,kFALSE,kFALSE);
    case kZNABG:           return ZDCTimeBGTrigger(aEsd,kASide);
    case kZNCBG:           return ZDCTimeBGTrigger(aEsd,kCSide);
    case kFMDA:            return FMDTrigger(aEsd, kASide);
    case kFMDC:            return FMDTrigger(aEsd, kCSide);
    case kTPCLaserWarmUp:  return IsLaserWarmUpTPCEvent(aEsd);
    case kTPCHVdip:        return IsHVdipTPCEvent(aEsd);
    case kIncompleteEvent: return IsIncompleteEvent(aEsd);
    case kEMCAL:           return EMCALCellsTrigger(aEsd);
    case kEmcalL0:         return EMCALTrigger(aEsd,kEmcalL0);
    case kEmcalL1GammaHigh:return EMCALTrigger(aEsd,kEmcalL1GammaHigh);
    case kEmcalL1GammaLow: return EMCALTrigger(aEsd,kEmcalL1GammaLow);
    case kEmcalL1JetHigh:  return EMCALTrigger(aEsd,kEmcalL1JetHigh);
    case kEmcalL1JetLow:   return EMCALTrigger(aEsd,kEmcalL1JetLow);
    case kTRDHCO:          return TRDTrigger(aEsd,kTRDHCO);
    case kTRDHJT:          return TRDTrigger(aEsd,kTRDHJT);
    case kTRDHSE:          return TRDTrigger(aEsd,kTRDHSE);
    case kTRDHQU:          return TRDTrigger(aEsd,kTRDHQU);
    case kTRDHEE:          return TRDTrigger(aEsd,kTRDHEE);
    case kCentral: {
      if (!aEsd->GetVZEROData()) { AliWarning("V0 centrality trigger bits were not filled!"); return kFALSE; }
      if (!aEsd->GetVZEROData()->TestBit(AliESDVZERO::kTriggerChargeBitsFilled)) return kFALSE;
      return aEsd->GetVZEROData()->GetTriggerBits() & (1<<AliESDVZERO::kCTA2andCTC2);
    }
    case kSemiCentral: {
      if (!aEsd->GetVZEROData()) { AliWarning("V0 centrality trigger bits were not filled!"); return kFALSE; }
      if (!aEsd->GetVZEROData()->TestBit(AliESDVZERO::kTriggerChargeBitsFilled)) return kFALSE;
      return aEsd->GetVZEROData()->GetTriggerBits() & (1<<AliESDVZERO::kCTA1andCTC1);
    }
    default: AliFatal(Form("Trigger type %d not implemented", triggerNoFlags));
  }
  return kFALSE;
}


//-------------------------------------------------------------------------------------------------
Bool_t AliTriggerAnalysis::IsOfflineTriggerFired(const AliESDEvent* aEsd, Trigger trigger){
  // checks if an event has been triggered "offline"
  UInt_t triggerNoFlags = (UInt_t) trigger % (UInt_t) kStartOfFlags;
  
  Bool_t decision = kFALSE;
  switch (triggerNoFlags) {
    case kAcceptAll:        decision = kTRUE; break;
    case kMB1:              decision = SPDGFOTrigger(aEsd, 0) ||  V0Trigger(aEsd, kASide, kFALSE) == kV0BB || V0Trigger(aEsd, kCSide, kFALSE) == kV0BB;  break;
    case kMB2:              decision = SPDGFOTrigger(aEsd, 0) && (V0Trigger(aEsd, kASide, kFALSE) == kV0BB || V0Trigger(aEsd, kCSide, kFALSE) == kV0BB); break;
    case kMB3:              decision = SPDGFOTrigger(aEsd, 0) &&  V0Trigger(aEsd, kASide, kFALSE) == kV0BB && V0Trigger(aEsd, kCSide, kFALSE) == kV0BB;  break;
    case kSPDGFO:           decision = SPDGFOTrigger(aEsd, 0); break;
    case kSPDGFOBits:       decision = SPDGFOTrigger(aEsd, 1); break;
    case kADA:              decision = ADTrigger(aEsd, kASide, kFALSE) == kADBB; break;
    case kADC:              decision = ADTrigger(aEsd, kCSide, kFALSE) == kADBB; break;
    case kADABG:            decision = ADTrigger(aEsd, kASide, kFALSE) == kADBG; break;
    case kADCBG:            decision = ADTrigger(aEsd, kCSide, kFALSE) == kADBG; break;
    case kV0A:              decision = V0Trigger(aEsd, kASide, kFALSE) == kV0BB; break;
    case kV0C:              decision = V0Trigger(aEsd, kCSide, kFALSE) == kV0BB; break;
    case kV0OR:             decision = V0Trigger(aEsd, kASide, kFALSE) == kV0BB || V0Trigger(aEsd, kCSide, kFALSE) == kV0BB; break;
    case kV0AND:            decision = V0Trigger(aEsd, kASide, kFALSE) == kV0BB && V0Trigger(aEsd, kCSide, kFALSE) == kV0BB; break;
    case kV0ABG:            decision = V0Trigger(aEsd, kASide, kFALSE) == kV0BG; break;
    case kV0CBG:            decision = V0Trigger(aEsd, kCSide, kFALSE) == kV0BG; break;
    case kZDC:              decision = ZDCTrigger(aEsd, kASide) || ZDCTrigger(aEsd, kCentralBarrel) || ZDCTrigger(aEsd, kCSide); break;
    case kZDCA:             decision = ZDCTrigger(aEsd, kASide); break;
    case kZDCC:             decision = ZDCTrigger(aEsd, kCSide); break;
    case kZNA:              decision = ZDCTDCTrigger(aEsd,kASide,kTRUE,kFALSE,kFALSE); break;
    case kZNC:              decision = ZDCTDCTrigger(aEsd,kCSide,kTRUE,kFALSE,kFALSE); break;
    case kZNABG:            decision = ZDCTimeBGTrigger(aEsd,kASide); break;
    case kZNCBG:            decision = ZDCTimeBGTrigger(aEsd,kCSide); break;
    case kFMDA:             decision = FMDTrigger(aEsd, kASide); break;
    case kFMDC:             decision = FMDTrigger(aEsd, kCSide); break;
    case kEMCAL:            decision = EMCALCellsTrigger(aEsd); break;
    case kEmcalL0:          decision = EMCALTrigger(aEsd,kEmcalL0);          break;
    case kEmcalL1GammaHigh: decision = EMCALTrigger(aEsd,kEmcalL1GammaHigh); break;
    case kEmcalL1GammaLow:  decision = EMCALTrigger(aEsd,kEmcalL1GammaLow);  break;
    case kEmcalL1JetHigh:   decision = EMCALTrigger(aEsd,kEmcalL1JetHigh);   break;
    case kEmcalL1JetLow:    decision = EMCALTrigger(aEsd,kEmcalL1JetLow);    break;
    case kTRDHCO:           decision = TRDTrigger(aEsd,kTRDHCO); break;
    case kTRDHJT:           decision = TRDTrigger(aEsd,kTRDHJT); break;
    case kTRDHSE:           decision = TRDTrigger(aEsd,kTRDHSE); break;
    case kTRDHQU:           decision = TRDTrigger(aEsd,kTRDHQU); break;
    case kTRDHEE:           decision = TRDTrigger(aEsd,kTRDHEE); break;
    case kNSD1:             decision = SPDFiredChips(aEsd, 0) >= 5 || (V0Trigger(aEsd, kASide, kFALSE) == kV0BB && V0Trigger(aEsd, kCSide, kFALSE) == kV0BB); break;
    case kFPANY:            decision |= SPDGFOTrigger(aEsd, 0); 
                            decision |= V0Trigger(aEsd, kASide, kFALSE) == kV0BB;
                            decision |= V0Trigger(aEsd, kCSide, kFALSE) == kV0BB;
                            decision |= ZDCTrigger(aEsd, kASide);
                            decision |= ZDCTrigger(aEsd, kCentralBarrel);
                            decision |= ZDCTrigger(aEsd, kCSide);
                            decision |= FMDTrigger(aEsd, kASide);
                            decision |= FMDTrigger(aEsd, kCSide);
                            break; 
    case kMB1Prime:         decision |= SPDGFOTrigger(aEsd, 0) && V0Trigger(aEsd, kASide, kFALSE) == kV0BB;
                            decision |= SPDGFOTrigger(aEsd, 0) && V0Trigger(aEsd, kCSide, kFALSE) == kV0BB;
                            decision |= V0Trigger(aEsd, kASide, kFALSE) == kV0BB && V0Trigger(aEsd, kCSide, kFALSE) == kV0BB;
                            break; 
    default:                AliFatal(Form("Trigger type %d not implemented", triggerNoFlags));
  }
  
  // hadron-level requirement: SPD tracklets
  if (decision && (trigger & kOneParticle))  {
    decision = kFALSE;
    
    const AliESDVertex* vertex = aEsd->GetPrimaryVertexSPD();
    const AliMultiplicity* mult = aEsd->GetMultiplicity();
    
    if (mult && vertex && vertex->GetNContributors() > 0 && (!vertex->IsFromVertexerZ() || vertex->GetDispersion() < 0.02) && TMath::Abs(vertex->GetZ()) < 5.5)    {
      for (Int_t i=0; i<mult->GetNumberOfTracklets(); ++i) {
        if (TMath::Abs(mult->GetEta(i)) < 1) {
          decision = kTRUE;
          break;
        }
      }
    }
  }
  
  // hadron level requirement: TPC tracks
  if (decision && (trigger & kOneTrack)) {
    decision = kFALSE;
    const AliESDVertex* vertex =0x0;
    vertex = aEsd->GetPrimaryVertexTracks();
    if (!vertex || vertex->GetNContributors() <= 0) vertex = aEsd->GetPrimaryVertexSPD();
    Float_t ptmin, ptmax;
    fEsdTrackCuts->GetPtRange(ptmin,ptmax);
    AliDebug(3, Form("ptmin = %f, ptmax = %f\n",ptmin, ptmax));
    
    if (vertex && vertex->GetNContributors() > 0 && (!vertex->IsFromVertexerZ() || vertex->GetDispersion() < 0.02) && TMath::Abs(vertex->GetZ()) < 10.) {
      AliDebug(3,Form("Check on the vertex passed\n"));
      for (Int_t i=0; i<aEsd->GetNumberOfTracks(); ++i){
        if (fTPCOnly == kFALSE){
          if (fEsdTrackCuts->AcceptTrack(aEsd->GetTrack(i))){
            AliDebug(2, Form("pt of track = %f --> check passed\n",aEsd->GetTrack(i)->Pt()));
            decision = kTRUE;
            break;
          }
        }
        else {
          // TPC only tracks
          AliESDtrack *tpcTrack = fEsdTrackCuts->GetTPCOnlyTrack((AliESDEvent*)aEsd, i);
          if (!tpcTrack){
            AliDebug(3,Form("track %d is NOT a TPC track",i));
            continue;
          }
          else{
            AliDebug(3,Form("track %d IS a TPC track",i));
            if (!(fEsdTrackCuts->AcceptTrack(tpcTrack))) {
              AliDebug(2, Form("TPC track %d NOT ACCEPTED, pt = %f, eta = %f",i,tpcTrack->Pt(), tpcTrack->Eta()));
              delete tpcTrack; tpcTrack = 0x0;
              continue;
            } // end if the TPC track is not accepted
            else{
              AliDebug(2, Form("TPC track %d ACCEPTED, pt = %f, eta = %f",i,tpcTrack->Pt(), tpcTrack->Eta()));
              decision = kTRUE;
              delete tpcTrack; tpcTrack = 0x0;
              break;
            } // end if the TPC track is accepted
          } // end if it is a TPC track
        } // end if you are looking at TPC only tracks			
      } // end loop on tracks
    } // end check on vertex
    else{
      AliDebug(4,Form("Check on the vertex not passed\n"));
      for (Int_t i=0; i<aEsd->GetNumberOfTracks(); ++i){
        if (fEsdTrackCuts->AcceptTrack(aEsd->GetTrack(i))){
          AliDebug(4,Form("pt of track = %f --> check would be passed if the vertex was ok\n",aEsd->GetTrack(i)->Pt()));
          break;
        }
      }
    }
    if (!decision) AliDebug(3,("Check for kOneTrack NOT passed\n"));
  }
  
  return decision;
}



//-------------------------------------------------------------------------------------------------
Int_t AliTriggerAnalysis::SPDFiredChips(const AliESDEvent* aEsd, Int_t origin, Bool_t fillHists, Int_t layer){
  // returns the number of fired chips in the SPD
  //
  // origin = 0 --> aEsd->GetMultiplicity()->GetNumberOfFiredChips() (filled from clusters)
  // origin = 1 --> aEsd->GetMultiplicity()->TestFastOrFiredChips() (from hardware bits)
  // layer  = 0 --> both layers
  // layer  = 1 --> inner
  // layer  = 2 --> outer
  
  const AliMultiplicity* mult = aEsd->GetMultiplicity();
  if (!mult) { AliFatal("AliMultiplicity not available"); return 0; }
  
  if (origin == 0) {
    if (layer == 0) return mult->GetNumberOfFiredChips(0) + mult->GetNumberOfFiredChips(1);
    return mult->GetNumberOfFiredChips(layer-1); 
  }
  
  if (origin == 1) {
    Int_t nChips = 0;
    Int_t firstChip = 0;
    Int_t lastChip  = 1200;
    if(layer == 1) lastChip  = 400;
    if(layer == 2) firstChip = 400;

    for (Int_t i=firstChip; i<lastChip; i++) {
      if (mult->TestFastOrFiredChips(i)) {
        // efficiency simulation (if enabled)
        if (fSPDGFOEfficiency) if (gRandom->Uniform() > fSPDGFOEfficiency->GetBinContent(i+1)) continue;
        nChips++;
        if (fillHists) fHistFiredBitsSPD->Fill(i);
      }
    }
    return nChips;
  }
  
  return -1;
}


//-------------------------------------------------------------------------------------------------
Int_t AliTriggerAnalysis::SSDClusters(const AliVEvent* event){ 
  return event->GetNumberOfITSClusters(4)+event->GetNumberOfITSClusters(5); 
}


//-------------------------------------------------------------------------------------------------
AliTriggerAnalysis::ADDecision AliTriggerAnalysis::ADTrigger(const AliESDEvent* aEsd, AliceSide side, Bool_t online, Bool_t fillHists){
  // See comments on V0Trigger
  
  AliESDAD* esdAD = aEsd->GetADData();
  if (!esdAD) { AliError("AliESDAD not available");  return kADInvalid; }
  if (side != kASide && side != kCSide) return kADInvalid;
  
  AliDebug(2,Form("In ADTrigger: %f %f",esdAD->GetADATime(),esdAD->GetADCTime()));

  if (fillHists) {
    if (fHistAD) fHistAD->Fill(esdAD->GetADCTime()-esdAD->GetADATime(),esdAD->GetADATime()+esdAD->GetADCTime());
    if (side == kASide && fHistADA) fHistADA->Fill(esdAD->GetADATime());
    if (side == kCSide && fHistADC) fHistADC->Fill(esdAD->GetADCTime());
  }
  
  if (online) {
    UShort_t bits = esdAD->GetTriggerBits();
    if (side==kASide && (bits & 1<<12)) return kADBB;
    if (side==kCSide && (bits & 1<<13)) return kADBB;
    if (side==kASide && (bits & 1<< 3)) return kADBG;
    if (side==kCSide && (bits & 1<< 5)) return kADBG;
    return kADEmpty;
  } else {
    if      (side == kASide) return (ADDecision) esdAD->GetADADecision();
    else if (side == kCSide) return (ADDecision) esdAD->GetADCDecision();
  }
  
  return kADEmpty;
}



//-------------------------------------------------------------------------------------------------
AliTriggerAnalysis::V0Decision AliTriggerAnalysis::V0Trigger(const AliESDEvent* aEsd, AliceSide side, Bool_t online, Bool_t fillHists){
  // Returns the V0 trigger decision in V0A | V0C
  //
  // Returns kV0Fake if the calculated average time is in a window where neither BB nor BG is expected. 
  // The rate of such triggers can be used to estimate the background. Note that the rate has to be 
  // rescaled with the size of the windows (numerical values see below in the code)
  //
  // argument 'online' is used as a switch between online and offline trigger algorithms
  //
  // Based on an algorithm by Cvetan Cheshkov
  
  AliESDVZERO* esdV0 = aEsd->GetVZEROData();
  if (!esdV0) { AliError("AliESDVZERO not available");  return kV0Invalid; }
  if (side != kASide && side != kCSide) return kV0Invalid;
  
  AliDebug(2,Form("In V0Trigger: %f %f",esdV0->GetV0ATime(),esdV0->GetV0CTime()));
  
  Int_t begin = (side == kASide) ? 32 :  0;
  Int_t end   = (side == kASide) ? 64 : 32;
  
  if (esdV0->TestBit(AliESDVZERO::kDecisionFilled)) {
    if (online) {
      if (esdV0->TestBit(AliESDVZERO::kOnlineBitsFilled)) {
        for (Int_t i = begin; i < end; ++i) if (esdV0->GetBBFlag(i)) return kV0BB;
        for (Int_t i = begin; i < end; ++i) if (esdV0->GetBGFlag(i)) return kV0BG;
        return kV0Empty;
      }
      else {
        AliWarning("V0 online trigger analysis is not yet available!");
        return kV0BB;
      }
    }
    else {
      if (fillHists) {
        if (side == kASide && fHistV0A) fHistV0A->Fill(esdV0->GetV0ATime());
        if (side == kCSide && fHistV0C) fHistV0C->Fill(esdV0->GetV0CTime());
      }
      if      (side == kASide) return (V0Decision) esdV0->GetV0ADecision();
      else if (side == kCSide) return (V0Decision) esdV0->GetV0CDecision();
    }
  }
  
  Float_t time = 0;
  Float_t weight = 0;
  if (fMC) {
    Int_t runRange;
    if      (aEsd->GetRunNumber() <= 104803) runRange = 0;
    else if (aEsd->GetRunNumber() <= 104876) runRange = 1;
    else runRange = 2;
    
    Float_t factors[3][64] = {
      /*104792-104803*/ {4.6,5.9,6.3,6.0,4.7,5.9,4.9,5.4,4.8,4.1,4.9,4.6,4.5,5.5,5.1,5.8,4.3,4.0,4.0,3.3,3.1,2.9,3.0,5.6,3.3,4.9,3.9,5.3,4.1,4.4,3.9,5.5,5.7,9.5,5.1,5.3,6.6,7.1,8.9,4.4,4.1,5.9,9.0,4.5,4.1,6.0,4.7,7.1,4.2,4.7,3.9,6.3,5.9,4.8,4.7,4.5,4.7,5.4,5.8,5.0,5.1,5.9,5.3,3.6},
      /*104841-104876*/ {4.6,4.8,4.9,4.8,4.3,4.9,4.4,4.5,4.6,5.0,4.7,4.6,4.7,4.6,4.6,5.5,4.7,4.5,4.7,5.0,6.5,7.6,5.3,4.9,5.5,4.8,4.6,4.9,4.5,4.5,4.6,4.9,5.7,9.8,4.9,5.2,7.1,7.1,8.1,4.4,4.0,6.0,8.3,4.6,4.2,5.6,4.6,6.4,4.4,4.7,4.5,6.5,6.0,4.7,4.5,4.4,4.8,5.5,5.9,5.3,5.0,5.7,5.1,3.6},
      /*104890-104892*/ {4.7,5.2,4.8,5.0,4.4,5.0,4.4,4.6,4.6,4.5,4.4,4.6,4.5,4.6,4.8,5.5,4.8,4.5,4.4,4.3,5.4,7.7,5.6,5.0,5.4,4.3,4.5,4.8,4.5,4.5,4.6,5.3,5.7,9.6,4.9,5.4,6.1,7.2,8.6,4.4,4.0,5.4,8.8,4.4,4.2,5.8,4.7,6.7,4.3,4.7,4.0,6.1,6.0,4.9,4.8,4.6,4.7,5.2,5.7,5.0,5.0,5.8,5.3,3.6}
    };
    Float_t dA = 77.4 - 11.0;
    Float_t dC = 77.4 - 2.9;
    // Time misalignment
    Float_t timeShift[64] = {0.477957 , 0.0889999 , 0.757669 , 0.205439 , 0.239666 , -0.183705 , 0.442873 , -0.281366 , 0.260976 , 0.788995 , 0.974758 , 0.548532 , 0.495023 , 0.868472 , 0.661167 , 0.358307 , 0.221243 , 0.530179 , 1.26696 , 1.33082 , 1.27086 , 1.77133 , 1.10253 , 0.634806 , 2.14838 , 1.50212 , 1.59253 , 1.66122 , 1.16957 , 1.52056 , 1.47791 , 1.81905 , -1.94123 , -1.29124 , -2.16045 , -1.78939 , -3.11111 , -1.87178 , -1.57671 , -1.70311 , -1.81208 , -1.94475 , -2.53058 , -1.7042 , -2.08109 , -1.84416 , -0.61073 , -1.77145 , 0.16999 , -0.0585339 , 0.00401133 , 0.397726 , 0.851111 , 0.264187 , 0.59573 , -0.158263 , 0.584362 , 1.20835 , 0.927573 , 1.13895 , 0.64648 , 2.18747 , 1.68909 , 0.451194};
    Float_t dA2 = 2.8, dC2 = 3.3;
    
    if (online) {
      for (Int_t i = begin; i < end; ++i) {
        Float_t tempAdc = esdV0->GetAdc(i)/factors[runRange][i];
        Float_t tempTime = (i >= 32) ? esdV0->GetTime(i)+dA+timeShift[i]+dA2 : esdV0->GetTime(i)+dC+timeShift[i]+dC2;
        if (esdV0->GetTime(i) >= 1e-6 && tempTime > fV0HwWinLow && tempTime < fV0HwWinHigh && tempAdc > fV0HwAdcThr) return kV0BB;
      }
      return kV0Empty;
    }
    else {
      for (Int_t i = begin; i < end; ++i) {
        Float_t tempAdc = esdV0->GetAdc(i)/factors[runRange][i];
        Float_t tempTime = (i >= 32) ? esdV0->GetTime(i)+dA : esdV0->GetTime(i)+dC;
        Float_t tempRawTime = (i >= 32) ? esdV0->GetTime(i)+dA+timeShift[i]+dA2 : esdV0->GetTime(i)+dC+timeShift[i]+dC2;
        if (esdV0->GetTime(i) >= 1e-6 && tempRawTime < 125.0 && tempAdc > fV0AdcThr) {
          weight += 1.0;
          time += tempTime;
        }
      }
    }
  }
  else {
    if (online) {
      for (Int_t i = begin; i < end; ++i) {
        if (esdV0->GetTime(i) >= 1e-6 && esdV0->GetTime(i) > fV0HwWinLow && esdV0->GetTime(i) < fV0HwWinHigh && esdV0->GetAdc(i) > fV0HwAdcThr) return kV0BB;
      }
      return kV0Empty;
    }
    else {
      for (Int_t i = begin; i < end; ++i) {
        if (esdV0->GetTime(i) > 1e-6 && esdV0->GetAdc(i) > fV0AdcThr) {
          Float_t correctedTime = V0CorrectLeadingTime(i, esdV0->GetTime(i), esdV0->GetAdc(i),aEsd->GetRunNumber());
          Float_t timeWeight = V0LeadingTimeWeight(esdV0->GetAdc(i));
          time += correctedTime*timeWeight;
          weight += timeWeight;
        }
      }
    }
  }
  
  if (weight > 0) time /= weight;
  time += fV0TimeOffset;
  
  if (fillHists) {
    if (side == kASide && fHistV0A) fHistV0A->Fill(time);
    if (side == kCSide && fHistV0C) fHistV0C->Fill(time);
  }
  
  if (side == kASide) {
    if (time > 68 && time < 100)  return kV0BB;
    if (time > 54 && time < 57.5) return kV0BG;
    if (time > 57.5 && time < 68) return kV0Fake;
  }
  if (side == kCSide) {
    if (time > 75.5 && time < 100) return kV0BB;
    if (time > 69.5 && time < 73)  return kV0BG; 
    if (time > 55 && time < 69.5)  return kV0Fake;
  }
  
  return kV0Empty;
}


//-------------------------------------------------------------------------------------------------
Float_t AliTriggerAnalysis::V0CorrectLeadingTime(Int_t i, Float_t time, Float_t adc, Int_t runNumber) const {
  // Correct for slewing and align the channels
  // Authors: Cvetan Cheshkov / Raphael Tieulent
  if (time == 0) return 0;

  // Time alignment
  Float_t timeShift[64] = {0.477957 , 0.0889999 , 0.757669 , 0.205439 , 0.239666 , -0.183705 , 0.442873 , -0.281366 , 0.260976 , 0.788995 , 0.974758 , 0.548532 , 0.495023 , 0.868472 , 0.661167 , 0.358307 , 0.221243 , 0.530179 , 1.26696 , 1.33082 , 1.27086 , 1.77133 , 1.10253 , 0.634806 , 2.14838 , 1.50212 , 1.59253 , 1.66122 , 1.16957 , 1.52056 , 1.47791 , 1.81905 , -1.94123 , -1.29124 , -2.16045 , -1.78939 , -3.11111 , -1.87178 , -1.57671 , -1.70311 , -1.81208 , -1.94475 , -2.53058 , -1.7042 , -2.08109 , -1.84416 , -0.61073 , -1.77145 , 0.16999 , -0.0585339 , 0.00401133 , 0.397726 , 0.851111 , 0.264187 , 0.59573 , -0.158263 , 0.584362 , 1.20835 , 0.927573 , 1.13895 , 0.64648 , 2.18747 , 1.68909 , 0.451194};
  
  if(runNumber < 106031) time -= timeShift[i];
  
  // Slewing correction
  if (adc == 0) return time;
  
  Double_t p1 = 1.57345e+1;
  Double_t p2 =-4.25603e-1;
  
  if(runNumber >= 106031) adc *= (2.5/4.0);
  return (time - p1*TMath::Power(adc,p2));
}


//-------------------------------------------------------------------------------------------------
Float_t AliTriggerAnalysis::V0LeadingTimeWeight(Float_t adc) const {
  if (adc < 1e-6) return 0;
  
  Float_t p1 = 40.211;
  Float_t p2 =-4.25603e-1;
  Float_t p3 = 0.5646;
  
  return 1./(p1*p1*TMath::Power(adc,2.*(p2-1.))+p3*p3);
}


//-------------------------------------------------------------------------------------------------
Bool_t AliTriggerAnalysis::ZDCTDCTrigger(const AliESDEvent* aEsd, AliceSide side, Bool_t useZN, Bool_t useZP, Bool_t fillHists) const{
  // Returns if ZDC triggered, based on TDC information 
  
  AliESDZDC *esdZDC = aEsd->GetESDZDC();
  
  Bool_t zdcNA = kFALSE;
  Bool_t zdcNC = kFALSE;
  Bool_t zdcPA = kFALSE;
  Bool_t zdcPC = kFALSE;
  
  if (fMC) { // If it's MC, we use the energy
    Double_t minEnergy = 0;
    zdcNA = esdZDC->GetZDCN2Energy()>minEnergy;
    zdcNC = esdZDC->GetZDCN1Energy()>minEnergy;
    zdcPA = esdZDC->GetZDCP2Energy()>minEnergy;
    zdcPC = esdZDC->GetZDCP1Energy()>minEnergy;
  }
  else {
    Bool_t tdc[32] = {kFALSE};
    for(Int_t itdc=0; itdc<32; itdc++){
      for(Int_t i=0; i<4; i++) tdc[itdc] |= esdZDC->GetZDCTDCData(itdc, i)!=0;
      if(fillHists && tdc[itdc]) fHistTDCZDC->Fill(itdc);
    }
    zdcNA = tdc[esdZDC->GetZNATDCChannel()];
    zdcNC = tdc[esdZDC->GetZNCTDCChannel()];
    zdcPA = tdc[esdZDC->GetZPATDCChannel()];
    zdcPC = tdc[esdZDC->GetZPCTDCChannel()];
  }
  
  if (side == kASide) return ((useZP && zdcPA) || (useZN && zdcNA)); 
  if (side == kCSide) return ((useZP && zdcPC) || (useZN && zdcNC)); 
  return kFALSE;
}


//-------------------------------------------------------------------------------------------------
Bool_t AliTriggerAnalysis::ZDCTimeTrigger(const AliESDEvent *aEsd, Bool_t fillHists) const {
  // This method implements a selection based on the timing in both sides of zdcN
  // It can be used in order to eliminate parasitic collisions
  AliESDZDC *esdZDC = aEsd->GetESDZDC();
  if(fMC) {
    UInt_t esdFlag =  esdZDC->GetESDQuality();
    Bool_t znaFired  = (esdFlag & 0x01) == 0x01;
    Bool_t zncFired  = (esdFlag & 0x10) == 0x10;
    return znaFired | zncFired;
  }
  else {
    Int_t detChZNA  = esdZDC->GetZNATDCChannel();
    Int_t detChZNC  = esdZDC->GetZNCTDCChannel();
    Int_t detChGate = esdZDC->IsZDCTDCcablingSet() ? 20 : 14;
    
    for(Int_t i=0;i<4;++i) {
      if (esdZDC->GetZDCTDCData(detChZNC,i)==0) continue;
      Float_t tdcC = 0.025*(esdZDC->GetZDCTDCData(detChZNC,i)-esdZDC->GetZDCTDCData(detChGate,i)); 
      Float_t tdcCcorr = esdZDC->GetZDCTDCCorrected(detChZNC,i); 
      for(Int_t j=0;j<4;++j) {
        if (esdZDC->GetZDCTDCData(detChZNA,j)==0) continue;
        Float_t tdcA = 0.025*(esdZDC->GetZDCTDCData(detChZNA,j)-esdZDC->GetZDCTDCData(detChGate,j));
        Float_t tdcAcorr = esdZDC->GetZDCTDCCorrected(detChZNA,j);
        if(fillHists) {
          fHistTimeZDC->Fill(tdcC-tdcA,tdcC+tdcA);
          fHistTimeCorrZDC->Fill(tdcCcorr-tdcAcorr,tdcCcorr+tdcAcorr);
        }
        if (esdZDC->TestBit(AliESDZDC::kCorrectedTDCFilled)) {
          if (TMath::Power((tdcCcorr-tdcAcorr-fZDCCutRefDeltaCorr)/fZDCCutSigmaDeltaCorr,2.)+
              TMath::Power((tdcCcorr+tdcAcorr-fZDCCutRefSumCorr  )/fZDCCutSigmaSumCorr,2.) < 1.) return kTRUE;
        }
        else {
          if (TMath::Power((tdcC-tdcA-fZDCCutRefDelta)/fZDCCutSigmaDelta,2.)+
              TMath::Power((tdcC+tdcA-fZDCCutRefSum  )/fZDCCutSigmaSum,2.  )<1.0) return kTRUE;
        }
      }
    }
  }
  return kFALSE;
}


//-------------------------------------------------------------------------------------------------
Bool_t AliTriggerAnalysis::ZDCTimeBGTrigger(const AliESDEvent *aEsd, AliceSide side) const{
  // This method implements a selection based on the timing in zdcN
  // It can be used in order to flag background
  if(fMC) return kFALSE;
  
  AliESDZDC* zdcData = aEsd->GetESDZDC();
  Bool_t znabadhit = kFALSE;
  Bool_t zncbadhit = kFALSE;
  
  Float_t tdcCcorr=999, tdcAcorr=999;
  
  Int_t detChZNA  = zdcData->GetZNATDCChannel();
  Int_t detChZNC  = zdcData->GetZNCTDCChannel();

  for(Int_t i = 0; i < 4; ++i) {
    if (zdcData->GetZDCTDCData(detChZNC,i)==0) continue;
    tdcCcorr = TMath::Abs(zdcData->GetZDCTDCCorrected(detChZNC,i));
    if(tdcCcorr<fZDCCutZNCTimeCorrMax && tdcCcorr>fZDCCutZNCTimeCorrMin) zncbadhit = kTRUE;
  }
  for(Int_t i = 0; i < 4; ++i) {
    if (zdcData->GetZDCTDCData(detChZNA,i)==0) continue;
    tdcAcorr = TMath::Abs(zdcData->GetZDCTDCCorrected(detChZNA,i));
    if(tdcAcorr<fZDCCutZNATimeCorrMax && tdcAcorr>fZDCCutZNATimeCorrMin) znabadhit = kTRUE;
  }
  
  if (side == kASide) return znabadhit;
  if (side == kCSide) return zncbadhit;

  return kFALSE;
}


//-------------------------------------------------------------------------------------------------
Bool_t AliTriggerAnalysis::ZDCTrigger(const AliESDEvent* aEsd, AliceSide side) const {
  // Returns if ZDC triggered
  
  AliESDZDC* zdcData = aEsd->GetESDZDC();
  UInt_t quality = zdcData->GetESDQuality();

  // from Nora's presentation, general first physics meeting 16.10.09
  static UInt_t zpc  = 0x20;
  static UInt_t znc  = 0x10;
  static UInt_t zem1 = 0x08;
  static UInt_t zem2 = 0x04;
  static UInt_t zpa  = 0x02;
  static UInt_t zna  = 0x01;
  
  if (side == kASide         && ((quality & zpa)  || (quality & zna ))) return kTRUE;
  if (side == kCentralBarrel && ((quality & zem1) || (quality & zem2))) return kTRUE;
  if (side == kCSide         && ((quality & zpc)  || (quality & znc ))) return kTRUE;
  
  return kFALSE;
}


//-------------------------------------------------------------------------------------------------
Int_t AliTriggerAnalysis::FMDHitCombinations(const AliESDEvent* aEsd, AliceSide side, Bool_t fillHists){
  // returns number of hit combinations above threshold
  // Authors: FMD team, Hans Dalsgaard (code merged from FMD/AliFMDOfflineTrigger)
  if (!fDoFMD) return -1;
  
  // Workaround for AliESDEvent::GetFMDData is not const!
  const AliESDFMD* fmdData = (const_cast<AliESDEvent*>(aEsd))->GetFMDData();
  if (!fmdData) {
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
          if (fillHists) fHistFMDSingle->Fill(mult);
          if (mult > fFMDLowCut)
            totalMult = totalMult + mult;
          else {
            if (totalMult > fFMDHitCut) triggers++;
            if (fillHists) fHistFMDSum->Fill(totalMult);
            totalMult = 0;
          }
        }
      }
    }
  }
  return triggers;
}


//-------------------------------------------------------------------------------------------------
Bool_t AliTriggerAnalysis::FMDTrigger(const AliESDEvent* aEsd, AliceSide side){
  // Returns if the FMD triggered
  // Authors: FMD team, Hans Dalsgaard (code merged from FMD/AliFMDOfflineTrigger)
  return FMDHitCombinations(aEsd, side, kFALSE);
}


//-------------------------------------------------------------------------------------------------
AliTriggerAnalysis::T0Decision AliTriggerAnalysis::T0Trigger(const AliESDEvent* aEsd, Bool_t online, Bool_t fillHists){
  // Returns the T0 TVDC trigger decision
  //  
  // argument 'online' is used as a switch between online and offline trigger algorithms
  // in online mode return 0TVX 
  // in offline mode in addition check pile-up and background :
  // pile-up readed from ESD: check if TVDC (0TVX module name) has more 1 hit;
  // backgroud flag readed from ESD : check in given time interval OrA and OrC were correct but TVDC not  
  // 
  // Based on an algorithm by Alla Maevskaya 
  
  const AliESDTZERO* esdT0 = aEsd->GetESDTZERO();
  if (!esdT0) {
    AliError("AliESDTZERO not available");
    return kT0Invalid;
  }

  Float_t  tvdc0 = esdT0->GetTVDC(0);
  if(fillHists) fHistT0->Fill(tvdc0);
  
  if (online) {
    if( aEsd->GetHeader()->IsTriggerInputFired("0TVX")) return kT0BB;
  }
  else {
    if (esdT0->GetPileupFlag()) return kT0DecPileup;
    if (esdT0->GetBackgroundFlag()) return kT0DecBG;
    if (tvdc0>-5 && tvdc0<5 && tvdc0!=0) return kT0BB; 
  }
  
  if (fMC) if(esdT0->GetT0zVertex()>-12.3 && esdT0->GetT0zVertex() < 10.3) return kT0BB; 
  return kT0Empty;
}

//-------------------------------------------------------------------------------------------------
Bool_t AliTriggerAnalysis::EMCALCellsTrigger(const AliESDEvent *aEsd){
  //
  // Returns the EMCAL trigger decision
  // so far only implemented for LHC11a data
  // see http://alisoft.cern.ch/viewvc/trunk/PWGGA/EMCALTasks/AliEmcalPhysicsSelection.cxx?view=markup&root=AliRoot Revision 56136
  //
  
  Bool_t isFired = kTRUE;
  const Int_t runNumber = aEsd->GetRunNumber();
  
  /*
  // Load EMCAL branches from the manager
  AliAnalysisManager *am = AliAnalysisManager::GetAnalysisManager();
  am->LoadBranch("EMCALCells.");
  am->LoadBranch("CaloClusters");
   */
  
  // Get EMCAL cells
  AliVCaloCells *cells = aEsd->GetEMCALCells();
  const Short_t nCells = cells->GetNumberOfCells();
  
  // count cells above threshold per sm
  Int_t nCellCount[10] = {0,0,0,0,0,0,0,0,0,0};
  for(Int_t iCell=0; iCell<nCells; ++iCell) {
    Short_t cellId = cells->GetCellNumber(iCell);
    Double_t cellE = cells->GetCellAmplitude(cellId);
    Int_t sm       = cellId / (24*48);
    if (cellE>0.1)
      ++nCellCount[sm];
  }
  
  // Trigger decision for LHC11a
  Bool_t isLedEvent = kFALSE;
  if ((runNumber>=144871) && (runNumber<=146860)) {
    if (nCellCount[4] > 100)
      isLedEvent = kTRUE;
    else {
      if ((runNumber>=146858) && (runNumber<=146860)) {
        if (nCellCount[3]>=35)
          isLedEvent = kTRUE;
      }
    }
  }
  
  if (isLedEvent) {
    isFired = kFALSE;
  }
  
  return isFired;
}


//----------------------------------------------------------------------------------------------------
Bool_t AliTriggerAnalysis::TRDTrigger(const AliESDEvent *esd, Trigger trigger){
  // evaluate the TRD trigger conditions,
  // so far HCO, HSE, HQU, HJT, HEE
  if(trigger!=kTRDHCO && trigger!=kTRDHJT && trigger!=kTRDHSE && trigger!=kTRDHQU && trigger!=kTRDHEE) {
    AliWarning("Beware you are erroneously trying to use this function (wrong trigger)");
    return kFALSE;
  }
  
  Int_t nTrdTracks = esd->GetNumberOfTrdTracks();
  if (nTrdTracks<=0) return kFALSE;
  if      (trigger==kTRDHCO) return kTRUE;
  else if (trigger!=kTRDHJT) {
    for (Int_t iTrack = 0; iTrack < nTrdTracks; ++iTrack) {
      AliESDTrdTrack *trdTrack = esd->GetTrdTrack(iTrack);
      if (!trdTrack) continue;
      // for the electron triggers we only consider matched tracks
      if(trigger==kTRDHQU) if (TMath::Abs(trdTrack->Pt())>fTRDptHQU && trdTrack->GetPID()>fTRDpidHQU) return kTRUE;
      if(trigger==kTRDHSE) if (TMath::Abs(trdTrack->Pt())>fTRDptHSE && trdTrack->GetPID()>fTRDpidHSE) return kTRUE; 
      if(trigger==kTRDHEE) if (TMath::Abs(trdTrack->Pt())>fTRDptHSE && trdTrack->GetPID()>fTRDpidHSE && trdTrack->GetSector()>=fTRDminSectorHEE && trdTrack->GetSector()<=fTRDmaxSectorHEE) return kTRUE;
    }
  } 
  else if (trigger==kTRDHJT) {
    Int_t nTracks[90] = { 0 }; // stack-wise counted number of tracks above pt threshold
    for (Int_t iTrack = 0; iTrack < nTrdTracks; ++iTrack) {
      AliESDTrdTrack *trdTrack = esd->GetTrdTrack(iTrack);    
      if (!trdTrack) continue;
      Int_t globalStack = 5*trdTrack->GetSector() + trdTrack->GetStack();
      // stack-wise counting of tracks above pt threshold for jet trigger
      if (TMath::Abs(trdTrack->GetPt()) >= fTRDptHJT) ++nTracks[globalStack];
    }
    // check if HJT condition is fulfilled in any stack
    for (Int_t iStack = 0; iStack < 90; iStack++) if (nTracks[iStack] >= fTRDnHJT) return kTRUE;
  }
  return kFALSE;
}


//-------------------------------------------------------------------------------------------------
Bool_t AliTriggerAnalysis::EMCALTrigger(const AliVEvent* event, Trigger trigger){
  AliVCaloTrigger* emcalTrigger = event->GetCaloTrigger("EMCAL");
  if (!emcalTrigger) return kFALSE;
  Int_t emcalTriggerBits = 0;
  emcalTrigger->GetTriggerBits(emcalTriggerBits);
  if      (trigger==kEmcalL0         ) { return emcalTriggerBits & 1<<0; }
  else if (trigger==kEmcalL1GammaHigh) { return emcalTriggerBits & 1<<1; }
  else if (trigger==kEmcalL1GammaLow ) { return emcalTriggerBits & 1<<2; }
  else if (trigger==kEmcalL1JetHigh  ) { return emcalTriggerBits & 1<<3; }
  else if (trigger==kEmcalL1JetLow   ) { return emcalTriggerBits & 1<<4; }
  else {
    AliWarning("Beware you are erroneously trying to use this function (wrong trigger)");
    return kFALSE;
  }
  return kFALSE;
}


//-------------------------------------------------------------------------------------------------
Bool_t AliTriggerAnalysis::IsSPDClusterVsTrackletBG(const AliVEvent* event, Bool_t fillHists){
  // rejects BG based on the cluster vs tracklet correlation
  // returns true if the event is BG
  const AliVMultiplicity* mult = event->GetMultiplicity();
  if (!mult) { AliFatal("No multiplicity object"); return 0; }
  Int_t ntracklet   = mult->GetNumberOfTracklets();
  Int_t spdClusters = event->GetNumberOfITSClusters(0) + event->GetNumberOfITSClusters(1);
  if(fillHists) fHistSPDClsVsTrk->Fill(ntracklet,spdClusters);
  return spdClusters > Float_t(fASPDCvsTCut) + Float_t(ntracklet)*fBSPDCvsTCut;
}


//-------------------------------------------------------------------------------------------------
Bool_t AliTriggerAnalysis::IsLaserWarmUpTPCEvent(const AliESDEvent* esd){
  // This function flags noisy TPC events which can happen during laser warm-up.
  Int_t trackCounter = 0;
  for (Int_t i=0; i<esd->GetNumberOfTracks(); ++i) {
    AliESDtrack *track = esd->GetTrack(i);
    if (!track) continue;
    if (track->GetTPCNcls() < 30) continue;
    if (TMath::Abs(track->Eta()) > 0.005) continue;
    if (track->Pt() < 4) continue;
    if (track->GetKinkIndex(0) > 0) continue;
    UInt_t status = track->GetStatus();
    if ((status&AliESDtrack::kITSrefit)==AliESDtrack::kITSrefit) continue; // explicitly ask for tracks without ITS refit
    if ((status&AliESDtrack::kTPCrefit)!=AliESDtrack::kTPCrefit) continue;
    if (track->GetTPCsignal() > 10) continue;          // explicitly ask for tracks without dE/dx
    trackCounter++;
  }
  if (trackCounter > 15) return kTRUE;
  return kFALSE;
}


//-------------------------------------------------------------------------------------------------
Bool_t AliTriggerAnalysis::IsHVdipTPCEvent(const AliESDEvent* esd) {
  // This function flags events in which the TPC chamber HV is not at its nominal value
  if (fMC) return kFALSE; // there are no dip events in MC
  if (!esd->IsDetectorOn(AliDAQ::kTPC)) return kTRUE;
  return kFALSE;
}


//-------------------------------------------------------------------------------------------------
Bool_t AliTriggerAnalysis::IsIncompleteEvent(const AliESDEvent* esd){
  // Check whether the event is incomplete 
  // (due to DAQ-HLT issues, it could be only part of the event was saved)
  if (fMC) return kFALSE; // there are no incomplete events on MC
  if ((esd->GetEventType() == 7) &&
      (esd->GetDAQDetectorPattern() & (1<<4)) &&
      !(esd->GetDAQAttributes() & (1<<7))) return kTRUE;
  return kFALSE;
}


//-------------------------------------------------------------------------------------------------
Long64_t AliTriggerAnalysis::Merge(TCollection* list){
  // Merge a list of objects with this (needed for PROOF).
  // Returns the number of merged objects (including this).
  if (!list) return 0;
  if (list->IsEmpty()) return 1;
  TIterator* iter = list->MakeIterator();
  TObject* obj;
  
  // collections of all histograms
  const Int_t nHists = 17;
  TList collections[nHists];
  
  Int_t count = 0;
  while ((obj = iter->Next())) {
    AliTriggerAnalysis* entry = dynamic_cast<AliTriggerAnalysis*> (obj);
    if (entry == 0) continue;
    Int_t n = 0;
    collections[n++].Add(entry->fHistAD);
    collections[n++].Add(entry->fHistADA);
    collections[n++].Add(entry->fHistADC);
    collections[n++].Add(entry->fHistV0A);
    collections[n++].Add(entry->fHistV0C);
    collections[n++].Add(entry->fHistZDC);
    collections[n++].Add(entry->fHistTDCZDC);
    collections[n++].Add(entry->fHistTimeZDC);
    collections[n++].Add(entry->fHistTimeCorrZDC);
    collections[n++].Add(entry->fHistFMDA);
    collections[n++].Add(entry->fHistFMDC);
    collections[n++].Add(entry->fHistFMDSingle);
    collections[n++].Add(entry->fHistFMDSum);
    collections[n++].Add(entry->fHistBitsSPD);
    collections[n++].Add(entry->fHistFiredBitsSPD);
    collections[n++].Add(entry->fHistSPDClsVsTrk);
    collections[n++].Add(entry->fHistT0);
    
    // merge fTriggerClasses
    TIterator* iter2 = entry->fTriggerClasses->MakeIterator();
    TObjString* obj2 = 0;
    while ((obj2 = dynamic_cast<TObjString*> (iter2->Next()))) {
      TParameter<Long64_t>* param2 = static_cast<TParameter<Long64_t>*> (entry->fTriggerClasses->GetValue(obj2));
      TParameter<Long64_t>* param1 = dynamic_cast<TParameter<Long64_t>*> (fTriggerClasses->GetValue(obj2));
      if (param1) { param1->SetVal(param1->GetVal() + param2->GetVal()); }
      else {        
        param1 = dynamic_cast<TParameter<Long64_t>*> (param2->Clone());
        fTriggerClasses->Add(new TObjString(obj2->String()), param1);
      }
    }
    delete iter2;
    count++;
  }
  
  Int_t n = 0;
  fHistAD->Merge(&collections[n++]);
  fHistADA->Merge(&collections[n++]);
  fHistADC->Merge(&collections[n++]);
  fHistV0A->Merge(&collections[n++]);
  fHistV0C->Merge(&collections[n++]);
  fHistZDC->Merge(&collections[n++]);
  fHistTDCZDC->Merge(&collections[n++]);
  if (fHistTimeZDC)     fHistTimeZDC->Merge(&collections[n++]);     else n++;
  if (fHistTimeCorrZDC) fHistTimeCorrZDC->Merge(&collections[n++]); else n++;
  fHistFMDA->Merge(&collections[n++]);
  fHistFMDC->Merge(&collections[n++]);
  fHistFMDSingle->Merge(&collections[n++]);
  fHistFMDSum->Merge(&collections[n++]);
  fHistBitsSPD->Merge(&collections[n++]);
  fHistFiredBitsSPD->Merge(&collections[n++]);
  fHistSPDClsVsTrk->Merge(&collections[n++]);
  fHistT0->Merge(&collections[n++]);
  delete iter;
  return count+1;
}


//-------------------------------------------------------------------------------------------------
void AliTriggerAnalysis::FillHistograms(const AliESDEvent* aEsd){
  // fills the histograms with the info from the ESD
  
  fHistBitsSPD->Fill(SPDFiredChips(aEsd, 0), SPDFiredChips(aEsd, 1, kTRUE));
  
  ADTrigger(aEsd, kASide, kFALSE, kTRUE);
  ADTrigger(aEsd, kCSide, kFALSE, kTRUE);
  V0Trigger(aEsd, kASide, kFALSE, kTRUE);
  V0Trigger(aEsd, kCSide, kFALSE, kTRUE);
  T0Trigger(aEsd, kFALSE, kTRUE);
  ZDCTDCTrigger(aEsd,kASide,kFALSE,kFALSE,kTRUE);
  ZDCTimeTrigger(aEsd,kTRUE);
  IsSPDClusterVsTrackletBG(aEsd, kTRUE);
  
  AliESDZDC* zdcData = aEsd->GetESDZDC();
  if (zdcData)  {
    UInt_t quality = zdcData->GetESDQuality();
    
    // from Nora's presentation, general first physics meeting 16.10.09
    static UInt_t zpc  = 0x20;
    static UInt_t znc  = 0x10;
    static UInt_t zem1 = 0x08;
    static UInt_t zem2 = 0x04;
    static UInt_t zpa  = 0x02;
    static UInt_t zna  = 0x01;
    
    fHistZDC->Fill(1, (quality & zna)  ? 1 : 0);
    fHistZDC->Fill(2, (quality & zpa)  ? 1 : 0);
    fHistZDC->Fill(3, (quality & zem2) ? 1 : 0);
    fHistZDC->Fill(4, (quality & zem1) ? 1 : 0);
    fHistZDC->Fill(5, (quality & znc)  ? 1 : 0);
    fHistZDC->Fill(6, (quality & zpc)  ? 1 : 0);
  }
  else {
    fHistZDC->Fill(-1);
    AliError("AliESDZDC not available");
  }
  
  if (fDoFMD) {
    fHistFMDA->Fill(FMDHitCombinations(aEsd, kASide, kTRUE));
    fHistFMDC->Fill(FMDHitCombinations(aEsd, kCSide, kTRUE));
  }
}


//-------------------------------------------------------------------------------------------------
void AliTriggerAnalysis::SaveHistograms() const {
  // write histograms to current directory
  if (fHistBitsSPD)      fHistBitsSPD->Write();
  if (fHistFiredBitsSPD) fHistFiredBitsSPD->Write();
  if (fHistAD )          fHistAD->Write();
  if (fHistADA)          fHistADA->Write();
  if (fHistADC)          fHistADC->Write();
  if (fHistV0A)          fHistV0A->Write();
  if (fHistV0C)          fHistV0C->Write();
  if (fHistZDC)          fHistZDC->Write();
  if (fHistTDCZDC)       fHistTDCZDC->Write();
  if (fHistTimeZDC)      fHistTimeZDC->Write();
  if (fHistTimeCorrZDC)  fHistTimeCorrZDC->Write();
  if (fHistFMDA)         fHistFMDA->Write();
  if (fHistFMDC)         fHistFMDC->Write();
  if (fHistFMDSingle)    fHistFMDSingle->Write();
  if (fHistFMDSum)       fHistFMDSum->Write();
  if (fSPDGFOEfficiency) fSPDGFOEfficiency->Write("fSPDGFOEfficiency");
  if (fHistSPDClsVsTrk)  fHistSPDClsVsTrk->Write("fHistSPDClsVsTrk");
  if (fHistT0)           fHistT0->Write("fHistT0");
  fTriggerClasses->Write("fTriggerClasses", TObject::kSingleKey);
}


//-------------------------------------------------------------------------------------------------
void AliTriggerAnalysis::FillTriggerClasses(const AliESDEvent* aEsd){
  // fills trigger classes map
  TParameter<Long64_t>* count = dynamic_cast<TParameter<Long64_t>*> (fTriggerClasses->GetValue(aEsd->GetFiredTriggerClasses().Data()));
  if (!count) {
    count = new TParameter<Long64_t>(aEsd->GetFiredTriggerClasses(), 0);
    fTriggerClasses->Add(new TObjString(aEsd->GetFiredTriggerClasses().Data()), count);
  }
  count->SetVal(count->GetVal() + 1);
}


//-------------------------------------------------------------------------------------------------
void AliTriggerAnalysis::PrintTriggerClasses() const {
  // print trigger classes
  
  Printf("Trigger Classes:");
  Printf("Event count for trigger combinations:");
  TMap singleTrigger;
  singleTrigger.SetOwner();
  TIterator* iter = fTriggerClasses->MakeIterator();
  TObjString* obj = 0;
  while ((obj = dynamic_cast<TObjString*> (iter->Next()))) {
    TParameter<Long64_t>* param = static_cast<TParameter<Long64_t>*> (fTriggerClasses->GetValue(obj));
    Printf(" %s: %ld triggers", obj->String().Data(), (Long_t)param->GetVal());
    TObjArray* tokens = obj->String().Tokenize(" ");
    for (Int_t i=0; i<tokens->GetEntries(); i++) {
      TParameter<Long64_t>* count = dynamic_cast<TParameter<Long64_t>*> (singleTrigger.GetValue(((TObjString*) tokens->At(i))->String().Data()));
      if (!count) {
        count = new TParameter<Long64_t>(((TObjString*) tokens->At(i))->String().Data(), 0);
        singleTrigger.Add(new TObjString(((TObjString*) tokens->At(i))->String().Data()), count);
      }
      count->SetVal(count->GetVal() + param->GetVal());
    }
    delete tokens;
  }
  delete iter;
  
  Printf("Event count for single trigger:");
  iter = singleTrigger.MakeIterator();
  while ((obj = dynamic_cast<TObjString*> (iter->Next()))) {
    TParameter<Long64_t>* param = static_cast<TParameter<Long64_t>*> (singleTrigger.GetValue(obj));
    Printf("  %s: %ld triggers", obj->String().Data(), (Long_t)param->GetVal());
  }
  delete iter;
  singleTrigger.DeleteAll();
}
