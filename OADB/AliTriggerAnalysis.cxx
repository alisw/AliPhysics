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
// This class provides function to check if events have been triggered based 
// on the data in ESD and AODs. The trigger bits, trigger class inputs and 
// only the data (offline trigger) can be used
// Origin: Jan Fiete Grosse-Oetringhaus, CERN
// Current support and development: Evgeny Kryshen, PNPI
//-------------------------------------------------------------------------

#include "TF1.h"
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
#include "AliAODTZERO.h"
#include "AliAODEvent.h"
ClassImp(AliTriggerAnalysis)

AliTriggerAnalysis::AliTriggerAnalysis(TString name) :
AliOADBTriggerAnalysis(name.Data()),
fSPDGFOEfficiency(0),
fDoFMD(kFALSE),
fHistList(0),
fHistFiredBitsSPD(0),
fHistSPDClsVsTkl(0),
fHistV0MOnVsOf(0),
fHistSPDOnVsOf(0),
fHistV0MOn(0),
fHistV0MOfAll(0),
fHistV0MOfAcc(0),
fHistSPDOnOuter(0),
fHistV0C3vs012(0),
fHistVIRvsBCmod4pup(0),
fHistVIRvsBCmod4acc(0),
fHistBBAflags(0),
fHistBBCflags(0),
fHistBGAflags(0),
fHistBGCflags(0),
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
fMC(kFALSE)
{
  // constructor
}

void AliTriggerAnalysis::SetParameters(AliOADBTriggerAnalysis* oadb){
  fZDCCutRefSumCorr     = oadb->GetZDCCutRefSumCorr();
  fZDCCutRefDeltaCorr   = oadb->GetZDCCutRefDeltaCorr();
  fZDCCutSigmaSumCorr   = oadb->GetZDCCutSigmaSumCorr();
  fZDCCutSigmaDeltaCorr = oadb->GetZDCCutSigmaDeltaCorr();
  fZDCCutZNATimeCorrMax = oadb->GetZDCCutZNATimeCorrMax();
  fZDCCutZNATimeCorrMin = oadb->GetZDCCutZNATimeCorrMin();
  fZDCCutZNCTimeCorrMax = oadb->GetZDCCutZNCTimeCorrMax();
  fZDCCutZNCTimeCorrMin = oadb->GetZDCCutZNCTimeCorrMin();
  fV0MOnVsOfA           = oadb->GetV0MOnVsOfA();
  fV0MOnVsOfB           = oadb->GetV0MOnVsOfB();
  fSPDOnVsOfA           = oadb->GetSPDOnVsOfA();
  fSPDOnVsOfB           = oadb->GetSPDOnVsOfB();
  fV0CasymA             = oadb->GetV0CasymA();
  fV0CasymB             = oadb->GetV0CasymB();
  fNBCsPast             = oadb->GetNBCsPast();
  fNBCsFuture           = oadb->GetNBCsFuture();
  fVIRBBAflags          = oadb->GetVIRBBAflags();
  fVIRBBCflags          = oadb->GetVIRBBCflags();
  fVIRBGAflags          = oadb->GetVIRBGAflags();
  fVIRBGCflags          = oadb->GSetVIRBGCflags();
  fVHMBBAflags          = oadb->GetVHMBBAflags();
  fVHMBBCflags          = oadb->GetVHMBBCflags();
  fVHMBGAflags          = oadb->GetVHMBGAflags();
  fVHMBGCflags          = oadb->GetVHMBGCflags();
  fV0MOnThreshold       = oadb->GetV0MOnThreshold();
  fV0MOfThreshold       = oadb->GetV0MOfThreshold();
  fSPDGFOThreshold      = oadb->GetSPDGFOThreshhold();
  fSH1OuterThreshold    = oadb->GetSH1OuterThreshold();
  fSH2OuterThreshold    = oadb->GetSH2OuterThreshold();
  fFMDLowCut            = oadb->GetFMDLowThreshold();
  fFMDHitCut            = oadb->GetFMDHitThreshold();
}

//-------------------------------------------------------------------------------------------------
AliTriggerAnalysis::~AliTriggerAnalysis(){
  // destructor
  if (fHistList)           { fHistList->Delete();          delete fHistList;       fHistList       = 0; }
  if (fTriggerClasses)     { fTriggerClasses->DeleteAll(); delete fTriggerClasses; fTriggerClasses = 0; }
}


//-------------------------------------------------------------------------------------------------
void AliTriggerAnalysis::EnableHistograms(Bool_t isLowFlux){
  // creates the monitoring histograms 
  // dynamical range of histograms can be adapted for pp and pPb via isLowFlux flag)
  // TODO check limits for FMD
  
  Int_t nBinsX = isLowFlux ?  400 :  300;
  Int_t nBinsY = isLowFlux ? 4000 : 1000;
  Float_t xMax = isLowFlux ?  400 : 2999.5;
  Float_t yMax = isLowFlux ? 4000 : 9999.5;
  
  // do not add these hists to the directory
  Bool_t oldStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);
  fHistFiredBitsSPD   = new TH1F("fHistFiredBitsSPD", "SPD GFO Hardware;chip number;events", 1200, -0.5, 1199.5);
  fHistSPDClsVsTkl    = new TH2F("fHistSPDClsVsTkl", "SPD Clusters vs Tracklets; n tracklets; n clusters", nBinsX, -0.5, xMax, nBinsY, -0.5, yMax);
  fHistV0MOnVsOf      = new TH2F("fHistV0MOnvsOf",";Offline V0M;Online V0M",1000,0,2000,1000,0,10000);
  fHistSPDOnVsOf      = new TH2F("fHistV0MOnvsOf",";Offline FOR;Online FOR",800,0,800,800,0,800);
  fHistVIRvsBCmod4pup = new TH2F("fHistVIRvsBCmod4pup","VIR vs BC%4 for events identified as SPD or V0 pileup;VIR;BC%4",21,-10.5,10.5,4,-0.5,3.5);
  fHistVIRvsBCmod4acc = new TH2F("fHistVIRvsBCmod4acc","VIR vs BC%4 for accepted events;VIR;BC%4",21,-10.5,10.5,4,-0.5,3.5);
  fHistV0C3vs012      = new TH2F("fHistV0C3vs012",";V0C012 multiplicity;V0C3 multiplicity",800,0,800,300,0,300);
  fHistBBAflags       = new TH1F("fHistBBAflags",";BBA flags;",33,-0.5,32.5);
  fHistBBCflags       = new TH1F("fHistBBCflags",";BBC flags;",33,-0.5,32.5);
  fHistBGAflags       = new TH1F("fHistBGAflags",";BGA flags;",33,-0.5,32.5);
  fHistBGCflags       = new TH1F("fHistBGCflags",";BGC flags;",33,-0.5,32.5);
  fHistV0MOn          = new TH1F("fHistV0MOn",";Online V0M;",10000,0,10000);
  fHistV0MOfAll       = new TH1F("fHistV0MOfAll","All;Offline V0M;",2000,0,2000);
  fHistV0MOfAcc       = new TH1F("fHistV0MOfAcc","Accepted;Offline V0M;",2000,0,2000);
  fHistSPDOnOuter     = new TH1F("fHistSPDOnOuter","Online outer FO chips",800,0,800);
  fHistAD             = new TH2F("fHistAD", "ADC+ADA vs ADC+ADA;ADC-ADA time (ns);ADC+ADA time (ns)", 300, -150, 150, 300, -50, 250);
  fHistADA            = new TH1F("fHistADA", "ADA;mean time (ns);events", 2000, -100, 100);
  fHistADC            = new TH1F("fHistADC", "ADC;mean time (ns);events", 2000, -100, 100);
  fHistV0A            = new TH1F("fHistV0A", "V0A;mean time (ns);events", 400, -100, 100);
  fHistV0C            = new TH1F("fHistV0C", "V0C;mean time (ns);events", 400, -100, 100);
  fHistZDC            = new TH1F("fHistZDC", "ZDC;trigger bits;events", 8, -1.5, 6.5);
  fHistTDCZDC         = new TH1F("fHistTDCZDC", "ZDC;TDC bits;events", 32, -0.5, 32-0.5);
  fHistTimeZDC        = new TH2F("fHistTimeZDC", "ZDC;TDC timing C-A;TDC timing C+A", 120,-30,30,120,-600,-540);
  fHistTimeCorrZDC    = new TH2F("fHistTimeCorrZDC", "ZDC;Corrected TDC timing C-A; Corrected TDC timing C+A", 120,-30,30,260,-100,30);
  fHistFMDA           = new TH1F("fHistFMDA", "FMDA;combinations above threshold;events", 102, -1.5, 100.5);
  fHistFMDC           = new TH1F("fHistFMDC", "FMDC;combinations above threshold;events", 102, -1.5, 100.5);
  fHistFMDSingle      = new TH1F("fHistFMDSingle", "FMD single;multiplicity value;counts", 1000, 0, 10);
  fHistFMDSum         = new TH1F("fHistFMDSum", "FMD sum;multiplicity value;counts", 1000, 0, 10);
  fHistT0             = new TH1F("fHistT0", "T0;time (ns);events", 100, -25, 25);
  
  TF1* fFuncV0MOnVsOf = new TF1("fFuncV0MOnVsOf","[0]+[1]*x",0,2000);
  fFuncV0MOnVsOf->SetParameters(fV0MOnVsOfA,fV0MOnVsOfB);
  fHistV0MOnVsOf->GetListOfFunctions()->Add(fFuncV0MOnVsOf);

  TF1* fFuncSPDOnVsOf = new TF1("fFuncSPDOnVsOf","[0]+[1]*x",0,800);
  fFuncSPDOnVsOf->SetParameters(fSPDOnVsOfA,fSPDOnVsOfB);
  fHistSPDOnVsOf->GetListOfFunctions()->Add(fFuncSPDOnVsOf);

  TF1* fFuncV0C3vs012 = new TF1("fFuncV0C3vs012","[0]+[1]*x",0,800);
  fFuncV0C3vs012->SetParameters(fV0CasymA,fV0CasymB);
  fHistV0C3vs012->GetListOfFunctions()->Add(fFuncV0C3vs012);
  
  fHistList = new TList();
  fHistList->Add(fHistV0A);
  fHistList->Add(fHistV0C);
  fHistList->Add(fHistFiredBitsSPD);
  fHistList->Add(fHistSPDClsVsTkl);
  fHistList->Add(fHistV0MOnVsOf);
  fHistList->Add(fHistSPDOnVsOf);
  fHistList->Add(fHistVIRvsBCmod4pup);
  fHistList->Add(fHistVIRvsBCmod4acc);
  fHistList->Add(fHistV0C3vs012);
  fHistList->Add(fHistBBAflags);
  fHistList->Add(fHistBBCflags);
  fHistList->Add(fHistBGAflags);
  fHistList->Add(fHistBGCflags);
  fHistList->Add(fHistV0MOn);
  fHistList->Add(fHistV0MOfAll);
  fHistList->Add(fHistV0MOfAcc);
  fHistList->Add(fHistSPDOnOuter);
  fHistList->Add(fHistAD);
  fHistList->Add(fHistADA);
  fHistList->Add(fHistADC);
  fHistList->Add(fHistV0A);
  fHistList->Add(fHistV0C);
  fHistList->Add(fHistZDC);
  fHistList->Add(fHistTDCZDC);
  fHistList->Add(fHistTimeZDC);
  fHistList->Add(fHistTimeCorrZDC);
  fHistList->Add(fHistFMDA);
  fHistList->Add(fHistFMDC);
  fHistList->Add(fHistFMDSingle);
  fHistList->Add(fHistFMDSum);
  fHistList->Add(fHistT0);
  fHistList->SetOwner();
  
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
    case kV0MOnVsOfPileup: str = "V0M on-cs-of pileup";       break;
    case kSPDOnVsOfPileup: str = "SPD on-cs-of pileup";       break;
    case kV0PFPileup:      str = "V0 PF pileup";              break;
    case kV0Casym:         str = "V0C012 vs V0C3 asymmetry";  break;
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
    case kVHM :            str = "VHM";                       break;
    case kV0M :            str = "V0M";                       break;
    case kSH1 :            str = "SH1";                       break;
    case kSH2 :            str = "SH2";                       break;
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
  return str;
}


//-------------------------------------------------------------------------------------------------
Bool_t AliTriggerAnalysis::IsTriggerFired(const AliVEvent* event, Trigger trigger){
  // checks if an event has been triggered
  if (trigger & kOfflineFlag) return IsOfflineTriggerFired(event, trigger);
  return IsTriggerBitFired(event, trigger);
}


//-------------------------------------------------------------------------------------------------
Bool_t AliTriggerAnalysis::IsTriggerBitFired(const AliVEvent* event, ULong64_t tclass) const {
  // Checks if corresponding bit in mask is on
  // TODO Why we need this function
  ULong64_t trigmask = event->GetTriggerMask();
  return (trigmask & (1ull << (tclass-1)));
}


//-------------------------------------------------------------------------------------------------
Int_t AliTriggerAnalysis::EvaluateTrigger(const AliVEvent* event, Trigger trigger){
  // evaluates a given trigger
  // trigger combinations are not supported, for that see IsOfflineTriggerFired

  UInt_t triggerNoFlags = (UInt_t) trigger % (UInt_t) kStartOfFlags;
  Bool_t offline = trigger & kOfflineFlag;
  
  if (!offline) {
    if ( triggerNoFlags==kT0BG
      || triggerNoFlags==kT0Pileup
      || triggerNoFlags==kSPDClsVsTrkBG
      || triggerNoFlags==kV0MOnVsOfPileup
      || triggerNoFlags==kSPDOnVsOfPileup
      || triggerNoFlags==kV0PFPileup
      || triggerNoFlags==kV0Casym
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
        ||triggerNoFlags==kVHM
        ||triggerNoFlags==kSH1
        ||triggerNoFlags==kSH2
        ||triggerNoFlags==kCentral
        ||triggerNoFlags==kSemiCentral
      ) AliFatal(Form("Offline trigger not available for trigger %d", triggerNoFlags));
  }
  
  switch (triggerNoFlags) {
    case kCTPV0A:          return event->GetHeader()->IsTriggerInputFired("V0A");
    case kCTPV0C:          return event->GetHeader()->IsTriggerInputFired("V0C");
    case kSPDGFO:          return SPDFiredChips(event, !offline, kFALSE, 0); 
    case kSPDGFOL0:        return SPDFiredChips(event, !offline, kFALSE, 1);
    case kSPDGFOL1:        return SPDFiredChips(event, !offline, kFALSE, 2);
    case kSPDClsVsTrkBG:   return IsSPDClusterVsTrackletBG(event);
    case kV0MOnVsOfPileup: return IsV0MOnVsOfPileup(event);
    case kSPDOnVsOfPileup: return IsSPDOnVsOfPileup(event);
    case kV0PFPileup:      return IsV0PFPileup(event);
    case kV0Casym:         return IsV0Casym(event);
    case kADA:             return ADTrigger(event, kASide, !offline) == kADBB; 
    case kADC:             return ADTrigger(event, kCSide, !offline) == kADBB;
    case kADABG:           return ADTrigger(event, kASide, !offline) == kADBG;
    case kADCBG:           return ADTrigger(event, kCSide, !offline) == kADBG;
    case kV0A:             return V0Trigger(event, kASide, !offline) == kV0BB; 
    case kV0C:             return V0Trigger(event, kCSide, !offline) == kV0BB;
    case kV0ABG:           return V0Trigger(event, kASide, !offline) == kV0BG;
    case kV0CBG:           return V0Trigger(event, kCSide, !offline) == kV0BG;
    case kVHM:             return VHMTrigger(event,!offline);
    case kV0M:             return V0MTrigger(event,!offline);
    case kSH1:             return SH1Trigger(event);
    case kSH2:             return SH2Trigger(event);
    case kT0:              return T0Trigger(event, !offline) == kT0BB;
    case kT0BG:            return T0Trigger(event, !offline) == kT0DecBG;
    case kT0Pileup:        return T0Trigger(event, !offline) == kT0DecPileup;
    case kZDCA:            return ZDCTrigger(event, kASide);
    case kZDCC:            return ZDCTrigger(event, kCSide);
    case kZDCTDCA:         return ZDCTDCTrigger(event, kASide);
    case kZDCTDCC:         return ZDCTDCTrigger(event, kCSide);
    case kZDCTime:         return ZDCTimeTrigger(event);
    case kZNA:             return ZDCTDCTrigger(event,kASide,kTRUE,kFALSE,kFALSE);
    case kZNC:             return ZDCTDCTrigger(event,kCSide,kTRUE,kFALSE,kFALSE);
    case kZNABG:           return ZDCTimeBGTrigger(event,kASide);
    case kZNCBG:           return ZDCTimeBGTrigger(event,kCSide);
    case kFMDA:            return FMDTrigger(event, kASide);
    case kFMDC:            return FMDTrigger(event, kCSide);
    case kTPCLaserWarmUp:  return IsLaserWarmUpTPCEvent(event);
    case kTPCHVdip:        return IsHVdipTPCEvent(event);
    case kIncompleteEvent: return IsIncompleteEvent(event);
    case kEMCAL:           return EMCALCellsTrigger(event);
    case kEmcalL0:         return EMCALTrigger(event,kEmcalL0);
    case kEmcalL1GammaHigh:return EMCALTrigger(event,kEmcalL1GammaHigh);
    case kEmcalL1GammaLow: return EMCALTrigger(event,kEmcalL1GammaLow);
    case kEmcalL1JetHigh:  return EMCALTrigger(event,kEmcalL1JetHigh);
    case kEmcalL1JetLow:   return EMCALTrigger(event,kEmcalL1JetLow);
    case kTRDHCO:          return TRDTrigger(event,kTRDHCO);
    case kTRDHJT:          return TRDTrigger(event,kTRDHJT);
    case kTRDHSE:          return TRDTrigger(event,kTRDHSE);
    case kTRDHQU:          return TRDTrigger(event,kTRDHQU);
    case kTRDHEE:          return TRDTrigger(event,kTRDHEE);
    case kCentral: {
      if (!event->GetVZEROData()) { AliWarning("V0 centrality trigger bits were not filled!"); return kFALSE; }
      if (!event->GetVZEROData()->TestBit(AliVVZERO::kTriggerChargeBitsFilled)) return kFALSE;
      return event->GetVZEROData()->GetTriggerBits() & (1<<AliVVZERO::kCTA2andCTC2);
    }
    case kSemiCentral: {
      if (!event->GetVZEROData()) { AliWarning("V0 centrality trigger bits were not filled!"); return kFALSE; }
      if (!event->GetVZEROData()->TestBit(AliVVZERO::kTriggerChargeBitsFilled)) return kFALSE;
      return event->GetVZEROData()->GetTriggerBits() & (1<<AliVVZERO::kCTA1andCTC1);
    }
    default: AliFatal(Form("Trigger type %d not implemented", triggerNoFlags));
  }
  return kFALSE;
}


//-------------------------------------------------------------------------------------------------
Bool_t AliTriggerAnalysis::IsOfflineTriggerFired(const AliVEvent* event, Trigger trigger){
  // checks if an event has been triggered "offline"
  UInt_t triggerNoFlags = (UInt_t) trigger % (UInt_t) kStartOfFlags;
  if (trigger & kOneParticle) AliError("AliTriggerAnalysis::kOneParticle functionality is obsolete");
  if (trigger & kOneTrack)    AliError("AliTriggerAnalysis::kOneTrack functionality is obsolete");

  Bool_t decision = kFALSE;
  switch (triggerNoFlags) {
    case kAcceptAll:        return kTRUE; 
    case kMB1:              return SPDGFOTrigger(event, 0) ||  V0Trigger(event, kASide, kFALSE) == kV0BB || V0Trigger(event, kCSide, kFALSE) == kV0BB;
    case kMB2:              return SPDGFOTrigger(event, 0) && (V0Trigger(event, kASide, kFALSE) == kV0BB || V0Trigger(event, kCSide, kFALSE) == kV0BB);
    case kMB3:              return SPDGFOTrigger(event, 0) &&  V0Trigger(event, kASide, kFALSE) == kV0BB && V0Trigger(event, kCSide, kFALSE) == kV0BB;
    case kSPDGFO:           return SPDGFOTrigger(event, 0);
    case kSPDGFOBits:       return SPDGFOTrigger(event, 1);
    case kV0MOnVsOfPileup:  return IsV0MOnVsOfPileup(event);
    case kSPDOnVsOfPileup:  return IsSPDOnVsOfPileup(event);
    case kV0PFPileup:       return IsV0PFPileup(event);
    case kV0Casym:          return IsV0Casym(event);
    case kADA:              return ADTrigger(event, kASide, kFALSE) == kADBB;
    case kADC:              return ADTrigger(event, kCSide, kFALSE) == kADBB;
    case kADABG:            return ADTrigger(event, kASide, kFALSE) == kADBG;
    case kADCBG:            return ADTrigger(event, kCSide, kFALSE) == kADBG;
    case kV0A:              return V0Trigger(event, kASide, kFALSE) == kV0BB;
    case kV0C:              return V0Trigger(event, kCSide, kFALSE) == kV0BB;
    case kV0OR:             return V0Trigger(event, kASide, kFALSE) == kV0BB || V0Trigger(event, kCSide, kFALSE) == kV0BB;
    case kV0AND:            return V0Trigger(event, kASide, kFALSE) == kV0BB && V0Trigger(event, kCSide, kFALSE) == kV0BB;
    case kV0ABG:            return V0Trigger(event, kASide, kFALSE) == kV0BG;
    case kV0CBG:            return V0Trigger(event, kCSide, kFALSE) == kV0BG;
    case kV0M:              return V0MTrigger(event,kFALSE);
    case kZDC:              return ZDCTrigger(event, kASide) || ZDCTrigger(event, kCentralBarrel) || ZDCTrigger(event, kCSide);
    case kZDCA:             return ZDCTrigger(event, kASide);
    case kZDCC:             return ZDCTrigger(event, kCSide);
    case kZNA:              return ZDCTDCTrigger(event,kASide,kTRUE,kFALSE,kFALSE);
    case kZNC:              return ZDCTDCTrigger(event,kCSide,kTRUE,kFALSE,kFALSE);
    case kZNABG:            return ZDCTimeBGTrigger(event,kASide);
    case kZNCBG:            return ZDCTimeBGTrigger(event,kCSide);
    case kFMDA:             return FMDTrigger(event, kASide);
    case kFMDC:             return FMDTrigger(event, kCSide);
    case kEMCAL:            return EMCALCellsTrigger(event);
    case kEmcalL0:          return EMCALTrigger(event,kEmcalL0);
    case kEmcalL1GammaHigh: return EMCALTrigger(event,kEmcalL1GammaHigh);
    case kEmcalL1GammaLow:  return EMCALTrigger(event,kEmcalL1GammaLow);
    case kEmcalL1JetHigh:   return EMCALTrigger(event,kEmcalL1JetHigh);
    case kEmcalL1JetLow:    return EMCALTrigger(event,kEmcalL1JetLow);
    case kTRDHCO:           return TRDTrigger(event,kTRDHCO);
    case kTRDHJT:           return TRDTrigger(event,kTRDHJT);
    case kTRDHSE:           return TRDTrigger(event,kTRDHSE);
    case kTRDHQU:           return TRDTrigger(event,kTRDHQU);
    case kTRDHEE:           return TRDTrigger(event,kTRDHEE);
    case kNSD1:             return SPDFiredChips(event, 0) >= 5 || (V0Trigger(event, kASide, kFALSE) == kV0BB && V0Trigger(event, kCSide, kFALSE) == kV0BB);
    case kFPANY:            decision |= SPDGFOTrigger(event, 0); 
                            decision |= V0Trigger(event, kASide, kFALSE) == kV0BB;
                            decision |= V0Trigger(event, kCSide, kFALSE) == kV0BB;
                            decision |= ZDCTrigger(event, kASide);
                            decision |= ZDCTrigger(event, kCentralBarrel);
                            decision |= ZDCTrigger(event, kCSide);
                            decision |= FMDTrigger(event, kASide);
                            decision |= FMDTrigger(event, kCSide);
                            return decision; 
    case kMB1Prime:         decision |= SPDGFOTrigger(event, 0) && V0Trigger(event, kASide, kFALSE) == kV0BB;
                            decision |= SPDGFOTrigger(event, 0) && V0Trigger(event, kCSide, kFALSE) == kV0BB;
                            decision |= V0Trigger(event, kASide, kFALSE) == kV0BB && V0Trigger(event, kCSide, kFALSE) == kV0BB;
                            return decision;
    default:                AliFatal(Form("Trigger type %d not implemented", triggerNoFlags));
  }
  
  return kFALSE;
}


//-------------------------------------------------------------------------------------------------
Int_t AliTriggerAnalysis::SPDFiredChips(const AliVEvent* event, Int_t origin, Bool_t fillHists, Int_t layer){
  // returns the number of fired chips in the SPD
  //
  // origin = 0 --> event->GetMultiplicity()->GetNumberOfFiredChips() (filled from clusters)
  // origin = 1 --> event->GetMultiplicity()->TestFastOrFiredChips() (from hardware bits)
  // layer  = 0 --> both layers
  // layer  = 1 --> inner
  // layer  = 2 --> outer
  
  const AliVMultiplicity* mult = event->GetMultiplicity();
  if (!mult) { 
    AliError("AliVMultiplicity not available"); 
    return -1; 
  }
  
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
AliTriggerAnalysis::ADDecision AliTriggerAnalysis::ADTrigger(const AliVEvent* event, AliceSide side, Bool_t online, Bool_t fillHists){
  // Returns the AD trigger decision 
  // argument 'online' is used as a switch between online and offline trigger algorithms
  
  const AliVAD* ad = event->GetADData();
  if (!ad) { 
    // print error only for runs from >=2015
    if (event->GetRunNumber()>=208505) AliError("AliVAD not available");
    return kADInvalid; 
  }
  if (side != kASide && side != kCSide) {
    AliError("Invalid AD side argument");
    return kADInvalid;
  }
  
  AliDebug(2,Form("In ADTrigger: %f %f",ad->GetADATime(),ad->GetADCTime()));

  if (online) {
    UShort_t bits = ad->GetTriggerBits();
    if (side==kASide && (bits & 1<<12)) return kADBB;
    if (side==kCSide && (bits & 1<<13)) return kADBB;
    if (side==kASide && (bits & 1<< 3)) return kADBG;
    if (side==kCSide && (bits & 1<< 5)) return kADBG;
  } else {
    if (fillHists) {
      if (fHistAD) fHistAD->Fill(ad->GetADCTime()-ad->GetADATime(),ad->GetADATime()+ad->GetADCTime());
      if (side == kASide && fHistADA) fHistADA->Fill(ad->GetADATime());
      if (side == kCSide && fHistADC) fHistADC->Fill(ad->GetADCTime());
    }
    if      (side == kASide) return (ADDecision) ad->GetADADecision();
    else if (side == kCSide) return (ADDecision) ad->GetADCDecision();
  }
  
  return kADEmpty;
}


//-------------------------------------------------------------------------------------------------
AliTriggerAnalysis::V0Decision AliTriggerAnalysis::V0Trigger(const AliVEvent* event, AliceSide side, Bool_t online, Bool_t fillHists){
  // Returns the V0 trigger decision 
  // argument 'online' is used as a switch between online and offline trigger algorithms
  
  const AliVVZERO* vzero = event->GetVZEROData();
  if (!vzero) { 
    AliError("AliVVZERO not available");  
    return kV0Invalid; 
  }
  if (!vzero->TestBit(AliVVZERO::kDecisionFilled)) {
    AliError("V0 decisions not filled");
    return kV0Invalid;
  }
  if (!vzero->TestBit(AliVVZERO::kOnlineBitsFilled)) {
    AliError("V0 online trigger bits not filled");
  }
  if (side != kASide && side != kCSide) {
    AliError("Invalid V0 side argument");
    return kV0Invalid;
  }
  
  AliDebug(2,Form("In V0Trigger: %f %f",vzero->GetV0ATime(),vzero->GetV0CTime()));
  
  if (online) {
    Int_t begin = (side == kASide) ? 32 :  0;
    Int_t end   = (side == kASide) ? 64 : 32;
    for (Int_t i=begin; i<end; i++) if (vzero->GetBBFlag(i)) return kV0BB;
    for (Int_t i=begin; i<end; i++) if (vzero->GetBGFlag(i)) return kV0BG;
  } else {
    if (fillHists) {
      if (side == kASide && fHistV0A) fHistV0A->Fill(vzero->GetV0ATime());
      if (side == kCSide && fHistV0C) fHistV0C->Fill(vzero->GetV0CTime());
    }
    if      (side == kASide) return (V0Decision) vzero->GetV0ADecision();
    else if (side == kCSide) return (V0Decision) vzero->GetV0CDecision();
  }
  
  return kV0Empty;
}


//-------------------------------------------------------------------------------------------------
Bool_t AliTriggerAnalysis::ZDCTDCTrigger(const AliVEvent* event, AliceSide side, Bool_t useZN, Bool_t useZP, Bool_t fillHists) const{
  // Returns if ZDC triggered, based on TDC information 
  if (event->GetDataLayoutType()!=AliVEvent::kESD) {
//    AliError("ZDCTDCTrigger method implemented for ESDs only");
    return kFALSE;
  }
  const AliESDEvent* aEsd = dynamic_cast<const AliESDEvent*>(event);

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
Bool_t AliTriggerAnalysis::ZDCTimeTrigger(const AliVEvent* event, Bool_t fillHists) const {
  // This method implements a selection based on the timing in both sides of zdcN
  // It can be used in order to eliminate parasitic collisions
  if (event->GetDataLayoutType()!=AliVEvent::kESD) {
//    AliError("ZDCTimeTrigger method implemented for ESDs only");
    return kFALSE;
  }
  const AliESDEvent* aEsd = dynamic_cast<const AliESDEvent*>(event);

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
    
    if (aEsd->GetRunNumber()>=245726 && aEsd->GetRunNumber()<=245793) detChZNA = 10; // use  timing from the common ZNA PMT
    
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
        if (esdZDC->TestBit(AliESDZDC::kCorrectedTDCFilled) && detChZNA == 10) {
          if (TMath::Power((tdcCcorr-tdcAcorr-123.1)/2.2,2.)+
              TMath::Power((tdcCcorr+tdcAcorr+123.1)/2.2,2.) < 1.) return kTRUE;
        } else if (esdZDC->TestBit(AliESDZDC::kCorrectedTDCFilled) && detChZNA!=10) {
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
Bool_t AliTriggerAnalysis::ZDCTimeBGTrigger(const AliVEvent* event, AliceSide side) const{
  // This method implements a selection based on the timing in zdcN
  // It can be used in order to flag background
  if(fMC) return kFALSE;
  
  if (event->GetDataLayoutType()!=AliVEvent::kESD) {
//    AliError("ZDCTimeBGTrigger method implemented for ESDs only");
    return kFALSE;
  }
  const AliESDEvent* aEsd = dynamic_cast<const AliESDEvent*>(event);

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
Bool_t AliTriggerAnalysis::ZDCTrigger(const AliVEvent* event, AliceSide side) const {
  // Returns if ZDC triggered

  if (event->GetDataLayoutType()!=AliVEvent::kESD) {
    AliError("ZDCTrigger method implemented for ESDs only");
    return kFALSE;
  }
  const AliESDEvent* aEsd = dynamic_cast<const AliESDEvent*>(event);
  
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
Bool_t AliTriggerAnalysis::FMDTrigger(const AliVEvent* event, AliceSide side){
  // Returns if the FMD triggered
  // Authors: FMD team, Hans Dalsgaard (code merged from FMD/AliFMDOfflineTrigger)
  if (event->GetDataLayoutType()!=AliVEvent::kESD) return 0;
  return FMDHitCombinations((AliESDEvent*) event, side, kFALSE);
}


//-------------------------------------------------------------------------------------------------
AliTriggerAnalysis::T0Decision AliTriggerAnalysis::T0Trigger(const AliVEvent* event, Bool_t online, Bool_t fillHists){
  // Returns the T0 TVDC trigger decision
  //  
  // argument 'online' is used as a switch between online and offline trigger algorithms
  // in online mode return 0TVX 
  // in offline mode in addition check pile-up and background :
  // pile-up read from ESD: check if TVDC (0TVX module name) has more 1 hit;
  // background flag read from ESD : check in given time interval OrA and OrC were correct but TVDC not
  // 
  // Based on an algorithm by Alla Maevskaya
  // TODO: implement online and offline selection in AOD
  // TODO: read vtx thresholds from OCDB
  
  if (event->GetDataLayoutType()==AliVEvent::kAOD) {
    // AOD analysis
    const AliAODTZERO* tzero = dynamic_cast<const AliAODEvent*>(event)->GetTZEROData();
    if (!tzero) {
      AliError("AliAODTZERO not available");
      return kT0Invalid;
    }
    if (fMC) if(tzero->GetT0zVertex()>-12.3 && tzero->GetT0zVertex() < 10.3) return kT0BB;
  } 
  else if (event->GetDataLayoutType()==AliVEvent::kESD) {
    // ESD analysis
    const AliESDTZERO* tzero = dynamic_cast<const AliESDEvent*>(event)->GetESDTZERO();
    if (!tzero) {
      AliError("AliESDTZERO not available");
      return kT0Invalid;
    }
    if (online) {
      if (event->GetHeader()->IsTriggerInputFired("0TVX")) return kT0BB;
    } else {
      Float_t tvdc0 = tzero->GetTVDC(0);
      if(fillHists) fHistT0->Fill(tvdc0);
      if (tzero->GetPileupFlag()) return kT0DecPileup;
      if (tzero->GetBackgroundFlag()) return kT0DecBG;
      if (tvdc0>-5 && tvdc0<5 && tvdc0!=0) return kT0BB;
    }
    if (fMC) if (tzero->GetT0zVertex()>-12.3 && tzero->GetT0zVertex() < 10.3) return kT0BB;
  }
  
  return kT0Empty;
}

//-------------------------------------------------------------------------------------------------
Bool_t AliTriggerAnalysis::EMCALCellsTrigger(const AliVEvent* event){
  //
  // Returns the EMCAL trigger decision
  // so far only implemented for LHC11a data
  // see http://alisoft.cern.ch/viewvc/trunk/PWGGA/EMCALTasks/AliEmcalPhysicsSelection.cxx?view=markup&root=AliRoot Revision 56136
  //
  
  Bool_t isFired = kTRUE;
  const Int_t runNumber = event->GetRunNumber();
  
  // Get EMCAL cells
  AliVCaloCells *cells = event->GetEMCALCells();
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
Bool_t AliTriggerAnalysis::TRDTrigger(const AliVEvent* event, Trigger trigger){
  // evaluate the TRD trigger conditions,
  // so far HCO, HSE, HQU, HJT, HEE
  if(trigger!=kTRDHCO && trigger!=kTRDHJT && trigger!=kTRDHSE && trigger!=kTRDHQU && trigger!=kTRDHEE) {
    AliWarning("Beware you are erroneously trying to use this function (wrong trigger)");
    return kFALSE;
  }
  
  Int_t nTrdTracks = event->GetNumberOfTrdTracks();
  if (nTrdTracks<=0) return kFALSE;
  if      (trigger==kTRDHCO) return kTRUE;
  else if (trigger!=kTRDHJT) {
    for (Int_t iTrack = 0; iTrack < nTrdTracks; ++iTrack) {
      AliVTrdTrack* trdTrack = event->GetTrdTrack(iTrack);
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
      AliVTrdTrack *trdTrack = event->GetTrdTrack(iTrack);    
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
  if (!mult) { 
    AliError("No multiplicity object"); 
    return kFALSE; 
  }
  Int_t nTkl = mult->GetNumberOfTracklets();
  Int_t nCls = event->GetNumberOfITSClusters(0) + event->GetNumberOfITSClusters(1);
  if (fillHists) fHistSPDClsVsTkl->Fill(nTkl,nCls);
  return nCls > fSPDClsVsTklA + nTkl*fSPDClsVsTklB;
}


//-------------------------------------------------------------------------------------------------
Bool_t AliTriggerAnalysis::IsLaserWarmUpTPCEvent(const AliVEvent* event){
  // This function flags noisy TPC events which can happen during laser warm-up.
  Int_t trackCounter = 0;
  for (Int_t i=0; i<event->GetNumberOfTracks(); i++) {
    AliVTrack *track = dynamic_cast<AliVTrack*>(event->GetTrack(i));
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
Bool_t AliTriggerAnalysis::IsHVdipTPCEvent(const AliVEvent* event) {
  // This function flags events in which the TPC chamber HV is not at its nominal value
  if (fMC) return kFALSE; // there are no dip events in MC
  
  if (event->GetDataLayoutType()!=AliVEvent::kESD) {
    AliError("IsHVdipTPCEvent method implemented for ESDs only");
    return kFALSE;
  }
  const AliESDEvent* aEsd = dynamic_cast<const AliESDEvent*>(event);

  if (!aEsd->IsDetectorOn(AliDAQ::kTPC)) return kTRUE;
  return kFALSE;
}


//-------------------------------------------------------------------------------------------------
Bool_t AliTriggerAnalysis::IsIncompleteEvent(const AliVEvent* event){
  // Check whether the event is incomplete 
  // (due to DAQ-HLT issues, it could be only part of the event was saved)
  if (fMC) return kFALSE; // there are no incomplete events on MC
  return const_cast<AliVEvent*>(event)->IsIncompleteDAQ(); 
}

//-------------------------------------------------------------------------------------------------
Bool_t AliTriggerAnalysis::IsV0MOnVsOfPileup(const AliVEvent* event, Bool_t fillHists){
  if (fMC) return kFALSE;
  AliVVZERO* vzero = event->GetVZEROData();
  if (!vzero) {
    AliError("AliVVZERO not available");
    return kFALSE;
  }
  Float_t on = vzero->GetTriggerChargeA()+vzero->GetTriggerChargeC();
  Float_t of = vzero->GetMTotV0A()+vzero->GetMTotV0C();
  if (fillHists) fHistV0MOnVsOf->Fill(of,on);
  return (on < fV0MOnVsOfA + fV0MOnVsOfB*of);
}

//-------------------------------------------------------------------------------------------------
Bool_t AliTriggerAnalysis::IsSPDOnVsOfPileup(const AliVEvent* event, Bool_t fillHists){
  if (fMC) return kFALSE;
  AliVMultiplicity* mult = event->GetMultiplicity();
  if (!mult) {
    AliError("AliVMultiplicity not available");
    return kFALSE;
  }
  TBits onMap = mult->GetFastOrFiredChips();
  TBits ofMap = mult->GetFiredChipMap();
  Int_t on = onMap.CountBits(400);
  Int_t of = ofMap.CountBits(400);
  if (fillHists) fHistSPDOnVsOf->Fill(of,on);
  return (on < fSPDOnVsOfA + fSPDOnVsOfB*of);
}


//-------------------------------------------------------------------------------------------------
Bool_t AliTriggerAnalysis::IsV0PFPileup(const AliVEvent* event, Bool_t fillHists){
  if (fMC) return kFALSE;
  AliVVZERO* vzero = event->GetVZEROData();
  if (!vzero) {
    AliError("AliVVZERO not available");
    return kFALSE;
  }

  Bool_t vir[21] = {0};
  UChar_t bcMod4 = event->GetBunchCrossNumber()%4;

  for (Int_t bc=0;bc<=20;bc++) {
    UChar_t nBBA=0;
    UChar_t nBBC=0;
    UChar_t nBGA=0;
    UChar_t nBGC=0;
    if (fVIRBBAflags<33) for (Int_t i=0;i<32;i++) nBBA+=vzero->GetPFBBFlag(i+32,bc);
    if (fVIRBBCflags<33) for (Int_t i=0;i<32;i++) nBBC+=vzero->GetPFBBFlag(i   ,bc);
    if (fVIRBGAflags<33) for (Int_t i=0;i<32;i++) nBGA+=vzero->GetPFBGFlag(i+32,bc);
    if (fVIRBGCflags<33) for (Int_t i=0;i<32;i++) nBGC+=vzero->GetPFBGFlag(i   ,bc);
    vir[bc] |= nBBA>=fVIRBBAflags;
    vir[bc] |= nBBC>=fVIRBBCflags;
    vir[bc] |= nBGA>=fVIRBGAflags;
    vir[bc] |= nBGC>=fVIRBGCflags;
    if (fillHists) {
      if (bc==10) continue;
      if (!vir[bc]) continue;
      if (!IsSPDOnVsOfPileup(event) && !IsV0MOnVsOfPileup(event)) continue;
      fHistVIRvsBCmod4pup->Fill(10-bc,bcMod4);
    }
  }
  
  // clock index is counting from future to past
  Int_t bcMin = 10 - fNBCsFuture + bcMod4;
  Int_t bcMax = 10 + fNBCsPast   + bcMod4;
  for (Int_t bc=bcMin;bc<=bcMax;bc++) {
    if (bc==10) continue; // skip current bc
    if (vir[bc]) return kTRUE;
  }

  if (fillHists) {
    for (Int_t bc=0;bc<=20;bc++) {
      if (bc==10) continue;
      if (!vir[bc]) continue;
      fHistVIRvsBCmod4acc->Fill(10-bc,bcMod4);
    }
  }
  return kFALSE;
}


//-------------------------------------------------------------------------------------------------
Bool_t AliTriggerAnalysis::IsV0Casym(const AliVEvent* event, Bool_t fillHists){
  if (fMC) return kFALSE;
  AliVVZERO* vzero = event->GetVZEROData();
  if (!vzero) {
    AliError("AliVVZERO not available");
    return kFALSE;
  }
  Float_t multV0C012 = vzero->GetMRingV0C(0)+vzero->GetMRingV0C(1)+vzero->GetMRingV0C(2);
  Float_t multV0C3   = vzero->GetMRingV0C(3);
  
  if (fillHists) fHistV0C3vs012->Fill(multV0C012,multV0C3);
  return (multV0C3 < fV0CasymA + fV0CasymB*multV0C012);
}

//-------------------------------------------------------------------------------------------------
Bool_t AliTriggerAnalysis::VHMTrigger(const AliVEvent* event, Bool_t fillHists){
  AliVVZERO* vzero = event->GetVZEROData();
  if (!vzero) {
    AliError("AliVVZERO not available");
    return kFALSE;
  }
  Int_t nBBA = 0;
  Int_t nBBC = 0;
  Int_t nBGA = 0;
  Int_t nBGC = 0;
  
  for (Int_t i=0;i<32;i++) {
    nBBA += vzero->GetBBFlag(i+32);
    nBBC += vzero->GetBBFlag(i   );
    nBGA += vzero->GetBGFlag(i+32);
    nBGC += vzero->GetBGFlag(i   );
  }
  if (fillHists){
    fHistBBAflags->Fill(nBBA);
    fHistBBCflags->Fill(nBBC);
    fHistBGAflags->Fill(nBGA);
    fHistBGCflags->Fill(nBGC);
  }
  
  Bool_t vhm = 1;
  vhm *= nBBA>=fVHMBBAflags;
  vhm *= nBBC>=fVHMBBCflags;
  vhm *= nBGA<=fVHMBGAflags;
  vhm *= nBGC<=fVHMBGCflags;
  
  return vhm;
}


//-------------------------------------------------------------------------------------------------
Bool_t AliTriggerAnalysis::V0MTrigger(const AliVEvent* event, Bool_t online, Bool_t fillHists){
  AliVVZERO* vzero = event->GetVZEROData();
  if (!vzero) {
    AliError("AliVVZERO not available");
    return kFALSE;
  }
  Int_t   on = vzero->GetTriggerChargeA()+vzero->GetTriggerChargeC();
  Float_t of = vzero->GetMTotV0A()+vzero->GetMTotV0C();

  if (fillHists) {
    fHistV0MOn->Fill(on);
    fHistV0MOfAll->Fill(of);
    if (of>=fV0MOfThreshold) fHistV0MOfAcc->Fill(of);
  }
  
  return online ? on>=fV0MOnThreshold: of>=fV0MOfThreshold;
}


//-------------------------------------------------------------------------------------------------
Bool_t AliTriggerAnalysis::SH1Trigger(const AliVEvent* event, Bool_t fillHists){
  if (fMC) return kFALSE;
  AliVMultiplicity* mult = event->GetMultiplicity();
  if (!mult) {
    AliError("AliVMultiplicity not available");
    return kFALSE;
  }
  TBits onMap = mult->GetFastOrFiredChips();
  Int_t on = onMap.CountBits(400);
  if (fillHists) fHistSPDOnOuter->Fill(on);
  return (on>=fSH1OuterThreshold);
}


//-------------------------------------------------------------------------------------------------
Bool_t AliTriggerAnalysis::SH2Trigger(const AliVEvent* event, Bool_t fillHists){
  if (fMC) return kFALSE;
  AliVMultiplicity* mult = event->GetMultiplicity();
  if (!mult) {
    AliError("AliVMultiplicity not available");
    return kFALSE;
  }
  TBits onMap = mult->GetFastOrFiredChips();
  Int_t on = onMap.CountBits(400);
  if (fillHists) fHistSPDOnOuter->Fill(on);
  return (on>=fSH2OuterThreshold);
}


//-------------------------------------------------------------------------------------------------
Long64_t AliTriggerAnalysis::Merge(TCollection* list){
  // Merge a list of objects with this (needed for PROOF).
  // Returns the number of merged objects (including this).
  if (!list) return 0;
  if (list->IsEmpty()) return 1;
  TIterator* iter = list->MakeIterator();
  TObject* obj;
  TList histListCollection; 
  Int_t count = 0;
  while ((obj = iter->Next())) {
    AliTriggerAnalysis* entry = dynamic_cast<AliTriggerAnalysis*> (obj);
    if (entry == 0) continue;
    histListCollection.Add(entry->fHistList);
    
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
  fHistList->Merge(&histListCollection);
  delete iter;
  return count+1;
}


//-------------------------------------------------------------------------------------------------
void AliTriggerAnalysis::FillHistograms(const AliVEvent* event){
  SPDFiredChips(event,1,kTRUE,0);
  Bool_t decisionADA       = ADTrigger(event, kASide, kFALSE, kTRUE);
  Bool_t decisionADC       = ADTrigger(event, kCSide, kFALSE, kTRUE);
  Bool_t decisionV0A       = V0Trigger(event, kASide, kFALSE, kTRUE);
  Bool_t decisionV0C       = V0Trigger(event, kCSide, kFALSE, kTRUE);
  Bool_t isZDCTDCTrigger   = ZDCTDCTrigger(event,kASide,kFALSE,kFALSE,kTRUE);
  Bool_t isSPDClsVsTklBG   = IsSPDClusterVsTrackletBG(event, kTRUE);
  Bool_t isV0MOnVsOfPileup = IsV0MOnVsOfPileup(event, kTRUE);
  Bool_t isSPDOnVsOfPileup = IsSPDOnVsOfPileup(event, kTRUE);
  Bool_t isV0PFPileup      = IsV0PFPileup(event,kTRUE);
  Bool_t isV0Casym         = IsV0Casym(event,kTRUE);
  Bool_t isVHMTrigger      = VHMTrigger(event,kTRUE);
  Bool_t isV0MTrigger      = V0MTrigger(event,kFALSE,kTRUE);
  Bool_t isSH1Trigger      = SH1Trigger(event,kTRUE);
//  TODO: Adjust for AOD
//  AliESDZDC* zdcData = event->GetESDZDC();
//  if (zdcData)  {
//    UInt_t quality = zdcData->GetESDQuality();
//    
//    // from Nora's presentation, general first physics meeting 16.10.09
//    static UInt_t zpc  = 0x20;
//    static UInt_t znc  = 0x10;
//    static UInt_t zem1 = 0x08;
//    static UInt_t zem2 = 0x04;
//    static UInt_t zpa  = 0x02;
//    static UInt_t zna  = 0x01;
//    
//    fHistZDC->Fill(1, (quality & zna)  ? 1 : 0);
//    fHistZDC->Fill(2, (quality & zpa)  ? 1 : 0);
//    fHistZDC->Fill(3, (quality & zem2) ? 1 : 0);
//    fHistZDC->Fill(4, (quality & zem1) ? 1 : 0);
//    fHistZDC->Fill(5, (quality & znc)  ? 1 : 0);
//    fHistZDC->Fill(6, (quality & zpc)  ? 1 : 0);
//  }
//  else {
//    fHistZDC->Fill(-1);
//    AliError("AliESDZDC not available");
//  }
  
//  if (fDoFMD) {
//    fHistFMDA->Fill(FMDHitCombinations(event, kASide, kTRUE));
//    fHistFMDC->Fill(FMDHitCombinations(event, kCSide, kTRUE));
//  }
}


//-------------------------------------------------------------------------------------------------
void AliTriggerAnalysis::SaveHistograms() const {
  // write histograms to current directory
  if (fSPDGFOEfficiency)   fSPDGFOEfficiency->Write();
  fTriggerClasses->Write("fTriggerClasses", TObject::kSingleKey);
  fHistList->Write("histos",TObject::kSingleKey);
}


//-------------------------------------------------------------------------------------------------
void AliTriggerAnalysis::FillTriggerClasses(const AliVEvent* event){
  // fills trigger classes map
  TParameter<Long64_t>* count = dynamic_cast<TParameter<Long64_t>*> (fTriggerClasses->GetValue(event->GetFiredTriggerClasses().Data()));
  if (!count) {
    count = new TParameter<Long64_t>(event->GetFiredTriggerClasses(), 0);
    fTriggerClasses->Add(new TObjString(event->GetFiredTriggerClasses().Data()), count);
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

void AliTriggerAnalysis::Browse(TBrowser *b){
   // Browse this object.
   // If b=0, there is no Browse call TObject::Browse(0) instead.
   //         This means TObject::Inspect() will be invoked indirectly


  static TObjString * strZDCCutRefSumCorr     =0;    
  static TObjString * strZDCCutRefDeltaCorr   =0;  
  static TObjString * strZDCCutSigmaSumCorr   =0;  
  static TObjString * strZDCCutSigmaDeltaCorr =0;
  static TObjString * strZDCCutZNATimeCorrMin =0;
  static TObjString * strZDCCutZNATimeCorrMax =0;
  static TObjString * strZDCCutZNCTimeCorrMin =0;
  static TObjString * strZDCCutZNCTimeCorrMax =0;

  if(strZDCCutRefSumCorr     ) delete strZDCCutRefSumCorr     ;
  if(strZDCCutRefDeltaCorr   ) delete strZDCCutRefDeltaCorr   ;
  if(strZDCCutSigmaSumCorr   ) delete strZDCCutSigmaSumCorr   ;
  if(strZDCCutSigmaDeltaCorr ) delete strZDCCutSigmaDeltaCorr ;
  if(strZDCCutZNATimeCorrMin ) delete strZDCCutZNATimeCorrMin ;
  if(strZDCCutZNATimeCorrMax ) delete strZDCCutZNATimeCorrMax ;
  if(strZDCCutZNCTimeCorrMin ) delete strZDCCutZNCTimeCorrMin ;
  if(strZDCCutZNCTimeCorrMax ) delete strZDCCutZNCTimeCorrMax ;
  
  strZDCCutRefSumCorr     = new TObjString(Form("ZDCCutRefSumCorr     %f", fZDCCutRefSumCorr    )); 
  strZDCCutRefDeltaCorr   = new TObjString(Form("ZDCCutRefDeltaCorr   %f", fZDCCutRefDeltaCorr  )); 
  strZDCCutSigmaSumCorr   = new TObjString(Form("ZDCCutSigmaSumCorr   %f", fZDCCutSigmaSumCorr  )); 
  strZDCCutSigmaDeltaCorr = new TObjString(Form("ZDCCutSigmaDeltaCorr %f", fZDCCutSigmaDeltaCorr)); 
  strZDCCutZNATimeCorrMin = new TObjString(Form("ZDCCutZNATimeCorrMin %f", fZDCCutZNATimeCorrMin));
  strZDCCutZNATimeCorrMax = new TObjString(Form("ZDCCutZNATimeCorrMax %f", fZDCCutZNATimeCorrMax));
  strZDCCutZNCTimeCorrMin = new TObjString(Form("ZDCCutZNCTimeCorrMin %f", fZDCCutZNCTimeCorrMin));
  strZDCCutZNCTimeCorrMax = new TObjString(Form("ZDCCutZNCTimeCorrMax %f", fZDCCutZNCTimeCorrMax));

  if (b) {
    // Creates a folder for each beam type containing the list of corresponding bx ids
    b->Add(strZDCCutRefSumCorr    );
    b->Add(strZDCCutRefDeltaCorr  );
    b->Add(strZDCCutSigmaSumCorr  );
    b->Add(strZDCCutSigmaDeltaCorr);
    b->Add(strZDCCutZNATimeCorrMin);
    b->Add(strZDCCutZNATimeCorrMax);
    b->Add(strZDCCutZNCTimeCorrMin);
    b->Add(strZDCCutZNCTimeCorrMax);
    b->Add(fHistList);
  }     
  else
    TObject::Browse(b);
}
