#ifndef ALITRDCALDCSPTRCBB_H
#define ALITRDCALDCSPTRCBB_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliTRDCalDCSPTRCbb.h 18952 2007-06-08 11:36:12Z cblume $ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  TRD calibration class for TRD GTU configuration parameters               //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "TNamed.h"

class TString;

class AliTRDCalDCSPTRCbb : public TNamed {

 public:

  AliTRDCalDCSPTRCbb();
  AliTRDCalDCSPTRCbb(const char *name, const char *title);
  AliTRDCalDCSPTRCbb(const AliTRDCalDCSPTRCbb &);
  virtual ~AliTRDCalDCSPTRCbb() { };

  TString GetControlBoxSide()                         { return fSide;                         }
  Int_t   GetControlBoxPrimary()                      { return fPrimary;                      }
  
  void    SetControlBoxSide(TString bs) const        { fSide = bs;                           }
  void    SetControlBoxPrimary(Int_t bp) const        { fPrimary = bp;                        }
  UInt_t  GetPreToNextCycles() const                  { return fPreToNextCycles;                         }
  UInt_t  GetL0ToNextCycles() const                   { return fL0ToNextCycles;                         }
  UInt_t  GetL1ToNextCycles() const                   { return fL1ToNextCycles;                         }
  UInt_t  GetClkTtcexShifts() const                   { return fClkTtcexShifts;                         }
  UInt_t  GetDisableGtuBusy() const                   { return fDisableGtuBusy;                         }
  UInt_t  GetDisableSorEorBusy() const                { return fDisableSorEorBusy;                         }
  UInt_t  GetScintEn() const                          { return fScintEn;                         }
  UInt_t  GetLutAsPre() const                         { return fLutAsPre;                         }
  UInt_t  GetCtpAsPre() const                         { return fCtpAsPre;                         }
  UInt_t  GetL0En() const                             { return fL0En;                         }
  UInt_t  GetL1En() const                             { return fL1En;                         }
  UInt_t  GetLutL0ToCtpEn() const                     { return fLutL0ToCtpEn;                         }
  UInt_t  GetLutPreToCtpEn() const                    { return fLutPreToCtpEn;                         }
  UInt_t  GetSoftwTrigToCtpEn() const                 { return fSoftwTrigToCtpEn;                         }
  UInt_t  GetL0Autogenerate() const                   { return fL0Autogenerate;                         }
  UInt_t  GetL1Autogenerate() const                   { return fL1Autogenerate;                         }
  UInt_t  GetChannelBDisable() const                  { return fChannelBDisable;                         }
  UInt_t  GetTtcexClkDisable() const                  { return fTtcexClkDisable;                         }
  UInt_t  GetPreAsL0En() const                        { return fPreAsL0En;                         }
  UInt_t  GetNoTriggerToSm() const                    { return fNoTriggerToSm;                         }
  UInt_t  GetSoftwTrigEn() const                      { return fSoftwTrigEn;                         }
  UInt_t  GetSoftwTrigPattern() const                 { return fSoftwTrigPattern;                         }
  UInt_t  GetSychrInput() const                       { return fSychrInput;                         }
  UInt_t  GetRandomTriggerThr() const                 { return fRandomTriggerThr;                         }
  UInt_t  GetPtChDelayCbA() const                     { return fPtChDelayCbA;                         }
  UInt_t  GetPtChDelayCbC() const                     { return fPtChDelayCbC;                         }
  UInt_t  GetPtChDelayGtu() const                     { return fPtChDelayGtu;                         }
  UInt_t  GetPtChDelayT0() const                      { return fPtChDelayT0;                         }
  UInt_t  GetPtChDelayT0() const                      { return fPtChDelayT0;                         }
  UInt_t  GetPtChDelayT1() const                      { return fPtChDelayT1;                         }
  UInt_t  GetPtChDelayT2() const                      { return fPtChDelayT2;                         }
  UInt_t  GetPtChDelayT3() const                      { return fPtChDelayT3;                         }
  UInt_t  GetPtChDelayT4() const                      { return fPtChDelayT4;                         }
  UInt_t  GetPtChDelayT5() const                      { return fPtChDelayT5;                         }
  UInt_t  GetPtChDelayT6() const                      { return fPtChDelayT6;                         }
  UInt_t  GetPtChDelayT7() const                      { return fPtChDelayT7;                         }
  UInt_t  GetClkLb() const                            { return fClkLb;                         }
  UInt_t  GetClkHb() const                            { return fClkHb;                         }
  UInt_t  GetPulseToSmLb() const                      { return fPulseToSmLb;                         }
  UInt_t  GetPulseToSmHb() const                      { return fPulseToSmHb;                         }
  UInt_t  GetDeadTimeLb() const                       { return fDeadTimeLb;                         }
  UInt_t  GetDeadTimeHb() const                       { return fDeadTimeHb;                         }
  UInt_t  GetPrePulsesLb() const                      { return fPrePulsesLb;                         }
  UInt_t  GetPrePulsesHb() const                      { return fPrePulsesHb;                         }
  UInt_t  GetL0CtpLb() const                          { return fL0CtpLb;                         }
  UInt_t  GetL0CtpHb() const                          { return fL0CtpHb;                         }
  UInt_t  GetL1CtpLb() const                          { return fL1CtpLb;                         }
  UInt_t  GetL1CtpHb() const                          { return fL1CtpHb;                         }
  UInt_t  GetL0LutLb() const                          { return fL0LutLb;                         }
  UInt_t  GetL0LutHb() const                          { return fL0LutHb;                         }
  UInt_t  GetNotTinToCtpPreLb() const                 { return fNotTinToCtpPreLb;                         }
  UInt_t  GetNotTinToCtpPreHb() const                 { return fNotTinToCtpPreHb;                         }
  UInt_t  GetPatternMatchCbALb() const                { return fPatternMatchCbALb;                         }
  UInt_t  GetPatternMatchCbAHb() const                { return fPatternMatchCbAHb;                         }
  UInt_t  GetPatternMatchCbCLb() const                { return fPatternMatchCbCLb;                         }
  UInt_t  GetPatternMatchCbCHb() const                { return fPatternMatchCbCHb;                         }
  UInt_t  GetPatternMatchCbTlmuLb() const             { return fPatternMatchCbTlmuLb;                         }
  UInt_t  GetPatternMatchCbTlmuHb() const             { return fPatternMatchCbTlmuHb;                         }
  UInt_t  GetTriggerSynchrLb() const                  { return fTriggerSynchrLb;                         }
  UInt_t  GetTriggerSynchrHb() const                  { return fTriggerSynchrHb;                         }
  UInt_t  GetTriggerToSmLb() const                    { return fTriggerToSmLb;                         }
  UInt_t  GetTriggerToSmHb() const                    { return fTriggerToSmHb;                         }
  UInt_t  GetPreLutLb() const                         { return fPreLutLb;                         }
  UInt_t  GetPreLutHb() const                         { return fPreLutHb;                         }
  UInt_t  GetMissingPreLb() const                     { return fMissingPreLb;                         }
  UInt_t  GetMissingPreHb() const                     { return fMissingPreHb;                         }
  UInt_t  GetUnnecessaryPreLb() const                 { return fUnnecessaryPreLb;                         }
  UInt_t  GetUnnecessaryPreHb() const                 { return fUnnecessaryPreHb;                         }
  
  
  void    SetPreToNextCycles(UInt_t ar)               { fPreToNextCycles = ar;                           }
  void    SetL0ToNextCycles(UInt_t ar)                { fL0ToNextCycles = ar;                           }
  void    SetL1ToNextCycles(UInt_t ar)                { fL1ToNextCycles = ar;                           }
  void    SetClkTtcexShifts(UInt_t ar)                { fClkTtcexShifts = ar;                           }
  void    SetDisableGtuBusy(UInt_t ar)                { fDisableGtuBusy = ar;                           }
  void    SetDisableSorEorBusy(UInt_t ar)             { fDisableSorEorBusy = ar;                           }
  void    SetScintEn(UInt_t ar)                       { fScintEn = ar;                           }
  void    SetLutAsPre(UInt_t ar)                      { fLutAsPre = ar;                           }
  void    SetCtpAsPre(UInt_t ar)                      { fCtpAsPre = ar;                           }
  void    SetL0En(UInt_t ar)                          { fL0En = ar;                           }
  void    SetL1En(UInt_t ar)                          { fL1En = ar;                           }
  void    SetLutL0ToCtpEn(UInt_t ar)                  { fLutL0ToCtpEn = ar;                           }
  void    SetLutPreToCtpEn(UInt_t ar)                 { fLutPreToCtpEn = ar;                           }
  void    SetSoftwTrigToCtpEn(UInt_t ar)              { fSoftwTrigToCtpEn = ar;                           }
  void    SetL0Autogenerate(UInt_t ar)                { fL0Autogenerate = ar;                           }
  void    SetL1Autogenerate(UInt_t ar)                { fL1Autogenerate = ar;                           }
  void    SetChannelBDisable(UInt_t ar)               { fChannelBDisable = ar;                           }
  void    SetTtcexClkDisable(UInt_t ar)               { fTtcexClkDisable = ar;                           }
  void    SetPreAsL0En(UInt_t ar)                     { fPreAsL0En = ar;                           }
  void    SetNoTriggerToSm(UInt_t ar)                 { fNoTriggerToSm = ar;                           }
  void    SetSoftwTrigEn(UInt_t ar)                   { fSoftwTrigEn = ar;                           }
  void    SetSoftwTrigPattern(UInt_t ar)              { fSoftwTrigPattern = ar;                           }
  void    SetSychrInput(UInt_t ar)                    { fSychrInput = ar;                           }
  void    SetRandomTriggerThr(UInt_t ar)              { fRandomTriggerThr = ar;                           }
  void    SetPtChDelayCbA(UInt_t ar)                  { fPtChDelayCbA = ar;                           }
  void    SetPtChDelayCbC(UInt_t ar)                  { fPtChDelayCbC = ar;                           }
  void    SetPtChDelayGtu(UInt_t ar)                  { fPtChDelayGtu = ar;                           }
  void    SetPtChDelayT0(UInt_t ar)                   { fPtChDelayT0 = ar;                           }
  void    SetPtChDelayT0(UInt_t ar)                   { fPtChDelayT0 = ar;                           }
  void    SetPtChDelayT1(UInt_t ar)                   { fPtChDelayT1 = ar;                           }
  void    SetPtChDelayT2(UInt_t ar)                   { fPtChDelayT2 = ar;                           }
  void    SetPtChDelayT3(UInt_t ar)                   { fPtChDelayT3 = ar;                           }
  void    SetPtChDelayT4(UInt_t ar)                   { fPtChDelayT4 = ar;                           }
  void    SetPtChDelayT5(UInt_t ar)                   { fPtChDelayT5 = ar;                           }
  void    SetPtChDelayT6(UInt_t ar)                   { fPtChDelayT6 = ar;                           }
  void    SetPtChDelayT7(UInt_t ar)                   { fPtChDelayT7 = ar;                           }
  void    SetClkLb(UInt_t ar)                         { fClkLb = ar;                           }
  void    SetClkHb(UInt_t ar)                         { fClkHb = ar;                           }
  void    SetPulseToSmLb(UInt_t ar)                   { fPulseToSmLb = ar;                           }
  void    SetPulseToSmHb(UInt_t ar)                   { fPulseToSmHb = ar;                           }
  void    SetDeadTimeLb(UInt_t ar)                    { fDeadTimeLb = ar;                           }
  void    SetDeadTimeHb(UInt_t ar)                    { fDeadTimeHb = ar;                           }
  void    SetPrePulsesLb(UInt_t ar)                   { fPrePulsesLb = ar;                           }
  void    SetPrePulsesHb(UInt_t ar)                   { fPrePulsesHb = ar;                           }
  void    SetL0CtpLb(UInt_t ar)                       { fL0CtpLb = ar;                           }
  void    SetL0CtpHb(UInt_t ar)                       { fL0CtpHb = ar;                           }
  void    SetL1CtpLb(UInt_t ar)                       { fL1CtpLb = ar;                           }
  void    SetL1CtpHb(UInt_t ar)                       { fL1CtpHb = ar;                           }
  void    SetL0LutLb(UInt_t ar)                       { fL0LutLb = ar;                           }
  void    SetL0LutHb(UInt_t ar)                       { fL0LutHb = ar;                           }
  void    SetNotTinToCtpPreLb(UInt_t ar)              { fNotTinToCtpPreLb = ar;                           }
  void    SetNotTinToCtpPreHb(UInt_t ar)              { fNotTinToCtpPreHb = ar;                           }
  void    SetPatternMatchCbALb(UInt_t ar)             { fPatternMatchCbALb = ar;                           }
  void    SetPatternMatchCbAHb(UInt_t ar)             { fPatternMatchCbAHb = ar;                           }
  void    SetPatternMatchCbCLb(UInt_t ar)             { fPatternMatchCbCLb = ar;                           }
  void    SetPatternMatchCbCHb(UInt_t ar)             { fPatternMatchCbCHb = ar;                           }
  void    SetPatternMatchCbTlmuLb(UInt_t ar)          { fPatternMatchCbTlmuLb = ar;                           }
  void    SetPatternMatchCbTlmuHb(UInt_t ar)          { fPatternMatchCbTlmuHb = ar;                           }
  void    SetTriggerSynchrLb(UInt_t ar)               { fTriggerSynchrLb = ar;                           }
  void    SetTriggerSynchrHb(UInt_t ar)               { fTriggerSynchrHb = ar;                           }
  void    SetTriggerToSmLb(UInt_t ar)                 { fTriggerToSmLb = ar;                           }
  void    SetTriggerToSmHb(UInt_t ar)                 { fTriggerToSmHb = ar;                           }
  void    SetPreLutLb(UInt_t ar)                      { fPreLutLb = ar;                           }
  void    SetPreLutHb(UInt_t ar)                      { fPreLutHb = ar;                           }
  void    SetMissingPreLb(UInt_t ar)                  { fMissingPreLb = ar;                           }
  void    SetMissingPreHb(UInt_t ar)                  { fMissingPreHb = ar;                           }
  void    SetUnnecessaryPreLb(UInt_t ar)              { fUnnecessaryPreLb = ar;                           }
  void    SetUnnecessaryPreHb(UInt_t ar)              { fUnnecessaryPreHb = ar;                           }
  
  

 protected:
  TString fSide; // side of the control box, either A, B or C 
  TInt_t  fPrimary; // 1 if its the primary control box, 2 for backup

  UInt_t  fPreToNextCycles; // value from the PreToNextCycles tag in the pt box type CbB
  UInt_t  fL0ToNextCycles; // value from the L0ToNextCycles tag in the pt box type CbB
  UInt_t  fL1ToNextCycles; // value from the L1ToNextCycles tag in the pt box type CbB
  UInt_t  fClkTtcexShifts; // value from the ClkTtcexShifts tag in the pt box type CbB
  UInt_t  fDisableGtuBusy; // value from the DisableGtuBusy tag in the pt box type CbB
  UInt_t  fDisableSorEorBusy; // value from the DisableSorEorBusy tag in the pt box type CbB
  UInt_t  fScintEn; // value from the ScintEn tag in the pt box type CbB
  UInt_t  fLutAsPre; // value from the LutAsPre tag in the pt box type CbB
  UInt_t  fCtpAsPre; // value from the CtpAsPre tag in the pt box type CbB
  UInt_t  fL0En; // value from the L0En tag in the pt box type CbB
  UInt_t  fL1En; // value from the L1En tag in the pt box type CbB
  UInt_t  fLutL0ToCtpEn; // value from the LutL0ToCtpEn tag in the pt box type CbB
  UInt_t  fLutPreToCtpEn; // value from the LutPreToCtpEn tag in the pt box type CbB
  UInt_t  fSoftwTrigToCtpEn; // value from the SoftwTrigToCtpEn tag in the pt box type CbB
  UInt_t  fL0Autogenerate; // value from the L0Autogenerate tag in the pt box type CbB
  UInt_t  fL1Autogenerate; // value from the L1Autogenerate tag in the pt box type CbB
  UInt_t  fChannelBDisable; // value from the ChannelBDisable tag in the pt box type CbB
  UInt_t  fTtcexClkDisable; // value from the TtcexClkDisable tag in the pt box type CbB
  UInt_t  fPreAsL0En; // value from the PreAsL0En tag in the pt box type CbB
  UInt_t  fNoTriggerToSm; // value from the NoTriggerToSm tag in the pt box type CbB
  UInt_t  fSoftwTrigEn; // value from the SoftwTrigEn tag in the pt box type CbB
  UInt_t  fSoftwTrigPattern; // value from the SoftwTrigPattern tag in the pt box type CbB
  UInt_t  fSychrInput; // value from the SychrInput tag in the pt box type CbB
  UInt_t  fRandomTriggerThr; // value from the RandomTriggerThr tag in the pt box type CbB
  UInt_t  fPtChDelayCbA; // value from the PtChDelayCbA tag in the pt box type CbB
  UInt_t  fPtChDelayCbC; // value from the PtChDelayCbC tag in the pt box type CbB
  UInt_t  fPtChDelayGtu; // value from the PtChDelayGtu tag in the pt box type CbB
  UInt_t  fPtChDelayT0; // value from the PtChDelayT0 tag in the pt box type CbB
  UInt_t  fPtChDelayT0; // value from the PtChDelayT0 tag in the pt box type CbB
  UInt_t  fPtChDelayT1; // value from the PtChDelayT1 tag in the pt box type CbB
  UInt_t  fPtChDelayT2; // value from the PtChDelayT2 tag in the pt box type CbB
  UInt_t  fPtChDelayT3; // value from the PtChDelayT3 tag in the pt box type CbB
  UInt_t  fPtChDelayT4; // value from the PtChDelayT4 tag in the pt box type CbB
  UInt_t  fPtChDelayT5; // value from the PtChDelayT5 tag in the pt box type CbB
  UInt_t  fPtChDelayT6; // value from the PtChDelayT6 tag in the pt box type CbB
  UInt_t  fPtChDelayT7; // value from the PtChDelayT7 tag in the pt box type CbB
  UInt_t  fClkLb; // value from the Clk low bit tag in the pt box type CbB
  UInt_t  fClkHb; // value from the Clk high bit tag in the pt box type CbB
  UInt_t  fPulseToSmLb; // value from the PulseToSm low bit tag in the pt box type CbB
  UInt_t  fPulseToSmHb; // value from the PulseToSm high bit tag in the pt box type CbB
  UInt_t  fDeadTimeLb; // value from the DeadTime low bit tag in the pt box type CbB
  UInt_t  fDeadTimeHb; // value from the DeadTime high bit tag in the pt box type CbB
  UInt_t  fPrePulsesLb; // value from the PrePulses low bit tag in the pt box type CbB
  UInt_t  fPrePulsesHb; // value from the PrePulses high bit tag in the pt box type CbB
  UInt_t  fL0CtpLb; // value from the L0Ctp low bit tag in the pt box type CbB
  UInt_t  fL0CtpHb; // value from the L0Ctp high bit tag in the pt box type CbB
  UInt_t  fL1CtpLb; // value from the L1Ctp low bit tag in the pt box type CbB
  UInt_t  fL1CtpHb; // value from the L1Ctp high bit tag in the pt box type CbB
  UInt_t  fL0LutLb; // value from the L0Lut low bit tag in the pt box type CbB
  UInt_t  fL0LutHb; // value from the L0Lut high bit tag in the pt box type CbB
  UInt_t  fNotTinToCtpPreLb; // value from the NotTinToCtpPre low bit tag in the pt box type CbB
  UInt_t  fNotTinToCtpPreHb; // value from the NotTinToCtpPre high bit tag in the pt box type CbB
  UInt_t  fPatternMatchCbALb; // value from the PatternMatchCbA low bit tag in the pt box type CbB
  UInt_t  fPatternMatchCbAHb; // value from the PatternMatchCbA high bit tag in the pt box type CbB
  UInt_t  fPatternMatchCbCLb; // value from the PatternMatchCbC low bit tag in the pt box type CbB
  UInt_t  fPatternMatchCbCHb; // value from the PatternMatchCbC high bit tag in the pt box type CbB
  UInt_t  fPatternMatchCbTlmuLb; // value from the PatternMatchCbTlmu low bit tag in the pt box type CbB
  UInt_t  fPatternMatchCbTlmuHb; // value from the PatternMatchCbTlmu high bit tag in the pt box type CbB
  UInt_t  fTriggerSynchrLb; // value from the TriggerSynchr low bit tag in the pt box type CbB
  UInt_t  fTriggerSynchrHb; // value from the TriggerSynchr high bit tag in the pt box type CbB
  UInt_t  fTriggerToSmLb; // value from the TriggerToSm low bit tag in the pt box type CbB
  UInt_t  fTriggerToSmHb; // value from the TriggerToSm high bit tag in the pt box type CbB
  UInt_t  fPreLutLb; // value from the PreLut low bit tag in the pt box type CbB
  UInt_t  fPreLutHb; // value from the PreLut high bit tag in the pt box type CbB
  UInt_t  fMissingPreLb; // value from the MissingPre low bit tag in the pt box type CbB
  UInt_t  fMissingPreHb; // value from the MissingPre high bit tag in the pt box type CbB
  UInt_t  fUnnecessaryPreLb; // value from the UnnecessaryPre low bit tag in the pt box type CbB
  UInt_t  fUnnecessaryPreHb; // value from the UnnecessaryPre high bit tag in the pt box type CbB


  ClassDef(AliTRDCalDCSPTRCbb,1)      //  TRD calibration class for TRD GTU parameters

};
#endif
