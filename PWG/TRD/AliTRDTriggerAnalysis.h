/* Copyright(c) 2013, ALICE Experiment at CERN, All rights reserved.      *
 * See cxx source for full Copyright notice                               */

// evaluate TRD trigger conditions,
// potentially with hardened conditions to remove
// triggers caused by conversions of low-pt photons
// at large radii
//
// Author: Jochen Klein <jochen.klein@cern.ch>

#ifndef ALITRDTRIGGERANALYSIS_H
#define ALITRDTRIGGERANALYSIS_H

class AliVEvent;

class AliTRDTriggerAnalysis : public TObject
{
public:
  AliTRDTriggerAnalysis();
  ~AliTRDTriggerAnalysis();

  enum TRDTrigger_t { kHCO = 0, kHJT, kHSE, kHQU, kHEE, kHlast };

  enum JetTriggerMode_t { kHJTDefault = 0, kHJTWindowZPhi };

  void ResetTriggers();
  Bool_t CalcTriggers(const AliVEvent* event);

  Bool_t IsFired(TRDTrigger_t trg) const {
    Obsolete("IsFired(...) is deprecated, use CheckCondition instead",
	     "now", "asap");
    return CheckCondition(trg);
  }

  Bool_t HasTriggeredConfirmed(TRDTrigger_t trg) const {
    return (HasTriggered(trg) && CheckCondition(trg));
  }
  Bool_t HasTriggered(TRDTrigger_t trg) const {
    return (fTriggerClasses & (1 << trg));
  }
  Bool_t HasFired(TRDTrigger_t trg) const {
    return (fTriggerInputs & (1 << trg));
  }
  Bool_t CheckCondition(TRDTrigger_t trg) const {
    return (fTriggerFlags[2 * trg] | fTriggerFlags[2 * trg + 1]);
  }
  Bool_t CheckCondition(TRDTrigger_t trg, Int_t stack) const {
    Int_t idx = 2 * trg + (stack / 64);
    Int_t bit = stack % 64;
    return (fTriggerFlags[idx] & (1ULL << bit));
  }

  Bool_t CheckTrgFlags(Int_t bit, Int_t sector) const {
    return (fTriggerContribs[sector] & (1 << bit));
  }

  void SetRequireMatch(Bool_t val) { fRequireMatch = val; }
  Bool_t GetRequireMatch() const { return fRequireMatch; }

  void SetRequireMatchElectron(Bool_t val) { fRequireMatchElectron = val; }
  Bool_t GetRequireMatchElectron() const { return fRequireMatchElectron; }

  void SetRequireInTime(Bool_t val) { fRequireInTime = val; }
  Bool_t GetRequireInTime() const { return fRequireInTime; }

  void SetJetTriggerMode(Int_t mode) { fJetTriggerMode = mode; }
  Int_t GetJetTriggerMode() const { return fJetTriggerMode; }

  void SetVerbosity(UChar_t val) { fVerbosity = val; }
  UChar_t GetVerbosity() const { return fVerbosity; }

protected:
  void MarkClass(TRDTrigger_t trg) { fTriggerClasses |= (1 << trg); }
  void MarkInput(TRDTrigger_t trg) { fTriggerInputs |= (1 << trg); }
   void MarkCondition(TRDTrigger_t trg, Int_t stack)
  {
    Int_t idx = 2 * trg + (stack / 64);
    Int_t bit = stack % 64;
    fTriggerFlags[idx] |= (1ULL << bit);
  }


  static const Int_t fgkNstacks = 90; // no. of TRD stacks (global)
  ULong64_t fTriggerFlags[2 * kHlast];   //! internal representation of condition checks
  UChar_t fTriggerInputs;  //! internal representation of trigger inputs
  UChar_t fTriggerClasses; //! internal representation of trigger classes

  // configuration
  UChar_t fVerbosity;      // verbosity level
  Bool_t fRequireMatch;	   // require a matched global track
			   // for all conditions
  Bool_t fRequireMatchElectron;	// require a matched global track
				// for the electron conditions
  Bool_t fRequireInTime;	// require the tracks to be in time
  Int_t  fJetTriggerMode;       // select mode for jet trigger
				// 0: default (stack-wise counting, as hw)
				// 1: count in overlapping windows of stack size

  // trigger thresholds
  UChar_t fTRDlayerMaskEl;      // mask for tracklet requirements
  UChar_t fTRDnTrackletsEl;     // min. number of tracklets
  Float_t fTRDptHSE;            // pt threshold for HSE trigger
  UChar_t fTRDpidHSE;           // PID threshold for HSE trigger
  Float_t fTRDptHQU;            // pt threshold for HQU trigger
  UChar_t fTRDpidHQU;           // PID threshold for HQU trigger
  Float_t fTRDptHEE;            // pt threshold for HEE trigger
  UChar_t fTRDpidHEE;           // PID threshold for HEE trigger
  UChar_t fTRDminSectorHEE;     // min sector for HEE trigger
  UChar_t fTRDmaxSectorHEE;     // max sector for HEE trigger
  Float_t fTRDptHJT;            // pt threshold for HJT trigger
  UChar_t fTRDnHJT;             // no of track threshold for HJT trigger

  UInt_t fTriggerContribs[18]; // temporary for debugging !!!

  ClassDef(AliTRDTriggerAnalysis, 1);
};

#endif
