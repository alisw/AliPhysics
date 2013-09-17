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

  enum TRDTrigger_t { kHCO = 0, kHJT, kHSE, kHQU, kHEE };

  void ResetTriggers() { fTriggerFlags = fTriggerInputs = fTriggerClasses = 0; }
  Bool_t CalcTriggers(const AliVEvent* event);

  Bool_t IsFired(TRDTrigger_t trg) const { return (fTriggerFlags & (1 << trg)); }

protected:
  void Fire(TRDTrigger_t trg) { fTriggerFlags |= (1 << trg); }

  // configuration
  UChar_t fTriggerFlags;   // internal representation of trigger decisions
  UChar_t fTriggerInputs;  // internal representation of trigger decisions
  UChar_t fTriggerClasses; // internal representation of trigger decisions
  Bool_t fRequireMatch;	   // require a matched global track
			   // for all conditions
  Bool_t fRequireMatchElectron;	// require a matched global track
				// for the electron conditions
  Bool_t fRequireInTime;	// require the tracks to be in time

  // trigger thresholds
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

  ClassDef(AliTRDTriggerAnalysis, 1);
};

#endif
