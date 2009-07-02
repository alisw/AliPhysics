#ifndef ALITRIGGERIR_H
#define ALITRIGGERIR_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

///////////////////////////////////////////////////////////////////////////////
//
//  Class represents CTP interaction record
//
//  The definition of the IR follows the ALICE internal note:
//  ALICE-INT-2002-010
//  The CTP raw-data payload will contain IRs within +- 1 orbit
//  around the triggered event.
//  The same IRs are propagated to the ESD (and AOD).
//
//  cvetan.cheshkov@cern.ch 10/07/2008
//
///////////////////////////////////////////////////////////////////////////////

#include <TObject.h>

class AliTriggerIR : public TObject {

public:
  AliTriggerIR();
  AliTriggerIR(UInt_t orbit, UInt_t nwords, UInt_t *words, Bool_t incomplete = kFALSE, Bool_t transerr = kFALSE);

  AliTriggerIR(const AliTriggerIR &rec);
  AliTriggerIR& operator= (const AliTriggerIR& rec);

  virtual   ~AliTriggerIR();

  //  Setters
  void SetOrbit(UInt_t orbit) {fOrbit=orbit;}
  void SetIncomplete(Bool_t flag) {fIncomplete=flag;}
  void SetTransErr(Bool_t flag) {fTransErr=flag;}

  //  Getters
  UInt_t GetOrbit() const {return fOrbit;}
  UInt_t GetNWord() const {return fNWord;}
  Bool_t* GetInt1s() const {return fInt1;}
  Bool_t* GetInt2s() const {return fInt2;}
  UShort_t* GetBCs() const {return fBC;}
  Bool_t GetIncomplete() const {return fIncomplete;}
  Bool_t GetTransErr() const {return fTransErr;}
  virtual void   Print( const Option_t* opt ="" ) const;

private:
  UInt_t    fOrbit;        // Orbit number
  UInt_t    fNWord;        // Number of recorded interaction signals
  Bool_t   *fInt1;         //[fNWord] signals for interaction 1
  Bool_t   *fInt2;         //[fNWord] signals for interaction 2
  UShort_t *fBC;           //[fNWord] bunch-crossing number
  Bool_t    fIncomplete;   // flag which says if the IR is incomplete or not
  Bool_t    fTransErr;     // flag which says if there was a transmission error (gap) or not

  ClassDef( AliTriggerIR, 1 )  // Trigger Interaction Record (one per orbit)
};

#endif
