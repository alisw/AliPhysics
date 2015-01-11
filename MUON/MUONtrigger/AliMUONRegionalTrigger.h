#ifndef ALIMUONREGIONALTRIGGER_H
#define ALIMUONREGIONALTRIGGER_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
// Revision of includes 07/05/2004

/// \ingroup trigger
/// \class AliMUONRegionalTrigger
/// \brief Reconstructed regional Trigger object
//  Author Ch. Finck

#include <TObject.h>

class AliMUONRegionalTrigger : public TObject {
 public:
  AliMUONRegionalTrigger();
  AliMUONRegionalTrigger(const AliMUONRegionalTrigger& rhs); // copy constructor !
  virtual ~AliMUONRegionalTrigger();
  AliMUONRegionalTrigger& operator=(const AliMUONRegionalTrigger& rhs); 

  // getter methods
  /// Return regional id 
  Int_t    GetId()        const {return fId;}
  /// Return local output Lpt & Hpt
  UInt_t   GetLocalOutput(Int_t n) const {return fLocalOutput[n];}
  /// Return local mask
  UShort_t GetLocalMask() const {return fLocalMask;}
  /// Return regional output (single muon:2, unlike sign:1, like sign:1)
  Char_t   GetOutput()    const {return fOutput;}

  // setter methods
  /// Set regional id 
  void  SetId(Int_t d)           {fId = d;}
  /// Set local output Lpt & Hpt
  void  SetLocalOutput(UInt_t local, Int_t n) {fLocalOutput[n] = local;}
  /// Set local mask
  void  SetLocalMask(UShort_t m) {fLocalMask = m;}
  /// Set regional output (single muon:2, unlike sign:1, like sign:1)
  void  SetOutput(Char_t o)      {fOutput = o;}

  virtual void Print(Option_t* opt="") const;
  
private:
  Int_t    fId;              ///< regional id 
  UInt_t   fLocalOutput[2];  ///< local output Lpt & Hpt
  UShort_t fLocalMask;       ///< local mask
  UChar_t  fOutput;          ///< regional output (single muon:2, unlike sign:1, like sign:1) 


  ClassDef(AliMUONRegionalTrigger,1)  // reconstructed regional Trigger object
};
#endif






