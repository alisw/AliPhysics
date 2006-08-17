#ifndef ALIMUONREGIONALTRIGGER_H
#define ALIMUONREGIONALTRIGGER_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
// Revision of includes 07/05/2004

/// \ingroup base
/// \class AliMUONRegionalTrigger
/// \brief Reconstructed regional Trigger object
//  Author Ch. Finck

#include <TObject.h>

class AliMUONRegionalTrigger : public TObject {
 public:
  AliMUONRegionalTrigger();
  AliMUONRegionalTrigger(const AliMUONRegionalTrigger& rhs); // copy constructor !
  virtual ~AliMUONRegionalTrigger(){;}
  AliMUONRegionalTrigger& operator=(const AliMUONRegionalTrigger& rhs); 

  // getter methods
  Int_t    GetId()        const {return fId;}
  UInt_t   GetLocalOutput(Int_t n) const {return fLocalOutput[n];}
  UShort_t GetLocalMask() const {return fLocalMask;}
  Char_t   GetOutput()    const {return fOutput;}

  // setter methods
  void  SetId(Int_t d)           {fId = d;}
  void  SetLocalOutput(UInt_t local, Int_t n) {fLocalOutput[n] = local;}
  void  SetLocalMask(UShort_t m) {fLocalMask = m;}
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






