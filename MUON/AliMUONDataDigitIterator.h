#ifndef ALIMUONDATADIGITITERATOR_H
#define ALIMUONDATADIGITITERATOR_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup base
/// \class AliMUONDataDigitIterator
/// \brief Iterator on digits (handled by AliMUONData).
/// 
/// \author Laurent Aphecetche

#ifndef ALIMUONVDATAITERATOR_H
#  include "AliMUONVDataIterator.h"
#endif

class AliMUONData;
class TClonesArray;

class AliMUONDataDigitIterator : public AliMUONVDataIterator
{
public:
  AliMUONDataDigitIterator(const AliMUONData* data, Int_t firstChamber, Int_t lastChamber);
  AliMUONDataDigitIterator(const AliMUONDataDigitIterator& rhs);
  AliMUONDataDigitIterator& operator=(const AliMUONDataDigitIterator& rhs);
  virtual ~AliMUONDataDigitIterator() {}
    
  TObject* Next();
  
  void Reset(); 
  
  Bool_t Remove();
  
private:
    void CopyTo(AliMUONDataDigitIterator& destination) const;
  
private:
    const AliMUONData* fData; //! Pointer to data accessor
  Int_t fFirstChamber; //! First chamber to iterate on
  Int_t fLastChamber; //! Last chamber to iterate on      
  TClonesArray* fDigits; //! Digits of the current chamber
  Int_t fCurrentDigit; //! Current position within fDigits array
  Int_t fCurrentChamber; //! Current chamber
  
  ClassDef(AliMUONDataDigitIterator,0)
};      

#endif
