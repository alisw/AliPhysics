#ifndef ALIMUONLOCALTRIGGER_H
#define ALIMUONLOCALTRIGGER_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
// Revision of includes 07/05/2004

#include <TObject.h>
#include <TArrayI.h>

class AliMUONLocalTrigger : public TObject {
 public:
  AliMUONLocalTrigger();
  AliMUONLocalTrigger(const AliMUONLocalTrigger& rhs); // copy constructor !
  AliMUONLocalTrigger(const Int_t* localtr, const TArrayI& digits);
  virtual ~AliMUONLocalTrigger(){;}
  AliMUONLocalTrigger& operator=(const AliMUONLocalTrigger& rhs); 

  // getter methods
  Int_t LoCircuit() const {return fLoCircuit;}
  Int_t LoStripX() const {return fLoStripX;}   
  Int_t LoDev() const {return fLoDev;}
  Int_t LoStripY() const {return fLoStripY;}
  Int_t LoLpt() const {return fLoLpt;}
  Int_t LoHpt() const {return fLoHpt;}
  Int_t LoApt() const {return fLoApt;}

   
  UShort_t GetX1Pattern() const {return fX1Pattern;}
  UShort_t GetX2Pattern() const {return fX2Pattern;}
  UShort_t GetX3Pattern() const {return fX3Pattern;}
  UShort_t GetX4Pattern() const {return fX4Pattern;}

  UShort_t GetY1Pattern() const {return fY1Pattern;}
  UShort_t GetY2Pattern() const {return fY2Pattern;}
  UShort_t GetY3Pattern() const {return fY3Pattern;}
  UShort_t GetY4Pattern() const {return fY4Pattern;}

  Char_t GetLoDecision();

  // setter methods
  void SetLoCircuit(Int_t loCir) {fLoCircuit = loCir;}
  void SetLoStripX(Int_t loStrX) {fLoStripX = loStrX;}   
  void SetLoDev(Int_t loDev)     {fLoDev = loDev;}
  void SetLoStripY(Int_t loStrY) {fLoStripY = loStrY;}
  void SetLoLpt(Int_t loLpt)     {fLoLpt = loLpt;}
  void SetLoHpt(Int_t loHpt)     {fLoHpt = loHpt;}
  void SetLoApt(Int_t loApt)     {fLoApt = loApt;}

  void SetX1Pattern(UShort_t pat) {fX1Pattern = pat;}
  void SetX2Pattern(UShort_t pat) {fX2Pattern = pat;}
  void SetX3Pattern(UShort_t pat) {fX3Pattern = pat;}
  void SetX4Pattern(UShort_t pat) {fX4Pattern = pat;}

  void SetY1Pattern(UShort_t pat) {fY1Pattern = pat;}
  void SetY2Pattern(UShort_t pat) {fY2Pattern = pat;}
  void SetY3Pattern(UShort_t pat) {fY3Pattern = pat;}
  void SetY4Pattern(UShort_t pat) {fY4Pattern = pat;}

  // data link
  Int_t NumberOfDigits() const { return fDigits.GetSize(); }
  Int_t GetDigitNumber(Int_t i) const { return fDigits[i]; }
  void GetDigit(Int_t i, Int_t& chamber, Int_t& cathode, Int_t& digit) const;

  static Int_t EncodeDigitNumber(Int_t chamber, Int_t cathode, Int_t digit);
  static void DecodeDigitNumber(Int_t digitnumber, Int_t& chamber, Int_t& cathode, Int_t& digit);

private:
  Int_t fLoCircuit; // circuit number 
  Int_t fLoStripX;  // X strip in MT11 
  Int_t fLoDev;     // deviation 
  Int_t fLoStripY;  // Y strip in MT11 
  Int_t fLoLpt;     // Low pt  0 : nothing, 1 : Minus, 2 : Plus, 3 : Undef
  Int_t fLoHpt;     // High pt 0 : nothing, 1 : Minus, 2 : Plus, 3 : Undef
  Int_t fLoApt;     // All pt  0 : nothing, 1 : Minus, 2 : Plus, 3 : Undef

  UShort_t fX1Pattern; // X and Y strip pattern for each chamber
  UShort_t fX2Pattern;
  UShort_t fX3Pattern;
  UShort_t fX4Pattern;
 
  UShort_t fY1Pattern;
  UShort_t fY2Pattern;
  UShort_t fY3Pattern;
  UShort_t fY4Pattern;

  Char_t fLoDecision; // local decision word (4 bits)

  TArrayI fDigits;    // List of digit numbers from which this object was created.

  ClassDef(AliMUONLocalTrigger,2)  // reconstructed Local Trigger object
};
#endif






