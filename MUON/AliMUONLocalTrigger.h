#ifndef ALIMUONLOCALTRIGGER_H
#define ALIMUONLOCALTRIGGER_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
// Revision of includes 07/05/2004

#include <TObject.h>

class AliMUONLocalTrigger : public TObject {
 public:
  AliMUONLocalTrigger();
  AliMUONLocalTrigger(const AliMUONLocalTrigger& rhs); // copy constructor !
  AliMUONLocalTrigger(Int_t* localtr);
  virtual ~AliMUONLocalTrigger(){;}
  AliMUONLocalTrigger& operator=(const AliMUONLocalTrigger& rhs); 

  Int_t LoCircuit() const {return fLoCircuit;}; 
  Int_t LoStripX() const {return fLoStripX;};    
  Int_t LoDev() const {return fLoDev;};     
  Int_t LoStripY() const {return fLoStripY;};  
  Int_t LoLpt() const {return fLoLpt;};     
  Int_t LoHpt() const {return fLoHpt;};     
  Int_t LoApt() const {return fLoApt;}; 

   
  UShort_t GetX1Pattern() const {return fX1Pattern;}
  UShort_t GetX2Pattern() const {return fX2Pattern;}
  UShort_t GetX3Pattern() const {return fX3Pattern;}
  UShort_t GetX4Pattern() const {return fX4Pattern;}

  UShort_t GetY1Pattern() const {return fY1Pattern;}
  UShort_t GetY2Pattern() const {return fY2Pattern;}
  UShort_t GetY3Pattern() const {return fY3Pattern;}
  UShort_t GetY4Pattern() const {return fY4Pattern;}

  Char_t GetLoDecision();

  ClassDef(AliMUONLocalTrigger,2)  // reconstructed Local Trigger object

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
};
#endif






