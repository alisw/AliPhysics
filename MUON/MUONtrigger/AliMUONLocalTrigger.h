#ifndef ALIMUONLOCALTRIGGER_H
#define ALIMUONLOCALTRIGGER_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
// Revision of includes 07/05/2004

/// \ingroup trigger
/// \class AliMUONLocalTrigger
/// \brief Reconstructed Local Trigger object
//  Author Ph. Crochet

#include <TObject.h>
#include "AliMUONRawStreamTriggerHP.h"

class AliMUONLocalStruct;
class TArrayS;

class AliMUONLocalTrigger : public TObject {
 public:
  AliMUONLocalTrigger();
  AliMUONLocalTrigger(const AliMUONLocalTrigger& rhs); // copy constructor !
  virtual ~AliMUONLocalTrigger();
  AliMUONLocalTrigger& operator=(const AliMUONLocalTrigger& rhs); 

  // getter methods
  //
        /// Return Circuit number
  Int_t LoCircuit() const {return fLoCircuit;}
        /// Return X strip in MT11
  Int_t LoStripX() const {return fLoStripX;}   
        /// Return Deviation
  Int_t LoDev() const {return fLoDev;}
        /// Return Sign of Deviation
  Int_t LoSdev() const {return fLoSdev;}
        /// Return Trig Y
  Int_t LoTrigY() const {return fLoTrigY;}
        /// Return Y strip in MT11
  Int_t LoStripY() const {return fLoStripY;}
        /// Return Low pt
  Int_t LoLpt() const {return fLoLpt;}
        /// Return High p
  Int_t LoHpt() const {return fLoHpt;}

           /// Return X strip pattern for chamber 11
  UShort_t GetX1Pattern() const {return fX1Pattern;}
           /// Return X strip pattern for chamber 12  
  UShort_t GetX2Pattern() const {return fX2Pattern;}
           /// Return X strip pattern for chamber 21 
  UShort_t GetX3Pattern() const {return fX3Pattern;}
           /// Return X strip pattern for chamber 22
  UShort_t GetX4Pattern() const {return fX4Pattern;}
           /// return X pattern array
   void    GetXPattern(TArrayS& array) const;


           /// Return Y strip pattern for chamber 11 
  UShort_t GetY1Pattern() const {return fY1Pattern;}
           /// Return Y strip pattern for chamber 12
  UShort_t GetY2Pattern() const {return fY2Pattern;}
           /// Return Y strip pattern for chamber 21
  UShort_t GetY3Pattern() const {return fY3Pattern;}
           /// Return Y strip pattern for chamber 22
  UShort_t GetY4Pattern() const {return fY4Pattern;}
           /// return Y pattern array
   void    GetYPattern(TArrayS& array) const;

  Char_t GetLoDecision() const;
  
  Int_t GetDeviation() const;
  void  SetDeviation(Int_t deviation);

  /// Trigger response X strips
  Bool_t IsTrigX() const;
  /// Trigger response Y strips
  Bool_t IsTrigY() const;

  /// Hit pattern from the re-calculated trigger response after removing chambers one-by-one
  UShort_t GetHitPatternFromResponse() const {return fHitPatternFromResponse; }

  /// Triggers from the re-calculated trigger response after removing chambers one-by-one
  UChar_t GetTriggerWithoutChamber() const {return fTriggerWithoutChamber; }

  // setter methods
  //
           /// Set Circuit number
  void SetLoCircuit(Int_t loCir) {fLoCircuit = loCir;}
           /// Set X strip in MT11
  void SetLoStripX(Int_t loStrX) {fLoStripX = loStrX;}   
           /// Set Deviation
  void SetLoDev(Int_t loDev)     {fLoDev = loDev;}
             /// Set Sign of Deviation
  void SetLoSdev(Int_t loSdev)   {fLoSdev = loSdev;}
           /// Set Trig Y
  void SetLoTrigY(Int_t loTrigY) {fLoTrigY = loTrigY;}
           /// Set Y strip in MT11
  void SetLoStripY(Int_t loStrY) {fLoStripY = loStrY;}
           /// Set Low pt
  void SetLoLpt(Int_t loLpt)     {fLoLpt = loLpt;}
           /// Set High pt
  void SetLoHpt(Int_t loHpt)     {fLoHpt = loHpt;}
   
           /// Set X strip pattern for chamber 11
  void SetX1Pattern(UShort_t pat) {fX1Pattern = pat;}
           /// Set X strip pattern for chamber 12
  void SetX2Pattern(UShort_t pat) {fX2Pattern = pat;}
           /// Set X strip pattern for chamber 21
  void SetX3Pattern(UShort_t pat) {fX3Pattern = pat;}
           /// Set X strip pattern for chamber 22
  void SetX4Pattern(UShort_t pat) {fX4Pattern = pat;}

           /// Set Y strip pattern for chamber 11
  void SetY1Pattern(UShort_t pat) {fY1Pattern = pat;}
           /// Set Y strip pattern for chamber 12
  void SetY2Pattern(UShort_t pat) {fY2Pattern = pat;}
           /// Set Y strip pattern for chamber 21
  void SetY3Pattern(UShort_t pat) {fY3Pattern = pat;}
           /// Set Y strip pattern for chamber 22
  void SetY4Pattern(UShort_t pat) {fY4Pattern = pat;}

  void SetLocalStruct(Int_t loCircuit, AliMUONLocalStruct& localStruct);
  void SetLocalStruct(Int_t loCircuit, const AliMUONRawStreamTriggerHP::AliLocalStruct& localStruct);

  Bool_t IsNull() const;
  
  virtual void Print(Option_t* opt="") const;
  
  virtual const char* GetName() const;

  /// The board would provide a trigger even after removing chamber ich [0,3]
  void SetTriggerWithoutChamber(Int_t ich){ fTriggerWithoutChamber |= 1 << (3 - ich); }

  /// Trigger algorithm did NOT find hit in the specified chamber and cathode
  void SetNoHitInPlane(Int_t icath, Int_t ich){ fHitPatternFromResponse &= ~(1 << (7 - 4*icath - ich)); }

private:
  Int_t fLoCircuit; ///< Circuit number 
  Int_t fLoStripX;  ///< X strip in MT11 
  Int_t fLoDev;     ///< Deviation
  Int_t fLoSdev;    ///< Sign of Deviation 
  Int_t fLoTrigY;   ///< Trig Y
  Int_t fLoStripY;  ///< Y strip in MT11 
  Int_t fLoLpt;     ///< Low pt  0 : nothing, 1 : Minus, 2 : Plus, 3 : Undef
  Int_t fLoHpt;     ///< High pt 0 : nothing, 1 : Minus, 2 : Plus, 3 : Undef
  
  UShort_t fX1Pattern; ///< X strip pattern for chamber 11
  UShort_t fX2Pattern; ///< X strip pattern for chamber 12
  UShort_t fX3Pattern; ///< X strip pattern for chamber 21
  UShort_t fX4Pattern; ///< X strip pattern for chamber 22

  UShort_t fY1Pattern; ///< Y strip pattern for chamber 11
  UShort_t fY2Pattern; ///< Y strip pattern for chamber 12
  UShort_t fY3Pattern; ///< Y strip pattern for chamber 21
  UShort_t fY4Pattern; ///< Y strip pattern for chamber 22

  UShort_t fHitPatternFromResponse; ///<  Fired plane according to re-computed repsonse
  UChar_t fTriggerWithoutChamber; ///< Pattern of triggers after chamber removal

  ClassDef(AliMUONLocalTrigger,5)  // reconstructed Local Trigger object
};
#endif






