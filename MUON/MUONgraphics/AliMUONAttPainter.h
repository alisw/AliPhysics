#ifndef ALIMUONATTPAINTER_H
#define ALIMUONATTPAINTER_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup graphics
/// \class AliMUONAttPainter
/// \brief Basic attributes shared by all painters
/// 
// Author Laurent Aphecetche, Subatech

#ifndef ROOT_TObject
#  include "TObject.h"
#endif
#ifndef ROOT_TString
#  include "TString.h"
#endif

class AliMUONAttPainter : public TObject
{
public:
  
  /// Internal status bits
  enum EBits {
    kIsCathode0         = BIT(14),
    kIsCathode1         = BIT(15),
    kIsBendingPlane     = BIT(16),
    kIsNonBendingPlane  = BIT(17),
    kIsFrontView        = BIT(18),
    kIsBackView         = BIT(19),
    kIsCathodeAndPlaneMutuallyExclusive = BIT(20),
    kIsValid            = BIT(21),
    kIsSinglePainter    = BIT(22),
    kIsCathodeAndPlaneDisabled = BIT(23)
  };
  
  AliMUONAttPainter();
  virtual ~AliMUONAttPainter();
  
  /// Return our name
  virtual const char* GetName() const { return Name(); }
  
  TString Name() const;
  
  TString CathodeName() const;
  
  TString ViewPointName() const;

  TString PlaneName() const;
  
  /// Whether cathode & plane are disabled
  Bool_t IsCathodeAndPlaneDisabled() const { return TestBit(kIsCathodeAndPlaneDisabled); }
  
  /// Whether we are representing bending plane
  Bool_t IsBendingPlane() const { return TestBit(kIsBendingPlane); }
  
  /// Whether we are representing cathode 0
  Bool_t IsCathode0() const { return TestBit(kIsCathode0); }
  
  /// Whether we are representing cathode 1
  Bool_t IsCathode1() const { return TestBit(kIsCathode1); }
  
  /// Whether we can select both cathode and plane
  Bool_t IsCathodeAndPlaneMutuallyExclusive() const { return TestBit(kIsCathodeAndPlaneMutuallyExclusive); }
  
  /// Whether we are defined by cathode
  Bool_t IsCathodeDefined() const { return IsCathode0() || IsCathode1(); }
  
  /// Whether we are representing non bending plane
  Bool_t IsNonBendingPlane() const { return TestBit(kIsNonBendingPlane); }
  
  /// Whether we are defined by plane
  Bool_t IsPlaneDefined() const { return IsBendingPlane() || IsNonBendingPlane(); }
  
  /// Whether we are valid
  Bool_t IsValid() const { return TestBit(kIsValid); }
  
  void Invert();
  
  /// Set our cathode states
  void SetCathode(Bool_t cath0, Bool_t cath1) { SetBit(kIsCathode0,cath0); SetBit(kIsCathode1,cath1); }
  
  /// Set our plane states
  void SetPlane(Bool_t bending, Bool_t nonBending) { SetBit(kIsBendingPlane,bending); SetBit(kIsNonBendingPlane,nonBending); }
  
  /// Set single status
  void SetSingle(Bool_t value) { SetBit(kIsSinglePainter,value); }
  
  /// Whether the painter is to be represented from front (as seen from IP)
  Bool_t IsFrontView() const { return TestBit(kIsFrontView); }
  
  /// Whether the painter is to be represented from back (as seen from IP)
  Bool_t IsBackView() const { return TestBit(kIsBackView); }
  
  /// Set view point
  void SetViewPoint(Bool_t front, Bool_t back) { SetBit(kIsFrontView,front); SetBit(kIsBackView,back); }
  
  /// Set mutually exclusive flag
  void SetCathodeAndPlaneMutuallyExclusive(Bool_t value) { SetBit(kIsCathodeAndPlaneMutuallyExclusive,value); }
  
  /// Set valid flag
  void SetValid(Bool_t value) { SetBit(kIsValid,value); }
  
  /// Whether we represent attributes of a single painter (if false, means it's a painter group)
  Bool_t IsSinglePainter() const { return TestBit(kIsSinglePainter); }
  
  /// Set cathode & plane disable flag
  void SetCathodeAndPlaneDisabled(Bool_t value) { SetBit(kIsCathodeAndPlaneDisabled,value); }

  void Print(Option_t* opt="") const;
  
  ClassDef(AliMUONAttPainter,1) // Basic attributes of painters
};

#endif
