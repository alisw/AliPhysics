/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$
// $MpId: AliMpIntPair.h,v 1.4 2005/08/26 15:43:36 ivana Exp $

/// \ingroup basic
/// \class AliMpIntPair
/// \brief A pair of integers.
///
/// The pair created by the default constructor is in invalide state,
/// setting one of values changes the state to valid.
///
/// Authors: David Guez, Ivana Hrivnacova; IPN Orsay

#ifndef ALI_MP_INT_PAIR_H
#define ALI_MP_INT_PAIR_H

#include <TObject.h>

class AliMpIntPair : public TObject
{
 public:
  AliMpIntPair(Int_t ix,Int_t iy);
  AliMpIntPair(Int_t ix,Int_t iy, Bool_t validity);
  AliMpIntPair(const AliMpIntPair& right);
  AliMpIntPair();
  virtual ~AliMpIntPair();

  // operators  
  Bool_t operator <  (const AliMpIntPair& pos2) const;
  Bool_t operator == (const AliMpIntPair& pos2) const;
  Bool_t operator != (const AliMpIntPair& pos2) const;
  AliMpIntPair& operator = (const AliMpIntPair& src) ;
  void operator += (const AliMpIntPair& op);
  void operator -= (const AliMpIntPair& op);

  // static get methods
  static AliMpIntPair Invalid() {return AliMpIntPair();}

  // get methods
  Int_t  GetFirst() const  {return fFirst;}
  Int_t  GetSecond() const {return fSecond;}
  Bool_t IsValid() const   {return fValidity;}

  // set methods
  void SetFirst(Int_t ix)  {fFirst=ix; fValidity=true; }
  void SetSecond(Int_t iy) {fSecond=iy; fValidity=true;}
  
 private:
  // data members
  Int_t   fFirst;    // position along x
  Int_t   fSecond;   // position along y
  Bool_t  fValidity; // validity

  ClassDef(AliMpIntPair,1) //utility class for the motif type
};

AliMpIntPair operator + (const AliMpIntPair& op1,const AliMpIntPair& op2);
AliMpIntPair operator - (const AliMpIntPair& op1,const AliMpIntPair& op2);
AliMpIntPair operator * (const AliMpIntPair& op1,const AliMpIntPair& op2);
ostream& operator << (ostream &stream,const AliMpIntPair& op);

#endif //ALI_MP_INT_PAIR_H
