// $Id$
// Category: basic
//
// Class AliMpIntPair
// --------------
// Class that defines the pair of integers.
// The pair created by the default constructor is in invalide state,
// setting one of values changes the state to valid.
//
// Included in AliRoot: 2003/05/02
// Authors: David Guez, Ivana Hrivnacova; IPN Orsay

#include <Riostream.h>

#include "AliMpIntPair.h"

ClassImp(AliMpIntPair)


///////////////////////////////////////////////////
//
// This class is a replacement for the standard STL
// pair<int,int> class, which can not be handed
// by the normal ROOT automatic streamer
// (at least in the ROOT version 3.03/03)
//
///////////////////////////////////////////////////


//_____________________________________________________________________________
AliMpIntPair::AliMpIntPair(Int_t ix,Int_t iy)
  : TObject(),
    fFirst(ix),
    fSecond(iy),
    fValidity(true) {
//
}

//_____________________________________________________________________________
AliMpIntPair::AliMpIntPair(Int_t ix,Int_t iy, Bool_t validity)
  : TObject(),
    fFirst(ix),
    fSecond(iy),
    fValidity(validity) {
//
}

//_____________________________________________________________________________
AliMpIntPair::AliMpIntPair()
  : TObject(),
    //fFirst(9999),
    //fSecond(9999),
    fFirst(0),
    fSecond(0),
    fValidity(false) {
//
}
//_____________________________________________________________________________
AliMpIntPair::AliMpIntPair(const AliMpIntPair& src):
  TObject(src),
  fFirst(src.fFirst),
  fSecond(src.fSecond),
  fValidity(src.fValidity)
{

}

//_____________________________________________________________________________
AliMpIntPair::~AliMpIntPair() {
//
}

//_____________________________________________________________________________
Bool_t AliMpIntPair::operator< (const AliMpIntPair& pos2) const
{
  // fFirst prior to fSecond
  if (fFirst<pos2.fFirst) return kTRUE;
  if (fFirst>pos2.fFirst) return kFALSE;
  if (fSecond<pos2.fSecond) return kTRUE;
  return kFALSE;
}

//_____________________________________________________________________________
Bool_t AliMpIntPair::operator== (const AliMpIntPair& pos2) const
{
  // are this and pos2 equals?
  
  // one valid, one invalid
  if (fValidity != pos2.fValidity) return false;
  
  // both invalid
  if (!fValidity) return true;
  
  // both valid
  return (fFirst==pos2.fFirst) && (fSecond==pos2.fSecond);
}

//_____________________________________________________________________________
Bool_t AliMpIntPair::operator!= (const AliMpIntPair& pos2) const
{
  // are this and pos2 equals?
  return !(*this == pos2);
}

//_____________________________________________________________________________
AliMpIntPair& AliMpIntPair::operator=(const AliMpIntPair& src) 
{
  // check assignement to self
  if (this == &src) return *this;

  // base class assignement
  TObject::operator=(src);

  // assignement operator
  fFirst = src.fFirst;
  fSecond = src.fSecond;
  fValidity = src.fValidity;
  
  return *this;
}

//_____________________________________________________________________________
void AliMpIntPair::operator += (const AliMpIntPair& op)
{
  // incrementation operator
  fFirst += op.fFirst;
  fSecond += op.fSecond;
  
  // operation only on valid pairs
  fValidity = fValidity && op.fValidity;
}
//_____________________________________________________________________________
void AliMpIntPair::operator -= (const AliMpIntPair& op)
{
  // decrementation operator
  fFirst -= op.fFirst;
  fSecond -= op.fSecond;

  // operation only on valid pairs
  fValidity = fValidity && op.fValidity;
}

//_____________________________________________________________________________
AliMpIntPair operator-(const AliMpIntPair& op1,const AliMpIntPair& op2)
{
  return AliMpIntPair(op1.GetFirst()-op2.GetFirst(),
                  op1.GetSecond()-op2.GetSecond(),
		  op1.IsValid() && op2.IsValid());
}
//_____________________________________________________________________________
AliMpIntPair operator+(const AliMpIntPair& op1,const AliMpIntPair& op2)
{
  return AliMpIntPair(op1.GetFirst()+op2.GetFirst(),
                  op1.GetSecond()+op2.GetSecond(),
		  op1.IsValid() && op2.IsValid());
}
//_____________________________________________________________________________
AliMpIntPair operator*(const AliMpIntPair& op1,const AliMpIntPair& op2)
{
  return AliMpIntPair(op1.GetFirst()*op2.GetFirst(),
                  op1.GetSecond()*op2.GetSecond(),
		  op1.IsValid() && op2.IsValid());
}
//_____________________________________________________________________________
ostream& operator<< (ostream &stream,const AliMpIntPair& op)
{
  if (op.IsValid()) {
    stream << '(';
    stream << op.GetFirst()<<','<<op.GetSecond()<<')';
    return stream;
  }  
  else { 
    stream << "AliMpIntPair::Invalid";
    return stream;
  }  
}

