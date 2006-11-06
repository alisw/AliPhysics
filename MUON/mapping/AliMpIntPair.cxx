/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

// $Id$
// $MpId: AliMpIntPair.cxx,v 1.7 2006/05/24 13:58:29 ivana Exp $
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

#include "AliMpIntPair.h"

#include "AliLog.h"

#include <Riostream.h>

/// \cond CLASSIMP
ClassImp(AliMpIntPair)
/// \endcond


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
    fValidity(true) 
{
/// Standard constructor
}

//_____________________________________________________________________________
AliMpIntPair::AliMpIntPair(Int_t ix,Int_t iy, Bool_t validity)
  : TObject(),
    fFirst(ix),
    fSecond(iy),
    fValidity(validity) 
{
/// Standard constructor with validity argument
}

//_____________________________________________________________________________
AliMpIntPair::AliMpIntPair()
  : TObject(),
    //fFirst(9999),
    //fSecond(9999),
    fFirst(0),
    fSecond(0),
    fValidity(false) 
{
/// Default constructor
}

//_____________________________________________________________________________
AliMpIntPair::AliMpIntPair(const AliMpIntPair& src):
  TObject(src),
  fFirst(src.fFirst),
  fSecond(src.fSecond),
  fValidity(src.fValidity)
{
/// Copy constructor
}

//_____________________________________________________________________________
AliMpIntPair::~AliMpIntPair() 
{
/// Destructor
}

//_____________________________________________________________________________
Bool_t AliMpIntPair::operator< (const AliMpIntPair& pos2) const
{
/// Less operator

  // fFirst prior to fSecond
  if (fFirst<pos2.fFirst) return kTRUE;
  if (fFirst>pos2.fFirst) return kFALSE;
  if (fSecond<pos2.fSecond) return kTRUE;
  return kFALSE;
}

//_____________________________________________________________________________
Bool_t AliMpIntPair::operator== (const AliMpIntPair& pos2) const
{
/// Equality operator

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
/// Non-equality operator

  // are this and pos2 equals?
  return !(*this == pos2);
}

//_____________________________________________________________________________
AliMpIntPair& AliMpIntPair::operator=(const AliMpIntPair& src) 
{
/// Assignment operator

  // check assignment to self
  if (this == &src) return *this;

  // base class assignment
  TObject::operator=(src);

  // assignment operator
  fFirst = src.fFirst;
  fSecond = src.fSecond;
  fValidity = src.fValidity;
  
  return *this;
}
//_____________________________________________________________________________
Int_t AliMpIntPair::Compare(const TObject* obj) const
{
/// Compare using operator <

  const AliMpIntPair* pair = dynamic_cast<const AliMpIntPair*>(obj);
  if ( !pair ) {
    AliErrorStream() << "Wrong object type." << endl;
    return -1;
  }  

  return ( *this < *pair ) ? -1 : 1;
}
//_____________________________________________________________________________
void AliMpIntPair::operator += (const AliMpIntPair& op)
{
/// Incrementation operator

  fFirst += op.fFirst;
  fSecond += op.fSecond;
  
  // operation only on valid pairs
  fValidity = fValidity && op.fValidity;
}
//_____________________________________________________________________________
void AliMpIntPair::operator -= (const AliMpIntPair& op)
{
/// Decrementation operator

  fFirst -= op.fFirst;
  fSecond -= op.fSecond;

  // operation only on valid pairs
  fValidity = fValidity && op.fValidity;
}

//_____________________________________________________________________________
AliMpIntPair operator-(const AliMpIntPair& op1,const AliMpIntPair& op2)
{
/// Substraction operator

  return AliMpIntPair(op1.GetFirst()-op2.GetFirst(),
                  op1.GetSecond()-op2.GetSecond(),
		  op1.IsValid() && op2.IsValid());
}
//_____________________________________________________________________________
AliMpIntPair operator+(const AliMpIntPair& op1,const AliMpIntPair& op2)
{
/// Addition operator

  return AliMpIntPair(op1.GetFirst()+op2.GetFirst(),
                  op1.GetSecond()+op2.GetSecond(),
		  op1.IsValid() && op2.IsValid());
}
//_____________________________________________________________________________
AliMpIntPair operator*(const AliMpIntPair& op1,const AliMpIntPair& op2)
{
/// Multiplication operator

  return AliMpIntPair(op1.GetFirst()*op2.GetFirst(),
                  op1.GetSecond()*op2.GetSecond(),
		  op1.IsValid() && op2.IsValid());
}
//_____________________________________________________________________________
ostream& operator<< (ostream &stream,const AliMpIntPair& op)
{
/// Output streaming

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

