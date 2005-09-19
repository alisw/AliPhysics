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
// $MpId: AliMpPad.cxx,v 1.6 2005/08/26 15:43:36 ivana Exp $
// Category: basic
//
// Class AliMpPad
// ---------------
// Class which encapsuate all informations about a pad
// Included in AliRoot: 2003/05/02
// Authors: David Guez, Ivana Hrivnacova; IPN Orsay

#include <Riostream.h>

#include "AliMpPad.h"

ClassImp(AliMpPad)

//////////////////////////////////////////////////////////
//
// This class encapsulate all the information about a pad
//
//////////////////////////////////////////////////////////

//
// foreign operators
//

//_____________________________________________________________________________
Bool_t operator==(const TVector2& v1,const TVector2& v2)
{
return v1.X()==v2.X() && v1.Y()==v2.Y();
}


//_____________________________________________________________________________
ostream& operator<<(ostream& out,const TVector2& v)
{
  out << '(' << v.X() << ',' << v.Y() << ')';
  return out; 
}

//_____________________________________________________________________________
AliMpPad::AliMpPad(const AliMpIntPair& location,const AliMpIntPair& indices,
                   const TVector2& position,const TVector2& dimensions,
                   Bool_t validity)
 : TObject(),
   fLocation(location),
   fIndices(indices),
   fPosition(position),
   fDimensions(dimensions),
   fValidity(validity)
{
/// Standard constructor                                                   \n
/// Be carefull : this constructor doesn't check the validity of
/// the correspondance between location and indices.
/// By default, validity is set true.
/// It is aimed to be used by MSegmentation methods, and never from outside....
}


//_____________________________________________________________________________
AliMpPad::AliMpPad()
  : TObject(),
    fLocation(AliMpIntPair::Invalid()),
    fIndices(AliMpIntPair::Invalid()),
    fPosition(-1.,-1.),
    fDimensions(0.,0.),
    fValidity(false) 
{
/// Default constructor - creates pad in invalid state
}


//_____________________________________________________________________________
AliMpPad::AliMpPad(const AliMpPad& src)
  : TObject(src)
{
/// Copy constructor

 *this = src;
}

//_____________________________________________________________________________
AliMpPad::~AliMpPad() 
{
/// Destructor
}

//_____________________________________________________________________________
AliMpPad& AliMpPad::operator = (const AliMpPad& src) 
{
/// Assignment operator
 
  // check assignment to self
  if (this == &src) return *this;

  // base class assignment
  TObject::operator=(src);

  // assignment operator
  fLocation   = src.fLocation;
  fIndices    = src.fIndices;
  fPosition.Set(src.fPosition);
  fDimensions.Set(src.fDimensions);
  fValidity = src.fValidity;

  return *this;
}

//_____________________________________________________________________________
Bool_t AliMpPad::operator == (const AliMpPad& pos2) const
{
/// Equality operator

  // are this and pos2 equals?

  // one valid, one invalid
  if (fValidity != pos2.fValidity) return false;
  
  // both invalid
  if (!fValidity) return true;
  
  // both valid
  return    (fLocation==pos2.fLocation) && (fIndices   ==pos2.fIndices   )
         && (fPosition==pos2.fPosition) && (fDimensions==pos2.fDimensions);
}
//_____________________________________________________________________________
Bool_t AliMpPad::operator!= (const AliMpPad& pos2) const
{
/// Non-equality operator

  // are this and pos2 equals?
  return !(*this==pos2);
}

//_____________________________________________________________________________
ostream& operator<< (ostream &out, const AliMpPad& op)
{
/// Output streaming

  if (op.IsValid()) {
    out << "Pad: Location " << op.GetLocation() 
        << "  Indices "     << op.GetIndices() 
	<< "  Position "    << op.Position()
        << "  Dimensions "  << op.Dimensions();
    return out;
  }
  else {
    out << "Pad::Invalid";
    return out;
  }  
}

//_____________________________________________________________________________
Bool_t operator < (const AliMpPad& left, const AliMpPad& right)
{
/// Less operator

  return left.GetIndices()<right.GetIndices();
}

//_____________________________________________________________________________
void AliMpPad::Print(const char* /*option*/) const
{
/// Prints all pad data.

  if (fValidity) {
    cout << "Indices: " << fIndices << "; "
         << " Location: " << fLocation << "; "
         << " Position: " << fPosition.X() << " " << fPosition.Y() << "; "
         << " Dimensions: " << fDimensions.X() << " " << fDimensions.Y() 
         << endl;
  }
  else {	 
    cout << "Pad::Invalid " << endl;
  }  
}

