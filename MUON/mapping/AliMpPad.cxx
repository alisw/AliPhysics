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
// $MpId: AliMpPad.cxx,v 1.9 2006/05/24 13:58:29 ivana Exp $
// Category: basic

//-----------------------------------------------------------------------------
// Class AliMpPad
// ---------------
// Class which encapsuate all informations about a pad
// Included in AliRoot: 2003/05/02
// Authors: David Guez, Ivana Hrivnacova; IPN Orsay
// root [0] .x testSectorAreaIterator.C
// Real time 0:00:56, CP time 36.270
//-----------------------------------------------------------------------------

#include "AliMpPad.h"
#include "AliLog.h"

#include <TClonesArray.h>
#include <Riostream.h>

/// \cond CLASSIMP
ClassImp(AliMpPad)
/// \endcond

const Int_t  AliMpPad::fgkMaxNofLocations = 6;

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
   fLocations(0),
   fNofLocations(0),
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
    fLocations(0),
    fNofLocations(0),
    fLocation(AliMpIntPair::Invalid()),
    fIndices(AliMpIntPair::Invalid()),
    fPosition(-1.,-1.),
    fDimensions(0.,0.),
    fValidity(false) 
{
/// Default constructor - creates pad in invalid state
}


//_____________________________________________________________________________
AliMpPad::AliMpPad(const AliMpPad& rhs)
  : TObject(rhs),
    fLocations(0),
    fNofLocations(0),
    fLocation(AliMpIntPair::Invalid()),
    fIndices(AliMpIntPair::Invalid()),
    fPosition(-1.,-1.),
    fDimensions(0.,0.),
    fValidity(false) 
{
/// Copy constructor

 *this = rhs;
}

//_____________________________________________________________________________
AliMpPad::~AliMpPad() 
{
/// Destructor

  delete [] fLocations;
}

//_____________________________________________________________________________
AliMpPad& AliMpPad::operator = (const AliMpPad& rhs) 
{
/// Assignment operator
 
  // check assignment to self
  if (this == &rhs) return *this;

  // base class assignment
  TObject::operator=(rhs);

  // assignment operator
  fLocation   = rhs.fLocation;
  fIndices    = rhs.fIndices;
  fPosition.Set(rhs.fPosition);
  fDimensions.Set(rhs.fDimensions);
  fValidity = rhs.fValidity;
  
  fLocations = 0;
  fNofLocations = rhs.fNofLocations;
  if ( rhs.GetNofLocations() ) {
    fLocations = new AliMpIntPair[fgkMaxNofLocations];
    for ( UInt_t i=0; i<rhs.fNofLocations; i++ )
      fLocations[i] = rhs.fLocations[i];
  }  			

  return *this;
}

//_____________________________________________________________________________
Bool_t AliMpPad::operator == (const AliMpPad& rhs) const
{
/// Equality operator

  // are this and rhs equals?

  // one valid, one invalid
  if (fValidity != rhs.fValidity) return false;
  
  // both invalid
  if (!fValidity) return true;
  
  // both valid
  Bool_t sameLocations = true;
  
  if (rhs.GetNofLocations()) {
    for (Int_t i=0; i<rhs.GetNofLocations(); i++) 
      if ( GetLocation(i) != rhs.GetLocation(i) )
        sameLocations = false;
  }
  
  return    (fLocation   == rhs.fLocation) 
         && (fIndices    == rhs.fIndices)
         && (fPosition   == rhs.fPosition) 
	 && (fDimensions == rhs.fDimensions)
	 && sameLocations;
}
//_____________________________________________________________________________
Bool_t AliMpPad::operator != (const AliMpPad& rhs) const
{
/// Non-equality operator

  // are this and rhs equals?
  return !(*this==rhs);
}

//_____________________________________________________________________________
Bool_t operator < (const AliMpPad& left, const AliMpPad& right)
{
/// Less operator

  return left.GetIndices()<right.GetIndices();
}

//_____________________________________________________________________________
Bool_t AliMpPad::AddLocation(const AliMpIntPair& location, Bool_t warn)
{
/// Add location to the collection if not yet present and
/// if collection is not yet full                                           \n
/// Return false and optionally give a warning if location is not 
/// added. 

  // Check maximum number limit
  if ( GetNofLocations() == fgkMaxNofLocations ) {
    if (warn) {
      AliWarningStream() << "Cannot add location: "
                         << location
			 << "  Maximum number has been reached." << endl;
    }
    return false;
  }  			 

  // Check if location is present
  if ( HasLocation(location) ) {
    if (warn) {
      AliWarningStream() << "Cannot add location: "
                         << location
			 << "  Location is already present." << endl;
    }
    return false;
  } 
  
  // Add location
  if ( ! fLocations)
    fLocations = new AliMpIntPair[fgkMaxNofLocations];
  
  fLocations[fNofLocations++] = location;
  return true;
}

//_____________________________________________________________________________
void AliMpPad::PrintOn(ostream& out) const
{
/// Prints all pad data.

  if ( !fValidity ) {
    out << "Pad::Invalid";
    return;
  }  

  out << "Pad: Location " << fLocation 
      << "  Indices "     << fIndices
      << "  Position "    << fPosition
      << "  Dimensions "  << fDimensions;

  if ( GetNofLocations() ) {
    out << endl;
    out << "     Other locations: ";

    for (Int_t i=0; i<GetNofLocations(); i++) 
        out << GetLocation(i) << "  ";
  }
}

//_____________________________________________________________________________
void AliMpPad::Print(const char* /*option*/) const
{
/// Prints all pad data.

  PrintOn(cout);
  cout << endl;
}

//_____________________________________________________________________________
Int_t  AliMpPad::GetNofLocations() const
{
/// Return number of other locations associated with this pad

  if (!fLocations) return 0;
  
  return fNofLocations;
}  
  

//_____________________________________________________________________________
AliMpIntPair AliMpPad::GetLocation(Int_t i) const
{
/// Return i-th other location associated with this pad

  if ( !fLocations || i<0 || i>=GetNofLocations() ) 
    return AliMpIntPair::Invalid();

  return fLocations[i];
}  

//_____________________________________________________________________________
Bool_t AliMpPad::HasLocation(const AliMpIntPair& location) const
{
/// Return true if given location is present either as fLocation
/// or in the collectio

  if (fLocation == location) return true;

  for ( Int_t i=0; i<GetNofLocations(); i++ ) {
    if ( GetLocation(i) == location ) return true;
  }
    
  return false;
}      

//_____________________________________________________________________________
ostream& operator<< (ostream &out, const AliMpPad& pad)
{
/// Output streaming

  pad.PrintOn(out);

  return out;
}

