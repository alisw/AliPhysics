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
#include "AliMpEncodePair.h"
#include "AliLog.h"

#include <TClonesArray.h>
#include <Riostream.h>

/// \cond CLASSIMP
ClassImp(AliMpPad)
/// \endcond

const Int_t  AliMpPad::fgkMaxNofLocations = 6;

//_____________________________________________________________________________
AliMpPad::AliMpPad(Int_t manuId, Int_t channel,
                   Int_t ix, Int_t iy,
                   Double_t x,  Double_t y, 
                   Double_t dx,  Double_t dy, 
                   Bool_t validity)
 : TObject(),
   fNofLocations(0),
   fLLocations(0),
   fLLocation(AliMp::Pair(manuId, channel)),
   fLIndices(AliMp::Pair(ix, iy)),
   fPositionX(x),
   fPositionY(y),
   fDimensionX(dx),
   fDimensionY(dy),
   fValidity(validity)
{
/// Standard constructor                                                   \n
/// Be carefull : this constructor doesn't check the validity of
/// the correspondance between location and indices.
/// By default, validity is set true.
/// It is aimed to be used by MSegmentation methods, and never from outside....
}

//_____________________________________________________________________________
AliMpPad::AliMpPad(Int_t manuId, Int_t channel,
                   MpPair_t indices,
                   Double_t x,  Double_t y, 
                   Double_t dx,  Double_t dy, 
                   Bool_t validity)
 : TObject(),
   fNofLocations(0),
   fLLocations(0),
   fLLocation(AliMp::Pair(manuId, channel)),
   fLIndices(indices),
   fPositionX(x),
   fPositionY(y),
   fDimensionX(dx),
   fDimensionY(dy),
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
    fNofLocations(0),
    fLLocations(0),
    fLLocation(0),
    fLIndices(0),
    fPositionX(-1.),
    fPositionY(-1.),
    fDimensionX(-1.),
    fDimensionY(-1.),
    fValidity(false) 
{
/// Default constructor - creates pad in invalid state
}

//_____________________________________________________________________________
AliMpPad::AliMpPad(const AliMpPad& rhs)
  : TObject(),
    fNofLocations(0),
    fLLocations(0),
    fLLocation(0),
    fLIndices(0),
    fPositionX(-1.),
    fPositionY(-1.),
    fDimensionX(-1.),
    fDimensionY(-1.),
    fValidity(false) 
{
/// Copy constructor

 *this = rhs;
}

//_____________________________________________________________________________
AliMpPad::~AliMpPad() 
{
/// Destructor

  delete [] fLLocations;
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
  fLLocation   = rhs.fLLocation;
  fLIndices    = rhs.fLIndices;
  fPositionX   = rhs.fPositionX;
  fPositionY   = rhs.fPositionY;
  fDimensionX  = rhs.fDimensionX;
  fDimensionY  = rhs.fDimensionY;
  fValidity = rhs.fValidity;
  
  delete [] fLLocations;
  fLLocations = 0;
  fNofLocations = rhs.fNofLocations;
  if ( rhs.GetNofLocations() ) {
    fLLocations = new MpPair_t[fgkMaxNofLocations];
    for ( UInt_t i=0; i<rhs.fNofLocations; i++ )
      fLLocations[i] = rhs.fLLocations[i];
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
  
  return    ( fLLocation  == rhs.fLLocation ) 
         && ( fLIndices   == rhs.fLIndices )
         && ( fPositionX  == rhs.fPositionX ) 
         && ( fPositionY  == rhs.fPositionY ) 
	 && ( fDimensionX == rhs.fDimensionX )
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

  if ( left.GetIx() < right.GetIx() ) return kTRUE;
  if ( left.GetIx() > right.GetIx() ) return kFALSE;
  if ( left.GetIy() < right.GetIy() ) return kTRUE;
  return kFALSE;
}

//_____________________________________________________________________________
Bool_t AliMpPad::AddLocation(Int_t localBoardId, Int_t localBoardChannel, 
                             Bool_t warn)
{
/// Add location to the collection if not yet present and
/// if collection is not yet full                                           \n
/// Return false and optionally give a warning if location is not 
/// added. 

  // Check maximum number limit
  if ( GetNofLocations() == fgkMaxNofLocations ) {
    if (warn) {
      AliWarningStream() << "Cannot add location: ("
                         << localBoardId << "," << localBoardChannel << ")."
       	                 << "  Maximum number has been reached." << endl;
    }
    return false;
  }  			 

  // Check if location is present
  if ( HasLocation(localBoardId, localBoardChannel) ) {
    if (warn) {
      AliWarningStream() << "Cannot add location: "
                         << localBoardId << "," << localBoardChannel << ")."
                         << "  Location is already present." << endl;
    }
    return false;
  } 
  
  // Add location
  if ( ! fLLocations)
    fLLocations = new MpPair_t[fgkMaxNofLocations];
  
  fLLocations[fNofLocations++] 
    = AliMp::Pair(localBoardId, localBoardChannel);

  return true;
}

//_____________________________________________________________________________
Int_t  AliMpPad::GetManuId() const
{
/// Return pad manu Id 

  return AliMp::PairFirst(fLLocation);
}  

//_____________________________________________________________________________
Int_t  AliMpPad::GetManuChannel() const
{
/// Return pad manu channel

  return AliMp::PairSecond(fLLocation);
}  

//_____________________________________________________________________________
Int_t  AliMpPad::GetIx() const
{
/// Return pad index ix

  return AliMp::PairFirst(fLIndices);
}  

//_____________________________________________________________________________
Int_t  AliMpPad::GetIy() const
{
/// Return pad index iy

  return AliMp::PairSecond(fLIndices);
}  

//_____________________________________________________________________________
void AliMpPad::PrintOn(ostream& out) const
{
/// Prints all pad data.

  if ( !fValidity ) {
    out << "Pad::Invalid";
    return;
  }  

  out << "Pad: Location ";
  AliMp::PairPut(out, fLLocation)
      << "  Indices ";     
  AliMp::PairPut(out,fLIndices)
      << "  Position "
      << "(" << fPositionX << "," << fPositionY << ")"
      << "  Dimensions "
      << "(" << fDimensionX << "," << fDimensionY << ")";

  if ( GetNofLocations() ) {
    out << endl;
    out << "     Other locations: ";

    for (Int_t i=0; i<GetNofLocations(); i++) 
        AliMp::PairPut(out,GetLocation(i)) << "  ";
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

  if (!fLLocations) return 0;
  
  return fNofLocations;
}  
  

//_____________________________________________________________________________
MpPair_t AliMpPad::GetLocation(Int_t i) const
{
/// Return i-th other location associated with this pad

  if ( !fLLocations || i<0 || i>=GetNofLocations() ) 
    return 0;

  return fLLocations[i];
}  

//_____________________________________________________________________________
Int_t AliMpPad::GetLocalBoardId(Int_t i) const
{
/// Return i-th other local board Id associated with this pad

  if ( !fLLocations || i<0 || i>=GetNofLocations() ) 
    return 0;

  return AliMp::PairFirst(fLLocations[i]);
}  

//_____________________________________________________________________________
Int_t AliMpPad::GetLocalBoardChannel(Int_t i) const
{
/// Return i-th other local board channel associated with this pad

  if ( !fLLocations || i<0 || i>=GetNofLocations() ) 
    return 0;

  return AliMp::PairSecond(fLLocations[i]);
}  

//_____________________________________________________________________________
Bool_t AliMpPad::HasLocation(Int_t localBoardId, Int_t localBoardChannel) const
{
/// Return true if given location is present either as fLLocation
/// or in the collectio

  MpPair_t location = AliMp::Pair(localBoardId, localBoardChannel);

  if (fLLocation == location) return true;

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

