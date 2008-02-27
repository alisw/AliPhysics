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

//-----------------------------------------------------------------------------
/// \class AliMUONVDigit
///
/// This is the base class of a MUON digit that most client code should deal with.
/// There should be no reason to have to use a concrete class in most cases.
///
/// All digits have basic features, like :
///
/// - a way to identify it : detection element, electronics card and
///   channel, cathode. Note that some static methods exists to compact 
///   those 4 informations into a single 4 bytes integer (stored in the
///   fUniqueID data member present in all TObjects). 
///
/// - its charge
///
/// - a set of boolean methods to indicate whether the digit has been calibrated, etc...
///
/// In addition, if HasMCInformation is true, the digit store also the list
/// of MC tracks that contributed to its charge
///
/// Also, if HasGeometryInformation is true, the digit knows the position and
/// the (half) dimensions (in cm) of the pad it corresponds to.
///
/// \author Laurent Aphecetche, Subatech
//-----------------------------------------------------------------------------

#include "AliMUONVDigit.h"

#include <Riostream.h>
#include <TClass.h>

/// \cond CLASSIMP
ClassImp(AliMUONVDigit)
/// \endcond

//_____________________________________________________________________________
AliMUONVDigit::AliMUONVDigit(Int_t detElemId, Int_t eCardId,
                             Int_t eCardChannel, Int_t cathode)
: TObject() 
{
  /// Normal constructor for trigger digits
  SetUniqueID(BuildUniqueID(detElemId,eCardId,eCardChannel,cathode));
}

//_____________________________________________________________________________
AliMUONVDigit::AliMUONVDigit()
: TObject() 
{
  /// Default ctor
}

//_____________________________________________________________________________
AliMUONVDigit::~AliMUONVDigit()
{
  /// dtor
}

//_____________________________________________________________________________
Bool_t 
AliMUONVDigit::IsEqual(const TObject* object) const
{
  /// Whether we're equal to object. 
  /// WARNING : only based on our identifiers (de,manu,channel,cathode), not our
  /// content (i.e. charge, status...)
  
  const AliMUONVDigit* d = static_cast<const AliMUONVDigit*>(object);
    
  return ( DetElemId() == d->DetElemId() &&
           Cathode() == d->Cathode() &&
           ManuId() == d->ManuId() &&
           ManuChannel() == d->ManuChannel() );
}

//_____________________________________________________________________________
Int_t 
AliMUONVDigit::Compare(const TObject* object) const
{
  /// Compare two digits, trying to get as complete an order as possible.
  /// We sort by DE, then by charge, then by manu, etc...
  ///
  const AliMUONVDigit* d = static_cast<const AliMUONVDigit*>(object);
  
  if ( DetElemId() > d->DetElemId() ) 
  {
    return 1;
  }
  else if ( DetElemId() < d->DetElemId() )
  {
    return -1;
  }
  else
  {
    if ( Charge() > d->Charge() )
    {
      return 1;
    }
    else if ( Charge() < d->Charge() )
    {
      return -1;
    }
    else
    {
      if ( ManuId() > d->ManuId() )
      {
        return 1;
      }
      else if ( ManuId() < d->ManuId() )
      {
        return -1;
      }
      else
      {
        if ( ManuChannel() > d->ManuChannel() )
        {
          return 1;
        }
        else if ( ManuChannel() < d->ManuChannel() )
        {
          return -1;
        }
      }
    }
  }
  return 0;
}

//_____________________________________________________________________________
UInt_t 
AliMUONVDigit::BuildUniqueID(Int_t detElemId, Int_t manuId, 
                             Int_t manuChannel, Int_t cathode)
{
  /// Build a single integer with id information
  return ( ( detElemId ) | ( manuId << 12 ) | ( manuChannel << 24 )
                | ( cathode << 30 ) );
}

//_____________________________________________________________________________
Int_t
AliMUONVDigit::DetElemId(UInt_t uniqueID)
{
  /// Return detection element id part of the uniqueID
  return uniqueID & 0xFFF;
}

//_____________________________________________________________________________
Int_t
AliMUONVDigit::ManuChannel(UInt_t uniqueID)
{
  /// Return manuChannel part of the uniqueID
  return ( uniqueID & 0x3F000000 ) >> 24;
}

//_____________________________________________________________________________
Int_t
AliMUONVDigit::ManuId(UInt_t uniqueID)
{
  /// Return manuId part of the uniqueID
  return ( uniqueID & 0xFFF000 ) >> 12;
}

//_____________________________________________________________________________
Int_t
AliMUONVDigit::Cathode(UInt_t uniqueID)
{
  /// Return the cathode part of the uniqueID
  return ( uniqueID & 0x40000000 ) >> 30;
}

//_____________________________________________________________________________
void
AliMUONVDigit::DecodeUniqueID(UInt_t uniqueID,
                              Int_t& detElemId, Int_t& manuId, 
                              Int_t& manuChannel, Int_t& cathode)
{
  /// Unpack uniqueID into 4 elements
  detElemId = DetElemId(uniqueID);
  manuId = ManuId(uniqueID);
  manuChannel = ManuChannel(uniqueID);
  cathode = Cathode(uniqueID);
}

//_____________________________________________________________________________
const char*
AliMUONVDigit::GetName() const
{
  /// Return the name of this digit, composed of its id parts.
  return Form("DE%04d-%04d-%02d-%d",
              DetElemId(),ManuId(),ManuChannel(),Cathode());
}

//_____________________________________________________________________________
void
AliMUONVDigit::Print(Option_t* opt) const
{
  /// Dump to screen.
  /// If opt=="tracks", info on tracks are printed too.
  
  cout << Form("<%s>: ID %12u DE %4d Cath %d (Ix,Iy)=(%3d,%3d) (Manu,Channel)=(%4d,%2d)"
               ", Charge=%7.2f",
               ClassName(),GetUniqueID(),
               DetElemId(),Cathode(),PadX(),PadY(),ManuId(),ManuChannel(),Charge());  
  if ( IsSaturated() ) 
  {
    cout << "(S)";
  }
  else
  {
    cout << "   ";
  }
  
  if ( IsCalibrated() )
  {
    cout << "(C)";
  }
  else
  {
    cout << "   ";
  }

  if ( IsUsed() )
  {
    cout << "(U)";
  }
  else
  {
    cout << "   ";
  }
  
  cout << Form(" ADC=%4d StatusMap=%04x",ADC(),StatusMap());

  TString options(opt);
  options.ToLower();
  if ( options.Contains("tracks") && HasMCInformation() )
  {
    cout << " Hit " << setw(3) << Hit();
    Int_t ntracks = Ntracks();
    if (ntracks) 
    {
      cout << " Tracks : " << setw(2) << ntracks;
      for ( Int_t i = 0; i < ntracks; ++i )
      {
        cout << " Track(" << i << ")=" << setw(3) << Track(i)
        << " Charge(" << i << ")=" << setw(5) << TrackCharge(i);
      }
    }
    else
    {
      cout << " no track info.";
    }
  }
  cout << endl;  
}
