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
* about the suitability of this software for any purpeateose. It is      *
* provided "as is" without express or implied warranty.                  *
**************************************************************************/

// $Id$
// $MpId: AliMpTrigger.cxx,v 1.4 2006/05/24 13:58:52 ivana Exp $

#include "AliMpTrigger.h"
#include "AliMpSlatSegmentation.h"
#include "AliMpConstants.h"
#include "AliLog.h"
#include "AliMpSlat.h"

#include "Riostream.h"
#include "TArrayI.h"
#include "TObjArray.h"

/// 
/// \class AliMpTrigger
/// 
/// A trigger 'slat' object. 
/// It is to be viewed as a superposition of  
/// virtual layers of AliMpSlat objects. The need for more than one layer  
/// arise from the fact that a given local board deals with strips  
/// located in different detelem. So a given strip (pad) can have several  
/// "locations".
///
/// \author Laurent Aphecetche

/// \cond CLASSIMP
ClassImp(AliMpTrigger)
/// \endcond

namespace
{
  Bool_t IsEqual(Double_t a, Double_t b, Double_t precision)
{
    if (b)
    {
      Double_t diff = TMath::Abs(b-a)/TMath::Abs(b);
      if ( diff < precision ) 
      {
        return kTRUE;
      }
    }
    else
    {
      if ( !a ) return kTRUE;
    }
    return kFALSE;
}
}

//_____________________________________________________________________________
AliMpTrigger::AliMpTrigger()
: TObject(), 
  fId(""), 
  fPlaneType(kNonBendingPlane), 
  fSlats(0),
  fSlatSegmentations(0),
  fMaxNofPadsY(0),
  fDX(0), 
  fDY(0)
{
  // default ctor

  AliDebugStream(1) << "this = " << this << endl;

  fSlats.SetOwner(kFALSE);
  fSlatSegmentations.SetOwner(kTRUE);
}

//_____________________________________________________________________________
AliMpTrigger::AliMpTrigger(const char* slatType, AliMpPlaneType bendingOrNot)
    :  TObject(), 
       fId(slatType), 
       fPlaneType(bendingOrNot), 
       fSlats(0),
       fSlatSegmentations(0),
       fMaxNofPadsY(0), 
       fDX(0), 
       fDY(0)
{
  // normal ctor

  AliDebugStream(1) << "this = " << this << endl;

  fSlats.SetOwner(kFALSE);
  fSlatSegmentations.SetOwner(kTRUE);
}

//_____________________________________________________________________________
AliMpTrigger::~AliMpTrigger()
{
  // dtor
  AliDebugStream(1) << "this = " << this << endl;

  fSlatSegmentations.Delete();
}

//_____________________________________________________________________________
Bool_t
AliMpTrigger::AdoptLayer(AliMpSlat* slat)
{
  // Adopt (i.e. we become owner of that pointer) a slat, as 
  // a layer of this trigger slat.

  AliDebug(2,Form("%s is adopting %s ",
                  GetID(),slat->GetID()));

  // Check that we keep our size constant.
    
  if ( GetSize() > 0 && 
       ( !::IsEqual(slat->DX(),fDX,AliMpConstants::LengthTolerance()) || 
         !::IsEqual(slat->DY(),fDY,AliMpConstants::LengthTolerance()) )
     )
  {
    AliError(Form("In %s trying to add a layer (%e,%e) (layer #%d) "
                  "of a different size than mine (%e,%e)",
                  GetID(),slat->DX(),slat->DY(),fSlats.GetEntries(),
                  fDX,fDY));
    return kFALSE;
  }
  fSlats.Add(slat);
  Bool_t owner(kTRUE);
  // the slat segmentation will be the owner of the slat, and will delete
  // it when it'll be deleted itself
  fSlatSegmentations.Add(new AliMpSlatSegmentation(slat,owner));
  fMaxNofPadsY = std::max(slat->GetMaxNofPadsY(),fMaxNofPadsY);
  fDX = std::max(fDX,slat->DX());
  fDY = std::max(fDY,slat->DY());
  return kTRUE;
}

//_____________________________________________________________________________
TVector2
AliMpTrigger::Dimensions() const
{
  // Returns the dimensions (half-sizes) of that slat (cm)
  return TVector2(DX(),DY());
}

//_____________________________________________________________________________
Double_t
AliMpTrigger::DX() const
{
  // Returns the half-size in X (cm)
  return fDX;
}

//_____________________________________________________________________________
Double_t
AliMpTrigger::DY() const
{
  // Returns the half-size in Y (cm)
  return fDY;
}

//_____________________________________________________________________________
void 
AliMpTrigger::GetAllLocalBoardNumbers(TArrayI& lbn) const
{
  // Fills lbn with the local board numbers we're dealing with
  Int_t n(0);
  for ( Int_t i = 0; i < GetSize(); ++i )
  {
    n += GetLayer(i)->GetNofElectronicCards();
  }
  
  lbn.Set(n);

  Int_t index(0);
  
  for ( Int_t i = 0; i < GetSize(); ++i )
  {
    TArrayI slbn;
    GetLayer(i)->GetAllMotifPositionsIDs(slbn);
    for ( Int_t j = 0; j < slbn.GetSize(); ++j )
    {
      lbn[index] = slbn[j];
      ++index;
    }
  }
}

//_____________________________________________________________________________
const char*
AliMpTrigger::GetID() const
{
  // returns the id of this slat
  return fId.Data();
}

//_____________________________________________________________________________
const char*
AliMpTrigger::GetName() const
{
  // returns the name (=id+bending/non-bending) of this slat
  TString name(GetID());
  if ( fPlaneType == kBendingPlane )
  {
    name += ".Bending";
  }
  else if ( fPlaneType == kNonBendingPlane )
  {
    name += ".NonBending";
  }
  else
  {
    name += ".Invalid";
  }
  return name.Data();
}

//_____________________________________________________________________________
AliMpSlat*
AliMpTrigger::GetLayer(int layer) const
{
  // Returns a given layer
  if ( IsLayerValid(layer) )
  {
    return (AliMpSlat*)fSlats.At(layer);
  }
  return 0;
}

//_____________________________________________________________________________
AliMpVSegmentation*
AliMpTrigger::GetLayerSegmentation(int layer) const
{
  // Returns a given layer
  if ( IsLayerValid(layer) )
  {
    return (AliMpSlatSegmentation*)fSlatSegmentations.At(layer);
  }
  return 0;
}

//_____________________________________________________________________________
Int_t
AliMpTrigger::GetNofPadsX() const
{
  // Returns the number of pad in x direction
  if ( !GetSize() ) return -1;
  if ( GetLayer(0) )
  {
    return GetLayer(0)->GetNofPadsX();
  }
  return -1;
}

//_____________________________________________________________________________
Int_t
AliMpTrigger::GetMaxNofPadsY() const
{
  // Maximum number of pads in y direction
  return fMaxNofPadsY;
}

//_____________________________________________________________________________
Int_t
AliMpTrigger::GetSize() const
{
  // Number of layers
  return fSlats.GetEntriesFast();
}

//_____________________________________________________________________________
Bool_t
AliMpTrigger::IsLayerValid(int layer) const
{
  // Whether a given layer index is valid or not
  if ( layer >= 0 && layer < GetSize() )
  {
    return kTRUE;
  }
  return kFALSE;
}

//_____________________________________________________________________________
AliMpPlaneType
AliMpTrigger::PlaneType() const
{
  // Bending or not
  return fPlaneType;
}

//_____________________________________________________________________________
TVector2
AliMpTrigger::Position() const
{
  // Slat position (cm)
  return TVector2(DX(),DY());
}

//_____________________________________________________________________________
void
AliMpTrigger::Print(Option_t* opt) const
{
  // Dump on screen
  cout << "AliMpTrigger::" << GetID();
  if ( GetSize() == 0 )
  {
    cout << " Empty";
  }
  else if ( GetSize() > 1 )
  {
    cout << " Number of layers : " << GetSize();
  }
  else 
  {
    cout << " One layer";
  }
  cout << endl;
  for ( Int_t i = 0; i < GetSize(); ++i ) 
  {
    cout << "   ";
    GetLayer(i)->Print(opt);
  }
}

//_____________________________________________________________________________
//_____________________________________________________________________________
//_____________________________________________________________________________
//_____________________________________________________________________________
