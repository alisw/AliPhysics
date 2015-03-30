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

#include "AliMUONPainterPadStore.h"

#include "AliMUONCalibParamND.h"
#include "AliMUON2DMap.h"
#include "AliMUONVStore.h"
#include "AliMUONVDigit.h"
#include "AliLog.h"
#include <Riostream.h>
#include <TArrayI.h>
#include <TVector2.h>

///\class AliMUONPainterPadStore
///
/// Container for pads
///
///\author Laurent Aphecetche, Subatech

using std::cout;
using std::endl;
///\cond CLASSIMP
ClassImp(AliMUONPainterPadStore)
///\endcond

//_____________________________________________________________________________
AliMUONPainterPadStore::AliMUONPainterPadStore(TRootIOCtor* /*dummy*/) : TObject(),
fPadStore(0x0)
{
  /// ctor
}

//_____________________________________________________________________________
AliMUONPainterPadStore::AliMUONPainterPadStore() : TObject(),
  fPadStore(new AliMUON2DMap(kTRUE))
{
    /// ctor
}

//_____________________________________________________________________________
AliMUONPainterPadStore::~AliMUONPainterPadStore()
{
  /// dtor
  delete fPadStore;
}

//_____________________________________________________________________________
Int_t
AliMUONPainterPadStore::FindPadID(const TArrayI& pads, Double_t x, Double_t y) const
{
  /// Find, in array of pads, the one which contains (x,y). Returns -1 if not
  /// found
  
  for ( Int_t i = 0; i < pads.GetSize(); ++i ) 
  {
    Int_t id = pads.At(i);
    
    TVector2 position;
    TVector2 dimensions;
    
    GetPadGeometry(id,position,dimensions);
    
    TVector2 bl(position-dimensions);
    TVector2 ur(position+dimensions);    
    if ( bl.X() <= x && ur.X() >= x && bl.Y() <= y && ur.Y() >= y ) 
    {
      return id;
    }
  }
  return -1;
}


//_____________________________________________________________________________
AliMUONVCalibParam*
AliMUONPainterPadStore::Get(Int_t detElemId, Int_t manuId) const
{
  /// Get the pad container for a given manu
  
  AliMUONVCalibParam* param = 
  static_cast<AliMUONVCalibParam*>(fPadStore->FindObject(detElemId,manuId));
  
  if (!param)
  {
    param = new AliMUONCalibParamND(4,64,detElemId,manuId,-1.0);
    fPadStore->Add(param);
  }
  
  return param;
}

//_____________________________________________________________________________
void
AliMUONPainterPadStore::GetBoundaries(const TArrayI& pads,
                                      Double_t& xmin,
                                      Double_t& ymin,
                                      Double_t& xmax,
                                      Double_t& ymax) const
{
  /// Get the area covered by an array of pads
  
  xmin=ymin=1E9;
  xmax=ymax=-1E9;
  
  for ( Int_t i = 0; i < pads.GetSize(); ++i ) 
  {
    Int_t id = pads.At(i);
    
    TVector2 position;
    TVector2 dimensions;
    
    GetPadGeometry(id,position,dimensions);
    
    TVector2 bl(position-dimensions);
    TVector2 ur(position+dimensions);
    xmin = TMath::Min(xmin,bl.X());
    ymin = TMath::Min(ymin,bl.Y());
    xmax = TMath::Max(xmax,ur.X());
    ymax = TMath::Max(ymax,ur.Y());
  }     
}

//_____________________________________________________________________________
void
AliMUONPainterPadStore::GetPadGeometry(Int_t padId, 
                                       TVector2& position,
                                       TVector2& dimensions) const
{
  /// Get the geomtry of one pad
  
  if ( padId < 0 ) 
  {
    AliError(Form("padId is < 0 : %d",padId));
    position.Set(0.0,0.0);
    dimensions.Set(-1.0,-1.0);
    return;
  }
  
  Int_t detElemId = AliMUONVDigit::DetElemId(padId);
  Int_t manuId = AliMUONVDigit::ManuId(padId);
  Int_t manuChannel = AliMUONVDigit::ManuChannel(padId);
  
  AliMUONVCalibParam* param = 
    static_cast<AliMUONVCalibParam*>(fPadStore->FindObject(detElemId,manuId));
  
  if (!param)
  {
    AliError(Form("Could not find object DE %d manu %d",detElemId,manuId));
    position.Set(0.0,0.0);
    dimensions.Set(-1.0,-1.0);
    return;
  }
  
  position.Set(param->ValueAsDouble(manuChannel,0),
               param->ValueAsDouble(manuChannel,1));
  
  dimensions.Set(param->ValueAsDouble(manuChannel,2),
                 param->ValueAsDouble(manuChannel,3));
  
}

//_____________________________________________________________________________
Int_t
AliMUONPainterPadStore::GetSize() const
{
  /// Get the number of pads we handle
  
  TIter next(fPadStore->CreateIterator());
  AliMUONVCalibParam* param;
  Int_t n(0);
  
  while ( ( param = static_cast<AliMUONVCalibParam*>(next()) ) )
  {
    for ( Int_t i = 0; i < param->Size(); ++i ) 
    {
      if ( param->ValueAsDouble(i,2) >= 0 && param->ValueAsDouble(i,3) >= 0 ) 
      {
        ++n;
      }
    }
  }
  
  return n;
}

//_____________________________________________________________________________
void
AliMUONPainterPadStore::PrintPads(const TArrayI& pads) const
{
  /// Printout
  cout << "n=" << pads.GetSize() << endl;
  
  for ( Int_t i = 0; i < pads.GetSize(); ++i ) 
  {
    Int_t id = pads.At(i);
    TVector2 position, dimensions;
    GetPadGeometry(id,position,dimensions);
    cout << Form("i %4d DE %4d ManuID %4d ManuChannel %2d (X,Y)=(%7.3f,%7.3f)"
                 " (DX,DY)=(%7.3f,%7.3f)",
                 i,
                 AliMUONVDigit::DetElemId(id),
                 AliMUONVDigit::ManuId(id),
                 AliMUONVDigit::ManuChannel(id),
                 position.X(),position.Y(),
                 dimensions.X(),dimensions.Y()) << endl;
  }
}

