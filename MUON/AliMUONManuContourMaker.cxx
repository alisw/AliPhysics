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

/// \class AliMUONManuContourMaker
///
/// Maker of manu contours. 
///
/// Make use of the AliMUONContourMaker class, but this one contains
/// specific things for MUON (as the mapping, for instance), hence its
/// separation from AliMUONContourMaker.
///
/// This class requires that the mapping is loaded before anything can be done.
///
/// \author Laurent Aphecetche, Subatech

#include "AliMUONManuContourMaker.h"

#include "AliCodeTimer.h"
#include "AliLog.h"
#include "AliMUONContour.h"
#include "AliMUONContourMaker.h"
#include "AliMUONPolygon.h"
#include "AliMpCathodType.h"
#include "AliMpConnection.h"
#include "AliMpConstants.h"
#include "AliMpDDLStore.h"
#include "AliMpDEManager.h"
#include "AliMpIntPair.h"
#include "AliMpMotifPosition.h"
#include "AliMpMotifType.h"
#include "AliMpManuIterator.h"
#include "AliMpPlaneType.h"
#include "AliMpSegmentation.h"
#include "AliMpUID.h"
#include "AliMpVMotif.h"
#include "AliMpVSegmentation.h"
#include "TGeoMatrix.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TString.h"
#include "TVector2.h"

///\cond CLASSIMP
ClassImp(AliMUONManuContourMaker)
///\endcond

//_____________________________________________________________________________
AliMUONManuContourMaker::AliMUONManuContourMaker(AliMpExMap* deTransformations)
: TObject(), fDETransformations(deTransformations), fLocalManuContours(222,1)
{
/// Standard constructor

  fLocalManuContours.SetOwnerKeyValue(kTRUE,kTRUE);  
}

//_____________________________________________________________________________
AliMUONManuContourMaker::~AliMUONManuContourMaker()
{
/// Destructor
}

//_____________________________________________________________________________
AliMUONContour* 
AliMUONManuContourMaker::CreateManuContour(Int_t detElemId, Int_t manuId, const char* name) const
{
  /// Create the contour of a given manu (global coordinates)
 
  AliCodeTimerAuto("");
  
  TString sname(name);
  
  if ( sname.Length()==0 )
  {
    sname = ManuPathName(detElemId,manuId);
  }
  
  const AliMpVSegmentation* seg = AliMpSegmentation::Instance()->GetMpSegmentationByElectronics(detElemId,manuId);
  const AliMpMotifPosition* motifPos = seg->MotifPosition(manuId);

  AliMUONContour* contour = CreateMotifContour(*motifPos);
  
  if (!contour)
  {
    AliError(Form("Could not build contour %s",sname.Data()));
    return 0x0;
  }
  
  contour->SetName(sname.Data());
  
  contour->Offset(motifPos->GetPositionX()-seg->GetPositionX(),
                  motifPos->GetPositionY()-seg->GetPositionY());
  
  TGeoHMatrix* matrix = 0x0;
  
  if ( fDETransformations ) 
  {
    matrix = static_cast<TGeoHMatrix*>(fDETransformations->GetValue(detElemId));
    if ( matrix ) contour->Transform(*matrix);
  }
  
  return contour;
}


//_____________________________________________________________________________
AliMUONContour* 
AliMUONManuContourMaker::CreateMotifContour(const AliMpMotifPosition& motifPosition) const
{
  /// Create the contour of a given MOTIF (i.e. local coordinates only).
  
  AliCodeTimerAuto("");
  
  TString mpName(NameIt(motifPosition));
  
  AliMUONContour* contour = static_cast<AliMUONContour*>(fLocalManuContours.GetValue(mpName.Data()));
  
  if ( contour ) 
  {
    // if we have already done the job, just have to clone it and we are done
    return static_cast<AliMUONContour*>(contour->Clone());
  }
  
  TObjArray polygons(AliMpConstants::ManuNofChannels()); // array of AliMUONPolygon objects
  polygons.SetOwner(kTRUE);
  
  AliMpVMotif* motif = motifPosition.GetMotif();
  
  AliMpMotifType* motifType = motif->GetMotifType();
  
  if ( motifType->IsFull() ) 
  {
    // motif is a simple rectangle. No need to loop over pads, we can
    // compute the contour right here and now.
    polygons.Add(new AliMUONPolygon(0.0,0.0,motif->DimensionX(),motif->DimensionY()));
  }
  else
  {
    for ( Int_t i = 0; i <= AliMpConstants::ManuNofChannels(); ++i ) 
    {
      AliMpConnection* connection = motifType->FindConnectionByGassiNum(i);
      
      if ( connection ) 
      {
        Int_t ix = connection->GetLocalIx();
        Int_t iy = connection->GetLocalIy();
        
        Double_t x,y,dx,dy;
        
        motif->GetPadDimensionsByIndices(ix,iy,dx,dy);
        motif->PadPositionLocal(ix,iy,x,y);
        
        AliMUONPolygon* pol = new AliMUONPolygon(x,y,dx,dy);
        polygons.Add(pol);
      }
    }
  }
  
  AliMUONContourMaker maker;
  
  contour = maker.CreateContour(polygons);
  
  if (!contour || !contour->IsValid() ) 
  {
    AliError(Form("Failed to properly create contour %s contour = %p",mpName.Data(),contour));
    if ( contour ) 
    {
      AliError(Form("nofVertices=%d area.isvalid=%d",contour->NumberOfVertices(),contour->Area().IsValid()));
      StdoutToAliError(contour->Area().Print(););
    }
    delete contour;
    return 0x0;
  }
  
  {
    AliCodeTimerAuto("localmanucontour.add");
    fLocalManuContours.Add(new TObjString(mpName),contour);
  }
  
  return static_cast<AliMUONContour*>(contour->Clone());
}

//_____________________________________________________________________________
TObjArray* 
AliMUONManuContourMaker::GenerateManuContours(Bool_t stopAtError)
{
  /// Generate the contours for all the manus, taking into account the given transformation
  /// (to go from local to global). That transformation need not be the real one (i.e.
  /// it can be an "exploded" one to ease visualization).
  
  AliCodeTimerAuto("");
  
  TObjArray* manuContours = new TObjArray;
  
  manuContours->SetOwner(kTRUE);
  
  AliMpManuIterator it;
  Int_t detElemId, manuId;
  Int_t nmanus(0);
  Int_t nok(0);

  while ( it.Next(detElemId,manuId) ) 
  {
    ++nmanus;
    AliMUONContour* contour = CreateManuContour(detElemId,manuId);
    if (contour)
    {
      manuContours->Add(contour);
    }
    else
    {
      if ( stopAtError )
      {
        break;
      }
    }
    ++nok;
  }
  
  AliInfo(Form("%d manus. %d contours successfully created",nmanus,nok));
  
  return manuContours;
}

//_____________________________________________________________________________
TString
AliMUONManuContourMaker::NameIt(const AliMpMotifPosition& motifPosition) const
{
  /// Get the name of an AliMpMotifPosition
  
  AliMpVMotif* motif = motifPosition.GetMotif();
  TString name(Form("%s",motif->GetID().Data()));
  
  for ( Int_t i = 0; i < motif->GetNofPadDimensions(); ++i )
  {
    name += Form("/%7.3f-%7.3f:",motif->GetPadDimensionX(i),motif->GetPadDimensionY(i));
  }
  
  return name;
}

//_____________________________________________________________________________
TString
AliMUONManuContourMaker::ManuPathName(Int_t detElemId, Int_t manuId, Bool_t withCathodeName)
{
  /// Get the name of a manu
  
  AliMp::PlaneType planeType;
  if ( manuId & AliMpConstants::ManuMask(AliMp::kNonBendingPlane) )
  {
    planeType = AliMp::kNonBendingPlane;
  }
  else
  {
    planeType = AliMp::kBendingPlane;
  }
  AliMp::CathodType cathodeType = AliMpDEManager::GetCathod(detElemId,planeType);
  
  Int_t chamberId = AliMpDEManager::GetChamberId(detElemId);
  Int_t stationId = chamberId/2;
  
  Int_t busPatchId = AliMpDDLStore::Instance()->GetBusPatchId(detElemId, manuId);
  
  AliMpUID id(cathodeType,stationId,chamberId,detElemId,busPatchId,manuId);
  
  if ( withCathodeName ) return id.PathName();
  
  TString name(id.PathName());
  
  name.ReplaceAll("Cathode0/","");
  name.ReplaceAll("Cathode1/","");
  
  return name;
}




