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

#include "AliMUONChamberPainter.h"

#include "AliMUONDEPainter.h"
#include "AliMUONContour.h"
#include "AliMUONPainterHelper.h"
#include "AliMUONVCalibParam.h"
#include "AliMUONVTrackerData.h"
#include "AliMpConstants.h"
#include "AliMpDEIterator.h"
#include "AliMpDEManager.h"
#include "AliMpPlaneType.h"
#include "AliMpSegmentation.h"
#include "AliMpStationType.h"
#include "AliMpVSegmentation.h"
#include "AliMUONObjectPair.h"
#include "AliLog.h"
#include <Riostream.h>
#include <TObjString.h>
#include <TArrayI.h>
#include <cassert>
#include <float.h>

/// \class AliMUONChamberPainter
///
/// Painter for one plane/cathode of one chamber
///
/// \author Laurent Aphecetche, Subatech

///\cond CLASSIMP
ClassImp(AliMUONChamberPainter)
///\endcond

//_____________________________________________________________________________
AliMUONChamberPainter::AliMUONChamberPainter()
: AliMUONVPainter(),
fChamberId(-1)
{
  /// ctor
}

//_____________________________________________________________________________
AliMUONChamberPainter::AliMUONChamberPainter(TRootIOCtor* ioCtor)
: AliMUONVPainter(ioCtor),
fChamberId(-1)
{
  /// ctor
}

//_____________________________________________________________________________
AliMUONChamberPainter::AliMUONChamberPainter(const AliMUONAttPainter& att, 
                                             Int_t chamberId)
: AliMUONVPainter("Chamber"),
fChamberId(chamberId)
{
  /// ctor

  AliMUONPainterHelper* h = AliMUONPainterHelper::Instance(); // to be sure mapping is loaded...

  AliMUONAttPainter chAtt(att);
  
  chAtt.SetCathodeAndPlaneMutuallyExclusive(kTRUE);
  
  SetAttributes(chAtt);
  
  SetID(chamberId,-1);
  SetName(h->ChamberName(fChamberId).Data());
  SetPathName(h->ChamberPathName(fChamberId).Data());
              
  AliMpDEIterator deIt;
  
  deIt.First(fChamberId);
  
  AliMUONContour* contour = h->GetContour(ContourName());
  TObjArray contourArray;
  
  while (!deIt.IsDone())
  {
    Int_t detElemId = deIt.CurrentDEId();

    AliMUONAttPainter deAtt(att);
    
    if ( att.IsCathodeDefined() ) 
    {
      deAtt.SetCathode(kFALSE,kFALSE);
      AliMp::PlaneType planeType;
      
      if ( att.IsCathode0() ) planeType = AliMpDEManager::GetPlaneType(detElemId,AliMp::kCath0);
      else planeType = AliMpDEManager::GetPlaneType(detElemId,AliMp::kCath1);
      
      Bool_t bending = ( planeType == AliMp::kBendingPlane );
      
      deAtt.SetPlane(bending,!bending);
    }
    
    assert(deAtt.IsPlaneDefined());
    
    AliMUONVPainter* painter = new AliMUONDEPainter(deAtt,detElemId);

    Add(painter);
    
    if (!contour)
    {
      contourArray.Add(painter->Contour());
    }
    
    deIt.Next();
  }
  
  Double_t xmin(1E9), xmax(-1E9), ymin(1E9), ymax(-1E9);
  TIter next(Children());
  AliMUONVPainter* painter;
  
  while ( ( painter = static_cast<AliMUONVPainter*>(next()) ) )
  {
    const AliMpArea& area = painter->Area();
    xmin = TMath::Min(xmin,area.LeftBorder());
    xmax = TMath::Max(xmax,area.RightBorder());
    ymin = TMath::Min(ymin,area.DownBorder());
    ymax = TMath::Max(ymax,area.UpBorder());
  }
  
  if ( contourArray.GetLast() >= 0 ) 
  {
    contour = h->MergeContours(contourArray,ContourName());
  }
  
  SetContour(contour);    
}

//_____________________________________________________________________________
AliMUONChamberPainter::AliMUONChamberPainter(const AliMUONChamberPainter& rhs)
: AliMUONVPainter(rhs),
fChamberId(rhs.fChamberId)
{
  /// copy ctor
  rhs.Copy(*this);
}

//_____________________________________________________________________________
AliMUONChamberPainter&
AliMUONChamberPainter::operator=(const AliMUONChamberPainter& rhs)
{
  /// assignment operator
  if ( this != &rhs )
  {
    rhs.Copy(*this);
  }
  return *this;
}

//_____________________________________________________________________________
AliMUONChamberPainter::~AliMUONChamberPainter()
{
  /// dtor
}

//_____________________________________________________________________________
void 
AliMUONChamberPainter::ComputeDataRange(const AliMUONVTrackerData& data, Int_t dataIndex, 
                                        Double_t& dataMin, Double_t& dataMax) const
{
  /// Compute data range spanned by this (plane of that) chamber
  dataMin = dataMax = data.Chamber(fChamberId, dataIndex);
}


//_____________________________________________________________________________
void
AliMUONChamberPainter::Copy(TObject& object) const
{
  /// Copy this to object
  AliMUONVPainter::Copy((AliMUONVPainter&)(object));
  ((AliMUONChamberPainter&)(object)).fChamberId = fChamberId;
}

//_____________________________________________________________________________
TString
AliMUONChamberPainter::Describe(const AliMUONVTrackerData& data, Int_t dataIndex,
                           Double_t, Double_t)
{
  /// Describe data at this chamber
  
  if (!data.HasChamber(fChamberId)) return "";
  
  Double_t value = data.Chamber(fChamberId,dataIndex);
  
  return AliMUONPainterHelper::Instance()->FormatValue(data.DimensionName(dataIndex).Data(),value);
}

//_____________________________________________________________________________
Bool_t 
AliMUONChamberPainter::IsIncluded() const
{
  /// whether this chamber is included in the readout or not
  return ( InteractiveReadOutConfig()->Chamber(fChamberId) > 0 );
}

//_____________________________________________________________________________
void
AliMUONChamberPainter::PaintArea(const AliMUONVTrackerData& data, Int_t dataIndex,
                                 Double_t min, Double_t max)
{
  /// Paint area of this chamber according to data
  
  if (!data.HasChamber(fChamberId)) return;
  
  Double_t value = data.Chamber(fChamberId,dataIndex);
  
  if ( value >= AliMUONVCalibParam::InvalidFloatValue() ) return;
  
  Int_t color = AliMUONPainterHelper::Instance()->ColorFromValue(value,min,max);
  
  PaintArea(color);
}

//_____________________________________________________________________________
AliMUONAttPainter 
AliMUONChamberPainter::Validate(const AliMUONAttPainter& attributes) const
{
  /// Normalize attributes
  
  AliMUONAttPainter norm(attributes);
  
  // A chamber painter must be either cathode defined or plane defined
  
  if ( norm.IsCathodeDefined() && norm.IsPlaneDefined() ) 
  {
    norm.SetValid(kFALSE);
  }

  if ( !norm.IsCathodeDefined() && !norm.IsPlaneDefined() ) 
  {
    norm.SetValid(kFALSE);
  }
  
  return norm;
}
