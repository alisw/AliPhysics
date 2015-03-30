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

#include "AliMUONBusPatchPainter.h"

#include "AliMUONManuPainter.h"
#include "AliMUONContour.h"
#include "AliMUONPainterHelper.h"
#include "AliMUONVCalibParam.h"
#include "AliMUONVTrackerData.h"
#include "AliMpBusPatch.h"
#include "AliMpConstants.h"
#include "AliMpDDLStore.h"
#include "AliMpDEManager.h"
#include "AliMpPlaneType.h"
#include "AliLog.h"
#include <TObjArray.h>
#include <TString.h>
#include <float.h>

/// \class AliMUONBusPatchPainter
///
/// Painter for one bus patch. Actually possibly for only part of one
/// buspatch (the part that is on the plane/cathode requested when
/// creating the painter)
///
/// \author Laurent Aphecetche, Subatech

///\cond CLASSIMP
ClassImp(AliMUONBusPatchPainter)
///\endcond

//_____________________________________________________________________________
AliMUONBusPatchPainter::AliMUONBusPatchPainter()
: AliMUONVPainter(),
fBusPatchId(-1)
{
  /// default ctor
}

//_____________________________________________________________________________
AliMUONBusPatchPainter::AliMUONBusPatchPainter(TRootIOCtor* ioCtor)
: AliMUONVPainter(ioCtor),
fBusPatchId(-1)
{
  /// default streaming ctor
}

//_____________________________________________________________________________
AliMUONBusPatchPainter::AliMUONBusPatchPainter(const AliMUONAttPainter& att, 
                                               Int_t busPatchId)
: AliMUONVPainter("BUSPATCH"),
fBusPatchId(busPatchId)
{
  /// normal ctor
  /// WARNING : the construction of this object can fail.
  /// You MUST check the IsValid() method afterwards (real world would
  /// be to use exception, but well, whether we should use exceptions
  /// in aliroot is still unclear to me.
  
  SetAttributes(Validate(att));
  
  AliMp::PlaneType planeType = ( Attributes().IsBendingPlane() ? AliMp::kBendingPlane : AliMp::kNonBendingPlane );

  Int_t detElemId = AliMpDDLStore::Instance()->GetDEfromBus(busPatchId);

  AliMUONPainterHelper* h = AliMUONPainterHelper::Instance();
  
  SetID(busPatchId,-1);
  SetName(h->BusPatchName(busPatchId));
  SetPathName(h->BusPatchPathName(busPatchId));
  
  AliMpBusPatch* busPatch = AliMpDDLStore::Instance()->GetBusPatch(fBusPatchId);
  
  Int_t mask = AliMpConstants::ManuMask(AliMp::kNonBendingPlane);
  
  AliMUONContour* bpContour = h->GetContour(ContourName());
  
  Double_t xmin(FLT_MAX), ymin(FLT_MAX), xmax(-FLT_MAX), ymax(-FLT_MAX);
  
  TObjArray contours;
  
  Int_t nmanus(0);
  
  for ( Int_t i = 0; i < busPatch->GetNofManus(); ++i ) 
  {
    Int_t manuId = busPatch->GetManuId(i);
    
    Bool_t correctPlane(kTRUE);
    
    if ( planeType == AliMp::kNonBendingPlane ) 
    {
      if ( ( manuId & mask ) == 0 ) correctPlane = kFALSE;
    }
    else
    {
      if ( ( manuId & mask ) == mask ) correctPlane = kFALSE;
    }

    if (!correctPlane) continue;
    
    ++nmanus;
    
    AliMUONVPainter* painter = new AliMUONManuPainter(Attributes(),
                                                      busPatch->GetDEId(),
                                                      manuId);
    
    Add(painter);
    
    const AliMpArea& area = painter->Area();
    
    xmin = TMath::Min(xmin,area.LeftBorder());
    ymin = TMath::Min(ymin,area.DownBorder());
    xmax = TMath::Max(xmax,area.RightBorder());
    ymax = TMath::Max(ymax,area.UpBorder());
    
    if (!bpContour)
    {
      contours.Add(painter->Contour());
    }
  }
  
  if ( !nmanus )
  {
    Invalidate();
    return;
  }
    
  if (!bpContour)
  {
    bpContour = h->MergeContours(contours,ContourName());
    if (!bpContour)
    {
      AliError("Could not merge those contours");
      StdoutToAliError(contours.Print(););
    }
  }
  
  SetContour(bpContour);
}

//_____________________________________________________________________________
void 
AliMUONBusPatchPainter::ComputeDataRange(const AliMUONVTrackerData& data, Int_t dataIndex, 
                                         Double_t& dataMin, Double_t& dataMax) const
{
  /// Compute the data range spanned by this bus patch (on this cathode or plane)
  dataMin = dataMax = data.BusPatch(fBusPatchId, dataIndex);
}

//_____________________________________________________________________________
AliMUONBusPatchPainter::AliMUONBusPatchPainter(const AliMUONBusPatchPainter& rhs)
: AliMUONVPainter(rhs), fBusPatchId(-1)
{
  /// Copy ctor
  rhs.Copy(*this);
}

//_____________________________________________________________________________
AliMUONBusPatchPainter&
AliMUONBusPatchPainter::operator=(const AliMUONBusPatchPainter& rhs)
{
  /// Assignment operator
  if ( this != &rhs ) 
  {
    rhs.Copy(*this);
  }
  return *this;
}

//_____________________________________________________________________________
AliMUONBusPatchPainter::~AliMUONBusPatchPainter()
{
  /// dtor
}

//_____________________________________________________________________________
void
AliMUONBusPatchPainter::Copy(TObject& object) const
{
  /// Copy this to object
  AliMUONVPainter::Copy((AliMUONVPainter&)(object));
  ((AliMUONBusPatchPainter&)(object)).fBusPatchId = fBusPatchId;
}

//_____________________________________________________________________________
Bool_t
AliMUONBusPatchPainter::IsIncluded() const
{
  /// whether this bus patch is included in the readout or not
  return ( InteractiveReadOutConfig()->BusPatch(fBusPatchId) > 0 );
}

//_____________________________________________________________________________
TString
AliMUONBusPatchPainter::Describe(const AliMUONVTrackerData& data, Int_t dataIndex, 
                                 Double_t, Double_t)
{
  /// Text about data
  
  if (!data.HasBusPatch(fBusPatchId)) return "";
  
  Double_t value = data.BusPatch(fBusPatchId,dataIndex);
  
  return AliMUONPainterHelper::Instance()->FormatValue(data.DimensionName(dataIndex).Data(),value);
}

//_____________________________________________________________________________
void
AliMUONBusPatchPainter::PaintArea(const AliMUONVTrackerData& data, Int_t dataIndex,
                                  Double_t min, Double_t max)
{
  /// Paint area of this buspatch according to the data
  
  if (!data.HasBusPatch(fBusPatchId)) return;
  
  Double_t value = data.BusPatch(fBusPatchId,dataIndex);
  
  if ( value >= AliMUONVCalibParam::InvalidFloatValue() ) return;
  
  Int_t color = AliMUONPainterHelper::Instance()->ColorFromValue(value,min,max);
  
  PaintArea(color);
}

//_____________________________________________________________________________
AliMUONAttPainter 
AliMUONBusPatchPainter::Validate(const AliMUONAttPainter& attributes) const
{
  /// Normalize attributes
  
  // we invalidate the attributes, if we have no manu in the requested plane
  // and we cross-check that both cathode and plane are up-to-date
  
  AliMUONAttPainter norm(attributes);
  
  Int_t detElemId = AliMpDDLStore::Instance()->GetDEfromBus(fBusPatchId);

  if (!norm.IsValid()) return norm;
  
  if ( !norm.IsCathodeDefined() )
  {
    AliMp::PlaneType planeType = ( norm.IsBendingPlane() ? AliMp::kBendingPlane : AliMp::kNonBendingPlane );
  
    AliMp::CathodType cathodeType = AliMpDEManager::GetCathod(detElemId,planeType);
    
    Bool_t cath0 = ( cathodeType == AliMp::kCath0 );
    
    norm.SetCathode(cath0,!cath0);
  }
  else if ( !norm.IsPlaneDefined() )
  {
    AliMp::CathodType cathodeType = ( norm.IsCathode0() ? AliMp::kCath0 : AliMp::kCath1 );
    
    AliMp::PlaneType planeType = AliMpDEManager::GetPlaneType(detElemId,cathodeType);
    
    Bool_t bending = ( planeType == AliMp::kBendingPlane );

    norm.SetPlane(bending,!bending);    
  }
  
  AliMpBusPatch* busPatch = AliMpDDLStore::Instance()->GetBusPatch(fBusPatchId);
  
  Int_t mask = AliMpConstants::ManuMask(AliMp::kNonBendingPlane);
  
  Int_t nb(0);
  Int_t b(0);
  
  for ( Int_t i = 0; i < busPatch->GetNofManus(); ++i ) 
  {
    Int_t manuId = busPatch->GetManuId(i);
    
    if ( manuId & mask  ) ++nb;
    else ++b;
  }
  
  if ( norm.IsBendingPlane() && !b ) norm.SetValid(kFALSE);
  if ( norm.IsNonBendingPlane() && !nb ) norm.SetValid(kFALSE);
  
  return norm;
}

