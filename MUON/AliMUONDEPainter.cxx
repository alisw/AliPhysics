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

#include "AliMUONDEPainter.h"

#include "AliMUONBusPatchPainter.h"
#include "AliMUONGeometryTransformer.h"
#include "AliMUONPCBPainter.h"
#include "AliMUONContour.h"
#include "AliMUONPainterHelper.h"
#include "AliMUONVCalibParam.h"
#include "AliMUONVTrackerData.h"
#include "AliMUONObjectPair.h"
#include "AliMpDDLStore.h"
#include "AliMpDEManager.h"
#include "AliMpDetElement.h"
#include "AliMpPCB.h"
#include "AliMpSector.h"
#include "AliMpSlat.h"
#include "AliLog.h"
#include <TObjString.h>

/// \class AliMUONDEPainter
///
/// Painter for one detection element
///
/// It draws a given plane (bending or non bending) of a given detection element
///
/// \author Laurent Aphecetche, Subatech

///\cond CLASSIMP
ClassImp(AliMUONDEPainter)
///\endcond

//_____________________________________________________________________________
AliMUONDEPainter::AliMUONDEPainter()
: AliMUONVPainter(),
fDetElemId(-1)
{
  /// default ctor
}

//_____________________________________________________________________________
AliMUONDEPainter::AliMUONDEPainter(TRootIOCtor* ioCtor)
: AliMUONVPainter(ioCtor),
fDetElemId(-1)
{
  /// default streaming ctor
}

//_____________________________________________________________________________
AliMUONDEPainter::AliMUONDEPainter(const AliMUONAttPainter& att, Int_t detElemId)
: AliMUONVPainter("DE"),
fDetElemId(detElemId)
{
  /// normal ctor

  AliMUONAttPainter deAtt(att);
  
  if ( att.IsCathodeDefined() )
  {
    AliMp::CathodType cathodType = ( att.IsCathode0() ? AliMp::kCath0 : AliMp::kCath1 ) ;
    
    Bool_t cath0 = ( cathodType == AliMp::kCath0 ) ;
    
    AliMp::PlaneType planeType = AliMpDEManager::GetPlaneType(detElemId,cathodType);
    
    Bool_t bending = (  planeType == AliMp::kBendingPlane );
    
    deAtt.SetCathode(cath0,!cath0);
    deAtt.SetPlane(bending,!bending);

  }
  
  if ( att.IsPlaneDefined() ) 
  {  
    AliMp::PlaneType planeType = ( att.IsBendingPlane() ? AliMp::kBendingPlane : AliMp::kNonBendingPlane );
  
    Bool_t bending = ( planeType == AliMp::kBendingPlane );
  
    Bool_t cath0 = ( AliMpDEManager::GetCathod(detElemId,planeType) == AliMp::kCath0 );
    
    deAtt.SetCathode(cath0,!cath0);
    deAtt.SetPlane(bending,!bending);

  }
  
  deAtt.SetCathodeAndPlaneMutuallyExclusive(kFALSE);
                                            
  SetAttributes(deAtt);
  
  AliMUONPainterHelper* h = AliMUONPainterHelper::Instance();
  
  SetID(detElemId,-1);
  SetName(h->DEName(fDetElemId).Data());
  SetPathName(h->DEPathName(fDetElemId).Data());
              
  AliMp::PlaneType planeType = ( Attributes().IsBendingPlane() ? AliMp::kBendingPlane : AliMp::kNonBendingPlane );
  
  Double_t x,y,z;
  
  if ( AliMpDEManager::GetStationType(DetElemId()) == AliMp::kStation345 ) 
  {
    const AliMpSlat* slat = h->GetSlat(DetElemId(),planeType);
  
    for ( Int_t i = 0; i < slat->GetSize(); ++i ) 
    {
      Add(new AliMUONPCBPainter(Attributes(),DetElemId(),i));
    }
    
    AliMUONPainterHelper::Instance()->Local2Global(fDetElemId,0.0,0.0,0.0,x,y,z);    
  }
  else if ( AliMpDEManager::GetStationType(DetElemId()) != AliMp::kStationTrigger )
  {
    const AliMpSector* sector = h->GetSector(DetElemId(),planeType);

    Double_t xl(sector->GetDimensionX());
    Double_t yl(sector->GetDimensionY());
    
    h->Local2Global(fDetElemId,xl,yl,0.0,x,y,z);
  }
  else
  {
    AliFatal("Not implemented for trigger !!!");
  }
  
  AliMUONContour* contour = h->GetContour(ContourName());
  
  TObjArray contourArray;
    
  AliMpDetElement* de = AliMpDDLStore::Instance()->GetDetElement(fDetElemId);
  
  AliDebug(1,Form("de %p n %d",de,de->GetNofBusPatches()));
  
  for ( Int_t i = 0; i < de->GetNofBusPatches(); ++i ) 
  {
    AliMUONBusPatchPainter* painter = new AliMUONBusPatchPainter(Attributes(),de->GetBusPatchId(i));
                                                                 
    if ( !painter->IsValid() ) 
    {
      AliDebug(1,Form("Skipping BP %d which seem to have no manu in plane %s",
                   de->GetBusPatchId(i),
                   Attributes().IsBendingPlane() ? "bending" : "non bending"));
      delete painter;
      continue;
    }
    else
    {
      Add(painter);
    }
    
    if ( !contour ) 
    {
      contourArray.Add(painter->Contour());
    }
  }
  
  if (!contour)
  {
    contour = h->MergeContours(contourArray,ContourName());
    if (!contour)
    {
      AliError(Form("%s : could not merge those contours",Name().Data()));
      StdoutToAliError(contourArray.Print(););
    }
  }
  
  SetContour(contour);  
}

//_____________________________________________________________________________
AliMUONDEPainter::AliMUONDEPainter(const AliMUONDEPainter& rhs):
AliMUONVPainter(rhs), fDetElemId(-1)
{
  /// copy ctor
  rhs.Copy(*this);
}

//_____________________________________________________________________________
AliMUONDEPainter& 
AliMUONDEPainter::operator=(const AliMUONDEPainter& rhs)
{
  /// assignment operator
  if ( this != &rhs ) 
  {
    rhs.Copy(*this);
  }
  return *this;
}

//_____________________________________________________________________________
AliMUONDEPainter::~AliMUONDEPainter()
{
  /// dtor = nop
}

//_____________________________________________________________________________
void 
AliMUONDEPainter::ComputeDataRange(const AliMUONVTrackerData& data, Int_t dataIndex, 
                                         Double_t& dataMin, Double_t& dataMax) const
{
  /// Compute the data range spanned by this detection element
  dataMin = dataMax = data.DetectionElement(fDetElemId, dataIndex);
}

//_____________________________________________________________________________
void
AliMUONDEPainter::Copy(TObject& object) const
{
  /// Copy this to object
  AliMUONVPainter::Copy((AliMUONVPainter&)(object));
  ((AliMUONDEPainter&)(object)).fDetElemId = fDetElemId;
}

//_____________________________________________________________________________
TString
AliMUONDEPainter::Describe(const AliMUONVTrackerData& data, Int_t dataIndex,
                             Double_t, Double_t)
{
  /// Describe data at this detection element
  
  if (!data.HasDetectionElement(fDetElemId)) return "";
  
  Double_t value = data.DetectionElement(fDetElemId,dataIndex);
  
  return AliMUONPainterHelper::Instance()->FormatValue(data.DimensionName(dataIndex).Data(),value);
}

//_____________________________________________________________________________
void
AliMUONDEPainter::FillManuList(TObjArray& manuList) const
{
  /// Fill (append to) manu list
  TIter next(Children());
  AliMUONVPainter* p;
  
  while ( ( p = static_cast<AliMUONVPainter*>(next()) ) )
  {
    if ( p->IsA() == AliMUONBusPatchPainter::Class() )
    {
      // Only consider bus patch painters (and not PCB ones),
      // in order not to double count some manus
      p->FillManuList(manuList);
    }
  }
}

//_____________________________________________________________________________
Bool_t
AliMUONDEPainter::IsIncluded() const
{
  /// whether this detection element is included in the readout or not
  return ( InteractiveReadOutConfig()->DetectionElement(fDetElemId) > 0 );
}

//_____________________________________________________________________________
void
AliMUONDEPainter::PaintArea(const AliMUONVTrackerData& data, Int_t dataIndex,
                                  Double_t min, Double_t max)
{
  /// Paint the area of this detection element
  
  if (!data.HasDetectionElement(fDetElemId)) return;
  
  Double_t value = data.DetectionElement(fDetElemId,dataIndex);
  
  if ( value >= AliMUONVCalibParam::InvalidFloatValue() ) return;
  
  Int_t color = AliMUONPainterHelper::Instance()->ColorFromValue(value,min,max);
  
  PaintArea(color);
}

//_____________________________________________________________________________
AliMUONAttPainter 
AliMUONDEPainter::Validate(const AliMUONAttPainter& attributes) const
{
  /// Normalize attributes
  
  AliMUONAttPainter norm(attributes);
  
  norm.SetCathodeAndPlaneMutuallyExclusive(kFALSE);
  
  if ( norm.IsCathodeDefined() && !norm.IsPlaneDefined() )
  {
    // only cathode known : derive the plane
    
    AliMp::CathodType cathodType = ( norm.IsCathode0() ? AliMp::kCath0 : AliMp::kCath1 );
    
    AliMp::PlaneType planeType = AliMpDEManager::GetPlaneType(fDetElemId,cathodType);

    Bool_t bending = ( planeType == AliMp::kBendingPlane ) ;
    
    norm.SetPlane(bending,!bending);
  }
  
  else if ( !norm.IsCathodeDefined() && norm.IsPlaneDefined() )
  {
    // only plane is known : derive the cathode
    
    AliMp::PlaneType planeType = ( norm.IsBendingPlane() ? AliMp::kBendingPlane : AliMp::kNonBendingPlane );
        
    Bool_t cath0 = ( AliMpDEManager::GetCathod(fDetElemId,planeType) == AliMp::kCath0 );
    
    norm.SetCathode(cath0,!cath0);
  }  
  else    
  {
    // check that both information are compatible
    
    AliMp::PlaneType planeType = ( norm.IsBendingPlane() ? AliMp::kBendingPlane : AliMp::kNonBendingPlane );

    AliMp::CathodType cathode = AliMpDEManager::GetCathod(fDetElemId,planeType);
    
    if ( (cathode == AliMp::kCath0 && norm.IsCathode1()) ||
         (cathode == AliMp::kCath1 && norm.IsCathode0()) ) 
    {
      norm.SetValid(kFALSE);
    }
  }
  
  return norm;
}

