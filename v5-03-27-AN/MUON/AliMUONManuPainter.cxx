
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

#include "AliMUONManuPainter.h"

#include "AliMpDCSNamer.h"
#include "AliLog.h"
#include "AliMUONContour.h"
#include "AliMUONManuContourMaker.h"
#include "AliMUONManuPadPainter.h"
#include "AliMUONPainterHelper.h"
#include "AliMUONVCalibParam.h"
#include "AliMUONVDigit.h"
#include "AliMUONVTrackerData.h"
#include "AliMpDEManager.h"
#include "AliMpManuUID.h"
#include "AliMpMotifPosition.h"
#include "AliMpMotifType.h"
#include "AliMpSlat.h"
#include "AliMpStationType.h"
#include "AliMpVMotif.h"
#include "AliMpVPadIterator.h"
#include <TArrayI.h>
#include <float.h>

///\class AliMUONManuPainter
///
/// Implementation of AliMUONVPainter for one manu (not the pads, only the manu
/// itself).
///
///\author Laurent Aphecetche, Subatech

///\cond CLASSIMP
ClassImp(AliMUONManuPainter)
///\endcond

//_____________________________________________________________________________
AliMUONManuPainter::AliMUONManuPainter(TRootIOCtor* ioCtor)
: AliMUONVPainter(ioCtor),
fDetElemId(-1),
fManuId(-1)
{
  /// ctor
}

//_____________________________________________________________________________
AliMUONManuPainter::AliMUONManuPainter()
: AliMUONVPainter(),
fDetElemId(-1),
fManuId(-1)
{
  /// ctor
}

//_____________________________________________________________________________
AliMUONManuPainter::AliMUONManuPainter(const AliMUONAttPainter& att,
                                       Int_t detElemId,
                                       Int_t manuId)
: AliMUONVPainter("MANU"),
  fDetElemId(detElemId),
  fManuId(manuId)
{
    /// ctor
    
  SetAttributes(att);
    
  AliMUONPainterHelper* h = AliMUONPainterHelper::Instance();
    
  SetID(detElemId,manuId);
  SetName(h->ManuName(manuId));
  SetPathName(h->ManuPathName(detElemId,manuId));
  
  
  AliMp::StationType stationType = AliMpDEManager::GetStationType(detElemId);
  
  if ( stationType == AliMp::kStationTrigger ) 
  {
    AliError("Hu ho. Not supposed to be used for trigger !");
    Invalidate();
    return;    
  }
    
  TString name = AliMUONManuContourMaker::ManuPathName(detElemId, manuId);

  AliMUONContour* contour = h->GetContour(name.Data());
  
  if (!contour)
  {
    AliError(Form("Could not get manuId %04d from DE %04d (name=%s)",manuId,detElemId,name.Data()));
  }
  
  SetContour(contour);
  
  Add(new AliMUONManuPadPainter(*this,fDetElemId,fManuId));
}

//_____________________________________________________________________________
AliMUONManuPainter::AliMUONManuPainter(const AliMUONManuPainter& rhs)
: AliMUONVPainter(rhs), fDetElemId(-1), fManuId(-1)
{
  /// copy ctor
  rhs.Copy(*this);
}

//_____________________________________________________________________________
AliMUONManuPainter&
AliMUONManuPainter::operator=(const AliMUONManuPainter& rhs)
{
  /// assignment operator
  if ( this != &rhs ) 
  {
    rhs.Copy(*this);
  }
  return *this;
}

//_____________________________________________________________________________
AliMUONManuPainter::~AliMUONManuPainter()
{
  /// dtor
}

//_____________________________________________________________________________
void 
AliMUONManuPainter::ComputeDataRange(const AliMUONVTrackerData& data, Int_t dataIndex, 
                                     Double_t& dataMin, Double_t& dataMax) const
{
  /// Compute data range spanned by this manu
  dataMin = dataMax = data.Manu(fDetElemId, fManuId, dataIndex);
}


//_____________________________________________________________________________
void
AliMUONManuPainter::Copy(TObject& object) const
{
  /// copyy this to object
  AliMUONVPainter::Copy((AliMUONVPainter&)(object));
  ((AliMUONManuPainter&)(object)).fDetElemId = fDetElemId;
  ((AliMUONManuPainter&)(object)).fManuId = fManuId;
}

//_____________________________________________________________________________
TString
AliMUONManuPainter::Describe(const AliMUONVTrackerData& data, Int_t dataIndex,
                             Double_t, Double_t)
{
  /// Describe data at this manu
  
  if (!data.HasManu(fDetElemId,fManuId)) return "";
  
  Double_t value = data.Manu(fDetElemId,fManuId,dataIndex);
  
  TString rv = AliMUONPainterHelper::Instance()->FormatValue(data.DimensionName(dataIndex).Data(),value);
  
  if ( TString(data.GetName()).Contains("HV") )
  {
    rv += "\n";
    
    AliMpDCSNamer hvNamer("TRACKER");
    
    if ( AliMpDEManager::GetStationType(fDetElemId) == AliMp::kStation12 )
    {
      Int_t sector = hvNamer.ManuId2Sector(fDetElemId,fManuId);

      rv += hvNamer.DCSChannelName(fDetElemId,sector);
    }
    else
    {
      rv += hvNamer.DCSChannelName(fDetElemId);
    }
  }
  
  return rv;
}

//_____________________________________________________________________________
void
AliMUONManuPainter::FillManuList(TObjArray& manuList) const
{
  /// Append our manu to the list
  manuList.Add(new AliMpManuUID(fDetElemId,fManuId));
}

//_____________________________________________________________________________
Bool_t
AliMUONManuPainter::IsIncluded() const
{
  /// whether this manu is included in the readout or not
  return ( InteractiveReadOutConfig()->Manu(fDetElemId,fManuId) > 0 );
}

//_____________________________________________________________________________
void
AliMUONManuPainter::PaintArea(const AliMUONVTrackerData& data, Int_t dataIndex,
                              Double_t min, Double_t max)
{
  /// Paint area of this manu according to data
  
  if (!data.HasManu(fDetElemId,fManuId)) return;

  Double_t value = data.Manu(fDetElemId,fManuId,dataIndex);
  
  if ( value >= AliMUONVCalibParam::InvalidFloatValue() ) return;
  
  Int_t color = AliMUONPainterHelper::Instance()->ColorFromValue(value,min,max);
  
  PaintArea(color);
}

//_____________________________________________________________________________
AliMUONAttPainter
AliMUONManuPainter::Validate(const AliMUONAttPainter& attributes) const
{
  /// Normalize the attributes
  
  /// check that cathode and plane are up-to-date, and that they are legal
  
  AliMUONAttPainter norm(attributes);

  norm.SetValid(kFALSE);
  
  if ( norm.IsCathodeDefined() ) 
  {
    if ( norm.IsCathode0() != Attributes().IsCathode0() ) return norm;
  }
  
  if ( norm.IsPlaneDefined() ) 
  {
    if ( norm.IsBendingPlane() != Attributes().IsBendingPlane() ) return norm;
  }
  
  norm.SetValid(kTRUE);
  
  if ( norm.IsCathodeDefined() && !norm.IsPlaneDefined() ) 
  {
    // derive plane from cathode
    
    AliMp::CathodType cathode = ( norm.IsCathode0() ? AliMp::kCath0 : AliMp::kCath1 ) ;

    AliMp::PlaneType planeType = AliMpDEManager::GetPlaneType(fDetElemId,cathode);
    
    Bool_t bending = ( planeType == AliMp::kBendingPlane );
    
    norm.SetPlane(bending,!bending);    
  }
  else if ( norm.IsPlaneDefined() && !norm.IsCathodeDefined() )
  {
    // derive cathode from plane

    AliMp::PlaneType planeType = ( norm.IsBendingPlane() ? AliMp::kBendingPlane : AliMp::kNonBendingPlane );

    AliMp::CathodType cathode = AliMpDEManager::GetCathod(fDetElemId,planeType);
          
    Bool_t cath0 = ( cathode == AliMp::kCath0 );
    
    norm.SetCathode(cath0,!cath0);    
  }
  else if ( norm.IsPlaneDefined() && norm.IsCathodeDefined() ) 
  {
    // check that cathode and plane matches
    
    AliMp::CathodType cathode = ( norm.IsCathode0() ? AliMp::kCath0 : AliMp::kCath1 ) ;

    AliMp::PlaneType planeType = AliMpDEManager::GetPlaneType(fDetElemId,cathode);
    
    Bool_t bending = ( planeType == AliMp::kBendingPlane );
    
    if ( bending != norm.IsBendingPlane() ) 
    {
      norm.SetValid(kFALSE);
    }
  }
  
  return norm;
}


