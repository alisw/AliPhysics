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

/// \class AliMUONPCBPainter
///
/// Implementation of AliMUONVPainter for slat's PCBs
/// 
/// \author Laurent Aphecetche, Subatech

#include "AliMUONPCBPainter.h"

#include "AliMUONManuPainter.h"
#include "AliMUONContour.h"
#include "AliMUONPainterHelper.h"
#include "AliMUONVCalibParam.h"
#include "AliMUONVTrackerData.h"
#include "AliMpDEManager.h"
#include "AliMpMotifPosition.h"
#include "AliMpPCB.h"
#include "AliMpPlaneType.h"
#include "AliMpSlat.h"
#include "AliLog.h"

/// \cond CLASSIMP
ClassImp(AliMUONPCBPainter)
/// \endcond

//_____________________________________________________________________________
AliMUONPCBPainter::AliMUONPCBPainter(TRootIOCtor* ioCtor)
: AliMUONVPainter(ioCtor),
fDetElemId(-1),
fPCBIndex(-1)
{
  /// root streaming ctor
}

//_____________________________________________________________________________
AliMUONPCBPainter::AliMUONPCBPainter()
: AliMUONVPainter(),
fDetElemId(-1),
fPCBIndex(-1)
{
  /// empty ctor
}

//_____________________________________________________________________________
AliMUONPCBPainter::AliMUONPCBPainter(const AliMUONAttPainter& att,
                                     Int_t detElemId,
                                     Int_t pcbNumber)
: AliMUONVPainter("PCB"),
  fDetElemId(detElemId),
  fPCBIndex(pcbNumber)
{
  /// Ctor
  
  SetAttributes(att);
  
  AliMUONPainterHelper* h = AliMUONPainterHelper::Instance();
  
  AliMp::PlaneType planeType = ( att.IsBendingPlane() ? AliMp::kBendingPlane : AliMp::kNonBendingPlane );
  
  const AliMpSlat* slat = AliMUONPainterHelper::Instance()->GetSlat(fDetElemId,planeType);
  
  SetID(detElemId,pcbNumber);
  SetName(h->PCBName(pcbNumber));
  SetPathName(h->PCBPathName(detElemId,pcbNumber));
  
  AliMpPCB* pcb = slat->GetPCB(fPCBIndex);
  
  AliMUONContour* contour = h->GetContour(ContourName());
  TObjArray contourArray;
  
  for ( Int_t imp = 0 ; imp < pcb->GetSize(); ++imp ) 
  {
    AliMpMotifPosition* mp = pcb->GetMotifPosition(imp);
    AliMUONVPainter* painter = new AliMUONManuPainter(Attributes(),fDetElemId,mp->GetID());
    Add(painter);
    if (!contour)
    {
      contourArray.Add(painter->Contour());
    }
  }
  
  Double_t x,y,z;
  
  h->Local2Global(fDetElemId,
                  pcb->X()-slat->GetPositionX(),
                  pcb->Y()-slat->GetPositionY(),
                  0.0,
                  x,y,z);
  
  if (!contour)
  {
    contour = h->MergeContours(contourArray,ContourName());
  }
  
  SetContour(contour);
}

//_____________________________________________________________________________
AliMUONPCBPainter::~AliMUONPCBPainter()
{
  /// dtor
}

//_____________________________________________________________________________
AliMUONPCBPainter::AliMUONPCBPainter(const AliMUONPCBPainter& rhs)
: AliMUONVPainter(rhs),
  fDetElemId(-1),
  fPCBIndex(-1)
{
    /// copy ctor
    rhs.Copy(*this);
}

//_____________________________________________________________________________
AliMUONPCBPainter&
AliMUONPCBPainter::operator=(const AliMUONPCBPainter& rhs)
{
  /// assignment operator
  if ( this != &rhs )
  {
    rhs.Copy(*this);
  }
  return *this;
}

//_____________________________________________________________________________
void 
AliMUONPCBPainter::ComputeDataRange(const AliMUONVTrackerData& data, Int_t dataIndex, 
                                   Double_t& dataMin, Double_t& dataMax) const
{
  /// Compute the min and max of this PCB data 
  dataMin = dataMax = data.PCB(fDetElemId, fPCBIndex,dataIndex);
}

//_____________________________________________________________________________
void
AliMUONPCBPainter::Copy(TObject& object) const
{
  /// Copy this to object
  AliMUONVPainter::Copy((AliMUONVPainter&)(object));
  ((AliMUONPCBPainter&)(object)).fDetElemId = fDetElemId;
  ((AliMUONPCBPainter&)(object)).fPCBIndex = fPCBIndex;
}

//_____________________________________________________________________________
TString
AliMUONPCBPainter::Describe(const AliMUONVTrackerData& data, Int_t dataIndex,
                                Double_t, Double_t)
{
  /// Describe data at this PCB
  
  if (!data.HasPCB(fDetElemId,fPCBIndex)) return "";
  
  Double_t value = data.PCB(fDetElemId,fPCBIndex,dataIndex);
  
  return AliMUONPainterHelper::Instance()->FormatValue(data.DimensionName(dataIndex).Data(),value);
}

//_____________________________________________________________________________
Bool_t
AliMUONPCBPainter::IsIncluded() const
{
  /// Whether this PCB is included in the read out or not
  return ( InteractiveReadOutConfig()->PCB(fDetElemId,fPCBIndex) > 0 );
}

//_____________________________________________________________________________
void
AliMUONPCBPainter::PaintArea(const AliMUONVTrackerData& data, Int_t dataIndex,
                            Double_t min, Double_t max)
{
  /// Fill the contour of this PCB with a color depending of the value of the data
  if (!data.HasPCB(fDetElemId,fPCBIndex)) return;
  
  Double_t value = data.PCB(fDetElemId,fPCBIndex,dataIndex);
  
  if ( value >= AliMUONVCalibParam::InvalidFloatValue() ) return;
  
  Int_t color = AliMUONPainterHelper::Instance()->ColorFromValue(value,min,max);
  
  PaintArea(color);
}

