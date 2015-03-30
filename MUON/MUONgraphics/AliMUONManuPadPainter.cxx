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

#include "AliMUONManuPadPainter.h"

#include "AliLog.h"
#include "AliMUONPainterGroup.h"
#include "AliMUONPainterHelper.h"
#include "AliMUONTrackerDataHistogrammer.h"
#include "AliMUONVCalibParam.h"
#include "AliMUONVDigit.h"
#include "AliMUONVTrackerData.h"
#include "AliMpConnection.h"
#include "AliMpConstants.h"
#include "AliMpDDLStore.h"
#include "AliMpPad.h"
#include "AliMpSegmentation.h"
#include "AliMpVSegmentation.h"
#include <TCanvas.h>
#include <TH1.h>
#include <TVirtualPad.h>
#include <TVirtualX.h>
#include <float.h>

///\class AliMUONManuPadPainter
///
/// Painter for the pads of one manu
///
///\author Laurent Aphecetche, Subatech

///\cond CLASSIMP
ClassImp(AliMUONManuPadPainter)
///\endcond

//_____________________________________________________________________________
AliMUONManuPadPainter::AliMUONManuPadPainter()
: AliMUONVPainter(), 
fDetElemId(-1), 
fManuId(-1),
fLineColorBck(-1),
fLineWidthBck(-1),
fFillColorBck(-1),
fFillStyleBck(-1)
{
  /// ctor
}

//_____________________________________________________________________________
AliMUONManuPadPainter::AliMUONManuPadPainter(TRootIOCtor* ioCtor)
: AliMUONVPainter(ioCtor), 
fDetElemId(-1), 
fManuId(-1),
fLineColorBck(-1),
fLineWidthBck(-1),
fFillColorBck(-1),
fFillStyleBck(-1)
{
  /// ctor
}

//_____________________________________________________________________________
AliMUONManuPadPainter::AliMUONManuPadPainter(const AliMUONVPainter& mother,
                                             Int_t detElemId,
                                             Int_t manuId)
: AliMUONVPainter("PAD"),
fDetElemId(detElemId),
fManuId(manuId),
fLineColorBck(-1),
fLineWidthBck(-1),
fFillColorBck(-1),
fFillStyleBck(-1)
{
  /// ctor
  SetID(detElemId,manuId);
  SetName(Form("%s/PAD",mother.GetName()));
  SetPathName(Form("%s/PAD",mother.PathName().Data()));
  SetContour(mother.Contour());
}

//_____________________________________________________________________________
AliMUONManuPadPainter::~AliMUONManuPadPainter()
{
  /// dtor
}

//_____________________________________________________________________________
void
AliMUONManuPadPainter::BackupStyle()
{
  /// Remember line and fill style values
  
  fFillStyleBck = gVirtualX->GetFillStyle();
  fFillColorBck = gVirtualX->GetFillColor();
  fLineColorBck = gVirtualX->GetLineColor();
  fLineWidthBck = gVirtualX->GetLineWidth();
}

//_____________________________________________________________________________
void
AliMUONManuPadPainter::ComputeDataRange(const AliMUONVTrackerData& data,
                                        Int_t dataIndex,
                                        Double_t& dataMin, Double_t& dataMax) const
{
  /// Compute data range spanned by this manu pads
  const AliMpVSegmentation* seg = AliMpSegmentation::Instance()->GetMpSegmentationByElectronics(fDetElemId,fManuId);
  
  dataMin = FLT_MAX;
  dataMax = -FLT_MAX;
  
  for ( Int_t manuChannel = 0; manuChannel < AliMpConstants::ManuNofChannels(); 
        ++manuChannel )
  {
    if ( seg->HasPadByLocation(fManuId,manuChannel) )
    {
      Double_t value = data.Channel(fDetElemId, fManuId, manuChannel, dataIndex);
      dataMin = TMath::Min(value,dataMin);
      dataMax = TMath::Max(value,dataMax);
    }
  }
}

//_____________________________________________________________________________
char* 
AliMUONManuPadPainter::GetObjectInfo(Int_t px, Int_t py) const
{
  /// Get object info
  return const_cast<char*>(Form("%s:%d:%d",GetName(),px,py));
}

//_____________________________________________________________________________
TString
AliMUONManuPadPainter::NameAtPosition(Double_t x, Double_t y) const
{
  /// Specific name, dependent on the position within this painter
  
  TString name("invalid");

  AliMpPad pad = PadByPosition(x,y);
  
  if ( pad.IsValid() )
  {  
    name = Form("%s%d",PathName().Data(),pad.GetManuChannel());
  }
  
  return name;
}

//_____________________________________________________________________________
Bool_t
AliMUONManuPadPainter::IsIncluded() const
{
  /// whether this manu is included in the readout or not
  return ( InteractiveReadOutConfig()->Manu(fDetElemId,fManuId) > 0 );
}

//_____________________________________________________________________________
TString 
AliMUONManuPadPainter::Describe(const AliMUONVTrackerData& data, Int_t dataIndex,
                                Double_t x, Double_t y)
{
  /// Describe data at given location
  
  if ( ! data.HasManu(fDetElemId,fManuId) ) return "";

  AliMpPad pad = PadByPosition(x,y);
  
  if ( pad.IsValid() ) 
  {
    Double_t value = data.Channel(fDetElemId,fManuId,pad.GetManuChannel(),dataIndex);
  
    return AliMUONPainterHelper::Instance()->FormatValue(data.DimensionName(dataIndex).Data(),value);
  }
  else
  {
    return "invalid";
  }
}

//_____________________________________________________________________________
void 
AliMUONManuPadPainter::DrawHistogramClone(Double_t* values) const
{
  /// Draw histogram for pad at (values[0],values[1])
  
  if ( !values ) return;

  AliMUONPainterGroup* group = Master()->PlotterGroup();
  
  if ( !group ) return; // no data to histogram in this painter
    
  AliMpPad pad = PadByPosition(values[0],values[1]);
  
  AliMUONVTrackerData* data = group->Data();
  
  AliMUONTrackerDataHistogrammer tdh(*data,0,-1);

  fHistogram = tdh.CreateChannelHisto(fDetElemId, fManuId,pad.GetManuChannel());

  if (fHistogram) 
  {
    new TCanvas();
    fHistogram->Draw();
  }  
}

//_____________________________________________________________________________
void
AliMUONManuPadPainter::PaintArea(const AliMUONVTrackerData& data,
                                 Int_t dataIndex,
                                 Double_t min,
                                 Double_t max)
{
    /// Paint area of this manu pads according to the data
  
  if ( !gPad ) return;
  
  if ( ! data.HasManu(fDetElemId,fManuId) ) return;
  
  AliMUONPainterHelper* h = AliMUONPainterHelper::Instance();
  
  BackupStyle();
    
  gVirtualX->SetLineColor(-1);
  gVirtualX->SetFillStyle(1001);
  
  const AliMpVSegmentation* seg = AliMpSegmentation::Instance()->GetMpSegmentationByElectronics(fDetElemId,fManuId);
  
  for ( Int_t i = 0; i < AliMpConstants::ManuNofChannels(); ++i ) 
  {    
    AliMpPad pad = seg->PadByLocation(fManuId,i,kFALSE);
    
    if ( pad.IsValid() )
    {
      Double_t value = data.Channel(fDetElemId,fManuId,i,dataIndex);
      
      if ( value >= AliMUONVCalibParam::InvalidFloatValue() ) continue;
      
      Int_t color = h->ColorFromValue(value,min,max);
      
      gVirtualX->SetFillColor(color);
      
      PaintPad(pad);
      
    }
  }
  
  RestoreStyle();
}
                      
//_____________________________________________________________________________
void
AliMUONManuPadPainter::PaintPad(const AliMpPad& pad) const
{
  Double_t blx = pad.GetPositionX()-pad.GetDimensionX();
  Double_t bly = pad.GetPositionY()-pad.GetDimensionY();
  
  Double_t urx = pad.GetPositionX()+pad.GetDimensionX();
  Double_t ury = pad.GetPositionY()+pad.GetDimensionY();
  
  Double_t xe1,ye1,xe2,ye2,z;

  AliMUONPainterHelper::Instance()->Local2Global(fDetElemId,blx,bly,0,xe1,ye1,z);
  AliMUONPainterHelper::Instance()->Local2Global(fDetElemId,urx,ury,0,xe2,ye2,z);
  
  gPad->PaintBox(xe1,ye1,xe2,ye2);
}

//_____________________________________________________________________________
void
AliMUONManuPadPainter::PaintOutline(Int_t color, Int_t, Double_t x, Double_t y)
{
  /// Paint the outline of our pads

  if ( !gPad ) return;
  
  Int_t lineColor = color >= 0 ? color : GetLineColor();
  
  if ( lineColor > 0 )
  {
    BackupStyle();
    
    gVirtualX->SetLineColor(lineColor);
    gVirtualX->SetFillStyle(0);
    
    if ( x < FLT_MAX && y < FLT_MAX ) 
    {
      // find pad to be drawn
      AliMpPad pad = PadByPosition(x,y);
      
      PaintPad(pad);        
    }
    else
    {
      const AliMpVSegmentation* seg = AliMpSegmentation::Instance()->GetMpSegmentationByElectronics(fDetElemId,fManuId);

      for ( Int_t i = 0; i < AliMpConstants::ManuNofChannels(); ++i ) 
      {    
        AliMpPad pad = seg->PadByLocation(fManuId,i,kFALSE);
        
        if (pad.IsValid()) PaintPad(pad);
      
        PaintPad(pad);        
      }
    }
    RestoreStyle();
  }  
  
}

//_____________________________________________________________________________
AliMpPad
AliMUONManuPadPainter::PadByPosition(Double_t x, Double_t y) const
{
  /// Find the pad at given exploded-position (x,y)
  
  const AliMpVSegmentation* seg = AliMpSegmentation::Instance()->GetMpSegmentationByElectronics(fDetElemId,fManuId);
  
  Double_t xg,yg,zg;
  Double_t xl,yl,zl;
  
  AliMUONPainterHelper::Instance()->Local2Global(fDetElemId,0.0,0.0,0.0,xg,yg,zg); // to get zg
  
  AliMUONPainterHelper::Instance()->Global2Local(fDetElemId,x,y,zg,xl,yl,zl);
  
  return seg->PadByPosition(xl,yl);
}

//_____________________________________________________________________________
void
AliMUONManuPadPainter::RestoreStyle()
{
  /// Restore line and fill style values
  
  gVirtualX->SetFillStyle(fFillStyleBck);
  gVirtualX->SetFillColor(fFillColorBck);
  gVirtualX->SetLineColor(fLineColorBck);
  gVirtualX->SetLineWidth(fLineWidthBck);
}


