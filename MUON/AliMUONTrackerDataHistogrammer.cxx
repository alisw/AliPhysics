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

#include "AliMUONTrackerDataHistogrammer.h"

#include "AliMUONPainterGroup.h"
#include "AliMUONSparseHisto.h"
#include "AliMUONVPainter.h"
#include "AliMUONVTrackerData.h"
#include "AliMpBusPatch.h"
#include "AliMpConstants.h"
#include "AliMpDDLStore.h"
#include "AliMpDEIterator.h"
#include "AliMpDetElement.h"
#include "AliMpManuUID.h"
#include <TH1.h>
#include <TObjArray.h>

///\class AliMUONTrackerDataHistogrammer
///
/// Class to generate histograms from AliMUONVTrackerData 
/// (and AliMUONVPainter) objects
///
/// \author Laurent Aphecetche, Subatech
///

///\cond CLASSIMP
ClassImp(AliMUONTrackerDataHistogrammer)
///\endcond CLASSIMP

//_____________________________________________________________________________
AliMUONTrackerDataHistogrammer::AliMUONTrackerDataHistogrammer(const AliMUONVTrackerData& data,
                                                               Int_t dim)
: TObject(),
fData(data),
fDim(dim)
{
  /// ctor
}

//_____________________________________________________________________________
AliMUONTrackerDataHistogrammer::~AliMUONTrackerDataHistogrammer()
{
  /// dtor
}

//_____________________________________________________________________________
void
AliMUONTrackerDataHistogrammer::Add(TH1& h, const AliMUONSparseHisto& sh) const
{
  /// Add sparse histo content to histogram.
  
  Double_t entries(h.GetEntries());
  
  for ( Int_t i = 0; i < sh.GetNbins(); ++i ) 
  {
    Int_t count = sh.GetBinContent(i);
    
    h.Fill(sh.GetBinCenter(i),count);
    
    entries += count;
  }
  
  h.SetEntries(entries);
  
  if (sh.HasUnderflow()) h.SetBinContent(0,1);
  if (sh.HasOverflow()) h.SetBinContent(h.GetNbinsX()+1,1);
}

//_____________________________________________________________________________
void
AliMUONTrackerDataHistogrammer::AddBusPatchHisto(TH1& h, Int_t busPatchId) const
{
  /// Add data from one bus patch to the histogram
  
  if ( fData.HasBusPatch(busPatchId ) )
  {
    AliMpBusPatch* busPatch = AliMpDDLStore::Instance()->GetBusPatch(busPatchId);
    for ( Int_t i = 0; i < busPatch->GetNofManus(); ++i ) 
    {
      Int_t manuId = busPatch->GetManuId(i);
      AddManuHisto(h,busPatch->GetDEId(),manuId);
    }
  }
}
//_____________________________________________________________________________
void
AliMUONTrackerDataHistogrammer::AddDEHisto(TH1& h, Int_t detElemId) const
{
  /// Add data from one detection element to the histogram
  
  if ( fData.HasDetectionElement(detElemId) )
  {
    AliMpDetElement* de = AliMpDDLStore::Instance()->GetDetElement(detElemId);
    for ( Int_t i = 0; i < de->GetNofBusPatches(); ++ i ) 
    {
      Int_t busPatchId = de->GetBusPatchId(i);
      AddBusPatchHisto(h,busPatchId);
    }
  }
}

//_____________________________________________________________________________
void
AliMUONTrackerDataHistogrammer::AddManuHisto(TH1& h, Int_t detElemId, Int_t manuId) const
{
  /// Add data from a given manu to histogram
  
  if ( fData.HasManu(detElemId,manuId) )
  {
    for ( Int_t i = 0; i < AliMpConstants::ManuNofChannels(); ++i ) 
    {
      if ( fData.HasChannel(detElemId,manuId,i) )
      {
        AliMUONSparseHisto* sh = fData.GetChannelSparseHisto(detElemId,manuId,i);
        
        if ( sh ) 
        {      
          Add(h,*sh);
        }
      }
    }
  }
}


//_____________________________________________________________________________
TH1*
AliMUONTrackerDataHistogrammer::CreateBusPatchHisto(Int_t busPatchId) const
{
  /// Create histogram of a given bus patch. Note that in order
  /// to keep memory footprint as low as possible, you should delete
  /// the returned pointer as soon as possible...
  
  TH1* h(0x0);
  
  if ( fData.HasBusPatch(busPatchId) && fData.IsHistogrammed(fDim)) 
  {
    h = CreateHisto(Form("BP%04d_%d",busPatchId,fDim));
    if ( h ) AddBusPatchHisto(*h,busPatchId);
  }
  
  return h;
}  


//_____________________________________________________________________________
TH1*
AliMUONTrackerDataHistogrammer::CreateChamberHisto(Int_t chamberId) const
{
  /// Create histogram of a given chamber. Note that in order
  /// to keep memory footprint as low as possible, you should delete
  /// the returned pointer as soon as possible...
  
  TH1* h(0x0);
  
  if ( fData.HasChamber(chamberId) && fData.IsHistogrammed(fDim))
  {
    h = CreateHisto(Form("CHAMBER%02d_%d",chamberId,fDim));
    if ( h ) 
    {
      AliMpDEIterator it;
      it.First(chamberId);
      while ( !it.IsDone() )
      {
        Int_t detElemId = it.CurrentDEId();
        AddDEHisto(*h,detElemId);
        it.Next();
      }
    }
  }
  
  return h;
}

//_____________________________________________________________________________
TH1*
AliMUONTrackerDataHistogrammer::CreateChannelHisto(Int_t detElemId, Int_t manuId, 
                                                   Int_t manuChannel) const
{
  /// Create histogram of a given channel. Note that in order
  /// to keep memory footprint as low as possible, you should delete
  /// the returned pointer as soon as possible...
  
  if ( fData.HasChannel(detElemId, manuId, manuChannel) && fData.IsHistogrammed(fDim) )
  {
    AliMUONSparseHisto* sh = fData.GetChannelSparseHisto(detElemId,manuId,manuChannel);
    
    if ( sh ) 
    {
      TH1* h = CreateHisto(Form("DE%04dMANU%04dCH%02d_%d",detElemId,manuId,manuChannel,fDim));
      if (h ) 
      {
        Add(*h,*sh);
      }
      return h;
    }
  }
  return 0x0;
}

//_____________________________________________________________________________
TH1*
AliMUONTrackerDataHistogrammer::CreateDEHisto(Int_t detElemId) const
{
  /// Create histogram of a given detection element. Note that in order
  /// to keep memory footprint as low as possible, you should delete
  /// the returned pointer as soon as possible...
  
  TH1* h(0x0);
  
  if ( fData.HasDetectionElement(detElemId) && fData.IsHistogrammed(fDim) ) 
  {
    h = CreateHisto(Form("DE%04d-%d",detElemId,fDim));
    if (h) AddDEHisto(*h,detElemId);
  }
  
  return h;
}

//_____________________________________________________________________________
TH1* 
AliMUONTrackerDataHistogrammer::CreateHisto(const AliMUONVPainter& painter)
{
  /// Create an histogram, from given dim of given data, 
  /// for all the channels handled by painter

  AliMUONPainterGroup* group = painter.Master()->PlotterGroup();
  
  if ( !group ) return 0x0; // no data to histogram in this painter
  
  AliMUONVTrackerData* data = group->Data();
  Int_t dim = data->InternalToExternal(group->DataIndex());
  
  AliMUONTrackerDataHistogrammer tdh(*data,dim);

  TObjArray manuArray;
  
  painter.FillManuList(manuArray);
  
  AliMpManuUID* mid;
  TIter next(&manuArray);
  
  TH1* histo = tdh.CreateHisto(Form("%s-%s",painter.PathName().Data(),painter.Attributes().GetName()));
  
  while ( ( mid = static_cast<AliMpManuUID*>(next()) ) )
  {
    TH1* h = tdh.CreateManuHisto(mid->DetElemId(),mid->ManuId());
    if ( h ) 
    {
      histo->Add(h);
    }
    delete h;
  }
  
  return histo;
}

//_____________________________________________________________________________
TH1*
AliMUONTrackerDataHistogrammer::CreateHisto(const char* name) const
{
  /// Create a single histogram

  Double_t xmin, xmax;
  
  fData.HistogramRange(xmin,xmax);
  
  TH1* h(0x0);
  
  if ( xmin != xmax ) 
  {
    h = new TH1I(name,Form("Data=%s Dim=%s",
                              fData.GetName(),
                              fData.ExternalDimensionName(fDim).Data()),
                    (1<<12),
                    xmin-0.5,
                    xmax-0.5);
    h->SetDirectory(0);
  }
  return h;
}


//_____________________________________________________________________________
TH1*
AliMUONTrackerDataHistogrammer::CreateManuHisto(Int_t detElemId, Int_t manuId) const
{
  /// Create histogram of a given manu. Note that in order
  /// to keep memory footprint as low as possible, you should delete
  /// the returned pointer as soon as possible...
  
  TH1* h(0x0);
  
  if ( fData.HasManu(detElemId, manuId) && fData.IsHistogrammed(fDim) ) 
  {
    h = CreateHisto(Form("DE%04dMANU%04d_%d",detElemId,manuId,fDim));
    if ( h ) AddManuHisto(*h,detElemId,manuId);
  }
  
  return h;
}

