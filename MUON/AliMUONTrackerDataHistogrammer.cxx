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

#include "AliLog.h"
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
#include <TClass.h>
#include <TH1.h>
#include <TObjArray.h>
#include <TROOT.h>
#include <TMath.h>

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
                                                               Int_t externalDim,
                                                               Int_t internalDim)
: TObject(),
fkData(data),
fExternalDim(externalDim),
fInternalDim(internalDim)
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
  
  if ( fkData.HasBusPatch(busPatchId ) )
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
  
  if ( fkData.HasDetectionElement(detElemId) )
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
  
  if ( fkData.HasManu(detElemId,manuId) )
  {
    if ( fkData.IsChannelLevelEnabled() )
    {
      for ( Int_t i = 0; i < AliMpConstants::ManuNofChannels(); ++i ) 
      {
        if ( fkData.HasChannel(detElemId,manuId,i) )
        {
          if ( IsInternalMode() ) 
          {
            h.Fill(fkData.Channel(detElemId,manuId,i,fInternalDim));
          }
          else
          {
            AliMUONSparseHisto* sh = fkData.GetChannelSparseHisto(detElemId,manuId,i,fExternalDim);
            
            if ( sh ) 
            {       
              Add(h,*sh);
            }
          }
        }
      }
    }
    else
    {
      if ( IsInternalMode() ) 
      {
        h.Fill(fkData.Manu(detElemId,manuId,fInternalDim));
      }
      else
      {
        AliMUONSparseHisto* sh = fkData.GetManuSparseHisto(detElemId,manuId,fExternalDim);
        if (sh)
        {
          Add(h,*sh);
        }
      }
    }
  }
}

//_____________________________________________________________________________
TH1*
AliMUONTrackerDataHistogrammer::CreateChannelHisto(Int_t detElemId, 
                                                   Int_t manuId, 
                                                   Int_t manuChannel) const
{
  /// Create histogram of a given channel. Note that in order
  /// to keep memory footprint as low as possible, you should delete
  /// the returned pointer as soon as possible...
  
  if ( fkData.HasChannel(detElemId, manuId, manuChannel) && fkData.IsHistogrammed(fExternalDim) )
  {
    AliMUONSparseHisto* sh = fkData.GetChannelSparseHisto(detElemId,manuId,manuChannel);
    
    if ( sh ) 
    {
      Int_t nbins((1<<sh->Nbits()));
      Double_t xmin,xmax;
      fkData.HistogramRange(xmin,xmax);
      
      TH1* h = CreateHisto(Form("DE%04dMANU%04dCH%02d",detElemId,manuId,manuChannel),
                           nbins,xmin,xmax);
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
AliMUONTrackerDataHistogrammer::CreateHisto(const char* name,
                                            Int_t nbins,
                                            Double_t xmin,
                                            Double_t xmax) const
{
  /// Create a single histogram
  
  TH1* h(0);
  
  if ( xmin < xmax ) 
  {
    h = new TH1F(name,name,nbins,xmin,xmax);
    h->SetDirectory(gROOT);
  }
	else
	{
		AliError(Form("Cannot create histo for name=%s nbins=%d xmin=%e xmax=%e",name,nbins,xmin,xmax));
	}
  return h;
}

//_____________________________________________________________________________
TH1* 
AliMUONTrackerDataHistogrammer::CreateHisto(const AliMUONVPainter& painter, 
																						Int_t externalDim,
																						Int_t internalDim)
{
  /// Create an histogram, from given dim of given data, 
  /// for all the channels handled by painter

  AliMUONPainterGroup* group = painter.Master()->PlotterGroup();
  
  if ( !group ) return 0x0; // no data to histogram in this painter
  
  AliMUONVTrackerData* data = group->Data();
  
  if ( externalDim >= data->ExternalDimension() )
  {
    AliErrorClass(Form("externalDim %d is out of bounds",externalDim));
    return 0x0;
  }

  if ( internalDim >= data->NumberOfDimensions() )
  {
    AliErrorClass(Form("internalDim %d is out of bounds",internalDim));
    return 0x0;
  }
  
  if ( internalDim < 0 && externalDim < 0 ) 
  {
    AliErrorClass("Both internal and external dim are < 0 !!!");
    return 0x0;
  }
  
  AliMUONTrackerDataHistogrammer tdh(*data,externalDim,internalDim);

  TObjArray manuArray;
  
  painter.FillManuList(manuArray);

  AliMpManuUID* mid;
  TIter next(&manuArray);

  TString basename(Form("%s-%s",painter.PathName().Data(),painter.Attributes().GetName()));
  TString ext;
  Int_t nbins((1<<12));
  Double_t xmin(0.0);
  Double_t xmax(0.0);
  
  if ( !tdh.IsInternalMode() ) 
  {
    data->HistogramRange(xmin,xmax);
    
    xmin -= 0.5;
    xmax -= 0.5;
    
    ext = data->ExternalDimensionName(externalDim).Data();
  }
  else
  {
    tdh.GetDataRange(manuArray,xmin,xmax);
    ext = data->DimensionName(internalDim).Data();
    nbins = 100;
  }
  
  TString name(Form("%s-%s",basename.Data(),ext.Data()));

  TH1* histo = tdh.CreateHisto(name.Data(),nbins,xmin,xmax);

  if ( histo ) 
  {
    while ( ( mid = static_cast<AliMpManuUID*>(next()) ) )
    {
      TH1* h = tdh.CreateManuHisto(mid->DetElemId(),mid->ManuId(),nbins,xmin,xmax);
      if ( h ) 
      {
        histo->Add(h);
      }
      delete h;
    }
  }
  else
  {
    AliErrorClass(Form("Could not create histo for painter %s (%p) data %s (%p) external dim %d internal dim %d",
                       painter.PathName().Data(),&painter,
                       data->GetName(),data,externalDim,internalDim));
  }
  
  if (histo) histo->SetDirectory(gROOT);
  
  return histo;
}

//_____________________________________________________________________________
TH1*
AliMUONTrackerDataHistogrammer::CreateManuHisto(Int_t detElemId, Int_t manuId,
                                                Int_t nbins,
                                                Double_t xmin,
                                                Double_t xmax) const
{
  /// Create histogram of a given manu. Note that in order
  /// to keep memory footprint as low as possible, you should delete
  /// the returned pointer as soon as possible...
  
  TH1* h(0x0);
  
  if ( !fkData.HasManu(detElemId,manuId) ) return 0x0;
  
  if ( ( fExternalDim >= 0 && fkData.IsHistogrammed(fExternalDim) ) ||
       ( fInternalDim >= 0 && fInternalDim < fkData.NumberOfDimensions() ) )
  {
    h = CreateHisto(Form("DE%04dMANU%04d",detElemId,manuId),
                    nbins,xmin,xmax);
    if ( h ) AddManuHisto(*h,detElemId,manuId);
  }
  
  return h;
}

//_____________________________________________________________________________
void
AliMUONTrackerDataHistogrammer::GetDataRange(const TObjArray& manuArray, 
                                             Double_t& xmin, Double_t& xmax) const
{
  /// Get data range (in case of InternalMode() only) spanned by the manus in
  /// manuArray
  
  xmin = FLT_MAX;
  xmax = -FLT_MAX;
  
  if (!IsInternalMode())
  {
    AliError("Cannot use this method for external mode !");
  }

  AliMpManuUID* mid;
  TIter next(&manuArray);
  
  while ( ( mid = static_cast<AliMpManuUID*>(next()) ) )
  {
    Int_t detElemId = mid->DetElemId();
    Int_t manuId = mid->ManuId();
    
    for ( Int_t i = 0; i < AliMpConstants::ManuNofChannels(); ++i ) 
    {
      if ( fkData.HasChannel(detElemId,manuId,i) ) 
      {
        Double_t value = fkData.Channel(detElemId,manuId,i,fInternalDim);
				
				if ( ! TMath::Finite(value) )
				{
					AliError(Form("Got a NaN for DE %d manu %d ch %d",detElemId,manuId,i));
				}
				else
				{
					xmin = TMath::Min(xmin,value);
					xmax = TMath::Max(xmax,value);
				}
      }
    }
  }

}

