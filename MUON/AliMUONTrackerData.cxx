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

#include "AliMUONTrackerData.h"

#include "AliCodeTimer.h"
#include "AliLog.h"
#include "AliMUON1DArray.h"
#include "AliMUON1DMap.h"
#include "AliMUON2DMap.h"
#include "AliMUONCalibParamND.h"
#include "AliMUONSparseHisto.h"
#include "AliMUONVStore.h"
#include "AliMpBusPatch.h"
#include "AliMpConstants.h"
#include "AliMpDDLStore.h"
#include "AliMpDEIterator.h"
#include "AliMpDEManager.h"
#include "AliMpDetElement.h"
#include "AliMpHVNamer.h"
#include "AliMpManuIterator.h"
#include <Riostream.h>
#include <TH1.h>
#include <TMath.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <TString.h>
#include <TVector2.h>
#include <float.h>

/// \class AliMUONTrackerData
///
/// Implementation of AliMUONVTrackerData class
///
/// \author Laurent Aphecetche, Subatech

///\cond CLASSIMP
ClassImp(AliMUONTrackerData)
///\endcond

const Int_t AliMUONTrackerData::fgkExtraDimension = 2;
const Int_t AliMUONTrackerData::fgkVirtualExtraDimension = 1;

//_____________________________________________________________________________
AliMUONTrackerData::AliMUONTrackerData(const char* name, const char* title,
                                       Int_t dimension,
                                       Bool_t runnable)
: AliMUONVTrackerData(name,title),
fChannelValues(0x0),
fManuValues(0x0),
fBusPatchValues(0x0),
fDEValues(0x0),
fChamberValues(0x0),
fPCBValues(0x0),
fDimension(dimension*2+fgkExtraDimension),
fNevents(0x0),
fDimensionNames(new TObjArray(fDimension+fgkVirtualExtraDimension)),
fExternalDimensionNames(new TObjArray(dimension)),
fExternalDimension(dimension),
fIsRunnable(runnable),
fHistogramming(new Int_t[fExternalDimension]),
fChannelHistos(0x0)
{  
  /// ctor
  memset(fHistogramming,0,sizeof(Int_t)); // histogramming is off by default. Use SetHistogramDimension to turn it on.
  fExternalDimensionNames->SetOwner(kTRUE);
  fDimensionNames->SetOwner(kTRUE);  
  fDimensionNames->AddAt(new TObjString("occ"),IndexOfOccupancyDimension());
  fDimensionNames->AddAt(new TObjString("N"),IndexOfNumberDimension());
  fDimensionNames->AddAt(new TObjString("n"),NumberOfDimensions()-fgkVirtualExtraDimension);
  Clear();
}

//_____________________________________________________________________________
AliMUONTrackerData::~AliMUONTrackerData()
{
  /// dtor
  delete fChannelValues;
  delete fManuValues;
  delete fBusPatchValues;
  delete fDEValues;
  delete fChamberValues;
  delete fPCBValues;
  delete fDimensionNames;
  delete fExternalDimensionNames;
  delete[] fHistogramming;
  delete fChannelHistos;
}

//_____________________________________________________________________________
Bool_t
AliMUONTrackerData::Add(const AliMUONVStore& store)
{
  /// Add the given external store to our internal store
  
  AliCodeTimerAuto(GetName());
    
  ++fNevents;
  NumberOfEventsChanged();
  
  if (!fChannelValues)
  {
    fChannelValues = store.Create();
    fManuValues = store.Create();
    fPCBValues = store.Create();
    fBusPatchValues = new AliMUON1DMap;
    fDEValues = new AliMUON1DMap;
    fChamberValues = new AliMUON1DArray;
  }
  
  TIter next(store.CreateIterator());
  AliMUONVCalibParam* external;
  
  while ( ( external = static_cast<AliMUONVCalibParam*>(next()) ) )
  {
    if ( external->Dimension() != ExternalDimension() )
    {
      AliError(Form("Incompatible dimensions %d vs %d",
                    external->Dimension(),ExternalDimension()));
      return kFALSE;
    }
    

    AliMUONVCalibParam* chamber, *de, *busPatch, *pcb, *manu, *channel;
    AliMpDetElement* mpde;
    
    Int_t manuId = GetParts(external,chamber,de,busPatch,pcb,manu,channel,mpde);
    
    Int_t detElemId = mpde->GetId();
    
    Double_t value[] = { 0.0, 0.0 };
    
    Int_t nch = mpde->NofChannelsInManu(manuId);
    
    for ( Int_t i = 0; i < external->Size(); ++i ) 
    {
      Bool_t existingChannel =  ( nch == AliMpConstants::ManuNofChannels() ? kTRUE
                                                                           : mpde->IsConnectedChannel(manuId,i));
      // note we only use IsConnectedChannel method when really needed, as
      // it costs (some) CPU time...
      
      if ( existingChannel ) 
      {
        Bool_t validChannel(kFALSE);
        
        for ( Int_t j = 0; j < external->Dimension(); ++j )
        {
          Double_t vext = external->IsDoublePrecision() ? 
            external->ValueAsDoubleFast(i,j) :
            external->ValueAsFloatFast(i,j);
          
          if ( vext >= AliMUONVCalibParam::InvalidFloatValue() ) continue;
          
          validChannel = kTRUE;
                    
          Int_t ix = External2Internal(j);
          
          value[0] = vext;
          value[1] = vext*vext;
          
          if ( IsHistogrammed(j) )
          {
            FillChannel(detElemId,manuId,i,j,vext);
          }
          
          for ( Int_t k = 0; k < 2; ++k ) 
          {
            channel->SetValueAsDoubleFast(i,ix+k,channel->ValueAsDoubleFast(i,ix+k)+value[k]);

            manu->SetValueAsDoubleFast(0,ix+k,manu->ValueAsDoubleFast(0,ix+k)+value[k]);            
            
            busPatch->SetValueAsDoubleFast(0,ix+k,busPatch->ValueAsDoubleFast(0,ix+k)+value[k]);
          
            de->SetValueAsDoubleFast(0,ix+k,de->ValueAsDoubleFast(0,ix+k)+value[k]);
          
            chamber->SetValueAsDoubleFast(0,ix+k,chamber->ValueAsDoubleFast(0,ix+k)+value[k]);
          
            if ( pcb ) 
            {
              pcb->SetValueAsDoubleFast(0,ix+k,pcb->ValueAsDoubleFast(0,ix+k)+value[k]);
            }
          }
        }
        
        if ( validChannel )
        {
          channel->SetValueAsDoubleFast(i,IndexOfOccupancyDimension(),
                                        channel->ValueAsDoubleFast(i,IndexOfOccupancyDimension())+1.0);
          manu->SetValueAsDoubleFast(0,IndexOfOccupancyDimension(),
                                                 manu->ValueAsDoubleFast(0,IndexOfOccupancyDimension())+1.0);        
          busPatch->SetValueAsDoubleFast(0,IndexOfOccupancyDimension(),
                                                         busPatch->ValueAsDoubleFast(0,IndexOfOccupancyDimension())+1.0);        
          de->SetValueAsDoubleFast(0,IndexOfOccupancyDimension(),
                                             de->ValueAsDoubleFast(0,IndexOfOccupancyDimension())+1.0);        
          chamber->SetValueAsDoubleFast(0,IndexOfOccupancyDimension(),
                                                       chamber->ValueAsDoubleFast(0,IndexOfOccupancyDimension())+1.0); 
          if ( pcb ) 
          {
            pcb->SetValueAsDoubleFast(0,IndexOfOccupancyDimension(),
                                      pcb->ValueAsDoubleFast(0,IndexOfOccupancyDimension())+1.0);        
          }
        }
      }
    }
  }
  
  return kTRUE;
}

//_____________________________________________________________________________
Double_t 
AliMUONTrackerData::BusPatch(Int_t busPatchId, Int_t dim) const
{
  /// Return the value of a given buspatch for a given dimension
  /// or 0 if not existing
  AliMUONVCalibParam* param = BusPatchParam(busPatchId);
  return param ? Value(*param,0,dim) : 0.0;
}

//_____________________________________________________________________________
AliMUONVCalibParam* 
AliMUONTrackerData::BusPatchParam(Int_t busPatchId, Bool_t create) const
{
  /// Return (if it exist), the VCalibParam for a given busPatch
  
  AliMUONVCalibParam* busPatch = fBusPatchValues ? static_cast<AliMUONVCalibParam*>
    (fBusPatchValues->FindObject(busPatchId)) : 0x0;
  
  if (!busPatch && create && fBusPatchValues)
  {
    busPatch = CreateBusPatchParam(busPatchId);
    fBusPatchValues->Add(busPatch);
  }
  
  return busPatch;
}

//_____________________________________________________________________________
AliMUONVCalibParam* 
AliMUONTrackerData::CreateBusPatchParam(Int_t busPatchId) const
{
  /// Create storage for one bus patch
  
  AliCodeTimerAuto("");
  
  AliMpBusPatch* bp = AliMpDDLStore::Instance()->GetBusPatch(busPatchId);
  
  if (!bp)
  {
    AliError(Form("Got an invalid buspatchId = %d",busPatchId));
    return 0x0;
  }
  
  AliMUONVCalibParam* busPatch = new AliMUONCalibParamND(Dimension(),1,busPatchId,0,0.0);
  
  // set the number of channels in that buspatch
  
  Int_t nchannels(0);
  
  Int_t detElemId = AliMpDDLStore::Instance()->GetDEfromBus(busPatchId);
  
  AliMpDetElement* de = AliMpDDLStore::Instance()->GetDetElement(detElemId);
  
  for ( Int_t i = 0; i < bp->GetNofManus(); ++i ) 
  {
    Int_t manuId = bp->GetManuId(i);
    nchannels += de->NofChannelsInManu(manuId);
  }
  
  busPatch->SetValueAsDouble(0,IndexOfNumberDimension(),nchannels);
  
  return busPatch;
}

//_____________________________________________________________________________
Double_t 
AliMUONTrackerData::Chamber(Int_t chamberId, Int_t dim) const
{
  /// Return the value fo a given chamber for a given dimension,
  /// or zero if not existing
  AliMUONVCalibParam* param = ChamberParam(chamberId);
  return param ? Value(*param,0,dim) : 0.0;
}

//_____________________________________________________________________________
AliMUONVCalibParam* 
AliMUONTrackerData::ChamberParam(Int_t chamberId, Bool_t create) const
{
  /// Return (if it exist) the VCalibParam for a given chamber
  
  AliMUONVCalibParam* chamber =  fChamberValues ? static_cast<AliMUONVCalibParam*>
  (fChamberValues->FindObject(chamberId)) : 0x0;
  
  if (!chamber && create && fChamberValues)
  {
    chamber = CreateChamberParam(chamberId);
    fChamberValues->Add(chamber);
  }
    
  return chamber;
}

//_____________________________________________________________________________
AliMUONVCalibParam* 
AliMUONTrackerData::CreateChamberParam(Int_t chamberId) const
{
  /// Create storage for one chamber
  
  AliCodeTimerAuto("");
  
  AliMUONVCalibParam* chamber = new AliMUONCalibParamND(Dimension(),1,chamberId,0,0.0);
  
  // set the number of channels in that chamber
  
  Int_t nchannels(0);
  
  AliMpDEIterator it;
  
  it.First(chamberId);
  
  while ( !it.IsDone() )
  {        
    AliMpDetElement* det = it.CurrentDE();
    
    for ( Int_t i = 0; i < det->GetNofBusPatches(); ++i ) 
    {
      Int_t busPatchId = det->GetBusPatchId(i);
      AliMpBusPatch* bp = AliMpDDLStore::Instance()->GetBusPatch(busPatchId);
      for ( Int_t j = 0; j < bp->GetNofManus(); ++j ) 
      {
        Int_t manuId = bp->GetManuId(j);
        nchannels += det->NofChannelsInManu(manuId);
      }        
    }
    
    it.Next();
  }
  
  chamber->SetValueAsDouble(0,IndexOfNumberDimension(),nchannels);
  
  return chamber;
}

//_____________________________________________________________________________
Double_t 
AliMUONTrackerData::Channel(Int_t detElemId, Int_t manuId, 
                            Int_t manuChannel, Int_t dim) const
{
  /// Return the value for a given channel for a given dimension
  
  AliMUONVCalibParam* param = ChannelParam(detElemId,manuId);
  
  return param ? Value(*param,manuChannel,dim) : 0.0;
}

//_____________________________________________________________________________
AliMUONVCalibParam* 
AliMUONTrackerData::ChannelParam(Int_t detElemId, Int_t manuId,
                                 AliMUONVCalibParam* external) const
{
  /// Return (if it exist) the VCalibParam for a given manu
  
  AliMUONVCalibParam* param = fChannelValues ? static_cast<AliMUONVCalibParam*>
    (fChannelValues->FindObject(detElemId,manuId)) : 0x0 ;
  
  if (!param && external && fChannelValues)
  {
    param = CreateDouble(*external);
    fChannelValues->Add(param);
  }
  
  return param;
}

//_____________________________________________________________________________
void 
AliMUONTrackerData::Clear(Option_t*)
{
  /// Clear all the values
  if ( fChannelValues ) fChannelValues->Clear();
  if ( fManuValues ) fManuValues->Clear();
  if ( fBusPatchValues) fBusPatchValues->Clear();
  if ( fPCBValues ) fPCBValues->Clear();
  if ( fDEValues) fDEValues->Clear();
  if ( fChamberValues ) fChamberValues->Clear();
  if ( fChannelHistos ) fChannelHistos->Clear();
  fNevents = 0;
  NumberOfEventsChanged();
}

//_____________________________________________________________________________
Double_t 
AliMUONTrackerData::Count(Int_t detElemId, Int_t manuId, 
                          Int_t manuChannel) const
{
  /// Return the number of times a given channel had data
  
  return Channel(detElemId,manuId,manuChannel,IndexOfNumberDimension());
}

//_____________________________________________________________________________
AliMUONVCalibParam*
AliMUONTrackerData::CreateDouble(const AliMUONVCalibParam& param) const
{
  /// Create a double version of VCalibParam, for internal use
  
  AliCodeTimerAuto("");
  
  AliMUONVCalibParam* c = new AliMUONCalibParamND(Dimension(),
                                                  param.Size(),
                                                  param.ID0(),
                                                  param.ID1(),
                                                  0.0);
  
  for ( Int_t i = 0; i < c->Size(); ++i ) 
  {
    c->SetValueAsDouble(i,IndexOfNumberDimension(),1.0);
  }
  
  return c;
}

//_____________________________________________________________________________
Double_t 
AliMUONTrackerData::DetectionElement(Int_t detElemId, Int_t dim) const
{
  /// Return the value for a given detection element for a given dimension
  AliMUONVCalibParam* param = DetectionElementParam(detElemId);
  return param ? Value(*param,0,dim) : 0.0;

}

//_____________________________________________________________________________
AliMUONVCalibParam* 
AliMUONTrackerData::DetectionElementParam(Int_t detElemId, Bool_t create) const
{
  /// Return (if it exist) the VCalibParam for a given detection element
  
  AliMUONVCalibParam* de = fDEValues ? static_cast<AliMUONVCalibParam*>
    (fDEValues->FindObject(detElemId)) : 0x0 ;
  
  if (!de && create && fDEValues)
  {
    de = CreateDetectionElementParam(detElemId);
    fDEValues->Add(de);
  }
  
  return de;
  
}

//_____________________________________________________________________________
AliMUONVCalibParam* 
AliMUONTrackerData::CreateDetectionElementParam(Int_t detElemId) const
{
  /// Create storage for one detection element
  
  AliCodeTimerAuto("");
  
  AliMUONVCalibParam*  de = new AliMUONCalibParamND(Dimension(),1,detElemId,0,0.0);
  
  AliMpDetElement* det = AliMpDDLStore::Instance()->GetDetElement(detElemId);
  Int_t nchannels(0);
  
  for ( Int_t i = 0; i < det->GetNofBusPatches(); ++i ) 
  {
    Int_t busPatchId = det->GetBusPatchId(i);
    AliMpBusPatch* bp = AliMpDDLStore::Instance()->GetBusPatch(busPatchId);
    for ( Int_t j = 0; j < bp->GetNofManus(); ++j ) 
    {
      Int_t manuId = bp->GetManuId(j);
      nchannels += det->NofChannelsInManu(manuId);
    }        
  }
  
  de->SetValueAsDouble(0,IndexOfNumberDimension(),nchannels);
  
  return de;
}

//_____________________________________________________________________________
TString 
AliMUONTrackerData::DimensionName(Int_t dim) const
{
  /// Get the name of a given dimension
  TObjString* value = static_cast<TObjString*>(fDimensionNames->At(dim));
  if ( value ) 
  {
    return value->String();
  }
  else
  {
    return TString("Invalid");
  }  
}

//_____________________________________________________________________________
TString 
AliMUONTrackerData::ExternalDimensionName(Int_t dim) const
{
  /// Get the name of a given external dimension
  
  TObjString* value = static_cast<TObjString*>(fExternalDimensionNames->At(dim));
  if ( value ) 
  {
    return value->String();
  }
  else
  {
    return TString("Invalid");
  }  
}

//_____________________________________________________________________________
void
AliMUONTrackerData::FillChannel(Int_t detElemId, Int_t manuId, Int_t manuChannel,
                                Int_t dim, Double_t value)
{
  /// Fill histogram of a given channel
  
  AliMUONSparseHisto* h = GetChannelHisto(detElemId, manuId, manuChannel,dim);
 
  h->Fill(static_cast<Int_t>(TMath::Nint(value)));
}

//_____________________________________________________________________________
TH1*
AliMUONTrackerData::CreateChannelHisto(Int_t detElemId, Int_t manuId, 
                                       Int_t manuChannel, Int_t dim)
{
  /// Create histogram of a given channel. Note that in order
  /// to keep memory footprint as low as possible, you should delete
  /// the returned pointer as soon as possible...

  if ( HasChannel(detElemId, manuId, manuChannel) && IsHistogrammed(dim) )
  {
    AliMUONSparseHisto* sh = GetChannelHisto(detElemId,manuId,manuChannel,dim);
  
    if ( sh ) 
    {
      TH1* h = new TH1I(Form("DE%04dMANU%04dCH%02d_%d",detElemId,manuId,manuChannel,dim),
                        Form("Data=%s Dim=%s",GetName(),ExternalDimensionName(dim).Data()),
                        4096,-0.5,4095.5);
      
      Add(*h,*sh);
      return h;
    }
  }
  return 0x0;
}

//_____________________________________________________________________________
void
AliMUONTrackerData::Add(TH1& h, const AliMUONSparseHisto& sh)
{
  /// Add sparse histo content to histogram.
  
  Double_t entries(h.GetEntries());
  
  for ( Int_t i = 0; i < sh.GetNbins(); ++i ) 
  {
    Int_t x = sh.GetBinContent(i);
    Int_t adc, count;
    sh.Decode(x,adc,count);
    h.Fill(adc,count);
    entries += count;
  }
  
  h.SetEntries(entries);
}

//_____________________________________________________________________________
TH1*
AliMUONTrackerData::CreateHisto(const char* name, Int_t dim) const
{
  /// Create a single histogram
 
  return new TH1I(name,Form("Data=%s Dim=%s",GetName(),ExternalDimensionName(dim).Data()),
                  4096,-0.5,4095.5);
}

//_____________________________________________________________________________
TH1*
AliMUONTrackerData::CreateBusPatchHisto(Int_t busPatchId, Int_t dim)
{
  /// Create histogram of a given bus patch. Note that in order
  /// to keep memory footprint as low as possible, you should delete
  /// the returned pointer as soon as possible...
  
  TH1* h(0x0);
  
  if ( HasBusPatch(busPatchId) && IsHistogrammed(dim)) 
  {
    h = CreateHisto(Form("BP%04d_%d",busPatchId,dim),dim);
    AddBusPatchHisto(*h,busPatchId,dim);
  }
  
  return h;
}  
 
//_____________________________________________________________________________
void
AliMUONTrackerData::AddBusPatchHisto(TH1& h, Int_t busPatchId, Int_t dim)
{
  /// Add data from one bus patch to the histogram
      
  if ( HasBusPatch(busPatchId ) )
  {
    AliMpBusPatch* busPatch = AliMpDDLStore::Instance()->GetBusPatch(busPatchId);
    for ( Int_t i = 0; i < busPatch->GetNofManus(); ++i ) 
    {
      Int_t manuId = busPatch->GetManuId(i);
      AddManuHisto(h,busPatch->GetDEId(),manuId,dim);
    }
  }
}


//_____________________________________________________________________________
TH1*
AliMUONTrackerData::CreateDEHisto(Int_t detElemId, Int_t dim)
{
  /// Create histogram of a given detection element. Note that in order
  /// to keep memory footprint as low as possible, you should delete
  /// the returned pointer as soon as possible...
  
  TH1* h(0x0);
  
  if ( HasDetectionElement(detElemId) && IsHistogrammed(dim) ) 
  {
    h = CreateHisto(Form("DE%04d-%d",detElemId,dim),dim);
    AddDEHisto(*h,detElemId,dim);
  }
  
  return h;
}

//_____________________________________________________________________________
void
AliMUONTrackerData::AddDEHisto(TH1& h, Int_t detElemId, Int_t dim)
{
  /// Add data from one detection element to the histogram
  
  if ( HasDetectionElement(detElemId) )
  {
    AliMpDetElement* de = AliMpDDLStore::Instance()->GetDetElement(detElemId);
    for ( Int_t i = 0; i < de->GetNofBusPatches(); ++ i ) 
    {
      Int_t busPatchId = de->GetBusPatchId(i);
      AddBusPatchHisto(h,busPatchId,dim);
    }
  }
}

//_____________________________________________________________________________
TH1*
AliMUONTrackerData::CreateManuHisto(Int_t detElemId, Int_t manuId, Int_t dim)
{
  /// Create histogram of a given manu. Note that in order
  /// to keep memory footprint as low as possible, you should delete
  /// the returned pointer as soon as possible...
  
  TH1* h(0x0);
  
  if ( HasManu(detElemId, manuId) && IsHistogrammed(dim) ) 
  {
    h = CreateHisto(Form("DE%04dMANU%04d_%d",detElemId,manuId,dim),dim);
    AddManuHisto(*h,detElemId,manuId,dim);
  }
  
  return h;
}

//_____________________________________________________________________________
void
AliMUONTrackerData::AddManuHisto(TH1& h, Int_t detElemId, Int_t manuId, Int_t dim)
{
  /// Add data from a given manu to histogram
  
  if ( HasManu(detElemId,manuId) )
  {
    for ( Int_t i = 0; i < AliMpConstants::ManuNofChannels(); ++i ) 
    {
      if ( HasChannel(detElemId,manuId,i) )
      {
        AliMUONSparseHisto* sh = GetChannelHisto(detElemId,manuId,i,dim);
      
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
AliMUONTrackerData::CreatePCBHisto(Int_t /*detElemId*/, Int_t /*pcbIndex*/, Int_t /*dim*/)
{
  /// Create histogram of a given PCB. Note that in order
  /// to keep memory footprint as low as possible, you should delete
  /// the returned pointer as soon as possible...
  
 // TH1* h(0x0);
//  
//  if ( HasPCB(detElemId, pcbIndex) && IsHistogrammed(dim)) 
//  {
//    h = CreateHisto(Form("DE%04dPCB1d_%d",detElemId,pcbIndex,dim),dim);
//  }
//  
//  return h;

  AliWarning("Not implemented (is it needed ?)");
  return 0x0;
}

//_____________________________________________________________________________
TH1*
AliMUONTrackerData::CreateChamberHisto(Int_t chamberId, Int_t dim)
{
  /// Create histogram of a given chamber. Note that in order
  /// to keep memory footprint as low as possible, you should delete
  /// the returned pointer as soon as possible...
  
  TH1* h(0x0);
  
  if ( HasChamber(chamberId) && IsHistogrammed(dim))
  {
    h = CreateHisto(Form("CHAMBER%02d_%d",chamberId,dim),dim);
    AliMpDEIterator it;
    it.First(chamberId);
    while ( !it.IsDone() )
    {
      Int_t detElemId = it.CurrentDEId();
      AddDEHisto(*h,detElemId,dim);
      it.Next();
    }
  }
  
  return h;
}

//_____________________________________________________________________________
AliMUONSparseHisto*
AliMUONTrackerData::GetChannelHisto(Int_t detElemId, Int_t manuId, 
                                    Int_t manuChannel, Int_t dim)
{
  /// Get histogram of a given channel
  
  if (!fChannelHistos) fChannelHistos = new AliMUON2DMap(kTRUE);
  
  TObjArray* dimArray = static_cast<TObjArray*>(fChannelHistos->FindObject(detElemId,manuId));
  if (!dimArray)
  {
    dimArray = new TObjArray(fExternalDimension);
    dimArray->SetUniqueID( ( manuId << 16 ) | detElemId );
    fChannelHistos->Add(dimArray);
  }
  
  TObjArray* channels = static_cast<TObjArray*>(dimArray->UncheckedAt(dim));
  if (!channels)
  {
    channels = new TObjArray(AliMpConstants::ManuNofChannels());
    dimArray->AddAt(channels,dim);
  }
  
  AliMUONSparseHisto* h = static_cast<AliMUONSparseHisto*>(channels->UncheckedAt(manuChannel));
  if (!h)
  {
    h = new AliMUONSparseHisto;
    h->SetUniqueID(( manuChannel << 16 ) | dim);
    channels->AddAt(h,manuChannel);
  }
  
  return h;
 
// below is an alternate implementation, using a 1DMap, which *seems* to be
// slightly SLOWER.
//
//  AliMUON1DMap* m = static_cast<AliMUON1DMap*>(fChannelHistos->FindObject(detElemId,manuId));
//  if (!m)
//  {
//    m = new AliMUON1DMap(kFALSE);
//    m->SetUniqueID( ( manuId << 16 ) | detElemId );
//    fChannelHistos->Add(m);
//  }
//  
//  UInt_t uid = ( manuChannel << 16 ) | dim;
//  
//  AliMUONSparseHisto* h = static_cast<AliMUONSparseHisto*>(m->FindObject(uid));
//  if (!h)
//  {
//    h = new AliMUONSparseHisto;
//    
//    h->SetUniqueID(uid);
//    
//    m->Add(h);
//  }

  return h;
}

//_____________________________________________________________________________
Int_t
AliMUONTrackerData::GetParts(AliMUONVCalibParam* external,
                             AliMUONVCalibParam*& chamber,
                             AliMUONVCalibParam*& de,
                             AliMUONVCalibParam*& busPatch,
                             AliMUONVCalibParam*& pcb,
                             AliMUONVCalibParam*& manu,
                             AliMUONVCalibParam*& channel,
                             AliMpDetElement*& mpde)
{
  /// Get containers at all levels
 
  AliMpDDLStore* ddlStore = AliMpDDLStore::Instance();
  
  Int_t detElemId = external->ID0();
  
  Int_t chamberId = AliMpDEManager::GetChamberId(detElemId);
  
  Int_t manuId = external->ID1();
  
  mpde = ddlStore->GetDetElement(detElemId);
  
  Int_t busPatchId = ddlStore->GetBusPatchId(detElemId,manuId);
  
  Int_t pcbIndex = -1;
  
  AliMp::StationType stationType = mpde->GetStationType();
  
  if ( stationType == AliMp::kStation345 ) 
  {
    AliMpHVNamer namer;
    pcbIndex = namer.ManuId2PCBIndex(detElemId,manuId);
  }
  
  channel = ChannelParam(detElemId,manuId,external);
  
  manu = ManuParam(detElemId,manuId,kTRUE);
  
  busPatch = BusPatchParam(busPatchId,kTRUE);
  
  de = DetectionElementParam(detElemId,kTRUE);
  
  chamber = ChamberParam(chamberId,kTRUE);
  
  pcb = 0x0;
  
  if ( pcbIndex >= 0 ) 
  {
    pcb = PCBParam(detElemId,pcbIndex,kTRUE);
  }
  
  return manuId;
}

//_____________________________________________________________________________
Bool_t 
AliMUONTrackerData::HasBusPatch(Int_t busPatchId) const
{
  /// Whether we have data for a given buspatch
  return ( BusPatchParam(busPatchId) != 0 );
}

//_____________________________________________________________________________
Bool_t 
AliMUONTrackerData::HasChamber(Int_t chamberId) const
{
  /// Whether we have data for a given chamber
  return ( ChamberParam(chamberId) != 0 );
}

//_____________________________________________________________________________
Bool_t 
AliMUONTrackerData::HasDetectionElement(Int_t detElemId) const
{
  /// Whether we have data for a given detection element
  return ( DetectionElementParam(detElemId) != 0 );
}

//_____________________________________________________________________________
Bool_t
AliMUONTrackerData::HasManu(Int_t detElemId, Int_t manuId) const
{
  /// Whether we have data for a given manu
  return ( ManuParam(detElemId,manuId) != 0 ); 
}

//_____________________________________________________________________________
Bool_t
AliMUONTrackerData::HasPCB(Int_t detElemId, Int_t pcbIndex) const
{
  /// Whether we have data for a given pcb
  return ( PCBParam(detElemId,pcbIndex) != 0 ); 
}

//_____________________________________________________________________________
Double_t 
AliMUONTrackerData::Manu(Int_t detElemId, Int_t manuId, Int_t dim) const
{
  /// Return the value for a given manu and a given dimension
  
  AliMUONVCalibParam* param = ManuParam(detElemId,manuId);
  return param ? Value(*param,0,dim) : 0.0;
}

//_____________________________________________________________________________
AliMUONVCalibParam* 
AliMUONTrackerData::ManuParam(Int_t detElemId, Int_t manuId, Bool_t create) const
{
  /// Get the VCalibParam for a given manu
  
  AliMUONVCalibParam* manu = fManuValues ? static_cast<AliMUONVCalibParam*>
    (fManuValues->FindObject(detElemId,manuId)) : 0x0 ;
  
  if (!manu && create && fManuValues)
  {
    manu = CreateManuParam(detElemId,manuId);
    fManuValues->Add(manu);
  }
  
  return manu;
}

//_____________________________________________________________________________
AliMUONVCalibParam* 
AliMUONTrackerData::CreateManuParam(Int_t detElemId, Int_t manuId) const
{
  /// Create storage for one manu
  
  AliCodeTimerAuto("");
  
  AliMUONVCalibParam* manu = new AliMUONCalibParamND(Dimension(),1,detElemId,manuId,0.0);
  
  // set the number of channels in that manu
  
  AliMpDetElement* de = AliMpDDLStore::Instance()->GetDetElement(detElemId);
  
  manu->SetValueAsDouble(0,IndexOfNumberDimension(),de->NofChannelsInManu(manuId));
  
  return manu;
}

//_____________________________________________________________________________
Int_t 
AliMUONTrackerData::NumberOfDimensions() const
{
  /// Number of dimensions we're dealing with
  
  return fDimension + fgkVirtualExtraDimension; 
}

//_____________________________________________________________________________
Double_t 
AliMUONTrackerData::PCB(Int_t detElemId, Int_t pcbIndex, Int_t dim) const
{
  /// Return the value of a given pcb for a given dimension

  AliMUONVCalibParam* param = PCBParam(detElemId,pcbIndex);
  
  return param ? Value(*param,pcbIndex,dim) : 0.0;
}

//_____________________________________________________________________________
AliMUONVCalibParam* 
AliMUONTrackerData::PCBParam(Int_t detElemId, Int_t pcbIndex, Bool_t create) const
{
  /// Return (if it exist) the VCalibParam for a given pcb

  AliMUONVCalibParam* pcb =  fPCBValues ? static_cast<AliMUONVCalibParam*>
    (fPCBValues->FindObject(detElemId,pcbIndex)) : 0x0 ;
  
  if (create && fPCBValues && !pcb)
  {
    pcb = CreatePCBParam(detElemId,pcbIndex);
    fPCBValues->Add(pcb);
  }
  
  return pcb;
}

//_____________________________________________________________________________
AliMUONVCalibParam* 
AliMUONTrackerData::CreatePCBParam(Int_t detElemId, Int_t pcbIndex) const
{
  /// Create storage for one PCB (station345 only)
  
  AliCodeTimerAuto("");
  
  AliMpHVNamer namer;
  
  AliMUONVCalibParam* pcb = new AliMUONCalibParamND(Dimension(),
                                                    namer.NumberOfPCBs(detElemId),
                                                    detElemId,
                                                    pcbIndex,
                                                    0.0);
  return pcb;
}

//_____________________________________________________________________________
void 
AliMUONTrackerData::Print(Option_t* wildcard, Option_t* opt) const
{
  /// Printout
  
  TNamed::Print(opt);
  
  if ( fIsRunnable ) 
  {
    cout << " Nevents=" << fNevents << endl;
  }

  for ( Int_t i = 0; i <= fExternalDimensionNames->GetLast(); ++i ) 
  {
    TObjString* name = static_cast<TObjString*>(fExternalDimensionNames->At(i));
    cout << Form("External Dimension %2d Name %s %s",i,
                 ( name ? name->String().Data() : "null"),
                 ( IsHistogrammed(i) ? "(histogrammed)" : "")) << endl;
  }
  
  for ( Int_t i = 0; i <= fDimensionNames->GetLast(); ++i ) 
  {
    TObjString* name = static_cast<TObjString*>(fDimensionNames->At(i));
    cout << Form("Internal Dimension %2d Name %s",i,
                 ( name ? name->String().Data() : "null")) << endl;
  }
    
  TString sopt(opt);
  sopt.ToUpper();
  
  if ( sopt.Contains("CHANNEL") && fChannelValues ) 
  {
    fChannelValues->Print(wildcard,opt);
  }

  if ( sopt.Contains("MANU") && fManuValues ) 
  {
    fManuValues->Print(wildcard,opt);
  }

  if ( sopt.Contains("BUSPATCH") && fBusPatchValues ) 
  {
    fBusPatchValues->Print(wildcard,opt);
  }

  if ( sopt.Contains("DE") && fDEValues ) 
  {
    fDEValues->Print(wildcard,opt);
  }

  if ( sopt.Contains("CHAMBER") && fChamberValues ) 
  {
    fChamberValues->Print(wildcard,opt);
  }
  
}

//_____________________________________________________________________________
void
AliMUONTrackerData::SetDimensionName(Int_t index, const char* name)
{  
  /// Set the name of a given dimension

  if ( index >= fExternalDimension ) 
  {
    AliError(Form("Index out of bounds : %d / %d",index,fExternalDimension));
    return;
  }
  
  Int_t ix = External2Internal(index);
  
  const char* prefix[] = { "mean", "sigma" };
  
  for ( Int_t i = 0; i < 2; ++i ) 
  {
    Int_t j = ix+i;
    
    SetInternalDimensionName(j,Form("%s of %s",prefix[i],name));
  }
  
  SetExternalDimensionName(index,name);
}

//_____________________________________________________________________________
void 
AliMUONTrackerData::SetHistogramDimension(Int_t index, Bool_t value)
{
  /// decide to make histos for a given dimension
  if ( index >= ExternalDimension() ) 
  {
    AliError(Form("Index out of bounds : %d / %d",index,ExternalDimension()));
    return;
  }
  
  fHistogramming[index] = value;
}

//_____________________________________________________________________________
void 
AliMUONTrackerData::SetInternalDimensionName(Int_t index, const char* value)
{
  /// Set the name of a given internal dimension
  if ( index >= fDimension ) 
  {
    AliError(Form("Index out of bounds : %d / %d",index,fDimension));
    return;
  }
  
  TObjString* ovalue = static_cast<TObjString*>(fDimensionNames->At(index));
    
  if ( ovalue ) 
  {
    fDimensionNames->Remove(ovalue);
    delete ovalue;
  }
  fDimensionNames->AddAt(new TObjString(value),index);
}

//_____________________________________________________________________________
void 
AliMUONTrackerData::SetExternalDimensionName(Int_t index, const char* value)
{
  /// Set the name of a given external dimension
  if ( index >= fExternalDimension ) 
  {
    AliError(Form("Index out of bounds : %d / %d",index,fExternalDimension));
    return;
  }
  
  TObjString* ovalue = static_cast<TObjString*>(fExternalDimensionNames->At(index));
  
  if ( ovalue ) 
  {
    fExternalDimensionNames->Remove(ovalue);
    delete ovalue;
  }
  fExternalDimensionNames->AddAt(new TObjString(value),index);
}

//_____________________________________________________________________________
Double_t 
AliMUONTrackerData::Value(const AliMUONVCalibParam& param, Int_t i, Int_t dim) const
{
  /// Compute the value for a given dim, using the internal information we have
  /// Basically we're converting sum of weights and sum of squares of weights
  /// into means and sigmas, and number of events into occupancy number.

  Double_t n = param.ValueAsDouble(i,IndexOfNumberDimension());
  
  if ( dim == IndexOfNumberDimension() ) return n; // the number of channels in any given element does not depend on the number of events
  
  Double_t occ = param.ValueAsDouble(i,IndexOfOccupancyDimension());

  if ( dim >= fDimension ) 
  {
    return occ;
  }
  
  if ( dim == IndexOfOccupancyDimension() ) return occ/n/NumberOfEvents();
  
  Double_t value = param.ValueAsDouble(i,dim);
  
  if ( value >= AliMUONVCalibParam::InvalidFloatValue() ) return AliMUONVCalibParam::InvalidFloatValue();
  
  if ( TMath::Even(dim) ) 
  {
    return value/occ;
  }
  else
  {
    Double_t n = occ;
    
    Double_t mean = param.ValueAsDouble(i,dim-1)/n;
    
    return  TMath::Sqrt(TMath::Abs((value-n*mean*mean)/(n-1.0)));
  }
  
  AliError("Why am I here ?");
  return 0.0;
}

