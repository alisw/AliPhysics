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

#include "AliMUONCalibParamND.h"
#include "AliMUONVStore.h"
#include "AliMpBusPatch.h"
#include "AliMpDDLStore.h"
#include "AliMpDEManager.h"
#include "AliMpDEIterator.h"
#include "AliMpDetElement.h"
#include "AliMpHVNamer.h"
#include "AliCodeTimer.h"
#include "AliLog.h"
#include <Riostream.h>
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
fExternalDimension(dimension),
fIsRunnable(runnable)
{  
  /// ctor
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
}

//_____________________________________________________________________________
Bool_t
AliMUONTrackerData::Add(const AliMUONVStore& store)
{
  /// We first convert the external store to a temporary internal store
  /// with more dimension (2*store's dimension)
  
  AliCodeTimerAuto(GetName())
  
  Int_t ndim(NumberOfDimensions()-fgkExtraDimension-fgkVirtualExtraDimension); 
  
  AliMUONVStore* istore = store.Create();
  
  TIter next(store.CreateIterator());
  AliMUONVCalibParam* external;
  
  AliCodeTimerStart("from external to internal store");
  
  while ( ( external = static_cast<AliMUONVCalibParam*>(next()) ) )
  {
    Int_t detElemId = external->ID0();
    Int_t manuId = external->ID1();
    
    AliMpDetElement* de = AliMpDDLStore::Instance()->GetDetElement(detElemId);
    
    AliMUONVCalibParam* internal = static_cast<AliMUONVCalibParam*>
      (istore->FindObject(detElemId,manuId));
    
    if (!internal)
    {
      internal = new AliMUONCalibParamND(ndim,external->Size(),
                                         detElemId, manuId, 
                                         0.0);
      istore->Add(internal);
    }
    
    for ( Int_t i = 0; i < external->Size(); ++i ) 
    {
      Bool_t connectPad = de->IsConnectedChannel(manuId,i);
      
      if (!connectPad) continue;
      
      for ( Int_t j = 0; j < external->Dimension(); ++j )
      {
        Int_t ix = External2Internal(j);
        
        Double_t vext = external->IsDoublePrecision() ? 
          external->ValueAsDouble(i,j) :
          external->ValueAsFloat(i,j);
        
        Double_t sumw = internal->ValueAsDouble(i,ix) + vext;
        Double_t sumw2 = internal->ValueAsDouble(i,ix+1) + vext*vext;
        
        internal->SetValueAsFloat(i,ix,sumw);
        internal->SetValueAsFloat(i,ix+1,sumw2);
      }
    }
  }
  
  AliCodeTimerStop("from external to internal store");
  
  /// and we add this internal store to what we already have
  
  InternalAdd(*istore);
  
  /// delete the temporary internal store.
  AliCodeTimerStart("delete");
  delete istore;
  AliCodeTimerStop("delete");
  
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
AliMUONTrackerData::BusPatchParam(Int_t busPatchId) const
{
  /// Return (if it exist), the VCalibParam for a given busPatch
  return fBusPatchValues ? static_cast<AliMUONVCalibParam*>
  (fBusPatchValues->FindObject(busPatchId)) : 0x0;
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
AliMUONTrackerData::ChamberParam(Int_t chamberId) const
{
  /// Return (if it exist) the VCalibParam for a given chamber
  return fChamberValues ? static_cast<AliMUONVCalibParam*>
  (fChamberValues->FindObject(chamberId)) : 0x0;
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
AliMUONTrackerData::ChannelParam(Int_t detElemId, Int_t manuId) const
{
  /// Return (if it exist) the VCalibParam for a given manu
  return fChannelValues ? static_cast<AliMUONVCalibParam*>
  (fChannelValues->FindObject(detElemId,manuId)) : 0x0 ;
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
  if ( fChamberValues) fChamberValues->Clear();
  fNevents = 0;
  NumberOfEventsChanged();
}

//_____________________________________________________________________________
Double_t 
AliMUONTrackerData::Count(Int_t detElemId, Int_t manuId, 
                          Int_t manuChannel) const
{
  /// Return the number of times a given channel had data
  AliMUONVCalibParam* param = static_cast<AliMUONVCalibParam*>
  (fChannelValues->FindObject(detElemId,manuId));
  
  if ( !param ) return 0.0;
  
  return param->ValueAsDouble(manuChannel,IndexOfOccupancyDimension());
}

//_____________________________________________________________________________
AliMUONVCalibParam*
AliMUONTrackerData::CreateDouble(const AliMUONVCalibParam& param) const
{
  /// Create a double version of VCalibParam, for internal use
  AliMUONVCalibParam* c = new AliMUONCalibParamND(param.Dimension()+fgkExtraDimension,
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
AliMUONTrackerData::DetectionElementParam(Int_t detElemId) const
{
  /// Return (if it exist) the VCalibParam for a given detection element
  return fDEValues ? static_cast<AliMUONVCalibParam*>
  (fDEValues->FindObject(detElemId)) : 0x0 ;
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
Bool_t 
AliMUONTrackerData::InternalAdd(const AliMUONVStore& store)
{
  /// Add the given store to our internal store
  /// Store must be of dimension = fDimension-1
  
  AliCodeTimerAuto(GetName());
  
  AliMpDDLStore* ddlStore = AliMpDDLStore::Instance();
  
  ++fNevents;
  NumberOfEventsChanged();
  
  if (!fChannelValues)
  {
    fChannelValues = store.Create();
    fManuValues = store.Create();
    fBusPatchValues = store.Create();
    fDEValues = store.Create();
    fChamberValues = store.Create();
    fPCBValues = store.Create();
  }
  
  TIter next(store.CreateIterator());
  AliMUONVCalibParam* external;
  
  AliMpHVNamer namer;
  
  while ( ( external = static_cast<AliMUONVCalibParam*>(next()) ) )
  {
    if ( external->Dimension() != fDimension-fgkExtraDimension ) 
    {
      AliError(Form("Incompatible dimensions %d vs %d",
                    external->Dimension(),fDimension-fgkExtraDimension));
      return kFALSE;
    }
    
    Int_t detElemId = external->ID0();
    
    AliMp::StationType stationType = AliMpDEManager::GetStationType(detElemId);
    
    Int_t chamberId = AliMpDEManager::GetChamberId(detElemId);
    
    Int_t manuId = external->ID1();
    
    AliMpDetElement* mpde = ddlStore->GetDetElement(detElemId);

    Int_t busPatchId = ddlStore->GetBusPatchId(detElemId,manuId);
    
    Int_t pcbIndex = -1;
    
    if ( stationType == AliMp::kStation345 ) 
    {
      pcbIndex = namer.ManuId2PCBIndex(detElemId,manuId);
    }

    AliMUONVCalibParam* channel = ChannelParam(detElemId,manuId);
    if (!channel)
    {
      channel = CreateDouble(*external);
      fChannelValues->Add(channel);
    }

    AliMUONVCalibParam* manu = ManuParam(detElemId,manuId);
    if (!manu)
    {
      manu = new AliMUONCalibParamND(external->Dimension()+fgkExtraDimension,
                                     1,
                                     detElemId,
                                     manuId,
                                     0.0);
      
      // set the number of channels in that manu
      
      AliMpDetElement* de = AliMpDDLStore::Instance()->GetDetElement(detElemId);
      
      manu->SetValueAsDouble(0,IndexOfNumberDimension(),de->NofChannelsInManu(manuId));
      
      fManuValues->Add(manu);
    }
    
    AliMUONVCalibParam* busPatch = BusPatchParam(busPatchId);
    if (!busPatch)
    {
      AliMpBusPatch* bp = AliMpDDLStore::Instance()->GetBusPatch(busPatchId);

      if (!bp)
      {
        AliError(Form("Got an invalid buspatchId = %d",busPatchId));
        continue;
      }
      
      busPatch = new AliMUONCalibParamND(external->Dimension()+fgkExtraDimension,
                                         1,
                                         busPatchId,
                                         0,
                                         0.0);
      
      // set the number of channels in that buspatch

      Int_t nchannels(0);
      
      AliMpDetElement* de = AliMpDDLStore::Instance()->GetDetElement(detElemId);

      for ( Int_t i = 0; i < bp->GetNofManus(); ++i ) 
      {
        Int_t manuId = bp->GetManuId(i);
        nchannels += de->NofChannelsInManu(manuId);
      }

      busPatch->SetValueAsDouble(0,IndexOfNumberDimension(),nchannels);
      
      fBusPatchValues->Add(busPatch);
    }

    AliMUONVCalibParam* de = DetectionElementParam(detElemId);
    if (!de)
    {
      de = new AliMUONCalibParamND(external->Dimension()+fgkExtraDimension,
                                         1,
                                         detElemId,
                                         0,
                                         0.0);
      
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
      
      fDEValues->Add(de);
    }

    AliMUONVCalibParam* chamber = ChamberParam(chamberId);
    if (!chamber)
    {
      chamber = new AliMUONCalibParamND(external->Dimension()+fgkExtraDimension,
                                   1,
                                   chamberId,
                                   0,
                                   0.0);

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
      
      fChamberValues->Add(chamber);
    }

    AliMUONVCalibParam* pcb = 0x0;
    
    if ( pcbIndex >= 0 ) 
    {
      pcb = PCBParam(detElemId,pcbIndex);
      if (!pcb)
      {
        pcb = new AliMUONCalibParamND(external->Dimension()+fgkExtraDimension,
                                        namer.NumberOfPCBs(detElemId),
                                        detElemId,
                                        pcbIndex,
                                        0.0);
        fPCBValues->Add(pcb);
      }
    }
    
    for ( Int_t i = 0; i < external->Size(); ++i ) 
    {
      Bool_t existingChannel = mpde->IsConnectedChannel(manuId,i);
      
      if ( existingChannel ) 
      {
        Bool_t validChannel(kFALSE);
        
        for ( Int_t j = 0; j < external->Dimension(); ++j )
        {
          if ( external->ValueAsFloat(i,j) >= AliMUONVCalibParam::InvalidFloatValue() ) continue;
          
          validChannel = kTRUE;
          
          Double_t vext = external->IsDoublePrecision() ? 
            external->ValueAsDouble(i,j) :
            external->ValueAsFloat(i,j);
          
          Double_t value = channel->ValueAsDouble(i,j) + vext;
          
          channel->SetValueAsDouble(i,j,value);
          
          manu->SetValueAsDouble(0,j,manu->ValueAsDouble(0,j)+vext);  
          
          busPatch->SetValueAsDouble(0,j,busPatch->ValueAsDouble(0,j)+vext);

          de->SetValueAsDouble(0,j,de->ValueAsDouble(0,j)+vext);

          chamber->SetValueAsDouble(0,j,chamber->ValueAsDouble(0,j)+vext);

          if ( pcb ) 
          {
            pcb->SetValueAsDouble(pcbIndex,j,pcb->ValueAsDouble(pcbIndex,j)+vext);
          }
        }
        
        if ( validChannel )
        {
          channel->SetValueAsDouble(i,IndexOfOccupancyDimension(),
                                  channel->ValueAsDouble(i,IndexOfOccupancyDimension())+1.0);
          manu->SetValueAsDouble(0,IndexOfOccupancyDimension(),
                               manu->ValueAsDouble(0,IndexOfOccupancyDimension())+1.0);        
          busPatch->SetValueAsDouble(0,IndexOfOccupancyDimension(),
                                 busPatch->ValueAsDouble(0,IndexOfOccupancyDimension())+1.0);        
          de->SetValueAsDouble(0,IndexOfOccupancyDimension(),
                                     de->ValueAsDouble(0,IndexOfOccupancyDimension())+1.0);        
          chamber->SetValueAsDouble(0,IndexOfOccupancyDimension(),
                               chamber->ValueAsDouble(0,IndexOfOccupancyDimension())+1.0); 
          if ( pcb ) 
          {
            pcb->SetValueAsDouble(pcbIndex,IndexOfOccupancyDimension(),
                                    pcb->ValueAsDouble(pcbIndex,IndexOfOccupancyDimension())+1.0);        
          }
        }
      }
    }
  }

  return kTRUE;
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
AliMUONTrackerData::ManuParam(Int_t detElemId, Int_t manuId) const
{
  /// Get the VCalibParam for a given manu
  return fManuValues ? static_cast<AliMUONVCalibParam*>
  (fManuValues->FindObject(detElemId,manuId)) : 0x0 ;
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
AliMUONTrackerData::PCBParam(Int_t detElemId, Int_t pcbIndex) const
{
  /// Return (if it exist) the VCalibParam for a given pcb
  return fPCBValues ? static_cast<AliMUONVCalibParam*>
  (fPCBValues->FindObject(detElemId,pcbIndex)) : 0x0 ;
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
  
  for ( Int_t i = 0; i <= fDimensionNames->GetLast(); ++i ) 
  {
    TObjString* name = static_cast<TObjString*>(fDimensionNames->At(i));
    cout << Form("Dimension %2d Name %s",i,
                 ( name ? name->String().Data() : "null")) << endl;
  }
  
  cout << Form("External Dimensions = %d",fExternalDimension) << endl;  

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
    Double_t sumw = param.ValueAsDouble(i,dim-1);
    Double_t mean = sumw/n;
    
    return  TMath::Sqrt(TMath::Abs(value/occ - mean*mean));
  }
  
  AliError("Why am I here ?");
  return 0.0;
}

