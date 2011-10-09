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
#include "AliDAQ.h"
#include "AliLog.h"
#include "AliMUON1DArray.h"
#include "AliMUON1DMap.h"
#include "AliMUON2DMap.h"
#include "AliMUONCalibParamND.h"
#include "AliMUONRejectList.h"
#include "AliMUONSparseHisto.h"
#include "AliMUONVStore.h"
#include "AliMpBusPatch.h"
#include "AliMpConstants.h"
#include "AliMpCDB.h"
#include "AliMpDDLStore.h"
#include "AliMpManuStore.h"
#include "AliMpDEIterator.h"
#include "AliMpDEManager.h"
#include "AliMpDetElement.h"
#include "AliMpDCSNamer.h"
#include "AliMpManuIterator.h"
#include "AliMpEncodePair.h"
#include "AliMpSegmentation.h"
#include <Riostream.h>
#include <TClass.h>
#include <TMath.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <TString.h>
#include <TTimeStamp.h>
#include <TVector2.h>
#include <cassert>
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
                                       Bool_t issingleevent)
: AliMUONVTrackerData(name,title),
fIsSingleEvent(issingleevent),
fChannelValues(0x0),
fManuValues(0x0),
fBusPatchValues(0x0),
fDEValues(0x0),
fChamberValues(0x0),
fPCBValues(0x0),
fDimension(External2Internal(dimension)+fgkExtraDimension),
fNevents(0),
fDimensionNames(new TObjArray(fDimension+fgkVirtualExtraDimension)),
fExternalDimensionNames(new TObjArray(dimension)),
fExternalDimension(dimension),
fHistogramming(new Int_t[fExternalDimension]),
fHistos(0x0),
fXmin(0.0),
fXmax(0.0),
fIsChannelLevelEnabled(kTRUE),
fIsManuLevelEnabled(kTRUE),
fIsBustPatchLevelEnabled(kTRUE),
fIsPCBLevelEnabled(kTRUE),
fNofDDLs(0),
fNofEventsPerDDL(0x0)
{  
  /// ctor
  memset(fHistogramming,0,fExternalDimension*sizeof(Int_t)); // histogramming is off by default. Use MakeHistogramForDimension to turn it on.
  fExternalDimensionNames->SetOwner(kTRUE);
  fDimensionNames->SetOwner(kTRUE);  
  fDimensionNames->AddAt(new TObjString("occ"),IndexOfOccupancyDimension());
  fDimensionNames->AddAt(new TObjString("N"),IndexOfNumberDimension());
  fDimensionNames->AddAt(new TObjString("n"),NumberOfDimensions()-fgkVirtualExtraDimension);
  Clear();
}

//_____________________________________________________________________________
AliMUONTrackerData::AliMUONTrackerData(const char* name, const char* title,
                                       const AliMUONVStore& manuValues)
: AliMUONVTrackerData(name,title),
fIsSingleEvent(kFALSE),
fChannelValues(0x0),
fManuValues(0x0),
fBusPatchValues(0x0),
fDEValues(0x0),
fChamberValues(0x0),
fPCBValues(0x0),
fDimension(0),
fNevents(0),
fDimensionNames(0x0),
fExternalDimensionNames(0x0),
fExternalDimension(0),
fHistogramming(0x0),
fHistos(0x0),
fXmin(0.0),
fXmax(0.0),
fIsChannelLevelEnabled(kFALSE),
fIsManuLevelEnabled(kTRUE),
fIsBustPatchLevelEnabled(kTRUE),
fIsPCBLevelEnabled(kTRUE),
fNofDDLs(0),
fNofEventsPerDDL(0x0)
{  
  /// ctor with pre-computed values at the manu level
  /// In this case, we force fIsChannelLevelEnabled = kFALSE
  /// ctor
  
  if (manuValues.GetSize()==0)
  {
    AliFatal("Cannot create a tracker data from nothing in that case !");
  }
  
  if ( !AliMpDDLStore::Instance(kFALSE) && !AliMpManuStore::Instance(kFALSE) )
  {
    AliError("Cannot work without (full) mapping");
    return;
  }
  
  TIter next(manuValues.CreateIterator());
  AliMUONVCalibParam* m = static_cast<AliMUONVCalibParam*>(next());
  
  Int_t dimension = ( m->Dimension() - fgkExtraDimension - fgkVirtualExtraDimension ) / 2;
  
  fDimension = External2Internal(dimension)+fgkExtraDimension;
  
  fDimensionNames = new TObjArray(fDimension+fgkVirtualExtraDimension);
  fExternalDimensionNames = new TObjArray(dimension);
  fExternalDimension = dimension;
  fHistogramming = new Int_t[fExternalDimension];
  memset(fHistogramming,0,fExternalDimension*sizeof(Int_t)); // histogramming is off by default. Use MakeHistogramForDimension to turn it on.

  fExternalDimensionNames->SetOwner(kTRUE);
  fDimensionNames->SetOwner(kTRUE);  
  fDimensionNames->AddAt(new TObjString("occ"),IndexOfOccupancyDimension());
  fDimensionNames->AddAt(new TObjString("N"),IndexOfNumberDimension());
  fDimensionNames->AddAt(new TObjString("n"),NumberOfDimensions()-fgkVirtualExtraDimension);
  Clear();
  TArrayI nevents(AliDAQ::NumberOfDdls("MUONTRK"));
  AssertStores();
  
  next.Reset();
  AliMUONVCalibParam* external;
  
  while ( ( external = static_cast<AliMUONVCalibParam*>(next()) ) )
  {
    Int_t detElemId, manuId;
    
    GetDEManu(*external,detElemId,manuId);
    
    AliMUONVCalibParam* chamber(0x0);
    AliMUONVCalibParam* de(0x0);
    AliMUONVCalibParam* busPatch(0x0);
    AliMUONVCalibParam* pcb(0x0);
    AliMUONVCalibParam* manu(0x0);
    AliMUONVCalibParam* channel(0x0);
    AliMpDetElement* mpde(0x0);
    
    AliMUONVCalibParam* wec = new AliMUONCalibParamND(external->Dimension()-1,1,detElemId,manuId,0);
    // as external, but without event count
    wec->SetValueAsDouble(0,0,external->ValueAsDouble(0,0));
    wec->SetValueAsDouble(0,1,external->ValueAsDouble(0,1));
    wec->SetValueAsDouble(0,2,external->ValueAsDouble(0,2));
    wec->SetValueAsDouble(0,3,external->ValueAsDouble(0,3));
    
    Int_t mid = GetParts(wec,chamber,de,busPatch,pcb,manu,channel,mpde);

    if ( manuId != mid ) 
    {
      AliError(Form("Something is wrong for DE %5d : manuId = %d vs mid = %d",detElemId,manuId,mid));
      continue;
    }
    
    if ( manuId < 0 ) 
    {
      AliError("Got a < 0 manuId. Should not happen here !");
      continue;
    }
    
    assert(channel==0x0);
    
    Int_t n1 = manu->ValueAsInt(0,IndexOfNumberDimension());
    Int_t n2 = external->ValueAsInt(0,IndexOfNumberDimension());
    if  ( n1 != n2 )
    {
      AliError(Form("Incoherent number of manu channels for DE %5d MANU %5d : %d vs %d",
                    detElemId,manuId,n1,n2));
    }
    
    Int_t nevt = external->ValueAsInt(0,4);
    
    Int_t busPatchId = AliMpDDLStore::Instance()->GetBusPatchId(detElemId,manuId);

    Int_t ddl = AliMpDDLStore::Instance()->GetDDLfromBus(busPatchId);

    if ( nevents[ddl] == 0 ) 
    {
      nevents[ddl] = nevt;
    }
    else
    {
      if ( nevents.At(ddl) != nevt ) 
      {
        AliError(Form("Nevt mismatch for DE %5d MANU %5d DDL %d : %d vs %d",
                      detElemId,manuId,ddl,nevents.At(ddl),nevt));
        continue;
      }
    }
    
    for ( Int_t i = 0; i < wec->Dimension()-1; ++i ) 
    {
      manu->SetValueAsDouble(0,i,manu->ValueAsDouble(0,i) + wec->ValueAsDouble(0,i));
      
      busPatch->SetValueAsDouble(0,i,busPatch->ValueAsDouble(0,i) + wec->ValueAsDouble(0,i));
      
      de->SetValueAsDouble(0,i,de->ValueAsDouble(0,i) + wec->ValueAsDouble(0,i));
      
      chamber->SetValueAsDouble(0,i,chamber->ValueAsDouble(0,i) + wec->ValueAsDouble(0,i));
    }    
    delete wec;
  }  
  
  UpdateNumberOfEvents(&nevents);

}

//_____________________________________________________________________________
AliMUONTrackerData::AliMUONTrackerData(const char* name, const char* title,
                                       const AliMUONVStore& deOrBpValues, Int_t val)
: AliMUONVTrackerData(name,title),
fIsSingleEvent(kFALSE),
fChannelValues(0x0),
fManuValues(0x0),
fBusPatchValues(0x0),
fDEValues(0x0),
fChamberValues(0x0),
fPCBValues(0x0),
fDimension(0),
fNevents(0),
fDimensionNames(0x0),
fExternalDimensionNames(0x0),
fExternalDimension(0),
fHistogramming(0x0),
fHistos(0x0),
fXmin(0.0),
fXmax(0.0),
fIsChannelLevelEnabled(kFALSE),
fIsManuLevelEnabled(kFALSE),
fIsBustPatchLevelEnabled(kFALSE),
fIsPCBLevelEnabled(kFALSE),
fNofDDLs(0),
fNofEventsPerDDL(0x0)
{  
  /// ctor with values at the detection element OR bus patch level
  /// In this case, we force fIsChannelLevelEnabled = fIsManuLevelEnabled = kFALSE
  /// ctor
  
  if (deOrBpValues.GetSize()==0)
  {
    AliFatal("Cannot create a tracker data from nothing in that case !");
  }
  
  if ( !AliMpDDLStore::Instance(kFALSE) && !AliMpManuStore::Instance(kFALSE) )
  {
    AliError("Cannot work without (full) mapping");
    return;
  }
  
  if ( val == 1 ) 
  {
    BuildFromDEStore(deOrBpValues);
  }
  else if ( val == 2 )
  {
    BuildFromBPStore(deOrBpValues);
  }
  else
  {
    AliFatal("Wrong parameter. Must be 1 or 2");
  }
  
  
}

//_____________________________________________________________________________
void AliMUONTrackerData::BuildFromDEStore(const AliMUONVStore& deValues)
{
  /// Fill internals from a store of values at the detection element level
  
  TIter next(deValues.CreateIterator());
  AliMUONVCalibParam* m = static_cast<AliMUONVCalibParam*>(next());
  
  Int_t dimension = ( m->Dimension() - fgkExtraDimension - fgkVirtualExtraDimension ) / 2;
  
  fDimension = External2Internal(dimension)+fgkExtraDimension;
  
  fDimensionNames = new TObjArray(fDimension+fgkVirtualExtraDimension);
  fExternalDimensionNames = new TObjArray(dimension);
  fExternalDimension = dimension;
  fHistogramming = new Int_t[fExternalDimension];
  memset(fHistogramming,0,fExternalDimension*sizeof(Int_t)); // histogramming is off by default. Use MakeHistogramForDimension to turn it on.
  
  fExternalDimensionNames->SetOwner(kTRUE);
  fDimensionNames->SetOwner(kTRUE);  
  fDimensionNames->AddAt(new TObjString("occ"),IndexOfOccupancyDimension());
  fDimensionNames->AddAt(new TObjString("N"),IndexOfNumberDimension());
  fDimensionNames->AddAt(new TObjString("n"),NumberOfDimensions()-fgkVirtualExtraDimension);
  Clear();
  TArrayI nevents(AliDAQ::NumberOfDdls("MUONTRK"));
  AssertStores();
  
  next.Reset();
  AliMUONVCalibParam* external;
  
  while ( ( external = static_cast<AliMUONVCalibParam*>(next()) ) )
  {
    Int_t detElemId = external->ID0();
    
    AliMpDetElement* mpde = AliMpDDLStore::Instance()->GetDetElement(detElemId,kFALSE);
    
    if (!mpde)
    {
      AliError(Form("Got an invalid DE (%d) from external store",detElemId));
      continue;
    }
    
    Int_t chamberId = AliMpDEManager::GetChamberId(detElemId);
    AliMUONVCalibParam* chamber = ChamberParam(chamberId,kTRUE);
    AliMUONVCalibParam* de = DetectionElementParam(detElemId,kTRUE);
    
    AliMUONVCalibParam* wec = new AliMUONCalibParamND(external->Dimension()-1,1,detElemId,0,0);
    // as external, but without event count
    wec->SetValueAsDouble(0,0,external->ValueAsDouble(0,0));
    wec->SetValueAsDouble(0,1,external->ValueAsDouble(0,1));
    wec->SetValueAsDouble(0,2,external->ValueAsDouble(0,2));
    wec->SetValueAsDouble(0,3,external->ValueAsDouble(0,3));
    
    Int_t n1 = de->ValueAsInt(0,IndexOfNumberDimension());
    Int_t n2 = external->ValueAsInt(0,IndexOfNumberDimension());
    if  ( n1 != n2 )
    {
      AliError(Form("Incoherent number of dimensions for DE%d : %d vs %d",
                    detElemId,n1,n2));
      continue;
    }
    
    Int_t nevt = external->ValueAsInt(0,4);
    
    Int_t ddl = mpde->GetDdlId();
    
    if ( nevents[ddl] == 0 ) 
    {
      nevents[ddl] = nevt;
    }
    else
    {
      if ( nevents.At(ddl) != nevt ) 
      {
        AliError(Form("Nevt mismatch for DE %5d  DDL %d : %d vs %d",
                      detElemId,ddl,nevents.At(ddl),nevt));
        continue;
      }
    }
    
    for ( Int_t i = 0; i < wec->Dimension()-1; ++i ) 
    {
      de->SetValueAsDouble(0,i,de->ValueAsDouble(0,i) + wec->ValueAsDouble(0,i));
      
      chamber->SetValueAsDouble(0,i,chamber->ValueAsDouble(0,i) + wec->ValueAsDouble(0,i));
    }    
    delete wec;
  }  
  
  UpdateNumberOfEvents(&nevents);
}

//_____________________________________________________________________________
void AliMUONTrackerData::BuildFromBPStore(const AliMUONVStore& bpValues)
{
  /// Fill internals from a store of values at the bus patch level
  
  fIsBustPatchLevelEnabled = kTRUE;
  
  TIter next(bpValues.CreateIterator());
  AliMUONVCalibParam* m = static_cast<AliMUONVCalibParam*>(next());
  
  Int_t dimension = ( m->Dimension() - fgkExtraDimension - fgkVirtualExtraDimension ) / 2;
  
  fDimension = External2Internal(dimension)+fgkExtraDimension;
  
  fDimensionNames = new TObjArray(fDimension+fgkVirtualExtraDimension);
  fExternalDimensionNames = new TObjArray(dimension);
  fExternalDimension = dimension;
  fHistogramming = new Int_t[fExternalDimension];
  memset(fHistogramming,0,fExternalDimension*sizeof(Int_t)); // histogramming is off by default. Use MakeHistogramForDimension to turn it on.
  
  fExternalDimensionNames->SetOwner(kTRUE);
  fDimensionNames->SetOwner(kTRUE);  
  fDimensionNames->AddAt(new TObjString("occ"),IndexOfOccupancyDimension());
  fDimensionNames->AddAt(new TObjString("N"),IndexOfNumberDimension());
  fDimensionNames->AddAt(new TObjString("n"),NumberOfDimensions()-fgkVirtualExtraDimension);
  Clear();
  TArrayI nevents(AliDAQ::NumberOfDdls("MUONTRK"));
  AssertStores();
  
  next.Reset();
  AliMUONVCalibParam* external;
  
  while ( ( external = static_cast<AliMUONVCalibParam*>(next()) ) )
  {
    Int_t busPatchId = external->ID0();
    
    AliMpBusPatch* mpbp = AliMpDDLStore::Instance()->GetBusPatch(busPatchId,kFALSE);
    
    if (!mpbp)
    {
      AliError(Form("Got an invalid buspatchId (%d) from external store",busPatchId));
      continue;
    }
    
    Int_t detElemId = AliMpDDLStore::Instance()->GetDEfromBus(busPatchId);
    Int_t chamberId = AliMpDEManager::GetChamberId(detElemId);
    AliMUONVCalibParam* chamber = ChamberParam(chamberId,kTRUE);
    AliMUONVCalibParam* de = DetectionElementParam(detElemId,kTRUE);
    AliMUONVCalibParam* bp = BusPatchParam(busPatchId,kTRUE);
    
    AliMUONVCalibParam* wec = new AliMUONCalibParamND(external->Dimension()-1,1,busPatchId,0,0);
    // as external, but without event count
    wec->SetValueAsDouble(0,0,external->ValueAsDouble(0,0));
    wec->SetValueAsDouble(0,1,external->ValueAsDouble(0,1));
    wec->SetValueAsDouble(0,2,external->ValueAsDouble(0,2));
    wec->SetValueAsDouble(0,3,external->ValueAsDouble(0,3));
    
    Int_t n1 = bp->ValueAsInt(0,IndexOfNumberDimension());
    Int_t n2 = external->ValueAsInt(0,IndexOfNumberDimension());
    if  ( n1 != n2 )
    {
      AliError(Form("Incoherent number of dimensions for BP%d : %d vs %d",
                    busPatchId,n1,n2));
      continue;
    }
    
    Int_t nevt = external->ValueAsInt(0,4);
    
    Int_t ddl = mpbp->GetDdlId();
    
    if ( nevents[ddl] == 0 ) 
    {
      nevents[ddl] = nevt;
    }
    else
    {
      if ( nevents.At(ddl) != nevt ) 
      {
        AliError(Form("Nevt mismatch for BP %5d  DDL %d : %d vs %d",
                      busPatchId,ddl,nevents.At(ddl),nevt));
        continue;
      }
    }
    
    for ( Int_t i = 0; i < wec->Dimension()-1; ++i ) 
    {
      bp->SetValueAsDouble(0,i,bp->ValueAsDouble(0,i) + wec->ValueAsDouble(0,i));
      
      de->SetValueAsDouble(0,i,de->ValueAsDouble(0,i) + wec->ValueAsDouble(0,i));
      
      chamber->SetValueAsDouble(0,i,chamber->ValueAsDouble(0,i) + wec->ValueAsDouble(0,i));
    }    
    delete wec;
  }  
  
  UpdateNumberOfEvents(&nevents);  
}

//_____________________________________________________________________________
AliMUONTrackerData::AliMUONTrackerData(const char* name, const char* title,
                                       const AliMUONRejectList& rejectList)
: AliMUONVTrackerData(name,title),
fIsSingleEvent(kFALSE),
fChannelValues(0x0),
fManuValues(0x0),
fBusPatchValues(0x0),
fDEValues(0x0),
fChamberValues(0x0),
fPCBValues(0x0),
fDimension(External2Internal(1)+fgkExtraDimension),
fNevents(0),
fDimensionNames(new TObjArray(fDimension+fgkVirtualExtraDimension)),
fExternalDimensionNames(new TObjArray(1)),
fExternalDimension(1),
fHistogramming(new Int_t[fExternalDimension]),
fHistos(0x0),
fXmin(0.0),
fXmax(0.0),
fIsChannelLevelEnabled(kTRUE),
fIsManuLevelEnabled(kTRUE),
fIsBustPatchLevelEnabled(kTRUE),
fIsPCBLevelEnabled(kFALSE),
fNofDDLs(0),
fNofEventsPerDDL(0x0)
{  
  
   /// ctor with values from the given rejectlist
    
  if ( !AliMpDDLStore::Instance(kFALSE) && !AliMpManuStore::Instance(kFALSE) )
  {
    AliError("Cannot work without (full) mapping");
    return;
  }
  
  memset(fHistogramming,0,fExternalDimension*sizeof(Int_t)); // histogramming is off by default. Use MakeHistogramForDimension to turn it on.
  
  fExternalDimensionNames->SetOwner(kTRUE);
  fDimensionNames->SetOwner(kTRUE);  
  fDimensionNames->AddAt(new TObjString("occ"),IndexOfOccupancyDimension());
  fDimensionNames->AddAt(new TObjString("N"),IndexOfNumberDimension());
  fDimensionNames->AddAt(new TObjString("n"),NumberOfDimensions()-fgkVirtualExtraDimension);
  Clear();
  TArrayI nevents(AliDAQ::NumberOfDdls("MUONTRK"));
  AssertStores();
  
  
  for ( Int_t chamberId = 0; chamberId < AliMpConstants::NofChambers(); ++chamberId ) 
  {
//    AliMUONVCalibParam* chamber = ChamberParam(chamberId,kTRUE);

    // FIXME : update the chamber value ?
    
    AliMpDEIterator deit;
  
    deit.First(chamberId);

    while ( !deit.IsDone() ) 
    {
      AliMpDetElement* mpde = deit.CurrentDE();
      
      Int_t detElemId = mpde->GetId();

      AliMUONVCalibParam* de = DetectionElementParam(detElemId,kTRUE);

      DispatchValue(*de,0,rejectList.DetectionElementProbability(detElemId),0.0,mpde->NofChannels());
      
      for ( Int_t iBusPatch = 0; iBusPatch < mpde->GetNofBusPatches(); ++iBusPatch ) 
      {
        Int_t busPatchId = mpde->GetBusPatchId(iBusPatch);
        
        AliMUONVCalibParam* bp = BusPatchParam(busPatchId,kTRUE);

        AliMpBusPatch* mpbp = AliMpDDLStore::Instance()->GetBusPatch(busPatchId);

        Int_t nch(0);
        
        for ( Int_t iManu = 0 ;iManu < mpbp->GetNofManus(); ++iManu )
        {
          Int_t manuId = mpbp->GetManuId(iManu);
          
          nch += mpde->NofChannelsInManu(manuId);
          
          AliMUONVCalibParam* manu = ManuParam(detElemId,manuId,kTRUE);
          
          DispatchValue(*manu,0,rejectList.ManuProbability(detElemId,manuId),0.0,mpde->NofChannelsInManu(manuId));
          
          AliMUONVCalibParam* c = ChannelParam(detElemId,manuId);
          
          if (!c)
          {
            c = new AliMUONCalibParamND(Dimension(),
                                        AliMpConstants::ManuNofChannels(),
                                        detElemId,
                                        manuId,
                                        0.0);
            fChannelValues->Add(c);
          }
          
          for ( Int_t manuChannel = 0; manuChannel < AliMpConstants::ManuNofChannels(); ++manuChannel )
          {
            DispatchValue(*c,manuChannel,rejectList.ChannelProbability(detElemId,manuId,manuChannel),0.0,1);
          }
        }

        DispatchValue(*bp,0,rejectList.BusPatchProbability(busPatchId),0.0,nch);

      }
      
      deit.Next();
    }
 
  }
  
  
  SetDimensionName(0,"RejectProba");

  UpdateNumberOfEvents(0x0);
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
  delete fHistos;
  delete[] fNofEventsPerDDL;
}

//_____________________________________________________________________________
Bool_t
AliMUONTrackerData::Add(const AliMUONVStore& store, TArrayI* nevents)
{
  /// Add the given external store to our internal store
  return InternalAdd(store,nevents,kFALSE);
}

//_____________________________________________________________________________
Bool_t 
AliMUONTrackerData::Add(const AliMUONTrackerData& data)
{
  /// Add data to *this
  // We do this by looping on all VCalibParam stored in the various containers,
  // and simply adding the values there.
  // Same thing for the number of events per DDL.
  // Same thing for sparsehistograms, if we have some.

  // First cross check we have compatible objects.
  
  AliCodeTimerAuto("",0);
  
  if ( fIsChannelLevelEnabled != data.fIsChannelLevelEnabled ) 
  {
    AliError("Incompatible IsChannelLevelEnabled status");
    return kFALSE;
  }
  
  if ( fIsManuLevelEnabled != data.fIsManuLevelEnabled ) 
  {
    AliError("Incompatible IsManuLevelEnabled status");
    return kFALSE;
  }
  
  if ( fIsSingleEvent != data.fIsSingleEvent ) 
  {
    AliError("Incompatible IsSingleEvent status");
    return kFALSE;
  }
  
  if ( fDimension != data.fDimension || fExternalDimension != data.fExternalDimension ) 
  {
    AliError("Incompatible dimensions");
    return kFALSE;
  }

  if ( fNofDDLs != data.fNofDDLs ) 
  {
    AliError("Incompatible number of Ddls");
    return kFALSE;
  }
  
  if ( ( !fHistogramming && data.fHistogramming ) || ( fHistogramming && !data.fHistogramming ) 
      || fXmin != data.fXmin || fXmax != data.fXmax ) 
  {
    AliError(Form("Incompatible histogramming (%p vs %p) (xmax = %e vs %e ; xmin = %e vs %e)",
             fHistogramming,data.fHistogramming,fXmax,data.fXmax,fXmin,data.fXmin));
    return kFALSE;
  }

  if ( fHistogramming )
  {
    for ( Int_t i = 0; i < fExternalDimension; ++i ) 
    {
      if ( fHistogramming[i] != data.fHistogramming[i] )
      {
        AliError(Form("Incompatible histogramming for external dimension %d",i));
        return kFALSE;
      }
    }
  }
  
  // OK. Seems we have compatible objects, so we can proceed with the actual
  // merging...
  
  if ( data.fChannelValues ) 
  {
    Add2D(*(data.fChannelValues),*fChannelValues);
  }
  
  if ( data.fManuValues ) 
  {
    Add2D(*(data.fManuValues),*fManuValues);
  }
  
  if ( data.fPCBValues ) 
  {
    Add2D(*(data.fPCBValues),*fPCBValues);
  }
  
  if ( data.fBusPatchValues ) 
  {
    Add1D(*(data.fBusPatchValues),*fBusPatchValues);
  }
  
  if ( data.fDEValues ) 
  {
    Add1D(*(data.fDEValues),*fDEValues);
  }
  
  if ( data.fChamberValues ) 
  {
    Add1D(*(data.fChamberValues),*fChamberValues);
  }

  for ( Int_t i = 0; i < fNofDDLs; ++i ) 
  {
    fNofEventsPerDDL[i] += data.fNofEventsPerDDL[i];
  }
  
  if ( data.fHistos ) 
  {
    TIter nexthisto(data.fHistos->CreateIterator());
    AliMUONVStore* store;
    while ( ( store = static_cast<AliMUONVStore*>(nexthisto()) ) )
    {
      TIter ns(store->CreateIterator());
      AliMUONSparseHisto* h;
      while ( ( h = static_cast<AliMUONSparseHisto*>(ns()) ) )
      {
        AliMUONVStore* thisStore = static_cast<AliMUONVStore*>(fHistos->FindObject(store->GetUniqueID()));
        
        if (!thisStore)
        {
          thisStore = store->Create();
          thisStore->SetUniqueID(store->GetUniqueID());
          fHistos->Add(thisStore);
        }
        
        AliMUONSparseHisto* mine = static_cast<AliMUONSparseHisto*>(thisStore->FindObject(h->GetUniqueID()));
        
        if (!mine) 
        {
          thisStore->Add(h);
        }
        else
        {
          mine->Add(*h);
        }
      }
    }
  }
  
  for ( Int_t i = 0 ; i < fNofDDLs; ++i ) 
  {
    fNevents = TMath::Max(fNevents,fNofEventsPerDDL[i]);
  }
  
  return kTRUE;
}

//_____________________________________________________________________________
void 
AliMUONTrackerData::Add2D(const AliMUONVStore& src, AliMUONVStore& dest) const
{
  /// Add one 2d store to another
  
  TIter next(src.CreateIterator());
  AliMUONVCalibParam* p;
  
  while ( ( p = static_cast<AliMUONVCalibParam*>(next()) ) )
  {
    AliMUONVCalibParam* a = static_cast<AliMUONVCalibParam*>(dest.FindObject(p->ID0(),p->ID1()));
    
    if (!a)
    {
      dest.Add(static_cast<AliMUONVCalibParam*>(p->Clone()));
    }
    else
    {
      AddCalibParams(*p,*a);
    }
  }
}

//_____________________________________________________________________________
void 
AliMUONTrackerData::Add1D(const AliMUONVStore& src, AliMUONVStore& dest) const
{
  /// Add one 1d store to another
  
  TIter next(src.CreateIterator());
  AliMUONVCalibParam* p;
  
  while ( ( p = static_cast<AliMUONVCalibParam*>(next()) ) )
  {
    AliMUONVCalibParam* a = static_cast<AliMUONVCalibParam*>(dest.FindObject(p->GetUniqueID()));
    
    if (!a)
    {
      dest.Add(static_cast<AliMUONVCalibParam*>(p->Clone()));
    }
    else
    {
      AddCalibParams(*p,*a);
    }
  }
}

//_____________________________________________________________________________
void
AliMUONTrackerData::AddCalibParams(const AliMUONVCalibParam& src, AliMUONVCalibParam& dest) const
{
  /// Add src to dest
  for ( Int_t i = 0; i < src.Size(); ++i ) 
  {
    for ( Int_t j = 0; j < src.Dimension(); ++j )
    {
      dest.SetValueAsFloat(i,j,src.ValueAsFloat(i,j));
    }
  }
}

//_____________________________________________________________________________
Bool_t
AliMUONTrackerData::Replace(const AliMUONVStore& store)
{
  /// Replace our values by values from the given external store
  Bool_t rv = InternalAdd(store,0x0,kTRUE);
  AliMUONVTrackerData::Replace(store);
  return rv;
}

//_____________________________________________________________________________
Bool_t
AliMUONTrackerData::UpdateNumberOfEvents(TArrayI* nevents)
{
  /// Update the number of events
  
  if (!fNofDDLs)
  {
    fNofDDLs = AliDAQ::NumberOfDdls("MUONTRK");
    fNofEventsPerDDL = new Int_t[fNofDDLs];
    for ( Int_t i = 0; i < fNofDDLs; ++i ) 
    {
      fNofEventsPerDDL[i] = 0;
    }
  }
  
  if (nevents)
  {
    if (nevents->GetSize() != fNofDDLs ) 
    {
      AliError(Form("nof of ddl per event array size is incorrect : got %d, expecting %d",
                    nevents->GetSize(),fNofDDLs));
      return kFALSE;
    }
    
    for ( Int_t i = 0 ; i < fNofDDLs; ++i ) 
    {
      fNofEventsPerDDL[i] += nevents->At(i);
      fNevents = TMath::Max(fNevents,fNofEventsPerDDL[i]);
    }
  }
  else
  {
    for ( Int_t i = 0 ; i < fNofDDLs; ++i ) 
    {
      ++fNofEventsPerDDL[i];
      fNevents = TMath::Max(fNevents,fNofEventsPerDDL[i]);
    }
  }
  return kTRUE;
}

//_____________________________________________________________________________
void
AliMUONTrackerData::AssertStores()
{
  /// Insure our stores are allocated
  
  if (!fChamberValues)
  {
    Int_t numberOfBusPatches(0);
    Int_t numberOfDEs(0);
    
    // get number of bus patches and number of detection element
    // to initialize fBusPatchValues and fDEValues below
    
    if (!AliMpDDLStore::Instance(false))
    {
      AliMpCDB::LoadAll();
    }
    
    TIter next(AliMpDDLStore::Instance()->CreateBusPatchIterator());
    while ( next() ) ++numberOfBusPatches;
    AliMpDEIterator deIt;
    deIt.First();
    while (!deIt.IsDone())
    {
      ++numberOfDEs;
      deIt.Next();
    }
    
    if ( fIsChannelLevelEnabled ) 
    {
      fChannelValues = new AliMUON2DMap(kTRUE);
    }
    if  ( fIsManuLevelEnabled ) 
    {
      fManuValues = new AliMUON2DMap(kTRUE);
    }
    if ( fIsPCBLevelEnabled )
    {
      fPCBValues = new AliMUON2DMap(kFALSE);
    }
    if ( fIsBustPatchLevelEnabled )
    {
      fBusPatchValues = new AliMUON1DMap(numberOfBusPatches);
    }
    fDEValues = new AliMUON1DMap(numberOfDEs);
    fChamberValues = new AliMUON1DArray;
  }
}

//_____________________________________________________________________________
Bool_t
AliMUONTrackerData::InternalAdd(const AliMUONVStore& store, TArrayI* nevents, Bool_t replace)
{
  /// Add the given external store to our internal store
  
  AliCodeTimerAuto(GetName(),0);
    
  if ( !replace)
  {
    if ( IsSingleEvent() && NumberOfEvents(-1) == 1 ) 
    {
      AliError(Form("%s is supposed to be single event only",GetName()));
      return kFALSE;
    }  
  }

  UpdateNumberOfEvents(nevents);
  
  AssertStores();
  
  TIter next(store.CreateIterator());
  AliMUONVCalibParam* external;
  
  Int_t nk(2);
  
  if ( IsSingleEvent() ) nk = 1;

  while ( ( external = static_cast<AliMUONVCalibParam*>(next()) ) )
  {
    if ( external->Dimension() != ExternalDimension() )
    {
      AliError(Form("Incompatible dimensions %d vs %d",
                    external->Dimension(),ExternalDimension()));
      return kFALSE;
    }
    
    AliMUONVCalibParam* chamber(0x0);
    AliMUONVCalibParam* de(0x0);
    AliMUONVCalibParam* busPatch(0x0);
    AliMUONVCalibParam* pcb(0x0);
    AliMUONVCalibParam* manu(0x0);
    AliMUONVCalibParam* channel(0x0);
    AliMpDetElement* mpde(0x0);
    
    Int_t manuId = GetParts(external,chamber,de,busPatch,pcb,manu,channel,mpde);
    
    if ( manuId < 0 ) continue;
    
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
            FillHisto(detElemId,manuId,i,j,vext);
          }
          
          for ( Int_t k = 0; k < nk; ++k ) 
          {
            Double_t e = ( replace && channel ) ? channel->ValueAsDoubleFast(i,ix+k) : 0.0;
            
						if ( channel ) 
						{
							channel->SetValueAsDoubleFast(i,ix+k,channel->ValueAsDoubleFast(i,ix+k)-e+value[k]);
						}
						
            if (manu)
            {
              manu->SetValueAsDoubleFast(0,ix+k,manu->ValueAsDoubleFast(0,ix+k)-e+value[k]);            
            }
            
            busPatch->SetValueAsDoubleFast(0,ix+k,busPatch->ValueAsDoubleFast(0,ix+k)-e+value[k]);
            
            de->SetValueAsDoubleFast(0,ix+k,de->ValueAsDoubleFast(0,ix+k)-e+value[k]);
            
            chamber->SetValueAsDoubleFast(0,ix+k,chamber->ValueAsDoubleFast(0,ix+k)-e+value[k]);
            
            if ( pcb ) 
            {
              pcb->SetValueAsDoubleFast(0,ix+k,pcb->ValueAsDoubleFast(0,ix+k)-e+value[k]);
            }
          }
        }
        
        if ( validChannel && !replace )
        {
					if ( channel ) 
					{
						channel->SetValueAsDoubleFast(i,IndexOfOccupancyDimension(),
																					channel->ValueAsDoubleFast(i,IndexOfOccupancyDimension())+1.0);
					}
					
          if (manu)
          {
            manu->SetValueAsDoubleFast(0,IndexOfOccupancyDimension(),
                                       manu->ValueAsDoubleFast(0,IndexOfOccupancyDimension())+1.0);        
          }
          
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
  
  NumberOfEventsChanged();
  
  return kTRUE;
}

//_____________________________________________________________________________
Double_t 
AliMUONTrackerData::BusPatch(Int_t busPatchId, Int_t dim) const
{
  /// Return the value of a given buspatch for a given dimension
  /// or 0 if not existing
  AliMUONVCalibParam* param = BusPatchParam(busPatchId);
  return param ? Value(*param,0,dim,DdlIdFromBusPatchId(busPatchId)) : 0.0;
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
  
  AliCodeTimerAuto("",0);
  
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
  
  // FIXME: is the Value() correct wrt to number of events in the case of
  // chamber ? Or should we do something custom at the chamber level 
  // (as it spans several ddls) ?
  
  AliMUONVCalibParam* param = ChamberParam(chamberId);
  return param ? Value(*param,0,dim,DdlIdFromChamberId(chamberId)) : 0.0;
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
  
  AliCodeTimerAuto("",0);
  
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
  
  return param ? Value(*param,manuChannel,dim,DdlIdFromDetElemId(detElemId)) : 0.0;
}

//_____________________________________________________________________________
AliMUONVCalibParam* 
AliMUONTrackerData::ChannelParam(Int_t detElemId, Int_t manuId,
                                 const AliMUONVCalibParam* external) const
{
  /// Return (if it exist) the VCalibParam for a given manu
  
  AliMUONVCalibParam* param = fChannelValues ? static_cast<AliMUONVCalibParam*>
    (fChannelValues->FindObject(detElemId,manuId)) : 0x0 ;
  
  if (!param && external && fChannelValues)
  {
    param = CreateDouble(*external,detElemId,manuId);
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
  if ( fHistos ) fHistos->Clear();
  for ( Int_t i = 0; i < fNofDDLs; ++i ) 
  {
    fNofEventsPerDDL[i] = 0;
  }
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
AliMUONTrackerData::CreateDouble(const AliMUONVCalibParam& param, 
                                 Int_t detElemId, Int_t manuId) const
{
  /// Create a double version of VCalibParam, for internal use
  
  AliCodeTimerAuto("",0);
  
  AliMUONVCalibParam* c = new AliMUONCalibParamND(Dimension(),
                                                  param.Size(),
                                                  detElemId,
                                                  manuId,
                                                  0.0);
  
  AliMpDetElement* de = AliMpDDLStore::Instance()->GetDetElement(detElemId,manuId);
  
  for ( Int_t i = 0; i < c->Size(); ++i ) 
  {
    Double_t value(0.0);
    
    if ( de->IsConnectedChannel(manuId,i) ) value = 1.0;
      
    c->SetValueAsDouble(i,IndexOfNumberDimension(),value);
  }
  
  return c;
}

//_____________________________________________________________________________
Double_t 
AliMUONTrackerData::DetectionElement(Int_t detElemId, Int_t dim) const
{
  /// Return the value for a given detection element for a given dimension
  AliMUONVCalibParam* param = DetectionElementParam(detElemId);
  return param ? Value(*param,0,dim,DdlIdFromDetElemId(detElemId)) : 0.0;

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
  
  AliCodeTimerAuto("",0);
  
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
Int_t 
AliMUONTrackerData::DdlIdFromBusPatchId(Int_t buspatchid) const
{
  /// Get the "local" ddlid (0..19) of a given buspatch
  AliMpBusPatch* bp = AliMpDDLStore::Instance()->GetBusPatch(buspatchid);
  if (bp)
  {
    return bp->GetDdlId();
  }
  return -1;
}

//_____________________________________________________________________________
Int_t 
AliMUONTrackerData::DdlIdFromDetElemId(Int_t detelemid) const
{
  /// Get the "local" ddlid (0..19) of a given detection element
  AliMpDetElement* de = AliMpDDLStore::Instance()->GetDetElement(detelemid);
  if (de)
  {
    return de->GetDdlId();
  }
  return -1;
}

//_____________________________________________________________________________
Int_t 
AliMUONTrackerData::DdlIdFromChamberId(Int_t chamberid) const
{
  /// Get the "local" ddlid (0..19) of a given chamber
  /// This has no real meaning (as there are several ddls per chamber),
  /// so we take the ddlid where we got the max number of events
  
  AliMpDEIterator it;
  
  it.First(chamberid);
  Int_t n(0);
  Int_t d(-1);
  
  while (!it.IsDone())
  {
    Int_t detElemId = it.CurrentDEId();
    Int_t ddlId = DdlIdFromDetElemId(detElemId);
    if ( NumberOfEvents(ddlId) > n ) 
    {
      n = NumberOfEvents(ddlId);
      d = ddlId;
    }
    it.Next();
  }
  
  return d;
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
void 
AliMUONTrackerData::DisableChannelLevel()
{ 
  /// Disable the storing of data at channel level
  
  delete fChannelValues;
  fChannelValues = 0x0;
  fIsChannelLevelEnabled = kFALSE; 
}

//_____________________________________________________________________________
void 
AliMUONTrackerData::DisableManuLevel()
{ 
  /// Disable the storing of data at manu level (and below)
  
  DisableChannelLevel();
  delete fManuValues;
  fManuValues = 0x0;
  fIsManuLevelEnabled = kFALSE; 
}

//_____________________________________________________________________________
Int_t 
AliMUONTrackerData::External2Internal(Int_t index) const 
{
  /// From external to internal dimension
  return IsSingleEvent() ? index : index*2;
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
AliMUONTrackerData::FillHisto(Int_t detElemId, Int_t manuId, Int_t manuChannel,
                              Int_t dim, Double_t value)
{
  /// Fill histogram of a given channel
  
  AliMUONSparseHisto* h(0x0);
  
	if ( fIsChannelLevelEnabled ) 
	{
		h = GetChannelSparseHisto(detElemId, manuId, manuChannel,dim);
  }
  else if ( fIsManuLevelEnabled ) 
  {
    h = GetManuSparseHisto(detElemId,manuId,dim);
  }
  
  AliDebug(1,Form("DE %04d MANU %04d CH %02d dim %d value %e h %p",detElemId,manuId,manuChannel,dim,value,h));
  
  if (h)
  {
		h->Fill(static_cast<Int_t>(TMath::Nint(value)));
	}
}

//_____________________________________________________________________________
AliMUONSparseHisto*
AliMUONTrackerData::GetManuSparseHisto(Int_t detElemId, Int_t manuId, 
                                       Int_t dim) const
{
  /// Get histogram of a given manu
  
  if (!fHistos) return 0x0;
  
  AliMUON1DArray* m = static_cast<AliMUON1DArray*>(fHistos->FindObject(detElemId,manuId));
  if (!m) return 0x0;
  
  AliMUONSparseHisto* h = static_cast<AliMUONSparseHisto*>(m->FindObject(dim));
  
  return h;
}

//_____________________________________________________________________________
AliMUONSparseHisto*
AliMUONTrackerData::GetManuSparseHisto(Int_t detElemId, Int_t manuId, Int_t dim)
{
  /// Get histogram of a given manu. Create it if necessary
  
  if (!fHistos) fHistos = new AliMUON2DMap(kTRUE);
  
  AliMUON1DArray* m = static_cast<AliMUON1DArray*>(fHistos->FindObject(detElemId,manuId));
  if (!m)
  {
    m = new AliMUON1DArray(NumberOfDimensions());
    m->SetUniqueID( ( manuId << 16 ) | detElemId );
    fHistos->Add(m);
  }
    
  AliMUONSparseHisto* h = static_cast<AliMUONSparseHisto*>(m->FindObject(dim));
  if (!h)
  {
    h = new AliMUONSparseHisto(fXmin,fXmax);
    
    h->SetUniqueID(dim);
    
    m->Add(h);
  }
  
   return h;
}

//_____________________________________________________________________________
AliMUONSparseHisto*
AliMUONTrackerData::GetChannelSparseHisto(Int_t detElemId, Int_t manuId, 
                                          Int_t manuChannel, Int_t dim) const
{
  /// Get histogram of a given channel
  
  if (!fHistos) return 0x0;
  
  AliMUON1DMap* m = static_cast<AliMUON1DMap*>(fHistos->FindObject(detElemId,manuId));
  if (!m) return 0x0;
  
  UInt_t uid = ( manuChannel << 16 ) | dim;
  
  AliMUONSparseHisto* h = static_cast<AliMUONSparseHisto*>(m->FindObject(uid));
  
  return h;
}

//_____________________________________________________________________________
AliMUONSparseHisto*
AliMUONTrackerData::GetChannelSparseHisto(Int_t detElemId, Int_t manuId, 
                                          Int_t manuChannel, Int_t dim)
{
  /// Get histogram of a given channel. Create it if necessary
  
  if (!fHistos) fHistos = new AliMUON2DMap(kTRUE);
  
  AliMUON1DMap* m = static_cast<AliMUON1DMap*>(fHistos->FindObject(detElemId,manuId));
  if (!m)
  {
    m = new AliMUON1DMap(AliMpConstants::ManuNofChannels()); // start with only 1 dim
    m->SetUniqueID( ( manuId << 16 ) | detElemId );
    fHistos->Add(m);
  }
  
  UInt_t uid = ( manuChannel << 16 ) | dim;
  
  AliMUONSparseHisto* h = static_cast<AliMUONSparseHisto*>(m->FindObject(uid));
  if (!h)
  {
    h = new AliMUONSparseHisto(fXmin,fXmax);
    
    h->SetUniqueID(uid);
    
    m->Add(h);
  }

  return h;
}

//_____________________________________________________________________________
void
AliMUONTrackerData::GetDEManu(const AliMUONVCalibParam& param,
                              Int_t& detElemId, Int_t& manuId) const
{
  /// Tries to get (detElemId,manuId) of param
  
  // Load mapping manu store
  if ( ! AliMpCDB::LoadManuStore() ) {
    AliError("Could not access manu store from OCDB !");
    return;
  }

  if ( param.ID1() <= 0 ) 
  {
    // we (probably) get a manu serial number
    Int_t serial = param.ID0();
    MpPair_t pair = AliMpManuStore::Instance()->GetDetElemIdManu(serial);
    detElemId = AliMp::PairFirst(pair);
    manuId = AliMp::PairSecond(pair);
    if ( !detElemId ) 
    {
      AliDebug(1,Form("DE %d manuId %d from serial %d is not correct !",
                      detElemId,manuId,serial));
    }
  }
  else
  {
    // we get a (de,manu) pair
    detElemId = param.ID0();
    manuId = param.ID1();
  }
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
 
  chamber = de = busPatch = pcb = manu = channel = 0x0;
  mpde = 0x0;
  
  AliMpDDLStore* ddlStore = AliMpDDLStore::Instance();
  
  Int_t detElemId;
  Int_t manuId;
  
  GetDEManu(*external,detElemId,manuId);
  
  mpde = ddlStore->GetDetElement(detElemId,kFALSE);

  if (!mpde) // can happen if reading e.g. capacitances store where we have data for non-connected manus
  {
    return -1;
  }
  
  // explicitely check that de,manu is correct
  const AliMpVSegmentation* mpseg = AliMpSegmentation::Instance()->GetMpSegmentationByElectronics(detElemId, manuId,kFALSE);
  
  if (!mpseg)
  {
    return -1;
  }
  
  Int_t chamberId = AliMpDEManager::GetChamberId(detElemId);
    
  Int_t busPatchId = ddlStore->GetBusPatchId(detElemId,manuId);
  
  if ( busPatchId <= 0 )
  {
    return -1;
  }
  
  Int_t pcbIndex = -1;
  
  AliMp::StationType stationType = mpde->GetStationType();
  
  if ( stationType == AliMp::kStation345 ) 
  {
    AliMpDCSNamer namer("TRACKER");
    pcbIndex = namer.ManuId2PCBIndex(detElemId,manuId);
  }
  
  if ( fIsChannelLevelEnabled ) 
  {
    channel = ChannelParam(detElemId,manuId,external);
  }
  
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
  return param ? Value(*param,0,dim,DdlIdFromDetElemId(detElemId)) : 0.0;
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
  
  AliCodeTimerAuto("",0);
  
  AliMUONVCalibParam* manu = new AliMUONCalibParamND(Dimension(),1,detElemId,manuId,0.0);
  
  // set the number of channels in that manu
  
  AliMpDetElement* de = AliMpDDLStore::Instance()->GetDetElement(detElemId);
  
  manu->SetValueAsDouble(0,IndexOfNumberDimension(),de->NofChannelsInManu(manuId));
  
  return manu;
}

//_____________________________________________________________________________
Long64_t
AliMUONTrackerData::Merge(TCollection* list)
{
  /// Merge this with a list of AliMUONVTrackerData

  if (!list) return 0;
  
  if ( list->IsEmpty() ) return NumberOfEvents(-1);
  
  TIter next(list);
  const TObject* o(0x0);
  
  while ( ( o = next() ) )
  {
    const AliMUONTrackerData* data = dynamic_cast<const AliMUONTrackerData*>(o);
    if (!data)
    {
      AliError(Form("Object named %s is not an AliMUONTrackerData ! Skipping it",
                    o->GetName()));
    }
    else
    {
      Bool_t ok = Add(*data);
      if (!ok)
      {
        AliError("Got incompatible objects");
      }
    }
  }
  
  return NumberOfEvents(-1);
}

//_____________________________________________________________________________
Int_t 
AliMUONTrackerData::NumberOfDimensions() const
{
  /// Number of dimensions we're dealing with
  
  return fDimension + fgkVirtualExtraDimension; 
}

//_____________________________________________________________________________
Int_t
AliMUONTrackerData::NumberOfEvents(Int_t ddlNumber) const
{
  /// Get the number of events we've seen for a given DDL, or the max
  /// in case ddlNumber<0

  Int_t n(0);
  
  if ( fNofEventsPerDDL && ddlNumber >= 0 && ddlNumber < fNofDDLs )
  {
    n = fNofEventsPerDDL[ddlNumber];
  }
  else
  {
    // get the max number of events
    return fNevents;
  }
  
  return n;
}

//_____________________________________________________________________________
Double_t 
AliMUONTrackerData::PCB(Int_t detElemId, Int_t pcbIndex, Int_t dim) const
{
  /// Return the value of a given pcb for a given dimension

  AliMUONVCalibParam* param = PCBParam(detElemId,pcbIndex);
  
  return param ? Value(*param,0,dim,DdlIdFromDetElemId(detElemId)) : 0.0;
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
  
  AliCodeTimerAuto("",0);
  
  AliMpDCSNamer namer("TRACKER");
  
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
  
  if ( !fIsSingleEvent ) 
  {
    for ( Int_t i = 0; i < fNofDDLs; ++i ) 
    {
      cout << Form("DDL %04d Nevents=%10d",AliDAQ::DdlID("MUONTRK",i),fNofEventsPerDDL[i]) << endl;
    }
  }

	if ( !fIsChannelLevelEnabled ) 
	{
		cout << "Is not storing data at the channel level" << endl;
	}

  if ( !fIsManuLevelEnabled ) 
	{
		cout << "Is not storing data at the manu level" << endl;
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
  
  if ( sopt.Contains("CHANNEL") )
  {
    if ( fIsChannelLevelEnabled ) 
    {      
      if ( fChannelValues ) fChannelValues->Print(wildcard,opt);
    }
    else
    {
      AliWarning("You requested channel values, but they were not stored !");
    }
  }

  if ( sopt.Contains("MANU") )
  {
    if ( fIsManuLevelEnabled ) 
    {
      if ( fManuValues ) fManuValues->Print(wildcard,opt);
    }
    else
    {
      AliWarning("You requested manu values, but they were not stored !");
    }
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
    AliError(Form("%s : dimension %s : Index out of bounds : %d / %d",
                  GetName(),
                  name,index,fExternalDimension));
    return;
  }
  
  Int_t ix = External2Internal(index);
  
  if ( !IsSingleEvent() ) 
  {
    const char* prefix[] = { "mean", "sigma" };
  
    for ( Int_t i = 0; i < 2; ++i ) 
    {
      Int_t j = ix+i;
    
      SetInternalDimensionName(j,Form("%s of %s",prefix[i],name));
    }
  }
  else
  {
    SetInternalDimensionName(index,name);
  }
  
  SetExternalDimensionName(index,name);
}

//_____________________________________________________________________________
void 
AliMUONTrackerData::MakeHistogramForDimension(Int_t index, Bool_t value, Double_t xmin, Double_t xmax)
{
  /// decide to make histos for a given dimension
  if ( index >= ExternalDimension() ) 
  {
    AliError(Form("Index out of bounds : %d / %d",index,ExternalDimension()));
    return;
  }
  
  AliWarning(Form("Will %s make histogram for data %s index %d : that might ressemble a memory leak depending on the input data",
                  value ? "":"not", GetName(),index));
  fHistogramming[index] = value;
  fXmin = xmin;
  fXmax = xmax;
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
AliMUONTrackerData::Value(const AliMUONVCalibParam& param, Int_t i, 
                          Int_t dim, Int_t ddlId) const
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
  
  if ( dim == IndexOfOccupancyDimension() ) 
  {
    if ( ddlId < 0 ) AliError("Got a negative ddl id !");
    return occ/n/NumberOfEvents(ddlId);
  }
  
  Double_t value = param.ValueAsDouble(i,dim);
  
  if ( value >= AliMUONVCalibParam::InvalidFloatValue() ) return AliMUONVCalibParam::InvalidFloatValue();
  
  if ( TMath::Even(dim) || IsSingleEvent() ) 
  {
		Double_t x = value/occ;
		
		return ( TMath::Finite(x) ? x : 0.0 ) ;
  }
  else
  {
    Double_t nn = occ;
    
    if ( nn > 1.0 ) 
    {
      Double_t mean = param.ValueAsDouble(i,dim-1)/nn;
    
      return TMath::Sqrt(TMath::Abs((value-nn*mean*mean)/(nn-1.0)));
    }
    else
    {
      return 0.0;
    }
  }
  
  AliError("Why am I here ?");
  return 0.0;
}

//_____________________________________________________________________________
void 
AliMUONTrackerData::Streamer(TBuffer &R__b)
{
  /// Customized streamer                                                    
  
  if (R__b.IsReading()) {
    AliMUONTrackerData::Class()->ReadBuffer(R__b, this);
    if ( !fNofDDLs )
    {
      // backward compatible mode : we set number of events
      // per DDL to the total number of events (the only information
      // we had before version 7 of that class)
      delete[] fNofEventsPerDDL;
      fNofDDLs=20;
      fNofEventsPerDDL = new Int_t[fNofDDLs];
      for ( Int_t i = 0; i < fNofDDLs; ++i ) 
      {
        fNofEventsPerDDL[i] = fNevents;
      }
    }
  } 
  else {
    AliMUONTrackerData::Class()->WriteBuffer(R__b, this);
  }
}

//_____________________________________________________________________________
Bool_t 
AliMUONTrackerData::ExportAsASCIIOccupancyFile(const char* filename, Int_t runNumber) const
{
  /// Export only the occupancy part, in a format compatible with what
  /// the online occupancy DA is writing
  
  if ( ! AliMpDDLStore::Instance(kFALSE) )
  {
    AliError("Mapping not loaded. Cannot work");
    return kFALSE;
  }
  
  if (!fManuValues)
  {
    AliError("No manu values. Cannot work");
    return kFALSE;
  }
  
  ofstream out(filename);
  
  if (out.bad())
  {
    AliError(Form("Cannot create file %s",filename));
    return kFALSE;
  }
  
  out << "//===========================================================================" << endl;
  out << "//  Hit counter exported from $Id$" << endl;
  out << "//===========================================================================" << endl;
  out << "//" << endl;
  out << "//       * Run Number          : " << runNumber << endl;
  out << "//       * File Creation Date  : " << TTimeStamp().AsString("l") << endl;
  out << "//---------------------------------------------------------------------------" << endl;
  out << "//  BP   MANU  SUM_N  NEVENTS" << endl;
  out << "//---------------------------------------------------------------------------" << endl;
  
  TIter next(fManuValues->CreateIterator());
  AliMUONVCalibParam* manu;
  
  while ( ( manu = static_cast<AliMUONVCalibParam*>(next()) ) )
  {
    Int_t detElemId = manu->ID0();
    Int_t manuId = manu->ID1();
    Int_t busPatchId = AliMpDDLStore::Instance()->GetBusPatchId(detElemId,manuId);
    Int_t ddl = AliMpDDLStore::Instance()->GetDDLfromBus( busPatchId);
    if ( busPatchId < 0 || ddl < 0 )
    {
      AliError(Form("Got invalid (DE,manu,bp,ddl)=(%d,%d,%d,%d). Skipping it",detElemId,manuId,busPatchId,ddl));
      continue;
    }

    Int_t nevents = fNofEventsPerDDL[ddl];
    
    out << Form("%5d %5d %10d %10d",busPatchId,manuId,manu->ValueAsInt(0,IndexOfOccupancyDimension()),nevents) << endl;
  }
  
  out.close();
  return kTRUE;
}
  
//_____________________________________________________________________________
void AliMUONTrackerData::DispatchValue(AliMUONVCalibParam& param, 
                                       Int_t index,
                                       Double_t y, 
                                       Double_t ey,
                                       Int_t nchannels)
{
  /// fills the calibparam with a single value
  
  Double_t sumn = 1000.0; // or any value strictly above 1
  Double_t sumw = sumn*y;
  Double_t sumw2 = (sumn-1)*ey*ey+sumw*sumw/sumn;
  
  param.SetValueAsDouble(index,0,sumw);
  param.SetValueAsDouble(index,1,sumw2);
  param.SetValueAsDouble(index,2,sumn);
  param.SetValueAsDouble(index,3,nchannels);
  
}
