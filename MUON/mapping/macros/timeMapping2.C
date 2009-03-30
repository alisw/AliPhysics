// $Id$

/// An alternative macro to time the PadBy*** methods of AliMpVSegmentation 
/// implementation(s) which can handle AliMpPad not derived from TObject.
///
/// By L. Aphecetche, Subatech
/// Modified by I. Hrivnacova, IPN Orsay

#if !defined(__CINT__) || defined(__MAKECINT__)

#include "AliMpVSegmentation.h"
#include "AliMpCDB.h"
#include "AliMpSegmentation.h"
#include "AliMpPad.h"
#include "AliCodeTimer.h"
#include "AliMpDEManager.h"
#include "AliMpConstants.h"
#include "AliMpManuIterator.h"
#include "AliMpDEIterator.h"
#include "AliMpCathodType.h"
#include "AliMpStationType.h"
#include "AliMpVPadIterator.h"

#include "AliSysInfo.h"

#include <TObjArray.h>
#include <TVector2.h>
#include <TTree.h>

#include <vector>

// The line below should be commented if you want to try this macro
// on revision before 31082 (where AliMpVSegmentation did not have the HasPadBy...
// methods).

#define HASPAD

//______________________________________________________________________________
Int_t StationId(Int_t detElemId)
{
  switch ( 1 + AliMpDEManager::GetChamberId(detElemId) / 2 )
  {
    case 1:
    case 2:
      return 12;
      break;
    case 3:
    case 4:
    case 5:
      return 345;
      break;
    default:
      return -1;
  }
}

//______________________________________________________________________________
void ByPosition(const AliMpVSegmentation* seg, Int_t detElemId, 
                const std::vector<AliMpPad>& pads)
{
  /// Time the PadByPosition method
  
  Int_t stationId = StationId(detElemId);

  AliCodeTimerAutoGeneral(Form("PadByPosition-St%d",stationId));

  std::vector<AliMpPad>::const_iterator it; 
  for ( it = pads.begin(); it != pads.end(); it++ ) 
  {
    seg->PadByPosition(it->Position(),kFALSE);
  }
}

//______________________________________________________________________________
void ByIndices(const AliMpVSegmentation* seg, Int_t detElemId)
{
  /// Time the (Has)PadByIndices method
  
  Int_t stationId = StationId(detElemId);
  {
    AliCodeTimerAutoGeneral(Form("PadByIndices-St%d",stationId));
    
    for ( Int_t ix = 0; ix < seg->MaxPadIndexX(); ++ix )
    {
      for ( Int_t iy = 0; iy < seg->MaxPadIndexY(); ++iy )
      {
        seg->PadByIndices(AliMpIntPair(ix,iy),kFALSE);
      }
    }
  }
  
#ifdef HASPAD        
  {
    AliCodeTimerAutoGeneral(Form("HasPadByIndices-St%d",stationId));
    
    for ( Int_t ix = 0; ix < seg->MaxPadIndexX(); ++ix )
    {
      for ( Int_t iy = 0; iy < seg->MaxPadIndexY(); ++iy )
      {
        seg->HasPadByIndices(AliMpIntPair(ix,iy));
      }
    }
  }
#endif        
}

//______________________________________________________________________________
void ByLocation(const AliMpVSegmentation* seg, Int_t detElemId, Int_t manuId)
{
  /// Time the (Has)PadByLocation method
  
  Int_t stationId = StationId(detElemId);
  {
    AliCodeTimerAutoGeneral(Form("PadByLocation-St%d",stationId));
  
    for ( Int_t manuChannel = 0; manuChannel < AliMpConstants::ManuNofChannels(); ++manuChannel )
    {
      seg->PadByLocation(AliMpIntPair(manuId,manuChannel),kFALSE);
    }
  }

#ifdef HASPAD
  {
    AliCodeTimerAutoGeneral(Form("HasPadByLocation-St%d",stationId));
    
    for ( Int_t manuChannel = 0; manuChannel < AliMpConstants::ManuNofChannels(); ++manuChannel )
    {
      seg->HasPadByLocation(AliMpIntPair(manuId,manuChannel));
    }
  }
#endif      
  
}

//______________________________________________________________________________
void timeMapping2(Int_t nloop=1)
{
  AliCodeTimer::Instance()->Reset();

  {
    AliSysInfo::AddStamp("0");
    AliCodeTimerAutoGeneral("Load mapping");
    AliMpCDB::LoadDDLStore2();
    AliSysInfo::AddStamp("1");
    AliCodeTimer::Instance()->Print();
    TTree t;
    t.ReadFile("syswatch.log");
    t.Scan("pI.fMemResident:sname");
  }

  AliCodeTimer::Instance()->Reset();
  
  for ( Int_t i = 0; i < nloop; ++i )
  {
    AliMpManuIterator it;
    
    Int_t detElemId;
    Int_t manuId;
    
    while ( it.Next(detElemId,manuId) )
    {
      const AliMpVSegmentation* seg = AliMpSegmentation::Instance()->GetMpSegmentationByElectronics(detElemId,manuId);
      
      ByLocation(seg,detElemId,manuId);
    }
    
    AliMpDEIterator deit;
    
    deit.First();
    
    while (!deit.IsDone())
    {
      Int_t detElemId = deit.CurrentDEId();
      
      if ( AliMpDEManager::GetStationType(detElemId) != AliMp::kStationTrigger ) 
      {
        
        for ( Int_t cath = 0; cath < 2; ++cath )
        {
          const AliMpVSegmentation* seg = AliMpSegmentation::Instance()->GetMpSegmentation(detElemId,AliMp::GetCathodType(cath));

          ByIndices(seg,detElemId);
          
          //TObjArray pads;
          //pads.SetOwner(kTRUE);
          std::vector<AliMpPad> pads;

          AliMpVPadIterator* pit = seg->CreateIterator();
        
          pit->First();
        
          while (!pit->IsDone())
          {
            AliMpPad pad = pit->CurrentItem();
            //pads.Add(new AliMpPad(pad));
            pads.push_back(pad);
            pit->Next();
          }
          
          delete pit;
          
          ByPosition(seg, detElemId, pads);
        }
        
      }
      
      deit.Next();
    }
    
  }
  AliCodeTimer::Instance()->Print();
}

#endif
