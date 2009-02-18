/// 
/// Macro to time the PadBy*** methods of AliMpVSegmentation implementation(s)
/// 
/// The output should ressemble this (output from a MacBook Pro 2.33 GHz on Feb, 18th, 2009)
///
/// I-AliCDBManager::Init: AliEn classes enabled in Root. AliCDBGrid factory registered.
/// I-AliCDBManager::SetDefaultStorage: Setting Default storage to: local://$ALICE_ROOT/OCDB
/// W-AliCDBManager::Get: Run number explicitly set in query: CDB cache temporarily disabled!
/// AliMpDCSNamer
///   ManuId2PCBIndex  R:0.0645s C:0.0400s (9676 slices)
///   ManuId2Sector  R:0.0434s C:0.0200s (7152 slices)
/// AliMpDetElement
///   AddManu 
///     R:1.6670s C:1.6900s (16828 slices)
///    slat R:0.8474s C:0.8600s (9676 slices)
///    st12 R:0.6383s C:0.6000s (7152 slices)
/// AliMpFastSegmentation
///   AliMpFastSegmentation 
///     AliMpSectorSegmentation R:0.0929s C:0.0900s (4 slices)
///     AliMpSlatSegmentation R:0.1142s C:0.1400s (38 slices)
/// General
///   timeMapping Load mapping R:2.9965s C:2.9000s (1 slices) (1)
/// ************************************
/// *    Row   * pI.fMemRe *     sname *
/// ************************************
/// *        0 *        30 *         0 *
/// *        1 *        76 *         1 * (2)
/// ************************************
/// AliMpMotifMap
///   GetMotifPosition  R:0.0795s C:0.1400s (7152 slices)
/// General
///   ByIndices 
///     HasPadByIndices-St12 R:0.1787s C:0.1700s (32 slices) (3)
///     HasPadByIndices-St345 R:0.1579s C:0.1900s (280 slices)
///     PadByIndices-St12 R:0.5637s C:0.5900s (32 slices)
///     PadByIndices-St345 R:0.5379s C:0.5100s (280 slices)
///   ByLocation 
///    HasPadByLocation-St12 R:0.0906s C:0.0800s (7152 slices) (4)
///    HasPadByLocation-St345 R:0.1216s C:0.1100s (9676 slices)
///    PadByLocation-St12 R:0.4506s C:0.4300s (7152 slices)
///    PadByLocation-St345 R:0.5874s C:0.5900s (9676 slices)
///   ByPosition 
///     PadByPosition-St12 R:7.6133s C:7.5700s (32 slices) (5)
///     PadByPosition-St345 R:2.3484s C:2.4300s (280 slices)
/// 
/// Interesting points in the output are : 
///
/// (1) : this is the total time it takes to (create and) load the mapping
/// (2) : row 1 - row 0 indicates the memory the mapping takes
/// (3) : the *PadByIndices* are timed here
/// (4) : the *PadByLocation* are timed here
/// (5) : the *PadByPosition* are timed here.
///
/// 3-4-5 : please note that the HasPadBy... methods are always faster, so
/// if you do not need the pad itself, but just to know if it's there, use
/// those. 
/// Note also that currently the PadByPosition is by far the slowest of the
/// 3 methods (Indices,Location,Position).
///
/// L. Aphecetche, Subatech
///

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
#include "TObjArray.h"
#include "TVector2.h"
#include "AliMpVPadIterator.h"
#include "AliSysInfo.h"
#include <TTree.h>

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
void ByPosition(const AliMpVSegmentation* seg, Int_t detElemId, const TObjArray& pads)
{
  /// Time the PadByPosition method
  
  Int_t stationId = StationId(detElemId);

  AliCodeTimerAutoGeneral(Form("PadByPosition-St%d",stationId));

  TIter next(&pads);
  AliMpPad* pad;
  
  while ( ( pad = static_cast<AliMpPad*>(next()) ) )
  {
    seg->PadByPosition(pad->Position(),kFALSE);
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
void timeMapping(Int_t nloop=1)
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
          
          TObjArray pads;
          pads.SetOwner(kTRUE);

          AliMpVPadIterator* pit = seg->CreateIterator();
        
          pit->First();
        
          while (!pit->IsDone())
          {
            AliMpPad pad = pit->CurrentItem();
            pads.Add(new AliMpPad(pad));
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