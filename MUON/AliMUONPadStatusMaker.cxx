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

//-----------------------------------------------------------------------------
/// \class AliMUONPadStatusMaker
///
/// Make a 2DStore of pad statuses, using different sources of information,
/// like pedestal values, gain values, and HV values.
///
/// \author Laurent Aphecetche
//-----------------------------------------------------------------------------

#include "AliMUONPadStatusMaker.h"

#include "AliCDBEntry.h"
#include "AliCDBManager.h"
#include "AliDCSValue.h"
#include "AliLog.h"
#include "AliMUON2DMap.h"
#include "AliMUON2DStoreValidator.h"
#include "AliMUONCalibParamNI.h"
#include "AliMUONCalibrationData.h"
#include "AliMUONHVNamer.h"
#include "AliMUONVCalibParam.h"
#include "AliMpArea.h"
#include "AliMpConstants.h"
#include "AliMpDEIterator.h"
#include "AliMpDEManager.h"
#include "AliMpIntPair.h"
#include "AliMpManuList.h"
#include "AliMpMotifMap.h"
#include "AliMpMotifPosition.h"
#include "AliMpPCB.h"
#include "AliMpPad.h"
#include "AliMpSector.h"
#include "AliMpSectorSegmentation.h"
#include "AliMpSegmentation.h"
#include "AliMpSlat.h"
#include "AliMpSlatSegmentation.h"
#include "AliMpStationType.h"
#include "AliMpVPadIterator.h"
#include <Riostream.h>
#include <TMap.h>
#include <TStopwatch.h>
#include <TString.h>


/// \cond CLASSIMP
ClassImp(AliMUONPadStatusMaker)
/// \endcond

//_____________________________________________________________________________
AliMUONPadStatusMaker::AliMUONPadStatusMaker(const AliMUONCalibrationData& calibData)
: fCalibrationData(calibData),
  fPedMeanLimits(0,4095),
  fPedSigmaLimits(0,4095),
  fHVSt12Limits(0,5000),
  fHVSt345Limits(0,5000)
{
   /// ctor
}

//_____________________________________________________________________________
AliMUONPadStatusMaker::~AliMUONPadStatusMaker()
{
  /// dtor.
}

//_____________________________________________________________________________
AliMUONVStore*
AliMUONPadStatusMaker::Combine(const AliMUONVStore& store1,
                               const AliMUONVStore& store2,
                               Int_t binShift) const
{
  /// Combine two status containers into one, shifting store2 status bits
  /// to the left by binShift before making an OR with store1.
  
  TStopwatch timer;
  timer.Start(kTRUE);
  
  AliMUONVStore* combined = static_cast<AliMUONVStore*>(store1.Clone());
  
  TIter next(store1.CreateIterator());
  AliMUONVCalibParam* param1;
  
  while ( ( param1 = static_cast<AliMUONVCalibParam*>(next()) ) )
  {
    Int_t detElemId = param1->ID0();
    Int_t manuId = param1->ID1();
    AliMUONVCalibParam* param2 = static_cast<AliMUONVCalibParam*>(store2.FindObject(detElemId,manuId));
    if (!param2)
    {
      AliWarning(Form("Could not get statuses for store2 for DE %d ManuId %d. Marking as missing.",
                    detElemId,manuId));
      param2 = static_cast<AliMUONVCalibParam*>(param1->Clone());
      for ( Int_t manuChannel = 0; manuChannel < param2->Size(); ++manuChannel )
      {
        param2->SetValueAsInt(manuChannel,0,kMissing);
      }      
    }
    AliMUONVCalibParam* paramCombined = static_cast<AliMUONVCalibParam*>(combined->FindObject(detElemId,manuId));
    if (!paramCombined)
    {
      paramCombined = static_cast<AliMUONVCalibParam*>(param2->Clone());
      combined->Add(paramCombined);
    }
    
    for ( Int_t manuChannel = 0; manuChannel < param1->Size(); ++manuChannel )
    {
      if ( AliMpManuList::DoesChannelExist(detElemId, manuId, manuChannel) )
      {
        Int_t status1(param1->ValueAsInt(manuChannel));
        Int_t status2(param2->ValueAsInt(manuChannel));
        
        Int_t status = status1 | (status2 << binShift);
        
        paramCombined->SetValueAsInt(manuChannel,0,status);
      }
    }
  }
  
  AliInfo("Timer:");
  StdoutToAliInfo(timer.Print(););
  
  return combined;
}

//_____________________________________________________________________________
AliMUONVStore* 
AliMUONPadStatusMaker::GeneratePadStatus(Int_t value)
{
  /// Generate a "fake" store, with all (detElemId,manuId) present,
  /// and containing all the same value
  
  AliMUONVStore* store = new AliMUON2DMap(true);
  
  TList* list = AliMpManuList::ManuList();
  
  AliMpIntPair* pair;
  
  TIter next(list);
  
  while ( ( pair = static_cast<AliMpIntPair*>(next()) ) ) 
  {
    Int_t detElemId = pair->GetFirst();
    Int_t manuId = pair->GetSecond();
    AliMUONVCalibParam* param = new AliMUONCalibParamNI(1,AliMpConstants::ManuNofChannels(),detElemId,manuId,value);
    store->Add(param);
  }
  
  delete list;
  
  return store;
}

//_____________________________________________________________________________
Bool_t 
AliMUONPadStatusMaker::GetSt12Status(const TMap& hvMap,
                                     Int_t detElemId, Int_t sector,
                                     Bool_t& hvChannelTooLow,
                                     Bool_t& hvChannelTooHigh,
                                     Bool_t& hvChannelON) const
{
  /// Get HV status for one HV sector of St12
  
  /// For a given PCB in a given DE, get the HV status (both the channel
  /// and the switch).
  /// Returns false if hv switch changed during the run.
  
  Bool_t error = kFALSE;
  hvChannelTooLow = kFALSE;
  hvChannelTooHigh = kFALSE;
  hvChannelON = kTRUE;
  
  AliMUONHVNamer hvNamer;
  
  TString hvChannel(hvNamer.DCSHVChannelName(detElemId,sector));
  
  TPair* hvPair = static_cast<TPair*>(hvMap.FindObject(hvChannel.Data()));
  if (!hvPair)
  {
    AliError(Form("Did not find expected alias (%s) for DE %d",
                  hvChannel.Data(),detElemId));  
    error = kTRUE;
  }
  else
  {
    TObjArray* values = static_cast<TObjArray*>(hvPair->Value());
    if (!values)
    {
      AliError(Form("Could not get values for alias %s",hvChannel.Data()));
      error = kTRUE;
    }
    else
    {
      // find out min and max value, and makes a cut
      Float_t hvMin(1E9);
      Float_t hvMax(0);
      TIter next(values);
      AliDCSValue* val;
      
      while ( ( val = static_cast<AliDCSValue*>(next()) ) )
      {
        Float_t hv = val->GetFloat();
        hvMin = TMath::Min(hv,hvMin);
        hvMax = TMath::Max(hv,hvMax);
      }
      
      float lowThreshold = fHVSt12Limits.X();
      float highThreshold = fHVSt12Limits.Y();
            
      if ( hvMin < lowThreshold ) hvChannelTooLow = kTRUE;
      if ( hvMax > highThreshold ) hvChannelTooHigh = kTRUE;
      if ( hvMin < 1 ) hvChannelON = kFALSE;
    }
  }
  
  return error;
}

//_____________________________________________________________________________
Bool_t 
AliMUONPadStatusMaker::GetSt345Status(const TMap& hvMap,
                                      Int_t detElemId, Int_t pcbIndex,
                                      Bool_t& hvChannelTooLow,
                                      Bool_t& hvChannelTooHigh,
                                      Bool_t& hvChannelON,
                                      Bool_t& hvSwitchON) const
{
  /// For a given PCB in a given DE, get the HV status (both the channel
  /// and the switch).
  /// Returns false if something goes wrong (in particular if 
  /// hv switch changed during the run).
  
  Bool_t error = kFALSE;
  hvChannelTooLow = kFALSE;
  hvChannelTooHigh = kFALSE;
  hvSwitchON = kTRUE;
  hvChannelON = kTRUE;
  
  AliMUONHVNamer hvNamer;
  
  TString hvChannel(hvNamer.DCSHVChannelName(detElemId));
  
  TPair* hvPair = static_cast<TPair*>(hvMap.FindObject(hvChannel.Data()));
  if (!hvPair)
  {
    AliError(Form("Did not find expected alias (%s) for DE %d",
                  hvChannel.Data(),detElemId));  
    error = kTRUE;
  }
  else
  {
    TObjArray* values = static_cast<TObjArray*>(hvPair->Value());
    if (!values)
    {
      AliError(Form("Could not get values for alias %s",hvChannel.Data()));
      error = kTRUE;
    }
    else
    {
      // find out min and max value, and makes a cut
      Float_t hvMin(1E9);
      Float_t hvMax(0);
      TIter next(values);
      AliDCSValue* val;
      
      while ( ( val = static_cast<AliDCSValue*>(next()) ) )
      {
        Float_t hv = val->GetFloat();
        hvMin = TMath::Min(hv,hvMin);
        hvMax = TMath::Max(hv,hvMax);
      }

      float lowThreshold = fHVSt345Limits.X();
      float highThreshold = fHVSt345Limits.Y();

      if ( hvMin < lowThreshold ) hvChannelTooLow = kTRUE;
      if ( hvMax > highThreshold ) hvChannelTooHigh = kTRUE;
      if ( hvMin < 1 ) hvChannelON = kFALSE;
    }
  }
  
  TString hvSwitch(hvNamer.DCSHVSwitchName(detElemId,pcbIndex));
  TPair* switchPair = static_cast<TPair*>(hvMap.FindObject(hvSwitch.Data()));
  if (!switchPair)
  {
    AliError(Form("Did not find expected alias (%s) for DE %d PCB %d",
                  hvSwitch.Data(),detElemId,pcbIndex));
    error = kTRUE;
  }
  else
  {
    TObjArray* values = static_cast<TObjArray*>(switchPair->Value());
    if (!values)
    {    
      AliError(Form("Could not get values for alias %s",hvSwitch.Data()));
      error = kTRUE;
    }
    else
    {
      // we'll count the number of ON/OFF for this pad, to insure
      // consistency (i.e. if status changed during the run, we should
      // at least notify this fact ;-) and hope it's not the norm)
      Int_t nTrue(0);
      Int_t nFalse(0);
      TIter next(values);
      AliDCSValue* val;
      
      while ( ( val = static_cast<AliDCSValue*>(next()) ) )
      {
        if ( val->GetBool() )
        {
          ++nTrue;
        }
        else
        {
          ++nFalse;
        }
      }
      
      if ( (nTrue>0 && nFalse>0) )
      {
        AliWarning(Form("Status of HV Switch %s changed during this run nTrue=%d nFalse=%d! Will consider it OFF",
                        hvSwitch.Data(),nTrue,nFalse));
        error = kTRUE;
      }
      
      if ( nFalse ) hvSwitchON = kFALSE;
    }
  }
  return error;
}

//_____________________________________________________________________________
AliMUONVStore* 
AliMUONPadStatusMaker::MakeGainStatus(const AliMUONVStore& /*gainValues*/) const
{
  /// FIXME: to be implemented
  AliWarning("Not implemented yet");
  return 0x0;
}

//_____________________________________________________________________________
AliMUONVStore* 
AliMUONPadStatusMaker::MakeHVStatus(const TMap& hvValues) const
{
  /// Scrutinize HV values and deduce an HV status for each pad
  
  TStopwatch timerSt12;
  TStopwatch timerSt345;
  
  timerSt12.Start(kTRUE);
  timerSt12.Stop();
  timerSt345.Start(kTRUE);
  timerSt345.Stop();
  
  AliMUONHVNamer hvNamer;
  
  AliMpDEIterator deIt;
  
  deIt.First();
  
  AliMUONVStore* hv = new AliMUON2DMap(kTRUE);
  
  while ( !deIt.IsDone() )
  {
    Int_t detElemId = deIt.CurrentDEId();
    
    switch ( AliMpDEManager::GetStationType(detElemId) )
    {
      case AliMp::kStation1:
      case AliMp::kStation2:
        timerSt12.Start(kFALSE);
        for ( int sector = 0; sector < 3; ++sector)
        {
          AliDebug(1,Form("detElemId %5d sector %d",detElemId,sector));

          Bool_t hvChannelTooLow, hvChannelTooHigh, hvChannelON;
          Bool_t error = GetSt12Status(hvValues,
                                       detElemId,sector,
                                       hvChannelTooLow,hvChannelTooHigh,
                                       hvChannelON);
          Int_t status = 0;
          if ( error ) status |= kHVError;
          if ( hvChannelTooLow ) status |= kHVTooLow;
          if ( hvChannelTooHigh ) status |= kHVTooHigh; 
          if ( !hvChannelON ) status |= kHVChannelOFF;
          SetStatusSt12(*hv,detElemId,sector,status);
          
        }
          timerSt12.Stop();
        break;
      case AliMp::kStation345:
      {
        timerSt345.Start(kFALSE);
        for ( Int_t pcbIndex = 0; pcbIndex < hvNamer.NumberOfPCBs(detElemId); ++pcbIndex)
        {
          AliDebug(1,Form("detElemId %5d pcbIndex %d",detElemId,pcbIndex));
          Bool_t hvChannelTooLow, hvChannelTooHigh, hvChannelON,hvSwitchON;
          Bool_t error = GetSt345Status(hvValues,
                                        detElemId,pcbIndex,
                                        hvChannelTooLow,hvChannelTooHigh,
                                        hvChannelON,hvSwitchON);
          Int_t status = 0;
          if ( error ) status |= kHVError;
          if ( hvChannelTooLow ) status |= kHVTooLow;
          if ( hvChannelTooHigh ) status |= kHVTooHigh; 
          if ( !hvSwitchON ) status |= kHVSwitchOFF; 
          if ( !hvChannelON) status |= kHVChannelOFF;
          SetStatusSt345(*hv,detElemId,pcbIndex,status);
        }
        timerSt345.Stop();
      }
        break;
      default:
        break;
    }
    deIt.Next();
  }
  
  AliInfo("St12 timer:");
  StdoutToAliInfo(timerSt12.Print(););
  AliInfo("St345 timer:");
  StdoutToAliInfo(timerSt345.Print(););
  
  return hv;
}

//_____________________________________________________________________________
AliMUONVStore* 
AliMUONPadStatusMaker::MakePedestalStatus(const AliMUONVStore& pedValues) const
{
  /// Assign a pedestal status to each pad
  
  TStopwatch timer;
  
  timer.Start(kTRUE);
  
  AliMUONVStore* pedStatuses = new AliMUON2DMap(kTRUE);
  
  TIter next(pedValues.CreateIterator());
  AliMUONVCalibParam* pedestals;
  Int_t nofManus(0);
  
  while ( ( pedestals = static_cast<AliMUONVCalibParam*>(next() ) ) )
  {
    Int_t detElemId = pedestals->ID0();
    Int_t manuId = pedestals->ID1();
    ++nofManus;
    for ( Int_t manuChannel = 0; manuChannel < pedestals->Size(); ++manuChannel )
    {
      Int_t status(0);
      if ( AliMpManuList::DoesChannelExist(detElemId, manuId, manuChannel) )
      {
        Float_t pedMean = pedestals->ValueAsFloat(manuChannel,0);
        Float_t pedSigma = pedestals->ValueAsFloat(manuChannel,1);
        if ( pedMean < fPedMeanLimits.X() ) status |= kPedMeanTooLow;
        if ( pedMean > fPedMeanLimits.Y() ) status |= kPedMeanTooHigh;
        if ( pedSigma < fPedSigmaLimits.X() ) status |= kPedSigmaTooLow;
        if ( pedSigma > fPedSigmaLimits.Y() ) status |= kPedSigmaTooHigh;
        if ( pedMean == 0 ) status |= kPedMeanZero;
        
        AliMUONVCalibParam* vStatus = 
          static_cast<AliMUONVCalibParam*>(pedStatuses->FindObject(detElemId,manuId));
        if ( !vStatus ) 
        {
          vStatus = new AliMUONCalibParamNI(1,AliMpConstants::ManuNofChannels(),detElemId,manuId,0);
          pedStatuses->Add(vStatus);
        }
        vStatus->SetValueAsInt(manuChannel,0,status);
      }
    }
  }
  
  AliInfo(Form("%d manus checked in :",nofManus));
  StdoutToAliInfo(timer.Print(););
  return pedStatuses;  
}

//_____________________________________________________________________________
AliMUONVStore* 
AliMUONPadStatusMaker::MakeStatus() const
{
  /// Read ped, gains and hv values from CDB, apply some Q&A and produces
  /// a combined status for each pad.

  TMap* hvValues = fCalibrationData.HV();
  AliMUONVStore* hvStatus(0x0);
  
  if (!hvValues)
  {
    AliError("Could not get HV values from CDB. Will create dummy ones and mark those as missing");
    hvStatus = GeneratePadStatus(kHVMissing);
  }
  else
  {
    hvStatus = MakeHVStatus(*hvValues);
  }
  
  AliMUONVStore* pedValues = fCalibrationData.Pedestals();
  AliMUONVStore* pedStatus(0x0);

  if (!pedValues)
  {
    AliError("Could not get pedestals values from CDB. Will create dummy ones and mark those as missing");
    pedStatus = GeneratePadStatus(kPedMissing);
  }
  else
  {
    pedStatus = MakePedestalStatus(*pedValues);
  }
  
  // FIXME: should do the same for gains as for hv and ped.    
  
  AliMUONVStore* status = Combine(*hvStatus,*pedStatus,8);
  
  delete hvStatus;
  delete pedStatus;
  
  // Insure we get all channels there (some or even all can be bad, but they
  // must be there somehow).
  
  AliMUON2DStoreValidator validator;
        
  TObjArray* a = validator.Validate(*status);
    
  if (a) 
  {
    // this should not happen.
    AliError("Status store not complete. Crash to follow soon...");
    StdoutToAliError(a->Print(););
    AliFatal("this should not happen at all!");
    delete status;
    status = 0x0;
  }
    
  return status;
}

//_____________________________________________________________________________
void
AliMUONPadStatusMaker::SetStatusSt12(AliMUONVStore& hvStatus,
                                     Int_t detElemId, 
                                     Int_t isector,
                                     Int_t status) const
{
  /// Flag all pads of detElemId (for St12) as bad.
  
  // FIXME: need a way to iterator on pads over a given HV sector for St12... 
  // we currently suppose that one sector is about a third of the chamber...
  // FIXME !! This has to be checked very carefully...
  
  const AliMp::CathodType kCathodes[] = { AliMp::kCath0, AliMp::kCath1 };
  
  for ( Int_t icathode = 0; icathode < 2; ++icathode )
  {
    const AliMpSectorSegmentation* seg = 
    static_cast<const AliMpSectorSegmentation*>(AliMpSegmentation::Instance()->GetMpSegmentation(detElemId,kCathodes[icathode]));
    const AliMpSector* sector = seg->GetSector();
    AliMpMotifMap* mMap = sector->GetMotifMap();
    TArrayI a;
    
    mMap->GetAllMotifPositionsIDs(a);
    
    TVector2 dim = seg->Dimensions();
    Double_t x = dim.X()*2;
    Double_t xmin = isector*x/3.0;
    Double_t xmax = xmin + x/3.0;   
    
    for ( Int_t i = 0; i < a.GetSize(); ++i ) 
    {
      AliMpMotifPosition* pos = mMap->FindMotifPosition(a[i]);
      Int_t manuId = pos->GetID();
      TVector2 position = pos->Position();
      if ( position.X() >= xmin && position.X() <= xmax) 
      {
        AliMUONVCalibParam* dead =
        static_cast<AliMUONVCalibParam*>(hvStatus.FindObject(detElemId,manuId));
        if (!dead)
        {
          dead = new AliMUONCalibParamNI(1,AliMpConstants::ManuNofChannels(),detElemId,manuId,status);
          hvStatus.Add(dead);
        }        
        else
        {
          // FIXME: this should really not happen, if we'd know really the
          // relationship between manuId and HV sector...
          // For the time being, let's leave it like that, for testing
          // purposes only. For production, this will have to be fixed.
          AliWarning("Please fixme.");
        }
      }
    }
  }  
}

//_____________________________________________________________________________
void
AliMUONPadStatusMaker::SetStatusSt345(AliMUONVStore& hvStatus,
                                      Int_t detElemId, Int_t pcbIndex,
                                      Int_t status) const
{
  /// Flag all pads of pcbIndex-th PCB of detElemId (for St345) as bad.
  
  const AliMp::CathodType kCathodes[] = { AliMp::kCath0, AliMp::kCath1 };
  
  for ( Int_t icathode = 0; icathode < 2; ++icathode )
  {
    const AliMpSlatSegmentation* seg = static_cast<const AliMpSlatSegmentation*>
    (AliMpSegmentation::Instance()->GetMpSegmentation(detElemId,kCathodes[icathode]));
    const AliMpSlat* slat = seg->Slat();
    const AliMpPCB* pcb = slat->GetPCB(pcbIndex);

    for ( Int_t i = 0; i < pcb->GetSize(); ++i ) 
    {
      AliMpMotifPosition* pos = pcb->GetMotifPosition(i);
      Int_t manuId = pos->GetID();
      AliMUONVCalibParam* dead = 
        static_cast<AliMUONVCalibParam*>(hvStatus.FindObject(detElemId,manuId));
      if (dead)
      {
        AliError(Form("dead is not null as expected from DE %d manuId %d",
                      detElemId,manuId));
      }
      if (!dead)
      {
        dead = new AliMUONCalibParamNI(1,AliMpConstants::ManuNofChannels(),detElemId,manuId,status);
        hvStatus.Add(dead);
      }
    }    
  }
}



