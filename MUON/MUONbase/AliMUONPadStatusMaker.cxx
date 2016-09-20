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
/// like pedestal values, LV values, and HV values.
///
/// \author Laurent Aphecetche
//-----------------------------------------------------------------------------

#include "AliMUONPadStatusMaker.h"

#include "AliMUON2DMap.h"
#include "AliMUON2DStoreValidator.h"
#include "AliMUONCalibParamNI.h"
#include "AliMUONCalibrationData.h"
#include "AliMUONLogger.h"
#include "AliMUONRecoParam.h"
#include "AliMUONStringIntMap.h"
#include "AliMUONTrackerData.h"
#include "AliMUONVCalibParam.h"

#include "AliMpArea.h"
#include "AliMpArrayI.h"
#include "AliMpBusPatch.h"
#include "AliMpCDB.h"
#include "AliMpConstants.h"
#include "AliMpDDLStore.h"
#include "AliMpDEManager.h"
#include "AliMpDetElement.h"
#include "AliMpDCSNamer.h"
#include "AliMpManuIterator.h"
#include "AliMpManuUID.h"
#include "AliMpStationType.h"

#include "AliCDBEntry.h"
#include "AliCDBManager.h"
#include "AliCodeTimer.h"
#include "AliDCSValue.h"
#include "AliLog.h"

#include <Riostream.h>
#include <TArrayI.h>
#include <TExMap.h>
#include <TFile.h>
#include <TKey.h>
#include <TMap.h>
#include <TROOT.h>
#include <TString.h>
#include <TSystem.h>

using std::cout;
using std::endl;
/// \cond CLASSIMP
ClassImp(AliMUONPadStatusMaker)
/// \endcond

//_____________________________________________________________________________
AliMUONPadStatusMaker::AliMUONPadStatusMaker(const AliMUONCalibrationData& calibData)
: fkCalibrationData(calibData),
fPedMeanLimits(0,4095),
fPedSigmaLimits(0,4095),
fManuOccupancyLimits(0,1.0),
fBuspatchOccupancyLimits(0,1.0),
fDEOccupancyLimits(0,1.0),
fStatus(new AliMUON2DMap(true)),
fHV(0x0),
fLV(0x0),
fPedestals(calibData.Pedestals()),
fTrackerData(0x0),
fConfig(calibData.Config())
{
  /// ctor
  SetHVLimit(-1,0.0);
}

//_____________________________________________________________________________
AliMUONPadStatusMaker::~AliMUONPadStatusMaker()
{
  /// dtor.

  delete fStatus;
  delete fHV;
  delete fLV;
  delete fTrackerData;
}

//_____________________________________________________________________________
TString
AliMUONPadStatusMaker::AsString(Int_t status)
{
  /// return a human readable version of the integer status

  if ( status == 0 )
  {
    return "Brave New World";
  }

  Int_t pedStatus;
  Int_t lvStatus;
  Int_t hvStatus;
  Int_t occStatus;

  DecodeStatus(status,pedStatus,hvStatus,lvStatus,occStatus);

  TString s;

  if ( pedStatus & kPedMeanZero ) s += "& Ped Mean is Zero ";
  if ( pedStatus & kPedMeanTooLow ) s += "& Ped Mean Too Low ";
  if ( pedStatus & kPedMeanTooHigh ) s += "& Ped Mean Too High ";
  if ( pedStatus & kPedSigmaTooLow ) s += "& Ped Sigma Too Low ";
  if ( pedStatus & kPedSigmaTooHigh ) s += "& Ped Sigma Too High ";
  if ( pedStatus & kPedMissing ) s += "& Ped is missing ";

	if ( lvStatus & kLVTooLow ) s+="& LV is Too Low ";
	if ( lvStatus & kLVMissing ) s+="& LV is missing ";

	if ( hvStatus & kHVError ) s+="& HV is on error ";
	if ( hvStatus & kHVTooLow ) s+="& HV is Too Low ";
	if ( hvStatus & kHVTooHigh ) s+="& HV is Too High ";
	if ( hvStatus & kHVChannelOFF ) s+="& HV has channel OFF ";
	if ( hvStatus & kHVSwitchOFF ) s+="& HV has switch OFF ";
	if ( hvStatus & kHVMissing ) s+="& HV is missing ";

  if ( occStatus & kManuOccupancyTooHigh ) s+="& manu occupancy too high ";
  if ( occStatus & kManuOccupancyTooLow ) s+="& manu occupancy too low ";
  if ( occStatus & kBusPatchOccupancyTooHigh ) s+="& bus patch occupancy too high ";
  if ( occStatus & kBusPatchOccupancyTooLow ) s+="& bus patch occupancy too low ";
  if ( occStatus & kDEOccupancyTooHigh ) s+="& DE occupancy too high ";
  if ( occStatus & kDEOccupancyTooLow ) s+="& DE occupancy too low ";

  if ( occStatus & kBusPatchRemovedByPAR ) s+="& BusPatch removed during PAR";

  if ( s[0] == '&' ) s[0] = ' ';

  return s;
}

//_____________________________________________________________________________
TString
AliMUONPadStatusMaker::AsCondition(Int_t mask)
{
  /// return a human readable version of the mask's equivalent condition

  TString s(AsString(mask));

  s.ReplaceAll("&","|");

  return s;
}

//_____________________________________________________________________________
Int_t
AliMUONPadStatusMaker::BuildStatus(Int_t pedStatus,
                                   Int_t hvStatus,
                                   Int_t lvStatus,
                                   Int_t occStatus)
{
  /// Build a complete status from specific parts (ped,hv,lv)

  return ( hvStatus & 0xFF ) | ( ( pedStatus & 0xFF ) << 8 ) |
  ( ( lvStatus & 0xFF ) << 16 ) |
  ( ( occStatus & 0xFF ) << 24 ) ;
}

//_____________________________________________________________________________
Int_t AliMUONPadStatusMaker::CheckConfigConsistencyWithPedestalInformation(Int_t detElemId,
                                                                           Int_t manuId) const
{
  /// Check the consistency between the information from the MUON/Calib/Config and
  /// MUON/Calib/Pedestals objects.

  AliMUONVCalibParam* pedestals = static_cast<AliMUONVCalibParam*>(fPedestals->FindObject(detElemId,manuId));

  AliMUONVCalibParam* config = static_cast<AliMUONVCalibParam*>(fConfig->FindObject(detElemId,manuId));

  if ( pedestals == 0 && config == 0 )
  {
    /// manu missing both in config and pedestal run : that is expected
    return 0;
  }

  if ( config == 0 && pedestals )
  {
    // a manu present in the pedestal run disappeared in the configuration
    // that is happening if we removed a bus patch _during_ the run and then
    // issued a PAR (Pause And Reconfigure) to change the readout configuration
    //
    // So, that's normal if all the manus of the same buspatch are in the same case.
    // Let's check that...
    AliMpBusPatch* busPatch = AliMpDDLStore::Instance()->GetBusPatch(detElemId,manuId);
    Int_t n = busPatch->GetNofManus();
    Int_t missing(0);
    for ( Int_t i = 0; i < n; ++i )
    {
      Int_t manu = busPatch->GetManuId(i);
      if ( fConfig->FindObject(detElemId,manuId) == 0x0 ) ++missing;
    }
    if ( missing != n )
    {
      AliError(Form("Got an inconsistent state between config and pedestal information for DE %4d MANU %4d BUSPATCH %4d : not all the manus from this bus patch are missing in the configuration ? ",detElemId,manuId,busPatch->GetId()));
      return -1;
    }
    return 1;
  }

  if ( pedestals == 0 && config != 0 )
  {
    AliError(Form("Got an inconsistent state between config and pedestal information for DE %4d MANU %4d : got a configuration but no pedestal values ???",detElemId,manuId));
    return -2;
  }

  return 0;
}

//_____________________________________________________________________________
void
AliMUONPadStatusMaker::DecodeStatus(Int_t status,
                                    Int_t& pedStatus,
                                    Int_t& hvStatus,
                                    Int_t& lvStatus,
                                    Int_t& occStatus)
{
  /// Decode complete status into specific parts (ped,hv,lv)

  occStatus = ( status & 0xFF000000 ) >> 24;
  lvStatus = ( status & 0xFF0000 ) >> 16;
  pedStatus = ( status & 0xFF00 ) >> 8;
  hvStatus = (status & 0xFF);
}

//_____________________________________________________________________________
Bool_t
AliMUONPadStatusMaker::HVSt12Status(Int_t detElemId, Int_t sector,
                                    Bool_t& hvChannelTooLow,
                                    Bool_t& hvChannelTooHigh,
                                    Bool_t& hvChannelON) const
{
  /// Get HV status for one HV sector of St12

  /// For a given PCB in a given DE, get the HV status (both the channel
  /// and the switch).
  /// Returns false if hv switch changed during the run.

  AliCodeTimerAuto("",0)

  if (!fHV) return kFALSE;

  Bool_t error = kFALSE;
  hvChannelTooLow = kFALSE;
  hvChannelTooHigh = kFALSE;
  hvChannelON = kTRUE;

  Int_t chamberId = AliMpDEManager::GetChamberId(detElemId);

  AliMpDCSNamer hvNamer("TRACKER");

  TString hvChannel(hvNamer.DCSAliasName(detElemId,sector));

  TMap* hvMap = fkCalibrationData.HV();
  TPair* hvPair = static_cast<TPair*>(hvMap->FindObject(hvChannel.Data()));
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
      // find out min value, and makes a cut
      Float_t hvMin(1E9);
      TIter next(values);
      AliDCSValue* val;

      while ( ( val = static_cast<AliDCSValue*>(next()) ) )
      {
        Float_t hv = val->GetFloat();
        hvMin = TMath::Min(hv,hvMin);
      }

      float lowThreshold = HVLimit(chamberId);

      if ( hvMin < lowThreshold ) hvChannelTooLow = kTRUE;
      if ( hvMin < hvNamer.TrackerHVOFF() ) hvChannelON = kFALSE;
    }
  }

  return error;
}

//_____________________________________________________________________________
Float_t
AliMUONPadStatusMaker::SwitchValue(const TObjArray& dcsArray)
{
  /// Loop over the dcs value for a single switch to decide whether
  /// we should consider it on or off

  // we'll count the number of ON/OFF for this pad, to insure
  // consistency (i.e. if status changed during the run, we should
  // at least notify this fact ;-) and hope it's not the norm)
  Int_t nTrue(0);
  Int_t nFalse(0);
  TIter next(&dcsArray);
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
    // change of state during the run, consider it off
    return 0.0;
  }

  if ( nFalse )
  {
    /// switch = FALSE means the HV was flowding up to the PCB.
    /// i.e. switch = FALSE = ON
    return 1.0;
  }

  return 0.0;
}

//_____________________________________________________________________________
Bool_t
AliMUONPadStatusMaker::HVSt345Status(Int_t detElemId, Int_t pcbIndex,
                                     Bool_t& hvChannelTooLow,
                                     Bool_t& hvChannelTooHigh,
                                     Bool_t& hvChannelON,
                                     Bool_t& hvSwitchON) const
{
  /// For a given PCB in a given DE, get the HV status (both the channel
  /// and the switch).
  /// Returns false if something goes wrong (in particular if
  /// hv switch changed during the run).

  AliCodeTimerAuto("",0)

  if (!fHV) return kFALSE;

  Bool_t error = kFALSE;
  hvChannelTooLow = kFALSE;
  hvChannelTooHigh = kFALSE;
  hvSwitchON = kTRUE;
  hvChannelON = kTRUE;

  AliMpDCSNamer hvNamer("TRACKER");

  Int_t chamberId = AliMpDEManager::GetChamberId(detElemId);

  TString hvChannel(hvNamer.DCSAliasName(detElemId));

  TMap* hvMap = fkCalibrationData.HV();

  TPair* hvPair = static_cast<TPair*>(hvMap->FindObject(hvChannel.Data()));
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
      // find out min value, and makes a cut
      Float_t hvMin(1E9);
      TIter next(values);
      AliDCSValue* val;

      while ( ( val = static_cast<AliDCSValue*>(next()) ) )
      {
        Float_t hv = val->GetFloat();
        hvMin = TMath::Min(hv,hvMin);
      }

      float lowThreshold = HVLimit(chamberId);

      if ( hvMin < lowThreshold ) hvChannelTooLow = kTRUE;
      if ( hvMin < hvNamer.TrackerHVOFF() ) hvChannelON = kFALSE;
    }
  }

  TString hvSwitch(hvNamer.DCSSwitchAliasName(detElemId,pcbIndex));
  TPair* switchPair = static_cast<TPair*>(hvMap->FindObject(hvSwitch.Data()));
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
      Float_t sv = SwitchValue(*values);
      if ( sv < 0.99 ) hvSwitchON = kFALSE;
    }
  }
  return error;
}

//_____________________________________________________________________________
Int_t
AliMUONPadStatusMaker::HVStatus(Int_t detElemId, Int_t manuId) const
{
  /// Get HV status of one manu

  AliCodeTimerAuto("",0)

  if ( !InternalHV() ) return kMissing;

  Long_t lint = InternalHV()->GetValue(AliMpManuUID::BuildUniqueID(detElemId,manuId));

  if ( lint )
  {
    return (Int_t)(lint - 1);
  }

  Int_t status(0);

  AliMpDCSNamer hvNamer("TRACKER");

  switch ( AliMpDEManager::GetStationType(detElemId) )
  {
    case AliMp::kStation12:
    {
      int sector = hvNamer.ManuId2Sector(detElemId,manuId);
      if ( sector >= 0 )
      {
        Bool_t hvChannelTooLow, hvChannelTooHigh, hvChannelON;
        Bool_t error = HVSt12Status(detElemId,sector,
                                    hvChannelTooLow,
                                    hvChannelTooHigh,
                                    hvChannelON);
        if ( error ) status |= kHVError;
        if ( hvChannelTooLow ) status |= kHVTooLow;
        if ( hvChannelTooHigh ) status |= kHVTooHigh;
        if ( !hvChannelON ) status |= kHVChannelOFF;
        // assign this status to all the other manus handled by the same HV channel
        SetHVStatus(detElemId,sector,status);
      }
    }
      break;
    case AliMp::kStation345:
    {
      int pcbIndex = hvNamer.ManuId2PCBIndex(detElemId,manuId);
      if ( pcbIndex >= 0 )
      {
        Bool_t hvChannelTooLow, hvChannelTooHigh, hvChannelON,hvSwitchON;
        Bool_t error = HVSt345Status(detElemId,pcbIndex,
                                     hvChannelTooLow,hvChannelTooHigh,
                                     hvChannelON,hvSwitchON);
        if ( error ) status |= kHVError;
        if ( hvChannelTooLow ) status |= kHVTooLow;
        if ( hvChannelTooHigh ) status |= kHVTooHigh;
        if ( !hvSwitchON ) status |= kHVSwitchOFF;
        if ( !hvChannelON) status |= kHVChannelOFF;
        // assign this status to all the other manus handled by the same HV channel
        SetHVStatus(detElemId,pcbIndex,status);
      }
    }
      break;
    default:
      break;
  }

  return status;
}

//_____________________________________________________________________________
Int_t AliMUONPadStatusMaker::LVStatus(Int_t detElemId, Int_t manuId) const
{
  /// Get LV status of one detection element

  AliCodeTimerAuto("",0)

  if ( !InternalLV() ) return kMissing;

  AliMp::PlaneType planeType = AliMp::kBendingPlane;

  AliMp::StationType stationType = AliMpDEManager::GetStationType(detElemId);

  if ( stationType == AliMp::kStation12 && ( manuId & AliMpConstants::ManuMask(AliMp::kNonBendingPlane) ) ) {
    // for St12 we need to know the plane type, for St345 we don't care
    planeType = AliMp::kNonBendingPlane;
  }

  Long_t lint = InternalLV()->GetValue(AliMpManuUID::BuildUniqueID(detElemId,planeType));

  if ( lint )
  {
    return (Int_t)(lint - 1);
  }

  Int_t status(0);
  Bool_t error(kFALSE);

  AliMpDCSNamer lvNamer("TRACKER");

  Int_t voltageType[] = { -1, 0, 1 };

  const float lowThreshold[] = { 1.5, 2.3, 1.5 };

  // we loop on the 3 voltages, and if any of them is below nominal,
  // we consider the LV group as off

  TMap* lvMap = fkCalibrationData.LV();

  Bool_t lvChannelON = kTRUE;

  for ( Int_t i = 0; i < 3; ++i )
  {
    TString lvGroup(lvNamer.DCSMCHLVAliasName(detElemId,voltageType[i],planeType));

    TPair* lvPair = static_cast<TPair*>(lvMap->FindObject(lvGroup.Data()));
    if (!lvPair)
    {
      AliError(Form("Did not find expected alias (%s) for DE %d",
                  lvGroup.Data(),detElemId));
      error = kTRUE;
    }
    else
    {
      TObjArray* values = static_cast<TObjArray*>(lvPair->Value());
      if (!values)
      {
        AliError(Form("Could not get values for alias %s",lvGroup.Data()));
        error = kTRUE;
      }
      else
      {
        // find out min value, and makes a cut
        Float_t lvMin(1E9);
        TIter next(values);
        AliDCSValue* val;

        while ( ( val = static_cast<AliDCSValue*>(next()) ) )
        {
          Float_t lv = val->GetFloat();
          lvMin = TMath::Min(lv,lvMin);
        }

        if ( lvMin < lowThreshold[i] ) lvChannelON = kFALSE;
      }
    }
  }

  if (!lvChannelON) {
    status |= kLVTooLow;
  }

  InternalLV()->Add(AliMpManuUID::BuildUniqueID(detElemId,planeType),status+1);

  return status;
}

//_____________________________________________________________________________
AliMUONVCalibParam*
AliMUONPadStatusMaker::Neighbours(Int_t detElemId, Int_t manuId) const
{
  /// Get the neighbours parameters for a given manu
  AliMUONVStore* neighbourStore = fkCalibrationData.Neighbours();
  return static_cast<AliMUONVCalibParam*>(neighbourStore->FindObject(detElemId,manuId));
}

//_____________________________________________________________________________
AliMUONVStore*
AliMUONPadStatusMaker::NeighboursStore() const
{
  /// Return the store containing all the neighbours
  return fkCalibrationData.Neighbours();
}

//_____________________________________________________________________________
AliMUONVCalibParam*
AliMUONPadStatusMaker::ComputeStatus(Int_t detElemId, Int_t manuId) const
{
  /// Compute the status of a given manu, using all available information,
  /// i.e. pedestals, LV, and HV

  AliMUONVCalibParam* param = new AliMUONCalibParamNI(1,AliMpConstants::ManuNofChannels(),detElemId,manuId,-1);
  fStatus->Add(param);

  AliMUONVCalibParam* pedestals = static_cast<AliMUONVCalibParam*>(fPedestals->FindObject(detElemId,manuId));

  Int_t hvStatus = HVStatus(detElemId,manuId);
  Int_t lvStatus = LVStatus(detElemId,manuId);

  Int_t occStatus = OccupancyStatus(detElemId,manuId);

  for ( Int_t manuChannel = 0; manuChannel < param->Size(); ++manuChannel )
  {
    Int_t pedStatus(0);

    if (pedestals)
    {
      Float_t pedMean = pedestals->ValueAsFloatFast(manuChannel,0);
      Float_t pedSigma = pedestals->ValueAsFloatFast(manuChannel,1);
      if ( pedMean < fPedMeanLimits.X() ) pedStatus |= kPedMeanTooLow;
      else if ( pedMean > fPedMeanLimits.Y() ) pedStatus |= kPedMeanTooHigh;
      if ( pedSigma < fPedSigmaLimits.X() ) pedStatus |= kPedSigmaTooLow;
      else if ( pedSigma > fPedSigmaLimits.Y() ) pedStatus |= kPedSigmaTooHigh;
      if ( pedMean == 0 ) pedStatus |= kPedMeanZero;
    }
    else
    {
      pedStatus = kPedMissing;
    }

    Int_t status = BuildStatus(pedStatus,hvStatus,lvStatus,occStatus);

    param->SetValueAsIntFast(manuChannel,0,status);
  }

  return param;
}

//_____________________________________________________________________________
Int_t
AliMUONPadStatusMaker::OccupancyStatus(Int_t detElemId, Int_t manuId) const
{
  /// Get the "other" status for a given manu

  Int_t rv(0);

  if ( InternalTrackerData() )
  {
    const Int_t occIndex = 2;

    Double_t occ = InternalTrackerData()->DetectionElement(detElemId,occIndex);

    if ( occ <= fDEOccupancyLimits.X() )
    {
      rv |= kDEOccupancyTooLow;
    }
    else if ( occ > fDEOccupancyLimits.Y() )
    {
      rv |= kDEOccupancyTooHigh;
    }

    Int_t busPatchId = AliMpDDLStore::Instance()->GetBusPatchId(detElemId,manuId);

    occ = InternalTrackerData()->BusPatch(busPatchId,occIndex);

    if ( occ <= fBuspatchOccupancyLimits.X() )
    {
      rv |= kBusPatchOccupancyTooLow;
    }
    else if ( occ > fBuspatchOccupancyLimits.Y() )
    {
      rv |= kBusPatchOccupancyTooHigh;
    }

    occ = InternalTrackerData()->Manu(detElemId,manuId,occIndex);

    if ( occ <= fManuOccupancyLimits.X() )
    {
      rv |= kManuOccupancyTooLow;
    }
    else if ( occ > fManuOccupancyLimits.Y() )
    {
      rv |= kManuOccupancyTooHigh;
    }
  }

  Int_t config = CheckConfigConsistencyWithPedestalInformation(detElemId,manuId);

  if (config==1)
  {
    Int_t bpid = AliMpDDLStore::Instance()->GetBusPatchId(detElemId,manuId);

    AliWarning(Form("BusPatchRemovedByPAR : BP %4d",bpid));
    rv |= kBusPatchRemovedByPAR;
  }

  return rv;
}

//_____________________________________________________________________________
AliMUONVCalibParam*
AliMUONPadStatusMaker::PadStatus(Int_t detElemId, Int_t manuId) const
{
  /// Get the status container for a given manu

  AliMUONVCalibParam* param = static_cast<AliMUONVCalibParam*>(fStatus->FindObject(detElemId,manuId));
  if (!param)
  {
    // not already there, so compute it now
    AliCodeTimerAuto("ComputeStatus",0);
    param = ComputeStatus(detElemId,manuId);
  }
  return param;
}

//_____________________________________________________________________________
Int_t
AliMUONPadStatusMaker::PadStatus(Int_t detElemId, Int_t manuId, Int_t manuChannel) const
{
  /// Get the status for a given channel

  AliMUONVCalibParam* param = static_cast<AliMUONVCalibParam*>(fStatus->FindObject(detElemId,manuId));
  if (!param)
  {
    // not already there, so compute it now
    param = ComputeStatus(detElemId,manuId);
  }
  return param->ValueAsInt(manuChannel,0);
}

//_____________________________________________________________________________
void
AliMUONPadStatusMaker::SetHVStatus(Int_t detElemId, Int_t index, Int_t status) const
{
  /// Assign status to all manus in a given HV "zone" (defined by index, meaning
  /// is different thing from St12 and St345)

  AliCodeTimerAuto("",0)

  AliMpDetElement* de = AliMpDDLStore::Instance()->GetDetElement(detElemId);

  const AliMpArrayI* manus = de->ManusForHV(index);

  for ( Int_t i = 0; i < manus->GetSize(); ++ i )
  {
    Int_t manuId = manus->GetValue(i);
    InternalHV()->Add(AliMpManuUID::BuildUniqueID(detElemId,manuId),status + 1);
  }
}

//_____________________________________________________________________________
Double_t
AliMUONPadStatusMaker::HVLimit(Int_t chamberId) const
{
  /// Get HV limit for a given chamber
  if ( chamberId >=0 && chamberId < 10 )
  {
    return fHVLimit[chamberId];
  }
  return 0.0;
}

//_____________________________________________________________________________
void
AliMUONPadStatusMaker::SetHVLimit(Int_t chamberId, Double_t hv)
{
  /// Set hv limit for a given chamber (or all if chamberId==-1)

  if ( chamberId == -1 )
  {
    for ( Int_t i = 0; i < 10; ++i )
    {
      fHVLimit[i] = hv;
    }
  }
  else if ( chamberId >= 0 && chamberId < 10 )
  {
    fHVLimit[chamberId]=hv;
  }
  else
  {
    AliError(Form("chamberId=%d is invalid",chamberId));
  }
}

//_____________________________________________________________________________
void
AliMUONPadStatusMaker::SetLimits(const AliMUONRecoParam& recoParams)
{
  /// Set the limits from the recoparam

  for ( int i = 0; i < 10; ++i )
  {
    SetHVLimit(i,recoParams.HVLimit(i));
  }

  SetPedMeanLimits(recoParams.PedMeanLowLimit(),recoParams.PedMeanHighLimit());
  SetPedSigmaLimits(recoParams.PedSigmaLowLimit(),recoParams.PedSigmaHighLimit());

  SetManuOccupancyLimits(recoParams.ManuOccupancyLowLimit(),recoParams.ManuOccupancyHighLimit());
  SetBuspatchOccupancyLimits(recoParams.BuspatchOccupancyLowLimit(),recoParams.BuspatchOccupancyHighLimit());
  SetDEOccupancyLimits(recoParams.DEOccupancyLowLimit(),recoParams.DEOccupancyHighLimit());
}

//_____________________________________________________________________________
void
AliMUONPadStatusMaker::Report(UInt_t mask)
{
  /// Report the number of bad pads, according to the mask,
  /// and the various reasons why they are bad (with occurence rates)

  AliInfo("");
  AliCodeTimerAuto("",0);

  AliMUONLogger log(1064008);

  Int_t nBadPads(0);
  Int_t nPads(0);

  AliMpManuIterator it;

  Int_t detElemId, manuId;

  while ( it.Next(detElemId,manuId) )
  {
    AliMpDetElement* de = AliMpDDLStore::Instance()->GetDetElement(detElemId);

    for ( Int_t i = 0; i < AliMpConstants::ManuNofChannels(); ++i )
    {
      if ( de->IsConnectedChannel(manuId,i) )
      {
        ++nPads;

        Int_t status = PadStatus(detElemId,manuId,i);

        if ( mask && ( status & mask) ) // note that if mask == 0, all pads are good...
        {
          ++nBadPads;
          log.Log(AsString(status));
        }
      }
    }
  }

  if (!nPads)
  {
    AliError("Got no pad from the iterator ?! That's not normal. Please check !");
    return;
  }

  TString msg;
  Int_t ntimes;

  cout << Form("According to mask %x (human readable form below) %6d pads are bad (over a total of %6d, i.e. %7.2f %%)",
               mask,nBadPads,nPads,nBadPads*100.0/nPads) << endl;
  cout << AliMUONPadStatusMaker::AsCondition(mask) << endl;
  cout << "--------" << endl;

  while ( log.Next(msg,ntimes) )
  {
    cout << Form("The message (%120s) occured %15d times (%7.4f %%)",msg.Data(),ntimes,ntimes*100.0/nPads) << endl;
  }

  TMap* hvMap = CalibrationData().HV();

  std::cout << "Map UniqueID = " << hvMap->GetUniqueID() << std::endl;
}

//_____________________________________________________________________________
AliMUONVTrackerData* AliMUONPadStatusMaker::InternalTrackerData() const
{
    if (!fTrackerData && fkCalibrationData.OccupancyMap()) 
    {
        fTrackerData = new AliMUONTrackerData("OCC","OCC",*(fkCalibrationData.OccupancyMap()));
    }
    return fTrackerData;
}

//_____________________________________________________________________________
TExMap* AliMUONPadStatusMaker::InternalHV() const
{
    if (!fHV && fkCalibrationData.HV())
    {
        fHV = new TExMap;
    }
    return fHV;
}

//_____________________________________________________________________________
TExMap* AliMUONPadStatusMaker::InternalLV() const
{
    if (!fLV && fkCalibrationData.LV())
    {
        fLV = new TExMap;
    }
    return fLV;
}
