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

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Class providing the calibration parameters by accessing the CDB           //
//                                                                           //
// Request an instance with AliTRDcalibDB::Instance()                        //
// If a new event is processed set the event number with SetRun              //
// Then request the calibration data                                         // 
//                                                                           //
// Author:                                                                   //
//   Jan Fiete Grosse-Oetringhaus (Jan.Fiete.Grosse-Oetringhaus@cern.ch)     //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <TROOT.h>
#include <TClonesArray.h>
#include <TObjArray.h>

#include "AliCDBManager.h"
#include "AliCDBEntry.h"
#include "AliLog.h"

#include "AliTRDPIDReference.h"
#include "AliTRDPIDResponseObject.h"
#include "AliTRDcalibDB.h"
#include "AliTRDtrapConfig.h"
#include "AliTRDtrapConfigHandler.h"
#include "AliTRDCommonParam.h"
#include "AliTRDgeometry.h"

#include "AliTRDCalROC.h"
#include "AliTRDCalPad.h"
#include "AliTRDCalDet.h"
#include "AliTRDCalDCS.h"
#include "AliTRDCalDCSFEE.h"
#include "AliTRDCalDCSv2.h"
#include "AliTRDCalDCSFEEv2.h"
#include "AliTRDCalPID.h"
#include "AliTRDCalMonitoring.h"
#include "AliTRDCalChamberStatus.h"
#include "AliTRDCalPadStatus.h"
#include "AliTRDCalSingleChamberStatus.h"
#include "AliTRDCalTrkAttach.h"
#include "AliTRDCalOnlineGainTable.h"

ClassImp(AliTRDcalibDB)

AliTRDcalibDB *AliTRDcalibDB::fgInstance   = 0;
Bool_t         AliTRDcalibDB::fgTerminated = kFALSE;

//_ singleton implementation __________________________________________________
AliTRDcalibDB* AliTRDcalibDB::Instance()
{
  //
  // Singleton implementation
  // Returns an instance of this class, it is created if neccessary
  //

  if (fgTerminated != kFALSE) {
    return 0;
  }

  if (fgInstance == 0) {
    fgInstance = new AliTRDcalibDB();
  }

  return fgInstance;

}

//_____________________________________________________________________________
void AliTRDcalibDB::Terminate()
{
  //
  // Singleton implementation
  // Deletes the instance of this class and sets the terminated flag,
  // instances cannot be requested anymore
  // This function can be called several times.
  //
  
  fgTerminated = kTRUE;
  
  if (fgInstance != 0) {
    delete fgInstance;
    fgInstance = 0;
  }

}

//_____________________________________________________________________________
AliTRDcalibDB::AliTRDcalibDB()
  :TObject()
  ,fRun(-1)
  ,fPRFsmp(0)
  ,fPRFbin(0)
  ,fPRFlo(0)
  ,fPRFhi(0)
  ,fPRFwid(0)
  ,fPRFpad(0)
  ,fPIDResponse(NULL)
  ,fOnlineGainTableID(0)
  ,fTrapConfig(0x0)
  ,fTrapConfigName("")
  ,fTrapConfigVersion("")
{
  //
  // Default constructor
  //
  // TODO Default runnumber is set to 0, this should be changed later
  //      to an invalid value (e.g. -1) to prevent
  // TODO invalid calibration data to be used.
  //

  fRun = AliCDBManager::Instance()->GetRun();

  for (Int_t i = 0; i < kCDBCacheSize; ++i) {
    fCDBCache[i]   = 0;
    fCDBEntries[i] = 0;
  }
  
  // Create the sampled PRF
  SamplePRF();
  
}

//_____________________________________________________________________________
AliTRDcalibDB::AliTRDcalibDB(const AliTRDcalibDB &c)
  :TObject(c)
  ,fRun(-1)
  ,fPRFsmp(0)
  ,fPRFbin(0)
  ,fPRFlo(0)
  ,fPRFhi(0)
  ,fPRFwid(0)
  ,fPRFpad(0)
  ,fPIDResponse(NULL)
  ,fOnlineGainTableID(0)
  ,fTrapConfig(0x0)
  ,fTrapConfigName("")
  ,fTrapConfigVersion("")
{
  //
  // Copy constructor (not that it make any sense for a singleton...)
  //

  for (Int_t i = 0; i < kCDBCacheSize; ++i) {
    fCDBCache[i]   = 0;
    fCDBEntries[i] = 0;
  }
  
  // Create the sampled PRF
  SamplePRF();

}

//_____________________________________________________________________________
AliTRDcalibDB &AliTRDcalibDB::operator=(const AliTRDcalibDB &c) 
{
  //
  // Assignment operator (same as above ...)
  //

  if (this != &c) {
    AliFatal("No assignment operator defined");
  }

  return *this;

}

//_____________________________________________________________________________
AliTRDcalibDB::~AliTRDcalibDB() 
{
  //
  // destructor
  //
  
  if (fPRFsmp) {
    delete [] fPRFsmp;
    fPRFsmp = 0;
  }

  if (fPIDResponse) {
    delete fPIDResponse;
    fPIDResponse = 0x0;
  }

  Invalidate();
  fgInstance   = 0;
}

//_caching functions____________________________________________________________
const TObject *AliTRDcalibDB::GetCachedCDBObject(Int_t id)
{
  //
  // Retrieves a cdb object with the given id. The objects are cached as
  // long as the run number is not changed.
  //
  // Put together the available objects here by using the lines
  //   a) For usual calibration objects:
  //      case kID<Name> : 
  //        return CacheCDBEntry(kID<Name>,"TRD/Calib/<Path>"); 
  //        break;
  //      See function CacheCDBEntry for details.
  //   and
  //   b) For calibration data which depends on two objects: One containing 
  //      a value per detector and one the local fluctuations per pad:
  //      case kID<Name> :
  //        return CacheMergeCDBEntry(kID<Name>,"TRD/Calib/<padPath>","TRD/Calib/<chamberPath>"); 
  //        break;
  //      See function CacheMergeCDBEntry for details.
  //
    
  switch (id) {

    // Parameters defined per pad and chamber
    case kIDVdriftPad : 
      return CacheCDBEntry(kIDVdriftPad         ,"TRD/Calib/LocalVdrift"); 
      break;
    case kIDVdriftChamber : 
      return CacheCDBEntry(kIDVdriftChamber     ,"TRD/Calib/ChamberVdrift"); 
      break;
    case kIDExBChamber : 
      return CacheCDBEntry(kIDExBChamber        ,"TRD/Calib/ChamberExB"); 
      break;
    case kIDT0Pad : 
      return CacheCDBEntry(kIDT0Pad             ,"TRD/Calib/LocalT0"); 
      break;
    case kIDT0Chamber : 
      return CacheCDBEntry(kIDT0Chamber         ,"TRD/Calib/ChamberT0"); 
      break;
    case kIDGainFactorPad : 
      return CacheCDBEntry(kIDGainFactorPad     ,"TRD/Calib/LocalGainFactor"); 
      break;
    case kIDGainFactorChamber : 
      return CacheCDBEntry(kIDGainFactorChamber ,"TRD/Calib/ChamberGainFactor"); 
      break;

    case kIDOnlineGainFactor : 
      switch(GetOnlineGainTableID()) {
        case 0:
          // For testing purposes only !!!
          AliInfo("No gain table name from OCDB. Use default table!");
          return CacheCDBEntry(kIDOnlineGainFactor  ,"TRD/Calib/Krypton_2011-01"); 
          break;
        case 1:
	  // Online gain table ID 1
          return CacheCDBEntry(kIDOnlineGainFactor  ,"TRD/Calib/Krypton_2011-01"); 
          break;
        case 2:
	  // Online gain table ID 2
          return CacheCDBEntry(kIDOnlineGainFactor  ,"TRD/Calib/Gaintbl_Uniform_FGAN0_2011-01"); 
          break;
        case 3:
	  // Online gain table ID 3
          return CacheCDBEntry(kIDOnlineGainFactor  ,"TRD/Calib/Gaintbl_Uniform_FGAN8_2011-01"); 
          break;
        case 4:
	  // Online gain table ID 4
          return CacheCDBEntry(kIDOnlineGainFactor  ,"TRD/Calib/Krypton_2011-02"); 
          break;
        case 5:
	  // Online gain table ID 5
          return CacheCDBEntry(kIDOnlineGainFactor  ,"TRD/Calib/Krypton_2011-03"); 
          break;
        case 6:
	  // Online gain table ID 6
          return CacheCDBEntry(kIDOnlineGainFactor  ,"TRD/Calib/Gaintbl_Uniform_FGAN0_2012-01"); 
          break;
        case 7:
	  // Online gain table ID 7
          return CacheCDBEntry(kIDOnlineGainFactor  ,"TRD/Calib/Gaintbl_Uniform_FGAN8_2012-01");
          break; 
        case 8:
	  // Online gain table ID 8
          return CacheCDBEntry(kIDOnlineGainFactor  ,"TRD/Calib/Krypton_2012-01"); 
          break;
        case 9:
	  // Online gain table ID 9
          return CacheCDBEntry(kIDOnlineGainFactor  ,"TRD/Calib/Krypton_2015-01"); 
          break;
        case 10:
	  // Online gain table ID 10
          return CacheCDBEntry(kIDOnlineGainFactor  ,"TRD/Calib/Krypton_2015-02"); 
          break;
      default:
	AliError(Form("unknown gaintable requested with ID"));
      }
      break;

    case kIDNoiseChamber : 
      return CacheCDBEntry(kIDNoiseChamber      ,"TRD/Calib/DetNoise"); 
      break;
    case kIDNoisePad : 
      return CacheCDBEntry(kIDNoisePad          ,"TRD/Calib/PadNoise"); 
      break;

    // Parameters defined per pad
    case kIDPRFWidth : 
      return CacheCDBEntry(kIDPRFWidth          ,"TRD/Calib/PRFWidth"); 
      break;

    // Status values
    case kIDChamberStatus : 
      return CacheCDBEntry(kIDChamberStatus     ,"TRD/Calib/ChamberStatus"); 
      break;
    case kIDPadStatus : 
      return CacheCDBEntry(kIDPadStatus         ,"TRD/Calib/PadStatus"); 
      break;

    // Global parameters
    case kIDMonitoringData : 
      return CacheCDBEntry(kIDMonitoringData    ,"TRD/Calib/MonitoringData"); 
      break;
    case kIDFEE : 
      return CacheCDBEntry(kIDFEE               ,"TRD/Calib/FEE"); 
      break;
    case kIDTrapConfig :
      return CacheCDBEntry(kIDTrapConfig        ,"TRD/Calib/TrapConfig"); 
      break;
    case kIDDCS :
      return CacheCDBEntry(kIDDCS               ,"TRD/Calib/DCS");
      break;
    case kIDPIDNN : 
      return CacheCDBEntry(kIDPIDNN             ,"TRD/Calib/PIDNN");
      break;
    case kIDPIDLQ : 
      return CacheCDBEntry(kIDPIDLQ             ,"TRD/Calib/PIDLQ"); 
      break;
    case kIDPIDLQ1D:
      return CacheCDBEntry(kIDPIDLQ1D           ,"TRD/Calib/PIDLQ1D");
      break;
    case kIDRecoParam : 
      return CacheCDBEntry(kIDRecoParam         ,"TRD/Calib/RecoParam"); 
      break;
    case kIDAttach : 
      return CacheCDBEntry(kIDAttach            ,"TRD/Calib/TrkAttach"); 
      break;
    case kIDPHQ :
      return CacheCDBEntry(kIDPHQ               ,"TRD/Calib/PHQ");
      break;
  }

  return 0;

}

//_____________________________________________________________________________
AliCDBEntry *AliTRDcalibDB::GetCDBEntry(const char *cdbPath)
{
  // 
  // Retrieves an entry with path <cdbPath> from the CDB.
  //
    
  AliCDBEntry *entry = AliCDBManager::Instance()->Get(cdbPath,fRun);
  if (!entry) { 
    AliError(Form("Failed to get entry: %s",cdbPath));
    return 0; 
  }
  
  return entry;

}

//_____________________________________________________________________________
const TObject *AliTRDcalibDB::CacheCDBEntry(Int_t id, const char *cdbPath)
{
  //
  // Caches the entry <id> with cdb path <cdbPath>
  //

  if (!fCDBCache[id]) {
    fCDBEntries[id] = GetCDBEntry(cdbPath);
    if (fCDBEntries[id]) {
      fCDBCache[id] = fCDBEntries[id]->GetObject();
      if (id == kIDOnlineGainFactor)
	AliInfo(Form("loaded gain table: %s", fCDBEntries[id]->GetId().GetAliCDBPath().GetPath().Data()));
    }
  } 
  
  return fCDBCache[id];

}

//_____________________________________________________________________________
void AliTRDcalibDB::SetRun(Long64_t run)
{
  //
  // Sets current run number. Calibration data is read from the corresponding file.
  // When the run number changes the caching is invalidated.
  //

  if (fRun == run) {
    return;
  }

  fRun = run;

  Invalidate();

}

//_____________________________________________________________________________
void AliTRDcalibDB::Invalidate()
{
  //
  // Invalidates cache (when run number is changed).
  //
  
  for (Int_t i = 0; i < kCDBCacheSize; ++i) {
    if (fCDBEntries[i]) {
      if (AliCDBManager::Instance()->GetCacheFlag() == kFALSE) {
        if ((fCDBEntries[i]->IsOwner() == kFALSE) && 
            (fCDBCache[i])) {
          delete fCDBCache[i];
	}
        delete fCDBEntries[i];
      }
      fCDBEntries[i] = 0;
      fCDBCache[i]   = 0;
    }
  }

  fOnlineGainTableID = 0;

}

//_____________________________________________________________________________
Float_t AliTRDcalibDB::GetNoise(Int_t det, Int_t col, Int_t row)
{
  //
  // Returns the noise level in ADC counts for the given pad.
  //

  const AliTRDCalPad *calPad     = dynamic_cast<const AliTRDCalPad *> 
                                   (GetCachedCDBObject(kIDNoisePad));
  if (!calPad) {
    return -1;
  }

  AliTRDCalROC       *roc        = calPad->GetCalROC(det);
  if (!roc) {
    return -1;
  }

  const AliTRDCalDet *calChamber = dynamic_cast<const AliTRDCalDet *> 
                                   (GetCachedCDBObject(kIDNoiseChamber));
  if (!calChamber) {
    return -1;
  }

  return calChamber->GetValue(det) * roc->GetValue(col,row);

}
 
//_____________________________________________________________________________
AliTRDCalROC *AliTRDcalibDB::GetNoiseROC(Int_t det)
{
  //
  // Returns the Vdrift calibration object for a given ROC
  // containing one number per pad 
  //
  
  const AliTRDCalPad     *calPad     = dynamic_cast<const AliTRDCalPad *> 
                                       (GetCachedCDBObject(kIDNoisePad));
  if (!calPad) {
    return 0;
  }

  AliTRDCalROC           *roc        = calPad->GetCalROC(det);
  if (!roc) {
    return 0;
  }
  else {
    return roc;
  }

}

//_____________________________________________________________________________
const AliTRDCalDet *AliTRDcalibDB::GetNoiseDet()
{
  //
  // Returns the Vdrift calibration object
  // containing one number per detector
  //
  
  const AliTRDCalDet     *calChamber = dynamic_cast<const AliTRDCalDet *> 
                                       (GetCachedCDBObject(kIDNoiseChamber));
  if (!calChamber) {
    return 0;
  }
  else {
    return calChamber;
  }

}

//_____________________________________________________________________________
Float_t AliTRDcalibDB::GetVdrift(Int_t det, Int_t col, Int_t row)
{
  //
  // Returns the drift velocity for the given pad.
  //

  const AliTRDCalPad *calPad     = dynamic_cast<const AliTRDCalPad *> 
                                   (GetCachedCDBObject(kIDVdriftPad));
  if (!calPad) {
    return -1;
  }

  AliTRDCalROC       *roc        = calPad->GetCalROC(det);
  if (!roc) {
    return -1;
  }

  const AliTRDCalDet *calChamber = dynamic_cast<const AliTRDCalDet *> 
                                   (GetCachedCDBObject(kIDVdriftChamber));
  if (!calChamber) {
    return -1;
  }

  return calChamber->GetValue(det) * roc->GetValue(col,row);

}
 
//_____________________________________________________________________________
AliTRDCalROC *AliTRDcalibDB::GetVdriftROC(Int_t det)
{
  //
  // Returns the Vdrift calibration object for a given ROC
  // containing one number per pad 
  //
  
  const AliTRDCalPad     *calPad     = dynamic_cast<const AliTRDCalPad *> 
                                       (GetCachedCDBObject(kIDVdriftPad));
  if (!calPad) {
    return 0;
  }

  AliTRDCalROC           *roc        = calPad->GetCalROC(det);
  if (!roc) {
    return 0;
  }
  else {
    return roc;
  }

}

//_____________________________________________________________________________
const AliTRDCalDet *AliTRDcalibDB::GetVdriftDet()
{
  //
  // Returns the Vdrift calibration object
  // containing one number per detector
  //
  
  const AliTRDCalDet     *calChamber = dynamic_cast<const AliTRDCalDet *> 
                                       (GetCachedCDBObject(kIDVdriftChamber));
  if (!calChamber) {
    return 0;
  }
  else {
    return calChamber;
  }

}

//_____________________________________________________________________________
TObjArray * AliTRDcalibDB::GetPHQ()
{
  //
  //return PHQ calibration object
  //
  TObjArray *arr = (TObjArray *) (GetCachedCDBObject(kIDPHQ));
  return arr;
}

//_____________________________________________________________________________
Float_t AliTRDcalibDB::GetVdriftAverage(Int_t det)
{
  //
  // Returns the average drift velocity for the given detector
  //

  const AliTRDCalDet *calDet     = dynamic_cast<const AliTRDCalDet *> 
                                   (GetCachedCDBObject(kIDVdriftChamber));
  if (!calDet) {
    return -1;
  }

  return calDet->GetValue(det);

}
//_____________________________________________________________________________
const AliTRDCalDet *AliTRDcalibDB::GetExBDet()
{
  //
  // Returns the exB calibration object
  // containing one number per detector
  //
  
  const AliTRDCalDet     *calChamber = dynamic_cast<const AliTRDCalDet *> (GetCachedCDBObject(kIDExBChamber));
  
  Double_t meanexb = 100.0;
  if (calChamber) meanexb = calChamber->GetMean();
  //printf("mean %f\n",meanexb);  

  if ((!calChamber) || (meanexb > 70.0)) {
    
    const AliTRDCalDet     *calChambervdrift = dynamic_cast<const AliTRDCalDet *> 
	(GetCachedCDBObject(kIDVdriftChamber));
      if (!calChambervdrift) {
	return 0;
      }
      else {
	AliTRDCalDet *calDetExB = new AliTRDCalDet("lorentz angle tan","lorentz angle tan (detector value)");
	for(Int_t k = 0; k < 540; k++){
	  calDetExB->SetValue(k,AliTRDCommonParam::Instance()->GetOmegaTau(calChambervdrift->GetValue(k)));
	}
	return calDetExB;
      }
  }
  else return calChamber;

}
//_____________________________________________________________________________
Float_t AliTRDcalibDB::GetT0(Int_t det, Int_t col, Int_t row)
{
  //
  // Returns t0 for the given pad.
  //
  
  const AliTRDCalPad     *calPad     = dynamic_cast<const AliTRDCalPad *> 
                                       (GetCachedCDBObject(kIDT0Pad));
  if (!calPad) {
    return -1;
  }

  AliTRDCalROC           *roc        = calPad->GetCalROC(det);
  if (!roc) {
    return -1;
  }

  const AliTRDCalDet     *calChamber = dynamic_cast<const AliTRDCalDet *> 
                                       (GetCachedCDBObject(kIDT0Chamber));
  if (!calChamber) {
    return -1;
  }

  return calChamber->GetValue(det) + roc->GetValue(col,row);

}
 
//_____________________________________________________________________________
AliTRDCalROC *AliTRDcalibDB::GetT0ROC(Int_t det)
{
  //
  // Returns the t0 calibration object for a given ROC
  // containing one number per pad 
  //
  
  const AliTRDCalPad     *calPad     = dynamic_cast<const AliTRDCalPad *> 
                                       (GetCachedCDBObject(kIDT0Pad));
  if (!calPad) {
    return 0;
  }

  AliTRDCalROC           *roc        = calPad->GetCalROC(det);
  if (!roc) {
    return 0;
  }
  else {
    return roc;
  }

}

//_____________________________________________________________________________
const AliTRDCalDet *AliTRDcalibDB::GetT0Det()
{
  //
  // Returns the t0 calibration object
  // containing one number per detector
  //
  
  const AliTRDCalDet     *calChamber = dynamic_cast<const AliTRDCalDet *> 
                                       (GetCachedCDBObject(kIDT0Chamber));
  if (!calChamber) {
    return 0;
  }
  else {
    return calChamber;
  }

}

//_____________________________________________________________________________
Float_t AliTRDcalibDB::GetT0Average(Int_t det)
{
  //
  // Returns the average t0 for the given detector
  //

  const AliTRDCalPad *calPad     = dynamic_cast<const AliTRDCalPad *> 
                                   (GetCachedCDBObject(kIDT0Pad));
  if (!calPad) {
    return -1;
  }

  AliTRDCalROC       *roc        = calPad->GetCalROC(det);
  if (!roc) {
    return -1;
  }

  const AliTRDCalDet *calDet     = dynamic_cast<const AliTRDCalDet *> 
                                   (GetCachedCDBObject(kIDT0Chamber));
  if (!calDet) {
    return -1;
  }

  Double_t sum = 0.0; 
  for (Int_t channel = 0; channel < roc->GetNchannels(); ++channel) {
    sum += roc->GetValue(channel);
  }
  sum /= roc->GetNchannels();
  sum += calDet->GetValue(det);
  return sum;

}

//_____________________________________________________________________________
AliTRDCalOnlineGainTableROC* AliTRDcalibDB::GetOnlineGainTableROC(Int_t det)
{
  //
  // Returns the online gain factor table for a given ROC.
  //

  if (!HasOnlineFilterGain()) {
    return 0x0;
  }
  
  const AliTRDCalOnlineGainTable *calOnline 
     = dynamic_cast<const AliTRDCalOnlineGainTable *> 
                                   (GetCachedCDBObject(kIDOnlineGainFactor));
  if (!calOnline) {
    return 0x0;
  }

  return calOnline->GetGainTableROC(det);

}

//_____________________________________________________________________________
Float_t AliTRDcalibDB::GetOnlineGainFactor(Int_t det, Int_t col, Int_t row)
{
  //
  // Returns the online gain factor for the given pad.
  //

  if (!HasOnlineFilterGain()) {
    return -1;
  }
  
  const AliTRDCalOnlineGainTable *calOnline 
     = dynamic_cast<const AliTRDCalOnlineGainTable *> 
                                   (GetCachedCDBObject(kIDOnlineGainFactor));
  if (!calOnline) {
    return -1;
  }

  return calOnline->GetGainCorrectionFactor(det,row,col);

}

//_____________________________________________________________________________
AliTRDCalROC *AliTRDcalibDB::GetGainFactorROC(Int_t det)
{
  //
  // Returns the gain factor calibration object for a given ROC
  //
  
  const AliTRDCalPad *calPad     = dynamic_cast<const AliTRDCalPad *> 
                                   (GetCachedCDBObject(kIDGainFactorPad));
  if (!calPad) {
    return 0;
  }

  AliTRDCalROC       *roc        = calPad->GetCalROC(det);
  if (!roc) {
    return 0;
  }
  else {
    return roc;
  }

}

//_____________________________________________________________________________
Float_t AliTRDcalibDB::GetGainFactor(Int_t det, Int_t col, Int_t row)
{
  //
  // Returns the gain factor for the given pad.
  //
  
  const AliTRDCalPad *calPad     = dynamic_cast<const AliTRDCalPad *> 
                                   (GetCachedCDBObject(kIDGainFactorPad));
  if (!calPad) {
    return -1;
  }

  AliTRDCalROC       *roc        = calPad->GetCalROC(det);
  if (!roc) {
    return -1;
  }

  const AliTRDCalDet *calChamber = dynamic_cast<const AliTRDCalDet *> 
                                   (GetCachedCDBObject(kIDGainFactorChamber));
  if (!calChamber) {
    return -1;
  }

  return calChamber->GetValue(det) * roc->GetValue(col,row);

}

//_____________________________________________________________________________
const AliTRDCalDet *AliTRDcalibDB::GetGainFactorDet()
{
  //
  // Returns the gain factor calibration object
  // containing one number per detector
  //

  const AliTRDCalDet *calChamber = dynamic_cast<const AliTRDCalDet *> 
                                   (GetCachedCDBObject(kIDGainFactorChamber));
  if (!calChamber) {
    return 0;
  }
  else {
    return calChamber;
  }

}

//_____________________________________________________________________________
Float_t AliTRDcalibDB::GetGainFactorAverage(Int_t det)
{
  //
  // Returns the average gain factor for the given detector
  //

  const AliTRDCalDet *calDet     = dynamic_cast<const AliTRDCalDet *> 
                                   (GetCachedCDBObject(kIDGainFactorChamber));
  if (!calDet) {
    return -1;
  }

  return calDet->GetValue(det);

}

//_____________________________________________________________________________
AliTRDCalROC *AliTRDcalibDB::GetPRFROC(Int_t det)
{
  //
  // Returns the PRF calibration object for a given ROC
  // containing one number per pad 
  //
  
  const AliTRDCalPad     *calPad     = dynamic_cast<const AliTRDCalPad *> 
                                       (GetCachedCDBObject(kIDPRFWidth));
  if (!calPad) {
    return 0;
  }

  AliTRDCalROC           *roc        = calPad->GetCalROC(det);
  if (!roc) {
    return 0;
  }
  else {
    return roc;
  }

}

//_____________________________________________________________________________
Float_t AliTRDcalibDB::GetPRFWidth(Int_t det, Int_t col, Int_t row)
{
  //
  // Returns the PRF width for the given pad.
  //
  
  const AliTRDCalPad *calPad     = dynamic_cast<const AliTRDCalPad *> 
                                   (GetCachedCDBObject(kIDPRFWidth));
  if (!calPad) {
    return -1;
  }

  AliTRDCalROC       *roc        = calPad->GetCalROC(det);
  if (!roc) {
    return -1;
  }

  return roc->GetValue(col,row);

}
  
//_____________________________________________________________________________
Int_t AliTRDcalibDB::ExtractTimeBinsFromString(TString tbstr)
{
  // Check if there is any content in the string first
  if (tbstr.Length() == 0) {
    AliErrorClass("Parameter for number of timebins is empty!");
    return -1;
  }

  // Check if we have the correct config parameter
  TString tbident  = "tb";
  TString tbsubstr = tbstr(0,2);
  if (!tbsubstr.EqualTo(tbident)) {
    AliErrorClass(Form("Parameter for number of timebins is corrupted (%s)!", tbstr.Data()));
    return -1;
  }

  tbstr.Remove(0,2);
  // check if there is more than a number
  if (!tbstr.IsDigit()) {
    AliErrorClass(Form("Parameter for number of timebins is corrupted (%s)!", tbstr.Data()));
    return -1;
  }

  return tbstr.Atoi();

}

//_____________________________________________________________________________
Int_t AliTRDcalibDB::GetNumberOfTimeBinsDCSBoard(){
  // This is an old way to extract the number of time bins,
  // only to be used as a fallback to see whether there's
  // a patched OCDB object!
  
  Int_t nTB=-1;
  // Get the array with SOR and EOR
  const TObjArray *dcsArr = dynamic_cast<const TObjArray *>(GetCachedCDBObject(kIDDCS));
  if(!dcsArr)
    return -1; // Error
  // Check CalDCS version
  Int_t calver = 0;
  if (!strcmp(dcsArr->At(0)->ClassName(),"AliTRDCalDCS")) calver = 1;
  else if (!strcmp(dcsArr->At(0)->ClassName(),"AliTRDCalDCSv2")) calver = 2;
  // CalDCS version 1
  if (calver == 1) {
    // DCS object
    AliTRDCalDCS *calSOR = dynamic_cast<AliTRDCalDCS *>(dcsArr->At(0));
    // SOR mandantory
    if(!calSOR)
      return -1; // Error
    else
      nTB=calSOR->GetGlobalNumberOfTimeBins();
    AliTRDCalDCS *calEOR = dynamic_cast<AliTRDCalDCS *>(dcsArr->At(1));
    if(calEOR && calEOR->GetGlobalNumberOfTimeBins()!=nTB)
      return -2; // Mixed
  }
  else if (calver == 2) {
    // DCS object
    AliTRDCalDCSv2 *calSOR = dynamic_cast<AliTRDCalDCSv2 *>(dcsArr->At(0));
    // SOR mandantory
    if(!calSOR)
      return -1; // Error
    else
      nTB=calSOR->GetGlobalNumberOfTimeBins();
    AliTRDCalDCSv2 *calEOR = dynamic_cast<AliTRDCalDCSv2 *>(dcsArr->At(1));
    if(calEOR && calEOR->GetGlobalNumberOfTimeBins()!=nTB)
      return -2; // Mixed
  }
  else{
    // No version 1 nor version 2 object
    return -1; // Error
  }

  // All well
  return nTB;
}
//_____________________________________________________________________________
Int_t AliTRDcalibDB::GetNumberOfTimeBinsDCS()
{
  //
  // Returns number of time bins from the DCS
  //

  // Initialize with values indicating no information
  TString cfg="";
  Int_t nTB = -1;
  // Extract the global configuration name
  GetGlobalConfigurationByChamber(cfg,kTimebin);
  // Extract the number of time bins from
  // the global configuration 
  if(cfg.Length()>0){
    nTB = ExtractTimeBinsFromString(cfg);
  }
  if(nTB>0){
    // All well, we could extract the number
    // of time bins from the configuration name
    return nTB;
  }
  else{
    // No number of time bins from config name.
    // No board responded or similar.
    // We patched some OCDB entries with 
    // only the global number of time bins set.
    // We should get these here
    nTB=GetNumberOfTimeBinsDCSBoard();
    if(nTB>0){
      AliWarning("Using old method for number of time bins."
		 " This is probably a patched OCDB entry");
      return nTB;
    }
    else{
      // Error
      AliError("No number of time bins either"
	       " from config name or patched OCDB entry");
      return -1;
    }
  }

  
  // // Get the corresponding parameter
  // TString cfgstr = "", cfgname = "";
  // GetGlobalConfiguration(cfgname);
  // if(cfgname.Length()==0)
  //   return -1;
  // GetDCSConfigParOption(cfgname, kTimebin, 0, cfgstr);
  // if(cfgstr.Length()==0)
  //   return -1;
  // return ExtractTimeBinsFromString(cfgstr);

}

//_____________________________________________________________________________
void AliTRDcalibDB::GetFilterType(TString &filterType)
{
  //
  // Returns the filter type
  //

  TString cfgname = "";
  GetGlobalConfiguration(cfgname);
  GetDCSConfigParOption(cfgname, kFltrSet, 0, filterType);

}

//_____________________________________________________________________________
Int_t AliTRDcalibDB::GetOnlineGainTableID()
{
  //
  // Get the gain table ID from the DCS
  //

  if (fOnlineGainTableID > 0) {
    return fOnlineGainTableID;
  }

  const TObjArray *dcsArr = dynamic_cast<const TObjArray *>(GetCachedCDBObject(kIDDCS));
  if (!dcsArr){
    return -1;
  }

  Int_t esor   = 0; // Take SOR
  Int_t calver = 0; // Check CalDCS version
  if (!strcmp(dcsArr->At(0)->ClassName(),"AliTRDCalDCS"))   calver = 1;
  if (!strcmp(dcsArr->At(0)->ClassName(),"AliTRDCalDCSv2")) calver = 2;

  if      (calver == 1) {

    // No data for old DCS object available, anyway
    return -1;

  } 
  else if (calver == 2) {

    // DCSv2 object
    const AliTRDCalDCSv2 *calDCSv2 = dynamic_cast<const AliTRDCalDCSv2 *>(dcsArr->At(esor));
    if(!calDCSv2){
      return -1;
    }

    TString tableName = "";
    for (Int_t i = 0; i < 540; i++) {
      const AliTRDCalDCSFEEv2 *calDCSFEEv2 = calDCSv2->GetCalDCSFEEObj(i);
      if (!calDCSFEEv2) {
        continue;
      }
      const TString tableNameTmp = calDCSFEEv2->GetGainTableName();
      if (tableNameTmp.Length() > 0) {
        if ((tableName.Length() > 0) &&
            (tableName != tableNameTmp)) {
          AliFatal(Form("Inconsistent gain table names! %s - %s"
                       ,tableName.Data(),tableNameTmp.Data()));
          continue; // maybe return -1;
        }
        tableName = tableNameTmp; // this contains the first entry found
      }
    }

    if      (tableName.CompareTo("Krypton_2011-01")               == 0)
      fOnlineGainTableID = 1;
    else if (tableName.CompareTo("Gaintbl_Uniform_FGAN0_2011-01") == 0)
      fOnlineGainTableID = 2;
    else if (tableName.CompareTo("Gaintbl_Uniform_FGAN8_2011-01") == 0)
      fOnlineGainTableID = 3;
    else if (tableName.CompareTo("Krypton_2011-02")               == 0)
      fOnlineGainTableID = 4;
    else if (tableName.CompareTo("Krypton_2011-03")               == 0)
      fOnlineGainTableID = 5;
    else if (tableName.CompareTo("Gaintbl_Uniform_FGAN0_2012-01") == 0)
      fOnlineGainTableID = 6;
    else if (tableName.CompareTo("Gaintbl_Uniform_FGAN8_2012-01") == 0)
      fOnlineGainTableID = 7;
    else if (tableName.CompareTo("Krypton_2012-01")               == 0)
      fOnlineGainTableID = 8;
    else if (tableName.CompareTo("Krypton_2015-01")               == 0)
      fOnlineGainTableID = 9;
    else if (tableName.CompareTo("Krypton_2015-02")               == 0)
      fOnlineGainTableID = 10;
    else
      AliFatal(Form("unknown gaintable <%s> requested", tableName.Data()));

    AliInfo(Form("looking for gaintable: %s (id %i)",
		 tableName.Data(), fOnlineGainTableID));

    if (fOnlineGainTableID > 0)
      return fOnlineGainTableID;
  } 
  else {

    AliError("NO DCS/DCSv2 OCDB entry found!");
    return -1;

  }

  return -1;

}

//_____________________________________________________________________________
void AliTRDcalibDB::GetGlobalConfiguration(TString &config)
{
  //
  // Get Configuration from the DCS
  //

  const TObjArray *dcsArr = dynamic_cast<const TObjArray *>(GetCachedCDBObject(kIDDCS));
  if(!dcsArr){
    AliError("No DCS CDB Object available!");
    config = "";
    return;
  }

  Int_t idSOR = 0, idEOR=1; // The index of SOR and EOR
  Bool_t hasSOR = (dcsArr->At(idSOR));
  Bool_t hasEOR = (dcsArr->At(idEOR));
  TString cfgSOR = "", cfgEOR = ""; // The configuration at SOR/EOR

  // The SOR object is mandatory
  if (!hasSOR) {
    AliError("NO SOR object found in CDB file!");
    config = "";
    return;
  }
  if (!hasEOR) AliWarning("NO EOR object found in CDB file!");

  // Check CalDCS version
  Int_t calver = 0;
  if (!strcmp(dcsArr->At(idSOR)->ClassName(),"AliTRDCalDCS")) calver = 1;
  else if (!strcmp(dcsArr->At(idSOR)->ClassName(),"AliTRDCalDCSv2")) calver = 2;

  // Get the configuration strings
  if (calver == 1) {
    // DCS object
    const AliTRDCalDCS *calSOR = dynamic_cast<const AliTRDCalDCS *>(dcsArr->At(idSOR));
    cfgSOR = calSOR->GetGlobalConfigName();
    if (hasEOR) {
      const AliTRDCalDCS *calEOR = dynamic_cast<const AliTRDCalDCS *>(dcsArr->At(idEOR));
      cfgEOR = calEOR->GetGlobalConfigName();
    }
  } 
  else if (calver == 2) {
    // DCSv2 object
    const AliTRDCalDCSv2 *calv2SOR = dynamic_cast<const AliTRDCalDCSv2 *>(dcsArr->At(idSOR));
    cfgSOR = calv2SOR->GetGlobalConfigName();
    if (hasEOR) {
      const AliTRDCalDCSv2 *calv2EOR = dynamic_cast<const AliTRDCalDCSv2 *>(dcsArr->At(idEOR));
      cfgEOR = calv2EOR->GetGlobalConfigName();
    }
  } 
  else {
    AliError("NO DCS/DCSv2 OCDB entry found!");
    config = "";
    return;
  }

  // If there is no EOR entry, return the SOR value
  if (!hasEOR || cfgEOR.Length()==0) {
    config = cfgSOR;
    return;
  }

  // Check if the configuration is the same for both
  if (cfgSOR.EqualTo(cfgEOR)) {
    config = cfgSOR;
    return;
  }

  // When both SOR and EOR have an entry but are different, the config is not defined
  AliError("Inconsistent configuration at start and end of run found!");
  config = "";
  return;

}
//_____________________________________________________________________________
void AliTRDcalibDB::GetGlobalConfigurationByChamber(TString &config,Int_t par, Int_t opt)
{
  //
  // Get Configuration from the DCS
  //

  // par is the enumeration from kFltrSet = 1 ... kAddOpti 6
  //   example:
  //      cf_p_nozs_tb30_csmtrk_ptrg
  //   par   1   2    3    4     5
  //      kFltrSet  kTimebin   kTrigSet
  //          kReadout  kTrkMode
  // opt is 0 for the value itself and 1,2,3 for the first, second option
  // opt has standard value of 0 in header file


  // The point is that we extract the parameter/option for 
  // each of the 540 chambers and compare only this 
  // parameter/option for a global value, not the full
  // configuration name with all parameters

  // Get the array with SOR and EOR
  const TObjArray *dcsArr = dynamic_cast<const TObjArray *>(GetCachedCDBObject(kIDDCS));
  if(!dcsArr){
    AliError("No DCS CDB Object available!");
    config = "";
    return;
  }

  // Check CalDCS version
  Int_t calver = 0;
  if (!strcmp(dcsArr->At(0)->ClassName(),"AliTRDCalDCS")) calver = 1;
  else if (!strcmp(dcsArr->At(0)->ClassName(),"AliTRDCalDCSv2")) calver = 2;
  //
  // Get the configuration strings
  //
  // CalDCS version 1
  if (calver == 1) {

    // Loop over SOR/EOR
    for(Int_t iSEOR=0;iSEOR<2;iSEOR++){

      // DCS object
      AliTRDCalDCS *cal = dynamic_cast<AliTRDCalDCS *>(dcsArr->At(iSEOR));
      if(!cal){
	// SOR mandantory
	if(iSEOR==0){
	  AliError("NO SOR object found in CDB file!");
	  config = "";
	  return;
	}
	else{
	  AliWarning("NO EOR object found in CDB file!");
	  continue;
	}
      } // Check DCS object is there
      
      // Loop over 540 chambers
      {
	Int_t iDet=0;
	// Find the first chamber in SOR
	if(iSEOR==0){
	  // Loop until we find the first chamber
	  for(;iDet<kNdet;iDet++){
	    const AliTRDCalDCSFEE *DCSFEEObj = cal->GetCalDCSFEEObj(iDet);
	    // Check it's an installed chamber that responded
	    if(DCSFEEObj && !DCSFEEObj->GetStatusBit()) {

	      // We ask for a valid configuration name starting with cf_
	      if(DCSFEEObj->GetConfigName().BeginsWith("cf_")){
		// Set the parameter of the first good ROC as global,
		// we take only the requested part!
		GetDCSConfigParOption(DCSFEEObj->GetConfigName(), par, opt, config);
		// Increase counter to not look at this chamber again
		// and break loop
		iDet++;
		break;
	      } // End: valid config name
	    } // End: first good chamber
	  } // End: find the first chamber
	} // End: is SOR

	// Now compare the other chambers with the first one
	for(;iDet<kNdet;iDet++){
	  const AliTRDCalDCSFEE *DCSFEEObj = cal->GetCalDCSFEEObj(iDet);
	  // Check it's an installed chamber that responded
	  if(DCSFEEObj && !DCSFEEObj->GetStatusBit()) {
	    
	    // We ask for a valid configuration name starting with cf_
	    if(DCSFEEObj->GetConfigName().BeginsWith("cf_")){

	      // Get the parameter/option of this chamber
	      TString tmpcfg;
	      GetDCSConfigParOption(DCSFEEObj->GetConfigName(), par, opt, tmpcfg);
	      // Compare with the global value
	      if(config.CompareTo(tmpcfg)){
		// Two cases mixed or changed during run
		if(iSEOR==0){
		  AliError("Mixed DCS configuration for different chambers!");
		  config="mixed";
		  return;
		} else {
		  AliError("Inconsistent configuration at start and end of run found!");
		  config = "";
		  return;
		}
	      }
	    }// Valid config name 
	  }// Good chamber
	} // Second half loop 
      } // Loop over 540 chambers
    } // Loop over SOR / EOR 
  } // End calver 1
  //
  // CalDCS version 2
  //
  else if (calver == 2) {

    // Loop over SOR/EOR
    for(Int_t iSEOR=0;iSEOR<2;iSEOR++){

      // DCS object
      AliTRDCalDCSv2 *cal = dynamic_cast<AliTRDCalDCSv2 *>(dcsArr->At(iSEOR));
      if(!cal){
	// SOR mandantory
	if(iSEOR==0){
	  AliError("NO SOR object found in CDB file!");
	  config = "";
	  return;
	}
	else{
	  AliWarning("NO EOR object found in CDB file!");
	  continue;
	}
      } // Check DCS object is there
      
      // Loop over 540 chambers
      {
	Int_t iDet=0;
	// Find the first chamber in SOR
	if(iSEOR==0){
	  // Loop until we find the first chamber
	  for(;iDet<kNdet;iDet++){
	    const AliTRDCalDCSFEEv2 *DCSFEEObj = cal->GetCalDCSFEEObj(iDet);
	    // Check it's an installed chamber that responded
	    if(DCSFEEObj && !DCSFEEObj->GetStatusBit()) {
	      // Check for a valid config name
	      if(DCSFEEObj->GetConfigName().BeginsWith("cf_")){

		// Set the parameter of the first good ROC as global,
		// we take only the requested part!
		GetDCSConfigParOption(DCSFEEObj->GetConfigName(), par, opt, config);
		// Increase counter to not look at this chamber again
		// and break loop
		iDet++;
		break;
	      } // End: valid config name
	    } // End: first good chamber
	  } // End: find the first chamber
	} // End: is SOR

	// Now compare the other chambers with the first one
	for(;iDet<kNdet;iDet++){
	  const AliTRDCalDCSFEEv2 *DCSFEEObj = cal->GetCalDCSFEEObj(iDet);
	  // Check it's an installed chamber that responded
	  if(DCSFEEObj && !DCSFEEObj->GetStatusBit()) {

	    // Check for a valid config name
	    if(DCSFEEObj->GetConfigName().BeginsWith("cf_")){

	      // Get the parameter/option of this chamber
	      TString tmpcfg;
	      GetDCSConfigParOption(DCSFEEObj->GetConfigName(), par, opt, tmpcfg);
	      // Compare with the global value
	      if(config.CompareTo(tmpcfg)){
		// Two cases mixed or changed during run
		if(iSEOR==0){
		  AliError("Mixed DCS configuration for different chambers!");
		  config="mixed";
		  return;
		} else {
		  AliError("Inconsistent configuration at start and end of run found!");
		  config = "";
		  return;
		}
	      }
	    } // End: valid config name
	  } // End: Good chamber
	} // End: Second half loop 
      } // Loop over 540 chambers
    } // Loop over SOR / EOR 
  } // End calver 2
  else {
    AliError("NO DCS/DCSv2 OCDB entry found!");
    config = "";
    return;
  }

}
//_____________________________________________________________________________
void AliTRDcalibDB::GetGlobalConfigurationVersion(TString &version)
{
  //
  // Get Version of Configuration from the DCS
  //

  const TObjArray *dcsArr = dynamic_cast<const TObjArray *>(GetCachedCDBObject(kIDDCS));
  if(!dcsArr){
    version = "";
    return;
  }

  Int_t esor   = 0; // Take SOR
  Int_t calver = 0; // Check CalDCS version
  if (!strcmp(dcsArr->At(0)->ClassName(),"AliTRDCalDCS"))   calver = 1;
  if (!strcmp(dcsArr->At(0)->ClassName(),"AliTRDCalDCSv2")) calver = 2;

  if      (calver == 1) {

    // DCS object
    const AliTRDCalDCS   *calDCS   = dynamic_cast<const AliTRDCalDCS *>(dcsArr->At(esor));
    if(!calDCS){
      version = "";
      return;
    } 
    version = calDCS->GetGlobalConfigVersion();

  } 
  else if (calver == 2) {

    // DCSv2 object
    const AliTRDCalDCSv2 *calDCSv2 = dynamic_cast<const AliTRDCalDCSv2 *>(dcsArr->At(esor));
    if(!calDCSv2){
      version = "";
      return;
    } 
    version = calDCSv2->GetGlobalConfigVersion();

  } 
  else {

    AliError("NO DCS/DCSv2 OCDB entry found!");

  }

}

//_____________________________________________________________________________
Int_t AliTRDcalibDB::GetNumberOfParsDCS(TString cname, Char_t delimiter)
{
  // Get the number of configuration parameters from the DCS config

  //AliInfo(Form("\"%s\" tokenized by \"%c\"", cname.Data(), delimiter));
  if(!cname.Length()) return -1;  // -1 for the "cf"
  Int_t nconf(0);
  for(Int_t ich(1); ich<cname.Length()-1; ich++){ if(cname[ich]==delimiter) nconf++;}
  return nconf;
}

//_____________________________________________________________________________
Int_t AliTRDcalibDB::GetNumberOfOptsDCS(TString cname, Int_t cfgType)
{
  // Get the number of options of a given configuration parameter from DCS

  Char_t cdelim = '_', // define the delimiters
         odelim = '-';
  Int_t nconfig = GetNumberOfParsDCS(cname, cdelim);

  // protect
  if ((nconfig == -1) || (nconfig < cfgType)) {
    AliErrorClass("Not enough parameters in DCS configuration name!");
    return 0;
  }

  TObjArray *carr = cname.Tokenize(cdelim);
  Int_t nopt = GetNumberOfParsDCS(((TObjString*)carr->At(cfgType))->GetString(), odelim);
  carr->Delete(); delete carr;
  return nopt;
}

//_____________________________________________________________________________
void AliTRDcalibDB::GetDCSConfigParOption(TString cname, Int_t cfgType, Int_t option, TString &cfgo)
{
  //
  // Get a configuration (see enum in header file) or the options of a configuration
  // option == 0 returns the configuration itself
  // option >  0 returns the optional parameter Nr. (option) of the configuration (cfgType)
  //

  Char_t cdelim = '_', // define the delimiters
         odelim = '-';

  Int_t nconfig = GetNumberOfParsDCS(cname, cdelim);
  // protect
  if (nconfig == -1) {
    AliErrorClass("DCS configuration name empty!");
    cfgo = "";
    return;
  } else if (nconfig < cfgType) {
    AliErrorClass(Form("Not enough parameters in DCS configuration name!"
		  " Name %s",cname.Data()));
    cfgo = "";
    return;
  }

  TObjArray *carr = cname.Tokenize(cdelim);
  TString cfgString(((TObjString*)carr->At(cfgType))->GetString());
  Int_t noptions = GetNumberOfParsDCS(cfgString, odelim);
  // protect
  if (noptions < option) {
    AliErrorClass(Form("Not enough options in DCS configuration name!"
		  " Name %s",cname.Data()));
    cfgo = "";
    carr->Delete(); delete carr;
    return;
  }
  TObjArray *oarr = cfgString.Tokenize(odelim);
  cfgo = ((TObjString*)oarr->At(option))->GetString();
  carr->Delete(); delete carr;
  oarr->Delete(); delete oarr;
  return;

}

//_____________________________________________________________________________
Bool_t AliTRDcalibDB::HasOnlineFilterPedestal()
{
  //
  // Checks whether pedestal filter was applied online
  //

  TString filterconfig;
  GetFilterType(filterconfig);

  return filterconfig.Contains("p");

}

//_____________________________________________________________________________
Bool_t AliTRDcalibDB::HasOnlineFilterGain()
{
  //
  // Checks whether online gain filter was applied
  //

  TString filterconfig;
  GetFilterType(filterconfig);

  return filterconfig.Contains("g");

}

//_____________________________________________________________________________
Bool_t AliTRDcalibDB::HasOnlineTailCancellation()
{
  //
  // Checks whether online tail cancellation was applied
  //

  TString filterconfig;
  GetFilterType(filterconfig);

  return filterconfig.Contains("t");

}

//_____________________________________________________________________________
Char_t AliTRDcalibDB::GetPadStatus(Int_t det, Int_t col, Int_t row)
{
  //
  // Returns the status of the given pad
  //

  const AliTRDCalPadStatus *cal  = dynamic_cast<const AliTRDCalPadStatus *> 
                                   (GetCachedCDBObject(kIDPadStatus));
  if (!cal) {
    return -1;
  }

  const AliTRDCalSingleChamberStatus *roc = cal->GetCalROC(det);
  if (!roc) {
    return -1;
  }

  return roc->GetStatus(col,row);

}

//_____________________________________________________________________________
AliTRDCalSingleChamberStatus* AliTRDcalibDB::GetPadStatusROC(Int_t det)
{
  //
  // Returns the pad status calibration object for a given ROC
  //

  const AliTRDCalPadStatus *cal  = dynamic_cast<const AliTRDCalPadStatus *> 
                                   (GetCachedCDBObject(kIDPadStatus));
  if (!cal) {
    return 0;
  }

  AliTRDCalSingleChamberStatus *roc = cal->GetCalROC(det);
  if (!roc) {
    return 0;
  }
  else {
    return roc;
  }

}

//_____________________________________________________________________________
Char_t AliTRDcalibDB::GetChamberStatus(Int_t det)
{
  //
  // Returns the status of the given chamber
  //

  const AliTRDCalChamberStatus *cal = dynamic_cast<const AliTRDCalChamberStatus *> 
                                      (GetCachedCDBObject(kIDChamberStatus));
  if (!cal) {
    return -1;
  }

  return cal->GetStatus(det);

}

//_____________________________________________________________________________
AliTRDrecoParam* AliTRDcalibDB::GetRecoParam(Int_t */*eventtype*/)
{
  //
  // Returns the TRD reconstruction parameters from the OCDB
  //

  const TClonesArray *recos = dynamic_cast<const TClonesArray*>(GetCachedCDBObject(kIDRecoParam));
  if (!recos) return 0x0;

  // calculate entry based on event type info
  Int_t n = 0; //f(eventtype[0], eventtype[1], ....)

  return (AliTRDrecoParam *) recos->UncheckedAt(n);

}

//_____________________________________________________________________________
Bool_t AliTRDcalibDB::IsPadMasked(Int_t det, Int_t col, Int_t row)
{
  //
  // Returns status, see name of functions for details ;-)
  //

  const AliTRDCalPadStatus         *cal = dynamic_cast<const AliTRDCalPadStatus *> 
                                          (GetCachedCDBObject(kIDPadStatus));
  if (!cal) {
    return -1;
  }

  return cal->IsMasked(det,col,row);

}

//_____________________________________________________________________________
Bool_t AliTRDcalibDB::IsPadBridgedLeft(Int_t det, Int_t col, Int_t row)
{
  //
  // Returns status, see name of functions for details ;-)
  //

  const AliTRDCalPadStatus         *cal = dynamic_cast<const AliTRDCalPadStatus *> 
                                          (GetCachedCDBObject(kIDPadStatus));
  if (!cal) {
    return -1;
  }

  return cal->IsBridgedLeft(det,col,row);

}

//_____________________________________________________________________________
Bool_t AliTRDcalibDB::IsPadBridgedRight(Int_t det, Int_t col, Int_t row)
{
  //
  // Returns status, see name of functions for details ;-)
  //

  const AliTRDCalPadStatus         * cal = dynamic_cast<const AliTRDCalPadStatus *> 
                                           (GetCachedCDBObject(kIDPadStatus));
  if (!cal) {
    return -1;
  }

  return cal->IsBridgedRight(det,col,row);

}

//_____________________________________________________________________________
Bool_t AliTRDcalibDB::IsPadNotConnected(Int_t det, Int_t col, Int_t row)
{
  //
  // Returns status, see name of functions for details ;-)
  //

  const AliTRDCalPadStatus         * cal = dynamic_cast<const AliTRDCalPadStatus *> 
                                           (GetCachedCDBObject(kIDPadStatus));
  if (!cal) {
    return -1;
  }

  return cal->IsNotConnected(det,col,row);

}

//_____________________________________________________________________________
Bool_t AliTRDcalibDB::IsChamberGood(Int_t det)
{
  //
  // Returns status, see name of functions for details ;-)
  //

  const AliTRDCalChamberStatus     * cal = dynamic_cast<const AliTRDCalChamberStatus *> 
                                           (GetCachedCDBObject(kIDChamberStatus));
  if (!cal) {
    return -1;
  }

  return cal->IsGood(det);

}

//_____________________________________________________________________________
Bool_t AliTRDcalibDB::IsChamberNoData(Int_t det)
{
  //
  // Returns status, see name of functions for details ;-)
  //

  const AliTRDCalChamberStatus     * cal = dynamic_cast<const AliTRDCalChamberStatus *> 
                                           (GetCachedCDBObject(kIDChamberStatus));
  if (!cal) {
    return -1;
  }

  return cal->IsNoData(det);

}

//_____________________________________________________________________________
Bool_t AliTRDcalibDB::IsHalfChamberNoData(Int_t det, Int_t side)
{
  //
  // Returns status, see name of functions for details ;-)
  //

  const AliTRDCalChamberStatus     * cal = dynamic_cast<const AliTRDCalChamberStatus *> 
                                           (GetCachedCDBObject(kIDChamberStatus));
  if (!cal) {
    return -1;
  }

  return side > 0 ? cal->IsNoDataSideB(det) : cal->IsNoDataSideA(det);

}

//_____________________________________________________________________________
Bool_t AliTRDcalibDB::IsChamberBadCalibrated(Int_t det)
{
  //
  // Returns status, see name of functions for details ;-)
  //

  const AliTRDCalChamberStatus     * cal = dynamic_cast<const AliTRDCalChamberStatus *> 
                                           (GetCachedCDBObject(kIDChamberStatus));
  if (!cal) {
    return -1;
  }

  return cal->IsBadCalibrated(det);

}

//_____________________________________________________________________________
Bool_t AliTRDcalibDB::IsChamberNotCalibrated(Int_t det)
{
  //
  // Returns status, see name of functions for details ;-)
  //

  const AliTRDCalChamberStatus     * cal = dynamic_cast<const AliTRDCalChamberStatus *> 
                                           (GetCachedCDBObject(kIDChamberStatus));
  if (!cal) {
    return -1;
  }

  return cal->IsNotCalibrated(det);

}

//_____________________________________________________________________________
const AliTRDCalPID *AliTRDcalibDB::GetPIDObject(AliTRDpidUtil::ETRDPIDMethod method)
{
  //
  // Returns the object storing the distributions for PID with likelihood
  //

  switch(method) {
  case AliTRDpidUtil::kLQ: 
    return dynamic_cast<const AliTRDCalPID *>(GetCachedCDBObject(kIDPIDLQ));
  case AliTRDpidUtil::kNN: 
    return dynamic_cast<const AliTRDCalPID *>(GetCachedCDBObject(kIDPIDNN));
  case AliTRDpidUtil::kESD:
    return 0x0; // To avoid compiler warnings 
  }

  return 0x0;

}

//_____________________________________________________________________________
AliTRDPIDResponse *AliTRDcalibDB::GetPIDResponse(AliTRDPIDResponse::ETRDPIDMethod /*method*/)
{
  //
  // Returns the PID response object for 1D-LQ
  //

  if (!fPIDResponse) {
    AliDebug(2, "Setting new PID response.");

    fPIDResponse = new AliTRDPIDResponse;

    // Load Reference Histos from OCDB
//    if(method == AliTRDPIDResponse::kLQ1D){
    //fPIDResponse->SetPIDmethod(AliTRDPIDResponse::kLQ1D);
    const TObjArray *references = dynamic_cast<const TObjArray *>(GetCachedCDBObject(kIDPIDLQ1D));

    TIter refs(references);
    TObject *obj = NULL;
    AliTRDPIDReference *ref = NULL;
    Bool_t hasReference = kFALSE;
    while ((obj = refs())){
      if ((ref = dynamic_cast<AliTRDPIDReference *>(obj))){
        AliDebug(2, "Setting new PID response object.");
        TDirectory *bkpdir = gDirectory;
        gROOT->cd();
        AliTRDPIDResponseObject *ro = new AliTRDPIDResponseObject;
        ro->SetPIDReference(ref);
        fPIDResponse->SetPIDResponseObject(ro);
        hasReference = kTRUE;
        gDirectory = bkpdir;
        break;
      }
    }

    if (!hasReference) {
      AliError("Reference histograms not found in the OCDB");
    }
  }

//  }

  return fPIDResponse;

}

//_____________________________________________________________________________
const AliTRDCalTrkAttach* AliTRDcalibDB::GetAttachObject()
{
  //
  // Returns the object storing likelihood distributions for cluster to track attachment
  //

  return dynamic_cast<const AliTRDCalTrkAttach*>(GetCachedCDBObject(kIDAttach));

}

//_____________________________________________________________________________
const AliTRDCalMonitoring *AliTRDcalibDB::GetMonitoringObject()
{
  //
  // Returns the object storing the monitoring data
  //

  return dynamic_cast<const AliTRDCalMonitoring *> 
         (GetCachedCDBObject(kIDMonitoringData));
   
}

//_____________________________________________________________________________
void AliTRDcalibDB::SamplePRF()
{
  //
  // Samples the pad response function (should maybe go somewhere else ...)
  //

  const Int_t kPRFbin = 61;

  Float_t prf[kNlayer][kPRFbin] = { 
                   {2.9037e-02, 3.3608e-02, 3.9020e-02, 4.5292e-02,
                    5.2694e-02, 6.1362e-02, 7.1461e-02, 8.3362e-02,
                    9.7063e-02, 1.1307e-01, 1.3140e-01, 1.5235e-01,
                    1.7623e-01, 2.0290e-01, 2.3294e-01, 2.6586e-01,
                    3.0177e-01, 3.4028e-01, 3.8077e-01, 4.2267e-01,
                    4.6493e-01, 5.0657e-01, 5.4655e-01, 5.8397e-01,
                    6.1767e-01, 6.4744e-01, 6.7212e-01, 6.9188e-01,
                    7.0627e-01, 7.1499e-01, 7.1851e-01, 7.1499e-01,
                    7.0627e-01, 6.9188e-01, 6.7212e-01, 6.4744e-01,
                    6.1767e-01, 5.8397e-01, 5.4655e-01, 5.0657e-01,
                    4.6493e-01, 4.2267e-01, 3.8077e-01, 3.4028e-01,
                    3.0177e-01, 2.6586e-01, 2.3294e-01, 2.0290e-01,
                    1.7623e-01, 1.5235e-01, 1.3140e-01, 1.1307e-01,
                    9.7063e-02, 8.3362e-02, 7.1461e-02, 6.1362e-02,
                    5.2694e-02, 4.5292e-02, 3.9020e-02, 3.3608e-02,
                    2.9037e-02},
                   {2.5478e-02, 2.9695e-02, 3.4655e-02, 4.0454e-02,
                    4.7342e-02, 5.5487e-02, 6.5038e-02, 7.6378e-02,
                    8.9696e-02, 1.0516e-01, 1.2327e-01, 1.4415e-01,
                    1.6794e-01, 1.9516e-01, 2.2573e-01, 2.5959e-01,
                    2.9694e-01, 3.3719e-01, 3.7978e-01, 4.2407e-01,
                    4.6889e-01, 5.1322e-01, 5.5569e-01, 5.9535e-01,
                    6.3141e-01, 6.6259e-01, 6.8882e-01, 7.0983e-01,
                    7.2471e-01, 7.3398e-01, 7.3761e-01, 7.3398e-01,
                    7.2471e-01, 7.0983e-01, 6.8882e-01, 6.6259e-01,
                    6.3141e-01, 5.9535e-01, 5.5569e-01, 5.1322e-01,
                    4.6889e-01, 4.2407e-01, 3.7978e-01, 3.3719e-01,
                    2.9694e-01, 2.5959e-01, 2.2573e-01, 1.9516e-01,
                    1.6794e-01, 1.4415e-01, 1.2327e-01, 1.0516e-01,
                    8.9696e-02, 7.6378e-02, 6.5038e-02, 5.5487e-02,
                    4.7342e-02, 4.0454e-02, 3.4655e-02, 2.9695e-02,
                    2.5478e-02},
                   {2.2363e-02, 2.6233e-02, 3.0782e-02, 3.6140e-02,
                    4.2535e-02, 5.0157e-02, 5.9197e-02, 6.9900e-02,
                    8.2707e-02, 9.7811e-02, 1.1548e-01, 1.3601e-01,
                    1.5998e-01, 1.8739e-01, 2.1840e-01, 2.5318e-01,
                    2.9182e-01, 3.3373e-01, 3.7837e-01, 4.2498e-01,
                    4.7235e-01, 5.1918e-01, 5.6426e-01, 6.0621e-01,
                    6.4399e-01, 6.7700e-01, 7.0472e-01, 7.2637e-01,
                    7.4206e-01, 7.5179e-01, 7.5551e-01, 7.5179e-01,
                    7.4206e-01, 7.2637e-01, 7.0472e-01, 6.7700e-01,
                    6.4399e-01, 6.0621e-01, 5.6426e-01, 5.1918e-01,
                    4.7235e-01, 4.2498e-01, 3.7837e-01, 3.3373e-01,
                    2.9182e-01, 2.5318e-01, 2.1840e-01, 1.8739e-01,
                    1.5998e-01, 1.3601e-01, 1.1548e-01, 9.7811e-02,
                    8.2707e-02, 6.9900e-02, 5.9197e-02, 5.0157e-02,
                    4.2535e-02, 3.6140e-02, 3.0782e-02, 2.6233e-02,
                    2.2363e-02},
                   {1.9635e-02, 2.3167e-02, 2.7343e-02, 3.2293e-02,
                    3.8224e-02, 4.5335e-02, 5.3849e-02, 6.4039e-02,
                    7.6210e-02, 9.0739e-02, 1.0805e-01, 1.2841e-01,
                    1.5216e-01, 1.7960e-01, 2.1099e-01, 2.4671e-01,
                    2.8647e-01, 3.2996e-01, 3.7660e-01, 4.2547e-01,
                    4.7536e-01, 5.2473e-01, 5.7215e-01, 6.1632e-01,
                    6.5616e-01, 6.9075e-01, 7.1939e-01, 7.4199e-01,
                    7.5838e-01, 7.6848e-01, 7.7227e-01, 7.6848e-01,
                    7.5838e-01, 7.4199e-01, 7.1939e-01, 6.9075e-01,
                    6.5616e-01, 6.1632e-01, 5.7215e-01, 5.2473e-01,
                    4.7536e-01, 4.2547e-01, 3.7660e-01, 3.2996e-01,
                    2.8647e-01, 2.4671e-01, 2.1099e-01, 1.7960e-01,
                    1.5216e-01, 1.2841e-01, 1.0805e-01, 9.0739e-02,
                    7.6210e-02, 6.4039e-02, 5.3849e-02, 4.5335e-02,
                    3.8224e-02, 3.2293e-02, 2.7343e-02, 2.3167e-02,
                    1.9635e-02},
                   {1.7224e-02, 2.0450e-02, 2.4286e-02, 2.8860e-02,
                    3.4357e-02, 4.0979e-02, 4.8966e-02, 5.8612e-02,
                    7.0253e-02, 8.4257e-02, 1.0102e-01, 1.2094e-01,
                    1.4442e-01, 1.7196e-01, 2.0381e-01, 2.4013e-01,
                    2.8093e-01, 3.2594e-01, 3.7450e-01, 4.2563e-01,
                    4.7796e-01, 5.2991e-01, 5.7974e-01, 6.2599e-01,
                    6.6750e-01, 7.0344e-01, 7.3329e-01, 7.5676e-01,
                    7.7371e-01, 7.8410e-01, 7.8793e-01, 7.8410e-01,
                    7.7371e-01, 7.5676e-01, 7.3329e-01, 7.0344e-01,
                    6.6750e-01, 6.2599e-01, 5.7974e-01, 5.2991e-01,
                    4.7796e-01, 4.2563e-01, 3.7450e-01, 3.2594e-01,
                    2.8093e-01, 2.4013e-01, 2.0381e-01, 1.7196e-01,
                    1.4442e-01, 1.2094e-01, 1.0102e-01, 8.4257e-02,
                    7.0253e-02, 5.8612e-02, 4.8966e-02, 4.0979e-02,
                    3.4357e-02, 2.8860e-02, 2.4286e-02, 2.0450e-02,
                    1.7224e-02},
                   {1.5096e-02, 1.8041e-02, 2.1566e-02, 2.5793e-02,
                    3.0886e-02, 3.7044e-02, 4.4515e-02, 5.3604e-02,
                    6.4668e-02, 7.8109e-02, 9.4364e-02, 1.1389e-01,
                    1.3716e-01, 1.6461e-01, 1.9663e-01, 2.3350e-01,
                    2.7527e-01, 3.2170e-01, 3.7214e-01, 4.2549e-01,
                    4.8024e-01, 5.3460e-01, 5.8677e-01, 6.3512e-01,
                    6.7838e-01, 7.1569e-01, 7.4655e-01, 7.7071e-01,
                    7.8810e-01, 7.9871e-01, 8.0255e-01, 7.9871e-01,
                    7.8810e-01, 7.7071e-01, 7.4655e-01, 7.1569e-01,
                    6.7838e-01, 6.3512e-01, 5.8677e-01, 5.3460e-01,
                    4.8024e-01, 4.2549e-01, 3.7214e-01, 3.2170e-01,
                    2.7527e-01, 2.3350e-01, 1.9663e-01, 1.6461e-01,
                    1.3716e-01, 1.1389e-01, 9.4364e-02, 7.8109e-02,
                    6.4668e-02, 5.3604e-02, 4.4515e-02, 3.7044e-02,
                    3.0886e-02, 2.5793e-02, 2.1566e-02, 1.8041e-02,
                    1.5096e-02}};

  // More sampling precision with linear interpolation
  fPRFlo  = -1.5;
  fPRFhi  =  1.5;
  Float_t pad[kPRFbin];
  Int_t   sPRFbin = kPRFbin;  
  Float_t sPRFwid = (fPRFhi - fPRFlo) / ((Float_t) sPRFbin);
  for (Int_t iPad = 0; iPad < sPRFbin; iPad++) {
    pad[iPad] = ((Float_t) iPad + 0.5) * sPRFwid + fPRFlo;
  }
  fPRFbin = 500;  
  fPRFwid = (fPRFhi - fPRFlo) / ((Float_t) fPRFbin);
  fPRFpad = ((Int_t) (1.0 / fPRFwid));

  if (fPRFsmp) delete [] fPRFsmp;
  fPRFsmp = new Float_t[kNlayer*fPRFbin];

  Int_t   ipos1;
  Int_t   ipos2;
  Float_t diff;

  for (Int_t iLayer = 0; iLayer < kNlayer; iLayer++) {

    for (Int_t iBin = 0; iBin < fPRFbin; iBin++) {

      Float_t bin = (((Float_t) iBin) + 0.5) * fPRFwid + fPRFlo;
      ipos1 = ipos2 = 0;
      diff  = 0;
      do {
        diff = bin - pad[ipos2++];
      } while ((diff > 0) && (ipos2 < kPRFbin));
      if      (ipos2 == kPRFbin) {
        fPRFsmp[iLayer*fPRFbin+iBin] = prf[iLayer][ipos2-1];
      }
      else if (ipos2 == 1) {
        fPRFsmp[iLayer*fPRFbin+iBin] = prf[iLayer][ipos2-1];
      }
      else {
        ipos2--;
        if (ipos2 >= kPRFbin) ipos2 = kPRFbin - 1;
        ipos1 = ipos2 - 1;
        fPRFsmp[iLayer*fPRFbin+iBin] = prf[iLayer][ipos2] 
                                     + diff * (prf[iLayer][ipos2] - prf[iLayer][ipos1]) 
                                     / sPRFwid;
      }

    }
  } 

}

//_____________________________________________________________________________
Int_t AliTRDcalibDB::PadResponse(Double_t signal, Double_t dist
                               , Int_t layer, Double_t *pad) const
{
  //
  // Applies the pad response
  // So far this is the fixed parametrization and should be replaced by
  // something dependent on calibration values
  //

  Int_t iBin  = ((Int_t) ((-dist - fPRFlo) / fPRFwid));
  Int_t iOff  = layer * fPRFbin;

  Int_t iBin0 = iBin - fPRFpad + iOff;
  Int_t iBin1 = iBin           + iOff;
  Int_t iBin2 = iBin + fPRFpad + iOff;

  pad[0] = 0.0;
  pad[1] = 0.0;
  pad[2] = 0.0;
  if ((iBin1 >= 0) && (iBin1 < (fPRFbin*kNlayer))) {

    if (iBin0 >= 0) {
      pad[0] = signal * fPRFsmp[iBin0];
    }
    pad[1] = signal * fPRFsmp[iBin1];
    if (iBin2 < (fPRFbin*kNlayer)) {
      pad[2] = signal * fPRFsmp[iBin2];
    }

    return 1;

  }
  else {

    return 0;

  }

}


AliTRDtrapConfig* AliTRDcalibDB::GetTrapConfig()
{
  // return an existing TRAPconfig or load it from the OCDB
  // in case of failure, a default TRAPconfig is created

  if (fTrapConfig)
    return fTrapConfig;
  else {
    if ((fTrapConfigName.Length() <= 0) || (fTrapConfigVersion.Length() <= 0)) {
      // query the configuration to be used
      this->GetGlobalConfiguration(fTrapConfigName);
      this->GetGlobalConfigurationVersion(fTrapConfigVersion);
    }

    // try to load the requested configuration
    this->LoadTrapConfig(fTrapConfigName, fTrapConfigVersion);

    // if we still don't have a valid TRAPconfig, we give up
    if (!fTrapConfig)
      AliFatal("Requested TRAP configuration not found!");

    AliInfo(Form("using TRAPconfig \"%s\"", fTrapConfig->GetTitle()));

    // we still have to load the gain tables
    // if the gain filter is active
    if (HasOnlineFilterGain()) {
      const Int_t nDets = AliTRDgeometry::Ndet();
      const Int_t nMcms = AliTRDgeometry::MCMmax();
      const Int_t nChs  = AliTRDgeometry::ADCmax();

      // gain factors are per MCM
      // allocate the registers accordingly
      for (Int_t ch = 0; ch < nChs; ++ch) {
      	AliTRDtrapConfig::TrapReg_t regFGAN = (AliTRDtrapConfig::TrapReg_t) (AliTRDtrapConfig::kFGA0 + ch);
      	AliTRDtrapConfig::TrapReg_t regFGFN = (AliTRDtrapConfig::TrapReg_t) (AliTRDtrapConfig::kFGF0 + ch);
      	fTrapConfig->SetTrapRegAlloc(regFGAN, AliTRDtrapConfig::kAllocByMCM);
      	fTrapConfig->SetTrapRegAlloc(regFGFN, AliTRDtrapConfig::kAllocByMCM);
      }

      for (Int_t iDet = 0; iDet < nDets; ++iDet) {
      	AliTRDCalOnlineGainTableROC *gainTbl = GetOnlineGainTableROC(iDet);
      	if (gainTbl) {
	  const Int_t nRobs = AliTRDgeometry::GetStack(iDet) == 2 ?
	    AliTRDgeometry::ROBmaxC0() : AliTRDgeometry::ROBmaxC1();
      	  for (Int_t rob = 0; rob < nRobs; ++rob) {
	    for (Int_t mcm = 0; mcm < nMcms; ++mcm) {
	      AliTRDCalOnlineGainTableMCM *gainTblMCM = gainTbl->GetGainTableMCM(rob, mcm);

	      // set ADC reference voltage
	      fTrapConfig->SetTrapReg(AliTRDtrapConfig::kADCDAC, gainTblMCM->GetAdcdac(), iDet, rob, mcm);

	      // set constants channel-wise
	      for (Int_t ch = 0; ch < nChs; ++ch) {
		AliTRDtrapConfig::TrapReg_t regFGAN = (AliTRDtrapConfig::TrapReg_t) (AliTRDtrapConfig::kFGA0 + ch);
		AliTRDtrapConfig::TrapReg_t regFGFN = (AliTRDtrapConfig::TrapReg_t) (AliTRDtrapConfig::kFGF0 + ch);
		fTrapConfig->SetTrapReg(regFGAN, gainTblMCM->GetFGAN(ch), iDet, rob, mcm);
		fTrapConfig->SetTrapReg(regFGFN, gainTblMCM->GetFGFN(ch), iDet, rob, mcm);
	      }
      	    }
      	  }
      	}
      }
    }

    return fTrapConfig;
  }
}


AliTRDtrapConfig* AliTRDcalibDB::LoadTrapConfig(const TString &name, const TString &version)
{
  // try to load the specified configuration from the OCDB
  // if it fails, or it does not exist, return null

  AliInfo(Form("looking for TRAPconfig \"%s.%s\"", name.Data(), version.Data()));

  const AliTRDCalTrapConfig *caltrap = dynamic_cast<const AliTRDCalTrapConfig*> (GetCachedCDBObject(kIDTrapConfig));

  if (caltrap) {
    TString configName(name);
    configName.Append(".");
    configName.Append(version);
    fTrapConfig = caltrap->Get(configName);
  }
  else {
    fTrapConfig = 0x0;
    AliError("No TRAPconfig entry found!");
  }

  return fTrapConfig;
}
