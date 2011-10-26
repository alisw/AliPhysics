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

#include <TClonesArray.h>
#include <TObjArray.h>

#include "AliCDBManager.h"
#include "AliCDBEntry.h"
#include "AliLog.h"

#include "AliTRDPIDReference.h"
#include "AliTRDcalibDB.h"
#include "AliTRDCommonParam.h"

#include "Cal/AliTRDCalROC.h"
#include "Cal/AliTRDCalPad.h"
#include "Cal/AliTRDCalDet.h"
#include "Cal/AliTRDCalDCS.h"
#include "Cal/AliTRDCalDCSv2.h"
#include "Cal/AliTRDCalDCSFEEv2.h"
#include "Cal/AliTRDCalPID.h"
#include "Cal/AliTRDCalMonitoring.h"
#include "Cal/AliTRDCalChamberStatus.h"
#include "Cal/AliTRDCalPadStatus.h"
#include "Cal/AliTRDCalSingleChamberStatus.h"
#include "Cal/AliTRDCalTrkAttach.h"
#include "Cal/AliTRDCalOnlineGainTable.h"

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
{
  //
  // Default constructor
  //
  // TODO Default runnumber is set to 0, this should be changed later
  //      to an invalid value (e.g. -1) to prevent
  // TODO invalid calibration data to be used.
  //

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
	  //printf("vdrift %f and exb %f\n",calChambervdrift->GetValue(k),calDetExB->GetValue(k));
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
AliTRDCalOnlineGainTableROC* AliTRDcalibDB::GetOnlineGainTableROC(Int_t det)
{
  //
  // Returns the online gain factor table for a given ROC.
  //
  
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
Int_t AliTRDcalibDB::GetNumberOfTimeBinsDCS()
{
  //
  // Returns Number of time bins from the DCS
  //

  Int_t nMixed = -2; // not the same number for all chambers
  Int_t nUndef = -1; // default value - has not been set!
  Int_t nTbSor = nUndef;
  Int_t nTbEor = nUndef;
  Int_t calver = 0; // Check CalDCS version

  const TObjArray *dcsArr = dynamic_cast<const TObjArray *>(GetCachedCDBObject(kIDDCS));
  if (!dcsArr) {
    AliError("No DCS object found!");
    return nUndef;
  }

  if (!strcmp(dcsArr->At(0)->ClassName(),"AliTRDCalDCS"))   calver = 1;
  if (!strcmp(dcsArr->At(0)->ClassName(),"AliTRDCalDCSv2")) calver = 2;

  if (calver == 1) {
    // DCS object
    const AliTRDCalDCS *calDCSsor = dynamic_cast<const AliTRDCalDCS *>(dcsArr->At(0));
    const AliTRDCalDCS *calDCSeor = dynamic_cast<const AliTRDCalDCS *>(dcsArr->At(1));
    if (!calDCSsor) {
      // the SOR file is mandatory
      AliError("NO SOR AliTRDCalDCS object found in CDB file!");
      return nUndef;
    }
    if (!calDCSeor) {
      // this can happen if the run is shorter than a couple of seconds.
      AliWarning("NO EOR AliTRDCalDCS object found in CDB file.");
    }

    // get the numbers
    nTbSor = calDCSsor->GetGlobalNumberOfTimeBins();
    if (calDCSeor) nTbEor = calDCSeor->GetGlobalNumberOfTimeBins();

  } else if (calver == 2) {
    // DCSv2 object
    const AliTRDCalDCSv2 *calDCSsorv2 = dynamic_cast<const AliTRDCalDCSv2 *>(dcsArr->At(0));
    const AliTRDCalDCSv2 *calDCSeorv2 = dynamic_cast<const AliTRDCalDCSv2 *>(dcsArr->At(1));
    if (!calDCSsorv2) {
      // the SOR file is mandatory
      AliError("NO SOR AliTRDCalDCSv2 object found in CDB file!");
      return nUndef;
    }
    if (!calDCSeorv2) {
      // this can happen if the run is shorter than a couple of seconds.
      AliWarning("NO EOR AliTRDCalDCSv2 object found in CDB file.");
    }

    // get the numbers
    nTbSor = calDCSsorv2->GetGlobalNumberOfTimeBins();
    if (calDCSeorv2) nTbEor = calDCSeorv2->GetGlobalNumberOfTimeBins();

  } else AliError("NO DCS/DCSv2 OCDB entry found!");

  // if they're the same return the value
  // -2 means mixed, -1: no data, >= 0: good number of time bins
  if (nTbSor == nTbEor) return nTbSor;

  // if they're differing:
  if (nTbSor == nMixed || nTbEor == nMixed) {
    AliWarning("Inconsistent number of time bins found!");
    return nMixed;
  }

  // one is undefined, the other ok -> return that one
  if (nTbSor == nUndef) return nTbEor;
  if (nTbEor == nUndef) return nTbSor;

  // only remains: two different numbers >= 0
  return nMixed;

}

//_____________________________________________________________________________
void AliTRDcalibDB::GetFilterType(TString &filterType)
{
  //
  // Returns the filter type
  //

  const TObjArray *dcsArr = dynamic_cast<const TObjArray *>(GetCachedCDBObject(kIDDCS));
  if(!dcsArr){
    filterType = "";
    return;
  }

  Int_t esor   = 0; // Take SOR
  Int_t calver = 0; // Check CalDCS version
  if (!strcmp(dcsArr->At(0)->ClassName(),"AliTRDCalDCS"))   calver = 1;
  if (!strcmp(dcsArr->At(0)->ClassName(),"AliTRDCalDCSv2")) calver = 2;

  if      (calver == 1) {

    // DCS object
    const AliTRDCalDCS *calDCS = dynamic_cast<const AliTRDCalDCS *>(dcsArr->At(esor));
    if(!calDCS){
      filterType = "";
      return;
    } 
    filterType = calDCS->GetGlobalFilterType();

  } 
  else if (calver == 2) {

    // DCSv2 object
    const AliTRDCalDCSv2 *calDCSv2 = dynamic_cast<const AliTRDCalDCSv2 *>(dcsArr->At(esor));
    if(!calDCSv2){
      filterType = "";
      return;
    } 
    filterType = calDCSv2->GetGlobalFilterType();

  } 
  else {

    AliError("NO DCS/DCSv2 OCDB entry found!");

  }

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
      const AliTRDCalDCSFEEv2 *calDCSFEEv2 = calDCSv2->GetCalDCSFEEObj(0);
      tableName = calDCSFEEv2->GetGainTableName();
      if (tableName.Length() > 0) {
        break;
      }
    }
    if (tableName.CompareTo("Krypton_2011-01") == 0) {
      fOnlineGainTableID = 1;
      return fOnlineGainTableID;
    }

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
    config = "";
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
      config = "";
      return;
    } 
    config = calDCS->GetGlobalConfigName();

  } 
  else if (calver == 2) {

    // DCSv2 object
    const AliTRDCalDCSv2 *calDCSv2 = dynamic_cast<const AliTRDCalDCSv2 *>(dcsArr->At(esor));
    if(!calDCSv2){
      config = "";
      return;
    } 
    config = calDCSv2->GetGlobalConfigName();

  } 
  else {

    AliError("NO DCS/DCSv2 OCDB entry found!");

  }

}

//_____________________________________________________________________________
Bool_t AliTRDcalibDB::HasOnlineFilterPedestal()
{
  //
  // Checks whether pedestal filter was applied online
  //

  TString cname;

  // Temporary: Get the filter config from the configuration name
  GetGlobalConfiguration(cname);
  TString filterconfig = cname(cname.First("_") + 1, cname.First("-") - cname.First("_") - 1);

  // TString filterconfig;
  //GetFilterType(filterconfig);

  return filterconfig.Contains("p");

}

//_____________________________________________________________________________
Bool_t AliTRDcalibDB::HasOnlineFilterGain()
{
  //
  // Checks whether online gain filter was applied
  //

  TString cname;

  // Temporary: Get the filter config from the configuration name
  GetGlobalConfiguration(cname);
  TString filterconfig = cname(cname.First("_") + 1, cname.First("-") - cname.First("_") - 1);

  //TString filterconfig;
  //GetFilterType(filterconfig);

  return filterconfig.Contains("g");

}

//_____________________________________________________________________________
Bool_t AliTRDcalibDB::HasOnlineTailCancellation()
{
  //
  // Checks whether online tail cancellation was applied
  //

  TString cname;

  // Temporary: Get the filter config from the configuration name
  GetGlobalConfiguration(cname);
  TString filterconfig = cname(cname.First("_") + 1, cname.First("-") - cname.First("_") - 1);

  //TString filterconfig;
  //GetFilterType(filterconfig);

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
AliTRDPIDResponse *AliTRDcalibDB::GetPIDResponse(AliTRDPIDResponse::ETRDPIDMethod method)
{
  //
  // Returns the PID response object for 1D-LQ
  //

  if (!fPIDResponse) {

    fPIDResponse = new AliTRDPIDResponse;

    // Load Reference Histos from OCDB
    fPIDResponse->SetPIDmethod(method);
    const TObjArray *references = dynamic_cast<const TObjArray *>(GetCachedCDBObject(kIDPIDLQ1D));

    TIter refs(references);
    TObject *obj = NULL;
    AliTRDPIDReference *ref = NULL;
    Bool_t hasReference = kFALSE;
    while ((obj = refs())){
      if ((ref = dynamic_cast<AliTRDPIDReference *>(obj))){
        fPIDResponse->Load(ref);
        hasReference = kTRUE;
        break;
      }
    }

    if (!hasReference) {
      AliError("Reference histograms not found in the OCDB");
    }

  }

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


