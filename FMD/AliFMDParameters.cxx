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

//____________________________________________________________________
//                                                                          
// Forward Multiplicity Detector based on Silicon wafers. 
//
// This class is a singleton that handles various parameters of
// the FMD detectors.  
// Eventually, this class will use the Conditions DB to get the
// various parameters, which code can then request from here.
//                                                       
#include "AliLog.h"		   // ALILOG_H
#include "AliFMDParameters.h"	   // ALIFMDPARAMETERS_H
#include "AliFMDGeometry.h"	   // ALIFMDGEOMETRY_H
#include "AliFMDRing.h"	           // ALIFMDRING_H
#include "AliFMDCalibGain.h"       // ALIFMDCALIBGAIN_H
#include "AliFMDCalibPedestal.h"   // ALIFMDCALIBPEDESTAL_H
#include "AliFMDCalibSampleRate.h" // ALIFMDCALIBPEDESTAL_H
#include "AliFMDAltroMapping.h"    // ALIFMDALTROMAPPING_H
#include <AliCDBManager.h>         // ALICDBMANAGER_H
#include <AliCDBEntry.h>           // ALICDBMANAGER_H
#include <Riostream.h>

//====================================================================
ClassImp(AliFMDParameters)
#if 0
  ; // This is here to keep Emacs for indenting the next line
#endif

//____________________________________________________________________
AliFMDParameters* AliFMDParameters::fgInstance = 0;

//____________________________________________________________________
const char* AliFMDParameters::fgkPulseGain	 = "FMD/Calib/PulseGain";
const char* AliFMDParameters::fgkPedestal	 = "FMD/Calib/Pedestal";
const char* AliFMDParameters::fgkDead	         = "FMD/Calib/Dead";
const char* AliFMDParameters::fgkSampleRate	 = "FMD/Calib/SampleRate";
const char* AliFMDParameters::fgkAltroMap	 = "FMD/Calib/AltroMap";
const char* AliFMDParameters::fgkZeroSuppression = "FMD/Calib/ZeroSuppression";


//____________________________________________________________________
AliFMDParameters* 
AliFMDParameters::Instance() 
{
  // Get static instance 
  if (!fgInstance) fgInstance = new AliFMDParameters;
  return fgInstance;
}

//____________________________________________________________________
AliFMDParameters::AliFMDParameters() 
  : fIsInit(kFALSE),
    fSiDeDxMip(1.664), 
    fFixedPulseGain(0), 
    fEdepMip(0),
    fZeroSuppression(0), 
    fSampleRate(0), 
    fPedestal(0), 
    fPulseGain(0), 
    fDeadMap(0), 
    fAltroMap(0)
{
  // Default constructor 
  SetVA1MipRange();
  SetAltroChannelSize();
  SetChannelsPerAltro();
  SetZeroSuppression();
  SetSampleRate();
  SetPedestal();
  SetPedestalWidth();
  SetPedestalFactor();
  SetThreshold();
}

//__________________________________________________________________
void
AliFMDParameters::Init()
{
  // Initialize the parameters manager.  We need to get stuff from the
  // CDB here. 
  if (fIsInit) return;
  InitPulseGain();
  InitPedestal();
  InitDeadMap();
  InitSampleRate();
  InitZeroSuppression();
  InitAltroMap();
  fIsInit = kTRUE;
  
}
//__________________________________________________________________
void
AliFMDParameters::InitPulseGain()
{
  AliCDBManager* cdb      = AliCDBManager::Instance();
  AliCDBEntry*   gain     = cdb->Get(fgkPulseGain);
  if (!gain) {
    AliWarning(Form("No %s found in CDB, perhaps you need to "
		    "use AliFMDCalibFaker?", fgkPulseGain));
    return;
  }
  
  AliDebug(1, Form("Got gain from CDB"));
  fPulseGain = dynamic_cast<AliFMDCalibGain*>(gain->GetObject());
  if (!fPulseGain) AliWarning("Invalid pulser gain object from CDB");
}
//__________________________________________________________________
void
AliFMDParameters::InitPedestal()
{
  AliCDBManager* cdb      = AliCDBManager::Instance();
  AliCDBEntry*   pedestal = cdb->Get(fgkPedestal);
  if (!pedestal) {
    AliWarning(Form("No %s found in CDB, perhaps you need to "
		    "use AliFMDCalibFaker?", fgkPedestal));
    return;
  }
  AliDebug(1, Form("Got pedestal from CDB"));
  fPedestal = dynamic_cast<AliFMDCalibPedestal*>(pedestal->GetObject());
  if (!fPedestal) AliWarning("Invalid pedestal object from CDB");
}

//__________________________________________________________________
void
AliFMDParameters::InitDeadMap()
{
  AliCDBManager* cdb      = AliCDBManager::Instance();
  AliCDBEntry*   deadMap  = cdb->Get(fgkDead);
  if (!deadMap) {
    AliWarning(Form("No %s found in CDB, perhaps you need to "
		    "use AliFMDCalibFaker?", fgkDead));
    return;
  }
  AliDebug(1, Form("Got dead map from CDB"));
  fDeadMap = dynamic_cast<AliFMDCalibDeadMap*>(deadMap->GetObject());
  if (!fDeadMap) AliWarning("Invalid dead map object from CDB");
}

//__________________________________________________________________
void
AliFMDParameters::InitZeroSuppression()
{
  AliCDBManager* cdb      = AliCDBManager::Instance();
  AliCDBEntry*   zeroSup  = cdb->Get(fgkZeroSuppression);
  if (!zeroSup) {
    AliWarning(Form("No %s found in CDB, perhaps you need to "
		    "use AliFMDCalibFaker?", fgkZeroSuppression));
    return;
  }
  AliDebug(1, Form("Got zero suppression from CDB"));
  fZeroSuppression = 
    dynamic_cast<AliFMDCalibZeroSuppression*>(zeroSup->GetObject());
  if (!fZeroSuppression)AliWarning("Invalid zero suppression object from CDB");
}

//__________________________________________________________________
void
AliFMDParameters::InitSampleRate()
{
  AliCDBManager* cdb      = AliCDBManager::Instance();
  AliCDBEntry*   sampRat  = cdb->Get(fgkSampleRate);
  if (!sampRat) {
    AliWarning(Form("No %s found in CDB, perhaps you need to "
		    "use AliFMDCalibFaker?", fgkSampleRate));
    return;
  }
  AliDebug(1, Form("Got zero suppression from CDB"));
  fSampleRate = dynamic_cast<AliFMDCalibSampleRate*>(sampRat->GetObject());
  if (!fSampleRate) AliWarning("Invalid zero suppression object from CDB");
}

//__________________________________________________________________
void
AliFMDParameters::InitAltroMap()
{
  AliCDBManager* cdb      = AliCDBManager::Instance();
  AliCDBEntry*   hwMap    = cdb->Get(fgkAltroMap);       
  if (!hwMap) {
    AliWarning(Form("No %s found in CDB, perhaps you need to "
		    "use AliFMDCalibFaker?", fgkAltroMap));
    fAltroMap = new AliFMDAltroMapping;
    return;
  }
  AliDebug(1, Form("Got ALTRO map from CDB"));
  fAltroMap = dynamic_cast<AliFMDAltroMapping*>(hwMap->GetObject());
  if (!fAltroMap) {
    AliWarning("Invalid ALTRO map object from CDB");
    fAltroMap = new AliFMDAltroMapping;
  }
}


//__________________________________________________________________
Float_t
AliFMDParameters::GetThreshold() const
{
  if (!fPulseGain) return fFixedThreshold;
  return fPulseGain->Threshold();
}

//__________________________________________________________________
Float_t
AliFMDParameters::GetPulseGain(UShort_t detector, Char_t ring, 
			       UShort_t sector, UShort_t strip) const
{
  // Returns the pulser calibrated gain for strip # strip in sector #
  // sector or ring id ring of detector # detector. 
  // 
  // For simulation, this is normally set to 
  // 
  //       VA1_MIP_Range 
  //    ------------------ * MIP_Energy_Loss
  //    ALTRO_channel_size
  // 
  if (!fPulseGain) { 
    if (fFixedPulseGain <= 0)
      fFixedPulseGain = fVA1MipRange * GetEdepMip() / fAltroChannelSize;
    return fFixedPulseGain;
  }  
  AliDebug(50, Form("pulse gain for FMD%d%c[%2d,%3d]=%f",
		    detector, ring, sector, strip,
		    fPulseGain->Value(detector, ring, sector, strip)));
  return fPulseGain->Value(detector, ring, sector, strip);
}

//__________________________________________________________________
Bool_t
AliFMDParameters::IsDead(UShort_t detector, Char_t ring, 
			 UShort_t sector, UShort_t strip) const
{
  if (!fDeadMap) return kFALSE;
  AliDebug(50, Form("Dead for FMD%d%c[%2d,%3d]=%s",
		    detector, ring, sector, strip,
		    fDeadMap->operator()(detector, ring, sector, strip) ? 
		    "no" : "yes"));
  return fDeadMap->operator()(detector, ring, sector, strip);
}

//__________________________________________________________________
UShort_t
AliFMDParameters::GetZeroSuppression(UShort_t detector, Char_t ring, 
				     UShort_t sector, UShort_t strip) const
{
  if (!fZeroSuppression) return fFixedZeroSuppression;
  // Need to map strip to ALTRO chip. 
  AliDebug(50, Form("zero sup. for FMD%d%c[%2d,%3d]=%f",
		    detector, ring, sector, strip,
		    fZeroSuppression->operator()(detector, ring, 
						 sector, strip)));
  return fZeroSuppression->operator()(detector, ring, sector, strip/128);
}

//__________________________________________________________________
UShort_t
AliFMDParameters::GetSampleRate(UShort_t ddl) const
{
  if (!fSampleRate) return fFixedSampleRate;
  // Need to map sector to digitizier card. 
  AliDebug(50, Form("Sample rate for %d=%d", ddl, fSampleRate->Rate(ddl)));
  return fSampleRate->Rate(ddl);
}

//__________________________________________________________________
Float_t
AliFMDParameters::GetPedestal(UShort_t detector, Char_t ring, 
			      UShort_t sector, UShort_t strip) const
{
  if (!fPedestal) return fFixedPedestal;
  AliDebug(50, Form("pedestal for FMD%d%c[%2d,%3d]=%f",
		    detector, ring, sector, strip,
		    fPedestal->Value(detector, ring, sector, strip)));
  return fPedestal->Value(detector, ring, sector, strip);
}

//__________________________________________________________________
Float_t
AliFMDParameters::GetPedestalWidth(UShort_t detector, Char_t ring, 
				   UShort_t sector, UShort_t strip) const
{
  if (!fPedestal) return fFixedPedestalWidth;
  AliDebug(50, Form("pedetal width for FMD%d%c[%2d,%3d]=%f",
		    detector, ring, sector, strip,
		    fPedestal->Width(detector, ring, sector, strip)));
  return fPedestal->Width(detector, ring, sector, strip);
}
  
//__________________________________________________________________
AliFMDAltroMapping*
AliFMDParameters::GetAltroMap() const
{
  return fAltroMap;
}


//__________________________________________________________________
Bool_t
AliFMDParameters::Hardware2Detector(UInt_t ddl, UInt_t addr, UShort_t& det,
				    Char_t& ring, UShort_t& sec, 
				    UShort_t& str) const
{
  if (!fAltroMap) return kFALSE;
  return fAltroMap->Hardware2Detector(ddl, addr, det, ring, sec, str);
}

//__________________________________________________________________
Bool_t
AliFMDParameters::Detector2Hardware(UShort_t det, Char_t ring, UShort_t sec, 
				    UShort_t str, UInt_t& ddl, 
				    UInt_t& addr) const			      
{
  if (!fAltroMap) return kFALSE;
  return fAltroMap->Detector2Hardware(det, ring, sec, str, ddl, addr);
}


//__________________________________________________________________
Float_t
AliFMDParameters::GetEdepMip() const 
{ 
  // Get energy deposited by a MIP in the silicon sensors
  if (fEdepMip <= 0){
    AliFMDGeometry* fmd = AliFMDGeometry::Instance();
    fEdepMip = (fSiDeDxMip 
		* fmd->GetRing('I')->GetSiThickness() 
		* fmd->GetSiDensity());
  }
  return fEdepMip;
}


  
  
  
//____________________________________________________________________
//
// EOF
//
