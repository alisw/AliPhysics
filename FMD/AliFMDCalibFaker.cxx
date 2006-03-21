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
// This task creates fake calibrations. Which calibration, depends on
// the bit mask passed to the constructor, or added by `AddCalib'.
//
// The default is to write all calibration parameters to a local
// storage `local://cdb' which is a directory in the current
// directory. 
//                                                       
#include "AliLog.h"		   // ALILOG_H
#include "AliFMDCalibFaker.h"      // ALIFMDCALIBFAKER_H
#include "AliFMDCalibGain.h"       // ALIFMDCALIBGAIN_H
#include "AliFMDCalibPedestal.h"   // ALIFMDCALIBPEDESTAL_H
#include "AliFMDCalibSampleRate.h" // ALIFMDCALIBPEDESTAL_H
#include "AliFMDAltroMapping.h"    // ALIFMDALTROMAPPING_H
#include <AliCDBManager.h>         // ALICDBMANAGER_H
#include <AliCDBEntry.h>           // ALICDBMANAGER_H
#include <Riostream.h>
#include <TSystem.h>
#include <TMath.h>
#include <TROOT.h>
#include <TRandom.h>

//====================================================================
ClassImp(AliFMDCalibFaker)
#if 0
  ; // This is here to keep Emacs for indenting the next line
#endif

//____________________________________________________________________
AliFMDCalibFaker::AliFMDCalibFaker(Int_t mask, const char* loc) 
  : TTask("FMDCalibFaker", loc),
    fMask(mask),
    fGain(-1),
    fThresholdFactor(.1),
    fThreshold(-1),
    fPedestalMin(20),
    fPedestalMax(30), 
    fDeadChance(0),
    fRate(1),
    fZeroThreshold(0),
    fRunMin(0),
    fRunMax(10)
{
  // Default constructor 
}


#define MAKE_META(meta) \
  do { \
    meta = new AliCDBMetaData; \
    meta->SetResponsible(gSystem->GetUserInfo()->fRealName.Data()); \
    meta->SetAliRootVersion(gROOT->GetVersion()); \
    meta->SetBeamPeriod(1); \
    meta->SetComment("Dummy data for testing"); } while (false);
 

//__________________________________________________________________
void
AliFMDCalibFaker::Exec(Option_t*)
{
  // Make the objects. 
  AliCDBManager*     cdb      = AliCDBManager::Instance();
  AliFMDParameters*  param    = AliFMDParameters::Instance();
  Float_t            maxADC   = param->GetAltroChannelSize();
  TObjArray          cleanup;

  if (GetTitle()) cdb->SetDefaultStorage(GetTitle());
    
  AliCDBMetaData* meta = 0;
  if (TESTBIT(fMask, kPulseGain)) {
    if (fGain <= 0) 
      fGain      = (param->GetVA1MipRange() * param->GetEdepMip() / maxADC);
    fThreshold = fThresholdFactor * param->GetEdepMip();
    AliFMDCalibGain* gain = MakePulseGain();
    AliCDBId         id(AliFMDParameters::fgkPulseGain, fRunMin, fRunMax);
    MAKE_META(meta);
    meta->SetProperty("key1", gain);
    cdb->Put(gain, id, meta);
    cleanup.Add(gain);
    cleanup.Add(meta);
  }  
  if (TESTBIT(fMask, kPedestal)) {
    fPedestalMin = TMath::Max(TMath::Min(fPedestalMin, maxADC), 0.F);
    fPedestalMax = TMath::Max(TMath::Min(fPedestalMax, maxADC), fPedestalMin);
    AliFMDCalibPedestal* pedestal = MakePedestal();
    AliCDBId             id(AliFMDParameters::fgkPedestal, fRunMin, fRunMax);
    MAKE_META(meta);
    meta->SetProperty("key1", pedestal);
    cdb->Put(pedestal, id, meta);
    cleanup.Add(pedestal);
    cleanup.Add(meta);
  }
  if (TESTBIT(fMask, kDeadMap)) {
    fDeadChance = TMath::Max(TMath::Min(fDeadChance, 1.F), 0.F);
    AliFMDCalibDeadMap* deadMap = MakeDeadMap();
    AliCDBId            id(AliFMDParameters::fgkDead, fRunMin, fRunMax);
    MAKE_META(meta);
    meta->SetProperty("key1", deadMap);
    cdb->Put(deadMap, id, meta);
    cleanup.Add(deadMap);
    cleanup.Add(meta);
  }
  if (TESTBIT(fMask, kZeroSuppression)) {
    fZeroThreshold = TMath::Min(fZeroThreshold, UShort_t(maxADC));
    AliFMDCalibZeroSuppression* zeroSup = MakeZeroSuppression();
    AliCDBId                    id(AliFMDParameters::fgkZeroSuppression, 
				   fRunMin, fRunMax);
    MAKE_META(meta);
    meta->SetProperty("key1", zeroSup);
    cdb->Put(zeroSup, id, meta);
    cleanup.Add(zeroSup);
    cleanup.Add(meta);
  }
  if (TESTBIT(fMask, kSampleRate)) {
    fRate = TMath::Max(TMath::Min(fRate, UShort_t(8)), UShort_t(1));
    AliFMDCalibSampleRate* rate = MakeSampleRate();
    AliCDBId               id(AliFMDParameters::fgkSampleRate,fRunMin,fRunMax);
    MAKE_META(meta);
    meta->SetProperty("key1", rate);
    cdb->Put(rate, id, meta);
    cleanup.Add(rate);
    cleanup.Add(meta);
  }
  if (TESTBIT(fMask, kAltroMap)) {
    AliFMDAltroMapping* altroMap = MakeAltroMap();
    AliCDBId            id(AliFMDParameters::fgkAltroMap, fRunMin, fRunMax);
    MAKE_META(meta);
    meta->SetProperty("key1", altroMap);
    cdb->Put(altroMap, id, meta);
    cleanup.Add(altroMap);
    cleanup.Add(meta);
  }
  cdb->Destroy();
  cleanup.Delete();
}


//__________________________________________________________________
AliFMDCalibGain*
AliFMDCalibFaker::MakePulseGain()
{
  // Make the actual data
  AliFMDCalibGain*  gain  = new AliFMDCalibGain;
  // Set threshold 
  gain->Set(fThreshold);
  for (UShort_t det = 1; det <= 3; det++) {
    Char_t rings[] = { 'I', (det == 1 ? '\0' : 'O'), '\0' };
    for (Char_t* ring = rings; *ring != '\0'; ring++) {
      UShort_t nSec = ( *ring == 'I' ? 20  :  40 );
      UShort_t nStr = ( *ring == 'I' ? 512 : 256 );
      for (UShort_t sec = 0; sec < nSec; sec++) {
        for (UShort_t str = 0; str < nStr; str++) {
          gain->Set(det, *ring, sec, str,
                    gRandom->Gaus(fGain, .01 * fGain));
        }
      }
    }
  }
  return gain;
}

//__________________________________________________________________
AliFMDCalibPedestal*
AliFMDCalibFaker::MakePedestal()
{
  // Make the actual data
  AliFMDCalibPedestal*  pedestal  = new AliFMDCalibPedestal;
  for (UShort_t det = 1; det <= 3; det++) {
    Char_t rings[] = { 'I', (det == 1 ? '\0' : 'O'), '\0' };
    for (Char_t* ring = rings; *ring != '\0'; ring++) {
      UShort_t nSec = ( *ring == 'I' ? 20  :  40 );
      UShort_t nStr = ( *ring == 'I' ? 512 : 256 );
      for (UShort_t sec = 0; sec < nSec; sec++) {
        for (UShort_t str = 0; str < nStr; str++) {
          pedestal->Set(det, *ring, sec, str,
			gRandom->Uniform(fPedestalMin, fPedestalMax), 1.5);
        }
      }
    }
  }
  return pedestal;
}

//__________________________________________________________________
AliFMDCalibDeadMap*
AliFMDCalibFaker::MakeDeadMap()
{
  // Make the actual data
  AliFMDCalibDeadMap*  deadmap  = new AliFMDCalibDeadMap;
  for (UShort_t det = 1; det <= 3; det++) {
    Char_t rings[] = { 'I', (det == 1 ? '\0' : 'O'), '\0' };
    for (Char_t* ring = rings; *ring != '\0'; ring++) {
      UShort_t nSec = ( *ring == 'I' ? 20  :  40 );
      UShort_t nStr = ( *ring == 'I' ? 512 : 256 );
      for (UShort_t sec = 0; sec < nSec; sec++) {
        for (UShort_t str = 0; str < nStr; str++) {
          deadmap->operator()(det, *ring, sec, str) = 
	    gRandom->Uniform(0, 1) < fDeadChance;
        }
      }
    }
  }
  return deadmap;
}

//__________________________________________________________________
AliFMDCalibZeroSuppression*
AliFMDCalibFaker::MakeZeroSuppression()
{
  // Make the actual data
  AliFMDCalibZeroSuppression*  zs  = new AliFMDCalibZeroSuppression;
  for (UShort_t det = 1; det <= 3; det++) {
    Char_t rings[] = { 'I', (det == 1 ? '\0' : 'O'), '\0' };
    for (Char_t* ring = rings; *ring != '\0'; ring++) {
      UShort_t nSec = ( *ring == 'I' ? 20  :  40 );
      UShort_t nStr = ( *ring == 'I' ? 512 : 256 );
      for (UShort_t sec = 0; sec < nSec; sec++) {
        for (UShort_t str = 0; str < nStr; str++) {
          zs->operator()(det, *ring, sec, str) =  fZeroThreshold;
        }
      }
    }
  }
  return zs;
}

//__________________________________________________________________
AliFMDCalibSampleRate*
AliFMDCalibFaker::MakeSampleRate()
{
  AliFMDCalibSampleRate*  sampleRate  = new AliFMDCalibSampleRate;
  for (int i = 0; i < 3; i++)
    sampleRate->Set(AliFMDParameters::kBaseDDL+i, fRate);
  return sampleRate;
}

//__________________________________________________________________
AliFMDAltroMapping*
AliFMDCalibFaker::MakeAltroMap()
{
  AliFMDAltroMapping*  m  = new AliFMDAltroMapping;
  return m;
}
  
  
  
//____________________________________________________________________
//
// EOF
//
