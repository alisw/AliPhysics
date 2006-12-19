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
/** @file    AliFMDCalibFaker.cxx
    @author  Christian Holm Christensen <cholm@nbi.dk>
    @date    Sun Mar 26 18:29:21 2006
    @brief   Make fake calibration data 
    @ingroup FMD_util
*/
//____________________________________________________________________
//                                                                          
// Forward Multiplicity Detector based on Silicon wafers. 
//
// This task creates fake calibrations. Which calibration, depends on
// the bit mask passed to the constructor, or added by `AddCalib'.
//
// The default is to write all calibration parameters to a local
// storage `local://$ALICE_ROOT' which is were the sources live (sigh!
// - why oh why do we need to shit where we eat - it's just not
// healty).
//                                                       
#include "AliLog.h"		   // ALILOG_H
#include "AliFMDCalibFaker.h"      // ALIFMDCALIBFAKER_H
#include "AliFMDCalibGain.h"       // ALIFMDCALIBGAIN_H
#include "AliFMDCalibPedestal.h"   // ALIFMDCALIBPEDESTAL_H
#include "AliFMDCalibSampleRate.h" // ALIFMDCALIBPEDESTAL_H
#include "AliFMDAltroMapping.h"    // ALIFMDALTROMAPPING_H
#include "AliFMDCalibStripRange.h" // ALIFMDCALIBSTRIPRANGE_H
#include <AliCDBManager.h>         // ALICDBMANAGER_H
#include <AliCDBEntry.h>           // ALICDBMANAGER_H
//#include <Riostream.h>
#include <TSystem.h>
// #include <TMath.h>
#include <TROOT.h>
#include <TRandom.h>
#include <TF1.h>

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
    fRunMax(10),
    fStripMin(0),
    fStripMax(127)
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
  Float_t            maxADC   = 1.F*param->GetAltroChannelSize();
  TObjArray          cleanup;

  if (GetTitle() && GetTitle()[0] != '\0') { 
    AliInfo(Form("Setting default storage to '%s'", GetTitle()));
    cdb->SetDefaultStorage(GetTitle());
  }
  
    
  AliCDBMetaData* meta = 0;
  if (TESTBIT(fMask, kPulseGain)) {
    // Info("Exec","Default gain to %f = %d * %f / %d", 
    //      (param->GetVA1MipRange() * param->GetEdepMip() / maxADC),
    //      param->GetVA1MipRange(), param->GetEdepMip(), Int_t(maxADC));
    if (fGain <= 0) {
      fGain = (param->GetVA1MipRange() * param->GetEdepMip() / maxADC);
    }
    fThreshold = fThresholdFactor * param->GetEdepMip();
    AliFMDCalibGain* gain = MakePulseGain();
    AliCDBId         id(AliFMDParameters::PulseGainPath(), fRunMin, fRunMax);
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
    AliCDBId             id(AliFMDParameters::PedestalPath(),fRunMin,fRunMax);
    MAKE_META(meta);
    meta->SetProperty("key1", pedestal);
    cdb->Put(pedestal, id, meta);
    cleanup.Add(pedestal);
    cleanup.Add(meta);
  }
  if (TESTBIT(fMask, kDeadMap)) {
    fDeadChance = TMath::Max(TMath::Min(fDeadChance, 1.F), 0.F);
    AliFMDCalibDeadMap* deadMap = MakeDeadMap();
    AliCDBId            id(AliFMDParameters::DeadPath(), fRunMin, fRunMax);
    MAKE_META(meta);
    meta->SetProperty("key1", deadMap);
    cdb->Put(deadMap, id, meta);
    cleanup.Add(deadMap);
    cleanup.Add(meta);
  }
  if (TESTBIT(fMask, kZeroSuppression)) {
    fZeroThreshold = TMath::Min(fZeroThreshold, UShort_t(maxADC));
    AliFMDCalibZeroSuppression* zeroSup = MakeZeroSuppression();
    AliCDBId                    id(AliFMDParameters::ZeroSuppressionPath(), 
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
    AliCDBId               id(AliFMDParameters::SampleRatePath(),
			      fRunMin,fRunMax);
    MAKE_META(meta);
    meta->SetProperty("key1", rate);
    cdb->Put(rate, id, meta);
    cleanup.Add(rate);
    cleanup.Add(meta);
  }
  if (TESTBIT(fMask, kStripRange)) {
    fRate = TMath::Max(TMath::Min(fRate, UShort_t(8)), UShort_t(1));
    AliFMDCalibStripRange* range = MakeStripRange();
    AliCDBId               id(AliFMDParameters::StripRangePath(),
			      fRunMin,fRunMax);
    MAKE_META(meta);
    meta->SetProperty("key1", range);
    cdb->Put(range, id, meta);
    cleanup.Add(range);
    cleanup.Add(meta);
  }
  if (TESTBIT(fMask, kAltroMap)) {
    AliFMDAltroMapping* altroMap = MakeAltroMap();
    AliCDBId            id(AliFMDParameters::AltroMapPath(), fRunMin, fRunMax);
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
AliFMDCalibFaker::MakePulseGain() const
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
Float_t 
AliFMDCalibFaker::MakeNoise(Char_t ring, UShort_t str) const
{
  const UShort_t innerN    = 512;
  const UShort_t outerN    = 256;
  const UShort_t innerCut  = 350;
  const UShort_t outerCut  = 190;
  const Float_t  innerBase =   1.2;
  const Float_t  outerBase =   2.1;
  const Float_t  innerInc  =   0.5;
  const Float_t  outerInc  =   0.8;
  Float_t cut, base, inc, n;
  switch (ring) {
  case 'I': 
    cut = innerCut; base = innerBase; inc = innerInc; n = innerN; break;
  case 'O': 
    cut = outerCut; base = outerBase; inc = outerInc; n = outerN; break;
  default:
    return -1;
  }
  Float_t bare = base + (str < cut ? 
			 str / cut * inc : 
			 inc  - (str - cut) / (n - cut) * inc);
  return bare + gRandom->Uniform(-.07, .07);
}
  
//__________________________________________________________________
AliFMDCalibPedestal*
AliFMDCalibFaker::MakePedestal() const
{
  // Make the actual data
  AliFMDCalibPedestal*  pedestal  = new AliFMDCalibPedestal;
  
  for (UShort_t det = 1; det <= 3; det++) {
    Char_t rings[] = { 'I', 'O', '\0' };
    for (Char_t* ring = rings; *ring != '\0'; ring++) {
      if (*ring == 'O' && det == 1) continue;
      UShort_t nSec = ( *ring == 'I' ? 20  :  40 );
      UShort_t nStr = ( *ring == 'I' ? 512 : 256 );
      for (UShort_t sec = 0; sec < nSec; sec++) {
        for (UShort_t str = 0; str < nStr; str++) {
	  Float_t noise = MakeNoise(*ring, str);
	  Float_t ped   = gRandom->Uniform(fPedestalMin, fPedestalMax);
          pedestal->Set(det, *ring, sec, str, ped, noise);
	}
      }
    }
  }
  return pedestal;
}

//__________________________________________________________________
AliFMDCalibDeadMap*
AliFMDCalibFaker::MakeDeadMap() const
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
AliFMDCalibFaker::MakeZeroSuppression() const
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
AliFMDCalibFaker::MakeSampleRate() const
{
  // Make sample rates 
  AliFMDCalibSampleRate*  sampleRate  = new AliFMDCalibSampleRate;
  for (UShort_t det = 1; det <= 3; det++) {
    Char_t rings[] = { 'I', (det == 1 ? '\0' : 'O'), '\0' };
    for (Char_t* ring = rings; *ring != '\0'; ring++) {
      UShort_t nSec = ( *ring == 'I' ? 20  :  40 );
      for (UShort_t sec = 0; sec < nSec; sec++) {
	sampleRate->Set(det, *ring, sec, 0, fRate);
      }
    }
  }
  return sampleRate;
}

//__________________________________________________________________
AliFMDCalibStripRange*
AliFMDCalibFaker::MakeStripRange() const
{
  // Make strip ranges 
  AliFMDCalibStripRange*  striprange  = new AliFMDCalibStripRange;
  for (UShort_t det = 1; det <= 3; det++) {
    Char_t rings[] = { 'I', (det == 1 ? '\0' : 'O'), '\0' };
    for (Char_t* ring = rings; *ring != '\0'; ring++) {
      UShort_t nSec = ( *ring == 'I' ? 20  :  40 );
      for (UShort_t sec = 0; sec < nSec; sec++) {
	striprange->Set(det, *ring, sec, 0, fStripMin, fStripMax);
      }
    }
  }
  return striprange;
}

//__________________________________________________________________
AliFMDAltroMapping*
AliFMDCalibFaker::MakeAltroMap() const
{
  // Make hardware mapping 
  AliFMDAltroMapping*  m  = new AliFMDAltroMapping;
  return m;
}
  
  
  
//____________________________________________________________________
//
// EOF
//
