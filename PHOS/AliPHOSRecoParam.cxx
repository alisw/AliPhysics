/**************************************************************************
 * Copyright(c) 2007, ALICE Experiment at CERN, All rights reserved.      *
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

// Base class for the PHOS reconstruction parameters.
// Do not use in the reconstruction; use derivative classes instead.
// Author: Boris Polichtchouk.

// --- AliRoot header files ---
#include "AliCDBManager.h"
#include "AliCDBEntry.h"
#include "AliLog.h"
#include "AliPHOSRecoParam.h"

ClassImp(AliPHOSRecoParam)

TObjArray* AliPHOSRecoParam::fgkMaps =0; //ALTRO mappings

//-----------------------------------------------------------------------------
AliPHOSRecoParam::AliPHOSRecoParam() :
  AliDetectorRecoParam(),
  fEMCClusteringThreshold(0.2),
  fEMCLocMaxCut(0.03),
  fEMCRawDigitThreshold(2),
  fEMCMinE(0.012),
  fEMCW0(4.5),
  fEMCSampleQualityCut(1.),
  fEMCEcoreRadius(3.),
  fEMCEcore2ESD(kFALSE),
  fEMCSubtractPedestals(kTRUE),
  fEMCUnfold(kTRUE),
  fEMCEnergyCorrectionOn(kTRUE),
  fEMCDecoderVersion(""),
  fGlobalAltroOffset(0),
  fCPVClusteringThreshold(0.0),
  fCPVLocMaxCut(0.03),
  fCPVMinE(0.0),
  fCPVW0(4.0),
  fCPVUnfold(kTRUE)
{
  //Default constructor.
}

//-----------------------------------------------------------------------------
AliPHOSRecoParam::AliPHOSRecoParam(const AliPHOSRecoParam& ):
  AliDetectorRecoParam(),
  fEMCClusteringThreshold(0.2),
  fEMCLocMaxCut(0.03),
  fEMCRawDigitThreshold(2),
  fEMCMinE(0.012),
  fEMCW0(4.5),
  fEMCSampleQualityCut(1.),
  fEMCEcoreRadius(3.),
  fEMCEcore2ESD(kFALSE),
  fEMCSubtractPedestals(kTRUE),
  fEMCUnfold(kTRUE),
  fEMCEnergyCorrectionOn(kTRUE),
  fEMCDecoderVersion(""),
  fGlobalAltroOffset(0),
  fCPVClusteringThreshold(0.0),
  fCPVLocMaxCut(0.03),
  fCPVMinE(0.0),
  fCPVW0(4.0),
  fCPVUnfold(kTRUE)
{
  //Copy constructor.
}

//-----------------------------------------------------------------------------
AliPHOSRecoParam& AliPHOSRecoParam::operator = (const AliPHOSRecoParam& recoParam)
{
  //Assignment operator.

  if(this != &recoParam) {
    fEMCClusteringThreshold = recoParam.fEMCClusteringThreshold;
    fEMCLocMaxCut           = recoParam.fEMCLocMaxCut;
    fEMCRawDigitThreshold   = recoParam.fEMCRawDigitThreshold;
    fEMCMinE                = recoParam.fEMCMinE;
    fEMCW0                  = recoParam.fEMCW0;
    fEMCSampleQualityCut    = recoParam.fEMCSampleQualityCut;
    fEMCEcoreRadius         = recoParam.fEMCEcoreRadius;
    fEMCEcore2ESD           = recoParam.fEMCEcore2ESD;
    fEMCSubtractPedestals   = recoParam.fEMCSubtractPedestals;
    fEMCUnfold              = recoParam.fEMCUnfold;
    fEMCEnergyCorrectionOn  = recoParam.fEMCEnergyCorrectionOn;
    fEMCDecoderVersion      = recoParam.fEMCDecoderVersion;
    fGlobalAltroOffset      = recoParam.fGlobalAltroOffset;
    fCPVClusteringThreshold = recoParam.fCPVClusteringThreshold;
    fCPVLocMaxCut           = recoParam.fCPVLocMaxCut;
    fCPVMinE                = recoParam.fCPVMinE;
    fCPVW0                  = recoParam.fCPVW0;
    fCPVUnfold              = recoParam.fCPVUnfold;
  }

  return *this;
}

//-----------------------------------------------------------------------------
void AliPHOSRecoParam::Print(Option_t * /*option*/) const
{
  AliDebug(2,Form("PHOS reconstruction parameters:\n"
		  "\tEMCClusteringThreshold = %f\n"
		  "\tEMCLocMaxCut           = %f\n"
		  "\tEMCRawDigitThreshold   = %f\n"
		  "\tEMCMinE                = %f\n"
		  "\tEMCW0                  = %f\n"
		  "\tEMCSampleQualityCut    = %f\n"
		  "\tEMCEcoreRadius         = %f\n"
		  "\tEMCEcore2ESD           = %d\n"
		  "\tEMCSubtractPedestals   = %d\n"
		  "\tEMCUnfold              = %d\n"
		  "\tEMCEnergyCorrectionOn  = %d\n"
		  "\tEMCDecoderVersion      = \"%s\"\n"
		  "\tGlobalAltroOffset      = %d",
		  fEMCClusteringThreshold,
		  fEMCLocMaxCut,
		  fEMCRawDigitThreshold,
		  fEMCMinE,
		  fEMCW0,
		  fEMCSampleQualityCut,
		  fEMCEcoreRadius,
		  fEMCEcore2ESD,
		  fEMCSubtractPedestals,
		  fEMCUnfold,
		  fEMCEnergyCorrectionOn,
		  fEMCDecoderVersion.Data(),
		  fGlobalAltroOffset));

}

//-----------------------------------------------------------------------------
AliPHOSRecoParam* AliPHOSRecoParam::GetDefaultParameters()
{
  //Default parameters for the reconstruction

  AliPHOSRecoParam* params = new AliPHOSRecoParam();
  return params;
}

//-----------------------------------------------------------------------------
const TObjArray* AliPHOSRecoParam::GetMappings()
{
  //Returns array of AliAltroMappings for RCU0..RCU3.
  //If not found, read it from OCDB.

  //Quick check as follows:
  //  root [0] AliCDBManager::Instance()->SetDefaultStorage("local://$ALICE_ROOT");
  //  root [1] AliCDBManager::Instance()->SetRun(1);
  //  root [2] TObjArray* maps = AliPHOSRecoParam::GetMappings();
  //  root [3] maps->Print();
  
  if(fgkMaps) return fgkMaps;
  
  AliCDBEntry* entry = AliCDBManager::Instance()->Get("PHOS/Calib/Mapping");
  if(entry)
    fgkMaps = (TObjArray*)entry->GetObject();
  
  return fgkMaps;
  
}
