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
///////////////////////////////////////////////////////////////////////////////
//
// FMD Reconstruction Parameters 
//
//
///////////////////////////////////////////////////////////////////////////////


#include "AliFMDRecoParam.h"

ClassImp(AliDetectorRecoParam)
#if 0 
; // Don't delete - for Emacs
#endif

//____________________________________________________________________
AliFMDRecoParam::AliFMDRecoParam(Float_t noiseFactor, 
				 Bool_t angleCorrect,
				 Bool_t sharingCorrect)
  : AliDetectorRecoParam(),
    fNoiseFactor(noiseFactor), 
    fAngleCorrect(angleCorrect), 
    fSharingCorrect(sharingCorrect)
{
  // Constructor
  SetName("FMD");
  SetTitle("FMD");
}

//____________________________________________________________________
AliFMDRecoParam*
AliFMDRecoParam::GetLowFluxParam()
{
  // 
  // Get low flux parameter
  //
  // Return:
  //    low flux parameters 
  //  
  AliFMDRecoParam* p = new AliFMDRecoParam(10, kTRUE, kFALSE);
  p->SetName("FMD_low_flux");
  p->SetTitle("FMD low flux");
  return p;
}
//____________________________________________________________________
AliFMDRecoParam*
AliFMDRecoParam::GetHighFluxParam()
{
  // 
  // Get high flux parameter
  //
  // Return:
  //    high flux parameters 
  //  
  AliFMDRecoParam* p = new AliFMDRecoParam(10, kTRUE, kFALSE);
  p->SetName("FMD_high_flux");
  p->SetTitle("FMD high flux");
  return p;
}
//____________________________________________________________________
AliFMDRecoParam*
AliFMDRecoParam::GetParam(AliRecoParam::EventSpecie_t specie)
{
  // 
  // Get parameters for a specific species 
  // 
  // Parameters:
  //    specie Species 
  // 
  // Return:
  //    Reconstruction paramters 
  //
  switch (specie) { 
  case AliRecoParam::kDefault: 
  case AliRecoParam::kCalib: 
  case AliRecoParam::kHighMult: return GetHighFluxParam(); 
  case AliRecoParam::kCosmic: 
  case AliRecoParam::kLowMult:  return GetLowFluxParam(); 
  }
  return new AliFMDRecoParam();
}


//____________________________________________________________________
//
//
// EOF
//
