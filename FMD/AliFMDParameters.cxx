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
//                                                       
#include "AliFMDParameters.h"	// ALIFMDPARAMETERS_H
#include "AliFMDGeometry.h"	// ALIFMDGEOMETRY_H
#include "AliFMDRing.h"	        // ALIFMDRING_H
#include "AliLog.h"		// ALILOG_H
#include <Riostream.h>

//====================================================================
ClassImp(AliFMDParameters)
#if 0
  ; // This is here to keep Emacs for indenting the next line
#endif

//____________________________________________________________________
AliFMDParameters* AliFMDParameters::fgInstance = 0;

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
  : fSiDeDxMip(1.664)
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
}



//__________________________________________________________________
Float_t
AliFMDParameters::GetEdepMip() const 
{ 
  // Get energy deposited by a MIP in the silicon sensors
  AliFMDGeometry* fmd = AliFMDGeometry::Instance();
  return (fSiDeDxMip 
	  * fmd->GetRing('I')->GetSiThickness() 
	  * fmd->GetSiDensity());
}

//____________________________________________________________________
//
// EOF
//
