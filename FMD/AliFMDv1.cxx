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
// Forward Multiplicity Detector based on Silicon wafers. This class
// contains the base procedures for the Forward Multiplicity detector
// Detector consists of 3 sub-detectors FMD1, FMD2, and FMD3, each of
// which has 1 or 2 rings of silicon sensors. 
//                                                       
// This class contains the detailed version of the FMD - that is, hits
// are produced during simulation. 
//                                                                           
// See also the class AliFMD for a more detailed explanation of the
// various componets. 
//
#include <TVirtualMC.h>		// ROOT_TVirtualMC
#include <AliRun.h>		// ALIRUN_H
#include <AliMC.h>		// ALIMC_H
#include <AliLog.h>		// ALILOG_H
#include "AliFMDv1.h"		// ALIFMDV1_H
#include "AliFMDSimulator.h"	// ALIFMDSIMULATOR_H
#include "AliFMDG3Simulator.h"	// ALIFMDG3SIMULATOR_H
#include "AliFMDGeoSimulator.h"	// ALIFMDGEOSIMULATOR_H

//____________________________________________________________________
ClassImp(AliFMDv1)
#if 0
  ; // This is here to keep Emacs for indenting the next line
#endif


//____________________________________________________________________
void 
AliFMDv1::StepManager()
{
  // Called for every step in the Forward Multiplicity Detector
  //
  // The message is deligated to AliFMDSimulator::Exec 
  // 
  if (!fSimulator) {
    AliFatal("No simulator object made");
    return;
  }
  fSimulator->Exec("");
}
//___________________________________________________________________
//
// EOF
//
