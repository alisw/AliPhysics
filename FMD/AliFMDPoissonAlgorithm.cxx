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
// Base class for FMD poisson algorithms. 
//
// Derived classes will implement various ways of reconstructing the
// charge particle multiplicity in the FMD.  
// 
#include "AliFMDPoissonAlgorithm.h"	// ALIFMDPOISSONALGORITHM_H
#include "AliFMDDigit.h"		// ALIFMDDIGIT_H

//____________________________________________________________________
ClassImp(AliFMDPoissonAlgorithm);

//____________________________________________________________________
AliFMDPoissonAlgorithm::AliFMDPoissonAlgorithm()
  : AliFMDReconstructionAlgorithm("Poisson", "Poisson")
{}


//____________________________________________________________________
void
AliFMDPoissonAlgorithm::Reset() 
{
  // Reset internal data 
  fTotal = 0;
  fEmpty.Reset(kFALSE);
}

//____________________________________________________________________
void
AliFMDPoissonAlgorithm::ProcessDigit(AliFMDDigit* digit, 
				     Float_t eta, Float_t phi, 
				     UShort_t count)
{
  // Process one digit. 
  // 
  // Parameters: 
  //    
  //   digit		Digit to process 
  //   ipZ		Z--coordinate of the primary interaction
  //                    vertex of this event 
  //
  if (!digit) return;
  fTotal++;
  if (charge < threshold) fEmpty(digit->Detector() - 1, 
				 digit->Ring(), 
				 digit->Setctor(), 
				 digit->Strip()) = kTRUE;
}


//____________________________________________________________________
// 
// EOF
//
