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
// Base class for FMD reconstruction algorithms. 
//
// Derived classes will implement various ways of reconstructing the
// charge particle multiplicity in the FMD.  
// 
#include "AliFMDReconstructionAlgorithm.h"	// ALIFMDRECONSTRUCTIONALGORITHM_H
#include "AliFMDDigit.h"			// ALIFMDDIGIT_H

//____________________________________________________________________
ClassImp(AliFMDReconstructionAlgorithm);

//____________________________________________________________________
AliFMDReconstructionAlgorithm::AliFMDReconstructionAlgorithm(const char* name, 
							     const char* title)
  : TNamed(name, title)
{}



//____________________________________________________________________
// 
// EOF
//
