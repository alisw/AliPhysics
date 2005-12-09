/**************************************************************************
 * Copyright(c) 2004, ALICE Experiment at CERN, All rights reserved. *
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
// Concrete implementation of AliFMDDetector 
//
// This implements the geometry for FMD1 
// FMD1 has only one ring, of type `inner'.  
// It is sitting at z=320.
// It is the FMD ring with highest eta
//
#include "AliFMD1.h"		// ALIFMD1_H 
#include "AliFMDRing.h"		// ALIFMDRING_H 


//====================================================================
ClassImp(AliFMD1)
#if 0
  ; // This is to keep Emacs from indenting the next line 
#endif 

//____________________________________________________________________
AliFMD1::AliFMD1(AliFMDRing* inner) 
  : AliFMDDetector(1, inner, 0)
{
  SetInnerZ(320);
}


//____________________________________________________________________
//
// EOF
//
