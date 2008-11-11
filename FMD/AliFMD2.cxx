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
/** @file    AliFMD2.cxx
    @author  Christian Holm Christensen <cholm@nbi.dk>
    @date    Sun Mar 26 18:25:51 2006
    @brief   Concrete implementation of AliFMDDetector for FMD2
*/
//____________________________________________________________________
//                                                                          
// Concrete implementation of AliFMDDetector 
//
// This implements the geometry for FMD2
// The FMD2 has two ring, one of both types. 
// FMD2 is mounted on the space-frame via 4 flanges
// Support is not fleshed ot yet. 
// Support will be simple compared to FMD3.
//
#include "AliFMD2.h"		// ALIFMD2_H 
#include "AliLog.h"
// #include "AliFMDRing.h"		// ALIFMDRING_H 

//====================================================================
ClassImp(AliFMD2)
#if 0
  ; // This is here to keep Emacs for indenting the next line
#endif

//____________________________________________________________________
AliFMD2::AliFMD2(AliFMDRing* inner, AliFMDRing* outer) 
  : AliFMDDetector(2, inner, outer)
{
  // Constructor 
  // SetInnerZ(83.4);
  // SetOuterZ(75.2);
  Double_t off = 0.414256-0.1963; // 2.35
  if (off != 0) 
    AliWarning(Form("Z position of FMD2 rings may be wrong by %fcm!", off));
  SetInnerZ(83.4+off);
  SetOuterZ(75.2+off);
}


//____________________________________________________________________
void
AliFMD2::Init() 
{
  // Initialize 
  AliFMDDetector::Init();
  SetInnerHoneyHighR(GetOuterHoneyHighR());
}

//____________________________________________________________________
//
// EOF
//
