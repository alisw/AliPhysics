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

//////////////////////////////////////////////////////////////////////////////
//                                                                          
// Concrete implementation of AliFMDSubDetector 
//
// This implements the geometry for FMD3
//
//////////////////////////////////////////////////////////////////////////////
#ifndef ALIFMD3_H
# include "AliFMD3.h"
#endif 
#ifndef ROOT_TVirtualMC
# include <TVirtualMC.h>
#endif
#ifndef ALILOG_H
# include "AliLog.h"
#endif

//____________________________________________________________________
ClassImp(AliFMD3);

//____________________________________________________________________
AliFMD3::AliFMD3() 
  : AliFMDSubDetector(3) 
{}

//____________________________________________________________________
AliFMD3::~AliFMD3() 
{}

//____________________________________________________________________
void 
AliFMD3::SetupGeometry(Int_t airId, Int_t kaptionId) 
{
  fInnerHoneyLowR  = fInner->GetLowR() + 1;
  fInnerHoneyHighR = fInner->GetHighR() + 1;
  fOuterHoneyLowR  = fOuter->GetLowR() + 1;
  fOuterHoneyHighR = fOuter->GetHighR() + 1;

  Double_t par[3];
  par[0] = fInner->GetLowR();
  par[1] = fOuterHoneyHighR;
  par[2] = fDz = (TMath::Abs(fInnerZ - fOuterZ) 
		  + fInner->GetSiThickness() 
		  + fInner->GetPrintboardThickness() 
		  + fInner->GetLegLength() 
		  + fInner->GetModuleSpacing() 
		  + fHoneycombThickness) / 2;
  fVolumeId = gMC->Gsvolu("FMD3", "TUBE", airId, par, 3);

  gMC->Matrix(fRotationId, 90, 0, 90, 90, 0, 0); 
  //0, 180, 90, 90, 180, 0);

  AliFMDSubDetector::SetupGeometry(airId, kaptionId);
}

//____________________________________________________________________
void 
AliFMD3::Geometry(const char* mother, Int_t pbRotId, 
		  Int_t idRotId, Double_t z) 
{
  z = fInnerZ - fDz;
  gMC->Gspos("FMD3", 1, mother, 0, 0, z, fRotationId);
  
  AliFMDSubDetector::Geometry("FMD3", pbRotId, idRotId, z);
}

  

//____________________________________________________________________
//
// EOF
//
