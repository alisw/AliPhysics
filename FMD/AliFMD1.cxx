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
// Concrete implementation of AliFMDSubDetector 
//
// This implements the geometry for FMD1 
//
#include "AliFMD1.h"		// ALIFMD1_H 
#include "AliFMDRing.h"		// ALIFMDRING_H 
#include "TVirtualMC.h"		// ROOT_TVirtualMC
#include "AliLog.h"		// ALILOG_H


//____________________________________________________________________
ClassImp(AliFMD1);

//____________________________________________________________________
AliFMD1::AliFMD1() 
  : AliFMDSubDetector(1) 
{
  // Default constructor for the FMD1 sub-detector 
}

//____________________________________________________________________
AliFMD1::~AliFMD1() 
{
  // Destructor - does nothing 
}

//____________________________________________________________________
void 
AliFMD1::SetupGeometry(Int_t airId, Int_t alId, Int_t /* cId */) 
{
  // Setup the FMD1 sub-detector geometry 
  // 
  // Parameters:
  // 
  //     airId         Id # of the Air medium 
  //     alId     Id # of the Aluminium medium 
  // 
  AliDebug(10, "\tDefining the geometry for FMD1");
  fInnerHoneyLowR  = fInner->GetLowR() + 1;
  fInnerHoneyHighR = fInner->GetHighR() + 1;
  fOuterHoneyLowR  = 0;
  fOuterHoneyHighR = 0;

  Double_t par[3];
  par[0] = fInner->GetLowR();
  par[1] = fInnerHoneyHighR;
  par[2] = fDz = (fInner->GetLegLength() 
		  + fInner->GetSiThickness() 
		  + fInner->GetPrintboardThickness() 
		  + fInner->GetModuleSpacing() 
		  + fHoneycombThickness) / 2;
  fVolumeId = gMC->Gsvolu("FMD1", "TUBE", airId, par, 3);

  // Rotate the full sub-detector 
  gMC->Matrix(fRotationId, 270, 180, 90, 90, 180, 0); 

  AliFMDSubDetector::SetupGeometry(airId, alId);
}

//____________________________________________________________________
void 
AliFMD1::Geometry(const char* mother, Int_t pbRotId, 
		  Int_t idRotId, Double_t z) 
{
  // Position the FMD1 sub-detector volume 
  // 
  // Parameters 
  //
  //     mother     name of the mother volume 
  //     pbRotId    Printboard roation matrix ID 
  //     idRotId    Identity rotation matrix ID 
  //     z          Z position (not really used here, but passed down)
  //
  // The Z passed in isn't used. 
  z = fInnerZ + fDz;
  AliDebug(10, Form("\tPutting FMD1 in %s at z=%lf cm", mother, z));
  gMC->Gspos("FMD1", 1, mother, 0, 0, z, fRotationId, "ONLY");

  AliFMDSubDetector::Geometry("FMD1", pbRotId, idRotId, z);
}

  

//____________________________________________________________________
//
// EOF
//
