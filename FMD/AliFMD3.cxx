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
// This implements the geometry for FMD3
//
#ifndef ROOT_TVirtualMC
# include <TVirtualMC.h>
#endif
#ifndef ROOT_TCONS
# include <TCONS.h>
#endif
#ifndef ROOT_TNode
# include <TNode.h>
#endif
#ifndef ROOT_TList
# include <TList.h>
#endif
#ifndef ALIFMD3_H
# include "AliFMD3.h"
#endif 
#ifndef ALILOG_H
# include "AliLog.h"
#endif
#ifndef ALIFMDRING_H
# include "AliFMDRing.h"
#endif 
#include <Riostream.h>

//____________________________________________________________________
ClassImp(AliFMD3);

//____________________________________________________________________
AliFMD3::AliFMD3() 
  : AliFMDSubDetector(3), 
    fVolumeId(0),
    fDz(0)
{
  // Default constructor for the FMD3 sub-detector 
}


//____________________________________________________________________
AliFMD3::~AliFMD3() 
{
  // Destructor - does nothing 
}


//____________________________________________________________________
void 
AliFMD3::SetupGeometry(Int_t airId, Int_t kaptionId) 
{
  // Setup the FMD3 sub-detector geometry 
  // 
  // Parameters:
  // 
  //     airId         Id # of the Air medium 
  //     kaptionId     Id # of the Aluminium medium 
  // 
  fInnerHoneyLowR  = fInner->GetLowR() + 1;
  fInnerHoneyHighR = fInner->GetHighR() + 1;
  fOuterHoneyLowR  = fOuter->GetLowR() + 1;
  fOuterHoneyHighR = fOuter->GetHighR() + 1;

  CalculateDz();
  Double_t par[3];
  par[0] = fInner->GetLowR();
  par[1] = fOuterHoneyHighR;
  par[2] = fDz;
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
  // Position the FMD3 sub-detector volume 
  // 
  // Parameters 
  //
  //     mother     name of the mother volume 
  //     pbRotId    Printboard roation matrix ID 
  //     idRotId    Identity rotation matrix ID 
  //     z          Z position (not really used here, but passed down)
  //
  z = fInnerZ - fDz;
  gMC->Gspos("FMD3", 1, mother, 0, 0, z, fRotationId);
  
  AliFMDSubDetector::Geometry("FMD3", pbRotId, idRotId, z);
}

  
//____________________________________________________________________
void 
AliFMD3::SimpleGeometry(TList* nodes, 
			TNode* mother, 
			Int_t colour, 
			Double_t zMother) 
{
  // We need to get the equation for the line that connects the 
  // outer circumfrences of the two rings, as  well as for the line
  // that connects the inner curcumfrences, so that we can project to
  // where the honey-comb actually ends. 
  // 
  // we have 
  //   
  //   y = a * x + b 
  //   b = y - a * x;
  // 
  // For the outer line, we have the two equations 
  // 
  //    fOuterHoneyHighR = a * x1 + b;
  //    fInnerHoneyHighR = a * x2 + b; 
  // 
  // where 
  // 
  //    x1 = (fOuterZ + fOuter->fSiThickness + fOuter->fPrintboardThickness 
  //          + fOuter->fLegLength + fModuleSpacing) 
  //       = fInner - fDz + fHoneycombThickness
  //    x2 = (fInnerZ + fInner->fSiThickness + fInner->fPrintboardThickness 
  //          + fInner->fLegLength + fModuleSpacing)
  // 
  // and 
  //
  //    a  = (fOuterHoneyHighR - fInnerHoneyHighR) / (x1 - x2)
  //    
  // 
  CalculateDz();
#if 1
  Double_t x1  = (fOuterZ - (fOuter->GetSiThickness() 
			     + fOuter->GetPrintboardThickness() 
			     + fOuter->GetLegLength() 
			     + fOuter->GetModuleSpacing()));
  Double_t x2  = (fInnerZ - (fInner->GetSiThickness() 
			     + fInner->GetPrintboardThickness() 
			     + fInner->GetLegLength() 
			     + fInner->GetModuleSpacing()));
  Double_t ao   = 0;
  Double_t ao1  = (fOuterHoneyHighR - fInnerHoneyHighR) / (x1 - x2);
  Double_t ao2  = ((fOuter->GetHighR() - fInner->GetHighR()) 
		   / (fOuterZ - fInnerZ));
  Double_t bo   = 0;
  if (ao2 > ao1) {
    cout << "Wafer determinds the size" << endl;
    ao  = ao2;
    bo  = fInner->GetHighR() - ao * fInnerZ;
  }
  else {
    ao = ao1;
    bo = fOuterHoneyHighR - ao * x1;
  }
  
  Double_t y1o = ao * (fInnerZ - 2 * fDz) + bo;
  Double_t y2o = ao * fInnerZ + bo;
#endif
  // We probably need to make a PCON here. 
  TShape* shape = new TCONS("FMD3", "FMD3", "", fDz, 
			    fOuter->GetLowR(),  y1o, /* fOuterHoneyHighR, */
			    fInner->GetLowR(),  y2o, /* fInnerHoneyHighR, */
			    0, 360);
  mother->cd();
  zMother = fInnerZ - fDz;  
  TNode* node = new TNode("FMD3", "FMD3", shape, 0, 0, zMother, 0);
  node->SetVisibility(0);
  nodes->Add(node);
  AliFMDSubDetector::SimpleGeometry(nodes, node, colour, zMother);
}

//____________________________________________________________________
void 
AliFMD3::CalculateDz() 
{
  if (fDz > 0) return;
  fDz = (TMath::Abs(fInnerZ - fOuterZ) 
	 + fOuter->GetSiThickness() 
	 + fOuter->GetPrintboardThickness() 
	 + fOuter->GetLegLength() 
	 + fOuter->GetModuleSpacing() 
	 + fHoneycombThickness) / 2;
}

//____________________________________________________________________
//
// EOF
//
