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
#include "TVirtualMC.h"		// ROOT_TVirtualMC
#include "AliFMD3Support.h"	// ALIFMD3SUPPORT_H 
#include "AliLog.h"		// ALILOG_H
#include <Riostream.h>		// ROOT_Riostream

//____________________________________________________________________
ClassImp(AliFMD3Support);

//____________________________________________________________________
const Char_t* AliFMD3Support::fgkNoseName   = "F3SN";
const Char_t* AliFMD3Support::fgkBackName   = "F3SB";
const Char_t* AliFMD3Support::fgkBeamName   = "F3SL";
const Char_t* AliFMD3Support::fgkFlangeName = "F3SF";

//____________________________________________________________________
AliFMD3Support::AliFMD3Support() 
  : fZ(0),
    fAlpha(-1),
    fNoseId(-1),
    fBeamId(-1),
    fBackId(-1), 
    fFlangeId(-1)
{
  // Default constructor for the support of FMD3 sub-detector 
  SetNoseZ();
  SetNoseLowR();
  SetNoseHighR();
  SetNoseLength();
  SetBackLowR();
  SetBackHighR();
  SetBackLength();
  SetBeamThickness();
  SetBeamWidth();
  SetConeLength();
  SetFlangeR();
  SetNBeam();
  SetNFlange();
}


//____________________________________________________________________
AliFMD3Support::~AliFMD3Support() 
{
  // Destructor - does nothing 
}


//____________________________________________________________________
void 
AliFMD3Support::SetupGeometry(Int_t    airId, 
			      Int_t    cId, 
			      Double_t innerZl, 
			      Double_t innerZh, 
			      Double_t innerRl, 
			      Double_t /* outerZl */, 
			      Double_t outerZh, 
			      Double_t outerRl)
{
  // Setup the FMD3 sub-detector geometry 
  // 
  // Parameters:
  // 
  //     airId         Id # of the Air medium 
  //     cId           Id # of the Carbon fibre medium 
  // 

  // Global stuff we need 
  Double_t zdist = fConeLength - fBackLength - fNoseLength;
  Double_t tdist = fBackHighR - fNoseHighR;
  Double_t beaml = TMath::Sqrt(zdist * zdist + tdist * tdist);
  Double_t theta = -180. * TMath::ATan2(tdist, zdist) / TMath::Pi();
  Double_t minZ  = TMath::Min(fNoseZ - fConeLength, outerZh);
  fAlpha         = tdist / zdist;
  fZ             = fNoseZ + (minZ - fNoseZ) / 2;
  AliDebug(30, Form("\tTheta = %lf", theta));
  
  const Char_t* mother = "FMD3";
  Double_t p[3 + 9 * 3];
  Double_t eps = 0;
  // ------------- Mother volume -------------------------------------
  // The planes should be defined with increasing Z, as it will become
  // invalid if not 
  // Global parameters 
  p[0]  = 0;
  p[1]  = 360;
  p[2]  = 8;
  // First plane (at back of back or outer ring)
  p[3]  = minZ - fZ - eps;
  p[4]  = outerRl - eps;				   
  p[5]  = fFlangeR + eps;
  // Second plane (at front of back, at end of flanges)
  p[6]  = fNoseZ - zdist - fNoseLength  - fZ  + eps;
  p[7]  = p[4];
  p[8]  = p[5];	     
  // Third plane (at front of back)
  p[9]  = p[6] - eps/2; 
  p[10] = p[7];
  p[11] = ConeR(p[9] + fZ)  + eps;		   
  // Fourth plane (at back of inner ring) 
  p[12] = innerZh - fZ - eps;	  
  p[13] = outerRl - eps;	  
  p[14] = ConeR(p[12] + fZ) + eps;
  // Fifth plane (at back of inner ring) 
  p[15] = p[12] - eps/2;	  
  p[16] = innerRl - eps;	 
  p[17] = ConeR(p[15] + fZ) + eps;
  // Sixth plane  (at front of inner ring)
  p[18] = innerZl - fZ + eps;	 
  p[19] = p[16];		  
  p[20] = ConeR(p[18] + fZ) + eps;
  // Seventh plane (at end of nose)
  p[21] = fNoseZ - fNoseLength  - fZ - eps;	
  p[22] = fNoseLowR - eps;			
  p[23] = ConeR(p[21] + fZ) + eps; // fNoseHighR;
  // Eight (and final) plane (at start of nose)
  p[24] = fNoseZ  - fZ + eps;
  p[25] = p[22];
  p[26] = fNoseHighR + eps;  

  // The volume 
  gMC->Gsvolu(mother, "PCON", airId, p, 27);
  
  // ------------- Support Structures --------------------------------
  fRotations.Set(fNBeam + fNFlange);
  Double_t par[3];
  
  // The nose 
  par[0]  = fNoseLowR;
  par[1]  = fNoseHighR;
  par[2]  = fNoseLength / 2;
  fNoseId = gMC->Gsvolu(fgkNoseName, "TUBE", cId, par, 3);

  // The Back 
  par[0]  = fBackLowR;
  par[1]  = fBackHighR;
  par[2]  = fBackLength / 2;
  fBackId = gMC->Gsvolu(fgkBackName, "TUBE", cId, par, 3);

  // The Beams 
  par[0]  = fBeamThickness  / 2;
  par[1]  = fBeamWidth / 2;
  par[2]  = beaml / 2;
  fBeamId = gMC->Gsvolu(fgkBeamName, "BOX", cId, par, 3);
  for (Int_t i = 0; i < fNBeam; i++) {
    // cout << "Making beam # " << i << endl;
    Double_t phi = 360. / fNBeam * i;
    Int_t    id;
    gMC->Matrix(id, 180 - theta, phi, 90, 90 + phi, theta, phi);
    fRotations[i] = id;
  }

  // The Flanges 
  par[0]    = (fFlangeR - fBackHighR) / 2;
  par[1]    = fBeamWidth / 2;
  par[2]    = fBackLength  / 2;
  fFlangeId = gMC->Gsvolu(fgkFlangeName, "BOX", cId, par, 3);
  for (Int_t i = 0; i < fNFlange; i++) {
    Double_t phi = 360. / fNFlange * i + 180. / fNFlange;
    Int_t id;
    gMC->Matrix(id, 90, phi, 90, 90+phi, 0, 0);
    fRotations[fNBeam + i] = id;
  }
}

//____________________________________________________________________
void 
AliFMD3Support::Geometry(const char* mother, Int_t idRotId, Double_t zTop) 
{
  // Position the FMD3 sub-detector volume 
  // 
  // Parameters 
  //
  //     mother     name of the mother volume 
  //     idRotId    Identity rotation matrix ID 
  //     z          Z position
  //

  // Common parameters 
  Double_t zdist = fConeLength - fBackLength - fNoseLength;
  Double_t tdist = fBackHighR - fNoseHighR;
  const Char_t* name = "FMD3";
  Double_t z = zTop;
  
  // Placing mother volume 
  AliDebug(10, Form("\tPutting %s in %s at z=%lf", name, mother, zTop));
  gMC->Gspos(name, 1, mother, 0, 0, zTop, idRotId, "ONLY");

  // Placing the nose 
  z       = fNoseZ - fNoseLength / 2 - fZ;
  AliDebug(10, Form("\tPutting %s in %s at z=%lf-%lf/2-%lf=%lf", 
		    fgkNoseName, name, fNoseZ, fNoseLength, fZ, z));
  gMC->Gspos(fgkNoseName, 1, name, 0., 0., z, idRotId, "");

  // Placing  the back 
  z       = fNoseZ - fNoseLength - zdist - fBackLength / 2 - fZ;
  AliDebug(10, Form("\tPutting %s in %s at z=%lf-%lf-%lf-%lf/2-%lf=%lf", 
		    fgkBackName, name, fNoseZ, fNoseLength, zdist, 
		    fBackLength, fZ, z));
  gMC->Gspos(fgkBackName, 1, name, 0., 0., z, idRotId, "");

  // Placing the beams 
  z          = fNoseZ - fNoseLength - zdist / 2 - fZ;
  Double_t r = fNoseHighR + tdist / 2;
  AliDebug(10, Form("\tPutting %s's in %s at z=%lf-%lf-%lf/2-%lf=%lf", 
		    fgkBeamName, name, fNoseZ, fNoseLength, zdist, fZ, z));
  for (Int_t i = 0; i < fNBeam; i++) {
    // cout << "Making beam # " << i << endl;
    Double_t phi = 360. / fNBeam * i;
    gMC->Gspos(fgkBeamName, i, name, 
	       r * TMath::Cos(TMath::Pi() / 180 * phi), 
	       r * TMath::Sin(TMath::Pi() / 180 * phi), 
	       z, fRotations[i], "");
  }

  // Placing the flanges 
  r         = fBackHighR + (fFlangeR - fBackHighR) / 2;
  z         = fNoseZ - fNoseLength - zdist - fBackLength / 2 - fZ;
  AliDebug(10, Form("\tPutting %s in %s at z=%lf-%lf-%lf-%lf/2-%lf=%lf", 
		    fgkFlangeName, name, fNoseZ, fNoseLength, zdist, 
		    fBackLength, fZ, z));
  for (Int_t i = 0; i < fNFlange; i++) {
    Double_t phi = 360. / fNFlange * i + 180. / fNFlange;
    gMC->Gspos(fgkFlangeName, i, name, 
	       r * TMath::Cos(TMath::Pi() / 180 * phi), 
	       r * TMath::Sin(TMath::Pi() / 180 * phi), 
	       z, fRotations[fNBeam + i], "");
  }

}

//____________________________________________________________________
Double_t
AliFMD3Support::ConeR(Double_t z, Option_t* opt) const
{
  // Calculate the cone radius at Z 
  if (fAlpha < 0) {
    Warning("ConeR", "alpha not set: %lf", fAlpha);
    return -1;
  }
  if (z > fNoseZ) { 
    Warning("ConeR", "z=%lf is before start of cone %lf", z, fNoseZ);
    return -1;
  }
  Double_t e = fBeamThickness / TMath::Cos(TMath::ATan(fAlpha));
  if (opt[0] == 'I' || opt[1] == 'i') e *= -1;
  if (z > fNoseZ - fNoseLength) return fNoseHighR + e; 
  if (z < fNoseZ - fConeLength + fBackLength) return fBackHighR + e; 
  Double_t r = fNoseHighR + fAlpha * TMath::Abs(z - fNoseZ + fNoseLength) + e;
  return r;
}

//____________________________________________________________________
void 
AliFMD3Support::Gsatt() const
{
  // Set drawing attributes for the FMD3 Support 
  gMC->Gsatt(fgkNoseName, "SEEN", 1);
  gMC->Gsatt(fgkBeamName, "SEEN", 1);
  gMC->Gsatt(fgkBackName, "SEEN", 1);
  gMC->Gsatt(fgkFlangeName, "SEEN", 1);
}


//____________________________________________________________________
//
// EOF
//
