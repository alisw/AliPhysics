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

/*
$Log$
*/

#include "AliMevSimParticle.h"


ClassImp(AliMevSimParticle)

   
///////////////////////////////////////////////////////////////////////////////////////

AliMevSimParticle::AliMevSimParticle()
  : TMevSimPartTypeParams() {

}

///////////////////////////////////////////////////////////////////////////////////////

AliMevSimParticle::AliMevSimParticle(PDG_t pdg, Int_t multmean, Int_t multvc, 
		  Float_t tempmean, Float_t tempstdev, Float_t sigmamean,
		  Float_t sigmastdev, Float_t expvelmean, Float_t expvelstdev)

  : TMevSimPartTypeParams(0, multmean, multvc, tempmean, tempstdev, 
			  sigmamean, sigmastdev, expvelmean, expvelstdev)  {


  // Calculate geant ID from pdg
  fConv = new TMevSimConverter();
  fPdg = pdg;
  if (fConv) fGPid = fConv->IdFromPDG(pdg);  

}

///////////////////////////////////////////////////////////////////////////////////////

AliMevSimParticle::~AliMevSimParticle() {
}

///////////////////////////////////////////////////////////////////////////////////////

void  AliMevSimParticle::SetPDG(PDG_t pdg) {

  fPdg = pdg;
  fGPid = fConv->IdFromPDG(pdg);
}

///////////////////////////////////////////////////////////////////////////////////////

PDG_t AliMevSimParticle::GetPDG() {
  
  fPdg = (PDG_t)fConv->PDGFromId(fGPid);
  return fPdg;
}

///////////////////////////////////////////////////////////////////////////////////////


  

