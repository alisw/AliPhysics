/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *;
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id$ */

////////////////////////////////////////////////
//  RecPoints classes for set ZDC             //
//  This class reconstructs the space         //
//  points from digits                        //
//  for the ZDC calorimeter                   //
////////////////////////////////////////////////


#include "AliZDCReco.h"

ClassImp(AliZDCReco)
  

//_____________________________________________________________________________
AliZDCReco::AliZDCReco() :
	
  TObject(),
  fZN1energy(0),
  fZP1energy(0),
  fZDC1energy(0),
  fZN2energy(0),
  fZP2energy(0),
  fZDC2energy(0),
  fZEMenergy(0),
  fNDetSpecNLeft(0),
  fNDetSpecPLeft(0),
  fNDetSpecNRight(0),
  fNDetSpecPRight(0),
  fNTrueSpecN(0),
  fNTrueSpecP(0),
  fNTrueSpec(0),
  fNPart(0),
  fImpPar(0)

{ 
  //
  // Default constructor
  //
}
  

//_____________________________________________________________________________
AliZDCReco::AliZDCReco(Float_t ezn1, Float_t ezp1, Float_t ezdc1, Float_t ezem,
	Float_t ezn2, Float_t ezp2, Float_t ezdc2,
     	Int_t detspnLeft, Int_t detsppLeft, Int_t detspnRight, Int_t detsppRight,
	Int_t trspn, Int_t trspp, Int_t trsp, Int_t part,
	Float_t b) :
	
  TObject(),
  fZN1energy(ezn1),
  fZP1energy(ezp1),
  fZDC1energy(ezdc1),
  fZN2energy(ezn2),
  fZP2energy(ezp2),
  fZDC2energy(ezdc2),
  fZEMenergy(ezem),
  fNDetSpecNLeft(detspnLeft),
  fNDetSpecPLeft(detsppLeft),
  fNDetSpecNRight(detspnRight),
  fNDetSpecPRight(detsppRight),
  fNTrueSpecN(trspn),
  fNTrueSpecP(trspp),
  fNTrueSpec(trsp),
  fNPart(part),
  fImpPar(b)

{ 
  //
  // Constructor
  //
  
}

//______________________________________________________________________________
AliZDCReco::AliZDCReco(const AliZDCReco &oldreco) :

  TObject()
{
  // Copy constructor

  fZN1energy = oldreco.GetZN1energy();
  fZP1energy = oldreco.GetZP1energy();            
  fZDC1energy = oldreco.GetZDC1energy();	  
  fZN2energy = oldreco.GetZN2energy();	   
  fZP2energy = oldreco.GetZP2energy(); 	   
  fZDC2energy = oldreco.GetZDC2energy();	  
  fZEMenergy = oldreco.GetZEMenergy();	   
  fNDetSpecNLeft = oldreco.GetNDetSpecNLeft();	
  fNDetSpecPLeft = oldreco.GetNDetSpecPLeft();	
  fNDetSpecNRight = oldreco.GetNDetSpecNRight();	
  fNDetSpecPRight = oldreco.GetNDetSpecPRight();	
  fNTrueSpecN = oldreco.GetNTrueSpecN();	  
  fNTrueSpecP = oldreco.GetNTrueSpecP();	  
  fNTrueSpec = oldreco.GetNTrueSpec();	  
  fNPart = oldreco.GetNPart();			 
  fImpPar = oldreco.GetImpPar();			 
}

//______________________________________________________________________________
void AliZDCReco::Print(Option_t *) const {
  //
  // Printing Reconstruction Parameters
  //
  printf("	---   Reconstruction -> EZN = %f TeV, EZP = %f TeV, EZDC = %f TeV,"
	 " EZEM = %f GeV \n 		NDetSpecN = %d, NDetSpecP = %d, Nspecn = %d,"
	 " Nspecp = %d, Npart = %d, b = %f fm.\n ", 
	 fZN1energy,fZP1energy,fZDC1energy,fZEMenergy,fNDetSpecNLeft,
	 fNDetSpecPLeft,fNTrueSpecN,fNTrueSpecP,fNPart,fImpPar);
}
