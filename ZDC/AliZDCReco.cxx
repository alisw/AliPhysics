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
AliZDCReco::AliZDCReco(Float_t ezn, Float_t ezp, Float_t ezdc, Float_t ezem,
     Int_t detspn, Int_t detspp, Int_t trspn, Int_t trspp, Int_t trsp, Int_t part, Float_t b)
{ 
  //
  // Standard constructor
  //
  fZNenergy  = ezn;
  fZPenergy  = ezp;
  fZDCenergy = ezdc;
  fZEMenergy = ezem;
  fNDetSpecN = detspn;
  fNDetSpecP = detspp;
  fNTrueSpecN = trspn;
  fNTrueSpecP = trspp;
  fNTrueSpec = trsp;
  fNPart     = part;
  fImpPar    = b;
  
}

//______________________________________________________________________________
void AliZDCReco::Print(Option_t *) const {
  //
  // Printing Reconstruction Parameters
  //
  printf("	---   Reconstruction -> EZN = %f TeV, EZP = %f TeV, EZDC = %f TeV,"
	 " EZEM = %f GeV \n 		NDetSpecN = %d, NDetSpecP = %d, Nspecn = %d,"
	 " Nspecp = %d, Npart = %d, b = %f fm.\n ", 
	 fZNenergy,fZPenergy,fZDCenergy,fZEMenergy,fNDetSpecN,fNDetSpecP,
	 fNTrueSpecN,fNTrueSpecP,fNPart,fImpPar);
}
