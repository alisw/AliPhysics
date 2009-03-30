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
  //
  fNDetSpecNLeft(0),
  fNDetSpecPLeft(0),
  fNDetSpecNRight(0),
  fNDetSpecPRight(0),
  fNTrueSpecNLeft(0),
  fNTrueSpecPLeft(0),
  fNTrueSpecLeft(0),
  fNTrueSpecNRight(0),
  fNTrueSpecPRight(0),
  fNTrueSpecRight(0),
  fNPartLeft(0),
  fNPartRight(0),
  fImpPar(0)

{ 
  //
  // Default constructor
  //
  for(Int_t i=0; i<10; i++){
     fZN1EnTow[i] = fZP1EnTow[i] = fZN2EnTow[i] = fZP2EnTow[i] = 0.;
     if(i<2){
       fZN1Energy[i] = fZP1Energy[i] = fZN2Energy[i] = fZP2Energy[i] =  0.;
       fZEM1signal[i] = fZEM2signal[i] = 0.;
       fPMRef1[i] = fPMRef2[i] = 0.;
     }
  }
}
  

//_____________________________________________________________________________
AliZDCReco::AliZDCReco(Float_t* ezn1, Float_t* ezp1, Float_t* ezn2, Float_t* ezp2,  
	     Float_t* ezn1tow, Float_t* ezp1tow,
	     Float_t* ezn2tow, Float_t* ezp2tow, 
	     Float_t* ezem1, Float_t* ezem2, 
	     Float_t* ref1, Float_t* ref2, 
	     //	   
	     Int_t detspnLeft,  Int_t detsppLeft, Int_t detspnRight, Int_t detsppRight,  
	     Int_t trspnLeft, Int_t trsppLeft, Int_t trspLeft, 
	     Int_t trspnRight, Int_t trsppRight, Int_t trspRight,
	     Int_t partLeft, Int_t partRight, Float_t b) :
	
  TObject(),
  //
  fNDetSpecNLeft(detspnLeft),
  fNDetSpecPLeft(detsppLeft),
  fNDetSpecNRight(detspnRight),
  fNDetSpecPRight(detsppRight),
  fNTrueSpecNLeft(trspnLeft),
  fNTrueSpecPLeft(trsppLeft),
  fNTrueSpecLeft(trspLeft),
  fNTrueSpecNRight(trspnRight),
  fNTrueSpecPRight(trsppRight),
  fNTrueSpecRight(trspRight),
  fNPartLeft(partLeft),
  fNPartRight(partRight),
  fImpPar(b)

{ 
  //
  // Constructor
  //
  for(Int_t j=0; j<10; j++){
     fZN1EnTow[j] =  ezn1tow[j];
     fZP1EnTow[j] =  ezp1tow[j];
     fZN2EnTow[j] =  ezn2tow[j];
     fZP2EnTow[j] =  ezp2tow[j];
     if(j<2){
       fZN1Energy[j] = ezn1[j];
       fZP1Energy[j] = ezp1[j];
       fZN2Energy[j] = ezn2[j];
       fZP2Energy[j] = ezp2[j];
       fZEM1signal[j] = ezem1[j];
       fZEM2signal[j] = ezem2[j];
       fPMRef1[j] = ref1[j];
       fPMRef2[j] = ref2[j];
     }
  }
  
}

//______________________________________________________________________________
AliZDCReco::AliZDCReco(const AliZDCReco &oldreco) :
TObject(),
fNDetSpecNLeft(oldreco.GetNDetSpecNLeft()),
fNDetSpecPLeft(oldreco.GetNDetSpecPLeft()),
fNDetSpecNRight(oldreco.GetNDetSpecNRight()),        
fNDetSpecPRight(oldreco.GetNDetSpecPRight()),       
fNTrueSpecNLeft(oldreco.GetNTrueSpecNLeft()), 	
fNTrueSpecPLeft(oldreco.GetNTrueSpecPLeft()), 	
fNTrueSpecLeft(oldreco.GetNTrueSpecLeft()),
fNTrueSpecNRight(oldreco.GetNTrueSpecNRight()),	
fNTrueSpecPRight(oldreco.GetNTrueSpecPRight()),	
fNTrueSpecRight(oldreco.GetNTrueSpecRight()),  	
fNPartLeft(oldreco.GetNPartLeft()),		       
fNPartRight(oldreco.GetNPartRight()),  		       
fImpPar(oldreco.GetImpPar())      
{
  // Copy constructor

  fZN1Energy[0]  = oldreco.GetZN1HREnergy();
  fZP1Energy[0]  = oldreco.GetZP1HREnergy();		
  fZN2Energy[0]  = oldreco.GetZN2HREnergy();	     
  fZP2Energy[0]  = oldreco.GetZP2HREnergy();	 
  //    
  fZN1Energy[1]  = oldreco.GetZN1LREnergy();
  fZP1Energy[1]  = oldreco.GetZP1LREnergy();	       
  fZN2Energy[1]  = oldreco.GetZN2LREnergy();	    
  fZP2Energy[1]  = oldreco.GetZP2LREnergy();	    
  //
  for(Int_t i=0; i<5; i++){	  
     fZN1EnTow[i]  = oldreco.GetZN1HREnTow(i);
     fZP1EnTow[i]  = oldreco.GetZP1HREnTow(i);
     fZN2EnTow[i]  = oldreco.GetZN2HREnTow(i);
     fZP2EnTow[i]  = oldreco.GetZP2HREnTow(i);
     fZN1EnTow[i+5]  = oldreco.GetZN1LREnTow(i);
     fZP1EnTow[i+5]  = oldreco.GetZP1LREnTow(i);
     fZN2EnTow[i+5]  = oldreco.GetZN2LREnTow(i);
     fZP2EnTow[i+5]  = oldreco.GetZP2LREnTow(i);
  }
  fZEM1signal[0] = oldreco.GetZEM1HRsignal();
  fZEM1signal[1] = oldreco.GetZEM1LRsignal();
  fZEM2signal[0] = oldreco.GetZEM2HRsignal();
  fZEM2signal[1] = oldreco.GetZEM2LRsignal();
  fPMRef1[0] = oldreco.GetPMRef1HRsignal();
  fPMRef1[1] = oldreco.GetPMRef1LRsignal();
  fPMRef2[0] = oldreco.GetPMRef2HRsignal();
  fPMRef2[1] = oldreco.GetPMRef2LRsignal();
}

//______________________________________________________________________________
void AliZDCReco::Print(Option_t *) const {
  //
  // Printing Reconstruction Parameters
  //
  printf(" \t ---   Reconstruction -> EZN1 = %f TeV, EZP1 = %f TeV,  EZEM1 = %f GeV ,  EZEM2 = %f GeV \n "		
	 "EZN2 = %f TeV, EZP2 = %f TeV \n"
	 " \t NDetSpecNLeft = %d, NDetSpecPLeft = %d, NspecnLeft = %d,"
	 " NspecpLeft = %d, NpartLeft = %d"
	 " \t NDetSpecNRight = %d, NDetSpecPRight = %d, NspecnRight = %d,"
	 " NspecpRight = %d, NpartRight = %d"
	 " \t b = %f fm\n ", 
	 fZN1Energy[0]/1000.,fZP1Energy[0]/1000.,fZEM1signal[0]/1000.,fZEM2signal[0]/1000., 
	 fZN2Energy[0]/1000., fZP2Energy[0]/1000.,
	 fNDetSpecNLeft,fNDetSpecPLeft,fNTrueSpecNLeft,fNTrueSpecPLeft,fNPartLeft,
	 fNDetSpecNRight,fNDetSpecPRight,fNTrueSpecNRight,fNTrueSpecPRight,fNPartRight,
	 fImpPar);
}
