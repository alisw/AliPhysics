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
  fNDetSpecNSideA(0),
  fNDetSpecPSideA(0),
  fNDetSpecNSideC(0),
  fNDetSpecPSideC(0),
  fNTrueSpectators(0),
  fNTrueSpecSideA(0),
  fNTrueSpecSideC(0),
  fNParticipants(0),
  fNPartSideA(0),
  fNPartSideC(0),
  fImpParameter(0),
  fImpParSideA(0),
  fImpParSideC(0),
  fRecoFlag(0x0),
  fEnergyFlag(kFALSE),
  fIsScalerOn(kFALSE)
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
  
  for(Int_t i=0; i<32; i++){
    fZDCScaler[i] = 0;
    for(Int_t ij=0; ij<4; ij++) fZDCTDCData[i][ij] = 0;
  }
}
  

//_____________________________________________________________________________
AliZDCReco::AliZDCReco(
     Float_t* ezn1, Float_t* ezp1, Float_t* ezn2, Float_t* ezp2,  
     Float_t* ezn1tow, Float_t* ezp1tow, Float_t* ezn2tow, Float_t* ezp2tow, 
     Float_t* ezem1, Float_t* ezem2, Float_t* ref1, Float_t* ref2, 
     //    
     Int_t detspnSideA,  Int_t detsppSideA, Int_t detspnSideC, Int_t detsppSideC,  
     Int_t trsp, Int_t trspSideA,Int_t trspSideC,
     Int_t npart, Int_t npartSideA, Int_t npartSideC, 
     Float_t b, Float_t bSideA, Float_t bSideC,
     UInt_t recoFlag, Bool_t energyFlag, Bool_t scalerOn, 
     UInt_t* scaler, Int_t tdcData[32][4]) :
	
  TObject(),
  //
  fNDetSpecNSideA(detspnSideA),
  fNDetSpecPSideA(detsppSideA),
  fNDetSpecNSideC(detspnSideC),
  fNDetSpecPSideC(detsppSideC),
  fNTrueSpectators(trsp),
  fNTrueSpecSideA(trspSideA),
  fNTrueSpecSideC(trspSideC),
  fNParticipants(npart),
  fNPartSideA(npartSideA),
  fNPartSideC(npartSideC),
  fImpParameter(b),
  fImpParSideA(bSideA),
  fImpParSideC(bSideC),
  fRecoFlag(recoFlag),
  fEnergyFlag(energyFlag),
  fIsScalerOn(scalerOn)
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
  for(Int_t j=0; j<32; j++){
    fZDCScaler[j] = scaler[j];
    for(Int_t y=0; y<4; y++) fZDCTDCData[j][y] = tdcData[j][y];
  }
}

//______________________________________________________________________________
AliZDCReco::AliZDCReco(const AliZDCReco &oldreco) :
TObject(),
fNDetSpecNSideA(oldreco.GetNDetSpecNSideA()),
fNDetSpecPSideA(oldreco.GetNDetSpecPSideA()),
fNDetSpecNSideC(oldreco.GetNDetSpecNSideC()),        
fNDetSpecPSideC(oldreco.GetNDetSpecPSideC()),       
fNTrueSpectators(oldreco.GetNTrueSpectators()),
fNTrueSpecSideA(oldreco.GetNTrueSpecSideA()),
fNTrueSpecSideC(oldreco.GetNTrueSpecSideC()),  	
fNParticipants(oldreco.GetNParticipants()),		       
fNPartSideA(oldreco.GetNPartSideA()),		       
fNPartSideC(oldreco.GetNPartSideC()),  		       
fImpParameter(oldreco.GetImpParameter()),      
fImpParSideA(oldreco.GetImpParSideA()),      
fImpParSideC(oldreco.GetImpParSideC()),
fRecoFlag(oldreco.GetRecoFlag()),
fEnergyFlag(oldreco.GetEnergyFlag()),
fIsScalerOn(oldreco.IsScalerOn())    
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
  for(Int_t j=0; j<32; j++){
    fZDCScaler[j] = oldreco.GetZDCScaler(j);
    for(Int_t y=0; y<4; y++) fZDCTDCData[j][y] = oldreco.GetZDCTDCData(j, y);
  }
}

//______________________________________________________________________________
void AliZDCReco::Print(Option_t *) const {
  //
  // Printing Reconstruction Parameters
  //
  printf(" ****************** AliZDCReco object ******************\n"
  	 "       ---------------   side A ---------------\n"
	 " E_ZN = %1.2f TeV, E_ZP = %1.2f TeV, "
	 " E_ZEM1 =  %1.2f TeV,  E_ZEM2 = %1.2f TeV\n "
	 " N_spec_n = %d, N_spec_p = %d,"
	 " N_part = %d, b = %1.4f fm\n"
  	 "       ---------------   side C ---------------\n"
	 " E_ZN = %1.2f TeV, E_ZP = %1.2f TeV, "
	 " N_spec_n = %d, N_spec_p = %d,"
	 " N_part = %d, b = %1.4f fm\n"
         " *******************************************************\n",
	 fZN2Energy[0]/1000., fZP2Energy[0]/1000.,
	 fZEM1signal[0]/1000.,fZEM2signal[0]/1000.,
	 fNDetSpecNSideA,fNDetSpecPSideA, fNPartSideA,fImpParSideA,
	 fZN1Energy[0]/1000.,fZP1Energy[0]/1000.,	
         fNDetSpecNSideC,fNDetSpecPSideC,fNPartSideC,fImpParSideC);
	 
}
