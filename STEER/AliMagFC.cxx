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

/* $Id$ */

//-------------------------------------------------------------------------
//     Constant magnetic field class
//     Used by AliRun class
//     Author:
//-------------------------------------------------------------------------

#include <stdlib.h>

#include "AliLog.h"
#include "AliMagFC.h"

ClassImp(AliMagFC)

//________________________________________
AliMagFC::AliMagFC(const char *name, const char *title, Int_t integ, 
		   Float_t factor, Float_t fmax)
    : AliMagF(name,title,integ,factor,fmax),
      fCompensator(kFALSE)
{
  // 
  // Standard constructor
  //
  fType = kConst;
  fMap  = 1;
}

//________________________________________
void AliMagFC::Field(Float_t *x, Float_t *b) const
{
  //
  // Method to return the field in a point
  //
  b[0]=b[1]=b[2]=0;
  if(fMap==1) {
    if(TMath::Abs(x[2])<700 && x[0]*x[0]+(x[1]+30)*(x[1]+30) < 560*560) {
      b[2]=2;
    } else {
      if ( -725 >= x[2] && x[2] >= -1225 ) {
	Float_t dz = TMath::Abs(-975-x[2])*0.01;
	b[0] = - (1-0.1*dz*dz)*7;
      }
      else {
	  ZDCField(x, b);
      }
    }
    if(fFactor!=1) {
	b[0]*=fFactor;
	b[1]*=fFactor;
	b[2]*=fFactor;
    }
  } else {
      AliFatal(Form("Invalid field map for constant field %d",fMap));
  }
}


void AliMagFC::ZDCField(Float_t *x, Float_t *b) const
{
//This is the ZDC part
    Float_t rad2 = x[0] * x[0] + x[1] * x[1];

    if (fCompensator && (x[2] > 919. && x[2] < 1231.) && rad2 < 16.) {
      // Compensator magnet at z = 1075 m 
	b[0] = 10.9;
	b[1] = 0.;
	b[2] = 0.;
	return;
    }
    
    
    if(x[2] < kCORBEG2 && x[2] > kCOREND2){
	if(rad2<kCOR2RA2){
	    b[0] = - kFCORN2;
	}
    }
    else if(x[2] < kZ1BEG && x[2] > kZ1END){  
	if(rad2<kZ1RA2){
	    b[0] =  kG1*x[1];
	    b[1] =  kG1*x[0];
	}
    }
    else if(x[2] < kZ2BEG && x[2] > kZ2END){  
	if(rad2<kZ2RA2){
	    b[0] = -kG1*x[1];
	    b[1] = -kG1*x[0];
	}
    }
    else if(x[2] < kZ3BEG && x[2] > kZ3END){  
	if(rad2<kZ3RA2){
	    b[0] = -kG1*x[1];
	    b[1] = -kG1*x[0];
	}
    }
    else if(x[2] < kZ4BEG && x[2] > kZ4END){  
	if(rad2<kZ4RA2){
	    b[0] =  kG1*x[1];
	    b[1] =  kG1*x[0];
	}
    }
    else if(x[2] < kD1BEG && x[2] > kD1END){ 
	if(rad2<kD1RA2){
	    b[1] = -kFDIP;
	}
    }
    else if(x[2] < kD2BEG && x[2] > kD2END){
	if(((x[0]-kXCEN1D2)*(x[0]-kXCEN1D2)+(x[1]-kYCEN1D2)*(x[1]-kYCEN1D2))<kD2RA2
	   || ((x[0]-kXCEN2D2)*(x[0]-kXCEN2D2)+(x[1]-kYCEN2D2)*(x[1]-kYCEN2D2))<kD2RA2){
	    b[1] = kFDIP;
	}
    }
    
// *************************** LEFT LINE ***********************************************
    
    if(x[2] > kCORBEG2l && x[2] < kCOREND2l){
	if(rad2<kCOR2RA2l){
	    b[0] = kFCORN2l;
	}
    }          
    else if(x[2] > kZ1BEGl && x[2] < kZ1ENDl){  
	if(rad2<kZ1RA2l){
	// First quadrupole of inner triplet de-focussing in x-direction
	    b[0] = -kG1l*x[1];
	    b[1] = -kG1l*x[0];
	}
    }
    else if(x[2] > kZ2BEGl && x[2] < kZ2ENDl){  
	if(rad2<kZ2RA2l){
	    b[0] = kG1l*x[1];
	    b[1] = kG1l*x[0];
	}
    }
    else if(x[2] > kZ3BEGl && x[2] < kZ3ENDl){  
	if(rad2<kZ3RA2l){
	    b[0] = kG1l*x[1];
	    b[1] = kG1l*x[0];
	}
    }
    else if(x[2] > kZ4BEGl && x[2] < kZ4ENDl){  
	if(rad2<kZ4RA2l){
	    b[0] = -kG1l*x[1];
	    b[1] = -kG1l*x[0];
	}
    }
    else if(x[2] > kD1BEGl && x[2] < kD1ENDl){ 
	if(rad2<kD1RA2l){
	    b[1] = kFDIPl;
	}
    }
    else if(x[2] > kD2BEGl && x[2] < kD2ENDl){
	if(((x[0]-kXCEN1D2l)*(x[0]-kXCEN1D2l)+(x[1]-kYCEN1D2l)*(x[1]-kYCEN1D2l))<kD2RA2l
	   || ((x[0]-kXCEN2D2l)*(x[0]-kXCEN2D2l)+(x[1]-kYCEN2D2l)*(x[1]-kYCEN2D2l))<kD2RA2l){
	    b[1] = -kFDIPl;
	}
    }
    
}
