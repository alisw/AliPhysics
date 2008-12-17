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
AliMagFC::AliMagFC()
    :AliMagF(),
    fCompensator(kFALSE),
    fBeamType(kBeamTypepp),
    fBeamEnergy(0),
    fQuadGradient(0),
    fDipoleField(0),
    fCCorrField(0), 
    fACorr1Field(0),
    fACorr2Field(0)
{
  // 
  // Default constructor
  //
}

//________________________________________
AliMagFC::AliMagFC(const char *name, const char *title, Int_t integ, 
		   Float_t factor, Float_t fmax)
    : AliMagF(name,title,integ,factor,fmax),
    fCompensator(kFALSE),
    fBeamType(kBeamTypepp), 
    fBeamEnergy(7000.),
    fQuadGradient(0),
    fDipoleField(0),
    fCCorrField(0), 
    fACorr1Field(0),
    fACorr2Field(0)

{
  // 
  // Standard constructor
  //
  fType = kConst;
  fMap  = 1;
  
}

//________________________________________
void AliMagFC::Field(const float *x, float *b) const
{
  //
  // Method to return the field in a point
  //
  b[0]=b[1]=b[2]=0;
  if(fMap==1) {
    if(TMath::Abs(x[2])<700 && x[0]*x[0]+(x[1]+30)*(x[1]+30) < 560*560) {
      b[2]=2;
    } 
    else {
      if(-725 >= x[2] && x[2] >= -1225 ){
	Float_t dz = TMath::Abs(-975-x[2])*0.01;
	b[0] = - (1-0.1*dz*dz)*7;
	if(fFactor!=1) {
	    b[0]*=fFactor;
	    b[1]*=fFactor;
	    b[2]*=fFactor;
	}
      }
      else {
	  ZDCField(x, b);
      }
    }

  } 
  else {
      AliFatal(Form("Invalid field map for constant field %d",fMap));
  }
}
//________________________________________
void AliMagFC::Field(const double *x, double *b) const
{
  //
  // Method to return the field in a point
  //
  b[0]=b[1]=b[2]=0;
  if(fMap==1) {
    if(TMath::Abs(x[2])<700 && x[0]*x[0]+(x[1]+30)*(x[1]+30) < 560*560) {
      b[2]=2;
    } 
    else {
      if(-725 >= x[2] && x[2] >= -1225 ){
	Float_t dz = TMath::Abs(-975-x[2])*0.01;
	b[0] = - (1-0.1*dz*dz)*7;
	if(fFactor!=1) {
	    b[0]*=fFactor;
	    b[1]*=fFactor;
	    b[2]*=fFactor;
	}
      }
      else {
	  ZDCField(x, b);
      }
    }

  } 
  else {
      AliFatal(Form("Invalid field map for constant field %d",fMap));
  }
}

//___________________________________________________
void AliMagFC::ZDCField(const float *x, float *b) const
{
  // ---- This is the ZDC part
  
  float rad2 = x[0] * x[0] + x[1] * x[1];
  static Bool_t init = kFALSE;

  if (! init) {
      init = kTRUE;
      //////////////////////////////////////////////////////////////////////
      // ---- Magnetic field values (according to beam type and energy) ----
      if(fBeamType==kBeamTypepp && fBeamEnergy == 5000.){
	  // p-p @ 5+5 TeV
	  fQuadGradient = 15.7145;
	  fDipoleField  = 27.0558;
	  // SIDE C
	  fCCorrField   = 9.7017;
	  // SIDE A
	  fACorr1Field  = -13.2143;
	  fACorr2Field  = -11.9909;
      } else if (fBeamType == kBeamTypepp && fBeamEnergy == 450.) {
	  // p-p 0.45+0.45 TeV
	  Float_t const kEnergyRatio = fBeamEnergy / 7000.;
	  
	  fQuadGradient = 22.0002 * kEnergyRatio;
	  fDipoleField  = 37.8781 * kEnergyRatio;
	  // SIDE C
	  fCCorrField   =  9.6908;
	  // SIDE A
	  fACorr1Field  = -13.2014;
	  fACorr2Field  = -9.6908;
      } else if ((fBeamType == kBeamTypepp && fBeamEnergy == 7000.) ||
		 (fBeamType == kBeamTypeAA))
      {
	  // Pb-Pb @ 2.7+2.7 TeV or p-p @ 7+7 TeV
	  fQuadGradient = 22.0002;
	  fDipoleField  = 37.8781;
	  // SIDE C
	  fCCorrField   = 9.6908;
	  // SIDE A
	  fACorr1Field  = -13.2014;
	  fACorr2Field  = -9.6908;
      }
  }
  
  
  // SIDE C **************************************************
  if(x[2]<0.){  
    if(x[2] < kCCorrBegin && x[2] > kCCorrEnd && rad2 < kCCorrSqRadius){
	if (fFactor != 0.) {
	    b[0] = fCCorrField;
	    b[1] = 0.;
	    b[2] = 0.;
	} 
    }
    else if(x[2] < kCQ1Begin && x[2] > kCQ1End && rad2 < kCQ1SqRadius){
	b[0] = fQuadGradient*x[1];
	b[1] = fQuadGradient*x[0];
	b[2] = 0.;
    }
    else if(x[2] < kCQ2Begin && x[2] > kCQ2End && rad2 < kCQ2SqRadius){
	b[0] = -fQuadGradient*x[1];
	b[1] = -fQuadGradient*x[0];
	b[2] = 0.;
    }
    else if(x[2] < kCQ3Begin && x[2] > kCQ3End && rad2 < kCQ3SqRadius){
	b[0] = -fQuadGradient*x[1];
	b[1] = -fQuadGradient*x[0];
	b[2] = 0.;
    }
    else if(x[2] < kCQ4Begin && x[2] > kCQ4End && rad2 < kCQ4SqRadius){
	b[0] = fQuadGradient*x[1];
	b[1] = fQuadGradient*x[0];
	b[2] = 0.;
    }
    else if(x[2] < kCD1Begin && x[2] > kCD1End && rad2 < kCD1SqRadius){
	b[1] = fDipoleField;
	b[2] = 0.;
	b[2] = 0.;
    }
    else if(x[2] < kCD2Begin && x[2] > kCD2End){
	if(((x[0]-kCD2XCentre1)*(x[0]-kCD2XCentre1)+(x[1]*x[1]))<kCD2SqRadius
	   || ((x[0]-kCD2XCentre2)*(x[0]-kCD2XCentre2)+(x[1]*x[1]))<kCD2SqRadius){
	  b[1] = -fDipoleField;
	  b[2] = 0.;
	  b[2] = 0.;
	}
    }
  }
  
  // SIDE A **************************************************
  else{        
    if(fCompensator && (x[2] > kACorr1Begin && x[2] < kACorr1End) && rad2 < kCCorr1SqRadius) {
      // Compensator magnet at z = 1075 m 
	if (fFactor != 0.) {
	    b[0] = fACorr1Field;
	    b[1] = 0.;
	    b[2] = 0.;
	}
	return;
    }
    
    if(x[2] > kACorr2Begin && x[2] < kACorr2End && rad2 < kCCorr2SqRadius){
	if (fFactor != 0.) {
	    b[0] = fACorr2Field;
	    b[1] = 0.;
	    b[2] = 0.;
	}
    }          
    else if(x[2] > kAQ1Begin && x[2] < kAQ1End && rad2 < kAQ1SqRadius){
	// First quadrupole of inner triplet de-focussing in x-direction
	b[0] = -fQuadGradient*x[1];
	b[1] = -fQuadGradient*x[0];
	b[2] = 0.;
    }
    else if(x[2] > kAQ2Begin && x[2] < kAQ2End && rad2 < kAQ2SqRadius){
	b[0] = fQuadGradient*x[1];
	b[1] = fQuadGradient*x[0];
	b[2] = 0.;
    }
    else if(x[2] > kAQ3Begin && x[2] < kAQ3End && rad2 < kAQ3SqRadius){
	b[0] = fQuadGradient*x[1];
	b[1] = fQuadGradient*x[0];
	b[2] = 0.;
    }
    else if(x[2] > kAQ4Begin && x[2] < kAQ4End && rad2 < kAQ4SqRadius){
	b[0] = -fQuadGradient*x[1];
	b[1] = -fQuadGradient*x[0];
	b[2] = 0.;
    }
    else if(x[2] > kAD1Begin && x[2] < kAD1End && rad2 < kAD1SqRadius){
	b[0] = 0.;
	b[1] = -fDipoleField;
	b[2] = 0.;
    }
    else if(x[2] > kAD2Begin && x[2] < kAD2End){
	if(((x[0]-kAD2XCentre1)*(x[0]-kAD2XCentre1)+(x[1]*x[1])) < kAD2SqRadius
	   || ((x[0]-kAD2XCentre2)*(x[0]-kAD2XCentre2)+(x[1]*x[1])) < kAD2SqRadius){
	    b[1] = fDipoleField;
	}
    }
  }
}

void AliMagFC::ZDCField(const double *x, double *b) const
{
  // ---- This is the ZDC part
  
  double rad2 = x[0] * x[0] + x[1] * x[1];
  static Bool_t init = kFALSE;

  if (! init) {
      init = kTRUE;
      //////////////////////////////////////////////////////////////////////
      // ---- Magnetic field values (according to beam type and energy) ----
      if(fBeamType==kBeamTypepp && fBeamEnergy == 5000.){
	  // p-p @ 5+5 TeV
	  fQuadGradient = 15.7145;
	  fDipoleField  = 27.0558;
	  // SIDE C
	  fCCorrField   = 9.7017;
	  // SIDE A
	  fACorr1Field  = -13.2143;
	  fACorr2Field  = -11.9909;
      } else if (fBeamType == kBeamTypepp && fBeamEnergy == 450.) {
	  // p-p 0.45+0.45 TeV
	  Float_t const kEnergyRatio = fBeamEnergy / 7000.;
	  
	  fQuadGradient = 22.0002 * kEnergyRatio;
	  fDipoleField  = 37.8781 * kEnergyRatio;
	  // SIDE C
	  fCCorrField   =  9.6908;
	  // SIDE A
	  fACorr1Field  = -13.2014;
	  fACorr2Field  = -9.6908;
      } else if ((fBeamType == kBeamTypepp && fBeamEnergy == 7000.) ||
		 (fBeamType == kBeamTypeAA))
      {
	  // Pb-Pb @ 2.7+2.7 TeV or p-p @ 7+7 TeV
	  fQuadGradient = 22.0002;
	  fDipoleField  = 37.8781;
	  // SIDE C
	  fCCorrField   = 9.6908;
	  // SIDE A
	  fACorr1Field  = -13.2014;
	  fACorr2Field  = -9.6908;
      }
  }
  
  
  // SIDE C **************************************************
  if(x[2]<0.){  
    if(x[2] < kCCorrBegin && x[2] > kCCorrEnd && rad2 < kCCorrSqRadius){
	if (fFactor != 0.) {
	    b[0] = fCCorrField;
	    b[1] = 0.;
	    b[2] = 0.;
	} 
    }
    else if(x[2] < kCQ1Begin && x[2] > kCQ1End && rad2 < kCQ1SqRadius){
	b[0] = fQuadGradient*x[1];
	b[1] = fQuadGradient*x[0];
	b[2] = 0.;
    }
    else if(x[2] < kCQ2Begin && x[2] > kCQ2End && rad2 < kCQ2SqRadius){
	b[0] = -fQuadGradient*x[1];
	b[1] = -fQuadGradient*x[0];
	b[2] = 0.;
    }
    else if(x[2] < kCQ3Begin && x[2] > kCQ3End && rad2 < kCQ3SqRadius){
	b[0] = -fQuadGradient*x[1];
	b[1] = -fQuadGradient*x[0];
	b[2] = 0.;
    }
    else if(x[2] < kCQ4Begin && x[2] > kCQ4End && rad2 < kCQ4SqRadius){
	b[0] = fQuadGradient*x[1];
	b[1] = fQuadGradient*x[0];
	b[2] = 0.;
    }
    else if(x[2] < kCD1Begin && x[2] > kCD1End && rad2 < kCD1SqRadius){
	b[1] = fDipoleField;
	b[2] = 0.;
	b[2] = 0.;
    }
    else if(x[2] < kCD2Begin && x[2] > kCD2End){
	if(((x[0]-kCD2XCentre1)*(x[0]-kCD2XCentre1)+(x[1]*x[1]))<kCD2SqRadius
	   || ((x[0]-kCD2XCentre2)*(x[0]-kCD2XCentre2)+(x[1]*x[1]))<kCD2SqRadius){
	  b[1] = -fDipoleField;
	  b[2] = 0.;
	  b[2] = 0.;
	}
    }
  }
  
  // SIDE A **************************************************
  else{        
    if(fCompensator && (x[2] > kACorr1Begin && x[2] < kACorr1End) && rad2 < kCCorr1SqRadius) {
      // Compensator magnet at z = 1075 m 
	if (fFactor != 0.) {
	    b[0] = fACorr1Field;
	    b[1] = 0.;
	    b[2] = 0.;
	}
	return;
    }
    
    if(x[2] > kACorr2Begin && x[2] < kACorr2End && rad2 < kCCorr2SqRadius){
	if (fFactor != 0.) {
	    b[0] = fACorr2Field;
	    b[1] = 0.;
	    b[2] = 0.;
	}
    }          
    else if(x[2] > kAQ1Begin && x[2] < kAQ1End && rad2 < kAQ1SqRadius){
	// First quadrupole of inner triplet de-focussing in x-direction
	b[0] = -fQuadGradient*x[1];
	b[1] = -fQuadGradient*x[0];
	b[2] = 0.;
    }
    else if(x[2] > kAQ2Begin && x[2] < kAQ2End && rad2 < kAQ2SqRadius){
	b[0] = fQuadGradient*x[1];
	b[1] = fQuadGradient*x[0];
	b[2] = 0.;
    }
    else if(x[2] > kAQ3Begin && x[2] < kAQ3End && rad2 < kAQ3SqRadius){
	b[0] = fQuadGradient*x[1];
	b[1] = fQuadGradient*x[0];
	b[2] = 0.;
    }
    else if(x[2] > kAQ4Begin && x[2] < kAQ4End && rad2 < kAQ4SqRadius){
	b[0] = -fQuadGradient*x[1];
	b[1] = -fQuadGradient*x[0];
	b[2] = 0.;
    }
    else if(x[2] > kAD1Begin && x[2] < kAD1End && rad2 < kAD1SqRadius){
	b[0] = 0.;
	b[1] = -fDipoleField;
	b[2] = 0.;
    }
    else if(x[2] > kAD2Begin && x[2] < kAD2End){
	if(((x[0]-kAD2XCentre1)*(x[0]-kAD2XCentre1)+(x[1]*x[1])) < kAD2SqRadius
	   || ((x[0]-kAD2XCentre2)*(x[0]-kAD2XCentre2)+(x[1]*x[1])) < kAD2SqRadius){
	    b[1] = fDipoleField;
	}
    }
  }
}

