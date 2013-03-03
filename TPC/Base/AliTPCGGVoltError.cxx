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

////////////////////////////////////////////////////////////////////////////
// AliTPCGGVoltError class                                                //
////////////////////////////////////////////////////////////////////////////


#include "AliMagF.h"
#include "TGeoGlobalMagField.h"
#include "AliTPCcalibDB.h"
#include "AliTPCParam.h"
#include "AliLog.h"

#include "AliTPCGGVoltError.h"
#include <TMath.h>

AliTPCGGVoltError::AliTPCGGVoltError()
  : AliTPCCorrection("GGVoltError","GatingGrid (GG) Voltage Error"),
    fC0(0.),fC1(0.),
    fDeltaVGGA(0.),fDeltaVGGC(0.),
    fInitLookUp(kFALSE)
{
  //
  // default constructor
  //
}

AliTPCGGVoltError::~AliTPCGGVoltError() {
  //
  // default destructor
  //
}

void AliTPCGGVoltError::Init() {
  //
  // Init function
  //
  AliMagF* magF= (AliMagF*)TGeoGlobalMagField::Instance()->GetField();
  if (!magF) AliError("Magneticd field - not initialized");
  Double_t bzField = magF->SolenoidField()/10.; //field in T
  AliTPCParam *param= AliTPCcalibDB::Instance()->GetParameters();
  if (!param) AliError("Parameters - not initialized");
  Double_t vdrift = param->GetDriftV()/1000000.; // [cm/us]   // From dataBase: to be updated: per second (ideally)
  Double_t ezField = 400; // [V/cm]   // to be updated: never (hopefully)
  Double_t wt = -10.0 * (bzField*10) * vdrift / ezField ; 
  //
  SetOmegaTauT1T2(wt,fT1,fT2);
  InitGGVoltErrorDistortion();
  //SetDeltaVGGA(0.0);//  ideally from the database
  //SetDeltaVGGC(0.0);//  ideally from the database
}

void AliTPCGGVoltError::Update(const TTimeStamp &/*timeStamp*/) {
  //
  // Update function
  //
  AliMagF* magF= (AliMagF*)TGeoGlobalMagField::Instance()->GetField();
  if (!magF) AliError("Magneticd field - not initialized");
  Double_t bzField = magF->SolenoidField()/10.; //field in T
  AliTPCParam *param= AliTPCcalibDB::Instance()->GetParameters();
  if (!param) AliError("Parameters - not initialized");
  Double_t vdrift = param->GetDriftV()/1000000.; // [cm/us]   // From dataBase: to be updated: per second (ideally)
  Double_t ezField = 400; // [V/cm]   // to be updated: never (hopefully)
  Double_t wt = -10.0 * (bzField*10) * vdrift / ezField ; 

  SetOmegaTauT1T2(wt,fT1,fT2);
  //  InitGGVoltErrorDistortion(); // not necessary in here since the Voltage should not change!
}



void AliTPCGGVoltError::GetCorrection(const Float_t x[],const Short_t roc,Float_t dx[]) {

  //
  // Gated Grid Voltage Error
  //
  // Calculates the effect of having an incorrect voltage on the A or C end plate Gated Grids.
  //
  // Electrostatic Equations from StarNote SN0253 by Howard Wieman.
  //
  
  if (!fInitLookUp) AliError("Lookup table was not initialized! You should do InitGGVoltErrorDistortion() ...");
  
  Int_t   order     = 1 ;               // FIXME: hardcoded? Linear interpolation = 1, Quadratic = 2         
 
  Double_t intEr, intEphi ;
  Double_t r, phi, z ;
  Int_t    sign ;

  Double_t deltaVGG;
  
  r   = TMath::Sqrt( x[0]*x[0] + x[1]*x[1] );
  phi = TMath::ATan2(x[1],x[0]);
  if ( phi < 0 ) phi += TMath::TwoPi();                   // Table uses phi from 0 to 2*Pi
  z   = x[2] ;

  if ( (roc%36) < 18 ) {
    sign =  1; 
    deltaVGG = fDeltaVGGA;           // (TPC End A)
  } else {
    sign = -1;                       // (TPC End C)
    deltaVGG = fDeltaVGGC; 
  }

  if ( sign==1  && z <  fgkZOffSet ) z =  fgkZOffSet;    // Protect against discontinuity at CE
  if ( sign==-1 && z > -fgkZOffSet ) z = -fgkZOffSet;    // Protect against discontinuity at CE

  Interpolate2DEdistortion( order, r, z, fGGVoltErrorER, intEr );
  intEphi = 0.0;  // Efield is symmetric in phi

  // Calculate distorted position
  if ( r > 0.0 ) {
    phi =  phi + deltaVGG*( fC0*intEphi - fC1*intEr ) / r;      
    r   =  r   + deltaVGG*( fC0*intEr   + fC1*intEphi );  
  }
  
  // Calculate correction in cartesian coordinates
  dx[0] = r * TMath::Cos(phi) - x[0];
  dx[1] = r * TMath::Sin(phi) - x[1]; 
  dx[2] = 0.; // z distortion not implemented (1st order distortions) - see e.g. AliTPCBoundaryVoltError-class



}


Float_t AliTPCGGVoltError::GetIntErOverEz(const Float_t x[],const Short_t roc) {
  //
  // This function is purely for calibration purposes
  // Calculates the integral (int Er/Ez dz) for the setted GG voltage offset 
  // 

  if (!fInitLookUp) AliError("Lookup table was not initialized! You should do InitGGVoltErrorDistortion() ...");

  Int_t   order     = 1 ;     // FIXME: so far hardcoded? Linear interpolation = 1, Quadratic = 2         
  
  Double_t intEr;
  Double_t r, phi, z ;
  Int_t    sign ;
  
  Double_t deltaVGG;
  
  r   = TMath::Sqrt( x[0]*x[0] + x[1]*x[1] );
  phi = TMath::ATan2(x[1],x[0]);
  if ( phi < 0 ) phi += TMath::TwoPi();        // Table uses phi from 0 to 2*Pi
  z   = x[2] ;

  if ( (roc%36) < 18 ) {
    sign =  1; 
    deltaVGG = fDeltaVGGA;           // (TPC End A)
  } else {
    sign = -1;                       // (TPC End C)
    deltaVGG = fDeltaVGGC; 
  }

  if ( sign==1  && z <  fgkZOffSet ) z =  fgkZOffSet;    // Protect against discontinuity at CE
  if ( sign==-1 && z > -fgkZOffSet ) z = -fgkZOffSet;    // Protect against discontinuity at CE

  Interpolate2DEdistortion(order, r, z, fGGVoltErrorER, intEr );

  return (intEr*deltaVGG);

}

void AliTPCGGVoltError::InitGGVoltErrorDistortion() {
  //
  // Initialization of the Lookup table which contains the solutions of the GG Error problem
  //

  Double_t r,z;
  Int_t nterms = 100 ;
  for ( Int_t i = 0 ; i < kNZ ; ++i ) {
    z = fgkZList[i] ;
    for ( Int_t j = 0 ; j < kNR ; ++j ) {
      r = fgkRList[j] ;
      fGGVoltErrorER[i][j] = 0.0 ; 	    
      Double_t intz = 0.0 ;
      for ( Int_t n = 1 ; n < nterms ; ++n ) {
	Double_t k    =  n * TMath::Pi() / fgkTPCZ0 ;
	Double_t ein  =  0 ;                    // Error potential on the IFC
	Double_t eout =  0 ;                    // Error potential on the OFC
	if ( z < 0 ) {
	  ein   =  -2.0 / ( k * (fgkCathodeV - fgkGG) ) ;       
	  eout  =  -2.0 / ( k * (fgkCathodeV - fgkGG) ) ;       
	}
	if ( z == 0 ) continue ;
	if ( z > 0 ) {
	  ein   =  -2.0 / ( k * (fgkCathodeV - fgkGG) ) ;       
	  eout  =  -2.0 / ( k * (fgkCathodeV - fgkGG) ) ;       
	}
	Double_t an   =  ein  * TMath::BesselK0( k*fgkOFCRadius ) - eout * TMath::BesselK0( k*fgkIFCRadius ) ;
	Double_t bn   =  eout * TMath::BesselI0( k*fgkIFCRadius ) - ein  * TMath::BesselI0( k*fgkOFCRadius ) ;
	Double_t numerator =
	  an * TMath::BesselI1( k*r ) - bn * TMath::BesselK1( k*r ) ;
	Double_t denominator =
	  TMath::BesselK0( k*fgkOFCRadius ) * TMath::BesselI0( k*fgkIFCRadius ) -
	  TMath::BesselK0( k*fgkIFCRadius ) * TMath::BesselI0( k*fgkOFCRadius ) ;
	Double_t zterm = TMath::Cos( k*(fgkTPCZ0-TMath::Abs(z)) ) - 1 ;
	intz += zterm * numerator / denominator ;
	// Assume series converges, break if small terms
	if ( n>10 && TMath::Abs(intz)*1.e-10 > TMath::Abs(numerator/denominator) ) break;   
      }
      fGGVoltErrorER[i][j] = (Double_t) intz ;

    }
  }
  
  fInitLookUp = kTRUE;
}



void AliTPCGGVoltError::Print(const Option_t* option) const {
  //
  // Print function to check the settings (e.g. voltage offsets)
  // option=="a" prints the C0 and C1 coefficents for calibration purposes
  //

  TString opt = option; opt.ToLower();
  printf("%s\n",GetTitle());
  printf(" - GG Voltage offset: A-side: %3.1f V, C-side: %3.1f V \n",fDeltaVGGA,fDeltaVGGC);  
  if (opt.Contains("a")) { // Print all details
    printf(" - T1: %1.4f, T2: %1.4f \n",fT1,fT2);
    printf(" - C1: %1.4f, C0: %1.4f \n",fC1,fC0);
  }    

  if (!fInitLookUp) AliError("Lookup table was not initialized! You should do InitGGVoltErrorDistortion() ...");

  
}
