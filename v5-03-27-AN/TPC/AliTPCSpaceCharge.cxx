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

// _________________________________________________________________
//
// Begin_Html
//   <h2> AliTPCSpaceCharge class  </h2>    
//   The class calculates the space point distortions due to a rotational 
//   symmetric space charge distribution with the TPC drift volume. 
//   <p>
//   The class uses the PoissonRelaxation2D to calculate the resulting 
//   electrical field inhomogeneities in the (r,z)-plane. Then, the 
//   Langevin-integral formalism is used to calculate the space point distortions. 
//   <p>
//   The class assumes, that the distortions scales linearly with the magnitude 
//   of the space charge distribution $\rho(r,z)$. The in here assumed distribution is 
//   $$\rho(r,z) = \frac{(A-B\,z)}{r^2} $$ wherein the factors A and B scale with the
//   event multiplicity and the interaction rate.
//   <p>
//   The scaling factor can be set via the function SetCorrectionFactor. An example of 
//   the shape of the distortions is given below.
// End_Html
//
// Begin_Macro(source)
//   {
//   gROOT->SetStyle("Plain"); gStyle->SetPalette(1);
//   TCanvas *c2 = new TCanvas("cAliTPCSpaceCharge","cAliTPCSpaceCharge",500,300); 
//   AliTPCSpaceCharge sc;
//   sc.SetOmegaTauT1T2(-0.32,1,1); // B=0.5 Tesla
//   sc.SetCorrectionFactor(0.0015);
//   sc.CreateHistoDRinZR(0.)->Draw("surf2");
//   return c2;
//   } 
// End_Macro
//
// Begin_Html
//   <p>
//   Date: 23/08/2010 <br>
// Authors: Jim Thomas, Stefan Rossegger 
// End_Html 
// _________________________________________________________________



#include "AliMagF.h"
#include "TGeoGlobalMagField.h"
#include "AliTPCcalibDB.h"
#include "AliTPCParam.h"
#include "AliLog.h"
#include "TMatrixD.h"

#include "TMath.h"
#include "AliTPCROC.h"
#include "AliTPCSpaceCharge.h"

ClassImp(AliTPCSpaceCharge)

AliTPCSpaceCharge::AliTPCSpaceCharge()
  : AliTPCCorrection("SpaceCharge2D","Space Charge 2D"),
    fC0(0.),fC1(0.),fCorrectionFactor(0.001),
    fInitLookUp(kFALSE)
{
  //
  // default constructor
  //
 
}

AliTPCSpaceCharge::~AliTPCSpaceCharge() {
  //
  // default destructor
  //
}



void AliTPCSpaceCharge::Init() {
  //
  // Initialization funtion
  //
  
  AliMagF* magF= (AliMagF*)TGeoGlobalMagField::Instance()->GetField();
  if (!magF) AliError("Magneticd field - not initialized");
  Double_t bzField = magF->SolenoidField()/10.; //field in T
  AliTPCParam *param= AliTPCcalibDB::Instance()->GetParameters();
  if (!param) AliError("Parameters - not initialized");
  Double_t vdrift = param->GetDriftV()/1000000.; // [cm/us]   // From dataBase: to be updated: per second (ideally)
  Double_t ezField = 400; // [V/cm]   // to be updated: never (hopefully)
  Double_t wt = -10.0 * (bzField*10) * vdrift / ezField ; 
  // Correction Terms for effective omegaTau; obtained by a laser calibration run
  SetOmegaTauT1T2(wt,fT1,fT2);

  InitSpaceChargeDistortion(); // fill the look up table
}

void AliTPCSpaceCharge::Update(const TTimeStamp &/*timeStamp*/) {
  //
  // Update function 
  //
  AliMagF* magF= (AliMagF*)TGeoGlobalMagField::Instance()->GetField();
  if (!magF) AliError("Magneticd field - not initialized");
  Double_t bzField = magF->SolenoidField()/10.; //field in T
  AliTPCParam *param= AliTPCcalibDB::Instance()->GetParameters();
  if (!param) AliError("Parameters - not initialized");
  Double_t vdrift = param->GetDriftV()/1000000.; // [cm/us]  // From dataBase: to be updated: per second (ideally)
  Double_t ezField = 400; // [V/cm]   // to be updated: never (hopefully)
  Double_t wt = -10.0 * (bzField*10) * vdrift / ezField ; 
  // Correction Terms for effective omegaTau; obtained by a laser calibration run
  SetOmegaTauT1T2(wt,fT1,fT2);

  //  SetCorrectionFactor(1.); // should come from some database

}



void AliTPCSpaceCharge::GetCorrection(const Float_t x[],const Short_t roc,Float_t dx[]) {
  //
  // Calculates the correction due the Space Charge effect within the TPC drift volume
  //   

  if (!fInitLookUp) {
    AliInfo("Lookup table was not initialized! Perform the inizialisation now ...");
    InitSpaceChargeDistortion();
  }
  Int_t   order     = 1 ;    // FIXME: hardcoded? Linear interpolation = 1, Quadratic = 2         
                        
  Double_t intEr, intEphi, intdEz;
  Double_t r, phi, z ;
  Int_t    sign;

  r      =  TMath::Sqrt( x[0]*x[0] + x[1]*x[1] ) ;
  phi    =  TMath::ATan2(x[1],x[0]) ;
  if ( phi < 0 ) phi += TMath::TwoPi() ;                   // Table uses phi from 0 to 2*Pi
  z      =  x[2] ;                                         // Create temporary copy of x[2]

  if ( (roc%36) < 18 ) {
    sign =  1;       // (TPC A side)
  } else {
    sign = -1;       // (TPC C side)
  }
  
  if ( sign==1  && z <  fgkZOffSet ) z =  fgkZOffSet;    // Protect against discontinuity at CE
  if ( sign==-1 && z > -fgkZOffSet ) z = -fgkZOffSet;    // Protect against discontinuity at CE
  

  if ( (sign==1 && z<0) || (sign==-1 && z>0) ) // just a consistency check
    AliError("ROC number does not correspond to z coordinate! Calculation of distortions is most likely wrong!");

  // Efield is symmetric in phi - 2D calculation
  intEphi = 0.0; 
  // Get the E field integrals
  Interpolate2DEdistortion( order, r, z, fLookUpErOverEz, intEr );
  // Get DeltaEz field integral
  Interpolate2DEdistortion( order, r, z, fLookUpDeltaEz, intdEz );
  
 
  // Calculate distorted position
  if ( r > 0.0 ) {
    phi =  phi + fCorrectionFactor *( fC0*intEphi - fC1*intEr ) / r;      
    r   =  r   + fCorrectionFactor *( fC0*intEr   + fC1*intEphi );  
  }
  Double_t dz = intdEz*fCorrectionFactor;
 
  // Calculate correction in cartesian coordinates
  dx[0] = - (r * TMath::Cos(phi) - x[0]);
  dx[1] = - (r * TMath::Sin(phi) - x[1]); 
  dx[2] = - dz;  // z distortion - (internally scaled with driftvelocity dependency 
                 // on the Ez field 

}

void AliTPCSpaceCharge::InitSpaceChargeDistortion() {
  //
  // Initialization of the Lookup table which contains the solutions of the 
  // poisson problem
  //

  const Float_t  gridSizeR   =  (fgkOFCRadius-fgkIFCRadius) / (kRows-1) ;
  const Float_t  gridSizeZ   =  fgkTPCZ0 / (kColumns-1) ;

  TMatrixD voltArray(kRows,kColumns);        // dummy boundary vectors
  TMatrixD chargeDensity(kRows,kColumns);    // charge
  TMatrixD arrayErOverEz(kRows,kColumns);    // solution in Er
  TMatrixD arrayDeltaEz(kRows,kColumns);    // solution in Ez

  Double_t  rList[kRows], zedList[kColumns] ;
  
  // Fill arrays with initial conditions.  V on the boundary and ChargeDensity in the volume.      
  for ( Int_t j = 0 ; j < kColumns ; j++ ) {
    Double_t zed = j*gridSizeZ ;
    zedList[j] = zed ;
    for ( Int_t i = 0 ; i < kRows ; i++ )  {
      Double_t radius = fgkIFCRadius + i*gridSizeR ;
      rList[i]           = radius ;
      voltArray(i,j)        = 0;  // Initialize voltArray to zero - not used in this class
      chargeDensity(i,j)     = 0;  // Initialize ChargeDensity to zero
    }
  }      

  // Fill the initial conditions
  for ( Int_t j = 1 ; j < kColumns-1 ; j++ ) {
    Double_t zed = j*gridSizeZ ;
    for ( Int_t i = 1 ; i < kRows-1 ; i++ ) { 
      Double_t radius = fgkIFCRadius + i*gridSizeR ;

      Double_t zterm = (fgkTPCZ0-zed) * (fgkOFCRadius*fgkOFCRadius - fgkIFCRadius*fgkIFCRadius) / fgkTPCZ0 ;
      // for 1/R**2 charge density in the TPC; then integrated in Z due to drifting ions
      chargeDensity(i,j) = zterm / ( TMath::Log(fgkOFCRadius/fgkIFCRadius) * ( radius*radius ) ) ;              
    }
  }


  // Solve the electrosatic problem in 2D 

  PoissonRelaxation2D( voltArray, chargeDensity, arrayErOverEz, arrayDeltaEz, kRows, kColumns, kIterations ) ;
  
  //Interpolate results onto standard grid for Electric Fields
  Int_t ilow=0, jlow=0 ;
  Double_t z,r;
  Float_t saveEr[2], saveEz[2] ;	      
  for ( Int_t i = 0 ; i < kNZ ; ++i )  {
    z = TMath::Abs( fgkZList[i] ) ; // assume symmetric behaviour on A and C side
    for ( Int_t j = 0 ; j < kNR ; ++j ) {

      // Linear interpolation !!
      r = fgkRList[j] ;
      Search( kRows,   rList, r, ilow ) ;          // Note switch - R in rows and Z in columns
      Search( kColumns, zedList, z, jlow ) ;
      if ( ilow < 0 ) ilow = 0 ;                   // check if out of range
      if ( jlow < 0 ) jlow = 0 ;   
      if ( ilow + 1  >=  kRows - 1 ) ilow =  kRows - 2 ;	      
      if ( jlow + 1  >=  kColumns - 1 ) jlow =  kColumns - 2 ; 
    
      saveEr[0] = arrayErOverEz(ilow,jlow) + 
	(arrayErOverEz(ilow,jlow+1)-arrayErOverEz(ilow,jlow))*(z-zedList[jlow])/gridSizeZ ;
      saveEr[1] = arrayErOverEz(ilow+1,jlow) + 
	(arrayErOverEz(ilow+1,jlow+1)-arrayErOverEz(ilow+1,jlow))*(z-zedList[jlow])/gridSizeZ ;
      saveEz[0] = arrayDeltaEz(ilow,jlow) + 
	(arrayDeltaEz(ilow,jlow+1)-arrayDeltaEz(ilow,jlow))*(z-zedList[jlow])/gridSizeZ ;
      saveEz[1] = arrayDeltaEz(ilow+1,jlow) + 
	(arrayDeltaEz(ilow+1,jlow+1)-arrayDeltaEz(ilow+1,jlow))*(z-zedList[jlow])/gridSizeZ ;

      
      fLookUpErOverEz[i][j] = saveEr[0] + (saveEr[1]-saveEr[0])*(r-rList[ilow])/gridSizeR ;
      fLookUpDeltaEz[i][j]  = saveEz[0] + (saveEz[1]-saveEz[0])*(r-rList[ilow])/gridSizeR ;

      if (fgkZList[i]<0)  fLookUpDeltaEz[i][j] *= -1; // C side is negative z
    }
  }
  
  fInitLookUp = kTRUE;

}

void AliTPCSpaceCharge::Print(const Option_t* option) const {
  //
  // Print function to check the settings of the boundary vectors
  // option=="a" prints the C0 and C1 coefficents for calibration purposes
  //

  TString opt = option; opt.ToLower();
  printf("%s\n",GetTitle());
  printf(" - Space Charge effects assuming a radial symmetric z over r^2 SC-distribution.\n");
  printf("   SC correction factor: %f \n",fCorrectionFactor);

  if (opt.Contains("a")) { // Print all details
    printf(" - T1: %1.4f, T2: %1.4f \n",fT1,fT2);
    printf(" - C1: %1.4f, C0: %1.4f \n",fC1,fC0);
  } 
   
  if (!fInitLookUp) AliError("Lookup table was not initialized! You should do InitSpaceChargeDistortion() ...");

}
