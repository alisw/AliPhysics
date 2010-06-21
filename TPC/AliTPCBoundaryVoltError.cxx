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

//////////////////////////////////////////////////////////////////////////////
//                                                                          //
// AliTPCBoundaryVoltError class                                            //
// The class calculates the space point distortions due to residual voltage //
// errors on the main boundaries of the TPC. For example, the inner vessel  //
// of the TPC is shifted by a certain amount, whereas the ROCs on the A side//
// and the ROCs on the C side follow this mechanical shift (at the inner    //
// vessel) in z direction (see example below). This example is commonly     //
// named "conical deformation" of the TPC field cage.                       //
//                                                                          //
// The class allows "effective Omega Tau" corrections.                      // 
//                                                                          //
// NOTE: This class is not  capable of calculation z distortions due to     //
//       drift velocity change in dependence of the electric field!!!       //
//                                                                          //
// date: 01/06/2010                                                         //
// Authors: Jim Thomas, Stefan Rossegger                                    //
//                                                                          //
// Example usage (e.g +1mm shift of "conical deformation")                  //
//  AliTPCBoundaryVoltError bve;                                            //
//  Float_t boundA[8] = {-40,-40,-40,0,0,0,0,-40}; // voltages A-side       //
//  Float_t boundC[6] = { 40, 40, 40,0,0,0}; // voltages C-side             //
//  bve.SetBoundariesA(boundA);                                             //
//  bve.SetBoundariesC(boundC);                                             //
//  bve.SetOmegaTauT1T2(0.32,1.,1.); // values ideally from OCDB            //
//  // initialization of the look up                                        //
//  bve.InitBoundaryVoltErrorDistortion();                                  // 
//  // plot dRPhi distortions ...                                           //
//  bve.CreateHistoDRPhiinZR(1.,100,100)->Draw("surf2");                    //
//////////////////////////////////////////////////////////////////////////////

#include "AliMagF.h"
#include "TGeoGlobalMagField.h"
#include "AliTPCcalibDB.h"
#include "AliTPCParam.h"
#include "AliLog.h"
#include "TMatrixD.h"

#include "TMath.h"
#include "AliTPCROC.h"
#include "AliTPCBoundaryVoltError.h"

ClassImp(AliTPCBoundaryVoltError)

AliTPCBoundaryVoltError::AliTPCBoundaryVoltError()
  : AliTPCCorrection("BoundaryVoltError","Boundary Voltage Error"),
    fC0(0.),fC1(0.),
    fInitLookUp(kFALSE)
{
  //
  // default constructor
  //
  for (Int_t i=0; i<8; i++){
    fBoundariesA[i]= 0;  
    if (i<6) fBoundariesC[i]= 0;
  }
}

AliTPCBoundaryVoltError::~AliTPCBoundaryVoltError() {
  //
  // default destructor
  //
}



void AliTPCBoundaryVoltError::Init() {
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

  InitBoundaryVoltErrorDistortion();
}

void AliTPCBoundaryVoltError::Update(const TTimeStamp &/*timeStamp*/) {
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


}



void AliTPCBoundaryVoltError::GetCorrection(const Float_t x[],const Short_t roc,Float_t dx[]) {
  //
  // Calculates the correction due e.g. residual voltage errors on the TPC boundaries
  //   

  if (!fInitLookUp) AliError("Lookup table was not initialized! You should do InitBoundaryVoltErrorDistortion() ...");

  Int_t   order     = 1 ;               // FIXME: hardcoded? Linear interpolation = 1, Quadratic = 2         
                                        // note that the poisson solution was linearly mirroed on this grid!
  Double_t intEr, intEphi ;
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
  

  intEphi = 0.0;  // Efield is symmetric in phi - 2D calculation

  if ( (sign==1 && z<0) || (sign==-1 && z>0) ) // just a consistency check
    AliError("ROC number does not correspond to z coordinate! Calculation of distortions is most likely wrong!");

  // Get the E field integral
  Interpolate2DEdistortion( order, r, z, fLookUpErOverEz, intEr );
  
  // Calculate distorted position
  if ( r > 0.0 ) {
    phi =  phi + ( fC0*intEphi - fC1*intEr ) / r;      
    r   =  r   + ( fC0*intEr   + fC1*intEphi );  
  }
  
  // Calculate correction in cartesian coordinates
  dx[0] = r * TMath::Cos(phi) - x[0];
  dx[1] = r * TMath::Sin(phi) - x[1]; 
  dx[2] = 0.; // z distortion not implemented (1st order distortions)

}

void AliTPCBoundaryVoltError::InitBoundaryVoltErrorDistortion() {
  //
  // Initialization of the Lookup table which contains the solutions of the 
  // Dirichlet boundary problem
  //

  const Float_t  gridSizeR   =  (fgkOFCRadius-fgkIFCRadius) / (kRows-1) ;
  const Float_t  gridSizeZ   =  fgkTPCZ0 / (kColumns-1) ;

  TMatrixD voltArrayA(kRows,kColumns), voltArrayC(kRows,kColumns); // boundary vectors
  TMatrixD chargeDensity(kRows,kColumns);                              // dummy charge
  TMatrixD arrayErOverEzA(kRows,kColumns), arrayErOverEzC(kRows,kColumns); // solution

  Double_t  rList[kRows], zedList[kColumns] ;
  
  // Fill arrays with initial conditions.  V on the boundary and ChargeDensity in the volume.      
  for ( Int_t j = 0 ; j < kColumns ; j++ ) {
    Double_t zed = j*gridSizeZ ;
    zedList[j] = zed ;
    for ( Int_t i = 0 ; i < kRows ; i++ )  {
      Double_t radius = fgkIFCRadius + i*gridSizeR ;
      rList[i]           = radius ;
      voltArrayA(i,j)        = 0;  // Initialize voltArrayA to zero
      voltArrayC(i,j)        = 0;  // Initialize voltArrayC to zero
      chargeDensity(i,j)     = 0;  // Initialize ChargeDensity to zero - not used in this class
    }
  }      


  // check if boundary values are the same for the C side (for later, saving some calculation time)

  Int_t symmetry = -1; // assume that  A and C side are identical (but anti-symmetric!) // e.g conical deformation
  Int_t sVec[8];

  // check if boundaries are different (regardless the sign)
  for (Int_t i=0; i<8; i++) { 
    if ((TMath::Abs(fBoundariesA[i]) - TMath::Abs(fBoundariesC[i])) > 1e-5) symmetry = 0;  
    sVec[i] = (Int_t) (TMath::Sign((Float_t)1.,fBoundariesA[i])*TMath::Sign((Float_t)1.,fBoundariesC[i])); // == -1 for anti-symmetry
  }
  if (symmetry==-1) { // still the same values?
    // check the kind of symmetry , if even ...
    if (sVec[0]==1 && sVec[1]==1 && sVec[2]==1 && sVec[3]==1 && sVec[4]==1 && sVec[5]==1 && sVec[6]==1 && sVec[7]==1 ) 
      symmetry =  1;
    else if (sVec[0]==-1 && sVec[1]==-1 && sVec[2]==-1 && sVec[3]==-1 && sVec[4]==-1 && sVec[5]==-1 && sVec[6]==-1 && sVec[7]==-1 ) 
      symmetry = -1;
    else
      symmetry =  0; // some of the values differ in the sign -> neither symmetric nor antisymmetric
  }



  // Solve the electrosatic problem in 2D 

  // Fill the complete Boundary vectors
  // Start at IFC at CE and work anti-clockwise through IFC, ROC, OFC, and CE (clockwise for C side)
  for ( Int_t j = 0 ; j < kColumns ; j++ ) {
    Double_t zed = j*gridSizeZ ;
    for ( Int_t i = 0 ; i < kRows ; i++ ) { 
      Double_t radius = fgkIFCRadius + i*gridSizeR ;

      // A side boundary vectors
      if ( i == 0 ) voltArrayA(i,j) += zed   *((fBoundariesA[1]-fBoundariesA[0])/((kColumns-1)*gridSizeZ))
	+ fBoundariesA[0] ; // IFC
      if ( j == kColumns-1 ) voltArrayA(i,j) += radius*((fBoundariesA[3]-fBoundariesA[2])/((kRows-1)*gridSizeR+fgkIFCRadius))
	+ fBoundariesA[2] ; // ROC
      if ( i == kRows-1 ) voltArrayA(i,j) += zed   *((fBoundariesA[4]-fBoundariesA[5])/((kColumns-1)*gridSizeZ))
	+ fBoundariesA[5] ; // OFC
      if ( j == 0 ) voltArrayA(i,j) += radius*((fBoundariesA[6]-fBoundariesA[7])/((kRows-1)*gridSizeR+fgkIFCRadius))
	+ fBoundariesA[7] ; // CE

      if (symmetry==0) {
	// C side boundary vectors
	if ( i == 0 ) voltArrayC(i,j) += zed   *((fBoundariesC[1]-fBoundariesC[0])/((kColumns-1)*gridSizeZ))
	  + fBoundariesC[0] ; // IFC
	if ( j == kColumns-1 ) voltArrayC(i,j) += radius*((fBoundariesC[3]-fBoundariesC[2])/((kRows-1)*gridSizeR+fgkIFCRadius))
	  + fBoundariesC[2] ; // ROC
	if ( i == kRows-1 ) voltArrayC(i,j) += zed   *((fBoundariesC[4]-fBoundariesC[5])/((kColumns-1)*gridSizeZ))
	  + fBoundariesC[5] ; // OFC
	if ( j == 0 ) voltArrayC(i,j) += radius*((fBoundariesC[6]-fBoundariesC[7])/((kRows-1)*gridSizeR+fgkIFCRadius))
	  + fBoundariesC[7] ; // CE

      }
    }
  }

  voltArrayA(0,0)               *= 0.5 ; // Use average boundary condition at corner
  voltArrayA(kRows-1,0)         *= 0.5 ; // Use average boundary condition at corner
  voltArrayA(0,kColumns-1)      *= 0.5 ; // Use average boundary condition at corner
  voltArrayA(kRows-1,kColumns-1)*= 0.5 ; // Use average boundary condition at corner

  if (symmetry==0) {
    voltArrayC(0,0)               *= 0.5 ; // Use average boundary condition at corner
    voltArrayC(kRows-1,0)         *= 0.5 ; // Use average boundary condition at corner
    voltArrayC(0,kColumns-1)      *= 0.5 ; // Use average boundary condition at corner
    voltArrayC(kRows-1,kColumns-1)*= 0.5 ; // Use average boundary condition at corner
  }


  // always solve the problem on the A side
  PoissonRelaxation2D( voltArrayA, chargeDensity, arrayErOverEzA, kRows, kColumns, kIterations ) ;

  if (symmetry!=0) { // A and C side are the same ("anti-symmetric" or "symmetric")
    for ( Int_t j = 0 ; j < kColumns ; j++ ) {
      for ( Int_t i = 0 ; i < kRows ; i++ ) { 
	arrayErOverEzC(i,j) = symmetry*arrayErOverEzA(i,j);
      }
    }
  } else if (symmetry==0) { // A and C side are different - Solve the problem on the C side too
    PoissonRelaxation2D( voltArrayC, chargeDensity, arrayErOverEzC, kRows, kColumns, kIterations ) ;
  }

  //Interpolate results onto standard grid for Electric Fields
  Int_t ilow=0, jlow=0 ;
  Double_t z,r;
  Float_t saveEr[2] ;	      
  for ( Int_t i = 0 ; i < kNZ ; ++i )  {
    z = TMath::Abs( fgkZList[i] ) ;
    for ( Int_t j = 0 ; j < kNR ; ++j ) {
      // Linear interpolation !!
      r = fgkRList[j] ;
	Search( kRows,   rList, r, ilow ) ;          // Note switch - R in rows and Z in columns
	Search( kColumns, zedList, z, jlow ) ;
	if ( ilow < 0 ) ilow = 0 ;                   // check if out of range
	if ( jlow < 0 ) jlow = 0 ;   
	if ( ilow + 1  >=  kRows - 1 ) ilow =  kRows - 2 ;	      
	if ( jlow + 1  >=  kColumns - 1 ) jlow =  kColumns - 2 ; 

	if (fgkZList[i]>0) {         // A side solution
	  saveEr[0] = arrayErOverEzA(ilow,jlow) + 
	    (arrayErOverEzA(ilow,jlow+1)-arrayErOverEzA(ilow,jlow))*(z-zedList[jlow])/gridSizeZ ;
	  saveEr[1] = arrayErOverEzA(ilow+1,jlow) + 
	    (arrayErOverEzA(ilow+1,jlow+1)-arrayErOverEzA(ilow+1,jlow))*(z-zedList[jlow])/gridSizeZ ;
	} else if (fgkZList[i]<0) {  // C side solution
	  saveEr[0] = arrayErOverEzC(ilow,jlow) + 
	    (arrayErOverEzC(ilow,jlow+1)-arrayErOverEzC(ilow,jlow))*(z-zedList[jlow])/gridSizeZ ;
	  saveEr[1] = arrayErOverEzC(ilow+1,jlow) + 
	    (arrayErOverEzC(ilow+1,jlow+1)-arrayErOverEzC(ilow+1,jlow))*(z-zedList[jlow])/gridSizeZ ;
	} else {
	  AliWarning("Field calculation at z=0 (CE) is not allowed!");
	  saveEr[0]=0; saveEr[1]=0;
	}
	fLookUpErOverEz[i][j] = saveEr[0] + (saveEr[1]-saveEr[0])*(r-rList[ilow])/gridSizeR ;
      }
  }
  
  /* delete [] saveEr;
     delete [] sVec;
     delete [] rList;
     delete [] zedList;
  */

  fInitLookUp = kTRUE;

}

void AliTPCBoundaryVoltError::Print(const Option_t* option) const {
  //
  // Print function to check the settings of the boundary vectors
  // option=="a" prints the C0 and C1 coefficents for calibration purposes
  //

  TString opt = option; opt.ToLower();
  printf("%s\n",GetTitle());
  printf(" - Voltage settings (on the TPC boundaries) - linearly interpolated\n");
  printf("  : A-side (anti-clockwise)\n"); 
  printf("     (0,1):\t IFC (CE) : %3.1f V \t IFC (ROC): %3.1f V \n",fBoundariesA[0],fBoundariesA[1]);
  printf("     (2,3):\t ROC (IFC): %3.1f V \t ROC (OFC): %3.1f V \n",fBoundariesA[2],fBoundariesA[3]);
  printf("     (4,5):\t OFC (ROC): %3.1f V \t OFC (CE) : %3.1f V \n",fBoundariesA[4],fBoundariesA[5]);
  printf("     (6,7):\t CE  (OFC): %3.1f V \t CE  (IFC): %3.1f V \n",fBoundariesA[6],fBoundariesA[7]);
  printf("  : C-side (clockwise)\n"); 
  printf("     (0,1):\t IFC (CE) : %3.1f V \t IFC (ROC): %3.1f V \n",fBoundariesC[0],fBoundariesC[1]);
  printf("     (2,3):\t ROC (IFC): %3.1f V \t ROC (OFC): %3.1f V \n",fBoundariesC[2],fBoundariesC[3]);
  printf("     (4,5):\t OFC (ROC): %3.1f V \t OFC (CE) : %3.1f V \n",fBoundariesC[4],fBoundariesC[5]);
  printf("     (6,7):\t CE  (OFC): %3.1f V \t CE  (IFC): %3.1f V \n",fBoundariesC[6],fBoundariesC[7]);

  // Check wether the settings of the Central Electrode agree (on the A and C side)
  // Note: they have to be antisymmetric!
  if (( TMath::Abs(fBoundariesA[6]+fBoundariesC[6])>1e-5) || ( TMath::Abs(fBoundariesA[7]+fBoundariesC[7])>1e-5 ) ){
    AliWarning("Boundary parameters for the Central Electrode (CE) are not anti-symmetric! HOW DID YOU MANAGE THAT?");
    AliWarning("Congratulations, you just splitted the Central Electrode of the TPC!");
    AliWarning("Non-physical settings of the boundary parameter at the Central Electrode");
  }

  if (opt.Contains("a")) { // Print all details
    printf(" - T1: %1.4f, T2: %1.4f \n",fT1,fT2);
    printf(" - C1: %1.4f, C0: %1.4f \n",fC1,fC0);
  } 
   
  if (!fInitLookUp) AliError("Lookup table was not initialized! You should do InitBoundaryVoltErrorDistortion() ...");

}


void AliTPCBoundaryVoltError::SetBoundariesA(Float_t boundariesA[8]){
  //
  // set voltage errors on the TPC boundaries - A side 
  //
  // Start at IFC at the Central electrode and work anti-clockwise (clockwise for C side) through 
  // IFC, ROC, OFC, and CE. The boundary conditions are currently defined to be a linear 
  // interpolation between pairs of numbers in the Boundary (e.g. fBoundariesA) vector.  
  // The first pair of numbers represent the beginning and end of the Inner Field cage, etc.
  // The unit of the error potential vector is [Volt], whereas 1mm shift of the IFC would 
  // correspond to ~ 40 V
  // 
  // Note: The setting for the CE will be passed to the C side!
  
  for (Int_t i=0; i<8; i++) {
    fBoundariesA[i]= boundariesA[i];  
    if (i>5) fBoundariesC[i]= -boundariesA[i]; // setting for the CE is passed to C side
  }

}
void AliTPCBoundaryVoltError::SetBoundariesC(Float_t boundariesC[6]){
  //
  // set voltage errors on the TPC boundaries - A side 
  //
  // Start at IFC at the Central electrode and work clockwise (for C side) through 
  // IFC, ROC and OFC. The boundary conditions are currently defined to be a linear 
  // interpolation between pairs of numbers in the Boundary (e.g. fBoundariesC) vector.  
  // The first pair of numbers represent the beginning and end of the Inner Field cage, etc.
  // The unit of the error potential vector is [Volt], whereas 1mm shift of the IFC would 
  // correspond to ~ 40 V
  // 
  // Note: The setting for the CE will be taken from the A side (pos 6 and 7)!

  for (Int_t i=0; i<6; i++) {
    fBoundariesC[i]= boundariesC[i];  
  }

}
