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
//   <h2>AliTPCFCVoltError3D class   </h2>       
//   The class calculates the space point distortions due to residual voltage errors 
//   on the Field Cage (FC) boundaries (rods and strips) of the TPC in 3D. It uses 
//   the Poisson relaxation technique in three dimension as implemented in the parent class. 
//   <p>
//   Although the calculation is performed in 3D, the calculation time was significantly 
//   reduced by using certain symmetry conditions within the calculation.
//   <p>
//   The input parameters can be set via the functions (e.g.) SetRodVoltShift(rod,dV) where 
//   rod is the number of the rod on which the voltage offset dV is set. The corresponding 
//   shift in z direction would be $dz=dV/40$ with an opposite sign for the C side. The 
//   rods are numbered in anti-clock-wise direction starting at $\phi=0$. Rods in the IFC 
//   are from 0 to 17, rods on the OFC go from 18 to 35. <br>
//   This convention is following the offline numbering scheme of the ROCs.
//   <p>
//   Different misalignment scenarios can be modeled: 
//   <ul>
//   <li> A rotated clip scenario is only possible at two positions (Rod 11 at IFC, rod 3(+18) at OFC) 
//        and can be set via SetRotatedClipVolt. The key feature is that at the mentioned positions,
//        the strip-ends were combined. At these positions, a systematic offset of one strip-end in
//        respect to the other is possible. 
//   <li> A normal rod offset, where the strips follow the movement of the rod, can be set via 
//        SetRodVoltShift. It has a anti-mirrored signature in terms of distortions when compared 
//        to the rotated clip. This misalignment is possible at each single rod of the FC.
//   <li> A simple rod offset, where the strips do not follow the shift, results in an even more 
//        localized distortion close to the rod. The difference to a rod shift, where the strips follow,
//        is more dominant on the OFC. This effect can be set via the function SetCopperRodShift.
//   </ul>
// End_Html
//
// Begin_Macro(source)
//   {
//   gROOT->SetStyle("Plain"); gStyle->SetPalette(1);
//   TCanvas *c2 = new TCanvas("cAliTPCVoltError3D","cAliTPCVoltError3D",500,450); 
//   AliTPCFCVoltError3D fc;
//   fc.SetOmegaTauT1T2(0,1,1); 
//   fc.SetRotatedClipVoltA(0,40);
//   fc.SetRodVoltShiftA(3,40); 
//   fc.SetCopperRodShiftA(7+18,40);
//   fc.SetRodVoltShiftA(15+18,40); 
//   fc.CreateHistoDRPhiinXY(10)->Draw("cont4z");
//   return c2;
//   } 
// End_Macro
//
// Begin_Html
//   <p>
//   Date: 08/08/2010  <br>
//   Authors: Jim Thomas, Stefan Rossegger  
// End_Html 
// _________________________________________________________________


#include "AliMagF.h"
#include "TGeoGlobalMagField.h"
#include "AliTPCcalibDB.h"
#include "AliTPCParam.h"
#include "AliLog.h"
#include "TMatrixD.h"
#include "TMatrixF.h"

#include "TMath.h"
#include "AliTPCROC.h"
#include "AliTPCFCVoltError3D.h"

ClassImp(AliTPCFCVoltError3D)

AliTPCFCVoltError3D::AliTPCFCVoltError3D()
  : AliTPCCorrection("FieldCageVoltErrors","FieldCage (Rods) Voltage Errors"),
    fC0(0.),fC1(0.),
    fInitLookUp(kFALSE)
{
  //
  // default constructor
  //

  // flags for filled 'basic lookup tables'
  for (Int_t i=0; i<6; i++){
    fInitLookUpBasic[i]= kFALSE;  
  }

  // Boundary settings 
  for (Int_t i=0; i<36; i++){
    fRodVoltShiftA[i] = 0;  
    fRodVoltShiftC[i] = 0;  
  }
  for (Int_t i=0; i<2; i++){
    fRotatedClipVoltA[i] = 0;  
    fRotatedClipVoltC[i] = 0;  
  }
  // 
  for (Int_t i=0; i<36; i++){
    fCopperRodShiftA[i] = 0;  
    fCopperRodShiftC[i] = 0;  
  }

  // Array which will contain the solution according to the setted boundary conditions
  // it represents a sum up of the 4 basic look up tables (created individually)
  // see InitFCVoltError3D() function
  for ( Int_t k = 0 ; k < kNPhi ; k++ ) {
    fLookUpErOverEz[k]   =  new TMatrixF(kNR,kNZ);  
    fLookUpEphiOverEz[k] =  new TMatrixF(kNR,kNZ);
    fLookUpDeltaEz[k]    =  new TMatrixF(kNR,kNZ);   
  }
  
  for ( Int_t k = 0 ; k < kPhiSlices ; k++ ) {
    fLookUpBasic1ErOverEz[k]   = 0;
    fLookUpBasic1EphiOverEz[k] = 0; 
    fLookUpBasic1DeltaEz[k]    = 0;

    fLookUpBasic2ErOverEz[k]   = 0;
    fLookUpBasic2EphiOverEz[k] = 0; 
    fLookUpBasic2DeltaEz[k]    = 0;

    fLookUpBasic3ErOverEz[k]   = 0;
    fLookUpBasic3EphiOverEz[k] = 0; 
    fLookUpBasic3DeltaEz[k]    = 0;

    fLookUpBasic4ErOverEz[k]   = 0;
    fLookUpBasic4EphiOverEz[k] = 0; 
    fLookUpBasic4DeltaEz[k]    = 0;
    
    fLookUpBasic5ErOverEz[k]   = 0;
    fLookUpBasic5EphiOverEz[k] = 0; 
    fLookUpBasic5DeltaEz[k]    = 0;

    fLookUpBasic6ErOverEz[k]   = 0;
    fLookUpBasic6EphiOverEz[k] = 0; 
    fLookUpBasic6DeltaEz[k]    = 0;
  }

}

AliTPCFCVoltError3D::~AliTPCFCVoltError3D() {
  //
  // destructor
  //
  
  for ( Int_t k = 0 ; k < kNPhi ; k++ ) {
    delete fLookUpErOverEz[k];
    delete fLookUpEphiOverEz[k];
    delete fLookUpDeltaEz[k];
  }

  for ( Int_t k = 0 ; k < kPhiSlices ; k++ ) {
    delete fLookUpBasic1ErOverEz[k];  // does nothing if pointer is zero!
    delete fLookUpBasic1EphiOverEz[k]; 
    delete fLookUpBasic1DeltaEz[k]; 

    delete fLookUpBasic2ErOverEz[k];  // does nothing if pointer is zero!
    delete fLookUpBasic2EphiOverEz[k]; 
    delete fLookUpBasic2DeltaEz[k]; 
    
    delete fLookUpBasic3ErOverEz[k];  // does nothing if pointer is zero!
    delete fLookUpBasic3EphiOverEz[k]; 
    delete fLookUpBasic3DeltaEz[k]; 

    delete fLookUpBasic4ErOverEz[k];  // does nothing if pointer is zero!
    delete fLookUpBasic4EphiOverEz[k]; 
    delete fLookUpBasic4DeltaEz[k]; 

    delete fLookUpBasic5ErOverEz[k];  // does nothing if pointer is zero!
    delete fLookUpBasic5EphiOverEz[k]; 
    delete fLookUpBasic5DeltaEz[k]; 

    delete fLookUpBasic6ErOverEz[k];  // does nothing if pointer is zero!
    delete fLookUpBasic6EphiOverEz[k]; 
    delete fLookUpBasic6DeltaEz[k]; 

  }
}

void AliTPCFCVoltError3D::Init() {
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

  if (!fInitLookUp) InitFCVoltError3D();
}

void AliTPCFCVoltError3D::Update(const TTimeStamp &/*timeStamp*/) {
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



void AliTPCFCVoltError3D::GetCorrection(const Float_t x[],const Short_t roc,Float_t dx[]) {
  //
  // Calculates the correction due e.g. residual voltage errors on the TPC boundaries
  //   
  const Double_t kEpsilon=Double_t(FLT_MIN);

  if (!fInitLookUp) {
    AliInfo("Lookup table was not initialized! Perform the inizialisation now ...");
    InitFCVoltError3D();
  }

  static Bool_t forceInit=kTRUE; //temporary needed for back compatibility with old OCDB entries
  if (forceInit &&fLookUpErOverEz[0]){
    if (TMath::Abs(fLookUpErOverEz[0]->Sum())<kEpsilon){//temporary needed for back compatibility with old OCDB entries
      ForceInitFCVoltError3D();
    }
    forceInit=kFALSE;
  }


  Int_t   order     = 1 ;               // FIXME: hardcoded? Linear interpolation = 1, Quadratic = 2         
                                        // note that the poisson solution was linearly mirroed on this grid!
  Float_t intEr, intEphi, intDeltaEz;
  Float_t r, phi, z ;
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

  // Get the Er and Ephi field integrals plus the integral over DeltaEz 
  intEr      = Interpolate3DTable(order, r, z, phi, kNR, kNZ, kNPhi, 
				  fgkRList, fgkZList, fgkPhiList, fLookUpErOverEz  );
  intEphi    = Interpolate3DTable(order, r, z, phi, kNR, kNZ, kNPhi, 
				  fgkRList, fgkZList, fgkPhiList, fLookUpEphiOverEz  );
  intDeltaEz = Interpolate3DTable(order, r, z, phi, kNR, kNZ, kNPhi, 
				  fgkRList, fgkZList, fgkPhiList, fLookUpDeltaEz  );

  //  printf("%lf %lf %lf\n",intEr,intEphi,intDeltaEz);

  // Calculate distorted position
  if ( r > 0.0 ) {
    phi =  phi + ( fC0*intEphi - fC1*intEr ) / r;      
    r   =  r   + ( fC0*intEr   + fC1*intEphi );  
  }
  
  // Calculate correction in cartesian coordinates
  dx[0] = r * TMath::Cos(phi) - x[0];
  dx[1] = r * TMath::Sin(phi) - x[1]; 
  dx[2] = intDeltaEz;  // z distortion - (internally scaled with driftvelocity dependency 
                       // on the Ez field plus the actual ROC misalignment (if set TRUE)

}

void AliTPCFCVoltError3D::InitFCVoltError3D() {
  //
  // Initialization of the Lookup table which contains the solutions of the 
  // Dirichlet boundary problem
  // Calculation of the single 3D-Poisson solver is done just if needed
  // (see basic lookup tables in header file)
  //

  const Int_t   order       =    1  ;  // Linear interpolation = 1, Quadratic = 2  
  const Float_t gridSizeR   =  (fgkOFCRadius-fgkIFCRadius) / (kRows-1) ;
  const Float_t gridSizeZ   =  fgkTPCZ0 / (kColumns-1) ;
  const Float_t gridSizePhi =  TMath::TwoPi() / ( 18.0 * kPhiSlicesPerSector);

  // temporary arrays to create the boundary conditions
  TMatrixD *arrayofArrayV[kPhiSlices], *arrayofCharge[kPhiSlices] ; 
  for ( Int_t k = 0 ; k < kPhiSlices ; k++ ) {
    arrayofArrayV[k]     =   new TMatrixD(kRows,kColumns) ;
    arrayofCharge[k]     =   new TMatrixD(kRows,kColumns) ;
  }
  Double_t  innerList[kPhiSlices], outerList[kPhiSlices] ;
  
  // list of point as used in the poisson relation and the interpolation (during sum up)
  Double_t  rlist[kRows], zedlist[kColumns] , philist[kPhiSlices];
  for ( Int_t k = 0 ; k < kPhiSlices ; k++ ) {
    philist[k] =  gridSizePhi * k;
    for ( Int_t i = 0 ; i < kRows ; i++ )    {
      rlist[i] = fgkIFCRadius + i*gridSizeR ;
      for ( Int_t j = 0 ; j < kColumns ; j++ ) { // Fill Vmatrix with Boundary Conditions
	zedlist[j]  = j * gridSizeZ ;
      }
    }
  }

  // ==========================================================================
  // Solve Poisson's equation in 3D cylindrical coordinates by relaxation technique
  // Allow for different size grid spacing in R and Z directions
  
  const Int_t   symmetry[6] =    {1,1,-1,-1,1,1}; // shifted rod (2x), rotated clip (2x), only rod shift on OFC (1x)

  // check which lookup tables are needed ---------------------------------

  Bool_t needTable[6] ={kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE};

  // Shifted rods & strips
  for ( Int_t rod = 0 ; rod < 18 ; rod++ ) {
    if (fRodVoltShiftA[rod]!=0) needTable[0]=kTRUE;
    if (fRodVoltShiftC[rod]!=0) needTable[0]=kTRUE;
  }
  for ( Int_t rod = 18 ; rod < 36 ; rod++ ) {
    if (fRodVoltShiftA[rod]!=0) needTable[1]=kTRUE;
    if (fRodVoltShiftC[rod]!=0) needTable[1]=kTRUE;
  }
  // Rotated clips
  if (fRotatedClipVoltA[0]!=0 || fRotatedClipVoltC[0]!=0) needTable[2]=kTRUE;
  if (fRotatedClipVoltA[1]!=0 || fRotatedClipVoltC[1]!=0) needTable[3]=kTRUE;
 
  // shifted Copper rods 
  for ( Int_t rod = 0 ; rod < 18 ; rod++ ) {
    if (fCopperRodShiftA[rod]!=0) needTable[4]=kTRUE;
    if (fCopperRodShiftC[rod]!=0) needTable[4]=kTRUE;
  }
  // shifted Copper rods 
  for ( Int_t rod = 18; rod < 36 ; rod++ ) {
    if (fCopperRodShiftA[rod]!=0) needTable[5]=kTRUE;
    if (fCopperRodShiftC[rod]!=0) needTable[5]=kTRUE;
  }


  // reserve the arrays for the basic lookup tables ----------------------
  if (needTable[0] && !fInitLookUpBasic[0]) {
    for ( Int_t k = 0 ; k < kPhiSlices ; k++ ) {   // Possibly make an extra table to be used for 0 == 360
      fLookUpBasic1ErOverEz[k]   =   new TMatrixD(kRows,kColumns);
      fLookUpBasic1EphiOverEz[k] =   new TMatrixD(kRows,kColumns);
      fLookUpBasic1DeltaEz[k]    =   new TMatrixD(kRows,kColumns);
      // will be deleted in destructor
    }
  }
  if (needTable[1] && !fInitLookUpBasic[1]) {
    for ( Int_t k = 0 ; k < kPhiSlices ; k++ ) {   // Possibly make an extra table to be used for 0 == 360
      fLookUpBasic2ErOverEz[k]   =   new TMatrixD(kRows,kColumns);
      fLookUpBasic2EphiOverEz[k] =   new TMatrixD(kRows,kColumns);
      fLookUpBasic2DeltaEz[k]    =   new TMatrixD(kRows,kColumns);
      // will be deleted in destructor
    }
  }
  if (needTable[2] && !fInitLookUpBasic[2]) {
    for ( Int_t k = 0 ; k < kPhiSlices ; k++ ) {   // Possibly make an extra table to be used for 0 == 360
      fLookUpBasic3ErOverEz[k]   =   new TMatrixD(kRows,kColumns);
      fLookUpBasic3EphiOverEz[k] =   new TMatrixD(kRows,kColumns);
      fLookUpBasic3DeltaEz[k]    =   new TMatrixD(kRows,kColumns);
      // will be deleted in destructor
    }
  }
  if (needTable[3] && !fInitLookUpBasic[3]) {
    for ( Int_t k = 0 ; k < kPhiSlices ; k++ ) {   // Possibly make an extra table to be used for 0 == 360
      fLookUpBasic4ErOverEz[k]   =   new TMatrixD(kRows,kColumns);
      fLookUpBasic4EphiOverEz[k] =   new TMatrixD(kRows,kColumns);
      fLookUpBasic4DeltaEz[k]    =   new TMatrixD(kRows,kColumns);
      // will be deleted in destructor
    }
  }
  if (needTable[4] && !fInitLookUpBasic[4]) {
    for ( Int_t k = 0 ; k < kPhiSlices ; k++ ) {   // Possibly make an extra table to be used for 0 == 360
      fLookUpBasic5ErOverEz[k]   =   new TMatrixD(kRows,kColumns);
      fLookUpBasic5EphiOverEz[k] =   new TMatrixD(kRows,kColumns);
      fLookUpBasic5DeltaEz[k]    =   new TMatrixD(kRows,kColumns);
      // will be deleted in destructor
    }
  }
  if (needTable[5] && !fInitLookUpBasic[5]) {
    for ( Int_t k = 0 ; k < kPhiSlices ; k++ ) {   // Possibly make an extra table to be used for 0 == 360
      fLookUpBasic6ErOverEz[k]   =   new TMatrixD(kRows,kColumns);
      fLookUpBasic6EphiOverEz[k] =   new TMatrixD(kRows,kColumns);
      fLookUpBasic6DeltaEz[k]    =   new TMatrixD(kRows,kColumns);
      // will be deleted in destructor
    }
  }
 
  // Set bondaries and solve Poisson's equation --------------------------
 
  for (Int_t look=0; look<6; look++) {
   
    Float_t inner18[18] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 } ;  
    Float_t outer18[18] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 } ; 
  
    if (needTable[look] && !fInitLookUpBasic[look]) {

      // Specify which rod is the reference; other distortions calculated by summing over multiple rotations of refrence
      // Note: the interpolation later on depends on it!! Do not change if not really needed!
      if (look==0) {
	AliInfo(Form("IFC - ROD&Strip shift : Filling table (~ %d sec)",3*(int)(kPhiSlices/5)));
	inner18[0] = 1;  
      } else if (look==1) {
	AliInfo(Form("OFC - ROD&Strip shift : Filling table (~ %d sec)",3*(int)(kPhiSlices/5)));
	outer18[0] = 1;  
      } else if (look==2) {
	AliInfo(Form("IFC - Clip rot. : Filling table (~ %d sec)",3*(int)(kPhiSlices/5)));
	inner18[0] = 1; 
      } else if (look==3) {
	AliInfo(Form("OFC - Clip rot. : Filling table (~ %d sec)",3*(int)(kPhiSlices/5)));
	outer18[0] = 1;  
      } else if (look==4) {
	AliInfo(Form("IFC - CopperRod shift : Filling table (~ %d sec)",3*(int)(kPhiSlices/5)));
	inner18[0] = 1;  
      } else if (look==5) {
	AliInfo(Form("OFC - CopperRod shift : Filling table (~ %d sec)",3*(int)(kPhiSlices/5)));
	outer18[0] = 1;  
      }
      // Build potentials on boundary stripes for a rotated clip (SYMMETRY==-1) or a shifted rod (SYMMETRY==1)
      if (look<4) {
	// in these cases, the strips follow the actual rod shift (linear interpolation between the rods)
	for ( Int_t k = 0 ; k < kPhiSlices ; k++ ) {
	  Int_t slices = kPhiSlicesPerSector ;
	  Int_t ipoint = k/slices  ;  
	  innerList[k] = (((ipoint+1)*slices-k)*inner18[ipoint]-(k-ipoint*slices)*inner18[ipoint+1])/slices ;
	  outerList[k] = (((ipoint+1)*slices-k)*outer18[ipoint]-(k-ipoint*slices)*outer18[ipoint+1])/slices ;
	  if ( (k%slices) == 0 && symmetry[look] == -1 ) { innerList[k] = 0.0 ; outerList[k] = 0.0 ; } 
	  // above, force Zero if Anti-Sym
	} 
      } else if (look==4){
	// in this case, we assume that the strips stay at the correct position, but the rods move
	// the distortion is then much more localized around the rod (indicated by real data)
	for ( Int_t k = 0 ; k < kPhiSlices ; k++ )
	  innerList[k] = outerList[k] = 0;
	innerList[0] = inner18[0];     // point at rod
	innerList[0] = inner18[0]/4*3;  // point close to rod (educated guess)
      } else if (look==5){
	// in this case, we assume that the strips stay at the correct position, but the copper plated OFC-rods move
	// the distortion is then much more localized around the rod (indicated by real data)
	for ( Int_t k = 0 ; k < kPhiSlices ; k++ )
	  innerList[k] = outerList[k] = 0;
	outerList[0] = outer18[0];   // point at rod
	outerList[1] = outer18[0]/4; // point close to rod (educated-`guessed` scaling)
      }

      // Fill arrays with initial conditions.  V on the boundary and Charge in the volume.
      for ( Int_t k = 0 ; k < kPhiSlices ; k++ ) {
	TMatrixD &arrayV    =  *arrayofArrayV[k] ;
	TMatrixD &charge    =  *arrayofCharge[k] ;
       	for ( Int_t i = 0 ; i < kRows ; i++ )    {
	  for ( Int_t j = 0 ; j < kColumns ; j++ ) { // Fill Vmatrix with Boundary Conditions
	    arrayV(i,j) = 0.0 ; 
	    charge(i,j) = 0.0 ;
	    if ( i == 0 )       arrayV(i,j) = innerList[k] ; 
	    if ( i == kRows-1 ) arrayV(i,j) = outerList[k] ; 
	  }
	}      
	// no charge in the volume
	for ( Int_t i = 1 ; i < kRows-1 ; i++ )  { 
	  for ( Int_t j = 1 ; j < kColumns-1 ; j++ ) {
	    charge(i,j)  =  0.0 ;
	  }
	}
      }      
     
      // Solve Poisson's equation in 3D cylindrical coordinates by relaxation technique
      // Allow for different size grid spacing in R and Z directions
      if (look==0) {
	PoissonRelaxation3D( arrayofArrayV, arrayofCharge, 
			     fLookUpBasic1ErOverEz, fLookUpBasic1EphiOverEz, fLookUpBasic1DeltaEz,
			     kRows, kColumns, kPhiSlices, gridSizePhi, kIterations, symmetry[0]) ;
	AliInfo("IFC - ROD&Strip shift : done ");
      } else if (look==1) {
	PoissonRelaxation3D( arrayofArrayV, arrayofCharge, 
			     fLookUpBasic2ErOverEz, fLookUpBasic2EphiOverEz, fLookUpBasic2DeltaEz,
			     kRows, kColumns, kPhiSlices, gridSizePhi, kIterations, symmetry[1]) ;
	
	AliInfo("OFC - ROD&Strip shift : done ");
      } else if (look==2) {
	PoissonRelaxation3D( arrayofArrayV, arrayofCharge, 
			     fLookUpBasic3ErOverEz, fLookUpBasic3EphiOverEz, fLookUpBasic3DeltaEz,
			     kRows, kColumns, kPhiSlices, gridSizePhi, kIterations, symmetry[2]) ;
	AliInfo("IFC - Clip rot. : done ");
      } else if (look==3) {
	PoissonRelaxation3D( arrayofArrayV, arrayofCharge, 
			     fLookUpBasic4ErOverEz, fLookUpBasic4EphiOverEz, fLookUpBasic4DeltaEz,
			     kRows, kColumns, kPhiSlices, gridSizePhi, kIterations, symmetry[3]) ;
	AliInfo("OFC - Clip rot. : done ");
      } else if (look==4) {
	PoissonRelaxation3D( arrayofArrayV, arrayofCharge, 
			     fLookUpBasic5ErOverEz, fLookUpBasic5EphiOverEz, fLookUpBasic5DeltaEz,
			     kRows, kColumns, kPhiSlices, gridSizePhi, kIterations, symmetry[4]) ;
	AliInfo("IFC - CopperRod shift : done ");
      } else if (look==5) {
	PoissonRelaxation3D( arrayofArrayV, arrayofCharge, 
			     fLookUpBasic6ErOverEz, fLookUpBasic6EphiOverEz, fLookUpBasic6DeltaEz,
			     kRows, kColumns, kPhiSlices, gridSizePhi, kIterations, symmetry[5]) ;
	AliInfo("OFC - CopperRod shift : done ");
      }
      
      fInitLookUpBasic[look] = kTRUE;
    }
  }
  

  // =============================================================================
  // Create the final lookup table with corresponding ROD shifts or clip rotations

  Float_t erIntegralSum   = 0.0 ;
  Float_t ephiIntegralSum = 0.0 ;
  Float_t ezIntegralSum   = 0.0 ;

  Double_t phiPrime     = 0. ;
  Double_t erIntegral   = 0. ;
  Double_t ephiIntegral = 0. ;
  Double_t ezIntegral   = 0. ;


  AliInfo("Calculation of combined Look-up Table");

  // Interpolate and sum the Basic lookup tables onto the standard grid
  Double_t  r, phi, z ;
  for ( Int_t k = 0 ; k < kNPhi ; k++ ) {
    phi = fgkPhiList[k] ;

    TMatrixF &erOverEz   =  *fLookUpErOverEz[k]  ;
    TMatrixF &ephiOverEz =  *fLookUpEphiOverEz[k];
    TMatrixF &deltaEz    =  *fLookUpDeltaEz[k]   ;

    for ( Int_t i = 0 ; i < kNZ ; i++ ) {
      z = TMath::Abs(fgkZList[i]) ;  // Symmetric solution in Z that depends only on ABS(Z)
      for ( Int_t j = 0 ; j < kNR ; j++ ) { 
	r = fgkRList[j] ;
	// Interpolate basicLookup tables; once for each rod, then sum the results
	
	erIntegralSum   = 0.0 ;
	ephiIntegralSum = 0.0 ;
	ezIntegralSum   = 0.0 ;
 
	// SHIFTED RODS =========================================================

	// IFC ROD SHIFTS +++++++++++++++++++++++++++++
	for ( Int_t rod = 0 ; rod < 18 ; rod++ ) {
	  
	  erIntegral = ephiIntegral = ezIntegral = 0.0 ;
	  
	  if ( fRodVoltShiftA[rod] == 0 && fgkZList[i] > 0) continue ;
	  if ( fRodVoltShiftC[rod] == 0 && fgkZList[i] < 0) continue ;

	  // Rotate to a position where we have maps and prepare to sum
	  phiPrime =  phi - rod*kPhiSlicesPerSector*gridSizePhi ;  

	  if ( phiPrime < 0 ) phiPrime += TMath::TwoPi() ;   // Stay in range from 0 to TwoPi    

	  if ( (phiPrime >= 0) && (phiPrime <= kPhiSlices*gridSizePhi) ) {
	    
	    erIntegral   = Interpolate3DTable(order, r, z, phiPrime, kRows, kColumns, kPhiSlices, 
					      rlist, zedlist, philist, fLookUpBasic1ErOverEz  );
	    ephiIntegral = Interpolate3DTable(order, r, z, phiPrime, kRows, kColumns, kPhiSlices,
					      rlist, zedlist, philist, fLookUpBasic1EphiOverEz);
	    ezIntegral   = Interpolate3DTable(order, r, z, phiPrime, kRows, kColumns, kPhiSlices,
					      rlist, zedlist, philist, fLookUpBasic1DeltaEz   );
	    
	  }  else if ( (phiPrime <= TMath::TwoPi()) && (phiPrime >= (TMath::TwoPi()-kPhiSlices*gridSizePhi)) ){
	    
	    phiPrime     = TMath::TwoPi() - phiPrime ;
	    
	    erIntegral   = Interpolate3DTable(order, r, z, phiPrime, kRows, kColumns, kPhiSlices, 
					      rlist, zedlist, philist, fLookUpBasic1ErOverEz  );
	    ephiIntegral = Interpolate3DTable(order, r, z, phiPrime, kRows, kColumns, kPhiSlices,
					      rlist, zedlist, philist, fLookUpBasic1EphiOverEz);
	    ezIntegral   = Interpolate3DTable(order, r, z, phiPrime, kRows, kColumns, kPhiSlices,
					      rlist, zedlist, philist, fLookUpBasic1DeltaEz   );
	    
	    // Flip symmetry of phi integral for symmetric boundary conditions
	    if ( symmetry[0] ==  1 ) ephiIntegral *= -1  ;    
	    // Flip symmetry of r integral if anti-symmetric boundary conditions 
	    if ( symmetry[0] == -1 ) erIntegral   *= -1  ;    

	  }

	  if ( fgkZList[i] > 0 ) {
	    erIntegralSum   += fRodVoltShiftA[rod]*erIntegral   ;
	    ephiIntegralSum += fRodVoltShiftA[rod]*ephiIntegral ;
	    ezIntegralSum   += fRodVoltShiftA[rod]*ezIntegral;
	  } else if ( fgkZList[i] < 0 ) {
	    erIntegralSum   += fRodVoltShiftC[rod]*erIntegral   ;
	    ephiIntegralSum += fRodVoltShiftC[rod]*ephiIntegral ;
	    ezIntegralSum   -= fRodVoltShiftC[rod]*ezIntegral;
	  }
	}

	// OFC ROD SHIFTS +++++++++++++++++++++++++++++
	for ( Int_t rod = 18 ; rod < 36 ; rod++ ) {

	  erIntegral = ephiIntegral = ezIntegral = 0.0 ;
	  
	  if ( fRodVoltShiftA[rod] == 0 && fgkZList[i] > 0) continue ;
	  if ( fRodVoltShiftC[rod] == 0 && fgkZList[i] < 0) continue ;

	  // Rotate to a position where we have maps and prepare to sum
	  phiPrime =  phi - (rod-18)*kPhiSlicesPerSector*gridSizePhi ;  
		  
	  if ( phiPrime < 0 ) phiPrime += TMath::TwoPi() ;   // Stay in range from 0 to TwoPi    

	  if ( (phiPrime >= 0) && (phiPrime <= kPhiSlices*gridSizePhi) ) {
	    
	    erIntegral   = Interpolate3DTable(order, r, z, phiPrime, kRows, kColumns, kPhiSlices, 
					      rlist, zedlist, philist, fLookUpBasic2ErOverEz  );
	    ephiIntegral = Interpolate3DTable(order, r, z, phiPrime, kRows, kColumns, kPhiSlices,
					      rlist, zedlist, philist, fLookUpBasic2EphiOverEz);
	    ezIntegral   = Interpolate3DTable(order, r, z, phiPrime, kRows, kColumns, kPhiSlices,
					      rlist, zedlist, philist, fLookUpBasic2DeltaEz   );
	    
	  }  else if ( (phiPrime <= TMath::TwoPi()) && (phiPrime >= (TMath::TwoPi()-kPhiSlices*gridSizePhi)) ){
	    
	    phiPrime     = TMath::TwoPi() - phiPrime ;
	    
	    erIntegral   = Interpolate3DTable(order, r, z, phiPrime, kRows, kColumns, kPhiSlices, 
					      rlist, zedlist, philist, fLookUpBasic2ErOverEz  );
	    ephiIntegral = Interpolate3DTable(order, r, z, phiPrime, kRows, kColumns, kPhiSlices,
					      rlist, zedlist, philist, fLookUpBasic2EphiOverEz);
	    ezIntegral   = Interpolate3DTable(order, r, z, phiPrime, kRows, kColumns, kPhiSlices,
					      rlist, zedlist, philist, fLookUpBasic2DeltaEz   );
	    
	    // Flip symmetry of phi integral for symmetric boundary conditions
	    if ( symmetry[1] ==  1 ) ephiIntegral *= -1  ;    
	    // Flip symmetry of r integral if anti-symmetric boundary conditions 
	    if ( symmetry[1] == -1 ) erIntegral   *= -1  ;    

	  }

	  if ( fgkZList[i] > 0 ) {
	    erIntegralSum   += fRodVoltShiftA[rod]*erIntegral   ;
	    ephiIntegralSum += fRodVoltShiftA[rod]*ephiIntegral ;
	    ezIntegralSum   += fRodVoltShiftA[rod]*ezIntegral;
	  } else if ( fgkZList[i] < 0 ) {
	    erIntegralSum   += fRodVoltShiftC[rod]*erIntegral   ;
	    ephiIntegralSum += fRodVoltShiftC[rod]*ephiIntegral ;
	    ezIntegralSum   -= fRodVoltShiftC[rod]*ezIntegral;
	  }

	} // rod loop - shited ROD


	// Rotated clips =========================================================

	Int_t rodIFC = 11; // resistor rod positions, IFC
	Int_t rodOFC =  3; // resistor rod positions, OFC
	// just one rod on IFC and OFC

	// IFC ROTATED CLIP +++++++++++++++++++++++++++++
	for ( Int_t rod = rodIFC ; rod < rodIFC+1 ; rod++ ) { // loop over 1 to keep the "ignore"-functionality

	  erIntegral = ephiIntegral = ezIntegral = 0.0 ;
	  if ( fRotatedClipVoltA[0] == 0 && fgkZList[i] > 0) continue ;
	  if ( fRotatedClipVoltC[0] == 0 && fgkZList[i] < 0) continue ;

	  // Rotate to a position where we have maps and prepare to sum
	  phiPrime =  phi - rod*kPhiSlicesPerSector*gridSizePhi ;  
	  
	  if ( phiPrime < 0 ) phiPrime += TMath::TwoPi() ;   // Stay in range from 0 to TwoPi    
	  
	  if ( (phiPrime >= 0) && (phiPrime <= kPhiSlices*gridSizePhi) ) {
	    
	    erIntegral   = Interpolate3DTable(order, r, z, phiPrime, kRows, kColumns, kPhiSlices, 
					      rlist, zedlist, philist, fLookUpBasic3ErOverEz  );
	    ephiIntegral = Interpolate3DTable(order, r, z, phiPrime, kRows, kColumns, kPhiSlices,
					      rlist, zedlist, philist, fLookUpBasic3EphiOverEz);
	    ezIntegral   = Interpolate3DTable(order, r, z, phiPrime, kRows, kColumns, kPhiSlices,
					      rlist, zedlist, philist, fLookUpBasic3DeltaEz   );
	    
	  }  else if ( (phiPrime <= TMath::TwoPi()) && (phiPrime >= (TMath::TwoPi()-kPhiSlices*gridSizePhi)) ){
	    
	    phiPrime     = TMath::TwoPi() - phiPrime ;
	    
	    erIntegral   = Interpolate3DTable(order, r, z, phiPrime, kRows, kColumns, kPhiSlices, 
					      rlist, zedlist, philist, fLookUpBasic3ErOverEz  );
	    ephiIntegral = Interpolate3DTable(order, r, z, phiPrime, kRows, kColumns, kPhiSlices,
					      rlist, zedlist, philist, fLookUpBasic3EphiOverEz);
	    ezIntegral   = Interpolate3DTable(order, r, z, phiPrime, kRows, kColumns, kPhiSlices,
					      rlist, zedlist, philist, fLookUpBasic3DeltaEz   );
	    
	    // Flip symmetry of phi integral for symmetric boundary conditions
	    if ( symmetry[2] ==  1 ) ephiIntegral *= -1  ;    
	    // Flip symmetry of r integral if anti-symmetric boundary conditions 
	    if ( symmetry[2] == -1 ) erIntegral   *= -1  ;    
	    
	  }
	  
	  if ( fgkZList[i] > 0 ) {
	    erIntegralSum   += fRotatedClipVoltA[0]*erIntegral   ;
	    ephiIntegralSum += fRotatedClipVoltA[0]*ephiIntegral ;
	    ezIntegralSum   += fRotatedClipVoltA[0]*ezIntegral;
	  } else if ( fgkZList[i] < 0 ) {
	    erIntegralSum   += fRotatedClipVoltC[0]*erIntegral   ;
	    ephiIntegralSum += fRotatedClipVoltC[0]*ephiIntegral ;
	    ezIntegralSum   -= fRotatedClipVoltC[0]*ezIntegral;
	  }
	}

	// OFC: ROTATED CLIP  +++++++++++++++++++++++++++++
	for ( Int_t rod = rodOFC ; rod < rodOFC+1 ; rod++ ) { // loop over 1 to keep the "ignore"-functionality
	  
	  erIntegral = ephiIntegral = ezIntegral = 0.0 ;
	  
	  if ( fRotatedClipVoltA[1] == 0 && fgkZList[i] > 0) continue ;
	  if ( fRotatedClipVoltC[1] == 0 && fgkZList[i] < 0) continue ;

	  // Rotate to a position where we have maps and prepare to sum
	  phiPrime =  phi - rod*kPhiSlicesPerSector*gridSizePhi ;  
	  
	  
	  if ( phiPrime < 0 ) phiPrime += TMath::TwoPi() ;   // Stay in range from 0 to TwoPi    
	  
	  if ( (phiPrime >= 0) && (phiPrime <= kPhiSlices*gridSizePhi) ) {
	    
	    erIntegral   = Interpolate3DTable(order, r, z, phiPrime, kRows, kColumns, kPhiSlices, 
					      rlist, zedlist, philist, fLookUpBasic4ErOverEz  );
	    ephiIntegral = Interpolate3DTable(order, r, z, phiPrime, kRows, kColumns, kPhiSlices,
					      rlist, zedlist, philist, fLookUpBasic4EphiOverEz);
	    ezIntegral   = Interpolate3DTable(order, r, z, phiPrime, kRows, kColumns, kPhiSlices,
					      rlist, zedlist, philist, fLookUpBasic4DeltaEz   );
	    
	  }  else if ( (phiPrime <= TMath::TwoPi()) && (phiPrime >= (TMath::TwoPi()-kPhiSlices*gridSizePhi)) ){
	    
	    phiPrime     = TMath::TwoPi() - phiPrime ;
	    
	    erIntegral   = Interpolate3DTable(order, r, z, phiPrime, kRows, kColumns, kPhiSlices, 
					      rlist, zedlist, philist, fLookUpBasic4ErOverEz  );
	    ephiIntegral = Interpolate3DTable(order, r, z, phiPrime, kRows, kColumns, kPhiSlices,
					      rlist, zedlist, philist, fLookUpBasic4EphiOverEz);
	    ezIntegral   = Interpolate3DTable(order, r, z, phiPrime, kRows, kColumns, kPhiSlices,
					      rlist, zedlist, philist, fLookUpBasic4DeltaEz   );
	    
	    // Flip symmetry of phi integral for symmetric boundary conditions
	    if ( symmetry[3] ==  1 ) ephiIntegral *= -1  ;    
	    // Flip symmetry of r integral if anti-symmetric boundary conditions 
	    if ( symmetry[3] == -1 ) erIntegral   *= -1  ;    
	    
	  }
	  
	  if ( fgkZList[i] > 0 ) {
	    erIntegralSum   += fRotatedClipVoltA[1]*erIntegral   ;
	    ephiIntegralSum += fRotatedClipVoltA[1]*ephiIntegral ;
	    ezIntegralSum   += fRotatedClipVoltA[1]*ezIntegral;
	  } else if ( fgkZList[i] < 0 ) {
	    erIntegralSum   += fRotatedClipVoltC[1]*erIntegral   ;
	    ephiIntegralSum += fRotatedClipVoltC[1]*ephiIntegral ;
	    ezIntegralSum   -= fRotatedClipVoltC[1]*ezIntegral;
	  }
	}

	// IFC Cooper rod shift  +++++++++++++++++++++++++++++
	for ( Int_t rod = 0 ; rod < 18 ; rod++ ) {
	  
	  erIntegral = ephiIntegral = ezIntegral = 0.0 ;
	  
	  if ( fCopperRodShiftA[rod] == 0 && fgkZList[i] > 0) continue ;
	  if ( fCopperRodShiftC[rod] == 0 && fgkZList[i] < 0) continue ;

	  // Rotate to a position where we have maps and prepare to sum
	  phiPrime =  phi - rod*kPhiSlicesPerSector*gridSizePhi ;  

	  if ( phiPrime < 0 ) phiPrime += TMath::TwoPi() ;   // Stay in range from 0 to TwoPi    

	  if ( (phiPrime >= 0) && (phiPrime <= kPhiSlices*gridSizePhi) ) {
	    
	    erIntegral   = Interpolate3DTable(order, r, z, phiPrime, kRows, kColumns, kPhiSlices, 
					      rlist, zedlist, philist, fLookUpBasic5ErOverEz  );
	    ephiIntegral = Interpolate3DTable(order, r, z, phiPrime, kRows, kColumns, kPhiSlices,
					      rlist, zedlist, philist, fLookUpBasic5EphiOverEz);
	    ezIntegral   = Interpolate3DTable(order, r, z, phiPrime, kRows, kColumns, kPhiSlices,
					      rlist, zedlist, philist, fLookUpBasic5DeltaEz   );
	    
	  }  else if ( (phiPrime <= TMath::TwoPi()) && (phiPrime >= (TMath::TwoPi()-kPhiSlices*gridSizePhi)) ){
	    
	    phiPrime     = TMath::TwoPi() - phiPrime ;
	    
	    erIntegral   = Interpolate3DTable(order, r, z, phiPrime, kRows, kColumns, kPhiSlices, 
					      rlist, zedlist, philist, fLookUpBasic5ErOverEz  );
	    ephiIntegral = Interpolate3DTable(order, r, z, phiPrime, kRows, kColumns, kPhiSlices,
					      rlist, zedlist, philist, fLookUpBasic5EphiOverEz);
	    ezIntegral   = Interpolate3DTable(order, r, z, phiPrime, kRows, kColumns, kPhiSlices,
					      rlist, zedlist, philist, fLookUpBasic5DeltaEz   );
	    
	    // Flip symmetry of phi integral for symmetric boundary conditions
	    if ( symmetry[4] ==  1 ) ephiIntegral *= -1  ;    
	    // Flip symmetry of r integral if anti-symmetric boundary conditions 
	    if ( symmetry[4] == -1 ) erIntegral   *= -1  ;    

	  }

	  if ( fgkZList[i] > 0 ) {
	    erIntegralSum   += fCopperRodShiftA[rod]*erIntegral   ;
	    ephiIntegralSum += fCopperRodShiftA[rod]*ephiIntegral ;
	    ezIntegralSum   += fCopperRodShiftA[rod]*ezIntegral;
	  } else if ( fgkZList[i] < 0 ) {
	    erIntegralSum   += fCopperRodShiftC[rod]*erIntegral   ;
	    ephiIntegralSum += fCopperRodShiftC[rod]*ephiIntegral ;
	    ezIntegralSum   -= fCopperRodShiftC[rod]*ezIntegral;
	  }
	}

	// OFC Cooper rod shift  +++++++++++++++++++++++++++++
	for ( Int_t rod = 18 ; rod < 36 ; rod++ ) {
	  
	  erIntegral = ephiIntegral = ezIntegral = 0.0 ;
	  
	  if ( fCopperRodShiftA[rod] == 0 && fgkZList[i] > 0) continue ;
	  if ( fCopperRodShiftC[rod] == 0 && fgkZList[i] < 0) continue ;

	  // Rotate to a position where we have maps and prepare to sum
	  phiPrime =  phi - (rod-18)*kPhiSlicesPerSector*gridSizePhi ;  

	  if ( phiPrime < 0 ) phiPrime += TMath::TwoPi() ;   // Stay in range from 0 to TwoPi    

	  if ( (phiPrime >= 0) && (phiPrime <= kPhiSlices*gridSizePhi) ) {
	    
	    erIntegral   = Interpolate3DTable(order, r, z, phiPrime, kRows, kColumns, kPhiSlices, 
					      rlist, zedlist, philist, fLookUpBasic6ErOverEz  );
	    ephiIntegral = Interpolate3DTable(order, r, z, phiPrime, kRows, kColumns, kPhiSlices,
					      rlist, zedlist, philist, fLookUpBasic6EphiOverEz);
	    ezIntegral   = Interpolate3DTable(order, r, z, phiPrime, kRows, kColumns, kPhiSlices,
					      rlist, zedlist, philist, fLookUpBasic6DeltaEz   );
	    
	  }  else if ( (phiPrime <= TMath::TwoPi()) && (phiPrime >= (TMath::TwoPi()-kPhiSlices*gridSizePhi)) ){
	    
	    phiPrime     = TMath::TwoPi() - phiPrime ;
	    
	    erIntegral   = Interpolate3DTable(order, r, z, phiPrime, kRows, kColumns, kPhiSlices, 
					      rlist, zedlist, philist, fLookUpBasic6ErOverEz  );
	    ephiIntegral = Interpolate3DTable(order, r, z, phiPrime, kRows, kColumns, kPhiSlices,
					      rlist, zedlist, philist, fLookUpBasic6EphiOverEz);
	    ezIntegral   = Interpolate3DTable(order, r, z, phiPrime, kRows, kColumns, kPhiSlices,
					      rlist, zedlist, philist, fLookUpBasic6DeltaEz   );
	    
	    // Flip symmetry of phi integral for symmetric boundary conditions
	    if ( symmetry[5] ==  1 ) ephiIntegral *= -1  ;    
	    // Flip symmetry of r integral if anti-symmetric boundary conditions 
	    if ( symmetry[5] == -1 ) erIntegral   *= -1  ;    

	  }

	  if ( fgkZList[i] > 0 ) {
	    erIntegralSum   += fCopperRodShiftA[rod]*erIntegral   ;
	    ephiIntegralSum += fCopperRodShiftA[rod]*ephiIntegral ;
	    ezIntegralSum   += fCopperRodShiftA[rod]*ezIntegral;
	  } else if ( fgkZList[i] < 0 ) {
	    erIntegralSum   += fCopperRodShiftC[rod]*erIntegral   ;
	    ephiIntegralSum += fCopperRodShiftC[rod]*ephiIntegral ;
	    ezIntegralSum   -= fCopperRodShiftC[rod]*ezIntegral;
	  }
	}

	// put the sum into the final lookup table
	erOverEz(j,i)   =  erIntegralSum;
	ephiOverEz(j,i) =  ephiIntegralSum;
	deltaEz(j,i)    =  ezIntegralSum;
	
	//	if (j==1) printf("%lf %lf %lf: %lf %lf %lf\n",r, phi, z, erIntegralSum,ephiIntegralSum,ezIntegralSum);
 
      } // endo r loop
    } // end of z loop
  } // end of phi loop


  // clear the temporary arrays lists
  for ( Int_t k = 0 ; k < kPhiSlices ; k++ )  {
    delete arrayofArrayV[k];
    delete arrayofCharge[k];
  }
 
  AliInfo(" done");
  fInitLookUp = kTRUE;

}

void AliTPCFCVoltError3D::Print(const Option_t* option) const {
  //
  // Print function to check the settings of the Rod shifts and the rotated clips
  // option=="a" prints the C0 and C1 coefficents for calibration purposes
  //

  TString opt = option; opt.ToLower();
  printf("%s\n",GetTitle());
  printf(" - ROD shifts  (residual voltage settings): 40V correspond to 1mm of shift\n");
  Int_t count = 0;
  printf("  : A-side - (Rod & Strips) \n     "); 
  for (Int_t i=0; i<36;i++) {
    if (fRodVoltShiftA[i]!=0) {
      printf("Rod%2d:%3.1fV ",i,fRodVoltShiftA[i]);
      count++;
    }
    if (count==0&&i==35) 
      printf("-> all at 0 V\n");
    else if (i==35){
      printf("\n");
      count=0;
    }
  } 
  printf("  : C-side - (Rod & Strips) \n     "); 
  for (Int_t i=0; i<36;i++) {
    if (fRodVoltShiftC[i]!=0) {
      printf("Rod%2d:%3.1fV ",i,fRodVoltShiftC[i]);
      count++;
    }
    if (count==0&&i==35) 
      printf("-> all at 0 V\n");
    else if (i==35){
      printf("\n");
      count=0;
    }
  } 
 
  printf(" - Rotated clips (residual voltage settings): 40V correspond to 1mm of 'offset'\n");
  if (fRotatedClipVoltA[0]!=0) {   printf("     A side (IFC): HV Rod: %3.1f V \n",fRotatedClipVoltA[0]); count++; }
  if (fRotatedClipVoltA[1]!=0) {   printf("     A side (OFC): HV Rod: %3.1f V \n",fRotatedClipVoltA[1]); count++; }
  if (fRotatedClipVoltC[0]!=0) {   printf("     C side (IFC): HV Rod: %3.1f V \n",fRotatedClipVoltC[0]); count++; }
  if (fRotatedClipVoltC[1]!=0) {   printf("     C side (OFC): HV Rod: %3.1f V \n",fRotatedClipVoltC[1]); count++; }
  if (count==0) 
    printf("     -> no rotated clips \n");

  count=0;
  printf(" - Copper ROD shifts (without strips):\n");
  printf("  : A-side - (Copper Rod shift) \n     "); 
  for (Int_t i=0; i<36;i++) {
    if (fCopperRodShiftA[i]!=0) {
      printf("Rod%2d:%3.1fV ",i,fCopperRodShiftA[i]);
      count++;
    }
    if (count==0&&i==35) 
      printf("-> all at 0 V\n");
    else if (i==35){
      printf("\n");
      count=0;
    }
  } 
  printf("  : C-side - (Copper Rod shift) \n     "); 
  for (Int_t i=0; i<36;i++) {
    if (fCopperRodShiftC[i]!=0) {
      printf("Rod%2d:%3.1fV ",i,fCopperRodShiftC[i]);
      count++;
    }
    if (count==0&&i==35) 
      printf("-> all at 0 V\n");
    else if (i==35){
      printf("\n");
      count=0;
    }
  } 

  if (opt.Contains("a")) { // Print all details
    printf(" - T1: %1.4f, T2: %1.4f \n",fT1,fT2);
    printf(" - C1: %1.4f, C0: %1.4f \n",fC1,fC0);
  }

  if (!fInitLookUp) AliError("Lookup table was not initialized! You should do InitFCVoltError3D() ...");

}
