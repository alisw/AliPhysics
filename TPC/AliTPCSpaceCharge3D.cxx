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
// AliTPCSpaceCharge3D class                                                //
// The class calculates the space point distortions due to a space charge   //
// effect ....                                                              //
// Method of calculation:                                                   //
// The analytical solution for the poisson problem in 3D (cylindrical coord)//
// is used in form of look up tables. PieceOfCake (POC) voxel were pre-     //
// calculated and can be sumed up (weighted) according to the what is needed// 
//                                                                          //
// The class allows "effective Omega Tau" corrections.                      // 
//                                                                          //
// NOTE: This class is not  capable of calculation z distortions due to     //
//       drift velocity change in dependence of the electric field!!!       //
//                                                                          //
// date: 19/06/2010                                                         //
// Authors: Stefan Rossegger                                                //
//                                                                          //
// Example usage:                                                           //
//////////////////////////////////////////////////////////////////////////////

#include "AliMagF.h"
#include "TGeoGlobalMagField.h"
#include "AliTPCcalibDB.h"
#include "AliTPCParam.h"
#include "AliLog.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TFile.h"
#include "TVector.h"
#include "TMatrix.h"
#include "TMatrixD.h"

#include "TMath.h"
#include "AliTPCROC.h"
#include "AliTPCSpaceCharge3D.h"

ClassImp(AliTPCSpaceCharge3D)

AliTPCSpaceCharge3D::AliTPCSpaceCharge3D()
  : AliTPCCorrection("SpaceCharge3D","Space Charge - 3D"),
    fC0(0.),fC1(0.),
    fCorrectionFactor(1.),
    fInitLookUp(kFALSE),
    fSCDataFileName((char*)"$(ALICE_ROOT)/TPC/Calib/maps/sc_3D_distribution_Sim.root"),
    fSCLookUpPOCsFileName((char*)"$(ALICE_ROOT)/TPC/Calib/maps/sc_3D_raw_18-18-26_17p-18p-25p-MN30.root")
{
  //
  // default constructor
  //

  // Array which will contain the solution according to the setted charge density distribution
  // see InitSpaceCharge3DDistortion() function
  for ( Int_t k = 0 ; k < kNPhi ; k++ ) {
    fLookUpErOverEz[k]   =  new TMatrixD(kNR,kNZ);  
    fLookUpEphiOverEz[k] =  new TMatrixD(kNR,kNZ);
    fLookUpDeltaEz[k]    =  new TMatrixD(kNR,kNZ); 
    fSCdensityDistribution[k] = new TMatrixD(kNR,kNZ);
  }

  printf("%s\n",fSCDataFileName);
  printf("%s\n",fSCLookUpPOCsFileName);
  SetSCDataFileName(fSCDataFileName);

}

AliTPCSpaceCharge3D::~AliTPCSpaceCharge3D() {
  //
  // default destructor
  //
 
  for ( Int_t k = 0 ; k < kNPhi ; k++ ) {
    delete fLookUpErOverEz[k];
    delete fLookUpEphiOverEz[k];
    delete fLookUpDeltaEz[k];
    delete fSCdensityDistribution[k];
  }

  //  delete fSCdensityDistribution;
}



void AliTPCSpaceCharge3D::Init() {
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

  InitSpaceCharge3DDistortion(); // fill the look up table
}

void AliTPCSpaceCharge3D::Update(const TTimeStamp &/*timeStamp*/) {
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



void AliTPCSpaceCharge3D::GetCorrection(const Float_t x[],const Short_t roc,Float_t dx[]) {
  //
  // Calculates the correction due the Space Charge effect within the TPC drift volume
  //   

  if (!fInitLookUp) {
    AliInfo("Lookup table was not initialized! Performing the inizialisation now ...");
    InitSpaceCharge3DDistortion();
  }

  Int_t   order     = 1 ;    // FIXME: hardcoded? Linear interpolation = 1, Quadratic = 2         
                        
  Double_t intEr, intEphi, intdEz ;
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

  // Get the Er and Ephi field integrals plus the integral over DeltaEz 
  intEr      = Interpolate3DTable(order, r, z, phi, kNR, kNZ, kNPhi, 
				  fgkRList, fgkZList, fgkPhiList, fLookUpErOverEz  );
  intEphi    = Interpolate3DTable(order, r, z, phi, kNR, kNZ, kNPhi, 
				  fgkRList, fgkZList, fgkPhiList, fLookUpEphiOverEz);
  intdEz = Interpolate3DTable(order, r, z, phi, kNR, kNZ, kNPhi, 
				  fgkRList, fgkZList, fgkPhiList, fLookUpDeltaEz   );

  // Calculate distorted position
  if ( r > 0.0 ) {
    phi =  phi + fCorrectionFactor *( fC0*intEphi - fC1*intEr ) / r;      
    r   =  r   + fCorrectionFactor *( fC0*intEr   + fC1*intEphi );  
  }
  Double_t dz = intdEz * fCorrectionFactor * fgkdvdE;
 
  // Calculate correction in cartesian coordinates
  dx[0] = r * TMath::Cos(phi) - x[0];
  dx[1] = r * TMath::Sin(phi) - x[1]; 
  dx[2] = dz;  // z distortion - (scaled with driftvelocity dependency on the Ez field and the overall scaling factor)

}


void AliTPCSpaceCharge3D::InitSpaceCharge3DDistortion() {
  //
  // Initialization of the Lookup table which contains the solutions of the 
  // "space charge" (poisson) problem 
  //
  // The sum-up uses a look-up table which contains different discretized Space charge fields 
  // in order to calculate the corresponding field deviations due to a given (discretized)
  // space charge distribution ....
  //
  // Method of calculation: Weighted sum up of the different fields within the look up table
  //

  if (fInitLookUp) {
    AliInfo("Lookup table was already initialized!");
    //    return;
  }

  AliInfo("Preparation of the weighted look-up table");
   
  TFile *f = new TFile(fSCLookUpPOCsFileName);
  if (!f) { 
    AliError("Precalculated POC-looup-table could not be found");
    return;
  }

  // units are in [m]
  TVector *gridf  = (TVector*) f->Get("constants"); 
  TVector &grid = *gridf;
  TMatrix *coordf  = (TMatrix*) f->Get("coordinates");
  TMatrix &coord = *coordf;
  TMatrix *coordPOCf  = (TMatrix*) f->Get("POCcoord");
  TMatrix &coordPOC = *coordPOCf;
  
  Bool_t flagRadSym = 0;
  if (grid(1)==1 && grid(4)==1) {
    AliInfo("LOOK UP TABLE IS RADIAL SYMETTRIC - Field in Phi is ZERO");
    flagRadSym=1;
  }

  Int_t rows      = (Int_t)grid(0);   // number of points in r direction 
  Int_t phiSlices = (Int_t)grid(1);   // number of points in phi         
  Int_t columns   = (Int_t)grid(2);   // number of points in z direction 

  const Float_t gridSizeR   =  grid(6)*100;    // unit in [cm]
  const Float_t gridSizePhi   =  grid(7);  // unit in [rad]
  const Float_t gridSizeZ =  grid(8)*100;      // unit in [cm]
 
  // temporary matrices needed for the calculation
  TMatrixD *arrayofErA[phiSlices], *arrayofEphiA[phiSlices], *arrayofdEzA[phiSlices];
  TMatrixD *arrayofErC[phiSlices], *arrayofEphiC[phiSlices], *arrayofdEzC[phiSlices];

  TMatrixD *arrayofEroverEzA[phiSlices], *arrayofEphioverEzA[phiSlices], *arrayofDeltaEzA[phiSlices];
  TMatrixD *arrayofEroverEzC[phiSlices], *arrayofEphioverEzC[phiSlices], *arrayofDeltaEzC[phiSlices];

 
  for ( Int_t k = 0 ; k < phiSlices ; k++ ) {
   
    arrayofErA[k]   =   new TMatrixD(rows,columns) ;
    arrayofEphiA[k] =   new TMatrixD(rows,columns) ; // zero if radial symmetric
    arrayofdEzA[k]  =   new TMatrixD(rows,columns) ;
    arrayofErC[k]   =   new TMatrixD(rows,columns) ;
    arrayofEphiC[k] =   new TMatrixD(rows,columns) ; // zero if radial symmetric
    arrayofdEzC[k]  =   new TMatrixD(rows,columns) ;

    arrayofEroverEzA[k]   =   new TMatrixD(rows,columns) ;
    arrayofEphioverEzA[k] =   new TMatrixD(rows,columns) ; // zero if radial symmetric
    arrayofDeltaEzA[k]    =   new TMatrixD(rows,columns) ;
    arrayofEroverEzC[k]   =   new TMatrixD(rows,columns) ;
    arrayofEphioverEzC[k] =   new TMatrixD(rows,columns) ; // zero if radial symmetric
    arrayofDeltaEzC[k]    =   new TMatrixD(rows,columns) ;

    // Set the values to zero the lookup tables
    // not necessary, it is done in the constructor of TMatrix - code deleted

  }
 
  // list of points as used in the interpolation (during sum up)
  Double_t  rlist[rows], zedlist[columns] , philist[phiSlices];
  for ( Int_t k = 0 ; k < phiSlices ; k++ ) {
    philist[k] =  gridSizePhi * k;
    for ( Int_t i = 0 ; i < rows ; i++ )    {
      rlist[i] = fgkIFCRadius + i*gridSizeR ;
      for ( Int_t j = 0 ; j < columns ; j++ ) { 
	zedlist[j]  = j * gridSizeZ ;
      }
    }
  } // only done once
  
  
  // Prepare summation Vectors
  //  Int_t numPosInd = (Int_t)(grid(0)*grid(1)*grid(2)); // Number of fulcrums in the look-up table
  //  Int_t numPOC = (Int_t)(grid(3)*grid(4)*grid(5));    // Number of POC conf. in the look-up table
  
  TTree *treePOC = (TTree*)f->Get("POCall");

  TVector *bEr  = 0;   TVector *bEphi= 0;   TVector *bEz  = 0;
  
  treePOC->SetBranchAddress("Er",&bEr);
  if (!flagRadSym) treePOC->SetBranchAddress("Ephi",&bEphi);
  treePOC->SetBranchAddress("Ez",&bEz);

  /*  TBranch *b2Er = treePOC->GetBranch("Er");
      TBranch *b2Ephi =0;
      if (!flagRadSym) b2Ephi = treePOC->GetBranch("Ephi");
      TBranch *b2Ez = treePOC->GetBranch("Ez");
  */

  // Read the complete tree and do a weighted sum-up over the POC configurations
  // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  Int_t treeNumPOC = (Int_t)treePOC->GetEntries(); // Number of POC conf. in the look-up table
  Int_t ipC = 0; // POC Conf. counter (note: different to the POC number in the tree!)

  for (Int_t itreepC=0; itreepC<treeNumPOC; itreepC++) { // loop over POC configurations in tree
  
    // if (!flagRadSym) cout<<ipC<<endl;

    treePOC->GetEntry(itreepC);

    /*    b2Er->GetEntry(itreepC); 
	  if (!flagRadSym) b2Ephi->GetEntry(itreepC); 
	  b2Ez->GetEntry(itreepC); 
    */

    // center of the POC voxel in [meter]
    Double_t r0 = coordPOC(ipC,0);
    Double_t phi0 = coordPOC(ipC,1);
    Double_t z0 = coordPOC(ipC,2);

    ipC++; // POC conf. counter

    // weights (charge density) at POC position on the A and C side (in C/m^3/e0)
    // note: coordinates are in [cm]
    Double_t weightA = GetSpaceChargeDensity(r0*100,phi0, z0*100); 
    Double_t weightC = GetSpaceChargeDensity(r0*100,phi0,-z0*100);
  
    // Summing up the vector components according to their weight

    Int_t ip = 0;
    for ( Int_t j = 0 ; j < columns ; j++ ) { 
      for ( Int_t i = 0 ; i < rows ; i++ )    {
	for ( Int_t k = 0 ; k < phiSlices ; k++ ) {
		  
	  // check wether the coordinates were screwed
	  if ((coord(0,ip)*100-rlist[i])>0.01 || 
	      (coord(1,ip)-philist[k])>0.01 || 
	      (coord(2,ip)*100-zedlist[j])>0.01) {
	    AliError("internal error: coordinate system was screwed during the sum-up");
	    printf("lookup: (r,phi,z)=(%lf,%lf,%lf)\n",coord(0,ip)*100,coord(1,ip),coord(2,ip)*100);
	    printf("sum-up: (r,phi,z)=(%lf,%lf,%lf)\n",rlist[i],philist[k],zedlist[j]);
	    AliError("Don't trust the results of the space charge calibration!");
	  }
	  
	  // unfortunately, the lookup tables were produced to be faster for phi symmetric charges
	  // This will be the most frequent usage (hopefully)
	  // That's why we have to do this here ...

	  TMatrixD &erA   =  *arrayofErA[k]  ;
	  TMatrixD &ephiA =  *arrayofEphiA[k];
	  TMatrixD &dEzA    =  *arrayofdEzA[k]   ;
   
	  TMatrixD &erC   =  *arrayofErC[k]  ;
	  TMatrixD &ephiC =  *arrayofEphiC[k];
	  TMatrixD &dEzC    =  *arrayofdEzC[k]   ;
   
	  // Sum up - Efield values in [V/m] -> transition to [V/cm]
	  erA(i,j) += ((*bEr)(ip)) * weightA /100;
	  erC(i,j) += ((*bEr)(ip)) * weightC /100;
	  if (!flagRadSym) {
	    ephiA(i,j) += ((*bEphi)(ip)) * weightA; // [V/rad]
	    ephiC(i,j) += ((*bEphi)(ip)) * weightC; // [V/rad]
	  }
	  dEzA(i,j)  += ((*bEz)(ip)) * weightA /100;
	  dEzC(i,j)  += ((*bEz)(ip)) * weightC /100;

	  // increase the counter
	  ip++;
	}
      }
    }
    
    // Rotation and summation in the rest of the dPhiSteps
    // which were not stored in the tree due to storage & symmetry reasons

    Int_t phiPoints = (Int_t) grid(1);
    Int_t phiPOC    = (Int_t) grid(4);
    
    // printf("%d %d\n",phiPOC,flagRadSym);
    
    for (Int_t phiiC = 1; phiiC<phiPOC; phiiC++) { // just used for non-radial symetric table
      
      r0 = coordPOC(ipC,0);
      phi0 = coordPOC(ipC,1);
      z0 = coordPOC(ipC,2);
      
      ipC++; // POC conf. counter
      
      // weights (charge density) at POC position on the A and C side (in C/m^3/e0)
      // note: coordinates are in [cm]
      weightA = GetSpaceChargeDensity(r0*100,phi0, z0*100); 
      weightC = GetSpaceChargeDensity(r0*100,phi0,-z0*100);
      
      
      // Summing up the vector components according to their weight
      ip = 0;
      for ( Int_t j = 0 ; j < columns ; j++ ) { 
	for ( Int_t i = 0 ; i < rows ; i++ )    {
	  for ( Int_t k = 0 ; k < phiSlices ; k++ ) {
	    
	    // Note: rotating the coordinated during the sum up
	    
	    Int_t rotVal = (phiPoints/phiPOC)*phiiC;
	    Int_t ipR = -1;
	    
	    if ((ip%phiPoints)>=rotVal) {
	      ipR = ip-rotVal;
	    } else {
	      ipR = ip+(phiPoints-rotVal);
	    }
	    
	    // unfortunately, the lookup tables were produced to be faster for phi symmetric charges
	    // This will be the most frequent usage (hopefully)
	    // That's why we have to do this here and not outside the loop ...
	    
	    TMatrixD &erA   =  *arrayofErA[k]  ;
	    TMatrixD &ephiA =  *arrayofEphiA[k];
	    TMatrixD &dEzA  =  *arrayofdEzA[k]   ;
	    
	    TMatrixD &erC   =  *arrayofErC[k]  ;
	    TMatrixD &ephiC =  *arrayofEphiC[k];
	    TMatrixD &dEzC  =  *arrayofdEzC[k]   ;
       
	    // Sum up - Efield values in [V/m] -> transition to [V/cm]
	    erA(i,j) += ((*bEr)(ipR)) * weightA /100;
	    erC(i,j) += ((*bEr)(ipR)) * weightC /100;
	    if (!flagRadSym) {
	      ephiA(i,j) += ((*bEphi)(ipR)) * weightA; // [V/rad]
	      ephiC(i,j) += ((*bEphi)(ipR)) * weightC; // [V/rad]
	    }
	    dEzA(i,j)  += ((*bEz)(ipR)) * weightA /100;
	    dEzC(i,j)  += ((*bEz)(ipR)) * weightC /100;

	    // increase the counter
	    ip++;
	  }
	}
      } // end coordinate loop

    } // end phi-POC summation (phiiC)
   

    // printf("POC: (r,phi,z) = (%lf %lf %lf) | weight(A,C): %03.1lf %03.1lf\n",r0,phi0,z0, weightA, weightC);

  } // end POC loop
 
  AliInfo("Division and integration");

  // -------------------------------------------------------------------------------
  // Division by the Ez (drift) field and integration along z

  Double_t ezField = (fgkCathodeV-fgkGG)/fgkTPCZ0; // = Electric Field (V/cm) Magnitude ~ -400 V/cm;

  for ( Int_t k = 0 ; k < phiSlices ; k++ ) { // phi loop

    // matrices holding the solution - summation of POC charges
    TMatrixD &erA   =  *arrayofErA[k]  ;
    TMatrixD &ephiA =  *arrayofEphiA[k];
    TMatrixD &dezA  =  *arrayofdEzA[k]   ;
    TMatrixD &erC   =  *arrayofErC[k]  ;
    TMatrixD &ephiC =  *arrayofEphiC[k];
    TMatrixD &dezC  =  *arrayofdEzC[k]   ;

    // matrices which contain the integrated fields (divided by the drift field)
    TMatrixD &erOverEzA   =  *arrayofEroverEzA[k]  ;
    TMatrixD &ephiOverEzA =  *arrayofEphioverEzA[k];
    TMatrixD &deltaEzA    =  *arrayofDeltaEzA[k];
    TMatrixD &erOverEzC   =  *arrayofEroverEzC[k]  ;
    TMatrixD &ephiOverEzC =  *arrayofEphioverEzC[k];
    TMatrixD &deltaEzC    =  *arrayofDeltaEzC[k];    
    
    for ( Int_t i = 0 ; i < rows ; i++ )    { // r loop
      for ( Int_t j = columns-1 ; j >= 0 ; j-- ) {// z loop 
	// Count backwards to facilitate integration over Z

	Int_t index = 1 ; // Simpsons rule if N=odd.If N!=odd then add extra point by trapezoidal rule.  

	erOverEzA(i,j) = 0; ephiOverEzA(i,j) = 0;  deltaEzA(i,j) = 0;
	erOverEzC(i,j) = 0; ephiOverEzC(i,j) = 0;  deltaEzC(i,j) = 0;

	for ( Int_t m = j ; m < columns ; m++ ) { // integration

	  erOverEzA(i,j)   += index*(gridSizeZ/3.0)*erA(i,m)/(-1*ezField) ;
	  erOverEzC(i,j)   += index*(gridSizeZ/3.0)*erC(i,m)/(-1*ezField)  ;
	  if (!flagRadSym) {
	    ephiOverEzA(i,j) += index*(gridSizeZ/3.0)*ephiA(i,m)/(-1*ezField)  ;
	    ephiOverEzC(i,j) += index*(gridSizeZ/3.0)*ephiC(i,m)/(-1*ezField)  ;
	  }
	  deltaEzA(i,j)    += index*(gridSizeZ/3.0)*dezA(i,m) ;
	  deltaEzC(i,j)    += index*(gridSizeZ/3.0)*dezC(i,m) ;

	  if ( index != 4 )  index = 4; else index = 2 ;

	}

	if ( index == 4 ) {
	  erOverEzA(i,j)   -= (gridSizeZ/3.0)*erA(i,columns-1)/(-1*ezField) ;
	  erOverEzC(i,j)   -= (gridSizeZ/3.0)*erC(i,columns-1)/(-1*ezField) ;
	  if (!flagRadSym) {
	    ephiOverEzA(i,j) -= (gridSizeZ/3.0)*ephiA(i,columns-1)/(-1*ezField) ;
	    ephiOverEzC(i,j) -= (gridSizeZ/3.0)*ephiC(i,columns-1)/(-1*ezField) ;
	  }
	  deltaEzA(i,j)    -= (gridSizeZ/3.0)*dezA(i,columns-1) ;
	  deltaEzC(i,j)    -= (gridSizeZ/3.0)*dezC(i,columns-1) ;
	}
	if ( index == 2 ) {
	  erOverEzA(i,j)   += (gridSizeZ/3.0)*(0.5*erA(i,columns-2)-2.5*erA(i,columns-1))/(-1*ezField) ;
	  erOverEzC(i,j)   += (gridSizeZ/3.0)*(0.5*erC(i,columns-2)-2.5*erC(i,columns-1))/(-1*ezField) ;
	  if (!flagRadSym) {
	    ephiOverEzA(i,j) += (gridSizeZ/3.0)*(0.5*ephiA(i,columns-2)-2.5*ephiA(i,columns-1))/(-1*ezField) ;
	    ephiOverEzC(i,j) += (gridSizeZ/3.0)*(0.5*ephiC(i,columns-2)-2.5*ephiC(i,columns-1))/(-1*ezField) ;
	  }
	  deltaEzA(i,j)    += (gridSizeZ/3.0)*(0.5*dezA(i,columns-2)-2.5*dezA(i,columns-1)) ;
	  deltaEzC(i,j)    += (gridSizeZ/3.0)*(0.5*dezC(i,columns-2)-2.5*dezC(i,columns-1)) ;
	}
	if ( j == columns-2 ) {
	  erOverEzA(i,j)   = (gridSizeZ/3.0)*(1.5*erA(i,columns-2)+1.5*erA(i,columns-1))/(-1*ezField) ;
	  erOverEzC(i,j)   = (gridSizeZ/3.0)*(1.5*erC(i,columns-2)+1.5*erC(i,columns-1))/(-1*ezField) ;
	  if (!flagRadSym) {
	    ephiOverEzA(i,j) = (gridSizeZ/3.0)*(1.5*ephiA(i,columns-2)+1.5*ephiA(i,columns-1))/(-1*ezField) ;
	    ephiOverEzC(i,j) = (gridSizeZ/3.0)*(1.5*ephiC(i,columns-2)+1.5*ephiC(i,columns-1))/(-1*ezField) ;
	  }
	  deltaEzA(i,j)    = (gridSizeZ/3.0)*(1.5*dezA(i,columns-2)+1.5*dezA(i,columns-1)) ;
	  deltaEzC(i,j)    = (gridSizeZ/3.0)*(1.5*dezC(i,columns-2)+1.5*dezC(i,columns-1)) ;
	}
	if ( j == columns-1 ) {
	  erOverEzA(i,j)   = 0;    erOverEzC(i,j)   = 0;
	  if (!flagRadSym) {
	    ephiOverEzA(i,j) = 0;  ephiOverEzC(i,j) = 0;
	  }
	  deltaEzA(i,j)    = 0;    deltaEzC(i,j)    = 0;
	}
      }
    }

  }
  
  AliInfo("Interpolation to Standard grid");

  // -------------------------------------------------------------------------------
  // Interpolate results onto the standard grid which is used for all AliTPCCorrections

  const Int_t order  = 1  ;  // Linear interpolation = 1, Quadratic = 2  

  Double_t  r, phi, z ;
  for ( Int_t k = 0 ; k < kNPhi ; k++ ) {
    phi = fgkPhiList[k] ;
	
    TMatrixD &erOverEz   =  *fLookUpErOverEz[k]  ;
    TMatrixD &ephiOverEz =  *fLookUpEphiOverEz[k];
    TMatrixD &deltaEz    =  *fLookUpDeltaEz[k]   ;
	
    for ( Int_t j = 0 ; j < kNZ ; j++ ) {

      z = TMath::Abs(fgkZList[j]) ;  // z position is symmetric
    
      for ( Int_t i = 0 ; i < kNR ; i++ ) { 
	r = fgkRList[i] ;

	// Interpolate Lookup tables onto standard grid
	if (fgkZList[j]>0) {
	  erOverEz(i,j)   = Interpolate3DTable(order, r, z, phi, rows, columns, phiSlices, 
					       rlist, zedlist, philist, arrayofEroverEzA  );
	  ephiOverEz(i,j) = Interpolate3DTable(order, r, z, phi, rows, columns, phiSlices,
					       rlist, zedlist, philist, arrayofEphioverEzA);
	  deltaEz(i,j)    = Interpolate3DTable(order, r, z, phi, rows, columns, phiSlices,
					       rlist, zedlist, philist, arrayofDeltaEzA   );
	} else {
	  erOverEz(i,j)   = Interpolate3DTable(order, r, z, phi, rows, columns, phiSlices, 
					       rlist, zedlist, philist, arrayofEroverEzC  );
	  ephiOverEz(i,j) = Interpolate3DTable(order, r, z, phi, rows, columns, phiSlices,
					       rlist, zedlist, philist, arrayofEphioverEzC);
	  deltaEz(i,j)  = - Interpolate3DTable(order, r, z, phi, rows, columns, phiSlices,
					       rlist, zedlist, philist, arrayofDeltaEzC   );
	  // negative coordinate system on C side
	}

      } // end r loop
    } // end z loop
  } // end phi loop

  
  // clear the temporary arrays lists
  for ( Int_t k = 0 ; k < phiSlices ; k++ )  {

    delete arrayofErA[k];  
    delete arrayofEphiA[k];
    delete arrayofdEzA[k];
    delete arrayofErC[k];  
    delete arrayofEphiC[k];
    delete arrayofdEzC[k];

    delete arrayofEroverEzA[k];  
    delete arrayofEphioverEzA[k];
    delete arrayofDeltaEzA[k];
    delete arrayofEroverEzC[k];  
    delete arrayofEphioverEzC[k];
    delete arrayofDeltaEzC[k];

  }


  fInitLookUp = kTRUE;

}


void AliTPCSpaceCharge3D::SetSCDataFileName(char *const fname) {
  //
  // Set & load the Space charge density distribution from a file 
  // (linear interpolation onto a standard grid)
  //

  fSCDataFileName = fname;

  TFile *f = new TFile(fSCDataFileName,"READ");
  if (!f) { 
    AliError(Form("File %s, which should contain the space charge distribution, could not be found",
		  fSCDataFileName));
    return;
  }
 
  TH3F *density = (TH3F*) f->Get("SpaceChargeDistribution");
  if (!density) { 
    AliError(Form("The indicated file (%s) does not contain a histogram called %s",
		  fSCDataFileName,"SpaceChargeDistribution"));
    return;
  }
 

  Double_t  r, phi, z ;
  for ( Int_t k = 0 ; k < kNPhi ; k++ ) {
    phi = fgkPhiList[k] ;
    TMatrixD &scdensity   =  *fSCdensityDistribution[k]  ;
    for ( Int_t j = 0 ; j < kNZ ; j++ ) {
      z = TMath::Abs(fgkZList[j]) ;  // z position is symmetric
      for ( Int_t i = 0 ; i < kNR ; i++ ) { 
	r = fgkRList[i] ;
	scdensity(i,j) = density->Interpolate(r,phi,z); // quite slow
	//  	printf("%lf %lf %lf: %lf\n",r,phi,z,scdensity(i,j));
      }
    }
  }

 
  f->Close();

  fInitLookUp = kFALSE;

  
}


Float_t  AliTPCSpaceCharge3D::GetSpaceChargeDensity(Float_t r, Float_t phi, Float_t z) {
  //
  // returns the (input) space charge density at a given point according 
  // Note: input in [cm], output in [C/m^3/e0] !!
  //

  if (!fSCdensityDistribution) {
    printf("Scheiss is nuul - argg\n");
    return 0.;
  }

  while (phi<0) phi += TMath::TwoPi();
  while (phi>TMath::TwoPi()) phi -= TMath::TwoPi();


  // Float_t sc =fSCdensityDistribution->Interpolate(r0,phi0,z0);
  Int_t order = 1; //
  Float_t sc = Interpolate3DTable(order, r, z, phi, kNR, kNZ, kNPhi, 
				  fgkRList, fgkZList, fgkPhiList, fSCdensityDistribution );

  //  printf("%lf %lf %lf: %lf\n",r,phi,z,sc);
  
  return sc;
}


TH2F * AliTPCSpaceCharge3D::CreateHistoSCinXY(Float_t z, Int_t nx, Int_t ny) {
  //
  // return a simple histogramm containing the space charge distribution (input for the calculation)
  //

  TH2F *h=CreateTH2F("spaceCharge",GetTitle(),"x [cm]","y [cm]","#rho_{sc} [C/m^{3}/e_{0}]",
		     nx,-250.,250.,ny,-250.,250.);

  for (Int_t iy=1;iy<=ny;++iy) {
    Double_t yp = h->GetYaxis()->GetBinCenter(iy);
    for (Int_t ix=1;ix<=nx;++ix) {
      Double_t xp = h->GetXaxis()->GetBinCenter(ix);
    
      Float_t r = TMath::Sqrt(xp*xp+yp*yp);
      Float_t phi = TMath::ATan2(yp,xp); 
   
      if (85.<=r && r<=250.) {
	Float_t sc = GetSpaceChargeDensity(r,phi,z)/fgke0; // in [C/m^3/e0]
	h->SetBinContent(ix,iy,sc); 
      } else {
	h->SetBinContent(ix,iy,0.);
      }
    }
  }
  
  return h;
} 

TH2F * AliTPCSpaceCharge3D::CreateHistoSCinZR(Float_t phi, Int_t nz, Int_t nr) {
  //
  // return a simple histogramm containing the space charge distribution (input for the calculation)
  //

  TH2F *h=CreateTH2F("spaceCharge",GetTitle(),"z [cm]","r [cm]","#rho_{sc} [C/m^{3}/e_{0}]",
		     nz,-250.,250.,nr,85.,250.);

  for (Int_t ir=1;ir<=nr;++ir) {
    Float_t r = h->GetYaxis()->GetBinCenter(ir);
    for (Int_t iz=1;iz<=nz;++iz) {
      Float_t z = h->GetXaxis()->GetBinCenter(iz);
      Float_t sc = GetSpaceChargeDensity(r,phi,z)/fgke0; // in [C/m^3/e0]
      h->SetBinContent(iz,ir,sc);
    }
  }

  return h;
} 

void AliTPCSpaceCharge3D::WriteChargeDistributionToFile(const char* fname) {
  //
  // Example on how to write a Space charge distribution into a File
  //  (see below: estimate from scaling STAR measurements to Alice)
  //

  TFile *f = new TFile(fname,"RECREATE");

  // some grid, not too course
  Int_t nr = 72;
  Int_t nphi = 180;
  Int_t nz = 250;

  Double_t dr = (fgkOFCRadius-fgkIFCRadius)/(nr+1);
  Double_t dphi = TMath::TwoPi()/(nphi+1);
  Double_t dz = 500./(nz+1);
  Double_t safty = 0.; // due to a root bug which does not interpolate the boundary (first and last bin) correctly

  TH3F *histo = new TH3F("charge","charge",
			 nr,fgkIFCRadius-dr-safty,fgkOFCRadius+dr+safty,
			 nphi,0-dphi-safty,TMath::TwoPi()+dphi+safty,
			 nz,-250-dz-safty,250+dz+safty);
 
  for (Int_t ir=1;ir<=nr;++ir) {
    Double_t rp = histo->GetXaxis()->GetBinCenter(ir);
    for (Int_t iphi=1;iphi<=nphi;++iphi) {
      Double_t phip = histo->GetYaxis()->GetBinCenter(iphi);
      for (Int_t iz=1;iz<=nz;++iz) {
	Double_t zp = histo->GetZaxis()->GetBinCenter(iz);


	zp = TMath::Abs(zp); // symmetric on both sides of the TPC
	
	// recalculation to meter
	Double_t lZ = 2.5; // approx. TPC drift length
	Double_t rpM = rp/100.; // in [m]
	Double_t zpM = zp/100.; // in [m]
	
	// setting of mb multiplicity and Interaction rate
	Double_t multiplicity = 950;
	Double_t intRate = 8000;

	// calculation of "scaled" parameters
	Double_t a = multiplicity*intRate/79175;
	Double_t b = a/lZ;
	
	Double_t charge = (a - b * zpM)/(rpM*rpM); // charge in [C/m^3/e0]
	
	// some 'arbitrary' GG leaks
	if (  (phip<(TMath::Pi()/9*3) && phip>(TMath::Pi()/9*2) ) ) { //||
	  //	      (phip<(TMath::Pi()/9*7) && phip>(TMath::Pi()/9*6) ) ||
	  //	      (phip<(TMath::Pi()/9*16) && phip>(TMath::Pi()/9*13) ) ){
	  if (rp>130&&rp<135) 
	    charge = 300; 
	}

	//	printf("%lf %lf %lf: %lf\n",rp,phip,zp,charge);
 
	charge = charge*fgke0; // [C/m^3]

	histo->SetBinContent(ir,iphi,iz,charge); 
      }
    }
  }


  histo->Write("SpaceChargeDistribution");
  f->Close();
  
}


void AliTPCSpaceCharge3D::Print(const Option_t* option) const {
  //
  // Print function to check the settings of the boundary vectors
  // option=="a" prints the C0 and C1 coefficents for calibration purposes
  //

  TString opt = option; opt.ToLower();
  printf("%s\n",GetTitle());
  printf(" - Space Charge effect with arbitrary 3D charge density (as input).\n");
  printf("   SC correction factor: %f \n",fCorrectionFactor);

  if (opt.Contains("a")) { // Print all details
    printf(" - T1: %1.4f, T2: %1.4f \n",fT1,fT2);
    printf(" - C1: %1.4f, C0: %1.4f \n",fC1,fC0);
  } 
   
  if (!fInitLookUp) AliError("Lookup table was not initialized! You should do InitSpaceCharge3DDistortion() ...");

}
