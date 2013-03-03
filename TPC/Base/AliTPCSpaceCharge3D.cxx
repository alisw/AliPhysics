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
//   <h2>  AliTPCSpaceCharge3D class   </h2>    
//   The class calculates the space point distortions due to an arbitrary space
//   charge distribution in 3D. 
//   <p>
//   The method of calculation is based on the analytical solution for the Poisson 
//   problem in 3D (cylindrical coordinates). The solution is used in form of 
//   look up tables, where the pre calculated solutions for different voxel 
//   positions are stored. These voxel solutions can be summed up according 
//   to the weight of the position of the applied space charge distribution.
//   Further details can be found in \cite[chap.5]{PhD-thesis_S.Rossegger}.
//   <p>
//   The class also allows a simple scaling of the resulting distortions
//   via the function SetCorrectionFactor. This speeds up the calculation 
//   time under the assumption, that the distortions scales linearly with the 
//   magnitude of the space charge distribution $\rho(r,z)$ and the shape stays 
//   the same at higher luminosities.
//   <p>
//   In contrast to the implementation in 2D (see the class AliTPCSpaceChargeabove), 
//   the input charge distribution can be of arbitrary character. An example on how 
//   to produce a corresponding charge distribution can be found in the function 
//   WriteChargeDistributionToFile. In there, a $\rho(r,z) = (A-B\,z)/r^2$, 
//   with slightly different magnitude on the A and C side (due to the muon absorber),
//   is superpositioned with a few leaking wires at arbitrary positions. 
// End_Html
//
// Begin_Macro(source)
//   {
//   gROOT->SetStyle("Plain"); gStyle->SetPalette(1);
//   TCanvas *c2 = new TCanvas("cAliTPCSpaceCharge3D","cAliTPCSpaceCharge3D",500,400); 
//   AliTPCSpaceCharge3D sc;
//   sc.WriteChargeDistributionToFile("SC_zr2_GGleaks.root");
//   sc.SetSCDataFileName("SC_zr2_GGleaks.root");
//   sc.SetOmegaTauT1T2(0,1,1); // B=0
//   sc.InitSpaceCharge3DDistortion();
//   sc.CreateHistoDRinXY(15,300,300)->Draw("colz");
//   return c2;
//   } 
// End_Macro
//
// Begin_Html
//   <p>
//   Date: 19/06/2010  <br>                                                       
//   Authors: Stefan Rossegger                                                
// End_Html 
// _________________________________________________________________


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
    fSCDataFileName(""),
    fSCLookUpPOCsFileName3D(""),
    fSCLookUpPOCsFileNameRZ(""),
    fSCLookUpPOCsFileNameRPhi(""),
    fSCdensityInRZ(0),
    fSCdensityInRPhiA(0), 
    fSCdensityInRPhiC(0)
{
  //
  // default constructor
  //

  // Array which will contain the solution according to the setted charge density distribution
  // see InitSpaceCharge3DDistortion() function
  for ( Int_t k = 0 ; k < kNPhi ; k++ ) {
    fLookUpErOverEz[k]   =  new TMatrixF(kNR,kNZ);  
    fLookUpEphiOverEz[k] =  new TMatrixF(kNR,kNZ);
    fLookUpDeltaEz[k]    =  new TMatrixF(kNR,kNZ); 
    fSCdensityDistribution[k] = new TMatrixF(kNR,kNZ);
  }
  fSCdensityInRZ   = new TMatrixD(kNR,kNZ);
  fSCdensityInRPhiA = new TMatrixD(kNR,kNPhi);
  fSCdensityInRPhiC = new TMatrixD(kNR,kNPhi);

  // location of the precalculated look up tables

  fSCLookUpPOCsFileName3D="$(ALICE_ROOT)/TPC/Calib/maps/sc_3D_raw_18-18-26_17p-18p-25p-MN30.root"; // rough estimate
  fSCLookUpPOCsFileNameRZ="$(ALICE_ROOT)/TPC/Calib/maps/sc_radSym_35-01-51_34p-01p-50p_MN60.root";
  fSCLookUpPOCsFileNameRPhi="$(ALICE_ROOT)/TPC/Calib/maps/sc_cChInZ_35-144-26_34p-18p-01p-MN30.root";
  //  fSCLookUpPOCsFileNameRPhi="$(ALICE_ROOT)/TPC/Calib/maps/sc_cChInZ_35-36-26_34p-18p-01p-MN40.root";
 


  // standard location of the space charge distibution ... can be changes
  fSCDataFileName="$(ALICE_ROOT)/TPC/Calib/maps/sc_3D_distribution_Sim.root";

  //  SetSCDataFileName(fSCDataFileName.Data()); // should be done by the user


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
  delete fSCdensityInRZ;
  delete fSCdensityInRPhiA;
  delete fSCdensityInRPhiC;

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
                        
  Float_t intEr, intEphi, intdEz ;
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
  dx[0] = - (r * TMath::Cos(phi) - x[0]);
  dx[1] = - (r * TMath::Sin(phi) - x[1]); 
  dx[2] = - dz;  // z distortion - (scaled with driftvelocity dependency on the Ez field and the overall scaling factor)

}

void AliTPCSpaceCharge3D::InitSpaceCharge3DDistortion() {
 //
  // Initialization of the Lookup table which contains the solutions of the 
  // "space charge" (poisson) problem - Faster and more accureate
  //
  // Method: Weighted sum-up of the different fields within the look up table 
  // but using two lookup tables with higher granularity in the (r,z) and the (rphi)- plane to emulate
  // more realistic space charges. (r,z) from primary ionisation. (rphi) for possible Gating leaks

  if (fInitLookUp) {
    AliInfo("Lookup table was already initialized!  Doing it again anyway ...");
    //    return;
  }
  
  // ------------------------------------------------------------------------------------------------------
  // step 1: lookup table in rz, fine grid, radial symetric, to emulate primary ionization

  AliInfo("Step 1: Preparation of the weighted look-up tables.");
   
  // lookup table in rz, fine grid

  TFile *fZR = new TFile(fSCLookUpPOCsFileNameRZ.Data(),"READ");
  if ( !fZR ) {
    AliError("Precalculated POC-looup-table in ZR could not be found");
    return;
  } 

  // units are in [m]
  TVector *gridf1  = (TVector*) fZR->Get("constants"); 
  TVector &grid1 = *gridf1;
  TMatrix *coordf1  = (TMatrix*) fZR->Get("coordinates");
  TMatrix &coord1 = *coordf1;
  TMatrix *coordPOCf1  = (TMatrix*) fZR->Get("POCcoord");
  TMatrix &coordPOC1 = *coordPOCf1;
  
  Int_t rows      = (Int_t)grid1(0);   // number of points in r direction  - from RZ or RPhi table
  Int_t phiSlices = (Int_t)grid1(1);   // number of points in phi          - from RPhi table
  Int_t columns   = (Int_t)grid1(2);   // number of points in z direction  - from RZ table

  Float_t gridSizeR   =  (fgkOFCRadius-fgkIFCRadius)/(rows-1);   // unit in [cm]
  Float_t gridSizeZ   =  fgkTPCZ0/(columns-1);                  // unit in [cm]

  // temporary matrices needed for the calculation // for rotational symmetric RZ table, phislices is 1
  
  TMatrixD *arrayofErA[kNPhiSlices], *arrayofdEzA[kNPhiSlices]; 
  TMatrixD *arrayofErC[kNPhiSlices], *arrayofdEzC[kNPhiSlices]; 

  TMatrixD *arrayofEroverEzA[kNPhiSlices], *arrayofDeltaEzA[kNPhiSlices];
  TMatrixD *arrayofEroverEzC[kNPhiSlices], *arrayofDeltaEzC[kNPhiSlices];

 
  for ( Int_t k = 0 ; k < phiSlices ; k++ ) {
   
    arrayofErA[k]   =   new TMatrixD(rows,columns) ;
    arrayofdEzA[k]  =   new TMatrixD(rows,columns) ;
    arrayofErC[k]   =   new TMatrixD(rows,columns) ;
    arrayofdEzC[k]  =   new TMatrixD(rows,columns) ;

    arrayofEroverEzA[k]   =   new TMatrixD(rows,columns) ;
    arrayofDeltaEzA[k]    =   new TMatrixD(rows,columns) ;
    arrayofEroverEzC[k]   =   new TMatrixD(rows,columns) ;
    arrayofDeltaEzC[k]    =   new TMatrixD(rows,columns) ;

    // zero initialization not necessary, it is done in the constructor of TMatrix 

  }
 
  // list of points as used during sum up
  Double_t  rlist1[kNRows], zedlist1[kNColumns];// , philist1[phiSlices];
  for ( Int_t i = 0 ; i < rows ; i++ )    {
    rlist1[i] = fgkIFCRadius + i*gridSizeR ;
    for ( Int_t j = 0 ; j < columns ; j++ ) { 
      zedlist1[j]  = j * gridSizeZ ;
    }
  }
  
  TTree *treePOC = (TTree*)fZR->Get("POCall");

  TVector *bEr  = 0;   //TVector *bEphi= 0;
  TVector *bEz  = 0;
  
  treePOC->SetBranchAddress("Er",&bEr);
  treePOC->SetBranchAddress("Ez",&bEz);


  // Read the complete tree and do a weighted sum-up over the POC configurations
  // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  Int_t treeNumPOC = (Int_t)treePOC->GetEntries(); // Number of POC conf. in the look-up table
  Int_t ipC = 0; // POC Conf. counter (note: different to the POC number in the tree!)

  for (Int_t itreepC=0; itreepC<treeNumPOC; itreepC++) { // ------------- loop over POC configurations in tree
  
    treePOC->GetEntry(itreepC);

    // center of the POC voxel in [meter]
    Double_t r0 = coordPOC1(ipC,0);
    Double_t phi0 = coordPOC1(ipC,1);
    Double_t z0 = coordPOC1(ipC,2);

    ipC++; // POC configuration counter

    // weights (charge density) at POC position on the A and C side (in C/m^3/e0)
    // note: coordinates are in [cm]
    Double_t weightA = GetSpaceChargeDensity(r0*100,phi0, z0*100, 1);  // partial load in r,z
    Double_t weightC = GetSpaceChargeDensity(r0*100,phi0,-z0*100, 1);  // partial load in r,z
  
    // Summing up the vector components according to their weight

    Int_t ip = 0;
    for ( Int_t j = 0 ; j < columns ; j++ ) { 
      for ( Int_t i = 0 ; i < rows ; i++ )    {
	for ( Int_t k = 0 ; k < phiSlices ; k++ ) {
		  
	  // check wether the coordinates were screwed
	  if (TMath::Abs((coord1(0,ip)*100-rlist1[i]))>1 || 
	      TMath::Abs((coord1(2,ip)*100-zedlist1[j])>1)) { 
	    AliError("internal error: coordinate system was screwed during the sum-up");
	    printf("sum-up: (r,z)=(%f,%f)\n",rlist1[i],zedlist1[j]);
	    printf("lookup: (r,z)=(%f,%f)\n",coord1(0,ip)*100,coord1(2,ip)*100);
	    AliError("Don't trust the results of the space charge calculation!");
	  }
	  
	  // unfortunately, the lookup tables were produced to be faster for phi symmetric charges
	  // This will be the most frequent usage (hopefully)
	  // That's why we have to do this here ...

	  TMatrixD &erA   =  *arrayofErA[k]  ;
	  TMatrixD &dEzA  =  *arrayofdEzA[k]   ;
   
	  TMatrixD &erC   =  *arrayofErC[k]  ;
	  TMatrixD &dEzC  =  *arrayofdEzC[k]   ;
   
	  // Sum up - Efield values in [V/m] -> transition to [V/cm]
	  erA(i,j) += ((*bEr)(ip)) * weightA /100;
	  erC(i,j) += ((*bEr)(ip)) * weightC /100;
	  dEzA(i,j)  += ((*bEz)(ip)) * weightA /100;
	  dEzC(i,j)  += ((*bEz)(ip)) * weightC /100;

	  // increase the counter
	  ip++;
	}
      }
    } // end coordinate loop
  } // end POC loop


  // -------------------------------------------------------------------------------
  // Division by the Ez (drift) field and integration along z

  //  AliInfo("Step 1: Division and integration");

  Double_t ezField = (fgkCathodeV-fgkGG)/fgkTPCZ0; // = Electric Field (V/cm) Magnitude ~ -400 V/cm;

  for ( Int_t k = 0 ; k < phiSlices ; k++ ) { // phi loop

    // matrices holding the solution - summation of POC charges // see above
    TMatrixD &erA   =  *arrayofErA[k]  ;
    TMatrixD &dezA  =  *arrayofdEzA[k]   ;
    TMatrixD &erC   =  *arrayofErC[k]  ;
    TMatrixD &dezC  =  *arrayofdEzC[k]   ;

    // matrices which will contain the integrated fields (divided by the drift field)
    TMatrixD &erOverEzA   =  *arrayofEroverEzA[k]  ;
    TMatrixD &deltaEzA    =  *arrayofDeltaEzA[k];
    TMatrixD &erOverEzC   =  *arrayofEroverEzC[k]  ;
    TMatrixD &deltaEzC    =  *arrayofDeltaEzC[k];    
    
    for ( Int_t i = 0 ; i < rows ; i++ )    { // r loop
      for ( Int_t j = columns-1 ; j >= 0 ; j-- ) {// z loop 
	// Count backwards to facilitate integration over Z

	Int_t index = 1 ; // Simpsons rule if N=odd.If N!=odd then add extra point 
	                  // by trapezoidal rule.  

	erOverEzA(i,j) = 0; //ephiOverEzA(i,j) = 0;
	deltaEzA(i,j) = 0;
	erOverEzC(i,j) = 0; //ephiOverEzC(i,j) = 0; 
	deltaEzC(i,j) = 0;

	for ( Int_t m = j ; m < columns ; m++ ) { // integration

	  erOverEzA(i,j)   += index*(gridSizeZ/3.0)*erA(i,m)/(-1*ezField) ;
	  erOverEzC(i,j)   += index*(gridSizeZ/3.0)*erC(i,m)/(-1*ezField)  ;
	  deltaEzA(i,j)    += index*(gridSizeZ/3.0)*dezA(i,m)/(-1) ;
	  deltaEzC(i,j)    += index*(gridSizeZ/3.0)*dezC(i,m)/(-1) ;

	  if ( index != 4 )  index = 4; else index = 2 ;

	}

	if ( index == 4 ) {
	  erOverEzA(i,j)   -= (gridSizeZ/3.0)*erA(i,columns-1)/(-1*ezField) ;
	  erOverEzC(i,j)   -= (gridSizeZ/3.0)*erC(i,columns-1)/(-1*ezField) ;
	  deltaEzA(i,j)    -= (gridSizeZ/3.0)*dezA(i,columns-1)/(-1) ;
	  deltaEzC(i,j)    -= (gridSizeZ/3.0)*dezC(i,columns-1)/(-1) ;
	}
	if ( index == 2 ) {
	  erOverEzA(i,j)   += (gridSizeZ/3.0)*(0.5*erA(i,columns-2)-2.5*erA(i,columns-1))/(-1*ezField) ;
	  erOverEzC(i,j)   += (gridSizeZ/3.0)*(0.5*erC(i,columns-2)-2.5*erC(i,columns-1))/(-1*ezField) ;
	  deltaEzA(i,j)    += (gridSizeZ/3.0)*(0.5*dezA(i,columns-2)-2.5*dezA(i,columns-1))/(-1) ;
	  deltaEzC(i,j)    += (gridSizeZ/3.0)*(0.5*dezC(i,columns-2)-2.5*dezC(i,columns-1))/(-1) ;
	}
	if ( j == columns-2 ) {
	  erOverEzA(i,j)   = (gridSizeZ/3.0)*(1.5*erA(i,columns-2)+1.5*erA(i,columns-1))/(-1*ezField) ;
	  erOverEzC(i,j)   = (gridSizeZ/3.0)*(1.5*erC(i,columns-2)+1.5*erC(i,columns-1))/(-1*ezField) ;
	  deltaEzA(i,j)    = (gridSizeZ/3.0)*(1.5*dezA(i,columns-2)+1.5*dezA(i,columns-1))/(-1) ;
	  deltaEzC(i,j)    = (gridSizeZ/3.0)*(1.5*dezC(i,columns-2)+1.5*dezC(i,columns-1))/(-1) ;
	}
	if ( j == columns-1 ) {
	  erOverEzA(i,j)   = 0;   
	  erOverEzC(i,j)   = 0;
	  deltaEzA(i,j)    = 0;  
	  deltaEzC(i,j)    = 0;
	}
      }
    }

  }
  
  //  AliInfo("Step 1: Interpolation to Standard grid");

  // -------------------------------------------------------------------------------
  // Interpolate results onto the standard grid which is used for all AliTPCCorrections classes

  const Int_t order  = 1  ;  // Linear interpolation = 1, Quadratic = 2  

  Double_t  r, z;//phi, z ;
  for ( Int_t k = 0 ; k < kNPhi ; k++ ) {
    //    phi = fgkPhiList[k] ;
	
    // final lookup table
    TMatrixF &erOverEzFinal   =  *fLookUpErOverEz[k]  ;
    TMatrixF &deltaEzFinal    =  *fLookUpDeltaEz[k]   ;
	
    // calculated and integrated tables - just one phi slice
    TMatrixD &erOverEzA   =  *arrayofEroverEzA[0]  ;
    TMatrixD &deltaEzA    =  *arrayofDeltaEzA[0];
    TMatrixD &erOverEzC   =  *arrayofEroverEzC[0]  ;
    TMatrixD &deltaEzC    =  *arrayofDeltaEzC[0];    
    
    
    for ( Int_t j = 0 ; j < kNZ ; j++ ) {

      z = TMath::Abs(fgkZList[j]) ;  // z position is symmetric
    
      for ( Int_t i = 0 ; i < kNR ; i++ ) { 
	r = fgkRList[i] ;

	// Interpolate Lookup tables onto standard grid
	if (fgkZList[j]>0) {
	  erOverEzFinal(i,j)   = Interpolate2DTable(order, r, z, rows, columns, rlist1, zedlist1, erOverEzA  );
	  deltaEzFinal(i,j)    = Interpolate2DTable(order, r, z, rows, columns, rlist1, zedlist1, deltaEzA   );
	} else {
	  erOverEzFinal(i,j)   = Interpolate2DTable(order, r, z, rows, columns, rlist1, zedlist1, erOverEzC  );
	  deltaEzFinal(i,j)    = - Interpolate2DTable(order, r, z, rows, columns, rlist1, zedlist1, deltaEzC   );
	  // negative coordinate system on C side
	}

      } // end r loop
    } // end z loop
  } // end phi loop

  
  // clear the temporary arrays lists
  for ( Int_t k = 0 ; k < phiSlices ; k++ )  {

    delete arrayofErA[k];  
    delete arrayofdEzA[k];
    delete arrayofErC[k];  
    delete arrayofdEzC[k];

    delete arrayofEroverEzA[k];  
    delete arrayofDeltaEzA[k];
    delete arrayofEroverEzC[k];  
    delete arrayofDeltaEzC[k];

  }

  fZR->Close();

  // ------------------------------------------------------------------------------------------------------
  // Step 2: Load and sum up lookup table in rphi, fine grid, to emulate for example a GG leak
 
  //  AliInfo("Step 2: Preparation of the weighted look-up table");
 
  TFile *fRPhi = new TFile(fSCLookUpPOCsFileNameRPhi.Data(),"READ");
  if ( !fRPhi ) {
    AliError("Precalculated POC-looup-table in RPhi could not be found");
    return;
  } 

  // units are in [m]
  TVector *gridf2  = (TVector*) fRPhi->Get("constants"); 
  TVector &grid2 = *gridf2;
  TMatrix *coordf2  = (TMatrix*) fRPhi->Get("coordinates");
  TMatrix &coord2 = *coordf2;
  TMatrix *coordPOCf2  = (TMatrix*) fRPhi->Get("POCcoord");
  TMatrix &coordPOC2 = *coordPOCf2;

  rows      = (Int_t)grid2(0);   // number of points in r direction   
  phiSlices = (Int_t)grid2(1);   // number of points in phi           
  columns   = (Int_t)grid2(2);   // number of points in z direction   

  gridSizeR   =  (fgkOFCRadius-fgkIFCRadius)/(rows-1);   // unit in [cm]
  Float_t gridSizePhi =  TMath::TwoPi()/phiSlices;         // unit in [rad]
  gridSizeZ   =  fgkTPCZ0/(columns-1);                  // unit in [cm]
 
  // list of points as used during sum up
  Double_t  rlist2[kNRows], philist2[kNPhiSlices], zedlist2[kNColumns]; 
  for ( Int_t k = 0 ; k < phiSlices ; k++ ) {
    philist2[k] =  gridSizePhi * k;
    for ( Int_t i = 0 ; i < rows ; i++ )    {
      rlist2[i] = fgkIFCRadius + i*gridSizeR ;
      for ( Int_t j = 0 ; j < columns ; j++ ) { 
	zedlist2[j]  = j * gridSizeZ ;
      }
    }
  } // only done once
 
  // temporary matrices needed for the calculation 

  TMatrixD *arrayofErA2[kNPhiSlices], *arrayofEphiA2[kNPhiSlices], *arrayofdEzA2[kNPhiSlices];
  TMatrixD *arrayofErC2[kNPhiSlices], *arrayofEphiC2[kNPhiSlices], *arrayofdEzC2[kNPhiSlices]; 

  TMatrixD *arrayofEroverEzA2[kNPhiSlices], *arrayofEphioverEzA2[kNPhiSlices], *arrayofDeltaEzA2[kNPhiSlices]; 
  TMatrixD *arrayofEroverEzC2[kNPhiSlices], *arrayofEphioverEzC2[kNPhiSlices], *arrayofDeltaEzC2[kNPhiSlices]; 

 
  for ( Int_t k = 0 ; k < phiSlices ; k++ ) {
   
    arrayofErA2[k]   =   new TMatrixD(rows,columns) ;
    arrayofEphiA2[k] =   new TMatrixD(rows,columns) ;
    arrayofdEzA2[k]  =   new TMatrixD(rows,columns) ;
    arrayofErC2[k]   =   new TMatrixD(rows,columns) ;
    arrayofEphiC2[k] =   new TMatrixD(rows,columns) ;
    arrayofdEzC2[k]  =   new TMatrixD(rows,columns) ;

    arrayofEroverEzA2[k]   =   new TMatrixD(rows,columns) ;
    arrayofEphioverEzA2[k] =   new TMatrixD(rows,columns) ; 
    arrayofDeltaEzA2[k]    =   new TMatrixD(rows,columns) ;
    arrayofEroverEzC2[k]   =   new TMatrixD(rows,columns) ;
    arrayofEphioverEzC2[k] =   new TMatrixD(rows,columns) ; 
    arrayofDeltaEzC2[k]    =   new TMatrixD(rows,columns) ;

    // zero initialization not necessary, it is done in the constructor of TMatrix 

  }
 
    
  treePOC = (TTree*)fRPhi->Get("POCall");

  //  TVector *bEr  = 0;   // done above
  TVector *bEphi= 0;
  //  TVector *bEz  = 0;   // done above

  treePOC->SetBranchAddress("Er",&bEr);
  treePOC->SetBranchAddress("Ephi",&bEphi);
  treePOC->SetBranchAddress("Ez",&bEz);

  // Read the complete tree and do a weighted sum-up over the POC configurations
  // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  treeNumPOC = (Int_t)treePOC->GetEntries(); // Number of POC conf. in the look-up table
  ipC = 0; // POC Conf. counter (note: different to the POC number in the tree!)

  for (Int_t itreepC=0; itreepC<treeNumPOC; itreepC++) { // ------------- loop over POC configurations in tree
  
    treePOC->GetEntry(itreepC);

    // center of the POC voxel in [meter]
    Double_t r0 = coordPOC2(ipC,0);
    Double_t phi0 = coordPOC2(ipC,1);
    //    Double_t z0 = coordPOC2(ipC,2);

     // weights (charge density) at POC position on the A and C side (in C/m^3/e0)
    // note: coordinates are in [cm]
    Double_t weightA = GetSpaceChargeDensity(r0*100,phi0, 0.499, 2);  // partial load in r,phi
    Double_t weightC = GetSpaceChargeDensity(r0*100,phi0,-0.499, 2);  // partial load in r,phi
  
    //    printf("-----\n%f %f : %e %e\n",r0,phi0,weightA,weightC);

    // Summing up the vector components according to their weight

    Int_t ip = 0;
    for ( Int_t j = 0 ; j < columns ; j++ ) { 
      for ( Int_t i = 0 ; i < rows ; i++ )    {
	for ( Int_t k = 0 ; k < phiSlices ; k++ ) {
		  
	  // check wether the coordinates were screwed
	  if (TMath::Abs((coord2(0,ip)*100-rlist2[i]))>1 || 
	      TMath::Abs((coord2(1,ip)-philist2[k]))>1 || 
	      TMath::Abs((coord2(2,ip)*100-zedlist2[j]))>1) { 
	    AliError("internal error: coordinate system was screwed during the sum-up");
	    printf("lookup: (r,phi,z)=(%f,%f,%f)\n",coord2(0,ip)*100,coord2(1,ip),coord2(2,ip)*100);
	    printf("sum-up: (r,phi,z)=(%f,%f,%f)\n",rlist2[i],philist2[k],zedlist2[j]);
	    AliError("Don't trust the results of the space charge calculation!");
	  }
	  
	  // unfortunately, the lookup tables were produced to be faster for phi symmetric charges
	  // This will be the most frequent usage (hopefully)
	  // That's why we have to do this here ...

	  TMatrixD &erA   =  *arrayofErA2[k]  ;
	  TMatrixD &ephiA =  *arrayofEphiA2[k];
	  TMatrixD &dEzA  =  *arrayofdEzA2[k] ;
   
	  TMatrixD &erC   =  *arrayofErC2[k]  ;
	  TMatrixD &ephiC =  *arrayofEphiC2[k];
	  TMatrixD &dEzC  =  *arrayofdEzC2[k]   ;
   
	  // Sum up - Efield values in [V/m] -> transition to [V/cm]
	  erA(i,j) += ((*bEr)(ip)) * weightA /100;
	  erC(i,j) += ((*bEr)(ip)) * weightC /100;
	  ephiA(i,j) += ((*bEphi)(ip)) * weightA/100; // [V/rad]
	  ephiC(i,j) += ((*bEphi)(ip)) * weightC/100; // [V/rad]
	  dEzA(i,j)  += ((*bEz)(ip)) * weightA /100;
	  dEzC(i,j)  += ((*bEz)(ip)) * weightC /100;

	  // increase the counter
	  ip++;
	}
      }
    } // end coordinate loop
    
    
    // Rotation and summation in the rest of the dPhiSteps
    // which were not stored in the this tree due to storage & symmetry reasons

    
    Int_t phiPoints = (Int_t) grid2(1);
    Int_t phiPOC    = (Int_t) grid2(4);
    
    //   printf("%d %d\n",phiPOC,flagRadSym);
    
    for (Int_t phiiC = 1; phiiC<phiPOC; phiiC++) { // just used for non-radial symetric table 
      
      Double_t phi0R = phiiC*phi0*2 + phi0; // rotate further

      // weights (charge density) at POC position on the A and C side (in C/m^3/e0)
      // note: coordinates are in [cm] // ecxept z
      weightA = GetSpaceChargeDensity(r0*100,phi0R, 0.499, 2);  // partial load in r,phi
      weightC = GetSpaceChargeDensity(r0*100,phi0R,-0.499, 2);  // partial load in r,phi
    
      // printf("%f %f : %e %e\n",r0,phi0R,weightA,weightC);
         
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
	    // This will be the most frequent usage 
	    // That's why we have to do this here and not outside the loop ...
	    
	    TMatrixD &erA   =  *arrayofErA2[k]  ;
	    TMatrixD &ephiA =  *arrayofEphiA2[k];
	    TMatrixD &dEzA  =  *arrayofdEzA2[k]   ;
	    
	    TMatrixD &erC   =  *arrayofErC2[k]  ;
	    TMatrixD &ephiC =  *arrayofEphiC2[k];
	    TMatrixD &dEzC  =  *arrayofdEzC2[k]   ;
       
	    // Sum up - Efield values in [V/m] -> transition to [V/cm]
	    erA(i,j) += ((*bEr)(ipR)) * weightA /100;
	    erC(i,j) += ((*bEr)(ipR)) * weightC /100;
	    ephiA(i,j) += ((*bEphi)(ipR)) * weightA/100; // [V/rad]
	    ephiC(i,j) += ((*bEphi)(ipR)) * weightC/100; // [V/rad]
	    dEzA(i,j)  += ((*bEz)(ipR)) * weightA /100;
	    dEzC(i,j)  += ((*bEz)(ipR)) * weightC /100;

	    // increase the counter
	    ip++;
	  }
	}
      } // end coordinate loop

    } // end phi-POC summation (phiiC)

    ipC++; // POC configuration counter

    //   printf("POC: (r,phi,z) = (%f %f %f) | weight(A,C): %03.1lf %03.1lf\n",r0,phi0,z0, weightA, weightC);
    
  }




  // -------------------------------------------------------------------------------
  // Division by the Ez (drift) field and integration along z

  //  AliInfo("Step 2: Division and integration");


  for ( Int_t k = 0 ; k < phiSlices ; k++ ) { // phi loop

    // matrices holding the solution - summation of POC charges // see above
    TMatrixD &erA   =  *arrayofErA2[k]  ;
    TMatrixD &ephiA =  *arrayofEphiA2[k];
    TMatrixD &dezA  =  *arrayofdEzA2[k]   ;
    TMatrixD &erC   =  *arrayofErC2[k]  ;
    TMatrixD &ephiC =  *arrayofEphiC2[k];
    TMatrixD &dezC  =  *arrayofdEzC2[k]   ;

    // matrices which will contain the integrated fields (divided by the drift field)
    TMatrixD &erOverEzA   =  *arrayofEroverEzA2[k]  ;
    TMatrixD &ephiOverEzA =  *arrayofEphioverEzA2[k];
    TMatrixD &deltaEzA    =  *arrayofDeltaEzA2[k];
    TMatrixD &erOverEzC   =  *arrayofEroverEzC2[k]  ;
    TMatrixD &ephiOverEzC =  *arrayofEphioverEzC2[k];
    TMatrixD &deltaEzC    =  *arrayofDeltaEzC2[k];    
    
    for ( Int_t i = 0 ; i < rows ; i++ )    { // r loop
      for ( Int_t j = columns-1 ; j >= 0 ; j-- ) {// z loop 
	// Count backwards to facilitate integration over Z

	Int_t index = 1 ; // Simpsons rule if N=odd.If N!=odd then add extra point by trapezoidal rule.  

	erOverEzA(i,j) = 0; 
	ephiOverEzA(i,j) = 0;
	deltaEzA(i,j) = 0;
	erOverEzC(i,j) = 0; 
	ephiOverEzC(i,j) = 0; 
	deltaEzC(i,j) = 0;

	for ( Int_t m = j ; m < columns ; m++ ) { // integration

	  erOverEzA(i,j)   += index*(gridSizeZ/3.0)*erA(i,m)/(-1*ezField) ;
	  erOverEzC(i,j)   += index*(gridSizeZ/3.0)*erC(i,m)/(-1*ezField)  ;
	  ephiOverEzA(i,j) += index*(gridSizeZ/3.0)*ephiA(i,m)/(-1*ezField)  ;
	  ephiOverEzC(i,j) += index*(gridSizeZ/3.0)*ephiC(i,m)/(-1*ezField)  ;
	  deltaEzA(i,j)    += index*(gridSizeZ/3.0)*dezA(i,m)/(-1) ;
	  deltaEzC(i,j)    += index*(gridSizeZ/3.0)*dezC(i,m)/(-1) ;

	  if ( index != 4 )  index = 4; else index = 2 ;

	}

	if ( index == 4 ) {
	  erOverEzA(i,j)   -= (gridSizeZ/3.0)*erA(i,columns-1)/(-1*ezField) ;
	  erOverEzC(i,j)   -= (gridSizeZ/3.0)*erC(i,columns-1)/(-1*ezField) ;
	  ephiOverEzA(i,j) -= (gridSizeZ/3.0)*ephiA(i,columns-1)/(-1*ezField) ;
	  ephiOverEzC(i,j) -= (gridSizeZ/3.0)*ephiC(i,columns-1)/(-1*ezField) ;
	  deltaEzA(i,j)    -= (gridSizeZ/3.0)*dezA(i,columns-1)/(-1) ;
	  deltaEzC(i,j)    -= (gridSizeZ/3.0)*dezC(i,columns-1)/(-1) ;
	}
	if ( index == 2 ) {
	  erOverEzA(i,j)   += (gridSizeZ/3.0)*(0.5*erA(i,columns-2)-2.5*erA(i,columns-1))/(-1*ezField) ;
	  erOverEzC(i,j)   += (gridSizeZ/3.0)*(0.5*erC(i,columns-2)-2.5*erC(i,columns-1))/(-1*ezField) ;
	  ephiOverEzA(i,j) += (gridSizeZ/3.0)*(0.5*ephiA(i,columns-2)-2.5*ephiA(i,columns-1))/(-1*ezField) ;
	  ephiOverEzC(i,j) += (gridSizeZ/3.0)*(0.5*ephiC(i,columns-2)-2.5*ephiC(i,columns-1))/(-1*ezField) ;
	  deltaEzA(i,j)    += (gridSizeZ/3.0)*(0.5*dezA(i,columns-2)-2.5*dezA(i,columns-1))/(-1) ;
	  deltaEzC(i,j)    += (gridSizeZ/3.0)*(0.5*dezC(i,columns-2)-2.5*dezC(i,columns-1))/(-1) ;
	}
	if ( j == columns-2 ) {
	  erOverEzA(i,j)   = (gridSizeZ/3.0)*(1.5*erA(i,columns-2)+1.5*erA(i,columns-1))/(-1*ezField) ;
	  erOverEzC(i,j)   = (gridSizeZ/3.0)*(1.5*erC(i,columns-2)+1.5*erC(i,columns-1))/(-1*ezField) ;
	  ephiOverEzA(i,j) = (gridSizeZ/3.0)*(1.5*ephiA(i,columns-2)+1.5*ephiA(i,columns-1))/(-1*ezField) ;
	  ephiOverEzC(i,j) = (gridSizeZ/3.0)*(1.5*ephiC(i,columns-2)+1.5*ephiC(i,columns-1))/(-1*ezField) ;
	  deltaEzA(i,j)    = (gridSizeZ/3.0)*(1.5*dezA(i,columns-2)+1.5*dezA(i,columns-1))/(-1) ;
	  deltaEzC(i,j)    = (gridSizeZ/3.0)*(1.5*dezC(i,columns-2)+1.5*dezC(i,columns-1))/(-1) ;
	}
	if ( j == columns-1 ) {
	  erOverEzA(i,j)   = 0;   
	  erOverEzC(i,j)   = 0;
	  ephiOverEzA(i,j) = 0; 
	  ephiOverEzC(i,j) = 0;
	  deltaEzA(i,j)    = 0;  
	  deltaEzC(i,j)    = 0;
	}
      }
    }

  }
  
  AliInfo("Step 2: Interpolation to Standard grid");

  // -------------------------------------------------------------------------------
  // Interpolate results onto the standard grid which is used for all AliTPCCorrections classes


  for ( Int_t k = 0 ; k < kNPhi ; k++ ) {
    Double_t phi = fgkPhiList[k] ;
	
    // final lookup table
    TMatrixF &erOverEzFinal   =  *fLookUpErOverEz[k]  ;
    TMatrixF &ephiOverEzFinal =  *fLookUpEphiOverEz[k];
    TMatrixF &deltaEzFinal    =  *fLookUpDeltaEz[k]   ;
	
    for ( Int_t j = 0 ; j < kNZ ; j++ ) {

      z = TMath::Abs(fgkZList[j]) ;  // z position is symmetric
    
      for ( Int_t i = 0 ; i < kNR ; i++ ) { 
	r = fgkRList[i] ;

	// Interpolate Lookup tables onto standard grid
	if (fgkZList[j]>0) {
	  erOverEzFinal(i,j)   += Interpolate3DTable(order, r, z, phi, rows, columns, phiSlices, 
					       rlist2, zedlist2, philist2, arrayofEroverEzA2  );
	  ephiOverEzFinal(i,j) += Interpolate3DTable(order, r, z, phi, rows, columns, phiSlices,
					       rlist2, zedlist2, philist2, arrayofEphioverEzA2);
	  deltaEzFinal(i,j)    += Interpolate3DTable(order, r, z, phi, rows, columns, phiSlices,
					       rlist2, zedlist2, philist2, arrayofDeltaEzA2   );
	} else {
	  erOverEzFinal(i,j)   += Interpolate3DTable(order, r, z, phi, rows, columns, phiSlices, 
					       rlist2, zedlist2, philist2, arrayofEroverEzC2  );
	  ephiOverEzFinal(i,j) += Interpolate3DTable(order, r, z, phi, rows, columns, phiSlices,
					       rlist2, zedlist2, philist2, arrayofEphioverEzC2);
	  deltaEzFinal(i,j)  -=  Interpolate3DTable(order, r, z, phi, rows, columns, phiSlices,
					       rlist2, zedlist2, philist2, arrayofDeltaEzC2   );
	}

      } // end r loop
    } // end z loop
  } // end phi loop
  
  
  // clear the temporary arrays lists
  for ( Int_t k = 0 ; k < phiSlices ; k++ )  {

    delete arrayofErA2[k];  
    delete arrayofEphiA2[k];
    delete arrayofdEzA2[k];
    delete arrayofErC2[k];  
    delete arrayofEphiC2[k];
    delete arrayofdEzC2[k];

    delete arrayofEroverEzA2[k];  
    delete arrayofEphioverEzA2[k];
    delete arrayofDeltaEzA2[k];
    delete arrayofEroverEzC2[k];  
    delete arrayofEphioverEzC2[k];
    delete arrayofDeltaEzC2[k];

  }
 
  fRPhi->Close();
  
  // FINISHED

  fInitLookUp = kTRUE;

}

void AliTPCSpaceCharge3D::InitSpaceCharge3DDistortionCourse() {
  //
  // Initialization of the Lookup table which contains the solutions of the 
  // "space charge" (poisson) problem 
  //
  // The sum-up uses a look-up table which contains different discretized Space charge fields 
  // in order to calculate the corresponding field deviations due to a given (discretized)
  // space charge distribution ....
  //
  // Method of calculation: Weighted sum-up of the different fields within the look up table
  // Note: Full 3d version: Course and slow ...

  if (fInitLookUp) {
    AliInfo("Lookup table was already initialized!");
    //    return;
  }

  AliInfo("Preparation of the weighted look-up table");
   
  TFile *f = new TFile(fSCLookUpPOCsFileName3D.Data(),"READ");
  if ( !f ) {
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
    //   AliInfo("LOOK UP TABLE IS RADIAL SYMETTRIC - Field in Phi is ZERO");
    flagRadSym=1;
  }

  Int_t rows      = (Int_t)grid(0);   // number of points in r direction 
  Int_t phiSlices = (Int_t)grid(1);   // number of points in phi         
  Int_t columns   = (Int_t)grid(2);   // number of points in z direction 

  const Float_t gridSizeR   =  (fgkOFCRadius-fgkIFCRadius)/(rows-1);   // unit in [cm]
  const Float_t gridSizePhi =  TMath::TwoPi()/phiSlices;         // unit in [rad]
  const Float_t gridSizeZ   =  fgkTPCZ0/(columns-1);                  // unit in [cm]
 
  // temporary matrices needed for the calculation
  TMatrixD *arrayofErA[kNPhiSlices], *arrayofEphiA[kNPhiSlices], *arrayofdEzA[kNPhiSlices];
  TMatrixD *arrayofErC[kNPhiSlices], *arrayofEphiC[kNPhiSlices], *arrayofdEzC[kNPhiSlices];

  TMatrixD *arrayofEroverEzA[kNPhiSlices], *arrayofEphioverEzA[kNPhiSlices], *arrayofDeltaEzA[kNPhiSlices];
  TMatrixD *arrayofEroverEzC[kNPhiSlices], *arrayofEphioverEzC[kNPhiSlices], *arrayofDeltaEzC[kNPhiSlices];

 
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
  Double_t  rlist[kNRows], zedlist[kNColumns] , philist[kNPhiSlices];
  for ( Int_t k = 0 ; k < phiSlices ; k++ ) {
    philist[k] =  gridSizePhi * k;
    for ( Int_t i = 0 ; i < rows ; i++ )    {
      rlist[i] = fgkIFCRadius + i*gridSizeR ;
      for ( Int_t j = 0 ; j < columns ; j++ ) { 
	zedlist[j]  = j * gridSizeZ ;
      }
    }
  } // only done once
  
  
  TTree *treePOC = (TTree*)f->Get("POCall");

  TVector *bEr  = 0;   TVector *bEphi= 0;   TVector *bEz  = 0;
  
  treePOC->SetBranchAddress("Er",&bEr);
  if (!flagRadSym) treePOC->SetBranchAddress("Ephi",&bEphi);
  treePOC->SetBranchAddress("Ez",&bEz);


  // Read the complete tree and do a weighted sum-up over the POC configurations
  // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  Int_t treeNumPOC = (Int_t)treePOC->GetEntries(); // Number of POC conf. in the look-up table
  Int_t ipC = 0; // POC Conf. counter (note: different to the POC number in the tree!)

  for (Int_t itreepC=0; itreepC<treeNumPOC; itreepC++) { // ------------- loop over POC configurations in tree
  
    treePOC->GetEntry(itreepC);

    // center of the POC voxel in [meter]
    Double_t r0 = coordPOC(ipC,0);
    Double_t phi0 = coordPOC(ipC,1);
    Double_t z0 = coordPOC(ipC,2);

    ipC++; // POC configuration counter

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
	  if (TMath::Abs((coord(0,ip)*100-rlist[i]))>1 || 
	      TMath::Abs((coord(1,ip)-philist[k]))>1 || 
	      TMath::Abs((coord(2,ip)*100-zedlist[j]))>1) { 
	    AliError("internal error: coordinate system was screwed during the sum-up");
	    printf("lookup: (r,phi,z)=(%f,%f,%f)\n",coord(0,ip)*100,coord(1,ip),coord(2,ip)*100);
	    printf("sum-up: (r,phi,z)=(%f,%f,%f)\n",rlist[i],philist[k],zedlist[j]);
	    AliError("Don't trust the results of the space charge calculation!");
	  }
	  
	  // unfortunately, the lookup tables were produced to be faster for phi symmetric charges
	  // This will be the most frequent usage (hopefully)
	  // That's why we have to do this here ...

	  TMatrixD &erA   =  *arrayofErA[k]  ;
	  TMatrixD &ephiA =  *arrayofEphiA[k];
	  TMatrixD &dEzA  =  *arrayofdEzA[k]   ;
   
	  TMatrixD &erC   =  *arrayofErC[k]  ;
	  TMatrixD &ephiC =  *arrayofEphiC[k];
	  TMatrixD &dEzC  =  *arrayofdEzC[k]   ;
   
	  // Sum up - Efield values in [V/m] -> transition to [V/cm]
	  erA(i,j) += ((*bEr)(ip)) * weightA /100;
	  erC(i,j) += ((*bEr)(ip)) * weightC /100;
	  if (!flagRadSym) {
	    ephiA(i,j) += ((*bEphi)(ip)) * weightA/100; // [V/rad]
	    ephiC(i,j) += ((*bEphi)(ip)) * weightC/100; // [V/rad]
	  }
	  dEzA(i,j)  += ((*bEz)(ip)) * weightA /100;
	  dEzC(i,j)  += ((*bEz)(ip)) * weightC /100;

	  // increase the counter
	  ip++;
	}
      }
    } // end coordinate loop
    
    
    // Rotation and summation in the rest of the dPhiSteps
    // which were not stored in the this tree due to storage & symmetry reasons

    Int_t phiPoints = (Int_t) grid(1);
    Int_t phiPOC    = (Int_t) grid(4);
    
    //   printf("%d %d\n",phiPOC,flagRadSym);
    
    for (Int_t phiiC = 1; phiiC<phiPOC; phiiC++) { // just used for non-radial symetric table 
      
      r0 = coordPOC(ipC,0);
      phi0 = coordPOC(ipC,1);
      z0 = coordPOC(ipC,2);
      
      ipC++; // POC conf. counter
      
      // weights (charge density) at POC position on the A and C side (in C/m^3/e0)
      // note: coordinates are in [cm]
      weightA = GetSpaceChargeDensity(r0*100,phi0, z0*100); 
      weightC = GetSpaceChargeDensity(r0*100,phi0,-z0*100);
      
      //     printf("%f %f %f: %e %e\n",r0,phi0,z0,weightA,weightC);
       
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
	    // This will be the most frequent usage 
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
	      ephiA(i,j) += ((*bEphi)(ipR)) * weightA/100; // [V/rad]
	      ephiC(i,j) += ((*bEphi)(ipR)) * weightC/100; // [V/rad]
	    }
	    dEzA(i,j)  += ((*bEz)(ipR)) * weightA /100;
	    dEzC(i,j)  += ((*bEz)(ipR)) * weightC /100;

	    // increase the counter
	    ip++;
	  }
	}
      } // end coordinate loop

    } // end phi-POC summation (phiiC)
   

    // printf("POC: (r,phi,z) = (%f %f %f) | weight(A,C): %03.1lf %03.1lf\n",r0,phi0,z0, weightA, weightC);
    
  }



  // -------------------------------------------------------------------------------
  // Division by the Ez (drift) field and integration along z

  AliInfo("Division and integration");

  Double_t ezField = (fgkCathodeV-fgkGG)/fgkTPCZ0; // = Electric Field (V/cm) Magnitude ~ -400 V/cm;

  for ( Int_t k = 0 ; k < phiSlices ; k++ ) { // phi loop

    // matrices holding the solution - summation of POC charges // see above
    TMatrixD &erA   =  *arrayofErA[k]  ;
    TMatrixD &ephiA =  *arrayofEphiA[k];
    TMatrixD &dezA  =  *arrayofdEzA[k]   ;
    TMatrixD &erC   =  *arrayofErC[k]  ;
    TMatrixD &ephiC =  *arrayofEphiC[k];
    TMatrixD &dezC  =  *arrayofdEzC[k]   ;

    // matrices which will contain the integrated fields (divided by the drift field)
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
	  deltaEzA(i,j)    += index*(gridSizeZ/3.0)*dezA(i,m)/(-1) ;
	  deltaEzC(i,j)    += index*(gridSizeZ/3.0)*dezC(i,m)/(-1) ;

	  if ( index != 4 )  index = 4; else index = 2 ;

	}

	if ( index == 4 ) {
	  erOverEzA(i,j)   -= (gridSizeZ/3.0)*erA(i,columns-1)/(-1*ezField) ;
	  erOverEzC(i,j)   -= (gridSizeZ/3.0)*erC(i,columns-1)/(-1*ezField) ;
	  if (!flagRadSym) {
	    ephiOverEzA(i,j) -= (gridSizeZ/3.0)*ephiA(i,columns-1)/(-1*ezField) ;
	    ephiOverEzC(i,j) -= (gridSizeZ/3.0)*ephiC(i,columns-1)/(-1*ezField) ;
	  }
	  deltaEzA(i,j)    -= (gridSizeZ/3.0)*dezA(i,columns-1)/(-1) ;
	  deltaEzC(i,j)    -= (gridSizeZ/3.0)*dezC(i,columns-1)/(-1) ;
	}
	if ( index == 2 ) {
	  erOverEzA(i,j)   += (gridSizeZ/3.0)*(0.5*erA(i,columns-2)-2.5*erA(i,columns-1))/(-1*ezField) ;
	  erOverEzC(i,j)   += (gridSizeZ/3.0)*(0.5*erC(i,columns-2)-2.5*erC(i,columns-1))/(-1*ezField) ;
	  if (!flagRadSym) {
	    ephiOverEzA(i,j) += (gridSizeZ/3.0)*(0.5*ephiA(i,columns-2)-2.5*ephiA(i,columns-1))/(-1*ezField) ;
	    ephiOverEzC(i,j) += (gridSizeZ/3.0)*(0.5*ephiC(i,columns-2)-2.5*ephiC(i,columns-1))/(-1*ezField) ;
	  }
	  deltaEzA(i,j)    += (gridSizeZ/3.0)*(0.5*dezA(i,columns-2)-2.5*dezA(i,columns-1))/(-1) ;
	  deltaEzC(i,j)    += (gridSizeZ/3.0)*(0.5*dezC(i,columns-2)-2.5*dezC(i,columns-1))/(-1) ;
	}
	if ( j == columns-2 ) {
	  erOverEzA(i,j)   = (gridSizeZ/3.0)*(1.5*erA(i,columns-2)+1.5*erA(i,columns-1))/(-1*ezField) ;
	  erOverEzC(i,j)   = (gridSizeZ/3.0)*(1.5*erC(i,columns-2)+1.5*erC(i,columns-1))/(-1*ezField) ;
	  if (!flagRadSym) {
	    ephiOverEzA(i,j) = (gridSizeZ/3.0)*(1.5*ephiA(i,columns-2)+1.5*ephiA(i,columns-1))/(-1*ezField) ;
	    ephiOverEzC(i,j) = (gridSizeZ/3.0)*(1.5*ephiC(i,columns-2)+1.5*ephiC(i,columns-1))/(-1*ezField) ;
	  }
	  deltaEzA(i,j)    = (gridSizeZ/3.0)*(1.5*dezA(i,columns-2)+1.5*dezA(i,columns-1))/(-1) ;
	  deltaEzC(i,j)    = (gridSizeZ/3.0)*(1.5*dezC(i,columns-2)+1.5*dezC(i,columns-1))/(-1) ;
	}
	if ( j == columns-1 ) {
	  erOverEzA(i,j)   = 0;   
	  erOverEzC(i,j)   = 0;
	  if (!flagRadSym) {
	    ephiOverEzA(i,j) = 0; 
	    ephiOverEzC(i,j) = 0;
	  }
	  deltaEzA(i,j)    = 0;  
	  deltaEzC(i,j)    = 0;
	}
      }
    }

  }
  

 
  AliInfo("Interpolation to Standard grid");

  // -------------------------------------------------------------------------------
  // Interpolate results onto the standard grid which is used for all AliTPCCorrections classes

  const Int_t order  = 1  ;  // Linear interpolation = 1, Quadratic = 2  

  Double_t  r, phi, z ;
  for ( Int_t k = 0 ; k < kNPhi ; k++ ) {
    phi = fgkPhiList[k] ;
	
    TMatrixF &erOverEz   =  *fLookUpErOverEz[k]  ;
    TMatrixF &ephiOverEz =  *fLookUpEphiOverEz[k];
    TMatrixF &deltaEz    =  *fLookUpDeltaEz[k]   ;
	
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


void AliTPCSpaceCharge3D::SetSCDataFileName(TString fname) {
  //
  // Set & load the Space charge density distribution from a file 
  // (linear interpolation onto a standard grid)
  //

 
  fSCDataFileName = fname;

  TFile *f = new TFile(fSCDataFileName.Data(),"READ");
  if (!f) { 
    AliError(Form("File %s, which should contain the space charge distribution, could not be found",
		  fSCDataFileName.Data()));
    return;
  }
 
  TH2F *densityRZ = (TH2F*) f->Get("SpaceChargeInRZ");
  if (!densityRZ) { 
    AliError(Form("The indicated file (%s) does not contain a histogram called %s",
		  fSCDataFileName.Data(),"SpaceChargeInRZ"));
    return;
  }

  TH3F *densityRPhi = (TH3F*) f->Get("SpaceChargeInRPhi");
  if (!densityRPhi) { 
    AliError(Form("The indicated file (%s) does not contain a histogram called %s",
		  fSCDataFileName.Data(),"SpaceChargeInRPhi"));
    return;
  }
 

  Double_t  r, phi, z ;

  TMatrixD &scDensityInRZ   =  *fSCdensityInRZ;
  TMatrixD &scDensityInRPhiA   =  *fSCdensityInRPhiA;
  TMatrixD &scDensityInRPhiC   =  *fSCdensityInRPhiC;
  for ( Int_t k = 0 ; k < kNPhi ; k++ ) {
    phi = fgkPhiList[k] ;
    TMatrixF &scDensity   =  *fSCdensityDistribution[k]  ;
    for ( Int_t j = 0 ; j < kNZ ; j++ ) {
      z = fgkZList[j] ; 
      for ( Int_t i = 0 ; i < kNR ; i++ ) { 
	r = fgkRList[i] ;

	// partial load in (r,z)
	if (k==0) // do just once
	  scDensityInRZ(i,j) =  densityRZ->Interpolate(r,z);

	// partial load in (r,phi)
	if ( j==0 || j == kNZ/2 ) {
	  if (z>0) 
	    scDensityInRPhiA(i,k) =  densityRPhi->Interpolate(r,phi,0.499);  // A side
	  else 
	    scDensityInRPhiC(i,k) =  densityRPhi->Interpolate(r,phi,-0.499); // C side
	}

	// Full 3D configuration ...
	if (z>0) 
	   scDensity(i,j) = scDensityInRZ(i,j) + scDensityInRPhiA(i,k); 
	else
	   scDensity(i,j) = scDensityInRZ(i,j) + scDensityInRPhiC(i,k); 
      }
    }
  }

  f->Close();

  fInitLookUp = kFALSE;

  
}


Float_t  AliTPCSpaceCharge3D::GetSpaceChargeDensity(Float_t r, Float_t phi, Float_t z, Int_t mode) {
  //
  // returns the (input) space charge density at a given point according 
  // Note: input in [cm], output in [C/m^3/e0] !!
  //

  while (phi<0) phi += TMath::TwoPi();
  while (phi>TMath::TwoPi()) phi -= TMath::TwoPi();


  // Float_t sc =fSCdensityDistribution->Interpolate(r0,phi0,z0);
  Int_t order = 1; //
  Float_t sc = 0;

  if (mode == 0) { // return full load
    sc = Interpolate3DTable(order, r, z, phi, kNR, kNZ, kNPhi, 
			    fgkRList, fgkZList, fgkPhiList, fSCdensityDistribution );
 
  } else if (mode == 1) { // return partial load in (r,z)
    TMatrixD &scDensityInRZ   =  *fSCdensityInRZ;
    sc = Interpolate2DTable(order, r, z, kNR, kNZ, fgkRList, fgkZList, scDensityInRZ );
    
  } else if (mode == 2) { // return partial load in (r,phi)

    if (z>0) {
      TMatrixD &scDensityInRPhi   =  *fSCdensityInRPhiA;
      sc = Interpolate2DTable(order, r, phi, kNR, kNPhi, fgkRList, fgkPhiList, scDensityInRPhi );
    } else {
      TMatrixD &scDensityInRPhi   =  *fSCdensityInRPhiC;
      sc = Interpolate2DTable(order, r, phi, kNR, kNPhi, fgkRList, fgkPhiList, scDensityInRPhi );
    }

  } else {
    // should i give a warning?
    sc = 0;
  }
  
  //  printf("%f %f %f: %f\n",r,phi,z,sc);
  
  return sc;
}


TH2F * AliTPCSpaceCharge3D::CreateHistoSCinXY(Float_t z, Int_t nx, Int_t ny, Int_t mode) {
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
	Float_t sc = GetSpaceChargeDensity(r,phi,z,mode)/fgke0; // in [C/m^3/e0]
	h->SetBinContent(ix,iy,sc); 
      } else {
	h->SetBinContent(ix,iy,0.);
      }
    }
  }
  
  return h;
} 

TH2F * AliTPCSpaceCharge3D::CreateHistoSCinZR(Float_t phi, Int_t nz, Int_t nr,Int_t mode ) {
  //
  // return a simple histogramm containing the space charge distribution (input for the calculation)
  //

  TH2F *h=CreateTH2F("spaceCharge",GetTitle(),"z [cm]","r [cm]","#rho_{sc} [C/m^{3}/e_{0}]",
		     nz,-250.,250.,nr,85.,250.);

  for (Int_t ir=1;ir<=nr;++ir) {
    Float_t r = h->GetYaxis()->GetBinCenter(ir);
    for (Int_t iz=1;iz<=nz;++iz) {
      Float_t z = h->GetXaxis()->GetBinCenter(iz);
      Float_t sc = GetSpaceChargeDensity(r,phi,z,mode)/fgke0; // in [C/m^3/e0]
      h->SetBinContent(iz,ir,sc);
    }
  }

  return h;
} 

void AliTPCSpaceCharge3D::WriteChargeDistributionToFile(const char* fname) {
  //
  // Example on how to write a Space charge distribution into a File
  //  (see below: estimate from scaling STAR measurements to Alice)
  // Charge distribution is splitted into two (RZ and RPHI) in order to speed up
  // the needed calculation time
  //

  TFile *f = new TFile(fname,"RECREATE");
  
  // some grid, not too course
  Int_t nr = 350;
  Int_t nphi = 180;
  Int_t nz = 500;

  Double_t dr = (fgkOFCRadius-fgkIFCRadius)/(nr+1);
  Double_t dphi = TMath::TwoPi()/(nphi+1);
  Double_t dz = 500./(nz+1);
  Double_t safty = 0.; // due to a root bug which does not interpolate the boundary (first and last bin) correctly


  // Charge distribution in ZR (rotational symmetric) ------------------

  TH2F *histoZR = new TH2F("chargeZR","chargeZR",
			   nr,fgkIFCRadius-dr-safty,fgkOFCRadius+dr+safty,
			   nz,-250-dz-safty,250+dz+safty);
 
  for (Int_t ir=1;ir<=nr;++ir) {
    Double_t rp = histoZR->GetXaxis()->GetBinCenter(ir);
    for (Int_t iz=1;iz<=nz;++iz) {
      Double_t zp = histoZR->GetYaxis()->GetBinCenter(iz);
      
      // recalculation to meter
      Double_t lZ = 2.5; // approx. TPC drift length
      Double_t rpM = rp/100.; // in [m]
      Double_t zpM = TMath::Abs(zp/100.); // in [m]
      
      // setting of mb multiplicity and Interaction rate
      Double_t multiplicity = 950;
      Double_t intRate = 7800;

      // calculation of "scaled" parameters
      Double_t a = multiplicity*intRate/79175;
      Double_t b = a/lZ;
	
      Double_t charge = (a - b*zpM)/(rpM*rpM); // charge in [C/m^3/e0]
	
      charge = charge*fgke0; // [C/m^3]
	
      if (zp<0) charge *= 0.9; // e.g. slightly less on C side due to front absorber

      //  charge = 0; // for tests
      histoZR->SetBinContent(ir,iz,charge); 
    }
  }
    
  histoZR->Write("SpaceChargeInRZ");
  

  // Charge distribution in RPhi (e.g. Floating GG wire) ------------
  
  TH3F *histoRPhi = new TH3F("chargeRPhi","chargeRPhi",
			     nr,fgkIFCRadius-dr-safty,fgkOFCRadius+dr+safty,
			     nphi,0-dphi-safty,TMath::TwoPi()+dphi+safty,
			     2,-1,1); // z part - to allow A and C side differences
  
  // some 'arbitrary' GG leaks
  Int_t   nGGleaks = 5;
  Double_t secPosA[5]    = {3,6,6,11,13};         // sector
  Double_t radialPosA[5] = {125,100,160,200,230}; // radius in cm
  Double_t secPosC[5]    = {1,8,12,15,15};        // sector
  Double_t radialPosC[5] = {245,120,140,120,190}; // radius in cm

  for (Int_t ir=1;ir<=nr;++ir) {
    Double_t rp = histoRPhi->GetXaxis()->GetBinCenter(ir);
    for (Int_t iphi=1;iphi<=nphi;++iphi) {
      Double_t phip = histoRPhi->GetYaxis()->GetBinCenter(iphi);
      for (Int_t iz=1;iz<=2;++iz) {
	Double_t zp = histoRPhi->GetZaxis()->GetBinCenter(iz);
	
	Double_t charge = 0;
	
	for (Int_t igg = 0; igg<nGGleaks; igg++) { // loop over GG leaks
	  
	  // A side
	  Double_t secPos = secPosA[igg]; 
	  Double_t radialPos = radialPosA[igg];

	  if (zp<0) { // C side
	    secPos = secPosC[igg]; 
	    radialPos = radialPosC[igg];
	  }	    

	  // some 'arbitrary' GG leaks
	  if (  (phip<(TMath::Pi()/9*(secPos+1)) && phip>(TMath::Pi()/9*secPos) ) ) { // sector slice
	    if ( rp>(radialPos-2.5) && rp<(radialPos+2.5))  // 5 cm slice
	      charge = 300; 
	  }
	  
	}		
	
	charge = charge*fgke0; // [C/m^3]

	histoRPhi->SetBinContent(ir,iphi,iz,charge); 
      }
    }
  }

  histoRPhi->Write("SpaceChargeInRPhi");

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
