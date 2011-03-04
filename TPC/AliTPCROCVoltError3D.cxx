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
//   <h2> AliTPCROCVoltError3D class   </h2>    
//   The class calculates the space point distortions due to z offsets of the 
//   ROCs via the residual voltage technique (Poisson relaxation) in 3D. 
//   Since the GG (part of the ROCs) represents the closure of the FC in z direction,
//   every misalignment in z produces not only dz distortions but also electrical 
//   field inhomogeneities throughout the volume, which produces additional dr and rd$\phi$ distortions.
//   <p>
//   Each ROC can be misaligned (in z direction) in three ways. A general z0 offset, 
//   an inclination along the x and an inclination along the y axis. The z-misalignment's
//   can be set via the function SetROCData(TMatrixD *mat) for each single chamber (ROC). 
//   The array size has to be (72,3) representing the 72 chambers according to the 
//   offline numbering scheme (IROC: roc$<$36; OROC: roc$\geq$36) and the three misalignment's
//   (see the source code for further details).
//   <p>
//   Internally, these z offsets (unit is cm)  are recalculated into residual voltage 
//   equivalents in order to make use of the relaxation technique. 
//   <p>
//   One has two possibilities when calculating the $dz$ distortions. The resulting 
//   distortions are purely due to the change of the drift velocity (along with the 
//   change of the drift field) when the SetROCDisplacement is FALSE. <br>
//   For this class, this is a rather unphysical setting and should be avoided. <br>
//   When the flag is set to TRUE, the corresponding offset in z is added to the dz 
//   calculation of the outer ROCs. <br>
//   For the Alice TPC gas, both effects are of similar magnitude. This means, if the 
//   drift length is sufficiently large, a z-offset of a chamber appears to have (approx.) 
//   twice the magnitude when one looks at the actual dz distortions.
//   <p>
//   In addition, this class allows a correction regarding the different arrival times 
//   of the electrons due to the geometrical difference of the inner and outer chambers.
//   The impact was simulated via Garfield. If the flag is set via the 
//   function SetElectronArrivalCorrection, the electron-arrival correction is added to the dz calculation.
// End_Html
//
// Begin_Macro(source) 
//   {
//   gROOT->SetStyle("Plain"); gStyle->SetPalette(1);
//   TCanvas *c2 = new TCanvas("c2","c2",500,400); 
//   AliTPCROCVoltError3D roc; 
//   roc.SetElectronArrivalCorrection(kFALSE);  // Correction for electron arrival offset, IROC vs OROC
//   roc.SetROCDisplacement(kTRUE);   // include the chamber offset in z when calculating the dz 
//   roc.SetOmegaTauT1T2(0,1,1); // B=0
//   roc.CreateHistoDZinXY(1.,300,300)->Draw("colz"); 
//   return c2;
//   } 
// End_Macro
//
// Begin_Html
//   <p>
//   Date: 08/08/2010    <br>                                                 
//   Authors: Jim Thomas, Stefan Rossegger                                
// End_Html 
// _________________________________________________________________


#include "AliMagF.h"
#include "TGeoGlobalMagField.h"
#include "AliTPCcalibDB.h"
#include "AliTPCParam.h"
#include "AliLog.h"
#include "TMatrixD.h"
#include "TFile.h"

#include "TMath.h"
#include "AliTPCROC.h"
#include "AliTPCROCVoltError3D.h"

ClassImp(AliTPCROCVoltError3D)

AliTPCROCVoltError3D::AliTPCROCVoltError3D()
  : AliTPCCorrection("ROCVoltErrors","ROC z alignment Errors"),
    fC0(0.),fC1(0.),
    fROCdisplacement(kTRUE),
    fElectronArrivalCorrection(kTRUE),
    fInitLookUp(kFALSE),
    fROCDataFileName(""),  
    fdzDataLinFit(0)
{
  //
  // default constructor
  //

  // Array which will contain the solution according to the setted boundary conditions
  // main input: z alignment of the Read Out chambers
  // see InitROCVoltError3D() function
  for ( Int_t k = 0 ; k < kNPhi ; k++ ) {
    fLookUpErOverEz[k]   =  new TMatrixF(kNR,kNZ);  
    fLookUpEphiOverEz[k] =  new TMatrixF(kNR,kNZ);
    fLookUpDeltaEz[k]    =  new TMatrixF(kNR,kNZ);   
  }
  fROCDataFileName="$ALICE_ROOT/TPC/Calib/maps/TPCROCdzSurvey.root";
  SetROCDataFileName(fROCDataFileName.Data()); // initialization of fdzDataLinFit is included

}

AliTPCROCVoltError3D::~AliTPCROCVoltError3D() {
  //
  // destructor
  //
  
  for ( Int_t k = 0 ; k < kNPhi ; k++ ) {
    delete fLookUpErOverEz[k];
    delete fLookUpEphiOverEz[k];
    delete fLookUpDeltaEz[k];
  }

  delete fdzDataLinFit;
}

void AliTPCROCVoltError3D::SetROCData(TMatrixD * matrix){
  //
  // Set a z alignment map of the chambers not via a file, but
  // directly via a TMatrix(72,3), where dz = p0 + p1*(lx-133.4) + p2*ly (all in cm)
  //
  if (!fdzDataLinFit) fdzDataLinFit=new TMatrixD(*matrix);
  else *fdzDataLinFit = *matrix;
}


void AliTPCROCVoltError3D::Init() {
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

  if (!fInitLookUp) InitROCVoltError3D();
}

void AliTPCROCVoltError3D::Update(const TTimeStamp &/*timeStamp*/) {
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

void  AliTPCROCVoltError3D::SetROCDataFileName(const char * fname) {
  //
  // Set / load the ROC data (linear fit of ROC misalignments)
  //

  fROCDataFileName = fname;
  
  TFile f(fROCDataFileName.Data(),"READ");
  TMatrixD *m = (TMatrixD*) f.Get("dzSurveyLinFitData");
  TMatrixD &mf = *m;

  // prepare some space

  if (fdzDataLinFit) delete fdzDataLinFit;
  fdzDataLinFit = new TMatrixD(72,3);
  TMatrixD &dataIntern = *fdzDataLinFit;
  
  for (Int_t iroc=0;iroc<72;iroc++) {
    dataIntern(iroc,0) = mf(iroc,0);  // z0 offset
    dataIntern(iroc,1) = mf(iroc,1);  // slope in x
    dataIntern(iroc,2) = mf(iroc,2);  // slope in y
  }

  f.Close();

  fInitLookUp = kFALSE;

}

void AliTPCROCVoltError3D::GetCorrection(const Float_t x[],const Short_t roc,Float_t dx[]) {
  //
  // Calculates the correction due e.g. residual voltage errors on the TPC boundaries
  //   
  const Double_t kEpsilon=Double_t(FLT_MIN);
  if (!fInitLookUp) {
    AliInfo("Lookup table was not initialized! Perform the inizialisation now ...");
    InitROCVoltError3D();
  }
  static Bool_t forceInit=kTRUE; //temporary needed for back compatibility with old OCDB entries
  if (forceInit&&fLookUpErOverEz[0]){
    if (fLookUpErOverEz[0]->Sum()<kEpsilon){//temporary needed for back compatibility with old OCDB entries
      ForceInitROCVoltError3D();
    }
    forceInit=kFALSE;
  }

  
  Int_t   order     = 1 ;    // FIXME: hardcoded? Linear interpolation = 1, Quadratic = 2         

  Float_t intEr, intEphi, intDeltaEz;
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
  intDeltaEz = Interpolate3DTable(order, r, z, phi, kNR, kNZ, kNPhi, 
				  fgkRList, fgkZList, fgkPhiList, fLookUpDeltaEz   );

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


  if (fElectronArrivalCorrection) {

    // correction for the OROC (in average, a 0.014usec longer drift time
    // due to different position of the anode wires) -> vd*dt -> 2.64*0.014 = 0.0369 cm
    // FIXME: correction are token from Magboltz simulations
    //        should be loaded from a database
 
    AliTPCROC * rocInfo = AliTPCROC::Instance();
    Double_t rCrossingROC  =  (rocInfo->GetPadRowRadii(0,62)+rocInfo->GetPadRowRadii(36,0))/2;
  
    if (r>rCrossingROC) {
      if (sign==1)
	dx[2] += 0.0369; // A side - negative correction
      else
	dx[2] -= 0.0369; // C side - positive correction
    }
    
  }

}

void AliTPCROCVoltError3D::InitROCVoltError3D() {
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
  TMatrixD *arrayofEroverEz[kPhiSlices], *arrayofEphioverEz[kPhiSlices], *arrayofDeltaEz[kPhiSlices] ; 

  for ( Int_t k = 0 ; k < kPhiSlices ; k++ ) {
    arrayofArrayV[k]     =   new TMatrixD(kRows,kColumns) ;
    arrayofCharge[k]     =   new TMatrixD(kRows,kColumns) ;
    arrayofEroverEz[k]   =   new TMatrixD(kRows,kColumns) ;
    arrayofEphioverEz[k] =   new TMatrixD(kRows,kColumns) ;
    arrayofDeltaEz[k]    =   new TMatrixD(kRows,kColumns) ;
  }
  
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
  
  const Int_t   symmetry = 0;
 
  // Set bondaries and solve Poisson's equation --------------------------
  
  if ( !fInitLookUp ) {
    
    AliInfo(Form("Solving the poisson equation (~ %d sec)",2*10*(int)(kPhiSlices/10)));
    
    for ( Int_t side = 0 ; side < 2 ; side++ ) {  // Solve Poisson3D twice; once for +Z and once for -Z
      
      for ( Int_t k = 0 ; k < kPhiSlices ; k++ )  {
	TMatrixD &arrayV    =  *arrayofArrayV[k] ;
	TMatrixD &charge    =  *arrayofCharge[k] ;
	
	//Fill arrays with initial conditions.  V on the boundary and Charge in the volume.
	for ( Int_t i = 0 ; i < kRows ; i++ ) {
	  for ( Int_t j = 0 ; j < kColumns ; j++ ) {  // Fill Vmatrix with Boundary Conditions
	    arrayV(i,j) = 0.0 ; 
	    charge(i,j) = 0.0 ;

	    Float_t radius0 = rlist[i] ;
	    Float_t phi0    = gridSizePhi * k ;
	    
	    // To avoid problems at sector boundaries, use an average of +- 1 degree from actual phi location
	    if ( j == (kColumns-1) ) {
	      arrayV(i,j) = 0.5*  ( GetROCVoltOffset( side, radius0, phi0+0.02 ) + GetROCVoltOffset( side, radius0, phi0-0.02 ) ) ;

	      if (side==1) // C side
		arrayV(i,j) = -arrayV(i,j); // minus sign on the C side to allow a consistent usage of global z when setting the boundaries
	    }
	  }
	}      
	
	for ( Int_t i = 1 ; i < kRows-1 ; i++ ) { 
	  for ( Int_t j = 1 ; j < kColumns-1 ; j++ ) {
	    charge(i,j)  =  0.0 ;
	  }
	}
      }      
      
      // Solve Poisson's equation in 3D cylindrical coordinates by relaxation technique
      // Allow for different size grid spacing in R and Z directions
      
      PoissonRelaxation3D( arrayofArrayV, arrayofCharge, 
			   arrayofEroverEz, arrayofEphioverEz, arrayofDeltaEz,
			   kRows, kColumns, kPhiSlices, gridSizePhi, kIterations, 
			   symmetry, fROCdisplacement) ;
      
      
      //Interpolate results onto a custom grid which is used just for these calculations.
      Double_t  r, phi, z ;
      for ( Int_t k = 0 ; k < kNPhi ; k++ ) {
	phi = fgkPhiList[k] ;
	
	TMatrixF &erOverEz   =  *fLookUpErOverEz[k]  ;
	TMatrixF &ephiOverEz =  *fLookUpEphiOverEz[k];
	TMatrixF &deltaEz    =  *fLookUpDeltaEz[k]   ;
	
	for ( Int_t j = 0 ; j < kNZ ; j++ ) {

	  z = TMath::Abs(fgkZList[j]) ;  // Symmetric solution in Z that depends only on ABS(Z)
  
	  if ( side == 0 &&  fgkZList[j] < 0 ) continue; // Skip rest of this loop if on the wrong side
	  if ( side == 1 &&  fgkZList[j] > 0 ) continue; // Skip rest of this loop if on the wrong side
	  
	  for ( Int_t i = 0 ; i < kNR ; i++ ) { 
	    r = fgkRList[i] ;

	    // Interpolate basicLookup tables; once for each rod, then sum the results
	    erOverEz(i,j)   = Interpolate3DTable(order, r, z, phi, kRows, kColumns, kPhiSlices, 
						 rlist, zedlist, philist, arrayofEroverEz  );
	    ephiOverEz(i,j) = Interpolate3DTable(order, r, z, phi, kRows, kColumns, kPhiSlices,
						 rlist, zedlist, philist, arrayofEphioverEz);
	    deltaEz(i,j)    = Interpolate3DTable(order, r, z, phi, kRows, kColumns, kPhiSlices,
						 rlist, zedlist, philist, arrayofDeltaEz  );

	    if (side == 1)  deltaEz(i,j) = -  deltaEz(i,j); // negative coordinate system on C side

	  } // end r loop
	}// end z loop
      }// end phi loop

      if ( side == 0 ) AliInfo(" A side done");
      if ( side == 1 ) AliInfo(" C side done");
    } // end side loop
  }
  
  // clear the temporary arrays lists
  for ( Int_t k = 0 ; k < kPhiSlices ; k++ )  {
    delete arrayofArrayV[k];
    delete arrayofCharge[k];
    delete arrayofEroverEz[k];  
    delete arrayofEphioverEz[k];
    delete arrayofDeltaEz[k];
  }
 

  fInitLookUp = kTRUE;

}


Float_t AliTPCROCVoltError3D::GetROCVoltOffset(Int_t side, Float_t r0, Float_t phi0) {
  // 
  // Returns the dz alignment data (in voltage equivalents) at 
  // the given position
  //

  Float_t xp = r0*TMath::Cos(phi0);
  Float_t yp = r0*TMath::Sin(phi0);
  
  // phi0 should be between 0 and 2pi 
  if (phi0<0) phi0+=TMath::TwoPi();
  Int_t roc = (Int_t)TMath::Floor((TMath::RadToDeg()*phi0)/20);
  if (side==1) roc+=18; // C side
  if (r0>132) roc+=36;  // OROC 
  
  // linear-plane data:  z = z0 + kx*lx + ky*ly (rotation in local coordinates)
  TMatrixD &fitData = *fdzDataLinFit;

  // local coordinates                                                          
  Double_t secAlpha = TMath::DegToRad()*(10.+20.*(((Int_t)roc)%18));
  Float_t lx = xp*TMath::Cos(secAlpha)+yp*TMath::Sin(secAlpha);
  Float_t ly = yp*TMath::Cos(secAlpha)-xp*TMath::Sin(secAlpha);

  // reference of rotation in lx is at the intersection between OROC and IROC
  // necessary, since the Fitprozedure is otherwise useless
  
  AliTPCROC * rocInfo = AliTPCROC::Instance();
  Double_t lxRef  = (rocInfo->GetPadRowRadii(0,62)+rocInfo->GetPadRowRadii(36,0))/2;
  
  Float_t dz = fitData(roc,0)+fitData(roc,1)*(lx-lxRef) + fitData(roc,2)*ly; // value in cm

  // aproximated Voltage-offset-aquivalent to the z misalignment
  // (linearly scaled with the z position)
  Double_t ezField = (fgkCathodeV-fgkGG)/fgkTPCZ0; // = ALICE Electric Field (V/cm) Magnitude ~ -400 V/cm; 
  Float_t voltOff = dz*ezField;            // values in "Volt equivalents"

  return voltOff;
}

TH2F * AliTPCROCVoltError3D::CreateHistoOfZAlignment(Int_t side, Int_t nx, Int_t ny) {
  //
  // return a simple histogramm containing the input to the poisson solver
  // (z positions of the Read-out chambers, linearly interpolated)

  char hname[100];
  if (side==0) snprintf(hname,100,"survey_dz_Aside");
  if (side==1) snprintf(hname,100,"survey_dz_Cside");

  TH2F *h = new TH2F(hname,hname,nx,-250.,250.,ny,-250.,250.);

  for (Int_t iy=1;iy<=ny;++iy) {
    Double_t yp = h->GetYaxis()->GetBinCenter(iy);
    for (Int_t ix=1;ix<=nx;++ix) {
      Double_t xp = h->GetXaxis()->GetBinCenter(ix);
    
      Float_t r0 = TMath::Sqrt(xp*xp+yp*yp);
      Float_t phi0 = TMath::ATan2(yp,xp); 
   
      Float_t dz = GetROCVoltOffset(side,r0,phi0); // in [volt]

      Double_t ezField = (fgkCathodeV-fgkGG)/fgkTPCZ0; // = ALICE Electric Field (V/cm) Magnitude ~ -400 V/cm; 
      dz = dz/ezField;    // in [cm]

      if (85.<=r0 && r0<=245.) {
	h->SetBinContent(ix,iy,dz); 
      } else {
	h->SetBinContent(ix,iy,0.);
      }
    }
  }
  
  h->GetXaxis()->SetTitle("x [cm]");
  h->GetYaxis()->SetTitle("y [cm]");
  h->GetZaxis()->SetTitle("dz [cm]");
  h->SetStats(0);
  //  h->DrawCopy("colz");

  return h;
} 

void AliTPCROCVoltError3D::Print(const Option_t* option) const {
  //
  // Print function to check the settings of the Rod shifts and the rotated clips
  // option=="a" prints the C0 and C1 coefficents for calibration purposes
  //

  TString opt = option; opt.ToLower();
  printf("%s\n",GetTitle());
  printf(" - z aligmnet of the TPC Read-Out chambers \n");
  printf("   (linear interpolation within the chamber:  dz = z0 + kx*(lx-133) + ky*ly [cm] ) \n");
  printf("   Info: Check the following data-file for more details: %s \n",fROCDataFileName.Data());

  if (opt.Contains("a")) { // Print all details
    TMatrixD &fitData = *fdzDataLinFit;
    printf(" A side:  IROC   ROCX=(z0,kx,ky): \n");
    for (Int_t roc = 0; roc<18; roc++) 
      printf("ROC%d:(%.2e,%.2e,%.2e) ",roc,fitData(roc,0),fitData(roc,1),fitData(roc,2));
    printf("\n A side:  OROC   ROCX=(z0,kx,ky): \n");
    for (Int_t roc = 36; roc<54; roc++) 
      printf("ROC%d:(%.2e,%.2e,%.2e) ",roc,fitData(roc,0),fitData(roc,1),fitData(roc,2));
    printf("\n C side:  IROC   ROCX=(z0,kx,ky): \n");
    for (Int_t roc = 18; roc<36; roc++) 
      printf("ROC%d:(%.2e,%.2e,%.2e) ",roc,fitData(roc,0),fitData(roc,1),fitData(roc,2));
    printf("\n C side:  OROC   ROCX=(z0,kx,ky): \n");
    for (Int_t roc = 54; roc<72; roc++) 
      printf("ROC%d:(%.2e,%.2e,%.2e) ",roc,fitData(roc,0),fitData(roc,1),fitData(roc,2));
    printf("\n\n");
    printf(" - T1: %1.4f, T2: %1.4f \n",fT1,fT2);
    printf(" - C1: %1.4f, C0: %1.4f \n",fC1,fC0);
  }

  if (!fInitLookUp) AliError("Lookup table was not initialized! You should do InitROCVoltError3D() ...");

}
