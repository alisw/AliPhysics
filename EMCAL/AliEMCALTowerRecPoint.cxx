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

//_________________________________________________________________________
//  RecPoint implementation for EMCAL-EMC 
//  An TowerRecPoint is a cluster of digits   
//*--
//*-- Author: Dmitri Peressounko (RRC KI & SUBATECH)


// --- ROOT system ---
#include "TMath.h" 

// --- Standard library ---

// --- AliRoot header files ---

#include "AliGenerator.h"
#include "AliEMCALGeometry.h"
#include "AliEMCALTowerRecPoint.h"
#include "AliEMCALGetter.h"

ClassImp(AliEMCALTowerRecPoint)

//____________________________________________________________________________
AliEMCALTowerRecPoint::AliEMCALTowerRecPoint() : AliEMCALRecPoint()
{
  // ctor

  fMulDigit   = 0 ;  
  fAmp   = 0. ;   
  fCoreEnergy = 0 ; 
  fEnergyList = 0 ;
  fTime = 0. ;
  fLocPos.SetX(0.)  ;      //Local position should be evaluated
}

//____________________________________________________________________________
AliEMCALTowerRecPoint::AliEMCALTowerRecPoint(const char * opt) : AliEMCALRecPoint(opt)
{
  // ctor
  
  fMulDigit   = 0 ;  
  fAmp   = 0. ;   
  fCoreEnergy = 0 ; 
  fEnergyList = 0 ;
  fTime = -1. ;
  fLocPos.SetX(1000000.)  ;      //Local position should be evaluated  
}

//____________________________________________________________________________
AliEMCALTowerRecPoint::~AliEMCALTowerRecPoint()
{
  // dtor

  if ( fEnergyList )
    delete[] fEnergyList ; 
}

//____________________________________________________________________________
void AliEMCALTowerRecPoint::AddDigit(AliEMCALDigit & digit, Float_t Energy)
{
  // Adds a digit to the RecPoint
  // and accumulates the total amplitude and the multiplicity 
  
  if(fEnergyList == 0)
    fEnergyList =  new Float_t[fMaxDigit]; 

  if ( fMulDigit >= fMaxDigit ) { // increase the size of the lists 
    fMaxDigit*=2 ; 
    Int_t * tempo = new ( Int_t[fMaxDigit] ) ; 
    Float_t * tempoE =  new ( Float_t[fMaxDigit] ) ;

    Int_t index ;     
    for ( index = 0 ; index < fMulDigit ; index++ ){
      tempo[index]  = fDigitsList[index] ;
      tempoE[index] = fEnergyList[index] ; 
    }
    
    delete [] fDigitsList ; 
    fDigitsList =  new ( Int_t[fMaxDigit] ) ;
 
    delete [] fEnergyList ;
    fEnergyList =  new ( Float_t[fMaxDigit] ) ;

    for ( index = 0 ; index < fMulDigit ; index++ ){
      fDigitsList[index] = tempo[index] ;
      fEnergyList[index] = tempoE[index] ; 
    }
 
    delete [] tempo ;
    delete [] tempoE ; 
  } // if
  
  fDigitsList[fMulDigit]   = digit.GetIndexInList()  ; 
  fEnergyList[fMulDigit]   = Energy ;
  fMulDigit++ ; 
  fAmp += Energy ; 

  // EvalEMCALMod(&digit) ;
}

//____________________________________________________________________________
Bool_t AliEMCALTowerRecPoint::AreNeighbours(AliEMCALDigit * digit1, AliEMCALDigit * digit2 ) const
{
  // Tells if (true) or not (false) two digits are neighbors
  
  Bool_t aren = kFALSE ;
  
  AliEMCALGeometry * phosgeom =  (AliEMCALGetter::Instance())->EMCALGeometry();

  Int_t relid1[3] ; 
  phosgeom->AbsToRelNumbering(digit1->GetId(), relid1) ; 

  Int_t relid2[3] ; 
  phosgeom->AbsToRelNumbering(digit2->GetId(), relid2) ; 
  
  Int_t rowdiff = TMath::Abs( relid1[1] - relid2[1] ) ;  
  Int_t coldiff = TMath::Abs( relid1[2] - relid2[2] ) ;  

  if (( coldiff <= 1 )  && ( rowdiff <= 1 ) && (coldiff + rowdiff > 0)) 
    aren = kTRUE ;
  
  return aren ;
}

//____________________________________________________________________________
Int_t AliEMCALTowerRecPoint::Compare(const TObject * obj) const
{
  // Compares two RecPoints according to their position in the EMCAL modules

  Float_t delta = 1 ; //Width of "Sorting row". If you change this 
                      //value (what is senseless) change as vell delta in
                      //AliEMCALTrackSegmentMakerv* and other RecPoints...
  Int_t rv ; 

  AliEMCALTowerRecPoint * clu = (AliEMCALTowerRecPoint *)obj ; 

 
  Int_t phosmod1 = GetEMCALArm() ;
  Int_t phosmod2 = clu->GetEMCALArm() ;

  TVector3 locpos1; 
  GetLocalPosition(locpos1) ;
  TVector3 locpos2;  
  clu->GetLocalPosition(locpos2) ;  

  if(phosmod1 == phosmod2 ) {
    Int_t rowdif = (Int_t)TMath::Ceil(locpos1.X()/delta)-(Int_t)TMath::Ceil(locpos2.X()/delta) ;
    if (rowdif> 0) 
      rv = 1 ;
     else if(rowdif < 0) 
       rv = -1 ;
    else if(locpos1.Z()>locpos2.Z()) 
      rv = -1 ;
    else 
      rv = 1 ; 
     }

  else {
    if(phosmod1 < phosmod2 ) 
      rv = -1 ;
    else 
      rv = 1 ;
  }

  return rv ; 
}
//______________________________________________________________________________
void AliEMCALTowerRecPoint::ExecuteEvent(Int_t /*event*/, Int_t, Int_t) const
{
  
  // Execute action corresponding to one event
  //  This member function is called when a AliEMCALRecPoint is clicked with the locator
  //
  //  If Left button is clicked on AliEMCALRecPoint, the digits are switched on    
  //  and switched off when the mouse button is released.
  
    
  //    AliEMCALGeometry * phosgeom =  (AliEMCALGetter::Instance())->EMCALGeometry();

//   static TGraph *  digitgraph = 0 ;
  
//   if (!gPad->IsEditable()) return;
  
//   TH2F * histo = 0 ;
//   TCanvas * histocanvas ; 

//   const TClonesArray * digits = gime->Digits() ;
  
//   switch (event) {
    
//   case kButton1Down: {
//     AliEMCALDigit * digit ;
//     Int_t iDigit;
//     Int_t relid[3] ;
    
//     const Int_t kMulDigit = AliEMCALTowerRecPoint::GetDigitsMultiplicity() ; 
//     Float_t * xi = new Float_t[kMulDigit] ; 
//     Float_t * zi = new Float_t[kMulDigit] ; 
    
//     // create the histogram for the single cluster 
//     // 1. gets histogram boundaries
//     Float_t ximax = -999. ; 
//     Float_t zimax = -999. ; 
//     Float_t ximin = 999. ; 
//     Float_t zimin = 999. ;
    
//     for(iDigit=0; iDigit<kMulDigit; iDigit++) {
//       digit = (AliEMCALDigit *) digits->At(fDigitsList[iDigit])  ;
//       emcalgeom->AbsToRelNumbering(digit->GetId(), relid) ;
//       emcalgeom->RelPosInModule(relid, xi[iDigit], zi[iDigit]);
//       if ( xi[iDigit] > ximax )
// 	ximax = xi[iDigit] ; 
//       if ( xi[iDigit] < ximin )
// 	ximin = xi[iDigit] ; 
//       if ( zi[iDigit] > zimax )
// 	zimax = zi[iDigit] ; 
//       if ( zi[iDigit] < zimin )
// 	zimin = zi[iDigit] ;     
//     }
//     ximax += emcalgeom->GetCrystalSize(0) / 2. ;
//     zimax += emcalgeom->GetCrystalSize(2) / 2. ;
//     ximin -= emcalgeom->GetCrystalSize(0) / 2. ;
//     zimin -= emcalgeom->GetCrystalSize(2) / 2. ;
//     Int_t xdim = (int)( (ximax - ximin ) / emcalgeom->GetCrystalSize(0) + 0.5  ) ; 
//     Int_t zdim = (int)( (zimax - zimin ) / emcalgeom->GetCrystalSize(2) + 0.5 ) ;
    
//     // 2. gets the histogram title
    
//     Text_t title[100] ; 
//     sprintf(title,"Energy=%1.2f GeV ; Digits ; %d ", GetEnergy(), GetDigitsMultiplicity()) ;
    
//     if (!histo) {
//       delete histo ; 
//       histo = 0 ; 
//     }
//     histo = new TH2F("cluster3D", title,  xdim, ximin, ximax, zdim, zimin, zimax)  ;
    
//     Float_t x, z ; 
//     for(iDigit=0; iDigit<kMulDigit; iDigit++) {
//       digit = (AliEMCALDigit *) digits->At(fDigitsList[iDigit])  ;
//       emcalgeom->AbsToRelNumbering(digit->GetId(), relid) ;
//       emcalgeom->RelPosInModule(relid, x, z);
//       histo->Fill(x, z, fEnergyList[iDigit] ) ;
//     }
    
//     if (!digitgraph) {
//       digitgraph = new TGraph(kMulDigit,xi,zi);
//       digitgraph-> SetMarkerStyle(5) ; 
//       digitgraph-> SetMarkerSize(1.) ;
//       digitgraph-> SetMarkerColor(1) ;
//       digitgraph-> Paint("P") ;
//     }
    
//     //    Print() ;
//     histocanvas = new TCanvas("cluster", "a single cluster", 600, 500) ; 
//     histocanvas->Draw() ; 
//     histo->Draw("lego1") ; 
    
//     delete[] xi ; 
//     delete[] zi ; 
    
//     break;
//   }
  
//   case kButton1Up: 
//     if (digitgraph) {
//       delete digitgraph  ;
//       digitgraph = 0 ;
//     }
//     break;
  
//    }
}

//____________________________________________________________________________
void  AliEMCALTowerRecPoint::EvalDispersion(Float_t logWeight,TClonesArray * digits)
{
  // Calculates the dispersion of the shower at the origine of the RecPoint

  Float_t d    = 0. ;
  Float_t wtot = 0. ;

  AliEMCALDigit * digit ;
 
  AliEMCALGeometry * emcalgeom = (AliEMCALGetter::Instance())->EMCALGeometry();
  

 // Calculates the center of gravity in the local EMCAL-module coordinates 
  
  Int_t iDigit;

  if (!fTheta || !fPhi ) 
    EvalGlobalPosition(logWeight, digits) ;
  
  const Float_t kDeg2Rad = TMath::DegToRad() ; 
    
  Float_t cyl_radius = 0 ;  

  if (IsInECA()) 
    cyl_radius = emcalgeom->GetIP2ECASection() ;
  else 
    Fatal("EvalDispersion", "Unexpected tower section!") ; 
  
  Float_t x =  fLocPos.X() ; 
  Float_t y =  fLocPos.Y() ; 
  Float_t z =  fLocPos.Z() ; 
  
  if (gDebug == 2) 
    printf("EvalDispersion: x,y,z = %f,%f,%f", x, y, z) ;

// Calculates the dispersion in coordinates 
  wtot = 0.;
  for(iDigit=0; iDigit < fMulDigit; iDigit++) {
    digit = (AliEMCALDigit *) digits->At(fDigitsList[iDigit])  ;
    Float_t thetai = 0. ;
    Float_t phii = 0.;
    emcalgeom->PosInAlice(digit->GetId(), thetai, phii);
    Float_t xi =  cyl_radius * TMath::Cos(phii * kDeg2Rad ) ;
    Float_t yi =  cyl_radius * TMath::Sin(phii * kDeg2Rad ) ; 
    Float_t zi =  cyl_radius / TMath::Tan(thetai * kDeg2Rad ) ; 

    if (gDebug == 2) 
      printf("EvalDispersion: id = %d, xi,yi,zi = %f,%f,%f", digit->GetId(), xi, yi, zi) ;

    Float_t w = TMath::Max(0.,logWeight+TMath::Log(fEnergyList[iDigit]/fAmp ) ) ;
    d += w * ( (xi-x)*(xi-x) + (zi-z)*(zi-z) ) ; 
    wtot+=w ;
  }
  
  if ( wtot > 0 ) 
    d /= wtot ;
  else 
    d = 0. ; 

  fDispersion = TMath::Sqrt(d) ;
 
}
//______________________________________________________________________________
void AliEMCALTowerRecPoint::EvalCoreEnergy(Float_t logWeight, TClonesArray * digits)
{
  // This function calculates energy in the core, 
  // i.e. within a radius rad = 3cm around the center. Beyond this radius
  // in accordance with shower profile the energy deposition 
  // should be less than 2%

  Float_t coreRadius = 10. ;

  AliEMCALDigit * digit ;
  Float_t wtot = 0. ;

  AliEMCALGeometry * emcalgeom = (AliEMCALGetter::Instance())->EMCALGeometry();    
  Int_t iDigit;

  if (!fTheta || !fPhi ) {
    for(iDigit=0; iDigit<fMulDigit; iDigit++) {
      digit = dynamic_cast<AliEMCALDigit *>(digits->At(fDigitsList[iDigit])) ;
      
      Float_t thetai ;
      Float_t phii ;
      emcalgeom->PosInAlice(digit->GetId(), thetai, phii);
      Float_t w = TMath::Max( 0., logWeight + TMath::Log( fEnergyList[iDigit] / fAmp ) ) ;
      fTheta = fTheta + thetai * w ;
      fPhi   += (phii * w );
      wtot  += w ;
    }
    
    if (wtot > 0 ) { 
      fTheta /= wtot ;
      fPhi   /= wtot ;
    } else { 
      fTheta = -1 ; 
      fPhi   = -1 ; 
    }
  }
  
  const Float_t kDeg2Rad = TMath::DegToRad() ; 
  
  Float_t cyl_radius = emcalgeom->GetIP2ECASection();
  Float_t x =  cyl_radius * TMath::Cos(fPhi * kDeg2Rad ) ;
  Float_t y =  cyl_radius * TMath::Cos(fPhi * kDeg2Rad ) ; 
  Float_t z =  cyl_radius * TMath::Tan(fTheta * kDeg2Rad ) ; 

  for(iDigit=0; iDigit < fMulDigit; iDigit++) {
    digit = (AliEMCALDigit *) ( digits->At(fDigitsList[iDigit]) ) ;
    Float_t thetai = 0. ;
    Float_t phii = 0. ;
    emcalgeom->PosInAlice(digit->GetId(), thetai, phii);
    
    Float_t xi =  cyl_radius * TMath::Cos(phii * kDeg2Rad ) ;
    Float_t yi =  cyl_radius * TMath::Sin(phii * kDeg2Rad ) ; 
    Float_t zi =  cyl_radius * TMath::Tan(thetai * kDeg2Rad ) ; 
     
    Float_t distance = TMath::Sqrt((xi-x)*(xi-x)+(yi-y)*(yi-y)+(zi-z)*(zi-z)) ;
    if(distance < coreRadius)
      fCoreEnergy += fEnergyList[iDigit] ;
  }
  
}

//____________________________________________________________________________
void  AliEMCALTowerRecPoint::EvalElipsAxis(Float_t logWeight,TClonesArray * digits)
{
  // Calculates the axis of the shower ellipsoid

  Double_t wtot = 0. ;
  Double_t x    = 0.;
  Double_t z    = 0.;
  Double_t dxx  = 0.;
  Double_t dzz  = 0.;
  Double_t dxz  = 0.;

  AliEMCALDigit * digit ;

  AliEMCALGeometry * emcalgeom = (AliEMCALGetter::Instance())->EMCALGeometry();

  Int_t iDigit;
  const Float_t kDeg2Rad = TMath::DegToRad() ; 
  
   Float_t cyl_radius = 0 ;  
  
  if (IsInECA()) 
    cyl_radius = emcalgeom->GetIP2ECASection() ;
  else 
    Fatal("EvalDispersion", "Unexpected tower section!") ; 

  for(iDigit=0; iDigit<fMulDigit; iDigit++) {
    digit = (AliEMCALDigit *) digits->At(fDigitsList[iDigit])  ;
    Float_t thetai = 0. ;
    Float_t phii = 0. ; 
    emcalgeom->PosInAlice(digit->GetId(), thetai, phii);
    Double_t w = TMath::Max(0.,logWeight+TMath::Log(fEnergyList[iDigit]/fAmp ) ) ;
    Float_t xi =  cyl_radius * TMath::Cos(fPhi * kDeg2Rad ) ;
    Float_t zi =  cyl_radius / TMath::Tan(fTheta * kDeg2Rad ) ; 
    dxx  += w * xi * xi ;
    x    += w * xi ;
    dzz  += w * zi * zi ;
    z    += w * zi ; 
    dxz  += w * xi * zi ; 
    wtot += w ;
  }
  if ( wtot > 0 ) { 
    dxx /= wtot ;
    x   /= wtot ;
    dxx -= x * x ;
    dzz /= wtot ;
    z   /= wtot ;
    dzz -= z * z ;
    dxz /= wtot ;
    dxz -= x * z ;

    
    //   //Apply correction due to non-perpendicular incidence
//   Double_t CosX ;
//   Double_t CosZ ;
//   AliEMCALGeometry * emcalgeom = (AliEMCALGetter::Instance())->EMCALGeometry();
  //   Double_t DistanceToIP= (Double_t ) emcalgeom->GetIPDistance() ;
  
//   CosX = DistanceToIP/TMath::Sqrt(DistanceToIP*DistanceToIP+x*x) ;
//   CosZ = DistanceToIP/TMath::Sqrt(DistanceToIP*DistanceToIP+z*z) ;

//   dxx = dxx/(CosX*CosX) ;
//   dzz = dzz/(CosZ*CosZ) ;
//   dxz = dxz/(CosX*CosZ) ;


    fLambda[0] =  0.5 * (dxx + dzz) + TMath::Sqrt( 0.25 * (dxx - dzz) * (dxx - dzz) + dxz * dxz )  ;
    if(fLambda[0] > 0)
      fLambda[0] = TMath::Sqrt(fLambda[0]) ;
    
    fLambda[1] =  0.5 * (dxx + dzz) - TMath::Sqrt( 0.25 * (dxx - dzz) * (dxx - dzz) + dxz * dxz )  ;
    if(fLambda[1] > 0) //To avoid exception if numerical errors lead to negative lambda.
      fLambda[1] = TMath::Sqrt(fLambda[1]) ;
    else
      fLambda[1]= 0. ;
  } else { 
    fLambda[0]= 0. ;
    fLambda[1]= 0. ;
  }
}

//____________________________________________________________________________
void AliEMCALTowerRecPoint::EvalAll(Float_t logWeight, TClonesArray * digits )
{
  // Evaluates all shower parameters

  AliEMCALRecPoint::EvalAll(logWeight,digits) ;
  EvalGlobalPosition(logWeight, digits) ;
  EvalElipsAxis(logWeight, digits) ;
  EvalDispersion(logWeight, digits) ;
  EvalCoreEnergy(logWeight, digits);
  EvalTime(digits) ;
}

//____________________________________________________________________________
void AliEMCALTowerRecPoint::EvalGlobalPosition(Float_t logWeight, TClonesArray * digits)
{
  // Calculates the center of gravity in the local EMCAL-module coordinates 
  Float_t wtot = 0. ;
 
  //  Int_t relid[3] ;
  
  AliEMCALDigit * digit ;
  AliEMCALGeometry * emcalgeom  =  (AliEMCALGetter::Instance())->EMCALGeometry();
  Int_t iDigit;

  for(iDigit=0; iDigit<fMulDigit; iDigit++) {
    digit = dynamic_cast<AliEMCALDigit *>(digits->At(fDigitsList[iDigit])) ;

    Float_t thetai ;
    Float_t phii ;
    emcalgeom->PosInAlice(digit->GetId(), thetai, phii);
    Float_t w = TMath::Max( 0., logWeight + TMath::Log( fEnergyList[iDigit] / fAmp ) ) ;
    fTheta = fTheta + thetai * w ;
    fPhi   += (phii * w );
    wtot  += w ;
  }

  if ( wtot > 0 ) { 
    fTheta /= wtot ;
    fPhi   /= wtot ;
  } else {
    fTheta = -1 ; 
    fPhi   = -1.; 
  }
  

  const Float_t kDeg2Rad = TMath::DegToRad() ; 
  
  Float_t cyl_radius = 0 ;  

  if (IsInECA()) 
    cyl_radius = emcalgeom->GetIP2ECASection() ;
  else 
    Fatal("EvalGlobalPosition", "Unexpected tower section!") ; 
  
  Float_t x =  cyl_radius * TMath::Cos(fPhi * kDeg2Rad ) ;
  Float_t y =  cyl_radius * TMath::Sin(fPhi * kDeg2Rad ) ; 
  Float_t z =  cyl_radius / TMath::Tan(fTheta * kDeg2Rad ) ; 
  
  fLocPos.SetX(x)  ;
  fLocPos.SetY(y) ;
  fLocPos.SetZ(z)  ;
    
  if (gDebug==2)
    printf("EvalGlobalPosition: x,y,z = %f,%f,%f", fLocPos.X(), fLocPos.Y(), fLocPos.Z()) ; 
  fLocPosM = 0 ;
}

//____________________________________________________________________________
Float_t AliEMCALTowerRecPoint::GetMaximalEnergy(void) const
{
  // Finds the maximum energy in the cluster
  
  Float_t menergy = 0. ;

  Int_t iDigit;

  for(iDigit=0; iDigit<fMulDigit; iDigit++) {
 
    if(fEnergyList[iDigit] > menergy) 
      menergy = fEnergyList[iDigit] ;
  }
  return menergy ;
}

//____________________________________________________________________________
Int_t AliEMCALTowerRecPoint::GetMultiplicityAtLevel(const Float_t H) const
{
  // Calculates the multiplicity of digits with energy larger than H*energy 
  
  Int_t multipl   = 0 ;
  Int_t iDigit ;
  for(iDigit=0; iDigit<fMulDigit; iDigit++) {

    if(fEnergyList[iDigit] > H * fAmp) 
      multipl++ ;
  }
  return multipl ;
}

//____________________________________________________________________________
Int_t  AliEMCALTowerRecPoint::GetNumberOfLocalMax(AliEMCALDigit **  maxAt, Float_t * maxAtEnergy,
					       Float_t locMaxCut,TClonesArray * digits) const
{ 
  // Calculates the number of local maxima in the cluster using fLocalMaxCut as the minimum
  // energy difference between two local maxima

  AliEMCALDigit * digit ;
  AliEMCALDigit * digitN ;
  
  Int_t iDigitN ;
  Int_t iDigit ;

  for(iDigit = 0; iDigit < fMulDigit; iDigit++)
    maxAt[iDigit] = (AliEMCALDigit*) digits->At(fDigitsList[iDigit])  ;
  
  for(iDigit = 0 ; iDigit < fMulDigit; iDigit++) {   
    if(maxAt[iDigit]) {
      digit = maxAt[iDigit] ;
          
      for(iDigitN = 0; iDigitN < fMulDigit; iDigitN++) {	
	digitN = (AliEMCALDigit *) digits->At(fDigitsList[iDigitN]) ; 
	
	if ( AreNeighbours(digit, digitN) ) {
	  if (fEnergyList[iDigit] > fEnergyList[iDigitN] ) {    
	    maxAt[iDigitN] = 0 ;
	    // but may be digit too is not local max ?
	    if(fEnergyList[iDigit] < fEnergyList[iDigitN] + locMaxCut) 
	      maxAt[iDigit] = 0 ;
	  }
	  else {
	    maxAt[iDigit] = 0 ;
	    // but may be digitN too is not local max ?
	    if(fEnergyList[iDigit] > fEnergyList[iDigitN] - locMaxCut) 
	      maxAt[iDigitN] = 0 ; 
	  } 
	} // if Areneighbours
      } // while digitN
    } // slot not empty
  } // while digit
  
  iDigitN = 0 ;
  for(iDigit = 0; iDigit < fMulDigit; iDigit++) { 
    if(maxAt[iDigit] ){
      maxAt[iDigitN] = maxAt[iDigit] ;
      maxAtEnergy[iDigitN] = fEnergyList[iDigit] ;
      iDigitN++ ; 
    }
  }
  return iDigitN ;
}
//____________________________________________________________________________
void AliEMCALTowerRecPoint::EvalTime(TClonesArray * digits){
  
  Float_t maxE = 0;
  Int_t maxAt = 0;
  for(Int_t idig=0; idig < fMulDigit; idig++){
    if(fEnergyList[idig] > maxE){
      maxE = fEnergyList[idig] ;
      maxAt = idig;
    }
  }
  fTime = ((AliEMCALDigit*) digits->At(fDigitsList[maxAt]))->GetTime() ;
  
}
//____________________________________________________________________________
void AliEMCALTowerRecPoint::Print(Option_t *) 
{
  // Print the list of digits belonging to the cluster
  
  printf("\n") ; 

  Int_t iDigit;
  printf("digits # = ");
  for(iDigit=0; iDigit<fMulDigit; iDigit++) {
     printf("%i ", fDigitsList[iDigit]); 
  } 
  
  printf("\nEnergies = ");
  for(iDigit=0; iDigit<fMulDigit; iDigit++) { 
    printf("%f ", fEnergyList[iDigit]); 
  }
  
   printf("\nPrimaries  ");
   for(iDigit = 0;iDigit < fMulTrack; iDigit++) {
     printf("%i ", fTracksList[iDigit]);
   }
   printf("\n       Multiplicity    = %i", fMulDigit);
   printf("\n       Cluster Energy  = %f", fAmp);
   printf("\n       Number of primaries %i", fMulTrack);
   printf("\n       Stored at position: %i", GetIndexInList()); 
}
 
//____________________________________________________________________________
const TVector3 AliEMCALTowerRecPoint::XYZInAlice(Float_t r, Float_t theta, Float_t phi) const 
{
  // spherical coordinates of recpoint in Alice reference frame

  if (gDebug == 2) 
    printf("XYZInAlice: r = %f, theta = %f, phi = %f", r, theta, phi) ; 

  if (theta == 9999. || phi == 9999. || r == 9999.) {
    TVector3  globalpos;  
    GetGlobalPosition(globalpos);
    phi   =  globalpos.X() * TMath::DegToRad() ; 
    r     =  globalpos.Y() ; 
    theta =  globalpos.Z() * TMath::DegToRad() ; 
  }
  else {
    theta *= TMath::DegToRad() ; 
    phi   *= TMath::DegToRad() ; 
  }
  
  Float_t y = r * TMath::Cos(phi) ; 
  Float_t x = r * TMath::Sin(phi) * TMath::Sin(theta) ; 
  Float_t z = r * TMath::Sin(phi) * TMath::Cos(theta) ; 
  
  TVector3 vec(z, x, y) ; 
  return vec ; 
}  
