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
//  RecPoint implementation for PHOS-EMC 
//  An EmcRecPoint is a cluster of digits   
//           
//*-- Author: Dmitri Peressounko (RRC KI & SUBATECH)


// --- ROOT system ---
#include "TPad.h"
#include "TH2.h"
#include "TMath.h" 
#include "TCanvas.h" 

// --- Standard library ---

#include <iostream.h> 

// --- AliRoot header files ---

 #include "AliGenerator.h"
#include "AliPHOSGeometry.h"
#include "AliPHOSEmcRecPoint.h"
#include "AliRun.h"

ClassImp(AliPHOSEmcRecPoint)

//____________________________________________________________________________
AliPHOSEmcRecPoint::AliPHOSEmcRecPoint() : AliPHOSRecPoint()
{
  // ctor

  fMulDigit   = 0 ;  
  fAmp   = 0. ;   
  fCoreEnergy = 0 ; 
  fEnergyList = 0 ;
  fGeom = AliPHOSGeometry::GetInstance() ;
  fLocPos.SetX(1000000.)  ;      //Local position should be evaluated
  
}

//____________________________________________________________________________
AliPHOSEmcRecPoint::~AliPHOSEmcRecPoint()
{
  // dtor
  if ( fEnergyList )
    delete[] fEnergyList ; 
}

//____________________________________________________________________________
void AliPHOSEmcRecPoint::AddDigit(AliPHOSDigit & digit, Float_t Energy)
{
  // Adds a digit to the RecPoint
  //  and accumulates the total amplitude and the multiplicity 
  
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

  EvalPHOSMod(&digit) ;
}

//____________________________________________________________________________
Bool_t AliPHOSEmcRecPoint::AreNeighbours(AliPHOSDigit * digit1, AliPHOSDigit * digit2 ) const
{
  // Tells if (true) or not (false) two digits are neighbors)
  
  Bool_t aren = kFALSE ;
  
  AliPHOSGeometry * phosgeom =  (AliPHOSGeometry *) fGeom ;
  Int_t relid1[4] ; 
  phosgeom->AbsToRelNumbering(digit1->GetId(), relid1) ; 

  Int_t relid2[4] ; 
  phosgeom->AbsToRelNumbering(digit2->GetId(), relid2) ; 
  
  Int_t rowdiff = TMath::Abs( relid1[2] - relid2[2] ) ;  
  Int_t coldiff = TMath::Abs( relid1[3] - relid2[3] ) ;  

  if (( coldiff <= 1 )  && ( rowdiff <= 1 ) && (coldiff + rowdiff > 0)) 
    aren = kTRUE ;
  
  return aren ;
}

//____________________________________________________________________________
Int_t AliPHOSEmcRecPoint::Compare(const TObject * obj) const
{
  // Compares two RecPoints according to their position in the PHOS modules

  Float_t delta = 1 ; //Width of "Sorting row". If you changibg this 
                      //value (what is senseless) change as vell delta in
                      //AliPHOSTrackSegmentMakerv* and other RecPoints...
  Int_t rv ; 

  AliPHOSEmcRecPoint * clu = (AliPHOSEmcRecPoint *)obj ; 

 
  Int_t phosmod1 = GetPHOSMod() ;
  Int_t phosmod2 = clu->GetPHOSMod() ;

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
void AliPHOSEmcRecPoint::ExecuteEvent(Int_t event, Int_t px, Int_t py)
{
//   Commented by Dmitri Peressounko: there is no possibility to ensure, 
//   that AliPHOSIndexToObject keeps the correct information.

//   // Execute action corresponding to one event
//   //  This member function is called when a AliPHOSRecPoint is clicked with the locator
//   //
//   //  If Left button is clicked on AliPHOSRecPoint, the digits are switched on    
//   //  and switched off when the mouse button is released.
//   //

//   //   static Int_t pxold, pyold;

//   AliPHOSIndexToObject * please =  AliPHOSIndexToObject::GetInstance() ; 
  
//   static TGraph *  digitgraph = 0 ;
  
//   if (!gPad->IsEditable()) return;
  
//   TH2F * histo = 0 ;
//   TCanvas * histocanvas ; 
  
//   switch (event) {
    
//   case kButton1Down: {
//     AliPHOSDigit * digit ;
//     AliPHOSGeometry * phosgeom =  (AliPHOSGeometry *) fGeom ;
//     Int_t iDigit;
//     Int_t relid[4] ;
    
//     const Int_t kMulDigit = AliPHOSEmcRecPoint::GetDigitsMultiplicity() ; 
//     Float_t * xi = new Float_t[kMulDigit] ; 
//     Float_t * zi = new Float_t[kMulDigit] ; 
    
//     // create the histogram for the single cluster 
//     // 1. gets histogram boundaries
//     Float_t ximax = -999. ; 
//     Float_t zimax = -999. ; 
//     Float_t ximin = 999. ; 
//     Float_t zimin = 999. ;
    
//     for(iDigit=0; iDigit<kMulDigit; iDigit++) {
//       digit = (AliPHOSDigit *) ( please->GimeDigit(fDigitsList[iDigit]) ) ;
//       phosgeom->AbsToRelNumbering(digit->GetId(), relid) ;
//       phosgeom->RelPosInModule(relid, xi[iDigit], zi[iDigit]);
//       if ( xi[iDigit] > ximax )
// 	ximax = xi[iDigit] ; 
//       if ( xi[iDigit] < ximin )
// 	ximin = xi[iDigit] ; 
//       if ( zi[iDigit] > zimax )
// 	zimax = zi[iDigit] ; 
//       if ( zi[iDigit] < zimin )
// 	zimin = zi[iDigit] ;     
//     }
//     ximax += phosgeom->GetCrystalSize(0) / 2. ;
//     zimax += phosgeom->GetCrystalSize(2) / 2. ;
//     ximin -= phosgeom->GetCrystalSize(0) / 2. ;
//     zimin -= phosgeom->GetCrystalSize(2) / 2. ;
//     Int_t xdim = (int)( (ximax - ximin ) / phosgeom->GetCrystalSize(0) + 0.5  ) ; 
//     Int_t zdim = (int)( (zimax - zimin ) / phosgeom->GetCrystalSize(2) + 0.5 ) ;
    
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
//       digit = (AliPHOSDigit *) ( please->GimeDigit(fDigitsList[iDigit]) ) ;
//       phosgeom->AbsToRelNumbering(digit->GetId(), relid) ;
//       phosgeom->RelPosInModule(relid, x, z);
//       histo->Fill(x, z, fEnergyList[iDigit] ) ;
//     }
    
//     if (!digitgraph) {
//       digitgraph = new TGraph(kMulDigit,xi,zi);
//       digitgraph-> SetMarkerStyle(5) ; 
//       digitgraph-> SetMarkerSize(1.) ;
//       digitgraph-> SetMarkerColor(1) ;
//       digitgraph-> Paint("P") ;
//     }
    
//     Print() ;
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
void  AliPHOSEmcRecPoint::EvalDispersion(Float_t logWeight,TClonesArray * digits)
{
  // Calculates the dispersion of the shower at the origine of the RecPoint

  Float_t d    = 0 ;
  Float_t wtot = 0 ;

  TVector3 locpos;
  GetLocalPosition(locpos);
  Float_t x = locpos.X() ;
  Float_t z = locpos.Z() ;

  AliPHOSDigit * digit ;
  AliPHOSGeometry * phosgeom =  (AliPHOSGeometry *) fGeom ;
  
  Int_t iDigit;
  for(iDigit=0; iDigit < fMulDigit; iDigit++) {
    digit = (AliPHOSDigit *) digits->At(fDigitsList[iDigit])  ;
    Int_t relid[4] ;
    Float_t xi ;
    Float_t zi ;
    phosgeom->AbsToRelNumbering(digit->GetId(), relid) ;
    phosgeom->RelPosInModule(relid, xi, zi);
    Float_t w = TMath::Max(0.,logWeight+TMath::Log(fEnergyList[iDigit]/fAmp ) ) ;
    d += w*((xi-x)*(xi-x) + (zi-z)*(zi-z) ) ; 
    wtot+=w ;
  }


  d /= wtot ;

  fDispersion = TMath::Sqrt(d) ;
}
//______________________________________________________________________________
void AliPHOSEmcRecPoint::EvalCoreEnergy(TClonesArray * digits)
{
  //This function calculates energy in the core, 
  //i.e. within radius rad = 3cm. Beyond this radius
  //in accoradnce with shower profile energy deposition 
  // should be less than 2%

  Float_t coreRadius = 3 ;

  TVector3 locpos;
  GetLocalPosition(locpos);
  Float_t x = locpos.X() ;
  Float_t z = locpos.Z() ;

  AliPHOSDigit * digit ;
  AliPHOSGeometry * phosgeom =  (AliPHOSGeometry *) fGeom ;
  
  Int_t iDigit;
  for(iDigit=0; iDigit < fMulDigit; iDigit++) {
    digit = (AliPHOSDigit *) ( digits->At(fDigitsList[iDigit]) ) ;
    Int_t relid[4] ;
    Float_t xi ;
    Float_t zi ;
    phosgeom->AbsToRelNumbering(digit->GetId(), relid) ;
    phosgeom->RelPosInModule(relid, xi, zi);    
    Float_t distance = TMath::Sqrt((xi-x)*(xi-x)+(zi-z)*(zi-z)) ;
    if(distance < coreRadius)
      fCoreEnergy += fEnergyList[iDigit] ;
  }

}

//____________________________________________________________________________
void  AliPHOSEmcRecPoint::EvalElipsAxis(Float_t logWeight,TClonesArray * digits)
{
  // Calculates the axis of the shower ellipsoid

  Double_t wtot = 0. ;
  Double_t x    = 0.;
  Double_t z    = 0.;
  Double_t dxx  = 0.;
  Double_t dzz  = 0.;
  Double_t dxz  = 0.;

  AliPHOSDigit * digit ;
  AliPHOSGeometry * phosgeom =  (AliPHOSGeometry *) fGeom ;
  Int_t iDigit;

  for(iDigit=0; iDigit<fMulDigit; iDigit++) {
    digit = (AliPHOSDigit *) digits->At(fDigitsList[iDigit])  ;
    Int_t relid[4] ;
    Float_t xi ;
    Float_t zi ;
    phosgeom->AbsToRelNumbering(digit->GetId(), relid) ;
    phosgeom->RelPosInModule(relid, xi, zi);
    Double_t w = TMath::Max(0.,logWeight+TMath::Log(fEnergyList[iDigit]/fAmp ) ) ;
    dxx  += w * xi * xi ;
    x    += w * xi ;
    dzz  += w * zi * zi ;
    z    += w * zi ; 
    dxz  += w * xi * zi ; 
    wtot += w ;
  }
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
//   Double_t DistanceToIP= (Double_t ) ((AliPHOSGeometry *) fGeom)->GetIPtoCrystalSurface() ;
  
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
}

//____________________________________________________________________________
void AliPHOSEmcRecPoint::EvalAll(Float_t logWeight, TClonesArray * digits )
{
  AliPHOSRecPoint::EvalAll(logWeight,digits) ;
  EvalLocalPosition(logWeight, digits) ;
  EvalElipsAxis(logWeight, digits) ;
  EvalDispersion(logWeight, digits) ;
  EvalCoreEnergy(digits);
}
//____________________________________________________________________________
void AliPHOSEmcRecPoint::EvalLocalPosition(Float_t logWeight, TClonesArray * digits)
{
  // Calculates the center of gravity in the local PHOS-module coordinates 

  Float_t wtot = 0. ;
 
  Int_t relid[4] ;

  Float_t x = 0. ;
  Float_t z = 0. ;
  
  AliPHOSDigit * digit ;

  AliPHOSGeometry * phosgeom =  (AliPHOSGeometry *) fGeom ;

  Int_t iDigit;

  for(iDigit=0; iDigit<fMulDigit; iDigit++) {
    digit = (AliPHOSDigit *) digits->At(fDigitsList[iDigit]) ;

    Float_t xi ;
    Float_t zi ;
    phosgeom->AbsToRelNumbering(digit->GetId(), relid) ;
    phosgeom->RelPosInModule(relid, xi, zi);
    Float_t w = TMath::Max( 0., logWeight + TMath::Log( fEnergyList[iDigit] / fAmp ) ) ;
    x    += xi * w ;
    z    += zi * w ;
    wtot += w ;

  }

  x /= wtot ;
  z /= wtot ;

  // Correction for the depth of the shower starting point (TDR p 127)  
  Float_t para = 0.925 ; 
  Float_t parb = 6.52 ; 

  Float_t xo,yo,zo ; //Coordinates of the origin
  gAlice->Generator()->GetOrigin(xo,yo,zo) ;

  Float_t phi = phosgeom->GetPHOSAngle(relid[0]) ;

  //Transform to the local ref.frame
  Float_t xoL,yoL ;
  xoL = xo*TMath::Cos(phi)-yo*TMath::Sin(phi) ;
  yoL = xo*TMath::Sin(phi)+yo*TMath::Cos(phi) ;
  
  Float_t radius = TMath::Sqrt((xoL-x)*(xoL-x)+
                               (phosgeom->GetIPtoCrystalSurface()-yoL)*(phosgeom->GetIPtoCrystalSurface()-yoL)+
                               (zo-z)*(zo-z));
 
  Float_t incidencephi = TMath::ATan((x-xoL ) / radius) ; 
  Float_t incidencetheta = TMath::ATan((z-zo) / radius) ;
 
  Float_t depthx =  ( para * TMath::Log(fAmp) + parb ) * TMath::Sin(incidencephi) ; 
  Float_t depthz =  ( para * TMath::Log(fAmp) + parb ) * TMath::Sin(incidencetheta) ; 
  

  fLocPos.SetX(x - depthx)  ;
  fLocPos.SetY(0.) ;
  fLocPos.SetZ(z - depthz)  ;

}

//____________________________________________________________________________
Float_t AliPHOSEmcRecPoint::GetMaximalEnergy(void) const
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
Int_t AliPHOSEmcRecPoint::GetMultiplicityAtLevel(const Float_t H) const
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
Int_t  AliPHOSEmcRecPoint::GetNumberOfLocalMax(Int_t *  maxAt, Float_t * maxAtEnergy,
					       Float_t locMaxCut,TClonesArray * digits) const
{ 
  // Calculates the number of local maxima in the cluster using fLocalMaxCut as the minimum
  //  energy difference between two local maxima

  AliPHOSDigit * digit ;
  AliPHOSDigit * digitN ;
  

  Int_t iDigitN ;
  Int_t iDigit ;

  for(iDigit = 0; iDigit < fMulDigit; iDigit++)
    maxAt[iDigit] = (Int_t) digits->At(fDigitsList[iDigit])  ;

  
  for(iDigit = 0 ; iDigit < fMulDigit; iDigit++) {   
    if(maxAt[iDigit] != -1) {
      digit = (AliPHOSDigit *) maxAt[iDigit] ;
          
      for(iDigitN = 0; iDigitN < fMulDigit; iDigitN++) {	
	digitN = (AliPHOSDigit *) digits->At(fDigitsList[iDigitN]) ; 
	
	if ( AreNeighbours(digit, digitN) ) {
	  if (fEnergyList[iDigit] > fEnergyList[iDigitN] ) {    
	    maxAt[iDigitN] = -1 ;
	    // but may be digit too is not local max ?
	    if(fEnergyList[iDigit] < fEnergyList[iDigitN] + locMaxCut) 
	      maxAt[iDigit] = -1 ;
	  }
	  else {
	    maxAt[iDigit] = -1 ;
	    // but may be digitN too is not local max ?
	    if(fEnergyList[iDigit] > fEnergyList[iDigitN] - locMaxCut) 
	      maxAt[iDigitN] = -1 ; 
	  } 
	} // if Areneighbours
      } // while digitN
    } // slot not empty
  } // while digit
  
  iDigitN = 0 ;
  for(iDigit = 0; iDigit < fMulDigit; iDigit++) { 
    if(maxAt[iDigit] != -1){
      maxAt[iDigitN] = maxAt[iDigit] ;
      maxAtEnergy[iDigitN] = fEnergyList[iDigit] ;
      iDigitN++ ; 
    }
  }
  return iDigitN ;
}



//____________________________________________________________________________
void AliPHOSEmcRecPoint::Print(Option_t * option) 
{
  // Print the list of digits belonging to the cluster
  
  cout << "AliPHOSEmcRecPoint: " << endl ;

  Int_t iDigit;
  cout << " digits # = " ;
  for(iDigit=0; iDigit<fMulDigit; iDigit++)
    cout << fDigitsList[iDigit] << "  " ;  
  cout << endl ;
  
  cout << " Energies = " ;
  for(iDigit=0; iDigit<fMulDigit; iDigit++) 
    cout  << fEnergyList[iDigit] << "  ";
  cout << endl ;
  
  cout << "       Multiplicity    = " << fMulDigit  << endl ;
  cout << "       Cluster Energy  = " << fAmp << endl ;
  cout << "       Stored at position " << GetIndexInList() << endl ; 
 
}
 
  
