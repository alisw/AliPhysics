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

//_________________________________________________________________________
// Rec Point in the PHOS EM calorimeter 
//*-- Author : Dmitri Peressounko RRC KI
//////////////////////////////////////////////////////////////////////////////

// --- ROOT system ---
#include "TPad.h"
#include "TH2.h"
#include "TMath.h" 
#include "TCanvas.h" 

// --- Standard library ---

#include <iostream>

// --- AliRoot header files ---

#include "AliPHOSGeometry.h"
#include "AliPHOSEmcRecPoint.h"
#include "AliRun.h"

ClassImp(AliPHOSEmcRecPoint)

//____________________________________________________________________________
AliPHOSEmcRecPoint::AliPHOSEmcRecPoint(Float_t W0, Float_t LocMaxCut)
  : AliPHOSRecPoint()
{
  // ctor

  fMulDigit   = 0 ;  
  fAmp   = 0. ;   
  fEnergyList = new Float_t[fMaxDigit]; 
  AliPHOSGeometry * PHOSGeom =  (AliPHOSGeometry *) fGeom ;
  fDelta     =  PHOSGeom->GetCrystalSize(0) ; 
  fW0        = W0 ;          
  fLocMaxCut = LocMaxCut ; 
  fLocPos.SetX(1000000.)  ;      //Local position should be evaluated
}

//____________________________________________________________________________
AliPHOSEmcRecPoint::~AliPHOSEmcRecPoint() 
{
  // dtor 
}

//____________________________________________________________________________
void AliPHOSEmcRecPoint::AddDigit(AliDigitNew & digit, Float_t Energy)
{
  // adds a digit to the digits list
  // and accumulates the total amplitude and the multiplicity 
  
  if ( fMulDigit >= fMaxDigit ) { // increase the size of the lists 
    fMaxDigit*=2 ; 
    int * tempo = new ( int[fMaxDigit] ) ; 
    Float_t * tempoE =  new ( Float_t[fMaxDigit] ) ;

    Int_t index ;     
    for ( index = 0 ; index < fMulDigit ; index++ ){
      tempo[index] = fDigitsList[index] ;
      tempoE[index] = fEnergyList[index] ; 
    }
    
    delete [] fDigitsList ; 
    fDigitsList =  new ( int[fMaxDigit] ) ;
 
    delete [] fEnergyList ;
    fEnergyList =  new ( Float_t[fMaxDigit] ) ;

    for ( index = 0 ; index < fMulDigit ; index++ ){
      fDigitsList[index] = tempo[index] ;
      fEnergyList[index] = tempoE[index] ; 
    }
 
    delete [] tempo ;
    delete [] tempoE ; 
  } // if
  
  fDigitsList[fMulDigit]   =  (int) &digit  ; 
  fEnergyList[fMulDigit++] = Energy ;
  fAmp += Energy ; 
}

//____________________________________________________________________________
Bool_t AliPHOSEmcRecPoint::AreNeighbours(AliPHOSDigit * digit1, AliPHOSDigit * digit2 ) 
{
  
  Bool_t aren = kFALSE ;
  
  AliPHOSGeometry * PHOSGeom =  (AliPHOSGeometry *) fGeom ;
  Int_t relid1[4] ; 
  PHOSGeom->AbsToRelNumbering(digit1->GetId(), relid1) ; 

  Int_t relid2[4] ; 
  PHOSGeom->AbsToRelNumbering(digit2->GetId(), relid2) ; 
  
  Int_t RowDiff = TMath::Abs( relid1[2] - relid2[2] ) ;  
  Int_t ColDiff = TMath::Abs( relid1[3] - relid2[3] ) ;  

  if (( ColDiff<=1 )  && ( RowDiff <= 1 ) && (ColDiff+RowDiff > 0)) 
    aren = kTRUE ;
  
  return aren ;
}

//____________________________________________________________________________
Int_t AliPHOSEmcRecPoint::Compare(TObject * obj)
{
  Int_t rv ; 

  AliPHOSEmcRecPoint * clu = (AliPHOSEmcRecPoint *)obj ; 

 
  Int_t PHOSMod1 = this->GetPHOSMod() ;
  Int_t PHOSMod2 = clu->GetPHOSMod() ;

  TVector3 LocPos1; 
  this->GetLocalPosition(LocPos1) ;
  TVector3 LocPos2;  
  clu->GetLocalPosition(LocPos2) ;  

  if(PHOSMod1 == PHOSMod2 ) {
    Int_t rowdif = (Int_t)TMath::Ceil(LocPos1.X()/fDelta)-(Int_t)TMath::Ceil(LocPos2.X()/fDelta) ;
    if (rowdif> 0) 
      rv = -1 ;
     else if(rowdif < 0) 
       rv = 1 ;
    else if(LocPos1.Z()>LocPos2.Z()) 
      rv = -1 ;
    else 
      rv = 1 ; 
     }

  else {
    if(PHOSMod1 < PHOSMod2 ) 
      rv = -1 ;
    else 
      rv = 1 ;
  }

  return rv ; 
}

//______________________________________________________________________________
void AliPHOSEmcRecPoint::ExecuteEvent(Int_t event, Int_t px, Int_t py)
{
  //Execute action corresponding to one event
  //  This member function is called when a AliPHOSRecPoint is clicked with the locator
  //
  //  If Left button is clicked on AliPHOSRecPoint, the digits are switched on    
  //  and switched off when the mouse button is released.
  //

  //   static Int_t pxold, pyold;

   static TGraph *  DigitGraph = 0 ;

   if (!gPad->IsEditable()) return;

   TH2F * Histo = 0 ;
   TCanvas * HistoCanvas ; 
   
   switch (event) {
   
   case kButton1Down: {
     AliPHOSDigit * digit ;
     AliPHOSGeometry * PHOSGeom =  (AliPHOSGeometry *) fGeom ;
     Int_t iDigit;
     Int_t relid[4] ;
     Float_t xi[fMulDigit] ;
     Float_t zi[fMulDigit] ;

     // create the histogram for the single cluster 
     // 1. gets histogram boundaries
     Float_t ximax = -999. ; 
     Float_t zimax = -999. ; 
     Float_t ximin = 999. ; 
     Float_t zimin = 999. ;
 
     for(iDigit=0; iDigit<fMulDigit; iDigit++) {
       digit = (AliPHOSDigit *) fDigitsList[iDigit];
       PHOSGeom->AbsToRelNumbering(digit->GetId(), relid) ;
       PHOSGeom->RelPosInModule(relid, xi[iDigit], zi[iDigit]);
       if ( xi[iDigit] > ximax )
	 ximax = xi[iDigit] ; 
       if ( xi[iDigit] < ximin )
	 ximin = xi[iDigit] ; 
       if ( zi[iDigit] > zimax )
	 zimax = zi[iDigit] ; 
       if ( zi[iDigit] < zimin )
	 zimin = zi[iDigit] ;     
     }
     ximax += PHOSGeom->GetCrystalSize(0) / 2. ;
     zimax += PHOSGeom->GetCrystalSize(2) / 2. ;
     ximin -= PHOSGeom->GetCrystalSize(0) / 2. ;
     zimin -= PHOSGeom->GetCrystalSize(2) / 2. ;
     Int_t xdim = (int)( (ximax - ximin ) / PHOSGeom->GetCrystalSize(0) + 0.5  ) ; 
     Int_t zdim = (int)( (zimax - zimin ) / PHOSGeom->GetCrystalSize(2) + 0.5 ) ;
 
     // 2. gets the histogram title

     Text_t title[100] ; 
     sprintf(title,"Energy=%1.2f GeV ; Digits ; %d ", GetEnergy(), GetDigitsMultiplicity()) ;
  
     if (!Histo) {
       delete Histo ; 
       Histo = 0 ; 
     }
     Histo = new TH2F("cluster3D", title,  xdim, ximin, ximax, zdim, zimin, zimax)  ;

     Float_t x, z ; 
     for(iDigit=0; iDigit<fMulDigit; iDigit++) {
       digit = (AliPHOSDigit *) fDigitsList[iDigit];
       PHOSGeom->AbsToRelNumbering(digit->GetId(), relid) ;
       PHOSGeom->RelPosInModule(relid, x, z);
       Histo->Fill(x, z, fEnergyList[iDigit] ) ;
     }

     if (!DigitGraph) {
       DigitGraph = new TGraph(fMulDigit,xi,zi);
       DigitGraph-> SetMarkerStyle(5) ; 
       DigitGraph-> SetMarkerSize(1.) ;
       DigitGraph-> SetMarkerColor(1) ;
       DigitGraph-> Paint("P") ;
     }

     Print() ;
     HistoCanvas = new TCanvas("cluser", "a single cluster", 600, 500) ; 
     HistoCanvas->Draw() ; 
     Histo->Draw("lego1") ; 

     break;
   }

   case kButton1Up: 
     if (DigitGraph) {
       delete DigitGraph  ;
       DigitGraph = 0 ;
     }
     break;
  
   }
}

//____________________________________________________________________________
Float_t  AliPHOSEmcRecPoint::GetDispersion() 
{
  Float_t D    = 0 ;
  Float_t wtot = 0 ;

  TVector3 LocPos;
  GetLocalPosition(LocPos);
  Float_t x = LocPos.X() ;
  Float_t z = LocPos.Z() ;
  //  Int_t i = GetPHOSMod() ;

  AliPHOSDigit * digit ;
  AliPHOSGeometry * PHOSGeom =  (AliPHOSGeometry *) fGeom ;
  
  Int_t iDigit;
  for(iDigit=0; iDigit<fMulDigit; iDigit++) {
    digit = (AliPHOSDigit *) fDigitsList[iDigit];
    Int_t relid[4] ;
    Float_t xi ;
    Float_t zi ;
    PHOSGeom->AbsToRelNumbering(digit->GetId(), relid) ;
    PHOSGeom->RelPosInModule(relid, xi, zi);
    Float_t w = TMath::Max(0.,fW0+TMath::Log(fEnergyList[iDigit]/fAmp ) ) ;
    D += w*((xi-x)*(xi-x) + (zi-z)*(zi-z) ) ; 
    wtot+=w ;
  }

  D /= wtot ;

  return TMath::Sqrt(D) ;
}

//____________________________________________________________________________
void  AliPHOSEmcRecPoint::GetElipsAxis(Float_t * lambda)
{
  Float_t wtot = 0. ;
  Float_t x    = 0.;
  Float_t z    = 0.;
  Float_t Dxx  = 0.;
  Float_t Dzz  = 0.;
  Float_t Dxz  = 0.;

  AliPHOSDigit * digit ;
  AliPHOSGeometry * PHOSGeom =  (AliPHOSGeometry *) fGeom ;
  Int_t iDigit;

  for(iDigit=0; iDigit<fMulDigit; iDigit++) {
    digit = (AliPHOSDigit *) fDigitsList[iDigit];
    Int_t relid[4] ;
    Float_t xi ;
    Float_t zi ;
    PHOSGeom->AbsToRelNumbering(digit->GetId(), relid) ;
    PHOSGeom->RelPosInModule(relid, xi, zi);
    Float_t w = TMath::Max(0.,fW0+TMath::Log(fEnergyList[iDigit]/fAmp ) ) ;
    Dxx  += w * xi * xi ;
    x    += w * xi ;
    Dzz  += w * zi * zi ;
    z    += w * zi ; 
    Dxz  += w * xi * zi ; 
    wtot += w ;
  }

  Dxx /= wtot ;
  x   /= wtot ;
  Dxx -= x * x ;
  Dzz /= wtot ;
  z   /= wtot ;
  Dzz -= z * z ;
  Dxz /= wtot ;
  Dxz -= x * z ;

  lambda[0] = TMath::Sqrt( 0.5 * (Dxx + Dzz) + TMath::Sqrt( 0.25 * (Dxx - Dzz) * (Dxx - Dzz) + Dxz * Dxz ) ) ;
  lambda[1] = TMath::Sqrt( 0.5 * (Dxx + Dzz) - TMath::Sqrt( 0.25 * (Dxx - Dzz) * (Dxx - Dzz) + Dxz * Dxz ) ) ;
}

//____________________________________________________________________________
Float_t AliPHOSEmcRecPoint::GetMaximalEnergy(void)
{
  Float_t menergy = 0. ;

  Int_t iDigit;

  for(iDigit=0; iDigit<fMulDigit; iDigit++) {
 
    if(fEnergyList[iDigit] > menergy) 
      menergy = fEnergyList[iDigit] ;
  }
  return menergy ;
}

//____________________________________________________________________________
Int_t AliPHOSEmcRecPoint::GetMultiplicityAtLevel(Float_t H) 
{
  Int_t multipl   = 0 ;
  Int_t iDigit ;
  for(iDigit=0; iDigit<fMulDigit; iDigit++) {

    if(fEnergyList[iDigit] > H * fAmp) 
      multipl++ ;
  }
  return multipl ;
}

//____________________________________________________________________________
Int_t  AliPHOSEmcRecPoint::GetNumberOfLocalMax(Int_t *  maxAt, Float_t * maxAtEnergy) 
{ 
  AliPHOSDigit * digit ;
  AliPHOSDigit * digitN ;
  

  Int_t iDigitN ;
  Int_t iDigit ;

  for(iDigit = 0; iDigit < fMulDigit; iDigit++){
    maxAt[iDigit] = fDigitsList[iDigit] ;
  }
  
  for(iDigit = 0 ; iDigit < fMulDigit; iDigit++) {   
    if(maxAt[iDigit] != -1) {
      digit = (AliPHOSDigit *) maxAt[iDigit] ;
         
      for(iDigitN = 0; iDigitN < fMulDigit; iDigitN++) {	
	digitN = (AliPHOSDigit *) fDigitsList[iDigitN] ; 
	
	if ( AreNeighbours(digit, digitN) ) {
	  if (fEnergyList[iDigit] > fEnergyList[iDigitN] ) {    
	    maxAt[iDigitN] = -1 ;
	    // but may be digit too is not local max ?
	    if(fEnergyList[iDigit] < fEnergyList[iDigitN] + fLocMaxCut) 
	      maxAt[iDigit] = -1 ;
	  }
	  else {
	    maxAt[iDigit] = -1 ;
	    // but may be digitN too is not local max ?
	    if(fEnergyList[iDigit] > fEnergyList[iDigitN] - fLocMaxCut) 
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
void AliPHOSEmcRecPoint::GetLocalPosition(TVector3 &LPos)
{
  if( fLocPos.X() < 1000000.) { // already evaluated
   LPos = fLocPos ;
   return ;
  }

  Float_t wtot = 0. ;
 
  Int_t relid[4] ;

  Float_t x = 0. ;
  Float_t z = 0. ;
  
  AliPHOSDigit * digit ;

  AliPHOSGeometry * PHOSGeom =  (AliPHOSGeometry *) fGeom ;

  Int_t iDigit;


  for(iDigit=0; iDigit<fMulDigit; iDigit++) {
    digit = (AliPHOSDigit *) fDigitsList[iDigit];

    Float_t xi ;
    Float_t zi ;
    PHOSGeom->AbsToRelNumbering(digit->GetId(), relid) ;
    PHOSGeom->RelPosInModule(relid, xi, zi);
    Float_t w = TMath::Max( 0., fW0 + TMath::Log( fEnergyList[iDigit] / fAmp ) ) ;
    x    += xi * w ;
    z    += zi * w ;
    wtot += w ;

  }

  x /= wtot ;
  z /= wtot ;
  fLocPos.SetX(x)  ;
  fLocPos.SetY(0.) ;
  fLocPos.SetZ(z)  ;

  LPos = fLocPos ;
}

// //____________________________________________________________________________
// AliPHOSEmcRecPoint& AliPHOSEmcRecPoint::operator = (AliPHOSEmcRecPoint Clu) 
// {
//   int * DL = Clu.GetDigitsList() ; 
 
//  if(fDigitsList) 
//     delete fDigitsList  ;

//   AliPHOSDigit * digit ;
 
//   Int_t iDigit;

//   for(iDigit=0; iDigit<fMulDigit; iDigit++) {
//     digit = (AliPHOSDigit *) DL[iDigit];
//     AddDigit(*digit) ;
//   }

//   fAmp       = Clu.GetTotalEnergy() ;
//   fGeom      = Clu.GetGeom() ;
//   TVector3 LocPos;
//   Clu.GetLocalPosition(LocPos) ;
//   fLocPos    = LocPos;
//   fMulDigit       = Clu.GetMultiplicity() ;
//   fMaxDigit  = Clu.GetMaximumMultiplicity() ;
//   fPHOSMod   = Clu.GetPHOSMod() ;
//   fW0        = Clu.GetLogWeightCut() ; 
//   fDelta     = Clu.GetDelta() ;
//   fLocMaxCut = Clu.GetLocMaxCut() ;
  
//   delete DL ; 
 
//   return *this ;
// }

//____________________________________________________________________________
void AliPHOSEmcRecPoint::Print(Option_t * option) 
{
  cout << "AliPHOSEmcRecPoint: " << endl ;

  AliPHOSDigit * digit ; 
  Int_t iDigit;
  AliPHOSGeometry * PHOSGeom =  (AliPHOSGeometry *) fGeom ;

  Float_t xi ;
  Float_t zi ;
  Int_t relid[4] ; 
 
  for(iDigit=0; iDigit<fMulDigit; iDigit++) {
    digit = (AliPHOSDigit *) fDigitsList[iDigit];
    PHOSGeom->AbsToRelNumbering(digit->GetId(), relid) ;
    PHOSGeom->RelPosInModule(relid, xi, zi);
    cout << " Id = " << digit->GetId() ;  
    cout << "   module  = " << relid[0] ;  
    cout << "   x  = " << xi ;  
    cout << "   z  = " << zi ;  
    cout << "   Energy = " << fEnergyList[iDigit] << endl ;
  }
  cout << "       Multiplicity    = " << fMulDigit  << endl ;
  cout << "       Cluster Energy  = " << fAmp << endl ;
  
}

//______________________________________________________________________________
void AliPHOSEmcRecPoint::Streamer(TBuffer &R__b)
{
   // Stream an object of class AliPHOSEmcRecPoint.

   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(); if (R__v) { }
      AliPHOSRecPoint::Streamer(R__b);
      R__b >> fDelta;
      R__b >> fLocMaxCut;
      R__b.ReadArray(fEnergyList);
      R__b >> fW0;
   } else {
      R__b.WriteVersion(AliPHOSEmcRecPoint::IsA());
      AliPHOSRecPoint::Streamer(R__b);
      R__b << fDelta;
      R__b << fLocMaxCut;
      R__b.WriteArray(fEnergyList, GetMaximumDigitMultiplicity() );
      R__b << fW0;
   }
}


