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
//  RecPoint implementation for PHOS-CPV
//  An CpvRecPoint is a cluster of digits   
//*-- Author: Yuri Kharlov
//  (after Dmitri Peressounko (RRC KI & SUBATECH))
//  30 October 2000 

// --- ROOT system ---
#include "TPad.h"
#include "TH2.h"
#include "TMath.h" 
#include "TCanvas.h" 

// --- Standard library ---

#include <iostream.h> 

// --- AliRoot header files ---

#include "AliPHOSGeometry.h"
#include "AliPHOSCpvRecPoint.h"
#include "AliPHOSPpsdRecPoint.h"
#include "AliRun.h"
#include "AliPHOSIndexToObject.h"

ClassImp(AliPHOSCpvRecPoint)

//____________________________________________________________________________
AliPHOSCpvRecPoint::AliPHOSCpvRecPoint(Float_t W0, Float_t LocMaxCut)
  : AliPHOSRecPoint()
{
  // ctor

  fMulDigit   = 0 ;  
  fAmp   = 0. ;   
  fEnergyList = new Float_t[fMaxDigit]; 
  AliPHOSGeometry * phosgeom =  (AliPHOSGeometry *) fGeom ;
  fDelta     =  phosgeom->GetCrystalSize(0) ; 
  fW0        = W0 ;          
  fLocMaxCut = LocMaxCut ; 
  fLocPos.SetX(1000000.)  ;      //Local position should be evaluated
  fLengX = -1;
  fLengZ = -1;
}

//____________________________________________________________________________
AliPHOSCpvRecPoint::~AliPHOSCpvRecPoint()
{
  // dtor
  if ( fEnergyList ) delete[] fEnergyList ; 
}

//____________________________________________________________________________
void AliPHOSCpvRecPoint::AddDigit(AliPHOSDigit & digit, Float_t Energy)
{
  // Adds a digit to the RecPoint
  //  and accumulates the total amplitude and the multiplicity 
  
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
}

//____________________________________________________________________________
Bool_t AliPHOSCpvRecPoint::AreNeighbours(AliPHOSDigit * digit1, AliPHOSDigit * digit2 ) 
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
Int_t AliPHOSCpvRecPoint::Compare(TObject * obj)
{
  // Compares two RecPoints according to their position in the PHOS modules

  Int_t rv ; 

  if( (strcmp(obj->ClassName() , "AliPHOSPpsdRecPoint" )) == 0)  // PPSD Rec Point
    {
      AliPHOSPpsdRecPoint * clu = (AliPHOSPpsdRecPoint *)obj ; 
      if(this->GetPHOSMod()  < clu->GetPHOSMod() ) 
	rv = -1 ;
      else 
	rv = 1 ;
      return rv ;
    }
  else
    {
      AliPHOSCpvRecPoint * clu  = (AliPHOSCpvRecPoint *) obj ; 
      
      Int_t phosmod1 = this->GetPHOSMod() ;
      Int_t phosmod2 = clu->GetPHOSMod() ;
      
      TVector3 locpos1; 
      this->GetLocalPosition(locpos1) ;
      TVector3 locpos2;  
      clu->GetLocalPosition(locpos2) ;  
      
      if(phosmod1 == phosmod2 ) {
	Int_t rowdif = (Int_t)TMath::Ceil(locpos1.X()/fDelta)-(Int_t)TMath::Ceil(locpos2.X()/fDelta) ;
	if (rowdif> 0) 
	  rv = -1 ;
	else if(rowdif < 0) 
	  rv = 1 ;
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
}

//______________________________________________________________________________
void AliPHOSCpvRecPoint::ExecuteEvent(Int_t event, Int_t px, Int_t py)
{
  // Execute action corresponding to one event
  //  This member function is called when a AliPHOSRecPoint is clicked with the locator
  //
  //  If Left button is clicked on AliPHOSRecPoint, the digits are switched on    
  //  and switched off when the mouse button is released.
  //

  //   static Int_t pxold, pyold;

  AliPHOSIndexToObject * please =  AliPHOSIndexToObject::GetInstance() ; 
  
  static TGraph *  digitgraph = 0 ;
  
  if (!gPad->IsEditable()) return;
  
  TH2F * histo = 0 ;
  TCanvas * histocanvas ; 
  
  switch (event) {
    
  case kButton1Down: {
    AliPHOSDigit * digit ;
    AliPHOSGeometry * phosgeom =  (AliPHOSGeometry *) fGeom ;
    Int_t iDigit;
    Int_t relid[4] ;
    
    const Int_t kMulDigit = AliPHOSCpvRecPoint::GetDigitsMultiplicity() ; 
    Float_t * xi = new Float_t[kMulDigit] ; 
    Float_t * zi = new Float_t[kMulDigit] ; 
    
    // create the histogram for the single cluster 
    // 1. gets histogram boundaries
    Float_t ximax = -999. ; 
    Float_t zimax = -999. ; 
    Float_t ximin = 999. ; 
    Float_t zimin = 999. ;
    
    for(iDigit=0; iDigit<kMulDigit; iDigit++) {
      digit = (AliPHOSDigit *) ( please->GimeDigit(fDigitsList[iDigit]) ) ;
      phosgeom->AbsToRelNumbering(digit->GetId(), relid) ;
      phosgeom->RelPosInModule(relid, xi[iDigit], zi[iDigit]);
      if ( xi[iDigit] > ximax )
	ximax = xi[iDigit] ; 
      if ( xi[iDigit] < ximin )
	ximin = xi[iDigit] ; 
      if ( zi[iDigit] > zimax )
	zimax = zi[iDigit] ; 
      if ( zi[iDigit] < zimin )
	zimin = zi[iDigit] ;     
    }
    ximax += phosgeom->GetCrystalSize(0) / 2. ;
    zimax += phosgeom->GetCrystalSize(2) / 2. ;
    ximin -= phosgeom->GetCrystalSize(0) / 2. ;
    zimin -= phosgeom->GetCrystalSize(2) / 2. ;
    Int_t xdim = (int)( (ximax - ximin ) / phosgeom->GetCrystalSize(0) + 0.5  ) ; 
    Int_t zdim = (int)( (zimax - zimin ) / phosgeom->GetCrystalSize(2) + 0.5 ) ;
    
    // 2. gets the histogram title
    
    Text_t title[100] ; 
    sprintf(title,"Energy=%1.2f GeV ; Digits ; %d ", GetEnergy(), GetDigitsMultiplicity()) ;
    
    if (!histo) {
      delete histo ; 
      histo = 0 ; 
    }
    histo = new TH2F("cluster3D", title,  xdim, ximin, ximax, zdim, zimin, zimax)  ;
    
    Float_t x, z ; 
    for(iDigit=0; iDigit<kMulDigit; iDigit++) {
      digit = (AliPHOSDigit *) ( please->GimeDigit(fDigitsList[iDigit]) ) ;
      phosgeom->AbsToRelNumbering(digit->GetId(), relid) ;
      phosgeom->RelPosInModule(relid, x, z);
      histo->Fill(x, z, fEnergyList[iDigit] ) ;
    }
    
    if (!digitgraph) {
      digitgraph = new TGraph(kMulDigit,xi,zi);
      digitgraph-> SetMarkerStyle(5) ; 
      digitgraph-> SetMarkerSize(1.) ;
      digitgraph-> SetMarkerColor(1) ;
      digitgraph-> Paint("P") ;
    }
    
    Print() ;
    histocanvas = new TCanvas("cluser", "a single cluster", 600, 500) ; 
    histocanvas->Draw() ; 
    histo->Draw("lego1") ; 
    
    delete[] xi ; 
    delete[] zi ; 
    
    break;
  }
  
  case kButton1Up: 
    if (digitgraph) {
      delete digitgraph  ;
      digitgraph = 0 ;
    }
    break;
  
   }
}

//____________________________________________________________________________
Float_t  AliPHOSCpvRecPoint::GetDispersion() 
{
  // Calculates the dispersion of the shower at the origine of the RecPoint

  AliPHOSIndexToObject * please =  AliPHOSIndexToObject::GetInstance() ; 

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
    digit = (AliPHOSDigit *) ( please->GimeDigit(fDigitsList[iDigit]) ) ;
    Int_t relid[4] ;
    Float_t xi ;
    Float_t zi ;
    phosgeom->AbsToRelNumbering(digit->GetId(), relid) ;
    phosgeom->RelPosInModule(relid, xi, zi);
    Float_t w = TMath::Max(0.,fW0+TMath::Log(fEnergyList[iDigit]/fAmp ) ) ;
    d += w*((xi-x)*(xi-x) + (zi-z)*(zi-z) ) ; 
    wtot+=w ;
  }

  d /= wtot ;

  return TMath::Sqrt(d) ;
}

//____________________________________________________________________________
void  AliPHOSCpvRecPoint::GetElipsAxis(Float_t * lambda)
{
  // Calculates the axis of the shower ellipsoid

  AliPHOSIndexToObject * please =  AliPHOSIndexToObject::GetInstance() ; 

  Float_t wtot = 0. ;
  Float_t x    = 0.;
  Float_t z    = 0.;
  Float_t dxx  = 0.;
  Float_t dzz  = 0.;
  Float_t dxz  = 0.;

  AliPHOSDigit * digit ;
  AliPHOSGeometry * phosgeom =  (AliPHOSGeometry *) fGeom ;
  Int_t iDigit;

  for(iDigit=0; iDigit<fMulDigit; iDigit++) {
    digit = (AliPHOSDigit *) ( please->GimeDigit(fDigitsList[iDigit]) ) ;
    Int_t relid[4] ;
    Float_t xi ;
    Float_t zi ;
    phosgeom->AbsToRelNumbering(digit->GetId(), relid) ;
    phosgeom->RelPosInModule(relid, xi, zi);
    Float_t w = TMath::Max(0.,fW0+TMath::Log(fEnergyList[iDigit]/fAmp ) ) ;
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

  Float_t tmp0 = 0.5 * (dxx + dzz) + TMath::Sqrt( 0.25 * (dxx - dzz) * (dxx - dzz) + dxz * dxz );
  Float_t tmp1 = 0.5 * (dxx + dzz) - TMath::Sqrt( 0.25 * (dxx - dzz) * (dxx - dzz) + dxz * dxz );
  if (tmp0>=0) lambda[0] = TMath::Sqrt(tmp0);
  else lambda[0] = -TMath::Sqrt(-tmp0);
  if (tmp1>=0) lambda[1] = TMath::Sqrt(tmp1);
  else lambda[1] = -TMath::Sqrt(-tmp1);
//    lambda[0] = TMath::Sqrt( 0.5 * (dxx + dzz) + TMath::Sqrt( 0.25 * (dxx - dzz) * (dxx - dzz) + dxz * dxz ) ) ;
//    lambda[1] = TMath::Sqrt( 0.5 * (dxx + dzz) - TMath::Sqrt( 0.25 * (dxx - dzz) * (dxx - dzz) + dxz * dxz ) ) ;
}

//____________________________________________________________________________
void AliPHOSCpvRecPoint::GetLocalPosition(TVector3 &LPos)
{
  // Calculates the center of gravity in the local PHOS-module coordinates 

  AliPHOSIndexToObject * please =  AliPHOSIndexToObject::GetInstance() ; 

  if( fLocPos.X() < 1000000.) { // already evaluated
   LPos = fLocPos ;
   return ;
  }

  Float_t wtot = 0. ;
 
  Int_t relid[4] ;

  Float_t x = 0. ;
  Float_t z = 0. ;
  
  AliPHOSDigit * digit ;

  AliPHOSGeometry * phosgeom =  (AliPHOSGeometry *) fGeom ;

  Int_t iDigit;

  for(iDigit=0; iDigit<fMulDigit; iDigit++) {
    digit = (AliPHOSDigit *) ( please->GimeDigit(fDigitsList[iDigit]) );

    Float_t xi ;
    Float_t zi ;
    phosgeom->AbsToRelNumbering(digit->GetId(), relid) ;
    phosgeom->RelPosInModule(relid, xi, zi);
    Float_t w = TMath::Max( 0., fW0 + TMath::Log( fEnergyList[iDigit] / fAmp ) ) ;
    x    += xi * w ;
    z    += zi * w ;
    wtot += w ;
  }

  if (wtot != 0) {
    x /= wtot ;
    z /= wtot ;
  } else {
    x = -1e6 ;
    z = -1e6 ;
    if (fMulDigit != 0) cout << "AliPHOSCpvRecPoint: too low log weight factor "
			     << "to evaluate cluster's center\n";
  }
  fLocPos.SetX(x)  ;
  fLocPos.SetY(0.) ;
  fLocPos.SetZ(z)  ;

  LPos = fLocPos ;
}

//____________________________________________________________________________
void AliPHOSCpvRecPoint::GetClusterLengths(Int_t &lengX, Int_t &lengZ)
{
  // Calculates the cluster lengths along X and Z axes
  // These characteristics are needed for CPV to tune
  // digitization+reconstruction to experimental data
  // Yuri Kharlov. 24 October 2000

  if( fLengX != -1 && fLengZ != -1 ) { // already evaluated
    lengX = fLengX ;
    lengZ = fLengZ ;
    return ;
  }

  AliPHOSIndexToObject * please =  AliPHOSIndexToObject::GetInstance() ; 

  Int_t relid[4] ;

  AliPHOSDigit * digit ;

  AliPHOSGeometry * phosgeom =  (AliPHOSGeometry *) fGeom ;

  const Int_t kMaxLeng=20;
  Int_t idX[kMaxLeng], idZ[kMaxLeng];
  lengX = 0;
  lengZ = 0;
  Bool_t dejavu;

  for(Int_t iDigit=0; iDigit<fMulDigit; iDigit++) {
    digit = (AliPHOSDigit *) ( please->GimeDigit(fDigitsList[iDigit]) );
    Int_t absId = digit->GetId();
    phosgeom->AbsToRelNumbering(absId, relid) ;

    Int_t i;
    dejavu=kFALSE;
    for (i=0; i<lengX; i++) if (relid[2]==idX[i]) { dejavu=kTRUE; break; }
    if (!dejavu) {
      idX[lengX]=relid[2];
      lengX++;
      lengX = TMath::Min(lengX,kMaxLeng);
    }

    dejavu=kFALSE;
    for (i=0; i<lengZ; i++) if (relid[3]==idZ[i]) { dejavu=kTRUE; break; }
    if (!dejavu) {
      idZ[lengZ]=relid[3];
      lengZ++;
      lengZ = TMath::Min(lengZ,kMaxLeng);
    }
  }
  fLengX = lengX;
  fLengZ = lengZ;
}

//____________________________________________________________________________
Float_t AliPHOSCpvRecPoint::GetMaximalEnergy(void)
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
Int_t AliPHOSCpvRecPoint::GetMultiplicityAtLevel(const Float_t H) 
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
Int_t  AliPHOSCpvRecPoint::GetNumberOfLocalMax(Int_t *  maxAt, Float_t * maxAtEnergy) 
{ 
  // Calculates the number of local maxima in the cluster using fLocalMaxCut as the minimum
  //  energy difference between two local maxima

  AliPHOSIndexToObject * please =  AliPHOSIndexToObject::GetInstance() ; 

  AliPHOSDigit * digit ;
  AliPHOSDigit * digitN ;
  

  Int_t iDigitN ;
  Int_t iDigit ;

  for(iDigit = 0; iDigit < fMulDigit; iDigit++){
    maxAt[iDigit] = (Int_t) ( please->GimeDigit(fDigitsList[iDigit]) ) ;
  }
  
  for(iDigit = 0 ; iDigit < fMulDigit; iDigit++) {   
    if(maxAt[iDigit] != -1) {
      digit = (AliPHOSDigit *) maxAt[iDigit] ;
          
      for(iDigitN = 0; iDigitN < fMulDigit; iDigitN++) {	
	digitN = (AliPHOSDigit *) ( please->GimeDigit(fDigitsList[iDigitN]) ) ; 
	
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
void AliPHOSCpvRecPoint::Print(Option_t * option) 
{
  // Print the list of digits belonging to the cluster
  
  cout << "AliPHOSCpvRecPoint: " << endl ;

  AliPHOSDigit * digit ; 
  Int_t iDigit;
  AliPHOSGeometry * phosgeom =  (AliPHOSGeometry *) fGeom ;

  Float_t xi ;
  Float_t zi ;
  Int_t relid[4] ; 
  
  AliPHOSIndexToObject * please =  AliPHOSIndexToObject::GetInstance() ; 

  for(iDigit=0; iDigit<fMulDigit; iDigit++) {
    digit = please->GimeDigit( fDigitsList[iDigit] ) ; 
    phosgeom->AbsToRelNumbering(digit->GetId(), relid) ;
    phosgeom->RelPosInModule(relid, xi, zi);
    cout << " Id = " << digit->GetId() ;  
    cout << "   module  = " << relid[0] ;  
    cout << "   x  = " << xi ;  
    cout << "   z  = " << zi ;  
    cout << "   Energy = " << fEnergyList[iDigit] << endl ;
  }
  cout << "       Multiplicity    = " << fMulDigit  << endl ;
  cout << "       Cluster Energy  = " << fAmp << endl ;
  cout << "       Stored at position " << GetIndexInList() << endl ; 
 
}
