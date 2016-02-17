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

/* History of cvs commits:
 *
 * $Log$
 * Revision 1.26  2007/06/18 07:02:44  kharlov
 * Change the signature of EvalLocalPosition() to obey the method virtuality from the parent class
 *
 * Revision 1.25  2007/03/06 06:47:28  kharlov
 * DP:Possibility to use actual vertex position added
 */

//_________________________________________________________________________
//  RecPoint implementation for PHOS-CPV
//  An CpvRecPoint is a cluster of digits   
//-- Author: Yuri Kharlov
//  (after Dmitri Peressounko (RRC KI & SUBATECH))
//  30 October 2000 

// --- ROOT system ---

#include "TMath.h" 
#include "TClonesArray.h" 

// --- Standard library ---

// --- AliRoot header files ---
#include "AliLog.h"
#include "AliPHOSGeometry.h" 
#include "AliPHOSDigit.h"
#include "AliPHOSCpvRecPoint.h"
#include "AliPHOSLoader.h"

ClassImp(AliPHOSCpvRecPoint)

//____________________________________________________________________________
AliPHOSCpvRecPoint::AliPHOSCpvRecPoint() : 
  AliPHOSEmcRecPoint(),
  fLengX(-1),
  fLengZ(-1)
{
  // ctor
}

//____________________________________________________________________________
AliPHOSCpvRecPoint::AliPHOSCpvRecPoint(const char * opt) : 
  AliPHOSEmcRecPoint(opt),
  fLengX(-1),
  fLengZ(-1)
{
   // ctor
}

//____________________________________________________________________________
AliPHOSCpvRecPoint::~AliPHOSCpvRecPoint()
{
  // dtor
}

//____________________________________________________________________________
Bool_t AliPHOSCpvRecPoint::AreNeighbours(AliPHOSDigit * digit1, AliPHOSDigit * digit2 ) const
{
  // Tells if (true) or not (false) two digits are neighbors)
  
  Bool_t aren = kFALSE ;
  
  AliPHOSGeometry * phosgeom =  AliPHOSGeometry::GetInstance() ;

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
Int_t AliPHOSCpvRecPoint::Compare(const TObject * obj) const
{
  // Compares two RecPoints according to their position in the PHOS modules

  Float_t delta = 1 ; //Width of "Sorting row". If you changibg this 
                      //value (what is senseless) change as vell delta in
                      //AliPHOSTrackSegmentMakerv* and other RecPoints...

  Int_t rv ; 

  AliPHOSCpvRecPoint * clu  = (AliPHOSCpvRecPoint *) obj ; 
  
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
void AliPHOSCpvRecPoint::ExecuteEvent(Int_t, Int_t, Int_t ) /*const*/
{
//   // Execute action corresponding to one event
//   //  This member function is called when a AliPHOSRecPoint is clicked with the locator
//   //
//   //  If Left button is clicked on AliPHOSRecPoint, the digits are switched on    
//   //  and switched off when the mouse button is released.
//   //

//   //   static Int_t pxold, pyold;

//   AliPHOSLoader * gime =  AliPHOSLoader::GetInstance() ; 
  
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
    
//     const Int_t kMulDigit = AliPHOSCpvRecPoint::GetDigitsMultiplicity() ; 
//     Float_t * xi = new Float_t[kMulDigit] ; 
//     Float_t * zi = new Float_t[kMulDigit] ; 
    
//     // create the histogram for the single cluster 
//     // 1. gets histogram boundaries
//     Float_t ximax = -999. ; 
//     Float_t zimax = -999. ; 
//     Float_t ximin = 999. ; 
//     Float_t zimin = 999. ;
    
//     for(iDigit=0; iDigit<kMulDigit; iDigit++) {
//       digit = (AliPHOSDigit *) ( gime->Digit(fDigitsList[iDigit]) ) ;
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
//       digit = (AliPHOSDigit *) ( gime->Digit(fDigitsList[iDigit]) ) ;
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
//     histocanvas = new TCanvas("cluser", "a single cluster", 600, 500) ; 
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
void AliPHOSCpvRecPoint::EvalAll(TClonesArray * digits)
{
  // Evaluate local coordinate assuming the vertex in (000) and no inclination
  AliPHOSEmcRecPoint::EvalAll(digits) ;
}
//____________________________________________________________________________
void AliPHOSCpvRecPoint::EvalAll(Float_t logWeight, TVector3 &vtx, TClonesArray * digits)
{
  // wraps other methods
  TVector3 vInc(0,1,0);
  AliPHOSEmcRecPoint::EvalAll(logWeight,vtx,digits) ;
  EvalLocalPosition(logWeight, vtx, digits,vInc) ;
  EvalClusterLengths(digits) ;
}
//____________________________________________________________________________
void AliPHOSCpvRecPoint::EvalLocalPosition(Float_t logWeight, TVector3 & /*vtx */, TClonesArray * digits, TVector3 &/* vInc */)
{
  // Calculates the center of gravity in the local PHOS-module coordinates 

  Float_t wtot = 0. ;
 
  Int_t relid[4] ;

  Float_t x = 0. ;
  Float_t z = 0. ;
  
  AliPHOSDigit * digit ;

  AliPHOSGeometry * phosgeom =  AliPHOSGeometry::GetInstance();
  
  Int_t iDigit;

  for(iDigit=0; iDigit<fMulDigit; iDigit++) {
    digit = (AliPHOSDigit *) digits->At(fDigitsList[iDigit]);

    Float_t xi ;
    Float_t zi ;
    phosgeom->AbsToRelNumbering(digit->GetId(), relid) ;
    phosgeom->RelPosInModule(relid, xi, zi);
    if (fAmp>0 && fEnergyList[iDigit]>0) {
      Float_t w = TMath::Max( 0., logWeight + TMath::Log( fEnergyList[iDigit] / fAmp ) ) ;
      x    += xi * w ;
      z    += zi * w ;
      wtot += w ;
    }
//    else
//      AliError(Form("Wrong energy %f and/or amplitude %f\n", fEnergyList[iDigit], fAmp));
  }

  if (wtot != 0) {
    x /= wtot ;
    z /= wtot ;
  } else {
    x = 999 ;
    z = 999 ;
//    if (fMulDigit != 0) 
//      AliWarning(Form("Too low log weight factor to evaluate cluster's center" )) ;
  }
  fLocPos.SetX(x)  ;
  fLocPos.SetY(0.) ;
  fLocPos.SetZ(z)  ;

}

//____________________________________________________________________________
void AliPHOSCpvRecPoint::EvalClusterLengths(TClonesArray * digits)
{
  //Modified 15.03.2001 by Dmitri Peressounko
  
  // Calculates the cluster lengths along X and Z axes
  // These characteristics are needed for CPV to tune
  // digitization+reconstruction to experimental data
  // Yuri Kharlov. 24 October 2000

  Int_t relid[4] ;

  AliPHOSDigit * digit ;

  AliPHOSGeometry * phosgeom =  AliPHOSGeometry::GetInstance();

  const Int_t kMaxLeng=20;
  Int_t idX[kMaxLeng], idZ[kMaxLeng];
  fLengX = 0;
  fLengZ = 0;
  Bool_t dejavu;

  for(Int_t iDigit=0; iDigit<fMulDigit; iDigit++) {
    digit = (AliPHOSDigit *) digits->At(fDigitsList[iDigit]) ;
    Int_t absId = digit->GetId();
    phosgeom->AbsToRelNumbering(absId, relid) ;

    Int_t i;
    dejavu=kFALSE;
    for (i=0; i<fLengX; i++) if (relid[2]==idX[i]) { dejavu=kTRUE; break; }
    if (!dejavu) {
      idX[fLengX]=relid[2];
      fLengX++;
      fLengX = TMath::Min(fLengX,kMaxLeng);
    }

    dejavu=kFALSE;
    for (i=0; i<fLengZ; i++) if (relid[3]==idZ[i]) { dejavu=kTRUE; break; }
    if (!dejavu) {
      idZ[fLengZ]=relid[3];
      fLengZ++;
      fLengZ = TMath::Min(fLengZ,kMaxLeng);
    }
  }
}

//____________________________________________________________________________
void AliPHOSCpvRecPoint::Print(const Option_t *) const
{
  // Print the list of digits belonging to the cluster
  
  TString message ; 
  message  =  "AliPHOSCpvRecPoint: " ;
  AliInfo(message.Data()) ; 
  
  Int_t iDigit;

  for(iDigit=0; iDigit<fMulDigit; iDigit++) 
    printf("\tDigit %d: id=%d, A=%f\n", iDigit, fDigitsList[iDigit], fEnergyList[iDigit]) ; 

  message  = "\tPad multiplicity    = %d\n" ;
  message += "\tCluster amplitude  = %f\n" ;
  message += "\tStored at position %d\n" ; 
 
  printf(message.Data(), fMulDigit, fAmp, GetIndexInList() ) ; 

}
