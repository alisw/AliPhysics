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
//  Base Class for PHOS Reconstructed Points  
//  Why should I put meaningless comments
//  just to satisfy
//  the code checker                
//*-- Author: Gines Martinez (SUBATECH)

// --- ROOT system ---
#include "TPad.h"
#include "TClonesArray.h"

// --- Standard library ---
#include <iostream.h>
#include <stdio.h>

// --- AliRoot header files ---

#include "AliPHOSGeometry.h"
#include "AliPHOSDigit.h"
#include "AliPHOSRecPoint.h"

ClassImp(AliPHOSRecPoint)


//____________________________________________________________________________
AliPHOSRecPoint::AliPHOSRecPoint()
  : AliRecPoint()
{
  // ctor

  fGeom = (AliPHOSGeometry*) AliPHOSGeometry::GetInstance() ;
  fMaxTrack = 20 ;
  fPHOSMod = 0;
}

//____________________________________________________________________________
Int_t AliPHOSRecPoint::DistancetoPrimitive(Int_t px, Int_t py)
{
  // Compute distance from point px,py to  a AliPHOSRecPoint considered as a Tmarker
  // Compute the closest distance of approach from point px,py to this marker.
  // The distance is computed in pixels units.

  TVector3 pos(0.,0.,0.) ;
  GetLocalPosition( pos) ;
  Float_t x =  pos.X() ;
  Float_t y =  pos.Z() ;
  const Int_t kMaxDiff = 10;
  Int_t pxm  = gPad->XtoAbsPixel(x);
  Int_t pym  = gPad->YtoAbsPixel(y);
  Int_t dist = (px-pxm)*(px-pxm) + (py-pym)*(py-pym);
  
  if (dist > kMaxDiff) return 9999;
  return dist;
}

//___________________________________________________________________________
 void AliPHOSRecPoint::Draw(Option_t *option)
 {
   // Draw this AliPHOSRecPoint with its current attributes
   
   AppendPad(option);
 }

//______________________________________________________________________________
void AliPHOSRecPoint::ExecuteEvent(Int_t event, Int_t px, Int_t py)
{
  // Execute action corresponding to one event
  // This member function is called when a AliPHOSRecPoint is clicked with the locator
  //
  // If Left button is clicked on AliPHOSRecPoint, the digits are switched on    
  // and switched off when the mouse button is released.

  //  static Int_t pxold, pyold;

  static TGraph *  digitgraph = 0 ;
  static TPaveText* clustertext = 0 ;
  
  if (!gPad->IsEditable()) return;
  
  switch (event) {
    
    
  case kButton1Down:{
    AliPHOSDigit * digit ;
    AliPHOSGeometry * phosgeom =  (AliPHOSGeometry *) fGeom ;
    Int_t iDigit;
    Int_t relid[4] ;
  
    const Int_t kMulDigit=AliPHOSRecPoint::GetDigitsMultiplicity() ;
    Float_t * xi = new Float_t [kMulDigit] ; 
    Float_t * zi = new Float_t [kMulDigit] ;
    
    for(iDigit = 0; iDigit < kMulDigit; iDigit++) {
      digit = (AliPHOSDigit *) fDigitsList[iDigit];
      phosgeom->AbsToRelNumbering(digit->GetId(), relid) ;
      phosgeom->RelPosInModule(relid, xi[iDigit], zi[iDigit]) ;
    }
    
    if (!digitgraph) {
      digitgraph = new TGraph(fMulDigit,xi,zi);
      digitgraph-> SetMarkerStyle(5) ; 
      digitgraph-> SetMarkerSize(1.) ;
      digitgraph-> SetMarkerColor(1) ;
      digitgraph-> Draw("P") ;
    }
    if (!clustertext) {
      
      TVector3 pos(0.,0.,0.) ;
      GetLocalPosition(pos) ;
      clustertext = new TPaveText(pos.X()-10,pos.Z()+10,pos.X()+50,pos.Z()+35,"") ;
      Text_t  line1[40] ;
      Text_t  line2[40] ;
      sprintf(line1,"Energy=%1.2f GeV",GetEnergy()) ;
      sprintf(line2,"%d Digits",GetDigitsMultiplicity()) ;
      clustertext ->AddText(line1) ;
      clustertext ->AddText(line2) ;
      clustertext ->Draw("");
    }
    gPad->Update() ; 
    Print() ;
    delete[] xi ; 
    delete[] zi ; 
   }
  
break;
  
  case kButton1Up:
    if (digitgraph) {
      delete digitgraph  ;
      digitgraph = 0 ;
    }
    if (clustertext) {
      delete clustertext ;
      clustertext = 0 ;
    }
    
    break;
    
  }
}
//____________________________________________________________________________
void AliPHOSRecPoint::EvalAll(Float_t logWeight,TClonesArray * digits) {
  //evaluates (if necessary) all RecPoint data members 

  EvalPrimaries(digits) ;
}
//____________________________________________________________________________
void AliPHOSRecPoint::EvalPHOSMod(AliPHOSDigit * digit) 
{
  // Returns the PHOS module in which the RecPoint is found

  if( fPHOSMod == 0){
  Int_t relid[4] ; 
  
  AliPHOSGeometry * phosgeom =  (AliPHOSGeometry *) fGeom ;

  phosgeom->AbsToRelNumbering(digit->GetId(), relid) ;
  fPHOSMod = relid[0];
  }
}

//______________________________________________________________________________
void  AliPHOSRecPoint::EvalPrimaries(TClonesArray * digits)
{
  // Constructs the list of primary particles (tracks) which have contributed to this RecPoint
  
  AliPHOSDigit * digit ;
  Int_t * tempo    = new Int_t[fMaxTrack] ;

  Int_t index ;  
  for ( index = 0 ; index < GetDigitsMultiplicity() ; index++ ) { // all digits
    digit = (AliPHOSDigit *) digits->At( fDigitsList[index] ) ; 
    Int_t nprimaries = digit->GetNprimary() ;
    Int_t * newprimaryarray = new Int_t[nprimaries] ;
    Int_t ii ; 
    for ( ii = 0 ; ii < nprimaries ; ii++)
      newprimaryarray[ii] = digit->GetPrimary(ii+1) ; 

    Int_t jndex ;
    for ( jndex = 0 ; jndex < nprimaries ; jndex++ ) { // all primaries in digit
      if ( fMulTrack > fMaxTrack ) {
	fMulTrack = - 1 ;
	cout << "AliPHOSRecPoint::GetNprimaries ERROR > increase fMaxTrack " << endl ;
	break ;
      }
      Int_t newprimary = newprimaryarray[jndex] ;
      Int_t kndex ;
      Bool_t already = kFALSE ;
      for ( kndex = 0 ; kndex < fMulTrack ; kndex++ ) { //check if not already stored
	if ( newprimary == tempo[kndex] ){
	  already = kTRUE ;
	  break ;
	}
      } // end of check
      if ( !already) { // store it
	tempo[fMulTrack] = newprimary ; 
	fMulTrack++ ;
      } // store it
    } // all primaries in digit
    delete newprimaryarray ; 
  } // all digits

  
  fTracksList = new Int_t[fMulTrack] ;
  for(index = 0; index < fMulTrack; index++)
   fTracksList[index] = tempo[index] ;
 
  delete tempo ;

}

//______________________________________________________________________________
void AliPHOSRecPoint::Paint(Option_t *)
{
  // Paint this ALiRecPoint as a TMarker  with its current attributes
  
  TVector3 pos(0.,0.,0.)  ;
  GetLocalPosition(pos)   ;
  Coord_t x = pos.X()     ;
  Coord_t y = pos.Z()     ;
  Color_t markercolor = 1 ;
  Size_t  markersize = 1. ;
  Style_t markerstyle = 5 ;
  
  if (!gPad->IsBatch()) {
    gVirtualX->SetMarkerColor(markercolor) ;
    gVirtualX->SetMarkerSize (markersize)  ;
    gVirtualX->SetMarkerStyle(markerstyle) ;
  }
  gPad->SetAttMarkerPS(markercolor,markerstyle,markersize) ;
  gPad->PaintPolyMarker(1,&x,&y,"") ;
}
