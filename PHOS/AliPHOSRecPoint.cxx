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
// PHOSRecPoint base class deriving from AliRecPoint
//*-- Author : Gines MARTINEZ  SUBATECH 
//////////////////////////////////////////////////////////////////////////////

// --- ROOT system ---
#include "TPad.h"

// --- Standard library ---
#include <iostream>
#include <cstdio>

// --- AliRoot header files ---

#include "AliPHOSGeometry.h"
#include "AliPHOSDigit.h"
#include "AliPHOSRecPoint.h"



ClassImp(AliPHOSRecPoint)


//____________________________________________________________________________
AliPHOSRecPoint::AliPHOSRecPoint()
  : AliRecPoint()
{
  fGeom =   AliPHOSGeometry::GetInstance() ;
  fPHOSMod = 0;
}

//____________________________________________________________________________
AliPHOSRecPoint::~AliPHOSRecPoint()
{
  // dtor
}
//____________________________________________________________________________
Int_t AliPHOSRecPoint::DistancetoPrimitive(Int_t px, Int_t py)
{
  //Compute distance from point px,py to  a AliPHOSRecPoint considered as a Tmarker
  //  Compute the closest distance of approach from point px,py to this marker.
  //  The distance is computed in pixels units.
  //

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
  //

  //  static Int_t pxold, pyold;

   static TGraph *  DigitGraph = 0 ;
   static TPaveText* ClusterText = 0 ;

   if (!gPad->IsEditable()) return;

   switch (event) {


   case kButton1Down:{
     AliPHOSDigit * digit ;
     AliPHOSGeometry * PHOSGeom =  (AliPHOSGeometry *) fGeom ;
     Int_t iDigit;
     Int_t relid[4] ;
     Float_t xi[fMulDigit] ;
     Float_t zi[fMulDigit] ;
 
     for(iDigit=0; iDigit<fMulDigit; iDigit++) {
       digit = (AliPHOSDigit *) fDigitsList[iDigit];
       PHOSGeom->AbsToRelNumbering(digit->GetId(), relid) ;
       PHOSGeom->RelPosInModule(relid, xi[iDigit], zi[iDigit]) ;
     }

     if (!DigitGraph) {
       DigitGraph = new TGraph(fMulDigit,xi,zi);
       DigitGraph-> SetMarkerStyle(5) ; 
       DigitGraph-> SetMarkerSize(1.) ;
       DigitGraph-> SetMarkerColor(1) ;
       DigitGraph-> Draw("P") ;
     }
     if (!ClusterText) {
  
       TVector3 pos(0.,0.,0.) ;
       GetLocalPosition(pos) ;
       ClusterText = new TPaveText(pos.X()-10,pos.Z()+10,pos.X()+50,pos.Z()+35,"") ;
       Text_t  line1[40] ;
       Text_t  line2[40] ;
       sprintf(line1,"Energy=%1.2f GeV",GetEnergy()) ;
       sprintf(line2,"%d Digits",GetDigitsMultiplicity()) ;
       ClusterText ->AddText(line1) ;
       ClusterText ->AddText(line2) ;
       ClusterText ->Draw("");
     }
     gPad->Update() ; 
     Print() ;
  }

     break;

   case kButton1Up:
     if (DigitGraph) {
       delete DigitGraph  ;
       DigitGraph = 0 ;
     }
     if (ClusterText) {
       delete ClusterText ;
       ClusterText = 0 ;
     }
     
     break;
     
   }
}


//____________________________________________________________________________
Int_t AliPHOSRecPoint::GetPHOSMod()
{ 
  if(fPHOSMod > 0) 
    return fPHOSMod ;

  Int_t relid[4] ;
  
  AliPHOSDigit * digit   ;
  digit = (AliPHOSDigit *) fDigitsList[0] ;
  AliPHOSGeometry * PHOSGeom =  (AliPHOSGeometry *) fGeom ;

  PHOSGeom->AbsToRelNumbering(digit->GetId(), relid) ;
  fPHOSMod = relid[0];
  return fPHOSMod ;
}

//______________________________________________________________________________
void AliPHOSRecPoint::Paint(Option_t *)
{
// Paint this ALiRecPoint as a TMarker  with its current attributes

   TVector3 pos(0.,0.,0.)  ;
   GetLocalPosition(pos)   ;
   Coord_t x = pos.X()     ;
   Coord_t y = pos.Z()     ;
   Color_t MarkerColor = 1 ;
   Size_t  MarkerSize = 1. ;
   Style_t MarkerStyle = 5 ;

   if (!gPad->IsBatch()) {
     gVirtualX->SetMarkerColor(MarkerColor) ;
     gVirtualX->SetMarkerSize (MarkerSize)  ;
     gVirtualX->SetMarkerStyle(MarkerStyle) ;
   }
   gPad->SetAttMarkerPS(MarkerColor,MarkerStyle,MarkerSize) ;
   gPad->PaintPolyMarker(1,&x,&y,"") ;
}
