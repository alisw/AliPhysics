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
// class of PHOS Sub Track
//*-- Author : Dmitri Peressounko RRC KI 
//////////////////////////////////////////////////////////////////////////////

// --- ROOT system ---
 
#include "TVector3.h"
#include "TPad.h"

// --- Standard library ---

#include <iostream>

// --- AliRoot header files ---

#include "AliPHOSTrackSegment.h" 
#include "AliPHOSv0.h"

ClassImp(AliPHOSTrackSegment)

//____________________________________________________________________________
AliPHOSTrackSegment::AliPHOSTrackSegment( AliPHOSEmcRecPoint * emc , AliPHOSPpsdRecPoint * ppsdRP1,
                  AliPHOSPpsdRecPoint * ppsdRP2  ) 
{     
  if( emc )   
    fEmcRecPoint =  emc ;

  if( ppsdRP1 )  
    fPpsdUp = ppsdRP1 ;

  if( ppsdRP2  ) 
    fPpsdLow = ppsdRP2 ;

  fCutOnDispersion = 1.5 ; 
}

//____________________________________________________________________________
AliPHOSTrackSegment::AliPHOSTrackSegment( const AliPHOSTrackSegment & ts) 
{
 ( (AliPHOSTrackSegment &)ts ).Copy(*this) ; 
}

//____________________________________________________________________________
AliPHOSTrackSegment::~AliPHOSTrackSegment() // dtor
{
//    fEmcRecPoint.Delete() ;   Not Owners !!!
//    fPpsdUp.Delete() ;
//    fPpsdLow.Delete() ;
}

//____________________________________________________________________________
void AliPHOSTrackSegment::Copy(TObject & obj) 
{
   TObject::Copy(obj) ;
   ( (AliPHOSTrackSegment &)obj ).fEmcRecPoint = fEmcRecPoint ; 
   ( (AliPHOSTrackSegment &)obj ).fPpsdLow     = fPpsdLow ; 
   ( (AliPHOSTrackSegment &)obj ).fPpsdUp      = fPpsdUp ; 
}
//____________________________________________________________________________
Int_t AliPHOSTrackSegment::DistancetoPrimitive(Int_t px, Int_t py)
{
  //Compute distance from point px,py to  a AliPHOSTrackSegment considered as a Tmarker
  //  Compute the closest distance of approach from point px,py to this marker.
  //  The distance is computed in pixels units.
  //
  Int_t div = 1 ;  
  TVector3 pos(0.,0.,0.) ;
 
  fEmcRecPoint->GetLocalPosition( pos) ;
  Float_t x =  pos.X() ;
  Float_t y =  pos.Z() ;
  if ( fPpsdLow ) {
    fPpsdLow->GetLocalPosition( pos ) ;
    x +=  pos.X() ;
    y +=  pos.Z() ;
    div++ ; 
  }
  if ( fPpsdUp ) {
    fPpsdUp->GetLocalPosition( pos ) ;
    x +=  pos.X() ;
    y +=  pos.Z() ;
    div++ ; 
  }
  x /= div ; 
  y /= div ; 

   const Int_t kMaxDiff = 10;
   Int_t pxm  = gPad->XtoAbsPixel(x);
   Int_t pym  = gPad->YtoAbsPixel(y);
   Int_t dist = (px-pxm)*(px-pxm) + (py-pym)*(py-pym);

   if (dist > kMaxDiff) return 9999;
   return dist;
}

//___________________________________________________________________________
 void AliPHOSTrackSegment::Draw(Option_t *option)
 {
   // Draw this AliPHOSTrackSegment with its current attribute

   AppendPad(option);
 }

//______________________________________________________________________________
void AliPHOSTrackSegment::ExecuteEvent(Int_t event, Int_t px, Int_t py)
{
  // Execute action corresponding to one event
  //  This member function is called when a AliPHOSTrackSegment is clicked with the locator
  //
  //  If Left button is clicked on AliPHOSRecPoint, the digits are switched on    
  //  and switched off when the mouse button is released.
  //
   static TPaveText* TrackSegmentText = 0 ;

   if (!gPad->IsEditable()) 
     return;

   switch (event) {

   case kButton1Down:{
    
     if (!TrackSegmentText) {
  
       TVector3 pos(0.,0.,0.) ;
       fEmcRecPoint->GetLocalPosition(pos) ;
       TrackSegmentText = new TPaveText(pos.X()-10,pos.Z()+10,pos.X()+5,pos.Z()+15,"") ;
       Text_t  line1[40] ;
       if (GetPartType() == GAMMA ) 
	 sprintf(line1,"PHOTON") ;
       if (GetPartType() == ELECTRON ) 
	 sprintf(line1,"ELECTRON") ;
       if (GetPartType() == NEUTRAL ) 
	 sprintf(line1,"NEUTRAL") ;
       if (GetPartType() == CHARGEDHADRON ) 
	 sprintf(line1,"CHARGED HADRON") ;

       TrackSegmentText ->AddText(line1) ;
       TrackSegmentText ->Draw("");
       gPad->Update() ; 
     }
   }

     break;

   case kButton1Up:
     if (TrackSegmentText) {
       delete TrackSegmentText ;
       TrackSegmentText = 0 ;
     }
     break;  
   }
}


//____________________________________________________________________________
Float_t AliPHOSTrackSegment::GetDistanceInPHOSPlane()
{
 
  TVector3 vecEmc ;
  fEmcRecPoint->GetLocalPosition(vecEmc) ;

  TVector3 vecPpsd ;
  if( fPpsdLow->GetMultiplicity() )  
    fPpsdLow->GetLocalPosition(vecPpsd)  ; 
  else { 
    vecPpsd.SetX(10000.) ;
  } 
  vecEmc -= vecPpsd ;

  Float_t R = vecEmc.Mag();;

  return R ;
}

//____________________________________________________________________________
TVector3 AliPHOSTrackSegment::GetMomentumDirection() 
{   
  TVector3 dir, tempo ; 
  TMatrix mdummy ;

  TVector3 posEmc ;
  fEmcRecPoint->GetGlobalPosition(posEmc, mdummy) ;

  TVector3 posPpsdl ;
  TVector3 posPpsdup ;
 
//   if( fPpsdLow ){
//     fPpsdLow->GetGlobalPosition(posPpsdl, mdummy) ; 
//     if( !fPpsdUp ) { // draw line trough 2 points
//       tempo = posEmc - posPpsdl;
//      }

//     else { // draw line through 3 points
//       fPpsdUp->GetGlobalPosition(posPpsdup, mdummy) ;
//       posPpsdl = 0.5 * ( posPpsdup + posPpsdl ) ; 
//       dir = posEmc - posPpsdl ;
//     }
//   }
//   else 
    tempo = posEmc ; 
    
  dir.SetX( tempo.X() ) ;  // assumes that a neutral comes from the vertex
  dir.SetY( tempo.Y() ) ;  
  dir.SetZ( -tempo.Z() ) ; 
  
  dir.SetMag(1.) ;
  return dir ;  
}

//____________________________________________________________________________
Int_t AliPHOSTrackSegment::GetPartType()  
{  
  // Returns 0 - gamma
  //         1 - e+, e-
  //         2 - neutral   
  //         3 - charged hadron

  Int_t PartType = GAMMA ;
                            
  if( fPpsdUp == 0 ) {     // Neutral

    if( fPpsdLow == 0 )    // Neutral  
      PartType = NEUTRAL ;   
    else                   // Gamma
      PartType = GAMMA ;               

  }

  else {             // Charged           

    if( fEmcRecPoint->GetDispersion() > fCutOnDispersion) 
      PartType = CHARGEDHADRON ;
    else  
      PartType = ELECTRON ;  
  
  }
  
  return   PartType ;                     

}

//____________________________________________________________________________
void AliPHOSTrackSegment::GetPosition( TVector3 & pos ) 
{  
  // Returns positions of hits
  TMatrix Dummy ;
  fEmcRecPoint->GetGlobalPosition(pos, Dummy) ;
}

//______________________________________________________________________________
void AliPHOSTrackSegment::Paint(Option_t *)
{
  //Paint this ALiPHOSTrackSegment as a TMarker  with its current attributes

   TVector3 posemc(999., 999., 999.) ;
   TVector3 posppsdl(999., 999., 999.) ;
   TVector3 posppsdu(999., 999., 999.) ;

   fEmcRecPoint->GetLocalPosition(posemc) ;
   if (fPpsdLow !=0 ) 
     fPpsdLow->GetLocalPosition(posppsdl) ;
     if (fPpsdUp !=0 ) 
       fPpsdUp->GetLocalPosition(posppsdu) ;

   Coord_t xemc   = posemc.X() ;
   Coord_t yemc   = posemc.Z() ;

   Coord_t yppsdl = posppsdl.Z() ;
   Coord_t xppsdl = posppsdl.X() ;

   Coord_t yppsdu = posppsdu.Z() ;
   Coord_t xppsdu = posppsdu.X() ;

   Color_t MarkerColor = 1 ;
   Size_t  MarkerSize  = 1.5 ;
   Style_t MarkerStyle = 20 ;

   if (!gPad->IsBatch()) {
     gVirtualX->SetMarkerColor(MarkerColor) ;
     gVirtualX->SetMarkerSize (MarkerSize)  ;
     gVirtualX->SetMarkerStyle(MarkerStyle) ;
   }
   gPad->SetAttMarkerPS(MarkerColor,MarkerStyle,MarkerSize) ;
   gPad->PaintPolyMarker(1, &xemc, &yemc, "") ;
   
   if (xppsdl != 999. && yppsdl != 999. ) {

     MarkerColor = 2 ;
     MarkerSize  = 1.25 ;
     MarkerStyle = 21 ;
     
     if (!gPad->IsBatch()) {
       gVirtualX->SetMarkerColor(MarkerColor) ;
       gVirtualX->SetMarkerSize (MarkerSize)  ;
       gVirtualX->SetMarkerStyle(MarkerStyle) ;
     }
     gPad->SetAttMarkerPS(MarkerColor,MarkerStyle,MarkerSize) ;
     gPad->PaintPolyMarker(1, &xppsdl, &yppsdl, "") ;
   }

    if (xppsdu != 999. && yppsdu != 999. ) {
  
      MarkerColor = 3 ;
      MarkerSize  = 1. ;
      MarkerStyle = 22 ;
      
      if (!gPad->IsBatch()) {
	gVirtualX->SetMarkerColor(MarkerColor) ;
	gVirtualX->SetMarkerSize (MarkerSize)  ;
	gVirtualX->SetMarkerStyle(MarkerStyle) ;
      }
      gPad->SetAttMarkerPS(MarkerColor,MarkerStyle,MarkerSize) ;
      gPad->PaintPolyMarker(1, &xppsdu, &yppsdu, "") ;
    }
}


//____________________________________________________________________________
void AliPHOSTrackSegment::Print()
{
  cout << "--------AliPHOSTrackSegment-------- "<<endl ;
  cout << "EMC Reconstructed Point: " << fEmcRecPoint << endl;
  
  TVector3 pos ;
  TMatrix Dummy ;  

  fEmcRecPoint->GetGlobalPosition( pos, Dummy ) ;
 
  cout << "    position " << pos.X() << "   " << pos.Y() << "  " << pos.Z() << "      Energy " << fEmcRecPoint->GetTotalEnergy() << endl ;
  cout << "PPSD Low Reconstructed Point: " << endl;
  
  if(fPpsdLow){
    fPpsdLow->GetGlobalPosition( pos , Dummy ) ;
    cout << "    position " << pos.X() << "   " << pos.Y() << "  " << pos.Z() << endl ;
  }

  cout << "PPSD Up Reconstructed Point: " << endl;
  
  if(fPpsdUp ){
    fPpsdUp->GetGlobalPosition( pos, Dummy ) ;
    cout << "    position " << pos.X() << "   " << pos.Y() << "  " << pos.Z()  << endl ;
  }

}

