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
AliPHOSTrackSegment::~AliPHOSTrackSegment() // dtor
{
//    fEmcRecPoint.Delete() ;   Not Owners !!!
//    fPpsdUp.Delete() ;
//    fPpsdLow.Delete() ;
}

//____________________________________________________________________________
Int_t AliPHOSTrackSegment::DistancetoPrimitive(Int_t px, Int_t py)
{
//*-*-*-*-*-*-*-*-*-*-*Compute distance from point px,py to  a AliPHOSTrackSegment considered as a Tmarker*-*-*-*-*-*
//*-*                  ===========================================
//  Compute the closest distance of approach from point px,py to this marker.
//  The distance is computed in pixels units.
//

   TVector3 pos(0.,0.,0.) ;
   fEmcRecPoint->GetLocalPosition( pos) ;
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
 void AliPHOSTrackSegment::Draw(Option_t *option)
 {
// //*-*-*-*-*-*-*-*-*-*-*Draw this AliPHOSTrackSegment with its current attributes*-*-*-*-*-*-*
// //*-*
   // assert(0==1);
  AppendPad(option);
 }

//______________________________________________________________________________
void AliPHOSTrackSegment::ExecuteEvent(Int_t event, Int_t px, Int_t py)
{
//*-*-*-*-*-*-*-*-*-*-*Execute action corresponding to one event*-*-*-*
//*-*                  =========================================
//  This member function is called when a AliPHOSTrackSegment is clicked with the locator
//
//  If Left button is clicked on AliPHOSRecPoint, the digits are switched on    
//  and switched off when the mouse button is released.
//
   static TPaveText* TrackSegmentText = 0 ;

   if (!gPad->IsEditable()) return;

   switch (event) {

   case kButton1Down:{
    
     if (!TrackSegmentText) {
  
       TVector3 pos(0.,0.,0.) ;
       fEmcRecPoint->GetLocalPosition(pos) ;
       TrackSegmentText = new TPaveText(pos.X()-10,pos.Z()+10,pos.X()+50,pos.Z()+35,"") ;
       Text_t  line1[40] ;
       if (GetPartType() == 0 ) sprintf(line1,"PHOTON") ;
       if (GetPartType() == 1 ) sprintf(line1,"NEUTRAL HADRON") ;
       if (GetPartType() == 2 ) sprintf(line1,"CHARGED HADRON") ;
       if (GetPartType() == 3 ) sprintf(line1,"ELECTRON") ;
       TrackSegmentText ->AddText(line1) ;
       TrackSegmentText ->Paint("");
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
Bool_t AliPHOSTrackSegment::GetMomentumDirection( TVector3 & dir ) 
{   
  // True if determined
  Bool_t ifdeterm = kTRUE ;

  if( fPpsdLow ){
    TMatrix mdummy ;
    if( fPpsdUp->GetMultiplicity() ) { // draw line trough 2 points
      TVector3 posEmc ;
      fEmcRecPoint->GetGlobalPosition(posEmc,mdummy) ;
      TVector3 posPpsd ;
      fPpsdLow->GetGlobalPosition(posPpsd,mdummy) ; 
      dir = posEmc - posPpsd ;
      dir.SetMag(1.) ;
    }

    else { // draw line through 3 pionts
      TVector3 posEmc ;
      fEmcRecPoint->GetGlobalPosition(posEmc,mdummy) ;
      TVector3 posPpsdl ;
      fPpsdLow->GetGlobalPosition(posPpsdl,mdummy) ; 
      TVector3 posPpsdup ;
      fPpsdUp->GetGlobalPosition(posPpsdup,mdummy) ;
      posPpsdl = 0.5*( posPpsdup+posPpsdl ) ; 
      dir = posEmc - posPpsdl ;
      dir.SetMag(1.) ;
    }
  
  }
  else 
    ifdeterm = kFALSE ;
 
  return ifdeterm ;
}

//____________________________________________________________________________
Int_t AliPHOSTrackSegment::GetPartType()  
{  
  // Returns 0 - gamma
  //         1 - e+, e-
  //         2 - neutral hadron  
  //         3 - charged hadron

  Int_t PartType =0;                            
  if( fPpsdUp ){     // Neutral

    if( fPpsdLow ) // Neutral hadron  
      PartType = 2 ;   
    else                                // Gamma
      PartType = 0 ;               

  }

  else {             // Charged           

    if( fEmcRecPoint->GetDispersion() > fCutOnDispersion) 
      PartType = 3 ;
    else  
      PartType = 1 ;  
  
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
//*-*-*-*-*-*-*-*-*-*-*Paint this ALiPHOSTrackSegment as a TMarker  with its current attributes*-*-*-*-*-*-*
//*-*                  =============================================
   TVector3 pos(0.,0.,0.)  ;
   fEmcRecPoint->GetLocalPosition(pos)   ;
   Coord_t x = pos.X()     ;
   Coord_t y = pos.Z()     ;
   Color_t MarkerColor = 1 ;
   Size_t  MarkerSize = 1. ;
   Style_t MarkerStyle = 29 ;

   if (GetPartType() == 0 ) MarkerStyle = 20 ;
   if (GetPartType() == 1 ) MarkerStyle = 21 ;
   if (GetPartType() == 2 ) MarkerStyle = 22 ;
   if (GetPartType() == 3 ) MarkerStyle = 23 ;

   if (!gPad->IsBatch()) {
     gVirtualX->SetMarkerColor(MarkerColor) ;
     gVirtualX->SetMarkerSize (MarkerSize)  ;
     gVirtualX->SetMarkerStyle(MarkerStyle) ;
   }
   gPad->SetAttMarkerPS(MarkerColor,MarkerStyle,MarkerSize) ;
   gPad->PaintPolyMarker(1,&x,&y,"") ;
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

