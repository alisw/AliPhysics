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
//  Track segment in PHOS
//  Can be : 1 EmcRecPoint
//           1 EmcRecPoint + 1 PPSD
//           1 EmcRecPoint + 1 PPSD + 1 PPSD     
//                  
//*-- Author:  Dmitri Peressounko (RRC KI & SUBATECH)

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
  // ctor

  if( emc )   
    fEmcRecPoint =  emc ;

  if( ppsdRP1 )  
    fPpsdUp = ppsdRP1 ;

  if( ppsdRP2  ) 
    fPpsdLow = ppsdRP2 ;

}

//____________________________________________________________________________
AliPHOSTrackSegment::AliPHOSTrackSegment( const AliPHOSTrackSegment & ts) 
{
  // Copy ctor

  ( (AliPHOSTrackSegment &)ts ).Copy(*this) ; 
}


//____________________________________________________________________________
void AliPHOSTrackSegment::Copy(TObject & obj) 
{
  // Copy of a track segment into another track segment

   TObject::Copy(obj) ;
   ( (AliPHOSTrackSegment &)obj ).fEmcRecPointId = fEmcRecPointId ; 
   ( (AliPHOSTrackSegment &)obj ).fPpsdLowId     = fPpsdLowId ; 
   ( (AliPHOSTrackSegment &)obj ).fPpsdUpId      = fPpsdUpId ; 
}
//____________________________________________________________________________
Int_t AliPHOSTrackSegment::DistancetoPrimitive(Int_t px, Int_t py)
{
  // Compute distance from point px,py to  a AliPHOSTrackSegment considered as a Tmarker
  // Compute the closest distance of approach from point px,py to this marker.
  // The distance is computed in pixels units.
  
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
 
  static TPaveText* textTS = 0 ;
  
  if (!gPad->IsEditable()) 
    return;
  
  switch (event) {
    
  case kButton1Down:{
    
    if (!textTS) {
      
      TVector3 pos(0.,0.,0.) ;
      fEmcRecPoint->GetLocalPosition(pos) ;
      textTS = new TPaveText(pos.X()-10,pos.Z()+10,pos.X()+5,pos.Z()+15,"") ;
      Text_t  line1[40] ;
      sprintf(line1,"See RecParticle for ID") ;
      textTS ->AddText(line1) ;
      textTS ->Draw("");
      gPad->Update() ; 
    }
  }
  
  break;
  
  case kButton1Up:
    if (textTS) {
      delete textTS ;
      textTS = 0 ;
    }
    break;  
  }
}

//____________________________________________________________________________
Float_t AliPHOSTrackSegment::GetDistanceInPHOSPlane()
{
  // Calculates the distance between the EMC RecPoint and PPSD RecPoint
  
  TVector3 vecEmc ;
  fEmcRecPoint->GetLocalPosition(vecEmc) ;
  
  TVector3 vecPpsd ;
  if( fPpsdLow->GetMultiplicity() )  
    fPpsdLow->GetLocalPosition(vecPpsd)  ; 
  else { 
    vecPpsd.SetX(10000.) ;
  } 
  vecEmc -= vecPpsd ;

  Float_t r = vecEmc.Mag();;

  return r ;
}

//____________________________________________________________________________
TVector3 AliPHOSTrackSegment::GetMomentumDirection() 
{ 
  // Calculates the momentum direction:
  //   1. if only a EMC RecPoint, direction is given by IP and this RecPoint
  //   2. if a EMC RecPoint and one PPSD RecPoint, direction is given by the line through the 2 recpoints 
  //   3. if a EMC RecPoint and two PPSD RecPoints, dirrection is given by the average line through 
  //      the 2 pairs of recpoints  
  // However because of the poor position resolution of PPSD the direction is always taken as if we were 
  //  in case 1.

  TVector3 dir(0,0,0) ; 
  TMatrix mdummy ;

  TVector3 posEmc ;
  fEmcRecPoint->GetGlobalPosition(posEmc, mdummy) ;

  // Correction for the depth of the shower starting point (TDR p 127) 

  Float_t energy = fEmcRecPoint->GetEnergy() ; 
  Float_t para = 0.925 ; 
  Float_t parb = 6.52 ; 

  TVector3 localpos ; 
  fEmcRecPoint->GetLocalPosition(localpos) ; 

  AliPHOSGeometry * geom = AliPHOSGeometry::GetInstance() ; 
  Float_t radius = geom->GetIPtoOuterCoverDistance() + geom->GetOuterBoxSize(1) ; 
  Float_t incidencephi = TMath::ATan(localpos.X() / radius) ; 
  Float_t incidencetheta = TMath::ATan(localpos.Z() / radius) ;
 
  Float_t depthx =  ( para * TMath::Log(energy) + parb ) * TMath::Sin(incidencephi) ; 
  Float_t depthz =  ( para * TMath::Log(energy) + parb ) * TMath::Sin(incidencetheta) ; 
  
  localpos.SetX(localpos.X() - depthx ) ;
  localpos.SetZ(localpos.Z() - depthz ) ;
 
  TVector3 emcglobalpos ;
  TMatrix  dummy ;

  fEmcRecPoint->GetGlobalPosition(emcglobalpos, dummy) ;

  emcglobalpos.SetX( emcglobalpos.X() - depthx ) ;  
  emcglobalpos.SetZ( emcglobalpos.Z() - depthz ) ;   
  
  // The following commeneted code becomes valid once the PPSD provides 
  // a reasonable position resolution, at least as good as EMC !
 
//   TVector3 ppsdlglobalpos ;
//   TVector3 ppsduglobalpos ;

//   if( fPpsdLow ){ // certainly a photon that has concerted
        
//     fPpsdLow->GetGlobalPosition(ppsdlglobalpos, mdummy) ; 
//     dir = emcglobalpos -  ppsdlglobalpos ; 
     
//     if( fPpsdUp ){ // nop looks like a charged
       
//        fPpsdUp->GetGlobalPosition(ppsduglobalpos, mdummy) ; 
//        dir = ( dir +  emcglobalpos -  ppsduglobalpos ) * 0.5 ; 
//      }
//   }
 
//   else { // looks like a neutral

    dir = emcglobalpos ;  
//  }

  dir.SetZ( -dir.Z() ) ;   // why ?  
  dir.SetMag(1.) ;
    
  return dir ;  
}

//____________________________________________________________________________
Int_t *  AliPHOSTrackSegment::GetPrimariesEmc(Int_t & number) 
{ 
  // Retrieves the primary particle(s) at the origin of the EMC RecPoint
  
  Int_t * rv = 0 ; 
  number = 0 ;
  if ( fEmcRecPoint )
    rv =  fEmcRecPoint->GetPrimaries(number) ; 

  return rv ; 
}

//____________________________________________________________________________
Int_t *  AliPHOSTrackSegment::GetPrimariesPpsdLow(Int_t & number) 
{ 
  // Retrieves the primary particle(s) at the origin of the lower PPSD RecPoint
  
  Int_t * rv = 0 ; 
  number = 0 ;
  if ( fPpsdLow )
    rv =  fPpsdLow->GetPrimaries(number) ; 

  return rv ; 
}

//____________________________________________________________________________
Int_t *  AliPHOSTrackSegment::GetPrimariesPpsdUp(Int_t & number) 
{ 
  // Retrieves the primary particle(s) at the origin of the upper PPSD  RecPoint
  
  Int_t * rv = 0 ; 
  number = 0 ;
  if ( fPpsdUp )
    rv =  fPpsdUp->GetPrimaries(number) ; 

  return rv ; 
}

//____________________________________________________________________________
void AliPHOSTrackSegment::GetPosition( TVector3 & pos ) 
{  
  // Returns position of the EMC RecPoint

  TMatrix dummy ;
  fEmcRecPoint->GetGlobalPosition(pos, dummy) ;
}


//______________________________________________________________________________
void AliPHOSTrackSegment::Paint(Option_t *)
{
  // Paint this ALiPHOSTrackSegment as a TMarker  with its current attributes

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
  
  Color_t markercolor = 1 ;
  Size_t  markersize  = 1.5 ;
  Style_t markerstyle = 20 ;
  
  if (!gPad->IsBatch()) {
    gVirtualX->SetMarkerColor(markercolor) ;
    gVirtualX->SetMarkerSize (markersize)  ;
    gVirtualX->SetMarkerStyle(markerstyle) ;
  }
  gPad->SetAttMarkerPS(markercolor,markerstyle,markersize) ;
  gPad->PaintPolyMarker(1, &xemc, &yemc, "") ;
  
  if (xppsdl != 999. && yppsdl != 999. ) {
    
    markercolor = 2 ;
    markersize  = 1.25 ;
    markerstyle = 21 ;
    
    if (!gPad->IsBatch()) {
      gVirtualX->SetMarkerColor(markercolor) ;
      gVirtualX->SetMarkerSize (markersize)  ;
      gVirtualX->SetMarkerStyle(markerstyle) ;
    }
    gPad->SetAttMarkerPS(markercolor,markerstyle,markersize) ;
    gPad->PaintPolyMarker(1, &xppsdl, &yppsdl, "") ;
  }
  
  if (xppsdu != 999. && yppsdu != 999. ) {
    
    markercolor = 3 ;
    markersize  = 1. ;
    markerstyle = 22 ;
    
    if (!gPad->IsBatch()) {
      gVirtualX->SetMarkerColor(markercolor) ;
      gVirtualX->SetMarkerSize (markersize)  ;
      gVirtualX->SetMarkerStyle(markerstyle) ;
    }
    gPad->SetAttMarkerPS(markercolor,markerstyle,markersize) ;
    gPad->PaintPolyMarker(1, &xppsdu, &yppsdu, "") ;
  }
}

//____________________________________________________________________________
void AliPHOSTrackSegment::Print()
{
  // Print all information on this track Segment
  
  cout << "--------AliPHOSTrackSegment-------- "<<endl ;
  cout << "EMC Reconstructed Point: " << fEmcRecPoint << endl;
  
  TVector3 pos ;
  TMatrix dummy ;  

  fEmcRecPoint->GetGlobalPosition( pos, dummy ) ;
 
  cout << "    position " << pos.X() << "   " << pos.Y() << "  " << pos.Z() << "      Energy " << fEmcRecPoint->GetTotalEnergy() << endl ;
  cout << "PPSD Low Reconstructed Point: " << endl;
  
  if(fPpsdLow){
    fPpsdLow->GetGlobalPosition( pos , dummy ) ;
    cout << "    position " << pos.X() << "   " << pos.Y() << "  " << pos.Z() << endl ;
  }

  cout << "PPSD Up Reconstructed Point: " << endl;
  
  if(fPpsdUp ){
    fPpsdUp->GetGlobalPosition( pos, dummy ) ;
    cout << "    position " << pos.X() << "   " << pos.Y() << "  " << pos.Z()  << endl ;
  }

}

