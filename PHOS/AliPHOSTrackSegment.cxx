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

#include <iostream.h>

// --- AliRoot header files ---

#include "AliPHOSTrackSegment.h" 
#include "AliPHOSv0.h"
#include "AliPHOSIndexToObject.h"

ClassImp(AliPHOSTrackSegment)

//____________________________________________________________________________
AliPHOSTrackSegment::AliPHOSTrackSegment( AliPHOSEmcRecPoint * emc , 
					  AliPHOSPpsdRecPoint * ppsdrp1,
					  AliPHOSPpsdRecPoint * ppsdrp2  ) 
{
  // ctor

  if( emc )   
    fEmcRecPoint =  emc->GetIndexInList() ;
  else 
    fEmcRecPoint = -1 ;

  if( ppsdrp1 )  
    fPpsdUpRecPoint = ppsdrp1->GetIndexInList() ;
 else 
    fPpsdUpRecPoint = -1 ;

  if( ppsdrp2  ) 
    fPpsdLowRecPoint = ppsdrp2->GetIndexInList() ;
  else 
    fPpsdLowRecPoint = -1 ;

  fIndexInList = -1 ;
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
   ( (AliPHOSTrackSegment &)obj ).fEmcRecPoint     = fEmcRecPoint ; 
   ( (AliPHOSTrackSegment &)obj ).fPpsdLowRecPoint = fPpsdLowRecPoint ; 
   ( (AliPHOSTrackSegment &)obj ).fPpsdUpRecPoint  = fPpsdUpRecPoint ; 
   ( (AliPHOSTrackSegment &)obj ).fIndexInList     = fIndexInList ; 
}
//____________________________________________________________________________
Int_t AliPHOSTrackSegment::DistancetoPrimitive(Int_t px, Int_t py)
{
  // Compute distance from point px,py to  a AliPHOSTrackSegment considered as a Tmarker
  // Compute the closest distance of approach from point px,py to this marker.
  // The distance is computed in pixels units.
  
  Int_t div = 1 ;  
  Int_t dist = 9999 ; 
  
  TVector3 pos(0.,0.,0.) ;
  
  AliPHOSEmcRecPoint  * emcrp   = GetEmcRecPoint() ; 
  AliPHOSPpsdRecPoint * ppsdlrp = GetPpsdLowRecPoint() ; 
  AliPHOSPpsdRecPoint * ppsdurp = GetPpsdUpRecPoint() ; 
  
  if ( emcrp != 0 ) {
    emcrp->GetLocalPosition( pos) ;
    Float_t x =  pos.X() ;
    Float_t y =  pos.Z() ;
    if ( ppsdlrp != 0 ) {
      ppsdlrp->GetLocalPosition( pos ) ;
      x +=  pos.X() ;
      y +=  pos.Z() ;
      div++ ; 
    }
    if ( ppsdurp != 0 ) {
      ppsdurp->GetLocalPosition( pos ) ;
      x +=  pos.X() ;
      y +=  pos.Z() ;
      div++ ; 
    }
    x /= div ; 
    y /= div ; 

    const Int_t kMaxDiff = 10;
    Int_t pxm  = gPad->XtoAbsPixel(x);
    Int_t pym  = gPad->YtoAbsPixel(y);
    dist = (px-pxm)*(px-pxm) + (py-pym)*(py-pym);
    
    if (dist > kMaxDiff) return 9999;
  }
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
  
  AliPHOSEmcRecPoint  * emcrp   = GetEmcRecPoint() ; 

  if (!gPad->IsEditable()) 
    return;
  
  switch (event) {
    
  case kButton1Down:{
    
    if (!textTS) {
      
      TVector3 pos(0.,0.,0.) ;
      emcrp->GetLocalPosition(pos) ;
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
  
  AliPHOSEmcRecPoint  * emcrp   = GetEmcRecPoint() ; 
  AliPHOSPpsdRecPoint * ppsdlrp = GetPpsdLowRecPoint() ; 

  TVector3 vecEmc ;
  emcrp->GetLocalPosition(vecEmc) ;
  
  TVector3 vecPpsd ;
  if ( ppsdlrp !=0 ) {
    if( ppsdlrp->GetMultiplicity() )  
      ppsdlrp->GetLocalPosition(vecPpsd)  ; 
    else { 
      vecPpsd.SetX(10000.) ;
    } 
    vecEmc -= vecPpsd ;
  }
  Float_t r = vecEmc.Mag();;

  return r ;
}

//____________________________________________________________________________
AliPHOSEmcRecPoint * AliPHOSTrackSegment::GetEmcRecPoint() const 
{
  // get the EMC recpoint at the origin of this track
 
  AliPHOSIndexToObject * please =  AliPHOSIndexToObject::GetInstance() ;
  AliPHOSEmcRecPoint * rv = 0 ;
  if (  fEmcRecPoint > -1 )
    rv = (AliPHOSEmcRecPoint *)please->GimeRecPoint( fEmcRecPoint, TString("emc") );
  
  return rv ;

}
  
//____________________________________________________________________________
 Float_t AliPHOSTrackSegment::GetEnergy()
{ 
  // Returns energy in EMC
  
  AliPHOSEmcRecPoint  * emcrp   = GetEmcRecPoint() ; 
  
  return emcrp->GetEnergy() ;
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


  AliPHOSEmcRecPoint  * emcrp   = GetEmcRecPoint() ; 
  // AliPHOSPpsdRecPoint * ppsdlrp = GetPpsdLowRecPoint() ; 
  // AliPHOSPpsdRecPoint * ppsdurp = GetPpsdUpRecPoint() ; 

  TVector3 dir(0,0,0) ; 
  TMatrix mdummy ;

  TVector3 posEmc ;
  emcrp->GetGlobalPosition(posEmc, mdummy) ;
 
  TVector3 emcglobalpos ;
  TMatrix  dummy ;

  emcrp->GetGlobalPosition(emcglobalpos, dummy) ;

  
// The following commeneted code becomes valid once the PPSD provides 
// a reasonable position resolution, at least as good as EMC ! 
//   TVector3 ppsdlglobalpos ;
//   TVector3 ppsduglobalpos ;
//   if( fPpsdLowRecPoint ){ // certainly a photon that has concerted
//     fPpsdLowRecPoint->GetGlobalPosition(ppsdlglobalpos, mdummy) ; 
//     dir = emcglobalpos -  ppsdlglobalpos ; 
//     if( fPpsdUpRecPoint ){ // not looks like a charged       
//        fPpsdUpRecPoint->GetGlobalPosition(ppsduglobalpos, mdummy) ; 
//        dir = ( dir +  emcglobalpos -  ppsduglobalpos ) * 0.5 ; 
//      }
//   }
//   else { // looks like a neutral
//    dir = emcglobalpos ;  
//  }

  dir = emcglobalpos ;  
  dir.SetZ( -dir.Z() ) ;   // why ?  
  dir.SetMag(1.) ;
    
  return dir ;  
}

//____________________________________________________________________________
Int_t AliPHOSTrackSegment:: GetPHOSMod(void) 
{
  // Returns the phos module which contains this track
 
  AliPHOSEmcRecPoint  * emcrp   = GetEmcRecPoint() ; 
  
  return emcrp->GetPHOSMod();  
}

//____________________________________________________________________________
AliPHOSPpsdRecPoint * AliPHOSTrackSegment::GetPpsdLowRecPoint() const 
{
  // Returns the lower PPSD rec point at the origin of this track
  
  AliPHOSIndexToObject * please =  AliPHOSIndexToObject::GetInstance() ;
  AliPHOSPpsdRecPoint * rv = 0 ;
  
  if ( fPpsdLowRecPoint > -1 )
    rv = (AliPHOSPpsdRecPoint *)please->GimeRecPoint( fPpsdLowRecPoint, TString("ppsd") ) ;
  
  return rv ; 
}

//____________________________________________________________________________
AliPHOSPpsdRecPoint * AliPHOSTrackSegment::GetPpsdUpRecPoint() const 
{
  // Returns the lower PPSD rec point at the origin of this track

  AliPHOSIndexToObject * please =  AliPHOSIndexToObject::GetInstance() ;
  AliPHOSPpsdRecPoint * rv = 0 ;
 
  if ( fPpsdUpRecPoint > -1 )
    rv =  (AliPHOSPpsdRecPoint *)please->GimeRecPoint( fPpsdUpRecPoint, TString("ppsd") ) ;

  return rv ;
}

//____________________________________________________________________________
Int_t *  AliPHOSTrackSegment::GetPrimariesEmc(Int_t & number) 
{ 
  // Retrieves the primary particle(s) at the origin of the EMC RecPoint
    
  AliPHOSEmcRecPoint  * emcrp   = GetEmcRecPoint() ; 

  Int_t * rv = 0 ; 
  number = 0 ;
  if ( emcrp )
    rv =  emcrp->GetPrimaries(number) ; 

  return rv ; 
}

//____________________________________________________________________________
Int_t *  AliPHOSTrackSegment::GetPrimariesPpsdLow(Int_t & number) 
{ 
  // Retrieves the primary particle(s) at the origin of the lower PPSD RecPoint
  
  AliPHOSPpsdRecPoint * ppsdlrp = GetPpsdLowRecPoint() ; 

  Int_t * rv = 0 ; 
  number = 0 ;
  if ( ppsdlrp )
    rv =  ppsdlrp->GetPrimaries(number) ; 

  return rv ; 
}

//____________________________________________________________________________
Int_t *  AliPHOSTrackSegment::GetPrimariesPpsdUp(Int_t & number) 
{ 
  // Retrieves the primary particle(s) at the origin of the upper PPSD  RecPoint
  
  AliPHOSPpsdRecPoint * ppsdurp = GetPpsdUpRecPoint() ; 

  Int_t * rv = 0 ; 
  number = 0 ;
  if ( ppsdurp )
    rv =  ppsdurp->GetPrimaries(number) ; 

  return rv ; 
}

//____________________________________________________________________________
void AliPHOSTrackSegment::GetPosition( TVector3 & pos ) 
{  
  // Returns position of the EMC RecPoint
  
  AliPHOSEmcRecPoint  * emcrp   = GetEmcRecPoint() ; 
 
  TMatrix dummy ;
  emcrp->GetGlobalPosition(pos, dummy) ;
}


//______________________________________________________________________________
void AliPHOSTrackSegment::Paint(Option_t *)
{
  // Paint this AliPHOSTrackSegment as a TMarker  with its current attributes

  AliPHOSEmcRecPoint  * emcrp   = GetEmcRecPoint() ; 
  AliPHOSPpsdRecPoint * ppsdlrp = GetPpsdLowRecPoint() ; 
  AliPHOSPpsdRecPoint * ppsdurp = GetPpsdUpRecPoint() ; 

  TVector3 posemc(999., 999., 999.) ;
  TVector3 posppsdl(999., 999., 999.) ;
  TVector3 posppsdu(999., 999., 999.) ;
  
  emcrp->GetLocalPosition(posemc) ;
  if (ppsdlrp !=0 ) 
    ppsdlrp->GetLocalPosition(posppsdl) ;
  if (ppsdurp !=0 ) 
    ppsdurp->GetLocalPosition(posppsdu) ;
  
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
void AliPHOSTrackSegment::Print(const char * opt)
{
  // Print all information on this track Segment
  
  AliPHOSEmcRecPoint  * emcrp   = GetEmcRecPoint() ; 
  AliPHOSPpsdRecPoint * ppsdlrp = GetPpsdLowRecPoint() ; 
  AliPHOSPpsdRecPoint * ppsdurp = GetPpsdUpRecPoint() ; 
  
  TVector3 pos ;
  TMatrix dummy ;  

  cout << "--------AliPHOSTrackSegment-------- "<<endl ;

  if ( emcrp != 0 ) {
    cout << "******** EMC Reconstructed Point: " << endl;
    emcrp->Print() ; 
    
    emcrp->GetGlobalPosition( pos, dummy ) ;
    
    cout << " Global position " << pos.X() << "   " << pos.Y() << "  " << pos.Z() << "      Energy " << emcrp->GetEnergy() << endl ;
  }
  
  if ( ppsdlrp != 0 ) {
    cout << "******** PPSD Low Reconstructed Point: " << endl;
    
    ppsdlrp->Print() ; 
    ppsdlrp->GetGlobalPosition( pos , dummy ) ;
    cout << "    position " << pos.X() << "   " << pos.Y() << "  " << pos.Z() << endl ;
  }

   if( ppsdurp != 0 ) {
     cout << "******** PPSD Up Reconstructed Point: " << endl;
     
     ppsdurp->Print() ; 
     ppsdurp->GetGlobalPosition( pos, dummy ) ;
     cout << "    position " << pos.X() << "   " << pos.Y() << "  " << pos.Z()  << endl ;
   }
   
}
//____________________________________________________________________________
void AliPHOSTrackSegment::SetPpsdUpRecPoint(AliPHOSPpsdRecPoint * PpsdUpRecPoint) 
{
  // gives an id from its position in the list
  if( PpsdUpRecPoint )  
    fPpsdUpRecPoint = PpsdUpRecPoint->GetIndexInList() ;
 else 
    fPpsdUpRecPoint = -1 ;
}

