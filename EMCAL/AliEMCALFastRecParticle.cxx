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
//  A  Particle modified by EMCAL response and produced by AliEMCALvFast
//*--
//  To become a general class of AliRoot ?    
//*--               
//*-- Author: Yves Schutz (SUBATECH)
//*--
/////////////////////////////////////////////////////////////////////////////

// --- ROOT system ---

// --- Standard library ---

// --- AliRoot header files ---

#include "AliEMCALFastRecParticle.h"
#include "TPad.h"
#include "TPaveText.h"

ClassImp(AliEMCALFastRecParticle) ; 

//____________________________________________________________________________
AliEMCALFastRecParticle::AliEMCALFastRecParticle() : TParticle()
{
  // ctor
  fType = 0 ; 
}

//____________________________________________________________________________
AliEMCALFastRecParticle::AliEMCALFastRecParticle(const AliEMCALFastRecParticle & rp)
  : TParticle(rp)
{
  // copy ctor

  fType        = rp.fType ;
  fPdgCode     = rp.fPdgCode;
  fStatusCode  = rp.fStatusCode;
  fMother[0]   = rp.fMother[0];
  fMother[1]   = rp.fMother[1];
  fDaughter[0] = rp.fDaughter[0];
  fDaughter[1] = rp.fDaughter[1];
  fWeight      = rp.fWeight;
  fCalcMass    = rp.fCalcMass;
  fPx          = rp.fPx;
  fPy          = rp.fPy;
  fPz          = rp.fPz;
  fE           = rp.fE;
  fVx          = rp.fVx;
  fVy          = rp.fVy;
  fVz          = rp.fVz;
  fVt          = rp.fVt;
  fPolarTheta  = rp.fPolarTheta;
  fPolarPhi    = rp.fPolarPhi;
  fParticlePDG = rp.fParticlePDG; 
}

//____________________________________________________________________________
 AliEMCALFastRecParticle::AliEMCALFastRecParticle(const TParticle & pp)
{
  // ctor from a TParticle (crummy?!)
 
  TParticle & pnoconst = (TParticle &)(pp) ;
  AliEMCALFastRecParticle & p = (AliEMCALFastRecParticle &)(pnoconst) ;
  fType        = 0  ;
  fPdgCode     = p.fPdgCode;
  fStatusCode  = p.fStatusCode;
  fMother[0]   = p.fMother[0];
  fMother[1]   = p.fMother[1];
  fDaughter[0] = p.fDaughter[0];
  fDaughter[1] = p.fDaughter[1];
  fWeight      = p.fWeight;
  fCalcMass    = p.fCalcMass;
  fPx          = p.fPx;
  fPy          = p.fPy;
  fPz          = p.fPz;
  fE           = p.fE;
  fVx          = p.fVx;
  fVy          = p.fVy;
  fVz          = p.fVz;
  fVt          = p.fVt;
  fPolarTheta  = p.fPolarTheta;
  fPolarPhi    = p.fPolarPhi;
  fParticlePDG = p.fParticlePDG; 

}

//____________________________________________________________________________
Int_t AliEMCALFastRecParticle::DistancetoPrimitive(Int_t px, Int_t py)
{
  //  Compute distance from point px,py to a AliEMCALFastRecParticle considered as a Tmarker
  //  Compute the closest distance of approach from point px,py to this marker.
  //  The distance is computed in pixels units.

  Double_t kRADDEG = 180. / TMath::Pi() ; 
  Coord_t x = Phi() * kRADDEG     ;
  Coord_t y = Theta() * kRADDEG     ;
  const Int_t kMaxDiff = 10;
  Int_t pxm  = gPad->XtoAbsPixel(x);
  Int_t pym  = gPad->YtoAbsPixel(y);
  Int_t dist = (px-pxm)*(px-pxm) + (py-pym)*(py-pym);
  
  if (dist > kMaxDiff) return 9999;
  return dist;
}

//___________________________________________________________________________
 void AliEMCALFastRecParticle::Draw(Option_t *option)
 {
   // Draw this AliEMCALFastRecParticle with its current attributes
    
   AppendPad(option);
 }

//______________________________________________________________________________
void AliEMCALFastRecParticle::ExecuteEvent(Int_t event, Int_t /*px*/, Int_t /*py*/)
{
  //  Execute action corresponding to one event
  //  This member function is called when a AliEMCALFastRecParticle is clicked with the locator
  
  if (!gPad->IsEditable()) 
    return ;
  
  static TPaveText * clustertext = 0 ; 
  
  switch (event) {
    
  case kButton1Down: {
    Double_t kRADDEG = 180. / TMath::Pi() ; 
    Coord_t x = Phi() * kRADDEG     ;
    Coord_t y = Theta() * kRADDEG     ;
    clustertext = new TPaveText(x-1, y+1, x+5, y+3, "") ;
    Text_t  line1[40] ;
    Text_t  line2[40] ;
    sprintf( line1, "PID: %s ", (const char*)Name() ) ;
    sprintf( line2, "ENERGY: %f ", Energy() ) ;
    clustertext ->AddText(line1) ;
    clustertext ->AddText(line2) ;
    clustertext ->Draw("");   
    gPad->Update() ; 
    break ;
  }
  
  case kButton1Up: {
    delete clustertext ; 
    clustertext = 0 ; 
    gPad->Update() ; 
    break ;
  }
  }
  
}

//____________________________________________________________________________
TString AliEMCALFastRecParticle::Name()const
{
  // Returns the name of the particle type (only valid if PIDv1 is employed)

  TString  name ; 
  
  if(fType == 127)
    name = "PHOTON_LOPU_HIEF" ;    //PCA = 001 TOF = 111 CPV = 111
  
  if(fType == 511)
    name = "PHOTON_HIPU_LOEF" ;    //PCA = 011 TOF = 111 CPV = 111
  
  if(fType == 255)
	name = "PHOTON_MED_PU_EF" ;    //PCA = 111 TOF = 111 CPV = 111
  
  if((fType == 383)||(fType == 447)) 
    name = "PHOTON_STRANGE" ;      //PCA = 101 or 110 TOF = 111 CPV = 111
  
  if(fType == 63)
    name = "NEUTRAL_FAST_HADRON" ; //PCA = 000 TOF = 111 CPV = 111
  
  if((fType == 504) || (fType == 505) ||(fType == 248)||(fType == 249)||(fType == 120)||(fType == 121))
    name = "CHARGED_FAST_EM" ;     //PCA = 111, 011 or 001 TOF =111 CPV = 000 or 001  
  
  if((fType == 56)||(fType == 57))
    name = "CHARGED_FAST_HADRON" ; //PCA = 000 TOF = 111 CPV = 000 or 001 
  
  if((fType < 8)&&(fType > 0))
    name = "NEUTRAL_SLOW_HADRON" ; //PCA = 000 TOF = 000 CPV = 001 or 011 or 111
  
  if((fType == 0))
    name = "CHARGED_SLOW_HADRON" ; //PCA = 000 TOF = 000 CPV = 000
  
  if((fType == 448) || (fType == 449) ||(fType == 192)||(fType == 193)||(fType == 64)||(fType == 64))
    name = "CHARGED_SLOW_EM" ;    //PCA = 111, 011 or 001 TOF =000 CPV = 000 or 001  
  
  return name ; 
}


//______________________________________________________________________________
void AliEMCALFastRecParticle::SetType(Int_t type) { 
  // sets the particle type 
  fType = type ; 
  
  if((type == 127) || (fType == 511) || (fType == 255) ||(fType == 383)||(fType == 447)){
    fPdgCode = 22 ; 
    return ;
  }
  
  if ((fType == 63)|| ((fType < 8)&&(fType > 0)) ){
    fPdgCode = 2112 ; 
    return ;
  }
  if ( ((fType == 504) || (fType == 505) ||(fType == 248)||(fType == 249)||(fType == 120)||(fType == 121)) ){
    fPdgCode = 11 ; 
    return ;
  }
  if ((fType == 448) || (fType == 449) ||(fType == 192)||(fType == 193)||(fType == 64)||(fType == 64)){
    fPdgCode = 13 ; 
    return ;
  }
  if((fType == 56)||(fType == 57)){
    fPdgCode = 211 ; 
    return ;
  }
  if (fType == 0){
    fPdgCode = 2212 ; 
    return ;
  }

}	    

//______________________________________________________________________________
void AliEMCALFastRecParticle::Paint(Option_t *)
{
  // Paint this ALiRecParticle in theta,phi coordinate as a TMarker  with its current attributes

  Double_t kRADDEG = 180. / TMath::Pi() ; 
  Coord_t x = Phi() * kRADDEG     ;
  Coord_t y = Theta() * kRADDEG     ;
  Color_t markercolor = 1 ;
  Size_t  markersize  = 1. ;
  Style_t markerstyle = 5 ;
  
  if (!gPad->IsBatch()) {
    gVirtualX->SetMarkerColor(markercolor) ;
    gVirtualX->SetMarkerSize (markersize)  ;
    gVirtualX->SetMarkerStyle(markerstyle) ;
  }
  gPad->SetAttMarkerPS(markercolor,markerstyle,markersize) ;
  gPad->PaintPolyMarker(1,&x,&y,"") ;
}

//____________________________________________________________________________
void AliEMCALFastRecParticle::Print(Option_t * /*opt*/)const
{
  // Print the type, energy and momentum of the reconstructed particle
  
  printf("Print Summary:") ; 
  printf("AliEMCALFastRecParticle > type is  %s\n", Name().Data()) ; 
  printf("                      Energy = %f\n", fE) ; 
  printf("                         Px     = %f\n", fPx) ; 
  printf("                         Py     = %f\n", fPy) ;
  printf("                         Pz     = %f\n", fPz) ; 
}
