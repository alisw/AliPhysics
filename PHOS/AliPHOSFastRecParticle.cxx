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
//  A  Particle modified by PHOS response and produced by AliPHOSvFast
//  To become a general class of AliRoot ?    
//               
//*-- Author: Yves Schutz (SUBATECH)

// --- ROOT system ---

// --- Standard library ---

#include <iostream.h>

// --- AliRoot header files ---

#include "AliPHOSFastRecParticle.h"
#include "TPad.h"
#include "TPaveText.h"

ClassImp(AliPHOSFastRecParticle)

//____________________________________________________________________________
 AliPHOSFastRecParticle::AliPHOSFastRecParticle(const AliPHOSFastRecParticle & rp)
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
 AliPHOSFastRecParticle::AliPHOSFastRecParticle(const TParticle & pp)
{
  // ctor from a TParticle (crummy?!)
 
  TParticle & pnoconst = (TParticle &)(pp) ;
  AliPHOSFastRecParticle & p = (AliPHOSFastRecParticle &)(pnoconst) ;

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
Int_t AliPHOSFastRecParticle::DistancetoPrimitive(Int_t px, Int_t py)
{
  //  Compute distance from point px,py to a AliPHOSFastRecParticle considered as a Tmarker
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
 void AliPHOSFastRecParticle::Draw(Option_t *option)
 {
   // Draw this AliPHOSFastRecParticle with its current attributes
    
   AppendPad(option);
 }

//______________________________________________________________________________
void AliPHOSFastRecParticle::ExecuteEvent(Int_t event, Int_t px, Int_t py)
{
  //  Execute action corresponding to one event
  //  This member function is called when a AliPHOSFastRecParticle is clicked with the locator
     
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
Int_t * AliPHOSFastRecParticle::GetPrimaries(Int_t & number) 
{
  // Retrieves the unique primary particle at the origine of the present reconstruced particle

  number = 1 ; 
  Int_t * list = new Int_t[1] ;
  list[0] = fPrimary ; 
  
  return list ;
}

//____________________________________________________________________________
TString AliPHOSFastRecParticle::Name()
{
  // Returns the name of the particle type
  
  TString  name ; 
  switch (fType) {
  case kGAMMA:
    name = "PHOTON" ;
    break ; 
   case kELECTRON:
     name = "ELECTRON" ;
    break ; 
   case kCHARGEDHA:
    name = "CHARGED_HA" ;
    break ; 
  case kNEUTRALHA:
    name = "NEUTRAL_HA" ; 
    break ; 
  case kNEUTRALEM:
    name = "NEUTRAL_EM" ; 
    break ; 
  case kGAMMAHA:
    name = "PHOTON_HA" ; 
    break ; 

  }
  return name ; 
}

//______________________________________________________________________________
void AliPHOSFastRecParticle::Paint(Option_t *)
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
void AliPHOSFastRecParticle::Print()
{
  // Print the typr, energy and momentum
  
  cout << "AliPHOSFastRecParticle > " << "type is  " << Name() << endl 
       << "                     " << "Energy = " << fE << endl 
       << "                     " << "Px     = " << fPx << endl 
       << "                     " << "Py     = " << fPy << endl 
       << "                     " << "Pz     = " << fPz << endl ; 
}
