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
// Reconstructed Particle
//*-- Y. Schutz:   SUBATECH 
//////////////////////////////////////////////////////////////////////////////

// --- ROOT system ---

// --- Standard library ---

// --- AliRoot header files ---

#include "AliPHOSRecParticle.h"
#include "TPad.h"

ClassImp(AliPHOSRecParticle)


//____________________________________________________________________________
 AliPHOSRecParticle::AliPHOSRecParticle(AliPHOSTrackSegment * ts)
{
  // ctor
 
  fPHOSTrackSegment = new AliPHOSTrackSegment(*ts) ; 
  fE                = ts->GetEnergy() ; 
  TVector3 momdir   = ts->GetMomentumDirection() ;
  fPx               = fE * momdir.X() ; 
  fPy               = fE * momdir.Y() ; 
  fPz               = fE * momdir.Z() ; 
  fType             = kUNDEFINED ;  
                           
}

//____________________________________________________________________________
 AliPHOSRecParticle::AliPHOSRecParticle(const AliPHOSRecParticle & rp)
{
  fPHOSTrackSegment = new AliPHOSTrackSegment( *( rp.GetPHOSTrackSegment()) ) ; 
  fType             = rp.fType ; 
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
 AliPHOSRecParticle::~AliPHOSRecParticle()
{
  if(!fPHOSTrackSegment) {
    delete fPHOSTrackSegment ;
    fPHOSTrackSegment = 0 ; 
  } 
}

//____________________________________________________________________________
Int_t AliPHOSRecParticle::DistancetoPrimitive(Int_t px, Int_t py)
{
  //  Compute distance from point px,py to a AliPHOSRecParticle considered as a Tmarker
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
 void AliPHOSRecParticle::Draw(Option_t *option)
 {
   // Draw this AliPHOSRecParticle with its current attributes
    
   AppendPad(option);
 }

//______________________________________________________________________________
void AliPHOSRecParticle::ExecuteEvent(Int_t event, Int_t px, Int_t py)
{
  //  Execute action corresponding to one event
  //  This member function is called when a AliPHOSRecParticle is clicked with the locator
  //
    
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
TString AliPHOSRecParticle::Name()
{
  TString  name ; 
  switch (fType) {
  case kGAMMA:
    name = "PHOTON" ;
    break ; 
   case kELECTRON:
     name = "ELECTRON" ;
    break ; 
  case kNEUTRAL:
    name = "NEUTRAL" ;
    break ; 
   case kCHARGEDHADRON:
    name = "CHARGED HADRON" ;
    break ; 
  case kNEUTRALHADRON:
    name = "NEUTRAL HADRON" ; 
    break ; 
  case kNEUTRALEM:
    name = "NEUTRAL EM" ; 
    break ; 
  case kGAMMAHADRON:
    name = "PHOTON HADRON" ; 
    break ; 

  }
  return name ; 
}

//______________________________________________________________________________
void AliPHOSRecParticle::Paint(Option_t *)
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
void AliPHOSRecParticle::Print()
{
  cout << "AliPHOSRecParticle > " << "type is  " << Name() << endl 
       << "                     " << "Energy = " << fE << endl 
       << "                     " << "Px     = " << fPx << endl 
       << "                     " << "Py     = " << fPy << endl 
       << "                     " << "Pz     = " << fPz << endl ; 
}
