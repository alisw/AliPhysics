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

ClassImp(AliEMCALFastRecParticle) 

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
Bool_t AliEMCALFastRecParticle::IsPhoton(TString purity) const
{
  //FIXME : just a blind copy from PHOS
  // Rec.Particle is a photon if it has a photon-like shape, fast and neutral
  // photon-like shape is defined with a purity "low", "medium" or "high"

  purity.ToLower();
  Bool_t photonLike = kFALSE;
  if      (purity == "low"   ) photonLike = TestPIDBit(6);
  else if (purity == "medium") photonLike = TestPIDBit(7);
  else if (purity == "high"  ) photonLike = TestPIDBit(8);
  if (photonLike                                   && //  photon by PCA
      (TestPIDBit(5)||TestPIDBit(4)||TestPIDBit(3))&& //  fast by TOF
      (TestPIDBit(2)||TestPIDBit(1)||TestPIDBit(0))&& //  neutral by CPV
      !TestPIDBit(14))                              //  no charged track
    return kTRUE ;
  else
    return kFALSE;
}

//____________________________________________________________________________
Bool_t AliEMCALFastRecParticle::IsPi0(TString purity) const
{
  //FIXME : just a blind copy from PHOS
  // Rec.Particle is a pi0 if it has a pi0-like shape, fast and neutral
  // pi0-like shape is defined with a purity "low", "medium" or "high"

  purity.ToLower();
  Bool_t pi0Like = kFALSE;
  if      (purity == "low"   ) pi0Like = TestPIDBit(9);
  else if (purity == "medium") pi0Like = TestPIDBit(10);
  else if (purity == "high"  ) pi0Like = TestPIDBit(11);
  else Error("IsPi0","Wrong purity type: %s",purity.Data());
  if (pi0Like                                      && //  pi0 by PCA
      (TestPIDBit(5)||TestPIDBit(4)||TestPIDBit(3))&& //  fast by TOF
      (TestPIDBit(2)||TestPIDBit(1)||TestPIDBit(0))&& //  neutral by CPV
      !TestPIDBit(14))                              //  no charged track
    return kTRUE ;
  else
    return kFALSE;
}

//____________________________________________________________________________
Bool_t AliEMCALFastRecParticle::IsElectron(TString purity) const
{
  //FIXME : just a blind copy from PHOS
  // Rec.Particle is an electron if it has a photon-like shape, fast and charged
  // photon-like shape is defined with a purity "low", "medium" or "high"

  purity.ToLower();
  Bool_t photonLike = kFALSE;
  if      (purity == "low"   ) photonLike = TestPIDBit(6);
  else if (purity == "medium") photonLike = TestPIDBit(7);
  else if (purity == "high"  ) photonLike = TestPIDBit(8);
  else Error("IsElectron","Wrong purity type: %s",purity.Data());
  
  if (photonLike                                   && //  photon by PCA
      (TestPIDBit(5)|| TestPIDBit(4)|| TestPIDBit(3))&& //  fast by TOF
      (!TestPIDBit(2)||!TestPIDBit(1)||!TestPIDBit(0))&& //  charged by CPV
      TestPIDBit(14))                                  //  no charged track
    return kTRUE ;
  else
    return kFALSE;
}

//____________________________________________________________________________
Bool_t AliEMCALFastRecParticle::IsEleCon(TString purity) const
{
  //FIXME : just a blind copy from PHOS
  // Rec.Particle is an electron if it has a photon-like shape, fast and charged
  // photon-like shape is defined with a purity "low", "medium" or "high"

  purity.ToLower();
  Bool_t photonLike = kFALSE;
  if      (purity == "low"   ) photonLike = TestPIDBit(6);
  else if (purity == "medium") photonLike = TestPIDBit(7);
  else if (purity == "high"  ) photonLike = TestPIDBit(8);
  else Error("IsElectron","Wrong purity type: %s",purity.Data());
  
  if (photonLike                                   && //  photon by PCA
      (TestPIDBit(5)|| TestPIDBit(4)|| TestPIDBit(3))&& //  fast by TOF
      (!TestPIDBit(2)||!TestPIDBit(1)||!TestPIDBit(0))&& //  charged by CPV
      !TestPIDBit(14))                                  //  no charged track
    return kTRUE ;
  else
    return kFALSE;
}

//____________________________________________________________________________
Bool_t AliEMCALFastRecParticle::IsHardPhoton() const
{
  //FIXME : just a blind copy from PHOS
  // Rec.Particle is a hard photon (E > 30 GeV) if its second moment M2x
  // corresponds to photons
  if (TestPIDBit(12) && !TestPIDBit(14))
    return kTRUE;
  else
    return kFALSE;
}

//____________________________________________________________________________
Bool_t AliEMCALFastRecParticle::IsHardPi0() const
{
  //FIXME : just a blind copy from EMCAL
  // Rec.Particle is a hard pi0 (E > 30 GeV) if its second moment M2x
  // corresponds to pi0
  if (TestPIDBit(13)&& !TestPIDBit(14))
    return kTRUE;
  else
    return kFALSE;
}

//____________________________________________________________________________
Bool_t AliEMCALFastRecParticle::IsHadron() const
{
  //FIXME : just a blind copy from EMCAL
  // Rec.Particle is an hadron if it does not look like
  // a low-purity photon nor low-purity pi0

  if ( !TestPIDBit(6) && !TestPIDBit(9) )             // not photon nor pi0
    return kTRUE ;
  else
    return kFALSE;
}

//____________________________________________________________________________
Bool_t AliEMCALFastRecParticle::IsChargedHadron() const
{
  //FIXME : just a blind copy from EMCAL
  // Rec.Particle is a charged hadron if it does not look like
  // a low-purity photon nor low-purity pi0 and is low-purity charged

  if ( !TestPIDBit(6) && !TestPIDBit(9) &&            // not photon nor pi0
       !TestPIDBit(2))                                // charged by CPV
    return kTRUE ;
  else
    return kFALSE;
}

//____________________________________________________________________________
Bool_t AliEMCALFastRecParticle::IsNeutralHadron() const
{
  //FIXME : just a blind copy from EMCAL
  // Rec.Particle is a neutral hadron if it does not look like
  // a low-purity photon nor low-purity pi0 and is high-purity neutral

  if ( !TestPIDBit(6) && !TestPIDBit(9) &&            // not photon nor pi0
        TestPIDBit(2))                                // neutral by CPV
    return kTRUE ;
  else
    return kFALSE;
}

//____________________________________________________________________________
Bool_t AliEMCALFastRecParticle::IsFastChargedHadron() const
{
  //FIXME : just a blind copy from EMCAL
  // Rec.Particle is a fast charged hadron if it does not look like
  // a low-purity photon nor low-purity pi0, is low-purity charged
  // and is high-purity fast

  if ( !TestPIDBit(6) && !TestPIDBit(9) &&            // not photon nor pi0
       !TestPIDBit(2) &&                              // charged by CPV
        TestPIDBit(5))                                // fast by TOF
    return kTRUE ;
  else
    return kFALSE;
}

//____________________________________________________________________________
Bool_t AliEMCALFastRecParticle::IsSlowChargedHadron() const
{
  //FIXME : just a blind copy from EMCAL
  // Rec.Particle is a slow neutral hadron if it does not look like
  // a low-purity photon nor low-purity pi0, is high-purity neutral
  // and is not high-purity fast

  if ( !TestPIDBit(6) && !TestPIDBit(9) &&            // not photon nor pi0
       !TestPIDBit(2) &&                              // charged by CPV
       !TestPIDBit(5))                                // slow by TOF
    return kTRUE ;
  else
    return kFALSE;

}

//____________________________________________________________________________
Bool_t AliEMCALFastRecParticle::IsFastNeutralHadron() const
{
  //FIXME : just a blind copy from EMCAL
  // Rec.Particle is a fast neutral hadron if it does not look like
  // a low-purity photon nor low-purity pi0, is high-purity neutral
  // and is high-purity fast

  if ( !TestPIDBit(6) && !TestPIDBit(9) &&            // not photon nor pi0
        TestPIDBit(2) &&                              // neutral by CPV
        TestPIDBit(5))                                // fast by TOF
    return kTRUE ;
  else
    return kFALSE;
}

//____________________________________________________________________________
Bool_t AliEMCALFastRecParticle::IsSlowNeutralHadron() const
{
  //FIXME : just a blind copy from EMCAL
  // Rec.Particle is a slow neutral hadron if it does not look like
  // a low-purity photon nor low-purity pi0, is high-purity neutral
  // and is not high-purity fast

  if ( !TestPIDBit(6) && !TestPIDBit(9) &&            // not photon nor pi0
        TestPIDBit(2) &&                              // neutral by CPV
       !TestPIDBit(5))                                // slow by TOF
    return kTRUE ;
  else
    return kFALSE;
}


//____________________________________________________________________________
TString AliEMCALFastRecParticle::Name()const
{
   // Returns the name of the particle type (only valid if PIDv1 is employed)

  TString  name ; 

  name = "Undefined particle" ;
  
  if      (IsPhoton("low"))
    name = "Photon low purity, ";
  else if (IsPhoton("medium"))
    name = "Photon medium purity, ";
  else if (IsPhoton("high"))
    name = "Photon high purity, ";

  if      (IsPi0("low"))
    name += "Pi0 low purity, ";
  else if (IsPi0("medium"))
    name += "Pi0 medium purity, ";
  else if (IsPi0("high"))
    name += "Pi0 high purity, ";

  if      (IsElectron("low"))
    name += "Electron low purity, ";
  else if (IsElectron("medium"))
    name += "Electron medium purity, ";
  else if (IsElectron("high"))
    name += "Electron high purity, ";

  if     (IsHadron()) {
    name = "hadron";
    if      (IsChargedHadron()) {
      name.Prepend("charged, ");
      if      (IsFastChargedHadron())
	name.Prepend("fast, ");
      else if (IsSlowChargedHadron())
	name.Prepend("slow, ");
    }
    else if (IsNeutralHadron()) {
      name.Prepend("neutral, ");
      if      (IsFastNeutralHadron())
	name.Prepend("fast, ");
      else if (IsSlowNeutralHadron())
	name.Prepend("slow, ");
    }
  }

  return name ; 
}

//______________________________________________________________________________
void AliEMCALFastRecParticle::SetType(Int_t type) { 
  // sets the particle type 
  // bit-mask of the particle type means the following:
  // bits 0,1,2   - neutral particle with low, medium and high purity
  // bits 3.4,5   - fast particle with low, medium and high purity
  // bits 6.7,8   - photon shower with low, medium and high purity
  // bits 9,10,11 - hard-pi0 shower with low, medium and high purity

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
