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
// Macros analyzing the ESD file
// Use Case : 
//          root> .L AnaESD.C++
//          root> ana() --> prints the objects stored in ESD
//                                              
// author  : Yves Schutz (CERN/SUBATECH)
// February 2004
//_________________________________________________________________________
#if !defined(__CINT__) || defined(__MAKECINT__)
#include "TFile.h"
#include "TMath.h"
#include "TH1D.h"
#include "AliPHOSGetter.h"
#include "AliPHOSGeometry.h"
#include "Riostream.h"
#include "AliESD.h"
#include "AliESDtrack.h"
#include "AliESDCaloTrack.h"
#include "AliEMCALRecParticle.h"
#include "AliPHOSRecParticle.h"
#include "AliKalmanTrack.h"
#endif

Double_t Match(TParticle * pp, AliESDtrack * cp) ; 
TH1D * heta = new TH1D("heta", "Eta correlation", 100, -2., 2.) ; 
TH1D * hphi = new TH1D("hphi", "Phi correlation", 360, 0., 360.) ; 
 
void Ana() 
{ 
  AliPHOSGetter * gime = AliPHOSGetter::Instance("galice.root") ; 
  Int_t nEvent = gime->MaxEvent() ;  
  Int_t event ; 
  AliESD * esd = 0 ;
  for (event = 0 ; event < nEvent; event++) {
    esd = gime->ESD(event) ; 
    //esd->Print();  
    Int_t caloindex ;
    // Calorimeter tracks 
    AliESDCaloTrack * ct ; 
    for (caloindex = 0 ; caloindex < esd->GetNumberOfCaloTracks() ; caloindex++) {
      // get the calorimeter type of particles (PHOS or EMCAL)
      ct = esd->GetCaloTrack(caloindex) ;
      TParticle * part = ct->GetRecParticle() ; 

      AliESDtrack * cp ; 
      Int_t cpindex ; 
      for (cpindex = 0 ; cpindex < esd->GetNumberOfTracks() ; cpindex++) {
	// get the charged tracks from central tracking
	cp = esd->GetTrack(cpindex) ;
	Double_t dist = Match(part, cp) ; 
	
	if (dist < 99999.) 
	  cout << "================ Distance = " << dist << endl ; 
      }
    }
  }
  //  heta->Draw() ; 
  hphi->Draw() ; 
}
Double_t Match(TParticle * part, AliESDtrack * cp) 
{
  // Calculates the distance (x,z) between  the particle detected by PHOS and 
  // the charged particle reconstructed by the global tracking 

  Double_t dist = 99999. ;
   
  AliPHOSRecParticle * pp  = dynamic_cast<AliPHOSRecParticle*>(part) ;  
  AliEMCALRecParticle * ep = dynamic_cast<AliEMCALRecParticle*>(part) ;  
   
  Int_t phN ; 
  Double_t phZ, phX ; 
  
  if (pp) { // it is a PHOS particle 
    
 //    cout << "PHOS particle # " << " pos (" 
// 	 << pp->GetPos().X() << ", " << pp->GetPos().Y() << ", " << pp->GetPos().Z() << ")" << endl ;
    
    AliPHOSGetter * gime = AliPHOSGetter::Instance() ; 
    gime->PHOSGeometry()->ImpactOnEmc(*pp, phN, phZ, phX) ; 
    Double_t xyzAtPHOS[3] ; 
    cp->GetOuterXYZ(xyzAtPHOS) ; 
    if ( (xyzAtPHOS[0] +  xyzAtPHOS[1] + xyzAtPHOS[2]) != 0.) { //it has reached PHOS
      //the next check are only if we want high quality tracks 
      //       ULong_t status = cp->GetStatus() ;  
      //       if ((status & AliESDtrack::kTRDput)==0) 
      // 	do not continue;
      //       if ((status & AliESDtrack::kTRDStop)!=0) 
      // 	do not continue;  
//       cout << "Charged particle # " << " pos (" 
// 	   << xyzAtPHOS[0] << ", " << xyzAtPHOS[1] << ", " << xyzAtPHOS[2] << ")" <<  endl ;     
      TVector3 poscp(xyzAtPHOS[0], xyzAtPHOS[1], xyzAtPHOS[2]) ;
      Int_t cpN ;
      Double_t cpZ,cpX ; 
      gime->PHOSGeometry()->ImpactOnEmc(poscp, cpN, cpZ, cpX) ; 
      if (cpN) {// we are inside the PHOS acceptance 
// 	cout << "Charged Matching 1: " << cpN << " " << cpZ << " " << cpX << endl ; 
// 	cout << "Charged Matching 2: " << phN << " " << phZ << " " << phX << endl ; 
	dist = TMath::Sqrt( (cpZ-phZ)*(cpZ-phZ) + (cpX-phX)*(cpX-phX)) ;  
      } 
      Double_t phTheta = pp->Theta() ; 
      Double_t phPhi   = pp->Phi() ;
      TParticle tempo ; 
      tempo.SetMomentum(xyzAtPHOS[0], xyzAtPHOS[1], xyzAtPHOS[2], 0.) ;  
      Double_t cpTheta = tempo.Theta() ; 
      Double_t cpPhi   = tempo.Phi() ;
      //cout << phTheta << " " << phPhi << " " << endl 
      //cout <<	 cpTheta << " " << cpPhi-phPhi << " " << endl ; 
      heta->Fill( (phTheta - cpTheta)*TMath::RadToDeg() ) ; 
      hphi->Fill( (phPhi - cpPhi)*TMath::RadToDeg() ) ; 
    }
  }
  
  if (ep) {
    //cout << "EMCAL particle # " << endl ; 
  }
  return dist ; 
}
 
