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
#include "AliPHOSGridFile.h"
#endif

void Match(TParticle * pp, AliESDtrack * cp, Double_t * dist) ; 
TH1D * heta = new TH1D("heta", "Eta correlation", 100, 0., 360.) ; 
TH1D * hphi = new TH1D("hphi", "Phi correlation", 100, 0., 360.) ; 
// Data Challenge identification
const TString kYear("2004") ; 
const TString kProd("02") ; 
const TString kVers("V4.01.Rev.00") ; 

Bool_t Ana(const TString type = "per5", const Int_t run = 1, const Int_t nOfEvt = 1) 
{ 
  Double_t dist[3] ; 
  // get the LFN file name in the Grid catalogue ; 
  AliPHOSGridFile lfn ; 
  if (!lfn.IsConnected()) 
    return kFALSE ; 
  lfn.SetPath(kYear, kProd, kVers, type) ;  
  lfn.SetRun(run) ; 

  //loop over the events 
  Int_t nevt, evt = 0 ; 
  for (nevt = 0 ; nevt < nOfEvt ; nevt++) { 
    evt++ ; 
    lfn.SetEvt(evt) ;
    TString fileName = lfn.GetLFN() ; 
    
    if (fileName.IsNull()) {
      nevt-- ; 
      continue ; 
    }

    printf(">>>>>>>>>>>> Processing %s-%s/%s/%s : run # %d event # %d \n", 
	   kYear.Data(), kProd.Data(), kVers.Data(), type.Data(), run, evt) ;  
    AliPHOSGetter * gime = AliPHOSGetter::Instance(fileName) ; 
    
    Int_t nEvent = gime->MaxEvent() ;  
    Int_t event ; 
    AliESD * esd = 0 ;
    for (event = 0 ; event < nEvent; event++) {
      esd = gime->ESD(event) ; 
      if (!esd) 
	return kFALSE ; 
      
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
	  Match(part, cp, dist) ; 
	}
	heta->Fill( dist[1] ) ; 
	hphi->Fill( dist[2] ) ; 
      }
    }
  }
  heta->Draw() ; 
  //hphi->Draw() ; 
  return kTRUE ; 
}

void Match(TParticle * part, AliESDtrack * cp, Double_t * dist) 
{
  // Calculates the distance (x,z) between  the particle detected by PHOS and 
  // the charged particle reconstructed by the global tracking 

   
  AliPHOSRecParticle * pp  = dynamic_cast<AliPHOSRecParticle*>(part) ;  
  AliEMCALRecParticle * ep = dynamic_cast<AliEMCALRecParticle*>(part) ;  
   
  Int_t phN ; 
  Double_t phZ, phX ; 
  
  if (pp) { // it is a PHOS particle 
    Double_t cpTheta,  cpPhi ;  
    Double_t phTheta,  phPhi ; 
    cpTheta = cpPhi = phTheta = phPhi = 0. ; 
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
	dist[0] = TMath::Sqrt( (cpZ-phZ)*(cpZ-phZ) + (cpX-phX)*(cpX-phX)) ;  
      } 
      phTheta = pp->Theta() ; 
      phPhi   = pp->Phi() ;
      TParticle tempo ; 
      tempo.SetMomentum(xyzAtPHOS[0], xyzAtPHOS[1], xyzAtPHOS[2], 0.) ;  
      cpTheta = tempo.Theta() ; 
      cpPhi   = tempo.Phi() ;
      //cout << phTheta << " " << phPhi << " " << endl 
      //cout <<	 cpTheta << " " << cpPhi-phPhi << " " << endl ; 
    }
    dist[1] = (phTheta - cpTheta)*TMath::RadToDeg() ; 
    dist[2] = (phPhi - cpPhi)*TMath::RadToDeg() ; 
  }
  
  if (ep) {
    //cout << "EMCAL particle # " << endl ; 
  }
}
