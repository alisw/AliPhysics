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
#include "TFile.h"
#include "TMath.h"
#include "AliPHOSGetter.h"
#include "Riostream.h"
#include "AliESD.h"
#include "AliESDtrack.h"
#include "AliESDCaloTrack.h"
#include "AliEMCALRecParticle.h"
#include "AliPHOSRecParticle.h"
#include "AliKalmanTrack.h"

void Ana() 
{
  AliPHOSGetter * gime = AliPHOSGetter::Instance("galice.root") ; 
  Int_t nEvent = gime->MaxEvent() ;  
  Int_t event ; 
  AliESD * esd = 0 ;
  for (event = 0 ; event < nEvent; event++) {
    esd = gime->ESD(event) ; 
    esd->Print();  
    Int_t index ;
    // Calorimeter tracks 
    AliESDCaloTrack * ct ; 
    AliPHOSRecParticle * pp ; 
    AliEMCALRecParticle * ep ; 
    for (index = 0 ; index < esd->GetNumberOfCaloTracks() ; index++) {
      ct = esd->GetCaloTrack(index) ;
      pp = dynamic_cast<AliPHOSRecParticle*>(ct->GetRecParticle()) ; 
      ep = dynamic_cast<AliEMCALRecParticle*>(ct->GetRecParticle()) ; 
      if (pp) { 
	TVector3 pos = pp->GetPos() ;
	cout << "PHOS particle # " << index << " pos " << 
	  TMath::Sqrt(pos.X()*pos.X() + pos.Y()*pos.Y() + pos.Z()*pos.Z() ) << endl ; 
      }
      if(ep) { 
 	TVector3 pos = ep->GetPos() ;
	cout << "EMCAL particle # " << index << " pos " << 
	  TMath::Sqrt(pos.X()*pos.X() + pos.Y()*pos.Y() + pos.Z()*pos.Z() ) << endl ; 
      }
    }
    //Charged tracks from central tracking
    AliESDtrack * cp ; 
    for (index = 0 ; index < esd->GetNumberOfTracks() ; index++) {
      cp = esd->GetTrack(index) ;
      ULong_t status = cp->GetStatus() ; 

      // check if the tracks comes out of TRD
      if ((status & AliESDtrack::kTRDout)==0) 
	continue;
      if ((status & AliESDtrack::kTRDStop)!=0) 
	continue;

      // Gets the Global coordinate of the track at the entrance of PHOS 
      Double_t xyzAtPHOS[3] ; 
      cp->GetOuterXYZ(xyzAtPHOS) ; 
      // cout << xyzAtPHOS[0] << " " << xyzAtPHOS[1] << " "<< xyzAtPHOS[2]  << endl ; 
      //cout << TMath::Sqrt(xyzAtPHOS[0]*xyzAtPHOS[0] + xyzAtPHOS[1]*xyzAtPHOS[1] + xyzAtPHOS[2]*xyzAtPHOS[2]) << endl ;  
     
      // Does the matching with PHOS/EMCAL
//       for (index = 0 ; index < esd->GetNumberOfCaloTracks() ; index++) {
// 	ct = esd->GetCaloTrack(index) ;
// 	pp = dynamic_cast<AliPHOSRecParticle*>(ct->GetRecParticle()) ; 
// 	ep = dynamic_cast<AliEMCALRecParticle*>(ct->GetRecParticle()) ; 
// 	if (pp) {
// 	  TVector3 pos = pp->GetPos() ; 
// 	}
//       }
    }
  }
}
