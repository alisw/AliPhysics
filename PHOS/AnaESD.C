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
#include "AliPHOSGetter.h"
#include "Riostream.h"
#include "AliESD.h"
#include "AliESDCaloTrack.h"
#include "AliEMCALRecParticle.h"
#include "AliPHOSRecParticle.h"

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
    AliESDCaloTrack * ct ; 
    AliPHOSRecParticle * pp ; 
    AliEMCALRecParticle * ep ;     
    for (index = 0 ; index < esd->GetNumberOfCaloTracks() ; index++) {
      ct = esd->GetCaloTrack(index) ;
      pp = dynamic_cast<AliPHOSRecParticle*>(ct->GetRecParticle()) ; 
      ep = dynamic_cast<AliEMCALRecParticle*>(ct->GetRecParticle()) ; 
      if (pp) 
	cout << "particle # " << index << " is from PHOS " << endl ; 
      if(ep) 
	cout << "particle # " << index << " is from EMCAL " << endl ; 
    }
  }
}
