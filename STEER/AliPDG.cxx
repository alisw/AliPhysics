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

// Class to encapsulate the ALICE updates to TDatabasePDG.h
// Can be used by TGeant3 and TGeant4
// It contains also the constants for the PDG particle IDs.
// Should evolve towards dynamical loading from external data base.
// Comments to: andreas.morsch@cern.ch 

#include "AliPDG.h"
#include "TDatabasePDG.h"

ClassImp(AliPDG)



void AliPDG::AddParticlesToPdgDataBase()
{

//
// Add particles to the PDG data base
//
  
  static Bool_t bAdded = kFALSE;
  // Check if already called
  if(bAdded)return;
  bAdded = true;

  TDatabasePDG *pdgDB = TDatabasePDG::Instance();
  const Int_t kspe=50000000;

  // PDG nuclear states are 10-digit numbers
  // 10LZZZAAAI e.g. deuteron is 
  // 1000010020
  const Int_t kion=1000000000;

  const Double_t kAu2Gev=0.9314943228;
  const Double_t khSlash = 1.0545726663e-27;
  const Double_t kErg2Gev = 1/1.6021773349e-3;
  const Double_t khShGev = khSlash*kErg2Gev;
  const Double_t kYear2Sec = 3600*24*365.25;


//
// Bottom mesons
// mass and life-time from PDG
//
  pdgDB->AddParticle("Upsilon(3S)","Upsilon(3S)",10.3552,kTRUE,
                     0,1,"Bottonium",200553);

// QCD diffractive states
  pdgDB->AddParticle("rho_diff0","rho_diff0",0,kTRUE,
		     0,0,"QCD diffr. state",9900110);
  pdgDB->AddParticle("pi_diffr+","pi_diffr+",0,kTRUE,
		     0,1,"QCD diffr. state",9900210);
  pdgDB->AddParticle("omega_di","omega_di",0,kTRUE,
		     0,0,"QCD diffr. state",9900220);
  pdgDB->AddParticle("phi_diff","phi_diff",0,kTRUE,
		     0,0,"QCD diffr. state",9900330);
  pdgDB->AddParticle("J/psi_di","J/psi_di",0,kTRUE,
		     0,0,"QCD diffr. state",9900440);
  pdgDB->AddParticle("n_diffr0","n_diffr0",0,kTRUE,
		     0,0,"QCD diffr. state",9902110);
  pdgDB->AddParticle("p_diffr+","p_diffr+",0,kTRUE,
		     0,1,"QCD diffr. state",9902210);

  // Some particles produced by HERWIG
  pdgDB->AddParticle("rho3(1690)0","rho(1690)0", 1.69, kTRUE,
		     0, 0,"rho", 117);
  pdgDB->AddParticle("pion2(1670)0","pion2(1670)0", 1.67, kTRUE,
		     0, 0,"pion", 10115);
  pdgDB->AddParticle("omega(1650)","omega(1650)", 1.65, kTRUE,
		     0, 0,"pion", 30223);

  // IONS
  //
  // Done by default now from Pythia6 table
  // Needed for other generators
  // So check if already defined


  Int_t ionCode = kion+10020;
  if(!pdgDB->GetParticle(ionCode)){
      pdgDB->AddParticle("Deuteron","Deuteron",2*kAu2Gev+8.071e-3,kTRUE,
			 0,1,"Ion",ionCode);
  }

  ionCode = kion+10030;
  if(!pdgDB->GetParticle(ionCode)){
    pdgDB->AddParticle("Triton","Triton",3*kAu2Gev+14.931e-3,kFALSE,
                     khShGev/(12.33*kYear2Sec),1,"Ion",ionCode);
  }

  ionCode = kion+20030;
  if(!pdgDB->GetParticle(ionCode)){
    pdgDB->AddParticle("HE3","HE3",3*kAu2Gev+14.931e-3,kFALSE,
                     0,2,"Ion",ionCode);
  }

  ionCode = kion+20040;
  if(!pdgDB->GetParticle(ionCode)){
    pdgDB->AddParticle("Alpha","Alpha",4*kAu2Gev+2.424e-3,kTRUE,
		       khShGev/(12.33*kYear2Sec),2,"Ion",ionCode);
  }


// Special particles
// 
  pdgDB->AddParticle("Cherenkov","Cherenkov",0,kFALSE,
		     0,0,"Special",kspe+50);
  pdgDB->AddParticle("FeedbackPhoton","FeedbackPhoton",0,kFALSE,
		     0,0,"Special",kspe+51);
  pdgDB->AddParticle("Lambda1520","Lambda1520",1.5195,kFALSE,
		     0.0156,0,"Resonance",3124);
  pdgDB->AddAntiParticle("Lambda1520bar",-3124);
}


