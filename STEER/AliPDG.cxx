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

/*
$Log$
Revision 1.2  2001/01/31 14:32:42  morsch
Some B mesons added

Revision 1.1  2000/12/21 16:48:39  morsch
AliPDG class, first commit.

*/

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

    TDatabasePDG *pdgDB = TDatabasePDG::Instance();

    const Int_t kion=10000000;
    const Int_t kspe=50000000;

    const Double_t kAu2Gev=0.9314943228;
    const Double_t khSlash = 1.0545726663e-27;
    const Double_t kErg2Gev = 1/1.6021773349e-3;
    const Double_t khShGev = khSlash*kErg2Gev;
    const Double_t kYear2Sec = 3600*24*365.25;
//
// Bottom mesons
// mass and life-time from PDG
//
// Done by default now from Pythia6 table!
//
//
// Ions 
//

  pdgDB->AddParticle("Deuteron","Deuteron",2*kAu2Gev+8.071e-3,kTRUE,
                     0,1,"Ion",kion+10020);
  pdgDB->AddParticle("Triton","Triton",3*kAu2Gev+14.931e-3,kFALSE,
                     khShGev/(12.33*kYear2Sec),1,"Ion",kion+10030);
  pdgDB->AddParticle("Alpha","Alpha",4*kAu2Gev+2.424e-3,kTRUE,
                     khShGev/(12.33*kYear2Sec),2,"Ion",kion+20040);
  pdgDB->AddParticle("HE3","HE3",3*kAu2Gev+14.931e-3,kFALSE,
                     0,2,"Ion",kion+20030);
// Special particles
// 
  pdgDB->AddParticle("Cherenkov","Cherenkov",0,kFALSE,
		     0,0,"Special",kspe+50);
  pdgDB->AddParticle("FeedbackPhoton","FeedbackPhoton",0,kFALSE,
		     0,0,"Special",kspe+51);

}


