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
    pdgDB->AddParticle("B(s)*0","B(s)*0",
		       5.4163, kTRUE, 0.047, +0.,"Meson",  533);
    
    pdgDB->AddParticle("B(s)*0 bar","B(s)*0 bar",
		       5.4163, kTRUE, 0.047, -0.,"Meson", -533);
    
// Charmed baryons
// 
// value for mass used by Hijing
    pdgDB->AddParticle("Sigma(c)*+","Sigma(c)*+",
		       2.4536, kTRUE, -1., +1.,"Baryon",  4214);
    
    pdgDB->AddParticle("Sigma(c)*-","Sigma(c)*-",
		       2.4536, kTRUE, -1., -1.,"Baryon", -4214);
// equivalent to 4312 ? Hijing uses m=2.55
    pdgDB->AddParticle("Xsi(c)0","Xsi(c)0",
		       2.4703, kTRUE, -1., +0.,"Baryon",  4132);
    
    pdgDB->AddParticle("Xsi(c)0 bar","Xsi(c)0 bar",
		       2.4703, kTRUE, -1., -0.,"Baryon", -4132);
// equivalent to 4322 ? Hijing uses m=2.55
    pdgDB->AddParticle("Xi(c)+","Xi(c)+",
		       2.4656, kFALSE, -1., +1.,"Baryon",  4232);
    
    pdgDB->AddParticle("Xi(c)-","Xi(c)-",
		       2.4656, kFALSE, -1., -1.,"Baryon", -4232);
// mass values from Hijing
    
    pdgDB->AddParticle("Xsi(c)*0","Xsi(c)*0",
		       2.63, kTRUE, -1., +0.,"Baryon",  4314);
    
    pdgDB->AddParticle("Xsi(c)*0 bar","Xsi(c)*0 bar",
		       2.63, kTRUE, -1., -0.,"Baryon", -4314);
    
    pdgDB->AddParticle("Xsi(c)*+","Xsi(c)*+",
		       2.63, kTRUE, -1., +1.,"Baryon",  4324);
    
    pdgDB->AddParticle("Xsi(c)*-","Xsi(c)*-",
		       2.63, kTRUE, -1., -1.,"Baryon", -4324);
    
// pdg mass value, Hijing uses m=2.73.
    pdgDB->AddParticle("Omega(c)0","Omega(c)0",
		       2.7040, kFALSE, khShGev/0.064e-12, +0.,"Baryon",  4332);
    
    pdgDB->AddParticle("Omega(c)0 bar","Omega(c)0 bar",
		       2.7040, kFALSE, khShGev/0.064e-12, -0.,"Baryon", -4332);
// mass value from Hijing
    pdgDB->AddParticle("Omega(c)*0","Omega(c)*0",
		       2.8000, kFALSE, -1., +0.,"Baryon",  4334);
    
    pdgDB->AddParticle("Omega(c)*0 bar","Omega(c)*0",
		       2.8000, kFALSE, -1., -0.,"Baryon", -4334);
    
    
// Xi(cc)
    
    pdgDB->AddParticle("Xsi(cc)+","Xsi(cc)+",
		     3.60, kTRUE, -1., +1.,"Baryon",  4412);
    
    pdgDB->AddParticle("Xsi(cc) bar-","Xsi(cc) bar-",
		       3.60, kTRUE, -1., -1.,"Baryon", -4412);
    
    pdgDB->AddParticle("Xsi*(cc)+","Xsi*(cc)+",
		       3.66, kTRUE, -1., +1.,"Baryon",  4414);

    pdgDB->AddParticle("Xsi*(cc) bar-","Xsi*(cc) bar-",
		       3.66, kTRUE, -1., -1.,"Baryon", -4414);
    
    
    pdgDB->AddParticle("Xsi(cc)++","Xsi(cc)++",
		       3.60, kTRUE, -1., +2.,"Baryon",  4422);
    
    pdgDB->AddParticle("Xsi(cc) bar--","Xsi(cc) bar--",
		       3.60, kTRUE, -1., -2.,"Baryon", -4422);
    
    
    pdgDB->AddParticle("Xsi*(cc)++","Xsi*(cc)++",
		       3.66, kTRUE, -1., +2.,"Baryon",  4424);
    
    pdgDB->AddParticle("Xsi*(cc) bar-","Xsi*(cc) bar-",
		       3.66, kTRUE, -1., -2.,"Baryon", -4424);
    
    pdgDB->AddParticle("Omega(cc)+","Omega(cc)+",
		       3.78, kTRUE, -1., +1.,"Baryon",  4432);
    
    pdgDB->AddParticle("Omega(cc) bar-","Omega(cc) bar-",
		       3.78, kTRUE, -1., -1.,"Baryon", -4432);
    
    pdgDB->AddParticle("Omega*(cc)+","Omega*(cc)+",
		       3.82, kTRUE, -1., +1.,"Baryon",  4434);
    
    pdgDB->AddParticle("Omega*(cc) bar-","Omega*(cc) bar-",
		       3.82, kTRUE, -1., -1.,"Baryon", -4434);
    
    
    pdgDB->AddParticle("Omega*(ccc)+","Omega*(cc)++",
		       4.91, kTRUE, -1., +2.,"Baryon",  4444);
    
    pdgDB->AddParticle("Omega*(ccc) bar--","Omega*(cc) bar--",
		       4.91, kTRUE, -1., -2.,"Baryon", -4444);
    
    

// Bottom baryons
//
// mass value from Hijing
    pdgDB->AddParticle("Sigma(b)*+","Sigma(b)*+",
		       5.8100, kFALSE, -1., +1.,"Baryon", 5224);
    
    pdgDB->AddParticle("Sigma(b)*-","Sigma(b)*-",
		       5.8100, kFALSE, -1., -1.,"Baryon", -5224);
    
    
    pdgDB->AddParticle("Xi(b)0","Xi(b)0",
		       5.8400, kFALSE, -1., +0.,"Baryon", 5232);
    
    pdgDB->AddParticle("Xi(b)0 bar","Xi(b)0 bar",
		       5.8100, kFALSE, -1., -0.,"Baryon", -5232);
    
// B(s) 
    pdgDB->AddParticle("Xi'(b)-","Xi'(b)-",
		       5.9600, kFALSE, -1., -1.,"Baryon", 5312);
    
    pdgDB->AddParticle("Xi'(b) bar+","Xi'(b) bar+",
		       5.9600, kFALSE, -1.,  1.,"Baryon", -5312);
    
    pdgDB->AddParticle("Xi*(b)-","Xi*(b)-",
		       5.9700, kFALSE, -1., -1.,"Baryon", 5314);
    
    pdgDB->AddParticle("Xi*(b) bar+","Xi*(b) bar+",
		       5.9700, kFALSE, -1.,  1.,"Baryon", -5314);
    
    pdgDB->AddParticle("Xi'(b)0","Xi'(b)0",
		       5.9600, kFALSE, -1., -0.,"Baryon", 5322);
    
    pdgDB->AddParticle("Xi'(b) bar0","Xi'(b) bar0",
		       5.9600, kFALSE, -1.,  0.,"Baryon", -5322);
    
    pdgDB->AddParticle("Xi*(b)0","Xi*(b)0",
		       5.9700, kFALSE, -1., -0.,"Baryon", 5324);
    
    pdgDB->AddParticle("Xi*(b) bar0","Xi*(b) bar0",
		       5.9700, kFALSE, -1.,  0.,"Baryon", -5324);
    
    pdgDB->AddParticle("Omega(b)-","Omega(b)-",
		       6.1200, kFALSE, -1., -1.,"Baryon", 5332);
    
    pdgDB->AddParticle("Omega(b) bar+","Omega(b) bar+",
		       6.1200, kFALSE, -1.,  1.,"Baryon", -5332);
    
    pdgDB->AddParticle("Omega*(b)-","Omega*(b)-",
		       6.1300, kFALSE, -1., -1.,"Baryon", 5334);
    
    pdgDB->AddParticle("Omega*(b) bar+","Omega*(b) bar+",
		       6.1300, kFALSE, -1.,  1.,"Baryon", -5334);


    pdgDB->AddParticle("Omega*(b)-","Omega*(b)-",
		       6.1300, kFALSE, -1., -1.,"Baryon", 5334);
    
    pdgDB->AddParticle("Omega*(b) bar+","Omega*(b) bar+",
		       6.1300, kFALSE, -1.,  1.,"Baryon", -5334);

// B(c) 

    pdgDB->AddParticle("Omega(bc)0","Omega(bc)0",
		       7.1900, kFALSE, -1., -0.,"Baryon", 5342);
    
    pdgDB->AddParticle("Omega(bc) bar0","Omega(bc) bar0",
		       7.1900, kFALSE, -1.,  0.,"Baryon", -5342);
    
    pdgDB->AddParticle("Xi'(bc)0","Xi'(bc)0",
		       7.0400, kFALSE, -1., -0.,"Baryon", 5412);
    
    pdgDB->AddParticle("Xi'(bc) bar0","Xi'(bc) bar0",
		       7.0400, kFALSE, -1.,  0.,"Baryon", -5412);
    
    pdgDB->AddParticle("Xi*(bc)0","Xi*(bc)0",
		       7.0500, kFALSE, -1., -0.,"Baryon", 5414);
    
    pdgDB->AddParticle("Xi*(bc) bar0","Xi*(bc) bar0",
		       7.0500, kFALSE, -1.,  0.,"Baryon", -5414);
    
    pdgDB->AddParticle("Xi'(bc)+","Xi'(bc)+",
		       7.0400, kFALSE, -1., +1.,"Baryon", 5422);
    
    pdgDB->AddParticle("Xi'(bc) bar-","Xi'(bc) bar-",
		       7.0400, kFALSE, -1., -1.,"Baryon", -5422);
    
    pdgDB->AddParticle("Xi*(bc)+","Xi*(bc)+",
		       7.0500, kFALSE, -1., +1.,"Baryon", 5424);
    
    pdgDB->AddParticle("Xi*(bc) bar-","Xi*(bc) bar-",
		       7.0500, kFALSE, -1., -1.,"Baryon", -5424);
    
    pdgDB->AddParticle("Omega'(bc)0","Omega'(bc)0",
		       7.2100, kFALSE, -1., -0.,"Baryon", 5432);
    
    pdgDB->AddParticle("Omega'(bc) bar0","Omega'(bc) bar0",
		       7.2100, kFALSE, -1.,  0.,"Baryon", -5432);
    
    pdgDB->AddParticle("Omega*(bc)0","Omega*(bc)0",
		       7.2200, kFALSE, -1., -0.,"Baryon", 5434);
    
    pdgDB->AddParticle("Omega*(bc) bar0","Omega*(bc) bar0",
		       7.2200, kFALSE, -1.,  0.,"Baryon", -5434);
// B(bcc)
    pdgDB->AddParticle("Omega(bcc)+","Omega(bcc)+",
		       8.3100, kFALSE, -1., +1.,"Baryon", 5442);
    
    pdgDB->AddParticle("Omega(bcc) bar-","Omega(bcc) bar-",
		       8.3100, kFALSE, -1., -1.,"Baryon", -5442);
    
    pdgDB->AddParticle("Omega*(bcc)+","Omega*(bcc)+",
		       8.3100, kFALSE, -1., +1.,"Baryon", 5444);
    
    pdgDB->AddParticle("Omega*(bcc) bar-","Omega*(bcc) bar-",
		       8.3100, kFALSE, -1., -1.,"Baryon", -5444);




// B(bb)

    pdgDB->AddParticle("Xsi(bb)-","Xsi(bb)-",
		       10.4200, kFALSE, -1., -1.,"Baryon", 5512);
    
    pdgDB->AddParticle("Xsi(bb) bar+","Xsi(bb) bar+",
		       10.4200, kFALSE, -1., +1.,"Baryon", -5512);
    
    pdgDB->AddParticle("Xsi*(bb)-","Xsi*(bb)-",
		       10.4400, kFALSE, -1., -1.,"Baryon", 5514);
    
    pdgDB->AddParticle("Xsi*(bb) bar+","Xsi*(bb) bar+",
		       10.4400, kFALSE, -1., +1.,"Baryon", -5514);
    
    pdgDB->AddParticle("Xsi(bb)0","Xsi(bb)0",
		       10.4200, kFALSE, -1., -0.,"Baryon", 5522);
    
    pdgDB->AddParticle("Xsi(bb) bar0","Xsi(bb) bar0",
		       10.4200, kFALSE, -1., +0.,"Baryon", -5522);
    
    pdgDB->AddParticle("Xsi*(bb)0","Xsi*(bb)0",
		       10.4400, kFALSE, -1., -0.,"Baryon", 5524);
    
    pdgDB->AddParticle("Xsi*(bb) bar0","Xsi*(bb) bar0",
		       10.4400, kFALSE, -1., +0.,"Baryon", -5524);
    
    pdgDB->AddParticle("Omega*(bb)-","Omega(bb)-",
		       10.6000, kFALSE, -1., -1.,"Baryon", 5532);
    
    pdgDB->AddParticle("Omega(bb) bar+","Omega(bb) bar+",
		       10.6000, kFALSE, -1., +1.,"Baryon", -5532);
    
    pdgDB->AddParticle("Omega*(bb)-","Omega*(bb)-",
		       10.6000, kFALSE, -1., -1.,"Baryon", 5534);
    
    pdgDB->AddParticle("Omega*(bb) bar+","Omega*(bb) bar+",
		       10.6000, kFALSE, -1., +1.,"Baryon", -5534);
    
// B(bbc)
    
    pdgDB->AddParticle("Omega(bbc)0","Omega(bbc)0",
		       11.7100, kFALSE, -1., -0.,"Baryon", 5542);
    
    pdgDB->AddParticle("Omega(bbc) bar0","Omega(bbc) bar0",
		       11.7100, kFALSE, -1., +0.,"Baryon", -5542);
    
    pdgDB->AddParticle("Omega*(bbc)0","Omega*(bbc)0",
		       11.7100, kFALSE, -1., -0.,"Baryon", 5544);

    pdgDB->AddParticle("Omega*(bbc) bar0","Omega*(bbc) bar0",
		     11.7100, kFALSE, -1., +0.,"Baryon", -5544);
// B(bbb)
    
    pdgDB->AddParticle("Omega*(bbb)-","Omega*(bbb)-",
		       15.1000, kFALSE, -1., -1.,"Baryon", 5544);

    pdgDB->AddParticle("Omega*(bbb) bar+","Omega*(bbb) bar+",
		       15.100, kFALSE, -1., +1.,"Baryon", -5544);
    
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

//
// special particles
// 
  pdgDB->AddParticle("Cherenkov","Cherenkov",0,kFALSE,
		     0,0,"Special",kspe+50);
  pdgDB->AddParticle("FeedbackPhoton","FeedbackPhoton",0,kFALSE,
		     0,0,"Special",kspe+51);
}
