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
//
// Realisation of AliGenReader to be used with AliGenExtFile
// It reads Pythia Jet events from a ntuple like event structure.
// The event format is defined in Init()
// NextEvent() is used to loop over events and NextParticle() to loop over particles.  
// Author: andreas.morsch@cern.ch
//
#include <TFile.h>
#include <TParticle.h>
#include <TDatabasePDG.h>
#include <TTree.h>

#include "AliGenReaderEcalJets.h"

ClassImp(AliGenReaderEcalJets)


AliGenReaderEcalJets::AliGenReaderEcalJets():
    fNcurrent(0),
    fNparticle(0),
    fNev(0),
    fNpart(0),
    fNjet(0),     
    fNsjet(0),
    fNpjet(0),
    fTreeNtuple(0) 
{
// Default constructor
    for (Int_t i = 0; i < 200; i++) {
	if (i < 2) {
	    fX[i]     = 0.;    
	    fXtyp[i]  = 0; 
	} else if (i < 10) {
	    fJet[i]   = 0.;  
	    fJeta[i]  = 0.; 
	    fJphi[i]  = 0.; 
	    fJset[i]  = 0.; 
	    fJseta[i] = 0.;
	    fJsphi[i] = 0.;
	    fJpet[i]  = 0.; 
	    fJpeta[i] = 0.;
	    fJpphi[i] = 0.;
	} else {
	    fXpt[i]  = 0.;
	    fXeta[i] = 0.;
	    fXphi[i] = 0.;
	    fXid[i]  = 0; 
	}
    }
}

 AliGenReaderEcalJets::AliGenReaderEcalJets(const AliGenReaderEcalJets &reader)
     :AliGenReader(reader),
      fNcurrent(0),
      fNparticle(0),
      fNev(0),
      fNpart(0),
      fNjet(0),     
      fNsjet(0),
      fNpjet(0),
      fTreeNtuple(0) 
{
    for (Int_t i = 0; i < 200; i++) {
	if (i < 2) {
	    fX[i]     = 0.;    
	    fXtyp[i]  = 0; 
	} else if (i < 10) {
	    fJet[i]   = 0.;  
	    fJeta[i]  = 0.; 
	    fJphi[i]  = 0.; 
	    fJset[i]  = 0.; 
	    fJseta[i] = 0.;
	    fJsphi[i] = 0.;
	    fJpet[i]  = 0.; 
	    fJpeta[i] = 0.;
	    fJpphi[i] = 0.;
	} else {
	    fXpt[i]  = 0.;
	    fXeta[i] = 0.;
	    fXphi[i] = 0.;
	    fXid[i]  = 0; 
	}
    }

    // Copy Constructor
    reader.Copy(*this);
}

void AliGenReaderEcalJets::Init() 
{
//
// reset the existing file environment and open a new root file if
// the pointer to the Fluka tree is null
    
    TFile *pFile=0;
    if (!pFile) {
	pFile = new TFile(fFileName);
	pFile->cd();
	printf("\n I have opened %s file \n", fFileName);
    }
// get the tree address in the Fluka boundary source file
    fTreeNtuple = (TTree*)gDirectory->Get("h1");
    TTree *h1=fTreeNtuple;
    h1->SetMakeClass(1);
//Set branch addresses
    h1->SetBranchAddress("nev",   &fNev);
    h1->SetBranchAddress("x",      fX);
    h1->SetBranchAddress("xtyp",   fXtyp);
    h1->SetBranchAddress("npart", &fNpart);
    h1->SetBranchAddress("xpt",    fXpt);
    h1->SetBranchAddress("xeta",   fXeta);
    h1->SetBranchAddress("xphi",   fXphi);
    h1->SetBranchAddress("xid",    fXid);
    h1->SetBranchAddress("njet",  &fNjet);
    h1->SetBranchAddress("jet",    fJet);
    h1->SetBranchAddress("jeta",   fJeta);
    h1->SetBranchAddress("jphi",   fJphi);
    h1->SetBranchAddress("nsjet", &fNsjet);
    h1->SetBranchAddress("jset",   fJset);
    h1->SetBranchAddress("jseta",  fJseta);
    h1->SetBranchAddress("jsphi",  fJsphi);
    h1->SetBranchAddress("npjet", &fNpjet);
    h1->SetBranchAddress("jpet",   fJpet);
    h1->SetBranchAddress("jpeta",  fJpeta);
    h1->SetBranchAddress("jpphi",  fJpphi);
}

Int_t AliGenReaderEcalJets::NextEvent() 
{
// Read the next event  
    Int_t nTracks=0, nread=0;
    
    TFile* pFile = fTreeNtuple->GetCurrentFile();
    pFile->cd();

    Int_t nentries = (Int_t) fTreeNtuple->GetEntries();
    if (fNcurrent < nentries) {
	Int_t nb = (Int_t)fTreeNtuple->GetEvent(fNcurrent);
	nread += nb;
	fNcurrent++;
	printf("\n Next event contains %d tracks! \n", fNpjet);
	nTracks    = fNpjet;
	fNparticle = 0;
	return nTracks;
    }
    return 0;
}

TParticle* AliGenReaderEcalJets::NextParticle() 
{
// Read the next particle

    Float_t p[4];
    Int_t    ipart  = fXid[fNparticle];
    Float_t  pt     = fXpt[fNparticle];
    Float_t  eta    = fXeta[fNparticle];
    Float_t  phi    = fXphi[fNparticle];
    Float_t  theta  = 2.*TMath::ATan(TMath::Exp(-eta));
    Double_t amass  = TDatabasePDG::Instance()->GetParticle(ipart)->Mass();

    p[0] = pt*TMath::Sin(phi);
    p[1] = pt*TMath::Cos(phi);      
    p[2] = pt/TMath::Cos(theta);
    p[3] = TMath::Sqrt(pt*pt+p[2]*p[2]+amass*amass);
    

    TParticle* particle = 
	new TParticle(ipart, 0, -1, -1, -1, -1, p[0], p[1], p[2], p[3], 
		      0., 0., 0., 0.);
    fNparticle++;
    return particle;
}



AliGenReaderEcalJets& AliGenReaderEcalJets::operator=(const  AliGenReaderEcalJets& rhs)
{
// Assignment operator
    rhs.Copy(*this);
    return (*this);
}

void AliGenReaderEcalJets::Copy(TObject&) const
{
    //
    // Copy 
    //
    Fatal("Copy","Not implemented!\n");
}







