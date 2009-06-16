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
// It reads Hijing events from a ntuple like event structure.
// The event format is defined in Init()
// NextEvent() is used to loop over events and NextParticle() to loop over particles.  
// Author: andreas.morsch@cern.ch
//
#include <TFile.h>
#include <TParticle.h>
#include <TDatabasePDG.h>
#include <TTree.h>

#include "AliGenReaderEcalHijing.h"

ClassImp(AliGenReaderEcalHijing)

AliGenReaderEcalHijing::AliGenReaderEcalHijing():
    fNcurrent(0),
    fNparticle(0),
    fTreeNtuple(0),
    fNjatt(0),
    fNahij(0),
    fNphij(0)
{
// Default constructor
}

AliGenReaderEcalHijing::AliGenReaderEcalHijing(const AliGenReaderEcalHijing &reader):
    AliGenReader(reader),
    fNcurrent(0),
    fNparticle(0),
    fTreeNtuple(0),
    fNjatt(0),
    fNahij(0),
    fNphij(0)
{
    // Copy constructor
    reader.Copy(*this);
}

void AliGenReaderEcalHijing::Init() 
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
    fTreeNtuple = (TTree*)gDirectory->Get("h2");
    TTree *h2=fTreeNtuple;
    h2->SetMakeClass(1);
//Set branch addresses
    h2->SetBranchAddress("njatt", &fNjatt);
    h2->SetBranchAddress("nahij", &fNahij);
    h2->SetBranchAddress("nphij", &fNphij);
    h2->SetBranchAddress("khij",   fKhij) ;
    h2->SetBranchAddress("pxhij",  fPxhij);
    h2->SetBranchAddress("pyhij",  fPyhij);
    h2->SetBranchAddress("pzhij",  fPzhij);
    h2->SetBranchAddress("ehij",   fEhij) ;
}

Int_t AliGenReaderEcalHijing::NextEvent() 
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
	printf("\n Next event contains %d tracks! \n", fNphij);
	nTracks    = fNphij;
	fNparticle = 0;
	return nTracks;
    }
    return 0;
}

TParticle* AliGenReaderEcalHijing::NextParticle() 
{
// Read the next particle

    Float_t p[4];
    Int_t ipart = fKhij[fNparticle];
    p[0] = fPxhij[fNparticle];
    p[1] = fPyhij[fNparticle];      
    p[2] = fPzhij[fNparticle];
    p[3] = fEhij[fNparticle];
    
    Double_t amass = TDatabasePDG::Instance()->GetParticle(ipart)->Mass();

    if(p[3] <= amass) {
	Warning("Generate","Particle %d  E = %f mass = %f \n",
		ipart, p[3], amass);
    } 
    TParticle* particle = 
	new TParticle(ipart, 0, -1, -1, -1, -1, p[0], p[1], p[2], p[3], 
		      0., 0., 0., 0.);
    fNparticle++;
    return particle;
}



AliGenReaderEcalHijing& AliGenReaderEcalHijing::operator=(const  AliGenReaderEcalHijing& rhs)
{
// Assignment operator
    rhs.Copy(*this);
    return (*this);
}

void AliGenReaderEcalHijing::Copy(TObject&) const
{
    //
    // Copy 
    //
    Fatal("Copy","Not implemented!\n");
}






