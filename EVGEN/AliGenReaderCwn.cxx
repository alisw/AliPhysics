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

// Read the old ALICE event format based on CW-ntuples
// http://consult.cern.ch/alice/Internal_Notes/1995/32/abstract
// .cwn file have to be converted to .root using h2root
// Use SetFileName(file) to read from "file" 
// Author: andreas.morsch@cern.ch

#include <TFile.h>
#include <TParticle.h>
#include <TDatabasePDG.h>
#include <TTree.h>
#include <TVirtualMC.h>

#include "AliGenReaderCwn.h"

ClassImp(AliGenReaderCwn)

AliGenReaderCwn::AliGenReaderCwn():
    fNcurrent(0),
    fNparticle(0),
    fNparticleMax(0),
    fTreeNtuple(0),
    fNihead(0),
    fNrhead(0),
    fIdpart(0),
    fTheta(0.),
    fPhi(0.),
    fP(0.),
    fE(0.)
{
// Default constructor
    Int_t i;
    for (i = 0; i <  6; i++) fRhead[i] = 0.;
    for (i = 0; i < 12; i++) fIhead[i] = 0;
}


AliGenReaderCwn::AliGenReaderCwn(const AliGenReaderCwn &reader):
    AliGenReader(reader),
    fNcurrent(0),
    fNparticle(0),
    fNparticleMax(0),
    fTreeNtuple(0),
    fNihead(0),
    fNrhead(0),
    fIdpart(0),
    fTheta(0.),
    fPhi(0.),
    fP(0.),
    fE(0.)
{
    // Copy constructor
    Int_t i;
    for (i = 0; i <  6; i++) fRhead[i] = 0.;
    for (i = 0; i < 12; i++) fIhead[i] = 0;
    reader.Copy(*this);
}


AliGenReaderCwn::~AliGenReaderCwn()
{
    delete fTreeNtuple;
}

void AliGenReaderCwn::Init() 
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
    fTreeNtuple = (TTree*)gDirectory->Get("h888");

    TTree *h2=fTreeNtuple;
//Set branch addresses
    h2->SetBranchAddress("Nihead",&fNihead);
    h2->SetBranchAddress("Ihead",fIhead);
    h2->SetBranchAddress("Nrhead",&fNrhead);
    h2->SetBranchAddress("Rhead",fRhead);
    h2->SetBranchAddress("Idpart",&fIdpart);
    h2->SetBranchAddress("Theta",&fTheta);
    h2->SetBranchAddress("Phi",&fPhi);
    h2->SetBranchAddress("P",&fP);
    h2->SetBranchAddress("E",&fE);
}

Int_t AliGenReaderCwn::NextEvent() 
{
// Read the next event  
    Int_t nTracks;
    fNparticle = 0;
    TFile* pFile = fTreeNtuple->GetCurrentFile();
    pFile->cd();

    Int_t nentries = (Int_t) fTreeNtuple->GetEntries();
    if (fNcurrent < nentries) {
	fNcurrent++;
	
	Int_t i5=fIhead[4];
	Int_t i6=fIhead[5];
	if (i5==0) {
	    printf("\n This should never happen !\n");
	    nTracks = 0;
	} else {
	    printf("\n Next event contains %d tracks! \n", i6);
	    nTracks = i6;
	}    
	fNparticleMax = nTracks;
	return nTracks;
    }

    return 0;
}

TParticle* AliGenReaderCwn::NextParticle() 
{
// Read next particle
//  
    Float_t prwn;
    Float_t p[4];
// Read the next particle
    if (fCode == kGEANT3) fIdpart=gMC->PDGFromId(fIdpart);
    Double_t amass = TDatabasePDG::Instance()->GetParticle(fIdpart)->Mass();
    if(fE<=amass) {
	Warning("Generate","Particle %d  E = %f mass = %f %f %f \n",
		fIdpart,fE,amass, fPhi, fTheta);
	prwn=0;
    } else {
	prwn=sqrt((fE+amass)*(fE-amass));
    }

    fTheta *= TMath::Pi()/180.;
    fPhi    = (fPhi-180)*TMath::Pi()/180.;      
    p[0] = prwn*TMath::Sin(fTheta)*TMath::Cos(fPhi);
    p[1] = prwn*TMath::Sin(fTheta)*TMath::Sin(fPhi);      
    p[2] = prwn*TMath::Cos(fTheta);
    p[3] = fE;
    TParticle* particle = new TParticle(fIdpart, 0, -1, -1, -1, -1, p[0], p[1], p[2], p[3], 0., 0., 0., 0.);
    fNcurrent++;
    fNparticle++;
    return particle;
}



AliGenReaderCwn& AliGenReaderCwn::operator=(const  AliGenReaderCwn& rhs)
{
// Assignment operator
    rhs.Copy(*this);
    return *this;
}

void AliGenReaderCwn::Copy(TObject&) const
{
    //
    // Copy 
    //
    Fatal("Copy","Not implemented!\n");
}




