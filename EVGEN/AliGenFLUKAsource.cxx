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
Revision 1.8  2000/02/14 14:49:38  morsch
Correct particle type for gamma and neutrons
More consistent calculation of momentum from kin. energy and mass

Revision 1.7  1999/11/03 17:43:20  fca
New version from G.Martinez & A.Morsch

Revision 1.6  1999/09/29 09:24:12  fca
Introduction of the Copyright and cvs Log

*/

#include "AliGenFLUKAsource.h"
#include "AliGenMUONlib.h"
#include "AliMC.h"
#include "AliRun.h"
#include "AliPDG.h"
#include <TDirectory.h>
#include <TFile.h>
#include <TTree.h>
#include <stdlib.h>
 ClassImp(AliGenFLUKAsource)
     AliGenFLUKAsource::AliGenFLUKAsource()
	 :AliGenerator(-1)
{
    //
    fName="FLUKA";
    fTitle="FLUKA Boundary Source";
    // Read in all particle types by default
    fIkine=6;
    // Set maximum admitted age of particles to 1.e-05 by default 
    fAgeMax=1.e-05;
    // Do not add weight
    fAddWeight=1.;
    // Shift the z-coordinate of the impact point by 4.5 cm only if it reads 
    // from  specific boundary source to the chamber (fZshift=4.5;),else there 
    // is no need to shift as it reads boundary source for the whole volume of 
    // the Muon Arm; the default value corresponds to boundary source for the
    // whole volume of the MUON Arm 
    fZshift=0;
    // Set the default file 
    fFileName="flukasource.root";

    fTreeFluka=0;
    fTreeChain = new TChain("h1");
//
//  Read all particles
    fNpart=-1;
    
}

AliGenFLUKAsource::AliGenFLUKAsource(Int_t npart)
    :AliGenerator(npart)
{
    //
    fName="FLUKA";
    fTitle="FLUKA Boundary Source";
    // Read in all particle types by default
    fIkine=6;
    // Set maximum admitted age of particles to 1.e-05 by default 
    fAgeMax=1.e-05;
    // Do not add weight
    fAddWeight=1.;
    // Shift the z-coordinate of the impact point by 4.5 cm only if it reads 
    // from  specific boundary source to the chamber (fZshift=4.5;),else there 
    // is no need to shift as it reads boundary source for the whole volume of 
    // the Muon Arm; the default value corresponds to boundary source for the
    // whole volume of the MUON Arm 
    fZshift=0;
    // Set the default file 
    fFileName="flukasource.root";

    fTreeFluka=0;
    fTreeChain = new TChain("h1"); 
}

//____________________________________________________________
AliGenFLUKAsource::~AliGenFLUKAsource()
{
 if (fTreeFluka) delete fTreeFluka;
 if (fTreeChain) delete fTreeChain;
 // if (fFileName)  delete fFileName;
}

void AliGenFLUKAsource::AddFile(const Text_t *filename)
{
    fTreeChain->Add(filename);
    
}


//____________________________________________________________
void AliGenFLUKAsource::FlukaInit() 
{
    TChain *h2=fTreeChain;
//Set branch addresses
    h2->SetBranchAddress("Ip",&Ip);
    h2->SetBranchAddress("Ipp",&Ipp);
    h2->SetBranchAddress("Xi",&Xi);
    h2->SetBranchAddress("Yi",&Yi);
    h2->SetBranchAddress("Zi",&Zi);
    h2->SetBranchAddress("Px",&Px);
    h2->SetBranchAddress("Py",&Py);
    h2->SetBranchAddress("Pz",&Pz);
    h2->SetBranchAddress("Ekin",&Ekin);
    h2->SetBranchAddress("Zv",&Zv);
    h2->SetBranchAddress("Rv",&Rv);
    h2->SetBranchAddress("Itra",&Itra);
    h2->SetBranchAddress("Igas",&Igas);
    h2->SetBranchAddress("Wgt",&Wgt);
    h2->SetBranchAddress("Etag",&Etag);
    h2->SetBranchAddress("Ptg",&Ptg);
    h2->SetBranchAddress("Age",&Age);
}

//____________________________________________________________
void AliGenFLUKAsource::Generate()
{ 
    AliMC* gMC = AliMC::GetMC();

    const Int_t ifluge[28]={kProton, kProtonBar, kElectron, kPositron,
			  kNuE, kNuEBar, kGamma, kNeutron, kNeutronBar,
			  kMuonPlus, kMuonMinus, kK0Long , kPiPlus, kPiMinus,
			  kKPlus, kKMinus, kLambda0, kLambda0Bar, kK0Short,
			  kSigmaMinus, kSigmaPlus, kSigma0, kPi0, kK0, kK0Bar,
			  0,kNuMu,kNuMuBar};
    Float_t polar[3]= {0,0,0};
  //
    Float_t origin[3];
    Float_t p[3];
    Float_t prwn;
    Float_t wgt, fwgt;
    Float_t phi;
    char name[100];
    Float_t amass, charge, tlife;
    Int_t itrtyp;
    Int_t iwgt;
    Int_t i, j, part, nt;
    static Int_t irwn=0;
    //
    Float_t random[2];
  //
    FlukaInit();
    TChain *h2=fTreeChain;
    Int_t nentries = (Int_t) h2->GetEntries();
    if (fNpart == -1) fNpart=Int_t(nentries*fFrac);
  
  // loop over number of particles
    Int_t nb=0;
    Int_t ev=gMC->CurrentEvent();
    for (i=0; i<fNpart;i++) {
	Int_t entry=fNpart*(ev-1)+i; 
	nb = (Int_t)h2->GetEvent(entry); 
	if (irwn > nentries) {
	    printf("No more entries in the FLUKA boundary source file\n");
	    TFile *File=0;
	    // Get AliRun object or create it 
	    if (!gAlice) {
		gAlice = (AliRun*)File->Get("gAlice");
		if (gAlice) printf("AliRun object found on file\n");
		if (!gAlice) gAlice = new AliRun("gAlice","Alice test program");
	    }
	    TTree *fAli=gAlice->TreeK();
	    if (fAli) File =fAli->GetCurrentFile();
	    File->cd();
	    printf("Generate - I'm out \n");
	    return;
	}   
	if (Ip > 28 || Ip < 0) {
	    irwn++;
	    continue;
	}
	
	if ((Ip != fIkine && fIkine != 6 && fIkine != 9 && fIkine != 10) || Age > fAgeMax){
	    irwn++;
	    continue;
	} else if (fIkine == 9) {
	    if (Ip == 7 || Ip == 8 || Age > fAgeMax) { 
		irwn++;
		continue;
	    }
	} else if (fIkine == 10) {
	    if (Ip == 8 || Age > fAgeMax) { 
		irwn++;
		continue;
	    }
	}
    

	irwn++;
//
// PDG code from FLUKA particle type (Ip)
	part=ifluge[int(Ip)-1];	
//
// Calculate momentum from kinetic energy and mass of the particle
	gMC->Gfpart(part, name, itrtyp,  
		    amass, charge, tlife); 
	prwn=Ekin*sqrt(1. + 2.*amass/Ekin);


	origin[0]=Yi;
	origin[1]=Xi;
	origin[2]=Zi;
	
	p[0]=Py*prwn;
	p[1]=Px*prwn;
	p[2]=Pz*prwn;

	//handle particle weight correctly
	wgt = (part == 13) ? Wgt*fAddWeight : Wgt;
	iwgt=Int_t(wgt);
	fwgt=wgt-Float_t(iwgt);
	gMC->Rndm(random,2);
	if (random[0] < fwgt) iwgt++;
	if (part==1 && iwgt>100) iwgt=100;
	Int_t nstack=0;
	for (j=0; j<iwgt; j++) {
	    gAlice->SetTrack(fTrackIt,-1,part,p,origin,polar,Age,"Primary",nt);
	    gMC->Rndm(random,2);
	    phi=2*random[1]*TMath::Pi();
	    Float_t pn1=p[0]*TMath::Sin(phi) - p[1]*TMath::Cos(phi);
	    Float_t pn2=p[0]*TMath::Cos(phi) + p[1]*TMath::Sin(phi);
	    p[0]=pn1;
	    p[1]=pn2;
	    Float_t on1=origin[0]*TMath::Sin(phi)-origin[1]*TMath::Cos(phi);
	    Float_t on2=origin[0]*TMath::Cos(phi)+origin[1]*TMath::Sin(phi);
	    origin[0]=on1;
	    origin[1]=on2;
	    nstack++;
	}
	if (nstack == 0) continue;
    }
    
    TFile *File=0;
// Get AliRun object or create it 
    if (!gAlice) {
	gAlice = (AliRun*)File->Get("gAlice");
	if (gAlice) printf("AliRun object found on file\n");
	if (!gAlice) gAlice = new AliRun("gAlice","Alice test program");
    }
    TTree *fAli=gAlice->TreeK();
    if (fAli) File =fAli->GetCurrentFile();
    File->cd();
}









