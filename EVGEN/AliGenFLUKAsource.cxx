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

// Read background particles from a FLUKA boundary source file
// This is a very special generator that works for background studies for the muon-spectrometer.
// The input files come from FLUKA simulations.
// Files can be chained. 
// Author: andreas.morsch@cern.ch

#include "AliGenFLUKAsource.h"
#include <stdlib.h>

#include <TDatabasePDG.h>
#include <TPDGCode.h>
#include <RVersion.h>
#include <TChain.h>
#include <TFile.h>
#include <TTree.h>
#include <TVirtualMC.h>

#include "AliRun.h"

ClassImp(AliGenFLUKAsource)

AliGenFLUKAsource::AliGenFLUKAsource()
    :AliGenerator(-1), 
     fIkine(6),
     fAgeMax(1.e-5), 
     fAddWeight(1.),
     fZshift(0.),
     fFrac(0.),
     fSourceId(-1),
     fFileName(0),
     fTreeChain(0),
     fTreeFluka(0),
     fIp(0.),
     fIpp(0.),
     fXi(0.),
     fYi(0.),
     fZi(0.),
     fPx(0.),
     fPy(0.),
     fPz(0.),
     fEkin(0.),
     fZv(0.),
     fRv(0.),
     fItra(0.),
     fIgas(0.),
     fWgt(0.),
     fEtag(0.),
     fPtg(0.),
     fAge(0.)
{
    // Constructor
    fName="FLUKA";
    fTitle="FLUKA Boundary Source";
    // Read in all particle types by default
    // Set maximum admitted age of particles to 1.e-05 by default 
    // Do not add weight
    // Shift the z-coordinate of the impact point by 4.5 cm only if it reads 
    // from  specific boundary source to the chamber (fZshift=4.5;),else there 
    // is no need to shift as it reads boundary source for the whole volume of 
    // the Muon Arm; the default value corresponds to boundary source for the
    // whole volume of the MUON Arm 
    // Set the default file 
    fTreeChain = new TChain("h1");
//
//  Read all particles
    fNpart=-1;
}

AliGenFLUKAsource::AliGenFLUKAsource(Int_t npart)
    :AliGenerator(npart), 
     fIkine(6),
     fAgeMax(1.e-5), 
     fAddWeight(1.),
     fZshift(0.),
     fFrac(0.),
     fSourceId(-1),
     fFileName(""),
     fTreeChain(new TChain("h1")),
     fTreeFluka(0),
     fIp(0.),
     fIpp(0.),
     fXi(0.),
     fYi(0.),
     fZi(0.),
     fPx(0.),
     fPy(0.),
     fPz(0.),
     fEkin(0.),
     fZv(0.),
     fRv(0.),
     fItra(0.),
     fIgas(0.),
     fWgt(0.),
     fEtag(0.),
     fPtg(0.),
     fAge(0.)
{
    // Constructor
    fName  = "FLUKA";
    fTitle = "FLUKA Boundary Source";
}

//____________________________________________________________
AliGenFLUKAsource::~AliGenFLUKAsource()
{
// Destructor
 if (fTreeFluka) delete fTreeFluka;
 if (fTreeChain) delete fTreeChain;
}

void AliGenFLUKAsource::AddFile(const Text_t *filename)
{
// Add a file to the chain
    fTreeChain->Add(filename);
    
}


//____________________________________________________________
void AliGenFLUKAsource::FlukaInit() 
{
// Set branch addresses of data entries
    TChain *h2=fTreeChain;
//Set branch addresses
    h2->SetBranchAddress("Ip",&fIp);
    h2->SetBranchAddress("Ipp",&fIpp);
    h2->SetBranchAddress("Xi",&fXi);
    h2->SetBranchAddress("Yi",&fYi);
    h2->SetBranchAddress("Zi",&fZi);
    h2->SetBranchAddress("Px",&fPx);
    h2->SetBranchAddress("Py",&fPy);
    h2->SetBranchAddress("Pz",&fPz);
    h2->SetBranchAddress("Ekin",&fEkin);
    h2->SetBranchAddress("Zv",&fZv);
    h2->SetBranchAddress("Rv",&fRv);
    h2->SetBranchAddress("Itra",&fItra);
    h2->SetBranchAddress("Igas",&fIgas);
    h2->SetBranchAddress("Wgt",&fWgt);
    h2->SetBranchAddress("Etag",&fEtag);
    h2->SetBranchAddress("Ptg",&fPtg);
    h2->SetBranchAddress("Age",&fAge);
}

//____________________________________________________________
void AliGenFLUKAsource::Generate()
{
// Generate one event 

    const Int_t kIfluge[28]={kProton, kProtonBar, kElectron, kPositron,
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
    Float_t amass;
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
	    TFile *pFile=0;
	    // Get AliRun object or create it 
	    if (!gAlice) {
		gAlice = (AliRun*)pFile->Get("gAlice");
		if (gAlice) printf("AliRun object found on file\n");
		if (!gAlice) gAlice = new AliRun("gAlice","Alice test program");
	    }
	    TTree *fAli=gAlice->TreeK();
	    if (fAli) pFile =fAli->GetCurrentFile();
	    pFile->cd();
	    printf("Generate - I'm out \n");
	    return;
	}   
	
	Int_t ifip = Int_t(fIp);
	

	if (fSourceId != -1 && fIgas !=fSourceId) {
	    irwn++;
	    continue;
	}
	
	if (ifip > 28 || ifip < 0) {
	    irwn++;
	    continue;
	}
	
	if ((ifip != fIkine && fIkine != kAll && fIkine != kCharged 
	     && fIkine != 10) || fAge > fAgeMax){
	    irwn++;
	    continue;
	} else if (fIkine == kCharged) {
	    if (ifip == 7 || ifip == 8 || fAge > fAgeMax) { 
		irwn++;
		continue;
	    }
	} else if (fIkine == kNoNeutron) {
	    if (ifip == 8 || fAge > fAgeMax) { 
		irwn++;
		continue;
	    }
	}
    

	irwn++;
//
// PDG code from FLUKA particle type (ifip)
	part=kIfluge[int(ifip)-1];	
//
// Calculate momentum from kinetic energy and mass of the particle
#if ROOT_VERSION_CODE > 197895
        amass = gMC->ParticleMass(part);
#else
	amass = (TDatabasePDG::Instance())->GetParticle(part)->Mass();
#endif
	prwn=fEkin*sqrt(1. + 2.*amass/fEkin);


	origin[0]=fYi;
	origin[1]=fXi;
	origin[2]=fZi;
	
	p[0]=fPy*prwn;
	p[1]=fPx*prwn;
	p[2]=fPz*prwn;

	//handle particle weight correctly
	wgt = (part == 13) ? fWgt*fAddWeight : fWgt;
	iwgt=Int_t(wgt);
	fwgt=wgt-Float_t(iwgt);
	Rndm(random,2);
	if (random[0] < fwgt) iwgt++;
	if (part==1 && iwgt>100) iwgt=100;
	Int_t nstack=0;
	for (j=0; j<iwgt; j++) {
	    PushTrack(fTrackIt,-1,part,p,origin,polar,fAge,kPPrimary,nt);
	    Rndm(random,2);
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
    
    TFile *pFile=0;
// Get AliRun object or create it 
    if (!gAlice) {
	gAlice = (AliRun*)pFile->Get("gAlice");
	if (gAlice) printf("AliRun object found on file\n");
	if (!gAlice) gAlice = new AliRun("gAlice","Alice test program");
    }
    TTree *fAli=gAlice->TreeK();
    if (fAli) pFile =fAli->GetCurrentFile();
    pFile->cd();
}






