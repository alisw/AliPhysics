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


// Class to read events from external (TNtupla) file
// Events -> neutron removal by EM dissociation of Pb nuclei
// Data from RELDIS code (by I. Pshenichov)

#include <TFile.h>
#include <TParticle.h>
#include <TTree.h>
#include <TVirtualMC.h>
#include <TDatabasePDG.h>
#include <TPDGCode.h>

#include "AliGenReaderEMD.h"

ClassImp(AliGenReaderEMD)


  // -----------------------------------------------------------------------------------
AliGenReaderEMD::AliGenReaderEMD() 
{
// Default constructor
    fStartEvent  = 0;
    fTreeNtuple  = 0;
    fIPSide      = 0;
    fPcToTrack = 0;
}

  // -----------------------------------------------------------------------------------
AliGenReaderEMD::~AliGenReaderEMD()
{
    delete fTreeNtuple;
}

// -----------------------------------------------------------------------------------
AliGenReaderEMD& AliGenReaderEMD::operator=(const  AliGenReaderEMD& rhs)
{
// Assignment operator
    rhs.Copy(*this);
    return *this;
}

// -----------------------------------------------------------------------------------
void AliGenReaderEMD::Copy(TObject&) const
{
    //
    // Copy 
    //
    Fatal("Copy","Not implemented!\n");
}

// -----------------------------------------------------------------------------------
void AliGenReaderEMD::Init() 
{
//
// Reset the existing file environment and open a new root file
    
    TFile *pFile=0;
    if (!pFile) {
	pFile = new TFile(fFileName);
	pFile->cd();
	printf("\n %s file opened\n\n", fFileName);
    }
    fTreeNtuple = (TTree*)gDirectory->Get("h2032");
    fNcurrent = fStartEvent;

    TTree *Ntu=fTreeNtuple;
    //
    // Set branch addresses
    // **** neutrons
    Ntu->SetBranchAddress("Nleft",&fNnLeft);
    Ntu->SetBranchAddress("Eleft",&fEnLeft);
    Ntu->SetBranchAddress("Pxl",  &fPxnLeft);
    Ntu->SetBranchAddress("Pyl",  &fPynLeft);
    Ntu->SetBranchAddress("Pzl",  &fPznLeft);
    Ntu->SetBranchAddress("Nright",&fNnRight);
    Ntu->SetBranchAddress("Eright",&fEnRight);
    Ntu->SetBranchAddress("Pxr",   &fPxnRight);
    Ntu->SetBranchAddress("Pyr",   &fPynRight);
    Ntu->SetBranchAddress("Pzr",   &fPznRight);
    // **** protons
    Ntu->SetBranchAddress("Nleft_p",&fNpLeft);
    Ntu->SetBranchAddress("Etaleft_p",&fEtapLeft);
    Ntu->SetBranchAddress("Pxl_p",  &fPxpLeft);
    Ntu->SetBranchAddress("Pyl_p",  &fPypLeft);
    Ntu->SetBranchAddress("Pzl_p",  &fPzpLeft);
    Ntu->SetBranchAddress("Nright_p",&fNpRight);
    Ntu->SetBranchAddress("Etaright_p",&fEtapRight);
    Ntu->SetBranchAddress("Pxr_p",   &fPxpRight);
    Ntu->SetBranchAddress("Pyr_p",   &fPypRight);
    Ntu->SetBranchAddress("Pzr_p",   &fPzpRight);
}

// -----------------------------------------------------------------------------------
Int_t AliGenReaderEMD::NextEvent() 
{
    // Read the next event  
    Int_t nTracks=0;
    
    TFile* pFile = fTreeNtuple->GetCurrentFile();
    pFile->cd();
    

    Int_t nentries = (Int_t) fTreeNtuple->GetEntries();
    if(fNcurrent < nentries) {
	fTreeNtuple->GetEvent(fNcurrent);
	printf("\n *** Event %d read ***\n\n",fNcurrent);
	fNcurrent++;
	//printf("\n \t \t LEFT-> %d neutrons and %d protons emitted\n", fNnLeft, fNpLeft);
	//printf("\t \t RIGHT-> %d neutrons and %d protons emitted\n\n", fNnRight, fNpRight);
	//
	// **** IPSIde 		=0->RIGHT side, =1->LEFT side  
	// #### fPcToTrack 	=0->neutrons, =1->protons
	if(fIPSide==0){
	  if(fPcToTrack==0){
	    printf("\n \t \t Tracking %d neutrons emitted on RIGHT side\n\n", fNnRight);
	    nTracks    = fNnRight;
	  }
	  else if(fPcToTrack==1){
	    printf("\n \t \t Tracking %d protons emitted on RIGHT side\n\n", fNpRight);
	    nTracks    = fNpRight;
	  }
	}
	else if(fIPSide==1){
	  if(fPcToTrack==0){
	    printf("\n \t \t Tracking %d neutrons emitted on LEFT side\n", fNnLeft);
	    nTracks    = fNnLeft;
	  }
	  else if(fPcToTrack==1){
	    printf("\n \t \t Tracking %d protons emitted on LEFT side\n", fNpLeft);
	    nTracks    = fNpLeft;
	  }
	}
	fNparticle = 0;
	return nTracks;
    }

    return 0;
}

// -----------------------------------------------------------------------------------
TParticle* AliGenReaderEMD::NextParticle() 
{
    // Read the next particle
    Float_t p[4];
    Int_t ipart = kNeutron;
    Double_t amass = TDatabasePDG::Instance()->GetParticle(kNeutron)->Mass();
    p[0] = fPxnRight[fNparticle];
    p[1] = fPynRight[fNparticle];
    p[2] = fPznRight[fNparticle];
    Float_t ptot = TMath::Sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2]);
    p[3] = TMath::Sqrt(ptot*ptot+amass*amass);
    //printf("\t Nucleon momentum: px = %f, py = %f, pz = %f\n", p[0],p[1],p[2]);
    
    if(p[3]<=amass) 
      Warning("Generate","Particle %d  E = %f mass = %f \n",ipart,p[3],amass);

    TParticle* particle = new TParticle(ipart, 0, -1, -1, -1, -1, 
    	p[0], p[1], p[2], p[3], 0., 0., 0., 0.);
    //fNcurrent++;
    fNparticle++;
    return particle;
}
