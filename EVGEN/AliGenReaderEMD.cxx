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
#include "AliStack.h"


ClassImp(AliGenReaderEMD)

AliGenReaderEMD::AliGenReaderEMD():
    fStartEvent(0),
    fNcurrent(0),  
    fNparticle(0), 
    fTreeNtuple(0),
    fIPSide(0),
    fPcToTrack(0),
    fNnAside(0),
    fEnAside(0),
    fNnCside(0),
    fEnCside(0),
    fNpAside(0),
    fEtapAside(0),
    fNpCside(0),
    fEtapCside(0)
{
// Default constructor
}

AliGenReaderEMD::AliGenReaderEMD(const AliGenReaderEMD &reader):
    AliGenReader(reader),
    fStartEvent(0),
    fNcurrent(0),  
    fNparticle(0), 
    fTreeNtuple(0),
    fIPSide(0),
    fPcToTrack(0),
    fNnAside(0),
    fEnAside(0),
    fNnCside(0),
    fEnCside(0),
    fNpAside(0),
    fEtapAside(0),
    fNpCside(0),
    fEtapCside(0)
{
    // Copy Constructor
    reader.Copy(*this);
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
    Ntu->SetBranchAddress("Nleft",&fNnAside);
    Ntu->SetBranchAddress("Eleft",&fEnAside);
    Ntu->SetBranchAddress("Pxl",  fPxnAside);
    Ntu->SetBranchAddress("Pyl",  fPynAside);
    Ntu->SetBranchAddress("Pzl",  fPznAside);
    Ntu->SetBranchAddress("Nright",&fNnCside);
    Ntu->SetBranchAddress("Eright",&fEnCside);
    Ntu->SetBranchAddress("Pxr",   fPxnCside);
    Ntu->SetBranchAddress("Pyr",   fPynCside);
    Ntu->SetBranchAddress("Pzr",   fPznCside);
    // **** protons
    Ntu->SetBranchAddress("Nleft_p",&fNpAside);
    Ntu->SetBranchAddress("Etaleft_p",&fEtapAside);
    Ntu->SetBranchAddress("Pxl_p",  fPxpAside);
    Ntu->SetBranchAddress("Pyl_p",  fPypAside);
    Ntu->SetBranchAddress("Pzl_p",  fPzpAside);
    Ntu->SetBranchAddress("Nright_p",&fNpCside);
    Ntu->SetBranchAddress("Etaright_p",&fEtapCside);
    Ntu->SetBranchAddress("Pxr_p",   fPxpCside);
    Ntu->SetBranchAddress("Pyr_p",   fPypCside);
    Ntu->SetBranchAddress("Pzr_p",   fPzpCside);
}

// -----------------------------------------------------------------------------------
Int_t AliGenReaderEMD::NextEvent() 
{
    // Read the next event  
    Int_t nTracks=0;
    fNparticle = 0;
    
    TFile* pFile = fTreeNtuple->GetCurrentFile();
    pFile->cd();
    

    Int_t nentries = (Int_t) fTreeNtuple->GetEntries();
    if(fNcurrent < nentries) {
	fTreeNtuple->GetEvent(fNcurrent);
	printf("\n *** Event %d read ***\n",fNcurrent);
	//printf("\t A side-> %d neutrons and %d protons emitted\n", fNnAside, fNpAside);
	//printf("\t C side-> %d neutrons and %d protons emitted\n\n", fNnCside, fNpCside);
	//
	// **** IPSIde = 0->both sides, 1->only Cside, 2->only A side  
	// #### fPcToTrack = 0->neutrons, 1->protons
	if(fPcToTrack==0){ // neutrons
	  if(fIPSide==0){
	    printf("\t Tracking %d+%d neutrons emitted on C+A side\n", fNnCside,fNnAside);
	    nTracks = fNnCside+fNnAside;
	  }
	  else if(fIPSide==1){ 
	    printf("\t Tracking %d neutrons emitted on C side\n", fNnCside);
	    nTracks    = fNnCside;
	  }
	  else if(fIPSide==2){
	    printf("\t Tracking %d neutrons emitted on A side\n", fNnAside);
	    nTracks    = fNnAside;
	  }
	}
	else if(fPcToTrack==1){ //protons
	  if(fIPSide==0){
	    printf("\t Tracking %d+%d protons emitted on C+A side\n", fNpCside,fNpAside);
	    nTracks    = fNpCside+fNpAside;
	  }
	  else if(fIPSide==1){
	    printf("\t Tracking %d protons emitted on C side\n", fNpCside);
	    nTracks    = fNpCside;
	  }
	  else if(fIPSide==2){
	    printf("\t Tracking %d protons emitted on A side\n", fNpAside);
	    nTracks    = fNpAside;
	  }
	}
	fNcurrent++;
	return nTracks;
    }

    return 0;
}

// -----------------------------------------------------------------------------------
TParticle* AliGenReaderEMD::NextParticle() 
{
    // Read the next particle
    Float_t p[4]={0.,0.,0.,0.};
    
    Int_t ipart=0;
    if(fPcToTrack==0) ipart = kNeutron;
    else  if(fPcToTrack==1) ipart = kProton;
    Double_t amass = TDatabasePDG::Instance()->GetParticle(ipart)->Mass();

    if(fPcToTrack==0){
      if(fIPSide==0){
        if(fNparticle<fNnCside){
          p[0] = fPxnCside[fNparticle];
          p[1] = fPynCside[fNparticle];
          p[2] = fPznCside[fNparticle];
	}
	else{
          p[0] = fPxnAside[fNparticle-fNnCside];
          p[1] = fPynAside[fNparticle-fNnCside];
          p[2] = fPznAside[fNparticle-fNnCside];
	}
      }
      else if(fIPSide==1){
        p[0] = fPxnCside[fNparticle];
        p[1] = fPynCside[fNparticle];
        p[2] = fPznCside[fNparticle];
      }
      else if(fIPSide==2){
        p[0] = fPxnAside[fNparticle-fNpCside];
        p[1] = fPynAside[fNparticle-fNpCside];
        p[2] = fPznAside[fNparticle-fNpCside];
      }
    }
    else if(fPcToTrack==1){
      if(fIPSide==0){
        if(fNparticle<fNnCside){
          p[0] = fPxpCside[fNparticle];
          p[1] = fPypCside[fNparticle];
          p[2] = fPzpCside[fNparticle];
	}
	else{
          p[0] = fPxpAside[fNparticle];
          p[1] = fPypAside[fNparticle];
          p[2] = fPzpAside[fNparticle];
	}
      }
      else if(fIPSide==1){
        p[0] = fPxpCside[fNparticle];
        p[1] = fPypCside[fNparticle];
        p[2] = fPzpCside[fNparticle];
      }
      else if(fIPSide==2){
        p[0] = fPxpAside[fNparticle];
        p[1] = fPypAside[fNparticle];
        p[2] = fPzpAside[fNparticle];
      }
    } 
    Float_t ptot = TMath::Sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2]);
    p[3] = TMath::Sqrt(ptot*ptot+amass*amass);
    //printf(" * pc %d: momentum (%f, %f, %f) \n", fNparticle, p[0],p[1],p[2]);
    
    if(p[3]<=amass) 
      Warning("Generate","Particle %d  E = %f mass = %f \n",ipart,p[3],amass);

    TParticle* particle = new TParticle(ipart, 0, -1, -1, -1, -1, 
    	p[0], p[1], p[2], p[3], 0., 0., 0., 0.);
    particle->SetBit(kTransportBit);
    fNparticle++;
    return particle;
}
