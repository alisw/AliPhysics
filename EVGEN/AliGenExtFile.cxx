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
Revision 1.10  2000/07/11 18:24:55  fca
Coding convention corrections + few minor bug fixes

Revision 1.9  2000/06/14 15:20:09  morsch
Include clean-up (IH)

Revision 1.8  2000/06/09 20:36:44  morsch
All coding rule violations except RS3 corrected

Revision 1.7  2000/02/16 14:56:27  morsch
Convert geant particle code into pdg code before putting particle on the stack.

Revision 1.6  1999/11/09 07:38:48  fca
Changes for compatibility with version 2.23 of ROOT

Revision 1.5  1999/09/29 09:24:12  fca
Introduction of the Copyright and cvs Log

*/

#include <iostream.h>

#include "AliGenExtFile.h"
#include "AliMC.h"
#include "AliRun.h"

#include <TDirectory.h>
#include <TDatabasePDG.h>
#include <TFile.h>
#include <stdlib.h>

 ClassImp(AliGenExtFile)
     AliGenExtFile::AliGenExtFile()
	 :AliGenerator(-1)
{
//  Constructor
    fName="ExtFile";
    fTitle="Primaries from ext. File";
    fFileName="dtujet93.root";
    fTreeNtuple=0;
    fNcurrent=0;
//
//  Read all particles
    fNpart=-1;
}

AliGenExtFile::AliGenExtFile(Int_t npart)
    :AliGenerator(npart)
{
//  Constructor
    fName="ExtFile";
    fTitle="Primaries from ext. File";
    fFileName="dtujet93.root";
    fTreeNtuple=0;
    fNcurrent=0;
}

AliGenExtFile::AliGenExtFile(const AliGenExtFile & ExtFile)
{
// copy constructor
}
//____________________________________________________________
AliGenExtFile::~AliGenExtFile()
{
// Destructor
    delete fTreeNtuple;
}

//____________________________________________________________
void AliGenExtFile::NtupleInit() 
{
//
// reset the existing file environment and open a new root file if
// the pointer to the Fluka tree is null
    
    TFile *pFile=0;
    if (fTreeNtuple==0) {
        if (!pFile) {
	    pFile = new TFile(fFileName);
	    pFile->cd();
	    cout<<"I have opened "<<fFileName<<" file "<<endl;
        }
// get the tree address in the Fluka boundary source file
	fTreeNtuple = (TTree*)gDirectory->Get("h888");
    } else {
        pFile = fTreeNtuple->GetCurrentFile();
        pFile->cd();
    }

    TTree *h2=fTreeNtuple;
//Set branch addresses
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


//____________________________________________________________
void AliGenExtFile::Generate()
{
// Generate particles

  Float_t polar[3]= {0,0,0};
  //
  Float_t origin[3]={0,0,0};
  Float_t p[3];
  Float_t random[6];
  Float_t prwn;
  Int_t i, j, nt, nTracks=0;
  //
  NtupleInit();
  TTree *h2=fTreeNtuple;
  Int_t nentries = (Int_t) h2->GetEntries();
  // loop over number of particles
  Int_t nb = (Int_t)h2->GetEvent(fNcurrent);
  Int_t i5=fIhead[4];
  Int_t i6=fIhead[5];

  for (j=0;j<3;j++) origin[j]=fOrigin[j];
  if(fVertexSmear==kPerTrack) {
    gMC->Rndm(random,6);
    for (j=0;j<3;j++) {
	origin[j]+=fOsigma[j]*TMath::Cos(2*random[2*j]*TMath::Pi())*
	    TMath::Sqrt(-2*TMath::Log(random[2*j+1]));
    }
  }

  if (fNcurrent >= nentries) {
      printf("\n No more entries !!! !\n");
      return;
  }
  
	  
  if (i5==0) {
      printf("\n This should never happen !\n");
  } else {
      printf("\n Next event contains %d tracks! \n", i6);
      nTracks=i6;
  }
  for (i=0; i<nTracks; i++) {
      fIdpart=gMC->PDGFromId(fIdpart);
      Double_t amass = TDatabasePDG::Instance()->GetParticle(fIdpart)->Mass();
      if(fE<=amass) {
	Warning("Generate","Particle %d no %d E = %f mass = %f\n",
		fIdpart,i,fE,amass);
	prwn=0;
      } else {
	prwn=sqrt((fE+amass)*(fE-amass));
      }

      fTheta *= TMath::Pi()/180.;
      fPhi    = (fPhi-180)*TMath::Pi()/180.;      
      if(fTheta<fThetaMin || fTheta>fThetaMax ||
	 fPhi<fPhiMin || fPhi>fPhiMax         ||
	 prwn<fPMin || prwn>fPMax)          
      {
	  ;
      } else {
	  p[0]=prwn*TMath::Sin(fTheta)*TMath::Cos(fPhi);
	  p[1]=prwn*TMath::Sin(fTheta)*TMath::Sin(fPhi);      
	  p[2]=prwn*TMath::Cos(fTheta);
	  
	  if(fVertexSmear==kPerTrack) {
	      gMC->Rndm(random,6);
	      for (j=0;j<3;j++) {
		  origin[j]=fOrigin[j]
		      +fOsigma[j]*TMath::Cos(2*random[2*j]*TMath::Pi())*
		      TMath::Sqrt(-2*TMath::Log(random[2*j+1]));
	      }
	  }
	  gAlice->SetTrack(fTrackIt,-1,fIdpart,p,origin,polar,0,"Primary",nt);
      }
      fNcurrent++;
      nb = (Int_t)h2->GetEvent(fNcurrent); 
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


AliGenExtFile& AliGenExtFile::operator=(const  AliGenExtFile& rhs)
{
// Assignment operator
    return *this;
}








