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
Revision 1.17  2001/11/09 09:12:58  morsch
Generalization by using AliGenReader object to read particles from file.

Revision 1.16  2001/07/27 17:09:35  morsch
Use local SetTrack, KeepTrack and SetHighWaterMark methods
to delegate either to local stack or to stack owned by AliRun.
(Piotr Skowronski, A.M.)

Revision 1.15  2001/01/23 13:29:37  morsch
Add method SetParticleCode and enum type Code_t to handle both PDG (new ntuples)
and GEANT3 codes (old ntuples) in input file.

Revision 1.14  2000/12/21 16:24:06  morsch
Coding convention clean-up

Revision 1.13  2000/11/30 07:12:50  alibrary
Introducing new Rndm and QA classes

Revision 1.12  2000/10/27 13:54:45  morsch
Remove explicite reference to input file from constuctor.

Revision 1.11  2000/10/02 21:28:06  fca
Removal of useless dependecies via forward declarations

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


// Event generator that using an instance of type AliGenReader
// reads particles from a file and applies cuts. 

#include <iostream.h>

#include "AliGenExtFile.h"
#include "AliRun.h"

#include <TParticle.h>
#include <TFile.h>
#include <TTree.h>


 ClassImp(AliGenExtFile)
     AliGenExtFile::AliGenExtFile()
	 :AliGenerator(-1)
{
//  Constructor
    fName="ExtFile";
    fTitle="Primaries from ext. File";
//
//  Read all particles
    fNpart  =- 1;
    fReader =  0;
}

AliGenExtFile::AliGenExtFile(Int_t npart)
    :AliGenerator(npart)
{
//  Constructor
    fName   = "ExtFile";
    fTitle  = "Primaries from ext. File";
    fReader = 0;
}

AliGenExtFile::AliGenExtFile(const AliGenExtFile & ExtFile)
{
// copy constructor
}
//____________________________________________________________
AliGenExtFile::~AliGenExtFile()
{
// Destructor
    delete fReader;
}

//___________________________________________________________
void AliGenExtFile::Init()
{
// Initialize
    if (fReader) fReader->Init();
}

    
void AliGenExtFile::Generate()
{
// Generate particles

  Float_t polar[3]  = {0,0,0};
  //
  Float_t origin[3] = {0,0,0};
  Float_t p[3];
  Float_t random[6];
  Int_t i, j, nt;
  //
  for (j=0;j<3;j++) origin[j]=fOrigin[j];
  if(fVertexSmear == kPerTrack) {
    Rndm(random,6);
    for (j = 0; j < 3; j++) {
	origin[j] += fOsigma[j]*TMath::Cos(2*random[2*j]*TMath::Pi())*
	    TMath::Sqrt(-2*TMath::Log(random[2*j+1]));
    }
  }

  Int_t nTracks = fReader->NextEvent(); 	
  if (nTracks == 0) {
      printf("\n No more events !!! !\n");
      return;
  }

  for (i = 0; i < nTracks; i++) {
      TParticle* iparticle = fReader->NextParticle();
      Double_t  theta = iparticle->Theta();
      Double_t  phi = iparticle->Phi();
      if (phi > TMath::Pi()) phi -= 2.*TMath::Pi();
      Double_t  pmom = iparticle->P();
      Double_t  pz   = iparticle->Pz();
      Double_t  e    = iparticle->Energy();
      Double_t  pt   = iparticle->Pt();
      Double_t  y;
      
      if ((e-pz) == 0) {
	  y = 20.;
      } else if ((e+pz) == 0.) {
	  y = -20.;
      } else {
	  y = 0.5*TMath::Log((e+pz)/(e-pz));	  
      }
       
      if(theta < fThetaMin || theta > fThetaMax ||
	 phi   < fPhiMin   || phi   > fPhiMax   ||
	 pmom  < fPMin     || pmom  > fPMax     ||
	 pt    < fPtMin    || pt    > fPtMax    ||
	 y     < fYMin     || y     > fYMax        )
      {
	  printf("\n Not selected %d %f %f %f %f %f", i, theta, phi, pmom, pt, y);
	  delete iparticle;
	  continue;
      }
      p[0] = iparticle->Px();
      p[1] = iparticle->Py();
      p[2] = iparticle->Pz();
      Int_t idpart = iparticle->GetPdgCode();
      if(fVertexSmear==kPerTrack) {
	  Rndm(random,6);
	  for (j = 0; j < 3; j++) {
	      origin[j]=fOrigin[j]
		  +fOsigma[j]*TMath::Cos(2*random[2*j]*TMath::Pi())*
		  TMath::Sqrt(-2*TMath::Log(random[2*j+1]));
	  }
      }
      SetTrack(fTrackIt,-1,idpart,p,origin,polar,0,kPPrimary,nt);
      delete iparticle;
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








