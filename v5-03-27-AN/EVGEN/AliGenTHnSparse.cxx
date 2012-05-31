/**************************************************************************
 * Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
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

//-----------------------------------------------------------------------
// Particle generator according to 4 correlated variables : here
// z, ptot, r, theta. The input is a THnSparse object included in
// the root file (path and name to be set via the SetTHnSparse method).
// This class is similar to AliGenFunction.
//-----------------------------------------------------------------------
// Author : X. Lopez - LPC Clermont (fr)
//-----------------------------------------------------------------------
/*
  Example for generation :
  	AliGenTHnSparse *gener = new AliGenTHnSparse();
	gener->SetNumberParticles(10);
	gener->SetPart(13,kTRUE); // for generating id 13 and -13
	gener->SetThnSparse("file_name","thn_name");
	gener->Init();
*/

#include <TRandom.h>
#include <TFile.h>
#include "THnSparse.h"

#include "AliGenTHnSparse.h"

// NEW
#include "AliRun.h"
#include "AliLog.h"
#include "AliGenEventHeader.h"
//

ClassImp(AliGenTHnSparse)

//_______________________________________________________________________
AliGenTHnSparse::AliGenTHnSparse():
  AliGenerator(),
  fHn(0),
  fFile(0),
  fIpart(0),
  fBoth(kFALSE)
{
    // Default constructor
    SetNumberParticles(1);
}

//_______________________________________________________________________
AliGenTHnSparse::AliGenTHnSparse(const AliGenTHnSparse& func):
  AliGenerator(),
  fHn(func.fHn),
  fFile(func.fFile),
  fIpart(func.fIpart),
  fBoth(func.fBoth)
{
    // Copy constructor
    SetNumberParticles(1);
}

//_______________________________________________________________________
AliGenTHnSparse & AliGenTHnSparse::operator=(const AliGenTHnSparse& func)
{
    // Assigment operator
    if(&func == this) return *this;
    fHn  = func.fHn;
    fFile  = func.fFile;
    fIpart  = func.fIpart;
    fBoth   = func.fBoth;
    
    return *this;
}

//_______________________________________________________________________
AliGenTHnSparse::~AliGenTHnSparse()
{
    // Destructor
    delete fFile;
}

//_______________________________________________________________________
void AliGenTHnSparse::Generate()
{
    // Generate Npart of id Ipart
    Int_t naccepted =0;
    
    Double_t rand[4]; //  z, ptot, r, theta
    Float_t pos[3], phi, ptot, theta, pt, z, r;
    Float_t mom[3];
    Int_t pdg = fIpart;
  
    for (Int_t ipart = 0; ipart < fNpart && naccepted<fNpart; ipart++) {

	fHn->GetRandom(rand);
	z=rand[0];
	ptot=rand[1];
	r=rand[2];
	theta=rand[3];

// Phi: same for position and momemtum
 
	phi=(-180+gRandom->Rndm()*360)*TMath::Pi()/180;

// position at production
	
	pos[0] = r*TMath::Cos(phi);
	pos[1] = r*TMath::Sin(phi);
	pos[2] = z;
	
// momentum at production

	pt     = ptot*TMath::Sin(theta);
	mom[0] = pt*TMath::Cos(phi); 
	mom[1] = pt*TMath::Sin(phi); 
	mom[2] = ptot*TMath::Cos(theta);

// propagation

	Float_t polarization[3] = {0,0,0};
	Int_t nt;

// Part and anti-part

	if(fBoth){
	    Double_t sign = gRandom->Rndm();
	    if(sign < 0.5) pdg = -fIpart;
	    else pdg = fIpart;
	}

	PushTrack(fTrackIt,-1,pdg,mom, pos, polarization,0,kPPrimary,nt);
	naccepted++;
  }

    AliGenEventHeader* header = new AliGenEventHeader("THn");
    gAlice->SetGenEventHeader(header);
    return;

}

//_______________________________________________________________________
void AliGenTHnSparse::Init()
{
    
    // Initialisation, check consistency of selected file
    printf("************ AliGenTHnSparse ****************\n");
    printf("*********************************************\n");
    if (!fHn){
	AliFatal("THnSparse file not specified");
    }

    return;
}

//_______________________________________________________________________
void AliGenTHnSparse::SetThnSparse(char *file_name, char *thn_name)
{

    // Open the file and get object
  fFile = new TFile(file_name);
  fHn = (THnSparseF*)(fFile->Get(thn_name));

}
