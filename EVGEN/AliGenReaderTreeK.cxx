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
Revision 1.4.4.1  2002/06/10 14:57:41  hristov
Merged with v3-08-02

Revision 1.5  2002/04/26 10:37:23  morsch
Method RewindEvent() added. (N. Carrer)

Revision 1.4  2002/03/22 08:25:33  morsch
TreeE connected correctly.

Revision 1.3  2001/12/12 11:21:37  morsch
Dummy copy constructor added.

Revision 1.2  2001/11/12 14:31:00  morsch
Memory leaks fixed. (M. Bondila)

Revision 1.1  2001/11/09 09:11:24  morsch
Realisation of AliGenReader that reads the kine tree (TreeK).

*/
#include <TFile.h>
#include <TTree.h>
#include <TParticle.h>

#include "AliGenReaderTreeK.h"
#include "AliStack.h"
#include "AliHeader.h"
#include "AliRun.h"

ClassImp(AliGenReaderTreeK);


AliGenReaderTreeK::AliGenReaderTreeK():AliGenReader() 
{
//  Default constructor
    fFileName       = NULL;
    fStack          = 0;
    fHeader         = 0;
    fNcurrent       = 0;
    fNparticle      = 0;
    fFile           = 0;
    fBaseFile       = 0;
    fTreeE          = 0;
}

AliGenReaderTreeK::AliGenReaderTreeK(const AliGenReaderTreeK &reader)
{
    ;
}


AliGenReaderTreeK::~AliGenReaderTreeK() 
{
// Destructor
    delete fTreeE;
}

void AliGenReaderTreeK::Init() 
{
// Initialization
// Connect base file and file to read from

    TTree *ali = gAlice->TreeE();
    if (ali) {
	fBaseFile = ali->GetCurrentFile();
    } else {
	printf("\n Warning: Basefile cannot be found !\n");
    }
    if (!fFile) fFile  = new TFile(fFileName);
}

Int_t AliGenReaderTreeK::NextEvent() 
{
//  Read the next event  
//  cd to file with old kine tree    
    if (!fBaseFile) Init();
    fFile->cd();
//  Connect header tree
    if (!fTreeE) fTreeE = (TTree*)gDirectory->Get("TE");
    if (fHeader) delete fHeader;
    fHeader = 0;
    fTreeE->SetBranchAddress("Header", &fHeader);
//  Get next event
    fTreeE->GetEntry(fNcurrent);
//  Connect Stack
    if (fStack) delete fStack;
    fStack = fHeader->Stack();
    fStack->GetEvent(fNcurrent);
//  cd back to base file
    fBaseFile->cd();
//
    fNcurrent++;
    fNparticle = 0;
    Int_t ntrack =  fStack->GetNtrack();
    printf("\n Next event contains %d particles", ntrack);
//    
    return  ntrack;
}

TParticle* AliGenReaderTreeK::NextParticle() 
{
//  Return next particle
    TParticle* part = fStack->Particle(fNparticle);
    fNparticle++;
    return part;
}

void AliGenReaderTreeK::RewindEvent()
{
  // Go back to the first particle of the event
  fNparticle = 0;
}


AliGenReaderTreeK& AliGenReaderTreeK::operator=(const  AliGenReaderTreeK& rhs)
{
// Assignment operator
    return *this;
}



