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
*/

#include <iostream.h>

#include "AliGenKine.h"
#include "AliMCProcess.h"
#include "AliMC.h"
#include "AliRun.h"
#include "AliStack.h"

#include <TROOT.h>
#include <TTree.h>
#include <TDirectory.h>
#include <TDatabasePDG.h>
#include <TFile.h>
#include <TH1.h>
#include <TArrayI.h>
#include <TClonesArray.h>
#include <TObjArray.h>
#include <TParticle.h>
#include <stdlib.h>

 ClassImp(AliGenKine)
     AliGenKine::AliGenKine()
	 :AliGenerator(-1)
{
//  Constructor
    fName           = "Kine";
    fTitle          = "Primaries from ext. File";
    fFileName       = NULL;
    fStack          = 0;
    fNcurrent       = 0;
//
//  Read all particles
    fNpart       = -1;
    fFile        =  0;
    fBaseFile    =  0;
}

AliGenKine::AliGenKine(Int_t npart)
    :AliGenerator(npart)
{
//  Constructor
    fName           = "Kine";
    fTitle          = "Primaries from ext. File";
    fFileName       = NULL;
    fStack          = 0;
    fNcurrent       = 0;
    fFile           = 0;
    fBaseFile       = 0;
}

//____________________________________________________________
AliGenKine::~AliGenKine()
{
// Destructor
    
}

//____________________________________________________________

//____________________________________________________________
void AliGenKine::Generate()
{
    Float_t polar[3], vpos[3], pmom[3];
//
//  Connect file and get next event  
//
    if (!fBaseFile) {
	TTree *ali = gAlice->TreeE();
	if (ali) fBaseFile = ali->GetCurrentFile();
    }

    if (!fFile) {
        fFile = new TFile(fFileName);
    }

//  cd to file with old kine tree    
    if (fStack) delete fStack;
    fStack = new AliStack(1000);
    fFile->cd();
//  Connect treeE
    TTree* treeE = (TTree*)gDirectory->Get("TE");
    treeE->SetBranchAddress("Stack", &fStack);
//  Get next event
    treeE->GetEntry(fNcurrent);
    fStack->GetEvent(fNcurrent);
//  cd back to base file
    fBaseFile->cd();
//
//  Read Particles
//
    Int_t ntr;

    for (Int_t i = 0; i < fStack->GetNtrack(); i++)
    {
	TParticle* part = fStack->Particle(i);

	
	Int_t pdg = part->GetPdgCode();
//	if (pdg == -1) continue;
	
	Int_t parent = part->GetFirstMother();
	Float_t tof  = part->T();
	
	vpos[0] = part->Vx();
	vpos[1] = part->Vy();	
	vpos[2] = part->Vz();

	pmom[0] = part->Px();
	pmom[1] = part->Py();
	pmom[2] = part->Pz();

	printf("\n %d %d %d %f", i, fStack->GetNtrack(), part->GetPdgCode(), pmom[2]);
	gAlice->SetTrack(fTrackIt, parent, pdg, pmom, vpos, polar, 
			 tof,  kPPrimary, ntr);
	gAlice->KeepTrack(ntr);
    }
    gAlice->SetHighWaterMark(ntr);
    fNcurrent++;
    
// cd back to output file

}








