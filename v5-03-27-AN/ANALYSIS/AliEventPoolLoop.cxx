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
 

// Realisation of an AliVEventPool which allows the user to
// run the analysis in a loop, i.e. passing several times over 
// the same event chain.
// Author Andreas Morsch
// andreas.morsch@cern.ch


#include "AliEventPoolLoop.h"
#include <TChain.h>
#include <TFile.h>
#include <TObjArray.h>

ClassImp(AliEventPoolLoop)


////////////////////////////////////////////////////////////////////////

AliEventPoolLoop::AliEventPoolLoop():
    AliVEventPool(),
    fMaxIterations(0),
    fNIteration(1),
    fChainClone(0)
{
  // Default constructor
}

AliEventPoolLoop::AliEventPoolLoop(Int_t nit):
    AliVEventPool(),
    fMaxIterations(nit),
    fNIteration(1),
    fChainClone(0)
{
  // Default constructor
}

AliEventPoolLoop::AliEventPoolLoop(const char* name, const char* title):
    AliVEventPool(name, title),
    fMaxIterations(0),
    fNIteration(1),
    fChainClone(0)
{
  // Constructor
}


AliEventPoolLoop::AliEventPoolLoop(const AliEventPoolLoop& obj):
    AliVEventPool(obj),
    fMaxIterations(obj.fMaxIterations),
    fNIteration(obj.fNIteration),
    fChainClone(0)
{
    // Copy constructor
}

AliEventPoolLoop& AliEventPoolLoop::operator=(const AliEventPoolLoop& other)
{
// Assignment operator
    AliVEventPool::operator=(other);
    fMaxIterations = other.fMaxIterations;
    fNIteration    = other.fNIteration;
    return *this;
}


void AliEventPoolLoop::Init()
{
// Initialisation

    fMaxIterations = 0;
    fNIteration    = 1;
}

TChain* AliEventPoolLoop::GetNextChain()
{
    // Get the next chain
    if (fNIteration > fMaxIterations) {
	return (0);
    } else {
	fNIteration++;
	/*
	if (fChainClone) {
	  fChainClone->Reset();
	  new (fChainClone) TChain(fChain->GetName());
	} else {
	  fChainClone = new TChain(fChain->GetName());
	}
	TObjArray* files = fChain->GetListOfFiles();
	Int_t n = files->GetEntriesFast();
	for (Int_t i = 0; i < n; i++) {
	  TFile* file = (TFile*) (files->At(i));
	  fChainClone->AddFile(file->GetTitle());
	  printf("Adding %s %s %s\n", fChainClone->GetName(), file->GetName(), file->GetTitle());

	}
	*/
	fChainClone = (TChain*) (fChain->Clone());
	return fChainClone;
    }
}

void  AliEventPoolLoop::GetCurrentBin(Float_t* /*bin*/)
{
    //
}

Int_t AliEventPoolLoop::GetDimension()
{
    //
    return (0);
}

