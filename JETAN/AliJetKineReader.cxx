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

// 
// Jet kinematics reader 
// MC reader for jet analysis
// Author: Andreas Morsch 
// andreas.morsch@cern.ch
//

// From root ...
#include <TClonesArray.h>
#include <TPDGCode.h>
#include <TParticle.h>
#include <TParticlePDG.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TSystem.h>
#include <TDatabasePDG.h>
#include <TRandom.h>
// From AliRoot ...
#include "AliJetKineReaderHeader.h"
#include "AliJetKineReader.h"
#include "AliRunLoader.h"
#include "AliStack.h"

ClassImp(AliJetKineReader)

AliJetKineReader::AliJetKineReader()
{
  // Constructor
  fReaderHeader = 0x0;
  fRunLoader    = 0x0;
  fMass         = 0;
  fPdgC         = 0;
}

//____________________________________________________________________________

AliJetKineReader::~AliJetKineReader()
{
  // Destructor
  delete fReaderHeader;
}

//____________________________________________________________________________

void AliJetKineReader::OpenInputFiles()
{
    // Opens the input file using the run loader
    const char* dirName = fReaderHeader->GetDirectory();
    char path[256];
    sprintf(path, "%s/galice.root",dirName);
    fRunLoader = AliRunLoader::Open(path);
    fRunLoader->LoadKinematics();
    fRunLoader->LoadHeader(); 
    
    Int_t nMax = fRunLoader->GetNumberOfEvents();
    printf("\nTotal number of events = %d", nMax);
    
  // set number of events in header
    if (fReaderHeader->GetLastEvent() == -1)
	fReaderHeader->SetLastEvent(nMax);
    else {
	Int_t nUsr = fReaderHeader->GetLastEvent();
	fReaderHeader->SetLastEvent(TMath::Min(nMax, nUsr));
    }
}

//____________________________________________________________________________

void AliJetKineReader::FillMomentumArray(Int_t event)
{
//
// Fill momentum array for event
//
    Int_t goodTrack = 0;
    // Clear array
    ClearArray();

    fRunLoader->GetEvent(event);
    AliStack* stack = fRunLoader->Stack();
    Int_t nt = stack->GetNprimary();
    // Get cuts set by user
    Double_t ptMin = ((AliJetKineReaderHeader*) fReaderHeader)->GetPtCut();
    fAliHeader = fRunLoader->GetHeader();
        
    // Loop over particles
    for (Int_t it = 0; it < nt; it++) {
	TParticle *part = stack->Particle(it);
	Int_t   status  = part->GetStatusCode();
	Int_t   pdg     = TMath::Abs(part->GetPdgCode());
	Float_t pt      = part->Pt(); 
	
	if (
	    (status != 1)            
	    || (pdg == 12 || pdg == 14) 
	    || (pt < ptMin)             
	    ) continue; 

	if (((AliJetKineReaderHeader*)fReaderHeader)->FastSimTPC()) {
	    Float_t charge = TDatabasePDG::Instance()->GetParticle(pdg)->Charge();
	      if (
		  (charge == 0)  
		  || (gRandom->Rndm() < 0.1) 
		  ) continue;
	}
	
	fMass = part->GetCalcMass();
	fPdgC = part->GetPdgCode();
	// Fill momentum array
	
	new ((*fMomentumArray)[goodTrack]) 
	    TLorentzVector(part->Px(),part->Py(),part->Pz(), part->Energy());
	goodTrack++;
  }
    printf("\nNumber of good tracks %d \n", goodTrack);
}


