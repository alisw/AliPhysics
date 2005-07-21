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
    // Get event from runloader
    fRunLoader->GetEvent(event);
    // Get the stack
    AliStack* stack = fRunLoader->Stack();
    // Number of primaries
    Int_t nt = stack->GetNprimary();
    // Get cuts set by user and header
    Double_t ptMin = ((AliJetKineReaderHeader*) fReaderHeader)->GetPtCut();
    fAliHeader = fRunLoader->GetHeader();
        
    // Loop over particles
    Int_t* flag  = new Int_t[nt];
    for (Int_t it = 0; it < nt; it++) {
	TParticle *part = stack->Particle(it);
	Int_t   status  = part->GetStatusCode();
	Int_t   pdg     = TMath::Abs(part->GetPdgCode());
	Float_t pt      = part->Pt(); 
	
	// Skip non-final state particles, neutrinos and particles with pt < pt_min 
	
	if (
	    (status != 1)            
	    || (pdg == 12 || pdg == 14 || pdg == 16) 
	    || (pt < ptMin)             
	    ) continue; 


	Float_t p       = part->P();
	Float_t p0      = p;
	Float_t eta     = part->Eta();
	Float_t phi     = part->Phi();

        // Fast simulation of TPC if requested
	if (((AliJetKineReaderHeader*)fReaderHeader)->FastSimTPC()) {
	    // Charged particles only
	    Float_t charge = TDatabasePDG::Instance()->GetParticle(pdg)->Charge();
	    if (charge == 0)               continue;
	    // Simulate efficiency
	    if (!Efficiency(p0, eta, phi)) continue;
	    // Simulate resolution
	    p = SmearMomentum(4, p0);
	} // Fast TPC
	
        // Fast simulation of EMCAL if requested
	if (((AliJetKineReaderHeader*)fReaderHeader)->FastSimEMCAL()) {
	    Float_t charge = TDatabasePDG::Instance()->GetParticle(pdg)->Charge();
	    // Charged particles
	    if (charge != 0) {
		// Simulate efficiency
		if (!Efficiency(p0, eta, phi)) continue;
		// Simulate resolution
		p = SmearMomentum(4, p0);
	    } // charged
	    // Neutral particles
	    // Exclude K0L, n, nbar
	    if (pdg == kNeutron || pdg == kK0Long) continue;
	} // Fast EMCAL
	
	fMass = part->GetCalcMass();
	fPdgC = part->GetPdgCode();
	// Fill momentum array
	Float_t r  = p/p0;
//	    	printf("smeared %13.3f %13.3f %13.3f\n", p0, p, r);
	
	Float_t px = r * part->Px(); 
	Float_t py = r * part->Py(); 
	Float_t pz = r * part->Pz();
	Float_t m  = part->GetMass();
	Float_t e  = TMath::Sqrt(px * px + py * py + pz * pz + m * m);
	
	new ((*fMomentumArray)[goodTrack]) 
	    TLorentzVector(px, py, pz, e);
	goodTrack++;
  }
    fSignalFlag.Set(goodTrack,flag);
    printf("\nNumber of good tracks in event %d= %d \n",event,goodTrack);
}


Float_t AliJetKineReader::SmearMomentum(Int_t ind, Float_t p)
{
//
//  The relative momentum error, i.e. 
//  (delta p)/p = sqrt (a**2 + (b*p)**2) * 10**-2,
//  where typically a = 0.75 and b = 0.16 - 0.24 depending on multiplicity
//  (the lower value is for dn/d(eta) about 2000, and the higher one for 8000)
//
//  If we include information from TRD b will be by a factor 2/3 smaller.
//
//  ind = 1: high multiplicity
//  ind = 2: low  multiplicity
//  ind = 3: high multiplicity + TRD
//  ind = 4: low  multiplicity + TRD

    Float_t pSmeared;
    Float_t a = 0.75;
    Float_t b = 0.12;

    if (ind == 1) b = 0.12;
    if (ind == 2) b = 0.08;
    if (ind == 3) b = 0.12;    
    if (ind == 4) b = 0.08;    
    
    Float_t sigma =  p * TMath::Sqrt(a * a + b * b * p * p)*0.01;
    pSmeared = p + gRandom->Gaus(0., sigma);
    return pSmeared;
}
Bool_t AliJetKineReader::Efficiency(Float_t p, Float_t /*eta*/, Float_t phi)
{
//
// Fast simulation of geometrical acceptance and tracking efficiency
//
//  Tracking
    Float_t eff = 0.99;
    if (p < 0.5) eff -= (0.5-p)*0.2/0.3;
// Geometry
    if (p > 0.5) {
	phi *= 180. / TMath::Pi();
	// Sector number 0 - 17
	Int_t isec  = Int_t(phi / 20.);
	// Sector centre
	Float_t phi0 = isec * 20. + 10.;
	Float_t phir = TMath::Abs(phi-phi0);
	// 2 deg of dead space
	if (phir > 9.) eff = 0.;
    }
    if (gRandom->Rndm() > eff) {
	return kFALSE;
    } else {
	return kTRUE;
    }
}

