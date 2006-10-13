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


//________________________________________________________________________
// Initial JetFinder input object
// --
//*-- Author: Mark Horner (LBL/UCT)
// --
// --
// --



#include "AliEMCALJetFinderInput.h"
#include "AliEMCALDigit.h"
#include "AliEMCALGeometry.h"
#include "AliEMCALParton.h"
class TClonesArray;


ClassImp(AliEMCALJetFinderInput)

AliEMCALJetFinderInput::AliEMCALJetFinderInput(): 
  fDigitsArray(0),fNDigits(0),fNMaxDigits(0),fTracksArray(0),fNTracks(0),fNMaxTracks(0),
  fPartonsArray(0),fNPartons(0),fNMaxPartons(0),fParticlesArray(0),
  fNParticles(0),fNMaxParticles(0),fDebug(0),fInitialised(0)
{
// Default constructor
// Defines all TClonesArrays  
// Creates full array of Digits - all other arrays will have objects created as necessary 

  //if (fDebug>0) Info("AliEMCALJetFinderInput","Beginning Constructor");	
  fDebug=0;
  fInitialised = kFALSE;
  fNDigits = 19152;     // This is the number of digits
  fNMaxDigits = fNDigits;
  fNMaxTracks = 3000;
  fNMaxParticles = 2000;
  fNMaxPartons = 4;
  fNTracks        = 0;
  fNPartons       = 0;
  fNParticles     = 0;

}


AliEMCALJetFinderInput::AliEMCALJetFinderInput(const AliEMCALJetFinderInput& ji)
  : TObject(ji), fDigitsArray(ji.fDigitsArray),fNDigits(ji.fNDigits),fNMaxDigits(ji.fNMaxDigits),
    fTracksArray(ji.fTracksArray),fNTracks(ji.fNTracks),fNMaxTracks(ji.fNMaxTracks),
    fPartonsArray(ji.fPartonsArray),fNPartons(ji.fNPartons),fNMaxPartons(ji.fNMaxPartons),
    fParticlesArray(ji.fParticlesArray),fNParticles(ji.fNParticles),fNMaxParticles(ji.fNMaxParticles),
    fDebug(ji.fDebug),fInitialised(ji.fInitialised)
{
  //copy ctor
}

void AliEMCALJetFinderInput::InitArrays()
{
if (fDebug>0) Info("AliEMCALJetFinderInput","Beginning InitArrays");
fNDigits = 19152;     // This is the number of digits
fNMaxDigits = fNDigits;
fNMaxTracks = 3000;
fNMaxParticles = 2000;
fNMaxPartons = 4;

fDigitsArray  = new TClonesArray ("AliEMCALDigit",fNDigits);     // This is the digits array for the EMCAL
fTracksArray  = new TClonesArray("TParticle",fNMaxTracks);
fPartonsArray = new TClonesArray("AliEMCALParton",fNMaxPartons);
fParticlesArray       = new TClonesArray ("TParticle",fNMaxParticles);

for (Int_t i = 0; i < 19152; i++)
{
	AliEMCALDigit tempdigit(0,0,i+1,0,1.0,-1);
	tempdigit.SetAmp(0);
//	 AliEMCALDigit(0,0,i+1,0,1.0,-1)
	new((*fDigitsArray)[i]) AliEMCALDigit(tempdigit);	// Default digit amplitude is zero
}
fNTracks	= 0;
fNPartons	= 0;
fNParticles	= 0;
fInitialised = kTRUE;
}


AliEMCALJetFinderInput::~AliEMCALJetFinderInput()
{
// Destructor -

if (fDebug>0) Info("~AliEMCALJetFinderInput","Beginning Destructor");	
if (fInitialised)
 {
       fDigitsArray->Delete();
       fTracksArray->Delete();
       fPartonsArray->Delete();
       fParticlesArray->Delete();
       delete fDigitsArray;
       delete fTracksArray;
       delete fPartonsArray;
       delete fParticlesArray;
}
		
}

void AliEMCALJetFinderInput::Reset(AliEMCALJetFinderResetType_t resettype)
{
// Clears the contents of all the arrays and returns to the state we were in 
// after the constructor
if (fDebug>1) Info("Reset","Beginning Reset");	
  
  if (!fInitialised) InitArrays();	

  if ( 	resettype == kResetData   ||	
	resettype == kResetDigits ||
	resettype == kResetAll   ){
		if (fInitialised)
 		{	
			fDigitsArray->Delete();			
	     		fTracksArray->Delete();
		 	fPartonsArray->Delete();
		        fParticlesArray->Delete();
		}
	  	for (Int_t i = 0; i < 19152; i++)
		{
			AliEMCALDigit tempdigit(0,0,i+1,0,1.0,-1);
			tempdigit.SetAmp(0);
			new((*fDigitsArray)[i]) AliEMCALDigit(tempdigit);	// Default digit amplitude is zero
		}
  }
  if (  resettype == kResetData   ||    
        resettype == kResetTracks ||
        resettype == kResetAll   ){
	  fNTracks = 0;
  }
 if (  resettype == kResetData    ||
       resettype == kResetPartons ||
       resettype == kResetAll   ){
	 fNPartons=0;
 }
 if (  resettype == kResetData      ||
       resettype == kResetParticles ||
       resettype == kResetAll   ){
	 fNParticles = 0;
 }

}
void AliEMCALJetFinderInput::AddEnergyToDigit(Int_t digitID,Int_t denergy)
{
// Adds energy to the digit specified - Note these are tower IDs and
// so start at 1!
	  	
if (fDebug>15) Info("AddEnergyToDigit","Beginning AddEnergyToDigit");
 if (!fInitialised) InitArrays();	

 Int_t amp = 0;
 AliEMCALDigit *mydigit = (AliEMCALDigit*)(*fDigitsArray)[digitID];
 amp = mydigit->GetAmp();
 mydigit->SetAmp(amp+denergy);
}	

void AliEMCALJetFinderInput::AddTrack(TParticle track)
{
// Adds a TParticle to the track array
	  	
if (fDebug>5) Info("AddTrack","Beginning AddTrack");
 
 	if (!fInitialised) InitArrays();	
	if (fNTracks < fNMaxTracks){  
		new((*fTracksArray)[fNTracks]) TParticle(track);
		fNTracks++;
		if (fDebug>5) Info("AddTracks","Added Tracks %i",fNTracks);	
	}else
	{
		Error("AddTracks","Cannot AddTrack - maximum exceeded");
	}

}

void AliEMCALJetFinderInput::AddTrack(TMCParticle *track)
{
// Adds a TParticle to the particle array
	
if (fDebug>5) Info("AddTrack","Beginning AddTrack");	

 if (!fInitialised) InitArrays();	
 if (fNTracks < fNMaxTracks){  
   new((*fTracksArray)[fNTracks]) 
     TParticle(track->GetKF(), track->GetKS(), track->GetParent(), 0,
	       track->GetFirstChild(), track->GetLastChild(),
	       track->GetPx(), track->GetPy(), track->GetPz(), 
	       track->GetEnergy(), 
	       track->GetVx(), track->GetVy(), track->GetVz(), 
	       track->GetTime());
   fNTracks++;
   if (fDebug>5) Info("AddTrack","Added Track %i",fNTracks);	
 } else {
   Error("AddTrack","Cannot AddTrack - maximum exceeded");
 }	
}


void AliEMCALJetFinderInput::AddParton(AliEMCALParton *parton)
{
// Adds an AliEMCALParton to the parton array 
	
if (fDebug>5) Info("AddParton","Beginning AddParton");

 	if (!fInitialised) InitArrays();	
	if (fNPartons < fNMaxPartons){  
		new((*fPartonsArray)[fNPartons]) AliEMCALParton(*parton);
		fNPartons++;
		if (fDebug>5) Info("AddParton","Added Parton %i",fNPartons);	
	}else
	{
		Error("AddParton","Cannot AddParton - maximum exceeded");
	}
}

void AliEMCALJetFinderInput::AddParticle(TParticle *particle)
{
// Adds a TParticle to the particle array
	
if (fDebug>5) Info("AddParticle","Beginning AddParticle");	

 	if (!fInitialised) InitArrays();	
	if (fNParticles < fNMaxParticles){  
		new((*fParticlesArray)[fNParticles]) TParticle(*particle);
		fNParticles++;
		if (fDebug>5) Info("AddParticle","Added Particle %i",fNParticles);	
	}else
	{
		Error("AddParticle","Cannot AddParticle - maximum exceeded");
	}

}

void AliEMCALJetFinderInput::AddParticle(TMCParticle *particle)
{
// Adds a TParticle to the particle array
	
if (fDebug>5) Info("AddParticle","Beginning AddParticle");	

 if (!fInitialised) InitArrays();	
 if (fNParticles < fNMaxParticles){  
   new((*fParticlesArray)[fNParticles]) 
     TParticle(particle->GetKF(), particle->GetKS(), particle->GetParent(), 0,
	       particle->GetFirstChild(), particle->GetLastChild(),
	       particle->GetPx(), particle->GetPy(), particle->GetPz(), 
	       particle->GetEnergy(), 
	       particle->GetVx(), particle->GetVy(), particle->GetVz(), 
	       particle->GetTime());
   fNParticles++;
   if (fDebug>5) Info("AddParticle","Added Particle %i",fNParticles);	
 } else {
   Error("AddParticle","Cannot AddParticle - maximum exceeded");
 }	
}


AliEMCALDigit* AliEMCALJetFinderInput::GetDigit(Int_t digitID)
{
// Returns a pointer to an AliEMCALDigit with the given ID
	
if (fDebug>15) Info("GetDigit","Beginning GetDigit");

    if (digitID >= fNDigits) return 0;
    return (AliEMCALDigit*)((*fDigitsArray)[digitID]);
}

TParticle* AliEMCALJetFinderInput::GetTrack(Int_t trackID)
{
// Returns a pointer to a TParticle specified by the ID
	
if (fDebug>15) Info("GetTrack","Beginning GetTrack with trackID %i",trackID);

    if (trackID >= fNTracks) return 0;
    return (TParticle*)((*fTracksArray)[trackID]);
}

AliEMCALParton* AliEMCALJetFinderInput::GetParton(Int_t partonID)
{
// Returns a pointer to an AliEMCALParton  specified by the ID
	
if (fDebug>15) Info("GetParton","Beginning GetParton");

    if (partonID >= fNPartons) return 0;
    return (AliEMCALParton*)((*fPartonsArray)[partonID]);
}

TParticle* AliEMCALJetFinderInput::GetParticle(Int_t particleID)
{
// Returns a pointer to a TParticle specified by the ID
	
if (fDebug>15) Info("GetParticle","Beginning GetParticle");

    if (particleID >= fNParticles) return 0;
    return (TParticle*)((*fParticlesArray)[particleID]);
}


	
