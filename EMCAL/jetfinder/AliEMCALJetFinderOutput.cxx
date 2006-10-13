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

//_________________________________________________________________________
//  Output object for jetfinder
// --
//*-- Author: 	Renan Cabrera (LBL)
//		Mark Horner (LBL/UCT)
// --


#include <stdio.h>
#include <TParticle.h>
#include <TTree.h>

//.....................
#include "AliEMCALJet.h"
#include "AliEMCALParton.h"
#include "AliEMCALJetFinderOutput.h"
#include "AliEMCALJetFinderInput.h"

ClassImp(AliEMCALJetFinderOutput)

//________________________________________________________________________
AliEMCALJetFinderOutput::AliEMCALJetFinderOutput() :
  fJetsArray(0),fPartonsArray(0),fNPartons(0),fNJets(0),
  fParticlesArray(0),fNParticles(0),fNMaxJets(0),fNMaxParticles(0),
  fNMaxPartons(0),fDebug(0),fInitialised(kFALSE)
{ 
	// Default constructor
	
	fNMaxJets=10;
	fNMaxParticles=2000;
	fNMaxPartons=4;
	fInitialised=kFALSE;
	fDebug=0;
	fNPartons=0;
	fNJets=0;    
	fNParticles=0;

	fJetsArray=0;
	fParticlesArray=0;
	fPartonsArray=0;

if (fDebug>0) Info("AliEMCALJetFinderOutput","Beginning Constructor");

} 

//________________________________________________________________________
AliEMCALJetFinderOutput::AliEMCALJetFinderOutput(const AliEMCALJetFinderOutput& jfo) 
  : TObject(jfo), fJetsArray(jfo.fJetsArray),fPartonsArray(jfo.fPartonsArray),fNPartons(jfo.fNPartons),
    fNJets(jfo.fNJets), fParticlesArray(jfo.fParticlesArray),fNParticles(jfo.fNParticles),
    fNMaxJets(jfo.fNMaxJets),fNMaxParticles(jfo.fNMaxParticles),fNMaxPartons(jfo.fNMaxPartons),
    fDebug(jfo.fDebug),fInitialised(jfo.fInitialised)
{ 
  //copy ctor
}

//________________________________________________________________________
void AliEMCALJetFinderOutput::InitArrays()
{
	// Initialise arrays - legacy from TClones days
if (fDebug>1) Info("AliEMCALJetFinderOutput","Beginning InitArrays");
 fParticlesArray=new TClonesArray("TParticle",fNMaxParticles);
 fPartonsArray=new TClonesArray("AliEMCALParton",fNMaxPartons);
 fJetsArray=new TClonesArray("AliEMCALJet",fNMaxJets);
 //fJetsArray->BypassStreamer(kFALSE);
 fInitialised=1;
}

//_______________________________________________________________________
AliEMCALJetFinderOutput::~AliEMCALJetFinderOutput()
{
	// Default destrucotr
if (fDebug>0) Info("~AliEMCALJetFinderOutput","Beginning Destructor");
 delete fParticlesArray;
 delete fPartonsArray;
 delete fJetsArray;
}

//_______________________________________________________________________
void AliEMCALJetFinderOutput::Reset(AliEMCALJetFinderResetType_t resettype)
{
// Reset stored data
	
if (fDebug>1) Info("Reset","Beginning Reset");
if (!fInitialised) InitArrays();
 if (	resettype == kResetAll ||
	resettype == kResetJets||
	resettype == kResetData ){
	 fNJets = 0;
	 if (fJetsArray)
	   fJetsArray->Clear();
 }
 if (   resettype == kResetAll ||
        resettype == kResetPartons||              
        resettype == kResetData ){
	 fNPartons = 0;
	 if (fPartonsArray)
	   fPartonsArray->Clear();
 }
 if (   resettype == kResetAll ||    
        resettype == kResetParticles||              
        resettype == kResetData ){
	 fNParticles = 0;
	 if (fParticlesArray)
	   fParticlesArray->Clear();
 }
}
//________________________________________________________________________
void AliEMCALJetFinderOutput::AddJet(AliEMCALJet* jet)
{
// Add a jet to the array
	
if (fDebug>1) Info("AddJet","Beginning AddJet");
if (!fInitialised) InitArrays();


  	if (fNJets < fNMaxJets){
		new(fJetsArray->AddrAt(fNJets))   AliEMCALJet( *jet );
		fNJets++;
	}else
	{
		Error("AddJet","Cannot AddJet - maximum exceeded");
                }
   
}


//_______________________________________________________________________
void AliEMCALJetFinderOutput::AddParton(AliEMCALParton* parton)
{
//Add a parton to the array
	
if (fDebug>1) Info("AddParton","Beginning AddParton");
if (!fInitialised) InitArrays();

	if (fNPartons < fNMaxPartons){
		new(fPartonsArray->AddrAt(fNPartons))  AliEMCALParton( *parton );
		fNPartons++;
	}else
	{
                Error("AddParton","Cannot AddParton - maximum exceeded");
	}
 
}

//_______________________________________________________________________
void AliEMCALJetFinderOutput::AddParticle(TParticle* particle)
{
//Add a particle tot he array
	
if (fDebug>1) Info("AddParticle","Beginning AddParticle");
if (!fInitialised) InitArrays();

	if (fNParticles < fNMaxParticles){
		new(fParticlesArray->AddrAt(fNParticles))  TParticle( *particle );
		fNParticles++;
	}else
	{
		Error("AddParticle","Cannot AddParticle - maximum exceeded");
                }
}

//______________________________________________________________________
AliEMCALJet* AliEMCALJetFinderOutput::GetJet(Int_t jetID)
{
	// return a jet
if (fDebug>1) Info("GetJet","Beginning GetJet");
	
  if (jetID >= fNJets) return 0;
  return (AliEMCALJet*)fJetsArray->At(jetID);
  
}

//______________________________________________________________________
AliEMCALParton* AliEMCALJetFinderOutput::GetParton(Int_t partonID)
{
	//return a parton
if (fDebug>1) Info("GetParton","Beginning GetParton");

  if (partonID >= fNPartons) return 0;
  return (AliEMCALParton*) fPartonsArray->At(partonID);
}

//______________________________________________________________________
TParticle* AliEMCALJetFinderOutput::GetParticle(Int_t particleID)
{
//return a particle
	
if (fDebug>1) Info("GetParticle","Beginning GetParticle");

  if (particleID >= fNParticles) return 0;
  return (TParticle*) fParticlesArray->At(particleID);

}

