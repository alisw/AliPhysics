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

//
// Class for a simple description of the LHC.
// The LHC is described by two beams,
// the interaction regions and the
// beam loss processes.
// Run paramters can be set in order to simulate the time evolution
// of emittance, number of particles per bunch and luminosity.
// Author: Andreas Morsch
// andreas.morsch@cern.ch

#include "AliLHC.h"
#include "AliLhcIRegion.h"
#include "AliLhcProcess.h"
#include "AliLhcBeam.h"

ClassImp(AliLHC)

AliLHC::AliLHC()
{
// Constructor
    fIRegions   = new TList;
    fProcesses  = new TList;
    fNRegions   = 0;
    fNProcesses = 0;
    fBeams = new TObjArray(2);
    //PH    (*fBeams)[0] = 0;
    //PH    (*fBeams)[1] = 0;    
    fBeams->AddAt(0,0);
    fBeams->AddAt(0,1);    
    fTime = 0;
    fTimeMax = 0;
    fTimeA = 0;
}

AliLHC::AliLHC(const AliLHC& lhc)
    : TObject(lhc)
{
// copy constructor
}

AliLHC::~AliLHC()
{
// Destructor
    delete fIRegions;
    delete fProcesses;
    delete fBeams;
}

void AliLHC::AddIRegion(AliLhcIRegion *region)
{
//
//  Add region to list   
     fIRegions->Add(region);
     fNRegions++;
 }

void AliLHC::AddProcess(AliLhcProcess *process)
{
//
//  Add process to list   
     fProcesses->Add(process);
     fNProcesses++;
 }

void AliLHC::SetBeams(AliLhcBeam *beam1, AliLhcBeam *beam2 )
{

//
//  Set the beams   

    (*fBeams)[0] = beam1;
    (*fBeams)[1] = beam2;
}

  void AliLHC::Init()
{
// Initialisation
    fNt    = 0;
    fNmax  = Int_t(fTimeMax/fTimeStep); 
    fTimeA = new Float_t[fNmax];
    //
    Beam(0)->Init();
    Beam(1)->Init();
  
    TIter next(fIRegions);
    AliLhcIRegion *region;
    //
    // Loop over generators and initialize
    while((region = (AliLhcIRegion*)next())) {
	region->Init();
	region->SetMonitor(fNmax);
    }
    
    Beam(0)->SetMonitor(fNmax);
    Beam(1)->SetMonitor(fNmax);

    TIter nextp(fProcesses);
    AliLhcProcess *process;
    //
    // Loop over generators and initialize
    while((process = (AliLhcProcess*)nextp())) {
	process->Init();
	process->SetMonitor(fNmax);
    }  
}

 void AliLHC::EvolveTime()
{
//
// Simulate Time Evolution
//
    while (fTime <= fTimeMax) {
	printf("\n Time: %f %f", fTime, fTimeStep);
	//
        //  Processes
        //	
	TIter next(fProcesses);
	AliLhcProcess *process;
	//
	// Evolve for each process
	while((process = (AliLhcProcess*)next())) {
	    process->Evolve(fTimeStep);
	    process->Record();
	}  
	//
	// Update and Monitoring
	//
	TIter nextregion(fIRegions);
	AliLhcIRegion *region;
	//
	while((region = (AliLhcIRegion*)nextregion())) {
	  printf("\n Region: %s, Luminosity %10.3e", 
		 region->GetName(), region->Luminosity());
	  region->Update();
	  region->Record();
	}
	Beam(0)->Record();
	fTimeA[fNt] = fTime/3600.;
	fTime+=fTimeStep;
	fNt++;
    }
}

void AliLHC::Evaluate()
{
  // Evaluation of the results
  TIter nextregion(fIRegions);
  AliLhcIRegion *region;
  //
  // Loop over generators and initialize
  while((region = (AliLhcIRegion*)nextregion())) {
    region->DrawPlots();
  }
  
  TIter next(fProcesses);
  AliLhcProcess *process;
  //
  // Evolve for each process
  while((process = (AliLhcProcess*)next())) {
    process->DrawPlots();
  }
  
  Beam(0)->DrawPlots();
}
   
AliLHC& AliLHC::operator=(const  AliLHC & /*rhs*/)
{
// Assignment operator
    return *this;
}


