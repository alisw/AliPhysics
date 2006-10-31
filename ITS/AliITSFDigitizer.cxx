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

//////////////////////////////////////////////////////////////////
// Class for fast reconstruction of recpoints                   //
//                                                              //
//                                                              //
//////////////////////////////////////////////////////////////////

#include <stdlib.h>
#include <TTree.h>

#include <AliRun.h>
#include <AliRunLoader.h>
#include <AliLoader.h>
#include <AliRunDigitizer.h>

#include "AliITSFDigitizer.h"
#include "AliITSgeom.h"
#include "AliITSsimulationFastPoints.h"

ClassImp(AliITSFDigitizer)

//______________________________________________________________________
AliITSFDigitizer::AliITSFDigitizer() : AliDigitizer(),
fITS(0),
fInit(kFALSE){
//
// Default constructor.
//
}
//______________________________________________________________________
AliITSFDigitizer::AliITSFDigitizer(AliRunDigitizer *mngr) : AliDigitizer(mngr),
fITS(0),
fInit(kFALSE){
//
// Standard constructor.
//
}
//______________________________________________________________________
AliITSFDigitizer::AliITSFDigitizer(const AliITSFDigitizer &rec):AliDigitizer(rec),
fITS(rec.fITS),
fInit(rec.fInit){
    // Copy constructor. 
  
}
//______________________________________________________________________
AliITSFDigitizer& AliITSFDigitizer::operator=(const AliITSFDigitizer& /*source*/){

    // Assignment operator. This is a function which is not allowed to be
    // done.
    Error("operator=","Assignment operator not allowed\n");
    return *this; 
}

//______________________________________________________________________
AliITSFDigitizer::~AliITSFDigitizer(){
//
// Default destructor. 
//
    fITS = 0; // don't delete fITS. Done else where.
}
//______________________________________________________________________
Bool_t AliITSFDigitizer::Init(){
//
// Initialization. 
// loads ITS and ITSgeom.
// Inputs:
//      none.
// Outputs:
//      none.

  
  fInit = kFALSE;
  if(!gAlice) {
    fITS      = 0;
    Warning("Init","gAlice not found");
    return fInit;
  }
  fITS = (AliITS *)(gAlice->GetDetector("ITS"));
  if(!fITS){
    Warning("Init","ITS not found");
    return fInit;
  } 
  if(!fITS->GetITSgeom()){
    Warning("Init","ITS geometry not found");
    return fInit;
  }
  return fInit = kTRUE;
}
////////////////////////////////////////////////////////////////////////
void AliITSFDigitizer::Exec(Option_t* opt){
//
// Main digitization function. 
// Inputs:
//      Option_t * opt  "deb" ... more verbose output 
//

    AliITSsimulationFastPoints *sim = new AliITSsimulationFastPoints();
    AliRunLoader* outrl = AliRunLoader::GetRunLoader(
                                        fManager->GetOutputFolderName());
    if (outrl == 0x0){
	Error("Exec","Can not find Run Loader in output folder.");
	return;
    }
    AliLoader* outgime = outrl->GetLoader("ITSLoader");
    if (outgime == 0x0){
	Error("Exec","Can not get TOF Loader from Output Run Loader.");
	return;
    }
    if(strstr(opt,"deb")){
	Info("Exec","sim=%p, outrl=%p, outgime=%p",sim,outrl,outgime);
    }
    TTree* outputTreeR = outgime->TreeR();
    if (outputTreeR == 0x0){
	outgime->MakeTree("R");
	outputTreeR = outgime->TreeR();
    }
    
    TClonesArray* recPoints = new TClonesArray("AliITSRecPoint",1000);
    TBranch* branch = outputTreeR->GetBranch("ITSRecPointsF");
    if(branch) branch->SetAddress(recPoints);
    else outputTreeR->Branch("ITSRecPointsF",&recPoints);
    Int_t nModules;
    fITS->InitModules(-1,nModules);

// load hits into modules
    for (Int_t iFile = 0; iFile < fManager->GetNinputs(); iFile++){
	AliRunLoader* rl = AliRunLoader::GetRunLoader(
                                       fManager->GetInputFolderName(iFile));
	if (rl == 0x0){
	    Error("Exec","Can not find Run Loader in input %d folder.",iFile);
	    return;
	}

	AliLoader* gime = rl->GetLoader("ITSLoader");
	if (gime == 0x0){
	    Error("Exec","Can not get TOF Loader from Input %d Run Loader.",
		  iFile);
	    return;
	}

	gime->LoadHits();
	fITS->FillModules(gime->TreeH(),fManager->GetMask(iFile));
	gime->UnloadHits();
    }
  
// transform hits to fast rec points

    AliITSgeom *geom = fITS->GetITSgeom();
    for(Int_t moduleIndex = 0; moduleIndex<geom->GetIndexMax(); moduleIndex++){
	sim->CreateFastRecPoints(moduleIndex,recPoints);
	outputTreeR->Fill();
	recPoints->Clear();
    }
    outrl->WriteRecPoints("OVERWRITE");
//  outputTreeR->AutoSave();
}
////////////////////////////////////////////////////////////////////////
