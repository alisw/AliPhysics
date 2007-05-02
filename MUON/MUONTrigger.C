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

// This macro is to be used to check the trigger algorithm w/o having to
// (re-)perform simulation and digitalization. 
// see full description in the REDAME file
// Author: P.Crochet (LPC)

#if !defined(__CINT__) || defined(__MAKECINT__)
#include "AliRun.h"
#include "AliMUON.h"
#include "AliMUONLoader.h"
#include "AliMUONData.h"
#include "AliMUONDigit.h"
#include "AliMUONTriggerElectronics.h"
#include "AliMUONCalibrationData.h"
#include "AliCDBManager.h"
#include <TClonesArray.h>
#endif
void MUONTrigger(char * FileNameSim="galice_sim.root", char * FileName="galice.root")
{
    // Creating Run Loader and openning file containing Digits 
    AliRunLoader * RunLoader = AliRunLoader::Open(FileName,"MUONLoader","UPDATE");
    if (RunLoader ==0x0) {
        printf(">>> Error : Error Opening %s file \n",FileName);
        return;
    }
    // Loading AliRun master
    if (RunLoader->GetAliRun() == 0x0) RunLoader->LoadgAlice();
    gAlice = RunLoader->GetAliRun();
    
    // Loading MUON subsystem
    AliLoader * MUONLoader = RunLoader->GetLoader("MUONLoader");
    MUONLoader->LoadDigits("READ");
    MUONLoader->LoadRecPoints("UPDATE"); // absolutely essential !!!
    
    Int_t nevents;
    nevents = RunLoader->GetNumberOfEvents();
    
    // Creating MUON data container
    AliMUONData* MUONData = new AliMUONData(MUONLoader,"MUON","MUON");
    
    // Creating MUONTriggerDecision
    TTask *TriggerProcessor;
    AliCDBManager* cdbManager = AliCDBManager::Instance();
    cdbManager->SetDefaultStorage("local://$ALICE_ROOT");
    Int_t runnumber = gAlice->GetRunNumber();
    AliMUONCalibrationData *CalibrationData = new AliMUONCalibrationData(runnumber);
    TriggerProcessor = new AliMUONTriggerElectronics(MUONData,CalibrationData);
    
    // Testing if Trigger has already been done
    RunLoader->GetEvent(0);    
    if (MUONLoader->TreeR()) {
	if (MUONData->IsTriggerBranchesInTree()) {
	    MUONLoader->UnloadRecPoints();
	    MUONLoader->LoadRecPoints("RECREATE");
	    printf("Recreating recpoints files\n");
	}
    }
    
    AliMUONDigit * mDigit;    
    Float_t digits[7];
    
    for(Int_t ievent = 0; ievent < nevents; ievent++) {
	printf(">>> Event %i out of %i \n",ievent,nevents);
	RunLoader->GetEvent(ievent);
	MUONData->SetTreeAddress("D");
	
	MUONData->GetDigits();
	for(Int_t ichamber=10; ichamber<14; ichamber++) {	    
	    Int_t ndigits = (Int_t) MUONData->Digits(ichamber)->GetEntriesFast();
	    for(Int_t idigit=0; idigit<ndigits; idigit++) {
		mDigit = static_cast<AliMUONDigit*>(MUONData->Digits(ichamber)->At(idigit));		 
		
		digits[0] = mDigit->PadX();
		digits[1] = mDigit->PadY();
		digits[2] = mDigit->Cathode();
		digits[3] = mDigit->Signal();
		digits[4] = mDigit->Physics();
		digits[5] = mDigit->Hit();
		digits[6] = mDigit->DetElemId();
		
//		printf("ichamber ix iy %d %d %d \n",ichamber,mDigit->PadX(),mDigit->PadY());
		
	    } // loop on digits
	} // loop on chambers
	
	
	if (MUONLoader->TreeR() == 0x0) {	
	    MUONLoader->MakeRecPointsContainer();
	} else {
	    if (MUONData->IsTriggerBranchesInTree()){ 
		if (ievent==0) MUONLoader->UnloadRecPoints();
		MUONLoader->MakeRecPointsContainer();
		cout << "Recreating RecPointsContainer and deleting previous ones" << "\n";
	    }
	}	
	
	MUONData->MakeBranch("TC");	
	MUONData->SetTreeAddress("TC");
	TriggerProcessor->ExecuteTask();

	MUONData->Fill("TC");
	MUONLoader->WriteRecPoints("OVERWRITE");  
        MUONData->ResetDigits();

    } // loop on events
    MUONLoader->UnloadDigits();
    MUONLoader->UnloadRecPoints();
}

