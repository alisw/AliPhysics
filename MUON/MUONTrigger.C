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

/// \ingroup macros
/// \file MUONTrigger.C
/// \brief This macro is to be used to check the trigger algorithm w/o having to
/// (re-)perform simulation and digitalization. 
///
/// See full description on the \ref README_trigger page.
///
/// \author P.Crochet (LPC)

#if !defined(__CINT__) || defined(__MAKECINT__)
#include "AliRun.h"
#include "AliMUON.h"
#include "AliMUONDigit.h"
#include "AliMUONTriggerElectronics.h"
#include "AliMUONCalibrationData.h"
#include "AliCDBManager.h"
#include "AliMUONDataInterface.h"
#include "AliMUONMCDataInterface.h"
#include "AliMUONVTriggerStore.h"
#include "AliMUONDigitStoreV1.h"
#include <TClonesArray.h>
#include "AliMpCDB.h"
#include <TFile.h>
#endif

void MUONTrigger(const char* filename)
{
    // Creating Run Loader and openning file containing Digits 
    AliRunLoader * RunLoader = AliRunLoader::Open(filename,"MUONLoader","UPDATE");
    if (RunLoader ==0x0) {
        printf(">>> Error : Error Opening %s file \n",filename);
        return;
    }
    // Loading AliRun master
    if (RunLoader->GetAliRun() == 0x0) RunLoader->LoadgAlice();
    gAlice = RunLoader->GetAliRun();
    
    // Loading MUON subsystem
    AliLoader* MUONLoader = RunLoader->GetDetectorLoader("MUON");
    MUONLoader->LoadDigits("READ");
    MUONLoader->LoadRecPoints("UPDATE"); // absolutely essential !!!    
    
    // Creating MUONTriggerDecision
    AliCDBManager* cdbManager = AliCDBManager::Instance();
    cdbManager->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
    
    Int_t runnumber = cdbManager->GetRun();
    AliMpCDB::LoadDDLStore();
    
    AliMUONCalibrationData *CalibrationData = new AliMUONCalibrationData(runnumber);
    AliMUONTriggerElectronics *TriggerProcessor	= new AliMUONTriggerElectronics(CalibrationData);
    
    Int_t nevents = RunLoader->GetNumberOfEvents();
    AliMUONVDigitStore* digitStore=0x0;
    AliMUONVTriggerStore* triggerStore=0x0;
    
    for(Int_t ievent = 0; ievent < nevents; ievent++) {
	printf(">>> Event %i out of %i \n",ievent,nevents);
	RunLoader->GetRunLoader()->GetEvent(ievent);
	
	MUONLoader->LoadRecPoints("update");
	MUONLoader->CleanRecPoints();
	MUONLoader->MakeRecPointsContainer();
	TTree* clustersTree = MUONLoader->TreeR();
	TFile* cfile = clustersTree->GetCurrentFile();
	if ( !cfile ) 
	{
	    cout << " could not find Cluster file " << endl;
	    return;
	}
	
	MUONLoader->LoadDigits("read");
	TTree* digitsTree = MUONLoader->TreeD();
	TFile* dfile = digitsTree->GetCurrentFile();
	if ( !dfile ) 
	{
	    cout << " could not find Digit file " << endl;
	    return;
	}
	
// here start reconstruction	
	if (!digitStore) digitStore = AliMUONVDigitStore::Create(*digitsTree);	
	if (!triggerStore) triggerStore = AliMUONVTriggerStore::Create(*digitsTree);
	// insure we start with empty stores
	if ( digitStore ) 
	{
	    digitStore->Clear(); 
	    Bool_t alone = ( triggerStore ? kFALSE : kTRUE );
	    Bool_t ok = digitStore->Connect(*digitsTree,alone);
	    if (!ok)
	    {
		cerr << "Could not connect digitStore to digitsTree \n";
		return;
	    }
	} else {
	    cerr << "digitStore does not exist " << "\n";
	    return;
	}
	
	digitsTree->GetEvent(0);
	
// process trigger response
	TriggerProcessor->Digits2Trigger(*digitStore,*triggerStore);
	
	//triggerStore->Print();
	
	Bool_t ok(kFALSE);
	if ( triggerStore ) {
	    ok = triggerStore->Connect(*clustersTree,kTRUE);
	    if (!ok)
	    {
		cerr << "Could not create triggerStore branches in TreeR " << "\n";
		return;
	    }
	} else {
	    cerr << "triggerStore does not exist " << "\n";
	    return;
	}

// fill TreeR
	clustersTree->Fill();
	MUONLoader->UnloadDigits();
	MUONLoader->WriteRecPoints("OVERWRITE");
	MUONLoader->UnloadRecPoints();

    }  // loop on events

}



