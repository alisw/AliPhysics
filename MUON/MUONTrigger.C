


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

// Macro MUONTrigger.C (TO BE COMPILED)
// for testing the C++ trigger code

#if !defined(__CINT__) || defined(__MAKECINT__)
#include "AliRun.h"
#include "AliMUON.h"
#include "AliMUONLoader.h"
#include "AliMUONData.h"
#include "AliMUONDigit.h"
#include "AliMUONTriggerDecision.h"
#endif
void MUONTrigger(char * FileName="galice.root", Int_t nevents=1, Int_t idebug=1)
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
//  AliMUONLoader *MUONLoader = (AliMUONLoader*) RunLoader->GetLoader("MUONLoader");
    AliLoader * MUONLoader = RunLoader->GetLoader("MUONLoader");
    MUONLoader->LoadDigits("READ");
    // Creating MUON data container
    AliMUONData* MUONData = new AliMUONData(MUONLoader,"MUON","MUON");
    
    // Creating MUONTriggerDecision
    AliMUONTriggerDecision* decision = new AliMUONTriggerDecision(MUONLoader ,idebug,MUONData);
    
    AliMUONDigit * mDigit;
    
    Int_t tracks[10];
    Int_t charges[10];
    Int_t digits[7];
    
    for(Int_t ievent = 0; ievent < nevents; ievent++) {
        printf(">>> Event %i out of %i \n",ievent,nevents);
        RunLoader->GetEvent(ievent);
        MUONData->SetTreeAddress("D");
        
	MUONData->GetDigits();
	for(Int_t ichamber=10; ichamber<14; ichamber++) {         
	    Int_t idigit, ndigits;
	    ndigits = (Int_t) MUONData->Digits(ichamber)->GetEntriesFast();
//            printf(">>> Chamber Cathode ndigits %d %d %d\n",ichamber,icathode,ndigits);             
	    for(idigit=0; idigit<ndigits; idigit++) {
		mDigit = static_cast<AliMUONDigit*>(MUONData->Digits(ichamber)->At(idigit));
		digits[0] = mDigit->PadX();
		digits[1] = mDigit->PadY();
		digits[2] = mDigit->Cathode();
		digits[3] = mDigit->Signal();
		digits[4] = mDigit->Physics();
		digits[5] = mDigit->Hit();
		digits[6] = mDigit->DetElemId();
		
		Int_t digitindex = 0 ;
                printf("ichamber ix iy %d %d %d \n",ichamber,mDigit->PadX(),mDigit->PadY());
		
		decision->AddDigit(ichamber, tracks, charges, digits, digitindex );
	    } // loop on digits
	} // loop on chambers
	//  } // loop on cathodes      
        MUONData->ResetDigits();
        decision->Trigger();
        decision->ClearDigits();
    } // loop on events
    MUONLoader->UnloadDigits();
}
