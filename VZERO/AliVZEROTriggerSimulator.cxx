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
// Class AliVZEROTriggerSimulator
// ------------------------------
//  Simulate the VZERO Trigger response
// Use FEE parameters stored in Database
// Can work on real data or in simulation
//

#include <TTree.h>
#include <TClonesArray.h>

#include "AliLog.h"
#include "AliCDBManager.h"
#include "AliCDBEntry.h"
#include "AliCDBStorage.h"
#include "AliCDBId.h"
#include "AliVZEROTriggerData.h"
#include "AliVZEROLogicalSignal.h"
#include "AliVZEROTriggerSimulator.h"
#include "AliVZEROdigit.h"

ClassImp(AliVZEROTriggerSimulator)

//_____________________________________________________________________________
AliVZEROTriggerSimulator::AliVZEROTriggerSimulator(TTree * digitsTree, TClonesArray* digits) : 
TObject(),fTriggerData(NULL),fDigitsTree(digitsTree),fDigits(digits),fTriggerWord(0)
{
	// constructor
	fTriggerData = LoadTriggerData();
	
	for(int i=0;i<64;i++) {
		fBBFlags[i] = fBGFlags[i] = kFALSE;
		fCharges[i] = 0;
	}
	GenerateBBWindows();
	GenerateBGWindows();
	for (int i=0; i<AliVZEROTriggerData::kNCIUBoards; i++) {
		fBBLatch[i] = new AliVZEROLogicalSignal(fTriggerData->GetLatchWin1(i),0); 
		fBGLatch[i] = new AliVZEROLogicalSignal(fTriggerData->GetLatchWin2(i),0); 
		fBBReset[i] = new AliVZEROLogicalSignal(fTriggerData->GetResetWin1(i),0);
		fBGReset[i] = new AliVZEROLogicalSignal(fTriggerData->GetResetWin2(i),0);		
	}
}
//_____________________________________________________________________________
AliVZEROTriggerSimulator::AliVZEROTriggerSimulator() : 
TObject(),fTriggerData(NULL),fDigitsTree(NULL),fDigits(NULL),fTriggerWord(0)
{
	// Default constructor
	fTriggerData = LoadTriggerData();

	for(int i=0;i<64;i++) {
		fBBFlags[i] = fBGFlags[i] = kFALSE;
		fCharges[i] = 0;
	}
	GenerateBBWindows();
	GenerateBGWindows();
	for (int i=0; i<AliVZEROTriggerData::kNCIUBoards; i++) {
		fBBLatch[i] = new AliVZEROLogicalSignal(fTriggerData->GetLatchWin1(i),0); 
		fBGLatch[i] = new AliVZEROLogicalSignal(fTriggerData->GetLatchWin2(i),0); 
		fBBReset[i] = new AliVZEROLogicalSignal(fTriggerData->GetResetWin1(i),0);
		fBGReset[i] = new AliVZEROLogicalSignal(fTriggerData->GetResetWin2(i),0);		
	}
}

//_____________________________________________________________________________
AliVZEROTriggerSimulator::~AliVZEROTriggerSimulator(){
// Destructor
	if(fBBGate) delete [] fBBGate;
	if(fBGGate) delete [] fBGGate;
	if(fBBLatch) delete [] fBBLatch;
	if(fBBReset) delete [] fBBReset;
	if(fBGLatch) delete [] fBGLatch;
	if(fBGReset) delete [] fBGReset;
}

//_____________________________________________________________________________
void AliVZEROTriggerSimulator::GenerateBBWindows() 
{
	// Generates the BB observation window
	for (int i=0; i<AliVZEROTriggerData::kNCIUBoards; i++) {
		AliVZEROLogicalSignal clk1BB(fTriggerData->GetClk1Win1(i),fTriggerData->GetDelayClk1Win1(i));
		AliVZEROLogicalSignal clk2BB(fTriggerData->GetClk2Win1(i),fTriggerData->GetDelayClk2Win1(i));
		fBBGate[i] = new AliVZEROLogicalSignal(clk1BB & clk2BB);
	}
}
//_____________________________________________________________________________
void AliVZEROTriggerSimulator::GenerateBGWindows() 
{
	// Generates the BG observation window
	for (int i=0; i<AliVZEROTriggerData::kNCIUBoards; i++) {
		AliVZEROLogicalSignal clk1BG(fTriggerData->GetClk1Win2(i),fTriggerData->GetDelayClk1Win2(i));
		AliVZEROLogicalSignal clk2BG(fTriggerData->GetClk2Win2(i),fTriggerData->GetDelayClk2Win2(i));
		fBGGate[i] = new AliVZEROLogicalSignal(clk1BG & clk2BG);
	}
}

//_____________________________________________________________________________
AliVZEROTriggerData * AliVZEROTriggerSimulator::LoadTriggerData() const 
{
	// Gets Trigger object for VZERO set
        AliDebug(1,"Loading Trigger parameters");
	AliCDBManager *man = AliCDBManager::Instance();
	
	
	AliCDBEntry *entry=0;
	
	entry = man->Get("VZERO/Trigger/Data");
	if(!entry){
		AliWarning("Load of calibration data from default storage failed!");
		AliWarning("Calibration data will be loaded from local storage ($ALICE_ROOT)");
		
		man->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
		entry = man->Get("VZERO/Trigger/Data",0);
	}
	
	// Retrieval of data in directory VZERO/Calib/Trigger:
	
	AliVZEROTriggerData *triggerData = NULL;
	
	if (entry) triggerData = (AliVZEROTriggerData*) entry->GetObject();
	if (!triggerData)  AliError("No Trigger data from database !");
	
	return triggerData;
}

//_____________________________________________________________________________
void AliVZEROTriggerSimulator::Run() {
	//AliInfo("Generating VZERO Triggers");
	
	// Loop over VZERO entries
	Int_t nEntries = (Int_t)fDigitsTree->GetEntries();
	for (Int_t ievt=0; ievt<nEntries; ievt++) {
		fDigitsTree->GetEvent(ievt);
		
		Int_t nDigits = fDigits->GetEntriesFast();
		
		for (Int_t iDigit=0; iDigit<nDigits; iDigit++) {
			AliVZEROdigit* digit = (AliVZEROdigit*)fDigits->At(iDigit);
			
			Int_t integrator = digit->Integrator();
			Int_t pmNumber   = digit->PMNumber();
			Int_t board   = pmNumber / 8;
			Int_t channel = pmNumber % 8;
			
			if(fTriggerData->GetEnableCharge(board,channel)) {
				fCharges[pmNumber] = digit->ADC();
				if(fTriggerData->GetPedestalSubtraction(board)) {
					if(fCharges[pmNumber]>=fTriggerData->GetPedestalCut(integrator,board,channel)){ 
						fCharges[pmNumber] -= fTriggerData->GetPedestal(integrator,board,channel);
					} else {
						fCharges[pmNumber] = 0;
					}
				}
			} else {
				fCharges[pmNumber] = 0;
			}
			
			Float_t time = (Float_t)digit->Time() / 10.; // digit->Time() in bin of 100 picoseconds. Divide by 10 to have it in nanoseconds
			time += fTriggerData->GetDelayHit(board,channel);
			
			//AliInfo(Form(" PM nb : %d ; ADC= %d ; TDC= %f  Enable Time %d charge %d",pmNumber,digit->ADC(),time,fTriggerData->GetEnableTiming(board,channel),fTriggerData->GetEnableCharge(board,channel)));
			fBBFlags[pmNumber] = fTriggerData->GetEnableTiming(board,channel) && fBBGate[board]->IsInCoincidence(time);
			fBGFlags[pmNumber] = fTriggerData->GetEnableTiming(board,channel) && fBGGate[board]->IsInCoincidence(time);
			
		} // end of loop over digits
	} // end of loop over events in digits tree
	
	Int_t nBBflagsV0A = 0;
	Int_t nBBflagsV0C = 0;
	Int_t nBGflagsV0A = 0;
	Int_t nBGflagsV0C = 0;
	Int_t chargeV0A   = 0;
	Int_t chargeV0C   = 0;

	for(int i=0;i<64;i++) {
		if(i<32) {
			nBBflagsV0C += fBBFlags[i]; 
			nBGflagsV0C += fBGFlags[i];
			chargeV0C += fCharges[i];
		} else {
			nBBflagsV0A += fBBFlags[i]; 
			nBGflagsV0A += fBGFlags[i];
			chargeV0A += fCharges[i];
		}
		//AliInfo(Form("Ch %d BB=%d BG=%d",i,fBBFlags[i],fBGFlags[i] )); 
	}
	
	// BBA
	if(nBBflagsV0A>=fTriggerData->GetBBAThreshold())  SetBBA();
	
	// BBC
	if(nBBflagsV0C>=fTriggerData->GetBBCThreshold())  SetBBC();

	// BBA_AND_BBC
	if(GetBBA() && GetBBC())  SetBBAandBBC();
	
	// BBA_OR_BBC
	if(GetBBA() || GetBBC()) SetBBAorBBC();

	// BGA
	if(nBGflagsV0A>=fTriggerData->GetBGAThreshold()) SetBGA();

	// BGC
	if(nBGflagsV0C>=fTriggerData->GetBGCThreshold()) SetBGC();
	
	// BGA_AND_BBC (Beam Gas from RB24 side)
	if(nBBflagsV0C>=fTriggerData->GetBBCForBGThreshold() && GetBGA()) SetBGAandBBC();
	
	// BGC_AND_BBA (Beam Gas from RB26 side)
	if(nBBflagsV0A>=fTriggerData->GetBBAForBGThreshold() && GetBGC()) SetBGCandBBA();

	// CTA1_AND_CTC1 (Centrality trigger 1)
	if(chargeV0A>=fTriggerData->GetCentralityV0AThrLow() && chargeV0C>=fTriggerData->GetCentralityV0CThrLow()) SetCTA1andCTC1();

	// CTA1_OR_CTC1 (Centrality trigger 1)
	if(chargeV0A>=fTriggerData->GetCentralityV0AThrLow() || chargeV0C>=fTriggerData->GetCentralityV0CThrLow()) SetCTA1orCTC1();
	
	// CTA2_AND_CTC2 (Centrality trigger 2)
	if(chargeV0A>=fTriggerData->GetCentralityV0AThrHigh() && chargeV0C>=fTriggerData->GetCentralityV0CThrHigh()) SetCTA2andCTC2();
	
	// CTA2_OR_CTC2 (Centrality trigger 2)
	if(chargeV0A>=fTriggerData->GetCentralityV0AThrHigh() || chargeV0C>=fTriggerData->GetCentralityV0CThrHigh()) SetCTA2orCTC2();
	
	// MTA_AND_MTC (Multiplicity Trigger)
	if((nBBflagsV0A<=fTriggerData->GetMultV0AThrHigh() && nBBflagsV0A>=fTriggerData->GetMultV0AThrLow())
	   && (nBBflagsV0C<=fTriggerData->GetMultV0CThrHigh() && nBBflagsV0C>=fTriggerData->GetMultV0CThrLow()) ) 
		SetMTAandMTC();
	
	// MTA_OR_MTC (Multiplicity Trigger)
	if((nBBflagsV0A<=fTriggerData->GetMultV0AThrHigh() && nBBflagsV0A>=fTriggerData->GetMultV0AThrLow())
	   || (nBBflagsV0C<=fTriggerData->GetMultV0CThrHigh() && nBBflagsV0C>=fTriggerData->GetMultV0CThrLow()) ) 
		SetMTAorMTC();
	
	// BGA_OR_BGC
	if(GetBGA() || GetBGC()) SetBGAorBGC();
	
	// (BGA and BBC) or (BGC and BBA) (Beam Gas from one of the two sides)
	if(GetBGAandBBC() || GetBGCandBBA()) SetBeamGas();

//	AliInfo(Form("BB Flags : V0A = %d  V0C = %d ",nBBflagsV0A, nBBflagsV0C )); 
//	AliInfo(Form("BG Flags : V0A = %d  V0C = %d ",nBGflagsV0A, nBGflagsV0C )); 
//	AliInfo(Form("Charges  : V0A = %d  V0C = %d ",chargeV0A, chargeV0C )); 
	
}



