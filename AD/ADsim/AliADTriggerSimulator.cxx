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
// Class AliADTriggerSimulator
// ------------------------------
//  Simulate the AD Trigger response
// Use FEE parameters stored in Database
// Can work on real data or in simulation
//

#include <TTree.h>
#include <TClonesArray.h>
#include <TParameter.h>

#include "AliLog.h"
#include "AliCDBManager.h"
#include "AliCDBEntry.h"
#include "AliCDBStorage.h"
#include "AliCDBId.h"
#include "AliADCalibData.h"
#include "AliADLogicalSignal.h"
#include "AliADTriggerSimulator.h"
#include "AliADdigit.h"
#include "AliADConst.h"
#include "AliCTPTimeParams.h"

ClassImp(AliADTriggerSimulator)

//_____________________________________________________________________________
AliADTriggerSimulator::AliADTriggerSimulator(TTree * digitsTree, TClonesArray* digits) : 
TObject(),fCalibData(NULL),fDigitsTree(digitsTree),fDigits(digits),fTriggerWord(0)
{
	// constructor
	fCalibData = LoadCalibData();
	LoadClockOffset();
	
	for(int i=0;i<16;i++) {
		fBBFlags[i] = fBGFlags[i] = kFALSE;
		fCharges[i] = 0.;
	}
	GenerateBBWindows();
	GenerateBGWindows();
	for (int i=0; i<kNCIUBoards; i++) {
		fBBLatch[i] = new AliADLogicalSignal(fCalibData->GetLatchWin1(i),0); 
		fBGLatch[i] = new AliADLogicalSignal(fCalibData->GetLatchWin2(i),0); 
		fBBReset[i] = new AliADLogicalSignal(fCalibData->GetResetWin1(i),0);
		fBGReset[i] = new AliADLogicalSignal(fCalibData->GetResetWin2(i),0);		
	}
}
//_____________________________________________________________________________
AliADTriggerSimulator::AliADTriggerSimulator() : 
TObject(),fCalibData(NULL),fDigitsTree(NULL),fDigits(NULL),fTriggerWord(0)
{
	// Default constructor
	fCalibData = LoadCalibData();
	LoadClockOffset();

	for(int i=0;i<16;i++) {
		fBBFlags[i] = fBGFlags[i] = kFALSE;
		fCharges[i] = 0;
	}
	GenerateBBWindows();
	GenerateBGWindows();
	for (int i=0; i<kNCIUBoards; i++) {
		fBBLatch[i] = new AliADLogicalSignal(fCalibData->GetLatchWin1(i),0); 
		fBGLatch[i] = new AliADLogicalSignal(fCalibData->GetLatchWin2(i),0); 
		fBBReset[i] = new AliADLogicalSignal(fCalibData->GetResetWin1(i),0);
		fBGReset[i] = new AliADLogicalSignal(fCalibData->GetResetWin2(i),0);		
	}
}

//_____________________________________________________________________________
AliADTriggerSimulator::~AliADTriggerSimulator(){
// Destructor
  for (Int_t i=0; i<kNCIUBoards; i++) {
    delete fBBGate[i];
    delete fBGGate[i];
    delete fBBLatch[i];
    delete fBBReset[i];
    delete fBGLatch[i];
    delete fBGReset[i];
  }
}

//_____________________________________________________________________________
void AliADTriggerSimulator::GenerateBBWindows() 
{
  // Generates the BB observation window
  // In case gates are open the windows are equal to 25ns
  if (AreGatesOpen()) {
	for (int i=0; i<kNCIUBoards; i++) {
	        fBBGate[i] = new AliADLogicalSignal();
		fBBGate[i]->SetStartTime(0.);
		fBBGate[i]->SetStopTime(25.0);
	}    
  }
  else {
	for (int i=0; i<kNCIUBoards; i++) {
		AliADLogicalSignal clk1BB(fCalibData->GetClk1Win1(i),fCalibData->GetDelayClk1Win1(i));
		AliADLogicalSignal clk2BB(fCalibData->GetClk2Win1(i),fCalibData->GetDelayClk2Win1(i));
		fBBGate[i] = new AliADLogicalSignal(clk1BB & clk2BB);
	}
  }
}
//_____________________________________________________________________________
void AliADTriggerSimulator::GenerateBGWindows() 
{
  // Generates the BG observation window
  // In case gates are open the windows are equal to 25ns
  if (AreGatesOpen()) {
	for (int i=0; i<kNCIUBoards; i++) {
	        fBGGate[i] = new AliADLogicalSignal();
		fBGGate[i]->SetStartTime(0.);
		fBGGate[i]->SetStopTime(25.0);
	}    
  }
  else {
	for (int i=0; i<kNCIUBoards; i++) {
		AliADLogicalSignal clk1BG(fCalibData->GetClk1Win2(i),fCalibData->GetDelayClk1Win2(i));
		AliADLogicalSignal clk2BG(fCalibData->GetClk2Win2(i),fCalibData->GetDelayClk2Win2(i));
		fBGGate[i] = new AliADLogicalSignal(clk1BG & clk2BG);
		// In AD-A we have a shift by -25ns which is controlled by
		// 'Delay Win2' = 7 instead of default 6.
		// The flag is not stored in OCDB so we have manually shift the
		// trigger windows
		if (i < 4) {
		  fBGGate[i]->SetStartTime(fBGGate[i]->GetStartTime()-25.0);
		  fBGGate[i]->SetStopTime(fBGGate[i]->GetStopTime()-25.0);
		}
	}
  }
}

//_____________________________________________________________________________
AliADCalibData * AliADTriggerSimulator::LoadCalibData() const 
{
	// Gets Trigger object for AD set
        AliDebug(1,"Loading Trigger parameters");
	AliCDBManager *man = AliCDBManager::Instance();
	
	
	AliCDBEntry *entry=0;
	
	entry = man->Get("AD/Calib/Data");
	if(!entry){
		AliFatal("Load of calibration data from default storage failed!");
		return NULL;
	}
	
	AliADCalibData *calibData = NULL;
	
	if (entry) calibData = (AliADCalibData*) entry->GetObject();
	if (!calibData)  AliError("No Trigger data from database !");
	
	return calibData;
}


//_____________________________________________________________________________
void AliADTriggerSimulator::LoadClockOffset()
{
  // This method is used in order to
  // retrieve the TDC clock offset including
  // roll-over, trig count and CTP L0->L1 delay

  AliCDBEntry *entry0 = AliCDBManager::Instance()->Get("AD/Calib/Data");
  if (!entry0) {
    AliFatal("AD Calib object is not found in OCDB !");
    return;
  }
  AliADCalibData *calibdata = (AliADCalibData*) entry0->GetObject();

  AliCDBEntry *entry = AliCDBManager::Instance()->Get("GRP/CTP/CTPtiming");
  if (!entry) {
    AliFatal("CTP timing parameters are not found in OCDB !");
    return;
  }
  AliCTPTimeParams *ctpParams = (AliCTPTimeParams*)entry->GetObject();
  Float_t l1Delay = (Float_t)ctpParams->GetDelayL1L0()*25.0;

  AliCDBEntry *entry1 = AliCDBManager::Instance()->Get("GRP/CTP/TimeAlign");
  if (!entry1) {
    AliFatal("CTP time-alignment is not found in OCDB !");
    return;
  }
  AliCTPTimeParams *ctpTimeAlign = (AliCTPTimeParams*)entry1->GetObject();
  l1Delay += ((Float_t)ctpTimeAlign->GetDelayL1L0()*25.0);

  for(Int_t board = 0; board < kNCIUBoards; ++board) {
    fClockOffset[board] = (((Float_t)calibdata->GetRollOver(board)-
			    (Float_t)calibdata->GetTriggerCountOffset(board))*25.0-
			   l1Delay+
			   kADOffset);
    AliDebug(1,Form("Board %d Offset %f",board,fClockOffset[board]));
  }
}

//_____________________________________________________________________________
void AliADTriggerSimulator::Run() {
	//AliInfo("Generating AD Triggers");
	
	// Loop over AD entries
	Int_t nEntries = (Int_t)fDigitsTree->GetEntries();
	for (Int_t ievt=0; ievt<nEntries; ievt++) {
		fDigitsTree->GetEvent(ievt);
		
		Int_t nDigits = fDigits->GetEntriesFast();
		
		for (Int_t iDigit=0; iDigit<nDigits; iDigit++) {
			AliADdigit* digit = (AliADdigit*)fDigits->At(iDigit);
			
			Int_t integrator = digit->Integrator();
			Int_t pmNumber   = digit->PMNumber();
			Int_t board   = AliADCalibData::GetBoardNumber(pmNumber);
			if (board < 0) continue;
			
			if(fCalibData->GetEnableCharge(pmNumber)) {
				fCharges[pmNumber] = digit->ChargeADC(kNClocks/2);
				if(fCalibData->GetPedestalSubtraction(board)) {
					if(fCharges[pmNumber]>=(Float_t) fCalibData->GetOnlinePedestalCut(integrator,pmNumber)){ 
						fCharges[pmNumber] -= (Float_t) fCalibData->GetOnlinePedestal(integrator,pmNumber);
					} else {
						fCharges[pmNumber] = 0.;
					}
				}
			} else {
				fCharges[pmNumber] = 0.;
			}
			
			Float_t time = digit->Time();
			time -= fClockOffset[board];

			AliDebug(10,Form(" Digit: %f %d %d %d %d %d %d %d %d",digit->Time(),
					 digit->ChargeADC(8),digit->ChargeADC(9),digit->ChargeADC(10),
					 digit->ChargeADC(11),digit->ChargeADC(12),digit->ChargeADC(13),
					 digit->ChargeADC(14),digit->ChargeADC(15)));
			AliDebug(10,Form(" PM nb : %d ; TDC= %f(%f)  Enable Time %d charge %d inCoin %d charge %f",
					 pmNumber,time,digit->Time(),
					 fCalibData->GetEnableTiming(pmNumber),fCalibData->GetEnableCharge(pmNumber),
					 fBBGate[board]->IsInCoincidence(time),fCharges[pmNumber]));
			fBBFlags[pmNumber] = fCalibData->GetEnableTiming(pmNumber) && fBBGate[board]->IsInCoincidence(time);
			fBGFlags[pmNumber] = fCalibData->GetEnableTiming(pmNumber) && fBGGate[board]->IsInCoincidence(time);
			// Store the BB and BG flags in the digits tree
			digit->SetBBflag(fBBFlags[pmNumber]);
			digit->SetBGflag(fBGFlags[pmNumber]);
			
		} // end of loop over digits
	} // end of loop over events in digits tree
	
	Int_t nBBflagsADA = 0;
	Int_t nBBflagsADC = 0;
	Int_t nBGflagsADA = 0;
	Int_t nBGflagsADC = 0;
	Float_t chargeADA   = 0.;
	Float_t chargeADC   = 0.;

	for(int i=0;i<16;i++) {
		if(i<8) {
			nBBflagsADC += fBBFlags[i]; 
			nBGflagsADC += fBGFlags[i];
			chargeADC += fCharges[i];
		} else {
			nBBflagsADA += fBBFlags[i]; 
			nBGflagsADA += fBGFlags[i];
			chargeADA += fCharges[i];
		}
		//AliInfo(Form("Ch %d BB=%d BG=%d",i,fBBFlags[i],fBGFlags[i] )); 
	}
	
	// BBA
	if(nBBflagsADA>=fCalibData->GetBBAThreshold())  SetBBA();
	
	// BBC
	if(nBBflagsADC>=fCalibData->GetBBCThreshold())  SetBBC();

	// BBA_AND_BBC
	if(GetBBA() && GetBBC())  SetBBAandBBC();
	
	// BBA_OR_BBC
	if(GetBBA() || GetBBC()) SetBBAorBBC();

	// BGA
	if(nBGflagsADA>=fCalibData->GetBGAThreshold()) SetBGA();

	// BGC
	if(nBGflagsADC>=fCalibData->GetBGCThreshold()) SetBGC();
	
	// BGA_AND_BBC (Beam Gas from RB24 side)
	if(nBBflagsADC>=fCalibData->GetBBCForBGThreshold() && GetBGA()) SetBGAandBBC();
	
	// BGC_AND_BBA (Beam Gas from RB26 side)
	if(nBBflagsADA>=fCalibData->GetBBAForBGThreshold() && GetBGC()) SetBGCandBBA();

	
	// MTA_AND_MTC (Multiplicity Trigger)
	if((nBBflagsADA<=fCalibData->GetMultADAThrHigh() && nBBflagsADA>=fCalibData->GetMultADAThrLow())
	   && (nBBflagsADC<=fCalibData->GetMultADCThrHigh() && nBBflagsADC>=fCalibData->GetMultADCThrLow()) ) 
		SetMTAandMTC();
	
	// MTA_OR_MTC (Multiplicity Trigger)
	if((nBBflagsADA<=fCalibData->GetMultADAThrHigh() && nBBflagsADA>=fCalibData->GetMultADAThrLow())
	   || (nBBflagsADC<=fCalibData->GetMultADCThrHigh() && nBBflagsADC>=fCalibData->GetMultADCThrLow()) ) 
		SetMTAorMTC();
	
	// BGA_OR_BGC
	if(GetBGA() || GetBGC()) SetBGAorBGC();
	
	// (BGA and BBC) or (BGC and BBA) (Beam Gas from one of the two sides)
	if(GetBGAandBBC() || GetBGCandBBA()) SetBeamGas();

//	AliInfo(Form("BB Flags : ADA = %d  ADC = %d ",nBBflagsADA, nBBflagsADC )); 
//	AliInfo(Form("BG Flags : ADA = %d  ADC = %d ",nBGflagsADA, nBGflagsADC )); 
//	AliInfo(Form("Charges  : ADA = %d  ADC = %d ",chargeADA, chargeADC )); 
	
}

//_____________________________________________________________________________
Bool_t AliADTriggerSimulator::AreGatesOpen() const {
  // The method check if the gates are suppossed to be open
  // (corresponding to 'Test Window' flag in DCS).
  // Since the flag is not stored in OCDB, we just check if
  // all the clock delays are 0 or not.
  // This rules should be followed when setting up the detector
  // at the level of DCS

  for (int i=0; i<kNCIUBoards; i++) {
    if (fCalibData->GetDelayClk1Win1(i)!=0 ||
	fCalibData->GetDelayClk2Win1(i)!=0 ||
	fCalibData->GetDelayClk1Win2(i)!=0 ||
	fCalibData->GetDelayClk2Win2(i)!=0)
      return kFALSE;
  }
  return kTRUE;
}

//_____________________________________________________________________________
void AliADTriggerSimulator::Print(Option_t* /* opt */) const
{
  // Prints the trigger windows as
  // initialized from the OCDB
  for (int i=0; i<kNCIUBoards; i++) {
    std::cout << "Board=" << i << "   BB (" << fBBGate[i]->GetStartTime() << " -> " << fBBGate[i]->GetStopTime() << ")   BG (" << fBGGate[i]->GetStartTime() << " -> " << fBGGate[i]->GetStopTime() << ")" << std::endl;
  }
  std::cout << std::endl;
}


