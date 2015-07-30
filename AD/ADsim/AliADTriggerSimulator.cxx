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
		fCharges[i] = 0;
		fTime[i] = 0;
	}
	GenerateBCMask();
	GenerateBBWindows();
	GenerateBGWindows();
			
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
		fTime[i] = 0;
	}
	GenerateBCMask();
	GenerateBBWindows();
	GenerateBGWindows();	
}

//_____________________________________________________________________________
AliADTriggerSimulator::~AliADTriggerSimulator(){
// Destructor
  for (Int_t i=0; i<kNCIUBoards; i++) {
    delete fBBGate[i];
    delete fBGGate[i];
    delete fBCMask[i];
  }
}

//_____________________________________________________________________________
void AliADTriggerSimulator::GenerateBCMask() 
{
  // Generates the BC mask
  // Gate at central clock 25ns wide
  for (Int_t i=0; i<kNCIUBoards; i++) {
  	  fBCMask[i] = new AliADLogicalSignal();
  	  fBCMask[i]->SetStartTime(fWindowOffset[i]);
  	  fBCMask[i]->SetStopTime(fWindowOffset[i]+25.0);
  }    

}


//_____________________________________________________________________________
void AliADTriggerSimulator::GenerateBBWindows() 
{
  // Generates the BB observation window
  // In case gates are open the windows are equal to 25ns
  if (AreGatesOpen()) {
	for (Int_t i=0; i<kNCIUBoards; i++) {
	        fBBGate[i] = new AliADLogicalSignal();
		fBBGate[i]->SetStartTime(fWindowOffset[i]);
		fBBGate[i]->SetStopTime(fWindowOffset[i]+25.0);
	}    
  }
  else {
	for (Int_t i=0; i<kNCIUBoards; i++) {
		AliADLogicalSignal clk1BB(fCalibData->GetClk1Win1(i),fCalibData->GetDelayClk1Win1(i),fCalibData->GetLatchWin1(i),fCalibData->GetResetWin1(i));
		AliADLogicalSignal clk2BB(fCalibData->GetClk2Win1(i),fCalibData->GetDelayClk2Win1(i),fCalibData->GetLatchWin1(i),fCalibData->GetResetWin1(i));
		fBBGate[i] = new AliADLogicalSignal(clk1BB & clk2BB);
		fBBGate[i]->SetStartTime(fBBGate[i]->GetStartTime()+fWindowOffset[i]);
		fBBGate[i]->SetStopTime(fBBGate[i]->GetStopTime()+fWindowOffset[i]);
					
	}
  }
}
//_____________________________________________________________________________
void AliADTriggerSimulator::GenerateBGWindows() 
{
  // Generates the BG observation window
  // In case gates are open the windows are equal to 25ns
  if (AreGatesOpen()) {
	for (Int_t i=0; i<kNCIUBoards; i++) {
	        fBGGate[i] = new AliADLogicalSignal();
		fBGGate[i]->SetStartTime(fWindowOffset[i]);
		fBGGate[i]->SetStopTime(fWindowOffset[i]+25.0);
	}    
  }
  else {
	for (Int_t i=0; i<kNCIUBoards; i++) {
		AliADLogicalSignal clk1BG(fCalibData->GetClk1Win2(i),fCalibData->GetDelayClk1Win2(i),fCalibData->GetLatchWin2(i),fCalibData->GetResetWin2(i));
		clk1BG.SetStartTime(clk1BG.GetStartTime()+2);
		clk1BG.SetStopTime(clk1BG.GetStopTime()+2);
		AliADLogicalSignal clk2BG(fCalibData->GetClk2Win2(i),fCalibData->GetDelayClk2Win2(i),fCalibData->GetLatchWin2(i),fCalibData->GetResetWin2(i));
		clk2BG.SetStartTime(clk2BG.GetStartTime()-2);
		clk2BG.SetStopTime(clk2BG.GetStopTime()-2);
		fBGGate[i] = new AliADLogicalSignal(clk1BG & clk2BG);
		fBGGate[i]->SetStartTime(fBGGate[i]->GetStartTime()+fWindowOffset[i]);
		fBGGate[i]->SetStopTime(fBGGate[i]->GetStopTime()+fWindowOffset[i]);	
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
  //Start of the central clock in HPTDC time
  for(Int_t board = 0; board < kNCIUBoards; ++board) {
    fWindowOffset[board] = (((Float_t)calibdata->GetRollOver(board)-
			    (Float_t)calibdata->GetTriggerCountOffset(board))*25.0
		             -l1Delay
		             -kADOffset);
    AliDebug(1,Form("Board %d Offset %f",board,fWindowOffset[board]));
  }
}
//_____________________________________________________________________________
void AliADTriggerSimulator::FillFlags(Bool_t *bbFlag, Bool_t *bgFlag, Float_t time[16]){

  for(Int_t i = 0; i<16; i++){
  	Int_t board   = AliADCalibData::GetBoardNumber(i);
	
	Bool_t inBCmask = fBCMask[board]->IsInCoincidence(time[i]);
	Bool_t inBBwindow = fBBGate[board]->IsInCoincidence(time[i]);
	Bool_t inBGwindow = fBGGate[board]->IsInCoincidence(time[i]);
	
	/*/
	Bool_t inBBwindow = kFALSE;
	Bool_t inBGwindow = kFALSE;
	
	AliADLogicalSignal fBBGateShifted;
	AliADLogicalSignal fBGGateShifted;
	for(Int_t j=-10; j<=10; j++){
		fBBGateShifted.SetStartTime(fBBGate[board]->GetStartTime() + 25*j);
		fBBGateShifted.SetStopTime(fBBGate[board]->GetStopTime() + 25*j);
		fBGGateShifted.SetStartTime(fBGGate[board]->GetStartTime() + 25*j);
		fBGGateShifted.SetStopTime(fBGGate[board]->GetStopTime() + 25*j);
		
		if(fBBGateShifted.IsInCoincidence(time[i])) inBBwindow = kTRUE;
		if(fBGGateShifted.IsInCoincidence(time[i])) inBGwindow = kTRUE;
		}
	/*/
	bbFlag[i] = inBBwindow;
	bgFlag[i] = inBGwindow;
	
	//AliInfo(Form("Ch %d Time=%.1f BCM=%d BB=%d BG=%d",i,time[i],inBCmask,bbFlag[i],bgFlag[i] ));
  	}
}
//_____________________________________________________________________________
void AliADTriggerSimulator::Run() {
//	AliInfo("Generating AD Triggers");
//	Print("");
	
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
				fCharges[pmNumber] = digit->ChargeADC(kADNClocks/2);
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
			
			fTime[pmNumber] = digit->Time();
			
		} // end of loop over digits
	} // end of loop over events in digits tree
	FillFlags(fBBFlags,fBGFlags,fTime);
	
	Int_t nBBflagsADA = 0;
	Int_t nBBflagsADC = 0;
	Int_t nBGflagsADA = 0;
	Int_t nBGflagsADC = 0;
	Float_t chargeADA   = 0.;
	Float_t chargeADC   = 0.;

	for(int i=0;i<16;i++) {
		if(i<8) chargeADC += fCharges[i];
		else chargeADA += fCharges[i];	
	}
	
	for(Int_t iChannel=0; iChannel<4; iChannel++) {//Loop over pairs of pads
    	    //Enable time is used to turn off the coincidence 
      	    if((!fCalibData->GetEnableTiming(iChannel) || fBBFlags[iChannel]) && (!fCalibData->GetEnableTiming(iChannel+4) || fBBFlags[iChannel+4])) nBBflagsADC++;
      	    if((!fCalibData->GetEnableTiming(iChannel) || fBGFlags[iChannel]) && (!fCalibData->GetEnableTiming(iChannel+4) || fBGFlags[iChannel+4])) nBGflagsADC++;

      	    if((!fCalibData->GetEnableTiming(iChannel+8) || fBBFlags[iChannel+8]) && (!fCalibData->GetEnableTiming(iChannel+12) || fBBFlags[iChannel+12])) nBBflagsADA++;
      	    if((!fCalibData->GetEnableTiming(iChannel+8) || fBGFlags[iChannel+8]) && (!fCalibData->GetEnableTiming(iChannel+12) || fBGFlags[iChannel+12])) nBGflagsADA++;
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

	//AliInfo(Form("BB Flags : ADA = %d  ADC = %d ",nBBflagsADA, nBBflagsADC )); 
	//AliInfo(Form("BG Flags : ADA = %d  ADC = %d ",nBGflagsADA, nBGflagsADC )); 
	//AliInfo(Form("Charges  : ADA = %d  ADC = %d ",chargeADA, chargeADC )); 
	
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


