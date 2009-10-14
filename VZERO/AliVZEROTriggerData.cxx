
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

// Class AliVZEROTriggerData
// -------------------------
// Retrieves and hold the FEE parameters
// The parameters are recieved from the shuttle 
// AliVZEROTriggerData is then used in the AliVZEROTriggerSimulator
//

#include <TObjString.h>
#include <TMap.h>

#include "AliLog.h"
#include "AliDCSValue.h"
#include "AliVZEROTriggerData.h"

ClassImp(AliVZEROTriggerData)
//________________________________________________________________

AliVZEROTriggerData::AliVZEROTriggerData() :
	TNamed(),
	fBBAThreshold(0),
	fBBCThreshold(0) ,  
	fBGAThreshold(0) ,  
	fBGCThreshold(0) ,  
	fBBAForBGThreshold(0) ,  
	fBBCForBGThreshold(0) ,  
	fCentralityVOAThrLow(0) ,  
	fCentralityVOAThrHigh(0) , 
	fCentralityVOCThrLow(0) ,  
	fCentralityVOCThrHigh(0) , 
	fMultV0AThrLow(0) ,  
	fMultV0AThrHigh(0) , 
	fMultV0CThrLow(0) ,  
	fMultV0CThrHigh(0),
	fRun(0),
	fStartTime(0),
	fEndTime(0),
	fIsProcessed(kFALSE)	

{
	// default constructor
	for(int i=0; i<kNCIUBoards ;i++) {
		fClk1Win1[i] = fClk1Win2[i] = 0;
		fDelayClk1Win1[i] = fDelayClk1Win2[i] = 0;
		fClk2Win1[i] = fClk2Win2[i] = 0;
		fDelayClk2Win1[i] = fDelayClk2Win2[i] = 0;
		fLatchWin1[i] = fLatchWin2[i] = 0;
		fResetWin1[i] = fResetWin2[i] = 0;
		fPedestalSubtraction[i] = kFALSE;
		fTimeResolution[i] = 25./256.; // Default time resolution
		fWidthResolution[i] = 0.4;     // Default time width resolution
	}
}
//________________________________________________________________
AliVZEROTriggerData::AliVZEROTriggerData(Int_t nRun, UInt_t startTime, UInt_t endTime) :
	TNamed(),
	fBBAThreshold(0),
	fBBCThreshold(0) ,  
	fBGAThreshold(0) ,  
	fBGCThreshold(0) ,  
	fBBAForBGThreshold(0) ,  
	fBBCForBGThreshold(0) ,  
	fCentralityVOAThrLow(0) ,  
	fCentralityVOAThrHigh(0) , 
	fCentralityVOCThrLow(0) ,  
	fCentralityVOCThrHigh(0) , 
	fMultV0AThrLow(0) ,  
	fMultV0AThrHigh(0) , 
	fMultV0CThrLow(0) ,  
	fMultV0CThrHigh(0),
	fRun(nRun),
	fStartTime(startTime),
	fEndTime(endTime),
	fIsProcessed(kFALSE)
{
	// Constructor
	for(int i=0; i<kNCIUBoards ;i++) {
		fClk1Win1[i] = fClk1Win2[i] = 0;
		fDelayClk1Win1[i] = fDelayClk1Win2[i] = 0;
		fClk2Win1[i] = fClk2Win2[i] = 0;
		fDelayClk2Win1[i] = fDelayClk2Win2[i] = 0;
		fLatchWin1[i] = fLatchWin2[i] = 0;
		fResetWin1[i] = fResetWin2[i] = 0;
		fPedestalSubtraction[i] = kFALSE;
		fTimeResolution[i] = 25./256.; // Default time resolution
		fWidthResolution[i] = 0.4;     // Default time width resolution
	}
	TString namst = "Trigger_FEE";
	SetName(namst.Data());
	SetTitle(namst.Data());
	
}

//________________________________________________________________
AliVZEROTriggerData::~AliVZEROTriggerData(){
	// destructor
}
//_____________________________________________________________________________
void AliVZEROTriggerData::FillData(AliVZERODataFEE * data){
	// Set all parameters from the data get by the shuttle
	TMap * params = data->GetParameters();
	TIter iter(params);	
	TObjString* aliasName;
	
	while ((  aliasName = (TObjString*) iter.Next() ))  {
		AliDCSValue* aValue = (AliDCSValue*) params->GetValue(aliasName);
		Float_t val;
		if(aValue) {
			val = aValue->GetFloat();
			//AliInfo(Form("%s : %f",aliasName->String().Data(), val));
			SetParameter(aliasName->String(),val);
		}
	}	
}

//_____________________________________________________________________________
void AliVZEROTriggerData::SetParameter(TString name, Float_t val){
	// Set given parameter
	
	Int_t iBoard = -1;
	Int_t iChannel = -1;

	TSeqCollection* nameSplit = name.Tokenize("/");
	TObjString * boardName = (TObjString *)nameSplit->At(2);
	if(!boardName->String().Contains("CCIU")) sscanf(boardName->String().Data(),"CIU%d",&iBoard);
	
	TString paramName = ((TObjString *)nameSplit->At(3))->String();
	Char_t channel[2] ; channel[1] = '\0';
	channel[0] = paramName[paramName.Sizeof()-2];
	sscanf(channel,"%d",&iChannel);
	
	if(name.Contains("DelayClk1Win1")) SetDelayClk1Win1((UShort_t) val,iBoard-1);
	else if(name.Contains("Clk1Win1")) SetClk1Win1((UShort_t) val,iBoard-1);
	else if(name.Contains("DelayClk1Win2")) SetDelayClk1Win2((UShort_t) val,iBoard-1);
	else if(name.Contains("Clk1Win2")) SetClk1Win2((UShort_t) val,iBoard-1);
	else if(name.Contains("DelayClk2Win1")) SetDelayClk2Win1((UShort_t) val,iBoard-1);
	else if(name.Contains("Clk2Win1")) SetClk2Win1((UShort_t) val,iBoard-1);
	else if(name.Contains("DelayClk2Win2")) SetDelayClk2Win2((UShort_t) val,iBoard-1);
	else if(name.Contains("Clk2Win2")) SetClk2Win2((UShort_t) val,iBoard-1);
	else if(name.Contains("LatchWin1")) SetLatchWin1((UShort_t) val,iBoard-1);
	else if(name.Contains("LatchWin2")) SetLatchWin2((UShort_t) val,iBoard-1);
	else if(name.Contains("ResetWin1")) SetResetWin1((UShort_t) val,iBoard-1);
	else if(name.Contains("ResetWin2")) SetResetWin2((UShort_t) val,iBoard-1);
	else if(name.Contains("PedestalSubtraction")) SetPedestalSubtraction((Bool_t) val,iBoard-1);
	else if(name.Contains("TimeResolution")) SetTimeResolution((UShort_t) val,iBoard-1);
	else if(name.Contains("WidthResolution")) SetWidthResolution((UShort_t) val,iBoard-1);
	else if(name.Contains("BBAThreshold")) SetBBAThreshold((UShort_t) val);
	else if(name.Contains("BBCThreshold")) SetBBCThreshold((UShort_t) val);
	else if(name.Contains("BGAThreshold")) SetBGAThreshold((UShort_t) val);
	else if(name.Contains("BGCThreshold")) SetBGCThreshold((UShort_t) val);
	else if(name.Contains("BBAForBGThreshold")) SetBBAForBGThreshold((UShort_t) val);
	else if(name.Contains("BBCForBGThreshold")) SetBBCForBGThreshold((UShort_t) val);
	else if(name.Contains("CentralityV0AThrLow")) SetCentralityV0AThrLow((UShort_t) val);
	else if(name.Contains("CentralityV0AThrHigh")) SetCentralityV0AThrHigh((UShort_t) val);
	else if(name.Contains("CentralityV0CThrLow")) SetCentralityV0CThrLow((UShort_t) val);
	else if(name.Contains("CentralityV0CThrHigh")) SetCentralityV0CThrHigh((UShort_t) val);
	else if(name.Contains("MultV0AThrLow")) SetMultV0AThrLow((UShort_t) val);
	else if(name.Contains("MultV0AThrHigh")) SetMultV0AThrHigh((UShort_t) val);
	else if(name.Contains("MultV0CThrLow")) SetMultV0CThrLow((UShort_t) val);
	else if(name.Contains("MultV0CThrHigh")) SetMultV0CThrHigh((UShort_t) val);
	else if(name.Contains("TriggerSelect")) SetTriggerSelected((UShort_t) val, iBoard -1);
	else if(name.Contains("EnableCharge")) SetEnableCharge((Bool_t) val, iBoard -1, iChannel-1);
	else if(name.Contains("EnableTiming")) SetEnableTiming((Bool_t) val, iBoard -1, iChannel-1);
	else if(name.Contains("DiscriThr")) SetDiscriThr((UShort_t) val, iBoard -1, iChannel-1);
	else if(name.Contains("DelayHit")) SetDelayHit((UShort_t) val, iBoard -1, iChannel-1);
	else if(name.Contains("PedOdd")) SetPedestal((UShort_t) val, 1, iBoard -1, iChannel-1);
	else if(name.Contains("PedEven")) SetPedestal((UShort_t) val, 0, iBoard -1, iChannel-1);
	else if(name.Contains("PedCutOdd")) SetPedestalCut((UShort_t) val, 1, iBoard -1, iChannel-1);
	else if(name.Contains("PedCutEven")) SetPedestalCut((UShort_t) val, 0, iBoard -1, iChannel-1);
	else AliError(Form("No Setter found for FEE parameter : %s",name.Data()));
}
//________________________________________________________________
void AliVZEROTriggerData::SetPedestalCut(UShort_t val,Int_t integrator, Int_t board, Int_t channel)
{
	// Set Pedestal Cut of individual channel 
	if(board<kNCIUBoards && channel<kNChannels) {
		if(integrator) fPedestalCutOdd[board][channel] = val;
		else fPedestalCutEven[board][channel] = val;
	} else AliError(Form("Impossible to write at : Board %d ; Channel %d",board,channel));
}
//________________________________________________________________
UShort_t AliVZEROTriggerData::GetPedestalCut(Int_t integrator, Int_t board, Int_t channel)
{
	// Get Pedestal Cut of individual channel 
	if(board<kNCIUBoards && channel<kNChannels) {
		if(integrator) return(fPedestalCutOdd[board][channel]);
		else return(fPedestalCutEven[board][channel]);
	}else AliError(Form("Impossible to read at : Board %d ; Channel %d",board,channel));
	return 0;
}
//________________________________________________________________
void AliVZEROTriggerData::SetPedestal(UShort_t val, Int_t integrator, Int_t board, Int_t channel)
{
	// Set Pedestal of individual channel 
	if(board<kNCIUBoards && channel<kNChannels) {
		if(integrator) fPedestalOdd[board][channel] = val;
		else fPedestalEven[board][channel] = val;
	} else AliError(Form("Impossible to write at : Board %d ; Channel %d ; Integrator %d ",board,channel,integrator));
}
//________________________________________________________________
UShort_t AliVZEROTriggerData::GetPedestal(Int_t integrator, Int_t board, Int_t channel)
{
	// Get Pedestal of individual channel 
	if(board<kNCIUBoards && channel<kNChannels) {
		if(integrator) return(fPedestalOdd[board][channel]);
		else return(fPedestalEven[board][channel]);
	} else AliError(Form("Impossible to read at : Board %d ; Channel %d",board,channel));
	return 0;
}
//________________________________________________________________
void AliVZEROTriggerData::SetDelayHit(UShort_t val,Int_t board, Int_t channel)
{
	// Set Delay of individual channel 
	if(board<kNCIUBoards && channel<kNChannels) fDelayHit[board][channel] = val;
	else AliError(Form("Impossible to write at : Board %d ; Channel %d",board,channel));
}
//________________________________________________________________
UShort_t AliVZEROTriggerData::GetDelayHit(Int_t board, Int_t channel)
{
	// Get Delay of individual channel 
	if(board<kNCIUBoards && channel<kNChannels) return(fDelayHit[board][channel]);
	else AliError(Form("Impossible to read at : Board %d ; Channel %d",board,channel));
	return 0;
}
//________________________________________________________________
void AliVZEROTriggerData::SetDiscriThr(UShort_t val,Int_t board, Int_t channel)
{
	// Set discriminator threshold
	if(board<kNCIUBoards && channel<kNChannels) fDiscriThr[board][channel] = val;
	else AliError(Form("Impossible to write at : Board %d ; Channel %d",board,channel));
}
//________________________________________________________________
UShort_t AliVZEROTriggerData::GetDiscriThr(Int_t board, Int_t channel)
{
	// Get discriminator threshold
	if(board<kNCIUBoards && channel<kNChannels) return(fDiscriThr[board][channel]);
	else AliError(Form("Impossible to read at : Board %d ; Channel %d",board,channel));
	return 0;
}
//________________________________________________________________
void AliVZEROTriggerData::SetEnableCharge(Bool_t val,Int_t board, Int_t channel)
{
	// Set the channels enabled for Charge triggers
	if(board<kNCIUBoards && channel<kNChannels) fEnableCharge[board][channel] = val;
	else AliError(Form("Impossible to write at : Board %d ; Channel %d",board,channel));
}
//________________________________________________________________
Bool_t AliVZEROTriggerData::GetEnableCharge(Int_t board, Int_t channel)
{
	// Get the channels enabled for Charge triggers
	if(board<kNCIUBoards && channel<kNChannels) return(fEnableCharge[board][channel]);
	else AliError(Form("Impossible to read at : Board %d ; Channel %d",board,channel));
	return kFALSE;
}
//________________________________________________________________
void AliVZEROTriggerData::SetEnableTiming(Bool_t val,Int_t board, Int_t channel)
{
	// Set the channels enabled for Timing triggers
	if(board<kNCIUBoards && channel<kNChannels) fEnableTiming[board][channel] = val;
	else AliError(Form("Impossible to write at : Board %d ; Channel %d",board,channel));
}
//________________________________________________________________
Bool_t AliVZEROTriggerData::GetEnableTiming(Int_t board, Int_t channel)
{
	// Get the channels enabled for Timing triggers
	if(board<kNCIUBoards && channel<kNChannels) return(fEnableTiming[board][channel]);
	else AliError(Form("Impossible to read at : Board %d ; Channel %d",board,channel));
	return kFALSE;
}
//________________________________________________________________
void AliVZEROTriggerData::SetTriggerSelected(UShort_t trigger, Int_t output)
{
	// Set the trigger selected on the outputs to CTP
	if(output<kNTriggerOutputs) fTriggerSelected[output] = trigger;
	else AliError(Form("Trigger output number %d not valid",output));
}

//________________________________________________________________
void AliVZEROTriggerData::SetClk1Win1(UShort_t* clks)
{
	// Set Win clock of BB
	if(clks) for(int t=0; t<kNCIUBoards; t++) SetClk1Win1(clks[t],t);
	else AliFatal("Profil Clock1 Win1 Not defined.");
}
//________________________________________________________________
void AliVZEROTriggerData::SetClk2Win1(UShort_t* clks)
{
	// Set Win clock of BB
	if(clks) for(int t=0; t<kNCIUBoards; t++) SetClk2Win1(clks[t],t);
	else AliFatal("Profil Clock2 Win1 Not defined.");
}
//________________________________________________________________
void AliVZEROTriggerData::SetClk1Win1(UShort_t clk, Int_t board)
{
	// Set Win clock of BB
	if(IsClkValid(clk) && (board<kNCIUBoards)) fClk1Win1[board] = clk;
	else AliError(Form("Profil Clock1 Win1 of board %d is not valid : %d",board,clk));
}
//________________________________________________________________
void AliVZEROTriggerData::SetClk2Win1(UShort_t clk, Int_t board)
{
	// Set Win clock of BB
	if(IsClkValid(clk) && (board<kNCIUBoards)) fClk2Win1[board] = clk;
	else AliError(Form("Profil Clock2 Win1 of board %d is not valid",board));
}
//________________________________________________________________
void AliVZEROTriggerData::SetClk1Win2(UShort_t* clks)
{
	// Set Win clock of BG
	if(clks) for(int t=0; t<kNCIUBoards; t++) SetClk1Win2(clks[t],t);
	else AliFatal("Profil Clock1 Win2 Not defined.");
}
//________________________________________________________________
void AliVZEROTriggerData::SetClk2Win2(UShort_t* clks)
{
	// Set Win clock of BG
	if(clks) for(int t=0; t<kNCIUBoards; t++) SetClk2Win2(clks[t],t);
	else AliFatal("Profil Clock2 Win2 Not defined.");
}
//________________________________________________________________
void AliVZEROTriggerData::SetClk1Win2(UShort_t clk, Int_t board)
{
	// Set Win clock of BG
	if(IsClkValid(clk) && (board<kNCIUBoards)) fClk1Win2[board] = clk;
	else AliError(Form("Profil Clock1 Win2 of board %d is not valid",board));
}
//________________________________________________________________
void AliVZEROTriggerData::SetClk2Win2(UShort_t clk, Int_t board)
{
	// Set Win clock of BG
	if(IsClkValid(clk) && (board<kNCIUBoards)) fClk2Win2[board] = clk;
	else AliError(Form("Profil Clock2 Win2 of board %d is not valid",board));
}
//________________________________________________________________
void AliVZEROTriggerData::SetDelayClk1Win1(UShort_t* delays)
{
	// Set Delay for Win clock of BB
	if(delays) for(int t=0; t<kNCIUBoards; t++) SetDelayClk1Win1(delays[t],t);
	else AliFatal("Profil Clock1 Win1 Delays Not defined.");
}
//________________________________________________________________
void AliVZEROTriggerData::SetDelayClk1Win1(UShort_t delay, Int_t board)
{
	// Set Delay for Win clock of BB
	if(delay>1023){
		AliWarning(Form("Profil Clock1 Win1 Delay of board %d should be less 1023 is currently %d. Truncated to the first 10 bits",board, delay));
		delay = delay & 0x3FF;
	}
	if(board<kNCIUBoards)	fDelayClk1Win1[board] = delay;
	else AliError("Trying to write out of the array");
}
//________________________________________________________________
void AliVZEROTriggerData::SetDelayClk2Win1(UShort_t* delays)
{
	// Set Delay for Win clock of BB
	if(delays) for(int t=0; t<kNCIUBoards; t++) SetDelayClk2Win1(delays[t],t);
	else AliFatal("Profil Clock2 Win1 Delays Not defined.");
}
//________________________________________________________________
void AliVZEROTriggerData::SetDelayClk2Win1(UShort_t delay, Int_t board)
{
	// Set Delay for Win clock of BB
	if(delay>1023){
		AliWarning(Form("Profil Clock2 Win1 Delay of board %d should be less 1023 is currently %d. Truncated to the first 10 bits",board, delay));
		delay = delay & 0x3FF;
	}
	if(board<kNCIUBoards)	fDelayClk2Win1[board] = delay;
	else AliError("Trying to write out of the array");
}
//________________________________________________________________
void AliVZEROTriggerData::SetDelayClk1Win2(UShort_t* delays)
{
	// Set Delay for Win clock of BG
	if(delays) for(int t=0; t<kNCIUBoards; t++) SetDelayClk1Win2(delays[t],t);
	else AliFatal("Profil Clock1 Win2 Delays Not defined.");
}
//________________________________________________________________
void AliVZEROTriggerData::SetDelayClk1Win2(UShort_t delay, Int_t board)
{
	// Set Delay for Win clock of BG
	if(delay>1023){
		AliWarning(Form("Profil Clock1 Win2 Delay of board %d should be less 1023 is currently %d. Truncated to the first 10 bits",board, delay));
		delay = delay & 0x3FF;
	}
	if(board<kNCIUBoards)	fDelayClk1Win2[board] = delay;
	else AliError("Trying to write out of the array");
}
//________________________________________________________________
void AliVZEROTriggerData::SetDelayClk2Win2(UShort_t* delays)
{
	// Set Delay for Win clock of BG
	if(delays) for(int t=0; t<kNCIUBoards; t++) SetDelayClk2Win2(delays[t],t);
	else AliFatal("Profil Clock2 Win2 Delays Not defined.");
}
//________________________________________________________________
void AliVZEROTriggerData::SetDelayClk2Win2(UShort_t delay, Int_t board)
{
	// Set Delay for Win clock of BG
	if(delay>1023){
		AliWarning(Form("Profil Clock2 Win2 Delay of board %d should be less 1023 is currently %d. Truncated to the first 10 bits",board, delay));
		delay = delay & 0x3FF;
	}
	if(board<kNCIUBoards)	fDelayClk2Win2[board] = delay;
	else AliError("Trying to write out of the array");
}
//________________________________________________________________
void AliVZEROTriggerData::SetLatchWin1(UShort_t *latchs){
	// Set Latch Win clock for BB
	if(latchs) for(int t=0; t<kNCIUBoards; t++) SetLatchWin1(latchs[t],t);
	else AliFatal("Latch Win1 profil Not defined.");
}
//________________________________________________________________
void AliVZEROTriggerData::SetLatchWin1(UShort_t latch, Int_t board)
{
	// Set Latch Win clock for BB
	if(IsClkValid(latch) && (board<kNCIUBoards)) fLatchWin1[board] = latch;
	else AliError(Form("Latch Win1 profil of board %d is not valid",board));
}
//________________________________________________________________
void AliVZEROTriggerData::SetLatchWin2(UShort_t *latchs){
	// Set Latch Win clock for BG
	if(latchs) for(int t=0; t<kNCIUBoards; t++) SetLatchWin2(latchs[t],t);
	else AliFatal("Latch Win2 profil Not defined.");
}
//________________________________________________________________
void AliVZEROTriggerData::SetLatchWin2(UShort_t latch, Int_t board)
{
	// Set Latch Win clock for BG
	if(IsClkValid(latch) && (board<kNCIUBoards)) fLatchWin2[board] = latch;
	else AliError(Form("Latch Win2 profil of board %d is not valid",board));
}
//________________________________________________________________
void AliVZEROTriggerData::SetResetWin1(UShort_t *resets){
	// Set Reset Win clock for BB
	if(resets) for(int t=0; t<kNCIUBoards; t++) SetResetWin1(resets[t],t);
	else AliFatal("Reset Win1 profil Not defined.");
}
//________________________________________________________________
void AliVZEROTriggerData::SetResetWin1(UShort_t reset, Int_t board)
{
	// Set Reset Win clock for BB
	if(IsClkValid(reset) && (board<kNCIUBoards)) fResetWin1[board] = reset;
	else AliError(Form("Reset Win1 profil of board %d is not valid",board));
}
//________________________________________________________________
void AliVZEROTriggerData::SetResetWin2(UShort_t *resets){
	// Set Reset Win clock for BG
	if(resets)  for(int t=0; t<kNCIUBoards; t++) SetResetWin2(resets[t],t);
	else AliFatal("Reset Win2 profil Not defined.");
}
//________________________________________________________________
void AliVZEROTriggerData::SetResetWin2(UShort_t reset, Int_t board)
{
	// Set Reset Win clock for BG
	if(IsClkValid(reset) && (board<kNCIUBoards)) fResetWin2[board] = reset;
	else AliError(Form("Reset Win2 profil of board %d is not valid",board));
}
//________________________________________________________________
void AliVZEROTriggerData::SetPedestalSubtraction(Bool_t *peds){
	// Set Pedestal Subtraction Parameter
	if(peds)  for(int t=0; t<kNCIUBoards; t++) SetPedestalSubtraction(peds[t],t);
	else AliFatal("Pedestal Subtraction Not defined.");
	
}
//________________________________________________________________
void AliVZEROTriggerData::SetPedestalSubtraction(Bool_t ped, Int_t board)
{
	// Set Pedestal Subtraction Parameter
	if((board<kNCIUBoards)) fPedestalSubtraction[board] = ped;
	else AliError(Form("Board %d is not valid",board));
}
//________________________________________________________________
void AliVZEROTriggerData::SetTimeResolution(UShort_t *resols){
	// Set Time Resolution of the TDC
	if(resols)  for(int t=0; t<kNCIUBoards; t++) SetTimeResolution(resols[t],t);
	else AliFatal("Time Resolution not defined.");
	
}
//________________________________________________________________
void AliVZEROTriggerData::SetTimeResolution(UShort_t resol, Int_t board)
{
	// Set Time Resolution of the TDC
	if((board<kNCIUBoards)) {
		switch(resol){
			case 0:
				fTimeResolution[board] = 25./256.;
				break;
			case 1:
				fTimeResolution[board] = 25./128.;
				break;
			case 2:
				fTimeResolution[board] = 25./64.;
				break;
			case 3:
				fTimeResolution[board] = 25./32.;
				break;
			case 4:
				fTimeResolution[board] = 25./16.;
				break;
			case 5:
				fTimeResolution[board] = 25./8.;
				break;
			case 6:
				fTimeResolution[board] = 6.25;
				break;
			case 7:
				fTimeResolution[board] = 12.5;
				break;
		}
	} else AliError(Form("Board %d is not valid",board));
}
//________________________________________________________________
void AliVZEROTriggerData::SetWidthResolution(UShort_t *resols){
	// Set Time Width Resolution of the TDC
	if(resols)  for(int t=0; t<kNCIUBoards; t++) SetWidthResolution(resols[t],t);
	else AliFatal("Width Resolution not defined.");
	
}
//________________________________________________________________
void AliVZEROTriggerData::SetWidthResolution(UShort_t resol, Int_t board)
{
	// Set Time Width Resolution of the TDC
	if((board<kNCIUBoards)){
		switch(resol){
			case 0:
				fWidthResolution[board] = 25./256.;
				break;
			case 1:
				fWidthResolution[board] = 25./128.;
				break;
			case 2:
				fWidthResolution[board] = 25./64.;
				break;
			case 3:
				fWidthResolution[board] = 25./32.;
				break;
			case 4:
				fWidthResolution[board] = 25./16.;
				break;
			case 5:
				fWidthResolution[board] = 25./8.;
				break;
			case 6:
				fWidthResolution[board] = 6.25;
				break;
			case 7:
				fWidthResolution[board] = 12.5;
				break;
			case 8:
				fWidthResolution[board] = 25.;
				break;
			case 9:
				fWidthResolution[board] = 50.;
				break;
			case 10:
				fWidthResolution[board] = 100.;
				break;
			case 11:
				fWidthResolution[board] = 200.;
				break;
			case 12:
				fWidthResolution[board] = 400.;
				break;
			case 13:
				fWidthResolution[board] = 800.;
				break;
				
		}
	}else AliError(Form("Board %d is not valid",board));
}

//________________________________________________________________
Bool_t	AliVZEROTriggerData::IsClkValid(UShort_t clock) const {
	// Check if the given clock has a valid profil.
	Bool_t word[5];
	Bool_t isValid = kTRUE;
	Short_t risingEdge = 0;
	Short_t fallingEdge = 0;
	for(int i=0 ; i<5 ; i++) word[i] = (clock >> i) & 0x1;
	
	if(word[0] != word[4]){
		if(word[4]) fallingEdge++;
		else risingEdge++;
	}	
	for(int i=1 ; i<5 ; i++){
		if(word[i] != word[i-1]) {
			if(word[i-1]) fallingEdge++;
			else risingEdge++;
		}
	}
	if((fallingEdge>1)||(risingEdge>1)) isValid = kFALSE;
	if(((risingEdge==0)&&(fallingEdge==0)) &&(!word[0]))  isValid = kFALSE;
	return isValid;
}



