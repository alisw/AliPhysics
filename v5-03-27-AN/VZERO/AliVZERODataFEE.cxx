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
// ---------------------
// Class AliVZERODataFEE
// ---------------------
// Used to process the TMap of DCS values comming from the shuttle.
// It stores into a TMap the FEE parameters for the given run number
//

#include <TTimeStamp.h>
#include <TObjString.h>
#include <TMap.h>

#include "AliDCSValue.h"
#include "AliLog.h"
#include "AliVZERODataFEE.h"

ClassImp(AliVZERODataFEE)

//_____________________________________________________________________________
AliVZERODataFEE::AliVZERODataFEE() :
	TObject(),
	fRun(0),
	fStartTime(0),
	fEndTime(0),
	fIsProcessed(kFALSE),
	fParameters(NULL)
{
	// Default constructor
}

//_____________________________________________________________________________
AliVZERODataFEE::AliVZERODataFEE(Int_t nRun, UInt_t startTime, UInt_t endTime) : 
	TObject(),
	fRun(nRun),
	fStartTime(startTime),
	fEndTime(endTime),
	fIsProcessed(kFALSE),
	fParameters(new TMap())
{
	// Constructor
	AliInfo(Form("\n\tRun %d \n\tStartTime %s \n\tEndTime %s", nRun,
				 TTimeStamp(startTime).AsString(),
				 TTimeStamp(endTime).AsString()));
	fParameters->SetOwnerValue();
	Init();
}

//_____________________________________________________________________________
AliVZERODataFEE::~AliVZERODataFEE()
{
	delete fParameters;
}

//_____________________________________________________________________________
void AliVZERODataFEE::Init(){
	// Initialization of DCS aliases
	int iAlias = 0;

	// CCIU Parameters
	
	fAliasNames[iAlias++] = "V00/FEE/CCIU/BBAThreshold";
	fAliasNames[iAlias++] = "V00/FEE/CCIU/BBCThreshold";
	fAliasNames[iAlias++] = "V00/FEE/CCIU/BGAThreshold";
	fAliasNames[iAlias++] = "V00/FEE/CCIU/BGCThreshold";
	fAliasNames[iAlias++] = "V00/FEE/CCIU/BBAForBGThreshold";
	fAliasNames[iAlias++] = "V00/FEE/CCIU/BBCForBGThreshold";
	fAliasNames[iAlias++] = "V00/FEE/CCIU/CentralityV0AThrLow";
	fAliasNames[iAlias++] = "V00/FEE/CCIU/CentralityV0AThrHigh";
	fAliasNames[iAlias++] = "V00/FEE/CCIU/CentralityV0CThrLow";
	fAliasNames[iAlias++] = "V00/FEE/CCIU/CentralityV0CThrHigh";
	fAliasNames[iAlias++] = "V00/FEE/CCIU/MultV0AThrLow";
	fAliasNames[iAlias++] = "V00/FEE/CCIU/MultV0AThrHigh";
	fAliasNames[iAlias++] = "V00/FEE/CCIU/MultV0CThrLow";
	fAliasNames[iAlias++] = "V00/FEE/CCIU/MultV0CThrHigh";
	for(int i=1;i<=5;i++) {
		fAliasNames[iAlias] = "V00/FEE/CCIU/TriggerSelect";
		fAliasNames[iAlias++] += Form("%d",i);
	}
	
	// CIU  Parameters
	
	for(int iCIU = 0; iCIU<8 ; iCIU++){
		for(int iParam=0; iParam<kNCIUParam;iParam++){
			fAliasNames[iAlias] = "V00/FEE/";
			fAliasNames[iAlias] += Form("CIU%d/",iCIU);
			
			fAliasNames[iAlias] += GetFEEParamName(iParam);
			iAlias++;
		}
		for(int iParam=kNCIUParam; iParam<kNCIUParam+kNChannelParam;iParam++){
			for(int iCh=1;iCh<=8;iCh++){
				fAliasNames[iAlias] = "V00/FEE/";
				fAliasNames[iAlias] += Form("CIU%d/",iCIU);
			
				fAliasNames[iAlias] += GetFEEParamName(iParam);
				fAliasNames[iAlias] += Form("%d",iCh);
			
				iAlias++;
			}
		}
		
	}

	if(iAlias!=kNAliases) 
		AliError(Form("Number of FEE Aliases defined not correct"));
	
}
//_____________________________________________________________________________
void AliVZERODataFEE::ProcessData(TMap& aliasMap){
// Process the TMap of DCS data
	
	if(!(fAliasNames[0])) Init();
	
	TObjArray *aliasArr;
	AliDCSValue* aValue;

	// starting loop on aliases
	for(int iAlias=0; iAlias<kNAliases; iAlias++){
		
		aliasArr = (TObjArray*) aliasMap.GetValue(fAliasNames[iAlias].Data());
		if(!aliasArr){
			AliError(Form("Alias %s not found!", fAliasNames[iAlias].Data()));
			return;
		}
				
		if(aliasArr->GetEntries()<1){
			AliError(Form("Alias %s has no entries!", fAliasNames[iAlias].Data()));
			continue;
		}
		
		TIter iterarray(aliasArr);
				
		AliDCSValue * lastVal = NULL;
		while((aValue = (AliDCSValue*) iterarray.Next())) lastVal = aValue; // Take only the last value
		
		fParameters->Add(new TObjString(fAliasNames[iAlias].Data()),lastVal);
		
	}
	
	fIsProcessed=kTRUE;
	
}
//_____________________________________________________________________________

TString AliVZERODataFEE::GetFEEParamName(Int_t iParam){
// Return the name of the FEE Parameter iParam
	TString result;
	if(iParam>kNCIUParam + kNChannelParam -1) {
		AliError(Form("Requesting FEE parameter number %d. Max parameter number is : %d",iParam,kNCIUParam + kNChannelParam-1));
		return result;
	}
	switch (iParam) {
		case 0: result = "Clk1Win1"; break;
		case 1: result = "Clk2Win1"; break;
		case 2: result = "Clk1Win2"; break;
		case 3: result = "Clk2Win2"; break;
		case 4: result = "DelayClk1Win1"; break;
		case 5: result = "DelayClk2Win1"; break;
		case 6: result = "DelayClk1Win2"; break;
		case 7: result = "DelayClk2Win2"; break;
		case 8: result = "LatchWin1"; break;
		case 9: result = "LatchWin2"; break;
		case 10: result = "ResetWin1"; break;
		case 11: result = "ResetWin2"; break;
		case 12: result = "PedestalSubtraction"; break;
		case 13: result = "EnableCharge"; break;
		case 14: result = "EnableTiming"; break;
		case 15: result = "PedEven"; break;
		case 16: result = "PedOdd"; break;
		case 17: result = "PedCutEven"; break;
		case 18: result = "PedCutOdd"; break;
		case 19: result = "DelayHit"; break;
		case 20: result = "DiscriThr"; break;
	}
	return result;
}

void AliVZERODataFEE::PrintAliases(){
	for(int i=0;i<kNAliases;i++) AliInfo(Form("%s",fAliasNames[i].Data()));
}

