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


#include "AliADDataDCS.h"

#include "AliDCSValue.h"
#include "AliLog.h"
#include "AliADConst.h"

#include <TGraph.h>
#include <TAxis.h>
#include <TCanvas.h>
#include <TTimeStamp.h>
#include <TMap.h>
#include <TString.h>
#include <TObjString.h>
#include <TH1F.h>

class TH2;
class AliCDBMetaData;
class TDatime;

// AliADDataDCS class
// main aim to introduce the aliases for the AD DCS
// data points to be then
// stored in the OCDB, and to process them. 
// ProcessData() method called by ADPreprocessor

ClassImp(AliADDataDCS)

//_____________________________________________________________________________
AliADDataDCS::AliADDataDCS():
	TObject(),
	fRun(0),
	fStartTime(0),
	fEndTime(0),
	fDaqStartTime(0),
	fDaqEndTime(0),
	fCtpStartTime(0),
	fCtpEndTime(0),
        fGraphs(NULL),
	fFEEParameters(NULL),
	fIsProcessed(kFALSE)

{
  // Default constructor
	for(int i=0;i<kNHvChannel;i++) {
		fDeadChannel[i] = kFALSE;
		fMeanHV[i]      = 100.0;
		fWidthHV[i]     = 0.0; 
		fHv[i]          = NULL;
	}
}

//_____________________________________________________________________________
AliADDataDCS::AliADDataDCS(Int_t nRun, UInt_t startTime, UInt_t endTime, UInt_t daqStartTime, UInt_t daqEndTime, UInt_t ctpStartTime, UInt_t ctpEndTime):
	TObject(),
	fRun(nRun),
	fStartTime(startTime),
	fEndTime(endTime),
	fDaqStartTime(daqStartTime),
	fDaqEndTime(daqEndTime),
	fCtpStartTime(ctpStartTime),
	fCtpEndTime(ctpEndTime),
	fGraphs(new TClonesArray("TGraph",kNGraphs)),
	fFEEParameters(new TMap()),
	fIsProcessed(kFALSE)

{

  // constructor with arguments
  	for(int i=0;i<kNHvChannel;i++) {
	 fDeadChannel[i] = kFALSE;        
	 fMeanHV[i]      = 100.0;
     fWidthHV[i]     = 0.0; 
	}
	AliInfo(Form("\n\tRun %d \n\tTime Created %s \n\tTime Completed %s \n\tDAQ start %s \n\tDAQ end %s \n\tCTP start %s \n\tCTP end %s   ", nRun,
		TTimeStamp(startTime).AsString(),
		TTimeStamp(endTime).AsString(),
		TTimeStamp(daqStartTime).AsString(),
		TTimeStamp(daqEndTime).AsString(),
		TTimeStamp(ctpStartTime).AsString(),
		TTimeStamp(ctpEndTime).AsString()
		));
	
	fFEEParameters->SetOwnerValue();
	Init();

}

//_____________________________________________________________________________
AliADDataDCS::~AliADDataDCS() {

  // destructor
  fGraphs->Clear("C");
  delete fFEEParameters;

}

//_____________________________________________________________________________
Bool_t AliADDataDCS::ProcessData(TMap& aliasMap){

  // method to process the data
  Bool_t success = kTRUE;

  if(!(fAliasNames[0])) Init();

  TObjArray *aliasArr;
  AliDCSValue* aValue;

  // starting loop on aliases
  for(int iAlias=0; iAlias<kNAliases; iAlias++){

    aliasArr = (TObjArray*) aliasMap.GetValue(fAliasNames[iAlias].Data());
    if(!aliasArr){
      AliError(Form("Missing data points for alias %s!", fAliasNames[iAlias].Data()));
      success = kFALSE;
      continue;
    }

    //Introduce(iAlias, aliasArr);
    
    if(aliasArr->GetEntries()<2){
      AliWarning(Form("Alias %s has just %d data points!",
		    fAliasNames[iAlias].Data(),aliasArr->GetEntries()));
    }
    
    TIter iterarray(aliasArr);
	
    if(iAlias<kNHvChannel){ // Treating HV values
    	Int_t nentries = aliasArr->GetEntries();

    	Double_t *times = new Double_t[nentries];
    	Double_t *values = new Double_t[nentries];

    	UInt_t iValue=0;
    	Float_t variation = 0.0;

    	while((aValue = (AliDCSValue*) iterarray.Next())) {
			UInt_t currentTime = aValue->GetTimeStamp();
			if(currentTime>fCtpEndTime) break;

   			values[iValue] = aValue->GetFloat();
   			times[iValue] = (Double_t) (currentTime);
			
			if(iValue>0) {
				if(values[iValue-1]>0.) variation = TMath::Abs(values[iValue]-values[iValue-1])/values[iValue-1];
				if(variation > 0.05) fDeadChannel[iAlias] = kTRUE;
				}
			fHv[iAlias]->Fill(values[iValue]);
   			iValue++;
    			}      
    	CreateGraph(iAlias, iValue, times, values); // fill graphs 

  	// calculate mean and rms of the first two histos
	Int_t iChannel	   = iAlias;	
	fMeanHV[iAlias]  = fHv[iAlias]->GetMean();
	fWidthHV[iAlias] = fHv[iAlias]->GetRMS();

    	delete[] values;
    	delete[] times;	
	}
      else if(iAlias >= kNHvChannel && iAlias<2*kNHvChannel){
        Int_t nentries = aliasArr->GetEntries();

    	Double_t *times = new Double_t[nentries];
    	Double_t *values = new Double_t[nentries];
	UInt_t iValue=0;
	
	while((aValue = (AliDCSValue*) iterarray.Next())) {
	               UInt_t currentTime = aValue->GetTimeStamp();
	               if(currentTime>fCtpEndTime) break;

   	               values[iValue] = aValue->GetFloat();
   	               times[iValue] = (Double_t) (currentTime);
		       iValue++;
	               }
	CreateGraph(iAlias, iValue, times, values); // fill graphs
        } 
      else { // Treating FEE Parameters
		AliDCSValue * lastVal = NULL;
		while((aValue = (AliDCSValue*) iterarray.Next())) lastVal = aValue; // Take only the last value
		fFEEParameters->Add(new TObjString(fAliasNames[iAlias].Data()),lastVal);
	}      
  }
  
  fIsProcessed=kTRUE;

  return success;
}

//_____________________________________________________________________________
void AliADDataDCS::Init(){

  // initialization of aliases and DCS data

  TString sindex;
  int iAlias = 0;
  
  for(int iPM = 0; iPM<16 ; iPM++){
  		fAliasNames[iAlias] = Form("AD0/HV/PM%d",iPM);
			
		fHv[iAlias] = new TH1F(fAliasNames[iAlias].Data(),fAliasNames[iAlias].Data(), 3000, kHvMin, kHvMax);
		fHv[iAlias]->GetXaxis()->SetTitle("Hv");
		iAlias++;
  }
  
  for(int iPM = 0; iPM<16 ; iPM++) fAliasNames[iAlias++] = Form("AD0/HV/Imon%d",iPM);
  
  // CCIU Parameters
	
  fAliasNames[iAlias++] = "AD0/FEE/CCIU/BBAThreshold";
  fAliasNames[iAlias++] = "AD0/FEE/CCIU/BBCThreshold";
  fAliasNames[iAlias++] = "AD0/FEE/CCIU/BGAThreshold";
  fAliasNames[iAlias++] = "AD0/FEE/CCIU/BGCThreshold";
  fAliasNames[iAlias++] = "AD0/FEE/CCIU/BBAForBGThreshold";
  fAliasNames[iAlias++] = "AD0/FEE/CCIU/BBCForBGThreshold";
  fAliasNames[iAlias++] = "AD0/FEE/CCIU/MultADAThrLow";
  fAliasNames[iAlias++] = "AD0/FEE/CCIU/MultADAThrHigh";
  fAliasNames[iAlias++] = "AD0/FEE/CCIU/MultADCThrLow";
  fAliasNames[iAlias++] = "AD0/FEE/CCIU/MultADCThrHigh";
  for(int i=1;i<=5;i++) {
    	fAliasNames[iAlias] = "AD0/FEE/CCIU/TriggerSelect";
    	fAliasNames[iAlias++] += Form("%d",i);
  }
  //Trigger parameters
  for(int iCIU = 0; iCIU<2 ; iCIU++){
  		fAliasNames[iAlias++] = Form("AD0/FEE/CIU%d/Clk1Win1",iCIU);
		fAliasNames[iAlias++] = Form("AD0/FEE/CIU%d/Clk1Win2",iCIU);
		fAliasNames[iAlias++] = Form("AD0/FEE/CIU%d/Clk2Win1",iCIU);
		fAliasNames[iAlias++] = Form("AD0/FEE/CIU%d/Clk2Win2",iCIU);
		fAliasNames[iAlias++] = Form("AD0/FEE/CIU%d/DelayClk1Win1",iCIU);
		fAliasNames[iAlias++] = Form("AD0/FEE/CIU%d/DelayClk1Win2",iCIU);
		fAliasNames[iAlias++] = Form("AD0/FEE/CIU%d/DelayClk2Win1",iCIU);
		fAliasNames[iAlias++] = Form("AD0/FEE/CIU%d/DelayClk2Win2",iCIU);
		fAliasNames[iAlias++] = Form("AD0/FEE/CIU%d/LatchWin1",iCIU);
		fAliasNames[iAlias++] = Form("AD0/FEE/CIU%d/LatchWin2",iCIU);
		fAliasNames[iAlias++] = Form("AD0/FEE/CIU%d/ResetWin1",iCIU);
		fAliasNames[iAlias++] = Form("AD0/FEE/CIU%d/ResetWin2",iCIU);
		fAliasNames[iAlias++] = Form("AD0/FEE/CIU%d/PedestalSubtraction",iCIU);
  }
  for(int iCIU = 0; iCIU<2 ; iCIU++){
	  for(int iCh=1;iCh<=8;iCh++){
	    	fAliasNames[iAlias++] = Form("AD0/FEE/CIU%d/EnableCharge%d",iCIU,iCh);
		fAliasNames[iAlias++] = Form("AD0/FEE/CIU%d/EnableTiming%d",iCIU,iCh);
		fAliasNames[iAlias++] = Form("AD0/FEE/CIU%d/PedEven%d",iCIU,iCh);
		fAliasNames[iAlias++] = Form("AD0/FEE/CIU%d/PedOdd%d",iCIU,iCh);
		fAliasNames[iAlias++] = Form("AD0/FEE/CIU%d/PedCutEven%d",iCIU,iCh);
		fAliasNames[iAlias++] = Form("AD0/FEE/CIU%d/PedCutOdd%d",iCIU,iCh);
	  }
  }
 // Time Resolution Parameters	
  for(int iCIU = 0; iCIU<2 ; iCIU++){
		fAliasNames[iAlias++] = Form("AD0/FEE/CIU%d/TimeResolution",iCIU);
		fAliasNames[iAlias++] = Form("AD0/FEE/CIU%d/WidthResolution",iCIU);
  }
  // HPTDC parameters
  for(int iCIU = 0; iCIU<2 ; iCIU++){
	  	fAliasNames[iAlias++] = Form("AD0/FEE/CIU%d/MatchWindow",iCIU);
	  	fAliasNames[iAlias++] = Form("AD0/FEE/CIU%d/SearchWindow",iCIU);
	  	fAliasNames[iAlias++] = Form("AD0/FEE/CIU%d/TriggerCountOffset",iCIU);
	  	fAliasNames[iAlias++] = Form("AD0/FEE/CIU%d/RollOver",iCIU);
  }

  for(int iCIU = 0; iCIU<2 ; iCIU++){
	  for(int iCh=1;iCh<=8;iCh++){
	    	fAliasNames[iAlias++] = Form("AD0/FEE/CIU%d/DelayHit%d",iCIU,iCh);
	  }
  }

  for(int iCIU = 0; iCIU<2 ; iCIU++){
	  for(int iCh=1;iCh<=8;iCh++){
	    	fAliasNames[iAlias++] = Form("AD0/FEE/CIU%d/DiscriThr%d",iCIU,iCh);
	  }
  }

  if(iAlias!=kNAliases) 
  	      AliError(Form("Number of DCS Aliases defined not correct"));

}

//_____________________________________________________________________________
void AliADDataDCS::Introduce(UInt_t numAlias, const TObjArray* aliasArr)const
{

  // method to introduce new aliases

  int entries=aliasArr->GetEntries();
  AliInfo(Form("************ Alias: %s **********",fAliasNames[numAlias].Data()));
  AliInfo(Form("    	%d DP values collected",entries));

}

//_____________________________________________________________________________
void AliADDataDCS::CreateGraph(int i, int dim, const Double_t *x, const Double_t *y)
{

   // Create graphics
   
   TGraph *gr = new((*fGraphs)[fGraphs->GetEntriesFast()]) TGraph(dim, x, y);

   gr->GetXaxis()->SetTimeDisplay(1);
   gr->SetTitle(fAliasNames[i].Data());

}


//_____________________________________________________________________________
void AliADDataDCS::Draw(const Option_t* /*option*/)
{
// Draw all histos and graphs

  if(!fIsProcessed) return;

  if(fGraphs->GetEntries()==0)  return;
  
  TCanvas *cHV = new TCanvas("AD0_HV","AD0_HV");
  cHV->Divide(4,4);
  
  for(int iPM = 0; iPM<16 ; iPM++){
  	cHV->cd(iPM+1);
  	((TGraph*) fGraphs->UncheckedAt(iPM))->SetMarkerStyle(20);
  	((TGraph*) fGraphs->UncheckedAt(iPM))->Draw("ALP");	  		
  	}
}

