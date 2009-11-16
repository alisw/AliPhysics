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


#include "AliVZERODataDCS.h"

#include "AliDCSValue.h"
#include "AliLog.h"

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

// AliVZERODataDCS class
// main aim to introduce the aliases for the VZERO DCS
// data points to be then
// stored in the OCDB, and to process them. 
// ProcessData() method called by VZEROPreprocessor

ClassImp(AliVZERODataDCS)

//_____________________________________________________________________________
AliVZERODataDCS::AliVZERODataDCS():
	TObject(),
	fRun(0),
	fStartTime(0),
	fEndTime(0),
    fGraphs("TGraph",kNGraphs),
	fFEEParameters(NULL),
	fIsProcessed(kFALSE)

{
  // Default constructor
	for(int i=0;i<kNHvChannel;i++) {
		fDeadChannel[i] = kFALSE;
		fMeanHV[i]      = 100.0;
		fWidthHV[i]     = 0.0; 
	}
}

//_____________________________________________________________________________
AliVZERODataDCS::AliVZERODataDCS(Int_t nRun, UInt_t startTime, UInt_t endTime):
	TObject(),
	fRun(nRun),
	fStartTime(startTime),
	fEndTime(endTime),
	fGraphs("TGraph",kNGraphs),
	fFEEParameters(new TMap()),
	fIsProcessed(kFALSE)

{

  // constructor with arguments
  	for(int i=0;i<kNHvChannel;i++) {
	 fDeadChannel[i] = kFALSE;        
	 fMeanHV[i]      = 100.0;
     fWidthHV[i]     = 0.0; 
	}
	AliInfo(Form("\n\tRun %d \n\tStartTime %s \n\tEndTime %s", nRun,
		TTimeStamp(startTime).AsString(),
		TTimeStamp(endTime).AsString()));
	
	fFEEParameters->SetOwnerValue();
	Init();

}

//_____________________________________________________________________________
AliVZERODataDCS::~AliVZERODataDCS() {

  // destructor
  fGraphs.Clear("C");
  delete fFEEParameters;

}

//_____________________________________________________________________________
void AliVZERODataDCS::ProcessData(TMap& aliasMap){

  // method to process the data

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

    //Introduce(iAlias, aliasArr);
    
    if(aliasArr->GetEntries()<2){
      AliWarning(Form("Alias %s has just %d entries!",
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
   			values[iValue] = aValue->GetFloat();
			if(iValue>0) {
				if(values[iValue-1]>0.) variation = TMath::Abs(values[iValue]-values[iValue-1])/values[iValue-1];
				if(variation > 0.10) fDeadChannel[GetOfflineChannel(iAlias)] = kTRUE;
			}
   			times[iValue] = (Double_t) (aValue->GetTimeStamp());
			fHv[iAlias]->Fill(values[iValue]);
			printf("%s %f Dead=%d\n",fAliasNames[iAlias].Data(),values[iValue],fDeadChannel[GetOfflineChannel(iAlias)]);
   			iValue++;
    	}      
    	CreateGraph(iAlias, aliasArr->GetEntries(), times, values); // fill graphs 

    	delete[] values;
    	delete[] times;	
	} else { // Treating FEE Parameters
		AliDCSValue * lastVal = NULL;
		while((aValue = (AliDCSValue*) iterarray.Next())) lastVal = aValue; // Take only the last value
		fFEEParameters->Add(new TObjString(fAliasNames[iAlias].Data()),lastVal);
	}      
  }
  
  	// calculate mean and rms of the first two histos
	// and convert index to aliroot channel
	for(int i=0;i<kNHvChannel;i++){
	    Int_t iChannel     = GetOfflineChannel(i);	
		fMeanHV[iChannel]  = fHv[i]->GetMean();
		fWidthHV[iChannel] = fHv[i]->GetRMS();
	}
    
  fIsProcessed=kTRUE;
}

//_____________________________________________________________________________
void AliVZERODataDCS::Init(){

  // initialization of aliases and DCS data

  TString sindex;
  int iAlias = 0;
  
  for(int iSide = 0; iSide<2 ; iSide++){
  	for(int iRing = 0; iRing<4 ; iRing++){
  		for(int iSector = 0; iSector<8 ; iSector++){
  			if(iSide == 0) fAliasNames[iAlias] = "V00/HV/V0A/SECTOR";
  			else fAliasNames[iAlias] = "V00/HV/V0C/SECTOR";
			sindex.Form("%d/RING%d",iSector,iRing);
			fAliasNames[iAlias] += sindex;
			
			fHv[iAlias] = new TH1F(fAliasNames[iAlias].Data(),fAliasNames[iAlias].Data(), 3000, kHvMin, kHvMax);
			fHv[iAlias]->GetXaxis()->SetTitle("Hv");
			iAlias++;
  		}
  	}
  }

 // Time Resolution Parameters
	
	for(int iCIU = 0; iCIU<8 ; iCIU++){
		fAliasNames[iAlias++] = Form("V00/FEE/CIU%d/TimeResolution",iCIU);
		fAliasNames[iAlias++] = Form("V00/FEE/CIU%d/WidthResolution",iCIU);
	}

  if(iAlias!=kNAliases) 
  	      AliError(Form("Number of DCS Aliases defined not correct"));

}

//_____________________________________________________________________________
void AliVZERODataDCS::Introduce(UInt_t numAlias, const TObjArray* aliasArr)const
{

  // method to introduce new aliases

  int entries=aliasArr->GetEntries();
  AliInfo(Form("************ Alias: %s **********",fAliasNames[numAlias].Data()));
  AliInfo(Form("    	%d DP values collected",entries));

}

//_____________________________________________________________________________
void AliVZERODataDCS::CreateGraph(int i, int dim, const Double_t *x, const Double_t *y)
{

   // Create graphics
   
   TGraph *gr = new(fGraphs[fGraphs.GetEntriesFast()]) TGraph(dim, x, y);

   gr->GetXaxis()->SetTimeDisplay(1);
   gr->SetTitle(fAliasNames[i].Data());

   AliInfo(Form("Array entries: %d",fGraphs.GetEntriesFast()));

}


//_____________________________________________________________________________
void AliVZERODataDCS::Draw(const Option_t* /*option*/)
{
// Draw all histos and graphs

  if(!fIsProcessed) return;

  if(fGraphs.GetEntries()==0)  return;
  
  TString canvasName;
  TCanvas *cHV[8];
  
  for(int iSide = 0 ;iSide<2;iSide++){
  	for(int iRing=0;iRing<4;iRing++){
  		if(iSide == 0)  canvasName = "V0A_Ring";
  		else  canvasName = "V0C_Ring";
  		canvasName += iRing;
  		int iCanvas = iSide*4 + iRing;
  		cHV[iCanvas] = new TCanvas(canvasName,canvasName);
  		cHV[iCanvas]->Divide(3,3);
  		for(int iSector=0;iSector<8;iSector++){
  			cHV[iCanvas]->cd(iSector+1);
  			int iChannel = iSide*32 + iRing*8 + iSector; 
  			((TGraph*) fGraphs.UncheckedAt(iChannel))->SetMarkerStyle(20);
  			((TGraph*) fGraphs.UncheckedAt(iChannel))->Draw("ALP");

  		}
  		  		
  	}
  }

}

