/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *  * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/*
$Log: AliT0DataDCS.cxx,v $
 Revision   2008/01/30 
Fetching data points from DCS, calculating mean and storing data to Reference DB 
 
 Version 1.1  2006/10
Preliminary test version (T.Malkiewicz)
*/

#include "AliT0DataDCS.h"

#include "AliCDBMetaData.h"
#include "AliDCSValue.h"
#include "AliLog.h"

#include <TTimeStamp.h>
#include <TObjString.h>
#include <TH2F.h>
#include <TProfile.h>
#include <TGraph.h>
#include <TDatime.h>
#include <TStyle.h>
#include <TCanvas.h>

// AliT0DataDCS class
// declaring DCS aliases for T0
// fetching T0 data points from DCS, 
// calculating mean values for the entire run
// and storing the result to Reference DB

ClassImp(AliT0DataDCS)

//---------------------------------------------------------------
AliT0DataDCS::AliT0DataDCS():
	TObject(),
	fRun(0),	
	fStartTime(0),
	fEndTime(0),
	fIsProcessed(kFALSE)
{
}

//---------------------------------------------------------------
AliT0DataDCS::AliT0DataDCS(Int_t nRun, UInt_t startTime, UInt_t endTime):
	TObject(),
	fRun(nRun),
	fStartTime(startTime),
	fEndTime(endTime),
	fIsProcessed(kFALSE)
{
	AliInfo(Form("\n\tRun %d \n\tStartTime %s \n\tEndTime %s", nRun,
	TTimeStamp(startTime).AsString(),
	TTimeStamp(endTime).AsString()));

	Init();

}

//---------------------------------------------------------------
AliT0DataDCS::~AliT0DataDCS() {

}

//---------------------------------------------------------------
Bool_t AliT0DataDCS::ProcessData(TMap& aliasMap)
{
 	        Float_t t0_scaler[32];
		Int_t aliasEntr[32];
		TObjArray *aliasArr;
        	for(int j=0; j<kNAliases; j++)
        	{
		  for (Int_t k=0;k<32;k++) 
		  {
		    t0_scaler[k]=0;
		    aliasEntr[k]=0;	
      		  }
		  aliasArr = (TObjArray*) aliasMap.GetValue(fAliasNames[j].Data());
                  if(!aliasArr)
                  {
                        AliError(Form("Alias %s not found!", fAliasNames[j].Data()));
                        continue;
                  }
                  if(aliasArr->GetEntries()<2)
                  {
                        AliError(Form("Alias %s has just %d entries!",
                                        fAliasNames[j].Data(),aliasArr->GetEntries()));
                        continue;
                  }
		  aliasEntr[j] = aliasArr->GetEntries();
		  for(int l=0; l<aliasEntr[j]; l++)
		  {
                  AliDCSValue *aValue=dynamic_cast<AliDCSValue*> (aliasArr->At(l));
		  t0_scaler[j]+= aValue->GetFloat();
	          }
		fScalerMean[j] = t0_scaler[j] / aliasEntr[j] ;
	        }
	fIsProcessed=kTRUE;
	return kTRUE;
}

//---------------------------------------------------------------
void AliT0DataDCS::Init()
{
	TString sindex;
	for(int i=0;i<kNAliases;i++)
	{
		fAliasNames[i] = "t00_ac_scaler_";
		 sindex.Form("%02d",i);
	        fAliasNames[i] += sindex;
	}

}

//---------------------------------------------------------------
void AliT0DataDCS::Introduce(UInt_t numAlias, const TObjArray* aliasArr){

	int entries=aliasArr->GetEntries();
	AliInfo(Form("************ Alias: %s **********",fAliasNames[numAlias].Data()));
	AliInfo(Form("    	%d DP values collected",entries));

}



