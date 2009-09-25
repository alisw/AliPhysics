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
// Class AliVZEROTrending
// ---------------------------
// 
//  class used in QA to publish variables evolution versus time in AMORE. 
//  These histo are the one which will be looked at by QA Shifter
// 
#include "TGraph.h"

#include "AliLog.h"
#include "AliVZEROTrending.h"

ClassImp(AliVZEROTrending)

//_____________________________________________________________________________
AliVZEROTrending::AliVZEROTrending() : TH1(), fNEntries(0)
{
	// Default constructor
}
//_____________________________________________________________________________
AliVZEROTrending::AliVZEROTrending(const char* name, const char* title) : TH1(), fNEntries(0)
{
	SetName(name);
	SetTitle(title);
}
//_____________________________________________________________________________
AliVZEROTrending::AliVZEROTrending(const AliVZEROTrending &trend) : 
	TH1(), fNEntries(trend.fNEntries)
{
	// Copy constructor
	SetName(trend.GetName());
	SetTitle(trend.GetTitle());
	for (int i = 0; i < kDataSize; i++) {
		fTime[i] = trend.fTime[i];
		for (int j = 0; j < 8; j++) {
			fData[j][i] = trend.fData[j][i];
		}
	}
}

//_____________________________________________________________________________
AliVZEROTrending::~AliVZEROTrending(){
	// Destructor
}
// -----------------------------------------------------------------			
void AliVZEROTrending::AddEntry(Double_t * data, UInt_t time)
{

	if(fNEntries<kDataSize){
		for (int i = 0; i < 8; i++)
		{
			fData[i][fNEntries] = data[i];
			fTime[fNEntries] = (double) time;
		}
		fNEntries++;	
	}else{

		for (int i = 0; i < kDataSize-1; i++){
			fTime[i] = fTime[i+1];
			for (int ich = 0; ich < 8; ich++){		
				fData[ich][i] = fData[ich][i+1];
			}	
		}
		for (int i = 0; i < 8; i++)
		{
			fData[i][fNEntries-1] = data[i];
			fTime[fNEntries-1] = (double) time;
		}
		
	}
	//printf("sizeof UInt_t Double_t %d %d\n",sizeof(UInt_t),sizeof(Double_t));
	//printf("Add Entry %d @ %f : %f %f %f %f %f %f %f %f \n",fNEntries,fTime[fNEntries-1], 
	//	data[0],data[1],data[2],data[3],data[4],data[5],data[6],data[7]);
}			
// -----------------------------------------------------------------			
void AliVZEROTrending::PrintEntry(UInt_t entry)
{

	if(entry>=fNEntries){
		AliError(Form("maximum entry is %d\n",fNEntries-1));
	}else{
		AliInfo(Form("Entry %d @ %f : %f %f %f %f %f %f %f %f \n",entry, fTime[entry],
			fData[0][entry],fData[1][entry],fData[2][entry],fData[3][entry],fData[4][entry],fData[5][entry],fData[6][entry],fData[7][entry]));

	}
}			
