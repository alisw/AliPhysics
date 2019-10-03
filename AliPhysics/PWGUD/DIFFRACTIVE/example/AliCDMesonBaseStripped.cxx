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
// AliCDMesonBaseStripped
//
//  Author:
//  Felix Reidt <Felix.Reidt@cern.ch>

#include "TH1F.h"
#include "TString.h"

#include "AliCDMesonBaseStripped.h"


//------------------------------------------------------------------------------
Int_t AliCDMesonBaseStripped::GetGapBin(TString tag, Int_t gapcg)
{
	//
	// retrieve gap topology for a given string
	//

	tag.ToUpper();

	Bool_t ka = kFALSE, kc = kFALSE;

	if(tag.Contains("V0")){
		ka = ka || (gapcg & kBitV0A);
		kc = kc || (gapcg & kBitV0C);
	}
	if(tag.Contains("FMD")){
		ka = ka || (gapcg & kBitFMDA);
		kc = kc || (gapcg & kBitFMDC);
	}
	if(tag.Contains("SPD")){
		ka = ka || (gapcg & kBitSPDA);
		kc = kc || (gapcg & kBitSPDC);
	}
	if(tag.Contains("TPC")){
		ka = ka || (gapcg & kBitTPCA);
		kc = kc || (gapcg & kBitTPCC);
	}
	if(ka && kc)
		return kBinNG;
	else{
		if(!ka && !kc)
			return kBinDG;
		else if(!kc)
			return kBinGC;
		else
			return kBinGA;
	}
}


//------------------------------------------------------------------------------
TH1F* AliCDMesonBaseStripped::GetHistStatsFlow()
{
	//
	// setup the stats flow histogram
	//

	TH1F *hist = new TH1F("c00_statsFlow", "",
	                      AliCDMesonBaseStripped::kBinLastValue,
	                      0, AliCDMesonBaseStripped::kBinLastValue);
	TAxis* axis = hist->GetXaxis();
	axis->SetBinLabel(AliCDMesonBaseStripped::kBinTotalInput+1, "total Input");
	axis->SetBinLabel(AliCDMesonBaseStripped::kBinGoodInput+1, "good ESDs");
	axis->SetBinLabel(AliCDMesonBaseStripped::kBinEventsAfterCuts+1,"after cuts");
	axis->SetBinLabel(AliCDMesonBaseStripped::kBinEventsWithOutPileUp+1,
	                  "w/o pile up");
	axis->SetBinLabel(AliCDMesonBaseStripped::kBinv0Gap+1, "with V0 DG gap");
	axis->SetBinLabel(AliCDMesonBaseStripped::kBinv0fmdGap+1,
	                  "with V0-FMD DG gap");
	axis->SetBinLabel(AliCDMesonBaseStripped::kBinv0fmdspdGap+1,
	                  "with V0-FMD-SPD DG gap");
	axis->SetBinLabel(AliCDMesonBaseStripped::kBinv0fmdspdtpcGap+1,
	                  "with V0-FMD-SPD-TPC DG gap");
	axis->SetBinLabel(AliCDMesonBaseStripped::kBinResidualTracks+1,
	                  "without residual tracks");
	axis->SetBinLabel(AliCDMesonBaseStripped::kBinResidualTracklets+1,
	                  "without residual tracklets");
	return hist;
}
