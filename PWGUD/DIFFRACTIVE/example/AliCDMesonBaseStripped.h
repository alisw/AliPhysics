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

#ifndef ALICDMESONBASESTRIPPED_H
#define ALICDMESONBASESTRIPPED_H

class TH1F;

class AliCDMesonBaseStripped
{
public:
	enum{
		// gap Conditions
		kBinDG = 1, // double gap
		kBinGC, // single gap c side
		kBinGA, // single gap a side
		kBinNG, // no gap

		// StatsFlow histogram entries
		// names for the bins are specified in the .cxx-file
		kBinTotalInput = 0,
		kBinGoodInput,
		kBinEventsAfterCuts,
		kBinEventsWithOutPileUp,
		kBinv0Gap,
		kBinv0fmdGap,
		kBinv0fmdspdGap,
		kBinv0fmdspdtpcGap,
		kBinResidualTracks,
		kBinResidualTracklets,
		kBinLastValue, // used to specify the correct histogram width


		// gap-condition bits used in AliAnalysisTaskCDMeson::fGapRun
		kBitBaseLine = (1<<0),

		kBitV0A  = (1<<1),
		kBitV0C  = (1<<2),
		kBitFMDA = (1<<3),
		kBitFMDC = (1<<4),

		kBitSPDA  = (1<<5),
		kBitSPDC  = (1<<6),
		kBitTPCA  = (1<<7),
		kBitTPCC  = (1<<8),
	};

	static Int_t GetGapBin(TString tag, Int_t gapcg);
	static TH1F* GetHistStatsFlow();
};

#endif
