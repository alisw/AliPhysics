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

// macro for a predefined set of track cuts
// author: Eulogio Serradilla <eulogio.serradilla@cern.ch>

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <Riostream.h>
#include <AliESDtrackCuts.h>
#include <TString.h>
#include "AliAnalysisTaskB2.h"
#endif

void AddNSigmaCuts(AliESDtrackCuts* trkCuts, Float_t maxNSigma=3)
{
//
// nsigma cuts
//
	trkCuts->SetRequireSigmaToVertex(kTRUE);
	trkCuts->SetMaxNsigmaToVertex(maxNSigma);
}

void AddDCACuts(AliESDtrackCuts* trkCuts, Float_t maxDCAxy=1.5, Float_t maxDCAz=2)
{
//
// DCA cuts
//
	trkCuts->SetMaxDCAToVertexXY(maxDCAxy);
	trkCuts->SetMaxDCAToVertexZ(maxDCAz);
}

AliESDtrackCuts* TrackCuts(AliAnalysisTaskB2* task, const TString& trksel, Double_t maxDCAxy, Double_t maxDCAz, Double_t maxNSigma, Bool_t xrows, Int_t minTPCnClsOrXRows, Double_t maxEta)
{
//
// Create an AliESDtrackCuts from a predefined set
//
	AliESDtrackCuts* trkCuts = new AliESDtrackCuts("AliESDtrackCuts");
	
	trkCuts->SetEtaRange(-maxEta, maxEta);
	trkCuts->SetPtRange(0.15, 100.);
	trkCuts->SetAcceptKinkDaughters(kFALSE);
	
	// ITS
	trkCuts->SetRequireITSRefit(kTRUE);
	trkCuts->SetMinNClustersITS(2);
	trkCuts->SetMaxChi2PerClusterITS(36);
	
	// TPC
	trkCuts->SetRequireTPCRefit(kTRUE);
	if(xrows)
	{
		trkCuts->SetMinNCrossedRowsTPC(minTPCnClsOrXRows);
		trkCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
	}
	else
	{
		trkCuts->SetMinNClustersTPC(minTPCnClsOrXRows);
	}
	trkCuts->SetMaxChi2PerClusterTPC(4.);
	//trkCuts->SetMaxChi2TPCConstrainedGlobal(36);
	
	TString tracksel = trksel;
	tracksel.ToLower();
	
	Int_t clusterCut = (xrows) ? 1 : 0;
	
	if(tracksel == "its_tpc_nsigma")
	{
		AddNSigmaCuts(trkCuts, maxNSigma);
	}
	else if(tracksel == "its_tpc_nsigma_spd1")
	{
		AddNSigmaCuts(trkCuts, maxNSigma);
		trkCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kFirst);
	}
	else if(tracksel == "its_tpc_nsigma_spd")
	{
		AddNSigmaCuts(trkCuts, maxNSigma);
		trkCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kAny);
	}
	else if(tracksel == "its_tpc_dca")
	{
		AddDCACuts(trkCuts,maxDCAxy,maxDCAz);
	}
	else if(tracksel == "its_tpc_dca_spd1")
	{
		AddDCACuts(trkCuts,maxDCAxy,maxDCAz);
		trkCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kFirst);
	}
	else if(tracksel == "its_tpc_dca_spd")
	{
		AddDCACuts(trkCuts,maxDCAxy,maxDCAz);
		trkCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kAny);
	}
	else if(tracksel == "its_tpc_tof_nsigma")
	{
		AddNSigmaCuts(trkCuts, maxNSigma);
		task->SetTOFmatch(1);
	}
	else if(tracksel == "its_tpc_tof_nsigma_spd1")
	{
		AddNSigmaCuts(trkCuts, maxNSigma);
		task->SetTOFmatch(1);
		trkCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kFirst);
	}
	else if(tracksel == "its_tpc_tof_nsigma_spd")
	{
		AddNSigmaCuts(trkCuts, maxNSigma);
		task->SetTOFmatch(1);
		trkCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kAny);
	}
	else if(tracksel == "its_tpc_tof_dca")
	{
		AddDCACuts(trkCuts,maxDCAxy,maxDCAz);
		task->SetTOFmatch(1);
	}
	else if(tracksel == "its_tpc_tof_dca_spd1")
	{
		AddDCACuts(trkCuts,maxDCAxy,maxDCAz);
		task->SetTOFmatch(1);
		trkCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kFirst);
	}
	else if(tracksel == "its_tpc_tof_dca_spd")
	{
		AddDCACuts(trkCuts,maxDCAxy,maxDCAz);
		task->SetTOFmatch(1);
		trkCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kAny);
	}
	else if(tracksel == "std_its_tpc_2009") // pp data 2009
	{
		delete trkCuts;
		return AliESDtrackCuts::GetStandardITSTPCTrackCuts2009(kFALSE);
	}
	else if(tracksel == "std_its_tpc_2010") // pp data 2010
	{
		delete trkCuts;
		return AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(kFALSE, clusterCut);
	}
	else if(tracksel == "std_its_tpc_2011") // pp data 2011
	{
		delete trkCuts;
		return AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kFALSE, clusterCut);
	}
	else
	{
		std::cerr << "Warning: no track selection criteria selected" << std::endl;
	}
	
	return trkCuts;
}
