/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */
//______________________________________________________________________________
// Basic histogram implimentation via AliJHistogramInterface.
// author:  D.J. Kim (dong.jo.kim@cern.ch)
// ALICE Group University of Jyvaskyla, Finland
// Note: Adapted for AliAnalysisTaskSE
//////////////////////////////////////////////////////////////////////////////
#ifndef ALIJFLOWHISTOS_H
#define ALIJFLOWHISTOS_H

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>

#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>
#include <TFile.h>
#include <TList.h>
#include <TLorentzVector.h>

#include <AliJHistogramInterface.h>


using namespace std;

class AliJFlowHistos : public AliJHistogramInterface
{

public:

	AliJFlowHistos(); //constructor
	AliJFlowHistos(const AliJFlowHistos& obj); // copy constructor
	virtual ~AliJFlowHistos();    //destructor
	AliJFlowHistos& operator=(const AliJFlowHistos& obj); // equal sign operator

	void CreateEventTrackHistos();
	#define NCENT 16 
	static Double_t CentBin[NCENT+1]; //8
	static Double_t pttJacek[74];
	static UInt_t NCentBin;
	static UInt_t NpttJacek;

	static int GetCentralityClass(Double_t);

	AliJHistManager * fHMG;//!
	AliJBin fHistCentBin;//! cent bin
	AliJBin fBin_DetSet;//! detector
	AliJBin fBin_hh;//! harmonics
	//===================================================
	// Event/Track histograms
	//===================================================
	AliJTH1D fh_pt;//! // for pt dist of tracks
	AliJTH1D fh_eta;//! // for eta dist of tracks
	AliJTH1D fh_phi;//! // for phi dist [ic][isub]
	AliJTH1D fh_EP;//! // for Q-Vector dist [ic][isub][ih]
	AliJTH1D fhEPCorrInHar;//! for 3subevt resolution [ic][isub][ih]
	//AliJTH2D fhEPCorr2D;//! EP correlations


};

#endif //ALIJFLOWHISTOS_H






















