/* Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice */

// Container class for histograms needed in the jT analysis.

//===========================================================
// AliJFlowHistos.h
//
// author: Marton Vargyas
//===========================================================

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

#include "AliJHistogramInterface.h"
#include "AliJFlowBaseTask.h" 


using namespace std;

class AliJFlowHistos : public AliJHistogramInterface
{

public:

	AliJFlowHistos(); //constructor
	AliJFlowHistos(const AliJFlowHistos& obj); // copy constructor
	virtual ~AliJFlowHistos();    //destructor
	AliJFlowHistos& operator=(const AliJFlowHistos& obj); // equal sign operator

	void CreateEventTrackHistos();
	#define CENTN 8
	static Double_t CentBin[CENTN+1]; //8
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


};

#endif //ALIJFLOWHISTOS_H






















