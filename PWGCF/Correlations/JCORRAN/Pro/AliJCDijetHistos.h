/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */
//______________________________________________________________________________
// Basic histogram implimentation via AliJHistogramInterface.
// author: O. Saarimaki, D.J. Kim (dong.jo.kim@cern.ch)
// ALICE Group University of Jyvaskyla, Finland
// Note: Adapted for AliAnalysisTaskSE
//////////////////////////////////////////////////////////////////////////////
#ifndef ALIJCDIJETHISTOS_H
#define ALIJCDIJETHISTOS_H

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

class AliJCDijetHistos : public AliJHistogramInterface
{

public:

	AliJCDijetHistos(); //constructor
	AliJCDijetHistos(const AliJCDijetHistos& obj); // copy constructor
	virtual ~AliJCDijetHistos();    //destructor
	AliJCDijetHistos& operator=(const AliJCDijetHistos& obj); // equal sign operator

	void CreateEventTrackHistos();
	#define NCENT 8
	static Double_t CentBin[NCENT+1]; //8
	static Double_t pttJacek[74+16];
	static UInt_t NCentBin;
	static UInt_t NpttJacek;

	static int GetCentralityClass(Double_t);

	AliJHistManager * fHMG;//!
	AliJBin fHistCentBin;//! cent bin
	//===================================================
	// Event/Track histograms
	//===================================================
	AliJTH1D fh_pt;     //! // for pt dist of tracks
	AliJTH1D fh_eta;    //! // for eta dist of tracks
	AliJTH2D fh_etaPhi; //! // for (eta,phi) dist of tracks
	AliJTH1D fh_phi;    //! // for phi dist of tracks

	AliJTH1D fh_jetPt;      //! // for pt dist of jets
	AliJTH1D fh_jetEta;     //! // for eta dist of jets
	AliJTH1D fh_jetPhi;     //! // for phi dist of jets
	AliJTH2D fh_jetEtaPhi;  //! // for (eta,phi) dist of jets
	AliJTH1D fh_rho;        //! // for event pt density
	AliJTH1D fh_rhom;       //! // for event mt density
	AliJTH1D fh_jetArea;    //! // for jet area spectrum
	AliJTH1D fh_jetAreaRho; //! // for jet area*pt spectrum
	AliJTH1D fh_corrJetPt;  //! // for pt dist of bg subtracted jets
	AliJTH1D fh_corrJetEta; //! // for eta dist of bg subtracted jets
	AliJTH1D fh_corrJetPhi; //! // for phi dist of bg subtracted jets

    AliJTH1D fh_dijetInvM;                //! // for dijet invariant mass
    AliJTH1D fh_dijetPtPair;              //! // for dijet pt
    AliJTH1D fh_dijetDeltaPhi;            //! // for dijet deltaPhi
    AliJTH1D fh_dijetInvMDeltaPhiCut;     //! // for dijet invariant mass after deltaPhi cut
    AliJTH1D fh_corrDijetInvM;            //! // for bg subtracted dijet invariant mass
    AliJTH1D fh_corrDijetPtPair;          //! // for bg subtracted dijet pt
    AliJTH1D fh_corrDijetDeltaPhi;        //! // for bg subtracted dijet deltaPhi
    AliJTH1D fh_corrDijetInvMDeltaPhiCut; //! // for bg subtracted dijet invariant mass after deltaPhi cut


};

#endif //ALIJCDIJETHISTOS_H
