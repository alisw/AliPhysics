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

// Add task B2
// author: Eulogio Serradilla <eulogio.serradilla@cern.ch>

AliAnalysisTaskB2* AddTaskB2(  const TString& species
                             , const TString& containername
                             , const TString& trksel
                             , Int_t pidProc
                             , const TString& periodname
                             , Bool_t   simulation       = kFALSE
                             , Bool_t   heavyIons        = kFALSE
                             , Double_t maxDCAxy         = 1
                             , Double_t maxDCAz          = 2
                             , Double_t maxEta           = 0.8
                             , Double_t maxY             = 0.5
                             , Bool_t ntrkMultTrigger    = 0
                             , Double_t minKNOmult       = -10
                             , Double_t maxKNOmult       = 10000
                             , Bool_t   V0AND            = kFALSE
                             , const TString& ztag       = ""
                             , Double_t maxVz            = 10
                             , Bool_t momentumCorr       = kFALSE
                             , const TString& binSize    = ""
                             , Bool_t xRowsTPC           = 0
                             , Int_t minTPCnClsOrXRows   = 70
                             , Double_t minCentrality    = 0
                             , Double_t maxCentrality    = 20
                             , Double_t minM2            = 2.
                             , Double_t maxM2            = 6.)
{
//
// Create, configure and add the analysis task to the analysis manager
//
	using namespace std;
	
	// sample config
	
	const Double_t kMaxVx     = 1.;
	const Double_t kMaxVy     = 1.;
	
	const Double_t kMaxNSigma = 3.;
	
	const Int_t kMaxNSigmaITS = 3;
	const Int_t kMaxNSigmaTPC = 3;
	const Int_t kMaxNSigmaTOF = 3;
	
	TString period = periodname;
	period.ToLower();
	
	// Get pointer to the existing analysis manager
	
	AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
	if (!mgr)
	{
		cerr << "AddTaskB2: no analysis manager to connect to" << endl;
		return 0;
	}
	
	if (!mgr->GetInputEventHandler())
	{
		cerr << "AddTaskB2: this task requires an input event handler" << endl;
		return 0;
	}
	
	// Create and configure the task
	
	AliAnalysisTaskB2* task = new AliAnalysisTaskB2(Form("B2.%s",containername.Data()));
	
	task->SetParticleSpecies(species);
	task->SetSimulation(simulation);
	task->SetHeavyIons(heavyIons);
	task->SetV0ANDtrigger(V0AND);
	
	task->SetMaxNSigmaITS(kMaxNSigmaITS);
	task->SetMaxNSigmaTPC(kMaxNSigmaTPC);
	task->SetMaxNSigmaTOF(kMaxNSigmaTOF);
	
	task->SetNtrkMultTrigger(ntrkMultTrigger);
	
	gROOT->LoadMacro("$ALICE_ROOT/PWGLF/SPECTRA/Nuclei/B2/macros/MeanNtrk.C");
	Double_t meanNtrk = MeanNtrk(period, maxEta, V0AND);
	
	task->SetMeanNtrk(meanNtrk);
	task->SetKNOmultInterval(minKNOmult, maxKNOmult);
	task->SetVertexXInterval(-kMaxVx, kMaxVx);
	task->SetVertexYInterval(-kMaxVy, kMaxVy);
	task->SetVertexZInterval(-maxVz, maxVz);
	
	task->SetEtaInterval(-maxEta, maxEta);
	task->SetRapidityInterval(-maxY, maxY);
	task->SetM2Interval(minM2, maxM2);
	
	task->SetCentralityInterval(minCentrality, maxCentrality);
	
	if(period=="lhc11a_wsdd" || period=="lhc11a_wosdd")
	{
		task->SetNoFastOnlyTrigger();
	}
	
	// momentum correction
	
	if(momentumCorr && species=="Deuteron")
	{
		gROOT->LoadMacro("$ALICE_ROOT/PWGLF/SPECTRA/Nuclei/B2/macros/MomentumCorrection.C");
		
		TProfile* pfx = MomentumCorrection(species);
		task->SetMomentumCorrectionProfile(pfx);
		if(pfx != 0) task->SetMomentumCorrection();
	}
	
	// histograms
	
	gROOT->LoadMacro("$ALICE_ROOT/PWGLF/SPECTRA/Nuclei/B2/macros/CreateHistograms.C");
	
	AliLnHistoMap* hMap = CreateHistograms(species, binSize, simulation, maxDCAxy, maxEta, maxY, heavyIons);
	
	task->SetHistogramMap(hMap);
	
	// track selection criteria
	
	gROOT->LoadMacro("$ALICE_ROOT/PWGLF/SPECTRA/Nuclei/B2/macros/TrackCuts.C");
	
	AliESDtrackCuts* trkCuts = TrackCuts(task, trksel, maxDCAxy, maxDCAz, kMaxNSigma, xRowsTPC, minTPCnClsOrXRows, maxEta);
	task->SetESDtrackCuts(trkCuts);
	
	// PID
	
	gROOT->LoadMacro("$ALICE_ROOT/PWGLF/SPECTRA/Nuclei/B2/macros/BetheBlochParams.C");
	
	Double_t bethe[5];
	BetheBlochParams(bethe, period);
	
	gROOT->LoadMacro("$ALICE_ROOT/PWGLF/SPECTRA/Nuclei/B2/macros/PriorProbabilities.C");
	
	Double_t prob[9];
	PriorProbabilities(prob, period, trksel, ztag);
	
	AliLnID* lnID = new AliLnID();
	
	lnID->SetTPCBetheBlochParams(bethe);
	lnID->SetPriorProbabilities(prob);
	lnID->SetPidProcedure(pidProc);
	
	if(!simulation) lnID->SetTPCChargeCorrection(2.3);
	
	task->SetPID(lnID);
	
	// Add task to the manager
	
	mgr->AddTask(task);
	
	// input and output containers
	
	AliAnalysisDataContainer* output = mgr->CreateContainer(containername.Data(), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:AliAnalysisTaskB2", AliAnalysisManager::GetCommonFileName()));
	
	mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
	mgr->ConnectOutput(task, 1, output);
	
	return task;
}
