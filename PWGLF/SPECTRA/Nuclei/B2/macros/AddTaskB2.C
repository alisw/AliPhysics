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

Double_t GetMeanNtrk(const TString& period);

AliAnalysisTaskB2* AddTaskB2(const TString& species,
                             const TString& outputfile,
                             const TString& trksel,
                             Int_t pidProc,
                             const TString& periodname,
                             Bool_t simulation=kFALSE,
                             Bool_t heavyIons=kFALSE,
                             Double_t maxDCAxy=1,
                             Double_t maxDCAz=2,
                             Double_t minKNOmult=-10,
                             Double_t maxKNOmult=10000,
                             Bool_t V0AND=kFALSE,
                             Double_t minCentrality=0,
                             Double_t maxCentrality=20)
{
//
// Create, configure and add the analysis task to the analysis manager
//
	using namespace std;
	
	// sample config
	
	const Double_t kMaxVx     = 1.;
	const Double_t kMaxVy     = 1.;
	const Double_t kMaxVz     = 10.;
	
	const Double_t kMaxY      = 0.5;
	const Double_t kMaxEta    = 0.8;
	const Double_t kMaxNSigma = 3.;
	
	const Int_t kMaxNSigmaITS = 3;
	const Int_t kMaxNSigmaTPC = 3;
	const Int_t kMaxNSigmaTOF = 3;
	const Int_t kMinTPCnCls   = 70;
	
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
	
	AliAnalysisTaskB2* task = new AliAnalysisTaskB2(Form("B2_%s",species.Data()));
	
	task->SetParticleSpecies(species);
	task->SetSimulation(simulation);
	task->SetHeavyIons(heavyIons);
	task->SetV0ANDtrigger(V0AND);
	
	task->SetMaxNSigmaITS(kMaxNSigmaITS);
	task->SetMaxNSigmaTPC(kMaxNSigmaTPC);
	task->SetMaxNSigmaTOF(kMaxNSigmaTOF);
	
	task->SetMeanNtrk(GetMeanNtrk(period));
	task->SetKNOmultInterval(minKNOmult, maxKNOmult);
	task->SetVertexXInterval(-kMaxVx, kMaxVx);
	task->SetVertexYInterval(-kMaxVy, kMaxVy);
	task->SetVertexZInterval(-kMaxVz, kMaxVz);
	
	task->SetCentralityInterval(minCentrality, maxCentrality);
	
	if(period=="lhc11a_wsdd" || period=="lhc11a_wosdd")
	{
		task->SetNoFastOnlyTrigger();
	}
	
	// histograms
	
	gROOT->LoadMacro("$ALICE_ROOT/PWGLF/SPECTRA/Nuclei/B2/macros/CreateHistograms.C");
	
	AliLnHistoMap* hMap = CreateHistograms(species, simulation, maxDCAxy, kMaxEta, kMaxY, heavyIons);
	
	task->SetHistogramMap(hMap);
	
	// track selection criteria
	
	gROOT->LoadMacro("$ALICE_ROOT/PWGLF/SPECTRA/Nuclei/B2/macros/TrackCuts.C");
	
	AliESDtrackCuts* trkCuts = TrackCuts(task, trksel, maxDCAxy, maxDCAz, kMaxNSigma, kMinTPCnCls, kMaxEta);
	task->SetESDtrackCuts(trkCuts);
	
	// PID
	
	gROOT->LoadMacro("$ALICE_ROOT/PWGLF/SPECTRA/Nuclei/B2/macros/BetheBlochParams.C");
	
	Double_t bethe[5];
	BetheBlochParams(bethe, period);
	
	gROOT->LoadMacro("$ALICE_ROOT/PWGLF/SPECTRA/Nuclei/B2/macros/PriorProbabilities.C");
	
	Double_t prob[9];
	PriorProbabilities(prob, period, trksel);
	
	AliLnID* lnID = new AliLnID();
	
	lnID->SetTPCBetheBlochParams(bethe);
	lnID->SetPriorProbabilities(prob);
	lnID->SetPidProcedure(pidProc);
	
	task->SetPID(lnID);
	
	// Add task to the manager
	
	mgr->AddTask(task);
	
	// input and output containers
	
	AliAnalysisDataContainer* output = mgr->CreateContainer(outputfile.Data(), AliLnHistoMap::Class(), AliAnalysisManager::kOutputContainer, outputfile.Data());
	
	mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
	mgr->ConnectOutput(task, 0, output);
	
	return task;
}

Double_t GetMeanNtrk(const TString& period)
{
//
// average track multiplicity <Ntrk> for the given period
//
	if(period =="lhc10c900")    return 3.70158; // pass3
	if(period =="lhc10b_pass2") return 9.86258;
	if(period =="lhc10c_pass2") return 9.61402;
	if(period =="lhc10b")       return 5.96104; // pass3
	if(period =="lhc10c")       return 5.94719; // pass3
	if(period =="lhc10d")       return 5.82333; // pass2
	if(period =="lhc10e")       return 5.89367; // pass2
	if(period =="lhc11a_wosdd") return 4.28597; // pass3
	if(period =="lhc11a_wsdd")  return 4.57960; // pass4
	
	// MC
	if(period =="lhc10e13")            return 3.17763;
	if(period =="lhc10f6a")            return 4.41362;
	if(period =="lhc10e21")            return 4.74991;
	if(period =="lhc11e3a_plus_wosdd") return 3.37669;
	if(period =="lhc11e3a_plus_wosdd") return 3.55973;
	
	return 1;
}
