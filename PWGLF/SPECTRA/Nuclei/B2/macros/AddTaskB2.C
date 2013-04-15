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

Double_t GetMeanNtrk(const TString& period, Double_t eta);
Double_t GetNSDMeanNtrk(const TString& period, Double_t eta);

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
                             , Double_t minKNOmult       = -10
                             , Double_t maxKNOmult       = 10000
                             , Bool_t   V0AND            = kFALSE
                             , const TString& ztag       = ""
                             , Double_t maxVz            = 10
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
	
	AliAnalysisTaskB2* task = new AliAnalysisTaskB2(Form("B2.%s",containername.Data()));
	
	task->SetParticleSpecies(species);
	task->SetSimulation(simulation);
	task->SetHeavyIons(heavyIons);
	task->SetV0ANDtrigger(V0AND);
	
	task->SetMaxNSigmaITS(kMaxNSigmaITS);
	task->SetMaxNSigmaTPC(kMaxNSigmaTPC);
	task->SetMaxNSigmaTOF(kMaxNSigmaTOF);
	
	Double_t meanNtrk = V0AND ? GetNSDMeanNtrk(period, maxEta) : GetMeanNtrk(period, maxEta);
	
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
	
	// histograms
	
	gROOT->LoadMacro("$ALICE_ROOT/PWGLF/SPECTRA/Nuclei/B2/macros/CreateHistograms.C");
	
	AliLnHistoMap* hMap = CreateHistograms(species, simulation, maxDCAxy, maxEta, maxY, heavyIons);
	
	task->SetHistogramMap(hMap);
	
	// track selection criteria
	
	gROOT->LoadMacro("$ALICE_ROOT/PWGLF/SPECTRA/Nuclei/B2/macros/TrackCuts.C");
	
	AliESDtrackCuts* trkCuts = TrackCuts(task, trksel, maxDCAxy, maxDCAz, kMaxNSigma, kMinTPCnCls, maxEta);
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

Double_t GetMeanNtrk(const TString& period, Double_t eta)
{
//
// average track multiplicity <Ntrk> for the given period and eta
//
	if(TMath::Abs(eta) > 0.51) // |eta|<0.8
	{
		if(period =="lhc10d")       return 9.47466; // pass2
	}
	else // |eta|<0.5
	{
		if(period =="lhc10c900")    return 3.70158; // pass3
		if(period =="lhc10b_pass2") return 9.86258;
		if(period =="lhc10c_pass2") return 9.61402;
		if(period =="lhc10b")       return 5.96104; // pass3
		if(period =="lhc10c")       return 5.94719; // pass3
		if(period =="lhc10d")       return 5.82333; // pass2
		if(period =="lhc10e")       return 5.89367; // pass2
		if(period =="lhc11a_wosdd") return 4.28597; // pass3
		if(period =="lhc11a_wsdd")  return 4.69927; // pass4
		
		// MC
		if(period =="lhc10e13")            return 3.13712;
		if(period =="lhc10f6a")            return 4.41362;
		if(period =="lhc10e21")            return 4.74991;
		if(period =="lhc11e3a_plus_wosdd") return 3.37669;
		if(period =="lhc11e3a_plus_wsdd")  return 3.47885;
		
		if(period =="lhc12a5a")            return 29.264;
		if(period =="lhc12a5bb")           return 31.0288;
		if(period =="lhc12a5bc")           return 30.6888;
		if(period =="lhc12a5bd")           return 30.3528;
		if(period =="lhc12a5be")           return 29.9859;
		if(period =="lhc12a5c_wsdd")       return 27.5981;
	}
	
	cerr << "Warning in GetMeanNtrk: no <Ntrk> for period " << period << " and |eta| < " << eta << endl;
	
	return 1;
}

Double_t GetNSDMeanNtrk(const TString& period, Double_t eta)
{
//
// average track multiplicity <Ntrk> for the given period and eta
// (NSD events)
//
	if(TMath::Abs(eta) > 0.51) // |eta|<0.8
	{
		if(period =="lhc10d")       return 9.77129; // pass2
	}
	else // |eta|<0.5
	{
		if(period =="lhc10c900")    return 3.84362; // pass3
		if(period =="lhc10b")       return 6.18470; // pass3
		if(period =="lhc10c")       return 6.16175; // pass3
		if(period =="lhc10d")       return 6.03108; // pass2
		if(period =="lhc10e")       return 6.10384; // pass2
		if(period =="lhc11a_wosdd") return 4.40312; // pass3
		if(period =="lhc11a_wsdd")  return 4.87609; // pass4
		
		// MC
		if(period =="lhc10e13")            return 3.33273;
		if(period =="lhc10f6a")            return 4.80771;
		if(period =="lhc10e21")            return 4.91967;
		if(period =="lhc11e3a_plus_wosdd") return 3.4774;
		if(period =="lhc11e3a_plus_wsdd")  return 3.6467;
		
		if(period =="lhc12a5a")            return 29.3414;
		if(period =="lhc12a5bb")           return 31.1514;
		if(period =="lhc12a5bc")           return 30.7877;
		if(period =="lhc12a5bd")           return 30.4706;
		if(period =="lhc12a5be")           return 30.1013;
		if(period =="lhc12a5c_wsdd")       return 27.6921;
	}
	
	cerr << "Warning in GetNSDMeanNtrk: no <Ntrk> for period " << period << " and |eta| < " << eta << endl;
	
	return 1;
}
