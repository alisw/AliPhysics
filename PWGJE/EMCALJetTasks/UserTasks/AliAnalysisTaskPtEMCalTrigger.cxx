/**************************************************************************
 * Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
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

/*
 * Analysis task of the pt analysis on EMCal-triggered events
 *
 *   Author: Markus Fasel
 */

#include <map>
#include <cstring>
#include <iostream>
#include <memory>
#include <vector>
#include <string>
#include <sstream>

#include <TDirectory.h>
#include <TH1.h>
#include <THashList.h>
#include <TKey.h>
#include <TMath.h>
#include <TObjArray.h>
#include <TString.h>

#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliESDtrack.h"
#include "AliESDVertex.h"

#include "AliEMCalHistoContainer.h"
#include "AliAnalysisTaskPtEMCalTrigger.h"

ClassImp(EMCalTriggerPtAnalysis::AliAnalysisTaskPtEMCalTrigger)

namespace EMCalTriggerPtAnalysis {

	//______________________________________________________________________________
	AliAnalysisTaskPtEMCalTrigger::AliAnalysisTaskPtEMCalTrigger():
                		AliAnalysisTaskSE(),
                		fResults(NULL),
                		fHistos(NULL),
                		fListTrackCuts(NULL)
	{
		/*
		 * Dummy constructor, initialising the values with default (NULL) values
		 */
	}

	//______________________________________________________________________________
	AliAnalysisTaskPtEMCalTrigger::AliAnalysisTaskPtEMCalTrigger(const char *name):
                		AliAnalysisTaskSE(name),
                		fResults(NULL),
                		fHistos(NULL),
                		fListTrackCuts(NULL)
	{
		/*
		 * Main constructor, setting default values for eta and zvertex cut
		 */
		DefineOutput(1, TList::Class());

		fListTrackCuts = new TList;
		fListTrackCuts->SetOwner(false);

		// Set default cuts
		fEtaRange.SetLimits(-0.8, 0.8);

	}

	//______________________________________________________________________________
	AliAnalysisTaskPtEMCalTrigger::~AliAnalysisTaskPtEMCalTrigger(){
		/*
		 * Destructor, deleting output
		 */
		//if(fTrackSelection) delete fTrackSelection;
		if(fHistos) delete fHistos;
		if(fListTrackCuts) delete fListTrackCuts;
	}

	//______________________________________________________________________________
	void AliAnalysisTaskPtEMCalTrigger::UserCreateOutputObjects(){
		/*
		 * Create the list of output objects and define the histograms.
		 * Also adding the track cuts to the list of histograms.
		 */
		fResults = new TList;
		fResults->SetOwner();

		fHistos = new AliEMCalHistoContainer("PtEMCalTriggerHistograms");
		fHistos->ReleaseOwner();

		std::map<std::string, std::string> triggerCombinations;
		const char *triggernames[6] = {"MinBias", "EMCJHigh", "EMCJLow", "EMCGHigh", "EMCGLow", "NoEMCal"};
		// Define axes for the trigger correlation histogram
		const TAxis *triggeraxis[5]; memset(triggeraxis, 0, sizeof(const TAxis *) * 5);
		const char *binlabels[2] = {"OFF", "ON"};
		TAxis mytrgaxis[5];
		for(int itrg = 0; itrg < 5; ++itrg){
			DefineAxis(mytrgaxis[itrg], triggernames[itrg], triggernames[itrg], 2, -0.5, 1.5, binlabels);
			triggeraxis[itrg] = mytrgaxis+itrg;
		}
		// Define names and titles for different triggers in the histogram container
		triggerCombinations.insert(std::pair<std::string,std::string>(triggernames[0], "min. bias events"));
		triggerCombinations.insert(std::pair<std::string,std::string>(triggernames[1], "jet-triggered events (high threshold)"));
		triggerCombinations.insert(std::pair<std::string,std::string>(triggernames[2], "jet-triggered events (low threshold)"));
		triggerCombinations.insert(std::pair<std::string,std::string>(triggernames[3], "jet-triggered events (high threshold)"));
		triggerCombinations.insert(std::pair<std::string,std::string>(triggernames[4], "jet-triggered events (low threshold)"));
		triggerCombinations.insert(std::pair<std::string,std::string>(triggernames[5], "non-EMCal-triggered events (low threshold)"));
		// Define axes for the pt histogram
		// Dimensions:
		// 1. pt
		// 2. eta
		// 3. phi
		// 4. vertex
		// 5. pileup (0 = all events, 1 = after pileup rejection)
		// 6. track cuts (0 = no cuts; 1 = after std cuts)
		TArrayD ptbinning, zvertexBinning, etabinning, pileupaxis(3);
		pileupaxis[0] = -0.5; pileupaxis[1] = 0.5; pileupaxis[2] = 1.5;
		CreateDefaultPtBinning(ptbinning);
		CreateDefaultZVertexBinning(zvertexBinning);
		CreateDefaultEtaBinning(etabinning);
		TAxis htrackaxes[6];
		DefineAxis(htrackaxes[0], "pt", "p_{t} (GeV/c)", ptbinning);
		DefineAxis(htrackaxes[1], "eta", "#eta", etabinning);
		DefineAxis(htrackaxes[2], "phi", "#phi", 100, 0, 2 * TMath::Pi());
		DefineAxis(htrackaxes[3], "zvertex", "z_{V} (cm)", zvertexBinning);
		DefineAxis(htrackaxes[4], "pileup", "Pileup rejection", 2, -0.5, 1.5);
		DefineAxis(htrackaxes[5], "trackcuts", "Track Cuts", (fListTrackCuts ? fListTrackCuts->GetEntries() : 0) + 1, -0.5, (fListTrackCuts ? fListTrackCuts->GetEntries() : 0) + 0.5);
		const TAxis *trackaxes[6];
		for(int iaxis = 0; iaxis < 6; ++iaxis) trackaxes[iaxis] = htrackaxes + iaxis;
		try{
			for(std::map<std::string,std::string>::iterator it = triggerCombinations.begin(); it != triggerCombinations.end(); ++it){
				const std::string name = it->first, &title = it->second;
				// Create event-based histogram
				fHistos->CreateTH2(Form("hEventHist%s", name.c_str()), Form("Event-based data for %s events; pileup rejection; z_{V} (cm)", title.c_str()), pileupaxis, zvertexBinning);
				// Create track-based histogram
				fHistos->CreateTHnSparse(Form("hTrackHist%s", name.c_str()), Form("Track-based data for %s events", title.c_str()), 6, trackaxes);
			}
			fHistos->CreateTHnSparse("hEventTriggers", "Trigger type per event", 5, triggeraxis);
		} catch (HistoContainerContentException &e){
			std::stringstream errormessage;
			errormessage << "Creation of histogram failed: " << e.what();
			AliError(errormessage.str().c_str());
		}
		fResults->Add(fHistos->GetListOfHistograms());
		if(fListTrackCuts && fListTrackCuts->GetEntries()){
			TIter cutIter(fListTrackCuts);
			AliESDtrackCuts *cutObject(NULL);
			while((cutObject = dynamic_cast<AliESDtrackCuts *>(cutIter()))){
				cutObject->DefineHistograms();
				fResults->Add(cutObject);
			}
		}
		PostData(1, fResults);
	}

	//______________________________________________________________________________
	void AliAnalysisTaskPtEMCalTrigger::UserExec(Option_t* /*option*/){
		/*
		 * Runs the event loop
		 *
		 * @param option: Additional options
		 */

		// Common checks: Have SPD vertex and primary vertex from tracks, and both need to have at least one contributor
		AliESDEvent *esd = static_cast<AliESDEvent *>(fInputEvent);
		const AliESDVertex *vtxTracks = esd->GetPrimaryVertex(),
				*vtxSPD = esd->GetPrimaryVertexSPD();
		if(!(vtxTracks && vtxSPD)) return;
		if(vtxTracks->GetNContributors() < 1 || vtxSPD->GetNContributors() < 1) return;

		double triggers[5]; memset(triggers, 0, sizeof(double) *5);
		if(fInputHandler->IsEventSelected() & AliVEvent::kINT7)
			triggers[0] = 1.;

		std::vector<std::string> triggerstrings;
		if(fInputHandler->IsEventSelected() & AliVEvent::kEMC7){
			// EMCal-triggered event, distinguish types
			TString trgstr(fInputEvent->GetFiredTriggerClasses());
			if(trgstr.Contains("EJ1")){
				triggerstrings.push_back("EMCJHigh");
				triggers[1] = 1;
			}
			if(trgstr.Contains("EJ2")){
				triggerstrings.push_back("EMCJLow");
				triggers[2] = 1;
			}
			if(trgstr.Contains("EG1")){
				triggerstrings.push_back("EMCGHigh");
				triggers[3] = 1;
			}
			if(trgstr.Contains("EG2")){
				triggerstrings.push_back("EMCGLow");
				triggers[4] = 1;
			}
		}
		try{
			fHistos->FillTHnSparse("hEventTriggers", triggers);
		} catch (HistoContainerContentException &e){
			std::stringstream errormessage;
			errormessage << "Filling of histogram failed: " << e.what();
			AliError(errormessage.str().c_str());
		}

		// apply event selection: Combine the Pileup cut from SPD with the other pA Vertex selection cuts.
		bool isPileupEvent = esd->IsPileupFromSPD();
		isPileupEvent = isPileupEvent || (TMath::Abs(vtxTracks->GetZ() - vtxSPD->GetZ()) > 0.5);
		double covSPD[6]; vtxSPD->GetCovarianceMatrix(covSPD);
		isPileupEvent = isPileupEvent || (TString(vtxSPD->GetTitle()).Contains("vertexer:Z") && TMath::Sqrt(covSPD[5]) > 0.25);

		// Fill event-based histogram
		const double &zv = vtxTracks->GetZ();
		if(triggers[0]) FillEventHist("MinBias", zv, isPileupEvent);
		if(!triggerstrings.size()) // Non-EMCal-triggered
			FillEventHist("NoEMCal", zv, isPileupEvent);
		else{
			// EMCal-triggered events
			for(std::vector<std::string>::iterator it = triggerstrings.begin(); it != triggerstrings.end(); ++it)
				FillEventHist(it->c_str(), zv, isPileupEvent);
		}

		AliESDtrack *track(NULL);
		// Loop over all tracks (No cuts applied)
		for(int itrk = 0; itrk < fInputEvent->GetNumberOfTracks(); ++itrk){
			track = dynamic_cast<AliESDtrack *>(fInputEvent->GetTrack(itrk));
			// first fill without pielup cut
			if(fEtaRange.IsInRange(track->Eta())) continue;
			if(triggers[0]) FillTrackHist("MinBias", track, zv, isPileupEvent, 0);
			if(!triggerstrings.size()) // Non-EMCal-triggered
				FillTrackHist("NoEMCal", track, zv, isPileupEvent, 0);
			else {
				// EMCal-triggered events
				for(std::vector<std::string>::iterator it = triggerstrings.begin(); it != triggerstrings.end(); ++it)
					FillTrackHist(it->c_str(), track, zv, isPileupEvent, 0);
			}
		}

		// Now apply track selection cuts
		// allow for several track selections to be tested at the same time
		// each track selection gets a different cut ID starting from 1
		// cut ID 0 is reserved for the case of no cuts
		if(fListTrackCuts && fListTrackCuts->GetEntries()){
			for(int icut = 0; icut < fListTrackCuts->GetEntries(); icut++){
				AliESDtrackCuts *trackSelection = static_cast<AliESDtrackCuts *>(fListTrackCuts->At(icut));
				std::auto_ptr<TObjArray> acceptedTracks(trackSelection->GetAcceptedTracks(esd));
				TIter trackIter(acceptedTracks.get());
				while((track = dynamic_cast<AliESDtrack *>(trackIter()))){
					if(!fEtaRange.IsInRange(track->Eta())) continue;
					if(triggers[0]) FillTrackHist("MinBias", track, zv, isPileupEvent, icut + 1);
					if(!triggerstrings.size()) // Non-EMCal-triggered
						FillTrackHist("NoEMCal", track, zv, isPileupEvent, icut + 1);
					else {
						// EMCal-triggered events
						for(std::vector<std::string>::iterator it = triggerstrings.begin(); it != triggerstrings.end(); ++it)
							FillTrackHist(it->c_str(), track, zv, isPileupEvent, icut + 1);
					}
				}
			}
		}
		PostData(1, fResults);
	}

	//______________________________________________________________________________
	void AliAnalysisTaskPtEMCalTrigger::CreateDefaultPtBinning(TArrayD &binning) const{
		/*
		 * Creating the default pt binning.
		 *
		 * @param binning: Array where to store the results.
		 */
		std::vector<double> mybinning;
		std::map<double,double> definitions;
		definitions.insert(std::pair<double,double>(2.5, 0.1));
		definitions.insert(std::pair<double,double>(7., 0.25));
		definitions.insert(std::pair<double,double>(15., 0.5));
		definitions.insert(std::pair<double,double>(25., 1.));
		definitions.insert(std::pair<double,double>(40., 2.5));
		definitions.insert(std::pair<double,double>(60., 5.));
		definitions.insert(std::pair<double,double>(100., 5.));
		double currentval = 0;
		for(std::map<double,double>::iterator id = definitions.begin(); id != definitions.end(); ++id){
			double limit = id->first, binwidth = id->second;
			while(currentval <= limit){
				currentval += binwidth;
				mybinning.push_back(currentval);
			}
		}
		binning.Set(mybinning.size());
		int ib = 0;
		for(std::vector<double>::iterator it = mybinning.begin(); it != mybinning.end(); ++it)
			binning[ib++] = *it;
	}

	//______________________________________________________________________________
	void AliAnalysisTaskPtEMCalTrigger::CreateDefaultZVertexBinning(TArrayD &binning) const {
		/*
		 * Creating default z-Vertex binning.
		 *
		 * @param binning: Array where to store the results.
		 */
		std::vector<double> mybinning;
		double currentval = -40;
		mybinning.push_back(currentval);
		while(currentval <= 40.){
			currentval += 0.1;
			mybinning.push_back(currentval);
		}
		binning.Set(mybinning.size());
		int ib = 0;
		for(std::vector<double>::iterator it = mybinning.begin(); it != mybinning.end(); ++it)
			binning[ib++] = *it;
	}

	//______________________________________________________________________________
	void AliAnalysisTaskPtEMCalTrigger::CreateDefaultEtaBinning(TArrayD& binning) const {
		/*
		 * Creating default z-Vertex binning.
		 *
		 * @param binning: Array where to store the results.
		 */
		std::vector<double> mybinning;
		double currentval = -0.8;
		mybinning.push_back(currentval);
		while(currentval <= 0.8){
			currentval += 0.05;
			mybinning.push_back(currentval);
		}
		binning.Set(mybinning.size());
		int ib = 0;
		for(std::vector<double>::iterator it = mybinning.begin(); it != mybinning.end(); ++it)
			binning[ib++] = *it;
	}

	//______________________________________________________________________________
	void AliAnalysisTaskPtEMCalTrigger::DefineAxis(TAxis& axis, const char* name,
			const char* title, const TArrayD& binning, const char** labels) {
		/*
		 * Define an axis with a given binning
		 *
		 * @param axis: Axis to be defined
		 * @param name: Name of the axis
		 * @param title: Title of the axis
		 * @param binning: axis binning
		 * @param labels (@optional): array of bin labels
		 */
		axis.Set(binning.GetSize()-1, binning.GetArray());
		axis.SetName(name);
		axis.SetTitle(title);
		if(labels){
			for(int ib = 1; ib <= axis.GetNbins(); ++ib)
				axis.SetBinLabel(ib, labels[ib]);
		}
	}

	//______________________________________________________________________________
	void AliAnalysisTaskPtEMCalTrigger::DefineAxis(TAxis& axis, const char* name,
			const char* title, int nbins, double min, double max,
			const char** labels) {
		/*
		 * Define an axis with number of bins from min to max
		 *
		 * @param axis: Axis to be defined
		 * @param name: Name of the axis
		 * @param title: Title of the axis
		 * @param nbins: Number of bins
		 * @param min: lower limit of the axis
		 * @param max: upper limit of the axis
		 * @param labels (@optional): array of bin labels
		 */
		axis.Set(nbins, min, max);
		axis.SetName(name);
		axis.SetTitle(title);
		if(labels){
			for(int ib = 1; ib <= axis.GetNbins(); ++ib)
				axis.SetBinLabel(ib, labels[ib]);
		}
	}

	//______________________________________________________________________________
	void AliAnalysisTaskPtEMCalTrigger::FillEventHist(const char* trigger,
			double vz, bool isPileup) {
		/*
		 * Fill event-based histogram
		 *
		 * @param trigger: name of the trigger configuration to be processed
		 * @param vz: z-position of the vertex
		 * @param isPileup: signalises if the event is flagged as pileup event
		 */
		char histname[1024];
		sprintf(histname, "hEventHist%s", trigger);
		try{
			fHistos->FillTH2(histname, vz, 0);
		} catch (HistoContainerContentException &e){
			std::stringstream errormessage;
			errormessage << "Filling of histogram failed: " << e.what();
			AliError(errormessage.str().c_str());
		}
		if(!isPileup){
			try{
				fHistos->FillTH2(histname, vz, 1);
			} catch(HistoContainerContentException &e){
				std::stringstream errormessage;
				errormessage << "Filling of histogram failed: " << e.what();
				AliError(errormessage.str().c_str());
			}
		}
	}

	//______________________________________________________________________________
	void AliAnalysisTaskPtEMCalTrigger::FillTrackHist(const char* trigger,
			const AliESDtrack* track, double vz, bool isPileup, int cut) {
		/*
		 * Fill track-based histogram with corresponding information
		 *
		 * @param trigger: name of the trigger
		 * @param track: ESD track selected
		 * @param vz: z-position of the vertex
		 * @param isPileup: flag event as pileup event
		 * @param cut: id of the cut (0 = no cut)
		 */
		double data[6] = {track->Pt(), track->Eta(), track->Phi(), vz, 0, cut};
		char histname[1024];
		sprintf(histname, "hTrackHist%s", trigger);
		try{
			fHistos->FillTHnSparse(histname, data);
		} catch (HistoContainerContentException &e){
			std::stringstream errormessage;
			errormessage << "Filling of histogram failed: " << e.what();
			AliError(errormessage.str().c_str());
		}
		if(!isPileup){
			data[4] = 1;
			try{
				fHistos->FillTHnSparse(histname, data);
			} catch (HistoContainerContentException &e){
				std::stringstream errormessage;
				errormessage << "Filling of histogram failed: " << e.what();
				AliError(errormessage.str().c_str());
			}
		}
	}

}

