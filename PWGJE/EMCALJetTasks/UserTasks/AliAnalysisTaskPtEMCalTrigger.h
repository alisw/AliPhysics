#ifndef ALIANALYSISTASKPTEMCALTRIGGER_H_
#define ALIANALYSISTASKPTEMCALTRIGGER_H_
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// Author: Markus Fasel

#include "AliAnalysisTaskSE.h"
#include "AliCutValueRange.h"
#include "AliESDtrackCuts.h"
#include <TList.h>

class TArrayD;
class Axis;
class AliESDtrack;

namespace EMCalTriggerPtAnalysis {
class AliEMCalHistoContainer;

class AliAnalysisTaskPtEMCalTrigger : public AliAnalysisTaskSE {
public:
	enum EEMCalTriggerType_t{
		kEMCalJetLow = 0,
		kEMCalJetHigh = 1,
		kEMCalGammaLow = 2,
		kEMCalGammaHigh = 3
	};
	AliAnalysisTaskPtEMCalTrigger();
	AliAnalysisTaskPtEMCalTrigger(const char *name);
	~AliAnalysisTaskPtEMCalTrigger();

	void UserCreateOutputObjects();
	void UserExec(Option_t* /*option*/);
	void Terminate(Option_t * /*option*/) {}

	void AddTrackCuts(AliESDtrackCuts *trackCuts) { fListTrackCuts->Add(trackCuts); }
	void SetEtaRange(double etamin, double etamax) { fEtaRange.SetLimits(etamin, etamax); }

private:
	AliAnalysisTaskPtEMCalTrigger(const AliAnalysisTaskPtEMCalTrigger &);
	AliAnalysisTaskPtEMCalTrigger &operator=(const AliAnalysisTaskPtEMCalTrigger &);
	void CreateDefaultPtBinning(TArrayD &binning) const;
	void CreateDefaultZVertexBinning(TArrayD &binning) const;
	void CreateDefaultEtaBinning(TArrayD &binning) const;
	void DefineAxis(TAxis &axis, const char *name, const char *title, const TArrayD &binning, const char **labels = NULL);
	void DefineAxis(TAxis &axis, const char *name, const char *title, int nbins, double min, double max, const char **labels = NULL);
	void FillEventHist(const char *trigger, double vz, bool isPileup);
	void FillTrackHist(const char *trigger, const AliESDtrack *track, double vz, bool isPileup, int cut);

	TList                         *fResults;              //! container for results
	AliEMCalHistoContainer        *fHistos;               //! Histogram container for the task
	TList 						  *fListTrackCuts;		  // List of track cuts

	// Cuts
	AliCutValueRange<double>      fEtaRange;              // Eta Selection Range

	ClassDef(AliAnalysisTaskPtEMCalTrigger, 1);           // Analysis of EMCal triggered events
};

}
#endif
