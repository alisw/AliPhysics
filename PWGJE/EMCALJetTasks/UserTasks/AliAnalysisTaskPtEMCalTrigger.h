#ifndef ALIANALYSISTASKPTEMCALTRIGGER_H_
#define ALIANALYSISTASKPTEMCALTRIGGER_H_
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// Author: Markus Fasel

#include "AliAnalysisTaskEmcal.h"
#include "AliCutValueRange.h"
#include "AliESDtrackCuts.h"
#include <TClonesArray.h>
#include <TList.h>

class TArrayD;
class Axis;
class AliESDtrack;
class AliVTrack;
class AliVParticle;

namespace EMCalTriggerPtAnalysis {
class AliEMCalHistoContainer;

class AliAnalysisTaskPtEMCalTrigger : public AliAnalysisTaskEmcal {
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

	virtual void UserCreateOutputObjects();
	virtual Bool_t Run();

	void AddESDTrackCuts(AliESDtrackCuts *trackCuts);
	void SetEtaRange(double etamin, double etamax) { fEtaRange.SetLimits(etamin, etamax); }
	void SetPtRange(double ptmin, double ptmax) { fPtRange.SetLimits(ptmin, ptmax); }
	void SetSwapEta() { fSwapEta = kTRUE; }
	void UseTriggersFromTriggerMaker() { fUseTriggersFromTriggerMaker = kTRUE; }

private:
	AliAnalysisTaskPtEMCalTrigger(const AliAnalysisTaskPtEMCalTrigger &);
	AliAnalysisTaskPtEMCalTrigger &operator=(const AliAnalysisTaskPtEMCalTrigger &);
	void CreateDefaultPtBinning(TArrayD &binning) const;
	void CreateDefaultZVertexBinning(TArrayD &binning) const;
	void CreateDefaultEtaBinning(TArrayD &binning) const;
	void DefineAxis(TAxis &axis, const char *name, const char *title, const TArrayD &binning, const char **labels = NULL);
	void DefineAxis(TAxis &axis, const char *name, const char *title, int nbins, double min, double max, const char **labels = NULL);
	void FillEventHist(const char *trigger, double vz, bool isPileup);
	void FillTrackHist(const char *trigger, const AliVTrack *track, double vz, bool isPileup, int cut, bool isMinBias);
	void FillClusterHist(const char *trigger, const AliVCluster *clust, bool isCalibrated, double vz, bool isPileup, bool isMinBias);
	void FillMCParticleHist(const AliVParticle * const part);
	bool IsTrueTrack(const AliVTrack *const) const;
	TString BuildTriggerString();
	const AliVVertex *GetSPDVertex() const;

	AliEMCalHistoContainer        *fHistos;               //! Histogram container for the task
	TList 						  *fListTrackCuts;		  // List of track cuts

	// Cuts
	AliCutValueRange<double>      fEtaRange;              // Eta Selection Range
	AliCutValueRange<double>	  fPtRange;				  // Pt Selection Range
	Bool_t						  fSwapEta;				  // Allow swapping of the eta sign in asymmetric collision systems
	Bool_t 						  fUseTriggersFromTriggerMaker; // Use trigger classes from trigger maker

	ClassDef(AliAnalysisTaskPtEMCalTrigger, 1);           // Analysis of EMCal triggered events
};

}
#endif
