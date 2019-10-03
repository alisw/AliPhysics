#ifndef ALIANALYSISTASKLINKTOMC_H
#define ALIANALYSISTASKLINKTOMC_H
/* This file is property of and copyright by the ALICE HLT Project        *
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

/* $Id$ */ 

///
/// @file   AliAnalysisTaskLinkToMC.h
/// @author Artur Szostak <artursz@iafrica.com>
/// @date   27 Oct 2008
/// @brief  Definition of the AliAnalysisTaskLinkToMC task to connect reconstructed tracks with MC tracks.
///

#include "AliAnalysisTaskSE.h"

class TList;
class TH2D;
class TMap;
class AliESDMuonTrack;
class AliMCParticle;

/**
 * \class AliAnalysisTaskLinkToMC
 * Connects reconstructed tracks found in the ESD to their corresponding Monte
 * Carlo tracks via the proximity of their hit coordinates.
 */
class AliAnalysisTaskLinkToMC : public AliAnalysisTaskSE
{
public:
	AliAnalysisTaskLinkToMC();
	AliAnalysisTaskLinkToMC(const char* name);

	virtual void UserCreateOutputObjects();
	virtual void UserExec(Option_t* option);
	virtual void Terminate(Option_t* option);
	
	/// Returns the absolute maximum distance in centimetres between the
	/// cluster and track reference that is allowed in the X direction.
	Double_t HardCutLimitX() const { return fHardCutLimitX; }
	
	/// Returns the absolute maximum distance in centimetres between the
	/// cluster and track reference that is allowed in the Y direction.
	Double_t HardCutLimitY() const { return fHardCutLimitY; }
	
	/// Returns the absolute maximum distance in centimetres between the
	/// cluster and track reference that is allowed in the Z direction.
	Double_t HardCutLimitZ() const { return fHardCutLimitZ; }
	
	void HardCutLimitX(Double_t value);
	void HardCutLimitY(Double_t value);
	void HardCutLimitZ(Double_t value);
	
	Double_t SigmaCut() const;
	
	/// Sets the maximum sigma value to allow for each X and Y dimension
	/// between clusters and track references.
	void SigmaCut(Double_t value) { fSigmaCut2 = value*value; }
	
	/// Returns the minimum number of clusters that must match between the
	/// ESD track clusters and track references.
	Int_t MinClusters() const { return fMinClusters; }
	
	/// Sets the minimum number of clusters that must match between the ESD
	/// track clusters and track references.
	void MinClusters(Int_t value) { fMinClusters = value; }
	
	/// Returns the minimum number of clusters that must match between the
	/// ESD track clusters and track references in stations 1 and 2.
	Int_t MinClustersInSt12() const { return fMinClustersInSt12; }
	
	/// Sets the minimum number of clusters that must match between the ESD
	/// track clusters and track references in stations 1 and 2.
	void MinClustersInSt12(Int_t value) { fMinClustersInSt12 = value; }
	
	/// Returns the minimum number of clusters that must match between the
	/// ESD track clusters and track references in stations 4 and 5.
	Int_t MinClustersInSt45() const { return fMinClustersInSt45; }
	
	/// Sets the minimum number of clusters that must match between the ESD
	/// track clusters and track references in stations 4 and 5.
	void MinClustersInSt45(Int_t value) { fMinClustersInSt45 = value; }
	
	/// Returns the flag indicating if histograms should be generated for
	/// checking the quality of track matching.
	Bool_t GenerateHistograms() const { return fMakeHists; }
	
	/// Set the flag indicating if histograms should be generated for checking
	/// the quality of track matching.
	void GenerateHistograms(Bool_t value) { fMakeHists = value; }
	
	/// Returns the flag indicating if the Z coordinate is to be used for
	/// track matching or not.
	Bool_t UseZCoordinate() const { return fUseZCoordinate; }
	
	/// Set the flag indicating if the Z coordinate should be used for track
	/// matching.
	void UseZCoordinate(Bool_t value) { fUseZCoordinate = value; }
	
	Bool_t ChamberMustMatch(Int_t chamber) const;
	void ChamberMustMatch(Int_t chamber, Bool_t value);
	Bool_t StationMustMatch(Int_t station) const;
	void StationMustMatch(Int_t station, Bool_t value);

private:

	/// Do not allow copying of this object. Method is not implemented.
	AliAnalysisTaskLinkToMC(const AliAnalysisTaskLinkToMC& obj);
	/// Do not allow copying of this object. Method is not implemented.
	AliAnalysisTaskLinkToMC& operator = (const AliAnalysisTaskLinkToMC& obj);
	
	void FindLinksToMC(TMap& links);
	void CalculateTrackCompatibility(
			AliESDMuonTrack* esdTrack, AliMCParticle* mcTrack,
			bool& tracksCompatible, Double_t& chi2, Int_t& noOfClusters
		);
	
	bool IsFindable(AliMCParticle* track) const;
	
	void FillHistsFromMC();
	void FillHistsFromLinks(TMap& links);
	void CreateAODTracks(TMap& links);

	TList* fHistos;  //! List of all histograms produced by this analysis task.
	TH2D* fFindableHist;  //! Histogram of reconstructible tracks (filled with MC information).
	TH2D* fFoundHistMC;  //! Histogram of reconstructed tracks (filled with MC information).
	TH2D* fFoundHist;  //! Histogram of reconstructed tracks (filled with reconstructed information).
	TH2D* fFakeHist;  //! Histogram of fake tracks (filled with reconstructed information).
	
	Double_t fHardCutLimitX; //! The absolute maximum difference between a cluster and reference track point in X direction in centimetres.
	Double_t fHardCutLimitY; //! The absolute maximum difference between a cluster and reference track point in Y direction in centimetres.
	Double_t fHardCutLimitZ; //! The absolute maximum difference between a cluster and reference track point in Z direction in centimetres.
	Double_t fVarX; //! The variance to use in the X direction if error value in the cluster was not set.
	Double_t fVarY; //! The variance to use in the Y direction if error value in the cluster was not set.
	Double_t fSigmaCut2; //! The cut on the chi2 X/Y direction component to make for each cluster fit.
	Int_t fMinClusters; //! Minimum number of clusters that must match between an ESD and MC track for it to be considered compatible.
	Int_t fMinClustersInSt12;  //! The minimum number of clusters that must match in stations 1 and 2.
	Int_t fMinClustersInSt45;  //! The minimum number of clusters that must match in stations 4 and 5.
	Bool_t fMakeHists;  //! Flag indicating if histogram filling and generation should be performed.
	Bool_t fUseZCoordinate;  //! Flag indicating if the Z coordinate should be used for comparrison.
	Bool_t fChamberMustMatch[14];  //! A flag for each chamber which must have a compatible hit matched.
	Bool_t fStationMustMatch[7];  //! A flag for each station which must have at least one compatible hit matched.

	ClassDef(AliAnalysisTaskLinkToMC, 0) // Analysis task to connect reconstructed tracks to their Monte Carlo counterparts.
};

#endif // ALIANALYSISTASKLINKTOMC_H

