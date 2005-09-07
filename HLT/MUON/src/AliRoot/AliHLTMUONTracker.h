// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com

#ifndef ALIHLTMUONTRACKER_H
#define ALIHLTMUONTRACKER_H

#include "AliTracker.h"
#include "AliLog.h"

#include "AliRoot/MicrodHLT.hpp"
#include "AliRoot/TriggerSource.hpp"
#include "AliRoot/ClusterSource.hpp"
#include "AliRoot/TrackSink.hpp"

class AliRunLoader;
class AliESD;

/* This class is a wrapper for the dHLT tracker implemented in the MicrodHLT
   object. It is used by the AliReconstruction framework to indirectly run the
   dHLT tracking algorithm.
 */

class AliHLTMUONTracker : public AliTracker
{
public:
	AliHLTMUONTracker(AliRunLoader* runloader);
	virtual ~AliHLTMUONTracker();

	// Inherited methods.
	virtual Int_t PropagateBack(AliESD* /*event*/) { return 0; };
	virtual Int_t RefitInward(AliESD* /*event*/) { return 0; };
	virtual Int_t LoadClusters(TTree* data);
	virtual void UnloadClusters();
	virtual AliCluster* GetCluster(Int_t /*i*/) const { return NULL; };
	virtual Int_t Clusters2Tracks(AliESD* event);

private:
	const AliMUONHLT::TriggerRecord*
	FindTriggerRecord(const AliMUONHLT::Track* track) const;

	void LeastSquaresFit(
			const Double_t x[4], const Double_t y[4],
			Double_t& m, Double_t& c
		) const;

	Double_t ComputeChi2(const AliMUONHLT::Track* track) const;

	AliMUONHLT::MicrodHLT*     fdHLT;      // dHLT tracker algorithm interface object.
	AliMUONHLT::TriggerSource* fTriggers;  // Trigger record input data object.
	AliMUONHLT::ClusterSource* fClusters;  // Reconstructed hit input data object (sorry about the poor object naming).
	AliMUONHLT::TrackSink*     fTracks;    // Track output data object.

	ClassDef(AliHLTMUONTracker, 1)  // dHLT tracker algorithm
};


#endif // ALIHLTMUONTRACKER_H

