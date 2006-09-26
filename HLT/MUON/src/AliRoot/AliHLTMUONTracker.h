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
	
	// Do not allow copying of this object.
	AliHLTMUONTracker(const AliHLTMUONTracker& /*object*/)
		: AliTracker(), fdHLT(NULL), fTriggers(NULL), fClusters(NULL),
		  fTracks(NULL)
	{}

	AliHLTMUONTracker& operator = (const AliHLTMUONTracker& /*object*/) { return *this; }


	const AliHLTMUONTriggerRecord*
	FindTriggerRecord(const AliHLTMUONTrack* track) const;

	void LeastSquaresFit(
			const Double_t x[4], const Double_t y[4],
			Double_t& m, Double_t& c
		) const;

	Double_t ComputeChi2(const AliHLTMUONTrack* track) const;

	AliHLTMUONMicrodHLT*     fdHLT;      // dHLT tracker algorithm interface object.
	AliHLTMUONTriggerSource* fTriggers;  // Trigger record input data object.
	AliHLTMUONClusterSource* fClusters;  // Reconstructed hit input data object (sorry about the poor object naming).
	AliHLTMUONTrackSink*     fTracks;    // Track output data object.

	ClassDef(AliHLTMUONTracker, 1)  // dHLT tracker algorithm
};


#endif // ALIHLTMUONTRACKER_H

