/*
 * AliEmcalMCTreeWriter.h
 *
 *  Created on: 03.10.2014
 *      Author: markusfasel
 */

#ifndef ALIEMCALMCTREEWRITER_H_
#define ALIEMCALMCTREEWRITER_H_

#include <cstring>
#include <TVector.h>
#include "AliAnalysisTaskEmcal.h"

class AliVCluster;
class AliVParticle;
class AliVTrack;
class TTree;
class vector;

class AliEmcalMCTreeWriter : public AliAnalysisTaskEmcal {
public:
	AliEmcalMCTreeWriter();
	AliEmcalMCTreeWriter(const char *name);
	virtual ~AliEmcalMCTreeWriter();

	virtual void UserCreateOutputObjects();
	virtual Bool_t Run();

protected:
	double GetClusterEnergy(AliVCluster * clust, Double_t *vpos, Double_t *evec) const;
	void GetShowerShape(AliVCluster *clust, Double_t *vector) const;
	void GetCellEnergies(AliVCluster *clust, TVectorD &cells, TVectorF &indices) const;
	void FindTracks(unsigned int label, std::vector<AliVTrack *> &tracks);
	void FindClusters(unsigned int label, std::vector<AliVCluster *> &clusters);
	bool AcceptParticle(const AliVParticle * const particle) const;

private:
	struct TrackInfo{
		Int_t pdg;
		Bool_t isPhysicalPrimary;
		Bool_t isUnique;
		Double_t Egen;
		Double_t E;
		Double_t pgen[3];
		Double_t prec[3];
		Double_t erec[3];
		Double_t showershape[2];		// Defined as M02 at first and M20 at second entry
		TVectorF fCellIndices;
		TVectorD fCellEnergies;

		TrackInfo(): pdg(0), isPhysicalPrimary(false), isUnique(true), Egen(0.), E(0.), fCellIndices(), fCellEnergies() {
			memset(pgen, 0, sizeof(Double_t) * 3);
			memset(prec, 0, sizeof(Double_t) * 3);
			memset(erec, 0, sizeof(Double_t) * 3);
			memset(showershape, 0, sizeof(Double_t) * 2);
		}
		void Reset(){
			pdg = 0;
			isPhysicalPrimary = false;
			isUnique = true;
			E = 0.;
			Egen = 0.;
			memset(pgen, 0, sizeof(Double_t) * 3);
			memset(prec, 0, sizeof(Double_t) * 3);
			memset(erec, 0, sizeof(Double_t) * 3);
			memset(showershape, 0, sizeof(Double_t) * 2);
			for(int i = 0; i < fCellIndices.GetNrows(); i++) fCellIndices[i] = 0;
			for(int i = 0; i < fCellEnergies.GetNrows(); i++) fCellEnergies[i] = 0;
		}
	};
	TTree *fOutputTree;						//! Output tree with tracks
	TrackInfo fOutputInfo;					// Track Info for the tree

	ClassDef(AliEmcalMCTreeWriter, 1)
};

#endif /* ALIEMCALMCTREEWRITER_H_ */
