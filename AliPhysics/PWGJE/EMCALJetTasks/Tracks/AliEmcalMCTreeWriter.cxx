/*
 * AliEmcalMCTreeWriter.cxx
 *
 *  Created on: 03.10.2014
 *      Author: markusfasel
 */

#include <iostream>
#include <string>
#include <TArrayI.h>
#include <TClonesArray.h>
#include <TLorentzVector.h>
#include <TPDGCode.h>
#include <TMath.h>
#include <TTree.h>
#include <TVector3.h>

#include "AliAODMCParticle.h"
#include "AliMCEvent.h"
#include "AliVParticle.h"
#include "AliVTrack.h"
#include "AliVCluster.h"

#include "AliEmcalMCTreeWriter.h"

AliEmcalMCTreeWriter::AliEmcalMCTreeWriter():
	AliAnalysisTaskEmcal(),
	fOutputTree(NULL),
	fOutputInfo()
{
	/*
	 * Dummy constructor
	 */
}

AliEmcalMCTreeWriter::AliEmcalMCTreeWriter(const char *name):
	AliAnalysisTaskEmcal(name,true),
	fOutputTree(NULL),
	fOutputInfo()
{
	/*
	 * Constructor
	 */
	DefineOutput(2,TTree::Class());
}

AliEmcalMCTreeWriter::~AliEmcalMCTreeWriter() {
	/*
	 * Destructor
	 */
	if(fOutputTree) delete fOutputTree;
}

void AliEmcalMCTreeWriter::UserCreateOutputObjects() {
	/*
	 * Create output tree, with two branches, one for the tracks matched and one for the clusters
	 */
	AliAnalysisTaskEmcal::UserCreateOutputObjects();

	// set the listnames
	std::string clusterlist(""), tracklist("");
	if(fIsEsd){
		tracklist = "ESDFilterTracks";
		//clusterlist = "EmcCaloClusters";
		clusterlist = "CaloClustersCorr";
	} else {
		tracklist = "AODFilterTracks";
		//clusterlist = "EmcCaloClusters";
		clusterlist = "CaloClustersCorr";
	}
	this->SetTracksName(tracklist.c_str());
	this->SetClusName(clusterlist.c_str());

	// Build the tree
	OpenFile(2);
	fOutputTree = new TTree("EMCalTree", "A tree with emcal information");
	fOutputTree->Branch("pdg", &fOutputInfo.pdg, "pdg/I");
	fOutputTree->Branch("energy", &fOutputInfo.E, "energy/D");
	fOutputTree->Branch("energyGen", &fOutputInfo.Egen, "energyGen/D");
	fOutputTree->Branch("isUnique", &fOutputInfo.isUnique, "isUnique/B");
	fOutputTree->Branch("isPhysicalPrimary", &fOutputInfo.isPhysicalPrimary, "isPhysicalPrimary/B");
	fOutputTree->Branch("pgx", &fOutputInfo.pgen[0], "pgx/D");
	fOutputTree->Branch("pgy", &fOutputInfo.pgen[1], "pgy/D");
	fOutputTree->Branch("pgz", &fOutputInfo.pgen[2], "pgz/D");
	fOutputTree->Branch("prx", &fOutputInfo.prec[0], "prx/D");
	fOutputTree->Branch("pry", &fOutputInfo.prec[1], "pry/D");
	fOutputTree->Branch("prz", &fOutputInfo.prec[2], "prz/D");
	fOutputTree->Branch("erx", &fOutputInfo.erec[0], "erx/D");
	fOutputTree->Branch("ery", &fOutputInfo.erec[1], "ery/D");
	fOutputTree->Branch("erz", &fOutputInfo.erec[2], "erz/D");
	fOutputTree->Branch("m02", &fOutputInfo.showershape[0], "m02/D");
	fOutputTree->Branch("m20", &fOutputInfo.showershape[1], "m20/D");
	fOutputTree->Branch("cells", "TVectorD", &fOutputInfo.fCellEnergies);
	fOutputTree->Branch("indices", "TVectorD", &fOutputInfo.fCellIndices);
	PostData(1, fOutput);
	PostData(2, fOutputTree);
}

Bool_t AliEmcalMCTreeWriter::Run() {
	/*
	 * Build the tree
	 */

	const AliVVertex *vertex = fInputEvent->GetPrimaryVertex();
	if(!vertex) return kFALSE;

	Double_t vpos[3];
	vertex->GetXYZ(vpos);
	TIter trackIter(fTracks);
	AliVTrack *rectrack(NULL);
	AliVParticle *particle(NULL);
	AliAODMCParticle *aodmc(NULL);
	AliVCluster *assoccluster(NULL);
	for(unsigned int ipart = 0; ipart < static_cast<unsigned int>(fMCEvent->GetNumberOfTracks()); ipart++){
		particle = fMCEvent->GetTrack(static_cast<int>(ipart));
		if(!AcceptParticle(particle)) continue;
		if(TMath::Abs(particle->Eta()) > 2) continue;		// Reject particles far outside the alice acceptance

		fOutputInfo.Reset();
		fOutputInfo.pdg = particle->PdgCode();
		if((aodmc = dynamic_cast<AliAODMCParticle *>(particle)))
			fOutputInfo.isPhysicalPrimary = aodmc->IsPhysicalPrimary();
		else
			fOutputInfo.isPhysicalPrimary = fMCEvent->IsPhysicalPrimary(static_cast<int>(ipart));
		particle->PxPyPz(fOutputInfo.pgen);
		fOutputInfo.Egen = particle->E();

		// First try to find a matching track
		std::vector<AliVTrack *> tracks;
		FindTracks(ipart, tracks);
		if(tracks.size() > 1) fOutputInfo.isUnique = false;
		for(std::vector<AliVTrack *>::iterator trackIter = tracks.begin(); trackIter != tracks.end(); trackIter++){
			rectrack = *trackIter;
			rectrack->GetPxPyPz(fOutputInfo.prec);
			memset(fOutputInfo.erec, 0, sizeof(Double_t) * 3);
			fOutputInfo.E = 0;
			assoccluster = dynamic_cast<AliVCluster *>(fCaloClusters->At(rectrack->GetEMCALcluster()));
			if(assoccluster){
				fOutputInfo.E = GetClusterEnergy(assoccluster, vpos, fOutputInfo.erec);
				GetShowerShape(assoccluster, fOutputInfo.showershape);
				GetCellEnergies(assoccluster, fOutputInfo.fCellEnergies, fOutputInfo.fCellIndices);
			}
			fOutputTree->Fill();
		}

		if(!tracks.size()){
			// No track found , neutral particle? Look for EMCal clusters
			std::vector<AliVCluster *> clusters;
			FindClusters(ipart, clusters);
			if(clusters.size() > 1) fOutputInfo.isUnique = false;
			for(std::vector<AliVCluster *>::iterator clustIter = clusters.begin(); clustIter != clusters.end(); clustIter++){
				assoccluster = *clustIter;
				fOutputInfo.E = GetClusterEnergy(assoccluster, vpos, fOutputInfo.erec);
				GetShowerShape(assoccluster, fOutputInfo.showershape);
				GetCellEnergies(assoccluster, fOutputInfo.fCellEnergies, fOutputInfo.fCellIndices);
				fOutputTree->Fill();
			}
		}
	}

	PostData(2, fOutputTree);
	return kTRUE;
}

double AliEmcalMCTreeWriter::GetClusterEnergy(AliVCluster* clust, Double_t *vpos, Double_t* evec) const {
	/*
	 * Get energy from a cluster. Fill also as directed quantity
	 */
	TLorentzVector clustervector;
	clust->GetMomentum(clustervector, vpos);
	TVector3 vec = clustervector.Vect().Unit();
	vec *= clustervector.E();
	if(evec)
		vec.GetXYZ(evec);
	return clustervector.E();
}

void AliEmcalMCTreeWriter::FindTracks(unsigned int label, std::vector<AliVTrack*>& tracks) {
	/*
	 * Reverse order finding: Find all tracks matched to a particle
	 */
	tracks.clear();
	TIter trackIter(fTracks);
	AliVTrack *track;
	while((track = dynamic_cast<AliVTrack *>(trackIter()))){
		if(static_cast<unsigned int>(TMath::Abs(track->GetLabel())) == label){
			// we have found a matching track
			tracks.push_back(track);
		}
	}
}

void AliEmcalMCTreeWriter::FindClusters(unsigned int label, std::vector<AliVCluster*>& clusters) {
	/*
	 * Reverse order: Find all clusters matched to a particle
	 */
	clusters.clear();
	TIter clusterIter(fCaloClusters);
	AliVCluster *cluster;
	while((cluster = dynamic_cast<AliVCluster *>(clusterIter()))){
		TArrayI labels;
		labels.Set(cluster->GetNLabels(), cluster->GetLabels());
		for(int *labIter = labels.GetArray(); labIter < labels.GetArray() + labels.GetSize(); labIter++){
			if(static_cast<unsigned int>(TMath::Abs(*labIter)) == label){
				// we have found a matching track
				clusters.push_back(cluster);
				break;
			}
		}
	}
}

void AliEmcalMCTreeWriter::GetShowerShape(AliVCluster* clust, Double_t* vector) const {
	vector[0] = clust->GetM02();
	vector[1] = clust->GetM20();
}

void AliEmcalMCTreeWriter::GetCellEnergies(AliVCluster* clust, TVectorD& cells, TVectorF &indices) const {
	cells.ResizeTo(clust->GetNCells());
	indices.ResizeTo(clust->GetNCells());
	UShort_t *cellIndices = clust->GetCellsAbsId();
	for(int icell = 0; icell < clust->GetNCells(); icell++){
		indices[icell] = cellIndices[icell];
		cells[icell] = fInputEvent->GetEMCALCells()->GetCellAmplitude(cellIndices[icell]);
	}
}

bool AliEmcalMCTreeWriter::AcceptParticle(const AliVParticle* const particle) const {
	/*
	 * accept particle in case it is a photon, pi0, charged pion, kaon, proton or electron
	 */
	 const int kAcceptPdg[7] = {kGamma, kPi0, kPiPlus, kKPlus, kProton, kElectron};
	 bool result = false;
	 for(const int *pdgiter = kAcceptPdg; pdgiter < kAcceptPdg + 7; pdgiter++){
		 if(TMath::Abs(particle->PdgCode()) == *pdgiter){
			 result = true;
			 break;
		 }
	 }
	 return result;
}
