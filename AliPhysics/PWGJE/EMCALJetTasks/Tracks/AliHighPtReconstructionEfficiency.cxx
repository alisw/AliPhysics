/**************************************************************************
 * Copyright(c) 1998-2015, ALICE Experiment at CERN, All rights reserved. *
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
#include <iostream>
#include <sstream>

#include <TFile.h>
#include <TKey.h>
#include <TList.h>
#include <TMath.h>
#include <TProfile.h>
#include <TNtuple.h>
#include <TString.h>
#include <TSystem.h>
#include <TTree.h>

#include "AliAnalysisManager.h"
#include "AliAODMCParticle.h"
#include "AliAODTrack.h"
#include "AliESDtrackCuts.h"
#include "AliESDtrack.h"
#include "AliLog.h"
#include "AliMCEvent.h"
#include "AliMCParticle.h"
#include "AliTrackReference.h"
#include "AliVEvent.h"
#include "AliVParticle.h"
#include "AliVTrack.h"

#include <fastjet/PseudoJet.hh>
#include <fastjet/JetDefinition.hh>
#include <fastjet/ClusterSequence.hh>

#include "AliReducedJetConstituent.h"
#include "AliReducedJetEvent.h"
#include "AliReducedJetInfo.h"
#include "AliReducedJetParticle.h"
#include "AliReducedMatchedTrack.h"

#include "AliHighPtReconstructionEfficiency.h"

/// \cond CLASSIMP
ClassImp(HighPtTracks::AliHighPtReconstructionEfficiency)
/// \endcond

namespace HighPtTracks {

/**
 * Dummy (I/O) constructor. Not to be used
 */
AliHighPtReconstructionEfficiency::AliHighPtReconstructionEfficiency() :
	AliAnalysisTaskSE(),
	fParticleMap(NULL),
	fJetTree(NULL),
	fJetEvent(NULL),
	fMaxEtaJets(0.5),
	fMaxEtaParticles(0.8),
	fMinPtParticles(0.15),
	fMaxDR(0.5),
	fCrossSection(0.),
	fNtrials(0),
	fPtHardBin(0),
	fTaskDebugMode(false)
{
  memset(fTrackCuts, 0, sizeof(AliESDtrackCuts *) * 2);
}

/**
 * Main constructor, initialises relevant settings for the task. Setting default values for
 * kinematic track selection
 *
 * \param name Name of the task
 */
AliHighPtReconstructionEfficiency::AliHighPtReconstructionEfficiency(const char* name):
	AliAnalysisTaskSE(name),
	fParticleMap(NULL),
	fJetTree(NULL),
	fJetEvent(NULL),
	fMaxEtaJets(0.5),
	fMaxEtaParticles(0.8),
	fMinPtParticles(0.15),
	fMaxDR(0.5),
	fCrossSection(0.),
	fNtrials(0),
	fPtHardBin(0),
	fTaskDebugMode(false)
{
	DefineOutput(1, TNtuple::Class());
  memset(fTrackCuts, 0, sizeof(AliESDtrackCuts *) * 2);
}

/**
 * Destructor, cleaning up memory allocated by the task.
 */
AliHighPtReconstructionEfficiency::~AliHighPtReconstructionEfficiency() {
  for(int icut = 0; icut < 2; icut++)
    if(fTrackCuts[icut]) delete fTrackCuts[icut];
	if(fParticleMap) delete fParticleMap;
}

/**
 * Preparing output objects and track selection: create a new TTree and
 * reduced event object, and set the branch address to the new reduced event.
 */
void AliHighPtReconstructionEfficiency::UserCreateOutputObjects() {
	OpenFile(1);

	fJetTree = new TTree("JetEvent","A tree with jet information");
	fJetEvent = new AliReducedJetEvent();
	fJetTree->Branch("JetEvent", "AliReducedJetEvent", fJetEvent, 128000,0);

	PostData(1, fJetTree);
}

/**
 * Function informing the user when a new file starts. Used to
 *
 * \return True if the function is executed successfully, false in case of errors
 */
bool AliHighPtReconstructionEfficiency::UserNotify(){
	TTree *tree = AliAnalysisManager::GetAnalysisManager()->GetTree();
	if (!tree) {
		AliError(Form("%s - UserNotify: No current tree!",GetName()));
		return kFALSE;
	}

	TFile *curfile = tree->GetCurrentFile();
	if (!curfile) {
		AliError(Form("%s - UserNotify: No current file!",GetName()));
		return kFALSE;
	}

	PythiaInfoFromFile(curfile->GetName(), fCrossSection, fNtrials, fPtHardBin);

	return kTRUE;
}

/**
 * Runs event loop and creates reconstructed event:
 *  -# Filter particles used for the jet finding
 *  -# Filter MC-true particles and create lookup of particles with matching reconstructed particle,
 *     based on the MC label
 *  -# Find jets
 *  -# Create Output event structure with output jets, and match the pre-filtered particles to the jet,
 *     select only those with a distance smaller than the distance cut
 *  -# Write out the event into a TTree
 *
 * \param option Opitional settings, not used here
 */
void AliHighPtReconstructionEfficiency::UserExec(Option_t* /*option*/) {
	TList listforjets;
	SelectParticlesForJetfinding(listforjets);

	// Find Jets
	AliVParticle *part(NULL);
	TIter partIter(&listforjets);
	std::vector<fastjet::PseudoJet> inputjets;
	while((part = dynamic_cast<AliVParticle *>(partIter()))){
		fastjet::PseudoJet inputpart(part->Px(), part->Py(), part->Pz(), part->E());
		inputpart.set_user_index(part->GetLabel());
		inputjets.push_back(inputpart);
	}
	fastjet::ClusterSequence jetfinder(inputjets, fastjet::JetDefinition(fastjet::antikt_algorithm, 0.4));
	std::vector<fastjet::PseudoJet> recjets = fastjet::sorted_by_pt(jetfinder.inclusive_jets());
	CreateRectrackLookup();

	new(fJetEvent) AliReducedJetEvent(fCrossSection, fNtrials, fPtHardBin);

	std::vector<AliReconstructedParticlePair> selectedParticles = SelectParticles();

	if(fTaskDebugMode){
		std::stringstream debugmessage;
		debugmessage << "Number of particles: gen[" << selectedParticles.size() << "], rec[" << fParticleMap->GetNumberOfParticles() << "]";
		AliDebug(1, debugmessage.str().c_str());

		int nparticlesGen = 0, nparticlesRec = 0;
		for(std::vector<AliReconstructedParticlePair>::iterator pairiter = selectedParticles.begin(); pairiter != selectedParticles.end(); ++pairiter){
			if(pairiter->GetMCTrueParticle()) nparticlesGen++;
			if(pairiter->GetRecTracks().GetNumberOfParticles()) nparticlesRec++;
		}
		AliDebug(1, Form("Among %d selected particles we find %d reconstructed particles.", nparticlesGen, nparticlesRec));
	}

	AliReducedJetInfo *recjet(NULL);
	for(std::vector<fastjet::PseudoJet>::iterator jetiter = recjets.begin(); jetiter != recjets.end(); ++jetiter){
		if(TMath::Abs(jetiter->eta()) > fMaxEtaJets) continue;
		recjet = new AliReducedJetInfo(jetiter->px(), jetiter->py(), jetiter->pz(), jetiter->E());
		ProcessJet(recjet, selectedParticles);
		ConvertConstituents(recjet, *jetiter);
		fJetEvent->AddReconstructedJet(recjet);
	}
	fJetTree->Fill();
	PostData(1, fJetTree);
}

/**
 * Checks if particle is accepted by a give cut type. Currently implemented cut types are the \f$ R_{AA} \f$
 * standard cuts and the hybrid track cuts.
 *
 * \param track
 * \param cuttype
 * \return True if particle survived the cut, false otherwise
 */
bool AliHighPtReconstructionEfficiency::IsSelected(const AliVTrack* const track, CutType_t cuttype) const {
	/*
	 * Select reconstructed particle by applying track cuts
	 */
	if(!fTrackCuts[cuttype]) return kTRUE;
	const AliESDtrack *inputtrack = dynamic_cast<const AliESDtrack *>(track);
	if(inputtrack){
		return fTrackCuts[cuttype]->AcceptTrack(inputtrack);
	} else {
		// AOD track: Use copy
		AliESDtrack copytrack(track);
		return fTrackCuts[cuttype]->AcceptTrack(&copytrack);
	}
}

/**
 * Selects particle as MC-true particle. At this level, the particle \f$ p_{t} \f$ and \f$ \eta \f$
 * are checked. Only charged particles are accepted.
 *
 * \param track Particle to be checked
 * \return True if particle is accepted as true particle according to the truth defintion
 */
bool AliHighPtReconstructionEfficiency::IsTrueSelected(const AliVParticle* const track) const {
	if(TMath::Abs(track->Eta()) > fMaxEtaParticles) return false;
	if(TMath::Abs(track->Pt()) < fMinPtParticles) return false;
	if(!track->Charge()) return false;
	return true;
}

/**
 * Select particles at generation level which are the input for the jet finder. Particles are
 * requested to be within the kinematic acceptance and physical primary. Both charged and neutral
 * particles are accepted for jet finding.
 *
 * \param particles Which are used as input for the jet finder, filled by this function.
 */
void AliHighPtReconstructionEfficiency::SelectParticlesForJetfinding(TList& particles) const {
	AliVParticle *part = 0;
	particles.Clear();
	for(int ipart = 0; ipart < fMCEvent->GetNumberOfTracks(); ipart++){
		part = fMCEvent->GetTrack(ipart);
		if(!IsPhysicalPrimary(part)) continue;
		if(TMath::Abs(part->Eta()) > fMaxEtaParticles) continue;
		if(TMath::Abs(part->Pt()) < fMinPtParticles) continue;
		particles.Add(part);
	}
}

/**
 * This function associates reconstructed particles with Monte-Carlo true particles. For this,
 * a lookup tables which maps all reconstructed particles to a label is used internally. The
 * lookup table is created when the event starts and belongs to the event. Only tracks which
 * are physical primary signal tracks and in the kinematic acceptance are selected. The result
 * is a vector of particle pairs, which store a pointer to the Monte-Carlo particle and a list of
 * reconstructed tracks associated.
 *
 * \return Vector of reconstructed particle sets
 */
std::vector<AliReconstructedParticlePair> AliHighPtReconstructionEfficiency::SelectParticles() const {
	AliVParticle *part(NULL);
	AliParticleList tracks;
	std::vector<AliReconstructedParticlePair> result;
	for(int ipart = 0; ipart <fMCEvent->GetNumberOfTracks(); ipart++){
		part = fMCEvent->GetTrack(ipart);
		if(!IsPhysicalPrimary(part)) continue;
		if(!IsTrueSelected(part)) continue;
		const AliParticleList *tmptracks = FindReconstructedParticleFast(ipart);
		if(tmptracks) tracks = *tmptracks;
		if(fTaskDebugMode){
			std::stringstream debugmessage;
			debugmessage << "Accepted particle with label " << ipart;
			if(tracks.GetNumberOfParticles())
				debugmessage << " - found reconstructed track";
			else
				debugmessage << " - did not find reconstructed track";
			AliDebug(1, debugmessage.str().c_str());
		}

		AliReconstructedParticlePair foundpair;
		foundpair.SetMCTrueParticle(part);
		foundpair.SetRecParticles(tracks);
		result.push_back(foundpair);
	}
	return result;
}

/**
 * Calculates the distance between a particle and the main jet axis.
 *
 * \param recjet Jet where the particle is associated to
 * \param inputtrack Paricle to check
 * \return Distance between jet main axis and particle
 */
double AliHighPtReconstructionEfficiency::GetDR(const fastjet::PseudoJet& recjet, const AliVParticle* const inputtrack) const {
	return recjet.delta_R(fastjet::PseudoJet(inputtrack->Px(), inputtrack->Py(), inputtrack->Pz(), inputtrack->E()));
}

/**
 * Fill reduced jet information with associated particles at generator and reconstruction level. Only
 * particles within a maximum distance to the jet are selected to be added to the jet structure. Also
 * all particles reconstructed with a minimum requirement (TPC refit and ITS refit) are added to a particle
 * at this step.
 *
 * \param recjet Jet to be filled
 * \param particles List of particle pairs containing particle at generator level and associated reconstruced tracks
 */
void AliHighPtReconstructionEfficiency::ProcessJet(AliReducedJetInfo * const recjet,
    const std::vector<AliReconstructedParticlePair>& particles) const
{
	double pvec[4];
	recjet->GetPxPyPxE(pvec[0], pvec[1], pvec[2], pvec[3]);
	fastjet::PseudoJet frecjet(pvec[0], pvec[1], pvec[2], pvec[3]);		// For distance to the main jet axis
	AliVParticle *mcpart(NULL);
	const AliVTrack *rectrack(NULL);
	for(std::vector<AliReconstructedParticlePair>::const_iterator partiter = particles.begin(); partiter != particles.end(); ++partiter){
		mcpart = partiter->GetMCTrueParticle();
		const AliParticleList &matchedTracks = partiter->GetRecTracks();
		double dr = GetDR(frecjet, mcpart);
		if(dr < fMaxDR){
			// Create new Particle and add it to the jet reconstructed jet
			AliReducedJetParticle * part = new AliReducedJetParticle(mcpart->Px(), mcpart->Py(), mcpart->Pz(), mcpart->E(), mcpart->PdgCode());
			part->SetNumberOfTPCtrackReferences(GetNumberOfTPCTrackReferences(mcpart));
			part->SetDistanceToMainJetAxis(dr);
			for(int itrk = 0; itrk < matchedTracks.GetNumberOfParticles(); itrk++){
			  rectrack = matchedTracks.GetParticle(itrk);
			  AliReducedMatchedTrack *matchedTrack = new AliReducedMatchedTrack(rectrack->Px(), rectrack->Py(), rectrack->Pz());
			  matchedTrack->SetNumberOfClustersTPC(rectrack->GetTPCNcls());
				if(rectrack->GetLabel() >= 0) matchedTrack->SetGoodTrackLabel(true);
				// FilterTrackCuts
				if(IsSelected(rectrack, kRJStandardCuts)) matchedTrack->SetSurvivedTrackCuts(AliReducedMatchedTrack::kHighPtStandard);
				if(IsSelected(rectrack, kRJHybridCuts)) matchedTrack->SetSurvivedTrackCuts(AliReducedMatchedTrack::kHighPtHybrid);
				part->AddMatchedTrack(matchedTrack);
			}
			recjet->AddParticleInCone(part);
		}
	}
}

/**
 * Creates a lookup table where particles are associated to MC labels and stores it in the event.
 */
void AliHighPtReconstructionEfficiency::CreateRectrackLookup() {
	if(fParticleMap) delete fParticleMap;
	fParticleMap = new AliParticleMap;
	for(int ipart = 0; ipart < fInputEvent->GetNumberOfTracks(); ipart++){
		AliVTrack *mytrack = dynamic_cast<AliVTrack *>(fInputEvent->GetTrack(ipart));
		// Only apply a minimum selection: TPC refit and ITS refit
		// check particle cuts on tree level
		if(!((mytrack->GetStatus() & AliVTrack::kITSrefit) && (mytrack->GetStatus() & AliVTrack::kTPCrefit))) continue;
		//if(!IsSelected(mytrack)) continue;
		fParticleMap->AddParticle(mytrack);
	}
}

/**
 * \brief Find reconstructed particles for a given label
 *
 * Find all reconstructed particles associated with a reconstruced track. Uses a lookup table which is
 * created at the start of the event.
 *
 * \param label Label to check
 * \return List of reconstructed particles (NULL if label is not found)
 */
const AliParticleList * AliHighPtReconstructionEfficiency::FindReconstructedParticleFast(int label) const {
  return fParticleMap->GetParticles(label);
}

/**
 * \brief Checks if the particle is a physical primary particle
 *
 * Function checks whether particle is a physical primary particle. Handles ESD and AOD MC particles.
 *
 * \param part Particle to check.
 * \return True if the particle is a physical primary particle, false otherwise
 */
bool AliHighPtReconstructionEfficiency::IsPhysicalPrimary(const AliVParticle* const part) const {
	const AliAODMCParticle *aodpart = dynamic_cast<const AliAODMCParticle *>(part);
	if(aodpart) return aodpart->IsPhysicalPrimary();	// AOD MC particle
	else if (dynamic_cast<const AliMCParticle *>(part)) return fMCEvent->IsPhysicalPrimary(part->GetLabel());		// ESD MC particle
	return false;
}

/**
 * Convert jet constituents at jet finder level to the reduced format and adds it to the reduced
 * reconstructed jet.
 *
 * \param recjet Reconstructed jet the consituent is associated to
 * \param inputjet A jet constituent to be added to the reconstructed jet
 */
void AliHighPtReconstructionEfficiency::ConvertConstituents(AliReducedJetInfo* const recjet, const fastjet::PseudoJet& inputjet) {
	std::vector<fastjet::PseudoJet> constituents = inputjet.constituents();
	for(std::vector<fastjet::PseudoJet>::const_iterator citer = constituents.begin(); citer != constituents.end(); ++citer){
		AliVParticle *mcpart = fMCEvent->GetTrack(citer->user_index());
		AliReducedJetConstituent *jetconst = new AliReducedJetConstituent(citer->px(), citer->py(), citer->px(), citer->E(), mcpart->PdgCode());
		recjet->AddConstituent(jetconst);
	}
}

/**
 * Extract cross section and number of trials from a root file created at generation level
 * and fills it to the fieds added in the parameter list. Function is called once a new file is
 * loaded via the function UserNotify. From AliAnalysisTaskEmcal.
 *
 * \param currFile File with PYTHIA hard informatino
 * \param fXsec Output field for cross section
 * \param fTrials Output field for number of trials
 * \param pthard Output field for generated \f$ p_{t} \f$ of the hard interaction
 * \return Success status (false in case of any error)
 */
bool AliHighPtReconstructionEfficiency::PythiaInfoFromFile(const char* currFile, double &fXsec, double &fTrials, int &pthard) const {
	TString file(currFile);
	fXsec = 0;
	fTrials = 1;

	if (file.Contains(".zip#")) {
		Ssiz_t pos1 = file.Index("root_archive",12,0,TString::kExact);
		Ssiz_t pos = file.Index("#",1,pos1,TString::kExact);
		Ssiz_t pos2 = file.Index(".root",5,TString::kExact);
		file.Replace(pos+1,pos2-pos1,"");
	} else {
		// not an archive take the basename....
		file.ReplaceAll(gSystem->BaseName(file.Data()),"");
	}

	// Get the pt hard bin
	TString strPthard(file);

	strPthard.Remove(strPthard.Last('/'));
	strPthard.Remove(strPthard.Last('/'));
	if (strPthard.Contains("AOD")) strPthard.Remove(strPthard.Last('/'));
	strPthard.Remove(0,strPthard.Last('/')+1);
	if (strPthard.IsDec())
		pthard = strPthard.Atoi();

	// problem that we cannot really test the existance of a file in a archive so we have to live with open error message from root
	TFile *fxsec = TFile::Open(Form("%s%s",file.Data(),"pyxsec.root"));

	if (!fxsec) {
		// next trial fetch the histgram file
		fxsec = TFile::Open(Form("%s%s",file.Data(),"pyxsec_hists.root"));
		if (!fxsec) {
			// not a severe condition but inciate that we have no information
			return kFALSE;
		} else {
			// find the tlist we want to be independtent of the name so use the Tkey
			TKey* key = (TKey*)fxsec->GetListOfKeys()->At(0);
			if (!key) {
				fxsec->Close();
				return kFALSE;
			}
			TList *list = dynamic_cast<TList*>(key->ReadObj());
			if (!list) {
				fxsec->Close();
				return kFALSE;
			}
			fXsec = ((TProfile*)list->FindObject("h1Xsec"))->GetBinContent(1);
			fTrials = ((TH1F*)list->FindObject("h1Trials"))->GetBinContent(1);
			fxsec->Close();
		}
	} else { // no tree pyxsec.root
		TTree *xtree = (TTree*)fxsec->Get("Xsection");
		if (!xtree) {
			fxsec->Close();
			return kFALSE;
		}
		UInt_t ntrials = 0;
		Double_t xsection = 0;
		xtree->SetBranchAddress("xsection",&xsection);
		xtree->SetBranchAddress("ntrials",&ntrials);
		xtree->GetEntry(0);
		fTrials = ntrials;
		fXsec = xsection;
		fxsec->Close();
	}
	return kTRUE;
 }

/**
 * Get the number of track references for a MC-true track
 *
 * \param trk the MC particle to check
 * \return The number of track references
 */
unsigned short AliHighPtReconstructionEfficiency::GetNumberOfTPCTrackReferences(AliVParticle* const trk) const {
  unsigned short nref = 0;
  if(trk->IsA() == AliMCParticle::Class()){
    AliMCParticle *part = dynamic_cast<AliMCParticle *>(trk);
    AliTrackReference *tref(NULL);
    for(int iref = 0; iref < part->GetNumberOfTrackReferences(); ++iref){
      tref = part->GetTrackReference(iref);
      if(tref->DetectorId() == AliTrackReference::kTPC) nref++;
    }
  }
  return nref;
}

}
