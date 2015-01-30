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
/*
 * Tree creation for tracking efficiency studies as function of the distance
 * to the main jet axis
 *
 * Author:
 *   Markus Fasel <markus.fasel@cern.ch>
 */
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

#include "AliHighPtReconstructionEfficiency.h"

ClassImp(HighPtTracks::AliHighPtReconstructionEfficiency)

namespace HighPtTracks {

AliHighPtReconstructionEfficiency::AliHighPtReconstructionEfficiency() :
	AliAnalysisTaskSE(),
	fTrackCuts(NULL),
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
}

AliHighPtReconstructionEfficiency::AliHighPtReconstructionEfficiency(const char* name):
	AliAnalysisTaskSE(name),
	fTrackCuts(NULL),
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
}

AliHighPtReconstructionEfficiency::~AliHighPtReconstructionEfficiency() {
	if(fTrackCuts) delete fTrackCuts;
	if(fParticleMap) delete fParticleMap;
}

void AliHighPtReconstructionEfficiency::UserCreateOutputObjects() {
	/*
	 * Prepare output objects and track selection
	 */
	OpenFile(1);

	fJetTree = new TTree("JetEvent","A tree with jet information");
	fJetEvent = new AliReducedJetEvent();
	fJetTree->Branch("JetEvent", "AliReducedJetEvent", fJetEvent, 128000,0);

	PostData(1, fJetTree);
}

bool AliHighPtReconstructionEfficiency::UserNotify(){
	// Read Pythia cross section, incorporated from AliAnalysisTaskEmcal

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

void AliHighPtReconstructionEfficiency::UserExec(Option_t*) {
	/*
	 * Event loop:
	 * 1. Filter particles used for the jet finding
	 * 2. Filter MC-true particles and create lookup of particles with matching reconstructed particle,
	 *    based on the MC label
	 * 3. Find jets
	 * 4. Create Output event structure with output jets, and match the pre-filtered particles to the jet,
	 *    select only those with a distance smaller than the distance cut
	 * 5. Write out the event into a TTree
	 */
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
			if(pairiter->GetRecTrack()) nparticlesRec++;
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

bool AliHighPtReconstructionEfficiency::IsSelected(const AliVTrack* const track) const {
	/*
	 * Select reconstructed particle by applying track cuts
	 */
	const AliESDtrack *inputtrack = dynamic_cast<const AliESDtrack *>(track);
	if(inputtrack){
		return fTrackCuts->AcceptTrack(inputtrack);
	} else {
		// AOD track: Use copy
		AliESDtrack copytrack(track);
		return fTrackCuts->AcceptTrack(&copytrack);
	}
}

bool AliHighPtReconstructionEfficiency::IsTrueSelected(const AliVParticle* const track) const {
	/*
	 * Accept particle as MC-true particle
	 */
	if(TMath::Abs(track->Eta()) > fMaxEtaParticles) return false;
	if(TMath::Abs(track->Pt()) < fMinPtParticles) return false;
	if(!track->Charge()) return false;
	return true;
}

void AliHighPtReconstructionEfficiency::SelectParticlesForJetfinding(TList& particles) const {
	/*
	 * Pre-fileter particles for the jet reconstruction
	 */
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

std::vector<AliReconstructedParticlePair> AliHighPtReconstructionEfficiency::SelectParticles() const {
	/*
	 * Select particles to be matched with the reconstructed jet
	 */
	const AliVParticle *part(NULL);
	AliVTrack *track(NULL);
	std::vector<AliReconstructedParticlePair> result;
	for(int ipart = 0; ipart <fMCEvent->GetNumberOfTracks(); ipart++){
		part = fMCEvent->GetTrack(ipart);
		if(!IsPhysicalPrimary(part)) continue;
		if(!IsTrueSelected(part)) continue;
		track = FindReconstructedParticleFast(ipart);
		if(fTaskDebugMode){
			std::stringstream debugmessage;
			debugmessage << "Accepted particle with label " << ipart;
			if(track)
				debugmessage << " - found reconstructed track";
			else
				debugmessage << " - did not find reconstructed track";
			AliDebug(1, debugmessage.str().c_str());
		}

		AliReconstructedParticlePair foundpair;
		foundpair.SetMCTrueParticle(part);
		foundpair.SetRecParticle(track);
		result.push_back(foundpair);
	}
	return result;
}

double AliHighPtReconstructionEfficiency::GetDR(const fastjet::PseudoJet& recjet, const AliVParticle* const inputtrack) const {
	/*
	 * Calculte distance to the main jet axis
	 */
	return recjet.delta_R(fastjet::PseudoJet(inputtrack->Px(), inputtrack->Py(), inputtrack->Pz(), inputtrack->E()));
}

AliVTrack* AliHighPtReconstructionEfficiency::FindReconstructedParticle(int label) const {
	/*
	 * Find matching reconstructed particle for a given MC label (slow version)
	 */
	AliVTrack *selected(NULL), *tmp(NULL);
	for(int itrk = 0; itrk < fInputEvent->GetNumberOfTracks(); itrk++){
		tmp = dynamic_cast<AliVTrack *>(fInputEvent->GetTrack(itrk));
		if(TMath::Abs(tmp->GetLabel()) == label){
			if(IsSelected(tmp)){
				selected = tmp;
				break;
			}
		}
	}
	return selected;
}

void AliHighPtReconstructionEfficiency::ProcessJet(
		AliReducedJetInfo * const recjet,
		const std::vector<AliReconstructedParticlePair>& particles
		) const
{
	/*
	 * Select particles with radius to the jet less than the maximum allowed radius and adds them to the reconstructed jet info
	 */
	double pvec[4];
	recjet->GetPxPyPxE(pvec[0], pvec[1], pvec[2], pvec[3]);
	fastjet::PseudoJet frecjet(pvec[0], pvec[1], pvec[2], pvec[3]);		// For distance to the main jet axis
	const AliVParticle *mcpart(NULL);
	const AliVTrack *rectrack(NULL);
	for(std::vector<AliReconstructedParticlePair>::const_iterator partiter = particles.begin(); partiter != particles.end(); ++partiter){
		mcpart = partiter->GetMCTrueParticle();
		rectrack = partiter->GetRecTrack();
		double dr = GetDR(frecjet, mcpart);
		if(dr < fMaxDR){
			// Create new Particle and add it to the jet reconstructed jet
			AliReducedJetParticle * part = new AliReducedJetParticle(mcpart->Px(), mcpart->Py(), mcpart->Pz(), mcpart->E(), mcpart->PdgCode(),rectrack ? true : false);
			part->SetDistanceToMainJetAxis(dr);
			if(rectrack){
				part->SetDeltaPt(TMath::Abs(mcpart->Pt()) - TMath::Abs(rectrack->Pt()));
				part->SetNumberOfClustersTPC(static_cast<unsigned char>(rectrack->GetTPCNcls()));
			}
			recjet->AddParticleInCone(part);
		}
	}
}

void AliHighPtReconstructionEfficiency::CreateRectrackLookup() {
	/*
	 * Create lookup table of particles and reconstructed tracks according to the MC label
	 */
	if(fParticleMap) delete fParticleMap;
	fParticleMap = new AliParticleMap;
	for(int ipart = 0; ipart < fInputEvent->GetNumberOfTracks(); ipart++){
		AliVTrack *mytrack = dynamic_cast<AliVTrack *>(fInputEvent->GetTrack(ipart));
		if(!IsSelected(mytrack)) continue;
		fParticleMap->AddParticle(mytrack);
	}
}

AliVTrack* AliHighPtReconstructionEfficiency::FindReconstructedParticleFast(int label) const {
	/*
	 * Find reconstructed particle for a given MC label using the lookup table (fast version)
	 */
	AliParticleList *particles = fParticleMap->GetParticles(label);
	if(!particles) return NULL;
	return particles->GetParticle(0);			// Always return the first occurrence
}

bool AliHighPtReconstructionEfficiency::IsPhysicalPrimary(const AliVParticle* const part) const {
	/*
	 * Method checking whether particle is physical primary, blind to ESD and AOD particles
	 */
	const AliAODMCParticle *aodpart = dynamic_cast<const AliAODMCParticle *>(part);
	if(aodpart) return aodpart->IsPhysicalPrimary();	// AOD MC particle
	else if (dynamic_cast<const AliMCParticle *>(part)) return fMCEvent->IsPhysicalPrimary(part->GetLabel());		// ESD MC particle
	return false;
}

void AliHighPtReconstructionEfficiency::ConvertConstituents(
		AliReducedJetInfo* const recjet, const fastjet::PseudoJet& inputjet) {
	/*
	 * Convert jet constituents into reduced format
	 */
	std::vector<fastjet::PseudoJet> constituents = inputjet.constituents();
	for(std::vector<fastjet::PseudoJet>::const_iterator citer = constituents.begin(); citer != constituents.end(); ++citer){
		AliVParticle *mcpart = fMCEvent->GetTrack(citer->user_index());
		AliReducedJetConstituent *jetconst = new AliReducedJetConstituent(citer->px(), citer->py(), citer->px(), citer->E(), mcpart->PdgCode());
		recjet->AddConstituent(jetconst);
	}
}

bool AliHighPtReconstructionEfficiency::PythiaInfoFromFile(const char* currFile, double &fXsec, double &fTrials, int &pthard) const {
	/*
	 * From AliAnalysisTaskEmcal: Get Pythia Hard information
	 */

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

}
