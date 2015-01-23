/**************************************************************************
 * Copyright(c) 1998-2014, ALICE Experiment at CERN, All rights reserved. *
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

#include "AliHighPtReconstructionEfficiency.h"

ClassImp(AliHighPtReconstructionEfficiency)

AliHighPtReconstructionEfficiency::AliHighPtReconstructionEfficiency() :
	AliAnalysisTaskSE(),
	fTrackCuts(NULL),
	fResults(NULL),
	fParticleMap(NULL),
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
	fResults(NULL),
	fParticleMap(NULL),
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
	OpenFile(1);

	fResults = new TNtuple("ParticlesInCone", "Particles related to a jet cone", "px:py:pz:energy:pdgcodee:conedist:jetpt:isreconstructed:crosssection:trials:pthardbin");

	// Create std track cuts
	fTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(true, 1);
	fTrackCuts->SetName("Standard Track cuts");
	fTrackCuts->SetMinNCrossedRowsTPC(120);
	fTrackCuts->SetMaxDCAToVertexXYPtDep("0.0182+0.0350/pt^1.01");

	PostData(1, fResults);
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
	TList listforjets;
	SelectParticlesForJetfinding(listforjets);
	std::vector<fastjet::PseudoJet> recjets = FindJets(listforjets);
	CreateRectrackLookup();

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

	for(std::vector<fastjet::PseudoJet>::iterator jetiter = recjets.begin(); jetiter != recjets.end(); ++jetiter){
		if(TMath::Abs(jetiter->eta()) > fMaxEtaJets) continue;
		ProcessJet(*jetiter, selectedParticles);
	}
	PostData(1, fResults);
}

bool AliHighPtReconstructionEfficiency::IsSelected(const AliVTrack* const track) const {
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
	if(TMath::Abs(track->Eta()) > fMaxEtaParticles) return false;
	if(TMath::Abs(track->Pt()) < fMinPtParticles) return false;
	if(!track->Charge()) return false;
	return true;
}

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

std::vector<fastjet::PseudoJet> AliHighPtReconstructionEfficiency::FindJets(const TList& inputparticles) const {
	AliVParticle *part(NULL);
	TIter partIter(&inputparticles);
	std::vector<fastjet::PseudoJet> inputjets;
	while((part = dynamic_cast<AliVParticle *>(partIter()))){
		inputjets.push_back(fastjet::PseudoJet(part->Px(), part->Py(), part->Pz(), part->E()));
	}
	fastjet::ClusterSequence jetfinder(inputjets, fastjet::JetDefinition(fastjet::antikt_algorithm, 0.4));
	return fastjet::sorted_by_pt(jetfinder.inclusive_jets());
}

std::vector<AliReconstructedParticlePair> AliHighPtReconstructionEfficiency::SelectParticles() const {
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
	return recjet.delta_R(fastjet::PseudoJet(inputtrack->Px(), inputtrack->Py(), inputtrack->Pz(), inputtrack->E()));
}

AliVTrack* AliHighPtReconstructionEfficiency::FindReconstructedParticle(int label) const {
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

void AliHighPtReconstructionEfficiency::ProcessJet(const fastjet::PseudoJet& recjet, const std::vector<AliReconstructedParticlePair>& particles) {
	const AliVParticle *mcpart(NULL);
	const AliVTrack *rectrack(NULL);
	for(std::vector<AliReconstructedParticlePair>::const_iterator partiter = particles.begin(); partiter != particles.end(); ++partiter){
		mcpart = partiter->GetMCTrueParticle();
		rectrack = partiter->GetRecTrack();
		double dr = GetDR(recjet, mcpart);
		if(dr < fMaxDR)
			fResults->Fill(mcpart->Px(), mcpart->Py(), mcpart->Pz(), mcpart->E(), mcpart->PdgCode(),dr,recjet.pt(), rectrack ? 1 : 0, fCrossSection, fNtrials, fPtHardBin);
	}
}

void AliHighPtReconstructionEfficiency::CreateRectrackLookup() {
	if(fParticleMap) delete fParticleMap;
	fParticleMap = new AliParticleMap;
	for(int ipart = 0; ipart < fInputEvent->GetNumberOfTracks(); ipart++){
		AliVTrack *mytrack = dynamic_cast<AliVTrack *>(fInputEvent->GetTrack(ipart));
		if(!IsSelected(mytrack)) continue;
		fParticleMap->AddParticle(mytrack);
	}
}

AliVTrack* AliHighPtReconstructionEfficiency::FindReconstructedParticleFast(int label) const {
	AliParticleList *particles = fParticleMap->GetParticles(label);
	if(!particles) return NULL;
	return particles->GetParticle(0);			// Always return the first occurrence
}

bool AliHighPtReconstructionEfficiency::IsPhysicalPrimary(const AliVParticle* const part) const {
	const AliAODMCParticle *aodpart = dynamic_cast<const AliAODMCParticle *>(part);
	if(aodpart) return aodpart->IsPhysicalPrimary();	// AOD MC particle
	else if (dynamic_cast<const AliMCParticle *>(part)) return fMCEvent->IsPhysicalPrimary(part->GetLabel());		// ESD MC particle
	return false;
}

bool AliHighPtReconstructionEfficiency::PythiaInfoFromFile(const char* currFile, double &fXsec, double &fTrials, int &pthard) const {
	//
	// From AliAnalysisTaskEmcal
	//

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

AliParticleMap::~AliParticleMap() {
	// Clean up all particle lists
	for(std::map<int, AliParticleList *>::iterator it = fParticles.begin(); it != fParticles.end(); ++it){
		delete it->second;
	}
}

void AliParticleMap::AddParticle(AliVTrack *track){
	int label = TMath::Abs(track->GetLabel());
	std::map<int, AliParticleList *>::iterator it = fParticles.find(label);
	if(it == fParticles.end()){			// not existing
		AliParticleList *nextparticle = new AliParticleList;
		nextparticle->AddParticle(track);
		fParticles.insert(std::pair<int, AliParticleList *>(label, nextparticle));
	} else {
		AliParticleList *mylist = it->second;
		mylist->AddParticle(track);
	}
}

AliParticleList* AliParticleMap::GetParticles(int label) const {
	AliParticleList *result = NULL, *content = NULL;
	std::map<int, AliParticleList *>::const_iterator found = fParticles.find(label);
	if(found != fParticles.end()){
		result = found->second;
	}
	return result;
}

void AliParticleMap::Print() const {
	for(std::map<int, AliParticleList *>::const_iterator it = fParticles.begin(); it != fParticles.end(); ++it){
		std::cout << "Particle with label " << it->first << ", number of reconstructed assigned: " << it->second->GetNumberOfParticles() << std::endl;
	}
}
