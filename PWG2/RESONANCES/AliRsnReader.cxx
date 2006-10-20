/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
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
 
//-------------------------------------------------------------------------
//                     Class AliRsnReader
//
//   Reader for conversion of ESD or Kinematics output into AliRsnEvent
//   .....
//   .....
//   .....
//   .....
//   .....
//   .....
//   .....
// 
// author: A. Pulvirenti             (email: alberto.pulvirenti@ct.infn.it)
//-------------------------------------------------------------------------

#include <Riostream.h>

#include <TH1.h>
#include <TH3.h>
#include <TFile.h>
#include <TTree.h>
#include <TArrayF.h>
#include <TParticle.h>
#include <TRandom.h>

#include "AliRun.h"
#include "AliESD.h"
#include "AliStack.h"
#include "AliHeader.h"
#include "AliESDtrack.h"
#include "AliRunLoader.h"
#include "AliGenPythiaEventHeader.h"

#include "AliRsnDaughter.h"
#include "AliRsnEvent.h"
#include "AliRsnReader.h"

ClassImp(AliRsnReader)

//--------------------------------------------------------------------------------------------------------
AliRsnReader::AliRsnReader()
{
//
// Constructor.
// Initializes all working parameters to default values:
// - ESD particle identification
// - rejection of non-ITS-refitted tracks
// - maximum distance allowed from primary vertex = 3.0 cm (beam pipe)
//
	fPIDMethod = kESDPID;
	Int_t i;
	for (i = 0; i < AliPID::kSPECIES; i++) fPrior[i] = 1.0;
	fPtLimit4PID = 4.0;
	fProbThreshold = 0.0;
	fMaxRadius = 3.0;      // the beam pipe

	fEvents = 0;
	for (i = AliPID::kElectron; i <= AliPID::kProton; i++) {
		fEffPos[i] = 0;
		fEffNeg[i] = 0;
	}
	fResPt = 0;
	fResLambda = 0;
	fResPhi = 0;
}
//--------------------------------------------------------------------------------------------------------
AliRsnReader::AliRsnReader(const AliRsnReader &copy) : TObject(copy)
{
//
// Copy constructor.
// Initializes all working parameters to same values of another simlar object.
// TTree data member is not created.
//
	fPIDMethod = copy.fPIDMethod;
	Int_t i;
	for (i = 0; i < AliPID::kSPECIES; i++) fPrior[i] = copy.fPrior[i];
	fPtLimit4PID = copy.fPtLimit4PID;
	fProbThreshold = copy.fProbThreshold;
	fMaxRadius = copy.fMaxRadius;

	fEvents = 0;
	for (i = AliPID::kElectron; i <= AliPID::kProton; i++) {
		fEffPos[i] = copy.fEffPos[i];
		fEffNeg[i] = copy.fEffNeg[i];
	}
	
	fResPt = copy.fResPt;
	fResLambda = copy.fResLambda;
	fResPhi = copy.fResPhi;
}
//--------------------------------------------------------------------------------------------------------
AliRsnReader& AliRsnReader::operator=(const AliRsnReader &copy)
{
//
// Assignment operator.
// Initializes all working parameters to same values of another simlar object.
// TTree data member is not created.
//
	fPIDMethod = copy.fPIDMethod;
	Int_t i;
	for (i = 0; i < AliPID::kSPECIES; i++) fPrior[i] = copy.fPrior[i];
	fPtLimit4PID = copy.fPtLimit4PID;
	fProbThreshold = copy.fProbThreshold;
	fMaxRadius = copy.fMaxRadius;

	fEvents = 0;
	for (i = AliPID::kElectron; i <= AliPID::kProton; i++) {
		fEffPos[i] = copy.fEffPos[i];
		fEffNeg[i] = copy.fEffNeg[i];
	}
	
	fResPt = copy.fResPt;
	fResLambda = copy.fResLambda;
	fResPhi = copy.fResPhi;
	
	return (*this);
}
//--------------------------------------------------------------------------------------------------------
void AliRsnReader::Clear(Option_t *option)
{
//
// Clear collection of filenames.
// If requested with the option "DELETELIST", 
// the collection object is also deleted.
//
	TString opt(option);
	
	if (!opt.CompareTo("TREE", TString::kIgnoreCase)) {
		fEvents->Reset();
		if (!opt.CompareTo("DELTREE", TString::kIgnoreCase)) {
			delete fEvents;
			fEvents = 0;
		}
	}
	
	Int_t i;
	for (i = AliPID::kElectron; i <= AliPID::kProton; i++) {
		fEffPos[i] = 0;
		fEffNeg[i] = 0;
	}
	fResPt = 0;
	fResLambda = 0;
	fResPhi = 0;
}
//--------------------------------------------------------------------------------------------------------
Bool_t AliRsnReader::EffSim(Int_t pdg, Double_t pt, Double_t eta, Double_t z)
{
//
// If efficiency histogram are present, they are used to simulate efficiency.
// Pt, Eta and Z are used to find reconstruction efficiency value to be used
// and PDG is used to select efficiency histogram for a given particle.
// An extraction is done, and it is supposed that particle must be accepted
// only when this function retunrs kTRUE (= successful extraction).
//
	// find particle sign from PDG code
	Int_t sign;
	if (TMath::Abs(pdg) >= 211) {
		if (pdg > 0) sign = 1; else sign = -1;
	}
	else {
		if (pdg > 0) sign = -1; else sign = 1;
	}
	
	// convert PDG code into one value in AliPID::kSPECIES
	// (if returned value is not a charged particle, procedure is skipped)
	Int_t index = FindType(pdg);
	TH3D *eff = 0;
	if (index >= AliPID::kElectron && index <= AliPID::kProton) {
		if (sign > 0) eff = fEffPos[index]; else eff = fEffNeg[index];
	}
	
	// if no efficiency histogram is defined, method returns a 'fail' result
	if (!eff) return kFALSE;
	
	// otherwise, a random extraction is made
	Int_t ibin = eff->FindBin(z, eta, pt);
	Double_t ref = (Double_t)eff->GetBinContent(ibin);
	Double_t ran = gRandom->Rndm();
	return (ran <= ref);
}
//--------------------------------------------------------------------------------------------------------
Double_t* AliRsnReader::GetPIDprobabilities(AliRsnDaughter track) const
{
//
// Computes PID probabilites starting from priors and weights
//
	Double_t *prob = new Double_t[AliPID::kSPECIES];
	
	Int_t i;
	
	// step 1 - compute the normalization factor
	Double_t sum = 0.0;
	for (i = 0; i < AliPID::kSPECIES; i++) {
		prob[i] = fPrior[i] * track.GetPIDweight(i);
		sum += prob[i];
	}
	if (sum <= 0.0) return 0;
	
	// step 2 - normalize PID weights by the relative prior probability
	for (Int_t i = 0; i < AliPID::kSPECIES; i++) {
		prob[i] /= sum;
	}
	
	return prob;
}
//--------------------------------------------------------------------------------------------------------
void AliRsnReader::Identify(AliRsnDaughter &track)
{
//
// Identifies a track according to method selected
//
	Bool_t doESDPID = (fPIDMethod == kESDPID);
	Bool_t keepRecSign = (fPIDMethod == kPerfectPIDwithRecSign);
	
	if (doESDPID) {
		// when doing ESD PID it is supposed that charge sign
		// comes from ESD track and it is not modified
		Double_t pt = track.GetPt();
		if (pt <= fPtLimit4PID) {
			Double_t *prob = GetPIDprobabilities(track);
			if (!prob) track.SetPDG(0);
			Int_t idx[AliPID::kSPECIES];
			TMath::Sort(AliPID::kSPECIES, prob, idx);
			Int_t imax = idx[0];
			Double_t maxprob = prob[imax];
			if (maxprob >= fProbThreshold) {
				track.SetPDG((UShort_t)AliPID::ParticleCode(imax));
			}
			delete [] prob;
		}
		else {
			track.SetPDG(0);
		}
	}
	else {
		Short_t truePDG = track.GetTruePDG();
		track.SetPDG((UShort_t)TMath::Abs(truePDG));
		if (!keepRecSign) {
			if (TMath::Abs(truePDG) <= 13) {
				if (truePDG > 0) track.SetSign(-1); else track.SetSign(1);
			}
			else {
				if (truePDG > 0) track.SetSign(1); else track.SetSign(-1);
			}
		}
	}
}
//--------------------------------------------------------------------------------------------------------
TTree* AliRsnReader::ReadParticles(const char *path, Option_t *option)
{
//
// Opens the file "galice.root" and kinematics in the path specified,
// loads Kinematics informations and fills and AliRsnEvent object
// with data coming from simulation.
// If required, an efficiency simulation can be done which cuts off some tracks
// depending on a MonteCarlo extraction weighted on the efficiency histograms
// passed to the class. 
// Moreover, a smearing can be done if requested, using the resolution histograms
// passed to the class.
// Allowed options:
// - "E" --> do efficiency simulation
// - "P" --> do momentum smearing
//
	fPIDMethod = kPerfectPID;

	TTree *events = new TTree("selection", "AliRsnEvents");
	TTree::SetBranchStyle(1);
	AliRsnEvent *event = new AliRsnEvent;
	TBranch *branch = events->Branch("events", "AliRsnEvent", &event, 32000, 1);
	branch->SetAutoDelete(kFALSE);
	
	// get path
	TString strPath(path);
	if (strPath.Length() < 1) return 0;
	strPath += '/';
	
	// evaluate options
	TString opt(option);
	opt.ToUpper();
	Bool_t doEffSim = (opt.Contains("E"));
	Bool_t doSmearing = (opt.Contains("P"));
	
	// initialize run loader
	AliRunLoader *runLoader = OpenRunLoader(path);
	if (!runLoader) return 0;
	Int_t nEvents = runLoader->GetNumberOfEvents();
	TArrayF vertex(3);
	
	// loop on events
	Int_t i, iev, index = 0, imum, nParticles, nStoredParticles = 0;
	Double_t vtot, px, py, pz, arg, eta;
	TParticle *particle = 0, *mum = 0;
	for (iev = 0; iev < nEvents; iev++) {
		
		// load given event
		cout << "\rEvent " << iev << " of " << nEvents << flush;
		runLoader->GetEvent(iev);
		
		// primary vertex
		runLoader->GetHeader()->GenEventHeader()->PrimaryVertex(vertex);
		//cout << "Process type: " << ((AliGenPythiaEventHeader*)runLoader->GetHeader()->GenEventHeader())->ProcessType() << endl;
		
		// kinematics
		AliStack *stack = runLoader->Stack();
		if (!stack) continue;
		
		// get number of particles to collect
		nParticles = stack->GetNtrack();
		if (!nParticles) return 0;
				
		// create new AliRsnEvent object
		event->Clear("DELETE");
		event->SetESD(kFALSE);
		event->SetPath(strPath);
		event->Init();
		
		// loop on tracks
		index = 0;
		for (i = 0; i < nParticles; i++) {
			
			// get particle
			particle = stack->Particle(i);
			if (!particle) continue;
			
			// check against maximum allowed distance from primary vertex
			vtot  = (particle->Vx() - vertex[0])*(particle->Vx() - vertex[0]);
			vtot += (particle->Vy() - vertex[1])*(particle->Vy() - vertex[1]);
			vtot += (particle->Vz() - vertex[2])*(particle->Vz() - vertex[2]);
			vtot  = TMath::Sqrt(vtot);
			if (vtot > fMaxRadius) continue;
			
			// efficiency selection
			if (doEffSim) {
				arg = 0.5 * TMath::ATan2(particle->Pz(), particle->Pt());
				if (arg > 0.) eta = -TMath::Log(arg); else continue;
				Bool_t result = EffSim(particle->GetPdgCode(), (Double_t)vertex[2], eta, (Double_t)particle->Pt());
				if (result == kFALSE) continue;
			}
			
			// smear kinematic parameters (if requested)
			px = particle->Px();
			py = particle->Py();
			pz = particle->Pz();
			if (doSmearing) SmearMomentum(px, py, pz);
									
			// add to collection of tracks
			nStoredParticles++;
			AliRsnDaughter *track = new AliRsnDaughter;
			if (!track->Adopt(particle)) {
				delete track;
				continue;
			}
			track->SetIndex(index);
			track->SetLabel(i);
			track->SetPxPyPz(px, py, pz);
			imum = particle->GetFirstMother();
			if (imum >= 0) {
				mum = stack->Particle(imum);
				track->SetMotherPDG( (Short_t)mum->GetPdgCode() );
			}
			Identify(*track);
			
			event->AddTrack(*track);
			delete track;
			
		} // end of loop over tracks
		
		// compute total multiplicity
		event->GetMultiplicity();

		// link to events tree and fill
		events->Fill();
	}
	
	runLoader->UnloadKinematics();
	delete runLoader;
	
	return events;
}
//--------------------------------------------------------------------------------------------------------
TTree* AliRsnReader::ReadTracks(const char *path, Option_t *option)
{
//
// Reads AliESD in a given path, with and "experimental" method.
// When using this method, no Kinematics information is assumed
// and particle identification is always done with the bayesian method.
// No Kinematics informations are stored.
// Allowed options (actually):
// - "R" : reject tracks which do not have the flag "kITSRefit" set to TRUE
// - "F" : reject 'fake' tracks (negative label)
//
	fPIDMethod = kESDPID;

	TTree *events = new TTree("selection", "AliRsnEvents");
	TTree::SetBranchStyle(1);
	AliRsnEvent *event = new AliRsnEvent;
	TBranch *branch = events->Branch("events", "AliRsnEvent", &event, 32000, 1);
	branch->SetAutoDelete(kFALSE);
	
	// get path
	TString strPath(path);
	if (strPath.Length() < 1) return 0;
	strPath += '/';
	
	// evaluate options
	TString opt(option);
	opt.ToUpper();
	Bool_t checkITSrefit = (opt.Contains("R"));
	Bool_t rejectFakes = (opt.Contains("F"));
	
	// opens ESD file
	TFile *fileESD = TFile::Open(Form("%s/AliESDs.root", strPath.Data()));
	if (!fileESD) return 0;
	if (fileESD->IsOpen() == kFALSE) return 0;
	TTree* treeESD = (TTree*)fileESD->Get("esdTree");
	AliESD *esd = 0;
	treeESD->SetBranchAddress("ESD", &esd);
	Int_t nev = (Int_t)treeESD->GetEntries();
	
	// loop on events
	Int_t i, nSelTracks = 0; 
	Double_t vertex[3];
	for (i = 0; i < nev; i++) {
	
		// message
		cout << "\rEvent " << i << flush;
		treeESD->GetEntry(i);
		
		// primary vertex
		vertex[0] = esd->GetVertex()->GetXv();
		vertex[1] = esd->GetVertex()->GetYv();
		vertex[2] = esd->GetVertex()->GetZv();
		
		// create new AliRsnEvent object
		event->Clear("DELETE");
		event->Init();
		event->SetPath(strPath);
		event->SetESD();
				
		// get number of tracks
		Int_t ntracks = esd->GetNumberOfTracks();
		if (!ntracks) continue;
		
		// store tracks from ESD
		Int_t index, label;
		Double_t vtot, v[3];
		for (index = 0; index < ntracks; index++) {
			
			// get track
			AliESDtrack *esdTrack = esd->GetTrack(index);
			
			// check against vertex constraint
			esdTrack->GetXYZ(v);
			vtot  = (v[0] - vertex[0])*(v[0] - vertex[0]);
			vtot += (v[1] - vertex[1])*(v[1] - vertex[1]);
			vtot += (v[2] - vertex[2])*(v[2] - vertex[2]);
			vtot  = TMath::Sqrt(vtot);
			if (vtot > fMaxRadius) continue;
			
			// check for ITS refit
			if (checkITSrefit) {
				if (!(esdTrack->GetStatus() & AliESDtrack::kITSrefit)) continue;
			}
			
			// check for fakes
			label = esdTrack->GetLabel();
			if (rejectFakes && (label < 0)) continue;
			
			// create AliRsnDaughter (and make Bayesian PID)
			AliRsnDaughter track;
			if (!track.Adopt(esdTrack, checkITSrefit)) continue;
			track.SetIndex(index);
			track.SetLabel(label);
			Identify(track);
			
			// store in TClonesArray
			event->AddTrack(track);
			nSelTracks++;
		}
		
		// compute total multiplicity
		event->GetMultiplicity();
	
		// link to events tree and fill
		events->Fill();
	}
	
	fileESD->Close();
	return events;
}
//--------------------------------------------------------------------------------------------------------
TTree* AliRsnReader::ReadTracksAndParticles(const char *path, Option_t *option)
{
//
// Reads AliESD in a given path, getting also informations from Kinematics.
// In this case, the PID method used is the one selected with apposite setter.
// Allowed options (actually):
// - "R" : reject tracks which do not have the flag "kITSRefit" set to TRUE
// - "F" : reject 'fake' tracks (negative label)
// - "M" : use 'true' momentum instead of reconstructed one
//
	TTree *events = new TTree("selection", "AliRsnEvents");
	TTree::SetBranchStyle(1);
	AliRsnEvent *event = new AliRsnEvent;
	TBranch *branch = events->Branch("events", "AliRsnEvent", &event, 32000, 1);
	branch->SetAutoDelete(kFALSE);
	
	// get path
	TString strPath(path);
	if (strPath.Length() < 1) return 0;
	
	// evaluate options
	TString opt(option);
	opt.ToUpper();
	Bool_t checkITSrefit = (opt.Contains("R"));
	Bool_t rejectFakes = (opt.Contains("F"));
	Bool_t copyMomentum = (opt.Contains("M"));
	
	// opens ESD file
	TFile *fileESD = TFile::Open(Form("%s/AliESDs.root", strPath.Data()));
	if (!fileESD) return 0;
	if (fileESD->IsOpen() == kFALSE) return 0;
	TTree* treeESD = (TTree*)fileESD->Get("esdTree");
	AliESD *esd = 0;
	treeESD->SetBranchAddress("ESD", &esd);
	Int_t nevRec = (Int_t)treeESD->GetEntries();
	
	// initialize run loader
	AliRunLoader *runLoader = OpenRunLoader(path);
	if (!runLoader) return 0;
	Int_t nevSim = runLoader->GetNumberOfEvents();
	
	// check number of reconstructed and simulated events
	if ( (nevSim != 0 && nevRec != 0) && (nevSim != nevRec) ) {
		cerr << "Count mismatch: sim = " << nevSim << " -- rec = " << nevRec << endl;
		return 0;
	}
	else if (nevSim == 0 && nevRec == 0) {
		cerr << "Count error: sim = rec = 0" << endl;
		return 0;
	}
	
	// loop on events
	Int_t i, procType, ntracks, nSelTracks = 0; 
	Double_t vertex[3];
	for (i = 0; i < nevRec; i++) {
	
		// get event
		cout << "\rEvent " << i << " " << flush;
		treeESD->GetEntry(i);
		runLoader->GetEvent(i);
				
		// reject event if it is diffractive
		procType = ((AliGenPythiaEventHeader*)runLoader->GetHeader()->GenEventHeader())->ProcessType();
		if (procType == 92 || procType == 93 || procType == 94) {
			cout << "Skipping diffractive event" << endl;
			continue;
		}
		
		// get particle stack
		AliStack *stack = runLoader->Stack();
				
		// primary vertex
		vertex[0] = esd->GetVertex()->GetXv();
		vertex[1] = esd->GetVertex()->GetYv();
		vertex[2] = esd->GetVertex()->GetZv();
		
		// multiplicity
		ntracks = esd->GetNumberOfTracks();
		if (!ntracks) {
			Warning("ReadTracksAndParticles", "Event %d has no tracks!", i);
			continue;
		}
		
		// create new AliRsnEvent object
		event->Clear("DELETE");
		event->Init();
		event->SetPath(strPath);
		event->SetESD();
				
		// store tracks from ESD
		Int_t index, label;
		Double_t vtot, v[3];
		for (index = 0; index < ntracks; index++) {
			
			// get track
			AliESDtrack *esdTrack = esd->GetTrack(index);
			
			// check against vertex constraint
			esdTrack->GetXYZ(v);
			vtot  = (v[0] - vertex[0])*(v[0] - vertex[0]);
			vtot += (v[1] - vertex[1])*(v[1] - vertex[1]);
			vtot += (v[2] - vertex[2])*(v[2] - vertex[2]);
			vtot  = TMath::Sqrt(vtot);
			if (vtot > fMaxRadius) continue;
			
			// check for fakes
			label = esdTrack->GetLabel();
			if (rejectFakes) {
				if (label < 0) continue;
			}
			
			// create AliRsnDaughter (and make Bayesian PID)
			AliRsnDaughter track;
			if (!track.Adopt(esdTrack, checkITSrefit)) continue;
			track.SetIndex(index);
			
			// retrieve particle and get Kine info
			TParticle *part = stack->Particle(TMath::Abs(label));
			track.SetTruePDG(part->GetPdgCode());
			Int_t mother = part->GetFirstMother();
			track.SetMother(mother);
			if (mother >= 0) {
				TParticle *mum = stack->Particle(mother);
				track.SetMotherPDG(mum->GetPdgCode());
			}
			if (copyMomentum) track.SetPxPyPz(part->Px(), part->Py(), part->Pz());
			
			// identification
			Identify(track);
			
			// store in TClonesArray
			event->AddTrack(track);
			nSelTracks++;
		}
		
		// compute total multiplicity
		event->GetMultiplicity();
	
		// link to events tree and fill
		events->Fill();
	}
	
	runLoader->UnloadKinematics();
	delete runLoader;
	fileESD->Close();
	
	return events;
}
//--------------------------------------------------------------------------------------------------------
void AliRsnReader::SetEfficiencyHistogram(AliPID::EParticleType type, TH3D *h, Bool_t pos)
{
//
// Sets efficiency histogram for a given particle species.
// If third argument is "true", histo is assigned to positive particles,
// otherwise it is assigned to negative particles.
//
	if (type >= AliPID::kElectron && type < AliPID::kPhoton) {
		if (pos) fEffPos[type] = h; else fEffNeg[type] = h;
	}
}
//--------------------------------------------------------------------------------------------------------
void AliRsnReader::SetPriorProbabilities(Double_t *prior)
{
//
// Set prior probabilities to be used in case of ESD PID.
//
	if (!prior) return;
	
	Int_t i = 0;
	for (i = 0; i < AliPID::kSPECIES; i++) fPrior[i] = prior[i];
}
//--------------------------------------------------------------------------------------------------------
void AliRsnReader::SetPriorProbability(AliPID::EParticleType type, Double_t value)
{
//
// Sets prior probability referred to a single particle species.
//
	if (type >= AliPID::kElectron && type < AliPID::kPhoton) fPrior[type] = value;
}
//--------------------------------------------------------------------------------------------------------
void AliRsnReader::SmearMomentum(Double_t &px, Double_t &py, Double_t &pz) 
{
//
// Use resolution histograms to do smearing of momentum
// (to introduce reconstruction effects)
//
	Double_t pt = TMath::Sqrt(px*px + py*py);
	Double_t lambda = TMath::ATan2(pz, pt);
	Double_t phi = TMath::ATan2(py, px);
	
	Int_t ibin;
	Double_t ref;
	if (fResPt) {
		ibin = fResPt->FindBin(pt);
		ref = (Double_t)fResPt->GetBinContent(ibin);
		pt *= 1.0 + ref;
	}
	if (fResLambda) {
		ibin = fResLambda->FindBin(pt);
		ref = (Double_t)fResLambda->GetBinContent(ibin);
		lambda *= 1.0 + ref;
	}
	if (fResPhi) {
		ibin = fResPhi->FindBin(pt);
		ref = (Double_t)fResPhi->GetBinContent(ibin);
		phi *= 1.0 + ref;
	}
	
	px = pt * TMath::Cos(phi);
	py = pt * TMath::Sin(phi);
	pz = pt * TMath::Tan(lambda);
}
//--------------------------------------------------------------------------------------------------------
AliPID::EParticleType AliRsnReader::FindType(Int_t pdg)
{
//
// Finds the enum value corresponding to a PDG code
//
	pdg = TMath::Abs(pdg);
	switch (pdg) {
		case   11: return AliPID::kElectron; break;
		case   13: return AliPID::kMuon; break;
		case  211: return AliPID::kPion; break;
		case  321: return AliPID::kKaon; break;
		case 2212: return AliPID::kProton; break;
		default  : return AliPID::kPhoton;
	}
}
//--------------------------------------------------------------------------------------------------------
AliRunLoader* AliRsnReader::OpenRunLoader(const char *path)
{
//
// Open the Run loader with events in a given path
//
	// clear gALICE
	if (gAlice) {
		delete gAlice;
		gAlice = 0;
	}
	
	// initialize run loader
	TString name(path);
	name += "/galice.root";
	AliRunLoader *runLoader = AliRunLoader::Open(name.Data());
	if (runLoader) {
		runLoader->LoadgAlice();
		gAlice = runLoader->GetAliRun();
		runLoader->LoadKinematics();
		runLoader->LoadHeader();
	}
	
	return runLoader;
}
