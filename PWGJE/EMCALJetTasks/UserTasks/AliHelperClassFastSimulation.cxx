#include <vector>

#include <TParticle.h>
#include <TArrayI.h>
#include <TRandom3.h>
#include <TObjArray.h>
#include <TF1.h>

#include "TDatime.h"

#include "AliEmcalJet.h"
#include "FJ_includes.h"
#include "AliFJWrapper.h"
#include "AliAODMCParticle.h"
#include "AliMCParticle.h"

#include "AliAnalysisTaskMTFPID.h"
#include "AliHelperClassFastSimulation.h"

ClassImp(AliHelperClassFastSimulation)

AliHelperClassFastSimulation::AliHelperClassFastSimulation(TF1** effFunctions, Double_t simEffFactor, Double_t simResFactor, Double_t simRes, Double_t jetMinPt, Double_t jetMaxEta, Double_t jetMinEta, AliFJWrapper* wrapper) {
	fEffFunctions = new TF1*[2*AliPID::kSPECIES];
	for (Int_t i=0;i<2*AliPID::kSPECIES;i++) {;
		fEffFunctions[i] = (TF1*)effFunctions[i]->Clone();
	}
	fFastSimEffFactor = simEffFactor;
	fFastSimResFactor = simResFactor;
	fFastSimRes = simRes;
	fJetMinPt = jetMinPt;
	fJetMaxEta = jetMaxEta;
	fJetMinEta = jetMinEta;
	
	if (wrapper)
		fWrapper = wrapper;
	else
		fWrapper = GetStandardJetWrapper();
	
	fRandom = new TRandom3();
  TDatime* date = new TDatime();
  fRandom->SetSeed(date->Get());
  
  delete date;
	
	fJetArray = new TObjArray();
	fParticlesArray = new TObjArray(50);
	fPDGCodeArray = new TArrayI(50);
}

AliHelperClassFastSimulation::~AliHelperClassFastSimulation() {
	try {
		delete fWrapper;
		fWrapper = 0x0;
	}
	catch (...) {
	}
	fJetArray->SetOwner(kTRUE);
	delete fJetArray;
	fJetArray = 0x0;
	fParticlesArray->SetOwner(kTRUE);
	delete fParticlesArray;
	fParticlesArray = 0x0;
	delete fEffFunctions;
	fEffFunctions = 0x0;
	delete fPDGCodeArray;
	fPDGCodeArray = 0x0;
}

void AliHelperClassFastSimulation::AddParticle(AliAODMCParticle* part) {
	AddInputVector(part->Pt(), part->Phi(), part->Eta(), part->M(), part->GetLabel(), part->Charge(), part->GetPdgCode());
}

void AliHelperClassFastSimulation::AddInputVector(Double_t pT, Double_t phi, Double_t eta, Double_t mass, Int_t index, Double_t charge, Int_t pdgCode) {
	Int_t mcID = AliAnalysisTaskMTFPID::PDGtoMCID(pdgCode);
	
	if (mcID == AliPID::kMuon || mcID >= AliPID::kSPECIES)
		return;
	
	Bool_t posCharge = charge > 0.0;

	Double_t xeff = fFastSimEffFactor * fEffFunctions[2*mcID+(Int_t)posCharge]->Eval(pT);
	Double_t x = fRandom->Rndm();
	
	if (x > xeff)
		return;

	Double_t smearedPt = 1.0/(fRandom->Gaus(1.0/pT,fFastSimResFactor * fFastSimRes)); 

	Double_t px = smearedPt * TMath::Cos(phi);
	Double_t py = smearedPt * TMath::Sin(phi);
	Double_t pz = smearedPt * TMath::SinH(eta);


	fWrapper->AddInputVector(px, py, pz, TMath::Sqrt(smearedPt * smearedPt + pz * pz + mass*mass), index);
	if (fPDGCodeArray->GetSize() < index)
		fPDGCodeArray->Set(index * 2);
		
	fPDGCodeArray->AddAt(pdgCode, index);
}

void AliHelperClassFastSimulation::Run() {
	fWrapper->Run();
	std::vector<fastjet::PseudoJet> jets = fWrapper->GetInclusiveJets();
	
	Int_t acceptedJetNumber = 0;
	for (UInt_t j=0;j<jets.size();++j) {
		if (jets[j].perp() < fJetMinPt || jets[j].eta() > fJetMaxEta || jets[j].eta() < fJetMinEta)
			continue;
		
		AliEmcalJet *jet = new AliEmcalJet(jets[j].perp(), jets[j].eta(), jets[j].phi(), jets[j].m());
		
		jet->SetLabel(j);
		fastjet::PseudoJet area(fWrapper->GetJetAreaVector(j));
		jet->SetArea(area.perp());
		jet->SetAreaEta(area.eta());
		jet->SetAreaPhi(area.phi());
		jet->SetAreaE(area.E());
//       jet->SetJetAcceptanceType(FindJetAcceptanceType(jet->Eta(), jet->Phi_0_2pi(), fRadius)); 
		
		std::vector<fastjet::PseudoJet> constituents(fWrapper->GetJetConstituents(j));
		
		TObjArray* partArray = new TObjArray();
		Int_t nOfParticles = 0;
		
		for (UInt_t ic = 0;ic<constituents.size();++ic) {
			Int_t uid = constituents[ic].user_index();
			
			if (uid == -1)    //Ghost particle
				continue;
			
			partArray->AddAtAndExpand(GetParticleFromConstituent(constituents[ic]), nOfParticles);
			nOfParticles++;
		}
		fParticlesArray->AddAtAndExpand(partArray, acceptedJetNumber);
		
		jet->SetNumberOfTracks(nOfParticles);
		fJetArray->AddAtAndExpand(jet, acceptedJetNumber);
		acceptedJetNumber++;
	}
}

AliEmcalJet* AliHelperClassFastSimulation::GetJet(Int_t nOfJet) {
	return (AliEmcalJet*)fJetArray->At(nOfJet);
}

AliAODMCParticle* AliHelperClassFastSimulation::GetTrackOfJet(Int_t nOfTrack, Int_t nOfJet) {
	return (AliAODMCParticle*)((TObjArray*)fParticlesArray->At(nOfJet))->At(nOfTrack);
}

AliAODMCParticle* AliHelperClassFastSimulation::GetParticleFromConstituent(fastjet::PseudoJet constituent) {
	Int_t label = constituent.user_index();
	
	TParticle* tpart = new TParticle(fPDGCodeArray->At(label), 0, 0, 0, 0, 0, constituent.px(), constituent.py(), constituent.pz(), constituent.E(), 0.0, 0.0, 0.0, 0.0);
	AliMCParticle* mcParticle = new AliMCParticle(tpart);
	AliAODMCParticle* part = new AliAODMCParticle(mcParticle, label);
	return part;
}

AliFJWrapper* AliHelperClassFastSimulation::GetStandardJetWrapper() {
	AliFJWrapper* wrapper = new AliFJWrapper("wrapper", "wrapper");
	wrapper->SetAreaType(fastjet::active_area_explicit_ghosts);
  wrapper->SetGhostArea(0.005);
  wrapper->SetR(0.4);
  //Currently not working, including AliEmcalJetTask fails
//         wrapper->SetAlgorithm(AliEmcalJetTask::ConvertToFJAlgo(jetContainer->GetJetAlgorithm()));
//         wrapper->SetRecombScheme(AliEmcalJetTask::ConvertToFJRecoScheme(jetContainer->GetRecombinationScheme()));
  wrapper->SetAlgorithm(fastjet::antikt_algorithm);
  wrapper->SetRecombScheme(fastjet::pt_scheme);
  wrapper->SetMaxRap(1);
	return wrapper;
}

UInt_t AliHelperClassFastSimulation::GetNJets() {
	if (fJetArray)
		return fJetArray->GetEntriesFast();
	else
		return 0;
}

UInt_t AliHelperClassFastSimulation::GetNParticlesOfJet(UInt_t jetNumber) {
	if (fParticlesArray && jetNumber <= GetNJets() - 1)
		return ((TObjArray*)fParticlesArray->At(jetNumber))->GetEntriesFast();
	else 
		return 0.0;
}
