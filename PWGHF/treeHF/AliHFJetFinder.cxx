/* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//*************************************************************************
// \class AliHFJetFinder
// \helper class to handle jet finding, matching and substructure
// \authors:
// N. Zardoshti, nima.zardoshti@cern.ch
/////////////////////////////////////////////////////////////

#include <cmath>
#include <limits>
#include "AliHFJetFinder.h"
#include "TMath.h"
#include "TRandom3.h"

/// \cond CLASSIMP
ClassImp(AliHFJetFinder);
/// \endcond


//________________________________________________________________
AliHFJetFinder::AliHFJetFinder():
  fMinJetPt(0.0),
  fJetRadius(0.4),
  fJetAlgorithm(JetAlgorithm::antikt),
  fJetRecombScheme(RecombScheme::e_scheme),
  fJetGhostArea(0.005),
  fJetAreaType(AreaType::active),
  fMinSubJetPt(0.0),
  fSubJetRadius(0.0),
  fSubJetAlgorithm(JetAlgorithm::ca),
  fSubJetRecombScheme(RecombScheme::e_scheme),
  fSoftDropZCut(0.1),
  fSoftDropBeta(0.0),
  fMinTrackPt(0.15),
  fMaxTrackPt(100.0),
  fMinTrackE(0.0),
  fMaxTrackE(1000.0),
  fMaxTrackEta(0.9),
  fMaxTrackPhi(10.0),
  fMinParticlePt(0.0),
  fMaxParticlePt(1000.0),
  fMaxParticleEta(0.9),
  fCharged(Charge::charged),
  fTrackingEfficiency(1.0),
  fDoJetSubstructure(false),
  fFastJetWrapper(0x0)
{
  //
  // Default constructor
  //
}

//________________________________________________________________
AliHFJetFinder::AliHFJetFinder(char *name):
  fMinJetPt(0.0),
  fJetRadius(0.4),
  fJetAlgorithm(JetAlgorithm::antikt),
  fJetRecombScheme(RecombScheme::e_scheme),
  fJetGhostArea(0.005),
  fJetAreaType(AreaType::active),
  fMinSubJetPt(0.0),
  fSubJetRadius(0.0),
  fSubJetAlgorithm(JetAlgorithm::ca),
  fSubJetRecombScheme(RecombScheme::e_scheme),
  fSoftDropZCut(0.1),
  fSoftDropBeta(0.0),
  fMinTrackPt(0.15),
  fMaxTrackPt(100.0),
  fMinTrackE(0.0),
  fMaxTrackE(1000.0),
  fMaxTrackEta(0.9),
  fMaxTrackPhi(0.9),
  fMinParticlePt(0.0),
  fMaxParticlePt(1000.0),
  fMaxParticleEta(0.9),
  fCharged(Charge::charged),
  fTrackingEfficiency(1.0),
  fDoJetSubstructure(false),
  fFastJetWrapper(0x0)
{
}


//________________________________________________________________
AliHFJetFinder::~AliHFJetFinder()
{
  //
  // Destructor
  //
  delete fFastJetWrapper;
}

//________________________________________________________________
//Set the jet finding parameters
void AliHFJetFinder::SetFJWrapper() 
{

  
  fFastJetWrapper = new AliFJWrapper("fFastJetWrapper","fFastJetWrapper");

  fFastJetWrapper->Clear();

  fFastJetWrapper->SetR(fJetRadius); 
  fFastJetWrapper->SetAlgorithm(JetAlgorithm(fJetAlgorithm));
  fFastJetWrapper->SetRecombScheme(RecombinationScheme(fJetRecombScheme));
  fFastJetWrapper->SetGhostArea(fJetGhostArea); 
  fFastJetWrapper->SetAreaType(AreaType(fJetAreaType)); 
}


//________________________________________________________________
//returns jet clustered with heavy flavour candidate
AliHFJet AliHFJetFinder::GetHFJet(TClonesArray *array, AliAODRecoDecayHF *cand, Double_t invmass){ 

  SetFJWrapper();
  AliHFJet hfjet;
  if (!cand) return hfjet;
  FindJets(array,cand, invmass);
  Int_t jet_index=Find_Candidate_Jet();
  if (jet_index==-1) return hfjet;
 
  std::vector<fastjet::PseudoJet> inclusive_jets = fFastJetWrapper->GetInclusiveJets();
  fastjet::PseudoJet jet = inclusive_jets[jet_index];
  if (jet.perp() < fMinJetPt) return hfjet;
  std::vector<fastjet::PseudoJet> constituents(fFastJetWrapper->GetJetConstituents(jet_index));

  SetJetVariables(hfjet, constituents, jet, 0, cand); 

  return hfjet;
}


//________________________________________________________________
//returns jet clustered with heavy flavour particle (MC)
AliHFJet AliHFJetFinder::GetHFMCJet(TClonesArray *array, AliAODMCParticle *mcpart){
  SetFJWrapper();
  AliHFJet hfjet;
  if (!mcpart) return hfjet;
  FindMCJets(array,mcpart);
  Int_t jet_index=Find_Candidate_Jet();
  if (jet_index==-1) return hfjet;

 
  std::vector<fastjet::PseudoJet> inclusive_jets = fFastJetWrapper->GetInclusiveJets();
  fastjet::PseudoJet jet = inclusive_jets[jet_index];
  if (jet.perp() < fMinJetPt) return hfjet;
  std::vector<fastjet::PseudoJet> constituents(fFastJetWrapper->GetJetConstituents(jet_index));

  SetMCJetVariables(hfjet, constituents, jet, 0, mcpart);
 
  return hfjet;
}


//________________________________________________________________
//returns vector of jets, including the jet with the heavy flavour candidate
std::vector<AliHFJet> AliHFJetFinder::GetHFJets(TClonesArray *array, AliAODRecoDecayHF *cand, Double_t invmass) {
  
  SetFJWrapper();
  std::vector<AliHFJet> hfjet_vec;
  hfjet_vec.clear();
  if (!cand) return hfjet_vec;
  FindJets(array, cand ,invmass);
  Int_t jet_index=Find_Candidate_Jet();
  if (jet_index==-1) return hfjet_vec;

  std::vector<fastjet::PseudoJet> inclusive_jets = fFastJetWrapper->GetInclusiveJets();
  
  for (Int_t i=0; i<inclusive_jets.size(); i++){
    fastjet::PseudoJet jet = inclusive_jets[i];
    if (jet.perp() < fMinJetPt) continue;
    std::vector<fastjet::PseudoJet> constituents(fFastJetWrapper->GetJetConstituents(i));

    AliHFJet hfjet;

    if (i==jet_index) SetJetVariables(hfjet, constituents, jet, i, cand); 
    else SetJetVariables(hfjet, constituents, jet, i, nullptr); 
    hfjet_vec.push_back(hfjet);
  }
  return hfjet_vec;
}


//________________________________________________________________
//returns vector of jets, including the jet with the heavy flavour particle (MC)
std::vector<AliHFJet> AliHFJetFinder::GetHFMCJets(TClonesArray *array, AliAODMCParticle *mcpart) {
  SetFJWrapper();
  std::vector<AliHFJet> hfjet_vec;
  hfjet_vec.clear();
  if (!mcpart) return hfjet_vec;
  FindMCJets(array, mcpart);
  Int_t jet_index=Find_Candidate_Jet();
  if (jet_index==-1) return hfjet_vec;

  std::vector<fastjet::PseudoJet> inclusive_jets = fFastJetWrapper->GetInclusiveJets();
  
  for (Int_t i=0; i<inclusive_jets.size(); i++){
    fastjet::PseudoJet jet = inclusive_jets[i];
    if (jet.perp() < fMinJetPt) continue;
    std::vector<fastjet::PseudoJet> constituents(fFastJetWrapper->GetJetConstituents(i));

    AliHFJet hfjet;
  
    if (i==jet_index) SetMCJetVariables(hfjet, constituents, jet, i, mcpart); 
    else SetMCJetVariables(hfjet, constituents, jet, i, nullptr);
    
    hfjet_vec.push_back(hfjet);
  }
  return hfjet_vec;
}



//________________________________________________________________
//returns vector of jets
std::vector<AliHFJet> AliHFJetFinder::GetJets(TClonesArray *array) {

  //Jet is clustered with heavy flavour meson and the corresponding variables are set
  SetFJWrapper();
  std::vector<AliHFJet> hfjet_vec;
  hfjet_vec.clear();
  FindJets(array,nullptr);

  std::vector<fastjet::PseudoJet> inclusive_jets = fFastJetWrapper->GetInclusiveJets();
  for (Int_t i=0; i<inclusive_jets.size(); i++){
    fastjet::PseudoJet jet = inclusive_jets[i];
    if (jet.perp() < fMinJetPt) continue;
    std::vector<fastjet::PseudoJet> constituents(fFastJetWrapper->GetJetConstituents(i));

    AliHFJet hfjet;

    SetJetVariables(hfjet, constituents, jet, i, nullptr);
  
    hfjet_vec.push_back(hfjet);
  }
  return hfjet_vec;
}


//________________________________________________________________
//returns vector of jets (MC)
std::vector<AliHFJet> AliHFJetFinder::GetMCJets(TClonesArray *array) {
  //Jet is clustered with heavy flavour meson and the corresponding variables are set
  SetFJWrapper();
  std::vector<AliHFJet> hfjet_vec;
  hfjet_vec.clear();
  FindMCJets(array,nullptr);

  std::vector<fastjet::PseudoJet> inclusive_jets = fFastJetWrapper->GetInclusiveJets();
  for (Int_t i=0; i<inclusive_jets.size(); i++){
    fastjet::PseudoJet jet = inclusive_jets[i];
    if (jet.perp() < fMinJetPt) continue;
    std::vector<fastjet::PseudoJet> constituents(fFastJetWrapper->GetJetConstituents(i));

    AliHFJet hfjet;

    SetMCJetVariables(hfjet, constituents, jet, i, nullptr);
  
    hfjet_vec.push_back(hfjet);
  }
  return hfjet_vec;
}



//________________________________________________________________
//Do jet finding, including replacing the daughters of the heavy flavour candidate with the candidiate itself, if needed
void AliHFJetFinder::FindJets(TClonesArray *array, AliAODRecoDecayHF *cand, Double_t invmass) {
  //Performs jet finding. Jets are stored in the fFastJetWrapper object

  std::vector<Int_t> daughter_vec;
  if (cand){
  
    daughter_vec.clear();

    AliVTrack *daughter;
    for (Int_t i = 0; i < cand->GetNDaughters(); i++) {
      daughter = dynamic_cast<AliVTrack *>(cand->GetDaughter(i));   
      if (!daughter) continue;
      daughter_vec.push_back(daughter->GetID());
    }
    AliTLorentzVector cand_lvec(0,0,0,0);
    cand_lvec.SetPtEtaPhiM(cand->Pt(), cand->Eta(), cand->Phi(), invmass);
    fFastJetWrapper->AddInputVector(cand_lvec.Px(), cand_lvec.Py(), cand_lvec.Pz(), cand_lvec.E(),0); 
  }

    
  bool isdaughter;
  AliAODTrack *track=NULL;
  for (Int_t i=0; i<array->GetEntriesFast(); i++) {
 
    track= dynamic_cast<AliAODTrack*>(array->At(i));
    if(!CheckTrack(track)) continue; 
    isdaughter=false;
    if (cand){
      for (Int_t j=0; j<daughter_vec.size(); j++){
	if (track->GetID()==daughter_vec[j]) isdaughter=true;
      }
      if(isdaughter) continue;
    }
    
    fFastJetWrapper->AddInputVector(track->Px(), track->Py(), track->Pz(), track->E(),i+100);
  }
  fFastJetWrapper->Run();
  //delete track;
}



//________________________________________________________________
//Do jet finding, including replacing the daughters of the heavy flavour particle with the cparticle itself, if needed (MC)
void AliHFJetFinder::FindMCJets(TClonesArray *array, AliAODMCParticle *mcpart) {
  //Performs jet finding. Jets are stored in the fFastJetWrapper object

  std::vector<Int_t> daughter_vec;
  if (mcpart){
  
    daughter_vec.clear();

    AliAODMCParticle *daughter;
    for (Int_t i = 0; i < mcpart->GetNDaughters(); i++) {
      daughter = dynamic_cast<AliAODMCParticle *>(array->At(mcpart->GetDaughterLabel(i)));   
      if (!daughter) continue;
      daughter_vec.push_back(daughter->GetLabel());
    }
    fFastJetWrapper->AddInputVector(mcpart->Px(), mcpart->Py(), mcpart->Pz(), mcpart->E(),0); 
  }

    
  bool isdaughter;
  AliAODMCParticle *particle=NULL;
  for (Int_t i=0; i<array->GetEntriesFast(); i++) {
 
    particle= dynamic_cast<AliAODMCParticle*>(array->At(i));
    if(!CheckParticle(particle)) continue; 
    isdaughter=false;
    if (mcpart){
      for (Int_t j=0; j<daughter_vec.size(); j++){
	if (particle->GetLabel()==daughter_vec[j]) isdaughter=true;
      }
      if(isdaughter) continue;
    }
    fFastJetWrapper->AddInputVector(particle->Px(), particle->Py(), particle->Pz(), particle->E(),i+100); 
  }
  fFastJetWrapper->Run();
  //delete track;
}


//________________________________________________________________
//Set the jet parameters in the AliHFJet object
void AliHFJetFinder::SetJetVariables(AliHFJet& hfjet, const std::vector<fastjet::PseudoJet>& constituents, const fastjet::PseudoJet& jet, Int_t jetID, AliAODRecoDecayHF *cand) {

  hfjet.fID=jetID;
  if (cand)hfjet.fHFMeson=1;
  else hfjet.fHFMeson=0;
  hfjet.fPt=jet.perp();
  hfjet.fEta=jet.pseudorapidity();
  hfjet.fPhi=jet.phi();
  if (cand) hfjet.fDeltaEta=jet.pseudorapidity()-cand->Eta();
  else hfjet.fDeltaEta=-99;
  if (cand) hfjet.fDeltaPhi=RelativePhi(jet.phi(),cand->Phi());
  else hfjet.fDeltaPhi=-99;
  if (cand) hfjet.fDeltaR=TMath::Sqrt(hfjet.fDeltaEta*hfjet.fDeltaEta + hfjet.fDeltaPhi*hfjet.fDeltaPhi);
  else hfjet.fDeltaR=-99;
  hfjet.fN=constituents.size();

  if (fDoJetSubstructure) SetJetSubstructureVariables(hfjet,constituents);
  
}


//________________________________________________________________
//Set the jet parameters in the AliHFJet object (MC)
void AliHFJetFinder::SetMCJetVariables(AliHFJet& hfjet, const std::vector<fastjet::PseudoJet>& constituents, const fastjet::PseudoJet& jet, Int_t jetID, AliAODMCParticle *mcpart) {

  hfjet.fID=jetID;
  if (mcpart)hfjet.fHFMeson=1;
  else hfjet.fHFMeson=0;
  hfjet.fPt=jet.perp();
  hfjet.fEta=jet.pseudorapidity();
  hfjet.fPhi=jet.phi();
  if (mcpart) hfjet.fDeltaEta=jet.pseudorapidity()-mcpart->Eta();
  else hfjet.fDeltaEta=-99;
  if (mcpart) hfjet.fDeltaPhi=RelativePhi(jet.phi(),mcpart->Phi());
  else hfjet.fDeltaPhi=-99;
  if (mcpart) hfjet.fDeltaR=TMath::Sqrt(hfjet.fDeltaEta*hfjet.fDeltaEta + hfjet.fDeltaPhi*hfjet.fDeltaPhi);
  else hfjet.fDeltaR=-99;
  hfjet.fN=constituents.size();

  if (fDoJetSubstructure) SetJetSubstructureVariables(hfjet,constituents);
  
}




//________________________________________________________________
//Set the jet substrucutre variables in the AliHFJet object (Data and MC)
void AliHFJetFinder::SetJetSubstructureVariables(AliHFJet& hfjet, const std::vector<fastjet::PseudoJet>& constituents) {

  Bool_t softdropped=kFALSE;
  Float_t Pt_jet=0.0;
  Float_t z=0.0;
  Float_t r=0.0;
  Float_t Nsd=0.0;
  Float_t Pt_mother=0.0;
  Float_t k0=0.0, k0_temp=0.0, Zk0=0.0, Rk0=0.0;
  Float_t k1=0.0, k1_temp=0.0, Zk1=0.0, Rk1=0.0;
  Float_t k2=0.0, k2_temp=0.0, Zk2=0.0, Rk2=0.0;
  Float_t kT=0.0, kT_temp=0.0, ZkT=0.0, RkT=0.0;

  if (fSubJetRadius==0.0) fSubJetRadius=fJetRadius*2.5;

  
  fastjet::JetDefinition subJet_definition(JetAlgorithm(fSubJetAlgorithm), fSubJetRadius,RecombinationScheme(fSubJetRecombScheme), fastjet::Best);

  
  try{
    fastjet::ClusterSequence cluster_sequence(constituents, subJet_definition);
    std::vector<fastjet::PseudoJet> reclustered_jet =  cluster_sequence.inclusive_jets(fMinSubJetPt);
    reclustered_jet = sorted_by_pt(reclustered_jet);
         
    fastjet::PseudoJet daughter_jet = reclustered_jet[0];
    fastjet::PseudoJet parent_subjet_1; 
    fastjet::PseudoJet parent_subjet_2;
    Pt_jet=daughter_jet.perp();

    while(daughter_jet.has_parents(parent_subjet_1,parent_subjet_2)){
      if(parent_subjet_1.perp() < parent_subjet_2.perp()) std::swap(parent_subjet_1,parent_subjet_2);
      z=parent_subjet_2.perp()/(parent_subjet_1.perp()+parent_subjet_2.perp());
      r=parent_subjet_1.delta_R(parent_subjet_2);
      Pt_mother=daughter_jet.perp();

      if (z >= fSoftDropZCut*TMath::Power(r/fJetRadius,fSoftDropBeta)){
	if(!softdropped){
	  hfjet.fZg = z;
	  hfjet.fRg = r;
	  hfjet.fPt_mother = Pt_mother;
	  softdropped=kTRUE;
	}
	Nsd++;
      }
      if (Pt_jet > 0.0){
	k0_temp=(1.0/Pt_jet)*z*(1-z)*Pt_mother*TMath::Power(r/fJetRadius,0.1);
	k1_temp=(1.0/Pt_jet)*z*(1-z)*Pt_mother*TMath::Power(r/fJetRadius,1);
	k2_temp=(1.0/Pt_jet)*z*(1-z)*Pt_mother*TMath::Power(r/fJetRadius,2);
	kT_temp=z*Pt_mother*TMath::Sin(r);
	if (k0_temp > k0){
	  k0=k0_temp;
	  Zk0 = z;
	  Rk0 = r;
	}
	if (k1_temp > k1){
	  k1=k1_temp;
	  Zk1 = z;
	  Rk1 = r;
	}
	if (k2_temp > k2){
	  k2=k2_temp;
	  Zk2 = z;
	  Rk2 = r;
	}
	if (kT_temp > kT){
	  kT=kT_temp;
	  ZkT = z;
	  RkT = r;
	}
      }
      
      daughter_jet=parent_subjet_1;
    }
    if (constituents.size() > 1){
      hfjet.fNsd = Nsd;
      hfjet.fk0=k0;
      hfjet.fZk0=Zk0;
      hfjet.fRk0=Rk0;
      hfjet.fk1=k1;
      hfjet.fZk1=Zk1;
      hfjet.fRk1=Rk1;
      hfjet.fk2=k2;
      hfjet.fZk2=Zk2;
      hfjet.fRk2=Rk2;
      hfjet.fkT=kT;
      hfjet.fZkT=ZkT;
      hfjet.fRkT=RkT;
    }

         
  } catch (fastjet::Error) { /*return -1;*/ }


  
}


//________________________________________________________________
//Apply quality check on tracks
Bool_t AliHFJetFinder::CheckTrack(AliAODTrack *track) { 
  if(!track) return false;

  TRandom3 Random;
  Random.SetSeed(0);
  Double_t Random_Number;
  Random_Number=Random.Rndm();
  if(Random_Number > fTrackingEfficiency) return false;

  AliTLorentzVector track_lvec(0,0,0,0);
  track_lvec.SetPtEtaPhiM(track->Pt(), track->Eta(), track->Phi(), 0.139); //set to mass of Pion
 
  if(track_lvec.Pt() > fMaxTrackPt) return false;
  if(track_lvec.Pt() < fMinTrackPt) return false;
  if(track_lvec.E() > fMaxTrackE) return false;
  if(track_lvec.E() < fMinTrackE) return false; 
  if(TMath::Abs(track_lvec.Eta()) > fMaxTrackEta) return false;
  if(TMath::Abs(track_lvec.Phi_0_2pi()) > fMaxTrackPhi) return false;

  if(track->Charge()==0) return false;
  if (!CheckFilterBits(track)) return false;
  
  return true;
}

//________________________________________________________________
//Apply filter bit selection on tracks
Bool_t AliHFJetFinder::CheckFilterBits(AliAODTrack *track) {

  if(!track->TestBit(BIT(4)) && !track->TestBit(BIT(9))) return false;
  if(!track->IsHybridGlobalConstrainedGlobal()) return false;
  return true;
}

//________________________________________________________________
//Apply quality check on particles (MC)
Bool_t AliHFJetFinder::CheckParticle(AliAODMCParticle *particle) {
  if(!particle) return false;
  if(!particle->IsPhysicalPrimary()) return false;
  if(particle->Pt() > fMaxParticlePt) return false;
  if(particle->Pt() < fMinParticlePt) return false;
  if(TMath::Abs(particle->Eta()) > fMaxParticleEta) return false;
  if(fCharged==Charge::charged && particle->Charge()==0) return false;
  if (fCharged==Charge::neutral && particle->Charge()!=0) return false;
  return true;
}



//________________________________________________________________
//Find the index of the jet containing the heavy flavour candidtae or particle (Data and MC)
Int_t AliHFJetFinder::Find_Candidate_Jet() {
  //Finds the label of the jet with the heavy flvaour meson
  std::vector<fastjet::PseudoJet> inclusive_jets = fFastJetWrapper->GetInclusiveJets(); 
  for (UInt_t i=0; i < inclusive_jets.size(); i++){
    std::vector<fastjet::PseudoJet> constituents(fFastJetWrapper->GetJetConstituents(i));
    for (UInt_t j = 0; j < constituents.size(); j++) { 
      if (constituents[j].user_index() == 0) {
	return i; 
      }
    }
  }
  return -1;
}


//________________________________________________________________________
//Transform the range of delta_phi from 0 to 2pi into -pi to pi
Float_t AliHFJetFinder::RelativePhi(Float_t phi_1, Float_t phi_2){

  if(phi_1 < -1*TMath::Pi()) phi_1 += (2*TMath::Pi());                                                             
  else if (phi_1 > TMath::Pi()) phi_1 -= (2*TMath::Pi());
  if(phi_2 < -1*TMath::Pi()) phi_2 += (2*TMath::Pi());
  else if (phi_2 > TMath::Pi()) phi_2 -= (2*TMath::Pi());
  Double_t deltaphi=phi_2-phi_1;
  if(deltaphi < -1*TMath::Pi()) deltaphi += (2*TMath::Pi());
  else if (deltaphi > TMath::Pi()) deltaphi -= (2*TMath::Pi());
  return deltaphi;
}
//________________________________________________________________________
//Returning fastjet jet or subjet algorithm
fastjet::JetFinder AliHFJetFinder::JetAlgorithm(Int_t jetalgo){
  if (jetalgo==JetAlgorithm::kt) return fastjet::kt_algorithm;
  else if (jetalgo==JetAlgorithm::ca) return fastjet::cambridge_algorithm; 
  else return fastjet::antikt_algorithm;
}
//________________________________________________________________________
//Returning fastjet jet or subjet recombination scheme
fastjet::RecombinationScheme AliHFJetFinder::RecombinationScheme(Int_t recombscheme){
  if (recombscheme==RecombScheme::pt_scheme) return fastjet::pt_scheme;
  else return fastjet::E_scheme;
}
//________________________________________________________________________
//Returning fastjet jet area type
fastjet::AreaType AliHFJetFinder::AreaType(Int_t area){
  if (area==AreaType::passive) return fastjet::passive_area;
  if (area==AreaType::voronoi) return fastjet::voronoi_area;
  else return fastjet::active_area;
}

