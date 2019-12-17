/*************************************************************************
* Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
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

//C++
#include <sstream>
#include <array>

// Root
#include <TClonesArray.h>
#include <TDatabasePDG.h>
#include <TParticlePDG.h>
#include <TVector3.h>
#include <THnSparse.h>
#include <TParticle.h>
#include <TMath.h>
#include <THashList.h>
#include <TFile.h>
#include <TRandom3.h>
#include <TGrid.h>
#include <TObjArray.h>

// Aliroot general
#include "AliLog.h"
#include "AliEMCALGeometry.h"
#include "AliAnalysisManager.h"
#include "AliVEventHandler.h"
#include "AliAnalysisDataSlot.h"
#include "AliAnalysisDataContainer.h"

// Aliroot HF
#include "AliAODRecoCascadeHF.h"
#include "AliAODRecoDecayHF2Prong.h"
#include "AliRDHFCutsD0toKpi.h"
#include "AliRDHFCutsDStartoKpipi.h"
#include "AliHFAODMCParticleContainer.h"
#include "AliHFTrackContainer.h"
#include "AliAnalysisVertexingHF.h"

// Aliroot EMCal jet framework
#include "AliEmcalJetTask.h"
#include "AliEmcalJet.h"
#include "AliJetContainer.h"
#include "AliParticleContainer.h"
#include "AliEmcalParticle.h"
#include "AliFJWrapper.h"
#include "AliRhoParameter.h"

#include "AliAnalysisTaskHFSubstructure.h"


#include "FJ_includes.h"





using std::cout;
using std::endl;


ClassImp(AliAnalysisTaskHFSubstructure)

//________________________________________________________________________
AliAnalysisTaskHFSubstructure::AliAnalysisTaskHFSubstructure() : 
AliAnalysisTaskEmcal("AliAnalysisTaskHFSubstructure", kTRUE),
  fJetShapeType(kData),
  fECandidateType(kD0toKpi),
  fIncludeInclusive(kFALSE),
  fIsBDecay(kFALSE),
  fRejectISR(kFALSE),
  fPromptReject(kFALSE),
  fAlienConnect(kFALSE),
  fBranchName(0),
  fCutsType(0),
  fCandidatePDG(421),
  fRejectedOrigin(0),
  fAcceptedDecay(kDecayD0toKpi),
  fJetRadius(0.4),
  fJetMinPt(0.0),
  fTrackingEfficiency(1.0),
  fCandidateArray(0),
  fAodEvent(0),
  fRDHFCuts(0),
  fFastJetWrapper(0),
  fFastJetWrapper_Truth(0),
  fShapesVar_Splittings_DeltaR(0),
  fShapesVar_Splittings_DeltaR_Truth(0),
  fShapesVar_Splittings_Zg(0),
  fShapesVar_Splittings_Zg_Truth(0),
  fShapesVar_Splittings_LeadingSubJetpT(0),
  fShapesVar_Splittings_LeadingSubJetpT_Truth(0),
  fShapesVar_Splittings_HardestSubJetD0(0),
  fShapesVar_Splittings_HardestSubJetD0_Truth(0),
  fShapesVar_Splittings_RadiatorE(0),
  fShapesVar_Splittings_RadiatorE_Truth(0),
  fShapesVar_Splittings_RadiatorpT(0),
  fShapesVar_Splittings_RadiatorpT_Truth(0),

  fTreeResponseMatrixAxis(0),
  fTreeSplittings(0),
  fhEvent(0x0)

{
  for(Int_t i=0;i<nVar;i++){
    fShapesVar[i]=0;
  }
  SetMakeGeneralHistograms(kTRUE);
  DefineOutput(1, TList::Class());
  DefineOutput(2, TTree::Class());
  DefineOutput(3, TTree::Class());
}

//________________________________________________________________________
AliAnalysisTaskHFSubstructure::AliAnalysisTaskHFSubstructure(const char *name) : 
  AliAnalysisTaskEmcal(name, kTRUE),
  fJetShapeType(kData),
  fECandidateType(kD0toKpi),
  fIncludeInclusive(kFALSE),
  fIsBDecay(kFALSE),
  fRejectISR(kFALSE),
  fPromptReject(kFALSE),
  fAlienConnect(kFALSE),
  fBranchName(0),
  fCutsType(0),
  fCandidatePDG(421),
  fRejectedOrigin(0),
  fAcceptedDecay(kDecayD0toKpi),
  fJetRadius(0.4),
  fJetMinPt(0.0),
  fTrackingEfficiency(1.0),
  fCandidateArray(0),
  fAodEvent(0),
  fRDHFCuts(0),
  fFastJetWrapper(0),
  fFastJetWrapper_Truth(0),
  fShapesVar_Splittings_DeltaR(0),
  fShapesVar_Splittings_DeltaR_Truth(0),
  fShapesVar_Splittings_Zg(0),
  fShapesVar_Splittings_Zg_Truth(0),
  fShapesVar_Splittings_LeadingSubJetpT(0),
  fShapesVar_Splittings_LeadingSubJetpT_Truth(0),
  fShapesVar_Splittings_HardestSubJetD0(0),
  fShapesVar_Splittings_HardestSubJetD0_Truth(0),
  fShapesVar_Splittings_RadiatorE(0),
  fShapesVar_Splittings_RadiatorE_Truth(0),
  fShapesVar_Splittings_RadiatorpT(0),
  fShapesVar_Splittings_RadiatorpT_Truth(0),

  fTreeResponseMatrixAxis(0),
  fTreeSplittings(0),
  fhEvent(0x0)
{
  // Standard constructor.
  for(Int_t i=0;i<nVar;i++){
    fShapesVar[i]=0;
  }
  SetMakeGeneralHistograms(kTRUE);
  DefineOutput(1, TList::Class());
  DefineOutput(2, TTree::Class());
  DefineOutput(3, TTree::Class());
}

//________________________________________________________________________
AliAnalysisTaskHFSubstructure::~AliAnalysisTaskHFSubstructure()
{
  // Destructor.
}

//________________________________________________________________________
void AliAnalysisTaskHFSubstructure::UserCreateOutputObjects()
{
  // Create user output.

  AliAnalysisTaskEmcal::UserCreateOutputObjects();

  Bool_t oldStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);
  TH1::AddDirectory(oldStatus);
  //create a tree used for the MC data and making a 4D response matrix

  
  
  const char* nameoutput = GetOutputSlot(2)->GetContainer()->GetName();
  fTreeResponseMatrixAxis = new TTree(nameoutput, nameoutput);
  TString *fShapesVarNames = new TString [nVar];
  fShapesVarNames[0] = "pT_Jet";
  fShapesVarNames[1] = "pT_Jet_Truth";
  fShapesVarNames[2] = "pT_D";
  fShapesVarNames[3] = "pT_D_Truth";
  fShapesVarNames[4] = "Inv_M_D";
  fShapesVarNames[5] = "Inv_M_D_Truth";
  fShapesVarNames[6] = "Flag_D";
  fShapesVarNames[7] = "Flag_D_Truth";
  fShapesVarNames[8] = "Prompt_PDG";
  fShapesVarNames[9] = "Prompt_PDG_Truth";
  // fShapesVarNames[10] = "NTracks";
  //fShapesVarNames[11] = "NTracks_Truth";
  //  fShapesVarNames[12] = "Eta_Jet";
  //  fShapesVarNames[13] = "Eta_Jet_Truth";
  //  fShapesVarNames[14] = "Eta_D";
  //  fShapesVarNames[15] = "Eta_D_Truth";
  // fShapesVarNames[16] = "Y_D";
  // fShapesVarNames[17] = "Y_D_Truth";
  
  for(Int_t ivar=0; ivar < nVar; ivar++){
    cout<<"looping over variables"<<endl;
    fTreeResponseMatrixAxis->Branch(fShapesVarNames[ivar].Data(), &fShapesVar[ivar], Form("%s/D", fShapesVarNames[ivar].Data()));
  }

  const char* nameoutput_Splittings = GetOutputSlot(3)->GetContainer()->GetName();
  fTreeSplittings = new TTree(nameoutput_Splittings, nameoutput_Splittings);
  TString *fShapesVarNames_Splittings=new TString[nVar_Splittings];
  fShapesVarNames_Splittings[0] = "DeltaR";
  fShapesVarNames_Splittings[1] = "DeltaR_Truth";
  fShapesVarNames_Splittings[2] = "Zg";
  fShapesVarNames_Splittings[3] = "Zg_Truth";
  // fShapesVarNames_Splittings[4] = "LeadingSubJetpT";
  //fShapesVarNames_Splittings[5] = "LeadingSubJetpT_Truth";
  fShapesVarNames_Splittings[4] = "HardestSubJetD0";
  fShapesVarNames_Splittings[5] = "HardestSubJetD0_Truth";
  fShapesVarNames_Splittings[6] = "RadiatorE";
  fShapesVarNames_Splittings[7] = "RadiatorE_Truth";
  fShapesVarNames_Splittings[8] = "RadiatorpT";
  fShapesVarNames_Splittings[9] = "RadiatorpT_Truth";
  fTreeSplittings->Branch(fShapesVarNames_Splittings[0].Data(), &fShapesVar_Splittings_DeltaR, 0,1);
  fTreeSplittings->Branch(fShapesVarNames_Splittings[1].Data(), &fShapesVar_Splittings_DeltaR_Truth, 0,1);
  fTreeSplittings->Branch(fShapesVarNames_Splittings[2].Data(), &fShapesVar_Splittings_Zg, 0,1);
  fTreeSplittings->Branch(fShapesVarNames_Splittings[3].Data(), &fShapesVar_Splittings_Zg_Truth, 0,1);
  // fTreeSplittings->Branch(fShapesVarNames_Splittings[4].Data(), &fShapesVar_Splittings_LeadingSubJetpT, 0,1);
  // fTreeSplittings->Branch(fShapesVarNames_Splittings[5].Data(), &fShapesVar_Splittings_LeadingSubJetpT_Truth, 0,1);
  fTreeSplittings->Branch(fShapesVarNames_Splittings[4].Data(), &fShapesVar_Splittings_HardestSubJetD0, 0,1);
  fTreeSplittings->Branch(fShapesVarNames_Splittings[5].Data(), &fShapesVar_Splittings_HardestSubJetD0_Truth, 0,1);
  fTreeSplittings->Branch(fShapesVarNames_Splittings[6].Data(), &fShapesVar_Splittings_RadiatorE, 0,1);
  fTreeSplittings->Branch(fShapesVarNames_Splittings[7].Data(), &fShapesVar_Splittings_RadiatorE_Truth, 0,1);
  fTreeSplittings->Branch(fShapesVarNames_Splittings[8].Data(), &fShapesVar_Splittings_RadiatorpT, 0,1);
  fTreeSplittings->Branch(fShapesVarNames_Splittings[9].Data(), &fShapesVar_Splittings_RadiatorpT_Truth, 0,1);

  fhEvent=new TH1D("fhEvent","fhEvent",40,-0.5,39.5);
  fOutput->Add(fhEvent);
  
  PostData(1,fOutput);
  PostData(2,fTreeResponseMatrixAxis);
  PostData(3,fTreeSplittings);
  // delete [] fShapesVarNames;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskHFSubstructure::Run()
{
  // Run analysis code here, if needed. It will be executed before FillHistograms().


  return kTRUE;
}
//________________________________________________________________________
Bool_t AliAnalysisTaskHFSubstructure::FillHistograms()
{

  Int_t Matching_AOD_deltaAODlevel = AliRDHFCuts::CheckMatchingAODdeltaAODevents();
  if (Matching_AOD_deltaAODlevel <= 0) return kTRUE;

  fhEvent->Fill(0); 
  //if(fCutsFileName.Contains("alien://") && fAlienConnect) TGrid::Connect("alien://");
  //TFile* Cuts_File = TFile::Open(fCutsFileName); 
  //TString cutsname="D0toKpiCuts";
  //if (fCutsType!="") cutsname += TString::Format("_%s", fCutsType.Data()); 
  //fRDHFCuts = dynamic_cast<AliRDHFCuts*>(fCutsFile->Get(cutsname));

  
  fAodEvent = dynamic_cast<AliAODEvent*>(fInputEvent);


  if(!fRDHFCuts->IsEventSelected(fAodEvent)) return kTRUE;
  fhEvent->Fill(1); 
  fFastJetWrapper=new AliFJWrapper("fastjetwrapper","fastjetwrapper");
  fFastJetWrapper->SetAreaType(fastjet::active_area); 
  fFastJetWrapper->SetGhostArea(0.005);  
  fFastJetWrapper->SetR(fJetRadius);
  fFastJetWrapper->SetAlgorithm(fastjet::antikt_algorithm);
  fFastJetWrapper->SetRecombScheme(static_cast<fastjet::RecombinationScheme>(0));

  fFastJetWrapper_Truth=new AliFJWrapper("fastjetwrapper_truth","fastjetwrapper_truth");

  if (fJetShapeType == kData || fJetShapeType == kDetSignal || fJetShapeType == kDetBackground || fJetShapeType == kDetReflection || fJetShapeType == kDet){

    TRandom3 Random;
    Random.SetSeed(0);
    Double_t Random_Number;
    
    fCandidateArray = dynamic_cast<TClonesArray*>(fAodEvent->GetList()->FindObject("D0toKpi"));

    AliHFAODMCParticleContainer *Particle_Container = NULL;
    if (fJetShapeType != kData) Particle_Container = (AliHFAODMCParticleContainer *) GetParticleContainer(1); 
    // if (!Particle_Container) continue;
       

    std::vector<AliAODRecoDecayHF2Prong*> D_Candidates_Vector;
    D_Candidates_Vector.clear();


    Int_t N_DMesons=0;
    //AliHFTrackContainer* Track_Container=(AliHFTrackContainer *) GetTrackContainer(0);
    AliHFTrackContainer* Track_Container=dynamic_cast<AliHFTrackContainer*>(GetTrackContainer(0));
    if (!Track_Container) return kTRUE;
    Track_Container->SetDMesonCandidate(NULL);
    for (Int_t i_D = 0; i_D < fCandidateArray->GetEntriesFast(); i_D++) {
      AliAODRecoDecayHF2Prong* D_Candidate = static_cast<AliAODRecoDecayHF2Prong*>(fCandidateArray->At(i_D));
      if (!D_Candidate) continue;
      if (!fRDHFCuts->IsInFiducialAcceptance(D_Candidate->Pt(), D_Candidate->Y(fCandidatePDG))) continue;
      

      
      Int_t Mass_Hypo_Type=fRDHFCuts->IsSelected(D_Candidate, AliRDHFCuts::kAll, fAodEvent);
      Int_t N_Mass_Hypotheses=1;
      if (Mass_Hypo_Type <= 0 || Mass_Hypo_Type>3) continue;
      else if (Mass_Hypo_Type ==3) N_Mass_Hypotheses=2;
      fhEvent->Fill(2);
      Int_t Matched_Truth_Particle_PDG=0;
      Int_t Is_Prompt_Correct_Quark_PDG=-1;
      if (fJetShapeType != kData){    

   
	const Int_t D_Candidtae_N_Daughters=2;
	Int_t D_Candidtae_Daughters_PDG[D_Candidtae_N_Daughters] = {211,321}; 
	Int_t D_Candidate_MatchedTruth_Label = D_Candidate->MatchToMC(fCandidatePDG, Particle_Container->GetArray(), D_Candidtae_N_Daughters, D_Candidtae_Daughters_PDG); 
	Bool_t Is_Prompt_Correct_Quark=kFALSE;
	if (D_Candidate_MatchedTruth_Label >= 0) {
	  AliAODMCParticle* Matched_Truth_Particle = static_cast<AliAODMCParticle*>(Particle_Container->GetArray()->At(D_Candidate_MatchedTruth_Label)); 
	  if (Matched_Truth_Particle) {
	    fhEvent->Fill(3); 
	   
	    Int_t Matched_Truth_Particle_Mother_Label=Matched_Truth_Particle->GetMother(); 
	    while (Matched_Truth_Particle_Mother_Label >= 0) {
	      AliAODMCParticle* Matched_Truth_Particle_Mother = static_cast<AliAODMCParticle*>(Particle_Container->GetArray()->At(Matched_Truth_Particle_Mother_Label)); 
	      if (Matched_Truth_Particle_Mother){
		Int_t Original_Quark_PDG=4; 
		if (fIsBDecay) Original_Quark_PDG=5; 
		if (TMath::Abs(Matched_Truth_Particle_Mother->GetPdgCode())==Original_Quark_PDG) Is_Prompt_Correct_Quark=kTRUE; 
		if (TMath::Abs(Matched_Truth_Particle_Mother->GetPdgCode()) == 4){
		  Is_Prompt_Correct_Quark_PDG=4;
		  fhEvent->Fill(4);
		  break;
		}
		if (TMath::Abs(Matched_Truth_Particle_Mother->GetPdgCode()) == 5){
		  Is_Prompt_Correct_Quark_PDG=5;
		  fhEvent->Fill(5);
		  break;
		}
		if (Matched_Truth_Particle_Mother_Label==Matched_Truth_Particle_Mother->GetMother()) break; 
	        Matched_Truth_Particle_Mother_Label=Matched_Truth_Particle_Mother->GetMother();
	      }
	      else break;
	      //delete Matched_Truth_Particle_Mother;
	    }
	    Matched_Truth_Particle_PDG = Matched_Truth_Particle->PdgCode();
	  }
	  //delete Matched_Truth_Particle;
	}
	//else continue; 
	if (fPromptReject && !Is_Prompt_Correct_Quark) continue; 

	//if (TMath::Abs(Matched_Truth_Particle_PDG)!=fCandidatePDG) continue; 

      }


      Double_t Inv_Mass_D=0.0;
      if (Mass_Hypo_Type==1){ 
	if (fJetShapeType==kData || fJetShapeType == kDet || (fJetShapeType==kDetSignal && Matched_Truth_Particle_PDG==fCandidatePDG) || (fJetShapeType==kDetBackground && Matched_Truth_Particle_PDG!=fCandidatePDG) || (fJetShapeType==kDetReflection && Matched_Truth_Particle_PDG==-fCandidatePDG)){
	  Inv_Mass_D=D_Candidate->InvMassD0();
	  fhEvent->Fill(6); 
	}
	else{
	  fhEvent->Fill(7); 
	  continue;
	}
      }

      if (Mass_Hypo_Type==2){ 
	if (fJetShapeType==kData || fJetShapeType == kDet || (fJetShapeType==kDetSignal && Matched_Truth_Particle_PDG==-fCandidatePDG) || (fJetShapeType==kDetBackground && Matched_Truth_Particle_PDG!=-fCandidatePDG) || (fJetShapeType==kDetReflection && Matched_Truth_Particle_PDG==fCandidatePDG)){
	  Inv_Mass_D=D_Candidate->InvMassD0bar();
	  fhEvent->Fill(8); 
	}
	else{
	  fhEvent->Fill(9); 
	  continue;
	}
      }

      
      
      for (Int_t i_Mass_Hypotheses=0; i_Mass_Hypotheses<N_Mass_Hypotheses; i_Mass_Hypotheses++){


	if (Mass_Hypo_Type==3){ 
	    if(i_Mass_Hypotheses==0){
	      if (fJetShapeType==kData || fJetShapeType == kDet || (fJetShapeType==kDetSignal && Matched_Truth_Particle_PDG==fCandidatePDG) || (fJetShapeType==kDetBackground && Matched_Truth_Particle_PDG!=fCandidatePDG) || (fJetShapeType==kDetReflection && Matched_Truth_Particle_PDG==-fCandidatePDG)){
		Inv_Mass_D=D_Candidate->InvMassD0();
		fhEvent->Fill(11); 
	      }
	      else{
		fhEvent->Fill(12); 
		continue;
	      }
	    }
	    if(i_Mass_Hypotheses==1){
	      if (fJetShapeType==kData || fJetShapeType == kDet || (fJetShapeType==kDetSignal && Matched_Truth_Particle_PDG==-fCandidatePDG) || (fJetShapeType==kDetBackground && Matched_Truth_Particle_PDG!=-fCandidatePDG) || (fJetShapeType==kDetReflection && Matched_Truth_Particle_PDG==fCandidatePDG)){
		Inv_Mass_D=D_Candidate->InvMassD0bar();
		fhEvent->Fill(13); 
	      }
	      else{
		fhEvent->Fill(14);
		continue;
	      }
	    }
	}
	//Random_Number=Random.Rndm();
	//if(Random_Number > fTrackingEfficiency*fTrackingEfficiency) continue; // here it shows that the D did not get reconstructed cause one of the daughters was missing...however should we do this before incase the same daughter is involved multiple times?
	fFastJetWrapper->Clear();
	AliTLorentzVector D_Candidate_LorentzVector(0,0,0,0);
	D_Candidate_LorentzVector.SetPtEtaPhiM(D_Candidate->Pt(), D_Candidate->Eta(), D_Candidate->Phi(), Inv_Mass_D);
	fFastJetWrapper->AddInputVector(D_Candidate_LorentzVector.Px(), D_Candidate_LorentzVector.Py(), D_Candidate_LorentzVector.Pz(), D_Candidate_LorentzVector.E(), 0); 
	Track_Container->SetDMesonCandidate(D_Candidate);
	AliAODTrack *Track = NULL;
	for (Int_t i_Track=0; i_Track<Track_Container->GetNTracks(); i_Track++){
	  Track = static_cast<AliAODTrack*>(Track_Container->GetAcceptParticle(i_Track));
	  if(!Track) continue;
	  if (Track->Pt() > 100.0 || TMath::Abs(Track->Eta()) > 0.9) continue;
	  Random_Number=Random.Rndm();
	  if(Random_Number > fTrackingEfficiency) continue;
	  fFastJetWrapper->AddInputVector(Track->Px(), Track->Py(), Track->Pz(), Track->E(),i_Track+100);
	}
	//	delete Track;

	fFastJetWrapper->Run();
	std::vector<fastjet::PseudoJet> Inclusive_Jets = fFastJetWrapper->GetInclusiveJets(); 
	for (UInt_t i_Jet=0; i_Jet < Inclusive_Jets.size(); i_Jet++){
	  Bool_t Is_D_Jet=kFALSE;
	  if (Inclusive_Jets[i_Jet].perp()<fJetMinPt) continue;
	  if (TMath::Abs(Inclusive_Jets[i_Jet].pseudorapidity()) > 0.9-fJetRadius) continue;
	  std::vector<fastjet::PseudoJet> Constituents(fFastJetWrapper->GetJetConstituents(i_Jet));
	  for (UInt_t i_Constituents = 0; i_Constituents < Constituents.size(); i_Constituents++) { 
	    if (Constituents[i_Constituents].user_index() == 0) Is_D_Jet = kTRUE; 	    
	  }
	  if (!Is_D_Jet) continue; 
	  fhEvent->Fill(10);
	  std::vector<Double_t> Splittings_Zg;
	  std::vector<Double_t> Splittings_DeltaR;
	  std::vector<Double_t> Splittings_LeadingSubJetpT;
	  std::vector<Double_t> Splittings_HardestSubJetD0;
	  std::vector<Double_t> Splittings_RadiatorE;
	  std::vector<Double_t> Splittings_RadiatorpT;

	  Bool_t Is_D_SubJet=kFALSE;
	fastjet::JetDefinition Jet_Definition(fastjet::cambridge_algorithm, fJetRadius*2.5,static_cast<fastjet::RecombinationScheme>(0), fastjet::Best); 
	try{
	  std::vector<fastjet::PseudoJet> Reclustered_Particles(fFastJetWrapper->GetJetConstituents(i_Jet));
	  fastjet::ClusterSequence Cluster_Sequence_CA(Reclustered_Particles, Jet_Definition);
	  std::vector<fastjet::PseudoJet> Reclustered_Jet =  Cluster_Sequence_CA.inclusive_jets(0.0);
	  Reclustered_Jet = sorted_by_pt(Reclustered_Jet); 
         

	  fastjet::PseudoJet Daughter_Jet = Reclustered_Jet[0];
	  fastjet::PseudoJet Parent_SubJet_1; 
	  fastjet::PseudoJet Parent_SubJet_2;  
    
	  while(Daughter_Jet.has_parents(Parent_SubJet_1,Parent_SubJet_2)){
	    if(Parent_SubJet_1.perp() < Parent_SubJet_2.perp()) std::swap(Parent_SubJet_1,Parent_SubJet_2);
	    // Splittings_LeadingSubJetpT.push_back(Parent_SubJet_1.perp());
	    vector < fastjet::PseudoJet > Hard_SubJet_Constituents = sorted_by_pt(Parent_SubJet_1.constituents());
	    Is_D_SubJet=kFALSE;
	    for(UInt_t j=0;j<Hard_SubJet_Constituents.size();j++){
	      if(Hard_SubJet_Constituents[j].user_index()==0) Is_D_SubJet=kTRUE;	      
	    }

	    if (!Is_D_SubJet) Splittings_HardestSubJetD0.push_back(1.0); 
	    else Splittings_HardestSubJetD0.push_back(2.0);  

	    // if(!Is_D_SubJet) std::swap(Parent_SubJet_1,Parent_SubJet_2);     
	    Splittings_DeltaR.push_back(Parent_SubJet_1.delta_R(Parent_SubJet_2));
	    Splittings_Zg.push_back(Parent_SubJet_2.perp()/(Parent_SubJet_1.perp()+Parent_SubJet_2.perp()));
	    Splittings_RadiatorE.push_back(Daughter_Jet.E());
	    Splittings_RadiatorpT.push_back(Daughter_Jet.perp());
	    Daughter_Jet=Parent_SubJet_1;
	  }

         
	} catch (fastjet::Error) { /*return -1;*/ }


	
	//	for (Int_t i_Mass_Hypotheses=0; i_Mass_Hypotheses<N_Mass_Hypotheses; i_Mass_Hypotheses++){ 

	      
	  D_Candidates_Vector.push_back(D_Candidate); 
	  if (i_Mass_Hypotheses==0) fhEvent->Fill(15); 
	  N_DMesons++;

	  Double_t Flag_D=-1.0;
	  if (Mass_Hypo_Type ==1) Flag_D=1.0; 
	  else if (Mass_Hypo_Type ==2) Flag_D=2.0; 
	  else if (Mass_Hypo_Type ==3 && i_Mass_Hypotheses==0) Flag_D=3.0; 
	  else if (Mass_Hypo_Type ==3 && i_Mass_Hypotheses==1) Flag_D=4.0; 

	
	  fShapesVar[0] = Inclusive_Jets[i_Jet].perp();
	  fShapesVar[1] = 0.0;
	  fShapesVar[2] = D_Candidate->Pt();
	  fShapesVar[3] = 0.0;
	  fShapesVar[4] = Inv_Mass_D;
	  fShapesVar[5] = 0.0;
	  fShapesVar[6] = Flag_D;
	  if (fJetShapeType == kData) fShapesVar[7] = 0.0;
	  else fShapesVar[7] = Matched_Truth_Particle_PDG;
	  fShapesVar[8] = Is_Prompt_Correct_Quark_PDG;
	  fShapesVar[9] = 0.0;
	  //  fShapesVar[10] = NTracks;
	  //fShapesVar[11] = 0.0;
	  // fShapesVar[12] = TMath::Abs(Inclusive_Jets[i_Jet].pseudorapidity());
	  // fShapesVar[13] = 0.0;
	  // fShapesVar[14] = Dmeson_Eta;
	  //fShapesVar[15] = 0.0;
	  //fShapesVar[16] = Dmeson_Y;
	  //fShapesVar[17] = 0.0;
	  


	  fShapesVar_Splittings_DeltaR.push_back(Splittings_DeltaR);
	  fShapesVar_Splittings_DeltaR_Truth.push_back(Splittings_DeltaR); 
	  fShapesVar_Splittings_Zg.push_back(Splittings_Zg);
	  fShapesVar_Splittings_Zg_Truth.push_back(Splittings_Zg); 
	  // fShapesVar_Splittings_LeadingSubJetpT.push_back(Splittings_LeadingSubJetpT);
	  //fShapesVar_Splittings_LeadingSubJetpT_Truth.push_back(Splittings_LeadingSubJetpT); 
	  fShapesVar_Splittings_HardestSubJetD0.push_back(Splittings_HardestSubJetD0);
	  fShapesVar_Splittings_HardestSubJetD0_Truth.push_back(Splittings_HardestSubJetD0);
	  fShapesVar_Splittings_RadiatorE.push_back(Splittings_RadiatorE);
	  fShapesVar_Splittings_RadiatorE_Truth.push_back(Splittings_RadiatorE);
	  fShapesVar_Splittings_RadiatorpT.push_back(Splittings_RadiatorpT);
	  fShapesVar_Splittings_RadiatorpT_Truth.push_back(Splittings_RadiatorpT); 

	  fTreeResponseMatrixAxis->Fill();
	  fTreeSplittings->Fill();
	  fShapesVar_Splittings_DeltaR.clear();
	  fShapesVar_Splittings_DeltaR_Truth.clear();
	  fShapesVar_Splittings_Zg.clear();
	  fShapesVar_Splittings_Zg_Truth.clear();
	  // fShapesVar_Splittings_LeadingSubJetpT.clear();
	  //fShapesVar_Splittings_LeadingSubJetpT_Truth.clear();
	  fShapesVar_Splittings_HardestSubJetD0.clear();
	  fShapesVar_Splittings_HardestSubJetD0_Truth.clear();
	  fShapesVar_Splittings_RadiatorE.clear();
	  fShapesVar_Splittings_RadiatorE_Truth.clear();
	  fShapesVar_Splittings_RadiatorpT.clear();
	  fShapesVar_Splittings_RadiatorpT_Truth.clear(); 
      }
      }
      //delete D_Candidate;
    }
    if(N_DMesons==0) fhEvent->Fill(16); 
    if(N_DMesons==1) fhEvent->Fill(17); 
    if(N_DMesons==2) fhEvent->Fill(18); 
    if(N_DMesons==3) fhEvent->Fill(19); 
    if(N_DMesons==4) fhEvent->Fill(20); 
    if(N_DMesons==5) fhEvent->Fill(21); 
    if(N_DMesons==6) fhEvent->Fill(22); 

    // delete Track_Container;
    // if (fJetShapeType != kData) delete Particle_Container;
  
  }



  if (fJetShapeType == kTrueDet){ 
    
    TRandom3 Random;
    Random.SetSeed(0);
    Double_t Random_Number;
    
    AliHFAODMCParticleContainer *Particle_Container = (AliHFAODMCParticleContainer*) GetParticleContainer(1); 
    Particle_Container->SetSpecialPDG(fCandidatePDG); 
    Particle_Container->SetRejectedOriginMap(fRejectedOrigin); 
    Particle_Container->SetAcceptedDecayMap(fAcceptedDecay);
    Particle_Container->SetRejectISR(fRejectISR);
    Particle_Container->SetCharge(AliParticleContainer::EChargeCut_t::kCharged);

    std::vector<fastjet::PseudoJet> Inclusive_Jets_Truth; 
    std::vector<std::pair<Int_t, Int_t>> Inclusive_Jets_Truth_Labels; 
    std::vector<Int_t> Unmatched_Truth_Level_D; 
    Int_t NMatched_DMeson_Jets=0;

    // Double_t NTracks=0;
    // Double_t NTracks_Truth=0;
    //Double_t Jet_Eta=0;
    //Double_t Jet_Eta_Truth=0;
    //Double_t Dmeson_Eta=0;
    //Double_t Dmeson_Eta_Truth=0;
    //Double_t Dmeson_Y=0;
    //Double_t Dmeson_Y_Truth=0;

    
    
    if (Particle_Container->IsSpecialPDGFound()){ 

      fhEvent->Fill(2); 

      fFastJetWrapper_Truth->SetAreaType(fastjet::active_area); 
      fFastJetWrapper_Truth->SetGhostArea(0.005);
      fFastJetWrapper_Truth->SetR(fJetRadius);
      fFastJetWrapper_Truth->SetAlgorithm(fastjet::antikt_algorithm);
      fFastJetWrapper_Truth->SetRecombScheme(static_cast<fastjet::RecombinationScheme>(0));
      fFastJetWrapper_Truth->Clear();
      
      AliAODMCParticle* Truth_Particle=NULL; 
      Int_t NTruthD=0; 
      for (Int_t i_Particle=0; i_Particle<Particle_Container->GetNParticles(); i_Particle++){ 
	Truth_Particle = static_cast<AliAODMCParticle*>(Particle_Container->GetAcceptMCParticle(i_Particle));
	if (!Truth_Particle) continue;
	//	if (TMath::Abs(Truth_Particle->Eta())>0.9) continue;
	if (TMath::Abs(Truth_Particle->PdgCode())==fCandidatePDG){
	  if (Truth_Particle->Pt() > 5.0){
	    if (TMath::Abs(Truth_Particle->Y()) > 0.8) continue;
	  }
	  else{
	    if(Truth_Particle->Y() < 0.2/15*Truth_Particle->Pt()*Truth_Particle->Pt()-1.9/15*Truth_Particle->Pt()-0.5 || Truth_Particle->Y() > -0.2/15*Truth_Particle->Pt()*Truth_Particle->Pt()+1.9/15*Truth_Particle->Pt()+0.5) continue;
	  }
	  std::pair<Int_t, Int_t> Inclusive_Jet_Truth_Labels;
	  Inclusive_Jet_Truth_Labels.first=Truth_Particle->GetLabel(); 
	  Inclusive_Jet_Truth_Labels.second=NTruthD; 
	  Inclusive_Jets_Truth_Labels.push_back(Inclusive_Jet_Truth_Labels);
	  Unmatched_Truth_Level_D.push_back(NTruthD);
	  fFastJetWrapper_Truth->AddInputVector(Truth_Particle->Px(), Truth_Particle->Py(), Truth_Particle->Pz(), Truth_Particle->E(),NTruthD); 
	  NTruthD++;
	}
	else{
	  if (TMath::Abs(Truth_Particle->Eta()) > 0.9) continue; 
	  fFastJetWrapper_Truth->AddInputVector(Truth_Particle->Px(), Truth_Particle->Py(), Truth_Particle->Pz(), Truth_Particle->E(),i_Particle+100);
	}
      }
      // delete Truth_Particle;
      fFastJetWrapper_Truth->Run();
      Inclusive_Jets_Truth = fFastJetWrapper_Truth->GetInclusiveJets(); 
      if (NTruthD==0) fhEvent->Fill(3); 
      if (NTruthD==1) fhEvent->Fill(4); 
      if (NTruthD==2) fhEvent->Fill(5); 
      if (NTruthD==3) fhEvent->Fill(6); 
      if (NTruthD==4) fhEvent->Fill(7); 
      if (NTruthD==5) fhEvent->Fill(8); 
      if (NTruthD==6) fhEvent->Fill(9); 
    }


  


    fCandidateArray = dynamic_cast<TClonesArray*>(fAodEvent->GetList()->FindObject("D0toKpi"));
    AliHFTrackContainer* Track_Container=(AliHFTrackContainer *) GetTrackContainer(0);
    if (!Track_Container) return kTRUE;
    Track_Container->SetDMesonCandidate(NULL);
    
    for (Int_t i_D = 0; i_D < fCandidateArray->GetEntriesFast(); i_D++) {  
      AliAODRecoDecayHF2Prong* D_Candidate = static_cast<AliAODRecoDecayHF2Prong*>(fCandidateArray->At(i_D)); 
      if (!D_Candidate) continue;
      if (!fRDHFCuts->IsInFiducialAcceptance(D_Candidate->Pt(), D_Candidate->Y(fCandidatePDG))) continue; 

      Int_t Mass_Hypo_Type=fRDHFCuts->IsSelected(D_Candidate, AliRDHFCuts::kAll, fAodEvent); 
      if (Mass_Hypo_Type <= 0 || Mass_Hypo_Type>3) continue; 
      fhEvent->Fill(10); 

      Int_t Matched_Truth_Particle_PDG=0;
      
   
      const Int_t D_Candidtae_N_Daughters=2;
      Int_t D_Candidtae_Daughters_PDG[D_Candidtae_N_Daughters] = {211,321}; 
      Int_t D_Candidate_MatchedTruth_Label = D_Candidate->MatchToMC(fCandidatePDG, Particle_Container->GetArray(), D_Candidtae_N_Daughters, D_Candidtae_Daughters_PDG); 
      Bool_t Is_Prompt_Correct_Quark=kFALSE;
      Int_t Is_Prompt_Correct_Quark_PDG=-1;
      AliAODMCParticle* Matched_Truth_Particle; 
      if (D_Candidate_MatchedTruth_Label >= 0) { 
	Matched_Truth_Particle = static_cast<AliAODMCParticle*>(Particle_Container->GetArray()->At(D_Candidate_MatchedTruth_Label)); 
	if (Matched_Truth_Particle) {
	  fhEvent->Fill(11); 
	  
	  Int_t Matched_Truth_Particle_Mother_Label=Matched_Truth_Particle->GetMother(); 
	  while (Matched_Truth_Particle_Mother_Label >= 0) {
	    AliAODMCParticle* Matched_Truth_Particle_Mother = static_cast<AliAODMCParticle*>(Particle_Container->GetArray()->At(Matched_Truth_Particle_Mother_Label)); 
	    if (Matched_Truth_Particle_Mother){
	      Int_t Original_Quark_PDG=4; 
	      if (fIsBDecay) Original_Quark_PDG=5; 
	      if (TMath::Abs(Matched_Truth_Particle_Mother->GetPdgCode()) == 4){
		fhEvent->Fill(12);		
		Is_Prompt_Correct_Quark_PDG=4;
		break;
	      }
	      if (TMath::Abs(Matched_Truth_Particle_Mother->GetPdgCode()) == 5){
		fhEvent->Fill(13);
		Is_Prompt_Correct_Quark_PDG=5;
		break;
	      }
	      if (TMath::Abs(Matched_Truth_Particle_Mother->GetPdgCode())==Original_Quark_PDG) Is_Prompt_Correct_Quark=kTRUE; 
	      if (Matched_Truth_Particle_Mother_Label==Matched_Truth_Particle_Mother->GetMother()) break; 
	      Matched_Truth_Particle_Mother_Label=Matched_Truth_Particle_Mother->GetMother();
	    }
	    else break;
	    // delete Matched_Truth_Particle_Mother;
	  }
	  Matched_Truth_Particle_PDG = Matched_Truth_Particle->PdgCode();
	}
      }
      else continue;
      if (fPromptReject && !Is_Prompt_Correct_Quark) continue; 

    
      if (Mass_Hypo_Type==1 && Matched_Truth_Particle_PDG!=fCandidatePDG){
	fhEvent->Fill(14); 
	continue;
      }
      else fhEvent->Fill(15); 
      if (Mass_Hypo_Type==2 && Matched_Truth_Particle_PDG!=-fCandidatePDG){
	fhEvent->Fill(16); 
	continue;
      }
      else fhEvent->Fill(17); 
      if (Mass_Hypo_Type==3 && TMath::Abs(Matched_Truth_Particle_PDG)!=fCandidatePDG){
	fhEvent->Fill(18); 
	continue;
      }
      else{
	fhEvent->Fill(19);
	if (Matched_Truth_Particle_PDG==fCandidatePDG) fhEvent->Fill(20); 
	if (Matched_Truth_Particle_PDG==-fCandidatePDG) fhEvent->Fill(21); 
      }


      Double_t Inv_Mass_D=0.0;
      Double_t Inv_Mass_D_Truth=0.0;

      if (Mass_Hypo_Type==1 || (Mass_Hypo_Type==3 && Matched_Truth_Particle_PDG==fCandidatePDG)){
	Inv_Mass_D=D_Candidate->InvMassD0();
	//Inv_Mass_D_Truth=Matched_Truth_Particle->InvMassD0();
	Inv_Mass_D_Truth=0.0;
      }
      if (Mass_Hypo_Type==2 || (Mass_Hypo_Type==3 && Matched_Truth_Particle_PDG==-fCandidatePDG)){
	Inv_Mass_D=D_Candidate->InvMassD0bar();
	//Inv_Mass_D_Truth=Matched_Truth_Particle->InvMassD0bar();
	Inv_Mass_D_Truth=0.0;
      }

      //Random_Number=Random.Rndm();
      //if(Random_Number > fTrackingEfficiency*fTrackingEfficiency) continue;
  
      fFastJetWrapper->Clear();
      AliTLorentzVector D_Candidate_LorentzVector(0,0,0,0);
      D_Candidate_LorentzVector.SetPtEtaPhiM(D_Candidate->Pt(), D_Candidate->Eta(), D_Candidate->Phi(), Inv_Mass_D);
      //if (TMath::Abs(D_Candidate->Eta())>0.9) continue;
      // Dmeson_Eta=TMath::Abs(D_Candidate->Eta());
      //Dmeson_Y=TMath::Abs(D_Candidate->Y(fCandidatePDG));
      fFastJetWrapper->AddInputVector(D_Candidate_LorentzVector.Px(), D_Candidate_LorentzVector.Py(), D_Candidate_LorentzVector.Pz(), D_Candidate_LorentzVector.E(), 0);

    
      if (!Track_Container) continue;
      Track_Container->SetDMesonCandidate(D_Candidate); 
      AliAODTrack *Track = NULL;
      for (Int_t i_Track=0; i_Track<Track_Container->GetNTracks(); i_Track++){ 
	Track = static_cast<AliAODTrack*>(Track_Container->GetAcceptParticle(i_Track));
	if(!Track) continue;
	if (Track->Pt() > 100.0 || TMath::Abs(Track->Eta())>0.9) continue;
	Random_Number=Random.Rndm();
	if(Random_Number > fTrackingEfficiency) continue;
	fFastJetWrapper->AddInputVector(Track->Px(), Track->Py(), Track->Pz(), Track->E(),i_Track+100); 
      }
      // delete Track;
      fFastJetWrapper->Run(); 

      
      std::vector<fastjet::PseudoJet> Inclusive_Jets = fFastJetWrapper->GetInclusiveJets(); 
      for (UInt_t i_Jet=0; i_Jet < Inclusive_Jets.size(); i_Jet++){ 
	Bool_t Is_D_Jet=kFALSE;
	if (Inclusive_Jets[i_Jet].perp()<fJetMinPt) continue;
	//Jet_Eta=TMath::Abs(Inclusive_Jets[i_Jet].pseudorapidity());
	if (TMath::Abs(Inclusive_Jets[i_Jet].pseudorapidity()) > 0.9-fJetRadius) continue;
	std::vector<fastjet::PseudoJet> Constituents(fFastJetWrapper->GetJetConstituents(i_Jet));
	//	NTracks=Constituents.size();
	for (UInt_t i_Constituents = 0; i_Constituents < Constituents.size(); i_Constituents++) { 
	  if (Constituents[i_Constituents].user_index() == 0) { 
	    Is_D_Jet = kTRUE; 
	  }
	}

      
	if (!Is_D_Jet) continue; 
	fhEvent->Fill(22); 
	
	
	Int_t i_Matched_D_Jet_Truth=-1;
	for (UInt_t k=0; k< Inclusive_Jets_Truth_Labels.size(); k++){
	  if(Inclusive_Jets_Truth_Labels[k].first==D_Candidate_MatchedTruth_Label) i_Matched_D_Jet_Truth=Inclusive_Jets_Truth_Labels[k].second; 
	}
	
	for (UInt_t i_Jet_Truth=0; i_Jet_Truth < Inclusive_Jets_Truth.size(); i_Jet_Truth++){
	  Bool_t Is_Jet_Truth_Matched=kFALSE;
	  if (TMath::Abs(Inclusive_Jets_Truth[i_Jet_Truth].pseudorapidity()) > 0.9-fJetRadius) continue;
	  //Jet_Eta_Truth=TMath::Abs(Inclusive_Jets_Truth[i_Jet_Truth].pseudorapidity());
	  std::vector<fastjet::PseudoJet> Constituents_Truth(fFastJetWrapper_Truth->GetJetConstituents(i_Jet_Truth));
	  // NTracks_Truth=Constituents_Truth.size();
	  for (UInt_t i_Constituents_Truth = 0; i_Constituents_Truth < Constituents_Truth.size(); i_Constituents_Truth++) { 
	    if (Constituents_Truth[i_Constituents_Truth].user_index() == i_Matched_D_Jet_Truth) { 
	      Is_Jet_Truth_Matched=kTRUE;
	      // Dmeson_Eta_Truth=TMath::Abs(Constituents_Truth[i_Constituents_Truth].pseudorapidity());
	      //Dmeson_Y_Truth=TMath::Abs(Constituents_Truth[i_Constituents_Truth].rapidity());
	      for (UInt_t i_Unmacthed_D=0; i_Unmacthed_D<Unmatched_Truth_Level_D.size(); i_Unmacthed_D++){ 
		if (Unmatched_Truth_Level_D[i_Unmacthed_D]==i_Matched_D_Jet_Truth) Unmatched_Truth_Level_D.erase(Unmatched_Truth_Level_D.begin()+i_Unmacthed_D); 
	      }
	    }
	  }
	  if (!Is_Jet_Truth_Matched) continue; 
	  fhEvent->Fill(23);
	  NMatched_DMeson_Jets++;
		

	
	  
	  std::vector<Double_t> Splittings_Zg;
	  std::vector<Double_t> Splittings_DeltaR;
	  // std::vector<Double_t> Splittings_LeadingSubJetpT;
	  std::vector<Double_t> Splittings_HardestSubJetD0;
	  std::vector<Double_t> Splittings_RadiatorE;
	  std::vector<Double_t> Splittings_RadiatorpT;

	  Bool_t Is_D_SubJet=kFALSE;
	  fastjet::JetDefinition Jet_Definition(fastjet::cambridge_algorithm, fJetRadius*2.5,static_cast<fastjet::RecombinationScheme>(0), fastjet::Best); 
	
	  try{
	    std::vector<fastjet::PseudoJet> Reclustered_Particles(fFastJetWrapper->GetJetConstituents(i_Jet));
	    fastjet::ClusterSequence Cluster_Sequence_CA(Reclustered_Particles, Jet_Definition);
	    std::vector<fastjet::PseudoJet> Reclustered_Jet =  Cluster_Sequence_CA.inclusive_jets(0.0);
	    Reclustered_Jet = sorted_by_pt(Reclustered_Jet);
         

	    fastjet::PseudoJet Daughter_Jet = Reclustered_Jet[0];
	    fastjet::PseudoJet Parent_SubJet_1; 
	    fastjet::PseudoJet Parent_SubJet_2;  
	  
	    while(Daughter_Jet.has_parents(Parent_SubJet_1,Parent_SubJet_2)){
	      if(Parent_SubJet_1.perp() < Parent_SubJet_2.perp()) std::swap(Parent_SubJet_1,Parent_SubJet_2);
	      //  Splittings_LeadingSubJetpT.push_back(Parent_SubJet_1.perp());
	      vector < fastjet::PseudoJet > Hard_SubJet_Constituents = sorted_by_pt(Parent_SubJet_1.constituents());
	      Is_D_SubJet=kFALSE;
	      for(UInt_t j=0;j<Hard_SubJet_Constituents.size();j++){
		if(Hard_SubJet_Constituents[j].user_index()==0) Is_D_SubJet=kTRUE;	      
	      }

	    
	      if (!Is_D_SubJet) Splittings_HardestSubJetD0.push_back(1.0); 
	      else Splittings_HardestSubJetD0.push_back(2.0);  

	      // if(!Is_D_SubJet) std::swap(Parent_SubJet_1,Parent_SubJet_2); 	    
	      Splittings_DeltaR.push_back(Parent_SubJet_1.delta_R(Parent_SubJet_2));
	      Splittings_Zg.push_back(Parent_SubJet_2.perp()/(Parent_SubJet_1.perp()+Parent_SubJet_2.perp()));
	      Splittings_RadiatorE.push_back(Daughter_Jet.E());
	      Splittings_RadiatorpT.push_back(Daughter_Jet.perp());
	      Daughter_Jet=Parent_SubJet_1;
	    }

         
	  } catch (fastjet::Error) { /*return -1;*/ }

	  std::vector<Double_t> Splittings_Zg_Truth;
	  std::vector<Double_t> Splittings_DeltaR_Truth;
	  //  std::vector<Double_t> Splittings_LeadingSubJetpT_Truth;
	  std::vector<Double_t> Splittings_HardestSubJetD0_Truth;
	  std::vector<Double_t> Splittings_RadiatorE_Truth;
	  std::vector<Double_t> Splittings_RadiatorpT_Truth;


	  Bool_t Is_D_SubJet_Truth=kFALSE;
	  fastjet::JetDefinition Jet_Definition_Truth(fastjet::cambridge_algorithm, fJetRadius*2.5,static_cast<fastjet::RecombinationScheme>(0), fastjet::Best); 
  
	  try{
	    std::vector<fastjet::PseudoJet> Reclustered_Particles_Truth(fFastJetWrapper_Truth->GetJetConstituents(i_Jet_Truth));
	    fastjet::ClusterSequence Cluster_Sequence_CA_Truth(Reclustered_Particles_Truth, Jet_Definition_Truth);
	    std::vector<fastjet::PseudoJet> Reclustered_Jet_Truth =  Cluster_Sequence_CA_Truth.inclusive_jets(0.0);
	    Reclustered_Jet_Truth = sorted_by_pt(Reclustered_Jet_Truth); 
         

	    fastjet::PseudoJet Daughter_Jet_Truth = Reclustered_Jet_Truth[0];
	    fastjet::PseudoJet Parent_SubJet_1_Truth; 
	    fastjet::PseudoJet Parent_SubJet_2_Truth;  
    
	    while(Daughter_Jet_Truth.has_parents(Parent_SubJet_1_Truth,Parent_SubJet_2_Truth)){
	      if(Parent_SubJet_1_Truth.perp() < Parent_SubJet_2_Truth.perp()) std::swap(Parent_SubJet_1_Truth,Parent_SubJet_2_Truth);
	      // Splittings_LeadingSubJetpT_Truth.push_back(Parent_SubJet_1_Truth.perp());
	      vector < fastjet::PseudoJet > Hard_SubJet_Constituents_Truth = sorted_by_pt(Parent_SubJet_1_Truth.constituents());
	      Is_D_SubJet_Truth=kFALSE;
	      for(UInt_t j=0;j<Hard_SubJet_Constituents_Truth.size();j++){
		if(Hard_SubJet_Constituents_Truth[j].user_index()==i_Matched_D_Jet_Truth) Is_D_SubJet_Truth=kTRUE; 	      
	      }

	    
	      if (!Is_D_SubJet_Truth) Splittings_HardestSubJetD0_Truth.push_back(1.0); 
	      else Splittings_HardestSubJetD0_Truth.push_back(2.0);  

	      // if(!Is_D_SubJet_Truth) std::swap(Parent_SubJet_1_Truth,Parent_SubJet_2_Truth);     
	      Splittings_DeltaR_Truth.push_back(Parent_SubJet_1_Truth.delta_R(Parent_SubJet_2_Truth));
	      Splittings_Zg_Truth.push_back(Parent_SubJet_2_Truth.perp()/(Parent_SubJet_1_Truth.perp()+Parent_SubJet_2_Truth.perp()));
	      Splittings_RadiatorE_Truth.push_back(Daughter_Jet_Truth.E());
	      Splittings_RadiatorpT_Truth.push_back(Daughter_Jet_Truth.perp());
	      Daughter_Jet_Truth=Parent_SubJet_1_Truth;
	    }

         
	  } catch (fastjet::Error) { /*return -1;*/ }


	  //set detector level flags
	  Double_t Flag_D=-1.0;
	  if (Mass_Hypo_Type ==1) Flag_D=1.0; 
	  else if (Mass_Hypo_Type ==2) Flag_D=2.0;
	  else if (Mass_Hypo_Type ==3 && Matched_Truth_Particle->GetPdgCode()==fCandidatePDG) Flag_D=3.0; 
	  else if (Mass_Hypo_Type ==3 && Matched_Truth_Particle->GetPdgCode()==-fCandidatePDG) Flag_D=4.0; 
	  
	  Double_t Flag_D_Truth=-1.0;
	  if (Matched_Truth_Particle->GetPdgCode()==fCandidatePDG) Flag_D_Truth=1.0;
	  if (Matched_Truth_Particle->GetPdgCode()==-fCandidatePDG) Flag_D_Truth=2.0;

	
	  fShapesVar[0] = Inclusive_Jets[i_Jet].perp();
	  fShapesVar[1] = Inclusive_Jets_Truth[i_Jet_Truth].perp();
	  fShapesVar[2] = D_Candidate->Pt();
	  fShapesVar[3] = Matched_Truth_Particle->Pt();
	  fShapesVar[4] = Inv_Mass_D;
	  fShapesVar[5] = Inv_Mass_D_Truth;
	  fShapesVar[6] = Flag_D;
	  fShapesVar[7] = Flag_D_Truth;
	  fShapesVar[8] = Is_Prompt_Correct_Quark_PDG;
	  fShapesVar[9] = 0.0;
	  // fShapesVar[10] = NTracks;
	  //fShapesVar[11] = NTracks_Truth;
	  //fShapesVar[12] = Jet_Eta;
	  //fShapesVar[13] = Jet_Eta_Truth;
	  //fShapesVar[14] = Dmeson_Eta;
	  //fShapesVar[15] = Dmeson_Eta_Truth;
	  //fShapesVar[16] = Dmeson_Y;
	  //fShapesVar[17] = Dmeson_Y_Truth;



	  fShapesVar_Splittings_DeltaR.push_back(Splittings_DeltaR);
	  fShapesVar_Splittings_DeltaR_Truth.push_back(Splittings_DeltaR_Truth); 
	  fShapesVar_Splittings_Zg.push_back(Splittings_Zg);
	  fShapesVar_Splittings_Zg_Truth.push_back(Splittings_Zg_Truth); 
	  //  fShapesVar_Splittings_LeadingSubJetpT.push_back(Splittings_LeadingSubJetpT);
	  //fShapesVar_Splittings_LeadingSubJetpT_Truth.push_back(Splittings_LeadingSubJetpT_Truth); 
	  fShapesVar_Splittings_HardestSubJetD0.push_back(Splittings_HardestSubJetD0);
	  fShapesVar_Splittings_HardestSubJetD0_Truth.push_back(Splittings_HardestSubJetD0_Truth);
	  fShapesVar_Splittings_RadiatorE.push_back(Splittings_RadiatorE);
	  fShapesVar_Splittings_RadiatorE_Truth.push_back(Splittings_RadiatorE_Truth);
	  fShapesVar_Splittings_RadiatorpT.push_back(Splittings_RadiatorpT);
	  fShapesVar_Splittings_RadiatorpT_Truth.push_back(Splittings_RadiatorpT_Truth); 

	
	  fTreeResponseMatrixAxis->Fill();
	  fTreeSplittings->Fill();

	  fShapesVar_Splittings_DeltaR.clear();
	  fShapesVar_Splittings_DeltaR_Truth.clear();
	  fShapesVar_Splittings_Zg.clear();
	  fShapesVar_Splittings_Zg_Truth.clear();
	  // fShapesVar_Splittings_LeadingSubJetpT.clear();
	  //fShapesVar_Splittings_LeadingSubJetpT_Truth.clear();
	  fShapesVar_Splittings_HardestSubJetD0.clear();
	  fShapesVar_Splittings_HardestSubJetD0_Truth.clear();
	  fShapesVar_Splittings_RadiatorE.clear();
	  fShapesVar_Splittings_RadiatorE_Truth.clear();
	  fShapesVar_Splittings_RadiatorpT.clear();
	  fShapesVar_Splittings_RadiatorpT_Truth.clear(); 
	}
      }
      // delete Matched_Truth_Particle;
      // delete D_Candidate;
    }
    // delete Track_Container;
    if(NMatched_DMeson_Jets==0) fhEvent->Fill(24); 
    if(NMatched_DMeson_Jets==1) fhEvent->Fill(25); 
    if(NMatched_DMeson_Jets==2) fhEvent->Fill(26); 
    if(NMatched_DMeson_Jets==3) fhEvent->Fill(27); 
    if(NMatched_DMeson_Jets==4) fhEvent->Fill(28); 
    if(NMatched_DMeson_Jets==5) fhEvent->Fill(29); 
    if(NMatched_DMeson_Jets==6) fhEvent->Fill(30); 
    
    if (fIncludeInclusive){ 

      if(Unmatched_Truth_Level_D.size()==0) fhEvent->Fill(31); 
      if(Unmatched_Truth_Level_D.size()==1) fhEvent->Fill(32); 
      if(Unmatched_Truth_Level_D.size()==2) fhEvent->Fill(33); 
      if(Unmatched_Truth_Level_D.size()==3) fhEvent->Fill(34); 
      if(Unmatched_Truth_Level_D.size()==4) fhEvent->Fill(35); 
      if(Unmatched_Truth_Level_D.size()==5) fhEvent->Fill(36); 
      if(Unmatched_Truth_Level_D.size()==6) fhEvent->Fill(37); 
      
      for (UInt_t i_Jet_Truth=0; i_Jet_Truth < Inclusive_Jets_Truth.size(); i_Jet_Truth++){ 
	AliAODMCParticle* Truth_D_Particle = NULL; 
	Bool_t Is_Unmatched_D=kFALSE;
	Int_t D_Meson_Matched_Index=-1;
	if (TMath::Abs(Inclusive_Jets_Truth[i_Jet_Truth].pseudorapidity()) > 0.9-fJetRadius) continue;
	//Jet_Eta_Truth=TMath::Abs(Inclusive_Jets_Truth[i_Jet_Truth].pseudorapidity());
	if (Inclusive_Jets_Truth[i_Jet_Truth].perp()<fJetMinPt) continue;
	std::vector<fastjet::PseudoJet> Constituents_Truth(fFastJetWrapper_Truth->GetJetConstituents(i_Jet_Truth));
	//NTracks_Truth=Constituents_Truth.size();
	for (UInt_t i_Constituents_Truth = 0; i_Constituents_Truth < Constituents_Truth.size(); i_Constituents_Truth++) { 
	  for (UInt_t i_Unmacthed_D=0; i_Unmacthed_D<Unmatched_Truth_Level_D.size(); i_Unmacthed_D++){ 
	    if(Constituents_Truth[i_Constituents_Truth].user_index() == Unmatched_Truth_Level_D[i_Unmacthed_D]){
	      Is_Unmatched_D=kTRUE;
	      D_Meson_Matched_Index=Constituents_Truth[i_Constituents_Truth].user_index();
	      for(UInt_t i_MC_Label=0; i_MC_Label<Inclusive_Jets_Truth_Labels.size(); i_MC_Label++){ 
		if (Inclusive_Jets_Truth_Labels[i_MC_Label].second==Constituents_Truth[i_Constituents_Truth].user_index()){
		  Truth_D_Particle=static_cast<AliAODMCParticle*>(Particle_Container->GetArray()->At(Inclusive_Jets_Truth_Labels[i_MC_Label].first));
		  // Dmeson_Eta_Truth=TMath::Abs(Truth_D_Particle->Eta());
		  //Dmeson_Y_Truth=TMath::Abs(Truth_D_Particle->Y());
		}
	      }
	    }
	  }
	}
	if (!Is_Unmatched_D) continue; 
	fhEvent->Fill(38);


	Bool_t Is_Prompt_Correct_Quark=kFALSE;
	Int_t Is_Prompt_Correct_Quark_PDG=-1;
	Int_t Truth_D_Particle_Mother_Label=Truth_D_Particle->GetMother(); 
	while (Truth_D_Particle_Mother_Label >= 0) {
	  AliAODMCParticle* Truth_D_Particle_Mother = static_cast<AliAODMCParticle*>(Particle_Container->GetArray()->At(Truth_D_Particle_Mother_Label)); 
	  if (Truth_D_Particle_Mother){
	    Int_t Original_Quark_PDG=4; 
	    if (fIsBDecay) Original_Quark_PDG=5; 
	    if (TMath::Abs(Truth_D_Particle_Mother->GetPdgCode()) == 4){
	      fhEvent->Fill(39);		
	      Is_Prompt_Correct_Quark_PDG=4;
	      break;
	    }
	    if (TMath::Abs(Truth_D_Particle_Mother->GetPdgCode()) == 5){
	      fhEvent->Fill(40);
	      Is_Prompt_Correct_Quark_PDG=5;
	      break;
	    }
	    if (TMath::Abs(Truth_D_Particle_Mother->GetPdgCode())==Original_Quark_PDG) Is_Prompt_Correct_Quark=kTRUE; 
	    if (Truth_D_Particle_Mother_Label==Truth_D_Particle_Mother->GetMother()) break; 
	    Truth_D_Particle_Mother_Label=Truth_D_Particle_Mother->GetMother();
	  }
	  else break;
	  // delete Truth_D_Particle_Mother;
	}

	if (fPromptReject && !Is_Prompt_Correct_Quark) continue; 


	Bool_t Is_D_SubJet_Truth=kFALSE;
	fastjet::JetDefinition Jet_Definition_Truth(fastjet::cambridge_algorithm, fJetRadius*2.5,static_cast<fastjet::RecombinationScheme>(0), fastjet::Best);

	std::vector<Double_t> Splittings_Zg_Truth;
	std::vector<Double_t> Splittings_DeltaR_Truth;
	//	std::vector<Double_t> Splittings_LeadingSubJetpT_Truth;
	std::vector<Double_t> Splittings_HardestSubJetD0_Truth;
	std::vector<Double_t> Splittings_RadiatorE_Truth;
	std::vector<Double_t> Splittings_RadiatorpT_Truth;

  
	try{
	  std::vector<fastjet::PseudoJet> Reclustered_Particles_Truth(fFastJetWrapper_Truth->GetJetConstituents(i_Jet_Truth));
	  fastjet::ClusterSequence Cluster_Sequence_CA_Truth(Reclustered_Particles_Truth, Jet_Definition_Truth);
	  std::vector<fastjet::PseudoJet> Reclustered_Jet_Truth =  Cluster_Sequence_CA_Truth.inclusive_jets(0.0);
	  Reclustered_Jet_Truth = sorted_by_pt(Reclustered_Jet_Truth); 

	  
	  fastjet::PseudoJet Daughter_Jet_Truth = Reclustered_Jet_Truth[0];
	  fastjet::PseudoJet Parent_SubJet_1_Truth; 
	  fastjet::PseudoJet Parent_SubJet_2_Truth;  
    
	  while(Daughter_Jet_Truth.has_parents(Parent_SubJet_1_Truth,Parent_SubJet_2_Truth)){
	    if(Parent_SubJet_1_Truth.perp() < Parent_SubJet_2_Truth.perp()) std::swap(Parent_SubJet_1_Truth,Parent_SubJet_2_Truth);
	    // Splittings_LeadingSubJetpT_Truth.push_back(Parent_SubJet_1_Truth.perp());
	    vector < fastjet::PseudoJet > Hard_SubJet_Constituents_Truth = sorted_by_pt(Parent_SubJet_1_Truth.constituents()); 
	    Is_D_SubJet_Truth=kFALSE;
	    for(UInt_t j=0;j<Hard_SubJet_Constituents_Truth.size();j++){
	      if(Hard_SubJet_Constituents_Truth[j].user_index()==D_Meson_Matched_Index) Is_D_SubJet_Truth=kTRUE;	      
	    }

	    
	    if (!Is_D_SubJet_Truth) Splittings_HardestSubJetD0_Truth.push_back(1.0); 
	    else Splittings_HardestSubJetD0_Truth.push_back(2.0);  

	    // if(!Is_D_SubJet_Truth) std::swap(Parent_SubJet_1_Truth,Parent_SubJet_2_Truth);     
	    Splittings_DeltaR_Truth.push_back(Parent_SubJet_1_Truth.delta_R(Parent_SubJet_2_Truth));
	    Splittings_Zg_Truth.push_back(Parent_SubJet_2_Truth.perp()/(Parent_SubJet_1_Truth.perp()+Parent_SubJet_2_Truth.perp()));
	    Splittings_RadiatorE_Truth.push_back(Daughter_Jet_Truth.E());
	    Splittings_RadiatorpT_Truth.push_back(Daughter_Jet_Truth.perp());
	    Daughter_Jet_Truth=Parent_SubJet_1_Truth;
	  }

         
	} catch (fastjet::Error) { /*return -1;*/ }
	
	Double_t Inv_Mass_D_Truth=0.0;
	Double_t Flag_D_Truth=-1.0;
	if (Truth_D_Particle->GetPdgCode()==fCandidatePDG){
	  //Inv_Mass_D_Truth=Truth_D_Particle->InvMassD0();
	  Inv_Mass_D_Truth=0.0;
	  Flag_D_Truth=3.0;
	}
	if (Truth_D_Particle->GetPdgCode()==-fCandidatePDG){
	  // Inv_Mass_D_Truth=Truth_D_Particle->InvMassD0bar();
	  Inv_Mass_D_Truth=0.0;
	  Flag_D_Truth=4.0;
	}

	
	fShapesVar[0] = 0.0;
	fShapesVar[1] = Inclusive_Jets_Truth[i_Jet_Truth].perp();
	fShapesVar[2] = 0.0;
	fShapesVar[3] = Truth_D_Particle->Pt();
	fShapesVar[4] = 0.0;
	fShapesVar[5] = Inv_Mass_D_Truth;
	fShapesVar[6] = 0.0;
	fShapesVar[7] = Flag_D_Truth;
	fShapesVar[8] = Is_Prompt_Correct_Quark_PDG;
	fShapesVar[9] = 0.0;
	//	fShapesVar[10] = 0.0;
	//	fShapesVar[11] = NTracks_Truth;
	//	fShapesVar[12] = 0.0;
	//	fShapesVar[13] = Jet_Eta_Truth;
	//	fShapesVar[14] = 0.0;
	//	fShapesVar[15] = Dmeson_Eta_Truth;
	//	fShapesVar[16] = 0.0;
	//	fShapesVar[17] = Dmeson_Y_Truth;

	fShapesVar_Splittings_DeltaR.push_back(Splittings_DeltaR_Truth);
	fShapesVar_Splittings_DeltaR_Truth.push_back(Splittings_DeltaR_Truth); 
	fShapesVar_Splittings_Zg.push_back(Splittings_Zg_Truth);
	fShapesVar_Splittings_Zg_Truth.push_back(Splittings_Zg_Truth); 
	//	fShapesVar_Splittings_LeadingSubJetpT.push_back(Splittings_LeadingSubJetpT_Truth);
	//	fShapesVar_Splittings_LeadingSubJetpT_Truth.push_back(Splittings_LeadingSubJetpT_Truth); 
	fShapesVar_Splittings_HardestSubJetD0.push_back(Splittings_HardestSubJetD0_Truth); 
	fShapesVar_Splittings_HardestSubJetD0_Truth.push_back(Splittings_HardestSubJetD0_Truth);
	fShapesVar_Splittings_RadiatorE.push_back(Splittings_RadiatorE_Truth);
	fShapesVar_Splittings_RadiatorE_Truth.push_back(Splittings_RadiatorE_Truth);
	fShapesVar_Splittings_RadiatorpT.push_back(Splittings_RadiatorpT_Truth);
	fShapesVar_Splittings_RadiatorpT_Truth.push_back(Splittings_RadiatorpT_Truth); 

	
	fTreeResponseMatrixAxis->Fill();
	fTreeSplittings->Fill();
	
	fShapesVar_Splittings_DeltaR.clear();
	fShapesVar_Splittings_DeltaR_Truth.clear();
	fShapesVar_Splittings_Zg.clear();
	fShapesVar_Splittings_Zg_Truth.clear(); 
	//	fShapesVar_Splittings_LeadingSubJetpT.clear();
	//	fShapesVar_Splittings_LeadingSubJetpT_Truth.clear();
	fShapesVar_Splittings_HardestSubJetD0.clear();
	fShapesVar_Splittings_HardestSubJetD0_Truth.clear();
	fShapesVar_Splittings_RadiatorE.clear();
	fShapesVar_Splittings_RadiatorE_Truth.clear();
	fShapesVar_Splittings_RadiatorpT.clear();
	fShapesVar_Splittings_RadiatorpT_Truth.clear();
	// delete Truth_D_Particle;		
      }	 
    }
    // delete Particle_Container;
  }


  if (fJetShapeType == kTrue){ 


    //truth level only jet finding

    Double_t NTracks_Truth=0;
    Double_t Jet_Eta_Truth=-5.0;
    Double_t Dmeson_Eta_Truth=-5.0;
    Double_t Dmeson_Y_Truth=-5.0;
    
    AliHFAODMCParticleContainer *Particle_Container = (AliHFAODMCParticleContainer*) GetParticleContainer(0);
    Particle_Container->SetSpecialPDG(fCandidatePDG); 
    Particle_Container->SetRejectedOriginMap(fRejectedOrigin); 
    Particle_Container->SetAcceptedDecayMap(fAcceptedDecay); 
    Particle_Container->SetRejectISR(fRejectISR);
    Particle_Container->SetCharge(AliParticleContainer::EChargeCut_t::kCharged);

    std::vector<fastjet::PseudoJet> Inclusive_Jets_Truth; 
    std::vector<std::pair<Int_t, Int_t>> Inclusive_Jets_Truth_Labels; 
    
    if (Particle_Container->IsSpecialPDGFound()){ 
      fhEvent->Fill(2); 

      fFastJetWrapper_Truth->SetAreaType(fastjet::active_area); 
      fFastJetWrapper_Truth->SetGhostArea(0.005);
      fFastJetWrapper_Truth->SetR(fJetRadius);
      fFastJetWrapper_Truth->SetAlgorithm(fastjet::antikt_algorithm);
      fFastJetWrapper_Truth->SetRecombScheme(static_cast<fastjet::RecombinationScheme>(1));
      fFastJetWrapper_Truth->Clear();
      
      AliAODMCParticle* Truth_Particle=NULL;
      Int_t NTruthD=0;
      for (Int_t i_Particle=0; i_Particle<Particle_Container->GetNParticles(); i_Particle++){ 
	Truth_Particle = static_cast<AliAODMCParticle*>(Particle_Container->GetAcceptMCParticle(i_Particle));
	if (!Truth_Particle) continue;
	if (TMath::Abs(Truth_Particle->GetPdgCode())==fCandidatePDG){ 
	  fhEvent->Fill(3); 	  
	  std::pair<Int_t, Int_t> Inclusive_Jet_Truth_Labels;
	  Inclusive_Jet_Truth_Labels.first=Truth_Particle->GetLabel(); 
	  Inclusive_Jet_Truth_Labels.second=NTruthD; 
	  Inclusive_Jets_Truth_Labels.push_back(Inclusive_Jet_Truth_Labels);
	  fFastJetWrapper_Truth->AddInputVector(Truth_Particle->Px(), Truth_Particle->Py(), Truth_Particle->Pz(), Truth_Particle->E(),NTruthD);
	  NTruthD++;
	}
	else fFastJetWrapper_Truth->AddInputVector(Truth_Particle->Px(), Truth_Particle->Py(), Truth_Particle->Pz(), Truth_Particle->E(),i_Particle+100); 
      }
      // delete Truth_Particle;
      fFastJetWrapper_Truth->Run();
      Inclusive_Jets_Truth = fFastJetWrapper_Truth->GetInclusiveJets(); 

      for (UInt_t i_Jet_Truth=0; i_Jet_Truth < Inclusive_Jets_Truth.size(); i_Jet_Truth++){ 
	AliAODMCParticle* Truth_D_Particle = NULL; 
	Bool_t Is_DJet_Truth=kFALSE;
	Int_t D_Meson_Matched_Index=-1;
	//if (TMath::Abs(Inclusive_Jets_Truth[i_Jet_Truth].pseudorapidity()) > 0.9-fJetRadius) continue;
	Jet_Eta_Truth=TMath::Abs(Inclusive_Jets_Truth[i_Jet_Truth].pseudorapidity());
	std::vector<fastjet::PseudoJet> Constituents_Truth(fFastJetWrapper_Truth->GetJetConstituents(i_Jet_Truth));
	NTracks_Truth=Constituents_Truth.size();
	for (UInt_t i_Constituents_Truth = 0; i_Constituents_Truth < Constituents_Truth.size(); i_Constituents_Truth++) { 
	  
	  if(Constituents_Truth[i_Constituents_Truth].user_index() >=0 && Constituents_Truth[i_Constituents_Truth].user_index() < NTruthD){
	    D_Meson_Matched_Index=Constituents_Truth[i_Constituents_Truth].user_index();
	    Is_DJet_Truth=kTRUE;
	    for(UInt_t i_MC_Label=0; i_MC_Label<Inclusive_Jets_Truth_Labels.size(); i_MC_Label++){
	      if (Inclusive_Jets_Truth_Labels[i_MC_Label].second==Constituents_Truth[i_Constituents_Truth].user_index()){
		Truth_D_Particle=static_cast<AliAODMCParticle*>(Particle_Container->GetArray()->At(Inclusive_Jets_Truth_Labels[i_MC_Label].first));
		Dmeson_Eta_Truth=TMath::Abs(Truth_D_Particle->Eta());
		Dmeson_Y_Truth=TMath::Abs(Truth_D_Particle->Y());
	      }
	    }
		
	  }
	}

	if (!Is_DJet_Truth) fhEvent->Fill(4); 
	if (Is_DJet_Truth) fhEvent->Fill(5);

	if (!fIncludeInclusive && !Is_DJet_Truth) continue; 

	Bool_t Is_Prompt_Correct_Quark=kFALSE;
	Int_t Is_Prompt_Correct_Quark_PDG=-1;
	Int_t Truth_D_Particle_Mother_Label=Truth_D_Particle->GetMother(); 
	while (Truth_D_Particle_Mother_Label >= 0) {
	  AliAODMCParticle* Truth_D_Particle_Mother = static_cast<AliAODMCParticle*>(Particle_Container->GetArray()->At(Truth_D_Particle_Mother_Label)); 
	  if (Truth_D_Particle_Mother){
	    Int_t Original_Quark_PDG=4; 
	    if (fIsBDecay) Original_Quark_PDG=5; 
	    if (TMath::Abs(Truth_D_Particle_Mother->GetPdgCode()) ==Original_Quark_PDG) Is_Prompt_Correct_Quark=kTRUE; 
	    if (TMath::Abs(Truth_D_Particle_Mother->GetPdgCode()) ==4){
	      fhEvent->Fill(6);
	      Is_Prompt_Correct_Quark_PDG=4;
	      break;
	    }
	    if (TMath::Abs(Truth_D_Particle_Mother->GetPdgCode()) ==5){
	      fhEvent->Fill(7);
	      Is_Prompt_Correct_Quark_PDG=5;
	      break;
	    }
	    if (Truth_D_Particle_Mother_Label==Truth_D_Particle_Mother->GetMother()) break; 
	    Truth_D_Particle_Mother_Label=Truth_D_Particle_Mother->GetMother();
	  }
	  else break;
	  // delete Truth_D_Particle_Mother;
	}

	if(fPromptReject && !Is_Prompt_Correct_Quark) continue;
	fhEvent->Fill(8); 

	Bool_t Is_D_SubJet_Truth=kFALSE;
	fastjet::JetDefinition Jet_Definition_Truth(fastjet::cambridge_algorithm, fJetRadius*2.5,static_cast<fastjet::RecombinationScheme>(1), fastjet::Best); 

	std::vector<Double_t> Splittings_Zg_Truth;
	std::vector<Double_t> Splittings_DeltaR_Truth;
	//	std::vector<Double_t> Splittings_LeadingSubJetpT_Truth;
	std::vector<Double_t> Splittings_HardestSubJetD0_Truth;
	std::vector<Double_t> Splittings_RadiatorE_Truth;
	std::vector<Double_t> Splittings_RadiatorpT_Truth;
	
	try{
	  std::vector<fastjet::PseudoJet> Reclustered_Particles_Truth(fFastJetWrapper_Truth->GetJetConstituents(i_Jet_Truth));
	  fastjet::ClusterSequence Cluster_Sequence_CA_Truth(Reclustered_Particles_Truth, Jet_Definition_Truth);
	  std::vector<fastjet::PseudoJet> Reclustered_Jet_Truth =  Cluster_Sequence_CA_Truth.inclusive_jets(0.0);
	  Reclustered_Jet_Truth = sorted_by_pt(Reclustered_Jet_Truth); 
         

	  fastjet::PseudoJet Daughter_Jet_Truth = Reclustered_Jet_Truth[0];
	  fastjet::PseudoJet Parent_SubJet_1_Truth; 
	  fastjet::PseudoJet Parent_SubJet_2_Truth;  
    
	  while(Daughter_Jet_Truth.has_parents(Parent_SubJet_1_Truth,Parent_SubJet_2_Truth)){
	    if(Parent_SubJet_1_Truth.perp() < Parent_SubJet_2_Truth.perp()) std::swap(Parent_SubJet_1_Truth,Parent_SubJet_2_Truth);
	    // Splittings_LeadingSubJetpT_Truth.push_back(Parent_SubJet_1_Truth.perp());
	    vector < fastjet::PseudoJet > Hard_SubJet_Constituents_Truth = sorted_by_pt(Parent_SubJet_1_Truth.constituents()); 
	    Is_D_SubJet_Truth=kFALSE;
	    for(UInt_t j=0;j<Hard_SubJet_Constituents_Truth.size();j++){
	      if(Hard_SubJet_Constituents_Truth[j].user_index()==D_Meson_Matched_Index) Is_D_SubJet_Truth=kTRUE;	      
	    }

	    if (fIncludeInclusive && !Is_DJet_Truth) Splittings_HardestSubJetD0_Truth.push_back(0.0); 
	    else if (Is_DJet_Truth && !Is_D_SubJet_Truth) Splittings_HardestSubJetD0_Truth.push_back(1.0); 
	    else if (Is_DJet_Truth && Is_D_SubJet_Truth)Splittings_HardestSubJetD0_Truth.push_back(2.0); 

	    // if(!Is_D_SubJet_Truth) std::swap(Parent_SubJet_1_Truth,Parent_SubJet_2_Truth); 	    
	    Splittings_DeltaR_Truth.push_back(Parent_SubJet_1_Truth.delta_R(Parent_SubJet_2_Truth));
	    Splittings_Zg_Truth.push_back(Parent_SubJet_2_Truth.perp()/(Parent_SubJet_1_Truth.perp()+Parent_SubJet_2_Truth.perp()));
	    Splittings_RadiatorE_Truth.push_back(Daughter_Jet_Truth.E());
	    Splittings_RadiatorpT_Truth.push_back(Daughter_Jet_Truth.perp());
	    Daughter_Jet_Truth=Parent_SubJet_1_Truth;
	  }

         
	} catch (fastjet::Error) { /*return -1;*/ }

	Double_t Inv_Mass_D_Truth=0.0;
	Double_t Flag_D_Truth=-1.0;
	Double_t D_Pt=-1.0;
	if (fIncludeInclusive && !Is_DJet_Truth){
	  Flag_D_Truth=0.0; 
	}
	else if (Is_DJet_Truth){
	  if (Truth_D_Particle->GetPdgCode()==fCandidatePDG){
	    //Inv_Mass_D_Truth=Truth_D_Particle->InvMassD0();
	    Inv_Mass_D_Truth=0.0;
	    Flag_D_Truth=1.0; 
	    D_Pt=Truth_D_Particle->Pt();
	    fhEvent->Fill(9);
	  }
	  if (Truth_D_Particle->GetPdgCode()==-fCandidatePDG){
	    // Inv_Mass_D_Truth=Truth_D_Particle->InvMassD0bar();
	    Inv_Mass_D_Truth=0.0;
	    Flag_D_Truth=2.0; 
	    D_Pt=Truth_D_Particle->Pt();
	    fhEvent->Fill(10); 
	  }
	}

	
	fShapesVar[0] = 0.0;
	fShapesVar[1] = Inclusive_Jets_Truth[i_Jet_Truth].perp();
	fShapesVar[2] = 0.0;
	fShapesVar[3] = D_Pt;
	fShapesVar[4] = 0.0;
	fShapesVar[5] = Inv_Mass_D_Truth;
	fShapesVar[6] = 0.0;
	fShapesVar[7] = Flag_D_Truth;
	fShapesVar[8] = 0.0;
	fShapesVar[9] = Is_Prompt_Correct_Quark_PDG;
	//	fShapesVar[10] = 0.0;
	//fShapesVar[11] = NTracks_Truth;
	//fShapesVar[12] = 0.0;
	//fShapesVar[13] = Jet_Eta_Truth;
	//fShapesVar[14] = 0.0;
	//	fShapesVar[15] = Dmeson_Eta_Truth;
	//	fShapesVar[16] = 0.0;
	//	fShapesVar[17] = Dmeson_Y_Truth;
	


	fShapesVar_Splittings_DeltaR.push_back(Splittings_DeltaR_Truth);
	fShapesVar_Splittings_DeltaR_Truth.push_back(Splittings_DeltaR_Truth); 
	fShapesVar_Splittings_Zg.push_back(Splittings_Zg_Truth);
	fShapesVar_Splittings_Zg_Truth.push_back(Splittings_Zg_Truth); 
	//	fShapesVar_Splittings_LeadingSubJetpT.push_back(Splittings_LeadingSubJetpT_Truth);
	//fShapesVar_Splittings_LeadingSubJetpT_Truth.push_back(Splittings_LeadingSubJetpT_Truth); 
	fShapesVar_Splittings_HardestSubJetD0.push_back(Splittings_HardestSubJetD0_Truth);
	fShapesVar_Splittings_HardestSubJetD0_Truth.push_back(Splittings_HardestSubJetD0_Truth);
	fShapesVar_Splittings_RadiatorE.push_back(Splittings_RadiatorE_Truth);
	fShapesVar_Splittings_RadiatorE_Truth.push_back(Splittings_RadiatorE_Truth);
	fShapesVar_Splittings_RadiatorpT.push_back(Splittings_RadiatorpT_Truth);
	fShapesVar_Splittings_RadiatorpT_Truth.push_back(Splittings_RadiatorpT_Truth); 

	
	fTreeResponseMatrixAxis->Fill();
	fTreeSplittings->Fill();

	fShapesVar_Splittings_DeltaR.clear();
	fShapesVar_Splittings_DeltaR_Truth.clear();
	fShapesVar_Splittings_Zg.clear();
	fShapesVar_Splittings_Zg_Truth.clear();
	//	fShapesVar_Splittings_LeadingSubJetpT.clear();
	//fShapesVar_Splittings_LeadingSubJetpT_Truth.clear();
	fShapesVar_Splittings_HardestSubJetD0.clear();
	fShapesVar_Splittings_HardestSubJetD0_Truth.clear();
	fShapesVar_Splittings_RadiatorE.clear();
	fShapesVar_Splittings_RadiatorE_Truth.clear();
	fShapesVar_Splittings_RadiatorpT.clear();
	fShapesVar_Splittings_RadiatorpT_Truth.clear();
	

	//	delete Truth_D_Particle;
      }
    if (NTruthD==0) fhEvent->Fill(11); 
    if (NTruthD==1) fhEvent->Fill(12); 
    if (NTruthD==2) fhEvent->Fill(13); 
    if (NTruthD==3) fhEvent->Fill(14); 
    if (NTruthD==4) fhEvent->Fill(15); 
    if (NTruthD==5) fhEvent->Fill(16); 
    if (NTruthD==6) fhEvent->Fill(17);
    }
    //delete Particle_Container; 
  }




  if (fJetShapeType == kDataInclusive){

    TRandom3 Random;
    Random.SetSeed(0);
    Double_t Random_Number;

    AliHFTrackContainer* Track_Container=dynamic_cast<AliHFTrackContainer*>(GetTrackContainer(0));
    if (!Track_Container) return kTRUE;
    //Track_Container->SetDMesonCandidate(NULL);
    fFastJetWrapper->Clear();
    // Double_t NTracks=0;
    //Double_t Jet_Eta=-5.0;
    //Double_t HardestTrack_Eta=-5.0;
    //Double_t HardestTrack_Y=-5.0;
    Double_t HardestTrack_Pt=-5.0;
    AliAODTrack *Track = NULL;
    for (Int_t i_Track=0; i_Track<Track_Container->GetNTracks(); i_Track++){
      Track = static_cast<AliAODTrack*>(Track_Container->GetAcceptParticle(i_Track));
      if(!Track) continue;
      if(Track->Pt() > 100.0 || TMath::Abs(Track->Eta()) > 0.9) continue;
      Random_Number=Random.Rndm();
      if(Random_Number > fTrackingEfficiency) continue;
      fFastJetWrapper->AddInputVector(Track->Px(), Track->Py(), Track->Pz(), Track->E(),i_Track+100);
    }
    fFastJetWrapper->Run();
    std::vector<fastjet::PseudoJet> Inclusive_Jets = fFastJetWrapper->GetInclusiveJets();
    for (UInt_t i_Jet=0; i_Jet < Inclusive_Jets.size(); i_Jet++){
      if (Inclusive_Jets[i_Jet].perp()<fJetMinPt) continue;
      if (TMath::Abs(Inclusive_Jets[i_Jet].pseudorapidity()) > 0.9-fJetRadius) continue;
      // Jet_Eta=TMath::Abs(Inclusive_Jets[i_Jet].pseudorapidity());
      fhEvent->Fill(23);

      HardestTrack_Pt=-5.0;
      std::vector<fastjet::PseudoJet> Constituents(fFastJetWrapper->GetJetConstituents(i_Jet));
      // NTracks=Constituents.size();
      for (UInt_t i_Constituents = 0; i_Constituents < Constituents.size(); i_Constituents++) { 
	if (Constituents[i_Constituents].perp() > HardestTrack_Pt){
	  HardestTrack_Pt=Constituents[i_Constituents].perp();
	  //HardestTrack_Eta=TMath::Abs(Constituents[i_Constituents].pseudorapidity());
	  //HardestTrack_Y=TMath::Abs(Constituents[i_Constituents].rapidity());
	}
      }

      
      std::vector<Double_t> Splittings_Zg;
      std::vector<Double_t> Splittings_DeltaR;
      // std::vector<Double_t> Splittings_LeadingSubJetpT;
      std::vector<Double_t> Splittings_HardestSubJetD0;
      std::vector<Double_t> Splittings_RadiatorE;
      std::vector<Double_t> Splittings_RadiatorpT;

      fastjet::JetDefinition Jet_Definition(fastjet::cambridge_algorithm, fJetRadius*2.5,static_cast<fastjet::RecombinationScheme>(0), fastjet::Best);
  
      try{
	std::vector<fastjet::PseudoJet> Reclustered_Particles(fFastJetWrapper->GetJetConstituents(i_Jet));
	fastjet::ClusterSequence Cluster_Sequence_CA(Reclustered_Particles, Jet_Definition);
	std::vector<fastjet::PseudoJet> Reclustered_Jet =  Cluster_Sequence_CA.inclusive_jets(0.0);
	Reclustered_Jet = sorted_by_pt(Reclustered_Jet);
         

	fastjet::PseudoJet Daughter_Jet = Reclustered_Jet[0];
	fastjet::PseudoJet Parent_SubJet_1; 
	fastjet::PseudoJet Parent_SubJet_2;

	while(Daughter_Jet.has_parents(Parent_SubJet_1,Parent_SubJet_2)){
	  if(Parent_SubJet_1.perp() < Parent_SubJet_2.perp()) std::swap(Parent_SubJet_1,Parent_SubJet_2);
	  Double_t Leading_Track_pT=0.0;
	  vector <fastjet::PseudoJet> Hard_SubJet_Constituents = sorted_by_pt(Parent_SubJet_1.constituents());
	  for (Int_t i_Leading_Track=0; i_Leading_Track < Hard_SubJet_Constituents.size(); i_Leading_Track++){
	    if (Hard_SubJet_Constituents[i_Leading_Track].perp() > Leading_Track_pT) Leading_Track_pT=Hard_SubJet_Constituents[i_Leading_Track].perp();
	  }
	  // Splittings_LeadingSubJetpT.push_back(Parent_SubJet_1.perp());
	  Splittings_HardestSubJetD0.push_back(Leading_Track_pT);  
	  Splittings_DeltaR.push_back(Parent_SubJet_1.delta_R(Parent_SubJet_2));
	  Splittings_Zg.push_back(Parent_SubJet_2.perp()/(Parent_SubJet_1.perp()+Parent_SubJet_2.perp()));
	  Splittings_RadiatorE.push_back(Daughter_Jet.E());
	  Splittings_RadiatorpT.push_back(Daughter_Jet.perp());
	  Daughter_Jet=Parent_SubJet_1;
	}

         
      } catch (fastjet::Error) { /*return -1;*/ }


	
      fShapesVar[0] = Inclusive_Jets[i_Jet].perp();
      fShapesVar[1] = 0.0;
      fShapesVar[2] = HardestTrack_Pt;
      fShapesVar[3] = 0.0;
      fShapesVar[4] = 0.0;
      fShapesVar[5] = 0.0;
      fShapesVar[6] = 0.0;
      fShapesVar[7] = 0.0;
      fShapesVar[8] = 0.0;
      fShapesVar[9] = 0.0;
      // fShapesVar[10] = NTracks;
      //fShapesVar[11] = 0.0;
      //fShapesVar[12] = Jet_Eta;
      //fShapesVar[13] = 0.0;
      //fShapesVar[14] = HardestTrack_Eta;
      //fShapesVar[15] = 0.0;
      //fShapesVar[16] = HardestTrack_Y;
      //fShapesVar[17] = 0.0;


      fShapesVar_Splittings_DeltaR.push_back(Splittings_DeltaR);
      fShapesVar_Splittings_DeltaR_Truth.push_back(Splittings_DeltaR); 
      fShapesVar_Splittings_Zg.push_back(Splittings_Zg);
      fShapesVar_Splittings_Zg_Truth.push_back(Splittings_Zg);
      // fShapesVar_Splittings_LeadingSubJetpT.push_back(Splittings_LeadingSubJetpT);
      //fShapesVar_Splittings_LeadingSubJetpT_Truth.push_back(Splittings_LeadingSubJetpT); 
      fShapesVar_Splittings_HardestSubJetD0.push_back(Splittings_HardestSubJetD0);
      fShapesVar_Splittings_HardestSubJetD0_Truth.push_back(Splittings_HardestSubJetD0);
      fShapesVar_Splittings_RadiatorE.push_back(Splittings_RadiatorE);
      fShapesVar_Splittings_RadiatorE_Truth.push_back(Splittings_RadiatorE);
      fShapesVar_Splittings_RadiatorpT.push_back(Splittings_RadiatorpT);
      fShapesVar_Splittings_RadiatorpT_Truth.push_back(Splittings_RadiatorpT); 


      fTreeResponseMatrixAxis->Fill();
      fTreeSplittings->Fill();

      fShapesVar_Splittings_DeltaR.clear();
      fShapesVar_Splittings_DeltaR_Truth.clear();
      fShapesVar_Splittings_Zg.clear();
      fShapesVar_Splittings_Zg_Truth.clear();
      // fShapesVar_Splittings_LeadingSubJetpT.clear();
      //fShapesVar_Splittings_LeadingSubJetpT_Truth.clear();
      fShapesVar_Splittings_HardestSubJetD0.clear();
      fShapesVar_Splittings_HardestSubJetD0_Truth.clear();
      fShapesVar_Splittings_RadiatorE.clear();
      fShapesVar_Splittings_RadiatorE_Truth.clear();
      fShapesVar_Splittings_RadiatorpT.clear();
      fShapesVar_Splittings_RadiatorpT_Truth.clear();

    }
  }











  

  delete fFastJetWrapper_Truth;
  delete fFastJetWrapper;
  
  return kTRUE;

 
}






//________________________________________________________________________
Bool_t AliAnalysisTaskHFSubstructure::RetrieveEventObjects() {
  //
  // retrieve event objects
  //
  if (!AliAnalysisTaskEmcal::RetrieveEventObjects())
    return kFALSE;

  return kTRUE;
}


//_______________________________________________________________________
void AliAnalysisTaskHFSubstructure::Terminate(Option_t *) 
{
  // Called once at the end of the analysis.
  
}
