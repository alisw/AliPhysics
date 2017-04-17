/**************************************************************************
 * Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appeuear in the supporting documentation. The authors make no claims   *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id$ */

//
//
//               Lc->pK0s analysis code
//
//  Input: AOD
//  Output: TTree or THnSparse (mass vs pT vs Centrality)
//
//  Cuts:
//  TTree: very loose cut
//  THnSparse: One THnSparse is created per cut. One cut is specified by
//  an array of bits, each bit corresponds to a cut in "Cut" function.
//  Use "AddCutStream" function to add a cut. 
//
//-------------------------------------------------------------------------
//
//                 Authors: Y.S Watanabe(a)
//  (a) CNS, the University of Tokyo
//  Contatcs: wyosuke@cns.s.u-tokyo.ac.jp
//-------------------------------------------------------------------------

#include <TSystem.h>
#include <TParticle.h>
#include <TParticlePDG.h>
#include <TH1F.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TRandom.h>
#include <THnSparse.h>
#include <TLorentzVector.h>
#include <TTree.h>
#include "TROOT.h"
#include <TDatabasePDG.h>
#include <AliAnalysisDataSlot.h>
#include <AliAnalysisDataContainer.h>
#include "AliMCEvent.h"
#include "AliAnalysisManager.h"
#include "AliAODMCHeader.h"
#include "AliAODHandler.h"
#include "AliLog.h"
#include "AliExternalTrackParam.h"
#include "AliAODVertex.h"
#include "AliESDVertex.h"
#include "AliAODRecoDecay.h"
#include "AliAODRecoDecayHF.h"
#include "AliAODRecoCascadeHF.h"
#include "AliESDtrack.h"
#include "AliAODTrack.h"
#include "AliAODv0.h"
#include "AliAODcascade.h"
#include "AliAODMCParticle.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisTaskSELc2pK0sfromAODtracks.h"
#include "AliPIDResponse.h"
#include "AliPIDCombined.h"
#include "AliTOFPIDResponse.h"
#include "AliAODPidHF.h"
#include "AliInputEventHandler.h"
#include "AliESDtrackCuts.h"
#include "AliNeutralTrackParam.h"
#include "AliKFParticle.h"
#include "AliKFVertex.h"
#include "AliExternalTrackParam.h"
#include "AliESDtrack.h"
#include "AliCentrality.h"
#include "AliVertexerTracks.h"
#include "AliNormalizationCounter.h"
#include <set>

using std::cout;
using std::endl;

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskSELc2pK0sfromAODtracks);
/// \endcond

//__________________________________________________________________________
AliAnalysisTaskSELc2pK0sfromAODtracks::AliAnalysisTaskSELc2pK0sfromAODtracks() : 
  AliAnalysisTaskSE(),
  fUseMCInfo(kFALSE),
  fOutput(0),
  fOutputAll(0),
  fListCuts(0),
  fCEvents(0),
  fHTrigger(0),
  fHCentrality(0),
  fHReactionPlane(0),
  fAnalCuts(0),
  fIsEventSelected(kFALSE),
  fWriteVariableTree(kFALSE),
  fWriteEachVariableTree(kFALSE),
  fWriteMCVariableTree(kFALSE),
  fVariablesTree(0),
  fProtonVariablesTree(0),
  fV0VariablesTree(0),
  fMCVariablesTree(0),
  fMCProtonVariablesTree(0),
  fMCV0VariablesTree(0),
  fReconstructPrimVert(kFALSE),
  fIsMB(kFALSE),
  fIsSemi(kFALSE),
  fIsCent(kFALSE),
  fIsINT7(kFALSE),
  fIsEMC7(kFALSE),
  fCandidateVariables(),
  fCandidateProtonVariables(),
  fCandidateV0Variables(),
  fCandidateMCVariables(),
  fCandidateMCProtonVariables(),
  fCandidateMCV0Variables(),
  fVtx1(0),
  fV1(0),
  fVtxZ(0),
  fBzkG(0),
  fCentrality(0),
  fReactionPlane(0),
  fRunNumber(0),
  fTriggerCheck(0),
  fUseCentralityV0M(kFALSE),
  fEvNumberCounter(0),
	fCounter(0),
	fHistonEvtvsRunNumber(0),
	fHistonProtonvsRunNumber(0),
	fHistonK0svsRunNumber(0),
  fHistoLcMCGen(0),
  fHistoLcK0SpMass(0),
  fHistoLcK0SpMassMix(0),
  fHistoLcK0SpMassCoarse(0),
  fHistoLcK0SpMassMixCoarse(0),
  fHistoK0spCorrelation(0),
  fHistoK0spCorrelationMix(0),
  fHistoK0spCorrelationMCS(0),
  fHistoLcK0SpMassMCS(0),
  fHistoLcK0SpPi0MassMCS(0),
  fHistoLcKPluspMass(0),
  fHistoLcKMinuspMass(0),
  fHistoLcKPluspMassMix(0),
  fHistoLcKMinuspMassMix(0),
  fHistoBachPt(0),
  fHistoBachPtMCS(0),
  fHistoBachPtMCGen(0),
  fHistoKaonPt(0),
  fHistoKaonPtMCS(0),
  fHistoKaonPtMCGen(0),
  fHistoK0sMassvsPt(0),
  fHistoK0sMassvsPtMCS(0),
  fHistoK0sMassvsPtMCGen(0),
  fHistod0Bach(0),
  fHistod0V0(0),
  fHistod0d0(0),
  fHistoV0CosPA(0),
  fHistoProbProton(0),
  fHistoDecayLength(0),
  fHistoK0SMass(0),
  fHistoMassTagV0Min(0),
  fHistoMassTagV0SameSignMin(0),
  fHistoResponseLcPt(0),
  fHistoResponseLcPt1(0),
  fHistoResponseLcPt2(0),
  fGTI(0),fGTIndex(0), fTrackBuffSize(19000),
  fDoEventMixing(0),
	fNumberOfEventsForMixing		(5),
	fNzVtxBins					(0), 
	fNCentBins					(0),
	fNRPBins					(0),
	fNOfPools(1),
	fEventBuffer(0x0),
	fEventInfo(0x0),
	fProtonTracks(0x0),
	fV0Tracks(0x0),
  fProtonCutVarsArray(0x0),
  fV0CutVarsArray(0x0)
{
  //
  // Default Constructor. 
  //
}

//___________________________________________________________________________
AliAnalysisTaskSELc2pK0sfromAODtracks::AliAnalysisTaskSELc2pK0sfromAODtracks(const Char_t* name,
									     AliRDHFCutsLctopK0sfromAODtracks* analCuts, 
									     Bool_t writeVariableTree) :
  AliAnalysisTaskSE(name),
  fUseMCInfo(kFALSE),
  fOutput(0),
  fOutputAll(0),
  fListCuts(0),
  fCEvents(0),
  fHTrigger(0),
  fHCentrality(0),
  fHReactionPlane(0),
  fAnalCuts(analCuts),
  fIsEventSelected(kFALSE),
  fWriteVariableTree(writeVariableTree),
  fWriteEachVariableTree(kFALSE),
  fWriteMCVariableTree(kFALSE),
  fVariablesTree(0),
  fProtonVariablesTree(0),
  fV0VariablesTree(0),
  fMCVariablesTree(0),
  fMCProtonVariablesTree(0),
  fMCV0VariablesTree(0),
  fReconstructPrimVert(kFALSE),
  fIsMB(kFALSE),
  fIsSemi(kFALSE),
  fIsCent(kFALSE),
  fIsINT7(kFALSE),
  fIsEMC7(kFALSE),
  fCandidateVariables(),
  fCandidateProtonVariables(),
  fCandidateV0Variables(),
  fCandidateMCVariables(),
  fCandidateMCProtonVariables(),
  fCandidateMCV0Variables(),
  fVtx1(0),
  fV1(0),
  fVtxZ(0),
  fBzkG(0),
  fCentrality(0),
  fReactionPlane(0),
  fRunNumber(0),
  fTriggerCheck(0),
  fUseCentralityV0M(kFALSE),
  fEvNumberCounter(0),
	fCounter(0),
	fHistonEvtvsRunNumber(0),
	fHistonProtonvsRunNumber(0),
	fHistonK0svsRunNumber(0),
  fHistoLcMCGen(0),
  fHistoLcK0SpMass(0),
  fHistoLcK0SpMassMix(0),
  fHistoLcK0SpMassCoarse(0),
  fHistoLcK0SpMassMixCoarse(0),
  fHistoK0spCorrelation(0),
  fHistoK0spCorrelationMix(0),
  fHistoK0spCorrelationMCS(0),
  fHistoLcK0SpMassMCS(0),
  fHistoLcK0SpPi0MassMCS(0),
  fHistoLcKPluspMass(0),
  fHistoLcKMinuspMass(0),
  fHistoLcKPluspMassMix(0),
  fHistoLcKMinuspMassMix(0),
  fHistoBachPt(0),
  fHistoBachPtMCS(0),
  fHistoBachPtMCGen(0),
  fHistoKaonPt(0),
  fHistoKaonPtMCS(0),
  fHistoKaonPtMCGen(0),
  fHistoK0sMassvsPt(0),
  fHistoK0sMassvsPtMCS(0),
  fHistoK0sMassvsPtMCGen(0),
  fHistod0Bach(0),
  fHistod0V0(0),
  fHistod0d0(0),
  fHistoV0CosPA(0),
  fHistoProbProton(0),
  fHistoDecayLength(0),
  fHistoK0SMass(0),
  fHistoMassTagV0Min(0),
  fHistoMassTagV0SameSignMin(0),
  fHistoResponseLcPt(0),
  fHistoResponseLcPt1(0),
  fHistoResponseLcPt2(0),
  fGTI(0),fGTIndex(0), fTrackBuffSize(19000),
  fDoEventMixing(0),
	fNumberOfEventsForMixing		(5),
	fNzVtxBins					(0), 
	fNCentBins					(0),
	fNRPBins					(0),
	fNOfPools(1),
	fEventBuffer(0x0),
	fEventInfo(0x0),
	fProtonTracks(0x0),
	fV0Tracks(0x0),
  fProtonCutVarsArray(0x0),
  fV0CutVarsArray(0x0)
{
  //
  // Constructor. Initialization of Inputs and Outputs
  //
  Info("AliAnalysisTaskSELc2pK0sfromAODtracks","Calling Constructor");

  DefineOutput(1,TList::Class());  //conters
  DefineOutput(2,TList::Class());
  DefineOutput(3,TList::Class());  //conters
  DefineOutput(4,TTree::Class());  //My private output
  DefineOutput(5,TTree::Class());  //My private output
  DefineOutput(6,TTree::Class());  //My private output
  DefineOutput(7,TTree::Class());  //My private output
  DefineOutput(8,AliNormalizationCounter::Class());
  DefineOutput(9,TTree::Class());  //My private output
  DefineOutput(10,TTree::Class());  //My private output
}

//___________________________________________________________________________
AliAnalysisTaskSELc2pK0sfromAODtracks::~AliAnalysisTaskSELc2pK0sfromAODtracks() {
  //
  // destructor
  //
  Info("~AliAnalysisTaskSELc2pK0sfromAODtracks","Calling Destructor");

  if (fOutput) {
    delete fOutput;
    fOutput = 0;
  }

  if (fOutputAll) {
    delete fOutputAll;
    fOutputAll = 0;
  }

  if (fListCuts) {
    delete fListCuts;
    fListCuts = 0;
  }


  if (fAnalCuts) {
    delete fAnalCuts;
    fAnalCuts = 0;
  }

  if (fVariablesTree) {
    delete fVariablesTree;
    fVariablesTree = 0;
  }

  if (fProtonVariablesTree) {
    delete fProtonVariablesTree;
    fProtonVariablesTree = 0;
  }
  if (fV0VariablesTree) {
    delete fV0VariablesTree;
    fV0VariablesTree = 0;
  }
  if (fMCVariablesTree) {
    delete fMCVariablesTree;
    fMCVariablesTree = 0;
  }
  if (fMCProtonVariablesTree) {
    delete fMCProtonVariablesTree;
    fMCProtonVariablesTree = 0;
  }
  if (fMCV0VariablesTree) {
    delete fMCV0VariablesTree;
    fMCV0VariablesTree = 0;
  }

	if(fCounter){
		delete fCounter;
		fCounter = 0;
	}

	if(fProtonTracks) fProtonTracks->Delete();
	delete fProtonTracks;
	if(fV0Tracks) fV0Tracks->Delete();
	delete fV0Tracks;
	if(fProtonCutVarsArray) fProtonCutVarsArray->Delete();
	delete fProtonCutVarsArray;
	if(fV0CutVarsArray) fV0CutVarsArray->Delete();
	delete fV0CutVarsArray;
  if(fEventBuffer){
    for(Int_t i=0; i<fNOfPools; i++) delete fEventBuffer[i];
    delete fEventBuffer;
  }
  delete fEventInfo;

  if (fGTI)
    delete[] fGTI;
  fGTI=0;
  if (fGTIndex)
    delete[] fGTIndex;
  fGTIndex=0;

}

//_________________________________________________
void AliAnalysisTaskSELc2pK0sfromAODtracks::Init() {
  //
  // Initialization
  //
  //

  fIsEventSelected=kFALSE;

  if (fDebug > 1) AliInfo("Init");

  fListCuts = new TList();
  fListCuts->SetOwner();
  fListCuts->SetName("ListCuts");
  fListCuts->Add(new AliRDHFCutsLctopK0sfromAODtracks(*fAnalCuts));
  PostData(2,fListCuts);

  return;
}

//_________________________________________________
void AliAnalysisTaskSELc2pK0sfromAODtracks::UserExec(Option_t *)
{
  //
  // UserExec
  //

  if (!fInputEvent) {
    AliError("NO EVENT FOUND!");
    return;
  }
  AliAODEvent* aodEvent = dynamic_cast<AliAODEvent*>(fInputEvent);
  fCEvents->Fill(1);
	fEvNumberCounter++;

  //------------------------------------------------
  // First check if the event has proper vertex and B
  //------------------------------------------------
  fBzkG = (Double_t)aodEvent->GetMagneticField(); 
  AliKFParticle::SetField(fBzkG);
  if (TMath::Abs(fBzkG)<0.001) {
    delete fV1;
    return;
  }
  fCEvents->Fill(2);

  fCounter->StoreEvent(aodEvent,fAnalCuts,fUseMCInfo);
  fIsEventSelected = fAnalCuts->IsEventSelected(aodEvent); 

  //------------------------------------------------
  // MC analysis setting
  //------------------------------------------------
  TClonesArray *mcArray = 0;
  AliAODMCHeader *mcHeader=0;
  if (fUseMCInfo) {
    // MC array need for maching
    mcArray = dynamic_cast<TClonesArray*>(aodEvent->FindListObject(AliAODMCParticle::StdBranchName()));
    if (!mcArray) {
      AliError("Could not find Monte-Carlo in AOD");
      return;
    }
    fCEvents->Fill(6); // in case of MC events
  
    // load MC header
    mcHeader = (AliAODMCHeader*)aodEvent->GetList()->FindObject(AliAODMCHeader::StdBranchName());
    if (!mcHeader) {
      AliError("AliAnalysisTaskSELc2pK0sfromAODtracks::UserExec: MC header branch not found!\n");
      return;
    }
    fCEvents->Fill(7); // in case of MC events
  
    Double_t zMCVertex = mcHeader->GetVtxZ();
    if (TMath::Abs(zMCVertex) > fAnalCuts->GetMaxVtxZ()) {
      AliDebug(2,Form("Event rejected: abs(zVtxMC)=%f > fAnalCuts->GetMaxVtxZ()=%f",zMCVertex,fAnalCuts->GetMaxVtxZ()));
      return;
    } else {
      fCEvents->Fill(17); // in case of MC events
    }
    if ((TMath::Abs(zMCVertex) < fAnalCuts->GetMaxVtxZ()) && (!fAnalCuts->IsEventRejectedDuePhysicsSelection()) && (!fAnalCuts->IsEventRejectedDueToTrigger())) {
			Bool_t selevt = MakeMCAnalysis(mcArray);
			if(!selevt) return;
		}
  }

  //------------------------------------------------
  // Event selection 
  //------------------------------------------------
  fVtx1 = (AliAODVertex*)aodEvent->GetPrimaryVertex();
  if (!fVtx1) return;

  Double_t pos[3],cov[6];
  fVtx1->GetXYZ(pos);
  fVtx1->GetCovarianceMatrix(cov);
  fV1 = new AliESDVertex(pos,cov,100.,100,fVtx1->GetName());
	fVtxZ = pos[2];

  Bool_t fIsTriggerNotOK = fAnalCuts->IsEventRejectedDueToTrigger();
  if(!fIsTriggerNotOK) fCEvents->Fill(3);
  fIsEventSelected = fAnalCuts->IsEventSelected(aodEvent); // better to initialize before CheckEventSelection call
  if(!fIsEventSelected) {
    delete fV1;
    return;
  }
  fCEvents->Fill(4);

  fIsMB=(((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected()&AliVEvent::kMB)==(AliVEvent::kMB);
  fIsSemi=(((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected()&AliVEvent::kSemiCentral)==(AliVEvent::kSemiCentral);
  fIsCent=(((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected()&AliVEvent::kCentral)==(AliVEvent::kCentral); 
  fIsINT7=(((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected()&AliVEvent::kINT7)==(AliVEvent::kINT7);  
  fIsEMC7=(((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected()&AliVEvent::kEMC7)==(AliVEvent::kEMC7);   
  fTriggerCheck = fIsMB+2*fIsSemi+4*fIsCent+8*fIsINT7+16*fIsEMC7;
  if(fIsMB) fHTrigger->Fill(1);
  if(fIsSemi) fHTrigger->Fill(2);
  if(fIsCent) fHTrigger->Fill(3);
  if(fIsINT7) fHTrigger->Fill(4);
  if(fIsEMC7) fHTrigger->Fill(5);
  if(fIsMB|fIsSemi|fIsCent) fHTrigger->Fill(7);
  if(fIsINT7|fIsEMC7) fHTrigger->Fill(8);
  if(fIsMB&fIsSemi) fHTrigger->Fill(10);
  if(fIsMB&fIsCent) fHTrigger->Fill(11);
  if(fIsINT7&fIsEMC7) fHTrigger->Fill(12);

	if(fUseCentralityV0M){
		AliCentrality *cent = aodEvent->GetCentrality();
		fCentrality = cent->GetCentralityPercentile("V0M");
	}else{
		fCentrality = 1.;
	}
	if(fCentrality<0.||fCentrality>100.-0.0000001) {
		delete fV1;
		return;
	}
  fHCentrality->Fill(fCentrality);
	fRunNumber = aodEvent->GetRunNumber();

	Int_t runnumber_offset = 0;
	Int_t runnumber = aodEvent->GetRunNumber();
	if(runnumber<=131000&&runnumber>=114000){
		runnumber_offset = 114000;//lhc10bcde
	}else if(runnumber<=196000&&runnumber>=195000){
		runnumber_offset = 195000;//lhc13bc
	}else if(runnumber<=170593&&runnumber>=167902){
		runnumber_offset = 167902;//lhc11h
	}
	fHistonEvtvsRunNumber->Fill(runnumber-runnumber_offset,1.);

  //------------------------------------------------
  // Check if the event has v0 candidate
  //------------------------------------------------
  //Int_t nv0 = aodEvent->GetNumberOfV0s();
  fCEvents->Fill(5);

  //------------------------------------------------
  // Main analysis done in this function
  //------------------------------------------------
  fAnalCuts->SetMagneticField(fBzkG);
  fAnalCuts->SetPrimaryVertex(pos);
  MakeAnalysis(aodEvent,mcArray);


  PostData(1,fOutput);
  PostData(3,fOutputAll);
  PostData(4,fVariablesTree);
  PostData(5,fProtonVariablesTree);
  PostData(6,fV0VariablesTree);
  PostData(7,fMCVariablesTree);
  PostData(8,fCounter);    
  PostData(9,fMCProtonVariablesTree);
  PostData(10,fMCV0VariablesTree);

  fIsEventSelected=kFALSE;

  delete fV1;
  return;
}

//________________________________________ terminate ___________________________
void AliAnalysisTaskSELc2pK0sfromAODtracks::Terminate(Option_t*)
{    
  // The Terminate() function is the last function to be called during
  // a query. It always runs on the client, it can be used to present
  // the results graphically or save the results to file.
  
  //AliInfo("Terminate","");
  AliAnalysisTaskSE::Terminate();
  
  fOutput = dynamic_cast<TList*> (GetOutputData(1));
  if (!fOutput) {     
    AliError("fOutput not available");
    return;
  }

  fOutputAll = dynamic_cast<TList*> (GetOutputData(3));
  if (!fOutputAll) {     
    AliError("fOutputAll not available");
    return;
  }

  return;
}

//___________________________________________________________________________
void AliAnalysisTaskSELc2pK0sfromAODtracks::UserCreateOutputObjects() 
{ 
  //
  // UserCreateOutputObject
  //
  //AliInfo(Form("CreateOutputObjects of task %s\n", GetName()));

  //------------------------------------------------
  // output object setting
  //------------------------------------------------
  fOutput = new TList();
  fOutput->SetOwner();
  fOutput->SetName("chist0");
  DefineGeneralHistograms(); // define general histograms
  PostData(1,fOutput);

  fOutputAll = new TList();
  fOutputAll->SetOwner();
  fOutputAll->SetName("anahisto");
  DefineAnalysisHistograms(); // define general histograms
  PostData(3,fOutputAll);

  DefineTreeVariables();
  PostData(4,fVariablesTree);

  DefineProtonTreeVariables();
  PostData(5,fProtonVariablesTree);

  DefineV0TreeVariables();
  PostData(6,fV0VariablesTree);

  DefineMCTreeVariables();
  PostData(7,fMCVariablesTree);

  DefineMCProtonTreeVariables();
  PostData(9,fMCProtonVariablesTree);

  DefineMCV0TreeVariables();
  PostData(10,fMCV0VariablesTree);

  //Counter for Normalization
  TString normName="NormalizationCounter";
  AliAnalysisDataContainer *cont = GetOutputSlot(8)->GetContainer();
  if(cont)normName=(TString)cont->GetName();
  fCounter = new AliNormalizationCounter(normName.Data());
  fCounter->Init();
  PostData(8,fCounter);

	if(fDoEventMixing){
		fProtonTracks = new TObjArray();
		fProtonTracks->SetOwner();
		fV0Tracks = new TObjArray();
		fV0Tracks->SetOwner();
		fProtonCutVarsArray = new TObjArray();
		fProtonCutVarsArray->SetOwner();
		fV0CutVarsArray = new TObjArray();
		fV0CutVarsArray->SetOwner();

		fNOfPools=fNCentBins*fNzVtxBins*fNRPBins;
		fEventBuffer = new TTree*[fNOfPools];
		for(Int_t i=0; i<fNOfPools; i++){
			fEventBuffer[i]=new TTree(Form("EventBuffer_%d",i), "Temporary buffer for event mixing");
			fEventBuffer[i]->Branch("zVertex", &fVtxZ);
			fEventBuffer[i]->Branch("centrality", &fCentrality);
			fEventBuffer[i]->Branch("reactionplane", &fReactionPlane);
			fEventBuffer[i]->Branch("eventInfo", "TObjString",&fEventInfo);
			fEventBuffer[i]->Branch("v1array", "TObjArray", &fV0Tracks);
			fEventBuffer[i]->Branch("v1varsarray", "TObjArray", &fV0CutVarsArray);
		}
	}

  fGTI = new AliAODTrack *[fTrackBuffSize]; // Array of pointers 
  fGTIndex = new Int_t [fTrackBuffSize]; // Array of index 

  return;
}

//-------------------------------------------------------------------------------
void AliAnalysisTaskSELc2pK0sfromAODtracks::MakeAnalysis
(
 AliAODEvent *aodEvent, TClonesArray *mcArray
 )
{
  //
  // Main Analysis part
  //
	if(fDoEventMixing){
		if(fProtonTracks) fProtonTracks->Delete();
		if(fV0Tracks) fV0Tracks->Delete();
		if(fProtonCutVarsArray) fProtonCutVarsArray->Delete();
		if(fV0CutVarsArray) fV0CutVarsArray->Delete();
	}

  ResetGlobalTrackReference();
  // ..and set it
  for (Int_t iTrack=0;iTrack<aodEvent->GetNumberOfTracks();iTrack++){
    // cast needed since the event now returns AliVTrack instead of AliAODTrack
    AliAODTrack *track = dynamic_cast<AliAODTrack *>(aodEvent->GetTrack(iTrack));
    if (!track) continue;

    // Store the reference of the global tracks
    StoreGlobalTrackReference(track,iTrack);
  }


  Int_t nV0s= aodEvent->GetNumberOfV0s();
  Int_t nTracks= aodEvent->GetNumberOfTracks();

  Int_t  seleTrkFlags[nTracks];
  Int_t nSeleTrks=0;
  SelectTrack(aodEvent,nTracks,nSeleTrks,seleTrkFlags,mcArray);

  Bool_t  seleV0Flags[nV0s];
  Int_t     nSeleV0=0;
  SelectV0(aodEvent,nV0s,nSeleV0,seleV0Flags,mcArray);

  //------------------------------------------------
  // V0 loop 
  //------------------------------------------------
  for (Int_t iv0 = 0; iv0<nV0s; iv0++) {
    if(!seleV0Flags[iv0]) continue;
    AliAODv0 *v0 = aodEvent->GetV0(iv0);
    if(!v0) continue;

    AliAODTrack *cptrack =  (AliAODTrack*)(v0->GetDaughter(0));
    AliAODTrack *cntrack =  (AliAODTrack*)(v0->GetDaughter(1));

    //------------------------------------------------
    // track loop 
    //------------------------------------------------
    for (Int_t itrk = 0; itrk<nTracks; itrk++) {
      if(seleTrkFlags[itrk]!=1) continue;
      AliAODTrack *trk = (AliAODTrack*)aodEvent->GetTrack(itrk);
      //if(trk->GetID()<0) continue;
      
      //TPC only track (BIT 7) does not have PID information 
      //In addition to that, TPC only tracks does not have good DCA resolution
      //(according to femtoscopy code)
      AliAODTrack *trkpid = 0;
      if(fAnalCuts->GetProdAODFilterBit()==7){
        trkpid = fGTI[-trk->GetID()-1];
      }else{
        trkpid = trk;
      }

      Int_t cpid = cptrack->GetID();
      Int_t cnid = cntrack->GetID();
      Int_t lpid = trkpid->GetID();
      if((cpid==lpid)||(cnid==lpid)) continue;

      if(!fAnalCuts->SelectWithRoughCuts(v0,trk)) continue;

      AliAODVertex *secVert = ReconstructSecondaryVertex(v0,trk,aodEvent);
      if(!secVert) continue;

      AliAODRecoCascadeHF *lcobj = MakeCascadeHF(v0,trk,trkpid,aodEvent,secVert);
      if(!lcobj) {
				continue;
      }

      FillROOTObjects(lcobj,v0,trk,trkpid,aodEvent,mcArray);

      lcobj->GetSecondaryVtx()->RemoveDaughters();
      lcobj->UnsetOwnPrimaryVtx();
      delete lcobj;lcobj=NULL;
      delete secVert;
    }
  }

  if(fDoEventMixing){
		fEventInfo->SetString(Form("Ev%d_esd%d_E%d_V%d",AliAnalysisManager::GetAnalysisManager()->GetNcalls(),((AliAODHeader*)aodEvent->GetHeader())->GetEventNumberESDFile(),fProtonTracks->GetEntries(),fV0Tracks->GetEntries()));
    Int_t ind=GetPoolIndex(fVtxZ,fCentrality,fReactionPlane);
    if(ind>=0 && ind<fNOfPools){
      if(fEventBuffer[ind]->GetEntries() >= fNumberOfEventsForMixing){
				DoEventMixingWithPools(ind);
				if(fEventBuffer[ind]->GetEntries() >= 20*fNumberOfEventsForMixing){
					ResetPool(ind);
				}
      }
      fEventBuffer[ind]->Fill();
    }
  }

}

////-------------------------------------------------------------------------------
void AliAnalysisTaskSELc2pK0sfromAODtracks::FillROOTObjects(AliAODRecoCascadeHF *lcobj, AliAODv0 *v0, AliAODTrack *trk, AliAODTrack *trkpid, AliAODEvent *aodEvent, TClonesArray *mcArray) 
{
  //
  // Fill histograms or tree depending on fWriteVariableTree 
  //
	if(!trk) return;
	if(!trkpid) return;
	if(!v0) return;

  Double_t mprPDG =  TDatabasePDG::Instance()->GetParticle(2212)->Mass();
  Double_t mk0sPDG =  TDatabasePDG::Instance()->GetParticle(310)->Mass();
  Double_t mlamPDG =  TDatabasePDG::Instance()->GetParticle(3122)->Mass();

	for(Int_t i=0;i<54;i++){
		fCandidateVariables[i] = -9999.;
	}

	Double_t pxp = trk->Px();
	Double_t pyp = trk->Py();
	Double_t pzp = trk->Pz();
	Double_t momp = sqrt(pxp*pxp+pyp*pyp+pzp*pzp);
	Double_t Ep = sqrt(momp*momp+mprPDG*mprPDG);

	Double_t pxv = v0->Px();
	Double_t pyv = v0->Py();
	Double_t pzv = v0->Pz();
	Double_t momv = sqrt(pxv*pxv+pyv*pyv+pzv*pzv);
	Double_t mv = v0->MassK0Short();
	Double_t Ev = sqrt(momv*momv+mv*mv);

	Double_t cosoa = (pxp*pxv+pyp*pyv+pzp*pzv)/momp/momv;

  fCandidateVariables[ 0] = lcobj->InvMassLctoK0sP();
  fCandidateVariables[ 1] = lcobj->Px();
  fCandidateVariables[ 2] = lcobj->Py();
  fCandidateVariables[ 3] = lcobj->Pz();
  fCandidateVariables[ 4] = v0->MassK0Short();
  fCandidateVariables[ 5] = lcobj->PxProng(0);
  fCandidateVariables[ 6] = lcobj->PyProng(0);
  fCandidateVariables[ 7] = lcobj->PzProng(0);
  fCandidateVariables[ 8] = lcobj->PxProng(1);
  fCandidateVariables[ 9] = lcobj->PyProng(1);
  fCandidateVariables[10] = lcobj->PzProng(1);
  fCandidateVariables[11] = fVtx1->GetX();
  fCandidateVariables[12] = fVtx1->GetY();
  fCandidateVariables[13] = fVtx1->GetZ();
  fCandidateVariables[14] = fCentrality;
  fCandidateVariables[15] = lcobj->DecayLengthXY();
  fCandidateVariables[16] = (Float_t) fAnalCuts->CalculateLcCosPAXY(lcobj);

  Double_t nSigmaTPCpr=-9999.;
  Double_t nSigmaTOFpr=-9999.;
  Double_t probProton=-9999.;
	if(fAnalCuts->GetIsUsePID())
	{
		nSigmaTPCpr = fAnalCuts->GetPidHF()->GetPidResponse()->NumberOfSigmasTPC(trkpid,AliPID::kProton);
		nSigmaTOFpr = fAnalCuts->GetPidHF()->GetPidResponse()->NumberOfSigmasTOF(trkpid,AliPID::kProton);
		if(fAnalCuts->GetPidHF()->GetUseCombined()){
			probProton = fAnalCuts->GetProtonProbabilityTPCTOF(trk);
		}
		fCandidateVariables[17] = nSigmaTPCpr;
		fCandidateVariables[18] = nSigmaTOFpr;
		fCandidateVariables[19] = probProton;
	}
	fCandidateVariables[20] = 0;
	fCandidateVariables[21] = lcobj->Getd0Prong(0);
	fCandidateVariables[22] = lcobj->Getd0Prong(1);

	Double_t pospx = v0->MomPosX();
	Double_t pospy = v0->MomPosY();
	Double_t pospz = v0->MomPosZ();
	Double_t pose = sqrt(pospx*pospx+pospy*pospy+pospz*pospz+0.000511*0.000511);
	Double_t negpx = v0->MomNegX();
	Double_t negpy = v0->MomNegY();
	Double_t negpz = v0->MomNegZ();
	Double_t nege = sqrt(negpx*negpx+negpy*negpy+negpz*negpz+0.000511*0.000511);
	Double_t massPhoton = sqrt(pow(pose+nege,2)-pow(pospx+negpx,2)-pow(pospy+negpy,2)-pow(pospz+negpz,2));
	fCandidateVariables[23] = v0->MassLambda();
	fCandidateVariables[24] = v0->MassAntiLambda();
	fCandidateVariables[25] = massPhoton;
	fCandidateVariables[26] = v0->DcaPosToPrimVertex();
	fCandidateVariables[27] = v0->DcaNegToPrimVertex();
	fCandidateVariables[28] = v0->DcaV0ToPrimVertex();
	fCandidateVariables[29] = v0->DcaV0Daughters();

  AliAODMCParticle *mclc = 0;
  AliAODMCParticle *mcpr = 0;
  AliAODMCParticle *mcv0 = 0;
  Int_t mclablc = 0;
	Int_t mcpdgpr_array[100];
	Int_t mcpdgv0_array[100];
	Int_t mclabelpr_array[100];
	Int_t mclabelv0_array[100];
	Int_t mcngen_pr=-9999;
	Int_t mcngen_v0=-9999;
	Bool_t islambdac = kFALSE;
	Bool_t islambdac3body = kFALSE;
	Double_t decayvertx_mc = -9999.;
	Double_t decayverty_mc = -9999.;
	Double_t decayvertz_mc = -9999.;
	Double_t genvertx_mc = -9999.;
	Double_t genverty_mc = -9999.;
	Double_t genvertz_mc = -9999.;
	if(fUseMCInfo && mcArray){
		mclablc =  MatchToMC(lcobj,mcArray,mcpdgpr_array, mcpdgv0_array,mclabelpr_array,mclabelv0_array,mcngen_pr,mcngen_v0);
    if(mclablc>-1){
      mclc = (AliAODMCParticle*) mcArray->At(mclablc);
			if(mclabelpr_array[0]>=0)
				mcpr = (AliAODMCParticle*) mcArray->At(mclabelpr_array[0]);
			if(mclabelv0_array[0]>=0)
				mcv0 = (AliAODMCParticle*) mcArray->At(mclabelv0_array[0]);
		}
		if(mclc){
			Int_t pdgcode = mclc->GetPdgCode();
			if(abs(pdgcode)==4122 && abs(mcpdgpr_array[1])==4122 && abs(mcpdgv0_array[1])==311 && abs(mcpdgv0_array[2])==4122 && mclc->GetNDaughters()==2){
				islambdac = kTRUE;
			}
			if(abs(pdgcode)==4122 && abs(mcpdgpr_array[1])==2214 && abs(mcpdgpr_array[2])==4122 && abs(mcpdgv0_array[1])==311 && abs(mcpdgv0_array[2])==4122 && mclc->GetNDaughters()==2){
				islambdac3body = kTRUE;
			}
			if(abs(pdgcode)==4122 && abs(mcpdgpr_array[1])==3222 && abs(mcpdgpr_array[2])==4122 && abs(mcpdgv0_array[1])==311 && abs(mcpdgv0_array[2])==4122 && mclc->GetNDaughters()==2){
				islambdac3body = kTRUE;
			}
			if(abs(pdgcode)==4122 && abs(mcpdgpr_array[1])==3222 && abs(mcpdgpr_array[2])==3224 && abs(mcpdgpr_array[3])==4122 && abs(mcpdgv0_array[1])==311 && abs(mcpdgv0_array[2])==4122 && mclc->GetNDaughters()==2){
				islambdac3body = kTRUE;
			}
			if(abs(pdgcode)==4122 && abs(mcpdgpr_array[1])==4122 && abs(mcpdgv0_array[1])==311 && abs(mcpdgv0_array[2])==313 && abs(mcpdgv0_array[3])==4122 && mclc->GetNDaughters()==2){
				islambdac3body = kTRUE;
			}
			genvertx_mc = mclc->Xv();
			genverty_mc = mclc->Yv();
			genvertz_mc = mclc->Zv();
		}
		if(mcpr){
			decayvertx_mc = mcpr->Xv();
			decayverty_mc = mcpr->Yv();
			decayvertz_mc = mcpr->Zv();
		}
	}
	fCandidateVariables[30] = (Float_t) islambdac;

  Double_t LcPx = lcobj->Px();
  Double_t LcPy = lcobj->Py();
  Double_t LcPt = TMath::Sqrt(LcPx*LcPx+LcPy*LcPy);

  Double_t d0z0[2],covd0z0[3];
  trk->PropagateToDCA(fVtx1,fBzkG,kVeryBig,d0z0,covd0z0);
  Double_t x0 = fVtx1->GetX();
  Double_t y0 = fVtx1->GetY();
  Double_t px0 = LcPx/LcPt;
  Double_t py0 = LcPy/LcPt;
  Double_t tx[3];
  trk->GetXYZ(tx);
  Double_t x1 = tx[0];
  Double_t y1 = tx[1];
  trk->GetPxPyPz(tx);
  Double_t px1 = tx[0];
  Double_t py1 = tx[1];
  Double_t pt1 = sqrt(px1*px1+py1*py1);
  px1 = px1/pt1;
  py1 = py1/pt1;

  Double_t dx = x0 - x1;
  Double_t dy = y0 - y1;

  Double_t Delta = -px0*py1+py0*px1;
  Double_t a0 = -9999.;
	if(Delta!=0)
	{
		a0 = 1./Delta * (py1 * dx - px1 * dy);
	}
  Double_t neovertx = x0 + a0 * px0;
  Double_t neoverty = y0 + a0 * py0;
  Double_t z0 = fVtx1->GetZ();
  Double_t neovertz = z0 + TMath::Abs(a0)*trk->Pz()/trk->Pt();

	fCandidateVariables[31] = neovertx;
	fCandidateVariables[32] = neoverty;
	fCandidateVariables[33] = neovertz;
	fCandidateVariables[34] = (Float_t) islambdac3body;
	fCandidateVariables[35] = decayvertx_mc;
	fCandidateVariables[36] = decayverty_mc;
	fCandidateVariables[37] = decayvertz_mc;
	if(mclc){
		fCandidateVariables[38] = mclc->Px();
		fCandidateVariables[39] = mclc->Py();
		fCandidateVariables[40] = mclc->Pz();
	}
	fCandidateVariables[41] = genvertx_mc;
	fCandidateVariables[42] = genverty_mc;
	fCandidateVariables[43] = genvertz_mc;

  Double_t x1_k0s = v0->DecayVertexV0X();
  Double_t y1_k0s = v0->DecayVertexV0Y();
  Double_t px1_k0s = v0->Px();
  Double_t py1_k0s = v0->Py();
  Double_t pt1_k0s = sqrt(px1_k0s*px1_k0s+py1_k0s*py1_k0s);
  px1_k0s = px1_k0s/pt1_k0s;
  py1_k0s = py1_k0s/pt1_k0s;

  Double_t dx_k0s = x0 - x1_k0s;
  Double_t dy_k0s = y0 - y1_k0s;

  Double_t Delta_k0s = -px0*py1_k0s+py0*px1_k0s;
  Double_t a0_k0s = -9999.;
	if(Delta_k0s!=0)
	{
		a0_k0s = 1./Delta_k0s * (py1_k0s * dx_k0s - px1_k0s * dy_k0s);
	}
  Double_t neovertx_k0s = x0 + a0_k0s * px0;
  Double_t neoverty_k0s = y0 + a0_k0s * py0;
	fCandidateVariables[44] = neovertx_k0s;
	fCandidateVariables[45] = neoverty_k0s;
	fCandidateVariables[46] = lcobj->GetDCA();
  if(mclc){
    fCandidateVariables[47] = mclc->GetPdgCode();
    fCandidateVariables[48] = mcpdgpr_array[0];
    fCandidateVariables[49] = mcpdgpr_array[1];
    fCandidateVariables[50] = mcpdgpr_array[2];
    fCandidateVariables[51] = mcpdgv0_array[0];
    fCandidateVariables[52] = mcpdgv0_array[1];
    fCandidateVariables[53] = mcpdgv0_array[2];
  }


  if(fWriteVariableTree)
    fVariablesTree->Fill();

	Double_t cont_correlation[5];
	Double_t phi_trig = trk->Phi();
	Double_t phi_assoc = v0->Phi();
	Double_t eta_trig = trk->Eta();
	Double_t eta_assoc = v0->Eta();
	Double_t deltaphi = phi_assoc-phi_trig;
	if(deltaphi<-M_PI/2.) deltaphi += 2 * M_PI;
	if(deltaphi>3*M_PI/2.) deltaphi -= 2 * M_PI;
	cont_correlation[0] =	deltaphi;
	cont_correlation[1] = eta_assoc-eta_trig;
	cont_correlation[2] = trk->Pt();
	cont_correlation[3] = v0->Pt();
	cont_correlation[4] = fCentrality;
	fHistoK0spCorrelation->Fill(cont_correlation);

	Double_t cont[3];
	cont[0] = lcobj->InvMassLctoK0sP();
	cont[1] = lcobj->Pt();
	cont[2] = fCentrality;
	if(cosoa>0.) fHistoLcK0SpMassCoarse->Fill(cont);

	if(fAnalCuts->IsSelected(lcobj,AliRDHFCuts::kAll))
	{
	    fHistoLcK0SpMass->Fill(cont);

      fHistod0Bach->Fill(lcobj->Getd0Prong(0));
      fHistod0V0->Fill(lcobj->Getd0Prong(1));
      fHistod0d0->Fill(lcobj->Getd0Prong(0)*lcobj->Getd0Prong(1));
      fHistoV0CosPA->Fill(lcobj->CosV0PointingAngle());
      fHistoProbProton->Fill(probProton);
      fHistoDecayLength->Fill(lcobj->DecayLengthXY()*(fAnalCuts->CalculateLcCosPAXY(lcobj)));
      fHistoK0SMass->Fill(v0->MassK0Short());

			if(fUseMCInfo){
				if(islambdac){
					fHistoLcK0SpMassMCS->Fill(cont);
					fHistoK0spCorrelationMCS->Fill(cont_correlation);
          fHistoResponseLcPt->Fill(mclc->Pt(),lcobj->Pt());
				}
			}
	}
	if(fUseMCInfo){
		if(islambdac3body){
			fHistoLcK0SpPi0MassMCS->Fill(cont);
		}
	}

  return;
}
////-------------------------------------------------------------------------------
void AliAnalysisTaskSELc2pK0sfromAODtracks::FillMixROOTObjects(TLorentzVector *trkp, TLorentzVector *v0, TVector *prvars, TVector *v0vars) 
{
  ///
  /// Fill histograms or tree depending on fWriteVariableTree
  ///
	if(!trkp) return;
	if(!v0) return;

	for(Int_t i=0;i<54;i++){
		fCandidateVariables[i] = -9999.;
	}

	Double_t pxp = trkp->Px();
	Double_t pyp = trkp->Py();
	Double_t pzp = trkp->Pz();
	Double_t momp = sqrt(pxp*pxp+pyp*pyp+pzp*pzp);
	Double_t Ep = sqrt(momp*momp+0.938272*0.938272);

	Double_t pxv = v0->Px();
	Double_t pyv = v0->Py();
	Double_t pzv = v0->Pz();
	Double_t momv = sqrt(pxv*pxv+pyv*pyv+pzv*pzv);
	Double_t mv = v0->M();
	Double_t Ev = sqrt(momv*momv+mv*mv);

	Double_t cosoa = (pxp*pxv+pyp*pyv+pzp*pzv)/momp/momv;

	Double_t pxsum = pxp + pxv;
	Double_t pysum = pyp + pyv;
	Double_t pzsum = pzp + pzv;
	Double_t Esum = Ep + Ev;
	Double_t mpk0s = sqrt(Esum*Esum-pxsum*pxsum-pysum*pysum-pzsum*pzsum);


	Double_t cont_correlation[5];
	Double_t phi_trig = trkp->Phi();
	Double_t phi_assoc = v0->Phi();
	Double_t eta_trig = trkp->Eta();
	Double_t eta_assoc = v0->Eta();
	Double_t deltaphi = phi_assoc-phi_trig;
	if(deltaphi<-M_PI/2.) deltaphi += 2 * M_PI;
	if(deltaphi>3*M_PI/2.) deltaphi -= 2 * M_PI;
	cont_correlation[0] =	deltaphi;
	cont_correlation[1] = eta_assoc-eta_trig;
	cont_correlation[2] = trkp->Pt();
	cont_correlation[3] = v0->Pt();
	cont_correlation[4] = fCentrality;
	fHistoK0spCorrelationMix->Fill(cont_correlation);

	fCandidateVariables[ 0] = mpk0s;
	fCandidateVariables[ 1] = pxsum;
	fCandidateVariables[ 2] = pysum;
	fCandidateVariables[ 3] = pzsum;
	fCandidateVariables[ 4] = v0->M();
	fCandidateVariables[ 5] = trkp->Px();
	fCandidateVariables[ 6] = trkp->Py();
	fCandidateVariables[ 7] = trkp->Pz();
	fCandidateVariables[ 8] = v0->Px();
	fCandidateVariables[ 9] = v0->Py();
	fCandidateVariables[10] = v0->Pz();
	fCandidateVariables[20] = 1;

  Double_t LcPx = pxsum;
  Double_t LcPy = pysum;
  Double_t LcPt = TMath::Sqrt(LcPx*LcPx+LcPy*LcPy);

  Double_t x0 = (*prvars)[5];
  Double_t y0 = (*prvars)[6];
  Double_t px0 = LcPx/LcPt;
  Double_t py0 = LcPy/LcPt;
  Double_t tx[3];
  Double_t x1 = (*prvars)[1];
  Double_t y1 = (*prvars)[2];
  Double_t px1 = (*prvars)[3];
  Double_t py1 = (*prvars)[4];

  Double_t dx = x0 - x1;
  Double_t dy = y0 - y1;

  Double_t Delta = -px0*py1+py0*px1;
  Double_t a0 = -9999.;
	if(Delta!=0)
	{
		a0 = 1./Delta * (py1 * dx - px1 * dy);
	}
  Double_t neovertx = x0 + a0 * px0;
  Double_t neoverty = y0 + a0 * py0;

	fCandidateVariables[31] = neovertx;
	fCandidateVariables[32] = neoverty;

  if(fWriteVariableTree)
    fVariablesTree->Fill();

  Double_t rdhfcutvars[7];
  rdhfcutvars[0] = mpk0s;
  rdhfcutvars[1] = sqrt(pxsum*pxsum+pysum*pysum);
  rdhfcutvars[2] = (*prvars)[7];
  rdhfcutvars[3] = (*prvars)[8];
  rdhfcutvars[4] = a0;
  rdhfcutvars[5] = (*prvars)[0];
  rdhfcutvars[6] = (*v0vars)[0];

	if(fAnalCuts->IsSelected(trkp,v0,rdhfcutvars,AliRDHFCuts::kAll))
	{
    Double_t cont[3];
    cont[0] = mpk0s;
    cont[1] = sqrt(pxsum*pxsum+pysum*pysum);
    cont[2] = fCentrality;
    fHistoLcK0SpMassMixCoarse->Fill(cont);
		fHistoLcK0SpMassMix->Fill(cont);
	}

  return;
}

////-------------------------------------------------------------------------------
void AliAnalysisTaskSELc2pK0sfromAODtracks::DefineTreeVariables() 
{
  //
  // Define tree variables
  //

  const char* nameoutput = GetOutputSlot(4)->GetContainer()->GetName();
  fVariablesTree = new TTree(nameoutput,"Candidates variables tree");
  Int_t nVar = 54;
  fCandidateVariables = new Float_t [nVar];
  TString * fCandidateVariableNames = new TString[nVar];

  fCandidateVariableNames[ 0]="InvMassLc2pK0s";
  fCandidateVariableNames[ 1]="LcPx";
  fCandidateVariableNames[ 2]="LcPy";
  fCandidateVariableNames[ 3]="LcPz";
  fCandidateVariableNames[ 4]="massK0Short";
  fCandidateVariableNames[ 5]="BachPx";
  fCandidateVariableNames[ 6]="BachPy";
  fCandidateVariableNames[ 7]="BachPz";
  fCandidateVariableNames[ 8]="V0Px";
  fCandidateVariableNames[ 9]="V0Py";
  fCandidateVariableNames[10]="V0Pz";
  fCandidateVariableNames[11]="PrimVertx";
  fCandidateVariableNames[12]="PrimVerty";
  fCandidateVariableNames[13]="PrimVertz";
  fCandidateVariableNames[14]="Centrality";
  fCandidateVariableNames[15]="DecayLengthXY";
  fCandidateVariableNames[16]="LcCosPAXY";
  fCandidateVariableNames[17]="nSigmaTPCpr";
  fCandidateVariableNames[18]="nSigmaTOFpr";
  fCandidateVariableNames[19]="probProton";
  fCandidateVariableNames[20]="Mixing";
  fCandidateVariableNames[21]="Bachd0";
  fCandidateVariableNames[22]="V0d0";
  fCandidateVariableNames[23]="massLambda";
  fCandidateVariableNames[24]="massAntiLambda";
  fCandidateVariableNames[25]="massPhoton";
  fCandidateVariableNames[26]="DcaPosToPrimVtx";
  fCandidateVariableNames[27]="DcaNegToPrimVtx";
  fCandidateVariableNames[28]="DcaV0ToPrimVtx";
  fCandidateVariableNames[29]="DcaV0Daughters";
  fCandidateVariableNames[30]="IsLambdac";
  fCandidateVariableNames[31]="SecVertX";
  fCandidateVariableNames[32]="SecVertY";
  fCandidateVariableNames[33]="SecVertZ";
  fCandidateVariableNames[34]="islambdac3body";
  fCandidateVariableNames[35]="SecVertXMC";
  fCandidateVariableNames[36]="SecVertYMC";
  fCandidateVariableNames[37]="SecVertZMC";
  fCandidateVariableNames[38]="LcPxMC";
  fCandidateVariableNames[39]="LcPyMC";
  fCandidateVariableNames[40]="LcPzMC";
  fCandidateVariableNames[41]="PrimVertXMC";
  fCandidateVariableNames[42]="PrimVertYMC";
  fCandidateVariableNames[43]="PrimVertZMC";
  fCandidateVariableNames[44]="SecVertXK0s";
  fCandidateVariableNames[45]="SecVertYK0s";
  fCandidateVariableNames[46]="DcaBachV0";
  fCandidateVariableNames[47]="MatchedPDG";
  fCandidateVariableNames[48]="ProtonPDG";
  fCandidateVariableNames[49]="MotherProtonPDG";
  fCandidateVariableNames[50]="GrMotherProtonPDG";
  fCandidateVariableNames[51]="V0PDG";
  fCandidateVariableNames[52]="MotherV0PDG";
  fCandidateVariableNames[53]="GrMotherV0PDG";

  for (Int_t ivar=0; ivar<nVar; ivar++) {
    fVariablesTree->Branch(fCandidateVariableNames[ivar].Data(),&fCandidateVariables[ivar],Form("%s/f",fCandidateVariableNames[ivar].Data()));
  }

  return;
}

////__________________________________________________________________________
void  AliAnalysisTaskSELc2pK0sfromAODtracks::DefineGeneralHistograms() {
  //
  // This is to define general histograms
  //

  fCEvents = new TH1F("fCEvents","conter",18,-0.5,17.5);
  fCEvents->SetStats(kTRUE);
  fCEvents->GetXaxis()->SetBinLabel(1,"X1");
  fCEvents->GetXaxis()->SetBinLabel(2,"Analyzed events");
  fCEvents->GetXaxis()->SetBinLabel(3,"AliAODVertex exists");
  fCEvents->GetXaxis()->SetBinLabel(4,"TriggerOK");
  fCEvents->GetXaxis()->SetBinLabel(5,"IsEventSelected");
  fCEvents->GetXaxis()->SetBinLabel(6,"CascadesHF exists");
  fCEvents->GetXaxis()->SetBinLabel(7,"MCarray exists");
  fCEvents->GetXaxis()->SetBinLabel(8,"MCheader exists");
  fCEvents->GetXaxis()->SetBinLabel(9,"triggerClass!=CINT1");
  fCEvents->GetXaxis()->SetBinLabel(10,"triggerMask!=kAnyINT");
  fCEvents->GetXaxis()->SetBinLabel(11,"triggerMask!=kAny");
  fCEvents->GetXaxis()->SetBinLabel(12,"vtxTitle.Contains(Z)");
  fCEvents->GetXaxis()->SetBinLabel(13,"vtxTitle.Contains(3D)");
  fCEvents->GetXaxis()->SetBinLabel(14,"vtxTitle.Doesn'tContain(Z-3D)");
  fCEvents->GetXaxis()->SetBinLabel(15,Form("zVtx<=%2.0fcm",fAnalCuts->GetMaxVtxZ()));
  fCEvents->GetXaxis()->SetBinLabel(16,"!IsEventSelected");
  fCEvents->GetXaxis()->SetBinLabel(17,"triggerMask!=kAnyINT || triggerClass!=CINT1");
  fCEvents->GetXaxis()->SetBinLabel(18,Form("zVtxMC<=%2.0fcm",fAnalCuts->GetMaxVtxZ()));
  //fCEvents->GetXaxis()->SetTitle("");
  fCEvents->GetYaxis()->SetTitle("counts");

  fHTrigger = new TH1F("fHTrigger","counter",18,-0.5,17.5);
  fHTrigger->SetStats(kTRUE);
  fHTrigger->GetXaxis()->SetBinLabel(1,"X1");
  fHTrigger->GetXaxis()->SetBinLabel(2,"kMB");
  fHTrigger->GetXaxis()->SetBinLabel(3,"kSemiCentral");
  fHTrigger->GetXaxis()->SetBinLabel(4,"kCentral");
  fHTrigger->GetXaxis()->SetBinLabel(5,"kINT7");
  fHTrigger->GetXaxis()->SetBinLabel(6,"kEMC7");
  //fHTrigger->GetXaxis()->SetBinLabel(7,"Space");
  fHTrigger->GetXaxis()->SetBinLabel(8,"kMB|kSemiCentral|kCentral");
  fHTrigger->GetXaxis()->SetBinLabel(9,"kINT7|kEMC7");
  fHTrigger->GetXaxis()->SetBinLabel(11,"kMB&kSemiCentral");
  fHTrigger->GetXaxis()->SetBinLabel(12,"kMB&kCentral");
  fHTrigger->GetXaxis()->SetBinLabel(13,"kINT7&kEMC7");

  fHCentrality = new TH1F("fHCentrality","conter",100,0.,100.);
  fHReactionPlane = new TH1F("fHReactionPlane","conter",100,-1.6,3.2);

  fOutput->Add(fCEvents);
  fOutput->Add(fHTrigger);
  fOutput->Add(fHCentrality);
  fOutput->Add(fHReactionPlane);

  return;
}
//__________________________________________________________________________
void  AliAnalysisTaskSELc2pK0sfromAODtracks::DefineAnalysisHistograms() 
{
  //
  // Define analyis histograms
  //
	
  //------------------------------------------------
  // Basic histogram
  //------------------------------------------------
  Int_t bins_base[3]=		{200			,20		,10};
  Double_t xmin_base[3]={2.286-0.5,0		,0.00};
  Double_t xmax_base[3]={2.286+0.5,20.	,100};
  fHistoLcK0SpMass = new THnSparseF("fHistoLcK0SpMass","",3,bins_base,xmin_base,xmax_base);
  fOutputAll->Add(fHistoLcK0SpMass);
  fHistoLcK0SpMassMix = new THnSparseF("fHistoLcK0SpMassMix","",3,bins_base,xmin_base,xmax_base);
  fOutputAll->Add(fHistoLcK0SpMassMix);
  fHistoLcK0SpMassMCS = new THnSparseF("fHistoLcK0SpMassMCS","",3,bins_base,xmin_base,xmax_base);
  fOutputAll->Add(fHistoLcK0SpMassMCS);
  fHistoLcK0SpPi0MassMCS = new THnSparseF("fHistoLcK0SpPi0MassMCS","",3,bins_base,xmin_base,xmax_base);
  fOutputAll->Add(fHistoLcK0SpPi0MassMCS);
  fHistoLcKPluspMass = new THnSparseF("fHistoLcKPluspMass","",3,bins_base,xmin_base,xmax_base);
  fOutputAll->Add(fHistoLcKPluspMass);
  fHistoLcKMinuspMass = new THnSparseF("fHistoLcKMinuspMass","",3,bins_base,xmin_base,xmax_base);
  fOutputAll->Add(fHistoLcKMinuspMass);
  fHistoLcKPluspMassMix = new THnSparseF("fHistoLcKPluspMassMix","",3,bins_base,xmin_base,xmax_base);
  fOutputAll->Add(fHistoLcKPluspMassMix);
  fHistoLcKMinuspMassMix = new THnSparseF("fHistoLcKMinuspMassMix","",3,bins_base,xmin_base,xmax_base);
  fOutputAll->Add(fHistoLcKMinuspMassMix);

  Int_t bins_lcmcgen[3]=	{100 ,20	,10};
  Double_t xmin_lcmcgen[3]={0.,-1.0	,0.0};
  Double_t xmax_lcmcgen[3]={20.,1.0	,100};
  fHistoLcMCGen = new THnSparseF("fHistoLcMCGen","",3,bins_lcmcgen,xmin_lcmcgen,xmax_lcmcgen);
  fOutputAll->Add(fHistoLcMCGen);

  Int_t bins_coarse[3]=		{160	,20		,10};
  Double_t xmin_coarse[3]={1.,0		,0.00};
  Double_t xmax_coarse[3]={5.8,20.	,100};
  fHistoLcK0SpMassCoarse = new THnSparseF("fHistoLcK0SpMassCoarse","",3,bins_coarse,xmin_coarse,xmax_coarse);
  fOutputAll->Add(fHistoLcK0SpMassCoarse);
  fHistoLcK0SpMassMixCoarse = new THnSparseF("fHistoLcK0SpMassMixCoarse","",3,bins_coarse,xmin_coarse,xmax_coarse);
  fOutputAll->Add(fHistoLcK0SpMassMixCoarse);

  Int_t bins_correlation[5]=	{40	,20,20,20,10};
  Double_t xmin_correlation[5]={-M_PI/2.,-2.4,0.,0.,0.0};
  Double_t xmax_correlation[5]={1.5*M_PI,2.4,10.,10.,100};
  fHistoK0spCorrelation = new THnSparseF("fHistoK0spCorrelation","",5,bins_correlation,xmin_correlation,xmax_correlation);
  fOutputAll->Add(fHistoK0spCorrelation);
  fHistoK0spCorrelationMix = new THnSparseF("fHistoK0spCorrelationMix","",5,bins_correlation,xmin_correlation,xmax_correlation);
  fOutputAll->Add(fHistoK0spCorrelationMix);
  fHistoK0spCorrelationMCS = new THnSparseF("fHistoK0spCorrelationMCS","",5,bins_correlation,xmin_correlation,xmax_correlation);
  fOutputAll->Add(fHistoK0spCorrelationMCS);


  //------------------------------------------------
  // checking histograms
  //------------------------------------------------
  fHistoBachPt = new TH2F("fHistoBachPt","Bachelor p_{T}",100,0.0,10.0,20,-1.,1.);
  fOutputAll->Add(fHistoBachPt);
  fHistoBachPtMCS = new TH2F("fHistoBachPtMCS","Bachelor p_{T}",100,0.0,10.0,20,-1.,1.);
  fOutputAll->Add(fHistoBachPtMCS);
  fHistoBachPtMCGen = new TH2F("fHistoBachPtMCGen","Bachelor p_{T}",100,0.0,10.0,20,-1.,1.);
  fOutputAll->Add(fHistoBachPtMCGen);

  fHistoKaonPt = new TH2F("fHistoKaonPt","Kaon p_{T}",100,0.0,10.0,20,-1.,1.);
  fOutputAll->Add(fHistoKaonPt);
  fHistoKaonPtMCS = new TH2F("fHistoKaonPtMCS","Kaon p_{T}",100,0.0,10.0,20,-1.,1.);
  fOutputAll->Add(fHistoKaonPtMCS);
  fHistoKaonPtMCGen = new TH2F("fHistoKaonPtMCGen","Kaon p_{T}",100,0.0,10.0,20,-1.,1.);
  fOutputAll->Add(fHistoKaonPtMCGen);

  fHistoK0sMassvsPt=new TH3F("fHistoK0sMassvsPt","K0s mass",100,0.497-0.05,0.497+0.05,20,0.,10.,20,-1.,1.);
  fOutputAll->Add(fHistoK0sMassvsPt);
  fHistoK0sMassvsPtMCS=new TH3F("fHistoK0sMassvsPtMCS","K0s mass",100,0.497-0.05,0.497+0.05,20,0.,10.,20,-1.,1.);
  fOutputAll->Add(fHistoK0sMassvsPtMCS);
  fHistoK0sMassvsPtMCGen=new TH3F("fHistoK0sMassvsPtMCGen","K0s mass",100,0.497-0.05,0.497+0.05,20,0.,10.,20,-1.,1.);
  fOutputAll->Add(fHistoK0sMassvsPtMCGen);

  fHistod0Bach = new TH1F("fHistod0Bach","Bachelor d_{0}",100,-0.5,0.5);
  fOutputAll->Add(fHistod0Bach);
  fHistod0V0 = new TH1F("fHistod0V0","V_{0} d_{0}",100,-0.5,0.5);
  fOutputAll->Add(fHistod0V0);
  fHistod0d0 = new TH1F("fHistod0d0","d_{0} d_{0}",100,-0.5,0.5);
  fOutputAll->Add(fHistod0d0);
  fHistoV0CosPA=new TH1F("fHistoV0CosPA","V0->Second vertex cospa",100,-1.,1.0);
  fOutputAll->Add(fHistoV0CosPA);
  fHistoProbProton=new TH1F("fHistoProbProton","ProbProton",100,0.,1.0);
  fOutputAll->Add(fHistoProbProton);
  fHistoDecayLength=new TH1F("fHistoDecayLength","Decay Length",100,-0.1,0.1);
  fOutputAll->Add(fHistoDecayLength);
  fHistoK0SMass=new TH1F("fHistoK0SMass","K0S mass",100,0.497-0.05,0.497+0.05);
  fOutputAll->Add(fHistoK0SMass);

  fHistonEvtvsRunNumber=new TH1F("fHistonEvtvsRunNumber","",20000,-0.5,19999.5);
  fOutputAll->Add(fHistonEvtvsRunNumber);
  fHistonProtonvsRunNumber=new TH1F("fHistonProtonvsRunNumber","",20000,-0.5,19999.5);
  fOutputAll->Add(fHistonProtonvsRunNumber);
  fHistonK0svsRunNumber=new TH1F("fHistonK0svsRunNumber","",20000,-0.5,19999.5);
  fOutputAll->Add(fHistonK0svsRunNumber);

  fHistoMassTagV0Min=new TH1F("fHistoMassTagV0Min","",1500,0,1.5);
  fOutputAll->Add(fHistoMassTagV0Min);
  fHistoMassTagV0SameSignMin=new TH1F("fHistoMassTagV0SameSignMin","",1500,0,1.5);
  fOutputAll->Add(fHistoMassTagV0SameSignMin);

  fHistoResponseLcPt = new TH2D("fHistoResponseLcPt","",100,0.,20.,100,0.,20.);
  fOutputAll->Add(fHistoResponseLcPt);
  fHistoResponseLcPt1 = new TH2D("fHistoResponseLcPt1","",100,0.,20.,100,0.,20.);
  fOutputAll->Add(fHistoResponseLcPt1);
  fHistoResponseLcPt2 = new TH2D("fHistoResponseLcPt2","",100,0.,20.,100,0.,20.);
  fOutputAll->Add(fHistoResponseLcPt2);

  return;
}

//________________________________________________________________________
AliAODRecoCascadeHF* AliAnalysisTaskSELc2pK0sfromAODtracks::MakeCascadeHF(AliAODv0 *v0, AliAODTrack *part, AliAODTrack *partpid, AliAODEvent * aod, AliAODVertex *secVert) 
{
  //
  // Create AliAODRecoCascadeHF object from the argument
  //

  if(!v0) return 0x0;
  if(!part) return 0x0;
  if(!partpid) return 0x0;
  if(!aod) return 0x0;

  //------------------------------------------------
  // PrimaryVertex
  //------------------------------------------------
  AliAODVertex *primVertexAOD;
  Bool_t unsetvtx = kFALSE;
  if(fReconstructPrimVert){
    primVertexAOD = CallPrimaryVertex(v0,part,aod);
    if(!primVertexAOD){
      primVertexAOD = fVtx1;
    }else{
      unsetvtx = kTRUE;
    }
  }else{
    primVertexAOD = fVtx1;
  }
  if(!primVertexAOD) return 0x0;
  Double_t posprim[3]; primVertexAOD->GetXYZ(posprim);

  //------------------------------------------------
  // DCA between tracks
  //------------------------------------------------
  AliESDtrack *esdtrack = new AliESDtrack((AliVTrack*)partpid);

  AliNeutralTrackParam *trackV0=NULL;
  const AliVTrack *trackVV0 = dynamic_cast<const AliVTrack*>(v0);
  if(trackVV0)  trackV0 = new AliNeutralTrackParam(trackVV0);

  Double_t xdummy, ydummy;
  Double_t dca = esdtrack->GetDCA(trackV0,fBzkG,xdummy,ydummy);


  //------------------------------------------------
  // Propagate all tracks to the secondary vertex and calculate momentum there
  //------------------------------------------------
	
  Double_t d0z0bach[2],covd0z0bach[3];
//  if(sqrt(pow(secVert->GetX(),2)+pow(secVert->GetY(),2))<1.){
//    part->PropagateToDCA(secVert,fBzkG,kVeryBig,d0z0bach,covd0z0bach);
//    trackV0->PropagateToDCA(secVert,fBzkG,kVeryBig);
//  }else{
//    part->PropagateToDCA(primVertexAOD,fBzkG,kVeryBig,d0z0bach,covd0z0bach);
//    trackV0->PropagateToDCA(primVertexAOD,fBzkG,kVeryBig);
//  }
  part->PropagateToDCA(primVertexAOD,fBzkG,kVeryBig,d0z0bach,covd0z0bach);
  trackV0->PropagateToDCA(primVertexAOD,fBzkG,kVeryBig);
  Double_t momv0_new[3]={-9999,-9999,-9999.};
  trackV0->GetPxPyPz(momv0_new);

  Double_t px[2],py[2],pz[2];
  px[0] = part->Px(); py[0] = part->Py(); pz[0] = part->Pz(); 
  px[1] = momv0_new[0]; py[1] = momv0_new[1]; pz[1] = momv0_new[2]; 

  //------------------------------------------------
  // d0
  //------------------------------------------------
  Double_t d0[3],d0err[3];

  part->PropagateToDCA(primVertexAOD,fBzkG,kVeryBig,d0z0bach,covd0z0bach);
  d0[0]= d0z0bach[0];
  d0err[0] = TMath::Sqrt(covd0z0bach[0]);

  Double_t d0z0v0[2],covd0z0v0[3];
  trackV0->PropagateToDCA(primVertexAOD,fBzkG,kVeryBig,d0z0v0,covd0z0v0);
  d0[1]= d0z0v0[0];
  d0err[1] = TMath::Sqrt(covd0z0v0[0]);

  //------------------------------------------------
  // Create AliAODRecoCascadeHF
  //------------------------------------------------
  Short_t charge = part->Charge();
  AliAODRecoCascadeHF *theCascade = new AliAODRecoCascadeHF(secVert,charge,px,py,pz,d0,d0err,dca);
  if(!theCascade)  
    {
      if(unsetvtx) delete primVertexAOD; primVertexAOD=NULL;
      if(esdtrack) delete esdtrack;
      if(trackV0) delete trackV0;
      return 0x0;
    }
  theCascade->SetOwnPrimaryVtx(primVertexAOD);
  UShort_t id[2]={(UShort_t)part->GetID(),(UShort_t)trackV0->GetID()};
  theCascade->SetProngIDs(2,id);

  theCascade->GetSecondaryVtx()->AddDaughter(partpid);
  theCascade->GetSecondaryVtx()->AddDaughter(v0);

  if(unsetvtx) delete primVertexAOD; primVertexAOD=NULL;
  if(esdtrack) delete esdtrack;
  if(trackV0) delete trackV0;

  return theCascade;
}

//________________________________________________________________________
AliAODVertex* AliAnalysisTaskSELc2pK0sfromAODtracks::CallPrimaryVertex(AliAODv0 *v0, AliAODTrack *trk, AliAODEvent* aod)
{
  //
  // Make an array of tracks which should not be used in primary vertex calculation and 
  // Call PrimaryVertex function
  //

  TObjArray *TrackArray = new TObjArray(3);
  
  AliESDtrack *cptrk1 = new AliESDtrack((AliVTrack*)trk);
  TrackArray->AddAt(cptrk1,0);
  
  AliESDtrack *cascptrack = new AliESDtrack((AliVTrack*)v0->GetDaughter(0));
  TrackArray->AddAt(cascptrack,1);
  AliESDtrack *cascntrack = new AliESDtrack((AliVTrack*)v0->GetDaughter(1));
  TrackArray->AddAt(cascntrack,2);
  
  AliAODVertex *newvert  = PrimaryVertex(TrackArray,aod);
  
  for(Int_t i=0;i<3;i++)
    {
      AliESDtrack *tesd = (AliESDtrack*)TrackArray->UncheckedAt(i);
      delete tesd;
    }
  TrackArray->Clear();
  delete TrackArray;
  
  return newvert;
}

//________________________________________________________________________
AliAODVertex* AliAnalysisTaskSELc2pK0sfromAODtracks::PrimaryVertex(const TObjArray *trkArray,
								   AliVEvent *event)
{
  //
  //Used only for pp
  //copied from AliAnalysisVertexingHF (except for the following 3 lines)
  //

  Bool_t fRecoPrimVtxSkippingTrks = kTRUE;
  Bool_t fRmTrksFromPrimVtx = kFALSE;

  AliESDVertex *vertexESD = 0;
  AliAODVertex *vertexAOD = 0;
  
  //vertexESD = new AliESDVertex(*fV1);
  

  if(!fRecoPrimVtxSkippingTrks && !fRmTrksFromPrimVtx) { 
    // primary vertex from the input event
    
    vertexESD = new AliESDVertex(*fV1);
    
  } else {
    // primary vertex specific to this candidate
    
    Int_t nTrks = trkArray->GetEntriesFast();
    AliVertexerTracks *vertexer = new AliVertexerTracks(event->GetMagneticField());
    
    if(fRecoPrimVtxSkippingTrks) { 
      // recalculating the vertex
      
      if(strstr(fV1->GetTitle(),"VertexerTracksWithConstraint")) {
	Float_t diamondcovxy[3];
	event->GetDiamondCovXY(diamondcovxy);
	Double_t pos[3]={event->GetDiamondX(),event->GetDiamondY(),0.};
	Double_t cov[6]={diamondcovxy[0],diamondcovxy[1],diamondcovxy[2],0.,0.,10.*10.};
	AliESDVertex *diamond = new AliESDVertex(pos,cov,1.,1);
	vertexer->SetVtxStart(diamond);
	delete diamond; diamond=NULL;
	if(strstr(fV1->GetTitle(),"VertexerTracksWithConstraintOnlyFitter")) 
	  vertexer->SetOnlyFitter();
      }
      Int_t skipped[1000];
      Int_t nTrksToSkip=0,id;
      AliExternalTrackParam *t = 0;
      for(Int_t i=0; i<nTrks; i++) {
	t = (AliExternalTrackParam*)trkArray->UncheckedAt(i);
	id = (Int_t)t->GetID();
	if(id<0) continue;
	skipped[nTrksToSkip++] = id;
      }
      // TEMPORARY FIX
      // For AOD, skip also tracks without covariance matrix
      Double_t covtest[21];
      for(Int_t j=0; j<event->GetNumberOfTracks(); j++) {
	AliVTrack *vtrack = (AliVTrack*)event->GetTrack(j);
	if(!vtrack->GetCovarianceXYZPxPyPz(covtest)) {
	  id = (Int_t)vtrack->GetID();
	  if(id<0) continue;
	  skipped[nTrksToSkip++] = id;
	}
      }
      for(Int_t ijk=nTrksToSkip; ijk<1000; ijk++) skipped[ijk]=-1;
      //
      vertexer->SetSkipTracks(nTrksToSkip,skipped);
      vertexESD = (AliESDVertex*)vertexer->FindPrimaryVertex(event); 
      
    } else if(fRmTrksFromPrimVtx && nTrks>0) { 
      // removing the prongs tracks
      
      TObjArray rmArray(nTrks);
      UShort_t *rmId = new UShort_t[nTrks];
      AliESDtrack *esdTrack = 0;
      AliESDtrack *t = 0;
      for(Int_t i=0; i<nTrks; i++) {
	t = (AliESDtrack*)trkArray->UncheckedAt(i);
	esdTrack = new AliESDtrack(*t);
	rmArray.AddLast(esdTrack);
	if(esdTrack->GetID()>=0) {
	  rmId[i]=(UShort_t)esdTrack->GetID();
	} else {
	  rmId[i]=9999;
	}
      }
      Float_t diamondxy[2]={static_cast<Float_t>(event->GetDiamondX()),static_cast<Float_t>(event->GetDiamondY())};
      vertexESD = vertexer->RemoveTracksFromVertex(fV1,&rmArray,rmId,diamondxy);
      delete [] rmId; rmId=NULL;
      rmArray.Delete();
      
    }
    
    delete vertexer; vertexer=NULL;
    if(!vertexESD) return vertexAOD;
    if(vertexESD->GetNContributors()<=0) { 
      //AliDebug(2,"vertexing failed"); 
      delete vertexESD; vertexESD=NULL;
      return vertexAOD;
    }
    
    
  }
  
  // convert to AliAODVertex
  Double_t pos[3],cov[6],chi2perNDF;
  vertexESD->GetXYZ(pos); // position
  vertexESD->GetCovMatrix(cov); //covariance matrix
  chi2perNDF = vertexESD->GetChi2toNDF();
  delete vertexESD; vertexESD=NULL;
  
  vertexAOD = new AliAODVertex(pos,cov,chi2perNDF);
  
  return vertexAOD;
}

//________________________________________________________________________
AliAODVertex* AliAnalysisTaskSELc2pK0sfromAODtracks::ReconstructSecondaryVertex(AliAODv0 *v0, AliAODTrack *part, AliAODEvent * aod) 
{
  //
  // Reconstruct secondary vertex from trkArray (Copied from AliAnalysisVertexingHF)
  //
  //------------------------------------------------
  // PrimaryVertex
  //------------------------------------------------
  AliAODVertex *primVertexAOD;
  Bool_t unsetvtx = kFALSE;
  if(fReconstructPrimVert){
    primVertexAOD = CallPrimaryVertex(v0,part,aod);
    if(!primVertexAOD){
      primVertexAOD = fVtx1;
    }else{
      unsetvtx = kTRUE;
    }
  }else{
    primVertexAOD = fVtx1;
  }
  if(!primVertexAOD) return 0x0;

  //------------------------------------------------
  // Secondary vertex
  //------------------------------------------------

  Double_t LcPx = part->Px()+v0->Px();
  Double_t LcPy = part->Py()+v0->Py();
  Double_t LcPt = TMath::Sqrt(LcPx*LcPx+LcPy*LcPy);

  Double_t d0z0[2],covd0z0[3];
  part->PropagateToDCA(primVertexAOD,fBzkG,kVeryBig,d0z0,covd0z0);
  Double_t x0 = primVertexAOD->GetX();
  Double_t y0 = primVertexAOD->GetY();
  Double_t px0 = LcPx/LcPt;
  Double_t py0 = LcPy/LcPt;
  Double_t tx[3];
  part->GetXYZ(tx);
  Double_t x1 = tx[0];
  Double_t y1 = tx[1];
  part->GetPxPyPz(tx);
  Double_t px1 = tx[0];
  Double_t py1 = tx[1];
  Double_t pt1 = sqrt(px1*px1+py1*py1);
  px1 = px1/pt1;
  py1 = py1/pt1;

  Double_t dx = x0 - x1;
  Double_t dy = y0 - y1;

  Double_t Delta = -px0*py1+py0*px1;
  Double_t a0 = -9999.;
  if(Delta!=0)
    {
      a0 = 1./Delta * (py1 * dx - px1 * dy);
    }
  Double_t neovertx = x0 + a0 * px0;
  Double_t neoverty = y0 + a0 * py0;
  Double_t z0 = primVertexAOD->GetZ();
  Double_t neovertz = z0 + TMath::Abs(a0)*part->Pz()/part->Pt();

  if(unsetvtx) delete primVertexAOD; primVertexAOD=NULL;

  Double_t pos[3],cov[6],chi2perNDF;
  pos[0]=neovertx;
  pos[1]=neoverty;
  pos[2]=neovertz;
  cov[0]=0.0;
  cov[1]=0.0;
  cov[2]=0.0;
  cov[3]=0.0;
  cov[4]=0.0;
  cov[5]=0.0;
  chi2perNDF=0.0;
  AliAODVertex *secVert = new AliAODVertex(pos,cov,chi2perNDF);
  if(!secVert){
    return 0x0;
  }
  return secVert;
}
//________________________________________________________________________
//AliAODVertex* AliAnalysisTaskSELc2pK0sfromAODtracks::ReconstructSecondaryVertex(AliAODv0 *v0, AliAODTrack *part, AliAODEvent * aod) 
//{
//  //
//  // Reconstruct secondary vertex from trkArray (Copied from AliAnalysisVertexingHF)
//	// Currently only returns Primary vertex (can we reconstruct secondary vertex from p - v0)
//  //
//	
//  AliAODVertex *primVertexAOD;
//  Bool_t unsetvtx = kFALSE;
//  if(fReconstructPrimVert){
//    primVertexAOD = CallPrimaryVertex(v0,part,aod);
//    if(!primVertexAOD){
//      primVertexAOD = fVtx1;
//    }else{
//      unsetvtx = kTRUE;
//    }
//  }else{
//    primVertexAOD = fVtx1;
//  }
//  if(!primVertexAOD) return 0x0;
//
//  AliESDVertex * vertexESD = new AliESDVertex(*fV1);
//
//  Double_t pos[3],cov[6],chi2perNDF;
//  vertexESD->GetXYZ(pos); // position
//  vertexESD->GetCovMatrix(cov); //covariance matrix
//  chi2perNDF = vertexESD->GetChi2toNDF();
//  delete vertexESD; vertexESD=NULL;
//  
//  AliAODVertex *secVert = new AliAODVertex(pos,cov,chi2perNDF);
//
//  return secVert;
//}
//________________________________________________________________________
void AliAnalysisTaskSELc2pK0sfromAODtracks::SelectTrack( const AliVEvent *event, Int_t trkEntries, Int_t &nSeleTrks,Int_t *seleFlags, TClonesArray *mcArray)
{
  //
  // Select good tracks using fAnalCuts (AliRDHFCuts object) and return the array of their ids
  //
  
  if(trkEntries==0) return;
  
  nSeleTrks=0;
  for(Int_t i=0; i<trkEntries; i++) {
    seleFlags[i] =0;
    
    AliVTrack *track;
    track = (AliVTrack*)event->GetTrack(i);
    
    //if(track->GetID()<0) continue;
    Double_t covtest[21];
    if(!track->GetCovarianceXYZPxPyPz(covtest)) continue;
    if(!fAnalCuts) continue;
    
    AliAODTrack *aodt = (AliAODTrack*)track;

    if(fAnalCuts->GetProdUseAODFilterBit()){
      Int_t filterbit = fAnalCuts->GetProdAODFilterBit();
      if(filterbit==7){
        if(!aodt->TestFilterBit(BIT(filterbit))) continue;
      }else{
        if(!aodt->TestFilterMask(BIT(filterbit))) continue;
      }
    }

    AliAODTrack *aodtpid = 0;
    if(fAnalCuts->GetProdAODFilterBit()==7){
      aodtpid = fGTI[-aodt->GetID()-1];
    }else{
      aodtpid = aodt;
    }

    if(fAnalCuts->SingleTrkCuts(aodt,aodtpid,fVtx1)){
      seleFlags[i]=1;
      nSeleTrks++;
			FillProtonROOTObjects(aodt,mcArray);

			Double_t minmass;
			Bool_t isv0 = fAnalCuts->TagV0(aodt,(AliAODEvent*)event,trkEntries,minmass);
			fHistoMassTagV0Min->Fill(minmass);
			if(isv0) seleFlags[i] = 0;

//			Double_t minmasslike;
//			fAnalCuts->TagV0SameSign(aodt,(AliAODEvent*)event,trkEntries,minmasslike);
//			fHistoMassTagV0SameSignMin->Fill(minmasslike);
    }
  } // end loop on tracks
}
//________________________________________________________________________
void AliAnalysisTaskSELc2pK0sfromAODtracks::SelectV0( const AliVEvent *event,Int_t nV0s,Int_t &nSeleV0, Bool_t *seleV0Flags, TClonesArray *mcArray)
{
  //
  // Select good V0 using fAnalCuts (AliRDHFCuts object) and return the array of their ids
  //

  nSeleV0 = 0;
  for(Int_t iv0=0;iv0<nV0s;iv0++)
    {
      seleV0Flags[iv0] = kFALSE;
      AliAODv0 *v0 = ((AliAODEvent*)event)->GetV0(iv0);

      if(!fAnalCuts) continue;
      if(fAnalCuts->SingleV0Cuts(v0,fVtx1)){
				seleV0Flags[iv0] = kTRUE;
				nSeleV0++;

				FillV0ROOTObjects(v0, mcArray);
      }
    }
}

////-------------------------------------------------------------------------------
void AliAnalysisTaskSELc2pK0sfromAODtracks::DefineProtonTreeVariables() 
{
  //
  /// Define proton tree variables
  //

  const char* nameoutput = GetOutputSlot(5)->GetContainer()->GetName();
  fProtonVariablesTree = new TTree(nameoutput,"proton variables tree");
  Int_t nVar = 3;
  fCandidateProtonVariables = new Float_t [nVar];
  TString * fCandidateVariableNames = new TString[nVar];

  fCandidateVariableNames[ 0]="PrPx";
  fCandidateVariableNames[ 1]="PrPy";
  fCandidateVariableNames[ 2]="PrPz";

  for (Int_t ivar=0; ivar<nVar; ivar++) {
    fProtonVariablesTree->Branch(fCandidateVariableNames[ivar].Data(),&fCandidateProtonVariables[ivar],Form("%s/f",fCandidateVariableNames[ivar].Data()));
  }

  return;
}
//________________________________________________________________________
void AliAnalysisTaskSELc2pK0sfromAODtracks::FillProtonROOTObjects(AliAODTrack *trk, TClonesArray *mcArray) 
{
  //
  /// Fill histograms or tree depending on fWriteVariableTree
  //

	if(!trk) return;

	TLorentzVector vpr;
	vpr.SetXYZM(trk->Px(),trk->Py(),trk->Pz(),0.938272);
  fHistoBachPt->Fill(vpr.Pt(),vpr.Rapidity());
	if(fDoEventMixing){
		fProtonTracks->AddLast(new TLorentzVector(trk->Px(),trk->Py(),trk->Pz(),trk->Charge()));

		Double_t d0z0[2],covd0z0[3];
		trk->PropagateToDCA(fVtx1,fBzkG,kVeryBig,d0z0,covd0z0);
		Double_t tx[3];
		trk->GetXYZ(tx);
		Double_t x1 = tx[0];
		Double_t y1 = tx[1];
		trk->GetPxPyPz(tx);
		Double_t px1 = tx[0];
		Double_t py1 = tx[1];
		Double_t pt1 = sqrt(px1*px1+py1*py1);
		px1 = px1/pt1;
		py1 = py1/pt1;

    TVector *varvec = new TVector(9);
    (*varvec)[0] = d0z0[0];
    (*varvec)[1] = x1;
    (*varvec)[2] = y1;
    (*varvec)[3] = px1;
    (*varvec)[4] = py1;
    (*varvec)[5] = fVtx1->GetX();
    (*varvec)[6] = fVtx1->GetY();
    (*varvec)[7] = fAnalCuts->GetPidHF()->GetPidResponse()->NumberOfSigmasTPC(trk,AliPID::kProton);
    (*varvec)[8] = fAnalCuts->GetPidHF()->GetPidResponse()->NumberOfSigmasTOF(trk,AliPID::kProton);
		fProtonCutVarsArray->AddLast(varvec);
	}

	Int_t pdgPr = -9999;
	if(fUseMCInfo)
	{
		Int_t labPr = trk->GetLabel();
		if(labPr<0) return;
		AliAODMCParticle *mcetrk = (AliAODMCParticle*)mcArray->At(labPr);
		if(!mcetrk) return;
		pdgPr = mcetrk->GetPdgCode();
		if(abs(pdgPr)!=2212) return;
		fHistoBachPtMCS->Fill(vpr.Pt(),vpr.Rapidity());
	}
}
////-------------------------------------------------------------------------------
void AliAnalysisTaskSELc2pK0sfromAODtracks::DefineV0TreeVariables() 
{
  //
  /// Define V0 tree variables
  //

  const char* nameoutput = GetOutputSlot(6)->GetContainer()->GetName();
  fV0VariablesTree = new TTree(nameoutput,"v0 variables tree");
  Int_t nVar = 3;
  fCandidateV0Variables = new Float_t [nVar];
  TString * fCandidateVariableNames = new TString[nVar];

  fCandidateVariableNames[ 0]="V0Px";
  fCandidateVariableNames[ 1]="V0Py";
  fCandidateVariableNames[ 2]="V0Pz";

  for (Int_t ivar=0; ivar<nVar; ivar++) {
    fV0VariablesTree->Branch(fCandidateVariableNames[ivar].Data(),&fCandidateV0Variables[ivar],Form("%s/f",fCandidateVariableNames[ivar].Data()));
  }

  return;
}
//________________________________________________________________________
void AliAnalysisTaskSELc2pK0sfromAODtracks::FillV0ROOTObjects(AliAODv0 *v0, TClonesArray *mcArray) 
{
  //
  /// Fill histograms or tree depending on fWriteVariableTree
  //
	if(!v0) return;

	TLorentzVector vk0;
	vk0.SetXYZM(v0->Px(),v0->Py(),v0->Pz(),v0->MassK0Short());
	fHistoK0sMassvsPt->Fill(v0->MassK0Short(),v0->Pt(),vk0.Rapidity());

	if(fDoEventMixing){
		Double_t ev0 = sqrt(v0->P2()+pow(v0->MassK0Short(),2));
		fV0Tracks->AddLast(new TLorentzVector(v0->Px(),v0->Py(),v0->Pz(),ev0));

    TVector *varvec = new TVector(1);
    (*varvec)[0] = v0->DcaV0ToPrimVertex();
		fV0CutVarsArray->AddLast(varvec);
	}

	Int_t v0pdgcode = -9999;
	if(fUseMCInfo)
	{
		Int_t pdgdgv0[2]={211,211};
		Int_t labV0 = v0->MatchToMC(310,mcArray,2,pdgdgv0); // the V0
		if(labV0<0) return;
		AliAODMCParticle *mcv0trk = (AliAODMCParticle*)mcArray->At(labV0);
		if(!mcv0trk) return;
		fHistoK0sMassvsPtMCS->Fill(v0->MassK0Short(),v0->Pt(),vk0.Rapidity());
	}
}
//_________________________________________________________________
Int_t AliAnalysisTaskSELc2pK0sfromAODtracks::GetPoolIndex(Double_t zvert, Double_t mult, Double_t rp){
	//
  // check in which of the pools the current event falls
	//

  Int_t theBinZ=TMath::BinarySearch(fNzVtxBins,fZvtxBins,zvert);
  if(theBinZ<0 || theBinZ>=fNzVtxBins) return -1;
  Int_t theBinM=TMath::BinarySearch(fNCentBins,fCentBins,mult);
  if(theBinM<0 || theBinM>=fNCentBins) return -1;
  Int_t theBinR=TMath::BinarySearch(fNRPBins,fRPBins,rp);
  if(theBinR<0 || theBinR>=fNRPBins) return -1;
  return fNRPBins*fNCentBins*theBinZ+fNRPBins*theBinM+theBinR;
}
//_________________________________________________________________
void AliAnalysisTaskSELc2pK0sfromAODtracks::ResetPool(Int_t poolIndex){
	//
  // delete the contets of the pool
	//
  if(poolIndex<0 || poolIndex>=fNOfPools) return;
  delete fEventBuffer[poolIndex];
  fEventBuffer[poolIndex]=new TTree(Form("EventBuffer_%d",poolIndex), "Temporary buffer for event mixing");

	fEventBuffer[poolIndex]->Branch("zVertex", &fVtxZ);
	fEventBuffer[poolIndex]->Branch("centrality", &fCentrality);
	fEventBuffer[poolIndex]->Branch("eventInfo", "TObjString",&fEventInfo);
	fEventBuffer[poolIndex]->Branch("v1array", "TObjArray", &fV0Tracks);
	fEventBuffer[poolIndex]->Branch("v1varsarray", "TObjArray", &fV0CutVarsArray);

  return;
}
//_________________________________________________________________
void AliAnalysisTaskSELc2pK0sfromAODtracks::DoEventMixingWithPools(Int_t poolIndex)
{
	//
  // perform mixed event analysis
	//

  if(poolIndex<0 || poolIndex>fNzVtxBins*fNCentBins*fNRPBins) return;
	if(fEventBuffer[poolIndex]->GetEntries()<fNumberOfEventsForMixing) return;

	Int_t nPro = fProtonTracks->GetEntries();
  Int_t nEvents=fEventBuffer[poolIndex]->GetEntries();
	Int_t nPro_test = fProtonCutVarsArray->GetEntries();
	if(nPro!=nPro_test){
		cout<<"Something wrong in mixing machinery"<<endl;
		exit(1);
	}

  TObjArray* v1array=0x0;
  TObjArray* v1varsarray=0x0;
  Float_t zVertex,cent;
  TObjString* eventInfo=0x0;
  fEventBuffer[poolIndex]->SetBranchAddress("eventInfo",&eventInfo);
  fEventBuffer[poolIndex]->SetBranchAddress("zVertex", &zVertex);
  fEventBuffer[poolIndex]->SetBranchAddress("centrality", &cent);
  fEventBuffer[poolIndex]->SetBranchAddress("v1array", &v1array);
  fEventBuffer[poolIndex]->SetBranchAddress("v1varsarray", &v1varsarray);
  for (Int_t i=0; i<nPro; i++)
  {
		TLorentzVector* trke=(TLorentzVector*) fProtonTracks->At(i);
    if(!trke)continue;
    TVector *prvarsarray = (TVector*)fProtonCutVarsArray->At(i);

		for(Int_t iEv=0; iEv<fNumberOfEventsForMixing; iEv++){
			fEventBuffer[poolIndex]->GetEvent(iEv + nEvents - fNumberOfEventsForMixing);
			Int_t nV01=v1array->GetEntries();
      for(Int_t iTr1=0; iTr1<nV01; iTr1++){
				TLorentzVector* v01=(TLorentzVector*)v1array->At(iTr1);
				if(!v01 ) continue;
				if(!fAnalCuts->SelectWithRoughCuts(v01,trke)) continue;

        TVector *v0varsarray = (TVector*) v1varsarray->At(iTr1);
				FillMixROOTObjects(trke,v01,prvarsarray,v0varsarray);
			}//v0 loop

		}//event loop
	}//track loop
}

////-------------------------------------------------------------------------------
void AliAnalysisTaskSELc2pK0sfromAODtracks::DefineMCTreeVariables() 
{
  ///
  /// Define electron tree variables
  ///

  const char* nameoutput = GetOutputSlot(7)->GetContainer()->GetName();
  fMCVariablesTree = new TTree(nameoutput,"MC variables tree");
  Int_t nVar = 1;
  fCandidateMCVariables = new Float_t [nVar];
  TString * fCandidateVariableNames = new TString[nVar];

  fCandidateVariableNames[ 0]="Centrality";

  for (Int_t ivar=0; ivar<nVar; ivar++) {
    fMCVariablesTree->Branch(fCandidateVariableNames[ivar].Data(),&fCandidateMCVariables[ivar],Form("%s/f",fCandidateVariableNames[ivar].Data()));
  }
  return;
}
////-------------------------------------------------------------------------------
void AliAnalysisTaskSELc2pK0sfromAODtracks::FillMCROOTObjects(AliAODMCParticle *mcpart, AliAODMCParticle *mcepart, AliAODMCParticle *mcv0part, Int_t decaytype) 
{
  //
  /// Fill histograms or tree depending on fWriteMCVariableTree
  //
	if(!mcpart) return;
	if(!mcepart) return;
	if(!mcv0part) return;


	if(decaytype==0){
		Double_t contmc[3];
		contmc[0] = mcpart->Pt();
		contmc[1] = mcpart->Y();
		contmc[2] = fCentrality;
		fHistoLcMCGen->Fill(contmc);
	}
}
////-------------------------------------------------------------------------------
void AliAnalysisTaskSELc2pK0sfromAODtracks::DefineMCProtonTreeVariables() 
{
  //
  // Define proton tree variables
  //

  const char* nameoutput = GetOutputSlot(9)->GetContainer()->GetName();
  fMCProtonVariablesTree = new TTree(nameoutput,"MC Proton variables tree");
  Int_t nVar = 1;
  fCandidateMCProtonVariables = new Float_t [nVar];
  TString * fCandidateVariableNames = new TString[nVar];

  fCandidateVariableNames[ 0]="Centrality";

  for (Int_t ivar=0; ivar<nVar; ivar++) {
    fMCProtonVariablesTree->Branch(fCandidateVariableNames[ivar].Data(),&fCandidateMCProtonVariables[ivar],Form("%s/f",fCandidateVariableNames[ivar].Data()));
  }
  return;
}
////-------------------------------------------------------------------------------
void AliAnalysisTaskSELc2pK0sfromAODtracks::FillMCProtonROOTObjects(AliAODMCParticle *mcepart, TClonesArray *mcArray) 
{
  //
  // Fill tree depending on fWriteMCVariableTree 
  //
	if(!mcepart) return;
}
////-------------------------------------------------------------------------------
void AliAnalysisTaskSELc2pK0sfromAODtracks::DefineMCV0TreeVariables() 
{
  //
  // Define Mc v0 tree variables
  //

  const char* nameoutput = GetOutputSlot(10)->GetContainer()->GetName();
  fMCV0VariablesTree = new TTree(nameoutput,"MC v0 variables tree");
  Int_t nVar = 1;
  fCandidateMCV0Variables = new Float_t [nVar];
  TString * fCandidateVariableNames = new TString[nVar];

  fCandidateVariableNames[ 0]="Centrality";

  for (Int_t ivar=0; ivar<nVar; ivar++) {
    fMCV0VariablesTree->Branch(fCandidateVariableNames[ivar].Data(),&fCandidateMCV0Variables[ivar],Form("%s/f",fCandidateVariableNames[ivar].Data()));
  }
  return;
}
////-------------------------------------------------------------------------------
void AliAnalysisTaskSELc2pK0sfromAODtracks::FillMCV0ROOTObjects(AliAODMCParticle *mcv0part, TClonesArray *mcArray) 
{
  //
  // Fill histograms or tree depending on fWriteMCVariableTree 
  //
	if(!mcv0part) return;
}
//_________________________________________________________________
Bool_t AliAnalysisTaskSELc2pK0sfromAODtracks::MakeMCAnalysis(TClonesArray *mcArray)
{
	//
  // Analyze AliAODmcparticle
	//
	Int_t nmcpart = mcArray->GetEntriesFast();

	for(Int_t i=0;i<nmcpart;i++)
	{
		AliAODMCParticle *mcpart = (AliAODMCParticle*) mcArray->At(i);
		if(TMath::Abs(mcpart->GetPdgCode())==4122){
			//cout<<"Lambdac"<<endl;
			Bool_t p_flag = kFALSE;
			Bool_t k0s_flag = kFALSE;
			AliAODMCParticle *mcepart = 0;
			AliAODMCParticle *mcv0part = 0;
			Int_t ndau = mcpart->GetLastDaughter()-mcpart->GetFirstDaughter()+1;
			if(ndau==2){
				for(Int_t idau=mcpart->GetFirstDaughter();idau<mcpart->GetLastDaughter()+1;idau++)
				{
					if(idau<0) break;
					AliAODMCParticle *mcdau = (AliAODMCParticle*) mcArray->At(idau);
					if(!mcdau) continue;
					if(TMath::Abs(mcdau->GetPdgCode())==2212){
						p_flag = kTRUE;
						mcepart = mcdau;
					}
					if(TMath::Abs(mcdau->GetPdgCode())==311){
						k0s_flag = kTRUE;
						mcv0part = mcdau;
					}
				}
			}

			Int_t decaytype = -9999;
			if(p_flag && k0s_flag) decaytype = 0;

			FillMCROOTObjects(mcpart,mcepart,mcv0part,decaytype);
		}

		if(TMath::Abs(mcpart->GetPdgCode())==2212 && mcpart->GetStatus()==1){
			AliESDtrackCuts *esdcuts = fAnalCuts->GetTrackCuts();
			Float_t etamin, etamax;
			esdcuts->GetEtaRange(etamin,etamax);
			if(fabs(mcpart->Eta())<etamax){
				fHistoBachPtMCGen->Fill(mcpart->Pt(),mcpart->Y());
			}
			FillMCProtonROOTObjects(mcpart, mcArray);
		}
		if(TMath::Abs(mcpart->GetPdgCode())==321 && mcpart->GetStatus()==1){
			AliESDtrackCuts *esdcuts = fAnalCuts->GetTrackCuts();
			Float_t etamin, etamax;
			esdcuts->GetEtaRange(etamin,etamax);
			if(fabs(mcpart->Eta())<etamax){
				fHistoKaonPtMCGen->Fill(mcpart->Pt(),mcpart->Y());
			}
		}
		if(TMath::Abs(mcpart->GetPdgCode())==310){
			Double_t etamin, etamax, rapmin, rapmax;
			fAnalCuts->GetProdV0EtaRange(etamin,etamax);
			fAnalCuts->GetProdV0RapRange(rapmin,rapmax);
			if((fabs(mcpart->Y())<rapmax) && (fabs(mcpart->Eta())<etamax)){
				fHistoK0sMassvsPtMCGen->Fill(0.497, mcpart->Pt(),mcpart->Y());
			}
			FillMCV0ROOTObjects(mcpart, mcArray);
		}
	}

	return kTRUE;
}
//________________________________________________________________________
Int_t AliAnalysisTaskSELc2pK0sfromAODtracks::MatchToMC(AliAODRecoCascadeHF *elobj, TClonesArray *mcArray, Int_t *pdgarray_pr, Int_t *pdgarray_v0, Int_t *labelarray_pr, Int_t *labelarray_v0,  Int_t &ngen_pr, Int_t &ngen_v0) 
{
  //
  // Match to MC
  //
	for(Int_t i=0;i<100;i++){
		pdgarray_pr[i] = -9999;
		labelarray_pr[i] = -9999;
		pdgarray_v0[i] = -9999;
		labelarray_v0[i] = -9999;
	}
	ngen_pr = 0;
	ngen_v0 = 0;

  AliVTrack *trk = dynamic_cast<AliVTrack*>(elobj->GetBachelor());
  if(!trk) return -1;
  Int_t labPr = trk->GetLabel();
	if(labPr<0) return -1;
	AliAODMCParticle *mcetrk = (AliAODMCParticle*)mcArray->At(labPr);
	if(!mcetrk) return -1;
	labelarray_pr[0] = labPr;
	pdgarray_pr[0] = mcetrk->GetPdgCode();
	ngen_pr ++;

  AliAODMCParticle *mcprimpr=0;
  mcprimpr = mcetrk;
  while(mcprimpr->GetMother()>=0) {
    Int_t labprim_pr=mcprimpr->GetMother();
    AliAODMCParticle *tmcprimpr = (AliAODMCParticle*)mcArray->At(labprim_pr);
    if(!tmcprimpr) {
			break;
    }

    mcprimpr = tmcprimpr;
		pdgarray_pr[ngen_pr] = mcprimpr->GetPdgCode();
		labelarray_pr[ngen_pr] = labprim_pr;
		ngen_pr ++;
		if(ngen_pr==100) break;
  }

  AliAODv0 *theV0 = dynamic_cast<AliAODv0*>(elobj->Getv0());
	if(!theV0) return -1;
	Int_t pdgdgv0[2]={211,211};
  Int_t labV0 = theV0->MatchToMC(310,mcArray,2,pdgdgv0); // the V0
	if(labV0<0) return -1;
	AliAODMCParticle *mcv0 = (AliAODMCParticle*)mcArray->At(labV0);
	if(!mcv0) return -1;
	labelarray_v0[0] = labV0;
	pdgarray_v0[0] = mcv0->GetPdgCode();
	ngen_v0 ++;

  AliAODMCParticle *mcprimv0=0;
  mcprimv0 = mcv0;
  while(mcprimv0->GetMother()>=0) {
    Int_t labprim_v0=mcprimv0->GetMother();
    AliAODMCParticle *tmcprimv0 = (AliAODMCParticle*)mcArray->At(labprim_v0);
    if(!tmcprimv0) {
			break;
    }

    mcprimv0 = tmcprimv0;
		pdgarray_v0[ngen_v0] = mcprimv0->GetPdgCode();
		labelarray_v0[ngen_v0] = labprim_v0;
		ngen_v0 ++;
		if(ngen_v0==100) break;
  }

	Bool_t same_flag = kFALSE;
	Int_t matchedlabel=-9999;
	for(Int_t iemc=0;iemc<ngen_pr;iemc++){
		for(Int_t ivmc=0;ivmc<ngen_v0;ivmc++){
			if(labelarray_pr[iemc]==labelarray_v0[ivmc]){
				same_flag = kTRUE;
				matchedlabel = labelarray_pr[iemc];
				break;
			}
		}
		if(same_flag) break;
	}

	return matchedlabel;

}

//________________________________________________________________________
void AliAnalysisTaskSELc2pK0sfromAODtracks::StoreGlobalTrackReference(AliAODTrack *track, Int_t index){
  //
  // Stores the pointer to the global track
  // copied from femtoscopy/k0analysis/plamanalysis
  //
  
  // Check that the id is positive
  if(track->GetID()<0){
    //    printf("Warning: track has negative ID: %d\n",track->GetID());
    return;
  }

  // Check id is not too big for buffer
  if(track->GetID()>=fTrackBuffSize){
    printf("Warning: track ID too big for buffer: ID: %d, buffer %d\n"
	   ,track->GetID(),fTrackBuffSize);
    return;
  }

  // Warn if we overwrite a track
  if(fGTI[track->GetID()]){
    // Seems like there are FilterMap 0 tracks
    // that have zero TPCNcls, don't store these!
    if( (!track->GetFilterMap()) &&
	(!track->GetTPCNcls())   )
      return;

    // Imagine the other way around,
    // the zero map zero clusters track
    // is stored and the good one wants 
    // to be added. We ommit the warning
    // and just overwrite the 'bad' track
    if( fGTI[track->GetID()]->GetFilterMap() ||
	fGTI[track->GetID()]->GetTPCNcls()   ){
      // If we come here, there's a problem
      printf("Warning! global track info already there!");
      printf("         TPCNcls track1 %u track2 %u",
	     (fGTI[track->GetID()])->GetTPCNcls(),track->GetTPCNcls());
      printf("         FilterMap track1 %u track2 %u\n",
	     (fGTI[track->GetID()])->GetFilterMap(),track->GetFilterMap());
    }
  } // Two tracks same id

  // // There are tracks with filter bit 0,
  // // do they have TPCNcls stored?
  // if(!track->GetFilterMap()){
  //   printf("Filter map is zero, TPCNcls: %u\n"
  // 	   ,track->GetTPCNcls());
  // }

  // Assign the pointer
  (fGTI[track->GetID()]) = track;
  (fGTIndex[track->GetID()]) = index;
}
//________________________________________________________________________
void AliAnalysisTaskSELc2pK0sfromAODtracks::ResetGlobalTrackReference(){
  // Sets all the pointers to zero. To be called at
  // the beginning or end of an event
  for(UShort_t i=0;i<fTrackBuffSize;i++){
    fGTI[i]=0;
    fGTIndex[i]=-9999;
  }
}
