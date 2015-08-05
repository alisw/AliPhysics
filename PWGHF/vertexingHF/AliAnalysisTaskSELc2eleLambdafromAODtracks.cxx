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
//               Lc->e Lambda  analysis code
//
//  Input: AOD
//  Output: TTree and/or THnSparse (mass vs pT vs Centrality)
//
//  Cuts:
//  TTree: SingleCuts on V0 and electron
//  THnSparse: In addition to that, IsSelected(obj, kCandidate) applied
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
#include <THnSparse.h>
#include <TLorentzVector.h>
#include <TTree.h>
#include "TROOT.h"
#include <TDatabasePDG.h>
#include <AliAnalysisDataSlot.h>
#include <AliAnalysisDataContainer.h>
#include "AliStack.h"
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
#include "AliAnalysisTaskSELc2eleLambdafromAODtracks.h"
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
#include "AliEventPoolManager.h"
#include "AliNormalizationCounter.h"

using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskSELc2eleLambdafromAODtracks)

//__________________________________________________________________________
AliAnalysisTaskSELc2eleLambdafromAODtracks::AliAnalysisTaskSELc2eleLambdafromAODtracks() : 
  AliAnalysisTaskSE(),
  fUseMCInfo(kFALSE),
  fOutput(0),
  fOutputAll(0),
  fListCuts(0),
  fCEvents(0),
  fHTrigger(0),
  fHCentrality(0),
  fAnalCuts(0),
  fIsEventSelected(kFALSE),
  fWriteVariableTree(kFALSE),
  fWriteEachVariableTree(kFALSE),
  fWriteMCVariableTree(kFALSE),
  fVariablesTree(0),
  fEleVariablesTree(0),
  fV0VariablesTree(0),
  fMCVariablesTree(0),
  fReconstructPrimVert(kFALSE),
  fIsMB(kFALSE),
  fIsSemi(kFALSE),
  fIsCent(kFALSE),
  fIsINT7(kFALSE),
  fIsEMC7(kFALSE),
  fCandidateVariables(),
  fCandidateEleVariables(),
  fCandidateV0Variables(),
  fCandidateMCVariables(),
  fVtx1(0),
  fV1(0),
  fVtxZ(0),
  fBzkG(0),
  fCentrality(0),
  fRunNumber(0),
  fTriggerCheck(0),
  fUseCentralityV0M(kFALSE),
  fEvNumberCounter(0),
  fHistoEleLambdaMass(0),
  fHistoEleLambdaMassRS(0),
  fHistoEleLambdaMassWS(0),
  fHistoEleLambdaMassRSMix(0),
  fHistoEleLambdaMassWSMix(0),
  fHistoEleLambdaMassvsElePtRS(0),
  fHistoEleLambdaMassvsElePtWS(0),
  fHistoEleLambdaMassvsElePtRSMix(0),
  fHistoEleLambdaMassvsElePtWSMix(0),
  fHistoElePtRS(0),
  fHistoElePtWS(0),
  fHistoElePtRSMix(0),
  fHistoElePtWSMix(0),
  fHistoEleLambdaMassMCS(0),
  fHistoEleLambdaMassMCGen(0),
  fHistoEleLambdaMassvsElePtMCS(0),
  fHistoEleLambdaMassvsElePtMCGen(0),
  fHistoElePtMCS(0),
  fHistoElePtMCGen(0),
  fHistoElePtvsEtaRS(0),
  fHistoElePtvsEtaWS(0),
  fHistoElePtvsEtaRSMix(0),
  fHistoElePtvsEtaWSMix(0),
  fHistoElePtvsEtaMCS(0),
  fHistoElePtvsEtaMCGen(0),
  fHistoElePtvsLambdaPtRS(0),
  fHistoElePtvsLambdaPtWS(0),
  fHistoElePtvsLambdaPtRSMix(0),
  fHistoElePtvsLambdaPtWSMix(0),
  fHistoElePtvsLambdaPtMCS(0),
  fHistoElePtvsLambdaPtvsLcPtMCS(0),
  fHistoElePtvsLambdaPtMCGen(0),
  fHistoElePtvsLambdaPtvsLcPtMCGen(0),
  fHistoElePtvsLambdaPtMCLcGen(0),
  fHistoElePtvsd0RS(0),
  fHistoElePtvsd0WS(0),
  fHistoElePtvsd0RSMix(0),
  fHistoElePtvsd0WSMix(0),
  fHistoElePtvsd0MCS(0),
  fHistoElePtvsd0PromptMCS(0),
  fHistoElePtvsd0BFeeddownMCS(0),
  fHistoEleLambdaMassFeeddownXic0MCS(0),
  fHistoEleLambdaMassFeeddownXic0MCGen(0),
  fHistoEleLambdaMassvsElePtFeeddownXic0MCS(0),
  fHistoEleLambdaMassvsElePtFeeddownXic0MCGen(0),
  fHistoElePtFeeddownXic0MCS(0),
  fHistoElePtFeeddownXic0MCGen(0),
  fHistoElePtvsEtaFeeddownXic0MCS(0),
  fHistoElePtvsEtaFeeddownXic0MCGen(0),
  fHistoElePtvsLambdaPtFeeddownXic0MCS(0),
  fHistoElePtvsLambdaPtFeeddownXic0MCGen(0),
  fHistoEleLambdaMassFeeddownXicPlusMCS(0),
  fHistoEleLambdaMassFeeddownXicPlusMCGen(0),
  fHistoEleLambdaMassvsElePtFeeddownXicPlusMCS(0),
  fHistoEleLambdaMassvsElePtFeeddownXicPlusMCGen(0),
  fHistoElePtFeeddownXicPlusMCS(0),
  fHistoElePtFeeddownXicPlusMCGen(0),
  fHistoElePtvsEtaFeeddownXicPlusMCS(0),
  fHistoElePtvsEtaFeeddownXicPlusMCGen(0),
  fHistoElePtvsLambdaPtFeeddownXicPlusMCS(0),
  fHistoElePtvsLambdaPtFeeddownXicPlusMCGen(0),
  fHistoBachPt(0),
  fHistoBachPtMCS(0),
  fHistoBachPtMCGen(0),
  fHistod0Bach(0),
  fHistoLambdaMassvsPt(0),
  fHistoLambdaMassvsPtMCS(0),
  fHistoLambdaMassvsPtMCGen(0),
  fHistoK0sMassvsPt(0),
  fHistoElectronTPCPID(0),
  fHistoElectronTOFPID(0),
  fHistoElectronTPCSelPID(0),
  fHistoElectronTOFSelPID(0),
  fHistoElectronTPCPIDSelTOF(0),
  fHistoElectronTPCPIDSelTOFSmallEta(0),
  fHistoElectronTPCPIDSelTOFLargeEta(0),
	fHistoElectronQovPtvsPhi(0),
	fHistoLambdaQovPtvsPhi(0),
	fCounter(0),
	fHistonEvtvsRunNumber(0),
	fHistonElevsRunNumber(0),
	fHistonLambdavsRunNumber(0),
  fDoEventMixing(0),
	fNumberOfEventsForMixing		(5),
	fNzVtxBins					(0), 
	fNCentBins					(0),
	fNOfPools(1),
	fEventBuffer(0x0),
	fEventInfo(0x0),
	fElectronTracks(0x0)
{
  //
  // Default Constructor. 
  //
	for(Int_t i=0;i<17;i++){
		fHistoElePtvsCutVarsRS[i] = 0;
		fHistoElePtvsCutVarsWS[i] = 0;
		fHistoElePtvsCutVarsMCS[i] = 0;
	}
	for(Int_t i=0;i<8;i++){
		fHistoElectronTPCPIDSelTOFEtaDep[i] = 0;
	}
}

//___________________________________________________________________________
AliAnalysisTaskSELc2eleLambdafromAODtracks::AliAnalysisTaskSELc2eleLambdafromAODtracks(const Char_t* name,
									     AliRDHFCutsLctoeleLambdafromAODtracks* analCuts, 
									     Bool_t writeVariableTree) :
  AliAnalysisTaskSE(name),
  fUseMCInfo(kFALSE),
  fOutput(0),
  fOutputAll(0),
  fListCuts(0),
  fCEvents(0),
  fHTrigger(0),
  fHCentrality(0),
  fAnalCuts(analCuts),
  fIsEventSelected(kFALSE),
  fWriteVariableTree(writeVariableTree),
  fWriteEachVariableTree(kFALSE),
  fWriteMCVariableTree(kFALSE),
  fVariablesTree(0),
  fEleVariablesTree(0),
  fV0VariablesTree(0),
  fMCVariablesTree(0),
  fReconstructPrimVert(kFALSE),
  fIsMB(kFALSE),
  fIsSemi(kFALSE),
  fIsCent(kFALSE),
  fIsINT7(kFALSE),
  fIsEMC7(kFALSE),
  fCandidateVariables(),
  fCandidateEleVariables(),
  fCandidateV0Variables(),
  fCandidateMCVariables(),
  fVtx1(0),
  fV1(0),
  fVtxZ(0),
  fBzkG(0),
  fCentrality(0),
  fRunNumber(0),
  fTriggerCheck(0),
  fUseCentralityV0M(kFALSE),
  fEvNumberCounter(0),
  fHistoEleLambdaMass(0),
  fHistoEleLambdaMassRS(0),
  fHistoEleLambdaMassWS(0),
  fHistoEleLambdaMassRSMix(0),
  fHistoEleLambdaMassWSMix(0),
  fHistoEleLambdaMassvsElePtRS(0),
  fHistoEleLambdaMassvsElePtWS(0),
  fHistoEleLambdaMassvsElePtRSMix(0),
  fHistoEleLambdaMassvsElePtWSMix(0),
  fHistoElePtRS(0),
  fHistoElePtWS(0),
  fHistoElePtRSMix(0),
  fHistoElePtWSMix(0),
  fHistoEleLambdaMassMCS(0),
  fHistoEleLambdaMassMCGen(0),
  fHistoEleLambdaMassvsElePtMCS(0),
  fHistoEleLambdaMassvsElePtMCGen(0),
  fHistoElePtMCS(0),
  fHistoElePtMCGen(0),
  fHistoElePtvsEtaRS(0),
  fHistoElePtvsEtaWS(0),
  fHistoElePtvsEtaRSMix(0),
  fHistoElePtvsEtaWSMix(0),
  fHistoElePtvsEtaMCS(0),
  fHistoElePtvsEtaMCGen(0),
  fHistoElePtvsLambdaPtRS(0),
  fHistoElePtvsLambdaPtWS(0),
  fHistoElePtvsLambdaPtRSMix(0),
  fHistoElePtvsLambdaPtWSMix(0),
  fHistoElePtvsLambdaPtMCS(0),
  fHistoElePtvsLambdaPtvsLcPtMCS(0),
  fHistoElePtvsLambdaPtMCGen(0),
  fHistoElePtvsLambdaPtvsLcPtMCGen(0),
  fHistoElePtvsLambdaPtMCLcGen(0),
  fHistoElePtvsd0RS(0),
  fHistoElePtvsd0WS(0),
  fHistoElePtvsd0RSMix(0),
  fHistoElePtvsd0WSMix(0),
  fHistoElePtvsd0MCS(0),
  fHistoElePtvsd0PromptMCS(0),
  fHistoElePtvsd0BFeeddownMCS(0),
  fHistoEleLambdaMassFeeddownXic0MCS(0),
  fHistoEleLambdaMassFeeddownXic0MCGen(0),
  fHistoEleLambdaMassvsElePtFeeddownXic0MCS(0),
  fHistoEleLambdaMassvsElePtFeeddownXic0MCGen(0),
  fHistoElePtFeeddownXic0MCS(0),
  fHistoElePtFeeddownXic0MCGen(0),
  fHistoElePtvsEtaFeeddownXic0MCS(0),
  fHistoElePtvsEtaFeeddownXic0MCGen(0),
  fHistoElePtvsLambdaPtFeeddownXic0MCS(0),
  fHistoElePtvsLambdaPtFeeddownXic0MCGen(0),
  fHistoEleLambdaMassFeeddownXicPlusMCS(0),
  fHistoEleLambdaMassFeeddownXicPlusMCGen(0),
  fHistoEleLambdaMassvsElePtFeeddownXicPlusMCS(0),
  fHistoEleLambdaMassvsElePtFeeddownXicPlusMCGen(0),
  fHistoElePtFeeddownXicPlusMCS(0),
  fHistoElePtFeeddownXicPlusMCGen(0),
  fHistoElePtvsEtaFeeddownXicPlusMCS(0),
  fHistoElePtvsEtaFeeddownXicPlusMCGen(0),
  fHistoElePtvsLambdaPtFeeddownXicPlusMCS(0),
  fHistoElePtvsLambdaPtFeeddownXicPlusMCGen(0),
  fHistoBachPt(0),
  fHistoBachPtMCS(0),
  fHistoBachPtMCGen(0),
  fHistod0Bach(0),
  fHistoLambdaMassvsPt(0),
  fHistoLambdaMassvsPtMCS(0),
  fHistoLambdaMassvsPtMCGen(0),
  fHistoK0sMassvsPt(0),
  fHistoElectronTPCPID(0),
  fHistoElectronTOFPID(0),
  fHistoElectronTPCSelPID(0),
  fHistoElectronTOFSelPID(0),
  fHistoElectronTPCPIDSelTOF(0),
  fHistoElectronTPCPIDSelTOFSmallEta(0),
  fHistoElectronTPCPIDSelTOFLargeEta(0),
	fHistoElectronQovPtvsPhi(0),
	fHistoLambdaQovPtvsPhi(0),
	fCounter(0),
	fHistonEvtvsRunNumber(0),
	fHistonElevsRunNumber(0),
	fHistonLambdavsRunNumber(0),
  fDoEventMixing(0),
	fNumberOfEventsForMixing		(5),
	fNzVtxBins					(0), 
	fNCentBins					(0),
	fNOfPools(1),
	fEventBuffer(0x0),
	fEventInfo(0x0),
	fElectronTracks(0x0)
{
  //
  // Constructor. Initialization of Inputs and Outputs
  //
  Info("AliAnalysisTaskSELc2eleLambdafromAODtracks","Calling Constructor");

	for(Int_t i=0;i<17;i++){
		fHistoElePtvsCutVarsRS[i] = 0;
		fHistoElePtvsCutVarsWS[i] = 0;
		fHistoElePtvsCutVarsMCS[i] = 0;
	}
	for(Int_t i=0;i<8;i++){
		fHistoElectronTPCPIDSelTOFEtaDep[i] = 0;
	}

  DefineOutput(1,TList::Class());  //conters
  DefineOutput(2,TList::Class());
  DefineOutput(3,TList::Class());  //conters
  DefineOutput(4,TTree::Class());  //My private output
  DefineOutput(5,TTree::Class());  //My private output
  DefineOutput(6,TTree::Class());  //My private output
  DefineOutput(7,TTree::Class());  //My private output
  DefineOutput(8,AliNormalizationCounter::Class());
}

//___________________________________________________________________________
AliAnalysisTaskSELc2eleLambdafromAODtracks::~AliAnalysisTaskSELc2eleLambdafromAODtracks() {
  //
  // destructor
  //
  Info("~AliAnalysisTaskSELc2eleLambdafromAODtracks","Calling Destructor");

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
  if (fEleVariablesTree) {
    delete fEleVariablesTree;
    fEleVariablesTree = 0;
  }
  if (fV0VariablesTree) {
    delete fV0VariablesTree;
    fV0VariablesTree = 0;
  }
  if (fMCVariablesTree) {
    delete fMCVariablesTree;
    fMCVariablesTree = 0;
  }
	if(fCounter){
		delete fCounter;
		fCounter = 0;
	}
}

//_________________________________________________
void AliAnalysisTaskSELc2eleLambdafromAODtracks::Init() {
  //
  // Initialization
  //
  //

  fIsEventSelected=kFALSE;

  if (fDebug > 1) AliInfo("Init");

  fListCuts = new TList();
  fListCuts->SetOwner();
  fListCuts->SetName("ListCuts");
  fListCuts->Add(new AliRDHFCutsLctoeleLambdafromAODtracks(*fAnalCuts));
  PostData(2,fListCuts);

  return;
}

//_________________________________________________
void AliAnalysisTaskSELc2eleLambdafromAODtracks::UserExec(Option_t *)
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

  //------------------------------------------------
  // First check if the event has proper B
  //------------------------------------------------
	
  fBzkG = (Double_t)aodEvent->GetMagneticField(); 
  AliKFParticle::SetField(fBzkG);
  if (TMath::Abs(fBzkG)<0.001) {
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
      AliError("AliAnalysisTaskSELc2eleLambdafromAODtracks::UserExec: MC header branch not found!\n");
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
			MakeMCAnalysis(mcArray);
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
	fEvNumberCounter++;

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
  MakeAnalysis(aodEvent,mcArray);


  PostData(1,fOutput);
  PostData(3,fOutputAll);
  PostData(4,fVariablesTree);
  PostData(5,fEleVariablesTree);
  PostData(6,fV0VariablesTree);
  PostData(7,fMCVariablesTree);
  PostData(8,fCounter);    

  fIsEventSelected=kFALSE;

  delete fV1;
  return;
}

//________________________________________ terminate ___________________________
void AliAnalysisTaskSELc2eleLambdafromAODtracks::Terminate(Option_t*)
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
void AliAnalysisTaskSELc2eleLambdafromAODtracks::UserCreateOutputObjects() 
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

  DefineEleTreeVariables();
  PostData(5,fEleVariablesTree);

  DefineV0TreeVariables();
  PostData(6,fV0VariablesTree);

  DefineMCTreeVariables();
  PostData(7,fMCVariablesTree);

  //Counter for Normalization
  TString normName="NormalizationCounter";
  AliAnalysisDataContainer *cont = GetOutputSlot(8)->GetContainer();
  if(cont)normName=(TString)cont->GetName();
  fCounter = new AliNormalizationCounter(normName.Data());
  fCounter->Init();
  PostData(8,fCounter);

	if(fDoEventMixing){
		fElectronTracks = new TObjArray();
		fElectronTracks->SetOwner();

		fNOfPools=fNCentBins*fNzVtxBins;
		fEventBuffer = new TTree*[fNOfPools];
		for(Int_t i=0; i<fNOfPools; i++){
			fEventBuffer[i]=new TTree(Form("EventBuffer_%d",i), "Temporary buffer for event mixing");
			fEventBuffer[i]->Branch("zVertex", &fVtxZ);
			fEventBuffer[i]->Branch("centrality", &fCentrality);
			fEventBuffer[i]->Branch("eventInfo", "TObjString",&fEventInfo);
			fEventBuffer[i]->Branch("earray", "TObjArray", &fElectronTracks);
		}
	}


  return;
}

//-------------------------------------------------------------------------------
void AliAnalysisTaskSELc2eleLambdafromAODtracks::MakeAnalysis
(
 AliAODEvent *aodEvent, TClonesArray *mcArray
 )
{
  //
  // Main Analysis part
  //
	if(fDoEventMixing && fElectronTracks) fElectronTracks->Delete();

  //------------------------------------------------
  // Select good track before hand to save time
  //------------------------------------------------

  Int_t nV0s= aodEvent->GetNumberOfV0s();
  Int_t nTracks= aodEvent->GetNumberOfTracks();

  Bool_t  seleTrkFlags[nTracks];
  Int_t nSeleTrks=0;
  SelectTrack(aodEvent,nTracks,nSeleTrks,seleTrkFlags,mcArray);

  Bool_t  seleV0Flags[nV0s];
  Int_t     nSeleV0=0;
  SelectV0(aodEvent,nV0s,nSeleV0,seleV0Flags,mcArray);

	Int_t runnumber_offset = 0;
	Int_t runnumber = aodEvent->GetRunNumber();
	if(runnumber<=131000&&runnumber>=114000){
		runnumber_offset = 114000;//lhc10bcde
	}else if(runnumber<=196000&&runnumber>=195000){
		runnumber_offset = 195000;//lhc13bc
	}else if(runnumber<=170593&&runnumber>=167902){
		runnumber_offset = 167902;//lhc11h
	}
	fHistonElevsRunNumber->Fill(runnumber-runnumber_offset,nSeleTrks);
	fHistonLambdavsRunNumber->Fill(runnumber-runnumber_offset,nSeleV0);

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
      if(!seleTrkFlags[itrk]) continue;
      AliAODTrack *trk = (AliAODTrack*)aodEvent->GetTrack(itrk);
      if(trk->GetID()<0) continue;

      Int_t cpid = cptrack->GetID();
      Int_t cnid = cntrack->GetID();
      Int_t lpid = trk->GetID();
      if((cpid==lpid)||(cnid==lpid)) continue;

      //if(!fAnalCuts->SelectWithRoughCuts(v0,trk)) continue;

      AliAODVertex *secVert = ReconstructSecondaryVertex(v0,trk,aodEvent);//Fake, prim vertex is just used as secondary vertex. place holder for future
      if(!secVert) continue;

      AliAODRecoCascadeHF *elobj = MakeCascadeHF(v0,trk,aodEvent,secVert,false);
      if(!elobj) {
	continue;
      }

      FillROOTObjects(elobj, v0,trk,mcArray,false);

      elobj->GetSecondaryVtx()->RemoveDaughters();
      elobj->UnsetOwnPrimaryVtx();
      delete elobj;elobj=NULL;
      delete secVert;
    }
  }

  if(fDoEventMixing){
		fEventInfo->SetString(Form("Ev%d_esd%d_E%d",AliAnalysisManager::GetAnalysisManager()->GetNcalls(),((AliAODHeader*)aodEvent->GetHeader())->GetEventNumberESDFile(),fElectronTracks->GetEntries()));
    Int_t ind=GetPoolIndex(fVtxZ,fCentrality);
    if(ind>=0 && ind<fNOfPools){
      if(fEventBuffer[ind]->GetEntries() >= fNumberOfEventsForMixing){
				DoEventMixingWithPools(ind,aodEvent,seleV0Flags);
				//ResetPool(ind);
      }
      fEventBuffer[ind]->Fill();
    }
  }
}
////-------------------------------------------------------------------------------
void AliAnalysisTaskSELc2eleLambdafromAODtracks::DefineTreeVariables() 
{
  //
  // Define tree variables
  //

  const char* nameoutput = GetOutputSlot(4)->GetContainer()->GetName();
  fVariablesTree = new TTree(nameoutput,"Candidates variables tree");
  Int_t nVar = 69;
  fCandidateVariables = new Float_t [nVar];
  TString * fCandidateVariableNames = new TString[nVar];

  fCandidateVariableNames[ 0]="Centrality";
  fCandidateVariableNames[ 1]="InvMassEleLambda";
  fCandidateVariableNames[ 2]="EleLambdaPt";
  fCandidateVariableNames[ 3]="EleLambdaPx";
  fCandidateVariableNames[ 4]="EleLambdaPy";
  fCandidateVariableNames[ 5]="EleLambdaPz";
  fCandidateVariableNames[ 6]="ElePx";
  fCandidateVariableNames[ 7]="ElePy";
  fCandidateVariableNames[ 8]="ElePz";
  fCandidateVariableNames[ 9]="V0Px";
  fCandidateVariableNames[10]="V0Py";
  fCandidateVariableNames[11]="V0Pz";
  fCandidateVariableNames[12]="AntiLambdaFlag";
  fCandidateVariableNames[13]="MassLambda";
  fCandidateVariableNames[14]="MassAntiLambda";
  fCandidateVariableNames[15]="Eled0";
  fCandidateVariableNames[16]="V0d0";
  fCandidateVariableNames[17]="nSigmaTPCele";
  fCandidateVariableNames[18]="nSigmaTOFele";
  fCandidateVariableNames[19]="nSigmaTPCv0pr";
  fCandidateVariableNames[20]="nSigmaTOFv0pr";
  fCandidateVariableNames[21]="EleCharge";
  fCandidateVariableNames[22]="ProtonPx";
  fCandidateVariableNames[23]="ProtonPy";
  fCandidateVariableNames[24]="ProtonPz";
  fCandidateVariableNames[25]="PiPx";
  fCandidateVariableNames[26]="PiPy";
  fCandidateVariableNames[27]="PiPz";
  fCandidateVariableNames[28]="mcpdglc";
  fCandidateVariableNames[29]="mclablc";
  fCandidateVariableNames[30]="mcpdgmomele";
  fCandidateVariableNames[31]="mcpdgmomv0";
  fCandidateVariableNames[32]="Mixing";
  fCandidateVariableNames[33]="mcpdgele";
  fCandidateVariableNames[34]="nSigmaTPCpr_etrk";
  fCandidateVariableNames[35]="nSigmaTOFpr_etrk";
  fCandidateVariableNames[36]="nSigmaTPCka_etrk";
  fCandidateVariableNames[37]="nSigmaTOFka_etrk";
  fCandidateVariableNames[38]="MassK0Short";
  fCandidateVariableNames[39]="mcpdggrmomele";
  fCandidateVariableNames[40]="mcpdggrmomv0";
  fCandidateVariableNames[41]="mcngenele";
  fCandidateVariableNames[42]="mcngenv0";
  fCandidateVariableNames[43]="mclcpx";
  fCandidateVariableNames[44]="mclcpy";
  fCandidateVariableNames[45]="mclcpz";
  fCandidateVariableNames[46]="mcelepx";
  fCandidateVariableNames[47]="mcelepy";
  fCandidateVariableNames[48]="mcelepz";
  fCandidateVariableNames[49]="mcv0px";
  fCandidateVariableNames[50]="mcv0py";
  fCandidateVariableNames[51]="mcv0pz";
  fCandidateVariableNames[52]="nSigmaTPCpi_etrk";
  fCandidateVariableNames[53]="nSigmaTOFpi_etrk";
  fCandidateVariableNames[54]="PrimVertx";
  fCandidateVariableNames[55]="PrimVerty";
  fCandidateVariableNames[56]="PrimVertz";
  fCandidateVariableNames[57]="V0Vertx";
  fCandidateVariableNames[58]="V0Verty";
  fCandidateVariableNames[59]="V0Vertz";

  fCandidateVariableNames[60]="DcaV0PrToPrimVertex";
  fCandidateVariableNames[61]="DcaV0PiToPrimVertex";
  fCandidateVariableNames[62]="DcaV0daughters";
  fCandidateVariableNames[63]="V0CosPointingAngle";
  fCandidateVariableNames[64]="V0ProperDecayLength";
  fCandidateVariableNames[65]="MassK0Short";

  fCandidateVariableNames[66]="nSigmaTPCv0pi";
  fCandidateVariableNames[67]="nSigmaTOFv0pi";

  fCandidateVariableNames[68]="EvNumber";

  for (Int_t ivar=0; ivar<nVar; ivar++) {
    fVariablesTree->Branch(fCandidateVariableNames[ivar].Data(),&fCandidateVariables[ivar],Form("%s/f",fCandidateVariableNames[ivar].Data()));
  }

  return;
}

////-------------------------------------------------------------------------------
void AliAnalysisTaskSELc2eleLambdafromAODtracks::FillROOTObjects(AliAODRecoCascadeHF *elobj, AliAODv0 *v0, AliAODTrack *trk, TClonesArray *mcArray, Bool_t mixing_flag) 
{
  //
  // Fill histograms or tree depending on fWriteVariableTree
  //
	if(!trk) return;
	if(!v0) return;

	for(Int_t i=0;i<67;i++){
		fCandidateVariables[i] = -9999.;
	}

	Bool_t anti_lambda_flag = kFALSE;
	if(fabs(v0->MassAntiLambda()-1.115683)<fAnalCuts->GetProdV0MassTolLambda()) anti_lambda_flag = kTRUE;

  AliAODTrack *cptrack =  (AliAODTrack*)(v0->GetDaughter(0));
  AliAODTrack *cntrack =  (AliAODTrack*)(v0->GetDaughter(1));
	if(cptrack->Charge()<0 && cntrack->Charge()>0){
		cptrack =  (AliAODTrack*)(v0->GetDaughter(1));
		cntrack =  (AliAODTrack*)(v0->GetDaughter(0));
	}

  fCandidateVariables[ 0] = fCentrality;
	UInt_t pdgdg[2]={11,3122};
  fCandidateVariables[ 1] = elobj->InvMass(2,pdgdg);
  fCandidateVariables[ 2] = elobj->Pt();
  fCandidateVariables[ 3] = elobj->Px();
  fCandidateVariables[ 4] = elobj->Py();
  fCandidateVariables[ 5] = elobj->Pz();
  fCandidateVariables[ 6] = elobj->PxProng(0);
  fCandidateVariables[ 7] = elobj->PyProng(0);
  fCandidateVariables[ 8] = elobj->PzProng(0);
  fCandidateVariables[ 9] = elobj->PxProng(1);
  fCandidateVariables[10] = elobj->PyProng(1);
  fCandidateVariables[11] = elobj->PzProng(1);
  fCandidateVariables[12] = anti_lambda_flag;
  fCandidateVariables[13] = v0->MassLambda();
  fCandidateVariables[14] = v0->MassAntiLambda();
  fCandidateVariables[15] = elobj->Getd0Prong(0);
  fCandidateVariables[16] = elobj->Getd0Prong(1);

  Double_t nSigmaTPCele=-9999.;
  Double_t nSigmaTOFele=-9999.;
  Double_t nSigmaTPCv0pr=-9999.;
  Double_t nSigmaTOFv0pr=-9999.;
  Double_t nSigmaTPCv0pi=-9999.;
  Double_t nSigmaTOFv0pi=-9999.;
  if(fAnalCuts->GetIsUsePID()&&!mixing_flag)
  {
		nSigmaTPCele = fAnalCuts->GetPidHF()->GetPidResponse()->NumberOfSigmasTPC(trk,AliPID::kElectron);
		nSigmaTOFele = fAnalCuts->GetPidHF()->GetPidResponse()->NumberOfSigmasTOF(trk,AliPID::kElectron);
    fCandidateVariables[17] = nSigmaTPCele;
    fCandidateVariables[18] = nSigmaTOFele;
  }

	if(fAnalCuts->GetUseLambdaPID()&&!mixing_flag)
	{
		if(anti_lambda_flag){
			nSigmaTPCv0pr = fAnalCuts->GetPidProton()->GetPidResponse()->NumberOfSigmasTPC(cntrack,AliPID::kProton);
			nSigmaTOFv0pr = fAnalCuts->GetPidProton()->GetPidResponse()->NumberOfSigmasTOF(cntrack,AliPID::kProton);
			nSigmaTPCv0pi = fAnalCuts->GetPidProton()->GetPidResponse()->NumberOfSigmasTPC(cptrack,AliPID::kPion);
			nSigmaTOFv0pi = fAnalCuts->GetPidProton()->GetPidResponse()->NumberOfSigmasTOF(cptrack,AliPID::kPion);
		}else{
			nSigmaTPCv0pr = fAnalCuts->GetPidProton()->GetPidResponse()->NumberOfSigmasTPC(cptrack,AliPID::kProton);
			nSigmaTOFv0pr = fAnalCuts->GetPidProton()->GetPidResponse()->NumberOfSigmasTOF(cptrack,AliPID::kProton);
			nSigmaTPCv0pi = fAnalCuts->GetPidProton()->GetPidResponse()->NumberOfSigmasTPC(cntrack,AliPID::kPion);
			nSigmaTOFv0pi = fAnalCuts->GetPidProton()->GetPidResponse()->NumberOfSigmasTOF(cntrack,AliPID::kPion);
		}
      fCandidateVariables[19] = nSigmaTPCv0pr;
      fCandidateVariables[20] = nSigmaTOFv0pr;
      fCandidateVariables[66] = nSigmaTPCv0pi;
      fCandidateVariables[67] = nSigmaTOFv0pi;
  }
  fCandidateVariables[21] = trk->Charge();

	if(anti_lambda_flag){
		fCandidateVariables[22] = cntrack->Px();
		fCandidateVariables[23] = cntrack->Py();
		fCandidateVariables[24] = cntrack->Pz();
		fCandidateVariables[25] = cptrack->Px();
		fCandidateVariables[26] = cptrack->Py();
		fCandidateVariables[27] = cptrack->Pz();
	}else{
		fCandidateVariables[22] = cptrack->Px();
		fCandidateVariables[23] = cptrack->Py();
		fCandidateVariables[24] = cptrack->Pz();
		fCandidateVariables[25] = cntrack->Px();
		fCandidateVariables[26] = cntrack->Py();
		fCandidateVariables[27] = cntrack->Pz();
	}

  AliAODMCParticle *mclc = 0;
  AliAODMCParticle *mcele = 0;
  AliAODMCParticle *mcv0 = 0;
  Int_t mclablc = 0;
	Int_t mcpdgele_array[100];
	Int_t mcpdgv0_array[100];
	Int_t mclabelele_array[100];
	Int_t mclabelv0_array[100];
	Int_t mcngen_ele=-9999;
	Int_t mcngen_v0=-9999;

	if(fUseMCInfo && mcArray){

		mclablc =  MatchToMC(elobj,mcArray,mcpdgele_array, mcpdgv0_array,mclabelele_array,mclabelv0_array,mcngen_ele,mcngen_v0);

    if(mclablc>-1){
      mclc = (AliAODMCParticle*) mcArray->At(mclablc);
			if(mclabelele_array[0]>=0)
				mcele = (AliAODMCParticle*) mcArray->At(mclabelele_array[0]);
			if(mclabelv0_array[0]>=0)
				mcv0 = (AliAODMCParticle*) mcArray->At(mclabelv0_array[0]);
			if(mclc){
				fCandidateVariables[28] = mclc->GetPdgCode();
				fCandidateVariables[29] = mclc->Label();
				fCandidateVariables[43] = mclc->Px();
				fCandidateVariables[44] = mclc->Py();
				fCandidateVariables[45] = mclc->Pz();
			}
			if(mcele){
				fCandidateVariables[46] = mcele->Px();
				fCandidateVariables[47] = mcele->Py();
				fCandidateVariables[48] = mcele->Pz();
			}
			if(mcv0){
				fCandidateVariables[49] = mcv0->Px();
				fCandidateVariables[50] = mcv0->Py();
				fCandidateVariables[51] = mcv0->Pz();
			}
		}
		fCandidateVariables[30] = mcpdgele_array[1];
		fCandidateVariables[31] = mcpdgv0_array[1];
		fCandidateVariables[33] = mcpdgele_array[0];
		fCandidateVariables[39] = mcpdgele_array[2];
		fCandidateVariables[40] = mcpdgv0_array[2];
		fCandidateVariables[41] = mcngen_ele;
		fCandidateVariables[42] = mcngen_v0;
	}
	fCandidateVariables[32] = mixing_flag;

  if(fAnalCuts->GetIsUsePID()&&!mixing_flag)
  {
		Double_t nSigmaTPCpr_etrk = fAnalCuts->GetPidHF()->GetPidResponse()->NumberOfSigmasTPC(trk,AliPID::kProton);
		Double_t nSigmaTOFpr_etrk = fAnalCuts->GetPidHF()->GetPidResponse()->NumberOfSigmasTOF(trk,AliPID::kProton);
		Double_t nSigmaTPCka_etrk = fAnalCuts->GetPidHF()->GetPidResponse()->NumberOfSigmasTPC(trk,AliPID::kKaon);
		Double_t nSigmaTOFka_etrk = fAnalCuts->GetPidHF()->GetPidResponse()->NumberOfSigmasTOF(trk,AliPID::kKaon);
		Double_t nSigmaTPCpi_etrk = fAnalCuts->GetPidHF()->GetPidResponse()->NumberOfSigmasTPC(trk,AliPID::kPion);
		Double_t nSigmaTOFpi_etrk = fAnalCuts->GetPidHF()->GetPidResponse()->NumberOfSigmasTOF(trk,AliPID::kPion);
    fCandidateVariables[34] = nSigmaTPCpr_etrk;
    fCandidateVariables[35] = nSigmaTOFpr_etrk;
    fCandidateVariables[36] = nSigmaTPCka_etrk;
    fCandidateVariables[37] = nSigmaTOFka_etrk;
    fCandidateVariables[52] = nSigmaTPCpi_etrk;
    fCandidateVariables[53] = nSigmaTOFpi_etrk;
  }
  fCandidateVariables[38] = v0->MassK0Short();

  fCandidateVariables[54] = fVtx1->GetX();
  fCandidateVariables[55] = fVtx1->GetY();
  fCandidateVariables[56] = fVtx1->GetZ();
  fCandidateVariables[57] = v0->DecayVertexV0X();
  fCandidateVariables[58] = v0->DecayVertexV0Y();
  fCandidateVariables[59] = v0->DecayVertexV0Z();

	Double_t lDcaPosToPrimVertex = v0->DcaPosToPrimVertex();
	Double_t lDcaNegToPrimVertex = v0->DcaNegToPrimVertex();
  if(!anti_lambda_flag){
		fCandidateVariables[60] = lDcaPosToPrimVertex;
		fCandidateVariables[61] = lDcaNegToPrimVertex;
  }else{
		fCandidateVariables[60] = lDcaNegToPrimVertex;
		fCandidateVariables[61] = lDcaPosToPrimVertex;
  }
	fCandidateVariables[62] = v0->DcaV0Daughters();
  Double_t posVtx[3] = {0.,0.,0.};
  fVtx1->GetXYZ(posVtx);
  fCandidateVariables[63] = v0->CosPointingAngle(posVtx); 
  Double_t ptotlam = TMath::Sqrt(pow(v0->Px(),2)+pow(v0->Py(),2)+pow(v0->Pz(),2));
  fCandidateVariables[64] = v0->DecayLengthV0(posVtx)*1.1157/ptotlam;
  fCandidateVariables[65] = v0->MassK0Short();

  fCandidateVariables[68] = fEvNumberCounter;


  if(fWriteVariableTree)
    fVariablesTree->Fill();

	if(fAnalCuts->IsSelected(elobj,AliRDHFCuts::kCandidate))
	{
		Double_t cont[3];
		cont[0] = elobj->InvMass(2,pdgdg);
		cont[1] = elobj->Pt();
		cont[2] = fCentrality;
		fHistoEleLambdaMass->Fill(cont);
		Double_t cont2[3];
		cont2[0] = elobj->InvMass(2,pdgdg);
		cont2[1] = trk->Pt();
		cont2[2] = fCentrality;
		Double_t cont_eleptvseta[3];
		cont_eleptvseta[0] = trk->Pt();
		cont_eleptvseta[1] = trk->Eta();
		cont_eleptvseta[2] = fCentrality;

		Double_t cont_eleptvslambdapt[3];
		cont_eleptvslambdapt[0] = trk->Pt();
		cont_eleptvslambdapt[1] = v0->Pt();
		cont_eleptvslambdapt[2] = fCentrality;

		Double_t cont_eleptvsd0[3];
		cont_eleptvsd0[0] = trk->Pt();
		cont_eleptvsd0[1] = elobj->Getd0Prong(0)*trk->Charge();
		cont_eleptvsd0[2] = fCentrality;


		if(mixing_flag){
			if((trk->Charge()>0 && !anti_lambda_flag) || (trk->Charge()<0 && anti_lambda_flag)){
				fHistoEleLambdaMassRSMix->Fill(cont);
				fHistoEleLambdaMassvsElePtRSMix->Fill(cont2);
				if(cont[0]<2.3){
					fHistoElePtRSMix->Fill(trk->Pt(),fCentrality);
					fHistoElePtvsEtaRSMix->Fill(cont_eleptvseta);
					fHistoElePtvsLambdaPtRSMix->Fill(cont_eleptvslambdapt);
					fHistoElePtvsd0RSMix->Fill(cont_eleptvsd0);
				}
			}else if((trk->Charge()<0 && !anti_lambda_flag) || (trk->Charge()>0 && anti_lambda_flag)){
				fHistoEleLambdaMassWSMix->Fill(cont);
				fHistoEleLambdaMassvsElePtWSMix->Fill(cont2);
				if(cont[0]<2.3){
					fHistoElePtWSMix->Fill(trk->Pt(),fCentrality);
					fHistoElePtvsEtaWSMix->Fill(cont_eleptvseta);
					fHistoElePtvsLambdaPtWSMix->Fill(cont_eleptvslambdapt);
					fHistoElePtvsd0WSMix->Fill(cont_eleptvsd0);
				}
			}
		}else{
			if((trk->Charge()>0 && !anti_lambda_flag) || (trk->Charge()<0 && anti_lambda_flag)){
				fHistoEleLambdaMassRS->Fill(cont);
				fHistoEleLambdaMassvsElePtRS->Fill(cont2);
				if(cont[0]<2.3){
					fHistoElePtRS->Fill(trk->Pt(),fCentrality);
					fHistoElePtvsEtaRS->Fill(cont_eleptvseta);
					fHistoElePtvsLambdaPtRS->Fill(cont_eleptvslambdapt);
					fHistoElePtvsd0RS->Fill(cont_eleptvsd0);

					for(Int_t ih=0;ih<17;ih++){
						Double_t cont_eleptvscutvars[3];
						cont_eleptvscutvars[0] = trk->Pt();
						cont_eleptvscutvars[2] = fCentrality;

						if(ih==0){
							cont_eleptvscutvars[1] = trk->GetTPCNcls();
						}else if(ih==1){
							cont_eleptvscutvars[1] = trk->GetTPCsignalN();
						}else if(ih==2){
							cont_eleptvscutvars[1] = nSigmaTPCele;
						}else if(ih==3){
							cont_eleptvscutvars[1] = nSigmaTOFele;
						}else if(ih==4){
							cont_eleptvscutvars[1] = trk->Eta();
						}else if(ih==5){
							cont_eleptvscutvars[1] = trk->GetITSNcls();
						}else if(ih==6){
							if(!anti_lambda_flag)
								cont_eleptvscutvars[1] = v0->MassLambda();
							else
								cont_eleptvscutvars[1] = v0->MassAntiLambda();
						}else if(ih==7){
							Double_t lPosV0[3];
							lPosV0[0] = v0->DecayVertexV0X();
							lPosV0[1] = v0->DecayVertexV0Y();
							lPosV0[2] = v0->DecayVertexV0Z();
							cont_eleptvscutvars[1] = TMath::Sqrt(lPosV0[0]*lPosV0[0]+lPosV0[1]*lPosV0[1]);
						}else if(ih==8){
							cont_eleptvscutvars[1] = v0->DcaV0Daughters();
						}else if(ih==9){
							if(!anti_lambda_flag)
								cont_eleptvscutvars[1] = v0->DcaPosToPrimVertex();
							else
								cont_eleptvscutvars[1] = v0->DcaNegToPrimVertex();
						}else if(ih==10){
							if(!anti_lambda_flag)
								cont_eleptvscutvars[1] = v0->DcaNegToPrimVertex();
							else
								cont_eleptvscutvars[1] = v0->DcaPosToPrimVertex();
						}else if(ih==11){
							cont_eleptvscutvars[1] =  v0->CosPointingAngle(posVtx);
						}else if(ih==12){
							cont_eleptvscutvars[1] =  v0->MassK0Short();
						}else if(ih==13){
							cont_eleptvscutvars[1] =  nSigmaTPCv0pr;
						}else if(ih==14){
							cont_eleptvscutvars[1] =  nSigmaTPCv0pi;
						}else if(ih==15){
							cont_eleptvscutvars[1] =  v0->Eta();
						}else if(ih==16){
							Double_t v0px = elobj->PxProng(1);
							Double_t v0py = elobj->PyProng(1);
							Double_t v0pz = elobj->PzProng(1);
							Double_t epx = elobj->PxProng(0);
							Double_t epy = elobj->PyProng(0);
							Double_t epz = elobj->PzProng(0);
							cont_eleptvscutvars[1] = acos((v0px*epx+v0py*epy+v0pz*epz)/sqrt(v0px*v0px+v0py*v0py+v0pz*v0pz)/sqrt(epx*epx+epy*epy+epz*epz));
						}else{
							cont_eleptvscutvars[1] = -9999.;
						}

						fHistoElePtvsCutVarsRS[ih]->Fill(cont_eleptvscutvars);
					}
				}
			}else if((trk->Charge()<0 && !anti_lambda_flag) || (trk->Charge()>0 && anti_lambda_flag)){
				fHistoEleLambdaMassWS->Fill(cont);
				fHistoEleLambdaMassvsElePtWS->Fill(cont2);
				if(cont[0]<2.3){
					fHistoElePtWS->Fill(trk->Pt(),fCentrality);
					fHistoElePtvsEtaWS->Fill(cont_eleptvseta);
					fHistoElePtvsLambdaPtWS->Fill(cont_eleptvslambdapt);
					fHistoElePtvsd0WS->Fill(cont_eleptvsd0);

					for(Int_t ih=0;ih<17;ih++){
						Double_t cont_eleptvscutvars[3];
						cont_eleptvscutvars[0] = trk->Pt();
						cont_eleptvscutvars[2] = fCentrality;

						if(ih==0){
							cont_eleptvscutvars[1] = trk->GetTPCNcls();
						}else if(ih==1){
							cont_eleptvscutvars[1] = trk->GetTPCsignalN();
						}else if(ih==2){
							cont_eleptvscutvars[1] = nSigmaTPCele;
						}else if(ih==3){
							cont_eleptvscutvars[1] = nSigmaTOFele;
						}else if(ih==4){
							cont_eleptvscutvars[1] = trk->Eta();
						}else if(ih==5){
							cont_eleptvscutvars[1] = trk->GetITSNcls();
						}else if(ih==6){
							if(!anti_lambda_flag)
								cont_eleptvscutvars[1] = v0->MassLambda();
							else
								cont_eleptvscutvars[1] = v0->MassAntiLambda();
						}else if(ih==7){
							Double_t lPosV0[3];
							lPosV0[0] = v0->DecayVertexV0X();
							lPosV0[1] = v0->DecayVertexV0Y();
							lPosV0[2] = v0->DecayVertexV0Z();
							cont_eleptvscutvars[1] = TMath::Sqrt(lPosV0[0]*lPosV0[0]+lPosV0[1]*lPosV0[1]);
						}else if(ih==8){
							cont_eleptvscutvars[1] = v0->DcaV0Daughters();
						}else if(ih==9){
							if(!anti_lambda_flag)
								cont_eleptvscutvars[1] = v0->DcaPosToPrimVertex();
							else
								cont_eleptvscutvars[1] = v0->DcaNegToPrimVertex();
						}else if(ih==10){
							if(!anti_lambda_flag)
								cont_eleptvscutvars[1] = v0->DcaNegToPrimVertex();
							else
								cont_eleptvscutvars[1] = v0->DcaPosToPrimVertex();
						}else if(ih==11){
							cont_eleptvscutvars[1] =  v0->CosPointingAngle(posVtx);
						}else if(ih==12){
							cont_eleptvscutvars[1] =  v0->MassK0Short();
						}else if(ih==13){
							cont_eleptvscutvars[1] =  nSigmaTPCv0pr;
						}else if(ih==14){
							cont_eleptvscutvars[1] =  nSigmaTPCv0pi;
						}else if(ih==15){
							cont_eleptvscutvars[1] =  v0->Eta();
						}else if(ih==16){
							Double_t v0px = elobj->PxProng(1);
							Double_t v0py = elobj->PyProng(1);
							Double_t v0pz = elobj->PzProng(1);
							Double_t epx = elobj->PxProng(0);
							Double_t epy = elobj->PyProng(0);
							Double_t epz = elobj->PzProng(0);
							cont_eleptvscutvars[1] = acos((v0px*epx+v0py*epy+v0pz*epz)/sqrt(v0px*v0px+v0py*v0py+v0pz*v0pz)/sqrt(epx*epx+epy*epy+epz*epz));
						}else{
							cont_eleptvscutvars[1] = -9999.;
						}

						fHistoElePtvsCutVarsWS[ih]->Fill(cont_eleptvscutvars);
					}
				}
			}
		}

		if(fUseMCInfo){
			if(mclc){
				Int_t pdgcode = mclc->GetPdgCode();
				if(abs(pdgcode)==4122 && abs(mcpdgele_array[1])==4122 && abs(mcpdgv0_array[1])==4122){
						fHistoEleLambdaMassMCS->Fill(cont);
						fHistoEleLambdaMassvsElePtMCS->Fill(cont2);
						if(cont[0]<2.3){
							fHistoElePtMCS->Fill(trk->Pt(),fCentrality);
							fHistoElePtvsEtaMCS->Fill(cont_eleptvseta);
							fHistoElePtvsLambdaPtMCS->Fill(cont_eleptvslambdapt);
							fHistoElePtvsd0MCS->Fill(cont_eleptvsd0);

							Double_t cont_eleptvslambdaptvslcpt[4];
							cont_eleptvslambdaptvslcpt[0] = cont_eleptvslambdapt[0];
							cont_eleptvslambdaptvslcpt[1] = cont_eleptvslambdapt[1];
							cont_eleptvslambdaptvslcpt[2] = mclc->Pt();
							cont_eleptvslambdaptvslcpt[3] = cont_eleptvslambdapt[2];
							fHistoElePtvsLambdaPtvsLcPtMCS->Fill(cont_eleptvslambdaptvslcpt);

							Int_t labmotherlc = mclc->GetMother();
							if(labmotherlc>=0){
								AliAODMCParticle *motherlc = (AliAODMCParticle*)mcArray->At(labmotherlc);
								Int_t pdgmotherlc = motherlc->GetPdgCode();
								if(TMath::Abs(pdgmotherlc)==511||TMath::Abs(pdgmotherlc)==521||TMath::Abs(pdgmotherlc)==5122||TMath::Abs(pdgmotherlc)==5132||TMath::Abs(pdgmotherlc)==5232||TMath::Abs(pdgmotherlc)==5332){
									fHistoElePtvsd0BFeeddownMCS->Fill(cont_eleptvsd0);
								}else{
									fHistoElePtvsd0PromptMCS->Fill(cont_eleptvsd0);
								}
							}else{
								fHistoElePtvsd0PromptMCS->Fill(cont_eleptvsd0);
							}

							for(Int_t ih=0;ih<17;ih++){
								Double_t cont_eleptvscutvars[3];
								cont_eleptvscutvars[0] = trk->Pt();
								cont_eleptvscutvars[2] = fCentrality;

								if(ih==0){
									cont_eleptvscutvars[1] = trk->GetTPCNcls();
								}else if(ih==1){
									cont_eleptvscutvars[1] = trk->GetTPCsignalN();
								}else if(ih==2){
									cont_eleptvscutvars[1] = nSigmaTPCele;
								}else if(ih==3){
									cont_eleptvscutvars[1] = nSigmaTOFele;
								}else if(ih==4){
									cont_eleptvscutvars[1] = trk->Eta();
								}else if(ih==5){
									cont_eleptvscutvars[1] = trk->GetITSNcls();
								}else if(ih==6){
									if(!anti_lambda_flag)
										cont_eleptvscutvars[1] = v0->MassLambda();
									else
										cont_eleptvscutvars[1] = v0->MassAntiLambda();
								}else if(ih==7){
									Double_t lPosV0[3];
									lPosV0[0] = v0->DecayVertexV0X();
									lPosV0[1] = v0->DecayVertexV0Y();
									lPosV0[2] = v0->DecayVertexV0Z();
									cont_eleptvscutvars[1] = TMath::Sqrt(lPosV0[0]*lPosV0[0]+lPosV0[1]*lPosV0[1]);
								}else if(ih==8){
									cont_eleptvscutvars[1] = v0->DcaV0Daughters();
								}else if(ih==9){
									if(!anti_lambda_flag)
										cont_eleptvscutvars[1] = v0->DcaPosToPrimVertex();
									else
										cont_eleptvscutvars[1] = v0->DcaNegToPrimVertex();
								}else if(ih==10){
									if(!anti_lambda_flag)
										cont_eleptvscutvars[1] = v0->DcaNegToPrimVertex();
									else
										cont_eleptvscutvars[1] = v0->DcaPosToPrimVertex();
								}else if(ih==11){
									cont_eleptvscutvars[1] =  v0->CosPointingAngle(posVtx);
								}else if(ih==12){
									cont_eleptvscutvars[1] =  v0->MassK0Short();
								}else if(ih==13){
									cont_eleptvscutvars[1] =  nSigmaTPCv0pr;
								}else if(ih==14){
									cont_eleptvscutvars[1] =  nSigmaTPCv0pi;
								}else if(ih==15){
									cont_eleptvscutvars[1] =  v0->Eta();
								}else if(ih==16){
									Double_t v0px = elobj->PxProng(1);
									Double_t v0py = elobj->PyProng(1);
									Double_t v0pz = elobj->PzProng(1);
									Double_t epx = elobj->PxProng(0);
									Double_t epy = elobj->PyProng(0);
									Double_t epz = elobj->PzProng(0);
									cont_eleptvscutvars[1] = acos((v0px*epx+v0py*epy+v0pz*epz)/sqrt(v0px*v0px+v0py*v0py+v0pz*v0pz)/sqrt(epx*epx+epy*epy+epz*epz));
								}else{
									cont_eleptvscutvars[1] = -9999.;
								}

								fHistoElePtvsCutVarsMCS[ih]->Fill(cont_eleptvscutvars);
							}
						}
				}
				if(abs(pdgcode)==4132 && abs(mcpdgele_array[1])==4132 && abs(mcpdgv0_array[1])==3312){
						fHistoEleLambdaMassFeeddownXic0MCS->Fill(cont);
						fHistoEleLambdaMassvsElePtFeeddownXic0MCS->Fill(cont2);
						if(cont[0]<2.3){
							fHistoElePtFeeddownXic0MCS->Fill(trk->Pt(),fCentrality);
							fHistoElePtvsEtaFeeddownXic0MCS->Fill(cont_eleptvseta);
							fHistoElePtvsLambdaPtFeeddownXic0MCS->Fill(cont_eleptvslambdapt);
						}
				}
				if(abs(pdgcode)==4232 && abs(mcpdgele_array[1])==4232 && abs(mcpdgv0_array[1])==3322){
						fHistoEleLambdaMassFeeddownXicPlusMCS->Fill(cont);
						fHistoEleLambdaMassvsElePtFeeddownXicPlusMCS->Fill(cont2);
						if(cont[0]<2.3){
							fHistoElePtFeeddownXicPlusMCS->Fill(trk->Pt(),fCentrality);
							fHistoElePtvsEtaFeeddownXicPlusMCS->Fill(cont_eleptvseta);
							fHistoElePtvsLambdaPtFeeddownXicPlusMCS->Fill(cont_eleptvslambdapt);
						}
				}
			}
		}
	}

  return;
}
////-------------------------------------------------------------------------------
void AliAnalysisTaskSELc2eleLambdafromAODtracks::DefineEleTreeVariables() 
{
  //
  // Define electron tree variables
  //

  const char* nameoutput = GetOutputSlot(5)->GetContainer()->GetName();
  fEleVariablesTree = new TTree(nameoutput,"electron variables tree");
  Int_t nVar = 20;
  fCandidateEleVariables = new Float_t [nVar];
  TString * fCandidateVariableNames = new TString[nVar];

  fCandidateVariableNames[ 0]="ElePx";
  fCandidateVariableNames[ 1]="ElePy";
  fCandidateVariableNames[ 2]="ElePz";
  fCandidateVariableNames[ 3]="TPCChi2overNDF";
  fCandidateVariableNames[ 4]="ITSNcls";
  fCandidateVariableNames[ 5]="TPCNcls";
  fCandidateVariableNames[ 6]="TPCNclsPID";
  fCandidateVariableNames[ 7]="TPCNclsRatio";
  fCandidateVariableNames[ 8]="d0R";
  fCandidateVariableNames[ 9]="d0Z";
  fCandidateVariableNames[10]="ITSClusterMap";
  fCandidateVariableNames[11]="nSigmaTPCele";
  fCandidateVariableNames[12]="nSigmaTOFele";
  fCandidateVariableNames[13]="nSigmaTPCpi";
  fCandidateVariableNames[14]="nSigmaTPCka";
  fCandidateVariableNames[15]="nSigmaTPCpr";
  fCandidateVariableNames[16]="EvNumber";
  fCandidateVariableNames[17]="EleCharge";
  fCandidateVariableNames[18]="Centrality";
  fCandidateVariableNames[19]="RunNumber";

  for (Int_t ivar=0; ivar<nVar; ivar++) {
    fEleVariablesTree->Branch(fCandidateVariableNames[ivar].Data(),&fCandidateEleVariables[ivar],Form("%s/f",fCandidateVariableNames[ivar].Data()));
  }

  return;
}
////-------------------------------------------------------------------------------
void AliAnalysisTaskSELc2eleLambdafromAODtracks::FillElectronROOTObjects(AliAODTrack *trk, TClonesArray *mcArray) 
{
  //
  // Fill histograms or tree depending on fWriteVariableTree 
  //

	if(!trk) return;

	fHistoBachPt->Fill(trk->Pt());
	fHistoElectronQovPtvsPhi->Fill(trk->Phi(),(Double_t)trk->Charge()/trk->Pt());

	if(fUseMCInfo)
	{
		Int_t labEle = trk->GetLabel();
		if(labEle>=0){
			AliAODMCParticle *mcetrk = (AliAODMCParticle*)mcArray->At(labEle);
			if(mcetrk){
				Int_t pdgEle = mcetrk->GetPdgCode();
				if(abs(pdgEle)==11){
					fHistoBachPtMCS->Fill(trk->Pt());
				}
			}
		}
	}

	if(!fWriteEachVariableTree) return;

	for(Int_t i=0;i<20;i++){
		fCandidateEleVariables[i] = -9999.;
	}

  fCandidateEleVariables[ 0] = trk->Px();
  fCandidateEleVariables[ 1] = trk->Py();
  fCandidateEleVariables[ 2] = trk->Pz();
  fCandidateEleVariables[ 3] = trk->Chi2perNDF();
  fCandidateEleVariables[ 4] = trk->GetITSNcls();
  fCandidateEleVariables[ 5] = trk->GetTPCncls();
  fCandidateEleVariables[ 6] = trk->GetTPCsignalN();
	if(trk->GetTPCNclsF()>0) 
		fCandidateEleVariables[ 7] = (Float_t)trk->GetTPCncls()/(Float_t)trk->GetTPCNclsF();

  Double_t d0z0[2],covd0z0[3];
  trk->PropagateToDCA(fVtx1,fBzkG,kVeryBig,d0z0,covd0z0);

  fCandidateEleVariables[ 8] = d0z0[0];
  fCandidateEleVariables[ 9] = d0z0[1];
	Int_t itsmap = trk->GetITSClusterMap();
	Int_t bit1 = 1;
	Int_t bit2 = 2;
	Bool_t spdfirst = (itsmap & bit1) == bit1;
	Bool_t spdsecond = (itsmap & bit2) == bit2;
  fCandidateEleVariables[10] = ((Int_t)spdfirst) + 2 * ((Int_t)spdsecond);

  if(fAnalCuts->GetIsUsePID())
  {
		Double_t nSigmaTPCele = fAnalCuts->GetPidHF()->GetPidResponse()->NumberOfSigmasTPC(trk,AliPID::kElectron);
		Double_t nSigmaTOFele = fAnalCuts->GetPidHF()->GetPidResponse()->NumberOfSigmasTOF(trk,AliPID::kElectron);
		Double_t nSigmaTPCpi_etrk = fAnalCuts->GetPidHF()->GetPidResponse()->NumberOfSigmasTPC(trk,AliPID::kPion);
		Double_t nSigmaTPCka_etrk = fAnalCuts->GetPidHF()->GetPidResponse()->NumberOfSigmasTPC(trk,AliPID::kKaon);
		Double_t nSigmaTPCpr_etrk = fAnalCuts->GetPidHF()->GetPidResponse()->NumberOfSigmasTPC(trk,AliPID::kProton);
    fCandidateEleVariables[11] = nSigmaTPCele;
    fCandidateEleVariables[12] = nSigmaTOFele;
    fCandidateEleVariables[13] = nSigmaTPCpi_etrk;
    fCandidateEleVariables[14] = nSigmaTPCka_etrk;
    fCandidateEleVariables[15] = nSigmaTPCpr_etrk;
  }
  fCandidateEleVariables[16] = fEvNumberCounter;
  fCandidateEleVariables[17] = trk->Charge();
  fCandidateEleVariables[18] = fCentrality;
  fCandidateEleVariables[19] = fRunNumber;

	fHistod0Bach->Fill(d0z0[0]);

	fEleVariablesTree->Fill();
}
////-------------------------------------------------------------------------------
void AliAnalysisTaskSELc2eleLambdafromAODtracks::DefineV0TreeVariables() 
{
  //
  // Define V0 tree variables
  //

  const char* nameoutput = GetOutputSlot(6)->GetContainer()->GetName();
  fV0VariablesTree = new TTree(nameoutput,"v0 variables tree");
  Int_t nVar = 27;
  fCandidateV0Variables = new Float_t [nVar];
  TString * fCandidateVariableNames = new TString[nVar];

  fCandidateVariableNames[ 0]="V0Px";
  fCandidateVariableNames[ 1]="V0Py";
  fCandidateVariableNames[ 2]="V0Pz";
  fCandidateVariableNames[ 3]="MassLambda";
  fCandidateVariableNames[ 4]="MassAntiLambda";
  fCandidateVariableNames[ 5]="ProtonPx";
  fCandidateVariableNames[ 6]="ProtonPy";
  fCandidateVariableNames[ 7]="ProtonPz";
  fCandidateVariableNames[ 8]="PionPx";
  fCandidateVariableNames[ 9]="PionPy";
  fCandidateVariableNames[10]="PionPz";
  fCandidateVariableNames[11]="RfidV0";
  fCandidateVariableNames[12]="DcaV0PrToPrimVertex";
  fCandidateVariableNames[13]="DcaV0PiToPrimVertex";
  fCandidateVariableNames[14]="DcaV0daughters";
  fCandidateVariableNames[15]="V0CosPointingAngle";
  fCandidateVariableNames[16]="V0ProperDecayLength";
  fCandidateVariableNames[17]="MassK0Short";
  fCandidateVariableNames[18]="nSigmaTPCpr";
  fCandidateVariableNames[19]="nSigmaTPCpi";
  fCandidateVariableNames[20]="TPCNCrossV0Pr";
  fCandidateVariableNames[21]="TPCNCrossV0Pi";
  fCandidateVariableNames[22]="TPCNCrossRatioV0Pr";
  fCandidateVariableNames[23]="TPCNCrossRatioV0Pi";
  fCandidateVariableNames[24]="EvNumber";
  fCandidateVariableNames[25]="Centrality";
  fCandidateVariableNames[26]="RunNumber";

  for (Int_t ivar=0; ivar<nVar; ivar++) {
    fV0VariablesTree->Branch(fCandidateVariableNames[ivar].Data(),&fCandidateV0Variables[ivar],Form("%s/f",fCandidateVariableNames[ivar].Data()));
  }

  return;
}

////-------------------------------------------------------------------------------
void AliAnalysisTaskSELc2eleLambdafromAODtracks::FillV0ROOTObjects(AliAODv0 *v0, TClonesArray *mcArray) 
{
  //
  // Fill histograms or tree depending on fWriteVariableTree 
  //
	if(!v0) return;

	if(TMath::Abs(v0->MassLambda()-1.1156)<fAnalCuts->GetProdV0MassTolLambda()){
		fHistoLambdaMassvsPt->Fill(v0->MassLambda(),v0->Pt());
	}
	if(TMath::Abs(v0->MassAntiLambda()-1.1156)<fAnalCuts->GetProdV0MassTolLambda()){
		fHistoLambdaMassvsPt->Fill(v0->MassAntiLambda(),v0->Pt());
	}
	fHistoK0sMassvsPt->Fill(v0->MassK0Short(),v0->Pt());

	Double_t momv0x = v0->MomV0X();
	Double_t momv0y = v0->MomV0Y();
	Double_t phi_alice = atan2(momv0y,momv0x);
	if(phi_alice<0.) phi_alice += 2 * M_PI;
	fHistoLambdaQovPtvsPhi->Fill(phi_alice,1./sqrt(momv0x*momv0x+momv0y*momv0y));

	if(fUseMCInfo)
	{
		Int_t pdgdgv0[2]={2212,211};
		Int_t labV0 = v0->MatchToMC(3122,mcArray,2,pdgdgv0); // the V0
		if(labV0>=0){
			if(TMath::Abs(v0->MassLambda()-1.1156)<0.03){
				fHistoLambdaMassvsPtMCS->Fill(v0->MassLambda(),v0->Pt());
			}
			if(TMath::Abs(v0->MassAntiLambda()-1.1156)<0.03){
				fHistoLambdaMassvsPtMCS->Fill(v0->MassAntiLambda(),v0->Pt());
			}
		}
	}

	if(!fWriteEachVariableTree) return;

	for(Int_t i=0;i<27;i++){
		fCandidateV0Variables[i] = -9999.;
	}

  AliAODTrack *cptrack =  (AliAODTrack*)(v0->GetDaughter(0));
  AliAODTrack *cntrack =  (AliAODTrack*)(v0->GetDaughter(1));
	if(!cptrack) return;
	if(!cntrack) return;
	if(cptrack->Charge()<0 && cntrack->Charge()>0){
		cptrack =  (AliAODTrack*)(v0->GetDaughter(1));
		cntrack =  (AliAODTrack*)(v0->GetDaughter(0));
	}
  Double_t mlamPDG   = TDatabasePDG::Instance()->GetParticle(3122)->Mass();

  fCandidateV0Variables[ 0] = v0->Px();
  fCandidateV0Variables[ 1] = v0->Py();
  fCandidateV0Variables[ 2] = v0->Pz();
  fCandidateV0Variables[ 3] = v0->MassLambda();
  fCandidateV0Variables[ 4] = v0->MassAntiLambda();

	Bool_t isparticle = kTRUE;
	if(fabs(v0->MassAntiLambda()-mlamPDG)<fAnalCuts->GetProdV0MassTolLambda()) isparticle=kFALSE;

	if(isparticle){
		fCandidateV0Variables[ 5] = cptrack->Px();
		fCandidateV0Variables[ 6] = cptrack->Py();
		fCandidateV0Variables[ 7] = cptrack->Pz();
		fCandidateV0Variables[ 8] = cntrack->Px();
		fCandidateV0Variables[ 9] = cntrack->Py();
		fCandidateV0Variables[10] = cntrack->Pz();
	}else{
		fCandidateV0Variables[ 5] = cntrack->Px();
		fCandidateV0Variables[ 6] = cntrack->Py();
		fCandidateV0Variables[ 7] = cntrack->Pz();
		fCandidateV0Variables[ 8] = cptrack->Px();
		fCandidateV0Variables[ 9] = cptrack->Py();
		fCandidateV0Variables[10] = cptrack->Pz();
	}

  Double_t lPosV0[3];
  lPosV0[0] = v0->DecayVertexV0X();
  lPosV0[1] = v0->DecayVertexV0Y();
  lPosV0[2] = v0->DecayVertexV0Z();
  Double_t decayvertV0 = TMath::Sqrt(lPosV0[0]*lPosV0[0]+lPosV0[1]*lPosV0[1]);
	fCandidateV0Variables[11] = decayvertV0;

	Double_t lDcaPosToPrimVertex = v0->DcaPosToPrimVertex();
	Double_t lDcaNegToPrimVertex = v0->DcaNegToPrimVertex();
  if(isparticle){
		fCandidateV0Variables[12] = lDcaPosToPrimVertex;
		fCandidateV0Variables[13] = lDcaNegToPrimVertex;
  }else{
		fCandidateV0Variables[12] = lDcaNegToPrimVertex;
		fCandidateV0Variables[13] = lDcaPosToPrimVertex;
  }
	fCandidateV0Variables[14] = v0->DcaV0Daughters();
  Double_t posVtx[3] = {0.,0.,0.};
  fVtx1->GetXYZ(posVtx);
  fCandidateV0Variables[15] = v0->CosPointingAngle(posVtx); 
  Double_t ptotlam = TMath::Sqrt(pow(v0->Px(),2)+pow(v0->Py(),2)+pow(v0->Pz(),2));
  fCandidateV0Variables[16] = v0->DecayLengthV0(posVtx)*mlamPDG/ptotlam;
  fCandidateV0Variables[17] = v0->MassK0Short();

  if(fAnalCuts->GetUseLambdaPID())
  {
		if(isparticle){
			Double_t nSigmaTPCv0pr = fAnalCuts->GetPidProton()->GetPidResponse()->NumberOfSigmasTPC(cptrack,AliPID::kProton);
			Double_t nSigmaTPCv0pi = fAnalCuts->GetPidPion()->GetPidResponse()->NumberOfSigmasTPC(cntrack,AliPID::kPion);
			fCandidateV0Variables[18] = nSigmaTPCv0pr;
			fCandidateV0Variables[19] = nSigmaTPCv0pi;
		}else{
			Double_t nSigmaTPCv0pr = fAnalCuts->GetPidProton()->GetPidResponse()->NumberOfSigmasTPC(cntrack,AliPID::kProton);
			Double_t nSigmaTPCv0pi = fAnalCuts->GetPidPion()->GetPidResponse()->NumberOfSigmasTPC(cptrack,AliPID::kPion);
			fCandidateV0Variables[18] = nSigmaTPCv0pr;
			fCandidateV0Variables[19] = nSigmaTPCv0pi;
		}
  }
	if(isparticle){
		fCandidateV0Variables[20] = cptrack->GetTPCClusterInfo(2,1);
		fCandidateV0Variables[21] = cntrack->GetTPCClusterInfo(2,1);
		if(cptrack->GetTPCNclsF()>0)
			fCandidateV0Variables[22] = (Float_t) cptrack->GetTPCClusterInfo(2,1)/(Float_t)cptrack->GetTPCNclsF();
		if(cntrack->GetTPCNclsF()>0)
			fCandidateV0Variables[23] =(Float_t)  cntrack->GetTPCClusterInfo(2,1)/(Float_t)cntrack->GetTPCNclsF();
	}else{
		fCandidateV0Variables[20] = cntrack->GetTPCClusterInfo(2,1);
		fCandidateV0Variables[21] = cptrack->GetTPCClusterInfo(2,1);
		if(cntrack->GetTPCNclsF()>0)
			fCandidateV0Variables[22] = (Float_t) cntrack->GetTPCClusterInfo(2,1)/(Float_t)cntrack->GetTPCNclsF();
		if(cptrack->GetTPCNclsF()>0)
			fCandidateV0Variables[23] = (Float_t) cptrack->GetTPCClusterInfo(2,1)/(Float_t)cptrack->GetTPCNclsF();
	}
	fCandidateV0Variables[24] = fEvNumberCounter;
	fCandidateV0Variables[25] = fCentrality;
	fCandidateV0Variables[26] = fRunNumber;


		fV0VariablesTree->Fill();
}
////-------------------------------------------------------------------------------
void AliAnalysisTaskSELc2eleLambdafromAODtracks::DefineMCTreeVariables() 
{
  //
  // Define electron tree variables
  //

  const char* nameoutput = GetOutputSlot(7)->GetContainer()->GetName();
  fMCVariablesTree = new TTree(nameoutput,"MC variables tree");
  Int_t nVar = 11;
  fCandidateMCVariables = new Float_t [nVar];
  TString * fCandidateVariableNames = new TString[nVar];

  fCandidateVariableNames[ 0]="Centrality";
  fCandidateVariableNames[ 1]="DecayType";
  fCandidateVariableNames[ 2]="LcPx";
  fCandidateVariableNames[ 3]="LcPy";
  fCandidateVariableNames[ 4]="LcPz";
  fCandidateVariableNames[ 5]="ElePx";
  fCandidateVariableNames[ 6]="ElePy";
  fCandidateVariableNames[ 7]="ElePz";
  fCandidateVariableNames[ 8]="V0Px";
  fCandidateVariableNames[ 9]="V0Py";
  fCandidateVariableNames[10]="V0Pz";

  for (Int_t ivar=0; ivar<nVar; ivar++) {
    fMCVariablesTree->Branch(fCandidateVariableNames[ivar].Data(),&fCandidateMCVariables[ivar],Form("%s/f",fCandidateVariableNames[ivar].Data()));
  }
  return;
}
////-------------------------------------------------------------------------------
void AliAnalysisTaskSELc2eleLambdafromAODtracks::FillMCROOTObjects(AliAODMCParticle *mcpart, AliAODMCParticle *mcepart, AliAODMCParticle *mcv0part, Int_t decaytype) 
{
  //
  // Fill histograms or tree depending on fWriteMCVariableTree 
  //
	if(!mcpart) return;
	if(!mcepart) return;
	if(!mcv0part) return;

	for(Int_t i=0;i<11;i++){
		fCandidateMCVariables[i] = -9999.;
	}

	fCandidateMCVariables[ 0] = fCentrality;
	fCandidateMCVariables[ 1] = (Float_t) decaytype;
	fCandidateMCVariables[ 2] = mcpart->Px();
	fCandidateMCVariables[ 3] = mcpart->Py();
	fCandidateMCVariables[ 4] = mcpart->Pz();
	fCandidateMCVariables[ 5] = mcepart->Px();
	fCandidateMCVariables[ 6] = mcepart->Py();
	fCandidateMCVariables[ 7] = mcepart->Pz();
	fCandidateMCVariables[ 8] = mcv0part->Px();
	fCandidateMCVariables[ 9] = mcv0part->Py();
	fCandidateMCVariables[10] = mcv0part->Pz();

	Double_t epx = mcepart->Px();
	Double_t epy = mcepart->Py();
	Double_t epz = mcepart->Pz();
	Double_t eE = sqrt(epx*epx+epy*epy+epz*epz+0.000511*0.000511);
	Double_t v0px = mcv0part->Px();
	Double_t v0py = mcv0part->Py();
	Double_t v0pz = mcv0part->Pz();
	Double_t v0E = sqrt(v0px*v0px+v0py*v0py+v0pz*v0pz+1.1157*1.1157);

	Double_t InvMassEleLambda = sqrt(pow(eE+v0E,2)-pow(epx+v0px,2)-pow(epy+v0py,2)-pow(epz+v0pz,2));

	Double_t cont[3];
	cont[0] = InvMassEleLambda;
	cont[1] = mcpart->Pt();
	cont[2] = fCentrality;
	Double_t cont2[3];
	cont2[0] = InvMassEleLambda;
	cont2[1] = mcepart->Pt();
	cont2[2] = fCentrality;
	Double_t cont_eleptvseta[3];
	cont_eleptvseta[0] = mcepart->Pt();
	cont_eleptvseta[1] = mcepart->Eta();
	cont_eleptvseta[2] = fCentrality;
	Double_t cont_eleptvslambdapt[3];
	cont_eleptvslambdapt[0] = mcepart->Pt();
	cont_eleptvslambdapt[1] = mcv0part->Pt();
	cont_eleptvslambdapt[2] = fCentrality;
	Double_t cont_eleptvslambdaptvslcpt[4];
	cont_eleptvslambdaptvslcpt[0] = mcepart->Pt();
	cont_eleptvslambdaptvslcpt[1] = mcv0part->Pt();
	cont_eleptvslambdaptvslcpt[2] = mcpart->Pt();
	cont_eleptvslambdaptvslcpt[3] = fCentrality;

	AliESDtrackCuts *esdcuts = fAnalCuts->GetTrackCuts();
	Float_t etamin, etamax;
	esdcuts->GetEtaRange(etamin,etamax);

	if(decaytype==0){
		fHistoEleLambdaMassMCGen->Fill(cont);
		if(fabs(mcepart->Eta())<etamax){
			fHistoEleLambdaMassvsElePtMCGen->Fill(cont2);
			if(InvMassEleLambda<2.3){
				fHistoElePtMCGen->Fill(mcepart->Pt(),fCentrality);
				fHistoElePtvsEtaMCGen->Fill(cont_eleptvseta);
				fHistoElePtvsLambdaPtMCGen->Fill(cont_eleptvslambdapt);
			}
		}
		if(fabs(mcpart->Y())<0.7){
			if(InvMassEleLambda<2.3){
				fHistoElePtvsLambdaPtMCLcGen->Fill(cont_eleptvslambdapt);
				fHistoElePtvsLambdaPtvsLcPtMCGen->Fill(cont_eleptvslambdaptvslcpt);
			}
		}
	}else if(decaytype==1){
		fHistoEleLambdaMassFeeddownXic0MCGen->Fill(cont);
		if(fabs(mcepart->Eta())<etamax){
			fHistoEleLambdaMassvsElePtFeeddownXic0MCGen->Fill(cont2);
			if(InvMassEleLambda<2.3){
				fHistoElePtFeeddownXic0MCGen->Fill(mcepart->Pt(),fCentrality);
				fHistoElePtvsEtaFeeddownXic0MCGen->Fill(cont_eleptvseta);
				fHistoElePtvsLambdaPtFeeddownXic0MCGen->Fill(cont_eleptvslambdapt);
			}
		}
	}else if(decaytype==2){
		fHistoEleLambdaMassFeeddownXicPlusMCGen->Fill(cont);
		if(fabs(mcepart->Eta())<etamax){
			fHistoEleLambdaMassvsElePtFeeddownXicPlusMCGen->Fill(cont2);
			if(InvMassEleLambda<2.3){
				fHistoElePtFeeddownXicPlusMCGen->Fill(mcepart->Pt(),fCentrality);
				fHistoElePtvsEtaFeeddownXicPlusMCGen->Fill(cont_eleptvseta);
				fHistoElePtvsLambdaPtFeeddownXicPlusMCGen->Fill(cont_eleptvslambdapt);
			}
		}
	}

	if(fWriteMCVariableTree)
		fMCVariablesTree->Fill();
}



////__________________________________________________________________________
void  AliAnalysisTaskSELc2eleLambdafromAODtracks::DefineGeneralHistograms() {
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


  fOutput->Add(fCEvents);
  fOutput->Add(fHTrigger);
  fOutput->Add(fHCentrality);

  return;
}
//__________________________________________________________________________
void  AliAnalysisTaskSELc2eleLambdafromAODtracks::DefineAnalysisHistograms() 
{
  //
  // Define analyis histograms
  //
	
  //------------------------------------------------
  // Basic histogram
  //------------------------------------------------
  Int_t bins_base[3]=		{16,100		,10};
  Double_t xmin_base[3]={1.1,0		,0.00};
  Double_t xmax_base[3]={3.1,10.	,100};
  fHistoEleLambdaMass = new THnSparseF("fHistoEleLambdaMass","",3,bins_base,xmin_base,xmax_base);
  fOutputAll->Add(fHistoEleLambdaMass);
  fHistoEleLambdaMassRS = new THnSparseF("fHistoEleLambdaMassRS","",3,bins_base,xmin_base,xmax_base);
  fOutputAll->Add(fHistoEleLambdaMassRS);
  fHistoEleLambdaMassWS = new THnSparseF("fHistoEleLambdaMassWS","",3,bins_base,xmin_base,xmax_base);
  fOutputAll->Add(fHistoEleLambdaMassWS);
  fHistoEleLambdaMassRSMix = new THnSparseF("fHistoEleLambdaMassRSMix","",3,bins_base,xmin_base,xmax_base);
  fOutputAll->Add(fHistoEleLambdaMassRSMix);
  fHistoEleLambdaMassWSMix = new THnSparseF("fHistoEleLambdaMassWSMix","",3,bins_base,xmin_base,xmax_base);
  fOutputAll->Add(fHistoEleLambdaMassWSMix);
  fHistoEleLambdaMassvsElePtRS = new THnSparseF("fHistoEleLambdaMassvsElePtRS","",3,bins_base,xmin_base,xmax_base);
  fOutputAll->Add(fHistoEleLambdaMassvsElePtRS);
  fHistoEleLambdaMassvsElePtWS = new THnSparseF("fHistoEleLambdaMassvsElePtWS","",3,bins_base,xmin_base,xmax_base);
  fOutputAll->Add(fHistoEleLambdaMassvsElePtWS);
  fHistoEleLambdaMassvsElePtRSMix = new THnSparseF("fHistoEleLambdaMassvsElePtRSMix","",3,bins_base,xmin_base,xmax_base);
  fOutputAll->Add(fHistoEleLambdaMassvsElePtRSMix);
  fHistoEleLambdaMassvsElePtWSMix = new THnSparseF("fHistoEleLambdaMassvsElePtWSMix","",3,bins_base,xmin_base,xmax_base);
  fOutputAll->Add(fHistoEleLambdaMassvsElePtWSMix);
  fHistoEleLambdaMassMCS = new THnSparseF("fHistoEleLambdaMassMCS","",3,bins_base,xmin_base,xmax_base);
  fOutputAll->Add(fHistoEleLambdaMassMCS);
  fHistoEleLambdaMassMCGen = new THnSparseF("fHistoEleLambdaMassMCGen","",3,bins_base,xmin_base,xmax_base);
  fOutputAll->Add(fHistoEleLambdaMassMCGen);
  fHistoEleLambdaMassvsElePtMCS = new THnSparseF("fHistoEleLambdaMassvsElePtMCS","",3,bins_base,xmin_base,xmax_base);
  fOutputAll->Add(fHistoEleLambdaMassvsElePtMCS);
  fHistoEleLambdaMassvsElePtMCGen = new THnSparseF("fHistoEleLambdaMassvsElePtMCGen","",3,bins_base,xmin_base,xmax_base);
  fOutputAll->Add(fHistoEleLambdaMassvsElePtMCGen);

  fHistoElePtRS = new TH2F("fHistoElePtRS","",100,0,10,10,0,100);
  fOutputAll->Add(fHistoElePtRS);
  fHistoElePtWS = new TH2F("fHistoElePtWS","",100,0.,10.,10,0,100);
  fOutputAll->Add(fHistoElePtWS);
  fHistoElePtRSMix = new TH2F("fHistoElePtRSMix","",100,0.,10.,10,0,100);
  fOutputAll->Add(fHistoElePtRSMix);
  fHistoElePtWSMix = new TH2F("fHistoElePtWSMix","",100,0.,10.,10,0,100);
  fOutputAll->Add(fHistoElePtWSMix);
  fHistoElePtMCS = new TH2F("fHistoElePtMCS","",100,0,10,10,0,100);
  fOutputAll->Add(fHistoElePtMCS);
  fHistoElePtMCGen = new TH2F("fHistoElePtMCGen","",100,0,10,10,0,100);
  fOutputAll->Add(fHistoElePtMCGen);

  Int_t bins_eleptvseta[3]=		{50,20	,10};
  Double_t xmin_eleptvseta[3]={0.,-1.	,0.0};
  Double_t xmax_eleptvseta[3]={5.,1.	,100};

  fHistoElePtvsEtaRS = new THnSparseF("fHistoElePtvsEtaRS","",3,bins_eleptvseta,xmin_eleptvseta,xmax_eleptvseta);
  fOutputAll->Add(fHistoElePtvsEtaRS);
  fHistoElePtvsEtaWS = new THnSparseF("fHistoElePtvsEtaWS","",3,bins_eleptvseta,xmin_eleptvseta,xmax_eleptvseta);
  fOutputAll->Add(fHistoElePtvsEtaWS);
  fHistoElePtvsEtaRSMix = new THnSparseF("fHistoElePtvsEtaRSMix","",3,bins_eleptvseta,xmin_eleptvseta,xmax_eleptvseta);
  fOutputAll->Add(fHistoElePtvsEtaRSMix);
  fHistoElePtvsEtaWSMix = new THnSparseF("fHistoElePtvsEtaWSMix","",3,bins_eleptvseta,xmin_eleptvseta,xmax_eleptvseta);
  fOutputAll->Add(fHistoElePtvsEtaWSMix);
  fHistoElePtvsEtaMCS = new THnSparseF("fHistoElePtvsEtaMCS","",3,bins_eleptvseta,xmin_eleptvseta,xmax_eleptvseta);
  fOutputAll->Add(fHistoElePtvsEtaMCS);
  fHistoElePtvsEtaMCGen = new THnSparseF("fHistoElePtvsEtaMCGen","",3,bins_eleptvseta,xmin_eleptvseta,xmax_eleptvseta);
  fOutputAll->Add(fHistoElePtvsEtaMCGen);

  Int_t bins_eleptvslambdapt[3]=	{50,20	,10};
  Double_t xmin_eleptvslambdapt[3]={0.,0.	,0.0};
  Double_t xmax_eleptvslambdapt[3]={5.,5.	,100};

  fHistoElePtvsLambdaPtRS = new THnSparseF("fHistoElePtvsLambdaPtRS","",3,bins_eleptvslambdapt,xmin_eleptvslambdapt,xmax_eleptvslambdapt);
  fOutputAll->Add(fHistoElePtvsLambdaPtRS);
  fHistoElePtvsLambdaPtWS = new THnSparseF("fHistoElePtvsLambdaPtWS","",3,bins_eleptvslambdapt,xmin_eleptvslambdapt,xmax_eleptvslambdapt);
  fOutputAll->Add(fHistoElePtvsLambdaPtWS);
  fHistoElePtvsLambdaPtRSMix = new THnSparseF("fHistoElePtvsLambdaPtRSMix","",3,bins_eleptvslambdapt,xmin_eleptvslambdapt,xmax_eleptvslambdapt);
  fOutputAll->Add(fHistoElePtvsLambdaPtRSMix);
  fHistoElePtvsLambdaPtWSMix = new THnSparseF("fHistoElePtvsLambdaPtWSMix","",3,bins_eleptvslambdapt,xmin_eleptvslambdapt,xmax_eleptvslambdapt);
  fOutputAll->Add(fHistoElePtvsLambdaPtWSMix);
  fHistoElePtvsLambdaPtMCS = new THnSparseF("fHistoElePtvsLambdaPtMCS","",3,bins_eleptvslambdapt,xmin_eleptvslambdapt,xmax_eleptvslambdapt);
  fOutputAll->Add(fHistoElePtvsLambdaPtMCS);
  fHistoElePtvsLambdaPtMCGen = new THnSparseF("fHistoElePtvsLambdaPtMCGen","",3,bins_eleptvslambdapt,xmin_eleptvslambdapt,xmax_eleptvslambdapt);
  fOutputAll->Add(fHistoElePtvsLambdaPtMCGen);
  fHistoElePtvsLambdaPtMCLcGen = new THnSparseF("fHistoElePtvsLambdaPtMCLcGen","",3,bins_eleptvslambdapt,xmin_eleptvslambdapt,xmax_eleptvslambdapt);
  fOutputAll->Add(fHistoElePtvsLambdaPtMCLcGen);

  Int_t bins_eleptvslambdaptvslcpt[4]=	{50,20,10,10};
  Double_t xmin_eleptvslambdaptvslcpt[4]={0.,0.,0.,0.0};
  Double_t xmax_eleptvslambdaptvslcpt[4]={5.,5.,10.,100};
  fHistoElePtvsLambdaPtvsLcPtMCS = new THnSparseF("fHistoElePtvsLambdaPtvsLcPtMCS","",4,bins_eleptvslambdaptvslcpt,xmin_eleptvslambdaptvslcpt,xmax_eleptvslambdaptvslcpt);
  fOutputAll->Add(fHistoElePtvsLambdaPtvsLcPtMCS);
  fHistoElePtvsLambdaPtvsLcPtMCGen = new THnSparseF("fHistoElePtvsLambdaPtvsLcPtMCGen","",4,bins_eleptvslambdaptvslcpt,xmin_eleptvslambdaptvslcpt,xmax_eleptvslambdaptvslcpt);
  fOutputAll->Add(fHistoElePtvsLambdaPtvsLcPtMCGen);

  Int_t bins_eleptvsd0[3]=	{50 ,50	,10};
  Double_t xmin_eleptvsd0[3]={0.,-0.2	,0.0};
  Double_t xmax_eleptvsd0[3]={5.,0.2	,100};

  fHistoElePtvsd0RS = new THnSparseF("fHistoElePtvsd0RS","",3,bins_eleptvsd0,xmin_eleptvsd0,xmax_eleptvsd0);
  fOutputAll->Add(fHistoElePtvsd0RS);
  fHistoElePtvsd0WS = new THnSparseF("fHistoElePtvsd0WS","",3,bins_eleptvsd0,xmin_eleptvsd0,xmax_eleptvsd0);
  fOutputAll->Add(fHistoElePtvsd0WS);
  fHistoElePtvsd0RSMix = new THnSparseF("fHistoElePtvsd0RSMix","",3,bins_eleptvsd0,xmin_eleptvsd0,xmax_eleptvsd0);
  fOutputAll->Add(fHistoElePtvsd0RSMix);
  fHistoElePtvsd0WSMix = new THnSparseF("fHistoElePtvsd0WSMix","",3,bins_eleptvsd0,xmin_eleptvsd0,xmax_eleptvsd0);
  fOutputAll->Add(fHistoElePtvsd0WSMix);
  fHistoElePtvsd0MCS = new THnSparseF("fHistoElePtvsd0MCS","",3,bins_eleptvsd0,xmin_eleptvsd0,xmax_eleptvsd0);
  fOutputAll->Add(fHistoElePtvsd0MCS);
  fHistoElePtvsd0PromptMCS = new THnSparseF("fHistoElePtvsd0PromptMCS","",3,bins_eleptvsd0,xmin_eleptvsd0,xmax_eleptvsd0);
  fOutputAll->Add(fHistoElePtvsd0PromptMCS);
  fHistoElePtvsd0BFeeddownMCS = new THnSparseF("fHistoElePtvsd0BFeeddownMCS","",3,bins_eleptvsd0,xmin_eleptvsd0,xmax_eleptvsd0);
  fOutputAll->Add(fHistoElePtvsd0BFeeddownMCS);


	//Feeddown from Xic0
  fHistoEleLambdaMassFeeddownXic0MCS = new THnSparseF("fHistoEleLambdaMassFeeddownXic0MCS","",3,bins_base,xmin_base,xmax_base);
  fOutputAll->Add(fHistoEleLambdaMassFeeddownXic0MCS);
  fHistoEleLambdaMassFeeddownXic0MCGen = new THnSparseF("fHistoEleLambdaMassFeeddownXic0MCGen","",3,bins_base,xmin_base,xmax_base);
  fOutputAll->Add(fHistoEleLambdaMassFeeddownXic0MCGen);
  fHistoEleLambdaMassvsElePtFeeddownXic0MCS = new THnSparseF("fHistoEleLambdaMassvsElePtFeeddownXic0MCS","",3,bins_base,xmin_base,xmax_base);
  fOutputAll->Add(fHistoEleLambdaMassvsElePtFeeddownXic0MCS);
  fHistoEleLambdaMassvsElePtFeeddownXic0MCGen = new THnSparseF("fHistoEleLambdaMassvsElePtFeeddownXic0MCGen","",3,bins_base,xmin_base,xmax_base);
  fOutputAll->Add(fHistoEleLambdaMassvsElePtFeeddownXic0MCGen);
  fHistoElePtFeeddownXic0MCS = new TH2F("fHistoElePtFeeddownXic0MCS","",100,0,10,10,0,100);
  fOutputAll->Add(fHistoElePtFeeddownXic0MCS);
  fHistoElePtFeeddownXic0MCGen = new TH2F("fHistoElePtFeeddownXic0MCGen","",100,0,10,10,0,100);
  fOutputAll->Add(fHistoElePtFeeddownXic0MCGen);
  fHistoElePtvsEtaFeeddownXic0MCS = new THnSparseF("fHistoElePtvsEtaFeeddownXic0MCS","",3,bins_eleptvseta,xmin_eleptvseta,xmax_eleptvseta);
  fOutputAll->Add(fHistoElePtvsEtaFeeddownXic0MCS);
  fHistoElePtvsEtaFeeddownXic0MCGen = new THnSparseF("fHistoElePtvsEtaFeeddownXic0MCGen","",3,bins_eleptvseta,xmin_eleptvseta,xmax_eleptvseta);
  fOutputAll->Add(fHistoElePtvsEtaFeeddownXic0MCGen);
  fHistoElePtvsLambdaPtFeeddownXic0MCS = new THnSparseF("fHistoElePtvsLambdaPtFeeddownXic0MCS","",3,bins_eleptvslambdapt,xmin_eleptvslambdapt,xmax_eleptvslambdapt);
  fOutputAll->Add(fHistoElePtvsLambdaPtFeeddownXic0MCS);
  fHistoElePtvsLambdaPtFeeddownXic0MCGen = new THnSparseF("fHistoElePtvsLambdaPtFeeddownXic0MCGen","",3,bins_eleptvslambdapt,xmin_eleptvslambdapt,xmax_eleptvslambdapt);
  fOutputAll->Add(fHistoElePtvsLambdaPtFeeddownXic0MCGen);

	//Feeddown from XicPlus
  fHistoEleLambdaMassFeeddownXicPlusMCS = new THnSparseF("fHistoEleLambdaMassFeeddownXicPlusMCS","",3,bins_base,xmin_base,xmax_base);
  fOutputAll->Add(fHistoEleLambdaMassFeeddownXicPlusMCS);
  fHistoEleLambdaMassFeeddownXicPlusMCGen = new THnSparseF("fHistoEleLambdaMassFeeddownXicPlusMCGen","",3,bins_base,xmin_base,xmax_base);
  fOutputAll->Add(fHistoEleLambdaMassFeeddownXicPlusMCGen);
  fHistoEleLambdaMassvsElePtFeeddownXicPlusMCS = new THnSparseF("fHistoEleLambdaMassvsElePtFeeddownXicPlusMCS","",3,bins_base,xmin_base,xmax_base);
  fOutputAll->Add(fHistoEleLambdaMassvsElePtFeeddownXicPlusMCS);
  fHistoEleLambdaMassvsElePtFeeddownXicPlusMCGen = new THnSparseF("fHistoEleLambdaMassvsElePtFeeddownXicPlusMCGen","",3,bins_base,xmin_base,xmax_base);
  fOutputAll->Add(fHistoEleLambdaMassvsElePtFeeddownXicPlusMCGen);
  fHistoElePtFeeddownXicPlusMCS = new TH2F("fHistoElePtFeeddownXicPlusMCS","",100,0,10,10,0,100);
  fOutputAll->Add(fHistoElePtFeeddownXicPlusMCS);
  fHistoElePtFeeddownXicPlusMCGen = new TH2F("fHistoElePtFeeddownXicPlusMCGen","",100,0,10,10,0,100);
  fOutputAll->Add(fHistoElePtFeeddownXicPlusMCGen);
  fHistoElePtvsEtaFeeddownXicPlusMCS = new THnSparseF("fHistoElePtvsEtaFeeddownXicPlusMCS","",3,bins_eleptvseta,xmin_eleptvseta,xmax_eleptvseta);
  fOutputAll->Add(fHistoElePtvsEtaFeeddownXicPlusMCS);
  fHistoElePtvsEtaFeeddownXicPlusMCGen = new THnSparseF("fHistoElePtvsEtaFeeddownXicPlusMCGen","",3,bins_eleptvseta,xmin_eleptvseta,xmax_eleptvseta);
  fOutputAll->Add(fHistoElePtvsEtaFeeddownXicPlusMCGen);
  fHistoElePtvsLambdaPtFeeddownXicPlusMCS = new THnSparseF("fHistoElePtvsLambdaPtFeeddownXicPlusMCS","",3,bins_eleptvslambdapt,xmin_eleptvslambdapt,xmax_eleptvslambdapt);
  fOutputAll->Add(fHistoElePtvsLambdaPtFeeddownXicPlusMCS);
  fHistoElePtvsLambdaPtFeeddownXicPlusMCGen = new THnSparseF("fHistoElePtvsLambdaPtFeeddownXicPlusMCGen","",3,bins_eleptvslambdapt,xmin_eleptvslambdapt,xmax_eleptvslambdapt);
  fOutputAll->Add(fHistoElePtvsLambdaPtFeeddownXicPlusMCGen);

  //------------------------------------------------
  // checking histograms
  //------------------------------------------------
  fHistoBachPt = new TH1F("fHistoBachPt","Bachelor p_{T}",100,0.0,5.0);
  fOutputAll->Add(fHistoBachPt);
  fHistoBachPtMCS = new TH1F("fHistoBachPtMCS","Bachelor p_{T}",100,0.0,5.0);
  fOutputAll->Add(fHistoBachPtMCS);
  fHistoBachPtMCGen = new TH1F("fHistoBachPtMCGen","Bachelor p_{T}",100,0.0,5.0);
  fOutputAll->Add(fHistoBachPtMCGen);
  fHistod0Bach = new TH1F("fHistod0Bach","Bachelor d_{0}",100,-0.5,0.5);
  fOutputAll->Add(fHistod0Bach);
  fHistoLambdaMassvsPt=new TH2F("fHistoLambdaMassvsPt","Lambda mass",100,1.116-0.05,1.116+0.05,20,0.,10.);
  fOutputAll->Add(fHistoLambdaMassvsPt);
  fHistoLambdaMassvsPtMCS=new TH2F("fHistoLambdaMassvsPtMCS","Lambda mass",100,1.116-0.05,1.116+0.05,20,0.,10.);
  fOutputAll->Add(fHistoLambdaMassvsPtMCS);
  fHistoLambdaMassvsPtMCGen=new TH2F("fHistoLambdaMassvsPtMCGen","Lambda mass",100,1.116-0.05,1.116+0.05,20,0.,10.);
  fOutputAll->Add(fHistoLambdaMassvsPtMCGen);
  fHistoK0sMassvsPt=new TH2F("fHistoK0sMassvsPt","K0s mass",100,0.497-0.05,0.497+0.05,20,0.,10.);
  fOutputAll->Add(fHistoK0sMassvsPt);

  fHistoElectronTPCPID=new TH2F("fHistoElectronTPCPID","",50,0.,5.,50,-20.,20.);
  fOutputAll->Add(fHistoElectronTPCPID);
  fHistoElectronTOFPID=new TH2F("fHistoElectronTOFPID","",50,0.,5.,50,-20.,20.);
  fOutputAll->Add(fHistoElectronTOFPID);
  fHistoElectronTPCSelPID=new TH2F("fHistoElectronTPCSelPID","",50,0.,5.,50,-20.,20.);
  fOutputAll->Add(fHistoElectronTPCSelPID);
  fHistoElectronTOFSelPID=new TH2F("fHistoElectronTOFSelPID","",50,0.,5.,50,-20.,20.);
  fOutputAll->Add(fHistoElectronTOFSelPID);
  fHistoElectronTPCPIDSelTOF=new TH2F("fHistoElectronTPCPIDSelTOF","",10,0.,5.,500,-10.,10.);
  fOutputAll->Add(fHistoElectronTPCPIDSelTOF);
  fHistoElectronTPCPIDSelTOFSmallEta=new TH2F("fHistoElectronTPCPIDSelTOFSmallEta","",10,0.,5.,500,-10.,10.);
  fOutputAll->Add(fHistoElectronTPCPIDSelTOFSmallEta);
  fHistoElectronTPCPIDSelTOFLargeEta=new TH2F("fHistoElectronTPCPIDSelTOFLargeEta","",10,0.,5.,500,-10.,10.);
  fOutputAll->Add(fHistoElectronTPCPIDSelTOFLargeEta);

	for(Int_t i=0;i<8;i++){
		fHistoElectronTPCPIDSelTOFEtaDep[i]=new TH2F(Form("fHistoElectronTPCPIDSelTOFEtaDep[%d]",i),"",10,0.,5.,500,-10.,10.);
		fOutputAll->Add(fHistoElectronTPCPIDSelTOFEtaDep[i]);
	}
  fHistoElectronQovPtvsPhi=new TH2F("fHistoElectronQovPtvsPhi","",70,0.,7.,50,-2.,2.);
  fOutputAll->Add(fHistoElectronQovPtvsPhi);
  fHistoLambdaQovPtvsPhi=new TH2F("fHistoLambdaQovPtvsPhi","",70,0.,7.,50,-2.,2.);
  fOutputAll->Add(fHistoLambdaQovPtvsPhi);

  fHistonEvtvsRunNumber=new TH1F("fHistonEvtvsRunNumber","",20000,-0.5,19999.5);
  fOutputAll->Add(fHistonEvtvsRunNumber);
  fHistonElevsRunNumber=new TH1F("fHistonElevsRunNumber","",20000,-0.5,19999.5);
  fOutputAll->Add(fHistonElevsRunNumber);
  fHistonLambdavsRunNumber=new TH1F("fHistonLambdavsRunNumber","",20000,-0.5,19999.5);
  fOutputAll->Add(fHistonLambdavsRunNumber);

	for(Int_t ih=0;ih<17;ih++){
		Int_t bins_eleptvscutvars[3];
		Double_t xmin_eleptvscutvars[3];
		Double_t xmax_eleptvscutvars[3];

		bins_eleptvscutvars[0] = 50;//electron pT bin
		xmin_eleptvscutvars[0] = 0.;
		xmax_eleptvscutvars[0] = 5.;
		bins_eleptvscutvars[2] = 10;//centrality bin
		xmin_eleptvscutvars[2] = 0.;
		xmax_eleptvscutvars[2] = 100.;

		if(ih==0 || ih==1){
			//0: TPC Ncluster 1: TPC ncluster PID
			bins_eleptvscutvars[1] = 40;
			xmin_eleptvscutvars[1] = 0.;
			xmax_eleptvscutvars[1] = 160.;
		}else if(ih==2 || ih==3){
			//2: nSigma(TPC,e) 3: nSigma(TOF,e)
			bins_eleptvscutvars[1] = 20;
			xmin_eleptvscutvars[1] = -5.;
			xmax_eleptvscutvars[1] = 5.;
		}else if(ih==4){
			//4: eta
			bins_eleptvscutvars[1] = 30;
			xmin_eleptvscutvars[1] = -1.5;
			xmax_eleptvscutvars[1] = 1.5;
		}else if(ih==5){
			//5: nITS cluster
			bins_eleptvscutvars[1] = 7;
			xmin_eleptvscutvars[1] = -0.5;
			xmax_eleptvscutvars[1] = 6.5;
		}else if(ih==6){
			//6: Lambda mass
			bins_eleptvscutvars[1] = 50;
			xmin_eleptvscutvars[1] = 1.1156-0.03;
			xmax_eleptvscutvars[1] = 1.1156+0.03;
		}else if(ih==7){
			//7: Rfid Lambda
			bins_eleptvscutvars[1] = 20;
			xmin_eleptvscutvars[1] = 0.;
			xmax_eleptvscutvars[1] = 5.;
		}else if(ih==8){
			//10: Dca V0
			bins_eleptvscutvars[1] = 20;
			xmin_eleptvscutvars[1] = 0.;
			xmax_eleptvscutvars[1] = 2.;
		}else if(ih==9 || ih==10 ){
			//9: DCA V0pr to prim 10: DCA V0pi to prim
			bins_eleptvscutvars[1] = 20;
			xmin_eleptvscutvars[1] = 0.;
			xmax_eleptvscutvars[1] = 0.5;
		}else if(ih==11){
			//11: CosPAv0
			bins_eleptvscutvars[1] = 20;
			xmin_eleptvscutvars[1] = 0.95;
			xmax_eleptvscutvars[1] = 1.0;
		}else if(ih==12){
			//12:K0s masss
			bins_eleptvscutvars[1] = 50;
			xmin_eleptvscutvars[1] = 0.497-0.03;
			xmax_eleptvscutvars[1] = 0.497+0.03;
		}else if(ih==13 || ih==14){
			//13: nSigmaTPC(pr), nSigma(pi)
			bins_eleptvscutvars[1] = 20;
			xmin_eleptvscutvars[1] = -5;
			xmax_eleptvscutvars[1] = 5;
		}else if(ih==15){
			//15: eta
			bins_eleptvscutvars[1] = 30;
			xmin_eleptvscutvars[1] = -1.5;
			xmax_eleptvscutvars[1] = 1.5;
		}else if(ih==16){
			//16: Opening angle
			bins_eleptvscutvars[1] = 20;
			xmin_eleptvscutvars[1] = 0.;
			xmax_eleptvscutvars[1] = 3.141592/2;
		}

		fHistoElePtvsCutVarsRS[ih] = new THnSparseF(Form("fHistoElePtvsCutVarsRS[%d]",ih),"",3,bins_eleptvscutvars,xmin_eleptvscutvars,xmax_eleptvscutvars);
		fOutputAll->Add(fHistoElePtvsCutVarsRS[ih]);
		fHistoElePtvsCutVarsWS[ih] = new THnSparseF(Form("fHistoElePtvsCutVarsWS[%d]",ih),"",3,bins_eleptvscutvars,xmin_eleptvscutvars,xmax_eleptvscutvars);
		fOutputAll->Add(fHistoElePtvsCutVarsWS[ih]);
		fHistoElePtvsCutVarsMCS[ih] = new THnSparseF(Form("fHistoElePtvsCutVarsMCS[%d]",ih),"",3,bins_eleptvscutvars,xmin_eleptvscutvars,xmax_eleptvscutvars);
		fOutputAll->Add(fHistoElePtvsCutVarsMCS[ih]);
	}

  return;
}

//________________________________________________________________________
AliAODRecoCascadeHF* AliAnalysisTaskSELc2eleLambdafromAODtracks::MakeCascadeHF(AliAODv0 *v0, AliAODTrack *part, AliAODEvent * aod, AliAODVertex *secVert, Bool_t mixing) 
{
  //
  // Create AliAODRecoCascadeHF object from the argument
  //

  if(!v0) return 0x0;
  if(!part) return 0x0;
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
  AliESDtrack *esdtrack = new AliESDtrack((AliVTrack*)part);

  AliNeutralTrackParam *trackV0=NULL;
  const AliVTrack *trackVV0 = dynamic_cast<const AliVTrack*>(v0);
  if(trackVV0)  trackV0 = new AliNeutralTrackParam(trackVV0);

  Double_t xdummy, ydummy;
  Double_t dca = esdtrack->GetDCA(trackV0,fBzkG,xdummy,ydummy);


  //------------------------------------------------
  // Propagate all tracks to the secondary vertex and calculate momentum there
  //------------------------------------------------
	
  Double_t d0z0bach[2],covd0z0bach[3];
  if(sqrt(pow(secVert->GetX(),2)+pow(secVert->GetY(),2))<1.){
    part->PropagateToDCA(secVert,fBzkG,kVeryBig,d0z0bach,covd0z0bach);
    trackV0->PropagateToDCA(secVert,fBzkG,kVeryBig);
  }else{
    part->PropagateToDCA(primVertexAOD,fBzkG,kVeryBig,d0z0bach,covd0z0bach);
    trackV0->PropagateToDCA(primVertexAOD,fBzkG,kVeryBig);
  }
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

	if(!mixing){
		//If I add daughters for mixing event, I had some error. 
		theCascade->GetSecondaryVtx()->AddDaughter(part);
		theCascade->GetSecondaryVtx()->AddDaughter(v0);
	}

  if(unsetvtx) delete primVertexAOD; primVertexAOD=NULL;
  if(esdtrack) delete esdtrack;
  if(trackV0) delete trackV0;

  return theCascade;
}

//________________________________________________________________________
AliAODVertex* AliAnalysisTaskSELc2eleLambdafromAODtracks::CallPrimaryVertex(AliAODv0 *v0, AliAODTrack *trk, AliAODEvent* aod)
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
AliAODVertex* AliAnalysisTaskSELc2eleLambdafromAODtracks::PrimaryVertex(const TObjArray *trkArray,
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
AliAODVertex* AliAnalysisTaskSELc2eleLambdafromAODtracks::ReconstructSecondaryVertex(AliAODv0 *v0, AliAODTrack *part, AliAODEvent * aod) 
{
  //
  // Reconstruct secondary vertex from trkArray (Copied from AliAnalysisVertexingHF)
	// Currently only returns Primary vertex (can we reconstruct secondary vertex from e - v0??)
  //
	
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

  AliESDVertex * vertexESD = new AliESDVertex(*fV1);

  Double_t pos[3],cov[6],chi2perNDF;
  vertexESD->GetXYZ(pos); // position
  vertexESD->GetCovMatrix(cov); //covariance matrix
  chi2perNDF = vertexESD->GetChi2toNDF();
  delete vertexESD; vertexESD=NULL;
  
  AliAODVertex *secVert = new AliAODVertex(pos,cov,chi2perNDF);

  return secVert;
}
//________________________________________________________________________
Int_t AliAnalysisTaskSELc2eleLambdafromAODtracks::MatchToMC(AliAODRecoCascadeHF *elobj, TClonesArray *mcArray, Int_t *pdgarray_ele, Int_t *pdgarray_v0, Int_t *labelarray_ele, Int_t *labelarray_v0,  Int_t &ngen_ele, Int_t &ngen_v0) 
{
  //
  // Match to MC
  //
	for(Int_t i=0;i<100;i++){
		pdgarray_ele[i] = -9999;
		labelarray_ele[i] = -9999;
		pdgarray_v0[i] = -9999;
		labelarray_v0[i] = -9999;
	}
	ngen_ele = 0;
	ngen_v0 = 0;

  AliVTrack *trk = dynamic_cast<AliVTrack*>(elobj->GetBachelor());
  if(!trk) return -1;
  Int_t labEle = trk->GetLabel();
	if(labEle<0) return -1;
	AliAODMCParticle *mcetrk = (AliAODMCParticle*)mcArray->At(labEle);
	if(!mcetrk) return -1;
	labelarray_ele[0] = labEle;
	pdgarray_ele[0] = mcetrk->GetPdgCode();
	ngen_ele ++;

  AliAODMCParticle *mcprimele=0;
  mcprimele = mcetrk;
  while(mcprimele->GetMother()>=0) {
    Int_t labprim_ele=mcprimele->GetMother();
    AliAODMCParticle *tmcprimele = (AliAODMCParticle*)mcArray->At(labprim_ele);
    if(!tmcprimele) {
			break;
    }

    mcprimele = tmcprimele;
		pdgarray_ele[ngen_ele] = mcprimele->GetPdgCode();
		labelarray_ele[ngen_ele] = labprim_ele;
		ngen_ele ++;
		if(ngen_ele==100) break;
  }

  AliAODv0 *theV0 = dynamic_cast<AliAODv0*>(elobj->Getv0());
	if(!theV0) return -1;
	Int_t pdgdgv0[2]={2212,211};
  Int_t labV0 = theV0->MatchToMC(3122,mcArray,2,pdgdgv0); // the V0
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
	for(Int_t iemc=0;iemc<ngen_ele;iemc++){
		for(Int_t ivmc=0;ivmc<ngen_v0;ivmc++){
			if(labelarray_ele[iemc]==labelarray_v0[ivmc]){
				same_flag = kTRUE;
				matchedlabel = labelarray_ele[iemc];
				break;
			}
		}
		if(same_flag) break;
	}

	return matchedlabel;

}
//________________________________________________________________________
void AliAnalysisTaskSELc2eleLambdafromAODtracks::SelectTrack( const AliVEvent *event, Int_t trkEntries, Int_t &nSeleTrks,Bool_t *seleFlags, TClonesArray *mcArray)
{
  //
  // Select good tracks using fAnalCuts (AliRDHFCuts object) and return the array of their ids
  //
  
  if(trkEntries==0) return;
  
  nSeleTrks=0;
  for(Int_t i=0; i<trkEntries; i++) {
    seleFlags[i] = kFALSE;
    
    AliVTrack *track;
    track = (AliVTrack*)event->GetTrack(i);
    
    if(track->GetID()<0) continue;
    Double_t covtest[21];
    if(!track->GetCovarianceXYZPxPyPz(covtest)) continue;
    
    AliAODTrack *aodt = (AliAODTrack*)track;
		Double_t nsigma_tpcele = -9999;
		Double_t nsigma_tofele = -9999;
		if(fAnalCuts->GetIsUsePID()){
			nsigma_tpcele = fAnalCuts->GetPidHF()->GetPidResponse()->NumberOfSigmasTPC(aodt,AliPID::kElectron);
			nsigma_tofele = fAnalCuts->GetPidHF()->GetPidResponse()->NumberOfSigmasTOF(aodt,AliPID::kElectron);
		}

    if(!fAnalCuts) continue;
    if(fAnalCuts->SingleTrkCutsNoPID(aodt,fVtx1)){
			fHistoElectronTPCPID->Fill(aodt->Pt(),nsigma_tpcele);
			fHistoElectronTOFPID->Fill(aodt->Pt(),nsigma_tofele);
			if(fabs(nsigma_tofele)<3.){
				fHistoElectronTPCPIDSelTOF->Fill(aodt->Pt(),nsigma_tpcele);
				Double_t eleeta = aodt->Eta();
				if(fabs(eleeta)<0.6)
					fHistoElectronTPCPIDSelTOFSmallEta->Fill(aodt->Pt(),nsigma_tpcele);
				if(fabs(eleeta)>0.6 && fabs(eleeta)<0.8)
					fHistoElectronTPCPIDSelTOFLargeEta->Fill(aodt->Pt(),nsigma_tpcele);
				if(eleeta>-0.8 && eleeta<-0.6){
					fHistoElectronTPCPIDSelTOFEtaDep[0]->Fill(aodt->Pt(),nsigma_tpcele);
				}else if(eleeta>-0.6&&eleeta<-0.4){
					fHistoElectronTPCPIDSelTOFEtaDep[1]->Fill(aodt->Pt(),nsigma_tpcele);
				}else if(eleeta>-0.4&&eleeta<-0.2){
					fHistoElectronTPCPIDSelTOFEtaDep[2]->Fill(aodt->Pt(),nsigma_tpcele);
				}else if(eleeta>-0.2&&eleeta<0.0){
					fHistoElectronTPCPIDSelTOFEtaDep[3]->Fill(aodt->Pt(),nsigma_tpcele);
				}else if(eleeta>0.0&&eleeta<0.2){
					fHistoElectronTPCPIDSelTOFEtaDep[4]->Fill(aodt->Pt(),nsigma_tpcele);
				}else if(eleeta>0.2&&eleeta<0.4){
					fHistoElectronTPCPIDSelTOFEtaDep[5]->Fill(aodt->Pt(),nsigma_tpcele);
				}else if(eleeta>0.4&&eleeta<0.6){
					fHistoElectronTPCPIDSelTOFEtaDep[6]->Fill(aodt->Pt(),nsigma_tpcele);
				}else if(eleeta>0.6&&eleeta<0.8){
					fHistoElectronTPCPIDSelTOFEtaDep[7]->Fill(aodt->Pt(),nsigma_tpcele);
				}
			}
		}
    if(fAnalCuts->SingleTrkCuts(aodt,fVtx1)){
      seleFlags[i]=kTRUE;
      nSeleTrks++;
			fHistoElectronTPCSelPID->Fill(aodt->Pt(),nsigma_tpcele);
			fHistoElectronTOFSelPID->Fill(aodt->Pt(),nsigma_tofele);
			FillElectronROOTObjects(aodt,mcArray);
			if(fDoEventMixing){
				fElectronTracks->AddLast(new AliAODTrack(*aodt));
			}
    }

  } // end loop on tracks
}
//________________________________________________________________________
void AliAnalysisTaskSELc2eleLambdafromAODtracks::SelectV0( const AliVEvent *event,Int_t nV0s,Int_t &nSeleV0, Bool_t *seleV0Flags, TClonesArray *mcArray)
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
//_________________________________________________________________
Int_t AliAnalysisTaskSELc2eleLambdafromAODtracks::GetPoolIndex(Double_t zvert, Double_t mult){
	//
  // check in which of the pools the current event falls
	//

  Int_t theBinZ=TMath::BinarySearch(fNzVtxBins,fZvtxBins,zvert);
  if(theBinZ<0 || theBinZ>=fNzVtxBins) return -1;
  Int_t theBinM=TMath::BinarySearch(fNCentBins,fCentBins,mult);
  if(theBinM<0 || theBinM>=fNCentBins) return -1;
  return fNCentBins*theBinZ+theBinM;
}
//_________________________________________________________________
void AliAnalysisTaskSELc2eleLambdafromAODtracks::ResetPool(Int_t poolIndex){
	//
  // delete the contets of the pool
	//
  if(poolIndex<0 || poolIndex>=fNOfPools) return;
  delete fEventBuffer[poolIndex];
  fEventBuffer[poolIndex]=new TTree(Form("EventBuffer_%d",poolIndex), "Temporary buffer for event mixing");

	fEventBuffer[poolIndex]->Branch("zVertex", &fVtxZ);
	fEventBuffer[poolIndex]->Branch("centrality", &fCentrality);
	fEventBuffer[poolIndex]->Branch("eventInfo", "TObjString",&fEventInfo);
	fEventBuffer[poolIndex]->Branch("earray", "TObjArray", &fElectronTracks);

  return;
}
//_________________________________________________________________
void AliAnalysisTaskSELc2eleLambdafromAODtracks::DoEventMixingWithPools(Int_t poolIndex,AliAODEvent *aodEvent, Bool_t *seleFlags)
{
	//
  // perform mixed event analysis
	//

  if(poolIndex<0 || poolIndex>fNzVtxBins*fNCentBins) return;
	if(fEventBuffer[poolIndex]->GetEntries()<fNumberOfEventsForMixing) return;

  Int_t nEvents=fEventBuffer[poolIndex]->GetEntries();

  TObjArray* earray=0x0;
  Float_t zVertex,cent;
  TObjString* eventInfo=0x0;
  fEventBuffer[poolIndex]->SetBranchAddress("earray", &earray);
  fEventBuffer[poolIndex]->SetBranchAddress("eventInfo",&eventInfo);
  fEventBuffer[poolIndex]->SetBranchAddress("zVertex", &zVertex);
  fEventBuffer[poolIndex]->SetBranchAddress("centrality", &cent);
  for (Int_t i=0; i<aodEvent->GetNumberOfV0s(); i++)
  {
    if(!seleFlags[i]) continue;
    AliAODv0* v0 = aodEvent->GetV0(i);
    if(!v0)continue;

		for(Int_t iEv=0; iEv<fNumberOfEventsForMixing; iEv++){
			fEventBuffer[poolIndex]->GetEvent(iEv + nEvents - fNumberOfEventsForMixing);
			TObjArray* earray1=(TObjArray*)earray->Clone();
			//Float_t zVertex1=zVertex;
			//Float_t mult1=cent;
			Int_t nElectrons=earray1->GetEntries();
			//Int_t evId1,esdId1,ne1;
			//sscanf((eventInfo->String()).Data(),"Ev%d_esd%d_K%d",&evId1,&esdId1,&ne1);
//			if(ne1!=nElectrons){ 
//				printf("AliAnalysisTaskSELc2eleLambdafromAODtracks::DoMixingWithPools ERROR: read event does not match to the stored one\n");
//				delete earray1;
//				continue;
//			}
      for(Int_t iTr1=0; iTr1<nElectrons; iTr1++){
				AliAODTrack* trk1=(AliAODTrack*)earray1->At(iTr1);
				if(!trk1) continue;

        //if(!fAnalCuts->SelectWithRoughCuts(v0,trk1)) continue;
				AliAODVertex *secVert = ReconstructSecondaryVertex(v0,trk1,aodEvent);//Fake, prim vertex is just used as secondary vertex. place holder for future
				if(!secVert) continue;

				AliAODRecoCascadeHF *elobj = MakeCascadeHF(v0,trk1,aodEvent,secVert,true);
				if(!elobj) {
						continue;
				}

				TClonesArray *fake = 0;
				FillROOTObjects(elobj,v0,trk1,fake,true);

				elobj->GetSecondaryVtx()->RemoveDaughters();
				elobj->UnsetOwnPrimaryVtx();
				delete elobj;elobj=NULL;
				delete secVert;
			}//track loop

			delete earray1;
		}//event loop
		
	}//v0 loop
}
//_________________________________________________________________
void AliAnalysisTaskSELc2eleLambdafromAODtracks::MakeMCAnalysis(TClonesArray *mcArray)
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
			Bool_t e_flag = kFALSE;
			Bool_t lam_flag = kFALSE;
			AliAODMCParticle *mcepart = 0;
			AliAODMCParticle *mcv0part = 0;
			for(Int_t idau=mcpart->GetFirstDaughter();idau<mcpart->GetLastDaughter()+1;idau++)
			{
				if(idau<0) break;
				AliAODMCParticle *mcdau = (AliAODMCParticle*) mcArray->At(idau);
				if(!mcdau) continue;
				if(TMath::Abs(mcdau->GetPdgCode())==11){
					e_flag = kTRUE;
					mcepart = mcdau;
				}
				if(TMath::Abs(mcdau->GetPdgCode())==3122){
					lam_flag = kTRUE;
					mcv0part = mcdau;
				}
			}

			Int_t decaytype = -9999;
			if(e_flag && lam_flag) decaytype = 0;

			FillMCROOTObjects(mcpart,mcepart,mcv0part,decaytype);
		}
		if(TMath::Abs(mcpart->GetPdgCode())==4132){
			//cout<<"Lambdac"<<endl;
			Bool_t e_flag = kFALSE;
			Bool_t xi_flag = kFALSE;
			Bool_t lam_flag = kFALSE;
			AliAODMCParticle *mcepart = 0;
			AliAODMCParticle *mccascpart = 0;
			AliAODMCParticle *mcv0part = 0;
			for(Int_t idau=mcpart->GetFirstDaughter();idau<mcpart->GetLastDaughter()+1;idau++)
			{
				if(idau<0) break;
				AliAODMCParticle *mcdau = (AliAODMCParticle*) mcArray->At(idau);
				if(!mcdau) continue;
				if(TMath::Abs(mcdau->GetPdgCode())==11){
					e_flag = kTRUE;
					mcepart = mcdau;
				}
				if(TMath::Abs(mcdau->GetPdgCode())==3312){
					xi_flag = kTRUE;
					mccascpart = mcdau;
					for(Int_t idauxi=mccascpart->GetFirstDaughter();idauxi<mccascpart->GetLastDaughter()+1;idauxi++)
					{
						if(idauxi<0) break;
						AliAODMCParticle *mcdauxi = (AliAODMCParticle*) mcArray->At(idauxi);
						if(!mcdauxi) continue;
						if(TMath::Abs(mcdauxi->GetPdgCode())==3122){
							lam_flag = kTRUE;
							mcv0part = mcdauxi;
						}
					}
				}
			}
			Int_t decaytype = -9999;
			if(e_flag && xi_flag && lam_flag) decaytype = 1;

			FillMCROOTObjects(mcpart,mcepart,mcv0part,decaytype);
		}
		if(TMath::Abs(mcpart->GetPdgCode())==4232){
			//cout<<"Lambdac"<<endl;
			Bool_t e_flag = kFALSE;
			Bool_t xi_flag = kFALSE;
			Bool_t lam_flag = kFALSE;
			AliAODMCParticle *mcepart = 0;
			AliAODMCParticle *mccascpart = 0;
			AliAODMCParticle *mcv0part = 0;
			for(Int_t idau=mcpart->GetFirstDaughter();idau<mcpart->GetLastDaughter()+1;idau++)
			{
				if(idau<0) break;
				AliAODMCParticle *mcdau = (AliAODMCParticle*) mcArray->At(idau);
				if(!mcdau) continue;
				if(TMath::Abs(mcdau->GetPdgCode())==11){
					e_flag = kTRUE;
					mcepart = mcdau;
				}
				if(TMath::Abs(mcdau->GetPdgCode())==3322){
					xi_flag = kTRUE;
					mccascpart = mcdau;
					for(Int_t idauxi=mccascpart->GetFirstDaughter();idauxi<mccascpart->GetLastDaughter()+1;idauxi++)
					{
						if(idauxi<0) break;
						AliAODMCParticle *mcdauxi = (AliAODMCParticle*) mcArray->At(idauxi);
						if(!mcdauxi) continue;
						if(TMath::Abs(mcdauxi->GetPdgCode())==3122){
							lam_flag = kTRUE;
							mcv0part = mcdauxi;
						}
					}
				}
			}
			Int_t decaytype = -9999;
			if(e_flag && xi_flag && lam_flag) decaytype = 2;

			FillMCROOTObjects(mcpart,mcepart,mcv0part,decaytype);
		}

		if(TMath::Abs(mcpart->GetPdgCode())==11 && mcpart->GetStatus()==1){
			AliESDtrackCuts *esdcuts = fAnalCuts->GetTrackCuts();
			Float_t etamin, etamax;
			esdcuts->GetEtaRange(etamin,etamax);
			if(fabs(mcpart->Eta())<etamax)
				fHistoBachPtMCGen->Fill(mcpart->Pt());
		}
		if(TMath::Abs(mcpart->GetPdgCode())==3122){
			Double_t etamin, etamax, rapmin, rapmax;
			fAnalCuts->GetProdV0EtaRange(etamin,etamax);
			fAnalCuts->GetProdV0RapRange(rapmin,rapmax);

			if((fabs(mcpart->Y())<rapmax) && (fabs(mcpart->Eta())<etamax))
				fHistoLambdaMassvsPtMCGen->Fill(1.115683, mcpart->Pt());
		}
	}
	return;
}
