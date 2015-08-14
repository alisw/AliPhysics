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
//               Xic->eXi analysis code
//
//  Input: AOD
//  Output: TTree or THnSparse (mass vs pT vs Centrality)
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
#include "AliAnalysisTaskSEXic2eleXifromAODtracks.h"
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

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskSEXic2eleXifromAODtracks);
/// \endcond

//__________________________________________________________________________
AliAnalysisTaskSEXic2eleXifromAODtracks::AliAnalysisTaskSEXic2eleXifromAODtracks() : 
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
  fCascVariablesTree(0),
  fMCVariablesTree(0),
  fReconstructPrimVert(kFALSE),
  fIsMB(kFALSE),
  fIsSemi(kFALSE),
  fIsCent(kFALSE),
  fIsINT7(kFALSE),
  fIsEMC7(kFALSE),
  fCandidateVariables(),
  fCandidateEleVariables(),
  fCandidateCascVariables(),
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
  fHistoEleXiMass(0),
  fHistoEleXiMassRS(0),
  fHistoEleXiMassWS(0),
  fHistoEleXiMassRSMix(0),
  fHistoEleXiMassWSMix(0),
  fHistoEleXiMassvsElePtRS(0),
  fHistoEleXiMassvsElePtWS(0),
  fHistoEleXiMassvsElePtRSMix(0),
  fHistoEleXiMassvsElePtWSMix(0),
  fHistoElePtRS(0),
  fHistoElePtWS(0),
  fHistoElePtRSMix(0),
  fHistoElePtWSMix(0),
  fHistoEleXiMassMCS(0),
  fHistoEleXiMassMCGen(0),
  fHistoEleXiMassvsElePtMCS(0),
  fHistoEleXiMassvsElePtMCGen(0),
  fHistoElePtMCS(0),
  fHistoElePtMCGen(0),
  fHistoElePtvsEtaRS(0),
  fHistoElePtvsEtaWS(0),
  fHistoElePtvsEtaRSMix(0),
  fHistoElePtvsEtaWSMix(0),
  fHistoElePtvsEtaMCS(0),
  fHistoElePtvsEtaMCGen(0),
  fHistoElePtvsXiPtRS(0),
  fHistoElePtvsXiPtWS(0),
  fHistoElePtvsXiPtRSMix(0),
  fHistoElePtvsXiPtWSMix(0),
  fHistoElePtvsXiPtMCS(0),
  fHistoElePtvsXiPtvsXicPtMCS(0),
  fHistoElePtvsXiPtMCGen(0),
  fHistoElePtvsXiPtvsXicPtMCGen(0),
  fHistoElePtvsXiPtMCXicGen(0),
  fHistoElePtvsd0RS(0),
  fHistoElePtvsd0WS(0),
  fHistoElePtvsd0RSMix(0),
  fHistoElePtvsd0WSMix(0),
  fHistoElePtvsd0MCS(0),
  fHistoElePtvsd0PromptMCS(0),
  fHistoElePtvsd0BFeeddownMCS(0),
  fHistoBachPt(0),
  fHistoBachPtMCS(0),
  fHistoBachPtMCGen(0),
  fHistod0Bach(0),
  fHistoXiMassvsPt(0),
  fHistoXiMassvsPtMCS(0),
  fHistoXiMassvsPtMCGen(0),
  fHistoOmegaMassvsPt(0),
  fHistoElectronTPCPID(0),
  fHistoElectronTOFPID(0),
  fHistoElectronTPCSelPID(0),
  fHistoElectronTOFSelPID(0),
  fHistoElectronTPCPIDSelTOF(0),
  fHistoElectronTPCPIDSelTOFSmallEta(0),
  fHistoElectronTPCPIDSelTOFLargeEta(0),
	fHistoElectronQovPtvsPhi(0),
	fHistoXiQovPtvsPhi(0),
	fCounter(0),
	fHistonEvtvsRunNumber(0),
	fHistonElevsRunNumber(0),
	fHistonXivsRunNumber(0),
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
	for(Int_t i=0;i<23;i++){
		fHistoElePtvsCutVarsRS[i] = 0;
		fHistoElePtvsCutVarsWS[i] = 0;
		fHistoElePtvsCutVarsMCS[i] = 0;
	}
	for(Int_t i=0;i<8;i++){
		fHistoElectronTPCPIDSelTOFEtaDep[i] = 0;
	}
}

//___________________________________________________________________________
AliAnalysisTaskSEXic2eleXifromAODtracks::AliAnalysisTaskSEXic2eleXifromAODtracks(const Char_t* name,
									     AliRDHFCutsXictoeleXifromAODtracks* analCuts, 
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
  fCascVariablesTree(0),
  fMCVariablesTree(0),
  fReconstructPrimVert(kFALSE),
  fIsMB(kFALSE),
  fIsSemi(kFALSE),
  fIsCent(kFALSE),
  fIsINT7(kFALSE),
  fIsEMC7(kFALSE),
  fCandidateVariables(),
  fCandidateEleVariables(),
  fCandidateCascVariables(),
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
  fHistoEleXiMass(0),
  fHistoEleXiMassRS(0),
  fHistoEleXiMassWS(0),
  fHistoEleXiMassRSMix(0),
  fHistoEleXiMassWSMix(0),
  fHistoEleXiMassvsElePtRS(0),
  fHistoEleXiMassvsElePtWS(0),
  fHistoEleXiMassvsElePtRSMix(0),
  fHistoEleXiMassvsElePtWSMix(0),
  fHistoElePtRS(0),
  fHistoElePtWS(0),
  fHistoElePtRSMix(0),
  fHistoElePtWSMix(0),
  fHistoEleXiMassMCS(0),
  fHistoEleXiMassMCGen(0),
  fHistoEleXiMassvsElePtMCS(0),
  fHistoEleXiMassvsElePtMCGen(0),
  fHistoElePtMCS(0),
  fHistoElePtMCGen(0),
  fHistoElePtvsEtaRS(0),
  fHistoElePtvsEtaWS(0),
  fHistoElePtvsEtaRSMix(0),
  fHistoElePtvsEtaWSMix(0),
  fHistoElePtvsEtaMCS(0),
  fHistoElePtvsEtaMCGen(0),
  fHistoElePtvsXiPtRS(0),
  fHistoElePtvsXiPtWS(0),
  fHistoElePtvsXiPtRSMix(0),
  fHistoElePtvsXiPtWSMix(0),
  fHistoElePtvsXiPtMCS(0),
  fHistoElePtvsXiPtvsXicPtMCS(0),
  fHistoElePtvsXiPtMCGen(0),
  fHistoElePtvsXiPtvsXicPtMCGen(0),
  fHistoElePtvsXiPtMCXicGen(0),
  fHistoElePtvsd0RS(0),
  fHistoElePtvsd0WS(0),
  fHistoElePtvsd0RSMix(0),
  fHistoElePtvsd0WSMix(0),
  fHistoElePtvsd0MCS(0),
  fHistoElePtvsd0PromptMCS(0),
  fHistoElePtvsd0BFeeddownMCS(0),
  fHistoBachPt(0),
  fHistoBachPtMCS(0),
  fHistoBachPtMCGen(0),
  fHistod0Bach(0),
  fHistoXiMassvsPt(0),
  fHistoXiMassvsPtMCS(0),
  fHistoXiMassvsPtMCGen(0),
  fHistoOmegaMassvsPt(0),
  fHistoElectronTPCPID(0),
  fHistoElectronTOFPID(0),
  fHistoElectronTPCSelPID(0),
  fHistoElectronTOFSelPID(0),
  fHistoElectronTPCPIDSelTOF(0),
  fHistoElectronTPCPIDSelTOFSmallEta(0),
  fHistoElectronTPCPIDSelTOFLargeEta(0),
	fHistoElectronQovPtvsPhi(0),
	fHistoXiQovPtvsPhi(0),
	fCounter(0),
	fHistonEvtvsRunNumber(0),
	fHistonElevsRunNumber(0),
	fHistonXivsRunNumber(0),
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
  Info("AliAnalysisTaskSEXic2eleXifromAODtracks","Calling Constructor");

	for(Int_t i=0;i<23;i++){
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
AliAnalysisTaskSEXic2eleXifromAODtracks::~AliAnalysisTaskSEXic2eleXifromAODtracks() {
  //
  // destructor
  //
  Info("~AliAnalysisTaskSEXic2eleXifromAODtracks","Calling Destructor");

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
  if (fCascVariablesTree) {
    delete fCascVariablesTree;
    fCascVariablesTree = 0;
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
void AliAnalysisTaskSEXic2eleXifromAODtracks::Init() {
  //
  // Initialization
  //
  //

  fIsEventSelected=kFALSE;

  if (fDebug > 1) AliInfo("Init");

  fListCuts = new TList();
  fListCuts->SetOwner();
  fListCuts->SetName("ListCuts");
  fListCuts->Add(new AliRDHFCutsXictoeleXifromAODtracks(*fAnalCuts));
  PostData(2,fListCuts);

  return;
}

//_________________________________________________
void AliAnalysisTaskSEXic2eleXifromAODtracks::UserExec(Option_t *)
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
      AliError("AliAnalysisTaskSEXic2eleXifromAODtracks::UserExec: MC header branch not found!\n");
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
  Bool_t fIsPhysSelNotOK = fAnalCuts->IsEventRejectedDuePhysicsSelection();
  Bool_t fIsNoVertex = fAnalCuts->IsEventRejectedDueToNotRecoVertex();
  if(!fIsTriggerNotOK && !fIsPhysSelNotOK && !fIsNoVertex && fabs(fVtx1->GetZ())<fAnalCuts->GetMaxVtxZ()) fCEvents->Fill(3);
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
  PostData(6,fCascVariablesTree);
  PostData(7,fMCVariablesTree);
  PostData(8,fCounter);    

  fIsEventSelected=kFALSE;

  delete fV1;
  return;
}

//________________________________________ terminate ___________________________
void AliAnalysisTaskSEXic2eleXifromAODtracks::Terminate(Option_t*)
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
void AliAnalysisTaskSEXic2eleXifromAODtracks::UserCreateOutputObjects() 
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

  DefineCascTreeVariables();
  PostData(6,fCascVariablesTree);

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
void AliAnalysisTaskSEXic2eleXifromAODtracks::MakeAnalysis
(
 AliAODEvent *aodEvent, TClonesArray *mcArray
 )
{
  //
  // Main Analysis part
  //
  //------------------------------------------------
  // Select good track before hand to save time
  //------------------------------------------------
	if(fDoEventMixing && fElectronTracks) fElectronTracks->Delete();

  Int_t nCascs= aodEvent->GetNumberOfCascades();
  Int_t nTracks= aodEvent->GetNumberOfTracks();

  Bool_t  seleTrkFlags[nTracks];
  Int_t nSeleTrks=0;
  SelectTrack(aodEvent,nTracks,nSeleTrks,seleTrkFlags,mcArray);

  Bool_t  seleCascFlags[nCascs];
  Int_t     nSeleCasc=0;
  SelectCascade(aodEvent,nCascs,nSeleCasc,seleCascFlags,mcArray);

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
	fHistonXivsRunNumber->Fill(runnumber-runnumber_offset,nSeleCasc);

  //------------------------------------------------
  // Cascade loop 
  //------------------------------------------------
  for (Int_t icasc = 0; icasc<nCascs; icasc++) {
    if(!seleCascFlags[icasc]) continue;
    AliAODcascade *casc = aodEvent->GetCascade(icasc);
    if(!casc) continue;

    AliAODTrack *cptrack =  (AliAODTrack*)(casc->GetDaughter(0));
    AliAODTrack *cntrack =  (AliAODTrack*)(casc->GetDaughter(1));
    AliAODTrack *cbtrack =  (AliAODTrack*)(casc->GetDecayVertexXi()->GetDaughter(0));

    //------------------------------------------------
    // track loop 
    //------------------------------------------------
    for (Int_t itrk = 0; itrk<nTracks; itrk++) {
      if(!seleTrkFlags[itrk]) continue;
      AliAODTrack *trk = (AliAODTrack*)aodEvent->GetTrack(itrk);
      if(trk->GetID()<0) continue;

      Int_t cpid = cptrack->GetID();
      Int_t cnid = cntrack->GetID();
      Int_t cbid = cbtrack->GetID();
      Int_t lpid = trk->GetID();
      if((cpid==lpid)||(cnid==lpid)||(cbid==lpid)) continue;

      //if(!fAnalCuts->SelectWithRoughCuts(v0,trk)) continue;

      AliAODVertex *secVert = ReconstructSecondaryVertex(casc,trk,aodEvent);//Fake, prim vertex is just used as secondary vertex. place holder for future
      if(!secVert) continue;

      AliAODRecoCascadeHF *exobj = MakeCascadeHF(casc,trk,aodEvent,secVert,false);
      if(!exobj) {
	continue;
      }

      FillROOTObjects(exobj, casc,trk,mcArray,false);

      exobj->GetSecondaryVtx()->RemoveDaughters();
      exobj->UnsetOwnPrimaryVtx();
      delete exobj;exobj=NULL;
      delete secVert;
    }
  }

  if(fDoEventMixing){
		fEventInfo->SetString(Form("Ev%d_esd%d_E%d",AliAnalysisManager::GetAnalysisManager()->GetNcalls(),((AliAODHeader*)aodEvent->GetHeader())->GetEventNumberESDFile(),fElectronTracks->GetEntries()));
    Int_t ind=GetPoolIndex(fVtxZ,fCentrality);
    if(ind>=0 && ind<fNOfPools){
      if(fEventBuffer[ind]->GetEntries() >= fNumberOfEventsForMixing){
				DoEventMixingWithPools(ind,aodEvent,seleCascFlags);
				//ResetPool(ind);
      }
      fEventBuffer[ind]->Fill();
    }
  }
}


////-------------------------------------------------------------------------------
void AliAnalysisTaskSEXic2eleXifromAODtracks::DefineTreeVariables() 
{
  //
  // Define tree variables
  //

  const char* nameoutput = GetOutputSlot(4)->GetContainer()->GetName();
  fVariablesTree = new TTree(nameoutput,"Candidates variables tree");
  Int_t nVar = 86;
  fCandidateVariables = new Float_t [nVar];
  TString * fCandidateVariableNames = new TString[nVar];

  fCandidateVariableNames[ 0]="Centrality";
  fCandidateVariableNames[ 1]="InvMassEleXi";
  fCandidateVariableNames[ 2]="EleXiPt";
  fCandidateVariableNames[ 3]="EleXiPx";
  fCandidateVariableNames[ 4]="EleXiPy";
  fCandidateVariableNames[ 5]="EleXiPz";
  fCandidateVariableNames[ 6]="ElePx";
  fCandidateVariableNames[ 7]="ElePy";
  fCandidateVariableNames[ 8]="ElePz";
  fCandidateVariableNames[ 9]="XiPx";
  fCandidateVariableNames[10]="XiPy";
  fCandidateVariableNames[11]="XiPz";
  fCandidateVariableNames[12]="XiCharge";
  fCandidateVariableNames[13]="MassXi";
  fCandidateVariableNames[14]="MassLambda";
  fCandidateVariableNames[15]="Eled0";
  fCandidateVariableNames[16]="Xid0";
  fCandidateVariableNames[17]="nSigmaTPCele";
  fCandidateVariableNames[18]="nSigmaTOFele";
  fCandidateVariableNames[19]="nSigmaTPCpr_etrk";
  fCandidateVariableNames[20]="nSigmaTOFpr_etrk";
  fCandidateVariableNames[21]="nSigmaTPCka_etrk";
  fCandidateVariableNames[22]="nSigmaTOFka_etrk";
  fCandidateVariableNames[23]="nSigmaTPCv0pr";
  fCandidateVariableNames[24]="nSigmaTOFv0pr";
  fCandidateVariableNames[25]="nSigmaTPCv0pi";
  fCandidateVariableNames[26]="nSigmaTOFv0pi";
  fCandidateVariableNames[27]="nSigmaTPCbachpi";
  fCandidateVariableNames[28]="nSigmaTOFbachpi";
  fCandidateVariableNames[29]="EleCharge";
  fCandidateVariableNames[30]="Mixing";
  fCandidateVariableNames[31]="DcaXiDaughters";
  fCandidateVariableNames[32]="DcaV0Daughters";
  fCandidateVariableNames[33]="DecayLengthXi";
  fCandidateVariableNames[34]="CosPointingAngleXi";
  fCandidateVariableNames[35]="DcaV0toPrimVertex";
  fCandidateVariableNames[36]="DcaPostoPrimVertex";
  fCandidateVariableNames[37]="DcaNegtoPrimVertex";
  fCandidateVariableNames[38]="DcaBachtoPrimVertex";
  fCandidateVariableNames[39]="DecayLengthV0";
  fCandidateVariableNames[40]="CosPointingAngleV0";

  fCandidateVariableNames[41]="mcpdgxic";
  fCandidateVariableNames[42]="mclabxic";
  fCandidateVariableNames[43]="mcxicpx";
  fCandidateVariableNames[44]="mcxicpy";
  fCandidateVariableNames[45]="mcxicpz";
  fCandidateVariableNames[46]="mcelepx";
  fCandidateVariableNames[47]="mcelepy";
  fCandidateVariableNames[48]="mcelepz";
  fCandidateVariableNames[49]="mccascpx";
  fCandidateVariableNames[50]="mccascpy";
  fCandidateVariableNames[51]="mccascpz";

  fCandidateVariableNames[52]="mcpdgele";
  fCandidateVariableNames[53]="mcpdgcasc";
  fCandidateVariableNames[54]="mcpdgmomele";
  fCandidateVariableNames[55]="mcpdgmomcasc";
  fCandidateVariableNames[56]="mcpdggrmomele";
  fCandidateVariableNames[57]="mcpdggrmomcasc";
  fCandidateVariableNames[58]="mcngenele";
  fCandidateVariableNames[59]="mcngencasc";

  fCandidateVariableNames[60]="nSigmaTPCpi_etrk";
  fCandidateVariableNames[61]="nSigmaTOFpi_etrk";

  fCandidateVariableNames[62]="V0PosPx";
  fCandidateVariableNames[63]="V0PosPy";
  fCandidateVariableNames[64]="V0PosPz";
  fCandidateVariableNames[65]="V0NegPx";
  fCandidateVariableNames[66]="V0NegPy";
  fCandidateVariableNames[67]="V0NegPz";
  fCandidateVariableNames[68]="V0VertX";
  fCandidateVariableNames[69]="V0VertY";
  fCandidateVariableNames[70]="V0VertZ";
  fCandidateVariableNames[71]="BachPx";
  fCandidateVariableNames[72]="BachPy";
  fCandidateVariableNames[73]="BachPz";
  fCandidateVariableNames[74]="XiVertX";
  fCandidateVariableNames[75]="XiVertY";
  fCandidateVariableNames[76]="XiVertZ";
  fCandidateVariableNames[77]="PrimVertX";
  fCandidateVariableNames[78]="PrimVertY";
  fCandidateVariableNames[79]="PrimVertZ";

  fCandidateVariableNames[80]="MassOmega";

	fCandidateVariableNames[81]= "EleITSMatch";
	fCandidateVariableNames[82]= "BachITSMatch";
	fCandidateVariableNames[83]= "V0PosITSMatch";
	fCandidateVariableNames[84]= "V0NegITSMatch";

  fCandidateVariableNames[85]="EvNumber";

  for (Int_t ivar=0; ivar<nVar; ivar++) {
    fVariablesTree->Branch(fCandidateVariableNames[ivar].Data(),&fCandidateVariables[ivar],Form("%s/f",fCandidateVariableNames[ivar].Data()));
  }

  return;
}

////-------------------------------------------------------------------------------
void AliAnalysisTaskSEXic2eleXifromAODtracks::FillROOTObjects(AliAODRecoCascadeHF *exobj, AliAODcascade *casc, AliAODTrack *trk, TClonesArray *mcArray, Bool_t mixing_flag) 
{
  //
  // Fill histograms or tree depending on fWriteVariableTree 
  //
	if(!trk) return;
	if(!casc) return;

	for(Int_t i=0;i<86;i++){
		fCandidateVariables[i] = -9999.;
	}


  AliAODTrack *cptrack =  (AliAODTrack*)(casc->GetDaughter(0));
  AliAODTrack *cntrack =  (AliAODTrack*)(casc->GetDaughter(1));
  AliAODTrack *cbtrack =  (AliAODTrack*)(casc->GetDecayVertexXi()->GetDaughter(0));
	if(cptrack->Charge()<0 && cntrack->Charge()>0){
		cptrack =   (AliAODTrack*)(casc->GetDaughter(1));
		cntrack =   (AliAODTrack*)(casc->GetDaughter(0));
	}
//  Double_t d0z0[2],covd0z0[3];
//  cptrack->PropagateToDCA(fVtx1,fBzkG,kVeryBig,d0z0,covd0z0);//propagate to primary vertex for debugging
//  cntrack->PropagateToDCA(fVtx1,fBzkG,kVeryBig,d0z0,covd0z0);//propagate to primary vertex for debugging
//  cbtrack->PropagateToDCA(fVtx1,fBzkG,kVeryBig,d0z0,covd0z0);//propagate to primary vertex for debugging

  Double_t posVtx[3] = {0.,0.,0.};
  fVtx1->GetXYZ(posVtx);


  fCandidateVariables[ 0] = fCentrality;
	UInt_t pdgdg[2]={11,3312};
  fCandidateVariables[ 1] = exobj->InvMass(2,pdgdg);
  fCandidateVariables[ 2] = exobj->Pt();
  fCandidateVariables[ 3] = exobj->Px();
  fCandidateVariables[ 4] = exobj->Py();
  fCandidateVariables[ 5] = exobj->Pz();
  fCandidateVariables[ 6] = exobj->PxProng(0);
  fCandidateVariables[ 7] = exobj->PyProng(0);
  fCandidateVariables[ 8] = exobj->PzProng(0);
  fCandidateVariables[ 9] = exobj->PxProng(1);
  fCandidateVariables[10] = exobj->PyProng(1);
  fCandidateVariables[11] = exobj->PzProng(1);
  fCandidateVariables[12] = casc->ChargeXi();
  fCandidateVariables[13] = casc->MassXi();
	if(casc->ChargeXi()<0)
		fCandidateVariables[14] = casc->MassLambda();
	else
		fCandidateVariables[14] = casc->MassAntiLambda();
  fCandidateVariables[15] = exobj->Getd0Prong(0);
  fCandidateVariables[16] = exobj->Getd0Prong(1);

	Double_t nSigmaTPCele = -9999.;
	Double_t nSigmaTOFele = -9999.;
  if(fAnalCuts->GetIsUsePID()&&!mixing_flag)
  {
		nSigmaTPCele = fAnalCuts->GetPidHF()->GetPidResponse()->NumberOfSigmasTPC(trk,AliPID::kElectron);
		nSigmaTOFele = fAnalCuts->GetPidHF()->GetPidResponse()->NumberOfSigmasTOF(trk,AliPID::kElectron);
    fCandidateVariables[17] = nSigmaTPCele;
    fCandidateVariables[18] = nSigmaTOFele;

		Double_t nSigmaTPCpr_etrk = fAnalCuts->GetPidHF()->GetPidResponse()->NumberOfSigmasTPC(trk,AliPID::kProton);
		Double_t nSigmaTOFpr_etrk = fAnalCuts->GetPidHF()->GetPidResponse()->NumberOfSigmasTOF(trk,AliPID::kProton);
		Double_t nSigmaTPCka_etrk = fAnalCuts->GetPidHF()->GetPidResponse()->NumberOfSigmasTPC(trk,AliPID::kKaon);
		Double_t nSigmaTOFka_etrk = fAnalCuts->GetPidHF()->GetPidResponse()->NumberOfSigmasTOF(trk,AliPID::kKaon);
		Double_t nSigmaTPCpi_etrk = fAnalCuts->GetPidHF()->GetPidResponse()->NumberOfSigmasTPC(trk,AliPID::kPion);
		Double_t nSigmaTOFpi_etrk = fAnalCuts->GetPidHF()->GetPidResponse()->NumberOfSigmasTOF(trk,AliPID::kPion);
    fCandidateVariables[19] = nSigmaTPCpr_etrk;
    fCandidateVariables[20] = nSigmaTOFpr_etrk;
    fCandidateVariables[21] = nSigmaTPCka_etrk;
    fCandidateVariables[22] = nSigmaTOFka_etrk;
    fCandidateVariables[60] = nSigmaTPCpi_etrk;
    fCandidateVariables[61] = nSigmaTOFpi_etrk;
  }

	Double_t nSigmaTPCv0pr=-9999.;
	Double_t nSigmaTOFv0pr=-9999.;
	Double_t nSigmaTPCv0pi=-9999.;
	Double_t nSigmaTOFv0pi=-9999.;
	Double_t nSigmaTPCbachpi=-9999.;
	Double_t nSigmaTOFbachpi=-9999.;
	if(fAnalCuts->GetUseCascadePID()&&!mixing_flag)
	{
		if(casc->ChargeXi()>0){
			nSigmaTPCv0pr = fAnalCuts->GetPidCascPr()->GetPidResponse()->NumberOfSigmasTPC(cntrack,AliPID::kProton);
			nSigmaTOFv0pr = fAnalCuts->GetPidCascPr()->GetPidResponse()->NumberOfSigmasTOF(cntrack,AliPID::kProton);
			nSigmaTPCv0pi = fAnalCuts->GetPidCascPi()->GetPidResponse()->NumberOfSigmasTPC(cptrack,AliPID::kPion);
			nSigmaTOFv0pi = fAnalCuts->GetPidCascPi()->GetPidResponse()->NumberOfSigmasTOF(cptrack,AliPID::kPion);
			nSigmaTPCbachpi = fAnalCuts->GetPidCascPi()->GetPidResponse()->NumberOfSigmasTPC(cbtrack,AliPID::kPion);
			nSigmaTOFbachpi = fAnalCuts->GetPidCascPi()->GetPidResponse()->NumberOfSigmasTOF(cbtrack,AliPID::kPion);
		}else{
			nSigmaTPCv0pr = fAnalCuts->GetPidCascPr()->GetPidResponse()->NumberOfSigmasTPC(cptrack,AliPID::kProton);
			nSigmaTOFv0pr = fAnalCuts->GetPidCascPr()->GetPidResponse()->NumberOfSigmasTOF(cptrack,AliPID::kProton);
			nSigmaTPCv0pi = fAnalCuts->GetPidCascPi()->GetPidResponse()->NumberOfSigmasTPC(cntrack,AliPID::kPion);
			nSigmaTOFv0pi = fAnalCuts->GetPidCascPi()->GetPidResponse()->NumberOfSigmasTOF(cntrack,AliPID::kPion);
			nSigmaTPCbachpi = fAnalCuts->GetPidCascPi()->GetPidResponse()->NumberOfSigmasTPC(cbtrack,AliPID::kPion);
			nSigmaTOFbachpi = fAnalCuts->GetPidCascPi()->GetPidResponse()->NumberOfSigmasTOF(cbtrack,AliPID::kPion);
		}
      fCandidateVariables[23] = nSigmaTPCv0pr;
      fCandidateVariables[24] = nSigmaTOFv0pr;
      fCandidateVariables[25] = nSigmaTPCv0pi;
      fCandidateVariables[26] = nSigmaTOFv0pi;
      fCandidateVariables[27] = nSigmaTPCbachpi;
      fCandidateVariables[28] = nSigmaTOFbachpi;
  }
  fCandidateVariables[29] = trk->Charge();
  fCandidateVariables[30] = (Float_t) mixing_flag;
  fCandidateVariables[31] = casc->DcaXiDaughters();
  fCandidateVariables[32] = casc->DcaV0Daughters();
  fCandidateVariables[33] = casc->DecayLengthXi(posVtx[0],posVtx[1],posVtx[2]);
  fCandidateVariables[34] = casc->CosPointingAngleXi(posVtx[0],posVtx[1],posVtx[2]);
  fCandidateVariables[35] = casc->DcaV0ToPrimVertex();
  fCandidateVariables[36] = casc->DcaPosToPrimVertex();
  fCandidateVariables[37] = casc->DcaNegToPrimVertex();
  fCandidateVariables[38] = casc->DcaBachToPrimVertex();
  fCandidateVariables[39] = casc->DecayLengthV0();
  fCandidateVariables[40] = casc->CosPointingAngle(casc->GetDecayVertexXi());

  AliAODMCParticle *mcxic = 0;
  AliAODMCParticle *mcele = 0;
  AliAODMCParticle *mccasc = 0;
  Int_t mclabxic = 0;
	Int_t mcpdgele_array[100];
	Int_t mcpdgcasc_array[100];
	Int_t mclabelele_array[100];
	Int_t mclabelcasc_array[100];
	Int_t mcngen_ele = -9999;
	Int_t mcngen_casc = -9999;

	if(fUseMCInfo && mcArray){

		mclabxic =  MatchToMC(exobj,mcArray,mcpdgele_array, mcpdgcasc_array,mclabelele_array,mclabelcasc_array,mcngen_ele,mcngen_casc);

    if(mclabxic>-1){
      mcxic = (AliAODMCParticle*) mcArray->At(mclabxic);
			if(mclabelele_array[0]>=0)
				mcele = (AliAODMCParticle*) mcArray->At(mclabelele_array[0]);
			if(mclabelcasc_array[0]>=0)
				mccasc = (AliAODMCParticle*) mcArray->At(mclabelcasc_array[0]);
			if(mcxic){
				fCandidateVariables[41] = mcxic->GetPdgCode();
				fCandidateVariables[42] = mcxic->Label();
				fCandidateVariables[43] = mcxic->Px();
				fCandidateVariables[44] = mcxic->Py();
				fCandidateVariables[45] = mcxic->Pz();
			}
			if(mcele){
				fCandidateVariables[46] = mcele->Px();
				fCandidateVariables[47] = mcele->Py();
				fCandidateVariables[48] = mcele->Pz();
			}
			if(mccasc){
				fCandidateVariables[49] = mccasc->Px();
				fCandidateVariables[50] = mccasc->Py();
				fCandidateVariables[51] = mccasc->Pz();
			}
		}
		fCandidateVariables[52] = mcpdgele_array[0];
		fCandidateVariables[53] = mcpdgcasc_array[0];
		fCandidateVariables[54] = mcpdgele_array[1];
		fCandidateVariables[55] = mcpdgcasc_array[1];
		fCandidateVariables[56] = mcpdgele_array[2];
		fCandidateVariables[57] = mcpdgcasc_array[2];
		fCandidateVariables[58] = mcngen_ele;
		fCandidateVariables[59] = mcngen_casc;
	}
	fCandidateVariables[62] = casc->MomPosX();
	fCandidateVariables[63] = casc->MomPosY();
	fCandidateVariables[64] = casc->MomPosZ();
	fCandidateVariables[65] = casc->MomNegX();
	fCandidateVariables[66] = casc->MomNegY();
	fCandidateVariables[67] = casc->MomNegZ();
	fCandidateVariables[68] = casc->DecayVertexV0X();
	fCandidateVariables[69] = casc->DecayVertexV0Y();
	fCandidateVariables[70] = casc->DecayVertexV0Z();
	fCandidateVariables[71] = casc->MomBachX();
	fCandidateVariables[72] = casc->MomBachY();
	fCandidateVariables[73] = casc->MomBachZ();
	fCandidateVariables[74] = casc->DecayVertexXiX();
	fCandidateVariables[75] = casc->DecayVertexXiY();
	fCandidateVariables[76] = casc->DecayVertexXiZ();
	fCandidateVariables[77] = fVtx1->GetX();
	fCandidateVariables[78] = fVtx1->GetY();
	fCandidateVariables[79] = fVtx1->GetZ();

	fCandidateVariables[80] = casc->MassOmega();

	if(trk) fCandidateVariables[81] = trk->GetITSClusterMap();
	if(cbtrack) fCandidateVariables[82] = cbtrack->GetITSClusterMap();
	if(cptrack) fCandidateVariables[83] = cptrack->GetITSClusterMap();
	if(cntrack) fCandidateVariables[84] = cntrack->GetITSClusterMap();

  fCandidateVariables[85] = fEvNumberCounter;

  if(fWriteVariableTree)
    fVariablesTree->Fill();

	if(fAnalCuts->IsSelected(exobj,AliRDHFCuts::kCandidate))
	{
		Double_t cont[3];
		cont[0] = exobj->InvMass(2,pdgdg);
		cont[1] = exobj->Pt();
		cont[2] = fCentrality;
		fHistoEleXiMass->Fill(cont);

		Double_t cont2[3];
		cont2[0] = exobj->InvMass(2,pdgdg);
		cont2[1] = trk->Pt();
		cont2[2] = fCentrality;

		Double_t cont_eleptvseta[3];
		cont_eleptvseta[0] = trk->Pt();
		cont_eleptvseta[1] = trk->Eta();
		cont_eleptvseta[2] = fCentrality;

		Double_t cont_eleptvsxipt[3];
		cont_eleptvsxipt[0] = trk->Pt();
		cont_eleptvsxipt[1] = sqrt(casc->MomXiX()*casc->MomXiX()+casc->MomXiY()*casc->MomXiY());
		cont_eleptvsxipt[2] = fCentrality;

		Double_t cont_eleptvsd0[3];
		cont_eleptvsd0[0] = trk->Pt();
		cont_eleptvsd0[1] = exobj->Getd0Prong(0)*trk->Charge();
		cont_eleptvsd0[2] = fCentrality;

		if(mixing_flag){
			if(trk->Charge()*casc->ChargeXi()<0){
				fHistoEleXiMassRSMix->Fill(cont);
				fHistoEleXiMassvsElePtRSMix->Fill(cont2);
				if(cont[0]<2.5){
					fHistoElePtRSMix->Fill(trk->Pt(),fCentrality);
					fHistoElePtvsEtaRSMix->Fill(cont_eleptvseta);
					fHistoElePtvsXiPtRSMix->Fill(cont_eleptvsxipt);
					fHistoElePtvsd0RSMix->Fill(cont_eleptvsd0);
				}
			}else{
				fHistoEleXiMassWSMix->Fill(cont);
				fHistoEleXiMassvsElePtWSMix->Fill(cont2);
				if(cont[0]<2.5){
					fHistoElePtWSMix->Fill(trk->Pt(),fCentrality);
					fHistoElePtvsEtaWSMix->Fill(cont_eleptvseta);
					fHistoElePtvsXiPtWSMix->Fill(cont_eleptvsxipt);
					fHistoElePtvsd0WSMix->Fill(cont_eleptvsd0);
				}
			}
		}else{
			if(trk->Charge()*casc->ChargeXi()<0){
				fHistoEleXiMassRS->Fill(cont);
				fHistoEleXiMassvsElePtRS->Fill(cont2);
				if(cont[0]<2.5){
					fHistoElePtRS->Fill(trk->Pt(),fCentrality);
					fHistoElePtvsEtaRS->Fill(cont_eleptvseta);
					fHistoElePtvsXiPtRS->Fill(cont_eleptvsxipt);
					fHistoElePtvsd0RS->Fill(cont_eleptvsd0);

					for(Int_t ih=0;ih<23;ih++){
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
							if(casc->ChargeXi()<0)
								cont_eleptvscutvars[1] = casc->MassLambda();
							else
								cont_eleptvscutvars[1] = casc->MassAntiLambda();
						}else if(ih==7){
							cont_eleptvscutvars[1] = casc->MassXi();
						}else if(ih==8){
							Double_t lPosV0[3];
							lPosV0[0] = casc->DecayVertexV0X();
							lPosV0[1] = casc->DecayVertexV0Y();
							lPosV0[2] = casc->DecayVertexV0Z();
							cont_eleptvscutvars[1] = TMath::Sqrt(lPosV0[0]*lPosV0[0]+lPosV0[1]*lPosV0[1]);
						}else if(ih==9){
							Double_t lPosXi[3];
							lPosXi[0] = casc->DecayVertexXiX();
							lPosXi[1] = casc->DecayVertexXiY();
							lPosXi[2] = casc->DecayVertexXiZ();
							cont_eleptvscutvars[1] = TMath::Sqrt(lPosXi[0]*lPosXi[0]+lPosXi[1]*lPosXi[1]);
						}else if(ih==10){
							cont_eleptvscutvars[1] = casc->DcaV0Daughters();
						}else if(ih==11){
							cont_eleptvscutvars[1] = casc->DcaXiDaughters();
						}else if(ih==12){
							cont_eleptvscutvars[1] = casc->DcaBachToPrimVertex();
						}else if(ih==13){
							if(casc->ChargeXi()<0.)
								cont_eleptvscutvars[1] = casc->DcaPosToPrimVertex();
							else
								cont_eleptvscutvars[1] = casc->DcaNegToPrimVertex();
						}else if(ih==14){
							if(casc->ChargeXi()>0.)
								cont_eleptvscutvars[1] = casc->DcaPosToPrimVertex();
							else
								cont_eleptvscutvars[1] = casc->DcaNegToPrimVertex();
						}else if(ih==15){
							cont_eleptvscutvars[1] =  casc->CosPointingAngle(casc->GetDecayVertexXi());
						}else if(ih==16){
							cont_eleptvscutvars[1] =  casc->CosPointingAngleXi(posVtx[0],posVtx[1],posVtx[2]);
						}else if(ih==17){
							cont_eleptvscutvars[1] =  nSigmaTPCv0pr;
						}else if(ih==18){
							cont_eleptvscutvars[1] =  nSigmaTPCv0pi;
						}else if(ih==19){
							cont_eleptvscutvars[1] =  nSigmaTPCbachpi;
						}else if(ih==20){
							cont_eleptvscutvars[1] =  casc->Eta();
						}else if(ih==21){
							cont_eleptvscutvars[1] =  0.5*TMath::Log((sqrt(casc->Ptot2Xi())+casc->MomXiZ())/(sqrt(casc->Ptot2Xi())-casc->MomXiZ()));
						}else if(ih==22){
							Double_t xipx = exobj->PxProng(1);
							Double_t xipy = exobj->PyProng(1);
							Double_t xipz = exobj->PzProng(1);
							Double_t epx = exobj->PxProng(0);
							Double_t epy = exobj->PyProng(0);
							Double_t epz = exobj->PzProng(0);
							cont_eleptvscutvars[1] = acos((xipx*epx+xipy*epy+xipz*epz)/sqrt(xipx*xipx+xipy*xipy+xipz*xipz)/sqrt(epx*epx+epy*epy+epz*epz));
						}else{
							cont_eleptvscutvars[1] = -9999.;
						}

						fHistoElePtvsCutVarsRS[ih]->Fill(cont_eleptvscutvars);
					}
				}
			}else{
				fHistoEleXiMassWS->Fill(cont);
				fHistoEleXiMassvsElePtWS->Fill(cont2);
				if(cont[0]<2.5){
					fHistoElePtWS->Fill(trk->Pt(),fCentrality);
					fHistoElePtvsEtaWS->Fill(cont_eleptvseta);
					fHistoElePtvsXiPtWS->Fill(cont_eleptvsxipt);
					fHistoElePtvsd0WS->Fill(cont_eleptvsd0);

					for(Int_t ih=0;ih<23;ih++){
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
							if(casc->ChargeXi()<0)
								cont_eleptvscutvars[1] = casc->MassLambda();
							else
								cont_eleptvscutvars[1] = casc->MassAntiLambda();
						}else if(ih==7){
							cont_eleptvscutvars[1] = casc->MassXi();
						}else if(ih==8){
							Double_t lPosV0[3];
							lPosV0[0] = casc->DecayVertexV0X();
							lPosV0[1] = casc->DecayVertexV0Y();
							lPosV0[2] = casc->DecayVertexV0Z();
							cont_eleptvscutvars[1] = TMath::Sqrt(lPosV0[0]*lPosV0[0]+lPosV0[1]*lPosV0[1]);
						}else if(ih==9){
							Double_t lPosXi[3];
							lPosXi[0] = casc->DecayVertexXiX();
							lPosXi[1] = casc->DecayVertexXiY();
							lPosXi[2] = casc->DecayVertexXiZ();
							cont_eleptvscutvars[1] = TMath::Sqrt(lPosXi[0]*lPosXi[0]+lPosXi[1]*lPosXi[1]);
						}else if(ih==10){
							cont_eleptvscutvars[1] = casc->DcaV0Daughters();
						}else if(ih==11){
							cont_eleptvscutvars[1] = casc->DcaXiDaughters();
						}else if(ih==12){
							cont_eleptvscutvars[1] = casc->DcaBachToPrimVertex();
						}else if(ih==13){
							if(casc->ChargeXi()<0.)
								cont_eleptvscutvars[1] = casc->DcaPosToPrimVertex();
							else
								cont_eleptvscutvars[1] = casc->DcaNegToPrimVertex();
						}else if(ih==14){
							if(casc->ChargeXi()>0.)
								cont_eleptvscutvars[1] = casc->DcaPosToPrimVertex();
							else
								cont_eleptvscutvars[1] = casc->DcaNegToPrimVertex();
						}else if(ih==15){
							cont_eleptvscutvars[1] =  casc->CosPointingAngle(casc->GetDecayVertexXi());
						}else if(ih==16){
							cont_eleptvscutvars[1] =  casc->CosPointingAngleXi(posVtx[0],posVtx[1],posVtx[2]);
						}else if(ih==17){
							cont_eleptvscutvars[1] =  nSigmaTPCv0pr;
						}else if(ih==18){
							cont_eleptvscutvars[1] =  nSigmaTPCv0pi;
						}else if(ih==19){
							cont_eleptvscutvars[1] =  nSigmaTPCbachpi;
						}else if(ih==20){
							cont_eleptvscutvars[1] =  casc->Eta();
						}else if(ih==21){
							cont_eleptvscutvars[1] =  0.5*TMath::Log((sqrt(casc->Ptot2Xi())+casc->MomXiZ())/(sqrt(casc->Ptot2Xi())-casc->MomXiZ()));
						}else if(ih==22){
							Double_t xipx = exobj->PxProng(1);
							Double_t xipy = exobj->PyProng(1);
							Double_t xipz = exobj->PzProng(1);
							Double_t epx = exobj->PxProng(0);
							Double_t epy = exobj->PyProng(0);
							Double_t epz = exobj->PzProng(0);
							cont_eleptvscutvars[1] = acos((xipx*epx+xipy*epy+xipz*epz)/sqrt(xipx*xipx+xipy*xipy+xipz*xipz)/sqrt(epx*epx+epy*epy+epz*epz));
						}else{
							cont_eleptvscutvars[1] = -9999.;
						}

						fHistoElePtvsCutVarsWS[ih]->Fill(cont_eleptvscutvars);
					}
				}
			}
		}

		if(fUseMCInfo){
			if(mcxic){
				Int_t pdgcode = mcxic->GetPdgCode();
				if(abs(pdgcode)==4132 && abs(mcpdgele_array[1])==4132 && abs(mcpdgcasc_array[1])==4132){
						fHistoEleXiMassMCS->Fill(cont);
						fHistoEleXiMassvsElePtMCS->Fill(cont2);
						if(cont[0]<2.5){
							fHistoElePtMCS->Fill(trk->Pt(),fCentrality);
							fHistoElePtvsEtaMCS->Fill(cont_eleptvseta);
							fHistoElePtvsXiPtMCS->Fill(cont_eleptvsxipt);
							fHistoElePtvsd0MCS->Fill(cont_eleptvsd0);

							Double_t cont_eleptvsxiptvsxicpt[4];
							cont_eleptvsxiptvsxicpt[0] = cont_eleptvsxipt[0];
							cont_eleptvsxiptvsxicpt[1] = cont_eleptvsxipt[1];
							cont_eleptvsxiptvsxicpt[2] = mcxic->Pt();
							cont_eleptvsxiptvsxicpt[3] = cont_eleptvsxipt[2];
							fHistoElePtvsXiPtvsXicPtMCS->Fill(cont_eleptvsxiptvsxicpt);

							Int_t labmotherxic = mcxic->GetMother();
							if(labmotherxic>=0){
								AliAODMCParticle *motherxic = (AliAODMCParticle*)mcArray->At(labmotherxic);
								Int_t pdgmotherxic = motherxic->GetPdgCode();
								if(TMath::Abs(pdgmotherxic)==511||TMath::Abs(pdgmotherxic)==521||TMath::Abs(pdgmotherxic)==5122||TMath::Abs(pdgmotherxic)==5132||TMath::Abs(pdgmotherxic)==5232||TMath::Abs(pdgmotherxic)==5332){
									fHistoElePtvsd0BFeeddownMCS->Fill(cont_eleptvsd0);
								}else{
									fHistoElePtvsd0PromptMCS->Fill(cont_eleptvsd0);
								}
							}else{
								fHistoElePtvsd0PromptMCS->Fill(cont_eleptvsd0);
							}

							for(Int_t ih=0;ih<23;ih++){
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
									if(casc->ChargeXi()<0)
										cont_eleptvscutvars[1] = casc->MassLambda();
									else
										cont_eleptvscutvars[1] = casc->MassAntiLambda();
								}else if(ih==7){
									cont_eleptvscutvars[1] = casc->MassXi();
								}else if(ih==8){
									Double_t lPosV0[3];
									lPosV0[0] = casc->DecayVertexV0X();
									lPosV0[1] = casc->DecayVertexV0Y();
									lPosV0[2] = casc->DecayVertexV0Z();
									cont_eleptvscutvars[1] = TMath::Sqrt(lPosV0[0]*lPosV0[0]+lPosV0[1]*lPosV0[1]);
								}else if(ih==9){
									Double_t lPosXi[3];
									lPosXi[0] = casc->DecayVertexXiX();
									lPosXi[1] = casc->DecayVertexXiY();
									lPosXi[2] = casc->DecayVertexXiZ();
									cont_eleptvscutvars[1] = TMath::Sqrt(lPosXi[0]*lPosXi[0]+lPosXi[1]*lPosXi[1]);
								}else if(ih==10){
									cont_eleptvscutvars[1] = casc->DcaV0Daughters();
								}else if(ih==11){
									cont_eleptvscutvars[1] = casc->DcaXiDaughters();
								}else if(ih==12){
									cont_eleptvscutvars[1] = casc->DcaBachToPrimVertex();
								}else if(ih==13){
									if(casc->ChargeXi()<0.)
										cont_eleptvscutvars[1] = casc->DcaPosToPrimVertex();
									else
										cont_eleptvscutvars[1] = casc->DcaNegToPrimVertex();
								}else if(ih==14){
									if(casc->ChargeXi()>0.)
										cont_eleptvscutvars[1] = casc->DcaPosToPrimVertex();
									else
										cont_eleptvscutvars[1] = casc->DcaNegToPrimVertex();
								}else if(ih==15){
									cont_eleptvscutvars[1] =  casc->CosPointingAngle(casc->GetDecayVertexXi());
								}else if(ih==16){
									cont_eleptvscutvars[1] =  casc->CosPointingAngleXi(posVtx[0],posVtx[1],posVtx[2]);
								}else if(ih==17){
									cont_eleptvscutvars[1] =  nSigmaTPCv0pr;
								}else if(ih==18){
									cont_eleptvscutvars[1] =  nSigmaTPCv0pi;
								}else if(ih==19){
									cont_eleptvscutvars[1] =  nSigmaTPCbachpi;
								}else if(ih==20){
									cont_eleptvscutvars[1] =  casc->Eta();
								}else if(ih==21){
									cont_eleptvscutvars[1] =  0.5*TMath::Log((sqrt(casc->Ptot2Xi())+casc->MomXiZ())/(sqrt(casc->Ptot2Xi())-casc->MomXiZ()));
								}else if(ih==22){
									Double_t xipx = exobj->PxProng(1);
									Double_t xipy = exobj->PyProng(1);
									Double_t xipz = exobj->PzProng(1);
									Double_t epx = exobj->PxProng(0);
									Double_t epy = exobj->PyProng(0);
									Double_t epz = exobj->PzProng(0);
									cont_eleptvscutvars[1] = acos((xipx*epx+xipy*epy+xipz*epz)/sqrt(xipx*xipx+xipy*xipy+xipz*xipz)/sqrt(epx*epx+epy*epy+epz*epz));
								}else{
									cont_eleptvscutvars[1] = -9999.;
								}

								fHistoElePtvsCutVarsMCS[ih]->Fill(cont_eleptvscutvars);
							}
						}
				}
			}
		}
	}

  return;
}

////-------------------------------------------------------------------------------
void AliAnalysisTaskSEXic2eleXifromAODtracks::DefineEleTreeVariables() 
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
void AliAnalysisTaskSEXic2eleXifromAODtracks::FillElectronROOTObjects(AliAODTrack *trk, TClonesArray *mcArray) 
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
void AliAnalysisTaskSEXic2eleXifromAODtracks::DefineCascTreeVariables() 
{
  //
  // Define V0 tree variables
  //

  const char* nameoutput = GetOutputSlot(6)->GetContainer()->GetName();
  fCascVariablesTree = new TTree(nameoutput,"cascade variables tree");
  Int_t nVar = 17;
  fCandidateCascVariables = new Float_t [nVar];
  TString * fCandidateVariableNames = new TString[nVar];

  fCandidateVariableNames[ 0]="Centrality";
  fCandidateVariableNames[ 1]="InvMassXi";
  fCandidateVariableNames[ 2]="XiPx";
  fCandidateVariableNames[ 3]="XiPy";
  fCandidateVariableNames[ 4]="XiPz";
  fCandidateVariableNames[ 5]="InvMassLambda";
  fCandidateVariableNames[ 6]="DcaXiDaughters";
  fCandidateVariableNames[ 7]="DcaV0Daughters";
  fCandidateVariableNames[ 8]="DecayLengthXi";
  fCandidateVariableNames[ 9]="CosPointingAngleXi";
  fCandidateVariableNames[10]="DcaV0toPrimVertex";
  fCandidateVariableNames[11]="DcaPostoPrimVertex";
  fCandidateVariableNames[12]="DcaNegtoPrimVertex";
  fCandidateVariableNames[13]="DcaBachtoPrimVertex";
  fCandidateVariableNames[14]="DecayLengthV0";
  fCandidateVariableNames[15]="CosPointingAngleV0";
  fCandidateVariableNames[16]="XiCharge";

  for (Int_t ivar=0; ivar<nVar; ivar++) {
    fCascVariablesTree->Branch(fCandidateVariableNames[ivar].Data(),&fCandidateCascVariables[ivar],Form("%s/f",fCandidateVariableNames[ivar].Data()));
  }

  return;
}

////-------------------------------------------------------------------------------
void AliAnalysisTaskSEXic2eleXifromAODtracks::FillCascROOTObjects(AliAODcascade *casc, TClonesArray *mcArray) 
{
  //
  // Fill histograms or tree depending on fWriteVariableTree (tree not implemented yet)
  //
	if(!casc) return;
	fHistoXiMassvsPt->Fill(casc->MassXi(),sqrt(casc->MomXiX()*casc->MomXiX()+casc->MomXiY()*casc->MomXiY()));
	fHistoOmegaMassvsPt->Fill(casc->MassOmega(),sqrt(casc->MomXiX()*casc->MomXiX()+casc->MomXiY()*casc->MomXiY()));
	Double_t momxix = casc->MomXiX();
	Double_t momxiy = casc->MomXiY();
	Double_t phi_alice = atan2(momxiy,momxix);
	if(phi_alice<0.) phi_alice += 2 * M_PI;
	fHistoXiQovPtvsPhi->Fill(phi_alice,(Double_t)casc->ChargeXi()/sqrt(momxix*momxix+momxiy*momxiy));

	if(fUseMCInfo){
		Int_t pdgDgcasc[2]={211,3122};
		Int_t pdgDgv0[2]={2212,211};
		Int_t labcasc = MatchToMCCascade(casc,3312,pdgDgcasc,pdgDgv0,mcArray); // the cascade
		if(labcasc>=0){
			fHistoXiMassvsPtMCS->Fill(casc->MassXi(),sqrt(casc->MomXiX()*casc->MomXiX()+casc->MomXiY()*casc->MomXiY()));
		}
	}

	if(!fWriteEachVariableTree) return;


	for(Int_t i=0;i<16;i++){
		fCandidateCascVariables[i] = -9999.;
	}
  Double_t posVtx[3] = {0.,0.,0.};
  fVtx1->GetXYZ(posVtx);

  fCandidateCascVariables[ 0] = fCentrality;
  fCandidateCascVariables[ 1] = casc->MassXi();
  fCandidateCascVariables[ 2] = casc->MomXiX();
  fCandidateCascVariables[ 3] = casc->MomXiY();
  fCandidateCascVariables[ 4] = casc->MomXiZ();
	if(casc->ChargeXi()<0)
		fCandidateCascVariables[ 5] = casc->MassLambda();
	else
		fCandidateCascVariables[ 5] = casc->MassAntiLambda();

  fCandidateCascVariables[ 6] = casc->DcaXiDaughters();
  fCandidateCascVariables[ 7] = casc->DcaV0Daughters();
  fCandidateCascVariables[ 8] = casc->DecayLengthXi(posVtx[0],posVtx[1],posVtx[2]);
  fCandidateCascVariables[ 9] = casc->CosPointingAngleXi(posVtx[0],posVtx[1],posVtx[2]);
  fCandidateCascVariables[10] = casc->DcaV0ToPrimVertex();
  fCandidateCascVariables[11] = casc->DcaPosToPrimVertex();
  fCandidateCascVariables[12] = casc->DcaNegToPrimVertex();
  fCandidateCascVariables[13] = casc->DcaBachToPrimVertex();
  fCandidateCascVariables[14] = casc->DecayLengthV0();
  fCandidateCascVariables[15] = casc->CosPointingAngle(casc->GetDecayVertexXi());
  fCandidateCascVariables[16] = casc->ChargeXi();


	fCascVariablesTree->Fill();
}
////-------------------------------------------------------------------------------
void AliAnalysisTaskSEXic2eleXifromAODtracks::DefineMCTreeVariables() 
{
  //
  // Define electron tree variables
  //

  const char* nameoutput = GetOutputSlot(7)->GetContainer()->GetName();
  fMCVariablesTree = new TTree(nameoutput,"MC variables tree");
  Int_t nVar = 14;
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
  fCandidateVariableNames[ 8]="CascPx";
  fCandidateVariableNames[ 9]="CascPy";
  fCandidateVariableNames[10]="CascPz";
  fCandidateVariableNames[11]="PdgCode";
  fCandidateVariableNames[12]="ElePdgCode";
  fCandidateVariableNames[13]="CascPdgCode";

  for (Int_t ivar=0; ivar<nVar; ivar++) {
    fMCVariablesTree->Branch(fCandidateVariableNames[ivar].Data(),&fCandidateMCVariables[ivar],Form("%s/f",fCandidateVariableNames[ivar].Data()));
  }
  return;
}
////-------------------------------------------------------------------------------
void AliAnalysisTaskSEXic2eleXifromAODtracks::FillMCROOTObjects(AliAODMCParticle *mcpart, AliAODMCParticle *mcepart, AliAODMCParticle *mccascpart, Int_t decaytype) 
{
  //
  // Fill histograms or tree depending on fWriteMCVariableTree 
  //
	if(!mcpart) return;
	if(!mcepart) return;
	if(!mccascpart) return;

	for(Int_t i=0;i<14;i++){
		fCandidateMCVariables[i] = -9999.;
	}

	fCandidateMCVariables[ 0] = fCentrality;
	fCandidateMCVariables[ 1] = decaytype;
	fCandidateMCVariables[ 2] = mcpart->Px();
	fCandidateMCVariables[ 3] = mcpart->Py();
	fCandidateMCVariables[ 4] = mcpart->Pz();
	fCandidateMCVariables[ 5] = mcepart->Px();
	fCandidateMCVariables[ 6] = mcepart->Py();
	fCandidateMCVariables[ 7] = mcepart->Pz();
	fCandidateMCVariables[ 8] = mccascpart->Px();
	fCandidateMCVariables[ 9] = mccascpart->Py();
	fCandidateMCVariables[10] = mccascpart->Pz();
	fCandidateMCVariables[11] = mcpart->GetPdgCode();
	fCandidateMCVariables[12] = mcepart->GetPdgCode();
	fCandidateMCVariables[13] = mccascpart->GetPdgCode();

	Double_t epx = mcepart->Px();
	Double_t epy = mcepart->Py();
	Double_t epz = mcepart->Pz();
	Double_t eE = sqrt(epx*epx+epy*epy+epz*epz+0.000511*0.000511);
	Double_t cascpx = mccascpart->Px();
	Double_t cascpy = mccascpart->Py();
	Double_t cascpz = mccascpart->Pz();
	Double_t cascE = sqrt(cascpx*cascpx+cascpy*cascpy+cascpz*cascpz+1.32171*1.32171);

	Double_t InvMassEleXi = sqrt(pow(eE+cascE,2)-pow(epx+cascpx,2)-pow(epy+cascpy,2)-pow(epz+cascpz,2));

	Double_t cont[3];
	cont[0] = InvMassEleXi;
	cont[1] = mcpart->Pt();
	cont[2] = fCentrality;
	Double_t cont2[3];
	cont2[0] = InvMassEleXi;
	cont2[1] = mcepart->Pt();
	cont2[2] = fCentrality;
	Double_t cont_eleptvseta[3];
	cont_eleptvseta[0] = mcepart->Pt();
	cont_eleptvseta[1] = mcepart->Eta();
	cont_eleptvseta[2] = fCentrality;
	Double_t cont_eleptvsxipt[3];
	cont_eleptvsxipt[0] = mcepart->Pt();
	cont_eleptvsxipt[1] = mccascpart->Pt();
	cont_eleptvsxipt[2] = fCentrality;
	Double_t cont_eleptvsxiptvsxicpt[4];
	cont_eleptvsxiptvsxicpt[0] = mcepart->Pt();
	cont_eleptvsxiptvsxicpt[1] = mccascpart->Pt();
	cont_eleptvsxiptvsxicpt[2] = mcpart->Pt();
	cont_eleptvsxiptvsxicpt[3] = fCentrality;

	AliESDtrackCuts *esdcuts = fAnalCuts->GetTrackCuts();
	Float_t etamin, etamax;
	esdcuts->GetEtaRange(etamin,etamax);

	if(decaytype==0){
		fHistoEleXiMassMCGen->Fill(cont);
		if(fabs(mcepart->Eta())<etamax){
			fHistoEleXiMassvsElePtMCGen->Fill(cont2);
			if(InvMassEleXi<2.5){
				fHistoElePtMCGen->Fill(mcepart->Pt(),fCentrality);
				fHistoElePtvsEtaMCGen->Fill(cont_eleptvseta);
				fHistoElePtvsXiPtMCGen->Fill(cont_eleptvsxipt);
			}
		}
		if(fabs(mcpart->Y())<0.7){
			if(InvMassEleXi<2.5){
				fHistoElePtvsXiPtMCXicGen->Fill(cont_eleptvsxipt);
				fHistoElePtvsXiPtvsXicPtMCGen->Fill(cont_eleptvsxiptvsxicpt);
			}
		}
	}

	if(fWriteMCVariableTree)
		fMCVariablesTree->Fill();
}

////__________________________________________________________________________
void  AliAnalysisTaskSEXic2eleXifromAODtracks::DefineGeneralHistograms() {
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
void  AliAnalysisTaskSEXic2eleXifromAODtracks::DefineAnalysisHistograms() 
{
  //
  // Define analyis histograms
  //
	
  //------------------------------------------------
  // Basic histogram
  //------------------------------------------------
  Int_t bins_base[3]=		{16	,100		,10};
  Double_t xmin_base[3]={1.3,0		,0.00};
  Double_t xmax_base[3]={3.3,10.	,100};

  fHistoEleXiMass = new THnSparseF("fHistoEleXiMass","",3,bins_base,xmin_base,xmax_base);
  fOutputAll->Add(fHistoEleXiMass);
  fHistoEleXiMassRS = new THnSparseF("fHistoEleXiMassRS","",3,bins_base,xmin_base,xmax_base);
  fOutputAll->Add(fHistoEleXiMassRS);
  fHistoEleXiMassWS = new THnSparseF("fHistoEleXiMassWS","",3,bins_base,xmin_base,xmax_base);
  fOutputAll->Add(fHistoEleXiMassWS);
  fHistoEleXiMassRSMix = new THnSparseF("fHistoEleXiMassRSMix","",3,bins_base,xmin_base,xmax_base);
  fOutputAll->Add(fHistoEleXiMassRSMix);
  fHistoEleXiMassWSMix = new THnSparseF("fHistoEleXiMassWSMix","",3,bins_base,xmin_base,xmax_base);
  fOutputAll->Add(fHistoEleXiMassWSMix);

  fHistoEleXiMassvsElePtRS = new THnSparseF("fHistoEleXiMassvsElePtRS","",3,bins_base,xmin_base,xmax_base);
  fOutputAll->Add(fHistoEleXiMassvsElePtRS);
  fHistoEleXiMassvsElePtWS = new THnSparseF("fHistoEleXiMassvsElePtWS","",3,bins_base,xmin_base,xmax_base);
  fOutputAll->Add(fHistoEleXiMassvsElePtWS);
  fHistoEleXiMassvsElePtRSMix = new THnSparseF("fHistoEleXiMassvsElePtRSMix","",3,bins_base,xmin_base,xmax_base);
  fOutputAll->Add(fHistoEleXiMassvsElePtRSMix);
  fHistoEleXiMassvsElePtWSMix = new THnSparseF("fHistoEleXiMassvsElePtWSMix","",3,bins_base,xmin_base,xmax_base);
  fOutputAll->Add(fHistoEleXiMassvsElePtWSMix);

  fHistoElePtRS=new TH2F("fHistoElePtRS","Right-sign e pt",100,0.,10.,10,0.,100.);
  fOutputAll->Add(fHistoElePtRS);
  fHistoElePtWS=new TH2F("fHistoElePtWS","Wrong-sign e pt",100,0.,10.,10,0.,100.);
  fOutputAll->Add(fHistoElePtWS);
  fHistoElePtRSMix=new TH2F("fHistoElePtRSMix","Right-sign e pt",100,0.,10.,10,0.,100.);
  fOutputAll->Add(fHistoElePtRSMix);
  fHistoElePtWSMix=new TH2F("fHistoElePtWSMix","Wrong-sign e pt",100,0.,10.,10,0.,100.);
  fOutputAll->Add(fHistoElePtWSMix);

  fHistoEleXiMassMCS = new THnSparseF("fHistoEleXiMassMCS","",3,bins_base,xmin_base,xmax_base);
  fOutputAll->Add(fHistoEleXiMassMCS);
  fHistoEleXiMassMCGen = new THnSparseF("fHistoEleXiMassMCGen","",3,bins_base,xmin_base,xmax_base);
  fOutputAll->Add(fHistoEleXiMassMCGen);
  fHistoEleXiMassvsElePtMCS = new THnSparseF("fHistoEleXiMassvsElePtMCS","",3,bins_base,xmin_base,xmax_base);
  fOutputAll->Add(fHistoEleXiMassvsElePtMCS);
  fHistoEleXiMassvsElePtMCGen = new THnSparseF("fHistoEleXiMassvsElePtMCGen","",3,bins_base,xmin_base,xmax_base);
  fOutputAll->Add(fHistoEleXiMassvsElePtMCGen);
  fHistoElePtMCS=new TH2F("fHistoElePtMCS","MC S e pt",100,0.,10.,10,0.,100.);
  fOutputAll->Add(fHistoElePtMCS);
  fHistoElePtMCGen=new TH2F("fHistoElePtMCGen","MC Gen e pt",100,0.,10.,10,0.,100.);
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

  Int_t bins_eleptvsxipt[3]=	{50,20	,10};
  Double_t xmin_eleptvsxipt[3]={0.,0.	,0.0};
  Double_t xmax_eleptvsxipt[3]={5.,5.	,100};

  fHistoElePtvsXiPtRS = new THnSparseF("fHistoElePtvsXiPtRS","",3,bins_eleptvsxipt,xmin_eleptvsxipt,xmax_eleptvsxipt);
  fOutputAll->Add(fHistoElePtvsXiPtRS);
  fHistoElePtvsXiPtWS = new THnSparseF("fHistoElePtvsXiPtWS","",3,bins_eleptvsxipt,xmin_eleptvsxipt,xmax_eleptvsxipt);
  fOutputAll->Add(fHistoElePtvsXiPtWS);
  fHistoElePtvsXiPtRSMix = new THnSparseF("fHistoElePtvsXiPtRSMix","",3,bins_eleptvsxipt,xmin_eleptvsxipt,xmax_eleptvsxipt);
  fOutputAll->Add(fHistoElePtvsXiPtRSMix);
  fHistoElePtvsXiPtWSMix = new THnSparseF("fHistoElePtvsXiPtWSMix","",3,bins_eleptvsxipt,xmin_eleptvsxipt,xmax_eleptvsxipt);
  fOutputAll->Add(fHistoElePtvsXiPtWSMix);
  fHistoElePtvsXiPtMCS = new THnSparseF("fHistoElePtvsXiPtMCS","",3,bins_eleptvsxipt,xmin_eleptvsxipt,xmax_eleptvsxipt);
  fOutputAll->Add(fHistoElePtvsXiPtMCS);
  fHistoElePtvsXiPtMCGen = new THnSparseF("fHistoElePtvsXiPtMCGen","",3,bins_eleptvsxipt,xmin_eleptvsxipt,xmax_eleptvsxipt);
  fOutputAll->Add(fHistoElePtvsXiPtMCGen);
  fHistoElePtvsXiPtMCXicGen = new THnSparseF("fHistoElePtvsXiPtMCXicGen","",3,bins_eleptvsxipt,xmin_eleptvsxipt,xmax_eleptvsxipt);
  fOutputAll->Add(fHistoElePtvsXiPtMCXicGen);

  Int_t bins_eleptvsxiptvsxicpt[4]=	{50,20,10,10};
  Double_t xmin_eleptvsxiptvsxicpt[4]={0.,0.,0.,0.0};
  Double_t xmax_eleptvsxiptvsxicpt[4]={5.,5.,10.,100};
  fHistoElePtvsXiPtvsXicPtMCS = new THnSparseF("fHistoElePtvsXiPtvsXicPtMCS","",4,bins_eleptvsxiptvsxicpt,xmin_eleptvsxiptvsxicpt,xmax_eleptvsxiptvsxicpt);
  fOutputAll->Add(fHistoElePtvsXiPtvsXicPtMCS);
  fHistoElePtvsXiPtvsXicPtMCGen = new THnSparseF("fHistoElePtvsXiPtvsXicPtMCGen","",4,bins_eleptvsxiptvsxicpt,xmin_eleptvsxiptvsxicpt,xmax_eleptvsxiptvsxicpt);
  fOutputAll->Add(fHistoElePtvsXiPtvsXicPtMCGen);

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
  fHistoXiMassvsPt=new TH2F("fHistoXiMassvsPt","Xi mass",100,1.32-0.05,1.32+0.05,20,0.,10.);
  fOutputAll->Add(fHistoXiMassvsPt);
  fHistoXiMassvsPtMCS=new TH2F("fHistoXiMassvsPtMCS","Xi mass",100,1.32-0.05,1.32+0.05,20,0.,10.);
  fOutputAll->Add(fHistoXiMassvsPtMCS);
  fHistoXiMassvsPtMCGen=new TH2F("fHistoXiMassvsPtMCGen","Xi mass",100,1.32-0.05,1.32+0.05,20,0.,10.);
  fOutputAll->Add(fHistoXiMassvsPtMCGen);
  fHistoOmegaMassvsPt=new TH2F("fHistoOmegaMassvsPt","Omega mass",100,1.67-0.05,1.67+0.05,20,0.,10.);
  fOutputAll->Add(fHistoOmegaMassvsPt);

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
  fHistoXiQovPtvsPhi=new TH2F("fHistoXiQovPtvsPhi","",70,0.,7.,50,-2.,2.);
  fOutputAll->Add(fHistoXiQovPtvsPhi);

  fHistonEvtvsRunNumber=new TH1F("fHistonEvtvsRunNumber","",20000,-0.5,19999.5);
  fOutputAll->Add(fHistonEvtvsRunNumber);
  fHistonElevsRunNumber=new TH1F("fHistonElevsRunNumber","",20000,-0.5,19999.5);
  fOutputAll->Add(fHistonElevsRunNumber);
  fHistonXivsRunNumber=new TH1F("fHistonXivsRunNumber","",20000,-0.5,19999.5);
  fOutputAll->Add(fHistonXivsRunNumber);

	for(Int_t ih=0;ih<23;ih++){
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
			//7: Xi mass
			bins_eleptvscutvars[1] = 50;
			xmin_eleptvscutvars[1] = 1.32-0.03;
			xmax_eleptvscutvars[1] = 1.32+0.03;
		}else if(ih==8 || ih==9){
			//8: Rfid Lambda, 9: Rfid Xi
			bins_eleptvscutvars[1] = 20;
			xmin_eleptvscutvars[1] = 0.;
			xmax_eleptvscutvars[1] = 5.;
		}else if(ih==10 || ih==11){
			//11: DCA Xi, 10: Dca V0
			bins_eleptvscutvars[1] = 20;
			xmin_eleptvscutvars[1] = 0.;
			xmax_eleptvscutvars[1] = 2.;
		}else if(ih==12 || ih==13 || ih==14){
			//12: DCA Bachto prim, 13: DCA V0pr to prim 14: DCA V0pi to prim
			bins_eleptvscutvars[1] = 20;
			xmin_eleptvscutvars[1] = 0.;
			xmax_eleptvscutvars[1] = 0.5;
		}else if(ih==15 || ih==16){
			//16: CosPAXi, 15: CosPAv0
			bins_eleptvscutvars[1] = 20;
			xmin_eleptvscutvars[1] = 0.95;
			xmax_eleptvscutvars[1] = 1.0;
		}else if(ih==17 || ih==18 || ih==19){
			//17: nSigma(TPC, bach)  18: nSigma(TPC, v0pr), 19: nSigma(TPC,v0pi)
			bins_eleptvscutvars[1] = 20;
			xmin_eleptvscutvars[1] = -5.;
			xmax_eleptvscutvars[1] = 5.;
		}else if(ih==20 || ih==21){
			//20: V0 eta 21:  Xi eta
			bins_eleptvscutvars[1] = 30;
			xmin_eleptvscutvars[1] = -1.5;
			xmax_eleptvscutvars[1] = 1.5;
		}else if(ih==22){
			//20: Opening angle
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
AliAODRecoCascadeHF* AliAnalysisTaskSEXic2eleXifromAODtracks::MakeCascadeHF(AliAODcascade *casc, AliAODTrack *part, AliAODEvent * aod, AliAODVertex *secVert, Bool_t mixing) 
{
  //
  // Create AliAODRecoCascadeHF object from the argument
  //

  if(!casc) return 0x0;
  if(!part) return 0x0;
  if(!aod) return 0x0;

  //------------------------------------------------
  // PrimaryVertex
  //------------------------------------------------
  AliAODVertex *primVertexAOD;
  Bool_t unsetvtx = kFALSE;
  if(fReconstructPrimVert){
    primVertexAOD = CallPrimaryVertex(casc,part,aod);
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

  Double_t xyz[3], pxpypz[3], cv[21]; Short_t sign;
  xyz[0]=casc->DecayVertexXiX();
  xyz[1]=casc->DecayVertexXiY();
  xyz[2]=casc->DecayVertexXiZ();
  pxpypz[0]=casc->MomXiX();
  pxpypz[1]=casc->MomXiY();
  pxpypz[2]=casc->MomXiZ();
  casc->GetCovarianceXYZPxPyPz(cv);
  sign=casc->ChargeXi();
  AliExternalTrackParam	*trackCasc = new AliExternalTrackParam(xyz,pxpypz,cv,sign);

  Double_t xdummy, ydummy;
  Double_t dca = esdtrack->GetDCA(trackCasc,fBzkG,xdummy,ydummy);


  //------------------------------------------------
  // Propagate all tracks to the secondary vertex and calculate momentum there
  //------------------------------------------------
	
  Double_t d0z0bach[2],covd0z0bach[3];
  if(sqrt(pow(secVert->GetX(),2)+pow(secVert->GetY(),2))<1.){
    part->PropagateToDCA(secVert,fBzkG,kVeryBig,d0z0bach,covd0z0bach);
    trackCasc->PropagateToDCA(secVert,fBzkG,kVeryBig);
  }else{
    part->PropagateToDCA(primVertexAOD,fBzkG,kVeryBig,d0z0bach,covd0z0bach);
    trackCasc->PropagateToDCA(primVertexAOD,fBzkG,kVeryBig);
  }
  Double_t momcasc_new[3]={-9999,-9999,-9999.};
  trackCasc->GetPxPyPz(momcasc_new);

  Double_t px[2],py[2],pz[2];
  px[0] = part->Px(); py[0] = part->Py(); pz[0] = part->Pz(); 
  px[1] = momcasc_new[0]; py[1] = momcasc_new[1]; pz[1] = momcasc_new[2]; 

  //------------------------------------------------
  // d0
  //------------------------------------------------
  Double_t d0[3],d0err[3];

  part->PropagateToDCA(primVertexAOD,fBzkG,kVeryBig,d0z0bach,covd0z0bach);
  d0[0]= d0z0bach[0];
  d0err[0] = TMath::Sqrt(covd0z0bach[0]);

  Double_t d0z0casc[2],covd0z0casc[3];
  trackCasc->PropagateToDCA(primVertexAOD,fBzkG,kVeryBig,d0z0casc,covd0z0casc);
  d0[1]= d0z0casc[0];
  d0err[1] = TMath::Sqrt(covd0z0casc[0]);

  //------------------------------------------------
  // Create AliAODRecoCascadeHF
  //------------------------------------------------
  Short_t charge = part->Charge();
  AliAODRecoCascadeHF *theCascade = new AliAODRecoCascadeHF(secVert,charge,px,py,pz,d0,d0err,dca);
  if(!theCascade)  
    {
      if(unsetvtx) delete primVertexAOD; primVertexAOD=NULL;
      if(esdtrack) delete esdtrack;
      if(trackCasc) delete trackCasc;
      return 0x0;
    }
  theCascade->SetOwnPrimaryVtx(primVertexAOD);
  UShort_t id[2]={(UShort_t)part->GetID(),(UShort_t)trackCasc->GetID()};
  theCascade->SetProngIDs(2,id);

	if(!mixing){
		theCascade->GetSecondaryVtx()->AddDaughter(part);
		theCascade->GetSecondaryVtx()->AddDaughter(casc);
	}

  if(unsetvtx) delete primVertexAOD; primVertexAOD=NULL;
  if(esdtrack) delete esdtrack;
  if(trackCasc) delete trackCasc;

  return theCascade;
}

//________________________________________________________________________
AliAODVertex* AliAnalysisTaskSEXic2eleXifromAODtracks::CallPrimaryVertex(AliAODcascade *casc, AliAODTrack *trk, AliAODEvent* aod)
{
  //
  // Make an array of tracks which should not be used in primary vertex calculation and 
  // Call PrimaryVertex function
  //

  TObjArray *TrackArray = new TObjArray(3);
  
  AliESDtrack *cptrk1 = new AliESDtrack((AliVTrack*)trk);
  TrackArray->AddAt(cptrk1,0);
  
  AliESDtrack *cascptrack = new AliESDtrack((AliVTrack*)casc->GetDaughter(0));
  TrackArray->AddAt(cascptrack,1);
  AliESDtrack *cascntrack = new AliESDtrack((AliVTrack*)casc->GetDaughter(1));
  TrackArray->AddAt(cascntrack,2);
  AliESDtrack *cascbtrack = new AliESDtrack((AliVTrack*)casc->GetDecayVertexXi()->GetDaughter(0));
  TrackArray->AddAt(cascbtrack,3);
  
  AliAODVertex *newvert  = PrimaryVertex(TrackArray,aod);
  
  for(Int_t i=0;i<4;i++)
    {
      AliESDtrack *tesd = (AliESDtrack*)TrackArray->UncheckedAt(i);
      delete tesd;
    }
  TrackArray->Clear();
  delete TrackArray;
  
  return newvert;
}

//________________________________________________________________________
AliAODVertex* AliAnalysisTaskSEXic2eleXifromAODtracks::PrimaryVertex(const TObjArray *trkArray,
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
AliAODVertex* AliAnalysisTaskSEXic2eleXifromAODtracks::ReconstructSecondaryVertex(AliAODcascade *casc, AliAODTrack *part, AliAODEvent * aod) 
{
  //
  // Reconstruct secondary vertex from trkArray (Copied from AliAnalysisVertexingHF)
  //
	
  AliAODVertex *primVertexAOD;
  Bool_t unsetvtx = kFALSE;
  if(fReconstructPrimVert){
    primVertexAOD = CallPrimaryVertex(casc,part,aod);
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
Int_t AliAnalysisTaskSEXic2eleXifromAODtracks::MatchToMC(AliAODRecoCascadeHF *exobj, TClonesArray *mcArray, Int_t *pdgarray_ele, Int_t *pdgarray_casc, Int_t *labelarray_ele, Int_t *labelarray_casc,  Int_t &ngen_ele, Int_t &ngen_casc) 
{
  //
  // Match to MC
  //
	for(Int_t i=0;i<100;i++){
		pdgarray_ele[i] = -9999;
		labelarray_ele[i] = -9999;
		pdgarray_casc[i] = -9999;
		labelarray_casc[i] = -9999;
	}
	ngen_ele = 0;
	ngen_casc = 0;

  AliVTrack *trk = dynamic_cast<AliVTrack*>(exobj->GetBachelor());
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

  AliAODcascade *theCascade = dynamic_cast<AliAODcascade*>(exobj->GetCascade());
	if(!theCascade) return -1;

	Int_t pdgDgcasc[2]={211,3122};
	Int_t pdgDgv0[2]={2212,211};
  Int_t labcasc = MatchToMCCascade(theCascade,3312,pdgDgcasc,pdgDgv0,mcArray); // the cascade
  if(labcasc<0) return -1;

	AliAODMCParticle *mccasc = (AliAODMCParticle*)mcArray->At(labcasc);
	if(!mccasc) return -1;
	labelarray_casc[0] = labcasc;
	pdgarray_casc[0] = mccasc->GetPdgCode();
	ngen_casc ++;

  AliAODMCParticle *mcprimcasc=0;
  mcprimcasc = mccasc;
  while(mcprimcasc->GetMother()>=0) {
    Int_t labprim_casc=mcprimcasc->GetMother();
    AliAODMCParticle *tmcprimcasc = (AliAODMCParticle*)mcArray->At(labprim_casc);
    if(!tmcprimcasc) {
			break;
    }

    mcprimcasc = tmcprimcasc;
		pdgarray_casc[ngen_casc] = mcprimcasc->GetPdgCode();
		labelarray_casc[ngen_casc] = labprim_casc;
		ngen_casc ++;
		if(ngen_casc==100) break;
  }

	Bool_t same_flag = kFALSE;
	Int_t matchedlabel=-9999;
	for(Int_t iemc=0;iemc<ngen_ele;iemc++){
		for(Int_t ivmc=0;ivmc<ngen_casc;ivmc++){
			if(labelarray_ele[iemc]==labelarray_casc[ivmc]){
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
Int_t AliAnalysisTaskSEXic2eleXifromAODtracks::MatchToMCCascade(AliAODcascade *theCascade, Int_t pdgabscasc, Int_t *pdgDgcasc, Int_t *pdgDgv0, TClonesArray *mcArray) const // the cascade
{
	//
	// Matching to MC of cascade
	//

	AliAODTrack *cptrack = (AliAODTrack*) theCascade->GetDaughter(0);
	if(!cptrack) return -1;
	Int_t label_p = cptrack->GetLabel();
	if(label_p<0) return -1;
	AliAODTrack *cntrack = (AliAODTrack*) theCascade->GetDaughter(1);
	if(!cntrack) return -1;
	Int_t label_n = cntrack->GetLabel();
	if(label_n<0) return -1;
	Int_t labv0 = theCascade->MatchToMC(pdgDgcasc[1],mcArray,2,pdgDgv0);
	if(labv0<0) return -1;
	AliAODMCParticle *mcpartv0= (AliAODMCParticle*) mcArray->At(labv0);

	AliAODTrack *cbtrack = (AliAODTrack*) theCascade->GetDecayVertexXi()->GetDaughter(0);
	if(!cbtrack) return -1;

	Int_t label_b = cbtrack->GetLabel();
	if(label_b<0) return -1;

	AliAODMCParticle *mcpartb= (AliAODMCParticle*) mcArray->At(label_b);
	Int_t pdgb = TMath::Abs(mcpartb->GetPdgCode());
	if(pdgb!=pdgDgcasc[0]) return -1;

	AliAODMCParticle *mcmotherv0=mcpartv0;
	Bool_t isFromXiv0 = kFALSE;
	Int_t labxiv0 = mcmotherv0->GetMother();
	if(labxiv0<0) return -1;
	mcmotherv0 =  (AliAODMCParticle*) mcArray->At(labxiv0);
	if(mcmotherv0){
		Int_t pdg = TMath::Abs(mcmotherv0 ->GetPdgCode());
		if(pdg==pdgabscasc){
			isFromXiv0 = kTRUE;
		}
	}
	if(!isFromXiv0) return -1;

	AliAODMCParticle *mcmotherb=mcpartb;
	Bool_t isFromXib = kFALSE;
	Int_t labxib = mcmotherb->GetMother();
	if(labxib<0) return -1;
	mcmotherb =  (AliAODMCParticle*) mcArray->At(labxib);
	if(mcmotherb){
		Int_t pdg = TMath::Abs(mcmotherb ->GetPdgCode());
		if(pdg==pdgabscasc){
			isFromXib = kTRUE;
		}
	}
	if(!isFromXib) return -1;

	if(labxiv0!=labxib) return -1;//Bachelor and V0 should come from the same Xi

	return labxib;
}
//________________________________________________________________________
void AliAnalysisTaskSEXic2eleXifromAODtracks::SelectTrack( const AliVEvent *event, Int_t trkEntries, Int_t &nSeleTrks,Bool_t *seleFlags, TClonesArray *mcArray)
{
  //
  // Select good tracks using fAnalCuts (AliRDHFCuts object) and return the array of their ids
  //
  
  //const Int_t entries = event->GetNumberOfTracks();
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
void AliAnalysisTaskSEXic2eleXifromAODtracks::SelectCascade( const AliVEvent *event,Int_t nCascs,Int_t &nSeleCasc, Bool_t *seleCascFlags, TClonesArray *mcArray)
{
  //
  // Select good Casc using fAnalCuts (AliRDHFCuts object) and return the array of their ids
  //

	Double_t primVtx[3];
	fVtx1->GetXYZ(primVtx);

  nSeleCasc = 0;
  for(Int_t icasc=0;icasc<nCascs;icasc++)
    {
      seleCascFlags[icasc] = kFALSE;
      AliAODcascade *casc = ((AliAODEvent*)event)->GetCascade(icasc);

      if(!fAnalCuts) continue;
      if(fAnalCuts->SingleCascadeCuts(casc,primVtx)){
				seleCascFlags[icasc] = kTRUE;
				nSeleCasc++;

//				AliAODTrack *cptrack =  (AliAODTrack*)(casc->GetDaughter(0));
//				AliAODTrack *cntrack =  (AliAODTrack*)(casc->GetDaughter(1));
//				AliAODTrack *cbtrack =  (AliAODTrack*)(casc->GetDecayVertexXi()->GetDaughter(0));
//				cout<<"Xi: "<<casc->MomXiX()<<" "<<casc->MomXiY()<<" "<<casc->MomXiZ()<<endl;
//				cout<<"V0: "<<casc->Px()<<" "<<casc->Py()<<" "<<casc->Pz()<<endl;
//				cout<<"from casc"<<endl;
//				cout<<casc->DecayVertexV0X()<<" "<<casc->DecayVertexV0Y()<<" "<<casc->DecayVertexV0Z()<<" "<<casc->MomPosX()<<" "<<casc->MomPosY()<<" "<<casc->MomPosZ()<<endl;
//				cout<<casc->DecayVertexV0X()<<" "<<casc->DecayVertexV0Y()<<" "<<casc->DecayVertexV0Z()<<" "<<casc->MomNegX()<<" "<<casc->MomNegY()<<" "<<casc->MomNegZ()<<endl;
//				cout<<casc->DecayVertexXiX()<<" "<<casc->DecayVertexXiY()<<" "<<casc->DecayVertexXiZ()<<" "<<casc->MomBachX()<<" "<<casc->MomBachY()<<" "<<casc->MomBachZ()<<endl;
//				cout<<"from track"<<endl;
//				cout<<cptrack->Xv()<<" "<<cptrack->Yv()<<" "<<cptrack->Zv()<<" "<<cptrack->Px()<<" "<<cptrack->Py()<<" "<<cptrack->Pz()<<endl;
//				cout<<cntrack->Xv()<<" "<<cntrack->Yv()<<" "<<cntrack->Zv()<<" "<<cntrack->Px()<<" "<<cntrack->Py()<<" "<<cntrack->Pz()<<endl;
//				cout<<cbtrack->Xv()<<" "<<cbtrack->Yv()<<" "<<cbtrack->Zv()<<" "<<cbtrack->Px()<<" "<<cbtrack->Py()<<" "<<cbtrack->Pz()<<endl;
//
//				cout<<"pos charge: "<<cptrack->Charge()<<endl;
//				cout<<"exact"<<endl;
//				Double_t xyz_cptrack[3];
//				cptrack->GetXYZAt(3.9,fBzkG,xyz_cptrack);
//				cout<<xyz_cptrack[0]<<" "<<xyz_cptrack[1]<<" "<<xyz_cptrack[2]<<endl;
//				cout<<"hand"<<endl;
//				cout<<"i am here1"<<endl;
//				Double_t v0vertr = 0.01*sqrt(casc->DecayVertexV0X()*casc->DecayVertexV0X()+casc->DecayVertexV0Y()*casc->DecayVertexV0Y());
//				cout<<"i am here2"<<endl;
//				Double_t v0pospt = sqrt(casc->MomPosX()*casc->MomPosX()+casc->MomPosY()*casc->MomPosY());
//				cout<<"i am here3"<<endl;
//				Double_t phi_atR = atan2(casc->DecayVertexV0Y(),casc->DecayVertexV0X())-asin((1.*0.3*0.5*v0vertr)/(2*v0pospt))+asin((1.*0.3*0.5*0.039)/(2*v0pospt));
//				cout<<"i am here4"<<endl;
//				//Double_t phi_atR = phi-asin((charge*0.3*0.5*R)/(2*pt));
//				if(phi_atR<-M_PI) phi_atR += 2 * M_PI;
//				if(phi_atR>M_PI) phi_atR -= 2 * M_PI;
//				cout<<3.9*cos(phi_atR)<<" "<<3.9*sin(phi_atR)<<endl;
//				cout<<endl;

				FillCascROOTObjects(casc,mcArray);
      }
    }
}
//_________________________________________________________________
Int_t AliAnalysisTaskSEXic2eleXifromAODtracks::GetPoolIndex(Double_t zvert, Double_t mult){
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
void AliAnalysisTaskSEXic2eleXifromAODtracks::ResetPool(Int_t poolIndex){
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
void AliAnalysisTaskSEXic2eleXifromAODtracks::DoEventMixingWithPools(Int_t poolIndex,AliAODEvent *aodEvent, Bool_t *seleFlags)
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
  for (Int_t i=0; i<aodEvent->GetNumberOfCascades(); i++)
  {
    if(!seleFlags[i]) continue;
    AliAODcascade* casc = aodEvent->GetCascade(i);
    if(!casc)continue;

		for(Int_t iEv=0; iEv<fNumberOfEventsForMixing; iEv++){
			fEventBuffer[poolIndex]->GetEvent(iEv + nEvents - fNumberOfEventsForMixing);
			TObjArray* earray1=(TObjArray*)earray->Clone();
			//Float_t zVertex1=zVertex;
			//Float_t mult1=cent;
			Int_t nElectrons=earray1->GetEntries();
			//Int_t evId1,esdId1,ne1;
			//sscanf((eventInfo->String()).Data(),"Ev%d_esd%d_K%d",&evId1,&esdId1,&ne1);
//			if(ne1!=nElectrons){ 
//				printf("AliAnalysisTaskSEXic2eleXifromAODtracks::DoMixingWithPools ERROR: read event does not match to the stored one\n");
//				delete earray1;
//				continue;
//			}
      for(Int_t iTr1=0; iTr1<nElectrons; iTr1++){
				AliAODTrack* trk1=(AliAODTrack*)earray1->At(iTr1);
				if(!trk1) continue;

				AliAODVertex *secVert = ReconstructSecondaryVertex(casc,trk1,aodEvent);//Fake, prim vertex is just used as secondary vertex. place holder for future
				if(!secVert) continue;

				AliAODRecoCascadeHF *exobj = MakeCascadeHF(casc,trk1,aodEvent,secVert,true);
				if(!exobj) {
						continue;
				}

				TClonesArray *fake = 0;
				FillROOTObjects(exobj,casc,trk1,fake,true);

				exobj->GetSecondaryVtx()->RemoveDaughters();
				exobj->UnsetOwnPrimaryVtx();
				delete exobj;exobj=NULL;
				delete secVert;
			}//track loop

			delete earray1;
		}//event loop
		
	}//v0 loop
}
//_________________________________________________________________
void AliAnalysisTaskSEXic2eleXifromAODtracks::MakeMCAnalysis(TClonesArray *mcArray)
{
	//
  // Analyze AliAODmcparticle
	//

	Int_t nmcpart = mcArray->GetEntriesFast();
	for(Int_t i=0;i<nmcpart;i++)
	{
		AliAODMCParticle *mcpart = (AliAODMCParticle*) mcArray->At(i);
		if(TMath::Abs(mcpart->GetPdgCode())==4132){
			Bool_t e_flag = kFALSE;
			Bool_t xi_flag = kFALSE;
			AliAODMCParticle *mcepart = 0;
			AliAODMCParticle *mccascpart = 0;
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
				}
			}

			Int_t decaytype = -9999;
			if(e_flag && xi_flag) decaytype = 0;

			FillMCROOTObjects(mcpart,mcepart,mccascpart,decaytype);
		}
		if(TMath::Abs(mcpart->GetPdgCode())==11 && mcpart->GetStatus()==1){
			AliESDtrackCuts *esdcuts = fAnalCuts->GetTrackCuts();
			Float_t etamin, etamax;
			esdcuts->GetEtaRange(etamin,etamax);
			if(fabs(mcpart->Eta())<etamax)
				fHistoBachPtMCGen->Fill(mcpart->Pt());
		}
		if(TMath::Abs(mcpart->GetPdgCode())==3312){
			Double_t etamin, etamax, rapmin, rapmax;
			fAnalCuts->GetProdCascEtaRange(etamin,etamax);
			fAnalCuts->GetProdCascRapRange(rapmin,rapmax);

			if((fabs(mcpart->Y())<rapmax) && (fabs(mcpart->Eta())<etamax))
				fHistoXiMassvsPtMCGen->Fill(1.32171, mcpart->Pt());
		}
	}
	return;
}
