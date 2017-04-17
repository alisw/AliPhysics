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
//               Omegac->eOmega analysis code
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
#include "AliAnalysisTaskSEOmegac2eleOmegafromAODtracks.h"
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
ClassImp(AliAnalysisTaskSEOmegac2eleOmegafromAODtracks);
/// \endcond

//__________________________________________________________________________
AliAnalysisTaskSEOmegac2eleOmegafromAODtracks::AliAnalysisTaskSEOmegac2eleOmegafromAODtracks() : 
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
  fTriggerCheck(0),
  fUseCentralityV0M(kFALSE),
  fEvNumberCounter(0),
  fHistoEleOmegaMass(0),
  fHistoEleOmegaMassRS(0),
  fHistoEleOmegaMassWS(0),
  fHistoEleOmegaMassRSMix(0),
  fHistoEleOmegaMassWSMix(0),
  fHistoEleOmegaMassvsElePtRS(0),
  fHistoEleOmegaMassvsElePtWS(0),
  fHistoEleOmegaMassvsElePtRSMix(0),
  fHistoEleOmegaMassvsElePtWSMix(0),
  fHistoElePtRS(0),
  fHistoElePtWS(0),
  fHistoElePtRSMix(0),
  fHistoElePtWSMix(0),
  fHistoEleOmegaMassMCS(0),
  fHistoEleOmegaMassMCGen(0),
  fHistoEleOmegaMassvsElePtMCS(0),
  fHistoEleOmegaMassvsElePtMCGen(0),
  fHistoElePtMCS(0),
  fHistoElePtMCGen(0),
  fHistoElePtvsEtaRS(0),
  fHistoElePtvsEtaWS(0),
  fHistoElePtvsEtaRSMix(0),
  fHistoElePtvsEtaWSMix(0),
  fHistoElePtvsEtaMCS(0),
  fHistoElePtvsEtaMCGen(0),
  fHistoElePtvsOmegaPtRS(0),
  fHistoElePtvsOmegaPtWS(0),
  fHistoElePtvsOmegaPtRSMix(0),
  fHistoElePtvsOmegaPtWSMix(0),
  fHistoElePtvsOmegaPtMCS(0),
  fHistoElePtvsOmegaPtMCGen(0),
  fHistoElePtvsd0RS(0),
  fHistoElePtvsd0WS(0),
  fHistoElePtvsd0RSMix(0),
  fHistoElePtvsd0WSMix(0),
  fHistoElePtvsd0MCS(0),
  fHistoBachPt(0),
  fHistoBachPtMCS(0),
  fHistoBachPtMCGen(0),
  fHistod0Bach(0),
  fHistoXiMassvsPt(0),
  fHistoOmegaMassvsPt(0),
  fHistoOmegaMassvsPtMCS(0),
  fHistoOmegaMassvsPtMCGen(0),
  fHistoElectronTPCPID(0),
  fHistoElectronTOFPID(0),
  fHistoElectronTPCSelPID(0),
  fHistoElectronTOFSelPID(0),
  fHistoElectronTPCPIDSelTOF(0),
  fHistoElectronTPCPIDSelTOFSmallEta(0),
  fHistoElectronTPCPIDSelTOFLargeEta(0),
	fCounter(0),
	fHistonEvtvsRunNumber(0),
	fHistonElevsRunNumber(0),
	fHistonOmegavsRunNumber(0),
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
}

//___________________________________________________________________________
AliAnalysisTaskSEOmegac2eleOmegafromAODtracks::AliAnalysisTaskSEOmegac2eleOmegafromAODtracks(const Char_t* name,
									     AliRDHFCutsOmegactoeleOmegafromAODtracks* analCuts, 
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
  fTriggerCheck(0),
  fUseCentralityV0M(kFALSE),
  fEvNumberCounter(0),
  fHistoEleOmegaMass(0),
  fHistoEleOmegaMassRS(0),
  fHistoEleOmegaMassWS(0),
  fHistoEleOmegaMassRSMix(0),
  fHistoEleOmegaMassWSMix(0),
  fHistoEleOmegaMassvsElePtRS(0),
  fHistoEleOmegaMassvsElePtWS(0),
  fHistoEleOmegaMassvsElePtRSMix(0),
  fHistoEleOmegaMassvsElePtWSMix(0),
  fHistoElePtRS(0),
  fHistoElePtWS(0),
  fHistoElePtRSMix(0),
  fHistoElePtWSMix(0),
  fHistoEleOmegaMassMCS(0),
  fHistoEleOmegaMassMCGen(0),
  fHistoEleOmegaMassvsElePtMCS(0),
  fHistoEleOmegaMassvsElePtMCGen(0),
  fHistoElePtMCS(0),
  fHistoElePtMCGen(0),
  fHistoElePtvsEtaRS(0),
  fHistoElePtvsEtaWS(0),
  fHistoElePtvsEtaRSMix(0),
  fHistoElePtvsEtaWSMix(0),
  fHistoElePtvsEtaMCS(0),
  fHistoElePtvsEtaMCGen(0),
  fHistoElePtvsOmegaPtRS(0),
  fHistoElePtvsOmegaPtWS(0),
  fHistoElePtvsOmegaPtRSMix(0),
  fHistoElePtvsOmegaPtWSMix(0),
  fHistoElePtvsOmegaPtMCS(0),
  fHistoElePtvsOmegaPtMCGen(0),
  fHistoElePtvsd0RS(0),
  fHistoElePtvsd0WS(0),
  fHistoElePtvsd0RSMix(0),
  fHistoElePtvsd0WSMix(0),
  fHistoElePtvsd0MCS(0),
  fHistoBachPt(0),
  fHistoBachPtMCS(0),
  fHistoBachPtMCGen(0),
  fHistod0Bach(0),
  fHistoXiMassvsPt(0),
  fHistoOmegaMassvsPt(0),
  fHistoOmegaMassvsPtMCS(0),
  fHistoOmegaMassvsPtMCGen(0),
  fHistoElectronTPCPID(0),
  fHistoElectronTOFPID(0),
  fHistoElectronTPCSelPID(0),
  fHistoElectronTOFSelPID(0),
  fHistoElectronTPCPIDSelTOF(0),
  fHistoElectronTPCPIDSelTOFSmallEta(0),
  fHistoElectronTPCPIDSelTOFLargeEta(0),
	fCounter(0),
	fHistonEvtvsRunNumber(0),
	fHistonElevsRunNumber(0),
	fHistonOmegavsRunNumber(0),
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
  Info("AliAnalysisTaskSEOmegac2eleOmegafromAODtracks","Calling Constructor");

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
AliAnalysisTaskSEOmegac2eleOmegafromAODtracks::~AliAnalysisTaskSEOmegac2eleOmegafromAODtracks() {
  //
  // destructor
  //
  Info("~AliAnalysisTaskSEOmegac2eleOmegafromAODtracks","Calling Destructor");

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


}

//_________________________________________________
void AliAnalysisTaskSEOmegac2eleOmegafromAODtracks::Init() {
  //
  // Initialization
  //
  //

  fIsEventSelected=kFALSE;

  if (fDebug > 1) AliInfo("Init");

  fListCuts = new TList();
  fListCuts->SetOwner();
  fListCuts->SetName("ListCuts");
  fListCuts->Add(new AliRDHFCutsOmegactoeleOmegafromAODtracks(*fAnalCuts));
  PostData(2,fListCuts);

  return;
}

//_________________________________________________
void AliAnalysisTaskSEOmegac2eleOmegafromAODtracks::UserExec(Option_t *)
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
  // First check if the event has proper vertex and B
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
      AliError("AliAnalysisTaskSEOmegac2eleOmegafromAODtracks::UserExec: MC header branch not found!\n");
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
	fEvNumberCounter++;

	Int_t runnumber_offset = 0;
	Int_t runnumber = aodEvent->GetRunNumber();
	if(runnumber<=117222&&runnumber>=114931){
		runnumber_offset = 114931;//lhc10b
	}else if(runnumber<=120829&&runnumber>=119159){
		runnumber_offset = 119159;//lhc10c
	}else if(runnumber<=126437&&runnumber>=122374){
		runnumber_offset = 122374;//lhc10d
	}else if(runnumber<=130840&&runnumber>=127712){
		runnumber_offset = 127712;//lhc10e
	}else if(runnumber<=195483&&runnumber>=195344){
		runnumber_offset = 195344;//lhc13b
	}else if(runnumber<=195677&&runnumber>=195529){
		runnumber_offset = 195529;//lhc13c
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
void AliAnalysisTaskSEOmegac2eleOmegafromAODtracks::Terminate(Option_t*)
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
void AliAnalysisTaskSEOmegac2eleOmegafromAODtracks::UserCreateOutputObjects() 
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
void AliAnalysisTaskSEOmegac2eleOmegafromAODtracks::MakeAnalysis
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
	if(runnumber<=117222&&runnumber>=114931){
		runnumber_offset = 114931;//lhc10b
	}else if(runnumber<=120829&&runnumber>=119159){
		runnumber_offset = 119159;//lhc10c
	}else if(runnumber<=126437&&runnumber>=122374){
		runnumber_offset = 122374;//lhc10d
	}else if(runnumber<=130840&&runnumber>=127712){
		runnumber_offset = 127712;//lhc10e
	}else if(runnumber<=195483&&runnumber>=195344){
		runnumber_offset = 195344;//lhc13b
	}else if(runnumber<=195677&&runnumber>=195529){
		runnumber_offset = 195529;//lhc13c
	}else if(runnumber<=170593&&runnumber>=167902){
		runnumber_offset = 167902;//lhc11h
	}
	fHistonElevsRunNumber->Fill(runnumber-runnumber_offset,nSeleTrks);
	fHistonOmegavsRunNumber->Fill(runnumber-runnumber_offset,nSeleCasc);

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
void AliAnalysisTaskSEOmegac2eleOmegafromAODtracks::DefineTreeVariables() 
{
  //
  // Define tree variables
  //

  const char* nameoutput = GetOutputSlot(4)->GetContainer()->GetName();
  fVariablesTree = new TTree(nameoutput,"Candidates variables tree");
  Int_t nVar = 64;
  fCandidateVariables = new Float_t [nVar];
  TString * fCandidateVariableNames = new TString[nVar];

  fCandidateVariableNames[ 0]="Centrality";
  fCandidateVariableNames[ 1]="InvMassEleOmega";
  fCandidateVariableNames[ 2]="EleOmegaPt";
  fCandidateVariableNames[ 3]="EleOmegaPx";
  fCandidateVariableNames[ 4]="EleOmegaPy";
  fCandidateVariableNames[ 5]="EleOmegaPz";
  fCandidateVariableNames[ 6]="ElePx";
  fCandidateVariableNames[ 7]="ElePy";
  fCandidateVariableNames[ 8]="ElePz";
  fCandidateVariableNames[ 9]="OmegaPx";
  fCandidateVariableNames[10]="OmegaPy";
  fCandidateVariableNames[11]="OmegaPz";
  fCandidateVariableNames[12]="OmegaCharge";
  fCandidateVariableNames[13]="MassOmega";
  fCandidateVariableNames[14]="MassLambda";
  fCandidateVariableNames[15]="Eled0";
  fCandidateVariableNames[16]="Omegad0";
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
  fCandidateVariableNames[27]="nSigmaTPCbachka";
  fCandidateVariableNames[28]="nSigmaTOFbachka";
  fCandidateVariableNames[29]="EleCharge";
  fCandidateVariableNames[30]="Mixing";
  fCandidateVariableNames[31]="DcaOmegaDaughters";
  fCandidateVariableNames[32]="DcaV0Daughters";
  fCandidateVariableNames[33]="DecayLengthXi";
  fCandidateVariableNames[34]="CosPointingAngleXi";
  fCandidateVariableNames[35]="DcaV0toPrimVertex";
  fCandidateVariableNames[36]="DcaPostoPrimVertex";
  fCandidateVariableNames[37]="DcaNegtoPrimVertex";
  fCandidateVariableNames[38]="DcaBachtoPrimVertex";
  fCandidateVariableNames[39]="DecayLengthV0";
  fCandidateVariableNames[40]="CosPointingAngleV0";

  fCandidateVariableNames[41]="mcpdgomegac";
  fCandidateVariableNames[42]="mclabomegac";
  fCandidateVariableNames[43]="mcomegacpx";
  fCandidateVariableNames[44]="mcomegacpy";
  fCandidateVariableNames[45]="mcomegacpz";
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
  fCandidateVariableNames[62]="MassXi";

  fCandidateVariableNames[63]="EvNumber";

  for (Int_t ivar=0; ivar<nVar; ivar++) {
    fVariablesTree->Branch(fCandidateVariableNames[ivar].Data(),&fCandidateVariables[ivar],Form("%s/f",fCandidateVariableNames[ivar].Data()));
  }

  return;
}

////-------------------------------------------------------------------------------
void AliAnalysisTaskSEOmegac2eleOmegafromAODtracks::FillROOTObjects(AliAODRecoCascadeHF *exobj, AliAODcascade *casc, AliAODTrack *trk, TClonesArray *mcArray, Bool_t mixing_flag) 
{
  //
  // Fill histograms or tree depending on fWriteVariableTree 
  //
	if(!trk) return;
	if(!casc) return;

	for(Int_t i=0;i<64;i++){
		fCandidateVariables[i] = -9999.;
	}


  AliAODTrack *cptrack =  (AliAODTrack*)(casc->GetDaughter(0));
  AliAODTrack *cntrack =  (AliAODTrack*)(casc->GetDaughter(1));
  AliAODTrack *cbtrack =  (AliAODTrack*)(casc->GetDecayVertexXi()->GetDaughter(0));
  Double_t posVtx[3] = {0.,0.,0.};
  fVtx1->GetXYZ(posVtx);


  fCandidateVariables[ 0] = fCentrality;
	UInt_t pdgdg[2]={11,3334};
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
  fCandidateVariables[13] = casc->MassOmega();
	if(casc->ChargeXi()<0)
		fCandidateVariables[14] = casc->MassLambda();
	else
		fCandidateVariables[14] = casc->MassAntiLambda();
  fCandidateVariables[15] = exobj->Getd0Prong(0);
  fCandidateVariables[16] = exobj->Getd0Prong(1);

  if(fAnalCuts->GetIsUsePID()&&!mixing_flag)
  {
		Double_t nSigmaTPCele = fAnalCuts->GetPidHF()->GetPidResponse()->NumberOfSigmasTPC(trk,AliPID::kElectron);
		Double_t nSigmaTOFele = fAnalCuts->GetPidHF()->GetPidResponse()->NumberOfSigmasTOF(trk,AliPID::kElectron);
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

	if(fAnalCuts->GetUseCascadePID()&&!mixing_flag)
	{
		Double_t nSigmaTPCv0pr=-9999.;
		Double_t nSigmaTOFv0pr=-9999.;
		Double_t nSigmaTPCv0pi=-9999.;
		Double_t nSigmaTOFv0pi=-9999.;
		Double_t nSigmaTPCbachka=-9999.;
		Double_t nSigmaTOFbachka=-9999.;
		if(casc->ChargeXi()>0){
			nSigmaTPCv0pr = fAnalCuts->GetPidCascPr()->GetPidResponse()->NumberOfSigmasTPC(cntrack,AliPID::kProton);
			nSigmaTOFv0pr = fAnalCuts->GetPidCascPr()->GetPidResponse()->NumberOfSigmasTOF(cntrack,AliPID::kProton);
			nSigmaTPCv0pi = fAnalCuts->GetPidCascPi()->GetPidResponse()->NumberOfSigmasTPC(cptrack,AliPID::kPion);
			nSigmaTOFv0pi = fAnalCuts->GetPidCascPi()->GetPidResponse()->NumberOfSigmasTOF(cptrack,AliPID::kPion);
			nSigmaTPCbachka = fAnalCuts->GetPidCascPi()->GetPidResponse()->NumberOfSigmasTPC(cbtrack,AliPID::kKaon);
			nSigmaTOFbachka = fAnalCuts->GetPidCascPi()->GetPidResponse()->NumberOfSigmasTOF(cbtrack,AliPID::kKaon);
		}else{
			nSigmaTPCv0pr = fAnalCuts->GetPidCascPr()->GetPidResponse()->NumberOfSigmasTPC(cptrack,AliPID::kProton);
			nSigmaTOFv0pr = fAnalCuts->GetPidCascPr()->GetPidResponse()->NumberOfSigmasTOF(cptrack,AliPID::kProton);
			nSigmaTPCv0pi = fAnalCuts->GetPidCascPi()->GetPidResponse()->NumberOfSigmasTPC(cntrack,AliPID::kPion);
			nSigmaTOFv0pi = fAnalCuts->GetPidCascPi()->GetPidResponse()->NumberOfSigmasTOF(cntrack,AliPID::kPion);
			nSigmaTPCbachka = fAnalCuts->GetPidCascPi()->GetPidResponse()->NumberOfSigmasTPC(cbtrack,AliPID::kKaon);
			nSigmaTOFbachka = fAnalCuts->GetPidCascPi()->GetPidResponse()->NumberOfSigmasTOF(cbtrack,AliPID::kKaon);
		}
      fCandidateVariables[23] = nSigmaTPCv0pr;
      fCandidateVariables[24] = nSigmaTOFv0pr;
      fCandidateVariables[25] = nSigmaTPCv0pi;
      fCandidateVariables[26] = nSigmaTOFv0pi;
      fCandidateVariables[27] = nSigmaTPCbachka;
      fCandidateVariables[28] = nSigmaTOFbachka;
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

  AliAODMCParticle *mcomegac = 0;
  AliAODMCParticle *mcele = 0;
  AliAODMCParticle *mccasc = 0;
  Int_t mclabomegac = 0;
	Int_t mcpdgele_array[100];
	Int_t mcpdgcasc_array[100];
	Int_t mclabelele_array[100];
	Int_t mclabelcasc_array[100];
	Int_t mcngen_ele = -9999;
	Int_t mcngen_casc = -9999;

	if(fUseMCInfo && mcArray){

		mclabomegac =  MatchToMC(exobj,mcArray,mcpdgele_array, mcpdgcasc_array,mclabelele_array,mclabelcasc_array,mcngen_ele,mcngen_casc);

    if(mclabomegac>-1){
      mcomegac = (AliAODMCParticle*) mcArray->At(mclabomegac);
			if(mclabelele_array[0]>=0)
				mcele = (AliAODMCParticle*) mcArray->At(mclabelele_array[0]);
			if(mclabelcasc_array[0]>=0)
				mccasc = (AliAODMCParticle*) mcArray->At(mclabelcasc_array[0]);
			if(mcomegac){
				fCandidateVariables[41] = mcomegac->GetPdgCode();
				fCandidateVariables[42] = mcomegac->Label();
				fCandidateVariables[43] = mcomegac->Px();
				fCandidateVariables[44] = mcomegac->Py();
				fCandidateVariables[45] = mcomegac->Pz();
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

  fCandidateVariables[62] = casc->MassXi();
  fCandidateVariables[63] = fEvNumberCounter;


  if(fWriteVariableTree)
    fVariablesTree->Fill();

	if(fAnalCuts->IsSelected(exobj,AliRDHFCuts::kCandidate))
	{
		Double_t cont[3];
		cont[0] = exobj->InvMass(2,pdgdg);
		cont[1] = exobj->Pt();
		cont[2] = fCentrality;
		fHistoEleOmegaMass->Fill(cont);

		Double_t cont2[3];
		cont2[0] = exobj->InvMass(2,pdgdg);
		cont2[1] = trk->Pt();
		cont2[2] = fCentrality;

		Double_t cont_eleptvseta[3];
		cont_eleptvseta[0] = trk->Pt();
		cont_eleptvseta[1] = trk->Eta();
		cont_eleptvseta[2] = fCentrality;

		Double_t cont_eleptvsomegapt[3];
		cont_eleptvsomegapt[0] = trk->Pt();
		cont_eleptvsomegapt[1] = sqrt(casc->MomXiX()*casc->MomXiX()+casc->MomXiY()*casc->MomXiY());
		cont_eleptvsomegapt[2] = fCentrality;

		Double_t cont_eleptvsd0[3];
		cont_eleptvsd0[0] = trk->Pt();
		cont_eleptvsd0[1] = exobj->Getd0Prong(0);
		cont_eleptvsd0[2] = fCentrality;

		if(mixing_flag){
			if(trk->Charge()*casc->ChargeXi()<0){
				fHistoEleOmegaMassRSMix->Fill(cont);
				fHistoEleOmegaMassvsElePtRSMix->Fill(cont2);
				if(cont[0]<2.7){
					fHistoElePtRSMix->Fill(trk->Pt(),fCentrality);
					fHistoElePtvsEtaRSMix->Fill(cont_eleptvseta);
					fHistoElePtvsOmegaPtRSMix->Fill(cont_eleptvsomegapt);
					fHistoElePtvsd0RSMix->Fill(cont_eleptvsd0);
				}
			}else{
				fHistoEleOmegaMassWSMix->Fill(cont);
				fHistoEleOmegaMassvsElePtWSMix->Fill(cont2);
				if(cont[0]<2.7){
					fHistoElePtWSMix->Fill(trk->Pt(),fCentrality);
					fHistoElePtvsEtaWSMix->Fill(cont_eleptvseta);
					fHistoElePtvsOmegaPtWSMix->Fill(cont_eleptvsomegapt);
					fHistoElePtvsd0WSMix->Fill(cont_eleptvsd0);
				}
			}
		}else{
			if(trk->Charge()*casc->ChargeXi()<0){
				fHistoEleOmegaMassRS->Fill(cont);
				fHistoEleOmegaMassvsElePtRS->Fill(cont2);
				if(cont[0]<2.7){
					fHistoElePtRS->Fill(trk->Pt(),fCentrality);
					fHistoElePtvsEtaRS->Fill(cont_eleptvseta);
					fHistoElePtvsOmegaPtRS->Fill(cont_eleptvsomegapt);
					fHistoElePtvsd0RS->Fill(cont_eleptvsd0);
				}
			}else{
				fHistoEleOmegaMassWS->Fill(cont);
				fHistoEleOmegaMassvsElePtWS->Fill(cont2);
				if(cont[0]<2.7){
					fHistoElePtWS->Fill(trk->Pt(),fCentrality);
					fHistoElePtvsEtaWS->Fill(cont_eleptvseta);
					fHistoElePtvsOmegaPtWS->Fill(cont_eleptvsomegapt);
					fHistoElePtvsd0WS->Fill(cont_eleptvsd0);
				}
			}
		}

		if(fUseMCInfo){
			if(mcomegac){
				Int_t pdgcode = mcomegac->GetPdgCode();
				if(abs(pdgcode)==4332 && abs(mcpdgele_array[1])==4332 && abs(mcpdgcasc_array[1])==4332){
						fHistoEleOmegaMassMCS->Fill(cont);
						fHistoEleOmegaMassvsElePtMCS->Fill(cont2);
						if(cont[0]<2.7){
							fHistoElePtMCS->Fill(trk->Pt(),fCentrality);
							fHistoElePtvsEtaMCS->Fill(cont_eleptvseta);
							fHistoElePtvsOmegaPtMCS->Fill(cont_eleptvsomegapt);
							fHistoElePtvsd0MCS->Fill(cont_eleptvsd0);
						}
				}
			}
		}
	}

  return;
}

////-------------------------------------------------------------------------------
void AliAnalysisTaskSEOmegac2eleOmegafromAODtracks::DefineEleTreeVariables() 
{
  //
  // Define electron tree variables
  //

  const char* nameoutput = GetOutputSlot(5)->GetContainer()->GetName();
  fEleVariablesTree = new TTree(nameoutput,"electron variables tree");
  Int_t nVar = 19;
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

  for (Int_t ivar=0; ivar<nVar; ivar++) {
    fEleVariablesTree->Branch(fCandidateVariableNames[ivar].Data(),&fCandidateEleVariables[ivar],Form("%s/f",fCandidateVariableNames[ivar].Data()));
  }

  return;
}
////-------------------------------------------------------------------------------
void AliAnalysisTaskSEOmegac2eleOmegafromAODtracks::FillElectronROOTObjects(AliAODTrack *trk, TClonesArray *mcArray) 
{
  //
  // Fill histograms or tree depending on fWriteVariableTree 
  //

	if(!trk) return;

	fHistoBachPt->Fill(trk->Pt());

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

	for(Int_t i=0;i<19;i++){
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

	fHistod0Bach->Fill(d0z0[0]);

	fEleVariablesTree->Fill();
}
////-------------------------------------------------------------------------------
void AliAnalysisTaskSEOmegac2eleOmegafromAODtracks::DefineCascTreeVariables() 
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
  fCandidateVariableNames[ 1]="InvMassOmega";
  fCandidateVariableNames[ 2]="OmegaPx";
  fCandidateVariableNames[ 3]="OmegaPy";
  fCandidateVariableNames[ 4]="OmegaPz";
  fCandidateVariableNames[ 5]="InvMassLambda";
  fCandidateVariableNames[ 6]="DcaOmegaDaughters";
  fCandidateVariableNames[ 7]="DcaV0Daughters";
  fCandidateVariableNames[ 8]="DecayLengthOmega";
  fCandidateVariableNames[ 9]="CosPointingAngleOmega";
  fCandidateVariableNames[10]="DcaV0toPrimVertex";
  fCandidateVariableNames[11]="DcaPostoPrimVertex";
  fCandidateVariableNames[12]="DcaNegtoPrimVertex";
  fCandidateVariableNames[13]="DcaBachtoPrimVertex";
  fCandidateVariableNames[14]="DecayLengthV0";
  fCandidateVariableNames[15]="CosPointingAngleV0";
  fCandidateVariableNames[16]="OmegaCharge";

  for (Int_t ivar=0; ivar<nVar; ivar++) {
    fCascVariablesTree->Branch(fCandidateVariableNames[ivar].Data(),&fCandidateCascVariables[ivar],Form("%s/f",fCandidateVariableNames[ivar].Data()));
  }

  return;
}

////-------------------------------------------------------------------------------
void AliAnalysisTaskSEOmegac2eleOmegafromAODtracks::FillCascROOTObjects(AliAODcascade *casc, TClonesArray *mcArray) 
{
  //
  // Fill histograms or tree depending on fWriteVariableTree (tree not implemented yet)
  //
	if(!casc) return;
	fHistoXiMassvsPt->Fill(casc->MassXi(),sqrt(casc->MomXiX()*casc->MomXiX()+casc->MomXiY()*casc->MomXiY()));
	fHistoOmegaMassvsPt->Fill(casc->MassOmega(),sqrt(casc->MomXiX()*casc->MomXiX()+casc->MomXiY()*casc->MomXiY()));

	if(fUseMCInfo){
		Int_t pdgDgcasc[2]={321,3122};
		Int_t pdgDgv0[2]={2212,211};
		Int_t labcasc = MatchToMCCascade(casc,3334,pdgDgcasc,pdgDgv0,mcArray); // the cascade
		if(labcasc>=0){
			fHistoOmegaMassvsPtMCS->Fill(casc->MassOmega(),sqrt(casc->MomXiX()*casc->MomXiX()+casc->MomXiY()*casc->MomXiY()));
		}
	}

	if(!fWriteEachVariableTree) return;


	for(Int_t i=0;i<16;i++){
		fCandidateCascVariables[i] = -9999.;
	}
  Double_t posVtx[3] = {0.,0.,0.};
  fVtx1->GetXYZ(posVtx);

  fCandidateCascVariables[ 0] = fCentrality;
  fCandidateCascVariables[ 1] = casc->MassOmega();
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
void AliAnalysisTaskSEOmegac2eleOmegafromAODtracks::DefineMCTreeVariables() 
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
void AliAnalysisTaskSEOmegac2eleOmegafromAODtracks::FillMCROOTObjects(AliAODMCParticle *mcpart, AliAODMCParticle *mcepart, AliAODMCParticle *mccascpart, Int_t decaytype) 
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
	Double_t cascE = sqrt(cascpx*cascpx+cascpy*cascpy+cascpz*cascpz+1.67245*1.67245);

	Double_t InvMassEleOmega = sqrt(pow(eE+cascE,2)-pow(epx+cascpx,2)-pow(epy+cascpy,2)-pow(epz+cascpz,2));

	Double_t cont[3];
	cont[0] = InvMassEleOmega;
	cont[1] = mcpart->Pt();
	cont[2] = fCentrality;
	Double_t cont2[3];
	cont2[0] = InvMassEleOmega;
	cont2[1] = mcepart->Pt();
	cont2[2] = fCentrality;
	Double_t cont_eleptvseta[3];
	cont_eleptvseta[0] = mcepart->Pt();
	cont_eleptvseta[1] = mcepart->Eta();
	cont_eleptvseta[2] = fCentrality;
	Double_t cont_eleptvsomegapt[3];
	cont_eleptvsomegapt[0] = mcepart->Pt();
	cont_eleptvsomegapt[1] = mccascpart->Pt();
	cont_eleptvsomegapt[2] = fCentrality;

	AliESDtrackCuts *esdcuts = fAnalCuts->GetTrackCuts();
	Float_t etamin, etamax;
	esdcuts->GetEtaRange(etamin,etamax);

	if(decaytype==0){
		fHistoEleOmegaMassMCGen->Fill(cont);
		if(fabs(mcepart->Eta())<etamax){
			fHistoEleOmegaMassvsElePtMCGen->Fill(cont2);
			if(InvMassEleOmega<2.7){
				fHistoElePtMCGen->Fill(mcepart->Pt(),fCentrality);
				fHistoElePtvsEtaMCGen->Fill(cont_eleptvseta);
				fHistoElePtvsOmegaPtMCGen->Fill(cont_eleptvsomegapt);
			}
		}
	}

	if(fWriteMCVariableTree)
		fMCVariablesTree->Fill();
}

////__________________________________________________________________________
void  AliAnalysisTaskSEOmegac2eleOmegafromAODtracks::DefineGeneralHistograms() {
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
void  AliAnalysisTaskSEOmegac2eleOmegafromAODtracks::DefineAnalysisHistograms() 
{
  //
  // Define analyis histograms
  //
	
  //------------------------------------------------
  // Basic histogram
  //------------------------------------------------
  Int_t bins_base[3]=		{16				,100		,10};
  Double_t xmin_base[3]={1.6,0		,0.00};
  Double_t xmax_base[3]={3.6,10.	,100};
  fHistoEleOmegaMass = new THnSparseF("fHistoEleOmegaMass","",3,bins_base,xmin_base,xmax_base);
  fOutputAll->Add(fHistoEleOmegaMass);
  fHistoEleOmegaMassRS = new THnSparseF("fHistoEleOmegaMassRS","",3,bins_base,xmin_base,xmax_base);
  fOutputAll->Add(fHistoEleOmegaMassRS);
  fHistoEleOmegaMassWS = new THnSparseF("fHistoEleOmegaMassWS","",3,bins_base,xmin_base,xmax_base);
  fOutputAll->Add(fHistoEleOmegaMassWS);
  fHistoEleOmegaMassRSMix = new THnSparseF("fHistoEleOmegaMassRSMix","",3,bins_base,xmin_base,xmax_base);
  fOutputAll->Add(fHistoEleOmegaMassRSMix);
  fHistoEleOmegaMassWSMix = new THnSparseF("fHistoEleOmegaMassWSMix","",3,bins_base,xmin_base,xmax_base);
  fOutputAll->Add(fHistoEleOmegaMassWSMix);

  fHistoEleOmegaMassvsElePtRS = new THnSparseF("fHistoEleOmegaMassvsElePtRS","",3,bins_base,xmin_base,xmax_base);
  fOutputAll->Add(fHistoEleOmegaMassvsElePtRS);
  fHistoEleOmegaMassvsElePtWS = new THnSparseF("fHistoEleOmegaMassvsElePtWS","",3,bins_base,xmin_base,xmax_base);
  fOutputAll->Add(fHistoEleOmegaMassvsElePtWS);
  fHistoEleOmegaMassvsElePtRSMix = new THnSparseF("fHistoEleOmegaMassvsElePtRSMix","",3,bins_base,xmin_base,xmax_base);
  fOutputAll->Add(fHistoEleOmegaMassvsElePtRSMix);
  fHistoEleOmegaMassvsElePtWSMix = new THnSparseF("fHistoEleOmegaMassvsElePtWSMix","",3,bins_base,xmin_base,xmax_base);
  fOutputAll->Add(fHistoEleOmegaMassvsElePtWSMix);

  fHistoElePtRS=new TH2F("fHistoElePtRS","Right-sign e pt",100,0.,10.,10,0.,100.);
  fOutputAll->Add(fHistoElePtRS);
  fHistoElePtWS=new TH2F("fHistoElePtWS","Wrong-sign e pt",100,0.,10.,10,0.,100.);
  fOutputAll->Add(fHistoElePtWS);
  fHistoElePtRSMix=new TH2F("fHistoElePtRSMix","Right-sign e pt",100,0.,10.,10,0.,100.);
  fOutputAll->Add(fHistoElePtRSMix);
  fHistoElePtWSMix=new TH2F("fHistoElePtWSMix","Wrong-sign e pt",100,0.,10.,10,0.,100.);
  fOutputAll->Add(fHistoElePtWSMix);

  fHistoEleOmegaMassMCS = new THnSparseF("fHistoEleOmegaMassMCS","",3,bins_base,xmin_base,xmax_base);
  fOutputAll->Add(fHistoEleOmegaMassMCS);
  fHistoEleOmegaMassMCGen = new THnSparseF("fHistoEleOmegaMassMCGen","",3,bins_base,xmin_base,xmax_base);
  fOutputAll->Add(fHistoEleOmegaMassMCGen);
  fHistoEleOmegaMassvsElePtMCS = new THnSparseF("fHistoEleOmegaMassvsElePtMCS","",3,bins_base,xmin_base,xmax_base);
  fOutputAll->Add(fHistoEleOmegaMassvsElePtMCS);
  fHistoEleOmegaMassvsElePtMCGen = new THnSparseF("fHistoEleOmegaMassvsElePtMCGen","",3,bins_base,xmin_base,xmax_base);
  fOutputAll->Add(fHistoEleOmegaMassvsElePtMCGen);
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

  Int_t bins_eleptvsomegapt[3]=	{50,20	,10};
  Double_t xmin_eleptvsomegapt[3]={0.,0.	,0.0};
  Double_t xmax_eleptvsomegapt[3]={5.,5.	,100};

  fHistoElePtvsOmegaPtRS = new THnSparseF("fHistoElePtvsOmegaPtRS","",3,bins_eleptvsomegapt,xmin_eleptvsomegapt,xmax_eleptvsomegapt);
  fOutputAll->Add(fHistoElePtvsOmegaPtRS);
  fHistoElePtvsOmegaPtWS = new THnSparseF("fHistoElePtvsOmegaPtWS","",3,bins_eleptvsomegapt,xmin_eleptvsomegapt,xmax_eleptvsomegapt);
  fOutputAll->Add(fHistoElePtvsOmegaPtWS);
  fHistoElePtvsOmegaPtRSMix = new THnSparseF("fHistoElePtvsOmegaPtRSMix","",3,bins_eleptvsomegapt,xmin_eleptvsomegapt,xmax_eleptvsomegapt);
  fOutputAll->Add(fHistoElePtvsOmegaPtRSMix);
  fHistoElePtvsOmegaPtWSMix = new THnSparseF("fHistoElePtvsOmegaPtWSMix","",3,bins_eleptvsomegapt,xmin_eleptvsomegapt,xmax_eleptvsomegapt);
  fOutputAll->Add(fHistoElePtvsOmegaPtWSMix);
  fHistoElePtvsOmegaPtMCS = new THnSparseF("fHistoElePtvsOmegaPtMCS","",3,bins_eleptvsomegapt,xmin_eleptvsomegapt,xmax_eleptvsomegapt);
  fOutputAll->Add(fHistoElePtvsOmegaPtMCS);
  fHistoElePtvsOmegaPtMCGen = new THnSparseF("fHistoElePtvsOmegaPtMCGen","",3,bins_eleptvsomegapt,xmin_eleptvsomegapt,xmax_eleptvsomegapt);
  fOutputAll->Add(fHistoElePtvsOmegaPtMCGen);

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
  fHistoOmegaMassvsPt=new TH2F("fHistoOmegaMassvsPt","Omega mass",100,1.67-0.05,1.67+0.05,20,0.,10.);
  fOutputAll->Add(fHistoOmegaMassvsPt);
  fHistoOmegaMassvsPtMCS=new TH2F("fHistoOmegaMassvsPtMCS","Omega mass",100,1.67-0.05,1.67+0.05,20,0.,10.);
  fOutputAll->Add(fHistoOmegaMassvsPtMCS);
  fHistoOmegaMassvsPtMCGen=new TH2F("fHistoOmegaMassvsPtMCGen","Omega mass",100,1.67-0.05,1.67+0.05,20,0.,10.);
  fOutputAll->Add(fHistoOmegaMassvsPtMCGen);

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

  fHistonEvtvsRunNumber=new TH1F("fHistonEvtvsRunNumber","",5000,-0.5,4999.5);
  fOutputAll->Add(fHistonEvtvsRunNumber);
  fHistonElevsRunNumber=new TH1F("fHistonElevsRunNumber","",5000,-0.5,4999.5);
  fOutputAll->Add(fHistonElevsRunNumber);
  fHistonOmegavsRunNumber=new TH1F("fHistonOmegavsRunNumber","",5000,-0.5,4999.5);
  fOutputAll->Add(fHistonOmegavsRunNumber);

  return;
}

//________________________________________________________________________
AliAODRecoCascadeHF* AliAnalysisTaskSEOmegac2eleOmegafromAODtracks::MakeCascadeHF(AliAODcascade *casc, AliAODTrack *part, AliAODEvent * aod, AliAODVertex *secVert, Bool_t mixing) 
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
AliAODVertex* AliAnalysisTaskSEOmegac2eleOmegafromAODtracks::CallPrimaryVertex(AliAODcascade *casc, AliAODTrack *trk, AliAODEvent* aod)
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
AliAODVertex* AliAnalysisTaskSEOmegac2eleOmegafromAODtracks::PrimaryVertex(const TObjArray *trkArray,
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
AliAODVertex* AliAnalysisTaskSEOmegac2eleOmegafromAODtracks::ReconstructSecondaryVertex(AliAODcascade *casc, AliAODTrack *part, AliAODEvent * aod) 
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
Int_t AliAnalysisTaskSEOmegac2eleOmegafromAODtracks::MatchToMC(AliAODRecoCascadeHF *exobj, TClonesArray *mcArray, Int_t *pdgarray_ele, Int_t *pdgarray_casc, Int_t *labelarray_ele, Int_t *labelarray_casc,  Int_t &ngen_ele, Int_t &ngen_casc) 
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

	Int_t pdgDgcasc[2]={321,3122};
	Int_t pdgDgv0[2]={2212,211};
  Int_t labcasc = MatchToMCCascade(theCascade,3334,pdgDgcasc,pdgDgv0,mcArray); // the cascade
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
Int_t AliAnalysisTaskSEOmegac2eleOmegafromAODtracks::MatchToMCCascade(AliAODcascade *theCascade, Int_t pdgabscasc, Int_t *pdgDgcasc, Int_t *pdgDgv0, TClonesArray *mcArray) const // the cascade
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
void AliAnalysisTaskSEOmegac2eleOmegafromAODtracks::SelectTrack( const AliVEvent *event, Int_t trkEntries, Int_t &nSeleTrks,Bool_t *seleFlags, TClonesArray *mcArray)
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
		fHistoElectronTPCPID->Fill(aodt->Pt(),nsigma_tpcele);
		fHistoElectronTOFPID->Fill(aodt->Pt(),nsigma_tofele);
		if(fabs(nsigma_tofele)<3.){
			fHistoElectronTPCPIDSelTOF->Fill(aodt->Pt(),nsigma_tpcele);
			Double_t eleeta = aodt->Eta();
			if(fabs(eleeta)<0.6)
				fHistoElectronTPCPIDSelTOFSmallEta->Fill(aodt->Pt(),nsigma_tpcele);
			if(fabs(eleeta)>0.6 && fabs(eleeta)<0.8)
				fHistoElectronTPCPIDSelTOFLargeEta->Fill(aodt->Pt(),nsigma_tpcele);
		}
  } // end loop on tracks
}
//________________________________________________________________________
void AliAnalysisTaskSEOmegac2eleOmegafromAODtracks::SelectCascade( const AliVEvent *event,Int_t nCascs,Int_t &nSeleCasc, Bool_t *seleCascFlags, TClonesArray *mcArray)
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

				FillCascROOTObjects(casc, mcArray);
      }
    }
}
//_________________________________________________________________
Int_t AliAnalysisTaskSEOmegac2eleOmegafromAODtracks::GetPoolIndex(Double_t zvert, Double_t mult){
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
void AliAnalysisTaskSEOmegac2eleOmegafromAODtracks::ResetPool(Int_t poolIndex){
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
void AliAnalysisTaskSEOmegac2eleOmegafromAODtracks::DoEventMixingWithPools(Int_t poolIndex,AliAODEvent *aodEvent, Bool_t *seleFlags)
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
//				printf("AliAnalysisTaskSEOmegac2eleOmegafromAODtracks::DoMixingWithPools ERROR: read event does not match to the stored one\n");
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
void AliAnalysisTaskSEOmegac2eleOmegafromAODtracks::MakeMCAnalysis(TClonesArray *mcArray)
{
	//
  // Analyze AliAODmcparticle
	//

	Int_t nmcpart = mcArray->GetEntriesFast();
	for(Int_t i=0;i<nmcpart;i++)
	{
		AliAODMCParticle *mcpart = (AliAODMCParticle*) mcArray->At(i);
		if(TMath::Abs(mcpart->GetPdgCode())==4332){
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
				if(TMath::Abs(mcdau->GetPdgCode())==3334){
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
		if(TMath::Abs(mcpart->GetPdgCode())==3334){
			Double_t etamin, etamax, rapmin, rapmax;
			fAnalCuts->GetProdCascEtaRange(etamin,etamax);
			fAnalCuts->GetProdCascRapRange(rapmin,rapmax);

			if((fabs(mcpart->Y())<rapmax) && (fabs(mcpart->Eta())<etamax))
				fHistoOmegaMassvsPtMCGen->Fill(1.67245, mcpart->Pt());
		}
	}
	return;
}
