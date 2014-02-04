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
//               Base class for Lc2V0 Analysis
//
//
//  The Lc spectra study is done 2D histograms:
//   cascade_invariantMass VS cascade_pT
//
//  Cuts have been centralized in AliRDHFCutsLctoV0 class
//
//-------------------------------------------------------------------------
//
//                 Authors: A.De Caro(a,b), P. Pagano(b)
//  (a) Centro 'E.Fermi' - Roma
//  (b) INFN and University of Salerno
//
//  Contatcs: decaro@sa.infn.it
//            paola.pagano@sa.infn.it
//-------------------------------------------------------------------------

#include <TSystem.h>
#include <TParticle.h>
#include <TParticlePDG.h>
#include <TH1F.h>
#include <TH1F.h>
#include <TH2F.h>
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
#include "AliAODRecoDecay.h"
#include "AliAODRecoDecayHF.h"
#include "AliAODRecoCascadeHF.h"
#include "AliAnalysisVertexingHF.h"
#include "AliESDtrack.h"
#include "AliAODTrack.h"
#include "AliAODv0.h"
#include "AliAODMCParticle.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisTaskSELc2V0bachelor.h"
#include "AliNormalizationCounter.h"
#include "AliAODPidHF.h"
#include "AliPIDResponse.h"
#include "AliTOFPIDResponse.h"
#include "AliInputEventHandler.h"
#include "AliESDtrackCuts.h"
#include "AliNeutralTrackParam.h"

using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskSELc2V0bachelor)

//__________________________________________________________________________
AliAnalysisTaskSELc2V0bachelor::AliAnalysisTaskSELc2V0bachelor() : AliAnalysisTaskSE(),
  fUseMCInfo(kFALSE),
  fOutput(0),
  fOutputAll(0),
  fOutputPIDBach(0),
  fCEvents(0),
  fPIDResponse(0),
  fIsK0SAnalysis(kFALSE),
  fCounter(0),
  fAnalCuts(0),
  fUseOnTheFlyV0(kFALSE),
  fIsEventSelected(kFALSE),
  fWriteVariableTree(kFALSE),
  fVariablesTree(0),
  fCandidateVariables(),
  fVtx1(0),
  fBzkG(0),
  fAdditionalChecks(kFALSE)
{
  //
  // Default ctor
  //
}
//___________________________________________________________________________
AliAnalysisTaskSELc2V0bachelor::AliAnalysisTaskSELc2V0bachelor(const Char_t* name,
							       AliRDHFCutsLctoV0* analCuts, Bool_t useOnTheFly,
							       Bool_t writeVariableTree, Bool_t additionalChecks) :
  AliAnalysisTaskSE(name),
  fUseMCInfo(kFALSE),
  fOutput(0),
  fOutputAll(0),
  fOutputPIDBach(0),
  fCEvents(0),
  fPIDResponse(0),
  fIsK0SAnalysis(kFALSE),
  fCounter(0),
  fAnalCuts(analCuts),
  fUseOnTheFlyV0(useOnTheFly),
  fIsEventSelected(kFALSE),
  fWriteVariableTree(writeVariableTree),
  fVariablesTree(0),
  fCandidateVariables(),
  fVtx1(0),
  fBzkG(0),
  fAdditionalChecks(additionalChecks)
{
  //
  // Constructor. Initialization of Inputs and Outputs
  //
  Info("AliAnalysisTaskSELc2V0bachelor","Calling Constructor");

  DefineOutput(1,TList::Class());  //conters
  DefineOutput(2,AliNormalizationCounter::Class());
  DefineOutput(3,AliRDHFCutsLctoV0::Class());
  if (!writeVariableTree) {
    DefineOutput(4,TList::Class());  //All Entries output
    DefineOutput(5,TList::Class());  //3sigma PID output
  } else {
    // Output slot #4 keeps a tree of the candidate variables after track selection
    DefineOutput(4,TTree::Class());  //My private output
  }

}

//___________________________________________________________________________
AliAnalysisTaskSELc2V0bachelor::~AliAnalysisTaskSELc2V0bachelor() {
  //
  // destructor
  //
  Info("~AliAnalysisTaskSELc2V0bachelor","Calling Destructor");
  
  if (fOutput) {
    delete fOutput;
    fOutput = 0;
  }

  if (fOutputAll) {
    delete fOutputAll;
    fOutputAll = 0;
  }

  if (fOutputPIDBach) {
    delete fOutputPIDBach;
    fOutputPIDBach = 0;
  }

  if (fPIDResponse) {
    delete  fPIDResponse;
  }

  if (fCounter) {
    delete fCounter;
    fCounter = 0;
  }

  if (fAnalCuts) {
    delete fAnalCuts;
    fAnalCuts = 0;
  }

  if (fVariablesTree) {
    delete fVariablesTree;
    fVariablesTree = 0;
  }

}
//_________________________________________________
void AliAnalysisTaskSELc2V0bachelor::Init() {
  //
  // Initialization
  //

  fIsEventSelected=kFALSE;

  if (fDebug > 1) AliInfo("Init");

  PostData(3,fAnalCuts);

  return;
}

//_________________________________________________
void AliAnalysisTaskSELc2V0bachelor::UserExec(Option_t *)
{
  // user exec
  if (!fInputEvent) {
    AliError("NO EVENT FOUND!");
    return;
  }

  AliAODEvent* aodEvent = dynamic_cast<AliAODEvent*>(fInputEvent);
  TClonesArray *arrayLctopKos=0;

  if (!aodEvent && AODEvent() && IsStandardAOD()) {
    // In case there is an AOD handler writing a standard AOD, use the AOD 
    // event in memory rather than the input (ESD) event.    
    aodEvent = dynamic_cast<AliAODEvent*> (AODEvent());
    // in this case the braches in the deltaAOD (AliAOD.VertexingHF.root)
    // have to taken from the AOD event hold by the AliAODExtension
    AliAODHandler* aodHandler = (AliAODHandler*) 
      ((AliAnalysisManager::GetAnalysisManager())->GetOutputEventHandler());

    if (aodHandler->GetExtensions()) {
      AliAODExtension *ext = (AliAODExtension*)aodHandler->GetExtensions()->FindObject("AliAOD.VertexingHF.root");
      AliAODEvent *aodFromExt = ext->GetAOD();
      arrayLctopKos=(TClonesArray*)aodFromExt->GetList()->FindObject("CascadesHF");
    }
  } else {
    arrayLctopKos=(TClonesArray*)aodEvent->GetList()->FindObject("CascadesHF");
  }

  fCEvents->Fill(1);

  if (fUseMCInfo)
    fAnalCuts->SetTriggerClass("");

  // AOD primary vertex
  fVtx1 = (AliAODVertex*)aodEvent->GetPrimaryVertex();
  if (!fVtx1) return;

  fIsEventSelected = fAnalCuts->IsEventSelected(aodEvent); // better to initialize before CheckEventSelection call

  CheckEventSelection(aodEvent);


  // fix for temporary bug in ESDfilter 
  fBzkG = (Double_t)aodEvent->GetMagneticField(); 
  if (TMath::Abs(fBzkG)<0.001) return;
  fCEvents->Fill(2);

  if (!arrayLctopKos) {
    AliInfo("Could not find array of HF cascades, skipping the event");
    return;
  } else {
    if (arrayLctopKos->GetEntriesFast()) {
      AliInfo(Form("Found %d cascades",arrayLctopKos->GetEntriesFast()));
    }
  }
  fCEvents->Fill(3);

  // mc analysis 
  TClonesArray *mcArray = 0;
  AliAODMCHeader *mcHeader=0;

  if (fUseMCInfo) {
    // MC array need for maching
    mcArray = dynamic_cast<TClonesArray*>(aodEvent->FindListObject(AliAODMCParticle::StdBranchName()));
    if (!mcArray) {
      AliError("Could not find Monte-Carlo in AOD");
      return;
    }
    fCEvents->Fill(4); // in case of MC events

    // load MC header
    mcHeader = (AliAODMCHeader*)aodEvent->GetList()->FindObject(AliAODMCHeader::StdBranchName());
    if (!mcHeader) {
      AliError("AliAnalysisTaskSELc2V0bachelor::UserExec: MC header branch not found!\n");
      return;
    }
    fCEvents->Fill(5); // in case of MC events
  }

  fCounter->StoreEvent(aodEvent,fAnalCuts,fUseMCInfo); // it is very important that it stays BEFORE any other event selection

  if (fVtx1->GetNContributors()>0) // this check is done in IsEventSelected
    fCEvents->Fill(6);

  if ( !fIsEventSelected ) return; // don't take into account not selected events 
  fCEvents->Fill(7);

  Int_t nSelectedAnal = 0;
  if (fIsK0SAnalysis) {
    MakeAnalysisForLc2prK0S(arrayLctopKos,mcArray, nSelectedAnal, fAnalCuts);

    if (nSelectedAnal)
      CheckEventSelectionWithCandidates(aodEvent);

  }

  fCounter->StoreCandidates(aodEvent,nSelectedAnal,kTRUE);
  fCounter->StoreCandidates(aodEvent,nSelectedAnal,kFALSE);

  PostData(1,fOutput);
  PostData(2,fCounter);
  if (!fWriteVariableTree) {
    PostData(4,fOutputAll);
    PostData(5,fOutputPIDBach);
  } else {
    PostData(4,fVariablesTree);
  }

  fIsEventSelected=kFALSE;

  return;
}

//________________________________________ terminate ___________________________
void AliAnalysisTaskSELc2V0bachelor::Terminate(Option_t*)
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
  
  //fCEvents = dynamic_cast<TH1F*>(fOutput->FindObject("fCEvents"));
  if (!fWriteVariableTree) {
    fOutputAll = dynamic_cast<TList*> (GetOutputData(4));
    if (!fOutputAll) {
      AliError("fOutputAll not available");
      return;
    }

    fOutputPIDBach = dynamic_cast<TList*> (GetOutputData(5));
    if (!fOutputPIDBach) {
      AliError("fOutputPIDBach not available");
      return;
    }
  }

  return;
}
//___________________________________________________________________________
void AliAnalysisTaskSELc2V0bachelor::UserCreateOutputObjects() { 
  // output
  AliInfo(Form("CreateOutputObjects of task %s\n", GetName()));
  
  AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
  fPIDResponse = inputHandler->GetPIDResponse();

  if (fAnalCuts->GetIsUsePID()){
    fAnalCuts->GetPidHF()->SetPidResponse(fPIDResponse);
    fAnalCuts->GetPidV0pos()->SetPidResponse(fPIDResponse);
    fAnalCuts->GetPidV0neg()->SetPidResponse(fPIDResponse);
    fAnalCuts->GetPidHF()->SetOldPid(kFALSE);
    fAnalCuts->GetPidV0pos()->SetOldPid(kFALSE);
    fAnalCuts->GetPidV0neg()->SetOldPid(kFALSE);
  }

  fOutput = new TList();
  fOutput->SetOwner();
  fOutput->SetName("chist0");
  DefineGeneralHistograms(); // define general histograms
  PostData(1,fOutput);

  fCounter = new AliNormalizationCounter("NormalizationCounter");
  fCounter->Init();
  PostData(2,fCounter);

  if (!fWriteVariableTree) {

    fOutputAll = new TList();
    fOutputAll->SetOwner();
    fOutputAll->SetName("listAll");

    fOutputPIDBach = new TList();
    fOutputPIDBach->SetOwner();
    fOutputPIDBach->SetName("listPIDBach");

    DefineAnalysisHistograms(); // define analysis histograms
  
    PostData(4,fOutputAll);
    PostData(5,fOutputPIDBach);
  }
  else {
    DefineTreeVariables();
    PostData(4,fVariablesTree);
  }

  return;
}

//-------------------------------------------------------------------------------
void AliAnalysisTaskSELc2V0bachelor::MakeAnalysisForLc2prK0S(TClonesArray *arrayLctopKos,
							     TClonesArray *mcArray,
							     Int_t &nSelectedAnal,
							     AliRDHFCutsLctoV0 *cutsAnal)
{

  // make the analysis

  Int_t pdgCand = 4122;
  Int_t pdgDgLctoV0bachelor[2]={2212,310}; // always 1st bachelor, 2nd V0
  Int_t pdgDgV0toDaughters[2]={211,211};

  // loop over cascades to search for candidates Lc->p+K0S
  Int_t nCascades= arrayLctopKos->GetEntriesFast();
  if (nCascades==0) {
    AliInfo("Could not find cascades, skipping the event");
    return;
  }

  for (Int_t iLctopK0S = 0; iLctopK0S<nCascades; iLctopK0S++) {

    ((TH1F*)(fOutput->FindObject("hCandidateSelection")))->Fill(0);

    // Lc candidates and K0S from Lc
    AliAODRecoCascadeHF* lcK0Spr = dynamic_cast<AliAODRecoCascadeHF*>(arrayLctopKos->At(iLctopK0S));
    if (!lcK0Spr) {
      AliDebug(2,Form("Cascade %d doens't exist, skipping",iLctopK0S));
      continue;
    }

    if (!lcK0Spr->GetSecondaryVtx()) {
      AliInfo("No secondary vertex"); // it will be done in AliRDHFCutsLctoV0::IsSelected
      continue;
    }

    if (lcK0Spr->GetNDaughters()!=2) {
      AliDebug(2,Form("Cascade %d has not 2 daughters (nDaughters=%d)",iLctopK0S,lcK0Spr->GetNDaughters())); // it will be done in AliRDHFCutsLctoV0::IsSelected
      continue;
    }

    AliAODv0 * v0part = dynamic_cast<AliAODv0*>(lcK0Spr->Getv0());
    AliAODTrack * bachPart = dynamic_cast<AliAODTrack*>(lcK0Spr->GetBachelor());
    if (!v0part || !bachPart) {
      AliDebug(2,Form("Cascade %d has no V0 or no bachelor object",iLctopK0S)); // it will be done in AliRDHFCutsLctoV0::IsSelected
      continue;
    }

    if (!v0part->GetSecondaryVtx()) {
      AliDebug(2,Form("No secondary vertex for V0 by cascade %d",iLctopK0S)); // it will be done in AliRDHFCutsLctoV0::IsSelected
      continue;
    }

    if (v0part->GetNDaughters()!=2) {
      AliDebug(2,Form("current V0 has not 2 daughters (onTheFly=%d, nDaughters=%d)",v0part->GetOnFlyStatus(),v0part->GetNDaughters())); // it will be done in AliRDHFCutsLctoV0::IsSelected
      continue;
    }

    AliAODTrack * v0Pos = dynamic_cast<AliAODTrack*>(lcK0Spr->Getv0PositiveTrack());
    AliAODTrack * v0Neg = dynamic_cast<AliAODTrack*>(lcK0Spr->Getv0NegativeTrack());
    if (!v0Neg || !v0Neg) {
      AliDebug(2,Form("V0 by cascade %d has no V0positive of V0negative object",iLctopK0S)); // it will be done in AliRDHFCutsLctoV0::IsSelected
      continue;
    }

    ((TH1F*)(fOutput->FindObject("hCandidateSelection")))->Fill(1);

    if (v0Pos->Charge() == v0Neg->Charge()) continue;

    ((TH1F*)(fOutput->FindObject("hCandidateSelection")))->Fill(2);

    Int_t isLc = 0;

    if (fUseMCInfo) {

      Int_t pdgCode=-2;

      // find associated MC particle for Lc -> p+K0 and K0S->pi+pi
      Int_t mcLabel = lcK0Spr->MatchToMC(pdgCand,pdgDgLctoV0bachelor[1],pdgDgLctoV0bachelor,pdgDgV0toDaughters,mcArray,kTRUE);
      if (mcLabel!=-1) {
	AliDebug(2,Form(" ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~cascade number %d (total cascade number = %d)", iLctopK0S,nCascades));

	AliAODMCParticle *partLc = dynamic_cast<AliAODMCParticle*>(mcArray->At(mcLabel));
	if (partLc) {
	  pdgCode = partLc->GetPdgCode();
	  if (pdgCode<0) AliDebug(2,Form(" ¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯ MClabel=%d ~~~~~~~~~~ pdgCode=%d", mcLabel, pdgCode));
	  pdgCode = TMath::Abs(pdgCode);
	  isLc = 1;
	}
      } else {
	AliDebug(2,Form("No candidate (cascade number %d -total cascade number = %d -)", iLctopK0S,nCascades));
	pdgCode=-1;
      }
    }

    FillLc2pK0Sspectrum(lcK0Spr, isLc,
			nSelectedAnal, cutsAnal,
			mcArray);
  }

  AliDebug(2, Form("Found %d Reco particles that are Lc!!", nSelectedAnal));

  return;
}

//________________________________________________________________________
void AliAnalysisTaskSELc2V0bachelor::FillLc2pK0Sspectrum(AliAODRecoCascadeHF *part,
							 Int_t isLc,
							 Int_t &nSelectedAnal,
							 AliRDHFCutsLctoV0 *cutsAnal,
							 TClonesArray *mcArray)
{
  //
  // Fill histos for Lc -> K0S+proton
  //

  TString fillthis="";

  AliAODTrack *bachelor = (AliAODTrack*)part->GetBachelor();
  Double_t momBach  = bachelor->P();

  AliAODv0 * v0part = (AliAODv0*)part->Getv0();
  Bool_t onFlyV0 = v0part->GetOnFlyStatus(); // on-the-flight V0s

  Bool_t areCutsUsingPID = cutsAnal->GetIsUsePID();
  cutsAnal->SetUsePID(kFALSE);
  Bool_t isInCascadeWindow = (((cutsAnal->IsSelectedSingleCut(part,AliRDHFCuts::kCandidate,0))&(AliRDHFCutsLctoV0::kLcToK0Spr))==(AliRDHFCutsLctoV0::kLcToK0Spr)); // cut on Lc->p+K0S invMass
  cutsAnal->SetUsePID(areCutsUsingPID);

  //if ( !( !onFlyV0 || (onFlyV0 && fUseOnTheFlyV0) ) ) return;
  if ( onFlyV0 && !fUseOnTheFlyV0 ) return;

  if (fAdditionalChecks) CheckCandidatesAtDifferentLevels(part,cutsAnal);

  if ( !( ( (cutsAnal->IsSelected(part,AliRDHFCuts::kTracks))&(AliRDHFCutsLctoV0::kLcToK0Spr)) == (AliRDHFCutsLctoV0::kLcToK0Spr) ) ) return;

  if ( !(cutsAnal->IsInFiducialAcceptance(part->Pt(),part->Y(4122))) ) return;

  if ( ( ( (cutsAnal->IsSelected(part,AliRDHFCuts::kAll))&(AliRDHFCutsLctoV0::kLcToK0Spr)) == (AliRDHFCutsLctoV0::kLcToK0Spr) ) ) nSelectedAnal++;

  // Fill candidate variable Tree (track selection, V0 invMass selection)
  if ( fWriteVariableTree ) {
    Double_t invmassK0S = v0part->MassK0Short();
    Double_t mk0sPDG = TDatabasePDG::Instance()->GetParticle(310)->Mass();
    if ( !onFlyV0 && isInCascadeWindow && part->CosV0PointingAngle()>0.99 && TMath::Abs(invmassK0S-mk0sPDG)<=0.05)
      FillTheTree(part,cutsAnal,mcArray,isLc);
    return;
  }

  cutsAnal->SetUsePID(kFALSE);
  Bool_t isCandidateSelectedCuts = (((cutsAnal->IsSelected(part,AliRDHFCuts::kCandidate))&(AliRDHFCutsLctoV0::kLcToK0Spr))==(AliRDHFCutsLctoV0::kLcToK0Spr)); // kinematic/topological cuts
  cutsAnal->SetUsePID(areCutsUsingPID);
  Bool_t isBachelorID = (((cutsAnal->IsSelected(part,AliRDHFCuts::kPID))&(AliRDHFCutsLctoV0::kLcToK0Spr))==(AliRDHFCutsLctoV0::kLcToK0Spr)); // ID x bachelor

  //if (((cutsAnal->IsSelected(part,AliRDHFCuts::kCandidate))&AliRDHFCutsLctoV0::kLcToK0Spr)==AliRDHFCutsLctoV0::kLcToK0Spr) {
  if (((cutsAnal->IsSelected(part,AliRDHFCuts::kAll))&AliRDHFCutsLctoV0::kLcToK0Spr)==AliRDHFCutsLctoV0::kLcToK0Spr) {
    if (fUseMCInfo && isLc && !fWriteVariableTree) {
      Int_t pdgCand1 = 4122;
      Int_t pdgDgLctoV0bachelor1[2]={2212,310};
      Int_t pdgDgV0toDaughters1[2]={211,211};
      Int_t mcLabel1=part->MatchToMC(pdgCand1,pdgDgLctoV0bachelor1[1],pdgDgLctoV0bachelor1,pdgDgV0toDaughters1,mcArray,kTRUE);
      AliDebug(2,Form(" Found true MC candidate: Lc->pK0S(%d) - onTheFly=%1d",mcLabel1,onFlyV0));
    }
  }

  Double_t nSigmaTPCpr=-999.;
  cutsAnal->GetPidHF()->GetnSigmaTPC(bachelor,4,nSigmaTPCpr);
  Double_t nSigmaTOFpr=-999.;
  cutsAnal->GetPidHF()->GetnSigmaTOF(bachelor,4,nSigmaTOFpr);

  Double_t nSigmaTPCpi=-999.;
  cutsAnal->GetPidHF()->GetnSigmaTPC(bachelor,2,nSigmaTPCpi);
  Double_t nSigmaTOFpi=-999.;
  cutsAnal->GetPidHF()->GetnSigmaTOF(bachelor,2,nSigmaTOFpi);

  Double_t nSigmaTPCka=-999.;
  cutsAnal->GetPidHF()->GetnSigmaTPC(bachelor,3,nSigmaTPCka);
  Double_t nSigmaTOFka=-999.;
  cutsAnal->GetPidHF()->GetnSigmaTOF(bachelor,3,nSigmaTOFka);

  //if (onFlyV0 && fUseOnTheFlyV0) {  
  if (onFlyV0) {  

    if (isCandidateSelectedCuts)
      FillAnalysisHistograms(part,isBachelorID,"");

  }
  //else if (!onFlyV0) {
  else {

    if (isCandidateSelectedCuts) {

      fillthis="histoprotonBachSigmaVspTOF";
      ((TH2F*)(fOutput->FindObject(fillthis)))->Fill(momBach,nSigmaTOFpr);
      fillthis="histoprotonBachSigmaVspTPC";
      ((TH2F*)(fOutput->FindObject(fillthis)))->Fill(momBach,nSigmaTPCpr);

      FillAnalysisHistograms(part,isBachelorID,"Offline");

    }

  }


  if (fUseMCInfo) {
    if (isLc==1) {
      //if (onFlyV0 && fUseOnTheFlyV0) {
      if (onFlyV0) {

	if (isCandidateSelectedCuts)
	  FillAnalysisHistograms(part,isBachelorID,"Sgn");
    
      }     
      //else if (!onFlyV0) {  
      else {  

 	if (isCandidateSelectedCuts) {

	  fillthis="histoprotonBachSigmaVspTOFsgn";
	  ((TH2F*)(fOutput->FindObject(fillthis)))->Fill(momBach,nSigmaTOFpr);
	  fillthis="histoprotonBachSigmaVspTPCsgn";
	  ((TH2F*)(fOutput->FindObject(fillthis)))->Fill(momBach,nSigmaTPCpr);

	  FillAnalysisHistograms(part,isBachelorID,"OfflineSgn");

	}

      }

    }// sgn
    else { // bkg
      //if (onFlyV0 && fUseOnTheFlyV0) {
      if (onFlyV0) {

	if (isCandidateSelectedCuts)
	  FillAnalysisHistograms(part,isBachelorID,"Bkg");

      }
      //else if (!onFlyV0) {
      else {

	if (isCandidateSelectedCuts) {

	  fillthis="histoprotonBachSigmaVspTOFbkg";
	  ((TH2F*)(fOutput->FindObject(fillthis)))->Fill(momBach,nSigmaTOFpr);
	  fillthis="histoprotonBachSigmaVspTPCbkg";
	  ((TH2F*)(fOutput->FindObject(fillthis)))->Fill(momBach,nSigmaTPCpr);

	  FillAnalysisHistograms(part,isBachelorID,"OfflineBkg");
	}

      }

    }
  } // if fUseMCInfo
 
  return;
}

//----------------------------------------------------
void AliAnalysisTaskSELc2V0bachelor::DefineK0SHistos()
{ 

  TString nameMass=" ", nameSgn=" ", nameBkg=" ";

  Double_t mLcPDG  = TDatabasePDG::Instance()->GetParticle(4122)->Mass();
  Double_t mK0SPDG = TDatabasePDG::Instance()->GetParticle(310)->Mass();

  if (fUseOnTheFlyV0) {

    // V0 invariant masses (on-the-fly)
    nameMass="histK0SMass";
    TH2F* spectrumK0SMass = new TH2F(nameMass.Data(),"K^{0}_{S} invariant mass VS p_{T}; M(#pi^{+}#pi^{-}) [GeV/c^{2}]; p_{T}(K^{0}_{S}) [GeV/c]; Entries",
				    1000,mK0SPDG-0.050,mK0SPDG+0.050,175,0.,35.);

    // Lc invariant masses (x K0S on-the-fly)
    nameMass="histLcMassByK0S";
    TH2F* spectrumLcMassByK0S = new TH2F(nameMass.Data(),"#Lambda_{c} invariant mass (by K^{0}_{S}) vs p_{T} ; m_{inv}(p-K^{0}_{S}) [GeV/c^{2}]; p_{T}(#Lambda_{c}) [GeV/c]",
					 1000,mLcPDG-0.250,mLcPDG+0.250,175,0.,35.);

    nameMass="histpK0Svsp";
    TH2F* momentumDistributionK0Svsp = new TH2F(nameMass.Data(),"#Lambda_{c}: p(K^{0}_{S}) vs p(p);  p_{p}; p_{K^{0}_{S}}  ",
						175,0.,35.,175,0.,35.);

    nameMass="histArmPodK0S";
    TH2F* armenterosPodK0S = new TH2F(nameMass.Data(),"K^{0}_{S}: Armenteros-Podolanski distribution; #frac{p_{L}^{+}-p_{L}^{-}}{p_{L}^{+}+p_{L}^{-}}; p_{T}^{+} [GeV/c]",
				      200,-1.,1.,300,0.,0.3);
 
    nameMass="histDCAtoPVvspK0S";
    TH2F *dcatoPVvspK0S = new TH2F(nameMass.Data(),"K^{0}_{S}: DCA to Primary Vertex vs K^{0}_{S} momentum ; p(K^{0}_{S}) [GeV/c]; DCA to Primary Vertex [n#sigmas]; Entries",
				   175,0.,35.,50,0.,5.);

    nameMass="histK0ScosPAwrtPVvspK0S";
    TH2F *cosPAwrtPVvspK0S = new TH2F(nameMass.Data(),"K^{0}_{S}: cosine of pointing angle wrt primary vertex vs K^{0}_{S} momentum ; p(K^{0}_{S}) [GeV/c]; cosine; Entries",
				      175,0.,35.,100,0.99,1.);

    TH2F* allspectrumK0SMass = (TH2F*)spectrumK0SMass->Clone(); 
    TH2F* allspectrumLcMassByK0S    = (TH2F*)spectrumLcMassByK0S->Clone(); 
    TH2F* allmomentumDistributionK0Svsp= (TH2F*)momentumDistributionK0Svsp->Clone(); 
    TH2F* alldcatoPVvspK0S=(TH2F*)dcatoPVvspK0S->Clone(); 
    TH2F* allcosV0PAwrtPVvspK0S=(TH2F*)cosPAwrtPVvspK0S->Clone(); 

    TH2F* pidBachspectrumK0SMass = (TH2F*)spectrumK0SMass->Clone(); 
    TH2F* pidBachspectrumLcMassByK0S    = (TH2F*)spectrumLcMassByK0S->Clone(); 
    TH2F* pidBachmomentumDistributionK0Svsp= (TH2F*)momentumDistributionK0Svsp->Clone(); 
    TH2F* pidBachdcatoPVvspK0S=(TH2F*)dcatoPVvspK0S->Clone(); 
    TH2F* pidBachcosV0PAwrtPVvspK0S=(TH2F*)cosPAwrtPVvspK0S->Clone(); 

    TH2F* allArmenterosPodK0S = (TH2F*)armenterosPodK0S->Clone();
    TH2F* pidBachArmenterosPodK0S = (TH2F*)armenterosPodK0S->Clone();

    fOutputAll->Add(allspectrumK0SMass);
    fOutputAll->Add(allspectrumLcMassByK0S);
    fOutputAll->Add(allmomentumDistributionK0Svsp); 
    fOutputAll->Add(allArmenterosPodK0S);
    fOutputAll->Add(alldcatoPVvspK0S);
    fOutputAll->Add(allcosV0PAwrtPVvspK0S);

    fOutputPIDBach->Add(pidBachspectrumK0SMass);
    fOutputPIDBach->Add(pidBachspectrumLcMassByK0S);
    fOutputPIDBach->Add(pidBachmomentumDistributionK0Svsp); 
    fOutputPIDBach->Add(pidBachArmenterosPodK0S);
    fOutputPIDBach->Add(pidBachdcatoPVvspK0S);
    fOutputPIDBach->Add(pidBachcosV0PAwrtPVvspK0S);
 
  }

  // V0 invariant masses (offline)
  nameMass="histK0SMassOffline";
  TH2F* spectrumK0SMassOffline = new TH2F(nameMass.Data(),"K^{0}_{S} invariant mass VS p_{T}; M(#pi^{+}#pi^{-}) [GeV/c^{2}]; p_{T}(K^{0}_{S}) [GeV/c]; Entries",
					 1000,mK0SPDG-0.050,mK0SPDG+0.050,175,0.,35.);

  // Lc invariant masses (x K0S offline)
  nameMass="histLcMassByK0SOffline";
  TH2F* spectrumLcMassOfflineByK0S = new TH2F(nameMass.Data(),"#Lambda_{c} invariant mass (by K^{0}_{S}) vs p_{T}; M(K^{0}_{S}p) [GeV/c^{2}]; p_{T} [GeV/c]",
					      1000,mLcPDG-0.250,mLcPDG+0.250,175,0.,35.);

  nameMass="histpK0SvspOffline";
  TH2F* momentumDistributionK0SvspOffline = new TH2F(nameMass.Data(),"#Lambda_{c}: p(K^{0}_{S}) vs p(p) - Offline ;  p_{p} [GeV/c]; p_{K^{0}_{S}} [GeV/c]",
						     175,0.,35.,175,0.,35.);

  nameMass="histArmPodK0SOffline";
  TH2F* armenterosPodK0SOff = new TH2F(nameMass.Data(),"K^{0}_{S}  Armenteros-Podolanski distribution - Offline; #frac{p_{L}^{+}-p_{L}^{-}}{p_{L}^{+}+p_{L}^{-}}; p_{T}^{+} [GeV/c]",
				       200,-1.,1.,300,0.,0.3);

  nameMass="histDCAtoPVvspK0SOffline";
  TH2F *dcatoPVvspK0SOffline = new TH2F(nameMass.Data(),"K^{0}_{S}: DCA to Primary Vertex vs  K^{0}_{S} invariant mass - Offline; p(K^{0}_{S}) [GeV/c]; DCA to Primary Vertex [n#sigmas]; Entries",
					175,0.,35.,50,0.,5.);

  nameMass="histK0ScosPAwrtPVvspK0SOffline";
  TH2F *cosPAwrtPVvspK0SOffline = new TH2F(nameMass.Data(),"K^{0}_{S}: cosine of pointing angle wrt primary vertex vs K^{0}_{S} momentum - Offline; p(K^{0}_{S}) [GeV/c]; cosine; Entries",
					   175,0.,35.,100,0.99,1.);



  TH2F* allspectrumK0SMassOffline = (TH2F*)spectrumK0SMassOffline->Clone(); 
  TH2F* allspectrumLcMassOfflineByK0S    = (TH2F*)spectrumLcMassOfflineByK0S->Clone(); 
  TH2F* allmomentumDistributionK0SvspOffline= (TH2F*)momentumDistributionK0SvspOffline->Clone(); 
  TH2F* alldcatoPVvspK0SOffline=(TH2F*)dcatoPVvspK0SOffline->Clone(); 
  TH2F* allcosPAwrtPVvspK0SOffline=(TH2F*)cosPAwrtPVvspK0SOffline->Clone(); 
  TH2F* allArmenterosPodK0SOff = (TH2F*)armenterosPodK0SOff->Clone();

  TH2F* pidBachspectrumK0SMassOffline = (TH2F*)spectrumK0SMassOffline->Clone(); 
  TH2F* pidBachspectrumLcMassOfflineByK0S    = (TH2F*)spectrumLcMassOfflineByK0S->Clone(); 
  TH2F* pidBachmomentumDistributionK0SvspOffline= (TH2F*)momentumDistributionK0SvspOffline->Clone(); 
  TH2F* pidBachdcatoPVvspK0SOffline=(TH2F*)dcatoPVvspK0SOffline->Clone(); 
  TH2F* pidBachcosPAwrtPVvspK0SOffline=(TH2F*)cosPAwrtPVvspK0SOffline->Clone(); 
  TH2F* pidBachArmenterosPodK0SOff = (TH2F*)armenterosPodK0SOff->Clone();


  fOutputAll->Add(allspectrumK0SMassOffline);
  fOutputAll->Add(allspectrumLcMassOfflineByK0S);
  fOutputAll->Add(allmomentumDistributionK0SvspOffline); 
  fOutputAll->Add(allArmenterosPodK0SOff);
  fOutputAll->Add(alldcatoPVvspK0SOffline);
  fOutputAll->Add(allcosPAwrtPVvspK0SOffline);

  fOutputPIDBach->Add(pidBachspectrumK0SMassOffline);
  fOutputPIDBach->Add(pidBachspectrumLcMassOfflineByK0S);
  fOutputPIDBach->Add(pidBachmomentumDistributionK0SvspOffline); 
  fOutputPIDBach->Add(pidBachArmenterosPodK0SOff);
  fOutputPIDBach->Add(pidBachdcatoPVvspK0SOffline);
  fOutputPIDBach->Add(pidBachcosPAwrtPVvspK0SOffline);

  if (fUseMCInfo) {

    if (fUseOnTheFlyV0) {

      nameSgn="histK0SMassSgn";
      nameBkg="histK0SMassBkg";
      TH2F* spectrumK0SMassSgn = new TH2F(nameSgn.Data(), "K^{0}_{S} - sgn: invariant mass VS p_{T} - MC; m_{inv}(#pi^{+}#pi^{-}) [GeV/c^{2}]; p_{T}(K^{0}_{S}) [GeV/c]; Entries",
					  1000,mK0SPDG-0.050,mK0SPDG+0.050,175,0.,35.);
      TH2F* spectrumK0SMassBkg = new TH2F(nameBkg.Data(), "K^{0}_{S} - bkg: invariant mass VS p_{T} - MC; m_{inv}(#pi^{+}#pi^{-}) [GeV/c^{2}]; p_{T}(K^{0}_{S}) [GeV/c]; Entries",
					  1000,mK0SPDG-0.050,mK0SPDG+0.050,175,0.,35.);

      TH2F* allspectrumK0SMassSgn = (TH2F*)spectrumK0SMassSgn->Clone(); 
      TH2F* allspectrumK0SMassBkg = (TH2F*) spectrumK0SMassBkg->Clone();  
      TH2F* pidBachspectrumK0SMassSgn = (TH2F*)spectrumK0SMassSgn->Clone(); 
      TH2F* pidBachspectrumK0SMassBkg = (TH2F*) spectrumK0SMassBkg->Clone();  
      fOutputAll->Add(allspectrumK0SMassSgn);
      fOutputAll->Add(allspectrumK0SMassBkg);
      fOutputPIDBach->Add(pidBachspectrumK0SMassSgn);
      fOutputPIDBach->Add(pidBachspectrumK0SMassBkg);


      nameSgn="histLcMassByK0SSgn";
      nameBkg="histLcMassByK0SBkg";
      TH2F* spectrumLcMassByK0SSgn = new TH2F(nameSgn.Data(), "#Lambda_{c} - sgn: invariant mass (by K^{0}_{S}) vs p_{T}  - MC; m_{inv}(p-K^{0}_{S}) [GeV/c^{2}];  p_{T}",
					      1000,mLcPDG-0.250,mLcPDG+0.250,175,0.,35.);
      TH2F* spectrumLcMassByK0SBkg = new TH2F(nameBkg.Data(), "#Lambda_{c} - bkg: invariant mass (by K^{0}_{S}) vs p_{T}  - MC; m_{inv}(p-K^{0}_{S}) [GeV/c^{2}]; p_{T}",
					      1000,mLcPDG-0.250,mLcPDG+0.250,175,0.,35.);

      TH2F* allspectrumLcMassByK0SSgn = (TH2F*)spectrumLcMassByK0SSgn->Clone(); 
      TH2F* allspectrumLcMassByK0SBkg = (TH2F*) spectrumLcMassByK0SBkg->Clone();  
      TH2F* pidBachspectrumLcMassByK0SSgn = (TH2F*)spectrumLcMassByK0SSgn->Clone(); 
      TH2F* pidBachspectrumLcMassByK0SBkg = (TH2F*) spectrumLcMassByK0SBkg->Clone();  
      fOutputAll->Add(allspectrumLcMassByK0SSgn);
      fOutputAll->Add(allspectrumLcMassByK0SBkg);
      fOutputPIDBach->Add(pidBachspectrumLcMassByK0SSgn);
      fOutputPIDBach->Add(pidBachspectrumLcMassByK0SBkg);

      nameSgn="histpK0SvspSgn";
      nameBkg="histpK0SvspBkg";
      TH2F* momentumDistributionK0SvspSgn= new TH2F(nameSgn.Data(),"#Lambda_{c} - sgn: K^{0}_{S} vs p Total Momentum Distribution - MC; p_{p}; p_{K^{0}_{S}}",
						    175,0.,35.,175,0.,35.);
      TH2F* momentumDistributionK0SvspBkg= new TH2F(nameBkg.Data(),"#Lambda_{c} - bkg: K^{0}_{S} vs p Total Momentum Distribution - MC; p_{p}; p_{K^{0}_{S}}",
						    175,0.,35.,175,0.,35.);

      TH2F* allmomentumDistributionK0SvspSgn= (TH2F*)momentumDistributionK0SvspSgn->Clone(); 
      TH2F* allmomentumDistributionK0SvspBkg= (TH2F*)momentumDistributionK0SvspBkg->Clone(); 
      TH2F* pidBachmomentumDistributionK0SvspSgn= (TH2F*)momentumDistributionK0SvspSgn->Clone(); 
      TH2F* pidBachmomentumDistributionK0SvspBkg= (TH2F*)momentumDistributionK0SvspBkg->Clone(); 
      fOutputAll->Add(allmomentumDistributionK0SvspSgn); 
      fOutputAll->Add(allmomentumDistributionK0SvspBkg); 
      fOutputPIDBach->Add(pidBachmomentumDistributionK0SvspSgn); 
      fOutputPIDBach->Add(pidBachmomentumDistributionK0SvspBkg); 


      // armenteros-podolanski plots K0S
      nameSgn="histArmPodK0SSgn";
      nameBkg="histArmPodK0SBkg";
      TH2F* armenterosPodK0SSgn = new TH2F(nameSgn.Data(),"K^{0}_{S}  Armenteros-Podolanski distribution (sgn); #frac{p_{L}^{+}-p_{L}^{-}}{p_{L}^{+}+p_{L}^{-}}; p_{T}^{+} [GeV/c]",
					   200,-1.,1.,300,0.,0.3);
      TH2F* armenterosPodK0SBkg = new TH2F(nameBkg.Data(),"K^{0}_{S}  Armenteros-Podolanski distribution (bkg); #frac{p_{L}^{+}-p_{L}^{-}}{p_{L}^{+}+p_{L}^{-}}; p_{T}^{+} [GeV/c]",
					   200,-1.,1.,300,0.,0.3);

      TH2F* allArmenterosPodK0SSgn = (TH2F*)armenterosPodK0SSgn->Clone();
      TH2F* allArmenterosPodK0SBkg = (TH2F*)armenterosPodK0SBkg->Clone();
      TH2F* pidBachArmenterosPodK0SSgn = (TH2F*)armenterosPodK0SSgn->Clone();
      TH2F* pidBachArmenterosPodK0SBkg = (TH2F*)armenterosPodK0SBkg->Clone();

      fOutputAll->Add(allArmenterosPodK0SSgn);
      fOutputAll->Add(allArmenterosPodK0SBkg);
      fOutputPIDBach->Add(pidBachArmenterosPodK0SSgn);
      fOutputPIDBach->Add(pidBachArmenterosPodK0SBkg);


      nameSgn="histDCAtoPVvspK0SSgn";
      nameBkg="histDCAtoPVvspK0SBkg";
      TH2F *dcatoPVvspK0SSgn=new TH2F(nameSgn.Data(),"K^{0}_{S} - sgn: DCA to Primary Vertex vs  K^{0}_{S} invariant mass (sgn); p(K^{0}_{S}) [GeV/c]; DCA to Primary Vertex [n#sigmas]; Entries",175,0.,35.,50,0.,5.);
      TH2F *dcatoPVvspK0SBkg=new TH2F(nameBkg.Data(),"K^{0}_{S} - bkg: DCA to Primary Vertex vs  K^{0}_{S} invariant mass (bkg); p(K^{0}_{S}) [GeV/c]; DCA to Primary Vertex [n#sigmas]; Entries",175,0.,35.,50,0.,5.);

      TH2F* alldcatoPVvspK0SSgn= (TH2F*)dcatoPVvspK0SSgn->Clone();
      TH2F* alldcatoPVvspK0SBkg= (TH2F*)dcatoPVvspK0SBkg->Clone();
      TH2F* pidBachdcatoPVvspK0SSgn= (TH2F*)dcatoPVvspK0SSgn->Clone();
      TH2F* pidBachdcatoPVvspK0SBkg= (TH2F*)dcatoPVvspK0SBkg->Clone();

      fOutputAll->Add(alldcatoPVvspK0SSgn);
      fOutputPIDBach->Add(pidBachdcatoPVvspK0SSgn);
      fOutputAll->Add(alldcatoPVvspK0SBkg);
      fOutputPIDBach->Add(pidBachdcatoPVvspK0SBkg);

    }


    nameSgn="histK0SMassOfflineSgn";
    nameBkg="histK0SMassOfflineBkg";
    TH2F* spectrumK0SMassOfflineSgn = new TH2F(nameSgn.Data(), "K^{0}_{S} - sgn: invariant mass VS p_{T} - MC; m_{inv}(#pi^{+}-#pi^{-}) [GeV/c^{2}]; p_{T}(K^{0}_{S}) [GeV/c]; Entries",
					      1000,mK0SPDG-0.050,mK0SPDG+0.050,175,0.,35.);
    TH2F* spectrumK0SMassOfflineBkg = new TH2F(nameBkg.Data(), "K^{0}_{S} - bkg: invariant mass VS p_{T} - MC; m_{inv}(#pi^{+}-#pi^{-}) [GeV/c^{2}]; p_{T}(K^{0}_{S}) [GeV/c]; Entries",
					      1000,mK0SPDG-0.050,mK0SPDG+0.050,175,0.,35.);

    TH2F* allspectrumK0SMassOfflineSgn = (TH2F*)spectrumK0SMassOfflineSgn->Clone(); 
    TH2F* allspectrumK0SMassOfflineBkg = (TH2F*) spectrumK0SMassOfflineBkg->Clone();  
    fOutputAll->Add(allspectrumK0SMassOfflineSgn);
    fOutputAll->Add(allspectrumK0SMassOfflineBkg);


    TH2F* pidBachspectrumK0SMassOfflineSgn = (TH2F*)spectrumK0SMassOfflineSgn->Clone(); 
    TH2F* pidBachspectrumK0SMassOfflineBkg = (TH2F*) spectrumK0SMassOfflineBkg->Clone();
    fOutputPIDBach->Add(pidBachspectrumK0SMassOfflineSgn);
    fOutputPIDBach->Add(pidBachspectrumK0SMassOfflineBkg);


    nameSgn="histLcMassByK0SOfflineSgn";
    nameBkg="histLcMassByK0SOfflineBkg";
    TH2F* spectrumLcMassOfflineByK0SSgn = new TH2F(nameSgn.Data(), "#Lambda_{c} - sgn: invariant mass (by K^{0}_{S})  vs p_{T} - MC; M(#Lambda_{c}) [GeV/c^{2}]; p_{T}",
						   1000,mLcPDG-0.250,mLcPDG+0.250,175,0.,35.);
    TH2F* spectrumLcMassOfflineByK0SBkg = new TH2F(nameBkg.Data(), "#Lambda_{c} - bkg: invariant mass (by K^{0}_{S})  vs p_{T} - MC; M(#Lambda_{c}) [GeV/c^{2}]; p_{T}",
						   1000,mLcPDG-0.250,mLcPDG+0.250,175,0.,35.);


    TH2F* allspectrumLcMassOfflineByK0SSgn = (TH2F*)spectrumLcMassOfflineByK0SSgn->Clone(); 
    TH2F* allspectrumLcMassOfflineByK0SBkg = (TH2F*) spectrumLcMassOfflineByK0SBkg->Clone();  
    TH2F* pidBachspectrumLcMassOfflineByK0SSgn = (TH2F*)spectrumLcMassOfflineByK0SSgn->Clone(); 
    TH2F* pidBachspectrumLcMassOfflineByK0SBkg = (TH2F*) spectrumLcMassOfflineByK0SBkg->Clone();  
    fOutputAll->Add(allspectrumLcMassOfflineByK0SSgn);
    fOutputAll->Add(allspectrumLcMassOfflineByK0SBkg);
    fOutputPIDBach->Add(pidBachspectrumLcMassOfflineByK0SSgn);
    fOutputPIDBach->Add(pidBachspectrumLcMassOfflineByK0SBkg);
  
 
    nameSgn="histpK0SvspOfflineSgn";
    nameBkg="histpK0SvspOfflineBkg";
    TH2F* momentumDistributionK0SvspOfflineSgn= new TH2F(nameSgn.Data(),"#Lambda_{c} - sgn: K^{0}_{S} vs p Total Momentum Distribution - Offline  - MC; p_{p}; p_{K^{0}_{S}}",
							 175,0.,35.,175,0.,35.);
    TH2F* momentumDistributionK0SvspOfflineBkg= new TH2F(nameBkg.Data(),"#Lambda_{c} - bkg: K^{0}_{S} vs p Total Momentum Distribution - Offline  - MC; p_{p}; p_{K^{0}_{S}}",
							 175,0.,35.,175,0.,35.);


    TH2F* allmomentumDistributionK0SvspOfflineSgn= (TH2F*)momentumDistributionK0SvspOfflineSgn->Clone(); 
    TH2F* allmomentumDistributionK0SvspOfflineBkg= (TH2F*)momentumDistributionK0SvspOfflineBkg->Clone(); 
    TH2F* pidBachmomentumDistributionK0SvspOfflineSgn= (TH2F*)momentumDistributionK0SvspOfflineSgn->Clone(); 
    TH2F* pidBachmomentumDistributionK0SvspOfflineBkg= (TH2F*)momentumDistributionK0SvspOfflineBkg->Clone(); 
    fOutputAll->Add(allmomentumDistributionK0SvspOfflineSgn); 
    fOutputAll->Add(allmomentumDistributionK0SvspOfflineBkg); 
    fOutputPIDBach->Add(pidBachmomentumDistributionK0SvspOfflineSgn); 
    fOutputPIDBach->Add(pidBachmomentumDistributionK0SvspOfflineBkg); 





    // armenteros-podolanski plots K0S (offline)
    nameSgn="histArmPodK0SOfflineSgn";
    nameBkg="histArmPodK0SOfflineBkg";
    TH2F* armenterosPodK0SOffSgn = new TH2F(nameSgn.Data(),"K^{0}_{S}  Armenteros-Podolanski distribution (sgn) -offline-; #frac{p_{L}^{+}-p_{L}^{-}}{p_{L}^{+}+p_{L}^{-}}; p_{T}^{+} [GeV/c]",
					    200,-1.,1.,300,0.,0.3);
    TH2F* armenterosPodK0SOffBkg = new TH2F(nameBkg.Data(),"K^{0}_{S}  Armenteros-Podolanski distribution (bkg) -offline-; #frac{p_{L}^{+}-p_{L}^{-}}{p_{L}^{+}+p_{L}^{-}}; p_{T}^{+} [GeV/c]",
					    200,-1.,1.,300,0.,0.3);
  

    TH2F* allArmenterosPodK0SOffSgn = (TH2F*)armenterosPodK0SOffSgn->Clone();
    TH2F* allArmenterosPodK0SOffBkg = (TH2F*)armenterosPodK0SOffBkg->Clone();
    TH2F* pidBachArmenterosPodK0SOffSgn = (TH2F*)armenterosPodK0SOffSgn->Clone();
    TH2F* pidBachArmenterosPodK0SOffBkg = (TH2F*)armenterosPodK0SOffBkg->Clone();


    fOutputAll->Add(allArmenterosPodK0SOffSgn);
    fOutputAll->Add(allArmenterosPodK0SOffBkg);
    fOutputPIDBach->Add(pidBachArmenterosPodK0SOffSgn);
    fOutputPIDBach->Add(pidBachArmenterosPodK0SOffBkg);


    nameSgn="histDCAtoPVvspK0SOfflineSgn";
    nameBkg="histDCAtoPVvspK0SOfflineBkg";
    TH2F *dcatoPVvspK0SOfflineSgn=new TH2F(nameSgn.Data(),"K^{0}_{S} -offline - (sgn): DCA to Primary Vertex vs  K^{0}_{S} invariant mass; p(K^{0}_{S}) [GeV/c]; DCA to Primary Vertex [n#sigmas]; Entries",
					   175,0.,35.,50,0.,5.);
    TH2F *dcatoPVvspK0SOfflineBkg=new TH2F(nameBkg.Data(),"K^{0}_{S} -offline - (bkg): DCA to Primary Vertex vs  K^{0}_{S} invariant mass; p(K^{0}_{S}) [GeV/c]; DCA to Primary Vertex [n#sigmas]; Entries",
					   175,0.,35.,50,0.,5.);
    

    TH2F* alldcatoPVvspK0SOfflineSgn= (TH2F*)dcatoPVvspK0SOfflineSgn->Clone();
    TH2F* pidBachdcatoPVvspK0SOfflineSgn= (TH2F*)dcatoPVvspK0SOfflineSgn->Clone();
    TH2F* alldcatoPVvspK0SOfflineBkg= (TH2F*)dcatoPVvspK0SOfflineBkg->Clone();
    TH2F* pidBachdcatoPVvspK0SOfflineBkg= (TH2F*)dcatoPVvspK0SOfflineBkg->Clone();



    fOutputAll->Add(alldcatoPVvspK0SOfflineSgn);
    fOutputPIDBach->Add(pidBachdcatoPVvspK0SOfflineSgn);
    fOutputAll->Add(alldcatoPVvspK0SOfflineBkg);
    fOutputPIDBach->Add(pidBachdcatoPVvspK0SOfflineBkg);

  }

  return;
}

//---------------------------
void AliAnalysisTaskSELc2V0bachelor::CheckEventSelection(AliAODEvent *aodEvent) {
  //
  // To fill control histograms
  //

  TClonesArray *arrayLctopKos=0;
  if (!aodEvent){
    if(AODEvent() && IsStandardAOD()) {
      // In case there is an AOD handler writing a standard AOD, use the AOD 
      // event in memory rather than the input (ESD) event.    
      aodEvent = dynamic_cast<AliAODEvent*> (AODEvent());
      // in this case the braches in the deltaAOD (AliAOD.VertexingHF.root)
      // have to taken from the AOD event hold by the AliAODExtension
      AliAODHandler* aodHandler = (AliAODHandler*) 
	((AliAnalysisManager::GetAnalysisManager())->GetOutputEventHandler());

      if (aodHandler->GetExtensions()) {
	AliAODExtension *ext = (AliAODExtension*)aodHandler->GetExtensions()->FindObject("AliAOD.VertexingHF.root");
	AliAODEvent *aodFromExt = ext->GetAOD();
	arrayLctopKos=(TClonesArray*)aodFromExt->GetList()->FindObject("CascadesHF");
      }
    }
  } else {
    arrayLctopKos=(TClonesArray*)aodEvent->GetList()->FindObject("CascadesHF");
  }

  Float_t zVertex = fVtx1->GetZ();
  TString titleVtx=fVtx1->GetTitle();

  if (TMath::Abs(fBzkG)>=0.001) {

    if (fAdditionalChecks) ((TH1F*)(fOutput->FindObject("hZ2")))->Fill(zVertex);

    if (arrayLctopKos) {

      if (fAdditionalChecks) ((TH1F*)(fOutput->FindObject("hZ3")))->Fill(zVertex);

      // mc analysis 
      TClonesArray *mcArray = 0;
      AliAODMCHeader *mcHeader=0;

      if (fUseMCInfo) {
	// MC array need for maching
	mcArray = dynamic_cast<TClonesArray*>(aodEvent->FindListObject(AliAODMCParticle::StdBranchName()));

	if (mcArray) {
	  // load MC header
	  if (fAdditionalChecks) ((TH1F*)(fOutput->FindObject("hZ4")))->Fill(zVertex);
	  mcHeader = (AliAODMCHeader*)aodEvent->GetList()->FindObject(AliAODMCHeader::StdBranchName());

	  if (mcHeader && fAdditionalChecks) ((TH1F*)(fOutput->FindObject("hZ5")))->Fill(zVertex);

	  // check on MC Lc Daughter
	  if (fAdditionalChecks) {
	    for (Int_t iii=0; iii<mcArray->GetEntries(); iii++)
	      SearchLcDaughter(mcArray,iii);
	  }

	}

      }

      if (fVtx1->GetNContributors()>0) {
	if (fAdditionalChecks) ((TH1F*)(fOutput->FindObject("hZ6")))->Fill(zVertex);

	TString firedTriggerClasses = aodEvent->GetFiredTriggerClasses(); // trigger class
	ULong64_t fTriggerMask=AliVEvent::kAnyINT;
	Bool_t check1 = kFALSE;
	if ( !fUseMCInfo && // don't do for MC...
	     (aodEvent->GetRunNumber()<136851 || aodEvent->GetRunNumber()>139517) ) { // ...and for PbPb 2010 data
	  if ( !(firedTriggerClasses.Contains("CINT1")) ) {
	    AliInfo(Form(" ======================== firedTriggerClasses.Data() = %s",firedTriggerClasses.Data()));
	    fCEvents->Fill(8);
	    if (fAdditionalChecks) ((TH1F*)(fOutput->FindObject("hZ8")))->Fill(zVertex);
	    check1 = kTRUE;
	  }
	}

	Bool_t isSelectedAAA = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & fTriggerMask);
	if (!isSelectedAAA) {
	  fCEvents->Fill(9);
	  if (fAdditionalChecks) ((TH1F*)(fOutput->FindObject("hZ9")))->Fill(zVertex);
	}

	if (!isSelectedAAA || check1) {
	  fCEvents->Fill(16);
	  if (fAdditionalChecks) ((TH1F*)(fOutput->FindObject("hZ16")))->Fill(zVertex);
	}

	fTriggerMask=AliVEvent::kAny;
	Bool_t isSelectedBBB = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & fTriggerMask);
	if (!isSelectedBBB) {
	  fCEvents->Fill(10);
	  if (fAdditionalChecks) ((TH1F*)(fOutput->FindObject("hZ10")))->Fill(zVertex);
	}

	if (titleVtx.Contains("Z")) {
	  fCEvents->Fill(11);
	  if (fAdditionalChecks) ((TH1F*)(fOutput->FindObject("hZ11")))->Fill(zVertex);
	}
	else if (titleVtx.Contains("3D")) {
	  fCEvents->Fill(12);
	  if (fAdditionalChecks) ((TH1F*)(fOutput->FindObject("hZ12")))->Fill(zVertex);
	} else {
	  fCEvents->Fill(13);
	  if (fAdditionalChecks) ((TH1F*)(fOutput->FindObject("hZ13")))->Fill(zVertex);
	}

	if (TMath::Abs(zVertex)<=fAnalCuts->GetMaxVtxZ()) {
	  fCEvents->Fill(14);
	  if (fAdditionalChecks) ((TH1F*)(fOutput->FindObject("hZ14")))->Fill(zVertex);
	}

	if ( fIsEventSelected ) {
	  if (fAdditionalChecks) ((TH1F*)(fOutput->FindObject("hZ7")))->Fill(zVertex);
	} else {
	  fCEvents->Fill(15);
	  if (fAdditionalChecks) ((TH1F*)(fOutput->FindObject("hZ15")))->Fill(zVertex);
	}

      } // nContributors>=1
    } // analysisArray exists
  } // magnetic field exists

  return;
}

//-----------------
void AliAnalysisTaskSELc2V0bachelor::CheckEventSelectionWithCandidates(AliAODEvent *aodEvent) {
  //
  // To fill control histograms
  //

  Float_t zVertex = fVtx1->GetZ();
  TString titleVtx=fVtx1->GetTitle();
  TString firedTriggerClasses = aodEvent->GetFiredTriggerClasses(); // trigger class
  ULong64_t fTriggerMask=AliVEvent::kAnyINT;

  ((TH1F*)(fOutput->FindObject("hEventsWithCandidates")))->Fill(6);
  if (fAdditionalChecks) ((TH1F*)(fOutput->FindObject("hZ6a")))->Fill(zVertex);

  Bool_t check1a = kFALSE;
  if ( !fUseMCInfo && // don't do for MC...
       (aodEvent->GetRunNumber()<136851 || aodEvent->GetRunNumber()>139517) ) { // ...and for PbPb 2010 data
    if ( !(firedTriggerClasses.Contains("CINT1")) ) {
      AliInfo(Form(" ======================== firedTriggerClasses.Data() = %s",firedTriggerClasses.Data()));
      ((TH1F*)(fOutput->FindObject("hEventsWithCandidates")))->Fill(8);
      if (fAdditionalChecks) ((TH1F*)(fOutput->FindObject("hZ8a")))->Fill(zVertex);
      check1a = kTRUE;
    }
  }

  Bool_t isSelectedAAAa = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & fTriggerMask);
  if (!isSelectedAAAa) {
    ((TH1F*)(fOutput->FindObject("hEventsWithCandidates")))->Fill(9);
    if (fAdditionalChecks) ((TH1F*)(fOutput->FindObject("hZ9a")))->Fill(zVertex);
  }

  if (!isSelectedAAAa || check1a) {
    ((TH1F*)(fOutput->FindObject("hEventsWithCandidates")))->Fill(16);
    if (fAdditionalChecks) ((TH1F*)(fOutput->FindObject("hZ16a")))->Fill(zVertex);
  }

  fTriggerMask=AliVEvent::kAny;
  Bool_t isSelectedBBBa = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & fTriggerMask);
  if (!isSelectedBBBa) {
    ((TH1F*)(fOutput->FindObject("hEventsWithCandidates")))->Fill(10);
    if (fAdditionalChecks) ((TH1F*)(fOutput->FindObject("hZ10a")))->Fill(zVertex);
  }

  if (titleVtx.Contains("Z")) {
    ((TH1F*)(fOutput->FindObject("hEventsWithCandidates")))->Fill(11);
    if (fAdditionalChecks) ((TH1F*)(fOutput->FindObject("hZ11a")))->Fill(zVertex);
  }
  else if (titleVtx.Contains("3D")) {
    ((TH1F*)(fOutput->FindObject("hEventsWithCandidates")))->Fill(12);
    if (fAdditionalChecks) ((TH1F*)(fOutput->FindObject("hZ12a")))->Fill(zVertex);
  } else {
    ((TH1F*)(fOutput->FindObject("hEventsWithCandidates")))->Fill(13);
    if (fAdditionalChecks) ((TH1F*)(fOutput->FindObject("hZ13a")))->Fill(zVertex);
  }

  if (TMath::Abs(zVertex)<=fAnalCuts->GetMaxVtxZ()) {
    ((TH1F*)(fOutput->FindObject("hEventsWithCandidates")))->Fill(14);
    if (fAdditionalChecks) ((TH1F*)(fOutput->FindObject("hZ14a")))->Fill(zVertex);
  }

  if ( fIsEventSelected ) {
    ((TH1F*)(fOutput->FindObject("hEventsWithCandidates")))->Fill(7);
    if (fAdditionalChecks) ((TH1F*)(fOutput->FindObject("hZ7a")))->Fill(zVertex);
  } else {
    ((TH1F*)(fOutput->FindObject("hEventsWithCandidates")))->Fill(15);
    if (fAdditionalChecks) ((TH1F*)(fOutput->FindObject("hZ15a")))->Fill(zVertex);
  }

  return;
}

//-----------------------------
Int_t AliAnalysisTaskSELc2V0bachelor::MatchToMC(AliAODRecoCascadeHF *lc2bacV0,
						Int_t *pdgDgLc2bacV0, Int_t *pdgDgV0,
						TClonesArray *mcArray) {
  //
  // This is now implemented in AliAODRecoCascadeHF
  //

  // bachelor
  AliAODTrack *bachelor = (AliAODTrack*)lc2bacV0->GetBachelor();
  if (!bachelor) return -1;
  Int_t labBachelor = TMath::Abs(bachelor->GetLabel());
  if (labBachelor<0) return -1;
  AliAODMCParticle *partBachelor = (AliAODMCParticle*)mcArray->At(labBachelor);
  if (!partBachelor) return -1;
  if (TMath::Abs(partBachelor->GetPdgCode())!=pdgDgLc2bacV0[0]) return -1;

  Int_t labBacMother = partBachelor->GetMother();
  if (labBacMother<0) return -1;
  AliAODMCParticle *partBacMother = (AliAODMCParticle*)mcArray->At(labBacMother);
  if (!partBacMother) return -1;
  if (TMath::Abs(partBacMother->GetPdgCode())!=4122) return -1;

  // V0
  AliAODTrack *posV0Daugh = (AliAODTrack*)lc2bacV0->Getv0PositiveTrack();
  AliAODTrack *negV0Daugh = (AliAODTrack*)lc2bacV0->Getv0NegativeTrack();
  if (!posV0Daugh || !negV0Daugh) return -1;

  Int_t labV0pos = TMath::Abs(posV0Daugh->GetLabel());
  Int_t labV0neg = TMath::Abs(negV0Daugh->GetLabel());
  if (labV0pos<0 || labV0neg<0) return -1;

  AliAODMCParticle *partV0pos = (AliAODMCParticle*)mcArray->At(labV0neg);
  AliAODMCParticle *partV0neg = (AliAODMCParticle*)mcArray->At(labV0pos);
  if (!partV0pos || !partV0neg) return -1;

  if ( ! ( (TMath::Abs(partV0pos->GetPdgCode())==pdgDgV0[0] &&
	    TMath::Abs(partV0neg->GetPdgCode())==pdgDgV0[1]) ||
	   (TMath::Abs(partV0pos->GetPdgCode())==pdgDgV0[1] &&
	    TMath::Abs(partV0neg->GetPdgCode())==pdgDgV0[0]) ) ) return -1;
  Int_t labV0posMother = partV0pos->GetMother();
  Int_t labV0negMother = partV0neg->GetMother();

  if (labV0posMother<0 || labV0negMother<0) return -1;
  if (labV0posMother!=labV0negMother) return -1;

  AliAODMCParticle *motherV0 = (AliAODMCParticle*)mcArray->At(labV0posMother);
  if (!motherV0) return-1;

  if (TMath::Abs(motherV0->GetPdgCode())!=pdgDgLc2bacV0[1]) return -1;
  Int_t labV0mother = motherV0->GetMother();
  if (labV0mother<0) return -1;
  AliAODMCParticle *gMotherV0 = (AliAODMCParticle*)mcArray->At(labV0mother);
  if (!gMotherV0) return-1;

  if ( !(pdgDgLc2bacV0[1]==310 && TMath::Abs(gMotherV0->GetPdgCode())==311) &&
       !(pdgDgLc2bacV0[1]==3122 && TMath::Abs(motherV0->GetPdgCode())==3122) ) return -1;

  if ( (pdgDgLc2bacV0[1]==310 && TMath::Abs(gMotherV0->GetPdgCode())==311) ) {
    Int_t labV0GMother = gMotherV0->GetMother();
    if (labV0GMother<0) return -1;
    AliAODMCParticle *ggMotherV0 = (AliAODMCParticle*)mcArray->At(labV0GMother);
    if (!ggMotherV0) return-1;

    if (TMath::Abs(ggMotherV0->GetPdgCode())!=4122) return -1;
    gMotherV0 = (AliAODMCParticle*)ggMotherV0;
    labV0mother=labV0GMother;
  }
  else if (pdgDgLc2bacV0[1]==3122 && TMath::Abs(motherV0->GetPdgCode())==3122) {
    if (TMath::Abs(gMotherV0->GetPdgCode())!=4122) return -1;
  }

  if (labBacMother!=labV0mother) {
    return -1;
  }

  return labBacMother;

}

//________________________________________________________________
Int_t AliAnalysisTaskSELc2V0bachelor::SearchLcDaughter(TClonesArray *arrayMC, Int_t iii) {
  //
  // This is to check Lc dinasty
  //

  Int_t indexToBeReturned=-999;

  Int_t pdgLc=4122;
  Int_t pdgLambda=3122;
  Int_t pdgV0=310;
  Int_t pdgK0=311;
  Int_t pdgBachelor=2212;
  Int_t pdgBachelorPi=211;

  TString fillthis="";
  fillthis="histMcStatLc";

  AliAODMCParticle *searchLc = dynamic_cast<AliAODMCParticle*>(arrayMC->At(iii));
  if(!searchLc) return -999;
  if (TMath::Abs(searchLc->GetPdgCode()) != pdgLc) return -999;

  ((TH1F*)(fOutput->FindObject(fillthis)))->Fill(0);
  indexToBeReturned = 0;

  ((TH1F*)(fOutput->FindObject(fillthis)))->Fill(TMath::Nint(searchLc->Charge()/3.)*1);
  indexToBeReturned = TMath::Nint(searchLc->Charge()/3.)*1;

  Int_t nDaughLc = searchLc->GetNDaughters();
  if (nDaughLc!=2) {
    ((TH1F*)(fOutput->FindObject(fillthis)))->Fill(TMath::Nint(searchLc->Charge()/3.)*10);
    indexToBeReturned = TMath::Nint(searchLc->Charge()/3.)*10;
    return indexToBeReturned;
  }

  Int_t index1=searchLc->GetDaughter(0);
  Int_t index2=searchLc->GetDaughter(1);
  if (index1<=0 || index2<=0) {
    return -999;
  }

  AliAODMCParticle *daugh1 = dynamic_cast<AliAODMCParticle*>(arrayMC->At(index1));
  AliAODMCParticle *daugh2 = dynamic_cast<AliAODMCParticle*>(arrayMC->At(index2));
  if (!daugh1 || !daugh2) return -999;

  Int_t daughPdg1 = TMath::Abs(daugh1->GetPdgCode());
  Int_t daughPdg2 = TMath::Abs(daugh2->GetPdgCode());
  if ( !( (daughPdg1==pdgBachelor && daughPdg2==pdgK0) ||
	  (daughPdg2==pdgBachelor && daughPdg1==pdgK0) ||
	  (daughPdg1==pdgLambda && daughPdg2==pdgBachelorPi) ||
	  (daughPdg2==pdgLambda && daughPdg1==pdgBachelorPi) ) ) {
    ((TH1F*)(fOutput->FindObject(fillthis)))->Fill(TMath::Nint(searchLc->Charge()/3.)*10);
    indexToBeReturned = TMath::Nint(searchLc->Charge()/3.)*10;
    return indexToBeReturned;
  }

  if (daughPdg1==pdgK0 || daughPdg1==pdgLambda) {
    index1=searchLc->GetDaughter(1);
    index2=searchLc->GetDaughter(0);
  }
  daugh1 = dynamic_cast<AliAODMCParticle*>(arrayMC->At(index1));
  daugh2 = dynamic_cast<AliAODMCParticle*>(arrayMC->At(index2));
  if (!daugh1 || !daugh2) return -999;

  daughPdg1=TMath::Abs(daugh1->GetPdgCode());
  daughPdg2=TMath::Abs(daugh2->GetPdgCode());

  if ( daughPdg1==pdgBachelor && daughPdg2==pdgK0 ) { // Lc+ -> p K0bar

    ((TH1F*)(fOutput->FindObject(fillthis)))->Fill(TMath::Nint(searchLc->Charge()/3.)*2);
    indexToBeReturned = TMath::Nint(searchLc->Charge()/3.)*2;

    Int_t nDaughK0 = daugh2->GetNDaughters();
    if (nDaughK0!=1) return -999;

    Int_t indexK0daugh=daugh2->GetDaughter(0);
    if (indexK0daugh<=0) return -999;

    AliAODMCParticle *daughK0 = dynamic_cast<AliAODMCParticle*>(arrayMC->At(indexK0daugh));
    if (!daughK0) return -999;

    Int_t daughK0Pdg=TMath::Abs(daughK0->GetPdgCode());
    if (daughK0Pdg!=pdgV0) {
      ((TH1F*)(fOutput->FindObject(fillthis)))->Fill(TMath::Nint(searchLc->Charge()/3.)*4); // K0L
      indexToBeReturned = TMath::Nint(searchLc->Charge()/3.)*4;
      return indexToBeReturned;
    }
    ((TH1F*)(fOutput->FindObject(fillthis)))->Fill(TMath::Nint(searchLc->Charge()/3.)*3); // K0S
    indexToBeReturned = TMath::Nint(searchLc->Charge()/3.)*3;

    Int_t nDaughK0S = daughK0->GetNDaughters();
    if (nDaughK0S!=2) {
      ((TH1F*)(fOutput->FindObject(fillthis)))->Fill(TMath::Nint(searchLc->Charge()/3.)*5); // other decays for K0S
      indexToBeReturned = TMath::Nint(searchLc->Charge()/3.)*5;
      return indexToBeReturned;
    }

    index1=daughK0->GetDaughter(0);
    index2=daughK0->GetDaughter(1);
    if(index1<=0 || index2<=0) {
      return -999;
    }

    AliAODMCParticle *daughK0S1 = dynamic_cast<AliAODMCParticle*>(arrayMC->At(index1));
    AliAODMCParticle *daughK0S2 = dynamic_cast<AliAODMCParticle*>(arrayMC->At(index2));
    if (!daughK0S1 || !daughK0S2) return -999;

    Int_t daughK0S1pdg=TMath::Abs(daughK0S1->GetPdgCode());
    Int_t daughK0S2pdg=TMath::Abs(daughK0S2->GetPdgCode());

    if ( daughK0S1pdg==211 && daughK0S2pdg==211 ) {
      ((TH1F*)(fOutput->FindObject(fillthis)))->Fill(TMath::Nint(searchLc->Charge()/3.)*6); // K0S -> pi+ pi-
      indexToBeReturned = TMath::Nint(searchLc->Charge()/3.)*6;
    } else {
      ((TH1F*)(fOutput->FindObject(fillthis)))->Fill(TMath::Nint(searchLc->Charge()/3.)*5); // other decays for K0S
      indexToBeReturned = TMath::Nint(searchLc->Charge()/3.)*5;
    }

  } //if (daughPdg1==pdgBachelor && daughPdg2==pdgK0)
  else if ( daughPdg1==pdgBachelorPi && daughPdg2==pdgLambda ) { // Lc+ -> pi+ Lambda

    ((TH1F*)(fOutput->FindObject(fillthis)))->Fill(TMath::Nint(searchLc->Charge()/3.)*7);
    indexToBeReturned = TMath::Nint(searchLc->Charge()/3.)*7;

    Int_t nDaughL = daugh2->GetNDaughters();
    if (nDaughL!=2) {
      ((TH1F*)(fOutput->FindObject(fillthis)))->Fill(TMath::Nint(searchLc->Charge()/3.)*8);
      indexToBeReturned = TMath::Nint(searchLc->Charge()/3.)*8;
      return indexToBeReturned;
    }

    index1=daugh2->GetDaughter(0);
    index2=daugh2->GetDaughter(1);
    if(index1<=0 || index2<=0) {
      return -999;
    }

    AliAODMCParticle *daughL1 = dynamic_cast<AliAODMCParticle*>(arrayMC->At(index1));
    AliAODMCParticle *daughL2 = dynamic_cast<AliAODMCParticle*>(arrayMC->At(index2));
    if (!daughL1 || !daughL2) return -999;

    Int_t daughL1pdg=TMath::Abs(daughL1->GetPdgCode());
    Int_t daughL2pdg=TMath::Abs(daughL2->GetPdgCode());
    if ( (daughL1pdg==211 && daughL2pdg==2212) ||
	 (daughL2pdg==211 && daughL1pdg==2212) ) {
      ((TH1F*)(fOutput->FindObject(fillthis)))->Fill(TMath::Nint(searchLc->Charge()/3.)*9);
      indexToBeReturned = TMath::Nint(searchLc->Charge()/3.)*9;
    } else {
      ((TH1F*)(fOutput->FindObject(fillthis)))->Fill(TMath::Nint(searchLc->Charge()/3.)*8);
      indexToBeReturned = TMath::Nint(searchLc->Charge()/3.)*8;
    }

  } //else if (daughPdg1==pdgBachelorPi && daughPdg2==pdgLambda)

  return indexToBeReturned;
}

//________________________________________________________________
void AliAnalysisTaskSELc2V0bachelor::FillArmPodDistribution(AliAODv0 *vZero,
							    TString histoTitle,
							    TList *histoList) {
  //
  // This is to fill Armenteros Podolanski plots
  //

  Double_t alpha = vZero->AlphaV0();
  Double_t qT    = vZero->PtArmV0();

  ((TH2F*)(histoList->FindObject(histoTitle)))->Fill(alpha,qT);

}

//-------------------------------------------------------------------------------
void AliAnalysisTaskSELc2V0bachelor::CheckCandidatesAtDifferentLevels(AliAODRecoCascadeHF *part, AliRDHFCutsLctoV0* cutsAnal) {
  //
  // This is to check candidates at different levels
  //

  Bool_t areCutsUsingPID = cutsAnal->GetIsUsePID();

  AliAODv0 * v0part = (AliAODv0*)part->Getv0();
  Bool_t onFlyV0 = v0part->GetOnFlyStatus(); // on-the-flight V0s

  AliAODTrack *bachelor = (AliAODTrack*)part->GetBachelor();

  if ( !onFlyV0 )
    ((TH1F*)(fOutput->FindObject("hCandidateSelection")))->Fill(3); // it counts number of candidates coming from offline V0s

  if ( cutsAnal->IsInFiducialAcceptance(part->Pt(),part->Y(4122)) )
    ((TH1F*)(fOutput->FindObject("hCandidateSelection")))->Fill(4);
  if ( (((cutsAnal->IsSelected(part,AliRDHFCuts::kTracks))&(AliRDHFCutsLctoV0::kLcToK0Spr))==(AliRDHFCutsLctoV0::kLcToK0Spr)) )
    ((TH1F*)(fOutput->FindObject("hCandidateSelection")))->Fill(5);
  cutsAnal->SetUsePID(kFALSE);
  if ( (((cutsAnal->IsSelected(part,AliRDHFCuts::kCandidate))&(AliRDHFCutsLctoV0::kLcToK0Spr))==(AliRDHFCutsLctoV0::kLcToK0Spr)) )
    ((TH1F*)(fOutput->FindObject("hCandidateSelection")))->Fill(6);
  cutsAnal->SetUsePID(areCutsUsingPID);
  if ( (((cutsAnal->IsSelected(part,AliRDHFCuts::kPID))&(AliRDHFCutsLctoV0::kLcToK0Spr))==(AliRDHFCutsLctoV0::kLcToK0Spr)) )
    ((TH1F*)(fOutput->FindObject("hCandidateSelection")))->Fill(7);
  if ( (((cutsAnal->IsSelected(part,AliRDHFCuts::kCandidate))&(AliRDHFCutsLctoV0::kLcToK0Spr))==(AliRDHFCutsLctoV0::kLcToK0Spr)) )
    ((TH1F*)(fOutput->FindObject("hCandidateSelection")))->Fill(8);
  if ( (((cutsAnal->IsSelected(part,AliRDHFCuts::kAll))&(AliRDHFCutsLctoV0::kLcToK0Spr))==(AliRDHFCutsLctoV0::kLcToK0Spr)) )
    ((TH1F*)(fOutput->FindObject("hCandidateSelection")))->Fill(9);

  if ( cutsAnal->IsInFiducialAcceptance(part->Pt(),part->Y(4122)) ) {

    if ( ( (cutsAnal->IsSelected(part,AliRDHFCuts::kTracks))&(AliRDHFCutsLctoV0::kLcToK0Spr)) == (AliRDHFCutsLctoV0::kLcToK0Spr) ) {

      Int_t aaa = cutsAnal->IsSelected(part,AliRDHFCuts::kTracks);
      if ( (aaa&AliRDHFCutsLctoV0::kLcToK0Spr)==AliRDHFCutsLctoV0::kLcToK0Spr ) {
	if ( ( (aaa&AliRDHFCutsLctoV0::kLcToLpi)==AliRDHFCutsLctoV0::kLcToLpi && bachelor->Charge()==-1)  ||
	     ( (aaa&AliRDHFCutsLctoV0::kLcToLBarpi)==AliRDHFCutsLctoV0::kLcToLBarpi && bachelor->Charge()==+1) )
	  ((TH1F*)(fOutput->FindObject("hSwitchOnCandidates1")))->Fill( -aaa );
	else
	  ((TH1F*)(fOutput->FindObject("hSwitchOnCandidates1")))->Fill( aaa );
      }

      cutsAnal->SetUsePID(kFALSE);
      aaa = cutsAnal->IsSelected(part,AliRDHFCuts::kCandidate);
      if ((aaa&AliRDHFCutsLctoV0::kLcToK0Spr)==AliRDHFCutsLctoV0::kLcToK0Spr) {
	if ( ( (aaa&AliRDHFCutsLctoV0::kLcToLpi)==AliRDHFCutsLctoV0::kLcToLpi && bachelor->Charge()==-1) ||
	     ( (aaa&AliRDHFCutsLctoV0::kLcToLBarpi)==AliRDHFCutsLctoV0::kLcToLBarpi && bachelor->Charge()==+1) )
	  ((TH1F*)(fOutput->FindObject("hSwitchOnCandidates2")))->Fill( -aaa );
	else
	  ((TH1F*)(fOutput->FindObject("hSwitchOnCandidates2")))->Fill( aaa );
      }
      cutsAnal->SetUsePID(areCutsUsingPID);

      aaa = cutsAnal->IsSelected(part,AliRDHFCuts::kPID);
      if ((aaa&AliRDHFCutsLctoV0::kLcToK0Spr)==AliRDHFCutsLctoV0::kLcToK0Spr) {
	if ( ( (aaa&AliRDHFCutsLctoV0::kLcToLpi)==AliRDHFCutsLctoV0::kLcToLpi && bachelor->Charge()==-1) ||
	     ( (aaa&AliRDHFCutsLctoV0::kLcToLBarpi)==AliRDHFCutsLctoV0::kLcToLBarpi && bachelor->Charge()==+1) )
	  ((TH1F*)(fOutput->FindObject("hSwitchOnCandidates3")))->Fill( -aaa );
	else
	  ((TH1F*)(fOutput->FindObject("hSwitchOnCandidates3")))->Fill( aaa );
      }

      aaa = cutsAnal->IsSelected(part,AliRDHFCuts::kAll);
      if ((aaa&AliRDHFCutsLctoV0::kLcToK0Spr)==AliRDHFCutsLctoV0::kLcToK0Spr) {
	if ( ( (aaa&AliRDHFCutsLctoV0::kLcToLpi)==AliRDHFCutsLctoV0::kLcToLpi && bachelor->Charge()==-1) ||
	     ( (aaa&AliRDHFCutsLctoV0::kLcToLBarpi)==AliRDHFCutsLctoV0::kLcToLBarpi && bachelor->Charge()==+1) )
	  ((TH1F*)(fOutput->FindObject("hSwitchOnCandidates4")))->Fill( -aaa );
	else
	  ((TH1F*)(fOutput->FindObject("hSwitchOnCandidates4")))->Fill( aaa );
      }

    }
  }

  return;
}

//-------------------------------------------------------------------------------
void AliAnalysisTaskSELc2V0bachelor::FillTheTree(AliAODRecoCascadeHF *part, AliRDHFCutsLctoV0 *cutsAnal, TClonesArray *mcArray, Int_t isLc) {
  //
  // This is to fill tree
  //

  Double_t mk0sPDG = TDatabasePDG::Instance()->GetParticle(310)->Mass();
  Double_t mLPDG   = TDatabasePDG::Instance()->GetParticle(3122)->Mass();

  Double_t invmassLc = part->InvMassLctoK0sP();
  Double_t invmassLc2Lpi = part->InvMassLctoLambdaPi();

  AliAODTrack *bachelor = (AliAODTrack*)part->GetBachelor();

  AliAODv0 * v0part = (AliAODv0*)part->Getv0();
  Double_t dcaV0ptp = v0part->GetDCA();
  Double_t invmassK0S = v0part->MassK0Short();
  Double_t invmassLambda = v0part->MassLambda();
  Double_t invmassLambdaBar = v0part->MassAntiLambda();

  Int_t isLc2LBarpi=0, isLc2Lpi=0;
  Int_t mcLabel = -1;
  Int_t isDp2K0Spi=0, isDs2K0SK=0;
  Int_t mcLabel2 = -1;
  Int_t mcLabel3 = -1;
  if (fUseMCInfo) {
    Int_t pdgCand = 4122;
    Int_t pdgDgLctoV0bachelor[2]={211,3122};
    Int_t pdgDgV0toDaughters[2]={2212,211};
    mcLabel = part->MatchToMC(pdgCand,pdgDgLctoV0bachelor[1],pdgDgLctoV0bachelor,pdgDgV0toDaughters,mcArray,kTRUE);
    if (mcLabel!=-1) {
      if (bachelor->Charge()==-1) isLc2LBarpi=1;
      if (bachelor->Charge()==+1) isLc2Lpi=1;
    }

    Int_t pdgCand2 = 411; // D+ -> pi+ K0S
    Int_t pdgCand3 = 431; // Ds+ -> K+ K0S
    Int_t pdgDgCand2[2]={211,310};
    Int_t pdgDgCand3[2]={321,310};
    pdgDgV0toDaughters[0]=211;
    pdgDgV0toDaughters[1]=211;
    mcLabel2 = part->MatchToMC(pdgCand2,pdgDgCand2[1],pdgDgCand2,pdgDgV0toDaughters,mcArray,kTRUE);
    mcLabel3 = part->MatchToMC(pdgCand3,pdgDgCand3[1],pdgDgCand3,pdgDgV0toDaughters,mcArray,kTRUE);
    if (mcLabel2!=-1) isDp2K0Spi=1;
    if (mcLabel3!=-1) isDs2K0SK=1;
  }

  Int_t isLcByMC = isLc+isLc2LBarpi*2+isLc2Lpi*4+isDp2K0Spi*8+isDs2K0SK*16;

  Int_t isK0S = 0;
  Int_t isLambda = 0;
  Int_t isLambdaBar = 0;
  Int_t isGamma = 0;
  if (fUseMCInfo) {
    Int_t pdgDg2prong[2] = {211, 211};
    Int_t labelK0S = v0part->MatchToMC(310,mcArray,2,pdgDg2prong);
    if (labelK0S>=0) isK0S = 1;

    pdgDg2prong[0] = 211;
    pdgDg2prong[1] = 2212;
    Int_t lambdaLabel = v0part->MatchToMC(3122,mcArray,2,pdgDg2prong);
    if (lambdaLabel>=0) {
      AliAODMCParticle *lambdaTrack = (AliAODMCParticle*)mcArray->At(lambdaLabel);
      if (lambdaTrack->GetPdgCode()==3122) isLambda = 1;
      else if (lambdaTrack->GetPdgCode()==-3122) isLambdaBar = 1;
    }

    pdgDg2prong[0] = 11;
    pdgDg2prong[1] = 11;
    Int_t gammaLabel = v0part->MatchToMC(22,mcArray,2,pdgDg2prong);
    if (gammaLabel>=0) {
      AliAODMCParticle *gammaTrack = (AliAODMCParticle*)mcArray->At(gammaLabel);
      if (gammaTrack->GetPdgCode()==22) isGamma = 1;
    }
  }

  Int_t isV0ByMC = isK0S+isLambdaBar*2+isLambda*4+isGamma*8;

  Int_t isBachelorSelected = (bachelor->TestFilterMask(BIT(4)))*1 + (!(bachelor->TestFilterMask(BIT(4))))*2;
  isBachelorSelected += (bachelor->GetLabel()<0)*4 + (bachelor->GetLabel()>=0)*8;
  if ( ( !(bachelor->HasPointOnITSLayer(0)) && !(bachelor->HasPointOnITSLayer(1)) ) )
    isBachelorSelected += 16;
  else {
    if ( bachelor->HasPointOnITSLayer(0) && !(bachelor->HasPointOnITSLayer(1)) )
      isBachelorSelected += 32;
    else if ( !(bachelor->HasPointOnITSLayer(0)) && bachelor->HasPointOnITSLayer(1) )
      isBachelorSelected += 64;
    else
      isBachelorSelected += 128;
  }

  AliAODTrack *v0pos = (AliAODTrack*)part->Getv0PositiveTrack();
  AliAODTrack *v0neg = (AliAODTrack*)part->Getv0NegativeTrack(); 

  Int_t areV0daughtersSelected = (v0pos->TestFilterMask(BIT(4)))*1 + (!(v0pos->TestFilterMask(BIT(4))))*2;
  areV0daughtersSelected += (v0pos->GetLabel()<0)*4 + (v0pos->GetLabel()>=0)*8;
  areV0daughtersSelected += (v0neg->TestFilterMask(BIT(4)))*16 + (!(v0neg->TestFilterMask(BIT(4))))*32;
  areV0daughtersSelected += (v0neg->GetLabel()<0)*64 + (v0neg->GetLabel()>=0)*128;

  Double_t nSigmaITSpr=-999.;
  cutsAnal->GetPidHF()->GetnSigmaITS(bachelor,4,nSigmaITSpr);
  Double_t nSigmaTPCpr=-999.;
  cutsAnal->GetPidHF()->GetnSigmaTPC(bachelor,4,nSigmaTPCpr);
  Double_t nSigmaTOFpr=-999.;
  cutsAnal->GetPidHF()->GetnSigmaTOF(bachelor,4,nSigmaTOFpr);

  Double_t nSigmaITSpi=-999.;
  cutsAnal->GetPidHF()->GetnSigmaITS(bachelor,2,nSigmaITSpi);
  Double_t nSigmaTPCpi=-999.;
  cutsAnal->GetPidHF()->GetnSigmaTPC(bachelor,2,nSigmaTPCpi);
  Double_t nSigmaTOFpi=-999.;
  cutsAnal->GetPidHF()->GetnSigmaTOF(bachelor,2,nSigmaTOFpi);

  Double_t nSigmaITSka=-999.;
  cutsAnal->GetPidHF()->GetnSigmaITS(bachelor,3,nSigmaITSka);
  Double_t nSigmaTPCka=-999.;
  cutsAnal->GetPidHF()->GetnSigmaTPC(bachelor,3,nSigmaTPCka);
  Double_t nSigmaTOFka=-999.;
  cutsAnal->GetPidHF()->GetnSigmaTOF(bachelor,3,nSigmaTOFka);


  Int_t flagToCheckCandidate = 1*(TMath::Abs(invmassK0S-mk0sPDG)<=0.050);
  flagToCheckCandidate+=2*((TMath::Abs(invmassLambdaBar-mLPDG)<=0.050) && (bachelor->Charge()==-1));
  flagToCheckCandidate+=4*((TMath::Abs(invmassLambda-mLPDG)<=0.050) && (bachelor->Charge()==+1));
  flagToCheckCandidate+=8*((TMath::Abs(invmassLambdaBar-mLPDG)<=0.050) && (bachelor->Charge()==+1));
  flagToCheckCandidate+=16*((TMath::Abs(invmassLambda-mLPDG)<=0.050) && (bachelor->Charge()==-1));

  /*
  Bool_t areCutsUsingPID = cutsAnal->GetIsUsePID();
  cutsAnal->SetUsePID(kFALSE);
  Int_t aaa = cutsAnal->IsSelected(part,AliRDHFCuts::kCandidate);
  if ( (aaa&AliRDHFCutsLctoV0::kLcToK0Spr)==AliRDHFCutsLctoV0::kLcToK0Spr ) {
    if ( aaa==AliRDHFCutsLctoV0::kLcToK0Spr ) {
      flagToCheckCandidate = aaa; // Lc->K0S+p OK
    } else {
      if ( (aaa&AliRDHFCutsLctoV0::kLcToLpi)==AliRDHFCutsLctoV0::kLcToLpi ) {
	if (bachelor->Charge()==+1)
	  flagToCheckCandidate = aaa; // Lc->Lambda+pi+
	else if (bachelor->Charge()==-1)
	  flagToCheckCandidate =-aaa; // Lambda+pi- AS Lc->K0S+p candidate
      }
      if ( (aaa&AliRDHFCutsLctoV0::kLcToLBarpi)==AliRDHFCutsLctoV0::kLcToLBarpi ) {
	if (bachelor->Charge()==-1)
	  flagToCheckCandidate = aaa; // Lc->LambdaBar+pi-
	else if (bachelor->Charge()==+1)
	  flagToCheckCandidate =-aaa; // LambdaBar+pi+ AS Lc->K0S+p candidate
      }
    }
  } else {
    //if ( aaa==AliRDHFCutsLctoV0::kLcToK0Spr ) {
    //flagToCheckCandidate = -10-(AliRDHFCutsLctoV0::kLcToK0Spr); // NEVER
    //} else {
      if ( (aaa&AliRDHFCutsLctoV0::kLcToLpi)==AliRDHFCutsLctoV0::kLcToLpi ) {
	if (bachelor->Charge()==+1)
	  flagToCheckCandidate = aaa; // Lc->Lambda+pi+ OK
	else if (bachelor->Charge()==-1)
	  flagToCheckCandidate =-aaa; // Lambda+pi- AS Lc->Lambda+pi+ candidate
      }
      if ( (aaa&AliRDHFCutsLctoV0::kLcToLBarpi)==AliRDHFCutsLctoV0::kLcToLBarpi ) {
	if (bachelor->Charge()==-1)
	  flagToCheckCandidate = aaa; // Lc->LambdaBar+pi- OK
	else if (bachelor->Charge()==+1)
	  flagToCheckCandidate =-aaa; // LambdaBar+pi+ AS Lc->LambdaBar+pi- candidate
      }
      //}
  }
  cutsAnal->SetUsePID(areCutsUsingPID);
  */

  fCandidateVariables[ 0] = fUseMCInfo+isLcByMC; // 0: real data; 1: bkg; 2: Lc->K0S+p; 3: Lc->LambdaBar+pbar; 5: Lc->Lambda+p
  fCandidateVariables[ 1] = fUseMCInfo+isV0ByMC; // 0: real data; 1: bkg; 2: K0S->pi+pi; 3: LambdaBar->pbar+pi+; 5: Lambda->p+pi-
  fCandidateVariables[ 2] = isBachelorSelected;
  fCandidateVariables[ 3] = areV0daughtersSelected;
  fCandidateVariables[ 4] = flagToCheckCandidate;
  fCandidateVariables[ 5] = invmassLc;
  fCandidateVariables[ 6] = invmassLc2Lpi;
  fCandidateVariables[ 7] = part->InvMass2Prongs(0,1,211,310); // D+ -> pi+ K0S
  fCandidateVariables[ 8] = part->InvMass2Prongs(0,1,321,310); // D+S -> K+ K0S
  fCandidateVariables[ 9] = invmassK0S;
  fCandidateVariables[10] = invmassLambda;
  fCandidateVariables[11] = invmassLambdaBar;
  fCandidateVariables[12] = v0part->InvMass2Prongs(0,1,11,11);
  fCandidateVariables[13] = part->GetDCA();
  fCandidateVariables[14] = dcaV0ptp;
  fCandidateVariables[15] = part->Getd0Prong(0);
  fCandidateVariables[16] = part->Getd0Prong(1);
  fCandidateVariables[17] = v0part->Getd0Prong(0);
  fCandidateVariables[18] = v0part->Getd0Prong(1);
  fCandidateVariables[19] = part->CosPointingAngle();
  fCandidateVariables[20] = part->CosV0PointingAngle();
  fCandidateVariables[21] = v0part->RadiusSecVtx();
  fCandidateVariables[22] = nSigmaITSpr;
  fCandidateVariables[23] = nSigmaITSpi;
  fCandidateVariables[24] = nSigmaITSka;
  fCandidateVariables[25] = nSigmaTPCpr;
  fCandidateVariables[26] = nSigmaTPCpi;
  fCandidateVariables[27] = nSigmaTPCka;
  fCandidateVariables[28] = nSigmaTOFpr;
  fCandidateVariables[29] = nSigmaTOFpi;
  fCandidateVariables[30] = nSigmaTOFka;
  fCandidateVariables[31] = part->Y(4122);
  fCandidateVariables[32] = bachelor->Eta();
  fCandidateVariables[33] = v0pos->Eta();
  fCandidateVariables[34] = v0neg->Eta();
  fCandidateVariables[35] = part->P();
  fCandidateVariables[36] = part->Pt();
  fCandidateVariables[37] = v0part->P();
  fCandidateVariables[38] = v0part->Pt();
  fCandidateVariables[39] = bachelor->P();
  fCandidateVariables[40] = bachelor->Pt();
  fCandidateVariables[41] = v0pos->P();
  fCandidateVariables[42] = v0pos->Pt();
  fCandidateVariables[43] = v0neg->P();
  fCandidateVariables[44] = v0neg->Pt();
  fCandidateVariables[45] = part->DecayLength();
  fCandidateVariables[46] = part->DecayLengthV0();
  fCandidateVariables[47] = part->CosPointingAngleXY();
  fCandidateVariables[48] = part->CosV0PointingAngleXY();
  fCandidateVariables[49] = part->DecayLengthXY();
  fCandidateVariables[50] = part->DecayLengthXYV0();
  fCandidateVariables[51] = part->NormalizedDecayLength();
  fCandidateVariables[52] = part->NormalizedV0DecayLength();
  fCandidateVariables[53] = part->NormalizedDecayLengthXY();
  fCandidateVariables[54] = part->NormalizedV0DecayLengthXY();
  Double_t xVtxLc=0, yVtxLc=0, zVtxLc=0;
  Double_t xLcMC=0,yLcMC=0,zLcMC=0;
  Double_t pxVtxBachelor=0, pyVtxBachelor=0, pzVtxBachelor=0;
  Double_t dcaForLc = PropagateToDCA(v0part,bachelor,fBzkG, xVtxLc, yVtxLc, zVtxLc, pxVtxBachelor, pyVtxBachelor, pzVtxBachelor);
  if (isLc) {
    Int_t pdgCand0 = 4122;
    Int_t pdgDgLctoV0bachelor0[2]={2212,310};
    Int_t pdgDgV0toDaughters0[2]={211,211};
    Int_t mcLabel0 = part->MatchToMC(pdgCand0,pdgDgLctoV0bachelor0[1],pdgDgLctoV0bachelor0,pdgDgV0toDaughters0,mcArray,kTRUE);
    AliAODMCParticle *partLc = dynamic_cast<AliAODMCParticle*>(mcArray->At(mcLabel0));
    if(partLc){
      AliAODMCParticle *partLcDaug0 = dynamic_cast<AliAODMCParticle*>(mcArray->At(partLc->GetDaughter(0)));
      if(partLcDaug0){
	xLcMC=partLcDaug0->Xv(), yLcMC=partLcDaug0->Yv(), zLcMC=partLcDaug0->Zv();
      }
    }
  } else if (isLc2LBarpi || isLc2Lpi) {
    AliAODMCParticle *partLc = dynamic_cast<AliAODMCParticle*>(mcArray->At(mcLabel));
    AliAODMCParticle *partLcDaug0 = dynamic_cast<AliAODMCParticle*>(mcArray->At(partLc->GetDaughter(0)));
    xLcMC=partLcDaug0->Xv(), yLcMC=partLcDaug0->Yv(), zLcMC=partLcDaug0->Zv();
  } else if (isDp2K0Spi) {
    AliAODMCParticle *partLc = dynamic_cast<AliAODMCParticle*>(mcArray->At(mcLabel2));
    AliAODMCParticle *partLcDaug0 = dynamic_cast<AliAODMCParticle*>(mcArray->At(partLc->GetDaughter(0)));
    xLcMC=partLcDaug0->Xv(), yLcMC=partLcDaug0->Yv(), zLcMC=partLcDaug0->Zv();
  } else if (isDs2K0SK) {
    AliAODMCParticle *partLc = dynamic_cast<AliAODMCParticle*>(mcArray->At(mcLabel3));
    AliAODMCParticle *partLcDaug0 = dynamic_cast<AliAODMCParticle*>(mcArray->At(partLc->GetDaughter(0)));
    xLcMC=partLcDaug0->Xv(), yLcMC=partLcDaug0->Yv(), zLcMC=partLcDaug0->Zv();
  }
  fCandidateVariables[55]=dcaForLc;
  fCandidateVariables[56]=part->GetSecVtxX();
  fCandidateVariables[57]=part->GetSecVtxY();
  fCandidateVariables[58]=part->GetSecVtxZ();
  fCandidateVariables[59]=xVtxLc;
  fCandidateVariables[60]=yVtxLc;
  fCandidateVariables[61]=zVtxLc;
  fCandidateVariables[62]=xLcMC;
  fCandidateVariables[63]=yLcMC;
  fCandidateVariables[64]=zLcMC;
  fCandidateVariables[65]=bachelor->Px();
  fCandidateVariables[66]=bachelor->Py();
  fCandidateVariables[67]=pxVtxBachelor;
  fCandidateVariables[68]=pyVtxBachelor;
  fCandidateVariables[69]=v0part->Px();
  fCandidateVariables[70]=v0part->Py();
  fCandidateVariables[71]=fVtx1->GetX();
  fCandidateVariables[72]=fVtx1->GetY();
  fCandidateVariables[73]=fVtx1->GetZ();
  fCandidateVariables[74]=part->CosThetaStar(0,4122,2212,310);
  fCandidateVariables[75]=part->CosThetaStar(1,4122,2212,310);
  fCandidateVariables[76]=v0part->Eta();
  fCandidateVariables[77]=v0part->Y(310);
  fCandidateVariables[78]=pzVtxBachelor;
  fCandidateVariables[79]=v0part->Pz();

  //fCandidateVariables[65] = bachelor->Px();
  //fCandidateVariables[66] = bachelor->Py();
  //fCandidateVariables[67] = bachelor->Pz();
  //fCandidateVariables[68] = v0pos->Px();
  //fCandidateVariables[69] = v0pos->Py();
  //fCandidateVariables[70] = v0pos->Pz();
  //fCandidateVariables[71] = v0neg->Px();
  //fCandidateVariables[72] = v0neg->Py();
  //fCandidateVariables[73] = v0neg->Pz();
  //fCandidateVariables[74] = part->PxProng(0);
  //fCandidateVariables[75] = part->PyProng(0);
  //fCandidateVariables[76] = part->PzProng(0);
  //fCandidateVariables[77] = part->PxProng(1);
  //fCandidateVariables[78] = part->PyProng(1);
  //fCandidateVariables[79] = part->PzProng(1);
  //fCandidateVariables[80] = v0part->PxProng(0);
  //fCandidateVariables[81] = v0part->PyProng(0);
  //fCandidateVariables[82] = v0part->PzProng(0);
  //fCandidateVariables[83] = v0part->PxProng(1);
  //fCandidateVariables[84] = v0part->PyProng(1);
  //fCandidateVariables[85] = v0part->PzProng(1);
  //fCandidateVariables[86] = part->QtProng(0);
  //fCandidateVariables[87] = part->Alpha();

  fVariablesTree->Fill();

  return;
}

//-------------------------------------------------------------------------------
void AliAnalysisTaskSELc2V0bachelor::DefineTreeVariables() {
  //
  // This is to define tree variables
  //

  const char* nameoutput = GetOutputSlot(4)->GetContainer()->GetName();
  fVariablesTree = new TTree(nameoutput,"Candidates variables tree");
  Int_t nVar = 80;
  fCandidateVariables = new Float_t [nVar];
  TString * fCandidateVariableNames = new TString[nVar];
  fCandidateVariableNames[ 0]="isLcByMC";
  fCandidateVariableNames[ 1]="isV0ByMC";
  fCandidateVariableNames[ 2]="flagToCheckBachelor";
  fCandidateVariableNames[ 3]="flagToCheckV0daughters";
  fCandidateVariableNames[ 4]="flagToCheckCandidate";
  fCandidateVariableNames[ 5]="massLc2K0Sp";
  fCandidateVariableNames[ 6]="massLc2Lambdapi";
  fCandidateVariableNames[ 7]="massD2K0Spi"; // D+ -> pi+ K0S
  fCandidateVariableNames[ 8]="massDS2K0SK"; // D+S -> K+ K0S
  fCandidateVariableNames[ 9]="massK0S";
  fCandidateVariableNames[10]="massLambda";
  fCandidateVariableNames[11]="massLambdaBar";
  fCandidateVariableNames[12]="massGamma";
  fCandidateVariableNames[13]="dcaLcptp"; // DCA Lc prong-to-prong
  fCandidateVariableNames[14]="dcaV0ptp";
  fCandidateVariableNames[15]="tImpParBach";
  fCandidateVariableNames[16]="tImpParV0";
  fCandidateVariableNames[17]="dcaV0postoPV";
  fCandidateVariableNames[18]="dcaV0negtoPV";
  fCandidateVariableNames[19]="cosPALc";
  fCandidateVariableNames[20]="cosPAK0S";
  fCandidateVariableNames[21]="rhoV0";
  fCandidateVariableNames[22]="nSigmaITSpr";
  fCandidateVariableNames[23]="nSigmaITSpi";
  fCandidateVariableNames[24]="nSigmaITSka";
  fCandidateVariableNames[25]="nSigmaTPCpr";
  fCandidateVariableNames[26]="nSigmaTPCpi";
  fCandidateVariableNames[27]="nSigmaTPCka";
  fCandidateVariableNames[28]="nSigmaTOFpr";
  fCandidateVariableNames[29]="nSigmaTOFpi";
  fCandidateVariableNames[30]="nSigmaTOFka";
  fCandidateVariableNames[31]="yLc";
  fCandidateVariableNames[32]="etaBach"; // etaBachelor
  fCandidateVariableNames[33]="etaV0pos"; // etaV0pos
  fCandidateVariableNames[34]="etaV0neg"; // etaV0neg
  fCandidateVariableNames[35]="LcP"; // @ DCA
  fCandidateVariableNames[36]="LcPt"; // @ DCA
  fCandidateVariableNames[37]="v0P"; // @ V0 DCA
  fCandidateVariableNames[38]="v0Pt"; // @ V0 DCA
  fCandidateVariableNames[39]="bachelorP"; // @ prim vtx
  fCandidateVariableNames[40]="bachelorPt"; // @ prim vtx
  fCandidateVariableNames[41]="V0positiveP"; // @ prim vtx
  fCandidateVariableNames[42]="V0positivePt"; // @ prim vtx
  fCandidateVariableNames[43]="V0negativeP"; // @ prim vtx
  fCandidateVariableNames[44]="V0negativePt"; // @ prim vtx
  fCandidateVariableNames[45]="decayLengthLc";
  fCandidateVariableNames[46]="decayLengthV0";
  fCandidateVariableNames[47]="cosPALcXY"; // cosPA XY x Lc
  fCandidateVariableNames[48]="cosPAV0XY"; // cosPA XY x V0
  fCandidateVariableNames[49]="decayLengthLcXY"; // decay length XY x Lc
  fCandidateVariableNames[50]="decayLengthV0XY"; // decay length XY x V0
  fCandidateVariableNames[51]="normalizedDecayLengthLc"; // normalized decay length x Lc
  fCandidateVariableNames[52]="normalizedDecayLengthV0"; // normalized decay length x V0
  fCandidateVariableNames[53]="normalizedDecayLengthXYLc"; // normalized decay length XY x Lc
  fCandidateVariableNames[54]="normalizedDecayLengthXYV0"; // normalized decay length XY x V0
  fCandidateVariableNames[55]="newLcDCA";
  fCandidateVariableNames[56]="xVtxLcBad";
  fCandidateVariableNames[57]="yVtxLcBad";
  fCandidateVariableNames[58]="zVtxLcBad";
  fCandidateVariableNames[59]="xVtxLcGood";
  fCandidateVariableNames[60]="yVtxLcGood";
  fCandidateVariableNames[61]="zVtxLcGood";
  fCandidateVariableNames[62]="xVtxLcMC";
  fCandidateVariableNames[63]="yVtxLcMC";
  fCandidateVariableNames[64]="zVtxLcMC";
  fCandidateVariableNames[65]="pxVtxBachelorBad";
  fCandidateVariableNames[66]="pyVtxBachelorBad";
  fCandidateVariableNames[67]="pxVtxBachelorGood";
  fCandidateVariableNames[68]="pyVtxBachelorGood";
  fCandidateVariableNames[69]="pxVtxV0";
  fCandidateVariableNames[70]="pyVtxV0";
  fCandidateVariableNames[71]="xPvtx";
  fCandidateVariableNames[72]="yPvtx";
  fCandidateVariableNames[73]="zPvtx";
  fCandidateVariableNames[74]="cosThetaStarBachelor";
  fCandidateVariableNames[75]="cosThetaStarV0";
  fCandidateVariableNames[76]="etaV0";
  fCandidateVariableNames[77]="yV0";
  fCandidateVariableNames[78]="pzVtxBachelorGood";
  fCandidateVariableNames[79]="pzVtxV0";

  //fCandidateVariableNames[65]="bachelorPx";
  //fCandidateVariableNames[66]="bachelorPy";
  //fCandidateVariableNames[67]="bachelorPz";
  //fCandidateVariableNames[68]="V0positivePx";
  //fCandidateVariableNames[69]="V0positivePy";
  //fCandidateVariableNames[70]="V0positivePz";
  //fCandidateVariableNames[71]="V0negativePx";
  //fCandidateVariableNames[72]="V0negativePy";
  //fCandidateVariableNames[73]="V0negativePz";
  //fCandidateVariableNames[74]="bachelorPxDCA";
  //fCandidateVariableNames[75]="bachelorPyDCA";
  //fCandidateVariableNames[76]="bachelorPzDCA";
  //fCandidateVariableNames[77]="v0PxDCA";
  //fCandidateVariableNames[78]="v0PyDCA";
  //fCandidateVariableNames[79]="v0PzDCA";
  //fCandidateVariableNames[80]="V0positivePxDCA";
  //fCandidateVariableNames[81]="V0positivePyDCA";
  //fCandidateVariableNames[82]="V0positivePzDCA";
  //fCandidateVariableNames[83]="V0negativePxDCA";
  //fCandidateVariableNames[84]="V0negativePyDCA";
  //fCandidateVariableNames[85]="V0negativePzDCA";
  //fCandidateVariableNames[86]="qtLc";
  //fCandidateVariableNames[87]="alphaLc";

  for (Int_t ivar=0; ivar<nVar; ivar++) {
    fVariablesTree->Branch(fCandidateVariableNames[ivar].Data(),&fCandidateVariables[ivar],Form("%s/f",fCandidateVariableNames[ivar].Data()));
  }

  return;
}

//__________________________________________________________________________
void  AliAnalysisTaskSELc2V0bachelor::DefineGeneralHistograms() {
  //
  // This is to define general histograms
  //

  fCEvents = new TH1F("fCEvents","conter",17,0,17);
  fCEvents->SetStats(kTRUE);
  fCEvents->GetXaxis()->SetBinLabel(1,"X1");
  fCEvents->GetXaxis()->SetBinLabel(2,"Analyzed events");
  fCEvents->GetXaxis()->SetBinLabel(3,"AliAODVertex exists");
  fCEvents->GetXaxis()->SetBinLabel(4,"CascadesHF exists");
  fCEvents->GetXaxis()->SetBinLabel(5,"MCarray exists");
  fCEvents->GetXaxis()->SetBinLabel(6,"MCheader exists");
  fCEvents->GetXaxis()->SetBinLabel(7,"GetNContributors()>0");
  fCEvents->GetXaxis()->SetBinLabel(8,"IsEventSelected");
  fCEvents->GetXaxis()->SetBinLabel(9,"triggerClass!=CINT1");
  fCEvents->GetXaxis()->SetBinLabel(10,"triggerMask!=kAnyINT");
  fCEvents->GetXaxis()->SetBinLabel(11,"triggerMask!=kAny");
  fCEvents->GetXaxis()->SetBinLabel(12,"vtxTitle.Contains(Z)");
  fCEvents->GetXaxis()->SetBinLabel(13,"vtxTitle.Contains(3D)");
  fCEvents->GetXaxis()->SetBinLabel(14,"vtxTitle.Doesn'tContain(Z-3D)");
  fCEvents->GetXaxis()->SetBinLabel(15,Form("zVtx<=%2.0fcm",fAnalCuts->GetMaxVtxZ()));
  fCEvents->GetXaxis()->SetBinLabel(16,"!IsEventSelected");
  fCEvents->GetXaxis()->SetBinLabel(17,"triggerMask!=kAnyINT || triggerClass!=CINT1");
  //fCEvents->GetXaxis()->SetTitle("");
  fCEvents->GetYaxis()->SetTitle("counts");

  fOutput->Add(fCEvents);
  TString fillthis="";

  if (fUseMCInfo && fAdditionalChecks) {
    fillthis="histMcStatLc";
    TH1F* mcStatisticLc = new TH1F(fillthis.Data(),"#Lambda_{c} generated and their decays",21,-10.5,10.5);
    fOutput->Add(mcStatisticLc);
  }

  //fillthis="histopionV0SigmaVspTOF";
  //TH2F *hpionV0SigmaVspTOF=new TH2F(fillthis.Data(),fillthis.Data(),300,0.,30.,100,-5.,5.);
  fillthis="histoprotonBachSigmaVspTOF";
  TH2F *hprotonBachSigmaVspTOF=new TH2F(fillthis.Data(),fillthis.Data(),300,0.,30.,100,-5.,5.);

  //fOutput->Add(hpionV0SigmaVspTOF);
  fOutput->Add(hprotonBachSigmaVspTOF);

  //fillthis="histopionV0SigmaVspTPC";
  //TH2F *hpionV0SigmaVspTPC=new TH2F(fillthis.Data(),fillthis.Data(),300,0.,30.,100,-5.,5.);
  fillthis="histoprotonBachSigmaVspTPC";
  TH2F *hprotonBachSigmaVspTPC=new TH2F(fillthis.Data(),fillthis.Data(),300,0.,30.,100,-5.,5.);

  //fOutput->Add(hpionV0SigmaVspTPC);
  fOutput->Add(hprotonBachSigmaVspTPC);

  if (fUseMCInfo) {

    //fillthis="histopionV0SigmaVspTOFsgn";
    //TH2F *hpionV0SigmaVspTOFsgn=new TH2F(fillthis.Data(),fillthis.Data(),300,0.,30.,100,-5.,5.);
    fillthis="histoprotonBachSigmaVspTOFsgn";
    TH2F *hprotonBachSigmaVspTOFsgn=new TH2F(fillthis.Data(),fillthis.Data(),300,0.,30.,100,-5.,5.);

    //fOutput->Add(hpionV0SigmaVspTOFsgn);
    fOutput->Add(hprotonBachSigmaVspTOFsgn);

    //fillthis="histopionV0SigmaVspTPCsgn";
    //TH2F *hpionV0SigmaVspTPCsgn=new TH2F(fillthis.Data(),fillthis.Data(),300,0.,30.,100,-5.,5.);
    fillthis="histoprotonBachSigmaVspTPCsgn";
    TH2F *hprotonBachSigmaVspTPCsgn=new TH2F(fillthis.Data(),fillthis.Data(),300,0.,30.,100,-5.,5.);

    //fOutput->Add(hpionV0SigmaVspTPCsgn);
    fOutput->Add(hprotonBachSigmaVspTPCsgn);


    //fillthis="histopionV0SigmaVspTOFbkg";
    //TH2F *hpionV0SigmaVspTOFbkg=new TH2F(fillthis.Data(),fillthis.Data(),300,0.,30.,100,-5.,5.);
    fillthis="histoprotonBachSigmaVspTOFbkg";
    TH2F *hprotonBachSigmaVspTOFbkg=new TH2F(fillthis.Data(),fillthis.Data(),300,0.,30.,100,-5.,5.);

    //fOutput->Add(hpionV0SigmaVspTOFbkg);
    fOutput->Add(hprotonBachSigmaVspTOFbkg);

    //fillthis="histopionV0SigmaVspTPCbkg";
    //TH2F *hpionV0SigmaVspTPCbkg=new TH2F(fillthis.Data(),fillthis.Data(),300,0.,30.,100,-5.,5.);
    fillthis="histoprotonBachSigmaVspTPCbkg";
    TH2F *hprotonBachSigmaVspTPCbkg=new TH2F(fillthis.Data(),fillthis.Data(),300,0.,30.,100,-5.,5.);

    //fOutput->Add(hpionV0SigmaVspTPCbkg);
    fOutput->Add(hprotonBachSigmaVspTPCbkg);

  }

  if (fAdditionalChecks) {

    TH1F *hZ2 = new TH1F("hZ2","",100,-50.,50.);
    fOutput->Add(hZ2);
    TH1F *hZ3 = new TH1F("hZ3","",100,-50.,50.);
    fOutput->Add(hZ3);
    TH1F *hZ4 = new TH1F("hZ4","",100,-50.,50.);
    fOutput->Add(hZ4);
    TH1F *hZ5 = new TH1F("hZ5","",100,-50.,50.);
    fOutput->Add(hZ5);
    TH1F *hZ6 = new TH1F("hZ6","",100,-50.,50.);
    fOutput->Add(hZ6);
    TH1F *hZ7 = new TH1F("hZ7","",100,-50.,50.);
    fOutput->Add(hZ7);
    TH1F *hZ8 = new TH1F("hZ8","",100,-50.,50.);
    fOutput->Add(hZ8);
    TH1F *hZ9 = new TH1F("hZ9","",100,-50.,50.);
    fOutput->Add(hZ9);
    TH1F *hZ10 = new TH1F("hZ10","",100,-50.,50.);
    fOutput->Add(hZ10);
    TH1F *hZ11 = new TH1F("hZ11","",100,-50.,50.);
    fOutput->Add(hZ11);
    TH1F *hZ12 = new TH1F("hZ12","",100,-50.,50.);
    fOutput->Add(hZ12);
    TH1F *hZ13 = new TH1F("hZ13","",100,-50.,50.);
    fOutput->Add(hZ13);
    TH1F *hZ14 = new TH1F("hZ14","",100,-50.,50.);
    fOutput->Add(hZ14);
    TH1F *hZ15 = new TH1F("hZ15","",100,-50.,50.);
    fOutput->Add(hZ15);
    TH1F *hZ16 = new TH1F("hZ16","",100,-50.,50.);
    fOutput->Add(hZ16);
  }

  TH1F *hCandidateSelection = new TH1F("hCandidateSelection","",10,-0.5,9.5);
  hCandidateSelection->GetXaxis()->SetBinLabel(1,"IsEventSelected");
  hCandidateSelection->GetXaxis()->SetBinLabel(2,"IsSecondaryVtx");
  hCandidateSelection->GetXaxis()->SetBinLabel(3,"V0toPosNeg");
  hCandidateSelection->GetXaxis()->SetBinLabel(4,"offlineV0");
  hCandidateSelection->GetXaxis()->SetBinLabel(5,"isInFiducialAcceptance");
  hCandidateSelection->GetXaxis()->SetBinLabel(6,"analCuts::kTracks");
  hCandidateSelection->GetXaxis()->SetBinLabel(7,"analCuts::kCandidateNoPID");
  hCandidateSelection->GetXaxis()->SetBinLabel(8,"analCuts::kPID");
  hCandidateSelection->GetXaxis()->SetBinLabel(9,"analCuts::kCandidateWithPID");
  hCandidateSelection->GetXaxis()->SetBinLabel(10,"analCuts::kAll");
  fOutput->Add(hCandidateSelection);

  TH1F *hEventsWithCandidates = new TH1F("hEventsWithCandidates","conter",11,5.5,16.5);
  hEventsWithCandidates->GetXaxis()->SetBinLabel(1,"GetNContributors()>0");
  hEventsWithCandidates->GetXaxis()->SetBinLabel(2,"IsEventSelected");
  hEventsWithCandidates->GetXaxis()->SetBinLabel(3,"triggerClass!=CINT1");
  hEventsWithCandidates->GetXaxis()->SetBinLabel(4,"triggerMask!=kAnyINT");
  hEventsWithCandidates->GetXaxis()->SetBinLabel(5,"triggerMask!=kAny");
  hEventsWithCandidates->GetXaxis()->SetBinLabel(6,"vtxTitle.Contains(Z)");
  hEventsWithCandidates->GetXaxis()->SetBinLabel(7,"vtxTitle.Contains(3D)");
  hEventsWithCandidates->GetXaxis()->SetBinLabel(8,"vtxTitle.Doesn'tContain(Z-3D)");
  hEventsWithCandidates->GetXaxis()->SetBinLabel(9,Form("zVtx<=%2.0fcm",fAnalCuts->GetMaxVtxZ()));
  hEventsWithCandidates->GetXaxis()->SetBinLabel(10,"!IsEventSelected");
  hEventsWithCandidates->GetXaxis()->SetBinLabel(11,"triggerMask!=kAnyINT || triggerClass!=CINT1");
  fOutput->Add(hEventsWithCandidates);

  if (fAdditionalChecks) {

    TH1F *hZ6a = new TH1F("hZ6a","",100,-50.,50.);
    fOutput->Add(hZ6a);
    TH1F *hZ7a = new TH1F("hZ7a","",100,-50.,50.);
    fOutput->Add(hZ7a);
    TH1F *hZ8a = new TH1F("hZ8a","",100,-50.,50.);
    fOutput->Add(hZ8a);
    TH1F *hZ9a = new TH1F("hZ9a","",100,-50.,50.);
    fOutput->Add(hZ9a);
    TH1F *hZ10a = new TH1F("hZ10a","",100,-50.,50.);
    fOutput->Add(hZ10a);
    TH1F *hZ11a = new TH1F("hZ11a","",100,-50.,50.);
    fOutput->Add(hZ11a);
    TH1F *hZ12a = new TH1F("hZ12a","",100,-50.,50.);
    fOutput->Add(hZ12a);
    TH1F *hZ13a = new TH1F("hZ13a","",100,-50.,50.);
    fOutput->Add(hZ13a);
    TH1F *hZ14a = new TH1F("hZ14a","",100,-50.,50.);
    fOutput->Add(hZ14a);
    TH1F *hZ15a = new TH1F("hZ15a","",100,-50.,50.);
    fOutput->Add(hZ15a);
    TH1F *hZ16a = new TH1F("hZ16a","",100,-50.,50.);
    fOutput->Add(hZ16a);
  }

  TH1F *hSwitchOnCandidates1 = new TH1F("hSwitchOnCandidates1","",15,-7.5,7.5);
  fOutput->Add(hSwitchOnCandidates1);
  TH1F *hSwitchOnCandidates2 = new TH1F("hSwitchOnCandidates2","",15,-7.5,7.5);
  fOutput->Add(hSwitchOnCandidates2);
  TH1F *hSwitchOnCandidates3 = new TH1F("hSwitchOnCandidates3","",15,-7.5,7.5);
  fOutput->Add(hSwitchOnCandidates3);
  TH1F *hSwitchOnCandidates4 = new TH1F("hSwitchOnCandidates4","",15,-7.5,7.5);
  fOutput->Add(hSwitchOnCandidates4);

  return;
}

//________________________________________________________________________
void  AliAnalysisTaskSELc2V0bachelor::DefineAnalysisHistograms() {
  //
  // This is to define analysis histograms
  //

  if (fIsK0SAnalysis) DefineK0SHistos();// hK0S histos declarations

  return;
}

//________________________________________________________________________
void  AliAnalysisTaskSELc2V0bachelor::FillAnalysisHistograms(AliAODRecoCascadeHF *part, Bool_t isBachelorID, TString appendthis) {
  //
  // This is to fill analysis histograms
  //

  TString fillthis="";

  Double_t invmassLc = part->InvMassLctoK0sP();
  Double_t lambdacpt = part->Pt();

  AliAODTrack *bachelor = (AliAODTrack*)part->GetBachelor();
  Double_t momBach  = bachelor->P();

  AliAODv0 * v0part = (AliAODv0*)part->Getv0();
  Double_t momK0S = v0part->P();
  Double_t ptK0S = v0part->Pt();
  Double_t dcaV0ptp = v0part->GetDCA();
  Double_t invmassK0S = v0part->MassK0Short();

    fillthis="histK0SMass"+appendthis;
    cout << fillthis << endl;
    ((TH2F*)(fOutputAll->FindObject(fillthis)))->Fill(invmassK0S,ptK0S);
    if (isBachelorID)  ((TH2F*)(fOutputPIDBach->FindObject(fillthis)))->Fill(invmassK0S,ptK0S);

    fillthis="histpK0Svsp"+appendthis;
    cout << fillthis << endl;
    ((TH2F*)(fOutputAll->FindObject(fillthis)))->Fill(momBach,momK0S);
    if (isBachelorID)  ((TH2F*)(fOutputPIDBach->FindObject(fillthis)))->Fill(momBach,momK0S);

    fillthis="histDCAtoPVvspK0S"+appendthis;
    cout << fillthis << endl;
    ((TH2F*)(fOutputAll->FindObject(fillthis)))->Fill(momK0S,dcaV0ptp);
    if (isBachelorID)  ((TH2F*)(fOutputPIDBach->FindObject(fillthis)))->Fill(momK0S,dcaV0ptp);

    fillthis="histArmPodK0S"+appendthis;
    cout << fillthis << endl;
    FillArmPodDistribution(v0part,fillthis,fOutputAll);
    if (isBachelorID) FillArmPodDistribution(v0part,fillthis,fOutputPIDBach);

    fillthis="histLcMassByK0S"+appendthis;
    cout << fillthis << endl;
    ((TH2F*)(fOutputAll->FindObject(fillthis)))->Fill(invmassLc,lambdacpt);
    if (isBachelorID)((TH2F*)(fOutputPIDBach->FindObject(fillthis)))->Fill(invmassLc,lambdacpt);

    return;
}
//---------------------------
Double_t AliAnalysisTaskSELc2V0bachelor::PropagateToDCA(AliAODv0 *v, AliAODTrack *bachelor, Double_t b,
							Double_t &xVtxLc, Double_t &yVtxLc, Double_t &zVtxLc,
							Double_t &pxVtxBachelor, Double_t &pyVtxBachelor, Double_t &pzVtxBachelor) {
  //--------------------------------------------------------------------
  // This function returns the DCA between the V0 and the track
  // This is a copy of AliCascadeVertexer::PropagateToDCA(...) method
  //--------------------------------------------------------------------

  // Get AliExternalTrackParam out of the AliAODTracks                                                      
  Double_t xyz[3], pxpypz[3], cv[21]; Short_t sign;
  bachelor->PxPyPz(pxpypz);
  bachelor->XvYvZv(xyz);
  bachelor->GetCovarianceXYZPxPyPz(cv);
  sign=bachelor->Charge();
  AliExternalTrackParam *t = new AliExternalTrackParam(xyz,pxpypz,cv,sign);

  Double_t alpha=t->GetAlpha(), cs1=TMath::Cos(alpha), sn1=TMath::Sin(alpha);
  //Double_t alpha = GetAlpha(xyz,pxpypz), cs1=TMath::Cos(alpha), sn1=TMath::Sin(alpha);

  // position and momentum of bachelor
  Double_t x1=xyz[0], y1=xyz[1], z1=xyz[2];
  Double_t px1=pxpypz[0], py1=pxpypz[1], pz1=pxpypz[2];

  // position and momentum of V0
  Double_t x2=v->DecayVertexV0X(),
    y2=v->DecayVertexV0Y(),
    z2=v->DecayVertexV0Z();
  Double_t px2=v->Px(),
    py2=v->Py(),
    pz2=v->Pz();

  /*
  AliAODTrack *trackP = (AliAODTrack*) v->GetDaughter(0);
  //Double_t pxpypzP[3]; trackP->PxPyPz(pxpypzP);
  //Double_t xyzP[3]; trackP->XvYvZv(xyzP);
  Double_t cvP[21]; trackP->GetCovarianceXYZPxPyPz(cvP);
  //Short_t signP=trackP->Charge();
  //AliExternalTrackParam *tP = new AliExternalTrackParam(xyzP,pxpypzP,cvP,signP);

  // Get AliExternalTrackParam out of the AliAODTrack
  AliAODTrack *trackN = (AliAODTrack*) v->GetDaughter(1);
  //Double_t pxpypzN[3]; trackN->PxPyPz(pxpypzN);
  //Double_t xyzN[3]; trackN->XvYvZv(xyzN);
  Double_t cvN[21]; trackN->GetCovarianceXYZPxPyPz(cvN);
  //Short_t signN=trackN->Charge();
  //AliExternalTrackParam *tN = new AliExternalTrackParam(xyzN,pxpypzN,cvN,signN);

  Double_t xyzV0[3]={x2,y2,z2};
  Double_t pxpypzV0[3]={px2,py2,pz2};
  Double_t cvV0[21]; for (Int_t ii=0; ii<21; ii++) cvV0[ii]=cvP[ii]+cvN[ii];
  AliNeutralTrackParam *trackV0 = new AliNeutralTrackParam(xyzV0,pxpypzV0,cvV0,0);
  */
 
  // calculation dca
  Double_t dd= Det(x2-x1,y2-y1,z2-z1,px1,py1,pz1,px2,py2,pz2);
  Double_t ax= Det(py1,pz1,py2,pz2);
  Double_t ay=-Det(px1,pz1,px2,pz2);
  Double_t az= Det(px1,py1,px2,py2);

  Double_t dca=TMath::Abs(dd)/TMath::Sqrt(ax*ax + ay*ay + az*az);

  // bachelor point @ the DCA
  Double_t t1 = Det(x2-x1,y2-y1,z2-z1,px2,py2,pz2,ax,ay,az)/
    Det(px1,py1,pz1,px2,py2,pz2,ax,ay,az);
  x1 += px1*t1; y1 += py1*t1; z1 += pz1*t1;

  //propagate track to the point of DCA
  Double_t rho1=x1*cs1 + y1*sn1;
  if (!t->PropagateTo(rho1,b)) {
    Error("PropagateToDCA","Propagation failed !");
    delete t; t=NULL;
    return 1.e+33;
  }

  Double_t pBachelorDCA[3]; t->GetPxPyPz(pBachelorDCA);
  pxVtxBachelor=pBachelorDCA[0], pyVtxBachelor=pBachelorDCA[1], pzVtxBachelor=pBachelorDCA[2];

  delete t; t=NULL;

  // V0 point @ the DCA
  Double_t t2 = Det(x1-x2,y1-y2,z1-z2,px1,py1,pz1,ax,ay,az)/
    Det(px2,py2,pz2,px1,py1,pz1,ax,ay,az);
  x2 += px2*t2; y2 += py2*t2; z2 += pz2*t2;


  // Lc decay vtx
  xVtxLc = 0.5*(x1+x2);
  yVtxLc = 0.5*(y1+y2);
  zVtxLc = 0.5*(z1+z2);
  
  return dca;

}

//---------------------------
Double_t AliAnalysisTaskSELc2V0bachelor::GetAlpha(Double_t xyz[3],Double_t pxpypz[3])
{
  //
  // To estimate alpha according to what done in the AliExternalTrackParam::Set(...) method
  //

  Double_t alpha = 0.;

  const double kSafe = 1e-5;
  Double_t radPos2 = xyz[0]*xyz[0]+xyz[1]*xyz[1];
  Double_t radMax  = 45.; // approximately ITS outer radius
  if (radPos2 < radMax*radMax) { // inside the ITS
    alpha = TMath::ATan2(pxpypz[1],pxpypz[0]);
  } else { // outside the ITS
    Float_t phiPos = TMath::Pi()+TMath::ATan2(-xyz[1], -xyz[0]);
     alpha =
       TMath::DegToRad()*(20*((((Int_t)(phiPos*TMath::RadToDeg()))/20))+10);
  }

  Double_t cs=TMath::Cos(alpha), sn=TMath::Sin(alpha);
  // protection: avoid alpha being too close to 0 or +-pi/2
  if (TMath::Abs(sn)<2*kSafe) {
    if (alpha>0) alpha += alpha< TMath::Pi()/2. ?  2*kSafe : -2*kSafe;
    else         alpha += alpha>-TMath::Pi()/2. ? -2*kSafe :  2*kSafe;
    cs=TMath::Cos(alpha);
    sn=TMath::Sin(alpha);
  }
  else if (TMath::Abs(cs)<2*kSafe) {
    if (alpha>0) alpha += alpha> TMath::Pi()/2. ? 2*kSafe : -2*kSafe;
    else         alpha += alpha>-TMath::Pi()/2. ? 2*kSafe : -2*kSafe;
    cs=TMath::Cos(alpha);
    sn=TMath::Sin(alpha);
  }


  return alpha;
}

//---------------------------
Double_t AliAnalysisTaskSELc2V0bachelor::Det(Double_t a00, Double_t a01,
					     Double_t a10, Double_t a11) const {
  //--------------------------------------------------------------------
  // This function calculates locally a 2x2 determinant.
  // This is a copy of the AliCascadeVertexer::Det(...) method
  //--------------------------------------------------------------------
  return a00*a11 - a01*a10;
}

//---------------------------
Double_t AliAnalysisTaskSELc2V0bachelor::Det(Double_t a00,Double_t a01,Double_t a02,
					     Double_t a10,Double_t a11,Double_t a12,
					     Double_t a20,Double_t a21,Double_t a22) const {
  //--------------------------------------------------------------------
  // This function calculates locally a 3x3 determinant
  // This is a copy of the AliCascadeVertexer::Det(...) method
  //--------------------------------------------------------------------
  return  a00*Det(a11,a12,a21,a22)-a01*Det(a10,a12,a20,a22)+a02*Det(a10,a11,a20,a21);
}

