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
#include <TH2F.h>
#include <TH3F.h>
#include <THnSparse.h>
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
#include "AliInputEventHandler.h"
#include "AliESDtrackCuts.h"
#include "AliNeutralTrackParam.h"
#include "AliVertexingHFUtils.h"

using std::cout;
using std::endl;

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskSELc2V0bachelor);
/// \endcond

//__________________________________________________________________________
AliAnalysisTaskSELc2V0bachelor::AliAnalysisTaskSELc2V0bachelor() : AliAnalysisTaskSE(),
  fUseMCInfo(kFALSE),
  fOutput(0),
  fOutputAll(0),
  fOutputPIDBach(0),
  fCEvents(0),
  fCounter(0),
  fAnalCuts(0),
  fUseOnTheFlyV0(kFALSE),
  fAODProtection(1),
  fIsEventSelected(kFALSE),
  fWriteVariableTree(kFALSE),
  fVariablesTree(0),
  fCandidateVariables(),
  fVtx1(0),
  fBzkG(0),
  fAdditionalChecks(kFALSE),
  fTrackRotation(kFALSE),
  fOutputPIDBachTR(0),
  fMinAngleForRot(5*TMath::Pi()/6),
  fMaxAngleForRot(7*TMath::Pi()/6),
  fMinMass(0),
  fMaxMass(0),
  fNRotations(9),
  fPtMinToFillTheTree(0.),
  fPtMaxToFillTheTree(999.),
  fUseTPCPIDtoFillTree(kFALSE),
  fSign(2),
  fCheckOrigin(kFALSE),
  fReconstructSecVtx(kFALSE)
{
  //
  /// Default ctor
  //

  Double_t mLcPDG = TDatabasePDG::Instance()->GetParticle(4122)->Mass();
  fMinMass=mLcPDG-0.250;
  fMaxMass=mLcPDG+0.250;

}
//___________________________________________________________________________
AliAnalysisTaskSELc2V0bachelor::AliAnalysisTaskSELc2V0bachelor(const Char_t* name,
								   AliRDHFCutsLctoV0* analCuts, Bool_t useOnTheFly,
								   Bool_t writeVariableTree, Bool_t additionalChecks, Bool_t trackRotation, Bool_t useTPCpid, Char_t sign, Bool_t origin) :
  AliAnalysisTaskSE(name),
  fUseMCInfo(kFALSE),
  fOutput(0),
  fOutputAll(0),
  fOutputPIDBach(0),
  fCEvents(0),
  fCounter(0),
  fAnalCuts(analCuts),
  fUseOnTheFlyV0(useOnTheFly),
  fAODProtection(1),
  fIsEventSelected(kFALSE),
  fWriteVariableTree(writeVariableTree),
  fVariablesTree(0),
  fCandidateVariables(),
  fVtx1(0),
  fBzkG(0),
  fAdditionalChecks(additionalChecks),
  fTrackRotation(trackRotation),
  fOutputPIDBachTR(0),
  fMinAngleForRot(5*TMath::Pi()/6),
  fMaxAngleForRot(7*TMath::Pi()/6),
  fMinMass(0),
  fMaxMass(0),
  fNRotations(9),
  fPtMinToFillTheTree(0.),
  fPtMaxToFillTheTree(999.),
  fUseTPCPIDtoFillTree(useTPCpid),
  fSign(sign),
  fCheckOrigin(origin),
  fReconstructSecVtx(kFALSE)
{
  //
  /// Constructor. Initialization of Inputs and Outputs
  //
  Info("AliAnalysisTaskSELc2V0bachelor","Calling Constructor");

  if (fWriteVariableTree && fTrackRotation) {
    AliInfo(Form("You cannot initialize fWriteVariableTree=%d and fTrackRotation=%d => fTrackRotation=0",fWriteVariableTree,fTrackRotation));
    fTrackRotation=kFALSE;
  }

  Double_t mLcPDG = TDatabasePDG::Instance()->GetParticle(4122)->Mass();
  fMinMass=mLcPDG-0.250;
  fMaxMass=mLcPDG+0.250;

  DefineOutput(1,TList::Class());  //conters
  DefineOutput(2,AliNormalizationCounter::Class());
  DefineOutput(3,AliRDHFCutsLctoV0::Class());
  if (!writeVariableTree) {
    DefineOutput(4,TList::Class());  //All Entries output
    DefineOutput(5,TList::Class());  //3sigma PID output
    if (trackRotation) {
      DefineOutput(6,TList::Class());  //All Entries output
    }
  } else {
    // Output slot #4 keeps a tree of the candidate variables after track selection
    DefineOutput(4,TTree::Class());  //My private output
  }

  if (fWriteVariableTree) fSign=2;

}

//___________________________________________________________________________
AliAnalysisTaskSELc2V0bachelor::~AliAnalysisTaskSELc2V0bachelor() {
  //
  /// destructor
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

  if (fOutputPIDBachTR) {
    delete fOutputPIDBachTR;
    fOutputPIDBachTR = 0;
  }

}
//_________________________________________________
void AliAnalysisTaskSELc2V0bachelor::Init() {
  //
  /// Initialization
  //

  fIsEventSelected=kFALSE;

  if (fDebug > 1) AliInfo("Init");

  PostData(3,fAnalCuts);

  return;
}

//_________________________________________________
void AliAnalysisTaskSELc2V0bachelor::UserExec(Option_t *)
{
  /// user exec
  if (!fInputEvent) {
    AliError("NO EVENT FOUND!");
    return;
  }

  AliAODEvent* aodEvent = dynamic_cast<AliAODEvent*>(fInputEvent);

  if(fAODProtection>=0)
  {
    //   Protection against different number of events in the AOD and deltaAOD
    //   In case of discrepancy the event is rejected.
    Int_t matchingAODdeltaAODlevel = AliRDHFCuts::CheckMatchingAODdeltaAODevents();
    if (matchingAODdeltaAODlevel<0 || (matchingAODdeltaAODlevel==0 && fAODProtection==1)) {
      // AOD/deltaAOD trees have different number of entries || TProcessID do not match while it was required
      fCEvents->Fill(19);
      return;
    }
  }

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
  fAnalCuts->SetMagneticField(fBzkG);
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

    Double_t zMCVertex = mcHeader->GetVtxZ();
    if (TMath::Abs(zMCVertex) > fAnalCuts->GetMaxVtxZ()) {
      AliDebug(2,Form("Event rejected: abs(zVtxMC)=%f > fAnalCuts->GetMaxVtxZ()=%f",zMCVertex,fAnalCuts->GetMaxVtxZ()));
      return;
    } else {
      fCEvents->Fill(17); // in case of MC events
    }
  }

  Int_t runnumber = aodEvent->GetRunNumber();
  if (aodEvent->GetTriggerMask() == 0 && (runnumber >= 195344 && runnumber <= 195677)){
    AliDebug(3,"Event rejected because of null trigger mask");
    return;
  }

  fCounter->StoreEvent(aodEvent,fAnalCuts,fUseMCInfo); // it is very important that it stays BEFORE any other event selection

  if (fVtx1->GetNContributors()>0) // this check is done in IsEventSelected
    fCEvents->Fill(6);

  if ( !fIsEventSelected ) return; // don't take into account not selected events
  fCEvents->Fill(7);

  Int_t nSelectedAnal = 0;
  MakeAnalysisForLc2prK0S(aodEvent,arrayLctopKos,mcArray, nSelectedAnal, fAnalCuts);

  if (nSelectedAnal)
    CheckEventSelectionWithCandidates(aodEvent);

  fCounter->StoreCandidates(aodEvent,nSelectedAnal,kTRUE);
  fCounter->StoreCandidates(aodEvent,nSelectedAnal,kFALSE);

  PostData(1,fOutput);
  PostData(2,fCounter);
  if (!fWriteVariableTree) {
    PostData(4,fOutputAll);
    PostData(5,fOutputPIDBach);
    if (fTrackRotation)
      PostData(6,fOutputPIDBachTR);
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

    if (fTrackRotation) {
      fOutputPIDBachTR = dynamic_cast<TList*> (GetOutputData(6));
      if (!fOutputPIDBachTR) {
	AliError("fOutputPIDBachTR not available");
	return;
      }
    }

  } else {
    fVariablesTree = dynamic_cast<TTree*> (GetOutputData(4));
    if (!fVariablesTree) {
      AliError("fVariablesTree not available");
      return;
    }
  }

  return;
}
//___________________________________________________________________________
void AliAnalysisTaskSELc2V0bachelor::UserCreateOutputObjects() {
  /// output
  AliInfo(Form("CreateOutputObjects of task %s\n", GetName()));

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

    if (fTrackRotation) {
      fOutputPIDBachTR = new TList();
      fOutputPIDBachTR->SetOwner();
      fOutputPIDBachTR->SetName("listPIDBachTR");
    }

    DefineK0SHistos(); // define analysis histograms
    if (fUseMCInfo && fCheckOrigin) DefineSignalHistosSeparatedPerOrigin(); // define analysis histograms for SNG separated for origin

    PostData(4,fOutputAll);
    PostData(5,fOutputPIDBach);

    if (fTrackRotation)
      PostData(6,fOutputPIDBachTR);

  }
  else {
    DefineTreeVariables();
    PostData(4,fVariablesTree);
  }

  return;
}

//-------------------------------------------------------------------------------
void AliAnalysisTaskSELc2V0bachelor::MakeAnalysisForLc2prK0S(AliAODEvent *aodEvent,TClonesArray *arrayLctopKos,
							     TClonesArray *mcArray,
							     Int_t &nSelectedAnal,
							     AliRDHFCutsLctoV0 *cutsAnal)
{

  /// make the analysis

  Int_t pdgCand = 4122;
  Int_t pdgDgLctoV0bachelor[2]={2212,310}; // always 1st bachelor, 2nd V0
  Int_t pdgDgV0toDaughters[2]={211,211};

  // loop over cascades to search for candidates Lc->p+K0S
  Int_t nCascades= arrayLctopKos->GetEntriesFast();
  if (nCascades==0) {
    AliInfo("Could not find cascades, skipping the event");
    return;
  }
  AliAnalysisVertexingHF *vHF = new AliAnalysisVertexingHF();
  for (Int_t iLctopK0S = 0; iLctopK0S<nCascades; iLctopK0S++) {

    ((TH1F*)(fOutput->FindObject("hCandidateSelection")))->Fill(0);

    // Lc candidates and K0S from Lc
    AliAODRecoCascadeHF* lcK0Spr = dynamic_cast<AliAODRecoCascadeHF*>(arrayLctopKos->At(iLctopK0S));
    if (!lcK0Spr) {
      AliDebug(2,Form("Cascade %d doens't exist, skipping",iLctopK0S));
      continue;
    }

    if (!(lcK0Spr->CheckCascadeFlags())) {
      AliDebug(2,Form("Cascade %d is not flagged as Lc candidate",iLctopK0S));
      continue;
    }

    Bool_t unsetvtx=kFALSE;
    if (!lcK0Spr->GetOwnPrimaryVtx()) {
      lcK0Spr->SetOwnPrimaryVtx(fVtx1);
      unsetvtx=kTRUE;
    }

    if(!vHF->FillRecoCasc(aodEvent,lcK0Spr,kFALSE)){//Fill the data members of the candidate only if they are empty.
       fCEvents->Fill(18);//monitor how often this fails
    continue;
    }
    if (fReconstructSecVtx) {
      if (!(vHF->RecoSecondaryVertexForCascades(aodEvent, lcK0Spr))) {
	continue;
      }
    }
    if (!lcK0Spr->GetSecondaryVtx()) {
      AliInfo("No secondary vertex"); // it will be done in AliRDHFCutsLctoV0::IsSelected
      continue;
    }

    if (lcK0Spr->GetNDaughters()!=2) {
      AliDebug(2,Form("Cascade %d has not 2 daughters (nDaughters=%d)",iLctopK0S,lcK0Spr->GetNDaughters())); // it will be done in AliRDHFCutsLctoV0::IsSelected
      continue;
    }

    if ( (fSign == 0 && lcK0Spr->Charge()<0) ||
	 (fSign == 1 && lcK0Spr->Charge()>0) ) {
      AliDebug(2,Form("Charge of the cascade %d is different with respect to the required one",iLctopK0S)); //
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
    if (!v0Neg || !v0Pos) {
      AliDebug(2,Form("V0 by cascade %d has no V0positive of V0negative object",iLctopK0S)); // it will be done in AliRDHFCutsLctoV0::IsSelected
      continue;
    }

    ((TH1F*)(fOutput->FindObject("hCandidateSelection")))->Fill(1);

    if (v0Pos->Charge() == v0Neg->Charge()) continue;

    ((TH1F*)(fOutput->FindObject("hCandidateSelection")))->Fill(2);

    Int_t isLc = 0;
    Int_t originLc = -1; // -1 -> RD; 0 -> BKG; 1 -> form c; 2 -> from b; 3 -> not from quark

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
	  AliVertexingHFUtils *util = new AliVertexingHFUtils();
	  Int_t pdgMom = util->CheckOrigin(mcArray,partLc,kFALSE);
	  if (pdgMom == 4) {
	    originLc=1; // from c
	  } else if (pdgMom == 5) {
	    originLc=2; // form b
	  } else {
	    Int_t isThereaQuark=util->CheckOrigin(mcArray,partLc,kTRUE);
	    if (isThereaQuark<=0) originLc=3; // not from q
	  }
	  delete util;
	}
      } else {
	AliDebug(2,Form("No MC candidate (cascade number %d -total cascade number = %d -)", iLctopK0S,nCascades));
	pdgCode=-1;
      }
    } else {
      originLc=-1; // real data
    }

    FillLc2pK0Sspectrum(lcK0Spr, isLc,
			nSelectedAnal, cutsAnal,
			mcArray, originLc);

    if (unsetvtx) lcK0Spr->UnsetOwnPrimaryVtx();

  }

  delete vHF;
  AliDebug(2, Form("Found %d Reco particles that are Lc!!", nSelectedAnal));

  return;
}

//________________________________________________________________________
void AliAnalysisTaskSELc2V0bachelor::FillLc2pK0Sspectrum(AliAODRecoCascadeHF *part,
							   Int_t isLc,
							   Int_t &nSelectedAnal,
							   AliRDHFCutsLctoV0 *cutsAnal,
							   TClonesArray *mcArray,
							   Int_t originLc)
{
  //
  /// Fill histos for Lc -> K0S+proton
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

  if ( onFlyV0 && !fUseOnTheFlyV0 ) return;

  if (fAdditionalChecks) CheckCandidatesAtDifferentLevels(part,cutsAnal);

  // track rotation
  if (fTrackRotation) {
    if (onFlyV0) {
      TrackRotation(cutsAnal,part,"");
    }
    else {
      TrackRotation(cutsAnal,part,"Offline");
    }
    if (fUseMCInfo) {
      if (isLc==1) {
	if (onFlyV0) {
	  TrackRotation(cutsAnal,part,"Sgn");
	}
	else {
	  TrackRotation(cutsAnal,part,"OfflineSgn");
	}
      }// sgn
      else { // bkg
	if (onFlyV0) {
	  TrackRotation(cutsAnal,part,"Bkg");
	}
	else {
	  TrackRotation(cutsAnal,part,"OfflineBkg");
	}
      }
    } // if fUseMCInfo
  } // if fTrackRotation




  if ( !(cutsAnal->IsInFiducialAcceptance(part->Pt(),part->Y(4122))) ) return;

  if ( !( ( (cutsAnal->IsSelected(part,AliRDHFCuts::kTracks))&(AliRDHFCutsLctoV0::kLcToK0Spr)) == (AliRDHFCutsLctoV0::kLcToK0Spr) ) ) return;

  if ( ( ( (cutsAnal->IsSelected(part,AliRDHFCuts::kAll))&(AliRDHFCutsLctoV0::kLcToK0Spr)) == (AliRDHFCutsLctoV0::kLcToK0Spr) ) ) nSelectedAnal++;

  // Fill candidate variable Tree (track selection, V0 invMass selection)
  if ( fWriteVariableTree ) {
    Double_t invmassK0S = v0part->MassK0Short();
    Double_t mk0sPDG = TDatabasePDG::Instance()->GetParticle(310)->Mass();
    Double_t nSigmaTPCpr=-999.;
    cutsAnal->GetPidHF()->GetnSigmaTPC(bachelor,4,nSigmaTPCpr);
    if ( !onFlyV0 && isInCascadeWindow &&
	 part->CosV0PointingAngle()>0.99 && TMath::Abs(invmassK0S-mk0sPDG)<=0.05 &&
	 part->Pt()>=fPtMinToFillTheTree && part->Pt()<fPtMaxToFillTheTree &&
	 (!fUseTPCPIDtoFillTree || (fUseTPCPIDtoFillTree && TMath::Abs(nSigmaTPCpr)<3.)))
      FillTheTree(part,cutsAnal,mcArray,isLc,originLc);
    return;
  }

  cutsAnal->SetUsePID(kFALSE);
  Bool_t isCandidateSelectedCuts = (((cutsAnal->IsSelected(part,AliRDHFCuts::kCandidate))&(AliRDHFCutsLctoV0::kLcToK0Spr))==(AliRDHFCutsLctoV0::kLcToK0Spr)); // kinematic/topological cuts
  cutsAnal->SetUsePID(areCutsUsingPID);
  Bool_t isBachelorID = (((cutsAnal->IsSelected(part,AliRDHFCuts::kPID))&(AliRDHFCutsLctoV0::kLcToK0Spr))==(AliRDHFCutsLctoV0::kLcToK0Spr)); // ID x bachelor

  //if (((cutsAnal->IsSelected(part,AliRDHFCuts::kCandidate))&(AliRDHFCutsLctoV0::kLcToK0Spr))==(AliRDHFCutsLctoV0::kLcToK0Spr)) {
  if (((cutsAnal->IsSelected(part,AliRDHFCuts::kAll))&(AliRDHFCutsLctoV0::kLcToK0Spr))==(AliRDHFCutsLctoV0::kLcToK0Spr)) {
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

  if (onFlyV0) {

    fillthis="histArmPodK0S";
    FillArmPodDistribution(v0part,fillthis,isCandidateSelectedCuts,isBachelorID);

    fillthis="histArmPodLc";
    FillArmPodDistribution(part,fillthis,isCandidateSelectedCuts,isBachelorID);

    //if (isCandidateSelectedCuts) {
    FillAnalysisHistograms(part,cutsAnal,"");
    //}
  }
  else {

    fillthis="histArmPodK0SOffline";
    FillArmPodDistribution(v0part,fillthis,isCandidateSelectedCuts,isBachelorID);

    fillthis="histArmPodLcOffline";
    FillArmPodDistribution(part,fillthis,isCandidateSelectedCuts,isBachelorID);

    FillAnalysisHistograms(part,cutsAnal,"Offline");
    if (isCandidateSelectedCuts) {
      fillthis="histoprotonBachSigmaVspTOF";
      ((TH2F*)(fOutput->FindObject(fillthis)))->Fill(momBach,nSigmaTOFpr);
      fillthis="histoprotonBachSigmaVspTPC";
      ((TH2F*)(fOutput->FindObject(fillthis)))->Fill(momBach,nSigmaTPCpr);
      //FillAnalysisHistograms(part,cutsAnal,"Offline");
    }
  }
  if (fUseMCInfo) {
    if (isLc==1) {
      if (onFlyV0) {

	fillthis="histArmPodK0SSgn";
	FillArmPodDistribution(v0part,fillthis,isCandidateSelectedCuts,isBachelorID);

	fillthis="histArmPodLcSgn";
	FillArmPodDistribution(part,fillthis,isCandidateSelectedCuts,isBachelorID);

	//if (isCandidateSelectedCuts) {
	FillAnalysisHistograms(part,cutsAnal,"Sgn");
	//}
	if (fCheckOrigin) {
	  switch (originLc) {
	  case 1:
	    FillAnalysisHistograms(part,cutsAnal,"SgnC");
	    break;
	  case 2:
	    FillAnalysisHistograms(part,cutsAnal,"SgnB");
	    break;
	  case 3:
	    FillAnalysisHistograms(part,cutsAnal,"SgnNoQ");
	    break;
	  }
	}
      }
      else {

	fillthis="histArmPodK0SOfflineSgn";
	FillArmPodDistribution(v0part,fillthis,isCandidateSelectedCuts,isBachelorID);

	fillthis="histArmPodLcOfflineSgn";
	FillArmPodDistribution(part,fillthis,isCandidateSelectedCuts,isBachelorID);

 	if (isCandidateSelectedCuts) {
	  fillthis="histoprotonBachSigmaVspTOFsgn";
	  ((TH2F*)(fOutput->FindObject(fillthis)))->Fill(momBach,nSigmaTOFpr);
	  fillthis="histoprotonBachSigmaVspTPCsgn";
	  ((TH2F*)(fOutput->FindObject(fillthis)))->Fill(momBach,nSigmaTPCpr);
	  //FillAnalysisHistograms(part,cutsAnal,"OfflineSgn");
	}

	if (fCheckOrigin) {
	  switch (originLc) {
	  case 1:
	    FillAnalysisHistograms(part,cutsAnal,"OfflineSgnC");
	    break;
	  case 2:
	    FillAnalysisHistograms(part,cutsAnal,"OfflineSgnB");
	    break;
	  case 3:
	    FillAnalysisHistograms(part,cutsAnal,"OfflineSgnNoQ");
	    break;
	  }
	}
	FillAnalysisHistograms(part,cutsAnal,"OfflineSgn");

      }
    }// sgn
    else { // bkg
      if (onFlyV0) {

	fillthis="histArmPodK0SBkg";
	FillArmPodDistribution(v0part,fillthis,isCandidateSelectedCuts,isBachelorID);

	fillthis="histArmPodLcBkg";
	FillArmPodDistribution(part,fillthis,isCandidateSelectedCuts,isBachelorID);

	//if (isCandidateSelectedCuts) {
	FillAnalysisHistograms(part,cutsAnal,"Bkg");
	//}
      }
      else {

	fillthis="histArmPodK0SOfflineBkg";
	FillArmPodDistribution(v0part,fillthis,isCandidateSelectedCuts,isBachelorID);

	fillthis="histArmPodLcOfflineBkg";
	FillArmPodDistribution(part,fillthis,isCandidateSelectedCuts,isBachelorID);

	FillAnalysisHistograms(part,cutsAnal,"OfflineBkg");
	if (isCandidateSelectedCuts) {
	  fillthis="histoprotonBachSigmaVspTOFbkg";
	  ((TH2F*)(fOutput->FindObject(fillthis)))->Fill(momBach,nSigmaTOFpr);
	  fillthis="histoprotonBachSigmaVspTPCbkg";
	  ((TH2F*)(fOutput->FindObject(fillthis)))->Fill(momBach,nSigmaTPCpr);
	  //FillAnalysisHistograms(part,cutsAnal,"OfflineBkg");
	}
      }
    }
  } // if fUseMCInfo

  return;
}

//----------------------------------------------------
void AliAnalysisTaskSELc2V0bachelor::DefineK0SHistos()
{

  Double_t mLcPDG  = TDatabasePDG::Instance()->GetParticle(4122)->Mass();
  Double_t mK0SPDG = TDatabasePDG::Instance()->GetParticle(310)->Mass();
  Double_t mMinLambdaPDG  = TDatabasePDG::Instance()->GetParticle(2212)->Mass()+
    TDatabasePDG::Instance()->GetParticle(211)->Mass();

  TString nameHisto=" ", nameHistoSgn=" ", nameHistoBkg=" ";
  TString titleHisto=" ", titleHistoSgn=" ", titleHistoBkg=" ";

  // pt(Lc)
  Double_t *binLimpTLc=new Double_t[11+1]; // 11 pT(Lc) bins
  binLimpTLc[ 0]= 0.;
  binLimpTLc[ 1]= 1.;
  binLimpTLc[ 2]= 2.;
  binLimpTLc[ 3]= 3.;
  binLimpTLc[ 4]= 4.;
  binLimpTLc[ 5]= 5.;
  binLimpTLc[ 6]= 6.;
  binLimpTLc[ 7]= 8.;
  binLimpTLc[ 8]=12.;
  binLimpTLc[ 9]=17.;
  binLimpTLc[10]=25.;
  binLimpTLc[11]=35.;

  // pt(prong)
  Double_t *binLimpTprong=new Double_t[41+1]; // 41 pT(prong) bins
  binLimpTprong[ 0]= 0.0;
  binLimpTprong[ 1]= 0.1;
  binLimpTprong[ 2]= 0.2;
  binLimpTprong[ 3]= 0.3;
  binLimpTprong[ 4]= 0.4;
  binLimpTprong[ 5]= 0.5;
  binLimpTprong[ 6]= 0.6;
  binLimpTprong[ 7]= 0.7;
  binLimpTprong[ 8]= 0.8;
  binLimpTprong[ 9]= 0.9;
  binLimpTprong[10]= 1.0;
  binLimpTprong[11]= 1.2;
  binLimpTprong[12]= 1.4;
  binLimpTprong[13]= 1.6;
  binLimpTprong[14]= 1.8;
  binLimpTprong[15]= 2.0;
  binLimpTprong[16]= 2.2;
  binLimpTprong[17]= 2.4;
  binLimpTprong[18]= 2.6;
  binLimpTprong[19]= 2.8;
  binLimpTprong[20]= 3.0;
  binLimpTprong[21]= 3.5;
  binLimpTprong[22]= 4.0;
  binLimpTprong[23]= 4.5;
  binLimpTprong[24]= 5.0;
  binLimpTprong[25]= 5.5;
  binLimpTprong[26]= 6.0;
  binLimpTprong[27]= 6.5;
  binLimpTprong[28]= 7.0;
  binLimpTprong[29]= 7.5;
  binLimpTprong[30]= 8.0;
  binLimpTprong[31]= 9.0;
  binLimpTprong[32]=10.0;
  binLimpTprong[33]=11.0;
  binLimpTprong[34]=12.0;
  binLimpTprong[35]=13.0;
  binLimpTprong[36]=14.0;
  binLimpTprong[37]=15.0;
  binLimpTprong[38]=20.0;
  binLimpTprong[39]=25.0;
  binLimpTprong[40]=30.0;
  binLimpTprong[41]=35.0;

  if (fUseOnTheFlyV0) {

    // V0 invariant masses (on-the-fly)
    nameHisto="histK0SMass";
    titleHisto="K^{0}_{S} invariant mass VS p_{T}; p_{T}(#Lambda_{c}) [GeV/c]; m_{inv}(#pi^{+},#pi^{-}) [GeV/c^{2}]; Entries";
    TH2F* spectrumK0SMass = new TH2F(nameHisto.Data(),titleHisto.Data(),11,binLimpTLc,1000,mK0SPDG-0.050,mK0SPDG+0.050);

    // Lc invariant masses (x K0S on-the-fly)
    nameHisto="histLcMassByK0S";
    titleHisto="#Lambda_{c} invariant mass (by K^{0}_{S}) vs p_{T}; m_{inv}(p,K^{0}_{S}) [GeV/c^{2}]; p_{T}(#Lambda_{c}) [GeV/c]";
    TH2F* spectrumLcMassByK0S = new TH2F(nameHisto.Data(),titleHisto.Data(),1000,mLcPDG-0.250,mLcPDG+0.250,11,binLimpTLc);

    nameHisto="histpK0Svsp";
    titleHisto="p(K^{0}_{S}) vs p(p); p(p) [GeV/c]; p(K^{0}_{S}) [GeV/c]";
    TH2F* momentumDistributionK0Svsp = new TH2F(nameHisto.Data(),titleHisto.Data(),41,binLimpTprong,41,binLimpTprong);

    nameHisto="histArmPodK0S";
    titleHisto="V0-candidate Armenteros-Podolanski distribution; #frac{p_{L}^{+}-p_{L}^{-}}{p_{L}^{+}+p_{L}^{-}}; p_{T}^{+} [GeV/c]";
    TH2F* armenterosPodK0S = new TH2F(nameHisto.Data(),titleHisto.Data(),200,-1.,1.,300,0.,0.3);

    nameHisto="histArmPodLc";
    titleHisto="#Lambda_{c}-candidate Armenteros-Podolanski distribution; #frac{p_{L}^{+}-p_{L}^{-}}{p_{L}^{+}+p_{L}^{-}}; p_{T}^{+} [GeV/c]";
    TH2F* armenterosPodLc = new TH2F(nameHisto.Data(),titleHisto.Data(),200,-4.,4.,800,0.,1.6);

    TH2F* allspectrumK0SMass = (TH2F*)spectrumK0SMass->Clone();
    TH2F* allspectrumLcMassByK0S = (TH2F*)spectrumLcMassByK0S->Clone();
    TH2F* allmomentumDistributionK0Svsp = (TH2F*)momentumDistributionK0Svsp->Clone();
    TH2F* allArmenterosPodK0S = (TH2F*)armenterosPodK0S->Clone();
    TH2F* allArmenterosPodLc = (TH2F*)armenterosPodLc->Clone();

    TH2F* pidBachspectrumK0SMass = (TH2F*)spectrumK0SMass->Clone();
    TH2F* pidBachspectrumLcMassByK0S = (TH2F*)spectrumLcMassByK0S->Clone();
    TH2F* pidBachmomentumDistributionK0Svsp = (TH2F*)momentumDistributionK0Svsp->Clone();
    TH2F* pidBachArmenterosPodK0S = (TH2F*)armenterosPodK0S->Clone();
    TH2F* pidBachArmenterosPodLc = (TH2F*)armenterosPodLc->Clone();

    fOutputAll->Add(allspectrumK0SMass);
    fOutputAll->Add(allspectrumLcMassByK0S);
    fOutputAll->Add(allmomentumDistributionK0Svsp);
    fOutputAll->Add(allArmenterosPodK0S);
    fOutputAll->Add(allArmenterosPodLc);

    fOutputPIDBach->Add(pidBachspectrumK0SMass);
    fOutputPIDBach->Add(pidBachspectrumLcMassByK0S);
    fOutputPIDBach->Add(pidBachmomentumDistributionK0Svsp);
    fOutputPIDBach->Add(pidBachArmenterosPodK0S);
    fOutputPIDBach->Add(pidBachArmenterosPodLc);

    nameHisto="histArmPodK0S0";
    titleHisto="V0-candidate Armenteros-Podolanski distribution; #frac{p_{L}^{+}-p_{L}^{-}}{p_{L}^{+}+p_{L}^{-}}; p_{T}^{+} [GeV/c]";
    TH2F* armenterosPodK0S0 = new TH2F(nameHisto.Data(),titleHisto.Data(),200,-1.,1.,300,0.,0.3);
    nameHisto="histArmPodLc0";
    titleHisto="#Lambda_{c}-candidate Armenteros-Podolanski distribution; #frac{p_{L}^{+}-p_{L}^{-}}{p_{L}^{+}+p_{L}^{-}}; p_{T}^{+} [GeV/c]";
    TH2F* armenterosPodLc0 = new TH2F(nameHisto.Data(),titleHisto.Data(),200,-4.,4.,800,0.,1.6);
    fOutputAll->Add(armenterosPodK0S0);
    fOutputAll->Add(armenterosPodLc0);


    if (fTrackRotation) {
      TH2F* pidBachTRspectrumLcMassByK0S = (TH2F*)spectrumLcMassByK0S->Clone();
      fOutputPIDBachTR->Add(pidBachTRspectrumLcMassByK0S);
    }



    nameHisto="histptK0S";
    titleHisto="p_{T}(K^{0}_{S}) vs p_{T}(#Lambda_{c}); p_{T}(#Lambda_{c}) [GeV/c]; p_{T}(K^{0}_{S}) [GeV/c]; Entries";
    TH2F* ptK0S = new TH2F(nameHisto.Data(),titleHisto.Data(),11,binLimpTLc,41,binLimpTprong);

    nameHisto="histptP";
    titleHisto="p_{T}(p) vs p_{T}(#Lambda_{c}); p_{T}(#Lambda_{c}) [GeV/c]; p_{T}(p) [GeV/c]; Entries";
    TH2F* ptP = new TH2F(nameHisto.Data(),titleHisto.Data(),11,binLimpTLc,41,binLimpTprong);

    nameHisto="histptPip";
    titleHisto="p_{T}(#pi^{+}) vs p_{T}(#Lambda_{c}); p_{T}(#Lambda_{c}) [GeV/c]; p_{T}(#pi^{+}) [GeV/c]; Entries";
    TH2F* ptPiP = new TH2F(nameHisto.Data(),titleHisto.Data(),11,binLimpTLc,41,binLimpTprong);

    nameHisto="histptPim";
    titleHisto="p_{T}(#pi^{-}) vs p_{T}(#Lambda_{c}); p_{T}(#Lambda_{c}) [GeV/c]; p_{T}(#pi^{-}) [GeV/c]; Entries";
    TH2F* ptPiM = new TH2F(nameHisto.Data(),titleHisto.Data(),11,binLimpTLc,41,binLimpTprong);

    nameHisto="histLambdaMass";
    titleHisto="m_{inv}(p,#pi^{-}) vs p_{T}(#Lambda_{c}); p_{T}(#Lambda_{c}) [GeV/c]; m_{inv}(p,#pi^{-}) [GeV/c^{2}]; Entries";
    TH2F* massLambda = new TH2F(nameHisto.Data(),titleHisto.Data(),11,binLimpTLc,1000,mMinLambdaPDG,mMinLambdaPDG+0.5);

    nameHisto="histLambdaBarMass";
    titleHisto="m_{inv}(#bar{p},#pi^{+}) vs p_{T}(#Lambda_{c}); p_{T}(#Lambda_{c}) [GeV/c]; m_{inv}(#bar{p},#pi^{+}) [GeV/c^{2}]; Entries";
    TH2F* massLambdaBar = new TH2F(nameHisto.Data(),titleHisto.Data(),11,binLimpTLc,1000,mMinLambdaPDG,mMinLambdaPDG+0.5);

    nameHisto="histGammaMass";
    titleHisto="m_{inv}(e^{+},e^{-}) vs p_{T}(#Lambda_{c}); p_{T}(#Lambda_{c}) [GeV/c]; m_{inv}(e^{+},e^{-}) [GeV/c^{2}]; Entries";
    TH2F* massGamma = new TH2F(nameHisto.Data(),titleHisto.Data(),11,binLimpTLc,100,0.,1.);

    nameHisto="histD0K0S";
    titleHisto="d_{0}(K^{0}_{S}) vs p_{T}(#Lambda_{c}); p_{T}(#Lambda_{c}) [GeV/c]; d_{0}(K^{0}_{S}) [#sigmas]; Entries";
    TH2F* d0K0S = new TH2F(nameHisto.Data(),titleHisto.Data(),11,binLimpTLc,1000,-1.,1.);

    nameHisto="histD0P";
    titleHisto="d_{0}(p) vs p_{T}(#Lambda_{c}); p_{T}(#Lambda_{c}) [GeV/c]; d_{0}(p) [cm]; Entries";
    TH2F* d0P = new TH2F(nameHisto.Data(),titleHisto.Data(),11,binLimpTLc,1000,-1.,1.);

    nameHisto="histCosPAK0S";
    titleHisto="K^{0}_{S} cosine of pointing angle wrt primary vertex vs p_{T}(#Lambda_{c}); p_{T}(#Lambda_{c}) [GeV/c]; cosine; Entries";
    TH2F *cosPAK0S = new TH2F(nameHisto.Data(),titleHisto.Data(),41,binLimpTprong,100,0.99,1.);

    nameHisto="histCosThetaProtonCMS";
    titleHisto="cosien of proton emission angle in Lc rest frame; p_{T}(#Lambda_{c}) [GeV/c]; cosine; Entries";
    TH2F *cosThePr = new TH2F(nameHisto.Data(),titleHisto.Data(),41,binLimpTprong,100,-1.,1.);

    nameHisto="histResignedD0";
    titleHisto="Proton d0 with different sign convention; p_{T}(#Lambda_{c}) [GeV/c]; d0 [cm]; Entries";
    TH2F *resignedD0 = new TH2F(nameHisto.Data(),titleHisto.Data(),41,binLimpTprong,100,-0.1,0.1);

    TH2F* allptK0S = (TH2F*)ptK0S->Clone();
    TH2F* allptP = (TH2F*)ptP->Clone();
    TH2F* allptPiP = (TH2F*)ptPiP->Clone();
    TH2F* allptPiM = (TH2F*)ptPiM->Clone();
    TH2F* allmassLambda = (TH2F*)massLambda->Clone();
    TH2F* allmassLambdaBar = (TH2F*)massLambdaBar->Clone();
    TH2F* allmassGamma = (TH2F*)massGamma->Clone();
    TH2F* alld0K0S = (TH2F*)d0K0S->Clone();
    TH2F* alld0P = (TH2F*)d0P->Clone();
    TH2F* allcosPAK0S = (TH2F*)cosPAK0S->Clone();
    TH2F* allcosThePr = (TH2F*)cosThePr->Clone();
    TH2F* allresignedD0 = (TH2F*)resignedD0->Clone();

    TH2F* pidptK0S = (TH2F*)ptK0S->Clone();
    TH2F* pidptP = (TH2F*)ptP->Clone();
    TH2F* pidptPiP = (TH2F*)ptPiP->Clone();
    TH2F* pidptPiM = (TH2F*)ptPiM->Clone();
    TH2F* pidmassLambda = (TH2F*)massLambda->Clone();
    TH2F* pidmassLambdaBar = (TH2F*)massLambdaBar->Clone();
    TH2F* pidmassGamma = (TH2F*)massGamma->Clone();
    TH2F* pidd0K0S = (TH2F*)d0K0S->Clone();
    TH2F* pidd0P = (TH2F*)d0P->Clone();
    TH2F* pidcosPAK0S = (TH2F*)cosPAK0S->Clone();
    TH2F* pidcosThePr = (TH2F*)cosThePr->Clone();
    TH2F* pidresignedD0 = (TH2F*)resignedD0->Clone();

    fOutputAll->Add(allptK0S);
    fOutputAll->Add(allptP);
    fOutputAll->Add(allptPiP);
    fOutputAll->Add(allptPiM);
    fOutputAll->Add(allmassLambda);
    fOutputAll->Add(allmassLambdaBar);
    fOutputAll->Add(allmassGamma);
    fOutputAll->Add(alld0K0S);
    fOutputAll->Add(alld0P);
    fOutputAll->Add(allcosPAK0S);
    fOutputAll->Add(allcosThePr);
    fOutputAll->Add(allresignedD0);

    fOutputPIDBach->Add(pidptK0S);
    fOutputPIDBach->Add(pidptP);
    fOutputPIDBach->Add(pidptPiP);
    fOutputPIDBach->Add(pidptPiM);
    fOutputPIDBach->Add(pidmassLambda);
    fOutputPIDBach->Add(pidmassLambdaBar);
    fOutputPIDBach->Add(pidmassGamma);
    fOutputPIDBach->Add(pidd0K0S);
    fOutputPIDBach->Add(pidd0P);
    fOutputPIDBach->Add(pidcosPAK0S);
    fOutputPIDBach->Add(pidcosThePr);
    fOutputPIDBach->Add(pidresignedD0);

    if (fTrackRotation) {

      TH2F* pidTRptK0S = (TH2F*)ptK0S->Clone();
      TH2F* pidTRptP = (TH2F*)ptP->Clone();
      TH2F* pidTRptPiP = (TH2F*)ptPiP->Clone();
      TH2F* pidTRptPiM = (TH2F*)ptPiM->Clone();
      TH2F* pidTRmassLambda = (TH2F*)massLambda->Clone();
      TH2F* pidTRmassLambdaBar = (TH2F*)massLambdaBar->Clone();
      TH2F* pidTRmassGamma = (TH2F*)massGamma->Clone();
      TH2F* pidTRcosPAK0S = (TH2F*)cosPAK0S->Clone();
      TH2F* pidTRcosThePr = (TH2F*)cosThePr->Clone();
      TH2F* pidTRresignedD0 = (TH2F*)resignedD0->Clone();
      fOutputPIDBachTR->Add(pidTRptK0S);
      fOutputPIDBachTR->Add(pidTRptP);
      fOutputPIDBachTR->Add(pidTRptPiP);
      fOutputPIDBachTR->Add(pidTRptPiM);
      fOutputPIDBachTR->Add(pidTRmassLambda);
      fOutputPIDBachTR->Add(pidTRmassLambdaBar);
      fOutputPIDBachTR->Add(pidTRmassGamma);
      fOutputPIDBachTR->Add(pidTRcosPAK0S);
      fOutputPIDBachTR->Add(pidTRcosThePr);
      fOutputPIDBachTR->Add(pidTRresignedD0);

    }

  }

  // V0 invariant masses (offline)
  nameHisto="histK0SMassOffline";
  titleHisto="K^{0}_{S} invariant mass VS p_{T}; p_{T}(#Lambda_{c}) [GeV/c]; m_{inv}(#pi^{+},#pi^{-}) [GeV/c^{2}]; Entries";
  TH2F* spectrumK0SMassOffline = new TH2F(nameHisto.Data(),titleHisto.Data(),11,binLimpTLc,1000,mK0SPDG-0.050,mK0SPDG+0.050);

  // Lc invariant masses (x K0S offline)
  nameHisto="histLcMassByK0SOffline";
  titleHisto="#Lambda_{c} invariant mass (by K^{0}_{S}) vs p_{T}; m_{inv}(p,K^{0}_{S}) [GeV/c^{2}]; p_{T}(#Lambda_{c}) [GeV/c]";
  TH2F* spectrumLcMassOfflineByK0S = new TH2F(nameHisto.Data(),titleHisto.Data(),1000,mLcPDG-0.250,mLcPDG+0.250,11,binLimpTLc);

  nameHisto="histpK0SvspOffline";
  titleHisto="p(K^{0}_{S}) vs p(p); p(p) [GeV/c]; p(K^{0}_{S}) [GeV/c]";
  TH2F* momentumDistributionK0SvspOffline = new TH2F(nameHisto.Data(),titleHisto.Data(),41,binLimpTprong,41,binLimpTprong);

  nameHisto="histArmPodK0SOffline";
  titleHisto="V0-candidate Armenteros-Podolanski distribution; #frac{p_{L}^{+}-p_{L}^{-}}{p_{L}^{+}+p_{L}^{-}}; p_{T}^{+} [GeV/c]";
  TH2F* armenterosPodK0SOffline = new TH2F(nameHisto.Data(),titleHisto.Data(),200,-1.,1.,300,0.,0.3);

  nameHisto="histArmPodLcOffline";
  titleHisto="#Lambda_{c}-candidate Armenteros-Podolanski distribution; #frac{p_{L}^{+}-p_{L}^{-}}{p_{L}^{+}+p_{L}^{-}}; p_{T}^{+} [GeV/c]";
  TH2F* armenterosPodLcOffline = new TH2F(nameHisto.Data(),titleHisto.Data(),200,-4.,4.,800,0.,1.6);

  TH2F* allspectrumK0SMassOffline = (TH2F*)spectrumK0SMassOffline->Clone();
  TH2F* allspectrumLcMassOfflineByK0S = (TH2F*)spectrumLcMassOfflineByK0S->Clone();
  TH2F* allmomentumDistributionK0SvspOffline = (TH2F*)momentumDistributionK0SvspOffline->Clone();
  TH2F* allArmenterosPodK0SOffline = (TH2F*)armenterosPodK0SOffline->Clone();
  TH2F* allArmenterosPodLcOffline = (TH2F*)armenterosPodLcOffline->Clone();

  TH2F* pidBachspectrumK0SMassOffline = (TH2F*)spectrumK0SMassOffline->Clone();
  TH2F* pidBachspectrumLcMassOfflineByK0S = (TH2F*)spectrumLcMassOfflineByK0S->Clone();
  TH2F* pidBachmomentumDistributionK0SvspOffline = (TH2F*)momentumDistributionK0SvspOffline->Clone();
  TH2F* pidBachArmenterosPodK0SOffline = (TH2F*)armenterosPodK0SOffline->Clone();
  TH2F* pidBachArmenterosPodLcOffline = (TH2F*)armenterosPodLcOffline->Clone();

  fOutputAll->Add(allspectrumK0SMassOffline);
  fOutputAll->Add(allspectrumLcMassOfflineByK0S);
  fOutputAll->Add(allmomentumDistributionK0SvspOffline);
  fOutputAll->Add(allArmenterosPodK0SOffline);
  fOutputAll->Add(allArmenterosPodLcOffline);

  fOutputPIDBach->Add(pidBachspectrumK0SMassOffline);
  fOutputPIDBach->Add(pidBachspectrumLcMassOfflineByK0S);
  fOutputPIDBach->Add(pidBachmomentumDistributionK0SvspOffline);
  fOutputPIDBach->Add(pidBachArmenterosPodK0SOffline);
  fOutputPIDBach->Add(pidBachArmenterosPodLcOffline);

  nameHisto="histArmPodK0SOffline0";
  titleHisto="V0-candidate Armenteros-Podolanski distribution; #frac{p_{L}^{+}-p_{L}^{-}}{p_{L}^{+}+p_{L}^{-}}; p_{T}^{+} [GeV/c]";
  TH2F* armenterosPodK0SOffline0 = new TH2F(nameHisto.Data(),titleHisto.Data(),200,-1.,1.,300,0.,0.3);
  nameHisto="histArmPodLcOffline0";
  titleHisto="#Lambda_{c}-candidate Armenteros-Podolanski distribution; #frac{p_{L}^{+}-p_{L}^{-}}{p_{L}^{+}+p_{L}^{-}}; p_{T}^{+} [GeV/c]";
  TH2F* armenterosPodLcOffline0 = new TH2F(nameHisto.Data(),titleHisto.Data(),200,-4.,4.,800,0.,1.6);
  fOutputAll->Add(armenterosPodK0SOffline0);
  fOutputAll->Add(armenterosPodLcOffline0);

  if (fTrackRotation) {
    TH2F* pidBachTRspectrumLcMassOfflineByK0S = (TH2F*)spectrumLcMassOfflineByK0S->Clone();
    fOutputPIDBachTR->Add(pidBachTRspectrumLcMassOfflineByK0S);
  }




  nameHisto="histptK0SOffline";
  titleHisto="p_{T}(K^{0}_{S}) vs p_{T}(#Lambda_{c}); p_{T}(#Lambda_{c}) [GeV/c]; p_{T}(K^{0}_{S}) [GeV/c]; Entries";
  TH2F* ptK0SOffline = new TH2F(nameHisto.Data(),titleHisto.Data(),11,binLimpTLc,41,binLimpTprong);

  nameHisto="histptPOffline";
  titleHisto="p_{T}(p) vs p_{T}(#Lambda_{c}); p_{T}(#Lambda_{c}) [GeV/c]; p_{T}(p) [GeV/c]; Entries";
  TH2F* ptPOffline = new TH2F(nameHisto.Data(),titleHisto.Data(),11,binLimpTLc,41,binLimpTprong);

  nameHisto="histptPipOffline";
  titleHisto="p_{T}(#pi^{+}) vs p_{T}(#Lambda_{c}); p_{T}(#Lambda_{c}) [GeV/c]; p_{T}(#pi^{+}) [GeV/c]; Entries";
  TH2F* ptPiPOffline = new TH2F(nameHisto.Data(),titleHisto.Data(),11,binLimpTLc,41,binLimpTprong);

  nameHisto="histptPimOffline";
  titleHisto="p_{T}(#pi^{-}) vs p_{T}(#Lambda_{c}); p_{T}(#Lambda_{c}) [GeV/c]; p_{T}(#pi^{-}) [GeV/c]; Entries";
  TH2F* ptPiMOffline = new TH2F(nameHisto.Data(),titleHisto.Data(),11,binLimpTLc,41,binLimpTprong);

  nameHisto="histLambdaMassOffline";
  titleHisto="m_{inv}(p,#pi^{-}) vs p_{T}(#Lambda_{c}); p_{T}(#Lambda_{c}) [GeV/c]; m_{inv}(p,#pi^{-}) [GeV/c^{2}]; Entries";
  TH2F* massLambdaOffline = new TH2F(nameHisto.Data(),titleHisto.Data(),11,binLimpTLc,1000,mMinLambdaPDG,mMinLambdaPDG+0.5);

  nameHisto="histLambdaBarMassOffline";
  titleHisto="m_{inv}(#bar{p},#pi^{+}) vs p_{T}(#Lambda_{c}); p_{T}(#Lambda_{c}) [GeV/c]; m_{inv}(#bar{p},#pi^{+}) [GeV/c^{2}]; Entries";
  TH2F* massLambdaBarOffline = new TH2F(nameHisto.Data(),titleHisto.Data(),11,binLimpTLc,1000,mMinLambdaPDG,mMinLambdaPDG+0.5);

  nameHisto="histGammaMassOffline";
  titleHisto="m_{inv}(e^{+},e^{-}) vs p_{T}(#Lambda_{c}); p_{T}(#Lambda_{c}) [GeV/c]; m_{inv}(e^{+},e^{-}) [GeV/c^{2}]; Entries";
  TH2F* massGammaOffline = new TH2F(nameHisto.Data(),titleHisto.Data(),11,binLimpTLc,100,0.,1.);

  nameHisto="histD0K0SOffline";
  titleHisto="d_{0}(K^{0}_{S}) vs p_{T}(#Lambda_{c}); p_{T}(#Lambda_{c}) [GeV/c]; d_{0}(K^{0}_{S}) [#sigmas]; Entries";
  TH2F* d0K0SOffline = new TH2F(nameHisto.Data(),titleHisto.Data(),11,binLimpTLc,1000,-1.,1.);

  nameHisto="histD0POffline";
  titleHisto="d_{0}(p) vs p_{T}(#Lambda_{c}); p_{T}(#Lambda_{c}) [GeV/c]; d_{0}(p) [cm]; Entries";
  TH2F* d0POffline = new TH2F(nameHisto.Data(),titleHisto.Data(),11,binLimpTLc,1000,-1.,1.);

  nameHisto="histCosPAK0SOffline";
  titleHisto="K^{0}_{S} cosine of pointing angle wrt primary vertex vs p_{T}(#Lambda_{c}); p_{T}(#Lambda_{c}) [GeV/c]; cosine; Entries";
  TH2F *cosPAK0SOffline = new TH2F(nameHisto.Data(),titleHisto.Data(),41,binLimpTprong,100,0.99,1.);

  nameHisto="histCosThetaProtonCMSOffline";
  titleHisto="cosien of proton emission angle in Lc rest frame; p_{T}(#Lambda_{c}) [GeV/c]; cosine; Entries";
  TH2F *cosThePrOffline = new TH2F(nameHisto.Data(),titleHisto.Data(),41,binLimpTprong,100,-1.,1.);

  nameHisto="histResignedD0Offline";
  titleHisto="Proton d0 with different sign convention; p_{T}(#Lambda_{c}) [GeV/c]; d0 [cm]; Entries";
  TH2F *resignedD0Offline = new TH2F(nameHisto.Data(),titleHisto.Data(),41,binLimpTprong,100,-0.1,0.1);


  TH2F* allptK0SOffline = (TH2F*)ptK0SOffline->Clone();
  TH2F* allptPOffline = (TH2F*)ptPOffline->Clone();
  TH2F* allptPiPOffline = (TH2F*)ptPiPOffline->Clone();
  TH2F* allptPiMOffline = (TH2F*)ptPiMOffline->Clone();
  TH2F* allmassLambdaOffline = (TH2F*)massLambdaOffline->Clone();
  TH2F* allmassLambdaBarOffline = (TH2F*)massLambdaBarOffline->Clone();
  TH2F* allmassGammaOffline = (TH2F*)massGammaOffline->Clone();
  TH2F* alld0K0SOffline = (TH2F*)d0K0SOffline->Clone();
  TH2F* alld0POffline = (TH2F*)d0POffline->Clone();
  TH2F* allcosPAK0SOffline = (TH2F*)cosPAK0SOffline->Clone();
  TH2F* allcosThePrOffline = (TH2F*)cosThePrOffline->Clone();
  TH2F* allresignedD0Offline = (TH2F*)resignedD0Offline->Clone();

  TH2F* pidptK0SOffline = (TH2F*)ptK0SOffline->Clone();
  TH2F* pidptPOffline = (TH2F*)ptPOffline->Clone();
  TH2F* pidptPiPOffline = (TH2F*)ptPiPOffline->Clone();
  TH2F* pidptPiMOffline = (TH2F*)ptPiMOffline->Clone();
  TH2F* pidmassLambdaOffline = (TH2F*)massLambdaOffline->Clone();
  TH2F* pidmassLambdaBarOffline = (TH2F*)massLambdaBarOffline->Clone();
  TH2F* pidmassGammaOffline = (TH2F*)massGammaOffline->Clone();
  TH2F* pidd0K0SOffline = (TH2F*)d0K0SOffline->Clone();
  TH2F* pidd0POffline = (TH2F*)d0POffline->Clone();
  TH2F* pidcosPAK0SOffline = (TH2F*)cosPAK0SOffline->Clone();
  TH2F* pidcosThePrOffline = (TH2F*)cosThePrOffline->Clone();
  TH2F* pidresignedD0Offline = (TH2F*)resignedD0Offline->Clone();

  fOutputAll->Add(allptK0SOffline);
  fOutputAll->Add(allptPOffline);
  fOutputAll->Add(allptPiPOffline);
  fOutputAll->Add(allptPiMOffline);
  fOutputAll->Add(allmassLambdaOffline);
  fOutputAll->Add(allmassLambdaBarOffline);
  fOutputAll->Add(allmassGammaOffline);
  fOutputAll->Add(alld0K0SOffline);
  fOutputAll->Add(alld0POffline);
  fOutputAll->Add(allcosPAK0SOffline);
  fOutputAll->Add(allcosThePrOffline);
  fOutputAll->Add(allresignedD0Offline);

  fOutputPIDBach->Add(pidptK0SOffline);
  fOutputPIDBach->Add(pidptPOffline);
  fOutputPIDBach->Add(pidptPiPOffline);
  fOutputPIDBach->Add(pidptPiMOffline);
  fOutputPIDBach->Add(pidmassLambdaOffline);
  fOutputPIDBach->Add(pidmassLambdaBarOffline);
  fOutputPIDBach->Add(pidmassGammaOffline);
  fOutputPIDBach->Add(pidd0K0SOffline);
  fOutputPIDBach->Add(pidd0POffline);
  fOutputPIDBach->Add(pidcosPAK0SOffline);
  fOutputPIDBach->Add(pidcosThePrOffline);
  fOutputPIDBach->Add(pidresignedD0Offline);

  if (fTrackRotation) {

    TH2F* pidTRptK0SOffline = (TH2F*)ptK0SOffline->Clone();
    TH2F* pidTRptPOffline = (TH2F*)ptPOffline->Clone();
    TH2F* pidTRptPiPOffline = (TH2F*)ptPiPOffline->Clone();
    TH2F* pidTRptPiMOffline = (TH2F*)ptPiMOffline->Clone();
    TH2F* pidTRmassLambdaOffline = (TH2F*)massLambdaOffline->Clone();
    TH2F* pidTRmassLambdaBarOffline = (TH2F*)massLambdaBarOffline->Clone();
    TH2F* pidTRmassGammaOffline = (TH2F*)massGammaOffline->Clone();
    TH2F* pidTRcosPAK0SOffline = (TH2F*)cosPAK0SOffline->Clone();
    TH2F* pidTRcosThePrOffline = (TH2F*)cosThePrOffline->Clone();
    TH2F* pidTRresignedD0Offline = (TH2F*)resignedD0Offline->Clone();
    fOutputPIDBachTR->Add(pidTRptK0SOffline);
    fOutputPIDBachTR->Add(pidTRptPOffline);
    fOutputPIDBachTR->Add(pidTRptPiPOffline);
    fOutputPIDBachTR->Add(pidTRptPiMOffline);
    fOutputPIDBachTR->Add(pidTRmassLambdaOffline);
    fOutputPIDBachTR->Add(pidTRmassLambdaBarOffline);
    fOutputPIDBachTR->Add(pidTRmassGammaOffline);
    fOutputPIDBachTR->Add(pidTRcosPAK0SOffline);
    fOutputPIDBachTR->Add(pidTRcosThePrOffline);
    fOutputPIDBachTR->Add(pidTRresignedD0Offline);

  }





  if (fUseMCInfo) {

    if (fUseOnTheFlyV0) {

      nameHistoSgn="histK0SMassSgn";
      nameHistoBkg="histK0SMassBkg";
      titleHistoSgn="K^{0}_{S} - sgn: invariant mass VS p_{T} - MC; p_{T}(#Lambda_{c}) [GeV/c]; m_{inv}(#pi^{+},#pi^{-}) [GeV/c^{2}]; Entries";
      titleHistoBkg="K^{0}_{S} - bkg: invariant mass VS p_{T} - MC; p_{T}(#Lambda_{c}) [GeV/c]; m_{inv}(#pi^{+},#pi^{-}) [GeV/c^{2}]; Entries";
      TH2F* spectrumK0SMassSgn = new TH2F(nameHistoSgn.Data(),titleHistoSgn.Data(),11,binLimpTLc,1000,mK0SPDG-0.050,mK0SPDG+0.050);
      TH2F* spectrumK0SMassBkg = new TH2F(nameHistoBkg.Data(),titleHistoBkg.Data(),11,binLimpTLc,1000,mK0SPDG-0.050,mK0SPDG+0.050);

      nameHistoSgn="histLcMassByK0SSgn";
      nameHistoBkg="histLcMassByK0SBkg";
      titleHistoSgn="#Lambda_{c} - sgn: invariant mass (by K^{0}_{S}) vs p_{T} - MC; m_{inv}(p,K^{0}_{S}) [GeV/c^{2}]; p_{T}(#Lambda_{c}) [GeV/c]";
      titleHistoBkg="#Lambda_{c} - bkg: invariant mass (by K^{0}_{S}) vs p_{T} - MC; m_{inv}(p,K^{0}_{S}) [GeV/c^{2}]; p_{T}(#Lambda_{c}) [GeV/c]";
      TH2F* spectrumLcMassByK0SSgn = new TH2F(nameHistoSgn.Data(),titleHistoSgn.Data(),1000,mLcPDG-0.250,mLcPDG+0.250,11,binLimpTLc);
      TH2F* spectrumLcMassByK0SBkg = new TH2F(nameHistoBkg.Data(),titleHistoBkg.Data(),1000,mLcPDG-0.250,mLcPDG+0.250,11,binLimpTLc);

      nameHistoSgn="histpK0SvspSgn";
      nameHistoBkg="histpK0SvspBkg";
      titleHistoSgn="#Lambda_{c} - sgn: K^{0}_{S} vs p Total Momentum Distribution - MC; p(p) [GeV/c]; p(K^{0}_{S}) [GeV/c]";
      titleHistoBkg="#Lambda_{c} - bkg: K^{0}_{S} vs p Total Momentum Distribution - MC; p(p) [GeV/c]; p(K^{0}_{S}) [GeV/c]";
      TH2F* momentumDistributionK0SvspSgn = new TH2F(nameHistoSgn.Data(),titleHistoSgn.Data(),41,binLimpTprong,41,binLimpTprong);
      TH2F* momentumDistributionK0SvspBkg = new TH2F(nameHistoBkg.Data(),titleHistoBkg.Data(),41,binLimpTprong,41,binLimpTprong);

      // armenteros-podolanski plots K0S
      nameHistoSgn="histArmPodK0SSgn";
      nameHistoBkg="histArmPodK0SBkg";
      titleHistoSgn="V0-candidate Armenteros-Podolanski distribution (sgn); #frac{p_{L}^{+}-p_{L}^{-}}{p_{L}^{+}+p_{L}^{-}}; p_{T}^{+} [GeV/c]";
      titleHistoBkg="V0-candidate Armenteros-Podolanski distribution (bkg); #frac{p_{L}^{+}-p_{L}^{-}}{p_{L}^{+}+p_{L}^{-}}; p_{T}^{+} [GeV/c]";
      TH2F* armenterosPodK0SSgn = new TH2F(nameHistoSgn.Data(),titleHistoSgn.Data(),200,-1.,1.,300,0.,0.3);
      TH2F* armenterosPodK0SBkg = new TH2F(nameHistoBkg.Data(),titleHistoBkg.Data(),200,-1.,1.,300,0.,0.3);

      nameHistoSgn="histArmPodLcSgn";
      nameHistoBkg="histArmPodLcBkg";
      titleHistoSgn="#Lambda_{c}-candidate Armenteros-Podolanski distribution (sgn); #frac{p_{L}^{+}-p_{L}^{-}}{p_{L}^{+}+p_{L}^{-}}; p_{T}^{+} [GeV/c]";
      titleHistoBkg="#Lambda_{c}-candidate Armenteros-Podolanski distribution (bkg); #frac{p_{L}^{+}-p_{L}^{-}}{p_{L}^{+}+p_{L}^{-}}; p_{T}^{+} [GeV/c]";
      TH2F* armenterosPodLcSgn = new TH2F(nameHistoSgn.Data(),titleHistoSgn.Data(),200,-4.,4.,800,0.,1.6);
      TH2F* armenterosPodLcBkg = new TH2F(nameHistoBkg.Data(),titleHistoBkg.Data(),200,-4.,4.,800,0.,1.6);

      TH2F* allspectrumK0SMassSgn = (TH2F*)spectrumK0SMassSgn->Clone();
      TH2F* allspectrumK0SMassBkg = (TH2F*)spectrumK0SMassBkg->Clone();
      TH2F* allspectrumLcMassByK0SSgn = (TH2F*)spectrumLcMassByK0SSgn->Clone();
      TH2F* allspectrumLcMassByK0SBkg = (TH2F*)spectrumLcMassByK0SBkg->Clone();
      TH2F* allmomentumDistributionK0SvspSgn = (TH2F*)momentumDistributionK0SvspSgn->Clone();
      TH2F* allmomentumDistributionK0SvspBkg = (TH2F*)momentumDistributionK0SvspBkg->Clone();
      TH2F* allArmenterosPodK0SSgn = (TH2F*)armenterosPodK0SSgn->Clone();
      TH2F* allArmenterosPodK0SBkg = (TH2F*)armenterosPodK0SBkg->Clone();
      TH2F* allArmenterosPodLcSgn = (TH2F*)armenterosPodLcSgn->Clone();
      TH2F* allArmenterosPodLcBkg = (TH2F*)armenterosPodLcBkg->Clone();

      TH2F* pidBachspectrumK0SMassSgn = (TH2F*)spectrumK0SMassSgn->Clone();
      TH2F* pidBachspectrumK0SMassBkg = (TH2F*)spectrumK0SMassBkg->Clone();
      TH2F* pidBachspectrumLcMassByK0SSgn = (TH2F*)spectrumLcMassByK0SSgn->Clone();
      TH2F* pidBachspectrumLcMassByK0SBkg = (TH2F*)spectrumLcMassByK0SBkg->Clone();
      TH2F* pidBachmomentumDistributionK0SvspSgn = (TH2F*)momentumDistributionK0SvspSgn->Clone();
      TH2F* pidBachmomentumDistributionK0SvspBkg = (TH2F*)momentumDistributionK0SvspBkg->Clone();
      TH2F* pidBachArmenterosPodK0SSgn = (TH2F*)armenterosPodK0SSgn->Clone();
      TH2F* pidBachArmenterosPodK0SBkg = (TH2F*)armenterosPodK0SBkg->Clone();
      TH2F* pidBachArmenterosPodLcSgn = (TH2F*)armenterosPodLcSgn->Clone();
      TH2F* pidBachArmenterosPodLcBkg = (TH2F*)armenterosPodLcBkg->Clone();

      fOutputAll->Add(allspectrumK0SMassSgn);
      fOutputAll->Add(allspectrumK0SMassBkg);
      fOutputAll->Add(allspectrumLcMassByK0SSgn);
      fOutputAll->Add(allspectrumLcMassByK0SBkg);
      fOutputAll->Add(allmomentumDistributionK0SvspSgn);
      fOutputAll->Add(allmomentumDistributionK0SvspBkg);
      fOutputAll->Add(allArmenterosPodK0SSgn);
      fOutputAll->Add(allArmenterosPodK0SBkg);
      fOutputAll->Add(allArmenterosPodLcSgn);
      fOutputAll->Add(allArmenterosPodLcBkg);

      fOutputPIDBach->Add(pidBachspectrumK0SMassSgn);
      fOutputPIDBach->Add(pidBachspectrumK0SMassBkg);
      fOutputPIDBach->Add(pidBachspectrumLcMassByK0SSgn);
      fOutputPIDBach->Add(pidBachspectrumLcMassByK0SBkg);
      fOutputPIDBach->Add(pidBachmomentumDistributionK0SvspSgn);
      fOutputPIDBach->Add(pidBachmomentumDistributionK0SvspBkg);
      fOutputPIDBach->Add(pidBachArmenterosPodK0SSgn);
      fOutputPIDBach->Add(pidBachArmenterosPodK0SBkg);
      fOutputPIDBach->Add(pidBachArmenterosPodLcSgn);
      fOutputPIDBach->Add(pidBachArmenterosPodLcBkg);

      nameHistoSgn="histArmPodK0SSgn0";
      nameHistoBkg="histArmPodK0SBkg0";
      titleHistoSgn="V0-candidate Armenteros-Podolanski distribution (sgn); #frac{p_{L}^{+}-p_{L}^{-}}{p_{L}^{+}+p_{L}^{-}}; p_{T}^{+} [GeV/c]";
      titleHistoBkg="V0-candidate Armenteros-Podolanski distribution (bkg); #frac{p_{L}^{+}-p_{L}^{-}}{p_{L}^{+}+p_{L}^{-}}; p_{T}^{+} [GeV/c]";
      TH2F* armenterosPodK0SSgn0 = new TH2F(nameHistoSgn.Data(),titleHistoSgn.Data(),200,-1.,1.,300,0.,0.3);
      TH2F* armenterosPodK0SBkg0 = new TH2F(nameHistoBkg.Data(),titleHistoBkg.Data(),200,-1.,1.,300,0.,0.3);
      fOutputAll->Add(armenterosPodK0SSgn0);
      fOutputAll->Add(armenterosPodK0SBkg0);
      nameHistoSgn="histArmPodLcSgn0";
      nameHistoBkg="histArmPodLcBkg0";
      titleHistoSgn="#Lambda_{c}-candidate Armenteros-Podolanski distribution (sgn); #frac{p_{L}^{+}-p_{L}^{-}}{p_{L}^{+}+p_{L}^{-}}; p_{T}^{+} [GeV/c]";
      titleHistoBkg="#Lambda_{c}-candidate Armenteros-Podolanski distribution (bkg); #frac{p_{L}^{+}-p_{L}^{-}}{p_{L}^{+}+p_{L}^{-}}; p_{T}^{+} [GeV/c]";
      TH2F* armenterosPodLcSgn0 = new TH2F(nameHistoSgn.Data(),titleHistoSgn.Data(),200,-4.,4.,800,0.,1.6);
      TH2F* armenterosPodLcBkg0 = new TH2F(nameHistoBkg.Data(),titleHistoBkg.Data(),200,-4.,4.,800,0.,1.6);
      fOutputAll->Add(armenterosPodLcSgn0);
      fOutputAll->Add(armenterosPodLcBkg0);

      if (fTrackRotation) {
	TH2F* pidBachTRspectrumLcMassByK0SSgn = (TH2F*)spectrumLcMassByK0SSgn->Clone();
	TH2F* pidBachTRspectrumLcMassByK0SBkg = (TH2F*)spectrumLcMassByK0SBkg->Clone();
	fOutputPIDBachTR->Add(pidBachTRspectrumLcMassByK0SSgn);
	fOutputPIDBachTR->Add(pidBachTRspectrumLcMassByK0SBkg);
      }



      nameHistoSgn="histptK0SSgn";
      nameHistoBkg="histptK0SBkg";
      titleHistoSgn="p_{T}(K^{0}_{S}) vs p_{T}(#Lambda_{c}); p_{T}(#Lambda_{c}) [GeV/c]; p_{T}(K^{0}_{S}) [GeV/c]; Entries";
      titleHistoBkg="p_{T}(K^{0}_{S}) vs p_{T}(#Lambda_{c}); p_{T}(#Lambda_{c}) [GeV/c]; p_{T}(K^{0}_{S}) [GeV/c]; Entries";
      TH2F* ptK0SSgn = new TH2F(nameHistoSgn.Data(),titleHistoSgn.Data(),11,binLimpTLc,41,binLimpTprong);
      TH2F* ptK0SBkg = new TH2F(nameHistoBkg.Data(),titleHistoBkg.Data(),11,binLimpTLc,41,binLimpTprong);

      nameHistoSgn="histptPSgn";
      nameHistoBkg="histptPBkg";
      titleHistoSgn="p_{T}(p) vs p_{T}(#Lambda_{c}); p_{T}(#Lambda_{c}) [GeV/c]; p_{T}(p) [GeV/c]; Entries";
      titleHistoBkg="p_{T}(p) vs p_{T}(#Lambda_{c}); p_{T}(#Lambda_{c}) [GeV/c]; p_{T}(p) [GeV/c]; Entries";
      TH2F* ptPSgn = new TH2F(nameHistoSgn.Data(),titleHistoSgn.Data(),11,binLimpTLc,41,binLimpTprong);
      TH2F* ptPBkg = new TH2F(nameHistoBkg.Data(),titleHistoBkg.Data(),11,binLimpTLc,41,binLimpTprong);

      nameHistoSgn="histptPipSgn";
      nameHistoBkg="histptPipBkg";
      titleHistoSgn="p_{T}(#pi^{+}) vs p_{T}(#Lambda_{c}); p_{T}(#Lambda_{c}) [GeV/c]; p_{T}(#pi^{+}) [GeV/c]; Entries";
      titleHistoBkg="p_{T}(#pi^{+}) vs p_{T}(#Lambda_{c}); p_{T}(#Lambda_{c}) [GeV/c]; p_{T}(#pi^{+}) [GeV/c]; Entries";
      TH2F* ptPiPSgn = new TH2F(nameHistoSgn.Data(),titleHistoSgn.Data(),11,binLimpTLc,41,binLimpTprong);
      TH2F* ptPiPBkg = new TH2F(nameHistoBkg.Data(),titleHistoBkg.Data(),11,binLimpTLc,41,binLimpTprong);

      nameHistoSgn="histptPimSgn";
      nameHistoBkg="histptPimBkg";
      titleHistoSgn="p_{T}(#pi^{-}) vs p_{T}(#Lambda_{c}); p_{T}(#Lambda_{c}) [GeV/c]; p_{T}(#pi^{-}) [GeV/c]; Entries";
      titleHistoBkg="p_{T}(#pi^{-}) vs p_{T}(#Lambda_{c}); p_{T}(#Lambda_{c}) [GeV/c]; p_{T}(#pi^{-}) [GeV/c]; Entries";
      TH2F* ptPiMSgn = new TH2F(nameHistoSgn.Data(),titleHistoSgn.Data(),11,binLimpTLc,41,binLimpTprong);
      TH2F* ptPiMBkg = new TH2F(nameHistoBkg.Data(),titleHistoBkg.Data(),11,binLimpTLc,41,binLimpTprong);

      nameHistoSgn="histLambdaMassSgn";
      nameHistoBkg="histLambdaMassBkg";
      titleHistoSgn="m_{inv}(p,#pi^{-}) vs p_{T}(#Lambda_{c}); p_{T}(#Lambda_{c}) [GeV/c]; m_{inv}(p,#pi^{-}) [GeV/c^{2}]; Entries";
      titleHistoBkg="m_{inv}(p,#pi^{-}) vs p_{T}(#Lambda_{c}); p_{T}(#Lambda_{c}) [GeV/c]; m_{inv}(p,#pi^{-}) [GeV/c^{2}]; Entries";
      TH2F* massLambdaSgn = new TH2F(nameHistoSgn.Data(),titleHistoSgn.Data(),11,binLimpTLc,1000,mMinLambdaPDG,mMinLambdaPDG+0.5);
      TH2F* massLambdaBkg = new TH2F(nameHistoBkg.Data(),titleHistoBkg.Data(),11,binLimpTLc,1000,mMinLambdaPDG,mMinLambdaPDG+0.5);

      nameHistoSgn="histLambdaBarMassSgn";
      nameHistoBkg="histLambdaBarMassBkg";
      titleHistoSgn="m_{inv}(#bar{p},#pi^{+}) vs p_{T}(#Lambda_{c}); p_{T}(#Lambda_{c}) [GeV/c]; m_{inv}(#bar{p},#pi^{+}) [GeV/c^{2}]; Entries";
      titleHistoBkg="m_{inv}(#bar{p},#pi^{+}) vs p_{T}(#Lambda_{c}); p_{T}(#Lambda_{c}) [GeV/c]; m_{inv}(#bar{p},#pi^{+}) [GeV/c^{2}]; Entries";
      TH2F* massLambdaBarSgn = new TH2F(nameHistoSgn.Data(),titleHistoSgn.Data(),11,binLimpTLc,1000,mMinLambdaPDG,mMinLambdaPDG+0.5);
      TH2F* massLambdaBarBkg = new TH2F(nameHistoBkg.Data(),titleHistoBkg.Data(),11,binLimpTLc,1000,mMinLambdaPDG,mMinLambdaPDG+0.5);

      nameHistoSgn="histGammaMassSgn";
      nameHistoBkg="histGammaMassBkg";
      titleHistoSgn="m_{inv}(e^{+},e^{-}) vs p_{T}(#Lambda_{c}); p_{T}(#Lambda_{c}) [GeV/c]; m_{inv}(e^{+},e^{-}) [GeV/c^{2}]; Entries";
      titleHistoBkg="m_{inv}(e^{+},e^{-}) vs p_{T}(#Lambda_{c}); p_{T}(#Lambda_{c}) [GeV/c]; m_{inv}(e^{+},e^{-}) [GeV/c^{2}]; Entries";
      TH2F* massGammaSgn = new TH2F(nameHistoSgn.Data(),titleHistoSgn.Data(),11,binLimpTLc,100,0.,1.);
      TH2F* massGammaBkg = new TH2F(nameHistoBkg.Data(),titleHistoBkg.Data(),11,binLimpTLc,100,0.,1.);

      nameHistoSgn="histD0K0SSgn";
      nameHistoBkg="histD0K0SBkg";
      titleHistoSgn="d_{0}(K^{0}_{S}) vs p_{T}(#Lambda_{c}); p_{T}(#Lambda_{c}) [GeV/c]; d_{0}(K^{0}_{S}) [#sigmas]; Entries";
      titleHistoBkg="d_{0}(K^{0}_{S}) vs p_{T}(#Lambda_{c}); p_{T}(#Lambda_{c}) [GeV/c]; d_{0}(K^{0}_{S}) [#sigmas]; Entries";
      TH2F* d0K0SSgn = new TH2F(nameHistoSgn.Data(),titleHistoSgn.Data(),11,binLimpTLc,1000,-1.,1.);
      TH2F* d0K0SBkg = new TH2F(nameHistoBkg.Data(),titleHistoBkg.Data(),11,binLimpTLc,1000,-1.,1.);

      nameHistoSgn="histD0PSgn";
      nameHistoBkg="histD0PBkg";
      titleHistoSgn="d_{0}(p) vs p_{T}(#Lambda_{c}); p_{T}(#Lambda_{c}) [GeV/c]; d_{0}(p) [cm]; Entries";
      titleHistoBkg="d_{0}(p) vs p_{T}(#Lambda_{c}); p_{T}(#Lambda_{c}) [GeV/c]; d_{0}(p) [cm]; Entries";
      TH2F* d0PSgn = new TH2F(nameHistoSgn.Data(),titleHistoSgn.Data(),11,binLimpTLc,1000,-1.,1.);
      TH2F* d0PBkg = new TH2F(nameHistoBkg.Data(),titleHistoBkg.Data(),11,binLimpTLc,1000,-1.,1.);

      nameHistoSgn="histCosPAK0SSgn";
      nameHistoBkg="histCosPAK0SBkg";
      titleHistoSgn="K^{0}_{S} cosine of pointing angle wrt primary vertex vs p_{T}(#Lambda_{c}); p_{T}(#Lambda_{c}) [GeV/c]; cosine; Entries";
      titleHistoBkg="K^{0}_{S} cosine of pointing angle wrt primary vertex vs p_{T}(#Lambda_{c}); p_{T}(#Lambda_{c}) [GeV/c]; cosine; Entries";
      TH2F *cosPAK0SSgn = new TH2F(nameHistoSgn.Data(),titleHistoSgn.Data(),41,binLimpTprong,100,0.99,1.);
      TH2F *cosPAK0SBkg = new TH2F(nameHistoBkg.Data(),titleHistoBkg.Data(),41,binLimpTprong,100,0.99,1.);

      nameHistoSgn="histCosThetaProtonCMSSgn";
      nameHistoBkg="histCosThetaProtonCMSBkg";
      titleHistoSgn="cosien of proton emission angle in Lc rest frame; p_{T}(#Lambda_{c}) [GeV/c]; cosine; Entries";
      titleHistoBkg="cosien of proton emission angle in Lc rest frame; p_{T}(#Lambda_{c}) [GeV/c]; cosine; Entries";
      TH2F *cosThePrSgn = new TH2F(nameHistoSgn.Data(),titleHistoSgn.Data(),41,binLimpTprong,100,-1.,1.);
      TH2F *cosThePrBkg = new TH2F(nameHistoBkg.Data(),titleHistoBkg.Data(),41,binLimpTprong,100,-1.,1.);


      nameHistoSgn="histResignedD0Sgn";
      nameHistoBkg="histResignedD0Bkg";
      titleHistoSgn="Proton d0 with different sign convention; p_{T}(#Lambda_{c}) [GeV/c]; d0 [cm]; Entries";
      titleHistoBkg="Proton d0 with different sign convention; p_{T}(#Lambda_{c}) [GeV/c]; d0 [cm]; Entries";
      TH2F *resignedD0Sgn = new TH2F(nameHistoSgn.Data(),titleHistoSgn.Data(),41,binLimpTprong,100,-0.1,0.1);
      TH2F *resignedD0Bkg = new TH2F(nameHistoBkg.Data(),titleHistoBkg.Data(),41,binLimpTprong,100,-0.1,0.1);

      TH2F* allptK0SSgn = (TH2F*)ptK0SSgn->Clone();
      TH2F* allptK0SBkg = (TH2F*)ptK0SBkg->Clone();
      TH2F* allptPSgn = (TH2F*)ptPSgn->Clone();
      TH2F* allptPBkg = (TH2F*)ptPBkg->Clone();
      TH2F* allptPiPSgn = (TH2F*)ptPiPSgn->Clone();
      TH2F* allptPiPBkg = (TH2F*)ptPiPBkg->Clone();
      TH2F* allptPiMSgn = (TH2F*)ptPiMSgn->Clone();
      TH2F* allptPiMBkg = (TH2F*)ptPiMBkg->Clone();
      TH2F* allmassLambdaSgn = (TH2F*)massLambdaSgn->Clone();
      TH2F* allmassLambdaBkg = (TH2F*)massLambdaBkg->Clone();
      TH2F* allmassLambdaBarSgn = (TH2F*)massLambdaBarSgn->Clone();
      TH2F* allmassLambdaBarBkg = (TH2F*)massLambdaBarBkg->Clone();
      TH2F* allmassGammaSgn = (TH2F*)massGammaSgn->Clone();
      TH2F* allmassGammaBkg = (TH2F*)massGammaBkg->Clone();
      TH2F* alld0K0SSgn = (TH2F*)d0K0SSgn->Clone();
      TH2F* alld0K0SBkg = (TH2F*)d0K0SBkg->Clone();
      TH2F* alld0PSgn = (TH2F*)d0PSgn->Clone();
      TH2F* alld0PBkg = (TH2F*)d0PBkg->Clone();
      TH2F* allcosPAK0SSgn = (TH2F*)cosPAK0SSgn->Clone();
      TH2F* allcosPAK0SBkg = (TH2F*)cosPAK0SBkg->Clone();
      TH2F* allcosThePrSgn = (TH2F*)cosThePrSgn->Clone();
      TH2F* allcosThePrBkg = (TH2F*)cosThePrBkg->Clone();
      TH2F* allresignedD0Sgn = (TH2F*)resignedD0Sgn->Clone();
      TH2F* allresignedD0Bkg = (TH2F*)resignedD0Bkg->Clone();

      TH2F* pidptK0SSgn = (TH2F*)ptK0SSgn->Clone();
      TH2F* pidptK0SBkg = (TH2F*)ptK0SBkg->Clone();
      TH2F* pidptPSgn = (TH2F*)ptPSgn->Clone();
      TH2F* pidptPBkg = (TH2F*)ptPBkg->Clone();
      TH2F* pidptPiPSgn = (TH2F*)ptPiPSgn->Clone();
      TH2F* pidptPiPBkg = (TH2F*)ptPiPBkg->Clone();
      TH2F* pidptPiMSgn = (TH2F*)ptPiMSgn->Clone();
      TH2F* pidptPiMBkg = (TH2F*)ptPiMBkg->Clone();
      TH2F* pidmassLambdaSgn = (TH2F*)massLambdaSgn->Clone();
      TH2F* pidmassLambdaBkg = (TH2F*)massLambdaBkg->Clone();
      TH2F* pidmassLambdaBarSgn = (TH2F*)massLambdaBarSgn->Clone();
      TH2F* pidmassLambdaBarBkg = (TH2F*)massLambdaBarBkg->Clone();
      TH2F* pidmassGammaSgn = (TH2F*)massGammaSgn->Clone();
      TH2F* pidmassGammaBkg = (TH2F*)massGammaBkg->Clone();
      TH2F* pidd0K0SSgn = (TH2F*)d0K0SSgn->Clone();
      TH2F* pidd0K0SBkg = (TH2F*)d0K0SBkg->Clone();
      TH2F* pidd0PSgn = (TH2F*)d0PSgn->Clone();
      TH2F* pidd0PBkg = (TH2F*)d0PBkg->Clone();
      TH2F* pidcosPAK0SSgn = (TH2F*)cosPAK0SSgn->Clone();
      TH2F* pidcosPAK0SBkg = (TH2F*)cosPAK0SBkg->Clone();
      TH2F* pidcosThePrSgn = (TH2F*)cosThePrSgn->Clone();
      TH2F* pidcosThePrBkg = (TH2F*)cosThePrBkg->Clone();
      TH2F* pidresignedD0Sgn = (TH2F*)resignedD0Sgn->Clone();
      TH2F* pidresignedD0Bkg = (TH2F*)resignedD0Bkg->Clone();

      fOutputAll->Add(allptK0SSgn);
      fOutputAll->Add(allptK0SBkg);
      fOutputAll->Add(allptPSgn);
      fOutputAll->Add(allptPBkg);
      fOutputAll->Add(allptPiPSgn);
      fOutputAll->Add(allptPiPBkg);
      fOutputAll->Add(allptPiMSgn);
      fOutputAll->Add(allptPiMBkg);
      fOutputAll->Add(allmassLambdaSgn);
      fOutputAll->Add(allmassLambdaBkg);
      fOutputAll->Add(allmassLambdaBarSgn);
      fOutputAll->Add(allmassLambdaBarBkg);
      fOutputAll->Add(allmassGammaSgn);
      fOutputAll->Add(allmassGammaBkg);
      fOutputAll->Add(alld0K0SSgn);
      fOutputAll->Add(alld0K0SBkg);
      fOutputAll->Add(alld0PSgn);
      fOutputAll->Add(alld0PBkg);
      fOutputAll->Add(allcosPAK0SSgn);
      fOutputAll->Add(allcosPAK0SBkg);
      fOutputAll->Add(allcosThePrSgn);
      fOutputAll->Add(allcosThePrBkg);
      fOutputAll->Add(allresignedD0Sgn);
      fOutputAll->Add(allresignedD0Bkg);

      fOutputPIDBach->Add(pidptK0SSgn);
      fOutputPIDBach->Add(pidptK0SBkg);
      fOutputPIDBach->Add(pidptPSgn);
      fOutputPIDBach->Add(pidptPBkg);
      fOutputPIDBach->Add(pidptPiPSgn);
      fOutputPIDBach->Add(pidptPiPBkg);
      fOutputPIDBach->Add(pidptPiMSgn);
      fOutputPIDBach->Add(pidptPiMBkg);
      fOutputPIDBach->Add(pidmassLambdaSgn);
      fOutputPIDBach->Add(pidmassLambdaBkg);
      fOutputPIDBach->Add(pidmassLambdaBarSgn);
      fOutputPIDBach->Add(pidmassLambdaBarBkg);
      fOutputPIDBach->Add(pidmassGammaSgn);
      fOutputPIDBach->Add(pidmassGammaBkg);
      fOutputPIDBach->Add(pidd0K0SSgn);
      fOutputPIDBach->Add(pidd0K0SBkg);
      fOutputPIDBach->Add(pidd0PSgn);
      fOutputPIDBach->Add(pidd0PBkg);
      fOutputPIDBach->Add(pidcosPAK0SSgn);
      fOutputPIDBach->Add(pidcosPAK0SBkg);
      fOutputPIDBach->Add(pidcosThePrSgn);
      fOutputPIDBach->Add(pidcosThePrBkg);
      fOutputPIDBach->Add(pidresignedD0Sgn);
      fOutputPIDBach->Add(pidresignedD0Bkg);

      if (fTrackRotation) {

	TH2F* pidTRptK0SSgn = (TH2F*)ptK0SSgn->Clone();
	TH2F* pidTRptK0SBkg = (TH2F*)ptK0SBkg->Clone();
	TH2F* pidTRptPSgn = (TH2F*)ptPSgn->Clone();
	TH2F* pidTRptPBkg = (TH2F*)ptPBkg->Clone();
	TH2F* pidTRptPiPSgn = (TH2F*)ptPiPSgn->Clone();
	TH2F* pidTRptPiPBkg = (TH2F*)ptPiPBkg->Clone();
	TH2F* pidTRptPiMSgn = (TH2F*)ptPiMSgn->Clone();
	TH2F* pidTRptPiMBkg = (TH2F*)ptPiMBkg->Clone();
	TH2F* pidTRmassLambdaSgn = (TH2F*)massLambdaSgn->Clone();
	TH2F* pidTRmassLambdaBkg = (TH2F*)massLambdaBkg->Clone();
	TH2F* pidTRmassLambdaBarSgn = (TH2F*)massLambdaBarSgn->Clone();
	TH2F* pidTRmassLambdaBarBkg = (TH2F*)massLambdaBarBkg->Clone();
	TH2F* pidTRmassGammaSgn = (TH2F*)massGammaSgn->Clone();
	TH2F* pidTRmassGammaBkg = (TH2F*)massGammaBkg->Clone();
	TH2F* pidTRcosPAK0SSgn = (TH2F*)cosPAK0SSgn->Clone();
	TH2F* pidTRcosPAK0SBkg = (TH2F*)cosPAK0SBkg->Clone();
	TH2F* pidTRcosThePrSgn = (TH2F*)cosThePrSgn->Clone();
	TH2F* pidTRcosThePrBkg = (TH2F*)cosThePrBkg->Clone();
	TH2F* pidTRresignedD0Sgn = (TH2F*)resignedD0Sgn->Clone();
	TH2F* pidTRresignedD0Bkg = (TH2F*)resignedD0Bkg->Clone();
	fOutputPIDBachTR->Add(pidTRptK0SSgn);
	fOutputPIDBachTR->Add(pidTRptK0SBkg);
	fOutputPIDBachTR->Add(pidTRptPSgn);
	fOutputPIDBachTR->Add(pidTRptPBkg);
	fOutputPIDBachTR->Add(pidTRptPiPSgn);
	fOutputPIDBachTR->Add(pidTRptPiPBkg);
	fOutputPIDBachTR->Add(pidTRptPiMSgn);
	fOutputPIDBachTR->Add(pidTRptPiMBkg);
	fOutputPIDBachTR->Add(pidTRmassLambdaSgn);
	fOutputPIDBachTR->Add(pidTRmassLambdaBkg);
	fOutputPIDBachTR->Add(pidTRmassLambdaBarSgn);
	fOutputPIDBachTR->Add(pidTRmassLambdaBarBkg);
	fOutputPIDBachTR->Add(pidTRmassGammaSgn);
	fOutputPIDBachTR->Add(pidTRmassGammaBkg);
	fOutputPIDBachTR->Add(pidTRcosPAK0SSgn);
	fOutputPIDBachTR->Add(pidTRcosPAK0SBkg);
	fOutputPIDBachTR->Add(pidTRcosThePrSgn);
	fOutputPIDBachTR->Add(pidTRcosThePrBkg);
	fOutputPIDBachTR->Add(pidTRresignedD0Sgn);
	fOutputPIDBachTR->Add(pidTRresignedD0Bkg);

      }


    } // useOnTheFly


    nameHistoSgn="histK0SMassOfflineSgn";
    nameHistoBkg="histK0SMassOfflineBkg";
    titleHistoSgn="K^{0}_{S} - sgn: invariant mass VS p_{T} - MC; p_{T}(#Lambda_{c}) [GeV/c]; m_{inv}(#pi^{+},#pi^{-}) [GeV/c^{2}]; Entries";
    titleHistoBkg="K^{0}_{S} - bkg: invariant mass VS p_{T} - MC; p_{T}(#Lambda_{c}) [GeV/c]; m_{inv}(#pi^{+},#pi^{-}) [GeV/c^{2}]; Entries";
    TH2F* spectrumK0SMassOfflineSgn = new TH2F(nameHistoSgn.Data(),titleHistoSgn.Data(),11,binLimpTLc,1000,mK0SPDG-0.050,mK0SPDG+0.050);
    TH2F* spectrumK0SMassOfflineBkg = new TH2F(nameHistoBkg.Data(),titleHistoBkg.Data(),11,binLimpTLc,1000,mK0SPDG-0.050,mK0SPDG+0.050);

    nameHistoSgn="histLcMassByK0SOfflineSgn";
    nameHistoBkg="histLcMassByK0SOfflineBkg";
    titleHistoSgn="#Lambda_{c} - sgn: invariant mass (by K^{0}_{S})  vs p_{T} - MC; m_{inv}(p,K^{0}_{S}) [GeV/c^{2}]; p_{T}(#Lambda_{c}) [GeV/c]";
    titleHistoBkg="#Lambda_{c} - bkg: invariant mass (by K^{0}_{S})  vs p_{T} - MC; m_{inv}(p,K^{0}_{S}) [GeV/c^{2}]; p_{T}(#Lambda_{c}) [GeV/c]";
    TH2F* spectrumLcMassOfflineByK0SSgn = new TH2F(nameHistoSgn.Data(),titleHistoSgn.Data(),1000,mLcPDG-0.250,mLcPDG+0.250,11,binLimpTLc);
    TH2F* spectrumLcMassOfflineByK0SBkg = new TH2F(nameHistoBkg.Data(),titleHistoBkg.Data(),1000,mLcPDG-0.250,mLcPDG+0.250,11,binLimpTLc);

    nameHistoSgn="histpK0SvspOfflineSgn";
    nameHistoBkg="histpK0SvspOfflineBkg";
    titleHistoSgn="#Lambda_{c} - sgn: K^{0}_{S} vs p Total Momentum Distribution - Offline - MC; p(p) [GeV/c]; p(K^{0}_{S}) [GeV/c]";
    titleHistoBkg="#Lambda_{c} - bkg: K^{0}_{S} vs p Total Momentum Distribution - Offline - MC; p(p) [GeV/c]; p(K^{0}_{S}) [GeV/c]";
    TH2F* momentumDistributionK0SvspOfflineSgn = new TH2F(nameHistoSgn.Data(),titleHistoSgn.Data(),41,binLimpTprong,41,binLimpTprong);
    TH2F* momentumDistributionK0SvspOfflineBkg = new TH2F(nameHistoBkg.Data(),titleHistoBkg.Data(),41,binLimpTprong,41,binLimpTprong);

    // armenteros-podolanski plots K0S (offline)
    nameHistoSgn="histArmPodK0SOfflineSgn";
    nameHistoBkg="histArmPodK0SOfflineBkg";
    titleHistoSgn="V0-candidate Armenteros-Podolanski distribution (sgn) -offline-; #frac{p_{L}^{+}-p_{L}^{-}}{p_{L}^{+}+p_{L}^{-}}; p_{T}^{+} [GeV/c]";
    titleHistoBkg="V0-candidate Armenteros-Podolanski distribution (bkg) -offline-; #frac{p_{L}^{+}-p_{L}^{-}}{p_{L}^{+}+p_{L}^{-}}; p_{T}^{+} [GeV/c]";
    TH2F* armenterosPodK0SOfflineSgn = new TH2F(nameHistoSgn.Data(),titleHistoSgn.Data(),200,-1.,1.,300,0.,0.3);
    TH2F* armenterosPodK0SOfflineBkg = new TH2F(nameHistoBkg.Data(),titleHistoBkg.Data(),200,-1.,1.,300,0.,0.3);

    nameHistoSgn="histArmPodLcOfflineSgn";
    nameHistoBkg="histArmPodLcOfflineBkg";
    titleHistoSgn="#Lambda_{c}-candidate Armenteros-Podolanski distribution (sgn) -offline-; #frac{p_{L}^{+}-p_{L}^{-}}{p_{L}^{+}+p_{L}^{-}}; p_{T}^{+} [GeV/c]";
    titleHistoBkg="#Lambda_{c}-candidate Armenteros-Podolanski distribution (bkg) -offline-; #frac{p_{L}^{+}-p_{L}^{-}}{p_{L}^{+}+p_{L}^{-}}; p_{T}^{+} [GeV/c]";
    TH2F* armenterosPodLcOfflineSgn = new TH2F(nameHistoSgn.Data(),titleHistoSgn.Data(),200,-4.,4.,800,0.,1.6);
    TH2F* armenterosPodLcOfflineBkg = new TH2F(nameHistoBkg.Data(),titleHistoBkg.Data(),200,-4.,4.,800,0.,1.6);


    TH2F* allspectrumK0SMassOfflineSgn = (TH2F*)spectrumK0SMassOfflineSgn->Clone();
    TH2F* allspectrumK0SMassOfflineBkg = (TH2F*) spectrumK0SMassOfflineBkg->Clone();
    TH2F* allspectrumLcMassOfflineByK0SSgn = (TH2F*)spectrumLcMassOfflineByK0SSgn->Clone();
    TH2F* allspectrumLcMassOfflineByK0SBkg = (TH2F*) spectrumLcMassOfflineByK0SBkg->Clone();
    TH2F* allmomentumDistributionK0SvspOfflineSgn = (TH2F*)momentumDistributionK0SvspOfflineSgn->Clone();
    TH2F* allmomentumDistributionK0SvspOfflineBkg = (TH2F*)momentumDistributionK0SvspOfflineBkg->Clone();
    TH2F* allArmenterosPodK0SOfflineSgn = (TH2F*)armenterosPodK0SOfflineSgn->Clone();
    TH2F* allArmenterosPodK0SOfflineBkg = (TH2F*)armenterosPodK0SOfflineBkg->Clone();
    TH2F* allArmenterosPodLcOfflineSgn = (TH2F*)armenterosPodLcOfflineSgn->Clone();
    TH2F* allArmenterosPodLcOfflineBkg = (TH2F*)armenterosPodLcOfflineBkg->Clone();

    TH2F* pidBachspectrumLcMassOfflineByK0SSgn = (TH2F*)spectrumLcMassOfflineByK0SSgn->Clone();
    TH2F* pidBachspectrumLcMassOfflineByK0SBkg = (TH2F*) spectrumLcMassOfflineByK0SBkg->Clone();
    TH2F* pidBachspectrumK0SMassOfflineSgn = (TH2F*)spectrumK0SMassOfflineSgn->Clone();
    TH2F* pidBachspectrumK0SMassOfflineBkg = (TH2F*) spectrumK0SMassOfflineBkg->Clone();
    TH2F* pidBachmomentumDistributionK0SvspOfflineSgn = (TH2F*)momentumDistributionK0SvspOfflineSgn->Clone();
    TH2F* pidBachmomentumDistributionK0SvspOfflineBkg = (TH2F*)momentumDistributionK0SvspOfflineBkg->Clone();
    TH2F* pidBachArmenterosPodK0SOfflineSgn = (TH2F*)armenterosPodK0SOfflineSgn->Clone();
    TH2F* pidBachArmenterosPodK0SOfflineBkg = (TH2F*)armenterosPodK0SOfflineBkg->Clone();
    TH2F* pidBachArmenterosPodLcOfflineSgn = (TH2F*)armenterosPodLcOfflineSgn->Clone();
    TH2F* pidBachArmenterosPodLcOfflineBkg = (TH2F*)armenterosPodLcOfflineBkg->Clone();

    fOutputAll->Add(allspectrumK0SMassOfflineSgn);
    fOutputAll->Add(allspectrumK0SMassOfflineBkg);
    fOutputAll->Add(allspectrumLcMassOfflineByK0SSgn);
    fOutputAll->Add(allspectrumLcMassOfflineByK0SBkg);
    fOutputAll->Add(allmomentumDistributionK0SvspOfflineSgn);
    fOutputAll->Add(allmomentumDistributionK0SvspOfflineBkg);
    fOutputAll->Add(allArmenterosPodK0SOfflineSgn);
    fOutputAll->Add(allArmenterosPodK0SOfflineBkg);
    fOutputAll->Add(allArmenterosPodLcOfflineSgn);
    fOutputAll->Add(allArmenterosPodLcOfflineBkg);

    fOutputPIDBach->Add(pidBachspectrumK0SMassOfflineSgn);
    fOutputPIDBach->Add(pidBachspectrumK0SMassOfflineBkg);
    fOutputPIDBach->Add(pidBachspectrumLcMassOfflineByK0SSgn);
    fOutputPIDBach->Add(pidBachspectrumLcMassOfflineByK0SBkg);
    fOutputPIDBach->Add(pidBachmomentumDistributionK0SvspOfflineSgn);
    fOutputPIDBach->Add(pidBachmomentumDistributionK0SvspOfflineBkg);
    fOutputPIDBach->Add(pidBachArmenterosPodK0SOfflineSgn);
    fOutputPIDBach->Add(pidBachArmenterosPodK0SOfflineBkg);
    fOutputPIDBach->Add(pidBachArmenterosPodLcOfflineSgn);
    fOutputPIDBach->Add(pidBachArmenterosPodLcOfflineBkg);

    nameHistoSgn="histArmPodK0SOfflineSgn0";
    nameHistoBkg="histArmPodK0SOfflineBkg0";
    titleHistoSgn="V0-candidate Armenteros-Podolanski distribution (sgn) -offline-; #frac{p_{L}^{+}-p_{L}^{-}}{p_{L}^{+}+p_{L}^{-}}; p_{T}^{+} [GeV/c]";
    titleHistoBkg="V0-candidate Armenteros-Podolanski distribution (bkg) -offline-; #frac{p_{L}^{+}-p_{L}^{-}}{p_{L}^{+}+p_{L}^{-}}; p_{T}^{+} [GeV/c]";
    TH2F* armenterosPodK0SOfflineSgn0 = new TH2F(nameHistoSgn.Data(),titleHistoSgn.Data(),200,-1.,1.,300,0.,0.3);
    TH2F* armenterosPodK0SOfflineBkg0 = new TH2F(nameHistoBkg.Data(),titleHistoBkg.Data(),200,-1.,1.,300,0.,0.3);
    nameHistoSgn="histArmPodLcOfflineSgn0";
    nameHistoBkg="histArmPodLcOfflineBkg0";
    titleHistoSgn="#Lambda_{c}-candidate Armenteros-Podolanski distribution (sgn) -offline-; #frac{p_{L}^{+}-p_{L}^{-}}{p_{L}^{+}+p_{L}^{-}}; p_{T}^{+} [GeV/c]";
    titleHistoBkg="#Lambda_{c}-candidate Armenteros-Podolanski distribution (bkg) -offline-; #frac{p_{L}^{+}-p_{L}^{-}}{p_{L}^{+}+p_{L}^{-}}; p_{T}^{+} [GeV/c]";
    TH2F* armenterosPodLcOfflineSgn0 = new TH2F(nameHistoSgn.Data(),titleHistoSgn.Data(),200,-4.,4.,800,0.,1.6);
    TH2F* armenterosPodLcOfflineBkg0 = new TH2F(nameHistoBkg.Data(),titleHistoBkg.Data(),200,-4.,4.,800,0.,1.6);
    fOutputAll->Add(armenterosPodK0SOfflineSgn0);
    fOutputAll->Add(armenterosPodK0SOfflineBkg0);
    fOutputAll->Add(armenterosPodLcOfflineSgn0);
    fOutputAll->Add(armenterosPodLcOfflineBkg0);

    if (fTrackRotation) {
      TH2F* pidBachTRspectrumLcMassOfflineByK0SSgn = (TH2F*)spectrumLcMassOfflineByK0SSgn->Clone();
      TH2F* pidBachTRspectrumLcMassOfflineByK0SBkg = (TH2F*) spectrumLcMassOfflineByK0SBkg->Clone();
      fOutputPIDBachTR->Add(pidBachTRspectrumLcMassOfflineByK0SSgn);
      fOutputPIDBachTR->Add(pidBachTRspectrumLcMassOfflineByK0SBkg);
    }




    nameHistoSgn="histptK0SOfflineSgn";
    nameHistoBkg="histptK0SOfflineBkg";
    titleHistoSgn="p_{T}(K^{0}_{S}) vs p_{T}(#Lambda_{c}); p_{T}(#Lambda_{c}) [GeV/c]; p_{T}(K^{0}_{S}) [GeV/c]; Entries";
    titleHistoBkg="p_{T}(K^{0}_{S}) vs p_{T}(#Lambda_{c}); p_{T}(#Lambda_{c}) [GeV/c]; p_{T}(K^{0}_{S}) [GeV/c]; Entries";
    TH2F* ptK0SOfflineSgn = new TH2F(nameHistoSgn.Data(),titleHistoSgn.Data(),11,binLimpTLc,41,binLimpTprong);
    TH2F* ptK0SOfflineBkg = new TH2F(nameHistoBkg.Data(),titleHistoBkg.Data(),11,binLimpTLc,41,binLimpTprong);

    nameHistoSgn="histptPOfflineSgn";
    nameHistoBkg="histptPOfflineBkg";
    titleHistoSgn="p_{T}(p) vs p_{T}(#Lambda_{c}); p_{T}(#Lambda_{c}) [GeV/c]; p_{T}(p) [GeV/c]; Entries";
    titleHistoBkg="p_{T}(p) vs p_{T}(#Lambda_{c}); p_{T}(#Lambda_{c}) [GeV/c]; p_{T}(p) [GeV/c]; Entries";
    TH2F* ptPOfflineSgn = new TH2F(nameHistoSgn.Data(),titleHistoSgn.Data(),11,binLimpTLc,41,binLimpTprong);
    TH2F* ptPOfflineBkg = new TH2F(nameHistoBkg.Data(),titleHistoBkg.Data(),11,binLimpTLc,41,binLimpTprong);

    nameHistoSgn="histptPipOfflineSgn";
    nameHistoBkg="histptPipOfflineBkg";
    titleHistoSgn="p_{T}(#pi^{+}) vs p_{T}(#Lambda_{c}); p_{T}(#Lambda_{c}) [GeV/c]; p_{T}(#pi^{+}) [GeV/c]; Entries";
    titleHistoBkg="p_{T}(#pi^{+}) vs p_{T}(#Lambda_{c}); p_{T}(#Lambda_{c}) [GeV/c]; p_{T}(#pi^{+}) [GeV/c]; Entries";
    TH2F* ptPiPOfflineSgn = new TH2F(nameHistoSgn.Data(),titleHistoSgn.Data(),11,binLimpTLc,41,binLimpTprong);
    TH2F* ptPiPOfflineBkg = new TH2F(nameHistoBkg.Data(),titleHistoBkg.Data(),11,binLimpTLc,41,binLimpTprong);

    nameHistoSgn="histptPimOfflineSgn";
    nameHistoBkg="histptPimOfflineBkg";
    titleHistoSgn="p_{T}(#pi^{-}) vs p_{T}(#Lambda_{c}); p_{T}(#Lambda_{c}) [GeV/c]; p_{T}(#pi^{-}) [GeV/c]; Entries";
    titleHistoBkg="p_{T}(#pi^{-}) vs p_{T}(#Lambda_{c}); p_{T}(#Lambda_{c}) [GeV/c]; p_{T}(#pi^{-}) [GeV/c]; Entries";
    TH2F* ptPiMOfflineSgn = new TH2F(nameHistoSgn.Data(),titleHistoSgn.Data(),11,binLimpTLc,41,binLimpTprong);
    TH2F* ptPiMOfflineBkg = new TH2F(nameHistoBkg.Data(),titleHistoBkg.Data(),11,binLimpTLc,41,binLimpTprong);

    nameHistoSgn="histLambdaMassOfflineSgn";
    nameHistoBkg="histLambdaMassOfflineBkg";
    titleHistoSgn="m_{inv}(p,#pi^{-}) vs p_{T}(#Lambda_{c}); p_{T}(#Lambda_{c}) [GeV/c]; m_{inv}(p,#pi^{-}) [GeV/c^{2}]; Entries";
    titleHistoBkg="m_{inv}(p,#pi^{-}) vs p_{T}(#Lambda_{c}); p_{T}(#Lambda_{c}) [GeV/c]; m_{inv}(p,#pi^{-}) [GeV/c^{2}]; Entries";
    TH2F* massLambdaOfflineSgn = new TH2F(nameHistoSgn.Data(),titleHistoSgn.Data(),11,binLimpTLc,1000,mMinLambdaPDG,mMinLambdaPDG+0.5);
    TH2F* massLambdaOfflineBkg = new TH2F(nameHistoBkg.Data(),titleHistoBkg.Data(),11,binLimpTLc,1000,mMinLambdaPDG,mMinLambdaPDG+0.5);

    nameHistoSgn="histLambdaBarMassOfflineSgn";
    nameHistoBkg="histLambdaBarMassOfflineBkg";
    titleHistoSgn="m_{inv}(#bar{p},#pi^{+}) vs p_{T}(#Lambda_{c}); p_{T}(#Lambda_{c}) [GeV/c]; m_{inv}(#bar{p},#pi^{+}) [GeV/c^{2}]; Entries";
    titleHistoBkg="m_{inv}(#bar{p},#pi^{+}) vs p_{T}(#Lambda_{c}); p_{T}(#Lambda_{c}) [GeV/c]; m_{inv}(#bar{p},#pi^{+}) [GeV/c^{2}]; Entries";
    TH2F* massLambdaBarOfflineSgn = new TH2F(nameHistoSgn.Data(),titleHistoSgn.Data(),11,binLimpTLc,1000,mMinLambdaPDG,mMinLambdaPDG+0.5);
    TH2F* massLambdaBarOfflineBkg = new TH2F(nameHistoBkg.Data(),titleHistoBkg.Data(),11,binLimpTLc,1000,mMinLambdaPDG,mMinLambdaPDG+0.5);

    nameHistoSgn="histGammaMassOfflineSgn";
    nameHistoBkg="histGammaMassOfflineBkg";
    titleHistoSgn="m_{inv}(e^{+},e^{-}) vs p_{T}(#Lambda_{c}); p_{T}(#Lambda_{c}) [GeV/c]; m_{inv}(e^{+},e^{-}) [GeV/c^{2}]; Entries";
    titleHistoBkg="m_{inv}(e^{+},e^{-}) vs p_{T}(#Lambda_{c}); p_{T}(#Lambda_{c}) [GeV/c]; m_{inv}(e^{+},e^{-}) [GeV/c^{2}]; Entries";
    TH2F* massGammaOfflineSgn = new TH2F(nameHistoSgn.Data(),titleHistoSgn.Data(),11,binLimpTLc,100,0.,1.);
    TH2F* massGammaOfflineBkg = new TH2F(nameHistoBkg.Data(),titleHistoBkg.Data(),11,binLimpTLc,100,0.,1.);

    nameHistoSgn="histD0K0SOfflineSgn";
    nameHistoBkg="histD0K0SOfflineBkg";
    titleHistoSgn="d_{0}(K^{0}_{S}) vs p_{T}(#Lambda_{c}); p_{T}(#Lambda_{c}) [GeV/c]; d_{0}(K^{0}_{S}) [#sigmas]; Entries";
    titleHistoBkg="d_{0}(K^{0}_{S}) vs p_{T}(#Lambda_{c}); p_{T}(#Lambda_{c}) [GeV/c]; d_{0}(K^{0}_{S}) [#sigmas]; Entries";
    TH2F* d0K0SOfflineSgn = new TH2F(nameHistoSgn.Data(),titleHistoSgn.Data(),11,binLimpTLc,1000,-1.,1.);
    TH2F* d0K0SOfflineBkg = new TH2F(nameHistoBkg.Data(),titleHistoBkg.Data(),11,binLimpTLc,1000,-1.,1.);

    nameHistoSgn="histD0POfflineSgn";
    nameHistoBkg="histD0POfflineBkg";
    titleHistoSgn="d_{0}(p) vs p_{T}(#Lambda_{c}); p_{T}(#Lambda_{c}) [GeV/c]; d_{0}(p) [cm]; Entries";
    titleHistoBkg="d_{0}(p) vs p_{T}(#Lambda_{c}); p_{T}(#Lambda_{c}) [GeV/c]; d_{0}(p) [cm]; Entries";
    TH2F* d0POfflineSgn = new TH2F(nameHistoSgn.Data(),titleHistoSgn.Data(),11,binLimpTLc,1000,-1.,1.);
    TH2F* d0POfflineBkg = new TH2F(nameHistoBkg.Data(),titleHistoBkg.Data(),11,binLimpTLc,1000,-1.,1.);

    nameHistoSgn="histCosPAK0SOfflineSgn";
    nameHistoBkg="histCosPAK0SOfflineBkg";
    titleHistoSgn="K^{0}_{S} cosine of pointing angle wrt primary vertex vs p_{T}(#Lambda_{c}); p_{T}(#Lambda_{c}) [GeV/c]; cosine; Entries";
    titleHistoBkg="K^{0}_{S} cosine of pointing angle wrt primary vertex vs p_{T}(#Lambda_{c}); p_{T}(#Lambda_{c}) [GeV/c]; cosine; Entries";
    TH2F *cosPAK0SOfflineSgn = new TH2F(nameHistoSgn.Data(),titleHistoSgn.Data(),41,binLimpTprong,100,0.99,1.);
    TH2F *cosPAK0SOfflineBkg = new TH2F(nameHistoBkg.Data(),titleHistoBkg.Data(),41,binLimpTprong,100,0.99,1.);

    nameHistoSgn="histCosThetaProtonCMSOfflineSgn";
    nameHistoBkg="histCosThetaProtonCMSOfflineBkg";
    titleHistoSgn="cosien of proton emission angle in Lc rest frame; p_{T}(#Lambda_{c}) [GeV/c]; cosine; Entries";
    titleHistoBkg="cosien of proton emission angle in Lc rest frame; p_{T}(#Lambda_{c}) [GeV/c]; cosine; Entries";
    TH2F *cosThePrOfflineSgn = new TH2F(nameHistoSgn.Data(),titleHistoSgn.Data(),41,binLimpTprong,100,-1.,1.);
    TH2F *cosThePrOfflineBkg = new TH2F(nameHistoBkg.Data(),titleHistoBkg.Data(),41,binLimpTprong,100,-1.,1.);


    nameHistoSgn="histResignedD0OfflineSgn";
    nameHistoBkg="histResignedD0OfflineBkg";
    titleHistoSgn="Proton d0 with different sign convention; p_{T}(#Lambda_{c}) [GeV/c]; d0 [cm]; Entries";
    titleHistoBkg="Proton d0 with different sign convention; p_{T}(#Lambda_{c}) [GeV/c]; d0 [cm]; Entries";
    TH2F *resignedD0OfflineSgn = new TH2F(nameHistoSgn.Data(),titleHistoSgn.Data(),41,binLimpTprong,100,-0.1,0.1);
    TH2F *resignedD0OfflineBkg = new TH2F(nameHistoBkg.Data(),titleHistoBkg.Data(),41,binLimpTprong,100,-0.1,0.1);

    TH2F* allptK0SOfflineSgn = (TH2F*)ptK0SOfflineSgn->Clone();
    TH2F* allptK0SOfflineBkg = (TH2F*)ptK0SOfflineBkg->Clone();
    TH2F* allptPOfflineSgn = (TH2F*)ptPOfflineSgn->Clone();
    TH2F* allptPOfflineBkg = (TH2F*)ptPOfflineBkg->Clone();
    TH2F* allptPiPOfflineSgn = (TH2F*)ptPiPOfflineSgn->Clone();
    TH2F* allptPiPOfflineBkg = (TH2F*)ptPiPOfflineBkg->Clone();
    TH2F* allptPiMOfflineSgn = (TH2F*)ptPiMOfflineSgn->Clone();
    TH2F* allptPiMOfflineBkg = (TH2F*)ptPiMOfflineBkg->Clone();
    TH2F* allmassLambdaOfflineSgn = (TH2F*)massLambdaOfflineSgn->Clone();
    TH2F* allmassLambdaOfflineBkg = (TH2F*)massLambdaOfflineBkg->Clone();
    TH2F* allmassLambdaBarOfflineSgn = (TH2F*)massLambdaBarOfflineSgn->Clone();
    TH2F* allmassLambdaBarOfflineBkg = (TH2F*)massLambdaBarOfflineBkg->Clone();
    TH2F* allmassGammaOfflineSgn = (TH2F*)massGammaOfflineSgn->Clone();
    TH2F* allmassGammaOfflineBkg = (TH2F*)massGammaOfflineBkg->Clone();
    TH2F* alld0K0SOfflineSgn = (TH2F*)d0K0SOfflineSgn->Clone();
    TH2F* alld0K0SOfflineBkg = (TH2F*)d0K0SOfflineBkg->Clone();
    TH2F* alld0POfflineSgn = (TH2F*)d0POfflineSgn->Clone();
    TH2F* alld0POfflineBkg = (TH2F*)d0POfflineBkg->Clone();
    TH2F* allcosPAK0SOfflineSgn = (TH2F*)cosPAK0SOfflineSgn->Clone();
    TH2F* allcosPAK0SOfflineBkg = (TH2F*)cosPAK0SOfflineBkg->Clone();
    TH2F* allcosThePrOfflineSgn = (TH2F*)cosThePrOfflineSgn->Clone();
    TH2F* allcosThePrOfflineBkg = (TH2F*)cosThePrOfflineBkg->Clone();
    TH2F* allresignedD0OfflineSgn = (TH2F*)resignedD0OfflineSgn->Clone();
    TH2F* allresignedD0OfflineBkg = (TH2F*)resignedD0OfflineBkg->Clone();

    TH2F* pidptK0SOfflineSgn = (TH2F*)ptK0SOfflineSgn->Clone();
    TH2F* pidptK0SOfflineBkg = (TH2F*)ptK0SOfflineBkg->Clone();
    TH2F* pidptPOfflineSgn = (TH2F*)ptPOfflineSgn->Clone();
    TH2F* pidptPOfflineBkg = (TH2F*)ptPOfflineBkg->Clone();
    TH2F* pidptPiPOfflineSgn = (TH2F*)ptPiPOfflineSgn->Clone();
    TH2F* pidptPiPOfflineBkg = (TH2F*)ptPiPOfflineBkg->Clone();
    TH2F* pidptPiMOfflineSgn = (TH2F*)ptPiMOfflineSgn->Clone();
    TH2F* pidptPiMOfflineBkg = (TH2F*)ptPiMOfflineBkg->Clone();
    TH2F* pidmassLambdaOfflineSgn = (TH2F*)massLambdaOfflineSgn->Clone();
    TH2F* pidmassLambdaOfflineBkg = (TH2F*)massLambdaOfflineBkg->Clone();
    TH2F* pidmassLambdaBarOfflineSgn = (TH2F*)massLambdaBarOfflineSgn->Clone();
    TH2F* pidmassLambdaBarOfflineBkg = (TH2F*)massLambdaBarOfflineBkg->Clone();
    TH2F* pidmassGammaOfflineSgn = (TH2F*)massGammaOfflineSgn->Clone();
    TH2F* pidmassGammaOfflineBkg = (TH2F*)massGammaOfflineBkg->Clone();
    TH2F* pidd0K0SOfflineSgn = (TH2F*)d0K0SOfflineSgn->Clone();
    TH2F* pidd0K0SOfflineBkg = (TH2F*)d0K0SOfflineBkg->Clone();
    TH2F* pidd0POfflineSgn = (TH2F*)d0POfflineSgn->Clone();
    TH2F* pidd0POfflineBkg = (TH2F*)d0POfflineBkg->Clone();
    TH2F* pidcosPAK0SOfflineSgn = (TH2F*)cosPAK0SOfflineSgn->Clone();
    TH2F* pidcosPAK0SOfflineBkg = (TH2F*)cosPAK0SOfflineBkg->Clone();
    TH2F* pidcosThePrOfflineSgn = (TH2F*)cosThePrOfflineSgn->Clone();
    TH2F* pidcosThePrOfflineBkg = (TH2F*)cosThePrOfflineBkg->Clone();
    TH2F* pidresignedD0OfflineSgn = (TH2F*)resignedD0OfflineSgn->Clone();
    TH2F* pidresignedD0OfflineBkg = (TH2F*)resignedD0OfflineBkg->Clone();

    fOutputAll->Add(allptK0SOfflineSgn);
    fOutputAll->Add(allptK0SOfflineBkg);
    fOutputAll->Add(allptPOfflineSgn);
    fOutputAll->Add(allptPOfflineBkg);
    fOutputAll->Add(allptPiPOfflineSgn);
    fOutputAll->Add(allptPiPOfflineBkg);
    fOutputAll->Add(allptPiMOfflineSgn);
    fOutputAll->Add(allptPiMOfflineBkg);
    fOutputAll->Add(allmassLambdaOfflineSgn);
    fOutputAll->Add(allmassLambdaOfflineBkg);
    fOutputAll->Add(allmassLambdaBarOfflineSgn);
    fOutputAll->Add(allmassLambdaBarOfflineBkg);
    fOutputAll->Add(allmassGammaOfflineSgn);
    fOutputAll->Add(allmassGammaOfflineBkg);
    fOutputAll->Add(alld0K0SOfflineSgn);
    fOutputAll->Add(alld0K0SOfflineBkg);
    fOutputAll->Add(alld0POfflineSgn);
    fOutputAll->Add(alld0POfflineBkg);
    fOutputAll->Add(allcosPAK0SOfflineSgn);
    fOutputAll->Add(allcosPAK0SOfflineBkg);
    fOutputAll->Add(allcosThePrOfflineSgn);
    fOutputAll->Add(allcosThePrOfflineBkg);
    fOutputAll->Add(allresignedD0OfflineSgn);
    fOutputAll->Add(allresignedD0OfflineBkg);

    fOutputPIDBach->Add(pidptK0SOfflineSgn);
    fOutputPIDBach->Add(pidptK0SOfflineBkg);
    fOutputPIDBach->Add(pidptPOfflineSgn);
    fOutputPIDBach->Add(pidptPOfflineBkg);
    fOutputPIDBach->Add(pidptPiPOfflineSgn);
    fOutputPIDBach->Add(pidptPiPOfflineBkg);
    fOutputPIDBach->Add(pidptPiMOfflineSgn);
    fOutputPIDBach->Add(pidptPiMOfflineBkg);
    fOutputPIDBach->Add(pidmassLambdaOfflineSgn);
    fOutputPIDBach->Add(pidmassLambdaOfflineBkg);
    fOutputPIDBach->Add(pidmassLambdaBarOfflineSgn);
    fOutputPIDBach->Add(pidmassLambdaBarOfflineBkg);
    fOutputPIDBach->Add(pidmassGammaOfflineSgn);
    fOutputPIDBach->Add(pidmassGammaOfflineBkg);
    fOutputPIDBach->Add(pidd0K0SOfflineSgn);
    fOutputPIDBach->Add(pidd0K0SOfflineBkg);
    fOutputPIDBach->Add(pidd0POfflineSgn);
    fOutputPIDBach->Add(pidd0POfflineBkg);
    fOutputPIDBach->Add(pidcosPAK0SOfflineSgn);
    fOutputPIDBach->Add(pidcosPAK0SOfflineBkg);
    fOutputPIDBach->Add(pidcosThePrOfflineSgn);
    fOutputPIDBach->Add(pidcosThePrOfflineBkg);
    fOutputPIDBach->Add(pidresignedD0OfflineSgn);
    fOutputPIDBach->Add(pidresignedD0OfflineBkg);

    if (fTrackRotation) {

      TH2F* pidTRptK0SOfflineSgn = (TH2F*)ptK0SOfflineSgn->Clone();
      TH2F* pidTRptK0SOfflineBkg = (TH2F*)ptK0SOfflineBkg->Clone();
      TH2F* pidTRptPOfflineSgn = (TH2F*)ptPOfflineSgn->Clone();
      TH2F* pidTRptPOfflineBkg = (TH2F*)ptPOfflineBkg->Clone();
      TH2F* pidTRptPiPOfflineSgn = (TH2F*)ptPiPOfflineSgn->Clone();
      TH2F* pidTRptPiPOfflineBkg = (TH2F*)ptPiPOfflineBkg->Clone();
      TH2F* pidTRptPiMOfflineSgn = (TH2F*)ptPiMOfflineSgn->Clone();
      TH2F* pidTRptPiMOfflineBkg = (TH2F*)ptPiMOfflineBkg->Clone();
      TH2F* pidTRmassLambdaOfflineSgn = (TH2F*)massLambdaOfflineSgn->Clone();
      TH2F* pidTRmassLambdaOfflineBkg = (TH2F*)massLambdaOfflineBkg->Clone();
      TH2F* pidTRmassLambdaBarOfflineSgn = (TH2F*)massLambdaBarOfflineSgn->Clone();
      TH2F* pidTRmassLambdaBarOfflineBkg = (TH2F*)massLambdaBarOfflineBkg->Clone();
      TH2F* pidTRmassGammaOfflineSgn = (TH2F*)massGammaOfflineSgn->Clone();
      TH2F* pidTRmassGammaOfflineBkg = (TH2F*)massGammaOfflineBkg->Clone();
      TH2F* pidTRcosPAK0SOfflineSgn = (TH2F*)cosPAK0SOfflineSgn->Clone();
      TH2F* pidTRcosPAK0SOfflineBkg = (TH2F*)cosPAK0SOfflineBkg->Clone();
      TH2F* pidTRcosThePrOfflineSgn = (TH2F*)cosThePrOfflineSgn->Clone();
      TH2F* pidTRcosThePrOfflineBkg = (TH2F*)cosThePrOfflineBkg->Clone();
      TH2F* pidTRresignedD0OfflineSgn = (TH2F*)resignedD0OfflineSgn->Clone();
      TH2F* pidTRresignedD0OfflineBkg = (TH2F*)resignedD0OfflineBkg->Clone();
      fOutputPIDBachTR->Add(pidTRptK0SOfflineSgn);
      fOutputPIDBachTR->Add(pidTRptK0SOfflineBkg);
      fOutputPIDBachTR->Add(pidTRptPOfflineSgn);
      fOutputPIDBachTR->Add(pidTRptPOfflineBkg);
      fOutputPIDBachTR->Add(pidTRptPiPOfflineSgn);
      fOutputPIDBachTR->Add(pidTRptPiPOfflineBkg);
      fOutputPIDBachTR->Add(pidTRptPiMOfflineSgn);
      fOutputPIDBachTR->Add(pidTRptPiMOfflineBkg);
      fOutputPIDBachTR->Add(pidTRmassLambdaOfflineSgn);
      fOutputPIDBachTR->Add(pidTRmassLambdaOfflineBkg);
      fOutputPIDBachTR->Add(pidTRmassLambdaBarOfflineSgn);
      fOutputPIDBachTR->Add(pidTRmassLambdaBarOfflineBkg);
      fOutputPIDBachTR->Add(pidTRmassGammaOfflineSgn);
      fOutputPIDBachTR->Add(pidTRmassGammaOfflineBkg);
      fOutputPIDBachTR->Add(pidTRcosPAK0SOfflineSgn);
      fOutputPIDBachTR->Add(pidTRcosPAK0SOfflineBkg);
      fOutputPIDBachTR->Add(pidTRcosThePrOfflineSgn);
      fOutputPIDBachTR->Add(pidTRcosThePrOfflineBkg);
      fOutputPIDBachTR->Add(pidTRresignedD0OfflineSgn);
      fOutputPIDBachTR->Add(pidTRresignedD0OfflineBkg);

    }

  } // useMCinfo


  if (fTrackRotation) {

    TH3F *phiVSthetaVSpt = new TH3F("phiVSthetaVSpt","",35,0.,35.,360,0.,2.*TMath::Pi(),100,40.*TMath::DegToRad(),140.*TMath::DegToRad());
    TH3F *phiVSthetaVSptRot = new TH3F("phiVSthetaVSptRot","",35,0.,35.,360,0.,2.*TMath::Pi(),100,40.*TMath::DegToRad(),140.*TMath::DegToRad());
    TH3F *phiVSthetaVSptOffline = new TH3F("phiVSthetaVSptOffline","",35,0.,35.,360,0.,2.*TMath::Pi(),100,40.*TMath::DegToRad(),140.*TMath::DegToRad());
    TH3F *phiVSthetaVSptRotOffline = new TH3F("phiVSthetaVSptRotOffline","",35,0.,35.,360,0.,2.*TMath::Pi(),100,40.*TMath::DegToRad(),140.*TMath::DegToRad());
    fOutputPIDBachTR->Add(phiVSthetaVSpt);
    fOutputPIDBachTR->Add(phiVSthetaVSptRot);
    fOutputPIDBachTR->Add(phiVSthetaVSptOffline);
    fOutputPIDBachTR->Add(phiVSthetaVSptRotOffline);

    TH1F *hNormRotated=new TH1F("hNormRotated","",fNRotations+1,-0.5,fNRotations+0.5);
    TH1F *hNormRotatedOffline=new TH1F("hNormRotatedOffline","",fNRotations+1,-0.5,fNRotations+0.5);
    /*
      hNormRotated->Sumw2();
      hNormRotatedOffline->Sumw2();

      hNormRotated->SetMinimum(0);
      hNormRotatedOffline->SetMinimum(0);
    */

    fOutputPIDBachTR->Add(hNormRotated);
    fOutputPIDBachTR->Add(hNormRotatedOffline);

    if (fUseMCInfo) {

      TH3F *phiVSthetaVSptSgn = new TH3F("phiVSthetaVSptSgn","",35,0.,35.,360,0.,2.*TMath::Pi(),100,40.*TMath::DegToRad(),140.*TMath::DegToRad());
      TH3F *phiVSthetaVSptRotSgn = new TH3F("phiVSthetaVSptRotSgn","",35,0.,35.,360,0.,2.*TMath::Pi(),100,40.*TMath::DegToRad(),140.*TMath::DegToRad());
      TH3F *phiVSthetaVSptOfflineSgn = new TH3F("phiVSthetaVSptOfflineSgn","",35,0.,35.,360,0.,2.*TMath::Pi(),100,40.*TMath::DegToRad(),140.*TMath::DegToRad());
      TH3F *phiVSthetaVSptRotOfflineSgn = new TH3F("phiVSthetaVSptRotOfflineSgn","",35,0.,35.,360,0.,2.*TMath::Pi(),100,40.*TMath::DegToRad(),140.*TMath::DegToRad());
      fOutputPIDBachTR->Add(phiVSthetaVSptSgn);
      fOutputPIDBachTR->Add(phiVSthetaVSptRotSgn);
      fOutputPIDBachTR->Add(phiVSthetaVSptOfflineSgn);
      fOutputPIDBachTR->Add(phiVSthetaVSptRotOfflineSgn);

      TH3F *phiVSthetaVSptBkg = new TH3F("phiVSthetaVSptBkg","",35,0.,35.,360,0.,2.*TMath::Pi(),100,40.*TMath::DegToRad(),140.*TMath::DegToRad());
      TH3F *phiVSthetaVSptRotBkg = new TH3F("phiVSthetaVSptRotBkg","",35,0.,35.,360,0.,2.*TMath::Pi(),100,40.*TMath::DegToRad(),140.*TMath::DegToRad());
      TH3F *phiVSthetaVSptOfflineBkg = new TH3F("phiVSthetaVSptOfflineBkg","",35,0.,35.,360,0.,2.*TMath::Pi(),100,40.*TMath::DegToRad(),140.*TMath::DegToRad());
      TH3F *phiVSthetaVSptRotOfflineBkg = new TH3F("phiVSthetaVSptRotOfflineBkg","",35,0.,35.,360,0.,2.*TMath::Pi(),100,40.*TMath::DegToRad(),140.*TMath::DegToRad());
      fOutputPIDBachTR->Add(phiVSthetaVSptBkg);
      fOutputPIDBachTR->Add(phiVSthetaVSptRotBkg);
      fOutputPIDBachTR->Add(phiVSthetaVSptOfflineBkg);
      fOutputPIDBachTR->Add(phiVSthetaVSptRotOfflineBkg);

      TH1F *hNormRotatedSgn=new TH1F("hNormRotatedSgn","",fNRotations+1,-0.5,fNRotations+0.5);
      TH1F *hNormRotatedOfflineSgn=new TH1F("hNormRotatedOfflineSgn","",fNRotations+1,-0.5,fNRotations+0.5);
      TH1F *hNormRotatedBkg=new TH1F("hNormRotatedBkg","",fNRotations+1,-0.5,fNRotations+0.5);
      TH1F *hNormRotatedOfflineBkg=new TH1F("hNormRotatedOfflineBkg","",fNRotations+1,-0.5,fNRotations+0.5);
      /*
	hNormRotatedSgn->Sumw2();
	hNormRotatedOfflineSgn->Sumw2();
	hNormRotatedBkg->Sumw2();
	hNormRotatedOfflineBkg->Sumw2();

	hNormRotatedSgn->SetMinimum(0);
	hNormRotatedOfflineSgn->SetMinimum(0);
	hNormRotatedBkg->SetMinimum(0);
	hNormRotatedOfflineBkg->SetMinimum(0);
      */

      fOutputPIDBachTR->Add(hNormRotatedSgn);
      fOutputPIDBachTR->Add(hNormRotatedOfflineSgn);
      fOutputPIDBachTR->Add(hNormRotatedBkg);
      fOutputPIDBachTR->Add(hNormRotatedOfflineBkg);

    }

    Int_t nMassBins=(Int_t)(fMaxMass*1000.-fMinMass*1000.);
    Double_t maxm=fMinMass+nMassBins*0.001;
    TH3F *hMassVsPtVsY=new TH3F("hMassVsPtVsY","",nMassBins,fMinMass,maxm,15,0.,15.,20,-1.,1.);
    TH3F *hMassVsPtVsYOffline=new TH3F("hMassVsPtVsYOffline","",nMassBins,fMinMass,maxm,15,0.,15.,20,-1.,1.);
    /*
      hMassVsPtVsY->Sumw2();
      hMassVsPtVsYOffline->Sumw2();

      hMassVsPtVsY->SetMinimum(0);
      hMassVsPtVsYOffline->SetMinimum(0);
    */

    fOutputPIDBachTR->Add(hMassVsPtVsY);
    fOutputPIDBachTR->Add(hMassVsPtVsYOffline);

    if (fUseMCInfo) {

      TH3F *hMassVsPtVsYSgn=new TH3F("hMassVsPtVsYSgn","",nMassBins,fMinMass,maxm,15,0.,15.,20,-1.,1.);
      TH3F *hMassVsPtVsYOfflineSgn=new TH3F("hMassVsPtVsYOfflineSgn","",nMassBins,fMinMass,maxm,15,0.,15.,20,-1.,1.);
      TH3F *hMassVsPtVsYBkg=new TH3F("hMassVsPtVsYBkg","",nMassBins,fMinMass,maxm,15,0.,15.,20,-1.,1.);
      TH3F *hMassVsPtVsYOfflineBkg=new TH3F("hMassVsPtVsYOfflineBkg","",nMassBins,fMinMass,maxm,15,0.,15.,20,-1.,1.);

      /*
	hMassVsPtVsYSgn->Sumw2();
	hMassVsPtVsYOfflineSgn->Sumw2();
	hMassVsPtVsYBkg->Sumw2();
	hMassVsPtVsYOfflineBkg->Sumw2();

	hMassVsPtVsYSgn->SetMinimum(0);
	hMassVsPtVsYOfflineSgn->SetMinimum(0);
	hMassVsPtVsYBkg->SetMinimum(0);
	hMassVsPtVsYOfflineBkg->SetMinimum(0);
      */

      fOutputPIDBachTR->Add(hMassVsPtVsYSgn);
      fOutputPIDBachTR->Add(hMassVsPtVsYOfflineSgn);
      fOutputPIDBachTR->Add(hMassVsPtVsYBkg);
      fOutputPIDBachTR->Add(hMassVsPtVsYOfflineBkg);

    }

    TH3F *hMassVsPtVsYRot=new TH3F("hMassVsPtVsYRot","",nMassBins,fMinMass,maxm,15,0.,15.,20,-1.,1.);
    TH3F *hMassVsPtVsYRotOffline=new TH3F("hMassVsPtVsYRotOffline","",nMassBins,fMinMass,maxm,15,0.,15.,20,-1.,1.);
    /*
      hMassVsPtVsYRot->Sumw2();
      hMassVsPtVsYRotOffline->Sumw2();

      hMassVsPtVsYRot->SetMinimum(0);
      hMassVsPtVsYRotOffline->SetMinimum(0);
    */

    fOutputPIDBachTR->Add(hMassVsPtVsYRot);
    fOutputPIDBachTR->Add(hMassVsPtVsYRotOffline);

    if (fUseMCInfo) {

      TH3F *hMassVsPtVsYRotSgn=new TH3F("hMassVsPtVsYRotSgn","",nMassBins,fMinMass,maxm,15,0.,15.,20,-1.,1.);
      TH3F *hMassVsPtVsYRotOfflineSgn=new TH3F("hMassVsPtVsYRotOfflineSgn","",nMassBins,fMinMass,maxm,15,0.,15.,20,-1.,1.);
      TH3F *hMassVsPtVsYRotBkg=new TH3F("hMassVsPtVsYRotBkg","",nMassBins,fMinMass,maxm,15,0.,15.,20,-1.,1.);
      TH3F *hMassVsPtVsYRotOfflineBkg=new TH3F("hMassVsPtVsYRotOfflineBkg","",nMassBins,fMinMass,maxm,15,0.,15.,20,-1.,1.);
      /*
	hMassVsPtVsYRotSgn->Sumw2();
	hMassVsPtVsYRotOfflineSgn->Sumw2();
	hMassVsPtVsYRotBkg->Sumw2();
	hMassVsPtVsYRotOfflineBkg->Sumw2();

	hMassVsPtVsYRotSgn->SetMinimum(0);
	hMassVsPtVsYRotOfflineSgn->SetMinimum(0);
	hMassVsPtVsYRotBkg->SetMinimum(0);
	hMassVsPtVsYRotOfflineBkg->SetMinimum(0);
      */

      fOutputPIDBachTR->Add(hMassVsPtVsYRotSgn);
      fOutputPIDBachTR->Add(hMassVsPtVsYRotOfflineSgn);
      fOutputPIDBachTR->Add(hMassVsPtVsYRotBkg);
      fOutputPIDBachTR->Add(hMassVsPtVsYRotOfflineBkg);

    }

    TH1F *hDeltaMass=new TH1F("hDeltaMass","",100,-0.4,0.4);
    TH1F *hDeltaMassOffline=new TH1F("hDeltaMassOffline","",100,-0.4,0.4);
    /*
      hDeltaMass->Sumw2();
      hDeltaMassOffline->Sumw2();

      hDeltaMass->SetMinimum(0);
      hDeltaMassOffline->SetMinimum(0);
    */

    fOutputPIDBachTR->Add(hDeltaMass);
    fOutputPIDBachTR->Add(hDeltaMassOffline);

    if (fUseMCInfo) {

      TH1F *hDeltaMassSgn=new TH1F("hDeltaMassSgn","",100,-0.4,0.4);
      TH1F *hDeltaMassOfflineSgn=new TH1F("hDeltaMassOfflineSgn","",100,-0.4,0.4);
      TH1F *hDeltaMassBkg=new TH1F("hDeltaMassBkg","",100,-0.4,0.4);
      TH1F *hDeltaMassOfflineBkg=new TH1F("hDeltaMassOfflineBkg","",100,-0.4,0.4);
      /*
	hDeltaMassSgn->Sumw2();
	hDeltaMassOfflineSgn->Sumw2();
	hDeltaMassBkg->Sumw2();
	hDeltaMassOfflineBkg->Sumw2();

	hDeltaMassSgn->SetMinimum(0);
	hDeltaMassOfflineSgn->SetMinimum(0);
	hDeltaMassBkg->SetMinimum(0);
	hDeltaMassOfflineBkg->SetMinimum(0);
      */

      fOutputPIDBachTR->Add(hDeltaMassSgn);
      fOutputPIDBachTR->Add(hDeltaMassOfflineSgn);
      fOutputPIDBachTR->Add(hDeltaMassBkg);
      fOutputPIDBachTR->Add(hDeltaMassOfflineBkg);

    }

    /*
      Int_t binSparseDMassRot[5]={nMassBins,100,24,40,20};
      Double_t edgeLowSparseDMassRot[5]={fMinMass,-0.4,0.,-4.,0};
      Double_t edgeHighSparseDMassRot[5]={maxm,0.4,12.,4.,3.14};
      THnSparse *hDeltaMassFullAnalysis=new THnSparseF("hDeltaMassFullAnalysis","hDeltaMassFullAnalysis;inv mass (GeV/c);#Delta inv mass (GeV/c); p_{T}^{#Lambda_{c}} (GeV/c); #Delta p_{T} (GeV/c); daughter angle (2prongs) (rad);",5,binSparseDMassRot,edgeLowSparseDMassRot,edgeHighSparseDMassRot);
      THnSparse *hDeltaMassFullAnalysisOffline=new THnSparseF("fDeltaMassFullAnalysisOffline","hDeltaMassFullAnalysisOffline;inv mass (GeV/c);#Delta inv mass (GeV/c); p_{T}^{#Lambda_{c}} (GeV/c); #Delta p_{T} (GeV/c); daughter angle (2prongs) (rad);",5,binSparseDMassRot,edgeLowSparseDMassRot,edgeHighSparseDMassRot);

      fOutputPIDBachTR->Add(hDeltaMassFullAnalysis);
      fOutputPIDBachTR->Add(hDeltaMassFullAnalysisOffline);

      if (fUseMCInfo) {

      THnSparse *hDeltaMassFullAnalysisSgn=new THnSparseF("hDeltaMassFullAnalysisSgn","hDeltaMassFullAnalysisSgn;inv mass (GeV/c);#Delta inv mass (GeV/c); p_{T}^{#Lambda_{c}} (GeV/c); #Delta p_{T} (GeV/c); daughter angle (2prongs) (rad);",5,binSparseDMassRot,edgeLowSparseDMassRot,edgeHighSparseDMassRot);
      THnSparse *hDeltaMassFullAnalysisOfflineSgn=new THnSparseF("fDeltaMassFullAnalysisOfflineSgn","hDeltaMassFullAnalysisOfflineSgn;inv mass (GeV/c);#Delta inv mass (GeV/c); p_{T}^{#Lambda_{c}} (GeV/c); #Delta p_{T} (GeV/c); daughter angle (2prongs) (rad);",5,binSparseDMassRot,edgeLowSparseDMassRot,edgeHighSparseDMassRot);

      THnSparse *hDeltaMassFullAnalysisBkg=new THnSparseF("hDeltaMassFullAnalysisBkg","hDeltaMassFullAnalysisBkg;inv mass (GeV/c);#Delta inv mass (GeV/c); p_{T}^{#Lambda_{c}} (GeV/c); #Delta p_{T} (GeV/c); daughter angle (2prongs) (rad);",5,binSparseDMassRot,edgeLowSparseDMassRot,edgeHighSparseDMassRot);
      THnSparse *hDeltaMassFullAnalysisOfflineBkg=new THnSparseF("fDeltaMassFullAnalysisOfflineBkg","hDeltaMassFullAnalysisOfflineBkg;inv mass (GeV/c);#Delta inv mass (GeV/c); p_{T}^{#Lambda_{c}} (GeV/c); #Delta p_{T} (GeV/c); daughter angle (2prongs) (rad);",5,binSparseDMassRot,edgeLowSparseDMassRot,edgeHighSparseDMassRot);

      fOutputPIDBachTR->Add(hDeltaMassFullAnalysisSgn);
      fOutputPIDBachTR->Add(hDeltaMassFullAnalysisOfflineSgn);
      fOutputPIDBachTR->Add(hDeltaMassFullAnalysisBkg);
      fOutputPIDBachTR->Add(hDeltaMassFullAnalysisOfflineBkg);

      }
    */

  }

  /*
    fOutputAll->Print();
    fOutputPIDBach->Print();
    if (fTrackRotation) fOutputPIDBachTR->Print();
  */
  return;
}

//---------------------------
void AliAnalysisTaskSELc2V0bachelor::CheckEventSelection(AliAODEvent *aodEvent) {
  //
  /// To fill control histograms
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
  ///
  /// To fill control histograms
  ///

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
  /// This is now implemented in AliAODRecoCascadeHF
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
  /// This is to check Lc dinasty
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
void AliAnalysisTaskSELc2V0bachelor::FillArmPodDistribution(AliAODRecoDecay *vZero,
							      TString histoTitle,
							      Bool_t isCandidateSelectedCuts,
							      Bool_t isBachelorID) {
  //
  /// This is to fill Armenteros Podolanski plots
  //

  Double_t alpha = vZero->Alpha();//AlphaV0();
  Double_t qT    = vZero->QtProng();//PtArmV0();

  ((TH2F*)(fOutputAll->FindObject(histoTitle+"0")))->Fill(alpha,qT);
  if (isCandidateSelectedCuts) {
    ((TH2F*)(fOutputAll->FindObject(histoTitle)))->Fill(alpha,qT);
    if (isBachelorID)
      ((TH2F*)(fOutputPIDBach->FindObject(histoTitle)))->Fill(alpha,qT);
  }

}

//-------------------------------------------------------------------------------
void AliAnalysisTaskSELc2V0bachelor::CheckCandidatesAtDifferentLevels(AliAODRecoCascadeHF *part, AliRDHFCutsLctoV0* cutsAnal) {
  //
  /// This is to check candidates at different levels
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
	if ( ( (aaa&AliRDHFCutsLctoV0::kLcToLpi)==AliRDHFCutsLctoV0::kLcToLpi && bachelor->Charge()<0)  ||
	     ( (aaa&AliRDHFCutsLctoV0::kLcToLBarpi)==AliRDHFCutsLctoV0::kLcToLBarpi && bachelor->Charge()>0) )
	  ((TH1F*)(fOutput->FindObject("hSwitchOnCandidates1")))->Fill( -aaa );
	else
	  ((TH1F*)(fOutput->FindObject("hSwitchOnCandidates1")))->Fill( aaa );
      }

      cutsAnal->SetUsePID(kFALSE);
      aaa = cutsAnal->IsSelected(part,AliRDHFCuts::kCandidate);
      if ((aaa&AliRDHFCutsLctoV0::kLcToK0Spr)==AliRDHFCutsLctoV0::kLcToK0Spr) {
	if ( ( (aaa&AliRDHFCutsLctoV0::kLcToLpi)==AliRDHFCutsLctoV0::kLcToLpi && bachelor->Charge()<0) ||
	     ( (aaa&AliRDHFCutsLctoV0::kLcToLBarpi)==AliRDHFCutsLctoV0::kLcToLBarpi && bachelor->Charge()>0) )
	  ((TH1F*)(fOutput->FindObject("hSwitchOnCandidates2")))->Fill( -aaa );
	else
	  ((TH1F*)(fOutput->FindObject("hSwitchOnCandidates2")))->Fill( aaa );
      }
      cutsAnal->SetUsePID(areCutsUsingPID);

      aaa = cutsAnal->IsSelected(part,AliRDHFCuts::kPID);
      if ((aaa&AliRDHFCutsLctoV0::kLcToK0Spr)==AliRDHFCutsLctoV0::kLcToK0Spr) {
	if ( ( (aaa&AliRDHFCutsLctoV0::kLcToLpi)==AliRDHFCutsLctoV0::kLcToLpi && bachelor->Charge()<0) ||
	     ( (aaa&AliRDHFCutsLctoV0::kLcToLBarpi)==AliRDHFCutsLctoV0::kLcToLBarpi && bachelor->Charge()>0) )
	  ((TH1F*)(fOutput->FindObject("hSwitchOnCandidates3")))->Fill( -aaa );
	else
	  ((TH1F*)(fOutput->FindObject("hSwitchOnCandidates3")))->Fill( aaa );
      }

      aaa = cutsAnal->IsSelected(part,AliRDHFCuts::kAll);
      if ((aaa&AliRDHFCutsLctoV0::kLcToK0Spr)==AliRDHFCutsLctoV0::kLcToK0Spr) {
	if ( ( (aaa&AliRDHFCutsLctoV0::kLcToLpi)==AliRDHFCutsLctoV0::kLcToLpi && bachelor->Charge()<0) ||
	     ( (aaa&AliRDHFCutsLctoV0::kLcToLBarpi)==AliRDHFCutsLctoV0::kLcToLBarpi && bachelor->Charge()>0) )
	  ((TH1F*)(fOutput->FindObject("hSwitchOnCandidates4")))->Fill( -aaa );
	else
	  ((TH1F*)(fOutput->FindObject("hSwitchOnCandidates4")))->Fill( aaa );
      }

    }
  }

  return;
}

//-------------------------------------------------------------------------------
void AliAnalysisTaskSELc2V0bachelor::FillTheTree(AliAODRecoCascadeHF *part, AliRDHFCutsLctoV0 *cutsAnal, TClonesArray *mcArray, Int_t isLc, Int_t checkLcOrigin) {
  //
  /// This is to fill tree
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
  Int_t isKstar12K0Spi=0, isKstar22K0Spi=0;
  Int_t mcLabel4 = -1;
  Int_t mcLabel5 = -1;
  Double_t ptCandByMC = 0.;//fmcPartCandidate->Pt();
  Double_t yCandByMC  = 0.;//fmcPartCandidate->Y();
  if (fUseMCInfo) {
    if (isLc) {
      Int_t pdgCand0 = 4122;
      Int_t pdgDgLctoV0bachelor0[2]={2212,310};
      Int_t pdgDgV0toDaughters0[2]={211,211};
      Int_t mcLabelLc2pK0S = part->MatchToMC(pdgCand0,pdgDgLctoV0bachelor0[1],pdgDgLctoV0bachelor0,pdgDgV0toDaughters0,mcArray,kTRUE);
      AliAODMCParticle *lambdaCpartMC = (AliAODMCParticle*)mcArray->At(mcLabelLc2pK0S);
      if (lambdaCpartMC) {
	ptCandByMC = lambdaCpartMC->Pt();
	yCandByMC  = lambdaCpartMC->Y();
      }
    }

    Int_t pdgCand = 4122;
    Int_t pdgDgLctoV0bachelor[2]={211,3122};
    Int_t pdgDgV0toDaughters[2]={2212,211};
    mcLabel = part->MatchToMC(pdgCand,pdgDgLctoV0bachelor[1],pdgDgLctoV0bachelor,pdgDgV0toDaughters,mcArray,kTRUE);
    if (mcLabel!=-1) {
      if (bachelor->Charge()<0) isLc2LBarpi=1;
      if (bachelor->Charge()>0) isLc2Lpi=1;
      AliAODMCParticle *lambdaCpartMC = (AliAODMCParticle*)mcArray->At(mcLabel);
      if (lambdaCpartMC) {
	ptCandByMC = lambdaCpartMC->Pt();
	yCandByMC  = lambdaCpartMC->Y();
      }
    }

    Int_t pdgCand2 = 411; // D+ -> pi+ K0S
    Int_t pdgCand3 = 431; // Ds+ -> K+ K0S
    Int_t pdgDgCand2[2]={211,310};
    Int_t pdgDgCand3[2]={321,310};
    pdgDgV0toDaughters[0]=211;
    pdgDgV0toDaughters[1]=211;
    mcLabel2 = part->MatchToMC(pdgCand2,pdgDgCand2[1],pdgDgCand2,pdgDgV0toDaughters,mcArray,kTRUE);
    mcLabel3 = part->MatchToMC(pdgCand3,pdgDgCand3[1],pdgDgCand3,pdgDgV0toDaughters,mcArray,kTRUE);
    if (mcLabel2!=-1) {
      isDp2K0Spi=1;
      AliAODMCParticle *lambdaCpartMC = (AliAODMCParticle*)mcArray->At(mcLabel2);
      if (lambdaCpartMC) {
	ptCandByMC = lambdaCpartMC->Pt();
	yCandByMC  = lambdaCpartMC->Y();
      }
    }
    if (mcLabel3!=-1) {
      isDs2K0SK=1;
      AliAODMCParticle *lambdaCpartMC = (AliAODMCParticle*)mcArray->At(mcLabel3);
      if (lambdaCpartMC) {
	ptCandByMC = lambdaCpartMC->Pt();
	yCandByMC  = lambdaCpartMC->Y();
      }
    }

    Int_t pdgCand4 = 313; // K*(892)+ -> pi+ K0S
    Int_t pdgCand5 = 325; // K*(1430)+ -> pi+ K0S
    Int_t pdgDgCand4[2]={211,310};
    Int_t pdgDgCand5[2]={211,310};
    pdgDgV0toDaughters[0]=211;
    pdgDgV0toDaughters[1]=211;
    mcLabel4 = part->MatchToMC(pdgCand4,pdgDgCand4[1],pdgDgCand4,pdgDgV0toDaughters,mcArray,kTRUE);
    mcLabel5 = part->MatchToMC(pdgCand5,pdgDgCand5[1],pdgDgCand5,pdgDgV0toDaughters,mcArray,kTRUE);
    if (mcLabel4!=-1) {
      isKstar12K0Spi=1;
      AliAODMCParticle *lambdaCpartMC = (AliAODMCParticle*)mcArray->At(mcLabel4);
      if (lambdaCpartMC) {
	ptCandByMC = lambdaCpartMC->Pt();
	yCandByMC  = lambdaCpartMC->Y();
      }
    }
    if (mcLabel5!=-1) {
      isKstar22K0Spi=1;
      AliAODMCParticle *lambdaCpartMC = (AliAODMCParticle*)mcArray->At(mcLabel5);
      if (lambdaCpartMC) {
	ptCandByMC = lambdaCpartMC->Pt();
	yCandByMC  = lambdaCpartMC->Y();
      }
    }
  }

  Int_t isLcByMC = isLc+isLc2LBarpi*2+isLc2Lpi*4+isDp2K0Spi*8+isDs2K0SK*16+isKstar12K0Spi*32+isKstar22K0Spi*64;

  Bool_t isMCparticleInFiducialAcceptance = kTRUE;
  if (isLc || isLc2LBarpi || isLc2Lpi || isDp2K0Spi || isDs2K0SK || isKstar12K0Spi || isKstar22K0Spi) {
    isMCparticleInFiducialAcceptance = cutsAnal->IsInFiducialAcceptance(ptCandByMC,yCandByMC);
  }

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
  areV0daughtersSelected += (v0pos->HasPointOnITSLayer(0))*16;
  areV0daughtersSelected += (v0pos->HasPointOnITSLayer(1))*32;
  areV0daughtersSelected += (v0pos->HasPointOnITSLayer(2))*64;
  areV0daughtersSelected += (v0pos->HasPointOnITSLayer(3))*128;
  areV0daughtersSelected += (v0pos->HasPointOnITSLayer(4))*256;
  areV0daughtersSelected += (v0pos->HasPointOnITSLayer(5))*512;

  areV0daughtersSelected += (v0neg->TestFilterMask(BIT(4)))*1024 + (!(v0neg->TestFilterMask(BIT(4))))*2048;
  areV0daughtersSelected += (v0neg->GetLabel()<0)*4096 + (v0neg->GetLabel()>=0)*8192;
  areV0daughtersSelected += (v0neg->HasPointOnITSLayer(0))*16384;
  areV0daughtersSelected += (v0neg->HasPointOnITSLayer(1))*32768;
  areV0daughtersSelected += (v0neg->HasPointOnITSLayer(2))*65536;
  areV0daughtersSelected += (v0neg->HasPointOnITSLayer(3))*131072;
  areV0daughtersSelected += (v0neg->HasPointOnITSLayer(4))*262144;
  areV0daughtersSelected += (v0neg->HasPointOnITSLayer(5))*524288;

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
  flagToCheckCandidate+=2*((TMath::Abs(invmassLambdaBar-mLPDG)<=0.050) && (bachelor->Charge()<0));
  flagToCheckCandidate+=4*((TMath::Abs(invmassLambda-mLPDG)<=0.050) && (bachelor->Charge()>0));
  flagToCheckCandidate+=8*((TMath::Abs(invmassLambdaBar-mLPDG)<=0.050) && (bachelor->Charge()>0));
  flagToCheckCandidate+=16*((TMath::Abs(invmassLambda-mLPDG)<=0.050) && (bachelor->Charge()<0));

  fCandidateVariables[ 0] = fUseMCInfo+isLcByMC; // 0: real data; 1: bkg; 2: Lc->K0S+p; 3: Lc->LambdaBar+pbar; 5: Lc->Lambda+p; 9: D+->K0S+pi; 17: Ds+->K0S+K; 33: K*+->K0S+pi; 65: K*+->K0S+K
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

  Double_t dcaForLc=0.;
  if (fAdditionalChecks) {
    Double_t xVtxLc=0, yVtxLc=0, zVtxLc=0;
    Double_t xLcMC=0,yLcMC=0,zLcMC=0;
    Double_t pxVtxBachelor=0, pyVtxBachelor=0, pzVtxBachelor=0;
    dcaForLc = PropagateToDCA(v0part,bachelor,fBzkG, xVtxLc, yVtxLc, zVtxLc, pxVtxBachelor, pyVtxBachelor, pzVtxBachelor);
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
    } else if (isKstar12K0Spi) {
      AliAODMCParticle *partLc = dynamic_cast<AliAODMCParticle*>(mcArray->At(mcLabel4));
      AliAODMCParticle *partLcDaug0 = dynamic_cast<AliAODMCParticle*>(mcArray->At(partLc->GetDaughter(0)));
      xLcMC=partLcDaug0->Xv(), yLcMC=partLcDaug0->Yv(), zLcMC=partLcDaug0->Zv();
    } else if (isKstar22K0Spi) {
      AliAODMCParticle *partLc = dynamic_cast<AliAODMCParticle*>(mcArray->At(mcLabel5));
      AliAODMCParticle *partLcDaug0 = dynamic_cast<AliAODMCParticle*>(mcArray->At(partLc->GetDaughter(0)));
      xLcMC=partLcDaug0->Xv(), yLcMC=partLcDaug0->Yv(), zLcMC=partLcDaug0->Zv();
    }

    fCandidateVariables[75]=part->GetSecVtxX();
    fCandidateVariables[76]=part->GetSecVtxY();
    fCandidateVariables[77]=part->GetSecVtxZ();
    fCandidateVariables[78]=xVtxLc;
    fCandidateVariables[79]=yVtxLc;
    fCandidateVariables[80]=zVtxLc;
    fCandidateVariables[81]=xLcMC;
    fCandidateVariables[82]=yLcMC;
    fCandidateVariables[83]=zLcMC;
    fCandidateVariables[84]=bachelor->Px();
    fCandidateVariables[85]=bachelor->Py();
    fCandidateVariables[86]=pxVtxBachelor;
    fCandidateVariables[87]=pyVtxBachelor;
    fCandidateVariables[88]=pzVtxBachelor;
    fCandidateVariables[89]=v0part->Px();
    fCandidateVariables[90]=v0part->Py();
    fCandidateVariables[91]=v0part->Pz();
    fCandidateVariables[92]=fVtx1->GetX();
    fCandidateVariables[93]=fVtx1->GetY();
    fCandidateVariables[94]=fVtx1->GetZ();
  }

  fCandidateVariables[55]=dcaForLc;

  fCandidateVariables[56]=part->CosThetaStar(0,4122,2212,310);
  fCandidateVariables[57]=part->CosThetaStar(1,4122,2212,310);
  fCandidateVariables[58]=v0part->Eta();
  fCandidateVariables[59]=v0part->Y(310);
  fCandidateVariables[60]=bachelor->Charge();
  fCandidateVariables[61]=isMCparticleInFiducialAcceptance;

  fCandidateVariables[62] = part->InvMass2Prongs(0,1,211,310); // Kstar( 892)+ -> pi+K0S
  fCandidateVariables[63] = part->InvMass2Prongs(0,1,321,310); // Kstar(1430)+ -> pi+K0S

  fCandidateVariables[64]=-1;
  fCandidateVariables[65]=-1;
  fCandidateVariables[66]=-1;
  fCandidateVariables[67]=-1;
  fCandidateVariables[68]=-1;
  if (fUseMCInfo) {
    if (bachelor->GetLabel()!=-1) {
      AliAODMCParticle *partBachelor = dynamic_cast<AliAODMCParticle*>(mcArray->At(TMath::Abs(bachelor->GetLabel())));
      if(partBachelor) fCandidateVariables[64]=partBachelor->GetPdgCode();
    }
    if (bachelor->GetLabel()!=-1 &&
	v0pos->GetLabel()!=-1 &&
	v0neg->GetLabel()!=-1) {
      const Int_t ndg=3;
      Int_t dgLabels[ndg]={TMath::Abs(bachelor->GetLabel()),
			   TMath::Abs(v0pos->GetLabel()),
			   TMath::Abs(v0neg->GetLabel())};
      Int_t ndgCk=0;
      Int_t *pdgDg=0;
      Int_t absLabelMother=-1;
      Int_t nDauCand=-1;
      fCandidateVariables[65]=SearchForCommonMother(mcArray,
						    dgLabels,ndg,ndgCk,pdgDg,absLabelMother,nDauCand);
    }
    if (v0pos->GetLabel()!=-1) {
      AliAODMCParticle *part1 = dynamic_cast<AliAODMCParticle*>(mcArray->At(TMath::Abs(v0pos->GetLabel())));
      if(part1) fCandidateVariables[66]=part1->GetPdgCode();
    }
    if (v0neg->GetLabel()!=-1) {
      AliAODMCParticle *part2 = dynamic_cast<AliAODMCParticle*>(mcArray->At(TMath::Abs(v0neg->GetLabel())));
      if(part2) fCandidateVariables[67]=part2->GetPdgCode();
    }
    if (v0pos->GetLabel()!=-1 &&
	v0neg->GetLabel()!=-1) {
      const Int_t ndg=2;
      Int_t dgLabels[ndg]={TMath::Abs(v0pos->GetLabel()),
			   TMath::Abs(v0neg->GetLabel())};
      Int_t ndgCk=0;
      Int_t *pdgDg=0;
      Int_t absLabelMother=-1;
      Int_t nDauCand=-1;
      fCandidateVariables[68]=SearchForCommonMother(mcArray,
						    dgLabels,ndg,ndgCk,pdgDg,absLabelMother,nDauCand);
    }
  }

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler *inputHandler=(AliInputEventHandler*)mgr->GetInputEventHandler();
  AliPIDResponse *pidResponse=inputHandler->GetPIDResponse();
  fCandidateVariables[69]=pidResponse->GetTOFResponse().GetStartTimeMask(bachelor->P());

  AliPIDCombined *objectPIDCombined=new AliPIDCombined;
  objectPIDCombined->SetDefaultTPCPriors();
  objectPIDCombined->SetDetectorMask(AliPIDResponse::kDetTPC+AliPIDResponse::kDetTOF);

  Double_t probTPCTOF[AliPID::kSPECIES]={-1.};
  UInt_t detUsed = objectPIDCombined->ComputeProbabilities(bachelor, pidResponse, probTPCTOF);

  Double_t probProton = -1.;
  //  Double_t probPion = -1.;
  //  Double_t probKaon = -1.;
  if (detUsed == (UInt_t)objectPIDCombined->GetDetectorMask() ) {
    AliDebug(2, Form("We have found the detector mask for TOF + TPC: probProton will be set to %f", probTPCTOF[AliPID::kProton]));
    probProton = probTPCTOF[AliPID::kProton];
    // probPion = probTPCTOF[AliPID::kPion];
    // probKaon = probTPCTOF[AliPID::kKaon];
  }
  else { // if you don't have both TOF and TPC, try only TPC
    objectPIDCombined->SetDetectorMask(AliPIDResponse::kDetTPC);
    AliDebug(2, "We did not find the detector mask for TOF + TPC, let's see only TPC");
    detUsed = objectPIDCombined->ComputeProbabilities(bachelor, pidResponse, probTPCTOF);
    AliDebug(2,Form(" detUsed (TPC case) = %d", detUsed));
    if (detUsed == (UInt_t)objectPIDCombined->GetDetectorMask()) {
      probProton = probTPCTOF[AliPID::kProton];
      // probPion = probTPCTOF[AliPID::kPion];
      // probKaon = probTPCTOF[AliPID::kKaon];
      AliDebug(2, Form("TPC only worked: probProton will be set to %f", probTPCTOF[AliPID::kProton]));
    }
    else {
      AliDebug(2, "Only TPC did not work...");
    }
    // resetting mask to ask for both TPC+TOF
    objectPIDCombined->SetDetectorMask(AliPIDResponse::kDetTPC+AliPIDResponse::kDetTOF);
  }
  AliDebug(2, Form("probProton = %f", probProton));

  // now we get the TPC and TOF single PID probabilities (only for Proton, or the tree will explode :) )
  Double_t probProtonTPC = -1.;
  Double_t probProtonTOF = -1.;
  Double_t pidTPC[AliPID::kSPECIES]={-1.};
  Double_t pidTOF[AliPID::kSPECIES]={-1.};
  Int_t respTPC = pidResponse->ComputePIDProbability(AliPIDResponse::kDetTPC, bachelor, AliPID::kSPECIES, pidTPC);
  Int_t respTOF = pidResponse->ComputePIDProbability(AliPIDResponse::kDetTOF, bachelor, AliPID::kSPECIES, pidTOF);
  if (respTPC == AliPIDResponse::kDetPidOk) probProtonTPC = pidTPC[AliPID::kProton];
  if (respTOF == AliPIDResponse::kDetPidOk) probProtonTOF = pidTOF[AliPID::kProton];

  fCandidateVariables[70]=probProton;
  fCandidateVariables[71]=probProtonTPC;
  fCandidateVariables[72]=probProtonTOF;

  fCandidateVariables[73]=checkLcOrigin;

  fCandidateVariables[74]=v0part->QtProng();

  delete objectPIDCombined;

  fVariablesTree->Fill();

  return;
}

//-------------------------------------------------------------------------------
void AliAnalysisTaskSELc2V0bachelor::DefineTreeVariables() {
  //
  /// This is to define tree variables
  //

  const char* nameoutput = GetOutputSlot(4)->GetContainer()->GetName();
  fVariablesTree = new TTree(nameoutput,"Candidates variables tree");
  Int_t nVar = 75;
  if (fAdditionalChecks) nVar = 95;
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

  fCandidateVariableNames[56]="cosThetaStarBachelor";
  fCandidateVariableNames[57]="cosThetaStarV0";
  fCandidateVariableNames[58]="etaV0";
  fCandidateVariableNames[59]="yV0";
  fCandidateVariableNames[60]="bachelorCharge";
  fCandidateVariableNames[61]="isMCparticleInFiducialAcceptance";

  fCandidateVariableNames[62]="massKstar12K0Spi"; // Kstar( 892)+ -> pi+ K0S
  fCandidateVariableNames[63]="massKstar22K0Spi"; // Kstar(1430)+ -> pi+ K0S
  fCandidateVariableNames[64]="pdgBachelor"; // pdg MC bachelor
  fCandidateVariableNames[65]="pdgCandidate"; // pdg MC candidate recovered via new method
  fCandidateVariableNames[66]="pdgV0pos"; // pdg MC V0 positive
  fCandidateVariableNames[67]="pdgV0neg"; // pdg MC V0 negative
  fCandidateVariableNames[68]="pdgV0Candidate"; // pdg MC V0candidate recovered via new method
  fCandidateVariableNames[69]="startTimeMask"; // start time mask

  fCandidateVariableNames[70]="combinedProtonProb";
  fCandidateVariableNames[71]="TPCProtonProb";
  fCandidateVariableNames[72]="TOFProtonProb";
  fCandidateVariableNames[73]="checkLcOrigin";

  fCandidateVariableNames[74]="qtProng0V0";

  if (fAdditionalChecks) {
    fCandidateVariableNames[75]="xVtxLcBad";
    fCandidateVariableNames[76]="yVtxLcBad";
    fCandidateVariableNames[77]="zVtxLcBad";
    fCandidateVariableNames[78]="xVtxLcGood";
    fCandidateVariableNames[79]="yVtxLcGood";
    fCandidateVariableNames[80]="zVtxLcGood";
    fCandidateVariableNames[81]="xVtxLcMC";
    fCandidateVariableNames[82]="yVtxLcMC";
    fCandidateVariableNames[83]="zVtxLcMC";
    fCandidateVariableNames[84]="pxVtxBachelorBad";
    fCandidateVariableNames[85]="pyVtxBachelorBad";
    fCandidateVariableNames[86]="pxVtxBachelorGood";
    fCandidateVariableNames[87]="pyVtxBachelorGood";
    fCandidateVariableNames[88]="pzVtxBachelorGood";
    fCandidateVariableNames[89]="pxVtxV0";
    fCandidateVariableNames[90]="pyVtxV0";
    fCandidateVariableNames[91]="pzVtxV0";
    fCandidateVariableNames[92]="xPvtx";
    fCandidateVariableNames[93]="yPvtx";
    fCandidateVariableNames[94]="zPvtx";
  }

  for (Int_t ivar=0; ivar<nVar; ivar++) {
    fVariablesTree->Branch(fCandidateVariableNames[ivar].Data(),&fCandidateVariables[ivar],Form("%s/f",fCandidateVariableNames[ivar].Data()));
  }

  return;
}

//__________________________________________________________________________
void  AliAnalysisTaskSELc2V0bachelor::DefineGeneralHistograms() {
  //
  /// This is to define general histograms
  //

  fCEvents = new TH1F("fCEvents","conter",20,0,20);
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
  fCEvents->GetXaxis()->SetBinLabel(18,Form("zVtxMC<=%2.0fcm",fAnalCuts->GetMaxVtxZ()));
  fCEvents->GetXaxis()->SetBinLabel(19,"Re-Fill Fail");
  fCEvents->GetXaxis()->SetBinLabel(20,"AOD Mismatch");
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
void  AliAnalysisTaskSELc2V0bachelor::FillAnalysisHistograms(AliAODRecoCascadeHF *part, AliRDHFCutsLctoV0 *cutsAnal, TString appendthis) {
  //
  // This is to fill analysis histograms
  //

  TString fillthis="";

  Bool_t isBachelorID = (((cutsAnal->IsSelected(part,AliRDHFCuts::kPID))&(AliRDHFCutsLctoV0::kLcToK0Spr))==(AliRDHFCutsLctoV0::kLcToK0Spr)); // ID x bachelor

  Bool_t areCutsUsingPID = cutsAnal->GetIsUsePID();
  cutsAnal->SetUsePID(kFALSE);

  Double_t invmassLc = part->InvMassLctoK0sP();
  Double_t lambdacpt = part->Pt();

  AliAODTrack *bachelor = (AliAODTrack*)part->GetBachelor();
  Double_t momBach = bachelor->P();
  Double_t ptBach = bachelor->Pt();

  AliAODv0 *v0part = (AliAODv0*)part->Getv0();
  Double_t momK0S = v0part->P();
  Double_t ptK0S = v0part->Pt();
  //Double_t dcaV0ptp = v0part->GetDCA();
  Double_t invmassK0S = v0part->MassK0Short();

  AliAODTrack *v0pos = (AliAODTrack*)part->Getv0PositiveTrack();
  Double_t ptV0pos = v0pos->Pt();
  AliAODTrack *v0neg = (AliAODTrack*)part->Getv0NegativeTrack();
  Double_t ptV0neg = v0neg->Pt();

  if (!appendthis.Contains("SgnC") && !appendthis.Contains("SgnB") && !appendthis.Contains("SgnNoQ")) {
    fillthis="histpK0Svsp"+appendthis;
    //cout << fillthis << endl;
    if ( ((cutsAnal->IsSelected(part,AliRDHFCuts::kCandidate))&(AliRDHFCutsLctoV0::kLcToK0Spr)) == (AliRDHFCutsLctoV0::kLcToK0Spr) ) {
      ((TH2F*)(fOutputAll->FindObject(fillthis)))->Fill(momBach,momK0S);
      if (isBachelorID)  ((TH2F*)(fOutputPIDBach->FindObject(fillthis)))->Fill(momBach,momK0S);
    }
  }

  fillthis="histLcMassByK0S"+appendthis;
  //cout << fillthis << endl;
  if ( ((cutsAnal->IsSelected(part,AliRDHFCuts::kCandidate))&(AliRDHFCutsLctoV0::kLcToK0Spr)) == (AliRDHFCutsLctoV0::kLcToK0Spr) ) {
    ((TH2F*)(fOutputAll->FindObject(fillthis)))->Fill(invmassLc,lambdacpt);
    if (isBachelorID)((TH2F*)(fOutputPIDBach->FindObject(fillthis)))->Fill(invmassLc,lambdacpt);
  }

  if (!appendthis.Contains("SgnC") && !appendthis.Contains("SgnB") && !appendthis.Contains("SgnNoQ")) {
    fillthis="histK0SMass"+appendthis;
    //    cout << fillthis << endl;
    cutsAnal->SetExcludedCut(2);
    if ( ((cutsAnal->IsSelected(part,AliRDHFCuts::kCandidate))&(AliRDHFCutsLctoV0::kLcToK0Spr)) == (AliRDHFCutsLctoV0::kLcToK0Spr) ) {
      ((TH2F*)(fOutputAll->FindObject(fillthis)))->Fill(lambdacpt,invmassK0S);
      if (isBachelorID)  ((TH2F*)(fOutputPIDBach->FindObject(fillthis)))->Fill(lambdacpt,invmassK0S);
    }
    cutsAnal->SetExcludedCut(-1);
  }

  fillthis="histptK0S"+appendthis;
  //cout << fillthis << endl;
  cutsAnal->SetExcludedCut(15);
  if ( ((cutsAnal->IsSelected(part,AliRDHFCuts::kCandidate))&(AliRDHFCutsLctoV0::kLcToK0Spr)) == (AliRDHFCutsLctoV0::kLcToK0Spr) ) {
    ((TH2F*)(fOutputAll->FindObject(fillthis)))->Fill(lambdacpt,ptK0S);
    if (isBachelorID)((TH2F*)(fOutputPIDBach->FindObject(fillthis)))->Fill(lambdacpt,ptK0S);
  }

  fillthis="histptP"+appendthis;
  //cout << fillthis << endl;
  cutsAnal->SetExcludedCut(4);
  if ( ((cutsAnal->IsSelected(part,AliRDHFCuts::kCandidate))&(AliRDHFCutsLctoV0::kLcToK0Spr)) == (AliRDHFCutsLctoV0::kLcToK0Spr) ) {
    ((TH2F*)(fOutputAll->FindObject(fillthis)))->Fill(lambdacpt,ptBach);
    if (isBachelorID)((TH2F*)(fOutputPIDBach->FindObject(fillthis)))->Fill(lambdacpt,ptBach);
  }

  fillthis="histptPip"+appendthis;
  //cout << fillthis << endl;
  cutsAnal->SetExcludedCut(5);
  if ( ((cutsAnal->IsSelected(part,AliRDHFCuts::kCandidate))&(AliRDHFCutsLctoV0::kLcToK0Spr)) == (AliRDHFCutsLctoV0::kLcToK0Spr) ) {
    ((TH2F*)(fOutputAll->FindObject(fillthis)))->Fill(lambdacpt,ptV0pos);
    if (isBachelorID)((TH2F*)(fOutputPIDBach->FindObject(fillthis)))->Fill(lambdacpt,ptV0pos);
  }

  fillthis="histptPim"+appendthis;
  //cout << fillthis << endl;
  cutsAnal->SetExcludedCut(6);
  if ( ((cutsAnal->IsSelected(part,AliRDHFCuts::kCandidate))&(AliRDHFCutsLctoV0::kLcToK0Spr)) == (AliRDHFCutsLctoV0::kLcToK0Spr) ) {
    ((TH2F*)(fOutputAll->FindObject(fillthis)))->Fill(lambdacpt,ptV0neg);
    if (isBachelorID)((TH2F*)(fOutputPIDBach->FindObject(fillthis)))->Fill(lambdacpt,ptV0neg);
  }

  if (!appendthis.Contains("SgnC") && !appendthis.Contains("SgnB") && !appendthis.Contains("SgnNoQ")) {
    fillthis="histLambdaMass"+appendthis;
    //    cout << fillthis << endl;
    cutsAnal->SetExcludedCut(13);
    if ( ((cutsAnal->IsSelected(part,AliRDHFCuts::kCandidate))&(AliRDHFCutsLctoV0::kLcToK0Spr)) == (AliRDHFCutsLctoV0::kLcToK0Spr) ) {
      ((TH2F*)(fOutputAll->FindObject(fillthis)))->Fill(lambdacpt,v0part->MassLambda());
      if (isBachelorID)((TH2F*)(fOutputPIDBach->FindObject(fillthis)))->Fill(lambdacpt,v0part->MassLambda());
    }

    fillthis="histLambdaBarMass"+appendthis;
    //    cout << fillthis << endl;
    cutsAnal->SetExcludedCut(13);
    if ( ((cutsAnal->IsSelected(part,AliRDHFCuts::kCandidate))&(AliRDHFCutsLctoV0::kLcToK0Spr)) == (AliRDHFCutsLctoV0::kLcToK0Spr) ) {
      ((TH2F*)(fOutputAll->FindObject(fillthis)))->Fill(lambdacpt,v0part->MassAntiLambda());
      if (isBachelorID)((TH2F*)(fOutputPIDBach->FindObject(fillthis)))->Fill(lambdacpt,v0part->MassAntiLambda());
    }

    fillthis="histGammaMass"+appendthis;
    //    cout << fillthis << endl;
    cutsAnal->SetExcludedCut(14);
    if ( ((cutsAnal->IsSelected(part,AliRDHFCuts::kCandidate))&(AliRDHFCutsLctoV0::kLcToK0Spr)) == (AliRDHFCutsLctoV0::kLcToK0Spr) ) {
      ((TH2F*)(fOutputAll->FindObject(fillthis)))->Fill(lambdacpt,v0part->InvMass2Prongs(0,1,11,11));
      if (isBachelorID)((TH2F*)(fOutputPIDBach->FindObject(fillthis)))->Fill(lambdacpt,v0part->InvMass2Prongs(0,1,11,11));
    }
  }

  fillthis="histD0K0S"+appendthis;
  //cout << fillthis << endl;
  cutsAnal->SetExcludedCut(11);
  if ( ((cutsAnal->IsSelected(part,AliRDHFCuts::kCandidate))&(AliRDHFCutsLctoV0::kLcToK0Spr)) == (AliRDHFCutsLctoV0::kLcToK0Spr) ) {
    ((TH2F*)(fOutputAll->FindObject(fillthis)))->Fill(lambdacpt,part->Getd0Prong(1));
    if (isBachelorID)((TH2F*)(fOutputPIDBach->FindObject(fillthis)))->Fill(lambdacpt,part->Getd0Prong(1));
  }

  fillthis="histD0P"+appendthis;
  //cout << fillthis << endl;
  cutsAnal->SetExcludedCut(10);
  if ( ((cutsAnal->IsSelected(part,AliRDHFCuts::kCandidate))&(AliRDHFCutsLctoV0::kLcToK0Spr)) == (AliRDHFCutsLctoV0::kLcToK0Spr) ) {
    ((TH2F*)(fOutputAll->FindObject(fillthis)))->Fill(lambdacpt,part->Getd0Prong(0));
    if (isBachelorID)((TH2F*)(fOutputPIDBach->FindObject(fillthis)))->Fill(lambdacpt,part->Getd0Prong(0));
  }

  fillthis="histCosPAK0S"+appendthis;
  //cout << fillthis << endl;
  cutsAnal->SetExcludedCut(9);
  if ( ((cutsAnal->IsSelected(part,AliRDHFCuts::kCandidate))&(AliRDHFCutsLctoV0::kLcToK0Spr)) == (AliRDHFCutsLctoV0::kLcToK0Spr) ) {
    ((TH2F*)(fOutputAll->FindObject(fillthis)))->Fill(lambdacpt,part->CosV0PointingAngle());
    if (isBachelorID)((TH2F*)(fOutputPIDBach->FindObject(fillthis)))->Fill(lambdacpt,part->CosV0PointingAngle());
  }

  fillthis="histCosThetaProtonCMS"+appendthis;
  //cout << fillthis << endl;
  cutsAnal->SetExcludedCut(16);
  if ( ((cutsAnal->IsSelected(part,AliRDHFCuts::kCandidate))&(AliRDHFCutsLctoV0::kLcToK0Spr)) == (AliRDHFCutsLctoV0::kLcToK0Spr) ) {
    ((TH2F*)(fOutputAll->FindObject(fillthis)))->Fill(lambdacpt,cutsAnal->GetProtonEmissionAngleCMS(part));
    if (isBachelorID)((TH2F*)(fOutputPIDBach->FindObject(fillthis)))->Fill(lambdacpt,cutsAnal->GetProtonEmissionAngleCMS(part));
  }

  fillthis="histResignedD0"+appendthis;
  //cout << fillthis << endl;
  cutsAnal->SetExcludedCut(18);
  if ( ((cutsAnal->IsSelected(part,AliRDHFCuts::kCandidate))&(AliRDHFCutsLctoV0::kLcToK0Spr)) == (AliRDHFCutsLctoV0::kLcToK0Spr) ) {
    ((TH2F*)(fOutputAll->FindObject(fillthis)))->Fill(lambdacpt,cutsAnal->GetReSignedd0(part));
    if (isBachelorID)((TH2F*)(fOutputPIDBach->FindObject(fillthis)))->Fill(lambdacpt,cutsAnal->GetReSignedd0(part));
  }

  cutsAnal->SetExcludedCut(-1);

  cutsAnal->SetUsePID(areCutsUsingPID);

  return;
}
//---------------------------
Double_t AliAnalysisTaskSELc2V0bachelor::PropagateToDCA(AliAODv0 *v, AliAODTrack *bachelor, Double_t b,
							  Double_t &xVtxLc, Double_t &yVtxLc, Double_t &zVtxLc,
							  Double_t &pxVtxBachelor, Double_t &pyVtxBachelor, Double_t &pzVtxBachelor) {
  //--------------------------------------------------------------------
  /// This function returns the DCA between the V0 and the track
  /// This is a copy of AliCascadeVertexer::PropagateToDCA(...) method
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
  /// To estimate alpha according to what done in the AliExternalTrackParam::Set(...) method
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
  /// This function calculates locally a 2x2 determinant.
  /// This is a copy of the AliCascadeVertexer::Det(...) method
  //--------------------------------------------------------------------
  return a00*a11 - a01*a10;
}

//---------------------------
Double_t AliAnalysisTaskSELc2V0bachelor::Det(Double_t a00,Double_t a01,Double_t a02,
					       Double_t a10,Double_t a11,Double_t a12,
					       Double_t a20,Double_t a21,Double_t a22) const {
  //--------------------------------------------------------------------
  /// This function calculates locally a 3x3 determinant
  /// This is a copy of the AliCascadeVertexer::Det(...) method
  //--------------------------------------------------------------------
  return  a00*Det(a11,a12,a21,a22)-a01*Det(a10,a12,a20,a22)+a02*Det(a10,a11,a20,a21);
}

//----------------------------------------------------------------------------
Int_t AliAnalysisTaskSELc2V0bachelor::MatchToMClabelC(AliAODRecoCascadeHF *candidate,
							TClonesArray *mcArray)
{
  //
  /// Check if this candidate is matched to a MC signal  Lc -> p K0S + X
  /// If no, return -1
  /// If yes, return label (>=0) of the AliAODMCParticle
  //

  AliAODv0 *theV0 = dynamic_cast<AliAODv0*>(candidate->Getv0()); // the V0
  AliVTrack *trk = dynamic_cast<AliVTrack*>(candidate->GetBachelor()); // the bachelor
  if (!trk || !theV0) return -1;

  if (trk->GetLabel()==-1) return -1;
  Int_t bachLabels = TMath::Abs(trk->GetLabel());
  AliAODMCParticle*bachelorMC = dynamic_cast<AliAODMCParticle*>(mcArray->At(bachLabels));
  if (!bachelorMC) return -1;
  if (TMath::Abs(bachelorMC->GetPdgCode())!=2212) return -1;
  Int_t indexMotherBach = bachelorMC->GetMother();
  if (indexMotherBach==-1) return -1;

  Int_t pdgDg2prong[2] = {211,211};
  Int_t lab2Prong = theV0->MatchToMC(310,mcArray,2,pdgDg2prong); // the V0
  if(lab2Prong<0) return -1;
  AliAODMCParticle*partK0S = dynamic_cast<AliAODMCParticle*>(mcArray->At(lab2Prong));
  if (!partK0S) return -1;
  Int_t indexMotherK0S = partK0S->GetMother();
  if (indexMotherK0S==-1) return -1;
  AliAODMCParticle*partK0 = dynamic_cast<AliAODMCParticle*>(mcArray->At(indexMotherK0S));
  if (!partK0) return -1;
  Int_t indexMotherK0 = partK0->GetMother();
  if (indexMotherK0==-1) return -1;

  if (indexMotherBach!=indexMotherK0) return -1; // p e K0S sono fratelli

  AliAODMCParticle*partLc = dynamic_cast<AliAODMCParticle*>(mcArray->At(indexMotherK0));
  if (!partLc) return -1;
  Int_t ndg2 = partLc->GetDaughter(1)-partLc->GetDaughter(0)+1;
  if (ndg2==2) return -1;

  TString stringaCheck = Form(">>>>>>>> %d -> ",partLc->GetPdgCode());
  for(Int_t ii=0; ii<ndg2; ii++) {
    AliAODMCParticle* partDau=(AliAODMCParticle*)(mcArray->At(partLc->GetDaughter(0)+ii));
    stringaCheck.Append(Form("  %d",partDau->GetPdgCode()));
  }
  //printf("%s \n",stringaCheck.Data());

  return indexMotherBach;

}
//--------------------------------------------------------------------------
Int_t AliAnalysisTaskSELc2V0bachelor::SearchForCommonMother(TClonesArray *mcArray,
							      Int_t dgLabels[10],Int_t ndg,
							      Int_t &ndgCk, Int_t *pdgDg, Int_t &absLabelMother, Int_t &nDauCand) const
{
  ///
  /// Check if this candidate is matched to a MC signal
  /// If no, return 0
  /// If yes, return pdgCode of particle
  ///

  Int_t lab=-1,labMother=-1,pdgMother=0;
  AliAODMCParticle *part=0;
  AliAODMCParticle *mother=0;

  // loop on daughter labels
  TArrayI **labelMother = new TArrayI*[ndg];
  for(Int_t i=0; i<ndg; i++) labelMother[i] = new TArrayI(0);
  for(Int_t i=0; i<ndg; i++) {
    lab = TMath::Abs(dgLabels[i]);
    if(lab<0) {
      AliDebug(2,Form("daughter with negative label %d",lab));
      delete [] labelMother;
      return 0;
    }
    part = (AliAODMCParticle*)mcArray->At(lab);
    if(!part) {
      AliDebug(2,"no MC particle");
      delete [] labelMother;
      return 0;
    }

    mother = part;
    while(mother->GetMother()>=0) {
      labMother=mother->GetMother();
      mother = (AliAODMCParticle*)mcArray->At(labMother);
      if(!mother) {
	AliDebug(2,"no MC mother particle");
	break;
      }
      pdgMother = TMath::Abs(mother->GetPdgCode());
      if (pdgMother<10 || (pdgMother>18 && pdgMother<111)) {
	break;
      }
      labelMother[i]->Set(labelMother[i]->GetSize()+1);
      labelMother[i]->AddAt(labMother,labelMother[i]->GetSize()-1);
    }

  } // end loop on daughters


  TString stringaCheck;
  for(Int_t i=0; i<ndg; i++) {
    AliAODMCParticle*part0 = (AliAODMCParticle*)mcArray->At(TMath::Abs(dgLabels[i]));
    stringaCheck.Append(Form("part[%d]->GetLabel()=%d(%d) | ",i,dgLabels[i],part0->GetPdgCode()));
    stringaCheck.Append(Form("labelMother[%d] = ",i));
    for (Int_t jj=0;jj<labelMother[i]->GetSize(); jj++)
      stringaCheck.Append(Form("%d, ",labelMother[i]->At(jj)));
  }
  AliDebug(2,Form("%s \n",stringaCheck.Data()));
  Int_t pdgToBeReturned=0;

  TString stringaCheck2;
  ndgCk=ndg;
  pdgDg = new Int_t[ndgCk];
  for (Int_t index=1; index<ndg; index++) {
    Bool_t found=kFALSE;
    for (Int_t jj=0;jj<labelMother[index]->GetSize(); jj++) {
      for (Int_t ii=0;ii<labelMother[0]->GetSize(); ii++) {
	if (labelMother[0]->At(ii)==labelMother[index]->At(jj) &&
	    labelMother[0]->At(ii)!=0 && labelMother[0]->At(ii)!=1 && !found) {
	  mother = (AliAODMCParticle*)mcArray->At(labelMother[0]->At(ii));
	  pdgToBeReturned=mother->GetPdgCode();
	  absLabelMother=labelMother[0]->At(ii);
	  AliDebug(2,Form("FOUND label for the mother of this candidate: %d (PDG=%d)\n",labelMother[0]->At(ii),pdgToBeReturned));
	  //mother->Print();
	  nDauCand=mother->GetNDaughters();
	  found = kTRUE;
	  AliAODMCParticle *partMC = (AliAODMCParticle*)mcArray->At(dgLabels[0]);
	  pdgDg[0]=partMC->GetPdgCode();
	  partMC = (AliAODMCParticle*)mcArray->At(dgLabels[index]);
	  pdgDg[index]=partMC->GetPdgCode();
	  if (index==1) stringaCheck2.Append(Form("found daughters -> %d(%d)",dgLabels[0],pdgDg[0]));
	  stringaCheck2.Append(Form(" %d(%d)",dgLabels[index],pdgDg[index]));
	  break;
	}
      }
      if (found) break;
    }
  }
  stringaCheck2.Prepend(Form("Ecco quanto trovato: %d(%d) with %d daughters; ",absLabelMother,pdgToBeReturned,nDauCand));
  AliDebug(2,Form("%s \n",stringaCheck2.Data()));

  delete [] labelMother;
  delete [] pdgDg;

  return pdgToBeReturned;

}

void AliAnalysisTaskSELc2V0bachelor::TrackRotation(AliRDHFCutsLctoV0 * cuts, AliAODRecoCascadeHF *part, TString appendthis)
{

  AliAODRecoCascadeHF *partCopy = new AliAODRecoCascadeHF(*part);

  Double_t px[2]={partCopy->PxProng(0),partCopy->PxProng(1)};
  Double_t py[2]={partCopy->PyProng(0),partCopy->PyProng(1)};
  Double_t pz[2]={partCopy->PzProng(0),partCopy->PzProng(1)};

  Double_t pt = partCopy->Pt();
  Int_t pdgD=4122;
  UInt_t pdgLc2pK0S[2]={2212,310};
  Double_t minv2 = partCopy->InvMass2(2,pdgLc2pK0S);
  Double_t mass=TMath::Sqrt(minv2);
  Double_t rapid = partCopy->Y(pdgD);

  TString fillthis;

  if ( ( ( (cuts->IsSelected(part,AliRDHFCuts::kAll))&(AliRDHFCutsLctoV0::kLcToK0Spr)) == (AliRDHFCutsLctoV0::kLcToK0Spr) ) ) {
    fillthis="hMassVsPtVsY"+appendthis;
    //cout << fillthis << endl;
    ((TH3F*)(fOutputPIDBachTR->FindObject(fillthis)))->Fill(mass,pt,rapid);

    fillthis="phiVSthetaVSpt"+appendthis;
    //cout << fillthis << endl;
    ((TH3F*)(fOutputPIDBachTR->FindObject(fillthis)))->Fill(pt,part->Phi(),part->Theta());
  }

  Int_t nRotated=0;
  Double_t massRot=0;// calculated later only if candidate is acceptable
  //  Double_t angleProngXY=TMath::ACos((px[0]*px[1]+py[0]*py[1])/TMath::Sqrt((px[0]*px[0]+py[0]*py[0])*(px[1]*px[1]+py[1]*py[1])));
  //  Double_t ptOrig=pt;
  Double_t rotStep=(fMaxAngleForRot-fMinAngleForRot)/(fNRotations-1); // -1 is to ensure that the last rotation is done with angle=fMaxAngleForRot

  for(Int_t irot=0; irot<fNRotations; irot++){
    Double_t phirot=fMinAngleForRot+rotStep*irot;
    Double_t tmpx=px[0];
    Double_t tmpy=py[0];
    px[0]=tmpx*TMath::Cos(phirot)-tmpy*TMath::Sin(phirot);
    py[0]=tmpx*TMath::Sin(phirot)+tmpy*TMath::Cos(phirot);
    partCopy->SetPxPyPzProngs(2,px,py,pz);
    pt = partCopy->Pt();
    minv2 = partCopy->InvMass2(2,pdgLc2pK0S);
    massRot=TMath::Sqrt(minv2);
    rapid = partCopy->Y(pdgD);
    //if(minv2>fMinMass*fMinMass && minv2<fMaxMass*fMaxMass){
    if ( cuts->IsInFiducialAcceptance(pt,partCopy->Y(4122)) ) {
      if ( ((cuts->IsSelected(partCopy,AliRDHFCuts::kAll))&(AliRDHFCutsLctoV0::kLcToK0Spr)) == (AliRDHFCutsLctoV0::kLcToK0Spr) ) {

	fillthis="histLcMassByK0S"+appendthis;
	//cout << fillthis << endl;
	((TH2F*)(fOutputPIDBachTR->FindObject(fillthis)))->Fill(massRot,pt);

	fillthis="hMassVsPtVsYRot"+appendthis;
	//cout << fillthis << endl;
	((TH3F*)(fOutputPIDBachTR->FindObject(fillthis)))->Fill(mass,pt,rapid);

	fillthis="phiVSthetaVSptRot"+appendthis;
	//cout << fillthis << endl;
	((TH3F*)(fOutputPIDBachTR->FindObject(fillthis)))->Fill(pt,partCopy->Phi(),partCopy->Theta());

	fillthis="hDeltaMass"+appendthis;
	//cout << fillthis << endl;
	((TH1F*)(fOutputPIDBachTR->FindObject(fillthis)))->Fill(massRot-mass);
	//if(fFullAnalysis){
	//Double_t pointRot[5]={mass,massRot-mass,ptOrig,pt-ptOrig,angleProngXY};
	//fillthis="hDeltaMassFullAnalysis"+appendthis;
	////cout << fillthis << endl;
	//((THnSparse*)(fOutputPIDBachTR->FindObject(fillthis)))->Fill(pointRot);
	//}
	nRotated++;
	//}
      }

      // fill additional histos for track-rotated candidates
      fillthis="histptK0S"+appendthis;
      //cout << fillthis << endl;
      cuts->SetExcludedCut(15);
      if ( ((cuts->IsSelected(partCopy,AliRDHFCuts::kAll))&(AliRDHFCutsLctoV0::kLcToK0Spr)) == (AliRDHFCutsLctoV0::kLcToK0Spr) ) {
	((TH2F*)(fOutputPIDBachTR->FindObject(fillthis)))->Fill(pt,TMath::Sqrt(px[1]*px[1]+py[1]*py[1]));
      }

      fillthis="histptP"+appendthis;
      //cout << fillthis << endl;
      cuts->SetExcludedCut(4);
      if ( ((cuts->IsSelected(partCopy,AliRDHFCuts::kAll))&(AliRDHFCutsLctoV0::kLcToK0Spr)) == (AliRDHFCutsLctoV0::kLcToK0Spr) ) {
	((TH2F*)(fOutputPIDBachTR->FindObject(fillthis)))->Fill(pt,TMath::Sqrt(px[0]*px[0]+py[0]*py[0]));
      }

      fillthis="histptPip"+appendthis;
      //cout << fillthis << endl;
      cuts->SetExcludedCut(5);
      if ( ((cuts->IsSelected(partCopy,AliRDHFCuts::kAll))&(AliRDHFCutsLctoV0::kLcToK0Spr)) == (AliRDHFCutsLctoV0::kLcToK0Spr) ) {
	((TH2F*)(fOutputPIDBachTR->FindObject(fillthis)))->Fill(pt,(partCopy->Getv0PositiveTrack())->Pt());
      }

      fillthis="histptPim"+appendthis;
      //cout << fillthis << endl;
      cuts->SetExcludedCut(6);
      if ( ((cuts->IsSelected(partCopy,AliRDHFCuts::kAll))&(AliRDHFCutsLctoV0::kLcToK0Spr)) == (AliRDHFCutsLctoV0::kLcToK0Spr) ) {
	((TH2F*)(fOutputPIDBachTR->FindObject(fillthis)))->Fill(pt,(partCopy->Getv0NegativeTrack())->Pt());
      }

      fillthis="histLambdaMass"+appendthis;
      //cout << fillthis << endl;
      cuts->SetExcludedCut(13);
      if ( ((cuts->IsSelected(partCopy,AliRDHFCuts::kAll))&(AliRDHFCutsLctoV0::kLcToK0Spr)) == (AliRDHFCutsLctoV0::kLcToK0Spr) ) {
	((TH2F*)(fOutputPIDBachTR->FindObject(fillthis)))->Fill(pt,(partCopy->Getv0())->MassLambda());
      }

      fillthis="histLambdaBarMass"+appendthis;
      //cout << fillthis << endl;
      cuts->SetExcludedCut(13);
      if ( ((cuts->IsSelected(partCopy,AliRDHFCuts::kAll))&(AliRDHFCutsLctoV0::kLcToK0Spr)) == (AliRDHFCutsLctoV0::kLcToK0Spr) ) {
	((TH2F*)(fOutputPIDBachTR->FindObject(fillthis)))->Fill(pt,(partCopy->Getv0())->MassAntiLambda());
      }

      fillthis="histGammaMass"+appendthis;
      //cout << fillthis << endl;
      cuts->SetExcludedCut(14);
      if ( ((cuts->IsSelected(partCopy,AliRDHFCuts::kAll))&(AliRDHFCutsLctoV0::kLcToK0Spr)) == (AliRDHFCutsLctoV0::kLcToK0Spr) ) {
	((TH2F*)(fOutputPIDBachTR->FindObject(fillthis)))->Fill(pt,(partCopy->Getv0())->InvMass2Prongs(0,1,11,11));
      }

      fillthis="histCosPAK0S"+appendthis;
      //cout << fillthis << endl;
      cuts->SetExcludedCut(9);
      if ( ((cuts->IsSelected(partCopy,AliRDHFCuts::kAll))&(AliRDHFCutsLctoV0::kLcToK0Spr)) == (AliRDHFCutsLctoV0::kLcToK0Spr) ) {
	((TH2F*)(fOutputPIDBachTR->FindObject(fillthis)))->Fill(pt,partCopy->CosV0PointingAngle());
      }

      fillthis="histCosThetaProtonCMS"+appendthis;
      //cout << fillthis << endl;
      cuts->SetExcludedCut(16);
      if ( ((cuts->IsSelected(partCopy,AliRDHFCuts::kAll))&(AliRDHFCutsLctoV0::kLcToK0Spr)) == (AliRDHFCutsLctoV0::kLcToK0Spr) ) {
        ((TH2F*)(fOutputPIDBachTR->FindObject(fillthis)))->Fill(pt,cuts->GetProtonEmissionAngleCMS(partCopy));
      }

      fillthis="histResignedD0"+appendthis;
      //cout << fillthis << endl;
      cuts->SetExcludedCut(18);
      if ( ((cuts->IsSelected(partCopy,AliRDHFCuts::kAll))&(AliRDHFCutsLctoV0::kLcToK0Spr)) == (AliRDHFCutsLctoV0::kLcToK0Spr) ) {
        ((TH2F*)(fOutputPIDBach->FindObject(fillthis)))->Fill(pt,cuts->GetReSignedd0(partCopy));
      }
      cuts->SetExcludedCut(-1);

    } // isInFiducialAcceptance

    px[0]=tmpx;
    py[0]=tmpy;
  }
  fillthis="hNormRotated"+appendthis;
  //cout << fillthis << endl;
  ((TH1F*)(fOutputPIDBachTR->FindObject(fillthis)))->Fill(nRotated);

  delete partCopy;

  return;

}

//----------------------------------------------------
void AliAnalysisTaskSELc2V0bachelor::DefineSignalHistosSeparatedPerOrigin()
{
  //
  // Define analysis histograms for SNG separated for origin (from c, from b and from no-quark)
  //

  if (!fUseMCInfo) return;
  if (!fCheckOrigin) return;

  Double_t mLcPDG  = TDatabasePDG::Instance()->GetParticle(4122)->Mass();
  Double_t mK0SPDG = TDatabasePDG::Instance()->GetParticle(310)->Mass();
  Double_t mMinLambdaPDG  = TDatabasePDG::Instance()->GetParticle(2212)->Mass()+
    TDatabasePDG::Instance()->GetParticle(211)->Mass();

  TString nameHistoSgnC=" ", nameHistoSgnB=" ", nameHistoSgnNoQ=" ";
  TString titleHistoSgnC=" ", titleHistoSgnB=" ", titleHistoSgnNoQ=" ";

  // pt(Lc)
  Double_t *binLimpTLc=new Double_t[11+1]; // 11 pT(Lc) bins
  binLimpTLc[ 0]= 0.;
  binLimpTLc[ 1]= 1.;
  binLimpTLc[ 2]= 2.;
  binLimpTLc[ 3]= 3.;
  binLimpTLc[ 4]= 4.;
  binLimpTLc[ 5]= 5.;
  binLimpTLc[ 6]= 6.;
  binLimpTLc[ 7]= 8.;
  binLimpTLc[ 8]=12.;
  binLimpTLc[ 9]=17.;
  binLimpTLc[10]=25.;
  binLimpTLc[11]=35.;

  // pt(prong)
  Double_t *binLimpTprong=new Double_t[41+1]; // 41 pT(prong) bins
  binLimpTprong[ 0]= 0.0;
  binLimpTprong[ 1]= 0.1;
  binLimpTprong[ 2]= 0.2;
  binLimpTprong[ 3]= 0.3;
  binLimpTprong[ 4]= 0.4;
  binLimpTprong[ 5]= 0.5;
  binLimpTprong[ 6]= 0.6;
  binLimpTprong[ 7]= 0.7;
  binLimpTprong[ 8]= 0.8;
  binLimpTprong[ 9]= 0.9;
  binLimpTprong[10]= 1.0;
  binLimpTprong[11]= 1.2;
  binLimpTprong[12]= 1.4;
  binLimpTprong[13]= 1.6;
  binLimpTprong[14]= 1.8;
  binLimpTprong[15]= 2.0;
  binLimpTprong[16]= 2.2;
  binLimpTprong[17]= 2.4;
  binLimpTprong[18]= 2.6;
  binLimpTprong[19]= 2.8;
  binLimpTprong[20]= 3.0;
  binLimpTprong[21]= 3.5;
  binLimpTprong[22]= 4.0;
  binLimpTprong[23]= 4.5;
  binLimpTprong[24]= 5.0;
  binLimpTprong[25]= 5.5;
  binLimpTprong[26]= 6.0;
  binLimpTprong[27]= 6.5;
  binLimpTprong[28]= 7.0;
  binLimpTprong[29]= 7.5;
  binLimpTprong[30]= 8.0;
  binLimpTprong[31]= 9.0;
  binLimpTprong[32]=10.0;
  binLimpTprong[33]=11.0;
  binLimpTprong[34]=12.0;
  binLimpTprong[35]=13.0;
  binLimpTprong[36]=14.0;
  binLimpTprong[37]=15.0;
  binLimpTprong[38]=20.0;
  binLimpTprong[39]=25.0;
  binLimpTprong[40]=30.0;
  binLimpTprong[41]=35.0;

  if (fUseOnTheFlyV0) {

    nameHistoSgnC="histLcMassByK0SSgnC";
    nameHistoSgnB="histLcMassByK0SSgnB";
    nameHistoSgnNoQ="histLcMassByK0SSgnNoQ";
    titleHistoSgnC="#Lambda_{c} #leftarrow c - sgn: invariant mass (by K^{0}_{S})  vs p_{T} - MC; m_{inv}(p,K^{0}_{S}) [GeV/c^{2}]; p_{T}(#Lambda_{c}) [GeV/c]";
    titleHistoSgnB="#Lambda_{c} #leftarrow b - sgn: invariant mass (by K^{0}_{S})  vs p_{T} - MC; m_{inv}(p,K^{0}_{S}) [GeV/c^{2}]; p_{T}(#Lambda_{c}) [GeV/c]";
    titleHistoSgnNoQ="#Lambda_{c} #leftarrow no quark - sgn: invariant mass (by K^{0}_{S})  vs p_{T} - MC; m_{inv}(p,K^{0}_{S}) [GeV/c^{2}]; p_{T}(#Lambda_{c}) [GeV/c]";
    TH2F* spectrumLcMassByK0SSgnC = new TH2F(nameHistoSgnC.Data(),titleHistoSgnC.Data(),1000,mLcPDG-0.250,mLcPDG+0.250,11,binLimpTLc);
    TH2F* spectrumLcMassByK0SSgnB = new TH2F(nameHistoSgnB.Data(),titleHistoSgnB.Data(),1000,mLcPDG-0.250,mLcPDG+0.250,11,binLimpTLc);
    TH2F* spectrumLcMassByK0SSgnNoQ = new TH2F(nameHistoSgnNoQ.Data(),titleHistoSgnNoQ.Data(),1000,mLcPDG-0.250,mLcPDG+0.250,11,binLimpTLc);

    TH2F* allspectrumLcMassByK0SSgnC = (TH2F*)spectrumLcMassByK0SSgnC->Clone();
    TH2F* allspectrumLcMassByK0SSgnB = (TH2F*) spectrumLcMassByK0SSgnB->Clone();
    TH2F* allspectrumLcMassByK0SSgnNoQ = (TH2F*) spectrumLcMassByK0SSgnNoQ->Clone();
    TH2F* pidBachspectrumLcMassByK0SSgnC = (TH2F*)spectrumLcMassByK0SSgnC->Clone();
    TH2F* pidBachspectrumLcMassByK0SSgnB = (TH2F*) spectrumLcMassByK0SSgnB->Clone();
    TH2F* pidBachspectrumLcMassByK0SSgnNoQ = (TH2F*) spectrumLcMassByK0SSgnNoQ->Clone();

    fOutputAll->Add(allspectrumLcMassByK0SSgnC);
    fOutputAll->Add(allspectrumLcMassByK0SSgnB);
    fOutputAll->Add(allspectrumLcMassByK0SSgnNoQ);
    fOutputPIDBach->Add(pidBachspectrumLcMassByK0SSgnC);
    fOutputPIDBach->Add(pidBachspectrumLcMassByK0SSgnB);
    fOutputPIDBach->Add(pidBachspectrumLcMassByK0SSgnNoQ);

    nameHistoSgnC="histptK0SSgnC";
    nameHistoSgnB="histptK0SSgnB";
    nameHistoSgnNoQ="histptK0SSgnNoQ";
    titleHistoSgnC="p_{T}(K^{0}_{S}) vs p_{T}(#Lambda_{c}); p_{T}(#Lambda_{c}) [GeV/c]; p_{T}(K^{0}_{S}) [GeV/c]; Entries";
    titleHistoSgnB="p_{T}(K^{0}_{S}) vs p_{T}(#Lambda_{c}); p_{T}(#Lambda_{c}) [GeV/c]; p_{T}(K^{0}_{S}) [GeV/c]; Entries";
    titleHistoSgnNoQ="p_{T}(K^{0}_{S}) vs p_{T}(#Lambda_{c}); p_{T}(#Lambda_{c}) [GeV/c]; p_{T}(K^{0}_{S}) [GeV/c]; Entries";
    TH2F* ptK0SSgnC = new TH2F(nameHistoSgnC.Data(),titleHistoSgnC.Data(),11,binLimpTLc,41,binLimpTprong);
    TH2F* ptK0SSgnB = new TH2F(nameHistoSgnB.Data(),titleHistoSgnB.Data(),11,binLimpTLc,41,binLimpTprong);
    TH2F* ptK0SSgnNoQ = new TH2F(nameHistoSgnNoQ.Data(),titleHistoSgnNoQ.Data(),11,binLimpTLc,41,binLimpTprong);

    nameHistoSgnC="histptPSgnC";
    nameHistoSgnB="histptPSgnB";
    nameHistoSgnNoQ="histptPSgnNoQ";
    titleHistoSgnC="p_{T}(p) vs p_{T}(#Lambda_{c}); p_{T}(#Lambda_{c}) [GeV/c]; p_{T}(p) [GeV/c]; Entries";
    titleHistoSgnB="p_{T}(p) vs p_{T}(#Lambda_{c}); p_{T}(#Lambda_{c}) [GeV/c]; p_{T}(p) [GeV/c]; Entries";
    titleHistoSgnNoQ="p_{T}(p) vs p_{T}(#Lambda_{c}); p_{T}(#Lambda_{c}) [GeV/c]; p_{T}(p) [GeV/c]; Entries";
    TH2F* ptPSgnC = new TH2F(nameHistoSgnC.Data(),titleHistoSgnC.Data(),11,binLimpTLc,41,binLimpTprong);
    TH2F* ptPSgnB = new TH2F(nameHistoSgnB.Data(),titleHistoSgnB.Data(),11,binLimpTLc,41,binLimpTprong);
    TH2F* ptPSgnNoQ = new TH2F(nameHistoSgnNoQ.Data(),titleHistoSgnNoQ.Data(),11,binLimpTLc,41,binLimpTprong);

    nameHistoSgnC="histptPipSgnC";
    nameHistoSgnB="histptPipSgnB";
    nameHistoSgnNoQ="histptPipSgnNoQ";
    titleHistoSgnC="p_{T}(#pi^{+}) vs p_{T}(#Lambda_{c}); p_{T}(#Lambda_{c}) [GeV/c]; p_{T}(#pi^{+}) [GeV/c]; Entries";
    titleHistoSgnB="p_{T}(#pi^{+}) vs p_{T}(#Lambda_{c}); p_{T}(#Lambda_{c}) [GeV/c]; p_{T}(#pi^{+}) [GeV/c]; Entries";
    titleHistoSgnNoQ="p_{T}(#pi^{+}) vs p_{T}(#Lambda_{c}); p_{T}(#Lambda_{c}) [GeV/c]; p_{T}(#pi^{+}) [GeV/c]; Entries";
    TH2F* ptPiPSgnC = new TH2F(nameHistoSgnC.Data(),titleHistoSgnC.Data(),11,binLimpTLc,41,binLimpTprong);
    TH2F* ptPiPSgnB = new TH2F(nameHistoSgnB.Data(),titleHistoSgnB.Data(),11,binLimpTLc,41,binLimpTprong);
    TH2F* ptPiPSgnNoQ = new TH2F(nameHistoSgnNoQ.Data(),titleHistoSgnB.Data(),11,binLimpTLc,41,binLimpTprong);

    nameHistoSgnC="histptPimSgnC";
    nameHistoSgnB="histptPimSgnB";
    nameHistoSgnNoQ="histptPimSgnNoQ";
    titleHistoSgnC="p_{T}(#pi^{-}) vs p_{T}(#Lambda_{c}); p_{T}(#Lambda_{c}) [GeV/c]; p_{T}(#pi^{-}) [GeV/c]; Entries";
    titleHistoSgnB="p_{T}(#pi^{-}) vs p_{T}(#Lambda_{c}); p_{T}(#Lambda_{c}) [GeV/c]; p_{T}(#pi^{-}) [GeV/c]; Entries";
    titleHistoSgnNoQ="p_{T}(#pi^{-}) vs p_{T}(#Lambda_{c}); p_{T}(#Lambda_{c}) [GeV/c]; p_{T}(#pi^{-}) [GeV/c]; Entries";
    TH2F* ptPiMSgnC = new TH2F(nameHistoSgnC.Data(),titleHistoSgnC.Data(),11,binLimpTLc,41,binLimpTprong);
    TH2F* ptPiMSgnB = new TH2F(nameHistoSgnB.Data(),titleHistoSgnB.Data(),11,binLimpTLc,41,binLimpTprong);
    TH2F* ptPiMSgnNoQ = new TH2F(nameHistoSgnNoQ.Data(),titleHistoSgnB.Data(),11,binLimpTLc,41,binLimpTprong);

    nameHistoSgnC="histD0K0SSgnC";
    nameHistoSgnB="histD0K0SSgnB";
    nameHistoSgnNoQ="histD0K0SSgnNoQ";
    titleHistoSgnC="d_{0}(K^{0}_{S}) vs p_{T}(#Lambda_{c}); p_{T}(#Lambda_{c}) [GeV/c]; d_{0}(K^{0}_{S}) [#sigmas]; Entries";
    titleHistoSgnB="d_{0}(K^{0}_{S}) vs p_{T}(#Lambda_{c}); p_{T}(#Lambda_{c}) [GeV/c]; d_{0}(K^{0}_{S}) [#sigmas]; Entries";
    titleHistoSgnNoQ="d_{0}(K^{0}_{S}) vs p_{T}(#Lambda_{c}); p_{T}(#Lambda_{c}) [GeV/c]; d_{0}(K^{0}_{S}) [#sigmas]; Entries";
    TH2F* d0K0SSgnC = new TH2F(nameHistoSgnC.Data(),titleHistoSgnC.Data(),11,binLimpTLc,1000,-1.,1.);
    TH2F* d0K0SSgnB = new TH2F(nameHistoSgnB.Data(),titleHistoSgnB.Data(),11,binLimpTLc,1000,-1.,1.);
    TH2F* d0K0SSgnNoQ = new TH2F(nameHistoSgnNoQ.Data(),titleHistoSgnNoQ.Data(),11,binLimpTLc,1000,-1.,1.);

    nameHistoSgnC="histD0PSgnC";
    nameHistoSgnB="histD0PSgnB";
    nameHistoSgnNoQ="histD0PSgnNoQ";
    titleHistoSgnC="d_{0}(p) vs p_{T}(#Lambda_{c}); p_{T}(#Lambda_{c}) [GeV/c]; d_{0}(p) [cm]; Entries";
    titleHistoSgnB="d_{0}(p) vs p_{T}(#Lambda_{c}); p_{T}(#Lambda_{c}) [GeV/c]; d_{0}(p) [cm]; Entries";
    titleHistoSgnNoQ="d_{0}(p) vs p_{T}(#Lambda_{c}); p_{T}(#Lambda_{c}) [GeV/c]; d_{0}(p) [cm]; Entries";
    TH2F* d0PSgnC = new TH2F(nameHistoSgnC.Data(),titleHistoSgnC.Data(),11,binLimpTLc,1000,-1.,1.);
    TH2F* d0PSgnB = new TH2F(nameHistoSgnB.Data(),titleHistoSgnB.Data(),11,binLimpTLc,1000,-1.,1.);
    TH2F* d0PSgnNoQ = new TH2F(nameHistoSgnNoQ.Data(),titleHistoSgnNoQ.Data(),11,binLimpTLc,1000,-1.,1.);

    nameHistoSgnC="histCosPAK0SSgnC";
    nameHistoSgnB="histCosPAK0SSgnB";
    nameHistoSgnNoQ="histCosPAK0SSgnNoQ";
    titleHistoSgnC="K^{0}_{S} cosine of pointing angle wrt primary vertex vs p_{T}(#Lambda_{c}); p_{T}(#Lambda_{c}) [GeV/c]; cosine; Entries";
    titleHistoSgnB="K^{0}_{S} cosine of pointing angle wrt primary vertex vs p_{T}(#Lambda_{c}); p_{T}(#Lambda_{c}) [GeV/c]; cosine; Entries";
    titleHistoSgnNoQ="K^{0}_{S} cosine of pointing angle wrt primary vertex vs p_{T}(#Lambda_{c}); p_{T}(#Lambda_{c}) [GeV/c]; cosine; Entries";
    TH2F *cosPAK0SSgnC = new TH2F(nameHistoSgnC.Data(),titleHistoSgnC.Data(),41,binLimpTprong,100,0.99,1.);
    TH2F *cosPAK0SSgnB = new TH2F(nameHistoSgnB.Data(),titleHistoSgnB.Data(),41,binLimpTprong,100,0.99,1.);
    TH2F *cosPAK0SSgnNoQ = new TH2F(nameHistoSgnNoQ.Data(),titleHistoSgnNoQ.Data(),41,binLimpTprong,100,0.99,1.);

    nameHistoSgnC="histCosThetaProtonCMSSgnC";
    nameHistoSgnB="histCosThetaProtonCMSSgnB";
    nameHistoSgnNoQ="histCosThetaProtonCMSSgnNoQ";
    titleHistoSgnC="cosien of proton emission angle in Lc rest frame; p_{T}(#Lambda_{c}) [GeV/c]; cosine; Entries";
    titleHistoSgnB="cosien of proton emission angle in Lc rest frame; p_{T}(#Lambda_{c}) [GeV/c]; cosine; Entries";
    titleHistoSgnNoQ="cosien of proton emission angle in Lc rest frame; p_{T}(#Lambda_{c}) [GeV/c]; cosine; Entries";
    TH2F *cosThePrSgnC = new TH2F(nameHistoSgnC.Data(),titleHistoSgnC.Data(),41,binLimpTprong,100,-1.,1.);
    TH2F *cosThePrSgnB = new TH2F(nameHistoSgnB.Data(),titleHistoSgnB.Data(),41,binLimpTprong,100,-1.,1.);
    TH2F *cosThePrSgnNoQ = new TH2F(nameHistoSgnNoQ.Data(),titleHistoSgnNoQ.Data(),41,binLimpTprong,100,-1.,1.);

    nameHistoSgnC="histResignedD0SgnC";
    nameHistoSgnB="histResignedD0SgnB";
    nameHistoSgnNoQ="histResignedD0SgnNoQ";
    titleHistoSgnC="Proton d0 with different sign convention; p_{T}(#Lambda_{c}) [GeV/c]; d0 [cm]; Entries";
    titleHistoSgnB="Proton d0 with different sign convention; p_{T}(#Lambda_{c}) [GeV/c]; d0 [cm]; Entries";
    titleHistoSgnNoQ="Proton d0 with different sign convention; p_{T}(#Lambda_{c}) [GeV/c]; d0 [cm]; Entries";
    TH2F *resignedD0SgnC = new TH2F(nameHistoSgnC.Data(),titleHistoSgnC.Data(),41,binLimpTprong,100,-0.1,0.1);
    TH2F *resignedD0SgnB = new TH2F(nameHistoSgnB.Data(),titleHistoSgnB.Data(),41,binLimpTprong,100,-0.1,0.1);
    TH2F *resignedD0SgnNoQ = new TH2F(nameHistoSgnNoQ.Data(),titleHistoSgnNoQ.Data(),41,binLimpTprong,100,-0.1,0.1);

    TH2F* allptK0SSgnC = (TH2F*)ptK0SSgnC->Clone();
    TH2F* allptK0SSgnB = (TH2F*)ptK0SSgnB->Clone();
    TH2F* allptK0SSgnNoQ = (TH2F*)ptK0SSgnNoQ->Clone();
    TH2F* allptPSgnC = (TH2F*)ptPSgnC->Clone();
    TH2F* allptPSgnB = (TH2F*)ptPSgnB->Clone();
    TH2F* allptPSgnNoQ = (TH2F*)ptPSgnNoQ->Clone();
    TH2F* allptPiPSgnC = (TH2F*)ptPiPSgnC->Clone();
    TH2F* allptPiPSgnB = (TH2F*)ptPiPSgnB->Clone();
    TH2F* allptPiPSgnNoQ = (TH2F*)ptPiPSgnNoQ->Clone();
    TH2F* allptPiMSgnC = (TH2F*)ptPiMSgnC->Clone();
    TH2F* allptPiMSgnB = (TH2F*)ptPiMSgnB->Clone();
    TH2F* allptPiMSgnNoQ = (TH2F*)ptPiMSgnNoQ->Clone();
    TH2F* alld0K0SSgnC = (TH2F*)d0K0SSgnC->Clone();
    TH2F* alld0K0SSgnB = (TH2F*)d0K0SSgnB->Clone();
    TH2F* alld0K0SSgnNoQ = (TH2F*)d0K0SSgnNoQ->Clone();
    TH2F* alld0PSgnC = (TH2F*)d0PSgnC->Clone();
    TH2F* alld0PSgnB = (TH2F*)d0PSgnB->Clone();
    TH2F* alld0PSgnNoQ = (TH2F*)d0PSgnNoQ->Clone();
    TH2F* allcosPAK0SSgnC = (TH2F*)cosPAK0SSgnC->Clone();
    TH2F* allcosPAK0SSgnB = (TH2F*)cosPAK0SSgnB->Clone();
    TH2F* allcosPAK0SSgnNoQ = (TH2F*)cosPAK0SSgnNoQ->Clone();
    TH2F* allcosThePrSgnC = (TH2F*)cosThePrSgnC->Clone();
    TH2F* allcosThePrSgnB = (TH2F*)cosThePrSgnB->Clone();
    TH2F* allcosThePrSgnNoQ = (TH2F*)cosThePrSgnNoQ->Clone();
    TH2F* allresignedD0SgnC = (TH2F*)resignedD0SgnC->Clone();
    TH2F* allresignedD0SgnB = (TH2F*)resignedD0SgnB->Clone();
    TH2F* allresignedD0SgnNoQ = (TH2F*)resignedD0SgnNoQ->Clone();

    TH2F* pidptK0SSgnC = (TH2F*)ptK0SSgnC->Clone();
    TH2F* pidptK0SSgnB = (TH2F*)ptK0SSgnB->Clone();
    TH2F* pidptK0SSgnNoQ = (TH2F*)ptK0SSgnNoQ->Clone();
    TH2F* pidptPSgnC = (TH2F*)ptPSgnC->Clone();
    TH2F* pidptPSgnB = (TH2F*)ptPSgnB->Clone();
    TH2F* pidptPSgnNoQ = (TH2F*)ptPSgnNoQ->Clone();
    TH2F* pidptPiPSgnC = (TH2F*)ptPiPSgnC->Clone();
    TH2F* pidptPiPSgnB = (TH2F*)ptPiPSgnB->Clone();
    TH2F* pidptPiPSgnNoQ = (TH2F*)ptPiPSgnNoQ->Clone();
    TH2F* pidptPiMSgnC = (TH2F*)ptPiMSgnC->Clone();
    TH2F* pidptPiMSgnB = (TH2F*)ptPiMSgnB->Clone();
    TH2F* pidptPiMSgnNoQ = (TH2F*)ptPiMSgnNoQ->Clone();
    TH2F* pidd0K0SSgnC = (TH2F*)d0K0SSgnC->Clone();
    TH2F* pidd0K0SSgnB = (TH2F*)d0K0SSgnB->Clone();
    TH2F* pidd0K0SSgnNoQ = (TH2F*)d0K0SSgnNoQ->Clone();
    TH2F* pidd0PSgnC = (TH2F*)d0PSgnC->Clone();
    TH2F* pidd0PSgnB = (TH2F*)d0PSgnB->Clone();
    TH2F* pidd0PSgnNoQ = (TH2F*)d0PSgnNoQ->Clone();
    TH2F* pidcosPAK0SSgnC = (TH2F*)cosPAK0SSgnC->Clone();
    TH2F* pidcosPAK0SSgnB = (TH2F*)cosPAK0SSgnB->Clone();
    TH2F* pidcosPAK0SSgnNoQ = (TH2F*)cosPAK0SSgnNoQ->Clone();
    TH2F* pidcosThePrSgnC = (TH2F*)cosThePrSgnC->Clone();
    TH2F* pidcosThePrSgnB = (TH2F*)cosThePrSgnB->Clone();
    TH2F* pidcosThePrSgnNoQ = (TH2F*)cosThePrSgnNoQ->Clone();
    TH2F* pidresignedD0SgnC = (TH2F*)resignedD0SgnC->Clone();
    TH2F* pidresignedD0SgnB = (TH2F*)resignedD0SgnB->Clone();
    TH2F* pidresignedD0SgnNoQ = (TH2F*)resignedD0SgnNoQ->Clone();

    fOutputAll->Add(allptK0SSgnC);
    fOutputAll->Add(allptK0SSgnB);
    fOutputAll->Add(allptK0SSgnNoQ);
    fOutputAll->Add(allptPSgnC);
    fOutputAll->Add(allptPSgnB);
    fOutputAll->Add(allptPSgnNoQ);
    fOutputAll->Add(allptPiPSgnC);
    fOutputAll->Add(allptPiPSgnB);
    fOutputAll->Add(allptPiPSgnNoQ);
    fOutputAll->Add(allptPiMSgnC);
    fOutputAll->Add(allptPiMSgnB);
    fOutputAll->Add(allptPiMSgnNoQ);
    fOutputAll->Add(alld0K0SSgnC);
    fOutputAll->Add(alld0K0SSgnB);
    fOutputAll->Add(alld0K0SSgnNoQ);
    fOutputAll->Add(alld0PSgnC);
    fOutputAll->Add(alld0PSgnB);
    fOutputAll->Add(alld0PSgnNoQ);
    fOutputAll->Add(allcosPAK0SSgnC);
    fOutputAll->Add(allcosPAK0SSgnB);
    fOutputAll->Add(allcosPAK0SSgnNoQ);
    fOutputAll->Add(allcosThePrSgnC);
    fOutputAll->Add(allcosThePrSgnB);
    fOutputAll->Add(allcosThePrSgnNoQ);
    fOutputAll->Add(allresignedD0SgnC);
    fOutputAll->Add(allresignedD0SgnB);
    fOutputAll->Add(allresignedD0SgnNoQ);

    fOutputPIDBach->Add(pidptK0SSgnC);
    fOutputPIDBach->Add(pidptK0SSgnB);
    fOutputPIDBach->Add(pidptK0SSgnNoQ);
    fOutputPIDBach->Add(pidptPSgnC);
    fOutputPIDBach->Add(pidptPSgnB);
    fOutputPIDBach->Add(pidptPSgnNoQ);
    fOutputPIDBach->Add(pidptPiPSgnC);
    fOutputPIDBach->Add(pidptPiPSgnB);
    fOutputPIDBach->Add(pidptPiPSgnNoQ);
    fOutputPIDBach->Add(pidptPiMSgnC);
    fOutputPIDBach->Add(pidptPiMSgnB);
    fOutputPIDBach->Add(pidptPiMSgnNoQ);
    fOutputPIDBach->Add(pidd0K0SSgnC);
    fOutputPIDBach->Add(pidd0K0SSgnB);
    fOutputPIDBach->Add(pidd0K0SSgnNoQ);
    fOutputPIDBach->Add(pidd0PSgnC);
    fOutputPIDBach->Add(pidd0PSgnB);
    fOutputPIDBach->Add(pidd0PSgnNoQ);
    fOutputPIDBach->Add(pidcosPAK0SSgnC);
    fOutputPIDBach->Add(pidcosPAK0SSgnB);
    fOutputPIDBach->Add(pidcosPAK0SSgnNoQ);
    fOutputPIDBach->Add(pidcosThePrSgnC);
    fOutputPIDBach->Add(pidcosThePrSgnB);
    fOutputPIDBach->Add(pidcosThePrSgnNoQ);
    fOutputPIDBach->Add(pidresignedD0SgnC);
    fOutputPIDBach->Add(pidresignedD0SgnB);
    fOutputPIDBach->Add(pidresignedD0SgnNoQ);

  }

  nameHistoSgnC="histLcMassByK0SOfflineSgnC";
  nameHistoSgnB="histLcMassByK0SOfflineSgnB";
  nameHistoSgnNoQ="histLcMassByK0SOfflineSgnNoQ";
  titleHistoSgnC="#Lambda_{c} #leftarrow c - sgn: invariant mass (by K^{0}_{S})  vs p_{T} - MC; m_{inv}(p,K^{0}_{S}) [GeV/c^{2}]; p_{T}(#Lambda_{c}) [GeV/c]";
  titleHistoSgnB="#Lambda_{c} #leftarrow b - sgn: invariant mass (by K^{0}_{S})  vs p_{T} - MC; m_{inv}(p,K^{0}_{S}) [GeV/c^{2}]; p_{T}(#Lambda_{c}) [GeV/c]";
  titleHistoSgnNoQ="#Lambda_{c} #leftarrow no quark - sgn: invariant mass (by K^{0}_{S})  vs p_{T} - MC; m_{inv}(p,K^{0}_{S}) [GeV/c^{2}]; p_{T}(#Lambda_{c}) [GeV/c]";
  TH2F* spectrumLcMassOfflineByK0SSgnC = new TH2F(nameHistoSgnC.Data(),titleHistoSgnC.Data(),1000,mLcPDG-0.250,mLcPDG+0.250,11,binLimpTLc);
  TH2F* spectrumLcMassOfflineByK0SSgnB = new TH2F(nameHistoSgnB.Data(),titleHistoSgnB.Data(),1000,mLcPDG-0.250,mLcPDG+0.250,11,binLimpTLc);
  TH2F* spectrumLcMassOfflineByK0SSgnNoQ = new TH2F(nameHistoSgnNoQ.Data(),titleHistoSgnNoQ.Data(),1000,mLcPDG-0.250,mLcPDG+0.250,11,binLimpTLc);

  TH2F* allspectrumLcMassOfflineByK0SSgnC = (TH2F*)spectrumLcMassOfflineByK0SSgnC->Clone();
  TH2F* allspectrumLcMassOfflineByK0SSgnB = (TH2F*) spectrumLcMassOfflineByK0SSgnB->Clone();
  TH2F* allspectrumLcMassOfflineByK0SSgnNoQ = (TH2F*) spectrumLcMassOfflineByK0SSgnNoQ->Clone();
  TH2F* pidBachspectrumLcMassOfflineByK0SSgnC = (TH2F*)spectrumLcMassOfflineByK0SSgnC->Clone();
  TH2F* pidBachspectrumLcMassOfflineByK0SSgnB = (TH2F*) spectrumLcMassOfflineByK0SSgnB->Clone();
  TH2F* pidBachspectrumLcMassOfflineByK0SSgnNoQ = (TH2F*) spectrumLcMassOfflineByK0SSgnNoQ->Clone();
  fOutputAll->Add(allspectrumLcMassOfflineByK0SSgnC);
  fOutputAll->Add(allspectrumLcMassOfflineByK0SSgnB);
  fOutputAll->Add(allspectrumLcMassOfflineByK0SSgnNoQ);
  fOutputPIDBach->Add(pidBachspectrumLcMassOfflineByK0SSgnC);
  fOutputPIDBach->Add(pidBachspectrumLcMassOfflineByK0SSgnB);
  fOutputPIDBach->Add(pidBachspectrumLcMassOfflineByK0SSgnNoQ);

  nameHistoSgnC="histptK0SOfflineSgnC";
  nameHistoSgnB="histptK0SOfflineSgnB";
  nameHistoSgnNoQ="histptK0SOfflineSgnNoQ";
  titleHistoSgnC="p_{T}(K^{0}_{S}) vs p_{T}(#Lambda_{c}); p_{T}(#Lambda_{c}) [GeV/c]; p_{T}(K^{0}_{S}) [GeV/c]; Entries";
  titleHistoSgnB="p_{T}(K^{0}_{S}) vs p_{T}(#Lambda_{c}); p_{T}(#Lambda_{c}) [GeV/c]; p_{T}(K^{0}_{S}) [GeV/c]; Entries";
  titleHistoSgnNoQ="p_{T}(K^{0}_{S}) vs p_{T}(#Lambda_{c}); p_{T}(#Lambda_{c}) [GeV/c]; p_{T}(K^{0}_{S}) [GeV/c]; Entries";
  TH2F* ptK0SOfflineSgnC = new TH2F(nameHistoSgnC.Data(),titleHistoSgnC.Data(),11,binLimpTLc,41,binLimpTprong);
  TH2F* ptK0SOfflineSgnB = new TH2F(nameHistoSgnB.Data(),titleHistoSgnB.Data(),11,binLimpTLc,41,binLimpTprong);
  TH2F* ptK0SOfflineSgnNoQ = new TH2F(nameHistoSgnNoQ.Data(),titleHistoSgnNoQ.Data(),11,binLimpTLc,41,binLimpTprong);

  nameHistoSgnC="histptPOfflineSgnC";
  nameHistoSgnB="histptPOfflineSgnB";
  nameHistoSgnNoQ="histptPOfflineSgnNoQ";
  titleHistoSgnC="p_{T}(p) vs p_{T}(#Lambda_{c}); p_{T}(#Lambda_{c}) [GeV/c]; p_{T}(p) [GeV/c]; Entries";
  titleHistoSgnB="p_{T}(p) vs p_{T}(#Lambda_{c}); p_{T}(#Lambda_{c}) [GeV/c]; p_{T}(p) [GeV/c]; Entries";
  titleHistoSgnNoQ="p_{T}(p) vs p_{T}(#Lambda_{c}); p_{T}(#Lambda_{c}) [GeV/c]; p_{T}(p) [GeV/c]; Entries";
  TH2F* ptPOfflineSgnC = new TH2F(nameHistoSgnC.Data(),titleHistoSgnC.Data(),11,binLimpTLc,41,binLimpTprong);
  TH2F* ptPOfflineSgnB = new TH2F(nameHistoSgnB.Data(),titleHistoSgnB.Data(),11,binLimpTLc,41,binLimpTprong);
  TH2F* ptPOfflineSgnNoQ = new TH2F(nameHistoSgnNoQ.Data(),titleHistoSgnNoQ.Data(),11,binLimpTLc,41,binLimpTprong);

  nameHistoSgnC="histptPipOfflineSgnC";
  nameHistoSgnB="histptPipOfflineSgnB";
  nameHistoSgnNoQ="histptPipOfflineSgnNoQ";
  titleHistoSgnC="p_{T}(#pi^{+}) vs p_{T}(#Lambda_{c}); p_{T}(#Lambda_{c}) [GeV/c]; p_{T}(#pi^{+}) [GeV/c]; Entries";
  titleHistoSgnB="p_{T}(#pi^{+}) vs p_{T}(#Lambda_{c}); p_{T}(#Lambda_{c}) [GeV/c]; p_{T}(#pi^{+}) [GeV/c]; Entries";
  titleHistoSgnNoQ="p_{T}(#pi^{+}) vs p_{T}(#Lambda_{c}); p_{T}(#Lambda_{c}) [GeV/c]; p_{T}(#pi^{+}) [GeV/c]; Entries";
  TH2F* ptPiPOfflineSgnC = new TH2F(nameHistoSgnC.Data(),titleHistoSgnC.Data(),11,binLimpTLc,41,binLimpTprong);
  TH2F* ptPiPOfflineSgnB = new TH2F(nameHistoSgnB.Data(),titleHistoSgnB.Data(),11,binLimpTLc,41,binLimpTprong);
  TH2F* ptPiPOfflineSgnNoQ = new TH2F(nameHistoSgnNoQ.Data(),titleHistoSgnB.Data(),11,binLimpTLc,41,binLimpTprong);

  nameHistoSgnC="histptPimOfflineSgnC";
  nameHistoSgnB="histptPimOfflineSgnB";
  nameHistoSgnNoQ="histptPimOfflineSgnNoQ";
  titleHistoSgnC="p_{T}(#pi^{-}) vs p_{T}(#Lambda_{c}); p_{T}(#Lambda_{c}) [GeV/c]; p_{T}(#pi^{-}) [GeV/c]; Entries";
  titleHistoSgnB="p_{T}(#pi^{-}) vs p_{T}(#Lambda_{c}); p_{T}(#Lambda_{c}) [GeV/c]; p_{T}(#pi^{-}) [GeV/c]; Entries";
  titleHistoSgnNoQ="p_{T}(#pi^{-}) vs p_{T}(#Lambda_{c}); p_{T}(#Lambda_{c}) [GeV/c]; p_{T}(#pi^{-}) [GeV/c]; Entries";
  TH2F* ptPiMOfflineSgnC = new TH2F(nameHistoSgnC.Data(),titleHistoSgnC.Data(),11,binLimpTLc,41,binLimpTprong);
  TH2F* ptPiMOfflineSgnB = new TH2F(nameHistoSgnB.Data(),titleHistoSgnB.Data(),11,binLimpTLc,41,binLimpTprong);
  TH2F* ptPiMOfflineSgnNoQ = new TH2F(nameHistoSgnNoQ.Data(),titleHistoSgnB.Data(),11,binLimpTLc,41,binLimpTprong);

  nameHistoSgnC="histD0K0SOfflineSgnC";
  nameHistoSgnB="histD0K0SOfflineSgnB";
  nameHistoSgnNoQ="histD0K0SOfflineSgnNoQ";
  titleHistoSgnC="d_{0}(K^{0}_{S}) vs p_{T}(#Lambda_{c}); p_{T}(#Lambda_{c}) [GeV/c]; d_{0}(K^{0}_{S}) [#sigmas]; Entries";
  titleHistoSgnB="d_{0}(K^{0}_{S}) vs p_{T}(#Lambda_{c}); p_{T}(#Lambda_{c}) [GeV/c]; d_{0}(K^{0}_{S}) [#sigmas]; Entries";
  titleHistoSgnNoQ="d_{0}(K^{0}_{S}) vs p_{T}(#Lambda_{c}); p_{T}(#Lambda_{c}) [GeV/c]; d_{0}(K^{0}_{S}) [#sigmas]; Entries";
  TH2F* d0K0SOfflineSgnC = new TH2F(nameHistoSgnC.Data(),titleHistoSgnC.Data(),11,binLimpTLc,1000,-1.,1.);
  TH2F* d0K0SOfflineSgnB = new TH2F(nameHistoSgnB.Data(),titleHistoSgnB.Data(),11,binLimpTLc,1000,-1.,1.);
  TH2F* d0K0SOfflineSgnNoQ = new TH2F(nameHistoSgnNoQ.Data(),titleHistoSgnNoQ.Data(),11,binLimpTLc,1000,-1.,1.);

  nameHistoSgnC="histD0POfflineSgnC";
  nameHistoSgnB="histD0POfflineSgnB";
  nameHistoSgnNoQ="histD0POfflineSgnNoQ";
  titleHistoSgnC="d_{0}(p) vs p_{T}(#Lambda_{c}); p_{T}(#Lambda_{c}) [GeV/c]; d_{0}(p) [cm]; Entries";
  titleHistoSgnB="d_{0}(p) vs p_{T}(#Lambda_{c}); p_{T}(#Lambda_{c}) [GeV/c]; d_{0}(p) [cm]; Entries";
  titleHistoSgnNoQ="d_{0}(p) vs p_{T}(#Lambda_{c}); p_{T}(#Lambda_{c}) [GeV/c]; d_{0}(p) [cm]; Entries";
  TH2F* d0POfflineSgnC = new TH2F(nameHistoSgnC.Data(),titleHistoSgnC.Data(),11,binLimpTLc,1000,-1.,1.);
  TH2F* d0POfflineSgnB = new TH2F(nameHistoSgnB.Data(),titleHistoSgnB.Data(),11,binLimpTLc,1000,-1.,1.);
  TH2F* d0POfflineSgnNoQ = new TH2F(nameHistoSgnNoQ.Data(),titleHistoSgnNoQ.Data(),11,binLimpTLc,1000,-1.,1.);

  nameHistoSgnC="histCosPAK0SOfflineSgnC";
  nameHistoSgnB="histCosPAK0SOfflineSgnB";
  nameHistoSgnNoQ="histCosPAK0SOfflineSgnNoQ";
  titleHistoSgnC="K^{0}_{S} cosine of pointing angle wrt primary vertex vs p_{T}(#Lambda_{c}); p_{T}(#Lambda_{c}) [GeV/c]; cosine; Entries";
  titleHistoSgnB="K^{0}_{S} cosine of pointing angle wrt primary vertex vs p_{T}(#Lambda_{c}); p_{T}(#Lambda_{c}) [GeV/c]; cosine; Entries";
  titleHistoSgnNoQ="K^{0}_{S} cosine of pointing angle wrt primary vertex vs p_{T}(#Lambda_{c}); p_{T}(#Lambda_{c}) [GeV/c]; cosine; Entries";
  TH2F *cosPAK0SOfflineSgnC = new TH2F(nameHistoSgnC.Data(),titleHistoSgnC.Data(),41,binLimpTprong,100,0.99,1.);
  TH2F *cosPAK0SOfflineSgnB = new TH2F(nameHistoSgnB.Data(),titleHistoSgnB.Data(),41,binLimpTprong,100,0.99,1.);
  TH2F *cosPAK0SOfflineSgnNoQ = new TH2F(nameHistoSgnNoQ.Data(),titleHistoSgnNoQ.Data(),41,binLimpTprong,100,0.99,1.);

  nameHistoSgnC="histCosThetaProtonCMSOfflineSgnC";
  nameHistoSgnB="histCosThetaProtonCMSOfflineSgnB";
  nameHistoSgnNoQ="histCosThetaProtonCMSOfflineSgnNoQ";
  titleHistoSgnC="cosien of proton emission angle in Lc rest frame; p_{T}(#Lambda_{c}) [GeV/c]; cosine; Entries";
  titleHistoSgnB="cosien of proton emission angle in Lc rest frame; p_{T}(#Lambda_{c}) [GeV/c]; cosine; Entries";
  titleHistoSgnNoQ="cosien of proton emission angle in Lc rest frame; p_{T}(#Lambda_{c}) [GeV/c]; cosine; Entries";
  TH2F *cosThePrOfflineSgnC = new TH2F(nameHistoSgnC.Data(),titleHistoSgnC.Data(),41,binLimpTprong,100,-1.,1.);
  TH2F *cosThePrOfflineSgnB = new TH2F(nameHistoSgnB.Data(),titleHistoSgnB.Data(),41,binLimpTprong,100,-1.,1.);
  TH2F *cosThePrOfflineSgnNoQ = new TH2F(nameHistoSgnNoQ.Data(),titleHistoSgnNoQ.Data(),41,binLimpTprong,100,-1.,1.);

  nameHistoSgnC="histResignedD0OfflineSgnC";
  nameHistoSgnB="histResignedD0OfflineSgnB";
  nameHistoSgnNoQ="histResignedD0OfflineSgnNoQ";
  titleHistoSgnC="Proton d0 with different sign convention; p_{T}(#Lambda_{c}) [GeV/c]; d0 [cm]; Entries";
  titleHistoSgnB="Proton d0 with different sign convention; p_{T}(#Lambda_{c}) [GeV/c]; d0 [cm]; Entries";
  titleHistoSgnNoQ="Proton d0 with different sign convention; p_{T}(#Lambda_{c}) [GeV/c]; d0 [cm]; Entries";
  TH2F *resignedD0OfflineSgnC = new TH2F(nameHistoSgnC.Data(),titleHistoSgnC.Data(),41,binLimpTprong,100,-0.1,0.1);
  TH2F *resignedD0OfflineSgnB = new TH2F(nameHistoSgnB.Data(),titleHistoSgnB.Data(),41,binLimpTprong,100,-0.1,0.1);
  TH2F *resignedD0OfflineSgnNoQ = new TH2F(nameHistoSgnNoQ.Data(),titleHistoSgnNoQ.Data(),41,binLimpTprong,100,-0.1,0.1);

  TH2F* allptK0SOfflineSgnC = (TH2F*)ptK0SOfflineSgnC->Clone();
  TH2F* allptK0SOfflineSgnB = (TH2F*)ptK0SOfflineSgnB->Clone();
  TH2F* allptK0SOfflineSgnNoQ = (TH2F*)ptK0SOfflineSgnNoQ->Clone();
  TH2F* allptPOfflineSgnC = (TH2F*)ptPOfflineSgnC->Clone();
  TH2F* allptPOfflineSgnB = (TH2F*)ptPOfflineSgnB->Clone();
  TH2F* allptPOfflineSgnNoQ = (TH2F*)ptPOfflineSgnNoQ->Clone();
  TH2F* allptPiPOfflineSgnC = (TH2F*)ptPiPOfflineSgnC->Clone();
  TH2F* allptPiPOfflineSgnB = (TH2F*)ptPiPOfflineSgnB->Clone();
  TH2F* allptPiPOfflineSgnNoQ = (TH2F*)ptPiPOfflineSgnNoQ->Clone();
  TH2F* allptPiMOfflineSgnC = (TH2F*)ptPiMOfflineSgnC->Clone();
  TH2F* allptPiMOfflineSgnB = (TH2F*)ptPiMOfflineSgnB->Clone();
  TH2F* allptPiMOfflineSgnNoQ = (TH2F*)ptPiMOfflineSgnNoQ->Clone();
  TH2F* alld0K0SOfflineSgnC = (TH2F*)d0K0SOfflineSgnC->Clone();
  TH2F* alld0K0SOfflineSgnB = (TH2F*)d0K0SOfflineSgnB->Clone();
  TH2F* alld0K0SOfflineSgnNoQ = (TH2F*)d0K0SOfflineSgnNoQ->Clone();
  TH2F* alld0POfflineSgnC = (TH2F*)d0POfflineSgnC->Clone();
  TH2F* alld0POfflineSgnB = (TH2F*)d0POfflineSgnB->Clone();
  TH2F* alld0POfflineSgnNoQ = (TH2F*)d0POfflineSgnNoQ->Clone();
  TH2F* allcosPAK0SOfflineSgnC = (TH2F*)cosPAK0SOfflineSgnC->Clone();
  TH2F* allcosPAK0SOfflineSgnB = (TH2F*)cosPAK0SOfflineSgnB->Clone();
  TH2F* allcosPAK0SOfflineSgnNoQ = (TH2F*)cosPAK0SOfflineSgnNoQ->Clone();
  TH2F* allcosThePrOfflineSgnC = (TH2F*)cosThePrOfflineSgnC->Clone();
  TH2F* allcosThePrOfflineSgnB = (TH2F*)cosThePrOfflineSgnB->Clone();
  TH2F* allcosThePrOfflineSgnNoQ = (TH2F*)cosThePrOfflineSgnNoQ->Clone();
  TH2F* allresignedD0OfflineSgnC = (TH2F*)resignedD0OfflineSgnC->Clone();
  TH2F* allresignedD0OfflineSgnB = (TH2F*)resignedD0OfflineSgnB->Clone();
  TH2F* allresignedD0OfflineSgnNoQ = (TH2F*)resignedD0OfflineSgnNoQ->Clone();

  TH2F* pidptK0SOfflineSgnC = (TH2F*)ptK0SOfflineSgnC->Clone();
  TH2F* pidptK0SOfflineSgnB = (TH2F*)ptK0SOfflineSgnB->Clone();
  TH2F* pidptK0SOfflineSgnNoQ = (TH2F*)ptK0SOfflineSgnNoQ->Clone();
  TH2F* pidptPOfflineSgnC = (TH2F*)ptPOfflineSgnC->Clone();
  TH2F* pidptPOfflineSgnB = (TH2F*)ptPOfflineSgnB->Clone();
  TH2F* pidptPOfflineSgnNoQ = (TH2F*)ptPOfflineSgnNoQ->Clone();
  TH2F* pidptPiPOfflineSgnC = (TH2F*)ptPiPOfflineSgnC->Clone();
  TH2F* pidptPiPOfflineSgnB = (TH2F*)ptPiPOfflineSgnB->Clone();
  TH2F* pidptPiPOfflineSgnNoQ = (TH2F*)ptPiPOfflineSgnNoQ->Clone();
  TH2F* pidptPiMOfflineSgnC = (TH2F*)ptPiMOfflineSgnC->Clone();
  TH2F* pidptPiMOfflineSgnB = (TH2F*)ptPiMOfflineSgnB->Clone();
  TH2F* pidptPiMOfflineSgnNoQ = (TH2F*)ptPiMOfflineSgnNoQ->Clone();
  TH2F* pidd0K0SOfflineSgnC = (TH2F*)d0K0SOfflineSgnC->Clone();
  TH2F* pidd0K0SOfflineSgnB = (TH2F*)d0K0SOfflineSgnB->Clone();
  TH2F* pidd0K0SOfflineSgnNoQ = (TH2F*)d0K0SOfflineSgnNoQ->Clone();
  TH2F* pidd0POfflineSgnC = (TH2F*)d0POfflineSgnC->Clone();
  TH2F* pidd0POfflineSgnB = (TH2F*)d0POfflineSgnB->Clone();
  TH2F* pidd0POfflineSgnNoQ = (TH2F*)d0POfflineSgnNoQ->Clone();
  TH2F* pidcosPAK0SOfflineSgnC = (TH2F*)cosPAK0SOfflineSgnC->Clone();
  TH2F* pidcosPAK0SOfflineSgnB = (TH2F*)cosPAK0SOfflineSgnB->Clone();
  TH2F* pidcosPAK0SOfflineSgnNoQ = (TH2F*)cosPAK0SOfflineSgnNoQ->Clone();
  TH2F* pidcosThePrOfflineSgnC = (TH2F*)cosThePrOfflineSgnC->Clone();
  TH2F* pidcosThePrOfflineSgnB = (TH2F*)cosThePrOfflineSgnB->Clone();
  TH2F* pidcosThePrOfflineSgnNoQ = (TH2F*)cosThePrOfflineSgnNoQ->Clone();
  TH2F* pidresignedD0OfflineSgnC = (TH2F*)resignedD0OfflineSgnC->Clone();
  TH2F* pidresignedD0OfflineSgnB = (TH2F*)resignedD0OfflineSgnB->Clone();
  TH2F* pidresignedD0OfflineSgnNoQ = (TH2F*)resignedD0OfflineSgnNoQ->Clone();

  fOutputAll->Add(allptK0SOfflineSgnC);
  fOutputAll->Add(allptK0SOfflineSgnB);
  fOutputAll->Add(allptK0SOfflineSgnNoQ);
  fOutputAll->Add(allptPOfflineSgnC);
  fOutputAll->Add(allptPOfflineSgnB);
  fOutputAll->Add(allptPOfflineSgnNoQ);
  fOutputAll->Add(allptPiPOfflineSgnC);
  fOutputAll->Add(allptPiPOfflineSgnB);
  fOutputAll->Add(allptPiPOfflineSgnNoQ);
  fOutputAll->Add(allptPiMOfflineSgnC);
  fOutputAll->Add(allptPiMOfflineSgnB);
  fOutputAll->Add(allptPiMOfflineSgnNoQ);
  fOutputAll->Add(alld0K0SOfflineSgnC);
  fOutputAll->Add(alld0K0SOfflineSgnB);
  fOutputAll->Add(alld0K0SOfflineSgnNoQ);
  fOutputAll->Add(alld0POfflineSgnC);
  fOutputAll->Add(alld0POfflineSgnB);
  fOutputAll->Add(alld0POfflineSgnNoQ);
  fOutputAll->Add(allcosPAK0SOfflineSgnC);
  fOutputAll->Add(allcosPAK0SOfflineSgnB);
  fOutputAll->Add(allcosPAK0SOfflineSgnNoQ);
  fOutputAll->Add(allcosThePrOfflineSgnC);
  fOutputAll->Add(allcosThePrOfflineSgnB);
  fOutputAll->Add(allcosThePrOfflineSgnNoQ);
  fOutputAll->Add(allresignedD0OfflineSgnC);
  fOutputAll->Add(allresignedD0OfflineSgnB);
  fOutputAll->Add(allresignedD0OfflineSgnNoQ);

  fOutputPIDBach->Add(pidptK0SOfflineSgnC);
  fOutputPIDBach->Add(pidptK0SOfflineSgnB);
  fOutputPIDBach->Add(pidptK0SOfflineSgnNoQ);
  fOutputPIDBach->Add(pidptPOfflineSgnC);
  fOutputPIDBach->Add(pidptPOfflineSgnB);
  fOutputPIDBach->Add(pidptPOfflineSgnNoQ);
  fOutputPIDBach->Add(pidptPiPOfflineSgnC);
  fOutputPIDBach->Add(pidptPiPOfflineSgnB);
  fOutputPIDBach->Add(pidptPiPOfflineSgnNoQ);
  fOutputPIDBach->Add(pidptPiMOfflineSgnC);
  fOutputPIDBach->Add(pidptPiMOfflineSgnB);
  fOutputPIDBach->Add(pidptPiMOfflineSgnNoQ);
  fOutputPIDBach->Add(pidd0K0SOfflineSgnC);
  fOutputPIDBach->Add(pidd0K0SOfflineSgnB);
  fOutputPIDBach->Add(pidd0K0SOfflineSgnNoQ);
  fOutputPIDBach->Add(pidd0POfflineSgnC);
  fOutputPIDBach->Add(pidd0POfflineSgnB);
  fOutputPIDBach->Add(pidd0POfflineSgnNoQ);
  fOutputPIDBach->Add(pidcosPAK0SOfflineSgnC);
  fOutputPIDBach->Add(pidcosPAK0SOfflineSgnB);
  fOutputPIDBach->Add(pidcosPAK0SOfflineSgnNoQ);
  fOutputPIDBach->Add(pidcosThePrOfflineSgnC);
  fOutputPIDBach->Add(pidcosThePrOfflineSgnB);
  fOutputPIDBach->Add(pidcosThePrOfflineSgnNoQ);
  fOutputPIDBach->Add(pidresignedD0OfflineSgnC);
  fOutputPIDBach->Add(pidresignedD0OfflineSgnB);
  fOutputPIDBach->Add(pidresignedD0OfflineSgnNoQ);

}
