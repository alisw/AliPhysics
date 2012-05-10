/**************************************************************************
* Copyright(c) 1998-2008,ALICE Experiment at CERN,All rights reserved.    *
*                                                                         *
* Author: The ALICE Off-line Project.                                     *
* Contributors are mentioned in the code where appropriate.               *
*                                                                         *
* Permission to use,copy,modify and distribute this software and its      *
* documentation strictly for non-commercial purposes is hereby granted    *
* without fee,provided that the above copyright notice appears in all     *
* copies and that both the copyright notice and this permission notice    *
* appear in the supporting documentation. The authors make no claims      *
* about the suitability of this software for any purpose. It is           *
* provided "as is" without express or implied warranty.                   *
**************************************************************************/

//==============================================================================
// FlowD2H main task:
// >> Select candidates and passes flowevents to the daughter tasks.
// >> The POIcuts are polymorphic based on the AliRDHFCuts class allowing the 
//    use of all charmed candidates reconstructed in the central barrel.
// Author: Carlos Perez (cperez@cern.ch)
//==============================================================================

/* $Id$ */

#include "TChain.h"

#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"

#include "AliCentrality.h"

#include "AliAODEvent.h"
#include "AliAODInputHandler.h"
#include "AliAODTrack.h"

#include "AliRDHFCuts.h"

#include "AliRDHFCutsD0toKpi.h"
#include "AliRDHFCutsDStartoKpipi.h"
#include "AliRDHFCutsDplustoKpipi.h"

#include "AliAODRecoDecayHF2Prong.h"
#include "AliAODRecoDecayHF3Prong.h"
#include "AliAODRecoCascadeHF.h"

#include "TMath.h"

#include "AliFlowCommonConstants.h"
#include "AliFlowEvent.h"
#include "AliFlowTrackCuts.h"
#include "AliFlowCandidateTrack.h"

#include "TObjArray.h"
#include "TList.h"
#include "TH1D.h"

#include "AliAnalysisTaskFlowD2H.h"

ClassImp(AliAnalysisTaskFlowD2H)

//_____________________________________________________________________________
AliAnalysisTaskFlowD2H::AliAnalysisTaskFlowD2H() :
  AliAnalysisTaskSE(), fTPCEvent(NULL), fVZEEvent(NULL), 
  fCutsTPC(NULL), fCutsVZE(NULL), fNoPOIs(NULL), fCutsPOI(NULL),
  fSource(0), fDebugV2(kFALSE), fMassBins(0), fMinMass(0.),
  fMaxMass(0.), fPtBinWidth(0), fHList(NULL), fEvent(NULL), 
  fCC(NULL), fRFPMTPC(NULL), fRFPPhiTPC(NULL), fCandidates(NULL)
{
// Default constructor
}

//_____________________________________________________________________________
AliAnalysisTaskFlowD2H::AliAnalysisTaskFlowD2H(const char *name,	
					       AliFlowTrackCuts *cutsTPC,
					       AliFlowTrackCuts *cutsVZE,
					       AliRDHFCuts *cutsPOIs,
					       Int_t specie) :
  AliAnalysisTaskSE(name), fTPCEvent(NULL), fVZEEvent(NULL), 
  fCutsTPC(cutsTPC), fCutsVZE(cutsVZE), fNoPOIs(NULL), fCutsPOI(cutsPOIs),
  fSource(specie), fDebugV2(kFALSE), fMassBins(0), fMinMass(0.),
  fMaxMass(0.), fPtBinWidth(0), fHList(NULL), fEvent(NULL),
  fCC(NULL), fRFPMTPC(NULL), fRFPPhiTPC(NULL), fCandidates(NULL)
{
// Standard constructor
  DefineInput( 0,TChain::Class());
  DefineOutput(1,TList::Class());
  DefineOutput(2,AliFlowEventSimple::Class());
  DefineOutput(3,AliFlowEventSimple::Class());
}

//_____________________________________________________________________________
AliAnalysisTaskFlowD2H::~AliAnalysisTaskFlowD2H(){
  // delete objects
  if(fTPCEvent) delete fTPCEvent;
  if(fVZEEvent) delete fVZEEvent;
  if(fCutsTPC) delete fCutsTPC;
  if(fCutsVZE) delete fCutsVZE;
  if(fNoPOIs) delete fNoPOIs;
  if(fCutsPOI) delete fCutsPOI;
  if(fHList) delete fHList;
  if(fCandidates) delete fCandidates;
}

//_____________________________________________________________________________
void AliAnalysisTaskFlowD2H::UserCreateOutputObjects(){
  // Define output objects + parameters
  if(fDebugV2)
    printf("DEBUGGER: Creating output\n");
  fHList = new TList();
  fHList->SetOwner();
  fEvent = new TH1D("Events","Events",3,0,3);
  fEvent->GetXaxis()->SetBinLabel(1,"REACHED");
  fEvent->GetXaxis()->SetBinLabel(2,"SELECTED");
  fEvent->GetXaxis()->SetBinLabel(3,"DELTA AOD REACHED");
  fHList->Add( fEvent );
  fCC = new TH1D("CentralityClass","Centrality Class",50,0,100);
  fHList->Add( fCC );
  fRFPMTPC = new TH1D("RFPMultiplicityTPC","RFP Multiplicity TPC",300,0,3000);
  fHList->Add( fRFPMTPC );
  fRFPPhiTPC = new TH1D("RFPPhiTPC","RFP Phi TPC",180,0,TMath::TwoPi());
  fHList->Add( fRFPPhiTPC );

  fCandidates = new TObjArray(100);
  fCandidates->SetOwner();

  AliFlowCommonConstants* cc = AliFlowCommonConstants::GetMaster();
  cc->SetNbinsMult(1);
  cc->SetNbinsPt(24/fPtBinWidth);
  cc->SetNbinsPhi(1);
  cc->SetNbinsEta(15);
  cc->SetNbinsQ(1);
  cc->SetNbinsMass( fMassBins );
  cc->SetMultMin(1);
  cc->SetMultMax(2);
  cc->SetPtMin(0);
  cc->SetPtMax(24);
  cc->SetPhiMin(0);
  cc->SetPhiMax(TMath::TwoPi());
  cc->SetEtaMin(-3.9);
  cc->SetEtaMax(+5.1);
  cc->SetQMin(0);
  cc->SetQMax(1);
  cc->SetMassMin( fMinMass );
  cc->SetMassMax( fMaxMass );

  fTPCEvent = new AliFlowEvent(3000);
  fVZEEvent = new AliFlowEvent(170);

  fNoPOIs = new AliFlowTrackCuts( "noPOIs" );
  fNoPOIs->SetParamType(AliFlowTrackCuts::kGlobal);
  fNoPOIs->SetPtRange(+1,-1);

  PostData(1,fHList);
  PostData(2,fTPCEvent);
  PostData(3,fVZEEvent);
}
//_____________________________________________________________________________
void AliAnalysisTaskFlowD2H::UserExec(Option_t *)
{
  // Do analysis + fIll histograms
  AliAODEvent *fAOD = dynamic_cast<AliAODEvent*>(InputEvent());

  if(!fAOD) return;
  fEvent->Fill( 0 );

  // floweventcuts::isselected() and alirdhfcuts::iseventselected() cut on the same 
  // values in the same way BUT the latter also loads the PIDresponse object from the 
  // event header!!! 
  //  if(!fEventCuts->IsSelected(fAOD)) return;
  if(!fCutsPOI->IsEventSelected(fAOD)) return;
  fEvent->Fill( 1 );

  fCC->Fill( fCutsPOI->GetCentrality(fAOD) );

  fCutsTPC->SetEvent( fAOD, MCEvent() );
  fCutsVZE->SetEvent( fAOD, MCEvent() );
  fNoPOIs->SetEvent( fAOD, MCEvent() );
  fTPCEvent->Fill( fCutsTPC, fNoPOIs );
  fVZEEvent->Fill( fCutsVZE, fNoPOIs );

  Int_t rfps = fTPCEvent->GetNumberOfRPs();
  fRFPMTPC->Fill( rfps );
  for(int iRPs=0; iRPs!=rfps; ++iRPs ) {
    AliFlowTrack *iRP = dynamic_cast<AliFlowTrack*>(fTPCEvent->GetTrack( iRPs ));
    if (!iRP) continue;
    fRFPPhiTPC->Fill( iRP->Phi() );
  }

  if (fDebugV2) printf("Event selected\n");
  fCandidates->SetLast(-1); // resets the array

  switch (fSource) {
    case (AliRDHFCuts::kD0toKpiCuts):
      FillD0toKpi(fAOD); break;
    case (AliRDHFCuts::kDstarCuts):
      FillDStartoKpipi(fAOD); break;
    case (AliRDHFCuts::kDplusCuts):
      FillDplustoKpipi(fAOD); break;
  }

  if(fDebugV2) printf("TPCevent %d | VZEevent %d\n", fTPCEvent->NumberOfTracks(), fVZEEvent->NumberOfTracks() );
  //inject candidates
  if (fDebugV2)  printf("I received %d candidates\n",fCandidates->GetEntriesFast());
  for(int iCand=0; iCand!=fCandidates->GetEntriesFast(); ++iCand ) {
    AliFlowCandidateTrack *cand = dynamic_cast<AliFlowCandidateTrack*>(fCandidates->At(iCand));
    if (!cand) continue;
    if (fDebugV2) printf(" >Checking at candidate %d with %d daughters: mass %f\n",iCand,cand->GetNDaughters(),cand->Mass());
    for(int iDau=0; iDau!=cand->GetNDaughters(); ++iDau) {
      if(fDebugV2) printf("  >Daughter %d with fID %d", iDau, cand->GetIDDaughter(iDau));
      for(int iRPs=0; iRPs!=fTPCEvent->NumberOfTracks(); ++iRPs ) {
	AliFlowTrack *iRP = dynamic_cast<AliFlowTrack*>(fTPCEvent->GetTrack( iRPs ));
	if (!iRP) continue;
	if( !iRP->InRPSelection() ) continue;
	if( cand->GetIDDaughter(iDau) == iRP->GetID() ) {
	  if(fDebugV2) printf(" was in RP set");
	  iRP->SetForRPSelection(kFALSE);
	  fTPCEvent->SetNumberOfRPs( fTPCEvent->GetNumberOfRPs() -1 );
	}
      }
      if(fDebugV2) printf("\n");
    }
    cand->SetForPOISelection(kTRUE);
    fTPCEvent->InsertTrack( ((AliFlowTrack*) cand) );
    fVZEEvent->InsertTrack( ((AliFlowTrack*) cand) );
  }
  if(fDebugV2) printf("TPCevent %d | VZEevent %d\n", fTPCEvent->NumberOfTracks(), fVZEEvent->NumberOfTracks() );

  PostData(1,fHList);
  PostData(2,fTPCEvent);
  PostData(3,fVZEEvent);
}
//______________________________________________________________________________
void AliAnalysisTaskFlowD2H::FillD0toKpi(const AliAODEvent *theAOD)
{
  // Fill D0->Kpi histos
  TList *listHF = (TList*) theAOD->GetList();
  if(!listHF) return;
  TClonesArray *listDzero = (TClonesArray*) listHF->FindObject("D0toKpi");
  if(!listDzero) return;
  int nEntries = listDzero->GetEntriesFast();
  if( !nEntries ) return;
  fEvent->Fill( 2 ); // EVERYTHING OKAY

  AliRDHFCutsD0toKpi *fCutsD0toKpi = (AliRDHFCutsD0toKpi*) fCutsPOI;
  if (fDebugV2) printf("  ᶫ%d candidates found. Looping...\n", nEntries);
  Int_t nIDs[2];
  for( int iEntry=0; iEntry!=nEntries; ++iEntry ) {
    AliAODRecoDecayHF2Prong *d0cand = (AliAODRecoDecayHF2Prong*) listDzero->UncheckedAt( iEntry );
    if( !d0cand ) continue;
    // APPLYING CUTS
    if( !d0cand->HasSelectionBit(AliRDHFCuts::kD0toKpiCuts) ) continue;
    if( !fCutsD0toKpi->IsInFiducialAcceptance(d0cand->Pt(),d0cand->Y(421)) )continue;
    Int_t topCut = fCutsD0toKpi->IsSelected( d0cand, AliRDHFCuts::kAll, NULL );
    if(!topCut) continue;
    // TO HANDLE AUTOCORRELATIONS
    nIDs[0] = ( (AliAODTrack*)d0cand->GetDaughter(0) )->GetID();
    nIDs[1] = ( (AliAODTrack*)d0cand->GetDaughter(1) )->GetID();
    if(topCut&1){
      Double_t massD0=d0cand->InvMassD0();
      MakeTrack(massD0, d0cand->Pt(), d0cand->Phi(), d0cand->Eta(), 2, nIDs);
      if(fDebugV2) printf("   ᶫInjecting D0 candidate\n");
    }
    if(topCut&2){
      Double_t massD0=d0cand->InvMassD0bar();
      MakeTrack(massD0, d0cand->Pt(), d0cand->Phi(), d0cand->Eta(), 2, nIDs);
    }
    if(fDebugV2) printf("   ᶫInjecting D0 candidate\n");
  }
}
//______________________________________________________________________________
void AliAnalysisTaskFlowD2H::FillDStartoKpipi(const AliAODEvent *theAOD )
{
  // Fills D* to Kpipi
  TList *listHF = (TList*) theAOD->GetList();
  if(!listHF) return;
  TClonesArray *listDstar = (TClonesArray*) listHF->FindObject("Dstar");
  if(!listDstar) return;
  int nEntries = listDstar->GetEntriesFast();
  if( !nEntries ) return;
  fEvent->Fill( 2 ); // EVERYTHING OKAY

  AliRDHFCutsDStartoKpipi *fCutsDStartoKpipi = (AliRDHFCutsDStartoKpipi*) fCutsPOI;
  if (fDebugV2) printf("  ᶫ%d candidates found. Looping...\n", nEntries);
  Int_t nIDs[3];
  for( int iEntry=0; iEntry!=nEntries; ++iEntry ) {
    AliAODRecoCascadeHF *dst = (AliAODRecoCascadeHF*) listDstar->UncheckedAt( iEntry );
    if( !dst ) continue;
    AliAODRecoDecayHF2Prong *d0cand = (AliAODRecoDecayHF2Prong*)dst->Get2Prong();
    if(!d0cand) continue;
    if( !d0cand->GetDaughter(0) ) continue;
    if( !d0cand->GetDaughter(1) ) continue;
    if( !dst->GetBachelor() ) continue;
    // APPLYING CUTS
    if( !dst->HasSelectionBit(AliRDHFCuts::kDstarCuts) ) continue;
    if( !fCutsDStartoKpipi->IsInFiducialAcceptance(dst->Pt(),dst->YDstar()) )continue;
    Int_t topCut=0;
    if(dst->Pt() > 5) {
      fCutsDStartoKpipi->SetUsePID(kFALSE);
      topCut = fCutsDStartoKpipi->IsSelected( dst, AliRDHFCuts::kAll );
      fCutsDStartoKpipi->SetUsePID(kTRUE);
    } else topCut = fCutsDStartoKpipi->IsSelected( dst, AliRDHFCuts::kAll ); 
    if(!topCut) continue;
    Double_t massDS=dst->DeltaInvMass();
    // TO HANDLE AUTOCORRELATIONS
    nIDs[0] = ((AliAODTrack*)d0cand->GetDaughter(0))->GetID();
    nIDs[1] = ((AliAODTrack*)d0cand->GetDaughter(1))->GetID();
    nIDs[2] = ((AliAODTrack*)dst->GetBachelor() )->GetID();
    // ADDING TRACK
    MakeTrack(massDS, dst->Pt(), dst->Phi(), dst->Eta(), 3, nIDs);
    if(fDebugV2) printf("   ᶫInjecting DStar candidate %d\n",iEntry);
  }
}
//______________________________________________________________________________
void AliAnalysisTaskFlowD2H::FillDplustoKpipi(const AliAODEvent *theAOD )
{
  // Fill D+ to Kpipi
  TList *listHF = (TList*) theAOD->GetList();
  if(!listHF) return;
  TClonesArray *listDplus = (TClonesArray*) listHF->FindObject("Charm3Prong");
  if(!listDplus) return;
  int nEntries = listDplus->GetEntriesFast();
  if( !nEntries ) return;
  fEvent->Fill( 2 ); // EVERYTHING OKAY

  AliRDHFCutsDplustoKpipi *fCutsDStartoKpipi = (AliRDHFCutsDplustoKpipi*) fCutsPOI;
  if (fDebugV2) printf("  ᶫ%d candidates found. Looping...\n", nEntries);
  Int_t nIDs[3];
  for( int iEntry=0; iEntry!=nEntries; ++iEntry ) {
    AliAODRecoDecayHF3Prong *dplu = (AliAODRecoDecayHF3Prong*) listDplus->UncheckedAt( iEntry );
    if( !dplu ) continue;
    // APPLYING CUTS
    if( !dplu->HasSelectionBit(AliRDHFCuts::kDplusCuts) ) continue;
    if( !fCutsDStartoKpipi->IsInFiducialAcceptance(dplu->Pt(),dplu->YDplus()) )continue;
    Int_t topCut = fCutsDStartoKpipi->IsSelected( dplu, AliRDHFCuts::kAll );
    if(!topCut) continue;
    Double_t massDp=dplu->InvMassDplus();
    // TO HANDLE AUTOCORRELATIONS
    nIDs[0] = ((AliAODTrack*)dplu->GetDaughter(0))->GetID();
    nIDs[1] = ((AliAODTrack*)dplu->GetDaughter(1))->GetID();
    nIDs[2] = ((AliAODTrack*)dplu->GetDaughter(2))->GetID();
    // ADDING TRACK
    MakeTrack(massDp, dplu->Pt(), dplu->Phi(), dplu->Eta(), 3, nIDs);
    if(fDebugV2) printf("   ᶫInjecting Dplus candidate %d\n",iEntry);
  }
}
//_______________________________________________________________________________
void AliAnalysisTaskFlowD2H::MakeTrack( Double_t mass, Double_t pt, Double_t phi, 
					Double_t eta, Int_t nDau, const Int_t *iID ) {
  // create track for flow tasks
  Bool_t overwrite = kTRUE;
  AliFlowCandidateTrack *oTrack = (static_cast<AliFlowCandidateTrack*> (fCandidates->At( fCandidates->GetLast()+1 )));
  if( !oTrack ) { // creates new
    oTrack = new AliFlowCandidateTrack();
    overwrite = kFALSE;
  } else { // overwrites
    oTrack->ClearMe();
  }
  oTrack->SetMass(mass);
  oTrack->SetPt(pt);
  oTrack->SetPhi(phi);
  oTrack->SetEta(eta);
  for(Int_t iDau=0; iDau!=nDau; ++iDau)
    oTrack->AddDaughter(iID[iDau]);
  oTrack->SetForPOISelection(kTRUE);
  oTrack->SetForRPSelection(kFALSE);
  if(overwrite) {
    fCandidates->SetLast( fCandidates->GetLast()+1 );
  } else {
    fCandidates->AddLast(oTrack);
  }
  return;
}

void AliAnalysisTaskFlowD2H::SetCommonConstants(Int_t massBins, Double_t minMass, Double_t maxMass, Int_t ptWidth) {
  fMassBins = massBins;
  fMinMass = minMass;
  fMaxMass = maxMass;
  fPtBinWidth = ptWidth;
}
