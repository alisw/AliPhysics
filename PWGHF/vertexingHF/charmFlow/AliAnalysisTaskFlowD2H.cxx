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
// >> Select candidates and passes the array to the daughter tasks.
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

#include "AliFlowEventCuts.h"
#include "AliFlowCandidateTrack.h"

#include "TObjArray.h"
#include "TList.h"
#include "TH1D.h"

#include "AliAnalysisTaskFlowD2H.h"

ClassImp(AliAnalysisTaskFlowD2H)

//_____________________________________________________________________________
AliAnalysisTaskFlowD2H::AliAnalysisTaskFlowD2H() :
AliAnalysisTaskSE(), fEventCuts(NULL), fCutsPOI(NULL), fSource(0),
  fDebugV2(kFALSE), fHList(NULL), fEvent(NULL), fCandidates(NULL)
{
// Default constructor
}

//_____________________________________________________________________________
AliAnalysisTaskFlowD2H::AliAnalysisTaskFlowD2H(const char *name,
					       AliFlowEventCuts *eventCuts,
					       AliRDHFCuts *cutsPOIs,
					       Int_t specie) :
  AliAnalysisTaskSE(name), fEventCuts(eventCuts), fCutsPOI(cutsPOIs),
  fSource(specie), fDebugV2(kFALSE), fHList(NULL), fEvent(NULL), fCandidates(NULL)
{
// Standard constructor
  DefineInput( 0,TChain::Class());
  DefineOutput(1,TList::Class());
  DefineOutput(2,TObjArray::Class());
}

//_____________________________________________________________________________
AliAnalysisTaskFlowD2H::~AliAnalysisTaskFlowD2H(){
  // delete objects
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

  fCandidates = new TObjArray(300);
  fCandidates->SetOwner();

  PostData(1,fHList);
  PostData(2,fCandidates);
}
//_____________________________________________________________________________
void AliAnalysisTaskFlowD2H::UserExec(Option_t *)
{
  // Do analysis + fIll histograms
  AliAODEvent *fAOD = dynamic_cast<AliAODEvent*>(InputEvent());

  if(!fAOD) return;
  fEvent->Fill( 0 );

  if(!fEventCuts->IsSelected(fAOD)) return;
  fEvent->Fill( 1 );

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

  if (fDebugV2) printf("Candidates inserted: %d\n", fCandidates->GetEntriesFast() );
  PostData(1,fHList);
  PostData(2,fCandidates);

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
    Double_t massD0=topCut>1?d0cand->InvMassD0bar():d0cand->InvMassD0();
    // TO HANDLE AUTOCORRELATIONS
    nIDs[0] = ( (AliAODTrack*)d0cand->GetDaughter(0) )->GetID();
    nIDs[1] = ( (AliAODTrack*)d0cand->GetDaughter(1) )->GetID();
    // ADDING TRACK
    MakeTrack(massD0, d0cand->Pt(), d0cand->Phi(), d0cand->Eta(), 2, nIDs);
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
