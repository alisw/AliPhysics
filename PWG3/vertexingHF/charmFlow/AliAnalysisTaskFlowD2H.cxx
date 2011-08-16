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
// >> Make flowEvent with RPcuts and POIcuts given in constructor and passes it
//    to the daughter tasks.
// >> The POIcuts are polymorphic based on the AliRDHFCuts class allowing the 
//    use of all charmed candidates reconstructed in the central barrel.
// Author: Carlos Perez (cperez@cern.ch)
//==============================================================================

/* $Id$ */

#include "TChain.h"
#include "TList.h"
#include "TH2D.h"
#include "TProfile.h"

#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"

#include "AliCentrality.h"

#include "AliAODEvent.h"
#include "AliAODInputHandler.h"
#include "AliAODTrack.h"

#include "TMath.h"
#include "TObjArray.h"

#include "AliFlowCandidateTrack.h"
#include "AliFlowTrackCuts.h"
#include "AliFlowEvent.h"
#include "AliFlowCommonConstants.h"

#include "AliRDHFCuts.h"
#include "AliRDHFCutsD0toKpi.h"
#include "AliRDHFCutsDStartoKpipi.h"
#include "AliRDHFCutsDplustoKpipi.h"

#include "AliAODRecoDecayHF2Prong.h"
#include "AliAODRecoDecayHF3Prong.h"
#include "AliAODRecoCascadeHF.h"

#include "AliAnalysisTaskFlowD2H.h"

ClassImp(AliAnalysisTaskFlowD2H)

//_____________________________________________________________________________
AliAnalysisTaskFlowD2H::AliAnalysisTaskFlowD2H() :
  AliAnalysisTaskSE(), fCutsRP(NULL), fCutsPOI(NULL), fSource(0),
  fDebug(kFALSE), fHList(NULL), fAnaCuts(NULL)
{
// Default constructor
  for(int i=0; i!=2; ++i)
    fEvent[i] = NULL;
  for(int i=0; i!=2; ++i)
    fMass[i] = NULL;
  for(int i=0; i!=2; ++i)
    fPOIEta[i] = 0;
  for(int i=0; i!=4; ++i)
    fFlowEta[i] = 0;
  for(int x=0; x!=2; ++x)
    fFlowPts[x] = 0;
  for(int x=0; x!=2; ++x)
    for(int m=0; m!=5; ++m)
      fFlowBands[x][m] = 0;
}

//_____________________________________________________________________________
AliAnalysisTaskFlowD2H::AliAnalysisTaskFlowD2H(const char *name,
    AliFlowTrackCuts *cutsRPs, AliRDHFCuts *cutsPOIs, Int_t specie) :
  AliAnalysisTaskSE(name), fCutsRP(cutsRPs), 
  fCutsPOI(cutsPOIs), fSource(specie),
  fDebug(kFALSE), fHList(NULL), fAnaCuts(NULL)
{
// Standard constructor
  for(int i=0; i!=2; ++i)
    fEvent[i] = NULL;
  for(int i=0; i!=2; ++i)
    fMass[i] = NULL;
  for(int i=0; i!=2; ++i)
    fPOIEta[i] = 0;
  for(int i=0; i!=4; ++i)
    fFlowEta[i] = 0;
  for(int x=0; x!=2; ++x)
    fFlowPts[x] = 0;
  for(int x=0; x!=2; ++x)
    for(int m=0; m!=5; ++m)
      fFlowBands[x][m] = 0;

  DefineInput( 0,TChain::Class());
  DefineOutput(1,TList::Class());
  DefineOutput(2,AliFlowEventSimple::Class()); // first band
  DefineOutput(3,AliFlowEventSimple::Class()); // second band
  DefineOutput(4,AliFlowEventSimple::Class()); // third band
  DefineOutput(5,AliFlowEventSimple::Class()); // fourth band
  DefineOutput(6,AliFlowEventSimple::Class()); // fifth band
}

//_____________________________________________________________________________
AliAnalysisTaskFlowD2H::~AliAnalysisTaskFlowD2H()
{
  if(fHList)
    delete fHList;
}

//_____________________________________________________________________________
void AliAnalysisTaskFlowD2H::UserCreateOutputObjects()
{
  if(fDebug)
    printf("DEBUGGER: Creating output\n");
  fHList = new TList();
  fHList->SetOwner();
  AddHistograms();
  PostData(1,fHList);

  AliFlowCommonConstants* cc = AliFlowCommonConstants::GetMaster();
  cc->SetNbinsMult(10000);
  cc->SetMultMin(0);
  cc->SetMultMax(10000);

  cc->SetNbinsPt(fFlowPts[1]-fFlowPts[0]);
  cc->SetPtMin(fFlowPts[0]);
  cc->SetPtMax(fFlowPts[1]);

  cc->SetNbinsPhi(180);
  cc->SetPhiMin(0.0);
  cc->SetPhiMax(TMath::TwoPi());

  cc->SetNbinsEta(200);
  cc->SetEtaMin(-5.0);
  cc->SetEtaMax(+5.0);

  cc->SetNbinsQ(500);
  cc->SetQMin(0.0);
  cc->SetQMax(3.0);
}
//_____________________________________________________________________________
void AliAnalysisTaskFlowD2H::AddHistograms()
{
  TList *tEvents = new TList();
    tEvents->SetName("Events");
    tEvents->SetOwner();
    for(int i=0; i!=2; ++i) {
      fEvent[i] = new TH2D(Form("Event%d",i),"Events;V0M;Arb",10,0,100,10,0,12000);
      tEvents->Add(fEvent[i]);
    }
    fAnaCuts=new TProfile("Cuts","Analysis Cuts", 10,0,10);
      fAnaCuts->Fill(0.5,fPOIEta[0],1);fAnaCuts->GetXaxis()->SetBinLabel(1,"ETAm");
      fAnaCuts->Fill(1.5,fPOIEta[1],1);fAnaCuts->GetXaxis()->SetBinLabel(2,"ETAM");
    tEvents->Add(fAnaCuts);
    fHList->Add(tEvents);
  TList *tCandidates;
    tCandidates = new TList();
    tCandidates->SetOwner();
    tCandidates->SetName(Form("Candidates%d",fSource));
    Double_t dBinMin=0, dBinMax=3;
    Int_t nBins=3;
    switch (fSource) {
      case (AliRDHFCuts::kD0toKpiCuts): // 360/72 (bw=5MeV)
        dBinMin = 1.695; dBinMax = 2.055; nBins=72; break;
      case (AliRDHFCuts::kDstarCuts): // 36/72 (bw=0.5MeV)
        dBinMin = 0.137; dBinMax = 0.173; nBins=72; break;
      case (AliRDHFCuts::kDplusCuts): // 360/72 (bw=5MeV)
        dBinMin = 1.695; dBinMax = 2.055; nBins=72; break;
    }
    for(int i=0; i!=2; ++i) {
      fMass[i] = new TH2D( Form("Mass%d",i),
                           Form("Mass%d;Mass [GeV];Pt [GeV]",i),
                           nBins, dBinMin, dBinMax,
                           fFlowPts[1]-fFlowPts[0], fFlowPts[0], fFlowPts[1] );
      tCandidates->Add(fMass[i]);
    }
    fHList->Add(tCandidates);
  if (fDebug) printf("DEBUGGER: Created histos for DMeson %d\n", fSource );
}
//_____________________________________________________________________________
void AliAnalysisTaskFlowD2H::NotifyRun()
{
}
//_____________________________________________________________________________
void AliAnalysisTaskFlowD2H::UserExec(Option_t *)
{
  AliAODEvent *fAOD = dynamic_cast<AliAODEvent*>(InputEvent());

  if(!fAOD)
    return;
  Int_t iMulti = fAOD->GetNTracks();  /// TO DO
  fEvent[0]->Fill( fCutsPOI->GetCentrality(fAOD), iMulti );
  if(!fCutsPOI->IsEventSelected(fAOD)) return;
  if(fCutsPOI->IsEventSelectedInCentrality(fAOD)>0) return;
  fEvent[1]->Fill( fCutsPOI->GetCentrality(fAOD), iMulti );

  if (fDebug) printf("Event selected\n");
  fCutsRP->SetEvent(fAOD,MCEvent());
  AliFlowTrackCuts* dummy = new AliFlowTrackCuts("null_cuts");
  dummy->SetParamType(AliFlowTrackCuts::kGlobal);
  dummy->SetPtRange(+1,-1); // select nothing QUICK
  dummy->SetEtaRange(+1,-1); // select nothing VZERO
  dummy->SetEvent(fAOD,MCEvent());
  AliFlowEvent *flowEvent[5];
  for(int r=0; r!=5; ++r) {
    flowEvent[r] = new AliFlowEvent(fCutsRP,dummy);
    flowEvent[r]->SetReferenceMultiplicity( iMulti );
    flowEvent[r]->DefineDeadZone(0,0,0,0);
    flowEvent[r]->TagSubeventsInEta( fFlowEta[0], fFlowEta[1],
                                     fFlowEta[2], fFlowEta[3] );
  }
  delete dummy;
  if (fDebug) printf(" ᶫFlowEvent has %d RPs\n", flowEvent[0]->NumberOfTracks() );

  switch (fSource) {
    case (AliRDHFCuts::kD0toKpiCuts):
      FillD0toKpi(fAOD,flowEvent); break;
    case (AliRDHFCuts::kDstarCuts):
      FillDStartoKpipi(fAOD,flowEvent); break;
    case (AliRDHFCuts::kDplusCuts):
      FillDplustoKpipi(fAOD,flowEvent); break;
  }

  for(int m=0; m!=5; ++m)
    PostData(2+m,flowEvent[m]);
  PostData(1,fHList);

}
//______________________________________________________________________________
void AliAnalysisTaskFlowD2H::FillD0toKpi(AliAODEvent *theAOD, 
                                         AliFlowEvent *theMB[5] )
{
  TList *listHF = (TList*) theAOD->GetList();
  if(!listHF) return;
  TClonesArray *listDzero = (TClonesArray*) listHF->FindObject("D0toKpi");
  if(!listDzero) return;
  int nEntries = listDzero->GetEntriesFast();
  if( !nEntries ) return;
  AliRDHFCutsD0toKpi *fCutsD0toKpi = (AliRDHFCutsD0toKpi*) fCutsPOI;
  if (fDebug) printf("  ᶫ%d candidates found. Looping...\n", nEntries);
  const Int_t ndght=2;
  Int_t nIDs[ndght];
  for( int iEntry=0; iEntry!=nEntries; ++iEntry ) {
    AliAODRecoDecayHF2Prong *D0 = (AliAODRecoDecayHF2Prong*) listDzero->UncheckedAt( iEntry );
    if( !D0 ) continue;
    if( !D0->HasSelectionBit(AliRDHFCuts::kD0toKpiCuts) ) continue;
    if( !fCutsD0toKpi->IsInFiducialAcceptance(D0->Pt(),D0->Y(421)) )continue;
    int topCut  = fCutsD0toKpi->IsSelected( D0, AliRDHFCuts::kAll, theAOD );
    int nLevel=topCut>0?1:0;
    if( (D0->Eta()<fPOIEta[0])||(D0->Eta()>fPOIEta[1]) )
      nLevel=0;
    Double_t MassD0=topCut>1?D0->InvMassD0bar():D0->InvMassD0();
    for(int h=0; h!=nLevel+1; ++h)
      fMass[h]->Fill(MassD0,D0->Pt());
    if( (D0->Pt()<fFlowPts[0]) || (D0->Pt()>fFlowPts[1]) ) continue;
    AliAODTrack* iT;
    for(Int_t i=0; i!=ndght; ++i) {
      iT = (AliAODTrack*)D0->GetDaughter(i);
      nIDs[i] = iT->GetID();
    }
    // Candidates Insertion (done in filling method: faster)
    if(nLevel)
      for(Int_t r=0; r!=5; ++r)
        if( (MassD0>=fFlowBands[0][r]) && (MassD0<fFlowBands[1][r]) ) {
          AliFlowCandidateTrack *sTrack = (AliFlowCandidateTrack*)
                                 MakeTrack(MassD0, D0->Pt(), D0->Phi(),
                                           D0->Eta(), ndght, nIDs);
          if(fDebug) printf("   ᶫInjecting D0 candidate on band %d \n", r);
          for(Int_t iDau=0; iDau!=ndght; ++iDau)
            for(Int_t iRPs=0; iRPs!=theMB[r]->NumberOfTracks(); ++iRPs) {
              AliFlowTrack *iRP = (AliFlowTrack*) (theMB[r]->GetTrack(iRPs));
              if(!iRP->InRPSelection()) continue;
              if( fabs(sTrack->GetIDDaughter(iDau)) == fabs(iRP->GetID()) ) {
                sTrack->SetDaughter(iDau,iRP);
                iRP->SetForRPSelection(kFALSE);
                if(fDebug) printf("    ᶫdaughter%d with fID %d was removed from this RP set\n", iDau, sTrack->GetIDDaughter(iDau));
              }
            }
          theMB[r]->AddTrack(sTrack);
        }
    // END of Candidates Insertion
  }
}
//______________________________________________________________________________
void AliAnalysisTaskFlowD2H::FillDStartoKpipi(AliAODEvent *theAOD, 
                                              AliFlowEvent *theMB[5] )
{
  TList *listHF = (TList*) theAOD->GetList();
  if(!listHF) return;
  TClonesArray *listDstar = (TClonesArray*) listHF->FindObject("Dstar");
  if(!listDstar) return;
  int nEntries = listDstar->GetEntriesFast();
  if( !nEntries ) return;
  AliRDHFCutsDStartoKpipi *fCutsDStartoKpipi = (AliRDHFCutsDStartoKpipi*) fCutsPOI;
  if (fDebug) printf("  ᶫ%d candidates found. Looping...\n", nEntries);
  const Int_t ndght=3;
  Int_t nIDs[ndght];
  for( int iEntry=0; iEntry!=nEntries; ++iEntry ) {
    AliAODRecoCascadeHF *DS = (AliAODRecoCascadeHF*) listDstar->UncheckedAt( iEntry );
    if( !DS ) continue;
    if( !DS->HasSelectionBit(AliRDHFCuts::kDstarCuts) ) continue;
    if( !fCutsDStartoKpipi->IsInFiducialAcceptance(DS->Pt(),DS->YDstar()) )continue;
    int topCut = fCutsDStartoKpipi->IsSelected( DS, AliRDHFCuts::kAll );
    int nLevel=topCut>0?1:0;
    if( (DS->Eta()<fPOIEta[0])||(DS->Eta()>fPOIEta[1]) )
      nLevel=0;
    Double_t MassDS=DS->DeltaInvMass();
    for(int h=0; h!=nLevel+1; ++h)
      fMass[h]->Fill(MassDS,DS->Pt());
    if( (DS->Pt()<fFlowPts[0]) || (DS->Pt()>fFlowPts[1]) ) continue;
    AliAODRecoDecayHF2Prong *D0 = (AliAODRecoDecayHF2Prong*)DS->Get2Prong();
    nIDs[0] = ((AliAODTrack*)D0->GetDaughter(0))->GetID();
    nIDs[1] = ((AliAODTrack*)D0->GetDaughter(1))->GetID();
    nIDs[2] = ((AliAODTrack*)DS->GetBachelor() )->GetID();
    // Candidates Insertion (done in filling method: faster)
    if(nLevel)
      for(Int_t r=0; r!=5; ++r)
        if( (MassDS>=fFlowBands[0][r]) && (MassDS<fFlowBands[1][r]) ) {
          AliFlowCandidateTrack *sTrack = (AliFlowCandidateTrack*)
                                 MakeTrack(MassDS, DS->Pt(), DS->Phi(),
                                           DS->Eta(), ndght, nIDs);
          if(fDebug) printf("   ᶫInjecting DStar candidate on band %d \n", r);
          for(Int_t iDau=0; iDau!=ndght; ++iDau)
            for(Int_t iRPs=0; iRPs!=theMB[r]->NumberOfTracks(); ++iRPs) {
              AliFlowTrack *iRP = (AliFlowTrack*) (theMB[r]->GetTrack(iRPs));
              if(!iRP->InRPSelection()) continue;
              if( fabs(sTrack->GetIDDaughter(iDau)) == fabs(iRP->GetID()) ) {
                sTrack->SetDaughter(iDau,iRP);
                iRP->SetForRPSelection(kFALSE);
                if(fDebug) printf("    ᶫdaughter%d with fID %d was removed from this RP set\n", iDau, sTrack->GetIDDaughter(iDau));
              }
            }
          theMB[r]->AddTrack(sTrack);
        }
    // END of Candidates Insertion
  }
}
//______________________________________________________________________________
void AliAnalysisTaskFlowD2H::FillDplustoKpipi(AliAODEvent *theAOD, 
                                              AliFlowEvent *theMB[5] )
{
  TList *listHF = (TList*) theAOD->GetList();
  if(!listHF) return;
  TClonesArray *listDplus = (TClonesArray*) listHF->FindObject("Charm3Prong");
  if(!listDplus) return;
  int nEntries = listDplus->GetEntriesFast();
  if( !nEntries ) return;
  AliRDHFCutsDplustoKpipi *fCutsDStartoKpipi = (AliRDHFCutsDplustoKpipi*) fCutsPOI;
  if (fDebug) printf("  ᶫ%d candidates found. Looping...\n", nEntries);
  const Int_t ndght=3;
  Int_t nIDs[ndght];
  for( int iEntry=0; iEntry!=nEntries; ++iEntry ) {
    AliAODRecoDecayHF3Prong *Dp = (AliAODRecoDecayHF3Prong*) listDplus->UncheckedAt( iEntry );
    if( !Dp ) continue;
    if( !Dp->HasSelectionBit(AliRDHFCuts::kDplusCuts) ) continue;
    if( !fCutsDStartoKpipi->IsInFiducialAcceptance(Dp->Pt(),Dp->YDplus()) )continue;
    int topCut = fCutsDStartoKpipi->IsSelected( Dp, AliRDHFCuts::kAll );
    int nLevel=topCut>0?1:0;
    if( (Dp->Eta()<fPOIEta[0])||(Dp->Eta()>fPOIEta[1]) )
      nLevel=0;
    Double_t MassDp=Dp->InvMassDplus();
    for(int h=0; h!=nLevel+1; ++h)
      fMass[h]->Fill(MassDp,Dp->Pt());
    if( (Dp->Pt()<fFlowPts[0]) || (Dp->Pt()>fFlowPts[1]) ) continue;
    nIDs[0] = ((AliAODTrack*)Dp->GetDaughter(0))->GetID();
    nIDs[1] = ((AliAODTrack*)Dp->GetDaughter(1))->GetID();
    nIDs[2] = ((AliAODTrack*)Dp->GetDaughter(2))->GetID();
    // Candidates Insertion (done in filling method: faster)
    if(nLevel)
      for(Int_t r=0; r!=5; ++r)
        if( (MassDp>=fFlowBands[0][r]) && (MassDp<fFlowBands[1][r]) ) {
          AliFlowCandidateTrack *sTrack = (AliFlowCandidateTrack*)
                                 MakeTrack(MassDp, Dp->Pt(), Dp->Phi(),
                                           Dp->Eta(), ndght, nIDs);
          if(fDebug) printf("   ᶫInjecting DStar candidate on band %d \n", r);
          for(Int_t iDau=0; iDau!=ndght; ++iDau)
            for(Int_t iRPs=0; iRPs!=theMB[r]->NumberOfTracks(); ++iRPs) {
              AliFlowTrack *iRP = (AliFlowTrack*) (theMB[r]->GetTrack(iRPs));
              if(!iRP->InRPSelection()) continue;
              if( fabs(sTrack->GetIDDaughter(iDau)) == fabs(iRP->GetID()) ) {
                sTrack->SetDaughter(iDau,iRP);
                iRP->SetForRPSelection(kFALSE);
                if(fDebug) printf("    ᶫdaughter%d with fID %d was removed from this RP set\n", iDau, sTrack->GetIDDaughter(iDau));
              }
            }
          theMB[r]->AddTrack(sTrack);
        }
    // END of Candidates Insertion
  }
}
//______________________________________________________________________________
AliFlowCandidateTrack* AliAnalysisTaskFlowD2H::MakeTrack( Double_t mass, 
                          Double_t pt, Double_t phi, Double_t eta, 
                          Int_t nDau, Int_t *iID ) {
  AliFlowCandidateTrack *sTrack = new AliFlowCandidateTrack();
  sTrack->SetMass(mass);
  sTrack->SetPt(pt);
  sTrack->SetPhi(phi);
  sTrack->SetEta(eta);
  for(Int_t iDau=0; iDau!=nDau; ++iDau)
    sTrack->AddDaughter(iID[iDau]);
  sTrack->SetForPOISelection(kTRUE);
  sTrack->SetForRPSelection(kFALSE);
  return sTrack;
}
//_____________________________________________________________________________
void AliAnalysisTaskFlowD2H::Terminate(Option_t *)
{

}

