/*************************************************************************
* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
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

/////////////////////////////////////////////////////
// AliAnalysisTaskFlowK0Candidates:
// Analysis task to select K0 candidates for flow analysis.
// Uses one AliESDtrackCuts object for both daughters and
// QA histograms to monitor the reconstruction.
// Author: Carlos Perez (cperez@cern.ch)
//////////////////////////////////////////////////////

#include "TChain.h"
#include "TList.h"
#include "TH1D.h"
#include "TVector3.h"

#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"

#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"

#include "AliAODEvent.h"
#include "AliAODInputHandler.h"
#include "AliAODTrack.h"

#include "AliCFManager.h"

#include "TMath.h"
#include "TObjArray.h"
#include "AliFlowCandidateTrack.h"
#include "AliFlowEventCuts.h"

#include "AliAnalysisTaskFlowK0Candidates.h"


ClassImp(AliAnalysisTaskFlowK0Candidates)

//_____________________________________________________________________________
AliAnalysisTaskFlowK0Candidates::AliAnalysisTaskFlowK0Candidates() :
  AliAnalysisTaskSE(),
  fCutsEvent(NULL),
  fCuts(NULL),
  fQAList(NULL),
  fMassMin(0),
  fMassMax(0),
  fDLcut(0),
  fEvent(NULL),
  fMulti(NULL)
{
  //
  for(int i=0; i!=4; ++i) {
    fMass[i] = NULL;
    fDCA[i]  = NULL;
    fDL[i]   = NULL;
    fCTP[i]  = NULL;
    fd0d0[i] = NULL;
    fPhi[i]  = NULL;
    fEta[i]  = NULL;
    fPt[i]   = NULL;
    fAPhi[i] = NULL;
    fAEta[i] = NULL;
    fAPt[i]  = NULL;
    fBPhi[i] = NULL;
    fBEta[i] = NULL;
    fBPt[i]  = NULL;
  }
}

//_____________________________________________________________________________
AliAnalysisTaskFlowK0Candidates::AliAnalysisTaskFlowK0Candidates(const char *name, AliFlowEventCuts *cutsEvent, AliESDtrackCuts *cuts, Double_t MassMin, Double_t MassMax) :
  AliAnalysisTaskSE(name),
  fCutsEvent(cutsEvent),
  fCuts(cuts),
  fQAList(NULL),
  fMassMin(MassMin),
  fMassMax(MassMax),
  fDLcut(0),
  fEvent(NULL),
  fMulti(NULL)
{
  //
  for(int i=0; i!=4; ++i) {
    fMass[i] = NULL;
    fDCA[i]  = NULL;
    fDL[i]   = NULL;
    fCTP[i]  = NULL;
    fd0d0[i] = NULL;
    fPhi[i]  = NULL;
    fEta[i]  = NULL;
    fPt[i]   = NULL;
    fAPhi[i] = NULL;
    fAEta[i] = NULL;
    fAPt[i]  = NULL;
    fBPhi[i] = NULL;
    fBEta[i] = NULL;
    fBPt[i]  = NULL;
  }
  DefineInput(  0, TChain::Class() );
  DefineOutput( 1, TObjArray::Class() );
  DefineOutput( 2, TList::Class() );
}

//_____________________________________________________________________________
AliAnalysisTaskFlowK0Candidates::~AliAnalysisTaskFlowK0Candidates()
{
  //
  if(fQAList) delete fQAList;
}

//_____________________________________________________________________________
void AliAnalysisTaskFlowK0Candidates::UserCreateOutputObjects()
{
  //
  fQAList = new TList();
  fQAList->SetOwner();
  AddQAEvents();
  AddQACandidates();
  PostData( 2, fQAList);
}
//_____________________________________________________________________________
void AliAnalysisTaskFlowK0Candidates::AddQAEvents()
{
  //
  TList *tQAEvents = new TList();
  tQAEvents->SetName("Events");
  tQAEvents->SetOwner();
  fEvent = new TH1D("Event", "Number of Events",    2, 0, 2);
  tQAEvents->Add(fEvent);
  fMulti = new TH1D("Multiplicity", "Multiplicity", 180, 0, 10000);
  tQAEvents->Add(fMulti);
  fQAList->Add(tQAEvents);
}  
//_____________________________________________________________________________
void AliAnalysisTaskFlowK0Candidates::AddQACandidates()
{
  //
  TList *tQACandidates[4];
  TList *tQADaughters[4];
  for(int i=0; i!=4; ++i) {
    tQACandidates[i] = new TList();
    tQACandidates[i]->SetOwner();
    tQACandidates[i]->SetName(Form("Candidates%d",i));
    fMass[i] = new TH1D( Form("Mass%i",i), "Mass;M_{#pi#pi} [GeV];Counts per MeV", 180, 0.41, 0.59); tQACandidates[i]->Add( fMass[i] );
    fDCA[i]  = new TH1D( Form("DCA%i" ,i), "DCA;[cm];Counts per 10 um",            180, 0.00, 0.18); tQACandidates[i]->Add( fDCA[i]  );
    fDL[i]   = new TH1D( Form("DL%i"  ,i), "Decay Length;[cm];Counts per 0.1 mm",  180, 0.00, 1.80); tQACandidates[i]->Add( fDL[i]   );
    fCTP[i]  = new TH1D( Form("CTP%i" ,i), "Cos#theta_{p}",                        180,-1.10, 1.10); tQACandidates[i]->Add( fCTP[i]  );
    fd0d0[i] = new TH1D( Form("d0d0%i",i), "d_{0}xd_{0};[cm^{2}];Cnts 0.01 mm^{2}",180,-0.009,0.009);tQACandidates[i]->Add( fd0d0[i] );
    fPhi[i]  = new TH1D( Form("Phi%i" ,i), "Phi;[rad];Counts per degree",     180,0,TMath::TwoPi()); tQACandidates[i]->Add( fPhi[i]  );
    fEta[i]  = new TH1D( Form("Eta%i" ,i), "Eta;;Counts per 0.04",                 180,-3.60, 3.60); tQACandidates[i]->Add( fEta[i]  );
    fPt[i]   = new TH1D( Form("Pt%i"  ,i), "Pt;[GeV];Counts per 0.1 GeV",          180, 0.00,18.00); tQACandidates[i]->Add( fPt[i]   );
    tQADaughters[i] = new TList();
    tQADaughters[i]->SetOwner();
    tQADaughters[i]->SetName(Form("Daughters%d",i));
    fAPhi[i] = new TH1D( Form("PhiBef%i",i), "Phi prePropagation;[rad];Cnts per degree",180,0,TMath::TwoPi()); tQADaughters[i]->Add( fAPhi[i] );
    fAEta[i] = new TH1D( Form("EtaBef%i",i), "Eta prePropagation;;Counts per 0.04"     ,180,-3.6,3.6); tQADaughters[i]->Add( fAEta[i] );
    fAPt[i]  = new TH1D( Form("PtBef%i" ,i), "Pt prePropagation;[GeV];Counts per 0.1 GeV",180, 0, 18); tQADaughters[i]->Add( fAPt[i]  );
    fBPhi[i] = new TH1D( Form("PhiAft%i",i), "Phi posPropagation;[rad];Cnts per degree",180,0,TMath::TwoPi()); tQADaughters[i]->Add( fBPhi[i] );
    fBPhi[i] = new TH1D( Form("EtaAft%i",i), "Eta posPropagation;;Counts per 0.04"     ,180,-3.6,3.6); tQADaughters[i]->Add( fBPhi[i] );
    fBPt[i]  = new TH1D( Form("PtAft%i" ,i), "Pt posPropagation;[GeV];Counts per 0.1 GeV",180, 0, 18); tQADaughters[i]->Add( fBPt[i]  );
    tQACandidates[i]->Add(tQADaughters[i]);
    fQAList->Add(tQACandidates[i]);
  }
}
//_____________________________________________________________________________
void AliAnalysisTaskFlowK0Candidates::NotifyRun()
{
  //
//  AliESDEvent *fESD = dynamic_cast<AliESDEvent*>(InputEvent());
//  if(!fESD) return;
}

//_____________________________________________________________________________
void AliAnalysisTaskFlowK0Candidates::UserExec(Option_t *)
{
  //
  AliESDEvent *fESD = dynamic_cast<AliESDEvent*>(InputEvent());
  AliAODEvent *fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
  if( (!fESD)&&(!fAOD) )
    return;
//  printf("\nEvent found");
  fEvent->Fill( 0 );
  if(fCutsEvent)
    if( !(fCutsEvent->IsSelected(InputEvent())) )
      return;
  fEvent->Fill( 1 );
  if(fESD)
    ReadFromESD(fESD);
  else if(fAOD)
    ReadFromAOD(fAOD);
}
//_____________________________________________________________________________
void AliAnalysisTaskFlowK0Candidates::ReadFromESD(const AliESDEvent *fESD)
{
  //
  TObjArray *gPOIselection =  new TObjArray();
  gPOIselection->SetOwner();
  int nTracks = fESD->GetNumberOfTracks();
  fMulti->Fill( nTracks );
  for(int i=0; i!=nTracks; ++i) {
    AliESDtrack *ioT = fESD->GetTrack(i);
    if(fCuts)
      if( !(fCuts->IsSelected(ioT)) )
        continue;
    // printf("\n first particle OK...");
    for(int j=i+1; j!=nTracks; ++j) {
      AliESDtrack *joT = fESD->GetTrack(j);
      if( (ioT->Charge()*joT->Charge()) > 0 )
        continue;
      if(fCuts)
        if( !(fCuts->IsSelected(joT)) )
          continue;
      // printf("\n second also...");
      AliESDtrack *iT = new AliESDtrack(*ioT);
      AliESDtrack *jT = new AliESDtrack(*joT);
      // getting distance of closest approach
      double DCA = iT->PropagateToDCA(jT,fESD->GetMagneticField());
      // printf(" propagated...");
      // getting decay length
      TVector3 vp, vv, vl;
      vp = TVector3( (iT->Xv()+jT->Xv())/2, (iT->Yv()+jT->Yv())/2, (iT->Zv()+jT->Zv())/2 ); // candidate position
      vv = TVector3( fESD->GetPrimaryVertex()->GetX(), fESD->GetPrimaryVertex()->GetY(), fESD->GetPrimaryVertex()->GetZ() ); // vertex position
      vl = vp - vv; // flight line
      double gDL;
      gDL  = vl.Mag();
      // getting cos pointing angle
      double gCTP;
      TVector3 vi, vj, vs;
      vi = TVector3( iT->Px(), iT->Py(), iT->Pz() );
      vj = TVector3( jT->Px(), jT->Py(), jT->Pz() );
      vs = vi + vj;
      gCTP = TMath::Cos( vl.Angle(vs) ); // Marta says: "Carlos, it is faster if you calculate this by hand!"
      double gD0D0;
      gD0D0 = iT->GetD(vv.X(),vv.Y(),fESD->GetMagneticField())*jT->GetD(vv.X(),vv.Y(),fESD->GetMagneticField());
      // getting invariant mass
      double sum12 = iT->P()*iT->P()+jT->P()*jT->P();
      double pro12 = iT->P()*iT->P()*jT->P()*jT->P();
      double gInvMass = 0;
      gInvMass += TMath::Power(0.13957018,2);
      gInvMass += sqrt( TMath::Power(0.13957018,4) + TMath::Power(0.13957018,2)*(sum12) + pro12 );
      gInvMass -= ( iT->Px()*jT->Px()+iT->Py()*jT->Py()+iT->Pz()*jT->Pz() );
      gInvMass *= 2;
      gInvMass = sqrt(gInvMass);
      // filtering
      int iLevel = 3;
      if( gCTP<0.8 )
        iLevel = 2;
      if( gDL<fDLcut ) // 0.5
        iLevel = 1;
      if( DCA>0.05 )
        iLevel = 0;
      // printf(" candidate at level %d...",iLevel);
      for(int h=0; h!=iLevel+1; ++h) {
        // printf(" %d",h);
        // candidate
        fDCA[h]->Fill( DCA );
        fDL[h]->Fill( gDL );
        fCTP[h]->Fill( gCTP );
        fd0d0[h]->Fill( gD0D0 );
        fMass[h]->Fill( gInvMass );
        fPhi[h]->Fill( vs.Phi()+TMath::Pi() );
        fEta[h]->Fill( vs.Eta() );
        fPt[h]->Fill( vs.Pt() );
        // daughters
        fAPhi[h]->Fill( iT->Phi() );
        fAEta[h]->Fill( iT->Eta() );
        fAPt[h]->Fill(  iT->Pt()  );
        fAPhi[h]->Fill( jT->Phi() );
        fAEta[h]->Fill( jT->Eta() );
        fAPt[h]->Fill(  jT->Pt()  );
        fBPhi[h]->Fill( ioT->Phi());
        fBPhi[h]->Fill( ioT->Eta());
        fBPt[h]->Fill(  ioT->Pt() );
        fBPhi[h]->Fill( joT->Phi());
        fBPhi[h]->Fill( joT->Eta());
        fBPt[h]->Fill(  joT->Pt() );
      }
      // building candidate
      if( (gInvMass>fMassMin) && (gInvMass<fMassMax) && (iLevel==3) ) {
        AliFlowCandidateTrack *sTrack = new AliFlowCandidateTrack();
        sTrack->SetMass( gInvMass );
        sTrack->SetPt( vs.Pt() );
        sTrack->SetPhi( vs.Phi()+TMath::Pi() );
        sTrack->SetEta( vs.Eta() );
        sTrack->SetCharge( 0 );
        sTrack->AddDaughter( iT->GetID() );
        sTrack->AddDaughter( jT->GetID() );
        gPOIselection->AddLast( sTrack );
      }
      delete iT;
      delete jT;
    } // endfor j
  } // endfor i
  PostData( 1, gPOIselection );
  PostData( 2, fQAList);
}

//_____________________________________________________________________________
void AliAnalysisTaskFlowK0Candidates::ReadFromAOD(const AliAODEvent *fAOD)
{
  //
  TObjArray *gPOIselection =  new TObjArray();
  gPOIselection->SetOwner();
  int nTracks = fAOD->GetNumberOfTracks();
  fMulti->Fill( nTracks );
  for(int i=0; i!=fAOD->GetNumberOfV0s(); ++i) {
    AliAODv0 *iV0 = (AliAODv0*) fAOD->GetV0(i);
    if ( iV0->GetOnFlyStatus() ) continue;
    if ( (iV0->ChargeProng(0)*iV0->ChargeProng(1))>0 ) continue;
    // getting distance of closest approach
    double DCA = iV0->DcaV0Daughters();
    // getting decay length
    double gDL;
    double vv[3];
    vv[0] = fAOD->GetPrimaryVertex()->GetX();
    vv[1] = fAOD->GetPrimaryVertex()->GetY();
    vv[2] = fAOD->GetPrimaryVertex()->GetZ();
    gDL = iV0->DecayLengthV0(vv);
    // cos pointing angle
    double gCTP;
    gCTP = iV0->CosPointingAngle(vv);
    double gD0D0;
    gD0D0 = iV0->Prodd0d0();
    // getting invariant mass
    double gInvMass = iV0->MassK0Short();
    // filtering
    int iLevel = 3;
    if( gCTP<0.8 )
      iLevel = 2;
    if( gDL<fDLcut ) // 0.5
      iLevel = 1;
    if( DCA>0.05 )
      iLevel = 0;
    for(int h=0; h!=iLevel+1; ++h) {
      // candidate
      fDCA[h]->Fill( DCA );
      fDL[h]->Fill( gDL );
      fCTP[h]->Fill( gCTP );
      fd0d0[h]->Fill( gD0D0 );
      fMass[h]->Fill( gInvMass );
      fPhi[h]->Fill( iV0->Phi() );
      fEta[h]->Fill( iV0->Eta() );
      fPt[h]->Fill( iV0->Pt() );
      // daughters
      fAPhi[h]->Fill( ( (AliAODTrack*) iV0->GetDaughter(0) )->Phi() );
      fAEta[h]->Fill( ( (AliAODTrack*) iV0->GetDaughter(0) )->Eta() );
      fAPt[h]->Fill(  ( (AliAODTrack*) iV0->GetDaughter(0) )->Pt()  );
      fAPhi[h]->Fill( ( (AliAODTrack*) iV0->GetDaughter(1) )->Phi() );
      fAEta[h]->Fill( ( (AliAODTrack*) iV0->GetDaughter(1) )->Eta() );
      fAPt[h]->Fill(  ( (AliAODTrack*) iV0->GetDaughter(1) )->Pt()  );
      fBPhi[h]->Fill( iV0->PhiProng(0) );
      fBPhi[h]->Fill( iV0->EtaProng(0) );
      fBPt[h]->Fill(  iV0->PtProng(0)  );
      fBPhi[h]->Fill( iV0->PhiProng(1) );
      fBPhi[h]->Fill( iV0->EtaProng(1) );
      fBPt[h]->Fill(  iV0->PtProng(1)  );
    }
    AliFlowCandidateTrack *sTrack = new AliFlowCandidateTrack();
    sTrack->SetMass( gInvMass );
    sTrack->SetPt( iV0->Pt() );
    sTrack->SetPhi( iV0->Phi() );
    sTrack->SetEta( iV0->Eta() );
    sTrack->SetCharge( 0 );
    sTrack->AddDaughter( iV0->GetPosID() );
    sTrack->AddDaughter( iV0->GetNegID() );
    gPOIselection->AddLast( sTrack );
  }
  PostData( 1, gPOIselection );
  PostData( 2, fQAList);
}

//_____________________________________________________________________________
void AliAnalysisTaskFlowK0Candidates::Terminate(Option_t *)
{
  //
}

