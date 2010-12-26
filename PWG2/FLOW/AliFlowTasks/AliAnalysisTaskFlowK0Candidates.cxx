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
  fList(NULL),
  fMass(NULL),
  fDCA(NULL),
  fDL(NULL),
  fCTP(NULL),
  fMassF(NULL),
  fDCAF(NULL),
  fDLF(NULL),
  fCTPF(NULL),
  fMassMin(0),
  fMassMax(0)
{
}

//_____________________________________________________________________________
AliAnalysisTaskFlowK0Candidates::AliAnalysisTaskFlowK0Candidates(const char *name, AliFlowEventCuts *cutsEvent, AliESDtrackCuts *cuts, Double_t MassMin, Double_t MassMax) :
  AliAnalysisTaskSE(name),
  fCutsEvent(cutsEvent),
  fCuts(cuts),
  fList(NULL),
  fMass(NULL),
  fDCA(NULL),
  fDL(NULL),
  fCTP(NULL),
  fMassF(NULL),
  fDCAF(NULL),
  fDLF(NULL),
  fCTPF(NULL),
  fMassMin(MassMin),
  fMassMax(MassMax)
{
  DefineInput(  0, TChain::Class() );
  DefineOutput( 1, TObjArray::Class() );
  DefineOutput( 2, TList::Class() );
}

//_____________________________________________________________________________
AliAnalysisTaskFlowK0Candidates::~AliAnalysisTaskFlowK0Candidates()
{
  if(fList) {
    fList->Delete();
    delete fList;
  }
}

//_____________________________________________________________________________
void AliAnalysisTaskFlowK0Candidates::UserCreateOutputObjects()
{
  fList = new TList();
  fList->SetOwner();
  fMass = new TH1D( "fMass","Mass;M_{#pi#pi} [GeV];Counts per MeV", 180, 0.41, 0.59); fList->Add( fMass );
  fDCA  = new TH1D( "fDCA" ,"DCA;[cm];Counts per 10 um",            180, 0.00, 0.18); fList->Add( fDCA  );
  fDL   = new TH1D( "fDL"  ,"Decay Length;[cm];Counts per 0.1 mm",  180, 0.00, 1.80); fList->Add( fDL   );
  fCTP  = new TH1D( "fCTP" ,"Cos#theta_{p}",                        180,-1.10, 1.10); fList->Add( fCTP  );

  fMassF= new TH1D( "fMassF","Mass after Cuts;M_{#pi#pi} [GeV];Counts per MeV", 180, 0.41, 0.59); fList->Add( fMassF );
  fDCAF = new TH1D( "fDCAF" ,"DCA after Cuts;[cm];Counts per 10 um",            180, 0.00, 0.18); fList->Add( fDCAF  );
  fDLF  = new TH1D( "fDLF"  ,"Decay Length after Cuts;[cm];Counts per 0.1 mm",  180, 0.00, 1.80); fList->Add( fDLF   );
  fCTPF = new TH1D( "fCTPF" ,"Cos#theta_{p} after Cuts",                        180,-1.10, 1.10); fList->Add( fCTPF  );

  PostData( 2, fList);
}


//_____________________________________________________________________________
void AliAnalysisTaskFlowK0Candidates::NotifyRun()
{
  AliESDEvent *fESD = dynamic_cast<AliESDEvent*>(InputEvent());
  if(!fESD) return;
}

//_____________________________________________________________________________
void AliAnalysisTaskFlowK0Candidates::UserExec(Option_t *)
{
  AliESDEvent *fESD = dynamic_cast<AliESDEvent*>(fInputEvent);
  if(!fESD) return;
  if(fCutsEvent)
    if( !(fCutsEvent->IsSelected(InputEvent())) ) return;

  TObjArray *POIselection =  new TObjArray();
  POIselection->SetOwner();

  int nTracks = fESD->GetNumberOfTracks();
  for(int i=0; i!=nTracks; ++i) {  // first particle
    AliESDtrack *ioT = fESD->GetTrack(i);
    if(fCuts)
      if( !(fCuts->IsSelected(ioT)) )
        continue;
    for(int j=i+1; j!=nTracks; ++j) {  // second particle
      AliESDtrack *joT = fESD->GetTrack(j);
      if( (ioT->Charge()*joT->Charge()) > 0 ) // only keeping opposite signs
        continue;
      if(fCuts)
        if( !(fCuts->IsSelected(joT)) )
          continue;
      AliESDtrack *iT = new AliESDtrack(*ioT);
      AliESDtrack *jT = new AliESDtrack(*joT);
      // getting distance of closest approach
      double DCA = iT->PropagateToDCA(jT,fESD->GetMagneticField());
      fDCA->Fill( DCA );
      // getting decay length
      double DL;
      double vx = fESD->GetPrimaryVertex()->GetX();
      double vy = fESD->GetPrimaryVertex()->GetY();
      double vz = fESD->GetPrimaryVertex()->GetZ();
      double px = (iT->Xv()+jT->Xv())/2;
      double py = (iT->Yv()+jT->Yv())/2;
      double pz = (iT->Zv()+jT->Zv())/2;
      DL  = (vx-px)*(vx-px);
      DL += (vy-py)*(vy-py);
      DL += (vz-pz)*(vz-pz);
      DL = sqrt(DL);
      fDL->Fill( DL );
      // cos pointing angle
      double CTP;
      TVector3 vi, vj, vs, vp;
      vi = TVector3( iT->Px(), iT->Py(), iT->Pz() );
      vj = TVector3( jT->Px(), jT->Py(), jT->Pz() );
      vs = vi + vj;
      vp = TVector3( px, py, pz );
      CTP = TMath::Cos( vp.Angle(vs) );
      fCTP->Fill(CTP);

      // getting invariant mass
      double sum12 = iT->P()*iT->P()+jT->P()*jT->P();
      double pro12 = iT->P()*iT->P()*jT->P()*jT->P();
      double InvMass = 0;
      InvMass += TMath::Power(0.13957018,2);
      InvMass += sqrt( TMath::Power(0.13957018,4) + TMath::Power(0.13957018,2)*(sum12) + pro12 );
      InvMass -= ( iT->Px()*jT->Px()+iT->Py()*jT->Py()+iT->Pz()*jT->Pz() );
      InvMass *= 2;
      InvMass = sqrt(InvMass);
      fMass->Fill( InvMass );

      // evaluating cuts
      bool bW   = (InvMass>fMassMin) && (InvMass<fMassMax);
      bool bCTP = (fabs(CTP)>0.95);
      bool bDCA = (DCA<0.05);
      bool bDL  = (DL>2.5);
      if( (!bW) || (!bDCA) || (!bDL) || (!bCTP) ) {
        delete iT;
        delete jT;
        continue;
      }
      fDCAF->Fill( DCA );
      fDLF->Fill( DL );
      fCTPF->Fill(CTP);
      fMassF->Fill( InvMass );

      AliFlowCandidateTrack *sTrack = new AliFlowCandidateTrack();
      //booking pt of selected candidates
      sTrack->SetMass( InvMass );
      sTrack->SetPt( vs.Pt() );
      sTrack->SetPhi( vs.Phi() );
      sTrack->SetEta( vs.Eta() );
      sTrack->SetCharge( 0 );
      sTrack->AddDaughter( i );
      sTrack->AddDaughter( j );
      //saving selected candidate
      POIselection->AddLast( sTrack );
      delete iT;
      delete jT;
    }
  }
  PostData( 1, POIselection );
  PostData( 2, fList);
}

//_____________________________________________________________________________
void AliAnalysisTaskFlowK0Candidates::Terminate(Option_t *)
{
}

