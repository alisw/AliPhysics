/*************************************************************************
* Copyright(c) 1998-2008,ALICE Experiment at CERN,All rights reserved. *
*                                                                        *
* Author: The ALICE Off-line Project.                                    *
* Contributors are mentioned in the code where appropriate.              *
*                                                                        *
* Permission to use,copy,modify and distribute this software and its   *
* documentation strictly for non-commercial purposes is hereby granted   *
* without fee,provided that the above copyright notice appears in all   *
* copies and that both the copyright notice and this permission notice   *
* appear in the supporting documentation. The authors make no claims     *
* about the suitability of this software for any purpose. It is          *
* provided "as is" without express or implied warranty.                  *
**************************************************************************/

/////////////////////////////////////////////////////
// AliAnalysisTaskFlowStrange:
// Analysis task to select K0/Lambda candidates for flow analysis.
// Uses one AliESDtrackCuts object for both daughters and
// QA histograms to monitor the reconstruction.
// Authors: Cristian Ivan (civan@cern.ch)
//          Carlos Perez (cperez@cern.ch)
//          Pawel Debski (pdebski@cern.ch)
//////////////////////////////////////////////////////

#include "TChain.h"
#include "TList.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TProfile.h"
#include "TVector3.h"

#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"

#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliESDVertex.h"
#include "AliESDv0.h"
#include "AliESDtrackCuts.h"

#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include "AliAODVertex.h"
#include "AliAODv0.h"

#include "TMath.h"
#include "TObjArray.h"
#include "AliFlowCandidateTrack.h"

#include "AliFlowTrackCuts.h"
#include "AliFlowEventCuts.h"
#include "AliFlowEvent.h"
#include "AliFlowCommonConstants.h"

#include "AliAnalysisTaskFlowStrange.h"

ClassImp(AliAnalysisTaskFlowStrange)

//=======================================================================
AliAnalysisTaskFlowStrange::AliAnalysisTaskFlowStrange() :
  AliAnalysisTaskSE(),
  fPIDResponse(NULL),
  fDebug(kFALSE),
  fSpecie(0),
  fMassBins(0),
  fMinMass(0.0),
  fMaxMass(0.0),
  fCutsEvent(NULL),
  fCutsRFPTPC(NULL),
  fCutsRFPVZE(NULL),
  fCutsPOI(NULL),
  fCutsDau(NULL),
  fFlowEventTPC(NULL),
  fFlowEventVZE(NULL),
  fCandidates(NULL),
  fQAList(NULL)
{
  //ctor
  for (Int_t i=0; i!=9; ++i)
    fV0Cuts[i] = 0;
}
//=======================================================================
AliAnalysisTaskFlowStrange::AliAnalysisTaskFlowStrange(const char *name,
						       AliFlowEventCuts *cutsEvent,
						       AliFlowTrackCuts *cutsRFPTPC,
						       AliFlowTrackCuts *cutsRFPVZE,
						       AliESDtrackCuts *cutsDau) :
  AliAnalysisTaskSE(name),
  fPIDResponse(NULL),
  fDebug(kFALSE),
  fSpecie(0),
  fMassBins(0),
  fMinMass(0.0),
  fMaxMass(0.0),
  fCutsEvent(cutsEvent),
  fCutsRFPTPC(cutsRFPTPC),
  fCutsRFPVZE(cutsRFPVZE),
  fCutsPOI(NULL),
  fCutsDau(cutsDau),
  fFlowEventTPC(NULL),
  fFlowEventVZE(NULL),
  fCandidates(NULL),
  fQAList(NULL)
{
  //ctor
  for (Int_t i=0; i!=9; ++i)
    fV0Cuts[i] = 0;

  DefineInput( 0,TChain::Class());
  DefineOutput(1,AliFlowEventSimple::Class()); // TPC object
  DefineOutput(2,AliFlowEventSimple::Class()); // VZE object
  DefineOutput(3,TList::Class());
}
//=======================================================================
AliAnalysisTaskFlowStrange::~AliAnalysisTaskFlowStrange()
{
  //dtor
  if (fQAList)       delete fQAList;
  if (fFlowEventTPC) delete fFlowEventTPC;
  if (fFlowEventVZE) delete fFlowEventVZE;
  if (fCandidates)   delete fCandidates;
  if (fCutsDau)      delete fCutsDau;
  if (fCutsPOI)      delete fCutsPOI;
  if (fCutsRFPTPC)   delete fCutsRFPTPC;
  if (fCutsRFPVZE)   delete fCutsRFPVZE;
}
//=======================================================================
void AliAnalysisTaskFlowStrange::UserCreateOutputObjects()
{
  //UserCreateOutputObjects
  fQAList=new TList();
  fQAList->SetOwner();
  AddQAEvents();
  AddQACandidates();

  AliFlowCommonConstants *cc = AliFlowCommonConstants::GetMaster();
  cc->SetNbinsMult(1); cc->SetMultMin(0);   cc->SetMultMax(1);
  cc->SetNbinsPt(60);  cc->SetPtMin(0.0);   cc->SetPtMax(12.0);
  cc->SetNbinsPhi(1);  cc->SetPhiMin(0.0);  cc->SetPhiMax(TMath::TwoPi());
  cc->SetNbinsEta(1);  cc->SetEtaMin(-2.0); cc->SetEtaMax(+2.0);
  cc->SetNbinsQ(1);    cc->SetQMin(0.0);    cc->SetQMax(1.0);
  cc->SetNbinsMass(fMassBins);
  cc->SetMassMin(fMinMass);
  cc->SetMassMax(fMaxMass);

  fCutsPOI = new AliFlowTrackCuts("dumb_cuts");
  fCutsPOI->SetParamType( fCutsRFPTPC->GetParamType() );
  fCutsPOI->SetPtRange(+1.0,-1.0);
  fCutsPOI->SetEtaRange(+1.0,-1.0);

  fFlowEventTPC = new AliFlowEvent(3000);
  fFlowEventVZE = new AliFlowEvent(1000);
  fCandidates = new TObjArray(100);
  fCandidates->SetOwner();

  AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
  fPIDResponse = inputHandler->GetPIDResponse();

  PostData(1,fFlowEventTPC);
  PostData(2,fFlowEventVZE);
  PostData(3,fQAList);

}
//=======================================================================
void AliAnalysisTaskFlowStrange::AddQAEvents()
{
  // function to add event qa
  TList *tQAEvents=new TList();
  tQAEvents->SetName("Event");
  tQAEvents->SetOwner();

  TH1D *tEvent = new TH1D("Events","Number of Events",3,0,3); tQAEvents->Add(tEvent);
  tEvent->GetXaxis()->SetBinLabel(1,"reached");
  tEvent->GetXaxis()->SetBinLabel(2,"selected");
  tEvent->GetXaxis()->SetBinLabel(3,"unexpected");

  TH1D *tTPCRFP = new TH1D("RFPTPC","TPC Reference Flow Particles;multiplicity",100,0,3000); tQAEvents->Add(tTPCRFP);
  TH1D *tVZERFP = new TH1D("RFPVZE","VZERO Reference Flow Particles;multiplicity",100,0,30000); tQAEvents->Add(tVZERFP);

  TProfile *tCuts = new TProfile("Cuts","Analysis Cuts",10,0,10);
  tCuts->Fill(0.5,fV0Cuts[0],1); tCuts->GetXaxis()->SetBinLabel(1,"dl");
  tCuts->Fill(1.5,fV0Cuts[1],1); tCuts->GetXaxis()->SetBinLabel(2,"dca");
  tCuts->Fill(2.5,fV0Cuts[2],1); tCuts->GetXaxis()->SetBinLabel(3,"ctp");
  tCuts->Fill(3.5,fV0Cuts[3],1); tCuts->GetXaxis()->SetBinLabel(4,"d0");
  tCuts->Fill(4.5,fV0Cuts[4],1); tCuts->GetXaxis()->SetBinLabel(5,"d0xd0");
  tCuts->Fill(5.5,fV0Cuts[5],1); tCuts->GetXaxis()->SetBinLabel(6,"qt");
  tCuts->Fill(6.5,fV0Cuts[6],1); tCuts->GetXaxis()->SetBinLabel(7,"min eta");
  tCuts->Fill(7.5,fV0Cuts[7],1); tCuts->GetXaxis()->SetBinLabel(8,"max eta");
  tCuts->Fill(8.5,fV0Cuts[8],1); tCuts->GetXaxis()->SetBinLabel(9,"use PID");
  tQAEvents->Add(tCuts);
  fQAList->Add(tQAEvents);
}
//=======================================================================
void AliAnalysisTaskFlowStrange::AddQACandidates()
{
  // function to add histogramming for candidates
  TList *tQACuts;
  TH2D *tDL, *tDCA, *tCTP, *tD0, *tD0D0;
  TH3D *tAP;

  tQACuts = new TList(); tQACuts->SetOwner(); tQACuts->SetName("QACutsBefore");
  tDL   = new TH2D("BefDL",  "DL;[cm];Pt [GeV]",  50,0,10,  24,0,12); tQACuts->Add(tDL);
  tDCA  = new TH2D("BefDCA", "DCA;[cm];Pt [GeV]", 50,0,1.001, 24,0,12); tQACuts->Add(tDCA);
  tCTP  = new TH2D("BefCTP", "CTP;;Pt [GeV]",     80,0.91,1.001, 24,0,12); tQACuts->Add(tCTP);
  tD0   = new TH2D("BefD0",  "D0;[cm];Pt [GeV]",  50,-1,+1, 24,0,12); tQACuts->Add(tD0);
  tD0D0 = new TH2D("BefD0D0","D0D0;[cm^{2}];Pt [GeV]", 100,-1,+1, 24,0,12); tQACuts->Add(tD0D0);
  tAP   = new TH3D("BefAP",  "AP;#alpha;q_{t}[GeV];Pt [GeV]", 80,-1,+1, 90,0,0.3, 24,0,12); tQACuts->Add(tAP);
  fQAList->Add(tQACuts);

  tQACuts = new TList(); tQACuts->SetOwner(); tQACuts->SetName("QACutsAfter");
  tAP   = new TH3D("AftAP","AP;#alpha;q_{t}[GeV];Pt [GeV]", 80,-1,+1, 90,0,0.3, 24,0,12); tQACuts->Add(tAP);
  TH1D *tPOI = new TH1D("POI","POIs;multiplicity",100,0,1000); tQACuts->Add(tPOI);
  fQAList->Add(tQACuts);

}

/*
void  AliAnalysisTaskFlowStrange::AddQAPID() {
 TList *tQAPID;
 tQAPID=new TList();
 tQAPID->SetOwner();
 tQAPID->SetName("QAPID");
 fdEdX_all_K0=new TH2D("fdEdX_all_K0","fdEdX_all;p_{T} [GeV];dEdx",150,0,15,200,0,20);
 fdEdX_all_L=new TH2D("fdEdX_all_L","fdEdX_all;p_{T} [GeV];dEdx",150,0,15,200,0,20);
 fdEdX_cut_K0=new TH2D("fdEdX_cut_K0","fdEdX_cut;p_{T} [GeV];dEdx",150,0,15,200,0,20);
 fdEdX_cut_L=new TH2D("fdEdX_cut_L","fdEdX_cut;p_{T} [GeV];dEdx",150,0,15,200,0,20);
 tQATracks->Add(fdEdX_all_K0);
 tQATracks->Add(fdEdX_all_L);
 tQATracks->Add(fdEdX_cut_K0);
 tQATracks->Add(fdEdX_cut_L);
 fQAList->Add(tQAPID);
}
*/
//=======================================================================
void AliAnalysisTaskFlowStrange::UserExec(Option_t *)
{
  // user exec
  AliESDEvent *tESD=dynamic_cast<AliESDEvent*>(InputEvent());
  AliAODEvent *tAOD=dynamic_cast<AliAODEvent*>(InputEvent());
  Bool_t acceptEvent=kFALSE;
  fCandidates->SetLast(-1);
  if(tESD) {
    ((TH1D*)((TList*)fQAList->FindObject("Event"))->FindObject("Events"))->Fill(0);
    if(fCutsEvent->IsSelected(tESD)) {
      acceptEvent=kTRUE;
      ReadFromESDv0(tESD);
    }
  } else if(tAOD) {
    ((TH1D*)((TList*)fQAList->FindObject("Event"))->FindObject("Events"))->Fill(0);
    if(fCutsEvent->IsSelected(tAOD)) {
      acceptEvent=kTRUE;
      ReadFromAODv0(tAOD);
    }
  }
  if(!acceptEvent) return;
  ((TH1D*)((TList*)fQAList->FindObject("Event"))->FindObject("Events"))->Fill(1);

  ((TH1D*)((TList*)fQAList->FindObject("Event"))->FindObject("RFPTPC"))->Fill( fFlowEventTPC->GetNumberOfRPs() );
  Double_t mult=0;
  for(Int_t i=0;i!=fFlowEventVZE->GetNumberOfRPs();++i) {
    AliFlowTrackSimple *pTrack = fFlowEventVZE->GetTrack(i);
    mult += pTrack->Weight();
  }
  ((TH1D*)((TList*)fQAList->FindObject("Event"))->FindObject("RFPVZE"))->Fill( mult );

  if(fDebug) printf("TPCevent %d | VZEevent %d\n",
		    fFlowEventTPC->NumberOfTracks(),
		    fFlowEventVZE->NumberOfTracks() );
  ((TH1D*)((TList*)fQAList->FindObject("QACutsAfter"))->FindObject("POI"))->Fill( fCandidates->GetEntriesFast() );
  AddCandidates();

  PostData(1,fFlowEventTPC);
  PostData(2,fFlowEventVZE);
  PostData(3,fQAList);
  return;
}
//=======================================================================
void AliAnalysisTaskFlowStrange::AddCandidates()
{
  // adds candidates to flow events (untaging if necessary)
  if(fDebug) printf("I received %d candidates\n",fCandidates->GetEntriesFast());
  for(int iCand=0; iCand!=fCandidates->GetEntriesFast(); ++iCand ) {
    AliFlowCandidateTrack *cand = dynamic_cast<AliFlowCandidateTrack*>(fCandidates->At(iCand));
    if(!cand) continue;
    if(fDebug) printf(" >Checking at candidate %d with %d daughters: mass %f\n",
		      iCand,cand->GetNDaughters(),cand->Mass());
    // untagging ===>
    for(int iDau=0; iDau!=cand->GetNDaughters(); ++iDau) {
      if(fDebug) printf("  >Daughter %d with fID %d", iDau, cand->GetIDDaughter(iDau));
      for(int iRPs=0; iRPs!=fFlowEventTPC->NumberOfTracks(); ++iRPs ) {
        AliFlowTrack *iRP = dynamic_cast<AliFlowTrack*>(fFlowEventTPC->GetTrack( iRPs ));
        if (!iRP) continue;
        if( !iRP->InRPSelection() ) continue;
        if( cand->GetIDDaughter(iDau) == iRP->GetID() ) {
          if(fDebug) printf(" was in RP set");
          iRP->SetForRPSelection(kFALSE);
          fFlowEventTPC->SetNumberOfRPs( fFlowEventTPC->GetNumberOfRPs() -1 );
        }
      }
      if(fDebug) printf("\n");
    }
    // <=== untagging
    cand->SetForPOISelection(kTRUE);
    fFlowEventTPC->InsertTrack( ((AliFlowTrack*) cand) );
    fFlowEventVZE->InsertTrack( ((AliFlowTrack*) cand) );
  }
  if(fDebug) printf("TPCevent %d | VZEevent %d\n",
		    fFlowEventTPC->NumberOfTracks(),
		    fFlowEventVZE->NumberOfTracks() );
}
//=======================================================================
void AliAnalysisTaskFlowStrange::ReadFromESDv0(AliESDEvent *tESD)
{
  fCutsRFPTPC->SetEvent(tESD,MCEvent());
  fCutsRFPVZE->SetEvent(tESD,MCEvent());
  fCutsPOI->SetEvent(tESD,MCEvent());
  fFlowEventTPC->Fill(fCutsRFPTPC,fCutsPOI);
  fFlowEventVZE->Fill(fCutsRFPVZE,fCutsPOI);

  Int_t nV0s = tESD->GetNumberOfV0s();
  AliESDv0 *myV0;
  Double_t dQT, dALPHA, dPT, dMASS=0.0;
  for (Int_t i=0; i!=nV0s; ++i) {
    myV0 = (AliESDv0*) tESD->GetV0(i);
    if(!myV0) continue;
    if(myV0->Pt()<0.1) continue; // skipping low momentum
    Int_t pass = PassesESDCuts(myV0,tESD);
    if(pass==0) continue;
    if(fSpecie==0) {
      dMASS = myV0->GetEffMass(2,2);
    } else {
      if(pass==1) dMASS = myV0->GetEffMass(2,4);
      if(pass==2) dMASS = myV0->GetEffMass(4,2);
      if(pass>2)  dMASS = myV0->GetEffMass(2,4);
    }
    dPT=myV0->Pt();
    dQT=myV0->PtArmV0();
    dALPHA=myV0->AlphaV0();
    ((TH3D*)((TList*)fQAList->FindObject("QACutsAfter"))->FindObject("AftAP"))  ->Fill(dALPHA,dQT,dPT);
    MakeTrack(dMASS, dPT, myV0->Phi(), myV0->Eta(), 
	      ((AliESDtrack*)tESD->GetTrack(myV0->GetIndex(0)))->GetID(),
	      ((AliESDtrack*)tESD->GetTrack(myV0->GetIndex(1)))->GetID());
  }
  return;
}
//=======================================================================
Int_t AliAnalysisTaskFlowStrange::PassesESDCuts(AliESDv0 *myV0, AliESDEvent *tESD)
{
  TVector3 vv = TVector3(tESD->GetPrimaryVertex()->GetX(),
			 tESD->GetPrimaryVertex()->GetY(),
			 tESD->GetPrimaryVertex()->GetZ()); // primary vertex
  if(myV0->GetOnFlyStatus() ) return 0;

  AliESDtrack *iT, *jT;

  // TESTING CHARGE
  int iPos, iNeg;
  iT=(AliESDtrack*) tESD->GetTrack( myV0->GetIndex(0) );
  if(iT->Charge()>0) {
    iPos = 0; iNeg = 1;
  } else {
    iPos = 1; iNeg = 0;
  }
  // END OF TEST

  iT=(AliESDtrack*) tESD->GetTrack( myV0->GetIndex(iPos) ); //positive daughter
  if(!fCutsDau->IsSelected(iT) ) return 0;
  jT=(AliESDtrack*) tESD->GetTrack( myV0->GetIndex(iNeg) ); //negative daughter
  if(!fCutsDau->IsSelected(jT) ) return 0;
  TVector3 vp,vl;
  vp=TVector3( myV0->Xv(), myV0->Yv(), myV0->Zv() ); // v0 vertex (same as GetXYZ)
  vl=vp-vv; // flight line
  Double_t dDL=vl.Mag();
  Double_t dDCA=myV0->GetDcaV0Daughters();
  Double_t dCTP=myV0->GetV0CosineOfPointingAngle();
  Double_t dD0P=iT->GetD(vv.X(),vv.Y(),tESD->GetMagneticField());
  Double_t dD0M=jT->GetD(vv.X(),vv.Y(),tESD->GetMagneticField());
  Double_t dD0D0=dD0P*dD0M;
  Double_t dQT=myV0->PtArmV0();
  Double_t dALPHA=myV0->AlphaV0();
  Double_t dPT=myV0->Pt();
  Double_t dETA=myV0->Eta();
  ((TH2D*)((TList*)fQAList->FindObject("QACutsBefore"))->FindObject("BefDL"))  ->Fill(dDL,dPT);
  ((TH2D*)((TList*)fQAList->FindObject("QACutsBefore"))->FindObject("BefDCA")) ->Fill(dDCA,dPT);
  ((TH2D*)((TList*)fQAList->FindObject("QACutsBefore"))->FindObject("BefCTP")) ->Fill(dCTP,dPT);
  ((TH2D*)((TList*)fQAList->FindObject("QACutsBefore"))->FindObject("BefD0"))  ->Fill(dD0M,dPT);
  ((TH2D*)((TList*)fQAList->FindObject("QACutsBefore"))->FindObject("BefD0"))  ->Fill(dD0P,dPT);
  ((TH2D*)((TList*)fQAList->FindObject("QACutsBefore"))->FindObject("BefD0D0"))->Fill(dD0D0,dPT);
  ((TH3D*)((TList*)fQAList->FindObject("QACutsBefore"))->FindObject("BefAP"))  ->Fill(dALPHA,dQT,dPT);
  Int_t passes = 1;
  if(dDL  <fV0Cuts[0]) passes = 0;
  if(dDCA >fV0Cuts[1]) passes = 0;
  if(dCTP <fV0Cuts[2]) passes = 0;
  if(TMath::Abs(dD0P) <fV0Cuts[3]) passes = 0;
  if(TMath::Abs(dD0M) <fV0Cuts[3]) passes = 0;
  if(dD0D0>fV0Cuts[4]) passes = 0;
  if(dETA <fV0Cuts[6]) passes = 0;
  if(dETA >fV0Cuts[7]) passes = 0;
  if(fSpecie==0) if(dQT<fV0Cuts[5]) passes = 0;

  if(passes&&fV0Cuts[8]) {
    switch(fSpecie) {
    case 0: // K0 PID
      if( (jT->GetInnerParam()->GetP()<15) &&
	  (TMath::Abs(fPIDResponse->NumberOfSigmasTPC(jT,AliPID::kPion))>3.) )
	passes = 0;
      if( (iT->GetInnerParam()->GetP()<15) &&
	  (TMath::Abs(fPIDResponse->NumberOfSigmasTPC(iT,AliPID::kPion))>3.) )
	passes = 0;
      break;
    case 1: // Lambda PID i==pos j==neg
      Int_t antilambda=2;
      Int_t lambda=1;
      // antilambda p- pi+
      if( (iT->GetInnerParam()->GetP()<15) && 
	  (TMath::Abs(fPIDResponse->NumberOfSigmasTPC(iT,AliPID::kPion))>3.) )
	antilambda = 0;
      if( (jT->GetInnerParam()->GetP()<15) && 
	  (TMath::Abs(fPIDResponse->NumberOfSigmasTPC(jT,AliPID::kProton))>3.) )
	antilambda = 0;
      // lambda p+ pi-
      if( (iT->GetInnerParam()->GetP()<15) && 
	  (TMath::Abs(fPIDResponse->NumberOfSigmasTPC(iT,AliPID::kProton))>3.) )
	lambda = 0;
      if( (jT->GetInnerParam()->GetP()<15) && 
	  (TMath::Abs(fPIDResponse->NumberOfSigmasTPC(jT,AliPID::kPion))>3.) )
	lambda = 0;
      passes=lambda+antilambda;
      break;
    }
  }
  return passes;
}
//=======================================================================
void AliAnalysisTaskFlowStrange::ReadFromAODv0(AliAODEvent *tAOD)
{
  fCutsRFPTPC->SetEvent(tAOD,MCEvent());
  fCutsRFPVZE->SetEvent(tAOD,MCEvent());
  fCutsPOI->SetEvent(tAOD,MCEvent());
  fFlowEventTPC->Fill(fCutsRFPTPC,fCutsPOI);
  fFlowEventVZE->Fill(fCutsRFPVZE,fCutsPOI);
  Int_t nV0s = tAOD->GetNumberOfV0s();
  AliAODv0 *myV0;
  Double_t dQT, dALPHA, dPT, dMASS=0.0;
  for (Int_t i=0; i!=nV0s; ++i) {
    myV0 = (AliAODv0*) tAOD->GetV0(i);
    if(!myV0) continue;
    if(myV0->Pt()<0.1) continue; // skipping low momentum
    Int_t pass = PassesAODCuts(myV0,tAOD);
    if(pass==0) continue;
    if(fSpecie==0) {
      dMASS = myV0->MassK0Short();
    } else {
      if(pass==1) dMASS = myV0->MassLambda();
      if(pass==2) dMASS = myV0->MassAntiLambda();
      if(pass>2)  dMASS = myV0->MassLambda();
    }
    dPT=myV0->Pt();
    dQT=myV0->PtArmV0();
    dALPHA=myV0->AlphaV0();
    ((TH3D*)((TList*)fQAList->FindObject("QACutsAfter"))->FindObject("AftAP"))  ->Fill(dALPHA,dQT,dPT);
    MakeTrack(dMASS, dPT, myV0->Phi(), myV0->Eta(),
              ((AliAODTrack*) myV0->GetDaughter(0))->GetID(),
              ((AliAODTrack*) myV0->GetDaughter(1))->GetID());
  }
  return;
}
//=======================================================================
Int_t AliAnalysisTaskFlowStrange::PassesAODCuts(AliAODv0 *myV0, AliAODEvent *tAOD)
{
  if (myV0->GetOnFlyStatus() ) return 0;

  //the following is needed in order to evualuate track-quality
  AliAODTrack *iT, *jT;
  AliAODVertex *vV0s = myV0->GetSecondaryVtx();
  Double_t pos[3],cov[6];
  vV0s->GetXYZ(pos);
  vV0s->GetCovarianceMatrix(cov);
  const AliESDVertex vESD(pos,cov,100.,100);


  // TESTING CHARGE
  int iPos, iNeg;
  iT=(AliAODTrack*) myV0->GetDaughter(0);
  if(iT->Charge()>0) {
    iPos = 0; iNeg = 1;
  } else {
    iPos = 1; iNeg = 0;
  }
  // END OF TEST

  iT=(AliAODTrack*) myV0->GetDaughter(iPos); // positive
  AliESDtrack ieT( iT );
  ieT.SetTPCClusterMap( iT->GetTPCClusterMap() );
  ieT.SetTPCSharedMap( iT->GetTPCSharedMap() );
  ieT.SetTPCPointsF( iT->GetTPCNclsF() );
  ieT.RelateToVertex(&vESD, tAOD->GetMagneticField(), 100);
  if (!fCutsDau->IsSelected( &ieT ) ) return 0;

  jT=(AliAODTrack*) myV0->GetDaughter(iNeg); // negative
  AliESDtrack jeT( jT );
  jeT.SetTPCClusterMap( jT->GetTPCClusterMap() );
  jeT.SetTPCSharedMap( jT->GetTPCSharedMap() );
  jeT.SetTPCPointsF( jT->GetTPCNclsF() );
  jeT.RelateToVertex(&vESD, tAOD->GetMagneticField(), 100);
  if (!fCutsDau->IsSelected( &jeT ) ) return 0;

  Double_t pvertex[3];
  pvertex[0]=tAOD->GetPrimaryVertex()->GetX();
  pvertex[1]=tAOD->GetPrimaryVertex()->GetY();
  pvertex[2]=tAOD->GetPrimaryVertex()->GetZ();
  Double_t dDL=myV0->DecayLengthV0( pvertex );
  Double_t dDCA=myV0->DcaV0Daughters();
  Double_t dCTP=myV0->CosPointingAngle( pvertex );
  Double_t dD0P=ieT.GetD(pvertex[0],pvertex[1],tAOD->GetMagneticField());
  Double_t dD0M=jeT.GetD(pvertex[0],pvertex[1],tAOD->GetMagneticField());
  Double_t dD0D0=dD0P*dD0M;
  //Double_t dRAD = myV0->RadiusV0;
  //Double_t dOpenAngle = myV0->OpenAngleV0();
  //Double_t dRapK=myV0->RapK0Short();
  //Double_t dRapL=myV0->RapLambda();
  //Double_t dDCAv0 = myV0->DcaV0ToPrimVertex();
  Double_t dQT=myV0->PtArmV0();
  Double_t dALPHA=myV0->AlphaV0();
  Double_t dPT=myV0->Pt();
  Double_t dETA=myV0->Eta();
  ((TH2D*)((TList*)fQAList->FindObject("QACutsBefore"))->FindObject("BefDL"))  ->Fill(dDL,dPT);
  ((TH2D*)((TList*)fQAList->FindObject("QACutsBefore"))->FindObject("BefDCA")) ->Fill(dDCA,dPT);
  ((TH2D*)((TList*)fQAList->FindObject("QACutsBefore"))->FindObject("BefCTP")) ->Fill(dCTP,dPT);
  ((TH2D*)((TList*)fQAList->FindObject("QACutsBefore"))->FindObject("BefD0"))  ->Fill(dD0M,dPT);
  ((TH2D*)((TList*)fQAList->FindObject("QACutsBefore"))->FindObject("BefD0"))  ->Fill(dD0P,dPT);
  ((TH2D*)((TList*)fQAList->FindObject("QACutsBefore"))->FindObject("BefD0D0"))->Fill(dD0D0,dPT);
  ((TH3D*)((TList*)fQAList->FindObject("QACutsBefore"))->FindObject("BefAP"))  ->Fill(dALPHA,dQT,dPT);
  Int_t passes = 1;
  if(dDL  <fV0Cuts[0]) passes = 0;
  if(dDCA >fV0Cuts[1]) passes = 0;
  if(dCTP <fV0Cuts[2]) passes = 0;
  if(TMath::Abs(dD0P) <fV0Cuts[3]) passes = 0;
  if(TMath::Abs(dD0M) <fV0Cuts[3]) passes = 0;
  if(dD0D0>fV0Cuts[4]) passes = 0;
  if(dETA <fV0Cuts[6]) passes = 0;
  if(dETA >fV0Cuts[7]) passes = 0;
  if(fSpecie==0) if(dQT<fV0Cuts[5]) passes = 0;

  if(passes&&fV0Cuts[8]) {
    switch(fSpecie) {
    case 0: // K0 PID
      if( (jT->GetTPCmomentum()<15) &&
	  (TMath::Abs(fPIDResponse->NumberOfSigmasTPC(jT,AliPID::kPion))>3.) )
	passes = 0;
      if( (iT->GetTPCmomentum()<15) &&
	  (TMath::Abs(fPIDResponse->NumberOfSigmasTPC(iT,AliPID::kPion))>3.) )
	passes = 0;
      break;
    case 1: // Lambda PID  i==pos j ==neg
      Int_t antilambda=2;
      Int_t lambda=1;
      // antilambda p- pi+
      if( (iT->GetTPCmomentum()<15) &&
	  (TMath::Abs(fPIDResponse->NumberOfSigmasTPC(iT,AliPID::kPion))>3.) )
	antilambda = 0;
      if( (jT->GetTPCmomentum()<15) &&
	  (TMath::Abs(fPIDResponse->NumberOfSigmasTPC(jT,AliPID::kProton))>3.) )
	antilambda = 0;
      // lambda p+ pi-
      if( (iT->GetTPCmomentum()<15) &&
	  (TMath::Abs(fPIDResponse->NumberOfSigmasTPC(iT,AliPID::kProton))>3.) )
	lambda = 0;
      if( (jT->GetTPCmomentum()<15) &&
	  (TMath::Abs(fPIDResponse->NumberOfSigmasTPC(jT,AliPID::kPion))>3.) )
	lambda = 0;
      passes=lambda+antilambda;
      break;
    }
  }
  return passes;
}
//=======================================================================
void AliAnalysisTaskFlowStrange::Terminate(Option_t *)
{
  //terminate
}
//=======================================================================
void AliAnalysisTaskFlowStrange::MakeTrack( Double_t mass, Double_t pt, Double_t phi,
					    Double_t eta, Int_t iid, Int_t jid )
{
  // create track for flow tasks
  if(fCandidates->GetLast()+1>=fCandidates->GetSize()) {
    fCandidates->Expand( 2*fCandidates->GetSize() );
  }
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
  oTrack->AddDaughter(iid);
  oTrack->AddDaughter(jid);
  oTrack->SetForPOISelection(kTRUE);
  oTrack->SetForRPSelection(kFALSE);
  if(overwrite) {
    fCandidates->SetLast( fCandidates->GetLast()+1 );
  } else {
    fCandidates->AddLast(oTrack);
  }
  return;
}
//=======================================================================
void AliAnalysisTaskFlowStrange::SetCommonConstants(Int_t massBins, Double_t minMass, Double_t maxMass)
{
  // setter for mass bins
  fMassBins = massBins;
  fMinMass = minMass;
  fMaxMass = maxMass;
}
//=======================================================================
void AliAnalysisTaskFlowStrange::SetCuts2010(int set) {
  // defines cuts to be used
  // fV0Cuts[9] dl dca ctp d0 d0d0 qt minEta maxEta PID
  switch(set) {
  case(0): // No cuts
    fV0Cuts[0] = -1e+6; fV0Cuts[1] = +1e+6; fV0Cuts[2] = -1e+6;
    fV0Cuts[3] = -1e+6; fV0Cuts[4] = +1e+6; fV0Cuts[5] = -1e+6;
    fV0Cuts[6] = -1e+6; fV0Cuts[7] = +1e+6; fV0Cuts[8] = 0;
  case(1): // Tight cuts
    fV0Cuts[0] = +0.5; fV0Cuts[1] = +0.5; fV0Cuts[2] = +0.998;
    fV0Cuts[3] = +0.1; fV0Cuts[4] = +0.0; fV0Cuts[5] = +0.105;
    fV0Cuts[6] = -0.8; fV0Cuts[7] = +0.8; fV0Cuts[8] = 0;
    break;
  case(2): // Tight cuts + PID
    fV0Cuts[0] = +0.5; fV0Cuts[1] = +0.5; fV0Cuts[2] = +0.998;
    fV0Cuts[3] = +0.1; fV0Cuts[4] = +0.0; fV0Cuts[5] = +0.105;
    fV0Cuts[6] = -0.8; fV0Cuts[7] = +0.8; fV0Cuts[8] = 1;
    break;
  }
}
