
//_____________________________________________________________________________
//    Class for semiforward UPC filter
//    Author: Jaroslav Adam
//
//    Fill structure of AliUPCEvent
//_____________________________________________________________________________

// c++ headers
#include <iostream>
#include <string.h>

// root headers
#include "TList.h"
#include "TH1I.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TParticle.h"
#include "TObjString.h"
#include "TFile.h"
#include "TArrayF.h"

// aliroot headers
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliMCEvent.h"
#include "AliAODVZERO.h"
#include "AliAODZDC.h"
#include "AliAODVertex.h"
#include "AliESDVZERO.h"
#include "AliESDZDC.h"
#include "AliESDVertex.h"
#include "AliAODTrack.h"
#include "AliMultiplicity.h"
#include "AliESDtrack.h"
#include "AliESDMuonTrack.h"
#include "AliAODMCParticle.h"
#include "AliMCParticle.h"
#include "AliTriggerAnalysis.h"
#include "AliMuonTrackCuts.h"
#include "AliESDtrackCuts.h"
#include "AliAODPid.h"
#include "AliGenEventHeader.h"
#include "AliAODMCHeader.h"
#include "AliPIDResponse.h"

// my headers
#include "AliAnalysisTaskUpcFilterSemiforward.h"
#include "AliUPCTrack.h"
#include "AliUPCMuonTrack.h"
#include "AliUPCEvent.h"

using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskUpcFilterSemiforward);

// task for upc semiforward filter
// jaroslav.adam@cern.ch

//_____________________________________________________________________________
AliAnalysisTaskUpcFilterSemiforward::AliAnalysisTaskUpcFilterSemiforward(const char *name)
 :AliAnalysisTaskSE(name),
  fIsESD(0), fIsMC(0), fMuonCuts(0x0), fTriggerAna(0x0), fCutsList(0x0), fPIDResponse(0x0),
  fHistList(0x0), fCounter(0x0), fTriggerCounter(0x0),
  fUPCEvent(0x0), fUPCTree(0x0)
{

  // Constructor

  DefineOutput(1, TTree::Class());
  DefineOutput(2, TList::Class());

}//AliAnalysisTaskUpcFilterSemiforward

//_____________________________________________________________________________
AliAnalysisTaskUpcFilterSemiforward::~AliAnalysisTaskUpcFilterSemiforward()
{
  // destructor

  if(fHistList) {delete fHistList; fHistList = 0x0;}
  if(fCounter) {delete fCounter; fCounter = 0x0;}
  if(fTriggerCounter) {delete fTriggerCounter; fTriggerCounter = 0x0;}
  if(fUPCEvent) {delete fUPCEvent; fUPCEvent = 0x0;}
  if(fUPCTree) {delete fUPCTree; fUPCTree = 0x0;}
  if(fMuonCuts) {delete fMuonCuts; fMuonCuts = 0x0;}
  if(fTriggerAna) {delete fTriggerAna; fTriggerAna = 0x0;}
  if(fCutsList) {delete[] fCutsList; fCutsList=0x0;}
  if(fPIDResponse) {delete fPIDResponse; fPIDResponse=0x0;}

}//~AliAnalysisTaskUpcFilterSemiforward

//_____________________________________________________________________________
void AliAnalysisTaskUpcFilterSemiforward::UserCreateOutputObjects()
{
  //muon track cuts
  fMuonCuts = new AliMuonTrackCuts("StdMuonCuts", "StdMuonCuts");
  fMuonCuts->SetFilterMask ( AliMuonTrackCuts::kMuPdca );
  fMuonCuts->Print("mask");
  fMuonCuts->SetAllowDefaultParams(kTRUE);
  fMuonCuts->SetPassName("muon_pass2");

  //trigger analysis for SPD FO fired chips in MC
  fTriggerAna = new AliTriggerAnalysis();
  fTriggerAna->SetAnalyzeMC( fIsMC );

  //PID response, input handler from AliAnalysisTaskSE
  //fPIDResponse = fInputHandler->GetPIDResponse();
  AliInputEventHandler* inputHandler = (AliInputEventHandler*) ( AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler() );
  fPIDResponse = inputHandler->GetPIDResponse();

  //ESD track cuts
  if( fIsESD ) {
    fCutsList = new AliESDtrackCuts*[32];
    for(UInt_t i=0; i<32; i++) fCutsList[i] = 0x0;

    // bits 0 - 10 set in $ALICE_ROOT/ANALYSIS/ESDfilter/macros/AddTaskESDFilter.C

    // Cuts on primary tracks
    AliESDtrackCuts* esdTrackCutsL = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts();
    fCutsList[0] = esdTrackCutsL;

    // standard cuts with very loose DCA
    AliESDtrackCuts* esdTrackCutsH = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kFALSE); 
    esdTrackCutsH->SetMaxDCAToVertexXY(2.4);
    esdTrackCutsH->SetMaxDCAToVertexZ(3.2);
    esdTrackCutsH->SetDCAToVertex2D(kTRUE);
    fCutsList[4] = esdTrackCutsH;

    // standard cuts with tight DCA cut
    AliESDtrackCuts *esdTrackCutsH2 = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011();
    fCutsList[5] = esdTrackCutsH2;

    // bits 11 - 31 free

    // selection by Daniele, UPC meeting, 15 July 2014
    // https://indico.cern.ch/event/330483/contribution/3/material/slides/0.pdf

    AliESDtrackCuts *cuts11 = new AliESDtrackCuts();
    cuts11->SetMinNClustersTPC(70);
    cuts11->SetRequireTPCRefit(kTRUE);
    cuts11->SetRequireITSRefit(kTRUE);
    cuts11->SetMaxChi2PerClusterTPC(4);
    cuts11->SetAcceptKinkDaughters(kFALSE);
    cuts11->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
    cuts11->SetMaxDCAToVertexXYPtDep("0.0105+0.0350/pt^1.01");
    cuts11->SetMaxDCAToVertexZ(2);
    cuts11->SetDCAToVertex2D(kFALSE);
    cuts11->SetRequireSigmaToVertex(kFALSE);
    cuts11->SetMaxChi2PerClusterITS(5.);
    fCutsList[11] = cuts11;

    AliESDtrackCuts *cuts12 = new AliESDtrackCuts();
    cuts12->SetMinNClustersTPC(70);
    cuts12->SetRequireTPCRefit(kTRUE);
    cuts12->SetRequireITSRefit(kTRUE);
    cuts12->SetMaxChi2PerClusterTPC(4);
    cuts12->SetAcceptKinkDaughters(kFALSE);
    cuts12->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
    cuts12->SetMaxDCAToVertexXYPtDep("0.0105+0.0350/pt^1.01");
    cuts12->SetMaxDCAToVertexZ(2);
    cuts12->SetDCAToVertex2D(kFALSE);
    cuts12->SetRequireSigmaToVertex(kFALSE);
    fCutsList[12] = cuts12;

    // setting in Evgeny's trees
    AliESDtrackCuts *cuts13 = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011();
    cuts13->SetMaxChi2TPCConstrainedGlobal(1e+10);
    fCutsList[13] = cuts13;

  }

  //output histograms
  fHistList = new TList();
  fHistList->SetOwner();
  fCounter = new TH1I("fCounter", "fCounter", 100, 1, 101);
  fHistList->Add(fCounter);
  fTriggerCounter = new TH2I("fTriggerCounter", "fTriggerCounter", 44000, 154000, 198000, AliUPCEvent::fgkNtrg+1, 0, AliUPCEvent::fgkNtrg+1);
  fHistList->Add(fTriggerCounter);

  //errors in counter start at 51
  if( !fPIDResponse ) fCounter->Fill( 51 ); // no PID response

  //output tree
  fUPCEvent  = new AliUPCEvent();
  fUPCEvent->SetIsMC( fIsMC );
  fUPCEvent->SetIsESD( fIsESD );

  TDirectory *pwd = gDirectory;
  OpenFile(1);
  fUPCTree = new TTree("fUPCTree", "fUPCTree");
  pwd->cd();
  fUPCTree->Branch("fUPCEvent", &fUPCEvent);

  PostData(1, fUPCTree);
  PostData(2, fHistList);

}//UserCreateOutputObjects

//_____________________________________________________________________________
void AliAnalysisTaskUpcFilterSemiforward::NotifyRun()
{

  //cout<<"AliAnalysisTaskUpcFilterSemiforward::NotifyRun()"<<endl;

  fMuonCuts->SetRun(fInputHandler); // use input handler from AliAnalysisTaskSE

  if( fIsESD ) { ((AliESDEvent*) InputEvent())->InitMagneticField(); }

}//NotifyRun

//_____________________________________________________________________________
void AliAnalysisTaskUpcFilterSemiforward::UserExec(Option_t *) 
{

  //cout<<"#################### Next event ##################"<<endl;

  // input event
  AliVEvent *vEvent = InputEvent();
  if(!vEvent) return;

  fUPCEvent->ClearEvent();

  fCounter->Fill( 1 ); // 1 = analyzed events

  // trigger
  TString trigger = vEvent->GetFiredTriggerClasses();
  Bool_t trgClasses[AliUPCEvent::fgkNtrg]; // array of fired trigger classes
  for(Int_t itrg=0; itrg<AliUPCEvent::fgkNtrg; itrg++) trgClasses[itrg] = kFALSE;

  trgClasses[1] = trigger.Contains("CMUP6-B"); // p-Pb FW
  trgClasses[2] = trigger.Contains("CMUP3-B"); // Pb-p FW
  trgClasses[3] = trigger.Contains("CMUP8-B"); // Pb-p FW

  trgClasses[4] = trigger.Contains("CMUP7-B"); // p-Pb SFW
  trgClasses[5] = trigger.Contains("CMUP5-B"); // Pb-p SFW
  trgClasses[6] = trigger.Contains("CMUP9-B"); // Pb-p SFW

  trgClasses[7] = trigger.Contains("CCUP7-B"); // CEN

  trgClasses[8] = trigger.Contains("CMUP7-ACE");
  trgClasses[9] = trigger.Contains("CMUP5-ACE");
  trgClasses[10]= trigger.Contains("CMUP9-ACE");

  trgClasses[11]= trigger.Contains("CMUP1-B"); // PbPb FW

  Bool_t isTrg = kFALSE;
  for(Int_t itrg=1; itrg<AliUPCEvent::fgkNtrg; itrg++) {
    if(!trgClasses[itrg]) continue;
    //trigger at itrg is fired
    fUPCEvent->SetTriggerClass( itrg , kTRUE );
    fTriggerCounter->Fill( vEvent->GetRunNumber() , itrg );
    isTrg = kTRUE;
  }
  if(!isTrg && !fIsMC) {PostData(2, fHistList); return;}
  //event passed the trigger

  fCounter->Fill( 2 ); // 2 = events after trigger (ESD and AOD)

  // ESD / AOD specific tasks: MC, SPD FO, L0 inputs, SPD tracklets, tracks, ZDC tdc
  Bool_t stat;
  if( fIsESD ) stat = RunESD();
  if(!fIsESD ) stat = RunAOD();
  if( !stat && !fIsMC ) {PostData(2, fHistList); return;}
  //event passed ESD / AOD specific selection

  fCounter->Fill( 3 ); // 3 = events after ESD / AOD specific part

  // input data
  const char *filnam = ((TTree*) GetInputData(0))->GetCurrentFile()->GetName();
  // reconstruction pass: -1 = unknown, 1 = pass1, 2 = pass2, counter indices: 9 (unknown), 11 (pass1), 12 (pass2)
  fUPCEvent->SetRecoPass( -1 );
  if( strstr(filnam,"/pass1/") ) {fUPCEvent->SetRecoPass( 1 ); fCounter->Fill( 11 );}
  if( strstr(filnam,"/pass2/") ) {fUPCEvent->SetRecoPass( 2 ); fCounter->Fill( 12 );}
  if( fUPCEvent->GetRecoPass() < 0 ) fCounter->Fill( 9 );
  fUPCEvent->SetInputFileName( filnam );
  fUPCEvent->SetEventNumber( ((TTree*) GetInputData(0))->GetTree()->GetReadEntry() );
  fUPCEvent->SetRunNumber( vEvent->GetRunNumber() );

  //VZERO
  AliVVZERO *dataVZERO = vEvent->GetVZEROData();
  if(!dataVZERO) {PostData(2, fHistList); return;}

  fUPCEvent->SetV0ADecision( dataVZERO->GetV0ADecision() );
  fUPCEvent->SetV0CDecision( dataVZERO->GetV0CDecision() );
  for (UInt_t iv=0; iv<32; iv++) {
    if( dataVZERO->BBTriggerV0C((Int_t)iv) ) fUPCEvent->SetBBtriggerV0Cmask(iv);
    if( dataVZERO->GetBBFlag((Int_t)iv) ) fUPCEvent->SetBBFlagV0Cmask(iv);
  }

  //ZDC
  AliVZDC *dataZDC = vEvent->GetZDCData();
  if(!dataZDC) {PostData(2, fHistList); return;}

  //energy in ZDC
  Double_t eZnc=0., eZpc=0., eZna=0., eZpa=0.;
  for(Int_t i=0; i<5; i++) {
    eZnc += dataZDC->GetZNCTowerEnergy()[i];
    eZpc += dataZDC->GetZPCTowerEnergy()[i];
    eZna += dataZDC->GetZNATowerEnergy()[i];
    eZpa += dataZDC->GetZPATowerEnergy()[i];
  }
  fUPCEvent->SetZNCEnergy( eZnc );
  fUPCEvent->SetZPCEnergy( eZpc );
  fUPCEvent->SetZNAEnergy( eZna );
  fUPCEvent->SetZPAEnergy( eZpa );

  //default primary vertex
  const AliVVertex *vtx = vEvent->GetPrimaryVertex();
  if(!vtx) {PostData(2, fHistList); return;}
  Double_t pos[3]; vtx->GetXYZ(pos);
  Double_t cov[6]; vtx->GetCovarianceMatrix(cov);
  fUPCEvent->SetPrimaryVertex(pos, vtx->GetChi2perNDF(), cov, vtx->GetNContributors());
  const char *vtxtitle = vtx->GetTitle();
  fUPCEvent->SetPrimaryVertexTitle( vtxtitle );

  fCounter->Fill( 4 ); // 4 = events written to the tree (ESD and AOD)

  fUPCTree ->Fill();
  PostData(1, fUPCTree);
  PostData(2, fHistList);

}//UserExec

//_____________________________________________________________________________
Bool_t AliAnalysisTaskUpcFilterSemiforward::RunAOD()
{
  //cout<<"#################### AOD event ##################"<<endl;

  //input AOD event
  AliAODEvent *aodEvent = (AliAODEvent*) InputEvent();
  if(!aodEvent) return kFALSE;

  fCounter->Fill( 22 ); // 22 = AOD analyzed events

  if(fIsMC) RunAODMC( (TClonesArray*) aodEvent->GetList()->FindObject(AliAODMCParticle::StdBranchName()),
                      (AliAODMCHeader*) aodEvent->FindListObject(AliAODMCHeader::StdBranchName())
                    );

  //nominal interaction point
  static AliAODVertex *vtx0 = 0x0;
  if( !vtx0 ) {
    // make a vertex at coordinates 0, 0, 0
    vtx0 = new AliAODVertex();
    vtx0->SetX(0.); vtx0->SetY(0.); vtx0->SetZ(0.);
  }

  //array of tracks for DCAs calculation
  //static TClonesArray *trackArray = 0x0;
  //if( !trackArray ) trackArray = new TClonesArray("AliAODTrack", 3); // number of DCAs

  Double_t pxpypz[3];
  UChar_t maskMan;
  Double_t xyzDca[2], cov[3];
  Float_t b[2], covF[3];
  Int_t nmun=0, ncen=0;
  // AOD tracks loop
  for(Int_t itr=0; itr<aodEvent->GetNumberOfTracks(); itr++) {
    AliAODTrack *trk = dynamic_cast<AliAODTrack *>(aodEvent->GetTrack(itr));
    if( !trk ) continue;

    //muon track
    if( trk->IsMuonTrack() ) {

      AliUPCMuonTrack *upcMuon = fUPCEvent->AddMuonTrack();
      fCounter->Fill( 31 ); // 31 = added muon track
      upcMuon->SetPtEtaPhi( trk->Pt(), trk->Eta(), trk->Phi() );
      upcMuon->SetCharge( trk->Charge() );
      upcMuon->SetMatchTrigger( trk->GetMatchTrigger() );
      upcMuon->SetRAtAbsorberEnd( trk->GetRAtAbsorberEnd() );
      upcMuon->SetChi2perNDF( trk->Chi2perNDF() );
      upcMuon->SetDCA( trk->DCA() );
      upcMuon->SetPxDCA( fMuonCuts->IsSelected(trk) );

      nmun++;
    }
    //central track
    else {

      maskMan = 0;
      if( trk->HasPointOnITSLayer(0) )                 { maskMan |= 1 << 0; }
      if( trk->HasPointOnITSLayer(1) )                 { maskMan |= 1 << 1; }
      if( trk->GetStatus() & AliESDtrack::kITSrefit )  { maskMan |= 1 << 2; }
      if( trk->GetStatus() & AliESDtrack::kTPCrefit )  { maskMan |= 1 << 3; }
      if( trk->GetTOFsignal() < 99999. && trk->GetTOFsignal() > 0. )  { maskMan |= 1 << 4; }

      //selection for at least one point in ITS and TPC refit
      //if( !(maskMan & (1 << 0)) && !(maskMan & (1 << 1)) ) continue;
      //if( !(maskMan & (1 << 3)) ) continue;

      //selection for its refit and tpc refit
      //if( !(maskMan & (1 << 2)) || !(maskMan & (1 << 3)) ) continue;

      //selection for tpc refit
      //if( !(maskMan & (1 << 3)) ) continue;

      //selection for at least one point in SPD and ITS refit and TPC refit
      //if( !(maskMan & (1 << 0)) && !(maskMan & (1 << 1)) ) continue;
      //if( !(maskMan & (1 << 2)) || !(maskMan & (1 << 3)) ) continue;

      //selection for at least one point in SPD and at least one TPC cluster
      if( !(maskMan & (1 << 0)) && !(maskMan & (1 << 1)) ) continue;
      if( trk->GetTPCNcls() == 0 ) continue;

      trk->GetPxPyPz(pxpypz);
      AliAODPid *apid = trk->GetDetPid();
      if(!apid) continue;

      AliUPCTrack *upcTrack = fUPCEvent->AddTrack();
      fCounter->Fill( 32 ); // 32 = added central track
      upcTrack->SetPxPyPz( pxpypz );
      upcTrack->SetMaskMan( maskMan );
      upcTrack->SetFilterMap( trk->GetFilterMap() );
      upcTrack->SetCharge( trk->Charge() );
      upcTrack->SetChi2perNDF( trk->Chi2perNDF() );
      upcTrack->SetTPCmomentum( apid->GetTPCmomentum() );
      upcTrack->SetTPCsignal( apid->GetTPCsignal() );
      upcTrack->SetTPCNcls( trk->GetTPCNcls() );
      upcTrack->SetTPCCrossedRows( (Float_t) trk->GetTPCNCrossedRows() );
      upcTrack->SetTPCNclsF( trk->GetTPCNclsF() );
      upcTrack->SetTPCNclsS( trk->GetTPCnclsS() );
      upcTrack->SetITSClusterMap( trk->GetITSClusterMap() );
      upcTrack->SetTOFsignal( trk->GetTOFsignal() );

      //array for DCA
      //trackArray->Clear("C");
      //new((*trackArray)[0]) AliAODTrack(*trk);
      //new((*trackArray)[1]) AliAODTrack(*trk);
      //new((*trackArray)[2]) AliAODTrack(*trk);

      //DCA to default primary vertex
      AliAODTrack *track1 = (AliAODTrack*) trk->Clone("track1");
      //AliAODTrack *track1 = (AliAODTrack*) trackArray->At( 0 );
      track1->PropagateToDCA(aodEvent->GetPrimaryVertex(), aodEvent->GetMagneticField(), 9999., xyzDca, cov);
      for(Int_t i=0; i<2; i++) {b[i] = (Float_t) xyzDca[i]; covF[i] = (Float_t) cov[i];}
      upcTrack->SetImpactParameters(b, covF);
      delete track1;

      //DCA to SPD vertex
      track1 = (AliAODTrack*) trk->Clone("track1");
      //track1 = (AliAODTrack*) trackArray->At( 1 );
      track1->PropagateToDCA(aodEvent->GetPrimaryVertexSPD(), aodEvent->GetMagneticField(), 9999., xyzDca, cov);
      for(Int_t i=0; i<2; i++) {b[i] = (Float_t) xyzDca[i]; covF[i] = (Float_t) cov[i];}
      upcTrack->SetImpactParametersSPD(b, covF);
      delete track1;

      //DCA to nominal interaction point
      track1 = (AliAODTrack*) trk->Clone("track1");
      //track1 = (AliAODTrack*) trackArray->At( 2 );
      track1->PropagateToDCA(vtx0, aodEvent->GetMagneticField(), 9999., xyzDca, cov);
      for(Int_t i=0; i<2; i++) {b[i] = (Float_t) xyzDca[i]; covF[i] = (Float_t) cov[i];}
      upcTrack->SetImpactParametersIP(b, covF);
      delete track1;

      ncen++;
    }
  }// AOD tracks loop

  //selection for at least one muon or central track
  //if( nmun + ncen < 1 ) return kFALSE;

  //selection for at least one muon and at least one central track
  //if( nmun < 1 || ncen < 1 ) return kFALSE;


  // L0 trigger inputs
  fUPCEvent->SetL0inputs( aodEvent->GetHeader()->GetL0TriggerInputs() );

  // Tracklets
  fUPCEvent->SetNumberOfTracklets( aodEvent->GetTracklets()->GetNumberOfTracklets() );

  // AOD ZDC for TDC
  AliAODZDC *dataZDCAOD = aodEvent->GetZDCData();
  if(!dataZDCAOD) {PostData(2, fHistList); return kFALSE;}

  Bool_t znctdc = kFALSE, znatdc = kFALSE;
  //if( dataZDCAOD->GetZNCTime() != 0. ) znctdc = kTRUE;
  //if( dataZDCAOD->GetZNATime() != 0. ) znatdc = kTRUE;
  if( TMath::Abs(dataZDCAOD->GetZNCTime()) > 20. ) znctdc = kTRUE;
  if( TMath::Abs(dataZDCAOD->GetZNATime()) > 20. ) znatdc = kTRUE;
  fUPCEvent->SetZNCtdc( znctdc );
  fUPCEvent->SetZNAtdc( znatdc );
  fUPCEvent->SetZPCtdc( kFALSE );
  fUPCEvent->SetZPAtdc( kFALSE );
  fUPCEvent->SetZNCTime( dataZDCAOD->GetZNCTime() );
  fUPCEvent->SetZNATime( dataZDCAOD->GetZNATime() );

  //SPD primary vertex in AOD
  AliAODVertex *vtx = aodEvent->GetPrimaryVertexSPD();
  if(!vtx) {PostData(2, fHistList); return kFALSE;}
  Double_t posVtx[3]; vtx->GetXYZ( posVtx );
  Double_t covVtx[6]; vtx->GetCovarianceMatrix( covVtx );
  fUPCEvent->SetPrimaryVertexSPD(posVtx, vtx->GetChi2perNDF(), covVtx, vtx->GetNContributors());
  const char *vtxtitle = vtx->GetTitle();
  fUPCEvent->SetPrimaryVertexSPDtitle( vtxtitle );

  return kTRUE;

}//RunAOD

//_____________________________________________________________________________
void AliAnalysisTaskUpcFilterSemiforward::RunAODMC(TClonesArray *arrayMC, AliAODMCHeader *headerMC)
{
  // run over AOD mc particles

  if( !arrayMC || !headerMC ) return;

  //generated vertex in AOD MC
  Double_t vtxD[3]; headerMC->GetVertex( vtxD );
  Float_t vtxF[3]; for(Int_t i=0; i<3; i++) vtxF[i] = (Float_t) vtxD[i];
  fUPCEvent->SetPrimaryVertexMC( vtxF );

  //loop over mc particles
  for(Int_t imc=0; imc<arrayMC->GetEntriesFast(); imc++) {
    AliAODMCParticle *aodmc = (AliAODMCParticle*) arrayMC->At(imc);
    if(!aodmc) continue;

    if(aodmc->GetMother() >= 0) continue;

    TParticle *part = fUPCEvent->AddMCParticle();
    part->SetMomentum(aodmc->Px(), aodmc->Py(), aodmc->Pz(), aodmc->E());
    part->SetProductionVertex(aodmc->Xv(), aodmc->Yv(), aodmc->Zv(), aodmc->T());
    part->SetFirstMother(aodmc->GetMother());
    part->SetLastDaughter(aodmc->GetNDaughters());
    part->SetPdgCode(aodmc->GetPdgCode());
    part->SetUniqueID(imc);
  }//loop over mc particles

}//RunAODMC

//_____________________________________________________________________________
Bool_t AliAnalysisTaskUpcFilterSemiforward::RunESD()
{
  //cout<<"#################### ESD event ##################"<<endl;

  // Input ESD event
  AliESDEvent *esdEvent = (AliESDEvent*) InputEvent();
  if(!esdEvent) return kFALSE;

  fCounter->Fill( 21 ); // 21 = ESD analyzed events

  if(fIsMC) {
    //SPD FO fired chips
    fUPCEvent->SetNSPDfiredInner( fTriggerAna->SPDFiredChips(esdEvent,1,kFALSE,1) );
    fUPCEvent->SetNSPDfiredOuter( fTriggerAna->SPDFiredChips(esdEvent,1,kFALSE,2) );
    fUPCEvent->SetFastOrFiredChips( (TBits*) &esdEvent->GetMultiplicity()->GetFastOrFiredChips() );

    RunESDMC();
  }

  static AliESDVertex *vtx0 = 0x0;
  if( !vtx0 ) {
    // make a vertex at coordinates 0, 0, 0
    vtx0 = new AliESDVertex();
    vtx0->SetXv(0.); vtx0->SetYv(0.); vtx0->SetZv(0.);
  }

  // ESD central tracks
  Double_t pxpypz[3];
  Float_t b[2]; Float_t cov[3];
  UChar_t maskMan;
  UInt_t filterMap;
  Int_t nmun=0, ncen=0;
  //ESD central tracks loop
  for(Int_t itr=0; itr<esdEvent->GetNumberOfTracks(); itr++) {
    AliESDtrack *eTrack = esdEvent->GetTrack(itr);
    if( !eTrack ) continue;

    // manual filter
    maskMan = 0;
    if( eTrack->HasPointOnITSLayer(0) )                 { maskMan |= 1 << 0; }
    if( eTrack->HasPointOnITSLayer(1) )                 { maskMan |= 1 << 1; }
    if( eTrack->GetStatus() & AliESDtrack::kITSrefit )  { maskMan |= 1 << 2; }
    if( eTrack->GetStatus() & AliESDtrack::kTPCrefit )  { maskMan |= 1 << 3; }
    if( eTrack->GetTOFsignal() < 99999. && eTrack->GetTOFsignal() > 0. )  { maskMan |= 1 << 4; }
    if( eTrack->GetKinkIndex(0) > 0 )                   { maskMan |= 1 << 5; } // bit set if track is kink candidate

    //selection for at least one point in ITS and TPC refit
    //if( !(maskMan & (1 << 0)) && !(maskMan & (1 << 1)) ) continue;
    //if( !(maskMan & (1 << 3)) ) continue;

    //selection for its refit and tpc refit
    //if( !(maskMan & (1 << 2)) || !(maskMan & (1 << 3)) ) continue;

    //selection for tpc refit
    //if( !(maskMan & (1 << 3)) ) continue;

    //selection for at least one point in SPD and ITS refit and TPC refit
    //if( !(maskMan & (1 << 0)) && !(maskMan & (1 << 1)) ) continue;
    //if( !(maskMan & (1 << 2)) || !(maskMan & (1 << 3)) ) continue;

    //selection for at least one point in SPD and at least one TPC cluster
    if( !(maskMan & (1 << 0)) && !(maskMan & (1 << 1)) ) continue;
    if( eTrack->GetTPCNcls() == 0 ) continue;

    //central track accepted to write to the UPC event

    // various configurations of AliESDtrackCuts
    filterMap = 0;
    for(UInt_t i=0; i<32; i++) {
      if( !fCutsList[i] ) continue;

      if( fCutsList[i]->AcceptTrack(eTrack) ) { filterMap |= 1 << i; }
    }

    eTrack->GetPxPyPz(pxpypz);

    //fill the UPC track
    AliUPCTrack *upcTrack = fUPCEvent->AddTrack();
    fCounter->Fill( 32 ); // 32 = added central track
    upcTrack->SetPxPyPz( pxpypz );
    upcTrack->SetMaskMan( maskMan );
    upcTrack->SetFilterMap( filterMap );
    upcTrack->SetCharge( eTrack->Charge() );
    upcTrack->SetChi2perNDF( eTrack->GetTPCchi2()/((Double_t) eTrack->GetTPCNcls()) ); // TPC chi2 used to fill this data member in esd
    upcTrack->SetTPCmomentum( eTrack->GetTPCmomentum() );
    upcTrack->SetTPCsignal( eTrack->GetTPCsignal() );
    upcTrack->SetTPCNcls( eTrack->GetTPCNcls() );
    upcTrack->SetTPCCrossedRows( eTrack->GetTPCCrossedRows() );
    upcTrack->SetTPCNclsF( eTrack->GetTPCNclsF() );
    upcTrack->SetTPCNclsS( eTrack->GetTPCnclsS() );
    upcTrack->SetITSchi2perNDF( eTrack->GetITSchi2()/((Double_t) eTrack->GetNcls(0)) );
    upcTrack->SetITSClusterMap( eTrack->GetITSClusterMap() );
    upcTrack->SetTOFsignal( eTrack->GetTOFsignal() );

    //TPC PID
    // AliPID::EParticleType = { kElectron = 0,  kMuon = 1,  kPion = 2,  kKaon = 3,  kProton = 4 .... }
    if( fPIDResponse && (maskMan & (1 << 3)) ) {
      for(Int_t itype=0; itype<=4; itype++) {
        upcTrack->SetNSigmasTPC( itype, fPIDResponse->NumberOfSigmasTPC( eTrack, (AliPID::EParticleType) itype ) );
      }
    }
    //TOF PID
    if( fPIDResponse && (maskMan & (1 << 4)) ) {
      for(Int_t itype=0; itype<=4; itype++) {
        upcTrack->SetNSigmasTOF( itype, fPIDResponse->NumberOfSigmasTOF( eTrack, (AliPID::EParticleType) itype ) );
      }
    }

    //DCA to default primary vertex
    eTrack->GetImpactParameters(b, cov);
    upcTrack->SetImpactParameters(b, cov);

    //DCA to SPD vertex
    AliESDtrack *track1 = (AliESDtrack*) eTrack->Clone("track1");
    track1->RelateToVertex( esdEvent->GetPrimaryVertexSPD(), esdEvent->GetMagneticField(), 9999. );
    track1->GetImpactParameters(b, cov);
    upcTrack->SetImpactParametersSPD(b, cov);
    delete track1;

    //DCA to nominal interaction point
    track1 = (AliESDtrack*) eTrack->Clone("track1");
    track1->RelateToVertex( vtx0, esdEvent->GetMagneticField(), 9999. );
    track1->GetImpactParameters(b, cov);
    upcTrack->SetImpactParametersIP(b, cov);
    delete track1;

    ncen++;
  } //ESD central tracks loop

  // ESD muon tracks
  //muon tracks loop
  for(Int_t itr=0; itr<esdEvent->GetNumberOfMuonTracks(); itr++) {
    AliESDMuonTrack *mTrack = esdEvent->GetMuonTrack(itr);
    if( !mTrack ) continue;

    AliUPCMuonTrack *upcMuon = fUPCEvent->AddMuonTrack();
    fCounter->Fill( 31 ); // 31 = added muon track
    upcMuon->SetPtEtaPhi( mTrack->Pt(), mTrack->Eta(), mTrack->Phi() );
    upcMuon->SetCharge( mTrack->Charge() );
    upcMuon->SetMatchTrigger( mTrack->GetMatchTrigger() );
    upcMuon->SetRAtAbsorberEnd( mTrack->GetRAtAbsorberEnd() );
    upcMuon->SetChi2perNDF( mTrack->GetChi2()/((Double_t) mTrack->GetNDF()) );
    upcMuon->SetDCA( mTrack->GetDCA() );
    upcMuon->SetPxDCA( fMuonCuts->IsSelected(mTrack) );

    nmun++;
  } //muon tracks loop

  //selection for at least one muon or central track
  //if( nmun + ncen < 1 ) return kFALSE;

  //selection for at least one muon and at least one central track
  //if( nmun < 1 || ncen < 1 ) return kFALSE;

  // L0 trigger inputs
  fUPCEvent->SetL0inputs( esdEvent->GetHeader()->GetL0TriggerInputs() );

  //Tracklets
  fUPCEvent->SetNumberOfTracklets( esdEvent->GetMultiplicity()->GetNumberOfTracklets() );

  // ESD ZDC for TDC
  AliESDZDC *dataZDCESD = esdEvent->GetESDZDC();
  if(!dataZDCESD) {PostData(2, fHistList); return kFALSE;}

  Bool_t znctdc = kFALSE, zpctdc = kFALSE, znatdc = kFALSE, zpatdc = kFALSE;
  Int_t tdcsum[4] = {0,0,0,0};
  Float_t znctime = 0., znatime = 0.;
  for(Int_t iz=0;iz<4;iz++) {
    for(Int_t i=0; i<4; i++) tdcsum[i] += dataZDCESD->GetZDCTDCData(10+i,iz);
    znctime += dataZDCESD->GetZDCTDCCorrected(10,iz);
    znatime += dataZDCESD->GetZDCTDCCorrected(12,iz);
  }
  if( tdcsum[0] != 0 ) znctdc = kTRUE;
  if( tdcsum[1] != 0 ) zpctdc = kTRUE;
  if( tdcsum[2] != 0 ) znatdc = kTRUE;
  if( tdcsum[3] != 0 ) zpatdc = kTRUE;
/*

    if( dataZDCESD->GetZDCTDCData(10,iz) ) znctdc = kTRUE;
    if( dataZDCESD->GetZDCTDCData(11,iz) ) zpctdc = kTRUE;
    if( dataZDCESD->GetZDCTDCData(12,iz) ) znatdc = kTRUE;
    if( dataZDCESD->GetZDCTDCData(13,iz) ) zpatdc = kTRUE;

*/
  fUPCEvent->SetZNCtdc( znctdc );
  fUPCEvent->SetZPCtdc( zpctdc );
  fUPCEvent->SetZNAtdc( znatdc );
  fUPCEvent->SetZPAtdc( zpatdc );
  fUPCEvent->SetZNCtdcData( tdcsum[0] );
  fUPCEvent->SetZPCtdcData( tdcsum[1] );
  fUPCEvent->SetZNAtdcData( tdcsum[2] );
  fUPCEvent->SetZPAtdcData( tdcsum[3] );
  fUPCEvent->SetZNCTime( znctime );
  fUPCEvent->SetZNATime( znatime );

  //SPD primary vertex in ESD
  const AliESDVertex *vtx = esdEvent->GetPrimaryVertexSPD();
  if(!vtx) {PostData(2, fHistList); return kFALSE;}
  Double_t posVtx[3]; vtx->GetXYZ( posVtx );
  Double_t covVtx[6]; vtx->GetCovarianceMatrix( covVtx );
  fUPCEvent->SetPrimaryVertexSPD(posVtx, vtx->GetChi2perNDF(), covVtx, vtx->GetNContributors());
  const char *vtxtitle = vtx->GetTitle();
  fUPCEvent->SetPrimaryVertexSPDtitle( vtxtitle );

  return kTRUE;

}//RunESD

//_____________________________________________________________________________
void AliAnalysisTaskUpcFilterSemiforward::RunESDMC()
{
  // ESD MC particles

  AliMCEvent *mcEvent = MCEvent();
  if(!mcEvent) {cout<<"Error: no ESD MC found"<<endl; return;}

  //generated vertex
  TArrayF vtx(3);
  mcEvent->GenEventHeader()->PrimaryVertex(vtx);
  fUPCEvent->SetPrimaryVertexMC( vtx.GetArray() );

  //loop over mc particles
  for(Int_t imc=0; imc<mcEvent->GetNumberOfTracks(); imc++) {
    AliMCParticle *esdmc = (AliMCParticle*) mcEvent->GetTrack(imc);
    if(!esdmc) continue;

    if(esdmc->GetMother() >= 0) continue;

    TParticle *part = fUPCEvent->AddMCParticle();
    part->SetMomentum(esdmc->Px(), esdmc->Py(), esdmc->Pz(), esdmc->E());
    part->SetProductionVertex(esdmc->Xv(), esdmc->Yv(), esdmc->Zv(), 0.);
    part->SetFirstMother(esdmc->GetMother());
    part->SetLastDaughter(esdmc->GetLastDaughter()-esdmc->GetFirstDaughter()+1);
    part->SetPdgCode(esdmc->PdgCode());
    part->SetUniqueID(imc);
  }//loop over mc particles

}//RunESDMC

//_____________________________________________________________________________
void AliAnalysisTaskUpcFilterSemiforward::Terminate(Option_t *) 
{

  cout<<"Analysis complete."<<endl;
}//Terminate






























