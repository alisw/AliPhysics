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
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

//-----------------------------------------------------------------
//         AliAnalysisTaskJetSpectraAOD class
//-----------------------------------------------------------------

#include "TChain.h"
#include "TTree.h"
#include "TLegend.h"
#include "TH1F.h"
#include "TH2F.h"
#include "THnSparse.h"
#include "TCanvas.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliVTrack.h"
#include "AliAODMCParticle.h"
#include "AliVParticle.h"
#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include "AliAODInputHandler.h"
#include "AliAnalysisTaskJetSpectraAOD.h"
#include "AliAnalysisTaskESDfilter.h"
#include "AliAnalysisDataContainer.h"
#include "AliCentrality.h"
#include "TProof.h"
#include "AliVEvent.h"
#include "AliStack.h"
#include <TMCProcess.h>

#include "AliSpectraAODTrackCuts.h"
#include "AliSpectraAODEventCuts.h"

//jet include
#include "AliAODHandler.h"
#include "AliAODJetEventBackground.h"
#include "AliAODJet.h"
#include "AliAnalysisHelperJetTasks.h"

#include <iostream>

using namespace std;

ClassImp(AliAnalysisTaskJetSpectraAOD) 

//________________________________________________________________________
AliAnalysisTaskJetSpectraAOD::AliAnalysisTaskJetSpectraAOD(const char *name) : AliAnalysisTaskSE(name),
  fAOD(0),
  fIsMC(0),
  fEventCuts(0),
  fTrackCuts(0),
  fVZEROside(0),
  fOutput(0),
  fAODJets(0),
  fJetBranchName(""),
  fListJets(0),
  fBackgroundBranch(""),
  fOfflineTrgMask(AliVEvent::kMB),
  fFilterMask(0),
  fJetPtMin(0x0),
  fJetEtaMin(0x0),
  fJetEtaMax(0x0),
  fLeadPtMin(0x0),
  fnCentBins(20),
  fnQvecBins(20),
  fnptLeadBins(4),
  fIsQvecCalibMode(0),
  fQvecUpperLim(100),
  fIsQvecCut(0),
  fQvecMin(0),
  fQvecMax(100),
  fHistEvtSelection(0x0),
  fDebug(0),
  fMinNcontributors(0),
  fRejectPileup(0),
  fR(0.4),
  fZvertexDiff(1),
  fZvertex(10.)
{
  // Default constructor
  
  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
  DefineOutput(2, AliSpectraAODEventCuts::Class());
  DefineOutput(3, AliSpectraAODTrackCuts::Class());
  
}

//________________________________________________________________________
AliAnalysisTaskJetSpectraAOD::~AliAnalysisTaskJetSpectraAOD()
{
   delete fListJets;
}
//________________________________________________________________________
//________________________________________________________________________
void AliAnalysisTaskJetSpectraAOD::UserCreateOutputObjects()
{
  // create output objects
  fOutput = new TList();
  fOutput->SetOwner();
  fOutput->SetName("chistpt");
  
  fListJets = new TList;
  
  if (!fTrackCuts) AliFatal("Track Cuts should be set in the steering macro");
  if (!fEventCuts) AliFatal("Event Cuts should be set in the steering macro");
    
  // binning common to all the THn
  const Double_t ptBins[] = {-50.,-45.,-40.,-35.,-30.,-25.,-20.,-15.,-10.,-5.,0.,5.,10.,15.,20.,25.,30.,35.,40.,45.,50.,55.,60.,65.,70.,75.,80.,85.,90.,100.,150.,200.};
  const Int_t nptBins=31;
  
  //dimensions of THnSparse for jets
  const Int_t nvarjet=5;
  //                                        pt_raw    pt_corr           cent             Q vec                pt_lead
  Int_t    binsHistRealJet[nvarjet] = {    nptBins,   nptBins,       fnCentBins,      fnQvecBins,       fnptLeadBins};
  Double_t xminHistRealJet[nvarjet] = {         0.,        0.,             0.,                0.,                  0.};
  Double_t xmaxHistRealJet[nvarjet] = {       200.,      200.,           100.,     fQvecUpperLim,                 20.};    
  THnSparseF* NSparseHistJet = new THnSparseF("NSparseHistJet","NSparseHistJet",nvarjet,binsHistRealJet,xminHistRealJet,xmaxHistRealJet);
  NSparseHistJet->GetAxis(0)->SetTitle("#it{p}_{T,raw}");
  NSparseHistJet->GetAxis(0)->SetName("pT_raw");
  NSparseHistJet->SetBinEdges(0,ptBins);
  NSparseHistJet->GetAxis(1)->SetTitle("#it{p}_{T,corr}");
  NSparseHistJet->GetAxis(1)->SetName("pT_corr");
  NSparseHistJet->SetBinEdges(1,ptBins);
  NSparseHistJet->GetAxis(2)->SetTitle(Form("%s cent",fEventCuts->GetCentralityMethod().Data()));
  NSparseHistJet->GetAxis(2)->SetName(Form("%s_cent",fEventCuts->GetCentralityMethod().Data()));
  NSparseHistJet->GetAxis(3)->SetTitle("Q vec");
  NSparseHistJet->GetAxis(3)->SetName("Q_vec");
  NSparseHistJet->GetAxis(4)->SetTitle("#it{p}_{T,lead}");
  NSparseHistJet->GetAxis(4)->SetName("pT_lead");
  fOutput->Add(NSparseHistJet);
  
  //dimensions of THnSparse for the normalization
  const Int_t nvarev=3;
  //                                             cent         Q vec         rho
  Int_t    binsHistRealEv[nvarev] = {    fnCentBins,      fnQvecBins,      100 };
  Double_t xminHistRealEv[nvarev] = {           0.,               0.,        0.};
  Double_t xmaxHistRealEv[nvarev] = {         100.,      fQvecUpperLim,    200.};
  THnSparseF* NSparseHistEv = new THnSparseF("NSparseHistEv","NSparseHistEv",nvarev,binsHistRealEv,xminHistRealEv,xmaxHistRealEv);
  NSparseHistEv->GetAxis(0)->SetTitle(Form("%s cent",fEventCuts->GetCentralityMethod().Data()));
  NSparseHistEv->GetAxis(0)->SetName(Form("%s_cent",fEventCuts->GetCentralityMethod().Data()));
  NSparseHistEv->GetAxis(1)->SetTitle("Q vec");
  NSparseHistEv->GetAxis(1)->SetName("Q_vec");
  NSparseHistEv->GetAxis(2)->SetTitle("rho");
  NSparseHistEv->GetAxis(2)->SetName("rho");
  fOutput->Add(NSparseHistEv);
  
//   TH1F* fHistTest = new TH1F("fHistTest","fHistTest",nptBins-1,ptBins);
//   fOutput->Add(fHistTest);

  fHistEvtSelection = new TH1I("fHistEvtSelection", "event selection", 6, -0.5, 7.5);
  fHistEvtSelection->GetXaxis()->SetBinLabel(1,"ACCEPTED");
  fHistEvtSelection->GetXaxis()->SetBinLabel(2,"events IN");
  fHistEvtSelection->GetXaxis()->SetBinLabel(3,"event selection (rejected)");
  fHistEvtSelection->GetXaxis()->SetBinLabel(4,"centrality (rejected)");
  fHistEvtSelection->GetXaxis()->SetBinLabel(5,"multiplicity (rejected)");
  fHistEvtSelection->GetXaxis()->SetBinLabel(6,"ESE (rejected)");
  fOutput->Add(fHistEvtSelection);
  
  PostData(1, fOutput  );
  PostData(2, fEventCuts);
  PostData(3, fTrackCuts);
  
}

//________________________________________________________________________
void AliAnalysisTaskJetSpectraAOD::UserExec(Option_t *)
{
  
  // check for jet branches
  if(!strlen(fJetBranchName.Data())){
    AliError("Jet branch name not set.");
    return;
  }
  
  // main event loop
  fAOD = dynamic_cast<AliAODEvent*>(fInputEvent);
  if (!fAOD) {
    AliWarning("ERROR: AliAODEvent not available \n");
    return;
  }
    
  if (strcmp(fAOD->ClassName(), "AliAODEvent"))
    {
      AliFatal("Not processing AODs");
    }
  
  TObject* outHandler = AliAnalysisManager::GetAnalysisManager()->GetOutputEventHandler();
  if( outHandler && outHandler->InheritsFrom("AliAODHandler") ) {
    fAODJets = ((AliAODHandler*)outHandler)->GetAOD();
  }

  /* -- event selection -- */
  fHistEvtSelection->Fill(1); // number of events before event selection

  //jet service task event selection.
  Bool_t selected=kTRUE;
  selected = AliAnalysisHelperJetTasks::Selected();
  if(!selected){
    // no selection by the service task, we continue
    PostData(1,fOutput);
    PostData(2, fEventCuts);
    PostData(3, fTrackCuts);
  return;}
  
  // physics selection: this is now redundant, all should appear as accepted after service task selection
  AliInputEventHandler* inputHandler = (AliInputEventHandler*)
  ((AliAnalysisManager::GetAnalysisManager())->GetInputEventHandler());
  //std::cout<<inputHandler->IsEventSelected()<<" "<<fOfflineTrgMask<<std::endl;
  if(!(inputHandler->IsEventSelected() & fOfflineTrgMask)){
    if(fDebug) Printf(" Trigger Selection: event REJECTED ... ");
    fHistEvtSelection->Fill(2);
    PostData(1,fOutput);
    PostData(2, fEventCuts);
    PostData(3, fTrackCuts);
    return;
  }
  
  AliAODVertex* primVtx = fAOD->GetPrimaryVertex();
  Int_t nTracksPrim = primVtx->GetNContributors();
  
  if (fDebug) Printf("%s:%d primary vertex selection: %d", (char*)__FILE__,__LINE__,nTracksPrim);
  if(nTracksPrim < fMinNcontributors){
    if (fDebug) Printf("%s:%d primary vertex selection: event REJECTED...",(char*)__FILE__,__LINE__); 
    fHistEvtSelection->Fill(3);
    PostData(1,fOutput);
    PostData(2, fEventCuts);
    PostData(3, fTrackCuts);
    return;
  }
  
  TString primVtxName(primVtx->GetName());

  if(primVtxName.CompareTo("TPCVertex",TString::kIgnoreCase) == 1){
    if (fDebug) Printf("%s:%d primary vertex selection: TPC vertex, event REJECTED...",(char*)__FILE__,__LINE__);
    fHistEvtSelection->Fill(4);
    PostData(1,fOutput);
    PostData(2, fEventCuts);
    PostData(3, fTrackCuts);
    return;
  }

  if(fRejectPileup && AliAnalysisHelperJetTasks::IsPileUp()){
    if (fDebug) Printf("%s:%d SPD pileup: event REJECTED...",(char*)__FILE__,__LINE__);
    fHistEvtSelection->Fill(5);
    PostData(1,fOutput);
    PostData(2, fEventCuts);
    PostData(3, fTrackCuts);
    return;
  }
    
  //cut on difference of the z-vertex
  if (fZvertexDiff || (fZvertex>0)) {
    
    Double_t vzPRI = +999;
    Double_t vzSPD = -999;
    
    const AliVVertex *pv = fAOD->GetPrimaryVertex(); // AOD
    if (pv && pv->GetNContributors()>0) {
      vzPRI = pv->GetZ();
    }
    
    const AliVVertex *sv = fAOD->GetPrimaryVertexSPD();
    if (sv && sv->GetNContributors()>0) {
      vzSPD = sv->GetZ();
    }
    
    Double_t  dvertex = TMath::Abs(vzPRI-vzSPD);

    if (fZvertexDiff && (dvertex>0.1)) {
      if (fDebug) Printf("%s:%d ZvertexDiff Event selection: event REJECTED...",(char*)__FILE__,__LINE__);
      return;}
    if ((fZvertex>0) && (TMath::Abs(vzPRI)>fZvertex)) {
      if (fDebug) Printf("%s:%d ZvertexDiff Event selection: event REJECTED...",(char*)__FILE__,__LINE__);
      return;}
  }
  
  //ESE event cuts
  if(!fEventCuts->IsSelected(fAOD,fTrackCuts)){
    if (fDebug) Printf("%s:%d ESE Event selection: event REJECTED...",(char*)__FILE__,__LINE__);
    fHistEvtSelection->Fill(6);
    PostData(1,fOutput);
    PostData(2, fEventCuts);
    PostData(3, fTrackCuts);
    return;
  }
  
  if (fDebug) Printf("%s:%d event ACCEPTED ...",(char*)__FILE__,__LINE__); 
  fHistEvtSelection->Fill(0.);
  // accepted events  
  // -- end event selection --
  
  
  Double_t Qvec=0.;//in case of MC we save space in the memory
  if(!fIsMC){
    if(fIsQvecCalibMode){
      if(fVZEROside==0)Qvec=fEventCuts->GetqV0A();
      else if (fVZEROside==1)Qvec=fEventCuts->GetqV0C();
    }
    else Qvec=fEventCuts->GetQvecPercentile(fVZEROside);
  }
  
  if(fIsQvecCut && (Qvec<fQvecMin || Qvec>fQvecMax) ) return;
  
  Double_t Cent=fEventCuts->GetCent();

  // get background
  AliAODJetEventBackground* externalBackground = 0;
  if(fAODJets && !externalBackground && fBackgroundBranch.Length()){
    externalBackground =  (AliAODJetEventBackground*)(fAODJets->FindListObject(fBackgroundBranch.Data()));
    if(!externalBackground)Printf("%s:%d Background branch not found %s",(char*)__FILE__,__LINE__,fBackgroundBranch.Data());;
  }
  
  Float_t rho = 0;
  if(externalBackground)rho = externalBackground->GetBackground(0);  //default schema
    if(rho==0) rho=-9999.; //rho value = 0 are non-physical -> removed from the distribution
  
  // fetch jets 
  TClonesArray *aodJets = dynamic_cast<TClonesArray*>(fAODJets->FindListObject(fJetBranchName.Data()));
  if(!aodJets){
    AliError(Form("no jet branch \"%s\" found, in the AODs are:", fJetBranchName.Data()));
    if(fAOD){
      Printf("Input AOD >>>>");
      fAOD->Print();    
    }
    return;
  }
  
  Bool_t kDeltaPt = kFALSE;
  if(fJetBranchName.Contains("RandomCone")){
//     cout<<"RANDOM CONES BRANCH!"<<endl;
    kDeltaPt = kTRUE;
  }
  
  
  fListJets->Clear();
  for (Int_t iJet = 0; iJet < aodJets->GetEntriesFast(); iJet++) {
    AliAODJet *jet = dynamic_cast<AliAODJet*>((*aodJets)[iJet]);
    if (jet) fListJets->Add(jet);
  }
  


  for(Int_t i=0; i<fListJets->GetEntries(); ++i){
    AliAODJet* jet = (AliAODJet*)(fListJets->At(i));
    
    if((jet->Eta()<fJetEtaMin)||(jet->Eta()>fJetEtaMax)) continue;
    
    Double_t ptJet   = jet->Pt();
    
    Double_t areaJet = jet->EffectiveAreaCharged();

    Double_t ptcorr = -999.;
    if(kDeltaPt) ptcorr = ptJet - (rho*TMath::Pi()*fR*fR);
    else ptcorr = ptJet-(rho*areaJet);
           
    Double_t varJet[5];
    varJet[0]=jet->Pt();
    varJet[1]=(Double_t)ptcorr;
    varJet[2]=(Double_t)Cent;
    varJet[3]=(Double_t)Qvec;
    
    if(!kDeltaPt){
      AliVParticle *leadTrack = LeadingTrackFromJetRefs(jet);
      if(fLeadPtMin>0 && leadTrack->Pt()<fLeadPtMin)continue;
    
      varJet[4]=(Double_t)leadTrack->Pt();
    }
    else varJet[4] = -1;
    
    ((THnSparseF*)fOutput->FindObject("NSparseHistJet"))->Fill(varJet);//jet loop
  }
  
  
//   //track loop
//   for (Int_t iTracks = 0; iTracks < fAOD->GetNumberOfTracks(); iTracks++) {  //FIXME loop on aod track... should be on jet tracks???
//     AliAODTrack* track = fAOD->GetTrack(iTracks);
// //     if(!(track->TestFilterBit(1024)))continue;
//     if (!fTrackCuts->IsSelected(track,kTRUE)) continue;
//       
//     TH1F* h=(TH1F*)fOutput->FindObject("fHistTest");h->Fill(track->Pt());
//     
//   } // end loop on tracks
  
  Double_t varEv[3];
  varEv[0]=Cent;
  varEv[1]=Qvec;
  varEv[2]=rho;
  ((THnSparseF*)fOutput->FindObject("NSparseHistEv"))->Fill(varEv);//event loop
  
  PostData(1,fOutput);
  PostData(2, fEventCuts);
  PostData(3, fTrackCuts);
  //Printf("............. end of Exec");
  
}

//_________________________________________________________________
AliVParticle *AliAnalysisTaskJetSpectraAOD::LeadingTrackFromJetRefs(AliAODJet* jet){
  if(!jet)return 0;
  TRefArray *refs = jet->GetRefTracks();
  if(!refs) return 0;
  AliVParticle *leading = 0;
  Float_t fMaxPt = 0;
  for(int i = 0;i<refs->GetEntriesFast();i++){
    AliVParticle *tmp = dynamic_cast<AliVParticle*>(refs->At(i));
    if(!tmp)continue;
    if(tmp->Pt()>fMaxPt){
      leading = tmp;
      fMaxPt = tmp->Pt();
    }
  }
  return leading;
}

//_________________________________________________________________
void   AliAnalysisTaskJetSpectraAOD::Terminate(Option_t *)
{
  // Terminate
}

//jet
void AliAnalysisTaskJetSpectraAOD::SetBranchNames(const TString &branch)
{
   fJetBranchName = branch;
}

void AliAnalysisTaskJetSpectraAOD::SetRecBackgroundBranch(const TString &bckbranch)
{
   fBackgroundBranch = bckbranch;
}
