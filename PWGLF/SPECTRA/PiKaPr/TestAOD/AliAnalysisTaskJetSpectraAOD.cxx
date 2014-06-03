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
  fFilterMask(0),
  fJetPtMin(0),
  fJetEtaMin(0x0),
  fJetEtaMax(0x0),
  fnCentBins(20),
  fnQvecBins(20),
  fIsQvecCalibMode(0),
  fQvecUpperLim(100),
  fIsQvecCut(0),
  fQvecMin(0),
  fQvecMax(100)
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
  Printf("\n\n\n\n\n\n In CreateOutput Object:");
  
  fOutput = new TList();
  fOutput->SetOwner();
  fOutput->SetName("chistpt");
  
  fListJets = new TList;
  
  if (!fTrackCuts) AliFatal("Track Cuts should be set in the steering macro");
  if (!fEventCuts) AliFatal("Event Cuts should be set in the steering macro");
    
  // binning common to all the THn
  const Double_t ptBins[] = {0.15,5.,10.,15.,20.,25.,30.,35.,40.,50.,75.,100.,150.,200.};
  const Int_t nptBins=13;
  
  //dimensions of THnSparse for jets
  const Int_t nvarjet=5;
  //                                        pt_raw    pt_corr           cent             Q vec            rho
  Int_t    binsHistRealJet[nvarjet] = {    nptBins,   nptBins,       fnCentBins,      fnQvecBins,          40.};
  Double_t xminHistRealJet[nvarjet] = {         0.,        0.,             0.,                0.,           0.};
  Double_t xmaxHistRealJet[nvarjet] = {       200.,      200.,           100.,     fQvecUpperLim,         200.};    
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
  NSparseHistJet->GetAxis(4)->SetTitle("rho");
  NSparseHistJet->GetAxis(4)->SetName("rho");
  fOutput->Add(NSparseHistJet);
  
  //dimensions of THnSparse for the normalization
  const Int_t nvarev=3;
  //                                             cent         Q vec         rho
  Int_t    binsHistRealEv[nvarev] = {    fnCentBins,      fnQvecBins,       40.};
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

// -- event selection --

  //event cuts
  if(!fEventCuts->IsSelected(fAOD,fTrackCuts))return;//event selection  //FIXME in our event selection need to put fAODJets?
  
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
  
  // accepted events  
  // -- end event selection --
  
  
  // get background
  AliAODJetEventBackground* externalBackground = 0;
  if(fAODJets && !externalBackground && fBackgroundBranch.Length()){
    externalBackground =  (AliAODJetEventBackground*)(fAODJets->FindListObject(fBackgroundBranch.Data()));
    if(!externalBackground)Printf("%s:%d Background branch not found %s",(char*)__FILE__,__LINE__,fBackgroundBranch.Data());;
  }
  
  Float_t rho = 0;
  if(externalBackground)rho = externalBackground->GetBackground(0);  //default schema
  
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

    Double_t ptcorr = ptJet-rho*areaJet;
           
    Double_t varJet[5];
    varJet[0]=jet->Pt();
    varJet[1]=(Double_t)ptcorr;
    varJet[2]=(Double_t)Cent;
    varJet[3]=(Double_t)Qvec;
    varJet[4]=(Double_t)rho;
    
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
void   AliAnalysisTaskJetSpectraAOD::Terminate(Option_t *)
{
  // Terminate
  printf("AliAnalysisTaskJetSpectraAOD: Terminate() \n");
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