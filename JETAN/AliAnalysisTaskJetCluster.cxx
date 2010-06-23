// **************************************
// Task used for the correction of determiantion of reconstructed jet spectra
// Compares input (gen) and output (rec) jets   
// *******************************************


/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
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

 
#include <TROOT.h>
#include <TRandom.h>
#include <TSystem.h>
#include <TInterpreter.h>
#include <TChain.h>
#include <TRandom.h>
#include <TFile.h>
#include <TKey.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TProfile.h>
#include <TList.h>
#include <TLorentzVector.h>
#include <TClonesArray.h>
#include  "TDatabasePDG.h"

#include "AliAnalysisTaskJetCluster.h"
#include "AliAnalysisManager.h"
#include "AliJetFinder.h"
#include "AliJetHeader.h"
#include "AliJetReader.h"
#include "AliJetReaderHeader.h"
#include "AliUA1JetHeaderV1.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliAODHandler.h"
#include "AliAODTrack.h"
#include "AliAODJet.h"
#include "AliAODMCParticle.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliStack.h"
#include "AliGenPythiaEventHeader.h"
#include "AliJetKineReaderHeader.h"
#include "AliGenCocktailEventHeader.h"
#include "AliInputEventHandler.h"


#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/AreaDefinition.hh"
#include "fastjet/JetDefinition.hh"
// get info on how fastjet was configured
#include "fastjet/config.h"


ClassImp(AliAnalysisTaskJetCluster)

AliAnalysisTaskJetCluster::AliAnalysisTaskJetCluster(): AliAnalysisTaskSE(),
  fAOD(0x0),
  fUseAODTrackInput(kFALSE),
  fUseAODMCInput(kFALSE),
  fUseGlobalSelection(kFALSE),
  fFilterMask(0),
  fTrackTypeRec(kTrackUndef),
  fTrackTypeGen(kTrackUndef),
  fAvgTrials(1),
  fExternalWeight(1),    
  fRecEtaWindow(0.5),
  fRparam(1.0), 
  fAlgorithm(fastjet::kt_algorithm),
  fStrategy(fastjet::Best),
  fRecombScheme(fastjet::BIpt_scheme),
  fAreaType(fastjet::active_area), 
  fh1Xsec(0x0),
  fh1Trials(0x0),
  fh1PtHard(0x0),
  fh1PtHardNoW(0x0),  
  fh1PtHardTrials(0x0),
  fh1NJetsRec(0x0),
  fh1NConstRec(0x0),
  fh1NConstLeadingRec(0x0),
  fh1PtJetsRecIn(0x0),
  fh1PtJetsLeadingRecIn(0x0),
  fh1PtJetConstRec(0x0),
  fh1PtJetConstLeadingRec(0x0),
  fh1PtTracksRecIn(0x0),
  fh1PtTracksLeadingRecIn(0x0),
  fh1NJetsRecRan(0x0),
  fh1NConstRecRan(0x0),
  fh1PtJetsLeadingRecInRan(0x0),
  fh1NConstLeadingRecRan(0x0),
  fh1PtJetsRecInRan(0x0),
  fh1PtTracksGenIn(0x0),
  fh1Nch(0x0),
  fh2NRecJetsPt(0x0),
  fh2NRecTracksPt(0x0),
  fh2NConstPt(0x0),
  fh2NConstLeadingPt(0x0),
  fh2JetPhiEta(0x0),
  fh2LeadingJetPhiEta(0x0),
  fh2JetEtaPt(0x0),
  fh2LeadingJetEtaPt(0x0),
  fh2TrackEtaPt(0x0),
  fh2LeadingTrackEtaPt(0x0),
  fh2JetsLeadingPhiEta(0x0),
  fh2JetsLeadingPhiPt(0x0),
  fh2TracksLeadingPhiEta(0x0),
  fh2TracksLeadingPhiPt(0x0),
  fh2TracksLeadingJetPhiPt(0x0),
  fh2JetsLeadingPhiPtW(0x0),
  fh2TracksLeadingPhiPtW(0x0),
  fh2TracksLeadingJetPhiPtW(0x0),
  fh2NRecJetsPtRan(0x0),
  fh2NConstPtRan(0x0),
  fh2NConstLeadingPtRan(0x0),
  fh2PtNch(0x0),
  fh2PtNchRan(0x0),
  fh2PtNchN(0x0),
  fh2PtNchNRan(0x0),
  fh2TracksLeadingJetPhiPtRan(0x0),
  fh2TracksLeadingJetPhiPtWRan(0x0),
  fHistList(0x0)  
{

}

AliAnalysisTaskJetCluster::AliAnalysisTaskJetCluster(const char* name):
  AliAnalysisTaskSE(name),
  fAOD(0x0),
  fUseAODTrackInput(kFALSE),
  fUseAODMCInput(kFALSE),
  fUseGlobalSelection(kFALSE),
  fFilterMask(0),
  fTrackTypeRec(kTrackUndef),
  fTrackTypeGen(kTrackUndef),
  fAvgTrials(1),
  fExternalWeight(1),    
  fRecEtaWindow(0.5),
  fRparam(1.0), 
  fAlgorithm(fastjet::kt_algorithm),
  fStrategy(fastjet::Best),
  fRecombScheme(fastjet::BIpt_scheme),
  fAreaType(fastjet::active_area), 
  fh1Xsec(0x0),
  fh1Trials(0x0),
  fh1PtHard(0x0),
  fh1PtHardNoW(0x0),  
  fh1PtHardTrials(0x0),
  fh1NJetsRec(0x0),
  fh1NConstRec(0x0),
  fh1NConstLeadingRec(0x0),
  fh1PtJetsRecIn(0x0),
  fh1PtJetsLeadingRecIn(0x0),
  fh1PtJetConstRec(0x0),
  fh1PtJetConstLeadingRec(0x0),
  fh1PtTracksRecIn(0x0),
  fh1PtTracksLeadingRecIn(0x0),
  fh1NJetsRecRan(0x0),
  fh1NConstRecRan(0x0),
  fh1PtJetsLeadingRecInRan(0x0),
  fh1NConstLeadingRecRan(0x0),
  fh1PtJetsRecInRan(0x0),
  fh1PtTracksGenIn(0x0),
  fh1Nch(0x0),
  fh2NRecJetsPt(0x0),
  fh2NRecTracksPt(0x0),
  fh2NConstPt(0x0),
  fh2NConstLeadingPt(0x0),
  fh2JetPhiEta(0x0),
  fh2LeadingJetPhiEta(0x0),
  fh2JetEtaPt(0x0),
  fh2LeadingJetEtaPt(0x0),
  fh2TrackEtaPt(0x0),
  fh2LeadingTrackEtaPt(0x0),
  fh2JetsLeadingPhiEta(0x0),
  fh2JetsLeadingPhiPt(0x0),
  fh2TracksLeadingPhiEta(0x0),
  fh2TracksLeadingPhiPt(0x0),
  fh2TracksLeadingJetPhiPt(0x0),
  fh2JetsLeadingPhiPtW(0x0),
  fh2TracksLeadingPhiPtW(0x0),
  fh2TracksLeadingJetPhiPtW(0x0),
  fh2NRecJetsPtRan(0x0),
  fh2NConstPtRan(0x0),
  fh2NConstLeadingPtRan(0x0),
  fh2PtNch(0x0),
  fh2PtNchRan(0x0),
  fh2PtNchN(0x0),
  fh2PtNchNRan(0x0),
  fh2TracksLeadingJetPhiPtRan(0x0),
  fh2TracksLeadingJetPhiPtWRan(0x0),
  fHistList(0x0)
{
  DefineOutput(1, TList::Class());  
}



Bool_t AliAnalysisTaskJetCluster::Notify()
{
  //
  // Implemented Notify() to read the cross sections
  // and number of trials from pyxsec.root
  // 
  return kTRUE;
}

void AliAnalysisTaskJetCluster::UserCreateOutputObjects()
{

  //
  // Create the output container
  //


  // Connect the AOD


  if (fDebug > 1) printf("AnalysisTaskJetCluster::UserCreateOutputObjects() \n");

  OpenFile(1);
  if(!fHistList)fHistList = new TList();

  Bool_t oldStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);

  //
  //  Histogram
    
  const Int_t nBinPt = 200;
  Double_t binLimitsPt[nBinPt+1];
  for(Int_t iPt = 0;iPt <= nBinPt;iPt++){
    if(iPt == 0){
      binLimitsPt[iPt] = 0.0;
    }
    else {// 1.0
      binLimitsPt[iPt] =  binLimitsPt[iPt-1] + 0.5;
    }
  }
  
  const Int_t nBinPhi = 90;
  Double_t binLimitsPhi[nBinPhi+1];
  for(Int_t iPhi = 0;iPhi<=nBinPhi;iPhi++){
    if(iPhi==0){
      binLimitsPhi[iPhi] = -1.*TMath::Pi();
    }
    else{
      binLimitsPhi[iPhi] = binLimitsPhi[iPhi-1] + 1/(Float_t)nBinPhi * TMath::Pi()*2;
    }
  }



  const Int_t nBinEta = 40;
  Double_t binLimitsEta[nBinEta+1];
  for(Int_t iEta = 0;iEta<=nBinEta;iEta++){
    if(iEta==0){
      binLimitsEta[iEta] = -2.0;
    }
    else{
      binLimitsEta[iEta] = binLimitsEta[iEta-1] + 0.1;
    }
  }

  const int nChMax = 100;

  fh1Xsec = new TProfile("fh1Xsec","xsec from pyxsec.root",1,0,1);
  fh1Xsec->GetXaxis()->SetBinLabel(1,"<#sigma>");

  fh1Trials = new TH1F("fh1Trials","trials root file",1,0,1);
  fh1Trials->GetXaxis()->SetBinLabel(1,"#sum{ntrials}");


  fh1NJetsRec = new TH1F("fh1NJetsRec","N reconstructed jets",120,-0.5,119.5);
  fh1NJetsRecRan = new TH1F("fh1NJetsRecRan","N reconstructed jets",120,-0.5,119.5);

  fh1NConstRec = new TH1F("fh1NConstRec","# jet constituents",120,-0.5,119.5);
  fh1NConstRecRan = new TH1F("fh1NConstRecRan","# jet constituents",120,-0.5,119.5);
  fh1NConstLeadingRec = new TH1F("fh1NConstLeadingRec","jet constituents",120,-0.5,119.5);
  fh1NConstLeadingRecRan = new TH1F("fh1NConstLeadingRecRan","jet constituents",120,-0.5,119.5);


  fh1PtHard = new TH1F("fh1PtHard","PYTHIA Pt hard;p_{T,hard}",nBinPt,binLimitsPt);
  fh1PtHardNoW = new TH1F("fh1PtHardNoW","PYTHIA Pt hard no weight;p_{T,hard}",nBinPt,binLimitsPt);
  fh1PtHardTrials = new TH1F("fh1PtHardTrials","PYTHIA Pt hard weight with trials;p_{T,hard}",nBinPt,binLimitsPt);

  fh1PtJetsRecIn  = new TH1F("fh1PtJetsRecIn","Rec jets P_T;p_{T} (GeV/c)",nBinPt,binLimitsPt);
  fh1PtJetsRecInRan  = new TH1F("fh1PtJetsRecInRan","Rec jets P_T;p_{T} (GeV/c)",nBinPt,binLimitsPt);
  fh1PtJetsLeadingRecIn = new TH1F("fh1PtJetsLeadingRecIn","Rec jets P_T;p_{T} (GeV/c)",nBinPt,binLimitsPt);
  fh1PtJetsLeadingRecInRan = new TH1F("fh1PtJetsLeadingRecInRan","Rec jets P_T;p_{T} (GeV/c)",nBinPt,binLimitsPt);
  fh1PtJetConstRec = new TH1F("fh1PtJetsConstRec","Rec jets constituents P_T;p_{T} (GeV/c)",nBinPt,binLimitsPt);
  fh1PtJetConstLeadingRec = new TH1F("fh1PtJetsConstLeadingRec","Rec jets constituents P_T;p_{T} (GeV/c)",nBinPt,binLimitsPt);
  fh1PtTracksRecIn  = new TH1F("fh1PtTracksRecIn","Rec tracks P_T #eta < 0.9;p_{T} (GeV/c)",nBinPt,binLimitsPt);
  fh1PtTracksLeadingRecIn  = new TH1F("fh1PtTracksLeadingRecIn","Rec tracks P_T #eta < 0.9;p_{T} (GeV/c)",nBinPt,binLimitsPt);
  fh1PtTracksGenIn  = new TH1F("fh1PtTracksGenIn","gen tracks P_T #eta < 0.9;p_{T} (GeV/c)",nBinPt,binLimitsPt);
  fh1Nch = new TH1F("fh1Nch","charged multiplicity; N_{ch}",nChMax,-0.5,nChMax-0.5);

  fh2NRecJetsPt = new TH2F("fh2NRecJetsPt","Number of jets above threshhold;p_{T,cut} (GeV/c);N_{jets}",nBinPt,binLimitsPt,50,-0.5,49.5);
  fh2NRecJetsPtRan = new TH2F("fh2NRecJetsPtRan","Number of jets above threshhold;p_{T,cut} (GeV/c);N_{jets}",nBinPt,binLimitsPt,50,-0.5,49.5);
  fh2NRecTracksPt = new TH2F("fh2NRecTracksPt","Number of tracks above threshhold;p_{T,cut} (GeV/c);N_{tracks}",nBinPt,binLimitsPt,50,-0.5,49.5);
  // 


  fh2NConstPt = new TH2F("fh2NConstPt","Number of constituents ;p_{T} (GeV/c);N",nBinPt,binLimitsPt,50,-0.5,49.5);
  fh2NConstLeadingPt = new TH2F("fh2NConstLeadingPt","Number of constituents ;p_{T} (GeV/c);N",nBinPt,binLimitsPt,50,-0.5,49.5);
  fh2NConstPtRan = new TH2F("fh2NConstPtRan","Number of constituents ;p_{T} (GeV/c);N",nBinPt,binLimitsPt,50,-0.5,49.5);
  fh2NConstLeadingPtRan = new TH2F("fh2NConstLeadingPtRan","Number of constituents ;p_{T} (GeV/c);N",nBinPt,binLimitsPt,50,-0.5,49.5);

  fh2PtNch = new TH2F("fh2PtNch","p_T of cluster vs. multiplicity; N_{ch};p_{T} (GeV/c);",nChMax,-0.5,nChMax-0.5,nBinPt,binLimitsPt);
  fh2PtNchRan = new TH2F("fh2PtNchRan","p_T of cluster vs. multiplicity ran; N_{ch};p_{T} (GeV/c);",nChMax,-0.5,nChMax-0.5,nBinPt,binLimitsPt);
  fh2PtNchN = new TH2F("fh2PtNchN","p_T of cluster vs. multiplicity N weighted; N_{ch};p_{T} (GeV/c);",nChMax,-0.5,nChMax-0.5,nBinPt,binLimitsPt);
  fh2PtNchNRan = new TH2F("fh2PtNchNRan","p_T of cluster vs. multiplicity N weighted ran; N_{ch};p_{T} (GeV/c);",nChMax,-0.5,nChMax-0.5,nBinPt,binLimitsPt);



  fh2JetPhiEta  = new TH2F("fh2JetPhiEta","eta vs phi all jets;#phi;#eta",
			   nBinPhi,0.,2.*TMath::Pi(),nBinEta,binLimitsEta);
  fh2LeadingJetPhiEta  = new TH2F("fh2LeadingJetPhiEta","eta vs phi leading jets;#phi;#eta",
				  nBinPhi,0.,2.*TMath::Pi(),nBinEta,binLimitsEta);

  fh2JetEtaPt  = new TH2F("fh2JetEtaPt","pt vs eta all jets;#eta;p_{T}",
			  nBinEta,binLimitsEta,nBinPt,binLimitsPt);
  fh2LeadingJetEtaPt  = new TH2F("fh2LeadingJetEtaPt","pT vs eta leading jets;#eta;p_{T}",
				 nBinEta,binLimitsEta,nBinPt,binLimitsPt);

  fh2TrackEtaPt  = new TH2F("fh2TrackEtaPt","pt vs eta all jets;#eta;p_{T}",
			  nBinEta,binLimitsEta,nBinPt,binLimitsPt);
  fh2LeadingTrackEtaPt  = new TH2F("fh2LeadingTrackEtaPt","pT vs eta leading jets;#eta;p_{T}",
				 nBinEta,binLimitsEta,nBinPt,binLimitsPt);



  fh2JetsLeadingPhiEta = new TH2F("fh2JetsLeadingPhiEta","delta eta vs delta phi to leading jet;#Delta#phi;#Delta#eta",
				nBinPhi,binLimitsPhi,nBinEta,binLimitsEta);
  fh2JetsLeadingPhiPt = new TH2F("fh2JetsLeadingPhiPt","leading p_T vs delta phi to leading jet;#Delta#phi;p_{T} (GeV/c)",
				nBinPhi,binLimitsPhi,nBinPt,binLimitsPt);
  fh2TracksLeadingPhiEta = new TH2F("fh2TracksLeadingPhiEta","delta eta vs delta phi to leading track;#Delta#phi;#Delta#eta",
				    nBinPhi,binLimitsPhi,nBinEta,binLimitsEta);
  fh2TracksLeadingPhiPt = new TH2F("fh2TracksLeadingPhiPt","leading p_T vs delta phi to leading jet;#Delta#phi;p_{T} (GeV/c)",
				 nBinPhi,binLimitsPhi,nBinPt,binLimitsPt);
  fh2TracksLeadingJetPhiPt = new TH2F("fh2TracksLeadingJetPhiPt","leading p_T vs delta phi to leading jet;#Delta#phi;p_{T} (GeV/c)",
				 nBinPhi,binLimitsPhi,nBinPt,binLimitsPt);
  fh2TracksLeadingJetPhiPtRan = new TH2F("fh2TracksLeadingJetPhiPtRan","leading p_T vs delta phi to leading jet;#Delta#phi;p_{T} (GeV/c)",
				 nBinPhi,binLimitsPhi,nBinPt,binLimitsPt);

  fh2JetsLeadingPhiPtW      = new TH2F("fh2JetsLeadingPhiPtW","leading p_T vs delta phi p_T weigted to leading jet;#Delta#phi;p_{T} (GeV/c)",
				nBinPhi,binLimitsPhi,nBinPt,binLimitsPt);
  fh2TracksLeadingPhiPtW    = new TH2F("fh2TracksLeadingPhiPtW","leading p_T vs delta phi to leading jet (p_T weighted);#Delta#phi;p_{T} (GeV/c)",
				    nBinPhi,binLimitsPhi,nBinPt,binLimitsPt);

  fh2TracksLeadingJetPhiPtW = new TH2F("fh2TracksLeadingJetPhiPtW","leading p_T vs delta phi to leading jet;#Delta#phi;p_{T} (GeV/c)",
				       nBinPhi,binLimitsPhi,nBinPt,binLimitsPt);
  fh2TracksLeadingJetPhiPtWRan = new TH2F("fh2TracksLeadingJetPhiPtWRan","leading p_T vs delta phi to leading jet;#Delta#phi;p_{T} (GeV/c)",
				       nBinPhi,binLimitsPhi,nBinPt,binLimitsPt);



  const Int_t saveLevel = 3; // large save level more histos
  if(saveLevel>0){
    fHistList->Add(fh1Xsec);
    fHistList->Add(fh1Trials);

    fHistList->Add(fh1NJetsRec);
    fHistList->Add(fh1NConstRec);
    fHistList->Add(fh1NConstLeadingRec);
    fHistList->Add(fh1PtJetsRecIn);
    fHistList->Add(fh1PtJetsLeadingRecIn);
    fHistList->Add(fh1PtTracksRecIn);
    fHistList->Add(fh1PtTracksLeadingRecIn);
    fHistList->Add(fh1PtJetConstRec);
    fHistList->Add(fh1PtJetConstLeadingRec);
    fHistList->Add(fh1NJetsRecRan);
    fHistList->Add(fh1NConstRecRan);
    fHistList->Add(fh1PtJetsLeadingRecInRan);
    fHistList->Add(fh1NConstLeadingRecRan);
    fHistList->Add(fh1PtJetsRecInRan);
    fHistList->Add(fh1Nch);
    fHistList->Add(fh2NRecJetsPt);
    fHistList->Add(fh2NRecTracksPt);
    fHistList->Add(fh2NConstPt);
    fHistList->Add(fh2NConstLeadingPt);
    fHistList->Add(fh2PtNch);
    fHistList->Add(fh2PtNchRan);
    fHistList->Add(fh2PtNchN);
    fHistList->Add(fh2PtNchNRan);
    fHistList->Add(fh2JetPhiEta);
    fHistList->Add(fh2LeadingJetPhiEta);
    fHistList->Add(fh2JetEtaPt);
    fHistList->Add(fh2LeadingJetEtaPt);
    fHistList->Add(fh2TrackEtaPt);
    fHistList->Add(fh2LeadingTrackEtaPt);
    fHistList->Add(fh2JetsLeadingPhiEta );
    fHistList->Add(fh2JetsLeadingPhiPt);
    fHistList->Add(fh2TracksLeadingPhiEta);
    fHistList->Add(fh2TracksLeadingPhiPt);
    fHistList->Add(fh2TracksLeadingJetPhiPt);
    fHistList->Add(fh2JetsLeadingPhiPtW);
    fHistList->Add(fh2TracksLeadingPhiPtW);
    fHistList->Add(fh2TracksLeadingJetPhiPtW);
    fHistList->Add(fh2NRecJetsPtRan);
    fHistList->Add(fh2NConstPtRan);
    fHistList->Add(fh2NConstLeadingPtRan);
    fHistList->Add(fh2TracksLeadingJetPhiPtRan);
    fHistList->Add(fh2TracksLeadingJetPhiPtWRan);
    }

  // =========== Switch on Sumw2 for all histos ===========
  for (Int_t i=0; i<fHistList->GetEntries(); ++i) {
    TH1 *h1 = dynamic_cast<TH1*>(fHistList->At(i));
    if (h1){
      h1->Sumw2();
      continue;
    }
    THnSparse *hn = dynamic_cast<THnSparse*>(fHistList->At(i));
    if(hn)hn->Sumw2();
  }
  TH1::AddDirectory(oldStatus);
}

void AliAnalysisTaskJetCluster::Init()
{
  //
  // Initialization
  //

  if (fDebug > 1) printf("AnalysisTaskJetCluster::Init() \n");

}

void AliAnalysisTaskJetCluster::UserExec(Option_t */*option*/)
{

  if(fUseGlobalSelection){
    // no selection by the service task, we continue
    if (fDebug > 1)Printf("Not selected %s:%d",(char*)__FILE__,__LINE__);
    PostData(1, fHistList);
    return;
  }

  //
  // Execute analysis for current event
  //
  AliESDEvent *fESD = 0;
  if(fUseAODTrackInput){    
    fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
    if(!fAOD){
      Printf("%s:%d AODEvent not found in Input Manager %d",(char*)__FILE__,__LINE__,fUseAODTrackInput);
      return;
    }
    // fethc the header
  }
  else{
    //  assume that the AOD is in the general output...
    fAOD  = AODEvent();
    if(!fAOD){
      Printf("%s:%d AODEvent not found in the Output",(char*)__FILE__,__LINE__);
      return;
    }
    if(fDebug>0){
      fESD = dynamic_cast<AliESDEvent*> (InputEvent());
    }
  }
  
  Bool_t selectEvent =  false;
  Bool_t physicsSelection = true; // fInputHandler->IsEventSelected();

  if(fAOD){
    const AliAODVertex *vtxAOD = fAOD->GetPrimaryVertex();
    if(vtxAOD->GetNContributors()>0){
	Float_t zvtx = vtxAOD->GetZ();
	Float_t yvtx = vtxAOD->GetY();
	Float_t xvtx = vtxAOD->GetX();
	Float_t r2   = yvtx*yvtx+xvtx*xvtx;  
	if(TMath::Abs(zvtx)<8.&&r2<1.){
	  if(physicsSelection)selectEvent = true;
	}
    }
  }
  if(!selectEvent){
    PostData(1, fHistList);
    return;
  }
  
  fh1Trials->Fill("#sum{ntrials}",1);
  

  if (fDebug > 10)Printf("%s:%d",(char*)__FILE__,__LINE__);

  // ==== General variables needed



  // we simply fetch the tracks/mc particles as a list of AliVParticles

  TList recParticles;
  TList genParticles;

  Int_t nT = GetListOfTracks(&recParticles,fTrackTypeRec);
  Float_t nCh = recParticles.GetEntries(); 
  fh1Nch->Fill(nCh);
  if(fDebug>2)Printf("%s:%d Selected Rec tracks: %d %d",(char*)__FILE__,__LINE__,nT,recParticles.GetEntries());
  nT = GetListOfTracks(&genParticles,fTrackTypeGen);
  if(fDebug>2)Printf("%s:%d Selected Gen tracks: %d %d",(char*)__FILE__,__LINE__,nT,genParticles.GetEntries());

  // find the jets....

  vector<fastjet::PseudoJet> inputParticlesRec;
  vector<fastjet::PseudoJet> inputParticlesRecRan;
  for(int i = 0; i < recParticles.GetEntries(); i++){
    AliVParticle *vp = (AliVParticle*)recParticles.At(i);
    // Carefull energy is not well determined in real data, should not matter for p_T scheme?
    // we talk total momentum here
    fastjet::PseudoJet jInp(vp->Px(),vp->Py(),vp->Pz(),vp->P());
    jInp.set_user_index(i);
    inputParticlesRec.push_back(jInp);

    // the randomized input changes eta and phi, but keeps the p_T
    Double_t pT = vp->Pt();
    Double_t eta = 1.8 * gRandom->Rndm() - 0.9;
    Double_t phi = 2.* TMath::Pi() * gRandom->Rndm();


    Double_t theta = 2.*TMath::ATan(TMath::Exp(-2.*eta));  
    Double_t pZ = pT/TMath::Tan(theta);

    Double_t pX = pT * TMath::Cos(phi);
    Double_t pY = pT * TMath::Sin(phi);
    Double_t p  = TMath::Sqrt(pT*pT+pZ*pZ); 

    fastjet::PseudoJet jInpRan(pX,pY,pZ,p);
    jInpRan.set_user_index(i);
    inputParticlesRecRan.push_back(jInpRan);

  }

  if(inputParticlesRec.size()==0){
    if(fDebug)Printf("%s:%d No input particles found, skipping event",(char*)__FILE__,__LINE__);
    PostData(1, fHistList);
    return;
  }
  
  // run fast jet
  fastjet::JetDefinition jetDef(fAlgorithm, fRparam, fRecombScheme, fStrategy);
  vector <fastjet::PseudoJet> inclusiveJets, sortedJets;
  fastjet::ClusterSequence clustSeq(inputParticlesRec, jetDef);
  
  inclusiveJets = clustSeq.inclusive_jets();
  sortedJets    = sorted_by_pt(inclusiveJets);
 
  if(fDebug)Printf("%s:%d Found %d jets",(char*)__FILE__,__LINE__, inclusiveJets.size());

  fh1NJetsRec->Fill(sortedJets.size());

 // loop over all jets an fill information, first on is the leading jet

 Int_t nRecOver = inclusiveJets.size();
 Int_t nRec     = inclusiveJets.size();
 if(inclusiveJets.size()>0){
   AliAODJet leadingJet (sortedJets[0].px(), sortedJets[0].py(), sortedJets[0].pz(), sortedJets[0].E());
    Float_t pt = leadingJet.Pt();

    Int_t iCount = 0;
    for(int i = 1;i <= fh2NRecJetsPt->GetNbinsX();i++){
      Float_t ptCut = fh2NRecJetsPt->GetXaxis()->GetBinCenter(i);
      while(pt<ptCut&&iCount<nRec){
	nRecOver--;
	iCount++;
	if(iCount<nRec){
	  pt = sortedJets[iCount].perp();
	}
      }
      if(nRecOver<=0)break;
      fh2NRecJetsPt->Fill(ptCut,nRecOver);
    }
    Float_t phi = leadingJet.Phi();
    if(phi<0)phi+=TMath::Pi()*2.;    
    Float_t eta = leadingJet.Eta();
    pt = leadingJet.Pt();

    // correlation of leading jet with tracks
    TIterator *recIter = recParticles.MakeIterator();
    AliVParticle *tmpRecTrack = (AliVParticle*)(recIter->Next());  

    recIter->Reset();
    while((tmpRecTrack = (AliVParticle*)(recIter->Next()))){
      Float_t tmpPt = tmpRecTrack->Pt();
      // correlation
      Float_t tmpPhi =  tmpRecTrack->Phi();     
      if(tmpPhi<0)tmpPhi+=TMath::Pi()*2.;    
      Float_t dPhi = phi - tmpPhi;
      if(dPhi>TMath::Pi())dPhi = dPhi - 2.*TMath::Pi();
      if(dPhi<(-1.*TMath::Pi()))dPhi = dPhi + 2.*TMath::Pi();      
      //      Float_t dEta = eta - tmpRecTrack->Eta();
      fh2TracksLeadingJetPhiPt->Fill(dPhi,pt);
      fh2TracksLeadingJetPhiPtW->Fill(dPhi,pt,tmpPt);
   }  



    for(int j = 0; j < nRec;j++){
      AliAODJet tmpRec (sortedJets[j].px(), sortedJets[j].py(), sortedJets[j].pz(), sortedJets[j].E());
      Float_t tmpPt = tmpRec.Pt();
      fh1PtJetsRecIn->Fill(tmpPt);
      // Fill Spectra with constituents
      vector<fastjet::PseudoJet> constituents = clustSeq.constituents(sortedJets[j]);
      fh1NConstRec->Fill(constituents.size());
      fh2PtNch->Fill(nCh,tmpPt);
      fh2PtNchN->Fill(nCh,tmpPt,constituents.size());
      fh2NConstPt->Fill(tmpPt,constituents.size());
      // loop over constiutents and fill spectrum
      for(unsigned int ic = 0; ic < constituents.size();ic++){
	AliVParticle *part = (AliVParticle*)recParticles.At(constituents[ic].user_index());
	fh1PtJetConstRec->Fill(part->Pt());
	if(j==0)fh1PtJetConstLeadingRec->Fill(part->Pt());
      }

      // correlation
      Float_t tmpPhi =  tmpRec.Phi();
      Float_t tmpEta =  tmpRec.Eta();
      if(tmpPhi<0)tmpPhi+=TMath::Pi()*2.;    

      if(j==0){
	fh1PtJetsLeadingRecIn->Fill(tmpPt);
	fh2LeadingJetPhiEta->Fill(tmpPhi,tmpEta);
	fh2LeadingJetEtaPt->Fill(tmpEta,tmpPt);
	fh1NConstLeadingRec->Fill(constituents.size());
	fh2NConstLeadingPt->Fill(tmpPt,constituents.size());
	continue;
      }
      fh2JetPhiEta->Fill(tmpRec.Phi(),tmpEta);
      fh2JetEtaPt->Fill(tmpEta,tmpPt);
      Float_t dPhi = phi - tmpPhi;
      if(dPhi>TMath::Pi())dPhi = dPhi - 2.*TMath::Pi();
      if(dPhi<(-1.*TMath::Pi()))dPhi = dPhi + 2.*TMath::Pi();      
      Float_t dEta = eta - tmpRec.Eta();
      fh2JetsLeadingPhiEta->Fill(dPhi,dEta);
      fh2JetsLeadingPhiPt->Fill(dPhi,pt);
      fh2JetsLeadingPhiPtW->Fill(dPhi,pt,tmpPt);
    }
    delete recIter;
 }
 

 // fill track information
 Int_t nTrackOver = recParticles.GetSize();
  // do the same for tracks and jets
 if(nTrackOver>0){
   TIterator *recIter = recParticles.MakeIterator();
   AliVParticle *tmpRec = (AliVParticle*)(recIter->Next());  
   Float_t pt = tmpRec->Pt();
   //    Printf("Leading track p_t %3.3E",pt);
   for(int i = 1;i <= fh2NRecTracksPt->GetNbinsX();i++){
     Float_t ptCut = fh2NRecTracksPt->GetXaxis()->GetBinCenter(i);
     while(pt<ptCut&&tmpRec){
       nTrackOver--;
       tmpRec = (AliVParticle*)(recIter->Next()); 
       if(tmpRec){
	 pt = tmpRec->Pt();
       }
     }
     if(nTrackOver<=0)break;
     fh2NRecTracksPt->Fill(ptCut,nTrackOver);
   }
   
   recIter->Reset();
   AliVParticle *leading = (AliVParticle*)recParticles.At(0);
   Float_t phi = leading->Phi();
   if(phi<0)phi+=TMath::Pi()*2.;    
   Float_t eta = leading->Eta();
   pt  = leading->Pt();
   while((tmpRec = (AliVParticle*)(recIter->Next()))){
     Float_t tmpPt = tmpRec->Pt();
     Float_t tmpEta = tmpRec->Eta();
     fh1PtTracksRecIn->Fill(tmpPt);
     fh2TrackEtaPt->Fill(tmpEta,tmpPt);
     if(tmpRec==leading){
       fh1PtTracksLeadingRecIn->Fill(tmpPt);
       fh2LeadingTrackEtaPt->Fill(tmpEta,tmpPt);
       continue;
     }
      // correlation
     Float_t tmpPhi =  tmpRec->Phi();
     
     if(tmpPhi<0)tmpPhi+=TMath::Pi()*2.;    
     Float_t dPhi = phi - tmpPhi;
     if(dPhi>TMath::Pi())dPhi = dPhi - 2.*TMath::Pi();
     if(dPhi<(-1.*TMath::Pi()))dPhi = dPhi + 2.*TMath::Pi();      
     Float_t dEta = eta - tmpRec->Eta();
     fh2TracksLeadingPhiEta->Fill(dPhi,dEta);
     fh2TracksLeadingPhiPt->Fill(dPhi,pt);
     fh2TracksLeadingPhiPtW->Fill(dPhi,pt,tmpPt);
   }  
   delete recIter;
 }

 // find the random jets
 vector <fastjet::PseudoJet> inclusiveJetsRan, sortedJetsRan;
 fastjet::ClusterSequence clustSeqRan(inputParticlesRecRan, jetDef);
  
 inclusiveJetsRan = clustSeqRan.inclusive_jets();
 sortedJetsRan    = sorted_by_pt(inclusiveJetsRan);

 // fill the jet information from random track

  fh1NJetsRecRan->Fill(sortedJetsRan.size());

 // loop over all jets an fill information, first on is the leading jet

 Int_t nRecOverRan = inclusiveJetsRan.size();
 Int_t nRecRan     = inclusiveJetsRan.size();
 if(inclusiveJetsRan.size()>0){
   AliAODJet leadingJet (sortedJetsRan[0].px(), sortedJetsRan[0].py(), sortedJetsRan[0].pz(), sortedJetsRan[0].E());
   Float_t pt = leadingJet.Pt();
   
   Int_t iCount = 0;
   for(int i = 1;i <= fh2NRecJetsPtRan->GetNbinsX();i++){
     Float_t ptCut = fh2NRecJetsPtRan->GetXaxis()->GetBinCenter(i);
      while(pt<ptCut&&iCount<nRecRan){
	nRecOverRan--;
	iCount++;
	if(iCount<nRecRan){
	  pt = sortedJetsRan[iCount].perp();
	}
      }
      if(nRecOverRan<=0)break;
      fh2NRecJetsPtRan->Fill(ptCut,nRecOverRan);
    }
    Float_t phi = leadingJet.Phi();
    if(phi<0)phi+=TMath::Pi()*2.;    
    pt = leadingJet.Pt();

    // correlation of leading jet with random tracks

    for(unsigned int ip = 0; ip < inputParticlesRecRan.size();ip++)
      { 
	Float_t tmpPt = inputParticlesRecRan[ip].perp();
	// correlation
	Float_t tmpPhi =  inputParticlesRecRan[ip].phi();
	if(tmpPhi<0)tmpPhi+=TMath::Pi()*2.;    
	Float_t dPhi = phi - tmpPhi;
	if(dPhi>TMath::Pi())dPhi = dPhi - 2.*TMath::Pi();
	if(dPhi<(-1.*TMath::Pi()))dPhi = dPhi + 2.*TMath::Pi();      
	//      Float_t dEta = eta - tmpRecTrack->Eta();
	fh2TracksLeadingJetPhiPtRan->Fill(dPhi,pt);
	fh2TracksLeadingJetPhiPtWRan->Fill(dPhi,pt,tmpPt);
      }  



    for(int j = 0; j < nRecRan;j++){
      AliAODJet tmpRec (sortedJetsRan[j].px(), sortedJetsRan[j].py(), sortedJetsRan[j].pz(), sortedJetsRan[j].E());
      Float_t tmpPt = tmpRec.Pt();
      fh1PtJetsRecInRan->Fill(tmpPt);
      // Fill Spectra with constituents
      vector<fastjet::PseudoJet> constituents = clustSeqRan.constituents(sortedJetsRan[j]);
      fh1NConstRecRan->Fill(constituents.size());
      fh2NConstPtRan->Fill(tmpPt,constituents.size());
      fh2PtNchRan->Fill(nCh,tmpPt);
      fh2PtNchNRan->Fill(nCh,tmpPt,constituents.size());
      // correlation
      Float_t tmpPhi =  tmpRec.Phi();
      if(tmpPhi<0)tmpPhi+=TMath::Pi()*2.;    

      if(j==0){
	fh1PtJetsLeadingRecInRan->Fill(tmpPt);
	fh1NConstLeadingRecRan->Fill(constituents.size());
	fh2NConstLeadingPtRan->Fill(tmpPt,constituents.size());
	continue;
      }
    }  
 }
 

 if (fDebug > 10)Printf("%s:%d",(char*)__FILE__,__LINE__);
 PostData(1, fHistList);
}

void AliAnalysisTaskJetCluster::Terminate(Option_t */*option*/)
{
// Terminate analysis
//
    if (fDebug > 1) printf("AnalysisJetCluster: Terminate() \n");
}


Int_t  AliAnalysisTaskJetCluster::GetListOfTracks(TList *list,Int_t type){

  if(fDebug>2)Printf("%s:%d Selecting tracks with %d",(char*)__FILE__,__LINE__,type);

  Int_t iCount = 0;
  if(type==kTrackAOD){
    AliAODEvent *aod = 0;
    if(fUseAODTrackInput)aod = dynamic_cast<AliAODEvent*>(InputEvent());
    else aod = AODEvent();
    if(!aod){
      return iCount;
    }
    for(int it = 0;it < aod->GetNumberOfTracks();++it){
      AliAODTrack *tr = aod->GetTrack(it);
      if((fFilterMask>0)&&!(tr->TestFilterBit(fFilterMask)))continue;
      if(TMath::Abs(tr->Eta())>0.9)continue;
      //      if(tr->Pt()<0.3)continue;
      list->Add(tr);
      iCount++;
    }
  }
  else if (type ==  kTrackKineAll||type == kTrackKineCharged){
    AliMCEvent* mcEvent = MCEvent();
    if(!mcEvent)return iCount;
    // we want to have alivpartilces so use get track
    for(int it = 0;it < mcEvent->GetNumberOfTracks();++it){
      if(!mcEvent->IsPhysicalPrimary(it))continue;
      AliMCParticle* part = (AliMCParticle*)mcEvent->GetTrack(it);
      if(type == kTrackKineAll){
	list->Add(part);
	iCount++;
      }
      else if(type == kTrackKineCharged){
	if(part->Particle()->GetPDG()->Charge()==0)continue;
	list->Add(part);
	iCount++;
      }
    }
  }
  else if (type == kTrackAODMCCharged || type == kTrackAODMCAll || type == kTrackAODMCChargedAcceptance) {
    AliAODEvent *aod = 0;
    if(fUseAODMCInput)aod = dynamic_cast<AliAODEvent*>(InputEvent());
    else aod = AODEvent();
    if(!aod)return iCount;
    TClonesArray *tca = dynamic_cast<TClonesArray*>(aod->FindListObject(AliAODMCParticle::StdBranchName()));
    if(!tca)return iCount;
    for(int it = 0;it < tca->GetEntriesFast();++it){
      AliAODMCParticle *part = dynamic_cast<AliAODMCParticle*>(tca->At(it));
      if(!part->IsPhysicalPrimary())continue;
      if(type == kTrackAODMCAll){
	list->Add(part);
	iCount++;
      }
      else if (type == kTrackAODMCCharged || type == kTrackAODMCChargedAcceptance ){
	if(part->Charge()==0)continue;
	if(kTrackAODMCCharged){
	  list->Add(part);
	}
	else {
	  if(TMath::Abs(part->Eta())>0.9)continue;
	  list->Add(part);
	}
	iCount++;
      }
    }
  }// AODMCparticle
  list->Sort();
  return iCount;

}

/*
Int_t AliAnalysisTaskJetCluster::AddParticlesFastJet(TList &particles,vector<fastjet::PseudoJet> &inputParticles){
  for(int i = 0; i < particles.GetEntries(); i++){
    AliVParticle *vp = (AliVParticle*)particles.At(i);
    // Carefull energy is not well determined in real data, should not matter for p_T scheme?
    fastjet::PseudoJet jInp(vp->Px(),vp->Py(),vp->Pz(),vp->E());
    jInp.set_user_index(i);
    inputParticles.push_back(jInp);
  }

  return 0;

}
*/
