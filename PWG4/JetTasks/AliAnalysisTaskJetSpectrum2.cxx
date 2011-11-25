// **************************************
// used for the correction of determiantion of reconstructed jet spectra
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
#include <TRandom3.h>
#include <TSystem.h>
#include <TInterpreter.h>
#include <TChain.h>
#include <TFile.h>
#include <TKey.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TList.h>
#include <TLorentzVector.h>
#include <TClonesArray.h>
#include <TRefArray.h>
#include  "TDatabasePDG.h"

#include "AliAnalysisTaskJetSpectrum2.h"
#include "AliAnalysisManager.h"
#include "AliJetFinder.h"
#include "AliTHn.h"
#include "AliJetHeader.h"
#include "AliJetReader.h"
#include "AliJetReaderHeader.h"
#include "AliUA1JetHeaderV1.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliAODHandler.h"
#include "AliAODTrack.h"
#include "AliAODJet.h"
#include "AliAODJetEventBackground.h"
#include "AliAODMCParticle.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliStack.h"
#include "AliGenPythiaEventHeader.h"
#include "AliJetKineReaderHeader.h"
#include "AliGenCocktailEventHeader.h"
#include "AliInputEventHandler.h"


#include "AliAnalysisHelperJetTasks.h"

ClassImp(AliAnalysisTaskJetSpectrum2)

AliAnalysisTaskJetSpectrum2::AliAnalysisTaskJetSpectrum2(): 
  AliAnalysisTaskSE(),
  fJetHeaderRec(0x0),
  fJetHeaderGen(0x0),
  fAODIn(0x0),
  fAODOut(0x0),
  fAODExtension(0x0),
  fhnJetContainer(0x0),
  fhnCorrelation(0x0),
  fhnEvent(0x0),
  f1PtScale(0x0),
  fBranchRec("jets"),
  fBranchGen(""),
  fBranchBkgRec(""), 
  fBranchBkgGen(""), 
  fNonStdFile(""),
  fRandomizer(0x0),
  fUseAODJetInput(kFALSE),
  fUseAODTrackInput(kFALSE),
  fUseAODMCInput(kFALSE),
  fUseGlobalSelection(kFALSE),
  fUseExternalWeightOnly(kFALSE),
  fLimitGenJetEta(kFALSE),
  fDoMatching(kFALSE),
  fNMatchJets(5),
  fNRPBins(3),
  fJetTriggerExcludeMask(AliAODJet::kHighTrackPtTriggered),
  fFilterMask(0),
  fEventSelectionMask(0),
  fAnalysisType(0),
  fTrackTypeRec(kTrackUndef),
  fTrackTypeGen(kTrackUndef),
  fEventClass(0),
  fRPMethod(0),
  fAvgTrials(1),
  fExternalWeight(1),    
  fJetRecEtaWindow(0.5),
  fTrackRecEtaWindow(0.5),
  fMinJetPt(0),
  fMinTrackPt(0.15),
  fDeltaPhiWindow(90./180.*TMath::Pi()),
  fCentrality(100),
  fRPAngle(0),
  fMultRec(0),
  fMultGen(0),
  fh1Xsec(0x0),
  fh1Trials(0x0),
  fh1PtHard(0x0),
  fh1PtHardNoW(0x0),  
  fh1PtHardTrials(0x0),
  fh1ZVtx(0x0),
  fh1RP(0x0),
  fh1Centrality(0x0),
  fh1TmpRho(0x0),
  fh2MultRec(0x0),
  fh2MultGen(0x0),
  fh2RPCentrality(0x0),
  fh2PtFGen(0x0),
  fh2RelPtFGen(0x0),
  fHistList(0x0)  
{

  for(int ij = 0;ij <kJetTypes;++ij){    
    fFlagJetType[ij] = 1; // default = on
    fh1NJets[ij] = 0;
    fh1SumPtTrack[ij] = 0;
    fh1PtJetsIn[ij] = 0;
    fh1PtJetsInRej[ij] = 0;
    fh1PtTracksIn[ij] = 0;
    fh1PtTracksInLow[ij] = 0;
    fh2NJetsPt[ij]  = 0;
    fh2NTracksPt[ij]  = 0;
    fp2MultRPPhiTrackPt[ij] = 0;
    fp2CentRPPhiTrackPt[ij] = 0;
    fhnJetPt[ij] = 0;
    fhnJetPtQA[ij] = 0;
    fhnTrackPt[ij] = 0;
    fhnTrackPtQA[ij] = 0;
    for(int i = 0;i <= kMaxJets;++i){
      fh2LTrackPtJetPt[ij][i] = 0;
      fh1PtIn[ij][i] = 0;
    }

    fh1DijetMinv[ij] = 0;      
    fh2DijetDeltaPhiPt[ij] = 0;  
    fh2DijetAsymPt[ij] = 0; 
    fh2DijetPt2vsPt1[ij] = 0;
    fh2DijetDifvsSum[ij] = 0;

  }
}

AliAnalysisTaskJetSpectrum2::AliAnalysisTaskJetSpectrum2(const char* name):
  AliAnalysisTaskSE(name),
  fJetHeaderRec(0x0),
  fJetHeaderGen(0x0),
  fAODIn(0x0),
  fAODOut(0x0),
  fAODExtension(0x0),
  fhnJetContainer(0x0),
  fhnCorrelation(0x0),
  fhnEvent(0x0),
  f1PtScale(0x0),
  fBranchRec("jets"),
  fBranchGen(""),
  fBranchBkgRec(""),
  fBranchBkgGen(""),
  fNonStdFile(""),
  fRandomizer(0x0),
  fUseAODJetInput(kFALSE),
  fUseAODTrackInput(kFALSE),
  fUseAODMCInput(kFALSE),
  fUseGlobalSelection(kFALSE),
  fUseExternalWeightOnly(kFALSE),
  fLimitGenJetEta(kFALSE),
  fDoMatching(kFALSE),
  fNMatchJets(5),
  fNRPBins(3),
  fJetTriggerExcludeMask(AliAODJet::kHighTrackPtTriggered),
  fFilterMask(0),
  fEventSelectionMask(0),
  fAnalysisType(0),
  fTrackTypeRec(kTrackUndef),
  fTrackTypeGen(kTrackUndef),
  fEventClass(0),
  fRPMethod(0),
  fAvgTrials(1),
  fExternalWeight(1),    
  fJetRecEtaWindow(0.5),
  fTrackRecEtaWindow(0.5),
  fMinJetPt(0),
  fMinTrackPt(0.15),
  fDeltaPhiWindow(90./180.*TMath::Pi()),
  fCentrality(100),
  fRPAngle(0),
  fMultRec(0),
  fMultGen(0),
  fh1Xsec(0x0),
  fh1Trials(0x0),
  fh1PtHard(0x0),
  fh1PtHardNoW(0x0),  
  fh1PtHardTrials(0x0),
  fh1ZVtx(0x0),
  fh1RP(0x0),
  fh1Centrality(0x0),
  fh1TmpRho(0x0),
  fh2MultRec(0x0),
  fh2MultGen(0x0),
  fh2RPCentrality(0x0),
  fh2PtFGen(0x0),
  fh2RelPtFGen(0x0),
  fHistList(0x0)
{

  for(int ij = 0;ij <kJetTypes;++ij){    
    fFlagJetType[ij] = 1; // default = on
    fh1NJets[ij] = 0;
    fh1SumPtTrack[ij] = 0;
    fh1PtJetsIn[ij] = 0;
    fh1PtJetsInRej[ij] = 0;
    fh1PtTracksIn[ij] = 0;
    fh1PtTracksInLow[ij] = 0;
    fh2NJetsPt[ij]  = 0;
    fh2NTracksPt[ij]  = 0;
    fp2MultRPPhiTrackPt[ij] = 0;
    fp2CentRPPhiTrackPt[ij] = 0;
    fhnJetPt[ij] = 0;
    fhnJetPtQA[ij] = 0;
    fhnTrackPt[ij] = 0;
    fhnTrackPtQA[ij] = 0;
    for(int i = 0;i <= kMaxJets;++i){
      fh2LTrackPtJetPt[ij][i] = 0;
      fh1PtIn[ij][i] = 0;
    }

    fh1DijetMinv[ij] = 0;      
    fh2DijetDeltaPhiPt[ij] = 0;  
    fh2DijetAsymPt[ij] = 0; 
    fh2DijetPt2vsPt1[ij] = 0;
    fh2DijetDifvsSum[ij] = 0;
  } 

  DefineOutput(1, TList::Class());  
}



Bool_t AliAnalysisTaskJetSpectrum2::Notify()
{



  //
  // Implemented Notify() to read the cross sections
  // and number of trials from pyxsec.root
  // 
  
  // Fetch the aod also from the input in,
  // have todo it in notify
  
  
  fAODIn = dynamic_cast<AliAODEvent*>(InputEvent());
  //  assume that the AOD is in the general output...
  fAODOut  = AODEvent();
  
  if(fNonStdFile.Length()!=0){
    // case that we have an AOD extension we need can fetch the jets from the extended output
    AliAODHandler *aodH = dynamic_cast<AliAODHandler*>(AliAnalysisManager::GetAnalysisManager()->GetOutputEventHandler());
    fAODExtension = (aodH?aodH->GetExtension(fNonStdFile.Data()):0);    
    if(!fAODExtension){
      if(fDebug>1)Printf("AODExtension found for %s",fNonStdFile.Data());
    }
  }


  TTree *tree = AliAnalysisManager::GetAnalysisManager()->GetTree();
  Float_t xsection = 0;
  Float_t ftrials  = 1;

  fAvgTrials = 1;
  if(tree){
    TFile *curfile = tree->GetCurrentFile();
    if (!curfile) {
      Error("Notify","No current file");
      return kFALSE;
    }
    if(!fh1Xsec||!fh1Trials){
      Printf("%s%d No Histogram fh1Xsec",(char*)__FILE__,__LINE__);
      return kFALSE;
    }
    AliAnalysisHelperJetTasks::PythiaInfoFromFile(curfile->GetName(),xsection,ftrials);
    fh1Xsec->Fill("<#sigma>",xsection);
    // construct a poor man average trials 
    Float_t nEntries = (Float_t)tree->GetTree()->GetEntries();
    if(ftrials>=nEntries && nEntries>0.)fAvgTrials = ftrials/nEntries;
  }  

  if(fDebug)Printf("Reading File %s",fInputHandler->GetTree()->GetCurrentFile()->GetName());

  return kTRUE;
}

void AliAnalysisTaskJetSpectrum2::UserCreateOutputObjects()
{

  
  // Connect the AOD

  if (fDebug > 1) printf("AnalysisTaskJetSpectrum2::UserCreateOutputObjects() \n");
  OpenFile(1);
  if(!fHistList)fHistList = new TList(); 
  PostData(1, fHistList); // post data in any case once

  if(!fRandomizer)fRandomizer = new TRandom3(0);

  fHistList->SetOwner(kTRUE);
  Bool_t oldStatus = TH1::AddDirectoryStatus(); 
  TH1::AddDirectory(kFALSE);


  
  // event npsparse cent, mult
  const Int_t nBinsSparse0 = 2;
  const Int_t nBins0[nBinsSparse0] = {     100, 500};
  const Double_t xmin0[nBinsSparse0]  = {    0,   0};
  const Double_t xmax0[nBinsSparse0]  = {  100,5000};
      

  fhnEvent = new THnSparseF("fhnEvent",";cent;mult",nBinsSparse0,
			    nBins0,xmin0,xmax0);
  fHistList->Add(fhnEvent);

  if(fDoMatching){
    MakeJetContainer();
    fHistList->Add(fhnCorrelation);
    fHistList->Add(fhnJetContainer);
  }
  //
  //  Histogram
    


  const Int_t nBinPt = 120;
  Double_t binLimitsPt[nBinPt+1];
  for(Int_t iPt = 0;iPt <= nBinPt;iPt++){
    if(iPt == 0){
      binLimitsPt[iPt] = -50.0;
    }
    else {// 1.0
      binLimitsPt[iPt] =  binLimitsPt[iPt-1] + 2.5;
    }
  }
  const Int_t nBinPhi = 90;
  Double_t binLimitsPhi[nBinPhi+1];
  for(Int_t iPhi = 0;iPhi<=nBinPhi;iPhi++){
    if(iPhi==0){
      binLimitsPhi[iPhi] = 0;
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


  fh1Xsec = new TProfile("fh1Xsec","xsec from pyxsec.root",1,0,1);
  fh1Xsec->GetXaxis()->SetBinLabel(1,"<#sigma>");
  fHistList->Add(fh1Xsec);
  fh1Trials = new TH1F("fh1Trials","trials root file",1,0,1);
  fh1Trials->GetXaxis()->SetBinLabel(1,"#sum{ntrials}");
  fHistList->Add(fh1Trials);
  fh1PtHard = new TH1F("fh1PtHard","PYTHIA Pt hard;p_{T,hard}",nBinPt,binLimitsPt);
  fHistList->Add(fh1PtHard);
  fh1PtHardNoW = new TH1F("fh1PtHardNoW","PYTHIA Pt hard no weight;p_{T,hard}",nBinPt,binLimitsPt);
  fHistList->Add(fh1PtHardNoW);
  fh1PtHardTrials = new TH1F("fh1PtHardTrials","PYTHIA Pt hard weight with trials;p_{T,hard}",nBinPt,binLimitsPt);
  fHistList->Add(fh1PtHardTrials);
  
  fh1ZVtx = new TH1F("fh1ZVtx","z vtx;z_{vtx} (cm)",400,-20,20);
  fHistList->Add(fh1ZVtx);


  fh1RP = new TH1F("fh1RP","RP;#Psi",440, -1.*TMath::Pi(), 2.*TMath::Pi());
  fHistList->Add(fh1RP);

  fh1Centrality = new TH1F("fh1Centrality","cent;cent (%)",103,-1,102);
  fHistList->Add(fh1Centrality);

  fh2MultRec = new TH2F("fh2MultRec","multiplicity rec;# tracks;# jetrefs",500,0,5000,500,0.,5000);
  fHistList->Add(fh2MultRec);
  fh2MultGen = new TH2F("fh2MultGen","multiplicity gen;# tracks;# jetrefs",400,0,5000,500,0.,5000);
  fHistList->Add(fh2MultGen);


  fh2RPCentrality = new TH2F("fh2RPCentrality" ,"Reaction Plane Angle" , 20, 0.,100., 180, 0, TMath::Pi());
  fHistList->Add(fh2RPCentrality);

  fh2PtFGen = new TH2F("fh2PtFGen",Form("%s vs. %s;p_{T,gen};p_{T,rec}",fBranchRec.Data(),fBranchGen.Data()),nBinPt,binLimitsPt,nBinPt,binLimitsPt);
  fHistList->Add(fh2PtFGen);

  fh2RelPtFGen = new TH2F("fh2RelPtFGen",";p_{T,gen};p_{T,rec}-p_{T,gen}/p_{T,Gen}",nBinPt,binLimitsPt,241,-2.41,2.41);
  fHistList->Add(fh2RelPtFGen);
  
  for(int ij = 0;ij <kJetTypes;++ij){    
    TString cAdd = "";
    TString cJetBranch = "";
    if(ij==kJetRec){
      cAdd = "Rec";
      cJetBranch = fBranchRec.Data();
    }
    else if (ij==kJetGen){
      cAdd = "Gen";
      cJetBranch = fBranchGen.Data();
    }
    else if (ij==kJetRecFull){
      cAdd = "RecFull";
      cJetBranch = fBranchRec.Data();
    }
    else if (ij==kJetGenFull){
      cAdd = "GenFull";
      cJetBranch = fBranchGen.Data();
    }

    if(cJetBranch.Length()==0)fFlagJetType[ij] = 0;
    if(!fFlagJetType[ij])continue;

    fh1NJets[ij] =new TH1F(Form("fh1N%sJets",cAdd.Data()),Form("N %s jets",cAdd.Data()),50,-0.5,49.5);
    fHistList->Add(fh1NJets[ij]);
    
    fh1PtJetsIn[ij]  = new TH1F(Form("fh1PtJets%sIn",cAdd.Data()),Form("%s jets p_T;p_{T} (GeV/c)",cAdd.Data()),nBinPt,binLimitsPt);
    fHistList->Add(fh1PtJetsIn[ij]);

    fh1PtJetsInRej[ij]  = new TH1F(Form("fh1PtJets%sInRej",cAdd.Data()),Form("%s jets p_T;p_{T} (GeV/c)",cAdd.Data()),nBinPt,binLimitsPt);
    fHistList->Add(fh1PtJetsInRej[ij]);
    
    fh1PtTracksIn[ij] = new TH1F(Form("fh1PtTracks%sIn",cAdd.Data()),Form("%s track p_T;p_{T} (GeV/c)",cAdd.Data()),nBinPt,binLimitsPt);
    fHistList->Add(fh1PtTracksIn[ij]);
    
    fh1PtTracksInLow[ij] = new TH1F(Form("fh1PtTracks%sInLow",cAdd.Data()),Form("%s track p_T;p_{T} (GeV/c)",cAdd.Data()),100,0.,5.);
    fHistList->Add(fh1PtTracksInLow[ij]);
    
    fh1SumPtTrack[ij] = new TH1F(Form("fh1SumPtTrack%s",cAdd.Data()),Form("Sum %s track p_T;p_{T} (GeV/c)",cAdd.Data()),1000,0.,3000.);
    fHistList->Add(fh1SumPtTrack[ij]);
    
    fh2NJetsPt[ij]  = new TH2F(Form("fh2N%sJetsPt",cAdd.Data()),Form("Number of %s jets above threshhold;p_{T,cut} (GeV/c);N_{jets}",cAdd.Data()),nBinPt,binLimitsPt,50,-0.5,49.5);
    fHistList->Add(fh2NJetsPt[ij]);
    
    fh2NTracksPt[ij]  = new TH2F(Form("fh2N%sTracksPt",cAdd.Data()),Form("Number of %s tracks above threshhold;p_{T,cut} (GeV/c);N_{tracks}",cAdd.Data()),nBinPt,binLimitsPt,1000,0.,5000);
    fHistList->Add(fh2NTracksPt[ij]);

    
    fp2MultRPPhiTrackPt[ij] = new TProfile2D(Form("fp2MultRPPhiTrackPt%s",cAdd.Data()),"RP phi vs Number of tracks;# tracks;#Delta#phi_{RP}; <p_{T}>",20,0,4000,181,-1./180.*TMath::Pi(),TMath::Pi(),"S");
    fHistList->Add(fp2MultRPPhiTrackPt[ij]);
    fp2CentRPPhiTrackPt[ij] = new TProfile2D(Form("fp2CentRPPhiTrackPt%s",cAdd.Data()),"RP phi vs cent;# cent;#Delta#phi_{RP}; <p_{T}>",10,0,100,181,-1./180.*TMath::Pi(),TMath::Pi(),"S");
    fHistList->Add(fp2CentRPPhiTrackPt[ij]);    

    // Bins:  Jet number: pTJet, cent, mult, RP, Area.   total bins = 4.5M
    const Int_t nBinsSparse1 = 6;
    const Int_t nBins1[nBinsSparse1] = {     kMaxJets+1,120, 10,  25,    fNRPBins, 10};
    const Double_t xmin1[nBinsSparse1]  = {        -0.5,-50,  0,   0,        -0.5, 0.};
    const Double_t xmax1[nBinsSparse1]  = {kMaxJets+0.5,250,100,5000,fNRPBins-0.5,1.0};
    
    fhnJetPt[ij] = new THnSparseF(Form("fhnJetPt%s",cAdd.Data()),";jet number;p_{T,jet};cent;# tracks;RP;area",nBinsSparse1,nBins1,xmin1,xmax1);
    fHistList->Add(fhnJetPt[ij]);
    
    // Bins:  Jet number: pTJet, cent, eta, phi, Area.   total bins = 9.72 M
    const Int_t nBinsSparse2 = 6;
      const Int_t nBins2[nBinsSparse2] = {     kMaxJets+1, 25,   5,  18,             360, 10};
      const Double_t xmin2[nBinsSparse2]  = {        -0.5,  0,   0,-0.9,              0,  0.};
      const Double_t xmax2[nBinsSparse2]  = {kMaxJets+0.5,250, 100, 0.9, 2.*TMath::Pi(),1.0};
      fhnJetPtQA[ij] = new THnSparseF(Form("fhnJetPtQA%s",cAdd.Data()),";jet number;p_{T,jet};cent;#eta;#phi;area",nBinsSparse2,nBins2,xmin2,xmax2);
      fHistList->Add(fhnJetPtQA[ij]);

      // Bins:track number  pTtrack, cent, mult, RP.   total bins = 224 k
      const Int_t nBinsSparse3 = 5;
      const Int_t nBins3[nBinsSparse3] = {       2,    100,     10,  20,    fNRPBins};
      const Double_t xmin3[nBinsSparse3]  = { -0.5,     0,   0,      0,        -0.5};
      const Double_t xmax3[nBinsSparse3]  = { 1.5,    200, 100,   4000,fNRPBins-0.5};  

      // change the binning ot the pT axis:
      Double_t *xPt3 = new Double_t[nBins3[1]+1];
      xPt3[0] = 0.;
      for(int i = 1; i<=nBins3[1];i++){
	if(xPt3[i-1]<2)xPt3[i] = xPt3[i-1] + 0.05; // 1 - 40
	else if(xPt3[i-1]<4)xPt3[i] = xPt3[i-1] + 0.2; // 41 - 50
	else if(xPt3[i-1]<10)xPt3[i] = xPt3[i-1] + 0.5; // 50 - 62
	else if(xPt3[i-1]<20)xPt3[i] = xPt3[i-1] +  1.; // 62 - 72
	else if(xPt3[i-1]<30)xPt3[i] = xPt3[i-1] + 2.5; // 74 - 78
	else xPt3[i] = xPt3[i-1] + 5.; // 78 - 100 = 140 
      }
      
      fhnTrackPt[ij] = new THnSparseF(Form("fhnTrackPt%s",cAdd.Data()),";track number;p_{T};cent;#tracks;RP",nBinsSparse3,nBins3,xmin3,xmax3);
      fhnTrackPt[ij]->SetBinEdges(1,xPt3);
      fHistList->Add(fhnTrackPt[ij]);
      delete [] xPt3;

      // Track QA bins track nr, pTrack, cent, eta, phi bins 5.4 M
      const Int_t nBinsSparse4 = 5;
      const Int_t nBins4[nBinsSparse4] =    {    2, 50,  10,  20, 360};
      const Double_t xmin4[nBinsSparse4]  = { -0.5,  0,   0, -1.0,   0.};
      const Double_t xmax4[nBinsSparse4]  = {  1.5,150, 100,  1.0,2.*TMath::Pi()};  

      // change the binning ot the pT axis:
      Double_t *xPt4 = new Double_t[nBins4[1]+1];
      xPt4[0] = 0.;
      for(int i = 1; i<=nBins4[1];i++){
	if(xPt4[i-1]<2)xPt4[i] = xPt4[i-1] + 0.1;
	else if(xPt4[i-1]<10)xPt4[i] = xPt4[i-1] + 0.5;
	else if(xPt4[i-1]<20)xPt4[i] = xPt4[i-1] +  1.;
	else if(xPt4[i-1]<30)xPt4[i] = xPt4[i-1] +  2.5;
	else xPt4[i] = xPt4[i-1] + 5.;
      }
      fhnTrackPtQA[ij] = new THnSparseF(Form("fhnTrackPtQA%s",cAdd.Data()),";track number;p_{T};cent;#eta;#phi",nBinsSparse4,nBins4,xmin4,xmax4);
      fhnTrackPtQA[ij]->SetBinEdges(1,xPt4);
      fHistList->Add(fhnTrackPtQA[ij]);
      delete [] xPt4;

      for(int i = 0;i <= kMaxJets;++i){
	fh1PtIn[ij][i] = new TH1F(Form("fh1Pt%sIn_j%d",cAdd.Data(),i),Form("%s p_T input ;p_{T}",cAdd.Data()),nBinPt,binLimitsPt);
	fHistList->Add(fh1PtIn[ij][i]);


	if(!fh1TmpRho)fh1TmpRho = new TH1F("fh1TmpRho","tmp histo for jet shape",40,0.,2);
	fh2LTrackPtJetPt[ij][i] = new TH2F(Form("fh2LTrackPtJetPt%s_j%d",cAdd.Data(),i),
					   Form("pt of leadin track within a jet vs jet %s;p_{T,lead in jet};p_{T.jet};",
						cAdd.Data()),
					   200,0.,200.,nBinPt,binLimitsPt);
	fHistList->Add(fh2LTrackPtJetPt[ij][i]);
      }


      fh1DijetMinv[ij]                = new TH1F(Form("fh1Dijet%sMinv",cAdd.Data()),"Dijet invariant mass;m_{JJ}",nBinPt,binLimitsPt);
      fHistList->Add(fh1DijetMinv[ij]);

      fh2DijetDeltaPhiPt[ij]       = new TH2F(Form("fh2Dijet%sDeltaPhiPt",cAdd.Data()),"Difference in the azimuthal angle;#Delta#phi;p_{T,2};Entries",180,0.,TMath::Pi(),nBinPt,binLimitsPt);
      fHistList->Add(fh2DijetDeltaPhiPt[ij]);

      fh2DijetAsymPt[ij]            = new TH2F(Form("fh2Dijet%sAsym",cAdd.Data()),"Pt asymmetry;#Deltap_{T}/(p_{T,1}+p_{T,2});p_{T,1};Entries",50,0.,1.,nBinPt,binLimitsPt);
      fHistList->Add(fh2DijetAsymPt[ij]);

      fh2DijetPt2vsPt1[ij]          = new TH2F(Form("fh2Dijet%sPt2vsPt1",cAdd.Data()),"Pt2 versus Pt1;p_{T,1} (GeV/c);p_{T,2} (GeV/c)",250,0.,250.,250,0.,250.);
      fHistList->Add(fh2DijetPt2vsPt1[ij]);
      fh2DijetDifvsSum[ij]         = new TH2F(Form("fh2Dijet%sDifvsSum",cAdd.Data()),"Pt difference vs Pt sum;p_{T,1}+p_{T,2} (GeV/c);#Deltap_{T} (GeV/c)",400,0.,400.,150,0.,150.);
      fHistList->Add( fh2DijetDifvsSum[ij]);
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

void AliAnalysisTaskJetSpectrum2::Init()
{
  //
  // Initialization
  //

  if (fDebug > 1) printf("AnalysisTaskJetSpectrum2::Init() \n");

}

void AliAnalysisTaskJetSpectrum2::UserExec(Option_t */*option*/){


  Bool_t selected = kTRUE;

  if(fUseGlobalSelection&&fEventSelectionMask==0){
    selected = AliAnalysisHelperJetTasks::Selected();
  }
  else if(fUseGlobalSelection&&fEventSelectionMask>0){
    selected = AliAnalysisHelperJetTasks::TestSelectInfo(fEventSelectionMask);
  }

  if(fEventClass>0){
    selected = selected&&(AliAnalysisHelperJetTasks::EventClass()==fEventClass);
  }

  if(!selected){
    // no selection by the service task, we continue
    if (fDebug > 1)Printf("Not selected %s:%d SelectInfo %d  Class %d",(char*)__FILE__,__LINE__, AliAnalysisHelperJetTasks::Selected(),AliAnalysisHelperJetTasks::EventClass());
    PostData(1, fHistList);
    return;
  }


  static AliAODEvent* aod = 0;

  // take all other information from the aod we take the tracks from
  if(!aod){
   if(fUseAODTrackInput)aod = fAODIn;
   else aod = fAODOut;
  }


  //
  // Execute analysis for current event
  //
  if (fDebug > 1)printf("Analysing event # %5d\n", (Int_t) fEntry);  
  AliAODHandler *aodH = dynamic_cast<AliAODHandler*>(AliAnalysisManager::GetAnalysisManager()->GetOutputEventHandler());
  if(!aodH){
    Printf("%s:%d no output aodHandler found Jet",(char*)__FILE__,__LINE__);
    return;
  }

  if (fDebug > 10)Printf("%s:%d",(char*)__FILE__,__LINE__);
  TClonesArray *aodRecJets = 0;

  if(fAODOut&&!aodRecJets){
    aodRecJets = dynamic_cast<TClonesArray*>(fAODOut->FindListObject(fBranchRec.Data()));
  }
  if(fAODExtension&&!aodRecJets){
    aodRecJets = dynamic_cast<TClonesArray*>(fAODExtension->GetAOD()->FindListObject(fBranchRec.Data()));
  }
  if(fAODIn&&!aodRecJets){
    aodRecJets = dynamic_cast<TClonesArray*>(fAODIn->FindListObject(fBranchRec.Data()));
  }



  if(!aodRecJets){
    if(fDebug){

	Printf("%s:%d no reconstructed Jet array with name %s in AOD",(char*)__FILE__,__LINE__,fBranchRec.Data());
	if(fAODIn){
	  Printf("Input AOD >>>>");
	  fAODIn->Print();    
	}
	if(fAODExtension){
	  Printf("AOD Extension >>>>");
	  fAODExtension->Print();
	}
	if(fAODOut){
	  Printf("Output AOD >>>>");
	  fAODOut->Print();    
	}
    }
	return;
  }

  TClonesArray *aodGenJets = 0;
  if(fBranchGen.Length()>0){
    if(fAODOut&&!aodGenJets){
      aodGenJets = dynamic_cast<TClonesArray*>(fAODOut->FindListObject(fBranchGen.Data()));
    }
    if(fAODExtension&&!aodGenJets){
      aodGenJets = dynamic_cast<TClonesArray*>(fAODExtension->GetAOD()->FindListObject(fBranchGen.Data()));
    }
    if(fAODIn&&!aodGenJets){
      aodGenJets = dynamic_cast<TClonesArray*>(fAODIn->FindListObject(fBranchGen.Data()));
    }

    if(!aodGenJets){
      Printf("%s:%d no generated Jet array with name %s in AOD",(char*)__FILE__,__LINE__,fBranchGen.Data());
      return;
    }
  }

  TClonesArray *aodBackRecJets = 0;
  if(fBranchBkgRec.Length()>0){
    if(fAODOut&&!aodBackRecJets){
      aodBackRecJets = dynamic_cast<TClonesArray*>(fAODOut->FindListObject(fBranchBkgRec.Data()));
    }
    if(fAODExtension&&!aodBackRecJets){
      aodBackRecJets = dynamic_cast<TClonesArray*>(fAODExtension->GetAOD()->FindListObject(fBranchBkgRec.Data()));
    }
    if(fAODIn&&!aodBackRecJets){
      aodBackRecJets = dynamic_cast<TClonesArray*>(fAODIn->FindListObject(fBranchBkgRec.Data()));
    }

    if(!aodBackRecJets){
      Printf("%s:%d no background rec Jet array with name %s in AOD",(char*)__FILE__,__LINE__,fBranchBkgRec.Data());
      return;
    }
  }


  TClonesArray *aodBackGenJets = 0;

  if(fBranchBkgGen.Length()>0){
    if(fAODOut&&!aodBackGenJets){
      aodBackGenJets = dynamic_cast<TClonesArray*>(fAODOut->FindListObject(fBranchBkgGen.Data()));
    }
    if(fAODExtension&&!aodBackGenJets){
      aodBackGenJets = dynamic_cast<TClonesArray*>(fAODExtension->GetAOD()->FindListObject(fBranchBkgGen.Data()));
    }
    if(fAODIn&&!aodBackGenJets){
      aodBackGenJets = dynamic_cast<TClonesArray*>(fAODIn->FindListObject(fBranchBkgGen.Data()));
    }

    if(!aodBackGenJets){
      Printf("%s:%d no background rec Jet array with name %s in AOD",(char*)__FILE__,__LINE__,fBranchBkgGen.Data());
      return;
    }
  }

 
  // new Scheme
  // first fill all the pure  histograms, i.e. full jets 
  // in case of the correaltion limit to the n-leading jets

  // reconstructed

  
  // generated


  // second fill the correlation histograms, here we limit to the n-th leading jet in case of the reconstructed



  TList genJetsList;         // full acceptance
  TList genJetsListCut;      // acceptance cut
  TList recJetsList;         // full acceptance
  TList recJetsListCut;      // acceptance cut


  GetListOfJets(&recJetsList,aodRecJets,0);
  GetListOfJets(&recJetsListCut,aodRecJets,1);

  if(fBranchGen.Length()>0){
    GetListOfJets(&genJetsList,aodGenJets,0);
    GetListOfJets(&genJetsListCut,aodGenJets,1);
  }

  Double_t eventW = 1;
  Double_t ptHard = 0; 
  Double_t nTrials = 1; // Trials for MC trigger 
  fh1Trials->Fill("#sum{ntrials}",fAvgTrials); 

  // Getting some global properties
  fCentrality = GetCentrality();
  if(fCentrality<=0)fCentrality = 0;
  fh1Centrality->Fill(fCentrality);


  if((fAnalysisType&kAnaMCESD)==kAnaMCESD){
    // this is the part we only use when we have MC information
    AliMCEvent* mcEvent = MCEvent();
    //    AliStack *pStack = 0; 
    if(!mcEvent){
      Printf("%s:%d no mcEvent",(char*)__FILE__,__LINE__);
      return;
    }
    AliGenPythiaEventHeader*  pythiaGenHeader = AliAnalysisHelperJetTasks::GetPythiaEventHeader(mcEvent);
    if(pythiaGenHeader){
      nTrials = pythiaGenHeader->Trials();
      ptHard  = pythiaGenHeader->GetPtHard();
      int iProcessType = pythiaGenHeader->ProcessType();
      // 11 f+f -> f+f
      // 12 f+barf -> f+barf
      // 13 f+barf -> g+g
      // 28 f+g -> f+g
      // 53 g+g -> f+barf
      // 68 g+g -> g+g
      if (fDebug > 10)Printf("%d iProcessType %d",__LINE__, iProcessType);
      if(fDebug>20)AliAnalysisHelperJetTasks::PrintStack(mcEvent);
    } 
  }// (fAnalysisType&kMCESD)==kMCESD)


  // we simply fetch the tracks/mc particles as a list of AliVParticles

  TList recParticles;
  TList genParticles;

  Int_t nT = GetListOfTracks(&recParticles,fTrackTypeRec);
  if(fDebug>2)Printf("%s:%d Selected Rec tracks: %d %d",(char*)__FILE__,__LINE__,nT,recParticles.GetEntries());
  nT = GetListOfTracks(&genParticles,fTrackTypeGen);
  if(fDebug>2)Printf("%s:%d Selected Gen tracks: %d %d",(char*)__FILE__,__LINE__,nT,genParticles.GetEntries());

  //  CalculateReactionPlaneAngle(&recParticles);
  fRPAngle = 0;
  
  if(fRPMethod==0)fRPAngle = aod->GetHeader()->GetEventplane();
  else if(fRPMethod==1||fRPMethod==2){
    fRPAngle = aod->GetHeader()->GetQTheta(fRPMethod);
  }
  fh1RP->Fill(fRPAngle);
  fh2RPCentrality->Fill(fCentrality,fRPAngle);
  // Event control and counting ...  
  // MC
  fh1PtHard->Fill(ptHard,eventW);
  fh1PtHardNoW->Fill(ptHard,1);
  fh1PtHardTrials->Fill(ptHard,nTrials);

  // Real
  if(aod->GetPrimaryVertex()){// No vtx for pure MC
    fh1ZVtx->Fill(aod->GetPrimaryVertex()->GetZ());
  }


  Int_t recMult1 = recParticles.GetEntries();
  Int_t genMult1 = genParticles.GetEntries();

  Int_t recMult2 = MultFromJetRefs(aodBackRecJets);
  Int_t genMult2 = MultFromJetRefs(aodBackGenJets);

  fh2MultRec->Fill(recMult1,recMult2);
  fh2MultGen->Fill(genMult1,genMult2);
  fMultRec = recMult1;
  if(fMultRec<=0)fMultRec = recMult2;
  fMultGen = genMult1;
  if(fMultGen<=0)fMultGen = genMult2;

  Double_t var0[2] = {0,};
  var0[0] = fCentrality;
  var0[1] = fMultRec;
  fhnEvent->Fill(var0);

  // the loops for rec and gen should be indentical... pass it to a separate
  // function ...
  // Jet Loop
  // Track Loop
  // Jet Jet Loop
  // Jet Track loop

  FillJetHistos(recJetsListCut,recParticles,kJetRec);
  FillJetHistos(recJetsList,recParticles,kJetRecFull);
  FillTrackHistos(recParticles,kJetRec);

  FillJetHistos(genJetsListCut,genParticles,kJetGen);
  FillJetHistos(genJetsList,genParticles,kJetGenFull);
  FillTrackHistos(genParticles,kJetGen);

  // Here follows the jet matching part
  // delegate to a separated method?

  if(fDoMatching){
    FillMatchHistos(recJetsList,genJetsList);
  }

  if (fDebug > 10)Printf("%s:%d",(char*)__FILE__,__LINE__);
  PostData(1, fHistList);
}

void AliAnalysisTaskJetSpectrum2::FillJetHistos(TList &jetsList,TList &particlesList,Int_t iType){

  if(iType>=kJetTypes){
    return;
  }
  if(!fFlagJetType[iType])return;

  Int_t refMult = fMultRec;
  if(iType==kJetGen||iType==kJetGenFull){
    refMult = fMultGen;
  }

  Int_t nJets = jetsList.GetEntries(); 
  fh1NJets[iType]->Fill(nJets);

  if(nJets<=0)return;
  
  Float_t ptOld = 1.E+32;
  Float_t pT0 = 0;
  Float_t pT1 = 0;
  Float_t phi0 = 0;
  Float_t phi1 = 0;
  Int_t ij0 = -1;
  Int_t ij1 = -1;

  Double_t var1[6] = {0,}; // jet number;p_{T,jet};cent;# tracks;RP;area
  var1[2] = fCentrality; 
  var1[3] = refMult;

  Double_t var2[6] = {0,}; // jet number;p_{T,jet};cent;#eta;#phi;area
  var2[2] = fCentrality;

  for(int ij = 0;ij < nJets;ij++){
    AliAODJet *jet = (AliAODJet*)jetsList.At(ij);
    Float_t ptJet = jet->Pt();
    if(ptJet<0.150)ptJet = jet->GetPtSubtracted(0);
    if(jet->Trigger()&fJetTriggerExcludeMask){
      fh1PtJetsInRej[iType]->Fill(ptJet);
      continue;
    }
    fh1PtJetsIn[iType]->Fill(ptJet);
    if(ptJet>ptOld){
      Printf("%s:%d Jets Type %d Not Sorted !! %d:%.3E %d:%.3E",(char*)__FILE__,__LINE__,iType,ij,ptJet,ij-1,ptOld);
    }
    ptOld = ptJet;
    
    // find the dijets assume sorting and acceptance cut...
    if(ij==0){
      ij0 = ij;
      pT0 = ptJet;
      phi0 = jet->Phi();
      if(phi0<0)phi0+=TMath::Pi()*2.;
    }
    else if(ptJet>pT1){
      // take only the backward hemisphere??                                                        
      phi1 = jet->Phi();
      if(phi1<0)phi1+=TMath::Pi()*2.;
      Float_t dPhi = phi1 - phi0;
      if(dPhi>TMath::Pi())dPhi = dPhi - 2.*TMath::Pi();
      if(dPhi<(-1.*TMath::Pi()))dPhi = dPhi + 2.*TMath::Pi();
      if(TMath::Abs(TMath::Pi()-dPhi)<fDeltaPhiWindow){
	ij1 = ij;
	pT1 = ptJet;
      }
    }
    // fill jet histos for kmax jets

      Float_t phiJet = jet->Phi();
      Float_t etaJet = jet->Eta();
      if(phiJet<0)phiJet+=TMath::Pi()*2.;    
      fh1TmpRho->Reset();
      if(ij<kMaxJets)fh1PtIn[iType][ij]->Fill(ptJet);

      fh1PtIn[iType][kMaxJets]->Fill(ptJet);
      // fill leading jets...
      AliVParticle *leadTrack = LeadingTrackFromJetRefs(jet);
      //      AliVParticle *leadTrack = LeadingTrackInCone(jet,&particlesList);
      Int_t phiBin = GetPhiBin(phiJet-fRPAngle);

      var1[1] = ptJet;
      var1[4] = phiBin;
      var1[5] = jet->EffectiveAreaCharged();

      var2[1] = ptJet;
      var2[3] = etaJet;
      var2[4] = phiJet;
      var2[5] = jet->EffectiveAreaCharged();
      if(ij<kMaxJets){
	if(leadTrack)fh2LTrackPtJetPt[iType][ij]->Fill(leadTrack->Pt(),ptJet);
	var1[0] = ij;
	var2[0] = ij;
	fhnJetPt[iType]->Fill(var1);
	fhnJetPtQA[iType]->Fill(var2);
      }
      var1[0] = kMaxJets;// fill for all jets
      var2[0] = kMaxJets;// fill for all jets
      fhnJetPt[iType]->Fill(var1);
      fhnJetPtQA[iType]->Fill(var2);
      if(leadTrack)fh2LTrackPtJetPt[iType][kMaxJets]->Fill(leadTrack->Pt(),ptJet);

      if(particlesList.GetSize()&&ij<kMaxJets){
	// Particles... correlated with jets...
	for(int it = 0;it<particlesList.GetEntries();++it){
	  AliVParticle *part = (AliVParticle*)particlesList.At(it);
	  Float_t deltaR = jet->DeltaR(part);
	  if(ptJet>0)fh1TmpRho->Fill(deltaR,part->Pt()/ptJet);
	}
	// fill the jet shapes
      }// if we have particles
  }// Jet Loop


  // do something with dijets...
  if(ij0>=0&&ij1>0){
    AliAODJet *jet0 = (AliAODJet*)jetsList.At(ij0);
    Double_t ptJet0 = jet0->Pt();
    Double_t phiJet0 = jet0->Phi();
    if(phiJet0<0)phiJet0+=TMath::Pi()*2.;    	

    AliAODJet *jet1 = (AliAODJet*)jetsList.At(ij1);
    Double_t ptJet1 = jet1->Pt();
    Double_t phiJet1 = jet1->Phi();
    if(phiJet1<0)phiJet1+=TMath::Pi()*2.;    	

    Float_t deltaPhi = phiJet0 - phiJet1;
    if(deltaPhi>TMath::Pi())deltaPhi = deltaPhi - 2.*TMath::Pi();
    if(deltaPhi<(-1.*TMath::Pi()))deltaPhi = deltaPhi + 2.*TMath::Pi();      
    deltaPhi = TMath::Abs(deltaPhi);
    fh2DijetDeltaPhiPt[iType]->Fill(deltaPhi,ptJet1);      

    Float_t asym = 9999;
    if((ptJet0+ptJet1)>0)asym = (ptJet0-ptJet1)/(ptJet0+ptJet1);
      fh2DijetAsymPt[iType]->Fill(asym,ptJet0);
      fh2DijetPt2vsPt1[iType]->Fill(ptJet0,ptJet1);        
      fh2DijetDifvsSum[iType]->Fill(ptJet0+ptJet1,ptJet0-ptJet1);        
      Float_t minv = 2.*(jet0->P()*jet1->P()-
			 jet0->Px()*jet1->Px()- 
			 jet0->Py()*jet1->Py()- 
			 jet0->Pz()*jet1->Pz());    // assume mass == 0;
      if(minv<0)minv=0; // prevent numerical instabilities
      minv = TMath::Sqrt(minv);
      fh1DijetMinv[iType]->Fill(minv);            
  }
  


  // count down the jets above thrueshold
  Int_t nOver = nJets;
  if(nOver>0){
    TIterator *jetIter = jetsList.MakeIterator();
    AliAODJet *tmpJet = (AliAODJet*)(jetIter->Next());  
    if(tmpJet){
      Float_t pt = tmpJet->Pt();
      for(int i = 1;i <= fh2NJetsPt[iType]->GetNbinsX();i++){
	Float_t ptCut = fh2NJetsPt[iType]->GetXaxis()->GetBinCenter(i);
	while(pt<ptCut&&tmpJet){
	  nOver--;
	  tmpJet = (AliAODJet*)(jetIter->Next()); 
	  if(tmpJet){
	    pt = tmpJet->Pt();
	  }
	}
	if(nOver<=0)break;
	fh2NJetsPt[iType]->Fill(ptCut,nOver);
      }
    }
    delete jetIter;
  }
}

void AliAnalysisTaskJetSpectrum2::FillTrackHistos(TList &particlesList,int iType){

  if(fFlagJetType[iType]<=0)return;
  Int_t refMult = fMultRec;
  if(iType==kJetGen||iType==kJetGenFull){
    refMult = fMultGen;

  }

  // 
  Double_t var3[5]; // track number;p_{T};cent;#tracks;RP
  var3[2] = fCentrality;
  var3[3] = refMult;
  Double_t var4[5]; // track number;p_{T};cent;#eta;#phi
  var4[2] = fCentrality;
  Int_t nTrackOver = particlesList.GetSize();
  // do the same for tracks and jets
  if(nTrackOver>0){
    TIterator *trackIter = particlesList.MakeIterator();
    AliVParticle *tmpTrack = (AliVParticle*)(trackIter->Next());  
    Float_t pt = tmpTrack->Pt();
    for(int i = 1;i <= fh2NTracksPt[iType]->GetNbinsX();i++){
      Float_t ptCut = fh2NTracksPt[iType]->GetXaxis()->GetBinCenter(i);
      while(pt<ptCut&&tmpTrack){
	nTrackOver--;
	tmpTrack = (AliVParticle*)(trackIter->Next()); 
	if(tmpTrack){
	  pt = tmpTrack->Pt();
	}
      }
      if(nTrackOver<=0)break;
      fh2NTracksPt[iType]->Fill(ptCut,nTrackOver);
    }

    
    trackIter->Reset();
    AliVParticle *leading = (AliVParticle*)particlesList.At(0);
    Float_t sumPt = 0;

    while((tmpTrack = (AliVParticle*)(trackIter->Next()))){
      Float_t tmpPt = tmpTrack->Pt();
      fh1PtTracksIn[iType]->Fill(tmpPt);
      fh1PtTracksInLow[iType]->Fill(tmpPt);

      sumPt += tmpPt;
      Float_t tmpPhi = tmpTrack->Phi();
      if(tmpPhi<0)tmpPhi+=TMath::Pi()*2.;    

      
      Float_t phiRP = tmpPhi-fRPAngle;
      if(phiRP>TMath::Pi())phiRP -= TMath::Pi();
      if(phiRP<0)phiRP += TMath::Pi();
      if(phiRP<0)phiRP += TMath::Pi();
      const float allPhi = -1./180.*TMath::Pi();

      if(tmpPt<100){
	fp2MultRPPhiTrackPt[iType]->Fill(refMult,phiRP,tmpPt);
	fp2MultRPPhiTrackPt[iType]->Fill(refMult,allPhi,tmpPt);
	
	fp2CentRPPhiTrackPt[iType]->Fill(fCentrality,phiRP,tmpPt);
	fp2CentRPPhiTrackPt[iType]->Fill(fCentrality,allPhi,tmpPt);
      }
      Int_t phiBin = GetPhiBin(tmpPhi-fRPAngle);
      var3[0] = 1;
      var3[1] = tmpPt;
      var3[4] = phiBin;
      
      var4[0] = 1;
      var4[1] = tmpPt;
      var4[3] = tmpTrack->Eta();
      var4[4] = tmpPhi;


      fhnTrackPt[iType]->Fill(var3);
      fhnTrackPtQA[iType]->Fill(var4);

      if(tmpTrack==leading){
	var3[0] = 0;
	var4[0] = 0;
	fhnTrackPt[iType]->Fill(var3);
	fhnTrackPtQA[iType]->Fill(var4);
	continue;
      }
    }  
    fh1SumPtTrack[iType]->Fill(sumPt);
    delete trackIter;
  }

}


void AliAnalysisTaskJetSpectrum2::FillMatchHistos(TList &recJetsList,TList &genJetsList){


  // Fill al the matching histos
  // It is important that the acceptances for the mathing are as large as possible
  // to avoid false matches on the edge of acceptance
  // therefore we add some extra matching jets as overhead

  static TArrayI aGenIndex(recJetsList.GetEntries());
  if(aGenIndex.GetSize()<recJetsList.GetEntries())aGenIndex.Set(recJetsList.GetEntries());

  static TArrayI aRecIndex(genJetsList.GetEntries());
  if(aRecIndex.GetSize()<genJetsList.GetEntries())aRecIndex.Set(genJetsList.GetEntries());

  if(fDebug){
    Printf("New Gens List %d rec index Array %d",genJetsList.GetEntries(),aRecIndex.GetSize());
    Printf("New Rec List %d gen indey Array %d",recJetsList.GetEntries(),aGenIndex.GetSize());
  }
  AliAnalysisHelperJetTasks::GetClosestJets(&genJetsList,TMath::Min((Int_t)fNMatchJets,(Int_t)genJetsList.GetEntries()),
					    &recJetsList,TMath::Min((Int_t)fNMatchJets,(Int_t)recJetsList.GetEntries()),
					    aGenIndex,aRecIndex,fDebug);

  if(fDebug){
    for(int i = 0;i< aGenIndex.GetSize();++i){ 
      if(aGenIndex[i]>=0)Printf("iGenFound: %d -> %d",i,aGenIndex[i]); 
    }
    for(int i = 0;i< aRecIndex.GetSize();++i){
      if(aRecIndex[i]>=0)Printf("iRecFound: %d -> %d",i,aRecIndex[i]); 
    }
  }

  Double_t container[6];

  // loop over generated jets
  // consider the 
  for(int ig = 0;ig < genJetsList.GetEntries();++ig){
    AliAODJet *genJet = (AliAODJet*)genJetsList.At(ig);
    Double_t ptGen = genJet->Pt();
    Double_t phiGen = genJet->Phi();
    if(phiGen<0)phiGen+=TMath::Pi()*2.;    
    Double_t etaGen = genJet->Eta();
    container[3] = ptGen;
    container[4] = etaGen;
    container[5] = phiGen;
    fhnJetContainer->Fill(&container[3],kStep0);
    if(JetSelected(genJet)){
      fhnJetContainer->Fill(&container[3],kStep1);
      Int_t ir = aRecIndex[ig];
      if(ir>=0&&ir<recJetsList.GetEntries()){   
	fhnJetContainer->Fill(&container[3],kStep2);
	AliAODJet* recJet = (AliAODJet*)recJetsList.At(ir); 
	if(JetSelected(recJet))fhnJetContainer->Fill(&container[3],kStep4);
	if(JetSelected(recJet))fhnJetContainer->Fill(&container[3],kStep3);
      }
    }
  }// loop over generated jets used for matching...



  // loop over reconstructed jets
  for(int ir = 0;ir < recJetsList.GetEntries();++ir){
    AliAODJet *recJet = (AliAODJet*)recJetsList.At(ir);
    Double_t etaRec = recJet->Eta();
    Double_t ptRec = recJet->Pt();
    Double_t phiRec = recJet->Phi();
    if(phiRec<0)phiRec+=TMath::Pi()*2.;    
    // do something with dijets...
    
    container[0] = ptRec;
    container[1] = etaRec;
    container[2] = phiRec;

    fhnJetContainer->Fill(container,kStep0+kMaxStep);
    if (fDebug > 10)Printf("%s:%d",(char*)__FILE__,__LINE__);
  
    if(JetSelected(recJet)){
      fhnJetContainer->Fill(container,kStep1+kMaxStep);
      // Fill Correlation
      Int_t ig = aGenIndex[ir];
      if(ig>=0 && ig<genJetsList.GetEntries()){
	fhnJetContainer->Fill(container,kStep2+kMaxStep);
	if (fDebug > 10)Printf("%s:%d ig = %d ir = %d",(char*)__FILE__,__LINE__,ig,ir);
	AliAODJet *genJet = (AliAODJet*)genJetsList.At(ig);
	Double_t ptGen  = genJet->Pt();
	Double_t phiGen = genJet->Phi();
	if(phiGen<0)phiGen+=TMath::Pi()*2.; 
	Double_t etaGen = genJet->Eta();
      
	container[3] = ptGen;
	container[4] = etaGen;
	container[5] = phiGen;
	// 
	// we accept only jets which are detected within a smaller window, to avoid ambigious pair association at the edges of the acceptance
	// 
	if(JetSelected(genJet))fhnJetContainer->Fill(container,kStep4+kMaxStep);
	fhnJetContainer->Fill(container,kStep3+kMaxStep);
	fhnCorrelation->Fill(container,0);
	if(ptGen>0){
	  Float_t delta = (ptRec-ptGen)/ptGen;
	  fh2RelPtFGen->Fill(ptGen,delta);
	  fh2PtFGen->Fill(ptGen,ptRec);
	}
      } 
    }// loop over reconstructed jets
  }
  if (fDebug > 10)Printf("%s:%d",(char*)__FILE__,__LINE__);
}


void AliAnalysisTaskJetSpectrum2::MakeJetContainer(){
  //
  // Create the particle container for the correction framework manager and 
  // link it
  //
  const Int_t kNvar   = 3 ; //number of variables on the grid:pt,eta, phi
  const Double_t kPtmin = 0.0, kPtmax = 250.; // we do not want to have empty bins at the beginning...
  const Double_t kEtamin = -3.0, kEtamax = 3.0;
  const Double_t kPhimin = 0., kPhimax = 2. * TMath::Pi();

  // can we neglect migration in eta and phi?
  // phi should be no problem since we cover full phi and are phi symmetric
  // eta migration is more difficult  due to needed acceptance correction
  // in limited eta range

  //arrays for the number of bins in each dimension
  Int_t iBin[kNvar];
  iBin[0] = 125; //bins in pt
  iBin[1] =  1; //bins in eta 
  iBin[2] = 1; // bins in phi

  //arrays for lower bounds :
  Double_t* binEdges[kNvar];
  for(Int_t ivar = 0; ivar < kNvar; ivar++)
    binEdges[ivar] = new Double_t[iBin[ivar] + 1];

  //values for bin lower bounds
  //  for(Int_t i=0; i<=iBin[0]; i++) binEdges[0][i]=(Double_t)TMath::Power(10,TMath::Log10(kPtmin) + (TMath::Log10(kPtmax)-TMath::Log10(kPtmin))/iBin[0]*(Double_t)i);  
  for(Int_t i=0; i<=iBin[0]; i++) binEdges[0][i]=(Double_t)kPtmin  + (kPtmax-kPtmin)/(Double_t)iBin[0]*(Double_t)i;
  for(Int_t i=0; i<=iBin[1]; i++) binEdges[1][i]=(Double_t)kEtamin  + (kEtamax-kEtamin)/iBin[1]*(Double_t)i;
  for(Int_t i=0; i<=iBin[2]; i++) binEdges[2][i]=(Double_t)kPhimin  + (kPhimax-kPhimin)/iBin[2]*(Double_t)i;


  fhnJetContainer = new AliTHn(Form("fahnJetContainer"),Form("AliTHn jet info"),kMaxStep*2,kNvar,iBin);
  for (int k=0; k<kNvar; k++) {
    fhnJetContainer->SetBinLimits(k,binEdges[k]);
  }

  //create correlation matrix for unfolding
  Int_t thnDim[2*kNvar];
  for (int k=0; k<kNvar; k++) {
    //first half  : reconstructed 
    //second half : MC
    thnDim[k]      = iBin[k];
    thnDim[k+kNvar] = iBin[k];
  }

  fhnCorrelation = new AliTHn("fahnCorrelation","AliTHn with correlations",1,2*kNvar,thnDim);
  for (int k=0; k<kNvar; k++) {
    fhnCorrelation->SetBinLimits(k,binEdges[k]);
    fhnCorrelation->SetBinLimits(k+kNvar,binEdges[k]);
  }

  for(Int_t ivar = 0; ivar < kNvar; ivar++)
    delete [] binEdges[ivar];


}

void AliAnalysisTaskJetSpectrum2::Terminate(Option_t */*option*/)
{
// Terminate analysis
//
    if (fDebug > 1) printf("AnalysisJetSpectrum2: Terminate() \n");
}


Int_t  AliAnalysisTaskJetSpectrum2::GetListOfTracks(TList *list,Int_t type){

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
      if(TMath::Abs(tr->Eta())>fTrackRecEtaWindow)continue;
      if(tr->Pt()<fMinTrackPt)continue;
      if(fDebug>0){
	if(tr->Pt()>20){
	  Printf("High pT track found in Event %d with p_T, %E",(int)Entry(),tr->Pt());
	  Printf("%s read event, %d",fInputHandler->GetTree()->GetCurrentFile()->GetName(),(Int_t)fInputHandler->GetTree()->GetReadEntry());
	  tr->Print();
	  //	tr->Dump();
	  AliESDEvent *fESD = dynamic_cast<AliESDEvent*> (InputEvent());
	  if(fESD){
	    AliESDtrack *esdTr = (AliESDtrack*)fESD->GetTrack(TMath::Abs(tr->GetID()+1));
	    esdTr->Print("");
	    // esdTr->Dump();
	  }
	}
      }
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
      if(part->Pt()<fMinTrackPt)continue;
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
      if(!part)continue;
      if(part->Pt()<fMinTrackPt)continue;
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
	  if(TMath::Abs(part->Eta())>fTrackRecEtaWindow)continue;
	  list->Add(part);
	}
	iCount++;
      }
    }
  }// AODMCparticle
  list->Sort();
  return iCount;

}


Float_t AliAnalysisTaskJetSpectrum2::GetCentrality(){
    AliAODEvent *aod = 0;
    if(fUseAODTrackInput)aod = dynamic_cast<AliAODEvent*>(InputEvent());
    else aod = AODEvent();
    if(!aod){
      return 101;
    }
    return aod->GetHeader()->GetCentrality();
}



Bool_t  AliAnalysisTaskJetSpectrum2::JetSelected(AliAODJet *jet){
  Bool_t selected = false;
  
  if(!jet)return selected;

  if(fabs(jet->Eta())<fJetRecEtaWindow&&jet->Pt()>fMinJetPt){
    selected = kTRUE;
  }
  return selected;

}

Int_t  AliAnalysisTaskJetSpectrum2::GetListOfJets(TList *list,TClonesArray* jarray,Int_t type){

  if(fDebug>2)Printf("%s:%d Selecting jets with cuts %d",(char*)__FILE__,__LINE__,type);
  Int_t iCount = 0;

  if(!jarray){
    Printf("%s:%d no Jet array",(char*)__FILE__,__LINE__);
    return iCount;
  }


  for(int ij=0;ij<jarray->GetEntries();ij++){
    AliAODJet* jet = (AliAODJet*)jarray->At(ij);
    if(!jet)continue;
    if(type==0){
      // No cut at all, main purpose here is sorting      
      list->Add(jet);
      iCount++;
    }
    else if(type == 1){
      // eta cut
      if(JetSelected(jet)){
	list->Add(jet);
	iCount++;
      }
    }
  }

  list->Sort();
  return iCount;

}


Int_t AliAnalysisTaskJetSpectrum2::MultFromJetRefs(TClonesArray* jets){
  if(!jets)return 0;

  Int_t refMult = 0;
  for(int ij = 0;ij < jets->GetEntries();++ij){
    AliAODJet* jet = (AliAODJet*)jets->At(ij);
    if(!jet)continue;
    TRefArray *refs = jet->GetRefTracks();
    if(!refs)continue;
    refMult += refs->GetEntries();
  }
  return refMult;

}


AliVParticle *AliAnalysisTaskJetSpectrum2::LeadingTrackFromJetRefs(AliAODJet* jet){
  if(!jet)return 0;
  TRefArray *refs = jet->GetRefTracks();
  if(!refs) return 0;
  AliVParticle *leading = 0;
  Float_t fMaxPt = 0;
  for(int i = 0;i<refs->GetEntries();i++){
    AliVParticle *tmp = dynamic_cast<AliVParticle*>(refs->At(i));
    if(!tmp)continue;
    if(tmp->Pt()>fMaxPt){
      leading = tmp;
      fMaxPt = tmp->Pt();
    }
  }
  return leading;
}


AliVParticle *AliAnalysisTaskJetSpectrum2::LeadingTrackInCone(AliAODJet* jet,TList *list,Float_t r){
  if(!jet)return 0;
  if(!list) return 0;
  AliVParticle *leading = 0;
  Float_t fMaxPt = 0;
  for(int i = 0;i<list->GetEntries();i++){
    AliVParticle *tmp = (AliVParticle*)(list->At(i));
    if(!tmp)continue;
    if(jet->DeltaR(tmp)>r)continue;
    if(tmp->Pt()>fMaxPt){
      leading = tmp;
      fMaxPt = tmp->Pt();
    }
  }
  return leading;
}

Int_t AliAnalysisTaskJetSpectrum2::GetPhiBin(Double_t phi)
{
    Int_t phibin=-1;
    if(!(TMath::Abs(phi)<=2*TMath::Pi())){AliError("phi w.r.t. RP out of defined range");return -1;}
    Double_t phiwrtrp=TMath::ACos(TMath::Abs(TMath::Cos(phi)));
    phibin=Int_t(fNRPBins*phiwrtrp/(0.5*TMath::Pi()));
    if(phibin<0||phibin>=fNRPBins){AliError("Phi Bin not defined");}
    return phibin;
}

