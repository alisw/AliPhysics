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

#include "AliAnalysisTaskJetSpectrum2.h"
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


#include "AliAnalysisHelperJetTasks.h"

ClassImp(AliAnalysisTaskJetSpectrum2)

AliAnalysisTaskJetSpectrum2::AliAnalysisTaskJetSpectrum2(): AliAnalysisTaskSE(),
							    fJetHeaderRec(0x0),
							    fJetHeaderGen(0x0),
							    fAOD(0x0),
							    fhnCorrelation(0x0),
  fhnCorrelationPhiZRec(0x0),
  f1PtScale(0x0),
							    fBranchRec("jets"),
							    fBranchGen(""),
							    fUseAODJetInput(kFALSE),
							    fUseAODTrackInput(kFALSE),
							    fUseAODMCInput(kFALSE),
							    fUseGlobalSelection(kFALSE),
							    fUseExternalWeightOnly(kFALSE),
							    fLimitGenJetEta(kFALSE),
							    fFilterMask(0),
  fEventSelectionMask(0),
							    fAnalysisType(0),
  fTrackTypeRec(kTrackUndef),
  fTrackTypeGen(kTrackUndef),
  fAvgTrials(1),
  fExternalWeight(1),    
							    fRecEtaWindow(0.5),
  fMinJetPt(0),
							    fDeltaPhiWindow(20./180.*TMath::Pi()),
  fh1Xsec(0x0),
  fh1Trials(0x0),
  fh1PtHard(0x0),
  fh1PtHardNoW(0x0),  
  fh1PtHardTrials(0x0),
  fh1NGenJets(0x0),
  fh1NRecJets(0x0),
  fh1PtTrackRec(0x0),   
  fh1SumPtTrackRec(0x0),  
  fh1SumPtTrackAreaRec(0x0),  
  fh1TmpRho(0x0),
  fh1PtJetsRecIn(0x0),
  fh1PtJetsLeadingRecIn(0x0),
  fh1PtTracksRecIn(0x0),
  fh1PtTracksLeadingRecIn(0x0),
  fh1PtTracksGenIn(0x0),
  fh2NRecJetsPt(0x0),
  fh2NRecTracksPt(0x0),
  fh2JetsLeadingPhiEta(0x0),
  fh2JetsLeadingPhiPt(0x0),
  fh2TracksLeadingPhiEta(0x0),
  fh2TracksLeadingPhiPt(0x0),
  fh2TracksLeadingJetPhiPt(0x0),
  fh2JetPtJetPhi(0x0),
  fh2TrackPtTrackPhi(0x0),
  fh2DijetDeltaPhiPt(0x0),      
  fh2DijetAsymPt(0x0),          
  fh2DijetAsymPtCut(0x0),       
  fh2DijetDeltaPhiDeltaEta(0x0),
  fh2DijetPt2vsPt1(0x0),        
  fh2DijetDifvsSum(0x0),        
  fh1DijetMinv(0x0),            
  fh1DijetMinvCut(0x0),         
  fHistList(0x0)  
{
  for(int i = 0;i < kMaxStep*2;++i){
    fhnJetContainer[i] = 0;
  }
  for(int i = 0;i < kMaxJets;++i){
    fh1PtRecIn[i] = fh1PtGenIn[i] = 0;
    
    fh2PhiPt[i] = 0;
    fh2PhiEta[i] = 0; 
    fh2RhoPtRec[i] = 0; 
    fh2RhoPtGen[i] = 0; 
    fh2PsiPtGen[i] = 0; 
    fh2PsiPtRec[i] = 0;
    fh2FragRec[i] = 0;
    fh2FragLnRec[i] = 0;
    fh2FragGen[i] = 0;
    fh2FragLnGen[i] = 0;
  }  

}

AliAnalysisTaskJetSpectrum2::AliAnalysisTaskJetSpectrum2(const char* name):
  AliAnalysisTaskSE(name),
  fJetHeaderRec(0x0),
  fJetHeaderGen(0x0),
  fAOD(0x0),
  fhnCorrelation(0x0),
  fhnCorrelationPhiZRec(0x0),
  f1PtScale(0x0),
  fBranchRec("jets"),
  fBranchGen(""),
  fUseAODJetInput(kFALSE),
  fUseAODTrackInput(kFALSE),
  fUseAODMCInput(kFALSE),
  fUseGlobalSelection(kFALSE),
  fUseExternalWeightOnly(kFALSE),
  fLimitGenJetEta(kFALSE),
  fFilterMask(0),
  fEventSelectionMask(0),
  fAnalysisType(0),
  fTrackTypeRec(kTrackUndef),
  fTrackTypeGen(kTrackUndef),
  fAvgTrials(1),
  fExternalWeight(1),    
  fRecEtaWindow(0.5),
  fMinJetPt(0),
  fDeltaPhiWindow(20./180.*TMath::Pi()),
  fh1Xsec(0x0),
  fh1Trials(0x0),
  fh1PtHard(0x0),
  fh1PtHardNoW(0x0),  
  fh1PtHardTrials(0x0),
  fh1NGenJets(0x0),
  fh1NRecJets(0x0),
  fh1PtTrackRec(0x0),   
  fh1SumPtTrackRec(0x0),  
  fh1SumPtTrackAreaRec(0x0),  
  fh1TmpRho(0x0),
  fh1PtJetsRecIn(0x0),
  fh1PtJetsLeadingRecIn(0x0),
  fh1PtTracksRecIn(0x0),
  fh1PtTracksLeadingRecIn(0x0),
  fh1PtTracksGenIn(0x0),
  fh2NRecJetsPt(0x0),
  fh2NRecTracksPt(0x0),
  fh2JetsLeadingPhiEta(0x0),
  fh2JetsLeadingPhiPt(0x0),
  fh2TracksLeadingPhiEta(0x0),
  fh2TracksLeadingPhiPt(0x0),
  fh2TracksLeadingJetPhiPt(0x0),
  fh2JetPtJetPhi(0x0),
  fh2TrackPtTrackPhi(0x0),
  fh2DijetDeltaPhiPt(0x0),      
  fh2DijetAsymPt(0x0),          
  fh2DijetAsymPtCut(0x0),       
  fh2DijetDeltaPhiDeltaEta(0x0),
  fh2DijetPt2vsPt1(0x0),        
  fh2DijetDifvsSum(0x0),        
  fh1DijetMinv(0x0),            
  fh1DijetMinvCut(0x0),         
  fHistList(0x0)
{

  for(int i = 0;i < kMaxStep*2;++i){
    fhnJetContainer[i] = 0;
  }  
  for(int i = 0;i < kMaxJets;++i){
    fh1PtRecIn[i] = fh1PtGenIn[i] = 0;
    
    fh2PhiPt[i] = 0;
    fh2PhiEta[i] = 0; 
    fh2RhoPtRec[i] = 0; 
    fh2RhoPtGen[i] = 0; 
    fh2PsiPtGen[i] = 0; 
    fh2PsiPtRec[i] = 0;
    fh2FragRec[i] = 0;
    fh2FragLnRec[i] = 0;
    fh2FragGen[i] = 0;
    fh2FragLnGen[i] = 0;
  }
  DefineOutput(1, TList::Class());  
}



Bool_t AliAnalysisTaskJetSpectrum2::Notify()
{
  //
  // Implemented Notify() to read the cross sections
  // and number of trials from pyxsec.root
  // 

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
  return kTRUE;
}

void AliAnalysisTaskJetSpectrum2::UserCreateOutputObjects()
{

  //
  // Create the output container
  //


  // Connect the AOD


  MakeJetContainer();


  if (fDebug > 1) printf("AnalysisTaskJetSpectrum2::UserCreateOutputObjects() \n");

  OpenFile(1);
  if(!fHistList)fHistList = new TList();
  fHistList->SetOwner(kTRUE);
  Bool_t oldStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);

  //
  //  Histogram
    
  const Int_t nBinPt = 240;
  Double_t binLimitsPt[nBinPt+1];
  for(Int_t iPt = 0;iPt <= nBinPt;iPt++){
    if(iPt == 0){
      binLimitsPt[iPt] = 0.0;
    }
    else {// 1.0
      binLimitsPt[iPt] =  binLimitsPt[iPt-1] + 1.0;
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


  const Int_t nBinPhi2 = 360;
  Double_t binLimitsPhi2[nBinPhi2+1];
  for(Int_t iPhi2 = 0;iPhi2<=nBinPhi2;iPhi2++){
    if(iPhi2==0){
      binLimitsPhi2[iPhi2] = 0.;
    }
    else{
      binLimitsPhi2[iPhi2] = binLimitsPhi2[iPhi2-1] + 1/(Float_t)nBinPhi2 * TMath::Pi()*2;
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

  const Int_t nBinFrag = 25;


  fh1Xsec = new TProfile("fh1Xsec","xsec from pyxsec.root",1,0,1);
  fh1Xsec->GetXaxis()->SetBinLabel(1,"<#sigma>");

  fh1Trials = new TH1F("fh1Trials","trials root file",1,0,1);
  fh1Trials->GetXaxis()->SetBinLabel(1,"#sum{ntrials}");

  fh1PtHard = new TH1F("fh1PtHard","PYTHIA Pt hard;p_{T,hard}",nBinPt,binLimitsPt);
  fh1PtHardNoW = new TH1F("fh1PtHardNoW","PYTHIA Pt hard no weight;p_{T,hard}",nBinPt,binLimitsPt);
  fh1PtHardTrials = new TH1F("fh1PtHardTrials","PYTHIA Pt hard weight with trials;p_{T,hard}",nBinPt,binLimitsPt);

  fh1NGenJets  = new TH1F("fh1NGenJets","N generated jets",20,-0.5,19.5);
  fh1NRecJets = new TH1F("fh1NRecJets","N reconstructed jets",20,-0.5,19.5);

  fh1PtTrackRec = new TH1F("fh1PtTrackRec","Rec track P_T #eta < 0.9;p_{T} (GeV/c)",nBinPt,binLimitsPt);
  fh1SumPtTrackRec = new TH1F("fh1SumPtTrackRec","Sum Rec track P_T #eta <0.9;p_{T,sum} (GeV/c)",nBinPt,binLimitsPt);
  fh1SumPtTrackAreaRec = new TH1F("fh1SumPtTrackAreaRec","Sum Rec track P_T #eta <0.9 / 1.8 * 2 * 0.4*0.4;p_{T,sum} (GeV/c)",nBinPt,binLimitsPt);
  
  fh1PtJetsRecIn  = new TH1F("fh1PtJetsRecIn","Rec jets P_T;p_{T} (GeV/c)",nBinPt,binLimitsPt);
  fh1PtJetsLeadingRecIn = new TH1F("fh1PtJetsLeadingRecIn","Rec jets P_T;p_{T} (GeV/c)",nBinPt,binLimitsPt);
  fh1PtTracksRecIn  = new TH1F("fh1PtTracksRecIn","Rec tracks P_T #eta < 0.9;p_{T} (GeV/c)",nBinPt,binLimitsPt);
  fh1PtTracksLeadingRecIn  = new TH1F("fh1PtTracksLeadingRecIn","Rec tracks P_T #eta < 0.9;p_{T} (GeV/c)",nBinPt,binLimitsPt);
  fh1PtTracksGenIn  = new TH1F("fh1PtTracksGenIn","gen tracks P_T #eta < 0.9;p_{T} (GeV/c)",nBinPt,binLimitsPt);

  fh2NRecJetsPt = new TH2F("fh2NRecJetsPt","Number of jets above threshhold;p_{T,cut} (GeV/c);N_{jets}",nBinPt,binLimitsPt,50,-0.5,49.5);
  fh2NRecTracksPt = new TH2F("fh2NRecTracksPt","Number of tracks above threshhold;p_{T,cut} (GeV/c);N_{tracks}",nBinPt,binLimitsPt,50,-0.5,49.5);
  // 

  fh2JetsLeadingPhiEta = new TH2F("fh2JetsLeadingPhiEta","delta eta vs delta phi to leading jet;#Delta#phi;#Delta#eta",
				nBinPhi,binLimitsPhi,nBinEta,binLimitsEta);
  fh2JetsLeadingPhiPt = new TH2F("fh2JetsLeadingPhiPt","leading p_T vs delta phi to leading jet;#Delta#phi;p_{T} (GeV/c)",
				nBinPhi,binLimitsPhi,nBinPt,binLimitsPt);
  fh2TracksLeadingPhiEta = new TH2F("fh2TracksLeadingPhiEta","delta eta vs delta phi to leading track;#Delta#phi;#Delta#eta",
				    nBinPhi,binLimitsPhi,nBinEta,binLimitsEta);
  fh2TracksLeadingPhiPt = new TH2F("fh2TracksLeadingPhiPt","leading p_T vs delta phi to leading jet;#Delta#phi;p_{T} (GeV/c)",
				 nBinPhi,binLimitsPhi,nBinPt,binLimitsPt);

  fh2JetPtJetPhi = new TH2F("fh2JetPtJetPhi","Reconstructed jet phi vs. pt",nBinPt,binLimitsPt,nBinPhi2,binLimitsPhi2);
  fh2TrackPtTrackPhi = new TH2F("fh2TrackPtTrackPhi","Reconstructed track phi vs. pt",nBinPt,binLimitsPt,nBinPhi2,binLimitsPhi2);


  for(int ij = 0;ij < kMaxJets;++ij){
    fh1PtRecIn[ij] = new TH1F(Form("fh1PtRecIn_j%d",ij),"rec p_T input ;p_{T,rec}",nBinPt,binLimitsPt);
    fh1PtGenIn[ij] = new TH1F(Form("fh1PtGenIn_j%d",ij),"found p_T input ;p_{T,gen}",nBinPt,binLimitsPt);

    fh2PhiPt[ij] =  new TH2F(Form("fh2PhiPtRec_j%d",ij),"Jet pt vs delta phi;#Delta#phi;p_{T,jet}",
			     nBinPhi,binLimitsPhi,nBinPt,binLimitsPt);

    fh2PhiEta[ij] =  new TH2F(Form("fh2PhiEtaRec_j%d",ij),"delta eta vs delta phi for jets;#Delta#phi;#Delta#eta",
			      nBinPhi,binLimitsPhi,nBinEta,binLimitsEta);

    fh2RhoPtRec[ij] =  new TH2F(Form("fh2RhoPtRec_j%d",ij),"jet shape rho for jets;r;p_{T,rec};",
				20,0.,1.,nBinPt,binLimitsPt);
    fh2PsiPtRec[ij] =  new TH2F(Form("fh2PsiPtRec_j%d",ij),"jet shape psi for jets;r;p_{T,rec};",
				20,0.,1.,nBinPt,binLimitsPt);

    fh2RhoPtGen[ij] =  new TH2F(Form("fh2RhoPtGen_j%d",ij),"jet shape rho for jets;r;p_{T,gen};",
				20,0.,1.,nBinPt,binLimitsPt);
    fh2PsiPtGen[ij] =  new TH2F(Form("fh2PsiPtGen_j%d",ij),"jet shape psi for jets;r;p_{T,gen};",
				20,0.,1.,nBinPt,binLimitsPt);
    if(!fh1TmpRho)fh1TmpRho = new TH1F("fh1TmpRho","tmp histo for jet shape",20,0.,1);


    fh2FragRec[ij] = new TH2F(Form("fh2FragRec_j%d",ij),"Jet Fragmentation;x=p_{T,i}/p_{T,jet};p_{T,jet};1/N_{jet}dN_{ch}/dx",
			   nBinFrag,0.,1.,nBinPt,binLimitsPt);
    fh2FragLnRec[ij] = new TH2F(Form("fh2FragLnRec_j%d",ij),"Jet Fragmentation Ln;#xi=ln(p_{T,jet}/p_{T,i});p_{T,jet}(GeV);1/N_{jet}dN_{ch}/d#xi",
			     nBinFrag,0.,10.,nBinPt,binLimitsPt);

    fh2FragGen[ij] = new TH2F(Form("fh2FragGen_j%d",ij),"Jet Fragmentation;x=p_{T,i}/p_{T,jet};p_{T,jet};1/N_{jet}dN_{ch}/dx",
			      nBinFrag,0.,1.,nBinPt,binLimitsPt);
    fh2FragLnGen[ij] = new TH2F(Form("fh2FragLnGen_j%d",ij),"Jet Fragmentation Ln;#xi=ln(p_{T,jet}/p_{T,i});p_{T,jet}(GeV);1/N_{jet}dN_{ch}/d#xi",
				nBinFrag,0.,10.,nBinPt,binLimitsPt);
  }

  // Dijet histograms

  fh2DijetDeltaPhiPt       = new TH2F("fh2DeltaPhiPt","Difference in the azimuthal angle;#Delta#phi;p_{T,1};Entries",180,0.,TMath::Pi(),nBinPt,binLimitsPt);
  fh2DijetAsymPt            = new TH2F("fh2DijetAsym","Pt asymmetry;#Deltap_{T}/(p_{T,1}+p_{T,2});p_{T,1};Entries",50,0.,1.,nBinPt,binLimitsPt);
  fh2DijetAsymPtCut         = new TH2F("fh2DijetAsymCut","Pt asymmetry after delta phi cut;#Deltap_{T}/(p_{T,1}+p_{T,2});p_{T,1};Entries",50,0.,1.,nBinPt,binLimitsPt);
  fh2DijetDeltaPhiDeltaEta = new TH2F("fh2DijetDeltaPhiDeltaEta","Difference in the azimuthal angle;#Delta#phi;Entries",180,0.,TMath::Pi(),20,-2.,2.);
  fh2DijetPt2vsPt1          = new TH2F("fh2DijetPt2vsPt1","Pt2 versus Pt1;p_{T,1} (GeV/c);p_{T,2} (GeV/c)",250,0.,250.,250,0.,250.);
  fh2DijetDifvsSum          = new TH2F("fh2DijetDifvsSum","Pt difference vs Pt sum;p_{T,1}+p_{T,2} (GeV/c);#Deltap_{T} (GeV/c)",400,0.,400.,150,0.,150.);
  fh1DijetMinv               = new TH1F("fh1DijetMinv","Dijet invariant mass;m_{JJ}",nBinPt,binLimitsPt);
  fh1DijetMinvCut           = new TH1F("fh1DijetMinvCut","Dijet invariant mass;m_{JJ}",nBinPt,binLimitsPt);




  const Int_t saveLevel = 3; // large save level more histos
  if(saveLevel>0){
    fHistList->Add(fh1Xsec);
    fHistList->Add(fh1Trials);
    fHistList->Add(fh1PtHard);
    fHistList->Add(fh1PtHardNoW);
    fHistList->Add(fh1PtHardTrials);
    if(fBranchGen.Length()>0){
      fHistList->Add(fh1NGenJets);
      fHistList->Add(fh1PtTracksGenIn);
    }
    fHistList->Add(fh1PtJetsRecIn);
    fHistList->Add(fh1PtJetsLeadingRecIn);
    fHistList->Add(fh1PtTracksRecIn);
    fHistList->Add(fh1PtTracksLeadingRecIn);
    fHistList->Add(fh1NRecJets);
    fHistList->Add(fh1PtTrackRec);
    fHistList->Add(fh1SumPtTrackRec);
    fHistList->Add(fh1SumPtTrackAreaRec);
    fHistList->Add(fh2NRecJetsPt);
    fHistList->Add(fh2NRecTracksPt);
    fHistList->Add(fh2JetsLeadingPhiEta );
    fHistList->Add(fh2JetsLeadingPhiPt );
    fHistList->Add(fh2TracksLeadingPhiEta);
    fHistList->Add(fh2TracksLeadingPhiPt);
    for(int i = 0;i<kMaxStep*2;++i)fHistList->Add(fhnJetContainer[i]);
    for(int ij = 0;ij<kMaxJets;++ij){
      fHistList->Add( fh1PtRecIn[ij]);

      if(fBranchGen.Length()>0){	
	fHistList->Add(fh1PtGenIn[ij]);
	fHistList->Add(fh2FragGen[ij]);
	fHistList->Add(fh2FragLnGen[ij]);
	fHistList->Add(fh2RhoPtGen[ij]);
	fHistList->Add(fh2PsiPtGen[ij]);
      }
      fHistList->Add( fh2PhiPt[ij]);
      fHistList->Add( fh2PhiEta[ij]);
      fHistList->Add( fh2RhoPtRec[ij]);
      fHistList->Add( fh2PsiPtRec[ij]);
      fHistList->Add( fh2FragRec[ij]);
      fHistList->Add( fh2FragLnRec[ij]);
    }
    fHistList->Add(fhnCorrelation);
    fHistList->Add(fhnCorrelationPhiZRec);
    fHistList->Add(fh2JetPtJetPhi);
    fHistList->Add(fh2TrackPtTrackPhi);

    fHistList->Add(fh2DijetDeltaPhiPt);       
    fHistList->Add(fh2DijetAsymPt);       
    fHistList->Add(fh2DijetAsymPtCut);               
    fHistList->Add(fh2DijetDeltaPhiDeltaEta);        
    fHistList->Add(fh2DijetPt2vsPt1);                
    fHistList->Add(fh2DijetDifvsSum);                
    fHistList->Add(fh1DijetMinv);                    
    fHistList->Add(fh1DijetMinvCut);                 
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

void AliAnalysisTaskJetSpectrum2::UserExec(Option_t */*option*/)
{

  Bool_t selected = kTRUE;

  if(fUseGlobalSelection&&fEventSelectionMask==0){
    selected = AliAnalysisHelperJetTasks::Selected();
  }
  else if(fUseGlobalSelection&&fEventSelectionMask>0){
    selected = AliAnalysisHelperJetTasks::TestSelectInfo(fEventSelectionMask);
  }

  if(!selected){
    // no selection by the service task, we continue
    if (fDebug > 1)Printf("Not selected %s:%d",(char*)__FILE__,__LINE__);
    PostData(1, fHistList);
    return;
  }


  //
  // Execute analysis for current event
  //
  AliESDEvent *fESD = 0;
  if(fUseAODJetInput){    
    fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
    if(!fAOD){
      Printf("%s:%d AODEvent not found in Input Manager %d",(char*)__FILE__,__LINE__,fUseAODJetInput);
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
  



  if (fDebug > 1)printf("Analysing event # %5d\n", (Int_t) fEntry);

  
  AliAODHandler *aodH = dynamic_cast<AliAODHandler*>(AliAnalysisManager::GetAnalysisManager()->GetOutputEventHandler());

  if(!aodH){
    Printf("%s:%d no output aodHandler found Jet",(char*)__FILE__,__LINE__);
    return;
  }

  if (fDebug > 10)Printf("%s:%d",(char*)__FILE__,__LINE__);
  TClonesArray *aodRecJets = dynamic_cast<TClonesArray*>(fAOD->FindListObject(fBranchRec.Data()));
  if(!aodRecJets){
    Printf("%s:%d no reconstructed Jet array with name %s in AOD",(char*)__FILE__,__LINE__,fBranchRec.Data());
    return;
  }

  // ==== General variables needed


  // We use statice array, not to fragment the memory
  AliAODJet genJets[kMaxJets];
  Int_t nGenJets = 0;
  AliAODJet recJets[kMaxJets];
  Int_t nRecJets = 0;
  ///////////////////////////


  Double_t eventW = 1;
  Double_t ptHard = 0; 
  Double_t nTrials = 1; // Trials for MC trigger 

  if(fUseExternalWeightOnly){
    eventW = fExternalWeight;
  }

  fh1Trials->Fill("#sum{ntrials}",fAvgTrials); 
  //  if(fDebug>0)aodH->SetFillAOD(kFALSE);
  if (fDebug > 10)Printf("%s:%d",(char*)__FILE__,__LINE__);
  if((fAnalysisType&kAnaMCESD)==kAnaMCESD){
    // this is the part we only use when we have MC information
    AliMCEvent* mcEvent = MCEvent();
    //    AliStack *pStack = 0; 
    if(!mcEvent){
      Printf("%s:%d no mcEvent",(char*)__FILE__,__LINE__);
      return;
    }
    AliGenPythiaEventHeader*  pythiaGenHeader = AliAnalysisHelperJetTasks::GetPythiaEventHeader(mcEvent);
    Int_t iCount = 0;  
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
      
      // fetch the pythia generated jets only to be used here
      Int_t nPythiaGenJets = pythiaGenHeader->NTriggerJets();
      AliAODJet pythiaGenJets[kMaxJets];
      for(int ip = 0;ip < nPythiaGenJets;++ip){
	if(iCount>=kMaxJets)continue;
	Float_t p[4];
	pythiaGenHeader->TriggerJet(ip,p);
	pythiaGenJets[iCount].SetPxPyPzE(p[0],p[1],p[2],p[3]);
	
	if(fBranchGen.Length()==0){
	  /*	
	  if(fLimitGenJetEta){
	    if(pythiaGenJets[iCount].Eta()>fJetHeaderRec->GetJetEtaMax()||
	       pythiaGenJets[iCount].Eta()<fJetHeaderRec->GetJetEtaMin())continue;
	  }
	  */
	  if(fMinJetPt>0&&pythiaGenJets[iCount].Pt()<fMinJetPt)continue;
	  // if we have MC particles and we do not read from the aod branch
	  // use the pythia jets
	  genJets[iCount].SetPxPyPzE(p[0],p[1],p[2],p[3]);
	  iCount++;
	}
      }
    }
    if(fBranchGen.Length()==0)nGenJets = iCount;    
  }// (fAnalysisType&kMCESD)==kMCESD)


  // we simply fetch the tracks/mc particles as a list of AliVParticles

  TList recParticles;
  TList genParticles;




  Int_t nT = GetListOfTracks(&recParticles,fTrackTypeRec);
  if(fDebug>2)Printf("%s:%d Selected Rec tracks: %d %d",(char*)__FILE__,__LINE__,nT,recParticles.GetEntries());
  nT = GetListOfTracks(&genParticles,fTrackTypeGen);
  if(fDebug>2)Printf("%s:%d Selected Gen tracks: %d %d",(char*)__FILE__,__LINE__,nT,genParticles.GetEntries());


  if (fDebug > 10)Printf("%s:%d",(char*)__FILE__,__LINE__);
  fh1PtHard->Fill(ptHard,eventW);
  fh1PtHardNoW->Fill(ptHard,1);
  fh1PtHardTrials->Fill(ptHard,nTrials);

  // If we set a second branch for the input jets fetch this 
  if(fBranchGen.Length()>0){
    TClonesArray *aodGenJets = dynamic_cast<TClonesArray*>(fAOD->FindListObject(fBranchGen.Data()));
    if(aodGenJets){
      Int_t iCount = 0;
      for(int ig = 0;ig < aodGenJets->GetEntries();++ig){
	if(iCount>=kMaxJets)continue;
	AliAODJet *tmp = dynamic_cast<AliAODJet*>(aodGenJets->At(ig));
	if(!tmp)continue;
	/*
	if(fLimitGenJetEta){
	  if(tmp->Eta()>fJetHeaderRec->GetJetEtaMax()||
	     tmp->Eta()<fJetHeaderRec->GetJetEtaMin())continue;
	}
	*/
	if(fMinJetPt>0&&tmp->Pt()<fMinJetPt)continue;
	genJets[iCount] = *tmp;
	iCount++;
      }
      nGenJets = iCount;
    }
    else{
      if(fDebug>1)Printf("%s:%d Generated jet branch %s not found",(char*)__FILE__,__LINE__,fBranchGen.Data());
      if(fDebug>2)fAOD->Print();
    }
  }

  fh1NGenJets->Fill(nGenJets);
  // We do not want to exceed the maximum jet number
  nGenJets = TMath::Min(nGenJets,kMaxJets);

  // Fetch the reconstructed jets...

  nRecJets = aodRecJets->GetEntries();

  nRecJets = aodRecJets->GetEntries();
  fh1NRecJets->Fill(nRecJets);

  // Do something with the all rec jets
  Int_t nRecOver = nRecJets;

  // check that the jets are sorted
  Float_t ptOld = 999999;
  for(int ir = 0;ir < nRecJets;ir++){
    AliAODJet *tmp = (AliAODJet*)(aodRecJets->At(ir));
    Float_t tmpPt = tmp->Pt();
    if(tmpPt>ptOld){
      Printf("%s:%d Jets Not Sorted %s !! %d:%.3E %d:%.3E",(char*)__FILE__,__LINE__,fBranchRec.Data(),ir,tmpPt,ir-1,ptOld);
    }
    ptOld = tmpPt;
  }


  if(nRecOver>0){
    TIterator *recIter = aodRecJets->MakeIterator();
    AliAODJet *tmpRec = (AliAODJet*)(recIter->Next());  
    Float_t pt = tmpRec->Pt();
    if(tmpRec){
      for(int i = 1;i <= fh2NRecJetsPt->GetNbinsX();i++){
	Float_t ptCut = fh2NRecJetsPt->GetXaxis()->GetBinCenter(i);
	while(pt<ptCut&&tmpRec){
	  nRecOver--;
	  tmpRec = (AliAODJet*)(recIter->Next()); 
	  if(tmpRec){
	    pt = tmpRec->Pt();
	  }
	}
	if(nRecOver<=0)break;
	fh2NRecJetsPt->Fill(ptCut,nRecOver);
      }
    }
    recIter->Reset();
    AliAODJet *leading = (AliAODJet*)aodRecJets->At(0);
    Float_t phi = leading->Phi();
    if(phi<0)phi+=TMath::Pi()*2.;    
    Float_t eta = leading->Eta();
    pt = leading->Pt();
    while((tmpRec = (AliAODJet*)(recIter->Next()))){
      Float_t tmpPt = tmpRec->Pt();
      fh1PtJetsRecIn->Fill(tmpPt);
      if(tmpRec==leading){
	fh1PtJetsLeadingRecIn->Fill(tmpPt);
	continue;
      }
      // correlation
      Float_t tmpPhi =  tmpRec->Phi();

      if(tmpPhi<0)tmpPhi+=TMath::Pi()*2.;    
      Float_t dPhi = phi - tmpRec->Phi();
      if(dPhi>TMath::Pi())dPhi = dPhi - 2.*TMath::Pi();
      if(dPhi<(-1.*TMath::Pi()))dPhi = dPhi + 2.*TMath::Pi();      
      Float_t dEta = eta - tmpRec->Eta();
      fh2JetsLeadingPhiEta->Fill(dPhi,dEta);
      fh2JetsLeadingPhiPt->Fill(dPhi,pt);
    }  
    delete recIter;
  }

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
      fh1PtTracksRecIn->Fill(tmpPt);
      if(tmpRec==leading){
	fh1PtTracksLeadingRecIn->Fill(tmpPt);
	continue;
      }
      // correlation
      Float_t tmpPhi =  tmpRec->Phi();

      if(tmpPhi<0)tmpPhi+=TMath::Pi()*2.;    
      Float_t dPhi = phi - tmpRec->Phi();
      if(dPhi>TMath::Pi())dPhi = dPhi - 2.*TMath::Pi();
      if(dPhi<(-1.*TMath::Pi()))dPhi = dPhi + 2.*TMath::Pi();      
      Float_t dEta = eta - tmpRec->Eta();
      fh2TracksLeadingPhiEta->Fill(dPhi,dEta);
      fh2TracksLeadingPhiPt->Fill(dPhi,pt);
    }  
    delete recIter;
  }
  
  if(genParticles.GetSize()){
    TIterator *genIter = genParticles.MakeIterator();
    AliVParticle *tmpGen = 0;
    while((tmpGen = (AliVParticle*)(genIter->Next()))){
      if(TMath::Abs(tmpGen->Eta())<0.9){
	Float_t tmpPt = tmpGen->Pt();
	fh1PtTracksGenIn->Fill(tmpPt);
      }  
    }
    delete genIter;
  }
  
  nRecJets = TMath::Min(nRecJets,kMaxJets);
  
  Int_t iCountRec = 0;
  for(int ir = 0;ir < nRecJets;++ir){
    AliAODJet *tmp = dynamic_cast<AliAODJet*>(aodRecJets->At(ir));
    if(!tmp)continue;
    if(tmp->Pt()<fMinJetPt)continue;
    recJets[ir] = *tmp;
    iCountRec++;
  }
  nRecJets = iCountRec;


  if (fDebug > 10)Printf("%s:%d",(char*)__FILE__,__LINE__);
  // Relate the jets
  Int_t iGenIndex[kMaxJets];    // Index of the generated jet for i-th rec -1 if none
  Int_t iRecIndex[kMaxJets];    // Index of the rec jet for i-th gen -1 if none
  

  for(int i = 0;i<kMaxJets;++i){    
    iGenIndex[i] = iRecIndex[i] = -1;
  }

  AliAnalysisHelperJetTasks::GetClosestJets(genJets,nGenJets,recJets,nRecJets,
					    iGenIndex,iRecIndex,fDebug);
  if (fDebug > 10)Printf("%s:%d",(char*)__FILE__,__LINE__);

  if(fDebug){
    for(int i = 0;i<kMaxJets;++i){
      if(iGenIndex[i]>=0)Printf("iGenFound: %d -> %d",i,iGenIndex[i]); 
      if(iRecIndex[i]>=0)Printf("iRecFound: %d -> %d",i,iRecIndex[i]); 
    }
  }




  Double_t container[6];
  Double_t containerPhiZ[6];

  // loop over generated jets

  // radius; tmp, get from the jet header itself
  // at some point, how todeal woht FastJet on MC?
  Float_t radiusGen =0.4;
  Float_t radiusRec =0.4;

  for(int ig = 0;ig < nGenJets;++ig){
    Double_t ptGen = genJets[ig].Pt();
    Double_t phiGen = genJets[ig].Phi();
    if(phiGen<0)phiGen+=TMath::Pi()*2.;    
    Double_t etaGen = genJets[ig].Eta();
    
    container[3] = ptGen;
    container[4] = etaGen;
    container[5] = phiGen;
    fhnJetContainer[kStep0]->Fill(&container[3],eventW);
    Int_t ir = iRecIndex[ig];

    if(TMath::Abs(etaGen)<fRecEtaWindow){
      fh1TmpRho->Reset();

      fhnJetContainer[kStep1]->Fill(&container[3],eventW);
      fh1PtGenIn[ig]->Fill(ptGen,eventW);
      // fill the fragmentation function
      for(int it = 0;it<genParticles.GetEntries();++it){
	AliVParticle *part = (AliVParticle*)genParticles.At(it);
	Float_t deltaR = genJets[ig].DeltaR(part);
	fh1TmpRho->Fill(deltaR,part->Pt()/ptGen);
	if(deltaR<radiusGen){
	  Float_t z = part->Pt()/ptGen;
	  Float_t lnz =  -1.*TMath::Log(z);
	  fh2FragGen[ig]->Fill(z,ptGen,eventW);
	  fh2FragLnGen[ig]->Fill(lnz,ptGen,eventW);
	}

      }
      Float_t rhoSum = 0;
      for(int ibx = 1;ibx<fh2RhoPtGen[ir]->GetNbinsX();ibx++){
	Float_t r = fh2RhoPtGen[ir]->GetXaxis()->GetBinCenter(ibx);
	Float_t rho = fh1TmpRho->GetBinContent(ibx);
	rhoSum += rho;
	fh2RhoPtGen[ig]->Fill(r,ptGen,rho);
	fh2PsiPtGen[ig]->Fill(r,ptGen,rhoSum);
      }
    }
    if(ir>=0&&ir<nRecJets){   
      fhnJetContainer[kStep2]->Fill(&container[3],eventW);
      Double_t etaRec = recJets[ir].Eta();
      if(TMath::Abs(etaRec)<fRecEtaWindow)fhnJetContainer[kStep4]->Fill(&container[3],eventW);
    }
  }// loop over generated jets
  
  Float_t sumPt = 0;
  for(int it = 0;it<recParticles.GetEntries();++it){
    AliVParticle *part = (AliVParticle*)recParticles.At(it);
    // fill sum pt and P_t of all paritles
    if(TMath::Abs(part->Eta())<0.9){
      Float_t pt = part->Pt();
      fh1PtTrackRec->Fill(pt,eventW);
      fh2TrackPtTrackPhi->Fill(pt,part->Phi());
      sumPt += pt;
    }
  }
  fh1SumPtTrackAreaRec->Fill(sumPt*0.4*0.4/(2.*1.8),eventW);
  fh1SumPtTrackRec->Fill(sumPt,eventW);

  
  // loop over reconstructed jets
  for(int ir = 0;ir < nRecJets;++ir){
    Double_t etaRec = recJets[ir].Eta();
    Double_t ptRec = recJets[ir].Pt();
    Double_t phiRec = recJets[ir].Phi();
    if(phiRec<0)phiRec+=TMath::Pi()*2.;    


  // do something with dijets...
  if(ir==1){
    Double_t etaRec1 = recJets[0].Eta();
    Double_t ptRec1 = recJets[0].Pt();
    Double_t phiRec1 = recJets[0].Phi();
    if(phiRec1<0)phiRec1+=TMath::Pi()*2.;    
    
  
    if(TMath::Abs(etaRec1)<fRecEtaWindow
      &&TMath::Abs(etaRec)<fRecEtaWindow){
  
      Float_t deltaPhi = phiRec1 - phiRec;
      
      if(deltaPhi>TMath::Pi())deltaPhi = deltaPhi - 2.*TMath::Pi();
      if(deltaPhi<(-1.*TMath::Pi()))deltaPhi = deltaPhi + 2.*TMath::Pi();      
      deltaPhi = TMath::Abs(deltaPhi);
      fh2DijetDeltaPhiPt->Fill(deltaPhi,ptRec1);      
      Float_t asym = (ptRec1-ptRec)/(ptRec1+ptRec);
      fh2DijetAsymPt->Fill(asym,ptRec1);
      fh2DijetDeltaPhiDeltaEta->Fill(deltaPhi,etaRec1-etaRec);
      fh2DijetPt2vsPt1->Fill(ptRec1,ptRec);        
      fh2DijetDifvsSum->Fill(ptRec1+ptRec,ptRec1-ptRec);        
      Float_t minv = 2.*(recJets[0].P()*recJets[1].P()-
			 recJets[0].Px()*recJets[1].Px()- 
			 recJets[0].Py()*recJets[1].Py()- 
			 recJets[0].Pz()*recJets[1].Py());
      minv = TMath::Sqrt(minv);
      // with mass == 0;
      
      fh1DijetMinv->Fill(minv);            
      if((TMath::Pi()-deltaPhi)<fDeltaPhiWindow){
	fh1DijetMinvCut->Fill(minv);         
	fh2DijetAsymPtCut->Fill(asym,ptRec1);      
      }
    }
  }


    container[0] = ptRec;
    container[1] = etaRec;
    container[2] = phiRec;
    containerPhiZ[0] = ptRec;
    containerPhiZ[1] = phiRec;
    if(ptRec>30.&&fDebug>0){
      // need to cast to int, otherwise the printf overwrites
      Printf("Jet found in Event %d with p_T, %E",(int)Entry(),ptRec);
      Printf("%s read event, %d",fInputHandler->GetTree()->GetCurrentFile()->GetName(),(int)fInputHandler->GetTree()->GetTree()->GetReadEntry());
      if(fESD)Printf("ESDEvent  GetEventNumberInFile(): %d",fESD->GetEventNumberInFile());
      //  aodH->SetFillAOD(kTRUE);
      fAOD->GetHeader()->Print();
      Printf("TriggerClasses: %s",fAOD->GetFiredTriggerClasses().Data());
      for(int it = 0;it < fAOD->GetNumberOfTracks();++it){
	AliAODTrack *tr = fAOD->GetTrack(it);
	if((fFilterMask>0)&&!(tr->TestFilterBit(fFilterMask)))continue;
	tr->Print();
	//	tr->Dump();
	if(fESD){
	  AliESDtrack *esdTr = (AliESDtrack*)fESD->GetTrack(tr->GetID());
	  esdTr->Print("");
	  // esdTr->Dump();
	}
      }
    }
  

    fhnJetContainer[kStep0+kMaxStep]->Fill(container,eventW);
    if (fDebug > 10)Printf("%s:%d",(char*)__FILE__,__LINE__);
  
    Float_t zLeading = -1;
    if(TMath::Abs(etaRec)<fRecEtaWindow){
      fh2JetPtJetPhi->Fill(ptRec,phiRec);
      fhnJetContainer[kStep1+kMaxStep]->Fill(container,eventW);
      fh1PtRecIn[ir]->Fill(ptRec,eventW);
      // fill the fragmentation function
      
      fh1TmpRho->Reset();

      for(int it = 0;it<recParticles.GetEntries();++it){
	AliVParticle *part = (AliVParticle*)recParticles.At(it);
	Float_t eta = part->Eta();
	if(TMath::Abs(eta)<0.9){
	  Float_t phi = part->Phi();
	  if(phi<0)phi+=TMath::Pi()*2.;    
	  Float_t dPhi = phi - phiRec;
	  Float_t dEta = eta - etaRec;
	  if(dPhi>TMath::Pi())dPhi = dPhi - 2.*TMath::Pi();
	  if(dPhi<(-1.*TMath::Pi()))dPhi = dPhi + 2.*TMath::Pi();      
	  fh2PhiPt[ir]->Fill(dPhi,ptRec,eventW);
	  fh2PhiEta[ir]->Fill(dPhi,dEta,eventW);
	}

	Float_t deltaR = recJets[ir].DeltaR(part);
	fh1TmpRho->Fill(deltaR,part->Pt()/ptRec);


	if(deltaR<radiusRec){
	  Float_t z = part->Pt()/ptRec;
	  if(z>zLeading)zLeading=z;
	  Float_t lnz =  -1.*TMath::Log(z);
	  fh2FragRec[ir]->Fill(z,ptRec,eventW);
	  fh2FragLnRec[ir]->Fill(lnz,ptRec,eventW);
	}
      }
      // fill the jet shapes
      Float_t rhoSum = 0;
      for(int ibx = 1;ibx<fh2RhoPtRec[ir]->GetNbinsX();ibx++){
	Float_t r = fh2RhoPtRec[ir]->GetXaxis()->GetBinCenter(ibx);
	Float_t rho = fh1TmpRho->GetBinContent(ibx);
	rhoSum += rho;
	fh2RhoPtRec[ir]->Fill(r,ptRec,rho);
	fh2PsiPtRec[ir]->Fill(r,ptRec,rhoSum);
      }
    }


    containerPhiZ[2] = zLeading;

    // Fill Correlation
    Int_t ig = iGenIndex[ir];
    if(ig>=0 && ig<nGenJets){
      fhnJetContainer[kStep2+kMaxStep]->Fill(container,eventW);
      if (fDebug > 10)Printf("%s:%d ig = %d ir = %d",(char*)__FILE__,__LINE__,ig,ir);
      Double_t ptGen  = genJets[ig].Pt();
      Double_t phiGen = genJets[ig].Phi();
      if(phiGen<0)phiGen+=TMath::Pi()*2.; 
      Double_t etaGen = genJets[ig].Eta();

      container[3] = ptGen;
      container[4] = etaGen;
      container[5] = phiGen;
      containerPhiZ[3] = ptGen;
      // 
      // we accept only jets which are detected within a smaller window, to avoid ambigious pair association at the edges of the acceptance
      // 

      if(TMath::Abs(etaGen)<fRecEtaWindow)fhnJetContainer[kStep4+kMaxStep]->Fill(container,eventW);
      if(TMath::Abs(etaRec)<fRecEtaWindow){
	fhnJetContainer[kStep3+kMaxStep]->Fill(container,eventW);
	fhnCorrelation->Fill(container);
	if(fhnCorrelationPhiZRec)fhnCorrelationPhiZRec->Fill(containerPhiZ);

      }// if etarec in window

    } 
    else{
      containerPhiZ[3] = 0;
      if(fhnCorrelationPhiZRec)fhnCorrelationPhiZRec->Fill(containerPhiZ);
    }
  }// loop over reconstructed jets


  if (fDebug > 10)Printf("%s:%d",(char*)__FILE__,__LINE__);
  if (fDebug > 10)Printf("%s:%d",(char*)__FILE__,__LINE__);
  PostData(1, fHistList);
}


void AliAnalysisTaskJetSpectrum2::MakeJetContainer(){
  //
  // Create the particle container for the correction framework manager and 
  // link it
  //
  const Int_t kNvar   = 3 ; //number of variables on the grid:pt,eta, phi
  const Double_t kPtmin = 0.0, kPtmax = 160.; // we do not want to have empty bins at the beginning...
  const Double_t kEtamin = -3.0, kEtamax = 3.0;
  const Double_t kPhimin = 0., kPhimax = 2. * TMath::Pi();
  const Double_t kZmin = 0., kZmax = 1;

  // can we neglect migration in eta and phi?
  // phi should be no problem since we cover full phi and are phi symmetric
  // eta migration is more difficult  due to needed acceptance correction
  // in limited eta range

  //arrays for the number of bins in each dimension
  Int_t iBin[kNvar];
  iBin[0] = 160; //bins in pt
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


  for(int i = 0;i < kMaxStep*2;++i){
    fhnJetContainer[i] = new THnSparseF(Form("fhnJetContainer%d",i),Form("THnSparse jet info %d"),kNvar,iBin);
    for (int k=0; k<kNvar; k++) {
      fhnJetContainer[i]->SetBinEdges(k,binEdges[k]);
    }
  }
  //create correlation matrix for unfolding
  Int_t thnDim[2*kNvar];
  for (int k=0; k<kNvar; k++) {
    //first half  : reconstructed 
    //second half : MC
    thnDim[k]      = iBin[k];
    thnDim[k+kNvar] = iBin[k];
  }

  fhnCorrelation = new THnSparseF("fhnCorrelation","THnSparse with correlations",2*kNvar,thnDim);
  for (int k=0; k<kNvar; k++) {
    fhnCorrelation->SetBinEdges(k,binEdges[k]);
    fhnCorrelation->SetBinEdges(k+kNvar,binEdges[k]);
  }
  fhnCorrelation->Sumw2();

  // for second correlation histogram


  const Int_t kNvarPhiZ   = 4; 
  //arrays for the number of bins in each dimension
  Int_t iBinPhiZ[kNvarPhiZ];
  iBinPhiZ[0] = 80; //bins in pt
  iBinPhiZ[1] = 72; //bins in phi 
  iBinPhiZ[2] = 20; // bins in Z
  iBinPhiZ[3] = 80; //bins in ptgen

  //arrays for lower bounds :
  Double_t* binEdgesPhiZ[kNvarPhiZ];
  for(Int_t ivar = 0; ivar < kNvarPhiZ; ivar++)
    binEdgesPhiZ[ivar] = new Double_t[iBinPhiZ[ivar] + 1];

  for(Int_t i=0; i<=iBinPhiZ[0]; i++) binEdgesPhiZ[0][i]=(Double_t)kPtmin  + (kPtmax-kPtmin)/(Double_t)iBinPhiZ[0]*(Double_t)i;
  for(Int_t i=0; i<=iBinPhiZ[1]; i++) binEdgesPhiZ[1][i]=(Double_t)kPhimin  + (kPhimax-kPhimin)/iBinPhiZ[1]*(Double_t)i;
  for(Int_t i=0; i<=iBinPhiZ[2]; i++) binEdgesPhiZ[2][i]=(Double_t)kZmin  + (kZmax-kZmin)/iBinPhiZ[2]*(Double_t)i;
  for(Int_t i=0; i<=iBinPhiZ[3]; i++) binEdgesPhiZ[3][i]=(Double_t)kPtmin  + (kPtmax-kPtmin)/(Double_t)iBinPhiZ[3]*(Double_t)i;

  fhnCorrelationPhiZRec = new THnSparseF("fhnCorrelationPhiZRec","THnSparse with correlations",kNvarPhiZ,iBinPhiZ);
  for (int k=0; k<kNvarPhiZ; k++) {
    fhnCorrelationPhiZRec->SetBinEdges(k,binEdgesPhiZ[k]);
  }
  fhnCorrelationPhiZRec->Sumw2();


  // Add a histogram for Fake jets

  for(Int_t ivar = 0; ivar < kNvar; ivar++)
    delete [] binEdges[ivar];

  for(Int_t ivar = 0; ivar < kNvarPhiZ; ivar++)
    delete [] binEdgesPhiZ[ivar];

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
      if(TMath::Abs(tr->Eta())>0.9)continue;
      if(fDebug>0){
	if(tr->Pt()>20){
	  Printf("High pT track found in Event %d with p_T, %E",(int)Entry(),tr->Pt());
	  Printf("%s read event, %d",fInputHandler->GetTree()->GetCurrentFile()->GetName(),fInputHandler->GetTree()->GetReadEntry());
	  tr->Print();
	  //	tr->Dump();
	  AliESDEvent *fESD = dynamic_cast<AliESDEvent*> (InputEvent());
	  if(fESD){
	    AliESDtrack *esdTr = (AliESDtrack*)fESD->GetTrack(tr->GetID());
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
