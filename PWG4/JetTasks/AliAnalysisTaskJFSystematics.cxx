// **************************************
// Task used for the systematic study of jet finders
//
// Compares input (gen) and output (rec) jets   
// gen can also be another jet finder on the rec level, matching is done in eta phi
//
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
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TProfile.h>
#include <TList.h>
#include <TKey.h>
#include <TLorentzVector.h>
#include <TClonesArray.h>
#include  "TDatabasePDG.h"

#include "AliAnalysisTaskJFSystematics.h"
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
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliStack.h"
#include "AliGenPythiaEventHeader.h"
#include "AliJetKineReaderHeader.h"
#include "AliGenCocktailEventHeader.h"
#include "AliInputEventHandler.h"


#include "AliAnalysisHelperJetTasks.h"

ClassImp(AliAnalysisTaskJFSystematics)

  const Int_t AliAnalysisTaskJFSystematics::fgkSysBins[AliAnalysisTaskJFSystematics::kSysTypes] = {0,AliAnalysisTaskJFSystematics::kMaxJets};
const char* AliAnalysisTaskJFSystematics::fgkSysName[AliAnalysisTaskJFSystematics::kSysTypes] = {"","j"};

AliAnalysisTaskJFSystematics::AliAnalysisTaskJFSystematics(): AliAnalysisTaskSE(),
  fJetHeaderRec(0x0),
  fJetHeaderGen(0x0),
  fAOD(0x0),
  fBranchRec("jets"),
  fBranchGen(""),
  fUseAODInput(kFALSE),
  fUseExternalWeightOnly(kFALSE),
  fLimitGenJetEta(kFALSE),
  fAnalysisType(0),
  fExternalWeight(1),    
  fRecEtaWindow(0.5),
  fAvgTrials(1),
  fh1Xsec(0x0),
  fh1Trials(0x0),
  fh1PtHard(0x0),
  fh1PtHardNoW(0x0),  
  fh1PtHardTrials(0x0),
  fh1NGenJets(0x0),
  fh1NRecJets(0x0),
  fh1PtRecIn(0x0),						     
  fh1PtRecOut(0x0),     
  fh1PtGenIn(0x0),      
  fh1PtGenOut(0x0),     
  fh2PtFGen(0x0),       
  fh2PhiFGen(0x0),      
  fh2EtaFGen(0x0),      
  fh2PtGenDeltaPhi(0x0),
  fh2PtGenDeltaEta(0x0),
  fh3RecInEtaPhiPt(0x0),   
  fh3RecOutEtaPhiPt(0x0),  
  fh3GenInEtaPhiPt(0x0),   
  fh3GenOutEtaPhiPt(0x0),  
  fhnCorrelation(0x0),     
  fHistList(0x0)
{
  // Default constructor
  /*
    for(int ij  = 0;ij<kMaxJets;++ij){
      fh1PtRecIn[ij] = fh1PtRecOut[ij] = fh1PtGenIn[ij] = fh1PtGenOut[ij] = 0;
      fh2PtFGen[ij] = fh2PhiFGen[ij] = fh2EtaFGen[ij] = fh2PtGenDeltaPhi[ij] =  fh2PtGenDeltaEta[ij] = 0;
      fh3RecInEtaPhiPt[ij] = fh3RecOutEtaPhiPt[ij] =fh3GenInEtaPhiPt[ij] =  fh3GenOutEtaPhiPt[ij] = 0;
      fhnCorrelation[ij] = 0;
    }
  */  
}

AliAnalysisTaskJFSystematics::AliAnalysisTaskJFSystematics(const char* name):
  AliAnalysisTaskSE(name),
  fJetHeaderRec(0x0),
  fJetHeaderGen(0x0),
  fAOD(0x0),
  fBranchRec("jets"),
  fBranchGen(""),
  fUseAODInput(kFALSE),
  fUseExternalWeightOnly(kFALSE),
  fLimitGenJetEta(kFALSE),
  fAnalysisType(0),
  fExternalWeight(1),    
  fRecEtaWindow(0.5),
  fAvgTrials(1),
  fh1Xsec(0x0),
  fh1Trials(0x0),
  fh1PtHard(0x0),
  fh1PtHardNoW(0x0),  
  fh1PtHardTrials(0x0),
  fh1NGenJets(0x0),
  fh1NRecJets(0x0),
  fh1PtRecIn(0x0),						     
  fh1PtRecOut(0x0),     
  fh1PtGenIn(0x0),      
  fh1PtGenOut(0x0),     
  fh2PtFGen(0x0),       
  fh2PhiFGen(0x0),      
  fh2EtaFGen(0x0),      
  fh2PtGenDeltaPhi(0x0),
  fh2PtGenDeltaEta(0x0),
  fh3RecInEtaPhiPt(0x0),   
  fh3RecOutEtaPhiPt(0x0),  
  fh3GenInEtaPhiPt(0x0),   
  fh3GenOutEtaPhiPt(0x0),  
  fhnCorrelation(0x0),     
  fHistList(0x0) 
{
  // Default constructor

  // Default constructor
  /*
  for(int ij  = 0;ij<kMaxJets;++ij){
    fh1PtRecIn[ij] = fh1PtRecOut[ij] = fh1PtGenIn[ij] = fh1PtGenOut[ij] = 0;
    fh2PtFGen[ij] = fh2PhiFGen[ij] = fh2EtaFGen[ij] = fh2PtGenDeltaPhi[ij] =  fh2PtGenDeltaEta[ij] = 0;
    fh3RecInEtaPhiPt[ij] = fh3RecOutEtaPhiPt[ij] =fh3GenInEtaPhiPt[ij] =  fh3GenOutEtaPhiPt[ij] = 0;
    fhnCorrelation[ij] = 0;
  } 
  */ 
  DefineOutput(1, TList::Class());  
}



Bool_t AliAnalysisTaskJFSystematics::Notify()
{
//
  // Implemented Notify() to read the cross sections
  // and number of trials from pyxsec.root
  // 
  TTree *tree = AliAnalysisManager::GetAnalysisManager()->GetTree();
  Float_t xsection = 0;
  Float_t ftrials  = 1;
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
    if(ftrials>=nEntries)fAvgTrials = ftrials/nEntries; 
  }
  return kTRUE;
}

void AliAnalysisTaskJFSystematics::UserCreateOutputObjects()
{

  //
  // Create the output container
  //


  if (fDebug > 1) printf("AnalysisTaskJFSystematics::UserCreateOutputObjects() \n");

  OpenFile(1);
  if(!fHistList)fHistList = new TList();

  Bool_t oldStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);

  //
  //  Histograms
  // 

  const Int_t nBinPt = 100;
  Double_t binLimitsPt[nBinPt+1];
  for(Int_t iPt = 0;iPt <= nBinPt;iPt++){
    if(iPt == 0){
      binLimitsPt[iPt] = 0.0;
    }
    else {// 1.0
      binLimitsPt[iPt] =  binLimitsPt[iPt-1] + 2.5;
    }
  }
  
  const Int_t nBinEta = 26;
  Double_t binLimitsEta[nBinEta+1] = {
    -1.6,-1.4,-1.2,-1.0,
    -0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0.0,
    0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 
    1.0, 1.2, 1.4, 1.6
  };


  const Int_t nBinPhi = 18;
  Double_t binLimitsPhi[nBinPhi+1];
  for(Int_t iPhi = 0;iPhi<=nBinPhi;iPhi++){
    if(iPhi==0){
      binLimitsPhi[iPhi] = 0;
    }
    else{
      binLimitsPhi[iPhi] = binLimitsPhi[iPhi-1] + 1/(Float_t)nBinPhi * TMath::Pi()*2;
    }
  }


  fh1Xsec = new TProfile("fh1Xsec","xsec from pyxsec.root",1,0,1);
  fh1Xsec->GetXaxis()->SetBinLabel(1,"<#sigma>");

  fh1Trials = new TH1F("fh1Trials","trials event header or pyxsec file",1,0,1);
  fh1Trials->GetXaxis()->SetBinLabel(1,"#sum{ntrials}");

  fh1PtHard = new TH1F("fh1PtHard","PYTHIA Pt hard;p_{T,hard}",nBinPt,binLimitsPt);
  fh1PtHardNoW = new TH1F("fh1PtHardNoW","PYTHIA Pt hard no weight;p_{T,hard}",nBinPt,binLimitsPt);
  fh1PtHardTrials = new TH1F("fh1PtHardTrials","PYTHIA Pt hard weight with trials;p_{T,hard}",nBinPt,binLimitsPt);
  fh1NGenJets  = new TH1F("fh1NGenJets","N generated jets",20,-0.5,19.5);
  fh1NRecJets = new TH1F("fh1NRecJets","N reconstructed jets",20,-0.5,19.5);

  // book the single histograms these we clone for various systematics
  // 
   fh1PtRecIn =  new TH1F("fh1PtRecIn","rec p_T input ;p_{T,rec}",nBinPt,binLimitsPt);
  fh1PtRecOut = new TH1F("fh1PtRecOut","rec p_T output jets;p_{T,rec}",nBinPt,binLimitsPt);
  fh1PtGenIn = new TH1F("fh1PtGenIn","found p_T input ;p_{T,gen}",nBinPt,binLimitsPt);
  fh1PtGenOut = new TH1F("fh1PtGenOut","found p_T output jets;p_{T,gen}",nBinPt,binLimitsPt);



  fh2PtFGen = new TH2F("fh2PtFGen","Pt Found vs. gen;p_{T,rec} (GeV/c);p_{T,gen} (GeV/c)",
			     nBinPt,binLimitsPt,nBinPt,binLimitsPt);
  fh2PhiFGen = new TH2F("fh2PhiFGen","#phi Found vs. gen;#phi_{rec};#phi_{gen}",
			nBinPhi,binLimitsPhi,nBinPhi,binLimitsPhi);
  fh2EtaFGen = new TH2F("fh2EtaFGen","#eta Found vs. gen;#eta_{rec};#eta_{gen}",
			nBinEta,binLimitsEta,nBinEta,binLimitsEta);
    
  fh2PtGenDeltaPhi = new TH2F("fh2PtGenDeltaPhi","delta phi vs. P_{T,gen};p_{T,gen} (GeV/c);#phi_{gen}-#phi_{rec}",
			      nBinPt,binLimitsPt,100,-1.0,1.0);
  fh2PtGenDeltaEta = new TH2F("fh2PtGenDeltaEta","delta eta vs. p_{T,gen};p_{T,gen} (GeV/c);#eta_{gen}-#eta_{rec}",
			      nBinPt,binLimitsPt,100,-1.0,1.0);


  fh3RecInEtaPhiPt = new TH3F("fh3RecInEtaPhiPt","Rec eta, phi, pt; #eta; #phi; p_{T,rec} (GeV/c)",
			      nBinEta,binLimitsEta,nBinPhi,binLimitsPhi,nBinPt,binLimitsPt);
  fh3RecOutEtaPhiPt = new TH3F("fh3RecOutEtaPhiPt","generated found jet Rec eta, phi, pt; #eta; #phi; p_{T,rec} (GeV/c)",
				   nBinEta,binLimitsEta,nBinPhi,binLimitsPhi,nBinPt,binLimitsPt);
  fh3GenInEtaPhiPt = new TH3F("fh3GenInEtaPhiPt","generated jet eta, phi, pt; #eta; #phi; p_{T,gen} (GeV/c)",
				  nBinEta,binLimitsEta,nBinPhi,binLimitsPhi,nBinPt,binLimitsPt);
  fh3GenOutEtaPhiPt = new TH3F("fh3GenOutEtaPhiPt","reconstructed found for Gen eta, phi, pt; #eta; #phi; p_{T,} (GeV/c)",
				   nBinEta,binLimitsEta,nBinPhi,binLimitsPhi,nBinPt,binLimitsPt);
  
  const int nbin[4] = {nBinPt,nBinPt,24,24};
  Double_t vLowEdge[4] = {0,0,-1.2,-1.2};
  Double_t vUpEdge[4] = {250,250,1.2,1.2};
  
  fhnCorrelation  = new THnSparseF("fhnCorrelation","Response Map", 4, nbin, vLowEdge, vUpEdge);



  fHistList->Add(fh1Xsec);
  fHistList->Add(fh1Trials);
  fHistList->Add(fh1PtHard);
  fHistList->Add(fh1PtHardNoW);
  fHistList->Add(fh1PtHardTrials);
  fHistList->Add(fh1NGenJets);
  fHistList->Add(fh1NRecJets);
  fHistList->Add(fh1PtRecIn);
  fHistList->Add(fh1PtRecOut);
  fHistList->Add(fh1PtGenIn);
  fHistList->Add(fh1PtGenOut);
  fHistList->Add(fh2PtFGen);
  fHistList->Add(fh2PhiFGen);
  fHistList->Add(fh2EtaFGen);
  fHistList->Add(fh2PtGenDeltaEta);
  fHistList->Add(fh2PtGenDeltaPhi);
  fHistList->Add(fh3RecOutEtaPhiPt);
  fHistList->Add(fh3GenOutEtaPhiPt);      
  fHistList->Add(fh3RecInEtaPhiPt);
  fHistList->Add(fh3GenInEtaPhiPt);
  fHistList->Add(fhnCorrelation);


  if(fAnalysisType==kSysJetOrder){
    // 
    for(int i = 0; i< fgkSysBins[kSysJetOrder];++i){
      TH1F *hTmp = (TH1F*)fh1PtRecIn->Clone(Form("%s_%s%d",fh1PtRecIn->GetName(),fgkSysName[kSysJetOrder],i));
      fHistList->Add(hTmp);
      hTmp = (TH1F*)fh1PtRecOut->Clone(Form("%s_%s%d",fh1PtRecOut->GetName(),fgkSysName[kSysJetOrder],i));
      fHistList->Add(hTmp);
      hTmp = (TH1F*)fh1PtGenIn->Clone(Form("%s_%s%d",fh1PtGenIn->GetName(),fgkSysName[kSysJetOrder],i));
      fHistList->Add(hTmp);
      hTmp = (TH1F*)fh1PtGenOut->Clone(Form("%s_%s%d",fh1PtGenOut->GetName(),fgkSysName[kSysJetOrder],i));
      fHistList->Add(hTmp);
      THnSparseF *hnTmp = (THnSparseF*)fhnCorrelation->Clone(Form("%s_%s%d",fhnCorrelation->GetName(),fgkSysName[kSysJetOrder],i));
      fHistList->Add(hnTmp);
    }
  }

  // =========== Switch on Sumw2 for all histos ===========
  for (Int_t i=0; i<fHistList->GetEntries(); ++i) {
    TH1 *h1 = dynamic_cast<TH1*>(fHistList->At(i));
    if (h1){
      // Printf("%s ",h1->GetName());
      h1->Sumw2();
      continue;
    }
    THnSparse *hn = dynamic_cast<THnSparse*>(fHistList->At(i));
    if(hn)hn->Sumw2();
  }

  TH1::AddDirectory(oldStatus);

  }

void AliAnalysisTaskJFSystematics::Init()
{
  //
  // Initialization
  //

  Printf(">>> AnalysisTaskJFSystematics::Init() debug level %d\n",fDebug);
  if (fDebug > 1) printf("AnalysisTaskJFSystematics::Init() \n");

}

void AliAnalysisTaskJFSystematics::UserExec(Option_t */*option*/)
{
  //
  // Execute analysis for current event
  //
 
  if(fUseAODInput){    
    fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
    if(!fAOD){
      Printf("%s:%d AODEvent not found in Input Manager %d",(char*)__FILE__,__LINE__,fUseAODInput);
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
  }
  
  if (fDebug > 1)printf("AliAnalysisTaskJFSystematics::Analysing event # %5d\n", (Int_t) fEntry);

  // ========= These pointers need to be valid in any case ======= 

  TClonesArray *aodRecJets = dynamic_cast<TClonesArray*>(fAOD->FindListObject(fBranchRec.Data()));
  if(!aodRecJets){
    Printf("%s:%d no reconstructed Jet array with name %s in AOD",(char*)__FILE__,__LINE__,fBranchRec.Data());
    return;
  }

  // We use static arrays, not to fragment the memory
  // 
  AliAODJet genJets[kMaxJets];
  Int_t nGenJets = 0;
  AliAODJet recJets[kMaxJets];
  Int_t nRecJets = 0;

  Double_t eventW = 1;
  Double_t ptHard = 0; 
  
  Double_t nTrials = 1; // Trials for MC trigger weigth for real data
  Int_t     iProcessType = 0;
  if(fUseExternalWeightOnly){
    eventW = fExternalWeight;
  }

  if (fDebug > 10)Printf("%s:%d",(char*)__FILE__,__LINE__);
  // this is the part where when we need to have MC information
  // we can also work on Reconstructed only when just comparing two JF
  AliMCEvent* mcEvent = MCEvent();
  if(!mcEvent){
    Printf("%s:%d no mcEvent",(char*)__FILE__,__LINE__);
  }
  else {
    AliGenPythiaEventHeader*  pythiaGenHeader = AliAnalysisHelperJetTasks::GetPythiaEventHeader(mcEvent);
    if(pythiaGenHeader){
      nTrials = pythiaGenHeader->Trials();
      ptHard  = pythiaGenHeader->GetPtHard();
      iProcessType = pythiaGenHeader->ProcessType();
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
      Int_t iCount = 0;
      for(int ip = 0;ip < nPythiaGenJets;++ip){
	if(iCount>=kMaxJets)continue;
	Float_t p[4];
	pythiaGenHeader->TriggerJet(ip,p);
	pythiaGenJets[iCount].SetPxPyPzE(p[0],p[1],p[2],p[3]);
	if(fLimitGenJetEta){
	  if(pythiaGenJets[iCount].Eta()>fJetHeaderRec->GetJetEtaMax()||
	     pythiaGenJets[iCount].Eta()<fJetHeaderRec->GetJetEtaMin())continue;
      }
	if(fBranchGen.Length()==0){
	  // if we have MC particles and we do not read from the aod branch
	  // use the pythia jets
	  genJets[iCount].SetPxPyPzE(p[0],p[1],p[2],p[3]);
	}
      iCount++;
      }
      if(fBranchGen.Length()==0)nGenJets = iCount;    
    }
  }// if we had the MCEvent

  if(nTrials==1&&fAvgTrials>1) fh1Trials->Fill("#sum{ntrials}",fAvgTrials); 
  else fh1Trials->Fill("#sum{ntrials}",nTrials); 

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
	if(fLimitGenJetEta){
	  if(tmp->Eta()>fJetHeaderRec->GetJetEtaMax()||
	     tmp->Eta()<fJetHeaderRec->GetJetEtaMin())continue;
	}
	genJets[iCount] = *tmp;
	iCount++;
      }
      nGenJets = iCount;
    }
    else{
      Printf("%s:%d Generated jet branch %s not found",(char*)__FILE__,__LINE__,fBranchGen.Data());
    }
  }

  fh1NGenJets->Fill(nGenJets);
  // We do not want to exceed the maximum jet number
  nGenJets = TMath::Min(nGenJets,kMaxJets);


  //
  // Fetch the reconstructed jets...
  //

  nRecJets = aodRecJets->GetEntries();
  fh1NRecJets->Fill(nRecJets);
  nRecJets = TMath::Min(nRecJets,kMaxJets);

  for(int ir = 0;ir < nRecJets;++ir){
    AliAODJet *tmp = dynamic_cast<AliAODJet*>(aodRecJets->At(ir));
    if(!tmp)continue;
    recJets[ir] = *tmp;
  }


  if (fDebug > 10)Printf("%s:%d",(char*)__FILE__,__LINE__);

  //
  // Relate the jets
  //
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

  //
  // Now the premliminaries are over, lets do the jet analysis
  //


  Double_t value[4]; // for the thnsparse 
  // loop over reconstructed jets
  for(int ir = 0;ir < nRecJets;++ir){
    Double_t ptRec = recJets[ir].Pt();
    Double_t phiRec = recJets[ir].Phi();
    if(phiRec<0)phiRec+=TMath::Pi()*2.;    
    Double_t etaRec = recJets[ir].Eta();
    if (fDebug > 10)Printf("%s:%d",(char*)__FILE__,__LINE__);
    fh1PtRecIn->Fill(ptRec,eventW);
    if(fAnalysisType==kSysJetOrder)((TH1F*)fHistList->FindObject(Form("fh1PtRecIn_%s%d",fgkSysName[kSysJetOrder],ir)))->Fill(ptRec,eventW);
    fh3RecInEtaPhiPt->Fill(etaRec,phiRec,ptRec,eventW);
    // Fill Correlation
    Int_t ig = iGenIndex[ir];
    if(ig>=0 && ig<nGenJets){
      if (fDebug > 10)Printf("%s:%d ig = %d ir = %d",(char*)__FILE__,__LINE__,ig,ir);
      fh1PtRecOut->Fill(ptRec,eventW);
      if(fAnalysisType==kSysJetOrder)((TH1F*)fHistList->FindObject(Form("fh1PtRecOut_%s%d",fgkSysName[kSysJetOrder],ir)))->Fill(ptRec,eventW);
      Double_t ptGen  = genJets[ig].Pt();
      Double_t phiGen = genJets[ig].Phi();
      if(phiGen<0)phiGen+=TMath::Pi()*2.; 
      Double_t etaGen = genJets[ig].Eta();

      fh3RecOutEtaPhiPt->Fill(etaRec,phiRec,ptRec,eventW);

      value[0] = ptRec;
      value[1] = ptGen;
      value[2] = etaRec;
      value[3] = etaGen;
      
      fhnCorrelation->Fill(value,eventW);
      if(fAnalysisType==kSysJetOrder)((THnSparseF*)fHistList->FindObject(Form("fhnCorrelation_%s%d",fgkSysName[kSysJetOrder],ir)))->Fill(value,eventW);
      // 
      // we accept only jets which are detected within a smaller window, to 
      // avoid ambigious pair association at the edges of the acceptance
      // 

      if(TMath::Abs(etaRec)<fRecEtaWindow){
	fh2PtFGen->Fill(ptRec,ptGen,eventW);
	fh2PhiFGen->Fill(phiRec,phiGen,eventW);
	fh2EtaFGen->Fill(etaRec,etaGen,eventW);
	fh2PtGenDeltaEta->Fill(ptGen,etaGen-etaRec,eventW);
	fh2PtGenDeltaPhi->Fill(ptGen,phiGen-phiRec,eventW);
      }// if etarec in window
    }//if ig valid
  }// loop over reconstructed jets
  

  // Now llop over generated jets

  if (fDebug > 10)Printf("%s:%d",(char*)__FILE__,__LINE__);
  for(int ig = 0;ig < nGenJets;++ig){
    Double_t ptGen = genJets[ig].Pt();
    // Fill Correlation
    Double_t phiGen = genJets[ig].Phi();
    if(phiGen<0)phiGen+=TMath::Pi()*2.;    
    Double_t etaGen = genJets[ig].Eta();
    fh3GenInEtaPhiPt->Fill(etaGen,phiGen,ptGen,eventW);
    fh1PtGenIn->Fill(ptGen,eventW);
    if (fDebug > 10)Printf("%s:%d",(char*)__FILE__,__LINE__);
    if(fAnalysisType==kSysJetOrder)((TH1F*)fHistList->FindObject(Form("fh1PtGenIn_%s%d",fgkSysName[kSysJetOrder],ig)))->Fill(ptGen,eventW);
    Int_t ir = iRecIndex[ig];
    if(ir>=0&&ir<nRecJets){
      fh1PtGenOut->Fill(ptGen,eventW);
      fh3GenOutEtaPhiPt->Fill(etaGen,phiGen,ptGen,eventW);
      if (fDebug > 10)Printf("%s:%d",(char*)__FILE__,__LINE__);
      if(fAnalysisType==kSysJetOrder)((TH1F*)fHistList->FindObject(Form("fh1PtGenOut_%s%d",fgkSysName[kSysJetOrder],ig)))->Fill(ptGen,eventW);
    }
  }// loop over reconstructed jets

  if (fDebug > 10)Printf("%s:%d",(char*)__FILE__,__LINE__);
  PostData(1, fHistList);
}

void AliAnalysisTaskJFSystematics::Terminate(Option_t */*option*/)
{
// Terminate analysis
//
    if (fDebug > 1) printf("AnalysisTaskJFSystematics: Terminate() \n");
    // Plot 


}
