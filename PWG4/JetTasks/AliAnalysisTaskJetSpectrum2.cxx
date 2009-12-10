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
  fBranchRec("jets"),
  fBranchGen(""),
  fUseAODInput(kFALSE),
  fUseGlobalSelection(kFALSE),
  fUseExternalWeightOnly(kFALSE),
  fLimitGenJetEta(kFALSE),
  fFilterMask(0),
  fAnalysisType(0),
  fTrackTypeRec(kTrackUndef),
  fTrackTypeGen(kTrackUndef),
  fAvgTrials(1),
  fExternalWeight(1),    
  fRecEtaWindow(0.5),
  fh1Xsec(0x0),
  fh1Trials(0x0),
  fh1PtHard(0x0),
  fh1PtHardNoW(0x0),  
  fh1PtHardTrials(0x0),
  fh1NGenJets(0x0),
  fh1NRecJets(0x0),
  fHistList(0x0)  
{
  for(int i = 0;i < kMaxStep*2;++i){
    fhnJetContainer[i] = 0;
  }
  for(int i = 0;i < kMaxJets;++i){
    fh2FragRec[i] = fh2FragLnRec[i] = fh2FragGen[i] = fh2FragLnGen[i] = 0;
    fh1PtRecIn[i] = fh1PtGenIn[i] = 0;
  }  

}

AliAnalysisTaskJetSpectrum2::AliAnalysisTaskJetSpectrum2(const char* name):
  AliAnalysisTaskSE(name),
  fJetHeaderRec(0x0),
  fJetHeaderGen(0x0),
  fAOD(0x0),
  fhnCorrelation(0x0),
  fBranchRec("jets"),
  fBranchGen(""),
  fUseAODInput(kFALSE),
  fUseGlobalSelection(kFALSE),
  fUseExternalWeightOnly(kFALSE),
  fLimitGenJetEta(kFALSE),
  fFilterMask(0),
  fAnalysisType(0),
  fTrackTypeRec(kTrackUndef),
  fTrackTypeGen(kTrackUndef),
  fAvgTrials(1),
  fExternalWeight(1),    
  fRecEtaWindow(0.5),
  fh1Xsec(0x0),
  fh1Trials(0x0),
  fh1PtHard(0x0),
  fh1PtHardNoW(0x0),  
  fh1PtHardTrials(0x0),
  fh1NGenJets(0x0),
  fh1NRecJets(0x0),
  fHistList(0x0)
{

  for(int i = 0;i < kMaxStep*2;++i){
    fhnJetContainer[i] = 0;
  }  
  for(int i = 0;i < kMaxJets;++i){
    fh2FragRec[i] = fh2FragLnRec[i] = fh2FragGen[i] = fh2FragLnGen[i] = 0;
    fh1PtRecIn[i] = fh1PtGenIn[i] = 0;
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
    if(ftrials>=nEntries)fAvgTrials = ftrials/nEntries; 
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

  Bool_t oldStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);

  //
  //  Histogram
    
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
  
  const Int_t nBinPhi = 30;
  Double_t binLimitsPhi[nBinPhi+1];
  for(Int_t iPhi = 0;iPhi<=nBinPhi;iPhi++){
    if(iPhi==0){
      binLimitsPhi[iPhi] = 0;
    }
    else{
      binLimitsPhi[iPhi] = binLimitsPhi[iPhi-1] + 1/(Float_t)nBinPhi * TMath::Pi()*2;
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

  // 
  for(int ij = 0;ij < kMaxJets;++ij){
    fh1PtRecIn[ij] = new TH1F(Form("fh1PtRecIn_j%d",ij),"rec p_T input ;p_{T,rec}",nBinPt,binLimitsPt);
    fh1PtGenIn[ij] = new TH1F(Form("fh1PtGenIn_j%d",ij),"found p_T input ;p_{T,gen}",nBinPt,binLimitsPt);

    fh2FragRec[ij] = new TH2F(Form("fh2FragRec_j%d",ij),"Jet Fragmentation;x=p_{T,i}/p_{T,jet};p_{T,jet};1/N_{jet}dN_{ch}/dx",
			   nBinFrag,0.,1.,nBinPt,binLimitsPt);
    fh2FragLnRec[ij] = new TH2F(Form("fh2FragLnRec_j%d",ij),"Jet Fragmentation Ln;#xi=ln(p_{T,jet}/p_{T,i});p_{T,jet}(GeV);1/N_{jet}dN_{ch}/d#xi",
			     nBinFrag,0.,10.,nBinPt,binLimitsPt);

    fh2FragGen[ij] = new TH2F(Form("fh2FragGen_j%d",ij),"Jet Fragmentation;x=p_{T,i}/p_{T,jet};p_{T,jet};1/N_{jet}dN_{ch}/dx",
			   nBinFrag,0.,1.,nBinPt,binLimitsPt);
    fh2FragLnGen[ij] = new TH2F(Form("fh2FragLnGen_j%d",ij),"Jet Fragmentation Ln;#xi=ln(p_{T,jet}/p_{T,i});p_{T,jet}(GeV);1/N_{jet}dN_{ch}/d#xi",
			     nBinFrag,0.,10.,nBinPt,binLimitsPt);
  }


  const Int_t saveLevel = 3; // large save level more histos
  if(saveLevel>0){
    fHistList->Add(fh1Xsec);
    fHistList->Add(fh1Trials);
    fHistList->Add(fh1PtHard);
    fHistList->Add(fh1PtHardNoW);
    fHistList->Add(fh1PtHardTrials);
    fHistList->Add(fh1NGenJets);
    fHistList->Add(fh1NRecJets);
    for(int i = 0;i<kMaxStep*2;++i)fHistList->Add(fhnJetContainer[i]);
    for(int ij = 0;ij<kMaxJets;++ij){
      fHistList->Add( fh1PtRecIn[ij]);
      fHistList->Add( fh1PtGenIn[ij]);
      fHistList->Add( fh2FragRec[ij]);
      fHistList->Add( fh2FragLnRec[ij]);
      fHistList->Add( fh2FragGen[ij]);
      fHistList->Add( fh2FragLnGen[ij]);
    }
    fHistList->Add(fhnCorrelation);
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

  Printf(">>> AnalysisTaskJetSpectrum2::Init() debug level %d\n",fDebug);
  if (fDebug > 1) printf("AnalysisTaskJetSpectrum2::Init() \n");

}

void AliAnalysisTaskJetSpectrum2::UserExec(Option_t */*option*/)
{

  if(! AliAnalysisHelperJetTasks::Selected()&&fUseGlobalSelection){
    // no selection by the service task, we continue
    if (fDebug > 10)Printf("%s:%d",(char*)__FILE__,__LINE__);
    PostData(1, fHistList);
    return;
  }

  //
  // Execute analysis for current event
  //
  AliESDEvent *fESD = 0;
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
    }
    if(fBranchGen.Length()==0)nGenJets = iCount;    
  }// (fAnalysisType&kMCESD)==kMCESD)


  // we simply fetch the tracks/mc particles as a list of AliVParticles

  TList recParticles;
  TList genParticles;

  GetListOfTracks(&recParticles,fTrackTypeRec);
  GetListOfTracks(&genParticles,fTrackTypeGen);

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
      fAOD->Print();
    }
  }

  fh1NGenJets->Fill(nGenJets);
  // We do not want to exceed the maximum jet number
  nGenJets = TMath::Min(nGenJets,kMaxJets);

  // Fetch the reconstructed jets...

  nRecJets = aodRecJets->GetEntries();

  


  fh1NRecJets->Fill(nRecJets);
  nRecJets = TMath::Min(nRecJets,kMaxJets);

  for(int ir = 0;ir < nRecJets;++ir){
    AliAODJet *tmp = dynamic_cast<AliAODJet*>(aodRecJets->At(ir));
    if(!tmp)continue;
    recJets[ir] = *tmp;
  }


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
      fhnJetContainer[kStep1]->Fill(&container[3],eventW);
      fh1PtGenIn[ig]->Fill(ptGen,eventW);
      // fill the fragmentation function
      for(int it = 0;it<genParticles.GetEntries();++it){
	AliVParticle *part = (AliVParticle*)genParticles.At(it);
	if(genJets[ig].DeltaR(part)<radiusGen){
	  Float_t z = part->Pt()/ptGen;
	  Float_t lnz =  -1.*TMath::Log(z);
	  fh2FragGen[ig]->Fill(z,ptGen,eventW);
	  fh2FragLnGen[ig]->Fill(lnz,ptGen,eventW);
	}
      }
      if(ir>=0&&ir<nRecJets){
	fhnJetContainer[kStep3]->Fill(&container[3],eventW);
      }
    }    

    if(ir>=0&&ir<nRecJets){   
      fhnJetContainer[kStep2]->Fill(&container[3],eventW);
      Double_t etaRec = recJets[ir].Eta();
      if(TMath::Abs(etaRec)<fRecEtaWindow)fhnJetContainer[kStep4]->Fill(&container[3],eventW);
    }
  }// loop over generated jets

  
  

  // loop over reconstructed jets
  for(int ir = 0;ir < nRecJets;++ir){
    Double_t etaRec = recJets[ir].Eta();
    Double_t ptRec = recJets[ir].Pt();
    Double_t phiRec = recJets[ir].Phi();
    if(phiRec<0)phiRec+=TMath::Pi()*2.;    
    container[0] = ptRec;
    container[1] = etaRec;
    container[2] = phiRec;

    if(ptRec>20.&&fDebug>0){
      // need to cast to int, otherwise the printf overwrites
      Printf("Jet found in Event %d with p_T, %E",(int)Entry(),ptRec);
      fAOD->GetHeader()->Print();
      for(int it = 0;it < fAOD->GetNumberOfTracks();++it){
	AliAODTrack *tr = fAOD->GetTrack(it);
	if((fFilterMask>0)&&!(tr->TestFilterBit(fFilterMask)))continue;
	tr->Print();
	tr->Dump();
	if(fESD){
	  AliESDtrack *esdTr = (AliESDtrack*)fESD->GetTrack(tr->GetID());
	  esdTr->Print("");
	  esdTr->Dump();
	}
      }
    }
  

    fhnJetContainer[kStep0+kMaxStep]->Fill(container,eventW);
    if (fDebug > 10)Printf("%s:%d",(char*)__FILE__,__LINE__);
    if(TMath::Abs(etaRec)<fRecEtaWindow){
      fhnJetContainer[kStep1+kMaxStep]->Fill(container,eventW);
      fh1PtRecIn[ir]->Fill(ptRec,eventW);
      // fill the fragmentation function
      for(int it = 0;it<recParticles.GetEntries();++it){
	AliVParticle *part = (AliVParticle*)recParticles.At(it);
	if(recJets[ir].DeltaR(part)<radiusRec){
	  Float_t z = part->Pt()/ptRec;
	  Float_t lnz =  -1.*TMath::Log(z);
	  fh2FragRec[ir]->Fill(z,ptRec,eventW);
	  fh2FragLnRec[ir]->Fill(lnz,ptRec,eventW);
	}
      }
    }


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

      // 
      // we accept only jets which are detected within a smaller window, to avoid ambigious pair association at the edges of the acceptance
      // 

      if(TMath::Abs(etaGen)<fRecEtaWindow)fhnJetContainer[kStep4+kMaxStep]->Fill(container,eventW);
      if(TMath::Abs(etaRec)<fRecEtaWindow){
	fhnJetContainer[kStep3+kMaxStep]->Fill(container,eventW);
	fhnCorrelation->Fill(container);
      }// if etarec in window

    } 
    ////////////////////////////////////////////////////  
    else{

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
  const Double_t kPtmin = 10.0, kPtmax = 260.; // we do not want to have empty bins at the beginning...
  const Double_t kEtamin = -3.0, kEtamax = 3.0;
  const Double_t kPhimin = 0., kPhimax = 2. * TMath::Pi();

  // can we neglect migration in eta and phi?
  // phi should be no problem since we cover full phi and are phi symmetric
  // eta migration is more difficult  due to needed acceptance correction
  // in limited eta range


  //arrays for the number of bins in each dimension
  Int_t iBin[kNvar];
  iBin[0] = 100; //bins in pt
  iBin[1] =  1; //bins in eta 
  iBin[2] = 1; // bins in phi

  //arrays for lower bounds :
  Double_t* binEdges[kNvar];
  for(Int_t ivar = 0; ivar < kNvar; ivar++)
    binEdges[ivar] = new Double_t[iBin[ivar] + 1];

  //values for bin lower bounds
  //  for(Int_t i=0; i<=iBin[0]; i++) binEdges[0][i]=(Double_t)TMath::Power(10,TMath::Log10(kPtmin) + (TMath::Log10(kPtmax)-TMath::Log10(kPtmin))/iBin[0]*(Double_t)i);  
  for(Int_t i=0; i<=iBin[0]; i++) binEdges[0][i]=(Double_t)kPtmin  + (kPtmax-kPtmin)/iBin[0]*(Double_t)i;
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

  // Add a histogram for Fake jets
  //  thnDim[3] = AliPID::kSPECIES;
  //  fFakeElectrons = new THnSparseF("fakeEkectrons", "Output for Fake Electrons", kNvar + 1, thnDim);
  // for(Int_t idim = 0; idim < kNvar; idim++)
  //  fFakeElectrons->SetBinEdges(idim, binEdges[idim]);
}

void AliAnalysisTaskJetSpectrum2::Terminate(Option_t */*option*/)
{
// Terminate analysis
//
    if (fDebug > 1) printf("AnalysisJetSpectrum2: Terminate() \n");
}


Int_t  AliAnalysisTaskJetSpectrum2::GetListOfTracks(TList *list,Int_t type){


  Int_t iCount = 0;
  if(type==kTrackAODIn||type==kTrackAODOut){
    AliAODEvent *aod = 0;
    if(type==kTrackAODIn)aod = dynamic_cast<AliAODEvent*>(InputEvent());
    else if (type==kTrackAODOut) aod = AODEvent();
    if(!aod){
      return iCount;
    }
    for(int it = 0;it < aod->GetNumberOfTracks();++it){
      AliAODTrack *tr = aod->GetTrack(it);
      if((fFilterMask>0)&&!(tr->TestFilterBit(fFilterMask)))continue;
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
  else if (type == kTrackAODMCCharged || type == kTrackAODMCAll ) {
    AliAODEvent *aod = dynamic_cast<AliAODEvent*>(InputEvent());
    if(!aod) aod = AODEvent();
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
      else if (type == kTrackAODMCCharged){
	if(part->Charge()==0)continue;
	list->Add(part);
	iCount++;
      }
    }
  }// AODMCparticle
  return iCount;

}
