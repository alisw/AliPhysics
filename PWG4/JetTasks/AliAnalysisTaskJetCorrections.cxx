// **************************************
// Task used for estimating a charged to neutral correction
// sona.pochybova@cern.ch
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


#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TList.h>
#include <TTree.h>
#include <TLorentzVector.h>
#include <TClonesArray.h>
#include <TRefArray.h>
#include <TVector3.h>
#include <TProfile.h>

#include "AliAnalysisTaskJetCorrections.h"
#include "AliAnalysisManager.h"
#include "AliAODEvent.h"
#include "AliAODVertex.h"
#include "AliAODHandler.h"
#include "AliAODTrack.h"
#include "AliAODJet.h"
#include "AliMCEvent.h"

#include "AliAnalysisHelperJetTasks.h"

//
//
// corrections to jet energy by sona 
// 
//

ClassImp(AliAnalysisTaskJetCorrections)

  AliAnalysisTaskJetCorrections::AliAnalysisTaskJetCorrections() : AliAnalysisTaskSE(),								   
								   fAOD(0x0),
								  
								   fBranchRec(""),
								   fBranchGen(""),
								   
								   fUseAODInput(kFALSE),

								   fR(0x0),
								   fList(0x0),
								   
								   fGlobVar(1),
								   fXsection(1),
								   
								   fhEGen(0x0),
								   fhERec(0x0),
								   
								   fhEGenRest(0x0),
								   fhERecRest(0x0),

								   fhEsumGenRest(0x0),
								   fhEsumRecRest(0x0),
								   
								   fhE2vsE1Gen(0x0),
								   fhE2vsE1Rec(0x0),
								   fhE2E1vsEsumGen(0x0),
								   fhE2E1vsEsumRec(0x0),
								   fhE2E1vsE1Gen(0x0),
								   fhE2E1vsE1Rec(0x0),
								   fhE2E1vsdPhiGen(0x0),
								   fhE2E1vsdPhiRec(0x0),

								   fhTrackBalance2(0x0),
								   fhTrackBalance3(0x0),

								   fhEt1Et22(0x0),
								   fhEt1Et23(0x0)

{
  //
  // ctor
  //
  for (Int_t i = 0; i < 3; i++)
    {
      fhECorrJet10[i] = 0;    
      fhECorrJet05[i] = 0;    
      fhECorrJet01[i] = 0;    
      fhECorrJet001[i] = 0;
      fhdEvsErec10[i] = 0;
      fhdEvsErec05[i] = 0;
      fhdEvsErec01[i] = 0;
      fhdEvsErec001[i] = 0;
      fhdPhidEta10[i] = 0;
      fhdPhidEta05[i] = 0;
      fhdPhidEta01[i] = 0;
      fhdPhidEta001[i] = 0;
      fhdPhidEtaPt10[i] = 0;
      fhdPhidEtaPt05[i] = 0;
      fhdPhidEtaPt01[i] = 0;
      fhdPhidEtaPt001[i] = 0;
    }
}

AliAnalysisTaskJetCorrections::AliAnalysisTaskJetCorrections(const char * name):
  AliAnalysisTaskSE(name),
  
  fAOD(0x0),
  
  fBranchRec(""),
  fBranchGen(""),
  
  fUseAODInput(kFALSE),

  fR(0x0),
  fList(0x0),

  fGlobVar(1),
  fXsection(1),

  fhEGen(0x0),
  fhERec(0x0),

  fhEGenRest(0x0),
  fhERecRest(0x0),
  
  fhEsumGenRest(0x0),
  fhEsumRecRest(0x0),
  
  fhE2vsE1Gen(0x0),
  fhE2vsE1Rec(0x0),
  fhE2E1vsEsumGen(0x0),
  fhE2E1vsEsumRec(0x0),
  fhE2E1vsE1Gen(0x0),
  fhE2E1vsE1Rec(0x0),
  fhE2E1vsdPhiGen(0x0),
  fhE2E1vsdPhiRec(0x0),

  fhTrackBalance2(0x0),
  fhTrackBalance3(0x0),

  fhEt1Et22(0x0),
  fhEt1Et23(0x0)
{
  //
  // ctor
  //
  for (Int_t i = 0; i < 3; i++)
    {
      fhECorrJet10[i] = 0;    
      fhECorrJet05[i] = 0;    
      fhECorrJet01[i] = 0;    
      fhECorrJet001[i] = 0;
      fhdEvsErec10[i] = 0;
      fhdEvsErec05[i] = 0;
      fhdEvsErec01[i] = 0;
      fhdEvsErec001[i] = 0;
      fhdPhidEta10[i] = 0;
      fhdPhidEta05[i] = 0;
      fhdPhidEta01[i] = 0;
      fhdPhidEta001[i] = 0;
      fhdPhidEtaPt10[i] = 0;
      fhdPhidEtaPt05[i] = 0;
      fhdPhidEtaPt01[i] = 0;
      fhdPhidEtaPt001[i] = 0;
    }
  DefineOutput(1, TList::Class());
}



Bool_t AliAnalysisTaskJetCorrections::Notify()
{
  //
  // Implemented Notify() to read the cross sections
  // and number of trials from pyxsec.root
  // 
  TTree *tree = AliAnalysisManager::GetAnalysisManager()->GetTree();
  UInt_t   ntrials  = 0;
  if(tree){
    TFile *curfile = tree->GetCurrentFile();
    if (!curfile) {
      Error("Notify","No current file");
      return kFALSE;
    }

    TString fileName(curfile->GetName());
    if(fileName.Contains("AliESDs.root")){
        fileName.ReplaceAll("AliESDs.root", "pyxsec.root");
    }
    else if(fileName.Contains("AliAOD.root")){
        fileName.ReplaceAll("AliAOD.root", "pyxsec.root");
    }
    else if(fileName.Contains("galice.root")){
        // for running with galice and kinematics alone...                      
        fileName.ReplaceAll("galice.root", "pyxsec.root");
    }
    TFile *fxsec = TFile::Open(fileName.Data());
    if(!fxsec){
      Printf("%s:%d %s not found in the Input",(char*)__FILE__,__LINE__,fileName.Data());
      // no a severe condition
      return kTRUE;
    }
    TTree *xtree = (TTree*)fxsec->Get("Xsection");
    if(!xtree){
      Printf("%s:%d tree not found in the pyxsec.root",(char*)__FILE__,__LINE__);
    }
    xtree->SetBranchAddress("xsection",&fXsection);
    xtree->SetBranchAddress("ntrials",&ntrials);
    xtree->GetEntry(0);
  }
  
  return kTRUE;
}


//___________________________________________________________________________________________________________________________________
void AliAnalysisTaskJetCorrections::UserCreateOutputObjects()
{
  //
  // Create the output container
  //
  //  Printf("Analysing event  %s :: # %5d\n", gSystem->pwd(), (Int_t) fEntry);

  if(fUseAODInput){
    fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
    if(!fAOD){
      Printf("%s:%d AODEvent not found in Input Manager %d",(char*)__FILE__,__LINE__,fUseAODInput);
      return;
    }    
  }
  else{
    //  assume that the AOD is in the general output...
    fAOD  = AODEvent();
    if(!fAOD){
      Printf("%s:%d AODEvent not found in the Output",(char*)__FILE__,__LINE__);
      return;
    }    
  }

  printf("AnalysisTaskJetSpectrum::UserCreateOutputObjects() \n");

  fList = new TList();

  fhEGen = new TH1F("EGen", "", 100, 0, 200);
  fhEGen->Sumw2();
  fList->Add(fhEGen);

  fhERec = new TH1F("ERec", "", 100, 0, 200);
  fhERec->Sumw2();
  fList->Add(fhERec);

  fhEGenRest = new TH1F("EGenRest", "", 100, 0, 200);
  fhEGenRest->Sumw2();
  fList->Add(fhEGenRest);

  fhERecRest = new TH1F("ERecRest", "", 100, 0, 200);
  fhERecRest->Sumw2();
  fList->Add(fhERecRest);

  fhEsumGenRest = new TH1F("EsumGenRest", "", 100, 0, 200);
  fhEsumGenRest->Sumw2();
  fList->Add(fhEsumGenRest);

  fhEsumRecRest = new TH1F("EsumRecRest", "", 100, 0, 200);
  fhEsumRecRest->Sumw2();
  fList->Add(fhEsumRecRest);

  fhE2vsE1Gen = new TH2F("E2vsE1Gen", "", 100, 0, 200, 100, 0, 200);
  fhE2vsE1Gen->Sumw2();
  fList->Add(fhE2vsE1Gen);

  fhE2vsE1Rec = new TH2F("E2vsE1Rec", "", 100, 0, 200, 100, 0, 200);
  fhE2vsE1Rec->Sumw2();
  fList->Add(fhE2vsE1Rec);

  fhE2E1vsEsumGen = new TH2F("E2E1vsEsumGen", "", 100, 0, 200, 25, 0, 1);
  fhE2E1vsEsumGen->Sumw2();
  fList->Add(fhE2E1vsEsumGen);

  fhE2E1vsEsumRec = new TH2F("E2E1vsEsumRec", "", 100, 0, 200, 25, 0, 1);
  fhE2E1vsEsumRec->Sumw2();
  fList->Add(fhE2E1vsEsumRec);

  fhE2E1vsE1Gen = new TH2F("E2E1vsE1Gen", "", 100, 0, 200, 25, 0, 1);
  fhE2E1vsE1Gen->Sumw2();
  fList->Add(fhE2E1vsE1Gen);

  fhE2E1vsE1Rec = new TH2F("E2E1vsE1Rec", "", 100, 0, 200, 25, 0, 1);
  fhE2E1vsE1Rec->Sumw2();
  fList->Add(fhE2E1vsE1Rec);

  fhE2E1vsdPhiGen =  new TH2F("E2E1vsdPhiGen", "", 64, -3.20, 3.20, 25, 0, 1);
  fList->Add(fhE2E1vsdPhiGen);

  fhE2E1vsdPhiRec =  new TH2F("E2E1vsdPhiRec", "", 64, -3.20, 3.20, 25, 0, 1);
  fList->Add(fhE2E1vsdPhiRec);
  
  fhTrackBalance2 = new TH2F("TrackBalance2", "", 60, 0, 30, 60, 0, 30);
  fhTrackBalance2->Sumw2();
  fList->Add(fhTrackBalance2);
  
  fhTrackBalance3 = new TH2F("TrackBalance3", "", 60, 0, 30, 60, 0, 30);
  fhTrackBalance3->Sumw2();
  fList->Add(fhTrackBalance3);

  fhEt1Et22 = new TH2F("Et1Et22", "", 100, 0, 50, 100, 0, 50);
  fhEt1Et22->Sumw2();
  fList->Add(fhEt1Et22);

  fhEt1Et23 = new TH2F("Et1Et23", "", 100, 0, 50, 100, 0, 50);
  fhEt1Et23->Sumw2();
  fList->Add(fhEt1Et23);

  for(Int_t i = 0; i < 3; i++)
    {
      fhECorrJet10[i] = new TProfile(Form("ECorrJet10%d", i+1), "", 100, 0, 200, 0, 10);
      fhECorrJet10[i]->SetXTitle("E_{rec} [GeV]");
      fhECorrJet10[i]->SetYTitle("C=E_{gen}/E_{rec}");
      fhECorrJet10[i]->Sumw2();

      fhECorrJet05[i] = new TProfile(Form("ECorrJet05%d", i+1), "", 100, 0, 200, 0, 10);
      fhECorrJet05[i]->SetXTitle("E_{rec} [GeV]");
      fhECorrJet05[i]->SetYTitle("C=E_{gen}/E_{rec}");
      fhECorrJet05[i]->Sumw2();

      fhECorrJet01[i] = new TProfile(Form("ECorrJet01%d", i+1), "", 100, 0, 200, 0, 10);
      fhECorrJet01[i]->SetXTitle("E_{rec} [GeV]");
      fhECorrJet01[i]->SetYTitle("C=E_{gen}/E_{rec}");
      fhECorrJet01[i]->Sumw2();

      fhECorrJet001[i] = new TProfile(Form("ECorrJet001%d", i+1), "", 100, 0, 200, 0, 10);
      fhECorrJet001[i]->SetXTitle("E_{rec} [GeV]");
      fhECorrJet001[i]->SetYTitle("C=E_{gen}/E_{rec}");
      fhECorrJet001[i]->Sumw2();

      fhdEvsErec10[i] = new TProfile(Form("dEvsErec10_%d", i+1),"", 100, 0, 200, -1, 10);
      fhdEvsErec10[i]->SetYTitle("|E_{rec}-E_{rec}|/E_{rec}");
      fhdEvsErec10[i]->SetXTitle("E_{rec} [GeV]");
      fhdEvsErec10[i]->Sumw2();

      fhdEvsErec05[i] = new TProfile(Form("dEvsErec05_%d", i+1),"", 100, 0, 200, -1, 10);
      fhdEvsErec05[i]->SetYTitle("|E_{rec}-E_{rec}|/E_{rec}");
      fhdEvsErec05[i]->SetXTitle("E_{rec} [GeV]");
      fhdEvsErec05[i]->Sumw2();

      fhdEvsErec01[i] = new TProfile(Form("dEvsErec01_%d", i+1),"", 100, 0, 200, -1, 10);
      fhdEvsErec01[i]->SetYTitle("|E_{rec}-E_{rec}|/E_{rec}");
      fhdEvsErec01[i]->SetXTitle("E_{rec} [GeV]");
      fhdEvsErec01[i]->Sumw2();

      fhdEvsErec001[i] = new TProfile(Form("dEvsErec001_%d", i+1),"", 100, 0, 200, -1, 10);
      fhdEvsErec001[i]->SetYTitle("|E_{rec}-E_{rec}|/E_{rec}");
      fhdEvsErec001[i]->SetXTitle("E_{rec} [GeV]");
      fhdEvsErec001[i]->Sumw2();

      fhdPhidEta10[i] = new TH2F(Form("dPhidEta10_%d", i+1), "", 63, (-1)*TMath::Pi(), TMath::Pi(), 18, -0.9, 0.9);
      fhdPhidEta10[i]->SetXTitle("#phi [rad]");
      fhdPhidEta10[i]->SetYTitle("#eta");
      fhdPhidEta10[i]->Sumw2();

      fhdPhidEta05[i] = new TH2F(Form("dPhidEta05_%d", i+1), "", 63, (-1)*TMath::Pi(), TMath::Pi(), 18, -0.9, 0.9);
      fhdPhidEta05[i]->SetXTitle("#phi [rad]");
      fhdPhidEta05[i]->SetYTitle("#eta");
      fhdPhidEta05[i]->Sumw2();

      fhdPhidEta01[i] = new TH2F(Form("dPhidEta01_%d", i+1), "", 63, (-1)*TMath::Pi(), TMath::Pi(), 18, -0.9, 0.9);
      fhdPhidEta01[i]->SetXTitle("#phi [rad]");
      fhdPhidEta01[i]->SetYTitle("#eta");
      fhdPhidEta01[i]->Sumw2();

      fhdPhidEta001[i] = new TH2F(Form("dPhidEta001_%d", i+1), "", 63, (-1)*TMath::Pi(), TMath::Pi(), 18, -0.9, 0.9);
      fhdPhidEta001[i]->SetXTitle("#phi [rad]");
      fhdPhidEta001[i]->SetYTitle("#eta");
      fhdPhidEta001[i]->Sumw2();

      fhdPhidEtaPt10[i] = new TH2F(Form("dPhidEtaPt10_%d", i+1), "", 63, (-1)*TMath::Pi(), TMath::Pi(), 18, -0.9, 0.9);
      fhdPhidEtaPt10[i]->SetXTitle("#phi [rad]");
      fhdPhidEtaPt10[i]->SetYTitle("#eta");
      fhdPhidEtaPt10[i]->Sumw2();

      fhdPhidEtaPt05[i] = new TH2F(Form("dPhidEtaPt05_%d", i+1), "", 63, (-1)*TMath::Pi(), TMath::Pi(), 18, -0.9, 0.9);
      fhdPhidEtaPt05[i]->SetXTitle("#phi [rad]");
      fhdPhidEtaPt05[i]->SetYTitle("#eta");
      fhdPhidEtaPt05[i]->Sumw2();

      fhdPhidEtaPt01[i] = new TH2F(Form("dPhidEtaPt01_%d", i+1), "", 63, (-1)*TMath::Pi(), TMath::Pi(), 18, -0.9, 0.9);
      fhdPhidEtaPt01[i]->SetXTitle("#phi [rad]");
      fhdPhidEtaPt01[i]->SetYTitle("#eta");
      fhdPhidEtaPt01[i]->Sumw2();

      fhdPhidEtaPt001[i] = new TH2F(Form("dPhidEtaPt001_%d", i+1), "", 63, (-1)*TMath::Pi(), TMath::Pi(), 18, -0.9, 0.9);
      fhdPhidEtaPt001[i]->SetXTitle("#phi [rad]");
      fhdPhidEtaPt001[i]->SetYTitle("#eta");
      fhdPhidEtaPt001[i]->Sumw2();

      fList->Add(fhECorrJet10[i]);
      fList->Add(fhECorrJet05[i]);
      fList->Add(fhECorrJet01[i]);
      fList->Add(fhECorrJet001[i]);
      fList->Add(fhdEvsErec10[i]);
      fList->Add(fhdEvsErec05[i]);
      fList->Add(fhdEvsErec01[i]);
      fList->Add(fhdEvsErec001[i]);
      fList->Add(fhdPhidEta10[i]);
      fList->Add(fhdPhidEta05[i]);
      fList->Add(fhdPhidEta01[i]);
      fList->Add(fhdPhidEta001[i]);
      fList->Add(fhdPhidEtaPt10[i]);
      fList->Add(fhdPhidEtaPt05[i]);
      fList->Add(fhdPhidEtaPt01[i]);
      fList->Add(fhdPhidEtaPt001[i]);
    }
    
  Printf("UserCreateOutputObjects finished\n");
}

//__________________________________________________________________________________________________________________________________________
void AliAnalysisTaskJetCorrections::Init()
{
  printf("AliAnalysisJetCut::Init() \n");
}

//____________________________________________________________________________________________________________________________________________
void AliAnalysisTaskJetCorrections::UserExec(Option_t * )
{
//  if (fDebug > 1) printf("Analysing event # %5d\n", (Int_t) fEntry);

 
  //create an AOD handler
  AliAODHandler *aodH = dynamic_cast<AliAODHandler*>(AliAnalysisManager::GetAnalysisManager()->GetOutputEventHandler());
  
  if(!aodH)
    {
      Printf("%s:%d no output aodHandler found Jet",(char*)__FILE__,__LINE__);
      return;
    }
  
  AliMCEvent* mcEvent =MCEvent();
  if(!mcEvent){
    Printf("%s:%d no mcEvent",(char*)__FILE__,__LINE__);
    return;
  }
  
  if (fDebug > 10)Printf("%s:%d",(char*)__FILE__,__LINE__);

 //primary vertex
  AliAODVertex * pvtx = dynamic_cast<AliAODVertex*>(fAOD->GetPrimaryVertex());
  if(!pvtx) return;

 AliAODJet genJets[kMaxJets];
 Int_t nGenJets = 0;
 
 AliAODJet recJets[kMaxJets];
 Int_t nRecJets = 0;

  //array of reconstructed jets from the AOD input
 TClonesArray *aodRecJets = dynamic_cast<TClonesArray*>(fAOD->FindListObject(fBranchRec.Data()));
 if(!aodRecJets){
   Printf("%s:%d no reconstructed Jet array with name %s in AOD",(char*)__FILE__,__LINE__,fBranchRec.Data());
   return;
 }
 
 // reconstructed jets
 nRecJets = aodRecJets->GetEntries(); 
 nRecJets = TMath::Min(nRecJets, kMaxJets);
 
 for(int ir = 0;ir < nRecJets;++ir)
   {
     AliAODJet *tmp = dynamic_cast<AliAODJet*>(aodRecJets->At(ir));
     if(!tmp)continue;
     recJets[ir] = *tmp;
    }
 
 // If we set a second branch for the input jets fetch this
 TClonesArray * aodGenJets = dynamic_cast<TClonesArray*>(fAOD->FindListObject(fBranchGen.Data()));
 
 if(!aodGenJets)
   {
     printf("NO MC jets branch with name %s  Found \n",fBranchGen.Data());
     return;
   }
 
 //   //Generated jets
 nGenJets = aodGenJets->GetEntries();
 nGenJets = TMath::Min(nGenJets, kMaxJets);
 
 for(Int_t ig =0 ; ig < nGenJets; ++ig)
   {
     AliAODJet * tmp = dynamic_cast<AliAODJet*>(aodGenJets->At(ig));
     if(!tmp)continue;
     genJets[ig] = * tmp;
   }

 Double_t eRec[kMaxJets];
 Double_t eGen[kMaxJets];
 
 Double_t eRecRest[kMaxJets];
 Double_t eGenRest[kMaxJets];

// AliAODJet jetRec[kMaxJets];
 AliAODJet jetGen[kMaxJets];
 
 Int_t idxRec[kMaxJets];
 Int_t idxGen[kMaxJets];

// Double_t EsumRec = 0;
 // Double_t EsumGen =0;
 
 TLorentzVector vRec[kMaxJets];
 TLorentzVector vGen[kMaxJets];

 TLorentzVector vsumRec;
 TLorentzVector vsumGen;
 
 TVector3 pRec[kMaxJets];
 TVector3 pGen[kMaxJets];
 
 // HISTOS FOR MC
 Int_t nGenSel = 0;
 Int_t counter = 0;
 Int_t tag = 0;
 
 AliAODJet selJets[kMaxJets];
 
 // loop for applying the separation cut
 for(Int_t i = 0; i < nGenJets; i++)
   {
     if(nGenJets == 1)
       {
	 selJets[nGenSel] = genJets[i];
	 nGenSel++;
       }
     else
       {
	 counter = 0;
	 tag = 0;
	 for(Int_t j = 0; j < nGenJets; j++)
	   {
	     if(i!=j)
	       {
		 Double_t dRij = genJets[i].DeltaR(&genJets[j]);
		 counter++;
		 if(dRij > 2*fR) tag++;
	       }
	   }
	 if(counter!=0)
	   {
	     if(tag/counter == 1)
	       {
		 selJets[nGenSel] = genJets[i];
		 nGenSel++;
	       }
	   }
       }
   }
 
 for (Int_t gj = 0; gj < nGenSel; gj++)
   {
     eGen[gj] = selJets[gj].E();
     fhEGen->Fill(eGen[gj], fXsection);
   }
 
 TMath::Sort(nGenSel, eGen, idxGen);
 for (Int_t ig = 0; ig < nGenSel; ig++)
   jetGen[ig] = selJets[idxGen[ig]];
 
 //rest frame MC jets
 for (Int_t i = 0; i < nGenSel; ++i)
   {
     vGen[i].SetPxPyPzE(jetGen[i].Px(), jetGen[i].Py(), jetGen[i].Pz(), jetGen[i].E());
     pGen[i].SetXYZ(vGen[i].Px(), vGen[i].Py(), vGen[i].Pz());
     vsumGen += vGen[i];
   }
 
 if(nGenSel > 1 && pGen[0].DeltaPhi(pGen[1]) > 2.8)
   {
     fhE2vsE1Gen->Fill(jetGen[0].E(), jetGen[1].E(), fXsection);
     fhE2E1vsEsumGen->Fill(jetGen[0].E()+jetGen[1].E(), TMath::Abs(jetGen[0].E()-jetGen[1].E())/jetGen[0].E(), fXsection);
     fhE2E1vsE1Gen->Fill(jetGen[0].E(),  TMath::Abs(jetGen[0].E()-jetGen[1].E())/jetGen[0].E(), fXsection);
     Double_t deltaPhi = (jetGen[0].Phi()-jetGen[1].Phi());
     if(deltaPhi > TMath::Pi()) deltaPhi = deltaPhi - 2.*TMath::Pi();
     if(deltaPhi < (-1.*TMath::Pi())) deltaPhi = deltaPhi + 2.*TMath::Pi();
     fhE2E1vsdPhiGen->Fill(deltaPhi,  TMath::Abs(jetGen[0].E()-jetGen[1].E())/jetGen[0].E(), fXsection);
   }
     
 Double_t fPxGen = vsumGen.Px();
 Double_t fPyGen = vsumGen.Py();
 Double_t fPzGen = vsumGen.Pz();
 Double_t fEGen = vsumGen.E();
 
 Double_t eSumGenRest = 0; 
 for (Int_t j = 0; j < nGenSel; j++)
   {
     vGen[j].Boost(-fPxGen/fEGen, -fPyGen/fEGen, -fPzGen/fEGen);
     eGenRest[j] = vGen[j].E();
     if(nGenSel > 1)
       fhEGenRest->Fill(eGenRest[j], fXsection);
     eSumGenRest += eGenRest[j];
   }
 
 if(nGenSel > 1)    
   fhEsumGenRest->Fill(eSumGenRest, fXsection);
 
 //END VARIABLES FOR MC JETS ---------------
 //   }
 
 //AOD JET VARIABLES
 Int_t nRecSel = 0;
 Int_t counter1 = 0;
 Int_t tag1 = 0;
 
 AliAODJet recSelJets[kMaxJets];
 
 for(Int_t i = 0; i < nRecJets; i++)
   {
     if(nRecJets == 1)
       {
	 recSelJets[nRecSel] = recJets[i];
	 nRecSel++;
       }
     else
       {
	 counter1 = 0;
	 tag1 = 0;
	 for(Int_t j = 0; j < nRecJets; j++)
	   {
	     if(i!=j)
	       {
		 Double_t dRij = recJets[i].DeltaR(&recJets[j]);
		 counter1++;
		 if(dRij > 2*fR) tag1++;
	       }
	   }
	 if(counter1!=0)
	   {
	     if(tag1/counter1 == 1)
	       {
		 recSelJets[nRecSel] = recJets[i];
		 nRecSel++;
	       }
	   }
       }
   } 
  
 if(nRecSel == 0) return;
 Printf("******NUMBER OF JETS AFTER DELTA R CUT : %d **********\n", nRecSel);
 //sort rec/gen jets by energy in C.M.S
 AliAODJet jetRecTmp[kMaxJets];
 Int_t nAccJets = 0;
 Double_t jetTrackPt[kTracks];
 TLorentzVector jetTrackTmp[kTracks];
 Int_t nTracks = 0;
 for (Int_t rj = 0; rj < nRecSel; rj++)
   {
     TRefArray * jetTracksAOD = dynamic_cast<TRefArray*>(recSelJets[rj].GetRefTracks());
     if(!jetTracksAOD) continue;
     if(jetTracksAOD->GetEntries() < 3) continue;
     Int_t nJetTracks = 0;
     for(Int_t j = 0; j < jetTracksAOD->GetEntries(); j++)
       {
	 AliAODTrack * track = dynamic_cast<AliAODTrack*>(jetTracksAOD->At(j));
	 if(!track) continue;
	 Double_t cv[21];
	 track->GetCovarianceXYZPxPyPz(cv);
	 if(cv[14] > 1000.) continue;
	 jetTrackPt[nTracks] = track->Pt();
	 jetTrackTmp[nTracks].SetPxPyPzE(track->Px(),track->Py(),track->Pz(),track->E());
	 nTracks++;
	 nJetTracks++;
       }
     if(nJetTracks < 4) continue;
     jetRecTmp[nAccJets] = recSelJets[rj];
     eRec[nAccJets] = recSelJets[rj].E();
     fhERec->Fill(eRec[nAccJets], fXsection);
     nAccJets++;
   }

 if(nAccJets == 0) return;
 if(nTracks == 0) return;

 Printf(" ************ Number of accepted jets : %d ************ \n", nAccJets);

 AliAODJet jetRecAcc[kMaxJets];
 TMath::Sort(nAccJets, eRec, idxRec);
 for (Int_t rj = 0; rj < nAccJets; rj++)
   jetRecAcc[rj] = jetRecTmp[idxRec[rj]];
 
 //rest frame for reconstructed jets
 for (Int_t i = 0; i < nAccJets; i++)
   {
     vRec[i].SetPxPyPzE(jetRecAcc[i].Px(), jetRecAcc[i].Py(), jetRecAcc[i].Pz(), jetRecAcc[i].E());
     pRec[i].SetXYZ(vRec[i].Px(), vRec[i].Py(), vRec[i].Pz());
     vsumRec += vRec[i];
   }

 //check balance of two leading hadrons, deltaPhi > 2.
 Int_t idxTrack[kTracks];
 TMath::Sort(nTracks, jetTrackPt, idxTrack);

 TLorentzVector jetTrack[kTracks];
 for(Int_t iTr = 0; iTr < nTracks; iTr++)
   jetTrack[iTr] = jetTrackTmp[idxTrack[iTr]];
 
 Int_t n = 1;
 while(jetTrack[0].DeltaPhi(jetTrack[n]) < 2.8)
   n++;

 Double_t et1 = 0;
 Double_t et2 = 0;
 for(Int_t iTr = 0; iTr < nTracks; iTr++)
   {
     if(TMath::Abs(jetTrack[0].DeltaPhi(jetTrack[iTr]) < 1.) && iTr != 0)
       et1 += jetTrack[iTr].Et();

     if(TMath::Abs(jetTrack[n].DeltaPhi(jetTrack[iTr]) < 1.) && iTr != n)
       et2 += jetTrack[iTr].Et();
   }

 if(nAccJets == 2)
   {
     fhTrackBalance2->Fill(jetTrack[0].Et(), jetTrack[n].Et());
     fhEt1Et22->Fill(et1, et2);
   }
 if(nAccJets == 3)
   {
     fhTrackBalance3->Fill(jetTrack[0].Et(), jetTrack[n].Et());
     fhEt1Et23->Fill(et1, et2);
   }
    
 if(nAccJets > 1 && pRec[0].DeltaPhi(pRec[1]) > 2.8)
   {
     fhE2vsE1Rec->Fill(jetRecAcc[0].E(), jetRecAcc[1].E(), fXsection);
     fhE2E1vsEsumRec->Fill(jetRecAcc[0].E()+jetRecAcc[1].E(), TMath::Abs(jetRecAcc[0].E()-jetRecAcc[1].E())/jetRecAcc[0].E(), fXsection);
     fhE2E1vsE1Rec->Fill(jetRecAcc[0].E(), TMath::Abs(jetRecAcc[0].E()-jetRecAcc[1].E())/jetRecAcc[0].E(), fXsection);
     Double_t deltaPhi = (jetRecAcc[0].Phi()-jetRecAcc[1].Phi());
     if(deltaPhi > TMath::Pi()) deltaPhi = deltaPhi - 2.*TMath::Pi();
     if(deltaPhi < (-1.*TMath::Pi())) deltaPhi = deltaPhi + 2.*TMath::Pi();
     fhE2E1vsdPhiRec->Fill(-1*deltaPhi, TMath::Abs(jetRecAcc[0].E()-jetRecAcc[1].E())/jetRecAcc[0].E(), fXsection);
   }
 
 Double_t fPx = vsumRec.Px();
 Double_t fPy = vsumRec.Py();
 Double_t fPz = vsumRec.Pz();
 Double_t fE = vsumRec.E();

 Double_t eSumRecRest = 0;
 for (Int_t j = 0; j < nAccJets; j++)
   {
     vRec[j].Boost(-fPx/fE, -fPy/fE, -fPz/fE);
     eRecRest[j] = vRec[j].E();
     if(nAccJets > 1)
       fhERecRest->Fill(eRecRest[j], fXsection);
     eSumRecRest += eRecRest[j];
   }
 if(nAccJets > 1)
   fhEsumRecRest->Fill(eSumRecRest, fXsection);
 
 // Relate the jets
 Int_t iGenIndex[kMaxJets];    // Index of the generated jet for i-th rec -1 if none
 Int_t iRecIndex[kMaxJets];    // Index of the rec jet for i-th gen -1 if none
 
 for(int i = 0;i<kMaxJets;++i){
   iGenIndex[i] = iRecIndex[i] = -1;
 }
 Double_t dR = 1.4;
 if(nAccJets == nGenSel)
   {
    AliAnalysisHelperJetTasks::GetClosestJets(jetGen,nGenSel,jetRecAcc,nAccJets,
						 iGenIndex,iRecIndex,0);
      
     for(int ir = 0;ir < nAccJets;ir++){
       Int_t ig = iGenIndex[ir];
       if(ig>=0&&ig<nGenSel){
	 Double_t dPhi = TMath::Abs(jetRecAcc[ir].Phi()-jetGen[ig].Phi());
	 if(dPhi > TMath::Pi()) dPhi = dPhi - 2.*TMath::Pi();
	 if(dPhi < (-1.*TMath::Pi())) dPhi = dPhi + 2.*TMath::Pi();
	 Double_t sigma = TMath::Abs(jetGen[ig].E()-jetRecAcc[ir].E())/jetGen[ig].E();
	 dR = jetRecAcc[ir].DeltaR(&jetGen[ig]);
	 if(dR < 2*fR && dR >= fR) 
	   {
	     switch(nAccJets)
	       {
	       case 1:
		 {
		   fhdEvsErec10[0]->Fill(jetRecAcc[ir].E(), sigma, fXsection);
		   fhECorrJet10[0]->Fill(jetRecAcc[ir].E(), jetGen[ig].E()/jetRecAcc[ir].E(), fXsection);
		   for(Int_t iTr = 0; iTr < nTracks; iTr++)
		     {
		       fhdPhidEtaPt10[0]->Fill(jetTrack[iTr].Phi(), jetTrack[iTr].Eta(), jetTrack[iTr].Pt());
		       fhdPhidEta10[0]->Fill(jetTrack[iTr].Phi(), jetTrack[ir].Eta());
		     }
		 }
		 break;
	       case 2:
		 {
		   fhdEvsErec10[1]->Fill(jetRecAcc[ir].E(), sigma, fXsection);
		   fhECorrJet10[1]->Fill(jetRecAcc[ir].E(), jetGen[ig].E()/jetRecAcc[ir].E(), fXsection);
		   for(Int_t iTr = 0; iTr < nTracks; iTr++)
		     {
		       fhdPhidEtaPt10[1]->Fill(jetTrack[iTr].Phi(), jetTrack[iTr].Eta(), jetTrack[iTr].Pt());
		       fhdPhidEta10[1]->Fill(jetTrack[iTr].Phi(), jetTrack[ir].Eta());
		     }
		 }
		 break;
	       case 3:
		 {
		   fhdEvsErec10[2]->Fill(jetRecAcc[ir].E(), sigma, fXsection);
		   fhECorrJet10[2]->Fill(jetRecAcc[ir].E(), jetGen[ig].E()/jetRecAcc[ir].E(), fXsection);
		   for(Int_t iTr = 0; iTr < nTracks; iTr++)
		     {
		       fhdPhidEtaPt10[2]->Fill(jetTrack[iTr].Phi(), jetTrack[iTr].Eta(), jetTrack[iTr].Pt());
		       fhdPhidEta10[2]->Fill(jetTrack[iTr].Phi(), jetTrack[ir].Eta());
		     }
		 }
		 break;
	       }
	   }
	 if(dR < fR && dR >= 0.1) 
	   {
	     switch(nAccJets)
	       {
	       case 1:
		 {
		   fhdEvsErec05[0]->Fill(jetRecAcc[ir].E(), sigma, fXsection);
		   fhECorrJet05[0]->Fill(jetRecAcc[ir].E(), jetGen[ig].E()/jetRecAcc[ir].E(), fXsection);
		   for(Int_t iTr = 0; iTr < nTracks; iTr++)
		     {
		       fhdPhidEtaPt05[0]->Fill(jetTrack[iTr].Phi(), jetTrack[iTr].Eta(), jetTrack[iTr].Pt());
		       fhdPhidEta05[0]->Fill(jetTrack[iTr].Phi(), jetTrack[ir].Eta());
		     }
		 }
		 break;
	       case 2:
		 {
		   fhdEvsErec05[1]->Fill(jetRecAcc[ir].E(), sigma, fXsection);
		   fhECorrJet05[1]->Fill(jetRecAcc[ir].E(), jetGen[ig].E()/jetRecAcc[ir].E(), fXsection);
		   for(Int_t iTr = 0; iTr < nTracks; iTr++)
		     {
		       fhdPhidEtaPt05[1]->Fill(jetTrack[iTr].Phi(), jetTrack[iTr].Eta(), jetTrack[iTr].Pt());
		       fhdPhidEta05[1]->Fill(jetTrack[iTr].Phi(), jetTrack[ir].Eta());
		     }
		 }
		 break;
	       case 3:
		 {
		   fhdEvsErec05[2]->Fill(jetRecAcc[ir].E(), sigma, fXsection);
		   fhECorrJet05[2]->Fill(jetRecAcc[ir].E(), jetGen[ig].E()/jetRecAcc[ir].E(), fXsection);
		   for(Int_t iTr = 0; iTr < nTracks; iTr++)
		     {
		       fhdPhidEtaPt05[2]->Fill(jetTrack[iTr].Phi(), jetTrack[iTr].Eta(), jetTrack[iTr].Pt());
		       fhdPhidEta05[2]->Fill(jetTrack[iTr].Phi(), jetTrack[ir].Eta());
		     }
		 }
		 break;
	       }
	   }
	 if(dR < 0.1 && dR >= 0.01)
	   {
	     switch(nAccJets)
	       {
	       case 1:
		 {
		   fhdEvsErec01[0]->Fill(jetRecAcc[ir].E(), sigma, fXsection);
		   fhECorrJet01[0]->Fill(jetRecAcc[ir].E(), jetGen[ig].E()/jetRecAcc[ir].E(), fXsection);
		   for(Int_t iTr = 0; iTr < nTracks; iTr++)
		     {
		       fhdPhidEtaPt01[0]->Fill(jetTrack[iTr].Phi(), jetTrack[iTr].Eta(), jetTrack[iTr].Pt());
		       fhdPhidEta01[0]->Fill(jetTrack[iTr].Phi(), jetTrack[ir].Eta());
		     }
		 }
		 break;
	       case 2:
		 {
		   fhdEvsErec01[1]->Fill(jetRecAcc[ir].E(), sigma, fXsection);
		   fhECorrJet01[1]->Fill(jetRecAcc[ir].E(), jetGen[ig].E()/jetRecAcc[ir].E(), fXsection);
		   for(Int_t iTr = 0; iTr < nTracks; iTr++)
		     {
		       fhdPhidEtaPt01[1]->Fill(jetTrack[iTr].Phi(), jetTrack[iTr].Eta(), jetTrack[iTr].Pt());
		       fhdPhidEta01[1]->Fill(jetTrack[iTr].Phi(), jetTrack[ir].Eta());
		     }
		 }
		 break;
	       case 3:
		 {
		   fhdEvsErec01[2]->Fill(jetRecAcc[ir].E(), sigma, fXsection);
		   fhECorrJet01[2]->Fill(jetRecAcc[ir].E(), jetGen[ig].E()/jetRecAcc[ir].E(), fXsection);
		   for(Int_t iTr = 0; iTr < nTracks; iTr++)
		     {
		       fhdPhidEtaPt01[2]->Fill(jetTrack[iTr].Phi(), jetTrack[iTr].Eta(), jetTrack[iTr].Pt());
		       fhdPhidEta01[2]->Fill(jetTrack[iTr].Phi(), jetTrack[ir].Eta());
		     }
		 }
		 break;
	       }
	   }
	 if(dR > 0.01) continue;
	 switch(nAccJets)
	   {
	   case 1:
	     {
	       fhECorrJet001[0]->Fill(jetRecAcc[ir].E(), jetGen[ig].E()/jetRecAcc[ir].E(), fXsection);
	       fhdEvsErec001[0]->Fill(jetRecAcc[ir].E(), sigma, fXsection);
	       for(Int_t iTr = 0; iTr < nTracks; iTr++)
		 {
		   fhdPhidEtaPt001[0]->Fill(jetTrack[iTr].Phi(), jetTrack[iTr].Eta(), jetTrack[iTr].Pt());
		   fhdPhidEta001[0]->Fill(jetTrack[iTr].Phi(), jetTrack[ir].Eta());
		 }
	     }
	     break;
	   case 2:
	     {
	       fhECorrJet001[1]->Fill(jetRecAcc[ir].E(), jetGen[ig].E()/jetRecAcc[ir].E(), fXsection);
	       fhdEvsErec001[1]->Fill(jetRecAcc[ir].E(), sigma, fXsection);
	       for(Int_t iTr = 0; iTr < nTracks; iTr++)
		 {
		   fhdPhidEtaPt001[1]->Fill(jetTrack[iTr].Phi(), jetTrack[iTr].Eta(), jetTrack[iTr].Pt());
		   fhdPhidEta001[1]->Fill(jetTrack[iTr].Phi(), jetTrack[ir].Eta());
		 }
	     }
	     break;
	   case 3:
	     {
	       fhECorrJet001[2]->Fill(jetRecAcc[ir].E(), jetGen[ig].E()/jetRecAcc[ir].E(), fXsection);
	       fhdEvsErec001[2]->Fill(jetRecAcc[ir].E(), sigma, fXsection);
	       for(Int_t iTr = 0; iTr < nTracks; iTr++)
		 {
		   fhdPhidEtaPt001[2]->Fill(jetTrack[iTr].Phi(), jetTrack[iTr].Eta(), jetTrack[iTr].Pt());
		   fhdPhidEta001[2]->Fill(jetTrack[iTr].Phi(), jetTrack[ir].Eta());
		 }
	     }
	     break;
	   }
       }
     }
     // loop over reconstructed jets
   }

 Printf("%s:%d",(char*)__FILE__,__LINE__);
 
 PostData(1, fList);
 
 Printf("%s:%d Data Posted",(char*)__FILE__,__LINE__);
  
}


//__________________________________________________________________________________________________________________________________________________
void AliAnalysisTaskJetCorrections::Terminate(Option_t *)
{
  printf("AnalysisJetCorrelation::Terminate()");

} 

//_______________________________________User defined functions_____________________________________________________________________________________

void AliAnalysisTaskJetCorrections::GetThrustAxis(TVector3 &n01, TVector3 * pTrack, const Int_t &nTracks)
{
  //
  // fetch the thrust axis
  //
  TVector3 psum;
  Double_t psum1 = 0;
  Double_t psum2 = 0;
  Double_t thrust[kTracks];
  Double_t th = -3;
  Double_t tpom = -1;
  Int_t j = 0;

  
  //  for(Int_t j = 0; j < nTracks; j++)
  while(TMath::Abs(th-tpom) > 10e-7 && j < nTracks)
    {
      th = tpom;
      psum.SetXYZ(0., 0., 0.);
      psum1 = 0;
      psum2 = 0;
      for(Int_t i = 0; i < nTracks; i++)
	{
	  psum1 += (TMath::Abs(n01.Dot(pTrack[i])));
	  psum2 += pTrack[i].Mag();
	  
	  if (n01.Dot(pTrack[i]) > 0) psum += pTrack[i];
	  if (n01.Dot(pTrack[i]) < 0) psum -= pTrack[i];
	}
      
      thrust[j] = psum1/psum2;
      
      //      if(th == thrust[j]) 
      //	break;
      
      tpom = thrust[j];
      
      n01 = psum.Unit();
      j++;
    }
}
//______________________________________________________________________________________________________

  
