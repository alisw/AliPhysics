
// ******
// Task for topological study of three jet events by sona
// *****
// *****
// *****


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

#include <cstdlib>
#include <iostream>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TList.h>
#include <TTree.h>
#include <TLorentzVector.h>
#include <TClonesArray.h>
#include <TRefArray.h>
#include <TVector3.h>
#include <TVectorD.h>
#include <TMatrixDSym.h>
#include <TMatrixDSymEigen.h>
#include <TProfile.h>

#include "AliAnalysisTaskThreeJets.h"
#include "AliAnalysisManager.h"
#include "AliAODEvent.h"
#include "AliAODVertex.h"
#include "AliAODHandler.h"
#include "AliAODTrack.h"
#include "AliAODJet.h"
#include "AliAODMCParticle.h"
//#include "AliGenPythiaEventHeader.h"
//#include "AliMCEvent.h"
//#include "AliStack.h"


#include "AliAnalysisHelperJetTasks.h"


ClassImp(AliAnalysisTaskThreeJets)

  AliAnalysisTaskThreeJets::AliAnalysisTaskThreeJets() : AliAnalysisTaskSE(),
						   
							 fAOD(0x0),
						 
							 fBranchRec(""),
							 fBranchGen(""),
						
							 fUseAODInput(kFALSE),

 							 fR(0x0),
							 fList(0x0),

							 fGlobVar(1),
							 fXsection(1),

							 fhX3X4Rec(0x0),
							 fhX3X4Gen(0x0),
 							 
						 	 fhMu35Rec(0x0),
							 fhMu34Rec(0x0),
							 fhMu45Rec(0x0),
    
							 fhMu35Gen(0x0),
							 fhMu34Gen(0x0),
							 fhMu45Gen(0x0),
  	
							 fhInOut(0x0),

							 fhThrustRec2(0x0),
							 fhThrustRec3(0x0),

							 fhThrustGen2(0x0),
							 fhThrustGen3(0x0),

							 fhCGen2(0x0),
							 fhCGen3(0x0),

							 fhSGen2(0x0),
							 fhSGen3(0x0),

							 fhAGen2(0x0),
							 fhAGen3(0x0),

							 fhCRec2(0x0),
							 fhCRec3(0x0),

							 fhSRec2(0x0),
							 fhSRec3(0x0),

							 fhARec2(0x0),
							 fhARec3(0x0),

							 fhX3(0x0),
							 fhX4(0x0),
							 fhX5(0x0),

							 fhXSec(0x0),

							 fhX3X4Rec60(0x0),
							 fhX3X4Rec60100(0x0),
							 fhX3X4Rec100(0x0),

							 fhX3X4Gen60(0x0),
							 fhX3X4Gen60100(0x0),
							 fhX3X4Gen100(0x0),

							 fhdPhiThrustGen(0x0),
							 fhdPhiThrustGenALL(0x0),

							 fhdPhiThrustRec(0x0),
							 fhdPhiThrustRecALL(0x0)

{

}

AliAnalysisTaskThreeJets::AliAnalysisTaskThreeJets(const char * name):
  AliAnalysisTaskSE(name),
  
  fAOD(0x0),
  
  fBranchRec(""),
  fBranchGen(""),
  
  fUseAODInput(kFALSE),

  fR(0x0),
  fList(0x0),

  fGlobVar(1),
  fXsection(1),

  fhX3X4Rec(0x0),
  fhX3X4Gen(0x0),
  
  fhMu35Rec(0x0),
  fhMu34Rec(0x0),
  fhMu45Rec(0x0),
  
  fhMu35Gen(0x0),
  fhMu34Gen(0x0),
  fhMu45Gen(0x0),
  
  fhInOut(0x0),
  
  fhThrustRec2(0x0),
  fhThrustRec3(0x0),
  
  fhThrustGen2(0x0),
  fhThrustGen3(0x0),
  
  fhCGen2(0x0),
  fhCGen3(0x0),
  
  fhSGen2(0x0),
  fhSGen3(0x0),

  fhAGen2(0x0),
  fhAGen3(0x0),
  
  fhCRec2(0x0),
  fhCRec3(0x0),
  
  fhSRec2(0x0),
  fhSRec3(0x0),
  
  fhARec2(0x0),
  fhARec3(0x0),
  
  fhX3(0x0),
  fhX4(0x0),
  fhX5(0x0),
  
  fhXSec(0x0),
  
  fhX3X4Rec60(0x0),
  fhX3X4Rec60100(0x0),
  fhX3X4Rec100(0x0),
  
  fhX3X4Gen60(0x0),
  fhX3X4Gen60100(0x0),
  fhX3X4Gen100(0x0),
  
  fhdPhiThrustGen(0x0),
  fhdPhiThrustGenALL(0x0),
  
  fhdPhiThrustRec(0x0),
  fhdPhiThrustRecALL(0x0)
{
  DefineOutput(1, TList::Class());
}



Bool_t AliAnalysisTaskThreeJets::Notify()
{
  //
  // Implemented Notify() to read the cross sections
  // and number of trials from pyxsec.root
  // 


  TTree *tree = AliAnalysisManager::GetAnalysisManager()->GetTree();
  Float_t xsection = 0;
  Float_t ftrials  = 1;

  Float_t fAvgTrials = 1;
  if(tree){
    TFile *curfile = tree->GetCurrentFile();
    if (!curfile) {
      Error("Notify","No current file");
      return kFALSE;
    }
    AliAnalysisHelperJetTasks::PythiaInfoFromFile(curfile->GetName(),xsection,ftrials);
    // construct a poor man average trials 
    Float_t nEntries = (Float_t)tree->GetTree()->GetEntries();
    if(ftrials>=nEntries)fAvgTrials = ftrials/nEntries; // CKB take this into account for normalisation
  }  

  if(xsection>0)fXsection  = xsection;

  return kTRUE;

}


//___________________________________________________________________________________________________________________________________
void AliAnalysisTaskThreeJets::UserCreateOutputObjects()
{
  //
  // Create the output container
  //
  //  Printf("Analysing event  %s :: # %5d\n", gSystem->pwd(), (Int_t) fEntry);

  printf("AnalysisTaskJetSpectrum::UserCreateOutputObjects() \n");

  fList = new TList();

  fhX3X4Gen = new TH2F("X3vsX4Gen", "", 22, 0.6, 1.02, 33, 0.4, 1.02);
  fhX3X4Gen->SetXTitle("X_{3}");
  fhX3X4Gen->SetYTitle("X_{4}");
  fhX3X4Gen->Sumw2();
  fList->Add(fhX3X4Gen);
     
  fhX3X4Rec = new TH2F("X3vsX4Rec", "", 22, 0.6, 1.02, 33, 0.4, 1.02);
  fhX3X4Rec->SetXTitle("X_{3}");
  fhX3X4Rec->SetYTitle("X_{4}");
  fhX3X4Rec->Sumw2();
  fList->Add(fhX3X4Rec);
     
  fhMu35Rec = new TH1F("Mu35Rec", "", 20,0.1,0.8);
  fhMu35Rec->Sumw2();
  fhMu35Rec->SetXTitle("#mu_{35}");
  fhMu35Rec->SetYTitle("#frac{dN}{d#mu_{35Rec}}");  
  fList->Add(fhMu35Rec);
	  
  fhMu34Rec = new TH1F("Mu34Rec", "", 20,0.5,1);
  fhMu34Rec->Sumw2();
  fhMu34Rec->SetXTitle("#mu_{34}");
  fhMu34Rec->SetYTitle("#frac{dN}{d#mu_{34}}");
  fList->Add(fhMu34Rec);
  
  fhMu45Rec = new TH1F("Mu45Rec", "", 20,0,0.65);
  fhMu45Rec->Sumw2();
  fhMu45Rec->SetXTitle("#mu_{45}");
  fhMu45Rec->SetYTitle("#frac{dN}{d#mu_{45}}");
  fList->Add(fhMu45Rec);
     
  fhMu35Gen = new TH1F("Mu35Gen", "", 20,0.1,0.8);
  fhMu35Gen->Sumw2();
  fhMu35Gen->SetXTitle("#mu_{35Gen}");
  fhMu35Gen->SetYTitle("#frac{dN}{d#mu_{35Gen}}");  
  fList->Add(fhMu35Gen);
	  
  fhMu34Gen = new TH1F("Mu34Gen", "", 20,0.5,1);
  fhMu34Gen->Sumw2();
  fhMu34Gen->SetXTitle("#mu_{34Gen}");
  fhMu34Gen->SetYTitle("#frac{dN}{d#mu_{34Gen}}");
  fList->Add(fhMu34Gen);
  
  fhMu45Gen = new TH1F("Mu45Gen", "", 20,0,0.65);
  fhMu45Gen->Sumw2();
  fhMu45Gen->SetXTitle("#mu_{45Gen}");
  fhMu45Gen->SetYTitle("#frac{dN}{d#mu_{45}}");
  fList->Add(fhMu45Gen);

  fhInOut = new TH1I("InOut", "", 6, 0, 6);
  fhInOut->SetXTitle("#RecJets_{GenJets=3}");
  fhInOut->SetYTitle("#Entries");
  fList->Add(fhInOut);

  fhThrustGen2 = new TH1F("ThrustGen2", "", 50, 0.5, 1);
  fhThrustGen2->Sumw2();
  fList->Add(fhThrustGen2);

  fhThrustGen3 = new TH1F("ThrustGen3", "", 50, 0.5, 1);
  fhThrustGen3->Sumw2();
  fList->Add(fhThrustGen3);

  fhThrustRec2 = new TH1F("ThrustRec2", "", 50, 0.5, 1);
  fhThrustRec2->Sumw2();
  fList->Add(fhThrustRec2);

  fhThrustRec3 = new TH1F("ThrustRec3", "", 50, 0.5, 1);
  fhThrustRec3->Sumw2();
  fList->Add(fhThrustRec3);

  fhCGen2 = new TH1F("CGen2", "", 100, 0, 1);
  fhCGen2->Sumw2();
  fList->Add(fhCGen2);

  fhCGen3 = new TH1F("CGen3", "", 100, 0, 1);
  fhCGen3->Sumw2();
  fList->Add(fhCGen3);

  fhCRec2 = new TH1F("CRec2", "", 100, 0, 1);
  fhCRec2->Sumw2();
  fList->Add(fhCRec2);

  fhCRec3 = new TH1F("CRec3", "", 100, 0, 1);
  fhCRec3->Sumw2();
  fList->Add(fhCRec3);

  fhSGen2 = new TH1F("SGen2", "", 100, 0, 1);
  fList->Add(fhSGen2);

  fhSGen3 = new TH1F("SGen3", "", 100, 0, 1);
  fList->Add(fhSGen3);

  fhSRec2 = new TH1F("SRec2", "", 100, 0, 1);
  fList->Add(fhSRec2);

  fhSRec3 = new TH1F("SRec3", "", 100, 0, 1);
  fList->Add(fhSRec3);

  fhAGen2 = new TH1F("AGen2", "", 50, 0, 0.5);
  fList->Add(fhAGen2);

  fhAGen3 = new TH1F("AGen3", "", 50, 0, 0.5);
  fList->Add(fhAGen3);

  fhARec2 = new TH1F("ARec2", "", 50, 0, 0.5);
  fList->Add(fhARec2);

  fhARec3 = new TH1F("ARec3", "", 50, 0, 0.5);
  fList->Add(fhARec3);

  fhX3 = new TH2F("X3", "", 22, 0.6, 1.02, 100, 0, 1);
  fhX3->SetYTitle("|X_{3}^{MC} - X_{3}^{AOD}|/X_{3}^{MC}");
  fhX3->SetXTitle("X_{3}");
  fhX3->Sumw2();
  fList->Add(fhX3);

  fhX4 = new TH2F("X4", "",33, 0.4, 1.02, 100, 0, 1);
  fhX4->SetYTitle("|X_{4}^{MC} - X_{4}^{AOD}|/X_{4}^{MC}");
  fhX4->SetXTitle("X_{4}");
  fhX4->Sumw2();
  fList->Add(fhX4);

  fhX5 = new TH2F("X5", "",100, 0., 1., 100, 0, 1);
  fhX5->SetYTitle("|X_{5}^{MC} - X_{5}^{AOD}|/X_{5}^{MC}");
  fhX5->SetXTitle("X_{5}");
  fhX5->Sumw2();
  fList->Add(fhX5);

  fhXSec = new TProfile("XSec", "", 200, 0, 200, 0, 1);
  fList->Add(fhXSec);

  fhX3X4Rec60 = new TH2F("X3vsX4Rec60", "", 22, 0.6, 1.02, 33, 0.4, 1.02);
  fhX3X4Rec60->SetXTitle("X_{3}");
  fhX3X4Rec60->SetYTitle("X_{4}");
  fhX3X4Rec60->Sumw2();
  fList->Add(fhX3X4Rec60);

  fhX3X4Rec60100 = new TH2F("X3vsX4Rec60100", "", 22, 0.6, 1.02, 33, 0.4, 1.02);
  fhX3X4Rec60100->SetXTitle("X_{3}");
  fhX3X4Rec60100->SetYTitle("X_{4}");
  fhX3X4Rec60100->Sumw2();
  fList->Add(fhX3X4Rec60100);

  fhX3X4Rec100 = new TH2F("X3vsX4Rec100", "", 22, 0.6, 1.02, 33, 0.4, 1.02);
  fhX3X4Rec100->SetXTitle("X_{3}");
  fhX3X4Rec100->SetYTitle("X_{4}");
  fhX3X4Rec100->Sumw2();
  fList->Add(fhX3X4Rec100);

  fhX3X4Gen60 = new TH2F("X3vsX4Gen60", "", 22, 0.6, 1.02, 33, 0.4, 1.02);
  fhX3X4Gen60->SetXTitle("X_{3}");
  fhX3X4Gen60->SetYTitle("X_{4}");
  fhX3X4Gen60->Sumw2();
  fList->Add(fhX3X4Gen60);

  fhX3X4Gen60100 = new TH2F("X3vsX4Gen60100", "", 22, 0.6, 1.02, 33, 0.4, 1.02);
  fhX3X4Gen60100->SetXTitle("X_{3}");
  fhX3X4Gen60100->SetYTitle("X_{4}");
  fhX3X4Gen60100->Sumw2();
  fList->Add(fhX3X4Gen60100);

  fhX3X4Gen100 = new TH2F("X3vsX4Gen100", "", 22, 0.6, 1.02, 33, 0.4, 1.02);
  fhX3X4Gen100->SetXTitle("X_{3}");
  fhX3X4Gen100->SetYTitle("X_{4}");
  fhX3X4Gen100->Sumw2();
  fList->Add(fhX3X4Gen100);

  fhdPhiThrustGen = new TH2F("dPhiThrustGen", "", 32, -1*TMath::Pi(), TMath::Pi(), 25, 0, 150);
  fhdPhiThrustGen->Sumw2();
  fList->Add(fhdPhiThrustGen);

  fhdPhiThrustGenALL = new TH2F("dPhiThrustGenALL", "", 32, -1*TMath::Pi(), TMath::Pi(), 25, 0,  150);
  fhdPhiThrustGenALL->Sumw2();
  fList->Add(fhdPhiThrustGenALL);

  fhdPhiThrustRec = new TH2F("dPhiThrustRec", "", 32, -1*TMath::Pi(), TMath::Pi(), 25, 0, 150);
  fhdPhiThrustRec->Sumw2();
  fList->Add(fhdPhiThrustRec);

  fhdPhiThrustRecALL = new TH2F("dPhiThrustRecALL", "", 32, -1*TMath::Pi(), TMath::Pi(), 25, 0, 150);
  fhdPhiThrustRecALL->Sumw2();
  fList->Add(fhdPhiThrustRecALL);
    
  if(fDebug)Printf("UserCreateOutputObjects finished\n");
}

//__________________________________________________________________________________________________________________________________________
void AliAnalysisTaskThreeJets::Init()
{
  printf("AliAnalysisJetCut::Init() \n");
}

//____________________________________________________________________________________________________________________________________________
void AliAnalysisTaskThreeJets::UserExec(Option_t * )
{
  if (fDebug > 1) printf("AliAnlysisTaskThreeJets::Analysing event # %5d\n", (Int_t) fEntry);

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

//   AliMCEvent* mcEvent =MCEvent();
//   if(!mcEvent){
//     Printf("%s:%d no mcEvent",(char*)__FILE__,__LINE__);
//     return;
//   }
  
  if (fDebug > 10)Printf("%s:%d",(char*)__FILE__,__LINE__);
  
  //primary vertex
  AliAODVertex * pvtx = dynamic_cast<AliAODVertex*>(fAOD->GetPrimaryVertex());
  if(!pvtx){
    if (fDebug > 1)     Printf("%s:%d AOD Vertex found",(char*)__FILE__,__LINE__);
    return;
  }
  
//   AliAODJet genJetsPythia[kMaxJets];
//   Int_t nPythiaGenJets = 0;
  
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
      printf("NO MC jets Found\n");
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
  
//   AliGenPythiaEventHeader*  pythiaGenHeader = AliAnalysisHelperJetTasks::GetPythiaEventHeader(mcEvent);
//   if(!pythiaGenHeader){
//     Printf("!!!NO GEN HEADER AVALABLE!!!");
//     return;
//   }
  
//   // Int_t ProcessType = pythiaGenHeader->ProcessType();
//   // if(ProcessType != 28) return;
//   nPythiaGenJets = pythiaGenHeader->NTriggerJets();
//   nPythiaGenJets = TMath::Min(nPythiaGenJets, kMaxJets);
   
//  fXsection = 1;

  Double_t eRec[kMaxJets];
  Double_t eGen[kMaxJets];
 
  Double_t eJetRec[kMaxJets];
  //  Double_t EJetGen[kMaxJets];

  AliAODJet jetRec[kMaxJets];
  AliAODJet jetGen[kMaxJets];

  Int_t idxRec[kMaxJets];
  Int_t idxGen[kMaxJets];
   
  Double_t xRec[kMaxJets]; 
  Double_t xGen[kMaxJets];

  Double_t eSumRec = 0;
  Double_t eSumGen = 0;
  
  TLorentzVector vRec[kMaxJets];
  TLorentzVector vRestRec[kMaxJets];
  
  TLorentzVector vGen[kMaxJets];
  TLorentzVector vRestGen[kMaxJets];

  TLorentzVector vsumRec;
  TLorentzVector vsumGen;

  TVector3 pRec[kMaxJets];
  TVector3 pGen[kMaxJets];

  TVector3 pTrack[kTracks];

  TVector3 pRestRec[kMaxJets];
  TVector3 pRestGen[kMaxJets];

  Double_t psumRestRec = 0;
  //  Double_t psumRestGen = 0;
//   //Pythia_________________________________________________________________________________________________________________

//   for(int ip = 0;ip < nPythiaGenJets;++ip)
//     {
//       if(ip>=kMaxJets)continue;
//       Float_t p[4];
//       pythiaGenHeader->TriggerJet(ip,p);
//       genJetsPythia[ip].SetPxPyPzE(p[0],p[1],p[2],p[3]);
//     }
  
//_________________________________________________________________________________________________________________________
 
 
//________histos for MC___________________________________________________________________________________________________________

  Int_t nGenSel = 0;
  Int_t counter = 0;
  Int_t tag = 0;

  AliAODJet selJets[kMaxJets];

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

  if(nGenSel == 0) return;

  for (Int_t gj = 0; gj < nGenSel; gj++)
    {
      eGen[gj] = selJets[gj].E();
    }
  
  TMath::Sort(nGenSel, eGen, idxGen);
  for (Int_t ig = 0; ig < nGenSel; ig++)
    {
      jetGen[ig] = selJets[idxGen[ig]];
    }

  fhXSec->Fill(jetGen[0].Pt(), fXsection);

//  AliStack * stack = mcEvent->Stack();

  TClonesArray *tca = dynamic_cast<TClonesArray*>(fAOD->FindListObject(AliAODMCParticle::StdBranchName()));
  if(!tca) return;
  Int_t nMCtracks = tca->GetEntries();
  Double_t * eTracksMC = new Double_t[kTracks];
  Double_t pTracksMC[kTracks];
  Int_t * idxTracksMC = new Int_t[kTracks];
  TLorentzVector jetTracksMC[kTracks];
  TLorentzVector jetTracksSortMC[kTracks];
  TVector3 pTrackMC[kTracks];
  TLorentzVector vTrackMCAll[kTracks];
  Double_t pTrackMCAll[kTracks];
  TLorentzVector vTrackMC[kTracks];
  TVector3 pTrackMCBoost[kTracks];
  Double_t eventShapes[4];

  Int_t nAccTr = 0;
  Int_t nInJet[kMaxJets];
  TLorentzVector inJetPartV[kMaxJets][kTracks];
  Int_t nAllTracksMC = 0;

  for(Int_t iTrack = 0; iTrack < nMCtracks; iTrack++)
    {
      //      TParticle * part = (TParticle*)stack->Particle(iTrack);
      AliAODMCParticle *part = dynamic_cast<AliAODMCParticle*>(tca->At(iTrack));
      if (!part) continue;
      if(!part->IsPhysicalPrimary())continue;      
      Double_t fEta = part->Eta();
      if(TMath::Abs(fEta) > .9) continue;
      vTrackMCAll[nAllTracksMC].SetPxPyPzE(part->Px(), part->Py(), part->Pz(), part->E());
      pTrackMCAll[nAllTracksMC] = part->Pt();
      nAllTracksMC++;
    }
  if(nAllTracksMC == 0) return;
  for(Int_t iJet = 0; iJet < nGenSel; iJet++)
    {
      Int_t nJetTracks = 0;
      for(Int_t i = 0; i < nAllTracksMC; i++)
	{
	  Double_t dPhi = (jetGen[iJet].Phi()-vTrackMCAll[i].Phi());
	  if(dPhi > TMath::Pi()) dPhi = dPhi - 2.*TMath::Pi();
	  if(dPhi < (-1.*TMath::Pi())) dPhi = dPhi + 2.*TMath::Pi();
	  Double_t dEta = (jetGen[iJet].Eta()-vTrackMCAll[i].Eta());
	  Double_t deltaR = TMath::Sqrt(dPhi*dPhi+dEta*dEta);
	  if(deltaR < fR && vTrackMCAll[i].Pt() > 1.5)
	    {
	      jetTracksMC[nAccTr] = vTrackMCAll[i];
	      eTracksMC[nAccTr] = vTrackMCAll[i].E();
	      pTracksMC[nAccTr] = vTrackMCAll[i].Pt();
	      inJetPartV[iJet][nJetTracks].SetPxPyPzE(vTrackMCAll[i].Px(), vTrackMCAll[i].Py(), vTrackMCAll[i].Pz(),vTrackMCAll[i].E());
	      nAccTr++;
	      nJetTracks++;
	    }
	}
      nInJet[iJet] = nJetTracks;
    }

  if(nAccTr == 0) return;
  if(fDebug)Printf("*********** Number of Jets : %d ***************\n", nGenSel);  
  Double_t pTav[kMaxJets];
  for(Int_t i = 0; i < nGenSel; i++)
    {
      Double_t pTsum = 0;
      if(fDebug)Printf("*********** Number of particles in Jet %d = %d *******************\n", i+3, nInJet[i]);
      for(Int_t iT = 0; iT < nInJet[i]; iT++)
	{
	  Double_t pt = inJetPartV[i][iT].Pt();
	  pTsum += pt;
	}
      pTav[i] = pTsum/nInJet[i]; 
    }
  
  TMath::Sort(nAllTracksMC, pTrackMCAll, idxTracksMC);
  for(Int_t i = 0; i < nAllTracksMC; i++)
    {
      jetTracksSortMC[i] = vTrackMCAll[idxTracksMC[i]];
      pTrackMC[i].SetXYZ(jetTracksSortMC[i].Px(), jetTracksSortMC[i].Py(), jetTracksSortMC[i].Pz());
      vTrackMC[i].SetPxPyPzE(jetTracksSortMC[i].Px(), jetTracksSortMC[i].Py(), jetTracksSortMC[i].Pz(), jetTracksSortMC[i].E());
    }

  TVector3 n01MC = pTrackMC[0].Unit();
  n01MC.SetZ(0.);

  //Thrust calculation, iterative method
  if(nGenSel > 1)
    { 
//       if(fGlobVar == 1)
// 	{
      if(fDebug)Printf("**************Shapes for MC*************");
      AliAnalysisHelperJetTasks::GetEventShapes(n01MC, pTrackMC, nAllTracksMC, eventShapes);
// 	}
      if(eventShapes[0] < 2/TMath::Pi()){
	Double_t eventShapesTest[4];
	TVector3 n01Test;
	Int_t rnd_max = nAllTracksMC;
	Int_t k = (rand()%rnd_max)+3;
	while(TMath::Abs(pTrackMC[k].X()) < 10e-5 && TMath::Abs(pTrackMC[k].Y()) < 10e-5){
	  k--;
	}
	n01Test = pTrackMC[k].Unit();
	n01Test.SetZ(0.);
	AliAnalysisHelperJetTasks::GetEventShapes(n01Test, pTrackMC, nAllTracksMC, eventShapesTest);
	eventShapes[0] = TMath::Max(eventShapes[0], eventShapesTest[0]);
	if(TMath::Abs(eventShapes[0]-eventShapesTest[0]) < 10e-7) n01MC = n01Test;
      }

      Double_t s = eventShapes[1];
      Double_t a = eventShapes[2];
      Double_t c = eventShapes[3];

      switch(nGenSel)
	{
	case 2:
	  {
	    fhAGen2->Fill(a);
	    fhSGen2->Fill(s);
	    fhCGen2->Fill(c);
	  }
	  break;
	case 3:
	  {
	    fhAGen3->Fill(a);
	    fhSGen3->Fill(s);
	    fhCGen3->Fill(c);
	  }
	  break;
	}
      Double_t thrust01MC = eventShapes[0];
      
      switch(nGenSel)
	{
	case 2:
	  fhThrustGen2->Fill(thrust01MC, fXsection);
	  break;
	case 3:
	  fhThrustGen3->Fill(thrust01MC, fXsection);
	  break;
	}
    }
  
  
  //rest frame MC jets
  for (Int_t i = 0; i < nGenSel; ++i)
    {
      vGen[i].SetPxPyPzE(jetGen[i].Px(), jetGen[i].Py(), jetGen[i].Pz(), jetGen[i].E());
      pGen[i].SetXYZ(vGen[i].Px(), vGen[i].Py(), vGen[i].Pz());
      vsumGen += vGen[i];
    }
  if(eventShapes[0] > 0.8 && nGenSel > 1)
    {
      for(Int_t i = 0; i < nGenSel; i++)
	fhdPhiThrustGen->Fill(n01MC.DeltaPhi(pGen[i]), jetGen[i].E());
    }
 
  if(eventShapes[0] <= 0.8 && nGenSel > 1)   
    {
      for(Int_t i = 0; i < nGenSel; i++)
	fhdPhiThrustGenALL->Fill(n01MC.DeltaPhi(pGen[i]), jetGen[i].E());
    }
      
  Double_t fPxGen = vsumGen.Px();
  Double_t fPyGen = vsumGen.Py();
  Double_t fPzGen = vsumGen.Pz();
  Double_t fEGen = vsumGen.E();

  Double_t eRestGen[kMaxJets];
  for (Int_t j = 0; j < nGenSel; j++)
    {
      vGen[j].Boost(-fPxGen/fEGen, -fPyGen/fEGen, -fPzGen/fEGen);
      eRestGen[j] = vGen[j].E();
    }

  for (Int_t j = 0; j < nAccTr; j++)
    {
      vTrackMC[j].Boost(-fPxGen/fEGen, -fPyGen/fEGen, -fPzGen/fEGen);
      pTrackMCBoost[j].SetXYZ(vTrackMC[j].Px(),vTrackMC[j].Py(),vTrackMC[j].Pz());
    }
      
  Int_t idxRestGen[kMaxJets];
  TMath::Sort(nGenSel, eRestGen, idxRestGen);
  for(Int_t j = 0; j < nGenSel; j++)
    {
      vRestGen[j] = vGen[idxRestGen[j]];
      eSumGen += vRestGen[j].E();
    }

  if (nGenSel == 3)
    {
      //      if(nInJet[0] < 3 || nInJet[1] < 3 || nInJet[2] < 3) return;
      //       if(pRestGen[1].DeltaPhi(pRestGen[2]) > 0.95 && pRestGen[1].DeltaPhi(pRestGen[2]) < 1.15)
      // 	{
      
      for(Int_t i = 0; i < nGenSel; i++)
	{
	  xGen[i] = 2*vRestGen[i].E()/eSumGen;
	}
      
      if(fDebug)Printf("***************** Values of Dalitz variables are : %f, %f, %f ****************\n", xGen[0], xGen[1], xGen[2]);
      
      if(fDebug)Printf("***************** fXSection = %f ******************\n", fXsection);
      if(eSumGen <= 60)
	fhX3X4Gen60->Fill(xGen[0], xGen[1], fXsection);
      
      if(eSumGen > 60 && eSumGen <= 100)
	fhX3X4Gen60100->Fill(xGen[0], xGen[1], fXsection);
      
      if(eSumGen > 100)
	fhX3X4Gen100->Fill(xGen[0], xGen[1], fXsection);
	  
      FillTopology(fhX3X4Gen, fhMu34Gen, fhMu45Gen, fhMu35Gen, xGen, pRestGen, fXsection);
    }
  


//_______________________________________________histos for MC_____________________________________________________


//_______________________________________________histos AOD________________________________________________________

//   Printf("Event Number : %d, Number of gen jets : %d ", fEntry, nGenJets);

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
  
  //sort rec/gen jets by energy in C.M.S
  for (Int_t rj = 0; rj < nRecSel; rj++)
    {
      eRec[rj] = recSelJets[rj].E();
    }
 
  Int_t nAODtracks = fAOD->GetNumberOfTracks();
  Int_t nTracks = 0; //tracks accepted in the whole event
  Int_t nTracksALL = 0;
  TLorentzVector jetTracks[kTracks];
  TLorentzVector jetTracksSort[kTracks];
  Double_t * eTracks = new Double_t[kTracks];
  Double_t pTracks[kTracks];
  Int_t * idxTracks = new Int_t[kTracks];
  Double_t eventShapesRec[4];
  Int_t jetMult[kMaxJets];
  //  TLorentzVector vTracksAll[kTracks];
  //  Double_t pTracksAll[kTracks];
  Int_t nAccJets = 0;
  AliAODJet jetRecAcc[kMaxJets];
  Int_t nJetTracks = 0;

  AliAODTrack jetTrack[kTracks];
  Double_t * cv = new Double_t[21];
  TMath::Sort(nRecSel, eRec, idxRec);
  for (Int_t rj = 0; rj < nRecSel; rj++)
    {
      nJetTracks = 0;
      eJetRec[rj] = eRec[idxRec[rj]];
      jetRec[rj] = recSelJets[idxRec[rj]];
      TRefArray * jetTracksAOD = dynamic_cast<TRefArray*>(jetRec[rj].GetRefTracks());
      if(!jetTracksAOD) continue;
      if(jetTracksAOD->GetEntries() < 3) continue;
      for(Int_t i = 0; i < jetTracksAOD->GetEntries(); i++)
	{
	  AliAODTrack * track = (AliAODTrack*)jetTracksAOD->At(i);
	  track->GetCovarianceXYZPxPyPz(cv);
	  if(cv[14] > 1000.) continue;
	  //	  jetTrack[nTracks] = *track;
	  //	  jetTracks[nTracks].SetPxPyPzE(jetTrack[nTracks].Px(), jetTrack[nTracks].Py(), jetTrack[nTracks].Pz(), jetTrack[nTracks].E());
	  //	  eTracks[nTracks] = jetTracks[nTracks].E();
	  //	  pTracks[nTracks] = jetTracks[nTracks].Pt();
	  nTracks++;
	  nJetTracks++;
	}
      if(nJetTracks < 3) continue;
      jetRecAcc[nAccJets] = jetRec[rj];
      jetMult[nAccJets] = jetTracksAOD->GetEntries();
      nAccJets++;
    }
 
  if (nAccJets == 0) return;

  for(Int_t i = 0; i < nAODtracks; i++)
    {
      AliAODTrack * track = dynamic_cast<AliAODTrack*>(fAOD->GetTrack(i));
      if(!track) continue;
      track->GetCovarianceXYZPxPyPz(cv);
      if(cv[14] > 1000.) continue;
      jetTrack[nTracksALL] = *track;
      jetTracks[nTracksALL].SetPxPyPzE(jetTrack[nTracksALL].Px(), jetTrack[nTracksALL].Py(), jetTrack[nTracksALL].Pz(), jetTrack[nTracksALL].E());
      eTracks[nTracksALL] = jetTracks[nTracksALL].E();
      pTracks[nTracksALL] = jetTracks[nTracksALL].Pt();
      nTracksALL++;
    }
      

  TLorentzVector vTrack[kTracks];
  TMath::Sort(nTracksALL, pTracks, idxTracks);
  for(Int_t i = 0; i < nTracksALL; i++)
    {
      jetTracksSort[i] = jetTracks[idxTracks[i]];
      pTrack[i].SetXYZ(jetTracksSort[i].Px(), jetTracksSort[i].Py(), jetTracksSort[i].Pz());
      vTrack[i].SetPxPyPzE(jetTracksSort[i].Px(), jetTracksSort[i].Py(), jetTracksSort[i].Pz(), jetTracksSort[i].E());
    }

  for (Int_t i = 0; i < nAccJets; ++i)
    {
      vRec[i].SetPxPyPzE(jetRecAcc[i].Px(), jetRecAcc[i].Py(), jetRecAcc[i].Pz(), jetRecAcc[i].E());
      pRec[i].SetXYZ(vRec[i].Px(), vRec[i].Py(), vRec[i].Pz());
      vsumRec += vRec[i];
    }

  //Thrust, iterative method, AODs
  TVector3 n01 = pTrack[0].Unit();
  n01.SetZ(0.);
  if(nAccJets > 1)
    {
//       if(fGlobVar == 1)
// 	{
      if(fDebug)Printf("*********Shapes for AOD********");
      AliAnalysisHelperJetTasks::GetEventShapes(n01, pTrack, nTracksALL, eventShapesRec);
      // 	}
      if(eventShapesRec[0] < 2/TMath::Pi()){
	Double_t eventShapesTest[4];
	TVector3 n01Test;
	Int_t rnd_max = nTracksALL;
	Int_t k = (rand()%rnd_max)+3;
	while(TMath::Abs(pTrack[k].X()) < 10e-5 && TMath::Abs(pTrack[k].Y()) < 10e-5){
	  k--;
	}
	
	n01Test = pTrack[k].Unit();
	n01Test.SetZ(0.);
	AliAnalysisHelperJetTasks::GetEventShapes(n01Test, pTrack, nTracksALL, eventShapesTest);
	eventShapesRec[0] = TMath::Max(eventShapesRec[0], eventShapesTest[0]);
	if(TMath::Abs(eventShapesRec[0]-eventShapesTest[0]) < 10e-7) n01 = n01Test;
      }
    
	
      //       fGlobVar = 0;      
      //       Double_t Max3 = TMath::Max(eventShapesRec0[0],eventShapesRec1[0]);
      //       Double_t Max4 = TMath::Max(eventShapesRec3[0],eventShapesRec2[0]);
      
      Double_t thrust = eventShapesRec[0];
      
      if(eventShapesRec[0] > 0.8)
	{
	  for(Int_t i = 0; i < nAccJets; i++)
	    fhdPhiThrustRec->Fill(n01.DeltaPhi(pRec[i]), jetRecAcc[i].E());
	  
	}
      
      if(eventShapesRec[0] <= 0.8)
	{
	  for(Int_t i = 0; i < nAccJets; i++)
	    fhdPhiThrustRecALL->Fill(n01.DeltaPhi(pRec[i]), jetRecAcc[i].E());
	}

      switch(nAccJets)
	{
	case 2:
	  fhThrustRec2->Fill(thrust, fXsection);
	  break;
	case 3:
	  fhThrustRec3->Fill(thrust, fXsection);
	  break;
	}
      
      switch(nAccJets)
	{
	case 2:
	  {
	    fhARec2->Fill(eventShapesRec[2], fXsection);
	    fhSRec2->Fill(eventShapesRec[1], fXsection);
	    fhCRec2->Fill(eventShapesRec[3], fXsection);
	  }
	  break;
	case 3:
	  {
	    fhARec3->Fill(eventShapesRec[2], fXsection);
	    fhSRec3->Fill(eventShapesRec[1], fXsection);
	    fhCRec3->Fill(eventShapesRec[3], fXsection);
	  }
	  break;
	}
    }

  //rest frame for reconstructed jets
  Double_t fPx = vsumRec.Px();
  Double_t fPy = vsumRec.Py();
  Double_t fPz = vsumRec.Pz();
  Double_t fE = vsumRec.E();
  	  	  
  TVector3 pTrackRest[kTracks];
  for(Int_t j = 0; j < nTracks; j++)
    {
      vTrack[j].Boost(-fPx/fE, -fPy/fE, -fPz/fE);
      pTrackRest[j].SetXYZ(vTrack[j].Px(), vTrack[j].Py(),vTrack[j].Pz());
    }
      
  Double_t eRestRec[kMaxJets];
  Int_t idxRestRec[kMaxJets];
  for (Int_t j = 0; j < nAccJets; j++)
    {
      vRec[j].Boost(-fPx/fE, -fPy/fE, -fPz/fE);
      eRestRec[j] = vRec[j].E();  
      eSumRec += vRec[j].E();
    }

  TMath::Sort(nAccJets, eRestRec, idxRestRec);
  for (Int_t i = 0; i < nAccJets; i++)
    {
      vRestRec[i] = vRec[idxRestRec[i]]; 
      pRestRec[i].SetXYZ(vRestRec[i].Px(), vRestRec[i].Py(), vRestRec[i].Pz());
      psumRestRec += pRestRec[i].Perp();
    } 

  if(nAccJets == 3)
    {
//       if(pRest[1].DeltaPhi(pRest[2]) > 0.95 && pRest[1].DeltaPhi(pRest[2]) < 1.15)
//  	{
	  fhInOut->Fill(nGenSel);
// 	  for(Int_t j = 0; j < nTracksALL; j++)
// 	    {
// 	      vTracksAll[j].Boost(-fPx/fE, -fPy/fE, -fPz/fE);
// 	      pTracksAll[j].SetXYZ(vTracksAll[j].Px(), vTracksAll[j].Py(),vTracksAll[j].Pz());
// 	      fhdPhiRec->Fill(pRest[0].DeltaPhi(pTracksAll[j]), pTracksAll[j].Perp(), fXsection);	      
// 	    }
	  //and the Dalitz variables and Energy distributions in the rest frame
	  for (Int_t i = 0; i < nAccJets; i++)
	    xRec[i] = 2*vRestRec[i].E()/eSumRec;
	  
	  if(eSumRec <= 60)
	    fhX3X4Rec60->Fill(xRec[0], xRec[1], fXsection);
	  
	  if(eSumRec > 60 && eSumRec <= 100)
	    fhX3X4Rec60100->Fill(xRec[0], xRec[1], fXsection);
	  
	  if(eSumRec > 100)
	    fhX3X4Rec100->Fill(xRec[0], xRec[1], fXsection);
	  
	  if(nAccJets == 3 && nAccJets == nGenJets)
	    {
	      fhX3->Fill(xGen[0], TMath::Abs(xGen[0]-xRec[0])/xGen[0], fXsection);
	      fhX4->Fill(xGen[1], TMath::Abs(xGen[1]-xRec[1])/xGen[1], fXsection);
	      fhX5->Fill(xGen[2], TMath::Abs(xGen[2]-xRec[2])/xGen[2], fXsection);
	    }
	  
	  FillTopology(fhX3X4Rec, fhMu34Rec, fhMu45Rec, fhMu35Rec, xRec, pRestRec, fXsection);
    }
  if(fDebug)Printf("%s:%d",(char*)__FILE__,__LINE__);
  
  PostData(1, fList);

  if(fDebug)Printf("%s:%d Data Posted",(char*)__FILE__,__LINE__);
  
}

//__________________________________________________________________________________________________________________________________________________
void AliAnalysisTaskThreeJets::Terminate(Option_t *)
{
  if(fDebug)printf(" AliAnalysisTaskThreeJets::Terminate()");

} 

//_______________________________________User defined functions_____________________________________________________________________________________
void AliAnalysisTaskThreeJets::FillTopology(TH2F * Dalitz, TH1F * fhMu34, TH1F * fhMu45, TH1F * fhMu35, Double_t * x, TVector3 * pRest, Double_t xsection)
{
  //
  // fill the topology histos
  //
  Dalitz->Fill(x[0], x[1], xsection);
  fhMu35->Fill(TMath::Sqrt(x[0]*x[2]*(1-(pRest[0].Unit()).Dot(pRest[2].Unit()))/2), xsection);
  fhMu34->Fill(TMath::Sqrt(x[0]*x[1]*(1-(pRest[0].Unit()).Dot(pRest[1].Unit()))/2), xsection);
  fhMu45->Fill(TMath::Sqrt(x[1]*x[2]*(1-(pRest[1].Unit()).Dot(pRest[2].Unit()))/2), xsection);
}

//_____________________________________________________________________________________________________________________________

Bool_t AliAnalysisTaskThreeJets::IsPrimChar(TParticle* aParticle, Int_t aTotalPrimaries, Bool_t adebug)
{
  //
  // this function checks if a particle from the event generator (i.e. among the nPrim particles in the stack)
  // shall be counted as a primary particle
  //
  // This function or a equivalent should be available in some common place of AliRoot
  //
  // WARNING: Call this function only for particles that are among the particles from the event generator!
  // --> stack->Particle(id) with id < stack->GetNprimary()

  // if the particle has a daughter primary, we do not want to count it
  if (aParticle->GetFirstDaughter() != -1 && aParticle->GetFirstDaughter() < aTotalPrimaries)
  {
    if (adebug)
      printf("Dropping particle because it has a daughter among the primaries.\n");
    return kFALSE;
  }

  Int_t pdgCode = TMath::Abs(aParticle->GetPdgCode());
  

  // skip quarks and gluon
  if (pdgCode <= 10 || pdgCode == 21)
  {
    if (adebug)
      printf("Dropping particle because it is a quark or gluon.\n");
    return kFALSE;
  }

  Int_t status = aParticle->GetStatusCode();
  // skip non final state particles..
  if(status!=1){
    if (adebug)
      printf("Dropping particle because it is not a final state particle.\n");
    return kFALSE;
  }

  if (strcmp(aParticle->GetName(),"XXX") == 0)
  {
    Printf("WARNING: There is a particle named XXX (pdg code %d).", pdgCode);
    return kFALSE;
  }

  TParticlePDG* pdgPart = aParticle->GetPDG();

  if (strcmp(pdgPart->ParticleClass(),"Unknown") == 0)
  {
    Printf("WARNING: There is a particle with an unknown particle class (pdg code %d).", pdgCode);
    return kFALSE;
  }

  if (pdgPart->Charge() == 0)
  {
    if (adebug)
      printf("Dropping particle because it is not charged.\n");
    return kTRUE;
  }

  return kTRUE;
}

//______________________________________________________________________________________________________


//__________________________________________________________________________________________________________________________

Double_t AliAnalysisTaskThreeJets::Exponent(Double_t x,const Double_t * const par) const
{
  return par[0]*TMath::Power(1/TMath::E(), TMath::Power(par[1]/x, par[2])+0.5*TMath::Power((x-par[3])/par[0], 2))+par[4]*x;
}

Double_t AliAnalysisTaskThreeJets::Exponent2(Double_t x,const Double_t * const par) const
{
  return par[0]*TMath::Power(1/TMath::E(), TMath::Power(par[1]/x, par[2]))+par[3]*x;
}

Double_t AliAnalysisTaskThreeJets::Gauss(Double_t x,const Double_t * const par) const
{
  return 1/(par[1])*TMath::Power(1/TMath::E(), 0.5*(x-par[0])*(x-par[0])/(par[1]*par[1]));
}

Double_t AliAnalysisTaskThreeJets::Total(Double_t x,const Double_t * const par) const
{
  return Exponent(x, par)+Gauss(x, &par[4]);
}

  
