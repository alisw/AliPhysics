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

/* $Id$ */
/*
 Analysis Task 
 for Dijet Analysis
 based on AOD
*/

#include "AliAnalysisTaskDiJets.h"
#include "AliAODEvent.h"
#include "AliAODJet.h"
#include "AliAODDiJet.h"
#include <TClonesArray.h>
#include <TLorentzVector.h>
#include <TStyle.h>


ClassImp(AliAnalysisTaskDiJets)

////////////////////////////////////////////////////////////////////////

AliAnalysisTaskDiJets::AliAnalysisTaskDiJets():
    AliAnalysisTaskSE(),
    fDiJets(0),
    fDiJetsIn(0),
    fUseAODInput(kFALSE),
    fFillAOD(kFALSE),
    fJetBranch("jets"),
    fAOD(0),
    fHistList(0),
    fH1DeltaPt(0),
    fH1DeltaPhi(0),
    fH1PhiImbal(0),
    fH1Asym(0),
    fH2Pt2vsPt1(0),
    fH2DifvsSum(0)
{
  // Default constructor
}

//----------------------------------------------------------------------
AliAnalysisTaskDiJets::AliAnalysisTaskDiJets(const char* name):
    AliAnalysisTaskSE(name),
    fDiJets(0),
    fDiJetsIn(0),
    fUseAODInput(kFALSE),
    fFillAOD(kFALSE),
    fJetBranch("jets"),
    fAOD(0),
    fHistList(0),
    fH1DeltaPt(0),
    fH1DeltaPhi(0),
    fH1PhiImbal(0),
    fH1Asym(0),
    fH2Pt2vsPt1(0),
    fH2DifvsSum(0)
{
  // Default constructor
    DefineOutput(1, TList::Class());  
}

//----------------------------------------------------------------------
void AliAnalysisTaskDiJets::UserCreateOutputObjects()
{
// Create the output container
//
    if (fDebug) printf("AnalysisTaskDiJets::CreateOutPutData() \n");
    fDiJets = new TClonesArray("AliAODDiJet", 0);
    if (fFillAOD){
      fDiJets->SetName(Form("dijets_%s",fJetBranch.Data()));
      AddAODBranch("TClonesArray", &fDiJets);
	}

    if (!fHistList) fHistList = new TList();
    fHistList->SetOwner();
    Float_t pi=TMath::Pi();
    gStyle->SetPalette(1);

    fH1DeltaPt  = new TH1F("DeltaPt","Difference between the jets' Pt;#Deltap_{T} (GeV/c);Entries",150,0.,150.);
    fH1DeltaPt->SetMarkerSize(0.6);
    fH1DeltaPt->SetMarkerColor(4);
    fH1DeltaPt->SetMarkerStyle(21);
    fH1DeltaPt->SetOption("E");

    fH1DeltaPhi = new TH1F("DeltaPhi","Difference in the azimuthal angle;#Delta#phi;Entries",100,0.,pi);
    fH1DeltaPhi->SetMarkerSize(0.6);
    fH1DeltaPhi->SetMarkerColor(4);
    fH1DeltaPhi->SetMarkerStyle(21);
    fH1DeltaPhi->SetOption("E");

    fH1PhiImbal = new TH1F("PhiImb","Phi imbalance;#phi;Entries",100,-pi,pi);
    fH1PhiImbal->SetMarkerSize(0.6);
    fH1PhiImbal->SetMarkerColor(4);
    fH1PhiImbal->SetMarkerStyle(21);
    fH1PhiImbal->SetOption("E");

    fH1Asym     = new TH1F("Asym","Pt asymmetry;#Deltap_{T}/(p_{T,1}+p_{T,2});Entries",50,0.,1.);
    fH1Asym->SetMarkerSize(0.6);
    fH1Asym->SetMarkerColor(4);
    fH1Asym->SetMarkerStyle(21);
    fH1Asym->SetOption("E");

    fH2Pt2vsPt1 = new TH2F("Pt2vsPt1","Pt2 versus Pt1;p_{T,1} (GeV/c);p_{T,2} (GeV/c)",250,0.,250.,250,0.,250.);
    fH2Pt2vsPt1->SetOption("cont0");

    fH2DifvsSum = new TH2F("DifvsSum","Pt difference vs Pt sum;p_{T,1}+p_{T,2} (GeV/c);#Deltap_{T} (GeV/c)",400,0.,400.,150,0.,150.);
    fH2DifvsSum->SetOption("cont0");

    fHistList->Add(fH1DeltaPt);
    fHistList->Add(fH1DeltaPhi);
    fHistList->Add(fH1PhiImbal);
    fHistList->Add(fH1Asym);
    fHistList->Add(fH2Pt2vsPt1);
    fHistList->Add(fH2DifvsSum);
}

//----------------------------------------------------------------------
void AliAnalysisTaskDiJets::Init()
{
    // Initialization
    if (fDebug) printf("AnalysisTaskDiJets::Init() \n");
}

//----------------------------------------------------------------------
void AliAnalysisTaskDiJets::UserExec(Option_t */*option*/)
{
// Execute analysis for current event
//
    if (fDiJets) fDiJets->Delete();

    if(fUseAODInput){
      fAOD = dynamic_cast<AliAODEvent*> (InputEvent());
      if(!fAOD){
        // We do not have an input AOD, look in the output
        if (fDebug) printf("%s:%d No AOD event in the input\n",(char*)__FILE__,__LINE__);
        return;
      }
    } else {
      fAOD = AODEvent();
      if(!fAOD){
        if (fDebug) printf("%s:%d AODEvent not found in the Output",(char*)__FILE__,__LINE__);
        return;
      }
    }

    TClonesArray* jets = (TClonesArray*) fAOD->FindListObject(fJetBranch.Data());
    // N.B. if we take the aod from the output this is always
    // empty and since it is the same as fDiJets 
    fDiJetsIn = (TClonesArray*) (fAOD->GetList()->FindObject("dijets"));

    if (fDiJetsIn) {
      if (fDebug) printf("Found %d dijets in old list \n", fDiJetsIn->GetEntries());
      AliAODJet* jj1, *jj2;
      AliAODDiJet* testJ;

      if (fDiJetsIn->GetEntries() > 0) {
        testJ = (AliAODDiJet*) (fDiJetsIn->At(0));
        jj1 = testJ->Jet(0);
        jj1->Print("");
        jj2 = testJ->Jet(1);
        jj2->Print("");
      }
    }

    Int_t nj = jets->GetEntriesFast();
    if (fDebug) printf("There are %5d jets in the event \n", nj);

    if (nj < 2){
      PostData(1, fHistList);
      return;
    }
    AliAODJet* jet1 = (AliAODJet*) (jets->At(0));
    TLorentzVector v1 = *(jet1->MomentumVector());
    AliAODJet* jet2 = (AliAODJet*) (jets->At(1));
    TLorentzVector v2 = *(jet2->MomentumVector());
    TLorentzVector v = v1 + v2;
    if (fDiJets) {
	Int_t ndi = fDiJets->GetEntriesFast();
	TClonesArray &lref = *fDiJets;
	new(lref[ndi]) AliAODDiJet(v);
	AliAODDiJet* dijet = (AliAODDiJet*) (fDiJets->At(ndi));
	dijet->SetJetRefs(jet1, jet2);
	fH1DeltaPhi->Fill(dijet->DeltaPhi());
	fH1PhiImbal->Fill(dijet->PhiImbalance());

    }
    
    fH1DeltaPt->Fill(jet1->Pt()-jet2->Pt());
    fH1Asym->Fill((jet1->Pt()-jet2->Pt())/(jet1->Pt()+jet2->Pt()));
    fH2Pt2vsPt1->Fill(jet1->Pt(),jet2->Pt());
    fH2DifvsSum->Fill(jet1->Pt()+jet2->Pt(),jet1->Pt()-jet2->Pt());

    PostData(1, fHistList);
    return;
}

//----------------------------------------------------------------------
void AliAnalysisTaskDiJets::Terminate(Option_t */*option*/)
{
// Terminate analysis
//
    if (fDebug) printf("AnalysisDiJets: Terminate() \n");
}

