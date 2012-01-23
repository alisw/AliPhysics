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
//
// Task for Efficiency studies
// Used for testing classes AliHFEcontainer and AliHFEfilter
// Creates Efficiency Histograms
//
// Author:
//   Markus Fasel <M.Fasel@gsi.de>
//
#include <TAxis.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TH1F.h>
#include <TH2.h>
#include <TIterator.h>
#include <TLegend.h>
#include <TObjArray.h>
#include <TPad.h>

#include "AliAnalysisManager.h"
#include "AliCFAcceptanceCuts.h"
#include "AliCFTrackIsPrimaryCuts.h"
#include "AliCFContainer.h"
#include "AliCFEffGrid.h"
#include "AliESDEvent.h"
#include "AliHFEcollection.h"
#include "AliHFEcontainer.h"
#include "AliHFEcutStep.h"
#include "AliHFEefficiency.h"
#include "AliHFEextraCuts.h"
#include "AliHFEtrackFilter.h"
#include "AliHFEtools.h"
#include "AliMCEvent.h"
#include "AliMCEventHandler.h"

ClassImp(AliHFEefficiency)

AliHFEefficiency::AliHFEefficiency():
  AliAnalysisTaskSE(),
  fFilter(NULL),
  fMCcut(NULL),
  fAcceptanceCuts(NULL),
  fEfficiency(NULL),
  fOutput(NULL),
  fCutTRD(kFALSE)
{
  //
  // Default constructor
  //
}

AliHFEefficiency::AliHFEefficiency(const Char_t *name):
  AliAnalysisTaskSE(name),
  fFilter(NULL),
  fMCcut(NULL),
  fAcceptanceCuts(NULL),
  fEfficiency(NULL),
  fOutput(NULL),
  fCutTRD(kFALSE)
{
  //
  // Main Constructor
  //
  DefineOutput(1, AliHFEcontainer::Class());
  DefineOutput(2, TList::Class());

  SetRunTerminate();
}

AliHFEefficiency::~AliHFEefficiency(){
  //
  // Destructor
  //
  if(fFilter) delete fFilter;
  if(fEfficiency) delete fEfficiency;
  if(fAcceptanceCuts) delete fAcceptanceCuts;
  if(fOutput) delete fOutput;
}

void AliHFEefficiency::UserCreateOutputObjects(){
  //
  // Create the output objects
  //
  fEfficiency = new AliHFEcontainer("Efficiency");
  fEfficiency->SetNumberOfVariables(4);
  // InitializeVariables
  fEfficiency->SetVariableName(0, "pt");
  fEfficiency->MakeLogarithmicBinning(0, 40, 0.1, 10);
  fEfficiency->SetVariableName(1, "Eta");
  fEfficiency->MakeLinearBinning(1, 20, -0.9, 0.9);
  fEfficiency->SetVariableName(2, "phi");
  fEfficiency->MakeLinearBinning(2 , 18, 0, 2 * TMath::Pi());
  fEfficiency->SetVariableName(3, "Charge");
  fEfficiency->MakeLinearBinning(3, 2, -1.1, 1.1);
  // New container
  fEfficiency->CreateContainer("MCFilter", "Efficiency for MC Studies", 2);
  AliCFContainer *cont = fEfficiency->GetCFContainer("MCFilter");
  cont->SetStepTitle(0, "MC Signal");
  cont->SetStepTitle(1, "MC Signal in acceptance");

  fFilter = new AliHFEtrackFilter("testfilter");
  fFilter->GenerateCutSteps();
  fFilter->InitCF(fEfficiency);
  fMCcut = fFilter->GetMCSignalCuts();
  if(fCutTRD){
    AliHFEcutStep *hfeTRD = fFilter->GetCutStep("HFETRD");
    AliHFEextraCuts *hfecuts = dynamic_cast<AliHFEextraCuts *>(hfeTRD->GetCut("HFETRDCuts"));
    if(hfecuts) hfecuts->SetMinTrackletsTRD(4);
  }
  AliHFEcutStep *hfeITS = fFilter->GetCutStep("HFEITS");
  AliHFEextraCuts *hfeitscuts = dynamic_cast<AliHFEextraCuts *>(hfeITS->GetCut("HFEPixelsCuts"));
  if(hfeitscuts)hfeitscuts->SetRequireITSpixel(AliHFEextraCuts::kFirst);

  AliHFEcutStep *primary = fFilter->GetCutStep("Primary");
  AliCFTrackIsPrimaryCuts *primcuts = dynamic_cast<AliCFTrackIsPrimaryCuts *>(primary->GetCut("PrimaryCuts"));
  if(primcuts){ 
    primcuts->SetMaxDCAToVertexXY(3);
    primcuts->SetMaxDCAToVertexZ(5);
  }

  fAcceptanceCuts = new AliCFAcceptanceCuts("Acceptance", "MC Acceptance Cuts");
  fAcceptanceCuts->SetMinNHitITS(3);
  fAcceptanceCuts->SetMinNHitTPC(2);
  fAcceptanceCuts->SetMinNHitTRD(0);

  // create additional histos for testing purpose
  fOutput = new AliHFEcollection("histos", "QA histograms");
  fOutput->CreateTH1F("hNtracks","Number of Tracks per Event", 100, 0, 100);
  fOutput->CreateTH1F("hPt","Pt of the tracks", 40, 0.1, 10, 0);
  fOutput->CreateTH1F("kinkIndex", "Kink Index", 400, -200, 200);
  fOutput->CreateTH1F("itspixel", "ITS PIXEL", 2, 0, 2);
  fOutput->CreateTH1F("mcmother", "Mother PDG", 1000, -500, 500);
  fOutput->CreateTH2F("ptres", "Pt Resolution", 40, 0.1, 10, 200, -1.5, 1.5, 0);
}

void AliHFEefficiency::UserExec(Option_t *){
  //
  // Event Loop
  // Filter track, fill Efficiency container
  //
  AliESDEvent *esdevent = dynamic_cast<AliESDEvent *>(fInputEvent);
  if(!esdevent){
    AliError("ESD Event required");
    return;
  }
  fEfficiency->NewEvent();
  fFilter->SetRecEvent(fInputEvent);
  if(fMCEvent){
    AliMCEventHandler *mcH = dynamic_cast<AliMCEventHandler *>(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
   if ( ! mcH ) {
     AliError("Cannot get MC truth event handler");
     return;
    }  
    if(mcH &&(!mcH->InitOk())) return;
    if(mcH &&(!mcH->TreeK())) return;
    if(mcH &&(!mcH->TreeTR())) return;
    fFilter->SetMC(fMCEvent);
    FilterMC();
  }
  fFilter->FilterTracks(esdevent);
  TObjArray *tracks = fFilter->GetFilteredTracks();
  TIterator *iter = tracks->MakeIterator();
  fOutput->Fill("hNtracks", tracks->GetEntriesFast());
  AliESDtrack *track = NULL;
  while((track = dynamic_cast<AliESDtrack *>(iter->Next()))){
    fOutput->Fill("hPt", track->Pt());  
    fOutput->Fill("kinkIndex", track->GetKinkIndex(0));
    if(TESTBIT(track->GetITSClusterMap(), 0))
      fOutput->Fill("itspixel",1);
    else
      fOutput->Fill("itspixel",0);
    if(fMCEvent){
      AliMCParticle *mctrack = dynamic_cast<AliMCParticle *>(fMCEvent->GetTrack(TMath::Abs(track->GetLabel())));
      if(mctrack){
        Int_t motherLabel = mctrack->Particle()->GetFirstMother();
        if(motherLabel){
          AliMCParticle *mother = dynamic_cast<AliMCParticle *>(fMCEvent->GetTrack(motherLabel));
          if(mother)fOutput->Fill("mcmother", mother->Particle()->GetPdgCode());
        }
        fOutput->Fill("ptres", mctrack->Pt(), (track->Pt() - mctrack->Pt())/mctrack->Pt());
      }
    }
  }
  delete iter;
  fFilter->Flush();
  PostData(1, fEfficiency);
  PostData(2, fOutput->GetList());
}

void AliHFEefficiency::Terminate(Option_t *){
  //
  // Evaluate Results
  //
  fEfficiency = dynamic_cast<AliHFEcontainer *>(GetOutputData(1));
  if(!fEfficiency){
    AliError("No Output data available");
    return;
  }

  if(!IsRunningTerminate()) return;
  PostProcess();

  TList *l = dynamic_cast<TList *>(GetOutputData(2));
  if(l) DrawPtResolution(l);
}

void AliHFEefficiency::FilterMC(){
  //
  // Cut MC tracks
  //
  Double_t cont[4] = {0., 0., 0., 0.};
  AliMCParticle *track = NULL;
  for(Int_t itrack = 0; itrack < fMCEvent->GetNumberOfTracks(); itrack++){
    track = dynamic_cast<AliMCParticle *>(fMCEvent->GetTrack(itrack));
    if(!track) continue;
    if(!fMCcut->IsSelected(track)) continue;
    cont[0] = track->Pt();
    cont[1] = track->Eta();
    cont[2] = track->Phi();
    cont[3] = track->Charge()/3;
    //AliMCParticle *mother = dynamic_cast<AliMCParticle *>(fMCEvent->GetTrack(track->Particle()->GetFirstMother()));
    //if(TMath::Abs(mother->Particle()->GetPdgCode()) != 443) continue;
    fEfficiency->FillCFContainer("MCFilter", 0, cont);
    if(!fAcceptanceCuts->IsSelected(track)) continue;
    fEfficiency->FillCFContainer("MCFilter", 1, cont);
  }
}

void AliHFEefficiency::Load(const char* filename){
  //
  // Load results for post processing
  //
  TFile *input = TFile::Open(filename);
  AliHFEcontainer *ctmp = dynamic_cast<AliHFEcontainer *>(input->Get("Efficiency"));
  if(ctmp) fEfficiency = dynamic_cast<AliHFEcontainer *>(ctmp->Clone());
  input->Close();
  delete input;
}

void AliHFEefficiency::PostProcess(){
  //
  // Calculate Efficiencies
  //
  fEfficiency->Print();

  // Prepare containers
  AliCFContainer *crec  = fEfficiency->MakeMergedCFContainer("merged", "Merged Container", "MCFilter:testfilter_container_signal");
  AliCFContainer *cmc   = fEfficiency->MakeMergedCFContainer("merged", "Merged Container", "MCFilter:testfilter_container_signalMC");
  AliCFEffGrid *gridRec = new AliCFEffGrid("effrec", "Efficiency Calculation using rec vars", *crec);
  AliCFEffGrid *gridMC  = new AliCFEffGrid("effmc", "Efficiency Calculation using mc vars", *cmc);

  TCanvas *ctmp = NULL;
  for(Int_t ivar = 0; ivar < crec->GetNVar(); ivar++){
    ctmp = new TCanvas(Form("cEffSignal%s", crec->GetVarTitle(ivar)), Form("Signal Efficiency (%s)", crec->GetVarTitle(ivar)));
    ctmp->cd();
    DrawSignalEfficiency(gridRec, crec, ivar);
    
    ctmp = new TCanvas(Form("cCutEff%s", crec->GetVarTitle(ivar)), Form("Cut Step Efficiency (%s)", crec->GetVarTitle(ivar)));
    ctmp->cd();
    DrawCutEfficiency(gridRec, crec, ivar);
  
    ctmp = new TCanvas(Form("cEffSignalMC%s", crec->GetVarTitle(ivar)), Form("Signal Efficiency (%s, MC)", crec->GetVarTitle(ivar)));
    ctmp->cd();
    DrawSignalEfficiency(gridMC, cmc, ivar);

    ctmp = new TCanvas(Form("cCutEffMC%s", crec->GetVarTitle(ivar)), Form("Cut Step Efficiency (%s, MC)", crec->GetVarTitle(ivar)));
    ctmp->cd();
    DrawCutEfficiency(gridMC, cmc, ivar);
  } 

  CalculatePTsmearing();

  delete gridRec;
  delete gridMC;
  delete crec;
  delete cmc;
}

void AliHFEefficiency::DrawSignalEfficiency(AliCFEffGrid *eff, AliCFContainer *cont, Int_t var){
  //
  // Efficiency of a cut step with respect to the signal
  //
  TH1 *hEff = NULL;
  Bool_t kFirst = kTRUE;
  TLegend *leg = new TLegend(0.5, 0.7, 0.89, 0.89);
  for(Int_t istep = 1; istep < cont->GetNStep(); istep++){
    eff->CalculateEfficiency(istep, 0);
    hEff = eff->Project(var);
    hEff->SetTitle(Form("Signal Efficiency (%s)", cont->GetVarTitle(var)));
    hEff->GetYaxis()->SetRangeUser(0., 1.6);
    hEff->GetYaxis()->SetTitle("Efficiency");
    hEff->SetStats(kFALSE);
    hEff->SetMarkerStyle(20+istep);
    hEff->SetMarkerColor(kBlue - 5 + istep);
    hEff->Draw(kFirst ? "ep" : "epsame");
    kFirst = kFALSE;
    leg->AddEntry(hEff, cont->GetStepTitle(istep), "p");
  }
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->Draw();
  gPad->Update();
}

void AliHFEefficiency::DrawCutEfficiency(AliCFEffGrid *eff, AliCFContainer *cont, Int_t var){
  //
  // Efficiency of a cut step with respect to the previous one
  //
  Bool_t kFirst = kTRUE;
  TH1 *hEff = NULL;
  TLegend *leg = new TLegend(0.5, 0.7, 0.89, 0.89);
  TString effTitle;
  Int_t start = 0, length = 0;
  for(Int_t istep = 1; istep < cont->GetNStep(); istep++){
    eff->CalculateEfficiency(istep, istep -1);
    hEff = eff->Project(var);
    effTitle = hEff->GetTitle();
    start = effTitle.Index("projected"); length = effTitle.Length() - start;
    effTitle.Remove(start - 1, length + 1);
    effTitle.ReplaceAll("Efficiency: ", "");
    hEff->SetTitle(Form("Cut Efficiency (%s)", cont->GetVarTitle(var)));
    hEff->GetYaxis()->SetRangeUser(0., 1.6);
    hEff->GetYaxis()->SetTitle("Efficiency");
    hEff->SetStats(kFALSE);
    hEff->SetMarkerStyle(20+istep);
    hEff->SetMarkerColor(kBlue - 5 + istep);
    hEff->Draw(kFirst ? "ep" : "epsame");
    kFirst = kFALSE;
    leg->AddEntry(hEff, effTitle.Data(), "p");
  }
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->Draw();
  gPad->Update();
}

void AliHFEefficiency::CalculatePTsmearing(){
  //
  // Efficiency of a cut step with respect to the same cut step filled with MC 
  // Information
  //
  AliCFContainer *cont = fEfficiency->MakeMergedCFContainer("merged", "Merged Container", "testfilter_container_signal:testfilter_container_signalMC");
  AliCFEffGrid *grid = new AliCFEffGrid("effc", "Efficiency Calculation", *cont);
  Bool_t kFirst = kTRUE;
  TH1 *hEff = NULL;
  TCanvas *c2 = new TCanvas("cSmear", "pT modification");
  c2->cd();
  TLegend *leg1 = new TLegend(0.5, 0.7, 0.89, 0.89);
  Int_t nStep = cont->GetNStep()/2;
  for(Int_t istep = 0; istep < nStep; istep++){
    grid->CalculateEfficiency(istep, istep + nStep);
    hEff = grid->Project(0);
    hEff->SetStats(kFALSE);
    hEff->SetMarkerStyle(20+istep);
    hEff->SetMarkerColor(kBlue - 5 + istep);
    hEff->Draw(kFirst ? "ep" : "epsame");
    kFirst = kFALSE;
    leg1->AddEntry(hEff, cont->GetStepTitle(istep), "p");
  }
  leg1->Draw();
  gPad->Update();

  delete cont;
  delete grid;
}

void AliHFEefficiency::DrawPtResolution(const TList * const l){
  //
  // Draw pt resolution
  //
  TH2 *h2 = dynamic_cast<TH2 *>(l->FindObject("ptres"));
  if(!h2) return;
  TGraphErrors *mean = new TGraphErrors(h2->GetNbinsX());
  TGraphErrors *rms = new TGraphErrors(h2->GetNbinsX());

  Double_t pt = 0., pterr = 0., mymean = 0., myrms = 0., mymeanerr = 0., myrmserr = 0.;
  TH1 *htmp = NULL;
  for(Int_t imom = 1; imom <= h2->GetXaxis()->GetLast(); imom++){
    pt = h2->GetXaxis()->GetBinCenter(imom);
    pterr = h2->GetXaxis()->GetBinWidth(imom)/2.;
    htmp = h2->ProjectionY("py", imom, imom);
    mymean = htmp->GetMean();
    myrms = htmp->GetRMS();
    mymeanerr = htmp->GetMeanError();
    myrmserr = htmp->GetRMSError();
    delete htmp;
    mean->SetPoint(imom -1, pt, mymean);
    rms->SetPoint(imom -1, pt, myrms);
    mean->SetPointError(imom-1,pterr,mymeanerr);
    rms->SetPointError(imom-1,pterr,myrmserr);
  }

  TCanvas *c1 = new TCanvas("cPtShift", "pT Shift");
  c1->cd();
  mean->SetMarkerStyle(22);
  mean->SetMarkerColor(kBlue);
  mean->Draw("ape");
  rms->SetMarkerStyle(22);
  rms->SetMarkerColor(kBlack);
  rms->Draw("pesame");
  TLegend *leg = new TLegend(0.5, 0.7, 0.89, 0.89);
  leg->AddEntry(mean, "Mean", "lp");
  leg->AddEntry(rms, "RMS", "lp");
  leg->Draw();
  gPad->Update();
} 

