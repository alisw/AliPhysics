#include <TFile.h>
#include <TH1.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TAxis.h>
#include <TStyle.h>
#include <TChain.h>
#include <TTree.h>
#include <TMath.h>
#include <TClonesArray.h>
#include <TLegend.h>

#include <AliXRDPROOFtoolkit.h>

#include "DebugClasses.C"
#include "my_tools.C"

#include <iostream>

using namespace std;

/*
  To run code:
  ============

  Info:
  * I did not recheck this code. For now I would just use the efficiency values I have.
  * The code could be made nicer. Esepcially some plots could be put in folders.

  Use AliRoot because of AliXRDPROOFtoolkit:
  gSystem->AddIncludePath("-I$ALICE_ROOT/TPC")
  gSystem->AddIncludePath("-I../lib")
  gSystem->AddIncludePath("-I../grid")
  gSystem->AddIncludePath("-I../macros")
  gROOT->SetMacroPath(".:../macros:../grid:../lib/")
  .L my_tools.C+
  .L DebugClasses.C+
  .L createEfficiency.C+

  Examples of visualization:
  DrawEfficiency("lhc10d_eff_pythia.root", "eff.root")

  // This is the correction I did for the low pT guys
  DrawCorrection("lhc10d_eff_pythia.root", "lhc10d_eff_phojet.root")



  CreateEff("lhc10d_mc_all.dat", 0, "lhc10d_eff_all.root")

  
*/

void DrawEfficiency(const Char_t* fileName, const Char_t* effFileName);
void DrawCorrection(const Char_t* fileNamePYTHIA, const Char_t* fileNamePHOJET);
TH1D* HistInvert(TH1D* hist);


void CreateEff(const Char_t* mcFileName, Int_t maxEvents, const Char_t* mcOutFileName,
	       Float_t centLow=-20, Float_t centHigh=-5)
{  
  gStyle->SetOptStat(0);
  //
  // Create output
  //
  TFile* outFile = new TFile(mcOutFileName, "RECREATE");

  const Int_t nPid = 7;
  TH1D* hMcIn[nPid]     = {0, 0, 0, 0, 0, 0, 0 };
  TH1D* hMcOut[nPid]    = {0, 0, 0, 0, 0, 0, 0 };
  TH1D* hMcSec[nPid]    = {0, 0, 0, 0, 0, 0, 0 };
  TH1D* hMcEff[nPid]     = {0, 0, 0, 0, 0, 0, 0 };
  TH1D* hMcInNeg[nPid]  = {0, 0, 0, 0, 0, 0, 0 };
  TH1D* hMcOutNeg[nPid] = {0, 0, 0, 0, 0, 0, 0 };
  TH1D* hMcSecNeg[nPid] = {0, 0, 0, 0, 0, 0, 0 };
  TH1D* hMcEffNeg[nPid]     = {0, 0, 0, 0, 0, 0, 0 };
  TH1D* hMcInPos[nPid]  = {0, 0, 0, 0, 0, 0, 0 };
  TH1D* hMcOutPos[nPid] = {0, 0, 0, 0, 0, 0, 0 };
  TH1D* hMcSecPos[nPid] = {0, 0, 0, 0, 0, 0, 0 };
  TH1D* hMcEffPos[nPid]     = {0, 0, 0, 0, 0, 0, 0 };

  Int_t color[nPid] = {1, 2, 3, 4, 5, 1, 1};

  const Int_t nPtBins = 68;
  Double_t xBins[nPtBins+1] = {0. ,  0.05, 0.1,  0.15, 0.2,  0.25, 0.3,  0.35, 0.4,  0.45,
			       0.5,  0.55, 0.6,  0.65, 0.7,  0.75, 0.8,  0.85, 0.9,  0.95,
			       1.0,  1.1 , 1.2,  1.3 , 1.4,  1.5 , 1.6,  1.7 , 1.8,  1.9 ,
			       2.0,  2.2 , 2.4,  2.6 , 2.8,  3.0 , 3.2,  3.4 , 3.6,  3.8 ,
			       4.0,  4.5 , 5.0,  5.5 , 6.0,  6.5 , 7.0,  8.0 , 9.0,  10.0,
			       11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 18.0, 20.0, 22.0, 24.0,
			       26.0, 28.0, 30.0, 32.0, 34.0, 36.0, 40.0, 45.0, 50.0 };
  
  for(Int_t pid = 0; pid < nPid; pid++) {
    
    hMcIn[pid] = new TH1D(Form("hIn%d", pid), Form("Efficiency (pid %d)", pid), 
			  nPtBins, xBins); 
    hMcInNeg[pid] = new TH1D(Form("hInNeg%d", pid), Form("Efficiency (pid %d, q < 0)", pid), 
			     nPtBins, xBins); 
    hMcInPos[pid] = new TH1D(Form("hInPos%d", pid), Form("Efficiency (pid %d, q < 0)", pid), 
			     nPtBins, xBins); 

    hMcIn[pid]->Sumw2();
    hMcIn[pid]->SetMarkerStyle(29);
    hMcIn[pid]->SetMarkerColor(color[pid]);
    hMcInNeg[pid]->Sumw2();
    hMcInNeg[pid]->SetMarkerStyle(24);
    hMcInNeg[pid]->SetMarkerColor(color[pid]);
    hMcInPos[pid]->Sumw2();
    hMcInPos[pid]->SetMarkerStyle(20);
    hMcInPos[pid]->SetMarkerColor(color[pid]);

    hMcOut[pid] = new TH1D(Form("hMcOut%d", pid), Form("MC out (pid %d)", pid), 
			   nPtBins, xBins); 
    hMcOutNeg[pid] = new TH1D(Form("hMcOutNeg%d", pid), Form("MC out (pid %d, q < 0)", pid), 
			      nPtBins, xBins); 
    hMcOutPos[pid] = new TH1D(Form("hMcOutPos%d", pid), Form("MC out (pid %d, q < 0)", pid), 
			      nPtBins, xBins); 

    hMcOut[pid]->Sumw2();
    hMcOut[pid]->SetMarkerStyle(29);
    hMcOut[pid]->SetMarkerColor(color[pid]);
    hMcOutNeg[pid]->Sumw2();
    hMcOutNeg[pid]->SetMarkerStyle(24);
    hMcOutNeg[pid]->SetMarkerColor(color[pid]);
    hMcOutPos[pid]->Sumw2();
    hMcOutPos[pid]->SetMarkerStyle(20);
    hMcOutPos[pid]->SetMarkerColor(color[pid]);

    hMcSec[pid] = new TH1D(Form("hSec%d", pid), Form("Secondaries (pid %d)", pid), 
			  nPtBins, xBins); 
    hMcSecNeg[pid] = new TH1D(Form("hSecNeg%d", pid), Form("Secondaries (pid %d, q < 0)", pid), 
			     nPtBins, xBins); 
    hMcSecPos[pid] = new TH1D(Form("hSecPos%d", pid), Form("Secondaries (pid %d, q < 0)", pid), 
			     nPtBins, xBins); 

    hMcSec[pid]->Sumw2();
    hMcSec[pid]->SetMarkerStyle(29);
    hMcSec[pid]->SetMarkerColor(color[pid]);
    hMcSecNeg[pid]->Sumw2();
    hMcSecNeg[pid]->SetMarkerStyle(24);
    hMcSecNeg[pid]->SetMarkerColor(color[pid]);
    hMcSecPos[pid]->Sumw2();
    hMcSecPos[pid]->SetMarkerStyle(20);
    hMcSecPos[pid]->SetMarkerColor(color[pid]);

    hMcEff[pid] = new TH1D(Form("hEff%d", pid), Form("Efficiency (pid %d)", pid), 
			  nPtBins, xBins); 
    hMcEffNeg[pid] = new TH1D(Form("hEffNeg%d", pid), Form("Efficiency (pid %d, q < 0)", pid), 
			     nPtBins, xBins); 
    hMcEffPos[pid] = new TH1D(Form("hEffPos%d", pid), Form("Efficiency (pid %d, q < 0)", pid), 
			     nPtBins, xBins); 

    hMcEff[pid]->Sumw2();
    hMcEff[pid]->SetMarkerStyle(29);
    hMcEff[pid]->SetMarkerColor(color[pid]);
    hMcEffNeg[pid]->Sumw2();
    hMcEffNeg[pid]->SetMarkerStyle(24);
    hMcEffNeg[pid]->SetMarkerColor(color[pid]);
    hMcEffPos[pid]->Sumw2();
    hMcEffPos[pid]->SetMarkerStyle(20);
    hMcEffPos[pid]->SetMarkerColor(color[pid]);
  }
      
  TTree* Tree   = 0;
  if(strstr(mcFileName, ".dat")) {
    
    AliXRDPROOFtoolkit tool;
    TChain* chain = tool.MakeChain(mcFileName,"tree", 0, 1000);
    chain->Lookup();
    Tree = chain;
  } else {
    TFile* mcFile = FindFileFresh(mcFileName);
    if(!mcFile)
      return;
    
    Tree = (TTree*)mcFile->Get("tree");
  }
  if(!Tree)
    return;
  
  
  DeDxEvent* event = 0;
  TClonesArray* trackArray = 0;
  TClonesArray* mcTrackArray = 0;
  Tree->SetBranchAddress("event", &event);
  Tree->SetBranchAddress("track"  , &trackArray);
  Tree->SetBranchAddress("trackMC"  , &mcTrackArray);
  Int_t nEvents = Tree->GetEntries();
  cout << "Number of events: " << nEvents << endl;
  
  if(maxEvents>0 && maxEvents < nEvents) {
    
    nEvents = maxEvents;
    cout << "N events was reduced to: " << maxEvents << endl;
  }
  
  Int_t currentRun = 0;

  for(Int_t n = 0; n < nEvents; n++) {
    
    Tree->GetEntry(n);
    
    if((n+1)%1000000==0)
      cout << "Event: " << n+1 << "/" << nEvents << endl;
    
    if(event->run != currentRun) {
      
      cout << "New run: " << event->run << endl;
      currentRun = event->run;
    }

    if(event->cent < centLow || event->cent > centHigh)
      continue;

    const Int_t nMcTracks = mcTrackArray->GetEntries();
      
    for(Int_t i = 0; i < nMcTracks; i++) {
	
      DeDxTrackMC* trackMC = (DeDxTrackMC*)mcTrackArray->At(i);

      // if(TMath::Abs(trackMC->pdgMC)==3312 || TMath::Abs(trackMC->pdgMC)==3334) 
      // 	continue; // Xi-!
	
      hMcIn[0]->Fill(trackMC->ptMC);
      hMcIn[trackMC->pidMC]->Fill(trackMC->ptMC);
	
      if(trackMC->qMC < 0) {
	  
	hMcInNeg[0]->Fill(trackMC->ptMC);
	hMcInNeg[trackMC->pidMC]->Fill(trackMC->ptMC);
      } else {
	  
	hMcInPos[0]->Fill(trackMC->ptMC);
	hMcInPos[trackMC->pidMC]->Fill(trackMC->ptMC);
      }
    }

    const Int_t nTracks = trackArray->GetEntries();
      
    for(Int_t i = 0; i < nTracks; i++) {
	
      DeDxTrack* track = (DeDxTrack*)trackArray->At(i);
	
      if(!(track->filter&1))
	continue;
      
      // if(TMath::Abs(track->mother)==3312 || TMath::Abs(track->mother)==3334) 
      // 	continue; // Xi+- or Omega+-!

      hMcOut[0]->Fill(track->pt);
      hMcOut[track->pid]->Fill(track->pt);
	
      if(track->q < 0) {
	  
	hMcOutNeg[0]->Fill(track->pt);
	hMcOutNeg[track->pid]->Fill(track->pt);
      } else {
	  
	hMcOutPos[0]->Fill(track->pt);
	hMcOutPos[track->pid]->Fill(track->pt);
      }
	
      if(track->primary==0) {
	hMcSec[0]->Fill(track->pt);
	hMcSec[track->pid]->Fill(track->pt);
	  
	if(track->q < 0) {
	    
	  hMcSecNeg[0]->Fill(track->pt);
	  hMcSecNeg[track->pid]->Fill(track->pt);
	} else {
	    
	  hMcSecPos[0]->Fill(track->pt);
	  hMcSecPos[track->pid]->Fill(track->pt);
	}
      }
    }
  }

  TH1D* hMcInPiKP = (TH1D*)hMcIn[1]->Clone("hMcInPiKP");
  hMcInPiKP->Add(hMcIn[2]);
  hMcInPiKP->Add(hMcIn[3]);
  hMcInPiKP->Divide(hMcIn[0]);

  for(Int_t pid = 0; pid < nPid; pid++) {
    
    hMcSec[pid]->Divide(hMcSec[pid], hMcOut[pid]); 
    hMcSecNeg[pid]->Divide(hMcSecNeg[pid], hMcOutNeg[pid]); 
    hMcSecPos[pid]->Divide(hMcSecPos[pid], hMcOutPos[pid]); 

    hMcEff[pid]   ->Divide(hMcOut[pid], hMcIn[pid]); 
    hMcEffNeg[pid]->Divide(hMcOutNeg[pid], hMcInNeg[pid]); 
    hMcEffPos[pid]->Divide(hMcOutPos[pid], hMcInPos[pid]); 
  }

  outFile->Write();
  outFile->Close();  
}

//_______________________________________________________________________
void DrawEfficiency(const Char_t* fileName, const Char_t* effFileName)
{
  TFile* file = FindFileFresh(fileName);
  if(!file)
    return;
  
  gStyle->SetOptStat(0);

  const Double_t ptmin  =  3.01;
  const Double_t ptmax  = 19.999;
  //  const Double_t ptmax  = 49.999;
  
  const Double_t effmin = 0.6;
  const Double_t effmax = 0.9;

  // const Double_t effmin = 0.0;
  // const Double_t effmax = 1.5;

  const Double_t secmin = 0.0;
  const Double_t secmax = 0.1;

  const Int_t nPid = 7;

  TH1D* hEff[nPid]     = {0, 0, 0, 0, 0, 0, 0 };
  TH1D* hEffNeg[nPid]  = {0, 0, 0, 0, 0, 0, 0 };
  TH1D* hEffPos[nPid]  = {0, 0, 0, 0, 0, 0, 0 };
  TH1D* hSec[nPid]    = {0, 0, 0, 0, 0, 0, 0 };
  TH1D* hSecNeg[nPid] = {0, 0, 0, 0, 0, 0, 0 };
  TH1D* hSecPos[nPid] = {0, 0, 0, 0, 0, 0, 0 };
  TH1D* hOut[nPid]    = {0, 0, 0, 0, 0, 0, 0 };

  for(Int_t pid = 0; pid < nPid; pid++) {

    hEff[pid] = (TH1D*)file->Get(Form("hEff%d", pid));

    hEffNeg[pid] = (TH1D*)file->Get(Form("hEffNeg%d", pid));

    hEffPos[pid] = (TH1D*)file->Get(Form("hEffPos%d", pid));

    hSec[pid] = (TH1D*)file->Get(Form("hSec%d", pid));

    hSecNeg[pid] = (TH1D*)file->Get(Form("hSecNeg%d", pid));

    hSecPos[pid] = (TH1D*)file->Get(Form("hSecPos%d", pid));

    hOut[pid] = (TH1D*)file->Get(Form("hMcOut%d", pid));

    // hEff[pid]->Rebin(4);
    // hEffNeg[pid]->Rebin(4);
    // hEffPos[pid]->Rebin(4);

    // hSec[pid]->Rebin(4);
    // hSecNeg[pid]->Rebin(4);
    // hSecPos[pid]->Rebin(4);

    // hOut[pid]->Rebin(4);

    // hEff[pid]->Scale(0.25);
    // hEffNeg[pid]->Scale(0.25);
    // hEffPos[pid]->Scale(0.25);

    // hSec[pid]->Scale(0.25);
    // hSecNeg[pid]->Scale(0.25);
    // hSecPos[pid]->Scale(0.25);

    // hOut[pid]->Scale(0.25);

    hEff[pid]->GetXaxis()->SetRangeUser(ptmin, ptmax);
    hEff[pid]->GetYaxis()->SetRangeUser(effmin, effmax);
    hEffNeg[pid]->GetXaxis()->SetRangeUser(ptmin, ptmax);
    hEffNeg[pid]->GetYaxis()->SetRangeUser(effmin, effmax);
    hEffPos[pid]->GetXaxis()->SetRangeUser(ptmin, ptmax);
    hEffPos[pid]->GetYaxis()->SetRangeUser(effmin, effmax);
    hSec[pid]->GetXaxis()->SetRangeUser(ptmin, ptmax);
    hSec[pid]->GetYaxis()->SetRangeUser(secmin, secmax);
    hSecNeg[pid]->GetXaxis()->SetRangeUser(ptmin, ptmax);
    hSecNeg[pid]->GetYaxis()->SetRangeUser(secmin, secmax);
    hSecPos[pid]->GetXaxis()->SetRangeUser(ptmin, ptmax);
    hSecPos[pid]->GetYaxis()->SetRangeUser(secmin, secmax);
    hOut[pid]->GetXaxis()->SetRangeUser(ptmin, ptmax);
  }

  TCanvas* cChEff = new TCanvas("cChEff", "All efficiency", 800, 300);
  cChEff->Clear();
  cChEff->Divide(2, 1);
  cChEff->cd(1);
  //  TF1* pionEff = new TF1("pionEff", "pol0", 0.0, 50.0);
  TF1* chEff = new TF1("chEff", "[0]*(1-[1]/x)", 0.0, 50.0);
  chEff->SetLineColor(1);
  chEff->SetParameters(0.7, 0.5);
  hEff[0]->Fit(chEff, "", "", ptmin, ptmax);
  
  cChEff->cd(2);
  hEffPos[0]->Draw();
  hEffNeg[0]->Draw("SAME");
  chEff->DrawCopy("SAME");
  cChEff->SaveAs("eff_charged.gif");

  TCanvas* cPionEff = new TCanvas("cPionEff", "Pion efficiency", 800, 300);
  cPionEff->Clear();
  cPionEff->Divide(2, 1);
  cPionEff->cd(1);
  //  TF1* pionEff = new TF1("pionEff", "pol0", 0.0, 50.0);
  TF1* pionEff = new TF1("pionEff", "[0]*(1-[1]/x)", 0.0, 50.0);
  pionEff->SetLineColor(2);
  pionEff->SetParameters(0.7, 0.5);
  hEff[1]->Fit(pionEff, "", "", ptmin, ptmax);
  
  cPionEff->cd(2);
  hEffPos[1]->Draw();
  hEffNeg[1]->Draw("SAME");
  pionEff->DrawCopy("SAME");
  cPionEff->SaveAs("eff_pion.gif");

  TCanvas* cKaonEff = new TCanvas("cKaonEff", "Kaon efficiency", 800, 300);
  cKaonEff->Clear();
  cKaonEff->Divide(2, 1);
  cKaonEff->cd(1);
  //  TF1* kaonEff = new TF1("kaonEff", "pol0", 0.0, 50.0);
  TF1* kaonEff = new TF1("kaonEff", "[0]*(1-[1]/x)", 0.0, 50.0);
  kaonEff->SetLineColor(3);
  kaonEff->SetParameters(0.7, 0.5);
  hEff[2]->Fit(kaonEff, "", "", ptmin, ptmax);
  
  cKaonEff->cd(2);
  hEffPos[2]->Draw();
  hEffNeg[2]->Draw("SAME");
  kaonEff->DrawCopy("SAME");
  cKaonEff->SaveAs("eff_kaon.gif");

  TCanvas* cProtonEff = new TCanvas("cProtonEff", "Proton efficiency", 800, 300);
  cProtonEff->Clear();
  cProtonEff->Divide(2, 1);
  cProtonEff->cd(1);
  //  TF1* protonEff = new TF1("protonEff", "pol0", 0.0, 50.0);
  TF1* protonEff = new TF1("protonEff", "[0]*(1-[1]/x)", 0.0, 50.0);
  protonEff->SetLineColor(4);
  protonEff->SetParameters(0.7, 0.5);
  hEff[3]->Fit(protonEff, "", "", ptmin, ptmax);
  
  cProtonEff->cd(2);
  hEffPos[3]->Draw();
  hEffNeg[3]->Draw("SAME");
  protonEff->DrawCopy("SAME");
  cProtonEff->SaveAs("eff_proton.gif");

  TCanvas* cAllEff = new TCanvas("cAllEff", "All efficiency", 400, 300);
  cAllEff->Clear();
  
  TH1D* hDummy = (TH1D*)hEff[1]->Clone("hDummy");
  hDummy->Reset();
  hDummy->SetTitle("Efficiency vs p_{T}; p_{T} [GeV/c]; Efficiency");
  hDummy->Draw();
  
  chEff->DrawCopy("SAME");
  pionEff->DrawCopy("SAME");
  kaonEff->DrawCopy("SAME");
  protonEff->DrawCopy("SAME");
  cAllEff->SaveAs("eff_all.gif");

  // TCanvas* cPionSec = new TCanvas("cPionSec", "Pion sec", 800, 300);
  // cPionSec->Clear();
  // cPionSec->Divide(2, 1);
  // cPionSec->cd(1);
  // TF1* pionSec = new TF1("pionSec", "pol0", 0.0, 50.0);
  // hSec[1]->Fit(pionSec, "", "", ptmin, ptmax);
  
  // cPionSec->cd(2);
  // hSecPos[1]->Draw();
  // hSecNeg[1]->Draw("SAME");
  // pionSec->Draw("SAME");
  // cPionSec->SaveAs("sec_pion.gif");

  // TCanvas* cKaonSec = new TCanvas("cKaonSec", "Kaon sec", 800, 300);
  // cKaonSec->Clear();
  // cKaonSec->Divide(2, 1);
  // cKaonSec->cd(1);
  // TF1* kaonSec = new TF1("kaonSec", "pol0", 0.0, 50.0);
  // hSec[2]->Fit(kaonSec, "", "", ptmin, ptmax);
  
  // cKaonSec->cd(2);
  // hSecPos[2]->Draw();
  // hSecNeg[2]->Draw("SAME");
  // kaonSec->Draw("SAME");
  // cKaonSec->SaveAs("sec_kaon.gif");

  // TCanvas* cProtonSec = new TCanvas("cProtonSec", "Proton sec", 800, 300);
  // cProtonSec->Clear();
  // cProtonSec->Divide(2, 1);
  // cProtonSec->cd(1);
  // TF1* protonSec = new TF1("protonSec", "pol0", 0.0, 50.0);
  // hSec[3]->Fit(protonSec, "", "", ptmin, ptmax);
  
  // cProtonSec->cd(2);
  // hSecPos[3]->Draw();
  // hSecNeg[3]->Draw("SAME");
  // protonSec->Draw("SAME");
  // cProtonSec->SaveAs("sec_proton.gif");


  TCanvas* cEffRatioPi = new TCanvas("cEffRatioPi", "eff pi / eff all", 400, 300);
  cEffRatioPi->Clear();
  //  TF1* pionEff = new TF1("pionEff", "pol0", 0.0, 50.0);
  TF1* effRatioPi = new TF1("effRatioPi", "pol0", 0.0, 50.0);
  effRatioPi->SetLineColor(6);
  effRatioPi->SetParameters(0.7, 0.5);
  TH1D* hEffRatioPi = (TH1D*)hEff[1]->Clone("hEffRatioPi");
  hEffRatioPi->SetTitle("; p_{T} [GeV/c]; #epsilon_{ch}/#epsilon_{pi}");
  hEffRatioPi->Divide(hEff[0], hEff[1], 1, 1, "B");
  //  hEffRatioPi->Divide(hEff[1], hEff[0]);
  hEffRatioPi->GetXaxis()->SetRangeUser(0.0, 19.99);
  hEffRatioPi->GetYaxis()->SetRangeUser(0.8, 1.1);
  hEffRatioPi->SetStats(kTRUE);
  hEffRatioPi->Fit(effRatioPi, "", "", ptmin, ptmax);
  cEffRatioPi->SaveAs("eff_ratio_pi.gif");
  cEffRatioPi->SaveAs("eff_ratio_pi.pdf");

  TCanvas* cEffRatioK = new TCanvas("cEffRatioK", "eff K / eff all", 400, 300);
  cEffRatioK->Clear();
  TF1* effRatioK = new TF1("effRatioK", "exp(-[1]*x)+[0]", 0.0, 50.0);
  effRatioK->SetParameters(1.0, 1.0);
  effRatioK->SetLineColor(6);
  TH1D* hEffChargedRebinned = (TH1D*)hEff[0]->Clone("hEffChargedRebinned");
  hEffChargedRebinned->Rebin(2);
  TH1D* hEffRatioK = (TH1D*)hEff[2]->Clone("hEffRatioK");
  hEffRatioK->Rebin(2);
  hEffRatioK->SetTitle("; p_{T} [GeV/c]; #epsilon_{ch}/#epsilon_{K}");
  hEffRatioK->Divide(hEffChargedRebinned, hEffRatioK, 1, 1, "B");
  //  hEffRatioK->Divide(hEff[1], hEff[0]);
  hEffRatioK->GetXaxis()->SetRangeUser(0.0, 19.99);
  hEffRatioK->GetYaxis()->SetRangeUser(0.8, 1.1);
  hEffRatioK->SetStats(kTRUE);
  hEffRatioK->Fit(effRatioK, "", "", ptmin, ptmax);
  cEffRatioK->SaveAs("eff_ratio_K.gif");
  cEffRatioK->SaveAs("eff_ratio_K.pdf");
 
  TCanvas* cEffRatioP = new TCanvas("cEffRatioP", "eff p / eff all", 400, 300);
  cEffRatioP->Clear();
  TF1* effRatioP = new TF1("effRatioP", "pol0", 0.0, 50.0);
  effRatioP->SetLineColor(6);
  effRatioP->SetParameters(0.7, 0.5);
  // TH1D* hEffChargedRebinned = (TH1D*)hEff[0]->Clone("hEffChargedRebinned");
  // hEffChargedRebinned->Rebin(2);
  TH1D* hEffRatioP = (TH1D*)hEff[3]->Clone("hEffRatioP");
  hEffRatioP->Rebin(2);
  hEffRatioP->SetTitle("; p_{T} [GeV/c]; #epsilon_{ch}/#epsilon_{p}");
  hEffRatioP->Divide(hEffChargedRebinned, hEffRatioP, 1, 1, "B");
  //  hEffRatioP->Divide(hEff[1], hEff[0]);
  hEffRatioP->GetXaxis()->SetRangeUser(0.0, 19.99);
  hEffRatioP->GetYaxis()->SetRangeUser(0.8, 1.1);
  hEffRatioP->SetStats(kTRUE);
  hEffRatioP->Fit(effRatioP, "", "", ptmin, ptmax);
  cEffRatioP->SaveAs("eff_ratio_p.gif");
  cEffRatioP->SaveAs("eff_ratio_p.pdf");


  // TCanvas* cMuonRatio = new TCanvas("cMuonRatio", "muon / all", 400, 300);
  // cMuonRatio->Clear();
  // //  TF1* pionMuon = new TF1("pionMuon", "pol0", 0.0, 50.0);
  // TF1* muonRatio = new TF1("muonRatio", "pol0", 0.0, 50.0);
  // muonRatio->SetLineColor(1);
  // muonRatio->SetParameter(0, 0.01);
  // TH1D* hMuonRatio = (TH1D*)hOut[5]->Clone("hMuonRatio");
  // hMuonRatio->SetTitle("; p_{T} [GeV/c]; muon/all");
  // hMuonRatio->Divide(hOut[5], hOut[0], 1, 1, "B");
  // hMuonRatio->GetYaxis()->SetRangeUser(0.0, 0.05);
  // hMuonRatio->SetStats(kTRUE);
  // hMuonRatio->Fit(muonRatio, "", "", ptmin, ptmax);
  // cMuonRatio->SaveAs("muon_ratio.gif");
  // cMuonRatio->SaveAs("muon_ratio.pdf");

  // TCanvas* cElectronRatio = new TCanvas("cElectronRatio", "electron / all", 400, 300);
  // cElectronRatio->Clear();
  // //  TF1* pionElectron = new TF1("pionElectron", "pol0", 0.0, 50.0);
  // TF1* electronRatio = new TF1("electronRatio", "pol0", 0.0, 50.0);
  // electronRatio->SetLineColor(1);
  // electronRatio->SetParameter(0, 0.01);
  // TH1D* hElectronRatio = (TH1D*)hOut[4]->Clone("hElectronRatio");
  // hElectronRatio->SetTitle("; p_{T} [GeV/c]; electron/all");
  // hElectronRatio->Divide(hOut[4], hOut[0], 1, 1, "B");
  // hElectronRatio->GetYaxis()->SetRangeUser(0.0, 0.05);
  // hElectronRatio->SetStats(kTRUE);
  // hElectronRatio->SetMarkerColor(6);
  // hElectronRatio->Fit(electronRatio, "", "", ptmin, ptmax);
  // cElectronRatio->SaveAs("electron_ratio.gif");
  // cElectronRatio->SaveAs("electron_ratio.pdf");
  

  TFile* effFile = new TFile(effFileName, "RECREATE");
  chEff->Write();
  pionEff->Write();
  kaonEff->Write();
  protonEff->Write();
  effFile->Close();
}


//_______________________________________________________________________
void DrawCorrection(const Char_t* fileNamePYTHIA, const Char_t* fileNamePHOJET)
{
  TFile* filePYTHIA = FindFileFresh(fileNamePYTHIA);
  if(!filePYTHIA)
    return;

  TFile* filePHOJET = FindFileFresh(fileNamePHOJET);
  if(!filePHOJET)
    return;
  
  gStyle->SetOptStat(0);

  const Double_t ptmin  =  0.0;
  const Double_t ptmax  = 4.999;
  //  const Double_t ptmax  = 49.999;
  
  const Double_t effmin = 0.95;
  const Double_t effmax = 1.12;

  const Int_t nPid = 7;

  TH1D* hInPYTHIA[nPid]     = {0, 0, 0, 0, 0, 0, 0 };
  TH1D* hInPHOJET[nPid]     = {0, 0, 0, 0, 0, 0, 0 };

  for(Int_t pid = 0; pid < nPid; pid++) {

    hInPYTHIA[pid] = (TH1D*)filePYTHIA->Get(Form("hIn%d", pid));
    hInPYTHIA[pid]->GetXaxis()->SetRangeUser(ptmin, ptmax);
    hInPYTHIA[pid]->GetYaxis()->SetRangeUser(effmin, effmax);

    hInPHOJET[pid] = (TH1D*)filePHOJET->Get(Form("hIn%d", pid));
    hInPHOJET[pid]->GetXaxis()->SetRangeUser(ptmin, ptmax);
    hInPHOJET[pid]->GetYaxis()->SetRangeUser(effmin, effmax);
  }

  hInPYTHIA[1]->Add(hInPYTHIA[2]);
  hInPYTHIA[1]->Add(hInPYTHIA[3]);
  hInPYTHIA[0]->Divide(hInPYTHIA[1], hInPYTHIA[0], 1, 1, "B");

  TH1D* histPYTHIA = HistInvert(hInPYTHIA[0]);
  histPYTHIA->SetLineColor(2);
  histPYTHIA->SetMarkerColor(2);

  hInPHOJET[1]->Add(hInPHOJET[2]);
  hInPHOJET[1]->Add(hInPHOJET[3]);
  hInPHOJET[0]->Divide(hInPHOJET[1], hInPHOJET[0], 1, 1, "B");

  TH1D* histPHOJET = HistInvert(hInPHOJET[0]);
  histPHOJET->SetLineColor(4);
  histPHOJET->SetMarkerColor(4);

  TCanvas* cChCorrection = new TCanvas("cChCorrection", "All efficiency", 400, 300);
  cChCorrection->Clear();

  cChCorrection->SetGridy();

  histPYTHIA->GetYaxis()->SetRangeUser(effmin, effmax);
  histPYTHIA->SetTitle("N_{ch}/(#pi + K + p) for primaries at generator level; p_{T} [GeV/c]; Ratio");
  histPYTHIA->Draw();

  histPHOJET->Draw("SAME");

  TLegend* legend = new TLegend(0.55, 0.22, 0.79, 0.42);    
  legend->SetBorderSize(0);
  legend->SetFillColor(0);
  legend->AddEntry(histPYTHIA, "PYTHIA", "P");
  legend->AddEntry(histPHOJET, "PHOJET", "P");
  legend->Draw();
  
  cChCorrection->SaveAs("charged_over_piKp.gif");

  TFile* file = new TFile("correction.root", "RECREATE");
  histPYTHIA->SetName("histPYTHIA");
  histPHOJET->SetName("histPHOJET");
  histPYTHIA->Write();
  histPHOJET->Write();
  file->Close();
}

//_______________________________________________________________________
TH1D* HistInvert(TH1D* hist)
{
  TH1D* histNew = new TH1D(*hist);
  histNew->Reset();
  histNew->SetName(Form("%s_inv", hist->GetName()));
  histNew->SetTitle(Form("%s (inverted)", hist->GetTitle()));
  
  const Int_t nBins = hist->GetNbinsX();
  
  for(Int_t i = 1; i <= nBins; i++) {
    
    if(hist->GetBinContent(i) == 0)
      continue;

    histNew->SetBinContent(i, 1.0/hist->GetBinContent(i));
    histNew->SetBinError(i, hist->GetBinError(i)/hist->GetBinContent(i)/hist->GetBinContent(i));
  }

  return histNew;
}
