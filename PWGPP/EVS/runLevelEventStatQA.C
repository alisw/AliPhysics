#ifndef __CINT__
#include "TTree.h"
#include "TFile.h"
#include "TMath.h"
#endif
#include "triggerInfo.C"

// TODO read number of bits from AliVEvent?
#define NBITS 29
TString bitNames[NBITS] = {
"kMB",
"kINT7",
"kMUON",
"kHighMult",
"kEMC1",
"kCINT5",
"kCMUS5",
"kMUSH7",
"kMUL7",
"kMUU7",
"kEMC7",
"kMUS7",
"kPHI1",
"kPHI7/kPHI8",
"kEMCEJE",
"kEMCEGA",
"kCentral",
"kSemiCentral",
"kDG5",
"kZED",
"kSPI7/kSPI8",
"kINT8",
"kMuonSingleLowPt8",
"kMuonSingleHighPt8",
"kMuonLikeLowPt8",
"kMuonUnlikeLowPt8",
"kMuonUnlikeLowPt0",
"kUserDefined",
"kTRD"
};


Int_t runLevelEventStatQA(TString qafilename="/data/alice/2010/LHC10b/000114783/pass4/QA_merge_archive.zip#event_stat.root", Int_t run=114783, TString ocdbStorage = "raw://"){
  printf("runLevelEventStatQA %s %i\n",qafilename.Data(),run);
  gStyle->SetOptStat(0);
  gStyle->SetLineScalePS(1.5);
  gStyle->SetPadBottomMargin(0.08);
  gStyle->SetPadRightMargin(0.02);
  gStyle->SetPadTopMargin(0.07);
  gStyle->SetPadLeftMargin(0.07);

  TFile* fin = new TFile(qafilename);
  fin->ls();
  fin->cd("trigger_histograms_+CINT1B-ABCE-NOPF-ALL,CINT1-B-NOPF-ALLNOTRD &1 *0");
  fin->ls();
//  TH1F* fHistV0A = (TH1F*) gDirectory->Get("fHistV0A");
//  fHistV0A->Draw();
//  return 0;
  TH2D* h = (TH2D*) fin->Get("fHistStatistics");
  if (!h) { printf("fHistStatistics not found\n"); return 1; }

  // tree variables
  Int_t all[NBITS]      = {0};
  Int_t accepted[NBITS] = {0};
  Int_t fill=0;
  Double_t duration=0;
  UInt_t l0b=0;
  Int_t nBCsPerOrbit=0;
  Double_t mu=0;
  Double_t lumi_seen=0;

  TString refClass="";
  Double_t refSigma=-1;
  if      (               run<=118501) { refSigma=  62.; refClass = "CINT1B-ABCE-NOPF-ALL";   } // pp_7.00: 62mb=54.3mb*1.15=sigma(VBAND)*R(INT1/VBAND) (Martino,2012-03-12,RunCond)
  else if (run>=118502 && run<=121040) { refSigma=  47.; refClass = "CINT1B-ABCE-NOPF-ALL";   } // pp_0.90: 47mb=52 mb *0.91=sigma(INEL)*R(INT1/INEL) (arxiv: 1208.4968, fig.10 + table 3)
  else if (run>=121041 && run<=126437) { refSigma=  62.; refClass = "CINT1B-ABCE-NOPF-ALL";   } // pp_7.00: 62mb=54.3mb*1.15=sigma(VBAND)*R(INT1/VBAND) (Martino,2012-03-12,RunCond)
  else if (run>=126438 && run<=127718) { refSigma=  62.; refClass = "CINT1-B-NOPF-ALLNOTRD";  } // pp_7.00: 62mb=54.3mb*1.15=sigma(VBAND)*R(INT1/VBAND) (Martino,2012-03-12,RunCond)
  else if (run>=127719 && run<=127730) { refSigma=  62.; refClass = "CINT1B-ABCE-NOPF-ALL";   } // pp_7.00: 62mb=54.3mb*1.15=sigma(VBAND)*R(INT1/VBAND) (Martino,2012-03-12,RunCond)
  else if (run>=127731 && run<=136848) { refSigma=  62.; refClass = "CINT1-B-NOPF-ALLNOTRD";  } // pp_7.00: 62mb=54.3mb*1.15=sigma(VBAND)*R(INT1/VBAND) (Martino,2012-03-12,RunCond)
  else if (run>=136849 && run<=139316) { refSigma=5970.; refClass = "C0SMH-B-NOPF-ALL";       } // PbPb_2.76: (Oyama,2011-05-20,RunCond)
  else if (run>=139328 && run<=139517) { refSigma=5970.; refClass = "C0SMH-B-NOPF-ALLNOTRD";  } // PbPb_2.76: (Oyama,2011-05-20,RunCond)
  else if (run>=145289 && run<=146860) { refSigma=  57.; refClass = "CINT1-B-NOPF-ALLNOTRD";  } // pp_2.76: 57mb=47.7mb*1.20=sigma(VBAND)*R(INT1/VBAND) (Martino,2012-03-12,RunCond)
  else if (run>=146808 && run<=146814) { refSigma=  57.; refClass = "CINT1-B-NOPF-ALL";       } // pp_2.76: 57mb=47.7mb*1.20=sigma(VBAND)*R(INT1/VBAND) (Martino,2012-03-12,RunCond)
  else if (run>=145815 && run<=146856) { refSigma=  57.; refClass = "CINT1-B-NOPF-ALLNOTRD";  } // pp_2.76: 57mb=47.7mb*1.20=sigma(VBAND)*R(INT1/VBAND) (Martino,2012-03-12,RunCond)
  else if (run>=146857 && run<=146857) { refSigma=  57.; refClass = "CINT1-B-NOPF-ALL";       } // pp_2.76: 57mb=47.7mb*1.20=sigma(VBAND)*R(INT1/VBAND) (Martino,2012-03-12,RunCond)
  else if (run>=146858 && run<=146860) { refSigma=  57.; refClass = "CINT1-B-NOPF-ALLNOTRD";  } // pp_2.76: 57mb=47.7mb*1.20=sigma(VBAND)*R(INT1/VBAND) (Martino,2012-03-12,RunCond)
  else if (run>=148370 && run<=157078) { refSigma=  54.; refClass = "CVBAND-B-NOPF-ALLNOTRD"; } // pp_7.00: 54.3mb (Martino,2012-03-12,RunCond)
  else if (run>=157079 && run<=165746) { refSigma=  24.; refClass = "C0TVX-B-NOPF-ALLNOTRD";  } // pp_7.00: 24mb=54.3mb*0.44=sigma(VBAND)*R(0TVX/VBAND) (Martino,2012-03-12,RunCond)
  else if (run>=166477 && run<=170593) { refSigma=4100.; refClass = "CVLN-B-NOPF-ALLNOTRD";   } // PbPb_2.76: (Martino,2013-03-15,RunCond)
  else if (run>=176658 && run<=177143) { refSigma=  25.; refClass = "C0TVX-B-NOPF-ALLNOTRD";  } // pp_8.00: (Artem, 2013-10-04,RunCond)
  else if (run>=177146 && run<=177147) { refSigma=  25.; refClass = "C0TVX-B-NOPF-CENTNOTRD"; } // pp_8.00: (Artem, 2013-10-04,RunCond)
  else if (run>=177148 && run<=177149) { refSigma=  25.; refClass = "C0TVX-B-NOPF-ALLNOTRD";  } // pp_8.00: (Artem, 2013-10-04,RunCond)
  else if (run>=177150 && run<=177506) { refSigma=  25.; refClass = "C0TVX-B-NOPF-CENTNOTRD"; } // pp_8.00: (Artem, 2013-10-04,RunCond)
  else if (run>=177580 && run<=178220) { refSigma=  25.; refClass = "C0TVX-B-NOPF-ALLNOTRD";  } // pp_8.00: (Artem, 2013-10-04,RunCond)
  else if (run>=179444 && run<=193692) { refSigma=  25.; refClass = "C0TVX-S-NOPF-ALLNOTRD";  } // pp_8.00: (Artem, 2013-10-04,RunCond)
  else if (run>=193693 && run<=193766) { refSigma=  25.; refClass = "C0TVX-B-NOPF-ALLNOTRD";  } // pp_8.00: (Artem, 2013-10-04,RunCond)
  else if (run>=195344 && run<=197388) { refSigma=1590.; refClass = "C0TVX-B-NOPF-ALLNOTRD";  } // pPb_5.02: arxiv:1405.1849
  else if (run>=197470 && run<=197692) { refSigma=  18.; refClass = "C0TVX-B-NOPF-ALLNOTRD";  } // pp_2.76: 18mb=47.7mb*0.39=sigma(VBAND)*R(0TVX/VBAND) (Martino,2012-03-12,RunCond)
  Double_t par[5] = {0};
  triggerInfo(run,refClass,ocdbStorage,par);
  fill         = TMath::Nint(par[0]);
  duration     = par[1];
  l0b          = TMath::Nint(par[2]);
  nBCsPerOrbit = TMath::Nint(par[3]);
  mu           = par[4];
  if (refSigma>0) lumi_seen = l0b/refSigma;

  TTree* t = new TTree("trending","tree of trending variables");
  t->Branch("run",&run);
  t->Branch("fill",&fill);
  t->Branch("bcs",&nBCsPerOrbit);
  t->Branch("duration",&duration);
  t->Branch("mu",&mu);
  t->Branch("l0b",&l0b);
  t->Branch("lumi_seen",&lumi_seen);
  t->Branch("all",&all,Form("all[%i]/I",NBITS));
  t->Branch("accepted",&accepted,Form("accepted[%i]/I",NBITS));

  TFile* fout = new TFile("trending.root","recreate");

  for (Int_t j=1;j<=h->GetNbinsY();j++){
    TString label = h->GetYaxis()->GetBinLabel(j);

    // skip background triggers
    // TODO introduce identifier to filter-out background triggers
    if (!label.Contains("-B-") && !label.Contains("-S-") && !(label.Contains("-ABCE-") && label.Contains("1B-"))) continue;

    printf("%s\n",label.Data());
    // Read mask
    // TODO think how to propagate mask with TBit aliases
    UInt_t mask = 0;
    TObjArray* array = label.Tokenize(" ");
    for (Int_t itoken=0;itoken<array->GetEntries();itoken++){
      TString token = array->At(itoken)->GetName();
      if (token[0]!='&') continue;
      token.Remove(0,1);
      mask = token.Atoi();
      break;
    }
    array->Delete();
    delete array;
    if (!mask) continue;
    // Fill all and accepted counters for active bits
    // TODO can we accidentally double count events?
    for (Int_t ibit=0;ibit<NBITS;ibit++) {
      if (!(mask & 1<<ibit)) continue; // to be changed with TBits
      all[ibit]      += Int_t(h->GetBinContent(1             ,j)); 
      accepted[ibit] += Int_t(h->GetBinContent(h->GetNbinsX(),j));
    }
    if (h->GetBinContent(1,j)<1) continue;
    TH1F* fHistV0A = (TH1F*) fin->Get(Form("trigger_histograms_%s/fHistV0A",label.Data()));
    TH1F* fHistV0C = (TH1F*) fin->Get(Form("trigger_histograms_%s/fHistV0C",label.Data()));
    TH1F* fHistFiredBitsSPD = (TH1F*) fin->Get(Form("trigger_histograms_%s/fHistFiredBitsSPD",label.Data()));
    TH2F* fHistBitsSPD      = (TH2F*) fin->Get(Form("trigger_histograms_%s/fHistBitsSPD",label.Data()));
    TH1F* fHistTDCZDC       = (TH1F*) fin->Get(Form("trigger_histograms_%s/fHistTDCZDC",label.Data()));
    TH2F* fHistTimeZDC      = (TH2F*) fin->Get(Form("trigger_histograms_%s/fHistTimeZDC",label.Data()));
    TH2F* fHistTimeCorrZDC  = (TH2F*) fin->Get(Form("trigger_histograms_%s/fHistTimeCorrZDC",label.Data()));

    const char* bitName = bitNames[j-1].Data();
    TCanvas* cHistV0A = new TCanvas(Form("cHistV0A",j),Form("cHistV0A_%s",bitName),1000,800);
    gPad->SetLogy();
    fHistV0A->SetTitle(Form("%s: V0A",bitName));
    fHistV0A->SetLineWidth(2);
    fHistV0A->SetLineColor(kBlue);
    fHistV0A->Draw();
    gPad->Print(Form("%s_histV0A.pdf",bitName));
    fHistV0A->Write(Form("%s_histV0A",bitName));

    TCanvas* cHistV0C = new TCanvas(Form("cHistV0C",j),Form("cHistV0C_%s",bitName),1000,800);
    gPad->SetLogy();
    fHistV0C->SetTitle(Form("%s: V0C",bitName));
    fHistV0C->SetLineWidth(2);
    fHistV0C->SetLineColor(kBlue);
    fHistV0C->Draw();
    gPad->Print(Form("%s_histV0C.pdf",bitName));
    fHistV0C->Write(Form("%s_histV0C",bitName));

    TCanvas* cHistFiredBitsSPD = new TCanvas(Form("cHistFiredBitsSPD",j),Form("cHistFiredBitsSPD_%s",bitName),1800,500);
    gPad->SetLogy();
    gPad->SetMargin(0.05,0.01,0.12,0.06);
    fHistFiredBitsSPD->SetTitle(Form("%s: hardware FO",bitName));
    fHistFiredBitsSPD->SetTitleFont(43);
    fHistFiredBitsSPD->SetTitleSize(25);
    fHistFiredBitsSPD->GetYaxis()->SetTitleFont(43);
    fHistFiredBitsSPD->GetXaxis()->SetLabelFont(43);
    fHistFiredBitsSPD->GetYaxis()->SetLabelFont(43);
    fHistFiredBitsSPD->GetYaxis()->SetTitleSize(25);
    fHistFiredBitsSPD->GetXaxis()->SetLabelSize(25);
    fHistFiredBitsSPD->GetYaxis()->SetLabelSize(25);
    fHistFiredBitsSPD->GetYaxis()->SetTickLength(0.01);
    fHistFiredBitsSPD->GetYaxis()->SetTitleOffset(0.5);
    fHistFiredBitsSPD->GetYaxis()->SetDecimals(1);
    fHistFiredBitsSPD->SetLineWidth(2);
    fHistFiredBitsSPD->SetLineColor(kBlue);
    fHistFiredBitsSPD->Draw();
    gPad->Print(Form("%s_histFiredBitsSPD.pdf",bitName));
    fHistFiredBitsSPD->Write(Form("%s_histFiredBitsSPD",bitName));

    TCanvas* cHistBitsSPD = new TCanvas(Form("cHistBitsSPD",j),Form("cHistBitsSPD_%s",bitName),800,800);
    gPad->SetLogz();
    gPad->SetMargin(0.12,0.12,0.10,0.06);
    fHistBitsSPD->SetTitle(Form("%s: hardware FO vs offline FO",bitName));
    fHistBitsSPD->GetXaxis()->SetTitleOffset(1.3);
    fHistBitsSPD->GetYaxis()->SetTitleOffset(1.6);
    fHistBitsSPD->Draw("colz");
    gPad->Print(Form("%s_histBitsSPD.pdf",bitName));
    fHistBitsSPD->Write(Form("%s_histBitsSPD",bitName));

    TCanvas* cHistTimeZDC = new TCanvas(Form("cHistTimeZDC",j),Form("cHistTimeZDC_%s",bitName),800,800);
    gPad->SetLogz();
    gPad->SetMargin(0.12,0.12,0.10,0.06);
    fHistTimeZDC->SetTitle(Form("%s: ZDC timing;TDC timing C-A;TDC timing C+A",bitName));
    fHistTimeZDC->GetXaxis()->SetTitleOffset(1.3);
    fHistTimeZDC->GetYaxis()->SetTitleOffset(1.6);
    fHistTimeZDC->Draw("colz");
    gPad->Print(Form("%s_histTimeZDC.pdf",bitName));
    fHistTimeZDC->Write(Form("%s_histTimeZDC",bitName));

    TCanvas* cHistTimeCorrZDC = new TCanvas(Form("cHistTimeCorrZDC",j),Form("cHistTimeCorrZDC_%s",bitName),800,800);
    gPad->SetLogz();
    gPad->SetMargin(0.12,0.12,0.10,0.06);
    fHistTimeCorrZDC->SetTitle(Form("%s: corrected ZDC timing;TDC timing C-A;TDC timing C+A",bitName));
    fHistTimeCorrZDC->GetXaxis()->SetTitleOffset(1.3);
    fHistTimeCorrZDC->GetYaxis()->SetTitleOffset(1.6);
    fHistTimeCorrZDC->Draw("colz");
    gPad->Print(Form("%s_histTimeCorrZDC.pdf",bitName));
    fHistTimeCorrZDC->Write(Form("%s_histTimeCorrZDC",bitName));
  }
  
  t->Fill();
  t->Write();
  fout->Close();
  return 0;
}
