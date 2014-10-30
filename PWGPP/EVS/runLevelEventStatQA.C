#ifndef __CINT__
#include "TTree.h"
#include "TFile.h"
#include "TMath.h"
#endif
#include "triggerInfo.C"

// TODO read number of bits from AliVEvent?
#define NBITS 29

Int_t runLevelEventStatQA(TString qafilename="/data/alice/2010/LHC10b/000114783/pass4/QA_merge_archive.zip#event_stat.root", Int_t run=114783, TString ocdbStorage = "raw://"){
  printf("runLevelEventStatQA %s %i\n",qafilename.Data(),run);
  TFile* fin = new TFile(qafilename);
  
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
  if      (              run<=126437) { refSigma=  62; refClass = "CINT1B-ABCE-NOPF-ALL";    }
  else if (run>126437 && run<=127718) { refSigma=  62; refClass = "CINT1-B-NOPF-ALLNOTRD";   }
  else if (run>127718 && run<=127730) { refSigma=  62; refClass = "CINT1B-ABCE-NOPF-ALL";    }
  else if (run>127730 && run<=136848) { refSigma=  62; refClass = "CINT1-B-NOPF-ALLNOTRD";   }
  else if (run>136848 && run<=139517) { refSigma=7640; refClass = "CMBACS2-B-NOPF-ALLNOTRD"; }
  else if (run>166476 && run<=170593) { refSigma=4100; refClass = "CVLN-B-NOPF-ALLNOTRD";    }
  else if (run>195144 && run<=197388) { refSigma=1590; refClass = "C0TVX-B-NOPF-ALLNOTRD";   }

  if (refSigma>0) {
    Double_t par[5] = {0};
    triggerInfo(run,refClass,ocdbStorage,par);
    fill         = TMath::Nint(par[0]);
    duration     = par[1];
    l0b          = TMath::Nint(par[2]);
    nBCsPerOrbit = TMath::Nint(par[3]);
    mu           = par[4];
    lumi_seen    = l0b/refSigma;
  }

  for (Int_t j=1;j<=h->GetNbinsY();j++){
    TString label = h->GetYaxis()->GetBinLabel(j);

    // skip background triggers
    // TODO introduce identifier to filter-out background triggers
    if (!label.Contains("-B-") && !label.Contains("-S-") && !(label.Contains("-ABCE-") && label.Contains("1B-"))) continue;

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
  }
  
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
  t->Fill();
  
  TFile* fout = new TFile("trending.root","recreate");
  t->Write();
  fout->Close();
  return 0;
}
