#if !defined (__CINT__) || (defined(__MAKECINT__))
#include <iostream>
#include "TClonesArray.h"
#include "TH1.h"
#include "AliPWGFunc.h"
//#include "LoadSpectraPiKPPbPb.C"
#include "AliParticleYield.h"
#include "TFile.h"
#include "TDatabasePDG.h"

#endif


void AddHistograms (TH1 * hdest, TH1 * hsource, Float_t scale = 1. , Int_t isLinear = 0) ;
TH1* PrintYieldAndError(TH1 * hsyst, TH1* hstat, Int_t ifunc = 0, Double_t massParticle=0.139) ;
TH1* ReweightSpectra(TClonesArray * histos, Double_t * weights, Int_t isLinear = 0, TString suffix = "") ;
void LoadStuff(Bool_t loadPbPb = kTRUE) ;
void AverageAndExtrapolate(TString what);



TH1 * hLambdaStat[7];
TH1 * hLambdaSyst[7];
TH1 * hK0SStat[7]   ;
TH1 * hK0SSyst[7]   ;



void AverageAndExtrapolate () {

  // MF COMMON LIBRARIES SHOULD NOT BE LOADED FOR THIS MACRO TO WORK
  LoadStuff(1); // True PbPb, False pPb

  //  PrintYieldAndError(hSpectraCentr_sys[0][kMyProton][0], hSpectraCentr_stat[0][kMyProton][0], 0, mass[2]);

  //  PrintYieldAndError(hLambdaSyst[0], hLambdaStat[0], 0, TDatabasePDG::Instance()->GetParticle("Lambda0")->Mass());
  // icentr, ipart, icharge

  //  AverageAndExtrapolate("pions_pos_0020");
  //  AverageAndExtrapolate("pions_pos_0010");
  // AverageAndExtrapolate("pions_neg_0010");
  //  AverageAndExtrapolate("pions_sum_0010");
  // AverageAndExtrapolate("kaons_pos_0010");
  // AverageAndExtrapolate("kaons_neg_0010");
  // AverageAndExtrapolate("protons_pos_0010");
  // AverageAndExtrapolate("protons_neg_0010");
  //  AverageAndExtrapolate("lambda_0010");
  //  AverageAndExtrapolate("k0s_0010");
  // AverageAndExtrapolate("pions_pos_6080");
  // AverageAndExtrapolate("pions_neg_6080");
  // AverageAndExtrapolate("kaons_pos_6080");
  // AverageAndExtrapolate("kaons_neg_6080");
  // AverageAndExtrapolate("protons_pos_6080");
  // AverageAndExtrapolate("protons_neg_6080");

  //  AverageAndExtrapolate("pions_pos_2040");
  AverageAndExtrapolate("pions_neg_2040");
  AverageAndExtrapolate("kaons_pos_2040");
  AverageAndExtrapolate("kaons_neg_2040");
  AverageAndExtrapolate("protons_pos_2040");
  AverageAndExtrapolate("protons_neg_2040");
  
}

TH1* PrintYieldAndError(TH1 * hsyst, TH1* hstat, Int_t ifunc, Double_t massParticle) {

  TF1 * f = 0;
  if (ifunc == 0) {
    f= BGBlastWave("fBlastWave",  massParticle);
    f->SetParameters(massParticle, 0.640, 0.097, 0.73, 100); // FIXME
  }
  else if (ifunc == 1) {
    f  = fm.GetTsallis(massParticle, 0.2, 11, 800, "fTsallisPion");
  }
  
  //  TH1 * h =  YieldMean(hstat, hsyst, f, 0.0, 100, 0.01, 0.1, "");
  TH1 * h =  YieldMean(hstat, hsyst, f, 0.0, 100);
  //  std::cout << "H " << h << std::endl;
  std::cout << "" << std::endl;
  std::cout << h->GetBinContent(1) << ", " << h->GetBinContent(2) << ", " << h->GetBinContent(3) << std::endl;
  std::cout << "" << std::endl;
  return h;


}


TH1* ReweightSpectra(TObjArray * histos, Double_t * weights, Int_t isLinear, TString suffix) {

  // sums a number of spectra with a given weight. Used to combine
  // several centrality bins.  The weights should add up to 1 and be
  // proportional to the width of the centrality bin for each
  // histogram
  // Suffix is added to the name of the histogram.
  // if linear = 1 errors are summed linearly rather than in quadrature


  TIter it(histos);
  TH1 * h = 0;
  TH1 * hsum = 0;
  Int_t ihisto = 0;
  while ((h = dynamic_cast<TH1*>(it.Next()))) {
    if(!h) {
      std::cout << "ERROR cannot get one of the histos!" << std::endl;
      return 0;
    }

    if (!hsum) {
      // First histogram, clone it
      hsum =  (TH1D*) h->Clone(TString(h->GetName())+suffix);
      hsum->Scale(weights[ihisto++]);
    }else {
      AddHistograms(hsum, h, weights[ihisto++], isLinear);
    }
  }
  return hsum;
}

void AddHistograms (TH1 * hdest, TH1 * hsource, Float_t scale, Int_t isLinear) {
  // THis method assumes that the 2 histos are consistent!!
  TH1 * hsourceLoc = (TH1*) hsource->Clone("hsourceLoc");
  hsourceLoc->Scale(scale);
  if(!isLinear) {
    hdest->Add(hsourceLoc);
  } 
  else {
    Int_t nbin = hsourceLoc->GetNbinsX();
    for(Int_t ibin = 0; ibin < nbin; ibin++){
      Double_t content = hdest->GetBinContent(ibin) + hsourceLoc->GetBinContent(ibin);
      Double_t error   = hdest->GetBinError  (ibin) + hsourceLoc->GetBinError  (ibin);
      hdest->SetBinContent(ibin,content);
      hdest->SetBinError(ibin, error);
    }
    

  }

  delete hsourceLoc;
  
  

}


void LoadStuff(Bool_t loadPbPb) {
  if(loadPbPb){
    gROOT->LoadMacro("LoadSpectraPiKPPbPb.C");
    LoadSpectraPiKPPbPb();
  }
  else {
    gROOT->LoadMacro("LoadSpectraPiKPProtonLead.C");
    LoadSpectraPiKPProtonLead();
  }
  LoadLibs();
  gROOT->LoadMacro("$ALICE_ROOT/PWGLF/SPECTRA/UTILS/YieldMean.C");
  gROOT->LoadMacro("$ALICE_ROOT/PWGLF/SPECTRA/UTILS/SpectraUtils.C");

  // Load Lambdas and K0s
  TFile * f = new TFile("k0s_lambda_final_spectra.root");
  const char * multTags[] = {  "0005",  "0510",  "1020",  "2040",  "4060",  "6080", "8090"};
  for (Int_t icentr = 0; icentr<7; icentr++) {
    hLambdaStat[icentr] = (TH1*) f->Get(Form("statonly_cent%s_Lambda",multTags[icentr]));
    hLambdaSyst[icentr] = (TH1*) f->Get(Form("statonly_cent%s_Lambda",multTags[icentr]));
    hK0SStat[icentr]    = (TH1*) f->Get(Form("systonly_cent%s_K0s",multTags[icentr]));
    hK0SSyst[icentr]    = (TH1*) f->Get(Form("statonly_cent%s_K0s",multTags[icentr]));

    // The bin 300-400 MeV was not used in the analysis
    hK0SStat[icentr]->SetBinContent(4,0);   
    hK0SSyst[icentr]->SetBinContent(4,0);
    hK0SStat[icentr]->SetBinError(4,0);   
    hK0SSyst[icentr]->SetBinError(4,0);
  }

}

void AverageAndExtrapolate(TString what) {

  TH1 * hstat =0;
  TH1 * hsyst =0;
  TH1 * hyieldmean =0;

  if(what == "pions_pos_0020"){ 

    TObjArray * arrStat = new TObjArray();
    arrStat->Add(hSpectraCentr_stat[0][kMyPion][kMyPos]);
    arrStat->Add(hSpectraCentr_stat[1][kMyPion][kMyPos]);
    arrStat->Add(hSpectraCentr_stat[2][kMyPion][kMyPos]);
    
    TObjArray * arrSyst = new TObjArray();
    arrSyst->Add(hSpectraCentr_sys [0][kMyPion][kMyPos]);
    arrSyst->Add(hSpectraCentr_sys [1][kMyPion][kMyPos]);
    arrSyst->Add(hSpectraCentr_sys [2][kMyPion][kMyPos]);
    
    Double_t weights[] = {0.25, 0.25, 0.5};
    //Double_t weights[] = {0.5, 0.5};
    
    TH1 * hstat = ReweightSpectra(arrStat, weights, 0, "_0020");
    TH1 * hsyst = ReweightSpectra(arrSyst, weights, 1, "_0020");

    hyieldmean = PrintYieldAndError(hsyst, hstat, 0, mass[0]);

    
  }
  if(what == "pions_pos_0010"){ 

    TObjArray * arrStat = new TObjArray();
    arrStat->Add(hSpectraCentr_stat[0][kMyPion][kMyPos]);
    arrStat->Add(hSpectraCentr_stat[1][kMyPion][kMyPos]);
    
    TObjArray * arrSyst = new TObjArray();
    arrSyst->Add(hSpectraCentr_sys [0][kMyPion][kMyPos]);
    arrSyst->Add(hSpectraCentr_sys [1][kMyPion][kMyPos]);
    
    Double_t weights[] = {0.5, 0.5};
    
    TH1 * hstat = ReweightSpectra(arrStat, weights, 0, "_0010");
    TH1 * hsyst = ReweightSpectra(arrSyst, weights, 1, "_0010");

    hyieldmean = PrintYieldAndError(hsyst, hstat, 0, mass[0]);

    
  }
  if(what == "pions_sum_0010"){ 

    TObjArray * arrStat = new TObjArray();
    arrStat->Add(hSpectraCentr_stat[0][kMyPion][kMyPos]);
    arrStat->Add(hSpectraCentr_stat[1][kMyPion][kMyPos]);
    arrStat->Add(hSpectraCentr_stat[0][kMyPion][kMyNeg]);
    arrStat->Add(hSpectraCentr_stat[1][kMyPion][kMyNeg]);
    
    TObjArray * arrSyst = new TObjArray();
    arrSyst->Add(hSpectraCentr_sys [0][kMyPion][kMyPos]);
    arrSyst->Add(hSpectraCentr_sys [1][kMyPion][kMyPos]);
    arrSyst->Add(hSpectraCentr_sys [0][kMyPion][kMyNeg]);
    arrSyst->Add(hSpectraCentr_sys [1][kMyPion][kMyNeg]);
    
    Double_t weights[] = {0.5, 0.5, 0.5, 0.5};
    
    TH1 * hstat = ReweightSpectra(arrStat, weights, 0, "_0010");
    TH1 * hsyst = ReweightSpectra(arrSyst, weights, 1, "_0010");

    hyieldmean = PrintYieldAndError(hsyst, hstat, 0, mass[0]);

    
  }
  if(what == "pions_neg_0010"){ 

    TObjArray * arrStat = new TObjArray();
    arrStat->Add(hSpectraCentr_stat[0][kMyPion][kMyNeg]);
    arrStat->Add(hSpectraCentr_stat[1][kMyPion][kMyNeg]);
    
    TObjArray * arrSyst = new TObjArray();
    arrSyst->Add(hSpectraCentr_sys [0][kMyPion][kMyNeg]);
    arrSyst->Add(hSpectraCentr_sys [1][kMyPion][kMyNeg]);
    
    Double_t weights[] = {0.5, 0.5};
    
    TH1 * hstat = ReweightSpectra(arrStat, weights, 0, "_0010");
    TH1 * hsyst = ReweightSpectra(arrSyst, weights, 1, "_0010");

    hyieldmean = PrintYieldAndError(hsyst, hstat, 0, mass[0]);

    
  }
  if(what == "protons_pos_0010"){ 

    TObjArray * arrStat = new TObjArray();
    arrStat->Add(hSpectraCentr_stat[0][kMyProton][kMyPos]);
    arrStat->Add(hSpectraCentr_stat[1][kMyProton][kMyPos]);
    
    TObjArray * arrSyst = new TObjArray();
    arrSyst->Add(hSpectraCentr_sys [0][kMyProton][kMyPos]);
    arrSyst->Add(hSpectraCentr_sys [1][kMyProton][kMyPos]);
    
    Double_t weights[] = {0.5, 0.5};
    
    TH1 * hstat = ReweightSpectra(arrStat, weights, 0, "_0010");
    TH1 * hsyst = ReweightSpectra(arrSyst, weights, 1, "_0010");

    hyieldmean = PrintYieldAndError(hsyst, hstat, 0, mass[2]);

    
  }
  if(what == "protons_neg_0010"){ 

    TObjArray * arrStat = new TObjArray();
    arrStat->Add(hSpectraCentr_stat[0][kMyProton][kMyNeg]);
    arrStat->Add(hSpectraCentr_stat[1][kMyProton][kMyNeg]);
    
    TObjArray * arrSyst = new TObjArray();
    arrSyst->Add(hSpectraCentr_sys [0][kMyProton][kMyNeg]);
    arrSyst->Add(hSpectraCentr_sys [1][kMyProton][kMyNeg]);
    
    Double_t weights[] = {0.5, 0.5};
    
    TH1 * hstat = ReweightSpectra(arrStat, weights, 0, "_0010");
    TH1 * hsyst = ReweightSpectra(arrSyst, weights, 1, "_0010");

    hyieldmean = PrintYieldAndError(hsyst, hstat, 0, mass[2]);

    
  }
    if(what == "kaons_pos_0010"){ 

    TObjArray * arrStat = new TObjArray();
    arrStat->Add(hSpectraCentr_stat[0][kMyKaon][kMyPos]);
    arrStat->Add(hSpectraCentr_stat[1][kMyKaon][kMyPos]);
    
    TObjArray * arrSyst = new TObjArray();
    arrSyst->Add(hSpectraCentr_sys [0][kMyKaon][kMyPos]);
    arrSyst->Add(hSpectraCentr_sys [1][kMyKaon][kMyPos]);
    
    Double_t weights[] = {0.5, 0.5};
    
    TH1 * hstat = ReweightSpectra(arrStat, weights, 0, "_0010");
    TH1 * hsyst = ReweightSpectra(arrSyst, weights, 1, "_0010");

    hyieldmean = PrintYieldAndError(hsyst, hstat, 0, mass[1]);

    
  }
  if(what == "kaons_neg_0010"){ 

    TObjArray * arrStat = new TObjArray();
    arrStat->Add(hSpectraCentr_stat[0][kMyKaon][kMyNeg]);
    arrStat->Add(hSpectraCentr_stat[1][kMyKaon][kMyNeg]);
    
    TObjArray * arrSyst = new TObjArray();
    arrSyst->Add(hSpectraCentr_sys [0][kMyKaon][kMyNeg]);
    arrSyst->Add(hSpectraCentr_sys [1][kMyKaon][kMyNeg]);
    
    Double_t weights[] = {0.5, 0.5};
    
    TH1 * hstat = ReweightSpectra(arrStat, weights, 0, "_0010");
    TH1 * hsyst = ReweightSpectra(arrSyst, weights, 1, "_0010");

    hyieldmean = PrintYieldAndError(hsyst, hstat, 0, mass[1]);

    
  }

  if(what == "pions_pos_6080"){ 

    TObjArray * arrStat = new TObjArray();
    arrStat->Add(hSpectraCentr_stat[7][kMyPion][kMyPos]);
    arrStat->Add(hSpectraCentr_stat[8][kMyPion][kMyPos]);
    
    TObjArray * arrSyst = new TObjArray();
    arrSyst->Add(hSpectraCentr_sys [7][kMyPion][kMyPos]);
    arrSyst->Add(hSpectraCentr_sys [8][kMyPion][kMyPos]);
    
    Double_t weights[] = {0.5, 0.5};
    
    TH1 * hstat = ReweightSpectra(arrStat, weights, 0, "_6080");
    TH1 * hsyst = ReweightSpectra(arrSyst, weights, 1, "_6080");

    hyieldmean = PrintYieldAndError(hsyst, hstat, 0, mass[0]);

    
  }
  if(what == "pions_neg_6080"){ 

    TObjArray * arrStat = new TObjArray();
    arrStat->Add(hSpectraCentr_stat[7][kMyPion][kMyNeg]);
    arrStat->Add(hSpectraCentr_stat[8][kMyPion][kMyNeg]);
    
    TObjArray * arrSyst = new TObjArray();
    arrSyst->Add(hSpectraCentr_sys [7][kMyPion][kMyNeg]);
    arrSyst->Add(hSpectraCentr_sys [8][kMyPion][kMyNeg]);
    
    Double_t weights[] = {0.5, 0.5};
    
    TH1 * hstat = ReweightSpectra(arrStat, weights, 0, "_6080");
    TH1 * hsyst = ReweightSpectra(arrSyst, weights, 1, "_6080");

    hyieldmean = PrintYieldAndError(hsyst, hstat, 0, mass[0]);

    
  }
  if(what == "protons_pos_6080"){ 

    TObjArray * arrStat = new TObjArray();
    arrStat->Add(hSpectraCentr_stat[7][kMyProton][kMyPos]);
    arrStat->Add(hSpectraCentr_stat[8][kMyProton][kMyPos]);
    
    TObjArray * arrSyst = new TObjArray();
    arrSyst->Add(hSpectraCentr_sys [7][kMyProton][kMyPos]);
    arrSyst->Add(hSpectraCentr_sys [8][kMyProton][kMyPos]);
    
    Double_t weights[] = {0.5, 0.5};
    
    TH1 * hstat = ReweightSpectra(arrStat, weights, 0, "_6080");
    TH1 * hsyst = ReweightSpectra(arrSyst, weights, 1, "_6080");

    hyieldmean = PrintYieldAndError(hsyst, hstat, 0, mass[2]);

    
  }
  if(what == "protons_neg_6080"){ 

    TObjArray * arrStat = new TObjArray();
    arrStat->Add(hSpectraCentr_stat[7][kMyProton][kMyNeg]);
    arrStat->Add(hSpectraCentr_stat[8][kMyProton][kMyNeg]);
    
    TObjArray * arrSyst = new TObjArray();
    arrSyst->Add(hSpectraCentr_sys [7][kMyProton][kMyNeg]);
    arrSyst->Add(hSpectraCentr_sys [8][kMyProton][kMyNeg]);
    
    Double_t weights[] = {0.5, 0.5};
    
    TH1 * hstat = ReweightSpectra(arrStat, weights, 0, "_6080");
    TH1 * hsyst = ReweightSpectra(arrSyst, weights, 1, "_6080");

    hyieldmean = PrintYieldAndError(hsyst, hstat, 0, mass[2]);

    
  }
    if(what == "kaons_pos_6080"){ 

    TObjArray * arrStat = new TObjArray();
    arrStat->Add(hSpectraCentr_stat[7][kMyKaon][kMyPos]);
    arrStat->Add(hSpectraCentr_stat[8][kMyKaon][kMyPos]);
    
    TObjArray * arrSyst = new TObjArray();
    arrSyst->Add(hSpectraCentr_sys [7][kMyKaon][kMyPos]);
    arrSyst->Add(hSpectraCentr_sys [8][kMyKaon][kMyPos]);
    
    Double_t weights[] = {0.5, 0.5};
    
    TH1 * hstat = ReweightSpectra(arrStat, weights, 0, "_6080");
    TH1 * hsyst = ReweightSpectra(arrSyst, weights, 1, "_6080");

    hyieldmean = PrintYieldAndError(hsyst, hstat, 0, mass[1]);

    
  }
  if(what == "kaons_neg_6080"){ 

    TObjArray * arrStat = new TObjArray();
    arrStat->Add(hSpectraCentr_stat[7][kMyKaon][kMyNeg]);
    arrStat->Add(hSpectraCentr_stat[8][kMyKaon][kMyNeg]);
    
    TObjArray * arrSyst = new TObjArray();
    arrSyst->Add(hSpectraCentr_sys [7][kMyKaon][kMyNeg]);
    arrSyst->Add(hSpectraCentr_sys [8][kMyKaon][kMyNeg]);
    
    Double_t weights[] = {0.5, 0.5};
    
    TH1 * hstat = ReweightSpectra(arrStat, weights, 0, "_6080");
    TH1 * hsyst = ReweightSpectra(arrSyst, weights, 1, "_6080");

    hyieldmean = PrintYieldAndError(hsyst, hstat, 0, mass[1]);

    
  }

  if(what == "pions_pos_2040"){ 

    TObjArray * arrStat = new TObjArray();
    arrStat->Add(hSpectraCentr_stat[3][kMyPion][kMyPos]);
    arrStat->Add(hSpectraCentr_stat[4][kMyPion][kMyPos]);
    
    TObjArray * arrSyst = new TObjArray();
    arrSyst->Add(hSpectraCentr_sys [3][kMyPion][kMyPos]);
    arrSyst->Add(hSpectraCentr_sys [4][kMyPion][kMyPos]);
    
    Double_t weights[] = {0.5, 0.5};
    
    TH1 * hstat = ReweightSpectra(arrStat, weights, 0, "_2040");
    TH1 * hsyst = ReweightSpectra(arrSyst, weights, 1, "_2040");

    hyieldmean = PrintYieldAndError(hsyst, hstat, 0, mass[0]);

    
  }
  if(what == "pions_neg_2040"){ 

    TObjArray * arrStat = new TObjArray();
    arrStat->Add(hSpectraCentr_stat[3][kMyPion][kMyNeg]);
    arrStat->Add(hSpectraCentr_stat[4][kMyPion][kMyNeg]);
    
    TObjArray * arrSyst = new TObjArray();
    arrSyst->Add(hSpectraCentr_sys [3][kMyPion][kMyNeg]);
    arrSyst->Add(hSpectraCentr_sys [4][kMyPion][kMyNeg]);
    
    Double_t weights[] = {0.5, 0.5};
    
    TH1 * hstat = ReweightSpectra(arrStat, weights, 0, "_2040");
    TH1 * hsyst = ReweightSpectra(arrSyst, weights, 1, "_2040");

    hyieldmean = PrintYieldAndError(hsyst, hstat, 0, mass[0]);

    
  }
    if(what == "kaons_pos_2040"){ 

    TObjArray * arrStat = new TObjArray();
    arrStat->Add(hSpectraCentr_stat[3][kMyKaon][kMyPos]);
    arrStat->Add(hSpectraCentr_stat[4][kMyKaon][kMyPos]);
    
    TObjArray * arrSyst = new TObjArray();
    arrSyst->Add(hSpectraCentr_sys [3][kMyKaon][kMyPos]);
    arrSyst->Add(hSpectraCentr_sys [4][kMyKaon][kMyPos]);
    
    Double_t weights[] = {0.5, 0.5};
    
    TH1 * hstat = ReweightSpectra(arrStat, weights, 0, "_2040");
    TH1 * hsyst = ReweightSpectra(arrSyst, weights, 1, "_2040");

    hyieldmean = PrintYieldAndError(hsyst, hstat, 0, mass[1]);

    
  }
  if(what == "kaons_neg_2040"){ 

    TObjArray * arrStat = new TObjArray();
    arrStat->Add(hSpectraCentr_stat[3][kMyKaon][kMyNeg]);
    arrStat->Add(hSpectraCentr_stat[4][kMyKaon][kMyNeg]);
    
    TObjArray * arrSyst = new TObjArray();
    arrSyst->Add(hSpectraCentr_sys [3][kMyKaon][kMyNeg]);
    arrSyst->Add(hSpectraCentr_sys [4][kMyKaon][kMyNeg]);
    
    Double_t weights[] = {0.5, 0.5};
    
    TH1 * hstat = ReweightSpectra(arrStat, weights, 0, "_2040");
    TH1 * hsyst = ReweightSpectra(arrSyst, weights, 1, "_2040");

    hyieldmean = PrintYieldAndError(hsyst, hstat, 0, mass[1]);

    
  }

  if(what == "protons_pos_2040"){ 

    TObjArray * arrStat = new TObjArray();
    arrStat->Add(hSpectraCentr_stat[3][kMyProton][kMyPos]);
    arrStat->Add(hSpectraCentr_stat[4][kMyProton][kMyPos]);
    
    TObjArray * arrSyst = new TObjArray();
    arrSyst->Add(hSpectraCentr_sys [3][kMyProton][kMyPos]);
    arrSyst->Add(hSpectraCentr_sys [4][kMyProton][kMyPos]);
    
    Double_t weights[] = {0.5, 0.5};
    
    TH1 * hstat = ReweightSpectra(arrStat, weights, 0, "_2040");
    TH1 * hsyst = ReweightSpectra(arrSyst, weights, 1, "_2040");

    hyieldmean = PrintYieldAndError(hsyst, hstat, 0, mass[2]);

    
  }
  if(what == "protons_neg_2040"){ 

    TObjArray * arrStat = new TObjArray();
    arrStat->Add(hSpectraCentr_stat[3][kMyProton][kMyNeg]);
    arrStat->Add(hSpectraCentr_stat[4][kMyProton][kMyNeg]);
    
    TObjArray * arrSyst = new TObjArray();
    arrSyst->Add(hSpectraCentr_sys [3][kMyProton][kMyNeg]);
    arrSyst->Add(hSpectraCentr_sys [4][kMyProton][kMyNeg]);
    
    Double_t weights[] = {0.5, 0.5};
    
    TH1 * hstat = ReweightSpectra(arrStat, weights, 0, "_2040");
    TH1 * hsyst = ReweightSpectra(arrSyst, weights, 1, "_2040");

    hyieldmean = PrintYieldAndError(hsyst, hstat, 0, mass[2]);

    
  }


  if(what == "lambda_0010"){ 

    TObjArray * arrStat = new TObjArray();
    arrStat->Add(hLambdaStat[0]);
    arrStat->Add(hLambdaStat[1]);
    
    TObjArray * arrSyst = new TObjArray();
    arrSyst->Add(hLambdaSyst[0]);
    arrSyst->Add(hLambdaSyst[1]);
    
    Double_t weights[] = {0.5, 0.5};
    
    TH1 * hstat = ReweightSpectra(arrStat, weights, 0, "_0010");
    TH1 * hsyst = ReweightSpectra(arrSyst, weights, 1, "_0010");

    hyieldmean = PrintYieldAndError(hsyst, hstat, 0, TDatabasePDG::Instance()->GetParticle("Lambda0")->Mass());

    
  }

  if(what == "k0s_0010"){ 

    TObjArray * arrStat = new TObjArray();
    arrStat->Add(hK0SStat[0]);
    arrStat->Add(hK0SStat[1]);
    
    TObjArray * arrSyst = new TObjArray();
    arrSyst->Add(hK0SSyst[0]);
    arrSyst->Add(hK0SSyst[1]);
    
    Double_t weights[] = {0.5, 0.5};
    
    TH1 * hstat = ReweightSpectra(arrStat, weights, 0, "_0010");
    TH1 * hsyst = ReweightSpectra(arrSyst, weights, 1, "_0010");

    hyieldmean = PrintYieldAndError(hsyst, hstat, 0, TDatabasePDG::Instance()->GetParticle("K_S0")->Mass());

    
  }



  
  TCanvas * c1 = new TCanvas("Averaging", "Averaging");
  c1->Divide(2,1);
  c1->cd(1);
  hstat->Draw();
  hsyst->Draw("same,e3");
  c1->cd(2);
  hyieldmean->Draw();
  

}


