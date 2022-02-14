# include "TH1F.h"
# include <iostream>
# include "TAxis.h"
# include "TCanvas.h"
# include "TLine.h"
# include "TLatex.h"
# include "TProfile.h"
# include "TProfile2D.h"
# include "TGraph.h"
# include "TF1.h"
# include "TList.h"
# include "TMath.h"
# include "TFile.h"
# include "TApplication.h"
# include "TStyle.h"
# include "TFractionFitter.h"
# include "TText.h"

int fits_pT_bins(){
  const char* names[6] = {"p", "#bar{p}", "#pi^{+}", "#pi^{-}", "K^{+}", "K^{-}"};
  const char* part_pairs[6] = {"PP","aPaP","PIpPIp","PImPIm", "KpKp", "KmKm"};
  const char* part_pairs_mc[6] = {"Proton", "ProtonMinus", "Pion", "PionMinus", "Kaon", "KaonMinus"};
  int pair_nb = 3;

  double pT_min = 0.5;
  double pT_max = 2.5;

  ofstream myfile;
  myfile.open(Form("output_fits_pT_%.1f_%.1f.txt", pT_min, pT_max));

  for(int iii = 0; iii < 4; iii++){

    pair_nb = iii;
  
    TH2D *mc_prim2D_XY, *mc_secMat2D_XY, *mc_secWeak2D_XY, *mc_prim2D_Z, *mc_secMat2D_Z, *mc_secWeak2D_Z;
    TH1D *mc_prim_XY, *mc_secMat_XY, *mc_secWeak_XY, *mc_prim_Z, *mc_secMat_Z, *mc_secWeak_Z;
    TH2D *data2D_XY, *data2D_Z;
    TH1D *data_XY, *data_Z, *data_XY_scaled;
    TH1D *hist1D_sm_prim_XY, *hist1D_sm_secMat_XY, *hist1D_sm_secWeak_XY, *hist1D_sm_prim_Z, *hist1D_sm_secMat_Z, *hist1D_sm_secWeak_Z;

    double fraction1, fraction2, fraction3;
  
    TFile *ifile_data = new TFile("./real/Mergeds.dir.root");
    data2D_XY = (TH2D*)ifile_data->Get(Form("DCARPtcutPass%stpcM0", part_pairs[pair_nb]));
    data2D_XY->GetYaxis()->SetRangeUser(pT_min, pT_max);
    data_XY = data2D_XY->ProjectionX();
    data_XY->Sumw2();
  
    // mirror bin - in case histos are not filled in the whole range
    int bin = 0;
    int mirror_bin = 2*data_XY->FindBin(0);
    for(int oo = 0; oo < data_XY->GetNbinsX()/2; oo++){
      data_XY->SetBinContent(bin+1, data_XY->GetBinContent(mirror_bin-2));
      data_XY->SetBinError(bin+1, data_XY->GetBinError(mirror_bin-2));
      bin++;
      mirror_bin--;
    }

    // reading MC data
    TFile *ifile_mc = new TFile("./MC/Mergeds.dir.root");
    mc_prim2D_XY = (TH2D*)ifile_mc->Get(Form("hPrimVsDCAXYM0%s", part_pairs_mc[pair_nb]))->Clone();
    mc_secMat2D_XY = (TH2D*)ifile_mc->Get(Form("hSecMatVsDCAXYM0%s", part_pairs_mc[pair_nb]))->Clone();
    mc_secWeak2D_XY = (TH2D*)ifile_mc->Get(Form("hSecWeakVsDCAXYM0%s", part_pairs_mc[pair_nb]))->Clone();


    mc_prim2D_XY->GetYaxis()->SetRangeUser(pT_min, pT_max);
    mc_secMat2D_XY->GetYaxis()->SetRangeUser(pT_min, pT_max);
    mc_secWeak2D_XY->GetYaxis()->SetRangeUser(pT_min, pT_max);

    mc_prim_XY = mc_prim2D_XY->ProjectionX();
    mc_secMat_XY = mc_secMat2D_XY->ProjectionX();
    mc_secWeak_XY = mc_secWeak2D_XY->ProjectionX();

    TObjArray *mc = new TObjArray(3);

    double scale_mc = 0;
    double scale_data = data_XY->Integral();

    hist1D_sm_prim_XY = new TH1D("histo_prim", "histo_prim", 400, -3, 3);
    hist1D_sm_secMat_XY = new TH1D("histo_secMat", "histo_secMat", 400, -3, 3);
    hist1D_sm_secWeak_XY = new TH1D("histo_secWeak", "histo_secWeak", 400, -3, 3);

    hist1D_sm_prim_XY = (TH1D*) mc_prim_XY->Clone();
    hist1D_sm_secMat_XY = (TH1D*) mc_secMat_XY->Clone();
    hist1D_sm_secWeak_XY = (TH1D*) mc_secWeak_XY->Clone();

    // mirrors again
    bin = hist1D_sm_prim_XY->GetNbinsX();
    mirror_bin = 1;
    for(int oo = 0; oo < hist1D_sm_prim_XY->GetNbinsX()/2; oo++){
      hist1D_sm_prim_XY->SetBinContent(bin-1, hist1D_sm_prim_XY->GetBinContent(mirror_bin));
      hist1D_sm_prim_XY->SetBinError(bin-1, hist1D_sm_prim_XY->GetBinError(mirror_bin));
      hist1D_sm_secMat_XY->SetBinContent(bin-1, hist1D_sm_secMat_XY->GetBinContent(mirror_bin));
      hist1D_sm_secMat_XY->SetBinError(bin-1, hist1D_sm_secMat_XY->GetBinError(mirror_bin));
      hist1D_sm_secWeak_XY->SetBinContent(bin-1, hist1D_sm_secWeak_XY->GetBinContent(mirror_bin));
      hist1D_sm_secWeak_XY->SetBinError(bin-1, hist1D_sm_secWeak_XY->GetBinError(mirror_bin));
      bin--;
      mirror_bin++;
    }

    scale_mc += hist1D_sm_prim_XY->Integral();
    scale_mc += hist1D_sm_secMat_XY->Integral();
    scale_mc += hist1D_sm_secWeak_XY->Integral();

    data_XY->Scale(scale_mc/scale_data);

    mc->Add(hist1D_sm_prim_XY);
    mc->Add(hist1D_sm_secMat_XY);
    mc->Add(hist1D_sm_secWeak_XY);

    TFractionFitter *fit = new TFractionFitter(data_XY, mc);

    // fine tuning
    if(pair_nb == 0){ // protons
      fit->Constrain(0, 0.5, 0.8); // prim
      fit->Constrain(1, 0.002, 100.0); // material
      fit->Constrain(2, 0.02, 100.0); // weak
    }
    else if(pair_nb == 1){ // anti-protons
      fit->Constrain(0, 0.6, 0.89); // prim
      fit->Constrain(1, 0.002, 100.0); // material
      fit->Constrain(2, 0.02, 100.0); // weak
    }
    else if(pair_nb == 2){ // pions   
      fit->Constrain(0, 0.89, 0.99); // prim
      fit->Constrain(1, 0.005, 1.0); // material
      fit->Constrain(2, 0.03, 1.0); // weak
    }
    else if(pair_nb == 3){ //pions minus
      fit->Constrain(0, 0.89, 0.99); // prim
      fit->Constrain(1, 0.002, 1.0); // material
      fit->Constrain(2, 0.02, 1.0); // weak

    }

    Int_t status = fit->Fit();
    cout << "status " << status << endl;
    double frac0, frac1, frac2;
    double err0, err1, err2;
    fit->GetResult(0, frac0, err0);
    fit->GetResult(1, frac1, err1);
    fit->GetResult(2, frac2, err2);

    // cout << "frac0, frac1, frac2 " << frac0 << " " << frac1 << " " << frac2 << endl << "sum: " << frac0+frac1+frac2 << endl;
  
    // plots
    TCanvas* cnv = new TCanvas(Form("DCA_XY_%d", pair_nb),Form("DCA_XY_%d", pair_nb),800,568);
    TPad* cXY = new TPad(Form("GraphsXY_%d",pair_nb),Form("GraphsXY_%d",pair_nb),0.01,0.05,0.95,0.95);
    gStyle->SetOptStat(000);
    gPad->SetLogy();

    if(status == 0){
  
      TH1D* MC0=(TH1D*)fit->GetMCPrediction(0);
      MC0->SetLineColor(kRed);
      TH1D* MC1=(TH1D*)fit->GetMCPrediction(1);
      MC1->SetLineColor(kBlue);
      MC1->SetTitle("DCA_{XY} distribution ");
      TH1D* MC2=(TH1D*)fit->GetMCPrediction(2);
      MC2->SetLineColor(kGreen+2);
      MC2->SetTitle("DCA_{XY} distribution ");
      hist1D_sm_prim_XY->SetLineColor(kRed);
      hist1D_sm_secMat_XY->SetLineColor(kBlue);
      hist1D_sm_secWeak_XY->SetLineColor(kGreen+2);

      data_XY->SetMarkerColor(kBlack);
      data_XY->SetMarkerStyle(21);
      data_XY->SetMarkerSize(0.7);
      data_XY->GetXaxis()->SetTitle("DCA_{XY} (cm)");
      data_XY->GetYaxis()->SetTitle("Entries");  

      data_XY->GetYaxis()->SetRangeUser(1, 50000000);
    
      data_XY->Draw("p");
      data_XY->SetTitle(Form("%s DCA_{XY} distribution", names[iii]));
      TH1D* result = (TH1D*) fit->GetPlot();
      result->SetLineColor(kMagenta);
      result->SetLineWidth(2);
      result->SetTitle("DCA_{XY} distribution ");
      result->Draw("hist same");
    
      hist1D_sm_prim_XY->Scale(frac0);
      hist1D_sm_secMat_XY->Scale(frac1);
      hist1D_sm_secWeak_XY->Scale(frac2);

      auto legend = new TLegend(0.6,0.7,0.88,0.88);
      legend->SetHeader("#sqrt{s} = 13 TeV, pp","C"); // option "C" allows to center the header
      legend->AddEntry(data_XY, "Data", "p");
      legend->AddEntry(result, "Combined fit", "l");
      legend->AddEntry(MC0,"Primaries","l");
      legend->AddEntry(MC1, "Sec. from material","l");
      legend->AddEntry(MC2,"Sec. form weak decays","l");
      // legend->AddEntry(hist1D_sm_prim_XY,"Primaries","l");
      // legend->AddEntry(hist1D_sm_secMat_XY, "Sec. from material","l");
      // legend->AddEntry(hist1D_sm_secWeak_XY,"Sec. form weak decays","l");
      legend->SetBorderSize(0);
      legend->Draw();

      hist1D_sm_prim_XY->SetTitle("DCA_{XY} distribution ");
      hist1D_sm_secMat_XY->SetTitle("DCA_{XY} distribution ");
      hist1D_sm_secWeak_XY->SetTitle("DCA_{XY} distribution ");
    
      // from template
      MC0->Draw("HISTsame");
      MC1->Draw("HISTsame");
      MC2->Draw("HISTsame");

      myfile << "********************************" << endl << part_pairs_mc[iii] << endl << "********************************" <<  endl;
      myfile << "Scale_mc: " << scale_mc << " Sum: " << MC0->Integral()+MC1->Integral()+MC2->Integral() << " data int.: " << data_XY->Integral() << endl;
      myfile << "Diff.: " << scale_mc-(MC0->Integral()+MC1->Integral()+MC2->Integral()) << endl;
      myfile << pair_nb << " particle DCA_XY"  << endl;
      double sum_MC = MC0->Integral()+MC1->Integral()+MC2->Integral();
      myfile << "Primaries: " << MC0->Integral()/sum_MC << endl;
      myfile << "Sec. mat.: " << MC1->Integral()/sum_MC << endl;
      myfile << "Sec. weak: " << MC2->Integral()/sum_MC << endl;

    }
  }
  
  return 0;
}
