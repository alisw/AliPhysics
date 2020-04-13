#include "src/Common.h"
#include "src/Utils.h"
using namespace utils;
#include "src/Plotting.h"
using namespace plotting;

#include <TFile.h>
#include <TDirectory.h>
#include <TLegend.h>
#include <TList.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TFitResult.h>
#include <TMatrixD.h>

#include <algorithm>

const char* kMultLab[kCentLength] = {"0-1","1-5","5-10","10-20","20-30","30-40","40-50","50-70","70-100","0-100"};
// const int kCentLength = 6; //Manuel
// const char* kMultLab[kCentLength] = {"0-5","5-10","10-20","20-40","40-100","0-100"};
const int kNspecies = 9;
const char* kSpeciesName[kNspecies] = {"Pi","Ka","Pr","Phi","Ks","K0","Lambda","Xi","Om"};
const char* kSpeciesLabels[kNspecies] = {"#pi", "k", "p", "#phi", "K_{s}", "K_{0}", "#Lambda", "#Xi", "#Omega"};
const float kSpeciesMasses[kNspecies] = {0.140, 0.494, 0.936};
const char* kMatter[2] = {"Pos","Neg"};
const char kSign[2] = {'+','-'};

void SignalLoss(int n_analysed_species=3, bool save_all_plots = false){

  TFile input(Form("%sSignalLoss.root",kBaseOutputDir.data()));
  TFile output(Form("%sSignalLoss_Correction.root",kBaseOutputDir.data()),"recreate");
  gStyle->SetOptStat(0);

  TCanvas* cSignalLoss[kCentLength][2];
  TCanvas* cTotSignalLoss[kCentLength];
  TCanvas* cRatioToPion[kNspecies][2];
  TCanvas* cTotRatioToPion[kNspecies];
  TH1F* hMeanSignalLoss[kCentLength][2];
  TH1F* hTotMeanSignalLoss[kCentLength];
  TH1F* hTotSystSignalLoss[kCentLength];
  TCanvas* cTotMeanSignalLoss = new TCanvas("cTotMeanSignalLoss","cTotMeanSignalLoss");
  TCanvas* cTotSpecies[kNspecies];
  TH1F* hSp[kCentLength][kNspecies][2];
  TH1F* hSpTot[kCentLength][kNspecies];
  TH1F* hRatioToPion[kCentLength][kNspecies][2];
  TH1F* hTotRatioToPion[kCentLength][kNspecies];
  TH1F* hDeutSignalLoss[kCentLength][2];
  TH1F* hRatioMatterAntimatter[kCentLength][3];
  TCanvas* cRatioMatterAntimatter[3];

  for(int iSpecies=0; iSpecies<n_analysed_species; iSpecies++){
    //Ratio to pion histograms
    cTotSpecies[iSpecies] = new TCanvas(Form("cTot_%s",kSpeciesName[iSpecies]),Form("cTot_%s",kSpeciesName[iSpecies]));
    cTotRatioToPion[iSpecies] = new TCanvas(Form("cTotRtP_%s",kSpeciesName[iSpecies]),Form("cTotRtP_%s",kSpeciesName[iSpecies]),3200,1800);
    for(int iS=0; iS<2; iS++){
      cRatioToPion[iSpecies][iS] = new TCanvas(Form("cRtP_%s%s",kSpeciesName[iSpecies],kMatter[iS]),Form("cRtP_%s%s",kSpeciesName[iSpecies],kMatter[iS]),3200,1800);
    }
  }

  TDirectory *ratio_dir = output.mkdir("ratios");
  TLegend* cent_leg = new TLegend(0.74,0.19,0.9,0.46);
  cent_leg->SetBorderSize(0);
  for(int iC=0; iC<kCentLength; iC++){

    cTotSignalLoss[iC] = new TCanvas(Form("cTotSignalLoss_%d",iC),Form("cTotSignalLoss_%d",iC));
    TDirectory *c_dir = output.mkdir(Form("%d",iC));

    for(int iS=0; iS<2; iS++){

      //deuteron signal loss
      hDeutSignalLoss[iC][iS] = new TH1F(Form("hDeutSignalLoss_%d_%s",iC,kMatter[iS]),Form("%s%%;#it{p}_{T} (GeV/#it{c}); #epsilon_{part}",kMultLab[iC]),kNPtBins,kPtBins);
      SetHistStyle(hDeutSignalLoss[iC][iS],kSpectraColors[iC],45);

      //Rebinned input histograms
      cSignalLoss[iC][iS] = new TCanvas(Form("cSL_%d_%s",iC,kMatter[iS]),Form("cSL_%d_%s",iC,kMatter[iS]),3200,1800);
      TLegend* leg = new TLegend(0.80,0.25,0.85,0.45);
      leg->SetBorderSize(0);
      //Mean correction
      hMeanSignalLoss[iC][iS] = new TH1F(Form("hMeanSignalLoss_%d_%s",iC,kMatter[iS]),Form("%s%%;#it{p}_{T} (GeV/#it{c});1/#epsilon_{part}",kMultLab[iC]),kNPtBins,kPtBins);
      hMeanSignalLoss[iC][iS]->GetYaxis()->SetRangeUser(0.7,1.3);

      hTotMeanSignalLoss[iC] = new TH1F(Form("hTotMeanSignalLoss_%d",iC),Form("%s%%;#it{p}_{T} (GeV/#it{c});1/#epsilon_{part}",kMultLab[iC]),kNPtBins,kPtBins);
      hTotMeanSignalLoss[iC]->GetYaxis()->SetRangeUser(0.7,1.3);

      hTotSystSignalLoss[iC] = new TH1F(Form("hTotSystSignalLoss_%d",iC),Form("%s%%;#it{p}_{T} (GeV/#it{c});Systematic uncertainty",kMultLab[iC]),kNPtBins,kPtBins);
      hTotSystSignalLoss[iC]->GetYaxis()->SetRangeUser(0.,1.);

      for(int iSpecies=0; iSpecies<n_analysed_species; iSpecies++){
        TList* lSp = (TList*)input.Get(Form("%s%s",kSpeciesName[iSpecies],kMatter[iS]));
        Requires(lSp,Form("list %s%s",kSpeciesName[iSpecies],kMatter[iS]));
        TH1F *hSp_tmp = (TH1F*)lSp->FindObject(Form("sgnLoss_%s%s_%s",kSpeciesName[iSpecies],kMatter[iS],kMultLab[iC]));
        Requires(hSp_tmp,Form("Missing %s",Form("histogram sgnLoss_%s%s_%s",kSpeciesName[iSpecies],kMatter[iS],kMultLab[iC])));
        hSp[iC][iSpecies][iS] = new TH1F(Form("h%s%s_%d",kSpeciesName[iSpecies],kMatter[iS],iC),Form("%s%%;#it{p}_{T} (GeV/#it{c}); 1/#epsilon_{part}",kMultLab[iC]),kNPtBins,kPtBins);
        //
        hRatioToPion[iC][iSpecies][iS] = new TH1F(Form("hRatioToPion_%s%s_%s",kSpeciesName[iSpecies],kMatter[iS],kMultLab[iC]),Form("; #it{p}_{T} (GeV/#it{c});%s^{%c} / %s^{%c}",kSpeciesLabels[iSpecies],kSign[iS],kSpeciesLabels[0],kSign[iS]),kNPtBins,kPtBins);
        hRatioToPion[iC][iSpecies][iS]->GetYaxis()->SetRangeUser(-2,2);

        int iSpBin = 1;
        float content_tmp[kNPtBins] = {0.};
        float err2_tmp[kNPtBins] = {0.};
        int count_tmp[kNPtBins] = {0};
        //Using the bins of the analysis
        for(int iB=1; iB<=hSp_tmp->GetNbinsX(); iB++){
          if(hSp_tmp->GetBinCenter(iB) < kPtBins[0]) continue;
          if(hSp_tmp->GetBinCenter(iB) > kCentPtLimits[iC]) break;
          if(hSp_tmp->GetBinCenter(iB) < kPtBins[iSpBin]){
            content_tmp[iSpBin-1] += hSp_tmp->GetBinContent(iB);
            err2_tmp[iSpBin-1] += Sq(hSp_tmp->GetBinError(iB));
            count_tmp[iSpBin-1]++;
          }
          else{
            iSpBin++;
            content_tmp[iSpBin-1] += hSp_tmp->GetBinContent(iB);
            err2_tmp[iSpBin-1] += Sq(hSp_tmp->GetBinError(iB));
            count_tmp[iSpBin-1]++;
          }
        }
        for(int i=0; i<kNPtBins; i++){
          if(hSp[iC][iSpecies][iS]->GetBinCenter(i+1)>kCentPtLimits[iC]) continue;
          hSp[iC][iSpecies][iS]->SetBinContent(i+1,content_tmp[i]/count_tmp[i]);
          hSp[iC][iSpecies][iS]->SetBinError(i+1,TMath::Sqrt(err2_tmp[i])/count_tmp[i]);
        }
        hSp[iC][iSpecies][iS]->SetDirectory(0);
        hSp[iC][iSpecies][iS]->GetYaxis()->SetRangeUser(0.7,1.3);
        SetHistStyle(hSp[iC][iSpecies][iS],kSpectraColors[iC],20+iSpecies);
        // creating unified histogram
        if(iS==1){
          hSpTot[iC][iSpecies] = (TH1F*)hSp[iC][iSpecies][0]->Clone(Form("hTot_%s_%d",kSpeciesName[iSpecies],iC));
          hSpTot[iC][iSpecies]->Add(hSp[iC][iSpecies][1]);
          hSpTot[iC][iSpecies]->Scale(0.5);
          hSpTot[iC][iSpecies]->GetYaxis()->SetRangeUser(0.7,1.3);
          hTotRatioToPion[iC][iSpecies]= (TH1F*) hSpTot[iC][iSpecies]->Clone(Form("hTotRatioToPion_%s_%s",kSpeciesName[iSpecies],kMultLab[iC]));
          hTotRatioToPion[iC][iSpecies]->Reset();
          hTotRatioToPion[iC][iSpecies]->SetTitle(Form("; #it{p}_{T} (GeV/#it{c});%s / %s",kSpeciesLabels[iSpecies],kSpeciesLabels[0]));
          hTotRatioToPion[iC][iSpecies]->GetYaxis()->SetRangeUser(-2,2);
          //hTotRatioToPion[iC][iSpecies] = new TH1F(Form("hTotRatioToPion_%s_%s",kSpeciesName[iSpecies],kMultLab[iC]),Form("; #it{p}_{T} (GeV/#it{c});%s / %s",kSpeciesLabels[iSpecies],kSpeciesLabels[0]),n_pt_bins,pt_bin_limits);
        }
        //Ratio antimatter-matter
        if(iS==1 && iSpecies<n_analysed_species){
          hRatioMatterAntimatter[iC][iSpecies] = (TH1F*)hSp[iC][iSpecies][1]->Clone(Form("hRatioAM_%s_%d",kSpeciesName[iSpecies],iC));
          hRatioMatterAntimatter[iC][iSpecies]->GetYaxis()->SetTitle(Form("%s^{%c}/%s^{%c}",kSpeciesLabels[iSpecies],kSign[1],kSpeciesLabels[iSpecies],kSign[0]));
          hRatioMatterAntimatter[iC][iSpecies]->Divide(hSp[iC][iSpecies][0]);
          hRatioMatterAntimatter[iC][iSpecies]->GetYaxis()->SetRangeUser(0.95,1.05);
          if(iC==0) cRatioMatterAntimatter[iSpecies] = new TCanvas(Form("cRatioAM_%s",kSpeciesName[iSpecies]),Form("cRatioAM_%s",kSpeciesName[iSpecies]),3200,2400);
        }
        leg->AddEntry(hSp[iC][iSpecies][iS],Form("%s^{%c}",kSpeciesLabels[iSpecies],kSign[iS]),"PE");
        cSignalLoss[iC][iS]->cd();
        if(iSpecies==0) hSp[iC][iSpecies][iS]->Draw();
        else hSp[iC][iSpecies][iS]->Draw("same");
      }
      //leg->AddEntry(hDeutSignalLoss[iC][iS],"d","PE");
      leg->Draw();
      //plotting pions
      c_dir->cd();
      if(save_all_plots) hSp[iC][0][iS]->Write();
      if(iS==1){
        hSpTot[iC][0]->Write();
        hRatioMatterAntimatter[iC][0]->Write();
      }
      //Plotting the ratio to pions of the correction
      for(int iSpecies=1; iSpecies<n_analysed_species; iSpecies++){
        hRatioToPion[iC][iSpecies][iS]->Divide(hSp[iC][iSpecies][iS],hSp[iC][0][iS]);
        SetHistStyle(hRatioToPion[iC][iSpecies][iS],kSpectraColors[iC]);
        hRatioToPion[iC][iSpecies][iS]->GetYaxis()->SetRangeUser(0.7,1.3);
        c_dir->cd();
        if(save_all_plots) hRatioToPion[iC][iSpecies][iS]->Write();
        if(iS==1){
          hTotRatioToPion[iC][iSpecies]->Divide(hSpTot[iC][iSpecies],hSpTot[iC][0]);
          SetHistStyle(hTotRatioToPion[iC][iSpecies],kSpectraColors[iC]);
          hTotRatioToPion[iC][iSpecies]->GetYaxis()->SetRangeUser(0.7,1.3);
          c_dir->cd();
          hTotRatioToPion[iC][iSpecies]->Write();
        }
        //SetHistStyle(hSp[iC][iSpecies][iS],kSpectraColors[iC]);
        if(save_all_plots) hSp[iC][iSpecies][iS]->Write();
        if(iS==1){
          hSpTot[iC][iSpecies]->Write();
          if(iSpecies<n_analysed_species) hRatioMatterAntimatter[iC][iSpecies]->Write();
        }
      }
      //Computing mean value of the correction
      for(int iB=1; iB<=kNPtBins; iB++){
        if(hSp[iC][0][iS]->GetBinCenter(iB)>kCentPtLimits[iC]) continue;
        vector<float> values;
        vector<float> errors;
        for(int iSpecies=0; iSpecies<n_analysed_species; iSpecies++){
          values.push_back(hSp[iC][iSpecies][iS]->GetBinContent(iB));
          errors.push_back(hSp[iC][iSpecies][iS]->GetBinError(iB));
        }
        hMeanSignalLoss[iC][iS]->SetBinContent(iB,TMath::Mean(values.begin(),values.end()));
        float error_val = 0.;
        for(auto p : errors){
          error_val += Sq(p);
        }
        error_val = TMath::Sqrt(error_val)/errors.size();
        hMeanSignalLoss[iC][iS]->SetBinError(iB,error_val);
      }
      c_dir->cd();
      SetHistStyle(hMeanSignalLoss[iC][iS],kSpectraColors[iC]);
      hMeanSignalLoss[iC][iS]->Write();
      //Computing mean value of the correction for both matter and anti-matter
      if(iS==1){
        for(int iB=1; iB<=kNPtBins; iB++){
          if(hSpTot[iC][0]->GetBinCenter(iB)>kCentPtLimits[iC]) continue;
          vector<float> values;
          vector<float> errors;
          for(int iSpecies=0; iSpecies<n_analysed_species; iSpecies++){
            values.push_back(hSpTot[iC][iSpecies]->GetBinContent(iB));
            errors.push_back(hSpTot[iC][iSpecies]->GetBinError(iB));
          }
          float val = TMath::Mean(values.begin(),values.end());
          hTotMeanSignalLoss[iC]->SetBinContent(iB,val);
          float error_val = 0.;
          error_val = (*std::max_element(values.begin(),values.end())-*std::min_element(values.begin(),values.end()))/2;
          hTotMeanSignalLoss[iC]->SetBinError(iB,error_val);
          hTotSystSignalLoss[iC]->SetBinContent(iB,error_val/val);
          hTotSystSignalLoss[iC]->SetBinError(iB,0);
        }
        c_dir->cd();
        SetHistStyle(hTotMeanSignalLoss[iC],kSpectraColors[iC]);
        hTotMeanSignalLoss[iC]->Write();
        hTotSystSignalLoss[iC]->Write();
        cTotMeanSignalLoss->cd();
        cent_leg->AddEntry(hTotMeanSignalLoss[iC],Form("%s %%",kMultLab[iC]),"PL");
        if(!iC) hTotMeanSignalLoss[iC]->Draw("");
        else  hTotMeanSignalLoss[iC]->Draw("SAME");
      }
    }
  }
  cTotMeanSignalLoss->cd();
  cent_leg->Draw();
  output.cd();
  cTotMeanSignalLoss->Write();

  for(int iSpecies=0; iSpecies<n_analysed_species; iSpecies++){
    TLegend* leg = new TLegend(0.74,0.19,0.9,0.46);
    leg->SetBorderSize(0);
    leg->SetHeader(Form("%s",kSpeciesLabels[iSpecies]));
    cTotSpecies[iSpecies]->cd();
    for(int iC=0; iC<kCentLength; iC++){
      if (iC==0) hSpTot[iC][iSpecies]->Draw();
      else hSpTot[iC][iSpecies]->Draw("same");
      leg->AddEntry(hSpTot[iC][iSpecies],Form("%s %%",kMultLab[iC]),"PE");
    }
    leg->Draw();
    output.cd();
    cTotSpecies[iSpecies]->Write();
  }

  for(int iC=0; iC<kCentLength; iC++){
    cTotSignalLoss[iC]->cd();
    TLegend leg(.74,0.19,0.9,0.46);
    leg.SetBorderSize(0);
    for(int iSpecies=0; iSpecies<3; iSpecies++){
      if(!iC) hSpTot[iC][iSpecies]->Draw();
      else  hSpTot[iC][iSpecies]->Draw("same");
      leg.AddEntry(hSpTot[iC][iSpecies],Form("%s",kSpeciesLabels[iSpecies]),"PR");
    }
    leg.Draw();
    cTotSignalLoss[iC]->Write();
  }

  //Plotting ratio

  ratio_dir->cd();
  for(int iS=0; iS<2; iS++){
    for(int iSpecies=1; iSpecies<n_analysed_species; iSpecies++){
      TLegend* leg_rtp = new TLegend(0.74,0.19,0.9,0.46);
      leg_rtp->SetBorderSize(0);
      for(int iC=0; iC<kCentLength; iC++){
        cRatioToPion[iSpecies][iS]->cd();
        hRatioToPion[iC][iSpecies][iS]->GetYaxis()->SetRangeUser(0.95,1.05);
        if(!iC) hRatioToPion[iC][iSpecies][iS]->Draw();
        else hRatioToPion[iC][iSpecies][iS]->Draw("same");
        leg_rtp->AddEntry(hRatioToPion[iC][iSpecies][iS],Form("%s %%",kMultLab[iC]),"PE");
        if(iS==1){
          cTotRatioToPion[iSpecies]->cd();
          hTotRatioToPion[iC][iSpecies]->GetYaxis()->SetRangeUser(0.95,1.05);
          if(!iC) hTotRatioToPion[iC][iSpecies]->Draw();
          else hTotRatioToPion[iC][iSpecies]->Draw("same");
          if(iSpecies<3){
            cRatioMatterAntimatter[iSpecies]->cd();
            if(!iC) hRatioMatterAntimatter[iC][iSpecies]->Draw();
            else hRatioMatterAntimatter[iC][iSpecies]->Draw("same");
            if(iSpecies==1){
              cRatioMatterAntimatter[0]->cd();
              if(!iC) hRatioMatterAntimatter[iC][0]->Draw();
              else hRatioMatterAntimatter[iC][0]->Draw("same");
            }
          }
        }
      }
      cRatioToPion[iSpecies][iS]->cd();
      leg_rtp->Draw();
      if(save_all_plots) cRatioToPion[iSpecies][iS]->Write();
      // tot ratio
      if(iS==1){
        cTotRatioToPion[iSpecies]->cd();
        leg_rtp->Draw();
        cTotRatioToPion[iSpecies]->Write();
        if(iSpecies<n_analysed_species){
          cRatioMatterAntimatter[iSpecies]->cd();
          leg_rtp->Draw();
          cRatioMatterAntimatter[iSpecies]->Write();
          if(iSpecies==1){
            cRatioMatterAntimatter[0]->cd();
            leg_rtp->Draw();
            cRatioMatterAntimatter[0]->Write();
          }
        }
      }
    }
  }

  //Mass scaling

  TH1F* Scaling[kNPtBins];
  TCanvas* cMassScaling[kCentLength][2];
  TH1F* hMassScaling[kCentLength][kNPtBins][2];
  for(int iC=0; iC<kCentLength; iC++){
    for(int iS=0; iS<2; iS++){
      cMassScaling[iC][iS] = new TCanvas(Form("cMassScaling_%d_%s",iC,kMatter[iS]),Form("cMassScaling_%d_%s",iC,kMatter[iS]));
      cMassScaling[iC][iS]->Divide(5,3);
      for(int iB=1; iB<=kNPtBins; iB++){
        if(hSp[iC][0][iS]->GetBinCenter(iB)>kCentPtLimits[iC]) break;
        cMassScaling[iC][iS]->cd(iB);
        hMassScaling[iC][iB][iS] = new TH1F(Form("hMassScaling_%d_%d_%c",iC,iB,kLetter[iS]),Form("%.1f < #it{p}_{T} #leq %.1f (GeV/#it{c});m (GeV/#it{c}^{2});1/#epsilon_{part}",kPtBins[iB-1],kPtBins[iB]),101,-0.05,1.05);
        for(int iSpecies=0; iSpecies<n_analysed_species; iSpecies++){
          hMassScaling[iC][iB][iS]->SetBinContent(hMassScaling[iC][iB][iS]->FindBin(kSpeciesMasses[iSpecies]),hSp[iC][iSpecies][iS]->GetBinContent(iB));
          hMassScaling[iC][iB][iS]->SetBinError(hMassScaling[iC][iB][iS]->FindBin(kSpeciesMasses[iSpecies]),hSp[iC][iSpecies][iS]->GetBinError(iB));
        }
        TF1 fitfunc("fitfunc","pol1");
        TFitResultPtr result = hMassScaling[iC][iB][iS]->Fit("fitfunc","SQ");
        TMatrixD cov_mat = result->GetCovarianceMatrix();
        float cov = cov_mat[0][1];

        hMassScaling[iC][iB][iS]->GetYaxis()->SetRangeUser(0.9*hMassScaling[iC][iB][iS]->GetBinContent(hMassScaling[iC][iB][iS]->FindBin(kSpeciesMasses[2])), 1.1*hMassScaling[iC][iB][iS]->GetMaximum());
        SetHistStyle(hMassScaling[iC][iB][iS],kBlack);
        hMassScaling[iC][iB][iS]->Draw("PE");
        output.cd(Form("%d",iC));
        if(save_all_plots) hMassScaling[iC][iB][iS]->Write();
        hDeutSignalLoss[iC][iS]->SetBinContent(iB,fitfunc.Eval(1.876)); // it is the deuteron mass
        float errore = TMath::Sqrt(Sq(1.876*fitfunc.GetParError(1))+Sq(fitfunc.GetParError(1))+2*cov*1.876);
        hDeutSignalLoss[iC][iS]->SetBinError(iB,errore);
      }
      output.cd(Form("%d",iC));
      cSignalLoss[iC][iS]->cd();
      hDeutSignalLoss[iC][iS]->Draw("same");
      output.cd(Form("%d",iC));
      cSignalLoss[iC][iS]->Write();
      hDeutSignalLoss[iC][iS]->Write();
      cMassScaling[iC][iS]->Write();
    }
  }
}
