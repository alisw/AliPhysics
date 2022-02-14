#include "src/Common.h"
#include "src/Utils.h"
using namespace utils;
#include <map>
#include <string>
#include <sstream>
#include <vector>

#include <TStyle.h>
#include <TFile.h>
#include <TList.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TF1.h>

using std::map;
using std::string;
using std::vector;

const char *base_dir = "/Users/lbariogl/cernbox/Deuterons13TeV/signal_loss/20180401_1009_LHC17d20a2_extra";

const char *prtName_ltx[] = { "#pi", "K", "p"};

const char *prtName_str[] = { "Pi", "Ka",  "Pr"};

const char *suffix[2] = {"Pos","Neg"};

const int kNCentBins = kCentLength;
const char* kMultLab[kNCentBins] ={"0-1","1-5","5-10","10-20","20-30","30-40","40-50","50-70","70-100","0-100"};
// const int kNCentBins = 6; //Manuel
// const char* kMultLab[kNCentBins] = {"0-5","5-10","10-20","20-40","40-100","0-100"};

const int kNspc = 3;

const unsigned int lowCentBin[kNCentBins] = { 2, 3, 4, 5, 6, 7, 8, 9, 11, 2 };
const unsigned int upCentBin[kNCentBins] = { 2, 3, 4, 5, 6, 7, 8, 10, 13, 13 };

// const unsigned int lowCentBin[kNCentBins] = {2,4,5,6,8,2};
// const unsigned int upCentBin[kNCentBins] = {3,4,5,7,13,13};


void GetSignalLossCorrection(){
  gStyle->SetOptStat(0);

  TFile input_file(Form("%s/AnalysisResults.root",base_dir),"read");
  if(!input_file.IsOpen()){
    printf("File not found\n");
    exit(1);
  }
  //TFile outoput_file(kSignalLossInput.data(),"recreate");
  TFile output_file("../results/ProvaSignalLoss.root","recreate");

  TList *list = (TList*)input_file.Get("signal_loss");
  if(!list){
    printf("Missing list: signal_loss\n");
    exit(1);
  }

  TH1F * fHistAccEvents = (TH1F*)list->FindObject("fHistAccEvents");
  Requires(fHistAccEvents,"fHistAccEvents");
  fHistAccEvents->SetDirectory(0);
  TH1F * fHistTrueINELgt0Events = (TH1F*)list->FindObject("fHistTrueINELgt0Events");
  Requires(fHistTrueINELgt0Events,"fHistTrueINELgt0Events");
  fHistTrueINELgt0Events->SetDirectory(0);

  output_file.cd();
  fHistAccEvents->Write();
  fHistTrueINELgt0Events->Write();
  TH1F* fHistEventLoss = (TH1F*) fHistAccEvents->Clone("fHistEventLoss");
  fHistEventLoss->Reset();
  fHistEventLoss->GetYaxis()->SetTitle("#epsilon_{evt}");
  fHistEventLoss->Divide(fHistAccEvents,fHistTrueINELgt0Events,1,1,"B");
  output_file.cd();
  fHistEventLoss->Write();

  int NormAcc[kCentLength];
  int NormTrue[kCentLength];
  float ratioNorm[kCentLength];
  for(int iC=0; iC<kCentLength; iC++){
    NormAcc[iC]=fHistAccEvents->Integral(lowCentBin[iC],upCentBin[iC]);
    NormTrue[iC]=fHistTrueINELgt0Events->Integral(lowCentBin[iC],upCentBin[iC]);
    ratioNorm[iC] = (float)NormAcc[iC]/(float)NormTrue[iC];
    printf("iC: %d NormAcc: %d NormTrue: %d ratioNorm: %f\n",iC,NormAcc[iC],NormTrue[iC],ratioNorm[iC]);
  }

  TH3F* fHistGenPartsAcc_tmp[2];
  TH3F* fHistGenPartsINELgt0_tmp[2];
  for(int iS=0; iS<2; iS++){
    fHistGenPartsAcc_tmp[iS] = (TH3F*)list->FindObject(Form("fHistGenPartsAcc%c",kLetter[iS]));
    Requires(fHistGenPartsAcc_tmp[iS],Form("fHistGenPartsAcc%c",kLetter[iS]));
    fHistGenPartsAcc_tmp[iS]->SetDirectory(0);
    fHistGenPartsINELgt0_tmp[iS] = (TH3F*)list->FindObject(Form("fHistGenPartsINELgt0%c",kLetter[iS]));
    Requires(fHistGenPartsINELgt0_tmp[iS],Form("fHistGenPartsINELgt0%c",kLetter[iS]));
    fHistGenPartsINELgt0_tmp[iS]->SetDirectory(0);
  }


  for(int iSpecies=0; iSpecies<kNspc; iSpecies++){
    for(int iMatter=0; iMatter<2; iMatter++){
      TList *out_list = new TList();
      out_list->SetName(Form("%s%s",prtName_str[iSpecies],suffix[iMatter]));
      out_list->SetOwner(true);
      for(int iCentrality=0; iCentrality<kNCentBins; iCentrality++){

        TH1F* fHistPtAccepted = (TH1F* )fHistGenPartsAcc_tmp[iMatter]->ProjectionY(Form("Accepted_%s%s_%s",prtName_str[iSpecies],suffix[iMatter],kMultLab[iCentrality]),lowCentBin[iCentrality],upCentBin[iCentrality],iSpecies+1,iSpecies+1);
        fHistPtAccepted->Sumw2();
        TH1F* fHistPtAcceptedNorm = (TH1F*) fHistPtAccepted->Clone("fHistPtAcceptedNorm");
        fHistPtAcceptedNorm->Scale(1./NormAcc[iCentrality]);

        TH1F* fHistPtTrueINELgt0 = (TH1F* )fHistGenPartsINELgt0_tmp[iMatter]->ProjectionY(Form("TrueINEL_%s%s_%s",prtName_str[iSpecies],suffix[iMatter],kMultLab[iCentrality]),lowCentBin[iCentrality],upCentBin[iCentrality],iSpecies+1,iSpecies+1);
        fHistPtTrueINELgt0->Sumw2();
        TH1F* fHistPtTrueINELgt0Norm = (TH1F*) fHistPtTrueINELgt0->Clone("fHistPtTrueINELgt0Norm");
        fHistPtTrueINELgt0Norm->Scale(1./NormTrue[iCentrality]);

        TH1D *fHistSigLoss = (TH1D*)fHistPtTrueINELgt0->Clone(Form("sgnLoss_%s%s_%s",prtName_str[iSpecies],suffix[iMatter],kMultLab[iCentrality]));
        fHistSigLoss->SetTitle(";#it{p}_{T} (GeV/#it{c});1/#epsilon_{part}");
        fHistSigLoss->Divide(fHistPtTrueINELgt0, fHistPtAccepted, 1, 1, "B");
        fHistSigLoss->SetLineColor(iSpecies + 1);

        TH1D *fHistSigLossNorm = (TH1D*)fHistPtTrueINELgt0Norm->Clone(Form("sgnLossNorm_%s%s_%s",prtName_str[iSpecies],suffix[iMatter],kMultLab[iCentrality]));
        fHistSigLossNorm->SetTitle(";#it{p}_{T} (GeV/#it{c});#epsilon_{evt}/#epsilon_{part}");
        fHistSigLossNorm->Divide(fHistPtTrueINELgt0Norm, fHistPtAcceptedNorm, 1, 1, "B");
        fHistSigLossNorm->SetLineColor(iSpecies + 1);

        TH1D* fHistCheck = (TH1D*)fHistPtAcceptedNorm->Clone(Form("fHistCheck_%s%s_%s",prtName_str[iSpecies],suffix[iMatter],kMultLab[iCentrality]));
        fHistCheck->Multiply(fHistSigLoss);
        fHistCheck->Divide(fHistPtTrueINELgt0Norm);

        out_list->Add(fHistSigLoss);
        out_list->Add(fHistSigLossNorm);
        out_list->Add(fHistCheck);
      }
      output_file.cd();
      out_list->Write("", TObject::kSingleKey);
    }
  }

}
