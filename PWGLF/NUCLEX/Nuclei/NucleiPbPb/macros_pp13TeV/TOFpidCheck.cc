#include "src/Common.h"
#include "src/Utils.h"
using namespace utils;

#include <TFile.h>
#include <TDirectory.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>

const char* dir_names[3] = {"nuclei_deuterons_","nuclei_deuterons_pid0","nuclei_deuterons_pid1"};
const char* leg_labels[3] = {"3.00 #sigma_{TPC}","3.25 #sigma_{TPC}","3.50 #sigma_{TPC}"};
const int vColors[3] = {kBlack, kRed, kBlue};


void TOFpidCheck() {
  gStyle->SetOptStat(0);
  /// Taking all the histograms from the MC file
  TFile signal_file(kSignalOutput.data());
  TFile output_file(Form("%sTOFpidCheck.root",kBaseOutputDir.data()),"recreate");
  TH1D* vRawTOF[3] = {nullptr};
  TH1D* vRawTPC[3] = {nullptr};
  for (int iS = 0; iS < 2; ++iS) {
    TDirectory* particle_dir = output_file.mkdir(kNames[iS].data());
    particle_dir->cd();
    for (int iC=0; iC<kCentLength; iC++) {
      TDirectory* cent_dir = particle_dir->mkdir(Form("%d",iC));
      cent_dir->cd();
      /// Getting raw signals
      TCanvas* can = new TCanvas(Form("%s_%c",kNames[iS].data(),iC));
      TLegend* l = new TLegend(0.5,0.5,0.7,0.7,Form("%s, %s",kNames[iS].data(),kRomanLabels[iC]),"brNDC");
      for(int i=0; i<3; i++){
        TH1D* rawTOF_tmp = (TH1D*)signal_file.Get(Form("%s/%s/Fits/hRawCounts%c%i",dir_names[i],kNames[iS].data(),kLetter[iS],iC));
        vRawTOF[i] = (TH1D*)rawTOF_tmp->Clone(Form("hRawTOF_%c_%i_pidvar%i",kLetter[iS],iC,i));
        vRawTOF[i]->SetLineColor(vColors[i]);
        vRawTOF[i]->SetMarkerColor(vColors[i]);
        l->AddEntry(vRawTOF[i],leg_labels[i],"PE");
        Requires(rawTOF_tmp,"Missing TOF raw counts");
        vRawTOF[i]->Write();
        can->cd();
        if(!i){
          vRawTOF[i]->Draw("PE");
        } else {
          vRawTOF[i]->Draw("SAME PE");
        }
      }
      l->Draw();
      cent_dir->cd();
      can->Write();
      if(iS==0 && iC==0){
        can->Print(Form("%splots/pid_check.pdf[",kBaseOutputDir.data()),"pdf");
      }
      can->Print(Form("%splots/pid_check.pdf",kBaseOutputDir.data()),"pdf");
      if(iS==1 && iC==kCentLength-1){
        can->Print(Form("%splots/pid_check.pdf]",kBaseOutputDir.data()),"pdf");
      }
    }
  }
  output_file.Close();
}
