#include "src/Common.h"
#include "src/Utils.h"
using namespace utils;
#include <TFile.h>
#include <RooPlot.h>
#include <TCanvas.h>
#include <TDirectory.h>
#include <iostream>

#include <string>
using std::string;

using namespace RooFit;

void Plotter(bool bTPC = false, bool short_mode = true){

  TFile file_in(kSignalOutput.data());
  string input_dir = (bTPC) ? "/TPConly" : "/Fits";

  string output_name = (bTPC) ? Form("%splotsTPC.root",kBaseOutputDir.data()) : Form("%splots.root",kBaseOutputDir.data());
  TFile file_out(output_name.data(),"RECREATE");

  string complete_output_pdf_string = kBaseOutputDir + "plots/all_plots_%s.pdf";
  const char* complete_output_pdf = complete_output_pdf_string.data();

  string complete_output_pdf_string_open = kBaseOutputDir + "plots/all_plots_%s.pdf[";
  const char* complete_output_pdf_open = complete_output_pdf_string_open.data();

  string complete_output_pdf_string_close = kBaseOutputDir + "plots/all_plots_%s.pdf]";
  const char* complete_output_pdf_close = complete_output_pdf_string_close.data();

  const char* kDetector = (bTPC) ? "TPC" : "TOF";

  TCanvas* fCanvas= nullptr;
  int iPad = 0;
  int nPads = 0;

  int counter = 0;

  for (auto list_key : *file_in.GetListOfKeys()) {
    /// Preliminary operation to read the list and create an output dir
    if(short_mode)
      if(counter>0) return;
    counter++;
    if (string(list_key->GetName()).find(kFilterListNames.data()) == string::npos) continue;

    TTList* list = (TTList*)file_in.Get(list_key->GetName());
    TDirectory* base_dir = file_out.mkdir(list_key->GetName());
    file_out.cd(list_key->GetName());
    for(int iS=0; iS<2; iS++){
      TDirectory* dir = base_dir->mkdir(kNames[iS].data());
      dir->cd();
      for(int iC=0; iC<kCentLength; iC++){
        dir->mkdir(Form("cent_%d",iC));
        dir->cd(Form("cent_%d",iC));
        iPad = 0;
        int nPads = 0;
        for(int iB = 1; iB <= kNPtBins; ++iB){
          float bin_center = (kPtBins[iB]+kPtBins[iB-1])/2;
          if(bTPC){
            if(bin_center>kTPCmaxPt) break;
          }
          else{
            if(bin_center<kTOFminPt) continue;
          }
          if(bin_center>kCentPtLimits[iC]) break;
          nPads++;
        }
        int page_counter = 0;
        for(int iB = 1; iB <= kNPtBins; ++iB){
          float bin_center = (kPtBins[iB]+kPtBins[iB-1])/2;
          if(bTPC){
            if(bin_center>kTPCmaxPt) break;
          }
          else{
            if(bin_center<kTOFminPt) continue;
          }
          if(bin_center>kCentPtLimits[iC]) break;
          if(iPad%6 == 0){
            if(fCanvas) delete fCanvas;
            fCanvas = new TCanvas(Form("Canvas_%d",iPad/6),Form("Canvas_%d",iPad/6),3200,2400);
            if(iS==0 && iC==0 && iPad==0){
              fCanvas->Print(Form(complete_output_pdf_open,kDetector));
            }
            fCanvas->Divide(2,2);
          }
          fCanvas->cd(iPad%4+1);
          string path = string(list_key->GetName()) + "/" + kNames[iS] + input_dir.data();
          RooPlot* fPlot = (bTPC) ? (RooPlot*)file_in.Get(Form("%s/C_%d/TPC_d%i_%i",path.data(),iC,iC,iB)) : (RooPlot*)file_in.Get(Form("%s/C_%d/d%i_%i",path.data(),iC,iC,iB));
          (bTPC) ? Requires(fPlot,Form("%s/TPC_d_%i_%i",path.data(),iC,iB)) : Requires(Form("%s/C_%d/d%i_%i",path.data(),iC,iC,iB));
          fPlot->Draw();
          if((iPad+1)%4 == 0){
            fCanvas->Print(Form("%splots/cent_%d_%c_%d.pdf",kBaseOutputDir.data(),iC,kLetter[iS],page_counter));
            fCanvas->Write();
            fCanvas->Print(Form(complete_output_pdf,kDetector));
            page_counter++;
            if((iPad+1)==nPads && iS==1 && iC==kCentLength-1 ){
              fCanvas->Print(Form(complete_output_pdf_close,kDetector));
            }
          }
          iPad++;
          if(iPad==nPads && (iPad)%4 != 0){
            fCanvas->Print(Form("%splots/cent_%d_%c_%d.pdf",kBaseOutputDir.data(),iC,kLetter[iS],page_counter));
            fCanvas->Write();
            fCanvas->Print(Form(complete_output_pdf,kDetector));
            if(iS==1 && iC==kCentLength-1)fCanvas->Print(Form(complete_output_pdf_close,kDetector));
          }
        }
      }
    }
  }
}
