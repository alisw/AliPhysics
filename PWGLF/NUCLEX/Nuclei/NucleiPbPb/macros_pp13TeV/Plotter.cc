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

  TCanvas* fCanvas= nullptr;
  int iPad = 0;
  int nPads = 0;

  const int kPtBinLimit[kCentLength] = {14,14,13,13,12,12,11,10,10,15};

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
        nPads = kPtBinLimit[iC] - 5 + 1;
        int page_counter = 0;
        for(int iB=1; iB<=kPtBinLimit[iC]; iB++){
          if(bTPC){
            if(iB>6) continue;
          }
          else{
            if(iB<5) continue;
          }
          if(iPad%6 == 0){
            if(fCanvas) delete fCanvas;
            fCanvas = new TCanvas(Form("Canvas_%d",iPad/6),Form("Canvas_%d",iPad/6),3200,2400);
            fCanvas->Divide(3,2);
          }
          fCanvas->cd(iPad%6+1);
          string path = string(list_key->GetName()) + "/" + kNames[iS] + input_dir.data();
          RooPlot* fPlot = (bTPC) ? (RooPlot*)file_in.Get(Form("%s/TPC_d_%i_%i",path.data(),iC,iB)) : (RooPlot*)file_in.Get(Form("%s/C_%d/d%i_%i",path.data(),iC,iC,iB));
          (bTPC) ? Requires(fPlot,Form("%s/TPC_d_%i_%i",path.data(),iC,iB)) : Requires(Form("%s/C_%d/d%i_%i",path.data(),iC,iC,iB));
          fPlot->Draw();
          if((iPad+1)%6 == 0){
            fCanvas->Print(Form("%splots/cent_%d_%c_%d.pdf",kBaseOutputDir.data(),iC,kLetter[iS],page_counter));
            fCanvas->Write();
            page_counter++;
          }
          iPad++;
          if(iPad==nPads && (iPad)%6 != 0){
            fCanvas->Print(Form("%splots/cent_%d_%c_%d.pdf",kBaseOutputDir.data(),iC,kLetter[iS],page_counter));
            fCanvas->Write();
          }
        }
      }
    }
  }
}
