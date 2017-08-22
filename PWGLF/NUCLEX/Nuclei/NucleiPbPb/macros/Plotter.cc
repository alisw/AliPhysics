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

string input_name = "results/signal.root";
string output_name = "results/plots.root";

string title = "ciaociao";

void Plotter(){

  TFile file_in(input_name.data());
  TFile file_out(output_name.data(),"RECREATE");

  TCanvas* fCanvas= nullptr;
  int iPad = 0;
  int nPads = 11;

  for (auto list_key : *file_in.GetListOfKeys()) {
    /// Preliminary operation to read the list and create an output dir
    if (string(list_key->GetName()).find(kFilterListNames.data()) == string::npos) continue;

    TTList* list = (TTList*)file_in.Get(list_key->GetName());
    TDirectory* base_dir = file_out.mkdir(list_key->GetName());
    file_out.cd(list_key->GetName());
    std::cout << string(list_key->GetName()).data() << std::endl;
    for(int iS=0; iS<2; iS++){
      std::cout << "   " << kNames[iS] << std::endl;
      TDirectory* dir = base_dir->mkdir(kNames[iS].data());
      dir->cd();
      for(int iC=0; iC<kNCentralityBins; iC++){
        std::cout << "      " << iC << std::endl;
        dir->mkdir(Form("cent_%d",iC));
        dir->cd(Form("cent_%d",iC));
        iPad=0;
        for(int iB=5; iB<=15; iB++){
          if(iPad%9 == 0){
            if(fCanvas) delete fCanvas;
            fCanvas = new TCanvas(Form("Canvas_%d",iPad/6),Form("Canvas_%d",iPad/6),3200,2400);
            fCanvas->Divide(3,2);
          }
          //if(iPad==0) fCanvas->Print(Form("cent_%d_%c.pdf[",iC,kLetter[iS]));
          fCanvas->cd(iPad%9+1);
          string path = string(list_key->GetName()) + "/" + kNames[iS] + "/TailTail";
          RooPlot* fPlot = (RooPlot*)file_in.Get(Form("%s/d%i_%i",path.data(),iC,iB));
          fPlot->Draw();
          if((iPad+1)%9 == 0){
            //fCanvas->Print(Form("cent_%d_%c.pdf",iC,kLetter[iS]));
            fCanvas->Write();
          }
          iPad++;
          if(iPad==nPads){
            //fCanvas->Print(Form("cent_%d_%c.pdf]",iC,kLetter[iS]));
            fCanvas->Write();
          }
        }
      }
    }
  }
}
