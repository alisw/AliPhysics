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

string output_name = "../results/plots.root";

void Plotter(bool short_mode = true){

  TFile file_in(kSignalOutput.data());
  TFile file_out(output_name.data(),"RECREATE");

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
    std::cout << string(list_key->GetName()).data() << std::endl;
    for(int iS=0; iS<2; iS++){
      std::cout << "   " << kNames[iS] << std::endl;
      TDirectory* dir = base_dir->mkdir(kNames[iS].data());
      dir->cd();
      for(int iC=0; iC<kCentLength; iC++){
        std::cout << "      " << iC << std::endl;
        dir->mkdir(Form("cent_%d",iC));
        dir->cd(Form("cent_%d",iC));
        iPad = 0;
        nPads = kPtBinLimit[iC] - 5 + 1;
        for(int iB=5; iB<=kPtBinLimit[iC]; iB++){
          if(iPad%6 == 0){
            if(fCanvas) delete fCanvas;
            fCanvas = new TCanvas(Form("Canvas_%d",iPad/6),Form("Canvas_%d",iPad/6),3200,2400);
            fCanvas->Divide(3,2);
          }
          //if(iPad==0) fCanvas->Print(Form("cent_%d_%c.pdf[",iC,kLetter[iS]));
          fCanvas->cd(iPad%6+1);
          string path = string(list_key->GetName()) + "/" + kNames[iS] + "/TailTail/C_" + to_string(iC) ;
          std::cout << "path " << path.data() << std::endl;
          RooPlot* fPlot = (RooPlot*)file_in.Get(Form("%s/d%i_%i",path.data(),iC,iB));
          Requires(fPlot,"RooPlot");
          fPlot->Draw();
          if((iPad+1)%6 == 0){
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
    break;
  }
}
