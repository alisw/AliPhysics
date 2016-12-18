#include <string>
using std::string;
#include <iostream>

#include <TFile.h>
#include <TH1D.h>
#include <TH3D.h>
#include <TList.h>
#include <TMath.h>

string names[] = {"#pi","k","p","n","d"};


///////// Functions in Utils.h moved here for AliPhysics version
template<class T> void Requires(T* obj, string msg = "") {
  if (!obj) {
    std::cout << "Missing object \"" << msg.data() << "\". Aborting." << std::endl;
    abort();
  }
}

template<typename T, typename L> void Get(T*& dest, L* from, string name) {
  dest = static_cast<T*>(from->Get(name.data()));
  Requires(dest, name);
}

template<typename T> void Get(T*& dest, TList* from, string name) {
  dest = static_cast<T*>(from->FindObject(name.data()));
  Requires(dest, name);
}

template<typename T, typename L> void Get(T*& dest, L& from, string name) {
  Get(dest,&from,name);
}
///////////

void AnalyseOutput() {
  TFile input_file("merged_output_nom.root");
  TList* input_list;
  Get(input_list, input_file, "fOutputList");
  TH1D* norm;
  TH3D* data;
  Get(norm, input_list, "fEventCounter");
  Get(data, input_list, "fPtSpectra");

  TFile output("outputfile.root","recreate");
  for (int iC = 1; iC <= data->GetNbinsX(); ++iC) {
    double normalisation = norm->GetBinContent(iC);
    if (normalisation < 1.f) continue;
    normalisation = 1. / normalisation;
    output.mkdir(Form("Cent%i",iC));
    output.cd(Form("Cent%i",iC));
    for (size_t iS = 1; iS <= 5; ++iS) {
      TH1D* pr0 = data->ProjectionZ(Form("%splus",names[iS - 1].data()),iC,iC,iS * 2 - 1, iS * 2 - 1);
      pr0->Scale(normalisation);
      TH1D* pr1 = data->ProjectionZ(Form("%sminus",names[iS - 1].data()),iC,iC,iS * 2, iS * 2);
      pr1->Scale(normalisation);
      TH1D* sum = (TH1D*)pr0->Clone(Form("%ssum",names[iS - 1].data()));
      sum->Add(pr1);
      sum->Scale(0.5);
      for (int iB = 1; iB <= pr0->GetNbinsX(); ++iB) {
        pr0->SetBinContent(iB, pr0->GetBinContent(iB) / (TMath::TwoPi() * pr0->GetBinCenter(iB))); 
        pr1->SetBinContent(iB, pr1->GetBinContent(iB) / (TMath::TwoPi() * pr0->GetBinCenter(iB))); 
        sum->SetBinContent(iB, sum->GetBinContent(iB) / (TMath::TwoPi() * pr0->GetBinCenter(iB))); 
      }
      pr0->Write();
      pr1->Write();
      sum->Write();
    }
  }
}
