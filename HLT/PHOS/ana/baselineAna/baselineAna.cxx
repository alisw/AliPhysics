
#include "TTree.h"
#include "TChain.h"
#include "TClonesArray.h"
#include "AliHLTPHOSBaseline.h"
#include <iostream>
#include "TH1F.h"
#include "TFile.h"

using namespace std;

Int_t baselineAna(const char*);

int main(int argc, const char** argv)
{
  baselineAna(argv[1]);
  return 0;
}

Int_t baselineAna(const char* runNb)
{
  char filepath [50];
  char outfile [50];

  bool finished = false;
  
  TChain *tree= new TChain("baselineTree");

  sprintf(filepath, "/tmp/phoshlt/analysis/data/baselines/run%s/*", runNb);
  sprintf(outfile, "/home/phoshlt/analysis/output/run%s/baselines.root", runNb);
  tree->Add(filepath);
 
  TClonesArray *baselineArray = new TClonesArray("AliHLTPHOSBaseline" , 100);
  tree->SetBranchAddress("Baselines", &baselineArray);
  TH1F *baselineHist = new TH1F("baselineHist", "Baselines", 100, 0, 99);
  AliHLTPHOSBaseline *baseline = 0;

  cout << endl << "Tree with " << tree->GetEntries() << " events...\n";
  for(int k = 0; k < tree->GetEntries(); k++)
    {
      tree->GetEntry(k);
      for(int j = 0; j < baselineArray->GetEntriesFast(); j++)
	{
	  baseline = (AliHLTPHOSBaseline*)baselineArray->At(j);
	  baselineHist->Fill(baseline->GetBaseline());
	}
    }
  
  
  cout << "Printing file...";

  
  TFile *file = new TFile(outfile,"recreate");
  baselineHist->Write();
  file->Close();
  delete file;
  file = 0;
  cout << "Done!\n";
  
  return 0;
}
  
 
