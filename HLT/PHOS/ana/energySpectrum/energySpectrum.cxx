
#include "TTree.h"
#include "TChain.h"
#include "TClonesArray.h"
#include "TH1F.h"
#include "TFile.h"
#include "AliHLTPHOSDebugRawDigit.h"
#include <iostream>

using namespace std;

Int_t energySpectrum(const char*);

int main(int argc, const char** argv)
{
  energySpectrum(argv[1]);
  return 0;
}

Int_t energySpectrum(const char* runNb)
{
  
  
  TH1F *spectrumHist = new TH1F("spectrumHist", "Digit energy spectrum", 200, 0, 200);
  TClonesArray *digArray = new TClonesArray("AliHLTPHOSDebugRawDigit" , 100);
  TChain *tree= new TChain("digitTree");
  AliHLTPHOSDebugRawDigit *digit = 0;

  char filepath [50];
  char outfile [50];

  sprintf(filepath, "/tmp/phoshlt/analysis/data/run%s/*", runNb);
  sprintf(outfile, "/home/phoshlt/analysis/output/run%s/spectrum.root", runNb);
  tree->Add(filepath);
  
  tree->SetBranchAddress("DebugRawDigit", &digArray);
  cout << endl << "Entries in tree: " << tree->GetEntries() << endl;

  for(int k = 0; k < tree->GetEntries(); k++)
    {
      tree->GetEntry(k);
      for(int j = 0; j < digArray->GetEntriesFast(); j++)
	{
	  digit = (AliHLTPHOSDebugRawDigit*)digArray->At(j);
	  if((digit->GetAmplitude() - digit->GetBaseline()) > 20 && digit->GetCrazyness())
	    spectrumHist->Fill(digit->GetAmplitude()- digit->GetBaseline());
	}
    }
  cout << "Printing file...";
  TFile *file = new TFile(outfile,"recreate");
  spectrumHist->Write();
  file->Close();
  delete file;
  file = 0;
  cout << "Done!\n";

  return 0;
}
  
  
