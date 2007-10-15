
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

Int_t energySpectrum(const char* filename)
{
  
  
  TH1F *spectrumHist = new TH1F("spectrumHist", "Digit energy spectrum", 200, 0, 200);
  TClonesArray *digArray = new TClonesArray("AliHLTPHOSDebugRawDigit" , 100);
  TChain *tree= new TChain("digitTree");
  AliHLTPHOSDebugRawDigit *digit = 0;

  tree->Add(filename);
  
  tree->SetBranchAddress("DebugRawDigit", &digArray);
  cout << endl << "Entries in tree: " << tree->GetEntries() << endl;

  for(int k = 0; k < tree->GetEntries(); k++)
    {
      tree->GetEntry(k);
      for(int j = 0; j < digArray->GetEntriesFast(); j++)
	{
	  digit = (AliHLTPHOSDebugRawDigit*)digArray->At(j);
	  spectrumHist->Fill(digit->GetAmplitude());
	}
    }
  
  spectrumHist->Draw();

  return 0;
}
  
  
