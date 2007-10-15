
#include "TTree.h"
#include "TChain.h"
#include "TClonesArray.h"
#include "AliHLTPHOSDebugRawDigit.h"
#include <iostream>
#include "TH1I.h"
#include "TH2I.h"
#include "TFile.h"

using namespace std;

Int_t sanityFix(const char*);

int main(int argc, const char** argv)
{
  sanityFix(argv[1]);
  return 0;
}

Int_t sanityFix(const char* runNb)
{
  char filepath [50];
  char outfile [50];

  TChain *tree= new TChain("digitTree");

  sprintf(filepath, "/tmp/phoshlt/analysis/data/run%s/*", runNb);
  sprintf(outfile, "/home/phoshlt/analysis/output/run%s/crazyness.root", runNb);
  tree->Add(filepath);
 
  TClonesArray *digArray = new TClonesArray("AliHLTPHOSDebugRawDigit" , 100);
  //tree->SetBranchAddress("digits", &digArray);
  tree->SetBranchAddress("DebugRawDigit", &digArray);
  TH1I *crazyHist = new TH1I("crazyHist", "Crazyness", 27, -6, 20);
  TH2I *crazyChannelHist = new TH2I("crazyChannelHist", "Crazy channels", 64, 0, 63, 56, 0, 55);
  AliHLTPHOSDebugRawDigit *digit = 0;

  cout << endl << "Tree with " << tree->GetEntries() << " events...\n";
  for(int k = 0; k < tree->GetEntries(); k++)
    {
      tree->GetEntry(k);
      for(int j = 0; j < digArray->GetEntriesFast(); j++)
	{
	  digit = (AliHLTPHOSDebugRawDigit*)digArray->At(j);
	  crazyHist->Fill(digit->GetCrazyness());
	  if(digit->GetCrazyness() == 0 && digit->GetGain() == 1)
	    {
	      crazyChannelHist->Fill(digit->GetX(), digit->GetZ());
	    }
	}
    }
  cout << "Printing file...";
  TFile *file = new TFile(outfile,"recreate");
  crazyHist->Write();
  crazyChannelHist->Write();
  file->Close();
  delete file;
  file = 0;
  cout << "Done!\n";
  
  return 0;
}
  
 
