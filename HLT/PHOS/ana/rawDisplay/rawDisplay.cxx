
#include "TTree.h"
#include "TChain.h"
#include "TClonesArray.h"
#include "AliHLTPHOSDebugRawDigit.h"
#include <iostream>
#include "TH2D.h"
#include "TFile.h"

using namespace std;

Int_t rawDisplay(const char*);

int main(int argc, const char** argv)
{
  rawDisplay(argv[1]);
  return 0;
}

Int_t rawDisplay(const char* runNb)
{
  char filepath [50];
  char outfile [50];

  bool finished = false;
  
  TChain *tree= new TChain("digitTree");

  sprintf(filepath, "/tmp/phoshlt/analysis/data/run%s/*", runNb);
  sprintf(outfile, "/home/phoshlt/analysis/output/run%s/rawHist.root", runNb);
  tree->Add(filepath);
 
  TClonesArray *digArray = new TClonesArray("AliHLTPHOSDebugRawDigit" , 100);
   tree->SetBranchAddress("DebugRawDigit", &digArray);
  TH1I *rawHist = new TH1I("rawHist", "Raw display", 70, 0, 69);
  AliHLTPHOSDebugRawDigit *digit = 0;

  cout << endl << "Tree with " << tree->GetEntries() << " events...\n";
  /*
  for(int k = 0; k < tree->GetEntries(); k++)
    {
      tree->GetEntry(k);
      for(int j = 0; j < digArray->GetEntriesFast(); j++)
	{
	  digit = (AliHLTPHOSDebugRawDigit*)digArray->At(j);
	  //if(digit->GetAmplitude() - digit->GetBaseline() > 200)
	    //	  if(digit->GetBaseline() && digit->GetCrazyness())
	  if(digit->GetAmplitude() < 30 && digit->GetAmplitude() > 25 )
	    {
	      UInt_t *data = digit->GetRawData();
	      for(Int_t i = 0; i < 70; i++)
		{
		  rawHist->SetBinContent(i, (digit->GetRawData())[i]);
		}
	      finished = true;
	      break;
	    }
	}
      if(finished) break;
      
      }*/

  tree->GetEntry(3);
  for(int j = 0; j < digArray->GetEntriesFast(); j++)
    {
      digit = (AliHLTPHOSDebugRawDigit*)digArray->At(j);
      if(digit->GetAmplitude() > 60)
	{
	  cout << digit->GetX() << " - " << digit->GetZ() << " - " << digit->GetAmplitude() << endl;
	}
      if(digit->GetX() == 57 && digit->GetZ() == 55 && digit->GetGain() == 0) 
	{
	  for(int a = 0; a < 70; a++)
	    {
	      cout << (digit->GetRawData())[a] << endl;
	      rawHist->SetBinContent(a, (digit->GetRawData())[a]);
	    }
	}
    }
  
  
  cout << "Printing file...";

  
  TFile *file = new TFile(outfile,"recreate");
  rawHist->Write();
  file->Close();
  delete file;
  file = 0;
  cout << "Done!\n";
  
  return 0;
}
  
 
