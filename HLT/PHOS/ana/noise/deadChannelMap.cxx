
#include "TTree.h"
#include "TChain.h"
#include "TClonesArray.h"
#include "AliHLTPHOSDebugRawDigit.h"
#include <iostream>
#include "TH2D.h"
#include "TFile.h"

using namespace std;

Int_t deadChannelMap(const char*);

int main(int argc, const char** argv)
{
  deadChannelMap(argv[1]);
  return 0;
}

Int_t deadChannelMap(const char* runNb)
{
  char filepath [50];
  char outfile [50];
  
  TChain *tree= new TChain("digitTree");

  sprintf(filepath, "/tmp/phoshlt/analysis/data/run%s/*", runNb);
  sprintf(outfile, "/home/phoshlt/analysis/output/run%s/deadMap.root", runNb);
  tree->Add(filepath);
 
  TClonesArray *digArray = new TClonesArray("AliHLTPHOSDebugRawDigit" , 100);
  //tree->SetBranchAddress("digits", &digArray);
  tree->SetBranchAddress("DebugRawDigit", &digArray);
  TH2D *deadMap = new TH2D("deadMap", "Dead Channel Map", 64, 0, 63, 56, 0, 55);
  AliHLTPHOSDebugRawDigit *digit = 0;

  cout << endl << "Chain with " << tree->GetEntries() << " events...\n";
  for(int k = 0; k < tree->GetEntries(); k++)
    {
      tree->GetEntry(k);
      //  cout << "     Event with " << digArray->GetEntriesFast() << " digits\n";
      for(int j = 0; j < digArray->GetEntriesFast(); j++)
	{
	  digit = (AliHLTPHOSDebugRawDigit*)digArray->At(j);
	  for(int i = 0; i < 70; i++)
	    {
	      if((digit->GetRawData())[i] > 55)
		{
		  if((digit->GetRawData())[i+1] > 55)
		    {
		      if((digit->GetRawData())[i+2] > 55)
			{
			  if((digit->GetRawData())[i+3] > 55)
			    {
			      deadMap->SetBinContent(digit->GetX(), digit->GetZ(), 10);
			      break;
			    }
			}
		    }
		}
	    }
	}
    }

  
  cout << "Printing file...";


  TFile *file = new TFile(outfile,"recreate");
  deadMap->Write();
  file->Close();
  delete file;
  file = 0;
  cout << "Done!\n";


  
  return 0;
}
  
  
