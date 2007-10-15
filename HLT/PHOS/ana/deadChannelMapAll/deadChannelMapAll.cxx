
#include "TTree.h"
#include "TChain.h"
#include "TClonesArray.h"
#include "AliHLTPHOSDebugRawDigit.h"
#include <iostream>
#include "TH2D.h"
#include "TFile.h"

using namespace std;

Int_t deadChannelMapAll();

int main(int argc, const char** argv)
{
  deadChannelMapAll();
  return 0;
}

Int_t deadChannelMapAll()
{
  char filepath [50];
  char outfile [50];
  
  TChain *tree= new TChain("digitTree");
  TChain *treeOld = new TChain("digitTree");
  TChain *treeOlder = new TChain("digitTree");

  sprintf(outfile, "/home/phoshlt/analysis/output/all/deadMapAll.root");
  tree->Add("/tmp/phoshlt/analysis/data/run7798/*");
  tree->Add("/tmp/phoshlt/analysis/data/run7799/*");
  tree->Add("/tmp/phoshlt/analysis/data/run7800/*");
  tree->Add("/tmp/phoshlt/analysis/data/run7801/*");
  tree->Add("/tmp/phoshlt/analysis/data/run7802/*");
  tree->Add("/tmp/phoshlt/analysis/data/run7804/*");
  tree->Add("/tmp/phoshlt/analysis/data/run7805/*");
  tree->Add("/tmp/phoshlt/analysis/data/run7806/*");
  
  // tree->Add("/tmp/phoshlt/analysis/data/run7775/*");
 //  tree->Add("/tmp/phoshlt/analysis/data/run7776/*");
  // tree->Add("/tmp/phoshlt/analysis/data/run7777/*");
 //  tree->Add("/tmp/phoshlt/analysis/data/run7778/*");
 //  tree->Add("/tmp/phoshlt/analysis/data/run7780/*");
 //  tree->Add("/tmp/phoshlt/analysis/data/run7781/*");
 //  tree->Add("/tmp/phoshlt/analysis/data/run7785/*");
 //  tree->Add("/tmp/phoshlt/analysis/data/run7786/*");
 //  tree->Add("/tmp/phoshlt/analysis/data/run7789/*");
 //  tree->Add("/tmp/phoshlt/analysis/data/run7790/*");
  // tree->Add("/tmp/phoshlt/analysis/data/run7791/*");
  // tree->Add("/tmp/phoshlt/analysis/data/run7792/*");
 //  tree->Add("/tmp/phoshlt/analysis/data/run7793/*");
 //  tree->Add("/tmp/phoshlt/analysis/data/run7794/*");
  // tree->Add("/tmp/phoshlt/analysis/data/run7795/*");
  // tree->Add("/tmp/phoshlt/analysis/data/run7796/*");
  // tree->Add("/tmp/phoshlt/analysis/data/run7797/*");
  // tree->Add("/home/phoshlt/analysis/data/run7798/*");
  

  
  // treeOld->Add("/tmp/phoshlt/analysis/data_old/run7762/*");
  // treeOld->Add("/tmp/phoshlt/analysis/data_old/run7765/*");
  // treeOld->Add("/tmp/phoshlt/analysis/data_old/run7767/*");
  // treeOld->Add("/tmp/phoshlt/analysis/data_old/run7768/*");
  // treeOld->Add("/tmp/phoshlt/analysis/data_old/run7769/*");
  // treeOld->Add("/tmp/phoshlt/analysis/data_old/run7771/*");
  
  // treeOlder->Add("/tmp/phoshlt/analysis/data_old/run7758/*");
  // treeOlder->Add("/tmp/phoshlt/analysis/data_old/run7760/*");
  // treeOlder->Add("/tmp/phoshlt/analysis/data_old/run7761/*");
  


  TClonesArray *digArray = new TClonesArray("AliHLTPHOSDebugRawDigit" , 100);
  TClonesArray *digArrayOld = new TClonesArray("AliHLTPHOSDebugRawDigit" , 100);
  TClonesArray *digArrayOlder = new TClonesArray("AliHLTPHOSDebugRawDigit" , 100);
  //tree->SetBranchAddress("digits", &digArray);
  tree->SetBranchAddress("DebugRawDigit", &digArray);
  treeOld->SetBranchAddress("DebugDigitsRaw", &digArrayOld);
  treeOlder->SetBranchAddress("digits", &digArrayOlder);

  TH2D *deadMapAll = new TH2D("deadMapAll", "Dead Channel Map All runs", 64, 0, 63, 56, 0, 55);
  AliHLTPHOSDebugRawDigit *digit = 0;

  cout << endl << "tree is a cain with " << tree->GetEntries() << " events...\n";
 //  cout << "treeOld is a cain with " << treeOld->GetEntries() << " events...\n";
 //  cout << "treeOlder is a cain with " << treeOlder->GetEntries() << " events...\n";
 //  cout << "Gives a total of " << tree->GetEntries() + treeOld->GetEntries() + treeOlder->GetEntries() << "events\n";

  for(int k = 0; k < tree->GetEntries(); k++)
    {
      tree->GetEntry(k);
      //  cout << "     Event with " << digArray->GetEntriesFast() << " digits\n";
      for(int j = 0; j < digArray->GetEntriesFast(); j++)
	{
	  digit = (AliHLTPHOSDebugRawDigit*)digArray->At(j);
	  for(int i = 0; i < 70; i++)
	    {
	      if((digit->GetRawData())[i] > 60)
		{
		  if((digit->GetRawData())[i+1] > 60)
		    {
		      if((digit->GetRawData())[i+2] > 60)
			{
			  if((digit->GetRawData())[i+3] > 60)
			    {
			      deadMapAll->SetBinContent(digit->GetX(), digit->GetZ(), 10);
			      break;
			    }
			}
		    }
		}
	    }
	}
    }
 for(int k = 0; k < treeOld->GetEntries(); k++)
    {
      treeOld->GetEntry(k);
      //  cout << "     Event with " << digArray->GetEntriesFast() << " digits\n";
      for(int j = 0; j < digArrayOld->GetEntriesFast(); j++)
	{
	  digit = (AliHLTPHOSDebugRawDigit*)digArrayOld->At(j);
	  for(int i = 0; i < 70; i++)
	    {
	      if((digit->GetRawData())[i] > 30)
		{
		  if((digit->GetRawData())[i+1] > 30)
		    {
		      if((digit->GetRawData())[i+2] > 30)
			{
			  if((digit->GetRawData())[i+3] > 30)
			    {
			      deadMapAll->SetBinContent(digit->GetX(), digit->GetZ(), 10);
			      break;
			    }
			}
		    }
		}
	    }
	}
    }
   for(int k = 0; k < treeOlder->GetEntries(); k++)
    {
      treeOlder->GetEntry(k);
      //  cout << "     Event with " << digArray->GetEntriesFast() << " digits\n";
      for(int j = 0; j < digArrayOlder->GetEntriesFast(); j++)
	{
	  digit = (AliHLTPHOSDebugRawDigit*)digArrayOlder->At(j);
	  for(int i = 0; i < 70; i++)
	    {
	      if((digit->GetRawData())[i] > 30)
		{
		  if((digit->GetRawData())[i+1] > 30)
		    {
		      if((digit->GetRawData())[i+2] > 30)
			{
			  if((digit->GetRawData())[i+3] > 30)
			    {
			      deadMapAll->SetBinContent(digit->GetX(), digit->GetZ(), 10);
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
  deadMapAll->Write();
  file->Close();
  delete file;
  file = 0;
  cout << "Done!\n";


  
  return 0;
}
  
  
