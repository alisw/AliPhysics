
#include "TTree.h"
#include "TChain.h"
#include "TClonesArray.h"
#include "AliHLTPHOSDebugRawDigit.h"
#include <iostream>

using namespace std;

Int_t anaRaw();

int main()
{
  anaRaw();
  return 0;
}

Int_t anaRaw()
{
  
  
  
  TClonesArray *digArray = new TClonesArray("AliHLTPHOSDebugRawDigit" , 100);
  TChain *tree= new TChain("digitTree");
  AliHLTPHOSDebugRawDigit *digit = 0;

  tree->Add("/tmp/phoshlt/analysis/data/run7800/*");
  
  tree->SetBranchAddress("DebugRawDigit", &digArray);
  cout << endl << "Entries in tree: " << tree->GetEntries() << endl;

  for(int k = 0; k < tree->GetEntries(); k++)
    {
      tree->GetEntry(k);
      for(int j = 0; j < digArray->GetEntriesFast(); j++)
	{
	  digit = (AliHLTPHOSDebugRawDigit*)digArray->At(j);
            
	  //  if((digit->GetAmplitude() - digit->GetBaseline()) < 50 && (digit->GetAmplitude() - digit->GetBaseline()) > 20)
	  if(digit->GetBaseline() > 70)
	    {
	      cout << "Digit energy: " << digit->GetAmplitude() << endl;
	      cout << "Baseline: " << digit->GetBaseline() << endl;
	      /*
	       if(digit->GetCrazyness())
		{
		  cout << "Crazyness: " << digit->GetCrazyness() << endl;
		  for(int i = 0; i < 70; i++)
		    { 
		      cout << (digit->GetRawData())[i] << endl;
		    }
		}
	      */
	    }	  
	}
    }
  
  return 0;
}
  
  
