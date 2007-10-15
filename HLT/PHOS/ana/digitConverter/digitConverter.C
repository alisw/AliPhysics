#include "TChain.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TFile.h"
#include "../digits/AliHLTPHOSDebugRawDigit.h"
#include "../digits/AliHLTPHOSDigit.h"
#include "../digits/AliHLTPHOSAltroConfig.h"
#include <iostream>

Int_t digitConverter(Int_t runNb, Int_t segNb)
{
  gROOT->ProcessLine(".L ../digits/AliHLTPHOSDigit.cxx+");
  gROOT->ProcessLine(".L ../digits/AliHLTPHOSDebugRawDigit.cxx+");
  char filename [50];
  sprintf(filename, "/tmp/phoshlt/aldaqpc019/hlt/run%d_digitTree_%d.root", runNb, segNb);

  char filepath [50];
 
  bool finished = false;
  
  TChain *tree= new TChain("digitTree");
  TTree *newTree = new TTree("digitTree", "Digits tree" );

  tree->Add(filename);
  TClonesArray *rawDigArray = new TClonesArray("AliHLTPHOSDebugRawDigit" , 100);
  TClonesArray *oldDigArray = new TClonesArray("AliHLTPHOSDigit" , 100);
  TClonesArray *newDigArray = new TClonesArray("AliHLTPHOSDigit" , 100);
  
  newTree->Branch("Digit", &newDigArray);
  tree->SetBranchAddress("DebugRawDigit", &digArray);
  AliHLTPHOSDebugRawDigit *digit = 0;
  AliHLTPHOSDigit *newDigit = 0;
  cout << endl << "Opening file: " << filename << endl;
  cout << "Converting tree with " << tree->GetEntries() << " events\n";
  
  for(int k = 0; k < tree->GetEntries(); k++)
    {
      tree->GetEntry(k);
      for(int j = 0; j < digArray->GetEntriesFast(); j++)
	{
	  digit = (AliHLTPHOSDebugRawDigit*)digArray->At(j);
	  newDigit = (AliHLTPHOSDigit*)newDigArray->New(j);
	  newDigit->SetX(digit->GetX());
	  newDigit->SetZ(digit->GetZ());
	  newDigit->SetAmplitude(digit->GetAmplitude());
	  newDigit->SetTime(digit->GetTime());
	  newDigit->SetEnergy(digit->GetEnergy());
	  newDigit->SetGain(digit->GetGain());
	  newDigit->SetRawData((Int_t*)(digit->GetRawData()));
	  if(digit->GetCrazyness() == 0)
	    newDigit->SetCrazyness(1);
	  if(digit->GetCrazyness() == 1)
	    newDigit->SetCrazyness(0);
	  newDigit->SetBaseline(digit->GetBaseline());
	}
      newTree->Fill();
      newDigArray->Clear();
    }
  cout << "Conversion done, writing file: " << filename << endl;;
  TFile *outfile = new TFile(filename,"recreate");
  newTree->Write();
  delete outfile;
  outfile = 0;
  cout << "Done!\n";
  return 0;
}
  
 
