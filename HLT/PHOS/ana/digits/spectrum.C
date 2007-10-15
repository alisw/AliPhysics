
Int_t spectrum(const char* filename)
{
  gROOT->ProcessLine(".L AliHLTPHOSDebugRawDigit.cxx++");
  
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
  
  
