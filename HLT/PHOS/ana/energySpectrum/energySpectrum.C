
Int_t energySpectrum(const char* runNb)
{
  gROOT->ProcessLine(".L ../digits/AliHLTPHOSAltroConfig.cxx++");
  gROOT->ProcessLine(".L ../digits/AliHLTPHOSDigit.cxx++");
  TH1F *spectrumHist = new TH1F("spectrumHist", "Digit energy spectrum", 200, 0, 200);
  TH1F *spectrumHistSingle = new TH1F("spectrumHist", "Digit energy spectrum", 200, 0, 200);
  TClonesArray *digArray = new TClonesArray("AliHLTPHOSDebugRawDigit" , 100);
  TChain *tree= new TChain("digitTree");
  AliHLTPHOSDebugRawDigit *digit = 0;

  char filepath [50];
  sprintf(filepath, "/tmp/phoshlt/analysis/data/run%s/*", runNb);
  tree->Add(filepath);
  
  tree->SetBranchAddress("DebugRawDigit", &digArray);
  cout << endl << "Entries in tree: " << tree->GetEntries() << endl;

  for(int k = 0; k < tree->GetEntries(); k++)
    {
      tree->GetEntry(k);
      for(int j = 0; j < digArray->GetEntriesFast(); j++)
	{
	  digit = (AliHLTPHOSDebugRawDigit*)digArray->At(j);
	  if(digit->GetAmplitude() > 20&& digit->GetCrazyness()==1 && digit->GetGain() == 1)
	    {
	    spectrumHist->Fill(digit->GetAmplitude());
	  if(digit->fX = 30 &&digit->fZ = 30)
	    {
	      spectrumHistSingle->Fill(digit->GetAmplitude());
	    }
	    }
	}
    }
  spectrumHist->Draw();

  return 0;
}
  
  
