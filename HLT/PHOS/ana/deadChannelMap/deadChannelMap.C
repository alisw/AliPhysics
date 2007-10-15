
Int_t deadChannelMap(const char* runNb)
{
  char filepath [50];
  gROOT->ProcessLine(".L ../digits/AliHLTPHOSAltroConfig.cxx++");
  gROOT->ProcessLine(".L ../digits/AliHLTPHOSConfig.cxx++");
  gROOT->ProcessLine(".L ../digits/AliHLTPHOSBase.cxx++");
  gROOT->ProcessLine(".L ../digits/AliHLTPHOSDigit.cxx++");
  TChain *tree= new TChain("digitTree");

  sprintf(filepath, "/tmp/phoshlt/aldaqpc019/hlt/run%s*", runNb);
  tree->Add(filepath);
 
  TClonesArray *digArray = new TClonesArray("AliHLTPHOSDigit" , 100);
  //tree->SetBranchAddress("digits", &digArray);
  tree->SetBranchAddress("Digit", &digArray);
  TH2D *deadMap = new TH2D("deadMap", "Dead Channel Map", 64, 0, 63, 56, 0, 55);
  AliHLTPHOSDigit *digit = 0;

  cout << endl << "Chain with " << tree->GetEntries() << " events...\n";
  for(int k = 0; k < tree->GetEntries(); k++)
    {
      tree->GetEntry(k);
      //  cout << "     Event with " << digArray->GetEntriesFast() << " digits\n";
      for(int j = 0; j < digArray->GetEntriesFast(); j++)
	{
	  digit = (AliHLTPHOSDigit*)digArray->At(j);
          if(digit->GetGain() == 1)
	  {
	  for(int i = 0; i < 70; i++)
	    {
	      if((digit->GetRawData())[i] > 300)
		{
		  if((digit->GetRawData())[i+1] > 300)
		    {
		      if((digit->GetRawData())[i+2] > 300)
			{
			  if((digit->GetRawData())[i+3] > 300)
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
    }

  
  deadMap->Draw();
  
  return 0;
}
  
  
