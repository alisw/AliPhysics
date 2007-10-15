Int_t rawDisplay(const char* runNb)
{
  gROOT->ProcessLine(".L ../digits/AliHLTPHOSAltroConfig.cxx++");
  gROOT->ProcessLine(".L ../digits/AliHLTPHOSConfig.cxx++");
  gROOT->ProcessLine(".L ../digits/AliHLTPHOSBase.cxx++");
  gROOT->ProcessLine(".L ../digits/AliHLTPHOSDigit.cxx++");
  char filepath [50];
  char outfile [50];
  Int_t *dataPtr = 0;
  fUpperBound = 40;
  fLowerBound = 20;
  fStartTime = 10;
  fZeroThreshold = 5;
  AliHLTPHOSDigit *digitPtr = 0;

  int flag = 0;
  int count =0 ;
  TChain *tree= new TChain("digitTree");

  sprintf(filepath, "/tmp/phoshlt/aldaqpc019/hlt/run%s_digitTree_0.root", runNb);
  cout << "Opening file: " << filepath << endl;
  tree->Add(filepath);
 
  TClonesArray *digArray = new TClonesArray("AliHLTPHOSDigit" , 100);
  tree->SetBranchAddress("Digit", &digArray);
  TH1I *rawHist = new TH1I("rawHist", "Raw display", 70, 0, 69);
  
  TH1I *rawHist2 = new TH1I("rawHist2", "Raw display", 70, 0, 69);
  cout << endl << "Tree with " << tree->GetEntries() << " events...\n";
  
  for(int k = 0; k < tree->GetEntries(); k++)
    {
      tree->GetEntry(k);
      for(int j = 0; j < digArray->GetEntriesFast(); j++)
	{
	  digit = (AliHLTPHOSDigit*)digArray->At(j);
	  //if(digit->GetAmplitude() - digit->GetBaseline() > 200)
	    //	  if(digit->GetBaseline() && digit->GetCrazyness())
	  if(digit->GetCrazyness() == 1 && digit->GetAmplitude() > 200 && digit->GetRawData()[50] < 10 )
	    {
	      count++;
	      UInt_t *data = digit->GetRawData();
	      for(Int_t i = 0; i < 70; i++)
		{
		  rawHist->SetBinContent(i, (digit->GetRawData())[i]);
		}
	      if(count == 2)
		{
		  cout << "x: "<< digit->GetX() << " - z: " << digit->GetZ() << endl;
		  for(int c = 0; c < digArray->GetEntriesFast(); c++)
		    {
		      d = (AliHLTPHOSDigit*)digArray->At(c);
		      if(d->GetX() == digit->GetX() && d->GetZ() == digit->GetZ() && digit->GetGain() != d->GetGain())
			{
			   for(Int_t b = 0; b < 70; b++)
			     {
			       rawHist2->SetBinContent(b, (d->GetRawData())[b]);
			     }
			}
		    }
		  finished = true;
		  break;
		}
	    }
	}
      if(finished) break;
    }
  /*
  for(int k = 0; k < tree->GetEntries(); k++)
    { 
      if(flag)
	break;
      tree->GetEntry(k);
      for(int j = 0; j < digArray->GetEntriesFast(); j++)
	{
	  bool IsMIP = true;
	  digitPtr = (AliHLTPHOSDigit*)digArray->At(j);
	  dataPtr = digitPtr->GetRawData();
	  if(digitPtr->GetCrazyness() != 0)
	     {
	       // cout << "Crazy!\n";
	       continue;
	     }
	  if(digitPtr->GetAmplitude() < fLowerBound || digitPtr->GetAmplitude() > fUpperBound)
	     {
	       //  cout << "Wrong amplitude!\n";
	       continue;
	     }
	   
	  for(Int_t time = (Int_t)(digitPtr->GetTime() - 2); time < (digitPtr->GetTime() - 3); time++)
	    {
	      if((Float_t)dataPtr[time] < (digitPtr->GetAmplitude() - (digitPtr->GetAmplitude())/10))
		{
		  //  cout << "Not flat enough!\n";
		  IsMIP = false;
		  break;
		}
	    }
	  if(!IsMIP)
	    continue;
	  for(Int_t sample = 0; sample < fStartTime; sample++)
	    {
	      if(dataPtr[sample] > fZeroThreshold || dataPtr[sample] < -fZeroThreshold)
		{
		  //  cout << endl << dataPtr[sample] << endl;
		  //cout << "Not stable baseline!\n";
		  IsMIP = false;
		  break;
		}
	    }
	  if(dataPtr[(Int_t)fStartTime + 5] < fZeroThreshold)
	    {
	      //      cout << "Not stable signal!\n";
	      IsMIP = false;
	    }
	  if(IsMIP)
	    {
	      cout << "Got MIP!\n";
//	      fMIPCountEvent++;
	      // fChannelHistPtr->Fill(digitPtr->fX, digitPtr->fZ);
	   
	      for(int a = 0; a < 70; a++)
		{
		  rawHist->SetBinContent(a, (digitPtr->GetRawData())[a]);
		}
	      flag++;
	      break;
	    }
	}
	}*/
  
  rawHist->Draw();
  
  TCanvas *can = new TCanvas("can");
  rawHist2->Draw();

  return 0;
}
  
 
