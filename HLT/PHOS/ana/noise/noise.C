
Int_t noise(int runNbHV, int runNbNoHV)
{
  char filepath [50];

  sprintf(filepath, "/opt/HLT-public/rundir/histo_outdata/run%d_RMSHistogram.root", runNbHV);
  TFile *infile = TFile::Open(filepath);
  
  TH2F *RMSHVHist = infile->Get("RMSHGMapHist");
  //infile->Close();
  
  sprintf(filepath, "/opt/HLT-public/rundir/histo_outdata/run%d_RMSHistogram.root", runNbNoHV);
  TFile *infile2 = TFile::Open(filepath);
  TH2F *RMSNoHVHist = infile2->Get("RMSHGMapHist");
  //  infile->Close();

  TH2F *DiffHist = new TH2F("diffHist", "Map", 64, 0 , 63, 56, 0, 55);
  
  for(Int_t x = 0; x < 64; x++)
    {
      for(Int_t z = 0; z < 56; z++)
	{
	  if(RMSNoHVHist->GetBinContent(x, z) > 3)
	    {
	      if((RMSNoHVHist->GetBinContent(x, z) - RMSHVHist->GetBinContent(x, z)) < 1 && (RMSNoHVHist->GetBinContent(x, z) - RMSHVHist->GetBinContent(x, z)) > -1)
		{
		  DiffHist->SetBinContent(x, z, 10);
		}
	    }
	}
    }
  
  TFile *outfile = TFile::Open("hist.root","RECREATE");
  DiffHist->Write();
  outfile->Close();
  //  TCanvas *c1 = new TCanvas("c1");
  //  RMSHVHist->Draw();
  //TCanvas *c2 = new TCanvas("c2");
  // RMSNoHVHist->Draw();
  // TCanvas *c3 = new TCanvas("c3");
  DiffHist->Draw();

  
  return 0;
}
  
  
