/*
  a very basic test macro for processing Pedestal events with AliCaloCalibPedestal
  (adopted from EMCALLEDda.C macro)
  Josh Hamblen
*/
#define AliDebugLevel() -1

void EMCALPEDda(const int runno = 476){
  
  int i, status;
  
  /* log start of process */
  printf("EMCAL DA started - %s\n",__FILE__);
  
  AliCaloCalibPedestal * calibPedestal = new 
    AliCaloCalibPedestal(AliCaloCalibPedestal::kEmCal); // pedestal and noise calibration
  
  int nevents=0;
  // Assume we are just interested in the 1st segment, _0.root below for fname*
  Char_t fname[256];
  sprintf(fname, "/local/data/Run_%09d.Seq_1A.Stream_0.root",runno);
  
  AliRawReader *rawReader = NULL;
  rawReader = new AliRawReaderRoot(fname);
  AliCaloRawStream *in = NULL; 
  in = new AliCaloRawStream(rawReader,"EMCAL");
  //in->SetOldRCUFormat(kTRUE);
  
  AliRawEventHeaderBase *aliHeader=NULL;
  int nev=0;
  /* read until EOF */
  while ( rawReader->NextEvent()) {
    nev++;
    
    aliHeader = (AliRawEventHeaderBase*) rawReader->GetEventHeader();
    calibPedestal->SetRunNumber( aliHeader->Get("RunNb") ); // just for fun; keep info on last run looked at
    
    // select physics and calibration events now (only calibration in future)
    if ( aliHeader->Get("Type") == AliRawEventHeaderBase::kPhysicsEvent || 
	 aliHeader->Get("Type") == AliRawEventHeaderBase::kCalibrationEvent  ) {
      
      nevents++;
      if(nevents%1000==0)cout<<"Event "<<nevents<<endl;
      //  Pedestal calibration
      calibPedestal->ProcessEvent(in);
    }
  } // loop over all events in file
    /* cleanup the reading handles */
  delete in;
  delete rawReader;    
  
  //
  // write results/histograms to rootfile
  //
  
  printf ("%d physics/calibration events processed.\n",nevents);
  
  // create output histograms and write to file
  Char_t outname[256];
  sprintf(outname, "PED_%09d.root",runno);
  
  calibPedestal->SaveHistograms(outname);

  printf("Wrote %s.\n",outname);
    
  // see if we can delete our analysis helper also
  delete calibPedestal;
  
  //
  // reopen the file that was just made to make the RMS hists
  //
  
  TFile *f=new TFile(outname,"update");
  
  // now take the profile histograms and plot the RMS from them 
  int ncols = hPedlowgain0->GetNbinsX();
  int nrows = hPedlowgain0->GetNbinsY();
  Double_t rms;
  
  TH2F *hPedRMSlowgain0 = new TH2F("hPedRMSlowgain0","Pedestal RMS, low gain, module 0",ncols,0,ncols,nrows,0,nrows); 
  TH2F *hPedRMShighgain0 = new TH2F("hPedRMShighgain0","Pedestal RMS, high gain, module 0",ncols,0,ncols,nrows,0,nrows);   
  TH2F *hSampleRMSlowgain0 = new TH2F("hSampleRMSlowgain0","All Samples RMS, low gain, module 0",ncols,0,ncols,nrows,0,nrows); 
  TH2F *hSampleRMShighgain0 = new TH2F("hSampleRMShighgain0","All Samples RMS, high gain, module 0",ncols,0,ncols,nrows,0,nrows);  

  for (int i=0; i < ncols; i++){
    for (int j=0; j < nrows; j++){
      
      rms = hPedlowgain0->GetBinError(i+1,j+1);
      hPedRMSlowgain0->Fill(i,j,rms);
      
      rms = hPedhighgain0->GetBinError(i+1,j+1);
      hPedRMShighgain0->Fill(i,j,rms);
      
      rms = hSamplelowgain0->GetBinError(i+1,j+1);
      hSampleRMSlowgain0->Fill(i,j,rms);
      
      rms = hSamplehighgain0->GetBinError(i+1,j+1);
      hSampleRMShighgain0->Fill(i,j,rms);
   }
  }
   
  // write the RMS hists to the file
  hPedRMSlowgain0->Write();
  hPedRMShighgain0->Write();
  hSampleRMSlowgain0->Write();
  hSampleRMShighgain0->Write();
  f->ls();

  // now draw the results
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);
  
  TCanvas *c1 = new TCanvas("c1","",800,600);
  c1->Divide(2,2);
  c1->cd(1);
  hPedlowgain0->Draw("colz");
  c1->cd(2);
  hPedRMSlowgain0->Draw("colz");
  c1->cd(3);
  hPedhighgain0->Draw("colz");
  c1->cd(4);
  hPedRMShighgain0->Draw("colz");
  c1->cd();

  TCanvas *c2 = new TCanvas("c2","",800,600);
  c2->Divide(2,2);
  c2->cd(1);
  hSamplelowgain0->Draw("colz");
  c2->cd(2);
  hSampleRMSlowgain0->Draw("colz");
  c2->cd(3);
  hSamplehighgain0->Draw("colz");
  c2->cd(4);
  hSampleRMShighgain0->Draw("colz");
  c2->cd();


}
