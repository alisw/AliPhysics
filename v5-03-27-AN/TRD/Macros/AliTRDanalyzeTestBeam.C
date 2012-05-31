
void AliTRDanalyzeTestBeam(Int_t run, Int_t begin, Int_t end) {
  
  gROOT->SetStyle("Plain");
  //gStyle->SetPadTopMargin(0.02);
  //gStyle->SetPadRightMargin(0.02);
  //gStyle->SetPadLeftMargin(0.1);

  gStyle->SetPalette(1);
  gStyle->SetOptStat(0);
  gStyle->SetOptDate();
  
  TGaxis::SetMaxDigits(3);

  TStopwatch st;
  st.Start();

  const Int_t N = 640;

  // declare histograms
  const int nbins = 100;
  const int nstart = 0;

  TH1D *mSi1L = new TH1D("mSi1L", ";number of Si1 fired pads", 50, nstart-0.5, nstart+nbins-0.5);
  TH1D *mSi2L = new TH1D("mS21L", ";number of Si2 fired pads", 50, nstart-0.5, nstart+nbins-0.5);
 
  TProfile *mSi1ChP = new TProfile("mSi1ChP",";Si1 pad number;signal", 1280, -0.5, 1279.5, 0, 200, "s");
  TProfile *mSi2ChP = new TProfile("mS21ChP",";Si2 pad number;signal", 1280, -0.5, 1279.5, 0, 200, "s");
  
  TH1D *mSi1N = new TH1D("mSi1N", "Noise Dist Si1;ADC", 100, 0, 50);
  TH1D *mSi2N = new TH1D("mSi2N", "Noise Dist Si2;ADC", 100, 0, 50);
  
  TH1D *mSiCh[4];
  mSiCh[0] = new TH1D("mSi1ChX", ";Si1X pad amplitude (ADC)", 250, -0.5, 249.5);
  mSiCh[1] = new TH1D("mSi1ChY", ";Si1Y pad amplitude (ADC)", 250, -0.5, 249.5);
  mSiCh[2] = new TH1D("mSi2ChX", ";Si2X pad amplitude (ADC)", 250, -0.5, 249.5);
  mSiCh[3] = new TH1D("mSi2ChY", ";Si2Y pad amplitude (ADC)", 250, -0.5, 249.5);
  
  TH1D *mSiFullCh[4];
  mSiFullCh[0] = new TH1D("mSi1fChX", "Si1X;max amplitude (ADC)", 300, -0.5, 299.5);
  mSiFullCh[1] = new TH1D("mSi1fChY", "Si1Y;max amplitude (ADC)", 300, -0.5, 299.5);
  mSiFullCh[2] = new TH1D("mSi2fChX", "Si2X;max amplitude (ADC)", 300, -0.5, 299.5);
  mSiFullCh[3] = new TH1D("mSi2fChY", "Si2Y;max amplitude (ADC)", 300, -0.5, 299.5);  

  TH2D *mPos[2];
  mPos[0] = new TH2D("posSi1", ";x Si1 (mm);y Si1 (mm)", 128, 0, 32, 128, 0, 32);
  mPos[1] = new TH2D("posSi2", ";x Si2 (mm);y Si2 (mm)", 128, 0, 32, 128, 0, 32);

  TH2D *mCor[2];
  mCor[0] = new TH2D("corX", ";Si1 X (mm);Si2 X (mm)", 128, 0, 32, 128, 0, 32);
  mCor[1] = new TH2D("corY", ";Si1 Y (mm);Si2 Y (mm)", 128, 0, 32, 128, 0, 32);

  TH2D *mChCor[2];
  mChCor[0] = new TH2D("ChCorSi1", ";Si1 amp X;Si1 amp Y", 100, 0, 200, 100, 0, 200);
  mChCor[1] = new TH2D("ChCorSi2", ";Si2 amp X;Si2 amp Y", 100, 0, 200, 100, 0, 200);
  
  gStyle->SetOptStat(11);
  TH1D *mPb   = new TH1D("mPb", ";Amp Pb (ADC)", 150, -0.5, 4499.5);
  TH1D *mCher = new TH1D("mCher", ";Amp Cherenkov (ADC)", 150, -0.5, 4499.5);
  TH2D *mPbCher = new TH2D("mPbCher", ";amp Cherenkov;amp Pb", 150, -0.5, 4499.5, 150, -0.5, 4599.5);
  // gStyle->SetOptStat(0);

  // needed by the AliTRDRawStreamTB
  AliCDBManager::Instance()->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  AliTRDcalibDB *calib = AliTRDcalibDB::Instance();
  calib->SetRun(0);

  // TRD data monitoring
  TH1D *mDet = new TH1D("mDet", ";chamber", 20, -0.5, 19.5);
  TH1D *mROB = new TH1D("mROB", ";rob", 20, -0.5, 19.5);
  TH1D *mTRDsig = new TH1D("mTRDsig", ";trdSignal", 100, -0.5, 299.5);
  
  //AliLog::SetClassDebugLevel("AliTRDRawStreamTB", 10);

  int counter = 0;
  // for(Int_t run = 365; run < 386; run++) {
  //for(Int_t run = 369; run < 382; run++) {
  //for(int run = 387; run < 389; run++) {

  // if (run == 389) continue;
    cout << run << endl;

    for(Int_t fn=begin; fn<end; fn++) {
      
      // connect to data 
      //const char *base="/Users/radomski/data/1GeV/";
      const char *base = "./";
      const char *filename = Form("%s/run%d_gdc_daq09.%03d.raw",base,run,fn);
      if (gSystem->AccessPathName(filename)) continue;
      cout << filename << endl;

      AliTRDtestBeam *data = new AliTRDtestBeam(filename);
      
      // process data
      while (data->NextEvent()) {
      
	if (!(counter%1000)) cout << "Event = " << counter << endl;
	counter++;
	
	/*
	AliTRDRawStreamTB *tb = data->GetTRDrawStream();
	while(tb->Next()) {
	  mROB->Fill(tb->GetROB());
	  mDet->Fill(tb->GetDet());
	  int *sig = tb->GetSignals();
	  mTRDsig->Fill(sig[0]);
	  mTRDsig->Fill(sig[1]);
	  mTRDsig->Fill(sig[2]);
	}
	delete tb;
	*/

	mCher->Fill(data->GetCher());
	mPb->Fill(data->GetPb());
	mPbCher->Fill(data->GetCher(), data->GetPb());

	mSi1L->Fill(data->GetNSi1());
	mSi2L->Fill(data->GetNSi2());
	
	for(int i=0; i<data->GetNSi1(); i++) {
	  Int_t q = data->GetSi1Charge(i);
	  Int_t a = data->GetSi1Address(i);
	  if (a == 0) continue; // noisy channels
	  mSi1ChP->Fill(a, q);
	  if (a < N) mSiCh[0]->Fill(q);
	  else mSiCh[1]->Fill(q);	       
	}
	
	for(int i=0; i<data->GetNSi2(); i++) {
	  Int_t q = data->GetSi2Charge(i);
	  Int_t a = data->GetSi2Address(i);
	  if (a == 0 || a == 1279) continue; // noisy channels
	  mSi2ChP->Fill(a, q);
	  if (a < N) mSiCh[2]->Fill(q);
	  else mSiCh[3]->Fill(q);	       
	}
	
	mSiFullCh[0]->Fill(data->GetQx(0));
	mSiFullCh[1]->Fill(data->GetQy(0));
	mSiFullCh[2]->Fill(data->GetQx(1));
	mSiFullCh[3]->Fill(data->GetQy(1));
	
	for(int k=0; k<2; k++)
	  mChCor[k]->Fill(data->GetQx(k), data->GetQy(k));
	
	/*
	  if (data->GetQx(0) < 20) continue;
	  if (data->GetQx(1) < 20) continue;
	  if (data->GetQy(0) < 20) continue;
	  if (data->GetQy(1) < 20) continue;
	*/
	
	for(int k=0; k<2; k++)
	  mPos[k]->Fill(data->GetX(k), data->GetY(k));
      
	mCor[0]->Fill(data->GetX(0), data->GetX(1));
	mCor[1]->Fill(data->GetY(0), data->GetY(1));
      }
      
      delete data;
    }
    //}

  // process histograms
  for(int i=1; i<1281; i++) mSi1N->Fill(mSi1ChP->GetBinError(i));
  for(int i=1; i<1281; i++) mSi2N->Fill(mSi2ChP->GetBinError(i));

  // display
  cout << "Number of Events = " << counter << endl;
  
  /**/
  TCanvas *c = new TCanvas("siliconSignal", "silicon signal");
  c->Divide(2,2, 0.01, 0.01);
  c->cd(1);
  mSi1L->Draw();
  c->cd(2);
  mSi2L->Draw();
  c->cd(3);
  mSi1ChP->Draw();
  c->cd(4);
  mSi2ChP->Draw();
  /* */

  /**/
  c = new TCanvas();
  c->Divide(2,2,0.01,0.01);
  c->cd(1);
  mSi1N->Draw();
  c->cd(2);
  mSi2N->Draw();
  /**/

  /**/
  // pads
  c = new TCanvas("siPads", "Silicon Pads");
  c->Divide(2,2, 0.01, 0.01);
  for(int i=0; i<4; i++) {
    c->cd(i+1);
    gPad->SetLogy();
    mSiCh[i]->Draw();
  }
  
  // clusters
  c = new TCanvas("siCluster", "silicon clusters");
  c->Divide(2,2, 0.01, 0.01);
  for(int i=0; i<4; i++) {
    c->cd(i+1);
    gPad->SetLogy();
    mSiFullCh[i]->Draw();
  }

  // position and correlation
  c = new TCanvas("siPosition", "reconstructed position");
  c->Divide(2,2, 0.01, 0.01);
  for(int i=0; i<2; i++) {
    c->cd(1+i);
    mPos[i]->Draw("col");
    c->cd(3+i);
    mCor[i]->Draw("col");  
  }
  
  c = new TCanvas("siCharge", "si charge correlation");
  c->Divide(2,2, 0.01, 0.01);
  c->cd(1);
  gPad->SetLogz();
  mChCor[0]->Draw("col");
  c->cd(2);
  gPad->SetLogz();
  mChCor[1]->Draw("col");
  /**/

  new TCanvas();
  gPad->SetLogy();
  mCher->Draw();
  
  // electron sample
  int bin = mCher->FindBin(500.);
  double ef = (mCher->Integral(bin, 151)/ mCher->GetSum());
  TLatex *l = new TLatex(2e3, 0.02*mCher->GetSum(), Form("Electron fraction = %.2f ", ef));
  l->Draw();
  
  new TCanvas();
  gPad->SetLogy();
  mPb->Draw();

  new TCanvas();
  gPad->SetLogz();
  mPbCher->Draw("colz");
  
  /*
  c = new TCanvas();
  c->Divide(2,2,0.01, 0.01);
  c->cd(1);
  gPad->SetLogy();
  mTRDsig->Draw();

  c->cd(2);
  mROB->Draw();

  c->cd(3);
  mDet->Draw();
  */
  
  st.Stop();
  st.Print();
}

//Int_t SelectEvent() {
//}
