void makeLaserDAOutput(Int_t run, Int_t nruns)
{
  TH1F* hAmp = new TH1F("hAmpLaser"," Laser amplitude ", 2400, 0.5, 40.5);  
  /* 
 Mips RUN 0.5 176418 0.6 176419 0.7 176420 0.8 176421 0.9 176422 1 176423 1.2 176424 1.4 176425 1.6 176426 1.8 176427 2.0 176428  3 176429 3.5 XXXX 4 176430 5 176431 6 176432 7 176433 8 176434 9 176435 10 176436  15 176437 20 176438 25 176439 
  */
  /*
  Int_t runs[27] = {175727, 175729, 175730, 175731, 175733,
		    175734, 175735, 175736, 175737, 175738,
		    175739, 175740, 175741, 175742, 175743,
		    175744, 175745, 175746, 175747, 175748,
		    175749, 175750, 175751, 175752, 175753 };
  */     
  Float_t mips[27]= { 0.5, 0.6, 0.7, 0.8, 0.9,
		      1,   1.2, 1.4, 1.6, 1.8,
		      2,   3,   4,  
		      5,   6,   7,    8,  9, 
		      10, 15,   20,  25,  30 }; 

  Float_t meanEstimate = 0;
 Int_t ik=0;
 TString histname1, histname2, histname3;
 TString buf1, buf2, buf3;
 TFile * hist = new TFile("DAoutput.root", "RECREATE");
 for (Int_t ii=0; ii<nruns; ii++)
   {
     hAmp->Fill(mips[ii]);
     cout<<mips[ii]<<endl; 
     TString filename = Form("caliblaser/histCalib%i.root",run+ii);
     TFile* f = new TFile(filename.Data());
     cout<<filename<<endl;
     for (Int_t iii=0; iii<24; iii++)
       {
	 histname1 = Form("hCFD%i_%i",iii+1,ii+1);
	 TH1F *  cfd = new TH1F(histname1.Data(),histname1.Data(),1000,2000,3000);
	 histname2 = Form("hQTC%i_%i",iii+1,ii+1);
	 TH1F *  qtc = new TH1F(histname2.Data(),histname2.Data(),2500,0,10000);
	 histname3 = Form("hLED%i_%i",iii+1,ii+1);
	 TH1F *  led = new TH1F(histname3.Data(),histname3.Data(),500,0, 1000);
	 
	 f->cd();
	 //	 f->ls();
	 if(iii<12 ) buf1=(Form("T0_C_%i_CFD",iii+1));
	 if(iii>11 ) buf1=(Form("T0_A_%i_CFD",iii+1-12));

	 TH1F * cfdfile = (TH1F*)f->Get(buf1);
	 TH1F * qtcfile =(TH1F*)f->Get(Form("QTC%i",iii+1));
	 TH1F * ledfile =(TH1F*)f->Get(Form("LEDminCFD%i",iii+1));
	 cfdfile->SetDirectory(0);
	 qtcfile->SetDirectory(0);
	 ledfile->SetDirectory(0);
	 cfd->Add(cfdfile,cfd,1,1);
	 qtc->Add(qtcfile,qtc,1,1);
	 led->Add(ledfile,led,1,1);
	 
	 if (iii == 0) {
	   meanestimate = cfd->GetBinCenter(cfd->GetMaximumBin());
	 }
	 cfd->GetXaxis()->SetRangeUser( meanestimate-500, meanestimate+500);
	 hist->cd();
	 cfd->Write();
	 led->Write();
	 qtc->Write();
	 cfdfile->Delete();
	 qtcfile->Delete();
	 ledfile->Delete();
	}
     f->Close();


    }
 hist->cd();
 hAmp->Write();

}
