///////////////////////////////////////////////
// Simple event display for channel spectra
// Input file should be raw DATE root tree file
///////////////////////////////////////////////

// for the readout: 1 FEE reads out LED ref. info for the whole SuperModule
// 1 CSP per StripModule
const int NSTRIPS = 24; // number of StripModules in SuperModule
const int TOTCHAN = NSTRIPS * 2;	// *2 since we have high gain and low gain
// we group the calibrations in 3 sets of 8 strips
const int NSETS = 3; 
const int NSTRIPS_IN_SET = 8;

// gamma2 fit function
double fitfun(double *x, double *par) {
  double Amp	= par[0];
  double Tmax	= par[1];
  double Tau	= par[2];
  double Ped	= par[3];
  double gammaN   = par[4];
  double t = 0;
  if(Tau) t = (x[0] - Tmax + Tau)/Tau;
  if(t<0) t = 0;
  
  // Since we have now set gammaN to 2, we replace the pow(t, 2) call in
  // double f = Amp * pow(t,gammaN) * exp(gammaN*(1-t)) + Ped;
  // with just t*t
  double f = Amp * t*t * exp(gammaN*(1-t)) + Ped;
  return f;
}

// main method
void LEDRef_evtdis(const int runno = 615,
		   const int gainv = 0,  /*0=low, 1=high*/
		   const int evtnum= -10,
		   int ymax=1023, // set the scale of plots
		   const int delay = 1)  // -1=no delay, wait for input, X>=0 => sleep aprox. X sec. after making plot
{
  // set ranges to plot
  const int strip_f = 0; // first
  const int strip_l = NSTRIPS - 1;
  
  const int nsamples = 65; // number of ADC time samples per channel and event
  
  const int saveplot = 0;
  const int numbering      = 1; // 0: no numbering, 1: nubering on each plot
  const int dofit = 0; // 0: no fit, 1: try to fit the spectra 
  const int debug    = 0;
  const float gammaN = 2;
  // end of setup    
  
  // Assume we are just interested in the 1st segment, _0.root below for fname*
  Char_t fname[256];
  sprintf(fname, "/local/data/Run_%09d.Seq_1A.Stream_0.root",runno);
  cout << "TOTCHAN " << TOTCHAN << endl;

  // set up a raw reader of the data
  AliRawReader *rawReader = NULL;
  rawReader = new AliRawReaderRoot(fname);
  AliCaloRawStream *in = NULL; 
  in = new AliCaloRawStream(rawReader,"EMCAL");

  // set up histograms
  TH1F *hfit[TOTCHAN];
  TF1 *f1[TOTCHAN];
  char ch_label[TOTCHAN][100];
  char buff1[100];
  char name[80];
  for(int i=0; i<TOTCHAN; i++) {
    sprintf(buff1,"hfit_%d",i);
    hfit[i] = new TH1F(buff1,"hfit", nsamples , -0.5, nsamples - 0.5);
    hfit[i]->SetDirectory(0);
    sprintf(name,"f1_%d",i);
    f1[i] = new TF1(name,fitfun,0,70,5);
    f1[i]->SetLineWidth(2);
    f1[i]->SetLineColor(2);

    //	int idx = istrip + NSTRIPS * gain; // encoding used later
    int gain = i / (NSTRIPS);
    int istrip = i % NSTRIPS;
    sprintf(ch_label[i], "Strip%02d", istrip);
  }
  
  TCanvas *cc1 = new TCanvas("cc1","3 columns of 8 strips each",600,800);
  int numcol = NSETS;
  int numrow = NSTRIPS_IN_SET;
  cc1->Divide(numcol, numrow);
  
  TText *t = new TText;
  t->SetTextSize(0.17);
  int clr[2] = {4,2}; // colors
  
  // figure out which events we should look at
  int firstevent = evtnum;
  int lastevent = evtnum;
  if (evtnum < 0) { // get a bunch of events
    firstevent = 0;
    lastevent = - evtnum;
  }
  if (evtnum == 0) { // get all events
    firstevent = 0;
    lastevent = 1000000;
  }
  
  Int_t iev =0;
  AliRawEventHeaderBase *aliHeader=NULL;    
  while ( rawReader->NextEvent() && iev < firstevent) {
    aliHeader = (AliRawEventHeaderBase*) rawReader->GetEventHeader();
    iev++;
  }
  
  // loop over selected events
  while ( rawReader->NextEvent() && iev <= lastevent) {
    aliHeader = (AliRawEventHeaderBase*) rawReader->GetEventHeader();
    int runNumber = aliHeader->Get("RunNb"); 
    
    cout << "Found run number " << runNumber << endl;
    
    // reset histograms
    for(int i=0; i<TOTCHAN; i++) {
      hfit[i]->Reset();
    }
    
    // get events (the "1" ensures that we actually select all events for now)
    if ( 1 || aliHeader->Get("Type") == AliRawEventHeaderBase::kPhysicsEvent ) {
      const UInt_t * evtId = aliHeader->GetP("Id");
      int evno_raw = (int) evtId[0];
      int timestamp = aliHeader->Get("Timestamp");
      
      cout << " evno " << evno_raw
	   << " size " << aliHeader->GetEventSize()
	   << " type " << aliHeader->Get("Type")
	   << " type name " << aliHeader->GetTypeName()
	   << " timestamp " << timestamp
	   << endl;
      
      /// process_event stream
      while ( in->Next() ) {
	
	int strip = in->GetColumn();
	int gain = in->GetRow();
	
	if (in->IsLEDMonData()) {
	  
	  int idx = strip + NSTRIPS*gain;
	  //cout << "hist idx " << idx << endl;
	  
	  if (idx < 0 || idx > TOTCHAN) { 
	    cout << "Hist idx out of range: " << idx << endl;
	  }
	  else { // reasonable range of idx
	    hfit[idx]->SetBinContent(in->GetTime(), in->GetSignal());
	  }

	} // LED Ref data only

      } // Raw data read
    
      // Next: let's actually plot the data..
      for (Int_t strip = strip_f; strip <= strip_l; strip++) {
	
	int idx = strip + NSTRIPS*gainv;
	
	// which set/column does the strip belong in
	int iset = strip / NSTRIPS_IN_SET; 	  
	int within_set = strip % NSTRIPS_IN_SET; 	  
	// on which pad should we plot it?
	int pad_id = (NSTRIPS_IN_SET-1-within_set)*NSETS + iset + 1;
	
	cout << "strip " << strip 
	     << ". set="<< iset << ", within_set=" << within_set
	     << ", pad=" << pad_id << endl;
	cc1->cd(pad_id);
	hfit[idx]->SetTitle("");
	hfit[idx]->SetFillColor(5);
	hfit[idx]->SetMaximum(ymax);
	hfit[idx]->SetMinimum(0);
	// we may or may not decide to fit the data
	if (dofit) {
	  f1[i]->SetParameter(0, 0); // initial guess; zero amplitude :=)
	  hfit[idx]->Fit(f1[i]);
	}
	hfit[idx]->Draw();
	if( numbering ) {
	  t->SetTextColor(clr[gainv]);
	  t->DrawTextNDC(0.65,0.65,ch_label[idx]);
	}
      }

      // add some extra text on the canvas
      // print a box showing run #, evt #, and timestamp
      cc1->cd();
      // first draw transparent pad
      TPad *trans = new TPad("trans","",0,0,1,1);
      trans->SetFillStyle(4000);
      trans->Draw();
      trans->cd();
      // then draw text
      TPaveText *label = new TPaveText(.2,.11,.8,.14,"NDC"); 
      //  label->Clear();
      label->SetBorderSize(1);
      label->SetFillColor(0);
      label->SetLineColor(clr[gainv]);
      label->SetTextColor(clr[gainv]);
      //label->SetFillStyle(0);
      TDatime d;
      d.Set(timestamp);
      sprintf(name,"Run %d, Event %d, Hist Max %d, %s",runno,iev,ymax,d.AsString());
      label->AddText(name);
      label->Draw();
      cc1->Update();
      cout << "Done" << endl;
      
      // some shenanigans to hold the plotting, if requested
      if (firstevent != lastevent) {
	if (delay == -1) {
	  // wait for character input before proceeding
	  cout << " enter y to proceed " << endl;
	  char dummy[2];
	  cin >> dummy;
	  cout << " read " << dummy << endl;
	  if (strcmp(dummy, "y")==0) {
	    cout << " ok, continuing with event " << iev+1 << endl;
	  }
	  else {
	    cout << " ok, exiting " << endl;
	    //exit(1);
	  }
	}
	else {
	  cout << "Sleeping for " << delay * 500 << endl;
	  gSystem->Sleep(delay * 500);
	}
      }

      // save plot, if setup/requested to do so
      char plotname[100];
      if (saveplot==1) {
	sprintf(plotname,"Run_%d_LEDRef_Ev%d_Gain%d_MaxHist%d.gif",
		runno,iev,gainv,ymax);  
	cout <<"SAVING plot:"<< plotname << endl;
	cc1->SaveAs(plotname);
      }

    } // event selection

  iev ++;
  } // event loop

}
