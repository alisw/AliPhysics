///////////////////////////////////////////////
// Simple event display for channel spectra
// Input file should be raw DATE root tree file
///////////////////////////////////////////////

// for the readout
const int NFEC    =  36;                        // Max 36 FrontEndCard per SM
const int NCHIP   =  4;                         // 4 ALTROs per FEC
const int NCHAN   = 16;                         // Channels per ALTRO
const int TOTCHAN = NFEC * NCHIP * NCHAN;	// Max uses ALTRO channels

const int NCOLS = 48; // per SM
const int NROWS = 24;
const int NMINI = 3; // number of ministrips, or T-cards, per strip

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
void RAW_evtdis(const int runno = 615,
		const int gainv = 1,  /*0=low, 1=high*/
		const int evtnum= 1,
		const int strip = 0,
		int ymax=200, // set the scale of plots
		const int delay = 1)  // -1=no delay, wait for input, X>=0 => sleep aprox. X sec. after making plot
{
  // set ranges to plot
  const int col_f = strip*2;
  const int col_l = col_f + 1;
  const int row_f = 0;
  const int row_l = NROWS - 1;

  const int nsamples = 65; // number of ADC time samples per channel and event

  const int saveplot = 0;
  const int numbering      = 1; // 0: no numbering, 1: nubering on each plot
  const int dofit = 0; // 0: no fit, 1: try to fit the spectra (not debugged) 
  const int debug    = 0;
  const float gammaN = 2;
  // end of setup    

  // Assume we are just interested in the 1st segment, _0.root below for fname*
  Char_t fname[256];
  sprintf(fname, "/local/data/Run_%09d.Seq_1A.Stream_0.root",runno);

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

    //	int idx = arow * NCOLS + acol + NCOLS*NROWS * gain; // encoding used later
    int gain = i / (NCOLS*NROWS);
    int row = (i % (NCOLS*NROWS))/NCOLS;
    int col = i % NCOLS;
    sprintf(ch_label[i], "Col%02d Row%02d", col, row);
  }

  int numcol = col_l-col_f+1;
  int numrow = row_l-row_f+1;

  TCanvas *cc[NMINI];
  for (int ic=0; ic<NMINI; ic++) {
    sprintf(buff1, "cc%d", ic);
    // use the usual StripX CY, X=0..23, Y=0..2 notation for T-cards (mini-strips)
    sprintf(name,"Strip %d MiniStrip(T-card) C%d", strip, ic);
    cc[ic] = new TCanvas(buff1,name, 400*ic, 10, 400,800);
    cc[ic]->Divide(numcol, numrow/NMINI);
    cc[ic]->SetTitle(name);
  }

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

	int acol = in->GetColumn();
	int arow = in->GetRow();
	int gain = in->GetCaloFlag();

	int idx = arow * NCOLS + acol + NCOLS*NROWS * gain;
	//cout << "hist idx " << idx << endl;
	
	if (idx < 0 || idx > TOTCHAN) { 
	  cout << "Hist idx out of range: " << idx << endl;
	}
	else { // reasonable range of idx (removes TRU and LEDMon data also)
	  hfit[idx]->SetBinContent(in->GetTime(), in->GetSignal());
	}

      } // Raw data read

      // Next: let's actually plot the data..
      int nrow_mini = NROWS/NMINI; // number of rows per T-card or mini-strip
      for (Int_t arow = row_f; arow <= row_l; arow++) {
	for (Int_t acol = col_f; acol <=  col_l; acol++) {
	  int idx = arow * NCOLS + acol + NCOLS*NROWS * gainv;

	  // which T-card does the tower belong to?
	  int mini = arow / (nrow_mini); 	  
	  // on which pad should we plot it?
	  int pad_id = ((row_l-arow)%nrow_mini) *(col_l-col_f+1) + acol - col_f + 1;
	  
	  cout << "row="<< arow << ", col=" << acol
	       << ", C" << mini << ", pad=" << pad_id << endl;
	  cc[mini]->cd(pad_id);
	  hfit[idx]->SetTitle("");
	  hfit[idx]->SetFillColor(5);
	  hfit[idx]->SetMaximum(ymax);
	  hfit[idx]->SetMinimum(0);
	  // we may or may not decide to fit the data
	  if (dofit) {
	    hfit[idx]->Fit(f1[i]);
	  }
	  hfit[idx]->Draw();
	  if( numbering ) {
	    t->SetTextColor(clr[gainv]);
	    t->DrawTextNDC(0.45,0.45,ch_label[idx]);
	  }
	}
      }

      // add some extra text on the canvases
      for (int ic=0; ic<NMINI; ic++) {
	// print a box showing run #, evt #, and timestamp
	cc[ic]->cd();
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
	cc[ic]->Update();
	cout << "Done" << endl;
      }

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

    } // event selection
    iev ++;
  } // event loop

  // save plot, if setup/requested to do so
  char plotname[100];
  if (saveplot==1) {
    for (int ic=0; ic<NMINI; ic++) {
      sprintf(plotname,"Run_%d_Ev%d_Strip%d_C%d_Gain%d_MaxHist%d.gif",
	      runno,iev,strip,ic,gainv,ymax);  
  
      cout <<"SAVING plot:"<< plotname << endl;
      cc[ic]->SaveAs(plotname);
    }
  }

}
