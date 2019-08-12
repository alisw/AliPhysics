const int kNSM = 20; // max # of SuperModules
const int kNumSens = 160; // max total

const int kNumSensPerSM = 8; // per SuperModule (nominally, not true for EMCal 1/3 SMs, a.k.a. SMs 10,11 & 18,19 seemingly...)

// init to strange values
const int kInitMin = 100;
const int kInitMax = 0;
// limits for acceptable values
const int kTooHigh = 35;
const int kTooLow = 15;
// not too far from SuperModule median
const float kMaxDiffMedian = 2;

const int debug = 0;

void scanAll(const char *listfile = "filelist.txt")
{

  int nruns = 0;
  ifstream fin(listfile);
  Int_t runno = 0;
  char fn[200];

  Int_t iSM = 0;
  Int_t nOk = 0;
  Int_t nPTot = 0;
  Int_t minStart = 0;
  Int_t maxEnd = 0;
  Float_t aveMin = 0;
  Float_t aveMax = 0;
  // output
  TFile destFile("all.root", "recreate");
  destFile.cd();
  TTree *tree = new TTree("tree","");
  tree->Branch("runno", &runno, "runno/I");
  tree->Branch("iSM", &iSM, "iSM/I");
  tree->Branch("nOk", &nOk, "nOk/I");
  tree->Branch("nPTot", &nPTot, "nPTot/I");
  tree->Branch("minStart", &minStart, "minStart/I");
  tree->Branch("maxEnd", &maxEnd, "maxEnd/I");
  tree->Branch("aveMin", &aveMin, "aveMin/F");
  tree->Branch("aveMax", &aveMax, "aveMax/F");

  // loop over input
  while ( fin.good() ) {
    fin >> runno >> fn;
    if ( fin.good() ) {
      cout << " filename fn " << fn << " runno " << runno << endl;
      nruns++;

      // get the file
      TFile *file0 = TFile::Open(fn);
      AliEMCALSensorTempArray *arr = (AliEMCALSensorTempArray *) AliCDBEntry->GetObject();
      if (debug>0) {
	file0->ls();
	arr->Print();
	cout << " NumSensors " << arr->NumSensors()
	     << " GetFirstIdDCS() " << arr->GetFirstIdDCS()
	     << " GetLastIdDCS() " << arr->GetLastIdDCS()
	     << endl;
      }
      file0->Close();

      // info for each sensor
      int np[kNumSens] = {0};
      double min[kNumSens] = {0};
      double max[kNumSens] = {0};
      int id = 0;

      minStart = 0x7FFFffff;
      maxEnd = 0;

      for (int isensor=0; isensor<kNumSens; isensor++) {
	AliEMCALSensorTemp *o = arr->GetSensor(isensor);
	if (o) {
	  iSM = o->GetSide() + o->GetSector()*2;
	  id = o->GetIdDCS();
	  int is = iSM*kNumSensPerSM + o->GetNum();
	  if (is != id) {
	    cout << " id mismatch: id " << id << " " << is << endl;
	  }
	  if (debug>1) {
	    o->Print();
	    cout << " side " << o->GetSide()
		 << " sector " << o->GetSector()
		 << " num " << o->GetNum()
		 << " startTime " << o->GetStartTime()
		 << " endTime " << o->GetEndTime()
		 << endl;
	  }
	  if ( minStart > o->GetStartTime() ) { minStart = o->GetStartTime(); }
	  if ( maxEnd < o->GetEndTime() ) { maxEnd = o->GetEndTime(); }

	  AliSplineFit *f = o->GetFit();
	  np[is] = 0;
	  min[is] = kInitMin;
	  max[is] = kInitMax;
	  
	  if (f) {
	    np[is] = f->GetKnots();
	    if (debug>1) {
	      cout << " np " << np[is] << endl;
	    }
	    Double_t *x = f->GetX();
	    Double_t *y0 = f->GetY0();
	    Double_t *y1 = f->GetY1();
	    for (int i=0; i<np[is]; i++) {
	      if (debug>1) {
		cout << " i " << i
		     << " x " << x[i]
		     << " y0 " << y0[i]
		     << " y1 " << y1[i]
		     << endl;
	      }
	      if (min[is]>y0[i]) min[is]=y0[i];
	      if (max[is]<y0[i]) max[is]=y0[i];
	    }
	  }
	}
      }
      
      for (iSM=0; iSM<kNSM; iSM++) {
	if (debug>1) {
	  cout << " iSM " << iSM << endl;
	}
	aveMin = 0;
	aveMax = 0;
	nOk = 0;
	nPTot = 0;

	double tmp[kNumSensPerSM] = {0};
	for (int is=0; is<kNumSensPerSM; is++) {
	  id = is + iSM*kNumSensPerSM;
	  if (np[id]>0 && max[id]<kTooHigh && min[id]>kTooLow) { // some OK readings
	    tmp[nOk] = min[id];
	    nOk++;
	  }
	}
	double median = TMath::Median(nOk, tmp);

	nOk = 0; // reset
	for (int is=0; is<kNumSensPerSM; is++) {
	  id = is + iSM*kNumSensPerSM;
	  if (np[id]>0 &&
	      TMath::Abs(median - max[id])<kMaxDiffMedian &&
	      TMath::Abs(median - min[id])<kMaxDiffMedian &&
	      max[id]<kTooHigh && min[id]>kTooLow) { // some OK readings
	    nOk++;
	    nPTot += np[id];
	    if (debug>1) {
	      printf("is %d np %d median %3.2f min %3.2f max %3.2f diff %3.2f\n",
		     is, np[id], median, min[id], max[id],
		     max[id] - min[id]);
	    }
	    aveMin += min[id];
	    aveMax += max[id];
	  }

	}

	if (nOk > 0) {
	  aveMin /= nOk;
	  aveMax /= nOk;

	  if (debug>0) {
	    printf("iSM %d average min %3.2f max %3.2f (max+min)/2 %3.2f\n",
		   iSM, aveMin, aveMax, (aveMin + aveMax)/2.);
	  }
	}

	tree->Fill();
      } // SM loop

    }
  }
  cout << " nruns " << nruns << endl;

  destFile.cd();
  tree->Write();
  destFile.Close();

  return;
}

