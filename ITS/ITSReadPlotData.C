Int_t ITSReadPlotData(char *filename = "galice.root", Int_t evNum = 0) {

	/*********************************************************************
	 *                                                                   *
	 *  Macro used to read hits, digits and recpoints for a module       *
	 *  It draws 3 TH2Fs where stores the 2-D coordinates of these       *
	 *  data for a specified module (for which the macro asks when       *
	 *  it starts, with some checks to avoid wrong detector coords.      *
	 *                                                                   *
	 *  Only a little 'experimental' warning:                            *
	 *  with all the tests I have made, I wasn't able to plot the        *
	 *  digits fo the 5th layer...                                       *
	 *  I skipped this problem with an alert to beware th user, while    *
	 *  in this situation the macro plots only recpoints and hits        *
	 *                                                                   *
	 *  Author: Alberto Pulvirenti                                       *
	 *                                                                   *
	 *********************************************************************/

	
	extern Int_t GetModuleHits (TObject* its, Int_t ID[0], Float_t*& X, Float_t*& Y, Float_t*& Z, Bool_t*& St);
	extern Int_t GetModuleRecPoints (TObject *its, Int_t ID[0], Float_t*& X, Float_t*& Z);
	extern Int_t GetModuleDigits(TObject *its, Int_t ID[0], Int_t dtype, Float_t*& X, Float_t*& Z);
	extern void  AssignCoords(TArrayI *ID);
	
	// First of all, here are put some variable declarations
	// that are useful in the following part:
	Int_t nparticles; // number of particles
	// ITS module coordinates [layer = 1, ladder = 2, det = 3] and absolute ID[0] of module [0]
	TArrayI ID(4);
	Int_t nmodules, dtype; // Total number of modules and module type (SSD, SPD, SDD)
	Float_t *x = 0, *y = 0, *z = 0; // Arrays where to store read coords
	Bool_t *St = 0; // Status of the track (hits only)

	// It's necessary to load the gAlice shared libs
	// if they aren't already stored in memory...
	if (gClassTable->GetID("AliRun") < 0) {
   	gROOT->LoadMacro("loadlibs.C");
		loadlibs();
	}
  // Anyway, this macro needs to read a gAlice file, so it
  // clears the gAlice object if there is already one in memory...
  else {
		if(gAlice){
			delete gAlice;
			gAlice = 0;
		}
	}

	// Now is opened the Root input file containing Geometry, Kine and Hits
  // by default its name must be "galice.root".
  // When the file is opened, its contens are shown.
	TFile *file = (TFile*)gROOT->GetListOfFiles()->FindObject(filename);
	if (!file) file = new TFile(filename);
		file->ls();
	
	// Then, the macro gets the AliRun object from file.
	// If this object is not present, an error occurs
	// and the execution is stopped.
	// Anyway, this operation needs some time,
	// don't worry about an apparent absence of output and actions...
	cout << "\nSearching in '" << filename << "' for an AliRun object ... " << flush;
	gAlice = (AliRun*)file->Get("gAlice");
	if (gAlice)
		cout << "FOUND!" << endl;
	else {
		cout<<"NOT FOUND! The Macro can't continue!" << endl;
		return 0;
	}
	
	// Then, the macro selects the event number specified. Default is 0.
	nparticles = gAlice->GetEvent(evNum);
	cout << "\nNumber of particles   = " << nparticles << endl;
	if (!nparticles) {
		cout << "With no particles I can't do much... Goodbye!" << endl;
		return 0;
	}
	
	// The next step is filling the ITS from the AliRun object.
	AliITS *ITS = (AliITS*)gAlice->GetModule("ITS");
  ITS->InitModules(-1, nmodules);
  cout << "Number of ITS modules = " << nmodules << endl;
	cout << "\nFilling modules (it takes a while, now)..." << flush;
	ITS->FillModules(0, 0, nmodules, " ", " ");
  cout << "DONE!" << endl;
	AliITSgeom *gm = ITS->GetITSgeom();
	AliITSDetType *det = ITS->DetType(dtype);	
	AliITSsegmentation *seg = det->GetSegmentationModel();	
	
	for(;;) {

    // Input phase.
    // The macro asks if the user wants to put a detector ID[0]
    // or prefers to input layer, ladder and detector.
    for (Int_t i = 0; i < 4; i++) ID[i] = -1;
    Int_t answ;
    do {
	   cout << "\nSelection modes:" << endl;
	  	cout << "1) by layer - ladder - detector numbers" << endl;
    	cout << "2) by unique detector ID" << endl;
	   cout << "0) exit macro" << endl;
	  	cout << "\nEnter your choice: ";
    	cin >> answ;
    } while (answ < 0 || answ > 2);
    switch (answ) {
    	case 0:
    		// input = 0 ---> EXIT
    		return;
    		break;
    	case 1:
    		// input = 1 ---> SELECTION BY COORD
    		do {
					cout << "\nLayer number [1-6, 0 to exit]: ";
					cin >> ID[1];
					if (!ID[1]) return 0;
				} while (ID[1] < 0 || ID[1] > 6);
				
				// Detector type: 0 --> SPD, 1 --> SDD, 2 --> SSD.
				// Layer 1,2 --> 0 / Layer 3,4 --> 1 / Layer 5,6 --> 2
				dtype = ID[1] / 3;
				
				// Once fixed the layer number, the macro calculates the max number
				// for ladder and detector from geometry, and accepts only suitable values.
				do {
					ID[2] = gm->GetNladders(ID[1]);
					cout << "Ladder number [1-" << ID[2] << ", 0 to exit]: ";
					cin >> ID[2];
					if (!ID[2]) return 0;
				} while (ID[2] < 0 || ID[2] > gm->GetNladders(ID[1]));
				do {
					ID[3] = gm->GetNdetectors(ID[1]);
					cout << "Detector number [1-" << ID[3] << ", 0 to exit]: ";
					cin >> ID[3];
					if (!ID[3]) return 0;
				} while (ID[3] < 0 || ID[3] > gm->GetNdetectors(ID[1]));
				break;
    	case 2:
    		// input = 2 ---> SELECTION BY ID[0]
    		do {
					ID[0] = gm->GetIndexMax();
					cout << "\n Detector ID number [0-" << ID[0] << ", -1 to exit]: ";
					cin >> ID[0];
					if (ID[0] == -1) return 0;
				} while (ID[0] < 0 || ID[0] > gm->GetIndexMax());
    	  break;
    };

    if (ID[0] == -1)
    		// If ID[0] = -1 the chioce was by coords, so it's necessary to assign the ID:
			ID[0] = gm->GetModuleIndex(ID[1], ID[2], ID[3]);
		else {	
			// Else we must get the layer, ladder and detector by the ID;
			// ...but first we must remember that the ID goest from 0 to NModules - 1!
			ID[0]--;
			ID[1] = ITS->GetModule(ID[0])->GetLayer();
			ID[2] = ITS->GetModule(ID[0])->GetLadder();
			ID[3] = ITS->GetModule(ID[0])->GetDet();
		}
				
		// Defines the histograms inside the `for' cycle, so they are destroyed at the end
		// of every read sequqnce, in order to mek another withour segmentation faults
		Text_t msg[250], xm = 0.0, ym = 0.0;
		switch (dtype) {
			case 0: xm = 1.5; zm = 7.0; break;
			case 1: xm = 7.5; zm = 8.0; break;
			case 2: xm = 7.5; zm = 4.5; break;
		}
		sprintf(msg, "Module index=%d lay=%d lad=%d det=%d", ID[0], ID[1], ID[2], ID[3]);
		TH2F *hhits = new TH2F("hhits", msg, 500, -xm, xm, 500, -zm, zm);     // Histogram of hits
		TH2F *hrecs = new TH2F("hrecs", msg, 500, -xm, xm, 500, -zm, zm);     // Histogram of recpoints
		TH2F *hdigits = new TH2F("hdigits", msg, 500, -xm, xm, 500, -zm, zm); // Histogram of digits
		
		cout << endl;
		
		// Reads hits...
		Int_t hits = GetModuleHits(ITS, ID[0], x, y, z, St);
		if (!hits) {
			cout << "No hits in module!" << endl;
			continue;
		}
		for (Int_t i = 0; i < hits; i++) if (!St[i]) hhits->Fill(x[i], z[i]);
		
		// Reads recpoints...
		Int_t recs = GetModuleRecPoints(ITS, ID[0], x, z);
		if (!recs) {
			cout << "No recpoints in module!" << endl;
			continue;
		}
		for (Int_t i = 0; i < recs; i++) hrecs->Fill(x[i], z[i]);
		
		// Reads digits...
		Int_t digits = GetModuleDigits(ITS, ID[0], dtype, x, z);
		if (!digits) {
			cout << "No digits in module!" << endl;
			//continue;
		}
		for (Int_t i = 0; i < digits; i++) hdigits->Fill(x[i], z[i]);

		// Draws read data...
		// HITS -------> red (2) crosses.
		// DIGITS -----> green (8) boxes.
		// REC-POINTS -> blue (4) St. Andrew's crosses.

		TCanvas *viewer = new TCanvas("viewer", "Module viewer canvas", 0, 0, 800, 800);
		viewer->cd();
		
		hdigits->SetMarkerStyle(25);
		hdigits->SetMarkerColor(8);
		hdigits->SetMarkerSize(2);
		hdigits->SetStats(kFALSE);
		hdigits->SetXTitle("Local X (cm)");
		hdigits->SetYTitle("Local Z (cm)");
		hdigits->Draw();
		
		hhits->SetMarkerStyle(5);
		hhits->SetMarkerColor(2);
		hhits->SetMarkerSize(3);
		hhits->SetStats(kFALSE);
		hhits->Draw("same");
	
		hrecs->SetMarkerStyle(2);
		hrecs->SetMarkerColor(4);
		hrecs->SetMarkerSize(3);
		hrecs->SetStats(kFALSE);
		hrecs->Draw("same");
		
		TLegend *legend = new TLegend(0.7, 0.8, 0.9, 0.9);
		legend->SetMargin(0.2);
		legend->AddEntry(hhits, "hits","P");
		legend->AddEntry(hrecs, "recpoints","P");
		legend->AddEntry(hdigits, "digits","P");
		legend->SetTextSizePixels(14);
		legend->Draw();
		
		viewer->Update();
	}
	
	cout << "Done. Goodbye" << endl;
	return;
}

Int_t GetModuleHits (TObject* its, Int_t ID, Float_t*& X, Float_t*& Y, Float_t*& Z, Bool_t*& St) {	
	// First of all, the macro selects the specified module,
	// then gets the array of hits in it and their number.
	AliITS *ITS = (AliITS*) its;
	AliITSmodule *module = ITS->GetModule(ID);
	TObjArray *hits_array = module->GetHits();
	Int_t hits_num = hits_array->GetEntriesFast();
	
	// Now, if this count returns 0, there's nothing to do,
	// while, if it doesn't, the first thing to do is dimensioning
	// the coordinate arrays, and then the loop can start.
	if (!hits_num)
		return 0;
	else {
		if (X) delete [] X;	
		if (Y) delete [] Y;
		if (Z) delete [] Z;
		if (St) delete [] St;
		X = new Float_t[hits_num];
		Y = new Float_t[hits_num];
		Z = new Float_t[hits_num];
		St = new Int_t[hits_num];
	}

	for (Int_t j = 0; j < hits_num; j++) {
		AliITShit *hit = (AliITShit*) hits_array->At(j);
		X[j]  = hit->GetXL();
		Y[j]  = hit->GetYL();
		Z[j]  = hit->GetZL();
		St[j] = hit->StatusEntering();
	}
	return hits_num;
}

Int_t GetModuleRecPoints (TObject *its, Int_t ID, Float_t*& X, Float_t*& Z) {
	
	// First of all, the macro selects the specified module,
	// then gets the array of recpoints in it and their number.
	AliITS *ITS = (AliITS*) its;
	TTree *TR = gAlice->TreeR();
	ITS->ResetRecPoints();
	TR->GetEvent(ID);
	TClonesArray *recs_array = ITS->RecPoints();
	Int_t recs_num = recs_array->GetEntries();

	// Now, if this count returns 0, there's nothing to do,
	// while, if it doesn't, the first thing to do is dimensioning
	// the coordinate and energy loss arrays, and then the loop can start.
	if (!recs_num)
		return 0;
	else {
		if (X) delete [] X;	
		if (Z) delete [] Z;
		X = new Float_t[recs_num];
		Z = new Float_t[recs_num];
	}
	for(Int_t j = 0; j < recs_num; j++) {
		AliITSRecPoint *recp = (AliITSRecPoint*)recs_array->At(j);
	   X[j] = recp->GetX();
		Z[j] = recp->GetZ();
	}
	return recs_num;	
}

Int_t GetModuleDigits(TObject *its, Int_t ID, Int_t dtype, Float_t*& X, Float_t*& Z) {

	// First of all, the macro selects the specified module,
	// then gets the array of recpoints in it and their number.
	AliITS *ITS = (AliITS*)its;
	TTree *TD = gAlice->TreeD();
	ITS->ResetDigits();
	TD->GetEvent(ID);
	TClonesArray *digits_array = ITS->DigitsAddress(dtype);
	AliITSgeom *gm = ITS->GetITSgeom();	
	AliITSDetType *det = ITS->DetType(dtype);	
	AliITSsegmentation *seg = det->GetSegmentationModel();	
	TArrayI ssdone(5000);  // used to match p and n side digits of strips
	TArrayI pair(5000);    // as above 	
	Int_t digits_num = digits_array->GetEntries();
	// Now, if this count returns 0, there's nothing to do,
	// while, if it doesn't, the first thing to do is dimensioning
	// the coordinate and energy loss arrays, and then the loop can start.

	if (!digits_num)
		return 0;
	else {
		if (X) delete [] X;			
		if (Z) delete [] Z;
		X = new Float_t[digits_num];		
		Z = new Float_t[digits_num];
	}
	
	// Get the coordinates of the module
	if (dtype == 2) {
		for (Int_t j = 0; j < digits_num; j++) {
			ssdone[j] = 0;			
			pair[j] = 0;
		}
	}
  for (Int_t j = 0; j < digits_num; j++) {
  	cout << j << endl;
		digit = (AliITSdigit*)digits_array->UncheckedAt(j);
		Int_t iz=digit->fCoord1;  // cell number z
		Int_t ix=digit->fCoord2;  // cell number x
    // Get local coordinates of the element (microns)
		if(dtype < 2)
    	seg->GetPadCxz(ix, iz, X[j], Z[j]);
    else {
			// SSD: if iz==0 ---> N side; if iz==1 P side
      if (ssdone[j] == 0) {
				ssdone[j]=1;
				pair[j]=-1;
				Bool_t pside = (iz == 1);
				Bool_t impaired = kTRUE;
				Int_t pstrip = 0;
				Int_t nstrip = 0;
				if (pside) pstrip = ix; else nstrip = ix;
				for (Int_t k = 0; k < digits_num; k++) {
					if (ssdone[k] == 0 && impaired) {
						AliITSdigitSSD *sdigit=(AliITSdigitSSD*)digits_array->UncheckedAt(k);
						if (sdigit->fCoord1 != iz && sdigit->GetTracks()[0] == digit->GetTracks()[0]) {
							ssdone[k]=2;
							pair[j]=k;
							if (pside) nstrip = sdigit->fCoord2; else pstrip = sdigit->fCoord2;
							impaired=kFALSE;
						}
					}
				}
        if (!impaired) seg->GetPadCxz(pstrip, nstrip, X[j], Z[j]);
			}
		}
		if (dtype == 0) {
			// !!!THIS CONVERSION TO HIT LRS SHOULD BE REMOVED AS SOON AS THE CODE IS FIXED
			X[j] = X[j]-seg->Dx() / 2.0;
			Z[j] = Z[j]-seg->Dz() / 2.0;
		}
		if (dtype != 1) {
			X[j] /= 10000.0;
			Z[j] /= 10000.0;
		}
	}
	return digits_num;
} 

