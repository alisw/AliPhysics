/******************************************************************************

  "ITSOccupancy.C"
  
  this macro calculates the mean occupancy of each ITS layer,       
  making also a distribution in function of the z-value of the           
  "fired" digits for each layer                                              
                                                                         
  Because of evaluation problems, it is necessary to permit the          
  analysis for both hits, digits and recpoints, so it's provided an     
  option data which must contain                                         
  - "H" for hits                                                         
  - "D" for digits                                                       
  - "R" for recpoints                                                    
  but it's possible to use several option in the same time               
  (like "HD", "HRD", and so on...)                                       
                                                                         
  to avoid the problem of drawing unrequested windows, the graph are     
  plotted only if it is explictly requested with the option "P", that    
  can be inserted together to the ones described above.                  
                                                                         
  It is also possible to compile this macro, to execute it faster.       
  To do this, it is necessary to:
  1) make sure that de "#define __COMPILED__" instruction below is UNcommented
     and the "#define __NOCOMPILED__" instruction is commented (remember that
	  if you comment or uncomment both these instructions, the macro can't work)
  2) type, at the root prompt, the instruction 
     gSystem->SetIncludePath("-I- -I$ALICE_ROOT/include -I$ALICE_ROOT/ITS -g")
	  to adjoin the ALICE include libraries to the main include directory
  3) load the function via the ".L ITSOccupancy.C++" statement (only the first
     time, because the first time the macro must be compiled). Frome the 
	  second time you use the macro, you must only load it via the 
	  ".L ITSOccupancy.C+" instruction (note the number of '+' signs in each case
  4) call the function with only its name, es: ITSOccupancy("HP", "galice.root")

 NOTE 1: option interpretation is case-insensitive ("h" eqv "H", etc.)    
 NOTE 2: if you insert a filename as argument, it must be inserted only the name
         (the macro attaches the extension ".root" to it, every time)

 Author: Alberto Pulvirenti                                             
******************************************************************************/

#define __COMPILED__
//#define __NOCOMPILED__

#ifdef __COMPILED__
	#include <iostream.h>
	#include <TROOT.h>
	#include <TArrayI.h>
	#include <TCanvas.h>
	#include <TClassTable.h>
	#include <TClonesArray.h>
	#include <TFile.h>
	#include <TObject.h>
	#include <TObjArray.h>
	#include <TTree.h>
	#include <TMath.h>
	#include <TString.h>
	#include <TH1.h>
	#include <AliRun.h>
	#include <AliITS.h>
	#include <AliITSgeom.h>
	#include <AliITSDetType.h>
	#include <AliITSRecPoint.h>
	#include <AliITSdigit.h>
	#include <AliITShit.h>
	#include <AliITSmodule.h> 
	#include <AliITSsegmentation.h>
	#include <AliITSsegmentationSPD.h> 
	#include <AliITSsegmentationSDD.h>
	#include <AliITSsegmentationSSD.h>
#endif

Int_t ITSOccupancy(char* opt, char *name = "galice", Int_t evNum = 0) {       

	extern void GetModuleDigits(TObject *its, Int_t ID, Int_t dtype, TH1D *counter);
	extern void GetModuleHits (TObject* its, Int_t ID, TH1D *counter);
	extern void GetModuleRecPoints (TObject *its, Int_t ID, TH1D *counter);
	extern Int_t DType(Int_t layer);
		 
	// Evaluating options
	TString option(opt);
	Bool_t HITS = (option.Contains("h") || option.Contains("H"));
	Bool_t DIGS = (option.Contains("d") || option.Contains("D"));
	Bool_t RECS = (option.Contains("r") || option.Contains("R"));
	Bool_t PLOT = (option.Contains("p") || option.Contains("P"));
	
	Text_t filename[100];
	char  gifH[100], gifD[100], gifR[100];
	strcpy (filename, name);
	strcpy (gifH, name);
	strcpy (gifD, name);
	strcpy (gifR, name);
	strcat (filename, ".root");
	strcat (gifH, "H.gif");
	strcat (gifD, "D.gif");
	strcat (gifR, "R.gif");
	
	if (!HITS && !DIGS && !RECS) {
		cout << "\n\n--- NO DATA-TYPE SPECIFIED!!! ---\n\n";
		return 0;
	}
	
	// --------------------------------
	//  1. ALICE objects setting phase
	// --------------------------------
	
	// Load the gAlice shared libs if not already in memory
#ifdef __NOCOMPILED__
	if (gClassTable->GetID("AliRun") < 0) {
    	gROOT->LoadMacro("loadlibs.C");
    	loadlibs();
  	}
#endif
   if(gAlice){
      delete gAlice;
	   gAlice=0;
   }
	// Open the Root input file containing the AliRun data 
	// (default name is "galice.root")
	TFile *file = (TFile*)gROOT->GetListOfFiles()->FindObject(filename);
	if (!file) file = new TFile(filename);
	file->ls();
	// Get the AliRun object from file (if there's none, the macro is stopped)
	cout << "\nSearching in '" << filename << "'" << endl;
	gAlice = (AliRun*)file->Get("gAlice");
	if (gAlice)
		cout << "Ok, I found an AliRun object..." << endl;
	else {
		cout<<"Sorry, there isn't any AliRun object here..." << endl;
		return 0;
	}
	// Select the event number specified. Default is 0.
	Int_t nparticles = gAlice->GetEvent(evNum);
	cout << "\nNumber of particles   = " << nparticles << endl;
	if (!nparticles) {
		cout << "With no particles I can't do much..." << endl;
		return 0;
	}
	// Initialize the ITS and its geometry
	AliITS *ITS = (AliITS*)gAlice->GetModule("ITS");
	AliITSgeom *gm = ITS->GetITSgeom();
	// Fill the AliITSmodule objects (only for hits)
	Int_t nmodules;
	if (HITS) {
	  	ITS->InitModules(-1, nmodules);
  		cout << "Number of ITS modules = " << nmodules << endl;
		cout << "\nFilling modules (it takes a while, now)..." << flush;
		ITS->FillModules(0, 0, nmodules, " ", " ");
	  	cout << "DONE!" << endl;	
	}

	// -------------------------------------------------
	//  B. Variables and objects definition and setting
	// -------------------------------------------------
	
	Double_t zmax[6] = {16.5, 16.5, 22.2, 29.7, 45.1, 50.8};  // max z measured (cm)
	Double_t binw[] = {2., 2., 2., 2., 5., 5.}; // bin width for TH1Ds (cm)
	
	Double_t tot[6];  // total number of channels 
	
	// Histos for plotting hits [1], digits [2] and points [3], and divisor [0]
	TH1D *hist[4][6];
	
	// Histograms initialization:
	// Using the AliITSsegmentation object, here is deduced
	// the maximum z that can be measurable for each layer, simply
	// obtaining the globar coordinates of the first digit of the first module 
	// of each layer (that is placed at the maximum z).
	for (Int_t i = 0; i < 6; i++) {
		AliITSDetType *det = ITS->DetType(DType(i));	
		AliITSsegmentation *seg = det->GetSegmentationModel();
		tot[i] = gm->GetNladders(i+1) * gm->GetNdetectors(i+1) * seg->GetNPads();
		cout.setf(ios::fixed);
		cout << tot[i] << endl;
		Text_t hname[7];
		for (Int_t j = 0; j < 4; j++) {
			hname[0] = 'h';
			hname[1] = 'i';
			hname[2] = 's';
			hname[3] = 't';
			hname[4] = '0' + j;
			hname[5] = '0' + i;
			hname[6] = '\0';
			Int_t nbins = (Int_t)(2. * zmax[i] / binw[i]);
			hist[j][i] = new TH1D(hname, "histo", nbins, -zmax[i], zmax[i]);
		}
	}
	// Here is put the value of the density for each bin of each layer's TH1D
	for (Int_t i = 0; i < 6; i++) {
		for (Int_t j = 0; j < hist[0][i]->GetNbinsX(); j++) {
			Double_t z = hist[0][i]->GetBinLowEdge(j);
			Double_t bincontent = tot[i] / (zmax[i] * 2.0) * hist[0][i]->GetBinWidth(j);
			hist[0][i]->Fill(z + 0.1, bincontent / 100.0);
			z += hist[0][i]->GetBinWidth(j);
		}
	}
	
	// --------------------------------------------------------------
	//  Counting the hits, digits and recpoints contained in the ITS
	// --------------------------------------------------------------

	for (Int_t L = 0; L < 6; L++) {
		
		// To avoid two nested loops, are calculated 
		// the ID of the first and last module of the L
		// (remember the L goest from 0 to 5 (not from 1 to 6)
		Int_t first, last;
		first = gm->GetModuleIndex(L + 1, 1, 1);
		last = gm->GetModuleIndex(L + 1, gm->GetNladders(L + 1), gm->GetNdetectors(L + 1));
				
		// Start loop on modules
		cout << "Examinating layer " << L + 1 << " ... " << flush;
		for (Int_t ID = first; ID <= last; ID++) {
			
			// Hits analysis (if requested)
			if (HITS)
				GetModuleHits(ITS, ID, hist[2][L]);
			
			// Digits analysis (if requested)
			if (DIGS)
				GetModuleDigits(ITS, ID, DType(L), hist[1][L]);
			
			// Recpoints analysis (if requested)
			if (RECS) 
				GetModuleRecPoints(ITS, ID, hist[3][L]);
		}
		
		hist[1][L]->Divide(hist[0][L]);
		hist[2][L]->Divide(hist[0][L]);
		hist[3][L]->Divide(hist[0][L]);
		cout << "DONE" << endl;
	}

	// Write on the screen the total means
	cout << endl;
	cout << "********* MEAN OCCUPANCIES *********" << endl;
	for (Int_t L = 0; L < 6; L++) {
		Double_t mean[3];
		cout << " LAYER " << L << ": " << endl;
		for (Int_t i = 0; i < 3; i++) {
			mean[i] = hist[i+1][L]->GetSumOfWeights() / hist[i+1][L]->GetNbinsX(); 
			Text_t title[50];
			sprintf(title, " - Layer %d, mean %4.2f - ", L + 1, mean[i]);
			hist[i+1][L]->SetTitle(title);
			cout.setf(ios::fixed);
			cout.precision(3);
			if (HITS && i == 1) {
				cout << "     hits occupancy = ";
				if (mean[i] < 10.0) cout << ' ';
				cout << mean[i] << "%" << endl;
			}
			else if (DIGS && i == 0) {
				cout << "   digits occupancy = ";
				if (mean[i] < 10.0) cout << ' ';
				cout << mean[i] << "%" << endl;
			}
			if (RECS && i == 2) {
				cout << "   recpts occupancy = ";
				if (mean[i] < 10.0) cout << ' ';
				cout << mean[i] << "%" << endl;
			}
		}
		cout << "------------------------------------" << endl;
	}
	
	// ----------------------
	//  Plots (if requested)
	// ----------------------
	if (!PLOT) return 1;
	
	TCanvas *view[4];
	if (HITS) {
		view[2] = new TCanvas("viewH", "Occupancy view (HITS)", 0, 0, 1050, 700);
		view[2]->Divide(3,2, 0.001, 0.001);
	}	
	if (DIGS) {
		view[1] = new TCanvas("viewD", "Occupancy view (DIGITS)", 20, 40, 1050, 700);
		view[1]->Divide(3,2, 0.001, 0.001);
	}
	if (RECS) {
		view[3] = new TCanvas("viewR", "Occupancy view (RECPOINTS)", 40, 60, 1050, 700);
		view[3]->Divide(3,2, 0.001, 0.001);
	}
	
	for (Int_t L = 0; L < 6; L++) {
		for (Int_t i = 1; i <= 3; i++) {
			if (DIGS && i == 1) 
				hist[i][L]->SetFillColor(kBlue);
			else if (HITS && i == 2) 
				hist[i][L]->SetFillColor(kRed);
			else if (RECS && i == 3) 
				hist[i][L]->SetFillColor(kGreen);
			
			if ((HITS && i == 2) || (DIGS && i == 1) || (RECS && i == 3)) {			
				view[i]->cd(L+1);			
				hist[i][L]->SetStats(kFALSE);
				hist[i][L]->Draw();
				hist[i][L]->GetXaxis()->SetTitle("zeta");
				hist[i][L]->GetYaxis()->SetTitle("%");
				view[i]->Update();
			}
		}
	}
	
	if (HITS) view[2]->SaveAs(gifH);
	if (DIGS) view[1]->SaveAs(gifD);
	if (RECS) view[3]->SaveAs(gifR);
	
	return 1;
}

void GetModuleDigits(TObject *its, Int_t ID, Int_t dtype, TH1D *counter) {
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
	Int_t digits_num = digits_array->GetEntries();
	
	// Get the coordinates of the module
  	for (Int_t j = 0; j < digits_num; j++) {
		Float_t lx[] = {0.0,0.0,0.0}, gx[] = {0.0,0.0,0.0};
		AliITSdigit *digit = (AliITSdigit*)digits_array->UncheckedAt(j);
		Int_t iz=digit->fCoord1;  // cell number z
		Int_t ix=digit->fCoord2;  // cell number x
    	// Get local coordinates of the element (microns)
    	seg->GetPadCxz(ix, iz, lx[0], lx[2]);
		if (dtype != 1) {
			// !!!THIS CONVERSION TO HIT LRS SHOULD BE REMOVED AS SOON AS THE CODE IS FIXED
			if (!dtype) {
				lx[0] -= seg->Dx() / 2.0;
				lx[2] -= seg->Dz() / 2.0;
			}
			lx[0] /= 10000.0; // convert from microns to cm (if not SDD)
			lx[2] /= 10000.0; // convert from microns to cm (if not SDD)
		}
		gm->LtoG(ID, lx, gx);
		counter->Fill(gx[2]);
	}
}

/*void GetModuleDigits(TObject *its, Int_t ID, Int_t dtype, TH1D *counter) {
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
	
	// Get the coordinates of the module
	if (dtype == 2) {
		for (Int_t j = 0; j < digits_num; j++) {
			ssdone[j] = 0;			
			pair[j] = 0;
		}
	}
  	for (Int_t j = 0; j < digits_num; j++) {
		Double_t lx[] = {0.0,0.0,0.0}, gx[] = {0.0,0.0,0.0};
		AliITSdigit *digit = (AliITSdigit*)digits_array->UncheckedAt(j);
		Int_t iz=digit->fCoord1;  // cell number z
		Int_t ix=digit->fCoord2;  // cell number x
    	// Get local coordinates of the element (microns)
		if(dtype < 2) {
    		seg->GetPadCxz(ix, iz, lx[0], lx[2]);
			if (dtype == 0) {
				// !!!THIS CONVERSION TO HIT LRS SHOULD BE REMOVED AS SOON AS THE CODE IS FIXED
				lx[0] -= seg->Dx() / 2.0;
				lx[2] -= seg->Dz() / 2.0;
				lx[0] /= 10000.0; // convert from microns to cm (SPD)
				lx[2] /= 10000.0; // convert from microns to cm (SPD)
			}
			gm->LtoG(ID, lx, gx);
			counter->Fill(gx[2]);
		}
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
        		if (!impaired) {
					seg->GetPadCxz(pstrip, nstrip, lx[0], lx[2]); 
					lx[0] /= 10000.0;
					lx[2] /= 10000.0;
					gm->LtoG(ID, lx, gx);
					counter->Fill(gx[2]);
				}
			}
		}
	}
}*/

void GetModuleHits (TObject* its, Int_t ID, TH1D *counter) {	
	// First of all, the macro selects the specified module,
	// then gets the array of hits in it and their number.
	AliITS *ITS = (AliITS*) its;
	AliITSmodule *module = ITS->GetModule(ID);
	TObjArray *hits_array = module->GetHits();
	Int_t hits_num = hits_array->GetEntriesFast();	
	// next, fills the counters with the hits coordinates' calcula
	for (Int_t j = 0; j < hits_num; j++) {
		AliITShit *hit = (AliITShit*) hits_array->At(j);
		Double_t Z = hit->GetZG();
		if (!hit->StatusEntering()) counter->Fill(Z);
	}
}

void GetModuleRecPoints (TObject *its, Int_t ID, TH1D *counter) {
	// First of all, the macro selects the specified module,
	// then gets the array of recpoints in it and their number.
	AliITS *ITS = (AliITS*) its;
	AliITSgeom *gm = ITS->GetITSgeom();	
	TTree *TR = gAlice->TreeR();
	if (!TR || TR->GetEntries() == 0) {
		cout << "The ITS wasn't clustered..." << endl;
		return;
	}
	ITS->ResetRecPoints();
	TR->GetEvent(ID);
	TClonesArray *recs_array = ITS->RecPoints();
	Int_t recs_num = recs_array->GetEntries();

	for(Int_t j = 0; j < recs_num; j++) {
		Double_t lx[3] = {0.,0.,0.}, gx[3] = {0.,0.,0.};
		AliITSRecPoint *recp = (AliITSRecPoint*)recs_array->At(j);
	   lx[0] = recp->GetX();
		lx[1] = 0.0;
		lx[2] = recp->GetZ();
		gm->LtoG(ID, lx, gx);
		counter->Fill(gx[2]);
	}
}

Int_t DType(Int_t layer) {
	if (layer == 0 || layer == 1)
		return 0;
	else if (layer == 2 || layer == 3)
		return 1;
	else 
		return 2;
}
