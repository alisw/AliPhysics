/******************************************************************************

  "AliITSOccupancy.C"
  
  this macro calculates the mean occupancy of each ITS layer, counting the
  "fired" digits of each module, and making two overall calculations:
     1) the calculus of the mean overall layer occupancy (as the ratio
	     between the total number of active digits and the total number of 
		  channels in the module;
	  2) a distribution of the occupancies as a funcion of z, to obtain some
		  kind of most significand data for this value, along the z axis
	
  The macro creates also the GIF and the EPS of the TCanvas where the 
  histograms are plotted. Muliple input files are supported (see the 
  input argument list below).
			  
                                                                        
  It is also possible (and recommended) to compile this macro, 
  to execute it faster. To do this, it is necessary to:
 
  1) type, at the root prompt, the instruction 
     gSystem->SetIncludePath("-I- -I$ALICE_ROOT/include -I$ALICE_ROOT/ITS -g")
	  to adjoin the ALICE include libraries to the main include directory
  2) load the function via the ".L AliITSOccupancy.C++" statement 
     (only the first time, because the first time the macro must be compiled). 
      From the second time you use the macro, you must only load it via the 
	  ".L AliITSOccupancy.C+" instruction (note the number of '+' signs in 
      each case
  3) call the function with only its name, e.g.: AliITSOccupancy()

 Author: Alberto Pulvirenti                                             
******************************************************************************/


#if !defined(__CINT__) || defined(__MAKECINT__)
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

TFile* AccessFile(TString inFile="galice.root", TString acctype="R");

Int_t AliITSOccupancy(TString FileHits="galice.root", TString FileDigits="galice.root", Int_t evNum = 0) {       
		 	
  // Open the Root input file containing the AliRun data and 
  // the Root output data file
  TFile *froot = AccessFile(FileHits);
  froot->ls();

  if(!(FileDigits.Data() == FileHits.Data())){
    gAlice->SetTreeDFileName(FileDigits);
  }
  gAlice->GetEvent(evNum);
  // Initialize the ITS and its geometry
  AliITS *ITS = (AliITS*)gAlice->GetModule("ITS");
  AliITSgeom *gm = ITS->GetITSgeom();

  // Variables and objects definition and setting
  Float_t zmax[] = {20.0,20.0,25.0,35.0,49.5,54.0}; // z edges for TH1F (cm)
  Int_t nbins[]  = {40,40,50,70,99,108};             // bins number for TH1Fs
	
  // Histos for plotting occupancy distributions
  TH1F *above[6], *below[6];
	
  for (Int_t lay = 0; lay < 6; lay++) {
    Int_t nlads = gm->GetNladders(lay+1);
    Int_t ndets = gm->GetNdetectors(lay+1);
    Int_t dtype = lay / 2;
    Int_t minID = gm->GetModuleIndex(lay+1, 1, 1);
    Int_t maxID = gm->GetModuleIndex(lay+1, nlads, ndets);
    Text_t ffname[20];
    sprintf(ffname, "h_%d", lay+1);
    below[lay] = new TH1F(ffname, "Total z distribution of digits", nbins[lay], -zmax[lay], zmax[lay]);
    cout << "Evaluating total channels number of layer " << lay+1 << endl;
    for (Int_t I = minID; I <= maxID; I++) {		
      AliITSDetType *det = ITS->DetType(dtype);
      AliITSsegmentation *seg = det->GetSegmentationModel();
      Int_t NX = seg->Npx();
      if(lay == 2 || lay == 3) NX = 192;
      Int_t NZ = seg->Npz();
      //			cout << "ID: " << I << ", Layer: " << lay+1 << ", NX = " << NX << ", NZ = " << NZ << endl;
      for (Int_t ix = 0; ix <= NX; ix++) {
        for (Int_t iz = 0; iz <= NZ; iz++) {
          Float_t lx[] = {0.0,0.0,0.0}, gx[] = {0.0,0.0,0.0};
          seg->DetToLocal(ix, iz, lx[0], lx[2]);
          gm->LtoG(I, lx, gx);
          below[lay]->Fill(gx[2]);
        }
      }
    }
    sprintf(ffname, "H_%d", lay+1);
    above[lay] = new TH1F(ffname, "histo", nbins[lay], -zmax[lay], zmax[lay]);
  }
		
  // Counting the hits, digits and recpoints contained in the ITS
  TTree *TD = gAlice->TreeD();
  ITS->ResetDigits();
  Float_t mean[6];
  Float_t average[6];

  for (Int_t L = 0; L < 6; L++) {
		
    cout << "Layer " << L + 1 << ": " << flush;				

    // To avoid two nested loops, are calculated 
    // the ID of the first and last module of the L
    // (remember the L goest from 0 to 5 (not from 1 to 6)
    Int_t first, last;
    first = gm->GetModuleIndex(L + 1, 1, 1);
    last = gm->GetModuleIndex(L + 1, gm->GetNladders(L + 1), gm->GetNdetectors(L + 1));
				
    // Start loop on modules
    for (Int_t ID = first; ID <= last; ID++) {
      Int_t dtype = L / 2;
      AliITSDetType *det = ITS->DetType(dtype);
      AliITSsegmentation *seg = det->GetSegmentationModel();
      if (dtype == 2) seg->SetLayer(L+1);
			
      TD->GetEvent(ID);
      TClonesArray *digits_array = ITS->DigitsAddress(dtype);
      Int_t digits_num = digits_array->GetEntries();
      // Get the coordinates of the module
      for (Int_t j = 0; j < digits_num; j++) {
        Float_t lx[] = {0.0,0.0,0.0}, gx[] = {0.0,0.0,0.0};
        AliITSdigit *digit = (AliITSdigit*)digits_array->UncheckedAt(j);
        Int_t iz=digit->fCoord1;  // cell number z
        Int_t ix=digit->fCoord2;  // cell number x
        // Get local coordinates of the element (microns)
        seg->DetToLocal(ix, iz, lx[0], lx[2]);
        gm->LtoG(ID, lx, gx);
        above[L]->Fill(gx[2]);
      }
    }
		
    Float_t occupied = above[L]->GetEntries();
    Float_t total = below[L]->GetEntries();
    cout << "Entries: " << occupied << ", " << total << endl;
    average[L] = 100.*occupied/total;
    above[L]->Divide(above[L], below[L], 100.0, 1.0);
    mean[L] = above[L]->GetSumOfWeights() / above[L]->GetNbinsX(); 
    cout.setf(ios::fixed);
    cout.precision(2);
    cout << " digits occupancy = " << mean[L] << "%" << endl;
    cout << " average digits occupancy = " << average[L] << "%" << endl;
  }
	
  TCanvas *view = new TCanvas("view", "Digits occupancy distributions", 600, 900);
  view->Divide(2, 3);
	
  for (Int_t I = 1; I < 7; I++) {
    view->cd(I);
    Text_t title[50];
    sprintf(title, "Layer %d: %4.2f%c", I, mean[I-1], '%');
    title[6] = (Text_t)I + '0';
    above[I-1]->SetTitle(title);
    above[I-1]->SetStats(kFALSE);
    above[I-1]->SetXTitle("z (cm)");
    above[I-1]->SetYTitle("%");
    above[I-1]->Draw();
    view->Update();
  }

  view->SaveAs("AliITSOccupancy_digit.gif");
  view->SaveAs("AliITSOccupancy_digit.eps");
  TFile *fOc = new TFile("AliITSOccupancy.root","recreate");
  fOc->cd();
  for (Int_t I = 0; I < 6; I++) above[I]->Write();
  fOc->Close();

  return 1;
}

//-------------------------------------------------------------------
TFile * AccessFile(TString FileName, TString acctype){

  // Function used to open the input file and fetch the AliRun object

  TFile *retfil = 0;
  TFile *file = (TFile*)gROOT->GetListOfFiles()->FindObject(FileName);
  if (file) {file->Close(); delete file; file = 0;}
  if(acctype.Contains("U")){
    file = new TFile(FileName,"update");
  }
  if(acctype.Contains("N") && !file){
    file = new TFile(FileName,"recreate");
  }
  if(!file) file = new TFile(FileName);   // default readonly
  if (!file->IsOpen()) {
	cerr<<"Can't open "<<FileName<<" !" << endl;
	return retfil;
  } 

  // Get AliRun object from file or return if not on file
  if (gAlice) {delete gAlice; gAlice = 0;}
  gAlice = (AliRun*)file->Get("gAlice");
  if (!gAlice) {
	cerr << "AliRun object not found on file"<< endl;
	return retfil;
  } 
  return file;
}
