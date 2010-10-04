// $Id$

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <Riostream.h>
#include <TFile.h>
#include <TKey.h>
#include <TTree.h>
#include <TMath.h>
#include <TChain.h>
#include <TH2F.h>
#include <TSystem.h>
#include <TTimeStamp.h>
#include <TStyle.h>
#include "TreeClasses.h"
#include "Utils.h"        // PlotUtils and Noti classes
#endif

Int_t nCuts = 2;
Bool_t quit = 0;
Bool_t printPDF = 0;
TFile* outFile = new TFile("out.root", "recreate");
const Double_t PI = TMath::Pi();
Int_t nEvents = -1; // Run over whole chain if < 0
MyHeader     *event   = 0;
TClonesArray *tracks  = 0;
TBranch      *hBranch = 0;
TBranch      *tBranch = 0;

// Event
//TH1F* hRun;		 
TH1F* hBx;		 // crossing times (ns)
TH1F* hNChips1;	 
TH1F* hNChips2;	 
TH1F* hNTracks;	 
TH1F* hNSelTracks;       // no kink daughters, NClTPC > 20.
TH1F* hNTracklets;	 
TH1F* hVx;             
TH1F* hVy;             
TH1F* hVz;               // from V0 (?)
TH1F* hVc;               // from TPC (?)
TH1F* hIsPileupSPD; 	 
TH1F* hNPileupSPD;	 
TH1F* hNPileup;	 
TH1F* hTrClassMask;	 // get names from constantin 
TH1F* hTrCluster;	 
TH1F* hVxSPD;          
TH1F* hVySPD;          
TH1F* hVzSPD;          
TH1F* hVcSPD;          
TH1F* hVxTPC;          
TH1F* hVyTPC;          
TH1F* hVzTPC;          
TH1F* hVcTPC;                     

// Track
TH1F* hSt;	     // status
TH1F* hC;	     // charge
TH1F* hPt;           
TH1F* hEta;          
TH1F* hPhi;          
TH1F* hNClTPC;	     
TH1F* hNClTPC1;	     // ?
TH1F* hNClTPCShared; 
TH1F* hNClITS;	     
TH1F* hChi2TPC;      
TH1F* hChi2TPC1;     
TH1F* hChi2ITS;      
TH1F* hD;            // Transverse DCA or "impact parameter"
TH1F* hZ;            // Z DCA
TH1F* hDTPC;         
TH1F* hZTPC;           

// Other
TH1F* hStBits;
TH1F* hL0Bits;
TH1F* hTrClBits;
// TH1F *hEta1, *hEta2, *hEta3;
// TH1F *hD1, *hD2, *hD3;
// TH1F *hZ1, *hZ2, *hZ3;
TH1F* hChi2ndfTPC;   

TH2F* hDZ;
TH2F* hDZTPC;

void BookHistos();
TObjArray* GetHistList(TFile& inFile, TString clname);
void SetHistProps();
void DrawHistos(TFile& inFile);
void FillEventInfo(const MyHeader& ev);
void FillTrackInfo(const MyPart& trk);

void DataChecker(const char* inFileNames = "../rootfiles/res_LHC10c_09212010/merged_*.root")
{
  // --- TTree/TChain input ---
  Noti *nme = new Noti;
  TChain *tree = new TChain("MyTree");
  tree->Add(inFileNames);
  tree->SetNotify(nme);
  
  Int_t nents = tree->GetEntries();
  cout << "Found " << nents << " entries!" << endl;
  if (nents<=0)
    return;
  if (nEvents>0 && nEvents<nents) { 
    nents=nEvents;
    cout << "Using " << nents << " entries!" << endl;
  }

  if (TClass::GetClass("MyPart"))
    TClass::GetClass("MyPart")->IgnoreTObjectStreamer();
  if (TClass::GetClass("MyTracklet"))
    TClass::GetClass("MyTracklet")->IgnoreTObjectStreamer();

  BookHistos();
  
  // --- Event loop ---
  for (Int_t i=0;i<nents;++i) {
    Int_t li = tree->LoadTree(i);
    if (nme->Notified()) {
      hBranch = tree->GetBranch("header");
      hBranch->SetAddress(&event);
      tBranch = tree->GetBranch("parts");
      tBranch->SetAddress(&tracks);
      nme->Reset();
    }
    
    hBranch->GetEntry(li);
    tBranch->GetEntry(li);

    FillEventInfo(*event);

    // --- Track loop ---
    Int_t ntracks = tracks->GetEntries();
    for (Int_t t=0;t<ntracks;++t) {
      MyPart *trk = (MyPart*)tracks->At(t);

      if (!event->fIsPileupSPD && trk->fPt > 0.4)
	FillTrackInfo(*trk);
    }
  }
  
  outFile->Write();
  SetHistProps();
  DrawHistos(*outFile);

  return;
}

void BookHistos()
{
  // From MyHeader
  // hRun		   = new TH1F("hRun", "", 1e3, 0, 5e5);	      	  	  
 hBx		   = new TH1F("hBx", "", 3000, 0, 3000);		      	  
 hNChips1	   = new TH1F("hNChips1", "", 150, 0, 150);	      
 hNChips2	   = new TH1F("hNChips2", "", 150, 0, 150);	      
 hNTracks	   = new TH1F("hNTracks", "", 500, 0, 500);	      
 hNSelTracks	   = new TH1F("hNSelTracks", "", 500, 0, 500);	      	  
 hNTracklets	   = new TH1F("hNTracklets", "", 100, 0, 100);	      	  
 hVx               = new TH1F("hVx", "", 100, -0.1, 0.0);               
 hVy               = new TH1F("hVy", "", 100, 0.15, 0.25);               
 hVz               = new TH1F("hVz", "", 40, -20, 20);               
 hVc               = new TH1F("hVc", "", 100, 0, 100);               
 hIsPileupSPD 	   = new TH1F("hIsPileupSPD", "", 3, 0, 3);      	  	  
 hNPileupSPD	   = new TH1F("hNPileupSPD", "", 10, 0, 10);	      	  
 hNPileup	   = new TH1F("hNPileup", "", 10, 0, 10);	      
 hTrCluster	   = new TH1F("hTrCluster", "", 10, 0, 10);	      	  
 hVxSPD            = new TH1F("hVxSPD", "", 100, -0.5, 0.5);            
 hVySPD            = new TH1F("hVySPD", "", 100, -0.5, 0.5);            
 hVzSPD            = new TH1F("hVzSPD", "", 100, -20, 20);            
 hVcSPD            = new TH1F("hVcSPD", "", 100, 0, 100);            
 hVxTPC            = new TH1F("hVxTPC", "", 100, -0.01, 0.01);             
 hVyTPC            = new TH1F("hVyTPC", "", 100, -0.01, 0.01);             
 hVzTPC            = new TH1F("hVzTPC", "", 100, -5, 5);            
 hVcTPC            = new TH1F("hVcTPC", "", 100, -1, 1);            
		  		      
 // Other evt histos
 hL0Bits           = new TH1F("hL0Bits", "", 32, 0, 32);   // UInt_t 
 hTrClBits         = new TH1F("hTrClBits", "", 64, 0, 64); // ULong64_t

 // From MyPart
 hC	     	   = new TH1F("hC", "", 5, -2.5, 2.5);	     	      	  
 hPt               = new TH1F("hPt", "", 200, 0, 20);               
 hEta              = new TH1F("hEta", "", 600, -3, 3);              
 hPhi              = new TH1F("hPhi", "", 100, 0, 2*PI);              
 hNClTPC	   = new TH1F("hNClTPC", "", 200, 0, 200);	      
 hNClTPC1	   = new TH1F("hNClTPC1", "", 200, 0, 200);	      
 hNClTPCShared     = new TH1F("hNClTPCShared", "", 200, 0, 200);     
 hNClITS	   = new TH1F("hNClITS", "", 20, 0, 20);	      
 hChi2TPC          = new TH1F("hChi2TPC", "", 100, 0, 100);          
 hChi2ndfTPC       = new TH1F("hChi2ndfTPC", "", 100, 0, 100); // AMA
 hChi2TPC1         = new TH1F("hChi2TPC1", "", 100, 0, 100);         
 hChi2ITS          = new TH1F("hChi2ITS", "", 100, 0, 100);          
 hD                = new TH1F("hD", "", 200, -0.1, 0.1);                
 // hD1               = new TH1F("hD1", "", 200, -0.1, 0.1);      // AMA
 // hD2               = new TH1F("hD2", "", 200, -0.1, 0.1);      // AMA    
 // hD3               = new TH1F("hD3", "", 200, -0.1, 0.1);      // AMA    
 hZ                = new TH1F("hZ", "", 200, -0.1, 0.1);                
 // hZ1               = new TH1F("hZ1", "", 200, -0.1, 0.1);      // AMA    
 // hZ2               = new TH1F("hZ2", "", 200, -0.1, 0.1);      // AMA    
 // hZ3               = new TH1F("hZ3", "", 200, -0.1, 0.1);      // AMA    
 hDTPC             = new TH1F("hDTPC", "", 200, -0.1, 0.1);                
 hZTPC             = new TH1F("hZTPC", "", 200, -0.1, 0.1);                

 // Other trk histos
 hStBits           = new TH1F("hStBits", "", 32, 0, 32); // ULong_t
 hDZ               = new TH2F("hDZ",    "", 100,-0.1,0.1, 100,-0.1,0.1);
 hDZTPC            = new TH2F("hDZTPC", "", 100,-0.1,0.1, 100,-0.1,0.1);

 // hEta1              = new TH1F("hEta1", "", 600, -3, 3);              
 // hEta2              = new TH1F("hEta2", "", 600, -3, 3);              
 // hEta3              = new TH1F("hEta3", "", 600, -3, 3);              

  return;
}

void SetHistProps()
{
  // MyHeader
  //  hRun		    ->SetTitle("Run number;Run;events");
  hBx		    ->SetTitle("Bx;bunch crossing time [ns];events");
  hNChips1	    ->SetTitle("NChips1;chips - SPD layer 1;events");
  hNChips2	    ->SetTitle("NChips2;chips - SPD layer 2;events");
  hNTracks	    ->SetTitle("NTracks;reconstructed tracks;events");
  hNSelTracks	    ->SetTitle("NSelTracks;selected tracks;events");
  hNTracklets	    ->SetTitle("NTracklets;tracklets;events");
  hVx               ->SetTitle("x-vertex;v_{x} [cm];events");
  hVy               ->SetTitle("y-vertex;v_{y} [cm];events");
  hVz               ->SetTitle("z-vertex;v_{z} [cm];events");
  hVc               ->SetTitle("Vertex contributors;tracks;events");
  hIsPileupSPD 	    ->SetTitle("IsPileupSPD;IsPileupSPD;events");
  hNPileupSPD	    ->SetTitle("NPileupSPD;pileup events;events");
  hNPileup	    ->SetTitle("NPileup;pileup events;events");
  hTrCluster	    ->SetTitle("Trigger cluster;trigger cluster;events");
  hVxSPD            ->SetTitle("SPD x-vertex;v_{x} [cm];events");
  hVySPD            ->SetTitle("SPD y-vertex;v_{y} [cm];events");
  hVzSPD            ->SetTitle("SPD z-vertex;v_{z} [cm];events");
  hVcSPD            ->SetTitle("SPD vertex contributors;tracks;events");
  hVxTPC            ->SetTitle("TPC x-vertex;v_{x} [cm];events");
  hVyTPC            ->SetTitle("TPC y-vertex;v_{y} [cm];events");
  hVzTPC            ->SetTitle("TPC z-vertex;v_{z} [cm];events");
  hVcTPC            ->SetTitle("TPC vertex contributors;tracks;events");

  // Other evt histos
  hL0Bits           ->SetTitle("L0 trigger bits;;events");
  hTrClBits         ->SetTitle("Trigger class bits;;events");

  // MyPart
  hC	     	    ->SetTitle("Charge;charge;tracks");
  hPt               ->SetTitle("Pt;p_{T} [GeV/c];tracks");
  hEta              ->SetTitle("Eta;#eta;tracks");
  // hEta1             ->SetTitle("Eta;#eta;tracks");
  // hEta2             ->SetTitle("Eta;#eta;tracks");
  // hEta3             ->SetTitle("Eta;#eta;tracks");
  hPhi              ->SetTitle("Phi;#phi;tracks");
  hNClTPC	    ->SetTitle("TPC clusters in track;TPC clusters;tracks");
  hNClTPC1	    ->SetTitle("TPC clusters in track (iter 1);TPC clusters;tracks");
  hNClTPCShared     ->SetTitle("NClTPCShared;ITS/TPC clusters;tracks");
  hNClITS	    ->SetTitle("ITS clusters;clusters;tracks");
  hChi2TPC          ->SetTitle("TPC #chi^{2};#chi^{2}_{TPC};tracks");
  hChi2TPC1         ->SetTitle("TPC #chi^{2} (iter 1);#chi^{2}_{TPC};tracks");
  hChi2ITS          ->SetTitle("ITS #chi^{2};#chi^{2}_{ITS};tracks");
  hD                ->SetTitle("Transverse DCA;x-y impact parameter [cm];tracks");
  // hD1               ->SetTitle("Transverse DCA;x-y impact parameter [cm];tracks");
  // hD2               ->SetTitle("Transverse DCA;x-y impact parameter [cm];tracks");
  // hD3               ->SetTitle("Transverse DCA;x-y impact parameter [cm];tracks");
  hZ                ->SetTitle("Longitudinal DCA;z impact parameter [cm];tracks");
  // hZ1               ->SetTitle("Longitudinal DCA;z impact parameter [cm];tracks");
  // hZ2               ->SetTitle("Longitudinal DCA;z impact parameter [cm];tracks");
  // hZ3               ->SetTitle("Longitudinal DCA;z impact parameter [cm];tracks");
  hDTPC             ->SetTitle("Transverse TPC DCA;x-y impact parameter [cm];tracks");
  hZTPC             ->SetTitle("Longitudinal TPC DCA;z impact parameter [cm];tracks");

  // Other trk histos
  hStBits           ->SetTitle("Track status flags;;tracks");
  hDZ               ->SetTitle("DCA;x-y impact parameter [cm];z impact parameter [cm]");
  hDZTPC            ->SetTitle("TPC DCA;x-y impact parameter [cm];z impact parameter [cm]");
  hChi2ndfTPC       ->SetTitle("#chi^{2}/N_{TPC clusters}/;#chi^{2};tracks");

  // Put this in TreeClasses.h (and don't change the ordering)!
  ULong_t flagValue[32] = 
    {
     MyPart::kITSin,
     MyPart::kITSout,
     MyPart::kITSrefit,
     MyPart::kITSpid,
     MyPart::kTPCin,
     MyPart::kTPCout,
     MyPart::kTPCrefit,
     MyPart::kTPCpid,
     MyPart::kTRDin,
     MyPart::kTRDout,
     MyPart::kTRDrefit,
     MyPart::kTRDpid,
     MyPart::kTOFin,
     MyPart::kTOFout,
     MyPart::kTOFrefit,
     MyPart::kTOFpid,
     MyPart::kHMPIDout,
     MyPart::kHMPIDpid,
     MyPart::kEMCALmatch,
     MyPart::kTRDbackup,
     0x0,     //20 missing
     MyPart::kPHOSmatch,
     0x0,     // 22-24 missing
     0x0,
     0x0,
     MyPart::kMultInV0,
     MyPart::kMultSec,
     MyPart::kGlobalMerge,
     MyPart::kITSpureSA,
     MyPart::kTRDStop,
     MyPart::kESDpid,
     MyPart::kTIME
    };

  // Put this in TreeClasses.h (and don't change the ordering)!
  TString flagLabel[32] =
    {
      "kITSin",
      "kITSout",
      "kITSrefit",
      "kITSpid",
      "kTPCin",
      "kTPCout",
      "kTPCrefit",
      "kTPCpid",
      "kTRDin",
      "kTRDout",
      "kTRDrefit",
      "kTRDpid",
      "kTOFin",
      "kTOFout",
      "kTOFrefit",
      "kTOFpid",
      "kHMPIDout",
      "kHMPIDpid",
      "kEMCALmatch",
      "kTRDbackup",
      "--",
      "kPHOSmatch",
      "--",
      "--",
      "--",
      "kMultInV0",
      "kMultSec",
      "kGlobalMerge",
      "kITSpureSA",
      "kTRDStop",
      "kESDpid",
      "kTIME"
    };

  if (0) {
    for (int i=0; i<28; i++) {
      int pos = TMath::Log2(flagValue[i]);
      cout << pos << " " << flagLabel[i].Data() << endl;
    }
  }
  
  for (Int_t n=0; n<32; ++n) {
    hStBits->GetXaxis()->SetBinLabel(n+1, flagLabel[n].Data());
  }
  return;
}

void FillEventInfo(const MyHeader& ev)
{
  //  hRun		  ->Fill(ev.fRun		);   
  hBx		  ->Fill(ev.fBx		   	);
  hNChips1	  ->Fill(ev.fNChips1	   	);
  hNChips2	  ->Fill(ev.fNChips2	   	);
  hNTracks	  ->Fill(ev.fNTracks	   	);
  hNSelTracks	  ->Fill(ev.fNSelTracks	   	);
  hNTracklets	  ->Fill(ev.fNTracklets	   	);
  hVx              ->Fill(ev.fVx               	);
  hVy              ->Fill(ev.fVy               	);
  hVz              ->Fill(ev.fVz               	);
  hVc              ->Fill(ev.fVc               	);
  hIsPileupSPD 	  ->Fill(ev.fIsPileupSPD 	);   
  hNPileupSPD	  ->Fill(ev.fNPileupSPD	   	);
  hNPileup	  ->Fill(ev.fNPileup	   	);
  hTrCluster	  ->Fill(ev.fTrCluster	   	);
  hVxSPD           ->Fill(ev.fVxSPD            	);
  hVySPD           ->Fill(ev.fVySPD            	);
  hVzSPD           ->Fill(ev.fVzSPD            	);
  hVcSPD           ->Fill(ev.fVcSPD            	);
  hVxTPC           ->Fill(ev.fVxTPC            	);
  hVyTPC           ->Fill(ev.fVyTPC            	);
  hVzTPC           ->Fill(ev.fVzTPC            	);
  hVcTPC           ->Fill(ev.fVcTPC            	);

  /*
    http://root.cern.ch/root/html/ListOfTypes.html  
    UChar_t Unsigned Character 1 byte (unsigned char)
    UInt_t Unsigned integer 4 bytes (unsigned int)
    ULong64_t Portable unsigned long integer 8 bytes
    ULong_t Unsigned long integer 8 bytes (unsigned long)
  */
  for (UInt_t n=0; n<32; n++) {
    ULong64_t uno = 1;
    Bool_t bit = ev.fL0 & uno<<n;
    hL0Bits->Fill(n+0.5, bit);
  }
  for (ULong64_t n=0; n<64; n++) {
    ULong64_t uno = 1;
    Bool_t bit = ev.fTrClassMask & uno<<n;
    hTrClBits->Fill(n+0.5, bit);
  }

  return;
}

void FillTrackInfo(const MyPart& trk)
{
  // N.B. a pt > 0.4 cut is applied outside this fn.

  // MyPart
  hC	     	 ->Fill(trk.fC);	     	      	  
  hPt            ->Fill(trk.fPt);               
  hEta           ->Fill(trk.fEta);              
  hPhi           ->Fill(trk.fPhi);              
  hNClTPC	 ->Fill(trk.fNClTPC);	      
  hNClTPC1	 ->Fill(trk.fNClTPC1);	      
  hNClTPCShared  ->Fill(trk.fNClTPCShared);     
  hNClITS        ->Fill(trk.fNClITS);	      
  hChi2TPC       ->Fill(trk.fChi2TPC);          
  hChi2TPC1      ->Fill(trk.fChi2TPC1);         
  hChi2ITS       ->Fill(trk.fChi2ITS);          
  hD             ->Fill(trk.fD);                
  hZ             ->Fill(trk.fZ);                
  hDTPC          ->Fill(trk.fDTPC);             
  hZTPC          ->Fill(trk.fZTPC);             

  Double_t chi2tpc = trk.fNClTPC > 0? trk.fChi2TPC/trk.fNClTPC : 999;
  hChi2ndfTPC    ->Fill(chi2tpc);
  hDZ            ->Fill(trk.fD, trk.fZ);                
  hDZTPC         ->Fill(trk.fDTPC, trk.fZTPC);

  // fSt is a ULong_t (8 bytes), but only the first 32 seem to make sense
  for (ULong_t n=0; n<32; ++n) {
    ULong64_t uno = 1;
    Bool_t bit = trk.fSt & uno<<n; // 0 or 1 * 2^n, downcasted
    hStBits->Fill(n+0.5, bit);
  }

  return;
}

TObjArray* GetHistList(TFile& inFile, TString clname)
{
  TObjArray* hList = new TObjArray();
  TIter next(inFile.GetListOfKeys()); 
  TKey *key;

  while ((key=(TKey*)next())) {
    TString className(key->GetClassName());
    TString keyName(key->GetName());
    if (0) printf("%10s %20s\n", className.Data(), keyName.Data());
    
    if (className.Contains(clname) && clname=="TH1F") {
      TH1F* h1 = (TH1F*)inFile.Get(keyName.Data());
      hList->Add(h1);
    }
    if (className.Contains(clname) && clname=="TH2F") {
      TH2F* h2 = (TH2F*)inFile.Get(keyName.Data());
      hList->Add(h2);
    }
  }
  return hList;
}

void DrawHistos(TFile& inFile)
{
  gStyle->SetTitleFontSize(0.06);
  gStyle->SetOptStat(0);

  TObjArray* cList = new TObjArray();
  TObjArray* h1List = GetHistList(inFile, "TH1F");
  TObjArray* h2List = GetHistList(inFile, "TH2F");

  Int_t fillCol=kYellow, lineCol=kBlack, mkrCol=kBlack, mkrSty=kFullCircle;
  Float_t mkrSize=1.0;
  PlotUtils::set_hist_props(h1List,lineCol,fillCol,mkrCol,mkrSty,mkrSize);

  Int_t nCanv = 0, nPrinted = 0;  
  TCanvas* c;

  // TH1Fs
  for (int i=0; i<h1List->GetEntries(); i++) {
    
    TH1F* h1 = (TH1F*)h1List->At(i);
    h1->SetStats(1);
    TString s(h1->GetName());
    TString options(h1->GetDrawOption());

    if (s.Contains("hEta2") || 
	s.Contains("hD2") || 
	s.Contains("hZ2")) {
      options.Append("same");
      h1->SetFillColor(kGreen);
    }
    if (s.Contains("hEta3") || 
	s.Contains("hD3") || 
	s.Contains("hZ3")) {
      options.Append("same");
      h1->SetFillColor(kAzure-9);
    }

    if (s.Contains("hChi2ndfTPC")) {
      h1->SetFillColor(kNone);
      h1->SetLineColor(kRed);
      h1->SetLineWidth(2);
      options.Append("same");
    }
     
    if (!options.Contains("same")) {
      TString ylabel(h1->GetYaxis()->GetTitle());
      TString prefix("c");
      if (ylabel=="events") prefix = "ev";
      if (ylabel=="tracks") prefix = "tr";
      char* title = Form("%s_%s", prefix.Data(), h1->GetName());
      cList->Add(new TCanvas(Form("c%i", nCanv++), title, 1));
    }  
    c = (TCanvas*)cList->Last(); 
    c->cd();
    
    if (s.Contains("hPt") ||
	s.Contains("NClTPC") ||
	s.Contains("Chi2") ) 
      c->SetLogy();
    if (s.Contains("hMult"))
      h1->SetStats(1);
    if (s.Contains("hStBits") || 
	s.Contains("hL0Bits") || 
	s.Contains("hTrClBits")) {
      h1->SetFillColor(kRed);
      h1->SetBarWidth(0.9 * h1->GetBinWidth(1));
      options.Append("hbar");
      c->SetLeftMargin(0.2);
    }
    if (!options.Contains("hbar"))
      PlotUtils::make_nice_axes(c, h1, 1.5);
    if (s.Contains("hVx") || s.Contains("hVy") || 
	s.Contains("hBx") || s.Contains("hTime") ||
	s.Contains("hEvNumber") || s.Contains("hC") || 
	s.Contains("hOrbit") || s=="hD" ||  s=="hZ")  
      h1->GetXaxis()->SetNdivisions(5); // unclutter

    printf("Drawing %20s opt: %s", s.Data(), options.Data());
    cout << endl;
    h1->Draw(options.Data());
    
  }

  // TH2Fs
  for (int i=0; i<h2List->GetEntries(); i++) {
    
    TH2F* h2 = (TH2F*)h2List->At(i);
    TString s(h2->GetName());
    TString options("colz");

    if (s.Contains("hDZ"))  
      h2->GetXaxis()->SetNdivisions(5); // unclutter
    
    cout << "Drawing " << s.Data() << endl;
    

    if (!options.Contains("same")) {
      TString ylabel(h2->GetYaxis()->GetTitle());
      TString prefix("c");
      if (ylabel=="events") prefix = "ev";
      if (ylabel=="tracks") prefix = "tr";
      char* title = Form("%s_%s", prefix.Data(), h2->GetName());
      cList->Add(new TCanvas(Form("c%i", nCanv++), title, 1));
    }  
    c = (TCanvas*)cList->Last(); 
    c->cd();

    PlotUtils::make_nice_axes(c, h2, 1.5);
    
    h2->Draw(options.Data());
  }
  
  if (printPDF) {
    TTimeStamp ts;
    UInt_t date = ts.GetDate();
    for (int i=0; i<cList->GetEntries(); i++) {
      TCanvas* c = (TCanvas*)cList->At(i); 
      c->Print(Form("plots/pngs/%s_%i.png", 
		    c->GetTitle(), date));
      if (nPrinted==0) {
	c->Print(Form("plots/datacheck_%i.ps%s", date, "["));
      }
      c->Print(Form("plots/datacheck_%i.ps%s",         
		    date, i==cList->GetEntries()-1 ? "]" : ""));
      nPrinted++;
    }
    TString base = Form("plots/datacheck_%i", date); 
    TString cmd = "ps2pdf " + base + ".ps " + base + ".pdf";             
    gSystem->Exec(cmd.Data()); 
  }

  return;
}
