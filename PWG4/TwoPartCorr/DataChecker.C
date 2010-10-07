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
#include <TLegend.h>
#include <vector>
#include "TreeClasses.h"
#include "Utils.h"
#endif

Bool_t printPDF = 1;
const Int_t nCuts = 3;
enum eCuts {CUT0=0, CUT1, CUT2}; // may split this into ev and trk
enum {EVT=0, TRK=1};
Int_t fillCol[2][nCuts] = {{kBlue, kYellow, kSpring-5}, 
			   {kAzure+1, kGreen, kRed}};
Int_t lineCol=kBlack, mkrCol=kBlack, mkrSty=kFullCircle;
Float_t mkrSize=1.0;
TLegend* leg[2];
TString cutDef[2][nCuts];
TFile* outFile = new TFile("out.root", "recreate");
const Double_t PI = TMath::Pi();
Int_t nEvents = -1; // Run over whole chain if < 0
MyHeader     *event   = 0;
TClonesArray *tracks  = 0;
TBranch      *hBranch = 0;
TBranch      *tBranch = 0;

// Event
vector<TH1F*> vhBx;		 // crossing times (ns)
vector<TH1F*> vhNChips1;	 
vector<TH1F*> vhNChips2;	 
vector<TH1F*> vhNTracks;	 
vector<TH1F*> vhNSelTracks;       // no kink daughters, NClTPC > 20.
vector<TH1F*> vhNTracklets;	 
vector<TH1F*> vhVx;             
vector<TH1F*> vhVy;             
vector<TH1F*> vhVz;               // from V0 (?)
vector<TH1F*> vhVc;               // from TPC (?)
vector<TH1F*> vhIsPileupSPD; 	 
vector<TH1F*> vhNPileupSPD;	 
vector<TH1F*> vhNPileup;	 
vector<TH1F*> vhTrClassMask;	 // get names from constantin 
vector<TH1F*> vhTrCluster;	 
vector<TH1F*> vhVxSPD;          
vector<TH1F*> vhVySPD;          
vector<TH1F*> vhVzSPD;          
vector<TH1F*> vhVcSPD;          
vector<TH1F*> vhVxTPC;          
vector<TH1F*> vhVyTPC;          
vector<TH1F*> vhVzTPC;          
vector<TH1F*> vhVcTPC;                     

// Track
vector<TH1F*> vhSt;	     // status
vector<TH1F*> vhC;	     // charge
vector<TH1F*> vhPt;           
vector<TH1F*> vhEta;          
vector<TH1F*> vhPhi;          
vector<TH1F*> vhNClTPC;	     
vector<TH1F*> vhNClTPC1;	     // ?
vector<TH1F*> vhNClTPCShared; 
vector<TH1F*> vhNClITS;	     
vector<TH1F*> vhChi2TPC;      
vector<TH1F*> vhChi2TPC1;     
vector<TH1F*> vhChi2ITS;      
vector<TH1F*> vhD;            // Transverse DCA or "impact parameter"
vector<TH1F*> vhZ;            // Z DCA
vector<TH1F*> vhDTPC;         
vector<TH1F*> vhZTPC;           

vector<TH1F*> vhStBits;
vector<TH1F*> vhL0Bits;
vector<TH1F*> vhTrClBits;
vector<TH1F*> vhChi2ndfTPC;   
vector<TH1F*> vhChi2ndfITS; // <---- Still need to make this  

vector<TH2F*> vhDZ;
vector<TH2F*> vhDZTPC;

void BookHistos();
TObjArray* GetHistList(TFile& inFile, TString clname);
void SetHistProps();
void DrawHistos(TFile& inFile);
void FillEventInfo(const MyHeader& ev, const Int_t k);
void FillTrackInfo(const MyPart& trk,  const Int_t k);

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
      // get bits from maps here and fill array of names
#if 0
TFile f("rootfiles/merged_run119037.root")
f.ls()
TrigClass_run119037_map->ls()
TrigClass_run119037_map->FindObject("CINT1B-ABCE-NOPF-ALL")
TPair *p = TrigClass_run119037_map->FindObject("CINT1B-ABCE-NOPF-ALL")
p->Value()->GetName()
#endif
    }
    
    hBranch->GetEntry(li);
    tBranch->GetEntry(li);

     Bool_t evCut0=1, evCut1=1, evCut2=1;
     Bool_t trCut0=1, trCut1=1, trCut2=1;
     evCut0 = 1;
     evCut1 = !event->fIsPileupSPD && fabs(event->fVz) < 10;
     evCut2 = event->fVc > 15;

    if (evCut0) {
      FillEventInfo(*event, CUT0);
      if (evCut1) {
	FillEventInfo(*event, CUT1);
	if (evCut2) {
	  FillEventInfo(*event, CUT2);
	}
      }
    }
    
    // --- Track loop ---
    Int_t ntracks = tracks->GetEntries();
    for (Int_t t=0;t<ntracks;++t) {
      MyPart *trk = (MyPart*)tracks->At(t);

     trCut0 = evCut1;
     trCut1 = trk->fNClTPC > 0;
     trCut2 = trk->fNClTPC > 70;
      
      if (trCut0) {
	FillTrackInfo(*trk, CUT0);
	if (trCut1) {
	  FillTrackInfo(*trk, CUT1);
	  if (trCut2) {
	    FillTrackInfo(*trk, CUT2);
	  }
	}
      }
    }
  }
  
  outFile->Write();
  SetHistProps();
  DrawHistos(*outFile);
  
  return;
}

void BookHistos()
{
  for (Int_t k=0; k<nCuts; k++) {
    // EVENT variables
    vhBx          .push_back(new TH1F(Form("EVThBx_%d",k),"", 3000, 0, 3000));
    vhNChips1	  .push_back(new TH1F(Form("EVThNChips1_%d",k),"", 150, 0, 150));	      
    vhNChips2	  .push_back(new TH1F(Form("EVThNChips2_%d",k),"", 150, 0, 150));	      
    vhNTracks	  .push_back(new TH1F(Form("EVThNTracks_%d",k),"", 500, 0, 500));	      
    vhNSelTracks  .push_back(new TH1F(Form("EVThNSelTracks_%d",k),"", 500, 0, 500));
    vhNTracklets  .push_back(new TH1F(Form("EVThNTracklets_%d",k),"", 100, 0, 100));
    vhVx          .push_back(new TH1F(Form("EVThVx_%d",k),"", 100, -0.1, 0.0));               
    vhVy          .push_back(new TH1F(Form("EVThVy_%d",k),"", 100, 0.15, 0.25));               
    vhVz          .push_back(new TH1F(Form("EVThVz_%d",k),"", 40, -20, 20));               
    vhVc          .push_back(new TH1F(Form("EVThVc_%d",k),"", 100, 0, 100));               
    vhIsPileupSPD .push_back(new TH1F(Form("EVThIsPileupSPD_%d",k),"", 3, 0, 3));  
    vhNPileupSPD  .push_back(new TH1F(Form("EVThNPileupSPD_%d",k),"", 10, 0, 10));
    vhNPileup	  .push_back(new TH1F(Form("EVThNPileup_%d",k),"", 10, 0, 10));	      
    vhTrCluster	  .push_back(new TH1F(Form("EVThTrCluster_%d",k),"", 10, 0, 10));
    vhVxSPD       .push_back(new TH1F(Form("EVThVxSPD_%d",k),"", 100, -0.5, 0.5));            
    vhVySPD       .push_back(new TH1F(Form("EVThVySPD_%d",k),"", 100, -0.5, 0.5));            
    vhVzSPD       .push_back(new TH1F(Form("EVThVzSPD_%d",k),"", 100, -20, 20));            
    vhVcSPD       .push_back(new TH1F(Form("EVThVcSPD_%d",k),"", 100, 0, 100));            
    vhVxTPC       .push_back(new TH1F(Form("EVThVxTPC_%d",k),"", 100, -0.01, 0.01));
    vhVyTPC       .push_back(new TH1F(Form("EVThVyTPC_%d",k),"", 100, -0.01, 0.01));
    vhVzTPC       .push_back(new TH1F(Form("EVThVzTPC_%d",k),"", 100, -5, 5));            
    vhVcTPC       .push_back(new TH1F(Form("EVThVcTPC_%d",k),"", 100, -1, 1));            
    vhL0Bits      .push_back(new TH1F(Form("EVThL0Bits_%d",k),"", 32, 0, 32));
    vhTrClBits    .push_back(new TH1F(Form("EVThTrClBits_%d",k),"", 64, 0, 64));
    
    // TRACK variables
    vhC	     	  .push_back(new TH1F(Form("TRKhC_%d",k),"", 5, -2.5, 2.5));
    vhPt          .push_back(new TH1F(Form("TRKhPt_%d",k),"", 200, 0, 20));               
    vhEta         .push_back(new TH1F(Form("TRKhEta_%d",k),"", 600, -3, 3));              
    vhPhi         .push_back(new TH1F(Form("TRKhPhi_%d",k),"", 100, -PI/6, 2*PI+PI/6));          
    vhNClTPC	  .push_back(new TH1F(Form("TRKhNClTPC_%d",k),"", 200, 0, 200));	      
    vhNClTPC1	  .push_back(new TH1F(Form("TRKhNClTPC1_%d",k),"", 200, 0, 200));	      
    vhNClTPCShared.push_back(new TH1F(Form("TRKhNClTPCShared_%d",k),"", 200, 0, 200));     
    vhNClITS	  .push_back(new TH1F(Form("TRKhNClITS_%d",k),"", 20, 0, 20));	      
    vhChi2TPC     .push_back(new TH1F(Form("TRKhChi2TPC_%d",k),"", 100, 0, 100));          
    vhChi2ndfTPC  .push_back(new TH1F(Form("TRKhChi2ndfTPC_%d",k),"", 100, 0, 100));
    vhChi2TPC1    .push_back(new TH1F(Form("TRKhChi2TPC1_%d",k),"", 100, 0, 100));         
    vhChi2ITS     .push_back(new TH1F(Form("TRKhChi2ITS_%d",k),"", 100, 0, 100));          
    vhD           .push_back(new TH1F(Form("TRKhD_%d",k),"", 200, -0.1, 0.1));                
    vhZ           .push_back(new TH1F(Form("TRKhZ_%d",k),"", 200, -0.1, 0.1));                
    vhDTPC        .push_back(new TH1F(Form("TRKhDTPC_%d",k),"", 200, -0.1, 0.1));
    vhZTPC        .push_back(new TH1F(Form("TRKhZTPC_%d",k),"", 200, -0.1, 0.1));
    vhStBits      .push_back(new TH1F(Form("TRKhStBits_%d",k),"", 32, 0, 32));
    vhDZ          .push_back(new TH2F(Form("TRKhDZ_%d",k),"", 100,-0.1,0.1, 100,-0.1,0.1));
    vhDZTPC       .push_back(new TH2F(Form("TRKhDZTPC_%d",k),"", 100,-0.1,0.1, 100,-0.1,0.1));
  }

 return;
}

void SetHistProps()
{
  for (Int_t k=0; k<nCuts; k++) {
    vhBx	  .at(k)->SetTitle("Bx;bunch crossing time [ns];events");
    vhNChips1	  .at(k)->SetTitle("NChips1;chips - SPD layer 1;events");
    vhNChips2	  .at(k)->SetTitle("NChips2;chips - SPD layer 2;events");
    vhNTracks	  .at(k)->SetTitle("NTracks;reconstructed tracks;events");
    vhNSelTracks  .at(k)->SetTitle("NSelTracks;selected tracks;events");
    vhNTracklets  .at(k)->SetTitle("NTracklets;tracklets;events");
    vhVx          .at(k)->SetTitle("x-vertex;v_{x} [cm];events");
    vhVy          .at(k)->SetTitle("y-vertex;v_{y} [cm];events");
    vhVz          .at(k)->SetTitle("z-vertex;v_{z} [cm];events");
    vhVc          .at(k)->SetTitle("Vertex contributors;tracks;events");
    vhIsPileupSPD .at(k)->SetTitle("IsPileupSPD;IsPileupSPD;events");
    vhNPileupSPD  .at(k)->SetTitle("NPileupSPD;pileup events;events");
    vhNPileup	  .at(k)->SetTitle("NPileup;pileup events;events");
    vhTrCluster	  .at(k)->SetTitle("Trigger cluster;trigger cluster;events");
    vhVxSPD       .at(k)->SetTitle("SPD x-vertex;v_{x} [cm];events");
    vhVySPD       .at(k)->SetTitle("SPD y-vertex;v_{y} [cm];events");
    vhVzSPD       .at(k)->SetTitle("SPD z-vertex;v_{z} [cm];events");
    vhVcSPD       .at(k)->SetTitle("SPD vertex contributors;tracks;events");
    vhVxTPC       .at(k)->SetTitle("TPC x-vertex;v_{x} [cm];events");
    vhVyTPC       .at(k)->SetTitle("TPC y-vertex;v_{y} [cm];events");
    vhVzTPC       .at(k)->SetTitle("TPC z-vertex;v_{z} [cm];events");
    vhVcTPC       .at(k)->SetTitle("TPC vertex contributors;tracks;events");
    vhL0Bits      .at(k)->SetTitle("L0 trigger bits;;events");
    vhTrClBits    .at(k)->SetTitle("Trigger class bits;;events");
    vhC           .at(k)->SetTitle("Charge;charge;tracks");
    vhPt          .at(k)->SetTitle("Pt;p_{T} [GeV/c];tracks");
    vhEta         .at(k)->SetTitle("Eta;#eta;tracks");
    vhPhi         .at(k)->SetTitle("Phi;#phi;tracks");
    vhNClTPC      .at(k)->SetTitle("TPC clusters in track;TPC clusters;tracks");
    vhNClTPC1     .at(k)->SetTitle("TPC clusters in track (iter 1);TPC clusters;tracks");
    vhNClTPCShared.at(k)->SetTitle("NClTPCShared;ITS/TPC clusters;tracks");
    vhNClITS      .at(k)->SetTitle("ITS clusters;clusters;tracks");
    vhChi2TPC     .at(k)->SetTitle("TPC #chi^{2};total #chi^{2}_{TPC};tracks");
    vhChi2TPC1    .at(k)->SetTitle("TPC #chi^{2} (iter 1);total #chi^{2}_{TPC};tracks");
    vhChi2ITS     .at(k)->SetTitle("ITS #chi^{2};#chi^{2}_{ITS};tracks");
    vhD           .at(k)->SetTitle("Transverse DCA;x-y impact parameter [cm];tracks");
    vhZ           .at(k)->SetTitle("Longitudinal DCA;z impact parameter [cm];tracks");
    vhDTPC        .at(k)->SetTitle("Transverse TPC DCA;x-y impact parameter [cm];tracks");
    vhZTPC        .at(k)->SetTitle("Longitudinal TPC DCA;z impact parameter [cm];tracks");
    vhStBits      .at(k)->SetTitle("Track status flags;;tracks");
    vhDZ          .at(k)->SetTitle("Impact parameter;x-y DCA [cm];z DCA [cm]");
    vhDZTPC       .at(k)->SetTitle("TPC Impact parameter;x-y DCA [cm];z DCA [cm]");
    vhChi2ndfTPC  .at(k)->SetTitle("#chi^{2} / N_{TPC clusters};#chi^{2}_{TPC} / N;tracks");
    
    for (Int_t n=0; n<32; ++n) {
      vhStBits.at(k)->GetXaxis()->SetBinLabel(n+1, flagLabel[n].Data());
    }
  }
  
  return;
}

void FillEventInfo(const MyHeader& ev, const Int_t k)
{
  if (k>=nCuts) {
    cout << "ERROR in FillEventInfo(): bad cut index " 
	 << k << endl;
    return;
  }
  vhBx         .at(k)->Fill(ev.fBx	   );
  vhNChips1    .at(k)->Fill(ev.fNChips1    );
  vhNChips2    .at(k)->Fill(ev.fNChips2    );
  vhNTracks    .at(k)->Fill(ev.fNTracks    );
  vhNSelTracks .at(k)->Fill(ev.fNSelTracks );
  vhNTracklets .at(k)->Fill(ev.fNTracklets );
  vhVx         .at(k)->Fill(ev.fVx         );
  vhVy         .at(k)->Fill(ev.fVy         );
  vhVz         .at(k)->Fill(ev.fVz         );
  vhVc         .at(k)->Fill(ev.fVc         );
  vhIsPileupSPD.at(k)->Fill(ev.fIsPileupSPD);   
  vhNPileupSPD .at(k)->Fill(ev.fNPileupSPD );
  vhNPileup    .at(k)->Fill(ev.fNPileup    );
  vhTrCluster  .at(k)->Fill(ev.fTrCluster  );
  vhVxSPD      .at(k)->Fill(ev.fVxSPD      );
  vhVySPD      .at(k)->Fill(ev.fVySPD      );
  vhVzSPD      .at(k)->Fill(ev.fVzSPD      );
  vhVcSPD      .at(k)->Fill(ev.fVcSPD      );
  vhVxTPC      .at(k)->Fill(ev.fVxTPC      );
  vhVyTPC      .at(k)->Fill(ev.fVyTPC      );
  vhVzTPC      .at(k)->Fill(ev.fVzTPC      );
  vhVcTPC      .at(k)->Fill(ev.fVcTPC      );

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
    vhL0Bits.at(k)->Fill(n+0.5, bit);
  }
  for (ULong64_t n=0; n<64; n++) {
    ULong64_t uno = 1;
    Bool_t bit = ev.fTrClassMask & uno<<n;
    vhTrClBits.at(k)->Fill(n+0.5, bit);
  }
  return;
}

void FillTrackInfo(const MyPart& trk, const Int_t k)
{
  // N.B. a pt > 0.4 cut is applied outside this fn.
  if (k>=nCuts) {
    cout << "ERROR in FillTrackInfo(): bad cut index " 
	 << k << endl;
    return;
  }

  Double_t chi2tpc = trk.fNClTPC > 0? trk.fChi2TPC/trk.fNClTPC : 999;
  vhC	     	.at(k)->Fill(trk.fC);	     	      	  
  vhPt          .at(k)->Fill(trk.fPt);               
  vhEta         .at(k)->Fill(trk.fEta);              
  vhPhi         .at(k)->Fill(trk.fPhi);              
  vhNClTPC	.at(k)->Fill(trk.fNClTPC);	      
  vhNClTPC1	.at(k)->Fill(trk.fNClTPC1);	      
  vhNClTPCShared.at(k)->Fill(trk.fNClTPCShared);     
  vhNClITS      .at(k)->Fill(trk.fNClITS);	      
  vhChi2TPC     .at(k)->Fill(trk.fChi2TPC);          
  vhChi2TPC1    .at(k)->Fill(trk.fChi2TPC1);         
  vhChi2ITS     .at(k)->Fill(trk.fChi2ITS);          
  vhD           .at(k)->Fill(trk.fD);                
  vhZ           .at(k)->Fill(trk.fZ);                
  vhDTPC        .at(k)->Fill(trk.fDTPC);             
  vhZTPC        .at(k)->Fill(trk.fZTPC);             
  vhChi2ndfTPC  .at(k)->Fill(chi2tpc);
  vhDZ          .at(k)->Fill(trk.fD, trk.fZ);                
  vhDZTPC       .at(k)->Fill(trk.fDTPC, trk.fZTPC);

  // fSt is a ULong_t (8 bytes), but only the first 32 seem to make sense
  for (ULong_t n=0; n<32; ++n) {
    ULong64_t uno = 1;
    Bool_t bit = trk.fSt & uno<<n; // 0 or 1 * 2^n, downcasted
    vhStBits.at(k)->Fill(n+0.5, bit);
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
      // Transpose from how histos are listed in file
      if (keyName.EndsWith("_0")) {
	hList->Add((TH1F*)inFile.Get(keyName.Data()));
	for (Int_t k=1; k<nCuts; k++) {
	  TString name(keyName.Replace(keyName.Length()-1, 1, Form("%d",k)));
	  //	  printf("%20s %20s\n", keyName.Data(), name.Data());
	  hList->Add((TH1F*)inFile.Get(name.Data()));
	}
      }
    }
    if (className.Contains(clname) && clname=="TH2F") {
      // Transpose from how histos are listed in file
      if (keyName.EndsWith("_0")) {
	hList->Add((TH2F*)inFile.Get(keyName.Data()));
	for (Int_t k=1; k<nCuts; k++) {
	  TString name(keyName.Replace(keyName.Length()-1, 1, Form("%d",k)));
	  //	  printf("%20s %20s\n", keyName.Data(), name.Data());
	  hList->Add((TH2F*)inFile.Get(name.Data()));
	}
      }
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

  PlotUtils::set_hist_props(h1List,lineCol,kNone,mkrCol,mkrSty,mkrSize);

  Int_t nCanv = 0, nPrinted = 0;  
  TCanvas* c;

  for (int i=0; i<2; i++) {
    leg[i] = new TLegend(0.6, 0.85, 0.99, 0.99);
    leg[i]->SetFillColor(kNone);
    leg[i]->SetBorderSize(0);
  }
  cutDef[EVT][0] = "no cuts";
  cutDef[EVT][1] = "!fIsPileupSPD && fabs(fVz) < 10";
  cutDef[EVT][2] = "fVc > 15";
  cutDef[TRK][0] = "!fIsPileupSPD && fabs(fVz) < 10";
  cutDef[TRK][1] = "fNClTPC > 0";
  cutDef[TRK][2] = "fNClTPC > 70";

  // TH1Fs
  for (int i=0; i<h1List->GetEntries(); i++) {
    
    TH1F* h1 = (TH1F*)h1List->At(i);
    h1->SetStats(1);
    TString s(h1->GetName());
    TString options(h1->GetDrawOption());
    
    Int_t T = s.Contains("EVT")? 0 : 1;

    int nLegEntries = leg[T]->GetListOfPrimitives()->GetEntries();

    if (s.Contains("_0")) {
      h1->SetFillColor(fillCol[T][CUT0]);
      if (nLegEntries < nCuts)
	leg[T]->AddEntry(h1, cutDef[T][CUT0].Data(), "f");
    }
    if (s.Contains("_1")) {
      h1->SetFillColor(fillCol[T][CUT1]);
      options.Append("same");
      if (nLegEntries < nCuts)
	leg[T]->AddEntry(h1, cutDef[T][CUT1].Data(), "f");
    }
    if (s.Contains("_2")) {
      h1->SetFillColor(fillCol[T][CUT2]);
      options.Append("same");
      if (nLegEntries < nCuts)
	leg[T]->AddEntry(h1, cutDef[T][CUT2].Data(), "f");
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
      //      h1->SetFillColor(kRed);
      h1->SetBarWidth(0.9 * h1->GetBinWidth(1));
      options.Append("hbar");
      c->SetLeftMargin(0.2);
    }
    if (!options.Contains("hbar"))
      PlotUtils::make_nice_axes(c, h1, 1.5);
    if (s.Contains("hVx") || s.Contains("hVy") || 
	s.Contains("hBx") || s.Contains("hC") || 
	s.Contains("hOrbit") || s.Contains("hD") ||  s.Contains("hZ") ) { 
      h1->GetXaxis()->SetNdivisions(5); // unclutter

      // if (s.Contains("_0")) 
      // 	h1->GetYaxis()->SetRangeUser(0.0, h1->GetYaxis()->GetXmax());
    }

    printf("Drawing %20s opt: %s", s.Data(), options.Data());
    cout << endl;
    h1->Draw(options.Data());
    leg[T]->Draw();
    
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
