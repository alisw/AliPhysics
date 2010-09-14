// $Id:$

#if !defined(__CINT__) || defined(__MAKECINT__)
#include "TreeClasses.h"
#include "EventPool.h"
#include "AutoCorr.h"
#endif

bool debug = true;

void runAutoCorr()
{
  const double PI = TMath::Pi();

  //TFile* f = new TFile("output.root", "read");
  //TList* l = (TList*)f->Get("output");
  //TTree* tree = (TTree*)l->At(11); 
  TChain *tree = new TChain("MyTree");
  tree->Add("res_LHC10e_09132010/mergedruns/merged_*.root");
  cout << "Entries " << tree->GetEntries() << endl;

  MyHeader* ev = 0;//new MyHeader();
  TClonesArray* trk = 0;//new TClonesArray("MyPart");
  AutoCorr* ac = new AutoCorr();

  TBranch* evBranch = tree->GetBranch("header");
  evBranch->SetAddress(&ev);
  
  TBranch* trBranch = tree->GetBranch("parts");
  trBranch->SetAddress(&trk);

  ////////////////////////
  // TEST DRIVE
  ////////////////////////
  Int_t poolsize = 20;
  const Int_t nMultBins = 3;
  Double_t multbin[nMultBins+1] = {180,230, 280, 500}; 
  const Int_t nZvtxBins = 3;
  Double_t zvtxbin[nZvtxBins+1] = {-6, -2, 2, 6};
  ac->InitEventPools(poolsize, nMultBins, multbin, nZvtxBins, zvtxbin);

  // TODO: encapsulate these in AutoCorr (pass arrays above to autocorr ctor)
  TH1F* hMultBins = new TH1F("hMultBins", "Event multiplicity binning", 
			     nMultBins, multbin);
  TH1F* hZvtxBins = new TH1F("hZvtxBins", "Event Z-vertex binning", 
			     nZvtxBins, zvtxbin);

  if (debug) {
    cout << "Mult. binning: ";
  for (int j=1; j<=hMultBins->GetNbinsX(); j++) 
    cout << hMultBins->GetBinLowEdge(j) << " ";
  cout << hMultBins->GetBinLowEdge(nMultBins) + 
    hMultBins->GetBinWidth(nMultBins) << " ";
  cout << endl;
    cout << "Z-vtx. binning: ";
  for (int j=1; j<=hZvtxBins->GetNbinsX(); j++) 
    cout << hZvtxBins->GetBinLowEdge(j) << " ";
  cout << hZvtxBins->GetBinLowEdge(nZvtxBins) + 
    hZvtxBins->GetBinWidth(nZvtxBins) << " ";
  cout << endl;
  }

  // More test examples  
  Short_t zmin = -5, zmax = 5;
  Int_t multmin = 60, multmax = 100;
  EventPool* testPool = new EventPool(poolsize); 
  EventPool* binnedPool = new EventPool(poolsize, multmin, multmax, zmin, zmax); 


  // TODO: Put histos in a manager
  // ---------------------------------
  Double_t etaMax = 2.2;

  TH2F* hSig[nMultBins];
  TH2F* hBkg[nMultBins];
  TH2F* hSB[nMultBins];
  TH1F* h1 = new TH1F("h1", "h1", 1000, -2, 2);
  TH1F* h2 = new TH1F("h2", "h2", 1000, -2, 2);
  TH1F* h3 = new TH1F("h3", "h3", 100, -PI/2, 3*PI/2);
  
  for (int iM=0; iM<nMultBins; iM++) {
    char* name;
    name = Form("hSig_%i", iM);
    hSig[iM] = new TH2F(name, name, 100, -PI/2, 3*PI/2, 100, -etaMax, etaMax);
    name = Form("hBkg_%i", iM);
    hBkg[iM] = new TH2F(name, name, 100, -PI/2, 3*PI/2, 100, -etaMax, etaMax);
    name = Form("hSB_%i", iM);
    hSB[iM] = new TH2F(name, name, 100, -PI/2, 3*PI/2, 100, -etaMax, etaMax);
  }
  
  TH2F* hDR = new TH2F("hDR", "hDR", 100, -0.01, 0.01, 100, -0.005, 0.005);
  hDR->SetTitle("#Delta#phi-#Delta#eta Distribution (no cut);#Delta#phi;#Delta#eta");
  TH2F* hDRcut = new TH2F("hDRcut", "hDRcut", 200, -0.2, 0.2, 200, -0.1, 0.1);
  hDRcut->SetTitle("#Delta#phi-#Delta#eta Distribution;#Delta#phi;#Delta#eta");
  h2->SetLineColor(kRed);
  h3->SetLineColor(kBlue);
  // ---------------------------------

  Int_t nEvents = tree->GetEntries();
  //nEvents = 10000;
  for (int n=0; n<nEvents; n++) {

    // GetEntry fills "header" and "parts" branches, reinitializing ev
    // and trk instances for each event.
    //tree->GetEntry(n); 
    evBranch->GetEntry(n);
    trBranch->GetEntry(n);

//    trk->Print("ls");
    if (n % 100 == 0) cout << n << endl;

    // tests...don't need
    testPool->UpdatePool(n, ev, trk);
    binnedPool->UpdatePool(n, ev, trk);

    ac->UpdatePools(n, ev, trk);

    int SignalEventMultBin = hMultBins->FindBin(ev->fNSelTracks) - 1;
    int SignalEventZvtxBin = hZvtxBins->FindBin(ev->fVz) - 1;
    EventPool* pool = ac->GetEventPool(SignalEventMultBin,
				       SignalEventZvtxBin);

    if (!pool) continue;     // No mixing pool in range specified by arrays
    if (!pool->IsPoolReady()) continue; // TODO: Figure out something better

    for (int i=0; i<ev->fNSelTracks-1; i++) {
      MyPart* t1 = (MyPart*)trk->At(i);
      if (!ac->IsTrackOk(t1)) continue;
      
      h1->Fill(t1->Eta());

      // --- Same-event correlation loop ---
      for (int j=i+1; j<ev->fNSelTracks; j++) {
	MyPart* t2 = (MyPart*)trk->At(j);
	Double_t deta = ac->DeltaEta(t1, t2);
	Double_t dphi = ac->DeltaPhi(t1, t2);
	hDR->Fill(dphi, deta); // No cut - contamination pk.

	if (ac->IsPairOk(t1, t2)) {
	  h2->Fill(deta);
	  h3->Fill(dphi);
	  hSig[SignalEventMultBin]->Fill(dphi, deta);
	  hDRcut->Fill(dphi, deta); // Show elliptical hole from cut 
	}
      }
      // ---

      // --- Event mixing loop---
      // TODO: change 100 to f*<mult> of bkg. event where f is some fixed fraction
      for (int j=i+1; j<0.2*ev->fNSelTracks; j++) { // TODO: update track selection
	MyPart* tbg = pool->GetRandomTrack();
	if (ac->IsMixedPairOk(t1, tbg)) {
	  Double_t deta_bg = ac->DeltaEta(t1, tbg);
	  Double_t dphi_bg = ac->DeltaPhi(t1, tbg);
	  hBkg[SignalEventMultBin]->Fill(dphi_bg, deta_bg);
	}
      }
      // ---      

    }
    
  }
  TCanvas* c1 = new TCanvas("c1", "c1", 1);
  c1->cd();
  h2->Draw();
  h1->Draw("same");

  TCanvas* c2 = new TCanvas("c2", "c2", 1);
  c2->cd();
  h3->Draw();

  TCanvas* c4 = new TCanvas("c4", "c4", 1);
  c4->Divide(2, 1, 0.001, 0.001);
  c4->cd(1);
  hDR->Draw("surf2");
  c4->cd(2);
  hDRcut->Draw("colz");

  TCanvas* cSig = new TCanvas("cSig", "cSig", 1);
  cSig->Divide(nMultBins, 1, 0.001, 0.001);
  for (int iM=0; iM<nMultBins; iM++) {
    cSig->cd(iM+1);
    hSig[iM]->Draw("surf1");
  }
  TCanvas* cBkg = new TCanvas("cBkg", "cBkg", 1);
  cBkg->Divide(nMultBins, 1, 0.001, 0.001);
  for (int iM=0; iM<nMultBins; iM++) {
    cBkg->cd(iM+1);
    hBkg[iM]->Draw("surf1");
  }
  
  TCanvas* cSB = new TCanvas("cSB", "cSB", 1);
  cSB->Divide(nMultBins, 1, 0.001, 0.001);
  for (int iM=0; iM<nMultBins; iM++) {
    cSB->cd(iM+1);
    hSB[iM]->Divide(hSig[iM], hBkg[iM]);
    hSB[iM]->Draw("surf1");
  }  
  return;
}

