// $Id$

// To run this:
// .x rootlogon.C
// .x runAutoCorr.C+

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TStopwatch.h>
#include <TString.h>
#include <TROOT.h>
#include "TreeClasses.h"
#include "EventPool.h"
#include "AutoCorr.h"
#endif

bool debug = true;

void runAutoCorr(const char* dataFile = 
		 "../rootfiles/res_LHC10e_09122010/mergedruns/merged_run130601.root", 
		 const char* outFileName = "output/130601.root")
{
  TFile* outFile = new TFile(outFileName, "recreate");
  const double PI = TMath::Pi();
  Int_t poolsize = 50;
  Int_t nMix = 5;
  const Int_t nMultBins = 2;
  Double_t multbin[nMultBins+1] = {100, 200, 400}; 
  const Int_t nZvtxBins = 2;
  Double_t zvtxbin[nZvtxBins+1] = {-15, 0, 15};
  const Int_t nPtBins = 3;
  Double_t ptbin[nPtBins+1] = {0.15, 0.50, 2.0, 10.};
  Double_t etaMin = -1.4; // for track cuts
  Double_t etaMax = 1.4;
  Double_t ptMin = 0.150; // for track cuts
  Double_t ptMax = 10.0;
  Double_t dEtaMax = 2.5; // for histogram limits
  TString sDataSet = "res_LHC10e_09122010";

  // Event cuts
  Int_t minVc = 25;
  Int_t maxNTracklets = 200;
  Double_t zMin = -15.0;
  Double_t zMax = 15.0;
  Int_t runNumber = -1;
  Int_t firstRunNumber = -1;

  TString sMultBins(""), sZvtxBins(""), sPtBins("");
  for (int i=0; i<nMultBins+1; i++) sMultBins.Append(Form("%d ", (Int_t)multbin[i]));
  for (int i=0; i<nZvtxBins+1; i++) sZvtxBins.Append(Form("%d ", (Int_t)zvtxbin[i]));
  for (int i=0; i<nPtBins+1;   i++)   sPtBins.Append(Form("%d ", (Int_t)ptbin[i]));

  // Document the applied cuts, get quick look with cuts->Print()
  TList* cuts = new TList();
  cuts->Add(new TNamed("data_sample", sDataSet.Data()));
  cuts->Add(new TNamed("mult_bins", sMultBins.Data()));
  cuts->Add(new TNamed("zvtx_bins", sZvtxBins.Data()));
  cuts->Add(new TNamed("mult_bins", sPtBins.Data()));
  cuts->Add(new TNamed("mass_cut", "No mass cut"));
  cuts->Add(new TNamed("eta_range", Form("%.1f < #eta < %.1f", etaMin, etaMax)));
  cuts->Add(new TNamed("z_vtx_range", Form("%.1f < z-vtx < %.1f", zMin, zMax)));
  cuts->Add(new TNamed("min_vc_ntracks", Form("No. vtx. contributors > %i", minVc)));
  cuts->Add(new TNamed("max_ntracklets", Form("No. trackets > %i", maxNTracklets)));
  cuts->Add(new TNamed("pool_size", Form("Event pool size = %i", poolsize )));
  cuts->Add(new TNamed("ntracks_mix", Form("Mixed tracks per real track = %i", nMix )));

  // Histos to hold binning info
  TH1F* hMultBins = new TH1F("hMultBins", "Event multiplicity binning", 
			     nMultBins, multbin);
  TH1F* hZvtxBins = new TH1F("hZvtxBins", "Event Z-vertex binning", 
			     nZvtxBins, zvtxbin);
  TH1F* hPtBins = new TH1F("hPtBins", "p_{T} binning", 
			     nPtBins, ptbin);

  TChain *tree = new TChain("MyTree");
  tree->Add(dataFile);
  cout << "Found " << tree->GetEntries() << " entries in tree!" << endl;

  MyHeader* ev = 0;
  TClonesArray* trk = 0;

  TBranch* evBranch = tree->GetBranch("header");
  evBranch->SetAddress(&ev);
  
  TBranch* trBranch = tree->GetBranch("parts");
  trBranch->SetAddress(&trk);

  AutoCorr* ac = new AutoCorr();
  ac->InitEventPools(poolsize, nMultBins, multbin, nZvtxBins, zvtxbin);

  if (debug) {
    cout << nMultBins << " x " << nZvtxBins 
	 << " event pool(s) to initialize." << endl;
    cout << "Mult. binning: ";
    for (int j=0; j<=nMultBins; j++) cout <<  multbin[j] << " ";
    cout << endl;
    cout << "Zvtx. binning: ";
    for (int j=0; j<=nZvtxBins; j++) cout <<  zvtxbin[j] << " ";
    cout << endl;
  }

  TH2F* hSig[nMultBins];
  TH2F* hBkg[nMultBins];
  TH2F* hSB[nMultBins];
  TH2F* hSigPt[nMultBins][nPtBins];
  TH2F* hBkgPt[nMultBins][nPtBins];
  TH2F* hSBPt[nMultBins][nPtBins];

  TH1F* hMult[nMultBins];
  TH1F* hEta = new TH1F("hEta", "hEta", 600, -3, 3);
  TH1F* hDEta = new TH1F("hDEta", "hDEta", 600, -3, 3);
  TH1F* hDPhi = new TH1F("hDPhi", "hDPhi", 60, -PI/2, 3*PI/2);
  TH1F* hMultAll = new TH1F("hMultAll", "hMultAll", 1000, 0, 1000);
  TH1F* hMultSel = new TH1F("hMultSel", "hMultSel", 1000, 0, 1000);
  TH1F* hPtAll = new TH1F("hPtAll", "hPtAll", 100, 0., 10.);
  TH1F* hPtSel = new TH1F("hPtSel", "hPtSel", 100, 0., 10.);
  TH1F* hVc = new TH1F("hVc", "Vertex contributor tracks", 1000, 0, 1000);
  TH1F* hNTracklets = new TH1F("hNTracklets", "Tracklets", 1000, 0, 1000);
  TH1F* hZvtx = new TH1F("hZvtx", "Event Z-vertex", 60, -30, 30);
  TH1F* hMass = new TH1F("hMass", "hMass", 1000, 0, 10);
  TH1F* hMassBg = new TH1F("hMassBg", "hMassBg", 1000, 0, 10);

  for (int iM=0; iM<nMultBins; iM++) {
    char* name;
    name = Form("hMult_%i", iM);
    hMult[iM] = new TH1F(name, name, 1000, 0, 1000);

    name = Form("hSig_%i", iM);
    hSig[iM] = new TH2F(name, name, 60, -PI/2, 3*PI/2, 60, -dEtaMax, dEtaMax);
    hSig[iM]->Sumw2();
    name = Form("hBkg_%i", iM);
    hBkg[iM] = new TH2F(name, name, 60, -PI/2, 3*PI/2, 60, -dEtaMax, dEtaMax);
    hBkg[iM]->Sumw2();
    name = Form("hSB_%i", iM);
    hSB[iM] = new TH2F(name, name, 60, -PI/2, 3*PI/2, 60, -dEtaMax, dEtaMax);
    hSB[iM]->Sumw2();

    for (int iPt=0; iPt<nPtBins; iPt++) {
      name = Form("hSigPt_%i_%i", iM, iPt);
      hSigPt[iM][iPt] = new TH2F(name, name, 60, -PI/2, 3*PI/2, 60, -dEtaMax, dEtaMax);
      hSigPt[iM][iPt]->Sumw2();
      name = Form("hBkgPt_%i_%i", iM, iPt);
      hBkgPt[iM][iPt] = new TH2F(name, name, 60, -PI/2, 3*PI/2, 60, -dEtaMax, dEtaMax);
      hBkgPt[iM][iPt]->Sumw2();
      name = Form("hSBPt_%i_%i", iM, iPt);
      hSBPt[iM][iPt] = new TH2F(name, name, 60, -PI/2, 3*PI/2, 60, -dEtaMax, dEtaMax);
      hSBPt[iM][iPt]->Sumw2();
    }
    
  }
  
  TH2F* hDR = new TH2F("hDR", "hDR", 100, -0.01, 0.01, 100, -0.005, 0.005);
  hDR->SetTitle("#Delta#phi-#Delta#eta Distribution (no cut);#Delta#phi;#Delta#eta");
  TH2F* hDRcut = new TH2F("hDRcut", "hDRcut", 100, -0.2, 0.2, 100, -0.1, 0.1);
  hDRcut->SetTitle("#Delta#phi-#Delta#eta Distribution;#Delta#phi;#Delta#eta");
  hDEta->SetLineColor(kRed);
  hDPhi->SetLineColor(kBlue);
  TH2F* hMixEff = new TH2F("hMixEff", "", 100, 0., 500., 11, -0.05, 1.05);
  hMixEff->SetTitle("Mixing efficiency vs. multiplicity;"
		    "multiplicity;"
		    "accepted/sampled tracks");
  // ---------------------------------

  Int_t nEvents = tree->GetEntries();
  Int_t nAcceptedEvents = 0;

  TStopwatch* watch = new TStopwatch();
  watch->Start();
  for (int n=0; n<nEvents; n++) {
    evBranch->GetEntry(n);
    if (!ac->IsEventOk(*ev, minVc, maxNTracklets, zMin, zMax))
      continue;

    trBranch->GetEntry(n);
    if (n % 100 == 0) cout << n << endl;

    hMultAll->Fill(ev->fNSelTracks);     // No cuts

    Int_t multBin = hMultBins->FindBin(ev->fNSelTracks) - 1;
    Int_t zvtxBin = hZvtxBins->FindBin(ev->fVz) - 1;
    if (!ac->InBounds(multBin, 0, nMultBins-1)) 
      continue;
    if (!ac->InBounds(zvtxBin, 0, nZvtxBins-1)) 
      continue;
    if (runNumber == -1) {
      runNumber = firstRunNumber = ev->fRun;
    }
    runNumber = ev->fRun;
    if (runNumber != firstRunNumber) {
      cout << "Warning: Run number has changed!" << endl; 
      firstRunNumber = runNumber;
    }
    nAcceptedEvents++;
    
    Int_t ntracks = trk->GetEntries();
    hMultSel->Fill(ntracks);            // After event selection
    hMult[multBin]->Fill(ntracks);
    hVc->Fill(ev->fVc);
    hZvtx->Fill(ev->fVz);
    hNTracklets->Fill(ev->fNTracklets);
    EventPool* pool = ac->GetEventPool(multBin, zvtxBin);

    if (pool->IsPoolReady()) {
      for (int i=0; i<ntracks; i++) {
	MyPart* t1 = (MyPart*)trk->At(i);
	hPtAll->Fill(t1->Pt());
	if (!ac->IsTrackOk(*t1, etaMin, etaMax, ptMin, ptMax))
	  continue;
	
	hEta->Fill(t1->Eta());
	hPtSel->Fill(t1->Pt());
	Int_t ptBin1 = hPtBins->FindBin(t1->Pt()) - 1;

	// --- Same-event correlation loop ---
	for (int j=i+1; j<ntracks; j++) {
	  MyPart* t2 = (MyPart*)trk->At(j);
	  if (!ac->IsTrackOk(*t2, etaMin, etaMax, ptMin, ptMax))
	    continue;
	  
	  Double_t deta = ac->DeltaEta(*t1, *t2);
	  Double_t dphi = ac->DeltaPhi(*t1, *t2);
	  hDR->Fill(dphi, deta); // No cut - contamination pk.
	  
	  if (ac->IsPairOk(*t1, *t2)) {

	    hDEta->Fill(deta);
	    hDPhi->Fill(dphi);
	    hSig[multBin]->Fill(dphi, deta);
	    hDRcut->Fill(dphi, deta); // Show elliptical hole from cut 
	    hMass->Fill(ac->InvMass(*t1, *t2));

	    Int_t ptBin2 = hPtBins->FindBin(t2->Pt()) - 1;
	    if (ptBin1==ptBin2)
	      hSigPt[multBin][ptBin1]->Fill(dphi, deta);
	  }
	
	  Int_t jMix = 0, jWatcher = 0, jLimit = 20*nMix;
	  for (jMix=0; jMix<nMix;) {
	    if (++jWatcher > jLimit) {
	      cout << "Warning: mix loop count > " << jLimit 
		   << ". Mixed with " << jMix << " tracks." 
		   << endl;
	      break;
	    }
	    MyPart* tbg = pool->GetRandomTrack();
	    if (!ac->IsTrackOk(*tbg, etaMin, etaMax, ptMin, ptMax))
	      continue;
	    
	    if (ac->IsMixedPairOk(*t1, *tbg)) {
	      hMassBg->Fill(ac->InvMass(*t1, *tbg));
	      Double_t deta_bg = ac->DeltaEta(*t1, *tbg);
	      Double_t dphi_bg = ac->DeltaPhi(*t1, *tbg);
	      hBkg[multBin]->Fill(dphi_bg, deta_bg);

	      Int_t ptBin3 = hPtBins->FindBin(tbg->Pt()) - 1;
	      if (ptBin1==ptBin3)
		hBkgPt[multBin][ptBin3]->Fill(dphi_bg, deta_bg);
	      
	      jMix++;
	    }
	  }
	  hMixEff->Fill(ntracks, (Double_t)jMix/jWatcher);
	}
      }
    }
    else {
      ac->UpdatePools(n, ev, trk);
    } 
  }
  watch->Stop();
  cout << "\nFinished event loop.\n" << endl;
  cout << nEvents << " events sampled from " << dataFile << endl;
  cout << nAcceptedEvents << " events selected." << endl;
  cout << "Used ROOT " << gROOT->GetVersion() << endl;
  cout << " CPU time:  " << watch->CpuTime()  << endl;
  cout << " Real time: " << watch->RealTime() << endl;

  cuts->Add(new TNamed("run_number", Form("Run number = %i", runNumber )));
  cuts->Add(new TNamed("n_events_file",Form("Events in file = %lld", tree->GetEntries())));
  cuts->Add(new TNamed("n_events_loop", Form("Events looped = %i", nEvents )));
  cuts->Add(new TNamed("n_events_kept", Form("Events kept = %i", nAcceptedEvents )));
  cuts->Add(new TNamed("cpu_time", Form("CPU time = %f", watch->RealTime() )));
  cuts->Add(new TNamed("real_time", Form("Real time = %f", watch->RealTime() )));
  cuts->Add(new TNamed("root_version", Form("ROOT version = %s", gROOT->GetVersion() )));
  cuts->Add(new TNamed("hostname", Form("Hostname = %s", gSystem->HostName() )));

  hMassBg->Scale(1./nMix);
  for (int iM=0; iM<nMultBins; iM++) {
    hBkg[iM]->Scale(1./nMix);
    hSB[iM]->Divide(hSig[iM], hBkg[iM]);  
    
    for (int iPt=0; iPt<nPtBins; iPt++) {
      hBkgPt[iM][iPt]->Scale(1./nMix);
      hSBPt[iM][iPt]->Divide(hSigPt[iM][iPt], hBkgPt[iM][iPt]);  
    }
  }
  
  outFile->Write();
  cuts->Write("cuts", TObject::kSingleKey);
  outFile->Close();
  return;
}

