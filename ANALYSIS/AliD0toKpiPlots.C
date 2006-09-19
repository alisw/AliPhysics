void AliD0toKpiPlots(const Char_t *inName="AliD0toKpi.root",
		     const Char_t *outName="D0histograms.root") {
  //--------------------------------------------------------------------------
  // This macro histograms many variables of D0->Kpi candidates
  //
  //     Andrea Dainese, andrea.dainese@lnl.infn.it
  //--------------------------------------------------------------------------

  gSystem->Load("libANALYSIS.so");

  // set of cuts
  Double_t D0Cuts[9] = {0.1,         // mass [GeV]
			1000000.,          // dca [micron]
			1.1,           // cosThetaStar
			0.,           // pT K [GeV/c]
			0.,           // pT Pi [GeV/c]    
			100000.,       // d0K upper [micron]
			100000.,       // d0Pi upper [micron]
			10000000000.,            // d0d0 [micron^2]
			-1.1};          // cosThetaPointing

  // number of events (for normalization)
  Bool_t normalize = kFALSE;
  Double_t  events = 1.;



  // define histograms
  TH1F *hptK = new TH1F("hptK","\"K\" p_{t} distribution",50,0,10);
  hptK->SetXTitle("p_{t} [GeV]");

  TH1F *hptPi = new TH1F("hptPi","\"#pi\" p_{t} distribution",50,0,10);
  hptPi->SetXTitle("p_{t} [GeV]");

  TH1F *hDCA = new TH1F("hDCA","DCA",50,0,1000);
  hDCA->SetXTitle("dca [#mu m]");

  TH1F *hptD0 = new TH1F("hptD0","D^{0} p_{t} distribution",40,0,40);
  hptD0->SetXTitle("p_{t} [GeV]");

  TH1F *hyD0 = new TH1F("hyD0","D^{0} rapidity distribution",50,-2,2);
  hyD0->SetXTitle("y");

  TH1F *hCPtaD0 = new TH1F("hCPtaD0","cosine of pointing angle distribution",100,-1,1);
  hCPtaD0->SetXTitle("cos #theta_{point}");

  TH1F *hCPtaXY = new TH1F("hCPtaXY","cosine of pointing angle in (x,y) plane",100,-1,1);
  hCPtaXY->SetXTitle("cos #theta_{point}");

  TH1F *hCts = new TH1F("hCts","cosine of decay angle",50,-1.2,1.2);
  hCts->SetXTitle("cos #theta^{*}");

  TH2F *hCtsVsPtK = new TH2F("hCtsVsPtK","cosine of decay angle VS \"K\" p_{t}",50,0,5,50,-1,1);
  hCtsVsPtK->SetYTitle("cos #theta^{*}");
  hCtsVsPtK->SetXTitle("p_{t} [GeV]");

  TH1F *hd0d0 = new TH1F("hd0d0","Product of the impact parameters",100,-100000,100000);
  hd0d0->SetXTitle("d_{0}^{K} #times d_{0}^{#pi} [#mu m^{2}]");

  TH1F *hd0K = new TH1F("hd0K","Impact parameter of \"K\"",100,-5000,5000);
  hd0K->SetXTitle("d_{0}^{K} [#mu m]");

  TH1F *hd0Pi = new TH1F("hd0Pi","Impact parameter of \"#pi\"",100,-5000,5000);
  hd0Pi->SetXTitle("d_{0}^{#pi} [#mu m]");

  TH2F *hCPtaVsd0d0 = new TH2F("hCPtaVsd0d0","cos #theta_{point} vs d_{0}^{K} #times d_{0}^{#pi}",100,-100000,100000,100,-1,1);
  hCPtaVsd0d0->SetXTitle("d_{0}^{K} #times d_{0}^{#pi} [#mu m^{2}]");
  hCPtaVsd0d0->SetYTitle("cos #theta_{point}");

  TH2F *hCPtaVsd0d0zoom = new TH2F("hCPtaVsd0d0zoom","cos #theta_{point} vs d_{0}^{K} #times d_{0}^{#pi}",100,-100000,0,100,.9,1);
  hCPtaVsd0d0zoom->SetXTitle("d_{0}^{K} #times d_{0}^{#pi} [#mu m^{2}]");
  hCPtaVsd0d0zoom->SetYTitle("cos #theta_{point}");

  TH2F *hd0d0VSptD0 = new TH2F("hd0d0VSptD0","d_{0}^{K} #times d_{0}^{#pi} VS D^{0} p_{t}",50,0,25,100,-120000,120000);
  hd0d0VSptD0->SetYTitle("d_{0}^{K} #times d_{0}^{#pi} [#mu m^{2}]");
  hd0d0VSptD0->SetXTitle("D^{0} p_{t} [GeV]");

  TH1F *hMass = new TH1F("hMass","Invariant mass distribution",50,1.765,1.965);
  hMass->SetXTitle("M[K,#pi] [GeV]");

  TH2F *hArm = new TH2F("hArm","Armenteros plot",50,-2,2,50,0,1);
  hArm->SetXTitle("#alpha");
  hArm->SetYTitle("q_{t}");

  // open input file and get tree
  TFile *inFile = TFile::Open(inName);

  TTree *treeD0 = (TTree*)inFile->Get("TreeD0");
  AliD0toKpi *D = 0; 
  treeD0->SetBranchAddress("D0toKpi",&D);
  Int_t entries = (Int_t)treeD0->GetEntries();

  printf("+++\n+++ Number of D0 in tree:  %d\n+++\n",entries);

  Double_t MD0,MD0bar,ctsD0,ctsD0bar,ctsPiD0,ctsPiD0bar;
  Double_t WgtD0,WgtD0bar;
  Double_t sampleABC=0.;
  Int_t okD0=0,okD0bar=0;
  Int_t nSel = 0;
  Int_t ptbin;

  // loop on D0
  for(Int_t i=0; i<entries; i++) {
    if(i%10000==0) printf(" candidate %d of %d\n",i,entries);

    // get event from tree
    treeD0->GetEvent(i);
    //--- select the PID strategy & compute weights
    //    D->ApplyPID("TOFparam_PbPb");
    //    D->ComputeWgts();
    // get weights for the three samples A+B+C 
    //    D->GetWgts(WgtD0,WgtD0bar,"ABC");
    WgtD0 = 1.; WgtD0bar = 1.;

    // normalize to 1 event
    if(normalize) { WgtD0 /= events; WgtD0bar /= events; }

    // check if candidate passes selection (as D0 or D0bar)
    D->Select(D0Cuts,okD0,okD0bar);

    // set weights to 0 if the candidate doesn't pass selection
    if(!okD0)    WgtD0=0.; 
    if(!okD0bar) WgtD0bar=0.;
    if(okD0 || okD0bar) nSel++;

    // count selected candidates
    sampleABC += WgtD0 + WgtD0bar;

    // inv mass and cosThetaStar
    D->InvMass(MD0,MD0bar);
    D->CosThetaStar(ctsD0,ctsD0bar);
    
    // fill histograms
    hptK->Fill(D->PtChild(1),WgtD0);
    hptK->Fill(D->PtChild(0),WgtD0bar);
    hptPi->Fill(D->PtChild(0),WgtD0);
    hptPi->Fill(D->PtChild(1),WgtD0bar);
    hd0K->Fill(D->Getd0Child(1),WgtD0);
    hd0K->Fill(D->Getd0Child(0),WgtD0bar);
    hd0Pi->Fill(D->Getd0Child(0),WgtD0);
    hd0Pi->Fill(D->Getd0Child(1),WgtD0bar);   
    hMass->Fill(MD0,WgtD0);
    hMass->Fill(MD0bar,WgtD0bar);
    hCts->Fill(ctsD0,WgtD0);
    hCts->Fill(ctsD0bar,WgtD0bar);
    hCtsVsPtK->Fill(D->PtChild(1),ctsD0,WgtD0);
    hCtsVsPtK->Fill(D->PtChild(0),ctsD0bar,WgtD0bar);
    hDCA->Fill(D->GetDCA(),WgtD0+WgtD0bar);
    hptD0->Fill(D->Pt(),WgtD0+WgtD0bar);
    hyD0->Fill(D->Rapidity(),WgtD0+WgtD0bar);
    hd0d0->Fill(D->ProdImpParams(),WgtD0+WgtD0bar);
    hCPtaD0->Fill(D->CosPointing(),WgtD0+WgtD0bar);
    hCPtaXY->Fill(D->CosPointingXY(),WgtD0+WgtD0bar);
    hCPtaVsd0d0->Fill(D->ProdImpParams(),D->CosPointing(),WgtD0+WgtD0bar);
    hd0d0VSptD0->Fill(D->Pt(),D->ProdImpParams(),WgtD0+WgtD0bar);
    hCPtaVsd0d0zoom->Fill(D->ProdImpParams(),D->CosPointing(),WgtD0+WgtD0bar);
    hArm->Fill(D->Alpha(),D->Qt(),WgtD0+WgtD0bar);
    
   
  } // end loop on D0 candidates

  inFile->Close();

  printf("\n\n --- Total number of candidates passing selection: %d\n\n --- Sum of weights sample A+B+C: %f\n\n",nSel,sampleABC); 
  
  // draw histograms
  TCanvas *c1 = new TCanvas("c1","pt K & pi",0,0,700,700);
  c1->SetLogy(); 
  hptK->Draw();
  hptPi->Draw("same");

  TCanvas *c2 = new TCanvas("c2","pt D0",0,0,700,700);
  c2->SetLogy(); 
  hptD0->Draw();

  TCanvas *c3 = new TCanvas("c3","rapidity D0",0,0,700,700);
  hyD0->Draw(); 

  TCanvas *c4 = new TCanvas("c4","pointing angle",0,0,700,700);
  hCPtaD0->Draw();

  TCanvas *c5 = new TCanvas("c5","d0 x d0",0,0,700,700);
  c5->SetLogy();
  hd0d0->Draw();

  TCanvas *c6 = new TCanvas("c6","pointing angle VS d0d0",0,0,700,700);
  c6->SetLogz();
  hCPtaVsd0d0->Draw("box");

  TCanvas *c7 = new TCanvas("c7","mass",0,0,700,700);
  hMass->Draw();

  TCanvas *c8 = new TCanvas("c8","armenteros",0,0,700,700);
  hArm->Draw("box");

  TCanvas *c9 = new TCanvas("c9","decay angle",0,0,700,700); 
  hCts->Draw();

  TCanvas *c10 = new TCanvas("c10","dca",0,0,700,700);
  c10->SetLogy();
  hDCA->Draw();

  TCanvas *c11 = new TCanvas("c11","d0 K & pi",0,0,700,700);
  c11->SetLogy();
  hd0K->Draw();
  hd0Pi->Draw("same");

  // write all histograms to file
  TFile *outFile = new TFile(outName,"recreate");
  hMass->Write();
  hDCA->Write(); 
  hCts->Write();
  hCtsVsPtK->Write();
  hArm->Write();
  hCPtaVsd0d0->Write();
  hd0d0VSptD0->Write();
  hCPtaVsd0d0zoom->Write();
  hd0d0->Write();
  hCPtaD0->Write();
  hCPtaXY->Write();
  hptK->Write();
  hptPi->Write();
  hptD0->Write();
  hyD0->Write();
  hd0K->Write();
  hd0Pi->Write();
  outFile->Close();
  
  return;
}




