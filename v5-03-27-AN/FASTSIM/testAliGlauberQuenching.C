void testAliGlauberQuenching() {

  AliFastGlauber g;
  g.Init(2);
  g.SetCentralityClass(0.00,0.10); // centrality (fraction of geometrical cross section)

  AliQuenchingWeights afq;
  afq.InitMult();
  afq.SetQTransport(1.);
  afq.SetECMethod(0);
  afq.SetLengthMax(8.);
  afq.SampleEnergyLoss();

  TCanvas *c = new TCanvas("cELD","Energy Loss Distribution",0,0,800,500);
  c->Divide(2,1);

  for(Int_t itype=1;itype<=2;itype++){
    c->cd(itype);
    gPad->SetLogy();
    Char_t name[100];
    Char_t hname[100];
    if(itype==1){
      sprintf(name,"Energy Loss Distribution - Quarks;E_{loss} (GeV);#"); 
      sprintf(hname,"hQuarks"); 
    } else {
      sprintf(name,"Energy Loss Distribution - Gluons;E_{loss} (GeV);#"); 
      sprintf(hname,"hGluons"); 
    }

    TH1F *h = new TH1F(hname,name,100,0,200);

    for(Int_t i=0;i<10000;i++){
      if(i % 100 == 0) cout << "." << flush;
      Double_t ell;
      g.GetLength(ell);
      Double_t loss=afq.GetELossRandom(itype,ell);
      h->Fill(loss);
    }

    h->SetStats(kTRUE);
    if(itype==1){
      h->SetLineColor(1);
      h->DrawCopy();
    }
    else {
      h->SetLineColor(2);
      h->DrawCopy();
    }
    delete h;
  }

  c->Update();
}

// second example using ell distribution
void testFastAliGlauberQuenching(Char_t *fname) {

  TFile f(fname);
  TH1F *hEll=(TH1F*)f.Get("hEll");
  if(!hEll) return;

  AliQuenchingWeights afq;
  afq.InitMult();
  afq.SetQTransport(1.);
  afq.SetECMethod(0);
  afq.SetLengthMax(5.);
  afq.SampleEnergyLoss();

  TCanvas *c = new TCanvas("cELD","Energy Loss Distribution",0,0,800,500);
  c->Divide(2,1);

  for(Int_t itype=1;itype<=2;itype++){
    c->cd(itype);
    gPad->SetLogy();
    Char_t name[200];
    if(itype==1)
      sprintf(name,"Energy Loss Distribution - Quarks;E_{loss} (GeV);#"); 
    else 
      sprintf(name,"Energy Loss Distribution - Gluons;E_{loss} (GeV);#"); 
    Char_t hname[200];
    if(itype==1)
      sprintf(hname,"hQuarks"); 
    else 
      sprintf(hname,"hGluons"); 

    TH1F *h = afq.ComputeELossHisto(itype,1.,hEll);

    h->SetStats(kTRUE);
    if(itype==1){
      h->SetLineColor(1);
      h->DrawCopy();
    }
    else {
      h->SetLineColor(2);
      h->DrawCopy();
    }
    delete h;
  }

  c->Update();
}
