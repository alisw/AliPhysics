computePriors(){
  char filename[100];
  sprintf(filename,"yourfile.root");
  TFile *f = new TFile(filename);
  TTree *t = f->Get("yourtree");

  Float_t mass[7] = {0.000511,0.139570,0.493677,0.938272,1.877837,2.817402,1.408054};

  Int_t nev = t->GetEntries();

  printf("ntracks = %i\n",nev);

  const Int_t maxstep = 10;

  TH1D *hpriors[7][maxstep];
  TH1D *hratio[7];

  for(Int_t i=0;i<7;i++){
    for(Int_t j=0;j<maxstep;j++){
      char histoname[100];
      sprintf(histoname,"priors%istep%i",i,j);
      hpriors[i][j] = new TH1D(histoname,histoname,10,0,10);
      if(j==0){
	sprintf(histoname,"ratio%i",i);
	hratio[i] = new TH1D(histoname,histoname,10,0,10);
	for(Int_t k=1;k <= 100;k++){
	  hpriors[i][j]->SetBinContent(k,1);
	}
      }
    }
  }

  for(Int_t isp=0;isp < 7;isp++) hratio[isp]->Sumw2();

  // loop variable definition 
  Float_t pt, eta, weight=1.0;
  Float_t ppi, pka, ppr, pel, pmu, pde, ptr, phe, ptot;

  for(Int_t j=1; j <maxstep;j++){
    for(Int_t isp=0;isp<7;isp++) hratio[isp]->Reset();
    printf("step %i\n",j);
    for(Int_t i=0;i<nev;i++){
      t->GetEvent(i);
      Int_t sw = t->GetLeaf("kTOF")->GetValue() * t->GetLeaf("kTPC")->GetValue(); 
      if(sw){
	pt = t->GetLeaf("pt")->GetValue();
	eta = t->GetLeaf("eta")->GetValue();
	//	weight = 1.0; // if you have fill the tree with weights
	if(pt > 0 && pt < 10){
	  ppi = t->GetLeaf("tofW")->GetValue(2)*t->GetLeaf("tpcW")->GetValue(2)*hpriors[1][j-1]->Interpolate(pt);
	  pka = t->GetLeaf("tofW")->GetValue(3)*t->GetLeaf("tpcW")->GetValue(3)*hpriors[2][j-1]->Interpolate(pt);
	  ppr = t->GetLeaf("tofW")->GetValue(4)*t->GetLeaf("tpcW")->GetValue(4)*hpriors[3][j-1]->Interpolate(pt);
	  pel = t->GetLeaf("tofW")->GetValue(0)*t->GetLeaf("tpcW")->GetValue(0)*hpriors[0][j-1]->Interpolate(pt);
	  pmu = t->GetLeaf("tofW")->GetValue(1)*t->GetLeaf("tpcW")->GetValue(1)*hpriors[0][j-1]->Interpolate(pt);
	  pde = t->GetLeaf("tofW")->GetValue(5)*t->GetLeaf("tpcW")->GetValue(5)*hpriors[4][j-1]->Interpolate(pt);
	  ptr = t->GetLeaf("tofW")->GetValue(6)*t->GetLeaf("tpcW")->GetValue(6)*hpriors[5][j-1]->Interpolate(pt);
	  phe = t->GetLeaf("tofW")->GetValue(7)*t->GetLeaf("tpcW")->GetValue(7)*hpriors[6][j-1]->Interpolate(pt);
	  ptot = ppi+pka+ppr+pel+pde+ptr+phe+pmu;

	  if(ptot > 0){
	    // fill ratio for |y|<0.5
	    if(TMath::Abs(eta2y(pt,mass[1],eta)) < 0.5) hratio[1]->Fill(pt,ppi/ptot*weight);
	    if(TMath::Abs(eta2y(pt,mass[2],eta)) < 0.5) hratio[2]->Fill(pt,pka/ptot*weight);
	    if(TMath::Abs(eta2y(pt,mass[3],eta)) < 0.5) hratio[3]->Fill(pt,ppr/ptot*weight);
	    if(TMath::Abs(eta2y(pt,mass[0],eta)) < 0.5) hratio[0]->Fill(pt,pel/ptot*weight);
	    if(TMath::Abs(eta2y(pt,mass[4],eta)) < 0.5) hratio[4]->Fill(pt,pde/ptot*weight);
	    if(TMath::Abs(eta2y(pt,mass[5],eta)) < 0.5) hratio[5]->Fill(pt,ptr/ptot*weight);
	    if(TMath::Abs(eta2y(pt,mass[6],eta)) < 0.5) hratio[6]->Fill(pt,phe/ptot*weight);

	    // fill priors assuming #eta cut during the fill of the tree
	    hpriors[1][j]->Fill(pt,ppi/ptot*weight);
	    hpriors[2][j]->Fill(pt,pka/ptot*weight);
	    hpriors[3][j]->Fill(pt,ppr/ptot*weight);
	    hpriors[0][j]->Fill(pt,pel/ptot*weight);
	    hpriors[4][j]->Fill(pt,pde/ptot*weight);
	    hpriors[5][j]->Fill(pt,ptr/ptot*weight);
	    hpriors[6][j]->Fill(pt,phe/ptot*weight);
	  }
	}
      }
    }

    for(Int_t isp=0;isp < 7;isp++) hpriors[isp][j]->Sumw2();

    // Normalize abundances to the pion one
    hpriors[0][j]->Divide(hpriors[1][j]);
    hpriors[2][j]->Divide(hpriors[1][j]);
    hpriors[3][j]->Divide(hpriors[1][j]);
    hpriors[4][j]->Divide(hpriors[1][j]);
    hpriors[5][j]->Divide(hpriors[1][j]);
    hpriors[6][j]->Divide(hpriors[1][j]);
    hpriors[1][j]->Divide(hpriors[1][j]);

    // make ratios w.r.t. pions (|y|<0.5)
    hratio[0]->Divide(hratio[1]);
    hratio[2]->Divide(hratio[1]);
    hratio[3]->Divide(hratio[1]);
    hratio[4]->Divide(hratio[1]);
    hratio[5]->Divide(hratio[1]);
    hratio[6]->Divide(hratio[1]);
    hratio[1]->Divide(hratio[1]);
  }

  // Write output
  sprintf(filename,"priors.root");
  TFile *fout = new TFile(filename,"RECREATE");
  for(Int_t j=maxstep-2; j <maxstep;j++){ // last 2 steps
    hpriors[0][j]->Write();
    hpriors[1][j]->Write();
    hpriors[2][j]->Write();
    hpriors[3][j]->Write();
    hpriors[4][j]->Write();
    hpriors[5][j]->Write();
    hpriors[6][j]->Write();
  }
  hratio[0]->Write();
  hratio[1]->Write();
  hratio[2]->Write();
  hratio[3]->Write();
  hratio[4]->Write();
  hratio[5]->Write();
  hratio[6]->Write();
  fout->Close();
}
Double_t y2eta(Double_t pt, Double_t m, Double_t y){
  Double_t mt=TMath::Sqrt(m*m+pt*pt);
  return TMath::ASinH(mt/pt*TMath::SinH(y));
}
Double_t eta2y(Double_t pt, Double_t m, Double_t eta){
  Double_t mt=TMath::Sqrt(m*m+pt*pt);
  return TMath::ASinH(pt/mt*TMath::SinH(eta));
}
