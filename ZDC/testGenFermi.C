void testGenFermi()
{
gROOT.Reset();
//
// ************* Parameters for AliGenHIJING generator **************
//
AliGenZDC *gener = new AliGenZDC();
gener->SetParticle(kNeutron);
gener->SetMomentum(2760.);
gener->SetDir(0,0,0,1);
gener->SetFermi(1);
gener->SetDiv(0,0,0);
gener->Init();
//
// ************* Creating canvas, pads & histograms **************
//
//TCanvas *c1 = new TCanvas("c1","Nucleon Momentum in LAB RS",0,10,580,700);
//pad11 = new TPad("pad11"," ",0.01,0.51,0.49,0.99);
//pad12 = new TPad("pad12"," ",0.51,0.51,0.99,0.99);
//pad13 = new TPad("pad13"," ",0.01,0.01,0.99,0.49);
//pad11->SetFillColor(18);
//pad12->SetFillColor(18);
//pad13->SetFillColor(18);
//pad11->Draw();
//pad12->Draw();
//pad13->Draw();

TCanvas *c2 = new TCanvas("c2","Nucleon Momentum Boosted with Fermi Momentum",600,10,600,700);
c2->SetFillColor(38);
pad21 = new TPad("pad21"," ",0.01,0.51,0.49,0.99);
pad22 = new TPad("pad22"," ",0.51,0.51,0.99,0.99);
pad23 = new TPad("pad23"," ",0.01,0.01,0.99,0.49);
pad21->SetFillColor(10);
pad22->SetFillColor(10);
pad23->SetFillColor(10);
pad21->Draw();
pad22->Draw();
pad23->Draw();

TCanvas *c3 = new TCanvas("c3","Fermi2Gaussian distributions",0,10,580,700);
c3->SetFillColor(38);
pad31 = new TPad("pad31"," ",0.01,0.51,0.99,0.99);
pad32 = new TPad("pad32"," ",0.01,0.01,0.99,0.49);
pad31->SetFillColor(10);
pad32->SetFillColor(10);
pad31->Draw();
pad32->Draw();

//TH1F *hpx  = new TH1F("hpx","Nucleon momentum Px",100,-100.,100.);
//TH1F *hpy  = new TH1F("hpy","Nucleon momentum Py",100,-100.,100.);
//TH1F *hpz  = new TH1F("hpz","Nucleon momentum Pz",100,-3000.,0.);

TH1F *hpbx = new TH1F("hpbx","Px boosted with Fermi momentum",100,-1.,1.);
TH1F *hpby = new TH1F("hpby","Py boosted with Fermi momentum",100,-1.,1.);
TH1F *hpbz = new TH1F("hpbz","Pz boosted with Fermi momentum",100,-8000.,0.);

TH1D *hdgp  = new TH1D("hdgp","Fermi Two Gaussian Distribution -> p",200,0.,200.);
TH1D *hdgn  = new TH1D("hdgn","Fermi Two Gaussian Distribution -> n",200,0.,200.);
//
// ************* Fermi Two Gaussian distributions **************
//
Double_t FermiDp[201], FermiDn[201];
for(Int_t i=0; i<=200; i++){
   FermiDp[i] = gener->GetFermi2p(i);
   FermiDn[i] = gener->GetFermi2n(i);
//     printf("	testGenZDC -> Fermi2p[%d] = %f, Fermi2n[%d] = %f\n",i,FermiDp[i],i,FermiDn[i]);
   hdgp->Fill((Axis_t)i,(Stat_t)FermiDp[i]);
   hdgn->Fill((Axis_t)i,(Stat_t)FermiDn[i]);
}
   
pad31->cd();
//pad31->GetFrame()->SetFillColor(10);
//pad31->GetFrame()->SetBorderMode(-1);
//pad31->GetFrame()->SetBorderSize(12);
hdgp->Draw();

pad32->cd();
//pad32->GetFrame()->SetFillColor(10);
//pad32->GetFrame()->SetBorderMode(-1);
//pad32->GetFrame()->SetBorderSize(12);
hdgn->Draw();
//
// ************* Generation of events **************
//
for(Int_t i=0; i<=10000; i++){
   gener->Generate();
  //
   // ************* Getting momenta **************
   //
//   Double_t px = gener->GetMomentum(0);
//   Double_t py = gener->GetMomentum(1);
//   Double_t pz = gener->GetMomentum(2);
   Double_t pboostx = gener->GetBoostMomentum(0);
   Double_t pboosty = gener->GetBoostMomentum(1);
   Double_t pboostz = gener->GetBoostMomentum(2);
//   printf ("	testGenZDC -> pz = %f pz_boost = %f\n",pz,pboostz);
   //
   // ************* Filling histograms **************
   //
//   hpx->Fill(px);
//   hpy->Fill(py);
//   hpz->Fill(pz);
   hpbx->Fill(pboostx);
   hpby->Fill(pboosty);
   hpbz->Fill(pboostz);
   //
   // ************* Drawing histograms **************
   //
//   pad11->cd();
//   pad11->GetFrame()->SetFillColor(10);
//   pad11->GetFrame()->SetBorderMode(-1);
//   pad11->GetFrame()->SetBorderSize(12);
//   hpx->Draw();

//   pad12->cd();
//   pad12->GetFrame()->SetFillColor(10);
//   pad12->GetFrame()->SetBorderMode(-1);
//   pad12->GetFrame()->SetBorderSize(12);
//   hpy->Draw();
   
//   pad13->cd();
//   pad13->GetFrame()->SetFillColor(10);
//   pad13->GetFrame()->SetBorderMode(-1);
//   pad13->GetFrame()->SetBorderSize(12);
//   hpz->Draw();
   
   pad21->cd();
//   pad21->GetFrame()->SetFillColor(10);
//   pad21->GetFrame()->SetBorderMode(-1);
//   pad21->GetFrame()->SetBorderSize(12);
   hpbx->Draw();
   
   pad22->cd();
//   pad22->GetFrame()->SetFillColor(10);
//   pad22->GetFrame()->SetBorderMode(-1);
//   pad22->GetFrame()->SetBorderSize(12);
   hpby->Draw();
   
   pad23->cd();
//   pad23->GetFrame()->SetFillColor(10);
//   pad23->GetFrame()->SetBorderMode(-1);
//   pad23->GetFrame()->SetBorderSize(12);
   hpbz->Draw();
}  
}
