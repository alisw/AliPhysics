void testGen()
{
gROOT.Reset();
//
// ************* Parameters for AliGenHIJING generator **************
//
AliGenZDC *gener = new AliGenZDC();
gener->SetParticle(kProton);
gener->SetMomentum(2760.);
gener->SetDir(0,0,0,1);
gener->SetFermi(1);
gener->SetDiv(0.000032,0.0001,2);
gener->Init();
//
// ************* Creating canvas, pads & histograms **************
//
TCanvas *c1 = new TCanvas("c1","Nucleon Momentum in LAB RS",0,10,580,700);
c1->SetFillColor(38);
pad11 = new TPad("pad11"," ",0.01,0.51,0.49,0.99);
pad12 = new TPad("pad12"," ",0.51,0.51,0.99,0.99);
pad13 = new TPad("pad13"," ",0.01,0.01,0.99,0.49);
pad11->SetFillColor(18);
pad12->SetFillColor(18);
pad13->SetFillColor(18);
pad11->Draw();
pad12->Draw();
pad13->Draw();

TCanvas *c2 = new TCanvas("c2","Nucleon Momentum with Fermi and Divergence",600,10,600,700);
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

TH1F *hpx  = new TH1F("hpx","Nucleon momentum Px",100,-1.,1.);
TH1F *hpy  = new TH1F("hpy","Nucleon momentum Py",100,-1.,1.);
TH1F *hpz  = new TH1F("hpz","Nucleon momentum Pz",100,2000.,3000.);

TH1F *hpbx = new TH1F("hpbx","Px with Fermi and Divergence",100,-1.,1.);
TH1F *hpby = new TH1F("hpby","Py with Fermi and Divergence",100,-1.,1.);
TH1F *hpbz = new TH1F("hpbz","Pz with Fermi and Divergence",100,0.,6000.);

// ************* Generation of events **************
//
for(Int_t i=0; i<=1000; i++){
   gener->Generate();
  //
   // ************* Getting momenta **************
   //
   Double_t px = gener->GetInMomentum(0);
   Double_t py = gener->GetInMomentum(1);
   Double_t pz = gener->GetInMomentum(2);
//   printf("Initial momentum -> px = %f, py = %f, pz = %f \n", px,py,pz);
   Double_t ptrackx = gener->GetTrackMomentum(0);
   Double_t ptracky = gener->GetTrackMomentum(1);
   Double_t ptrackz = gener->GetTrackMomentum(2);
//   printf("Track momentum -> px = %f, py = %f, pz = %f \n\n", ptrackx,ptracky,ptrackz);
   //
   // ************* Filling histograms **************
   //
   hpx->Fill(px);
   hpy->Fill(py);
   hpz->Fill(pz);
   hpbx->Fill(ptrackx);
   hpby->Fill(ptracky);
   hpbz->Fill(ptrackz);
   //
   // ************* Drawing histograms **************
   //
   pad11->cd();
//   pad11->GetFrame()->SetFillColor(10);
//   pad11->GetFrame()->SetBorderMode(-1);
//   pad11->GetFrame()->SetBorderSize(12);
   hpx->Draw();

   pad12->cd();
//   pad12->GetFrame()->SetFillColor(10);
//   pad12->GetFrame()->SetBorderMode(-1);
//   pad12->GetFrame()->SetBorderSize(12);
   hpy->Draw();
   
   pad13->cd();
//   pad13->GetFrame()->SetFillColor(10);
//   pad13->GetFrame()->SetBorderMode(-1);
//   pad13->GetFrame()->SetBorderSize(12);
   hpz->Draw();
   
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
