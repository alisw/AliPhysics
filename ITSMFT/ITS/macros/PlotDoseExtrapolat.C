Double_t OccupancyAtZAndR(Double_t z=0., Double_t r=2.2, Double_t occ0=2800, Double_t r0=3.7);
Double_t DoseAtZAndR(Double_t z=0., Double_t r=2.2, Double_t dose0=2800, Double_t z0=14.1, Double_t r0=3.7);
Double_t focc(Double_t* zr, Double_t* par);
Double_t fdosez(Double_t* zz, Double_t* par);
Double_t fdose(Double_t* zr, Double_t* par);

void PlotDoseExtrapolat(Bool_t zoom=kTRUE) {
//Double_t sigma=5.9; // PbPb 2010 run 
 Double_t sigma=7.94; // nominal parameter ( https://edms.cern.ch/file/445882/5/Vol_1_Chapter_21.pdf )
TFile *f = new TFile("prova.root","RECREATE");
//
//// Now the estimate of the occupancy
//
TF1 *occh = new TF1("Dose vs z at given r ",fdosez,-60.,60.,3);
//
Double_t r0=3.7;
Double_t z0=14.1;
Double_t dose0=2.2E+3 / 1.;
occh->SetParameter(0,1);  //  average occupancy at r0
occh->SetParameter(2,r0);  // r0
occh->SetParameter(1,z0);  // z0
Double_t tot=occh->Integral(-1.*z0,z0);
tot/=2*z0;
cout << tot << endl;
// now you can normalize
TF2 *occ; 
if (zoom) occ = new TF2("Dose vs z and r ",fdose,-20.,20.,1.5,5.,3); // zoomed
else occ = new TF2("Dose vs z and r ",fdose,-60.,60.,2.,50.,3); // not zoomed
TCanvas *c4=new TCanvas();
c4->SetRightMargin(0.15);
//occh->Draw();
//c3->SetLogy();
// 1 sigma
occ->SetParameter(0,dose0/tot);  //  occupancy at r0
occ->SetParameter(2,r0);  // r0
occ->SetParameter(1,z0);  // z0
occ->GetXaxis()->SetTitle("z (cm)");
occ->GetYaxis()->SetTitle("r (cm)");
occ->GetZaxis()->SetTitle("Dose (Gy)");
occ->GetZaxis()->SetTitleOffset(1.25);
occ->SetContour(20);
occ->Draw("colz");
occ->Write();

//
f->Close();
}


Double_t DoseAtZAndR(Double_t z, Double_t r, Double_t occ0, Double_t z0, Double_t r0) {
if(r<=0) {
  cout << "r should be positive " << endl;
 return 0;
}

Double_t occ=0;
cout << " occ0=" << occ0 << "  r0=" << r0 << "  r=" << r << " z=" << z << endl;
if (TMath::Abs(z)<0.0000001) occ=occ0;
else occ=occ0*r0/r * TMath::Sin(TMath::ATan(r/z));
cout << occ << endl;
return TMath::Abs(occ);
}


Double_t focc(Double_t* zr, Double_t* par) {
Double_t z=zr[0];
Double_t r=zr[1];
Double_t occ0=par[0];
Double_t r0=par[1];
return OccupancyAtZAndR(z,r,occ0,r0);
}

Double_t fdosez(Double_t* zz, Double_t* par) {
Double_t z=zz[0];
Double_t r=par[2];
Double_t pippo=par[0];
Double_t r0=par[2];
Double_t z0=par[1];
return DoseAtZAndR(z,r,pippo,z0,r0);
}


Double_t fdose(Double_t* zr, Double_t* par) {
Double_t z=zr[0];
Double_t r=zr[1];
Double_t dose0=par[0];
Double_t z0=par[1];
Double_t r0=par[2];
return DoseAtZAndR(z,r,dose0,z0,r0);
}

