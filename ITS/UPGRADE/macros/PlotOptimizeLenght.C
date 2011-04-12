//
// Macro to estimates the acceptance and the occupancy
//
Double_t EtaMaxAtGivenZvtx(Double_t r=3.9,Double_t z=14.1, Double_t zvtx=0);
Double_t OptimizeLenght(Double_t r0=2.2, Double_t eta=1.0, Double_t r1= 43.0, Double_t zmax=5.9);
Double_t fnc(Double_t* r, Double_t* par);
Double_t feta(Double_t* zz, Double_t* par);
Double_t zPipeFnc(Double_t* r, Double_t* par);
Double_t ZOfConicalPipeAtR(Double_t r=2.2, Double_t sig=5.9,Double_t thetadeg=9.0);
Double_t OccupancyAtZAndR(Double_t z=0., Double_t r=2.2, Double_t occ0=2800, Double_t r0=3.7);
Double_t focc(Double_t* zr, Double_t* par);

void PlotOptimizeLenght() {
//Double_t sigma=5.9; // PbPb 2010 run 
 Double_t sigma=7.94; // nominal parameter ( https://edms.cern.ch/file/445882/5/Vol_1_Chapter_21.pdf )
TFile *f = new TFile("prova.root","RECREATE");
TF1 *func1 = new TF1("z",fnc,1.,45.,3); 
TCanvas *c1=new TCanvas();
c1->SetGridx();
c1->SetGridy();
c1->SetLogy();
// 1 sigma
func1->SetParameter(0,1.0);   // eta 
func1->SetParameter(1,43.0);  // zeta
func1->SetParameter(2,sigma);   // sigma 
func1->Draw();
func1->GetXaxis()->SetTitle("radius (cm)");
func1->GetYaxis()->SetTitle("#pm z (cm)");
func1->Write();
// 2 sigma
TF1* func2=new TF1("z2",fnc,1.,45.,3);
func2->SetParameter(0,1.0);   // eta 
func2->SetParameter(1,43.0);  // zeta
func2->SetParameter(2,2*sigma);   // sigma
func2->SetLineStyle(10);
func2->Draw("SAME");
func2->Write();
// 0 sigma
TF1* func0=new TF1("z0",fnc,1.,45.,3);
func0->SetParameter(0,1.0);   // eta 
func0->SetParameter(1,43.0);  // zeta
func0->SetParameter(2,0*sigma);   // sigma
func0->SetLineStyle(2);
func0->Draw("SAME");
func0->Write();
//
//// now at eta=1.3
// 
// 0 sigma
TF1* func0_b=new TF1(*func0);
func0_b->SetParameter(0,1.3);   // eta
func0_b->SetLineColor(2);
func0_b->Draw("SAME");
func0_b->Write();
// 1 sigma
TF1* func1_b=new TF1(*func1);
func1_b->SetParameter(0,1.3);   // eta
func1_b->SetLineColor(2);
func1_b->Draw("SAME");
func1_b->Write();
// 2 sigma
TF1* func2_b=new TF1(*func2);
func2_b->SetParameter(0,1.3);   // eta
func2_b->SetLineColor(2);
func2_b->Draw("SAME");
func2_b->Write();
// Legend
TLatex* eta = new TLatex(37.,13.,"#eta = 1.0");
eta->Draw();
TLatex* eta_b = new TLatex(37.,20.,"#eta = 1.3");
eta_b->SetTextColor(2);
eta_b->Draw();
TLine* l2= new TLine(5.,76., 13.,76. );
TLine* l1= new TLine(5.,57., 13.,57.);
TLine* l0= new TLine(5.,41., 13.,41.);
l2->SetLineStyle(10);
l0->SetLineStyle(2);
l0->SetLineWidth(2);
l1->SetLineWidth(2);
l2->SetLineWidth(2);
TLatex* t2 = new TLatex(14., 76.,"zvtx = 2 #sigma");
TLatex* t1 = new TLatex(14., 57.,"zvtx = 1 #sigma");
TLatex* t0 = new TLatex(14., 41.,"zvtx = 0");
t2->Draw(); 
t1->Draw(); 
t0->Draw(); 
l0->Draw();
l1->Draw();
l2->Draw();
//
//// Now plot eta versus Z at fixed radius
//
TCanvas *c2=new TCanvas();
c1->SetGridx();
c1->SetGridy();
// 2.2 cm 
//
// 0 sigma
TF1* eta0=new TF1("eta0",feta,2.,30.,2);
eta0->SetParameter(0,2.2);  // radius
eta0->SetParameter(1,0.);   // sigma
eta0->SetLineStyle(2);
eta0->GetXaxis()->SetTitle("#pm z (cm)");
eta0->GetYaxis()->SetTitle("#eta_{max}");
eta0->SetMinimum(0.);
eta0->Draw();
eta0->Write();
// 1 sigma 
TF1 *eta1 = new TF1("eta1",feta,sigma,30.,2);
eta1->SetParameter(0,2.2);  // radius
eta1->SetParameter(1,sigma);   // sigma 
eta1->Write();
eta1->Draw("Same");
// 2 sigma
TF1* eta2=new TF1("eta2",feta,2*sigma,30.,2);
eta2->SetParameter(0,2.2);  // radius
eta2->SetParameter(1,2*sigma);   // sigma
eta2->SetLineStyle(10);
eta2->Draw("SAME");
eta2->Write();
// 2.5 cm 
//
TF1* eta0_b=new TF1(*eta0);
eta0_b->SetParameter(0,2.5);   // radius
eta0_b->SetLineColor(2);
eta0_b->Draw("SAME");
eta0_b->Write();
// 1 sigma
TF1* eta1_b=new TF1(*eta1);
eta1_b->SetParameter(0,2.5);   // radius
eta1_b->SetLineColor(2);
eta1_b->Draw("SAME");
eta1_b->Write();
// 2 sigma
TF1* eta2_b=new TF1(*eta2);
eta2_b->SetParameter(0,2.5);   // radius
eta2_b->SetLineColor(2);
eta2_b->Draw("SAME");
eta2_b->Write();
// 2.8 cm 
TF1* eta0_c=new TF1(*eta0);
eta0_c->SetParameter(0,2.8);   // radius
eta0_c->SetLineColor(4);
eta0_c->Draw("SAME");
eta0_c->Write();
// 1 sigma
TF1* eta1_c=new TF1(*eta1);
eta1_c->SetParameter(0,2.8);   // radius
eta1_c->SetLineColor(4);
eta1_c->Draw("SAME");
eta1_c->Write();
// 2 sigma
TF1* eta2_c=new TF1(*eta2);
eta2_c->SetParameter(0,2.8);   // radius
eta2_c->SetLineColor(4);
eta2_c->Draw("SAME");
eta2_c->Write();
//
//
// Legend
TLatex* rra = new TLatex(21.,1.3,"r = 2.2 cm");
rra->Draw();
TLatex* rrb = new TLatex(21.,1.0,"r = 2.5 cm");
rrb->SetTextColor(2);
rrb->Draw();
TLatex* rrc = new TLatex(21.,0.7,"r = 2.8 cm");
rrc->SetTextColor(4);
rrc->Draw();
TLine* ll2= new TLine(4.,3.2, 9.,3.2);
TLine* ll1= new TLine(4.,2.9, 9.,2.9);
TLine* ll0= new TLine(4.,2.6, 9.,2.6);
ll2->SetLineStyle(10);
ll0->SetLineStyle(2);
ll0->SetLineWidth(2);
ll1->SetLineWidth(2);
ll2->SetLineWidth(2);
TLatex* tt2 = new TLatex(10., 3.2,"zvtx = 2 #sigma");
TLatex* tt1 = new TLatex(10., 2.9,"zvtx = 1 #sigma");
TLatex* tt0 = new TLatex(10., 2.6,"zvtx = 0");
tt2->Draw();
tt1->Draw();
tt0->Draw();
ll0->Draw();
ll1->Draw();
ll2->Draw();
//
//// Now plot Boundary from Conical Beam Pipe
//
TF1 *bound1 = new TF1("z_r bound",zPipeFnc,2.,10.,2);
TCanvas *c3=new TCanvas();
c3->SetGridx();
c3->SetGridy();
//c3->SetLogy();
// 1 sigma
bound1->SetParameter(0,1.0*sigma);  // 1 sigma 
bound1->SetParameter(1,9.0);  // theta in deg
bound1->Draw();
bound1->GetXaxis()->SetTitle("radius (cm)");
bound1->GetYaxis()->SetTitle("z (cm)");
bound1->Write();
//
// 2 sigma
TF1* bound2=new TF1("z_r bound 2",zPipeFnc,2.,10.,2);
bound2->SetParameter(0,2.0*sigma);   // eta 
bound2->SetParameter(1,9.0);  // theta in deg
bound2->SetLineStyle(10);
bound2->Draw("SAME");
bound2->Write();
// 0 sigma
TF1* bound0=new TF1("z_r bound 0",zPipeFnc,2.,10.,2);
bound0->SetParameter(0,0.);   // eta 
bound0->SetParameter(1,9.0);  // theta in deg
bound0->SetLineStyle(2);
bound0->Draw("SAME");
bound0->Write();
// caption
//
TLine* hl2= new TLine(2.5,52., 4.,52. );
TLine* hl1= new TLine(2.5,46., 4.,46.);
TLine* hl0= new TLine(2.5,40., 4.,40.);
hl2->SetLineStyle(10);
hl0->SetLineStyle(2);
hl0->SetLineWidth(2);
hl1->SetLineWidth(2);
hl2->SetLineWidth(2);
TLatex* ht2 = new TLatex(4.2, 52.,"zvtx = - 2 #sigma");
TLatex* ht1 = new TLatex(4.2, 46.,"zvtx = - 1 #sigma");
TLatex* ht0 = new TLatex(4.2, 40.,"zvtx = 0");
ht2->Draw();
ht1->Draw();
ht0->Draw();
hl0->Draw();
hl1->Draw();
hl2->Draw();
//
//// Now the estimate of the occupancy
//
TF2 *occ = new TF2("Occupancy vs z and r ",focc,0.,10.,2.,5.,2);
TCanvas *c4=new TCanvas();
c4->SetRightMargin(0.15);
//c3->SetLogy();
// 1 sigma
occ->SetParameter(0,28);  //  occupancy at r0
occ->SetParameter(1,3.7);  // r0
occ->GetXaxis()->SetTitle("z (cm)");
occ->GetYaxis()->SetTitle("r (cm)");
occ->GetZaxis()->SetTitle("charged particle / cm^2");
occ->GetZaxis()->SetTitleOffset(1.25);
occ->SetContour(20);
occ->Draw("colz");
occ->Write();

//
f->Close();
}

Double_t feta(Double_t* zz, Double_t* par) {
Double_t z=zz[0];
Double_t r=par[0];
Double_t zvtx=par[1];
return EtaMaxAtGivenZvtx(r,z,zvtx);
}


Double_t fnc(Double_t* r, Double_t* par) {
Double_t r0=r[0];
Double_t eta=par[0];
Double_t r1=par[1];
Double_t zmax=par[2];
return OptimizeLenght(r0, eta, r1, zmax);
}

Double_t OptimizeLenght(Double_t r0, Double_t eta, Double_t r1, Double_t zmax) {
Double_t z0=0.;
Double_t theta= 2.0*TMath::ATan(TMath::Exp(-1*eta));
Double_t z1=r1/TMath::Tan(theta);
Double_t theta_max=TMath::ATan((z1-zmax)/r1);
  //printf("\n \n Diamond size x-y: %10.3e z: %10.3e\n \n", sigmaxy, sigmaz)
z0=zmax+r0*TMath::Tan(theta_max);
return z0;
}


Double_t EtaMaxAtGivenZvtx(Double_t r, Double_t z, Double_t zvtx) {
Double_t etamax=0;
if (r<=0) return 0.;
if (zvtx>=z) {
 cout << " z should be larger than zvtx = " << zvtx << endl;
 return 0.;
}
Double_t zeff= z - zvtx;
Double_t thetaHalf=r/zeff;
thetaHalf=0.5*TMath::ATan(thetaHalf);
etamax=-1.*TMath::Log(TMath::Tan(thetaHalf));
return etamax;
}

Double_t ZOfConicalPipeAtR(Double_t r, Double_t sig, Double_t thetadeg){
Double_t rpipe=2.0;
if (r<rpipe) {
 cout << " radius r less than beam pipe" << endl;
 return 0;
}
 Double_t a=TMath::Tan(thetadeg*TMath::Pi()/180.);
 Double_t b=sig*a;
 Double_t z=(r-b)/a;
 return z;
}

Double_t OccupancyAtZAndR(Double_t z, Double_t r, Double_t occ0, Double_t r0) {
if(r<=0) {
  cout << "r should be positive " << endl;
 return 0;
}
Double_t occ=0;
occ=occ0*r0/r * TMath::Sin(TMath::ATan(r/z));
return occ;
}

Double_t zPipeFnc(Double_t* r, Double_t* par) {
Double_t rr=r[0];
Double_t sig=par[0];
Double_t thetadeg=par[1];
return ZOfConicalPipeAtR(rr, sig, thetadeg);
}

Double_t focc(Double_t* zr, Double_t* par) {
Double_t z=zr[0];
Double_t r=zr[1];
Double_t occ0=par[0];
Double_t r0=par[1];
return OccupancyAtZAndR(z,r,occ0,r0);
}
