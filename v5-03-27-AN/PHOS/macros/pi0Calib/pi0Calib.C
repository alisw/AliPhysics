Double_t rpol2(Double_t *x, Double_t *par)
{
  Double_t rp2 = par[0] + par[1]*x[0] + par[2]*x[0]*x[0];
  return rp2;
}

TF1* InitPol2()
{
  TF1* f1 = new TF1("rpol2",rpol2,0.,300.,3); 
  f1->SetParNames("p0","p1","p2");

  f1->SetParameter("p0",1.);
  f1->SetParameter("p1",1.);
  f1->SetParameter("p2",1.);

  return f1;
}

Double_t gausPol2(Double_t *x, Double_t *par)
{
  //gaus+second order polynomial.
  
  Double_t gp2 = par[0] * TMath::Exp( -(x[0]-par[1])*(x[0]-par[1]) /
                                      (2*par[2]*par[2]) ) /
    + par[3] + par[4]*x[0] + par[5]*x[0]*x[0]; 
  
  return gp2;
}



void pi0Calib(const char* file="Sum_All_Emc.root")
{

  AliCDBManager::Instance()->SetDefaultStorage("local://$ALICE_ROOT");
  AliCDBManager::Instance()->SetSpecificStorage("PHOS/*","local://./");

  AliPHOSCalibData cdb(0);

  Int_t fMinEntr = 50; // min 50 pi0 canditates per cell
  Double_t xmin = 100; // fit range lower end
  Double_t xmax = 170; // fit range upper end

  const Float_t pi0m = 134.98; 
  Double_t pi0f;
  const Int_t fMod=2;
  char hname[128];
  TH1F* h1=0;
  Double_t cc_i, cc_ii;
  
  TH1F* hPi0Mass = new TH1F("pi0_mass","#pi^{0} mass in calibrated cells",100,0.,300.);
  
  TFile f(file);
  
  for(Int_t iMod=0; iMod<5; iMod++) {
    for(Int_t iX=0; iX<64; iX++) {
      for(Int_t iZ=0; iZ<56; iZ++) {

	sprintf(hname,"%d_%d_%d",iMod,iX,iZ);
	h1 = (TH1F*)f.Get(hname);
	if(!h1) continue;
	if(h1->GetEntries()<fMinEntr) continue;
	
	TF1* rp2 = InitPol2();
	rp2->SetRange(180,260);
// 	rp2->SetRange(50.,110.);
	h1->Fit(rp2,"RNQ");
	rp2->Print();
	
	h1->GetXaxis()->SetRange(17,210);
	Double_t max = h1->GetBinCenter(h1->GetMaximumBin()); // peak

	TF1* f1 = new TF1("gpol2",gausPol2,xmin,xmax,6);
	f1->SetParNames("Constant","Mean","Sigma","p0","p1","p2");

	printf("Trying to fit %s, %d entries, peak at %.3f\n",
	       h1->GetName(),h1->GetEntries(),max);

	f1->SetParameter("Constant",h1->GetMaximum());
	f1->SetParameter("Mean",max);
	f1->SetParameter("Sigma",1.);
	f1->SetParameter("p0",rp2->GetParameter("p0"));
	f1->SetParameter("p1",rp2->GetParameter("p1"));
	f1->SetParameter("p2",rp2->GetParameter("p2"));
 
	f1->SetLineColor(kBlue);
	h1->Fit(f1,"R+");
	
	pi0f = f1->GetParameter("Mean");
	hPi0Mass->Fill(pi0f);

	//CC from the previous iteration.
	cc_i = cdb.GetADCchannelEmc(iMod+1,iZ+1,iX+1);

	if(pi0f>100. && pi0f<200.) {
	  Double_t k = pi0m/pi0f;
	  cc_ii = cc_i*(1+k**2)/2.;
	  cdb.SetADCchannelEmc(iMod+1,iZ+1,iX+1,cc_ii);
	}
	
	printf("%s: %d entries, mean=%.3f, peak=%.3f, rms= %.3f. pi0 = %.3f,
             cc_i=%f, cc_ii=%f\n",
	       h1->GetName(),h1->GetEntries(),h1->GetMean(),max,h1->GetRMS(),
	       pi0f,cc_i,cc_ii);

      }
    }
  }
  
  TFile ff("pi0Calib.root","recreate");
  hPi0Mass->Write();
  
  //Writing new calibration coefficients (CCs) to OCDB
  
  AliCDBMetaData *md= new AliCDBMetaData();
  md->SetResponsible("Boris Polishchuk");
  
  cdb.WriteEmc(0,999999,md);
  cdb.WriteCpv(0,999999,md);
  
}
