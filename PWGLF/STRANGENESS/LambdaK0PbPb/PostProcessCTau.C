#if !defined(__CINT__) || defined(__MAKECINT__)
   #include <TMath.h>
   #include <TROOT.h>
   #include <Riostream.h>
   #include <TCanvas.h>

   #include <TString.h>

   #include <TFile.h>
   #include <TList.h>
   #include <TH1F.h>
   #include <TH1D.h>
   #include <TF2.h>
   #include <TFitResult.h>
   #include <TFitResultPtr.h>
   #include <TH2F.h>
   #include <TH3F.h>
#endif
 
extern TROOT *gROOT;
//extern TFile *gFile;


void PostProcessCTau() {
  //-----------------------------------------------------------------
  // PostProcess the histograms produced by AliAnalysisTaskCTau* 
  // Produce:
  // -- V0 particle spectra (not corrected for the feeddown)
  // -- Estimated (anti)Lambda feed-down spectra (from Xi's only)
  // -- V0 particle c*tau distributions (feeddown corrected) 
  //-----------------------------------------------------------------

  //TString name("cTau_0090"); //centrality
  //TString name("cTau_0005"); //centrality
  //TString name("cTau_2040"); //centrality
  TString name("cTau_4060"); //centrality
  //TString name("cTau_6080"); //centrality
  //TString name("cTau_8090"); //centrality


  //+++++ Real data histograms
  TFile::Open("new/real/all/" + name + ".root");
  TList *list=(TList *)gFile->Get(name); 

  TH2F *fK0sSi=(TH2F*)list->FindObject("fK0sSi");fK0sSi->Sumw2();
  TH2F *fLambdaSi=(TH2F*)list->FindObject("fLambdaSi");fLambdaSi->Sumw2();
  TH1F *fXiSiP=(TH1F*)list->FindObject("fXiSiP"); fXiSiP->Sumw2();
  TH2F *fLambdaBarSi=
        (TH2F*)list->FindObject("fLambdaBarSi"); fLambdaBarSi->Sumw2();
  TH1F *fXiBarSiP=(TH1F*)list->FindObject("fXiBarSiP"); fXiBarSiP->Sumw2();

  TH1F *fMult=(TH1F*)list->FindObject("fMult");
  //+++++

  const Double_t nEvents=fMult->Integral();

  //+++++ MC histograms
  TFile::Open("new/mc_b_plus/all/" + name + "_mc.root");
  TList *listmc=(TList *)gFile->Get(name+"_mc"); 

  TH2F *fK0sMC=(TH2F*)listmc->FindObject("fK0sMC"); fK0sMC->Sumw2();
  TH2F *fK0sAs=(TH2F*)listmc->FindObject("fK0sAs"); fK0sAs->Sumw2();

  TH2F *fLambdaMC=(TH2F*)listmc->FindObject("fLambdaMC"); fLambdaMC->Sumw2();
  TH2F *fLambdaAs=(TH2F*)listmc->FindObject("fLambdaAs"); fLambdaAs->Sumw2();

  TH2F *fLambdaBarMC=
        (TH2F*)listmc->FindObject("fLambdaBarMC"); fLambdaBarMC->Sumw2();
  TH2F *fLambdaBarAs=
        (TH2F*)listmc->FindObject("fLambdaBarAs"); fLambdaBarAs->Sumw2();

  TH1F *fXiSiPMC=(TH1F*)listmc->FindObject("fXiSiP");fXiSiPMC->Sumw2();
  TH3F *fLambdaFromXi=(TH3F*)listmc->FindObject("fLambdaFromXi");
  fLambdaFromXi->Sumw2();

  TH1F *fXiBarSiPMC=(TH1F*)listmc->FindObject("fXiBarSiP");fXiBarSiPMC->Sumw2();
  TH3F *fLambdaBarFromXiBar=(TH3F*)listmc->FindObject("fLambdaBarFromXiBar");
  fLambdaBarFromXiBar->Sumw2();
  //+++++


  Double_t myExp2(Double_t *v, Double_t *p);
  void Correct(TH1 *rw, const TH1 *as, const TH1 *mc);
  void Normalise(Double_t, Double_t, Double_t, Double_t, TH1 *h);

  const Int_t    iMax=3;   // Number of V0 particles
  const TString  pnam[]={"K^{0}_{S}", "#Lambda", "#bar{#Lambda}"};
  const Double_t brch[]={0.69, 0.64, 0.64};
  const Double_t mass[]={0.4977, 1.115, 1.115};
  const Double_t ctau[]={2.68, 7.89, 7.89};

  const TH2 *in[]={
    fK0sSi,fK0sAs,fK0sMC,
    fLambdaSi,fLambdaAs,fLambdaMC,
    fLambdaBarSi,fLambdaBarAs,fLambdaBarMC
  };
  TH1 *fd[]={
    0,0,0,
    fLambdaFromXi,fXiSiP,fXiSiPMC,
    fLambdaBarFromXiBar,fXiBarSiP,fXiBarSiPMC
  };
  Double_t wbx=fK0sSi->GetBinWidth(3);
  Int_t nbx=fLambdaFromXi->GetNbinsX();
  Int_t nby=fLambdaFromXi->GetNbinsY();
  Int_t nbz=fLambdaFromXi->GetNbinsZ();
 
  for (Int_t i=0; i<iMax; i++) {
      TH2 *cr=(TH2*)in[i*3 + 0]->Clone();
      Correct(cr, in[i*3 + 1], in[i*3 + 2]);

      //++++ pT spectrum
      TH1 *pt=cr->ProjectionX("_px",0,-1,"e");
      TString ptName = pnam[i] + " p_{T} spectrum";
      pt->SetTitle(ptName.Data());
      Normalise(brch[i], 2*0.5, wbx, nEvents, pt);      
      Double_t eipt=0., ipt=pt->IntegralAndError(1, nbx, eipt);

      new TCanvas; 
      pt->Draw(); 

      if (i>0) {
      //++++ feeddown
	 TH3 *fd3=(TH3*)fd[3*i + 0];
	 TH1 *rl =(TH1*)fd[3*i + 1];
	 TH1 *mc =(TH1*)fd[3*i + 2];
         rl->Divide(mc);
         
         for (Int_t k=1; k<=nbx; k++) {
             for (Int_t l=1; l<=nby; l++) {
                 for (Int_t m=1; m<=nbz; m++) {
                     Float_t c=fd3->GetBinContent(k,l,m);
                     c *= rl->GetBinContent(m);
                     fd3->SetBinContent(k,l,m,c);
		 }
	     }
	 }

         TH2 *fd2=(TH2*)fd3->Project3D("yxe");
         TH1 *fd1=fd2->ProjectionX("_px",0,-1,"e");
         Normalise(brch[i], 2*0.5, wbx, nEvents, fd1);
         Double_t eifd=0., ifd=fd1->IntegralAndError(1, nbx, eifd);

         Double_t f=ifd/ipt;
         Double_t ef=1/ipt*TMath::Sqrt(eifd*eifd + f*f*eipt*eipt);
         cout<<endl<<"Global FD correction: "<<f<<"+/-"<<ef<<endl;

         //new TCanvas();
         fd1->Draw("same");

         cr->Add(fd2,-1);
      } 
      //continue;
 
      //++++ c*tau
      TF2 *f2=new TF2("myexpo2",myExp2,0.,10.,0.,100.,1+1+1+nbx);
      f2->SetParameter(0, ctau[i]);
      f2->FixParameter(1, ctau[i]);
      f2->FixParameter(2, mass[i]);
      for (Int_t p=1+1+1; p<=nbx+1+1; p++)  
          f2->SetParameter(p,fK0sSi->GetBinContent(p,1));

      new TCanvas; 
      TFitResultPtr r=cr->Fit(f2,"SQ");

      Int_t status   = r;
      Double_t chi2  = r->Chi2()/r->Ndf();  
      Double_t slope = r->Parameter(0);  
      Double_t error = r->ParError(0);  
      cout<<endl<<pnam[i]<<"  \tstatus: "<<status<<"   chi2/ndf: "<<chi2<<
	"\tc*tau: "<<slope<<"+/-"<<error<<endl<<endl;
  }

  return;
} 

Double_t myExp2(Double_t *v, Double_t *p) {
   Double_t pt=v[0];
   Double_t lt=v[1];

   Double_t ctau=p[1];
   Double_t mass=p[2];
   Double_t ct=mass*lt/pt;
   if ((lt < 2) || (ct > 2.5*ctau)) {
      TF1::RejectPoint();
      return 0;        
   }

   Int_t i=pt/0.1 + 1 + 1;
   return p[i]*TMath::Exp(-ct/p[0]);
}

void Correct(TH1 *rw, const TH1 *as, const TH1 *mc) {
  TH1 *eff = (TH1*)as->Clone();
  eff->Divide(as,mc,1,1,"b");
  rw->Divide(eff);
  delete eff;
} 

void Normalise(Double_t br, Double_t yw, Double_t bw, Double_t ne, TH1 *pt) {
   pt->Scale(1/br);    // branching ratio
   pt->Scale(1/yw);    // rapidity window
   pt->Scale(1/bw);    // bin width 
   pt->Scale(1/ne);    // number of events
}
