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

  TString name("cTau_0090"); //centrality
  //TString name("cTau_0005"); //centrality
  //TString name("cTau_2040"); //centrality
  //TString name("cTau_4060"); //centrality
  //TString name("cTau_6080"); //centrality
  //TString name("cTau_8090"); //centrality


  //+++++ Real data histograms
  TFile::Open("new/real/all/" + name + ".root");
  TList *list=(TList *)gFile->Get(name); 

  TH2F *fK0sSi=(TH2F*)list->FindObject("fK0sSi");fK0sSi->Sumw2();
  TH2F *fLambdaSi=(TH2F*)list->FindObject("fLambdaSi");fLambdaSi->Sumw2();
  TH1F *fXiSiP=(TH1F*)list->FindObject("fXiSiP"); fXiSiP->Sumw2();

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

  TH1F *fXiSiPMC=(TH1F*)listmc->FindObject("fXiSiP");fXiSiPMC->Sumw2();
  TH3F *fLambdaFromXi=(TH3F*)listmc->FindObject("fLambdaFromXi");
  fLambdaFromXi->Sumw2();
  //+++++


  Double_t myExp2(Double_t *v, Double_t *p);
  void Correct(TH1 *rw, const TH1 *as, const TH1 *mc);

  const Int_t    iMax=2;   // Number of V0 particles
  const TString  pnam[]={"K^{0}_{s}", "#Lambda"};
  const Double_t brch[]={0.69, 0.64};
  const Double_t mass[]={0.4977, 1.115};
  const Double_t ctau[]={2.68, 7.89};

  const TH2 *in[]={
    fK0sSi,fK0sAs,fK0sMC,
    fLambdaSi,fLambdaAs,fLambdaMC
  };
  Int_t    nbx=fK0sSi->GetNbinsX();
  Double_t wbx=fK0sSi->GetBinWidth(3);

  for (Int_t i=0; i<iMax; i++) {
      TH2 *cr=(TH2*)in[i*3 + 0]->Clone();
      Correct(cr, in[i*3 + 1], in[i*3 + 2]);

      //++++ pT spectrum
      TH1 *pt=cr->ProjectionX("_px",0,-1,"e");
      TString ptName = pnam[i] + " p_{T} spectrum";
      pt->SetTitle(ptName.Data());

      pt->Scale(1/brch[i]);    // branching
      pt->Scale(1/(2*0.5));    // rapidity
      pt->Scale(1/wbx);        // bin width in pt 
      pt->Scale(1/nEvents);    // number of events
       
      new TCanvas; 
      pt->Draw(); 

 
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
      cout<<pnam[i]<<"  \tstatus: "<<status<<"\tchi2/ndf: "<<chi2<<
           "\tc*tau: "<<slope<<"+/-"<<error<<endl;
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

