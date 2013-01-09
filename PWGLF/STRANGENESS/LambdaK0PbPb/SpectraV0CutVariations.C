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

TFile *fAss;
TFile *fGen;
TFile *fRaw;

void Generated();
void Corrections(Float_t cMin, Float_t cMax, TString centr);
void RawYields(Float_t cMin, Float_t cMax, TString centr);
void Spectra(TString centr);

void SpectraV0CutVariations(Float_t cMin, Float_t cMax, TString centr) {

  fRaw=TFile::Open((centr+"/AliV0CutVariations.root").Data());
  fAss=TFile::Open((centr+"/AliV0CutVariationsMC.root").Data());
  fGen=TFile::Open("Generated.root");
  if (!fGen) {
     Generated();
     fGen=TFile::Open("Generated.root");
  }

  Corrections(cMin,cMax,centr);
  RawYields(cMin,cMax,centr);
  Spectra(centr);

  fAss->Close();
  fGen->Close();
  fRaw->Close();
}

TH1 *GetEfficiency(Float_t, Float_t, const Char_t *, const Char_t *);
TH1 *GetFeedDown(Float_t, Float_t, TString);

void Corrections(Float_t cmin, Float_t cmax, TString centr) {

  TString name;

  TH1 *effK0s=
  GetEfficiency(cmin, cmax, "fK0sAs","f3dHistPrimRawPtVsYVsMultK0Short");
  name="eff_K0s_";
  name+=centr;
  effK0s->SetName(name.Data());

  TH1 *effLambda=
  GetEfficiency(cmin, cmax, "fLambdaAs","f3dHistPrimRawPtVsYVsMultLambda");
  name="eff_Lambda_";
  name+=centr;
  effLambda->SetName(name.Data());
    TH1 *fdLambda=GetFeedDown(cmin, cmax, centr);
    name="fd_Lambda_";
    name+=centr;
    fdLambda->SetName(name.Data());


  TH1 *effLambdaBar=
  GetEfficiency(cmin,cmax,"fLambdaBarAs","f3dHistPrimRawPtVsYVsMultAntiLambda");
  name="eff_LambdaBar_";
  name+=centr;
  effLambdaBar->SetName(name.Data());

  TFile *f=TFile::Open("SpectraV0CutVariations.root","update");
    effK0s->Write();
    effLambda->Write();    fdLambda->Write();
    effLambdaBar->Write();
  f->Close();
}

TH1 *GetEfficiency(Float_t cMin, Float_t cMax, 
		   const Char_t *chis, const Char_t *znam) {
  Double_t xbins[]={
   0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,
   1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,
   2.2,2.4,2.6,2.8,3.0,3.2,3.4,3.6,3.8,4.0,
   4.5,5.0,5.5,6.5,8.0,10.0,12.0
  };
  Int_t nb=sizeof(xbins)/sizeof(Double_t);
  Int_t nb1=nb-1;

  // Numerator
  fAss->cd();
  TH2F *f2d=(TH2F*)gDirectory->Get(chis); f2d->Sumw2();
  TH1D *hAs=f2d->ProjectionX("hAs",0,-1,"e"); 

  // Denominator
  fGen->cd();
  TH3F *f3d = (TH3F*)gDirectory->Get(znam);
  f3d->Sumw2();

  TH1D *fpt = f3d->ProjectionX("fpt",
                  f3d->GetYaxis()->FindBin(-0.5+1e-2),
                  f3d->GetYaxis()->FindBin(+0.5-1e-2),
                  f3d->GetZaxis()->FindBin(cMin),
                  f3d->GetZaxis()->FindBin(cMax)
              );
   TH1 *hMc=fpt->Rebin(nb1,"hMc",xbins);
  
  //Efficiency 
  TH1 *eff = (TH1*)hAs->Clone();
  eff->Divide(hAs,hMc,1,1,"b");

  return eff;
}

Bool_t GetBinContentError(const TH1 *h, Double_t x, Double_t &c, Double_t &e) {
   Int_t i1=h->GetXaxis()->FindFixBin(x);
   if (i1 <=0) return kFALSE;
   Int_t nb=h->GetNbinsX();
   if (i1 >nb) return kFALSE;

   Double_t x1=h->GetBinCenter(i1);

   if (TMath::Abs(x1-x) < 1e-13) {
     c=h->GetBinContent(i1);
     e=h->GetBinError(i1);
     return kTRUE;
   }

   Int_t i2 = (x1 < x) ? i1+1 : i1-1;
   if (i2 <=0) i2=i1+1;
   if (i2 >nb) i2=i1-1;
   Double_t x2=h->GetBinCenter(i2);

   Double_t c1=h->GetBinContent(i1);
   Double_t c2=h->GetBinContent(i2);
   c = c1 + (x - x1)*(c2 - c1)/(x2 - x1);
   e = 0.5*(h->GetBinError(i1) + h->GetBinError(i2));

   return kTRUE;
}

TH1 *GetFeedDown(Float_t cMin, Float_t cMax, TString cent) {
  // Get the FD matrix
  fAss->cd();
  TH3F *f3d = (TH3F*)gDirectory->Get("fLambdaFromXi");
  f3d->Sumw2();
  TH2D *fdMat = (TH2D*)f3d->Project3D("zxe");
  const Int_t nbx=fdMat->GetNbinsX();
  const Int_t nby=fdMat->GetNbinsY();
  TAxis *xiPtAxis=fdMat->GetYaxis();

  // Re-normalise the FD matrix with MC Xi spectrum
  fGen->cd();
  f3d = (TH3F*)gDirectory->Get("f3dHistGenPtVsYVsMultXiMinus");
  f3d->Sumw2();

  TH1D *hMcXi = f3d->ProjectionX("fpt",
                  f3d->GetYaxis()->FindBin(-0.5+1e-2),
                  f3d->GetYaxis()->FindBin(+0.5-1e-2),
                  f3d->GetZaxis()->FindBin(cMin),
                  f3d->GetZaxis()->FindBin(cMax)
              );

  TH2D *h2McXi=new TH2D(*fdMat); 
  h2McXi->Reset();
  h2McXi->Sumw2();

  for (Int_t i=1; i<=nbx; i++) {
      for (Int_t j=1; j<=nby; j++) {
	  Double_t pt=xiPtAxis->GetBinCenter(j);
          Double_t c=0.,e=0.;
          if (!GetBinContentError(hMcXi,pt,c,e)) continue;
          h2McXi->SetBinContent(i,j,c);
          h2McXi->SetBinError(i,j,e);
      }
  }
  fdMat->Divide(h2McXi);
  delete h2McXi;


  // Multiply the re-normalised matrix with the real Xi spectrum
  TString fileName("systcorrectedpt_cent");
     //*** "Dictionary" ***
     if (cent=="0005") cent="0010";
     if (cent=="6080") cent="6090";
     if (cent=="8090") cent="6090";
  TFile *f=TFile::Open((fileName+cent+".root").Data());
     TH1F *hReXi=(TH1F*)gDirectory->Get("correctedpt_0");

     TH2D *h2ReXi=new TH2D(*fdMat); 
     h2ReXi->Reset();
     h2ReXi->Sumw2();

     for (Int_t i=1; i<=nbx; i++) {
         for (Int_t j=1; j<=nby; j++) {
	    Double_t pt=xiPtAxis->GetBinCenter(j);
            Double_t c=0.,e=0.;
            if (!GetBinContentError(hReXi,pt,c,e)) continue;
            h2ReXi->SetBinContent(i,j,c);
            h2ReXi->SetBinError(i,j,e);
         }
     }
     fdMat->Multiply(h2ReXi);
     delete h2ReXi;
     f->Close();

  TH1 *fd=fdMat->ProjectionX("_px",0,-1,"e");
  fd->Scale(0.1,"width");
  return fd;
}


void RawYields(Float_t cMin, Float_t cMax, TString centr) {
  TString name;
  TH2F *f2d=0;

  //+++ Number of events for normalisation
  TFile *file=TFile::Open("LHC10h_pass2/Merged.root");
    TList *v0list=(TList *)gFile->Get("PWGLFExtractV0_PP/clistV0");
    TH1F *fMult=(TH1F*)v0list->FindObject("fHistMultiplicity");
    Int_t i1=fMult->GetXaxis()->FindBin(cMin+1e-2);
    Int_t i2=fMult->GetXaxis()->FindBin(cMax+1e-2);
    Float_t nEvents=fMult->Integral(i1,i2);
  file->Close();

  fRaw->cd();

  name="raw_K0s_";
  name+=centr;  
  f2d=(TH2F*)gDirectory->Get("fK0sSi"); f2d->Sumw2();
  TH1D *rawK0s=f2d->ProjectionX(name,0,-1,"e");
  rawK0s->Scale(1/nEvents,"width");
 
  name="raw_Lambda_";
  name+=centr;  
  f2d=(TH2F*)gDirectory->Get("fLambdaSi"); f2d->Sumw2();
  TH1D *rawLambda=f2d->ProjectionX(name,0,-1,"e");
  rawLambda->Scale(1/nEvents,"width");
 
  name="raw_LambdaBar_";
  name+=centr;  
  f2d=(TH2F*)gDirectory->Get("fLambdaBarSi"); f2d->Sumw2();
  TH1D *rawLambdaBar=f2d->ProjectionX(name,0,-1,"e");
  rawLambdaBar->Scale(1/nEvents,"width");
 
  TFile *f=TFile::Open("SpectraV0CutVariations.root","update");
    rawK0s->Write();
    rawLambda->Write();
    rawLambdaBar->Write();
  f->Close();
}

void Spectra(TString centr) {
  TString name;
  TH1 *eff=0;
  TH1D *raw=0;
  TH1D *spe=0;

  TFile *f=TFile::Open("SpectraV0CutVariations.root","update");
    name="eff_K0s_";
    name+=centr;
    eff = (TH1*)gDirectory->Get(name.Data());
    name="raw_K0s_";
    name+=centr;
    raw = (TH1D*)gDirectory->Get(name.Data());
    spe = new TH1D(*raw);
    spe->Divide(eff);
    name="K0s_";
    name+=centr;
    spe->SetName(name.Data());
    spe->Write();  

    name="eff_Lambda_";
    name+=centr;
    eff = (TH1*)gDirectory->Get(name.Data());
    name="raw_Lambda_";
    name+=centr;
    raw = (TH1D*)gDirectory->Get(name.Data());
    spe = new TH1D(*raw);
    spe->Divide(eff);
    name="Lambda_";
    name+=centr;
    spe->SetName(name.Data());
    spe->Write();  

    name="eff_LambdaBar_";
    name+=centr;
    eff = (TH1*)gDirectory->Get(name.Data());
    name="raw_LambdaBar_";
    name+=centr;
    raw = (TH1D*)gDirectory->Get(name.Data());
    spe = new TH1D(*raw);
    spe->Divide(eff);
    name="LambdaBar_";
    name+=centr;
    spe->SetName(name.Data());
    spe->Write();  
  f->Close();
}

void Generated() {
  TList *v0listMC=0;
  TH3F *h3=0;
  //
  TFile::Open("LHC11a10b_plus/Merged.root");  
  v0listMC=(TList *)gFile->Get("PWGLFExtractPerformanceV0_PP_MC/clistV0MC");

  TH3F *k0s = 
    (TH3F*)v0listMC->FindObject("f3dHistPrimRawPtVsYVsMultK0Short");
  TH3F *k0s_nonInj = 
    (TH3F*)v0listMC->FindObject("f3dHistPrimRawPtVsYVsMultNonInjK0Short");

  TH3F *lam = 
    (TH3F*)v0listMC->FindObject("f3dHistPrimRawPtVsYVsMultLambda");
  TH3F *lam_nonInj = 
    (TH3F*)v0listMC->FindObject("f3dHistPrimRawPtVsYVsMultNonInjLambda");

  TH3F *alam = 
    (TH3F*)v0listMC->FindObject("f3dHistPrimRawPtVsYVsMultAntiLambda");
  TH3F *alam_nonInj = 
    (TH3F*)v0listMC->FindObject("f3dHistPrimRawPtVsYVsMultNonInjAntiLambda");

  TH3F *xiMinus = 
    (TH3F*)v0listMC->FindObject("f3dHistGenPtVsYVsMultXiMinus");
  TH3F *xiPlus = 
    (TH3F*)v0listMC->FindObject("f3dHistGenPtVsYVsMultXiPlus");


  //    
  TFile::Open("LHC11a10b_bis/Merged.root");  
  v0listMC=(TList *)gFile->Get("PWGLFExtractPerformanceV0_PP_MC/clistV0MC");

  h3 = (TH3F*)v0listMC->FindObject("f3dHistPrimRawPtVsYVsMultK0Short");
  k0s->Add(h3);
  h3 = (TH3F*)v0listMC->FindObject("f3dHistPrimRawPtVsYVsMultNonInjK0Short");
  k0s_nonInj->Add(h3); 

  h3 = (TH3F*)v0listMC->FindObject("f3dHistPrimRawPtVsYVsMultLambda");
  lam->Add(h3); 
  h3 = (TH3F*)v0listMC->FindObject("f3dHistPrimRawPtVsYVsMultNonInjLambda");
  lam_nonInj->Add(h3); 

  h3 = (TH3F*)v0listMC->FindObject("f3dHistPrimRawPtVsYVsMultAntiLambda");
  alam->Add(h3); 
  h3 = (TH3F*)v0listMC->FindObject("f3dHistPrimRawPtVsYVsMultNonInjAntiLambda");
  alam_nonInj->Add(h3); 

  h3 = (TH3F*)v0listMC->FindObject("f3dHistGenPtVsYVsMultXiMinus");
  xiMinus->Add(h3);
  h3 = (TH3F*)v0listMC->FindObject("f3dHistGenPtVsYVsMultXiPlus");
  xiPlus->Add(h3); 


  //
  TFile::Open("LHC11a10a_bis/Merged.root");
  v0listMC=(TList *)gFile->Get("PWGLFExtractPerformanceV0_PP_MC/clistV0MC");

  h3 = (TH3F*)v0listMC->FindObject("f3dHistPrimRawPtVsYVsMultK0Short");
  k0s_nonInj->Add(h3); 
  h3 = (TH3F*)v0listMC->FindObject("f3dHistPrimRawPtVsYVsMultLambda");
  lam_nonInj->Add(h3); 
  h3 = (TH3F*)v0listMC->FindObject("f3dHistPrimRawPtVsYVsMultAntiLambda");
  alam_nonInj->Add(h3); 
  h3 = (TH3F*)v0listMC->FindObject("f3dHistGenPtVsYVsMultXiMinus");
  xiMinus->Add(h3);
  h3 = (TH3F*)v0listMC->FindObject("f3dHistGenPtVsYVsMultXiPlus");
  xiPlus->Add(h3); 

  TFile *f=TFile::Open("Generated.root","new");
  k0s->Write(); k0s_nonInj->Write();
  lam->Write(); lam_nonInj->Write();
  alam->Write(); alam_nonInj->Write();
  xiMinus->Write();
  xiPlus->Write();
  f->Close();
}

