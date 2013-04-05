/* 
Example macro which shows how to use the AliBlastwaveFit code.
Prepared by F. Noferini (fnoferin@cern.ch)

It assume you have pi, K, proton spectra and v2 and it perform a global
blastwave fit in the selected pt range for the requested species.

Centrality is taken as an argument, while all the other options are set
in common variables.

Here: stat and syst errors considered.
*/

Bool_t kPiPlusSpectra=kTRUE; // request fit to pion spectra
Bool_t kKaPlusSpectra=kTRUE; // request fit to kaon spectra
Bool_t kPrPlusSpectra=kTRUE; // request fit to proton spectra
Bool_t kPiV2=kTRUE;          // request fit to pion v2
Bool_t kKaV2=kTRUE;          // request fit to kaon v2
Bool_t kPrV2=kTRUE;          // request fit to antiproton v2

Bool_t kCont = kFALSE;       // request the contour at the end (it needs time)

TH1D *hpiplus;
TH1D *hkaplus;
TH1D *hprplus;

TGraphAsymmErrors *gpiv2;
TGraphAsymmErrors *gkav2;
TGraphAsymmErrors *gprv2;

// particle masses
Float_t mPi = 1.39570000000000000e-01;
Float_t mKa = 4.93676999999999977e-01;
Float_t mPr = 9.38271999999999995e-01;

// centralities
Int_t cmin[] = {0,5,10,20,30,40,50,60,70,80};
Int_t cmax[] = {5,10,20,30,40,50,60,70,80,100};

// pt range for the fit
Float_t ptMinPi = 0.5;
Float_t ptMaxPi = 1.0;
Float_t ptMinKa = 0.2;
Float_t ptMaxKa = 1.5;
Float_t ptMinPr = 0.3;
Float_t ptMaxPr = 3.0;

Bool_t kLoaded = kFALSE;

fitblastwave(Int_t ic=2){

  if(! kLoaded) LoadLib();

  char name[200];

  // Get Spectra
  TFile *f = new TFile("SPECTRA_COMB_20120809_default.root");
  if(f){
    sprintf(name,"cent%i_pion_plus",ic);
    hpiplus = (TH1D *) f->Get(name);
    
    if(! hpiplus) kPiPlusSpectra = kFALSE; 
    
    sprintf(name,"cent%i_kaon_plus",ic);
    hkaplus = (TH1D *) f->Get(name);
    
    if(! hkaplus) kKaPlusSpectra = kFALSE; 
    
    sprintf(name,"cent%i_proton_plus",ic);
    hprplus = (TH1D *) f->Get(name);
 
    if(! hprplus) kPrPlusSpectra = kFALSE; 
  }

  TCanvas *csp = new TCanvas("cspectra","cspectra");
  csp->Divide(3,1);
  csp->cd(1)->SetLogy();
  if(hpiplus) hpiplus->Draw();
  csp->cd(2)->SetLogy();
  if(hkaplus) hkaplus->Draw();
  csp->cd(3)->SetLogy();
  if(hprplus) hprplus->Draw();

  // Get v2
  sprintf(name,"v2/v2SP_pion_%02i_%02i.txt",cmin[ic],cmax[ic]);
  gpiv2 = GetGraph(name);
  if(! gpiv2) kPiV2 = kFALSE;
  sprintf(name,"v2/v2SP_kaon_%02i_%02i.txt",cmin[ic],cmax[ic]);
  gkav2 = GetGraph(name);
  if(! gkav2) kKaV2 = kFALSE;
  sprintf(name,"v2/v2SP_antipr_%02i_%02i.txt",cmin[ic],cmax[ic]);
  gprv2 = GetGraph(name);
  if(! gprv2) kPrV2 = kFALSE;

  TCanvas *cv2 = new TCanvas("cv2","cv2");
  cv2->Divide(3,1);
  cv2->cd(1);
  if(gpiv2) gpiv2->Draw("AP");
  cv2->cd(2);
  if(gkav2) gkav2->Draw("AP");
  cv2->cd(3);
  if(gprv2) gprv2->Draw("AP");

  // initialize fitter
  AliBlastwaveFit2D *bwPi = new AliBlastwaveFit2D("pions",mPi);
  bwPi->SetMinPt(ptMinPi);
  bwPi->SetMaxPt(ptMaxPi);

  AliBlastwaveFit2D *bwKa = new AliBlastwaveFit2D("kaons",mKa);
  bwKa->SetMinPt(ptMinKa);
  bwKa->SetMaxPt(ptMaxKa);

  AliBlastwaveFit2D *bwPr = new AliBlastwaveFit2D("protons",mPr);
  bwPr->SetMinPt(ptMinPr);
  bwPr->SetMaxPt(ptMaxPr);

  if(kPiPlusSpectra) bwPi->SetSpectrumObj(hpiplus);
  if(kPiV2) bwPi->SetV2Obj(gpiv2);
  if(kKaPlusSpectra) bwKa->SetSpectrumObj(hkaplus);
  if(kKaV2) bwKa->SetV2Obj(gkav2);
  if(kPrPlusSpectra) bwPr->SetSpectrumObj(hprplus);
  if(kPrV2) bwPr->SetV2Obj(gprv2);

  AliBlastwaveFitter *fitter = new AliBlastwaveFitter("fitterBW");
  fitter->AddFitFunction(bwPi);
  fitter->AddFitFunction(bwKa);
  fitter->AddFitFunction(bwPr);

  // go to fit
  fitter->PrepareToFit(); // initialize the fitter object
  fitter->Fit();          // perform the fit (it will take some time)

  // Draw fit
  if(hpiplus){
    csp->cd(1);
    bwPi->GetSpectraFit()->Draw("SAME");
  }
  if(hkaplus){
    csp->cd(2);
    bwKa->GetSpectraFit()->Draw("SAME");
  }
  if(hprplus){
    csp->cd(3);
    bwPr->GetSpectraFit()->Draw("SAME");
  }

  if(gpiv2){
    cv2->cd(1);
    bwPi->GetV2Fit()->Draw("SAME");
  }
  if(gkav2){
    cv2->cd(2);
    bwKa->GetV2Fit()->Draw("SAME");
  }
  if(gprv2){
    cv2->cd(3);
    bwPr->GetV2Fit()->Draw("SAME");
  }

  // Print some outputs
  printf("Chi2 = %f\n",fitter->GetChi2());
  printf("N.D.G.F. = %f\n",fitter->GetNDGF());
  printf("<beta> = %f\n",bwPi->GetMeanBeta());

  // Do the contour if requested
  if(kCont){
    TCanvas *cCont = new TCanvas("cContour","cContour");
    cCont->Divide(3,1);
    cCont->cd(1);
    TGraph *gCont1 = fitter->DoContour(50,2,0,1); // rho and T
    gCont1->Draw("AP");
    cCont->cd(2);
    TGraph *gCont2 = fitter->DoContourBetaT(50,2,0,1); // beta and T
    gCont2->Draw("AP");
    cCont->cd(3);
    TGraph *gCont3 = fitter->DoContour(50,1,3,1); // s2 and rho2
    gCont3->Draw("AP");
  }
}

LoadLib(){
  gSystem->Load("libVMC.so");
  gSystem->Load("libPhysics.so");
  gSystem->Load("libTree.so");
  gSystem->Load("libMinuit.so");
  gSystem->Load("libSTEERBase.so");
  gSystem->Load("libPWGCFflowBW.so");
  kLoaded = kTRUE;
}


TGraphAsymmErrors *GetGraph(const char *filename){
  char name[200];

  sprintf(name,"cat %s|grep -v BinCenter >/tmp/tempov2",filename);

  system(name);

  Float_t x,y,statLow,statHigh,systLow,systHigh;

  Int_t np=0;
  Double_t pt[100],pterr[100],v2[100],v2errLow[100],v2errHigh[100];

  FILE *f = fopen("/tmp/tempov2","r");
  while(fscanf(f,"%f %f %f %f %f %f",&x,&y,&statHigh,&statLow,&systHigh,&systLow) == 6){
    pt[np] = x;
    pterr[np] = 0;
    v2[np] = y;
    v2errLow[np] = TMath::Sqrt(statLow*statLow + systLow*systLow);
    v2errHigh[np] = TMath::Sqrt(statHigh*statHigh + systHigh*systHigh);
    np++;
  }
  fclose(f);

  TGraphAsymmErrors *g = new TGraphAsymmErrors(np,pt,v2,pterr,pterr,v2errLow,v2errHigh);

  return g;
}


Float_t ComputePiV2int(Int_t ic=2){

  if(! kLoaded) LoadLib();

  char name[200];

  // Get Spectra
  TFile *f = new TFile("SPECTRA_COMB_20120809_default.root");
  if(f){
    sprintf(name,"cent%i_pion_plus",ic);
    hpiplus = (TH1D *) f->Get(name);
    
    if(! hpiplus) return 0.0; 
  }


  TCanvas *csp = new TCanvas("cpiSpectra","cpiSpectra");
  hpiplus->Draw();

  // Get v2
  sprintf(name,"v2/v2QC4_pion_%02i_%02i.txt",cmin[ic],cmax[ic]);
  gpiv2 = GetGraph(name);
  if(! gpiv2) return 0.0;
  

  TCanvas *cv2 = new TCanvas("cpiV2","cpiV2");
  gpiv2->Draw("AP");
  
  // initialize fitter
  AliBlastwaveFit2D *bwPi = new AliBlastwaveFit2D("pionsSp",mPi);
  bwPi->SetMinPt(0.1);
  bwPi->SetMaxPt(2.0);

  AliBlastwaveFit2D *bwPi2 = new AliBlastwaveFit2D("pionsV2",mPi);
  bwPi2->SetMinPt(0.2);
  bwPi2->SetMaxPt(1.0);

  AliBlastwaveFitSpectra *bwPi3 = new AliBlastwaveFitSpectra("pionsSpHP",mPi); // use only spectra function to avoid overwriting of other fit (parameters are static memebers)
  bwPi3->SetMinPt(1.5);
  bwPi3->SetMaxPt(3.0);

  bwPi->SetSpectrumObj(hpiplus);

  bwPi2->SetV2Obj(gpiv2);

  bwPi3->SetSpectrumObj(hpiplus);

  AliBlastwaveFitter *fitter = new AliBlastwaveFitter("fitterPion");
  fitter->AddFitFunction(bwPi);
  fitter->AddFitFunction(bwPi2);

  AliBlastwaveFitter *fitter2 = new AliBlastwaveFitter("fitterPion2");
  fitter2->AddFitFunction(bwPi3);

  // go to fit
  fitter->PrepareToFit(); // initialize the fitter object
  fitter->Fit();          // perform the fit (it will take some time)
 
  // Print some outputs
  printf("Chi2 = %f\n",fitter->GetChi2());
  printf("N.D.G.F. = %f\n",fitter->GetNDGF());
  printf("<beta> = %f\n",bwPi->GetMeanBeta());

  fitter2->PrepareToFit(); // initialize the fitter object
  fitter2->Fit();          // perform the fit (it will take some time)

  TF1 *fitSP = bwPi->GetSpectraFit();
  TF1 *fitV2 = bwPi2->GetV2Fit();
  TF1 *fitSP2 = bwPi3->GetSpectraFit();

  printf("2)Chi2 = %f\n",fitter2->GetChi2());
  printf("2)N.D.G.F. = %f\n",fitter2->GetNDGF());
  printf("2)<beta> = %f\n",bwPi3->GetMeanBeta());

  // Draw fit
  csp->cd();
  fitSP->Draw("SAME");
  fitSP2->Draw("SAME");
  cv2->cd();
  fitV2->Draw("SAME");

  fitSP->SetRange(0.0001,1.3);
  fitSP2->SetRange(1.5,6);
  fitV2->SetRange(0.0001,1);
  
  Float_t num = 0;
  Float_t den = 0;

  Float_t num1 = 0;
  Float_t den1 = 0;

  Float_t num2 = 0;
  Float_t den2 = 0;

  for(Int_t i=0;i<10;i++){ // form 0 to 0.2
    Float_t x1 = i*0.02;
    Float_t x2 = (i+1)*0.02;
    Float_t xm = (x1+x2)*0.5;

    Float_t yield = fitSP->Integral(x1,x2);
    Float_t v2 = fitV2->Eval(xm);

    if(xm > 0.2){
      den += yield;
      num += yield * v2;
      
      den1 += yield*(1 + 0.0202/xm + 0.03 * xm);
      num1 += yield*(1 + 0.0202/xm + 0.03 * xm) * v2 * (1 - 0.05);
      
      den2 += yield*(1 - 0.0202/xm - 0.03 * xm);
      num2 += yield*(1 - 0.0202/xm - 0.03 * xm) * v2 * (1 + 0.05);
    }
  }

  for(Int_t i=0; i < gpiv2->GetN();i++){
    Float_t x = gpiv2->GetX()[i];

    Float_t frSyst = 1 - 2*(x-0.2)/(3-0.2);

    if(x > 0.2 && x < 5){
      Float_t binwidth = 0.05;
      if(x < 3) binwidth = 0.05;
      else if(x < 4) binwidth = 0.1;
      else binwidth = 0.2;

      Float_t yield = hpiplus->Interpolate(x) * 2 * binwidth;
      if(x>3) yield = fitSP2->Integral(x-binwidth,x+binwidth);
      Float_t v2 = gpiv2->GetY()[i];
      Float_t v2err1 = gpiv2->GetEYlow()[i];
      Float_t v2err2 = gpiv2->GetEYhigh()[i];

      den += yield;
      num += yield * v2;

      den1 += yield*(1 + frSyst*(0.0202/x + 0.03 * x));
      num1 += yield*(1 + frSyst*(0.0202/x + 0.03 * x)) * v2 * (1 - v2err1);

      den2 += yield*(1 - frSyst*(0.0202/x + 0.03 * x));
      num2 += yield*(1 - frSyst*(0.0202/x + 0.03 * x)) * v2 * (1 + v2err2);

      printf("pt<%f) v2int = %f\n",x+binwidth,num/den);
    }
  }

  printf("Integrated flow for pions (0 < p_T < 6 GeV/c) = %f\n",num/den);
  printf("Syst. = %f\n",(num2/den2 - num1/den1)/2);
  printf("Yield = %f\n",den);
}

Float_t ComputeKaV2int(Int_t ic=2){

  if(! kLoaded) LoadLib();

  char name[200];

  // Get Spectra
  TFile *f = new TFile("SPECTRA_COMB_20120809_default.root");
  if(f){
    sprintf(name,"cent%i_kaon_plus",ic);
    hkaplus = (TH1D *) f->Get(name);
    
    if(! hkaplus) return 0.0; 
  }


  TCanvas *csp = new TCanvas("ckaSpectra","ckaSpectra");
  hkaplus->Draw();

  // Get v2
  sprintf(name,"v2/v2QC4_kaon_%02i_%02i.txt",cmin[ic],cmax[ic]);
  gkav2 = GetGraph(name);
  if(! gkav2) return 0.0;
  

  TCanvas *cv2 = new TCanvas("ckaV2","ckaV2");
  gkav2->Draw("AP");
  
  // initialize fitter
  AliBlastwaveFit2D *bwKa = new AliBlastwaveFit2D("kaonsSp",mKa);
  bwKa->SetMinPt(0.2);
  bwKa->SetMaxPt(1.3);

  AliBlastwaveFit2D *bwKa2 = new AliBlastwaveFit2D("kaonsV2",mKa);
  bwKa2->SetMinPt(0.25);
  bwKa2->SetMaxPt(1.0);


  AliBlastwaveFitSpectra *bwKa3 = new AliBlastwaveFitSpectra("kaonsSpHP",mKa); // use only spectra function to avoid overwriting of other fit (parameters are static memebers)
  bwKa3->SetMinPt(1.5);
  bwKa3->SetMaxPt(3.0);

  bwKa->SetSpectrumObj(hkaplus);

  bwKa2->SetV2Obj(gkav2);

  bwKa3->SetSpectrumObj(hkaplus);

  AliBlastwaveFitter *fitter = new AliBlastwaveFitter("fitterKaon");
  fitter->AddFitFunction(bwKa);
  fitter->AddFitFunction(bwKa2);

  AliBlastwaveFitter *fitter2 = new AliBlastwaveFitter("fitterKaon2");
  fitter2->AddFitFunction(bwKa3);

  // go to fit
  fitter->PrepareToFit(); // initialize the fitter object
  fitter->Fit();          // perform the fit (it will take some time)

  // Print some outputs
  printf("Chi2 = %f\n",fitter->GetChi2());
  printf("N.D.G.F. = %f\n",fitter->GetNDGF());
  printf("<beta> = %f\n",bwKa->GetMeanBeta());

  fitter2->PrepareToFit(); // initialize the fitter object
  fitter2->Fit();          // perform the fit (it will take some time)

  TF1 *fitSP = bwKa->GetSpectraFit();
  TF1 *fitV2 = bwKa2->GetV2Fit();
  TF1 *fitSP2 = bwKa3->GetSpectraFit();

  printf("2)Chi2 = %f\n",fitter2->GetChi2());
  printf("2)N.D.G.F. = %f\n",fitter2->GetNDGF());
  printf("2)<beta> = %f\n",bwKa3->GetMeanBeta());

  // Draw fit
  csp->cd();
  fitSP->Draw("SAME");
  fitSP2->Draw("SAME");
  cv2->cd();
  fitV2->Draw("SAME");

  fitSP->SetRange(0.0001,1.3);
  fitSP2->SetRange(1.5,6);
  fitV2->SetRange(0.0001,1);
  
  Float_t num = 0;
  Float_t den = 0;

  Float_t num1 = 0;
  Float_t den1 = 0;

  Float_t num2 = 0;
  Float_t den2 = 0;

  for(Int_t i=0;i<10;i++){ // form 0 to 0.2
    Float_t x1 = i*0.025;
    Float_t x2 = (i+1)*0.025;
    Float_t xm = (x1+x2)*0.5;

    Float_t yield = fitSP->Integral(x1,x2);
    Float_t v2 = fitV2->Eval(xm);

    if(xm > 0.2){
      den += yield;
      num += yield * v2;
      
      den1 += yield*(1 + 0.0215/xm + 0.05 * xm);
      num1 += yield*(1 + 0.0215/xm + 0.05 * xm) * v2 * (1 - 0.1);
      
      den2 += yield*(1 - 0.0215/xm - 0.05 * xm);
      num2 += yield*(1 - 0.0215/xm - 0.05 * xm) * v2 * (1 + 0.1);    
    }
  }

  for(Int_t i=0; i < gkav2->GetN();i++){
    Float_t x = gkav2->GetX()[i];
    Float_t frSyst = 1 - 2*(x-0.25)/(3-0.25);

    if(x > 0.25 && x < 5){
      Float_t binwidth = 0.05;
      if(x < 3) binwidth = 0.05;
      else if(x < 4) binwidth = 0.1;
      else binwidth = 0.2;

      Float_t yield = hkaplus->Interpolate(x) * 2 * binwidth;
      if(x>3) yield = fitSP2->Integral(x-binwidth,x+binwidth);
      Float_t v2 = gkav2->GetY()[i];
      Float_t v2err1 = gkav2->GetEYlow()[i];
      Float_t v2err2 = gkav2->GetEYhigh()[i];

      den += yield;
      num += yield * v2;

      den1 += yield*(1 + frSyst*(0.0215/x + 0.05 * x));
      num1 += yield*(1 + frSyst*(0.0215/x + 0.05 * x)) * v2 * (1 - v2err1);

      den2 += yield*(1 - frSyst*(0.0215/x + 0.05 * x));
      num2 += yield*(1 - frSyst*(0.0215/x + 0.05 * x)) * v2 * (1 + v2err2);

      printf("pt<%f) v2int = %f\n",x+binwidth,num/den);
    }
  }

  printf("Integrated flow for kaons (0 < p_T < 6 GeV/c) = %f\n",num/den);
  printf("Syst. = %f\n",(num2/den2 - num1/den1)/2);
  printf("Yield = %f\n",den);
}

Float_t ComputePrV2int(Int_t ic=2){

  if(! kLoaded) LoadLib();

  char name[200];

  // Get Spectra
  TFile *f = new TFile("SPECTRA_COMB_20120809_default.root");
  if(f){
    sprintf(name,"cent%i_proton_minus",ic);
    hprplus = (TH1D *) f->Get(name);
    
    if(! hprplus) return 0.0; 
  }


  TCanvas *csp = new TCanvas("cprSpectra","cprSpectra");
  hprplus->Draw();

  // Get v2
  sprintf(name,"v2/v2QC4_antipr_%02i_%02i.txt",cmin[ic],cmax[ic]);
  gprv2 = GetGraph(name);
  if(! gprv2) return 0.0;
  

  TCanvas *cv2 = new TCanvas("cprV2","cprV2");
  gprv2->Draw("AP");
  
  // initialize fitter
  AliBlastwaveFit2D *bwPr = new AliBlastwaveFit2D("antiprotonsSp",mPr);
  bwPr->SetMinPt(0.3);
  bwPr->SetMaxPt(1.2);

  AliBlastwaveFit2D *bwPr2 = new AliBlastwaveFit2D("antiprotonsV2",mPr);
  bwPr2->SetMinPt(0.3);
  bwPr2->SetMaxPt(2.0);

  AliBlastwaveFitSpectra *bwPr3 = new AliBlastwaveFitSpectra("antiprotonsSpHP",mPr); // use only spectra function to avoid overwriting of other fit (parameters are static memebers)
  bwPr3->SetMinPt(2.5);
  bwPr3->SetMaxPt(4.5);

  bwPr->SetSpectrumObj(hprplus);

  bwPr2->SetV2Obj(gprv2);

  bwPr3->SetSpectrumObj(hprplus);

  AliBlastwaveFitter *fitter = new AliBlastwaveFitter("fitterProton");
  fitter->AddFitFunction(bwPr);
  fitter->AddFitFunction(bwPr2);

  AliBlastwaveFitter *fitter2 = new AliBlastwaveFitter("fitterProton2");
  fitter2->AddFitFunction(bwPr3);

  // go to fit
  fitter->PrepareToFit(); // initialize the fitter object
  fitter->Fit();          // perform the fit (it will take some time)

  // Print some outputs
  printf("Chi2 = %f\n",fitter->GetChi2());
  printf("N.D.G.F. = %f\n",fitter->GetNDGF());
  printf("<beta> = %f\n",bwPr->GetMeanBeta());

  fitter2->PrepareToFit(); // initialize the fitter object
  fitter2->Fit();          // perform the fit (it will take some time)

  TF1 *fitSP = bwPr->GetSpectraFit();
  TF1 *fitV2 = bwPr2->GetV2Fit();
  TF1 *fitSP2 = bwPr3->GetSpectraFit();

  printf("2)Chi2 = %f\n",fitter2->GetChi2());
  printf("2)N.D.G.F. = %f\n",fitter2->GetNDGF());
  printf("2)<beta> = %f\n",bwPr3->GetMeanBeta());

  // Draw fit
  csp->cd();
  fitSP->Draw("SAME");
  fitSP2->Draw("SAME");
  cv2->cd();
  fitV2->Draw("SAME");

  fitSP->SetRange(0.0001,1.2);
  fitSP2->SetRange(2.5,6);
  fitV2->SetRange(0.0001,2);
  
  Float_t num = 0;
  Float_t den = 0;

  Float_t num1 = 0;
  Float_t den1 = 0;

  Float_t num2 = 0;
  Float_t den2 = 0;

  for(Int_t i=0;i<10;i++){ // form 0 to 0.2
    Float_t x1 = i*0.03;
    Float_t x2 = (i+1)*0.03;
    Float_t xm = (x1+x2)*0.5;

    Float_t yield = fitSP->Integral(x1,x2);
    Float_t v2 = fitV2->Eval(xm);

    if(xm > 0.2){
      den += yield;
      num += yield * v2;
      
      den1 += yield*(1 + 0.064/xm + 0.0083 * xm * xm);
      num1 += yield*(1 + 0.064/xm + 0.0083 * xm * xm) * v2 * (1 - 0.2);
      
      den2 += yield*(1 - 0.064/xm - 0.0083 * xm * xm);
      num2 += yield*(1 - 0.064/xm - 0.0083 * xm * xm) * v2 * (1 + 0.2);
    }
  }

  for(Int_t i=0; i < gprv2->GetN();i++){
    Float_t x = gprv2->GetX()[i];
    Float_t frSyst = 1 - 2*(x-0.3)/(4.5-0.3);

    if(x > 0.3 && x < 5){
      Float_t binwidth = 0.05;
      if(x < 3) binwidth = 0.05;
      else if(x < 4) binwidth = 0.1;
      else binwidth = 0.2;

      Float_t yield = hprplus->Interpolate(x) * 2 * binwidth;
      if(x>4.5) yield = fitSP2->Integral(x-binwidth,x+binwidth);
      Float_t v2 = gprv2->GetY()[i];
      Float_t v2err1 = gprv2->GetEYlow()[i];
      Float_t v2err2 = gprv2->GetEYhigh()[i];

      den += yield;
      num += yield * v2;

      den1 += yield*(1 + frSyst*(0.064/x + 0.0083 * x * x));
      num1 += yield*(1 + frSyst*(0.064/x + 0.0083 * x * x)) * v2 * (1 - v2err1);

      den2 += yield*(1 - frSyst*(0.064/x + 0.0083 * x * x));
      num2 += yield*(1 - frSyst*(0.064/x + 0.0083 * x * x)) * v2 * (1 + v2err2);

      printf("pt<%f) v2int = %f\n",x+binwidth,num/den);
    }
  }

  printf("Integrated flow for antiprotons (0 < p_T < 6 GeV/c) = %f\n",num/den);
  printf("Syst. = %f\n",(num2/den2 - num1/den1)/2);
  printf("Yield = %f\n",den);
}
