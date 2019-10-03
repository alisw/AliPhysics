#include"TF1.h"
#include"TH1D.h"
#include"TH2F.h"
#include"TMath.h"
#include"TSystem.h"
#include"TCanvas.h"
#include"TFile.h"
#include"TGraphErrors.h"
#include"AliPIDperfContainer.h"
#include"TRandom.h"

Int_t LoadLib();
void doeffPi(Int_t pos=1,Float_t prob=0.1,Float_t etaminkp=-0.8,Float_t etamaxkp=0.8);
TH2F *GetHistoPip(Float_t pt=1,Float_t ptM=1.1,Float_t pMinkp=0,Float_t pMinkn=0.,Float_t etaminkp=-0.8,Float_t etamaxkp=0.8);
TH2F *GetHistoPin(Float_t pt=1,Float_t ptM=1.1,Float_t pMinkn=0,Float_t pMinkp=0.,Float_t etaminkp=-0.8,Float_t etamaxkp=0.8);
void fit(TH1D *h,Float_t *a=NULL,char *opt="",char *opt2="",Float_t pt=1.5);
void AddHisto(TH2F *h1,TH2F *h2,Float_t w);

TObject* fContPid1;
TObject* fContPid2;
const Int_t nBinPid = 14; // pt,eta, ptPip, ptPin, PPip, PPin, TOF3sigmaPip, TOF3sigmaPin, isPhiTrue, nsigmaPip, nsigmaPin
// 0.985 < mass < 1.045 (60) and 0 < centrality < 100 (10)
Int_t binPid[nBinPid] = {1/*ptPhi*/,1/*EtaPi*/,20/*pt+*/,20/*pt-*/,5/*P+*/,5/*P-*/,2/*TOFmatch+*/,2/*TOFmatch-*/,2/*istrue*/,4/*Nsigma+*/,4/*Nsigma-*/,1/*DeltaPhi+*/,1/*DeltaPhi-*/,1/*Psi*/};
Float_t xmin[nBinPid] = {1,-0.8,0.3,0.3,0,0,-0.5,-0.5,-0.5,0,0,-TMath::Pi(),-TMath::Pi(),-TMath::Pi()/2};
Float_t xmax[nBinPid] = {5,0.8,4.3,4.3,1,1,1.5,1.5,1.5,7.5,7.5,TMath::Pi(),TMath::Pi(),TMath::Pi()/2};

TF1 *fsign;
TF1 *fall;
TF1 *fback;

Int_t ifunc=0;

Float_t fitmin = 0.3;
Float_t fitmax = 0.7;

Int_t cmin = 1; // min 1
Int_t cmax = 10;//max 10

Float_t weightS = -1.;

Int_t rebinsize = 1;

Int_t parplotted = 2;

Bool_t isMC = kFALSE; // don't change this (is set automatically)
Bool_t selectTrue = kTRUE; // put it to true to remove background (only for MC)
Bool_t keepTrue = kTRUE; // put it to false to fit only background (only for MC)

Bool_t kGoodMatch = kFALSE; // to check good matching

Bool_t kSigma2vs3 = kFALSE; // to check good matching

Bool_t require5sigma = kFALSE; // don't touch this flag

Bool_t bayesVsigma = kFALSE; // only to do checks

Bool_t kTOFmatch = kFALSE; // for combined PID requires TOF matching

Bool_t kOverAll = kFALSE;
Bool_t kOverAllTOFmatch = kFALSE;
Bool_t kOverAll2Sigma = kFALSE;

Bool_t kPid2Sigma = kFALSE;
Bool_t kPid3Sigma = kFALSE;

TH2F *hmatched;
TH2F *htracked;

Bool_t kLoaded=kFALSE;

Int_t LoadLib(){
  weightS = -1.;

  require5sigma = kFALSE;

  if(! kLoaded){
    gSystem->Load("libVMC");
    gSystem->Load("libPhysics");
    gSystem->Load("libTree");
    gSystem->Load("libMinuit");
    gSystem->Load("libSTEERBase");
    gSystem->Load("libANALYSIS");
    gSystem->Load("libAOD");
    gSystem->Load("libESD");
    gSystem->Load("libANALYSIS");
    gSystem->Load("libANALYSISalice");
    gSystem->Load("libCORRFW");
    gSystem->Load("libNetx");
    gSystem->Load("libPWGPPpid");

    TFile *f = new TFile("AnalysisResults.root");
    TList *l = (TList *) f->Get("contK0sBayes1");
    TList *l2 = (TList *) f->Get("contK0sBayes2");

    if(! (l && l2)) return 0;

    fContPid1 = (AliPIDperfContainer *) l->FindObject("contPID");
    fContPid2 = (AliPIDperfContainer *) l->FindObject("contPID2");
    hmatched = (TH2F *) l2->FindObject("hMatchPi"); 
    htracked = (TH2F *) l2->FindObject("hTrackingPi"); 
  }
  kLoaded = kTRUE;

  // check if MC
  Float_t x[] = {xmin[0]+0.001,xmin[1]+0.001,xmin[2]+0.001,xmin[3]+0.001,xmin[4]+0.001,xmin[5]+0.001,xmin[6]+0.001,xmin[7]+0.001,1/*trueMC*/,xmin[9],xmin[10],xmin[11],xmin[12],xmin[13]};
  Float_t x2[] = {xmax[0],xmax[1],xmax[2],xmax[3],xmax[4],xmax[5],xmax[6],xmax[7],xmax[8],xmax[9],xmax[10],xmax[11],xmax[12],xmax[13]};

  AliPIDperfContainer *tmp = (AliPIDperfContainer *) fContPid1;
  TH1D *h = tmp->GetQA(0, x, x2)->ProjectionX("checkMC");

  if(h->GetEntries()) isMC = kTRUE;
  else isMC=kFALSE;

  if(!isMC){
    selectTrue = kFALSE;
    keepTrue = kTRUE;
  }
  else{
    printf("MC truth found!!!!!!\nIt is MC!!!!!!");
  }

  fsign = new TF1("fsign","gaus(0) +0.5*[0]*TMath::Exp(-[3]*TMath::Abs(x-[1]))",fitmin,fitmax);
  fback = new TF1("fback","pol2",fitmin,fitmax);
  fall = new TF1("fall","gaus(0) +0.5*[0]*TMath::Exp(-[3]*TMath::Abs(x-[1])) + pol2(4)",fitmin,fitmax);

  fsign->SetLineColor(2);
  fback->SetLineColor(4);

  if(kSigma2vs3){
    kGoodMatch=kFALSE;
    kOverAll = 0;
  }

  if(bayesVsigma){
    kOverAll = 0;
    kGoodMatch=kFALSE;
    kSigma2vs3=kFALSE;
    kTOFmatch=kTRUE;
    weightS = -0.7;
  }

  return 1;
}

void doeffPi(Int_t pos,Float_t prob,Float_t etaminkp,Float_t etamaxkp){
  LoadLib();
  TH1D *hm = hmatched->ProjectionX("matchingPiEff",cmin,cmax);
  TH1D *ht = htracked->ProjectionX("tracking",cmin,cmax);

  hm->GetYaxis()->SetTitle("TOF matching eff.");
  hm->SetTitle("Using probability as weights");

  hm->Sumw2();
  ht->Sumw2();

  hm->Divide(hm,ht,1,1,"B");

  Int_t nptbin = binPid[2];
  Float_t minptbin = xmin[2];
  Float_t maxptbin = xmax[2];

  if(pos == 0){
    nptbin = binPid[3];
    minptbin = xmin[3];
    maxptbin = xmax[3];
  }

  if(prob > 0.1999|| kPid3Sigma ||kPid2Sigma){
    kGoodMatch = kFALSE;
    kSigma2vs3 = kFALSE;
//    if(! kOverAll) require5sigma = kTRUE;
    if(!isMC) weightS = -0.95;
  }

  TCanvas *c1 = new TCanvas();
  c1->Divide((nptbin+1)/2,2);
  TH2F *hh,*hh2;
  TH1D *h;
  char name[100];
  Float_t b[50][3];

  Double_t xx[50],yy[50];
  Double_t exx[50],eyy[50];

  for(Int_t i=0;i < nptbin;i++){
    c1->cd(i+1);//->SetLogy();
    Float_t ptmin = minptbin+(maxptbin-minptbin)/nptbin*(i);
    Float_t ptmax = minptbin+(maxptbin-minptbin)/nptbin*(i+1);

    xx[i] = (ptmin+ptmax)/2;
    exx[i] = (-ptmin+ptmax)/2;

    Float_t pp=0.1;
    if(prob < 0.2) pp = 0.;
    if(pos) hh=GetHistoPip(ptmin,ptmax,pp,0.0,etaminkp,etamaxkp);
    else hh=GetHistoPin(ptmin,ptmax,pp,0.0);
    sprintf(name,"TOF matched: %f < p_{T} < %f GeV/#it{c}",ptmin,ptmax);
    hh->SetTitle(name);
    sprintf(name,"hNoPid%i",i);
    
    pp=prob;
    if(prob < 0.2) pp = 0.1;
    if(pos) hh2=GetHistoPip(ptmin,ptmax,pp,0.0,etaminkp,etamaxkp);
    else hh2=GetHistoPin(ptmin,ptmax,pp,0.0);
    AddHisto(hh,hh2,weightS);

    h = hh->ProjectionX(name,cmin,cmax);
    h->RebinX(rebinsize);
    h->Draw("ERR");
    h->SetMarkerStyle(24);
    b[i][0]=-1;
    Int_t ntrial = 0;
    Float_t chi2 = 10000;
    while(ntrial < 3 && (chi2 > 20 + 1000*selectTrue)){
      fit(h,b[i],"WW","",xx[i]);
      c1->Update();
//       getchar();
      fit(h,b[i],"","",xx[i]);
      ntrial++;
      chi2 = b[i][2];
      printf("chi2 = %f\n",chi2);
      c1->Update();
//       getchar();
      
    }

    yy[i] = fall->GetParameter(parplotted);
    eyy[i] = fall->GetParError(parplotted);
  }

  TGraphErrors *gpar = new TGraphErrors(nptbin,xx,yy,exx,eyy);
  c1->cd(8);
//   gpar->Draw("AP");
  gpar->SetMarkerStyle(20);

  TCanvas *c2 = new TCanvas();
  c2->Divide((nptbin+1)/2,2);
  Float_t b2[50][3];

  for(Int_t i=0;i < nptbin;i++){
    c2->cd(i+1);
    Float_t ptmin = minptbin+(maxptbin-minptbin)/nptbin*(i);
    Float_t ptmax = minptbin+(maxptbin-minptbin)/nptbin*(i+1);

    Float_t pp=prob;
    if(prob < 0.2) pp = 0.1;
    if(pos) hh=GetHistoPip(ptmin,ptmax,pp,0.0,etaminkp,etamaxkp);
    else hh=GetHistoPin(ptmin,ptmax,pp,0.0);
    sprintf(name,"P_{TOF} > 0.8: %f < p_{T} < %f GeV/#it{c}",ptmin,ptmax);
    hh->SetTitle(name);
    sprintf(name,"hPid60_%i",i);
    h = hh->ProjectionX(name,cmin,cmax);
    h->RebinX(rebinsize);
    h->Draw("ERR");
    h->SetMarkerStyle(24);
    b2[i][0]=-1;
    Int_t ntrial = 0;
    Float_t chi2 = 10000;
    while(ntrial < 3 && (chi2 > 20 + 1000*selectTrue)){
      fit(h,b2[i],"WW","",xx[i]);
      fit(h,b2[i],"","",xx[i]);
      ntrial++;
      chi2 = b2[i][2];
      printf("chi2 = %f\n",chi2);
    }
    yy[i] = fall->GetParameter(parplotted);
    eyy[i] = fall->GetParError(parplotted);

  }

  TGraphErrors *gpar2 = new TGraphErrors(nptbin,xx,yy,exx,eyy);
  c2->cd(8);
//   gpar2->Draw("AP");
  gpar2->SetMarkerStyle(20);
  
  Double_t xpt[50],expt[50],eff[50],efferr[50];
  for(Int_t i=0;i<nptbin;i++){
    printf("%f +/- %f -  %f +/- %f\n",b[i][0],b[i][1],b2[i][0],b2[i][1]);

    Float_t ptmin = minptbin+(maxptbin-minptbin)/nptbin*(i);
    Float_t ptmax = minptbin+(maxptbin-minptbin)/nptbin*(i+1);

    xpt[i] = (ptmin+ptmax)/2;
    expt[i] = (-ptmin+ptmax)/2;
    eff[i] = b2[i][0]/(b[i][0]-b2[i][0]*weightS);

    //    b[i][0] = b[i][0]-b2[i][0]*weightS;

    //    efferr[i] = TMath::Sqrt(b[i][1]*b[i][1]/b[i][0]/b[i][0] + b2[i][1]*b2[i][1]/b2[i][0]/b2[i][0])*(b2[i][0]+b2[i][1])*(1+weightS*(b2[i][0]-b2[i][1])/b[i][0])/b[i][0];//*(1-eff[i]);//der2*der2*(b[i][1]*b[i][1] - b2[i][1]*b2[i][1]));

    efferr[i] = 1./(b[i][0]-b2[i][0]*weightS)/(b[i][0]-b2[i][0]*weightS)*TMath::Sqrt(b[i][0]*b[i][0]*b2[i][1]*b2[i][1] + b2[i][0]*b2[i][0]*b[i][1]*b[i][1]);

    if(TMath::Abs(efferr[i]) > 1)efferr[i]=1;
  }
  new TCanvas();
  TGraphErrors *geff = new TGraphErrors(nptbin,xpt,eff,expt,efferr);
  geff->Draw("AP");

  char flag[100];
  flag[0] = '\0';

  if(isMC){
    if(selectTrue) sprintf(flag,"true");
    else if(!keepTrue) sprintf(flag,"back");
  }

  Bool_t kWriteME = kFALSE;

  char flag2[100];
  flag2[0] = '\0';

  char etarange[100];
  sprintf(etarange,"_%.1f-%.1f_",etaminkp,etamaxkp);

  if(kGoodMatch)
    sprintf(flag2,"GM");

  if(bayesVsigma)
    sprintf(flag2,"BayesVsSigma");

  if(kSigma2vs3)
    sprintf(flag2,"Sigma2vs3");

  if(kOverAll)
    sprintf(flag2,"OverAll");
  if(kOverAllTOFmatch)
    sprintf(flag2,"OverAllTOF");
  if(kOverAll2Sigma)
    sprintf(flag2,"OverAll2sigma");

  if(kPid3Sigma)
    sprintf(flag2,"pid3sigma");
  if(kPid2Sigma)
    sprintf(flag2,"pid2sigma");


  if(pos){
    if(prob >=0.2) sprintf(name,"pionPos%sP%iEff%i_%i%s%s.root",etarange,Int_t(prob*100),(cmin-1)*10,cmax*10,flag,flag2);
    else{
      sprintf(name,"pionPos%sMatchEff%i_%i%s%s.root",etarange,(cmin-1)*10,cmax*10,flag,flag2);
      if(!(kOverAll || bayesVsigma || kGoodMatch || kSigma2vs3)) kWriteME = kTRUE;
    }
  }
  else{
    if(prob >=0.2) sprintf(name,"pionNeg%sP%iEff%i_%i%s%s.root",etarange,Int_t(prob*100),(cmin-1)*10,cmax*10,flag,flag2);
    else sprintf(name,"pionNeg%sMatchEff%i_%i%s%s.root",etarange,(cmin-1)*10,cmax*10,flag,flag2);
  }

  geff->SetTitle("#pi efficiency (from K^{0}_{s});p_{T} (GeV/#it{c};efficiency");
  TFile *fout = new TFile(name,"RECREATE");
  geff->Write();
  if(kWriteME) hm->Write();
  fout->Close();

  if(kWriteME) hm->Draw("SAME");
}

TH2F *GetHistoPip(Float_t pt,Float_t ptM,Float_t pMinkp,Float_t pMinkn,Float_t etaminkp,Float_t etamaxkp){

  Float_t x[] = {xmin[0]+0.001,etaminkp+0.001,pt+0.001,xmin[3]+0.001,pMinkp+0.001,pMinkn+0.001,(pMinkp>0.09 || kPid3Sigma||kPid2Sigma)+0.001,kTOFmatch+0.001,selectTrue,xmin[9],xmin[10],xmin[11],xmin[12],xmin[13]};
  Float_t x2[] = {xmax[0],etamaxkp-0.001,ptM-0.001,xmax[3],xmax[4],xmax[5],xmax[6],xmax[7],keepTrue,xmax[9],xmax[10],xmax[11],xmax[12],xmax[13]};

  if(kOverAll){
    x[6] = 0.0001;
    x2[9] = 5.9;
    if(pMinkp > 0.19) x2[9] = 4.9;
  }

  if(kOverAllTOFmatch && pMinkp > 0.19){
    x[6] = 1.0001;
    x2[9] = 4.9;
  }
  
  if(kOverAll2Sigma && pMinkp > 0.09){
    x2[9] = 2.;
    x[6] = 1.0001;
  }

  if(kGoodMatch){
    x[6] = 1.0001;
    if(pMinkp > 0)
      x2[9] = 4.9;
      
  }

  if(kTOFmatch){
    x[6] = 1.0001;
  }

  if(kSigma2vs3){
    x[6] = 1.0001;
    x2[9] = 3;
    if(pMinkp > 0)
      x2[9] = 2;
  }

  if(bayesVsigma){
    if(pMinkp > 0){
      x[4] = 0.2001;
      x2[9] = 5;
    }
    else{
      x2[9] = 3;
    }

    
  }

  if(require5sigma) x2[9] = 4.9;
  if(kPid3Sigma && pMinkp>0.09) x2[9] = 2.9;
  if(kPid2Sigma && pMinkp>0.09) x2[9] = 1.9;

  printf("max sigma = %f\n",x2[9]);

  AliPIDperfContainer *tmp = (AliPIDperfContainer *) fContPid1;

  TH2F *h = tmp->GetQA(0, x, x2);

  h->GetXaxis()->SetTitle("M_{K^{0}_{s}} (GeV/#it{c}^{2})");
  h->GetYaxis()->SetTitle("centrality [%]");

  return h;
}

TH2F *GetHistoPin(Float_t pt,Float_t ptM,Float_t pMinkn,Float_t pMinkp,Float_t etaminkp,Float_t etamaxkp){

  Float_t x[] = {xmin[0]+0.001,etaminkp+0.001,xmin[2]+0.001,pt+0.001,pMinkp+0.001,pMinkn+0.001,kTOFmatch+0.001,(pMinkn>0.09 || kPid3Sigma|| kPid2Sigma)+0.001,selectTrue,xmin[9],xmin[10],xmin[11],xmin[12],xmin[13]};
  Float_t x2[] = {xmax[0],etamaxkp-0.001,xmax[2],ptM-0.001,xmax[4],xmax[5],xmax[6],xmax[7],keepTrue,xmax[9],xmax[10],xmax[11],xmax[12],xmax[13]};

  if(kOverAll){
    x[7] = 0.0001;
    x2[10] = 5.9;
    if(pMinkn > 0.19) x2[10] = 4.9;
  }

  if(kOverAllTOFmatch && pMinkn > 0.19){
    x[7] = 1.0001;
    x2[10] = 4.9;
  }

  if(kOverAll2Sigma && pMinkn > 0.09){
    x2[10] = 2;
    x[7] = 1.0001;
  }

  if(kGoodMatch){
    x[7] = 1.0001;
    if(pMinkn > 0)
      x2[10] = 4.9;
      
  }

  if(kTOFmatch){
    x[7] = 1.0001;
  }

  if(kSigma2vs3){
    x[7] = 1.0001;
    x2[10] = 3;
    if(pMinkn > 0)
      x2[10] = 2;
  }
 
 if(bayesVsigma){
    if(pMinkn > 0){
      x[5] = 0.2001;
      x2[10] = 5;
    }
    else{
      x2[10] = 3;
    }    
  }

  if(require5sigma) x2[10] = 4.9;

  if(kPid3Sigma && pMinkn>0.09) x2[10] = 2.9;
  if(kPid2Sigma && pMinkn>0.09) x2[10] = 1.9;

  printf("max sigma = %f\n",x2[10]);

  AliPIDperfContainer *tmp = (AliPIDperfContainer *) fContPid2;

  TH2F *h = tmp->GetQA(0, x, x2);

  h->GetXaxis()->SetTitle("M_{K^{0}_{s}} (GeV/#it{c}^{2})");
  h->GetYaxis()->SetTitle("centrality [%]");

  return h;
}

void fit(TH1D *h,Float_t *a,char *opt,char *opt2,Float_t pt){
  if(h->Integral(1,h->GetNbinsX()) < 1){
    if(a){
      a[0]=0.001;
      a[1]=1;
    }
    return;
  }
  

 fall->SetParameter(0,100);
 fall->SetParameter(1,0.4971);
 fall->SetParameter(2,2.89748e-03);
 fall->FixParameter(3,230+30/pt);

 fall->SetParLimits(0,0.00001,1000000);
 fall->SetParLimits(1,0.4965,0.4985);
 fall->SetParLimits(2,0.0025,0.005);
 //fall->SetParLimits(3,200,350);

 fall->ReleaseParameter(4);
 fall->ReleaseParameter(5);

 if(selectTrue){
   fall->FixParameter(4,0);
   fall->FixParameter(5,0);
   fall->FixParameter(6,0);
 }
   fall->FixParameter(6,0);


 char namenew[100];
 sprintf(namenew,"%s_%i",h->GetName(),Int_t(gRandom->Rndm()*10000));
 TH1D *h2 = new TH1D(*h);
 h2->SetName(namenew);

// Float_t entries = h2->GetBinContent(h2->FindBin(0.497));
//  printf("entries under the peak = %f, pt = %f\n",entries,pt);
//  getchar();

 if(pt > 2.5){
   if(pt < 2.8) h2->RebinX(2);
   else if(pt < 3) h2->RebinX(4);
   else h2->RebinX(4);
 }

 h=h2;

 char name[100];
 TF1 *ftmp=fall;

 TF1 *ftmp2=new TF1(*fsign);
 sprintf(name,"fsign%i",ifunc);
 ftmp2->SetName(name);

 TF1 *ftmp3=new TF1(*fback);
 sprintf(name,"ftmp3%i",ifunc);
 ftmp3->SetName(name);

 ifunc++;

 h->Fit(ftmp,opt,opt2,fitmin,fitmax);
 h->Draw("ERR");

 ftmp2->SetParameter(0,ftmp->GetParameter(0));
 ftmp2->SetParameter(1,ftmp->GetParameter(1));
 ftmp2->SetParameter(2,ftmp->GetParameter(2));
 ftmp2->SetParameter(3,ftmp->GetParameter(3));
 ftmp2->Draw("SAME");
 ftmp3->SetParameter(0,ftmp->GetParameter(4));
 ftmp3->SetParameter(1,ftmp->GetParameter(5));
 ftmp3->SetParameter(2,ftmp->GetParameter(6));
 ftmp3->Draw("SAME");

 Float_t mean = ftmp->GetParameter(1);
 Float_t sigma = TMath::Abs(ftmp->GetParameter(2));

 Float_t signI = ftmp2->Integral(mean-10*sigma,mean+10*sigma)/h->GetBinWidth(1);
 Float_t backI = ftmp3->Integral(mean-3*sigma,mean+3*sigma)/h->GetBinWidth(1);

 if(signI < 0) signI = 0;
 if(backI < 1) backI = 1;

 Float_t errI = TMath::Sqrt(ftmp->GetParError(0)*ftmp->GetParError(0)/(0.001+ftmp->GetParameter(0))/(0.001+ftmp->GetParameter(0)));

 printf("signal(5 sigma) = %f +/- %f(fit) +/- %f(stat)\n",signI,errI*signI,TMath::Sqrt(signI));
 printf("backgr(3sigma) = %f\n",backI);
 printf("significance(3 sigma) = %f\n",signI/sqrt(signI+backI));

 if(a){
   a[0]=signI;
   a[1]=signI*errI*signI*errI + signI;
   a[1] = TMath::Sqrt(a[1]);
   if(ftmp->GetNDF()) a[2] = ftmp->GetChisquare()/ftmp->GetNDF();


   if(selectTrue){
     a[0] = h->Integral(1,h->GetNbinsX());
     a[1] = TMath::Sqrt(a[0]);
   }
 }
}

void AddHisto(TH2F *h1,TH2F *h2,Float_t w){
  Int_t nbinx = h1->GetNbinsX();
  Int_t nbiny = h1->GetNbinsY();

  for(Int_t i=1;i<=nbinx;i++){
    for(Int_t j=1;j<=nbiny;j++){
      Double_t val = h1->GetBinContent(i,j) + h2->GetBinContent(i,j)*w;
      Float_t err = TMath::Min(TMath::Sqrt(val),val);
      h1->SetBinContent(i,j,val);
      h1->SetBinError(i,j,err);
    }
  }
}
