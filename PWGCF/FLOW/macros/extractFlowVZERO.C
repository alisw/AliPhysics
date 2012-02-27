extractFlowVZERO(Int_t itech=0,Int_t ibin = 3,Int_t step = 0,Int_t arm=2,Float_t pTh = 0.9,Int_t charge=0,Int_t addbins=0){
  TCanvas *can1 = new TCanvas("cV0Acheck","cV0Acheck");

/*
read outVZEROv2.root or outVZEROv3.root files

itech=0 TPC stand alone PID and selection
itech=1 TPC stand alone PID but TOF track selection
itech=2 TPC&TOF PID (TOF strictly required)
itech=3 TPC|TOF PID (the same of before but using TPC stand alone when TOF not available)

ibin = centrality bin

step = species (0=all charges, 1=pi, 2=K, 3=p, 4=e, 5=d, 6=t, 7=3He

arm = armonic 2 (v2) or 3 (v3)

pTh = probability threshold = 0.6 , 0.8 , 0.9 , 0.95

charge = 1(pos), -1(neg), 0(both)

addbins = to merge more centrality bin (i.e. ibin = 2, addbins = 3 merge 2,3,4,5 centrality bins) - preliminary version!!!!!!

*/

  // load libraries
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libCORRFW.so");

  // get objects
  char filein[100];
  sprintf(filein,"outVZEROv%i.root",arm);
  TFile *f=new TFile(filein);
  sprintf(filein,"contVZEROv%i",arm);
  TList *l = f->Get(filein);
  l->ls();
  AliCFContainer *c = l->At(0);
  AliCFContainer *c3 = l->At(1);
  TProfile *p1 = l->At(2);
  TProfile *p2 = l->At(3);
  TProfile *p3 = l->At(4);
  Int_t iVar[2] = {1,3};

  Int_t iVarPsi[1] = {5};
  /* 0=centrality , 1= cos(arm*deltaphi) , 2=charge , 3=pt , 4=probability , 5=EP , 6=PIDmask*/
  Double_t iMin[7] = {ibin,-1,-1.5,0,pTh+0.00001,-3.14,0};
  Double_t iMax[7] = {ibin+addbins,1,1.5,20,1.1,3.14,1};

  char *iTechName[4] = {"tpc","tpcInTof","tof","tpcORtof"};

  Bool_t kOnlyTPC = kFALSE;

  if(itech==0 && (step!=0)){
      kOnlyTPC = kTRUE;
  }
  else if(itech==1 && (step!=0)){
      iMin[6] = 2;
      iMax[6] = 2;
  }
  else if(itech==2 || itech==1){
      iMin[6] = 1;
      iMax[6] = 1;
  }
  else{
      iMin[6] = 0;
      iMax[6] = 1;
  }

  if(charge==-1) iMax[2]=0;
  if(charge==1) iMin[2]=0;

  // EP distribution needed for NUA corrections
  TH2F *hPhiA = l->At(5);
  TH2F *hPhiC = l->At(6);

  TH1D *hPhiPA = hPhiA->ProjectionY("_py",iMin[0]+1,iMax[0]+1);
  TH1D *hPhiPC = hPhiC->ProjectionY("_py",iMin[0]+1,iMax[0]+1);

  // EP resolution variables
  Float_t res1=0,res2=0,res3=0; 
  Float_t eres1=0,eres2=0,eres3=0; 

  for(Int_t i=iMin[0];i<=iMax[0];i++){ // in case of more centrality bins weighting with the errors
    if(p1->GetBinError(i+1)){
      eres1 += 1./p1->GetBinError(i+1)/p1->GetBinError(i+1);
      res1 += p1->GetBinContent(i+1)/p1->GetBinError(i+1)/p1->GetBinError(i+1);      
    }
    if(p2->GetBinError(i+1)){
      eres2 += 1./p2->GetBinError(i+1)/p2->GetBinError(i+1);
      res2 += p2->GetBinContent(i+1)/p2->GetBinError(i+1)/p2->GetBinError(i+1);      
    }
    if(p3->GetBinError(i+1)){
      eres3 += 1./p3->GetBinError(i+1)/p3->GetBinError(i+1);
      res3 += p3->GetBinContent(i+1)/p3->GetBinError(i+1)/p3->GetBinError(i+1);      
    }

    if(eres1) res1 /= eres1;
    if(eres2) res2 /= eres2;
    if(eres3) res3 /= eres3;

    if(eres1) eres1 = TMath::Sqrt(1./eres1);
    if(eres1) eres2 = TMath::Sqrt(1./eres2);
    if(eres1) eres3 = TMath::Sqrt(1./eres3);
  }

  // NUA correction (fit to EP distribution)
  TF1 *fpol = new TF1("fPol","pol0");
  hPhiPA->Fit("fPol","","",-TMath::Pi()/arm,TMath::Pi()/arm);
  Float_t scalA = fPol->GetParameter(0);
  hPhiPC->Fit("fPol","","",-TMath::Pi()/arm,TMath::Pi()/arm);
  Float_t scalC = fPol->GetParameter(0);

  AliCFContainer *c2[20];
  AliCFContainer *c4[20];
  TH2D *h[20],*h2[20];

  AliCFContainer *c2bis[20];
  AliCFContainer *c4bis[20];

  AliCFContainer *cPsi2[20];
  AliCFContainer *cPsi4[20];
  TH1D *hPsi[20],*hPsi2[20];
  AliCFContainer *cPsi2bis[20];
  AliCFContainer *cPsi4bis[20];

  Float_t intA = 0;
  Float_t intC = 0;

  for(Int_t i=5;i<15;i++){
    printf("%i\n",i);
    iMin[5] = hPhiPA->GetBinCenter(i+1);
    iMax[5] = hPhiPA->GetBinCenter(i+1);
    if(!kOnlyTPC){
      c2[i] = c->MakeSlice(2,iVar,iMin,iMax);
      c4[i]= c3->MakeSlice(2,iVar,iMin,iMax);
    }
    else{ // merge to maskPID bins (TPC standalone without TOF + TPC stand alone with TOF)
      iMin[6] = 0;
      iMax[6] = 0;
      c2[i] = c->MakeSlice(2,iVar,iMin,iMax);
      c4[i]= c3->MakeSlice(2,iVar,iMin,iMax);

      iMin[6] = 2;
      iMax[6] = 2;
      c2bis[i] = c->MakeSlice(2,iVar,iMin,iMax);
      c4bis[i]= c3->MakeSlice(2,iVar,iMin,iMax);
   }


    if(!kOnlyTPC){
      cPsi2[i] = c->MakeSlice(1,iVarPsi,iMin,iMax);
      cPsi4[i]= c3->MakeSlice(1,iVarPsi,iMin,iMax);
    }
    else{ // merge to maskPID bins (TPC standalone without TOF + TPC stand alone with TOF)
      iMin[6] = 0;
      iMax[6] = 0;
      cPsi2[i] = c->MakeSlice(1,iVarPsi,iMin,iMax);
      cPsi4[i]= c3->MakeSlice(1,iVarPsi,iMin,iMax);

      iMin[6] = 2;
      iMax[6] = 2;
      cPsi2bis[i] = c->MakeSlice(1,iVarPsi,iMin,iMax);
      cPsi4bis[i]= c3->MakeSlice(1,iVarPsi,iMin,iMax);
    }

    h[i] = (TH2D *) c2[i]->Project(step,0,1);
    h2[i] = (TH2D *) c4[i]->Project(step,0,1);
    if(kOnlyTPC){
      TH2D *htemp = (TH2D *) c2bis[i]->Project(step,0,1);
      h[i]->Add(htemp);
      htemp = (TH2D *) c4bis[i]->Project(step,0,1);
      h2[i]->Add(htemp);
    }

    if(hPhiPA->GetBinContent(i)) h[i]->Scale(scalA/hPhiPA->GetBinContent(i+1)); // reweighting for NUA correction
    if(hPhiPC->GetBinContent(i)) h2[i]->Scale(scalC/hPhiPC->GetBinContent(i+1)); // reweighting for NUA correction

    hPsi[i] = (TH1D *) cPsi2[i]->Project(step,0);
    hPsi2[i] = (TH1D *) cPsi4[i]->Project(step,0);

    if(kOnlyTPC){ // merge to maskPID bins (TPC standalone without TOF + TPC stand alone with TOF)
      TH1D *htemp2 = (TH1D *) cPsi2bis[i]->Project(step,0);
      hPsi[i]->Add(htemp2);
      htemp2 = (TH1D *) cPsi4bis[i]->Project(step,0);
      hPsi2[i]->Add(htemp2);
    }

    if(hPhiPA->GetBinContent(i+1)) hPsi[i]->Scale(scalA/hPhiPA->GetBinContent(i+1)); // check of reweighting for NUA correction
    if(hPhiPC->GetBinContent(i+1)) hPsi2[i]->Scale(scalC/hPhiPC->GetBinContent(i+1)); // check of reweighting for NUA correction

    intA += hPsi[i]->Integral();
    intC += hPsi2[i]->Integral();

  }

//   hPhiPA->Scale(intA / (scalA)/10);
//   hPhiPC->Scale(intC / (scalC)/10);

  // NUA correction check V0A
  hPhiPA->Draw();
  hPsi[5]->Draw("SAME");
  hPsi[5]->Scale(scalA / intA * 10);
  for(Int_t i=6;i<15;i++){
    h[5]->Add(h[i]);
    h2[5]->Add(h2[i]);
    hPsi[i]->Draw("SAME");
    hPsi[i]->Scale(scalA / intA * 10);
  }

  // NUA correction check V0C
  TCanvas *can2 = new TCanvas("cV0Ccheck","cV0Ccheck");
  hPhiPC->Draw();
  hPsi2[5]->Draw("SAME");
  for(Int_t i=6;i<15;i++){
    h[5]->Add(h[i]);
    h2[5]->Add(h2[i]);
    hPsi2[i]->Draw("SAME");
    hPsi2[i]->Scale(scalC /intC * 10);
  }
  hPsi2[5]->Scale(scalC /intC * 10);

  // Flow for V0A and V0C separately
  TCanvas *can3 = new TCanvas("cFlowComp","cFlowComp");
  TProfile *pp = h[5]->ProfileY();
  pp->Draw();
  TProfile *pp2 = h2[5]->ProfileY();
  pp2->Draw("SAME");

  printf("nev (selected within 0-80%) = %i\n",p1->GetEntries());


  // correction for resoltion
  Float_t scaling = sqrt(res1*res3/res2);
  pp->Scale(1./scaling);
  printf("resolution V0A = %f\n",scaling);
  Float_t err1_2 = eres1*eres1/res1/res1/4 +
                   eres2*eres2/res2/res2/4 +
                   eres3*eres3/res3/res3/4;
  Float_t err2_2 = err1_2;
  err1_2 /= scaling*scaling;
  scaling = sqrt(res2*res3/res1);
  err2_2 /= scaling*scaling;
  pp2->Scale(1./scaling);
  printf("resolution V0C =%f\n",scaling);

  char title[100];
  sprintf(title,"VZERO EP;p_{t} (GeV/c);v_{%i}",arm);
  pp->SetTitle(title);

  // Average V0A-V0C
  TH1D *pAll = pp->ProjectionX();

  for(Int_t i=1;i <= pAll->GetNbinsX();i++){
       Float_t e1 = err1_2*pp->GetBinContent(i)*pp->GetBinContent(i) + pp->GetBinError(i)*pp->GetBinError(i);
       Float_t e2 = err2_2*pp2->GetBinContent(i)*pp2->GetBinContent(i) + pp2->GetBinError(i)*pp2->GetBinError(i);
       Float_t xval = 0,exval = 0;
       if(e1 >0 && e2>0){
         xval = (pp->GetBinContent(i)/e1 + pp2->GetBinContent(i)/e2)/(1/e1 + 1/e2);     
         exval = 1./sqrt(1/e1 + 1/e2);
       }
       pAll->SetBinContent(i,xval);
       pAll->SetBinError(i,exval);
  }
  // combined measurement
  TCanvas *can4 = new TCanvas("cFlowCombined","cFlowCombined");
  pAll->Draw();

  char name[100];
  sprintf(name,"out%i-%i_%i_%4.2f_%i%sv%i.root",iMin[0],iMax[0],step,pTh,charge,iTechName[itech],arm);
  TFile *fout = new TFile(name,"RECREATE");

  pAll->SetName("histo");
  pAll->Write();
  can1->Write();
  can2->Write();
  can3->Write();
  fout->Close();
}
