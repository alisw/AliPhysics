// E. Sicking, CERN
// macro for analysing output of AliAnalysisTaskHardSoft

void plotSignalBackground(){

  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);  

  Float_t piHalf = TMath::Pi()/2.;

  TFile *file0 = TFile::Open(Form("Clean.root"));
  TList *l1 = (TList *)file0->Get("HardSoftHists");


  //histograms with entries how often a trigger particle 
  //is found at given Nch
  TH1F *trigger = (TH1F*)l1->FindObject("fTriggerESD");
  TCanvas * c0= new TCanvas("c0", "c0", 100, 100, 620, 480);
  c0->cd();
  trigger->Draw();
  //  for(Int_t i=0;i<100;i++){
  //     cout << i << " : "<< trigger->GetBinContent(i+1) << endl;
  //   }
  

  //delta phi histograms for each Nch bin
  TH1D *dp[100];
  TF1 *f[100]; // fit function for background region around pi/2
  Float_t fAtPiHalf[100]; // value  f(pi/2)
  
  for(Int_t i=0;i<100;i++){
    dp[i]= (TH1D*)l1->FindObject(Form("fDPhiLeadingNchBinESD%02d",i));
    //normalize to number of triggers
    if(dp[i]->GetEntries()>0)dp[i]->Scale(1./trigger->GetBinContent(i+1)); 
   
    //define fit function
    f[i]= new TF1(Form("fun%2d",i), "pol2", TMath::Pi()/2-1, TMath::Pi()/2+1);
    f[i]->SetLineColor(kRed);
    fAtPiHalf[i]=0;
  
  }

  TGraphErrors *tge = new TGraphErrors(); // of side band = background
  tge->SetMarkerStyle(4);
  tge->SetMarkerColor(kRed);
  tge->SetLineColor(kRed);

  TGraphErrors *tgeSignal = new TGraphErrors(); //integral of dphi<0.7
  tgeSignal->SetMarkerStyle(4);
  tgeSignal->SetMarkerColor(kBlue);
  tgeSignal->SetLineColor(kBlue);

  TGraphErrors *tgeSignalWOB = new TGraphErrors();//integral of dphi<0.7 
                                                  //- background
  tgeSignalWOB->SetMarkerStyle(4);
  tgeSignalWOB->SetMarkerColor(kGreen);
  tgeSignalWOB->SetLineColor(kGreen);

  

  TCanvas * c= new TCanvas("c", "c", 150, 150, 820, 620);
  c->Divide(7,7);
  for(Int_t i=0;i<49;i++){
    c->cd(i+1);
    if(dp[i]->GetEntries()>0) {
      dp[i]->Fit(f[i],"R");
      dp[i]->SetMinimum(0);
    }

    // get value of fit function at pi/2 -> side band
    fAtPiHalf[i]=f[i]->Eval(piHalf,0,0);
    
    Double_t integral07=dp[i]->Integral(1,40);// bin 40 -> dphi=07
   
 

    Float_t radius=0.7;
    Float_t DeltaPhiBinWith = TMath::Pi()/180.;
    Float_t bkgSampleSize   = DeltaPhiBinWith*1.8;
    Float_t halfCircle      = TMath::Pi()*(radius)*(radius)/2;
    Float_t akzeptanz       =0.8;
    Float_t ratio= akzeptanz*halfCircle/bkgSampleSize;


    tge->SetPoint(i,i,fAtPiHalf[i] ); //side band
    tgeSignal->SetPoint(i,i, integral07); // integral
    tgeSignalWOB->SetPoint(i,i, integral07-fAtPiHalf[i]*40);//integral - bkg  
    
    
  }


  //define legend
  TLegend *legp;
  TF1 *fun1p;
  TF1 *fun2p;
  TF1 *fun3p;
  TF1 *fun4p;
  TF1 *fun5p;
  fun1p= new TF1("fun1p","gaus",-5.0,5.0);
  fun2p= new TF1("fun2p","gaus",-5.0,5.0);
  fun3p= new TF1("fun3p","gaus",-5.0,5.0);
  fun4p= new TF1("fun4p","gaus",-5.0,5.0);
  fun5p= new TF1("fun5p","gaus",-5.0,5.0);
  fun1p->SetMarkerColor(kBlue);
  fun1p->SetMarkerStyle(4);
  fun2p->SetMarkerColor(kGreen);
  fun2p->SetMarkerStyle(4);
  fun3p->SetMarkerColor(kMagenta);
  fun3p->SetMarkerStyle(4);
  fun4p->SetMarkerColor(kRed);
  fun4p->SetMarkerStyle(4);
  legp= new TLegend(0.1,0.7,0.75,0.9);
  legp->SetFillColor(kWhite);
  legp->AddEntry(fun1p,"integral dphi<0.7","p");   
  legp->AddEntry(fun4p,"background in side band","p");   
  legp->AddEntry(fun2p,"integral dphi<0.7 - background","p");   
  legp->AddEntry(fun5p,"Ntracks around trigger in R<0.7","l");   
  legp->AddEntry(fun3p,"Ntracks around trigger in R<0.7 - background","p");   


  //define histogram with name of axes
  TCanvas * c1= new TCanvas("c1", "c1", 200, 200, 620, 480);
  c1->cd();
  TH1F * histo = new TH1F("histo", "",100,-0.5,99.5 );
  histo->SetMaximum(10);
  histo->SetMinimum(0.01);
  histo->GetXaxis()->SetRangeUser(0,50);
  histo->SetXTitle("N_{charge}");
  histo->SetTitle("");
  histo->SetYTitle("");
  histo->Draw();
  legp->Draw();

  //draw output
  tgeSignal->Draw("p");
  tge->Draw("p");
  tgeSignalWOB->Draw("p");

  
  //number of associate particles within  R<0.7 around trigger
  TProfile *pSignal= (TProfile*)l1->FindObject("fNchAssInRESD");
  pSignal->Draw("same");

  TGraphErrors *t = new TGraphErrors();
  t->SetMarkerStyle(4);
  t->SetMarkerColor(kMagenta);
  t->SetLineColor(kMagenta);
  t->SetLineStyle(2);
  for(Int_t i=0;i<49;i++){
    t->SetPoint(i,i,pSignal->GetBinContent(i+1)-fAtPiHalf[i]*ratio);
  }

  t->Draw("p");

  //   TCanvas * c5= new TCanvas("c5", "c5", 300, 300, 620, 480);
  //   c5->cd();
  //   t->Draw("ap");
  
  //   TCanvas * c6= new TCanvas("c6", "c6", 400, 400, 620, 480);
  //   c6->cd();
  //   tge->Draw("ap");

}
