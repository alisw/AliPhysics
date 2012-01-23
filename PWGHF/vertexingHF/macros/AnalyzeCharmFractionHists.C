TCanvas *cCompareSimple,*cCompareCum;
TCanvas *cSubtraction,*cCompareSubtractSimple,*cCompareSubtractCum;
const Int_t timeWait=500;

void AnalyzeCharmFractionHists(Int_t ptbin,Int_t rebin){

  //
  // Macro for analyzing the impact parameter distribution of the D0 candidates.
  //  
  // andrea.rossi@ts.infn.it
  //
  //==========================================================================


  DrawHistos("d0D0_PureBack.root","d0D0NoMCSel_SideBand.root","hd0D0",ptbin,rebin);
  SubtractHist("d0D0_Signal.root","d0D0NoMCSel.root","d0D0NoMCSel_SideBand.root","hd0D0",ptbin,rebin);
  return;

}

TH1F* CumulativeHist(const TH1F *h1,TString nameHist=0x0){
  //ASSUME SAME BINNING
  TString histname="hCumulative";
  if(!nameHist.IsNull())histname.Append(nameHist.Data());
  Int_t nbins;
  Double_t minX,maxX,binwidth,cumul=0.,errcumul=0.;
  nbins=h1->GetNbinsX();
  minX=h1->GetBinCenter(1);
  maxX=h1->GetBinCenter(nbins);
  binwidth=h1->GetBinWidth(1);

  TH1F *h1copy=new TH1F();
  *h1copy=*h1;
  //h1copy->Sumw2();

  TH1F *hCumul=new TH1F(histname.Data(),histname.Data(),nbins,minX-binwidth/2.,maxX+binwidth/2.);


  for(Int_t j=1;j<=nbins;j++){
    cumul+=h1copy->GetBinContent(j);
    hCumul->SetBinContent(j,cumul);
  }
  hCumul->Sumw2();
  hCumul->Scale(1./h1copy->Integral(1,nbins));

  delete h1copy;
  return hCumul;
}


///////////////////// HISTOGRAMS COMPARISON   //////////////////////////////////
//              
//    Compare the desired histograms. The default case is the comparison
//    between the d0 distribution for background (from MC) D0 candidates with the
//    d0 distribution for "side bands" d0 candidates
//
///////////////////////////////////////////////////////////////////////////////////
             
void DrawHistos(Int_t pt,Int_t rebin=1){
  TString hist="hd0D0";
  DrawHistos("","",hist,pt,rebin);
  return;
}
void DrawHistos(TString file1,TString file2,TString hist,Int_t pt,Int_t rebin=1){
  
  if(pt>=0){
    hist.Append("pt");
    hist+=pt;
  }
  TString nameH=hist;
  Double_t integr1,integr2;
  if(file1.IsNull())file1="d0D0_PureBack.root";
  if(file2.IsNull())file2="d0D0NoMCSel_SideBand.root";
  TFile *fSide=TFile::Open(file2.Data());
  TFile *fmer=TFile::Open(file1.Data());
  TH1F *h1=(TH1F*)fmer->Get(hist.Data());
  nameH.Append("_Gaus");
  h1->SetName(nameH.Data());
  nameH=hist;
  TH1F *h2=(TH1F*)fSide->Get(hist.Data());
  nameH.Append("_SideBands");
  h2->SetName(nameH.Data());

  Bool_t satisfied=kFALSE;
  
  while(!satisfied){
    h1->Rebin(rebin);
    h2->Rebin(rebin);
    TH1F *hCumul1=CumulativeHist(h1,"_Gaus");
    hCumul1->SetDrawOption("E4");
    TH1F *hCumul2=CumulativeHist(h2,"_SideBands");
    hCumul2->SetDrawOption("E4");
    
    
    TH1F *hDivCumul=new TH1F();
    *hDivCumul=*hCumul2;
    //    hDivCumul->Sumw2();
    hDivCumul->Divide(hCumul1,hCumul2);
    hDivCumul->SetName("RatioCumulGausSdBands");
    hDivCumul->SetTitle("Ratio Cumulative Gaus Over SideBands");
    
    TH1F *hGausComp=new TH1F();
    *hGausComp=*h1;
    nameH=h1->GetName();
    nameH.Append("Compare");
    hGausComp->SetName(nameH.Data());
    hGausComp->SetLineColor(kRed);
    hGausComp->Sumw2();
    integr1=hGausComp->Integral();
    
    TH1F *hSideComp=new TH1F();
    *hSideComp=*h2;
    nameH=h2->GetName();
    nameH.Append("Compare");
    hSideComp->SetName(nameH.Data());
    hSideComp->SetLineColor(kBlue);
    hSideComp->Sumw2();
    integr2=hSideComp->Integral();
    if(integr2>integr1)hSideComp->Scale(integr1/integr2);
    else hGausComp->Scale(integr2/integr1);
    
    TH1F *hDivision=new TH1F();
    *hDivision=*h1;
    hDivision->SetName("RatioGausSdBands");
    hDivision->SetTitle("Ratio Gaus Over Sd Bands");
    hDivision->Sumw2();
    hDivision->Divide(hGausComp,hSideComp);
    
    
    cCompareSimple=new TCanvas("cCompareSimple","Comparison of Under Gaus and Under Side Band background",700,700);
    cCompareSimple->Divide(1,2);
    cCompareSimple->cd(1);
    hSideComp->Draw();
    hGausComp->Draw("Sames");
    cCompareSimple->Update();
    TPaveStats *p=hGausComp->FindObject("stats");
    p->SetTextColor(kRed);
    p=(TPaveStats*)hSideComp->FindObject("stats");
    p->SetTextColor(kBlue);
    cCompareSimple->cd(2);
    hDivision->Draw();
    cCompareSimple->Update();
    
    cCompareCum=new TCanvas("cCompareCum","Comparison of cumulative function",700,700);
    
    cCompareCum->Divide(1,2);
    
    cCompareCum->cd(1);
    
    hCumul1->SetLineColor(kRed);
    hCumul2->SetLineColor(kBlue);
    hCumul1->SetFillColor(kRed);
    hCumul2->SetFillColor(kBlue);
    hCumul1->Draw();
    hCumul2->Draw("Sames");
    cCompareCum->cd(2);
    hDivCumul->Draw();
    cCompareCum->Update();
    p=(TPaveStats*)hCumul1->FindObject("stats");
    p->SetTextColor(kRed);
    p=(TPaveStats*)hCumul2->FindObject("stats");
    p->SetTextColor(kBlue);
    cCompareCum->Update();
    // cCompareCum->cd(2);
    //hCumul2->Draw();
    //  cCompareCum->Update();
  
    for(Int_t timewait=0;timewait<timeWait;timewait++){
      cCompareCum->Update();
      cCompareSimple->Update();    
    }
    cout<<"Are you satisfied?"<<endl;
    cin>>satisfied;
    
    if(!satisfied){
     
    
      delete hCumul1;
      delete hCumul2;
      delete hDivCumul;
      delete cCompareCum;
      delete hSideComp;
      delete hGausComp;
      delete cCompareSimple;
      delete hDivision;
   
      cout<<"Rebin"<<endl;
      cin>>rebin;
      if(rebin<0)return;
    }
  }
  
  return;
  
}

///////////////////////////////////////   BACKGROUND SUBTRACTION EXAMPLE       ///////////////////
//
//
//   Performs: Signal Imp.Par Distr - B* Normalized Background Imp.Par Distr.(from Side Bands) 
//             and compare it with the MC Signal Imp. Par distribution. 
//             The amount of Background B is taken as the true one from the difference between 
//             the integrals of the "Observed Signal" and the MC Signal d0 distributions.
//           
//////////////////////////////////////////////////////////////////////////////////////////////

void SubtractHist(Int_t pt,Int_t rebin=1){
  SubtractHist("d0D0_Signal.root","d0D0NoMCSel.root","d0D0NoMCSel_SideBand.root","hd0D0",pt,rebin);
  
  return;
}

void SubtractHist(TString fileSignal,TString fileNoMC,TString fileNoMCSB,TString hist,Int_t pt,Int_t rebin){
  
  if(pt>=0){
    hist.Append("pt");
    hist+=pt;
  }
  
  TString nameH=hist;
  Double_t integr1,integr2;
  /*  if(fileNoMC.IsNull())fileNoMC="d0D0merged.root";
      if(fileNoMCSB.IsNull())fileNoMCSB="d0D0SideBmerged.root";*/
  if(gSystem->AccessPathName(fileSignal.Data())){
    Printf("Wrong signal file! \n");
    return;
  }
  if(gSystem->AccessPathName(fileNoMC.Data())){
    Printf("Wrong d0distr under inv mass peak file! \n");
    return;
  }
  if(gSystem->AccessPathName(fileNoMCSB.Data())){
    Printf("Wrong  d0distr Side band file! \n");
    return;
  }

  TFile *fSide=TFile::Open(fileNoMCSB.Data());
  TFile *fNoMCSignal=TFile::Open(fileNoMC.Data());
  
  TH1F *h1=(TH1F*)fNoMCSignal->Get(hist.Data());
  nameH.Append("_SelSignal");
  h1->SetName(nameH.Data());
  nameH=hist;
  h1->Sumw2();


  TH1F *h2=(TH1F*)fSide->Get(hist.Data());
  nameH.Append("_SideBands");
  h2->SetName(nameH.Data());
  h2->Sumw2();


  Double_t integrGaus,integrSB,integrSign;
  integrGaus=h1->Integral();
  integrSB=h2->Integral();
  

  TFile *fSignal=TFile::Open(fileSignal.Data());
  TH1F *hSign=(TH1F*)fSignal->Get(hist.Data());
  nameH=hist;
  nameH.Append("_MCSignal");
  hSign->SetName(nameH.Data());
  hSign->Sumw2();

  Double_t integrGaus,integrSB,integrSign;
  integrGaus=h1->Integral();
  integrSB=h2->Integral();
  integrSign=hSign->Integral();


   Bool_t satisfied=kFALSE;
   
   while(!satisfied){
     
     h1->Rebin(rebin);
     h2->Rebin(rebin);
     hSign->Rebin(rebin);
     
     TH1F *hGausSubtract=new TH1F();
     *hGausSubtract=*h1;
     nameH=hist;
     hist.Append("_Subtracted");
     hGausSubtract->SetName(hist.Data());
     //hGausSubtract->Sumw2();
     hGausSubtract->Add(h1,h2,1.,-1./integrSB*(integrGaus-integrSign));
     /*     hGausSubtract->Draw();
	    return;*/
     cSubtraction=new TCanvas("cSubtraction","cSubtraction",600,700);
     cSubtraction->cd();
     hGausSubtract->SetLineColor(kRed);
     hGausSubtract->Draw("E0"); 
     hSign->SetLineColor(kBlue);
     hSign->Draw("sames");
     cSubtraction->Update();
     TPaveStats *p=hGausSubtract->FindObject("stats");
     p->SetTextColor(kRed);
     cSubtraction->Update();
     p=(TPaveStats*)hSign->FindObject("stats");
     p->SetTextColor(kBlue);
     cSubtraction->Update();
     /*TString strText="Side Band integral: ";
       strText+=integrSB;
       TLatex *lat=new TLatex(500.,500.,strText.Data());// *text=p->AddText(strText.Data());
       lat->Draw();
       TPaveText *pvInfo = new TPaveText(72.74276,79.18782,92.84497,95.43147,"brNDC");
       pvInfo->SetName("pvInfo");
       pvInfo->SetBorderSize(2);
       pvInfo->SetFillColor(19);
       pvInfo->SetTextAlign(12);
       TText *text=p->AddText(strText.Data());
       pvInfo->Draw();
       hGausSubtract->GetListOfFunctions()->Add(pvInfo);
       //     pvInfo->SetParent(hGausSubtract->GetListOfFunctions());
       
       //     strText="MCSignal integral: ";
       //strText+=integrSign;
       text=p->AddText(strText.Data());
       strText="Sel Signal integral: ";
       strText+=integrGaus;
       text=p->AddText(strText.Data());
       cSubtraction->Modified();*/
  
     
     TH1F *hCumulGausSubtract=CumulativeHist(hGausSubtract,"_BackSubtr");
     hCumulGausSubtract->SetDrawOption("E4");
     TH1F *hCumulSign=CumulativeHist(hSign,"_MCSignal");
     hCumulSign->SetDrawOption("E4");
     
     
     TH1F *hDivCumul=new TH1F();
     *hDivCumul=*hCumulSign;
     hDivCumul->Divide(hCumulGausSubtract,hCumulSign);
     hDivCumul->SetName("RatioCumulBackSubtr_MCSignal");
     hDivCumul->SetTitle("Ratio Cumulative BackSubtracted Over MCSignal");
    
     TH1F *hGausSubComp=new TH1F();
     *hGausSubComp=*hGausSubtract;
     nameH=hGausSubtract->GetName();
     nameH.Append("Compare");
     hGausSubComp->SetName(nameH.Data());
     hGausSubComp->SetLineColor(kRed);
     integr1=hGausSubComp->Integral();
     
     TH1F *hSignComp=new TH1F();
     *hSignComp=*hSign;
     nameH=hSign->GetName();
     nameH.Append("Compare");
     hSignComp->SetName(nameH.Data());
     hSignComp->SetLineColor(kBlue);
     integr2=hSignComp->Integral();
     if(integr2>integr1)hSignComp->Scale(integr1/integr2);
     else hGausSubComp->Scale(integr2/integr1);
     
     TH1F *hDivision=new TH1F();
     *hDivision=*hGausSubtract;
     hDivision->SetName("RatioBackSubtr_Signal");
     hDivision->SetTitle("Ratio BackSubtracted Over MCSignal");
     hDivision->Divide(hGausSubComp,hSignComp);
     
     
     cCompareSubtractSimple=new TCanvas("cCompareSubtractSimple","Comparison of BackSubtracted and MCSignal",700,700);
     cCompareSubtractSimple->Divide(1,2);
     cCompareSubtractSimple->cd(1);
     hSignComp->Draw();
     hGausSubComp->Draw("Sames");
     cCompareSubtractSimple->cd(2);
     hDivision->SetLineColor(1);
     hDivision->Draw();
     cCompareSubtractSimple->Update();
     p=(TPaveStats*)hGausSubComp->FindObject("stats");
     p->SetTextColor(kRed);
     p=(TPaveStats*)hSignComp->FindObject("stats");
     p->SetTextColor(kBlue);
     cCompareSubtractSimple->Update();

     
     cCompareSubtractCum=new TCanvas("cCompareSubtractCumulative","Comparison of cumulative functions",700,700);
    
    cCompareSubtractCum->Divide(1,2);
    
    cCompareSubtractCum->cd(1);
    
    hCumulGausSubtract->SetLineColor(kRed);
    hCumulSign->SetLineColor(kBlue);
    hCumulGausSubtract->SetFillColor(kRed);
    hCumulSign->SetFillColor(kBlue);
    hCumulGausSubtract->Draw();
    hCumulSign->Draw("Sames");
    cCompareSubtractCum->cd(2);
    hDivCumul->SetLineColor(1);
    hDivCumul->Draw();

    cCompareSubtractCum->Update();
    p=(TPaveStats*)hCumulGausSubtract->FindObject("stats");
    p->SetTextColor(kRed);
    p=(TPaveStats*)hCumulSign->FindObject("stats");
    p->SetTextColor(kBlue);
    cCompareSubtractCum->Update();

    
    for(Int_t timewait=0;timewait<timeWait;timewait++){
      cCompareSubtractCum->Update();
      cCompareSubtractSimple->Update();  
      cSubtraction->Update();  
    }

    cout<<"Are you satisfied?"<<endl;
    cin>>satisfied;
    
    if(!satisfied){
      cout<<"Rebin"<<endl;
      cin>>rebin;
      delete hCumulGausSubtract;
      delete hCumulSign;
      delete hDivCumul;
      delete cCompareSubtractCum;
      delete hSignComp;
      delete hGausSubComp;
      delete cCompareSubtractSimple;
      delete hDivision;
    }
  }
  
   printf("Side Band integral: %d \n",integrSB);
   printf("Selected Signal integral: %d \n",integrGaus);
   printf("MC Signal integral: %d \n",integrSign);
   printf(" -> Background = %d \n",integrGaus-integrSign);
   
  return;
  
}

//////////////////////////////// FIT DISTRIBUTION WITH GAUSSIAN + EXP TAILS    ///////////////////
//
//            Fit a istogram with a gaus(mean,sigma) + a*exp(mean2,lambda) function
//
/////////////////////////////////////////////////////////////////////////////////////////////////
void DoFit(TString filename,TString histo,TString gausSide,Int_t ptbin,Int_t rebin=-1){

  TString  histname=histo;
  TString fileout=histo;

  if(ptbin>=0){
    histname.Append("pt");
    histname+=ptbin;
    fileout.Append("pt");
    fileout+=ptbin;
  }
  else fileout.Append("ptAll");
  
  TFile *fIN=TFile::Open(filename.Data());
  
  TH1F *h=(TH1F*)fIN->Get(histname.Data());

  Int_t ok=0;  
  Int_t binL=1,binR=h->GetNbinsX();
  
  while(ok!=1){
    if(rebin>0)h->Rebin(rebin);

    Double_t integr=h->Integral(binL,binR);
    Double_t binwidth=h->GetBinWidth(10);
    Double_t minX,maxX;
    Int_t nbin=h->GetNbinsX();
 
    
    TCanvas *c1=new TCanvas("c1","c1",700,600);
    c1->cd();
    //    h->SetFillColor(kRed);
    h->Draw("");
    //    c1->SetLogy();
    c1->Update();
  
    TF1 *f=CreateFunc(integr,binwidth);
    h->Fit(f->GetName(),"RI","",h->GetBinCenter(binL)-binwidth/2.,h->GetBinCenter(binR)+binwidth/2.);
    c1->Update();
    
    cout<<"Is it ok?"<<endl;
    cin>>ok;
    if(ok==1){
 
      fileout.Append("_");
      fileout.Append(gausSide.Data());
      fileout.Append("Fit.root");
      
      TFile *fOUT=new TFile(fileout.Data(),"RECREATE");
      fOUT->cd();
      c1->Write();
      h->Write();
      fOUT->Close();
      delete c1;
      delete f;
    }
    else if(ok==0){
      cout<<"rebin value= ";
      cin>>rebin;
      cout<<endl;
      cout<<"Change Interval?"<<endl;
      cin>>ok;
      if(ok){
	cout<<"Lower value= ";
	cin>>minX;
	cout<<endl;
	cout<<"Upper value= ";
	cin>>maxX;
	cout<<endl;
      
	binL=h->FindBin(minX);
	binR=h->FindBin(maxX);
      }
      ok=0;
      delete f;
      delete c1;
    }
    else {
      ok=1;
      delete f;
      delete c1;
    }
  }

  return;
}


TF1* CreateFunc(Double_t integral,Double_t binw){
  TF1 *func=new TF1("func","[5]*[6]*((1.-[2])*1./2./[1]*TMath::Exp(-TMath::Abs((x-[0])/[1]))+[2]/TMath::Sqrt(2*TMath::Pi())/[4]*TMath::Exp(-1*(x-[3])*(x-[3])/2./[4]/[4]))");

  func->FixParameter(5,integral);
  func->SetParName(5,"histInt");
  func->FixParameter(6,binw);
  func->SetParName(6,"binW");
  func->SetParLimits(0,-100,100.);
  func->SetParName(0,"expMean");
  func->SetParLimits(1,0.001,1000);
  func->SetParName(1,"expDecay");
  func->SetParLimits(2,0.00001,1.);
  func->SetParName(2,"fracGaus");
  func->SetParLimits(3,-100,100.);
  func->SetParName(3,"gausMean");
  func->SetParLimits(4,5.,200.);
  func->SetParName(4,"gausSigma");

  func->SetParameter(1,50.);
  func->SetParameter(4,50.);
  
  return func;

}

////////////////////////////////// CHI2 TEST: ////////////////////////
//                               IN PROGRESS
//
///////////////////////////////////////////////////////////////////////////////
Double_t GetChi2(const TH1F *h1,const TH1F *h2,TH1F *hchi2,Int_t &ndof){//ASSUME SAME BINNING
  
  Int_t int1=h1->Integral();
  Int_t int2=h2->Integral();
  Int_t nm,nM;
  TH1F **h=new TH1F*[2];
  TH1F *hchisq=new TH1F("hChi2","Chi2 histo",1000,0.,10.);
  
  if(int2>int1){
    nM=int2;
    nm=int1;
    h[0]=h2;
    h[1]=h1;
  }
  else {
    nM=int1;
    nm=int2;
    h[0]=h1;
    h[1]=h2;
  }
  ndof=0;
  Int_t nbin=h[0]->GetNbinsX();//ASSUME SAME BINNING
  
  Double_t chi2=0.,chi2Sum=0.;
  Double_t f1,fTrue;
  for(Int_t j=1;j<=nbin;j++){//THE BINNING IS NOT TAKEN INTO ACCOUNT
    
    if(h[1]->GetBinContent(j)==0||h[0]->GetBinContent(j)==0)continue;
    fTrue=h[0]->GetBinContent(j)/nM;
    chi2=(h[1]->GetBinContent(j)-fTrue*nm)*(h[1]->GetBinContent(j)-fTrue*nm)/h[1]->GetBinContent(j);
    hchisq->Fill(chi2);
    chi2Sum+=chi2;
    ndof++;
  }
  
  *hchi2=*hchisq;
  delete hchisq;

  return chi2Sum;
  

}

void TestChi2(){


  TH1F *h1=new TH1F("h1","h1",50,-10.,10.);
  TH1F *h2=new TH1F("h2","h2",50,-10.,10.);
  
  h1->FillRandom("gaus",50000);
  h2->FillRandom("gaus",80000);
  TH1F *hchi2=new TH1F();
  Int_t ndof;
  Double_t chi2=GetChi2(h1,h2,hchi2,ndof);
  hchi2->Draw();
  cout<<"The chi2 is "<<chi2;
}