void RP(){

  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);

  TFile *f = new TFile("../flow11h.root");

  if(!f->IsOpen()){
    printf("No such file\n");
    return 1;
  }

  char key[55] ;
//Get histograms 

  //Read resolution
  TH2F * hcos2AC = (TH2F*)f->Get("cos2AC");
  TH2F * hcos2V0AC = (TH2F*)f->Get("cos2V0AC");
  TH2F * hcos2V0ATPC = (TH2F*)f->Get("cos2V0ATPC");
  TH2F * hcos2V0CTPC = (TH2F*)f->Get("cos2V0CTPC");

//  TH1F * hresCh = (TH1F*)f2->Get("hEvPlResEta4");  

  TH1D * hres = 0 ;
  TH1D * hresF = 0 ;

  TH1D* rpT2  = new TH1D("rpTPC", "rpTPC 2 sub", hcos2AC->GetNbinsY(),0,100.);
  TH1D* rpT3 = new TH1D("rpTPC2", "rpTPC 3 sub", hcos2AC->GetNbinsY(),0,100.);
  TH1D* rpA  = new TH1D("rpV0A", "rpV0A", hcos2AC->GetNbinsY(),0,100.);
  TH1D* rpC  = new TH1D("rpV0C", "rpV0C", hcos2AC->GetNbinsY(),0,100.);

  Double_t resT2, resTA, resTC, resAC;

  for(int i=0;i<hcos2AC->GetNbinsY();i++){
    hres = hcos2AC->ProjectionX("res",i,i) ;
    resT2 = GetRes(hres);
    if(resT2>0)rpT2->SetBinContent(i, resT2);

    hresF = hcos2V0AC->ProjectionX("res2",i,i) ;
    resAC = GetCos(hresF);
    hres = hcos2V0ATPC->ProjectionX("res3",i,i) ;
    resTA = GetCos(hres);
    hres = hcos2V0CTPC->ProjectionX("res4",i,i) ;
    resTC = GetCos(hres);

    if(resAC>0)rpT3->SetBinContent(i, TMath::Sqrt(resTA*resTC/resAC));
    if(resTC>0)rpA->SetBinContent(i, TMath::Sqrt(resTA*resAC/resTC));
    if(resTA>0)rpC->SetBinContent(i, TMath::Sqrt(resTC*resAC/resTA));

//cout<<"i="<<i<<", rpA="<<TMath::Sqrt(resTA*resAC/resTC)<< endl;
  }

  rpT2->SetTitle("Reaction plane resolution with different ALICE detectors");
  rpT2->GetXaxis()->SetTitle("centrality percentage");
  rpT2->GetYaxis()->SetTitle("resolution");

  rpT2->SetMarkerStyle(20);
  rpT2->SetMarkerColor(2);
  rpT2->SetLineColor(2);
  rpT2->Draw("p");

  rpT3->SetMarkerStyle(21);
  rpT3->SetMarkerColor(2);
  rpT3->SetLineColor(2);
  rpT3->Draw("psame");

  rpA->SetMarkerStyle(22);
  rpA->SetMarkerColor(kMagenta);
  rpA->SetLineColor(kMagenta);
  rpA->Draw("psame");

  rpC->SetMarkerStyle(22);
  rpC->SetMarkerColor(kGreen);
  rpC->SetLineColor(kGreen);
  rpC->Draw("psame");

/*
  hresCh->SetMarkerStyle(21);
  hresCh->SetMarkerColor(kBlue);
  hresCh->SetLineColor(kBlue);
  hresCh->Draw("same");
*/

  TLegend * l = new TLegend(0.12,0.15,0.6,0.35) ;
  l->SetFillColor(0) ;
  l->SetTextSize(0.03) ;
  l->AddEntry(rpT2,"RPRes with TPC 2 subs","p") ;
  l->AddEntry(rpT3,"RPRes with TPC 3 subs","p") ;
  l->AddEntry(rpA,"RPRes with V0A","p") ;
  l->AddEntry(rpC,"RPRes with V0C","p") ;
//  l->AddEntry(hresCh,"RPRes TPC (Alex Dobrin)","p") ;
  l->Draw() ;

}

Double_t Ollitrault(Double_t chi){
  Double_t x = 0.25*chi*chi;
  Double_t resk1 = 0.626657 * chi * exp(-x) * (TMath::BesselI0((float)x) + TMath::BesselI1((float)x));
  Double_t resk2 = 0.626657 * chi * exp(-x) * (TMath::Sqrt(2./TMath::Pi()/x)*TMath::SinH(x) + TMath::Sqrt(2./TMath::Pi()/x)*(TMath::CosH(x) - TMath::SinH(x)/x));

  return resk1;
}

Double_t GetCos(TH1D* hres){
  Double_t resMean = hres->GetMean() ;
  cout<<resMean<<endl;
  if(resMean>0)return resMean ;
  else return 0. ;
}

Double_t GetRes(TH1D* hres){
  Double_t resMean = hres->GetMean() ;
  Double_t resOld=0 ;
  if(resMean>0)resOld=1./TMath::Sqrt(resMean) ;

  float rxn2=resMean;

  double chi = 2;
  for(int j=0; j<15; j++){
    double resSub = sqrt(rxn2);
    Double_t resEve = Ollitrault(chi);
    chi= (resEve < resSub) ? chi+1.0*pow(0.5, j) : chi-1.0*pow(0.5, j);
  }
  return Ollitrault(TMath::Sqrt(2)*chi);
}

