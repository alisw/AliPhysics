char inname[256];

//  const Int_t nbin=11 ;
//  Double_t xa[14]={1.,1.5,2.,2.5,3.,4.,5.,6.,8.,10.,14.,20.} ;
  const Int_t nbin=12 ;
  Double_t xa[13]={0.6,1.,1.5,2.,3.,4.,5.,7.,9.,11.,13.,16.,20.} ;
//  const Int_t nbin=11 ;
//  Double_t xa[12]={0.6,1.,1.5,2.,2.5,3.,4.,5.,7.,10.,15.,20.} ;
    

Double_t Ollitrault(Double_t chi){
  Double_t x = 0.25*chi*chi;
  Double_t resk1 = 0.626657 * chi * exp(-x) * (TMath::BesselI0((float)x) + TMath::BesselI1((float)x));
  Double_t resk2 = 0.626657 * chi * exp(-x) * (TMath::Sqrt(2./TMath::Pi()/x)*TMath::SinH(x) + TMath::Sqrt(2./TMath::Pi()/x)*(TMath::CosH(x) - TMath::SinH(x)/x));

  return resk1;
}

Double_t GetCos(TH2F* hcos2AC, Int_t cen){
  TH1D *hres = 0 ;
  TH2F * hPHOS = (TH2F*)gROOT->FindObjectAny("hCenPHOS");
  Double_t resMean = 0, weight=0. ;

  if(cen==21){
  TFile * f3 = new TFile("../flow11h_Apr14.root") ;
  hPHOS = (TH2F*)f3->Get("hCenPHOS");
  }

Int_t bin0, bin1;

if(cen==10){ bin0=1; bin1=3;}
if(cen==11){ bin0=5; bin1=9;}
if(cen==0){ bin0=1; bin1=2;}
if(cen==1){ bin0=2; bin1=3;}
if(cen==2){ bin0=3; bin1=5;}
if(cen==3){ bin0=5; bin1=7;}
if(cen==4){ bin0=7; bin1=9;}
if(cen==5){ bin0=9; bin1=11;}
if(cen==20){ bin0=3; bin1=11;} //10-50%
if(cen==21){ bin0=1; bin1=5;} //0-20%

  for(Int_t i=bin0; i<bin1; i++){
    TH1D * hwi = hPHOS->ProjectionX("wi",5*i-4,5*i) ; 
    Double_t wi = hwi->GetMean() ;
    hres = hcos2AC->ProjectionX("res",i,i) ; 
    resMean += wi*hres->GetMean() ;
    weight += wi ;
    delete hwi ;
    delete hres ;
  }
  if(weight>0.){
    resMean/=weight ;
    return resMean ;
  }
  else{
    return 0. ;
  }
}
Double_t GetRes(TH2F* hcos2AC, Int_t cen){

  TH1D *hres = 0 ;
  Double_t resMean = 0, weight=0. ;

Int_t bin0, bin1;

if(cen==10){ bin0=1; bin1=3;} //0-10 (Central)
if(cen==11){ bin0=5; bin1=9;} //20-40 (SemiCentral)
if(cen==0){ bin0=1; bin1=2;}
if(cen==1){ bin0=2; bin1=3;}
if(cen==2){ bin0=3; bin1=5;}
if(cen==3){ bin0=5; bin1=7;}
if(cen==4){ bin0=7; bin1=9;}
if(cen==5){ bin0=9; bin1=11;}
if(cen==20){ bin0=3; bin1=11;} //10-50%
if(cen==21){ bin0=1; bin1=5;} //0-20%

  TH1D *hres = 0 ;
  TH2F * hPHOS = (TH2F*)gROOT->FindObjectAny("hCenPHOS");
  Double_t resMean = 0, weight=0. ;
  if(cen==21){
  TFile * f3 = new TFile("../flow11h_Apr14.root") ;
  hPHOS = (TH2F*)f3->Get("hCenPHOS");
  }

  for(Int_t i=bin0; i<bin1; i++){
    TH1D * hwi = hPHOS->ProjectionX("wi",5*i-4,5*i) ; 
    Double_t wi = hwi->GetMean() ;
    hres = hcos2AC->ProjectionX("res",i,i) ; 
    resMean += wi*hres->GetMean() ;
    weight += wi ;
    delete hwi ;
    delete hres ;
  }
  if(weight>0.){
    resMean/=weight ;
  }
  else{
    resMean = 0. ;
  }

    
  double chi = 2;
  double resSub = TMath::Sqrt(resMean);
  for(int j=0; j<15; j++){
    double resEve = Ollitrault(chi);
    chi= (resEve < resSub) ? chi+1.0*pow(0.5, j) : chi-1.0*pow(0.5, j);
  }
  return Ollitrault(TMath::Sqrt(2.)*chi);
}

TH3F* GetRealMixed(Int_t cen, char w[255], char d[255], char kind[255], char kind2[255], char PID[255], Int_t Real=1){
  if(cen==10||cen==0||cen==1||cen==21) sprintf(inname,"../flow11hCentral_Apr14.root");
  else sprintf(inname,"../flow11hSemiCentral_Apr14.root");
  char inname2[255];
  if(cen==21) sprintf(inname2,"../flow11hSemiCentral_Apr14.root");

  TFile * f = new TFile(inname) ;
  if(cen==21) TFile * f2 = new TFile(inname2) ;

  char key[125] ;
  TH3F * h3DR =0 ;TH3F * h3DM=0;

  if(cen==10){//centrality 0-10 
    sprintf(key,"h%sMassPt%s%s%s_cen0",w,d,kind,PID) ;
    h3DR = (TH3F*)f->Get(key) ;
    sprintf(key,"hMiMassPt%s%s%s_cen0",d,kind2,PID) ;
    h3DM= (TH3F*)f->Get(key) ;

    sprintf(key,"h%sMassPt%s%s%s_cen1",w,d,kind,PID) ;
    TH3F *h3DRa = (TH3F*)f->Get(key) ;
    sprintf(key,"hMiMassPt%s%s%s_cen1",d,kind2,PID) ;
    TH3F *h3DMa= (TH3F*)f->Get(key) ;

    sprintf(key,"h%sMassPt%s%s%s_cen2",w,d,kind,PID) ;
    h3DRb = (TH3F*)f->Get(key) ;
    sprintf(key,"hMiMassPt%s%s%s_cen2",d,kind2,PID) ;
    h3DMb= (TH3F*)f->Get(key) ;
    sprintf(key,"h%sMassPt%s%s%s_cen3",w,d,kind,PID) ;
    h3DRc = (TH3F*)f->Get(key) ;
    sprintf(key,"hMiMassPt%s%s%s_cen3",d,kind2,PID) ;
    h3DMc= (TH3F*)f->Get(key) ;

    TH1D* hCen = ((TH2F*)f->Get("hCentrality"))->ProjectionX();
    Double_t e1 = hCen->Integral(6,8)/3.;
    Double_t e2 = hCen->Integral(9,9);
    Double_t e3 = hCen->Integral(10,10);

cout<<"w2="<<e1/e2<<", w3="<<e1/e3<<endl;

    h3DR->Add(h3DRa); delete h3DRa;
    h3DM->Add(h3DMa); delete h3DMa;
    h3DR->Add(h3DRb,e1/e2); delete h3DRb;
    h3DM->Add(h3DMb,e1/e2); delete h3DMb;
    h3DR->Add(h3DRc,e1/e3); delete h3DRc;
    h3DM->Add(h3DMc,e1/e3); delete h3DMc;
  }
  if(cen==21){//centrality 0-20 
    sprintf(key,"h%sMassPt%s%s%s_cen0",w,d,kind,PID) ;
    h3DR = (TH3F*)f->Get(key) ;
    sprintf(key,"hMiMassPt%s%s%s_cen0",d,kind2,PID) ;
    h3DM= (TH3F*)f->Get(key) ;

    sprintf(key,"h%sMassPt%s%s%s_cen1",w,d,kind,PID) ;
    TH3F *h3DRa = (TH3F*)f->Get(key) ;
    sprintf(key,"hMiMassPt%s%s%s_cen1",d,kind2,PID) ;
    TH3F *h3DMa= (TH3F*)f->Get(key) ;

    sprintf(key,"h%sMassPt%s%s%s_cen2",w,d,kind,PID) ;
    h3DRb = (TH3F*)f->Get(key) ;
    sprintf(key,"hMiMassPt%s%s%s_cen2",d,kind2,PID) ;
    h3DMb= (TH3F*)f->Get(key) ;
    sprintf(key,"h%sMassPt%s%s%s_cen3",w,d,kind,PID) ;
    h3DRc = (TH3F*)f->Get(key) ;
    sprintf(key,"hMiMassPt%s%s%s_cen3",d,kind2,PID) ;
    h3DMc= (TH3F*)f->Get(key) ;

    TH1D* hCen = ((TH2F*)f->Get("hCentrality"))->ProjectionX();
    Double_t e1 = hCen->Integral(6,8)/3.;
    Double_t e2 = hCen->Integral(9,9);
    Double_t e3 = hCen->Integral(10,10);

cout<<"w2="<<e1/e2<<", w3="<<e1/e3<<endl;

    h3DR->Add(h3DRa); delete h3DRa;
    h3DM->Add(h3DMa); delete h3DMa;
    h3DR->Add(h3DRb,e1/e2); delete h3DRb;
    h3DM->Add(h3DMb,e1/e2); delete h3DMb;
    h3DR->Add(h3DRc,e1/e3); delete h3DRc;
    h3DM->Add(h3DMc,e1/e3); delete h3DMc;

//now 10-20
    sprintf(key,"h%sMassPt%s%s%s_cen4",w,d,kind,PID) ;
    TH3F* h3DR1 = (TH3F*)f2->Get(key) ;
    sprintf(key,"hMiMassPt%s%s%s_cen4",d,kind2,PID) ;
    TH3F* h3DM1= (TH3F*)f2->Get(key) ;

    sprintf(key,"h%sMassPt%s%s%s_cen0",w,d,kind,PID) ;
    TH3F* h3DRa = (TH3F*)f2->Get(key) ;
    sprintf(key,"hMiMassPt%s%s%s_cen0",d,kind2,PID) ;
    TH3F* h3DMa= (TH3F*)f2->Get(key) ;
    sprintf(key,"h%sMassPt%s%s%s_cen1",w,d,kind,PID) ;
    TH3F* h3DRb = (TH3F*)f2->Get(key) ;
    sprintf(key,"hMiMassPt%s%s%s_cen1",d,kind2,PID) ;
    TH3F* h3DMb= (TH3F*)f2->Get(key) ;
    sprintf(key,"h%sMassPt%s%s%s_cen2",w,d,kind,PID) ;
    TH3F* h3DRc = (TH3F*)f2->Get(key) ;
    sprintf(key,"hMiMassPt%s%s%s_cen2",d,kind2,PID) ;
    TH3F* h3DMc= (TH3F*)f2->Get(key) ;
    sprintf(key,"h%sMassPt%s%s%s_cen3",w,d,kind,PID) ;
    TH3F* h3DRd = (TH3F*)f2->Get(key) ;
    sprintf(key,"hMiMassPt%s%s%s_cen3",d,kind2,PID) ;
    TH3F* h3DMd= (TH3F*)f2->Get(key) ;

    TH1D* hCen2 = ((TH2F*)f2->Get("hCentrality"))->ProjectionX();
    Double_t e1 = hCen2->Integral(11,11);
    Double_t e2 = hCen2->Integral(12,12);
    Double_t e3 = hCen2->Integral(13,13);
    Double_t e4 = hCen2->Integral(14,15)/2.;
    Double_t e5 = hCen2->Integral(16,20)/5.;

cout<<"w_11="<<e5/e1<<", w_12="<<e5/e2<<", w_13="<<e5/e3<<", w_14-15="<<e5/e4<<endl;

    h3DR->Add(h3DR1,e5/e1); delete h3DR1;
    h3DM->Add(h3DM1,e5/e1); delete h3DM1;
    h3DR->Add(h3DRa,e5/e1); delete h3DRa;
    h3DM->Add(h3DMa,e5/e1); delete h3DMa;
    h3DR->Add(h3DRb,e5/e2); delete h3DRb;
    h3DM->Add(h3DMb,e5/e2); delete h3DMb;
    h3DR->Add(h3DRc,e5/e3); delete h3DRc;
    h3DM->Add(h3DMc,e5/e3); delete h3DMc;
    h3DR->Add(h3DRd,e5/e4); delete h3DRd;
  }
  if(cen==11){//centrality 20-40 
    sprintf(key,"h%sMassPt%s%s%s_cen5",w,d,kind,PID) ;
    h3DR = (TH3F*)f->Get(key) ;
    sprintf(key,"hMiMassPt%s%s%s_cen5",d,kind2,PID) ;
    h3DM= (TH3F*)f->Get(key) ;
    
    sprintf(key,"h%sMassPt%s%s%s_cen6",w,d,kind,PID) ;
    TH3F *h3DRa = (TH3F*)f->Get(key) ;
    sprintf(key,"hMiMassPt%s%s%s_cen6",d,kind2,PID) ;
    TH3F *h3DMa= (TH3F*)f->Get(key) ;

    h3DR->Add(h3DRa); delete h3DRa;
    h3DM->Add(h3DMa); delete h3DMa;
  }
  if(cen==20){//centrality 10-50
    sprintf(key,"h%sMassPt%s%s%s_cen4",w,d,kind,PID) ;
    h3DR = (TH3F*)f->Get(key) ;
    sprintf(key,"hMiMassPt%s%s%s_cen4",d,kind2,PID) ;
    h3DM= (TH3F*)f->Get(key) ;

    sprintf(key,"h%sMassPt%s%s%s_cen5",w,d,kind,PID) ;
    TH3F *h3DRa = (TH3F*)f->Get(key) ;
    sprintf(key,"hMiMassPt%s%s%s_cen5",d,kind2,PID) ;
    TH3F *h3DMa= (TH3F*)f->Get(key) ;

    sprintf(key,"h%sMassPt%s%s%s_cen6",w,d,kind,PID) ;
    TH3F *h3DRb = (TH3F*)f->Get(key) ;
    sprintf(key,"hMiMassPt%s%s%s_cen6",d,kind2,PID) ;
    TH3F *h3DMb= (TH3F*)f->Get(key) ;

    sprintf(key,"h%sMassPt%s%s%s_cen7",w,d,kind,PID) ;
    TH3F *h3DRc = (TH3F*)f->Get(key) ;
    sprintf(key,"hMiMassPt%s%s%s_cen7",d,kind2,PID) ;
    TH3F *h3DMc= (TH3F*)f->Get(key) ;

    h3DR->Add(h3DRa); delete h3DRa;
    h3DM->Add(h3DMa); delete h3DMa;
    h3DR->Add(h3DRb); delete h3DRb;
    h3DM->Add(h3DMb); delete h3DMb;
    h3DR->Add(h3DRc); delete h3DRc;
    h3DM->Add(h3DMc); delete h3DMc;
  }
//small bins
  if(cen==0){//centrality 0-5 
    sprintf(key,"h%sMassPt%s%s%s_cen0",w,d,kind,PID) ;
    h3DR = (TH3F*)f->Get(key) ;
    sprintf(key,"hMiMassPt%s%s%s_cen0",d,kind2,PID) ;
    h3DM= (TH3F*)f->Get(key) ;
  }
  if(cen==1){//centrality 5-10   
    sprintf(key,"h%sMassPt%s%s%s_cen1",w,d,kind,PID) ;
    h3DR = (TH3F*)f->Get(key) ;
    sprintf(key,"hMiMassPt%s%s%s_cen1",d,kind2,PID) ;
    h3DM= (TH3F*)f->Get(key) ;
    sprintf(key,"h%sMassPt%s%s%s_cen2",w,d,kind,PID) ;
    h3DRa = (TH3F*)f->Get(key) ;
    sprintf(key,"hMiMassPt%s%s%s_cen2",d,kind2,PID) ;
    h3DMa= (TH3F*)f->Get(key) ;
    sprintf(key,"h%sMassPt%s%s%s_cen3",w,d,kind,PID) ;
    h3DRb = (TH3F*)f->Get(key) ;
    sprintf(key,"hMiMassPt%s%s%s_cen3",d,kind2,PID) ;
    h3DMb= (TH3F*)f->Get(key) ;

    TH1D* hCen = ((TH2F*)f->Get("hCentrality"))->ProjectionX();
    Double_t e1 = hCen->Integral(6,8)/3.;
    Double_t e2 = hCen->Integral(9,9);
    Double_t e3 = hCen->Integral(10,10);

cout<<"w2="<<e1/e2<<", w3="<<e1/e3<<endl;

    h3DR->Add(h3DRa,e1/e2); delete h3DRa;
    h3DM->Add(h3DMa,e1/e2); delete h3DMa;
    h3DR->Add(h3DRb,e1/e3); delete h3DRb;
    h3DM->Add(h3DMb,e1/e3); delete h3DMb;

/*
    h3DR->Add(h3DRa); delete h3DRa;
    h3DM->Add(h3DMa); delete h3DMa;
    h3DR->Add(h3DRb); delete h3DRb;
    h3DM->Add(h3DMb); delete h3DMb;
*/
  }
  if(cen==2){//centrality 10-20   
    sprintf(key,"h%sMassPt%s%s%s_cen4",w,d,kind,PID) ;
    h3DR = (TH3F*)f->Get(key) ;
    sprintf(key,"hMiMassPt%s%s%s_cen4",d,kind2,PID) ;
    h3DM= (TH3F*)f->Get(key) ;

    sprintf(key,"h%sMassPt%s%s%s_cen0",w,d,kind,PID) ;
    TH3F* h3DRa = (TH3F*)f->Get(key) ;
    sprintf(key,"hMiMassPt%s%s%s_cen0",d,kind2,PID) ;
    TH3F* h3DMa= (TH3F*)f->Get(key) ;
    sprintf(key,"h%sMassPt%s%s%s_cen1",w,d,kind,PID) ;
    TH3F* h3DRb = (TH3F*)f->Get(key) ;
    sprintf(key,"hMiMassPt%s%s%s_cen1",d,kind2,PID) ;
    TH3F* h3DMb= (TH3F*)f->Get(key) ;
    sprintf(key,"h%sMassPt%s%s%s_cen2",w,d,kind,PID) ;
    TH3F* h3DRc = (TH3F*)f->Get(key) ;
    sprintf(key,"hMiMassPt%s%s%s_cen2",d,kind2,PID) ;
    TH3F* h3DMc= (TH3F*)f->Get(key) ;
    sprintf(key,"h%sMassPt%s%s%s_cen3",w,d,kind,PID) ;
    TH3F* h3DRd = (TH3F*)f->Get(key) ;
    sprintf(key,"hMiMassPt%s%s%s_cen3",d,kind2,PID) ;
    TH3F* h3DMd= (TH3F*)f->Get(key) ;

    TH1D* hCen = ((TH2F*)f->Get("hCentrality"))->ProjectionX();
    Double_t e1 = hCen->Integral(11,11);
    Double_t e2 = hCen->Integral(12,12);
    Double_t e3 = hCen->Integral(13,13);
    Double_t e4 = hCen->Integral(14,15)/2.;
    Double_t e5 = hCen->Integral(16,20)/5.;

cout<<"w_11="<<e5/e1<<", w_12="<<e5/e2<<", w_13="<<e5/e3<<", w_14-15="<<e5/e4<<endl;

    h3DR->Add(h3DRa,e5/e1); delete h3DRa;
    h3DM->Add(h3DMa,e5/e1); delete h3DMa;
    h3DR->Add(h3DRb,e5/e2); delete h3DRb;
    h3DM->Add(h3DMb,e5/e2); delete h3DMb;
    h3DR->Add(h3DRc,e5/e3); delete h3DRc;
    h3DM->Add(h3DMc,e5/e3); delete h3DMc;
    h3DR->Add(h3DRd,e5/e4); delete h3DRd;
    h3DM->Add(h3DMd,e5/e4); delete h3DMd;
  }
  if(cen==3){//centrality 20-30   
    sprintf(key,"h%sMassPt%s%s%s_cen5",w,d,kind,PID) ;
    h3DR = (TH3F*)f->Get(key) ;
    sprintf(key,"hMiMassPt%s%s%s_cen5",d,kind2,PID) ;
    h3DM= (TH3F*)f->Get(key) ;
  }
  if(cen==4){//centrality 30-40   
    sprintf(key,"h%sMassPt%s%s%s_cen6",w,d,kind,PID) ;
    h3DR = (TH3F*)f->Get(key) ;
    sprintf(key,"hMiMassPt%s%s%s_cen6",d,kind2,PID) ;
    h3DM= (TH3F*)f->Get(key) ;
  }
  if(cen==5){//centrality 40-50   
    sprintf(key,"h%sMassPt%s%s%s_cen7",w,d,kind,PID) ;
    h3DR = (TH3F*)f->Get(key) ;
    sprintf(key,"hMiMassPt%s%s%s_cen7",d,kind2,PID) ;
    h3DM= (TH3F*)f->Get(key) ;
  }


  if(Real)return h3DR;
  else return h3DM;
}

TH1F * GetCen(){
cout<<inname<<endl;
  TFile * f = new TFile(inname) ;
  return (TH1F*)f->Get("hCentrality") ;
}

TH2F * GetCenPHOS(){
  TFile * f = new TFile(inname) ;
  return (TH2F*)f->Get("hCenPHOS") ;
}


