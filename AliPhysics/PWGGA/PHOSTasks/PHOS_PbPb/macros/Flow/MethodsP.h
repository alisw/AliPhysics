char inname[256];

//  const Int_t nbin=11 ;
//  Double_t xa[14]={1.,1.5,2.,2.5,3.,4.,5.,6.,8.,10.,14.,20.} ;
//  const Int_t nbin=12 ;
//  Double_t xa[13]={1.,1.5,2.,2.5,3.,4.,5.,7.,9.,11.,13.,16.,20.} ;

//const Int_t nbin=16 ;
//Double_t xa[17]={0.6,1,1.5,2,2.5,3,3.5,4,5,6,7,8,10,12,14,16,20};

//from PCM
const Int_t nbin=29 ;
Double_t xa[30]={0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,3.0,3.5,4.0,5,6,7,8,10,12,14,16,20};

//  const Int_t nbin=12 ;
//  Double_t xa[13]={0.6,1.,1.5,2.,2.5,3.,4.,5.,7.,10.,15.,20.} ;


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

Int_t bin0, bin1;

if(cen==10){ bin0=1; bin1=3;}
if(cen==11){ bin0=5; bin1=9;}
if(cen==0){ bin0=1; bin1=2;}
if(cen==1){ bin0=2; bin1=3;}
if(cen==2){ bin0=3; bin1=5;}
if(cen==3){ bin0=5; bin1=7;}
if(cen==4){ bin0=7; bin1=9;}
if(cen==5){ bin0=9; bin1=11;}

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

  TH1D *hres = 0 ;
  TH2F * hPHOS = (TH2F*)gROOT->FindObjectAny("hCenPHOS");
  Double_t resMean = 0, weight=0. ;

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

TH2F* GetRealMixed(Int_t cen, char w[255], char d[255], char kind[255], char kind2[255], char PID[255], Int_t Real=1){
  if(cen==10||cen==0||cen==1) sprintf(inname,"../flow11hCentral_Apr14.root");
  else sprintf(inname,"../flow11hSemiCentral_Apr14.root");
  TFile * f = new TFile(inname) ;
  char key[125] ;
  TH2F * h2DR =0 ;TH2F * h2DM=0;

  if(cen==10){//centrality 0-10 
    sprintf(key,"h%sPhotPhi%s%s%s_cen0",w,d,kind,PID) ;
//cout<<key<<endl;
    h2DR = (TH2F*)f->Get(key) ;

    sprintf(key,"h%sPhotPhi%s%s%s_cen1",w,d,kind,PID) ;
    TH2F *h2DRa = (TH2F*)f->Get(key) ;

    sprintf(key,"h%sPhotPhi%s%s%s_cen2",w,d,kind,PID) ;
    h2DRb = (TH2F*)f->Get(key) ;
    sprintf(key,"h%sPhotPhi%s%s%s_cen3",w,d,kind,PID) ;
    h2DRc = (TH2F*)f->Get(key) ;

    TH1D* hCen = ((TH2F*)f->Get("hCentrality"))->ProjectionX();
    Double_t e1 = hCen->Integral(6,8)/3.;
    Double_t e2 = hCen->Integral(9,9);
    Double_t e3 = hCen->Integral(10,10);

cout<<"w2="<<e1/e2<<", w3="<<e1/e3<<endl;

//    h2DR->Add(h2DRa); delete h2DRa;
    h2DRa->Add(h2DRb,e1/e2); delete h2DRb;
    h2DRa->Add(h2DRc,e1/e3); delete h2DRc;
    h2DR->Add(h2DRa); delete h2DRa;


//    h2DR->Add(h2DRa); delete h2DRa;
//    h2DR->Add(h2DRb); delete h2DRb;
//    h2DR->Add(h2DRc); delete h2DRc;

  }
  if(cen==11){//centrality 20-40 
    sprintf(key,"h%sPhotPhi%s%s%s_cen5",w,d,kind,PID) ;
    h2DR = (TH2F*)f->Get(key) ;
    
    sprintf(key,"h%sPhotPhi%s%s%s_cen6",w,d,kind,PID) ;
    TH2F *h2DRa = (TH2F*)f->Get(key) ;

    h2DR->Add(h2DRa); delete h2DRa;
  }
  if(cen==0){//centrality 0-5 
    sprintf(key,"h%sPhotPhi%s%s%s_cen0",w,d,kind,PID) ;
    h2DR = (TH2F*)f->Get(key) ;
  }
  if(cen==1){//centrality 5-10   
    sprintf(key,"h%sPhotPhi%s%s%s_cen1",w,d,kind,PID) ;
    h2DR = (TH2F*)f->Get(key) ;
    sprintf(key,"h%sPhotPhi%s%s%s_cen2",w,d,kind,PID) ;
    TH2F* h2DRa = (TH2F*)f->Get(key) ;
    sprintf(key,"h%sPhotPhi%s%s%s_cen3",w,d,kind,PID) ;
    TH2F* h2DRb = (TH2F*)f->Get(key) ;

    TH1D* hCen = ((TH2F*)f->Get("hCentrality"))->ProjectionX();
    Double_t e1 = hCen->Integral(6,8)/3.;
    Double_t e2 = hCen->Integral(9,9);
    Double_t e3 = hCen->Integral(10,10);

cout<<"w2="<<e1/e2<<", w3="<<e1/e3<<endl;

    h2DR->Add(h2DRa,e1/e2); delete h2DRa;
    h2DR->Add(h2DRb,e1/e3); delete h2DRb;

//    h2DR->Add(h2DRa); delete h2DRa;
//    h2DR->Add(h2DRb); delete h2DRb;
  }
  if(cen==2){//centrality 10-20   
    sprintf(key,"h%sPhotPhi%s%s%s_cen4",w,d,kind,PID) ;
    h2DR = (TH2F*)f->Get(key) ;

    sprintf(key,"h%sPhotPhi%s%s%s_cen0",w,d,kind,PID) ;
    TH2F* h2DRa = (TH2F*)f->Get(key) ;
    sprintf(key,"h%sPhotPhi%s%s%s_cen1",w,d,kind,PID) ;
    TH2F* h2DRb = (TH2F*)f->Get(key) ;
    sprintf(key,"h%sPhotPhi%s%s%s_cen2",w,d,kind,PID) ;
    TH2F* h2DRc = (TH2F*)f->Get(key) ;
    sprintf(key,"h%sPhotPhi%s%s%s_cen3",w,d,kind,PID) ;
    TH2F* h2DRd = (TH2F*)f->Get(key) ;


    TH1D* hCen = ((TH2F*)f->Get("hCentrality"))->ProjectionX();
    Double_t e1 = hCen->Integral(11,11);
    Double_t e2 = hCen->Integral(12,12);
    Double_t e3 = hCen->Integral(13,13);
    Double_t e4 = hCen->Integral(14,15)/2.;
    Double_t e5 = hCen->Integral(16,20)/5.;

cout<<"w_11="<<e5/e1<<", w_12="<<e5/e2<<", w_13="<<e5/e3<<", w_14-15="<<e5/e4<<endl;

    h2DR->Add(h2DRa,e5/e1); delete h2DRa;
    h2DR->Add(h2DRb,e5/e2); delete h2DRb;
    h2DR->Add(h2DRc,e5/e3); delete h2DRc;
    h2DR->Add(h2DRd,e5/e4); delete h2DRd;

  }
  if(cen==3){//centrality 20-30   
    sprintf(key,"h%sPhotPhi%s%s%s_cen5",w,d,kind,PID) ;
    h2DR = (TH2F*)f->Get(key) ;
  }
  if(cen==4){//centrality 30-40   
    sprintf(key,"h%sPhotPhi%s%s%s_cen6",w,d,kind,PID) ;
    h2DR = (TH2F*)f->Get(key) ;
  }
  if(cen==5){//centrality 40-50   
    sprintf(key,"h%sPhotPhi%s%s%s_cen7",w,d,kind,PID) ;
    h2DR = (TH2F*)f->Get(key) ;
  }

  if(Real)return h2DR;
  else return h2DM;
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


