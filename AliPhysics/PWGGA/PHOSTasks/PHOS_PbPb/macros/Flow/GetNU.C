int drawMB=0;

void GetNU(int side=0){

gStyle->SetOptStat(0);

  TFile * f = new TFile("flow11h_QA.root") ;

  TFile * f2 = new TFile("flow11hMB_QA.root") ;

  char what[20];
  if(side==0) sprintf(what,"phiRPV0Aflat") ;
  if(side==1) sprintf(what,"phiRPV0Cflat") ;
  if(side==2) sprintf(what,"phiRPflat") ;

  TH2F* h=f->Get(what);
  TH2F* h2=f2->Get(what);

  TH3D* hPHOSphi=f->Get("hPHOSphi");

  TF1 *fu=new TF1("fu","[0]*(1+2*[1]*TMath::Cos(2*x)+2*[2]*TMath::Sin(2*x))",1,5);

  TH1D* hcos=new TH1D("cos","cos",10,0.,50.);
  TH1D* hsin=new TH1D("sin","sin",10,0.,50.);
  TH1D* hcosMB=new TH1D("cosMB","cosMB",10,0.,50.);
  TH1D* hsinMB=new TH1D("sinMB","sinMB",10,0.,50.);
  TH1D* hsinv2=new TH1D("sinv2","sinv2",10,0.,50.);

  for(int cen=1;cen<=10;cen++){
  char name[255];
  sprintf(name,"phi_%d",cen);
  TH3D* hphi = hPHOSphi->Clone(name);
  hphi->GetXaxis()->SetRangeUser((cen-1)*5,(cen)*5);
  TH1D* hCos=hphi->Project3D("z");

  double meanC=0,meanS=0,errC=0,errS=0,mean=0;

  for(int i=0;i<hCos->GetXaxis()->GetNbins();i++){
    double bin=hCos->GetBinContent(i+1);
    double err=hCos->GetBinError(i+1);
    double phi=hCos->GetBinCenter(i+1);
    mean+=bin;
    meanC+=TMath::Cos(2*phi)*bin;
    meanS+=TMath::Sin(2*phi)*bin;
    errC+=TMath::Cos(2*phi)*err;
    errS+=TMath::Cos(2*phi)*err;
  }
  meanC/=mean;
  meanS/=mean;
  errC/=mean;
  errS/=mean;

  cout<<"i: "<<cen<<", cos: "<<meanC<<"+-"<<errC<<", sin: "<<meanS<<"+-"<<errS<<endl;

    char name[255];
    sprintf(name,"EP_cen%d",cen);
    TH1D * h1D = h->ProjectionX(name,cen,cen) ;
    sprintf(name,"MB_cen%d",cen);
    TH1D * h1DMB = h2->ProjectionX(name,cen,cen) ;

    h1D->Fit(fu,"q0");

    double sinv2 = meanS*fu->GetParameter(1) - meanC*fu->GetParameter(2);
    double err = meanS*fu->GetParError(1) + errS*fu->GetParameter(1);

    hsinv2->SetBinContent(cen,sinv2);
    hsinv2->SetBinError(cen,err);

    hcos->SetBinContent(cen,fu->GetParameter(1));
    hcos->SetBinError(cen,fu->GetParError(1));
    hsin->SetBinContent(cen,fu->GetParameter(2));
    hsin->SetBinError(cen,fu->GetParError(2));

    h1DMB->Fit(fu,"q0");

    hcosMB->SetBinContent(cen,fu->GetParameter(1));
    hcosMB->SetBinError(cen,fu->GetParError(1));
    hsinMB->SetBinContent(cen,fu->GetParameter(2));
    hsinMB->SetBinError(cen,fu->GetParError(2));
  }

  if(side==0)  hcos->SetTitle("V0A EP");
  if(side==1)  hcos->SetTitle("V0C EP");
  if(side==2)  hcos->SetTitle("TPC EP");

  hcos->SetMarkerStyle(20) ; 
  hcos->SetMarkerColor(2) ; 
  hcos->SetLineColor(2) ; 
  hcos->SetLineWidth(3);

  hcos->SetMinimum(-0.01) ;
  hcos->SetMaximum(0.01) ;
  hcos->SetXTitle("centrality") ;
  hcos->SetYTitle("<cos2#Psi> and <sin2#Psi") ;

  hsin->SetMarkerStyle(20) ;
  hsin->SetMarkerColor(3) ;
  hsin->SetLineColor(3) ;
  hsin->SetLineWidth(2);

  hcosMB->SetMarkerStyle(25) ;
  hcosMB->SetMarkerColor(2) ;
  hcosMB->SetLineColor(2) ;
  hcosMB->SetLineWidth(3);

  hsinMB->SetMarkerStyle(25) ;
  hsinMB->SetMarkerColor(3) ;
  hsinMB->SetLineColor(3) ;
  hsinMB->SetLineWidth(2);

  hsinv2->SetMarkerStyle(20) ;
  hsinv2->SetMarkerColor(4) ;
  hsinv2->SetLineColor(4) ;
  hsinv2->SetLineWidth(2);


  hcos->Draw("pl") ;
  hsin->Draw("pl same") ;
  hsinv2->Draw("pl same");
  if(drawMB){
    hcosMB->Draw("pl same") ;
    hsinMB->Draw("pl same") ;
  }

  TLegend * leg = new TLegend(0.6,0.7,0.9,0.9) ;
  leg->AddEntry(hcos,"cos2#Psi","p") ;
  leg->AddEntry(hsin,"sin2#Psi","p") ;
  leg->AddEntry(hsinv2,"sin2#phi*cos2#Psi - cos2#phi*sin2#Psi","p") ;

  if(drawMB){
    leg->AddEntry(hcosMB,"cos2#Psi, kMB","p") ;
    leg->AddEntry(hsinMB,"sin2#Psi, kMB","p") ;
  }
  leg->Draw() ;
}
