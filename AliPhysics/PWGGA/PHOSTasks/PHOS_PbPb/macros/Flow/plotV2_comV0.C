//#include "pointsForPionGroup.C"
#include "v2All.C"
#include "v2CMS.C"
#include "v2PHENIX.C"

Int_t drawch = 0;
Int_t drawall = 0; //V0A and V0C separate
Int_t drawsum = 1; //PHOS V0A+V0C
Int_t drawphenix = 0;
Int_t drawcms = 0;
Int_t drawpcm = 1;
Int_t drawwhdg = 0;
Int_t drawmerge = 1; //PHOS+PCM
Int_t merge = 1;
Int_t fit = 0;
Int_t diff = 2; //1 - PHOS+PCM - charged (fit), 2 - PHOS - PHOS+PCM, 3 - PCM - PHOS+PCM, 4 - PHOS - charged (fit), 5 - PCM - charged (fit)

//Compare pi0 v2 with V0 detector (method 1 or 2) with other measurements - charged pion v2, CMS pi0 v2, PHENIX pi0 v2
//Merge results from PHOS and PCM

void plotV2_comV0(Int_t cen=1, Int_t method=1){
  char key[255];
  gStyle->SetErrorX(0.2);

if(diff==1||diff==4||diff==5){
drawch=1;
fit=1;
}

  TGraphAsymmErrors * v2ch; // = (TGraphErrors*)f3->Get(cen3.Data()) ;

  if(cen==10) v2ch=v2Pion(0,10);
  if(cen==11) v2ch=v2Pion(20,40);
  if(cen==0) v2ch=v2Pion(0,5);
  if(cen==1) v2ch=v2Pion(5,10);
  if(cen==2) v2ch=v2Pion(10,20);
  if(cen==3) v2ch=v2Pion(20,30);
  if(cen==4) v2ch=v2Pion(30,40);
  if(cen==5) v2ch=v2Pion(40,50);
  if(cen==20) v2ch=v2Pion(10,50);
  if(cen==21) v2ch=v2Pion(0,20);

  v2ch->SetMarkerStyle(7) ;
  v2ch->SetMarkerColor(kOrange-2) ;
  v2ch->SetLineColor(kOrange-2) ;
  v2ch->SetFillColor(kOrange-2) ;
  v2ch->SetLineWidth(1) ;

  TGraphAsymmErrors * v2chAlex; // = (TGraphErrors*)f3->Get(cen3.Data()) ;

  if(cen==10) v2chAlex=v2PionAlex(0,10);
  if(cen==11) v2chAlex=v2PionAlex(20,40);
  if(cen==0) v2chAlex=v2PionAlex(0,5);
  if(cen==1) v2chAlex=v2PionAlex(5,10);
  if(cen==2) v2chAlex=v2PionAlex(10,20);
  if(cen==3) v2chAlex=v2PionAlex(20,30);
  if(cen==4) v2chAlex=v2PionAlex(30,40);
  if(cen==5) v2chAlex=v2PionAlex(40,50);
  if(cen==20) v2chAlex=v2PionAlex(10,50);
  if(cen==21) v2chAlex=v2PionAlex(0,20);

  v2chAlex->SetMarkerStyle(22) ;
  v2chAlex->SetMarkerColor(kYellow) ;
  v2chAlex->SetLineColor(kYellow) ;
  v2chAlex->SetFillColor(kYellow) ;
  v2chAlex->SetLineWidth(1) ;

  TList *coll = new TList();
  coll->Add(v2chAlex);

  TGraphAsymmErrors *v2chAll = v2ch->Clone("v2chAll");
  v2chAll->Merge(coll);

  TGraphErrors *v2CMS; // = (TGraphErrors*)f3->Get(cen3.Data()) ;

  TFile *fPCM = new TFile("PCMPi0v2_new.root") ;

  TGraphErrors *v2PCM, *v2PCMsys;

  if(cen==11){
    v2PCM = (TGraphErrors*)fPCM->Get("Pi0IMv2_V0_20-40");
    v2PCMsys = (TGraphErrors*)fPCM->Get("Pi0IMv2sys_V0_20-40");
  }
  if(cen==21){
    v2PCM = (TGraphErrors*)fPCM->Get("Pi0IMv2_V0_0-20");
    v2PCMsys = (TGraphErrors*)fPCM->Get("Pi0IMv2sys_V0_0-20");
  }

  if(cen==11||cen==21){
  v2PCM->SetMarkerStyle(29) ;
  v2PCM->SetMarkerColor(kBlue) ;
  v2PCM->SetMarkerSize(1.6) ;
  v2PCM->SetLineColor(kBlue) ;
  v2PCM->SetLineWidth(1) ;
  v2PCMsys->SetMarkerStyle(24) ;
  v2PCMsys->SetMarkerColor(kBlue) ;
  v2PCMsys->SetLineColor(kBlue) ;
  v2PCMsys->SetLineWidth(1) ;
  v2PCMsys->SetFillStyle(0);
  v2PCMsys->SetLineColor(kBlue) ;
  }
  else if(merge){ cout<<"Wrong centrality for merging PCM and PHOS (valid 11 and 21 only!)"<<endl; return;}

//  if(cen==10) v2CMS=v2CMS2030();
//  if(cen==11) v2CMS=
//  if(cen==0) v2CMS=v2PionAlex(0,5);
//  if(cen==1) v2CMS=v2PionAlex(5,10);
//  if(cen==2) v2CMS=v2PionAlex(10,20);
  if(cen==3) v2CMS=v2CMS2030();
  if(cen==4) v2CMS=v2CMS3040();
  if(cen==5) v2CMS=v2CMS4050();

  TGraphErrors * v2PHENIX; //

  if(cen==2) v2PHENIX=v2PHENIX1020();
  if(cen==3) v2PHENIX=v2PHENIX2030();
  if(cen==4) v2PHENIX=v2PHENIX3040();
  if(cen==5) v2PHENIX=v2PHENIX4050();
  if(cen==10) v2PHENIX=v2PHENIX010();

if(cen==3||cen==4||cen==5){
  v2CMS->SetMarkerStyle(23) ;
  v2CMS->SetMarkerColor(kViolet) ;
  v2CMS->SetLineColor(kViolet) ;
  v2CMS->SetLineWidth(1) ;
}
if(cen==2||cen==3||cen==4||cen==5||cen==10){
  v2PHENIX->SetMarkerStyle(23) ;
  v2PHENIX->SetMarkerColor(kBlue) ;
  v2PHENIX->SetLineColor(kBlue) ;
  v2PHENIX->SetLineWidth(1) ;
}

fstream fWHDG;
char wname[255];
if(cen==2) sprintf(wname,"WHDGlhcPi0v21020b.dat");
if(cen==3) sprintf(wname,"WHDGlhcPi0v22030b.dat");
if(cen==4) sprintf(wname,"WHDGlhcPi0v23040b.dat");
if(cen==5) sprintf(wname,"WHDGlhcPi0v24050b.dat");

if(cen>=2&&cen<=5){
fWHDG.open(wname,ios::in);

const int N=4;
Double_t xW[N],yW[N],exW[N],eyW[N];

int i=0;
while(!fWHDG.eof()){
  fWHDG>>xW[i]>>yW[i];
  exW[i]=0;
  eyW[i]=0;
  i++;
  if(i>=N) break;
}
fWHDG.close();

TGraphErrors * gWHDG = new TGraphErrors(N,xW,yW,exW,eyW);
gWHDG->SetLineColor(kGreen-5);
gWHDG->SetLineWidth(3);
gWHDG->SetMarkerStyle(20);
gWHDG->SetMarkerColor(kGreen-5);
}

  char name[255];
  sprintf(name,"v2_method%d_QA.root",method); 
  TFile *f = new TFile(name) ;

  char key[255] ;
  char kind[15];
  sprintf(kind,"Both2core"); //PID
  char w[25];
  sprintf(w,""); //weight - "" or "NW"

  sprintf(key,"v2stat_m%d_%s_%s_%d_%d",method,w,kind,cen,0) ;
  cout<<key<<endl;
  TH1D * v2stat0 = (TH1D*)f->Get(key) ;
  sprintf(key,"v2sys_m%d_%s_%s_%d_%d",method,w,kind,cen,0) ;
  TH1D * v2sys0 = (TH1D*)f->Get(key) ;

  sprintf(key,"v2stat_m%d_%s_%s_%d_%d",method,w,kind,cen,1) ;
  TH1D * v2stat1 = (TH1D*)f->Get(key) ;
  sprintf(key,"v2sys_m%d_%s_%s_%d_%d",method,w,kind,cen,1) ;
  TH1D * v2sys1 = (TH1D*)f->Get(key) ;

  v2sys0->SetFillStyle(0);
  v2sys0->SetLineColor(kBlue-9) ;

  v2sys1->SetFillStyle(0) ;
  v2sys1->SetLineColor(kGreen) ;

  TH1D * box = new TH1D("box","",100,0.,16.) ;
  if(cen==0||cen==1||cen==10) box->SetMinimum(-0.20001) ;
  else box->SetMinimum(-0.10001) ;
  box->SetMaximum(0.5) ;

/*  if(diff){
    box->SetMinimum(-0.10001) ;
    box->SetMaximum(0.1) ;
  }
*/

  box->SetLabelOffset(0.01,"Y");
  box->SetLabelOffset(0.01,"X");
  box->SetTitleOffset(0.5,"Y");
  box->SetTitleSize(0.05,"X");
  box->SetTitleSize(0.05,"Y");
  box->SetLabelSize(0.04,"X");
  box->SetLabelSize(0.04,"Y");
  box->SetNdivisions(505,"Y");
  box->GetXaxis()->SetTitleOffset(0.9) ;
  box->GetYaxis()->SetTitleOffset(0.9) ;

  box->SetXTitle("p_{T} (GeV/c)") ;
  box->SetYTitle("v_{2}{EP}") ;

  TCanvas * c= new TCanvas("a","v2") ;
  c->SetFillColor(0) ;
  c->SetFillStyle(0) ;
  c->SetFrameBorderMode(0);
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetTitleFillColor(kWhite);
  gStyle->SetPalette(1);
  gStyle->SetCanvasBorderMode(-1);
  gStyle->SetCanvasBorderSize(1);
  gStyle->SetFrameFillColor(0);
  gStyle->SetFrameFillStyle(0);
  gStyle->SetFrameBorderSize(1);

  Double_t x[100],y[100],ex[100],ey[100] ;

  for(Int_t i=0;i<v2stat0->GetNbinsX();i++){
    x[i]=v2stat0->GetXaxis()->GetBinCenter(i+1) ;
    y[i]=v2stat0->GetBinContent(i+1) ;
    ey[i]=v2stat0->GetBinError(i+1) ;
    ex[i]=0. ;
  }

  TGraphErrors * g = new TGraphErrors(v2stat0->GetNbinsX(),x,y,ex,ey) ;
  g->SetMarkerStyle(20) ;
  g->SetMarkerColor(4) ;
  g->SetLineColor(4) ;
  g->SetLineWidth(1) ;


 for(Int_t i=0;i<v2stat1->GetNbinsX();i++){
    x[i]=v2stat1->GetXaxis()->GetBinCenter(i+1) ;
    y[i]=v2stat1->GetBinContent(i+1) ;
    ey[i]=v2stat1->GetBinError(i+1) ;
    ex[i]=0. ;
  }

  TGraphErrors * g2 = new TGraphErrors(v2stat1->GetNbinsX(),x,y,ex,ey) ;
  g2->SetMarkerStyle(20) ;
  g2->SetMarkerColor(kBlack) ;
  g2->SetLineColor(kBlack) ;
  g2->SetLineWidth(1) ;

//merge V0A V0C
  TH1D *v2sumstat = v2stat0->Clone("v2sumstat");
  TH1D *v2sumsys = v2sys0->Clone("v2sumsys");
cout<<"merge v0a and v0c"<<endl;
  for(Int_t i=0;i<v2sumstat->GetNbinsX();i++){
    Double_t bin1 = v2stat1->GetBinContent(i+1) ;
    Double_t bin2 = v2stat0->GetBinContent(i+1) ;
    Double_t stat1 = v2stat1->GetBinError(i+1) ;
    Double_t stat2 = v2stat0->GetBinError(i+1) ;
    Double_t sys1 = v2sys1->GetBinError(i+1) ;
    Double_t sys2 = v2sys0->GetBinError(i+1) ;

    Double_t stat=0;
    if(stat1!=0&&stat2!=0)stat = 1./TMath::Sqrt(1./stat1/stat1 + 1./stat2/stat2);
    Double_t sys=0;
    if(sys1!=0&&sys2!=0) sys = 1./TMath::Sqrt(1./sys1/sys1 + 1./sys2/sys2);
    Double_t weight = 1, mean = bin1;
    if(sys1+stat1!=0 && sys2+stat2!=0){
      weight=(1./(sys1*sys1+stat1*stat1) + 1./(sys2*sys2+stat2*stat2));
      mean = bin1/(sys1*sys1+stat1*stat1) + bin2/(sys2*sys2+stat2*stat2);
    }
    mean /= weight;

//if(cen==0 && i==0){mean=-100.; stat=0; sys=0;}

    cout<<"bin1 = "<<bin1<<" bin2 = "<<bin2<<" mean = "<<mean<<" sys1 = "<<sys1<<" sys2 = "<<sys2<<" sys = "<<sys<<endl;

    v2sumstat->SetBinContent(i+1,mean);
    v2sumsys->SetBinContent(i+1,mean);
    v2sumstat->SetBinError(i+1,stat);
    v2sumsys->SetBinError(i+1,sys);
  }

  v2sumsys->SetFillStyle(0) ;
  v2sumsys->SetLineColor(kRed) ;

//merge PHOS and PCM
  if(cen==21 || cen==11){
cout<<"merge PHOS and PCM"<<endl;
  TH1D *v2mstat = v2stat0->Clone("v2mstat");
  TH1D *v2msys = v2sys0->Clone("v2msys");

  for(Int_t i=0;i<v2mstat->GetNbinsX();i++){
    Double_t bin1 = v2sumstat->GetBinContent(i+1) ;
    Double_t bin2,tmp;

    v2PCM->GetPoint(i,tmp,bin2);//GetBinContent(i+1) ;
    Double_t stat1 = v2sumstat->GetBinError(i+1) ;
    Double_t stat2 = v2PCM->GetErrorY(i);//->GetBinError(i+1) ;
    Double_t sys1 = v2sumsys->GetBinError(i+1) ;
    Double_t sys2 = v2PCMsys->GetErrorY(i);//GetBinError(i+1) ;

    Double_t stat=0;
    if(stat1!=0&&stat2!=0)stat = 1./TMath::Sqrt(1./stat1/stat1 + 1./stat2/stat2);
    Double_t sys=0;
    if(sys1!=0&&sys2!=0) sys = 1./TMath::Sqrt(1./sys1/sys1 + 1./sys2/sys2);
    Double_t weight = 1, mean = bin1;
    if(sys1+stat1!=0 && sys2+stat2!=0){
      weight=(1./(sys1*sys1+stat1*stat1) + 1./(sys2*sys2+stat2*stat2));
      mean = bin1/(sys1*sys1+stat1*stat1) + bin2/(sys2*sys2+stat2*stat2);
    }
    mean /= weight;

    if(i>7){ //no PCM
	mean=bin1;
	stat=stat1;
	sys=sys1;
    }
//if(cen==0 && i==0){mean=-100.; stat=0; sys=0;}

    cout<<"i = "<<i<<", bin1 = "<<bin1<<" bin2 = "<<bin2<<" mean = "<<mean<<endl;
    cout<<"sys1 = "<<sys1<<" sys2 = "<<sys2<<" sys = "<<sys<<" stat1 = "<<stat1<<" stat2 = "<<stat2<<" stat = "<<stat<<endl;

/*if(cen==21&&i==0){
    v2mstat->SetBinContent(i+1,0);
    v2msys->SetBinContent(i+1,0);
    v2mstat->SetBinError(i+1,0);
    v2msys->SetBinError(i+1,0);
    v2sumstat->SetBinContent(i+1,0);
    v2sumsys->SetBinContent(i+1,0);
    v2sumstat->SetBinError(i+1,0);
    v2sumsys->SetBinError(i+1,0);
}
else{
*/
    v2mstat->SetBinContent(i+1,mean);
    v2msys->SetBinContent(i+1,mean);
    v2mstat->SetBinError(i+1,stat);
    v2msys->SetBinError(i+1,sys);
//}
  }

  v2msys->SetFillStyle(0) ;
  v2msys->SetLineColor(kBlack) ;
  v2msys->SetLineWidth(1) ;

  for(Int_t i=0;i<v2mstat->GetNbinsX();i++){
    x[i]=v2mstat->GetXaxis()->GetBinCenter(i+1) ;
    y[i]=v2mstat->GetBinContent(i+1) ;
    ey[i]=v2mstat->GetBinError(i+1) ;
    ex[i]=0. ;
  }

  TGraphErrors * gm = new TGraphErrors(v2mstat->GetNbinsX(),x,y,ex,ey) ;
  gm->SetMarkerStyle(24) ;
  gm->SetMarkerColor(kBlack);
  gm->SetLineColor(kBlack) ;
  gm->SetLineWidth(1) ;
  gm->SetMarkerSize(1);
  }

  for(Int_t i=0;i<v2sumstat->GetNbinsX();i++){
    x[i]=v2sumstat->GetXaxis()->GetBinCenter(i+1) ;
    y[i]=v2sumstat->GetBinContent(i+1) ;
    ey[i]=v2sumstat->GetBinError(i+1) ;
    ex[i]=0. ;
  }

  TGraphErrors * gsum = new TGraphErrors(v2sumstat->GetNbinsX(),x,y,ex,ey) ;
  gsum->SetMarkerStyle(21) ;
  gsum->SetMarkerColor(kRed);
  gsum->SetLineColor(kRed) ;
  gsum->SetLineWidth(1) ;
  gsum->SetMarkerSize(1.5);

//DRAW ALL!!!
  box->Draw() ;

  if(drawch){
    v2ch->Draw("3");
    v2chAlex->Draw("3 same");
  }
//if(cen==3||cen==4||cen==5) v2CMS->Draw("p same");

  if(drawall){ 
    v2sys0->Draw("E2same");
    v2sys1->Draw("E2same");
    g->Draw("p") ;
    g2->Draw("p");
  }

  if(drawwhdg&&cen>=2&&cen<=5){
   gWHDG->Draw("pl");
  }

  if(fit){
//  TF1* ff=new TF1("ff","[0]*(1-exp(-[1]*x))/([2]+exp([3]+[4]*x+[5]*x*x))");
  TF1* ff=new TF1("ff","([0]+exp([1]+[2]*x+[3]*x*x))/([4]+exp([5]+[6]*x+[7]*x*x))",0,20);

  ff->SetParameters(0,0.05);
  ff->SetParameters(1,0.05);

  v2chAll->Fit(ff,"m");
  ff->Draw("same");

  Double_t *p = ff->GetParameters();
  cout<<p[0]<<"*(1-exp(-1*("<<p[1]<<")*x))/("<<p[2]<<"+exp("<<p[3]<<"+("<<p[4]<<")*x+("<<p[5]<<")*x*x))"<<endl;

  TGraphAsymmErrors* v2ch_cpy = v2ch->Clone("v2ch_cpy");
  TGraphAsymmErrors* v2ch_diff = v2ch->Clone("v2ch_diff");

  v2ch_cpy->Apply(ff);
  }

//PHOS sum V0A and V0C
  if(drawsum){
    gsum->Draw("p");
    v2sumsys->Draw("E2same");

  TFile fout("v2PHOS.root","update");
  char nname[100];
  char cname[100]="";
  if(cen==10) sprintf(cname,"cen010");
  if(cen==11) sprintf(cname,"cen2040");
  if(cen==20) sprintf(cname,"cen1050");

  if(cen==0) sprintf(cname,"cen05");
  if(cen==1) sprintf(cname,"cen510");
  if(cen==2) sprintf(cname,"cen1020");
  if(cen==3) sprintf(cname,"cen2030");
  if(cen==4) sprintf(cname,"cen3040");
  if(cen==5) sprintf(cname,"cen4050");

  sprintf(nname,"v2statV0_%s_11h",cname) ;
  gsum->SetName(nname) ;
  sprintf(nname,"v2sysV0_%s_11h",cname) ;
  v2sumsys->SetName(nname) ;

  gsum->Write(0,TObject::kOverwrite) ;
  v2sumsys->Write(0,TObject::kOverwrite) ;
  fout.Close() ;

  }

  if(fit){
  for(int i=0;i<v2ch_cpy->GetN();i++){
    Double_t px1,px2,py1,py2,pyPHOS,pyPCM,py,px;
    v2ch_cpy->GetPoint(i,px1,py1); //PHOS+PCM
    v2ch->GetPoint(i,px2,py2); //PHOS+PCM

    v2ch_diff->SetPoint(i,px1,py1-py2);
  }
  v2ch_diff->SetLineColor(kGray+1);
  v2ch_diff->SetMarkerColor(kGray+1);

  v2ch_diff->Draw("3");
  }

  if(drawpcm){
  if(cen==11||cen==21){
    v2PCMsys->Draw("E2same");
    v2PCM->Draw("p");
  }
  }

//  if(merge){
    if(diff){
      TGraphErrors *gtmp = gsum->Clone("gtmp");

//if(merge)gtmp=gm->Clone("gtmp");
//else gtmp=gsum->Clone("gtmp");

      TGraphErrors *gmnew;
      TH1D *v2msysnew;
      TGraphErrors *gmsysnew;

if(diff==1){
	gmnew=(TGraphErrors*)gm->Clone("gmnew");
	v2msysnew=(TH1D*)v2msys->Clone("v2msysnew");
}
else if(diff==2||diff==4){
        gmnew=(TGraphErrors*)gsum->Clone("gmnew");
        v2msysnew=(TH1D*)v2sumsys->Clone("v2msysnew");
}
else if(diff==3||diff==5){
        gmnew=(TGraphErrors*)v2PCM->Clone("gmnew");
        gmsysnew=(TGraphErrors*)v2PCMsys->Clone("v2msysnew");
}

//      TH1D *v2msysnew = v2msys->Clone("v2msysnew");

      if(fit)gtmp->Apply(ff);
      for(int i=0;i<gmnew->GetN();i++){
	Double_t px1=0,px2=0,py1=0,py2=0,pyPHOS=0,pyPCM=0,py=0,px=0;
	if(diff==1||diff==2||diff==3)gm->GetPoint(i,px1,py1); //PHOS+PCM
	if(fit)gtmp->GetPoint(i,px2,py2); //ch fit

	if(px1!=px2) cout<<"ERROR!"<<endl;
cout<<"points: "<<py1<<" "<<py2<<endl;

        py=py1-py2;
	gsum->GetPoint(i,px,pyPHOS);
	if(merge)v2PCM->GetPoint(i,px,pyPCM);
	 if(i>7) pyPCM=0;
        if(diff==2){ //PHOS - PHOS+PCM
                py=pyPHOS-py1;
        }
        if(diff==3){ //PCM - PHOS+PCM
                py=pyPCM-py1;
        }
        if(diff==4){ //PHOS - ch
                py=pyPHOS-py2;
        }
        if(diff==5){ //PCM - ch
                py=pyPCM-py2;
        }

	gmnew->SetPoint(i,px,py);
	if(diff==1||diff==2||diff==4)v2msysnew->SetBinContent(i+1,py);
        if(diff==3||diff==5)gmsysnew->SetPoint(i,px,py);

      }
    gmnew->SetLineColor(kGreen);
    gmnew->SetMarkerColor(kGreen);

    gmnew->Draw("p");
    if(diff==1||diff==2||diff==4){
	v2msysnew->SetLineColor(kGreen);
        v2msysnew->SetMarkerColor(kGreen);
	v2msysnew->Draw("E2same");
    }
    if(diff==3||diff==5){
        gmsysnew->SetLineColor(kGreen);
        gmsysnew->SetMarkerColor(kGreen);
	gmsysnew->Draw("E2same");
    }
//    }

  if(drawmerge){
    gm->Draw("p");
    v2msys->Draw("E2same");
  }

  if(diff!=0){
    TF1* ffuu=new TF1("ffuu","[0]");
    gmnew->Fit(ffuu,"m");

  TPaveText* t3=new TPaveText(0.12,0.08,0.57,0.18,"NDC"); //0.63,0.5,0.8,0.70,"NDC");
  t3->SetFillStyle(0);
  t3->SetBorderSize(0);
  t3->SetTextSize(0.04) ;

  char ffuuname[255];

  sprintf(ffuuname,"fit diff: %f #pm %f",ffuu->GetParameter(0),ffuu->GetParError(0));
  t3->AddText(0.,0.,ffuuname);
  t3->Draw();
  }
  }

  if(drawcms)
    if(cen==3||cen==4||cen==5) v2CMS->Draw("p same");
  if(drawphenix)
    if(cen==2||cen==3||cen==4||cen==5||cen==10) v2PHENIX->Draw("p same");

//TITLE

  TPaveText* t2=new TPaveText(0.12,0.88,0.57,0.98,"NDC"); //0.63,0.5,0.8,0.70,"NDC");
  t2->SetFillStyle(0);
  t2->SetBorderSize(0);
  t2->SetTextSize(0.04) ;

  char tname[255];

  sprintf(tname,"#pi^{0} v_{2}^{EP} centrality 0-10%%");
  if(cen==10) t2->AddText(0.,0.,tname);
  sprintf(tname,"#pi^{0} v_{2}^{EP} centrality 20-40%%");
  if(cen==11) t2->AddText(0.,0.,tname);
  sprintf(tname,"#pi^{0} v_{2}^{EP} centrality 0-5%%");
  if(cen==0) t2->AddText(0.,0.,tname);
  sprintf(tname,"#pi^{0} v_{2}^{EP} centrality 5-10%%");
  if(cen==1) t2->AddText(0.,0.,tname);
  sprintf(tname,"#pi^{0} v_{2}^{EP} centrality 10-20%%");
  if(cen==2) t2->AddText(0.,0.,tname);
  sprintf(tname,"#pi^{0} v_{2}^{EP} centrality 20-30%%");
  if(cen==3) t2->AddText(0.,0.,tname);
  sprintf(tname,"#pi^{0} v_{2}^{EP} centrality 30-40%%");
  if(cen==4) t2->AddText(0.,0.,tname);
  sprintf(tname,"#pi^{0} v_{2}^{EP} centrality 40-50%%");
  if(cen==5) t2->AddText(0.,0.,tname);
  sprintf(tname,"#pi^{0} v_{2}^{EP} centrality 10-50%%");
  if(cen==20) t2->AddText(0.,0.,tname);
  sprintf(tname,"#pi^{0} v_{2}^{EP} centrality 0-20%%");
  if(cen==21) t2->AddText(0.,0.,tname);

  t2->Draw();

  TLegend * l = new TLegend(0.15,0.68,0.95,0.88) ;
// TLegend * l = new TLegend(0.22,0.08,0.75,0.28) ;

  l->SetFillColor(0) ;
  l->SetTextSize(0.04) ;


  if(drawsum)
	if(method==1) l->AddEntry(gsum,"PHOS #pi^{0} dNdPhi method, VZERO EP 3 subs","p") ;
	else if(method==2) l->AddEntry(gsum,"PHOS #pi^{0} inv mass method, VZERO EP 3 subs","p") ;

  if(drawpcm)
  if(cen==11||cen==21)l->AddEntry(v2PCM,"PCM #pi^{0}, VZERO EP 3 subs","p") ;

  if(drawall){
    l->AddEntry(g,"PHOS #pi^{0}, V0A EP 3 subs","p") ;
    l->AddEntry(g2,"PHOS #pi^{0}, V0C EP 3 subs","p") ;
  }

  if(merge){
    l->AddEntry(gm,"PHOS+PCM #pi^{0}, V0 EP 3 subs","p") ;
    if(diff==1)l->AddEntry(gmnew,"PHOS+PCM #pi^{0} - #pi^{#pm}, V0 EP","p") ;
    if(diff==2)l->AddEntry(gmnew,"PHOS - PHOS+PCM #pi^{0}, V0 EP 3 subs","p") ;
    if(diff==3)l->AddEntry(gmnew,"PCM - PHOS+PCM #pi^{0}, V0 EP 3 subs","p") ;
    if(diff==4)l->AddEntry(gmnew,"PHOS #pi^{0} - #pi^{#pm}, V0 EP","p") ;
    if(diff==5)l->AddEntry(gmnew,"PCM #pi^{0} - #pi^{#pm}, V0 EP","p") ;
  }
  if(drawch){
    l->AddEntry(v2ch,"TOF #pi^{#pm}, V0 EP 3 subs","f") ;
    l->AddEntry(v2chAlex,"TPC #pi^{#pm} (v2Alex)","f") ;
  }
if(drawcms)
if(cen==3||cen==4||cen==5)l->AddEntry(v2CMS,"CMS #pi^{0}","p");
if(drawphenix)
if(cen==2||cen==3||cen==4||cen==5||cen==10) l->AddEntry(v2PHENIX,"PHENIX #pi^{0}, MPC+RXNin","p");
if(drawwhdg)
if(cen>=2&&cen<=5)l->AddEntry(gWHDG,"WHDG predictions for #pi^{0} v_{2} at LHC energies","p");

  l->Draw() ;

//if(cen==4) ALICEPreliminary2(c);
//  else ALICEPreliminary(c);
  
//  if(cen==6)TGraphErrors * gph = new TGraphErrors(17,pT2,v22,xbins,v2e2);
//  if(cen==3)TGraphErrors * gph = new TGraphErrors(17,pT4,v24,xbins,v2e4);
//  if(cen==4)TGraphErrors * gph = new TGraphErrors(17,pT6,v26,xbins,v2e6);

 if(diff==0){
    if(drawall)
    sprintf(name,"combined_%s_%s_%d_%d_V0EPall.gif",w,kind,cen,method) ;
    else
    sprintf(name,"combined_%s_%s_%d_%d_V0EP.gif",w,kind,cen,method) ;
 }
 else{
    if(drawall)
    sprintf(name,"combined_%s_%s_%d_%d_V0EPall_diff%d.gif",w,kind,cen,method,diff) ;
    else
    sprintf(name,"combined_%s_%s_%d_%d_V0EP_diff%d.gif",w,kind,cen,method,diff) ;
}
    c->SaveAs(name) ;

 if(diff==0){
    if(drawall)
    sprintf(name,"combined_%s_%s_%d_%d_V0EPall.eps",w,kind,cen,method) ;
    else
    sprintf(name,"combined_%s_%s_%d_%d_V0EP.eps",w,kind,cen,method) ;
 }
 else{
    if(drawall)
    sprintf(name,"combined_%s_%s_%d_%d_V0EPall_diff%d.eps",w,kind,cen,method,diff) ;
    else
    sprintf(name,"combined_%s_%s_%d_%d_V0EP_diff%d.eps",w,kind,cen,method,diff) ;
}

    c->SaveAs(name) ;


}
//-----------------------------------------------------------------------------
void ALICEPreliminary(TPad *c){

  TPad *myPadLogo = new TPad("myPadLogo", "Pad for ALICE Logo",0.69,0.70,0.84,0.89);//0.65,0.70,0.8,0.89);
  myPadLogo->SetFillColor(0); 
  myPadLogo->SetBorderMode(0);
  myPadLogo->SetBorderSize(2);
  myPadLogo->SetFrameBorderMode(0);
  myPadLogo->SetLeftMargin(0.0);
  myPadLogo->SetTopMargin(0.0);
  myPadLogo->SetBottomMargin(0.0);
  myPadLogo->SetRightMargin(0.0);
  myPadLogo->SetFillStyle(0);
  myPadLogo->Draw();
  myPadLogo->cd();

  TASImage *myAliceLogo = new TASImage("./alice_logo.png");
  myAliceLogo->Draw();

  c->cd() ;
  TPaveText* t1=new TPaveText(0.67,0.6,0.84,0.70,"NDC"); //0.63,0.5,0.8,0.70,"NDC");
  t1->SetFillStyle(0);
  t1->SetBorderSize(0);
  t1->SetTextSize(0.03) ;
  t1->AddText(0.,0.,"ALICE performance");
  t1->AddText(0.,0.,"#color[1]{18/05/2011}");
  t1->AddText(0.,0.,"Pb-Pb->#pi^{0}+X @ #sqrt{s_{NN}}=2.76 TeV");

//  t1->AddText(0.,0.,"9.8#times 10^{6} events");
//  t1->AddText(0.,0.,"PHOS");
  t1->SetTextColor(kRed);
  t1->SetTextFont(42);
  t1->Draw();
}

//-----------------------------------------------------------------------------
void ALICEPreliminary2(TPad *c){

  TPad *myPadLogo = new TPad("myPadLogo", "Pad for ALICE Logo",0.19,0.70,0.34,0.89);//0.65,0.70,0.8,0.89);

  myPadLogo->SetFillColor(0); 
  myPadLogo->SetBorderMode(0);
  myPadLogo->SetBorderSize(2);
  myPadLogo->SetFrameBorderMode(0);
  myPadLogo->SetLeftMargin(0.0);
  myPadLogo->SetTopMargin(0.0);
  myPadLogo->SetBottomMargin(0.0);
  myPadLogo->SetRightMargin(0.0);
  myPadLogo->SetFillStyle(0);
  myPadLogo->Draw();
  myPadLogo->cd();

  TASImage *myAliceLogo = new TASImage("./alice_logo.png");
  myAliceLogo->Draw();

  c->cd() ;
  TPaveText* t1=new TPaveText(0.17,0.6,0.34,0.70,"NDC");
  t1->SetFillStyle(0);
  t1->SetBorderSize(0);
  t1->SetTextSize(0.03) ;
  t1->AddText(0.,0.,"ALICE performance");
  t1->AddText(0.,0.,"#color[1]{18/05/2011}");
  t1->AddText(0.,0.,"PbPb->#pi^{0}+X @ #sqrt{s_{NN}}=2.76 TeV");

//  t1->AddText(0.,0.,"9.8#times 10^{6} events");
//  t1->AddText(0.,0.,"PHOS");
  t1->SetTextColor(kRed);
  t1->SetTextFont(42);
  t1->Draw();
}


