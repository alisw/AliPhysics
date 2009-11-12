#if !defined(__CINT__) || defined(__MAKECINT__)
#include <Riostream.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TFile.h>
#include <TLine.h>
#include <TNtuple.h>
#include <TStyle.h>
#include <TString.h>
#include <TSystem.h>
#include <TH1F.h>
#include <TGraphErrors.h>
#include <TF1.h>
#endif

void PlotVertexESD(TString vtxtype="SPD",
		   TString fname="Vertex.Performance.root",
		   TString ntname="fNtupleVertexESD",
		   Bool_t useztrue=kTRUE,
		   Int_t optgif=0) {
  //-------------------------------------------------------------------------
  // 
  // Plot output of AliAnalysisTaskVertexESD.
  // Origin: francesco.prino@to.infn.it
  //
  //-------------------------------------------------------------------------

  // ranges for ideal geometry
  /*
  Float_t rangeHPull=40.;
  Float_t rangeHRes=2500.;
  Int_t binsHRes=250;
  Float_t rangeGrAve = 100; // micron
  Float_t rangeGrRms = 900; // micron
  Float_t rangeGrPull = 5; 
  */

  // ranges for misaligned geometry
  Float_t rangeHPull=40.;
  Float_t rangeHRes=10000.*0.5;
  Int_t binsHRes=250*0.5;
  Float_t rangeGrAve = 10000; // micron
  Float_t rangeGrRms = 5000; // micron
  Float_t rangeGrPull = 15; 


  Float_t limcont[10]={0.,2.5,3.5,4.5,6.5,9.5,14.5,20.5,99999.};
  const Int_t nbinsmult=7;
  Float_t limmult[nbinsmult+1]={0.,10.,14.,24.,30.,44.,60.,99999.};
  //const Int_t nbinseff=11;
  //Float_t limeff[nbinseff+1]={0.,0.5,1.5,2.5,3.5,4.5,5.5,7.5,10.5,13.,15.,999999.};
  const Int_t nbinseff=9;
  Float_t limeff[nbinseff+1]={0.,0.5,1.5,2.5,3.5,4.5,5.5,7.5,10.5,999999.};
  const Int_t nbinz = 14;
  Float_t limitzt[nbinz+1]={-20,-18,-15.,-12.,-9.,-6.,-3.,0.,3.,6.,9.,12.,15.,18,20};
  //Float_t limitzt[nbinz+1]={-3,0.,3.,6.,9.,11.,13.,15,17,19,21};

  TH1F **hxm=new TH1F*[nbinsmult];
  TH1F **hym=new TH1F*[nbinsmult];
  TH1F **hzm=new TH1F*[nbinsmult];
  TH1F **hxpullm=new TH1F*[nbinsmult];
  TH1F **hypullm=new TH1F*[nbinsmult];
  TH1F **hzpullm=new TH1F*[nbinsmult];
  TH1F **hmult=new TH1F*[nbinsmult];

  TH1F **hxc=new TH1F*[8];
  TH1F **hyc=new TH1F*[8];
  TH1F **hzc=new TH1F*[8];
  TH1F **hxpullc=new TH1F*[8];
  TH1F **hypullc=new TH1F*[8];
  TH1F **hzpullc=new TH1F*[8];
  TH1F **hcont=new TH1F*[8];

  TH1F **hxz=new TH1F*[nbinz];
  TH1F **hyz=new TH1F*[nbinz];
  TH1F **hzz=new TH1F*[nbinz];
  TH1F **hz=new TH1F*[nbinz];

  TH1F **hmulteff=new TH1F*[nbinseff];
  TH1F **haux=new TH1F*[nbinseff];
  TH1F *htot=new TH1F("htot","",100,-1000,1000);

  TH1F **htrkseff=new TH1F*[nbinseff];
  TH1F **htrksaux=new TH1F*[nbinseff];
  TH1F **htrks3Daux=new TH1F*[nbinseff];

  Char_t hisnam[10];
  TCanvas *c1=new TCanvas("c1","Residuals",1000,700);
  c1->Divide(nbinsmult,3,0.001,0.001);
  TCanvas *c1z=new TCanvas("c1z","Residuals z",1000,700);
  c1z->Divide(nbinz,3,0.001,0.001);
  TCanvas *cp=new TCanvas("cp","Pulls",1000,700);
  cp->Divide(nbinsmult,3,0.001,0.001);

  TGraphErrors *gavexm=new TGraphErrors(0);
  TGraphErrors *gaveym=new TGraphErrors(0);
  TGraphErrors *gavezm=new TGraphErrors(0);
  TGraphErrors *grmsxm=new TGraphErrors(0);
  TGraphErrors *grmsym=new TGraphErrors(0);
  TGraphErrors *grmszm=new TGraphErrors(0);
  TGraphErrors *gavexmg=new TGraphErrors(0);
  TGraphErrors *gaveymg=new TGraphErrors(0);
  TGraphErrors *gavezmg=new TGraphErrors(0);
  TGraphErrors *grmsxmg=new TGraphErrors(0);
  TGraphErrors *grmsymg=new TGraphErrors(0);
  TGraphErrors *grmszmg=new TGraphErrors(0);
  TGraphErrors *gpullxm=new TGraphErrors(0);
  TGraphErrors *gpullym=new TGraphErrors(0);
  TGraphErrors *gpullzm=new TGraphErrors(0);
  TGraphErrors *gpullxmg=new TGraphErrors(0);
  TGraphErrors *gpullymg=new TGraphErrors(0);
  TGraphErrors *gpullzmg=new TGraphErrors(0);
  TGraphErrors *geffm=new TGraphErrors(0);
  TGraphErrors *gefftrks=new TGraphErrors(0);
  TGraphErrors *geff3Dtrks=new TGraphErrors(0);

  TGraphErrors *gavexc=new TGraphErrors(0);
  TGraphErrors *gaveyc=new TGraphErrors(0);
  TGraphErrors *gavezc=new TGraphErrors(0);
  TGraphErrors *grmsxc=new TGraphErrors(0);
  TGraphErrors *grmsyc=new TGraphErrors(0);
  TGraphErrors *grmszc=new TGraphErrors(0);
  TGraphErrors *gpullxc=new TGraphErrors(0);
  TGraphErrors *gpullyc=new TGraphErrors(0);
  TGraphErrors *gpullzc=new TGraphErrors(0);

  TGraph * gbeamxz=new TGraph(0);
  TGraph * gbeamyz=new TGraph(0);
  TGraph * gbeamxy=new TGraph(0);
  TH2F *hbeamxz = new TH2F("hbeamxz","",100,-14.,14.,100,-1.5,1.5);
  TH2F *hbeamyz = new TH2F("hbeamyz","",100,-14.,14.,100,-1.5,1.5);
  TH1F *hbeamx = new TH1F("hbeamx","",1000,-1.5,1.5);
  TH1F *hbeamy = new TH1F("hbeamy","",1000,-1.5,1.5);
  TH2F *hbeamxy = new TH2F("hbeamxy","",100,-1.5,1.5,100,-1.5,1.5);



  gavexm->SetName("gavexm");
  grmsxm->SetName("grmsxm");
  gavexmg->SetName("gavexmg");
  grmsxmg->SetName("grmsxmg");
  gpullxm->SetName("gpullxm");
  gpullxmg->SetName("gpullxmg");
  gaveym->SetName("gaveym");
  grmsym->SetName("grmsym");
  gaveymg->SetName("gaveymg");
  grmsymg->SetName("grmsymg");
  gpullym->SetName("gpullym");
  gpullymg->SetName("gpullymg");
  gavezm->SetName("gavezm");
  grmszm->SetName("grmszm");
  gavezmg->SetName("gavezmg");
  grmszmg->SetName("grmszmg");
  gpullzm->SetName("gpullzm");
  gpullzmg->SetName("gpullzmg");
  geffm->SetName("geffm");
  gefftrks->SetName("gefftrks");
  geff3Dtrks->SetName("geff3Dtrks");
  gavexc->SetName("gavexc");
  grmsxc->SetName("grmsxc");
  gpullxc->SetName("gpullxc");
  gaveyc->SetName("gaveyc");
  grmsyc->SetName("grmsyc");
  gpullyc->SetName("gpullyc");
  gavezc->SetName("gavezc");
  grmszc->SetName("grmszc");
  gpullzc->SetName("gpullzc");

  TGraphErrors *gavexz=new TGraphErrors(0);
  TGraphErrors *gaveyz=new TGraphErrors(0);
  TGraphErrors *gavezz=new TGraphErrors(0);
  TGraphErrors *grmsxz=new TGraphErrors(0);
  TGraphErrors *grmsyz=new TGraphErrors(0);
  TGraphErrors *grmszz=new TGraphErrors(0);
  TGraphErrors *gavexzg=new TGraphErrors(0);
  TGraphErrors *gaveyzg=new TGraphErrors(0);
  TGraphErrors *gavezzg=new TGraphErrors(0);
  TGraphErrors *grmsxzg=new TGraphErrors(0);
  TGraphErrors *grmsyzg=new TGraphErrors(0);
  TGraphErrors *grmszzg=new TGraphErrors(0);
  TGraphErrors *geffz=new TGraphErrors(0);

  gavexz->SetName("gavexz");
  grmsxz->SetName("grmsxz");
  gavexzg->SetName("gavexzg");
  grmsxzg->SetName("grmsxzg");
  gaveyz->SetName("gaveyz");
  grmsyz->SetName("grmsyz");
  gaveyzg->SetName("gaveyzg");
  grmsyzg->SetName("grmsyzg");
  gavezz->SetName("gavezz");
  grmszz->SetName("grmszz");
  gavezzg->SetName("gavezzg");
  grmszzg->SetName("grmszzg");
  geffz->SetName("geffz");

  TF1 *fitf;

  // histogram creation
  for(Int_t khis=0;khis<nbinseff;khis++){
    sprintf(hisnam,"hmeff%d",khis);
    hmulteff[khis]=new TH1F(hisnam,"",100,0.,200.);
    sprintf(hisnam,"htrkseff%d",khis);
    htrkseff[khis]=new TH1F(hisnam,"",100,0.,200.);
    sprintf(hisnam,"haux%d",khis);
    haux[khis]=new TH1F(hisnam,"",binsHRes,-rangeHRes,rangeHRes);
    sprintf(hisnam,"htrksaux%d",khis);
    htrksaux[khis]=new TH1F(hisnam,"",binsHRes,-rangeHRes,rangeHRes);
    sprintf(hisnam,"htrks3Daux%d",khis);
    htrks3Daux[khis]=new TH1F(hisnam,"",binsHRes,-rangeHRes,rangeHRes);
  }
  for(Int_t khis=0;khis<nbinsmult;khis++){
    sprintf(hisnam,"hxm%d",khis);
    hxm[khis]=new TH1F(hisnam,"",binsHRes,-rangeHRes,rangeHRes);
    sprintf(hisnam,"hym%d",khis);
    hym[khis]=new TH1F(hisnam,"",binsHRes,-rangeHRes,rangeHRes);
    sprintf(hisnam,"hzm%d",khis);
    hzm[khis]=new TH1F(hisnam,"",binsHRes,-rangeHRes,rangeHRes);
    sprintf(hisnam,"hxpm%d",khis);
    hxpullm[khis]=new TH1F(hisnam,"",100,-rangeHPull,rangeHPull);
    sprintf(hisnam,"hypm%d",khis);
    hypullm[khis]=new TH1F(hisnam,"",100,-rangeHPull,rangeHPull);
    sprintf(hisnam,"hzpm%d",khis);
    hzpullm[khis]=new TH1F(hisnam,"",100,-rangeHPull,rangeHPull);
    sprintf(hisnam,"hmult%d",khis);
    hmult[khis]=new TH1F(hisnam,"",100,0.,200.);
  }

  for(Int_t khis=0;khis<8;khis++){
    sprintf(hisnam,"hxc%d",khis);
    hxc[khis]=new TH1F(hisnam,"",binsHRes,-rangeHRes,rangeHRes);
    sprintf(hisnam,"hyc%d",khis);
    hyc[khis]=new TH1F(hisnam,"",binsHRes,-rangeHRes,rangeHRes);
    sprintf(hisnam,"hzc%d",khis);
    hzc[khis]=new TH1F(hisnam,"",binsHRes,-rangeHRes,rangeHRes);
    sprintf(hisnam,"hxpc%d",khis);
    hxpullc[khis]=new TH1F(hisnam,"",100,-rangeHPull,rangeHPull);
    sprintf(hisnam,"hypc%d",khis);
    hypullc[khis]=new TH1F(hisnam,"",100,-rangeHPull,rangeHPull);
    sprintf(hisnam,"hzpc%d",khis);
    hzpullc[khis]=new TH1F(hisnam,"",100,-rangeHPull,rangeHPull);
    sprintf(hisnam,"hcont%d",khis);
    hcont[khis]=new TH1F(hisnam,"",100,0.,200.);
  }

  for(Int_t khis=0;khis<nbinz;khis++){
    sprintf(hisnam,"hxz%d",khis);
    hxz[khis]=new TH1F(hisnam,"",binsHRes,-rangeHRes,rangeHRes);
    sprintf(hisnam,"hyz%d",khis);
    hyz[khis]=new TH1F(hisnam,"",binsHRes,-rangeHRes,rangeHRes);
    sprintf(hisnam,"hzz%d",khis);
    hzz[khis]=new TH1F(hisnam,"",binsHRes,-rangeHRes,rangeHRes);
    sprintf(hisnam,"hz%d",khis);
    hz[khis]=new TH1F(hisnam,"",100,-20.,20.);   
  }

  Float_t totev=0,totevtriggered=0,nvtx3D=0,nvtxZ=0;
  TFile *f=new TFile(fname.Data());
  TList *cOutput = (TList*)f->Get("cOutput");
  TNtuple *nt=(TNtuple*)cOutput->FindObject(ntname.Data());
  Int_t nnnev=nt->GetEntries();
  printf("Events = %d\n",nnnev);
  Float_t xVtx,xdiffVtx,xerrVtx;
  Float_t yVtx,ydiffVtx,yerrVtx;
  Float_t zVtx,zdiffVtx,zerrVtx;
  Float_t ntrklets,ncontrVtx,dndy,triggered,vtx3D;
  Float_t xtrue,ytrue,ztrue,zref;
  
  TString sxx="x"; sxx.Append(vtxtype.Data());
  nt->SetBranchAddress(sxx.Data(),&xVtx);
  TString syy="y"; syy.Append(vtxtype.Data());
  nt->SetBranchAddress(syy.Data(),&yVtx);
  TString szz="z"; szz.Append(vtxtype.Data());
  nt->SetBranchAddress(szz.Data(),&zVtx);
  
  TString xerr="xerr"; xerr.Append(vtxtype.Data());
  nt->SetBranchAddress(xerr.Data(),&xerrVtx);
  TString yerr="yerr"; yerr.Append(vtxtype.Data());
  nt->SetBranchAddress(yerr.Data(),&yerrVtx);
  TString zerr="zerr"; zerr.Append(vtxtype.Data());
  nt->SetBranchAddress(zerr.Data(),&zerrVtx);

  TString trkstitle;
  if(vtxtype.Contains("TPC")) {
    nt->SetBranchAddress("nTPCin",&ntrklets);
    trkstitle="TPC tracks pointing to beam pipe";
  } else {
    nt->SetBranchAddress("ntrklets",&ntrklets);
    trkstitle="SPD tracklets";
  }
  TString ntrks="ntrks"; ntrks.Append(vtxtype.Data());
  nt->SetBranchAddress(ntrks.Data(),&ncontrVtx);
  nt->SetBranchAddress("dndygen",&dndy);
  
  nt->SetBranchAddress("xtrue",&xtrue);
  nt->SetBranchAddress("ytrue",&ytrue);
  nt->SetBranchAddress("ztrue",&ztrue);

  nt->SetBranchAddress("triggered",&triggered);

  nt->SetBranchAddress("SPD3D",&vtx3D);


  // loop on events
  for(Int_t iev=0;iev<nnnev;iev++){
    nt->GetEvent(iev);

    xdiffVtx=10000.*(xVtx-xtrue);
    ydiffVtx=10000.*(yVtx-ytrue);
    zdiffVtx=10000.*(zVtx-ztrue);


    zref = (useztrue ? ztrue : zVtx);
    if(!vtxtype.Contains("SPD")) vtx3D=1.;
    if(triggered<0.5) continue; // not triggered
    totevtriggered += 1.;
    if(ncontrVtx>0) {
      htot->Fill(zdiffVtx);
      if(vtx3D>0.5) {
	nvtx3D += 1.;
      } else {
	nvtxZ += 1.;
      }
    }

    if(ncontrVtx>0 && vtx3D>0.5) {
      gbeamxz->SetPoint(gbeamxz->GetN(),zVtx,xVtx);
      gbeamyz->SetPoint(gbeamxz->GetN(),zVtx,yVtx);
      gbeamxy->SetPoint(gbeamxz->GetN(),xVtx,yVtx);
      
      hbeamx->Fill(xVtx);
      hbeamy->Fill(yVtx);
      hbeamxz->Fill(zVtx,xVtx);
      hbeamyz->Fill(zVtx,yVtx);
      hbeamxy->Fill(xVtx,yVtx);
    }

    for(Int_t khis=0;khis<nbinseff;khis++){
      if(dndy>=limeff[khis] && dndy<limeff[khis+1]){
	hmulteff[khis]->Fill(dndy);
	if(ncontrVtx>0) haux[khis]->Fill(zdiffVtx);
      }
      if(ntrklets>=limeff[khis] && ntrklets<limeff[khis+1]){
	htrkseff[khis]->Fill(ntrklets);
	if(ncontrVtx>0) htrksaux[khis]->Fill(zdiffVtx);
	if(ncontrVtx>0 && vtx3D>0.5) htrks3Daux[khis]->Fill(zdiffVtx);
      }
    }
    for(Int_t khis=0;khis<nbinsmult;khis++){
      if(ntrklets>=limmult[khis] && ntrklets<limmult[khis+1]){
	hmult[khis]->Fill(ntrklets);
	if(ncontrVtx>0){
	  if(vtx3D>0.5) hxm[khis]->Fill(xdiffVtx);
	  if(vtx3D>0.5) hym[khis]->Fill(ydiffVtx);
	  hzm[khis]->Fill(zdiffVtx);
	  if(vtx3D>0.5) hxpullm[khis]->Fill(xdiffVtx/10000./xerrVtx);
	  if(vtx3D>0.5) hypullm[khis]->Fill(ydiffVtx/10000./yerrVtx);
	  hzpullm[khis]->Fill(zdiffVtx/10000./zerrVtx);
	}
      }
    }
    for(Int_t khis=0;khis<8;khis++){
      if(ncontrVtx>=limcont[khis]&&ncontrVtx<limcont[khis+1]){
	hcont[khis]->Fill(ncontrVtx);
	if(ncontrVtx>0){	
	  if(vtx3D>0.5)hxc[khis]->Fill(xdiffVtx);
	  if(vtx3D>0.5)hyc[khis]->Fill(ydiffVtx);
	  hzc[khis]->Fill(zdiffVtx);
	  if(vtx3D>0.5)hxpullc[khis]->Fill(xdiffVtx/10000./xerrVtx);
	  if(vtx3D>0.5)hypullc[khis]->Fill(ydiffVtx/10000./yerrVtx);
	  hzpullc[khis]->Fill(zdiffVtx/10000./zerrVtx);
	}
      }
    }
    for(Int_t khis=0;khis<nbinz;khis++){
      if(zref>=limitzt[khis]&&zref<limitzt[khis+1]){
	hz[khis]->Fill(zref);
	if(ncontrVtx>0){
	  if(vtx3D>0.5)hxz[khis]->Fill(xdiffVtx);
	  if(vtx3D>0.5)hyz[khis]->Fill(ydiffVtx);
	  hzz[khis]->Fill(zdiffVtx);
	}
      }
    }
  }
  totev+=nnnev;

  if(totev==0){
    printf("Total number of events = 0\n");
    return;
  }
  for(Int_t khis=0;khis<nbinseff;khis++){
    Double_t x=hmulteff[khis]->GetMean();
    Double_t ex=hmulteff[khis]->GetRMS();;
    Float_t nEv=(Float_t)(hmulteff[khis]->GetEntries());
    Float_t trkeff=-1.;
    cout<<"Eff. dNch/dy bin "<<khis<<" ("<<limeff[khis]<<"-"<<limeff[khis+1]<<")  # Events ="<<nEv<<" with vertex ="<<haux[khis]->GetEntries()<<endl;
    if(nEv>0) trkeff=(Float_t)(haux[khis]->GetEntries())/nEv;
    geffm->SetPoint(khis,x,trkeff);
    geffm->SetPointError(khis,ex,0.);
  }
  for(Int_t khis=0;khis<nbinseff;khis++){
    Double_t x=htrkseff[khis]->GetMean();
    Double_t ex=htrkseff[khis]->GetRMS();;
    Float_t nEv=(Float_t)(htrkseff[khis]->GetEntries());
    Float_t trkeff=-1.;
    cout<<"Eff. trks bin "<<khis<<" ("<<limeff[khis]<<"-"<<limeff[khis+1]<<")  # Events ="<<nEv<<" with vertex ="<<htrksaux[khis]->GetEntries()<<endl;
    if(nEv>0) trkeff=(Float_t)(htrksaux[khis]->GetEntries())/nEv;
    gefftrks->SetPoint(khis,x,trkeff);
    gefftrks->SetPointError(khis,ex,0.);
    if(nEv>0) trkeff=(Float_t)(htrks3Daux[khis]->GetEntries())/nEv;
    geff3Dtrks->SetPoint(khis,x,trkeff);
    geff3Dtrks->SetPointError(khis,ex,0.);
  }

  for(Int_t khis=0;khis<nbinsmult;khis++){ 

    c1->cd(khis+1);
    if(hxm[khis]->GetEntries()<10) continue;
    hxm[khis]->Draw();
    hxm[khis]->Fit("gaus","Q0");
    fitf= hxm[khis]->GetFunction("gaus");
    Double_t avexg=fitf->GetParameter(1);
    Double_t eavexg=fitf->GetParError(1);
    Double_t rmsxg=fitf->GetParameter(2);
    Double_t ermsxg=fitf->GetParError(2);
    c1->cd(nbinsmult+khis+1);
    hym[khis]->Draw();
    hym[khis]->Fit("gaus","Q0");
    fitf= hym[khis]->GetFunction("gaus");
    Double_t aveyg=fitf->GetParameter(1);
    Double_t eaveyg=fitf->GetParError(1);
    Double_t rmsyg=fitf->GetParameter(2);
    Double_t ermsyg=fitf->GetParError(2);
    c1->cd(2*nbinsmult+khis+1);
    hzm[khis]->Draw();
    hzm[khis]->Fit("gaus","Q0");
    fitf= hzm[khis]->GetFunction("gaus");
    Double_t avezg=fitf->GetParameter(1);
    Double_t eavezg=fitf->GetParError(1);
    Double_t rmszg=fitf->GetParameter(2);
    Double_t ermszg=fitf->GetParError(2);

    cp->cd(khis+1);
    hxpullm[khis]->Draw();
    hxpullm[khis]->Fit("gaus","Q0");
    fitf= hxpullm[khis]->GetFunction("gaus");
    Double_t pullxg=fitf->GetParameter(2);
    Double_t epullxg=fitf->GetParError(2);
    cp->cd(nbinsmult+khis+1);
    hypullm[khis]->Draw();
    hypullm[khis]->Fit("gaus","Q0");
    fitf= hypullm[khis]->GetFunction("gaus");
    Double_t pullyg=fitf->GetParameter(2);
    Double_t epullyg=fitf->GetParError(2);
    cp->cd(2*nbinsmult+khis+1);
    hzpullm[khis]->Draw();
    hzpullm[khis]->Fit("gaus","Q0");
    fitf= hzpullm[khis]->GetFunction("gaus");
    Double_t pullzg=fitf->GetParameter(2);
    Double_t epullzg=fitf->GetParError(2);
 

    Double_t rmsxt=hxm[khis]->GetRMS();
    Double_t rmsyt=hym[khis]->GetRMS();
    Double_t rmszt=hzm[khis]->GetRMS();
    Double_t ermsxt=hxm[khis]->GetRMSError();
    Double_t ermsyt=hym[khis]->GetRMSError();
    Double_t ermszt=hzm[khis]->GetRMSError();
    Double_t avext=hxm[khis]->GetMean();
    Double_t aveyt=hym[khis]->GetMean();
    Double_t avezt=hzm[khis]->GetMean();
    Double_t eavext=hxm[khis]->GetMeanError();
    Double_t eaveyt=hym[khis]->GetMeanError();
    Double_t eavezt=hzm[khis]->GetMeanError();
    Double_t pullxt=hxpullm[khis]->GetRMS();
    Double_t pullyt=hypullm[khis]->GetRMS();
    Double_t pullzt=hzpullm[khis]->GetRMS();
    Double_t epullxt=hxpullm[khis]->GetRMSError();
    Double_t epullyt=hypullm[khis]->GetRMSError();
    Double_t epullzt=hzpullm[khis]->GetRMSError();

    Double_t x=hmult[khis]->GetMean();
    if(hmult[khis]->GetEntries()==0) x=-1;
    Double_t ex=hmult[khis]->GetRMS();;

    gavexm->SetPoint(khis,x,avext);
    gavexm->SetPointError(khis,ex,eavext);
    gaveym->SetPoint(khis,x,aveyt);
    gaveym->SetPointError(khis,ex,eaveyt);
    gavezm->SetPoint(khis,x,avezt);
    gavezm->SetPointError(khis,ex,eavezt);
    grmsxm->SetPoint(khis,x,rmsxt);
    grmsxm->SetPointError(khis,ex,ermsxt);
    grmsym->SetPoint(khis,x,rmsyt);
    grmsym->SetPointError(khis,ex,ermsyt);
    grmszm->SetPoint(khis,x,rmszt);
    grmszm->SetPointError(khis,ex,ermszt);

    gavexmg->SetPoint(khis,x,avexg);
    gavexmg->SetPointError(khis,ex,eavexg);
    gaveymg->SetPoint(khis,x,aveyg);
    gaveymg->SetPointError(khis,ex,eaveyg);
    gavezmg->SetPoint(khis,x,avezg);
    gavezmg->SetPointError(khis,ex,eavezg);
    grmsxmg->SetPoint(khis,x,rmsxg);
    grmsxmg->SetPointError(khis,ex,ermsxg);
    grmsymg->SetPoint(khis,x,rmsyg);
    grmsymg->SetPointError(khis,ex,ermsyg);
    grmszmg->SetPoint(khis,x,rmszg);
    grmszmg->SetPointError(khis,ex,ermszg);

    gpullxm->SetPoint(khis,x,pullxt);
    gpullxm->SetPointError(khis,ex,epullxt);
    gpullym->SetPoint(khis,x,pullyt);
    gpullym->SetPointError(khis,ex,epullyt);
    gpullzm->SetPoint(khis,x,pullzt);
    gpullzm->SetPointError(khis,ex,epullzt);

    gpullxmg->SetPoint(khis,x,pullxg);
    gpullxmg->SetPointError(khis,ex,epullxg);
    gpullymg->SetPoint(khis,x,pullyg);
    gpullymg->SetPointError(khis,ex,epullyg);
    gpullzmg->SetPoint(khis,x,pullzg);
    gpullzmg->SetPointError(khis,ex,epullzg);

    Float_t nEv=hmult[khis]->GetEntries();
    cout<<"Mult. bin "<<khis<<"  # Events ="<<nEv<<endl;
  }

  for(Int_t khis=0;khis<8;khis++){

    Double_t rmsxt=hxc[khis]->GetRMS();
    Double_t rmsyt=hyc[khis]->GetRMS();
    Double_t rmszt=hzc[khis]->GetRMS();
    Double_t ermsxt=hxc[khis]->GetRMSError();
    Double_t ermsyt=hyc[khis]->GetRMSError();
    Double_t ermszt=hzc[khis]->GetRMSError();
    Double_t avext=hxc[khis]->GetMean();
    Double_t aveyt=hyc[khis]->GetMean();
    Double_t avezt=hzc[khis]->GetMean();
    Double_t eavext=hxc[khis]->GetMeanError();
    Double_t eaveyt=hyc[khis]->GetMeanError();
    Double_t eavezt=hzc[khis]->GetMeanError();
    Double_t pullxt=hxpullc[khis]->GetRMS();
    Double_t pullyt=hypullc[khis]->GetRMS();
    Double_t pullzt=hzpullc[khis]->GetRMS();
    Double_t epullxt=hxpullc[khis]->GetRMSError();
    Double_t epullyt=hypullc[khis]->GetRMSError();
    Double_t epullzt=hzpullc[khis]->GetRMSError();

    Double_t x=hcont[khis]->GetMean();
    Double_t ex=hcont[khis]->GetRMS();;

    gavexc->SetPoint(khis,x,avext);
    gavexc->SetPointError(khis,ex,eavext);
    gaveyc->SetPoint(khis,x,aveyt);
    gaveyc->SetPointError(khis,ex,eaveyt);
    gavezc->SetPoint(khis,x,avezt);
    gavezc->SetPointError(khis,ex,eavezt);
    grmsxc->SetPoint(khis,x,rmsxt);
    grmsxc->SetPointError(khis,ex,ermsxt);
    grmsyc->SetPoint(khis,x,rmsyt);
    grmsyc->SetPointError(khis,ex,ermsyt);
    grmszc->SetPoint(khis,x,rmszt);
    grmszc->SetPointError(khis,ex,ermszt);

    gpullxc->SetPoint(khis,x,pullxt);
    gpullxc->SetPointError(khis,ex,epullxt);
    gpullyc->SetPoint(khis,x,pullyt);
    gpullyc->SetPointError(khis,ex,epullyt);
    gpullzc->SetPoint(khis,x,pullzt);
    gpullzc->SetPointError(khis,ex,epullzt);

    Float_t nEv=hcont[khis]->GetEntries();
    cout<<"Contrib. bin "<<khis<<"  # Events ="<<nEv<<endl;
  }

  for(Int_t khis=0; khis<nbinz; khis++){

    Double_t rmsxt=hxz[khis]->GetRMS();
    Double_t rmsyt=hyz[khis]->GetRMS();
    Double_t rmszt=hzz[khis]->GetRMS();
    Double_t ermsxt=hxz[khis]->GetRMSError();
    Double_t ermsyt=hyz[khis]->GetRMSError();
    Double_t ermszt=hzz[khis]->GetRMSError();
    Double_t avext=hxz[khis]->GetMean();
    Double_t aveyt=hyz[khis]->GetMean();
    Double_t avezt=hzz[khis]->GetMean();
    Double_t eavext=hxz[khis]->GetMeanError();
    Double_t eaveyt=hyz[khis]->GetMeanError();
    Double_t eavezt=hzz[khis]->GetMeanError();

    Float_t nEv=hz[khis]->GetEntries();
    Double_t x=-999.;
    if(nEv>0) x=hz[khis]->GetMean();
    Double_t ex=hz[khis]->GetRMS();;
    
    gavexz->SetPoint(khis,x,avext);
    gavexz->SetPointError(khis,ex,eavext);
    gaveyz->SetPoint(khis,x,aveyt);
    gaveyz->SetPointError(khis,ex,eaveyt);
    gavezz->SetPoint(khis,x,avezt);
    gavezz->SetPointError(khis,ex,eavezt);
    grmsxz->SetPoint(khis,x,rmsxt);
    grmsxz->SetPointError(khis,ex,ermsxt);
    grmsyz->SetPoint(khis,x,rmsyt);
    grmsyz->SetPointError(khis,ex,ermsyt);
    grmszz->SetPoint(khis,x,rmszt);
    grmszz->SetPointError(khis,ex,ermszt);

    c1z->cd(khis+1);
    hxz[khis]->Draw();
    hxz[khis]->Fit("gaus","Q0");
    fitf= hxz[khis]->GetFunction("gaus");
    Double_t avexg=fitf->GetParameter(1);
    Double_t eavexg=fitf->GetParError(1);
    Double_t rmsxg=fitf->GetParameter(2);
    Double_t ermsxg=fitf->GetParError(2);
    c1z->cd(1*nbinz+khis+1);
    hyz[khis]->Draw();
    hyz[khis]->Fit("gaus","Q0");
    fitf= hyz[khis]->GetFunction("gaus");
    Double_t aveyg=fitf->GetParameter(1);
    Double_t eaveyg=fitf->GetParError(1);
    Double_t rmsyg=fitf->GetParameter(2);
    Double_t ermsyg=fitf->GetParError(2);
    c1z->cd(2*nbinz+khis+1);
    hzz[khis]->Draw();
    hzz[khis]->Fit("gaus","Q0");
    fitf= hzz[khis]->GetFunction("gaus");
    Double_t avezg=fitf->GetParameter(1);
    Double_t eavezg=fitf->GetParError(1);
    Double_t rmszg=fitf->GetParameter(2);
    Double_t ermszg=fitf->GetParError(2);

    gavexzg->SetPoint(khis,x,avexg);
    gavexzg->SetPointError(khis,ex,eavexg);
    gaveyzg->SetPoint(khis,x,aveyg);
    gaveyzg->SetPointError(khis,ex,eaveyg);
    gavezzg->SetPoint(khis,x,avezg);
    gavezzg->SetPointError(khis,ex,eavezg);
    grmsxzg->SetPoint(khis,x,rmsxg);
    grmsxzg->SetPointError(khis,ex,ermsxg);
    grmsyzg->SetPoint(khis,x,rmsyg);
    grmsyzg->SetPointError(khis,ex,ermsyg);
    grmszzg->SetPoint(khis,x,rmszg);
    grmszzg->SetPointError(khis,ex,ermszg);

    Float_t zeff=-999.;
    if(nEv>0) zeff=hzz[khis]->GetEntries()/nEv;
    geffz->SetPoint(khis,x,zeff);
    geffz->SetPointError(khis,ex,0.);
    
    cout<<"Z bin "<<khis<<"  # Events ="<<nEv<<endl;
  }
 
  Double_t efftrk=htot->GetEntries()/totevtriggered;

  printf("EVENTS STATISTICS:\n Total: %d\n Triggered (MB1): %d\n Triggered and with vertex %d\n",totev,totevtriggered,htot->GetEntries());
  if(vtxtype.Contains("SPD")) printf("  %d with Vertexer3D, %d with VertexerZ\n",nvtx3D,nvtxZ);
  printf("Overall efficiency (for triggered evts) Vertexer%s = %f / %f = %f\n",vtxtype.Data(),htot->GetEntries(),totevtriggered,efftrk);

  TFile* in = new TFile("vert-graphs.root","recreate");
  gbeamxz->Write("gbeamxz");
  gbeamyz->Write("gbeamyz");
  grmsxm->Write();
  grmsxmg->Write();
  gpullxm->Write();
  gpullxmg->Write();
  gavexm->Write();
  gavexmg->Write();
  grmsym->Write();
  grmsymg->Write();
  gpullym->Write();
  gpullymg->Write();
  gaveym->Write();
  gaveymg->Write();
  grmszm->Write();
  grmszmg->Write();
  gpullzm->Write();
  gpullzmg->Write();
  gavezm->Write();
  gavezmg->Write();
  geffm->Write();
  gefftrks->Write();
  geff3Dtrks->Write();
  grmsxc->Write();
  gpullxc->Write();
  gavexc->Write();
  grmsyc->Write();
  gpullyc->Write();
  gaveyc->Write();
  grmszc->Write();
  gpullzc->Write();
  gavezc->Write();
  gavexz->Write();
  gaveyz->Write();
  gavezz->Write();
  grmsxz->Write();
  grmsyz->Write();
  grmszz->Write();
  gavexzg->Write();
  gaveyzg->Write();
  gavezzg->Write();
  grmsxzg->Write();
  grmsyzg->Write();
  grmszzg->Write();

  in->Close();

  gStyle->SetOptTitle(0);

  Char_t outgif[100];
  TLine *lin0=new TLine(0,0,60,0);
  lin0->SetLineStyle(2);

  TCanvas *cg1=new TCanvas("cg1","Histo mean");
  cg1->SetBottomMargin(0.14);
  cg1->SetTopMargin(0.08);
  cg1->SetLeftMargin(0.14);
  cg1->SetRightMargin(0.08);
  gavexm->SetMarkerStyle(20);
  gavexm->SetMarkerColor(1);
  gaveym->SetMarkerStyle(21);
  gaveym->SetMarkerColor(2);
  gavezm->SetMarkerStyle(22);
  gavezm->SetMarkerColor(4);
  TLegend *leg=new TLegend(0.18,0.70,0.25,0.90);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  TLegendEntry *ent=leg->AddEntry(gavexm,"x","P");
  ent->SetTextColor(1);
  ent=leg->AddEntry(gaveym,"y","P");
  ent->SetTextColor(2);
  ent=leg->AddEntry(gavezm,"z","P");
  ent->SetTextColor(4);
  gavezm->GetXaxis()->SetLimits(0.,60.);
  gavezm->SetMinimum(-rangeGrAve);
  gavezm->SetMaximum(rangeGrAve);
  gavezm->Draw("AP");
  gavezm->GetXaxis()->SetTitle(trkstitle.Data());
  gavezm->GetXaxis()->SetTitleSize(0.05);
  gavezm->GetYaxis()->SetTitle("<Pos_{found}-Pos_{true}> [#mum]");
  gavezm->GetYaxis()->SetTitleSize(0.05);
  gavexm->Draw("PSAME");
  gaveym->Draw("PSAME");
  leg->Draw();
  lin0->Draw();
  sprintf(outgif,"vert%s-ave-mult.gif",vtxtype.Data());
  if(optgif) cg1->SaveAs(outgif);


  TCanvas *cg2=new TCanvas("cg2","Histo RMS");
  cg2->SetBottomMargin(0.14);
  cg2->SetTopMargin(0.08);
  cg2->SetLeftMargin(0.14);
  cg2->SetRightMargin(0.08);
  grmsxm->SetMarkerStyle(20);
  grmsxm->SetMarkerColor(1);
  grmsym->SetMarkerStyle(21);
  grmsym->SetMarkerColor(2);
  grmszm->SetMarkerStyle(22);
  grmszm->SetMarkerColor(4);
  grmszm->SetMinimum(0);
  grmszm->SetMaximum(rangeGrRms);
  grmszm->GetXaxis()->SetLimits(0,60);
  grmszm->Draw("AP");
  grmszm->GetXaxis()->SetTitle(trkstitle.Data());
  grmszm->GetXaxis()->SetTitleSize(0.05);
  grmszm->GetYaxis()->SetTitle("Resolution [#mum]");
  grmszm->GetYaxis()->SetTitleSize(0.05);
  grmsym->Draw("PSAME");
  grmsxm->Draw("PSAME");
  grmsxmg->SetMarkerStyle(24);
  grmsxmg->SetMarkerColor(1);
  grmsymg->SetMarkerStyle(25);
  grmsymg->SetMarkerColor(2);
  grmszmg->SetMarkerStyle(26);
  grmszmg->SetMarkerColor(4);
  grmszmg->SetMinimum(0);
  grmszmg->SetMaximum(rangeGrRms);
  grmszmg->GetXaxis()->SetLimits(0,60);
  grmszmg->Draw("PSAME");
  grmszmg->GetXaxis()->SetTitle(trkstitle.Data());
  grmszmg->GetXaxis()->SetTitleSize(0.05);
  grmszmg->GetYaxis()->SetTitle("Resolution [#mum]");
  grmszmg->GetYaxis()->SetTitleSize(0.05);
  grmsymg->Draw("PSAME");
  grmsxmg->Draw("PSAME");
  leg->Draw();
  sprintf(outgif,"vert%s-rms-mult.gif",vtxtype.Data());
  if(optgif) cg2->SaveAs(outgif);




  TCanvas *cg3=new TCanvas("cg3","Efficiency vs dNch/dy");
  cg3->SetBottomMargin(0.14);
  cg3->SetTopMargin(0.08);
  cg3->SetLeftMargin(0.14);
  cg3->SetRightMargin(0.08);
  geffm->SetMarkerStyle(22);
  geffm->SetMarkerColor(1);
  geffm->GetXaxis()->SetLimits(0.,40.);
  geffm->SetMinimum(0.);
  geffm->SetMaximum(1.2);
  geffm->Draw("AP");
  geffm->GetXaxis()->SetTitle("MC dN_{ch}/dy in |y|<1");
  geffm->GetXaxis()->SetTitleSize(0.05);
  geffm->GetYaxis()->SetTitle("efficiency");
  geffm->GetYaxis()->SetTitleSize(0.05);
  sprintf(outgif,"vert%s-eff-mult.gif",vtxtype.Data());
  if(optgif) cg3->SaveAs(outgif);


  TCanvas *cg3b=new TCanvas("cg3b","Efficiency vs tracks");
  cg3b->SetBottomMargin(0.14);
  cg3b->SetTopMargin(0.08);
  cg3b->SetLeftMargin(0.14);
  cg3b->SetRightMargin(0.08);
  gefftrks->SetMarkerStyle(22);
  gefftrks->SetMarkerColor(1);
  gefftrks->GetXaxis()->SetLimits(0.,40.);
  gefftrks->SetMinimum(0.);
  gefftrks->SetMaximum(1.2);
  gefftrks->Draw("AP");
  gefftrks->GetXaxis()->SetTitle(trkstitle.Data());
  gefftrks->GetXaxis()->SetTitleSize(0.05);
  gefftrks->GetYaxis()->SetTitle("efficiency");
  gefftrks->GetYaxis()->SetTitleSize(0.05);
  if(vtxtype.Contains("SPD")) {
    geff3Dtrks->SetMarkerStyle(26);
    geff3Dtrks->SetMarkerColor(1);
    geff3Dtrks->Draw("P");
  }
  sprintf(outgif,"vert%s-eff-mult.gif",vtxtype.Data());
  if(optgif) cg3b->SaveAs(outgif);


  TCanvas *cg4=new TCanvas("cg4","Pulls");
  cg4->SetBottomMargin(0.14);
  cg4->SetTopMargin(0.08);
  cg4->SetLeftMargin(0.14);
  cg4->SetRightMargin(0.08);
  gpullxm->SetMarkerStyle(20);
  gpullxm->SetMarkerColor(1);
  gpullym->SetMarkerStyle(21);
  gpullym->SetMarkerColor(2);
  gpullzm->SetMarkerStyle(22);
  gpullzm->SetMarkerColor(4);
  gpullzm->GetXaxis()->SetLimits(0,60);
  gpullzm->SetMinimum(0);
  gpullzm->SetMaximum(rangeGrPull);
  gpullzm->Draw("AP");
  gpullzm->GetXaxis()->SetTitle(trkstitle.Data());
  gpullzm->GetXaxis()->SetTitleSize(0.05);
  gpullzm->GetYaxis()->SetTitle("PULL");
  gpullzm->GetYaxis()->SetTitleSize(0.05);
  gpullxm->Draw("PSAME");
  gpullym->Draw("PSAME");
  gpullxmg->SetMarkerStyle(24);
  gpullxmg->SetMarkerColor(1);
  gpullymg->SetMarkerStyle(25);
  gpullymg->SetMarkerColor(2);
  gpullzmg->SetMarkerStyle(26);
  gpullzmg->SetMarkerColor(4);
  gpullzmg->GetXaxis()->SetLimits(0,60);
  gpullzmg->SetMinimum(0);
  gpullzmg->SetMaximum(rangeGrPull);
  gpullzmg->Draw("PSAME");
  gpullzmg->GetXaxis()->SetTitle(trkstitle.Data());
  gpullzmg->GetXaxis()->SetTitleSize(0.05);
  gpullzmg->GetYaxis()->SetTitle("PULL");
  gpullzmg->GetYaxis()->SetTitleSize(0.05);
  gpullxmg->Draw("PSAME");
  gpullymg->Draw("PSAME");
  TLine *lin=new TLine(0,1,60,1);
  lin->SetLineStyle(2);
  lin->Draw();
  leg->Draw();
  sprintf(outgif,"vert%s-pull-mult.gif",vtxtype.Data());
  if(optgif) cg4->SaveAs(outgif);




  TCanvas *cz1=new TCanvas("cz1","Efficiency vs. Z");
  cz1->SetBottomMargin(0.14);
  cz1->SetTopMargin(0.08);
  cz1->SetLeftMargin(0.14);
  cz1->SetRightMargin(0.08);
  geffz->SetMarkerStyle(22);
  geffz->SetMarkerColor(1);
  geffz->GetXaxis()->SetLimits(-20,20.);
  geffz->SetMinimum(0.);
  geffz->SetMaximum(1.);
  geffz->Draw("AP");
  geffz->GetXaxis()->SetTitle("Z [cm]");
  geffz->GetXaxis()->SetTitleSize(0.05);
  geffz->GetYaxis()->SetTitle("efficiency");
  geffz->GetYaxis()->SetTitleSize(0.05);
  sprintf(outgif,"vert%s-eff-z.gif",vtxtype.Data());
  if(optgif) cz1->SaveAs(outgif);


  TLine *lin0z=new TLine(-20,0,20,0);
  lin0z->SetLineStyle(2);
  TCanvas *cz2=new TCanvas("cz2","Mean vs. Z");
  cz2->SetBottomMargin(0.14);
  cz2->SetTopMargin(0.08);
  cz2->SetLeftMargin(0.14);
  cz2->SetRightMargin(0.08);
  gavexz->SetMarkerStyle(20);
  gavexz->SetMarkerColor(1);
  gaveyz->SetMarkerStyle(21);
  gaveyz->SetMarkerColor(2);
  gavezz->SetMarkerStyle(22);
  gavezz->SetMarkerColor(4);
  gavezz->GetXaxis()->SetLimits(-20,20.);
  gavezz->SetMinimum(-rangeGrAve);
  gavezz->SetMaximum(rangeGrAve);
  gavezz->Draw("AP");
  gavezz->GetXaxis()->SetTitle("Z [cm]");
  gavezz->GetXaxis()->SetTitleSize(0.05);
  gavezz->GetYaxis()->SetTitle("<Pos_{found}-Pos_{true}> [#mum]");
  gavezz->GetYaxis()->SetTitleSize(0.05);
  gavexz->Draw("PSAME");
  gaveyz->Draw("PSAME");
  lin0z->Draw();
  gavexzg->SetMarkerStyle(24);
  gavexzg->SetMarkerColor(1);
  gavexzg->Draw("P");
  gaveyzg->SetMarkerStyle(25);
  gaveyzg->SetMarkerColor(2);
  gaveyzg->Draw("P");
  gavezzg->SetMarkerStyle(26);
  gavezzg->SetMarkerColor(4);
  gavezzg->Draw("P");
  leg->Draw();
  sprintf(outgif,"vert%s-ave-z.gif",vtxtype.Data());
  if(optgif) cz2->SaveAs(outgif);


  TCanvas *cz3=new TCanvas("cz3","Resolution vs. Z");
  cz3->SetBottomMargin(0.14);
  cz3->SetTopMargin(0.08);
  cz3->SetLeftMargin(0.14);
  cz3->SetRightMargin(0.08);
  grmsxz->SetMarkerStyle(20);
  grmsxz->SetMarkerColor(1);
  grmsyz->SetMarkerStyle(21);
  grmsyz->SetMarkerColor(2);
  grmszz->SetMarkerStyle(22);
  grmszz->SetMarkerColor(4);
  grmszz->SetMinimum(0);
  grmszz->SetMaximum(rangeGrRms);
  grmszz->GetXaxis()->SetLimits(-20,20);
  grmszz->Draw("AP");
  grmszz->GetXaxis()->SetTitle("Z [cm]");
  grmszz->GetXaxis()->SetTitleSize(0.05);
  grmszz->GetYaxis()->SetTitle("Resolution [#mum]");
  grmszz->GetYaxis()->SetTitleSize(0.05);
  grmsxz->Draw("PSAME");
  grmsyz->Draw("PSAME");
  leg->Draw();
  sprintf(outgif,"vert%s-rms-z.gif",vtxtype.Data());
  if(optgif) cz3->SaveAs(outgif);

  gStyle->SetPalette(1);
  TCanvas *cbeam = new TCanvas("cbeam","Beam Long",800,800);
  cbeam->Divide(2,2);
  cbeam->cd(1);
  hbeamx->GetYaxis()->SetTitle("X [cm]");
  hbeamx->Draw();
  cbeam->cd(3);
  hbeamx->GetYaxis()->SetTitle("Y [cm]");
  hbeamy->Draw();
  cbeam->cd(2);
  cbeam_2->SetLogz();
  //gbeamxz->SetMarkerStyle(7);
  hbeamxz->Draw("colz");
  hbeamxz->GetXaxis()->SetTitle("Z [cm]");
  hbeamxz->GetYaxis()->SetTitle("X [cm]");
  cbeam_1->Update();
  TPaveStats *st1=(TPaveStats*)hbeamxz->GetListOfFunctions()->FindObject("stats");
  st1->SetX1NDC(0.13);
  st1->SetX2NDC(0.33);
  cbeam_2->Modified();
  cbeam_2->Update();
  cbeam->cd(4);
  cbeam_4->SetLogz();
  //gbeamyz->SetMarkerStyle(7);
  hbeamyz->Draw("colz");
  hbeamyz->GetXaxis()->SetTitle("Z [cm]");
  hbeamyz->GetYaxis()->SetTitle("Y [cm]");
  cbeam_4->Update();
  TPaveStats *st2=(TPaveStats*)hbeamyz->GetListOfFunctions()->FindObject("stats");
  st2->SetX1NDC(0.13);
  st2->SetX2NDC(0.33);
  cbeam_4->Modified();
  cbeam_4->Update();
  cbeam->Update();

  TCanvas *cbeam2 = new TCanvas("cbeam2","Beam Transv",500,500);
  cbeam2->SetLogz();
  cbeam2->SetRightMargin(0.14);
  //gbeamxy->SetMarkerStyle(7);
  hbeamxy->Draw("colz");
  hbeamxy->GetXaxis()->SetTitle("X [cm]");
  hbeamxy->GetYaxis()->SetTitle("Y [cm]");
  cbeam2->Update();
  TPaveStats *st3=(TPaveStats*)hbeamxy->GetListOfFunctions()->FindObject("stats");
  st3->SetX1NDC(0.13);
  st3->SetX2NDC(0.33);
  cbeam2->Modified();
  cbeam2->Update();

  return;
}
//----------------------------------------------------------------------------
void ComputeVtxMean(TString vtxtype="TRK",
		    TString fname="Vertex.Performance.root",
		    TString ntname="fNtupleVertexESD",
		    Int_t nEventsToUse=10000,
		    Int_t mincontr=1) {
  //-----------------------------------------------------------------------
  // Compute weighted mean and cov. matrix from the ntuple
  //-----------------------------------------------------------------------
  gStyle->SetOptStat(0);
  //gStyle->SetOptFit(0);

  Double_t diamondx=0.0200.,diamondy=0.0200.,diamondz=7.5.;

  Double_t avx=0.;
  Double_t avy=0.;
  Double_t avz=0.;
  Double_t wgtavx=0.;
  Double_t wgtavy=0.;
  Double_t wgtavz=0.;
  Double_t sum=0.;
  Double_t sumwgtx=0.;
  Double_t sumwgty=0.;
  Double_t sumwgtz=0.;
  Double_t rmsx=0;
  Double_t rmsy=0;
  Double_t rmsz=0;
  Double_t varx=0.;
  Double_t vary=0.;
  Double_t varz=0.;
  Double_t covxy=0.;
  Double_t covxz=0.;
  Double_t covyz=0.;
  Double_t eavx,eavy,eavz,ewgtavx,ewgtavy,ewgtavz;

  TH1F* hx = new TH1F("hx","",200,-0.1,0.1);
  hx->SetXTitle("vertex x [#mu m]");
  hx->SetYTitle("events");
  TH1F* hy = new TH1F("hy","",200,-0.1,0.1);
  hy->SetXTitle("vertex y [#mu m]");
  hy->SetYTitle("events");
  TH1F* hz = new TH1F("hz","",200,-20,20);
  hz->SetXTitle("vertex z [cm]");
  hz->SetYTitle("events");


  TFile *f=new TFile(fname.Data());
  TList *cOutput = (TList*)f->Get("cOutput");
  TNtuple *nt=(TNtuple*)cOutput->FindObject(ntname.Data());
  Int_t nnnev=nt->GetEntries();
  printf("Events = %d\n",nnnev);
  Float_t xVtx,xdiffVtx,xerrVtx;
  Float_t yVtx,ydiffVtx,yerrVtx;
  Float_t zVtx,zdiffVtx,zerrVtx;
  Float_t ntrklets,ncontrVtx,dndy,triggered,vtx3D;
  Float_t ztrue,zref;
  
  TString sxx="x"; sxx.Append(vtxtype.Data());
  nt->SetBranchAddress(sxx.Data(),&xVtx);
  TString syy="y"; syy.Append(vtxtype.Data());
  nt->SetBranchAddress(syy.Data(),&yVtx);
  TString szz="z"; szz.Append(vtxtype.Data());
  nt->SetBranchAddress(szz.Data(),&zVtx);
  
  TString xerr="xerr"; xerr.Append(vtxtype.Data());
  nt->SetBranchAddress(xerr.Data(),&xerrVtx);
  TString yerr="yerr"; yerr.Append(vtxtype.Data());
  nt->SetBranchAddress(yerr.Data(),&yerrVtx);
  TString zerr="zerr"; zerr.Append(vtxtype.Data());
  nt->SetBranchAddress(zerr.Data(),&zerrVtx);

  TString trkstitle;
  if(vtxtype.Contains("TPC")) {
    nt->SetBranchAddress("nTPCin",&ntrklets);
    trkstitle="TPC tracks pointing to beam pipe";
  } else {
    nt->SetBranchAddress("ntrklets",&ntrklets);
    trkstitle="SPD tracklets";
  }
  TString ntrks="ntrks"; ntrks.Append(vtxtype.Data());
  nt->SetBranchAddress(ntrks.Data(),&ncontrVtx);
  nt->SetBranchAddress("dndygen",&dndy);
  
  nt->SetBranchAddress("ztrue",&ztrue);

  nt->SetBranchAddress("triggered",&triggered);

  nt->SetBranchAddress("SPD3D",&vtx3D);

  Int_t total=0;

  Int_t divider=(Int_t)(nnnev/nEventsToUse);
  // first loop on events
  for(Int_t iev=0;iev<nnnev;iev++) {
    if(iev%divider!=0) continue;
    total++;
    nt->GetEvent(iev);
    if(!vtxtype.Contains("SPD")) vtx3D=1.;
    //if(vtx3D<0.5) continue;
    if(triggered<0.5) continue; // not triggered
    if(ncontrVtx<=0) continue; // no vertex

    if(ncontrVtx<mincontr) continue;

    avx += xVtx;
    avy += yVtx;
    avz += zVtx;
    sum += 1.;
    wgtavx += xVtx/xerrVtx/xerrVtx;
    wgtavy += yVtx/yerrVtx/yerrVtx;
    wgtavz += zVtx/zerrVtx/zerrVtx;
    sumwgtx += 1./xerrVtx/xerrVtx;
    sumwgty += 1./yerrVtx/yerrVtx;
    sumwgtz += 1./zerrVtx/zerrVtx;
     
    hx->Fill(xVtx);
    hy->Fill(yVtx);
    hz->Fill(zVtx);
  }
  
  avx /= sum;
  avy /= sum;
  avz /= sum;
  wgtavx /= sumwgtx;
  wgtavy /= sumwgty;
  wgtavz /= sumwgtz;
  ewgtavx = 1./TMath::Sqrt(sumwgtx);
  ewgtavy = 1./TMath::Sqrt(sumwgty);
  ewgtavz = 1./TMath::Sqrt(sumwgtz);
  

  // second loop on events
  for(Int_t iev=0;iev<nnnev;iev++){
    if(iev%divider!=0) continue;
    nt->GetEvent(iev);
    if(!vtxtype.Contains("SPD")) vtx3D=1.;
    if(vtx3D<0.5) continue;
    if(triggered<0.5) continue; // not triggered
    if(ncontrVtx<=0) continue; // no vertex

    if(ncontrVtx<mincontr) continue;
  
    varx += (xVtx-avx)*(xVtx-avx);
    vary += (yVtx-avy)*(yVtx-avy);
    varz += (zVtx-avz)*(zVtx-avz);
    covxy += (xVtx-avx)*(yVtx-avy);
    covxz += (xVtx-avx)*(zVtx-avz);
    covyz += (yVtx-avy)*(zVtx-avz);
  }
  
  varx /= sum;
  vary /= sum;
  varz /= sum;
  covxy /= sum;
  covxz /= sum;
  covyz /= sum;
  rmsx = TMath::Sqrt(varx);
  rmsy = TMath::Sqrt(vary);
  rmsz = TMath::Sqrt(varz);
  eavx = rmsx/TMath::Sqrt(sum);
  eavy = rmsy/TMath::Sqrt(sum);
  eavz = rmsz/TMath::Sqrt(sum);
  

  printf("\n\nNumber of events: Total %d, Used %d\n",total,sum);
  printf("Minimum number of contributors: %d\n",mincontr);
  printf("Average:\n x = (%f +- %f) cm\n y = (%f +- %f) cm\n z = (%f +- %f) cm\n",avx,eavx,avy,eavy,avz,eavz);
  printf("Weighted Average:\n x = (%f +- %f) cm\n y = (%f +- %f) cm\n z = (%f +- %f) cm\n",wgtavx,ewgtavx,wgtavy,ewgtavy,wgtavz,ewgtavz);
  printf("RMS:\n x = %f cm\n y = %f cm\n z = %f cm\n",rmsx,rmsy,rmsz);

  TCanvas *c = new TCanvas("c","c",0,0,1000,500);
  c->Divide(3,1);
  c->cd(1);
  hx->Draw();
  TF1 *gx = new TF1("gx","gaus",-1000,1000);
  gx->SetLineColor(2);
  hx->Fit(gx,"Q");
  TF1 *gxx = (TF1*)gx->Clone("gxx");
  gxx->FixParameter(2,diamondx);
  gxx->SetLineStyle(2);
  gxx->Draw("same");
  c->cd(2);
  hy->Draw();
  TF1 *gy = new TF1("gy","gaus",-1000,1000);
  gy->SetLineColor(2);
  hy->Fit(gy,"Q");
  TF1 *gyy = (TF1*)gy->Clone("gyy");
  gyy->FixParameter(2,diamondy);
  gyy->SetLineStyle(2);
  gyy->Draw("same");
  c->cd(3);
  hz->Draw();
  TF1 *gz = new TF1("gz","gaus",-10,10);
  gz->SetLineColor(2);
  hz->Fit(gz,"Q");
  TF1 *gzz = (TF1*)gz->Clone("gzz");
  gzz->FixParameter(2,diamondz);
  gzz->SetLineStyle(2);
  gzz->Draw("same");


  return;
}
