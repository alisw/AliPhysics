/*
  
  .x ~/UliStyle.C
  .x ~/NimStyle.C
  .x ~/rootlogon.C
  TProof::Open("");
  gSystem->Load("libSTAT.so");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libTPCcalib");
  gSystem->Load("libSTAT.so");

  gSystem->AddIncludePath("-I$ALICE_ROOT/TPC/macros");  
  gROOT->LoadMacro("$ALICE_ROOT/TPC/macros/AliXRDPROOFtoolkit.cxx+");
  AliXRDPROOFtoolkit tool; 
  TChain * chain = tool.MakeChain("cosmic.txt","Track0",0,1000000);
  chain->Lookup();
  chain->SetProof(kTRUE);

  .L $ALICE_ROOT/TPC/CalibMacros/CosmicPerformance.C+
  chainCosmic=chain;
  MakeCuts()
  MakeAlias();
  Make1PtPlot();
  Draw1Pt();
  Draw1PtPull();

  MakeZPlot();
  DrawZ();
  DrawZPull();

  //
  PtResolPt();
  
  
*/

#include "TTree.h"
#include "TChain.h"
#include "TPad.h"
#include "TCanvas.h"
#include "TCut.h"
#include "TH1.h"
#include "TH2F.h"
#include "AliTPCcalibV0.h"
#include "AliExternalTrackParam.h"

TChain * chainCosmic=0;
Int_t kmicolors[10]={1,2,3,4,6,7,8,9,10,11};
Int_t kmimarkers[10]={21,22,23,24,25,26,27,28,29,30};

AliTPCcalibV0 v0;
Bool_t bUseCorrection =kFALSE;
//
// Global  cuts
//
TCut cutDiff[6];  // parameter diff cuts - 5 all
TCut cutPull[6];  // pull diff cuts      - 5 all
TCut cutGeomTPC;  // TPC geometrical cuts 
TCut cutN("cutN","min(Orig0.fTPCncls,Orig1.fTPCncls)>80");
TCut cutN120("cutN120","min(Orig0.fTPCncls,Orig1.fTPCncls)+abs(Tr1.fP[4]*10)>120");
TCut cutN130("cutN120","min(Orig0.fTPCncls,Orig1.fTPCncls)+abs(Tr1.fP[4]*20)>130");
TCut cutS("cutS","!crossI&&!crossO");
TCut cutRun("run<620600");
TCut cutAll;
TCut cutCustom=cutN130+cutS+cutRun;
TCut cut1Pt = "Tr0.fP[1]>0"; // Use only A side for Pt study
//
// Aliases
//
TString dP[5];      // delta of parameters
TString sP[5];      // sigma of parameters
TString pP[5];      // pull  of parameters
TString axisYDP[5];  // axis title
TString axisYPP[5];  // axis title
TString axisYDPm[5];  // axis title
TString axisYPPm[5];  // axis title
//
//
Float_t range[5] = {5,3,10,3,0.05};      // ranges for diff histogram

Float_t scale[5] = {10,10,1000,1000,1};  // scaling factor
//
// Z depend
TH2F * histoDPZ[5]={0,0,0,0,0};
TH2F * histoSPZ[5]={0,0,0,0,0};
TH2F * histoPPZ[5]={0,0,0,0,0};
TH1F * hmDPZ[5];
TH1F * hmSPZ[5];
TH1F * hmPPZ[5];
TH1F * hsDPZ[5];
TH1F * hsSPZ[5];
TH1F * hsPPZ[5];
//
// Pt depend
TH2F * histoDP1Pt[5]={0,0,0,0,0};
TH2F * histoSP1Pt[5]={0,0,0,0,0};
TH2F * histoPP1Pt[5]={0,0,0,0,0};
TH1F * hmDP1Pt[5];
TH1F * hmSP1Pt[5];
TH1F * hmPP1Pt[5];
TH1F * hsDP1Pt[5];
TH1F * hsSP1Pt[5];
TH1F * hsPP1Pt[5];

void MakeCutsParam(){
  //
  // Default selection criteria
  //
  cutDiff[0]="abs(Tr0.fP[0]+Tr1.fP[0])<3";
  cutDiff[1]="abs(Tr0.fP[1]-Tr1.fP[1])<15"; 
  cutDiff[2]="abs(Tr1.fAlpha-Tr0.fAlpha+pi)<0.1";
  cutDiff[3]="abs(Tr0.fP[3]+Tr1.fP[3])<0.1";     
  cutDiff[4]="abs(Tr0.fP[4]+Tr1.fP[4])<0.5";     
  for (Int_t i=0;i<5;i++) cutDiff[5]+=cutDiff[i];
  //
  cutPull[0]="abs(Tr0.fP[0]+Tr1.fP[0])/sqrt(Tr0.fC[0]+Tr1.fC[0])<10";
  cutPull[1]="1"; 
  cutPull[2]="((Tr1.fAlpha-Tr0.fAlpha+pi))/sqrt(Tr0.fC[5]+Tr1.fC[5])<10";
  cutPull[3]="1";     
  cutPull[4]="abs(Tr0.fP[4]+Tr1.fP[4])/sqrt(Tr0.fC[14]+Tr1.fC[14])<10";     
  for (Int_t i=0;i<5;i++) cutPull[5]+=cutPull[i];
}

void MakeGeomCuts(){
//
// Geometrical cut 
//
  TCut cutOx("Op1.fX>240&&Op0.fX>240");
  TCut cutOz("abs(Op1.fP[1])<240&&abs(Op0.fP[1])<240");
  TCut cutIz("abs(Ip1.fP[1])<240&&abs(Ip0.fP[1])<240");
  TCut cutX00("abs(x00)<70");
  TCut cutX10("abs(x10)<70");
  TCut cutT1P2("abs(Ip1.fP[2])<0.8");
  TCut cutT0P2("abs(Ip0.fP[2])<0.8");
  cutGeomTPC = cutOx+cutOz+cutIz+cutX00+cutX10+cutT1P2+cutT0P2;
}

void MakeCuts(){
  // make cuts all 
  MakeGeomCuts();
  MakeCutsParam();
  cutAll = cutDiff[5]+cutPull[5]+cutGeomTPC;
}

void MakeAlias(){
  AliExternalTrackParam p;
  dP[0]="(Tr0.fP[0]+Tr1.fP[0])";
  dP[1]="(Tr0.fP[1]-Tr1.fP[1])"; 
  dP[2]="(Tr1.fAlpha-Tr0.fAlpha+pi)";
  dP[3]="(Tr0.fP[3]+Tr1.fP[3])";     
  dP[4]="(Tr0.fP[4]+Tr1.fP[4])";     
  for (Int_t i=0;i<5;i++){
    sP[i]=Form("%f*sqrt(Tr0.fC[%d]+Tr1.fC[%d])",scale[i],p.GetIndex(i,i),p.GetIndex(i,i));
    pP[i]=Form("%s/sqrt(Tr0.fC[%d]+Tr1.fC[%d])",dP[i].Data(),scale[i],p.GetIndex(i,i),p.GetIndex(i,i));
    dP[i]+=Form("*%f",scale[i]);
  }
 
  axisYDP[0]="#sigma_{r#phi} (mm)";
  axisYDP[1]="#sigma_{z} (mm)";
  axisYDP[2]="#sigma_{#phi} (mrad)";
  axisYDP[3]="#sigma_{#theta} (mrad)";
  axisYDP[4]="#sigma_{1/pt} (1/GeV))";
  //
  axisYPP[0]="#sigma_{r#phi} (Unit)";
  axisYPP[1]="#sigma_{z} (Unit)";
  axisYPP[2]="#sigma_{#phi} (Unit)";
  axisYPP[3]="#sigma_{#theta} (Unit)";
  axisYPP[4]="#sigma_{1/pt} (Unit))";
  //
  axisYDPm[0]="#Delta_{r#phi} (mm)";
  axisYDPm[1]="#Delta_{z} (mm)";
  axisYDPm[2]="#Delta_{#phi} (mrad)";
  axisYDPm[3]="#Delta_{#theta} (mrad)";
  axisYDPm[4]="#Delta_{1/pt} (1/GeV))";
  //
  axisYPPm[0]="#Delta_{r#phi} (Unit)";
  axisYPPm[1]="#Delta_{z} (Unit)";
  axisYPPm[2]="#Delta_{#phi} (Unit)";
  axisYPPm[3]="#Delta_{#theta} (Unit)";
  axisYPPm[4]="#Delta_{1/pt} (Unit))";
}

void MakeZPlot(){
  
  for (Int_t i=0;i<5;i++){
    char hname[100];
    sprintf(hname,"dP%ivZ",i);
    //
    if ( histoDPZ[i]==0){
      histoDPZ[i]= new TH2F(hname, hname,10,-240,240,100,-range[i],range[i]);
      sprintf(hname,"sP%ivZ",i);
      histoSPZ[i]= new TH2F(hname, hname,10,-240,240,100,0,range[i]/3);
      sprintf(hname,"pP%ivZ",i);
      histoPPZ[i]= new TH2F(hname, hname,10,-240,240,100,-6,6);
    }
    //
    histoDPZ[i]->SetXTitle("z (cm)");
    histoSPZ[i]->SetXTitle("z (cm)");
    histoPPZ[i]->SetXTitle("z (cm)");
    histoDPZ[i]->SetYTitle(axisYDP[i]);
    histoSPZ[i]->SetYTitle(axisYDP[i]);
    histoPPZ[i]->SetYTitle(axisYPP[i]);
    chainCosmic->Draw((dP[i]+"/sqrt(2.):Tr0.fP[1]>>"+histoDPZ[i]->GetName()),cutAll+cutCustom+cutS);
    chainCosmic->Draw((sP[i]+"/sqrt(2.):Tr0.fP[1]>>"+histoSPZ[i]->GetName()),cutAll+cutCustom+cutS); 
    chainCosmic->Draw((pP[i]+":Tr0.fP[1]>>"+histoPPZ[i]->GetName()),cutAll+cutCustom+cutS); 
  }
  
  TObjArray array(3);
  {
    for (Int_t i=0;i<5;i++){
      histoDPZ[i]->FitSlicesY(0,0,-1,0,"QNR",&array);
      hmDPZ[i] = (TH1F*)((array.At(1))->Clone());
      hsDPZ[i] = (TH1F*)((array.At(2))->Clone());
      histoSPZ[i]->FitSlicesY(0,0,-1,0,"QNR",&array);
      hmSPZ[i] = (TH1F*)((array.At(1))->Clone());
      hsSPZ[i] = (TH1F*)((array.At(2))->Clone());
      histoPPZ[i]->FitSlicesY(0,0,-1,0,"QNR",&array);
      hmPPZ[i] = (TH1F*)((array.At(1))->Clone());
      hsPPZ[i] = (TH1F*)((array.At(2))->Clone());
      //
      hmDPZ[i]->SetYTitle(axisYDPm[i]);
      hmSPZ[i]->SetYTitle(axisYDPm[i]);
      hmPPZ[i]->SetYTitle(axisYPPm[i]);
      hsDPZ[i]->SetMinimum(0);
      hsDPZ[i]->SetYTitle(axisYDP[i]);
      hsSPZ[i]->SetYTitle(axisYDP[i]);
      hsPPZ[i]->SetYTitle(axisYPP[i]);
      //
      hmDPZ[i]->SetMarkerColor(kmicolors[1]);
      hmDPZ[i]->SetMarkerStyle(kmimarkers[1]);
      hsDPZ[i]->SetMarkerColor(kmicolors[2]);
      hsDPZ[i]->SetMarkerStyle(kmimarkers[2]);
      hmPPZ[i]->SetMarkerColor(kmicolors[1]);
      hmPPZ[i]->SetMarkerStyle(kmimarkers[1]);
      hsPPZ[i]->SetMarkerColor(kmicolors[2]);
      hsPPZ[i]->SetMarkerStyle(kmimarkers[2]);
    }
  }
}

void Make1PtPlot(){
  
  for (Int_t i=0;i<5;i++){
    char hname[100];
    sprintf(hname,"dP%iv1Pt",i);
    //
    if ( histoDP1Pt[i]==0){
      histoDP1Pt[i]= new TH2F(hname, hname,6,0.02,0.75,100,-range[i],range[i]);
      sprintf(hname,"sP%ivZ",i);
      histoSP1Pt[i]= new TH2F(hname, hname,6,0.02,0.75,100,0,range[i]/3);
      sprintf(hname,"pP%ivZ",i);
      histoPP1Pt[i]= new TH2F(hname, hname,6,0.02,0.75,100,-6,6);
    }
    //
    histoDP1Pt[i]->SetXTitle("#sqrt{1/p_{t} (GeV)}");
    histoSP1Pt[i]->SetXTitle("#sqrt{1/p_{t} (GeV)}");
    histoPP1Pt[i]->SetXTitle("#sqrt{1/p_{t} (GeV)}");
    histoDP1Pt[i]->SetYTitle(axisYDP[i]);
    histoSP1Pt[i]->SetYTitle(axisYDP[i]);
    histoPP1Pt[i]->SetYTitle(axisYPP[i]);
    chainCosmic->Draw((dP[i]+"/sqrt(2.):sqrt(0.5*abs(Tr0.fP[4]-Tr1.fP[4]))>>"+histoDP1Pt[i]->GetName()),cutAll+cutCustom+cutS+cut1Pt);
    chainCosmic->Draw((sP[i]+"/sqrt(2.):sqrt(0.5*abs(Tr0.fP[4]-Tr1.fP[4]))>>"+histoSP1Pt[i]->GetName()),cutAll+cutCustom+cutS+cut1Pt);
    chainCosmic->Draw((pP[i]+":sqrt(0.5*abs(Tr0.fP[4]-Tr1.fP[4]))>>"+histoPP1Pt[i]->GetName()),cutAll+cutCustom+cutS+cut1Pt);    
  }
  
  TObjArray array(3);
  {
    for (Int_t i=0;i<5;i++){
      histoDP1Pt[i]->FitSlicesY(0,0,-1,0,"QNR",&array);
      hmDP1Pt[i] = (TH1F*)((array.At(1))->Clone());
      hsDP1Pt[i] = (TH1F*)((array.At(2))->Clone());
      histoSP1Pt[i]->FitSlicesY(0,0,-1,0,"QNR",&array);
      hmSP1Pt[i] = (TH1F*)((array.At(1))->Clone());
      hsSP1Pt[i] = (TH1F*)((array.At(2))->Clone());
      histoPP1Pt[i]->FitSlicesY(0,0,-1,0,"QNR",&array);
      hmPP1Pt[i] = (TH1F*)((array.At(1))->Clone());
      hsPP1Pt[i] = (TH1F*)((array.At(2))->Clone());
      hmDP1Pt[i]->SetYTitle(axisYDPm[i]);
      hmSP1Pt[i]->SetYTitle(axisYDPm[i]);
      hmPP1Pt[i]->SetYTitle(axisYPPm[i]);
      hsDP1Pt[i]->SetMinimum(0);
      hsDP1Pt[i]->SetYTitle(axisYDP[i]);
      hsSP1Pt[i]->SetYTitle(axisYDP[i]);
      hsPP1Pt[i]->SetYTitle(axisYPP[i]);
      //
      hmDP1Pt[i]->SetMarkerColor(kmicolors[1]);
      hmDP1Pt[i]->SetMarkerStyle(kmimarkers[1]);
      hsDP1Pt[i]->SetMarkerColor(kmicolors[2]);
      hsDP1Pt[i]->SetMarkerStyle(kmimarkers[2]);
      hmPP1Pt[i]->SetMarkerColor(kmicolors[1]);
      hmPP1Pt[i]->SetMarkerStyle(kmimarkers[1]);
      hsPP1Pt[i]->SetMarkerColor(kmicolors[2]);
      hsPP1Pt[i]->SetMarkerStyle(kmimarkers[2]);


    }
  }
}

void DrawZ(){
  TCanvas *czd = new TCanvas("Z depend (abs)","Z depend (abs)",700,900);
  czd->Divide(2,5);
  for (Int_t i=0;i<5;i++){
    czd->cd(2*i+1);
    hmDPZ[i]->Draw("");
    czd->cd(2*i+2);
    hsDPZ[i]->Draw("");
  }
  czd->SaveAs("picResol/deltaPxZ.eps");
  czd->SaveAs("picResol/deltaPxZ.gif");
  czd->SaveAs("picResol/deltaPxZ.root");
}

void DrawZPull(){
  TCanvas *czp = new TCanvas("Z depend (Pull)","Z depend (Pull)",700,900);
  czp->Divide(2,5);
  for (Int_t i=0;i<5;i++){
    czp->cd(2*i+1);
    hmPPZ[i]->Draw("");
    czp->cd(2*i+2);
    hsPPZ[i]->Draw("");
  }
  czp->SaveAs("picResol/pullPxZ.eps");
  czp->SaveAs("picResol/pullPxZ.gif");
  czp->SaveAs("picResol/pullPxZ.root");

}


void Draw1Pt(){
  TCanvas *cpd = new TCanvas("1/Pt depend","1/Pt depend",700,900);
  cpd->Divide(2,5);
  for (Int_t i=0;i<5;i++){
    cpd->cd(2*i+1);
    hmDP1Pt[i]->Draw("");
    cpd->cd(2*i+2);
    hsDP1Pt[i]->Draw("");
  }
  cpd->SaveAs("picResol/deltaPx1Pt.eps");
  cpd->SaveAs("picResol/deltaPx1Pt.gif");
  cpd->SaveAs("picResol/deltaPx1Pt.root");
}
void Draw1PtPull(){
  TCanvas *cpp = new TCanvas("Pull 1/Pt","Pull 1/Pt",700,900);
  cpp->Divide(2,5);
  for (Int_t i=0;i<5;i++){
    cpp->cd(2*i+1);
    hmPP1Pt[i]->Draw("");
    cpp->cd(2*i+2);
    hsPP1Pt[i]->Draw("");
  }
  cpp->SaveAs("picResol/pullPx1Pt.eps");
  cpp->SaveAs("picResol/pullPx1Pt.gif");
  cpp->SaveAs("picResol/pullPx1Pt.root");
}


/*
//
//
//


void DrawPtSpectra(){
  TH1F * hisPt0 = new TH1F("hisPt0","hisPt0",50,0,100);
  TH1F * hisPtC = new TH1F("hisPtC","hisPtC",50,0,100);
  chainCosmic->Draw("Tr0.Pt()>>hisPt0",cutAll+cutCustom);
  chainCosmic->Draw("Tr0.Pt()>>hisPtC","abs(Tr0.fP[4])>3*sqrt(Tr0.fC[14])"+cutAll+cutCustom);
  //
  hisPt0->SetXTitle("p_{t} (GeV)");
  hisPt0->SetLineColor(kmicolors[1]);
  hisPtC->SetLineColor(kmicolors[2]);
  //
  hisPtC->Fit("exp");
  hisPt0->Draw("");
  hisPtC->Draw("same"); 
  TLegend * legend = new TLegend(.4,.7, .99, .99,"Cosmic p_{t} spectra");
  legend->AddEntry(hisPt0,"Raw spectra");
  legend->AddEntry(hisPtC,"Selection abs(p_{t})<3#sigma_{p_{t}}");
  legend->Draw();
  gPad->SaveAs("picSpectra/ptSpectra.eps");
  gPad->SaveAs("picSpectra/ptSpectra.gif");
  gPad->SaveAs("picSpectra/ptSpectra.root");
}



void InitCuts(){
  //
  // Init cuts
  //
  chainCosmic->Draw(">>listELP",cutAll,"entryList");
  TEntryList *elist = (TEntryList*)gDirectory->Get("listELP");
  chainCosmic->SetEntryList(elist);
  //
  chainCosmic->Draw(">>listELFit",cutAll+cuthpt+cutS+cutRun,"entryList");
  TEntryList *elistFit = (TEntryList*)gDirectory->Get("listELFit");
  chainCosmic->SetEntryList(elistFit);
}

void SetAlias(){
  //
  // Set aliases
  //
  chainCosmic->SetAlias("dP0","(Tr0.fP[0]+Tr1.fP[0])");
  chainCosmic->SetAlias("dP1","(Tr0.fP[1]-Tr1.fP[1])");
  chainCosmic->SetAlias("dP2","(Tr1.fAlpha-Tr0.fAlpha+pi)");
  chainCosmic->SetAlias("dP3","(Tr0.fP[3]+Tr1.fP[3])");
  chainCosmic->SetAlias("dP4","(Tr0.fP[4]+Tr1.fP[4])");
  //
  chainCosmic->SetAlias("sP0","sqrt(Tr0.fC[0]+Tr1.fC[0])");
  chainCosmic->SetAlias("sP1","sqrt(Tr0.fC[2]+Tr1.fC[2])");
  chainCosmic->SetAlias("sP2","sqrt(Tr0.fC[5]+Tr0.fC[5])");
  chainCosmic->SetAlias("sP3","sqrt(Tr0.fC[9]+Tr1.fC[9])");
  chainCosmic->SetAlias("sP4","sqrt(Tr0.fC[14]+Tr1.fC[14])");
  //
  chainCosmic->SetAlias("corrP0","0");
  chainCosmic->SetAlias("corrP1","0");
  chainCosmic->SetAlias("corrP2","0");
  chainCosmic->SetAlias("corrP3","0");
  chainCosmic->SetAlias("corrP4","0");
  //
  chainCosmic->SetAlias("dR","(1-abs(Tr0.fP[1]/250))");
  chainCosmic->SetAlias("side","(-1+(Tr0.fP[1]>0)*2)");
  chainCosmic->SetAlias("meanPt","((Tr0.Pt()+Tr1.Pt())/2.)");

}

void Correction(){
  //
  // Fit corrections
  //
  TStatToolkit toolkit;
  Double_t chi2=0;
  Int_t    npoints=0;
  TVectorD fitParam;
  TMatrixD covMatrix;
  //
  TString fstring="";
  fstring+="side++";
  //
  fstring+="dR++";
  fstring+="dR*dR++";
  fstring+="Tr0.fP[3]++";
  //
  fstring+="dR*side++";
  fstring+="dR*dR*side++";
  fstring+="Tr0.fP[3]*side++";
  //
  TString * strP0 = TStatToolkit::FitPlane(chainCosmic,"dP0", fstring.Data(),"1", chi2,npoints,fitParam,covMatrix);
  chainCosmic->SetAlias("corrP0",strP0->Data());

  TString * strP1 = TStatToolkit::FitPlane(chainCosmic,"dP1", fstring.Data(),"!crossI&&!crossO", chi2,npoints,fitParam,covMatrix);
  chainCosmic->SetAlias("corrP1",strP1->Data());

  TString * strP2 = TStatToolkit::FitPlane(chainCosmic,"dP2", fstring.Data(),"1", chi2,npoints,fitParam,covMatrix);
  chainCosmic->SetAlias("corrP2",strP2->Data());

  TString * strP3 = TStatToolkit::FitPlane(chainCosmic,"dP3", fstring.Data(),"!crossI&&!crossO", chi2,npoints,fitParam,covMatrix);
  chainCosmic->SetAlias("corrP3",strP3->Data());

  TString * strP4 = TStatToolkit::FitPlane(chainCosmic,"dP4", fstring.Data(),"1", chi2,npoints,fitParam,covMatrix);
  chainCosmic->SetAlias("corrP4",strP4->Data());
  bUseCorrection = kTRUE;
}



void DrawNClusterGeom(){
  //
  //
  //
  TH2F * hNd("hNd","hNd",50,0,190,100,10,160);
  chainCosmic->Draw("Orig0.fTPCncls:sqrt(x00^2+x01^2)>>hNd",cutOx+cutOz+cutIx+cutIz);
  hNd->FitSlicesY();
  hNd_1->SetXTitle("DCA_{r} (cm)");
  hNd_1->SetYTitle("Mean number of clusters");
  gPad->SaveAs("pic/NCl_Radius.eps");
  gPad->SaveAs("pic/NCl_Radius.gif");
  gPad->SaveAs("pic/NCl_Radius.root");


}


void PtResolPt(){

  TH2F * hdPtPt = new TH2F("hdPtPt","hdPtPt",20,0.5,30,100,-60,60);
  TH2F * hsPtPt = new TH2F("hsPtPt","hsPtPt",20,0.5,30,200,-0,30);
  TH2F * hdPtPtNoCor = new TH2F("hdPtPtNoCor","hdPtPtNoCor",20,0.5,30,100,-60,60);
  TH2F * hdPtPtCor = new TH2F("hdPtPtCor","hdPtPtCor",20,0.5,30,200,-60,60);
  //
  v0.BinLogX(hdPtPt);
  v0.BinLogX(hsPtPt);
  v0.BinLogX(hdPtPtNoCor);
  v0.BinLogX(hdPtPtCor);

  chainCosmic->Draw("100*((Tr0.Pt()-Tr1.Pt())/meanPt)/sqrt(2.):meanPt>>hdPtPt","side>0"+cutAll+cutRun+cutS+cutN120,"");
  chainCosmic->Draw("100*(sP4*meanPt)/sqrt(2.):meanPt>>hsPtPt","side>0"+cutAll+cutRun+cutS+cutN120,"");
  if (bUseCorrection) chainCosmic->Draw("100*((dP4-corrP4)*meanPt)/sqrt(2.):meanPt>>hdPtPtCorr","side>0"+cutAll+cutRun+cutS+cutN120,"");
  if (bUseCorrection) chainCosmic->Draw("100*(dP4*meanPt)/sqrt(2.):meanPt>>hdPtPtNoCorr","side>0"+cutAll+cutRun+cutS+cutN120,"");
  hdPtPt->FitSlicesY();
  hsPtPt->FitSlicesY();
  hdPtPt_2->Fit("pol1");
  {if (bUseCorrection){
     hdPtPtNoCor->FitSlicesY();
     hdPtPtCor->FitSlicesY();
  }}

  hdPtPt_2->SetXTitle("p_{t} (GeV)");
  hdPtPt_2->SetYTitle("#sigma_{p_{t}}/p_{t} (%)");
  hdPtPt_2->SetMinimum(0.5);
  hsPtPt_1->SetMinimum(0.5);
  hdPtPt_2->SetLineColor(kmicolors[1]);
  hdPtPt_2->SetMarkerStyle(kmimarkers[1]);

  hdPtPt_2->Draw();
  {if (bUseCorrection){
     hdPtPtNoCor_2->Draw("same");
     hdPtPtCor_2->Draw("same");
  }}
  gPad->SaveAs("picPerformance/SigmaPt_pt.gif");
  gPad->SaveAs("picPerformance/SigmaPt_pt.eps");
  gPad->SaveAs("picPerformance/SigmaPt_pt.root");
  hsPtPt_1->SetLineColor(kmicolors[2]);
  hsPtPt_1->SetMarkerStyle(kmimarkers[2]);
  hsPtPt_1->Draw("same");
  gPad->SaveAs("picPerformance/SigmaPt_ptLimit.gif");
  gPad->SaveAs("picPerformance/SigmaPt_ptLimit.eps");
  gPad->SaveAs("picPerformance/SigmaPt_ptLimit.root");

}





void PtResolN(){
  //
  //
  //
  TH2F * hdP4Ncl= new TH2F("hdp4Ncl","hdp4Ncl",5,80,160,100,-0.1,0.1);
  chainCosmic->Draw("(dP4-corrP4)/sqrt(2.):min(Orig0.fTPCncls,Orig1.fTPCncls)>>hdp4Ncl","side>0"+cuthpt+cutRun+cutS,"");
  hdp4Ncl->FitSlicesY();
  hdp4Ncl_2->SetXTitle("Number of clusters");
  hdp4Ncl_2->SetYTitle("#sigma 1/p_{t} (1/GeV)");
  hdp4Ncl_2->Draw();
  gPad->SaveAs("pic/SigmaP4_N.gif");
  gPad->SaveAs("pic/SigmaP4_N.eps");
  gPad->SaveAs("pic/SigmaP4_N.root");

  //
  //
  TH2F * hdP4PullNcl = new TH2F("hdP4PullNcl","hdP4PullNcl",5,80,160,100,-6.1,6.1);
  chainCosmic->Draw("(Tr1.fP[4]+Tr0.fP[4])/sqrt(Tr1.fC[14]+Tr0.fC[14]):min(Orig0.fTPCncls,Orig1.fTPCncls)>>hdP4PullNcl","side>0"+cuthpt+cutRun+cutS,"");
  hdP4PullNcl->FitSlicesY();
  hdP4PullNcl_2->SetXTitle("Number of clusters");
  hdP4PullNcl_2->SetYTitle("#sigma 1/p_{t} (Unit)");
  hdP4PullNcl_2->Draw();
  gPad->SaveAs("pic/PullP4_N.gif");
  gPad->SaveAs("pic/PullP4_N.eps");
  gPad->SaveAs("pic/PullP4_N.root");

}



void dEdxRatio(){
  TH2F hratioPt("hratioPt","hratioPt",10,0,10,100,0.6,1.4);
  chainCosmic->Draw("Seed0.CookdEdxNorm(0.01,0.65,1,0,159,0,kFALSE)/(0.3091*Seed0.CookdEdxNorm(0.01,0.65,0,0,159,0,kFALSE)):sqrt(Tr0.P())>>hratioPt","min(Orig0.fTPCncls,Orig1.fTPCncls)>60&&abs(mag)>0.1&&abs(Tr0.fP[3])>0.03","",10000000);
  //
  //
  //
  hratioPt->FitSlicesY();
  TH2F hratioP3("hratioP3","hratioP3",10,0.0,0.6,100,0.6,1.4);
  chainCosmic->Draw("Seed0.CookdEdxNorm(0.01,0.65,1,0,159,0,kFALSE)/(0.3091*Seed0.CookdEdxNorm(0.01,0.65,0,0,159,0,kFALSE)):abs(Tr0.fP[3])>>hratioP3","min(Orig0.fTPCncls,Orig1.fTPCncls)>60&&abs(mag)>0.1&&abs(Tr0.fP[3])>0.03","",10000000);
}




///////////////////////////////////////////////////////////////////////////
//
//
//
// RESOLUTION  as function of 1/pt
//
//
//
///////////////////////////////////////////////////////////////////////////

void P0resol1Pt(){
  //
  // P0 - Y -DCA resolution as function of the 1/pt
  //
  TH2F * hdP01Pt = new TH2F("hdP01Pt","hdP01Pt",6,0,1.0,100,-1.05,1.05);
  TH2F * hdP01PtNoCor = new TH2F("hdP01PtNoCor","hdP01PtNoCor",6,0,1.0,100,-1.05,1.05);
  chainCosmic->Draw("(dP0-corrP0)/sqrt(2.):sqrt(1/meanPt)>>hdP01Pt","side>0"+cutAll+cutRun+cutS+cutN120,"");
  chainCosmic->Draw("(dP0-0)/sqrt(2.):sqrt(1/meanPt)>>hdP01PtNoCor","side>0"+cutAll+cutRun+cutS+cutN120,"");
  hdP01Pt->FitSlicesY();
  hdP01PtNoCor->FitSlicesY();
  TH2F * hsP01Pt = new TH2F("hsP01Pt","hsP01Pt",6,0,1.0,100,-0.05,0.5);
  chainCosmic->Draw("(sP0)/sqrt(2.):sqrt(1/meanPt)>>hsP01Pt","side>0"+cutAll+cutRun+cutS+cutN120,"");
  TH1 *hsP01Pt_2 = hsP01Pt->ProfileX();

  hdP01Pt_2->SetMinimum(0);
  hdP01Pt_2->SetXTitle("#sqrt{1/p_{t} (1/GeV)}");
  hdP01Pt_2->SetYTitle("#sigma_{y} (cm)");
  hdP01Pt_2->SetMarkerStyle(kmimarkers[0]);  
  hsP01Pt_2->SetMarkerStyle(kmimarkers[2]);  
  hdP01PtNoCor_2->SetMarkerStyle(kmimarkers[1]);
  hdP01Pt_2->Draw();  
  hsP01Pt_2->Draw("same");  
  hdP01PtNoCor_2->Draw("same");

  gPad->SaveAs("picPerformance/SigmaP0_1pt.gif");
  gPad->SaveAs("picPerformance/SigmaP0_1pt.eps");
  gPad->SaveAs("picPerformance/SigmaP0_1pt.root");
  //
  TH2F * hPullP01Pt = new TH2F("hhPullP01Pt","hhPullP01Pt",6,0,1,50,-5.05,5.05);
  TH2F * hncPullP01Pt = new TH2F("hhncPullP01Pt","hhncPullP01Pt",6,0,1,50,-5.05,5.05);
  chainCosmic->Draw("(dP0-corrP0)/sP0:sqrt(1/meanPt)>>hhPullP01Pt","side>0"+cutAll+cutRun+cutS+cutN120,"");
  chainCosmic->Draw("(dP0-0)/sP0:sqrt(1/meanPt)>>hhncPullP01Pt","side>0"+cutAll+cutRun+cutS+cutN120,"");
  hPullP01Pt->FitSlicesY();
  hncPullP01Pt->FitSlicesY();
  hhPullP01Pt_2->SetMarkerStyle(kmimarkers[0]);
  hhncPullP01Pt_2->SetMarkerStyle(kmimarkers[1]);
  hhPullP01Pt_2->SetMinimum(0);
  hhPullP01Pt_2->SetXTitle("#sqrt{1/p_{t} (1/GeV)}");
  hhPullP01Pt_2->SetYTitle("#sigma_{y}(Unit) ");
  hhPullP01Pt_2->Draw();
  //hhncPullP01Pt_2->Draw("same");
  gPad->SaveAs("picPerformance/PullP0_1pt.gif");
  gPad->SaveAs("picPerformance/PullP0_1pt.eps");
  gPad->SaveAs("picPerformance/PullP0_1pt.root");
}



void P1resol1Pt(){
  //
  // P1 - Z -DCA resolution as function of the 1/pt
  //
  TH2F * hdP11Pt = new TH2F("hdP11Pt","hdP11Pt",6,0,1.0,100,-1.05,1.05);
  TH2F * hdP11PtNoCor = new TH2F("hdP11PtNoCor","hdP11PtNoCor",6,0,1.0,100,-1.05,1.05);
  chainCosmic->Draw("(dP1-corrP1)/sqrt(2.):sqrt(1/meanPt)>>hdP11Pt","side>0"+cutAll+cutRun+cutS+cutN120,"");
  chainCosmic->Draw("(dP1-0)/sqrt(2.):sqrt(1/meanPt)>>hdP11PtNoCor","side>0"+cutAll+cutRun+cutS+cutN120,"");
  hdP11Pt->FitSlicesY();
  hdP11PtNoCor->FitSlicesY();
  hdP11Pt_2->SetMinimum(0);
  hdP11Pt_2->SetXTitle("#sqrt{1/p_{t} (1/GeV)}");
  hdP11Pt_2->SetYTitle("#sigma_{z} (cm)");
  hdP11Pt_2->SetMarkerStyle(kmimarkers[0]);  
  hdP11PtNoCor_2->SetMarkerStyle(kmimarkers[1]);
  hdP11Pt_2NoCor->Draw();  
  //  hdP11PtNoCor_2->Draw("same");
  gPad->SaveAs("picPerformance/SigmaP1_1pt.gif");
  gPad->SaveAs("picPerformance/SigmaP1_1pt.eps");
  gPad->SaveAs("picPerformance/SigmaP1_1pt.root");

  //
  TH2F * hPullP11Pt = new TH2F("hhPullP11Pt","hhPullP11Pt",6,0,1,50,-5.05,5.05);
  TH2F * hncPullP11Pt = new TH2F("hhncPullP11Pt","hhncPullP11Pt",6,0,1,50,-5.05,5.05);
  chainCosmic->Draw("(dP1-corrP1)/sP1:sqrt(1/meanPt)>>hhPullP11Pt","side>0"+cutAll+cutRun+cutS+cutN120,"");
  chainCosmic->Draw("(dP1-0)/sP1:sqrt(1/meanPt)>>hhncPullP11Pt","side>0"+cutAll+cutRun+cutS+cutN120,"");
  hPullP11Pt->FitSlicesY();
  hncPullP11Pt->FitSlicesY();
  hhPullP11Pt_2->SetMarkerStyle(kmimarkers[0]);
  hhncPullP11Pt_2->SetMarkerStyle(kmimarkers[1]);
  hhPullP11Pt_2->SetMinimum(0);
  hhPullP11Pt_2->SetXTitle("#sqrt{1/p_{t} (1/GeV)}");
  hhPullP11Pt_2->SetYTitle("#sigma_{z}(Unit) ");
  hhPullP11Pt_2NoCor->Draw();
  //hhncPullP11Pt_2->Draw("same");
  gPad->SaveAs("picPerformance/PullP1_1pt.gif");
  gPad->SaveAs("picPerformance/PullP1_1pt.eps");
  gPad->SaveAs("picPerformance/PullP1_1pt.root");
}

void P2resol1Pt(){
  //
  // P2 - Z -DCA resolution as function of the 1/pt
  //
  TH2F * hdP21Pt = new TH2F("hdP21Pt","hdP21Pt",6,0,1.0,100,-20.05,20.05);
  TH2F * hdP21PtNoCor = new TH2F("hdP21PtNoCor","hdP21PtNoCor",6,0,1.0,100,-20.05,20.05);
  chainCosmic->Draw("1000*(dP2-corrP2)/sqrt(2.):sqrt(1/meanPt)>>hdP21Pt","side>0"+cutAll+cutRun+cutS+cutN120,"");
  chainCosmic->Draw("1000*(dP2-0)/sqrt(2.):sqrt(1/meanPt)>>hdP21PtNoCor","side>0"+cutAll+cutRun+cutS+cutN120,"");
  hdP21Pt->FitSlicesY();
  hdP21PtNoCor->FitSlicesY();
  hdP21Pt_2->SetMinimum(0);
  hdP21Pt_2->SetXTitle("#sqrt{1/p_{t} (1/GeV)}");
  hdP21Pt_2->SetYTitle("#sigma_{#Phi} (mrad)");
  hdP21Pt_2->SetMarkerStyle(kmimarkers[0]);  
  hdP21PtNoCor_2->SetMarkerStyle(kmimarkers[1]);
  hdP21Pt_2NoCor->Draw();  
  //  hdP21PtNoCor_2->Draw("same");
  gPad->SaveAs("picPerformance/SigmaP2_1pt.gif");
  gPad->SaveAs("picPerformance/SigmaP2_1pt.eps");
  gPad->SaveAs("picPerformance/SigmaP2_1pt.root");

  //
  TH2F * hPullP21Pt = new TH2F("hhPullP21Pt","hhPullP21Pt",6,0,1,50,-5.05,5.05);
  TH2F * hncPullP21Pt = new TH2F("hhncPullP21Pt","hhncPullP21Pt",6,0,1,50,-5.05,5.05);
  chainCosmic->Draw("(dP2-corrP2)/sP2:sqrt(1/meanPt)>>hhPullP21Pt","side>0"+cutAll+cutRun+cutS+cutN120,"");
  chainCosmic->Draw("(dP2-0)/sP2:sqrt(1/meanPt)>>hhncPullP21Pt","side>0"+cutAll+cutRun+cutS+cutN120,"");
  hPullP21Pt->FitSlicesY();
  hncPullP21Pt->FitSlicesY();
  hhPullP21Pt_2->SetMarkerStyle(kmimarkers[0]);
  hhncPullP21Pt_2->SetMarkerStyle(kmimarkers[1]);
  hhPullP21Pt_2->SetMinimum(0);
  hhPullP21Pt_2->SetXTitle("#sqrt{1/p_{t} (1/GeV)}");
  hhPullP21Pt_2->SetYTitle("#sigma_{#Phi}(Unit) ");
  hhPullP21Pt_2NoCor->Draw();
  //hhncPullP21Pt_2->Draw("same");
  gPad->SaveAs("picPerformance/PullP2_1pt.gif");
  gPad->SaveAs("picPerformance/PullP2_1pt.eps");
  gPad->SaveAs("picPerformance/PullP2_1pt.root");
}

void P3resol1Pt(){
  //
  // P3 - Z -DCA resolution as function of the 1/pt
  //
  TH2F * hdP31Pt = new TH2F("hdP31Pt","hdP31Pt",6,0,1.0,100,-5.05,5.05);
  TH2F * hdP31PtNoCor = new TH2F("hdP31PtNoCor","hdP31PtNoCor",6,0,1.0,100,-5.05,5.05);
  chainCosmic->Draw("1000*(dP3-corrP3)/sqrt(2.):sqrt(1/meanPt)>>hdP31Pt","side>0"+cutAll+cutRun+cutS+cutN120,"");
  chainCosmic->Draw("1000*(dP3-0)/sqrt(2.):sqrt(1/meanPt)>>hdP31PtNoCor","side>0"+cutAll+cutRun+cutS+cutN120,"");
  hdP31Pt->FitSlicesY();
  hdP31PtNoCor->FitSlicesY();
  hdP31Pt_2->SetMinimum(0);
  hdP31Pt_2->SetXTitle("#sqrt{1/p_{t} (1/GeV)}");
  hdP31Pt_2->SetYTitle("#sigma_{#Theta} (mrad)");
  hdP31Pt_2->SetMarkerStyle(kmimarkers[0]);  
  hdP31PtNoCor_2->SetMarkerStyle(kmimarkers[1]);
  hdP31Pt_2NoCor->Draw();  
  //  hdP31PtNoCor_2->Draw("same");
  gPad->SaveAs("picPerformance/SigmaP3_1pt.gif");
  gPad->SaveAs("picPerformance/SigmaP3_1pt.eps");
  gPad->SaveAs("picPerformance/SigmaP3_1pt.root");

  //
  TH2F * hPullP31Pt = new TH2F("hhPullP31Pt","hhPullP31Pt",6,0,1,50,-5.05,5.05);
  TH2F * hncPullP31Pt = new TH2F("hhncPullP31Pt","hhncPullP31Pt",6,0,1,50,-5.05,5.05);
  chainCosmic->Draw("(dP3-corrP3)/sP3:sqrt(1/meanPt)>>hhPullP31Pt","side>0"+cutAll+cutRun+cutS+cutN120,"");
  chainCosmic->Draw("(dP3-0)/sP3:sqrt(1/meanPt)>>hhncPullP31Pt","side>0"+cutAll+cutRun+cutS+cutN120,"");
  hPullP31Pt->FitSlicesY();
  hncPullP31Pt->FitSlicesY();
  hhPullP31Pt_2->SetMarkerStyle(kmimarkers[0]);
  hhncPullP31Pt_2->SetMarkerStyle(kmimarkers[1]);
  hhPullP31Pt_2->SetMinimum(0);
  hhPullP31Pt_2->SetXTitle("#sqrt{1/p_{t} (1/GeV)}");
  hhPullP31Pt_2->SetYTitle("#sigma_{#Theta}(Unit) ");
  hhPullP31Pt_2NoCor->Draw();
  //hhncPullP31Pt_2->Draw("same");
  gPad->SaveAs("picPerformance/PullP3_1pt.gif");
  gPad->SaveAs("picPerformance/PullP3_1pt.eps");
  gPad->SaveAs("picPerformance/PullP3_1pt.root");
}






void P4resol1Pt(){
  //
  // P4 - 1/Pt resolution as function of the 1/pt
  //
  TH2F * hdP41Pt = new TH2F("hdP41Pt","hdP41Pt",6,0,1.0,100,-0.05,0.05);
  TH2F * hdP41PtNoCor = new TH2F("hdP41PtNoCor","hdP41PtNoCor",6,0,1.0,100,-0.05,0.05);
  chainCosmic->Draw("(dP4-corrP4)/sqrt(2.):sqrt(1/meanPt)>>hdP41Pt","side>0"+cutAll+cutRun+cutS+cutN120,"");
  chainCosmic->Draw("(dP4-0)/sqrt(2.):sqrt(1/meanPt)>>hdP41PtNoCor","side>0"+cutAll+cutRun+cutS+cutN120,"");
  hdP41Pt->FitSlicesY();
  hdP41PtNoCor->FitSlicesY();
  hdP41Pt_2->SetMinimum(0);
  hdP41Pt_2->SetXTitle("#sqrt{1/p_{t} (1/GeV)}");
  hdP41Pt_2->SetYTitle("#sigma_{1/pt} (1/GeV)");
  hdP41Pt_2->SetMarkerStyle(kmimarkers[0]);  
  hdP41PtNoCor_2->SetMarkerStyle(kmimarkers[1]);
  hdP41Pt_2NoCor->Draw();  
  //hdP41PtNoCor_2->Draw("same");
  gPad->SaveAs("picPerformance/SigmaP4_1pt.gif");
  gPad->SaveAs("picPerformance/SigmaP4_1pt.eps");
  gPad->SaveAs("picPerformance/SigmaP4_1pt.root");

  //
  TH2F * hPullP41Pt = new TH2F("hhPullP41Pt","hhPullP41Pt",6,0,1,50,-5.05,5.05);
  TH2F * hncPullP41Pt = new TH2F("hhncPullP41Pt","hhncPullP41Pt",6,0,1,50,-5.05,5.05);
  chainCosmic->Draw("(dP4-corrP4)/sP4:sqrt(1/meanPt)>>hhPullP41Pt","side>0"+cutAll+cutRun+cutS+cutN120,"");
  chainCosmic->Draw("(dP4-0)/sP4:sqrt(1/meanPt)>>hhncPullP41Pt","side>0"+cutAll+cutRun+cutS+cutN120,"");
  hPullP41Pt->FitSlicesY();
  hncPullP41Pt->FitSlicesY();
  hhPullP41Pt_2->SetMarkerStyle(kmimarkers[0]);
  hhncPullP41Pt_2->SetMarkerStyle(kmimarkers[1]);
  hhPullP41Pt_2->SetMinimum(0);
  hhPullP41Pt_2->SetXTitle("#sqrt{1/p_{t}} (1/GeV)}");
  hhPullP41Pt_2->SetYTitle("#sigma_{1/pt} (Unit) ");
  hhPullP41Pt_2NoCor->Draw();
  //hhncPullP41Pt_2->Draw("same");
  gPad->SaveAs("picPerformance/PullP4_1pt.gif");
  gPad->SaveAs("picPerformance/PullP4_1pt.eps");
  gPad->SaveAs("picPerformance/PullP4_1pt.root");
}

//////////////////////////////////////////////////////
//
//
// RESOLUTION as function of Z
//
//
//////////////////////////////////////////////////////

void P0resolZ(){
  //
  //
  //
  TH2F * hdP0Z = new TH2F("hdP0Z","hdP0Z",10,-250,250,100,-3.05,3.05);
  TH2F * hdP0ZNoCor = new TH2F("hdP0ZNoCor","hdP0ZNoCor",10,-250,250,100,-3.05,3.05);
  chainCosmic->Draw("(dP0-corrP0)/sqrt(2.):x02>>hdP0Z",cutAll+cutRun+cutS+cutN120,"");
  chainCosmic->Draw("(dP0-0)/sqrt(2.):x02>>hdP0ZNoCor",cutAll+cutRun+cutS+cutN120,"");
  hdP0Z->FitSlicesY();
  hdP0ZNoCor->FitSlicesY();
  //
  TH2F * hsP0Z = new TH2F("hsP0Z","hsP0Z",10,-250,250,200,0.0,2);
  chainCosmic->Draw("(sP0)/sqrt(2.):x02>>hsP0Z",cutAll+cutRun+cutS+cutN120,"");  
  TH1 * hsP0Z_2 = hsP0Z->ProfileX();

  //
  hdP0Z_2->SetMinimum(0);
  hdP0Z_2->SetXTitle("Z position (cm)");
  hdP0Z_2->SetYTitle("#sigma_{y} (cm)");
  hdP0Z_2->SetMarkerStyle(kmimarkers[0]);  
  hdP0ZNoCor_2->SetMarkerStyle(kmimarkers[1]);
  hdP0Z_2->Draw();  
  hsP0Z_2->Draw("same");  
  hdP0ZNoCor_2->Draw("same");

  gPad->SaveAs("picPerformance/SigmaP0_z.gif");
  gPad->SaveAs("picPerformance/SigmaP0_z.eps");
  gPad->SaveAs("picPerformance/SigmaP0_z.root");
  //
  TH2F * hdPP0Z = new TH2F("hdPP0Z","hdPP0Z",8,-200,200,50,-5.05,5.05);
  TH2F * hncdPP0Z = new TH2F("hncdPP0Z","hncdPP0Z",8,-200,200,50,-5.05,5.05);
  chainCosmic->Draw("(dP0-corrP0)/sP0:x02>>hdPP0Z",cutAll+cutRun+cutS+cutN120,"");
  chainCosmic->Draw("(dP0-0)/sP0:x02>>hncdPP0Z",cutAll+cutRun+cutS+cutN120,"");
  hdPP0Z->FitSlicesY();
  hncdPP0Z->FitSlicesY();
  hdPP0Z_2->SetMarkerStyle(kmimarkers[0]);
  hncdPP0Z_2->SetMarkerStyle(kmimarkers[1]);
  hdPP0Z_2->SetMinimum(0);
  hdPP0Z_2->SetXTitle("Z position (cm)");
  hdPP0Z_2->SetYTitle("#sigma_{y}(Unit) ");
  hdPP0Z_2->Draw();
  hncdPP0Z_2->Draw("same");
  gPad->SaveAs("picPerformance/PullP0_z.gif");
  gPad->SaveAs("picPerformance/PullP0_z.eps");
  gPad->SaveAs("picPerformance/PullP0_z.root");
}

void P1resolZ(){
  //
  //
  //
  TH2F *hdP1Z = new TH2F("hdP1Z","hdP1Z",10,-250,250,100,-1.05,1.05);
  TH2F *hdP1ZNoCor=new TH2F("hdP1ZNoCor","hdP1ZNoCor",10,-250,250,100,-1.05,1.05);
  chainCosmic->Draw("(dP1-corrP1)/sqrt(2.):x02>>hdP1Z",cutAll+cutRun+cutS+cutN120,"");
  chainCosmic->Draw("(dP1-0)/sqrt(2.):x02>>hdP1ZNoCor",cutAll+cutRun+cutS+cutN120,"");
  TH2F * hsP1Z = new TH2F("hsP1Z","hsP1Z",10,-250,250,200,0.0,1);
  chainCosmic->Draw("(sP1)/sqrt(2.):x02>>hsP1Z",cutAll+cutRun+cutS+cutN120,"");  
  TH1 * hsP1Z_2 = hsP1Z->ProfileX();

  hdP1Z->FitSlicesY();
  hdP1ZNoCor->FitSlicesY();
  hdP1Z_2->SetMinimum(0);
  hdP1Z_2->SetXTitle("Z position (cm)");
  hdP1Z_2->SetYTitle("#sigma_{z} (cm)");
  hdP1Z_2->SetMarkerStyle(kmimarkers[0]);  
  hsP1Z_2->SetMarkerStyle(kmimarkers[2]);  
  hdP1ZNoCor_2->SetMarkerStyle(kmimarkers[1]);

  hdP1ZNoCor_2->SetMinimum(0);
  hdP1Z_2->Draw("");
  hdP1ZNoCor_2->Draw("same");
  hsP1Z_2->Draw("same");

  gPad->SaveAs("picPerformance/SigmaP1_z.gif");
  gPad->SaveAs("picPerformance/SigmaP1_z.eps");
  gPad->SaveAs("picPerformance/SigmaP1_z.root");
  //
  TH2F * hdPP1Z = new TH2F("hdPP1Z","hdPP1Z",8,-200,200,50,-5.05,5.05);
  TH2F * hncdPP1Z = new TH2F("hncdPP1Z","hncdPP1Z",8,-200,200,50,-5.05,5.05);
  chainCosmic->Draw("(dP1-corrP1)/sP1:x02>>hdPP1Z",cutAll+cutRun+cutS+cutN120,"");
  chainCosmic->Draw("(dP1-0)/sP1:x02>>hncdPP1Z",cutAll+cutRun+cutS+cutN120,"");
  hdPP1Z->FitSlicesY();
  hncdPP1Z->FitSlicesY();
  hdPP1Z_2->SetMarkerStyle(kmimarkers[0]);
  hncdPP1Z_2->SetMarkerStyle(kmimarkers[1]);
  hdPP1Z_2->SetMinimum(0);
  hncdPP1Z_2->SetMinimum(0);
  hncdPP1Z_2->SetXTitle("Z position (cm)");
  hncdPP1Z_2->SetYTitle("#sigma_{z} (Unit) ");
  hncdPP1Z_2->Draw("");
  hdPP1Z_2->Draw();
  gPad->SaveAs("picPerformance/PullP1_z.gif");
  gPad->SaveAs("picPerformance/PullP1_z.eps");
  gPad->SaveAs("picPerformance/PullP1_z.root");
}


void P2resolZ(){
  //
  //
  //
  TH2F * hdP2Z = new TH2F("hdP2Z","hdP2Z",10,-250,250,50,-20.0,20.0);
  chainCosmic->Draw("1000*(dP2-corrP2)/sqrt(2.):x02>>hdP2Z",cutAll+cutRun+cutS+cutN120,"");
  hdP2Z->FitSlicesY();
  hdP2Z_2->SetXTitle("Z position (cm)");
  hdP2Z_2->SetYTitle("#sigma_{#phi} (mrad)");
  hdP2Z_2->Draw();
  gPad->SaveAs("picPerformance/SigmaP2_z.gif");
  gPad->SaveAs("picPerformance/SigmaP2_z.eps");
  gPad->SaveAs("picPerformance/SigmaP2_z.root");

  //
  TH2F * hdPP2Z = new TH2F("hdPP2Z","hdPP2Z",8,-200,200,50,-5.05,5.05);
  chainCosmic->Draw("(dP2-corrP2)/sP2:x02>>hdPP2Z",cutAll+cutRun+cutS+cutN120,"");
  hdPP2Z->FitSlicesY();
  hdPP2Z_2->SetXTitle("Z position (cm)");
  hdPP2Z_2->SetYTitle("#sigma_{#Phi} (Unit) ");
  hdPP2Z_2->Draw();
  gPad->SaveAs("picPerformance/PullP2_z.gif");
  gPad->SaveAs("picPerformance/PullP2_z.eps");
  gPad->SaveAs("picPerformance/PullP2_z.root");

}

void P3resolZ(){
  //
  //
  //
  TH2F * hdP3Z= new TH2F("hdP3Z","hdP3Z",10,-250,250,50,-5,5);
  chainCosmic->Draw("1000*(dP3-corrP3)/sqrt(2.):x02>>hdP3Z",cutAll+cutRun+cutS+cutN120,"");
  hdP3Z->FitSlicesY();
  hdP3Z_2->SetMinimum(0);
  hdP3Z_2->SetXTitle("Z position (cm)");
  hdP3Z_2->SetYTitle("#sigma_{#Theta}  (mrad)");
  hdP3Z_2->Draw();
  gPad->SaveAs("picPerformance/SigmaP3_z.gif");
  gPad->SaveAs("picPerformance/SigmaP3_z.eps");
  gPad->SaveAs("picPerformance/SigmaP3_z.root");
  //
  TH2F * hdPP3Z=  new TH2F("hdPP3Z","hdPP3Z",8,-200,200,50,-5.05,5.05);
  chainCosmic->Draw("(dP3-corrP3)/sP3:x02>>hdPP3Z",cutAll+cutRun+cutS+cutN120,"");
  hdPP3Z->FitSlicesY();
  hdPP3Z->SetMinimum(0);
  hdPP3Z_2->SetXTitle("Z position (cm)");
  hdPP3Z_2->SetYTitle("#sigma_{#Theta} (Unit) ");
  hdPP3Z_2->Draw();
  //
  gPad->SaveAs("picPerformance/PullP3_z.gif");
  gPad->SaveAs("picPerformance/PullP3_z.eps");
  gPad->SaveAs("picPerformance/PullP3_z.root");
}



void P4resolZ(){
  //
  //
  //
  TH2F *hdP4Z = new TH2F("hdP4Z","hdP4Z",10,-250,250,100,-0.05,0.05);
  TH2F *hdP4ZNoCor=new TH2F("hdP4ZNoCor","hdP4ZNoCor",10,-250,250,100,-0.05,0.05);
  chainCosmic->Draw("(dP4-corrP4)/sqrt(2.):x02>>hdP4Z",cutAll+cutRun+cutS+cutN120,"");
  chainCosmic->Draw("(dP4-0)/sqrt(2.):x02>>hdP4ZNoCor",cutAll+cutRun+cutS+cutN120,"");
  hdP4Z->FitSlicesY();
  hdP4ZNoCor->FitSlicesY();
  hdP4Z_2->SetMinimum(0);
  hdP4Z_2->SetXTitle("Z position (cm)");
  hdP4Z_2->SetYTitle("#sigma_{1/pt} (1/GeV)");
  hdP4Z_2->SetMarkerStyle(kmimarkers[0]);  
  hdP4ZNoCor_2->SetMarkerStyle(kmimarkers[1]);
  hdP4Z_2->Draw();  
  //hdP4ZNoCor_2->Draw("same");
  gPad->SaveAs("picPerformance/SigmaP4_z.gif");
  gPad->SaveAs("picPerformance/SigmaP4_z.eps");
  gPad->SaveAs("picPerformance/SigmaP4_z.root");
  //
  TH2F * hdPP4Z = new TH2F("hdPP4Z","hdPP4Z",8,-200,200,50,-5.05,5.05);
  TH2F * hncdPP4Z = new TH2F("hncdPP4Z","hncdPP4Z",8,-200,200,50,-5.05,5.05);
  chainCosmic->Draw("(dP4-corrP4)/sP4:x02>>hdPP4Z",cutAll+cutRun+cutS+cutN120,"");
  chainCosmic->Draw("(dP4-0)/sP4:x02>>hncdPP4Z",cutAll+cutRun+cutS+cutN120,"");
  hdPP4Z->FitSlicesY();
  hncdPP4Z->FitSlicesY();
  hdPP4Z_2->SetMarkerStyle(kmimarkers[0]);
  hncdPP4Z_2->SetMarkerStyle(kmimarkers[1]);
  hdPP4Z_2->SetMinimum(0);
  hdPP4Z_2->SetXTitle("Z position (cm)");
  hdPP4Z_2->SetYTitle("#sigma_{1/pt} (Unit) ");
  hdPP4Z_2->Draw();
  //hncdPP4Z_2->Draw("same");
  gPad->SaveAs("picPerformance/PullP4_z.gif");
  gPad->SaveAs("picPerformance/PullP4_z.eps");
  gPad->SaveAs("picPerformance/PullP4_z.root");

}

*/
