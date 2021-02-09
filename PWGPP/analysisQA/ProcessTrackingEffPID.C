#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TH1.h>
#include <TH2.h>
#include <TAxis.h>
#include <TFile.h>
#include <TStyle.h>
#include <TMath.h>
#include <THnSparse.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include "AliCFContainer.h"
#include "AliPID.h"
#endif

void SetHistoStyle(TH1* h,  Int_t col, Int_t mar, Int_t opt);
int GetEmptyMarker(int mar);
int GetLighterCol(int col);
int GetDarkerCol(int col);
TH1D* Project(THnSparseF* hSparse,Int_t var, Double_t ptmin, Double_t ptmax);

Bool_t firstCall=kTRUE;
TH1D* hToSave[100];
Int_t nToSave=0;


void ProcessTrackingEffPID(TString filname="AnalysisResults.root",TString suffix=""){


  // Graphical setting for particle species
  TString partname[AliPID::kSPECIESC]={"e","#mu","#pi","K","p","d","t","3He","#alpha"};
  TString partnamepos[AliPID::kSPECIESC]={"e^{+}","#mu^{+}","#pi^{+}","K^{+}","p","d","t","3He","#alpha"};
  TString partnameneg[AliPID::kSPECIESC]={"e^{-}","#mu^{-}","#pi^{-}","K^{-}","#bar{p}","#bar{d}","#bar{t}","#bar{3He}","#bar{#alpha}"};
  Bool_t show[AliPID::kSPECIESC]={1,0,1,1,1,1,0,0,0};
  Int_t color[AliPID::kSPECIESC]={kGreen+2,kYellow+2,1,kRed+1,kBlue+1,kMagenta+1,kOrange+2,kCyan+2,kGray+1};
  Int_t marker[AliPID::kSPECIESC]={22,29,20,21,33,23,29,21,23};

  // Steps in efficiency calculations
  const Int_t nSteps=5;
  TString hName[nSteps]={"hGen","hGenEvSel","hReconstructed","hReconstructedPID","hReconstructedTOF"};
  const Int_t nRatios=5;
  Int_t hRatioNum[nRatios]={2,2,3,4,3};
  Int_t hRatioDen[nRatios]={0,1,1,1,4};
  TString ratioTit[nRatios]={"Reconstructed/Generated",
			     "Reconstructed/(Generated in selected events)",
			     "Reconstructed+TOF/(Generated in selected events)",
			     "Reconstructed+PID/(Generated in selected events)",
			     "Reconstructed+PID/(Reconstrcuted+TOF)"};
  Int_t markerStep[nRatios]={28,20,29,27,25};
  TString multVar="N_{tracklets}";
  
  // different projection variables
  TString varname[5]={"Eta","Phi","Pt","Mult","Zvert"};
  const Int_t nProjections=9;
  Double_t theVar[nProjections]={2,1,1,0,0,3,3,4,4};
  const Int_t nPtBins=3; // full pt, low pt, high pt
  Double_t ptLow[nPtBins]={-1.,0.3,3.};
  Double_t ptHigh[nPtBins]={999999.,0.5,10.};
  Int_t thePtBin[nProjections]={0,1,2,1,2,1,2,1,2};
  Double_t maxMult=10000.;

  Int_t nToShow=0;
  for (int iSpecies = 0; iSpecies < AliPID::kSPECIESC; iSpecies++) {
    if(show[iSpecies]) ++nToShow;
  }
  
  TFile* f = new TFile(filname.Data());
  TDirectoryFile* d = (TDirectoryFile*)f->Get("TrackEffPID");
  if(!d){
    printf("TDirectoryFile TrackEffPID does not exist\n");
    f->ls();
    return;
  }
  TList* l = (TList*)d->Get(Form("listTrackEffPID%s",suffix.Data()));
  if(!l){
    printf("TList listTrackEffPID%s does not exist\n",suffix.Data());
    d->ls();
    return;
  }
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);



      
  TH1D* hEff[AliPID::kSPECIESC][3][nProjections];
  TH1D* hGen[AliPID::kSPECIESC][3][nProjections];
  TH1D* hRec[AliPID::kSPECIESC][3][nProjections];
  TH1D* hEffRatio[AliPID::kSPECIESC][nProjections];

  TH1D* hPtStep[AliPID::kSPECIESC][nSteps];
  TH1D* hPtRatio[AliPID::kSPECIESC][nRatios];
  for (int iSpecies = 0; iSpecies < AliPID::kSPECIESC; iSpecies++) {
    for (int iStep=0; iStep<nSteps; iStep++){
      TString namepos=Form("%s_%s_pos",hName[iStep].Data(),AliPID::ParticleShortName(iSpecies));
      TString nameneg=Form("%s_%s_neg",hName[iStep].Data(),AliPID::ParticleShortName(iSpecies));
      THnSparseF* hSparseP=(THnSparseF*)l->FindObject(namepos.Data());
      THnSparseF* hSparseN=(THnSparseF*)l->FindObject(nameneg.Data());
      TH1D* htmp1=Project(hSparseP,2,-1.,9999.);
      TH1D* htmp2=Project(hSparseN,2,-1.,9999.);
      htmp1->Add(htmp2);
      hPtStep[iSpecies][iStep]=(TH1D*)htmp1->Clone(Form("%s_Pt_%s_posneg",hName[iStep].Data(),AliPID::ParticleShortName(iSpecies)));
      delete htmp1;
      delete htmp2;
    }
  }
  for (int iSpecies = 0; iSpecies < AliPID::kSPECIESC; iSpecies++) {
    for (int iRatio = 0; iRatio < nRatios; iRatio++) {
      Int_t theNum=hRatioNum[iRatio];
      Int_t theDen=hRatioDen[iRatio];
      hPtRatio[iSpecies][iRatio]=(TH1D*)hPtStep[iSpecies][theNum]->Clone(Form("hRatio%d_%d",theNum,theDen));
      hPtRatio[iSpecies][iRatio]->Divide(hPtStep[iSpecies][theNum],hPtStep[iSpecies][theDen],1,1,"B");
      SetHistoStyle(hPtRatio[iSpecies][iRatio],color[iSpecies],markerStep[iRatio],0);
    }
  }


  TString charge[2] = {"pos","neg"};
  for (int iSpecies = 0; iSpecies < AliPID::kSPECIESC; iSpecies++) {
    for (int iCharge = 0; iCharge < 2; ++iCharge) {
      TString nameSpGen=Form("hGenEvSel_%s_%s",AliPID::ParticleShortName(iSpecies),charge[iCharge].Data());
      TString nameSpRec=Form("hReconstructed_%s_%s",AliPID::ParticleShortName(iSpecies),charge[iCharge].Data());
      THnSparseF* hSpGen=(THnSparseF*)l->FindObject(nameSpGen.Data());
      THnSparseF* hSpRec=(THnSparseF*)l->FindObject(nameSpRec.Data());
      printf("%s %s \n",nameSpRec.Data(),nameSpGen.Data());
      maxMult=hSpRec->GetAxis(3)->GetXmax();
      TString tit3=hSpRec->GetAxis(3)->GetTitle();
      if(tit3.Contains("b (fm")) multVar="b (fm)";
      else if(tit3.Contains("in cone")) multVar="N tracks in cone";
      else multVar=tit3.Data();
      for(Int_t iP=0; iP<nProjections; iP++){
	Int_t iVar=theVar[iP];
	Int_t iPt=thePtBin[iP];
	Double_t ptmin=ptLow[iPt];
	Double_t ptmax=ptHigh[iPt];
	TString ptrange="";
	if(iPt==1) ptrange="_LowPt";
	if(iPt==2) ptrange="_HighPt";
	hGen[iSpecies][iCharge][iP]=Project(hSpGen,iVar,ptmin,ptmax);
	hGen[iSpecies][iCharge][iP]->SetName(Form("hGen%s_%s_%s%s",varname[iVar].Data(),AliPID::ParticleShortName(iSpecies),charge[iCharge].Data(),ptrange.Data()));
	hRec[iSpecies][iCharge][iP]=Project(hSpRec,iVar,ptmin,ptmax);
 	hRec[iSpecies][iCharge][iP]->SetName(Form("hRec%s_%s_%s%s",varname[iVar].Data(),AliPID::ParticleShortName(iSpecies),charge[iCharge].Data(),ptrange.Data()));
     }
    }
    for(Int_t iP=0; iP<nProjections; iP++){
      TString nameg=hGen[iSpecies][0][iP]->GetName();
      nameg.ReplaceAll("_pos","_posneg");
      hGen[iSpecies][2][iP]=(TH1D*)hGen[iSpecies][0][iP]->Clone(nameg.Data());
      hGen[iSpecies][2][iP]->Add(hGen[iSpecies][1][iP]);
      TString namer=hRec[iSpecies][0][iP]->GetName();
      namer.ReplaceAll("_pos","_posneg");
      hRec[iSpecies][2][iP]=(TH1D*)hRec[iSpecies][0][iP]->Clone(namer.Data());
      hRec[iSpecies][2][iP]->Add(hRec[iSpecies][1][iP]);
    }
  }
  for (int iSpecies = 0; iSpecies < AliPID::kSPECIESC; iSpecies++) {
    for (int iCharge = 0; iCharge < 3; ++iCharge) {
      for(Int_t iP=0; iP<nProjections; iP++){
	TString namee=hRec[iSpecies][iCharge][iP]->GetName();
	namee.ReplaceAll("hRec","hEff");
	hEff[iSpecies][iCharge][iP]=(TH1D*)hRec[iSpecies][iCharge][iP]->Clone(namee.Data());
	hEff[iSpecies][iCharge][iP]->Divide(hRec[iSpecies][iCharge][iP],hGen[iSpecies][iCharge][iP],1.,1.,"B");
	SetHistoStyle(hEff[iSpecies][iCharge][iP],color[iSpecies],marker[iSpecies],iCharge);
      }
    }
  }
  for (int iSpecies = 0; iSpecies < AliPID::kSPECIESC; iSpecies++) {
    for(Int_t iP=0; iP<nProjections; iP++){
      TString namer=hEff[iSpecies][0][iP]->GetName();
      namer.ReplaceAll("hEff","hRatioPosNegEff");
      namer.ReplaceAll("_pos","");
      hEffRatio[iSpecies][iP]=(TH1D*)hEff[iSpecies][0][iP]->Clone(namer.Data());
      hEffRatio[iSpecies][iP]->Divide(hEff[iSpecies][1][iP]);
      SetHistoStyle(hEffRatio[iSpecies][iP],kGray+1,marker[iSpecies],0);
    }
  }

  TH2F* hFramePt=new TH2F("hFramePt","",100,0.05,30.,100.,0.,1.2);
  hFramePt->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  hFramePt->GetYaxis()->SetTitle("Efficiency");
  hFramePt->GetYaxis()->SetTitleOffset(1.2);

  TH2F* hFrameRatio=new TH2F("hFrameRatio","",100,0.05,30.,100.,0.8,1.2);
  hFrameRatio->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  hFrameRatio->GetYaxis()->SetTitle("Ratio +/- efficiency");
  hFrameRatio->GetYaxis()->SetTitleOffset(1.2);

  TH2F* hFramePhi=new TH2F("hFramePhi","",100,0.,2.*TMath::Pi(),100.,0.,1.2);
  hFramePhi->GetXaxis()->SetTitle("#varphi (rad)");
  hFramePhi->GetYaxis()->SetTitle("Efficiency");
  hFramePhi->GetYaxis()->SetTitleOffset(1.2);

  TH2F* hFrameEta=new TH2F("hFrameEta","",100,-1.,1.,100.,0.,1.2);
  hFrameEta->GetXaxis()->SetTitle("#eta");
  hFrameEta->GetYaxis()->SetTitle("Efficiency");
  hFrameEta->GetYaxis()->SetTitleOffset(1.2);

  TH2F* hFrameZvert=new TH2F("hFrameZvert","",100,-12.,12.,100.,0.,1.2);
  hFrameZvert->GetXaxis()->SetTitle("z_{vertex} (cm)");
  hFrameZvert->GetYaxis()->SetTitle("Efficiency");
  hFrameZvert->GetYaxis()->SetTitleOffset(1.2);

  TH2F* hFrameMult=new TH2F("hFrameMult","",10000,0.,maxMult,100.,0.,1.2);
  hFrameMult->GetXaxis()->SetTitle(multVar.Data());
  hFrameMult->GetYaxis()->SetTitle("Efficiency");
  hFrameMult->GetYaxis()->SetTitleOffset(1.2);

  TLatex* tcuts=new TLatex(0.7,0.2,"#splitline{|#eta|<0.8}{|z_{vert}|<10 cm}");
  tcuts->SetNDC();
  tcuts->SetTextFont(43);
  tcuts->SetTextSize(26);

  TLatex* tsp[AliPID::kSPECIESC];
  for (int iSpecies = 0; iSpecies < AliPID::kSPECIESC; iSpecies++) {
    tsp[iSpecies]=new TLatex(0.18,0.82,AliPID::ParticleName(iSpecies));
    tsp[iSpecies]->SetNDC();
    tsp[iSpecies]->SetTextFont(43);
    tsp[iSpecies]->SetTextSize(24);
  }

  Int_t ny=1;
  if(nToShow>=3) ny=2;
  if(nToShow>=7) ny=3;
  Int_t nx=TMath::Ceil((Float_t)(nToShow+1)/(Float_t)ny);

  TCanvas* cest=new TCanvas("cest","EffStepByStep",1600,800);
  cest->Divide(nx,ny);
  Int_t thePad=1;
  for(Int_t iSp=0; iSp<AliPID::kSPECIESC; iSp++){
    if(!show[iSp]) continue;
    cest->cd(thePad);
    gPad->SetLeftMargin(0.12);
    gPad->SetRightMargin(0.08);
    gPad->SetTickx();
    gPad->SetTicky();
    gPad->SetLogx();
    hFramePt->Draw();
    for(Int_t iRatio=0; iRatio<nRatios; iRatio++){
      hPtRatio[iSp][iRatio]->Draw("same");
    }
    tsp[iSp]->Draw();
    ++thePad;
  }
  cest->cd(thePad);  
  TLegend* legra=new TLegend(0.1,0.27,0.9,0.73);
  legra->SetBorderSize(1);
  legra->SetMargin(0.15);
  for(Int_t iRatio=0; iRatio<nRatios; iRatio++){
    legra->AddEntry(hPtRatio[2][iRatio],ratioTit[iRatio].Data(),"P");
  }
  legra->Draw();

  ny=1;
  nx=nToShow;
  if(nToShow>5){
    ny=2;
    nx=TMath::Ceil((Float_t)nToShow/2.);
  }
  TCanvas* ceptch=new TCanvas("ceptch","EffVsPtByCharge",1600,800);
  ceptch->Divide(nx,ny,0,0);
  thePad=1;
  for(Int_t iSp=0; iSp<AliPID::kSPECIESC; iSp++){
    if(!show[iSp]) continue;
    ceptch->cd(thePad);
    TPad* thisPad=(TPad*)ceptch->GetPad(thePad);
    thisPad->Divide(1,2,0);
    thisPad->cd(1);
    gPad->SetLeftMargin(0.12);
    gPad->SetRightMargin(0.08);
    gPad->SetTickx();
    gPad->SetTicky();
    gPad->SetLogx();
    hFramePt->Draw();
    TLegend* legch=new TLegend(0.6,0.15,0.89,0.4);
    legch->SetTextFont(43);
    legch->SetTextSize(22);
    for(Int_t iP=0; iP<nProjections; iP++){
      if(theVar[iP]==2){
	hEff[iSp][0][iP]->Draw("same");
	hEff[iSp][1][iP]->Draw("same");
	legch->AddEntry(hEff[iSp][0][iP],partnamepos[iSp].Data(),"P");
	legch->AddEntry(hEff[iSp][1][iP],partnameneg[iSp].Data(),"P");
      }
    }
    legch->Draw();
    thisPad->cd(2);
    gPad->SetLeftMargin(0.12);
    gPad->SetRightMargin(0.08);
    gPad->SetTickx();
    gPad->SetTicky();
    gPad->SetLogx();
    hFrameRatio->Draw();
    for(Int_t iP=0; iP<nProjections; iP++){
      if(theVar[iP]==2){
	hEffRatio[iSp][iP]->Draw("same");
      }
    }
    tsp[iSp]->Draw();
    ++thePad;
  }

  TCanvas* cept=new TCanvas("cept","EffVsPt",1600,800);
  cept->Divide(2,1);
  cept->cd(1);
  gPad->SetLeftMargin(0.12);
  gPad->SetRightMargin(0.08);
  gPad->SetTickx();
  gPad->SetTicky();
  hFramePt->Draw();
  TLegend* legpt=new TLegend(0.7,0.3,0.89,0.55);
  cept->cd(2);
  gPad->SetLeftMargin(0.12);
  gPad->SetRightMargin(0.08);
  gPad->SetTickx();
  gPad->SetLogx();
  hFramePt->Draw();
  for(Int_t iP=0; iP<nProjections; iP++){
    if(theVar[iP]==2){
      for(Int_t iSp=0; iSp<AliPID::kSPECIESC; iSp++){
	if(!show[iSp]) continue;
	cept->cd(1);
	hEff[iSp][2][iP]->Draw("same");
	cept->cd(2);
	hEff[iSp][2][iP]->Draw("same");
	legpt->AddEntry(hEff[iSp][2][iP],partname[iSp].Data(),"P");
      }
    }
  }
  cept->cd(1);
  legpt->Draw();
  tcuts->Draw();
  cept->cd(2);
  legpt->Draw();
  tcuts->Draw();

  cept->SaveAs("EfficVsPt.png");

  TCanvas** cevar=new TCanvas*[4];
  cevar[0]=new TCanvas("cephi","EffVsPhi",1600,800);
  cevar[1]=new TCanvas("ceeta","EffVsEta",1600,800);
  cevar[2]=new TCanvas("cezv","EffVsZvert",1600,800);
  cevar[3]=new TCanvas("cemult","EffVsMult",1600,800);

  TLegend** legvar=new TLegend*[8];
  for(Int_t iv=0; iv<4; iv++){
    cevar[iv]->Divide(2,1);
    for(Int_t ip=1; ip<=2; ip++){
      cevar[iv]->cd(ip);
      gPad->SetLeftMargin(0.12);
      gPad->SetRightMargin(0.08);
      gPad->SetTickx();
      gPad->SetTicky();
      if(iv==0) hFramePhi->Draw();
      else if(iv==1) hFrameEta->Draw();
      else if(iv==2) hFrameZvert->Draw();
      else if(iv==3) hFrameMult->Draw();
      legvar[iv*2+(ip-1)]=new TLegend(0.16,0.73,0.87,0.89);
      legvar[iv*2+(ip-1)]->SetNColumns(2);
      legvar[iv*2+(ip-1)]->SetColumnSeparation(0.15);
      legvar[iv*2+(ip-1)]->SetMargin(0.15);
    }
  }

  for(Int_t iP=0; iP<nProjections; iP++){
    Int_t iPt=thePtBin[iP];
    Double_t ptmin=ptLow[iPt];
    Double_t ptmax=ptHigh[iPt];
    Int_t iCanv=-1;
    Int_t iPad=-1;
    if(theVar[iP]==1) iCanv=0;
    else if(theVar[iP]==0) iCanv=1;
    else if(theVar[iP]==4) iCanv=2;
    else if(theVar[iP]==3) iCanv=3;
    if(thePtBin[iP]==1) iPad=1;
    else if(thePtBin[iP]==2) iPad=2;
    if(iCanv<0 || iPad<0) continue;
    cevar[iCanv]->cd(iPad);
    for(Int_t iSp=0; iSp<AliPID::kSPECIESC; iSp++){
      if(!show[iSp]) continue;
      hEff[iSp][2][iP]->Draw("same");
      legvar[2*iCanv+iPad-1]->AddEntry(hEff[iSp][2][iP],Form("%s, %.1f<p_{T}<%.1f GeV/c",partname[iSp].Data(),ptmin,ptmax),"P");	
    }
  }
  for(Int_t iv=0; iv<4; iv++){
    cevar[iv]->cd(1);
    legvar[iv*2]->Draw();
    tcuts->Draw();
    cevar[iv]->cd(2);
    legvar[iv*2+1]->Draw();
    tcuts->Draw();
  }
  cevar[0]->SaveAs("EfficVsPhi.png");
  cevar[1]->SaveAs("EfficVsEta.png");
  cevar[2]->SaveAs("EfficVsZvert.png");
  cevar[3]->SaveAs("EfficVsMult.png");

  TFile* filout=new TFile("TrackingEffPID.root","recreate");
  for (int iSpecies = 0; iSpecies < AliPID::kSPECIESC; iSpecies++) {
    for (int iCharge = 0; iCharge < 3; ++iCharge) {
      for(Int_t iP=0; iP<nProjections; iP++){
	hEff[iSpecies][iCharge][iP]->Write();
      }
    }
  }
  filout->Close();
  
}

TH1D* Project(THnSparseF* hSparse, Int_t var, Double_t ptmin, Double_t ptmax){

  Int_t ndim=hSparse->GetNdimensions();

  for(Int_t j=0; j<ndim; j++){
    TAxis* ax=hSparse->GetAxis(j);
    Int_t nbins=ax->GetNbins();
    if(j!=var) ax->SetRange(0,nbins+1);
    else ax->SetRange(1,nbins);
  }

  if(var!=2 && ptmin>0 && ptmax<100){
    TAxis* ax2=hSparse->GetAxis(2);
    Int_t bl=ax2->FindBin(ptmin*1.0001);
    Int_t bu=ax2->FindBin(ptmax*0.9999);
    //    printf("Pt bins considered: %d-%d\n",bl,bu);
    ax2->SetRange(bl,bu);
  }

  //cut at |eta|<0.8 if projecting vs. other variables
  if(var!=0){
    TAxis* ax0=hSparse->GetAxis(0);
    Int_t bl=ax0->FindBin(-0.7999);
    Int_t bu=ax0->FindBin(0.7999);
    //    printf("Eta bins considered: %d-%d\n",bl,bu);
    ax0->SetRange(bl,bu);
  }

  // cut zvertex at 10 cm
  TAxis* ax4=hSparse->GetAxis(4);
  Int_t bl=ax4->FindBin(-9.9999);
  Int_t bu=ax4->FindBin(9.9999);
  ax4->SetRange(bl,bu);
  //  printf("zvertex bins considered: %d-%d\n",bl,bu);

  TH1D* hPr=hSparse->Projection(var);
  hPr->Sumw2();
  return hPr;
}


void SetHistoStyle(TH1* h,  Int_t col, Int_t mar, Int_t opt){
  if(opt==1) mar=GetEmptyMarker(mar);
  if(opt==2) col=GetLighterCol(col);
  h->SetMarkerStyle(mar);
  h->SetMarkerColor(col);
  h->SetLineColor(col);
}

int GetDarkerCol(int col){
  if(col==1) return kGray+2;
  if(col==2) return kRed+1;
  if(col==4) return kBlue+1;
  if(col==kGreen+1) return kGreen+2;
  if(col==6) return kMagenta+1;
  return col;
}

int GetLighterCol(int col){
  if(col==1) return kGray+2;
  if(col==kRed+1) return 2;
  if(col==kBlue+1) return 4;
  if(col==kGreen+2) return kGreen+1;
  if(col==kMagenta+1) return 6;
  if(col==kCyan+2) return kCyan;
  if(col==kOrange+2) return kOrange+1;
  if(col==kGray+1) return kGray;
  if(col==kYellow+2) return kYellow+1;
  return col;
}

int GetEmptyMarker(int mar){
  if(mar==20) return 24;
  if(mar==21) return 25;
  if(mar==22) return 26;
  if(mar==23) return 32;
  if(mar==29) return 30;
  if(mar==33) return 27;
  if(mar==34) return 28;
  return 24;
}
