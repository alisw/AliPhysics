#if !defined(__CINT__) || defined(__MAKECINT__)

#include "Riostream.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TF1.h"
#include "TGraph.h"
#include "TColor.h"
#include "TLegend.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TCutG.h"
#include "TString.h"
#include "TPaveStats.h"
#include "TLine.h"

#endif



//______________________________________________________________________________
const char *workdir = "$HOME/QA_hmpid";

// AliHMPIDParam *pParam;
const Double_t mass[3] = {0.13957,0.49370,0.93823};

// Run by Run Residual Histos:
TH1F *fHmpMipTrkDistPosX[7];
TH1F *fHmpMipTrkDistNegX[7];
TH1F *fHmpMipTrkDistPosY[7];
TH1F *fHmpMipTrkDistNegY[7];
//Run by Run MIP Landau Histos:
TH1F *fHmpMipCharge[7];
// Trend Residual Histos
TH1F *hgHmpXresidNegFromHistoMean[7];
TH1F *hgHmpXresidPosFromHistoMean[7];
TH1F *hgHmpYresidNegFromHistoMean[7];
TH1F *hgHmpYresidPosFromHistoMean[7];
// TREND MIP & Sigma plot.
TH1F *hgHmpMipCharge[7];

void CreateContainers(char * dataType, Int_t year, char *period, char *pass, Int_t number);

void Analyze_QA_run(char * dataType, Int_t year, char *period, char *pass, Int_t run);

void DrawLegend(char * dataType, Int_t year, char *period, char *pass, Int_t run);

void Make_Trending(char * dataType, Int_t year, char *period, char *pass, Int_t run, Int_t index);

void Print_Trending(char * dataType, Int_t year, char *period, char *pass, Int_t numberofruns);

void Save_Trending(char * dataType, Int_t year, char *period, char *pass);

//______________________________________________________________________________
//______________________________________________________________________________
//______________________________________________________________________________


Int_t makeHMPQA(char * dataType, Int_t year, char *period, char *pass){

  gStyle->SetOptStat(0);

  // check if data or simulation
  Bool_t isRealData = kTRUE;
  if(strcmp(dataType,"sim")==0) isRealData = kFALSE;

  ifstream pFile;

  pFile.open( Form("%s_%i_%s_%s.data",dataType, year, period, pass) );

  if(!pFile){cout<<"DataSet: "<<Form("%s_%i_%s_%s.data", dataType, year, period, pass)<<" File not present. Check file name!!!"<<endl; return 0;}

  Int_t number = 0;
  pFile>>number;
  cout<<number<<endl;

  // create output directories
  CreateContainers(dataType, year, period, pass, number);

  for(Int_t ir=0;ir<number;ir++){

    Int_t run=0;

    pFile>>run;

    cout<<run<<endl;

    // qa run...
    gStyle->SetOptStat(1);

    Analyze_QA_run(dataType, year, period, pass, run);

    //trending...
    gStyle->SetOptStat(0);
    Make_Trending(dataType, year, period, pass, run, ir);

  }

  Print_Trending(dataType, year, period, pass, number);

  Save_Trending(dataType, year, period, pass);

  return 0;
}

//______________________________________________________________________________
void CreateContainers(char * dataType, Int_t year, char *period, char *pass, Int_t number){

  gSystem->Exec(Form("mkdir -p %s/%s/%i/%s/%s", workdir, dataType, year, period, pass));

  gSystem->Exec(Form("mkdir -p %s/%s/%i/%s/%s/pdf", workdir, dataType, year, period, pass));

  gSystem->Exec(Form("mkdir -p %s/%s/%i/%s/%s/rootfiles", workdir, dataType, year, period, pass) );

  gSystem->Exec(Form("mkdir -p %s/%s/%i/%s/%s/trending", workdir, dataType, year, period, pass) );

  for(Int_t ich = 0 ; ich < 7; ich++){

    //TREND Residual Histos Container Objects:

    hgHmpXresidNegFromHistoMean[ich] = new TH1F(Form("hgHmpXresidNegFromHistoMean%d",ich),Form("Neg: Xresid From histo Mean;;Neg: Xresid From histo Mean (cm)"),number,0,number);
    hgHmpXresidNegFromHistoMean[ich]->GetXaxis()->SetLabelSize(0.02);
    hgHmpXresidNegFromHistoMean[ich]->Sumw2();
    hgHmpXresidNegFromHistoMean[ich]->SetMarkerStyle(24);    hgHmpXresidNegFromHistoMean[ich]->SetMarkerColor(ich+1);  hgHmpXresidNegFromHistoMean[ich]->SetLineColor(ich+1);

    hgHmpXresidPosFromHistoMean[ich] = new TH1F(Form("hgHmpXresidPosFromHistoMean%d",ich),Form("Pos: Xresid From histo Mean;;Pos: Xresid From histo Mean (cm)"),number,0,number);
    hgHmpXresidPosFromHistoMean[ich]->GetXaxis()->SetLabelSize(0.02);
    hgHmpXresidPosFromHistoMean[ich]->Sumw2();
    hgHmpXresidPosFromHistoMean[ich]->SetMarkerStyle(24);    hgHmpXresidPosFromHistoMean[ich]->SetMarkerColor(ich+1);  hgHmpXresidPosFromHistoMean[ich]->SetLineColor(ich+1);

    hgHmpYresidNegFromHistoMean[ich] = new TH1F(Form("hgHmpYresidNegFromHistoMean%d",ich),Form("Neg: Yresid From histo Mean;;Neg: Yresid From histo Mean (cm)"),number,0,number);
    hgHmpYresidNegFromHistoMean[ich]->GetXaxis()->SetLabelSize(0.02);
    hgHmpYresidNegFromHistoMean[ich]->Sumw2();
    hgHmpYresidNegFromHistoMean[ich]->SetMarkerStyle(24);    hgHmpYresidNegFromHistoMean[ich]->SetMarkerColor(ich+1);  hgHmpYresidNegFromHistoMean[ich]->SetLineColor(ich+1);

    hgHmpYresidPosFromHistoMean[ich] = new TH1F(Form("hgHmpYresidPosFromHistoMean%d",ich),Form("Pos: Yresid From histo Mean;;Pos: Yresid From histo Mean (cm)"),number,0,number);
    hgHmpYresidPosFromHistoMean[ich]->GetXaxis()->SetLabelSize(0.02);
    hgHmpYresidPosFromHistoMean[ich]->Sumw2();
    hgHmpYresidPosFromHistoMean[ich]->SetMarkerStyle(24);    hgHmpYresidPosFromHistoMean[ich]->SetMarkerColor(ich+1);  hgHmpYresidPosFromHistoMean[ich]->SetLineColor(ich+1);

    //TREND MIP Histos Container Objects:

    hgHmpMipCharge[ich] = new TH1F(Form("hgHmpMipChargeCh%d",ich),Form("MIP Charge (MPV/Sigma);;MIP charge (MPV/Sigma)"),number,0,number);
    hgHmpMipCharge[ich] -> Sumw2();
    hgHmpMipCharge[ich] -> GetXaxis()->SetLabelSize(0.02);

  }
}

//______________________________________________________________________________
void Analyze_QA_run(char * dataType, Int_t year, char *period, char *pass, Int_t run){

  //open QA file.
  TFile *fin_qa_results = TFile::Open( Form("%s/%s/%i/%s/%s/QAresults_%i.root",workdir, dataType, year, period, pass, run) );

  TList *list = (TList*)fin_qa_results->Get("HmpidQA/HmpidQA");

  if(!list) return;

  TCanvas * cmain = new TCanvas("cmain","cmain",1200,900);
  cmain->SaveAs(Form("%s/%s/%i/%s/%s/pdf/qa_run_%i.pdf[", workdir, dataType, year, period, pass, run) );

  // check if data or simulation
  Bool_t isRealData = kTRUE;

  if(strcmp(dataType,"sim")==0) isRealData = kFALSE;

  TFile *fout = new TFile(Form("%s/%s/%i/%s/%s/rootfiles/hmpid_qa_run_%i.root", workdir, dataType, year, period, pass, run),"RECREATE");

  for(Int_t ich = 0 ; ich < 7; ich++){

    // X/Y Residual Histos
    fHmpMipTrkDistPosX[ich] = (TH1F*) list->FindObject(Form("fHmpMipTrkDistPosX%d",ich));
    fHmpMipTrkDistPosX[ich]->GetXaxis()->SetRangeUser(-10.,10.);
    fHmpMipTrkDistPosX[ich]->SetLineColor(2);
    fHmpMipTrkDistPosX[ich]->SetMarkerColor(2);
    fHmpMipTrkDistPosX[ich]->SetYTitle("Entries");

    fHmpMipTrkDistNegX[ich] = (TH1F*) list->FindObject(Form("fHmpMipTrkDistNegX%d",ich));
    fHmpMipTrkDistNegX[ich]->GetXaxis()->SetRangeUser(-10.,10.);
    fHmpMipTrkDistNegX[ich]->SetLineColor(4);
    fHmpMipTrkDistNegX[ich]->SetMarkerColor(4);
    fHmpMipTrkDistNegX[ich]->SetYTitle("Entries");

    fHmpMipTrkDistPosY[ich] = (TH1F*) list->FindObject(Form("fHmpMipTrkDistPosY%d",ich));
    fHmpMipTrkDistPosY[ich]->GetXaxis()->SetRangeUser(-10.,10.);
    fHmpMipTrkDistPosY[ich]->SetLineColor(2);
    fHmpMipTrkDistPosY[ich]->SetMarkerColor(2);
    fHmpMipTrkDistPosY[ich]->SetYTitle("Entries");

    fHmpMipTrkDistNegY[ich] = (TH1F*) list->FindObject(Form("fHmpMipTrkDistNegY%d",ich));
    fHmpMipTrkDistNegY[ich]->GetXaxis()->SetRangeUser(-10.,10.);
    fHmpMipTrkDistNegY[ich]->SetLineColor(4);
    fHmpMipTrkDistNegY[ich]->SetMarkerColor(4);
    fHmpMipTrkDistNegY[ich]->SetYTitle("Entries");

    //MIP Landau fit Histos:
    fHmpMipCharge[ich] =  (TH1F*) list->FindObject(Form("fHmpMipCharge%d",ich));
    fHmpMipCharge[ich]->SetYTitle("Entries");
  }



  TCanvas *c1 = new TCanvas("c1","c1",1200,900);
  c1->Divide(3,3);
  gStyle->SetOptStat(0);
  c1->cd(1); fHmpMipTrkDistPosX[6]->Draw(); fHmpMipTrkDistNegX[6]->Draw("same");   c1->cd(2); fHmpMipTrkDistPosX[5]->Draw(); fHmpMipTrkDistNegX[5]->Draw("same");   c1->cd(3);  gPad->Clear(); DrawLegend(dataType, year, period, pass, run);
  c1->cd(4); fHmpMipTrkDistPosX[4]->Draw(); fHmpMipTrkDistNegX[4]->Draw("same");   c1->cd(5); fHmpMipTrkDistPosX[3]->Draw(); fHmpMipTrkDistNegX[3]->Draw("same");   c1->cd(6);fHmpMipTrkDistPosX[2]->Draw(); fHmpMipTrkDistNegX[2]->Draw("same");
                                                                                   c1->cd(8); fHmpMipTrkDistPosX[1]->Draw(); fHmpMipTrkDistNegX[1]->Draw("same");   c1->cd(9); fHmpMipTrkDistPosX[0]->Draw(); fHmpMipTrkDistNegX[0]->Draw("same");

  c1->cd(7);
  TLegend *lg = new TLegend(0.15,0.3,0.85,0.98);
  lg->AddEntry(fHmpMipTrkDistPosY[0],"Positive","lp");
  lg->AddEntry(fHmpMipTrkDistNegY[0],"Negative","lp");
  lg->Draw();

  c1->SaveAs(Form("%s/%s/%i/%s/%s/pdf/qa_run_%i.pdf", workdir, dataType, year, period, pass, run) );


  TCanvas *c2 = new TCanvas("c2","c2",1200,900);
  c2->Divide(3,3);
  gStyle->SetOptStat(0);
  c2->cd(1); fHmpMipTrkDistPosY[6]->Draw(); fHmpMipTrkDistNegY[6]->Draw("same");   c2->cd(2); fHmpMipTrkDistPosY[5]->Draw(); fHmpMipTrkDistNegY[5]->Draw("same");   c2->cd(3);  gPad->Clear(); DrawLegend(dataType, year, period, pass, run);
  c2->cd(4); fHmpMipTrkDistPosY[4]->Draw(); fHmpMipTrkDistNegY[4]->Draw("same");   c2->cd(5); fHmpMipTrkDistPosY[3]->Draw(); fHmpMipTrkDistNegY[3]->Draw("same");   c2->cd(6);fHmpMipTrkDistPosY[2]->Draw(); fHmpMipTrkDistNegY[2]->Draw("same");
                                                                                  c2->cd(8); fHmpMipTrkDistPosY[1]->Draw(); fHmpMipTrkDistNegY[1]->Draw("same");   c2->cd(9); fHmpMipTrkDistPosY[0]->Draw(); fHmpMipTrkDistNegY[0]->Draw("same");

  c2->cd(7);
  lg->Draw();

  c2->SaveAs(Form("%s/%s/%i/%s/%s/pdf/qa_run_%i.pdf", workdir, dataType, year, period, pass, run) );

  TCanvas *c3 = new TCanvas("c3","c3",1200,900);
  c3->SetLogz();
  TH2F *fHmpCkovPesd = (TH2F*) list->FindObject(Form("fHmpCkovPesd"));
  fHmpCkovPesd->SetAxisRange(0,0.8,"Y");
  fHmpCkovPesd->SetAxisRange(0,6.0,"X");
  fHmpCkovPesd->SetMinimum(2);
  gStyle->SetOptStat(0);
  fHmpCkovPesd->Draw("colz");

//   TF1 *fPion = new TF1("fPion","TMath::ACos(TMath::Sqrt(0.13957*0.13957+x*x)/(x*1.289))",0.2 , 6.);
//   fPion->SetLineColor(kBlack);
//   fPion->Draw("same");
//
//   TF1 *fKaon = new TF1("fKaon","TMath::ACos(TMath::Sqrt(0.49370*0.49370+x*x)/(x*1.2989))",0.7 , 6.);
//   fKaon->SetLineColor(kBlack);
//   fKaon->Draw("same");
//
//   TF1 *fProton = new TF1("fProton","TMath::ACos(TMath::Sqrt(0.93823*0.93823+x*x)/(x*1.289))",1.21 , 6);
//   fProton->SetLineColor(kBlack);
//   fProton->Draw("same");

  c3->SaveAs(Form("%s/%s/%i/%s/%s/pdf/qa_run_%i.pdf", workdir, dataType, year, period, pass, run) );
  gStyle->SetPalette(1);

  //MIP value distribution Function.

  //MIP CHARGE DISTRIBUTION RANGE:
  //To Change if needed
  Double_t mipFitRange[2];
  /*REAL DATA*/   if(isRealData == kTRUE) { mipFitRange[0]=300; mipFitRange[1]=1200; }
  /*MONTE CARLO*/ if(isRealData == kFALSE) { mipFitRange[0]=500; mipFitRange[1]=2000; }

  for(Int_t ich = 0 ; ich < 7; ich++){
    c1->cd();
    if ( fHmpMipCharge[ich]->GetEntries() < 200 ) {fHmpMipCharge[ich]->Rebin(30);} else { fHmpMipCharge[ich]->Rebin(20); }
    if ( fHmpMipCharge[ich]->GetEntries() > 360 ) {
      fHmpMipCharge[ich]->Fit("landau","Q","",mipFitRange[0],mipFitRange[1]);
      fHmpMipCharge[ich]->Fit("landau","Q","",mipFitRange[0],mipFitRange[1]);
    }
    fHmpMipCharge[ich]->SetAxisRange(200,3000,"X");
    fHmpMipCharge[ich]->Draw();
//    rootfiles->cd();
//    fHmpMipCharge[ich]->Write();
  }// end RICH loop.

//     rootfiles->Close();
  gStyle->SetOptFit(0);
  c1->cd();
  c1->Clear();
  c1->Divide(3,3);
  Int_t padMap[7]={9,8,6,5,4,2,1};   Double_t valneg = -1, valpos = -1;
  for(Int_t ich = 0 ; ich < 7; ich++){
    c1->cd(padMap[ich])->SetGridy();
    c1->cd(padMap[ich])->SetGridx();
    fHmpMipCharge[ich]->Draw(); c1->Update();
  }// end RICH loop.
  c1->SaveAs(Form("%s/%s/%i/%s/%s/pdf/qa_run_%i.pdf", workdir, dataType, year, period, pass, run) );


  cmain->SaveAs(Form("%s/%s/%i/%s/%s/pdf/qa_run_%i.pdf]", workdir, dataType, year, period, pass, run) );

  //Write Run by Run histos in a root file.
  fout->cd();
  fHmpCkovPesd->Write();
  for(Int_t ich = 0 ; ich < 7; ich++){
    // X/Y Residual Histos
    fHmpMipTrkDistPosX[ich]->Write();
    fHmpMipTrkDistNegX[ich]->Write();
    fHmpMipTrkDistPosY[ich]->Write();
    fHmpMipTrkDistNegY[ich]->Write();

    //MIP Landau fit Histos:
    fHmpMipCharge[ich]->Write();
  }

  fout->Close();
}

//______________________________________________________________________________
void DrawLegend(char * dataType, Int_t year, char *period, char *pass, Int_t run){

  //used to draw legend

  TLegend *pLeg=new TLegend(0.15,0.3,0.85,0.98);

  pLeg->AddEntry((TObject*)0, Form("%s",dataType), "");

  pLeg->AddEntry((TObject*)0, Form("%i",year), "");

  pLeg->AddEntry((TObject*)0, Form("%s",period), "");

  pLeg->AddEntry((TObject*)0, Form("%s",pass), "");

  pLeg->AddEntry((TObject*)0, Form("%i",run), "");

  pLeg->Draw();

}//DrawLegend()

//______________________________________________________________________________
void Make_Trending(char * dataType, Int_t year, char *period, char *pass, Int_t run, Int_t index){

  for(Int_t ich = 0 ; ich < 7; ich++){

    // X Residual Negative
    hgHmpXresidNegFromHistoMean[ich]->GetXaxis()->SetBinLabel(index+1,Form("%i",run));
    if(fHmpMipTrkDistNegX[ich]->GetEntries()<100.) { hgHmpXresidNegFromHistoMean[ich] -> SetBinContent(index+1, -999); }
    else { hgHmpXresidNegFromHistoMean[ich] -> SetBinContent(index+1, fHmpMipTrkDistNegX[ich]->GetMean()); }
    hgHmpXresidNegFromHistoMean[ich] -> SetBinError(index+1, 0.000000000000000001);

    // X Residual Positive
    hgHmpXresidPosFromHistoMean[ich]->GetXaxis()->SetBinLabel(index+1,Form("%i",run));
    if(fHmpMipTrkDistPosX[ich]->GetEntries()<100.) { hgHmpXresidPosFromHistoMean[ich] -> SetBinContent(index+1, -999); }
    else { hgHmpXresidPosFromHistoMean[ich] -> SetBinContent(index+1, fHmpMipTrkDistPosX[ich]->GetMean()); }
    hgHmpXresidPosFromHistoMean[ich] -> SetBinError(index+1, 0.000000000000000001);

    // Y Residual Negative
    hgHmpYresidNegFromHistoMean[ich]->GetXaxis()->SetBinLabel(index+1,Form("%i",run));
    if( fHmpMipTrkDistNegY[ich]->GetEntries()<100.) { hgHmpYresidNegFromHistoMean[ich] -> SetBinContent(index+1, -999); }
    else { hgHmpYresidNegFromHistoMean[ich] -> SetBinContent(index+1, fHmpMipTrkDistNegY[ich]->GetMean()); }
    hgHmpYresidNegFromHistoMean[ich] -> SetBinError(index+1, 0.000000000000000001);

    // Y Residual Positive
    hgHmpYresidPosFromHistoMean[ich]->GetXaxis()->SetBinLabel(index+1,Form("%i",run));
    if( fHmpMipTrkDistPosY[ich]->GetEntries()<100. ){hgHmpYresidPosFromHistoMean[ich] -> SetBinContent(index+1, -999);}
    else { hgHmpYresidPosFromHistoMean[ich] -> SetBinContent(index+1, fHmpMipTrkDistPosY[ich]->GetMean()); }
    hgHmpYresidPosFromHistoMean[ich] -> SetBinError(index+1, 0.000000000000000001);

    // MIP TREND HISTOS.
    hgHmpMipCharge[ich] -> GetXaxis()->SetBinLabel(index+1,Form("%d",run));
    Double_t lanMpv = 0 , lanSigma = 0 ;
    if(fHmpMipCharge[ich]->GetEntries() > 0  && fHmpMipCharge[ich] -> GetFunction("landau"))  {
      lanMpv = fHmpMipCharge[ich] -> GetFunction("landau")->GetParameter(1); //MPV
      lanSigma = fHmpMipCharge[ich] -> GetFunction("landau")->GetParameter(2); //sigma
      hgHmpMipCharge[ich]->SetBinContent(index+1,lanMpv);
      hgHmpMipCharge[ich]->SetBinError(index+1,lanSigma);
    } else {      hgHmpMipCharge[ich]->SetBinContent(index+1,-999); hgHmpMipCharge[ich]->SetBinError(index+1, 0.000000000000000001);}

  }// end RICH loop.
}

//______________________________________________________________________________
void Print_Trending(char * dataType, Int_t year, char *period, char *pass, Int_t numberofruns){

  TCanvas * cmain = new TCanvas("cmain","cmain",1200,800);
  cmain->SaveAs(Form("%s/%s/%i/%s/%s/trending/hmpid_qa_trend.pdf[", workdir, dataType, year, period, pass) );

  TLine *nl1 = new TLine(0.,1.5,numberofruns,1.5); nl1->SetLineStyle(2); nl1->SetLineColor(kRed); nl1->SetLineWidth(1.5); nl1->Draw("same");

  TLine *nl2 = new TLine(0.,-1.5,numberofruns,-1.5); nl2->SetLineStyle(2); nl2->SetLineColor(kRed); nl2->SetLineWidth(1.5); nl2->Draw("same");

  TCanvas* c1 = new TCanvas("c1","c1",1200,800);
  c1->Divide(1,2);

  c1->cd(1);

  TLegend *pLeg=new TLegend(0.85,0.99,0.99,0.8);

  TString psChambers[7]={"RICH0","RICH1","RICH2","RICH3","RICH4","RICH5","RICH6"};

  for(Int_t ich = 0 ; ich < 7; ich++){

    //__ X Residual Negative
    if( ich == 0 ){
      hgHmpXresidNegFromHistoMean[ich]->SetAxisRange(-3,3,"Y");
      hgHmpXresidNegFromHistoMean[ich]->Draw("p");
      pLeg->AddEntry(hgHmpXresidNegFromHistoMean[ich],Form("%s",psChambers[ich].Data()),"p");
    } else {
      hgHmpXresidNegFromHistoMean[ich]->SetAxisRange(-3,3,"Y");
      hgHmpXresidNegFromHistoMean[ich]->Draw("psame");
      pLeg->AddEntry(hgHmpXresidNegFromHistoMean[ich],Form("%s",psChambers[ich].Data()),"p");
    }
  }//end RICH loop

  nl1->Draw("same");

  nl2->Draw("same");

  pLeg->Draw("same");

  c1->cd(2);

  pLeg->Clear();

  for(Int_t ich = 0 ; ich < 7; ich++){

    //__ X Residual Positive
    if( ich == 0 ){
      hgHmpXresidPosFromHistoMean[ich]->SetAxisRange(-3,3,"Y");
      hgHmpXresidPosFromHistoMean[ich]->Draw("p");
      pLeg->AddEntry(hgHmpXresidPosFromHistoMean[ich],Form("%s",psChambers[ich].Data()),"p");
    } else {
      hgHmpXresidPosFromHistoMean[ich]->SetAxisRange(-3,3,"Y");
      hgHmpXresidPosFromHistoMean[ich]->Draw("psame");
      pLeg->AddEntry(hgHmpXresidPosFromHistoMean[ich],Form("%s",psChambers[ich].Data()),"p");
    }
  }//end RICH loop

  nl1->Draw("same");

  nl2->Draw("same");

  pLeg->Draw("same");

  c1->SaveAs(Form("%s/%s/%i/%s/%s/trending/hmpid_qa_trend.pdf", workdir, dataType, year, period, pass) );


  TCanvas* c2 = new TCanvas("c2","c2",1200,800);
  c2->Divide(1,2);

  c2->cd(1);

  pLeg->Clear();

  for(Int_t ich = 0 ; ich < 7; ich++){

    //__ Y Residual Negative
    if( ich == 0 ){
      hgHmpYresidNegFromHistoMean[ich]->SetAxisRange(-3,3,"Y");
      hgHmpYresidNegFromHistoMean[ich]->Draw("p");
      pLeg->AddEntry(hgHmpYresidNegFromHistoMean[ich],Form("%s",psChambers[ich].Data()),"p");
    } else {
      hgHmpYresidNegFromHistoMean[ich]->SetAxisRange(-3,3,"Y");
      hgHmpYresidNegFromHistoMean[ich]->Draw("psame");
      pLeg->AddEntry(hgHmpYresidNegFromHistoMean[ich],Form("%s",psChambers[ich].Data()),"p");
    }
  }//end RICH loop

  nl1->Draw("same");

  nl2->Draw("same");

  pLeg->Draw("same");

  c2->cd(2);
  pLeg->Clear();
  for(Int_t ich = 0 ; ich < 7; ich++){

    //__ Y Residual Positive.
    if( ich == 0 ){
      hgHmpYresidPosFromHistoMean[ich]->SetAxisRange(-3,3,"Y");
      hgHmpYresidPosFromHistoMean[ich]->Draw("p");
      pLeg->AddEntry(hgHmpYresidPosFromHistoMean[ich],Form("%s",psChambers[ich].Data()),"p");
    } else {
      hgHmpYresidPosFromHistoMean[ich]->SetAxisRange(-3,3,"Y");
      hgHmpYresidPosFromHistoMean[ich]->Draw("psame");
      pLeg->AddEntry(hgHmpYresidPosFromHistoMean[ich],Form("%s",psChambers[ich].Data()),"p");
    }
  }//end RICH loop

  nl1->Draw("same");

  nl2->Draw("same");

  pLeg->Draw("same");

  c2->SaveAs(Form("%s/%s/%i/%s/%s/trending/hmpid_qa_trend.pdf", workdir, dataType, year, period, pass) );

  //MIP TREND Histos
  TCanvas* c3 = new TCanvas("c3","c3",1200,800);

  pLeg->Clear();

  for(Int_t ich = 0 ; ich < 7; ich++){// RICH loop
    hgHmpMipCharge[ich]->SetMarkerStyle(24);    hgHmpMipCharge[ich]->SetMarkerColor(ich+1);  hgHmpMipCharge[ich]->SetLineColor(ich+1);
    if( ich == 0 ){
      hgHmpMipCharge[ich]->SetAxisRange(0,1000,"Y");
      hgHmpMipCharge[ich]->Draw();
      pLeg->AddEntry(hgHmpMipCharge[ich],Form("%s",psChambers[ich].Data()),"p");
    } else {
      hgHmpMipCharge[ich]->Draw("same");
      pLeg->AddEntry(hgHmpMipCharge[ich],Form("%s",psChambers[ich].Data()),"p");
    }
  }//end RICH loop

  pLeg->Draw("same");

  c3->SaveAs(Form("%s/%s/%i/%s/%s/trending/hmpid_qa_trend.pdf", workdir, dataType, year, period, pass) );


  cmain->SaveAs(Form("%s/%s/%i/%s/%s/trending/hmpid_qa_trend.pdf]", workdir, dataType, year, period, pass) );

}

//______________________________________________________________________________
void Save_Trending(char * dataType, Int_t year, char *period, char *pass){

  // Write files
  TFile *fout = new TFile(Form("%s/%s/%i/%s/%s/trending/hmpid_qa_trend.root", workdir, dataType, year, period, pass), "RECREATE");

  for(Int_t ich = 0 ; ich < 7; ich++) hgHmpXresidNegFromHistoMean[ich] -> Write();

  for(Int_t ich = 0 ; ich < 7; ich++) hgHmpXresidPosFromHistoMean[ich] -> Write();

  for(Int_t ich = 0 ; ich < 7; ich++) hgHmpYresidNegFromHistoMean[ich] -> Write();

  for(Int_t ich = 0 ; ich < 7; ich++) hgHmpYresidPosFromHistoMean[ich] -> Write();

  for(Int_t ich = 0 ; ich < 7; ich++) hgHmpMipCharge[ich] -> Write();

  fout->Close();

}
