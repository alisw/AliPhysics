#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TCanvas.h>
#include <TGrid.h>
#include <TFile.h>
#include <TList.h>
#include <TGridResult.h>
#include <TPaveStats.h>
#include <TGraph.h>
#include <TMath.h>
#include <TGraphErrors.h>
#include <TH1.h>
#include <TF1.h>
#include <TH2.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TLatex.h>
#include <TString.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TROOT.h>
#endif

void PlotGeneral(TFile* fildat, TFile* filMC, TCanvas**& clist, Int_t& cnum);
void PlotITSsa(TFile* fildat, TFile* filMC, TCanvas**& clist, Int_t& cnum);
void PlotSDD(TFile* fildat, TFile* filMC, TCanvas**& clist, Int_t& cnum);
void GetGainModuleLevelSSD(TFile* fildat, TFile* filMC, TCanvas**& clist, Int_t& cnum);
void VertexQAMacro(TFile *fildat, TFile* filMC, TCanvas **&clist, Int_t &cnum);
void SPDVertexPileup(TFile* fildat, TFile* filMC, TCanvas**& clist, Int_t& cnum);
void PlotSPD(TFile* fildat, TFile* filMC, TCanvas**& clist, Int_t& cnum);
void PlotSPDVtxPileup(TFile* fildat, TCanvas**& clist, Int_t& cnum);
Bool_t PlotITSTPCMatchingEff(TFile *f, TFile* filMC, TCanvas**& clist,Int_t& cnum);
void SetDrawAtt(Int_t markerstyle,Int_t markercolor,Int_t markersize,Int_t linecolor,Int_t linewidth,TH1 *h1);
Double_t LangausFun(Double_t *x, Double_t *par);
void SaveC(TFile &fout, TCanvas**& clist, Int_t cnum);
TString GetRunNumber();
void PlotPID(TFile* fildat,TCanvas**& clist,Int_t& cnum);

//  the run number is available to all the functions. Its value is set by AliITSQAchecks
  Int_t gRunNumber = 0;
  Int_t gRunNumberMC = 0;
  TString pdfFileNames="";


//_______________________________________________________________________
//void AliITSQAchecks_pp(TString option="local",             // local analysis version
//                              Int_t nRun=257912,
//                              TString filenamedata= "QAresults_barrel_257912.root",
//                              TString filenameMC="QAresults_barrel_257912.root",
//                              Int_t nRunMC=257912){

/*void AliITSQAchecks_pp(TString option="grid",                // grid analysis version
                                       Int_t nRun=244540,
//                                       //TString filenamedata="alien:///alice/cern.ch/user/g/germain/data/2016/LHC16c/000251255/quiet_beams_pass1/QAresults.root",
//                                       //TString filenameMC="alien:///alice/cern.ch/user/g/germain/data/2016/LHC16c/000251255/quiet_beams_pass1/QAresults.root",
                                       TString filenamedata="alien:///alice/data/2015/LHC15n/000244540/pass2/QAresults.root",
                                       TString filenameMC="alien:///alice/data/2015/LHC15n/000244540/pass2/QAresults.root",
                                       Int_t nRunMC=245829){
*/
void AliITSQAchecks(TString option="grid",                // grid analysis version
                         Int_t nRun=245952,
                         //                                  TString filenamedata="alien:///alice/data/2015/LHC15o/000245064/lowIR_standaloneITS/QAresults.root",
                         TString filenamedata="alien:///alice/data/2015/LHC15o/000245952/pass1/QAresults.root",
                         //                                  TString filenameMC="../Trending/MC/QAresults_245064.root", // merge parziale
                         TString filenameMC="alien:///alice/sim/2016/LHC16g1a/245952/QAresults.root",
                         Int_t nRunMC=245952){

/*void AliITSQAchecks(TString option="grid",                // grid analysis version
                         Int_t nRun=262425,
                         //                                  TString filenamedata="alien:///alice/data/2015/LHC15o/000245064/lowIR_standaloneITS/QAresults.root",
                         TString filenamedata="alien:///alice/data/2016/LHC16o/000262425/cpass1_pass1/QAresults_barrel.root",
                         //                                  TString filenameMC="../Trending/MC/QAresults_245064.root", // merge parziale
                         TString filenameMC="alien:///alice/data/2016/LHC16o/000262425/cpass1_pass1/QAresults_barrel.root",
                         Int_t nRunMC=262425){
*/
  // THIS MACRO SHOULD BE COMPILED. IT DOES NOT WORK WITH THE INTERPRETER
  // option:  "local" if filenamedata is the name of a local file
  //          "grid" if on alien
  //          the two options are mutually exclusive (no mixed local/Grid combination)
  // nRun: data run number (histos title)
  // filenamedata: total path of the data file on alien in "grid" mode, local file name in "local mode"
  // filenameMC: total path of the sim file on alien in "grid" mode, local file name in "local mode"
  // nRunMC:  sim run number for comparison (histos title)
  // Select here what you want to display
  // the complete selection string is
  // "general ITSSA SPD SDD SSD vertex ITSTPC"
  // Contact:  Elena Botta: botta@to.infn.it

/* $Id: AliITSQAchecks.C 53441 2011-12-06 16:30:40Z masera $ */

  gRunNumber = nRun; // data histos title
  gRunNumberMC = nRunMC; // sim histos title
  
// N.B.: PID option shpuld be used for ppass productions only
//    TString selection("pileup");
    
//    TString selection("general ITSSA SPD SDD SSD vertex ITSTPC pileup PID");
    TString selection("ITSTPC");
//    TString selection("ITSSA SPD vertex pileup");
 
  // TString selection("vertex"); 
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(1111);
    
  TFile *fildat, *filMC;

  if(option.Contains("local")){
    fildat=new TFile(filenamedata.Data());
    printf("Opened file %s\n",fildat->GetName());
    filMC=new TFile(filenameMC.Data());
    printf("Opened file %s\n",filMC->GetName());
    }
  else{
    TGrid::Connect("alien:");
    fildat=TFile::Open(filenamedata.Data());
    printf("Opened file %s\n",fildat->GetName());
    filMC=TFile::Open(filenameMC.Data());
    printf("Opened file %s\n",filMC->GetName());
  }
    

  TCanvas** clist;
  Int_t cnum;
  char rn[10];
//  cout<<gRunNumberMC<<endl;
  sprintf(rn,"%d",gRunNumber);
  TString strRN(rn);
  TString founame="Outfil"+strRN+".root";
  TFile fout(founame,"recreate");
  if(selection.Contains("general")){
    PlotGeneral(fildat,filMC,clist,cnum);
    printf("GENERAL - cnum = %d\n",cnum);
      SaveC(fout,clist,cnum);

  }
  if(selection.Contains("ITSSA")){
    PlotITSsa(fildat,filMC,clist,cnum);
    printf("ITSSA - cnum = %d\n",cnum);
    SaveC(fout,clist,cnum);
  }
  if(selection.Contains("SDD")){
    PlotSDD(fildat,filMC,clist,cnum);
    printf("SDD - cnum = %d\n",cnum);
    SaveC(fout,clist,cnum);
  }
  if(selection.Contains("SSD")){
    GetGainModuleLevelSSD(fildat,filMC,clist,cnum);
    printf("SSD - cnum = %d\n",cnum);
    SaveC(fout,clist,cnum);
  }
  if(selection.Contains("vertex")){
    VertexQAMacro(fildat,filMC,clist,cnum);
    printf("VERTEX - cnum = %d\n",cnum);
    SaveC(fout,clist,cnum);
    SPDVertexPileup(fildat,filMC,clist,cnum);
    printf("Pileup - cnum = %d\n",cnum);
    SaveC(fout,clist,cnum);
  }
  if(selection.Contains("SPD")){
    PlotSPD(fildat,filMC,clist,cnum);
    printf("SPD - cnum = %d\n",cnum);
    SaveC(fout,clist,cnum);
  }
    if(selection.Contains("ITSTPC")){
        PlotITSTPCMatchingEff(fildat,filMC,clist,cnum);
        printf("ITSTPC - cnum = %d\n",cnum);
        SaveC(fout,clist,cnum);
    }
    if(selection.Contains("pileup")){
        PlotSPDVtxPileup(fildat,clist,cnum);
        printf("SPD vtx pileup - cnum = %d\n",cnum);
        if(cnum>0) SaveC(fout,clist,cnum);
    }
    if(selection.Contains("PID")){
        PlotPID(fildat,clist,cnum);
        printf("PID - cnum = %d\n",cnum);
        SaveC(fout,clist,cnum);
    }

  fout.Close();

  // merge the pdf files
  TString command("gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=merged");
  command=command+strRN+".pdf "+pdfFileNames;
  gSystem->Exec(command.Data());
  printf(" Merging the pdf file:  %s \n",command.Data());
  
}

//_______________________________________________________________________
void PlotGeneral(TFile* fildat, TFile* filMC, TCanvas**& clist, Int_t& cnum){
  TDirectoryFile* df=(TDirectoryFile*)fildat->Get("SDD_Performance");
  if(!df){
    printf("SDD_Performance data MISSING -> Exit\n");
    return;
  }
  TList* l=(TList*)df->Get("coutputRP");
  if(!l){
    printf("coutputRP TList data MISSING -> Exit\n");
    return;
  }
  cnum=1; // number of canvases 
  clist= new TCanvas* [1];//array of pointers to TCanvases
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(1111);
  TH1F* hcllay=(TH1F*)l->FindObject("hCluInLay");
  TH1F* hev=(TH1F*)l->FindObject("hNEvents");
  Int_t nTotEvents=hev->GetBinContent(1);
  Int_t nTrigEvents=hev->GetBinContent(2);
  Int_t nEvents=nTotEvents;
  printf("---- Statistics ----\n");
  printf("Number of Events = %d\n",nTotEvents);
  if(nTrigEvents>0){ 
    printf("Number of Triggered Events = %d\n",nTrigEvents);
    nEvents=nTrigEvents;
  }else{
    printf("No request on the trigger done when running the task\n");
  }
  if(hcllay){
    Double_t norm=hcllay->GetBinContent(1);
    if(norm>0.){
      hcllay->Scale(1./norm);
      hcllay->SetTitle("");
      hcllay->GetXaxis()->SetRange(2,7);
      hcllay->SetMinimum(0.);
      hcllay->SetMaximum(1.1);
      hcllay->SetMarkerStyle(23);
      TString ctitle=GetRunNumber()+"General checks: PointPerLayer - DATA";
      TCanvas* ceffL=new TCanvas("ceffL",ctitle,1000,800);
      clist[0]=ceffL;
      // ceffL->Divide(1,2);
      // ceffL->cd(1);
      ceffL->SetGridy();
      hcllay->Draw(); 
      TLatex* tg=new TLatex(0.15,0.2,"Fraction of tracks with point in ITS layer");
      tg->SetTextSize(0.04);
      tg->SetNDC();
      tg->SetTextColor(1);
      tg->Draw();
      TString testo="Run "+GetRunNumber()+" -  DATA";
      TLatex* tg2 = new TLatex(0.15,0.85,testo.Data());
      tg2->SetTextSize(0.04);
      tg2->SetNDC();
      tg2->SetTextColor(2);
      tg2->Draw();
      hcllay->GetXaxis()->SetTitle("Layer");
      hcllay->GetYaxis()->SetTitle("Fraction of tracks with point in layer");
      ceffL->Update();
      ceffL->SaveAs("track_points_per_layer.pdf");
      pdfFileNames+=" track_points_per_layer.pdf";
    }
  }
// MC file
    TDirectoryFile* dfMC=(TDirectoryFile*)filMC->Get("SDD_Performance");
    if(!dfMC){
        printf("SDD_Performance MC MISSING -> Exit\n");
        return;
    }
    TList* lMC=(TList*)dfMC->Get("coutputRP");
    if(!lMC){
        printf("coutputRP TList MC MISSING -> Exit\n");
        return;
    }
    cnum=1; // number of canvases
    clist= new TCanvas* [1];//array of pointers to TCanvases
    gROOT->SetStyle("Plain");
    gStyle->SetOptStat(1111);
    TH1F* hcllayMC=(TH1F*)lMC->FindObject("hCluInLay");
    TH1F* hevMC=(TH1F*)lMC->FindObject("hNEvents");
    Int_t nTotEventsMC=hevMC->GetBinContent(1);
    Int_t nTrigEventsMC=hevMC->GetBinContent(2);
    Int_t nEventsMC=nTotEventsMC;
    printf("---- Statistics ----\n");
    printf("Number of MC Events = %d\n",nTotEventsMC);
    if(nTrigEvents>0){
        printf("Number of MC Triggered Events = %d\n",nTrigEventsMC);
        nEventsMC=nTrigEventsMC;
    }else{
        printf("No request on the trigger done when running the task\n");
    }
    if(hcllayMC){
        Double_t normMC=hcllayMC->GetBinContent(1);
        if(normMC>0.){
            hcllayMC->Scale(1./normMC);
            hcllayMC->SetTitle("");
            hcllayMC->GetXaxis()->SetRange(2,7);
            hcllayMC->SetMinimum(0.);
            hcllayMC->SetMaximum(1.1);
            hcllayMC->SetMarkerStyle(23);
            TString ctitleMC=GetRunNumber()+"General checks: PointPerLayer - MC";
            TCanvas* ceffLMC=new TCanvas("ceffLMC",ctitleMC,1000,800);
            clist[0]=ceffLMC;
            // ceffL->Divide(1,2);
            // ceffL->cd(1);
            ceffLMC->SetGridy();
            hcllayMC->Draw();
            TLatex* tgMC=new TLatex(0.15,0.2,"Fraction of tracks with point in ITS layer");
            tgMC->SetTextSize(0.04);
            tgMC->SetNDC();
            tgMC->SetTextColor(1);
            tgMC->Draw();
            TString testoMC="Run "+GetRunNumber()+" -  MC";
            TLatex* tg2MC = new TLatex(0.15,0.85,testoMC.Data());
            tg2MC->SetTextSize(0.04);
            tg2MC->SetNDC();
            tg2MC->SetTextColor(2);
            tg2MC->Draw();
            hcllayMC->GetXaxis()->SetTitle("Layer");
            hcllayMC->GetYaxis()->SetTitle("Fraction of tracks with point in layer");
            ceffLMC->Update();
            ceffLMC->SaveAs("track_points_per_layer_MC.pdf");
            pdfFileNames+=" track_points_per_layer_MC.pdf";
        }
    }

}


//_______________________________________________________________________
//////////////////////////////////////////////////////////////////////
/// ITSsa ////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
 void PlotITSsa(TFile* fildat, TFile* filMC, TCanvas**& clist, Int_t& cnum){
    TDirectoryFile* df=(TDirectoryFile*)fildat->Get("TracksITSsa");
    if(!df) df=(TDirectoryFile*)fildat->Get("ITSsaTracks"); // <<------
    if(!df){
      printf("ITSsa_Performance data MISSING -> Exit\n");
      return;
    }
 
    TList* l=(TList*)df->Get("clistITSsaTracks");
    if(!df){
      printf("clistITSsaTracks TList data MISSING -> Exit\n");
      return;
    }
    cnum=2+2+3*2; // number of canvases (+2 per i plot 2-d relativi al trigger e a tracks vs tracklets; +3*2 per i plot etaphi perlayer)
    clist= new TCanvas* [cnum];//array of pointers to TCanvases
     gROOT->SetStyle("Plain");
  gStyle->SetOptStat(1111);


  TH1F* hPtTPCITS=(TH1F*)l->FindObject("hPtTPCITS");
  TH1F* hPtITSsa=(TH1F*)l->FindObject("hPtITSsa");
  TH1F* hPtITSpureSA=(TH1F*)l->FindObject("hPtITSpureSA");

  TH2F* hEtaPhiTPCITS=(TH2F*)l->FindObject("hEtaPhiTPCITS");
  TH2F* hEtaPhiITSsa=(TH2F*)l->FindObject("hEtaPhiITSsa");
  TH2F* hEtaPhiITSpureSA=(TH2F*)l->FindObject("hEtaPhiITSpureSA");
  TH1F* hChi2TPCITS=(TH1F*)l->FindObject("hChi2TPCITS");
  TH1F* hChi2ITSsa=(TH1F*)l->FindObject("hChi2ITSsa");
  TH1F* hChi2ITSpureSA=(TH1F*)l->FindObject("hChi2ITSpureSA");

// Elena: gestire il check sull'esistenza degli histos
     TH1F* hNEvents=(TH1F*)l->FindObject("hNEvents");
     TH2F* hNEventsVsTrig=(TH2F*)l->FindObject("hNEventsVsTrig");
     TH2F* hPureSAtracksVsTracklets=(TH2F*)l->FindObject("hPureSAtracksVsTracklets");
     TH2F* hITSTPCtracksVsTracklets=(TH2F*)l->FindObject("hITSTPCtracksVsTracklets");
     TH2F* hEtaPhiTracksLay1TPCITS=(TH2F*)l->FindObject("hEtaPhiTracksLay1TPCITS");
     TH2F* hEtaPhiTracksLay2TPCITS=(TH2F*)l->FindObject("hEtaPhiTracksLay2TPCITS");
     TH2F* hEtaPhiTracksLay3TPCITS=(TH2F*)l->FindObject("hEtaPhiTracksLay3TPCITS");
     TH2F* hEtaPhiTracksLay4TPCITS=(TH2F*)l->FindObject("hEtaPhiTracksLay4TPCITS");
     TH2F* hEtaPhiTracksLay5TPCITS=(TH2F*)l->FindObject("hEtaPhiTracksLay5TPCITS");
     TH2F* hEtaPhiTracksLay6TPCITS=(TH2F*)l->FindObject("hEtaPhiTracksLay6TPCITS");
     TH2F* hEtaPhiTracksLay1ITSpureSA=(TH2F*)l->FindObject("hEtaPhiTracksLay1ITSpureSA");
     TH2F* hEtaPhiTracksLay2ITSpureSA=(TH2F*)l->FindObject("hEtaPhiTracksLay2ITSpureSA");
     TH2F* hEtaPhiTracksLay3ITSpureSA=(TH2F*)l->FindObject("hEtaPhiTracksLay3ITSpureSA");
     TH2F* hEtaPhiTracksLay4ITSpureSA=(TH2F*)l->FindObject("hEtaPhiTracksLay4ITSpureSA");
     TH2F* hEtaPhiTracksLay5ITSpureSA=(TH2F*)l->FindObject("hEtaPhiTracksLay5ITSpureSA");
     TH2F* hEtaPhiTracksLay6ITSpureSA=(TH2F*)l->FindObject("hEtaPhiTracksLay6ITSpureSA");
//
     
     // fraction of pureSA tracks with hits in ITS layers - DATA
     TH2F* hNcluPion2d = (TH2F*)l->FindObject("hCluInLayITSpureSAPion");
     TH1D* hNcluPion = hNcluPion2d->ProjectionY();
     Float_t fracTpi[6], efracTpi[6];
     for(Int_t iLay=0; iLay<6; iLay++){
           fracTpi[iLay]=hNcluPion->GetBinContent(iLay+2)/hNcluPion->GetBinContent(1);
           efracTpi[iLay]=TMath::Sqrt(fracTpi[iLay]*(1-fracTpi[iLay])/hNcluPion->GetBinContent(1));
     }

     TH1F* histoFracpureSA=new TH1F("histoFracpureSA","",7,0.5,7.5);
     histoFracpureSA->GetYaxis()->SetRangeUser(0.,1.1);
     histoFracpureSA->SetMarkerStyle(20);
     histoFracpureSA->GetXaxis()->SetTitle("Layer");
     histoFracpureSA->GetYaxis()->SetTitle("Fraction of pureSA pion tracks w/ points in layer");
     for(Int_t iLay=0; iLay<6; iLay++){
         histoFracpureSA->SetBinContent(iLay+1,fracTpi[iLay]);
         histoFracpureSA->SetBinError(iLay+1,efracTpi[iLay]);
     }
     
     
     // fraction of pureSA tracks with hits in ITS layers - DATA
     
  TH1F* hRatio=(TH1F*)hPtTPCITS->Clone("hRatio");
  TH1F* hRatio1=(TH1F*)hPtTPCITS->Clone("hRatio1");
   hRatio->Add(hPtITSsa);
  hRatio->Divide(hPtITSpureSA);
  hRatio->SetStats(0);
  hRatio1->Divide(hPtITSsa);
  hRatio1->SetStats(0);

     // partire da qui: gestire il clist[2]-->clist[4] con divisione 2,1, plottare, salvare su file pdf singolo e globale
     // partire da qui: gestire il clist[4]-->clist[10] con divisione 2,1, plottare, salvare su file pdf singolo e globale

     TString ctitle=GetRunNumber()+"ITS standalone: ITS configuration and trigger mask - DATA";
     TCanvas* cITSsa0=new TCanvas("cITSsa0",ctitle,1200,1200);
     clist[0]=cITSsa0;
     cITSsa0->Divide(1,2);
     cITSsa0->cd(1);
     hNEvents->Draw();
     cITSsa0->cd(2);
     if(hNEventsVsTrig)hNEventsVsTrig->Draw();
     cITSsa0->Update();
     cITSsa0->SaveAs("ITSsa0.pdf");
     pdfFileNames+=" ITSsa0.pdf";
     
     ctitle=GetRunNumber()+"ITS standalone: tracks vs tracklets - DATA";
     TCanvas* cITSsa0a=new TCanvas("cITSsa0a",ctitle,1200,1200);
     clist[1]=cITSsa0a;
     cITSsa0a->Divide(1,2);
     cITSsa0a->cd(1);
     if(hPureSAtracksVsTracklets){hPureSAtracksVsTracklets->GetXaxis()->SetTitle("nTracklets");hPureSAtracksVsTracklets->Draw("colz");}
     cITSsa0a->cd(2);
     if(hITSTPCtracksVsTracklets){hITSTPCtracksVsTracklets->GetXaxis()->SetTitle("nTracklets");hITSTPCtracksVsTracklets->Draw("colz");}
     cITSsa0a->Update();
     cITSsa0a->SaveAs("ITSsa0a.pdf");
     pdfFileNames+=" ITSsa0a.pdf";

// new
     TCanvas* cITSsa0b;
     ctitle=GetRunNumber()+"ITS standalone: eta-phi distribution for TPCITS tracks, SPD - DATA";
     if(hEtaPhiTracksLay1TPCITS){
     cITSsa0b=new TCanvas("cITSsa0b",ctitle,1200,1000);
     clist[2]=cITSsa0b;
     cITSsa0b->Divide(1,2);
     cITSsa0b->SetLogz(); // prova per vedere meglio i protoni, test sui pioni
     cITSsa0b->cd(1);
     if(hEtaPhiTracksLay1TPCITS){hEtaPhiTracksLay1TPCITS->GetXaxis()->SetTitle("Eta");hEtaPhiTracksLay1TPCITS->GetYaxis()->SetTitle("Phi");hEtaPhiTracksLay1TPCITS->Draw("colz");}
     cITSsa0b->cd(2);
     if(hEtaPhiTracksLay2TPCITS){hEtaPhiTracksLay2TPCITS->GetXaxis()->SetTitle("Eta");hEtaPhiTracksLay2TPCITS->GetYaxis()->SetTitle("Phi");hEtaPhiTracksLay2TPCITS->Draw("colz");}
     cITSsa0b->Update();
     cITSsa0b->SaveAs("TPCITS_etaphi_SPD.pdf");
     pdfFileNames+=" TPCITS_etaphi_SPD.pdf";
     }
     
     TCanvas* cITSsa0c;
     ctitle=GetRunNumber()+"ITS standalone: eta-phi distribution for TPCITS tracks, SDD - DATA";
     if(hEtaPhiTracksLay3TPCITS){
     cITSsa0c=new TCanvas("cITSsa0c",ctitle,1200,1000);
     clist[3]=cITSsa0c;
     cITSsa0c->Divide(1,2);
     cITSsa0c->cd(1);
     if(hEtaPhiTracksLay3TPCITS){hEtaPhiTracksLay3TPCITS->GetXaxis()->SetTitle("Eta");hEtaPhiTracksLay3TPCITS->GetYaxis()->SetTitle("Phi");hEtaPhiTracksLay3TPCITS->Draw("colz");}
     cITSsa0c->cd(2);
     if(hEtaPhiTracksLay4TPCITS){hEtaPhiTracksLay4TPCITS->GetXaxis()->SetTitle("Eta");hEtaPhiTracksLay4TPCITS->GetYaxis()->SetTitle("Phi");hEtaPhiTracksLay4TPCITS->Draw("colz");}
     cITSsa0c->Update();
     cITSsa0c->SaveAs("TPCITS_etaphi_SDD.pdf");
     pdfFileNames+=" TPCITS_etaphi_SDD.pdf";
     }
     
     TCanvas* cITSsa0d;
     ctitle=GetRunNumber()+"ITS standalone: eta-phi distribution for TPCITS tracks, SSD - DATA";
     if(hEtaPhiTracksLay5TPCITS){
     cITSsa0d=new TCanvas("cITSsa0d",ctitle,1200,1000);
     clist[4]=cITSsa0d;
     cITSsa0d->Divide(1,2);
     cITSsa0d->cd(1);
     if(hEtaPhiTracksLay5TPCITS){hEtaPhiTracksLay5TPCITS->GetXaxis()->SetTitle("Eta");hEtaPhiTracksLay5TPCITS->GetYaxis()->SetTitle("Phi");hEtaPhiTracksLay5TPCITS->Draw("colz");}
     cITSsa0d->cd(2);
     if(hEtaPhiTracksLay6TPCITS){hEtaPhiTracksLay6TPCITS->GetXaxis()->SetTitle("Eta");hEtaPhiTracksLay6TPCITS->GetYaxis()->SetTitle("Phi");hEtaPhiTracksLay6TPCITS->Draw("colz");}
     cITSsa0d->Update();
     cITSsa0d->SaveAs("TPCITS_etaphi_SSD.pdf");
     pdfFileNames+=" TPCITS_etaphi_SSD.pdf";
     }
     
     TCanvas* cITSsa0e;
     ctitle=GetRunNumber()+"ITS standalone: eta-phi distribution for ITSpureSA tracks, SPD - DATA";
     if(hEtaPhiTracksLay1ITSpureSA){
     cITSsa0e=new TCanvas("cITSsa0e",ctitle,1200,1000);
     clist[5]=cITSsa0e;
     cITSsa0e->Divide(1,2);
     cITSsa0e->cd(1);
     if(hEtaPhiTracksLay1ITSpureSA){hEtaPhiTracksLay1ITSpureSA->GetXaxis()->SetTitle("Eta");hEtaPhiTracksLay1ITSpureSA->GetYaxis()->SetTitle("Phi");hEtaPhiTracksLay1ITSpureSA->Draw("colz");}
     cITSsa0e->cd(2);
     if(hEtaPhiTracksLay2ITSpureSA){hEtaPhiTracksLay2ITSpureSA->GetXaxis()->SetTitle("Eta");hEtaPhiTracksLay2ITSpureSA->GetYaxis()->SetTitle("Phi");hEtaPhiTracksLay2ITSpureSA->Draw("colz");}
     cITSsa0e->Update();
     cITSsa0e->SaveAs("ITSpureSA_etaphi_SPD.pdf");
     pdfFileNames+=" ITSpureSA_etaphi_SPD.pdf";
     }

     TCanvas* cITSsa0f;
     ctitle=GetRunNumber()+"ITS standalone: eta-phi distribution for ITSpureSA tracks, SDD - DATA";
     if(hEtaPhiTracksLay3ITSpureSA){
     cITSsa0f=new TCanvas("cITSsa0f",ctitle,1200,1000);
     clist[6]=cITSsa0f;
     cITSsa0f->Divide(1,2);
     cITSsa0f->cd(1);
     if(hEtaPhiTracksLay3ITSpureSA){hEtaPhiTracksLay3ITSpureSA->GetXaxis()->SetTitle("Eta");hEtaPhiTracksLay3ITSpureSA->GetYaxis()->SetTitle("Phi");hEtaPhiTracksLay3ITSpureSA->Draw("colz");}
     cITSsa0f->cd(2);
     if(hEtaPhiTracksLay4ITSpureSA){hEtaPhiTracksLay4ITSpureSA->GetXaxis()->SetTitle("Eta");hEtaPhiTracksLay4ITSpureSA->GetYaxis()->SetTitle("Phi");hEtaPhiTracksLay4ITSpureSA->Draw("colz");}
     cITSsa0f->Update();
     cITSsa0f->SaveAs("ITSpureSA_etaphi_SDD.pdf");
     pdfFileNames+=" ITSpureSA_etaphi_SDD.pdf";
     }
     
     TCanvas* cITSsa0g;
     ctitle=GetRunNumber()+"ITS standalone: eta-phi distribution for ITSpureSA tracks, SSD - DATA";
     if(hEtaPhiTracksLay5ITSpureSA){
     cITSsa0g=new TCanvas("cITSsa0g",ctitle,1200,1000);
     clist[7]=cITSsa0g;
     cITSsa0g->Divide(1,2);
     cITSsa0g->cd(1);
     if(hEtaPhiTracksLay5ITSpureSA){hEtaPhiTracksLay5ITSpureSA->GetXaxis()->SetTitle("Eta");hEtaPhiTracksLay5ITSpureSA->GetYaxis()->SetTitle("Phi");hEtaPhiTracksLay5ITSpureSA->Draw("colz");}
     cITSsa0g->cd(2);
     if(hEtaPhiTracksLay6ITSpureSA){hEtaPhiTracksLay6ITSpureSA->GetXaxis()->SetTitle("Eta");hEtaPhiTracksLay6ITSpureSA->GetYaxis()->SetTitle("Phi");hEtaPhiTracksLay6ITSpureSA->Draw("colz");}
     cITSsa0g->Update();
     cITSsa0g->SaveAs("ITSpureSA_etaphi_SSD.pdf");
     pdfFileNames+=" ITSpureSA_etaphi_SSD.pdf";
     }
     
     // new

     ctitle=GetRunNumber()+"ITS standalone: performance vs Pt - DATA";
  TCanvas* cITSsa1=new TCanvas("cITSsa1",ctitle,1200,1200);
  clist[8]=cITSsa1;
  cITSsa1->Divide(1,3);
  cITSsa1->cd(1);
  hPtITSpureSA->Draw();
  hPtITSpureSA->GetXaxis()->SetTitle("Pt (GeV/c)");
  gPad->Update();
  TPaveStats *st1=(TPaveStats*)hPtITSpureSA->GetListOfFunctions()->FindObject("stats");
  st1->SetY1NDC(0.71);
  st1->SetY2NDC(0.9);
  hPtTPCITS->SetLineColor(2);
  hPtTPCITS->GetXaxis()->SetTitle("Pt (GeV/c)");
  hPtTPCITS->Draw("sames");
  //hPtTPCITS->Draw();
  gPad->Update();
  TPaveStats *st2=(TPaveStats*)hPtTPCITS->GetListOfFunctions()->FindObject("stats");
  st2->SetY1NDC(0.71);
  st2->SetY2NDC(0.9);
  st2->SetTextColor(2);

  hPtITSsa->SetLineColor(4);
  hPtITSsa->Draw("sames");
  gPad->Update();
  TPaveStats *st3=(TPaveStats*)hPtITSsa->GetListOfFunctions()->FindObject("stats");
  st3->SetY1NDC(0.51);
  st3->SetY2NDC(0.7);
  st3->SetTextColor(4);
  TLegend* leg=new TLegend(0.5,0.5,0.69,0.79);
  leg->SetFillColor(0);
     TLatex* ttype=new TLatex(0.35,0.75,"DATA");
     ttype->SetNDC();
     ttype->SetTextColor(2);
     ttype->SetTextSize(0.1);
     ttype->Draw();
  TLegendEntry* ent=leg->AddEntry(hPtTPCITS,"TPC+ITS","L");
  ent->SetTextColor(hPtTPCITS->GetLineColor());
  ent=leg->AddEntry(hPtITSsa,"ITSsa","L");
  ent->SetTextColor(hPtITSsa->GetLineColor());
   // to be used only with pp data (ITS pure SA)  
  ent=leg->AddEntry(hPtITSpureSA,"ITS pureSA","L");
  ent->SetTextColor(hPtITSpureSA->GetLineColor());
  leg->Draw();
  cITSsa1->cd(2);
   gPad->SetGridx();
   gPad->SetGridy();
   hRatio1->GetXaxis()->SetTitle("Pt (GeV/c)");
   hRatio->GetXaxis()->SetTitle("Pt (GeV/c)");
   hRatio1->GetYaxis()->SetTitle("TPCITS/ITSsa");
   hRatio->GetYaxis()->SetTitle("(TPCITS+ITSsa)/ITSpureSA");
   //   hRatio1->DrawCopy();
   hRatio->DrawCopy();
   //   TLatex* tratio=new TLatex(0.2,0.75,"TPC+ITS/ITSsa vs Pt");
   TLatex* tratio=new TLatex(0.2,0.75,"TPCITS+ITSsa/ITSpureSA vs Pt");
   tratio->SetNDC();
   tratio->SetTextColor(1);
  tratio->Draw();
  cITSsa1->cd(3);
  hChi2ITSpureSA->Scale(1./hChi2ITSpureSA->GetEntries());
  hChi2ITSsa->Scale(1./hChi2ITSsa->GetEntries());
  hChi2TPCITS->Scale(1./hChi2TPCITS->GetEntries());
  hChi2ITSpureSA->SetLineColor(1);
  hChi2ITSpureSA->Draw("");
  hChi2TPCITS->SetLineColor(2);
  hChi2TPCITS->Draw("sames");
  TLatex* tchi=new TLatex(0.25,0.85,"chi2 vs Pt");
  tchi->SetNDC();
  tchi->SetTextColor(1);
  tchi->Draw();
  gPad->Update();
  TPaveStats *stc2=(TPaveStats*)hChi2TPCITS->GetListOfFunctions()->FindObject("stats");
  stc2->SetY1NDC(0.71);
  stc2->SetY2NDC(0.9);
  stc2->SetTextColor(2);
  //  c2->Update();
  hChi2ITSsa->SetLineColor(4);
  hChi2ITSsa->Draw("sames");
  gPad->Update();
  TPaveStats *stc3=(TPaveStats*)hChi2ITSsa->GetListOfFunctions()->FindObject("stats");
  stc3->SetY1NDC(0.51);
  stc3->SetY2NDC(0.7);
  stc3->SetTextColor(4);
  leg->Draw();

  cITSsa1->Update();
  cITSsa1->SaveAs("ITSsa1.pdf");
  pdfFileNames+=" ITSsa1.pdf";
     
  gStyle->SetPalette(1);
  hEtaPhiITSpureSA->SetStats(0);
  hEtaPhiITSpureSA->SetTitle("ITS pureSA - DATA");
  hEtaPhiITSsa->SetStats(0);
  hEtaPhiITSsa->SetTitle("ITSsa tracks - DATA");
  hEtaPhiTPCITS->SetStats(0);
  hEtaPhiTPCITS->SetTitle("TPC+ITS tracks - DATA");
  ctitle=GetRunNumber()+"Eta-phi distribution for ITSsa and TPC+ITS tracks - DATA";
  TCanvas* cITSsa2=new TCanvas("cITSsa2",ctitle,1200,800);
  clist[9]=cITSsa2;
  cITSsa2->Divide(3,1); //for ITSpuresa
  //cITSsa2->Divide(2,1); //for no ITSpureSA
  cITSsa2->cd(1);
  //TPad* p1=new TPad("p1","Tracking: tracks distribution in EtaPhi",0,0.5,1.,1.);
  //p1->Divide(3,1);
  //p1->cd (1);
  hEtaPhiITSpureSA->Draw("colz");
  hEtaPhiITSpureSA->GetXaxis()->SetTitle("Eta");
  hEtaPhiITSpureSA->GetYaxis()->SetTitle("Phi");
  cITSsa2->cd(2);
  //p1->cd(2);
  hEtaPhiITSsa->Draw("colz");
  hEtaPhiITSsa->GetXaxis()->SetTitle("Eta");
  hEtaPhiITSsa->GetYaxis()->SetTitle("Phi");
   cITSsa2->cd(3);
   //p1->cd(3);
  hEtaPhiTPCITS->Draw("colz");
  hEtaPhiTPCITS->GetXaxis()->SetTitle("Eta");
  hEtaPhiTPCITS->GetYaxis()->SetTitle("Phi");
  //c4->cd(4);
  cITSsa2->SaveAs("ITSsa2.pdf");  
  pdfFileNames+=" ITSsa2.pdf";

// MC file
     TDirectoryFile* dfMC=(TDirectoryFile*)filMC->Get("TracksITSsa");
     if(!dfMC) dfMC=(TDirectoryFile*)filMC->Get("ITSsaTracks");
     if(!dfMC){
         printf("ITSsa_Performance MC MISSING -> Exit\n");
         return;
     }
     
     TList* lMC=(TList*)dfMC->Get("clistITSsaTracks");
     if(!dfMC){
         printf("clistITSsaTracks TList MC MISSING -> Exit\n");
         return;
     }
     cnum=2; // number of canvases
     clist= new TCanvas* [2+1];//array of pointers to TCanvases +1 per fraction of pion pureSA tracks with hit in ITS layers
     gROOT->SetStyle("Plain");
     gStyle->SetOptStat(1111);
     
     
     TH1F* hPtTPCITSMC=(TH1F*)lMC->FindObject("hPtTPCITS");
     TH1F* hPtITSsaMC=(TH1F*)lMC->FindObject("hPtITSsa");
     TH1F* hPtITSpureSAMC=(TH1F*)lMC->FindObject("hPtITSpureSA");
     
     TH2F* hEtaPhiTPCITSMC=(TH2F*)lMC->FindObject("hEtaPhiTPCITS");
     TH2F* hEtaPhiITSsaMC=(TH2F*)lMC->FindObject("hEtaPhiITSsa");
     TH2F* hEtaPhiITSpureSAMC=(TH2F*)lMC->FindObject("hEtaPhiITSpureSA");
     TH1F* hChi2TPCITSMC=(TH1F*)lMC->FindObject("hChi2TPCITS");
     TH1F* hChi2ITSsaMC=(TH1F*)lMC->FindObject("hChi2ITSsa");
     TH1F* hChi2ITSpureSAMC=(TH1F*)lMC->FindObject("hChi2ITSpureSA");
     
     TH1F* hRatioMC=(TH1F*)hPtTPCITSMC->Clone("hRatio");
     TH1F* hRatio1MC=(TH1F*)hPtTPCITSMC->Clone("hRatio1");
     hRatioMC->Add(hPtITSsaMC);
     hRatioMC->Divide(hPtITSpureSAMC);
     hRatioMC->SetStats(0);
     hRatio1MC->Divide(hPtITSsaMC);
     hRatio1MC->SetStats(0);
     
     // fraction of pureSA tracks with hits in ITS layers - MC
     TH2F* hNcluPion2dMC = (TH2F*)lMC->FindObject("hCluInLayITSpureSAPion");
     TH1D* hNcluPionMC = hNcluPion2dMC->ProjectionY();
     Float_t fracTpiMC[6], efracTpiMC[6];
     for(Int_t iLay=0; iLay<6; iLay++){
         fracTpiMC[iLay]=hNcluPionMC->GetBinContent(iLay+2)/hNcluPionMC->GetBinContent(1);
         efracTpiMC[iLay]=TMath::Sqrt(fracTpiMC[iLay]*(1-fracTpiMC[iLay])/hNcluPionMC->GetBinContent(1));
     }
     
     TH1F* histoFracpureSAMC=new TH1F("histoFracpureSAMC","",7,0.5,7.5);
     histoFracpureSAMC->GetYaxis()->SetRangeUser(0.,1.1);
     histoFracpureSAMC->SetMarkerStyle(20);
     histoFracpureSAMC->GetXaxis()->SetTitle("Layer");
     histoFracpureSAMC->GetYaxis()->SetTitle("Fraction of pureSA pion tracks w/ points in layer");

     for(Int_t iLay=0; iLay<6; iLay++){
         histoFracpureSAMC->SetBinContent(iLay+1,fracTpiMC[iLay]);
         histoFracpureSAMC->SetBinError(iLay+1,efracTpiMC[iLay]);
     }
     
     // fraction of pureSA tracks with hits in ITS layers - MC

     
     TString ctitleMC=GetRunNumber()+"ITS standalone: performance vs Pt - MC ";
     TCanvas* cITSsa1MC=new TCanvas("cITSsa1MC",ctitleMC,1200,1200);
     clist[0]=cITSsa1MC;
     cITSsa1MC->Divide(1,3);
     cITSsa1MC->cd(1);
     hPtITSpureSAMC->Draw();
     hPtITSpureSAMC->GetXaxis()->SetTitle("Pt (GeV/c)");
     gPad->Update();
     TPaveStats *st1MC=(TPaveStats*)hPtITSpureSAMC->GetListOfFunctions()->FindObject("stats");
     st1MC->SetY1NDC(0.71);
     st1MC->SetY2NDC(0.9);
     hPtTPCITSMC->SetLineColor(2);
     hPtTPCITSMC->GetXaxis()->SetTitle("Pt (GeV/c)");
     hPtTPCITSMC->Draw("sames");
     //hPtTPCITSMC->Draw();
     gPad->Update();
     TPaveStats *st2MC=(TPaveStats*)hPtTPCITSMC->GetListOfFunctions()->FindObject("stats");
     st2MC->SetY1NDC(0.71);
     st2MC->SetY2NDC(0.9);
     st2MC->SetTextColor(2);
     
     hPtITSsaMC->SetLineColor(4);
     hPtITSsaMC->Draw("sames");
     gPad->Update();
     TPaveStats *st3MC=(TPaveStats*)hPtITSsaMC->GetListOfFunctions()->FindObject("stats");
     st3MC->SetY1NDC(0.51);
     st3MC->SetY2NDC(0.7);
     st3MC->SetTextColor(4);
     TLegend* legMC=new TLegend(0.5,0.5,0.69,0.79);
     legMC->SetFillColor(0);
     TLatex* ttypeMC=new TLatex(0.35,0.75,"MC");
     ttypeMC->SetNDC();
     ttypeMC->SetTextColor(2);
     ttypeMC->SetTextSize(0.1);
     ttypeMC->Draw();
     TLegendEntry* entMC=legMC->AddEntry(hPtTPCITSMC,"TPC+ITS","L");
     entMC->SetTextColor(hPtTPCITSMC->GetLineColor());
     entMC=legMC->AddEntry(hPtITSsaMC,"ITSsa","L");
     entMC->SetTextColor(hPtITSsaMC->GetLineColor());
     // to be used only with pp data (ITS pure SA)
     entMC=legMC->AddEntry(hPtITSpureSAMC,"ITS pureSA","L");
     entMC->SetTextColor(hPtITSpureSAMC->GetLineColor());
     legMC->Draw();
     cITSsa1MC->cd(2);
     gPad->SetGridx();
     gPad->SetGridy();
     hRatio1MC->GetXaxis()->SetTitle("Pt (GeV/c)");
     hRatioMC->GetXaxis()->SetTitle("Pt (GeV/c)");
     hRatio1MC->GetYaxis()->SetTitle("TPCITS/ITSsa");
     hRatioMC->GetYaxis()->SetTitle("(TPCITS+ITSsa)/ITSpureSA");
     //   hRatio1->DrawCopy();
     hRatioMC->DrawCopy();
     //   TLatex* tratio=new TLatex(0.2,0.75,"TPC+ITS/ITSsa vs Pt");
     TLatex* tratioMC=new TLatex(0.2,0.75,"TPCITS+ITSsa/ITSpureSA vs Pt");
     tratioMC->SetNDC();
     tratioMC->SetTextColor(1);
     tratioMC->Draw();
     cITSsa1MC->cd(3);
     hChi2ITSpureSAMC->Scale(1./hChi2ITSpureSAMC->GetEntries());
     hChi2ITSsaMC->Scale(1./hChi2ITSsaMC->GetEntries());
     hChi2TPCITSMC->Scale(1./hChi2TPCITSMC->GetEntries());
     hChi2ITSpureSAMC->SetLineColor(1);
     hChi2ITSpureSAMC->Draw("");
     hChi2TPCITSMC->SetLineColor(2);
     hChi2TPCITSMC->Draw("sames");
     TLatex* tchiMC=new TLatex(0.25,0.85,"chi2 vs Pt");
     tchiMC->SetNDC();
     tchiMC->SetTextColor(1);
     tchiMC->Draw();
     gPad->Update();
     TPaveStats *stc2MC=(TPaveStats*)hChi2TPCITSMC->GetListOfFunctions()->FindObject("stats");
     stc2MC->SetY1NDC(0.71);
     stc2MC->SetY2NDC(0.9);
     stc2MC->SetTextColor(2);
     //  c2->Update();
     hChi2ITSsaMC->SetLineColor(4);
     hChi2ITSsaMC->Draw("sames");
     gPad->Update();
     TPaveStats *stc3MC=(TPaveStats*)hChi2ITSsaMC->GetListOfFunctions()->FindObject("stats");
     stc3MC->SetY1NDC(0.51);
     stc3MC->SetY2NDC(0.7);
     stc3MC->SetTextColor(4);
     legMC->Draw();
     
     cITSsa1MC->Update();
     cITSsa1MC->SaveAs("ITSsa1MC.pdf");
     pdfFileNames+=" ITSsa1MC.pdf";
     gStyle->SetPalette(1);
     hEtaPhiITSpureSAMC->SetStats(0);
     hEtaPhiITSpureSAMC->SetTitle("ITS pureSA - MC");
     hEtaPhiITSsaMC->SetStats(0);
     hEtaPhiITSsaMC->SetTitle("ITSsa tracks - MC");
     hEtaPhiTPCITSMC->SetStats(0);
     hEtaPhiTPCITSMC->SetTitle("TPC+ITS tracks - MC");
     ctitleMC=GetRunNumber()+"Eta-phi distribution for ITSsa and TPC+ITS tracks - MC";
     TCanvas* cITSsa2MC=new TCanvas("cITSsa2MC",ctitleMC,1200,800);
     clist[1]=cITSsa2MC;
     cITSsa2MC->Divide(3,1); //for ITSpuresa
     //cITSsa2->Divide(2,1); //for no ITSpureSA
     cITSsa2MC->cd(1);
     //TPad* p1=new TPad("p1","Tracking: tracks distribution in EtaPhi",0,0.5,1.,1.);
     //p1->Divide(3,1);
     //p1->cd (1);
     hEtaPhiITSpureSAMC->Draw("colz");
     hEtaPhiITSpureSAMC->GetXaxis()->SetTitle("Eta");
     hEtaPhiITSpureSAMC->GetYaxis()->SetTitle("Phi");
     cITSsa2MC->cd(2);
     //p1->cd(2);
     hEtaPhiITSsaMC->Draw("colz");
     hEtaPhiITSsaMC->GetXaxis()->SetTitle("Eta");
     hEtaPhiITSsaMC->GetYaxis()->SetTitle("Phi");
     cITSsa2MC->cd(3);
     //p1->cd(3);
     hEtaPhiTPCITSMC->Draw("colz");
     hEtaPhiTPCITSMC->GetXaxis()->SetTitle("Eta");
     hEtaPhiTPCITSMC->GetYaxis()->SetTitle("Phi");
     //c4->cd(4);
     cITSsa2MC->SaveAs("ITSsa2MC.pdf");
     pdfFileNames+=" ITSsa2MC.pdf";
 
      // fraction of pureSA tracks with hits in ITS layers - drawing
     TCanvas* cITSpureSAfrac=new TCanvas("cITSpureSAfrac",ctitleMC,1200,1200);
     clist[2]=cITSpureSAfrac;
     cITSpureSAfrac->Divide(1,2);
     cITSpureSAfrac->cd(1);
     histoFracpureSA->Draw();
     TLatex* ti=new TLatex(0.15,0.2,"Fraction of pureSA #pi tracks with point in ITS layer");
     ti->SetTextSize(0.04);
     ti->SetNDC();
     ti->SetTextColor(1);
     ti->Draw();
     TString testo2="Run "+GetRunNumber()+" -  DATA";
     TLatex* tg2i = new TLatex(0.15,0.85,testo2.Data());
     tg2i->SetTextSize(0.06);
     tg2i->SetNDC();
     tg2i->SetTextColor(2);
     tg2i->Draw();

     cITSpureSAfrac->cd(2);
     histoFracpureSAMC->Draw();
     ti->Draw();
     TString testo2mc="Run "+GetRunNumber()+" -  MC";
     TLatex* tg2imc = new TLatex(0.15,0.85,testo2mc.Data());
     tg2imc->SetTextSize(0.06);
     tg2imc->SetNDC();
     tg2imc->SetTextColor(2);
     tg2imc->Draw();

     cITSpureSAfrac->SaveAs("ITSpureSAFrac.pdf");
     pdfFileNames+=" ITSpureSAFrac.pdf";
     // fraction of pureSA tracks with hits in ITS layers - drawing
     
 }

//-----------------------------------------------------
// ///////////  Plot SDD ////////////
//_______________________________________________________________________
void SetDrawAtt(Int_t markerstyle,Int_t markercolor,Int_t markersize,Int_t linecolor,Int_t linewidth,TH1 *h1){ 


  h1->SetMarkerStyle(markerstyle);
  h1->SetMarkerColor(markercolor);
  h1->SetMarkerSize(markersize);
  h1->SetLineColor(linecolor);
  h1->SetLineWidth(linewidth);
}

//_______________________________________________________________________
Double_t LangausFun(Double_t *x, Double_t *par) {

  //Fit parameters:
  //par[0]=Width (scale) parameter of Landau density
  //par[1]=Most Probable (MP, location) parameter of Landau density
  //par[2]=Total area (integral -inf to inf, normalization constant)
  //par[3]=Width (sigma) of convoluted Gaussian function
  //
  //In the Landau distribution (represented by the CERNLIB approximation), 
  //the maximum is located at x=-0.22278298 with the location parameter=0.
  //This shift is corrected within this function, so that the actual
  //maximum is identical to the MP parameter.

  // Numeric constants
  Double_t invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
  Double_t mpshift  = -0.22278298;       // Landau maximum location

  // Control constants
  Double_t np = 100.0;      // number of convolution steps
  Double_t sc =   5.0;      // convolution extends to +-sc Gaussian sigmas

  // Variables
  Double_t xx;
  Double_t mpc;
  Double_t fland;
  Double_t sum = 0.0;
  Double_t xlow,xupp;
  Double_t step;
  Double_t i;


  // MP shift correction
  mpc = par[1] - mpshift * par[0]; 

  // Range of convolution integral
  xlow = x[0] - sc * par[3];
  xupp = x[0] + sc * par[3];

  step = (xupp-xlow) / np;

  // Convolution integral of Landau and Gaussian by sum
  for(i=1.0; i<=np/2; i++) {
    xx = xlow + (i-.5) * step;
    fland = TMath::Landau(xx,mpc,par[0]) / par[0];
    sum += fland * TMath::Gaus(x[0],xx,par[3]);

    xx = xupp - (i-.5) * step;
    fland = TMath::Landau(xx,mpc,par[0]) / par[0];
    sum += fland * TMath::Gaus(x[0],xx,par[3]);
  }

  return (par[2] * step * sum * invsq2pi / par[3]);

}
//_________________________________________________________________________
void PlotSDD(TFile* fildat, TFile* filMC, TCanvas**& clist, Int_t& cnum){
  TDirectoryFile* df=(TDirectoryFile*)fildat->Get("SDD_Performance");
    if(!df){
      printf("SDD_Performance data MISSING -> Exit\n");
      return;
    }
    TList* l=(TList*)df->Get("coutputRP");
    if(!df){
      printf("coutputRP TList data MISSING -> Exit\n");
      return;
    }
    cnum=1; // number of canvases 
    clist= new TCanvas* [cnum];//array of pointers to TCanvases
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(1111);
  TH1F* htimT=(TH1F*)l->FindObject("hDrTimTPAll");
  TH1F* htimTe=(TH1F*)l->FindObject("hDrTimTPExtra");
  TH1F* htimTne=(TH1F*)l->FindObject("hDrTimTPNoExtra");
  htimT->Rebin(4);
  htimTe->Rebin(4);
  htimTne->Rebin(4);
  htimT->SetLineWidth(2);
  htimTe->SetLineWidth(2);
  htimTne->SetLineWidth(2);  
  // TH1F* hev=(TH1F*)l->FindObject("hNEvents");
  // Int_t nTotEvents=hev->GetBinContent(2);
  // Int_t nTrigEvents=hev->GetBinContent(3);
  // Int_t nEvents=nTotEvents;
  // printf("---- Statistics ----\n");
  // printf("Number of Events = %d\n",nTotEvents);
  // if(nTrigEvents>0){ 
  //   printf("Number of Triggered Events = %d\n",nTrigEvents);
  //   nEvents=nTrigEvents;
  // }else{
  //   printf("No request on the trigger done when running the task\n");
  // }
  // if(hcllay){
  //   Double_t norm=hcllay->GetBinContent(1);
  //   if(norm>0.){
  //     hcllay->Scale(1./norm);
  //     hcllay->SetTitle("");
  //     hcllay->GetXaxis()->SetRange(2,7);
  //     hcllay->SetMinimum(0.);
  //     hcllay->SetMaximum(1.1);
  //     hcllay->SetMarkerStyle(23);
  //     TCanvas* ceffL=new TCanvas("ceffL","General: PointPerLayer",800,1000);
  //     clist[0]=ceffL;
  //     ceffL->Divide(1,2);
  //     ceffL->cd(1);
  //     ceffL->SetGridy();
  //     hcllay->Draw();
  //     hcllay->GetXaxis()->SetTitle("Layer");
  //     hcllay->GetYaxis()->SetTitle("Fraction of tracks with point in layer");
  //     ceffL->Update();
  //   }
  // }
  TH1F* hSigTim[8];
  TGraphErrors* gmpv=new TGraphErrors(0);
  TGraphErrors* gsigg=new TGraphErrors(0);
  TGraphErrors* gsigl=new TGraphErrors(0);
  gmpv->SetTitle("");
  gsigg->SetTitle("");
  gsigl->SetTitle("");
  Int_t iPoint=0;
  TF1 *lfun = new TF1("LangausFun",LangausFun,50.,300.,4);
  for(Int_t it=0; it<8; it++){
    hSigTim[it]=(TH1F*)l->FindObject(Form("hSigTimeInt%d",it));
    if(hSigTim[it]->GetEntries()>200){
      lfun->SetLineWidth(2);
      lfun->SetParameter(0,5.);
      lfun->SetParameter(1,80.);
      lfun->SetParameter(2,hSigTim[it]->GetEntries()/10.);
      lfun->SetParameter(3,10.);
      lfun->SetParLimits(3,0.,20);

      //      hSigTim[it]->Fit("LangausFun","QLR");
      hSigTim[it]->Fit("LangausFun","QON");
      hSigTim[it]->GetXaxis()->SetTitle(Form("dE/dx, time interval %d",it+1));
      hSigTim[it]->GetYaxis()->SetTitle("Events");
      Float_t mpv=lfun->GetParameter(1);
      Float_t empv=lfun->GetParError(1);
      Float_t sig=lfun->GetParameter(3);
      Float_t esig=lfun->GetParError(3);
      Float_t sigl=lfun->GetParameter(0);
      Float_t esigl=lfun->GetParError(0);
      gmpv->SetPoint(iPoint,(Float_t)it,mpv);
      gmpv->SetPointError(iPoint,0.,empv);
      gsigg->SetPoint(iPoint,(Float_t)it,sig);
      gsigg->SetPointError(iPoint,0.,esig);
      gsigl->SetPoint(iPoint,(Float_t)it,sigl);
      gsigl->SetPointError(iPoint,0.,esigl);
      ++iPoint;
      printf("Bin %d - MPV=%.3f  \t SigmaLandau=%.3f  \t SigmaGaus=%.3f\n",it,mpv,sigl,sig);
    }
  }
  TString ctitle=GetRunNumber()+"SDD: DriftTime - dE/dx - DATA";
  TCanvas* ctim=new TCanvas("ctim",ctitle,800,1000);
  clist[0]=ctim;
  ctim->Divide(1,2);
  ctim->cd(1);
  // htimT->Draw();
  // htimTe->SetLineColor(2);
  // htimTe->Draw("same");
  htimTne->SetLineColor(4);
  htimTne->Draw("");
  htimTne->GetXaxis()->SetTitle("Drift Time (ns)");
  htimTne->GetYaxis()->SetTitle("TrackPoints");
  htimTne->GetYaxis()->SetTitleOffset(1.2);
  htimTne->SetTitle("Drift Time from Track Points (ns) - DATA");
  // TLatex* ta=new TLatex(0.5,0.85,"All Clusters");
  // ta->SetNDC();
  // ta->SetTextColor(1);
  // ta->Draw();
  // TLatex* te=new TLatex(0.5,0.8,"Extra Clusters");
  // te->SetNDC();
  // te->SetTextColor(2);
  // te->Draw();
  //  TLatex* tn=new TLatex(0.3,0.3,"Non-Extra Clusters");
  TLatex* tn=new TLatex(0.3,0.3,"Clusters on SDD modules");
  tn->SetNDC();
  tn->SetTextColor(4);
  tn->Draw();
  TLine* tlin3=new TLine(450.,0.,450.,htimTne->GetMaximum());
    tlin3->SetLineColor(2);
    tlin3->SetLineWidth(2);
    tlin3->SetLineStyle(7);
    tlin3->Draw("same");
  TLine* tlin4=new TLine(620.,0.,620.,htimTne->GetMaximum());
    tlin4->SetLineColor(2);
    tlin4->SetLineWidth(2);
    tlin4->SetLineStyle(7);
    tlin4->Draw("same");
  TLatex* tlimit1=new TLatex(0.2,0.5,"Range for t0");
  tlimit1->SetNDC();
  tlimit1->SetTextColor(2);
  tlimit1->Draw();
  TLine* tlin5=new TLine(6200.,0.,6200.,htimTne->GetMaximum());
    tlin5->SetLineColor(2);
    tlin5->SetLineStyle(7);
    tlin5->SetLineWidth(2);
    tlin5->Draw("same");
  TLine* tlin6=new TLine(5150.,0.,5150.,htimTne->GetMaximum());
    tlin6->SetLineColor(2);
    tlin6->SetLineWidth(2);
    tlin6->SetLineStyle(7);
    tlin6->Draw("same");
  TLatex* tlimit2=new TLatex(0.6,0.5,"Range for falling edge");
  tlimit2->SetNDC();
  tlimit2->SetTextColor(2);
  tlimit2->Draw();

  //  ctim->Update();
  //  TCanvas* cpars=new TCanvas("cpars","Params",800,600);
  ctim->cd(2);
  gPad->SetLeftMargin(0.14);
  gPad->SetFrameLineWidth(2);
  gPad->SetTickx();
  gPad->SetTicky();
  gmpv->SetMarkerStyle(20);
  gmpv->SetMinimum(75);
  gmpv->SetMaximum(90);
  gmpv->GetXaxis()->SetLimits(-0.2,6.8);
  gmpv->Draw("AP");
  gmpv->GetXaxis()->SetTitle("Drift Time interval number");
  gmpv->GetYaxis()->SetTitle("Landau MPV (keV)");
  gmpv->GetXaxis()->SetTitleSize(0.05);
  gmpv->GetYaxis()->SetTitleSize(0.05);
  gmpv->GetYaxis()->SetTitleOffset(1.2);
  TLatex* tex=new TLatex(0.2,0.75,"dE/dx MPV vs Drift time interval - DATA");
  tex->SetNDC();
  tex->SetTextColor(1);
  tex->Draw();
 //  cpars->Update();
  ctim->Update();
  ctim->SaveAs("SDD.pdf");
  pdfFileNames+=" SDD.pdf";
    
// MC file
    
    TDirectoryFile* dfMC=(TDirectoryFile*)filMC->Get("SDD_Performance");
    if(!dfMC){
        printf("SDD_Performance MC MISSING -> Exit\n");
        return;
    }
    TList* lMC=(TList*)dfMC->Get("coutputRP");
    if(!dfMC){
        printf("coutputRP TList MC MISSING -> Exit\n");
        return;
    }
    cnum=1; // number of canvases
    clist= new TCanvas* [cnum];//array of pointers to TCanvases
    gROOT->SetStyle("Plain");
    gStyle->SetOptStat(1111);
    TH1F* htimTMC=(TH1F*)lMC->FindObject("hDrTimTPAll");
    TH1F* htimTeMC=(TH1F*)lMC->FindObject("hDrTimTPExtra");
    TH1F* htimTneMC=(TH1F*)lMC->FindObject("hDrTimTPNoExtra");
    htimTMC->Rebin(4);
    htimTeMC->Rebin(4);
    htimTneMC->Rebin(4);
    htimTMC->SetLineWidth(2);
    htimTeMC->SetLineWidth(2);
    htimTneMC->SetLineWidth(2);
    TH1F* hSigTimMC[8];
    TGraphErrors* gmpvMC=new TGraphErrors(0);
    TGraphErrors* gsiggMC=new TGraphErrors(0);
    TGraphErrors* gsiglMC=new TGraphErrors(0);
    gmpvMC->SetTitle("");
    gsiggMC->SetTitle("");
    gsiglMC->SetTitle("");
    Int_t iPointMC=0;
    TF1 *lfunMC = new TF1("LangausFunMC",LangausFun,50.,300.,4);
    for(Int_t it=0; it<8; it++){
        hSigTimMC[it]=(TH1F*)lMC->FindObject(Form("hSigTimeInt%d",it));
        if(hSigTimMC[it]->GetEntries()>200){
            lfunMC->SetLineWidth(2);
            lfunMC->SetParameter(0,5.);
            lfunMC->SetParameter(1,80.);
            lfunMC->SetParameter(2,hSigTim[it]->GetEntries()/10.);
            lfunMC->SetParameter(3,10.);
            lfunMC->SetParLimits(3,0.,20);
            
            //      hSigTimMC[it]->Fit("LangausFun","QLR");
            hSigTimMC[it]->Fit("LangausFunMC","QON");
            hSigTimMC[it]->GetXaxis()->SetTitle(Form("dE/dx, time interval %d",it+1));
            hSigTimMC[it]->GetYaxis()->SetTitle("Events");
            Float_t mpvMC=lfunMC->GetParameter(1);
            Float_t empvMC=lfunMC->GetParError(1);
            Float_t sigMC=lfunMC->GetParameter(3);
            Float_t esigMC=lfunMC->GetParError(3);
            Float_t siglMC=lfunMC->GetParameter(0);
            Float_t esiglMC=lfunMC->GetParError(0);
            gmpvMC->SetPoint(iPointMC,(Float_t)it,mpvMC);
            gmpvMC->SetPointError(iPointMC,0.,empvMC);
            gsiggMC->SetPoint(iPointMC,(Float_t)it,sigMC);
            gsiggMC->SetPointError(iPointMC,0.,esigMC);
            gsiglMC->SetPoint(iPointMC,(Float_t)it,siglMC);
            gsiglMC->SetPointError(iPointMC,0.,esiglMC);
            ++iPointMC;
            printf("Bin %d - MPV=%.3f  \t SigmaLandau=%.3f  \t SigmaGaus=%.3f\n",it,mpvMC,siglMC,sigMC);
        }
    }
    TString ctitleMC=GetRunNumber()+"SDD: DriftTime - dE/dx MC - MC";
    TCanvas* ctimMC=new TCanvas("ctimMC",ctitleMC,800,1000);
    clist[0]=ctimMC;
    ctimMC->Divide(1,2);
    ctimMC->cd(1);
    // htimTMC->Draw();
    // htimTeMC->SetLineColor(2);
    // htimTeMC->Draw("same");
    htimTneMC->SetLineColor(4);
    htimTneMC->Draw("");
    htimTneMC->GetXaxis()->SetTitle("Drift Time (ns)");
    htimTneMC->GetYaxis()->SetTitle("TrackPoints");
    htimTneMC->GetYaxis()->SetTitleOffset(1.2);
    htimTneMC->SetTitle("Drift Time from Track Points (ns) - MC");

    // TLatex* ta=new TLatex(0.5,0.85,"All Clusters");
    // ta->SetNDC();
    // ta->SetTextColor(1);
    // ta->Draw();
    // TLatex* te=new TLatex(0.5,0.8,"Extra Clusters");
    // te->SetNDC();
    // te->SetTextColor(2);
    // te->Draw();
    //  TLatex* tn=new TLatex(0.3,0.3,"Non-Extra Clusters");
    TLatex* tnMC=new TLatex(0.3,0.3,"Clusters on SDD modules");
    tnMC->SetNDC();
    tnMC->SetTextColor(4);
    tnMC->Draw();
    TLine* tlin3MC=new TLine(450.,0.,450.,htimTneMC->GetMaximum());
    tlin3MC->SetLineColor(2);
    tlin3MC->SetLineWidth(2);
    tlin3MC->SetLineStyle(7);
    tlin3MC->Draw("same");
    TLine* tlin4MC=new TLine(620.,0.,620.,htimTneMC->GetMaximum());
    tlin4MC->SetLineColor(2);
    tlin4MC->SetLineWidth(2);
    tlin4MC->SetLineStyle(7);
    tlin4MC->Draw("same");
    TLatex* tlimit1MC=new TLatex(0.2,0.5,"Range for t0");
    tlimit1MC->SetNDC();
    tlimit1MC->SetTextColor(2);
    tlimit1MC->Draw();
    TLine* tlin5MC=new TLine(6200.,0.,6200.,htimTneMC->GetMaximum());
    tlin5MC->SetLineColor(2);
    tlin5MC->SetLineStyle(7);
    tlin5MC->SetLineWidth(2);
    tlin5MC->Draw("same");
    TLine* tlin6MC=new TLine(5150.,0.,5150.,htimTneMC->GetMaximum());
    tlin6MC->SetLineColor(2);
    tlin6MC->SetLineWidth(2);
    tlin6MC->SetLineStyle(7);
    tlin6MC->Draw("same");
    TLatex* tlimit2MC=new TLatex(0.6,0.5,"Range for falling edge");
    tlimit2MC->SetNDC();
    tlimit2MC->SetTextColor(2);
    tlimit2MC->Draw();
    
    //  ctim->Update();
    //  TCanvas* cpars=new TCanvas("cpars","Params",800,600);
    ctimMC->cd(2);
    gPad->SetLeftMargin(0.14);
    gPad->SetFrameLineWidth(2);
    gPad->SetTickx();
    gPad->SetTicky();
    gmpvMC->SetMarkerStyle(20);
    gmpvMC->SetMinimum(75);
    gmpvMC->SetMaximum(90);
    gmpvMC->GetXaxis()->SetLimits(-0.2,6.8);
    gmpvMC->Draw("AP");
    gmpvMC->GetXaxis()->SetTitle("Drift Time interval number");
    gmpvMC->GetYaxis()->SetTitle("Landau MPV (keV)");
    gmpvMC->GetXaxis()->SetTitleSize(0.05);
    gmpvMC->GetYaxis()->SetTitleSize(0.05);
    gmpvMC->GetYaxis()->SetTitleOffset(1.2);
    TLatex* texMC=new TLatex(0.2,0.75,"dE/dx MPV vs Drift time interval - MC");
    texMC->SetNDC();
    texMC->SetTextColor(1);
    texMC->Draw();
    //  cpars->Update();
    ctimMC->Update();
    ctimMC->SaveAs("SDDMC.pdf");
    pdfFileNames+=" SDDMC.pdf";
   
} 

//_______________________________________________________________________
//////////////// SSD ///////////////////////
//_______________________________________________________________________
void GetGainModuleLevelSSD(TFile* fildat, TFile* filMC, TCanvas**& clist, Int_t& cnum)
{
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(1111);
  gStyle->SetPalette(1,0);
  cnum=1;
  clist=new TCanvas*[1]; 

  TDirectoryFile* df=(TDirectoryFile*)fildat->Get("PWGPPdEdxSSDQA");
    if(!df){
        printf("SSD_Performance MISSING -> Exit\n");
        cnum=-1;
        return;
    }

  TList* listin=(TList*)df->Get("SSDdEdxQA");
  if(!listin) {
        printf("SSDdEdxQA TList MISSING -> Exit\n");
      cnum=-1;
        return;
    }
  
  TH2F* fHistQ=0x0;
  fHistQ=(TH2F*)listin ->FindObject("QACharge");
  fHistQ->SetStats(111);
  fHistQ->SetTitle("SSD Charge vs module number - DATA");
  if(!fHistQ) return;
    
  TH2F* fHistCR=(TH2F*)listin ->FindObject("QAChargeRatio");
  fHistCR->SetStats(0);
  fHistCR->SetTitle("SSD Charge Ratio vs module number - DATA");

  if(!fHistCR) return;

  TH1F* fHistMPVs=new TH1F("SSD HistMPVS","HistMPVs - DATA;MPV;N",75,70,95);
	
  TH1F* fHistCRmean=new TH1F("SSD HistCRmean","HistCRmean - DATA;CRmean;N",200,-1,1);
	
  TH1F *fMPVGraph = new TH1F("SSD MPV","MPVgraph - DATA;Module number;MPV",1698,-0.5,1697.5);
  fMPVGraph->SetMarkerColor(kRed);
  fMPVGraph->SetMarkerSize(0.5);
  fMPVGraph->SetMarkerStyle(22);
  fMPVGraph->SetStats(111111);
  
  TH1F *fCRmeanGraph = new TH1F("SSD CRmeangraph","CRmeangraph - DATA;Module number;MPV",1698,-0.5,1697.5);
  fCRmeanGraph->SetMarkerColor(kBlue);
  fCRmeanGraph->SetMarkerSize(0.5);
  fCRmeanGraph->SetMarkerStyle(23);
  fCRmeanGraph->SetStats(111111);

  Float_t mpv[1698];
  Int_t ntofit=200;
 
  //  ofstream outfiletxtbad;
  //outfiletxtbad.open("QALHC11eCPass1_bis/badModules.txt");	
    for (int i =0;i<1698;i++)
    {
      //      cout<<i<<endl;
      TString tmpQ("Q");
      tmpQ+=i;
      TString tmpCR("CR");
      tmpCR+=i;
      TH1D* fHist1DCR= fHistCR->ProjectionY(tmpCR,i+1,i+1);
      Double_t mean=fHist1DCR->GetMean();
      if(!(TMath::Abs(mean)<1.0)||fHist1DCR->GetEntries()<10)
	  continue;
      fHistCRmean->Fill(mean);
      fCRmeanGraph->SetBinContent(i+1,mean);
      fCRmeanGraph->SetBinError(i+1,fHist1DCR->GetRMS());
      fCRmeanGraph->GetYaxis()->SetTitle("CR");
      TH1D* fHist1DQ=fHistQ->ProjectionY(tmpQ,i+1,i+1);
      //check bad modules
      if(fHist1DQ->GetEntries()<ntofit)
	{
	  //outfiletxtbad<<"Low statistic \t module= "<<i<<" netries="<<fHist1DQ->GetEntries()<<endl;
	  continue;
	}
      else
	{
	  tmpQ+="fit";
	  Float_t range=fHist1DQ->GetBinCenter(fHist1DQ->GetMaximumBin());
	  TF1 *f1 = new TF1(tmpQ,LangausFun,range*0.45,range*3.0,4);
	  f1->SetParameters(7.0,range,1.0,5.5);
	  Float_t normalization=fHist1DQ->GetEntries()*fHist1DQ->GetXaxis()->GetBinWidth(2)/f1->Integral(range*0.45,range*3.0);
	  f1->SetParameters(7.0,range,normalization,5.5);
	  //f1->SetParameters(7.0,range,fHist1DQ->GetMaximum(),5.5);
	  f1->SetParNames("sigma Landau","MPV","N","sigma Gaus");
	  f1->SetParLimits(0,2.0,100.0);
	  f1->SetParLimits(3,0.0,100.0);
	  if(fHist1DQ->Fit(tmpQ,"BRQON")==0)
	    {
	      mpv[i]=f1->GetParameter(1);
	      fHistMPVs->Fill(mpv[i]);	
	      fMPVGraph->SetBinContent(i+1,f1->GetParameter(1));
	      fMPVGraph->SetBinError(i+1,f1->GetParError(1));
	      if(mpv[i]<75.0)
		{
		  //outfiletxtbad<<"MPV lower than 75 \t module="<<i<<endl;
		}	
	      if(mpv[i]>100.0)
		{
		  // outfiletxtbad<<"MPV higher than 100 \t module="<<i<<endl;		  
		}
	      if(f1->GetParError(1)>1.0)
		{
		  //outfiletxtbad<<"MPV high error on MPV  \t module="<<i<<endl;				
		}
	    }
	  else
	    {
	      mpv[i]=1;
	      //outfiletxtbad<<"BAD FIT \t module="<<i<<endl;
	      continue;
	    }	
	}	
    }	
  
    TString ctitle=GetRunNumber()+"SSD Calibration 1 - DATA";
  TCanvas *c1SSD = new TCanvas("c1SSD",ctitle,1000,1000);
  clist[0]=c1SSD;
  c1SSD->Divide(2,3);
  c1SSD->cd(1);
  fHistQ->DrawCopy("colz");
  c1SSD->cd(2);
  fHistCR->DrawCopy("colz");

  //  TCanvas *c2SSD = new TCanvas("c2SSD","SSD Calibration 2",1000,1000);
  //clist[1]=c2SSD;
  //c2SSD->Divide(2,2);
  c1SSD->cd(3);
  fMPVGraph->DrawCopy();  
     TLine* tlin0=new TLine(0.,80.,1698.,80.);
    tlin0->SetLineColor(2);
    tlin0->SetLineWidth(2);
    tlin0->Draw("same");
     TLine* tlin01=new TLine(0.,90.,1698.,90.);
    tlin01->SetLineColor(2);
    tlin01->SetLineWidth(2);
    tlin01->Draw("same");
  c1SSD->cd(4); 
  fHistMPVs->DrawCopy();
  c1SSD->cd(5);
  fCRmeanGraph->DrawCopy();
    TLine* tlin1=new TLine(0.,0.2,1698.,0.2);
    tlin1->SetLineColor(2);
    tlin1->SetLineWidth(2);
    tlin1->Draw("same");
    TLine* tlin2=new TLine(0.,-0.2,1698.,-0.2);
    tlin2->SetLineColor(2);
    tlin2->SetLineWidth(2);
    tlin2->Draw("same");
  TLatex* ta1=new TLatex(0.2,0.8,"SSD Calibration");
  ta1->SetNDC();
  ta1->SetTextSize(0.05);
  ta1->SetTextColor(2);
  ta1->Draw("same");
  c1SSD->cd(6);
  fHistCRmean->DrawCopy();
  c1SSD->Update();
  c1SSD->SaveAs("SSD.pdf");
  pdfFileNames+=" SSD.pdf";
    
// MC file
    clist=new TCanvas*[1];
    
    TDirectoryFile* dfMC=(TDirectoryFile*)filMC->Get("PWGPPdEdxSSDQA");
    TList* listinMC=(TList*)dfMC->Get("SSDdEdxQA");
    if(!listinMC) return;
    TH2F* fHistQMC=0x0;
    fHistQMC=(TH2F*)listinMC ->FindObject("QACharge");
    fHistQMC->SetStats(111);
    fHistQMC->SetTitle("SSD Charge vs module number - MC");
    if(!fHistQMC) return;
    TH2F* fHistCRMC=(TH2F*)listinMC ->FindObject("QAChargeRatio");
    fHistCRMC->SetStats(0);
    fHistCRMC->SetTitle("SSD Charge Ratio vs module number - MC");
    
    if(!fHistCRMC) return;
    
    TH1F* fHistMPVsMC=new TH1F("SSD HistMPVS MC","HistMPVs - MC;MPV;N",75,70,95);
    
    TH1F* fHistCRmeanMC=new TH1F("SSD HistCRmean MC","HistCRmean - MC;CRmean;N",200,-1,1);
    
    TH1F *fMPVGraphMC = new TH1F("SSD MPV MC","MPVgraph - MC;Module number;MPV",1698,-0.5,1697.5);
    fMPVGraphMC->SetMarkerColor(kRed);
    fMPVGraphMC->SetMarkerSize(0.5);
    fMPVGraphMC->SetMarkerStyle(22);
    fMPVGraphMC->SetStats(111111);
    
    TH1F *fCRmeanGraphMC = new TH1F("SSD CRmeangraph MC","CRmeangraph - MC;Module number;MPV",1698,-0.5,1697.5);
    fCRmeanGraphMC->SetMarkerColor(kBlue);
    fCRmeanGraphMC->SetMarkerSize(0.5);
    fCRmeanGraphMC->SetMarkerStyle(23);
    fCRmeanGraphMC->SetStats(111111);
    
    Float_t mpvMC[1698];
    Int_t ntofitMC=200;
    
    //  ofstream outfiletxtbad;
    //outfiletxtbad.open("QALHC11eCPass1_bis/badModules.txt");
    for (int i =0;i<1698;i++)
    {
        //      cout<<i<<endl;
        TString tmpQMC("Q");
        tmpQMC+=i;
        TString tmpCRMC("CR");
        tmpCRMC+=i;
        TH1D* fHist1DCRMC= fHistCRMC->ProjectionY(tmpCRMC,i+1,i+1);
        Double_t meanMC=fHist1DCRMC->GetMean();
        if(!(TMath::Abs(meanMC)<1.0)||fHist1DCRMC->GetEntries()<10)
            continue;
        fHistCRmeanMC->Fill(meanMC);
        fCRmeanGraphMC->SetBinContent(i+1,meanMC);
        fCRmeanGraphMC->SetBinError(i+1,fHist1DCRMC->GetRMS());
        fCRmeanGraphMC->GetYaxis()->SetTitle("CR");
        TH1D* fHist1DQMC=fHistQMC->ProjectionY(tmpQMC,i+1,i+1);
        //check bad modules
        if(fHist1DQMC->GetEntries()<ntofit)
        {
            //outfiletxtbad<<"Low statistic \t module= "<<i<<" netries="<<fHist1DQ->GetEntries()<<endl;
            continue;
        }
        else
        {
            tmpQMC+="fit";
            Float_t rangeMC=fHist1DQMC->GetBinCenter(fHist1DQMC->GetMaximumBin());
            TF1 *f1MC = new TF1(tmpQMC,LangausFun,rangeMC*0.45,rangeMC*3.0,4);
            f1MC->SetParameters(7.0,rangeMC,1.0,5.5);
            Float_t normalizationMC=fHist1DQMC->GetEntries()*fHist1DQMC->GetXaxis()->GetBinWidth(2)/f1MC->Integral(rangeMC*0.45,rangeMC*3.0);
            f1MC->SetParameters(7.0,rangeMC,normalizationMC,5.5);
            //f1->SetParameters(7.0,range,fHist1DQ->GetMaximum(),5.5);
            f1MC->SetParNames("sigma Landau","MPV","N","sigma Gaus");
            f1MC->SetParLimits(0,2.0,100.0);
            f1MC->SetParLimits(3,0.0,100.0);
            if(fHist1DQMC->Fit(tmpQMC,"BRQON")==0)
            {
                mpvMC[i]=f1MC->GetParameter(1);
                fHistMPVsMC->Fill(mpvMC[i]);
                fMPVGraphMC->SetBinContent(i+1,f1MC->GetParameter(1));
                fMPVGraphMC->SetBinError(i+1,f1MC->GetParError(1));
                if(mpvMC[i]<75.0)
                {
                    //outfiletxtbad<<"MPV lower than 75 \t module="<<i<<endl;
                }
                if(mpvMC[i]>100.0)
                {
                    // outfiletxtbad<<"MPV higher than 100 \t module="<<i<<endl;
                }
                if(f1MC->GetParError(1)>1.0)
                {
                    //outfiletxtbad<<"MPV high error on MPV  \t module="<<i<<endl;
                }
            }
            else
            {
                mpvMC[i]=1;
                //outfiletxtbad<<"BAD FIT \t module="<<i<<endl;
                continue;
            }
        }
    }
    
    TString ctitleMC=GetRunNumber()+"SSD Calibration 1 - MC";
    TCanvas *c1SSDMC = new TCanvas("c1SSDMC",ctitleMC,1000,1000);
    clist[0]=c1SSDMC;
    c1SSDMC->Divide(2,3);
    c1SSDMC->cd(1);
    fHistQMC->DrawCopy("colz");
    c1SSDMC->cd(2);
    fHistCRMC->DrawCopy("colz");
    
    //  TCanvas *c2SSD = new TCanvas("c2SSD","SSD Calibration 2",1000,1000);
    //clist[1]=c2SSD;
    //c2SSD->Divide(2,2);
    c1SSDMC->cd(3);
    fMPVGraphMC->DrawCopy();
    TLine* tlin0MC=new TLine(0.,80.,1698.,80.);
    tlin0MC->SetLineColor(2);
    tlin0MC->SetLineWidth(2);
    tlin0MC->Draw("same");
    TLine* tlin01MC=new TLine(0.,90.,1698.,90.);
    tlin01MC->SetLineColor(2);
    tlin01MC->SetLineWidth(2);
    tlin01MC->Draw("same");
    c1SSDMC->cd(4);
    fHistMPVsMC->DrawCopy();
    c1SSDMC->cd(5);
    fCRmeanGraphMC->DrawCopy();
    TLine* tlin1MC=new TLine(0.,0.2,1698.,0.2);
    tlin1MC->SetLineColor(2);
    tlin1MC->SetLineWidth(2);
    tlin1MC->Draw("same");
    TLine* tlin2MC=new TLine(0.,-0.2,1698.,-0.2);
    tlin2MC->SetLineColor(2);
    tlin2MC->SetLineWidth(2);
    tlin2MC->Draw("same");
    TLatex* ta1MC=new TLatex(0.2,0.8,"SSD Calibration");
    ta1MC->SetNDC();
    ta1MC->SetTextSize(0.05);
    ta1MC->SetTextColor(2);
    ta1MC->Draw("same");
    c1SSDMC->cd(6);
    fHistCRmeanMC->DrawCopy();
    c1SSDMC->Update();
    c1SSDMC->SaveAs("SSDMC.pdf");
    pdfFileNames+=" SSDMC.pdf";

}

//_______________________________________________________________________
void VertexQAMacro(TFile *fildat, TFile *filMC, TCanvas **&clist, Int_t &cnum){

	TDirectoryFile *dir = (TDirectoryFile*)fildat->Get("Vertex_Performance");
	if(!dir){
		Printf("Vertex directory not found... check!");
        cnum=-1;
        return;
	}
	
	TList *lt = (TList*)dir->Get("cOutputVtxESD");
    if(!lt){
        Printf("cOutputVtxESD TList not found... check!");
        cnum=-1;
        return;
    }
    
	
	cnum = 1;
	clist = new TCanvas*[1];

	
	TH1F *xVtxSPD = (TH1F*)lt->FindObject("fhSPDVertexX");
	TH1F *yVtxSPD = (TH1F*)lt->FindObject("fhSPDVertexY");
    TH1F *zVtxSPD = (TH1F*)lt->FindObject("fhSPDVertexZ");
	
	TH1F *zVtxSPD_Zonly = (TH1F*)lt->FindObject("fhSPDVertexZonly");
	
	if(!zVtxSPD_Zonly){
		Printf("using SPD 3D histo, Zonly not available");
	    zVtxSPD_Zonly = (TH1F*)lt->FindObject("fhSPDVertexZ");	
	}
	
	TH1F *xVtxTRK = (TH1F*)lt->FindObject("fhTRKVertexX");
	TH1F *yVtxTRK = (TH1F*)lt->FindObject("fhTRKVertexY");
	TH1F *zVtxTRK = (TH1F*)lt->FindObject("fhTRKVertexZ");
	
	TH2F *hntrksSPDvsSPDcls = (TH2F*)lt->FindObject("fhntrksSPDvsSPDcls");
    TH2F *hntrksZvsSPDcls = (TH2F*)lt->FindObject("fhntrksZvsSPDcls");
	
	Bool_t histoCorelation = kTRUE;
	
	if(!hntrksZvsSPDcls){
		Printf("skipping the second part, no 2D histos available");
		histoCorelation=kFALSE;	
	}

	TString ctitle=GetRunNumber()+"TRKandSPD3DxVtx - DATA";
	TCanvas *TRK_SPD3D_Vtx = new TCanvas("TRKandSPD3DVtx",ctitle,1000,1000);
	TRK_SPD3D_Vtx->Divide(3,2);
	clist[0]=TRK_SPD3D_Vtx;
	gStyle->SetOptFit(111);

	TRK_SPD3D_Vtx->cd(1);
	xVtxSPD->SetMarkerStyle(20);
	xVtxSPD->SetLineWidth(3);
	xVtxSPD->SetMarkerColor(kBlue+2);
	TF1 *fx = new TF1("gaus", "gaus", -1, 1);
  	xVtxTRK->SetMarkerStyle(20);
	xVtxTRK->SetLineWidth(4);
	xVtxTRK->SetLineColor(2);
	xVtxTRK->Draw("PE");
	xVtxTRK->Fit("gaus", "M");
	xVtxSPD->Draw("PE SAME");
	xVtxTRK->GetXaxis()->SetRangeUser(-0.05, 0.15);
	xVtxSPD->GetXaxis()->SetRangeUser(-0.05, 0.15);
	
	TLatex* tVTX1=new TLatex(0.15,0.85,"VertexSPD - DATA");
    tVTX1->SetNDC();
    tVTX1->SetTextColor(kBlue+2);
    tVTX1->Draw();
	TLatex* tVTX2=new TLatex(0.15,0.8,"VertexTRK - DATA");
    tVTX2->SetNDC();
    tVTX2->SetTextColor(2);
    tVTX2->Draw();
	
	TRK_SPD3D_Vtx->cd(2);
	yVtxSPD->SetMarkerStyle(20);
	yVtxSPD->SetLineWidth(3);
	yVtxSPD->SetMarkerColor(kBlue+2);
	TF1 *fy = new TF1("gaus", "gaus", -1, 1);
	yVtxTRK->SetMarkerStyle(20);
	yVtxTRK->SetLineWidth(3);
	yVtxTRK->SetLineColor(2);
	yVtxTRK->Draw("PE");
	yVtxTRK->Fit("gaus", "M");
	yVtxSPD->Draw("PE SAME");
	yVtxTRK->GetXaxis()->SetRangeUser(-0.2, 0.6);
	yVtxSPD->GetXaxis()->SetRangeUser(-0.2, 0.6);

	TLatex* tVTX3=new TLatex(0.15,0.85,"VertexSPD - DATA");
  tVTX3->SetNDC();
  tVTX3->SetTextColor(kBlue+2);
  tVTX3->Draw();
	TLatex* tVTX4=new TLatex(0.15,0.8,"VertexTRK - DATA");
  tVTX4->SetNDC();
  tVTX4->SetTextColor(2);
  tVTX4->Draw();
 	
	
	
	TRK_SPD3D_Vtx->cd(3);
	
	TF1 *fz = new TF1("gaus", "gaus", -20, 20);
	zVtxTRK->SetMarkerStyle(20);
	zVtxTRK->SetLineWidth(3);
	zVtxTRK->SetMarkerColor(2);
	zVtxTRK->SetLineColor(2);
	zVtxTRK->Draw("PE");
	zVtxTRK->Fit("gaus", "M");
	zVtxSPD->SetMarkerStyle(20);
	zVtxSPD->SetLineWidth(1);
	zVtxSPD->SetLineColor(kBlue+2);
	zVtxSPD->SetMarkerColor(kBlue+2);
	zVtxSPD->SetMarkerSize(0.8);
	zVtxSPD->Draw("PE SAME");
	TLatex* tVTX5=new TLatex(0.15,0.85,"VertexSPD - DATA");
  tVTX5->SetNDC();
  tVTX5->SetTextColor(kBlue+2);
  tVTX5->Draw();
	TLatex* tVTX6=new TLatex(0.15,0.8,"VertexTRK - DATA");
  tVTX6->SetNDC();
  tVTX6->SetTextColor(2);
  tVTX6->Draw();

	
	
	//	TCanvas *corrContrSPDClusters = new TCanvas("corrContrSPDClusters", "corrContrSPDClusters");
	// corrContrSPDClusters->Divide(3,1);
	//clist[1]=corrContrSPDClusters;
		TRK_SPD3D_Vtx->cd(4);
		//	corrContrSPDClusters->cd(1);
	zVtxSPD_Zonly->SetLineWidth(3);
	zVtxSPD_Zonly->SetLineColor(kBlue+2);
	zVtxSPD_Zonly->Draw();
	TLatex* tVTX7=new TLatex(0.15,0.8,"Vertex Z only - DATA");
  tVTX7->SetNDC();
  tVTX7->SetTextColor(2);
  tVTX7->Draw();
	
	if(histoCorelation){
			TRK_SPD3D_Vtx->cd(5);
			//corrContrSPDClusters->cd(2);
		hntrksSPDvsSPDcls->SetMarkerStyle(20);
	        hntrksSPDvsSPDcls->GetXaxis()->SetRangeUser(0.,400);
		hntrksSPDvsSPDcls->GetYaxis()->SetRangeUser(0.,4500.);
        hntrksSPDvsSPDcls->SetTitle("ncontributors SPD3D vs number of cluster SPD - DATA");
		hntrksSPDvsSPDcls->Draw();
	
		TRK_SPD3D_Vtx->cd(6);
		//		corrContrSPDClusters->cd(3);
		hntrksZvsSPDcls->SetMarkerStyle(20);
		hntrksZvsSPDcls->GetXaxis()->SetRangeUser(0.,150.);
	        hntrksZvsSPDcls->GetYaxis()->SetRangeUser(0.,1000.);
        hntrksZvsSPDcls->SetTitle("ncontributors SPDZ vs number of cluster SPD - DATA");
		hntrksZvsSPDcls->Draw();
   }	
	TRK_SPD3D_Vtx->SaveAs("vertex.pdf");
	pdfFileNames+=" vertex.pdf";
	delete fx;
	delete fy;
	delete fz;

// MC file
    TDirectoryFile *dirMC = (TDirectoryFile*)filMC->Get("Vertex_Performance");
    if(!dirMC){
        Printf("Vertex directory MC not found... check!");
    }
    
    TList *ltMC = (TList*)dirMC->Get("cOutputVtxESD");
    
    cnum = 1;
    clist = new TCanvas*[1];
    
    
    TH1F *xVtxSPDMC = (TH1F*)ltMC->FindObject("fhSPDVertexX");
    TH1F *yVtxSPDMC = (TH1F*)ltMC->FindObject("fhSPDVertexY");
    TH1F *zVtxSPDMC = (TH1F*)ltMC->FindObject("fhSPDVertexZ");
    
    TH1F *zVtxSPD_ZonlyMC = (TH1F*)ltMC->FindObject("fhSPDVertexZonly");
    
    if(!zVtxSPD_ZonlyMC){
        Printf("using SPD 3D histo, Zonly not available");
        zVtxSPD_ZonlyMC = (TH1F*)ltMC->FindObject("fhSPDVertexZ");
    }
    
    TH1F *xVtxTRKMC = (TH1F*)ltMC->FindObject("fhTRKVertexX");
    TH1F *yVtxTRKMC = (TH1F*)ltMC->FindObject("fhTRKVertexY");
    TH1F *zVtxTRKMC = (TH1F*)ltMC->FindObject("fhTRKVertexZ");
    
    TH2F *hntrksSPDvsSPDclsMC = (TH2F*)ltMC->FindObject("fhntrksSPDvsSPDcls");
    TH2F *hntrksZvsSPDclsMC = (TH2F*)ltMC->FindObject("fhntrksZvsSPDcls");
    
    Bool_t histoCorelationMC = kTRUE;
    
    if(!hntrksZvsSPDclsMC){
        Printf("skipping the second part, no 2D histos available");
        histoCorelationMC=kFALSE;
    }
    
    TString ctitleMC=GetRunNumber()+"TRKandSPD3DxVtx - MC";
    TCanvas *TRK_SPD3D_VtxMC = new TCanvas("TRKandSPD3DVtxMC",ctitleMC,1000,1000);
    TRK_SPD3D_VtxMC->Divide(3,2);
    clist[0]=TRK_SPD3D_VtxMC;
    gStyle->SetOptFit(111);
    
    TRK_SPD3D_VtxMC->cd(1);
    xVtxSPDMC->SetMarkerStyle(20);
    xVtxSPDMC->SetLineWidth(3);
    xVtxSPDMC->SetMarkerColor(kBlue+2);
    TF1 *fxMC = new TF1("gaus", "gaus", -1, 1);
    xVtxTRKMC->SetMarkerStyle(20);
    xVtxTRKMC->SetLineWidth(4);
    xVtxTRKMC->SetLineColor(2);
    xVtxTRKMC->Draw("PE");
    xVtxTRKMC->Fit("gaus", "M");
    xVtxSPDMC->Draw("PE SAME");
    xVtxTRKMC->GetXaxis()->SetRangeUser(-0.05, 0.15);
    xVtxSPDMC->GetXaxis()->SetRangeUser(-0.05, 0.15);
    
    TLatex* tVTX1MC=new TLatex(0.15,0.85,"VertexSPD - MC");
    tVTX1MC->SetNDC();
    tVTX1MC->SetTextColor(kBlue+2);
    tVTX1MC->Draw();
    TLatex* tVTX2MC=new TLatex(0.15,0.8,"VertexTRK - MC");
    tVTX2MC->SetNDC();
    tVTX2MC->SetTextColor(2);
    tVTX2MC->Draw();
    
    TRK_SPD3D_VtxMC->cd(2);
    yVtxSPDMC->SetMarkerStyle(20);
    yVtxSPDMC->SetLineWidth(3);
    yVtxSPDMC->SetMarkerColor(kBlue+2);
    TF1 *fyMC = new TF1("gaus", "gaus", -1, 1);
    yVtxTRKMC->SetMarkerStyle(20);
    yVtxTRKMC->SetLineWidth(3);
    yVtxTRKMC->SetLineColor(2);
    yVtxTRKMC->Draw("PE");
    yVtxTRKMC->Fit("gaus", "M");
    yVtxSPDMC->Draw("PE SAME");
    yVtxTRKMC->GetXaxis()->SetRangeUser(-0.2, 0.6);
    yVtxSPDMC->GetXaxis()->SetRangeUser(-0.2, 0.6);
    
    TLatex* tVTX3MC=new TLatex(0.15,0.85,"VertexSPD - MC");
    tVTX3MC->SetNDC();
    tVTX3MC->SetTextColor(kBlue+2);
    tVTX3MC->Draw();
    TLatex* tVTX4MC=new TLatex(0.15,0.8,"VertexTRK - MC");
    tVTX4MC->SetNDC();
    tVTX4MC->SetTextColor(2);
    tVTX4MC->Draw();
    
    
    TRK_SPD3D_VtxMC->cd(3);
    
    TF1 *fzMC = new TF1("gaus", "gaus", -20, 20);
    zVtxTRKMC->SetMarkerStyle(20);
    zVtxTRKMC->SetLineWidth(3);
    zVtxTRKMC->SetMarkerColor(2);
    zVtxTRKMC->SetLineColor(2);
    zVtxTRKMC->Draw("PE");
    zVtxTRKMC->Fit("gaus", "M");
    zVtxSPDMC->SetMarkerStyle(20);
    zVtxSPDMC->SetLineWidth(1);
    zVtxSPDMC->SetLineColor(kBlue+2);
    zVtxSPDMC->SetMarkerColor(kBlue+2);
    zVtxSPDMC->SetMarkerSize(0.8);
    zVtxSPDMC->Draw("PE SAME");
    TLatex* tVTX5MC=new TLatex(0.15,0.85,"VertexSPD - MC");
    tVTX5MC->SetNDC();
    tVTX5MC->SetTextColor(kBlue+2);
    tVTX5MC->Draw();
    TLatex* tVTX6MC=new TLatex(0.15,0.8,"VertexTRK - MC");
    tVTX6MC->SetNDC();
    tVTX6MC->SetTextColor(2);
    tVTX6MC->Draw();
    
    
    //	TCanvas *corrContrSPDClusters = new TCanvas("corrContrSPDClusters", "corrContrSPDClusters");
    // corrContrSPDClusters->Divide(3,1);
    //clist[1]=corrContrSPDClusters;
    TRK_SPD3D_VtxMC->cd(4);
    //	corrContrSPDClusters->cd(1);
    zVtxSPD_ZonlyMC->SetLineWidth(3);
    zVtxSPD_ZonlyMC->SetLineColor(kBlue+2);
    zVtxSPD_ZonlyMC->Draw();
    TLatex* tVTX7MC=new TLatex(0.15,0.8,"Vertex Z only - MC");
    tVTX7MC->SetNDC();
    tVTX7MC->SetTextColor(2);
    tVTX7MC->Draw();
    
    if(histoCorelationMC){
        TRK_SPD3D_VtxMC->cd(5);
        //corrContrSPDClusters->cd(2);
        hntrksSPDvsSPDclsMC->SetMarkerStyle(20);
        hntrksSPDvsSPDclsMC->GetXaxis()->SetRangeUser(0.,400);
        hntrksSPDvsSPDclsMC->GetYaxis()->SetRangeUser(0.,4500.);
        hntrksSPDvsSPDclsMC->SetTitle("ncontributors SPD3D vs number of cluster SPD - MC");
        hntrksSPDvsSPDclsMC->Draw();
        
        TRK_SPD3D_VtxMC->cd(6);
        //		corrContrSPDClusters->cd(3);
        hntrksZvsSPDclsMC->SetMarkerStyle(20);
        hntrksZvsSPDclsMC->GetXaxis()->SetRangeUser(0.,150.);
        hntrksZvsSPDclsMC->GetYaxis()->SetRangeUser(0.,1000.);
        hntrksZvsSPDclsMC->SetTitle("ncontributors SPDZ vs number of cluster SPD - MC");
        hntrksZvsSPDclsMC->Draw();
    }	
    TRK_SPD3D_VtxMC->SaveAs("vertexMC.pdf");
    pdfFileNames+=" vertexMC.pdf";
    delete fxMC;
    delete fyMC;
    delete fzMC;

}

//_______________________________________________________________________
void SPDVertexPileup(TFile *fildat, TFile *filMC, TCanvas **&clist, Int_t &cnum){
    
    gROOT->SetStyle("Plain");
    gStyle->SetPalette(1);
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(111);
    clist= new TCanvas*[1];//array of pointers to TCanvases

    TDirectoryFile *vtxdata = (TDirectoryFile*)fildat->Get("Vertex_Performance");
    if(!vtxdata){
        Printf("Vertex directory not found... check!");
        cnum=-1;
        return;
    }
    TList *fListData = (TList*)vtxdata->Get("cOutputVtxESD");
    if(!fListData){
        Printf("cOutputVtxESD TList not found... check!");
        cnum=-1;
        return;
    }
    TH1F *zVtxSPDData = (TH1F*)fListData->FindObject("fhSPDVertexZ");
    TH1F *zVtxSPDpilData = (TH1F*)fListData->FindObject("fhSPDVertexZPile");
    
    TDirectoryFile *vtxMC = (TDirectoryFile*)filMC->Get("Vertex_Performance");
    if(!vtxMC){
        Printf("Vertex directory not found... check!");
        cnum=-1;
        return;
    }
    TList *fListMC = (TList*)vtxMC->Get("cOutputVtxESD");
    if(!fListMC){
        Printf("cOutputVtxESD TList not found... check!");
        cnum=-1;
        return;
    }
    
    cnum = 1; // number of canvases
    TH1F *zVtxSPDMC = (TH1F*)fListMC->FindObject("fhSPDVertexZ");
    TH1F *zVtxSPDpilMC = (TH1F*)fListMC->FindObject("fhSPDVertexZPile");

    Char_t message[50];
   
    TString ctitle = GetRunNumber()+"SPD vertex pileup";
    TCanvas *pileup = new TCanvas("pileup",ctitle,1200,800);
    clist[0]=pileup;
    pileup->Divide(2,2);
    pileup->cd(1);
    zVtxSPDData->SetTitle("SPDVertex z - DATA");
    zVtxSPDData->Draw();
    pileup->cd(2);
    zVtxSPDpilData->Draw();
    zVtxSPDpilData->SetTitle("SPDVertexPile z - DATA");
    Double_t ratio=zVtxSPDpilData->GetEntries()/zVtxSPDData->GetEntries();
    sprintf(message,"pileups/vertices = %f",ratio);
    TText *txt1 = new TText(0.5,0.8,message);
    txt1->SetNDC();
    txt1->Draw();
    
    pileup->cd(3);
    zVtxSPDMC->SetTitle("SPDVertex z - MC");
    zVtxSPDMC->Draw();
    pileup->cd(4);
    zVtxSPDpilMC->SetTitle("SPDVertexPile z - MC");
    zVtxSPDpilMC->Draw();
    ratio=zVtxSPDpilMC->GetEntries()/zVtxSPDMC->GetEntries();
    sprintf(message,"pileups/vertices = %f",ratio);
    TText *txt2 = new TText(0.5,0.8,message);
    txt2->SetNDC();
    txt2->Draw();

    pileup->SaveAs("SPD_pileup.pdf");
    pdfFileNames+=" SPD_pileup.pdf";
}

//_______________________________________________________________________
void PlotSPD(TFile *fildat, TFile *filMC, TCanvas **&clist, Int_t &cnum){

  gROOT->SetStyle("Plain");
  gStyle->SetPalette(1);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(111);
  clist= new TCanvas* [2];//array of pointers to TCanvases

   TDirectoryFile *spddata = (TDirectoryFile*)fildat->Get("SPD_Performance");
    if(!spddata){
        Printf("SPD_Performance directory not found... check!");
        cnum=-1;
        return;
    }
 
    spddata->cd();
   TList *fListData = (TList*)spddata->Get("coutput1");
    if(!fListData){
        Printf("cOutput1 TList not found... check!");
        cnum=-1;
        return;
    }
    
   TString fTitleData = "Data";
   TString fTitleMc = "MC";

    cnum=2; // number of canvases

  Double_t nevtsData = ((TH1I*)(fListData->FindObject("hEventsProcessed")))->GetEntries();
  printf("   #events in %s : %f \n",fTitleData.Data(),nevtsData);

  TDirectoryFile *spdmc = (TDirectoryFile*)filMC->Get("SPD_Performance");
  spdmc->cd();
  TList *fListMc = (TList*)spdmc->Get("coutput1");
  Double_t nevtsMc = ((TH1I*)(fListMc->FindObject("hEventsProcessed")))->GetEntries();
  printf("   #events in %s : %f \n",fTitleMc.Data(),nevtsMc);

  TH2F *trackData = (TH2F*)fListData->FindObject("hSPDphivsSPDeta");
  trackData->SetTitle(Form("%s %s",trackData->GetTitle(),fTitleData.Data()));
  TH1D *trackDataPhi = trackData->ProjectionY();
  if(!trackDataPhi) printf("NO 1 \n");

  trackDataPhi->SetTitle("Tracklets vs Phi");
  TH1D *trackDataEta = trackData->ProjectionX();
  if(!trackDataEta) printf("NO 2 \n");
  trackDataEta->SetTitle("Tracklets vs eta");

  TH1F etaData, phiData;
  trackDataEta->Copy(etaData);
  trackDataPhi->Copy(phiData);

  TH1F etaFrac, phiFrac, mcEta, mcPhi;
  trackDataEta->Copy(etaFrac);
  trackDataPhi->Copy(phiFrac);

  TH2F *trackMc = (TH2F*)fListMc->FindObject("hSPDphivsSPDeta");
    trackMc->SetTitle(Form("%s %s",trackMc->GetTitle(),fTitleMc.Data()));

  TString ctitle = GetRunNumber()+"tracklets";
  TCanvas *tracklets = new TCanvas("tracklets",ctitle,1200,600);
  clist[0]=tracklets;
  tracklets->Divide(2,1);
  tracklets->cd(1);
  tracklets->cd(1)->SetRightMargin(0.15);
  trackData->SetTitle(Form("Run %d - DATA",gRunNumber));
  trackData->DrawCopy("colz");
  tracklets->cd(2);
  tracklets->cd(2)->SetRightMargin(0.15);
  trackMc->SetTitle(Form("Run %d - MC",gRunNumberMC));
  TH1D *h = (TH1D*)trackMc->DrawCopy("colz");
  tracklets->SaveAs("SPDtracklets.pdf");
  pdfFileNames+=" SPDtracklets.pdf";

  TH1D *trackMcPhi = trackMc->ProjectionY();
  trackMcPhi->SetTitle(Form("%s",h->GetTitle()));
  TH1D *trackMcEta = trackMc->ProjectionX();
  trackMcEta->SetTitle(Form("%s",h->GetTitle()));

  TH1F etaMc, phiMc;
  trackMcEta->Copy(etaMc);
  trackMcPhi->Copy(phiMc);

  trackMcEta->Copy(mcEta);
  trackMcPhi->Copy(mcPhi);

  etaFrac.Scale(1./etaFrac.GetEntries());
  mcEta.Scale(1./mcEta.GetEntries());
  etaFrac.Add(&mcEta,-1);
  etaFrac.Divide(&mcEta);

  if(phiFrac.GetEntries()>0) phiFrac.Scale(1./phiFrac.GetEntries());
  if(mcPhi.GetEntries()>0) mcPhi.Scale(1./mcPhi.GetEntries());
  phiFrac.Add(&mcPhi,-1);
  if(mcPhi.GetEntries()>0) phiFrac.Divide(&mcPhi);

  ctitle = GetRunNumber()+"tracklets and ratios vs eta and phi";
  TCanvas *track = new TCanvas("track",ctitle,1200,1200);
  clist[1]=track;
  track->Divide(2,2);
  track->cd(1);
  phiData.SetLineColor(kRed);
  phiData.SetLineWidth(2);
  phiData.Scale(1./phiData.GetEntries());
  phiData.DrawCopy();
  if(phiMc.GetEntries()>0) phiMc.Scale(1./phiMc.GetEntries());
  TLatex* tphi=new TLatex(0.25,0.85,Form("Red = %d DATA; Blue = %d MC",gRunNumber,gRunNumberMC));
  tphi->SetNDC();
  tphi->SetTextSize(0.04);
   tphi->SetTextColor(1);
  tphi->Draw();  
  phiMc.DrawCopy("same");
  track->cd(2);
  etaData.SetLineColor(kRed);
  etaData.SetLineWidth(2);
  etaData.Scale(1./etaData.GetEntries());
  etaData.DrawCopy();
  if(etaMc.GetEntries()>0) etaMc.Scale(1./etaMc.GetEntries());
    TLatex* tphi2=new TLatex(0.15,0.35,Form("Red = %d DATA; Blue = %d MC",gRunNumber,gRunNumberMC));
    tphi2->SetNDC();
    tphi2->SetTextSize(0.04);
    tphi2->SetTextColor(1);
  tphi2->Draw();
  etaMc.DrawCopy("same");
  track->cd(3);
  phiFrac.SetLineColor(1);
  TLatex* tratio=new TLatex(0.2,0.85,Form("Run %d DATA / Run %d MC",gRunNumber,gRunNumberMC));
  tratio->SetNDC();
  tratio->SetTextColor(1);
  TLatex* tratio2=new TLatex(0.2,0.80,"(DATA - MC)/MC");
    tratio2->SetNDC();
    tratio2->SetTextColor(1);
  phiFrac.DrawCopy();
  tratio->Draw();
    tratio2->Draw();
  track->cd(4);
  etaFrac.SetLineColor(1);
  etaFrac.DrawCopy();  
  tratio->Draw();
    tratio2->Draw();
  track->SaveAs("SPD_eta_phi.pdf");
  pdfFileNames+=" SPD_eta_phi.pdf";
}

//_______________________________________________________________________
Bool_t PlotITSTPCMatchingEff(TFile *f, TFile *fMC, TCanvas**& clist,Int_t& cnum) {

    cnum=1;
  clist = new TCanvas*[1];

  //  clist = new TCanvas* [1];
  TString ctitle = GetRunNumber()+"ITS-TPC match - DATA";
  TCanvas* cITSTPCmatch = new TCanvas("cITSTPCmatch",ctitle,10,10,1200,600);
  clist[0]=cITSTPCmatch;
  // cITSTPCmatch->Divide(2,1);
  // cITSTPCmatch->cd(1);
  gPad->SetGridy();
  gPad->SetLogx();
  // cITSTPCmatch->cd(2);
  // gPad->SetGridy();
  // gPad->SetLogx();

  //  clist = cITSTPCmatch;

  if(!f) return kFALSE;

  TList *list=0;
  TList *listSPD=0;
  TDirectoryFile *dir=0;

  // count active SPD HSs
  dir=(TDirectoryFile*)f->GetDirectory("SPD_Performance");
  if(dir) listSPD = (TList*)dir->Get("coutput1");
  if(!dir) return kFALSE;

  Float_t spdFrac[2]={0.,0.};
  TH1F *hnHSsSPD=new TH1F("hnHSsSPD","Active HSs in SPD layers 1 and 2; layer; HSs",2,0.5,2.5);
  if(listSPD) {
    //listSPD->Print();
    TH1F *hFiredChip = (TH1F*)listSPD->FindObject("hFiredChip");
    Int_t nHSsInner=0,nHSsOuter=0;
    for(Int_t i=0;i<400;i++) if(hFiredChip->GetBinContent(i)>0) nHSsInner++;
    for(Int_t i=400;i<1200;i++) if(hFiredChip->GetBinContent(i)>0) nHSsOuter++;
    nHSsInner = (Int_t)(nHSsInner/10);
    nHSsOuter = (Int_t)(nHSsOuter/10);
    hnHSsSPD->SetBinContent(1,nHSsInner);
    hnHSsSPD->SetBinContent(2,nHSsOuter);
    spdFrac[0]=(Float_t)nHSsInner/40.;
    spdFrac[1]=(Float_t)nHSsOuter/80.;
  }
  TGraph *spdFrac0=new TGraph(1);
  spdFrac0->SetPoint(0,0.08,spdFrac[0]);
  spdFrac0->SetMarkerColor(1); spdFrac0->SetMarkerStyle(20);
  TGraph *spdFrac1=new TGraph(1);
  spdFrac1->SetPoint(0,0.08,spdFrac[1]);
  spdFrac1->SetMarkerColor(1); spdFrac1->SetMarkerStyle(24);
  TLegend *l2=new TLegend(0.1,0.62,0.5,0.93);
  l2->SetBorderSize(1);
  l2->AddEntry(spdFrac0,"Frac. active SPD0","p");
  l2->AddEntry(spdFrac1,"Frac. active SPD1","p");

  //
  // Efficiencies for CENTRAL
  //
  dir=(TDirectoryFile*)f->GetDirectory("ITS_Performance");
  if(dir) list = (TList*)dir->Get("cOutputITS");
  if(!list) return kFALSE;

  TH1F *fHistPtTPCInAcc = (TH1F*)list->FindObject("fHistPtTPCInAcc");
  TH1F *fHistPtITSMI6InAcc = (TH1F*)list->FindObject("fHistPtITSMI6InAcc");
  TH1F *fHistPtITSMI5InAcc = (TH1F*)list->FindObject("fHistPtITSMI5InAcc");
  TH1F *fHistPtITSMI4InAcc = (TH1F*)list->FindObject("fHistPtITSMI4InAcc");
  TH1F *fHistPtITSMI3InAcc = (TH1F*)list->FindObject("fHistPtITSMI3InAcc");
  TH1F *fHistPtITSMI2InAcc = (TH1F*)list->FindObject("fHistPtITSMI2InAcc");
  TH1F *fHistPtITSMISPDInAcc = (TH1F*)list->FindObject("fHistPtITSMISPDInAcc");
  TH1F *fHistPtITSMIoneSPDInAcc = (TH1F*)list->FindObject("fHistPtITSMIoneSPDInAcc");
  TH1F *fHistPtITSTPCsel = (TH1F*)list->FindObject("fHistPtITSTPCsel");
  TH1F *fHistPtITSMIge2InAcc = (TH1F*)fHistPtITSMI6InAcc->Clone("fHistPtITSMIge2InAcc");
  fHistPtITSMIge2InAcc->Add(fHistPtITSMI5InAcc);
  fHistPtITSMIge2InAcc->Add(fHistPtITSMI4InAcc);
  fHistPtITSMIge2InAcc->Add(fHistPtITSMI3InAcc);
  fHistPtITSMIge2InAcc->Add(fHistPtITSMI2InAcc);


  TLegend *l3=new TLegend(0.5,0.62,0.95,0.93);
  l3->SetBorderSize(1);
  //cITSTPCmatch->cd(1);
  fHistPtITSMIge2InAcc->SetTitle("Fraction of prolonged tracks with N ITS points - DATA");
  fHistPtITSMIge2InAcc->SetYTitle("ITS+TPC / TPC");
  fHistPtITSMIge2InAcc->Divide(fHistPtITSMIge2InAcc,fHistPtTPCInAcc,1,1,"B");
  fHistPtITSMIge2InAcc->SetMaximum(1.6);
  fHistPtITSMIge2InAcc->SetMinimum(0);
  fHistPtITSMIge2InAcc->GetXaxis()->SetRangeUser(0.1,30);
  fHistPtITSMIge2InAcc->Draw();
  l3->AddEntry(fHistPtITSMIge2InAcc,">=2 cls","l");
  fHistPtITSMI6InAcc->Divide(fHistPtITSMI6InAcc,fHistPtTPCInAcc,1,1,"B");
  fHistPtITSMI6InAcc->SetLineColor(2);
  l3->AddEntry(fHistPtITSMI6InAcc,"6 cls","l");
  fHistPtITSMI6InAcc->Draw("same");
  fHistPtITSMI5InAcc->Divide(fHistPtITSMI5InAcc,fHistPtTPCInAcc,1,1,"B");
  fHistPtITSMI5InAcc->SetLineColor(3);
  l3->AddEntry(fHistPtITSMI5InAcc,"5 cls","l");
  fHistPtITSMI5InAcc->Draw("same");
  fHistPtITSMI4InAcc->Divide(fHistPtITSMI4InAcc,fHistPtTPCInAcc,1,1,"B");
  fHistPtITSMI4InAcc->SetLineColor(4);
  l3->AddEntry(fHistPtITSMI4InAcc,"4 cls","l");
  fHistPtITSMI4InAcc->Draw("same");
  fHistPtITSMI3InAcc->Divide(fHistPtITSMI3InAcc,fHistPtTPCInAcc,1,1,"B");
  fHistPtITSMI3InAcc->SetLineColor(6);
  l3->AddEntry(fHistPtITSMI3InAcc,"3 cls","l");
  fHistPtITSMI3InAcc->Draw("same");
  fHistPtITSMI2InAcc->Divide(fHistPtITSMI2InAcc,fHistPtTPCInAcc,1,1,"B");
  fHistPtITSMI2InAcc->SetLineColor(7);
  l3->AddEntry(fHistPtITSMI2InAcc,"2 cls","l");
  fHistPtITSMI2InAcc->Draw("same");
  fHistPtITSMISPDInAcc->Divide(fHistPtITSMISPDInAcc,fHistPtTPCInAcc,1,1,"B");
  fHistPtITSMISPDInAcc->SetLineColor(9);
  l3->AddEntry(fHistPtITSMISPDInAcc,"2SPD + any","l");
  fHistPtITSMISPDInAcc->Draw("same");
  fHistPtITSMIoneSPDInAcc->Divide(fHistPtITSMIoneSPDInAcc,fHistPtTPCInAcc,1,1,"B");
  fHistPtITSMIoneSPDInAcc->SetLineColor(15);
  l3->AddEntry(fHistPtITSMIoneSPDInAcc,">=1SPD + any","l");
  fHistPtITSMIoneSPDInAcc->Draw("same");
  fHistPtITSTPCsel->Divide(fHistPtITSTPCsel,fHistPtTPCInAcc,1,1,"B");
  fHistPtITSTPCsel->SetLineColor(kAzure+1);
  l3->AddEntry(fHistPtITSTPCsel,">=1SPD + any + d_{0} cut","l");
  fHistPtITSTPCsel->Draw("same");
  fHistPtITSMIge2InAcc->Draw("same");
  l3->Draw();
  l2->Draw();
  spdFrac0->Draw("p");
  spdFrac1->Draw("p");

    // TOFITS
    
    ctitle = GetRunNumber()+"ITS-TOF match - DATA";
    TCanvas* cITSTOFmatch = new TCanvas("cITSTOFmatch",ctitle,10,10,1200,600);
    clist[1]=cITSTOFmatch;
    // cITSTPCmatch->Divide(2,1);
    // cITSTPCmatch->cd(1);
    gPad->SetGridy();
    gPad->SetLogx();
    
    TH1F *fHistPtTPCInAccTOF = (TH1F*)list->FindObject("fHistPtTPCInAccTOFbc0");
    TH1F *fHistPtITSMI6InAccTOF = (TH1F*)list->FindObject("fHistPtITSMI6InAccTOFbc0");
    TH1F *fHistPtITSMI5InAccTOF = (TH1F*)list->FindObject("fHistPtITSMI5InAccTOFbc0");
    TH1F *fHistPtITSMI4InAccTOF = (TH1F*)list->FindObject("fHistPtITSMI4InAccTOFbc0");
    TH1F *fHistPtITSMI3InAccTOF = (TH1F*)list->FindObject("fHistPtITSMI3InAccTOFbc0");
    TH1F *fHistPtITSMI2InAccTOF = (TH1F*)list->FindObject("fHistPtITSMI2InAccTOFbc0");
    TH1F *fHistPtITSMISPDInAccTOF = (TH1F*)list->FindObject("fHistPtITSMISPDInAccTOFbc0");
    TH1F *fHistPtITSMIoneSPDInAccTOF = (TH1F*)list->FindObject("fHistPtITSMIoneSPDInAccTOFbc0");
    TH1F *fHistPtITSTPCselTOF = (TH1F*)list->FindObject("fHistPtITSTPCselTOFbc0");
    TH1F *fHistPtITSMIge2InAccTOF = (TH1F*)fHistPtITSMI6InAccTOF->Clone("fHistPtITSMIge2InAccTOF");
    fHistPtITSMIge2InAccTOF->Add(fHistPtITSMI5InAccTOF);
    fHistPtITSMIge2InAccTOF->Add(fHistPtITSMI4InAccTOF);
    fHistPtITSMIge2InAccTOF->Add(fHistPtITSMI3InAccTOF);
    fHistPtITSMIge2InAccTOF->Add(fHistPtITSMI2InAccTOF);
    
    
    TLegend *l3t=new TLegend(0.5,0.62,0.95,0.93);
    l3t->SetBorderSize(1);
    //cITSTPCmatch->cd(1);
    fHistPtITSMIge2InAccTOF->SetTitle("Fraction of prolonged tracks (TOF) with N ITS points - DATA");
    fHistPtITSMIge2InAccTOF->SetYTitle("ITS+TOF / TPC");
    fHistPtITSMIge2InAccTOF->Divide(fHistPtITSMIge2InAccTOF,fHistPtTPCInAccTOF,1,1,"B");
    fHistPtITSMIge2InAccTOF->SetMaximum(1.6);
    fHistPtITSMIge2InAccTOF->SetMinimum(0);
    fHistPtITSMIge2InAccTOF->GetXaxis()->SetRangeUser(0.1,30);
    fHistPtITSMIge2InAccTOF->Draw();
    l3t->AddEntry(fHistPtITSMIge2InAccTOF,">=2 cls","l");
    fHistPtITSMI6InAccTOF->Divide(fHistPtITSMI6InAccTOF,fHistPtTPCInAccTOF,1,1,"B");
    fHistPtITSMI6InAccTOF->SetLineColor(2);
    l3t->AddEntry(fHistPtITSMI6InAccTOF,"6 cls","l");
    fHistPtITSMI6InAccTOF->Draw("same");
    fHistPtITSMI5InAccTOF->Divide(fHistPtITSMI5InAccTOF,fHistPtTPCInAccTOF,1,1,"B");
    fHistPtITSMI5InAccTOF->SetLineColor(3);
    l3t->AddEntry(fHistPtITSMI5InAccTOF,"5 cls","l");
    fHistPtITSMI5InAccTOF->Draw("same");
    fHistPtITSMI4InAccTOF->Divide(fHistPtITSMI4InAccTOF,fHistPtTPCInAccTOF,1,1,"B");
    fHistPtITSMI4InAccTOF->SetLineColor(4);
    l3t->AddEntry(fHistPtITSMI4InAccTOF,"4 cls","l");
    fHistPtITSMI4InAccTOF->Draw("same");
    fHistPtITSMI3InAccTOF->Divide(fHistPtITSMI3InAccTOF,fHistPtTPCInAccTOF,1,1,"B");
    fHistPtITSMI3InAccTOF->SetLineColor(6);
    l3t->AddEntry(fHistPtITSMI3InAccTOF,"3 cls","l");
    fHistPtITSMI3InAccTOF->Draw("same");
    fHistPtITSMI2InAccTOF->Divide(fHistPtITSMI2InAccTOF,fHistPtTPCInAccTOF,1,1,"B");
    fHistPtITSMI2InAccTOF->SetLineColor(7);
    l3t->AddEntry(fHistPtITSMI2InAccTOF,"2 cls","l");
    fHistPtITSMI2InAccTOF->Draw("same");
    fHistPtITSMISPDInAccTOF->Divide(fHistPtITSMISPDInAccTOF,fHistPtTPCInAccTOF,1,1,"B");
    fHistPtITSMISPDInAccTOF->SetLineColor(9);
    l3t->AddEntry(fHistPtITSMISPDInAccTOF,"2SPD + any","l");
    fHistPtITSMISPDInAccTOF->Draw("same");
    fHistPtITSMIoneSPDInAccTOF->Divide(fHistPtITSMIoneSPDInAccTOF,fHistPtTPCInAccTOF,1,1,"B");
    fHistPtITSMIoneSPDInAccTOF->SetLineColor(15);
    l3t->AddEntry(fHistPtITSMIoneSPDInAccTOF,">=1SPD + any","l");
    fHistPtITSMIoneSPDInAccTOF->Draw("same");
    fHistPtITSTPCselTOF->Divide(fHistPtITSTPCselTOF,fHistPtTPCInAccTOF,1,1,"B");
    fHistPtITSTPCselTOF->SetLineColor(kAzure+1);
    l3t->AddEntry(fHistPtITSTPCselTOF,">=1SPD + any + d_{0} cut","l");
    fHistPtITSTPCselTOF->Draw("same");
    fHistPtITSMIge2InAccTOF->Draw("same");
    l3t->Draw();
    l2->Draw();
    spdFrac0->Draw("p");
    spdFrac1->Draw("p");

  cITSTPCmatch->SaveAs("TPCITSmatching.pdf");
  pdfFileNames+=" TPCITSmatching.pdf";

// MC file
    clist = new TCanvas*[1];
    
    //  clist = new TCanvas* [1];
    TString ctitleMC = GetRunNumber()+"ITS-TPC match -- MC";
    TCanvas* cITSTPCmatchMC = new TCanvas("cITSTPCmatchMC",ctitleMC,10,10,1200,600);
    clist[0]=cITSTPCmatchMC;
    gPad->SetGridy();
    gPad->SetLogx();
    
    //  clist = cITSTPCmatch;
    
    if(!fMC) return kFALSE;
    
    TList *listMC=0;
    TList *listSPDMC=0;
    TDirectoryFile *dirMC=0;
    
    // count active SPD HSs
    dirMC=(TDirectoryFile*)fMC->GetDirectory("SPD_Performance");
    if(dirMC) listSPDMC = (TList*)dirMC->Get("coutput1");
    if(!dirMC) return kFALSE;
    
    Float_t spdFracMC[2]={0.,0.};
    TH1F *hnHSsSPDMC=new TH1F("hnHSsSPDMC","Active HSs in SPD layers 1 and 2; layer; HSs",2,0.5,2.5);
    if(listSPDMC) {
        //listSPD->Print();
        TH1F *hFiredChipMC = (TH1F*)listSPDMC->FindObject("hFiredChip");
        Int_t nHSsInnerMC=0,nHSsOuterMC=0;
        for(Int_t i=0;i<400;i++) if(hFiredChipMC->GetBinContent(i)>0) nHSsInnerMC++;
        for(Int_t i=400;i<1200;i++) if(hFiredChipMC->GetBinContent(i)>0) nHSsOuterMC++;
        nHSsInnerMC = (Int_t)(nHSsInnerMC/10);
        nHSsOuterMC = (Int_t)(nHSsOuterMC/10);
        hnHSsSPDMC->SetBinContent(1,nHSsInnerMC);
        hnHSsSPDMC->SetBinContent(2,nHSsOuterMC);
        spdFracMC[0]=(Float_t)nHSsInnerMC/40.;
        spdFracMC[1]=(Float_t)nHSsOuterMC/80.;
    }
    TGraph *spdFrac0MC=new TGraph(1);
    spdFrac0MC->SetPoint(0,0.08,spdFracMC[0]);
    spdFrac0MC->SetMarkerColor(1); spdFrac0MC->SetMarkerStyle(20);
    TGraph *spdFrac1MC=new TGraph(1);
    spdFrac1MC->SetPoint(0,0.08,spdFracMC[1]);
    spdFrac1MC->SetMarkerColor(1); spdFrac1MC->SetMarkerStyle(24);
    TLegend *l2MC=new TLegend(0.1,0.62,0.5,0.93);
    l2MC->SetBorderSize(1);
    l2MC->AddEntry(spdFrac0MC,"Frac. active SPD0","p");
    l2MC->AddEntry(spdFrac1MC,"Frac. active SPD1","p");
    

    //
    // Efficiencies for CENTRAL
    //
    dirMC=(TDirectoryFile*)fMC->GetDirectory("ITS_Performance");
    if(dirMC) listMC = (TList*)dirMC->Get("cOutputITS");
    if(!listMC) return kFALSE;

    TH1F *fHistPtTPCInAccMC = (TH1F*)listMC->FindObject("fHistPtTPCInAcc");
    TH1F *fHistPtITSMI6InAccMC = (TH1F*)listMC->FindObject("fHistPtITSMI6InAcc");
    TH1F *fHistPtITSMI5InAccMC = (TH1F*)listMC->FindObject("fHistPtITSMI5InAcc");
    TH1F *fHistPtITSMI4InAccMC = (TH1F*)listMC->FindObject("fHistPtITSMI4InAcc");
    TH1F *fHistPtITSMI3InAccMC = (TH1F*)listMC->FindObject("fHistPtITSMI3InAcc");
    TH1F *fHistPtITSMI2InAccMC = (TH1F*)listMC->FindObject("fHistPtITSMI2InAcc");
    TH1F *fHistPtITSMISPDInAccMC = (TH1F*)listMC->FindObject("fHistPtITSMISPDInAcc");
    TH1F *fHistPtITSMIoneSPDInAccMC = (TH1F*)listMC->FindObject("fHistPtITSMIoneSPDInAcc");
    TH1F *fHistPtITSTPCselMC = (TH1F*)listMC->FindObject("fHistPtITSTPCsel");
    TH1F *fHistPtITSMIge2InAccMC = (TH1F*)fHistPtITSMI6InAccMC->Clone("fHistPtITSMIge2InAccMC");
    fHistPtITSMIge2InAccMC->Add(fHistPtITSMI5InAccMC);
    fHistPtITSMIge2InAccMC->Add(fHistPtITSMI4InAccMC);
    fHistPtITSMIge2InAccMC->Add(fHistPtITSMI3InAccMC);
    fHistPtITSMIge2InAccMC->Add(fHistPtITSMI2InAccMC);
    
    
    TLegend *l3MC=new TLegend(0.5,0.62,0.95,0.93);
    l3MC->SetBorderSize(1);
    //cITSTPCmatch->cd(1);
    fHistPtITSMIge2InAccMC->SetTitle("Fraction of prolonged tracks with N ITS points - MC");
    fHistPtITSMIge2InAccMC->SetYTitle("ITS+TPC / TPC");
    fHistPtITSMIge2InAccMC->Divide(fHistPtITSMIge2InAccMC,fHistPtTPCInAccMC,1,1,"B");
    fHistPtITSMIge2InAccMC->SetMaximum(1.6);
    fHistPtITSMIge2InAccMC->SetMinimum(0);
    fHistPtITSMIge2InAccMC->GetXaxis()->SetRangeUser(0.1,30);
    fHistPtITSMIge2InAccMC->Draw();

    printf("accedo ai dati MC ME 1 \n");

    l3MC->AddEntry(fHistPtITSMIge2InAccMC,">=2 cls","l");
    fHistPtITSMI6InAccMC->Divide(fHistPtITSMI6InAccMC,fHistPtTPCInAccMC,1,1,"B");
    fHistPtITSMI6InAccMC->SetLineColor(2);
    l3MC->AddEntry(fHistPtITSMI6InAccMC,"6 cls","l");
    fHistPtITSMI6InAccMC->Draw("same");
    fHistPtITSMI5InAccMC->Divide(fHistPtITSMI5InAccMC,fHistPtTPCInAccMC,1,1,"B");
    fHistPtITSMI5InAccMC->SetLineColor(3);
    l3MC->AddEntry(fHistPtITSMI5InAccMC,"5 cls","l");
    fHistPtITSMI5InAccMC->Draw("same");
    fHistPtITSMI4InAccMC->Divide(fHistPtITSMI4InAccMC,fHistPtTPCInAccMC,1,1,"B");
    fHistPtITSMI4InAccMC->SetLineColor(4);
    l3MC->AddEntry(fHistPtITSMI4InAccMC,"4 cls","l");
    fHistPtITSMI4InAccMC->Draw("same");
    fHistPtITSMI3InAccMC->Divide(fHistPtITSMI3InAccMC,fHistPtTPCInAccMC,1,1,"B");
    fHistPtITSMI3InAccMC->SetLineColor(6);
    l3MC->AddEntry(fHistPtITSMI3InAccMC,"3 cls","l");
    fHistPtITSMI3InAccMC->Draw("same");
    fHistPtITSMI2InAccMC->Divide(fHistPtITSMI2InAccMC,fHistPtTPCInAccMC,1,1,"B");
    fHistPtITSMI2InAccMC->SetLineColor(7);
    l3MC->AddEntry(fHistPtITSMI2InAccMC,"2 cls","l");
    fHistPtITSMI2InAccMC->Draw("same");
    fHistPtITSMISPDInAccMC->Divide(fHistPtITSMISPDInAccMC,fHistPtTPCInAccMC,1,1,"B");
    fHistPtITSMISPDInAccMC->SetLineColor(9);
    l3MC->AddEntry(fHistPtITSMISPDInAccMC,"2SPD + any","l");
    fHistPtITSMISPDInAccMC->Draw("same");
    fHistPtITSMIoneSPDInAccMC->Divide(fHistPtITSMIoneSPDInAccMC,fHistPtTPCInAccMC,1,1,"B");
    fHistPtITSMIoneSPDInAccMC->SetLineColor(15);
    l3MC->AddEntry(fHistPtITSMIoneSPDInAccMC,">=1SPD + any","l");
    fHistPtITSMIoneSPDInAccMC->Draw("same");
    fHistPtITSTPCselMC->Divide(fHistPtITSTPCselMC,fHistPtTPCInAccMC,1,1,"B");
    fHistPtITSTPCselMC->SetLineColor(kAzure+1);
    l3MC->AddEntry(fHistPtITSTPCselMC,">=1SPD + any + d_{0} cut","l");
    fHistPtITSTPCselMC->Draw("same");
    fHistPtITSMIge2InAccMC->Draw("same");
    l3MC->Draw();
    l2MC->Draw();
    spdFrac0MC->Draw("p");
    spdFrac1MC->Draw("p");

    // TOFITS
    
    ctitle = GetRunNumber()+"ITS-TOF match - MC";
    TCanvas* cITSTOFmatchMC = new TCanvas("cITSTOFmatchMC",ctitle,10,10,1200,600);
    clist[1]=cITSTOFmatchMC;
    // cITSTPCmatch->Divide(2,1);
    // cITSTPCmatch->cd(1);
    gPad->SetGridy();
    gPad->SetLogx();
    
    TH1F *fHistPtTPCInAccTOFMC = (TH1F*)listMC->FindObject("fHistPtTPCInAccTOFbc0");
    TH1F *fHistPtITSMI6InAccTOFMC = (TH1F*)listMC->FindObject("fHistPtITSMI6InAccTOFbc0");
    TH1F *fHistPtITSMI5InAccTOFMC = (TH1F*)listMC->FindObject("fHistPtITSMI5InAccTOFbc0");
    TH1F *fHistPtITSMI4InAccTOFMC = (TH1F*)listMC->FindObject("fHistPtITSMI4InAccTOFbc0");
    TH1F *fHistPtITSMI3InAccTOFMC = (TH1F*)listMC->FindObject("fHistPtITSMI3InAccTOFbc0");
    TH1F *fHistPtITSMI2InAccTOFMC = (TH1F*)listMC->FindObject("fHistPtITSMI2InAccTOFbc0");
    TH1F *fHistPtITSMISPDInAccTOFMC = (TH1F*)listMC->FindObject("fHistPtITSMISPDInAccTOFbc0");
    TH1F *fHistPtITSMIoneSPDInAccTOFMC = (TH1F*)listMC->FindObject("fHistPtITSMIoneSPDInAccTOFbc0");
    TH1F *fHistPtITSTPCselTOFMC = (TH1F*)listMC->FindObject("fHistPtITSTPCselTOFbc0");
    TH1F *fHistPtITSMIge2InAccTOFMC = (TH1F*)fHistPtITSMI6InAccTOFMC->Clone("fHistPtITSMIge2InAccTOFMC");
    fHistPtITSMIge2InAccTOFMC->Add(fHistPtITSMI5InAccTOFMC);
    fHistPtITSMIge2InAccTOFMC->Add(fHistPtITSMI4InAccTOFMC);
    fHistPtITSMIge2InAccTOFMC->Add(fHistPtITSMI3InAccTOFMC);
    fHistPtITSMIge2InAccTOFMC->Add(fHistPtITSMI2InAccTOFMC);
    
    
    TLegend *l3tMC=new TLegend(0.5,0.62,0.95,0.93);
    l3tMC->SetBorderSize(1);
    //cITSTPCmatch->cd(1);
    fHistPtITSMIge2InAccTOFMC->SetTitle("Fraction of prolonged tracks (TOF) with N ITS points - MC");
    fHistPtITSMIge2InAccTOFMC->SetYTitle("ITS+TOF / TPC");
    fHistPtITSMIge2InAccTOFMC->Divide(fHistPtITSMIge2InAccTOFMC,fHistPtTPCInAccTOFMC,1,1,"B");
    fHistPtITSMIge2InAccTOFMC->SetMaximum(1.6);
    fHistPtITSMIge2InAccTOFMC->SetMinimum(0);
    fHistPtITSMIge2InAccTOFMC->GetXaxis()->SetRangeUser(0.1,30);
    fHistPtITSMIge2InAccTOFMC->Draw();
    l3tMC->AddEntry(fHistPtITSMIge2InAccTOFMC,">=2 cls","l");
    fHistPtITSMI6InAccTOFMC->Divide(fHistPtITSMI6InAccTOFMC,fHistPtTPCInAccTOFMC,1,1,"B");
    fHistPtITSMI6InAccTOFMC->SetLineColor(2);
    l3tMC->AddEntry(fHistPtITSMI6InAccTOFMC,"6 cls","l");
    fHistPtITSMI6InAccTOFMC->Draw("same");
    fHistPtITSMI5InAccTOFMC->Divide(fHistPtITSMI5InAccTOFMC,fHistPtTPCInAccTOFMC,1,1,"B");
    fHistPtITSMI5InAccTOFMC->SetLineColor(3);
    l3tMC->AddEntry(fHistPtITSMI5InAccTOFMC,"5 cls","l");
    fHistPtITSMI5InAccTOFMC->Draw("same");
    fHistPtITSMI4InAccTOFMC->Divide(fHistPtITSMI4InAccTOFMC,fHistPtTPCInAccTOFMC,1,1,"B");
    fHistPtITSMI4InAccTOFMC->SetLineColor(4);
    l3tMC->AddEntry(fHistPtITSMI4InAccTOFMC,"4 cls","l");
    fHistPtITSMI4InAccTOFMC->Draw("same");
    fHistPtITSMI3InAccTOFMC->Divide(fHistPtITSMI3InAccTOFMC,fHistPtTPCInAccTOFMC,1,1,"B");
    fHistPtITSMI3InAccTOFMC->SetLineColor(6);
    l3tMC->AddEntry(fHistPtITSMI3InAccTOFMC,"3 cls","l");
    fHistPtITSMI3InAccTOFMC->Draw("same");
    fHistPtITSMI2InAccTOFMC->Divide(fHistPtITSMI2InAccTOFMC,fHistPtTPCInAccTOFMC,1,1,"B");
    fHistPtITSMI2InAccTOFMC->SetLineColor(7);
    l3tMC->AddEntry(fHistPtITSMI2InAccTOFMC,"2 cls","l");
    fHistPtITSMI2InAccTOFMC->Draw("same");
    fHistPtITSMISPDInAccTOFMC->Divide(fHistPtITSMISPDInAccTOFMC,fHistPtTPCInAccTOFMC,1,1,"B");
    fHistPtITSMISPDInAccTOFMC->SetLineColor(9);
    l3tMC->AddEntry(fHistPtITSMISPDInAccTOFMC,"2SPD + any","l");
    fHistPtITSMISPDInAccTOFMC->Draw("same");
    fHistPtITSMIoneSPDInAccTOFMC->Divide(fHistPtITSMIoneSPDInAccTOFMC,fHistPtTPCInAccTOFMC,1,1,"B");
    fHistPtITSMIoneSPDInAccTOFMC->SetLineColor(15);
    l3tMC->AddEntry(fHistPtITSMIoneSPDInAccTOFMC,">=1SPD + any","l");
    fHistPtITSMIoneSPDInAccTOFMC->Draw("same");
    fHistPtITSTPCselTOFMC->Divide(fHistPtITSTPCselTOFMC,fHistPtTPCInAccTOFMC,1,1,"B");
    fHistPtITSTPCselTOFMC->SetLineColor(kAzure+1);
    l3tMC->AddEntry(fHistPtITSTPCselTOFMC,">=1SPD + any + d_{0} cut","l");
    fHistPtITSTPCselTOFMC->Draw("same");
    fHistPtITSMIge2InAccTOFMC->Draw("same");
    l3tMC->Draw();
    l2MC->Draw();
    spdFrac0MC->Draw("p");
    spdFrac1MC->Draw("p");
    
    
    cITSTPCmatchMC->SaveAs("TPCITSmatchingMC.pdf");
    pdfFileNames+=" TPCITSmatchingMC.pdf";

    cITSTOFmatch->SaveAs("TOFITSmatching.pdf");
    pdfFileNames+=" TOFITSmatching.pdf";
    cITSTOFmatchMC->SaveAs("TOFITSmatchingMC.pdf");
    pdfFileNames+=" TOFITSmatchingMC.pdf";

//
// MC/DATA efficiency ratios
//
    
    Int_t nPtBins=34;
    Float_t xPtBins[35]={0,0.025,0.05,0.075,0.1,0.125,0.15,0.175,0.2,0.225,0.25,0.275,0.3,0.325,0.35,0.375,0.4,0.5,0.6,0.7,0.8,1.0,1.5,2.,2.5,3,4,5,6,8,10,15,20,25,30};
    Float_t num[nPtBins], denom[nPtBins], ratio[nPtBins];
    Float_t ernum[nPtBins], erdenom[nPtBins], erratio[nPtBins];

// TPCITS ME ratios

    TH1F *ratio456data = (TH1F *)fHistPtITSMI6InAcc->Clone("ratio456data");
    ratio456data->Add(fHistPtITSMI5InAcc);
    ratio456data->Add(fHistPtITSMI4InAcc);

    TH1F *ratio456MC = (TH1F *)fHistPtITSMI6InAccMC->Clone("ratio456MC");
    ratio456MC->Add(fHistPtITSMI5InAccMC);
    ratio456MC->Add(fHistPtITSMI4InAccMC);

    TH1F *ratio456 = new TH1F("ratio456","ratio between MC and data TPCITS efficiency, tracks with>=4 hit; p_{t} [GeV/c]; MC eff/DATA eff",nPtBins,xPtBins);

    TH1F *ratio2any = new TH1F("ratio2any","ratio between MC and data TPCITS efficiency, tracks with 2SPD+any hit; p_{t} [GeV/c]; MC eff/DATA eff",nPtBins,xPtBins);

    
    
    for(Int_t cc=0;cc<nPtBins;cc++){                // 4-5-6 clusters
        num[cc]=ratio456MC->GetBinContent(cc);
        ernum[cc]=ratio456MC->GetBinError(cc);
        denom[cc]=ratio456data->GetBinContent(cc);
        erdenom[cc]=ratio456data->GetBinError(cc);
        if(denom[cc] > 0.){
            ratio[cc]=num[cc]/denom[cc];
            erratio[cc]=num[cc]/denom[cc] * TMath::Sqrt((ernum[cc]/num[cc])*(ernum[cc]/num[cc])+(erdenom[cc]/denom[cc])*(erdenom[cc]/denom[cc]));
//            erratio[cc]=TMath::Sqrt(TMath::Abs(num[cc]*(1-num[cc]/denom[cc])))/denom[cc];
            if(erratio[cc] == 0.) erratio[cc]=0.01;
        }
        else{
            ratio[cc]=0.0;
            erratio[cc]=0.01;
        }
        ratio456->SetBinContent(cc+1,ratio[cc]);
        ratio456->SetBinError(cc+1,erratio[cc]);
    }
    
    for(Int_t cc=0;cc<nPtBins;cc++){            // 2 SPD + any
        num[cc]=fHistPtITSMISPDInAccMC->GetBinContent(cc);
        ernum[cc]=fHistPtITSMISPDInAccMC->GetBinError(cc);
        denom[cc]=fHistPtITSMISPDInAcc->GetBinContent(cc);
        erdenom[cc]=fHistPtITSMISPDInAcc->GetBinError(cc);
        if(denom[cc] > 0.){
            ratio[cc]=num[cc]/denom[cc];
            erratio[cc]=num[cc]/denom[cc] * TMath::Sqrt((ernum[cc]/num[cc])*(ernum[cc]/num[cc])+(erdenom[cc]/denom[cc])*(erdenom[cc]/denom[cc]));
//            erratio[cc]=TMath::Sqrt(TMath::Abs(num[cc]*(1-num[cc]/denom[cc])))/denom[cc];
            if(erratio[cc] == 0.) erratio[cc]=0.01;
        }
        else{
            ratio[cc]=0.0;
            erratio[cc]=0.01;
        }
        ratio2any->SetBinContent(cc+1,ratio[cc]);
        ratio2any->SetBinError(cc+1,erratio[cc]);
    }
    
    // TOFITS ME ratios

    TH1F *ratio456dataTOF = (TH1F *)fHistPtITSMI6InAccTOF->Clone("ratio456dataTOF");
    ratio456dataTOF->Add(fHistPtITSMI5InAccTOF);
    ratio456dataTOF->Add(fHistPtITSMI4InAccTOF);
    
    TH1F *ratio456MCTOF = (TH1F *)fHistPtITSMI6InAccTOFMC->Clone("ratio456MCTOF");
    ratio456MCTOF->Add(fHistPtITSMI5InAccTOFMC);
    ratio456MCTOF->Add(fHistPtITSMI4InAccTOFMC);
    
    TH1F *ratio456TOF = new TH1F("ratio456","ratio of MC and data TOFITS efficiency, tracks with>=4 hit; p_{t} [GeV/c]; MC eff/DATA eff",nPtBins,xPtBins);

    TH1F *ratio2anyTOF = new TH1F("ratio2any","ratio of MC and data TOFITS efficiency, tracks with 2SPD+any hit; p_{t} [GeV/c]; MC eff/DATA eff",nPtBins,xPtBins);
    

    for(Int_t cc=0;cc<nPtBins;cc++){                // 4-5-6 clusters
        num[cc]=ratio456MCTOF->GetBinContent(cc);
        ernum[cc]=ratio456MCTOF->GetBinError(cc);
        denom[cc]=ratio456dataTOF->GetBinContent(cc);
        erdenom[cc]=ratio456dataTOF->GetBinError(cc);
        if(denom[cc] > 0.){
            ratio[cc]=num[cc]/denom[cc];
            erratio[cc]=num[cc]/denom[cc] * TMath::Sqrt((ernum[cc]/num[cc])*(ernum[cc]/num[cc])+(erdenom[cc]/denom[cc])*(erdenom[cc]/denom[cc]));
//            erratio[cc]=TMath::Sqrt(TMath::Abs(num[cc]*(1-num[cc]/denom[cc])))/denom[cc];
            if(erratio[cc] == 0.0) erratio[cc]=0.01;
        }
        else{
            ratio[cc]=0.0;
            erratio[cc]=0.01;
        }
        ratio456TOF->SetBinContent(cc+1,ratio[cc]);
        ratio456TOF->SetBinError(cc+1,erratio[cc]);
    }
    
    for(Int_t cc=0;cc<nPtBins;cc++){                // 2 SPD + any
        num[cc]=fHistPtITSMISPDInAccTOFMC->GetBinContent(cc);
        ernum[cc]=fHistPtITSMISPDInAccTOFMC->GetBinError(cc);
        denom[cc]=fHistPtITSMISPDInAccTOF->GetBinContent(cc);
        erdenom[cc]=fHistPtITSMISPDInAccTOF->GetBinError(cc);
        if(denom[cc] > 0.){
            ratio[cc]=num[cc]/denom[cc];
            erratio[cc]=num[cc]/denom[cc] * TMath::Sqrt((ernum[cc]/num[cc])*(ernum[cc]/num[cc])+(erdenom[cc]/denom[cc])*(erdenom[cc]/denom[cc]));
//            erratio[cc]=TMath::Sqrt(TMath::Abs(num[cc]*(1-num[cc]/denom[cc])))/denom[cc];
            if(erratio[cc] == 0.0) erratio[cc]=0.01;
        }
        else{
            ratio[cc]=0.0;
            erratio[cc]=0.01;
        }
        ratio2anyTOF->SetBinContent(cc+1,ratio[cc]);
        ratio2anyTOF->SetBinError(cc+1,erratio[cc]);
    }

    ctitle = GetRunNumber()+"ITS-TPC matchMC/matchDATA ratio";
    TCanvas* cMatchRatio = new TCanvas("cMatchRatio",ctitle,10,10,1000,1000);
    cMatchRatio->Divide(1,2);
    cMatchRatio->cd(1);
    ratio2any->SetLineColor(2);
    ratio2any->SetMaximum(2.0);
    ratio2any->SetMinimum(0);
    ratio2any->GetXaxis()->SetRangeUser(0.1,30);
    ratio2any->Draw();
    TLine *l1 = new TLine(0.,1.,30.,1.);
    l1->SetLineColor(kGreen+1);
    l1->Draw();
    cMatchRatio->cd(2);
    ratio456->SetLineColor(4);
    ratio456->SetMaximum(2.0);
    ratio456->SetMinimum(0);
    ratio456->GetXaxis()->SetRangeUser(0.1,30);
    ratio456->Draw();
    TLine *l1b = new TLine(0.,1.,30.,1.);
    l1b->SetLineColor(kOrange+1);
    l1b->Draw();
    cMatchRatio->Update();
    cMatchRatio->SaveAs("MatchRatio_TPCITS.pdf");
    pdfFileNames+=" MatchRatio_TPCITS.pdf";

    ctitle = GetRunNumber()+"TOF-TPC matchMC/matchDATA ratio";
    TCanvas* cMatchRatioTOF = new TCanvas("cMatchRatioTOF",ctitle,10,10,1000,1000);
    cMatchRatioTOF->Divide(1,2);
    cMatchRatioTOF->cd(1);
    ratio2anyTOF->SetLineColor(2);
    ratio2anyTOF->SetMaximum(2.0);
    ratio2anyTOF->SetMinimum(0);
    ratio2anyTOF->GetXaxis()->SetRangeUser(0.1,30);
    ratio2anyTOF->Draw();
    l1->Draw();
    cMatchRatioTOF->cd(2);
    ratio456TOF->SetLineColor(4);
    ratio456TOF->SetMaximum(2.0);
    ratio456TOF->SetMinimum(0);
    ratio456TOF->GetXaxis()->SetRangeUser(0.1,30);
    ratio456TOF->Draw();
    l1b->Draw();
    cMatchRatioTOF->Update();
    cMatchRatioTOF->SaveAs("MatchRatio_TOFITS.pdf");
    pdfFileNames+=" MatchRatio_TOFITS.pdf";

    return kTRUE;
}

//_______________________________________________________________________
void PlotSPDVtxPileup(TFile *fildat, TCanvas **&clist, Int_t &cnum){
    
    gROOT->SetStyle("Plain");
    gStyle->SetPalette(1);
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(111);
    cnum=2; // number of canvases
    clist= new TCanvas* [5];//array of pointers to TCanvases
    
    TDirectoryFile *pildata = (TDirectoryFile*)fildat->Get("CheckPileupQA");
    if(!pildata) {
        cout << "CheckPileupQA directory not found for DATA -- return" << endl;
        cnum =-1;
        return;
    }
    
    pildata->cd();
    TList *fListData = (TList*)pildata->Get("clistPileupSPDQA");
    if(!fListData) {
        cout << "clistPileupSPDQA TList not found for DATA -- return" << endl;
        cnum =-1;
        return;
    }
    
    TString fTitleData = "Data";
    TString fTitleMc = "MC";
    
//    TDirectoryFile *spdmc = (TDirectoryFile*)filMC->Get("CheckPileupQA");
//    spdmc->cd();
//    TList *fListMc = (TList*)spdmc->Get("clistPileupSPDQA");

    TH1F *PilVtxData = (TH1F*)fListData->FindObject("hNOfPileupVertSPD");
    PilVtxData->SetTitle(Form("%s %s",PilVtxData->GetTitle(),fTitleData.Data()));
    
    
    TString ctitle = GetRunNumber()+"pileup vertex";
    TCanvas *vertices = new TCanvas("pileupvtx",ctitle,1000,600);
    clist[0]=vertices;
    vertices->SetLogy();
    PilVtxData->Draw();
    PilVtxData->GetXaxis()->SetTitle("pileup vertices");
    TLatex* tp1=new TLatex(0.5,0.85,"pileup vertices multiplicity");
    tp1->SetNDC();
    tp1->Draw();
    vertices->SaveAs("Pileup_SPD_vtx.pdf");
    pdfFileNames+=" Pileup_SPD_vtx.pdf";

    TH1F *PilTrklData = (TH1F*)fListData->FindObject("hNtracklPilSPD");
    PilTrklData->SetTitle(Form("%s %s",PilTrklData->GetTitle(),fTitleData.Data()));
    Double_t npil=PilTrklData->GetEntries();

    TH1F *NoPilTrklData = (TH1F*)fListData->FindObject("hNtracklNoPilSPD");
    NoPilTrklData->SetTitle(Form("%s - %s",NoPilTrklData->GetTitle(),fTitleData.Data()));
    Double_t nnopil=NoPilTrklData->GetEntries();
    
    ctitle = GetRunNumber()+"tracklets";
    TCanvas *trackletspil = new TCanvas("trackletspil",ctitle,1000,600);
    clist[1]=trackletspil;
    trackletspil->SetLogy();
    PilTrklData->SetLineColor(2);
    NoPilTrklData->Draw();
    NoPilTrklData->GetXaxis()->SetTitle("tracklets");
    PilTrklData->Draw("same");
    TLatex* tphi=new TLatex(0.25,0.85, "Red = pileup; Black = no pileup");
    tphi->SetNDC();
    tphi->Draw();
    trackletspil->SaveAs("Pileup_Trkl.pdf");
    pdfFileNames+=" Pileup_Trkl.pdf";

    TH1F *PilCl1Data = (TH1F*)fListData->FindObject("hNCL1PilSPD");
    PilCl1Data->SetTitle(Form("%s - %s",PilCl1Data->GetTitle(),fTitleData.Data()));
    
    TH1F *NoPilCl1Data = (TH1F*)fListData->FindObject("hNCL1NoPilSPD");
    NoPilCl1Data->SetTitle(Form("%s - %s",NoPilCl1Data->GetTitle(),fTitleData.Data()));
    
    ctitle = GetRunNumber()+"clusters1";
    TCanvas *clusters1= new TCanvas("clusters1",ctitle,1000,600);
    clist[2]=clusters1;
    clusters1->SetLogy();
    PilCl1Data->SetLineColor(2);
    NoPilCl1Data->Draw();
    NoPilCl1Data->GetXaxis()->SetTitle("clusters");
    PilCl1Data->Draw("same");
    tphi->Draw();
    //    TLatex* tphi=new TLatex(0.25,0.85,Form("Red = %d DATA; Blue = %d MC",gRunNumber,gRunNumberMC));
    clusters1->SaveAs("Pileup_CL1.pdf");
    pdfFileNames+=" Pileup_CL1.pdf";

    TH1F *PilContrData = (TH1F*)fListData->FindObject("hContribPrimVertPilSPD");
    PilContrData->SetTitle("Number of contributors to primary vertex");
    
    TH1F *NoPilContrData = (TH1F*)fListData->FindObject("hContribPrimVertNoPilSPD");
    NoPilContrData->SetTitle(Form("%s - %s",NoPilContrData->GetTitle(),fTitleData.Data()));
    
    ctitle = GetRunNumber()+"contrib";
    TCanvas *contrib = new TCanvas("contrib",ctitle,1000,600);
    clist[3]=contrib;
    contrib->SetLogy();
    PilContrData->SetLineColor(2);
    NoPilContrData->Draw();
    NoPilContrData->GetXaxis()->SetTitle("SPD1 contributors");
    PilContrData->Draw("same");
    TLatex* tphi3=new TLatex(0.25,0.85,"Red = contr. prim. vertex pileup tagging");
    TLatex* tphi3s=new TLatex(0.25,0.78,"Black = contr. prim. vertex no pileup tagging");
    tphi3->SetNDC();
    tphi3s->SetNDC();
    tphi3->Draw(); tphi3s->Draw();
    clusters1->SaveAs("Pileup_contrib.pdf");
    pdfFileNames+=" Pileup_contrib.pdf";

    TH1F *PilZdiffPData = (TH1F*)fListData->FindObject("hZDiffFirstPilSPD");
    PilZdiffPData->SetTitle(Form("%s - %s",PilZdiffPData->GetTitle(),fTitleData.Data()));
    
    TH1F *PilZdiffSData = (TH1F*)fListData->FindObject("hZDiffSecondPilSPD");
    PilZdiffSData->SetTitle(Form("%s - %s",PilZdiffSData->GetTitle(),fTitleData.Data()));
    
    ctitle = GetRunNumber()+"zdiff";
    TCanvas *zdiff = new TCanvas("zdiff",ctitle,1000,600);
    clist[4]=zdiff;
    zdiff->SetLogy();
    PilZdiffPData->SetLineColor(2);
    PilZdiffPData->Draw();
    PilZdiffPData->GetXaxis()->SetTitle("zprim - zpileup (cm)");
    PilZdiffSData->Draw("same");
    TLatex* tz=new TLatex(0.15,0.85,"Red = first pileup vertex, Black = second pileup vertex");
    tz->SetNDC();
    tz->Draw();
    //    TLatex* tphi=new TLatex(0.25,0.85,Form("Red = %d DATA; Blue = %d MC",gRunNumber,gRunNumberMC));
    clusters1->SaveAs("Pileup_zdiff.pdf");
    pdfFileNames+=" Pileup_zdiff.pdf";
}



//_______________________________________________________________________

void PlotPID(TFile *fildat, TCanvas **&clist, Int_t &cnum){
    
    gROOT->SetStyle("Plain");
    gStyle->SetPalette(1);
//    gStyle->SetOptStat(0);
    gStyle->SetOptFit(111);
    cnum=2; // number of canvases
    clist= new TCanvas* [2];//array of pointers to TCanvases

    TString fTitleData = " - TPCITS tracks";
//    TString fTitleMc = "MC";
    
    TDirectoryFile *QAdata = (TDirectoryFile*)fildat->Get("PIDqa");
    if(!QAdata) {
        cout << "PIDqa directory not found for DATA -- return" << endl;
        cnum =-1;
        return;
    }
    TDirectoryFile *fListData = (TDirectoryFile*)QAdata->Get("PIDqa");
    if(!fListData) {
        cout << "PIDqa directory not found for DATA -- return" << endl;
        cnum =-1;
        return;
    }
   TList *fListITS = (TList*)fListData->FindObject("ITS");
    if(!fListITS) {
        cout << "ITS TList not found for DATA -- return" << endl;
        cnum =-1;
        return;
    }

    TH2F *nsigma_pi = (TH2F *)fListITS->FindObject("hNsigmaP_ITS_pion");
    nsigma_pi->SetTitle(Form("%s %s",nsigma_pi->GetTitle(),fTitleData.Data()));
    TH2F *nsigma_pro = (TH2F *)fListITS->FindObject("hNsigmaP_ITS_proton");
    nsigma_pro->SetTitle(Form("%s %s",nsigma_pro->GetTitle(),fTitleData.Data()));

    
    TString ctitle = GetRunNumber()+"ITS QA";
    TCanvas *its_qa = new TCanvas("its_qa",ctitle,1000,600);
    clist[0]=its_qa;
    its_qa->Divide(2,1);
    TVirtualPad *pad1 = its_qa->GetPad(1);
    pad1->SetLogz();
    its_qa->cd(1);
    nsigma_pi->GetXaxis()->SetRangeUser(0.0,3.0);
    nsigma_pi->Draw("colz");
    its_qa->cd(2);
//    nsigma_pro->GetXaxis()->SetRangeUser(0.0,3.0);
//    nsigma_pro->Draw("");
    TH1F *profX = (TH1F *)nsigma_pi->ProfileX();
    profX->Draw("");
    
    its_qa->SaveAs("ITS_pid_qa.pdf");
    pdfFileNames+=" ITS_pid_qa.pdf";

    TList *fListITS_pureSA = (TList*)fListData->FindObject("ITS_PureSA");
    fTitleData = " - ITSpureSA tracks";
    TH2F *nsigma_pi_pSA = (TH2F *)fListITS_pureSA->FindObject("hNsigmaP_ITSPureSA_pion");
    nsigma_pi_pSA->SetTitle(Form("%s %s",nsigma_pi_pSA->GetTitle(),fTitleData.Data()));
    TH2F *nsigma_pro_pSA = (TH2F *)fListITS_pureSA->FindObject("hNsigmaP_ITSPureSA_proton");
    nsigma_pro_pSA->SetTitle(Form("%s %s",nsigma_pro_pSA->GetTitle(),fTitleData.Data()));
    
    
    ctitle = GetRunNumber()+"ITSPureSA QA";
    TCanvas *itspSA_qa = new TCanvas("itspSA_qa",ctitle,1000,600);
    clist[1]=itspSA_qa;
    itspSA_qa->Divide(2,1);
    TVirtualPad *pad2 = itspSA_qa->GetPad(1);
    pad2->SetLogz();
    itspSA_qa->cd(1);
    nsigma_pi_pSA->GetXaxis()->SetRangeUser(0.0,3.0);
    nsigma_pi_pSA->Draw("colz");
    itspSA_qa->cd(2);
    nsigma_pro_pSA->GetXaxis()->SetRangeUser(0.0,3.0);
    nsigma_pro_pSA->Draw("");
    
    itspSA_qa->SaveAs("ITSpSA_pid_qa.pdf");
    pdfFileNames+=" ITSpSA_pid_qa.pdf";
    
}
    //_______________________________________________________________________

void SaveC(TFile &fout, TCanvas**& clist, Int_t cnum){
  TDirectory *current = gDirectory;
  fout.cd();
  for(Int_t i=0;i<cnum;i++)clist[i]->Write();
  delete[] clist;
  current->cd();
}

//_______________________________________________________________________
TString GetRunNumber(){
  // returns a string with the run number
  char rn[10];
  sprintf(rn,"%d  ",gRunNumber);
  TString str(rn);
  return str;
}

