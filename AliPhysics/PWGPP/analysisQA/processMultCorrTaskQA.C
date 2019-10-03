/*
//  processMultCorrTaskQA.C                                                  //
//  This macro process plots after the merge, corresponding to plots of      //
//  Multiplicity probabilitys and correlations for                           //
//  different V0 estimators.                                                 //
//___________________________________________________________________________//
 *  Author: Héctor Bello Martínez                                            *
 *  hector.bello.martinez@cern.ch                                            *
 *  for: ALICE collaboration                                                 *
//___Created by Hector on 14/06/15._________________________________________ */


#include <TFile.h>
#include <TList.h>
#include <TTree.h>
#include <TChain.h>
#include <TClonesArray.h>
#include <TObjString.h>
#include <TString.h>
#include <TROOT.h>
#include <TMath.h>
#include <TH2D.h>
#include <TProfile.h>
#include "TCanvas.h"
#include "TLegend.h"
#include <iostream>
#include <fstream>
#include <string>




void processMultCorrTaskQA(const char *filePath = "AnalysisResults.root",
                        TString suffix="eps",
                        const char *outfile="MultCorrtaskQA_output.root"){
    
    TFile *Resfile = TFile::Open(filePath, "read");
    if(!Resfile) {
        Printf("FATAL: File doesn't exist");
        return;
    }
    //TFile *Resfile = new TFile("AnalysisResults.root");
    
    TList* listPP =new TList();
    listPP= (TList*) Resfile->Get("PWGLF_PPVsMult/cList");
    if(!listPP) {
        Printf("FATAL: TList - PWGLF_PPVsMult/cList -  doesn't exist");
        return;
    }

    
    TH1F* Mult[5];
    TH1F* MultPerC[4];
    TH2D* MultCorr[5];
    TH2D* MultCorrPerC[4];
    TProfile* corrpfx[4];
    TProfile* corrpfxPerC[4];
    
    
    TH1F* EventSel = (TH1F*)listPP->FindObject("fHistEventCounter");
    TH1F* dNdeta = (TH1F*)listPP->FindObject("fdNdeta");
    TH2D* Modules = (TH2D*)listPP->FindObject("fModulesV0");
    
    Mult[0] = (TH1F*)listPP->FindObject("fHistRefMult08");
    Mult[1] = (TH1F*)listPP->FindObject("fHistRefMult05");
    Mult[2] = (TH1F*)listPP->FindObject("fHistV0Aamp");
    Mult[3] = (TH1F*)listPP->FindObject("fHistV0Camp");
    Mult[4] = (TH1F*)listPP->FindObject("fHistV0Mamp");
    MultPerC[1] = (TH1F*)listPP->FindObject("fHistV0A");
    MultPerC[2] = (TH1F*)listPP->FindObject("fHistV0C");
    MultPerC[3] = (TH1F*)listPP->FindObject("fHistV0M");
 
    TH2D* MultCorr[0]= (TH2D*)listPP->FindObject("fcorrRef05Ref08");
    TH2D* MultCorr[1]= (TH2D*)listPP->FindObject("fcorrV0AampRef08");
    TH2D* MultCorr[2]= (TH2D*)listPP->FindObject("fcorrV0CampRef08");
    TH2D* MultCorr[3]= (TH2D*)listPP->FindObject("fcorrV0MampRef08");
    TH2D* MultCorrPerC[1]= (TH2D*)listPP->FindObject("fcorrV0ARef08");
    TH2D* MultCorrPerC[2]= (TH2D*)listPP->FindObject("fcorrV0CRef08");
    TH2D* MultCorrPerC[3]= (TH2D*)listPP->FindObject("fcorrV0MRef08");

    TProfile* corrpfx[0]= (TProfile*)listPP->FindObject("fcorrRef05Ref08pfx");
    TProfile* corrpfx[1]= (TProfile*)listPP->FindObject("fcorrV0AampRef08pfx");
    TProfile* corrpfx[2]= (TProfile*)listPP->FindObject("fcorrV0CampRef08pfx");
    TProfile* corrpfx[3]= (TProfile*)listPP->FindObject("fcorrV0MampRef08pfx");
    TProfile* corrpfxPerC[1]= (TProfile*)listPP->FindObject("fcorrV0ARef08pfx");
    TProfile* corrpfxPerC[2]= (TProfile*)listPP->FindObject("fcorrV0CRef08pfx");
    TProfile* corrpfxPerC[3]= (TProfile*)listPP->FindObject("fcorrV0MRef08pfx");

    TProfile* HPcorrpfx[5];
    TProfile* HPcorrpfxPerC[5];
    for(Int_t i=0; i<4; i++){
        HPcorrpfx[i] = new TProfile(Form("HPcorrpfx%d",i),";;",700,0,700);
    }
    for(Int_t i=1; i<4; i++){
        HPcorrpfxPerC[i] = new TProfile(Form("HPcorrpfxPerC%d",i),";;",700,0,700);
    }
    //-------------------------------------------------------------
    Int_t etabins=100;
    //MBins{0,1, 2, 5, 8, 12, 17, 22, 27, 35, 45, 55, 65, 85, 100}
    //{[1, 3], [4, 6], [7, 9], [10, 14], [15, 19] [20, 24] [25, 29] [30, 39] [40, 49] [50, 59] [60, 69] [70, 99] [100, infty]
    // NBins[12]={0, 0.01, 0.1, 1, 5, 10, 15, 20, 30, 40, 50, 70,100};
    //Int_t V0Bins=12;
    //-------------------------------------------------------------
    TH1F* HistMult[5];
    TH1F* HistMultPerC[4];
    
    TH1D* hdeta= new TH1D("hdeta","(1/Nev)dN/d#eta distribution;#eta;(1/N_{ev})dN/d#eta",etabins,-1,1);
  for(Int_t i=0; i<5; i++){
     if(i==0)
    HistMult[i] = new TH1F( Form("HistMult_%d",i), "Multiplicity |#eta| < 0.8; Ref. Mult. |#eta| < 0.8; P(Nch)",200,0,200);
    else if(i==1){
    HistMult[i] = new TH1F( Form("HistMult_%d",i), "Multiplicity |#eta| < 0.5; Ref. Mult. |#eta| < 0.5; P(Nch)",200,0,200);
    //HistMultPerC[i] = new TH1F( Form("HistMultPerC_%d",i), "Multiplicity V0A; V0A Percentil; P(Nch)",V0Bins,NBins);}
    HistMultPerC[i] = new TH1F( Form("HistMultPerC_%d",i), "Multiplicity V0A; V0A Percentil; P(Nch)",100,0,100);}
    else if(i==2){
    HistMult[i] = new TH1F( Form("HistMult_%d",i), "Multiplicity V0A; V0A Amplitud; P(Nch)",700,0,700);
    HistMultPerC[2] = new TH1F( Form("HistMultPerC_%d",i), "Multiplicity V0C; V0C Percentil; P(Nch)",100,0,100);}
    else if(i==3){
    HistMult[i] = new TH1F( Form("HistMult_%d",i), "Multiplicity V0C; V0C Amplitud; P(Nch)",700,0,700);
    HistMultPerC[3] = new TH1F( Form("HistMultPerC_%d",i), "Multiplicity V0M; V0M Percentil; P(Nch)",100,0,100);}
    else if(i==4)
    HistMult[i] = new TH1F( Form("HistMult_%d",i), "Multiplicity V0M; V0M Amplitud; P(Nch)",700,0,700);
  }
    
    Double_t IntegMult[5]=0;
    Double_t IntegMultPerC[4]=0;
    Double_t contmult[5]=0;
    Double_t contmultPerC[4]=0;
    Double_t Errmult[5]=0;
    Double_t ErrmultPerC[4]=0;
    
    for(Int_t i=0; i<5; i++){
       IntegMult[i]=Mult[i]->Integral(0,700);
        for(Int_t in=1; in<700+1; in++){
            if(IntegMult[i]==0)continue;
            contmult[i]=Mult[i]->GetBinContent(in)/IntegMult[i];
            Errmult[i]=Mult[i]->GetBinError(in)/IntegMult[i];
            HistMult[i]->SetBinContent(in,contmult[i]);
            HistMult[i]->SetBinError(in,Errmult[i]);
        }
    }

    
   
    TCanvas* cCorrMult[5];
    for(Int_t i=0; i<4; i++){
        HPcorrpfx[i]=MultCorr[i]->ProfileX();
        cCorrMult[i] = new TCanvas(Form("cCorrMult_%d",i),Form("TH2CorrMult_%d",i),200,10,900,600);
        //cCorrMult[i]->SetLogz();
        gStyle->SetOptStat(1011);
        gStyle->SetTitleFontSize(0.09);
        MultCorr[i]->SetTitleSize(0.05,"X");
        MultCorr[i]->SetTitleSize(0.05,"Y");
        //HPcorrpfx[i]->SetLineColor("1");
        //TLegend* legstvspt = new TLegend(0.2, 0.1, .5, 0.4);
        //HistMult[i]->SetMarkerStyle(20);
        //HistMult[i]->SetMarkerColor(2);
        //legstvspt->AddEntry(stvspt[1], "Soft scatt. perugia2011", "p");
        MultCorr[i]->Draw("Colz");
        HPcorrpfx[i]->Draw("sames");
        //corrpfx[i]->Draw("sames");
        cCorrMult[i]->SaveAs(Form("fig_CorrMult_%d.%s",i,suffix.Data()));
    }
    
    TCanvas* cCorrMultPerC[5];
    for(Int_t i=1; i<4; i++){
        HPcorrpfxPerC[i]=MultCorrPerC[i]->ProfileX();
        cCorrMultPerC[i] = new TCanvas(Form("cCorrMultPerC_%d",i),Form("TH2CorrMultPerC_%d",i),200,10,900,600);
        //cCorrMultPerC[i]->SetLogz();
        gStyle->SetOptStat(1011);
        gStyle->SetTitleFontSize(0.09);
        MultCorrPerC[i]->SetTitleSize(0.05,"X");
        MultCorrPerC[i]->SetTitleSize(0.05,"Y");
        MultCorrPerC[i]->Draw("Colz");
        HPcorrpfxPerC[i]->Draw("sames");
        //corrpfxPerC[i]->Draw("sames");
        cCorrMultPerC[i]->SaveAs(Form("fig_CorrMultPerC_%d.%s",i,suffix.Data()));
    }

    
    TCanvas* cMult[4];
    for(Int_t i=0; i<5; i++){
    cMult[i] = new TCanvas(Form("cMult_%d",i),Form("THMult_%d",i),200,10,900,600);
    cMult[i]->SetLogy();
    gStyle->SetOptStat(1011);
    gStyle->SetTitleFontSize(0.09);
    HistMult[i]->SetTitleSize(0.06,"X");
    HistMult[i]->SetTitleSize(0.06,"Y");
    HistMult[i]->GetXaxis()->SetRangeUser(1,700);
    //TLegend* legstvspt = new TLegend(0.2, 0.1, .5, 0.4);
    //HistMult[i]->SetMarkerStyle(20);
    //HistMult[i]->SetMarkerColor(2);
    HistMult[i]->SetLineColor(2);
    //legstvspt->AddEntry(stvspt[1], "Soft scatt. perugia2011", "p");
    //HistMult[i]->Draw("E p");
    HistMult[i]->Draw("Hist");
    cMult[i]->SaveAs(Form("fig_Multiplicity_%d.%s",i,suffix.Data()));
    }
    //legstvspt->Draw();
 
    TCanvas* cMultPerC[4];
    for(Int_t i=1; i<4; i++){
        cMultPerC[i] = new TCanvas(Form("cMultPerC_%d",i),Form("THMultPerC_%d",i),200,10,900,600);
        gStyle->SetOptStat(1011);
        gStyle->SetTitleFontSize(0.09);
        MultPerC[i]->SetTitleSize(0.06,"X");
        MultPerC[i]->SetTitleSize(0.06,"Y");
        MultPerC[i]->GetYaxis()->SetTitle("Counts");
        //TLegend* legstvspt = new TLegend(0.2, 0.1, .5, 0.4);
        MultPerC[i]->SetMarkerStyle(20);
        MultPerC[i]->SetMarkerColor(8);
        MultPerC[i]->GetYaxis()->SetRangeUser(0,100);
         //MultPerC[i]->GetYaxis()->SetRangeUser(0.8,1.2);
        MultPerC[i]->Fit("pol0");
        //legstvspt->AddEntry(stvspt[1], "Soft scatt. perugia2011", "p");
        MultPerC[i]->Draw("E p");
        cMultPerC[i]->SaveAs(Form("fig_MultPerC_%d.%s",i,suffix.Data()));
    }

    
    Double_t IntegNev=dNdeta->Integral(-1,1);
    Double_t conteta=0;
    Double_t Erreta=0;
    for(Int_t in=0; in<etabins+1; in++){
        conteta=dNdeta->GetBinContent(in)/IntegNev;
        Erreta=dNdeta->GetBinError(in)/IntegNev;
        hdeta->SetBinContent(in,conteta);
        hdeta->SetBinError(in,Erreta);
    }
    
    TCanvas* cEvSel = new TCanvas("cEvSel","THEvSel",200,10,900,600);
    gStyle->SetOptStat(1011);
    gStyle->SetTitleFontSize(0.09);
    EventSel->SetTitleSize(0.06,"X");
    EventSel->SetTitleSize(0.06,"Y");
    //EventSel->SetFillColor(4);
    /*EventSel->GetYaxis()->SetTitle("counts");
    EventSel->GetXaxis()->SetTitle("Criteria");
    EventSel->GetXaxis()->SetRangeUser(0,7);
    EventSel->GetYaxis()->SetRangeUser(0.8,1.2);
    EventSel->SetTitle("Event Selection");
     */
    //TLegend* legstvspt = new TLegend(0.2, 0.1, .5, 0.4);
    //stvspt[1]->SetMarkerStyle(20);
    //stvspt[1]->SetMarkerColor(6);
    //legstvspt->AddEntry(stvspt[1], "Soft scatt. perugia2011", "p");
    EventSel->Draw();
    cEvSel->SaveAs(Form("fig_EvSel.%s",suffix.Data()));
    //legstvspt->Draw();
    
    TCanvas* ceta = new TCanvas("ceta","THETA",200,10,900,600);
    ceta->SetLogy();
    gStyle->SetOptStat(1011);
    gStyle->SetTitleFontSize(0.09);
    hdeta->SetTitleSize(0.06,"X");
    hdeta->SetTitleSize(0.06,"Y");
    hdeta->SetTitle("#eta distribution");
    //hdeta->GetYaxis()->SetTitle("(1/Nev)dN/d#eta");
    hdeta->GetYaxis()->SetRangeUser(0.001,1);
    hdeta->GetXaxis()->SetTitle("#eta (rads)");
    hdeta->Draw("Ep");
    ceta->SaveAs(Form("fig_eta.%s",suffix.Data()));
    
    TCanvas* cmodul = new TCanvas("cmodul","THmodul",200,10,900,600);
    cmodul->SetLogz();
    gStyle->SetOptStat(1011);
    gStyle->SetTitleFontSize(0.09);
    Modules->SetTitleSize(0.05,"X");
    Modules->SetTitleSize(0.05,"Y");
    Modules->GetYaxis()->SetTitle("Multiplicyty");
    Modules->GetXaxis()->SetTitle("V0 Modules");
    Modules->GetYaxis()->SetRangeUser(0,200);
    Modules->SetTitle("Multiplicity in V0 Modules");
    Modules->Draw("Colz");
    cmodul->SaveAs(Form("fig_V0Modules.%s",suffix.Data()));
   // TFile *fout=new TFile("ProcesedMultCorrTask.root","recreate");
    TFile *fout = TFile::Open(outfile,"UPDATE");
    fout->ls();
    
    TDirectoryFile *cdd = NULL;
    cdd = (TDirectoryFile*)fout->Get("RESULTS");
    if(!cdd) {
        Printf("Warning: RESULTS <dir> doesn't exist, creating a new one");
        cdd = (TDirectoryFile*)fout->mkdir("RESULTS");
    }
    cdd->cd();
    cdd->ls();
    //fout->cd();
                                           
    EventSel->Write();
    Modules->Write();
    hdeta->Write();
    for(Int_t i=1; i<4; i++){
        MultPerC[i]->Write();
        MultCorrPerC[i]->Write();
    }
    for(Int_t i=0; i<5; i++){
        HistMult[i]->Write();
    }
    for(Int_t i=0; i<4; i++){
        MultCorr[i]->Write();
    }
    //for(Int_t i=0; i<3; i++){
    
    //}
 
    fout->Close();
    


}

