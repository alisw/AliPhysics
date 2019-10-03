/*


gSystem->SetIncludePath("-I$ROOTSYS/include -I$ALICE_ROOT/ -I$ALICE_ROOT/include");
gSystem->Load("libANALYSIS");
gSystem->Load("libANALYSISalice");
.L MakeSensitivityPlots.C++g 
 .L MatchingEffPlots.C++g
MatchingEffPlots()


 */

#include "TH1D.h"
#include "TCanvas.h"
#include "TFile.h"
#include "THn.h"
#include "TList.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TLatex.h"
#include "MakeSensitivityPlots.C"

const int pid = 1;      // pid = 1 —> pions, pid = 2 —> kaons, pid = 3 —> protons, pid = 5 —> unidentified

//const char

void MatchingEffPlots(){
    
  //gROOT->LoadMacro("ll.C");
  //    gROOT->LoadMacro("MakeSensitivityPlots2.C++");
  
  // ll();
    
    TCanvas *c[3];
    TCanvas *c2[3];
    const Char_t *proj[3] = {"pT","eta","phi"};
    
    TH1D *data[4];
    TH1D *mc[4];
    
    TH1D *ratio[4];
    
    char dataset[5] = "bcde";

    for(Int_t i=0;i<3;i++){     //loop over the projections: pT, eta, phi
    
        c[i] = new TCanvas(Form("%s",proj[i]),Form("%s",proj[i]),1000,1000);
        gStyle->SetOptStat(0);
        c[i]->Divide(2,2);
    
        for(Int_t k=0;k<4;k++){     //loop over the periods: LHC10b,10c,10d,10e
            c[i]->cd(k+1);
            
            TString dataFile = Form("data_LHC10%c.root",dataset[k]);
            TString mcFile = Form("mc_LHC10%c.root",dataset[k]);
            
            data[k] = (TH1D*)GetITSTPCMatchingEff(dataFile/*dataFile.Data()*/,pid,i+1);
            mc[k]   = (TH1D*)GetITSTPCMatchingEff(mcFile/*mcFile.Data()*/,pid,i+1);
            data[k]->Draw("Ep");
            data[k]->SetTitle(Form("LHC10%c",dataset[k]));
            if (i==0) data[k]->GetXaxis()->SetTitle("#it{p}_{T}");
            if (i==1) data[k]->GetXaxis()->SetTitle("#eta");
            if (i==2) data[k]->GetXaxis()->SetTitle("#varphi");
        
            data[k]->SetLineColor(kRed-3);
            mc[k]->Draw("Ep,same");
            mc[k]->SetLineColor(kBlue-3);
            data[k]->SetMaximum(1);
        
            if(i==0){
                gPad->SetLogx();
            }
        
            if(i==1||i==0){
                data[k]->SetMinimum(.8);
            }
        
            TLegend *leg = new TLegend(0.15,0.15,.35,0.25);
            leg->AddEntry(data[k],"Data","f");
            leg->AddEntry(mc[k],"MC","f");
            leg->SetBorderSize(0);
            leg->SetFillColor(0);
            leg->Draw();
        }
        
        c2[i] = new TCanvas(Form("ratio%s",proj[i]),Form("ratio%s",proj[i]),1000,1000);
        gStyle->SetOptStat(0);
        c2[i]->Divide(2,2);
        
        for(Int_t k=0;k<4;k++){     //loop over the periods: LHC10b,10c,10d,10e
            c2[i]->cd(k+1);
            
            TString dataFile = Form("data_LHC10%c.root",dataset[k]);
            TString mcFile = Form("mc_LHC10%c.root",dataset[k]);
            
            data[k] = (TH1D*)GetITSTPCMatchingEff(dataFile.Data(),pid,i+1);
            mc[k]   = (TH1D*)GetITSTPCMatchingEff(mcFile.Data(),pid,i+1);
            ratio[k] = (TH1D*)data[k]->Clone();
            ratio[k]->Reset();
            ratio[k]->Divide(data[k],mc[k],1,1);
            ratio[k]->Draw("");
            ratio[k]->SetTitle(Form("LHC10%c",dataset[k]));
            if (i==0) ratio[k]->GetXaxis()->SetTitle("#it{p}_{T}");
            if (i==1) ratio[k]->GetXaxis()->SetTitle("#eta");
            if (i==2) ratio[k]->GetXaxis()->SetTitle("#varphi");
            
            ratio[k]->SetLineColor(kBlack);
            //ratio[k]->SetMaximum(1);
            
            if(i==0){
                gPad->SetLogx();
            }
            
            if(i==1||i==0){
                //data[k]->SetMinimum(.8);
            }
            
            TLegend *leg = new TLegend(0.15,0.15,.35,0.25);
            leg->AddEntry(ratio[k],"Data/MC","f");
            //leg->AddEntry(mc[k],"MC","f");
            leg->SetBorderSize(0);
            leg->SetFillColor(0);
            leg->Draw();
        }
    }

        
}
        
    /*    c[i]->cd(2);
        gStyle->SetOptStat(0);

        TH1D * data_10c = GetITSTPCMatchingEff("data_LHC10c.root",pid,i+1);
        TH1D * mc_10c   = GetITSTPCMatchingEff("mc_LHC10c.root",pid,i+1);
        data_10c->Draw("");
        data_10c->SetTitle("LHC10c");
        if (i==0) data_10c->GetXaxis()->SetTitle("#it{p}_{T}");
        if (i==1) data_10c->GetXaxis()->SetTitle("#eta");
        if (i==2) data_10c->GetXaxis()->SetTitle("#varphi");
        //data_10c->SetMarkerStyle(20);
        data_10c->SetLineColor(kRed-3);
        mc_10c->Draw("same");
        //mc_10c->SetMarkerStyle(20);
        mc_10c->SetLineColor(kBlue-3);
        data_10c->SetMaximum(1);
        
        if(i==0){
            gPad->SetLogx();
        }
        
        if(i==1||i==0){
            data_10c->SetMinimum(.8);
            
        }
        
        TLegend *leg = new TLegend(0.15,0.2,.35,0.3);
        leg->AddEntry(data_10c,"Data","f");
        leg->AddEntry(mc_10c,"MC","f");
        leg->SetBorderSize(0);
        gStyle->SetFillColor(0);
        leg->Draw();
        
 
        c[i]->cd(3);
        gStyle->SetOptStat(0);

        TH1D * data_10d = GetITSTPCMatchingEff("data_LHC10d.root",pid,i+1);
        TH1D * mc_10d   = GetITSTPCMatchingEff("mc_LHC10d.root",pid,i+1);
        data_10d->Draw("");
        data_10d->SetTitle("LHC10d");
        if (i==0) data_10d->GetXaxis()->SetTitle("#it{p}_{T}");
        if (i==1) data_10d->GetXaxis()->SetTitle("#eta");
        if (i==2) data_10d->GetXaxis()->SetTitle("#varphi");
        //data_10d->SetMarkerStyle(20);
        data_10d->SetLineColor(kRed-3);
        mc_10d->Draw("same");
        //mc_10d->SetMarkerStyle(20);
        mc_10d->SetLineColor(kBlue-3);
        data_10d->SetMaximum(1);
        
        if(i==0){
            gPad->SetLogx();
        }
        
        if(i==1||i==0){
            data_10d->SetMinimum(.8);
            
        }
        
        TLegend *leg = new TLegend(0.15,0.2,.35,0.3);
        leg->AddEntry(data_10d,"Data","f");
        leg->AddEntry(mc_10d,"MC","f");
        leg->SetBorderSize(0);
        gStyle->SetFillColor(0);
        leg->Draw();
        

        c[i]->cd(4);
        gStyle->SetOptStat(0);

        TH1D * data_10e = GetITSTPCMatchingEff("data_LHC10e.root",pid,i+1);
        TH1D * mc_10e   = GetITSTPCMatchingEff("mc_LHC10e.root",pid,i+1);
        data_10e->Draw("");
        data_10e->SetTitle("LHC10e");
        if (i==0) data_10e->GetXaxis()->SetTitle("#it{p}_{T}");
        if (i==1) data_10e->GetXaxis()->SetTitle("#eta");
        if (i==2) data_10e->GetXaxis()->SetTitle("#varphi");
        //data_10e->SetMarkerStyle(20);
        data_10e->SetLineColor(kRed-3);
        mc_10e->Draw("same");
        //mc_10e->SetMarkerStyle(20);
        mc_10e->SetLineColor(kBlue-3);
        data_10e->SetMaximum(1);
        
        if(i==0){
            gPad->SetLogx();
        }
        
        if(i==1 || i==0){
            data_10e->SetMinimum(.8);
            
        }
        
        TLegend *leg = new TLegend(0.15,0.2,.35,0.3);
        leg->AddEntry(data_10e,"Data","f");
        leg->AddEntry(mc_10e,"MC","f");
        leg->SetBorderSize(0);
        gStyle->SetFillColor(0);
        leg->Draw(); */
        

  /*  }
    
    for(Int_t i=0;i<3;i++){     //loop over the projections: pT, eta, phi

        c2[i] = new TCanvas(Form("ratio%s",proj[i]),Form("ratio%s",proj[i]),1000,1000);
        gStyle->SetOptStat(0);
        c2[i]->Divide(2,2);

    
        for(Int_t k=0;k<4;k++){     //loop over the periods: LHC10b,10c,10d,10e
            c2[i]->cd(k+1);
        
            TString dataFile = Form("data_LHC10%c.root",dataset[k]);
            TString mcFile = Form("mc_LHC10%c.root",dataset[k]);
        
            data[k] = (TH1D*)GetITSTPCMatchingEff(dataFile.Data(),pid,i+1);
            mc[k]   = (TH1D*)GetITSTPCMatchingEff(mcFile.Data(),pid,i+1);
            ratio[k] = (TH1D*)data[k]->Clone();
            ratio[k]->Reset();
            ratio[k]->Divide(data[k],mc[k],1,1);
            ratio[k]->Draw("");
            ratio[k]->SetTitle(Form("LHC10%c",dataset[k]));
            if (i==0) ratio[k]->GetXaxis()->SetTitle("#it{p}_{T}");
            if (i==1) ratio[k]->GetXaxis()->SetTitle("#eta");
            if (i==2) ratio[k]->GetXaxis()->SetTitle("#varphi");
        
            ratio[k]->SetLineColor(kBlack);
            //ratio[k]->SetMaximum(1);
        
            if(i==0){
                gPad->SetLogx();
            }
        
            if(i==1||i==0){
                //data[k]->SetMinimum(.8);
            }
        
            TLegend *leg = new TLegend(0.15,0.15,.35,0.25);
            leg->AddEntry(ratio[k],"Data/MC","f");
            //leg->AddEntry(mc[k],"MC","f");
            leg->SetBorderSize(0);
            leg->SetFillColor(0);
            leg->Draw();
        }

    
    //for identified: only vs pT
    
    }
}*/
