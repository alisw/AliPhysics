#include <TTree.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TFile.h>
#include <TString.h>
#include <TMath.h>
#include <TRandom.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TF3.h>
#include <TDatabasePDG.h>
#include <TParticlePDG.h>
#include <TF1.h>
#include <TStyle.h>
#include <TRandom.h>
#include <TBenchmark.h>
#include <TCanvas.h>
#include <TPaveText.h>
#include <TMultiGraph.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TGraphErrors.h>
#include <TLegendEntry.h>

#include "hCorrF.h"

Double_t fithist_flat(Double_t *x, Double_t *par);
Double_t nScale(TH1D *num, TH1D *den, Int_t b1, Int_t b2);


void mSetCanvas(TCanvas *cl) {
    cl->SetBorderMode(0);
    cl->SetBorderSize(2);
    cl->SetLeftMargin(0.2);
    cl->SetBottomMargin(0.2);
    cl->SetTopMargin(1.25);
    cl->SetRightMargin(0.15);
    cl->SetFrameBorderMode(0);
    cl->SetFrameBorderMode(0);
}

//============================================================================

void CorrF() {
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(1);
    gStyle->SetStatH(0.23);
    gStyle->SetLabelSize(0.08);
    gStyle->SetTitleSize(0.08);
    gStyle->SetTitleSize(0.08,"t");

    sum2cfs();
}
//============================================================================
void sum2cfs() {
    TCanvas *rawcfs[cfiles+5]; ///zmienic
    //ALICE data
    Int_t ic=0;
    TH1D *hc[cKt][cMu],*hn[cKt][cMu],*hd[cKt][cMu];
    TH1D *h1,*h2,*h3,*h4;
    const TString nfile[10]={"3/BC/final13bc_3","3/DE/final13de_3","1/BC/final13bc_1","1/DE/final13de_1" };
    
    for (Int_t i=0;i<cfiles;i++){
        rawcfs[i]=new TCanvas("rawcfs","rawcfs",1150,1650);
        rawcfs[i]->Divide(cKt,cMu);
        mSetCanvas(rawcfs[i]);
        TString filetemp= nfile[i]+".root";
        fexp=TFile::Open(filetemp);
        Int_t ic=0;
        for(Int_t m=0;m<cMu;m++){
            for(Int_t k=0;k<cKt;k++){
                ic++;
                rawcfs[i]->cd(ic);
                gPad->SetBottomMargin(0.2);
                gPad->SetLeftMargin(0.2);
                h1 = (TH1D*)fexp->Get(Form("NumcqinvKptpcM%ikT%i", m, k));
                h2 = (TH1D*)fexp->Get(Form("DencqinvKptpcM%ikT%i", m, k));
                
                hn[k][m]=new TH1D(*h1);
                hd[k][m]=new TH1D(*h2);
                hc[k][m]=new TH1D(*hn[k][m]);
                hc[k][m]->Divide(hd[k][m]);
                hc[k][m]->Scale(nScale(hn[k][m],hd[k][m],cbin1,cbin2));
                
                hcf0[i][k][m]=new TH1D(*hc[k][m]);

                hcf0[i][k][m]->GetXaxis()->SetRangeUser(0,0.5);
                hcf0[i][k][m]->SetMaximum(2.0);
                hcf0[i][k][m]->SetMinimum(0.5);
                hcf0[i][k][m]->SetTitle(cTitkT[k]+cTitMu[m]);
                hcf0[i][k][m]->SetMarkerStyle(gma[m]);
                hcf0[i][k][m]->SetMarkerColor(gcv[m]);
                hcf0[i][k][m]->SetLineColor(gcv[m]);
                
                
                hcf0[i][k][m]->GetXaxis()->SetTitle("q_{inv} (GeV/c)");
                hcf0[i][k][m]->GetXaxis()->SetLabelSize(0.07);
                hcf0[i][k][m]->GetXaxis()->SetTitleSize(0.08);
                hcf0[i][k][m]->GetXaxis()->SetTitleOffset(1.05);
                hcf0[i][k][m]->GetXaxis()->SetTickLength(0.065);
                hcf0[i][k][m]->GetXaxis()->SetNdivisions(505);
                hcf0[i][k][m]->GetYaxis()->SetTitle("CF");
                hcf0[i][k][m]->GetYaxis()->SetLabelSize(0.07);
                hcf0[i][k][m]->GetYaxis()->SetLabelOffset(0.01);
                hcf0[i][k][m]->GetYaxis()->SetTitleSize(0.08);
                hcf0[i][k][m]->GetYaxis()->SetTitleOffset(1.05);
                hcf0[i][k][m]->GetYaxis()->SetTickLength(0.025);
                hcf0[i][k][m]->GetYaxis()->SetNdivisions(505);


                hcf0[i][k][m]->Draw();
            }
        }
        filetemp=nfile[i]+".gif";
        filetemp.Remove(0,5);
        rawcfs[i]->SetTitle(filetemp);
    //    rawcfs[i]->Print(filetemp);
    }

    ic=0;
    rawcfs[cfiles+1]=new TCanvas("2","2",1150,1650);
    rawcfs[cfiles+1]->Divide(cKt,cMu);
    mSetCanvas(rawcfs[cfiles+1]);
     for(Int_t m=0;m<cMu;m++){
        for(Int_t k=0;k<cKt;k++){
            ic++;

            rawcfs[cfiles+1]->cd(ic);

           
            gPad->SetBottomMargin(0.2);
            gPad->SetLeftMargin(0.2);
            hc[k][m]= new TH1D(*hcf0[1][k][m]);
            hc[k][m]->Divide(hcf0[3][k][m]);

            hc[k][m]->Draw();
        }
    }
    rawcfs[cfiles+1]->Print("de_3_de_1.gif");

}
//============================================================================
Double_t nScale(TH1D *num, TH1D *den, Int_t b1, Int_t b2) {
  Double_t suma=num->Integral(b1,b2);
  Double_t sumb=den->Integral(b1,b2);
  if (sumb == 0) {
    cout << "Sorry Integral("<< den->GetName()<<")=0" <<endl;
    exit(1);
  }         
 return sumb/suma;
}
