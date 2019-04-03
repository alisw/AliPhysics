

//Macro to create tracking 3D (pt, eta, Zvtx) efficiency file for Dh Corr analysis
//Jitendra

#include <iostream>
#include "TSystem.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TDirectoryFile.h"
#include "THnSparse.h"
#include "TAxis.h"
#include "TH1D.h"
#include "TH3F.h"
#include "TH2F.h"

#if !defined(__CINT__) || defined(__MAKECINT__)
#include "AliCFContainer.h"
#include "AliCFGridSparse.h"
#endif

using namespace std;


void CheckFileQuality();
void DrawDhCorr_STE_3DEfficiency_in_pTEtaZvtx(){
    
    TFile* f = new TFile("AnalysisResults_pp2017_TrackEff.root");
    TDirectoryFile* d = (TDirectoryFile*)f->Get("PWGPP_CFSingleTrack");
    AliCFContainer *data = (AliCFContainer*) (d->Get("containerpp2017_pi"));
    
    AliCFGridSparse* gridSparsenum = (AliCFGridSparse*)data->GetGrid(1); // GenAcc
    THnSparse* numData = (THnSparse*)gridSparsenum->GetGrid();
    AliCFGridSparse* gridSparseden = (AliCFGridSparse*)data->GetGrid(6); // Reco
    THnSparse* denData = (THnSparse*)gridSparseden->GetGrid();


    TFile* f2 = new TFile("AnalysisResults_pp2017_TrackEff.root");
    TDirectoryFile* d2 = (TDirectoryFile*)f2->Get("PWGPP_CFSingleTrack");
    AliCFContainer *data2 = (AliCFContainer*) (d2->Get("containerpp2017_k"));
    
    AliCFGridSparse* gridSparsenum2 = (AliCFGridSparse*)data2->GetGrid(1); // GenAcc
    THnSparse* numData2 = (THnSparse*)gridSparsenum2->GetGrid();
    AliCFGridSparse* gridSparseden2 = (AliCFGridSparse*)data2->GetGrid(6); // Reco
    THnSparse* denData2 = (THnSparse*)gridSparseden2->GetGrid();
    
    
    TFile* f3 = new TFile("AnalysisResults_pp2017_TrackEff.root");
    TDirectoryFile* d3 = (TDirectoryFile*)f3->Get("PWGPP_CFSingleTrack");
    AliCFContainer *data3 = (AliCFContainer*) (d3->Get("containerpp2017_p"));
    
    AliCFGridSparse* gridSparsenum3 = (AliCFGridSparse*)data3->GetGrid(1); // GenAcc
    THnSparse* numData3 = (THnSparse*)gridSparsenum3->GetGrid();
    AliCFGridSparse* gridSparseden3 = (AliCFGridSparse*)data3->GetGrid(6); // Reco
    THnSparse* denData3 = (THnSparse*)gridSparseden3->GetGrid();
    
    
    TFile* f4 = new TFile("AnalysisResults_pp2017_TrackEff.root");
    TDirectoryFile* d4 = (TDirectoryFile*)f4->Get("PWGPP_CFSingleTrack");
    AliCFContainer *data4 = (AliCFContainer*) (d4->Get("containerpp2017_e"));
    
    AliCFGridSparse* gridSparsenum4 = (AliCFGridSparse*)data4->GetGrid(1); // GenAcc
    THnSparse* numData4 = (THnSparse*)gridSparsenum4->GetGrid();
    AliCFGridSparse* gridSparseden4 = (AliCFGridSparse*)data4->GetGrid(6); // Reco
    THnSparse* denData4 = (THnSparse*)gridSparseden4->GetGrid();
    
    
    TFile* f5 = new TFile("AnalysisResults_pp2017_TrackEff.root");
    TDirectoryFile* d5 = (TDirectoryFile*)f5->Get("PWGPP_CFSingleTrack");
    AliCFContainer *data5 = (AliCFContainer*) (d5->Get("containerpp2017_mu"));
    
    AliCFGridSparse* gridSparsenum5 = (AliCFGridSparse*)data5->GetGrid(1); // GenAcc
    THnSparse* numData5 = (THnSparse*)gridSparsenum5->GetGrid();
    AliCFGridSparse* gridSparseden5 = (AliCFGridSparse*)data5->GetGrid(6); // Reco
    THnSparse* denData5 = (THnSparse*)gridSparseden5->GetGrid();
    
           
    TH3D* hnum = (TH3D*)numData->Projection(0,1,4);
    TH3D* hden = (TH3D*)denData->Projection(0,1,4);
    
    TH3D* hnum2 = (TH3D*)numData2->Projection(0,1,4);
    TH3D* hden2 = (TH3D*)denData2->Projection(0,1,4);
    
    TH3D* hnum3 = (TH3D*)numData3->Projection(0,1,4);
    TH3D* hden3 = (TH3D*)denData3->Projection(0,1,4);
    
    TH3D* hnum4 = (TH3D*)numData4->Projection(0,1,4);
    TH3D* hden4 = (TH3D*)denData4->Projection(0,1,4);
    
    TH3D* hnum5 = (TH3D*)numData5->Projection(0,1,4);
    TH3D* hden5 = (TH3D*)denData5->Projection(0,1,4);                

    
    TCanvas *c = new TCanvas("c", "pT, Eta, Zvtx Efficiency distrubution", 400,900);
    c->Divide(1,2);

    TH3D* hnumTot = (TH3D*)hnum->Clone("hnumtot");
    hnumTot->Add(hnum2);
    hnumTot->Add(hnum3);
    hnumTot->Add(hnum4);
    hnumTot->Add(hnum5);
    TH3D* hdenTot = (TH3D*)hden->Clone("hdentot");
    hdenTot->Add(hden2);
    hdenTot->Add(hden3);
    hdenTot->Add(hden4);
    hdenTot->Add(hden5);
    
    TH3D* heff = (TH3D*)hdenTot->Clone("heff");
    heff->Divide(hdenTot,hnumTot,1,1,"B");
    
    c->cd(1);
    heff->Draw();
    //Rebinning
    const Int_t newLimitsBins = 26;
    Double_t* newLimits = new Double_t[newLimitsBins+1];
    newLimits[0]   = 0.30;
    newLimits[1]   = 0.40;
    newLimits[2]   = 0.50;
    newLimits[3]   = 0.60;
    newLimits[4]   = 0.70;
    newLimits[5]   = 0.80;
    newLimits[6]   = 0.90;
    newLimits[7]   = 1.00;
    newLimits[8]   = 1.20;
    newLimits[9]   = 1.40;
    newLimits[10]  = 1.60;
    newLimits[11]  = 1.80;
    newLimits[12]  = 2.00;
    newLimits[13]  = 2.33;
    newLimits[14]  = 2.66;
    newLimits[15]  = 3.00;
    newLimits[16]  = 3.50;
    newLimits[17]  = 4.00;
    newLimits[18]  = 4.40;
    newLimits[19]  = 4.80;
    newLimits[20]  = 5.20;
    newLimits[21]  = 6.00;
    newLimits[22]  = 7.00;
    newLimits[23]  = 8.00;
    newLimits[24]  = 10.00;
    newLimits[25]  = 16.00;
    newLimits[26]  = 24.00;
   
    THnSparse* RbnumData = (THnSparse*)numData->Clone("numNew");
    RbnumData->Reset();
    TAxis* axis = (TAxis*)RbnumData->GetAxis(0);
    axis->Set(newLimitsBins,newLimits);
    RbnumData->SetBinEdges(0,newLimits);
    RbnumData->RebinnedAdd(numData, 1);
    TH3D* Rbhnum = (TH3D*)RbnumData->Projection(0,1,4);
    gStyle->SetOptStat(kTRUE);
    
    THnSparse* RbdenData = (THnSparse*)denData->Clone("denNew");
    RbdenData->Reset();
    TAxis* axis2 = (TAxis*)RbdenData->GetAxis(0);
    axis2->Set(newLimitsBins,newLimits);
    RbdenData->SetBinEdges(0,newLimits);
    RbdenData->RebinnedAdd(denData, 1);
    TH3D* Rbhden = (TH3D*)RbdenData->Projection(0,1,4);
    gStyle->SetOptStat(kTRUE);
    
    THnSparse* RbnumData2 = (THnSparse*)numData2->Clone("numNew2");
    RbnumData2->Reset();
    TAxis* axis = (TAxis*)RbnumData2->GetAxis(0);
    axis->Set(newLimitsBins,newLimits);
    RbnumData2->SetBinEdges(0,newLimits);
    RbnumData2->RebinnedAdd(numData2, 1);
    TH3D* Rbhnum2 = (TH3D*)RbnumData2->Projection(0,1,4);
    gStyle->SetOptStat(kTRUE);
    
    THnSparse* RbdenData2 = (THnSparse*)denData2->Clone("denNew2");
    RbdenData2->Reset();
    TAxis* axis2 = (TAxis*)RbdenData2->GetAxis(0);
    axis2->Set(newLimitsBins,newLimits);
    RbdenData2->SetBinEdges(0,newLimits);
    RbdenData2->RebinnedAdd(denData2, 1);
    TH3D* Rbhden2 = (TH3D*)RbdenData2->Projection(0,1,4);
    gStyle->SetOptStat(kTRUE);
    
    THnSparse* RbnumData3 = (THnSparse*)numData3->Clone("numNew3");
    RbnumData3->Reset();
    TAxis* axis = (TAxis*)RbnumData3->GetAxis(0);
    axis->Set(newLimitsBins,newLimits);
    RbnumData3->SetBinEdges(0,newLimits);
    RbnumData3->RebinnedAdd(numData3, 1);
    TH3D* Rbhnum3 = (TH3D*)RbnumData3->Projection(0,1,4);
    gStyle->SetOptStat(kTRUE);
    
    THnSparse* RbdenData3 = (THnSparse*)denData3->Clone("denNew3");
    RbdenData3->Reset();
    TAxis* axis2 = (TAxis*)RbdenData3->GetAxis(0);
    axis2->Set(newLimitsBins,newLimits);
    RbdenData3->SetBinEdges(0,newLimits);
    RbdenData3->RebinnedAdd(denData3, 1);
    TH3D* Rbhden3 = (TH3D*)RbdenData3->Projection(0,1,4);
    gStyle->SetOptStat(kTRUE);
    
    THnSparse* RbnumData4 = (THnSparse*)numData4->Clone("numNew4");
    RbnumData4->Reset();
    TAxis* axis = (TAxis*)RbnumData4->GetAxis(0);
    axis->Set(newLimitsBins,newLimits);
    RbnumData4->SetBinEdges(0,newLimits);
    RbnumData4->RebinnedAdd(numData4, 1);
    TH3D* Rbhnum4 = (TH3D*)RbnumData4->Projection(0,1,4);
    gStyle->SetOptStat(kTRUE);
    
    THnSparse* RbdenData4 = (THnSparse*)denData4->Clone("denNew4");
    RbdenData4->Reset();
    TAxis* axis2 = (TAxis*)RbdenData4->GetAxis(0);
    axis2->Set(newLimitsBins,newLimits);
    RbdenData4->SetBinEdges(0,newLimits);
    RbdenData4->RebinnedAdd(denData4, 1);
    TH3D* Rbhden4 = (TH3D*)RbdenData4->Projection(0,1,4);
    gStyle->SetOptStat(kTRUE);
    
    THnSparse* RbnumData5 = (THnSparse*)numData5->Clone("numNew5");
    RbnumData5->Reset();
    TAxis* axis = (TAxis*)RbnumData5->GetAxis(0);
    axis->Set(newLimitsBins,newLimits);
    RbnumData5->SetBinEdges(0,newLimits);
    RbnumData5->RebinnedAdd(numData5, 1);
    TH3D* Rbhnum5 = (TH3D*)RbnumData5->Projection(0,1,4);
    gStyle->SetOptStat(kTRUE);
    
    THnSparse* RbdenData5 = (THnSparse*)denData5->Clone("denNew5");
    RbdenData5->Reset();
    TAxis* axis2 = (TAxis*)RbdenData5->GetAxis(0);
    axis2->Set(newLimitsBins,newLimits);
    RbdenData5->SetBinEdges(0,newLimits);
    RbdenData5->RebinnedAdd(denData5, 1);
    TH3D* Rbhden5 = (TH3D*)RbdenData5->Projection(0,1,4);
    gStyle->SetOptStat(kTRUE);
    
    TH3D* RbhdenTot = (TH3D*)Rbhden->Clone("Rbhdentot");
    RbhdenTot->Add(Rbhden2);
    RbhdenTot->Add(Rbhden3);
    RbhdenTot->Add(Rbhden4);
    RbhdenTot->Add(Rbhden5);
    TH3D* RbhnumTot = (TH3D*)Rbhnum->Clone("Rbhnumtot");
    RbhnumTot->Add(Rbhnum2);
    RbhnumTot->Add(Rbhnum3);
    RbhnumTot->Add(Rbhnum4);
    RbhnumTot->Add(Rbhnum5);
       
    TH3D* Rbheff = (TH3D*)RbhdenTot->Clone("heff_rebin");
    Rbheff->Divide(RbhdenTot,RbhnumTot,1,1,"B");
    gStyle->SetOptStat(kTRUE);
    c->cd(2);
    Rbheff->Draw();
    
    
    // Num bins checking
    Bool_t Rebintest = kTRUE;
    Float_t sum = 0, sumnew = 0;
    Float_t content = 0, content2 = 0;
    if(Rebintest){
        Printf("Original THnSparse <-------------");
        for (Int_t ibin = 1; ibin<=hnum->GetNbinsX(); ibin++){
            content=0;
            for (Int_t ibin2 = 1; ibin2<=hnum->GetNbinsY(); ibin2++){
                for (Int_t ibin3 = 1; ibin3<=hnum->GetNbinsZ(); ibin3++){
                    content+=hnum->GetBinContent(ibin,ibin2,ibin3);
                }
            }
            Printf("ibin = %d, low edge = %0.2f, content = %0.2f",ibin,hnum->GetXaxis()->GetBinLowEdge(ibin), content);
            sum+=content;
            
        }
        Printf("Rebinned THnSparse <-------------");
        for (Int_t ibin = 1; ibin<=Rbhnum->GetNbinsX(); ibin++){
            content2=0;
            for (Int_t ibin2 = 1; ibin2<=Rbhnum->GetNbinsY(); ibin2++){
                for (Int_t ibin3 = 1; ibin3<=Rbhnum->GetNbinsZ(); ibin3++){
                    content2+=Rbhnum->GetBinContent(ibin,ibin2,ibin3);
                }
            }
            Printf("ibin = %d, low edge = %0.2f, content = %0.2f",ibin,Rbhnum->GetXaxis()->GetBinLowEdge(ibin), content2);
            sumnew+=content2;
        }
        Printf("Sum_before = %f, Sum_After = %f, Checking difference.. =%f%s", sum, sumnew, (sum-sumnew)*100/sum,"%");
    }
    
    
    if(((sum-sumnew)*100/sum)<.01){
        TFile* fout = new TFile("STE3D_DhCorrPP_Pass4_03Dec15.root","RECREATE");
        c->Write();
        fout->Close();
        Printf( "1>>> rebin OKEY, STE3D_DhCorrPP_Pass4_03Dec15.root file created");
        CheckFileQuality();

    }
    else Printf(" ======> rebin Failed");
    
    
}


//funciton to check 3D file Quality (e.g. Empty bins or low/high efficiency)
void CheckFileQuality(){

    Printf( "\n2>>> Checking file Quality..");
    TFile fC("STE3D_DhCorrPP_Pass4_03Dec15.root");
    TCanvas *cPass2 = (TCanvas*)fC.Get("c");
    TH3D *hcEffPass2=(TH3D*)cPass2->FindObject("heff_rebin");
    
    cout << "Number of Pt bin = " <<hcEffPass2->GetNbinsX() << endl;
    cout << "Number of Eta bin = " <<hcEffPass2->GetNbinsY() << endl;
    cout << "Number of Zvtx bin = " <<hcEffPass2->GetNbinsZ() << endl;

    Printf( "\n**Suspected bins but can be due to less stats in bin");
    for(int i=1; i<=hcEffPass2->GetNbinsX(); i++){
        for(int j=1; j<=hcEffPass2->GetNbinsY(); j++){
            for(int k=1; k<=hcEffPass2->GetNbinsZ(); k++){
                if(hcEffPass2->GetBinContent(i,j,k) < 0.50 || hcEffPass2->GetBinContent(i,j,k) > 1.00){
                    cout <<"(Xbin = " << i <<", Ybin = "<<j<< ", Zbin = "<< k <<") -----> " <<  hcEffPass2->GetBinContent(i,j,k) << endl;
                }
            }
        }
    }
}


void DrawDhCorr_STE_1DEfficiency_in_pT(){
    
    TFile* f = new TFile("AnalysisResults_pp2017_TrackEff.root");
    TDirectoryFile* d = (TDirectoryFile*)f->Get("PWGPP_CFSingleTrack");
    AliCFContainer *data = (AliCFContainer*) (d->Get("containerpp2017_pi"));
    
    AliCFGridSparse* gridSparsenum = (AliCFGridSparse*)data->GetGrid(1); // GenAcc
    THnSparse* numData = (THnSparse*)gridSparsenum->GetGrid();
    AliCFGridSparse* gridSparseden = (AliCFGridSparse*)data->GetGrid(6); // Reco
    THnSparse* denData = (THnSparse*)gridSparseden->GetGrid();
    
    TFile* f2 = new TFile("AnalysisResults_pp2017_TrackEff.root");
    TDirectoryFile* d2 = (TDirectoryFile*)f2->Get("PWGPP_CFSingleTrack");
    AliCFContainer *data2 = (AliCFContainer*) (d2->Get("containerpp2017_k"));
    
    AliCFGridSparse* gridSparsenum2 = (AliCFGridSparse*)data2->GetGrid(1); // GenAcc
    THnSparse* numData2 = (THnSparse*)gridSparsenum2->GetGrid();
    AliCFGridSparse* gridSparseden2 = (AliCFGridSparse*)data2->GetGrid(6); // Reco
    THnSparse* denData2 = (THnSparse*)gridSparseden2->GetGrid();
    
    TFile* f3 = new TFile("AnalysisResults_pp2017_TrackEff.root");
    TDirectoryFile* d3 = (TDirectoryFile*)f3->Get("PWGPP_CFSingleTrack");
    AliCFContainer *data3 = (AliCFContainer*) (d3->Get("containerpp2017_p"));
    
    AliCFGridSparse* gridSparsenum3 = (AliCFGridSparse*)data3->GetGrid(1); // GenAcc
    THnSparse* numData3 = (THnSparse*)gridSparsenum3->GetGrid();
    AliCFGridSparse* gridSparseden3 = (AliCFGridSparse*)data3->GetGrid(6); // Reco
    THnSparse* denData3 = (THnSparse*)gridSparseden3->GetGrid();
    
    TFile* f4 = new TFile("AnalysisResults_pp2017_TrackEff.root");
    TDirectoryFile* d4 = (TDirectoryFile*)f4->Get("PWGPP_CFSingleTrack");
    AliCFContainer *data4 = (AliCFContainer*) (d4->Get("containerpp2017_e"));
    
    AliCFGridSparse* gridSparsenum4 = (AliCFGridSparse*)data4->GetGrid(1); // GenAcc
    THnSparse* numData4 = (THnSparse*)gridSparsenum4->GetGrid();
    AliCFGridSparse* gridSparseden4 = (AliCFGridSparse*)data4->GetGrid(6); // Reco
    THnSparse* denData4 = (THnSparse*)gridSparseden4->GetGrid();
    
    TFile* f5 = new TFile("AnalysisResults_pp2017_TrackEff.root");
    TDirectoryFile* d5 = (TDirectoryFile*)f5->Get("PWGPP_CFSingleTrack");
    AliCFContainer *data5 = (AliCFContainer*) (d5->Get("containerpp2017_mu"));
    
    AliCFGridSparse* gridSparsenum5 = (AliCFGridSparse*)data5->GetGrid(1); // GenAcc
    THnSparse* numData5 = (THnSparse*)gridSparsenum5->GetGrid();
    AliCFGridSparse* gridSparseden5 = (AliCFGridSparse*)data5->GetGrid(6); // Reco
    THnSparse* denData5 = (THnSparse*)gridSparseden5->GetGrid();
            
    TH1D* hnum = (TH1D*)numData->Projection(0);
    TH1D* hden = (TH1D*)denData->Projection(0);
    
    TH1D* hnum2 = (TH1D*)numData2->Projection(0);
    TH1D* hden2 = (TH1D*)denData2->Projection(0);
    
    TH1D* hnum3 = (TH1D*)numData3->Projection(0);
    TH1D* hden3 = (TH1D*)denData3->Projection(0);
    
    TH1D* hnum4 = (TH1D*)numData4->Projection(0);
    TH1D* hden4 = (TH1D*)denData4->Projection(0);
    
    TH1D* hnum5 = (TH1D*)numData5->Projection(0);
    TH1D* hden5 = (TH1D*)denData5->Projection(0);                    


    TCanvas *c = new TCanvas("c", "pT, Eta, Zvtx Efficiency distrubution", 800,900);
    c->Divide(1,2);
    
    TH1D* hnumTot = (TH1D*)hnum->Clone("hnumtot");
    hnumTot->Add(hnum2);
    hnumTot->Add(hnum3);
    hnumTot->Add(hnum4);
    hnumTot->Add(hnum5);
    TH1D* hdenTot = (TH1D*)hden->Clone("hdentot");
    hdenTot->Add(hden2);
    hdenTot->Add(hden3);
    hdenTot->Add(hden4);
    hdenTot->Add(hden5);
    TH1D* heff = (TH1D*)hdenTot->Clone("heff");
    heff->Divide(hdenTot,hnumTot,1,1,"B");
    
    c->cd(1);
    heff->SetMarkerStyle(21);
    heff->SetMarkerSize(0.8);
    heff->Draw();
    //Rebinning
    const Int_t newLimitsBins = 26;
    Double_t* newLimits = new Double_t[newLimitsBins+1];
    newLimits[0]   = 0.30;
    newLimits[1]   = 0.40;
    newLimits[2]   = 0.50;
    newLimits[3]   = 0.60;
    newLimits[4]   = 0.70;
    newLimits[5]   = 0.80;
    newLimits[6]   = 0.90;
    newLimits[7]   = 1.00;
    newLimits[8]   = 1.20;
    newLimits[9]   = 1.40;
    newLimits[10]  = 1.60;
    newLimits[11]  = 1.80;
    newLimits[12]  = 2.00;
    newLimits[13]  = 2.33;
    newLimits[14]  = 2.66;
    newLimits[15]  = 3.00;
    newLimits[16]  = 3.50;
    newLimits[17]  = 4.00;
    newLimits[18]  = 4.40;
    newLimits[19]  = 4.80;
    newLimits[20]  = 5.20;
    newLimits[21]  = 6.00;
    newLimits[22]  = 7.00;
    newLimits[23]  = 8.00;
    newLimits[24]  = 10.00;
    newLimits[25]  = 16.00;
    newLimits[26]  = 24.00;
    
    THnSparse* RbnumData = (THnSparse*)numData->Clone("numNew");
    RbnumData->Reset();
    TAxis* axis = (TAxis*)RbnumData->GetAxis(0);
    axis->Set(newLimitsBins,newLimits);
    RbnumData->SetBinEdges(0,newLimits);
    RbnumData->RebinnedAdd(numData, 1);
    TH1D* Rbhnum = (TH1D*)RbnumData->Projection(0);
    gStyle->SetOptStat(kTRUE);
    
    THnSparse* RbdenData = (THnSparse*)denData->Clone("denNew");
    RbdenData->Reset();
    TAxis* axis2 = (TAxis*)RbdenData->GetAxis(0);
    axis2->Set(newLimitsBins,newLimits);
    RbdenData->SetBinEdges(0,newLimits);
    RbdenData->RebinnedAdd(denData, 1);
    TH1D* Rbhden = (TH1D*)RbdenData->Projection(0);
    gStyle->SetOptStat(kTRUE);
    
    THnSparse* RbnumData2 = (THnSparse*)numData2->Clone("numNew2");
    RbnumData2->Reset();
    TAxis* axis = (TAxis*)RbnumData2->GetAxis(0);
    axis->Set(newLimitsBins,newLimits);
    RbnumData2->SetBinEdges(0,newLimits);
    RbnumData2->RebinnedAdd(numData2, 1);
    TH1D* Rbhnum2 = (TH1D*)RbnumData2->Projection(0);
    gStyle->SetOptStat(kTRUE);
    
    THnSparse* RbdenData2 = (THnSparse*)denData2->Clone("denNew2");
    RbdenData2->Reset();
    TAxis* axis2 = (TAxis*)RbdenData2->GetAxis(0);
    axis2->Set(newLimitsBins,newLimits);
    RbdenData2->SetBinEdges(0,newLimits);
    RbdenData2->RebinnedAdd(denData2, 1);
    TH1D* Rbhden2 = (TH1D*)RbdenData2->Projection(0);
    gStyle->SetOptStat(kTRUE);
    
    THnSparse* RbnumData3 = (THnSparse*)numData3->Clone("numNew3");
    RbnumData3->Reset();
    TAxis* axis = (TAxis*)RbnumData3->GetAxis(0);
    axis->Set(newLimitsBins,newLimits);
    RbnumData3->SetBinEdges(0,newLimits);
    RbnumData3->RebinnedAdd(numData3, 1);
    TH1D* Rbhnum3 = (TH1D*)RbnumData3->Projection(0);
    gStyle->SetOptStat(kTRUE);
    
    THnSparse* RbdenData3 = (THnSparse*)denData3->Clone("denNew3");
    RbdenData3->Reset();
    TAxis* axis2 = (TAxis*)RbdenData3->GetAxis(0);
    axis2->Set(newLimitsBins,newLimits);
    RbdenData3->SetBinEdges(0,newLimits);
    RbdenData3->RebinnedAdd(denData3, 1);
    TH1D* Rbhden3 = (TH1D*)RbdenData3->Projection(0);
    gStyle->SetOptStat(kTRUE);
    
    THnSparse* RbnumData4 = (THnSparse*)numData4->Clone("numNew4");
    RbnumData4->Reset();
    TAxis* axis = (TAxis*)RbnumData4->GetAxis(0);
    axis->Set(newLimitsBins,newLimits);
    RbnumData4->SetBinEdges(0,newLimits);
    RbnumData4->RebinnedAdd(numData4, 1);
    TH1D* Rbhnum4 = (TH1D*)RbnumData4->Projection(0);
    gStyle->SetOptStat(kTRUE);
    
    THnSparse* RbdenData4 = (THnSparse*)denData4->Clone("denNew4");
    RbdenData4->Reset();
    TAxis* axis2 = (TAxis*)RbdenData4->GetAxis(0);
    axis2->Set(newLimitsBins,newLimits);
    RbdenData4->SetBinEdges(0,newLimits);
    RbdenData4->RebinnedAdd(denData4, 1);
    TH1D* Rbhden4 = (TH1D*)RbdenData4->Projection(0);
    gStyle->SetOptStat(kTRUE);
    
    THnSparse* RbnumData5 = (THnSparse*)numData5->Clone("numNew5");
    RbnumData5->Reset();
    TAxis* axis = (TAxis*)RbnumData5->GetAxis(0);
    axis->Set(newLimitsBins,newLimits);
    RbnumData5->SetBinEdges(0,newLimits);
    RbnumData5->RebinnedAdd(numData5, 1);
    TH1D* Rbhnum5 = (TH1D*)RbnumData5->Projection(0);
    gStyle->SetOptStat(kTRUE);
    
    THnSparse* RbdenData5 = (THnSparse*)denData5->Clone("denNew5");
    RbdenData5->Reset();
    TAxis* axis2 = (TAxis*)RbdenData5->GetAxis(0);
    axis2->Set(newLimitsBins,newLimits);
    RbdenData5->SetBinEdges(0,newLimits);
    RbdenData5->RebinnedAdd(denData5, 1);
    TH1D* Rbhden5 = (TH1D*)RbdenData5->Projection(0);
    gStyle->SetOptStat(kTRUE);
    
    TH1D* RbhdenTot = (TH1D*)Rbhden->Clone("Rbhdentot");
    RbhdenTot->Add(Rbhden2);
    RbhdenTot->Add(Rbhden3);
    RbhdenTot->Add(Rbhden4);
    RbhdenTot->Add(Rbhden5);
    TH1D* RbhnumTot = (TH1D*)Rbhnum->Clone("Rbhnumtot");
    RbhnumTot->Add(Rbhnum2);
    RbhnumTot->Add(Rbhnum3);
    RbhnumTot->Add(Rbhnum4);
    RbhnumTot->Add(Rbhnum5);
    
    
    TH1D* Rbheff = (TH1D*)RbhdenTot->Clone("heff_rebin");
    Rbheff->Divide(RbhdenTot,RbhnumTot,1,1,"B");
    gStyle->SetOptStat(kTRUE);
    c->cd(2);
    Rbheff->SetMarkerStyle(21);
    Rbheff->SetMarkerSize(0.8);
    Rbheff->Draw();    
    
    // Num bins checking
    Bool_t Rebintest = kTRUE;
    Float_t sum = 0, sumnew = 0;
    Float_t content = 0, content2 = 0;
    if(Rebintest){
        Printf("Original THnSparse <-------------");
        for (Int_t ibin = 1; ibin<=hnum->GetNbinsX(); ibin++){
            content=0;
            for (Int_t ibin2 = 1; ibin2<=hnum->GetNbinsY(); ibin2++){
                for (Int_t ibin3 = 1; ibin3<=hnum->GetNbinsZ(); ibin3++){
                    content+=hnum->GetBinContent(ibin,ibin2,ibin3);
                }
            }
            Printf("ibin = %d, low edge = %0.2f, content = %0.2f",ibin,hnum->GetXaxis()->GetBinLowEdge(ibin), content);
            sum+=content;
            
        }
        Printf("Rebinned THnSparse <-------------");
        for (Int_t ibin = 1; ibin<=Rbhnum->GetNbinsX(); ibin++){
            content2=0;
            for (Int_t ibin2 = 1; ibin2<=Rbhnum->GetNbinsY(); ibin2++){
                for (Int_t ibin3 = 1; ibin3<=Rbhnum->GetNbinsZ(); ibin3++){
                    content2+=Rbhnum->GetBinContent(ibin,ibin2,ibin3);
                }
            }
            Printf("ibin = %d, low edge = %0.2f, content = %0.2f",ibin,Rbhnum->GetXaxis()->GetBinLowEdge(ibin), content2);
            sumnew+=content2;
        }
        Printf("Sum_before = %f, Sum_After = %f, Checking difference.. =%f%s", sum, sumnew, (sum-sumnew)*100/sum,"%");
    }
    
    
    if(((sum-sumnew)*100/sum)<.01){
        TFile* fout = new TFile("STE3D_DhCorrPP_Pass4_03Dec15.root","RECREATE");
        c->Write();
        fout->Close();
        Printf( "1>>> rebin OKEY, STE3D_DhCorrPP_Pass4_03Dec15.root file created");
        CheckFileQuality();

    }
    else Printf(" ======> rebin Failed");
    
    
}


/*
//additional lines 
 1. macro compliation
 root [0] gSystem->AddIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT -I$ALICE_ROOT/include -I$ALICE_PHYSICS -I$ALICE_PHYSICS/include -I$ALICE_ROOT/ITS -I$ALICE_ROOT/TPC -I$ALICE_ROOT/CONTAINERS -I$ALICE_ROOT/STEER/STEER -I$ALICE_ROOT/STEER/STEERBase -I$ALICE_ROOT/STEER/ESD -I$ALICE_ROOT/STEER/AOD -I$ALICE_ROOT/TRD -I$ALICE_PHYSICS/macros -I$ALICE_PHYSICS/ANALYSIS  -I$ALICE_PHYSICS/OADB -I$ALICE_PHYSICS/PWGHF/ -I$ALICE_PHYSICS/PWGHF/base -I$ALICE_PHYSICS/PWG/FLOW/Base -I$ALICE_PHYSICS/PWG/FLOW/Tasks -I$ALICE_PHYSICS/PWGHF/vertexingHF -I$ALICE_PHYSICS/PWGHF/correlationHF -I$ALICE_PHYSICS/PWGHF/correlationHF/macros  -g");
 root [1] gROOT->LoadMacro("$ALICE_PHYSICS/PWGHF/vertexingHF/macros/LoadLibraries.C")
 root [2] .L DrawDhCorr_STE_3DEfficiency_in_pTEtaZvtx.C++
 
 
 2. Getting histogram from file
 TFile* fAssoTracksEffMap = TFile::Open("STE3D_DhCorrPP_Pass4_03Dec15.root");
 TCanvas *cTrackEffMap = (TCanvas*)fAssoTracksEffMap->Get("c");
 TCanvas *c = new TCanvas ("c", "Eff histogram from root file", 500, 600);
 TH3D *hEffTrackMap = (TH3D*)cTrackEffMap->FindObject("heff_rebin");
 hEffTrackMap->Draw();
 */

void UserBins_DrawDhCorr_STE_1DEfficiency_in_pT(){
    
    TFile* f = new TFile("AnalysisResults_f2b_CENT_woSDD.root");
    TDirectoryFile* d = (TDirectoryFile*)f->Get("PWGPP_CFSingleTrack");
    AliCFContainer *data = (AliCFContainer*) (d->Get("containerpp2017_pi"));
    
    AliCFGridSparse* gridSparsenum = (AliCFGridSparse*)data->GetGrid(1); // GenAcc
    THnSparse* numData = (THnSparse*)gridSparsenum->GetGrid();
    AliCFGridSparse* gridSparseden = (AliCFGridSparse*)data->GetGrid(6); // Reco
    THnSparse* denData = (THnSparse*)gridSparseden->GetGrid();
    
    TFile* f2 = new TFile("AnalysisResults_f2b_CENT_woSDD.root");
    TDirectoryFile* d2 = (TDirectoryFile*)f2->Get("PWGPP_CFSingleTrack");
    AliCFContainer *data2 = (AliCFContainer*) (d2->Get("containerpp2017_k"));
    
    AliCFGridSparse* gridSparsenum2 = (AliCFGridSparse*)data2->GetGrid(1); // GenAcc
    THnSparse* numData2 = (THnSparse*)gridSparsenum2->GetGrid();
    AliCFGridSparse* gridSparseden2 = (AliCFGridSparse*)data2->GetGrid(6); // Reco
    THnSparse* denData2 = (THnSparse*)gridSparseden2->GetGrid();
    
    TFile* f3 = new TFile("AnalysisResults_f2b_CENT_woSDD.root");
    TDirectoryFile* d3 = (TDirectoryFile*)f3->Get("PWGPP_CFSingleTrack");
    AliCFContainer *data3 = (AliCFContainer*) (d3->Get("containerpp2017_p"));
    
    AliCFGridSparse* gridSparsenum3 = (AliCFGridSparse*)data3->GetGrid(1); // GenAcc
    THnSparse* numData3 = (THnSparse*)gridSparsenum3->GetGrid();
    AliCFGridSparse* gridSparseden3 = (AliCFGridSparse*)data3->GetGrid(6); // Reco
    THnSparse* denData3 = (THnSparse*)gridSparseden3->GetGrid();
    
    TFile* f4 = new TFile("AnalysisResults_f2b_CENT_woSDD.root");
    TDirectoryFile* d4 = (TDirectoryFile*)f4->Get("PWGPP_CFSingleTrack");
    AliCFContainer *data4 = (AliCFContainer*) (d4->Get("containerpp2017_e"));
    
    AliCFGridSparse* gridSparsenum4 = (AliCFGridSparse*)data4->GetGrid(1); // GenAcc
    THnSparse* numData4 = (THnSparse*)gridSparsenum4->GetGrid();
    AliCFGridSparse* gridSparseden4 = (AliCFGridSparse*)data4->GetGrid(6); // Reco
    THnSparse* denData4 = (THnSparse*)gridSparseden4->GetGrid();
    
    TFile* f5 = new TFile("AnalysisResults_f2b_CENT_woSDD.root");
    TDirectoryFile* d5 = (TDirectoryFile*)f5->Get("PWGPP_CFSingleTrack");
    AliCFContainer *data5 = (AliCFContainer*) (d5->Get("containerpp2017_mu"));
    
    AliCFGridSparse* gridSparsenum5 = (AliCFGridSparse*)data5->GetGrid(1); // GenAcc
    THnSparse* numData5 = (THnSparse*)gridSparsenum5->GetGrid();
    AliCFGridSparse* gridSparseden5 = (AliCFGridSparse*)data5->GetGrid(6); // Reco
    THnSparse* denData5 = (THnSparse*)gridSparseden5->GetGrid();
            
    TH1D* hnum = (TH1D*)numData->Projection(0);
    TH1D* hden = (TH1D*)denData->Projection(0);
    
    TH1D* hnum2 = (TH1D*)numData2->Projection(0);
    TH1D* hden2 = (TH1D*)denData2->Projection(0);
    
    TH1D* hnum3 = (TH1D*)numData3->Projection(0);
    TH1D* hden3 = (TH1D*)denData3->Projection(0);
    
    TH1D* hnum4 = (TH1D*)numData4->Projection(0);
    TH1D* hden4 = (TH1D*)denData4->Projection(0);
    
    TH1D* hnum5 = (TH1D*)numData5->Projection(0);
    TH1D* hden5 = (TH1D*)denData5->Projection(0);                    


    TCanvas *c = new TCanvas("c", "pT, Eta, Zvtx Efficiency distrubution", 400,900);
    c->Divide(1,2);
    
    TH1D* hnumTot = (TH1D*)hnum->Clone("hnumtot");
    hnumTot->Add(hnum2);
    hnumTot->Add(hnum3);
    hnumTot->Add(hnum4);
    hnumTot->Add(hnum5);
    TH1D* hdenTot = (TH1D*)hden->Clone("hdentot");
    hdenTot->Add(hden2);
    hdenTot->Add(hden3);
    hdenTot->Add(hden4);
    hdenTot->Add(hden5);
    TH1D* heff = (TH1D*)hdenTot->Clone("heff");
    heff->Divide(hdenTot,hnumTot,1,1,"B");
    
    c->cd(1);
    heff->Draw();
    //Rebinning
    const Int_t newLimitsBins = 4;
    Double_t* newLimits = new Double_t[newLimitsBins+1];
    newLimits[0]  = 0.30;
    newLimits[1]  = 1.00;
    newLimits[2]  = 2.00;
    newLimits[3]  = 3.00;
    newLimits[4]  = 24.00;
    
    THnSparse* RbnumData = (THnSparse*)numData->Clone("numNew");
    RbnumData->Reset();
    TAxis* axis = (TAxis*)RbnumData->GetAxis(0);
    axis->Set(newLimitsBins,newLimits);
    RbnumData->SetBinEdges(0,newLimits);
    RbnumData->RebinnedAdd(numData, 1);
    TH1D* Rbhnum = (TH1D*)RbnumData->Projection(0);
    gStyle->SetOptStat(kTRUE);
    
    THnSparse* RbdenData = (THnSparse*)denData->Clone("denNew");
    RbdenData->Reset();
    TAxis* axis2 = (TAxis*)RbdenData->GetAxis(0);
    axis2->Set(newLimitsBins,newLimits);
    RbdenData->SetBinEdges(0,newLimits);
    RbdenData->RebinnedAdd(denData, 1);
    TH1D* Rbhden = (TH1D*)RbdenData->Projection(0);
    gStyle->SetOptStat(kTRUE);
    
    THnSparse* RbnumData2 = (THnSparse*)numData2->Clone("numNew2");
    RbnumData2->Reset();
    TAxis* axis = (TAxis*)RbnumData2->GetAxis(0);
    axis->Set(newLimitsBins,newLimits);
    RbnumData2->SetBinEdges(0,newLimits);
    RbnumData2->RebinnedAdd(numData2, 1);
    TH1D* Rbhnum2 = (TH1D*)RbnumData2->Projection(0);
    gStyle->SetOptStat(kTRUE);
    
    THnSparse* RbdenData2 = (THnSparse*)denData2->Clone("denNew2");
    RbdenData2->Reset();
    TAxis* axis2 = (TAxis*)RbdenData2->GetAxis(0);
    axis2->Set(newLimitsBins,newLimits);
    RbdenData2->SetBinEdges(0,newLimits);
    RbdenData2->RebinnedAdd(denData2, 1);
    TH1D* Rbhden2 = (TH1D*)RbdenData2->Projection(0);
    gStyle->SetOptStat(kTRUE);
    
    THnSparse* RbnumData3 = (THnSparse*)numData3->Clone("numNew3");
    RbnumData3->Reset();
    TAxis* axis = (TAxis*)RbnumData3->GetAxis(0);
    axis->Set(newLimitsBins,newLimits);
    RbnumData3->SetBinEdges(0,newLimits);
    RbnumData3->RebinnedAdd(numData3, 1);
    TH1D* Rbhnum3 = (TH1D*)RbnumData3->Projection(0);
    gStyle->SetOptStat(kTRUE);
    
    THnSparse* RbdenData3 = (THnSparse*)denData3->Clone("denNew3");
    RbdenData3->Reset();
    TAxis* axis2 = (TAxis*)RbdenData3->GetAxis(0);
    axis2->Set(newLimitsBins,newLimits);
    RbdenData3->SetBinEdges(0,newLimits);
    RbdenData3->RebinnedAdd(denData3, 1);
    TH1D* Rbhden3 = (TH1D*)RbdenData3->Projection(0);
    gStyle->SetOptStat(kTRUE);
    
    THnSparse* RbnumData4 = (THnSparse*)numData4->Clone("numNew4");
    RbnumData4->Reset();
    TAxis* axis = (TAxis*)RbnumData4->GetAxis(0);
    axis->Set(newLimitsBins,newLimits);
    RbnumData4->SetBinEdges(0,newLimits);
    RbnumData4->RebinnedAdd(numData4, 1);
    TH1D* Rbhnum4 = (TH1D*)RbnumData4->Projection(0);
    gStyle->SetOptStat(kTRUE);
    
    THnSparse* RbdenData4 = (THnSparse*)denData4->Clone("denNew4");
    RbdenData4->Reset();
    TAxis* axis2 = (TAxis*)RbdenData4->GetAxis(0);
    axis2->Set(newLimitsBins,newLimits);
    RbdenData4->SetBinEdges(0,newLimits);
    RbdenData4->RebinnedAdd(denData4, 1);
    TH1D* Rbhden4 = (TH1D*)RbdenData4->Projection(0);
    gStyle->SetOptStat(kTRUE);
    
    THnSparse* RbnumData5 = (THnSparse*)numData5->Clone("numNew5");
    RbnumData5->Reset();
    TAxis* axis = (TAxis*)RbnumData5->GetAxis(0);
    axis->Set(newLimitsBins,newLimits);
    RbnumData5->SetBinEdges(0,newLimits);
    RbnumData5->RebinnedAdd(numData5, 1);
    TH1D* Rbhnum5 = (TH1D*)RbnumData5->Projection(0);
    gStyle->SetOptStat(kTRUE);
    
    THnSparse* RbdenData5 = (THnSparse*)denData5->Clone("denNew5");
    RbdenData5->Reset();
    TAxis* axis2 = (TAxis*)RbdenData5->GetAxis(0);
    axis2->Set(newLimitsBins,newLimits);
    RbdenData5->SetBinEdges(0,newLimits);
    RbdenData5->RebinnedAdd(denData5, 1);
    TH1D* Rbhden5 = (TH1D*)RbdenData5->Projection(0);
    gStyle->SetOptStat(kTRUE);
    
    TH1D* RbhdenTot = (TH1D*)Rbhden->Clone("Rbhdentot");
    RbhdenTot->Add(Rbhden2);
    RbhdenTot->Add(Rbhden3);
    RbhdenTot->Add(Rbhden4);
    RbhdenTot->Add(Rbhden5);
    TH1D* RbhnumTot = (TH1D*)Rbhnum->Clone("Rbhnumtot");
    RbhnumTot->Add(Rbhnum2);
    RbhnumTot->Add(Rbhnum3);
    RbhnumTot->Add(Rbhnum4);
    RbhnumTot->Add(Rbhnum5);
    
    
    TH1D* Rbheff = (TH1D*)RbhdenTot->Clone("heff_rebin");
    Rbheff->Divide(RbhdenTot,RbhnumTot,1,1,"B");
    gStyle->SetOptStat(kTRUE);
    c->cd(2);
    Rbheff->Draw();
    
    for(int i = 1; i<=newLimitsBins; i++) printf("Bin %d, Eff = %1.5f\n",i,Rbheff->GetBinContent(i));
    
    // Num bins checking
    Bool_t Rebintest = kTRUE;
    Float_t sum = 0, sumnew = 0;
    Float_t content = 0, content2 = 0;
    if(Rebintest){
        Printf("Original THnSparse <-------------");
        for (Int_t ibin = 1; ibin<=hnum->GetNbinsX(); ibin++){
            content=0;
            for (Int_t ibin2 = 1; ibin2<=hnum->GetNbinsY(); ibin2++){
                for (Int_t ibin3 = 1; ibin3<=hnum->GetNbinsZ(); ibin3++){
                    content+=hnum->GetBinContent(ibin,ibin2,ibin3);
                }
            }
            Printf("ibin = %d, low edge = %0.2f, content = %0.2f",ibin,hnum->GetXaxis()->GetBinLowEdge(ibin), content);
            sum+=content;
            
        }
        Printf("Rebinned THnSparse <-------------");
        for (Int_t ibin = 1; ibin<=Rbhnum->GetNbinsX(); ibin++){
            content2=0;
            for (Int_t ibin2 = 1; ibin2<=Rbhnum->GetNbinsY(); ibin2++){
                for (Int_t ibin3 = 1; ibin3<=Rbhnum->GetNbinsZ(); ibin3++){
                    content2+=Rbhnum->GetBinContent(ibin,ibin2,ibin3);
                }
            }
            Printf("ibin = %d, low edge = %0.2f, content = %0.2f",ibin,Rbhnum->GetXaxis()->GetBinLowEdge(ibin), content2);
            sumnew+=content2;
        }
        Printf("Sum_before = %f, Sum_After = %f, Checking difference.. =%f%s", sum, sumnew, (sum-sumnew)*100/sum,"%");
    }
    
    
    if(((sum-sumnew)*100/sum)<.01){
        TFile* fout = new TFile("STE3D_DhCorrPP_Pass4_03Dec15.root","RECREATE");
        c->Write();
        fout->Close();
        Printf( "1>>> rebin OKEY, STE3D_DhCorrPP_Pass4_03Dec15.root file created");
        CheckFileQuality();

    }
    else Printf(" ======> rebin Failed");

    
}


void DrawRelativeAbundances_VsPt(){
    
    TFile* ftot = new TFile("AnalysisResults_pp2017_TrackEff.root");
    TDirectoryFile* dtot = (TDirectoryFile*)ftot->Get("PWGPP_CFSingleTrack");
    AliCFContainer *datatot = (AliCFContainer*) (dtot->Get("containerpp2017"));
    
    AliCFGridSparse* gridSparsenum_tot = (AliCFGridSparse*)datatot->GetGrid(1); // GenAcc
    THnSparse* numData_tot = (THnSparse*)gridSparsenum_tot->GetGrid();

    TFile* f = new TFile("AnalysisResults_pp2017_TrackEff.root");
    TDirectoryFile* d = (TDirectoryFile*)f->Get("PWGPP_CFSingleTrack");
    AliCFContainer *data = (AliCFContainer*) (d->Get("containerpp2017_pi"));
    
    AliCFGridSparse* gridSparsenum = (AliCFGridSparse*)data->GetGrid(1); // GenAcc
    THnSparse* numData = (THnSparse*)gridSparsenum->GetGrid();
    
    TFile* f2 = new TFile("AnalysisResults_pp2017_TrackEff.root");
    TDirectoryFile* d2 = (TDirectoryFile*)f2->Get("PWGPP_CFSingleTrack");
    AliCFContainer *data2 = (AliCFContainer*) (d2->Get("containerpp2017_k"));
    
    AliCFGridSparse* gridSparsenum2 = (AliCFGridSparse*)data2->GetGrid(1); // GenAcc
    THnSparse* numData2 = (THnSparse*)gridSparsenum2->GetGrid();
    
    TFile* f3 = new TFile("AnalysisResults_pp2017_TrackEff.root");
    TDirectoryFile* d3 = (TDirectoryFile*)f3->Get("PWGPP_CFSingleTrack");
    AliCFContainer *data3 = (AliCFContainer*) (d3->Get("containerpp2017_p"));
    
    AliCFGridSparse* gridSparsenum3 = (AliCFGridSparse*)data3->GetGrid(1); // GenAcc
    THnSparse* numData3 = (THnSparse*)gridSparsenum3->GetGrid();
    
    TFile* f4 = new TFile("AnalysisResults_pp2017_TrackEff.root");
    TDirectoryFile* d4 = (TDirectoryFile*)f4->Get("PWGPP_CFSingleTrack");
    AliCFContainer *data4 = (AliCFContainer*) (d4->Get("containerpp2017_e"));
    
    AliCFGridSparse* gridSparsenum4 = (AliCFGridSparse*)data4->GetGrid(1); // GenAcc
    THnSparse* numData4 = (THnSparse*)gridSparsenum4->GetGrid();
    
    TFile* f5 = new TFile("AnalysisResults_pp2017_TrackEff.root");
    TDirectoryFile* d5 = (TDirectoryFile*)f5->Get("PWGPP_CFSingleTrack");
    AliCFContainer *data5 = (AliCFContainer*) (d5->Get("containerpp2017_mu"));
    
    AliCFGridSparse* gridSparsenum5 = (AliCFGridSparse*)data5->GetGrid(1); // GenAcc
    THnSparse* numData5 = (THnSparse*)gridSparsenum5->GetGrid();
            
    TH1D* hnum_total = (TH1D*)numData_tot->Projection(0);
    TH1D* hnum_pi = (TH1D*)numData->Projection(0); 
    TH1D* hnum_k  = (TH1D*)numData2->Projection(0);    
    TH1D* hnum_p  = (TH1D*)numData3->Projection(0);    
    TH1D* hnum_e  = (TH1D*)numData4->Projection(0);    
    TH1D* hnum_mu = (TH1D*)numData5->Projection(0);

    TCanvas *c = new TCanvas();
    hnum_pi->Divide(hnum_pi,hnum_total,1,1,"B");
    hnum_pi->SetMarkerColor(kRed);
    hnum_pi->SetLineColor(kRed);
    hnum_pi->SetMarkerStyle(21);
    hnum_pi->SetMarkerSize(0.8);
    hnum_pi->SetMinimum(0.0001);
    hnum_pi->SetMaximum(1.8);
    hnum_pi->Draw();
    hnum_k->Divide(hnum_k,hnum_total,1,1,"B");
    hnum_k->SetMarkerColor(kBlue);
    hnum_k->SetLineColor(kBlue);
    hnum_k->SetMarkerStyle(21);
    hnum_k->SetMarkerSize(0.8);    
    hnum_k->Draw("same");
    hnum_p->Divide(hnum_p,hnum_total,1,1,"B");
    hnum_p->SetMarkerColor(kGreen+1);
    hnum_p->SetLineColor(kGreen+1);
    hnum_p->SetMarkerStyle(21);
    hnum_p->SetMarkerSize(0.8);    
    hnum_p->Draw("same");
    hnum_e->Divide(hnum_e,hnum_total,1,1,"B");
    hnum_e->SetMarkerColor(kMagenta);
    hnum_e->SetLineColor(kMagenta);
    hnum_e->SetMarkerStyle(21);
    hnum_e->SetMarkerSize(0.8);    
    hnum_e->Draw("same");
    hnum_mu->Divide(hnum_mu,hnum_total,1,1,"B");
    hnum_mu->SetMarkerColor(kCyan+1);
    hnum_mu->SetLineColor(kCyan+1);
    hnum_mu->SetMarkerStyle(21);
    hnum_mu->SetMarkerSize(0.8);    
    hnum_mu->Draw("same");

    TLegend* leg=new TLegend(0.65,0.15,0.95,0.40);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->AddEntry(hnum_pi,"pions","epl");
    leg->AddEntry(hnum_k,"kaons","epl");
    leg->AddEntry(hnum_p,"protons","epl");
    leg->AddEntry(hnum_e,"electrons","epl");
    leg->AddEntry(hnum_mu,"muons","epl");
    leg->Draw(); 

    c->SetLogy();
    c->Draw();

}