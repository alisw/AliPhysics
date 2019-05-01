

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
    
    TFile* f = new TFile("AnalysisResults_TrackEffOnHFenr.root");
    TDirectoryFile* d = (TDirectoryFile*)f->Get("PWGPP_CFSingleTrack");
    AliCFContainer *data = (AliCFContainer*) (d->Get("containerMin2ITSCls"));
    
    AliCFGridSparse* gridSparsenum = (AliCFGridSparse*)data->GetGrid(1); // GenAcc
    THnSparse* numData = (THnSparse*)gridSparsenum->GetGrid();
    AliCFGridSparse* gridSparseden = (AliCFGridSparse*)data->GetGrid(6); // Reco
    THnSparse* denData = (THnSparse*)gridSparseden->GetGrid();
    
    TCanvas *c = new TCanvas("c", "pT, Eta, Zvtx Efficiency distrubution", 400,900);
    c->Divide(1,2);
    
    TH3D* hnum = (TH3D*)numData->Projection(0,1,4);
    TH3D* hden = (TH3D*)denData->Projection(0,1,4);
    TH3D* heff = (TH3D*)hden->Clone("heff");
    heff->Divide(hden,hnum,1,1,"B");
    
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
    newLimits[13]  = 2.25;
    newLimits[14]  = 2.50;
    newLimits[15]  = 2.75;
    newLimits[16]  = 3.00;
    newLimits[17]  = 3.50;
    newLimits[18]  = 4.00;
    newLimits[19]  = 4.50;
    newLimits[20]  = 5.50;
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
    
    TH3D* Rbheff = (TH3D*)Rbhden->Clone("heff_rebin");
    Rbheff->Divide(Rbhden,Rbhnum,1,1,"B");
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
    
    TFile* f = new TFile("AnalysisResults_f2a_CENT_woSDD.root");
    TDirectoryFile* d = (TDirectoryFile*)f->Get("PWGPP_CFSingleTrack");
    AliCFContainer *data = (AliCFContainer*) (d->Get("containerprelSQM"));
    
    AliCFGridSparse* gridSparsenum = (AliCFGridSparse*)data->GetGrid(1); // GenAcc
    THnSparse* numData = (THnSparse*)gridSparsenum->GetGrid();
    AliCFGridSparse* gridSparseden = (AliCFGridSparse*)data->GetGrid(6); // Reco
    THnSparse* denData = (THnSparse*)gridSparseden->GetGrid();
    
    TCanvas *c = new TCanvas("c", "pT, Eta, Zvtx Efficiency distrubution", 400,900);
    c->Divide(1,2);
    
    TH1D* hnum = (TH1D*)numData->Projection(0);
    TH1D* hden = (TH1D*)denData->Projection(0);
    TH1D* heff = (TH1D*)hden->Clone("heff");
    heff->Divide(hden,hnum,1,1,"B");
    
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
    newLimits[13]  = 2.25;
    newLimits[14]  = 2.50;
    newLimits[15]  = 2.75;
    newLimits[16]  = 3.00;
    newLimits[17]  = 3.50;
    newLimits[18]  = 4.00;
    newLimits[19]  = 4.50;
    newLimits[20]  = 5.50;
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
    
    TH1D* Rbheff = (TH1D*)Rbhden->Clone("heff_rebin");
    Rbheff->Divide(Rbhden,Rbhnum,1,1,"B");
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

