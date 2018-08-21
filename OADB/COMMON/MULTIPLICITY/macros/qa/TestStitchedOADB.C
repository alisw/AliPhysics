#if !defined (__CINT__) || (defined(__MAKECINT__))
#include <iostream>
#include "TROOT.h"
#include "AliMultSelectionCuts.h"
#include "AliOADBMultSelection.h"
#include "TH1F.h"
#include "TFile.h"
#include "AliOADBContainer.h"
#endif

void TestStitchedOADB(TString lPeriodName = "LHC18c",
                      Long_t         lRun = 285957 ) {
    // Code to test OADB in Offline Mode
    // Further developments needed!
    //
    // Todo: Loop on all existing estimators, add flatness histograms
    
    TString gLibs[] =   {"STEER",
        "ANALYSIS", "ANALYSISalice", "ANALYSIScalib", "OADB"};
    TString thislib = "lib";
    for(Int_t ilib = 0; ilib<4; ilib++){
        thislib="lib";
        thislib.Append(gLibs[ilib].Data());
        cout<<"Will load "<<thislib.Data()<<endl;
        gSystem->Load(thislib.Data());
    }
    gSystem->SetIncludePath("-I$ROOTSYS/include  -I$ALICE_ROOT/include -I$ALICE_PHYSICS/include");
    
    cout<<"Alive! "<<endl;
    
    TFile * f = new TFile (Form("OADB-%s.root", lPeriodName.Data()));
    AliOADBContainer * oadbContMS = (AliOADBContainer*) f->Get("MultSel");
    AliOADBMultSelection * oadbMultSelection = (AliOADBMultSelection*)oadbContMS->GetObject(lRun, "Default");
    
    oadbMultSelection->GetEventCuts()->Print();
    oadbMultSelection->GetMultSelection()->PrintInfo();

    Double_t anchor_percentile = ((AliMultEstimator*)oadbMultSelection->GetMultSelection()->GetEstimator("V0M"))->GetAnchorPercentile();
    
    //Attempt to open object for test run
    Long_t lTestRun = lRun;
    cout<<"Load OADB..."<<endl;
    AliOADBMultSelection * oadbMultSelection = (AliOADBMultSelection*)oadbContMS->GetObject(lTestRun, "Default");
    cout<<"OADB loaded for run "<<lTestRun<<endl;
    cout<<"Anchor percentile for this run: "<<anchor_percentile<<endl;
    

    //FIXME: Substitute with loop over several estimators, as needed
    TH1F *hCalibV0M = oadbMultSelection->GetCalibHisto(0);
    
    AliMultSelection *fSelection = (AliMultSelection*) oadbMultSelection->GetMultSelection();
    
    fSelection->PrintInfo();
    
    TCanvas *c1 = new TCanvas("c1", "",800,600);
    hCalibV0M->Draw(); 
    
    //
    // set percentile boundaries for the estimator histos
    // (based on what is implemented in the calibration)
    Double_t lDesiredBoundaries[1000];
    Long_t   lNDesiredBoundaries=0;
    lDesiredBoundaries[0] = 0.0;
    //From High To Low Multiplicity
    if(kFALSE) { // wide binning
        for( Int_t ib = 1; ib < 101; ib++) { // 100 bins  ] 0.0 , 0.1 ]
            lNDesiredBoundaries++;
            lDesiredBoundaries[lNDesiredBoundaries] = lDesiredBoundaries[lNDesiredBoundaries-1] + 0.01;
            // cout << "loop 1: " << lDesiredBoundaries[lNDesiredBoundaries] << endl;
        }
        for( Int_t ib = 1; ib < 91; ib++) { // 90 bins ] 1.0 , 10. ]
            lNDesiredBoundaries++;
            lDesiredBoundaries[lNDesiredBoundaries] = lDesiredBoundaries[lNDesiredBoundaries-1] + 0.1;
            // cout << "loop 2: " << lDesiredBoundaries[lNDesiredBoundaries] << endl;
        }
        for( Int_t ib = 1; ib < 91; ib++) { // 90 bins ] 10.0 , 100.0 ]
            lNDesiredBoundaries++;
            lDesiredBoundaries[lNDesiredBoundaries] = lDesiredBoundaries[lNDesiredBoundaries-1] + 1.0;
            // cout << "loop 3: " << lDesiredBoundaries[lNDesiredBoundaries] << endl;
        }
    }
    else { // default binning
        for( Int_t ib = 1; ib < 101; ib++) { // 100 bins  ] 0.0 , 0.1 ]
            lNDesiredBoundaries++;
            lDesiredBoundaries[lNDesiredBoundaries] = lDesiredBoundaries[lNDesiredBoundaries-1] + 0.001;
            // cout << "loop 1: " << lDesiredBoundaries[lNDesiredBoundaries] << endl;
        }
        for( Int_t ib = 1; ib < 91; ib++) { // 90 bins  ] 0.1 , 1.0 ]
            lNDesiredBoundaries++;
            lDesiredBoundaries[lNDesiredBoundaries] = lDesiredBoundaries[lNDesiredBoundaries-1] + 0.01;
            // cout << "loop 2: " << lDesiredBoundaries[lNDesiredBoundaries] << endl;
        }
        for( Int_t ib = 1; ib < 91; ib++) { // 90 bins ] 1.0 , 10. ]
            lNDesiredBoundaries++;
            lDesiredBoundaries[lNDesiredBoundaries] = lDesiredBoundaries[lNDesiredBoundaries-1] + 0.1;
            // cout << "loop 3: " << lDesiredBoundaries[lNDesiredBoundaries] << endl;
        }
        for( Int_t ib = 1; ib < 91; ib++) { // 90 bins ] 10.0 , 100.0 ]
            lNDesiredBoundaries++;
            lDesiredBoundaries[lNDesiredBoundaries] = lDesiredBoundaries[lNDesiredBoundaries-1] + 1.0;
            // cout << "loop 4: " << lDesiredBoundaries[lNDesiredBoundaries] << endl;
        }
    }
    // for(int i=0; i<lNDesiredBoundaries; ++i)
    // cout << i << "\t[ " << lDesiredBoundaries[i] << " , " << lDesiredBoundaries[i+1] << " ]\t" << lDesiredBoundaries[i+1]-lDesiredBoundaries[i] << endl;

    
    TH1F* hTestPercentileMB = new TH1F("hTestPercentileMB", "", lNDesiredBoundaries, lDesiredBoundaries);
    TH1F* hTestPercentileHM = new TH1F("hTestPercentileHM", "", lNDesiredBoundaries, lDesiredBoundaries);
    TH1F* hTestPercentileHM_anchored = new TH1F("hTestPercentileHM_anchored", "", lNDesiredBoundaries, lDesiredBoundaries);
    
    //Input file here!
    TString fInputFileNameMB = Form("../MB/files/AnalysisResults_%d.root", lRun); 
    TString fInputFileNameHM = Form("../VHM/files/AnalysisResults_%d.root", lRun); 
    
    cout<<" Offline Calibration Test "<<endl;
    cout<<" * Input File MB.....: "<<fInputFileNameMB.Data()<<endl;
    cout<<" * Input File VHM....: "<<fInputFileNameHM.Data()<<endl;
    cout<<endl;
    cout<<" Event Selection Peformed: "<<endl;
    oadbMultSelection->GetEventCuts()->Print();
    cout<<endl;

    //Declare minimal set of variables in memory
    Float_t fAmplitude_V0A = 0;
    Float_t fAmplitude_V0C = 0;
    Float_t fAmplitude_V0Apartial = 0;
    Float_t fAmplitude_V0Cpartial = 0;
    Float_t fAmplitude_V0AEq = 0;
    Float_t fAmplitude_V0CEq = 0;
    //Event Selection Variables
    Bool_t fEvSel_IsNotPileupInMultBins      = kFALSE ;
    Bool_t fEvSel_Triggered                  = kFALSE ;
    Bool_t fEvSel_INELgtZERO                 = kFALSE ;
    Bool_t fEvSel_PassesTrackletVsCluster    = kFALSE ;
    Bool_t fEvSel_HasNoInconsistentVertices  = kFALSE ;
    Float_t fEvSel_VtxZ                      = 10.0 ;
    UInt_t fEvSel_TriggerMask                = 0;
    Int_t fRefMultEta8;
    Int_t fRunNumber;
    Int_t fnContributors; 
    
    
    // STEP 1: Basic I/O
    cout<<"(1) Opening File MB"<<endl;
    
    //Open File
    TFile *fInputMB = TFile::Open( fInputFileNameMB.Data(), "READ");
    if(!fInputMB){
        AliWarningF("File %s not found!", fInputFileNameMB.Data() );
        return kFALSE;
    }
    //Locate TTree object
    TTree* fTreeMB = (TTree*)fInputMB->FindObjectAny("fTreeEvent");
    if(!fTreeMB){
        AliWarning("fTreeEvent object not found!" );
        return kFALSE;
    }
    //SetBranchAddresses
    fTreeMB->SetBranchAddress("fAmplitude_V0A",&fAmplitude_V0A);
    fTreeMB->SetBranchAddress("fAmplitude_V0C",&fAmplitude_V0C);
    fTreeMB->SetBranchAddress("fAmplitude_V0Apartial",&fAmplitude_V0Apartial);
    fTreeMB->SetBranchAddress("fAmplitude_V0Cpartial",&fAmplitude_V0Cpartial);
    fTreeMB->SetBranchAddress("fAmplitude_V0AEq",&fAmplitude_V0AEq);
    fTreeMB->SetBranchAddress("fAmplitude_V0CEq",&fAmplitude_V0CEq);
    
    fTreeMB->SetBranchAddress("fEvSel_IsNotPileupInMultBins",&fEvSel_IsNotPileupInMultBins);
    fTreeMB->SetBranchAddress("fEvSel_PassesTrackletVsCluster",&fEvSel_PassesTrackletVsCluster);
    fTreeMB->SetBranchAddress("fEvSel_HasNoInconsistentVertices",&fEvSel_HasNoInconsistentVertices);
    fTreeMB->SetBranchAddress("fEvSel_Triggered",&fEvSel_Triggered);
    fTreeMB->SetBranchAddress("fEvSel_INELgtZERO",&fEvSel_INELgtZERO);
    fTreeMB->SetBranchAddress("fEvSel_VtxZ",&fEvSel_VtxZ);
    
    fTreeMB->SetBranchAddress("fRunNumber",&fRunNumber);
    fTreeMB->SetBranchAddress("fRefMultEta8",&fRefMultEta8);
    fTreeMB->SetBranchAddress("fnContributors",&fnContributors);
    
    Long64_t lNEv = fTreeMB->GetEntries();
    cout<<"(1) File opened, event count is "<<lNEv<<endl;
    
    Long_t lSelectedMB = 0;
    Long_t lAllMB = 0;
    
    for(Long64_t iEv = 0; iEv<lNEv; iEv++) {

        fTreeMB->GetEntry(iEv); //Look at next event
        if ( lRun != fRunNumber ) continue; //skip if not desired run
        lAllMB++;

        if ( !fEvSel_Triggered ) continue; 
        if ( !fEvSel_IsNotPileupInMultBins ) continue; 
        if ( !fEvSel_PassesTrackletVsCluster ) continue; 
        if ( !fEvSel_INELgtZERO ) continue ;
        if ( !fEvSel_HasNoInconsistentVertices ) continue ; 
        if ( TMath::Abs(fEvSel_VtxZ) > 10. ) continue; 
        lSelectedMB++;

        Float_t lRawValue = ((fAmplitude_V0A)+(fAmplitude_V0C));
        Float_t lPercentile = hCalibV0M->GetBinContent( hCalibV0M->FindBin ( lRawValue ) );
        if( TMath::Abs( fEvSel_VtxZ )<10 ) hTestPercentileMB->Fill( lPercentile );
    }
    
    cout<<"All     : "<<lAllMB<<endl;
    cout<<"Selected: "<<lSelectedMB<<endl;
    cout<<endl;

    // STEP 2: Basic I/O
    cout<<"(2) Opening File VHM"<<endl;
    
    //Open File
    TFile *fInputHM = TFile::Open( fInputFileNameHM.Data(), "READ");
    if(!fInputHM){
        AliWarningF("File %s not found!", fInputFileNameHM.Data() );
        return kFALSE;
    }
    //Locate TTree object
    TTree* fTreeHM = (TTree*)fInputHM->FindObjectAny("fTreeEvent");
    if(!fTreeHM){
        AliWarning("fTreeEvent object not found!" );
        return kFALSE;
    }
    //SetBranchAddresses
    fTreeHM->SetBranchAddress("fAmplitude_V0A",&fAmplitude_V0A);
    fTreeHM->SetBranchAddress("fAmplitude_V0C",&fAmplitude_V0C);
    fTreeHM->SetBranchAddress("fAmplitude_V0Apartial",&fAmplitude_V0Apartial);
    fTreeHM->SetBranchAddress("fAmplitude_V0Cpartial",&fAmplitude_V0Cpartial);
    fTreeHM->SetBranchAddress("fAmplitude_V0AEq",&fAmplitude_V0AEq);
    fTreeHM->SetBranchAddress("fAmplitude_V0CEq",&fAmplitude_V0CEq);
    
    fTreeHM->SetBranchAddress("fEvSel_IsNotPileupInMultBins",&fEvSel_IsNotPileupInMultBins);
    fTreeHM->SetBranchAddress("fEvSel_PassesTrackletVsCluster",&fEvSel_PassesTrackletVsCluster);
    fTreeHM->SetBranchAddress("fEvSel_HasNoInconsistentVertices",&fEvSel_HasNoInconsistentVertices);
    fTreeHM->SetBranchAddress("fEvSel_Triggered",&fEvSel_Triggered);
    fTreeHM->SetBranchAddress("fEvSel_INELgtZERO",&fEvSel_INELgtZERO);
    fTreeHM->SetBranchAddress("fEvSel_VtxZ",&fEvSel_VtxZ);
    
    fTreeHM->SetBranchAddress("fRunNumber",&fRunNumber);
    fTreeHM->SetBranchAddress("fRefMultEta8",&fRefMultEta8);
    fTreeHM->SetBranchAddress("fnContributors",&fnContributors);
    
    lNEv = fTreeHM->GetEntries();
    cout<<"(2) File opened, event count is "<<lNEv<<endl;
    
    Long_t lSelectedHM = 0;
    Long_t lAllHM = 0;
    
    for(Long64_t iEv = 0; iEv<lNEv; iEv++) {

        fTreeHM->GetEntry(iEv); //Look at next event
        if ( lRun != fRunNumber ) continue; //skip if not desired run
        lAllHM++;

        if ( !fEvSel_Triggered ) continue; 
        if ( !fEvSel_IsNotPileupInMultBins ) continue; 
        if ( !fEvSel_PassesTrackletVsCluster ) continue; 
        if ( !fEvSel_INELgtZERO ) continue ;
        if ( !fEvSel_HasNoInconsistentVertices ) continue ; 
        if ( TMath::Abs(fEvSel_VtxZ) > 10. ) continue; 
        lSelectedHM++;

        Float_t lRawValue = ((fAmplitude_V0A)+(fAmplitude_V0C));
        Float_t lPercentile = hCalibV0M->GetBinContent( hCalibV0M->FindBin ( lRawValue ) );
        if( TMath::Abs( fEvSel_VtxZ )<10 ) hTestPercentileHM->Fill( lPercentile );
        if( TMath::Abs( fEvSel_VtxZ )<10 && lPercentile<anchor_percentile) hTestPercentileHM_anchored->Fill( lPercentile );
    }
    
    cout<<"All     : "<<lAllHM<<endl;
    cout<<"Selected: "<<lSelectedHM<<endl;
    cout<<endl;
    
    gStyle->SetOptStat(0); 
    gStyle->SetGridColor(kGray);

    hTestPercentileMB->Sumw2(); 
    hTestPercentileHM->Sumw2(); 
    hTestPercentileHM_anchored->Sumw2(); 
    
    Double_t normFactorMB = 100./((Double_t)hTestPercentileMB->Integral(hTestPercentileMB->FindBin(0.+1.e-9), hTestPercentileMB->FindBin(100.-1.e-9)));
    Double_t normFactorHM = 0.05/((Double_t)hTestPercentileHM->Integral(hTestPercentileHM->FindBin(0.+1.e-9), hTestPercentileHM->FindBin(0.05-1.e-9)));
    hTestPercentileMB->Scale( normFactorMB, "width"); 
    hTestPercentileHM->Scale( normFactorHM, "width"); 
    hTestPercentileHM_anchored->Scale( normFactorHM, "width"); 
    
    hTestPercentileMB->SetMarkerColor(kBlue); 
    hTestPercentileMB->SetLineColor(kBlue); 
    hTestPercentileMB->SetLineWidth(2);

    hTestPercentileHM->SetMarkerColor(kGreen); 
    hTestPercentileHM->SetLineColor(kGreen); 
    hTestPercentileHM->SetLineWidth(2);
    
    hTestPercentileHM_anchored->SetMarkerColor(kBlack); 
    hTestPercentileHM_anchored->SetLineColor(kBlack); 
    hTestPercentileHM_anchored->SetFillColor(kBlack);
    hTestPercentileHM_anchored->SetFillStyle(3004);
    hTestPercentileHM_anchored->SetLineStyle(2);
    hTestPercentileHM_anchored->SetLineWidth(1);
    
    //TCanvas *c2 = new TCanvas( "c2", "", 800,600);
    TCanvas *c2 = new TCanvas( "c2", "", 10, 10, 800, 600);
    c2->SetTicks(1,1); 
    c2->SetLeftMargin(0.15); 
    c2->SetBottomMargin(0.13); 
    c2->SetRightMargin(0.03); 
    c2->SetTopMargin(0.03); 
    //c2->SetLogx();
    //c2->SetLogy();
    c2->SetGridx();
    c2->SetGridy();
    hTestPercentileMB -> Draw("hist");
    hTestPercentileMB -> GetYaxis()->SetTitle("Normalized Counts"); 
    hTestPercentileMB -> GetYaxis()->SetTitleSize(0.055); 
    hTestPercentileMB -> GetYaxis()->SetTitleOffset(1.2); 
    hTestPercentileMB -> GetXaxis()->SetTitle("V0M Percentile"); 
    hTestPercentileMB -> GetXaxis()->SetTitleSize(0.055); 
    hTestPercentileMB -> GetXaxis()->SetRangeUser(0., 0.2); 
    hTestPercentileMB -> GetYaxis()->SetRangeUser(0., 2.0); 
    //
    hTestPercentileHM->Draw("hist same");
    //
    hTestPercentileHM_anchored->Draw("hist same");
    
    //hPos->SetMarkerColor(kGreen+1); 
    //hPos->SetLineColor(kGreen+1); 
    //hNeg->SetMarkerColor(kRed+1); 
    //hNeg->SetLineColor(kRed+1); 
    
    //hPos->Draw("same"); 
    //hNeg->Draw("same"); 
    
    TLegend* leg = new TLegend(0.55, 0.70, 0.9, 0.9);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.04); 

   leg->SetHeader(Form("Run: %d", lRun));
   leg->AddEntry(hTestPercentileMB, Form("%s MB", lPeriodName.Data()), "lp");
   leg->AddEntry(hTestPercentileHM, Form("%s VHM", lPeriodName.Data()), "lp");
   leg->AddEntry(hTestPercentileHM_anchored, Form("%s VHM (anchored)", lPeriodName.Data()), "lp");
   leg->Draw();
    
   c2->SaveAs(Form("checkStitching_%d.pdf", lRun)); 

   TCanvas *c3 = new TCanvas("c3", "", 820, 10, 800, 600);
   c3->SetTicks(1,1); 
   c3->SetLeftMargin(0.15); 
   c3->SetBottomMargin(0.13); 
   c3->SetRightMargin(0.03); 
   c3->SetTopMargin(0.03); 
   c3->SetGridx();
   c3->SetGridy();
   TH1D* h1 = (TH1D*)hTestPercentileMB->Clone("h1_copy");
   h1->GetXaxis()->SetRangeUser(0., 100.);
   h1->Draw("hist");
   TH1D* h2 = (TH1D*)hTestPercentileHM->Clone("h2_copy");
   h2->Draw("hist same");
   leg->Draw();

    
}

    

    
    
