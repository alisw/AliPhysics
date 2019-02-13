#if !defined (__CINT__) || (defined(__MAKECINT__))
#include <iostream>
#include "TROOT.h"
#include "AliMultSelectionCuts.h"
#include "AliOADBMultSelection.h"
#include "TH1F.h"
#include "TFile.h"
#include "AliOADBContainer.h"
#endif

// global pointer to calib histogram
TH1D *hCalib = 0x0;

void TestStitchedOADB(const Char_t* inputDir,
                      TString lPeriodName = "LHC18f",
                      Int_t          lRun = 287071 ) {
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
    //TH1F *hCalibV0M = oadbMultSelection->GetCalibHisto(0);
    hCalib = (TH1D*)oadbMultSelection->GetCalibHisto(0);
    
    AliMultSelection *fSelection = (AliMultSelection*) oadbMultSelection->GetMultSelection();
    
    fSelection->PrintInfo();
    
    TCanvas *c1 = new TCanvas("c1", "",800,600);
    hCalib->Draw(); 
    
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

    
    //Input file here!
    TString fInputFileName = Form("%s/AnalysisResults_%d.root", inputDir, lRun); 
    
    cout<<" Offline Calibration Test "<<endl;
    cout<<" * Input File .....: "<<fInputFileName.Data()<<endl;
    cout<<endl;
    cout<<" Event Selection Peformed: "<<endl;
    oadbMultSelection->GetEventCuts()->Print();
    cout<<endl;

    // MB part ********************************************************

    // STEP 1: Basic I/O
    cout<<"(1) Opening File MB"<<endl;
    
    //Open File
    TFile *fInput = TFile::Open( fInputFileName.Data(), "READ");
    if(!fInput){
        AliWarningF("File %s not found!", fInputFileName.Data() );
        return kFALSE;
    }
    //Locate TTree object
    TTree* fTree = (TTree*)fInput->FindObjectAny("fTreeEvent");
    if(!fTree){
        AliWarning("fTreeEvent object not found!" );
        return kFALSE;
    }

    TH1D* hTestPercentileMB = new TH1D("hTestPercentileMB", "", lNDesiredBoundaries, lDesiredBoundaries);

    Long64_t lNEv = fTree->GetEntries();
    cout<<"(1) File opened, event count is "<<lNEv<<endl;
    
    Long64_t lSelectedMB = fTree->Draw("get_percentile(fAmplitude_V0A+fAmplitude_V0C)>>hTestPercentileMB", 
                                         Form("fRunNumber==%d && fEvSel_Triggered && fEvSel_IsNotPileupInMultBins && fEvSel_PassesTrackletVsCluster && fEvSel_INELgtZERO && fEvSel_HasNoInconsistentVertices && TMath::Abs(fEvSel_VtxZ)<=10.0 && isSelectedMB(fEvSel_TriggerMask)", lRun),
                                         "goff");
    Long64_t lAllMB = lNEv;
    
    cout<<"All     : "<<lAllMB<<endl;
    cout<<"Selected: "<<lSelectedMB<<endl;
    cout<<endl;

    // VHM part *******************************************************
    // STEP 2: Basic I/O
    cout<<"(2) Opening File VHM"<<endl;
    
    //Open File
    TH1D* hTestPercentileHM = new TH1D("hTestPercentileHM", "", lNDesiredBoundaries, lDesiredBoundaries);
    TH1D* hTestPercentileHM_anchored = new TH1D("hTestPercentileHM_anchored", "", lNDesiredBoundaries, lDesiredBoundaries);
    
    lNEv = fTree->GetEntries();
    cout<<"(2) File opened, event count is "<<lNEv<<endl;
    
    Long64_t lSelectedHM = fTree->Draw(Form("get_percentile(fAmplitude_V0A+fAmplitude_V0C)>>hTestPercentileHM"), 
                                         Form("fRunNumber==%d && fEvSel_Triggered && fEvSel_IsNotPileupInMultBins && fEvSel_PassesTrackletVsCluster && fEvSel_INELgtZERO && fEvSel_HasNoInconsistentVertices && TMath::Abs(fEvSel_VtxZ)<=10.0 && isSelectedHM(fEvSel_TriggerMask)", lRun),
                                         "goff");
    Long64_t lSelectedHM_anchored = fTree->Draw(Form("get_percentile(fAmplitude_V0A+fAmplitude_V0C)>>hTestPercentileHM_anchored"), 
                                                Form("get_percentile(fAmplitude_V0A+fAmplitude_V0C)<%lf && fRunNumber==%d && fEvSel_Triggered && fEvSel_IsNotPileupInMultBins && fEvSel_PassesTrackletVsCluster && fEvSel_INELgtZERO && fEvSel_HasNoInconsistentVertices && TMath::Abs(fEvSel_VtxZ)<=10.0 && isSelectedHM(fEvSel_TriggerMask)", anchor_percentile, lRun),
                                                  "goff");
    Long64_t lAllHM = lNEv;
    
    cout<<"All     : "<<lAllHM<<endl;
    cout<<"Selected: "<<lSelectedHM<<" ("<<lSelectedHM_anchored<<" anchored)"<<endl;
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
    
   c2->SaveAs(Form("temp/checkStitching/checkStitching_%d.png", lRun)); 

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
   c3->SaveAs(Form("temp/checkStitching/checkStitchingZoomOut_%d.png", lRun)); 
}

// function to get V0M percentile
Double_t get_percentile(Double_t amplitude) 
{
    
    Double_t percentile = hCalib->GetBinContent( hCalib->FindBin( amplitude ) );

    return percentile;

}

Bool_t isSelectedMB(AliVEvent::EOfflineTriggerTypes trgMask) {
   
   Bool_t isSel = kFALSE;
   isSel = trgMask & AliVEvent::kINT7;
   
   return isSel;
}

Bool_t isSelectedHM(AliVEvent::EOfflineTriggerTypes trgMask) {
   
   Bool_t isSel = kFALSE;
   isSel = trgMask & AliVEvent::kHighMultV0;
   
   return isSel;
}
