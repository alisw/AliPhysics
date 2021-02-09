#if !defined (__CINT__) || (defined(__MAKECINT__))
#include <iostream>
#include "TROOT.h"
#include "AliMultSelectionCuts.h"
#include "AliOADBMultSelection.h"
#include "TH1F.h"
#include "TFile.h"
#include "AliOADBContainer.h"
#endif


void TestWrite() ;
void TestRead(Int_t run) ;
void TestUpdate(Int_t run) ;

void TestOADBMultSelPP(const Char_t* inputDir, TString lPeriodName = "LHC18f", Int_t run = 287784, const Char_t* chunkName="") {
    //Code to test OADB in Offline Mode
    //Further developments needed!
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
    
    TFile * f = new TFile (Form("OADB-%s-MB.root", lPeriodName.Data()));
    //TFile * f = new TFile (Form("/wrk15a/rderradi/alice/AliPhysics/OADB/COMMON/MULTIPLICITY/data/OADB-%s.root", lPeriodName.Data()));
    AliOADBContainer * oadbContMS = (AliOADBContainer*) f->Get("MultSel");
    AliOADBMultSelection * oadbMultSelection = (AliOADBMultSelection* )oadbContMS->GetObject(run, "Default");
    
    oadbMultSelection->GetEventCuts()->Print();
    oadbMultSelection->GetMultSelection()->PrintInfo();
    
    //Attempt to open object for test run
    Long_t lTestRun = run;
    cout<<"Load OADB..."<<endl;
    AliOADBMultSelection * oadbMultSelection = (AliOADBMultSelection* )oadbContMS->GetObject(lTestRun, "Default");
    cout<<"OADB loaded for run "<<lTestRun<<endl;
    

    //FIXME: Substitute with loop over several estimators, as needed
    TH1F *htesting = oadbMultSelection->GetCalibHisto(0);
    
    AliMultSelection *fSelection = (AliMultSelection*) oadbMultSelection->GetMultSelection();
    
    fSelection->PrintInfo();
        
    TCanvas *c1 = new TCanvas("c1", "",800,600);
    htesting->Draw(); 
    
    c1->SaveAs(Form("temp/oadbTests/htesting-%s-%d.png", lPeriodName.Data(), run));
    
    //Input file here!
    //TString fInputFileName = "AnalysisTest.root" ;
    TString fInputFileName = Form("%s/AnalysisResults_%d.root", inputDir, run);
    if(chunkName[0]!='\0') fInputFileName = Form("%s/AnalysisResults_%s.root", inputDir, chunkName);
    
    cout<<" Offline Calibration Test "<<endl;
    cout<<" * Input File.....: "<<fInputFileName.Data()<<endl;
    cout<<endl;
    cout<<" Event Selection Peformed: "<<endl;
    oadbMultSelection->GetEventCuts()->Print();
    cout<<endl;
    
    // STEP 1: Basic I/O
    cout<<"(1) Opening File"<<endl;
    
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
    
    //SetBranchAddresses
    fTree->SetBranchAddress("fAmplitude_V0A",&fAmplitude_V0A);
    fTree->SetBranchAddress("fAmplitude_V0C",&fAmplitude_V0C);
    fTree->SetBranchAddress("fAmplitude_V0Apartial",&fAmplitude_V0Apartial);
    fTree->SetBranchAddress("fAmplitude_V0Cpartial",&fAmplitude_V0Cpartial);
    fTree->SetBranchAddress("fAmplitude_V0AEq",&fAmplitude_V0AEq);
    fTree->SetBranchAddress("fAmplitude_V0CEq",&fAmplitude_V0CEq);
    
    fTree->SetBranchAddress("fEvSel_IsNotPileupInMultBins",&fEvSel_IsNotPileupInMultBins);
    fTree->SetBranchAddress("fEvSel_PassesTrackletVsCluster",&fEvSel_PassesTrackletVsCluster);
    fTree->SetBranchAddress("fEvSel_HasNoInconsistentVertices",&fEvSel_HasNoInconsistentVertices);
    fTree->SetBranchAddress("fEvSel_Triggered",&fEvSel_Triggered);
    fTree->SetBranchAddress("fEvSel_INELgtZERO",&fEvSel_INELgtZERO);
    fTree->SetBranchAddress("fEvSel_VtxZ",&fEvSel_VtxZ);
    fTree->SetBranchAddress("fEvSel_TriggerMask",&fEvSel_TriggerMask);
    
    fTree->SetBranchAddress("fRunNumber",&fRunNumber);
    fTree->SetBranchAddress("fRefMultEta8",&fRefMultEta8);
    fTree->SetBranchAddress("fnContributors",&fnContributors);
    
    Long64_t lNEv = fTree->GetEntries();
    cout<<"(1) File opened, event count is "<<lNEv<<endl;
    
    Float_t lDeterminedPercentile = -100;
        
    TH1F* hTestPercentile = new TH1F("hTestPercentile", "",100,0,100);
    TH1F* hPos = new TH1F("hPos", "",100,0,100);
    TH1F* hNeg = new TH1F("hNeg", "",100,0,100);
    
    Int_t lRun = run; //Select by hand
    Long_t lSelected = 0;
    Long_t lAll = 0;
    
    for(Long64_t iEv = 0; iEv<lNEv; iEv++) {
        //for(Long64_t iEv = 0; iEv<10000; iEv++) {

        fTree->GetEntry(iEv); //Look at next event
        if ( lRun != fRunNumber ) continue; //skip if not desired run
        lAll ++;

        if ( !fEvSel_Triggered ) continue; 
        if ( !fEvSel_IsNotPileupInMultBins ) continue; 
        if ( !fEvSel_PassesTrackletVsCluster ) continue; 
        if ( !fEvSel_INELgtZERO ) continue ;
        if ( !fEvSel_HasNoInconsistentVertices ) continue ; 
        if ( TMath::Abs(fEvSel_VtxZ) > 10. ) continue; 
        if ( !(fEvSel_TriggerMask & (AliVEvent::kINT7 | AliVEvent::kINT7inMUON)) ) continue;


        lSelected ++;

        Float_t lRawValue = ((fAmplitude_V0A)+(fAmplitude_V0C));

        lDeterminedPercentile = htesting-> GetBinContent( htesting ->FindBin ( lRawValue ) );

        if( TMath::Abs( fEvSel_VtxZ )<10 ) hTestPercentile -> Fill ( lDeterminedPercentile );
        if( 7 < fEvSel_VtxZ && fEvSel_VtxZ < 18 ) hPos -> Fill ( lDeterminedPercentile );
        if( -18 < fEvSel_VtxZ && fEvSel_VtxZ < -7 ) hNeg -> Fill ( lDeterminedPercentile );
    }
    
    cout<<"All     : "<<lAll<<endl;
    cout<<"Selected: "<<lSelected<<endl;
        
    gStyle->SetOptStat(0); 

    hTestPercentile->Sumw2(); 
    hNeg->Sumw2(); 
    hPos->Sumw2(); 
    
    hTestPercentile->Scale( 1./hTestPercentile->GetEntries()); 
    hNeg->Scale( 1./hNeg->GetEntries()); 
    hPos->Scale( 1./hPos->GetEntries()); 
    
    hTestPercentile->SetMarkerColor(kBlack); 
    hTestPercentile->SetLineColor(kBlack); 
    
    TCanvas *c2 = new TCanvas( "c2", "", 800,600);
    c2->SetTicks(1,1); 
    c2->SetLeftMargin(0.15); 
    c2->SetBottomMargin(0.13); 
    c2->SetRightMargin(0.03); 
    c2->SetTopMargin(0.03); 
    hTestPercentile -> Draw();
    hTestPercentile -> GetYaxis()->SetRangeUser(0.00850,0.013);
    hTestPercentile -> GetYaxis()->SetTitle("Fraction of Sample"); 
    hTestPercentile -> GetYaxis()->SetTitleSize(0.055); 
    hTestPercentile -> GetYaxis()->SetTitleOffset(1.2); 
    hTestPercentile -> GetXaxis()->SetTitle("V0M Percentile"); 
    hTestPercentile -> GetXaxis()->SetTitleSize(0.055); 
    
    hPos->SetMarkerColor(kGreen+1); 
    hPos->SetLineColor(kGreen+1); 
    hNeg->SetMarkerColor(kRed+1); 
    hNeg->SetLineColor(kRed+1); 
    
    hPos->Draw("same"); 
    hNeg->Draw("same"); 
    
    TLegend* leg = new TLegend(0.448, 0.7944, 0.928, 0.909);
    leg->SetTextSize(0.03); 

   leg->AddEntry(hTestPercentile, "Nominal Vertex, |z| < 10 cm", "lp");
   leg->AddEntry(hNeg, "Displaced Vertex, -18 cm < z < -7 cm", "lp");
   leg->AddEntry(hPos, "Displaced Vertex, 7 cm < z < 18 cm", "lp");
   leg->Draw();
 
   c2->SaveAs(Form("temp/oadbTests/percentDistrib-%s-%d.png", lPeriodName.Data(), run));
   //c2->SaveAs("testnewcorrection.pdf"); 
    
}
