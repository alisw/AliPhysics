#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1.h"
#include "TH2.h"
#include "TBranch.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TList.h"
#include "TStopwatch.h"
#include "TProfile.h"
#include "TCanvas.h"

#include <iostream>
#include <map>

void EstimatorQA(Char_t* inputFileName="AnalysisResults.root", Char_t* outputFileName="EstimatorQAPlots.pdf" ){
   
    TCanvas* runCanvas = new TCanvas("runCanvas","quantities as a function of run number",800,600);
    
    TStopwatch* timer = new TStopwatch();
    gStyle->SetPalette(55);
    // Macro to produce QA plots of multiplicity estimators in pp
    // Checks on run dependence, z vertex dependence etc. of the quantities used in
    // in forming the estimators
    // Plots of flatness of the produced estimators

    
    // To plot versus run we need to discover the list of runs, use them as labels, also cut one-by-one on them?
    // Something like CountEvents.C example can do that
    // not currently possible - TTreeReader is only available in ROOT v6 onwards
    
// OLD file handling section
    
    // open the file
  //  TFile *f = TFile::Open(inputFileName);
 //   if (f == 0) {
        // if we cannot open the file, print an error message and return immediatly
  //      std::cout << "Error: cannot open input file " << inputFileName << " " << std::endl;
//return;
//}

    // Is this name fixed or varying?
//    TDirectoryFile* treeDir = gROOT->FindObject("PWGLF_StrVsMult");
 //   TDirectoryFile* treeDir = (TDirectoryFile*)gROOT->FindObject("MultSelection");
//if (treeDir == 0) {
        // if cannot retrieve the internal ROOT directory return
 //       std::cout << "Error: get internal directory " << "PWGLF_StrVsMult" << " " << std::endl;
//return;
//}
    
//treeDir->cd();
    
    //TTree* evTree = (TTree*)gROOT->FindObject("fTreeEvent");
 
    //
    
    TChain* evTree = new TChain("MultSelection/fTreeEvent","");
                                                // String depends on TDirectory name and TTree name
    evTree->Add(inputFileName);
    // Add some error handling here!
    
    
    // Changes above to use TChain to add the files in a directory
    // if we call the TChain evTree then everything below is unchanged
    
    timer->Stop();
    timer->Print();
    timer->Start(kFALSE);
    
    // Procedure for getting the run numbers:
    // Create histogram with bin size of exactly 1
    // Make the bin edges 0.5, eg new TH1D("name","title",m,x+0.5,x+m+0.5)
    // for m bins starting at run number m
    // Create automatic histogram  first to find the rough range

    
    evTree->Draw("fRunNumber>>runCoarse","","N");
    TH1D* hRunCoarse = (TH1D*)gROOT->FindObject("runCoarse");
    
    std::cout << "Entries" << hRunCoarse->GetEntries() << std::endl;
    
    Int_t nStartRun, nEndRun, nRunRange, nCoursebins;
    nStartRun = hRunCoarse->GetBinLowEdge(1);
    nEndRun = hRunCoarse->GetBinLowEdge(1+hRunCoarse->GetNbinsX());
    nRunRange = nEndRun-nStartRun+1;
    TH1D* hRunFine = new TH1D("runFine","Run number",nRunRange,nStartRun-0.5,nEndRun+0.5);

    // Here using TCut might work better
    evTree->Draw("fRunNumber>>runFine","fEvSel_INELgtZERO&&fEvSel_HasAtLeastSPDVertex&&fEvSel_IsNotPileupInMultBins&&(fEvSel_nContributors>0)","N");
    
    timer->Stop();
    timer->Print();
    timer->Start(kFALSE);
    
    Int_t numberOfRuns = 0;
    Int_t thisRun;
    TList* listOfRuns = new TList();
    TH1D* hTest = new TH1D("test","a test",1,0,1);
    hTest->SetBit(TH1::kCanRebin);
    
    // Add a map (it means macro has to be compiled .x EstimatorQA.C)
    std::map<int, int> runmap;
    std::map<int, int> neventsmap;


    for (Int_t i=1; i<=hRunFine->GetEntries(); i++) {
        if (hRunFine->GetBinContent(i)>0) {
            thisRun = hRunFine->GetBinLowEdge(i)+0.5;
            std::cout << "Run " << thisRun << " has " << hRunFine->GetBinContent(i) << " entries" << std::endl;
            numberOfRuns++;
            listOfRuns->Add(new TObjString(Form("%i",thisRun)));
            runmap.insert( std::pair<int,int>(thisRun,numberOfRuns));
            neventsmap.insert( std::pair<int,int>(thisRun,hRunFine->GetBinContent(i)));
            hTest->Fill(Form("%i",thisRun),hRunFine->GetBinContent(i));
        }
    }

    timer->Stop();
    timer->Print();
    timer->Start(kFALSE);
    
    listOfRuns->Print();
    hTest->LabelsDeflate();
    hTest->Draw();
    
    // All the variables needed for tree access go here
    Float_t v0A, v0C, v0Apartial, v0Cpartial, v0M, v0AEq, v0CEq, v0MEq, vtxZ;
    Int_t refMultEta5, refmultEta8, runNumber, evSel_nContributors;
    Bool_t evSel_HasAtLeastSPDVertex, evSel_INELgtZERO, evSel_IsNotPileupInMultBins;
    
    evTree->SetBranchAddress("fAmplitude_V0A",&v0A);
    evTree->SetBranchAddress("fAmplitude_V0C",&v0C);
    //evTree->SetBranchAddress("fAmplitude_V0M",&v0M); // not in David's latest Trees
    evTree->SetBranchAddress("fAmplitude_V0Apartial",&v0Apartial);
    evTree->SetBranchAddress("fAmplitude_V0Cpartial",&v0Cpartial);
    evTree->SetBranchAddress("fAmplitude_V0AEq",&v0AEq);
    evTree->SetBranchAddress("fAmplitude_V0CEq",&v0CEq);
    //evTree->SetBranchAddress("fAmplitude_V0MEq",&v0MEq); // not in David's latest Trees
    evTree->SetBranchAddress("fRefMultEta5",&refMultEta5);
    evTree->SetBranchAddress("fRefMultEta8",&refmultEta8);
    evTree->SetBranchAddress("fEvSel_VtxZ",&vtxZ);
    evTree->SetBranchAddress("fEvSel_nContributors",&evSel_nContributors);
    evTree->SetBranchAddress("fEvSel_HasAtLeastSPDVertex",&evSel_HasAtLeastSPDVertex);
    evTree->SetBranchAddress("fEvSel_INELgtZERO",&evSel_INELgtZERO);
    evTree->SetBranchAddress("fEvSel_IsNotPileupInMultBins",&evSel_IsNotPileupInMultBins);
    evTree->SetBranchAddress("fRunNumber",&runNumber);

    // Histograms vs run
    
    TH2D* hVtxZvsRun = new TH2D("vtxZvsRun", "run number vs Vertex Z;;Z_{VTX} (cm)",numberOfRuns,0,numberOfRuns,200,-20.,20.);
    hVtxZvsRun->Sumw2();

    //
    Float_t maxV0 = 100;
    TH2D* hV0AvsRun = new TH2D("v0AvsRun", "run number vs V0A;;V0A amplitude",numberOfRuns,0,numberOfRuns,maxV0,0,maxV0);
    hV0AvsRun->Sumw2();
    TH2D* hV0CvsRun = new TH2D("v0CvsRun", "run number vs V0C;;V0C amplitude",numberOfRuns,0,numberOfRuns,maxV0,0,maxV0);
    TH2D* hV0MvsRun = new TH2D("v0MvsRun", "run number vs V0M;;V0M amplitude",numberOfRuns,0,numberOfRuns,maxV0,0,maxV0);
    TH2D* hV0ApartialvsRun = new TH2D("v0ApartialvsRun", "run number vs V0Apartial;;V0Apartial amplitude",numberOfRuns,0,numberOfRuns,maxV0,0,maxV0);
    TH2D* hV0CpartialvsRun = new TH2D("v0CpartialvsRun", "run number vs V0Cpartial;;V0Cpartial amplitude",numberOfRuns,0,numberOfRuns,maxV0,0,maxV0);
    TH2D* hV0AEqvsRun = new TH2D("v0AEqvsRun", "run number vs V0AEq;;V0AEq amplitude",numberOfRuns,0,numberOfRuns,maxV0,0,maxV0);
    TH2D* hV0CEqvsRun = new TH2D("v0CEqvsRun", "run number vs V0CEq;;V0CEq amplitude",numberOfRuns,0,numberOfRuns,maxV0,0,maxV0);
    TH2D* hV0MEqvsRun = new TH2D("v0MEqvsRun", "run number vs V0MEq;;V0MEq amplitude",numberOfRuns,0,numberOfRuns,maxV0,0,maxV0);
    
    // 2-D histograms, integrated over all runs
    TH2D* hVtxZvsV0A = new TH2D("vtxZvsV0A","Vertex Z vs V0A;;V0A amplitude;Z_{VTX} (cm)",200.,-20,20,maxV0,0,maxV0);
    TH2D* hVtxZvsV0C = new TH2D("vtxZvsV0C","Vertex Z vs V0C;;V0C amplitude;Z_{VTX} (cm)",200.,-20,20,maxV0,0,maxV0);
    
    //Monitoring histograms
    TH1D* hEventsPerRun = new TH1D("eventsPerRun","Number of events in each run",numberOfRuns,0,numberOfRuns);
    TH1D* hGoodEventsPerRun = new TH1D("GoodventsPerRun","Number of selected events in each run",numberOfRuns,0,numberOfRuns);
    
    Int_t r=0;
    for (std::map<int,int>::iterator it=runmap.begin(); it!=runmap.end(); it++) {
        r++;
        hEventsPerRun->GetXaxis()->SetBinLabel(it->second,Form("%i",it->first));
        hEventsPerRun->SetBinContent(r,neventsmap[it->first]);
        std::cout << "Weight is: " << 1./neventsmap[it->first] << endl;
    }

    hEventsPerRun->Draw();
    
    // Weight given when filling so that distributions are normalised to the total events in each run
    // Also there are events which don't satisfy the event selection criteria so must ensure criteria are
    // the same in this loop and in the earlier code to determine events per run
    // Not trivial since one uses the tree variables directly and here the local variable from SetBranchAddress
    Double_t runWeight;
    Int_t runIndex;
    Bool_t evSelected;
    
    for (Int_t i=1; i<evTree->GetEntries(); i++) {
        evTree->GetEntry(i);
        //Get the index once only per entry
        // Perhaps even take advantage of likely ordering within the Tree and only get the index if
        // the run number has changed?
        evSelected = evSel_INELgtZERO && evSel_HasAtLeastSPDVertex && (evSel_nContributors > 0) && evSel_IsNotPileupInMultBins;
        runWeight = 1./neventsmap[runNumber];
        //runWeight= 1.; // For testing, could make a switch
        runIndex = runmap[runNumber]-1.;
        if (evSelected) {
            
            hVtxZvsRun->Fill(runIndex,vtxZ,runWeight);
            hV0AvsRun->Fill(runIndex,v0A,runWeight);
            hVtxZvsV0A->Fill(vtxZ,v0A);
            
            hV0CvsRun->Fill(runIndex,v0C,runWeight);
            hVtxZvsV0C->Fill(vtxZ,v0C);
        }
        // Put all the histogram filling in the same loop
    }

    // Set the labels by iterating over the map
    for (std::map<int,int>::iterator it=runmap.begin(); it!=runmap.end(); it++) {
        hVtxZvsRun->GetXaxis()->SetBinLabel(it->second,Form("%i",it->first));
        hV0AvsRun->GetXaxis()->SetBinLabel(it->second,Form("%i",it->first));
        hV0CvsRun->GetXaxis()->SetBinLabel(it->second,Form("%i",it->first));
        // All "per run" histograms to go here. If there are very many if might be good to store them
        // in some container (list, arrary etc.) and iterate over that
    }
    
    TProfile* hVtxZvsRun_prof = hVtxZvsRun->ProfileX();
    hVtxZvsRun_prof->SetLineColor(kWhite);
    hVtxZvsRun_prof->SetLineWidth(2);

    TProfile* hV0AvsRun_prof = hV0AvsRun->ProfileX();
    hV0AvsRun_prof->SetLineColor(kWhite);
    hV0AvsRun_prof->SetLineWidth(2);

    TProfile* hV0CvsRun_prof = hV0CvsRun->ProfileX();
    hV0AvsRun_prof->SetLineColor(kWhite);
    hV0AvsRun_prof->SetLineWidth(2);
    
    // Can add something to compare the profile histograms with results from FitSlicesY
    // Simple in the case of Zvertex when the functional form og the fit (Gaussian) is known
    // Not so obvious for V0 distributions
    
    // Draw some stuff and save into pdf
    // note root conventions for opening and closing the pdf to get multiple pages
    
    hVtxZvsRun->Draw("COLZ");
    hVtxZvsRun_prof->Draw("SAME");
    
    runCanvas->SaveAs(Form("%s(",outputFileName),"RECREATE"); // NB opening '(' in filename
    
    hV0AvsRun->Draw("COLZ");
    hV0AvsRun_prof->Draw("SAME");
    
    runCanvas->SaveAs(Form("%s",outputFileName));
    
    hVtxZvsV0A->Draw("COLZ");
    
    runCanvas->SaveAs(Form("%s",outputFileName));
    
    hVtxZvsV0C->Draw("COLZ");

    
    runCanvas->SaveAs(Form("%s)",outputFileName));  // NB closing '(' in filename
    

    timer->Stop();
    timer->Print();
    
}