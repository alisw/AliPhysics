#include <TCanvas.h>
#include <TProfile2D.h>
#include <TH2F.h>
#include <TFile.h>
#include <TTree.h>
#include <Riostream.h>
#include <TStyle.h>
#include <TROOT.h>

/// Macro for readins trees produced with countOfflineRun2.C, macro for EMCAL raw data analyisis
/// Instructions:
/// per each run: root
/// .L ReadTreeCountRun2.C+
/// PlotTreeDDL()

// config for LHC Run 2 EMCal data (2015-): 
// 20 SuperModules total
const int kNSMEMCal = 12;
const int kNSMDCal = 8;
const int kNSM = 20;

// 40 optical fibers/DDLs total; 2 per SuperModule
// 2 of the 1/3 SMs (A5 and C12) actually don't have 2 DDLs, but space is reserved for them anyhow (the missing DDLs are # 21 and 39)
const int kNDDLEMCal = 24; // 0..23
const int kNDDLDCal = 16; // 24..39
const int kNDDL = 40; 
// There are also 2 DDLs (and 2 LDCs) with trigger=STU data; this macro ignores the STU DDLs

// 12 LDCs with EMCal/DCal data
const int kNLDCEMCal = 8; // aldaqpc137..144 a.k.a. ldc-EMCAL-C-0..3 ldc-EMCAL-A-0..3 ; 3 DDLs per LDC except for the last one   
const int kNLDCDCal = 5; // aldaqpc146..150 a.k.a. ldc-DCAL-0..4; 3 DDLs per LDCs
const int kNLDC = 13; 


// how many channels per DDL and SuperModule
const int kNCHAN = 1152; // FEE channels, per DDL
const int kNCHANTRU = 96; // channels per TRU
const int kNCHANLEDREF = 48; // channels for LED Mon Ref. 

const int kNCHANDDL0 = kNCHAN + kNCHANTRU + kNCHANLEDREF;
const int kNCHANDDL1 = kNCHAN + 2*kNCHANTRU;
 
const int kNCHANTot = kNCHANDDL0 + kNCHANDDL1; // per SuperModule

TTree* ReadTree(TString filename, TString treename);

TTree* ReadTree(TString filename, TString treename){
   TFile *fin = new TFile(filename);
   if(!fin->IsOpen()){
      Printf("File %s not found", filename.Data());
      return 0x0;
   }
   TTree *tree = (TTree*)fin->Get(treename);
   return tree;
   
}

void PlotTreeDDL(TString filename = "output_tree_15000214142011.10.root", Int_t ampCut = 10){
   TString treename = "treeDDL";
   TTree *tree = ReadTree(filename, treename);
   if(!tree) {
      Printf("Tree %s not found, exit", treename.Data());
      return;
   }
   
   Int_t evno = 0;
   UInt_t type = 0;
   UInt_t period = 0;
   UInt_t orbit = 0;
   UInt_t bc = 0;
   
   // variables for treeDDL
   Int_t iDDL = 0;
   Int_t nChan = 0;
   Int_t nSamp = 0;
   Int_t hwaddress[kNCHANDDL1] = {0};
   Int_t column[kNCHANDDL1] = {0};
   Int_t row[kNCHANDDL1] = {0};
   Int_t caloflag[kNCHANDDL1] = {0};
   Int_t nbunches[kNCHANDDL1] = {0};
   Int_t nsamples[kNCHANDDL1] = {0};
   Int_t min[kNCHANDDL1] = {0};
   Int_t max[kNCHANDDL1] = {0};
   Int_t timeAtMax[kNCHANDDL1] = {0};

   tree->SetBranchAddress("evno", &evno);
   tree->SetBranchAddress("type", &type);
   tree->SetBranchAddress("period", &period);
   tree->SetBranchAddress("orbit", &orbit);
   tree->SetBranchAddress("bc", &bc);
   tree->SetBranchAddress("iDDL", &iDDL);
   tree->SetBranchAddress("nChan", &nChan);
   tree->SetBranchAddress("nSamp", &nSamp);
   tree->SetBranchAddress("hwaddress", &hwaddress);
   tree->SetBranchAddress("column", &column);
   tree->SetBranchAddress("row", &row);
   tree->SetBranchAddress("caloflag", &caloflag); 
   tree->SetBranchAddress("nbunches", &nbunches); 
   tree->SetBranchAddress("nsamples", &nsamples); 
   tree->SetBranchAddress("min", &min);
   tree->SetBranchAddress("max", &max);
   tree->SetBranchAddress("timeAtMax", &timeAtMax);
   
   TProfile2D *hAmp[kNSM];
   TProfile2D *hTime[kNSM];
      
   for(Int_t ism = 0;ism<kNSM; ism++){
      hAmp[ism] = new TProfile2D(Form("hAmp%d-Acut%d", ism, ampCut),Form("hAmp%d;column;row;Amp (max - min > %d)", ism, ampCut), 121, -0.5, 47.5, 24, -0.5, 23.5);
      hTime[ism] = new TProfile2D(Form("hTime%d-Acut%d", ism, ampCut),Form("hTime%d;column;row;time (max - min > %d)", ism, ampCut), 121, -0.5, 47.5, 24, -0.5, 23.5);
   }
   for(Int_t ie=0; ie< tree->GetEntries(); ie++){
      tree->GetEntry(ie);
      Int_t iSM = iDDL/2;
      for (Int_t ic=0; ic<nChan; ic++) {
      	 //Printf("Channel %d", ic);
      if ( (max[ic]-min[ic]) > ampCut ) {
      	 //Printf("Amp > %d", ampCut);
	if (caloflag[ic] == 1) {
	   //Printf("Filling SM %d %d %d, %d",iSM, column[ic], row[ic], max[ic]-min[ic] );
	  hAmp[iSM]->Fill(column[ic], row[ic], max[ic]-min[ic]); 
	  hTime[iSM]->Fill(column[ic], row[ic], timeAtMax[ic]); 
	}
      }
    }
   }
   
   TCanvas *cTree2DAmp = new TCanvas(Form("cTree2DAmp-Acut%d",ampCut), Form("Amp (> %d)",ampCut), 1000,1000);
   cTree2DAmp -> Divide (5,4); //might need adjustment
   TCanvas *cTree2DTim = new TCanvas(Form("cTree2DTim-Acut%d", ampCut), Form("Time (amp > %d)",ampCut), 1000,1000);
   cTree2DTim -> Divide (5,4); //might need adjustment
   for(Int_t ism = 0; ism < kNSM; ism++){
      cTree2DAmp -> cd(ism+1);
      hAmp[ism] ->Draw("colz");
      cTree2DTim -> cd(ism+1);
      hTime[ism] ->Draw("colz");
      
   }

   cTree2DAmp->SaveAs(Form("%s.png", cTree2DAmp->GetName()));
   cTree2DTim->SaveAs(Form("%s.png", cTree2DTim->GetName()));
}


