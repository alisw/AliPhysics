#if !defined(__CINT__) || defined(__MAKECINT__)
#include "TString.h"
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TList.h"
#include "TChain.h"
#include "TObject.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TSystem.h"
#include "Riostream.h"

#include "AliITSgeomTGeo.h"
#include "AliLog.h"
#include "AliRawReaderRoot.h"
#include "AliITSRawStreamSSD.h"
#include "AliCDBManager.h"
#include "AliCDBEntry.h"
#endif

/*  $Id:  $    */

//============================================================//
static const Int_t fgkNumOfLDCs = 8;      //number of SSD LDCs
static const Int_t fgkNumOfDDLs = 16;      //number of SSD DDLs
static const Int_t fgkSSDMODULES = 1698;      //total number of SSD modules
static const Int_t fgkSSDLADDERSLAYER5 = 34; //ladders on layer 5
static const Int_t fgkSSDLADDERSLAYER6 = 38; //ladders on layer 6
static const Int_t fgkSSDMODULESPERLADDERLAYER5 = 22; //modules per ladder - layer 5                           
static const Int_t fgkSSDMODULESPERLADDERLAYER6 = 25; //modules per ladder - layer 6                          
static const Int_t fgkSSDMODULESLAYER5 = 748; //total number of SSD modules - layer5                           
static const Int_t fgkSSDMODULESLAYER6 = 950; //total number of SSD modules - layer6                         
static const Int_t fgkNumberOfPSideStrips = 768; //number of P-side strips
//============================================================//

TList *initCM();
void makeCM(const char* filename, Int_t nEvents, TList *list);

//__________________________________________________________//
void readSSDCommonMode(const char* filename = "raw.root",
		       Int_t nEvents = -1) {
  //Reads the CM pseudo-channels and produces the CM map for both layers 
  //and for p and n-side.
  gStyle->SetPalette(1,0);

  TList *list = initCM();
  //list->ls();
  Printf("CM histograms: %d",list->GetEntries());
  makeCM(filename,nEvents,list);
} 
 
//__________________________________________________________//
TList *initCM() {
  //Initializes the histograms and returns the TList object
  TList *list = new TList();

  Int_t gLayer = 0,gLadder = 0, gModule = 0;
  Int_t gHistCounter = 0;
  TString gTitle;
  TH1F *gHistSSDCMModule[2*fgkSSDMODULES]; 
  for(Int_t iModule = 500; iModule < fgkSSDMODULES + 500; iModule++) {
    AliITSgeomTGeo::GetModuleId(iModule,gLayer,gLadder,gModule);
    gTitle = "SSD_CM_PSide_Layer"; gTitle += gLayer;
    gTitle += "_Ladder"; gTitle += gLadder;
    gTitle += "_Module"; gTitle += gModule;
    gHistSSDCMModule[gHistCounter] = new TH1F(gTitle.Data(),gTitle.Data(),
					      100,-50.,50.);
    gHistSSDCMModule[gHistCounter]->GetXaxis()->SetTitle("CM");
    list->Add(gHistSSDCMModule[gHistCounter]);
    gHistCounter += 1;
  }
  for(Int_t iModule = 500; iModule < fgkSSDMODULES + 500; iModule++) {
    AliITSgeomTGeo::GetModuleId(iModule,gLayer,gLadder,gModule);
    gTitle = "SSD_CM_NSide_Layer"; gTitle += gLayer;
    gTitle += "_Ladder"; gTitle += gLadder;
    gTitle += "_Module"; gTitle += gModule;
    gHistSSDCMModule[gHistCounter] = new TH1F(gTitle.Data(),gTitle.Data(),
					      100,-50.,50.);
    gHistSSDCMModule[gHistCounter]->GetXaxis()->SetTitle("CM");
    list->Add(gHistSSDCMModule[gHistCounter]);
    gHistCounter += 1;
  }

  return list;
}

//__________________________________________________________//
void makeCM(const char* filename, Int_t nEvents, TList *list) {
  //Function to read the CM values
  Int_t gStripNumber = 0;
  Int_t gLayer = 0,gLadder = 0, gModule = 0;
  
  //==================================================//
  AliCDBManager *fCDBManager = AliCDBManager::Instance();
  fCDBManager->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  Int_t runNumber = atoi(gSystem->Getenv("DATE_RUN_NUMBER"));
  if(!runNumber) 
    Printf("DATE_RUN_NUMBER not defined!!!\n");
  
  fCDBManager->SetRun(runNumber);
  AliCDBEntry *geomGRP = fCDBManager->Get("GRP/Geometry/Data");
  if(!geomGRP) cout<<"GRP geometry not found!!!"<<endl;;
  //==================================================//
    
  //==================================================//
  TH2F *fHistPSideCMMapLayer5 = new TH2F("fHistPSideCMMapLayer5",
					 "Layer 5 - P side;N_{module};N_{ladder}",
					 22,1,23,
					 34,500,534);
  fHistPSideCMMapLayer5->GetXaxis()->SetTitleColor(1);
  fHistPSideCMMapLayer5->GetZaxis()->SetRangeUser(0.,20.);
  fHistPSideCMMapLayer5->SetStats(kFALSE);
  fHistPSideCMMapLayer5->GetYaxis()->SetTitleOffset(1.8);
  fHistPSideCMMapLayer5->GetXaxis()->SetNdivisions(22);
  fHistPSideCMMapLayer5->GetYaxis()->SetNdivisions(34);
  fHistPSideCMMapLayer5->GetXaxis()->SetLabelSize(0.03);
  fHistPSideCMMapLayer5->GetYaxis()->SetLabelSize(0.03);
  fHistPSideCMMapLayer5->GetZaxis()->SetTitleOffset(1.6);
  fHistPSideCMMapLayer5->GetZaxis()->SetTitle("RMS(CM) (p-side)");

  TH2F *fHistNSideCMMapLayer5 = new TH2F("fHistNSideCMMapLayer5",
					 "Layer 5 - N side;N_{module};N_{ladder}",
					 22,1,23,
					 34,500,534);
  fHistNSideCMMapLayer5->GetXaxis()->SetTitleColor(1);
  fHistNSideCMMapLayer5->GetZaxis()->SetRangeUser(0.,20.);
  fHistNSideCMMapLayer5->SetStats(kFALSE);
  fHistNSideCMMapLayer5->GetYaxis()->SetTitleOffset(1.8);
  fHistNSideCMMapLayer5->GetXaxis()->SetNdivisions(22);
  fHistNSideCMMapLayer5->GetYaxis()->SetNdivisions(34);
  fHistNSideCMMapLayer5->GetXaxis()->SetLabelSize(0.03);
  fHistNSideCMMapLayer5->GetYaxis()->SetLabelSize(0.03);
  fHistNSideCMMapLayer5->GetZaxis()->SetTitleOffset(1.6);
  fHistNSideCMMapLayer5->GetZaxis()->SetTitle("RMS(CM) (n-side)");
  
  TH2F *fHistPSideCMMapLayer6 = new TH2F("fHistPSideCMMapLayer6",
					 "Layer 6 - P side;N_{module};N_{ladder}",
					 25,1,26,
					 38,600,638);
  fHistPSideCMMapLayer6->GetXaxis()->SetTitleColor(1);
  fHistPSideCMMapLayer6->GetZaxis()->SetRangeUser(0.,20.);
  fHistPSideCMMapLayer6->SetStats(kFALSE);
  fHistPSideCMMapLayer6->GetYaxis()->SetTitleOffset(1.8);
  fHistPSideCMMapLayer6->GetXaxis()->SetNdivisions(25);
  fHistPSideCMMapLayer6->GetYaxis()->SetNdivisions(38);
  fHistPSideCMMapLayer6->GetXaxis()->SetLabelSize(0.03);
  fHistPSideCMMapLayer6->GetYaxis()->SetLabelSize(0.03);
  fHistPSideCMMapLayer6->GetZaxis()->SetTitleOffset(1.6);
  fHistPSideCMMapLayer6->GetZaxis()->SetTitle("RMS(CM) (p-side)");

  TH2F *fHistNSideCMMapLayer6 = new TH2F("fHistNSideCMMapLayer6",
					 "Layer 6 - N side;N_{module};N_{ladder}",
					 25,1,26,
					 38,600,638);
  fHistNSideCMMapLayer6->GetXaxis()->SetTitleColor(1);
  fHistNSideCMMapLayer6->GetZaxis()->SetRangeUser(0.,20.);
  fHistNSideCMMapLayer6->SetStats(kFALSE);
  fHistNSideCMMapLayer6->GetYaxis()->SetTitleOffset(1.8);
  fHistNSideCMMapLayer6->GetXaxis()->SetNdivisions(25);
  fHistNSideCMMapLayer6->GetYaxis()->SetNdivisions(38);
  fHistNSideCMMapLayer6->GetXaxis()->SetLabelSize(0.03);
  fHistNSideCMMapLayer6->GetYaxis()->SetLabelSize(0.03);
  fHistNSideCMMapLayer6->GetZaxis()->SetTitleOffset(1.6);
  fHistNSideCMMapLayer6->GetZaxis()->SetTitle("RMS(CM) (n-side)");
  //==================================================//

  TChain *chain = new TChain("RAW");
  chain->Add(filename);
  Int_t nTotalEvents = chain->GetEntries();
  if(nEvents == -1) nEvents = nTotalEvents;
  
  AliRawReaderRoot *rawReader = new AliRawReaderRoot(filename);
  Int_t iEvent = 0;
  Int_t fSSDEvent = 0;
  for(iEvent = 0; iEvent < nEvents; iEvent++) {
    cout<<"Event: "<<iEvent+1<<"/"<<nEvents<<endl;
    rawReader->Select("ITSSSD",-1,-1);  
    rawReader->Reset(); 
    rawReader->NextEvent();   
    
    if(rawReader->GetType() != 7) continue;
    
    fSSDEvent += 1;
    AliITSRawStreamSSD gSSDStream(rawReader);    
    while (gSSDStream.Next()) {
      if(gSSDStream.GetModuleID() < 0) continue;
      AliITSgeomTGeo::GetModuleId(gSSDStream.GetModuleID(),gLayer,gLadder,gModule);
      //if(gSSDStream.GetModuleID() != 500) continue;
      //Printf("Module id: %d - Layer: %d - Ladder: %d - Module: %d",gSSDStream.GetModuleID(),gLayer,gLadder,gModule);

      if(gSSDStream.GetStrip() >= 0) continue;
      gStripNumber = (gSSDStream.GetSideFlag() == 0) ? gSSDStream.GetStrip() : -gSSDStream.GetStrip() + 2*fgkNumberOfPSideStrips;
      //Printf("Module id: %d - Strip: %d - strip number: %d - Signal: %lf",gSSDStream.GetModuleID(),gSSDStream.GetStrip(),gStripNumber,signal);
      if(TMath::Abs(gSSDStream.GetStrip()) < 7)
	((TH1*)list->At(gSSDStream.GetModuleID()-500))->Fill(gSSDStream.GetSignal());
      if(TMath::Abs(gSSDStream.GetStrip()) > 6)
	((TH1*)list->At(1698+gSSDStream.GetModuleID()-500))->Fill(gSSDStream.GetSignal());
      
    }//streamer loop
  }//event loop
  
  //compute the rms of the CM values
  Double_t rmsPsideCM = 0.0, rmsNsideCM = 0.0;
  for(Int_t i = 0; i < 1698; i++) {
    rmsPsideCM = 0.0; rmsNsideCM = 0.0;
    AliITSgeomTGeo::GetModuleId(i+500,gLayer,gLadder,gModule);
    //Printf("%s - %s",((TH1*)list->At(i))->GetName(),
    //((TH1*)list->At(1698+i))->GetName());
    rmsPsideCM = ((TH1*)list->At(i))->GetRMS();
    rmsNsideCM = ((TH1*)list->At(1698+i))->GetRMS();
    if(rmsPsideCM == 0) rmsPsideCM = 0.001;
    if(rmsNsideCM == 0) rmsNsideCM = 0.001;
    //Printf("rmsPside: %lf - rmsNside: %lf",rmsPsideCM,rmsNsideCM);
    if(gLayer == 5) {
      fHistPSideCMMapLayer5->SetBinContent(gModule,gLadder,rmsPsideCM);
      fHistNSideCMMapLayer5->SetBinContent(gModule,gLadder,rmsNsideCM);
    }
     if(gLayer == 6) {
      fHistPSideCMMapLayer6->SetBinContent(gModule,gLadder,rmsPsideCM);
      fHistNSideCMMapLayer6->SetBinContent(gModule,gLadder,rmsNsideCM);
    }
    }

  TFile *foutput = TFile::Open("SSD.CM.root","recreate");
  list->Write();
  fHistPSideCMMapLayer5->Write(); fHistNSideCMMapLayer5->Write();
  fHistPSideCMMapLayer6->Write(); fHistNSideCMMapLayer6->Write();
  foutput->Close();
}
  
//__________________________________________________________//
void drawSSDCM(const char* filename = "SSD.CM.root") {
  gStyle->SetPalette(1,0);

  TFile *f = TFile::Open(filename);  

  TCanvas *c1 = new TCanvas("c1","CM values",0,0,700,650);
  c1->SetFillColor(10); c1->SetHighLightColor(10); c1->Divide(2,2);
  c1->cd(1);
  TH1F *fHistPSideCMMapLayer5 = dynamic_cast<TH1F *>(f->Get("fHistPSideCMMapLayer5"));
  fHistPSideCMMapLayer5->Draw("colz");
  c1->cd(2);
  TH1F *fHistNSideCMMapLayer5 = dynamic_cast<TH1F *>(f->Get("fHistNSideCMMapLayer5"));
  fHistNSideCMMapLayer5->Draw("colz");
  c1->cd(3);
  TH1F *fHistPSideCMMapLayer6 = dynamic_cast<TH1F *>(f->Get("fHistPSideCMMapLayer6"));
  fHistPSideCMMapLayer6->Draw("colz");
  c1->cd(4);
  TH1F *fHistNSideCMMapLayer6 = dynamic_cast<TH1F *>(f->Get("fHistNSideCMMapLayer6"));
  fHistNSideCMMapLayer6->Draw("colz");
}
