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
static const Int_t fgkNumOfChips = 6;      //number of SSD chips per module per side
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
  TH1F *gHistSSDCMModule[2*fgkSSDMODULES][fgkNumOfChips]; 
  for(Int_t iModule = 500; iModule < fgkSSDMODULES + 500; iModule++) {
    AliITSgeomTGeo::GetModuleId(iModule,gLayer,gLadder,gModule);
    for(Int_t iChip = 0; iChip < fgkNumOfChips; iChip++) {
      gTitle = "SSD_CM_PSide_Layer"; gTitle += gLayer;
      gTitle += "_Ladder"; gTitle += gLadder;
      gTitle += "_Module"; gTitle += gModule;
      gTitle += "_Chip"; gTitle += iChip+1;
      gHistSSDCMModule[gHistCounter][iChip] = new TH1F(gTitle.Data(),
						       gTitle.Data(),
						       100,-50.,50.);
      gHistSSDCMModule[gHistCounter][iChip]->GetXaxis()->SetTitle("CM");
      list->Add(gHistSSDCMModule[gHistCounter][iChip]);
    }
    gHistCounter += 1;
  }
  for(Int_t iModule = 500; iModule < fgkSSDMODULES + 500; iModule++) {
    AliITSgeomTGeo::GetModuleId(iModule,gLayer,gLadder,gModule);
    for(Int_t iChip = 0; iChip < fgkNumOfChips; iChip++) {
      gTitle = "SSD_CM_NSide_Layer"; gTitle += gLayer;
      gTitle += "_Ladder"; gTitle += gLadder;
      gTitle += "_Module"; gTitle += gModule;
      gTitle += "_Chip"; gTitle += iChip+1;
      gHistSSDCMModule[gHistCounter][iChip] = new TH1F(gTitle.Data(),
						       gTitle.Data(),
						       100,-50.,50.);
    gHistSSDCMModule[gHistCounter][iChip]->GetXaxis()->SetTitle("CM");
    list->Add(gHistSSDCMModule[gHistCounter][iChip]);
    }
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
  TH2F *fHistPSideMeanCMMapLayer5 = new TH2F("fHistPSideMeanCMMapLayer5",
					 "Layer 5 - P side;N_{module};N_{ladder}",
					 22,1,23,
					 34,500,534);
  fHistPSideMeanCMMapLayer5->GetXaxis()->SetTitleColor(1);
  fHistPSideMeanCMMapLayer5->GetZaxis()->SetRangeUser(0.,20.);
  fHistPSideMeanCMMapLayer5->SetStats(kFALSE);
  fHistPSideMeanCMMapLayer5->GetYaxis()->SetTitleOffset(1.8);
  fHistPSideMeanCMMapLayer5->GetXaxis()->SetNdivisions(22);
  fHistPSideMeanCMMapLayer5->GetYaxis()->SetNdivisions(34);
  fHistPSideMeanCMMapLayer5->GetXaxis()->SetLabelSize(0.03);
  fHistPSideMeanCMMapLayer5->GetYaxis()->SetLabelSize(0.03);
  fHistPSideMeanCMMapLayer5->GetZaxis()->SetTitleOffset(1.6);
  fHistPSideMeanCMMapLayer5->GetZaxis()->SetTitle("RMS(CM) (p-side)");

  TH2F *fHistNSideMeanCMMapLayer5 = new TH2F("fHistNSideMeanCMMapLayer5",
					 "Layer 5 - N side;N_{module};N_{ladder}",
					 22,1,23,
					 34,500,534);
  fHistNSideMeanCMMapLayer5->GetXaxis()->SetTitleColor(1);
  fHistNSideMeanCMMapLayer5->GetZaxis()->SetRangeUser(0.,20.);
  fHistNSideMeanCMMapLayer5->SetStats(kFALSE);
  fHistNSideMeanCMMapLayer5->GetYaxis()->SetTitleOffset(1.8);
  fHistNSideMeanCMMapLayer5->GetXaxis()->SetNdivisions(22);
  fHistNSideMeanCMMapLayer5->GetYaxis()->SetNdivisions(34);
  fHistNSideMeanCMMapLayer5->GetXaxis()->SetLabelSize(0.03);
  fHistNSideMeanCMMapLayer5->GetYaxis()->SetLabelSize(0.03);
  fHistNSideMeanCMMapLayer5->GetZaxis()->SetTitleOffset(1.6);
  fHistNSideMeanCMMapLayer5->GetZaxis()->SetTitle("RMS(CM) (n-side)");
  
  TH2F *fHistPSideMeanCMMapLayer6 = new TH2F("fHistPSideMeanCMMapLayer6",
					 "Layer 6 - P side;N_{module};N_{ladder}",
					 25,1,26,
					 38,600,638);
  fHistPSideMeanCMMapLayer6->GetXaxis()->SetTitleColor(1);
  fHistPSideMeanCMMapLayer6->GetZaxis()->SetRangeUser(0.,20.);
  fHistPSideMeanCMMapLayer6->SetStats(kFALSE);
  fHistPSideMeanCMMapLayer6->GetYaxis()->SetTitleOffset(1.8);
  fHistPSideMeanCMMapLayer6->GetXaxis()->SetNdivisions(25);
  fHistPSideMeanCMMapLayer6->GetYaxis()->SetNdivisions(38);
  fHistPSideMeanCMMapLayer6->GetXaxis()->SetLabelSize(0.03);
  fHistPSideMeanCMMapLayer6->GetYaxis()->SetLabelSize(0.03);
  fHistPSideMeanCMMapLayer6->GetZaxis()->SetTitleOffset(1.6);
  fHistPSideMeanCMMapLayer6->GetZaxis()->SetTitle("RMS(CM) (p-side)");

  TH2F *fHistNSideMeanCMMapLayer6 = new TH2F("fHistNSideMeanCMMapLayer6",
					 "Layer 6 - N side;N_{module};N_{ladder}",
					 25,1,26,
					 38,600,638);
  fHistNSideMeanCMMapLayer6->GetXaxis()->SetTitleColor(1);
  fHistNSideMeanCMMapLayer6->GetZaxis()->SetRangeUser(0.,20.);
  fHistNSideMeanCMMapLayer6->SetStats(kFALSE);
  fHistNSideMeanCMMapLayer6->GetYaxis()->SetTitleOffset(1.8);
  fHistNSideMeanCMMapLayer6->GetXaxis()->SetNdivisions(25);
  fHistNSideMeanCMMapLayer6->GetYaxis()->SetNdivisions(38);
  fHistNSideMeanCMMapLayer6->GetXaxis()->SetLabelSize(0.03);
  fHistNSideMeanCMMapLayer6->GetYaxis()->SetLabelSize(0.03);
  fHistNSideMeanCMMapLayer6->GetZaxis()->SetTitleOffset(1.6);
  fHistNSideMeanCMMapLayer6->GetZaxis()->SetTitle("RMS(CM) (n-side)");

  //____________________________________________________________//
  TH2F *fHistPSideRMSCMMapLayer5 = new TH2F("fHistPSideRMSCMMapLayer5",
					 "Layer 5 - P side;N_{module};N_{ladder}",
					 22,1,23,
					 34,500,534);
  fHistPSideRMSCMMapLayer5->GetXaxis()->SetTitleColor(1);
  fHistPSideRMSCMMapLayer5->GetZaxis()->SetRangeUser(0.,20.);
  fHistPSideRMSCMMapLayer5->SetStats(kFALSE);
  fHistPSideRMSCMMapLayer5->GetYaxis()->SetTitleOffset(1.8);
  fHistPSideRMSCMMapLayer5->GetXaxis()->SetNdivisions(22);
  fHistPSideRMSCMMapLayer5->GetYaxis()->SetNdivisions(34);
  fHistPSideRMSCMMapLayer5->GetXaxis()->SetLabelSize(0.03);
  fHistPSideRMSCMMapLayer5->GetYaxis()->SetLabelSize(0.03);
  fHistPSideRMSCMMapLayer5->GetZaxis()->SetTitleOffset(1.6);
  fHistPSideRMSCMMapLayer5->GetZaxis()->SetTitle("RMS(CM) (p-side)");

  TH2F *fHistNSideRMSCMMapLayer5 = new TH2F("fHistNSideRMSCMMapLayer5",
					 "Layer 5 - N side;N_{module};N_{ladder}",
					 22,1,23,
					 34,500,534);
  fHistNSideRMSCMMapLayer5->GetXaxis()->SetTitleColor(1);
  fHistNSideRMSCMMapLayer5->GetZaxis()->SetRangeUser(0.,20.);
  fHistNSideRMSCMMapLayer5->SetStats(kFALSE);
  fHistNSideRMSCMMapLayer5->GetYaxis()->SetTitleOffset(1.8);
  fHistNSideRMSCMMapLayer5->GetXaxis()->SetNdivisions(22);
  fHistNSideRMSCMMapLayer5->GetYaxis()->SetNdivisions(34);
  fHistNSideRMSCMMapLayer5->GetXaxis()->SetLabelSize(0.03);
  fHistNSideRMSCMMapLayer5->GetYaxis()->SetLabelSize(0.03);
  fHistNSideRMSCMMapLayer5->GetZaxis()->SetTitleOffset(1.6);
  fHistNSideRMSCMMapLayer5->GetZaxis()->SetTitle("RMS(CM) (n-side)");
  
  TH2F *fHistPSideRMSCMMapLayer6 = new TH2F("fHistPSideRMSCMMapLayer6",
					 "Layer 6 - P side;N_{module};N_{ladder}",
					 25,1,26,
					 38,600,638);
  fHistPSideRMSCMMapLayer6->GetXaxis()->SetTitleColor(1);
  fHistPSideRMSCMMapLayer6->GetZaxis()->SetRangeUser(0.,20.);
  fHistPSideRMSCMMapLayer6->SetStats(kFALSE);
  fHistPSideRMSCMMapLayer6->GetYaxis()->SetTitleOffset(1.8);
  fHistPSideRMSCMMapLayer6->GetXaxis()->SetNdivisions(25);
  fHistPSideRMSCMMapLayer6->GetYaxis()->SetNdivisions(38);
  fHistPSideRMSCMMapLayer6->GetXaxis()->SetLabelSize(0.03);
  fHistPSideRMSCMMapLayer6->GetYaxis()->SetLabelSize(0.03);
  fHistPSideRMSCMMapLayer6->GetZaxis()->SetTitleOffset(1.6);
  fHistPSideRMSCMMapLayer6->GetZaxis()->SetTitle("RMS(CM) (p-side)");

  TH2F *fHistNSideRMSCMMapLayer6 = new TH2F("fHistNSideRMSCMMapLayer6",
					 "Layer 6 - N side;N_{module};N_{ladder}",
					 25,1,26,
					 38,600,638);
  fHistNSideRMSCMMapLayer6->GetXaxis()->SetTitleColor(1);
  fHistNSideRMSCMMapLayer6->GetZaxis()->SetRangeUser(0.,20.);
  fHistNSideRMSCMMapLayer6->SetStats(kFALSE);
  fHistNSideRMSCMMapLayer6->GetYaxis()->SetTitleOffset(1.8);
  fHistNSideRMSCMMapLayer6->GetXaxis()->SetNdivisions(25);
  fHistNSideRMSCMMapLayer6->GetYaxis()->SetNdivisions(38);
  fHistNSideRMSCMMapLayer6->GetXaxis()->SetLabelSize(0.03);
  fHistNSideRMSCMMapLayer6->GetYaxis()->SetLabelSize(0.03);
  fHistNSideRMSCMMapLayer6->GetZaxis()->SetTitleOffset(1.6);
  fHistNSideRMSCMMapLayer6->GetZaxis()->SetTitle("RMS(CM) (n-side)");
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
	((TH1*)list->At((gSSDStream.GetModuleID()-500)*fgkNumOfChips+gSSDStream.GetStrip()))->Fill(gSSDStream.GetSignal());
      if(TMath::Abs(gSSDStream.GetStrip()) > 6)
	((TH1*)list->At(fgkSSDMODULES*fgkNumOfChips+(gSSDStream.GetModuleID()-500)*fgkNumOfChips+gSSDStream.GetStrip()-fgkNumOfChips))->Fill(gSSDStream.GetSignal());
    }//streamer loop
  }//event loop
  
  //compute the rms of the CM values
  TH1F *gHistCMDummy = new TH1F("gHistCMDummy","",100,-50.,50.);
  Double_t meanPsideCM = 0.0, meanNsideCM = 0.0;
  Double_t rmsPsideCM = 0.0, rmsNsideCM = 0.0;
  for(Int_t iModule = 0; iModule < fgkSSDMODULES; iModule++) {
      meanPsideCM = 0.0; meanNsideCM = 0.0;
    rmsPsideCM = 0.0; rmsNsideCM = 0.0;
    AliITSgeomTGeo::GetModuleId(iModule+500,gLayer,gLadder,gModule);
  
    gHistCMDummy->Clear();
    for(Int_t iChip = 0; iChip < fgkNumOfChips; iChip++) {
      //cout<<"Name: "<<dynamic_cast<TH1*>(list->At(iModule*fgkNumOfChips+iChip))->GetName()<<endl;
      gHistCMDummy->Add((TH1*)list->At(iModule*fgkNumOfChips+iChip));
    } 
    meanPsideCM = TMath::Abs(gHistCMDummy->GetMean());
    rmsPsideCM = gHistCMDummy->GetRMS();

    gHistCMDummy->Clear();
    for(Int_t iChip = 0; iChip < fgkNumOfChips; iChip++) {
      //cout<<"Name: "<<dynamic_cast<TH1*>(list->At(fgkSSDMODULES*fgkNumOfChips+iModule*fgkNumOfChips+iChip))->GetName()<<endl;
      gHistCMDummy->Add((TH1*)list->At(fgkSSDMODULES*fgkNumOfChips+iModule*fgkNumOfChips+iChip));
    }
    meanNsideCM = TMath::Abs(gHistCMDummy->GetMean());
    rmsNsideCM = gHistCMDummy->GetRMS();

    if(meanPsideCM == 0) meanPsideCM = 0.001;
    if(meanNsideCM == 0) meanNsideCM = 0.001;
    if(rmsPsideCM == 0) rmsPsideCM = 0.001;
    if(rmsNsideCM == 0) rmsNsideCM = 0.001;
 
    if(gLayer == 5) {
      fHistPSideMeanCMMapLayer5->SetBinContent(gModule,gLadder,meanPsideCM);
      fHistNSideMeanCMMapLayer5->SetBinContent(gModule,gLadder,meanNsideCM);
      fHistPSideRMSCMMapLayer5->SetBinContent(gModule,gLadder,rmsPsideCM);
      fHistNSideRMSCMMapLayer5->SetBinContent(gModule,gLadder,rmsNsideCM);
    }
     if(gLayer == 6) {
      fHistPSideMeanCMMapLayer6->SetBinContent(gModule,gLadder,meanPsideCM);
      fHistNSideMeanCMMapLayer6->SetBinContent(gModule,gLadder,meanNsideCM);
      fHistPSideRMSCMMapLayer6->SetBinContent(gModule,gLadder,rmsPsideCM);
      fHistNSideRMSCMMapLayer6->SetBinContent(gModule,gLadder,rmsNsideCM);
    }
    }

  TFile *foutput = TFile::Open("SSD.CM.root","recreate");
  list->Write();
  fHistPSideMeanCMMapLayer5->Write(); fHistNSideMeanCMMapLayer5->Write();
  fHistPSideMeanCMMapLayer6->Write(); fHistNSideMeanCMMapLayer6->Write();
  fHistPSideRMSCMMapLayer5->Write(); fHistNSideRMSCMMapLayer5->Write();
  fHistPSideRMSCMMapLayer6->Write(); fHistNSideRMSCMMapLayer6->Write();
  foutput->Close();
}
  
//__________________________________________________________//
void drawSSDCM(const char* filename = "SSD.CM.root") {
  gStyle->SetPalette(1,0);

  TFile *f = TFile::Open(filename);  

  TCanvas *c1 = new TCanvas("c1","Mean of CM values",0,0,700,650);
  c1->SetFillColor(10); c1->SetHighLightColor(10); c1->Divide(2,2);
  c1->cd(1);
  TH2F *fHistPSideMeanCMMapLayer5 = dynamic_cast<TH2F *>(f->Get("fHistPSideMeanCMMapLayer5"));
  fHistPSideMeanCMMapLayer5->Draw("colz");
  c1->cd(2);
  TH2F *fHistNSideMeanCMMapLayer5 = dynamic_cast<TH2F *>(f->Get("fHistNSideMeanCMMapLayer5"));
  fHistNSideMeanCMMapLayer5->Draw("colz");
  c1->cd(3);
  TH2F *fHistPSideMeanCMMapLayer6 = dynamic_cast<TH2F *>(f->Get("fHistPSideMeanCMMapLayer6"));
  fHistPSideMeanCMMapLayer6->Draw("colz");
  c1->cd(4);
  TH2F *fHistNSideMeanCMMapLayer6 = dynamic_cast<TH2F *>(f->Get("fHistNSideMeanCMMapLayer6"));
  fHistNSideMeanCMMapLayer6->Draw("colz");

  TCanvas *c2 = new TCanvas("c2","RMS of CM values",100,100,700,650);
  c2->SetFillColor(10); c2->SetHighLightColor(10); c2->Divide(2,2);
  c2->cd(1);
  TH2F *fHistPSideRMSCMMapLayer5 = dynamic_cast<TH2F *>(f->Get("fHistPSideRMSCMMapLayer5"));
  fHistPSideRMSCMMapLayer5->Draw("colz");
  c2->cd(2);
  TH2F *fHistNSideRMSCMMapLayer5 = dynamic_cast<TH2F *>(f->Get("fHistNSideRMSCMMapLayer5"));
  fHistNSideRMSCMMapLayer5->Draw("colz");
  c2->cd(3);
  TH2F *fHistPSideRMSCMMapLayer6 = dynamic_cast<TH2F *>(f->Get("fHistPSideRMSCMMapLayer6"));
  fHistPSideRMSCMMapLayer6->Draw("colz");
  c2->cd(4);
  TH2F *fHistNSideRMSCMMapLayer6 = dynamic_cast<TH2F *>(f->Get("fHistNSideRMSCMMapLayer6"));
  fHistNSideRMSCMMapLayer6->Draw("colz");
}
