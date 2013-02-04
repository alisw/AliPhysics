#include "TROOT.h"
#include "TChain.h"
#include "TString.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "AliRawReaderChain.h"
#include "AliCaloRawStreamV3.h"
#include "AliCDBManager.h"

#include <fstream>
#include <iostream>

#include "AliTRUPedestalAnalysis.h"
#include "AliTRUPedestalOutput.h"
#include <TStyle.h>

void addFilesToChain(const TString rawFileList, TChain* chain);
void printHists(AliTRUPedestalOutput* output);
void saveResults(AliTRUPedestalOutput* output, TString saveToFile);



void truPedestalAnalysis(TString rawFileList = "files.txt", TString saveToFile = "truOutput.root")
{
  // Raw Chain Initialization
  TChain *chain = new TChain("RAW");
  addFilesToChain(rawFileList, chain );
  AliRawReaderChain* rawChain = new AliRawReaderChain(chain);
  rawChain->Reset();

  // PHOS Raw Stream Initialization
  AliCaloRawStreamV3* phosRawStream = new AliCaloRawStreamV3(rawChain,"PHOS");

  // Analysis object Initialization
  AliTRUPedestalAnalysis* anaObj = new AliTRUPedestalAnalysis();

  // Loop over events in Chain
  UInt_t runNumber = -1;
  Int_t event_count = 0;
  while (rawChain->NextEvent())
  { // Print out event number:
    std::cout << "\r" << "event: " << ++event_count
  	      << "/"<< rawChain->GetNumberOfEvents() << " " << std::flush;
    if( rawChain->GetRunNumber() != runNumber ){
      // if new event number, update OCDB
      runNumber = rawChain->GetRunNumber();
      AliCDBManager::Instance()->SetRun(runNumber);
      Printf("New run number, current run number is: %d", runNumber);
    }

    // Process Event using analysis Object
    anaObj->ProcessEvent(phosRawStream);
  }

  // Save output to file, in form of single entry in tree
  saveResults(anaObj->GetOutput(), saveToFile);

  printHists(anaObj->GetOutput());


  // plotting
  //anaObj->GetOutput()->SaveAs("test.root");
  // cleanup
  //delete anaObj;

  //delete phosRawStream;
  //delete rawChain;
  //delete chain;
}

void addFilesToChain(const TString rawFileList, TChain* chain)
{
  // Open rawFileList and add files contained to chain
  std::ifstream inList;
  inList.open( rawFileList.Data() );
  std::string line;
  while( std::getline(inList, line ) ) {
    printf("Add file %s to chain\n", line.c_str());
    chain->Add(line.c_str());
 }
  inList.close();
}

void printHists(AliTRUPedestalOutput* output){
  gROOT->SetStyle("Plain");
  gStyle->SetPadRightMargin(0.11);
  
  TCanvas* canvas =  new TCanvas;
  TH1F* pedhist = output->GetPedestals();
  pedhist->GetXaxis()->SetRangeUser(512-10, 512+10);
  canvas->SetLogy();
  pedhist->DrawCopy();

  canvas =  new TCanvas;
  TH1F* rmshist = output->GetPedestalRMS();
  rmshist->GetXaxis()->SetRangeUser(0, 5);
  canvas->SetLogy();
  rmshist->DrawCopy();

  canvas = new TCanvas;
  TH1F* pedId  = output->GetPedestalsId();
  pedId->DrawCopy();

  canvas = new TCanvas;
  TH1I* samples = output->GetPedestalSamples();
  samples->DrawCopy();
  
  canvas = new TCanvas;
  canvas->Divide(3, 2);
  canvas->cd(1);
  TH2F* ped2d_2 = output->GetPedestals2d(2);
  ped2d_2->GetZaxis()->SetRangeUser(512-10, 512+10);
  ped2d_2->DrawCopy("colz");

  canvas->cd(2);
  TH2F* ped2d_3 = output->GetPedestals2d(3);
  ped2d_3->GetZaxis()->SetRangeUser(512-10, 512+10);
  ped2d_3->DrawCopy("colz");

  canvas->cd(3);
  TH2F* ped2d_4 = output->GetPedestals2d(4);
  ped2d_4->GetZaxis()->SetRangeUser(512-10, 512+10);
  ped2d_4->DrawCopy("colz");

  canvas->cd(4);
  TH2F* rms2d_2 = output->GetPedestalRMS2d(2);
  //rms2d_2->GetZaxis()->SetRangeUser(0, 5);
  rms2d_2->DrawCopy("colz");

  canvas->cd(5);
  TH2F* rms2d_3 = output->GetPedestalRMS2d(3);
  //rms2d_3->GetZaxis()->SetRangeUser(0, 5);
  rms2d_3->DrawCopy("colz");

  canvas->cd(6);
  TH2F* rms2d_4 = output->GetPedestalRMS2d(4);
  //rms2d_4->GetZaxis()->SetRangeUser(0, 5);
  rms2d_4->DrawCopy("colz");
}

void saveResults(AliTRUPedestalOutput* output, TString saveToFile)
{
  if( ! output )
    gROOT->Error("truPedestalAnalysis.C::saveResults: no ouput", "no msgfmt");

  TTree* tree = new TTree("pedestalTree", "Pedestal Analysis Tree");
  TBranch* branch = tree->Branch("pedOutput", &output);

  Printf("Filling:");
  branch->Print();;
  tree->Fill();

  tree->SaveAs(saveToFile);
}
