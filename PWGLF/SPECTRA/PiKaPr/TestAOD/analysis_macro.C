#include "AliAODTrack.h"
#include <iostream>
#include "TFile.h"
#include "TH1I.h"
#include "TH2.h"
#include "AliSpectraAODHistoManager.h"
#include "AliSpectraAODEventCuts.h"
#include "AliSpectraAODTrackCuts.h"
#include "TCanvas.h"
#include "TLegend.h"

using namespace std;

void analysis_macro(const char * dataFile = "Pt.AOD.1._data_ptcut.root", const char * mcFile ="Pt.AOD.1._MC.root") {
   
  // gSystem->Load("libTree.so");
  // gSystem->Load("libGeom.so");
  // gSystem->Load("libVMC.so");
  // gSystem->Load("libPhysics.so");
  // gSystem->Load("libSTEERBase.so");
  // gSystem->Load("libESD.so");
  // gSystem->Load("libAOD.so");
  // gSystem->Load("libANALYSIS.so");
  // gSystem->Load("libANALYSISalice.so");
  // gSystem->Load("libANALYSIS");
  // gSystem->Load("libANALYSISalice");
  // gROOT->LoadMacro("$ALICE_ROOT/PWG0/CreateESDChain.C");
  // gSystem->AddIncludePath("-I$ALICE_ROOT/include");
  // gStyle->SetPalette(1);
  // gStyle->SetFillColor(kWhite);

  // gROOT->LoadMacro("AliSpectraAODTrackCuts.cxx+g");
  // gROOT->LoadMacro("AliSpectraAODEventCuts.cxx+g");
  // gROOT->LoadMacro("AliSpectraAODHistoManager.cxx+g");
  // gROOT->LoadMacro("AliAnalysisTaskSpectraAOD.cxx+g");
 
  const char * histoname[] =
    {
      // names of histos we want to draw. please observe order of histos, elsewise name will not make sense
      "#frac{MC_{#sigma, rec}}{Data}, P^{+}",
      "#frac{MC_{#sigma, rec}}{Data}, K^{+}",
      "#frac{MC_{#sigma, rec}}{Data}, #pi^{+}",
      "#frac{MC_{#sigma, rec}}{Data}, P^{-}",
      "#frac{MC_{#sigma, rec}}{Data}, K^{-}",
      "#frac{MC_{#sigma, rec}}{Data}, #pi^{-}",      
          
      "#frac{MC_{#sigma, rec, prim}}{Data}, P^{+}",
      "#frac{MC_{#sigma, rec, prim}}{Data}, K^{+}",
      "#frac{MC_{#sigma, rec, prim}}{Data}, #pi^{+}",
      "#frac{MC_{#sigma, rec, prim}}{Data}, P^{-}",
      "#frac{MC_{#sigma, rec, prim}}{Data}, K^{-}",
      "#frac{MC_{#sigma, rec, prim}}{Data}, #pi^{-}",     
         
      "#frac{MC_{#sigma, rec, sec}}{Data}, P^{+}",
      "#frac{MC_{#sigma, rec, sec}}{Data}, K^{+}",
      "#frac{MC_{#sigma, rec, sec}}{Data}, #pi^{+}",
      "#frac{MC_{#sigma, rec, sec}}{Data}, P^{-}",
      "#frac{MC_{#sigma, rec, sec}}{Data}, K^{-}",
      "#frac{MC_{#sigma, rec, sec}}{Data}, #pi^{-}",     
         
      "#frac{MC_{true, gen}}{MC_{#sigma, rec}}, P^{+}",
      "#frac{MC_{true, gen}}{MC_{#sigma, rec}}, K^{+}",
      "#frac{MC_{true, gen}}{MC_{#sigma, rec}}, #pi^{+}",
      "#frac{MC_{true, gen}}{MC_{#sigma, rec}}, P^{-}",
      "#frac{MC_{true, gen}}{MC_{#sigma, rec}}, K^{-}",
      "#frac{MC_{true, gen}}{MC_{#sigma, rec}}, #pi^{-}",     
              
      "#frac{MC_{true, gen}}{MC_{#sigma, rec}} x Data, P^{+}",
      "#frac{MC_{true, gen}}{MC_{#sigma, rec}} x Data, K^{+}",
      "#frac{MC_{true, gen}}{MC_{#sigma, rec}} x Data, #pi^{+}",
      "#frac{MC_{true, gen}}{MC_{#sigma, rec}} x Data, P^{-}",
      "#frac{MC_{true, gen}}{MC_{#sigma, rec}} x Data, K^{-}",
      "#frac{MC_{true, gen}}{MC_{#sigma, rec}} x Data, #pi^{-}",     

      "Data, P^{+}",
      "Data, K^{+}",
      "Data, #pi^{+}",
      "Data, P^{-}",
      "Data, K^{-}",
      "Data, #pi^{-}",      

      "#frac{MC_{rec #sigma}}{MC_{true, rec}}, P^{+}",
      "#frac{MC_{rec #sigma}}{MC_{true, rec}}, K^{+}",
      "#frac{MC_{rec #sigma}}{MC_{true, rec}}, #pi^{+}",
      "#frac{MC_{rec #sigma}}{MC_{true, rec}}, P^{-}",
      "#frac{MC_{rec #sigma}}{MC_{true, rec}}, K^{-}",
      "#frac{MC_{rec #sigma}}{MC_{true, rec}}, #pi^{-}",     

      "#frac{MC_{rec #sigma sec}}{MC_{#sigma, rec}}, P^{+}",
      "#frac{MC_{rec #sigma sec}}{MC_{#sigma, rec}}, K^{+}",
      "#frac{MC_{rec #sigma sec}}{MC_{#sigma, rec}}, #pi^{+}",
      "#frac{MC_{rec #sigma sec}}{MC_{#sigma, rec}}, P^{-}",
      "#frac{MC_{rec #sigma sec}}{MC_{#sigma, rec}}, K^{-}",
      "#frac{MC_{rec #sigma sec}}{MC_{#sigma, rec}}, #pi^{-}",     
    };

  const char * histotitle[] =
    {
      // titles of histos we want to draw. please observe order of histos, elsewise name will not make sense
      "MC #sigma, rec / Data, P^{+}",
      "MC #sigma, rec / Data, K^{+}",
      "MC #sigma, rec / Data, #pi^{+}",
      "MC #sigma, rec / Data, P^{-}",
      "MC #sigma, rec / Data, K^{-}",
      "MC #sigma, rec / Data, #pi^{-}",     
   
      "MC #sigma, rec, prim / Data, P^{+}",
      "MC #sigma, rec, prim  / Data, K^{+}",
      "MC #sigma, rec, prim  / Data, #pi^{+}",
      "MC #sigma, rec, prim  / Data, P^{-}",
      "MC #sigma, rec, prim  / Data, K^{-}",
      "MC #sigma, rec, prim  / Data, #pi^{-}",      
    
      "MC #sigma, rec, sec / Data, P^{+}",
      "MC #sigma, rec, sec / Data, K^{+}",
      "MC #sigma, rec, sec / Data, #pi^{+}",
      "MC #sigma, rec, sec / Data, P^{-}",
      "MC #sigma, rec, sec / Data, K^{-}",
      "MC #sigma, rec, sec / Data, #pi^{-}",      
          
      "Correction factor, P^{+}",
      "Correction factor, K^{+}",
      "Correction factor, #pi^{+}",
      "Correction factor, P^{-}",
      "Correction factor, K^{-}",
      "Correction factor, #pi^{-}",      

      "True Data, P^{+}",
      "True Data, K^{+}",
      "True Data, #pi^{+}",
      "True Data, P^{-}",
      "True Data, K^{-}",
      "True Data, #pi^{-}",      

      "Raw Data, P^{+}",
      "Raw Data, K^{+}",
      "Raw Data, #pi^{+}",
      "Raw Data, P^{-}",
      "Raw Data, K^{-}",
      "Raw Data, #pi^{-}",      

      "Identificaiton quality, P^{+}",
      "Identificaiton quality, K^{+}",
      "Identificaiton quality, #pi^{+}",
      "Identificaiton quality, P^{-}",
      "Identificaiton quality, K^{-}",
      "Identificaiton quality, #pi^{-}",     

      "Contribution of secondaries, P^{+}",
      "Contribution of secondaries, K^{+}",
      "Contribution of secondaries, #pi^{+}",
      "Contribution of secondaries, P^{-}",
      "Contribution of secondaries, K^{-}",
      "Contribution of secondaries, #pi^{-}",     
    };

  // Open root MC file and get classes
  cout << "Analysis Macro" << endl;
  cout << "  > Reading MC data" << endl;
  TFile *_mc = TFile::Open(mcFile);
  AliSpectraAODHistoManager* hman_mc = (AliSpectraAODHistoManager*) _mc->Get("SpectraHistos");
  AliSpectraAODEventCuts* ecuts_mc = (AliSpectraAODEventCuts*) _mc->Get("Event Cuts");
  AliSpectraAODTrackCuts* tcuts_mc = (AliSpectraAODTrackCuts*) _mc->Get("Track Cuts");
  // print info about mc track and Event cuts
  cout << " -- Info about MC -- "<< endl;
  ecuts_mc->PrintCuts();
  tcuts_mc->PrintCuts();
  // get the mc histo's necessary to write the contamination histos from the rootfile
  TH1F* spectrahistos_mc[30];
  cout << " -- Reading and normalizing histograms -- " << endl;
  for( Int_t i = 0 ; i <= AliSpectraNameSpace::kNPtSpecies ; i++ ) { 
    spectrahistos_mc[i] = dynamic_cast<TH1F*>(hman_mc->GetHistogram((AliSpectraNameSpace::AODPtHist_t)i)); 
    Double_t events_mc =  1. / ecuts_mc->NumberOfEvents();
    if(!spectrahistos_mc[i]) {
      continue;
    }    
    spectrahistos_mc[i]->Scale(events_mc, "width");
  }
  cout << "  > Reading real data " << endl;
  // proceed likewise for data
  TFile *_data = TFile::Open(dataFile);
  AliSpectraAODHistoManager* hman_data = (AliSpectraAODHistoManager*) _data->Get("SpectraHistos");
  AliSpectraAODEventCuts* ecuts_data = (AliSpectraAODEventCuts*) _data->Get("Event Cuts");
  AliSpectraAODTrackCuts* tcuts_data = (AliSpectraAODTrackCuts*) _data->Get("Track Cuts");
  // print info about track and Event cuts
  cout << " -- Info about data -- " << endl;
  ecuts_data->PrintCuts();
  tcuts_data->PrintCuts();
  // get the histo's necessary to write the contamination histos from the rootfile
  TH1F* spectrahistos_data[30];
  cout << " -- Reading and normalizing histograms -- " << endl;
  for( Int_t i = 0 ; i <= AliSpectraNameSpace::kNPtSpecies ; i++ ) { 
    spectrahistos_data[i] = dynamic_cast<TH1F*>(hman_data->GetHistogram((AliSpectraNameSpace::AODPtHist_t)i));
    if(!spectrahistos_data[i]) {
      continue;
    }    
    Double_t events_data =  1. / ecuts_data->NumberOfEvents();
    //normalize the histos    
    spectrahistos_data[i]->Scale(events_data, "width");
  }
  //create output file and select it as the active file
  TFile * output = new TFile("analysis_output.root", "RECREATE");
  output->cd();
  // Write the scaled data and MC histos in the output file for convenience
  output->mkdir("Data");
  output->cd("Data");
  for(Int_t ihist = 0; ihist <= AliSpectraNameSpace::kNPtSpecies; ihist++){
      spectrahistos_data[ihist]->Write();
  }
  output->mkdir("MC");
  output->cd("MC");
  for(Int_t ihist = 0; ihist <= AliSpectraNameSpace::kNPtSpecies; ihist++){
      spectrahistos_mc[ihist]->Write();
  }
  output->mkdir("ANALYSIS");
  output->cd("ANALYSIS");
  // perform some analysis!
  //  next loop: MC_Sigma_Rec over Rec_Data
  TCanvas* c1 = new TCanvas("MC reconstructed over Data","MC reconstructed over Data");
  c1->Divide(3,2);
  TH1F* kRatioPtMCSigmaRecSigmaRec[6];
  for (Int_t i = 0 ; i < 6 ; i ++ )
    {
      kRatioPtMCSigmaRecSigmaRec[i] = (TH1F*)(spectrahistos_mc[i])->Clone(histotitle[i]);
      kRatioPtMCSigmaRecSigmaRec[i]->Divide(spectrahistos_data[i]);
      c1->cd(1+i);
      kRatioPtMCSigmaRecSigmaRec[i]->SetTitle(histoname[i]);
      kRatioPtMCSigmaRecSigmaRec[i]->GetYaxis()->SetTitle("#frac{d^{2} N}{dy dp_{T}} (c / GeV)");
      kRatioPtMCSigmaRecSigmaRec[i]->Draw();
      kRatioPtMCSigmaRecSigmaRec[i]->Write();
    }
  // next loop: MC_sigma_primaries over Rec_data
  TCanvas* c2 = new TCanvas("MC reconstructed primaries over Data", "MC reconstructed primaries over Data");
  c2->Divide(3,2);
  TH1F* kRatioPtMCSigmaRecPrimarySigmaRec[6];
  for (Int_t i = 0 ; i < 6 ; i ++ )
    {
      kRatioPtMCSigmaRecPrimarySigmaRec[i] = (TH1F*)(spectrahistos_mc[i+AliSpectraNameSpace::kHistPtRecSigmaPrimaryProtonPlus])->Clone(histotitle[i+6]);
      kRatioPtMCSigmaRecPrimarySigmaRec[i]->Divide(spectrahistos_data[i]);
      c2->cd(1+i);
      kRatioPtMCSigmaRecPrimarySigmaRec[i]->SetTitle(histoname[i+6]);
      kRatioPtMCSigmaRecPrimarySigmaRec[i]->GetYaxis()->SetTitle("#frac{d^{2} N}{dy dp_{T}} (c / GeV)");
      kRatioPtMCSigmaRecPrimarySigmaRec[i]->Draw();
      kRatioPtMCSigmaRecPrimarySigmaRec[i]->Write();
    }
  // next loop: MC_Sigma_secondaries over Rec_data
  TCanvas* c3 = new TCanvas("MC reconstructed secondaries over Data","MC reconstructed secondaries over Data");
  c3->Divide(3,2);
  TH1F* kRatioPtMCSigmaRecSigmaRecSecondary[6];
  for (Int_t i = 0 ; i < 6 ; i ++ )
    {
      kRatioPtMCSigmaRecSigmaRecSecondary[i] = (TH1F*)(spectrahistos_mc[i+AliSpectraNameSpace::kHistPtRecSigmaSecondaryProtonPlus])->Clone(histotitle[i+12]);
      kRatioPtMCSigmaRecSigmaRecSecondary[i]->Divide(spectrahistos_data[i]);
      c3->cd(1+i);
      kRatioPtMCSigmaRecSigmaRecSecondary[i]->SetTitle(histoname[12+i]);
      kRatioPtMCSigmaRecSigmaRecSecondary[i]->GetYaxis()->SetTitle("#frac{d^{2} N}{dy dp_{T}} (c / GeV)");
      kRatioPtMCSigmaRecSigmaRecSecondary[i]->Draw();
      kRatioPtMCSigmaRecSigmaRecSecondary[i]->Write();
    }
  //
  TCanvas* c4 = new TCanvas("MC true generated over MC true reconstructed","MC true generated over MC true reconstructed");
  c4->Divide(3,2);
  TH1F* kRatioPtMCTrueMCRec[6];
  for (Int_t i = 0 ; i < 6 ; i ++ )
    {
      kRatioPtMCTrueMCRec[i] = (TH1F*)(spectrahistos_mc[i+AliSpectraNameSpace::kHistPtGenTruePrimaryProtonPlus])->Clone(histotitle[i+18]);
      kRatioPtMCTrueMCRec[i]->Divide(spectrahistos_mc[i+AliSpectraNameSpace::kHistPtRecSigmaPrimaryProtonPlus]);
      c4->cd(1+i);
      kRatioPtMCTrueMCRec[i]->SetTitle(histoname[18+i]);
      kRatioPtMCTrueMCRec[i]->GetYaxis()->SetTitle("#frac{d^{2} N}{dy dp_{T}} (c / GeV)");
      kRatioPtMCTrueMCRec[i]->Draw();
      kRatioPtMCTrueMCRec[i]->Write();
    }
  
  TCanvas* c5 = new TCanvas("True Data","True Data");
  c5->Divide(3,2);
  TH1F* kRatioPtTrueData[6];
  for (Int_t i = 0 ; i < 6 ; i ++ )
    {
      kRatioPtTrueData[i] = (TH1F*)(spectrahistos_mc[i+AliSpectraNameSpace::kHistPtGenTruePrimaryProtonPlus])->Clone(histotitle[i+24]);
      kRatioPtTrueData[i]->Divide(spectrahistos_mc[i]); // Michele: Since we are correcting all contributions at once, this should be exaclty the same as in the data!
      kRatioPtTrueData[i]->Multiply(spectrahistos_data[i]);
      c5->cd(1+i);
      kRatioPtTrueData[i]->SetTitle(histoname[24+i]);
      kRatioPtTrueData[i]->GetYaxis()->SetTitle("#frac{d^{2} N}{dy dp_{T}} (c / GeV)");
      kRatioPtTrueData[i]->Draw();
      kRatioPtTrueData[i]->Write();
    }
 
  TCanvas* c6 = new TCanvas("Raw data", "Raw data");
  c6->Divide(3,2);
  TH1F* kRawData[6];
  for (Int_t i = 0 ; i < 6 ; i ++ )
    {
      kRawData[i] = (TH1F*)(spectrahistos_data[i+AliSpectraNameSpace::kHistPtRecSigmaProtonPlus])->Clone(histotitle[i+30]);
      c6->cd(1+i);
      kRawData[i]->SetTitle(histoname[30+i]);
      kRawData[i]->GetYaxis()->SetTitle("#frac{d^{2} N}{dy dp_{T}} (c / GeV)");
      kRawData[i]->Draw();
      kRawData[i]->Write();
    }
  
  TCanvas* c7 = new TCanvas("contamination", "contamination");
  c7->Divide(3,2);
  TH1F* kContamination[6];
  for (Int_t i = 0 ; i < 6 ; i ++ )
    {
      kContamination[i] = (TH1F*)(spectrahistos_mc[i+AliSpectraNameSpace::kHistPtRecSigmaProtonPlus])->Clone(histotitle[i+36]);
      kContamination[i]->Divide(spectrahistos_mc[i+AliSpectraNameSpace::kHistPtRecTrueProtonPlus]);
      c7->cd(1+i);
      kContamination[i]->SetTitle(histoname[36+i]);
      kContamination[i]->GetYaxis()->SetTitle("#frac{d^{2} N}{dy dp_{T}} (c / GeV)");
      kContamination[i]->Draw();
      kContamination[i]->Write();
    }
  
  TCanvas* c8 = new TCanvas("contamination secondary", "contamination secondary");
  c8->Divide(3,2);
  TH1F* kSecContamination[6];
  for (Int_t i = 0 ; i < 6 ; i ++ )
    {
      kSecContamination[i] = (TH1F*)(spectrahistos_mc[i+AliSpectraNameSpace::kHistPtRecSigmaSecondaryProtonPlus])->Clone(histotitle[i+42]);
      kSecContamination[i]->Divide(spectrahistos_mc[i+AliSpectraNameSpace::kHistPtRecSigmaProtonPlus]);
      c8->cd(1+i);
      kSecContamination[i]->SetTitle(histoname[42+i]);
      kSecContamination[i]->GetYaxis()->SetTitle("#frac{d^{2} N}{dy dp_{T}} (c / GeV)");
      kSecContamination[i]->Draw();
      kSecContamination[i]->Write();
    }
   
  TCanvas* c9 = new TCanvas("true data and generated", "true data and generated");
  c9->Divide(3,2);
  TH1F* ktrueandgenerated[6];
  for (Int_t i = 0 ; i < 6 ; i ++ )
    {
      // kSecContamination[i] = (TH1F*)(spectrahistos_mc[i+AliSpectraNameSpace::kHistPtRecSigmaSecondaryProtonPlus])->Clone(histotitle[i+42]);
      // kSecContamination[i]->Divide(spectrahistos_mc[i+AliSpectraNameSpace::kHistPtRecSigmaProtonPlus]);
      c9->cd(1+i);
      // kSecContamination[i]->SetTitle(histoname[42+i]);
      // kSecContamination[i]->GetYaxis()->SetTitle("#frac{d^{2} N}{dy dp_{T}} (c / GeV)");
      kRatioPtTrueData[i]->Draw();
      spectrahistos_mc[i+AliSpectraNameSpace::kHistPtGenTruePrimaryProtonPlus]->SetMarkerColor(kRed);
      spectrahistos_mc[i+AliSpectraNameSpace::kHistPtGenTruePrimaryProtonPlus]->SetMarkerStyle(kOpenCircle);
      spectrahistos_mc[i+AliSpectraNameSpace::kHistPtGenTruePrimaryProtonPlus]->SetLineColor(kRed);
      spectrahistos_mc[i+AliSpectraNameSpace::kHistPtGenTruePrimaryProtonPlus]->Draw("same");

      // legend
      TLegend * leg = new TLegend(0.6,0.64,0.89,0.79);  //coordinates are fractions
      leg->AddEntry(spectrahistos_mc[i+AliSpectraNameSpace::kHistPtGenTruePrimaryProtonPlus],"Generated","p"); 
      leg->AddEntry(kRatioPtTrueData[i],"Corrected","p"); 
      leg->Draw();
    }
   
  cout << "Analysis complete" << endl;
  cout << " - Please exit root to close and save output file " << endl;
}
