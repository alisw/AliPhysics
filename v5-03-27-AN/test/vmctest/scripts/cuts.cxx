/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

// $Id$

// Program to generate the plots with energy cuts and corresponding
// cuts in range for all materials per detector.
// The program uses a text file output from geant4_vmc
// (in development version only).
//
// By I.Hrivnacova, IPN Orsay

#include <Rtypes.h>
#include <TPaveText.h>
#include <TString.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TAxis.h>
#include <TLegend.h>

#include <map>
#include <vector>
#include <string>
#include <fstream>
#include <iostream>

using namespace std;

struct MaterialData {
  MaterialData(string name, 
               Double_t x1, Double_t x2, Double_t x3, 
               Double_t x4, Double_t x5, Double_t x6)
    : matName(name), 
      rangeGam(x1), 
      rangeEle(x2), 
      cutGam(x3), 
      cutEle(x4),
      cutGamLimit(x5), 
      cutEleLimit(x6) {}
  string matName;
  Double_t rangeGam;
  Double_t rangeEle;
  Double_t cutGam;
  Double_t cutEle;
  Double_t cutGamLimit;
  Double_t cutEleLimit;
};  
  
map<string, vector<MaterialData> > materialDataMap;  

void readCuts() 
{
  std::ifstream input("cuts.txt");
  while ( ! input.eof() ) {
    string matName;
    Int_t number;
    Double_t x1, x2, x3, x4, x5, x6;
    input >> number;
    if ( input.eof() ) break;
    input  >> matName >> x1 >> x2 >> x3 >> x4 >> x5 >> x6;
    // cout << ".. reading " << number << " " << matName << endl;
    
    if ( x1 > 1e10 ) x1 = 1e04;
    if ( x2 > 1e10 ) x2 = 1e04;
    if ( x3 > 1e10 ) x3 = 1e04;
    if ( x4 > 1e10 ) x4 = 1e04;
    if ( x5 > 1e10 ) x5 = 1e04;
    if ( x6 > 1e10 ) x6 = 1e04;
    
    string detName = matName;
    detName.erase(detName.find('_'), detName.size()-detName.find('_'));

    std::map<string, std::vector<MaterialData> >::iterator it
      = materialDataMap.find(detName);
    if ( it == materialDataMap.end() ) {
      materialDataMap.insert(
        std::pair<string, vector<MaterialData> >(detName, vector<MaterialData>()));
    }
    materialDataMap[detName].
      push_back(MaterialData(matName, x1, x2, x3, x4, x5, x6));
  }
}        
  
void plotCuts(const string& detName) 
{
  cout << "Processing " << detName << " ..." << endl;
  
  //if ( detName == "DIPO" ) return;
  //if ( detName == "HALL" ) return;
  //if ( detName == "SHIL" ) return;
  //if ( detName == "ZDC" ) return;

  TPaveText* paveText = new TPaveText(0.1, 0.1, 0.98, 0.98);
  paveText->SetTextAlign(11);

  std::vector<MaterialData> materialData 
    = materialDataMap[detName];

  Double_t matNumber[300];
  Double_t rangeGam[300];
  Double_t rangeEle[300];
  Double_t cutGam[300];
  Double_t cutEle[300];
  Double_t cutGamLimit[300];
  Double_t cutEleLimit[300];
  std::vector<string> matNames;
 
  for (UInt_t i=0; i<materialData.size(); i++ ) {
  
    matNumber[i] = i;
    rangeGam[i] = materialData[i].rangeGam;  
    rangeEle[i] = materialData[i].rangeEle;
    cutGam[i] = materialData[i].cutGam;
    cutEle[i] = materialData[i].cutEle;
    cutGamLimit[i] = materialData[i].cutGamLimit;
    cutEleLimit[i] = materialData[i].cutEleLimit;
   
    TString legend;
    legend += i;
    legend += "  ";
    legend += materialData[i].matName.data();    
    paveText->AddText(legend.Data());
  }
  
  TGraph* gr1 = new TGraph(materialData.size(), matNumber, rangeGam); 
  TGraph* gr2 = new TGraph(materialData.size(), matNumber, cutGam); 
  TGraph* gr3 = new TGraph(materialData.size(), matNumber, cutGamLimit); 
  TGraph* gr4 = new TGraph(materialData.size(), matNumber, rangeEle); 
  TGraph* gr5 = new TGraph(materialData.size(), matNumber, cutEle); 
  TGraph* gr6 = new TGraph(materialData.size(), matNumber, cutEleLimit); 

  // gamma range cut
  gr1->SetMarkerColor(kBlue);
  gr1->SetMarkerStyle(22);
  gr1->SetTitle("Range cut for gamma");
  gr1->GetXaxis()->SetTitle("Material number");
  gr1->GetYaxis()->SetTitle("Range [mm]");
  gr1->SetMinimum(1e-03);
  gr1->SetMaximum(1e+04);

  // gamma energy threshold
  gr2->SetMarkerColor(kBlue);
  gr2->SetMarkerStyle(22);
  gr2->SetTitle("Energy threshold for gamma");
  gr2->GetXaxis()->SetTitle("Material number");
  gr2->GetYaxis()->SetTitle("Energy [MeV]");
  gr2->SetMinimum(1e-04);
  gr2->SetMaximum(1e+04);

  // gamma user limit
  gr3->SetMarkerColor(kBlack);
  gr3->SetMarkerStyle(23);
  gr3->SetTitle("User limit for gamma (GUTGAM)");
  gr3->GetXaxis()->SetTitle("Material number");
  gr3->GetYaxis()->SetTitle("Energy [MeV]");
  gr3->SetMinimum(1e-04);
  gr3->SetMaximum(1e+04);

  // e- range cut
  gr4->SetMarkerColor(kRed);
  gr4->SetMarkerStyle(22);
  gr4->SetTitle("Range cut for e-");
  gr4->GetXaxis()->SetTitle("Material number");
  gr4->GetYaxis()->SetTitle("Range [mm]");
  gr4->GetYaxis()->SetRange(1e-03, 1e+04);
  gr4->SetMinimum(1e-03);
  gr4->SetMaximum(1e+04);

  // e- energy threshold
  gr5->SetMarkerColor(kRed);
  gr5->SetMarkerStyle(22);
  gr5->SetTitle("Energy threshold for e-");
  gr5->GetXaxis()->SetTitle("Material number");
  gr5->GetYaxis()->SetTitle("Energy [MeV]");
  gr5->SetMinimum(1e-04);
  gr5->SetMaximum(1e+04);

  // e- user limit
  gr6->SetMarkerColor(kBlack);
  gr6->SetMarkerStyle(23);
  gr6->SetTitle("User limit for e- (CUTELE)");
  gr6->GetXaxis()->SetTitle("Material number");
  gr6->GetYaxis()->SetTitle("Energy [MeV]");
  gr6->SetMinimum(1e-04);
  gr6->SetMaximum(1e+04);

  //TLegend* leg = new TLegend(0.1,0.7,0.48,0.9);
  TLegend* leg = new TLegend(0.1, 0.1, 0.98, 0.98);
  leg->SetHeader("Energy threshold/User cuts");
  leg->AddEntry(gr2,"Energy threshold for gamma","P");
  leg->AddEntry(gr3,"User cut for gamma (CUGAM)","P");
  leg->AddEntry(gr5,"Energy threshold for e-","P");
  leg->AddEntry(gr6,"User cut for e- (CUTELE)","P");

  TCanvas* canvas = new TCanvas("Range cuts", "Range cuts", 1200, 800);
  canvas->Divide(3,2);
  
  // Draw graph using the function Draw()
  canvas->cd(1);
  canvas->cd(1)->SetLogy();
  gr1->Draw("AP");
 
  canvas->cd(2);
  canvas->cd(2)->SetLogy();
  gr2->Draw("AP");

  //canvas->cd(3);
  //canvas->cd(3)->SetLogy();
  gr3->Draw("P");
 
  canvas->cd(3);
  leg->Draw();

  canvas->cd(4);
  canvas->cd(4)->SetLogy();
  gr4->Draw("AP");

  canvas->cd(5);
  canvas->cd(5)->SetLogy();
  gr5->Draw("AP");
 
  //canvas->cd(6);
  //canvas->cd(6)->SetLogy();
  gr6->Draw("P");

  canvas->cd(6);

  //TCanvas* canvas2 = new TCanvas("Materials", "Materials", 200, 800);
  //canvas2->cd();
  paveText->Draw();
  
  string outName = detName;
  outName += "_cuts.gif";
  canvas->SaveAs(outName.data());
  
  //string outName2 = detName;
  //outName2 += "_cuts_legend.gif";
  //canvas2->SaveAs(outName2.data());
}

int main() {
  
  readCuts();
  
  std::map<string, std::vector<MaterialData> >::iterator it;
  for ( it= materialDataMap.begin(); it != materialDataMap.end(); it++ )
    plotCuts(it->first);
    
  return 0;  
}  



