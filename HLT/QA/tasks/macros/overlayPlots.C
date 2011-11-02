// $Id$
// 
//  Macro to overlay the histograms produced by  
//  HLT/QA/tasks/macros/drawTHnSparse.C
//  
//  It assumes a txt file where the input is specified in 
//  the following format:
//  
//   file1 legend1
//   file2 legend2
//  //file3 legend3
//  //file4 legend4
//  ...
//  So it is possible to "comment out" a file by a // in the beginning of the name. While reading the 
//  the names of the input files, the macro skips the ones that have the // in front of them.
//
//  Run it by:
//  aliroot overlayPlots.C++ 
// 
//  @ingroup alihlt_qa
//  @author Kalliopi.Kanaki@ift.uib.no 


#if !defined(__CINT__) || defined(__MAKECINT__)
#include "TSystem.h"
#include "TROOT.h"
#include "TFile.h"
#include "TString.h"
#include "TList.h"
#include "THnSparse.h"
#include "TCanvas.h"
#include "TText.h"
#include "TPaveText.h"
#include "TPaveStats.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TPad.h"
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <vector>
using std::endl;
#endif

void printStats(TH1D *h1, TH1D *h2);
void defineYaxisMax(TH1D *h1, TH1D *h2);

void overlayPlots(TString plottype="track" /*or event*/,const char* option="HLT"/* or "OFF" */,  Bool_t bAddRunName=kTRUE, Bool_t bDrawNormalized=kTRUE, string fi="files.txt"){


  printf("test\n"); 
  gROOT->SetStyle("Plain");
  gStyle->SetPalette(1);
  gStyle->SetOptStat("emr");
  gStyle->SetTitleX(gStyle->GetPadLeftMargin());
 
  char filenames[100];
  sprintf(filenames,"%s",fi.c_str());
  ifstream infile;
  infile.open(filenames);
  if(!infile){
    printf("File %s does not exist", fi.data());
    return;
  }
  string c;
  string f;
  int nr_textfile = 0;

  vector<string> file;
  vector<string> cutnames;
  string plotname="";
  string t("//");

  while(1){
    infile >> f >> c;
    if(!infile.good()) break;
    if (f.compare(0, t.length(), t) == 0) continue;
    file.push_back(f);
    cutnames.push_back(c);
    nr_textfile++;
  }
  infile.close();
  cout << cutnames.size()<< endl;

  printf("Number of files: %d\n", nr_textfile);
  
  TCanvas *ca;
  TFile   *ff; 
  TPad    *pad; 
  TH1D    *g[nr_textfile];
  
  TCanvas *d = new TCanvas("d",Form("Compare %s distributions",option),1600,1000);
  d->Divide(4,2);
  //d->Divide(3,2);
  TLegend *l;
  if(plottype.Contains("event"))
    l = new TLegend(0,0,0.95,0.95);
  else
    l = new TLegend(0.35,0.3,0.89,0.8);
  l->SetFillColor(10);
  l->SetLineColor(10);
  
  char cut[100];  

  int nr_pads=0;
  if(plottype.Contains("event"))nr_pads=7;
  else nr_pads=9;

  for(int j=1; j<nr_pads; j++){ 
    for(int i=0; i<nr_textfile; i++){ 
       
      ff = TFile::Open(file[i].data());   
      if(!ff || ff->IsZombie()){
	printf("Non-existent, corrupted or zombie file %s\n", file[i].data());
	return;
      } 
      if(plottype.Contains("event"))
	ca  = (TCanvas*)ff->GetObjectUnchecked("can1");
      else		    
	ca  = (TCanvas*)ff->GetObjectUnchecked("can3");		    
      if(!ca){
	printf("Empty canvas in file %s.\n", file[i].data());
	continue;
      }
      if(plottype.Contains("event"))
	pad = (TPad*)ca->GetListOfPrimitives()->FindObject(Form("can1_%d",j));         	
      else
	pad = (TPad*)ca->GetListOfPrimitives()->FindObject(Form("can3_%d",j));         	
      if(!pad){
	printf("Empty pad in canvas %s.\n", ca->GetName());
	continue;	     
      }
      if(plottype.Contains("event"))
	g[i] =(TH1D*)pad->FindObject(Form("fEvent%s_proj_%d",option,j-1));
      else
	g[i] =(TH1D*)pad->FindObject(Form("fTrack%s_proj_%d",option,j-1));
      if(!g[i]){
	printf("Empty histogram for i=%d, file %s.\n", i, file[i].data());
	continue;
      }
        
      d->cd(j);
		
      if(i==0){
	g[i]->SetLineColor(kBlack);
	//	  defineYaxisMax(hlt, off);
	if(bDrawNormalized) 
	  g[i]->DrawNormalized();
	else
	  g[i]->Draw();
	if(option=="OFF"){
	  TPaveStats *st = (TPaveStats*)g[i]->FindObject("stats");
	  st->SetTextColor(kBlack);
	  d->Update();
	}
      }
      else { 
	if(i<4)
	  g[i]->SetLineColor(i+1);
	else
	  g[i]->SetLineColor(i+2);
	defineYaxisMax(g[0], g[i]);
	if(bDrawNormalized)  g[i]->DrawNormalized("sames");
	else g[i]->Draw("sames");

      }

      if(!bDrawNormalized)					 
	if(i>0) printStats(g[i-1], g[i]);
	
      ff->Close();
      sprintf( cut,"%s",cutnames[i].c_str() );
      if((j==2&&plottype.Contains("track")) || ((j==6&&plottype.Contains("event")))) {

	l->AddEntry(g[i],cut,"l");	    	
	//	  cutnames[i].resize(6);
	plotname+="_"+cutnames[i];
	//cout << "Adding Run: " << plotname <<endl;
	cout << "Adding Run: " << cutnames[i] <<endl;
      }
      else continue;
    }

    if(j==2 && plottype.Contains("track") ) l->Draw("same");
    if(plottype.Contains("event") && j==6) {
      d->cd(7);
      l->Draw();

    }
    d->Update();
  }

  d->Update(); 

  sprintf(filenames,"%s",plotname.c_str());
  cout << filenames << endl;
  if( bAddRunName){
    d->SaveAs(Form("overlay_%s_for%s.root",option,filenames));
    d->Print(Form("overlay_%s_for%s.png",option,filenames));
  }
  else{
    d->SaveAs(Form("overlay_%s.root",option));
    d->Print(Form("overlay_%s.png",option));
  }

  return;
}

void printStats(TH1D *h1, TH1D *h2){  
  gPad->Update();
  TPaveStats *st1 = (TPaveStats*)h1->FindObject("stats");
  st1->SetLineColor(0);

  gPad->Update();
  TPaveStats *st2 = (TPaveStats*)h2->FindObject("stats");
  st2->SetY2NDC(st1->GetY1NDC()-0.05);
  st2->SetY1NDC(st2->GetY2NDC()-TMath::Abs(st1->GetY1NDC()-st1->GetY2NDC()));
  st2->SetLineColor(0);
  st2->SetTextColor(h2->GetLineColor());
  st2->SetFillStyle(0);
  st2->Draw();
  return;
}

void defineYaxisMax(TH1D *h1, TH1D *h2){   
  //Y axis
  if(h1->GetMaximum() > h2->GetMaximum()) h2->SetMaximum(1.1*h1->GetMaximum());
  else h1->SetMaximum(1.1*h2->GetMaximum());
  
  h1->SetMinimum(0);
  h2->SetMinimum(0);
 
  // X axis  
  double xmin, xmax;  
  if(h1->GetBinLowEdge(1) > h2->GetBinLowEdge(1)) xmin = h1->GetBinLowEdge(1);
  else xmin = h2->GetBinLowEdge(1);
  if(h1->GetBinLowEdge(h1->GetNbinsX()+1) > h2->GetBinLowEdge(h1->GetNbinsX()+1)) xmax = h1->GetBinLowEdge(h1->GetNbinsX()+1);
  else xmax = h2->GetBinLowEdge(h2->GetNbinsX()+1);
  
  h2->SetAxisRange(xmin, xmax, "X");
  return;
}
