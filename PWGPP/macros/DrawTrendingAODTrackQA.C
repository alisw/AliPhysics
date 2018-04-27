#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TTree.h>
#include <TH1.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#endif

// Macro to draw the trending information from the output of PlotAODtrackQA.C

TH1F* CreateHisto(TString nam, Int_t tote);

Bool_t DrawTrendingAODTrackQA(TString mergedTrendFile = "trending.root", 
			    TString treename="trendingTrack"){

  TFile *fin = TFile::Open(mergedTrendFile.Data());
  if(!fin){
    printf("Cannot open file with AOD track QA trending: %s\n",mergedTrendFile.Data());
    return kFALSE;
  }
  TTree * ttree = (TTree*) fin->Get(treename.Data());
  if (!ttree){
    printf("Trending tree not found in file %s\n",mergedTrendFile.Data());
    return kFALSE;
  }
  Int_t nrun;
  ttree->SetBranchAddress("nrun",&nrun);
  TObjArray* lb=(TObjArray*)ttree->GetListOfBranches();
  Int_t nVars=lb->GetEntries();
  Int_t totEnt=ttree->GetEntries();
  Float_t* vectVars=new Float_t[nVars-1];
  Float_t* vectErrs=new Float_t[nVars-1];
  for(Int_t j=0; j<nVars-1; j++){
    vectVars[j]=-999.;
    vectErrs[j]=0.00001;
  }
  TH1F** htr=new TH1F*[nVars-1];
  Int_t kv=0;
  for(Int_t j=0; j<nVars; j++){
    TBranch* br=(TBranch*)lb->At(j);
    printf("Branch %d  %s\n",j,br->GetName());
    TString bnam=br->GetName();
    if(!bnam.Contains("nrun") && !bnam.BeginsWith("err")){
      ttree->SetBranchAddress(bnam,&vectVars[kv]);
      //    ttree->SetBranchAddress(Form("err%s",bnam.Data()),&vectErrs[kv]);
      htr[kv]=CreateHisto(bnam,totEnt);
      ++kv;
    }
  }
  Int_t totHisto=kv;

  for(Int_t je=0; je<totEnt; je++){
    ttree->GetEvent(je);
    printf(" Run %d\n",nrun);
    for(Int_t k=0; k<totHisto; k++){
      htr[k]->SetBinContent(je+1,vectVars[k]);
      if(vectErrs[k]>-900) htr[k]->SetBinError(je+1,vectErrs[k]);
      htr[k]->GetXaxis()->SetBinLabel(je+1,Form("%d",nrun));
    }
  }

  Int_t colors[12]={kRed+1,kRed-7,kOrange+1,kYellow+1,kGreen+1,kGreen,kCyan,kBlue+1,kMagenta,kMagenta+1,kGray+1,1};
  Int_t mstyl[12]={24,21,22,23,20,25,26,27,28,32,33,34};

  TCanvas* cFBfrac=new TCanvas("cFBfrac","Filt bit frac",1200,900);
  Bool_t first=kTRUE;
  TLegend* legfb=new TLegend(0.15,0.8,0.89,0.89);
  legfb->SetNColumns(6);
  for(Int_t k=0; k<totHisto; k++){
    TString hname=htr[k]->GetName();
    if(hname.Contains("RatioFB")){
      htr[k]->GetYaxis()->SetTitle("Ratio to #tracks passing FB0");
      htr[k]->SetMinimum(0);
      htr[k]->SetMaximum(1.2);
      Int_t theFB;
      sscanf(hname.Data(),"RatioFB%dFB0",&theFB);
       if(theFB>=0 && theFB<12){
	htr[k]->SetLineColor(colors[theFB]);
	htr[k]->SetMarkerColor(colors[theFB]);
	htr[k]->SetMarkerStyle(mstyl[theFB]);
	htr[k]->SetMarkerSize(0.7);
	htr[k]->SetLineWidth(2);
	legfb->AddEntry(htr[k],Form("FB%d",theFB),"P")->SetTextColor(colors[theFB]);
	if(!first){ 
	  htr[k]->Draw("same");
	}else{
	  htr[k]->Draw();
	  first=kFALSE;
	}
      }
    }
  }
  legfb->Draw();
  cFBfrac->SaveAs("TrendFracTrackPerFiltBit.png");
  cFBfrac->SaveAs("TrendFracTrackPerFiltBit.root");

  TCanvas* cFBpt=new TCanvas("cFBpt","Filt bit <pt>",1200,900);
  first=kTRUE;
  for(Int_t k=0; k<totHisto; k++){
    TString hname=htr[k]->GetName();
    if(hname.Contains("averPtFB")){
      htr[k]->GetYaxis()->SetTitle("<p_{T}> (GeV/c)");
      htr[k]->SetMinimum(0.3);
      htr[k]->SetMaximum(0.9);
      Int_t theFB;
      sscanf(hname.Data(),"averPtFB%d",&theFB);
       if(theFB>=0 && theFB<12){
	htr[k]->SetLineColor(colors[theFB]);
	htr[k]->SetMarkerColor(colors[theFB]);
	htr[k]->SetMarkerStyle(mstyl[theFB]);
	htr[k]->SetMarkerSize(0.7);
	htr[k]->SetLineWidth(2);
	if(!first){ 
	  htr[k]->Draw("same");
	}else{
	  htr[k]->Draw();
	  first=kFALSE;
	}
      }
    }
  }
  legfb->Draw();
  cFBpt->SaveAs("TrendAverPtPerFiltBit.png");
  cFBpt->SaveAs("TrendAverPtPerFiltBit.root");

  TCanvas* cME=new TCanvas("cME","MatchEff",1200,900);
  cME->Divide(1,3);
  Bool_t firstPlot[3]={kTRUE,kTRUE,kTRUE};
  TLegend* legME= new TLegend(0.2,0.75,0.8,0.89);
  legME->SetMargin(0.15);
  legME->SetNColumns(2);

  for(Int_t k=0; k<totHisto; k++){
    TString hname=htr[k]->GetName();
    if(hname.Contains("MatchEff")){
      htr[k]->GetYaxis()->SetTitle("Matching Efficiency");
      htr[k]->GetYaxis()->SetRangeUser(0.2,1.1);
      TString legTxt="ITS refit";
      if(hname.Contains("SPD")){
	htr[k]->SetMarkerColor(kBlue-7);
	htr[k]->SetLineColor(kBlue-7);	
	htr[k]->SetMarkerStyle(20);
	legTxt.Append(" + SPDany");
      }else{
	htr[k]->SetMarkerColor(1);
	htr[k]->SetLineColor(1);
	htr[k]->SetMarkerStyle(20);	
      }
      Int_t iPtBin=-1;
      if(hname.Contains("EtaNeg")){
	legTxt.Append(", #eta<0");
	htr[k]->SetMarkerStyle(25);	
      }else{
 	legTxt.Append(", #eta>0");
     }
      if(hname.Contains("350")){ 
	iPtBin=0;
	htr[k]->SetTitle("p_{T} = 0.35 GeV/c");
      }
      if(hname.Contains("1000")){
	iPtBin=1;
	htr[k]->SetTitle("p_{T} = 1 GeV/c");
       }
      if(hname.Contains("4000")){ 
	iPtBin=2;
	htr[k]->SetTitle("p_{T} = 4 GeV/c");
      }
      if(iPtBin==0) legME->AddEntry(htr[k],legTxt.Data(),"P");
      cME->cd(iPtBin+1);
      if(!firstPlot[iPtBin]){
	htr[k]->Draw("same");
      }else{
	htr[k]->Draw();
	firstPlot[iPtBin]=kFALSE;
      }
    }
  }
  cME->cd(1);
  legME->Draw();
  cME->SaveAs("TrendMatchingEff.png");
  cME->SaveAs("TrendMatchingEff.root");

  return kTRUE;

}


TH1F* CreateHisto(TString nam, Int_t tote){
  TH1F* h=new TH1F(nam.Data(),Form(" ; run ; %s",nam.Data()),tote,0.5,tote+0.5);
  h->SetStats(0);
  h->SetMarkerStyle(20);
  return h;
}
