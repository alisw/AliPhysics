#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TTree.h>
#include <TH1.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#endif

// Macro to draw the trending information from the output of PlotAODvertexQA.C

TH1F* CreateHisto(TString nam, Int_t tote);

Bool_t DrawTrendingAODVertexQA(TString mergedTrendFile = "trending.root", 
			       TString treename="trendingVert"){

  TFile *fin = TFile::Open(mergedTrendFile.Data());
  if(!fin){
    printf("Cannot open file with AOD vertex QA trending: %s\n",mergedTrendFile.Data());
    return kFALSE;
  }
  TTree * ttree = (TTree*) fin->Get(treename.Data());
  if (!ttree){
    printf("Trending tree %s not found in file %s\n",treename.Data(),mergedTrendFile.Data());
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
    vectErrs[j]=-999.;
  }
  TH1F** htr=new TH1F*[nVars-1];
  Int_t kv=0;
  for(Int_t j=0; j<nVars; j++){
    TBranch* br=(TBranch*)lb->At(j);
    printf("Branch %d  %s\n",j,br->GetName());
    TString bnam=br->GetName();
    if(!bnam.Contains("nrun") && !bnam.BeginsWith("e")){
      ttree->SetBranchAddress(bnam,&vectVars[kv]);
      if(bnam.Contains("mean") && !bnam.Contains("cont")){
	ttree->SetBranchAddress(Form("e%s",bnam.Data()),&vectErrs[kv]);
      }
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
      else htr[k]->SetBinError(je+1,0.00001);
      htr[k]->GetXaxis()->SetBinLabel(je+1,Form("%d",nrun));
    }
  }

  TCanvas* ctyp=new TCanvas("ctyp","VertexType",1400,900);
  ctyp->Divide(1,2);
  ctyp->cd(2);
  gPad->SetLogy();
  Bool_t first=kTRUE;
  TLegend* leg=new TLegend(0.18,0.18,0.89,0.28);
  leg->SetNColumns(4);
  for(Int_t k=0; k<totHisto; k++){
    TString hname=htr[k]->GetName();
    if(hname.Contains("fTrackV") ||
       hname.Contains("fSPD3D") ||
       hname.Contains("fSPDz") ||
       hname.Contains("fTPC") ||
       hname.Contains("fInvalid")){
      htr[k]->GetYaxis()->SetTitle("Fraction of events per primary vertex type");
      htr[k]->GetYaxis()->SetRangeUser(5e-5,1.1);
      TString legTxt="";
      if(hname.Contains("fTrackV")){
	htr[k]->SetMarkerColor(kBlue);
	htr[k]->SetLineColor(kBlue);	
	htr[k]->SetMarkerStyle(20);
	legTxt="Tracks";
      }else if(hname.Contains("fSPD3D")){
	htr[k]->SetMarkerColor(1);
	htr[k]->SetLineColor(1);
	htr[k]->SetMarkerStyle(21);
	legTxt="SPD-3D";
      }else if(hname.Contains("fSPDz")){
	htr[k]->SetMarkerColor(kGray+1);
	htr[k]->SetLineColor(kGray+1);
	htr[k]->SetMarkerStyle(22);
	legTxt="SPD-Z";
      }else if(hname.Contains("fTPC")){
	htr[k]->SetMarkerColor(kMagenta+1);
	htr[k]->SetLineColor(kMagenta+1);
	htr[k]->SetMarkerStyle(24);
	legTxt="TPC";
      }else if(hname.Contains("fInvalid")){
	htr[k]->SetMarkerColor(2);
	htr[k]->SetLineColor(2);
	htr[k]->SetMarkerStyle(23);
	legTxt="Invalid";
      }
      leg->AddEntry(htr[k],legTxt.Data(),"PL");
      if(!first){
	ctyp->cd(1);
	htr[k]->Draw("same");
	ctyp->cd(2);
	htr[k]->Draw("same");
      }else{
	ctyp->cd(1);
	htr[k]->Draw();
	ctyp->cd(2);
	htr[k]->Draw();
	first=kFALSE;
      }
    }
  }
  ctyp->cd(2);
  leg->Draw();
  ctyp->SaveAs("TrendVertexType.png");
  ctyp->SaveAs("TrendVertexType.root");

  TCanvas* cpos1=new TCanvas("cpos1","VertexPosition",1400,900);
  cpos1->Divide(1,2);
  TCanvas* cpos2=new TCanvas("cpos2","VertexPosition",1400,900);
  cpos2->Divide(1,2);

  for(Int_t k=0; k<totHisto; k++){
    TString hname=htr[k]->GetName();
    Bool_t draw=kFALSE;
    if(hname.Contains("meanx")){
      htr[k]->GetYaxis()->SetTitle("<x_{vert}> (cm)");
      cpos1->cd(1);
      draw=kTRUE;
    }
    else if(hname.Contains("meany")){
      htr[k]->GetYaxis()->SetTitle("<y_{vert}> (cm)");
      cpos1->cd(2);
      draw=kTRUE;
    }
    else if(hname.Contains("meanz")){
      htr[k]->GetYaxis()->SetTitle("<z_{vert}> (cm)");
      cpos2->cd(1);
      draw=kTRUE;
    }
    else if(hname.Contains("meancont")){
      htr[k]->GetYaxis()->SetTitle("<N_{contributors}>");
      cpos2->cd(2);
      draw=kTRUE;
    }
    if(draw) htr[k]->Draw();
  }
  cpos1->SaveAs("TrendVertexPosXY.png");
  cpos1->SaveAs("TrendVertexPosXY.root");
  cpos2->SaveAs("TrendVertexPosZContrib.png");
  cpos2->SaveAs("TrendVertexPosZContrib.root");

  first=kTRUE;
  TCanvas* cpil=new TCanvas("cpil","Pileup",1400,900);
  cpil->SetLogy();
  TLegend* legp=new TLegend(0.18,0.18,0.89,0.28);
  legp->SetNColumns(2);
  legp->SetColumnSeparation(0.2);
  for(Int_t k=0; k<totHisto; k++){
    TString hname=htr[k]->GetName();
    Bool_t draw=kFALSE;
    TString legTxt="";
    if(hname.Contains("pil")){
      htr[k]->GetYaxis()->SetTitle("Fraction of pile-up tagged events");
      htr[k]->GetYaxis()->SetRangeUser(7e-4,1.1);
      if(hname.Contains("fracpilSPD")){
	htr[k]->SetMarkerColor(1);
	htr[k]->SetLineColor(1);
	htr[k]->SetMarkerStyle(21);
	legTxt="SPD all tagged";
      }else if(hname.Contains("fracselpilSPD")){
	htr[k]->SetMarkerColor(1);
	htr[k]->SetLineColor(1);
	htr[k]->SetMarkerStyle(25);
	legTxt="SPD tagged + selected";
      }else if(hname.Contains("fracpilMV")){
	htr[k]->SetMarkerColor(kBlue);
	htr[k]->SetLineColor(kBlue);	
	htr[k]->SetMarkerStyle(20);
	legTxt="Track MultVert all tagged";
      }else if(hname.Contains("fracselpilMV")){
	htr[k]->SetMarkerColor(kBlue);
	htr[k]->SetLineColor(kBlue);	
	htr[k]->SetMarkerStyle(24);
	legTxt="Track MultVert tagged + selected";
      }
      legp->AddEntry(htr[k],legTxt.Data(),"PL");
      if(!first){
	htr[k]->Draw("same");
      }else{
	htr[k]->Draw();
	first=kFALSE;
      }
    }
  }
  legp->Draw();
  cpil->SaveAs("TrendPileupTag.png");
  cpil->SaveAs("TrendPileupTag.root");

  return kTRUE;

}


TH1F* CreateHisto(TString nam, Int_t tote){
  TH1F* h=new TH1F(nam.Data(),Form(" ; run ; %s",nam.Data()),tote,0.5,tote+0.5);
  h->SetStats(0);
  h->SetMarkerStyle(20);
  return h;
}
