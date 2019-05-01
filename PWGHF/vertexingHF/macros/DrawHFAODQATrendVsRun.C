#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TTree.h>
#include <TH1.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#endif

// Macro to draw the trending information from the output of readMCPerform.C

TH1F* CreateHisto(TString nam, Int_t tote);

Bool_t DrawHFAODQATrendVsRun(TString mergedTrendFile = "trending.root", 
			   TString treename="trendingHF"){

  TFile *fin = TFile::Open(mergedTrendFile.Data());
  if(!fin){
    printf("Cannot open file with HF QA trending: %s\n",mergedTrendFile.Data());
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
    if(!bnam.Contains("nrun")){
      ttree->SetBranchAddress(bnam,&vectVars[kv]);
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

  TCanvas* ccan=new TCanvas("ccan","Candidates",1400,900);
  ccan->Divide(3,2);
  for(Int_t k=0; k<totHisto; k++){
    TString hname=htr[k]->GetName();
    Bool_t draw=kFALSE;
    if(hname.Contains("nDzeroCandperEv")){
      htr[k]->GetYaxis()->SetTitle("N. D^{0} candidates / event");
      ccan->cd(1);
      draw=kTRUE;
    }
    else if(hname.Contains("nDplusCandperEv")){
      htr[k]->GetYaxis()->SetTitle("N. D^{+} candidates / event");
      ccan->cd(2);
      draw=kTRUE;
    }
    else if(hname.Contains("nDsCandperEv")){
      htr[k]->GetYaxis()->SetTitle("N. D_{s}^{+} candidates / event");
      ccan->cd(4);
      draw=kTRUE;
    }
    else if(hname.Contains("nDstarCandperEv")){
      htr[k]->GetYaxis()->SetTitle("N. D*^{+} candidates / event");
      ccan->cd(3);
      draw=kTRUE;
    }
    else if(hname.Contains("nLcCandperEv")){
      htr[k]->GetYaxis()->SetTitle("N. #Lambda_{c}^{+} candidates / event");
      ccan->cd(5);
      draw=kTRUE;
    }
    if(draw) htr[k]->Draw();
  }
  ccan->SaveAs("TrendHFCandidates.png");
  ccan->SaveAs("TrendHFCandidates.root");

  return kTRUE;
}

TH1F* CreateHisto(TString nam, Int_t tote){
  TH1F* h=new TH1F(nam.Data(),Form(" ; run ; %s",nam.Data()),tote,0.5,tote+0.5);
  h->SetStats(0);
  h->SetMarkerStyle(20);
  return h;
}
