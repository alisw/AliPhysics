#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TTree.h>
#include <TH1.h>
#include <TFile.h>
#include <TCanvas.h>
#endif

TH1F* CreateHisto(TString nam, Int_t tote);

Bool_t DrawTrendingPIDQA(TString mergedTrendFile = "trending.root"){
  TFile *fin = TFile::Open(mergedTrendFile.Data());
  if(!fin){
    printf("Cannot open file with PID QA trending: %s\n",mergedTrendFile.Data());
    return kFALSE;
  }
  TTree * ttree = (TTree*) fin->Get("trending");
  if (!ttree){
    printf("Trending tree not found in file %s\n",mergedTrendFile.Data());
    return kFALSE;
  }
  Int_t nrun;
  ttree->SetBranchAddress("nrun",&nrun);
  TObjArray* lb=(TObjArray*)ttree->GetListOfBranches();
  Int_t nVars=lb->GetEntries();
  Float_t* vect=new Float_t[nVars-1];
  TH1F** htr=new TH1F*[nVars-1];
  Int_t kv=0;
  Int_t totEnt=ttree->GetEntries();
  for(Int_t j=0; j<nVars; j++){
    TBranch* br=(TBranch*)lb->At(j);
    printf("Branch %d  %s\n",j,br->GetName());
    TString bnam=br->GetName();
    if(!bnam.Contains("nrun")){
      ttree->SetBranchAddress(bnam,&vect[kv]);
      htr[kv]=CreateHisto(bnam,totEnt);
      ++kv;
    }
  }
  for(Int_t je=0; je<totEnt; je++){
    ttree->GetEvent(je);
    printf(" Run %d\n",nrun);
    for(Int_t k=0; k<nVars-1; k++){
      htr[k]->SetBinContent(je+1,vect[k]);
      htr[k]->GetXaxis()->SetBinLabel(je+1,Form("%d",nrun));
    }
  }


  TCanvas* c=0x0;
  for(Int_t k=0; k<nVars-1; k++){
    if(k%4==0){
      if(c) c->SaveAs(Form("PIDqaTrend%d.png",k/4));
      c=new TCanvas(Form("c%d",k/4),Form("c%d",k/4),1500,600);
      c->Divide(2,2);
    }
    Int_t thePad=(k%4)+1;
    c->cd(thePad);
    htr[k]->Draw("");
    htr[k]->Draw("psame");
  }
  if(c) c->SaveAs(Form("PIDqaTrend%d.png",(nVars-1)/4));
  delete [] vect;

}


TH1F* CreateHisto(TString nam, Int_t tote){
  TH1F* h=new TH1F(nam.Data(),Form(" ; run ; %s",nam.Data()),tote,0.5,tote+0.5);
  h->SetStats(0);
  h->SetMarkerStyle(20);
  return h;
}
