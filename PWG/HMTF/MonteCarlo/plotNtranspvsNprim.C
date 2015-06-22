#if !defined (__CINT__) || (defined(__MAKECINT__))
#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "AliHeader.h"
#include "AliGenEventHeader.h"
#include "AliStack.h"
#include "TParticle.h"
#include "TH2F.h"
#include "TH1F.h"

#endif

void plotNtranspvsNprim(Int_t estimatescaling=1) {

  gSystem->Load("libpythia6_4_25.so");
  gSystem->Load("libpythia6.so");

  TFile* f = new TFile("1Mdummypy8/galice.root");
  TTree* t = (TTree*) f->Get("TE");
  AliHeader* header = new AliHeader;
  AliGenPythiaEventHeader * headPy  = 0;
  AliGenDPMjetEventHeader * headPho = 0;

  t->SetBranchAddress("Header",&header);
  TH2F* hTranspvsPrim = new TH2F("hTranspvsPrim","Total number of transported particles vs total primaries in stack",2000,0.,2000.,1000,0.,1000.);
  TH1F* hAllTransp = new TH1F("hAllTransp","Number of transported particles distribution",2000,0.,2000.);
  TH1F* hIntegrAllTransp = new TH1F("hIntegrAllTransp","Integral histogram of Number of transported particles distribution",2000,0.,2000.);
  TH1F* hScaleAllTransp = new TH1F("hScaleAllTransp","Scaled Number of transported particles distribution",2000,0.,2000.);

  hTranspvsPrim->Sumw2();
  hAllTransp->Sumw2();
  hIntegrAllTransp->Sumw2();
  hScaleAllTransp->Sumw2();

  for(Int_t i=0; i<t->GetEntries(); i++) {

    t->GetEntry(i);
    AliStack* stack = new AliStack*;
    stack=header->Stack();
    AliGenEventHeader * htmp = header->GenEventHeader();
    Double_t weight=htmp->EventWeight();
//   if(weight>40000) cout<<weight<<endl;
//    headPy =  (AliGenPythiaEventHeader*) htmp;

    Int_t nprim = stack->GetNprimary();
    Int_t ntransp = stack->GetNtransported();
//    TParticle* p = (TParticle*) stack->GetCurrentTrack();
//    if (p->Eta()>0.5) ||  (p->Eta()<-0.5)) cout<<"no"; 
    hTranspvsPrim->Fill(nprim,ntransp);
//        if(headPy->ProcessType() != 92 && headPy->ProcessType() != 93 && headPy->ProcessType() != 94) 
// if(headPy->ProcessType() == 68)
//        {
    hAllTransp->Fill(ntransp);
    if (estimatescaling==0) hScaleAllTransp->Fill(ntransp,weight);
//  cout<<headPy->ProcessType()<<endl;// != 92 && headPy->ProcessType() != 93 && headPy->ProcessType() != 94)\n";
//        }
//    cout<<ntransp<<" "<<nprim<<endl;
  }
  hAllTransp->Rebin(1);
  hAllTransp->Scale(1./hAllTransp->Integral());
  if (estimatescaling==0)  hScaleAllTransp->Scale(1./hScaleAllTransp->Integral());

  for (Int_t i=0; i<hAllTransp->GetNbinsX(); i++) {
  if (estimatescaling==0) hIntegrAllTransp->SetBinContent(i,hAllTransp->Integral(i,hAllTransp->GetNbinsX()));
  }

  Int_t nThresholds=9;
  Int_t* Thresholds = new Int_t[nThresholds];
  Double_t Scaling=1.;
  Double_t* Scalingarr = new Double_t[nThresholds];
//  Double_t limit=1.;
  for (Int_t i=0; i<nThresholds; i++) {
if (estimatescaling==1) Scaling = Scaling/10.;
//    limit=limit-Scaling;
//    Thresholds[i]=hIntegrAllTransp->GetBinLowEdge(hIntegrAllTransp->FindLastBinAbove(limit));//Scaling));

   Thresholds[i]=50*(i+1);
   if (estimatescaling==1) Scaling=hAllTransp->Integral((i==0)?0:hAllTransp->FindBin(Thresholds[i-1]),((i<nThresholds-1)?hAllTransp->FindBin(Thresholds[i]):hAllTransp->GetNbinsX()));//hAllTransp->Integral(i,hAllTransp->GetNbinsX());//hIntegrAllTransp->GetBinContent(Thresholds[i]);
//hIntegrAllTransp->SetBinContent(i,hScaleAllTransp->Integral((i==0)?0:hScaleAllTransp->FindBin(Thresholds[i-1]),((i<nThresholds-1)?hScaleAllTransp->FindBin(Thresholds[i]):hScaleAllTransp->GetNbinsX())));
    for (Int_t k=(i==0)?0:hAllTransp->FindBin(Thresholds[i-1]); k<((i<nThresholds-1)?hAllTransp->FindBin(Thresholds[i]):hAllTransp->GetNbinsX()); k++) {
      if (estimatescaling==1) hIntegrAllTransp->SetBinContent(k,hAllTransp->Integral(k,hAllTransp->GetNbinsX()));
//      else hIntegrAllTransp->SetBinContent(k,hAllTransp->Integral(k,hAllTransp->GetNbinsX()));
      if (estimatescaling==1) hScaleAllTransp->SetBinContent(k,(1./Scaling)*hAllTransp->GetBinContent(k));
    }

//    std::cout<<"Threshold for "<<Scaling<<" most central events: "<<Thresholds[i]<<"; Nevents="<<hAllTransp->Integral((i==0)?0:hAllTransp->FindBin(Thresholds[i-1]),((i<nThresholds-1)?hAllTransp->FindBin(Thresholds[i]):hAllTransp->GetNbinsX()))<<"; NeventsInt="<<hIntegrAllTransp->Integral((i==0)?0:hIntegrAllTransp->FindBin(Thresholds[i-1]),((i<nThresholds-1)?hIntegrAllTransp->FindBin(Thresholds[i]):hIntegrAllTransp->GetNbinsX()))<<"; NeventsScale="<<hScaleAllTransp->Integral((i==0)?0:hScaleAllTransp->FindBin(Thresholds[i-1]),((i<nThresholds-1)?hScaleAllTransp->FindBin(Thresholds[i]):hScaleAllTransp->GetNbinsX()))<<std::endl;
 std::cout<<"Threshold="<<Thresholds[i]<<"; Scaling="<<Scaling<<std::endl;
  Scalingarr[i]=Scaling;
  }

  std::cout<<"threshold=\"";
  for (Int_t i=0; i<nThresholds; i++) {
    std::cout<<Thresholds[i]<<((i<nThresholds-1)?",":"\"\n");
  }

  std::cout<<"scaling=\"";
  for (Int_t i=0; i<nThresholds; i++) {
//    std::cout<<Scalingarr[i]*1000000.<<((i<nThresholds-1)?",":"\"\n");
    Int_t Scaleint = (Int_t)(0.5+(Scalingarr[i]/Scalingarr[nThresholds-1]));
    Scalingarr[i]=Scaleint;
    std::cout<<Scaleint<<((i<nThresholds-1)?",":"\"\n");
  }

  //hTranspvsPrim->Print();

  TCanvas* c = new TCanvas("Canvas1","Canvas1");

  c->Divide(2,2);

  c->cd(1);
  hTranspvsPrim->SetOption("colz");
  hTranspvsPrim->Draw();
  hTranspvsPrim->SetDrawOption("colz");
  c->cd(2);
  hAllTransp->Draw();

  c->cd(3);
  hIntegrAllTransp->Draw();

  c->cd(4);
//  if (estimatescaling==1) 
hScaleAllTransp->Draw();

}
