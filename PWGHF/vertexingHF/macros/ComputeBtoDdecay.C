#if !defined(__CINT__) || defined(__MAKECINT__)
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TPythia6.h"
#include "AliTPythia8.h"
#include "AliDecayerPythia8.h"
#include "TPythia6Decayer.h"
#include "TCanvas.h"
#include "TClonesArray.h"
#include "TSystem.h"
#include "TRandom.h"
#include "TParticle.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TLorentzVector.h"
#include "TDatabasePDG.h"
#include "Riostream.h"
#include "TObjArray.h"
#include "TMCParticle.h"
#include "TStyle.h"
#include "TTree.h"
#include "TRandom3.h"
#endif

TH1D* ReadFONLL(TString filename, Int_t option=0);



void ComputeBtoDdecay(Int_t nGener=10000000,
		      Int_t pythiaver=8,
		      TString fileNameFONLL="FONLL-Bhadron-dsdpt-sqrts5020-50MeVbins.txt",
		      Int_t opt4ff=0,
		      Int_t optForNorm=0,
		      Bool_t writeTree=kFALSE){

  Int_t pdgD0=421;
  Int_t pdgDp=411;
  Int_t pdgDs=431;
  Int_t pdgLc=4122;

  Double_t fracB[4]={0.401,0.401,0.105,0.093};
  if(opt4ff==0){
    // ppbar fractions
    fracB[0]=0.340;
    fracB[1]=0.340;
    fracB[2]=0.101;
    fracB[3]=0.219;
  }else if(opt4ff==1){
    // e+e- fractions
    fracB[0]=0.412;
    fracB[1]=0.412;
    fracB[2]=0.088;
    fracB[3]=0.088;    
  }
  
  TVirtualMCDecayer* pdec=0x0;
  
  if(pythiaver==6){
    gSystem->Load("liblhapdf.so");      // Parton density functions
    gSystem->Load("libEGPythia6.so");   // TGenerator interface
    gSystem->Load("libpythia6.so");     // Pythia
    //    gSystem->Load("libAliPythia6.so");  // ALICE specific implementations
    pdec=new TPythia6Decayer();
  }else{
    gSystem->Load("liblhapdf.so");      // Parton density functions
    gSystem->Load("libpythia8.so");
    gSystem->Load("libAliPythia8.so");
    gSystem->Setenv("PYTHIA8DATA", gSystem->ExpandPathName("$ALICE_ROOT/PYTHIA8/pythia8/xmldoc"));
    gSystem->Setenv("LHAPDF",      gSystem->ExpandPathName("$ALICE_ROOT/LHAPDF"));
    gSystem->Setenv("LHAPATH",     gSystem->ExpandPathName("$ALICE_ROOT/LHAPDF/PDFsets"));
    pdec=new AliDecayerPythia8();
  }

  TDatabasePDG* db=TDatabasePDG::Instance();
  pdec->Init();

  TH1D * hBptDistr = ReadFONLL(fileNameFONLL.Data(),0);
  Double_t xsecb=0;
  for(Int_t i=1; i<=hBptDistr->GetNbinsX(); i++){
    xsecb+=(hBptDistr->GetBinContent(i)*hBptDistr->GetBinWidth(i));
  }
  printf("B hadron production cross section = %f pb\n",xsecb);
  hBptDistr->Scale(1.e-6); // convert to ub

  TH1F* hD0Origin=new TH1F("hD0Origin","D0 mother ; ; Entries",4,-0.5,3.5);
  hD0Origin->GetXaxis()->SetBinLabel(1,"B+");
  hD0Origin->GetXaxis()->SetBinLabel(2,"B0");
  hD0Origin->GetXaxis()->SetBinLabel(3,"B_s");
  hD0Origin->GetXaxis()->SetBinLabel(4,"Lb");
  TH1F* hDpOrigin=new TH1F("hDpOrigin","D+ mother ; ; Entries",4,-0.5,3.5);
  hDpOrigin->GetXaxis()->SetBinLabel(1,"B+");
  hDpOrigin->GetXaxis()->SetBinLabel(2,"B0");
  hDpOrigin->GetXaxis()->SetBinLabel(3,"B_s");
  hDpOrigin->GetXaxis()->SetBinLabel(4,"Lb");
  TH1F* hDsOrigin=new TH1F("hDsOrigin","D_{s} mother ; ; Entries",4,-0.5,3.5);
  hDsOrigin->GetXaxis()->SetBinLabel(1,"B+");
  hDsOrigin->GetXaxis()->SetBinLabel(2,"B0");
  hDsOrigin->GetXaxis()->SetBinLabel(3,"B_s");
  hDsOrigin->GetXaxis()->SetBinLabel(4,"Lb");
  TH1F* hLcOrigin=new TH1F("hLcOrigin","#Lambda_{c} mother ; ; Entries",4,-0.5,3.5);
  hLcOrigin->GetXaxis()->SetBinLabel(1,"B+");
  hLcOrigin->GetXaxis()->SetBinLabel(2,"B0");
  hLcOrigin->GetXaxis()->SetBinLabel(3,"B_s");
  hLcOrigin->GetXaxis()->SetBinLabel(4,"Lb");

  TH1F* hB0dau=new TH1F("hB0dau","",5,-1.5,3.5);
  hB0dau->GetXaxis()->SetBinLabel(1,"All B0");
  hB0dau->GetXaxis()->SetBinLabel(2,"D0");
  hB0dau->GetXaxis()->SetBinLabel(3,"D+");
  hB0dau->GetXaxis()->SetBinLabel(4,"D_s");
  hB0dau->GetXaxis()->SetBinLabel(5,"Lc");
  TH1F* hBpdau=new TH1F("hBplusdau","",5,-1.5,3.5);
  hBpdau->GetXaxis()->SetBinLabel(1,"All B+");
  hBpdau->GetXaxis()->SetBinLabel(2,"D0");
  hBpdau->GetXaxis()->SetBinLabel(3,"D+");
  hBpdau->GetXaxis()->SetBinLabel(4,"D_s");
  hBpdau->GetXaxis()->SetBinLabel(5,"Lc");
  TH1F* hBsdau=new TH1F("hBsdau","",5,-1.5,3.5);
  hBsdau->GetXaxis()->SetBinLabel(1,"All Bs");
  hBsdau->GetXaxis()->SetBinLabel(2,"D0");
  hBsdau->GetXaxis()->SetBinLabel(3,"D+");
  hBsdau->GetXaxis()->SetBinLabel(4,"D_s");
  hBsdau->GetXaxis()->SetBinLabel(5,"Lc");
  TH1F* hLbdau=new TH1F("hLbdau","",5,-1.5,3.5);
  hLbdau->GetXaxis()->SetBinLabel(1,"All Lb");
  hLbdau->GetXaxis()->SetBinLabel(2,"D0");
  hLbdau->GetXaxis()->SetBinLabel(3,"D+");
  hLbdau->GetXaxis()->SetBinLabel(4,"D_s");
  hLbdau->GetXaxis()->SetBinLabel(5,"Lc");
  
  TH1F* hD0pt=new TH1F("hD0pt","",1001,0.,50.05);
  TH1F* hDppt=new TH1F("hDpluspt","",1001,0.,50.05);
  TH1F* hDspt=new TH1F("hDspt","",1001,0.,50.05);
  TH1F* hLcpt=new TH1F("hLcpt","",1001,0.,50.05);

  TH2D* hD0PtByOrigin=new TH2D("hD0PtByOrigin","",4,-0.5,3.5,1001,0.,50.05);
  hD0PtByOrigin->GetXaxis()->SetBinLabel(1,"B+");
  hD0PtByOrigin->GetXaxis()->SetBinLabel(2,"B0");
  hD0PtByOrigin->GetXaxis()->SetBinLabel(3,"B_s");
  hD0PtByOrigin->GetXaxis()->SetBinLabel(4,"Lb");
  TH2D* hDpPtByOrigin=new TH2D("hDplusPtByOrigin","",4,-0.5,3.5,1001,0.,50.05);
  hDpPtByOrigin->GetXaxis()->SetBinLabel(1,"B+");
  hDpPtByOrigin->GetXaxis()->SetBinLabel(2,"B0");
  hDpPtByOrigin->GetXaxis()->SetBinLabel(3,"B_s");
  hDpPtByOrigin->GetXaxis()->SetBinLabel(4,"Lb");
  TH2D* hDsPtByOrigin=new TH2D("hDsPtByOrigin","",4,-0.5,3.5,1001,0.,50.05);
  hDsPtByOrigin->GetXaxis()->SetBinLabel(1,"B+");
  hDsPtByOrigin->GetXaxis()->SetBinLabel(2,"B0");
  hDsPtByOrigin->GetXaxis()->SetBinLabel(3,"B_s");
  hDsPtByOrigin->GetXaxis()->SetBinLabel(4,"Lb");
  TH2D* hLcPtByOrigin=new TH2D("hLcPtByOrigin","",4,-0.5,3.5,1001,0.,50.05);
  hLcPtByOrigin->GetXaxis()->SetBinLabel(1,"B+");
  hLcPtByOrigin->GetXaxis()->SetBinLabel(2,"B0");
  hLcPtByOrigin->GetXaxis()->SetBinLabel(3,"B_s");
  hLcPtByOrigin->GetXaxis()->SetBinLabel(4,"Lb");

  TH2F* hD0PtVsB0pt=new TH2F("hD0PtVsB0pt"," ; p_{T}(B0) ; p_{T}(D0)",100,0.,50.,100.,0.,50.);
  TH2F* hD0PtVsBppt=new TH2F("hD0PtVsBpluspt"," ; p_{T}(B+) ; p_{T}(D0)",100,0.,50.,100.,0.,50.);
  TH2F* hD0PtVsBspt=new TH2F("hD0PtVsBspt"," ; p_{T}(Bs) ; p_{T}(D0)",100,0.,50.,100.,0.,50.);
  TH2F* hD0PtVsLbpt=new TH2F("hD0PtVsLbpt"," ; p_{T}(Lb) ; p_{T}(D0)",100,0.,50.,100.,0.,50.);
  TH2F* hDpPtVsB0pt=new TH2F("hDplusPtVsB0pt"," ; p_{T}(B0) ; p_{T}(D+)",100,0.,50.,100.,0.,50.);
  TH2F* hDpPtVsBppt=new TH2F("hDplusPtVsBpluspt"," ; p_{T}(B+) ; p_{T}(D+)",100,0.,50.,100.,0.,50.);
  TH2F* hDpPtVsBspt=new TH2F("hDplusPtVsBspt"," ; p_{T}(Bs) ; p_{T}(D+)",100,0.,50.,100.,0.,50.);
  TH2F* hDpPtVsLbpt=new TH2F("hDplusPtVsLbpt"," ; p_{T}(Lb) ; p_{T}(D+)",100,0.,50.,100.,0.,50.);
  TH2F* hDsPtVsB0pt=new TH2F("hDsPtVsB0pt"," ; p_{T}(B0) ; p_{T}(Ds)",100,0.,50.,100.,0.,50.);
  TH2F* hDsPtVsBppt=new TH2F("hDsPtVsBpluspt"," ; p_{T}(B+) ; p_{T}(Ds)",100,0.,50.,100.,0.,50.);
  TH2F* hDsPtVsBspt=new TH2F("hDsPtVsBspt"," ; p_{T}(Bs) ; p_{T}(Ds)",100,0.,50.,100.,0.,50.);
  TH2F* hDsPtVsLbpt=new TH2F("hDsPtVsLbpt"," ; p_{T}(Lb) ; p_{T}(Ds)",100,0.,50.,100.,0.,50.);
  TH2F* hLcPtVsB0pt=new TH2F("hLcPtVsB0pt"," ; p_{T}(B0) ; p_{T}(Lc)",100,0.,50.,100.,0.,50.);
  TH2F* hLcPtVsBppt=new TH2F("hLcPtVsBpluspt"," ; p_{T}(B+) ; p_{T}(Lc)",100,0.,50.,100.,0.,50.);
  TH2F* hLcPtVsBspt=new TH2F("hLcPtVsBspt"," ; p_{T}(Bs) ; p_{T}(Lc)",100,0.,50.,100.,0.,50.);
  TH2F* hLcPtVsLbpt=new TH2F("hLcPtVsLbpt"," ; p_{T}(Lb) ; p_{T}(Lc)",100,0.,50.,100.,0.,50.);


  TTree* fTreeDecays = 0x0;
  Int_t pdgB = -9999;
  Double_t ptB = -1.;
  Double_t pB = -1.;
  Double_t yB = -1.;
  vector<float> arrptD;
  vector<float> arrpD;
  vector<float> arryD;
  vector<int> arrpdgD;
  Double_t norm = xsecb;

  if(writeTree){
    fTreeDecays = new TTree("fTreeDecays", "fTreeDecays");
    fTreeDecays->Branch("pdgB", &pdgB);
    fTreeDecays->Branch("ptB", &ptB);
    fTreeDecays->Branch("pB", &pB);
    fTreeDecays->Branch("yB", &yB);
    fTreeDecays->Branch("ptD", &arrptD);
    fTreeDecays->Branch("pD", &arrpD);
    fTreeDecays->Branch("yD", &arryD);
    fTreeDecays->Branch("pdgD", &arrpdgD);
    fTreeDecays->Branch("norm", &norm);
  }
  
  TRandom3* gener=new TRandom3(0);
  TClonesArray *array = new TClonesArray("TParticle",100);
  TLorentzVector* vec=new TLorentzVector();

  Double_t countB=0;
  
  for(Int_t itry=0; itry<nGener; itry++){
    if(itry%10000==0) printf("Particle %d\n",itry);

    Int_t iBin=1;
    Double_t value=gener->Rndm();
    TH1F* hdautofill=0x0;
    TH2F* hptD0tofill=0x0;
    TH2F* hptDptofill=0x0;
    TH2F* hptDstofill=0x0;
    TH2F* hptLctofill=0x0;
    if(value<fracB[0]){ 
      pdgB=511;
      iBin=1;
      hdautofill=hBpdau;
      hptD0tofill=hD0PtVsBppt;
      hptDptofill=hDpPtVsBppt;
      hptDstofill=hDsPtVsBppt;
      hptLctofill=hLcPtVsBppt;
   }else if(value<(fracB[0]+fracB[1])){ 
      pdgB=521;
      iBin=2;
      hdautofill=hB0dau;
      hptD0tofill=hD0PtVsB0pt;
      hptDptofill=hDpPtVsB0pt;
      hptDstofill=hDsPtVsB0pt;
      hptLctofill=hLcPtVsB0pt;
    }else if(value<(fracB[0]+fracB[1]+fracB[2])){ 
      pdgB=531;
      iBin=3;
      hdautofill=hBsdau;
      hptD0tofill=hD0PtVsBspt;
      hptDptofill=hDpPtVsBspt;
      hptDstofill=hDsPtVsBspt;
      hptLctofill=hLcPtVsBspt;
    }else{
      pdgB=5122;
      iBin=4;
      hdautofill=hLbdau;
      hptD0tofill=hD0PtVsLbpt;
      hptDptofill=hDpPtVsLbpt;
      hptDstofill=hDsPtVsLbpt;
      hptLctofill=hLcPtVsLbpt;
    }
    hdautofill->Fill(-1.);
    
    Double_t mass=db->GetParticle(pdgB)->Mass();
    ptB=hBptDistr->GetRandom();
    Double_t phiB=gener->Rndm()*2*TMath::Pi();
    yB=gener->Rndm()*2.-1.; // flat in -1<y<1
    Double_t px=ptB*TMath::Cos(phiB);
    Double_t py=ptB*TMath::Sin(phiB);
    Double_t mt=TMath::Sqrt(mass*mass+ptB*ptB);
    Double_t pz=mt*TMath::SinH(yB);
    pB=TMath::Sqrt(ptB*ptB+pz*pz);
    Double_t E=TMath::Sqrt(mass*mass+pB*pB);
    vec->SetPxPyPzE(px,py,pz,E);
    pdec->Decay(pdgB,vec);
    if(optForNorm==0) countB+=1.;
    else if(optForNorm==1 && TMath::Abs(yB)<0.5) countB+=1.;
    Int_t nentries = pdec->ImportParticles(array);
    //    TParticle* bmes=(TParticle*)array->At(0);

    for(int j=0; j<nentries; j++){
      TParticle * part = (TParticle*)array->At(j);
      Int_t pdgdau=TMath::Abs(part->GetPdgCode());
      Double_t ptD=-999;
      Double_t yD=-999;
      if(pdgdau==pdgD0 || pdgdau==pdgDp || pdgdau==pdgDs || pdgdau==pdgLc){
	ptD=part->Pt();
	yD=part->Y();
	arrptD.push_back(ptD);
	arrpD.push_back(part->P());
	arryD.push_back(yD);
	arrpdgD.push_back(pdgdau);
	if(optForNorm==0 || (optForNorm==1 && TMath::Abs(yD)<0.5)){
	  if(pdgdau==pdgD0){
	    hD0Origin->Fill(iBin-1);
	    hD0pt->Fill(ptD);
	    hdautofill->Fill(0.);
	    hD0PtByOrigin->Fill(iBin-1,ptD);
	    hptD0tofill->Fill(ptB,ptD);
	  }else if(pdgdau==pdgDp){
	    hDpOrigin->Fill(iBin-1);
	    hDppt->Fill(ptD);
	    hdautofill->Fill(1.);
	    hDpPtByOrigin->Fill(iBin-1,ptD);
	    hptDptofill->Fill(ptB,ptD);
	  }else if(pdgdau==pdgDs){
	    hDsOrigin->Fill(iBin-1);
	    hDspt->Fill(ptD);
	    hdautofill->Fill(2.);
	    hDsPtByOrigin->Fill(iBin-1,ptD);
	    hptDstofill->Fill(ptB,ptD);
	  }else if(pdgdau==pdgLc){
	    hLcOrigin->Fill(iBin-1);
	    hLcpt->Fill(ptD);
	    hdautofill->Fill(3.);
	    hLcPtByOrigin->Fill(iBin-1,ptD);
	    hptLctofill->Fill(ptB,ptD);
	  }
	}
      }
    }
    if(arrptD.size() == 0){
      arrptD.push_back(-1);
      arrpD.push_back(-1);
      arryD.push_back(-999);
      arrpdgD.push_back(-1);
    }
    if(fTreeDecays) fTreeDecays->Fill();
    arrptD.clear();
    arrpD.clear();
    arryD.clear();
    arrpdgD.clear();
    array->Clear();
  }
  
  delete vec;
  delete array;
  delete gener;
  
  hD0pt->Scale(xsecb/1e6/countB/hD0pt->GetBinWidth(1));
  hDppt->Scale(xsecb/1e6/countB/hDppt->GetBinWidth(1));
  hDspt->Scale(xsecb/1e6/countB/hDspt->GetBinWidth(1));
  hLcpt->Scale(xsecb/1e6/countB/hLcpt->GetBinWidth(1));
  hD0PtByOrigin->Scale(xsecb/1e6/countB/hD0PtByOrigin->GetYaxis()->GetBinWidth(1));
  hDpPtByOrigin->Scale(xsecb/1e6/countB/hDpPtByOrigin->GetYaxis()->GetBinWidth(1));
  hDsPtByOrigin->Scale(xsecb/1e6/countB/hDsPtByOrigin->GetYaxis()->GetBinWidth(1));
  hLcPtByOrigin->Scale(xsecb/1e6/countB/hLcPtByOrigin->GetYaxis()->GetBinWidth(1));

  printf("Cross sections for B = %f ub \n",xsecb/1e6);
  Double_t xsecD0=0;
  for(Int_t i=1; i<=hD0pt->GetNbinsX(); i++){
    xsecD0+=(hD0pt->GetBinContent(i)*hD0pt->GetBinWidth(i));
  }
  Double_t xsecDp=0;
  for(Int_t i=1; i<=hDppt->GetNbinsX(); i++){
    xsecDp+=(hDppt->GetBinContent(i)*hDppt->GetBinWidth(i));
  }
  Double_t xsecDs=0;
  for(Int_t i=1; i<=hDspt->GetNbinsX(); i++){
    xsecDs+=(hDspt->GetBinContent(i)*hDspt->GetBinWidth(i));
  }
  Double_t xsecLc=0;
  for(Int_t i=1; i<=hLcpt->GetNbinsX(); i++){
    xsecLc+=(hLcpt->GetBinContent(i)*hLcpt->GetBinWidth(i));
  }
  printf("Cross sections for D0 = %f ub  D+ = %f ub   Ds=%f ub  Lc=%f ub\n",xsecD0,xsecDp,xsecDs,xsecLc);

  hD0pt->SetStats(0);
  hDppt->SetStats(0);
  hDspt->SetStats(0);
  hLcpt->SetStats(0);
  hBptDistr->SetStats(0);
  hBptDistr->GetYaxis()->SetTitle("d#sigma/dp_{T} (#mub/GeV)");
  hBptDistr->GetXaxis()->SetTitle("p_{T} (GeV)");
  hD0pt->GetYaxis()->SetTitle("d#sigma/dp_{T} (#mub/GeV)");
  hD0pt->GetXaxis()->SetTitle("p_{T} (GeV)");
  hDppt->GetYaxis()->SetTitle("d#sigma/dp_{T} (#mub/GeV)");
  hDppt->GetXaxis()->SetTitle("p_{T} (GeV)");
  hDspt->GetYaxis()->SetTitle("d#sigma/dp_{T} (#mub/GeV)");
  hDspt->GetXaxis()->SetTitle("p_{T} (GeV)");
  hLcpt->GetYaxis()->SetTitle("d#sigma/dp_{T} (#mub/GeV)");
  hLcpt->GetXaxis()->SetTitle("p_{T} (GeV)");

  TH1F* hDsKKpipt=(TH1F*)hDspt->Clone("hDsKKpipt");
  hDsKKpipt->Scale(0.0227);
  hDsKKpipt->GetYaxis()->SetTitle("d#sigma/dp_{T}xBR (#mub/GeV)");
    
  TCanvas* c1=new TCanvas("c1","B mother",1500,900);
  c1->Divide(2,2);
  c1->cd(1);
  hD0Origin->Draw();
  c1->cd(2);
  hDpOrigin->Draw();
  c1->cd(3);
  hDsOrigin->Draw();
  c1->cd(4);
  hLcOrigin->Draw();

  
  TCanvas* c2=new TCanvas("c2","PtD vs PtB",1500,1000);
  c2->Divide(4,4);
  c2->cd(1);
  gPad->SetLogz();
  hD0PtVsB0pt->Draw("colz");
  c2->cd(2);
  gPad->SetLogz();
  hD0PtVsBppt->Draw("colz");
  c2->cd(3);
  gPad->SetLogz();
  hD0PtVsBspt->Draw("colz");
  c2->cd(4);
  gPad->SetLogz();
  hD0PtVsLbpt->Draw("colz");
  c2->cd(5);
  gPad->SetLogz();
  hDpPtVsB0pt->Draw("colz");
  c2->cd(6);
  gPad->SetLogz();
  hDpPtVsBppt->Draw("colz");
  c2->cd(7);
  gPad->SetLogz();
  hDpPtVsBspt->Draw("colz");
  c2->cd(8);
  gPad->SetLogz();
  hDpPtVsLbpt->Draw("colz");
  c2->cd(9);
  gPad->SetLogz();
  hDsPtVsB0pt->Draw("colz");
  c2->cd(10);
  gPad->SetLogz();
  hDsPtVsBppt->Draw("colz");
  c2->cd(11);
  gPad->SetLogz();
  hDsPtVsBspt->Draw("colz");
  c2->cd(12);
  gPad->SetLogz();
  hDsPtVsLbpt->Draw("colz");
  c2->cd(13);
  gPad->SetLogz();
  hLcPtVsB0pt->Draw("colz");
  c2->cd(14);
  gPad->SetLogz();
  hLcPtVsBppt->Draw("colz");
  c2->cd(15);
  gPad->SetLogz();
  hLcPtVsBspt->Draw("colz");
  c2->cd(16);
  gPad->SetLogz();
  hLcPtVsLbpt->Draw("colz");
  c2->SaveAs(Form("DecayKine_Pythia%d.png",pythiaver));
  
  TCanvas* c3=new TCanvas("c3","pt-diff xsec",900,800);
  gPad->SetLogy();
  hBptDistr->Draw();
  hD0pt->SetLineColor(2);
  hD0pt->Draw("same");
  hDppt->SetLineColor(4);
  hDppt->Draw("same");
  hDspt->SetLineColor(kGreen+1);
  hDspt->Draw("same");
  hLcpt->SetLineColor(kMagenta+1);
  hLcpt->Draw("same");
  TLegend* leg=new TLegend(0.5,0.6,0.89,0.89);
  leg->AddEntry(hBptDistr,"B hadron, FONLL","L")->SetTextColor(hBptDistr->GetLineColor());
  leg->AddEntry(hD0pt,Form("D^{0} #leftarrowB (FONLL+PYTHIA%d)",pythiaver),"L")->SetTextColor(hD0pt->GetLineColor());
  leg->AddEntry(hDppt,Form("D^{+} #leftarrowB (FONLL+PYTHIA%d)",pythiaver),"L")->SetTextColor(hDppt->GetLineColor());
  leg->AddEntry(hDspt,Form("D_{s}^{+} #leftarrowB (FONLL+PYTHIA%d)",pythiaver),"L")->SetTextColor(hDspt->GetLineColor());
  leg->AddEntry(hLcpt,Form("#Lambda_{c}^{+} #leftarrowB (FONLL+PYTHIA%d)",pythiaver),"L")->SetTextColor(hLcpt->GetLineColor());
  leg->Draw();
  c3->SaveAs(Form("XsecBandDfromB_FONLLPythia%d.png",pythiaver));

  TString outfilnam=Form("DfromBtest_FONLLPythia%d",pythiaver);
  if(opt4ff==0) outfilnam.Append("_FFppbar");
  else if(opt4ff==1) outfilnam.Append("_FFee");
  else outfilnam.Append("_FFold");
  if(optForNorm==1) outfilnam.Append("_yDcut");
  outfilnam.Append(".root");
  TFile* outfil=new TFile(outfilnam.Data(),"recreate");
  hBptDistr->Write();
  hD0pt->Write();
  hDppt->Write();
  hDspt->Write();
  hLcpt->Write();
  hDsKKpipt->Write();
  hB0dau->Write();
  hBpdau->Write();
  hBsdau->Write();
  hLbdau->Write();
  hD0PtByOrigin->Write();
  hDpPtByOrigin->Write();
  hDsPtByOrigin->Write();
  hLcPtByOrigin->Write();
  hD0PtVsB0pt->Write();
  hD0PtVsBppt->Write();
  hD0PtVsBspt->Write();
  hD0PtVsLbpt->Write();
  hDpPtVsB0pt->Write();
  hDpPtVsBppt->Write();
  hDpPtVsBspt->Write();
  hDpPtVsLbpt->Write();
  hDsPtVsB0pt->Write();
  hDsPtVsBppt->Write();
  hDsPtVsBspt->Write();
  hDsPtVsLbpt->Write();
  hLcPtVsB0pt->Write();
  hLcPtVsBppt->Write();
  hLcPtVsBspt->Write();
  hLcPtVsLbpt->Write();
  if(fTreeDecays) fTreeDecays->Write();
  outfil->Close();
}


//----------------------------------------------
TH1D* ReadFONLL(TString filename, Int_t option){
  FILE* infil=fopen(filename.Data(),"r");
  Char_t line[200];
  Char_t* rc;
  for(Int_t il=0; il<16; il++){
    rc=fgets(line,200,infil);
    if(strstr(line,"central")) break;
  }
  Float_t pt,csc,csmin,csmax,dum;
  Double_t ptmin=999,ptmax=0.;
  Int_t iPt=0;
  Double_t x[2000],y[2000];
  Bool_t ok;
  while(!feof(infil)){
    ok=fscanf(infil,"%f %f %f %f",&pt,&csc,&csmin,&csmax);
    for(Int_t i=0; i<12;i++) ok=fscanf(infil,"%f",&dum);
    if(feof(infil)) break;
    if(pt==0.) continue;
    if(pt<ptmin) ptmin=pt;
    if(pt>ptmax) ptmax=pt;
    x[iPt]=pt;
    if(option==1) y[iPt]=csmin;
    else if(option==2) y[iPt]=csmax;
    else  y[iPt]=csc;
    iPt++;
  }
  fclose(infil);
  Double_t binw=(ptmax-ptmin)/(iPt-1);
  TString hname="hfonllB";
  TH1D* hfonll=new TH1D(hname.Data(),"",iPt,ptmin-0.5*binw,ptmax+0.5*binw);
  printf("nPoints=%d  binw=%f ptmin=%f  ptmax=%f  firstbin=%f lastbin=%f\n",iPt,binw,hfonll->GetXaxis()->GetXmin(),hfonll->GetXaxis()->GetXmax(),x[0],x[iPt-1]);
  for(Int_t iBin=0; iBin<iPt; iBin++){
    hfonll->SetBinContent(iBin+1,y[iBin]);
  }
  return hfonll;

}
