#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TTree.h>
#include <TH1.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TLegend.h>
#endif

// Macro to draw the trending information from the output of PlotESDtrackQA.C

TH1F* CreateHisto(TString nam, Int_t tote);

Bool_t DrawTrendingESDTrackQA(TString mergedTrendFile = "trending.root"){
  TFile *fin = TFile::Open(mergedTrendFile.Data());
  if(!fin){
    printf("Cannot open file with ESD track QA trending: %s\n",mergedTrendFile.Data());
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
    if(!bnam.Contains("nrun") && !bnam.BeginsWith("err")){
      ttree->SetBranchAddress(bnam,&vectVars[kv]);
      ttree->SetBranchAddress(Form("err%s",bnam.Data()),&vectErrs[kv]);
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

  TCanvas* cMEp=new TCanvas("cMEp","MatchEffPosEta",900,900);
  cMEp->Divide(1,3);
  Bool_t firstPlotP[3]={kTRUE,kTRUE,kTRUE};
  TCanvas* cMEn=new TCanvas("cMEn","MatchEffNegEta",900,900);
  cMEn->Divide(1,3);
  Bool_t firstPlotN[3]={kTRUE,kTRUE,kTRUE};
  TLegend* legME= new TLegend(0.15,0.17,0.6,0.4);
  legME->SetMargin(0.15);
  for(Int_t k=0; k<totHisto; k++){
    TString hname=htr[k]->GetName();
    if(hname.Contains("MatchEff")){
      htr[k]->GetYaxis()->SetTitle("Matching Efficiency");
      htr[k]->GetYaxis()->SetRangeUser(0.7,1.);
      TString legTxt="ITS refit";
      if(hname.Contains("SPDany")){
	htr[k]->SetMarkerColor(kBlue-7);
	htr[k]->SetLineColor(kBlue-7);	
	legTxt.Append(" + SPDany");
      }
      if(hname.Contains("TOFbc")){ 
	htr[k]->SetMarkerStyle(25);
	legTxt.Append(" + TOFbc=0");
      }	
      Int_t iPtBin=-1;
      TString etaReg="#eta>0";
      if(hname.Contains("EtaNeg")) etaReg="#eta<0";
      if(hname.Contains("350")){ 
	iPtBin=0;
	htr[k]->SetTitle(Form("%s  p_{T} = 0.35 GeV/c",etaReg.Data()));
      }
      if(hname.Contains("1000")){
	iPtBin=1;
	htr[k]->SetTitle(Form("%s  p_{T} = 1 GeV/c",etaReg.Data()));
       }
      if(hname.Contains("4000")){ 
	iPtBin=2;
	htr[k]->SetTitle(Form("%s  p_{T} = 4 GeV/c",etaReg.Data()));
      }
      Int_t draw=-1;
      if(iPtBin>=0){
	if(hname.Contains("EtaPos")){
	  cMEp->cd(iPtBin+1);
	  draw=1;
	  if(firstPlotP[iPtBin]){
	    draw=0;
	    firstPlotP[iPtBin]=kFALSE;
	  }
	  if(iPtBin==0) legME->AddEntry(htr[k],legTxt.Data(),"P");
	}
	if(hname.Contains("EtaNeg")){
	  cMEn->cd(iPtBin+1);
	  draw=1;
	  if(firstPlotN[iPtBin]){
	    draw=0;
	    firstPlotN[iPtBin]=kFALSE;
	  }
	}
      }
      if(draw==0) htr[k]->Draw();
      else if(draw==1) htr[k]->Draw("same");
    }
  }
  cMEp->cd(1);
  legME->Draw();
  cMEn->cd(1);
  legME->Draw();
  cMEp->SaveAs("TrendMatchingEffEtaPos.png");
  cMEp->SaveAs("TrendMatchingEffEtaPos.root");
  cMEn->SaveAs("TrendMatchingEffEtaNeg.png");
  cMEn->SaveAs("TrendMatchingEffEtaNeg.root");

  TCanvas* cPN=new TCanvas("cPN","PosNegCharge",900,900);
  cPN->Divide(1,3);
  Bool_t firstPlotPN[3]={kTRUE,kTRUE,kTRUE};
  TLegend* legPN= new TLegend(0.15,0.17,0.6,0.4);
  legPN->SetMargin(0.15);
  for(Int_t k=0; k<totHisto; k++){
    TString hname=htr[k]->GetName();
    if(hname.Contains("PosNegCharge")){
      htr[k]->GetYaxis()->SetTitle("Positive/negative charge tracks");
      htr[k]->GetYaxis()->SetRangeUser(0.96,1.04);
      TString legTxt="TPC cuts";
      if(hname.Contains("SPDany")){
	htr[k]->SetMarkerColor(kBlue-7);
	htr[k]->SetLineColor(kBlue-7);	
	legTxt.Append(" + ITSrefit + SPDany");
      }
      if(hname.Contains("EtaNeg")){
	htr[k]->SetMarkerStyle(27);
	legTxt.Append(" #eta<0");
      }else{
	legTxt.Append(" #eta>0");
      }
      Int_t iPtBin=-1;
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
      Int_t draw=-1;
      if(iPtBin>=0){
	cPN->cd(iPtBin+1);
	draw=1;
	if(firstPlotPN[iPtBin]){
	  draw=0;
	  firstPlotPN[iPtBin]=kFALSE;
	}
	if(iPtBin==0) legPN->AddEntry(htr[k],legTxt.Data(),"P");
      }
      if(draw==0) htr[k]->Draw();
      else if(draw==1) htr[k]->Draw("same");
    }
  }
  cPN->cd(1);
  legPN->Draw();
  cPN->SaveAs("TrendPosNegCharge.png");
  cPN->SaveAs("TrendPosNegCharge.root");


  TCanvas* cBHpi=new TCanvas("cBHpi","FracBadHypPion",900,900);
  cBHpi->Divide(1,3);
  Bool_t firstPlotPion[3]={kTRUE,kTRUE,kTRUE};
  TCanvas* cBHp=new TCanvas("cBHp","FracBadHypProton",900,900);
  cBHp->Divide(1,3);
  Bool_t firstPlotProt[3]={kTRUE,kTRUE,kTRUE};
  TLegend* legBH= new TLegend(0.15,0.17,0.6,0.4);
  legBH->SetMargin(0.15);
  for(Int_t k=0; k<totHisto; k++){
    TString hname=htr[k]->GetName();
    if(hname.Contains("BadHyp")){
      htr[k]->GetYaxis()->SetTitle("Fraction tracked with wrong mass hypothesis");
      htr[k]->GetYaxis()->SetRangeUser(1e-4,1.);
      TString legTxt="TPC cuts";
      if(hname.Contains("ITSref")){
	htr[k]->SetMarkerColor(kRed+1);
	htr[k]->SetLineColor(kRed+1);	
	htr[k]->SetMarkerStyle(25);
	legTxt.Append(" + ITSrefit");
      }
      if(hname.Contains("SPDany")){
	htr[k]->SetMarkerColor(kBlue-7);
	htr[k]->SetLineColor(kBlue-7);	
	legTxt.Append(" + ITSrefit + SPDany");
      }
      Int_t iPtBin=-1;
      TString hadSpec="Pion";
      if(hname.Contains("Proton")) hadSpec="Proton";
      if(hname.Contains("350")){ 
	iPtBin=0;
	htr[k]->SetTitle(Form("%s  p_{T} = 0.35 GeV/c",hadSpec.Data()));
      }
      if(hname.Contains("600")){
	iPtBin=1;
	htr[k]->SetTitle(Form("%s  p_{T} = 0.6 GeV/c",hadSpec.Data()));
       }
      if(hname.Contains("900")){ 
	iPtBin=2;
	htr[k]->SetTitle(Form("%s  p_{T} = 0.9 GeV/c",hadSpec.Data()));
      }
      Int_t draw=-1;
      if(iPtBin>=0){
	if(hname.Contains("Pion")){
	  cBHpi->cd(iPtBin+1);
	  draw=1;
	  if(firstPlotPion[iPtBin]){
	    draw=0;
	    firstPlotPion[iPtBin]=kFALSE;
	  }
	  if(iPtBin==0) legBH->AddEntry(htr[k],legTxt.Data(),"P");
	}
	if(hname.Contains("Proton")){
	  cBHp->cd(iPtBin+1);
	  draw=1;
	  if(firstPlotProt[iPtBin]){
	    draw=0;
	    firstPlotProt[iPtBin]=kFALSE;
	  }
	}
      }
      if(draw==0){
	gPad->SetLogy();
	htr[k]->Draw();
      }
      else if(draw==1) htr[k]->Draw("same");
    }
  }
  cBHpi->cd(1);
  legBH->Draw();
  cBHp->cd(1);
  legBH->Draw();
  cBHpi->SaveAs("TrendFracBadPhyPion.png");
  cBHpi->SaveAs("TrendFracBadPhyPion.root");
  cBHp->SaveAs("TrendFracBadPhyProton.png");
  cBHp->SaveAs("TrendFracBadPhyProton.root");


  TCanvas* cRes=new TCanvas("cRes","PtResol",900,900);
  cRes->Divide(1,3);
  Bool_t firstPlotRes[3]={kTRUE,kTRUE,kTRUE};
  TLegend* legRes= new TLegend(0.15,0.6,0.6,0.87);
  legRes->SetMargin(0.15);
  for(Int_t k=0; k<totHisto; k++){
    TString hname=htr[k]->GetName();
    if(hname.Contains("PtResol")){
      htr[k]->GetYaxis()->SetTitle("p_{T} resolution from cov. matrix");
      htr[k]->GetYaxis()->SetRangeUser(0.,0.05);
      TString legTxt="";
      if(hname.Contains("Pt10000")){
	htr[k]->SetMarkerColor(kRed);
	htr[k]->SetLineColor(kRed);	
	legTxt.Append("p_{T} = 10 GeV/c");
      }else if(hname.Contains("Pt4000")){
	htr[k]->SetMarkerColor(kOrange+1);
	htr[k]->SetLineColor(kOrange+1);	
	htr[k]->SetMarkerStyle(25);
	legTxt.Append("p_{T} = 4 GeV/c");
      }else if(hname.Contains("Pt1000")){
	htr[k]->SetMarkerColor(kGreen+1);
	htr[k]->SetLineColor(kGreen+1);	
	htr[k]->SetMarkerStyle(33);
	legTxt.Append("p_{T} = 1 GeV/c");
      }else if(hname.Contains("Pt350")){
	htr[k]->SetMarkerColor(kBlue);
	htr[k]->SetLineColor(kBlue);	
	htr[k]->SetMarkerStyle(24);
	legTxt.Append("p_{T} = 0.35 GeV/c");
      }
      Int_t iPad=-1;
      if(hname.Contains("TPC")){ 
	iPad=0;
	htr[k]->SetTitle("TPC cuts");
      }
      if(hname.Contains("ITSref")){
	iPad=1;
	htr[k]->SetTitle("TPC cuts + ITS refit");
       }
      if(hname.Contains("SPDany")){ 
	iPad=2;
	htr[k]->SetTitle("TPC cuts + ITS refit + SPDany");
      }
      Int_t draw=-1;
      if(iPad>=0){
	cRes->cd(iPad+1);
	draw=1;
	if(firstPlotRes[iPad]){
	  draw=0;
	  firstPlotRes[iPad]=kFALSE;
	}
	if(iPad==0) legRes->AddEntry(htr[k],legTxt.Data(),"P");
      }
      if(draw==0) htr[k]->Draw();
      else if(draw==1) htr[k]->Draw("same");
    }
  }
  cRes->cd(1);
  legRes->Draw();
  cRes->SaveAs("TrendPtResol.png");
  cRes->SaveAs("TrendPtResol.root");

  TCanvas* cMass=new TCanvas("cMass","V0mass",900,900);
  cMass->Divide(1,3);
  Bool_t firstPlotMass[3]={kTRUE,kTRUE,kTRUE};
  TCanvas* cSigma=new TCanvas("cSigma","V0sigma",900,900);
  cSigma->Divide(1,3);
  Bool_t firstPlotSigma[3]={kTRUE,kTRUE,kTRUE};

  TLegend* legK0= new TLegend(0.15,0.6,0.6,0.87);
  legK0->SetMargin(0.15);
  TLegend* legL= new TLegend(0.15,0.6,0.6,0.87);
  legL->SetMargin(0.15);
  for(Int_t k=0; k<totHisto; k++){
    TString hname=htr[k]->GetName();
    if(hname.Contains("mass") || hname.Contains("sigma")){
      Int_t iPad=-1;
      TString legTxt="";
      if(hname.Contains("K0")){
	if(hname.Contains("mass")) htr[k]->GetYaxis()->SetRangeUser(0.49,0.51);
	if(hname.Contains("sigma")) htr[k]->GetYaxis()->SetRangeUser(2e-3,6e-3);
	htr[k]->SetTitle("K^{0}_{S}");
	iPad=0;
	if(hname.Contains("0Pt600")){
	  htr[k]->SetMarkerColor(kBlue);
	  htr[k]->SetLineColor(kBlue);	
	  htr[k]->SetMarkerStyle(24);
	  legTxt.Append("0<p_{T}<0.6 GeV/c");
	}else if(hname.Contains("600Pt1000")){
	  htr[k]->SetMarkerColor(kGreen+1);
	  htr[k]->SetLineColor(kGreen+1);	
	  htr[k]->SetMarkerStyle(27);
	  legTxt.Append("0.6<p_{T}<1 GeV/c");
	}else if(hname.Contains("1000Pt3000")){
	  htr[k]->SetMarkerColor(kOrange+1);
	  htr[k]->SetLineColor(kOrange+1);	
	  htr[k]->SetMarkerStyle(25);
	  legTxt.Append("1<p_{T}<3 GeV/c");
	}else if(hname.Contains("5000Pt5200")){
	  htr[k]->SetMarkerColor(kRed);
	  htr[k]->SetLineColor(kRed);
	  htr[k]->SetMarkerStyle(26);
	  legTxt.Append("5<p_{T}<5.2 GeV/c");
	}else{
	  legTxt.Append("All p_{T}");
	}
     }else if(hname.Contains("Lambda")){
	if(hname.Contains("mass")) htr[k]->GetYaxis()->SetRangeUser(1.105,1.125);
	if(hname.Contains("sigma")) htr[k]->GetYaxis()->SetRangeUser(8e-4,25e-4);
	if(hname.Contains("Lambdabar")){
	  iPad=2;
	  htr[k]->SetTitle("#bar{#Lambda}");
	}else{
	  iPad=1;
	  htr[k]->SetTitle("#Lambda");
	}
	if(hname.Contains("0Rad3")){
	  htr[k]->SetMarkerColor(kMagenta+1);
	  htr[k]->SetLineColor(kMagenta+1);	
	  htr[k]->SetMarkerStyle(24);
	  legTxt.Append("0<R<3 cm");
	}else if(hname.Contains("3Rad6")){
	  htr[k]->SetMarkerColor(kBlue+1);
	  htr[k]->SetLineColor(kBlue+1);	
	  htr[k]->SetMarkerStyle(27);
	  legTxt.Append("3<R<6 cm");
	}else if(hname.Contains("8Rad23")){
	  htr[k]->SetMarkerColor(kOrange-7);
	  htr[k]->SetLineColor(kOrange-7);	
	  htr[k]->SetMarkerStyle(25);
	  legTxt.Append("8<R<23 cm");
	}else if(hname.Contains("28Rad43")){
	  htr[k]->SetMarkerColor(kRed+1);
	  htr[k]->SetLineColor(kRed+1);
	  htr[k]->SetMarkerStyle(26);
	  legTxt.Append("28<R<43 cm");
	}else{
	  legTxt.Append("All radii");
	}
      }
      Int_t draw=-1;
      if(iPad>=0){
	if(hname.Contains("mass")){ 
	  htr[k]->GetYaxis()->SetTitle("Inv. Mass peak position (GeV/c^{2})");
 	  cMass->cd(iPad+1);
	  draw=1;
	  if(firstPlotMass[iPad]){
	    draw=0;
	    firstPlotMass[iPad]=kFALSE;
	  }
	  if(iPad==0) legK0->AddEntry(htr[k],legTxt.Data(),"P");
	  if(iPad==1) legL->AddEntry(htr[k],legTxt.Data(),"P");
	}else if(hname.Contains("sigma")){
	  htr[k]->GetYaxis()->SetTitle("Inv. Mass peak width (GeV/c^{2})");
	  cSigma->cd(iPad+1);
	  draw=1;
	  if(firstPlotSigma[iPad]){
	    draw=0;
	    firstPlotSigma[iPad]=kFALSE;
	  }
	}
      }
      if(draw==0) htr[k]->Draw();
      else if(draw==1) htr[k]->Draw("same");
    }
  }
  cMass->cd(1);
  legK0->Draw();
  cMass->cd(2);
  legL->Draw();
  cSigma->cd(1);
  legK0->Draw();
  cSigma->cd(2);
  legL->Draw();
  cMass->SaveAs("TrendV0Mass.png");
  cMass->SaveAs("TrendV0Mass.root");
  cSigma->SaveAs("TrendV0Sigma.png");
  cSigma->SaveAs("TrendV0Sigma.root");

  return kTRUE;
}


TH1F* CreateHisto(TString nam, Int_t tote){
  TH1F* h=new TH1F(nam.Data(),Form(" ; run ; %s",nam.Data()),tote,0.5,tote+0.5);
  h->SetStats(0);
  h->SetMarkerStyle(20);
  return h;
}
