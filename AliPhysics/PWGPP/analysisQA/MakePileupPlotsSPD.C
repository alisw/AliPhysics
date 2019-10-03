#include <TObjString.h>
#include <TString.h>
#include <iostream>
#include <fstream>
#include <string>

Double_t fgaussian(Double_t* x, Double_t* par){
  if(x[0]>-1. && x[0]<1.){
    TF1::RejectPoint();
    return 0;
  }
  Double_t nsig=(x[0]-par[1])/par[2];
  return par[0]*TMath::Exp(-0.5*nsig*nsig);
}


void MakePileupPlotsSPD(const char *filein = "AnalysisResults.root",
		     TString suffix = "SPDcont5",
		     TString imsuffix="jpg",
                     const char *outfile = "output.root"){
 
  gInterpreter->ExecuteMacro("$ALICE_PHYSICS/..//src/PWGHF/vertexingHF/macros/LoadLibraries.C");

  TString algo="SPD";

  TString text=Form("Algorithm: %s",algo.Data());
  if(suffix=="SPDcont5" && algo=="SPD") text.Append("-5contr");
  TLatex* tal=new TLatex(0.22,0.84,text.Data());
  tal->SetNDC();
  

  gStyle->SetOptTitle(0);

  TFile *fil = TFile::Open(filein, "read");
      if(!fil) {
            Printf("FATAL: File doesn't exist");
            return;
      }
 
 if(suffix=="SPDcont5")
 { const Char_t * outDir="TaskPileupFPrinoSPDcont5";}
  else {const Char_t* outDir="TaskPileupFPrinoSPDcont3";}
  CreateDir(outDir);
  CreateDir(Form("%s",outDir));

 // TFile* fil=new TFile("AnalysisResults.root");
  TDirectoryFile* df=(TDirectoryFile*)fil->Get(Form("CheckPileup_%s",suffix.Data()));
      if(!df) {
	Printf(Form("FATAL: TDirectory CheckPileup_%s doesn't exist",suffix.Data()));
            return;
      }
  
  AliCounterCollection* cnt=(AliCounterCollection*)df->Get("counters");
  TString runString=cnt->GetKeyWords("RUN");
  cout<<runString.Data()<<endl;
  TObjArray* runArray=runString.Tokenize(",");

  const Int_t nRuns=runArray->GetEntries();
  Int_t* runList=new Int_t[nRuns];
  for(Int_t i=0; i<nRuns; i++){
    TObjString* os=(TObjString*)runArray->At(i);
    TString numb=os->GetString();
    runList[i]=numb.Atoi();
    cout<< runList[i]<<endl;
  }

  TH1F* hEvTrig=new TH1F("hEvTrig","",nRuns,0.5,nRuns+0.5);
  TH1F* hEvPhysSel=new TH1F("hEvPhysSel","",nRuns,0.5,nRuns+0.5);
  TH1F* hEvSPDVert=new TH1F("hEvSPDVert","",nRuns,0.5,nRuns+0.5);
  TH1F* hEvSPDPil=new TH1F("hEvSPDPil","",nRuns,0.5,nRuns+0.5);
  TH1F* hEvMVPil=new TH1F("hEvMVPil","",nRuns,0.5,nRuns+0.5);

  TH1F* hFracPileupSPD=new TH1F("hFracPileupSPD","",nRuns,0.5,nRuns+0.5);
  TH1F* hFracPileupMV=new TH1F("hFracPileupMV","",nRuns,0.5,nRuns+0.5);

  for(Int_t irun=0; irun<nRuns; irun++){
    Int_t trig=cnt->GetSum(Form("Event:Triggered/RUN:%d",runList[irun]));
    Int_t ps=cnt->GetSum(Form("Event:PhysSel/RUN:%d",runList[irun]));
    Int_t spdvert=cnt->GetSum(Form("Event:SPDVert/RUN:%d",runList[irun]));
    Int_t pileupMV=cnt->GetSum(Form("Event:PileupMV/RUN:%d",runList[irun]));
    Int_t pileupSPD=cnt->GetSum(Form("Event:PileupSPD/RUN:%d",runList[irun]));
    printf("%d MV=%d(%f) SPD=%d(%f)\n",spdvert,pileupMV,(Double_t)pileupMV/(Double_t)spdvert,pileupSPD,(Double_t)pileupSPD/(Double_t)spdvert);
    hEvTrig->SetBinContent(irun+1,trig);
    hEvTrig->GetXaxis()->SetBinLabel(irun+1,Form("%d",runList[irun]));
    hEvPhysSel->SetBinContent(irun+1,ps);
    hEvPhysSel->GetXaxis()->SetBinLabel(irun+1,Form("%d",runList[irun]));
    hEvSPDVert->SetBinContent(irun+1,spdvert);
    hEvSPDVert->GetXaxis()->SetBinLabel(irun+1,Form("%d",runList[irun]));
    hEvSPDPil->SetBinContent(irun+1,pileupSPD);
    hEvSPDPil->GetXaxis()->SetBinLabel(irun+1,Form("%d",runList[irun]));
    hEvMVPil->SetBinContent(irun+1,pileupMV);
    hEvMVPil->GetXaxis()->SetBinLabel(irun+1,Form("%d",runList[irun]));
    hFracPileupSPD->SetBinContent(irun+1,(Double_t)pileupSPD/(Double_t)spdvert);
    hFracPileupSPD->SetBinError(irun+1,TMath::Sqrt((Double_t)pileupSPD)/(Double_t)spdvert);
    hFracPileupSPD->GetXaxis()->SetBinLabel(irun+1,Form("%d",runList[irun]));
    hFracPileupMV->SetBinContent(irun+1,(Double_t)pileupMV/(Double_t)spdvert);
    hFracPileupMV->SetBinError(irun+1,TMath::Sqrt((Double_t)pileupMV)/(Double_t)spdvert);
    hFracPileupMV->GetXaxis()->SetBinLabel(irun+1,Form("%d",runList[irun]));
  }

  TCanvas* c0=new TCanvas("c0","Run dependence",1400,700);
  c0->Divide(2,1);
  c0->cd(1);
  gPad->SetLogy();
  hEvTrig->SetStats(0);
  hEvTrig->SetMinimum(hEvSPDPil->GetMinimum()/100.);
  hEvTrig->GetXaxis()->SetTitle("Run number");
  hEvTrig->GetYaxis()->SetTitle("Events");
  hEvTrig->GetYaxis()->SetTitleOffset(1.3);
  hEvTrig->Draw("e");
  hEvPhysSel->SetLineColor(6);
  hEvPhysSel->Draw("esame");
  hEvSPDVert->SetLineColor(kGreen+2);
  hEvSPDVert->Draw("esame");
  hEvSPDPil->SetLineColor(2);
  hEvSPDPil->Draw("esame");
  hEvMVPil->SetLineColor(4);
  hEvMVPil->Draw("esame");
  TLegend* leg=new TLegend(0.15,0.15,0.6,0.35);
  leg->AddEntry(hEvTrig,"Triggered.","L")->SetTextColor(hEvTrig->GetLineColor());
  leg->AddEntry(hEvPhysSel,"Phys. Sel.","L")->SetTextColor(hEvPhysSel->GetLineColor());
  leg->AddEntry(hEvSPDVert,"With SPD vert","L")->SetTextColor(hEvSPDVert->GetLineColor());
  leg->AddEntry(hEvSPDPil,"SPD pileup","L")->SetTextColor(hEvSPDPil->GetLineColor());
  leg->AddEntry(hEvMVPil,"MultiVert pileup","L")->SetTextColor(hEvMVPil->GetLineColor());
  leg->Draw();
  c0->cd(2);
  gPad->SetLeftMargin(0.15);
  gPad->SetRightMargin(0.05);
  hFracPileupSPD->GetXaxis()->SetTitle("Run number");
  hFracPileupSPD->GetYaxis()->SetTitle("Fraction of events with pileup");
  hFracPileupSPD->GetYaxis()->SetTitleOffset(1.8);
  hFracPileupSPD->SetStats(0);
  hFracPileupSPD->SetMinimum(0);
  hFracPileupSPD->SetMaximum(1.2*TMath::Max(hFracPileupSPD->GetMaximum(),hFracPileupMV->GetMaximum()));
  hFracPileupSPD->SetLineColor(2);
  hFracPileupSPD->Draw("e");
  hFracPileupMV->SetLineColor(4);
  hFracPileupMV->Draw("esame");
  c0->SaveAs(Form("%s/fig_rundependence.%s",outDir,imsuffix.Data()));   

  gStyle->SetOptFit(1111);
  TList* lv=(TList*)df->Get(Form("clistPrimaryV_%s",suffix.Data()));
  lv->ls();
  TList* lp=(TList*)df->Get(Form("clistPileup%s_%s",algo.Data(),suffix.Data()));
  lp->ls();
  TH1F* hNOfPileupVert=(TH1F*)lp->FindObject(Form("hNOfPileupVert%s",algo.Data()));
  TCanvas* cn=new TCanvas("cn","Npilvert",800,800);
  cn->SetLogy();
  hNOfPileupVert->GetXaxis()->SetTitle("Number of pileup vertices");
  hNOfPileupVert->GetYaxis()->SetTitle("Events");
  hNOfPileupVert->GetYaxis()->SetTitleOffset(1.2);
  hNOfPileupVert->Draw();
  tal->Draw();
  cn->SaveAs(Form("%s/fig_Numpileupvtx.%s",outDir,imsuffix.Data()));

  TH1F* hZVertSPD=(TH1F*)lv->FindObject("hZVertSPD");
  TH1F* hZDiffTaggingPil=(TH1F*)lp->FindObject(Form("hZDiffTaggingPilZDiamcut%s",algo.Data()));
  hZDiffTaggingPil->Rebin(2);
  TCanvas* c1=new TCanvas("c1","Z diff",1400,700);
  c1->Divide(2,1);
  c1->cd(1);
  gPad->SetLeftMargin(0.15);
  gPad->SetRightMargin(0.05);
  hZVertSPD->GetXaxis()->SetTitle("z_{SPDvertex} (cm)");
  hZVertSPD->GetYaxis()->SetTitle("Events");
  hZVertSPD->GetYaxis()->SetTitleOffset(1.8);
  hZVertSPD->Fit("gaus","","",-10.,10.);
  TF1* fdiam=(TF1*)hZVertSPD->GetListOfFunctions()->FindObject("gaus");
  Double_t sig=fdiam->GetParameter(2);
  c1->cd(2);
  gPad->SetLeftMargin(0.15);
  gPad->SetRightMargin(0.05);
  TF1* fgausdz=new TF1("fgausdz",fgaussian,-10.,10.,3);
  fgausdz->SetParameter(0,2*hZDiffTaggingPil->Integral(hZDiffTaggingPil->FindBin(1.),hZDiffTaggingPil->FindBin(10.)));
  fgausdz->SetParameter(1,0.);
  fgausdz->SetParameter(2,fdiam->GetParameter(2)*TMath::Sqrt(2.));
  hZDiffTaggingPil->GetXaxis()->SetTitle("#Deltaz_{vertex} = z_{primary}-z_{pile} (cm)");
  hZDiffTaggingPil->GetYaxis()->SetTitle("Events");
  hZDiffTaggingPil->GetYaxis()->SetTitleOffset(1.8);
  hZDiffTaggingPil->Fit(fgausdz,"R");
  Double_t sigd=fgausdz->GetParameter(2);

  TF1* fgausdiff=new TF1("fgausdiff","gaus",-10.,10.);
  for(Int_t ip=0; ip<3; ip++) fgausdiff->SetParameter(ip,fgausdz->GetParameter(ip));
  // fgausdiff->SetParameter(0,fgausdz->GetParameter(0));
  fgausdiff->SetLineColor(2);
  fgausdiff->Draw("same");
  tal->Draw();
  //  hZDiffTaggingPil->Fit(fgausdiff);
  //  fgausdz->Draw("same");
  
  printf("\n-------------------------------------------------\n");
  printf("DeltaZ Gaussian fit sigma = %f (expected = %f)\n",sigd,sig*TMath::Sqrt(2));
  Double_t intFunc=fgausdiff->Integral(-10.,10.)/hZDiffTaggingPil->GetBinWidth(1);
  Double_t hentries=hZDiffTaggingPil->Integral(hZDiffTaggingPil->GetXaxis()->FindBin(-9.999),hZDiffTaggingPil->GetXaxis()->FindBin(9.999));
  printf("Integral of gaussian fit=%f   Integral of histo = %f\n",intFunc,hentries);
  printf("Efficiency = %f\n",hentries/intFunc);
  Double_t expEff=fdiam->Integral(0.8,10.)/fdiam->Integral(0.,10.);
  printf("Expected efficiency from integral of diamond profile = %f\n",expEff);
  printf("-------------------------------------------------\n");
 
  TH1F* hContribPrimVertPil=(TH1F*)lp->FindObject(Form("hContribPrimVertPil%s",algo.Data()));
  hContribPrimVertPil->SetLineColor(2);
  TH1F* hContribPrimVertNoPil=(TH1F*)lp->FindObject(Form("hContribPrimVertNoPil%s",algo.Data()));
  hContribPrimVertNoPil->SetLineColor(1);
  TH1F* hContribTaggingPil=(TH1F*)lp->FindObject(Form("hContribTaggingPil%s",algo.Data()));
  hContribTaggingPil->SetLineColor(4);
  hContribTaggingPil->SetLineWidth(2);
  TH1F* hContribFirstPil=(TH1F*)lp->FindObject(Form("hContribFirstPil%s",algo.Data()));
  hContribFirstPil->SetLineColor(kMagenta+2);
  TH1F* hContribSecondPil=(TH1F*)lp->FindObject(Form("hContribSecondPil%s",algo.Data()));
  hContribSecondPil->SetLineColor(6);
  c1->SaveAs(Form("%s/fig_zvertex.%s",outDir,imsuffix.Data()));   

  TCanvas* c3=new TCanvas("c3","",800,800);
  c3->SetLogy();
  hContribPrimVertNoPil->GetXaxis()->SetRangeUser(-0.5,120.5);
  hContribPrimVertNoPil->GetXaxis()->SetTitle("Vertex contributors");
  hContribPrimVertNoPil->GetYaxis()->SetTitle("Events");
  hContribPrimVertNoPil->GetYaxis()->SetTitleOffset(1.2);
  hContribPrimVertNoPil->Draw();
  c3->Update();
  TPaveStats* st=(TPaveStats*)hContribPrimVertNoPil->GetListOfFunctions()->FindObject("stats");
  st->SetTextColor(hContribPrimVertNoPil->GetLineColor());
  st->SetY1NDC(0.81);
  st->SetY2NDC(0.96);
  c3->Modified();
  hContribPrimVertPil->Draw("sames");
  c3->Update();
  st=(TPaveStats*)hContribPrimVertPil->GetListOfFunctions()->FindObject("stats");
  st->SetTextColor(hContribPrimVertPil->GetLineColor());
  st->SetY1NDC(0.65);
  st->SetY2NDC(0.80);
  c3->Modified();
  hContribTaggingPil->Draw("sames");
  c3->Update();
  st=(TPaveStats*)hContribTaggingPil->GetListOfFunctions()->FindObject("stats");
  st->SetTextColor(hContribTaggingPil->GetLineColor());
  st->SetY1NDC(0.49);
  st->SetY2NDC(0.64);
  c3->Modified();
  hContribFirstPil->Draw("sames");
  c3->Update();
  st=(TPaveStats*)hContribFirstPil->GetListOfFunctions()->FindObject("stats");
  st->SetTextColor(hContribFirstPil->GetLineColor());
  st->SetY1NDC(0.33);
  st->SetY2NDC(0.48);
  c3->Modified();
  hContribSecondPil->Draw("sames");
  c3->Update();
  st=(TPaveStats*)hContribSecondPil->GetListOfFunctions()->FindObject("stats");
  st->SetTextColor(hContribSecondPil->GetLineColor());
  st->SetY1NDC(0.17);
  st->SetY2NDC(0.32);
  c3->Modified();
  tal->Draw();
 
  TH1F* hNCL1Pil=(TH1F*)lp->FindObject(Form("hNCL1Pil%s",algo.Data()));
  hNCL1Pil->SetLineColor(2);
  TH1F* hNCL1NoPil=(TH1F*)lp->FindObject(Form("hNCL1NoPil%s",algo.Data()));
  hNCL1NoPil->SetLineColor(1);
  TH1F* hFracPilCL1=new TH1F(Form("hFracPilCL1%s",algo.Data()),"",hNCL1Pil->GetNbinsX(),hNCL1Pil->GetXaxis()->GetXmin(),hNCL1Pil->GetXaxis()->GetXmax());
  hFracPilCL1->SetStats(0);
  for(Int_t i=1; i<=hNCL1Pil->GetNbinsX(); i++){
    Double_t x=hNCL1Pil->GetBinContent(i);
    Double_t y=hNCL1NoPil->GetBinContent(i);
    if(x>0 || y>0){
      Double_t r=x/(x+y);
      Double_t er=TMath::Sqrt(r*(1-r)/(x+y));
      hFracPilCL1->SetBinContent(i,r);
      hFracPilCL1->SetBinError(i,er);
    }
  }
  c3->SaveAs(Form("%s/fig_contributions.%s",outDir,imsuffix.Data())); 

  TCanvas* c4=new TCanvas("c4","",1500,750);
  c4->Divide(2,1);
  c4->cd(1);
  gPad->SetLogy();
  //  hNCL1NoPil->GetXaxis()->SetRangeUser(-0.5,120.5);
  hNCL1NoPil->GetXaxis()->SetTitle("CL1 multiplicity");
  hNCL1NoPil->GetYaxis()->SetTitle("Events");
  hNCL1NoPil->GetYaxis()->SetTitleOffset(1.2);
  hNCL1NoPil->Draw();
  gPad->Update();
  st=(TPaveStats*)hNCL1NoPil->GetListOfFunctions()->FindObject("stats");
  st->SetTextColor(hNCL1NoPil->GetLineColor());
  st->SetY1NDC(0.81);
  st->SetY2NDC(0.96);
  gPad->Modified();
  hNCL1Pil->Draw("sames");
  gPad->Update();
  st=(TPaveStats*)hNCL1Pil->GetListOfFunctions()->FindObject("stats");
  st->SetTextColor(hNCL1Pil->GetLineColor());
  st->SetY1NDC(0.65);
  st->SetY2NDC(0.80);
  gPad->Modified();
  c4->cd(2);
  hFracPilCL1->GetXaxis()->SetTitle("CL1 multiplicity");
  hFracPilCL1->GetYaxis()->SetTitle("Fraction of events tagged as pileup");
  hFracPilCL1->GetYaxis()->SetTitleOffset(1.2);
  hFracPilCL1->Draw();
  tal->Draw();


  TH1F* hNtracklPil=(TH1F*)lp->FindObject(Form("hNtracklPil%s",algo.Data()));
  hNtracklPil->SetLineColor(2);
  TH1F* hNtracklNoPil=(TH1F*)lp->FindObject(Form("hNtracklNoPil%s",algo.Data()));
  hNtracklNoPil->SetLineColor(1);
  TH1F* hFracPilTrkl=new TH1F(Form("hFracPilTrkl%s",algo.Data()),"",hNtracklPil->GetNbinsX(),hNtracklPil->GetXaxis()->GetXmin(),hNtracklPil->GetXaxis()->GetXmax());
  hFracPilTrkl->SetStats(0);
  for(Int_t i=1; i<=hNtracklPil->GetNbinsX(); i++){
    Double_t x=hNtracklPil->GetBinContent(i);
    Double_t y=hNtracklNoPil->GetBinContent(i);
    if(x>0 || y>0){
      Double_t r=x/(x+y);
      Double_t er=TMath::Sqrt(r*(1-r)/(x+y));
      hFracPilTrkl->SetBinContent(i,r);
      hFracPilTrkl->SetBinError(i,er);
    }
  }
  c4->SaveAs(Form("%s/fig_CL1.%s",outDir,imsuffix.Data())); 

  TCanvas* c5=new TCanvas("c5","",1500,750);
  c5->Divide(2,1);
  c5->cd(1);
  gPad->SetLogy();
  hNtracklNoPil->GetXaxis()->SetRangeUser(-0.5,120.5);
  hNtracklNoPil->GetXaxis()->SetTitle("N_{tracklets}");
  hNtracklNoPil->GetYaxis()->SetTitle("Events");
  hNtracklNoPil->GetYaxis()->SetTitleOffset(1.2);
  hNtracklNoPil->Draw();
  gPad->Update();
  st=(TPaveStats*)hNtracklNoPil->GetListOfFunctions()->FindObject("stats");
  st->SetTextColor(hNtracklNoPil->GetLineColor());
  st->SetY1NDC(0.81);
  st->SetY2NDC(0.96);
  gPad->Modified();
  hNtracklPil->Draw("sames");
  gPad->Update();
  st=(TPaveStats*)hNtracklPil->GetListOfFunctions()->FindObject("stats");
  st->SetTextColor(hNtracklPil->GetLineColor());
  st->SetY1NDC(0.65);
  st->SetY2NDC(0.80);
  gPad->Modified();
  c5->cd(2);
  hFracPilTrkl->GetXaxis()->SetRangeUser(-0.5,120.5);
  hFracPilTrkl->GetXaxis()->SetTitle("N_{tracklets}");
  hFracPilTrkl->GetYaxis()->SetTitle("Fraction of events tagged as pileup");
  hFracPilTrkl->GetYaxis()->SetTitleOffset(1.2);
  hFracPilTrkl->Draw();
  tal->Draw();
  c5->SaveAs(Form("%s/fig_tracklets.%s",outDir,imsuffix.Data())); 

}
//--------------------------------
void CreateDir(const Char_t* dirName)
{
   TString pwd(gSystem->pwd());
   gSystem->cd(pwd.Data());
     
   if(gSystem->cd(dirName)) {
   gSystem->cd(pwd.Data());
   } else {
    gSystem->mkdir(dirName, kTRUE); // kTRUE means recursive
    }
}      
