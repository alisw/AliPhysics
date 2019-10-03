#if !defined (__CINT__) || defined (__CLING__)
#include <Riostream.h>
#include <TMath.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TColor.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TGaxis.h>
#include <TSpline.h>

#endif

//macro for the computation of q2 percentiles using splines 
//author: Fabrizio Grosa, fabrizio.grosa@to.infn.it ,INFN Torino

//*******************************************//
//                                           //
//    Main Function: SplineFitCumulativeq2   //
//                                           //
//*******************************************//

using namespace std;

enum {kTaskHFvn, kTaskCharmHadronvn};

const TString infilename = "./AnalysisResults.root";
const TString indirname = "PWGHF_D2H_HFvn_Dplus_3050VZERO_EvShapeEP";
const TString inlistname = "coutputvnDplus_3050VZERO_EvShapeEP";

const TString outfilename = "./Splines_3050.root";

Int_t SplineFitCumulativeq2(int task = kTaskCharmHadronvn) {
    
  gStyle->SetOptStat(0);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadBottomMargin(0.12);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  
  const Int_t nDet=6;
  TString detnames[nDet] = {"TPC","TPCNegEta","TPCPosEta","VZERO","VZEROA","VZEROC"};

  //input file
  TFile* infile = TFile::Open(infilename.Data(),"READ");
  if(!infile) {return 1;}
  TDirectoryFile* dir = (TDirectoryFile*)infile->Get(indirname.Data());
  if(!dir) {cerr << "TDirectoryFile name wrong! Exit." << endl; return 2;}
  TList* list = (TList*)dir->Get(inlistname.Data());
  if(!list) {cerr << "TList name wrong! Exit." << endl; return 3;}
  
  TString histoname = "";
  if(task==kTaskHFvn) 
    histoname = "hq2vsCentr";
  else if(task==kTaskCharmHadronvn) {
    histoname = "fHistqnVsCentr";
    for(Int_t iDet=3; iDet<nDet; iDet++) 
      detnames[iDet].ReplaceAll("VZERO","V0");
  }

  TH2F* hq2VsCentr[nDet];
  for(Int_t iDet=0; iDet<nDet; iDet++) {
    hq2VsCentr[iDet] = (TH2F*)list->FindObject(Form("%s%s",histoname.Data(),detnames[iDet].Data()));
    if(!hq2VsCentr[iDet]) {cerr << Form("%s%s not found! Exit.",histoname.Data(),detnames[iDet].Data()) << endl; return 4;}
  } 

  //compute the splines in centrality bins for each detector
  const Int_t nCentrBins = hq2VsCentr[0]->GetXaxis()->GetNbins();
  const Int_t nq2Bins = hq2VsCentr[0]->GetYaxis()->GetNbins();
  TH1F* hq2Int[nDet][nCentrBins];
  TSpline3* sq2Int[nDet][nCentrBins];
  TList* splinelist[nDet];

  cout << "\n\n" << endl;
  for(Int_t iDet=0; iDet<nDet; iDet++) {
    splinelist[iDet] = new TList();
    splinelist[iDet]->SetOwner(0);
    
    cout << Form("Fitting splines for qn%s",detnames[iDet].Data()) << endl;
    for(Int_t iCentr=0; iCentr<nCentrBins; iCentr++) {
      TH1F* hq2_CentrBin = (TH1F*)hq2VsCentr[iDet]->ProjectionY("hq2_CentrBin",iCentr+1,iCentr+1);
      Double_t q2integral=0;
      hq2Int[iDet][iCentr] = new TH1F(Form("hq2Int_%d_%s",iCentr,detnames[iDet].Data()),Form(";#it{q}_{2}^{%s};#it{q}_{2}^{%s} normalised integral",detnames[iDet].Data(),detnames[iDet].Data()),nq2Bins,0.,15.);
      for(Int_t iq2=0; iq2<nq2Bins; iq2++) {
        q2integral += hq2_CentrBin->GetBinContent(iq2+1);
        Double_t q2 = hq2_CentrBin->GetBinCenter(iq2+1);
        hq2Int[iDet][iCentr]->SetBinContent(iq2+1,q2integral/hq2_CentrBin->Integral()*100);
      }
      Double_t min_centr = hq2VsCentr[iDet]->GetXaxis()->GetBinLowEdge(iCentr+1);
      Double_t max_centr = hq2VsCentr[iDet]->GetXaxis()->GetBinLowEdge(iCentr+1)+hq2VsCentr[iDet]->GetXaxis()->GetBinWidth(iCentr+1);
      sq2Int[iDet][iCentr] = new TSpline3(hq2Int[iDet][iCentr]);
      sq2Int[iDet][iCentr]->SetName(Form("sq2Int_centr_%0.f_%0.f",min_centr,max_centr));
      splinelist[iDet]->Add(sq2Int[iDet][iCentr]);
    }
  }
    
  //save splines in a file  
  TFile outfile(outfilename.Data(),"recreate");
  for(Int_t iDet=0; iDet<nDet; iDet++) {
    splinelist[iDet]->Write(Form("SplineListq2%s",detnames[iDet].Data()),1);
  }
  outfile.Close();

  //Plot some examples for each detector
  TCanvas* cq2Int[nDet];
  TLegend* leg = new TLegend(0.6,0.2,0.8,0.4);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.035);

  for(Int_t iDet=0; iDet<nDet; iDet++) {
    cq2Int[iDet] = new TCanvas(Form("cq2Int%s",detnames[iDet].Data()),"",800,800);
    
    if(iDet==0) {
      leg->AddEntry(hq2Int[iDet][0],Form("Centr %0.f-%0.f(%%)",hq2VsCentr[iDet]->GetXaxis()->GetBinLowEdge(1),hq2VsCentr[iDet]->GetXaxis()->GetBinLowEdge(1)+hq2VsCentr[iDet]->GetXaxis()->GetBinWidth(1)),"p");
      leg->AddEntry(hq2Int[iDet][10],Form("Centr %0.f-%0.f(%%)",hq2VsCentr[iDet]->GetXaxis()->GetBinLowEdge(11),hq2VsCentr[iDet]->GetXaxis()->GetBinLowEdge(11)+hq2VsCentr[iDet]->GetXaxis()->GetBinWidth(11)),"p");
      leg->AddEntry(hq2Int[iDet][19],Form("Centr %0.f-%0.f(%%)",hq2VsCentr[iDet]->GetXaxis()->GetBinLowEdge(20),hq2VsCentr[iDet]->GetXaxis()->GetBinLowEdge(20)+hq2VsCentr[iDet]->GetXaxis()->GetBinWidth(20)),"p");
    }

    hq2Int[iDet][19]->SetMarkerStyle(kFullCircle);
    hq2Int[iDet][19]->SetMarkerColor(kGreen+2);
    hq2Int[iDet][19]->Draw("P");
    sq2Int[iDet][19]->SetLineColor(kGreen);
    sq2Int[iDet][19]->SetLineWidth(2);
    sq2Int[iDet][19]->Draw("same");
    hq2Int[iDet][10]->SetMarkerStyle(kFullSquare);
    hq2Int[iDet][10]->SetMarkerColor(kRed+2);
    hq2Int[iDet][10]->Draw("Psame");
    sq2Int[iDet][10]->SetLineColor(kRed);
    sq2Int[iDet][10]->SetLineWidth(2);
    sq2Int[iDet][10]->Draw("same");
    hq2Int[iDet][0]->SetMarkerStyle(kFullDiamond);
    hq2Int[iDet][0]->SetMarkerColor(kBlue);
    hq2Int[iDet][0]->Draw("Psame");
    sq2Int[iDet][0]->SetLineWidth(2);
    sq2Int[iDet][0]->SetLineColor(kCyan);
    sq2Int[iDet][0]->Draw("same");
    leg->Draw("same");
  }

  return 0;
}
