#include <stdio.h>
#include <TString.h>
#include <THnSparse.h>
#include <TFile.h>
#include <TH1D.h>
#include <TDirectoryFile.h>
#include <TList.h>
#include <TSpline.h>

//Demo macro to get percentiles from the Spherocity D mesons task output 
//Enter your multiplicity intervals in the main QuantileUsingSplines function and run it using ROOT
//The macro was tested using ROOT6 on a Ubuntu machine
//The obtained output needs to be fed to The AliAnalysisTaskSEDvsEventShapes in order to obtain spherocity percentiles
//For any questions text Marco Giacalone via Mattermost or email: mgiacalo@cern.ch

void Quantiles(int year = 19, TString fname = "AnalysisResult.root", TString fnamemc = "AnalysisResult.root", float minmul = 1., float maxmul = 200.) {

  //GenComp is for computed Spherocity, while Gen is for true spherocity

  TString baseSpline = Form("DSplines%d", year);
  if (year > 18 | year < 16)
    baseSpline = "DSplines";
  TString myfile = fnamemc.Data();//"MCGenPurpose/4380-2018/AnalysisResults.root"; //for now running only for MC
  TString myfiledata = fname.Data();//"Data/4794-2018/AnalysisResults.root";
  TFile *fil_1 = TFile::Open(myfile.Data(),"READ");
  TFile *fildata = TFile::Open(myfiledata.Data(),"READ");

  TDirectoryFile *ddata = (TDirectoryFile*)fildata->Get("PWG3_D2H_DEvtShape_DStarmgiacalo_sph_TPConly_S0unweighted");
  TList *listdata = (TList*)ddata->Get("coutputDStarmgiacalo_sph_TPConly_S0unweighted"); 

  TDirectoryFile *d_2 = (TDirectoryFile*)fil_1->Get("PWG3_D2H_DEvtShape_DStarmgiacaloMC_sph_TPConly_S0unweighted");
  TList *listMC = (TList*)d_2->Get("coutputDStarmgiacaloMC_sph_TPConly_S0unweighted");
  
  //MC Steps Splines

  THnSparse *thnMC = (THnSparse*)listMC->FindObject("hTotalSparseEvtShape");
  thnMC->GetAxis(1)->SetRangeUser(minmul,maxmul-1);
  TH1D *hMC=(TH1D*)thnMC->Projection(0);
  TH1D *hGen=(TH1D*)thnMC->Projection(2);
  const Int_t nBinsMC = hMC->GetXaxis()->GetNbins();
  const Int_t nBinsGen = hGen->GetXaxis()->GetNbins();

  //This is for Data Splines

  THnSparse *thnData = (THnSparse*)listdata->FindObject("hTotalSparseEvtShape");
  thnData->GetAxis(1)->SetRangeUser(minmul,maxmul-1);
  TH1D *hData=(TH1D*)thnData->Projection(0);
  const Int_t nBinsData = hData->GetXaxis()->GetNbins();

  const int numspline = 3;

  TH1F* splInt[numspline];
  
  TSpline3* spline[numspline];
  TList* splinelist = new TList();
  splinelist->SetOwner(0);

  //Int stands for interpolation

  for(Int_t k=0; k<numspline; k++){
    Double_t hINT=0;
    if(k==0){
      splInt[k] = new TH1F(Form("splIntMC%0.f%0.f",minmul,maxmul-1),Form("MCRecoSpherocity%0.f%0.f;Spherocity;Spherocity normalised integral",minmul,maxmul-1),nBinsMC,0.,1.);
      printf("Fitting Spherocity Splines MC Reco with Multiplicity [%0.f,%0.f]\n", minmul, maxmul-1);
      for(Int_t j=0; j<nBinsMC; j++){
        hINT+=hMC->GetBinContent(j+1);
        Double_t val = hMC->GetBinCenter(j+1);
        splInt[k]->SetBinContent(j+1,hINT/hMC->Integral()*100);
      }
      spline[k] = new TSpline3(splInt[k]);
      spline[k]->SetName(Form("%sMC",baseSpline.Data()));
      cout << "Spline at 0.6 for Reco is " << spline[k]->Eval(0.6) << endl;
      splinelist->Add(spline[k]);
    }
    else if(k==1){
      splInt[k] = new TH1F(Form("splIntGen%0.f%0.f",minmul,maxmul-1),Form("GenSpherocity%0.f%0.f;Spherocity;Spherocity normalised integral",minmul,maxmul-1),nBinsGen,0.,1.);
      printf("Fitting Spherocity Splines MC Gen with Multiplicity [%0.f,%0.f]\n", minmul, maxmul-1);
      for(Int_t j=0; j<nBinsGen; j++){
        hINT+=hGen->GetBinContent(j+1);
        Double_t val = hGen->GetBinCenter(j+1);
        splInt[k]->SetBinContent(j+1,hINT/hGen->Integral()*100);
      }
      spline[k] = new TSpline3(splInt[k]);
      spline[k]->SetName(Form("%sMCGen",baseSpline.Data()));
      cout << "Spline at 0.6 for Gen is " << spline[k]->Eval(0.6) << endl;
      splinelist->Add(spline[k]);
    }
    else if(k==2){
      splInt[k] = new TH1F(Form("splIntData%0.f%0.f",minmul,maxmul-1),Form("DataSpherocity%0.f%0.f;Spherocity;Spherocity normalised integral",minmul,maxmul-1),nBinsData,0.,1.);
      printf("Fitting Spherocity Splines Reco Data with Multiplicity [%0.f,%0.f]\n", minmul, maxmul-1);
      for(Int_t j=0; j<nBinsData; j++){
        hINT+=hData->GetBinContent(j+1);
        Double_t val = hData->GetBinCenter(j+1);
        splInt[k]->SetBinContent(j+1,hINT/hData->Integral()*100);
      }
      spline[k] = new TSpline3(splInt[k]);
      spline[k]->SetName(Form("%sData",baseSpline.Data()));
      cout << "Spline at 0.6 for Data is " << spline[k]->Eval(0.6) << endl;
      splinelist->Add(spline[k]);
    }
  }

  //Plot some examples
  TCanvas* cspheroInt[numspline];
  TCanvas* cspherohist[numspline];
  TLegend* leg = new TLegend(0.6,0.2,0.8,0.4);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.035);

  for(int k=0; k<numspline; k++){
    if(k==0)
      cspherohist[k] = new TCanvas(Form("cspheroHistMC_%0.f_%0.f",minmul,maxmul-1),"",800,800);
    else if(k==1)
      cspherohist[k] = new TCanvas(Form("cspheroHistGen_%0.f_%0.f",minmul,maxmul-1),"",800,800);
    else if(k==2)
      cspherohist[k] = new TCanvas(Form("cspheroHistData_%0.f_%0.f",minmul,maxmul-1),"",800,800);   
    leg->AddEntry(splInt[k],Form("Mult %0.f-%0.f%%",minmul,maxmul-1),"p"); 
    splInt[k]->Draw();
    TLine* line25 = new TLine(0.,25.,1., 25.);
    TLine* line50 = new TLine(0.,50.,1., 50.);
    TLine* line75 = new TLine(0.,75.,1., 75.);
    line25->SetLineColor(kRed);
    line50->SetLineColor(kRed);
    line75->SetLineColor(kRed);
    line25->Draw("same");
    line50->Draw("same");
    line75->Draw("same");
  }

  for(int k=0; k<numspline; k++){
    if(k==0)
      cspheroInt[k] = new TCanvas(Form("cspheroIntMC_%0.f_%0.f",minmul,maxmul-1),"",800,800);
    else if(k==1)
      cspheroInt[k] = new TCanvas(Form("cspheroIntGen_%0.f_%0.f",minmul,maxmul-1),"",800,800); 
    else if(k==2)
      cspheroInt[k] = new TCanvas(Form("cspheroIntData_%0.f_%0.f",minmul,maxmul-1),"",800,800); 
    leg->AddEntry(splInt[k],Form("Mult %0.f-%0.f%%",minmul,maxmul-1),"p"); 
    splInt[k]->SetMarkerStyle(kFullCircle);
    splInt[k]->SetMarkerColor(kGreen+2);
    splInt[k]->Draw("P");
    spline[k]->SetLineColor(kGreen);
    spline[k]->SetLineWidth(2);
    spline[k]->Draw("same");   
    cspheroInt[k]->Print(Form("%s.png", splInt[k]->GetTitle()));
  }

  //save splines in a file  
  TFile outfile("Splines.root","UPDATE");
  if(minmul == 1. && maxmul == 200.)
    splinelist->Write(Form("%s", baseSpline.Data()),1);
  else  
    splinelist->Write(Form("%s%0.f-%0.f", baseSpline.Data(),minmul, maxmul-1),1);
  /*cspheroInt[0]->Write();
  cspheroInt[1]->Write();
  cspherohist[0]->Write();
  cspherohist[1]->Write();*/
  outfile.Close();

/*
  TCanvas *c1 = new TCanvas("c1","demo quantiles",10,10,700,900);
  c1->Divide(1,2);
  TGraph *gr0 = new TGraph(nq,xq[0],yq[0]);
  c1->cd(1);
  gPad->SetGrid();
  gr0->SetMarkerStyle(20);
  gr0->Draw();

  // show the quantiles in the bottom pad
  c1->cd(2);
  gPad->SetGrid();
  TGraph *gr = new TGraph(nq,xq[1],yq[1]);
  gr->SetMarkerStyle(21);
  gr->Draw("alp");*/
}

void QuantileUsingSplines(TString splfile = "Splines.root"){
  TFile* splif = TFile::Open(splfile.Data(), "READ");
  if(splif){
    splif->Close();
    if(remove(splfile.Data()) != 0)
      cerr << "Failed to Delete file " << splfile.Data() << std::endl;
    else
      cout << "File " << splfile.Data() << " Deleted" << std::endl;   
  }
  TString fname[] = {"Data/4803-4804-2016/AnalysisResults.root", "Data/4805-2017/AnalysisResults_merged.root", "Data/4806-2018/AnalysisResults_merged.root", "Data/FullMerge.root"};
  TString fnamemc[] = {"MCGenPurpose/4390-2016/AnalysisResults.root", "MCGenPurpose/4392-2018/AnalysisResults.root", "MCGenPurpose/4393-2017/AnalysisResults.root", "MCGenPurpose/FullMerge.root"};
  const float mults[] = {20.,31.,82.};
  int years[] = {16, 17, 18, 19};
  for(int k = 3; k < 4; k++){
    for(int j = 0; j<=3; j++){
      if(j == 2)
        Quantiles(years[k], fname[k], fnamemc[k], mults[0],mults[j]);
      else if(j == 3)
        Quantiles(years[k], fname[k], fnamemc[k], 1., 200.);
      else
        Quantiles(years[k], fname[k], fnamemc[k], mults[j],mults[j+1]);
    }
  }
}

void comparison(){
  TFile * spl = new TFile("Splines.root", "READ");
  TCanvas * c[3];
  TLegend* leg[3];
  TString names[] = {"Data", "MC", "MCGen"};
  int years[] = {16, 17, 18, 19};
  int col[] = {kRed, kBlue, kGreen, kBlack};
  TList * li[4];
  TString baseSpline = "SpheroSpline";
  TSpline3 * splines[3];
  for(int z = 0; z < 3; z++){
    c[z] = new TCanvas(names[z].Data(),names[z].Data(), 1920, 1080);
    leg[z] = new TLegend(0.6,0.2,0.8,0.4);
    leg[z]->SetBorderSize(0);
    leg[z]->SetTextSize(0.035);
  }
  for(int i = 0; i < 4; i++){
    if(years[i] == 19)
      li[i] = (TList*)spl->Get(Form("%s", baseSpline.Data()));
    else 
      li[i] = (TList*)spl->Get(Form("%s%d", baseSpline.Data(), years[i]));
    for(int k = 0; k < 3; k++){
      c[k]->cd();
      if(years[i] == 19)
        splines[k] = (TSpline3*)li[i]->FindObject(Form("%s%s", baseSpline.Data(), names[k].Data()));
      else
        splines[k] = (TSpline3*)li[i]->FindObject(Form("%s%d%s", baseSpline.Data(), years[i], names[k].Data()));
      splines[k]->SetLineColor(col[i]);
      if(years[i] == 19)
        leg[k]->AddEntry(splines[k],Form("%s%sMerged", baseSpline.Data(), names[k].Data()),"l");  
      else
        leg[k]->AddEntry(splines[k],Form("%s%d%s", baseSpline.Data(), years[i], names[k].Data()),"l");
      if(years[i] == 16)
        splines[k]->Draw();
      else  
        splines[k]->Draw("SAME");  
      leg[k]->Draw();  
    }
  }
  c[0]->Print("DataSplines.png");
  c[1]->Print("MCSplines.png");
  c[2]->Print("MCGenSplines.png");
}
