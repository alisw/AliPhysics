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

void Quantiles(float minmul = 20., float maxmul = 82.) {

  //GenComp is for computed Spherocity, while Gen is for true spherocity

  TString baseSpline = "DStarSpline";
  TString myfile ="../MC/FullMerge.root";
  TString myfiledata ="../ExpData/FullMerge.root";
  TFile *fil_1 = TFile::Open(myfile.Data(),"READ");
  TFile *fildata = TFile::Open(myfiledata.Data(),"READ");

  TDirectoryFile *ddata = (TDirectoryFile*)fildata->Get("PWG3_D2H_DEvtShape_DStarmgiacalo_sph_TPConly_S0unweighted");
  TList *listdata = (TList*)ddata->Get("coutputDStarmgiacalo_sph_TPConly_S0unweighted"); 

  TDirectoryFile *d_2 = (TDirectoryFile*)fil_1->Get("PWG3_D2H_DEvtShape_DStarmgiacaloMC_sph_TPConly_S0unweighted");
  TList *list2 = (TList*)d_2->Get("coutputEffCorrDStarmgiacaloMC_sph_TPConly_S0unweighted"); 
  TList *listMC = (TList*)d_2->Get("coutputDStarmgiacaloMC_sph_TPConly_S0unweighted");
  
  //MC Steps Splines

  THnSparse *thnReco = (THnSparse*)list2->FindObject("hMCRecoPrompt");
  THnSparse *hMCAccGenPrompt =(THnSparse*)list2->FindObject("hMCAccGenPrompt");
  THnSparse *thnRecoFD = (THnSparse*)list2->FindObject("hMCRecoFeeddown");
  THnSparse *hMCAccGenFD =(THnSparse*)list2->FindObject("hMCAccGenBFeeddown");
  thnRecoFD->GetAxis(1)->SetRangeUser(minmul,maxmul-1);
  hMCAccGenFD->GetAxis(1)->SetRangeUser(minmul,maxmul-1);
  thnReco->GetAxis(1)->SetRangeUser(minmul,maxmul-1);
  hMCAccGenPrompt->GetAxis(1)->SetRangeUser(minmul,maxmul-1);
  TH1D *hReco=(TH1D*)thnReco->Projection(2);
  TH1D *hGen=(TH1D*)hMCAccGenPrompt->Projection(5);
  TH1D *hGenComp=(TH1D*)hMCAccGenPrompt->Projection(2);
  TH1D *hRecoFD=(TH1D*)thnRecoFD->Projection(2);
  TH1D *hGenFD=(TH1D*)hMCAccGenFD->Projection(2);
  const Int_t nBinsReco = hReco->GetXaxis()->GetNbins();
  const Int_t nBinsGen = hGen->GetXaxis()->GetNbins();
  const Int_t nBinsGenComp = hGenComp->GetXaxis()->GetNbins();
  const Int_t nBinsRecoFD = hRecoFD->GetXaxis()->GetNbins();
  const Int_t nBinsGenFD = hGenFD->GetXaxis()->GetNbins();

  //This is for Data Splines

  THnSparse *thnData = (THnSparse*)listdata->FindObject("hSparseEvtShape");
  thnData->GetAxis(3)->SetRangeUser(minmul,maxmul-1);
  TH1D *hData=(TH1D*)thnData->Projection(2);
  const Int_t nBinsData = hData->GetXaxis()->GetNbins();

  //This is for MC Splines
  
  THnSparse *thnMC = (THnSparse*)listMC->FindObject("hSparseEvtShape");
  thnMC->SetName("hSparseEvtShapeMC");
  thnMC->GetAxis(3)->SetRangeUser(minmul,maxmul-1);
  TH1D *hMC=(TH1D*)thnMC->Projection(2);
  const Int_t nBinsMC = hMC->GetXaxis()->GetNbins();

  THnSparse *thnMCPrompt = (THnSparse*)listMC->FindObject("hSparseEvtShapePrompt");
  thnMCPrompt->SetName("hSparseEvtShapeMCPrompt");
  thnMCPrompt->GetAxis(3)->SetRangeUser(minmul,maxmul-1);
  TH1D *hMCPrompt=(TH1D*)thnMCPrompt->Projection(2);
  const Int_t nBinsMCPrompt = hMCPrompt->GetXaxis()->GetNbins();

  THnSparse *thnMCFD = (THnSparse*)listMC->FindObject("hSparseEvtShapeFeeddown");
  thnMCFD->SetName("hSparseEvtShapeMCFeeddown");
  thnMCFD->GetAxis(3)->SetRangeUser(minmul,maxmul-1);
  TH1D *hMCFD=(TH1D*)thnMCFD->Projection(2);
  const Int_t nBinsMCFD = hMCFD->GetXaxis()->GetNbins();

  const int numspline = 9;

  TH1F* splInt[numspline];
  
  TSpline3* spline[numspline];
  TList* splinelist = new TList();
  splinelist->SetOwner(0);

  for(Int_t k=0; k<numspline; k++){
    Double_t hINT=0;
    if(k==0){
      splInt[k] = new TH1F(Form("splIntReco%0.f%0.f",minmul,maxmul-1),Form("Reco Spherocity;Spherocity;Spherocity normalised integral"),nBinsReco,0.,1.);
      printf("Fitting Spherocity Splines Reco with Multiplicity [%0.f,%0.f]\n", minmul, maxmul-1);
      for(Int_t j=0; j<nBinsReco; j++){
        hINT+=hReco->GetBinContent(j+1);
        Double_t val = hReco->GetBinCenter(j+1);
        splInt[k]->SetBinContent(j+1,hINT/hReco->Integral()*100);
      }
      spline[k] = new TSpline3(splInt[k]);
      spline[k]->SetName(Form("%sReco_%0.f_%0.f",baseSpline.Data(),minmul,maxmul-1));
      cout << "Spline at 0.6 for Reco is " << spline[k]->Eval(0.6) << endl;
      splinelist->Add(spline[k]);
    }
    else if(k==1){
      splInt[k] = new TH1F(Form("splIntGen%0.f%0.f",minmul,maxmul-1),Form("Gen Spherocity;Spherocity;Spherocity normalised integral"),nBinsGen,0.,1.);
      printf("Fitting Spherocity Splines Gen with Multiplicity [%0.f,%0.f]\n", minmul, maxmul-1);
      for(Int_t j=0; j<nBinsGen; j++){
        hINT+=hGen->GetBinContent(j+1);
        Double_t val = hGen->GetBinCenter(j+1);
        splInt[k]->SetBinContent(j+1,hINT/hGen->Integral()*100);
      }
      spline[k] = new TSpline3(splInt[k]);
      spline[k]->SetName(Form("%sGen_%0.f_%0.f",baseSpline.Data(),minmul,maxmul-1));
      cout << "Spline at 0.6 for Gen is " << spline[k]->Eval(0.6) << endl;
      splinelist->Add(spline[k]);
    }
    else if(k==2){
      splInt[k] = new TH1F(Form("splIntGenComp%0.f%0.f",minmul,maxmul-1),Form("GenComp Spherocity;Spherocity;Spherocity normalised integral"),nBinsGenComp,0.,1.);
      printf("Fitting Spherocity Splines GenComp with Multiplicity [%0.f,%0.f]\n", minmul, maxmul-1);
      for(Int_t j=0; j<nBinsGenComp; j++){
        hINT+=hGenComp->GetBinContent(j+1);
        Double_t val = hGenComp->GetBinCenter(j+1);
        splInt[k]->SetBinContent(j+1,hINT/hGenComp->Integral()*100);
      }
      spline[k] = new TSpline3(splInt[k]);
      spline[k]->SetName(Form("%sGenComp_%0.f_%0.f",baseSpline.Data(),minmul,maxmul-1));
      cout << "Spline at 0.6 for GenComp is " << spline[k]->Eval(0.6) << endl;
      splinelist->Add(spline[k]);
    }
    else if(k==3){
      splInt[k] = new TH1F(Form("splIntRecoFD%0.f%0.f",minmul,maxmul-1),Form("RecoFD Spherocity;Spherocity;Spherocity normalised integral"),nBinsRecoFD,0.,1.);
      printf("Fitting Spherocity Splines RecoFD with Multiplicity [%0.f,%0.f]\n", minmul, maxmul-1);
      for(Int_t j=0; j<nBinsRecoFD; j++){
        hINT+=hRecoFD->GetBinContent(j+1);
        Double_t val = hRecoFD->GetBinCenter(j+1);
        splInt[k]->SetBinContent(j+1,hINT/hRecoFD->Integral()*100);
      }
      spline[k] = new TSpline3(splInt[k]);
      spline[k]->SetName(Form("%sRecoFD_%0.f_%0.f",baseSpline.Data(),minmul,maxmul-1));
      cout << "Spline at 0.6 for RecoFD is " << spline[k]->Eval(0.6) << endl;
      splinelist->Add(spline[k]);
    }
    else if(k==4){
      splInt[k] = new TH1F(Form("splIntGenFD%0.f%0.f",minmul,maxmul-1),Form("GenFD Spherocity;Spherocity;Spherocity normalised integral"),nBinsGenFD,0.,1.);
      printf("Fitting Spherocity Splines GenFD with Multiplicity [%0.f,%0.f]\n", minmul, maxmul-1);
      for(Int_t j=0; j<nBinsGenFD; j++){
        hINT+=hGenFD->GetBinContent(j+1);
        Double_t val = hGenFD->GetBinCenter(j+1);
        splInt[k]->SetBinContent(j+1,hINT/hGenFD->Integral()*100);
      }
      spline[k] = new TSpline3(splInt[k]);
      spline[k]->SetName(Form("%sGenFD_%0.f_%0.f",baseSpline.Data(),minmul,maxmul-1));
      cout << "Spline at 0.6 for GenFD is " << spline[k]->Eval(0.6) << endl;
      splinelist->Add(spline[k]);
    }
    else if(k==5){
      splInt[k] = new TH1F(Form("splIntData%0.f%0.f",minmul,maxmul-1),Form("Data Spherocity;Spherocity;Spherocity normalised integral"),nBinsData,0.,1.);
      printf("Fitting Spherocity Splines Data with Multiplicity [%0.f,%0.f]\n", minmul, maxmul-1);
      for(Int_t j=0; j<nBinsData; j++){
        hINT+=hData->GetBinContent(j+1);
        Double_t val = hData->GetBinCenter(j+1);
        splInt[k]->SetBinContent(j+1,hINT/hData->Integral()*100);
      }
      spline[k] = new TSpline3(splInt[k]);
      spline[k]->SetName(Form("%sData_%0.f_%0.f",baseSpline.Data(),minmul,maxmul-1));
      cout << "Spline at 0.6 for Data is " << spline[k]->Eval(0.6) << endl;
      splinelist->Add(spline[k]);
    }
    else if(k==6){
      splInt[k] = new TH1F(Form("splIntMC%0.f%0.f",minmul,maxmul-1),Form("MC Spherocity;Spherocity;Spherocity normalised integral"),nBinsMC,0.,1.);
      printf("Fitting Spherocity Splines MC with Multiplicity [%0.f,%0.f]\n", minmul, maxmul-1);
      for(Int_t j=0; j<nBinsMC; j++){
        hINT+=hMC->GetBinContent(j+1);
        Double_t val = hMC->GetBinCenter(j+1);
        splInt[k]->SetBinContent(j+1,hINT/hMC->Integral()*100);
      }
      spline[k] = new TSpline3(splInt[k]);
      spline[k]->SetName(Form("%sMC_%0.f_%0.f",baseSpline.Data(),minmul,maxmul-1));
      cout << "Spline at 0.6 for MC is " << spline[k]->Eval(0.6) << endl;
      splinelist->Add(spline[k]);
    }
    else if(k==7){
      splInt[k] = new TH1F(Form("splIntMCPrompt%0.f%0.f",minmul,maxmul-1),Form("MCPrompt Spherocity;Spherocity;Spherocity normalised integral"),nBinsMCPrompt,0.,1.);
      printf("Fitting Spherocity Splines MCPrompt with Multiplicity [%0.f,%0.f]\n", minmul, maxmul-1);
      for(Int_t j=0; j<nBinsMC; j++){
        hINT+=hMCPrompt->GetBinContent(j+1);
        Double_t val = hMCPrompt->GetBinCenter(j+1);
        splInt[k]->SetBinContent(j+1,hINT/hMCPrompt->Integral()*100);
      }
      spline[k] = new TSpline3(splInt[k]);
      spline[k]->SetName(Form("%sMCPrompt_%0.f_%0.f",baseSpline.Data(),minmul,maxmul-1));
      cout << "Spline at 0.6 for MCPrompt is " << spline[k]->Eval(0.6) << endl;
      splinelist->Add(spline[k]);
    }
    else if(k==8){
      splInt[k] = new TH1F(Form("splIntMCFD%0.f%0.f",minmul,maxmul-1),Form("MCFD Spherocity;Spherocity;Spherocity normalised integral"),nBinsMCFD,0.,1.);
      printf("Fitting Spherocity Splines MCFD with Multiplicity [%0.f,%0.f]\n", minmul, maxmul-1);
      for(Int_t j=0; j<nBinsMC; j++){
        hINT+=hMCFD->GetBinContent(j+1);
        Double_t val = hMCFD->GetBinCenter(j+1);
        splInt[k]->SetBinContent(j+1,hINT/hMCFD->Integral()*100);
      }
      spline[k] = new TSpline3(splInt[k]);
      spline[k]->SetName(Form("%sMCFD_%0.f_%0.f",baseSpline.Data(),minmul,maxmul-1));
      cout << "Spline at 0.6 for MCFD is " << spline[k]->Eval(0.6) << endl;
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
      cspherohist[k] = new TCanvas(Form("cspheroHistReco_%0.f_%0.f",minmul,maxmul-1),"",800,800);
    else if(k==1)
      cspherohist[k] = new TCanvas(Form("cspheroHistGen_%0.f_%0.f",minmul,maxmul-1),"",800,800);
    else if(k==2)
      cspherohist[k] = new TCanvas(Form("cspheroHistGenComp_%0.f_%0.f",minmul,maxmul-1),"",800,800);   
    else if(k==3)
      cspherohist[k] = new TCanvas(Form("cspheroHistRecoFD_%0.f_%0.f",minmul,maxmul-1),"",800,800); 
    else if(k==4)
      cspherohist[k] = new TCanvas(Form("cspheroHistGenFD_%0.f_%0.f",minmul,maxmul-1),"",800,800);     
    else if(k==5)
      cspherohist[k] = new TCanvas(Form("cspheroHistData_%0.f_%0.f",minmul,maxmul-1),"",800,800);     
    else if(k==6)
      cspherohist[k] = new TCanvas(Form("cspheroHistMC_%0.f_%0.f",minmul,maxmul-1),"",800,800); 
    else if(k==7)
      cspherohist[k] = new TCanvas(Form("cspheroHistMCPrompt_%0.f_%0.f",minmul,maxmul-1),"",800,800); 
    else if(k==8)
      cspherohist[k] = new TCanvas(Form("cspheroHistMCFD_%0.f_%0.f",minmul,maxmul-1),"",800,800);     
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
      cspheroInt[k] = new TCanvas(Form("cspheroIntReco_%0.f_%0.f",minmul,maxmul-1),"",800,800);
    else if(k==1)
      cspheroInt[k] = new TCanvas(Form("cspheroIntGen_%0.f_%0.f",minmul,maxmul-1),"",800,800); 
    else if(k==2)
      cspheroInt[k] = new TCanvas(Form("cspheroIntGenComp_%0.f_%0.f",minmul,maxmul-1),"",800,800); 
    else if(k==3)
      cspheroInt[k] = new TCanvas(Form("cspheroIntRecoFD_%0.f_%0.f",minmul,maxmul-1),"",800,800);
    else if(k==4)
      cspheroInt[k] = new TCanvas(Form("cspheroIntGenFD_%0.f_%0.f",minmul,maxmul-1),"",800,800);   
    else if(k==5)
      cspheroInt[k] = new TCanvas(Form("cspheroIntData_%0.f_%0.f",minmul,maxmul-1),"",800,800);  
    else if(k==6)
      cspheroInt[k] = new TCanvas(Form("cspheroIntMC_%0.f_%0.f",minmul,maxmul-1),"",800,800);
    else if(k==7)
      cspheroInt[k] = new TCanvas(Form("cspheroIntMCPrompt_%0.f_%0.f",minmul,maxmul-1),"",800,800);
    else if(k==8)
      cspheroInt[k] = new TCanvas(Form("cspheroIntMCFD_%0.f_%0.f",minmul,maxmul-1),"",800,800);          
    leg->AddEntry(splInt[k],Form("Mult %0.f-%0.f%%",minmul,maxmul-1),"p"); 
    splInt[k]->SetMarkerStyle(kFullCircle);
    splInt[k]->SetMarkerColor(kGreen+2);
    splInt[k]->Draw("P");
    spline[k]->SetLineColor(kGreen);
    spline[k]->SetLineWidth(2);
    spline[k]->Draw("same");   
  }

  //save splines in a file  
  TFile outfile("Splines.root","UPDATE");
  splinelist->Write(Form("%s%0.f-%0.f", baseSpline.Data(), minmul, maxmul-1),1);
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
  const int nmult = 2;
  const float mults[nmult+1]{20.,31.,82.};
  for(int k=0; k<=nmult; k++){
    if(k==nmult)
      Quantiles(mults[0],mults[k]);
    else
      Quantiles(mults[k],mults[k+1]); 
  }
}
