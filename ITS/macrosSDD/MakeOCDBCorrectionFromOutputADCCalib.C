#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TH1F.h>
#include <TH2F.h>
#include <TF1.h>
#include <TPad.h>
#include <TGraphErrors.h>
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TGrid.h>
#include <TGridResult.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLatex.h>
#include <AliCDBEntry.h>
#include <AliITSresponseSDD.h>
#include <AliITSCalibrationSDD.h>
#include <AliCDBManager.h>

#endif

// Inputs: Output of MakeSDDCalib.C + AliITSresponseSDD Object used to produce the cpass1 (default Run192772_192779_v0_s0.root)
// Outputs: new ADCtoKeV and ADCvsDriftTimeHistos
// Launch as: aliroot MakeCorrectionFromOutputCalib.C+

void MakeOCDBCorrectionFromOutputADCCalib(Int_t run = 246994, TString dirname = ".",TString fperiod = "",TString OCDBname = "Run192772_192779_v0_s0") {
    
  TGrid::Connect("alien:");
  AliCDBManager * man = AliCDBManager::Instance();
  man->SetDefaultStorage("alien://folder=/alice/data/2015/OCDB");
  man->SetRun(run);
  AliCDBEntry *pedpul = man->Get("ITS/Calib/CalibSDD");
  TObjArray *calSDD = (TObjArray *)pedpul->GetObject();
 
  FILE* outtxt=fopen(Form("%s/outlargecorr_%d.txt",dirname.Data(),run),"w");

  printf("Opening file CalibResults.root in %s/",dirname.Data());
  TFile *fin = new TFile(Form("%s/SDDADCCalibResults_%d.root",dirname.Data(),run));
  if(fin) {
    TH1F *hmpvModpar0 = (TH1F*)fin->Get("hmpvModpar0");
    TH1F *hmpvModpar1 = (TH1F*)fin->Get("hmpvModpar1");
        
    // getting ADC2keV
    TFile* fr = TFile::Open(Form("%s/%s.root",dirname.Data(),OCDBname.Data()));
    AliCDBEntry* e = (AliCDBEntry*)fr->Get("AliCDBEntry");
    AliITSresponseSDD* r = (AliITSresponseSDD*)e->GetObject();
        
    TH1F* hADCtokeV       = new TH1F("hADCtokeV","",260,239.5,499.5);
    TH1F* hADCvsDriftTime = new TH1F("hADCvsDriftTime","",260,239.5,499.5);
    r->ls();
        
    for(Int_t iMod = 240; iMod < 500; iMod++){
      Float_t ak = r->GetADCtokeV(iMod);
      Float_t adcdrtime = r->GetADCvsDriftTime(iMod);
      printf("mod %d: adc->keV = %f, ADCvsDriftTime = %f\n ",iMod,ak,adcdrtime);
      hADCtokeV->SetBinContent(iMod-240+1,ak);
      hADCvsDriftTime->SetBinContent(iMod-240+1,adcdrtime);
      //calculate new ADC2KeV
      Double_t kprim = hmpvModpar0->GetBinContent(hmpvModpar0->FindBin(iMod));
      if(kprim==0)  kprim = 84;
      if(iMod==376) kprim = 84;
      Double_t Corr = kprim*ak/84;
      printf("mpv=%f | Corr(adc->keV*mpv/84) = %f\n",kprim,Corr); //84 is the expected MPV for a Minimum Ionizing Particle crossing 300 um of Silicon
      hmpvModpar0->SetBinContent(hmpvModpar0->FindBin(iMod),Corr);
      //Calculate new ADCvsDrTime
      Double_t sprim = hmpvModpar1->GetBinContent(hmpvModpar1->FindBin(iMod));
      if(iMod==376) sprim = 0;
      Double_t Corr2 = adcdrtime - (sprim*ak/2);
      hmpvModpar1->SetBinContent(hmpvModpar1->FindBin(iMod),Corr2);
      printf("ADCvsDriftTime=%f | Corr(ADCvsDriftTime_old-(ADCvsDriftTime_new*adc->keV/2)) = %f\n",adcdrtime,Corr2);
    }
        
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    hmpvModpar0->GetYaxis()->SetTitle("K^{*}=KK'/84");
    hmpvModpar1->GetYaxis()->SetTitle("S^{*}=S-S'K/2");
    hmpvModpar0->SetTitle("hADCtokeV - New");
    hmpvModpar0->SetName("hADCtokeV");
    hADCtokeV->SetTitle(Form("hADCtokeV - %s",OCDBname.Data()));
    hADCtokeV->SetLineColor(2);
    TCanvas *c0=new TCanvas("c0","c0",1400,900);
    c0->Divide(2,2);
    c0->cd(1);
    hmpvModpar0->DrawClone("HIST");
    hADCtokeV->DrawClone("HISTSAME");
    TLegend* leg0=new TLegend(0.14,0.75,0.89,0.89);
    leg0->SetMargin(0.15);
    leg0->AddEntry(hmpvModpar0,hmpvModpar0->GetTitle(),"L")->SetTextColor(hmpvModpar0->GetLineColor());
    leg0->AddEntry(hADCtokeV,hADCtokeV->GetTitle(),"L")->SetTextColor(hADCtokeV->GetLineColor());
    leg0->Draw();
    c0->cd(3);
    TH1F* hRatioADC=(TH1F*)hmpvModpar0->Clone("hRatioADC");
    hRatioADC->Divide(hmpvModpar0,hADCtokeV);
    for(Int_t j=1; j<=hRatioADC->GetNbinsX(); j++) hRatioADC->SetBinError(j,0.0001);
    hRatioADC->SetStats(0);
    hRatioADC->SetLineWidth(2);
    hRatioADC->SetMarkerStyle(7);
    hRatioADC->GetYaxis()->SetTitle("Ratio new/OCDB");
    hRatioADC->GetYaxis()->SetTitle("Ratio new/OCDB");
    hRatioADC->Draw();
    c0->cd(2);
    hmpvModpar1->SetTitle("hADCvsDriftTime - New");
    hmpvModpar1->SetName("hADCvsDriftTime");
    hADCvsDriftTime->SetTitle(Form("hADCvsDriftTime %s",OCDBname.Data()));//- Run166530_999999999_v4_s0");
    hADCvsDriftTime->SetLineColor(2);
    hmpvModpar1->DrawClone("HIST");
    hADCvsDriftTime->DrawClone("HISTSAME");
    TLegend* leg1=new TLegend(0.14,0.75,0.89,0.89);
    leg1->SetMargin(0.15);
    leg1->AddEntry(hmpvModpar1,hmpvModpar1->GetTitle(),"L")->SetTextColor(hmpvModpar1->GetLineColor());
    leg1->AddEntry(hADCvsDriftTime,hADCvsDriftTime->GetTitle(),"L")->SetTextColor(hADCvsDriftTime->GetLineColor());
    leg1->Draw();
    c0->cd(4);
    TH1F* hRatioAvsDT=(TH1F*)hmpvModpar1->Clone("hRatioAvsDT");
    hRatioAvsDT->Divide(hmpvModpar1,hADCvsDriftTime);
    for(Int_t j=1; j<=hRatioADC->GetNbinsX(); j++) hRatioAvsDT->SetBinError(j,0.0001);
    hRatioAvsDT->SetStats(0);
    hRatioAvsDT->SetMarkerStyle(7);
    hRatioAvsDT->SetLineWidth(2);
    hRatioAvsDT->GetYaxis()->SetTitle("Ratio new/OCDB");
    hRatioAvsDT->Draw();

    for(Int_t j=1; j<=hRatioADC->GetNbinsX(); j++){
      Double_t r = hRatioADC->GetBinContent(j);
      Int_t theMod=(Int_t)(hRatioADC->GetBinCenter(j)+0.1);
      if(r<0.96 || r>1.04){ 
	printf("Large Correction (%f) for module %d",r,theMod);
	AliITSCalibrationSDD *cal=(AliITSCalibrationSDD*)calSDD->At(theMod-240);
	if(cal->IsBad()) printf("  --> BAD module\n");
	else{ 
	  printf ("  **** Good module ****\n");
	  fprintf(outtxt,"Large Correction (%f) for module %d **** Good module with %d dead channels****\n",r,theMod,cal->GetDeadChannels());
	}
      }
    }

    TFile *out=new TFile(Form("%s/CorrectiondEdxSDD_%s_%s_%d.root",dirname.Data(),fperiod.Data(),OCDBname.Data(),run),"recreate");
    hmpvModpar0->Write();
    hmpvModpar1->Write();
    c0->SaveAs(Form("%s/CalibParams_run%d.eps",dirname.Data(),run));
    out->Close();
    delete out;


  } else printf("no input file!\n");
} //end main

