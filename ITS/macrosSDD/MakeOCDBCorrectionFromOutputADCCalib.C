#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TH1F.h>
#include <TH2F.h>
#include <TF1.h>
#include <TPad.h>
#include <TGraphErrors.h>
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TGrid.h>
#include <TGridResult.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLatex.h>
#include <AliCDBEntry.h>
#include <AliITSresponseSDD.h>

#endif

// Inputs: Output of MakeSDDCalib.C + AliITSresponseSDD Object used to produce the cpass1 (default Run192772_192779_v0_s0.root)
// Outputs: new ADCtoKeV and ADCvsDriftTimeHistos
// Launch as: aliroot MakeCorrectionFromOutputCalib.C+

void MakeOCDBCorrectionFromOutputADCCalib(TString fname=".",TString OCDBname="Run192772_192779_v0_s0")
{
  Printf("Opening file CalibResults.root in %s/",fname.Data());
  TFile *fin = new TFile(Form("%s/SDDADCCalibResults.root",fname.Data()));
    
  TH1F *hmpvModpar0=(TH1F*)fin->Get("hmpvModpar0");
  TH1F *hmpvModpar1=(TH1F*)fin->Get("hmpvModpar1");
  Float_t s=0.0101;//Default ADC vs drift time
  // getting ADC2keV
  TFile* fr=TFile::Open(Form("%s.root",OCDBname.Data()));
  AliCDBEntry* e=(AliCDBEntry*)fr->Get("AliCDBEntry");
  AliITSresponseSDD* r=(AliITSresponseSDD*)e->GetObject();
  TH1F* hADCtokeV=new TH1F("hADCtokeV","",260,239.5,499.5);
  TH1F* hADCvsDriftTime=new TH1F("hADCvsDriftTime","",260,239.5,499.5);
  r->ls();
    
  for(Int_t iMod=240; iMod<500; iMod++){
    Float_t ak=r->GetADCtokeV(iMod);
    Float_t adcdrtime=r->GetADCvsDriftTime(iMod);
    printf("mod %d\nadc->keV=%f ",iMod,ak);
    hADCtokeV->SetBinContent(iMod-240+1,ak);
    hADCvsDriftTime->SetBinContent(iMod-240+1,adcdrtime);
    //calculate new ADC2KeV
    Double_t kprim=hmpvModpar0->GetBinContent(hmpvModpar0->FindBin(iMod));
    if(kprim==0)kprim=84;
    if(iMod==376)kprim=84;
    Double_t Corr=kprim*ak/84;
    printf("mpv=%f | Corr(adc->keV*mpv/84)=%f\n",kprim,Corr);
    hmpvModpar0->SetBinContent(hmpvModpar0->FindBin(iMod),Corr);
    //Calculate new ADCvsDrTime
    Double_t sprim=hmpvModpar1->GetBinContent(hmpvModpar1->FindBin(iMod));
    if(iMod==376){
      sprim=0;
    }
    Double_t Corr2=adcdrtime-(sprim*ak/2);
    hmpvModpar1->SetBinContent(hmpvModpar1->FindBin(iMod),Corr2);
    printf("ADCvsDriftTime=%f | Corr(ADCvsDriftTime_old-(ADCvsDriftTime_new*adc->keV/2))=%f\n",adcdrtime,Corr2);
  }
    
  hmpvModpar0->GetYaxis()->SetTitle("K^{*}=KK'/84");
  hmpvModpar1->GetYaxis()->SetTitle("S^{*}=S-S'K/2");
  hmpvModpar0->SetTitle("hADCtokeV - New");
  hmpvModpar0->SetName("hADCtokeV");
  hADCtokeV->SetTitle("hADCtokeV - Run166530_999999999_v4_s0");
  hADCtokeV->SetLineColor(2);
  TCanvas *c0=new TCanvas("c0","c0");
  c0->cd(1);
  hmpvModpar0->DrawClone("HIST");
  hADCtokeV->DrawClone("HISTSAME");
  gPad->BuildLegend();
    
    
  hmpvModpar1->SetTitle("hADCvsDriftTime - New");
  hmpvModpar1->SetName("hADCvsDriftTime");
  hADCvsDriftTime->SetTitle("hADCvsDriftTime - Run166530_999999999_v4_s0");
  hADCvsDriftTime->SetLineColor(2);
  TCanvas *c1=new TCanvas("c1","c1");
  hmpvModpar1->DrawClone("HIST");
  hADCvsDriftTime->DrawClone("HISTSAME");
  c1->BuildLegend();
  TFile *out=new TFile(Form("%s/CorrectiondEdxSDD_%s_%s.root",fname.Data(),fname.Data(),OCDBname.Data()),"recreate");
  hmpvModpar0->Write();
  hmpvModpar1->Write();
  out->Close();
  delete out;
    
} //end main

