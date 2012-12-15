#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TH2F.h>
#include <TFile.h>
#include <TGrid.h>
#include <TCanvas.h>
#include <TStopwatch.h>
#include <TStyle.h>
#include <TProfile.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TLatex.h>
#include <TGrid.h>
#include "AliRawReaderDate.h"
#include "AliRawReaderRoot.h"
#include "AliITSRawStreamSDD.h"
#include "AliITSRawStreamSDDCompressed.h"
#endif


void CheckDataSizeSDD(TString datafil="12000188359004.16.root",
		      Int_t firstEv=0, 
		      Int_t lastEv=123456){

  if(datafil.Contains("alien:")) TGrid::Connect("alien:");
  printf("FILE: %s\n",datafil.Data());

  gStyle->SetOptStat(0);
  gStyle->SetTitleFont(42,"XY");
  gStyle->SetLabelFont(42,"XYZ");

  Double_t maxOcc=2500.;
  TH2F* hCellsOnMod=new TH2F("hCellsOnMod","",260,239.5,499.5,500,0.,maxOcc);
  TH2F* hCellsOnDDL=new TH2F("hCellsOnDDL","",24,-0.5,23.5,500,0.,maxOcc*12);
  TH1F* hMaxOccMod=new TH1F("hMaxOccMod","",260,239.5,499.5);
  TH1F* hMaxOccDDL=new TH1F("hMaxOccDDL","",24,-0.5,23.5);

  AliITSDDLModuleMapSDD* dmap=new AliITSDDLModuleMapSDD();
  dmap->SetJun09Map();


  Int_t iev=firstEv;
  AliRawReader *rd; 
  if(datafil.Contains(".root")){
    rd=new AliRawReaderRoot(datafil.Data(),iev);
  }else{
    rd=new AliRawReaderDate(datafil.Data(),iev);
  }

  Bool_t writtenoutput=kFALSE;
  Int_t countMod[260],countDDL[24];
  Int_t maxOccMod[260],maxOccDDL[24];
  for(Int_t im=0; im<260; im++) maxOccMod[im]=0;
  for(Int_t id=0; id<24; id++) maxOccDDL[id]=0;

  do{

    printf("Event # %d\n",iev);
    for(Int_t im=0; im<260; im++) countMod[im]=0;
    for(Int_t id=0; id<24; id++) countDDL[id]=0;
    rd->Reset();
    UChar_t cdhAttr=AliITSRawStreamSDD::ReadBlockAttributes(rd);
    UInt_t amSamplFreq=AliITSRawStreamSDD::ReadAMSamplFreqFromCDH(cdhAttr);
    AliITSRawStream* s=AliITSRawStreamSDD::CreateRawStreamSDD(rd,cdhAttr);
    if(!writtenoutput){
      printf("Use %s raw stream, sampling frequency %d MHz\n",s->ClassName(),amSamplFreq);
      writtenoutput=kTRUE;
    }

    while(s->Next()){
      if(s->IsCompletedModule()==kFALSE && s->IsCompletedDDL()==kFALSE){
	Int_t counts=s->GetSignal();
	if(counts>0){
	  Int_t iDDL=rd->GetDDLID();
	  ++countDDL[iDDL];
	  Int_t iMod=s->GetCarlosId();
	  Int_t modInd=dmap->GetModuleNumber(iDDL,iMod)-240;
	  if(modInd>=0 && modInd<260) ++countMod[modInd];
	}
      }
    }
    for(Int_t im=0; im<260; im++){
      hCellsOnMod->Fill(im+240,countMod[im]);
      if(countMod[im]>maxOccMod[im]) maxOccMod[im]=countMod[im];
    }
    for(Int_t id=0; id<24; id++){
      hCellsOnDDL->Fill(id,countDDL[id]);
      if(countDDL[id]>maxOccDDL[id]) maxOccDDL[id]=countDDL[id];
    }
    iev++;
    
  }while(rd->NextEvent()&&iev<=lastEv);

  for(Int_t im=0; im<260; im++) hMaxOccMod->SetBinContent(im+1,maxOccMod[im]);
  for(Int_t id=0; id<24; id++) hMaxOccDDL->SetBinContent(id+1,maxOccDDL[id]);

  TProfile* hAveOccDDL=hCellsOnDDL->ProfileX();
  TProfile* hAveOccMod=hCellsOnMod->ProfileX();
  TH1D* hOccDstDDL[24];
  Int_t totEv=0;
  for(Int_t i=0; i<24;i++){
    hOccDstDDL[i]=(TH1D*)hCellsOnDDL->ProjectionY(Form("hOccDstDDL%d",i),i+1,i+1);
    printf("DDL %d  Entries %.0f\n",i,hOccDstDDL[i]->Integral());
    totEv=hOccDstDDL[i]->Integral();
  }
  Double_t xmax=hMaxOccDDL->GetMaximum()*1.1;
  Int_t ddlForPlot[4]={9,15,17,21};

  TCanvas* c1= new TCanvas("c1","Module occupancy",1100,750);
  c1->Divide(2,2);
  c1->cd(1);
  gPad->SetLogz();
  hCellsOnMod->GetXaxis()->SetTitle("Module Id");
  hCellsOnMod->GetYaxis()->SetTitle("Number of cells on per event");
  hCellsOnMod->GetYaxis()->SetTitleOffset(1.3);
  hCellsOnMod->Draw("colz");
  c1->cd(3);
  hAveOccMod->GetXaxis()->SetTitle("Module Id");
  hAveOccMod->GetYaxis()->SetTitle("<Occupancy>");
  hAveOccMod->GetYaxis()->SetTitleOffset(1.3);
  hAveOccMod->Draw();
  c1->cd(4);
  hMaxOccMod->GetXaxis()->SetTitle("Module Id");
  hMaxOccMod->GetYaxis()->SetTitle("Maximum Occupancy");
  hMaxOccMod->GetYaxis()->SetTitleOffset(1.3);
  hMaxOccMod->Draw();

  TCanvas* c2= new TCanvas("c2","DDL occupancy",1100,750);
  c2->Divide(2,2);
  c2->cd(1);
  gPad->SetLogz();
  hCellsOnDDL->GetXaxis()->SetTitle("DDL Number");
  hCellsOnDDL->GetYaxis()->SetTitle("Number of cells on per event");
  hCellsOnDDL->GetYaxis()->SetTitleOffset(1.3);
  hCellsOnDDL->Draw("colz");
  TLatex* textev=new TLatex(0.15,0.8,Form("%d events",totEv));
  textev->SetNDC();
  textev->Draw();
  c2->cd(2);
  gPad->SetLogy();
  TLegend* leg=new TLegend(0.6,0.6,0.89,0.89);
  leg->SetFillColor(0);
  for(Int_t i=0; i<4; i++){
    Int_t iddl=ddlForPlot[i];
    hOccDstDDL[iddl]->GetXaxis()->SetTitle("Occupancy");
    hOccDstDDL[iddl]->GetXaxis()->SetRangeUser(0.,xmax);
    hOccDstDDL[iddl]->SetLineColor(i+1);
    if(i==0) hOccDstDDL[iddl]->Draw();
    else hOccDstDDL[iddl]->Draw("same");
    leg->AddEntry(hOccDstDDL[iddl],Form("DDL %d",iddl),"L")->SetTextColor(hOccDstDDL[iddl]->GetLineColor());
  }
  leg->Draw();
  c2->cd(3);
  hAveOccDDL->GetXaxis()->SetTitle("DDL Number");
  hAveOccDDL->GetYaxis()->SetTitle("<Occupancy>");
  hAveOccDDL->GetYaxis()->SetTitleOffset(1.3);
  hAveOccDDL->Draw();
  c2->cd(4);
  hMaxOccDDL->GetXaxis()->SetTitle("DDL Number");
  hMaxOccDDL->GetYaxis()->SetTitle("Maximum Occupancy");
  hMaxOccDDL->GetYaxis()->SetTitleOffset(1.3);
  hMaxOccDDL->Draw();

  TString outfilname=datafil.Data();
  if(outfilname.Contains("alien")){
    outfilname.Remove(0,outfilname.Length()-22);
  }
  outfilname.Prepend("DataSize_");
  printf("%s\n",outfilname.Data());
  TFile* outf=new TFile(outfilname.Data(),"recreate");
  hCellsOnMod->Write();
  hMaxOccMod->Write();
  hCellsOnDDL->Write();
  hMaxOccDDL->Write();
  outf->Close();
}
