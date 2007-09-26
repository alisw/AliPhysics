/**************************************************************************
 * Copyright(c) 2007-2009, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/
#include <TFile.h>
#include "AliITSOnlineSDDCMN.h"
#include "AliLog.h"
#include <TH2F.h>
#include <TMath.h>


///////////////////////////////////////////////////////////////////
//                                                               //
// Implementation of the class used for analysis of SDD noise    //
// corrected for common mode                                     //
// Origin: F.Prino, Torino, prino@to.infn.it                     //
//                                                               //
///////////////////////////////////////////////////////////////////


ClassImp(AliITSOnlineSDDCMN)
//______________________________________________________________________
  AliITSOnlineSDDCMN::AliITSOnlineSDDCMN():AliITSOnlineSDD(),fNEvents(0),fMinCorrNoise(0.),fMaxCorrNoise(0.),fNSigmaNoise(0.)
{
  // default constructor
  Reset();
  SetMinNoise();
  SetMaxNoise();
  SetNSigmaNoise();
}
//______________________________________________________________________
  AliITSOnlineSDDCMN::AliITSOnlineSDDCMN(Int_t mod, Int_t sid):AliITSOnlineSDD(mod,sid),fNEvents(0),fMinCorrNoise(0.),fMaxCorrNoise(0.),fNSigmaNoise(0.)
{
  // default constructor
  Reset();
  SetMinNoise();
  SetMaxNoise();
  SetNSigmaNoise();
}
//______________________________________________________________________
AliITSOnlineSDDCMN::~AliITSOnlineSDDCMN(){
  // Destructor
}
//______________________________________________________________________
void AliITSOnlineSDDCMN::Reset(){
  //
  fNEvents=0;
  for(Int_t i=0;i<fgkNAnodes;i++){
    fGoodAnode[i]=1;
    fBaseline[i]=0.;
    fRawNoise[i]=0.;
    fCMN[i]=0.;
    fSumCorrNoise[i]=0.;
  }
  ReadBaselines();
}
//______________________________________________________________________
void AliITSOnlineSDDCMN::ReadBaselines(){
  // assume baselines and good anodes are taken from previous run
  Char_t basfilnam[100];
  sprintf(basfilnam,"SDDbase_step1_mod%03d_sid%d.data",fModuleId,fSide);
  FILE* basf=fopen(basfilnam,"r");
  if(basf==0){
    AliWarning("Baselinefile not present, Set all baselines to 50\n");
    for(Int_t ian=0;ian<fgkNAnodes;ian++){ 
      fBaseline[ian]=50.;
      fGoodAnode[ian]=1;
    }
    return;
  }
  Int_t n,ok;
  Float_t base,rms,cmn,corrnoi;
  for(Int_t ian=0;ian<fgkNAnodes;ian++){
    fscanf(basf,"%d %d %f %f %f %f\n",&n,&ok,&base,&rms,&cmn,&corrnoi);
    fGoodAnode[ian]=ok;
    fBaseline[ian]=base;
    fRawNoise[ian]=rms;
    fCMN[ian]=cmn;
  }
  fclose(basf);
}
//______________________________________________________________________
void  AliITSOnlineSDDCMN::ValidateAnodes(){
  //
  for(Int_t ian=0;ian<fgkNAnodes;ian++){
    if(!fGoodAnode[ian]) continue;
    if(GetAnodeCorrNoise(ian)>fMaxCorrNoise || GetAnodeCorrNoise(ian)<fMinCorrNoise) fGoodAnode[ian]=0;
    if(GetAnodeCorrNoise(ian)>fNSigmaNoise*CalcMeanNoise()) fGoodAnode[ian]=0;
  }
}

//______________________________________________________________________
void AliITSOnlineSDDCMN::AddEvent(TH2F* hrawd){
  // 
  fNEvents++;
  const Int_t kTimeBins=fLastGoodTB-fFirstGoodTB+1;
  TH2F* hcorrd=new TH2F("hcorrd","",hrawd->GetNbinsX(),hrawd->GetXaxis()->GetXmin(),hrawd->GetXaxis()->GetXmax(),hrawd->GetNbinsY(),hrawd->GetYaxis()->GetXmin(),hrawd->GetYaxis()->GetXmax());
  for(Int_t itb=fFirstGoodTB;itb<=fLastGoodTB;itb++){
    Float_t sumEven=0., sumOdd=0.;
    Int_t countEven=0, countOdd=0;
    for(Int_t ian=0;ian<fgkNAnodes;ian+=2){
      if(!fGoodAnode[ian]) continue;
      sumEven+=hrawd->GetBinContent(itb+1,ian+1)-fBaseline[ian];
      countEven++;
    }
    for(Int_t ian=1;ian<fgkNAnodes;ian+=2){
      if(!fGoodAnode[ian]) continue;
      sumOdd+=hrawd->GetBinContent(itb+1,ian+1)-fBaseline[ian];
      countOdd++;
    }
    for(Int_t ian=0;ian<fgkNAnodes;ian++){
      if(!fGoodAnode[ian]) continue;
      Float_t meanN;
      if(ian%2==0) meanN=sumEven/(Float_t)countEven;
      else meanN=sumOdd/(Float_t)countOdd;
      Float_t cntCorr=hrawd->GetBinContent(itb+1,ian+1)-fCMN[ian]*meanN;
      hcorrd->SetBinContent(itb+1,ian+1,cntCorr);
    }
  }

  for(Int_t ian=0;ian<fgkNAnodes;ian++){
    if(!fGoodAnode[ian]) continue;
    Float_t sumQ=0.;
    for(Int_t itb=fFirstGoodTB;itb<=fLastGoodTB;itb++){
      sumQ+=TMath::Power(hcorrd->GetBinContent(itb+1,ian+1)-fBaseline[ian],2);      
    }
    fSumCorrNoise[ian]+=TMath::Sqrt(sumQ/(Float_t)kTimeBins);
  }
  delete hcorrd;
}
//______________________________________________________________________
Float_t AliITSOnlineSDDCMN::CalcMeanNoise() const{
  //
  Float_t meanns=0.;
  Int_t cnt=0;
  for(Int_t ian=0;ian<fgkNAnodes;ian++){
    if(!fGoodAnode[ian]) continue;  
    meanns+=GetAnodeCorrNoise(ian);
    cnt++;
  }
  if(cnt>0) meanns/=(Float_t)cnt;
  return meanns;
}
//______________________________________________________________________
void AliITSOnlineSDDCMN::WriteToASCII(){
  //
  Char_t outfilnam[100];
  sprintf(outfilnam,"SDDbase_step2_mod%03d_sid%d.data",fModuleId,fSide);
  FILE* outf=fopen(outfilnam,"w");
  for(Int_t ian=0;ian<fgkNAnodes;ian++){
    fprintf(outf,"%d %d %8.3f %8.3f %8.3f %8.3f\n",ian,IsAnodeGood(ian),GetAnodeBaseline(ian),GetAnodeRawNoise(ian),GetAnodeCommonMode(ian),GetAnodeCorrNoise(ian));
  }
  fclose(outf);  
}

//______________________________________________________________________
TH1F* AliITSOnlineSDDCMN::GetBaselineAnodeHisto() const {
  //
  Char_t hisnam[20];  
  sprintf(hisnam,"hb%03ds%d",fModuleId,fSide);
  TH1F* h=new TH1F(hisnam,"",256,-0.5,255.5);
  for(Int_t ian=0;ian<fgkNAnodes;ian++){
    h->SetBinContent(ian+1,GetAnodeBaseline(ian));
  }
  return h;
}
//______________________________________________________________________
TH1F* AliITSOnlineSDDCMN::GetRawNoiseAnodeHisto() const {
  //
  Char_t hisnam[20];  
  sprintf(hisnam,"hn%03ds%d",fModuleId,fSide);
  TH1F* h=new TH1F(hisnam,"",256,-0.5,255.5);
  for(Int_t ian=0;ian<fgkNAnodes;ian++){
    h->SetBinContent(ian+1,GetAnodeRawNoise(ian));
  }
  return h;
}
//______________________________________________________________________
TH1F* AliITSOnlineSDDCMN::GetCorrNoiseAnodeHisto() const {
  //
  Char_t hisnam[20];  
  sprintf(hisnam,"hc%03ds%d",fModuleId,fSide);
  TH1F* h=new TH1F(hisnam,"",256,-0.5,255.5);
  for(Int_t ian=0;ian<fgkNAnodes;ian++){
    h->SetBinContent(ian+1,GetAnodeCorrNoise(ian));
  }
  return h;
}
//______________________________________________________________________
TH1F* AliITSOnlineSDDCMN::GetBaselineHisto() const {
  //
  Char_t hisnam[20];  
  sprintf(hisnam,"hdb%03ds%d",fModuleId,fSide);
  TH1F* h=new TH1F(hisnam,"",100,0.,150.);
  for(Int_t ian=0;ian<fgkNAnodes;ian++){
    h->Fill(GetAnodeBaseline(ian));
  }
  return h;
}
//______________________________________________________________________
TH1F* AliITSOnlineSDDCMN::GetRawNoiseHisto() const {
  //
  Char_t hisnam[20];  
  sprintf(hisnam,"hdn%03ds%d",fModuleId,fSide);
  TH1F* h=new TH1F(hisnam,"",100,0.,8.);
  for(Int_t ian=0;ian<fgkNAnodes;ian++){
    h->Fill(GetAnodeRawNoise(ian));
  }
  return h;
}
//______________________________________________________________________
TH1F* AliITSOnlineSDDCMN::GetCorrNoiseHisto() const {
  //
  Char_t hisnam[20];  
  sprintf(hisnam,"hdc%03ds%d",fModuleId,fSide);
  TH1F* h=new TH1F(hisnam,"",100,0.,8.);
  for(Int_t ian=0;ian<fgkNAnodes;ian++){
    h->Fill(GetAnodeCorrNoise(ian));
  }
  return h;
}
//______________________________________________________________________
Bool_t AliITSOnlineSDDCMN::WriteToROOT(TFile *fil){
  //
  if(fil==0){ 
    AliWarning("Invalid pointer to ROOT file");
    return kFALSE;    
  }
  Char_t hisnam[20];
  fil->cd();
  sprintf(hisnam,"hgood%03ds%d",fModuleId,fSide);
  TH1F hgood(hisnam,"",256,-0.5,255.5);
  sprintf(hisnam,"hbase%03ds%d",fModuleId,fSide);
  TH1F hbase(hisnam,"",256,-0.5,255.5);
  sprintf(hisnam,"hnois%03ds%d",fModuleId,fSide);
  TH1F hnois(hisnam,"",256,-0.5,255.5);
  sprintf(hisnam,"hcmn%03ds%d",fModuleId,fSide);
  TH1F hcmn(hisnam,"",256,-0.5,255.5);
  sprintf(hisnam,"hcorn%03ds%d",fModuleId,fSide);
  TH1F hcorn(hisnam,"",256,-0.5,255.5);
  for(Int_t ian=0;ian<fgkNAnodes;ian++){
    hgood.SetBinContent(ian+1,float(IsAnodeGood(ian)));
    hbase.SetBinContent(ian+1,GetAnodeBaseline(ian));
    hnois.SetBinContent(ian+1,GetAnodeRawNoise(ian));
    hcmn.SetBinContent(ian+1,GetAnodeCommonMode(ian));
    hcorn.SetBinContent(ian+1,GetAnodeCorrNoise(ian));
  }
  hgood.Write();
  hbase.Write();
  hnois.Write();
  hcmn.Write();
  hcorn.Write();
  return kTRUE;
}
