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
  AliITSOnlineSDDCMN::AliITSOnlineSDDCMN():AliITSOnlineSDD(),fNEvents(0),fLowThreshold(0),fHighThreshold(0),fMinCorrNoise(0.),fMaxCorrNoise(0.),fNSigmaNoise(0.)
{
  // default constructor
  Reset();
  SetMinNoise();
  SetMaxNoise();
  SetNSigmaNoise();
}
//______________________________________________________________________
  AliITSOnlineSDDCMN::AliITSOnlineSDDCMN(Int_t nddl, Int_t ncarlos, Int_t sid):AliITSOnlineSDD(nddl,ncarlos,sid),fNEvents(0),fLowThreshold(0),fHighThreshold(0),fMinCorrNoise(0.),fMaxCorrNoise(0.),fNSigmaNoise(0.)
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
  // Reset counters
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
  TString basfilnam;
  basfilnam.Form("SDDbase_step1_ddl%02dc%02d_sid%d.data",fDDL,fCarlos,fSide);
  FILE* basf=fopen(basfilnam.Data(),"r");
  if(basf==0){
    AliWarning(Form("Baseline file not present (ddl %d  carlos %d side %d, Set all baselines to 50\n",fDDL,fCarlos,fSide));
    for(Int_t ian=0;ian<fgkNAnodes;ian++){ 
      fBaseline[ian]=50.;
      fEqBaseline[ian]=50;
      fOffsetBaseline[ian]=0;
      fGoodAnode[ian]=1;
    }
    return;
  }
  Int_t check = fscanf(basf,"%d\n",&fHighThreshold);
  if(check<1)AliError("Error while reading file with baselines");
  check = fscanf(basf,"%d\n",&fLowThreshold);
  if(check<1)AliError("Error while reading file with baselines");
  Int_t n,ok,eqbase,offbase;
  Float_t base,rms,cmn,corrnoi;
  for(Int_t ian=0;ian<fgkNAnodes;ian++){
    check = fscanf(basf,"%d %d %f %d %d %f %f %f\n",&n,&ok,&base,&eqbase,&offbase,&rms,&cmn,&corrnoi);
    if(check<1)AliError("Error while reading file with baselines");
    fGoodAnode[ian]=ok;
    fBaseline[ian]=base;
    fEqBaseline[ian]=eqbase;
    fOffsetBaseline[ian]=offbase;
    fRawNoise[ian]=rms;
    fCMN[ian]=cmn;
  }
  fclose(basf);
}
//______________________________________________________________________
void  AliITSOnlineSDDCMN::ValidateAnodes(){
  // Tag good/bad anodes
  for(Int_t ian=0;ian<fgkNAnodes;ian++){
    if(!fGoodAnode[ian]) continue;
    if(GetAnodeCorrNoise(ian)>fMaxCorrNoise || GetAnodeCorrNoise(ian)<fMinCorrNoise) fGoodAnode[ian]=0;
    if(GetAnodeCorrNoise(ian)>fNSigmaNoise*CalcMeanNoise()) fGoodAnode[ian]=0;
  }
}

//______________________________________________________________________
TH2F* AliITSOnlineSDDCMN::GetCleanEvent(TH2F* hrawd) const {
  // Fills an histogram with counts corrected for common mode noise

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
  return hcorrd;
}
//______________________________________________________________________
void AliITSOnlineSDDCMN::AddEvent(TH2F* hrawd){
  // analyzes one event and adds its ontribution to the various counters

  fNEvents++;
  TH2F* hcorrd=GetCleanEvent(hrawd);

  for(Int_t ian=0;ian<fgkNAnodes;ian++){
    if(!fGoodAnode[ian]) continue;
    Float_t sumQ=0.;
    Int_t cnt=0;
    for(Int_t itb=fFirstGoodTB;itb<=fLastGoodTB;itb++){
      Float_t cntdiff=hcorrd->GetBinContent(itb+1,ian+1)-fBaseline[ian];
      sumQ+=cntdiff*cntdiff;
      cnt++;    
    }
    if(cnt != 0)fSumCorrNoise[ian]+=TMath::Sqrt(sumQ/(Float_t)cnt);
  }
  delete hcorrd;
}
//______________________________________________________________________
Float_t AliITSOnlineSDDCMN::CalcMeanNoise() const{
  // compute average noise

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
  // writes parameters of each channel into an ASCII file 
  // to be then read by the PULSER DA (AliITSOnlineSDDTP)

  TString outfilnam;
  outfilnam.Form("SDDbase_step2_ddl%02dc%02d_sid%d.data",fDDL,fCarlos,fSide);
  FILE* outf=fopen(outfilnam.Data(),"w");
  fprintf(outf,"%d\n",fHighThreshold);
  fprintf(outf,"%d\n",fLowThreshold);
  for(Int_t ian=0;ian<fgkNAnodes;ian++){
    fprintf(outf,"%d %d %8.3f %d %d %8.3f %8.3f %8.3f\n",ian,IsAnodeGood(ian),GetAnodeBaseline(ian),GetAnodeEqualizedBaseline(ian),GetAnodeBaselineOffset(ian),GetAnodeRawNoise(ian),GetAnodeCommonMode(ian),GetAnodeCorrNoise(ian));
  }
  fclose(outf);  
}

//______________________________________________________________________
TH1F* AliITSOnlineSDDCMN::GetBaselineAnodeHisto() const {
  // produce histogram with baseline vs. anode number
  TString hisnam;  
  hisnam.Form("hbase%02dc%02ds%d",fDDL,fCarlos,fSide);
  TH1F* h=new TH1F(hisnam.Data(),"",256,-0.5,255.5);
  for(Int_t ian=0;ian<fgkNAnodes;ian++){
    h->SetBinContent(ian+1,GetAnodeBaseline(ian));
  }
  return h;
}
//______________________________________________________________________
TH1F* AliITSOnlineSDDCMN::GetRawNoiseAnodeHisto() const {
  // produce histogram with raw noise vs. anode number
  TString hisnam;  
  hisnam.Form("hnois%02dc%02ds%d",fDDL,fCarlos,fSide);
  TH1F* h=new TH1F(hisnam.Data(),"",256,-0.5,255.5);
  for(Int_t ian=0;ian<fgkNAnodes;ian++){
    h->SetBinContent(ian+1,GetAnodeRawNoise(ian));
  }
  return h;
}
//______________________________________________________________________
TH1F* AliITSOnlineSDDCMN::GetCorrNoiseAnodeHisto() const {
  // produce histogram with corrected noise vs. anode number
  TString hisnam;  
  hisnam.Form("hcorn%02dc%02ds%d",fDDL,fCarlos,fSide);
  TH1F* h=new TH1F(hisnam.Data(),"",256,-0.5,255.5);
  for(Int_t ian=0;ian<fgkNAnodes;ian++){
    h->SetBinContent(ian+1,GetAnodeCorrNoise(ian));
  }
  return h;
}
//______________________________________________________________________
TH1F* AliITSOnlineSDDCMN::GetCMNCoefAnodeHisto() const {
  // produce histogram with coefficients for common mode noise subtraction
  TString hisnam;  
  hisnam.Form("hcmn%02dc%02ds%d",fDDL,fCarlos,fSide);
  TH1F* h=new TH1F(hisnam.Data(),"",256,-0.5,255.5);
  for(Int_t ian=0;ian<fgkNAnodes;ian++){
    h->SetBinContent(ian+1,GetAnodeCommonMode(ian));
  }
  return h;
}
//______________________________________________________________________
TH1F* AliITSOnlineSDDCMN::GetStatusAnodeHisto() const {
  // produce histogram with status bit of each anode
  TString hisnam;  
  hisnam.Form("hgood%02dc%02ds%d",fDDL,fCarlos,fSide);
  TH1F* h=new TH1F(hisnam.Data(),"",256,-0.5,255.5);
  for(Int_t ian=0;ian<fgkNAnodes;ian++){
    h->SetBinContent(ian+1,float(IsAnodeGood(ian)));
  }
  return h;
}
//______________________________________________________________________
TH1F* AliITSOnlineSDDCMN::GetBaselineHisto() const {
  // produce histogram with baseline distribution
  TString hisnam;  
  hisnam.Form("hdbd%02dc%02ds%d",fDDL,fCarlos,fSide);
  TH1F* h=new TH1F(hisnam.Data(),"",100,0.,150.);
  for(Int_t ian=0;ian<fgkNAnodes;ian++){
    h->Fill(GetAnodeBaseline(ian));
  }
  return h;
}
//______________________________________________________________________
TH1F* AliITSOnlineSDDCMN::GetRawNoiseHisto() const {
  // produce histogram with raw noise distribution
  TString hisnam;  
  hisnam.Form("hdnd%02dc%02ds%d",fDDL,fCarlos,fSide);
  TH1F* h=new TH1F(hisnam.Data(),"",100,0.,8.);
  for(Int_t ian=0;ian<fgkNAnodes;ian++){
    h->Fill(GetAnodeRawNoise(ian));
  }
  return h;
}
//______________________________________________________________________
TH1F* AliITSOnlineSDDCMN::GetCorrNoiseHisto() const {
  // produce histogram with corrected noise distribution
  TString hisnam;  
  hisnam.Form("hdcd%02dc%02ds%d",fDDL,fCarlos,fSide);
  TH1F* h=new TH1F(hisnam.Data(),"",100,0.,8.);
  for(Int_t ian=0;ian<fgkNAnodes;ian++){
    h->Fill(GetAnodeCorrNoise(ian));
  }
  return h;
}
//______________________________________________________________________
Bool_t AliITSOnlineSDDCMN::WriteToROOT(TFile *fil){
  // writes output into a root file
  if(fil==0){ 
    AliWarning("Invalid pointer to ROOT file");
    return kFALSE;    
  }
  TString hisnam;
  fil->cd();
  hisnam.Form("hgood%02dc%02ds%d",fDDL,fCarlos,fSide);
  TH1F hgood(hisnam.Data(),"",256,-0.5,255.5);
  hisnam.Form("hbase%02dc%02ds%d",fDDL,fCarlos,fSide);
  TH1F hbase(hisnam.Data(),"",256,-0.5,255.5);
  hisnam.Form("hnois%02dc%02ds%d",fDDL,fCarlos,fSide);
  TH1F hnois(hisnam.Data(),"",256,-0.5,255.5);
  hisnam.Form("hcmn%02dc%02ds%d",fDDL,fCarlos,fSide);
  TH1F hcmn(hisnam.Data(),"",256,-0.5,255.5);
  hisnam.Form("hcorn%02dc%02ds%d",fDDL,fCarlos,fSide);
  TH1F hcorn(hisnam.Data(),"",256,-0.5,255.5);
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
