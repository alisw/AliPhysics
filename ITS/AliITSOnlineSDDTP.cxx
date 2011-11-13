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
#include "AliITSOnlineSDDTP.h"
#include "AliLog.h"
#include <TH2F.h>
#include <TMath.h>


///////////////////////////////////////////////////////////////////
//                                                               //
// Implemetation of the class SDD Test Pulse analysis            //
// Origin: F.Prino, Torino, prino@to.infn.it                     //
//                                                               //
///////////////////////////////////////////////////////////////////


ClassImp(AliITSOnlineSDDTP)

//______________________________________________________________________
AliITSOnlineSDDTP::AliITSOnlineSDDTP():AliITSOnlineSDD(),fDAC(0.),fLowThreshold(0),fHighThreshold(0),fNSigmaGain(0.),fNSigmaNoise(0.)
{
  // default constructor
  Reset();
  SetNSigmaGain();
  SetNSigmaNoise();
}
//______________________________________________________________________
AliITSOnlineSDDTP::AliITSOnlineSDDTP(Int_t nddl, Int_t ncarlos, Int_t sid, Float_t xDAC):AliITSOnlineSDD(nddl,ncarlos,sid),fDAC(xDAC),fLowThreshold(0),fHighThreshold(0),fNSigmaGain(0.),fNSigmaNoise(0.)
{
  // standard constructor
  Reset();
  SetNSigmaGain();
  SetNSigmaNoise();
}
//______________________________________________________________________
AliITSOnlineSDDTP::~AliITSOnlineSDDTP(){
  // Destructor
}
//______________________________________________________________________
void AliITSOnlineSDDTP::Reset(){
  // reset all counters
  for(Int_t i=0;i<fgkNAnodes;i++){
    fNEvents[i]=0;
    fGoodAnode[i]=1;
    fBaseline[i]=0.;
    fCMN[i]=0.;
    fRawNoise[i]=0.;
    fCorrNoise[i]=0.;
    fSumTPPeak[i]=0.;
    fTPPos[i]=0.;
  }
  ReadBaselines();
}

//______________________________________________________________________
void AliITSOnlineSDDTP::AddEvent(TH2F* hrawd){
  // analyzes current event and sum its contribution to the various counters
  for(Int_t ian=0;ian<fgkNAnodes;ian++){
    Float_t auxmax=0.;
    Int_t auxtb=0;
    if(!fGoodAnode[ian]) continue;
    for(Int_t itb=fFirstGoodTB;itb<=fLastGoodTB;itb++){
      Float_t cnt=hrawd->GetBinContent(itb+1,ian+1);
      if(cnt>auxmax){ 
	auxmax=cnt;
	auxtb=itb;
      }
    }
    if(auxmax>fBaseline[ian]+fNSigmaNoise*fRawNoise[ian]){
      fSumTPPeak[ian]+=auxmax-fBaseline[ian];
      fTPPos[ian]+=auxtb;
      fNEvents[ian]++;
    }
  }
}
//______________________________________________________________________
void AliITSOnlineSDDTP::ReadBaselines(){
  // assume baselines and good anodes are taken from previous run
  TString basfilnam;
  basfilnam.Form("SDDbase_step2_ddl%02dc%02d_sid%d.data",fDDL,fCarlos,fSide);
  FILE* basf=fopen(basfilnam.Data(),"r");
  if(basf==0){
    AliWarning(Form("Baseline file not present (ddl %d  carlos %d side %d, Set all baselines to 20",fDDL,fCarlos,fSide));
    for(Int_t ian=0;ian<fgkNAnodes;ian++){ 
      fBaseline[ian]=20.;
      fEqBaseline[ian]=20;
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
    fBaseline[ian]=base;
    fEqBaseline[ian]=eqbase;
    fOffsetBaseline[ian]=offbase;
    fGoodAnode[ian]=ok;
    fRawNoise[ian]=rms;
    fCMN[ian]=cmn;
    fCorrNoise[ian]=corrnoi;
  }
  fclose(basf);
}

//______________________________________________________________________
Bool_t AliITSOnlineSDDTP::IsModuleGood() const{
  //
  // Check if there is at least 1 good anode
  //
  for(Int_t ian=0;ian<fgkNAnodes;ian++){
    if(fGoodAnode[ian]) return kTRUE;
  }
  return kFALSE;
}
//______________________________________________________________________
void AliITSOnlineSDDTP::ValidateAnodes(){
  // tag good/bad channels
  Float_t meang,rmsg;
  StatGain(meang,rmsg);
  Float_t lowlim=meang-fNSigmaGain*rmsg;
  Float_t hilim=meang+fNSigmaGain*rmsg;

  for(Int_t ian=0;ian<fgkNAnodes;ian++){
    if(!fGoodAnode[ian]) continue;
    if(GetChannelGain(ian)<lowlim||GetChannelGain(ian)>hilim) fGoodAnode[ian]=0;
  }
}


//______________________________________________________________________
void AliITSOnlineSDDTP::StatGain(Float_t &mean, Float_t  &rms) const {
  // compute average gain and rms
  Float_t sum=0.,sumq=0.;
  Int_t cnt=0;
  for(Int_t ian=0;ian<fgkNAnodes;ian++){
    if(!fGoodAnode[ian]) continue;
    if(fNEvents[ian]==0) continue;
    Float_t chgain=GetChannelGain(ian);
    sum+=chgain;
    sumq+=chgain*chgain;
    cnt++;
  }
  if(cnt>0){ 
    mean=sum/(Float_t)cnt;
    Float_t variance=sumq/(Float_t)cnt-mean*mean;
    if(variance>0.) rms=TMath::Sqrt(variance);
    else rms = 0;
  }else{ 
    mean=0.;
    rms=0.;
  }
  return;
}

//______________________________________________________________________
void AliITSOnlineSDDTP::WriteToASCII(){
  // writes parameters of each channel into an ASCII file 
  // to be sent to FXS by the DA and processed by the SHUTTLE

  TString outfilnam;
  outfilnam.Form("SDDbase_ddl%02dc%02d_sid%d.data",fDDL,fCarlos,fSide);
  FILE* outf=fopen(outfilnam.Data(),"w");
  fprintf(outf,"%d %d %d\n",fCarlos,fSide,IsModuleGood());
  fprintf(outf,"%d\n",fHighThreshold);
  fprintf(outf,"%d\n",fLowThreshold);
  for(Int_t ian=0;ian<fgkNAnodes;ian++){
    fprintf(outf,"%d %d %8.3f %d %d %8.3f %8.3f %8.3f %8.3f\n",ian,IsAnodeGood(ian),GetAnodeBaseline(ian),GetAnodeEqualizedBaseline(ian),GetAnodeBaselineOffset(ian),GetAnodeRawNoise(ian),GetAnodeCommonMode(ian),GetAnodeCorrNoise(ian),GetChannelGain(ian));
  }
  fclose(outf);  
}
//______________________________________________________________________
TH1F* AliITSOnlineSDDTP::GetBaselineAnodeHisto() const {
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
TH1F* AliITSOnlineSDDTP::GetRawNoiseAnodeHisto() const {
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
TH1F* AliITSOnlineSDDTP::GetCorrNoiseAnodeHisto() const {
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
TH1F* AliITSOnlineSDDTP::GetCMNCoefAnodeHisto() const {
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
TH1F* AliITSOnlineSDDTP::GetStatusAnodeHisto() const {
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
TH1F* AliITSOnlineSDDTP::GetGainAnodeHisto() const {
  // produce histogram with gain vs. anode number
  TString hisnam;  
  hisnam.Form("hgain%02dc%02ds%d",fDDL,fCarlos,fSide);
  TH1F* h=new TH1F(hisnam.Data(),"",256,-0.5,255.5);
  for(Int_t ian=0;ian<fgkNAnodes;ian++){
    h->SetBinContent(ian+1,GetChannelGain(ian));
  }
  return h;
}
//______________________________________________________________________
Bool_t AliITSOnlineSDDTP::WriteToROOT(TFile *fil){
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
  hisnam.Form("hgain%02dc%02ds%d",fDDL,fCarlos,fSide);
  TH1F hgain(hisnam.Data(),"",256,-0.5,255.5);
  hisnam.Form("htptb%02dc%02ds%d",fDDL,fCarlos,fSide);
  TH1F htptb(hisnam.Data(),"",256,-0.5,255.5);
  for(Int_t ian=0;ian<fgkNAnodes;ian++){
    hgood.SetBinContent(ian+1,float(IsAnodeGood(ian)));
    hbase.SetBinContent(ian+1,GetAnodeBaseline(ian));
    hnois.SetBinContent(ian+1,GetAnodeRawNoise(ian));
    hcmn.SetBinContent(ian+1,GetAnodeCommonMode(ian));
    hcorn.SetBinContent(ian+1,GetAnodeCorrNoise(ian));
    hgain.SetBinContent(ian+1,GetChannelGain(ian));
    htptb.SetBinContent(ian+1,GetTimeBinTPPeak(ian));
  }
  hgood.Write();
  hbase.Write();
  hnois.Write();
  hcmn.Write();
  hcorn.Write();
  hgain.Write();
  htptb.Write();
  return kTRUE;
}

