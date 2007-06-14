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
#include "AliITSOnlineSDDBTP.h"
#include <TH2F.h>
#include <TMath.h>


///////////////////////////////////////////////////////////////////
//                                                               //
// Class used for SDD baseline, noise and gain analysis          //
// Origin: F.Prino, Torino, prino@to.infn.it                     //
//                                                               //
///////////////////////////////////////////////////////////////////


ClassImp(AliITSOnlineSDDBTP)
//______________________________________________________________________
  AliITSOnlineSDDBTP::AliITSOnlineSDDBTP():AliITSOnlineSDD(),fNBaseEvents(0),fNTPEvents(0),fMinBaseline(0.),fMaxBaseline(0.),fMinRawNoise(0.),fMaxRawNoise(0.),fNSigmaNoise(0.),fNSigmaGain(0.)
{
  // default constructor
  Reset();
  SetMinBaseline();
  SetMaxBaseline();
  SetMinRawNoise();
  SetMaxRawNoise();
  SetNSigmaNoise();
  SetNSigmaGain();
}
//______________________________________________________________________
  AliITSOnlineSDDBTP::AliITSOnlineSDDBTP(Int_t mod, Int_t sid):AliITSOnlineSDD(mod,sid),fNBaseEvents(0),fNTPEvents(0),fMinBaseline(0.),fMaxBaseline(0.),fMinRawNoise(0.),fMaxRawNoise(0.),fNSigmaNoise(0.),fNSigmaGain(0.)
{
  // default constructor
  Reset();
  SetMinBaseline();
  SetMaxBaseline();
  SetMinRawNoise();
  SetMaxRawNoise();
  SetNSigmaNoise();
  SetNSigmaGain();
}
//______________________________________________________________________
AliITSOnlineSDDBTP::~AliITSOnlineSDDBTP(){
  // Destructor
}
//______________________________________________________________________
void AliITSOnlineSDDBTP::Reset(){
  //
  fNBaseEvents=0;
  fNTPEvents=0;
  for(Int_t i=0;i<fgkNAnodes;i++){
    fGoodAnode[i]=1;
    fSumBaseline[i]=0.;
    fSumRawNoise[i]=0.;
    fSumCMN[i]=0.;
    fSumTPPeak[i]=0.;
    fTPPos[i]=0.;
  }
}
//______________________________________________________________________
void  AliITSOnlineSDDBTP::ValidateAnodes(){
  //
  for(Int_t ian=0;ian<fgkNAnodes;ian++){
    fGoodAnode[ian]=1;
    if(GetAnodeBaseline(ian)>fMaxBaseline || GetAnodeBaseline(ian)<fMinBaseline) fGoodAnode[ian]=0;
    if(GetAnodeRawNoise(ian)>fMaxRawNoise || GetAnodeRawNoise(ian)<fMinRawNoise) fGoodAnode[ian]=0;
    if(GetAnodeRawNoise(ian)>fNSigmaNoise*CalcMeanRawNoise()) fGoodAnode[ian]=0;
  }
  if(fNTPEvents>0){
    Float_t meang,rmsg;
    StatGain(meang,rmsg);
    Float_t lowlim=meang-fNSigmaGain*rmsg;
    Float_t hilim=meang+fNSigmaGain*rmsg;
    for(Int_t ian=0;ian<fgkNAnodes;ian++){
      if(!fGoodAnode[ian]) continue;
      if(GetChannelGain(ian)<lowlim||GetChannelGain(ian)>hilim) fGoodAnode[ian]=0;
    }
  }
}

//______________________________________________________________________
void AliITSOnlineSDDBTP::AddTPEvent(TH2F* hrawd, Float_t xDAC){
  // 
  if(fNBaseEvents==0) return;
  fNTPEvents++;
  Float_t tbmax=(Float_t)hrawd->GetNbinsX();
  for(Int_t ian=0;ian<fgkNAnodes;ian++){
    Float_t auxmax=0.;
    Int_t auxtb=0;
    if(!fGoodAnode[ian]) continue;
    for(Int_t itb=0;itb<tbmax;itb++){
      Float_t cnt=hrawd->GetBinContent(itb+1,ian+1);
      if(cnt>auxmax){ 
	auxmax=cnt;
	auxtb=itb;
      }
    }
    fSumTPPeak[ian]+=(auxmax-GetAnodeBaseline(ian))/xDAC;
    fTPPos[ian]+=auxtb;
  }
}
//______________________________________________________________________
void AliITSOnlineSDDBTP::AddBaseEvent(TH2F* hrawd){
  // 
  fNBaseEvents++;
  Float_t tbmax=(Float_t)hrawd->GetNbinsX();
  Float_t sum[fgkNAnodes];
  for(Int_t ian=0;ian<fgkNAnodes;ian++){
    Float_t sumQ=0.;
    sum[ian]=0.;
    for(Int_t itb=0;itb<tbmax;itb++){
      sum[ian]+=hrawd->GetBinContent(itb+1,ian+1);
      sumQ+=TMath::Power(hrawd->GetBinContent(itb+1,ian+1),2);      
    }
    sum[ian]/=tbmax;
    sumQ/=tbmax;
    fSumBaseline[ian]+=sum[ian];
    fSumRawNoise[ian]+=sumQ;
    if(fNBaseEvents==1) ValidateAnodes();
  }


  const Int_t kTbmax=int(tbmax);
  Float_t *cmnEven = new Float_t[kTbmax];
  Float_t *cmnOdd  = new Float_t[kTbmax];
  for(Int_t itb=0;itb<tbmax;itb++){
    Float_t sumEven=0., sumOdd=0.;
    Int_t countEven=0,countOdd=0;
    for(Int_t ian=0;ian<fgkNAnodes;ian+=2){
      if(!fGoodAnode[ian]) continue;
      sumEven+=hrawd->GetBinContent(itb+1,ian+1)-sum[ian];
      countEven++;
    }
    for(Int_t ian=1;ian<fgkNAnodes;ian+=2){
      if(!fGoodAnode[ian]) continue;
      sumOdd+=hrawd->GetBinContent(itb+1,ian+1)-sum[ian];
      countOdd++;
    }
    cmnEven[itb]=sumEven/countEven;
    cmnOdd[itb]=sumOdd/countOdd;
  }
  for(Int_t ian=0;ian<fgkNAnodes;ian++){
    Float_t num=0.,den=0.;
    if(!fGoodAnode[ian]) continue;
    for(Int_t itb=0;itb<tbmax;itb++){
      Float_t cmnCoef=cmnOdd[itb];
      if(ian%2==0) cmnCoef=cmnEven[itb];
      num+=(hrawd->GetBinContent(itb+1,ian+1)-sum[ian])*cmnCoef;
      den+=TMath::Power(cmnCoef,2);
    }
    if(den!=0) fSumCMN[ian]+=num/den;
  }
  delete [] cmnEven;
  delete [] cmnOdd;
}
//______________________________________________________________________
Float_t AliITSOnlineSDDBTP::CalcMeanRawNoise() const{
  //
  Float_t meanns=0.;
  Int_t cnt=0;
  for(Int_t ian=0;ian<fgkNAnodes;ian++){
    if(!fGoodAnode[ian]) continue;  
    meanns+=GetAnodeRawNoise(ian);
    cnt++;
  }
  if(cnt>0) meanns/=(Float_t)cnt;
  return meanns;
}
//______________________________________________________________________
void AliITSOnlineSDDBTP::StatGain(Float_t &mean, Float_t  &rms){
  //
  Float_t sum=0.,sumq=0.;
  Int_t cnt=0;
  for(Int_t ian=0;ian<fgkNAnodes;ian++){
    if(!fGoodAnode[ian]) continue;
    sum+=GetChannelGain(ian);
    sumq+=TMath::Power(GetChannelGain(ian),2);
    cnt++;
  }
  if(cnt>0){ 
    mean=sum/(Float_t)cnt;
    rms=TMath::Sqrt(sumq/(Float_t)cnt-mean*mean);
  }else{ 
    mean=0.;
    rms=0.;
  }
  return;
}
//______________________________________________________________________
void AliITSOnlineSDDBTP::WriteToASCII(){
  //
  Char_t outfilnam[100];
  sprintf(outfilnam,"SDDbase_mod%03d_sid%d.data",fModuleId,fSide);
  FILE* outf=fopen(outfilnam,"w");
  Float_t corrnoise=2.;
  for(Int_t ian=0;ian<fgkNAnodes;ian++){
    fprintf(outf,"%d %d %8.3f %8.3f %8.3f %8.3f %8.3f\n",ian,IsAnodeGood(ian),GetAnodeBaseline(ian),GetAnodeRawNoise(ian),GetAnodeCommonMode(ian),corrnoise,GetChannelGain(ian));
  }
  fclose(outf);  
}
