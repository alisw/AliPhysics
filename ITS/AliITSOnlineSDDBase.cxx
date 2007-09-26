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
#include "AliITSOnlineSDDBase.h"
#include <TH2F.h>
#include <TMath.h>


///////////////////////////////////////////////////////////////////
//                                                               //
// Implementation of the class used for SDD baselines            //
// and noise analysis                                            //
// Origin: F.Prino, Torino, prino@to.infn.it                     //
//                                                               //
///////////////////////////////////////////////////////////////////


ClassImp(AliITSOnlineSDDBase)
//______________________________________________________________________
  AliITSOnlineSDDBase::AliITSOnlineSDDBase():AliITSOnlineSDD(),fNEvents(0),fMinBaseline(0.),fMaxBaseline(0.),fMinRawNoise(0.),fMaxRawNoise(0.),fNSigmaNoise(0.)
{
  // default constructor
  Reset();
  SetMinBaseline();
  SetMaxBaseline();
  SetMinRawNoise();
  SetMaxRawNoise();
  SetNSigmaNoise();
}
//______________________________________________________________________
  AliITSOnlineSDDBase::AliITSOnlineSDDBase(Int_t mod, Int_t sid):AliITSOnlineSDD(mod,sid),fNEvents(0),fMinBaseline(0.),fMaxBaseline(0.),fMinRawNoise(0.),fMaxRawNoise(0.),fNSigmaNoise(0.)
{
  // default constructor
  Reset();
  SetMinBaseline();
  SetMaxBaseline();
  SetMinRawNoise();
  SetMaxRawNoise();
  SetNSigmaNoise();
}
//______________________________________________________________________
AliITSOnlineSDDBase::~AliITSOnlineSDDBase(){
  // Destructor
}
//______________________________________________________________________
void AliITSOnlineSDDBase::Reset(){
  //
  fNEvents=0;
  for(Int_t i=0;i<fgkNAnodes;i++){
    fGoodAnode[i]=1;
    fSumBaseline[i]=0.;
    fSumRawNoise[i]=0.;
    fSumCMN[i]=0.;
  }
}
//______________________________________________________________________
void  AliITSOnlineSDDBase::ValidateAnodes(){
  //
  for(Int_t ian=0;ian<fgkNAnodes;ian++){
    fGoodAnode[ian]=1;
    if(GetAnodeBaseline(ian)>fMaxBaseline || GetAnodeBaseline(ian)<fMinBaseline) fGoodAnode[ian]=0;
    if(GetAnodeRawNoise(ian)>fMaxRawNoise || GetAnodeRawNoise(ian)<fMinRawNoise) fGoodAnode[ian]=0;
    if(GetAnodeRawNoise(ian)>fNSigmaNoise*CalcMeanRawNoise()) fGoodAnode[ian]=0;
  }
}

//______________________________________________________________________
void AliITSOnlineSDDBase::AddEvent(TH2F* hrawd){
  // 
  fNEvents++;
  const Int_t kTimeBins=fLastGoodTB-fFirstGoodTB+1;
  Float_t sum[fgkNAnodes];
  for(Int_t ian=0;ian<fgkNAnodes;ian++){
    Float_t sumQ=0.;
    sum[ian]=0.;
    for(Int_t itb=fFirstGoodTB;itb<=fLastGoodTB;itb++){
      sum[ian]+=hrawd->GetBinContent(itb+1,ian+1);
      sumQ+=TMath::Power(hrawd->GetBinContent(itb+1,ian+1),2);      
    }
    sum[ian]/=(Float_t)kTimeBins;
    sumQ/=(Float_t)kTimeBins;
    fSumBaseline[ian]+=sum[ian];
    fSumRawNoise[ian]+=sumQ;
    if(fNEvents==1) ValidateAnodes();
  }


  Float_t *cmnEven = new Float_t[kTimeBins];
  Float_t *cmnOdd  = new Float_t[kTimeBins];
  for(Int_t itb=fFirstGoodTB;itb<=fLastGoodTB;itb++){
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
    if(countEven>0) cmnEven[itb]=sumEven/countEven;
    if(countOdd>0) cmnOdd[itb]=sumOdd/countOdd;
  }
  for(Int_t ian=0;ian<fgkNAnodes;ian++){
    Float_t num=0.,den=0.;
    if(!fGoodAnode[ian]) continue;
    for(Int_t itb=fFirstGoodTB;itb<=fLastGoodTB;itb++){
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
Float_t AliITSOnlineSDDBase::CalcMeanRawNoise() const{
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
void AliITSOnlineSDDBase::WriteToASCII(){
  //
  Char_t outfilnam[100];
  sprintf(outfilnam,"SDDbase_step1_mod%03d_sid%d.data",fModuleId,fSide);
  FILE* outf=fopen(outfilnam,"w");
  Float_t corrnoise=2.;
  for(Int_t ian=0;ian<fgkNAnodes;ian++){
    fprintf(outf,"%d %d %11.6f %11.6f %11.6f %11.6f\n",ian,IsAnodeGood(ian),GetAnodeBaseline(ian),GetAnodeRawNoise(ian),GetAnodeCommonMode(ian),corrnoise);
  }
  fclose(outf);  
}
