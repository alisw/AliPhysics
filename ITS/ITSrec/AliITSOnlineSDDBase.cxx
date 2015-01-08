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

/*  $Id$   */

const Int_t AliITSOnlineSDDBase::fgkMaxCorr=63; // 6 but correction

ClassImp(AliITSOnlineSDDBase)
//______________________________________________________________________
  AliITSOnlineSDDBase::AliITSOnlineSDDBase():AliITSOnlineSDD(),fNEvents(0),fMinBaseline(0.),fMaxBaseline(0.),fMinRawNoise(0.),fMaxRawNoise(0.),fNSigmaNoise(0.),fGoldenBaseline(0.),fLowThrFact(0.),fHighThrFact(0.)
{
  // default constructor
  Reset();
  SetMinBaseline();
  SetMaxBaseline();
  SetMinRawNoise();
  SetMaxRawNoise();
  SetNSigmaNoise();
  SetGoldenBaselineValue();
  SetZeroSuppThresholds();
}
//______________________________________________________________________
AliITSOnlineSDDBase::AliITSOnlineSDDBase(Int_t nddl, Int_t ncarlos, Int_t sid):AliITSOnlineSDD(nddl,ncarlos,sid),fNEvents(0),fMinBaseline(0.),fMaxBaseline(0.),fMinRawNoise(0.),fMaxRawNoise(0.),fNSigmaNoise(0.),fGoldenBaseline(0.),fLowThrFact(0.),fHighThrFact(0.)
{
  // default constructor
  Reset();
  SetMinBaseline();
  SetMaxBaseline();
  SetMinRawNoise();
  SetMaxRawNoise();
  SetNSigmaNoise();
  SetGoldenBaselineValue();
  SetZeroSuppThresholds();
}
//______________________________________________________________________
AliITSOnlineSDDBase::~AliITSOnlineSDDBase(){
  // Destructor
}
//______________________________________________________________________
void AliITSOnlineSDDBase::Reset(){
  // reset all counters
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
  // tag good/bad channels
  for(Int_t ian=0;ian<fgkNAnodes;ian++){
    fGoodAnode[ian]=1;
    Float_t basel=GetAnodeBaseline(ian);
    Float_t rawn=GetAnodeRawNoise(ian);
    Float_t ratio=0.;
    if(rawn>0) ratio=basel/rawn;
    if(basel>fMaxBaseline || basel<fMinBaseline) fGoodAnode[ian]=0;
    else if(rawn>fMaxRawNoise || rawn<fMinRawNoise) fGoodAnode[ian]=0;
    else if(rawn>fNSigmaNoise*CalcMeanRawNoise()) fGoodAnode[ian]=0;
    else if(ratio<3.) fGoodAnode[ian]=0;
  }
}

//______________________________________________________________________
void AliITSOnlineSDDBase::AddEvent(TH2F* hrawd){
  // analyzes one event and adds its ontribution to the various counters

  fNEvents++;
  const Int_t kTimeBins=fLastGoodTB+1;
  Float_t sum[fgkNAnodes];
  for(Int_t ian=0;ian<fgkNAnodes;ian++){
    Float_t sumQ=0.;
    sum[ian]=0.;
    Int_t cnt=0;
    for(Int_t itb=fFirstGoodTB;itb<=fLastGoodTB;itb++){
      Float_t cbin=hrawd->GetBinContent(itb+1,ian+1);
      sum[ian]+=cbin;
      sumQ+=cbin*cbin;
      cnt++;
    }
    if(cnt != 0){
      sum[ian]/=(Float_t)cnt;
      sumQ/=(Float_t)cnt;
    }
    fSumBaseline[ian]+=sum[ian];
    fSumRawNoise[ian]+=sumQ;
  }
  if(fNEvents==1) ValidateAnodes();

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
      den+=cmnCoef*cmnCoef;
    }
    if(den!=0) fSumCMN[ian]+=num/den;
  }

  delete [] cmnEven;
  delete [] cmnOdd;
}
//______________________________________________________________________
void AliITSOnlineSDDBase::GetMinAndMaxBaseline(Float_t &basMin, Float_t &basMax) const {
  // fills mininum and maximum baseline values
  basMin=1008.;
  basMax=0.;
  for(Int_t ian=0;ian<fgkNAnodes;ian++){
    if(!fGoodAnode[ian]) continue;
    Float_t bas=GetAnodeBaseline(ian);
    if(bas>0 && bas < basMin) basMin=bas;
    if(bas>0 && bas > basMax) basMax=bas;
  }
}
//______________________________________________________________________
Float_t AliITSOnlineSDDBase::GetMinimumBaseline() const {
  // returns anode with minum baseline value in hybrid
  Float_t basMin=1008.;
  for(Int_t ian=0;ian<fgkNAnodes;ian++){
    if(!fGoodAnode[ian]) continue;
    Float_t bas=GetAnodeBaseline(ian);
    if(bas>0 && bas < basMin) basMin=bas;
  }
  return basMin;
}
//______________________________________________________________________
Float_t AliITSOnlineSDDBase::CalcMeanRawNoise() const{
  // compute mean value of raw noise
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
  // writes parameters of each channel into an ASCII file 
  // to be then read in the successive step for common mode noise
  // correction (AliITSOnlineSDDCMN)

  TString outfilnam;
  Float_t basMin,basMax;
  GetMinAndMaxBaseline(basMin,basMax);
  Float_t finalVal=basMin;
  if(basMin>fGoldenBaseline && basMax<fGoldenBaseline+fgkMaxCorr) finalVal=fGoldenBaseline;
  if(basMax<basMin+fgkMaxCorr && basMax>fGoldenBaseline+fgkMaxCorr) finalVal=basMax-fgkMaxCorr;

  Float_t avNoise=CalcMeanRawNoise();
  Int_t thrL=(Int_t)(finalVal+fLowThrFact*avNoise+0.5);
  Int_t thrH=(Int_t)(finalVal+fHighThrFact*avNoise+0.5);
  if(CountGoodAnodes()==0) thrH=255;

  outfilnam.Form("SDDbase_step1_ddl%02dc%02d_sid%d.data",fDDL,fCarlos,fSide);
  FILE* outf=fopen(outfilnam.Data(),"w");
  fprintf(outf,"%d\n",thrH);
  fprintf(outf,"%d\n",thrL);
  Float_t corrnoise=2.;
  for(Int_t ian=0;ian<fgkNAnodes;ian++){
    Float_t bas=GetAnodeBaseline(ian);
    Int_t corr=(Int_t)(bas-finalVal+0.5);
    if(corr>fgkMaxCorr) corr=fgkMaxCorr; // only 6 bits in jtag for correction
    if(corr<0) corr=0; // avoid negative numbers
    fprintf(outf,"%d %d %11.6f %d %d %11.6f %11.6f %11.6f\n",ian,IsAnodeGood(ian),GetAnodeBaseline(ian),(Int_t)finalVal,corr,GetAnodeRawNoise(ian),GetAnodeCommonMode(ian),corrnoise);
  }
  fclose(outf);  
}
