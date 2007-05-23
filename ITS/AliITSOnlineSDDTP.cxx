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
#include "AliITSOnlineSDDTP.h"
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
AliITSOnlineSDDTP::AliITSOnlineSDDTP():AliITSOnlineSDD(),fNEvents(0),fDAQ(0.),fNSigmaGain(0.)
{
  // default constructor
  Reset();
  SetNSigmaGain();
}
//______________________________________________________________________
AliITSOnlineSDDTP::AliITSOnlineSDDTP(Int_t mod, Int_t sid, Float_t xDAQ):AliITSOnlineSDD(mod,sid),fNEvents(0),fDAQ(xDAQ),fNSigmaGain(0.)
{
  // standard constructor
  Reset();
  SetNSigmaGain();
}
//______________________________________________________________________
AliITSOnlineSDDTP::~AliITSOnlineSDDTP(){
  // Destructor
}
//______________________________________________________________________
void AliITSOnlineSDDTP::Reset(){
  fNEvents=0;
  for(Int_t i=0;i<fgkNAnodes;i++){
    fGoodAnode[i]=1;
    fBaseline[i]=0.;
    fSumTPPeak[i]=0.;
    fTPPos[i]=0.;
  }
  ReadBaselines();
}

//______________________________________________________________________
void AliITSOnlineSDDTP::AddEvent(TH2F* hrawd){
  // 
  fNEvents++;
  Double_t tbmax=(Double_t)hrawd->GetNbinsX();
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
    fSumTPPeak[ian]+=auxmax-fBaseline[ian];
    fTPPos[ian]+=auxtb;
  }
}
//______________________________________________________________________
void AliITSOnlineSDDTP::ReadBaselines(){
  // assume baselines and good anodes are taken from previous run
  Char_t basfilnam[100];
  sprintf(basfilnam,"SDDbase_step1_mod%03d_sid%d.data",fModuleId,fSide);
  FILE* basf=fopen(basfilnam,"r");
  if(basf==0){
    printf("Baselinefile not present, Set all baselines to 50\n");
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
    fBaseline[ian]=base;
    fGoodAnode[ian]=ok;
  }
  fclose(basf);
}

//______________________________________________________________________
void AliITSOnlineSDDTP::ValidateAnodes(){
  Float_t meang,rmsg;
  StatGain(meang,rmsg);
  printf("<gain>=%f,rms=%f\n",meang,rmsg);
  Float_t lowlim=meang-fNSigmaGain*rmsg;
  Float_t hilim=meang+fNSigmaGain*rmsg;

  for(Int_t ian=0;ian<fgkNAnodes;ian++){
    if(!fGoodAnode[ian]) continue;
    if(GetChannelGain(ian)<lowlim||GetChannelGain(ian)>hilim) fGoodAnode[ian]=0;
  }
}


//______________________________________________________________________
void AliITSOnlineSDDTP::StatGain(Float_t &mean, Float_t  &rms){
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
void AliITSOnlineSDDTP::WriteToFXS(){
  //
  Char_t basfilnam[100];
  sprintf(basfilnam,"SDDbase_step1_mod%03d_sid%d.data",fModuleId,fSide);
  FILE* basf=fopen(basfilnam,"r");
  Int_t n,ok;
  Float_t base,rms,cmn,corrnoi;
  Float_t noise[fgkNAnodes],cmncoef[fgkNAnodes],corrnoise[fgkNAnodes];
  for(Int_t ian=0;ian<fgkNAnodes;ian++){
    fscanf(basf,"%d %d %f %f %f %f\n",&n,&ok,&base,&rms,&cmn,&corrnoi);
    noise[ian]=rms;
    cmncoef[ian]=cmn;
    corrnoise[ian]=corrnoi;
  }
  fclose(basf);
  printf("Read All******************\n");
  Char_t outfilnam[100];
  sprintf(outfilnam,"SDDbase_mod%03d_sid%d.data",fModuleId,fSide);
  FILE* outf=fopen(outfilnam,"w");
  for(Int_t ian=0;ian<fgkNAnodes;ian++){
    fprintf(outf,"%d %d %8.3f %8.3f %8.3f %8.3f %8.3f\n",ian,IsAnodeGood(ian),fBaseline[ian], noise[ian],cmncoef[ian],corrnoise[ian],GetChannelGain(ian));
  }
  fclose(outf);  
}

