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
#include "AliITSOnlineSDDInjectors.h"
#include <TH2F.h>
#include <TGraphErrors.h>
#include <TMath.h>

/* $Id$ */

///////////////////////////////////////////////////////////////////
//                                                               //
// Implementation of the class used for SDD injector analysis    //
// Origin: F.Prino, Torino, prino@to.infn.it                     //
//                                                               //
///////////////////////////////////////////////////////////////////

ClassImp(AliITSOnlineSDDInjectors)

const Float_t AliITSOnlineSDDInjectors::fgkSaturation = 1008.;

//______________________________________________________________________
AliITSOnlineSDDInjectors::AliITSOnlineSDDInjectors():AliITSOnlineSDD(),fHisto(),fTbZero(0.),fParam(),fPolOrder(0),fMinDriftVel(0.),fMaxDriftVel(0.),fThreshold(0.),fTimeDiffTB()
{
  // default constructor
  SetMinDriftVel();
  SetMaxDriftVel();
  SetRangeLine1();
  SetRangeLine2();
  SetRangeLine3();
  SetPositions();
  SetPolOrder();
  SetThreshold();
  SetTimeDiffTB();
}
//______________________________________________________________________
AliITSOnlineSDDInjectors::AliITSOnlineSDDInjectors(Int_t mod, Int_t sid):AliITSOnlineSDD(mod,sid),fHisto(),fTbZero(0.),fParam(),fPolOrder(0),fMinDriftVel(0.),fMaxDriftVel(0.),fThreshold(0.),fTimeDiffTB()
{ 
// standard constructor
  SetMinDriftVel();
  SetMaxDriftVel();
  SetRangeLine1();
  SetRangeLine2();
  SetRangeLine3();
  SetPositions();
  SetPolOrder();
  SetThreshold();
  SetTimeDiffTB();
}
//______________________________________________________________________
AliITSOnlineSDDInjectors::~AliITSOnlineSDDInjectors(){
  // Destructor
  if(fHisto) delete fHisto;  
  if(fParam) delete [] fParam;
}
//______________________________________________________________________
void AliITSOnlineSDDInjectors::SetPositions(){
  // 
  Float_t kLinFromCenterUm[3]={31860.,17460.,660.};
  Float_t kAnodeFromCenterUm=35085;
  for(Int_t i=0;i<3;i++){
    fPosition[i]=kAnodeFromCenterUm-kLinFromCenterUm[i];
    fPosition[i]/=10000.; // from microns to cm
  }
}
//______________________________________________________________________
void AliITSOnlineSDDInjectors::Reset(){
  //
  for(Int_t i=0;i<kNInjectors;i++){ 
    fDriftVel[i]=0.;
    fSigmaDriftVel[i]=0.;
  }
  for(Int_t i=0;i<kNInjectors;i++){
    for(Int_t j=0;j<3;j++){
      fGoodInj[i][j]=0;
      fCentroid[i][j]=0.;
      fRMSCentroid[i][j]=0.;
    }
  }
}
//______________________________________________________________________
void AliITSOnlineSDDInjectors::AnalyzeEvent(TH2F* his){
  //
  Reset();
  fHisto=his;
  FindGoodInjectors();
  FindCentroids();
  CalcTimeBinZero();
  for(Int_t j=0;j<kNInjectors;j++) CalcDriftVelocity(j);
  FitDriftVelocityVsAnode();
}
//______________________________________________________________________
TGraphErrors* AliITSOnlineSDDInjectors::GetLineGraph(Int_t jlin) const{
  // 
  Float_t x[4],y[4],ex[4],ey[4];
  x[0]=0.;
  ex[0]=0.;
  y[0]=fTbZero;
  ey[0]=0.;
  for(Int_t i=0;i<3;i++){
    x[i+1]=fPosition[i];
    ex[i+1]=0.;
    y[i+1]=fCentroid[jlin][i];
    ey[i+1]=fRMSCentroid[jlin][i];
  }
  TGraphErrors *g=new TGraphErrors(4,x,y,ex,ey);
  return g;
}
//______________________________________________________________________
Float_t AliITSOnlineSDDInjectors::GetDriftCoordinate(Float_t cAnode, Float_t cTimeBin){
  //
  Float_t vel=0;
  for(Int_t i=0;i<=fPolOrder;i++) vel+=fParam[i]*TMath::Power(cAnode,(Float_t)i);
  return vel*(cTimeBin-(fTbZero-fTimeDiffTB))*25/1000.; 
}
//______________________________________________________________________
TGraphErrors* AliITSOnlineSDDInjectors::GetDriftVelocityGraph() const{
  // 
  Int_t ipt=0;
  TGraphErrors *g=new TGraphErrors(0);
  for(Int_t i=0;i<kNInjectors;i++){
    if(fDriftVel[i]>0){ 
      g->SetPoint(ipt,GetAnodeNumber(i),fDriftVel[i]);
      g->SetPointError(ipt,0,fSigmaDriftVel[i]);
      ipt++;
    }
  }
  return g;
}
//______________________________________________________________________
void AliITSOnlineSDDInjectors::CalcTimeBinZero(){
  //
  Float_t tzero=0.,intCont=0.;
  for(Int_t ian=0;ian<fgkNAnodes;ian++){
    for(Int_t itb=1;itb<fTbMin[0];itb++){
      Float_t cont=fHisto->GetBinContent(itb,ian+1);
      if(cont>fThreshold){
	tzero+=cont*float(itb);
	intCont+=cont;
      }
    }
  }
  if(intCont>0) fTbZero=tzero/intCont;
}
//______________________________________________________________________
void AliITSOnlineSDDInjectors::FitDriftVelocityVsAnode(){
  // fits the anode dependence of drift velocity with a polynomial
  const Int_t kNn=fPolOrder+1;
  Float_t **mat = new Float_t*[kNn];
  for(Int_t i=0; i < kNn; i++) mat[i] = new Float_t[kNn];
  Float_t *vect = new Float_t[kNn];
  for(Int_t k1=0;k1<kNn;k1++){
    vect[k1]=0;
    for(Int_t k2=0;k2<kNn;k2++){
      mat[k1][k2]=0;
      for(Int_t n=0; n<kNInjectors;n++){
	Float_t x=(Float_t)GetAnodeNumber(n);
	if(fDriftVel[n]>0) mat[k1][k2]+=TMath::Power(x,k1+k2)/TMath::Power(fSigmaDriftVel[n],2);
      }
    }
  }
  for(Int_t k1=0;k1<kNn;k1++){
    for(Int_t n=0; n<kNInjectors;n++){
      Float_t x=(Float_t)GetAnodeNumber(n);
      if(fDriftVel[n]>0) vect[k1]+=fDriftVel[n]*TMath::Power(x,k1)/TMath::Power(fSigmaDriftVel[n],2);
    }
  }
  Int_t *iPivot = new Int_t[kNn];
  Int_t *indxR = new Int_t[kNn];
  Int_t *indxC = new Int_t[kNn];
  for(Int_t i=0;i<kNn;i++) iPivot[i]=0;
  Int_t iCol=-1,iRow=-1;
  for(Int_t i=0;i<kNn;i++){
    Float_t big=0.;
    for(Int_t j=0;j<kNn;j++){
      if(iPivot[j]!=1){
	for(Int_t k=0;k<kNn;k++){
	   if(iPivot[k]==0){
	     if(TMath::Abs(mat[j][k])>=big){
	       big=TMath::Abs(mat[j][k]);
	       iRow=j;
	       iCol=k;
	     }
	  }
	}
      }
    }
    iPivot[iCol]++;
    Float_t aux;
    if(iRow!=iCol){
      for(Int_t l=0;l<kNn;l++){
	aux=mat[iRow][l];
	mat[iRow][l]=mat[iCol][l];
	mat[iCol][l]=aux;
      }
      aux=vect[iRow];
      vect[iRow]=vect[iCol];
      vect[iCol]=aux;
    }
    indxR[i]=iRow;
    indxC[i]=iCol;
    if(mat[iCol][iCol]==0) break;
    Float_t pivinv=1./mat[iCol][iCol];
    mat[iCol][iCol]=1;
    for(Int_t l=0;l<kNn;l++) mat[iCol][l]*=pivinv;
    vect[iCol]*=pivinv;
    for(Int_t m=0;m<kNn;m++){
      if(m!=iCol){
	aux=mat[m][iCol];
	mat[m][iCol]=0;
	for(Int_t n=0;n<kNn;n++) mat[m][n]-=mat[iCol][n]*aux;
	vect[m]-=vect[iCol]*aux;
      }
    }    
  }
  delete [] iPivot;
  delete [] indxR;
  delete [] indxC;

  if(fParam) delete [] fParam;
  fParam=new Float_t[kNn];
  for(Int_t i=0; i<kNn;i++)fParam[i]=vect[i];

  for(Int_t i=0; i < kNn; i++) delete [] mat[i];
  delete [] mat;
  delete [] vect;
}
//______________________________________________________________________
void AliITSOnlineSDDInjectors::CalcDriftVelocity(Int_t jlin){
  // 
  Float_t sumY=0,sumX=0,sumXX=0,sumYY=0.,sumXY=0,sumWEI=0.;
  Int_t npt=0;
  Float_t y[3],ey[3];
  Float_t tzero=0,erry=0;
  for(Int_t i=0;i<3;i++){ 
    y[i]=fCentroid[jlin][i];
    ey[i]=fRMSCentroid[jlin][i];
  }
  for(Int_t i=0;i<3;i++){
    if(fGoodInj[jlin][i]){
      sumY+=y[i]/ey[i]/ey[i];
      sumX+=fPosition[i]/ey[i]/ey[i];
      sumXX+=fPosition[i]*fPosition[i]/ey[i]/ey[i];
      sumYY+=y[i]*y[i]/ey[i]/ey[i];
      sumXY+=fPosition[i]*y[i]/ey[i]/ey[i];
      sumWEI+=1./ey[i]/ey[i];
      tzero=fTbZero/ey[i]/ey[i];
      erry=ey[i];
      npt++;
    }
  }
  Float_t vel=0,evel=0;
  if(npt>1){ 
    Float_t slope=(sumWEI*sumXY-sumY*sumX)/(sumWEI*sumXX-sumX*sumX);
    Float_t eslope=TMath::Sqrt(sumWEI/(sumWEI*sumXX-sumX*sumX));
    vel=1./slope*10000./25.;// micron/ns
    evel=eslope/slope/slope*10000./25.;// micron/ns
  }
  if(npt==1){
    Float_t slope=(sumY-tzero)/sumX;
    Float_t eslope=erry/sumX;
    vel=1./slope*10000./25.;// micron/ns    
    evel=eslope/slope/slope*10000./25.;// micron/ns
  }
  if(vel>fMaxDriftVel||vel<fMinDriftVel){ 
    vel=0.;
    evel=0.;
  }
  fDriftVel[jlin]=vel;
  fSigmaDriftVel[jlin]=evel;
}
//______________________________________________________________________
Int_t AliITSOnlineSDDInjectors::GetAnodeNumber(Int_t iInjLine) const{
  //
  Int_t ian=-1;
  if(iInjLine>32) return ian;
  if(!fSide){
    ian=iInjLine*8;
    if(iInjLine==32) ian--;
  }else{
    ian=iInjLine*8-1;
    if(iInjLine==0) ian=0;
  }
  return ian;
}

//______________________________________________________________________
void AliITSOnlineSDDInjectors::FindGoodInjectors(){
  // 
  for(Int_t iii=0;iii<kNInjectors;iii++){
    Int_t ian=GetAnodeNumber(iii);
    for(Int_t ninj=0;ninj<3;ninj++){
      for(Int_t jjj=fTbMin[ninj];jjj<fTbMax[ninj];jjj++){
	Float_t c1=fHisto->GetBinContent(jjj,ian+1);
	Float_t c2=fHisto->GetBinContent(jjj+1,ian+1);
	Float_t c3=fHisto->GetBinContent(jjj+2,ian+1);
	if(c1>fThreshold && c2>fThreshold && c3>fThreshold){ 
	  fGoodInj[iii][ninj]=1;
	  break;
	}
      }
      //      for(Int_t jjj=fTbMin[ninj];jjj<fTbMax[ninj];jjj++){
      //      	Float_t c1=fHisto->GetBinContent(jjj,ian+1);
      //      	if(c1>=fgkSaturation){
      //      	  fGoodInj[iii][ninj]=0;
      //	  break;
      //	}
      //      }
    }
  }
}
//______________________________________________________________________
void AliITSOnlineSDDInjectors::FindCentroids(){
  // 
  for(Int_t iii=0;iii<kNInjectors;iii++){
    Int_t ian=GetAnodeNumber(iii);
    for(Int_t ninj=0;ninj<3;ninj++){
      if(!fGoodInj[iii][ninj]) continue;
      Float_t maxcont=0;
      Int_t ilmax=-1;
      for(Int_t jjj=fTbMin[ninj];jjj<fTbMax[ninj];jjj++){
	Float_t cont=fHisto->GetBinContent(jjj,ian+1);
	if(cont>maxcont){
	  maxcont=cont;
	  ilmax=jjj;
	}
      }
      Float_t intCont=0;
      Int_t jjj=ilmax;
      while(1){
	Float_t cont=fHisto->GetBinContent(jjj,ian+1);
	if(cont<fThreshold) break;
	if(cont<fgkSaturation){
	  fCentroid[iii][ninj]+=cont*(Float_t)jjj;
	  fRMSCentroid[iii][ninj]+=cont*TMath::Power((Float_t)jjj,2);
	  intCont+=cont;
	}
	jjj--;
      }
      jjj=ilmax+1;
      while(1){
	Float_t cont=fHisto->GetBinContent(jjj,ian+1);
	if(cont<fThreshold) break;
	if(cont<fgkSaturation){
	  fCentroid[iii][ninj]+=cont*float(jjj);
	  fRMSCentroid[iii][ninj]+=cont*TMath::Power((Float_t)jjj,2);
	  intCont+=cont;
	}
	jjj++;
      }
      if(intCont>0){ 
	fCentroid[iii][ninj]/=intCont;
	fRMSCentroid[iii][ninj]=TMath::Sqrt(fRMSCentroid[iii][ninj]/intCont-fCentroid[iii][ninj]*fCentroid[iii][ninj]);
      }
      else{ 
	fCentroid[iii][ninj]=0.;
	fRMSCentroid[iii][ninj]=0.;
	fGoodInj[iii][ninj]=0;
      }
      if(fRMSCentroid[iii][ninj]==0) fGoodInj[iii][ninj]=0;
    }
  }
}
//______________________________________________________________________
void AliITSOnlineSDDInjectors::PrintInjMap(){
  //
  for(Int_t iii=0;iii<kNInjectors;iii++){
    printf("Line%d-Anode%d: %d %d %d\n",iii,GetAnodeNumber(iii),fGoodInj[iii][0],fGoodInj[iii][1],fGoodInj[iii][2]);
  }
}
//______________________________________________________________________
void AliITSOnlineSDDInjectors::PrintCentroids(){
  //
  for(Int_t iii=0;iii<kNInjectors;iii++){
    printf("Line%d-Anode%d: %f+-%f %f+-%f %f+-%f\n",iii,GetAnodeNumber(iii),fCentroid[iii][0],fRMSCentroid[iii][0],fCentroid[iii][1],fRMSCentroid[iii][1],fCentroid[iii][2],fRMSCentroid[iii][2]);
  }
}
//______________________________________________________________________
void AliITSOnlineSDDInjectors::WriteToASCII(Int_t evNumb, UInt_t timeStamp, Int_t optAppend){
  //
  Char_t outfilnam[100];
  sprintf(outfilnam,"SDDinj_mod%03d_sid%d.data",fModuleId,fSide);
  FILE* outf;
  if(optAppend==0) outf=fopen(outfilnam,"w");
  else outf=fopen(outfilnam,"a");
  fprintf(outf,"%d   %d   ",evNumb,timeStamp);
  for(Int_t ic=0;ic<fPolOrder+1;ic++){
    fprintf(outf,"%G ",fParam[ic]);
  }
  fprintf(outf,"\n");
  fclose(outf);  
}
