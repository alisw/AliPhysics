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
const Float_t AliITSOnlineSDDInjectors::fgkDefaultLThreshold = 5.;
const Float_t AliITSOnlineSDDInjectors::fgkDefaultHThreshold = 25.;
const Float_t AliITSOnlineSDDInjectors::fgkDefaultMinSpeed = 5.5;
const Float_t AliITSOnlineSDDInjectors::fgkDefaultMaxSpeed = 9.0;
const Float_t AliITSOnlineSDDInjectors::fgkDefaultMaxErr = 1.5;
const Int_t   AliITSOnlineSDDInjectors::fgkDefaultPolOrder = 3;
const Float_t AliITSOnlineSDDInjectors::fgkDefaultTimeStep = 50.;
const UShort_t AliITSOnlineSDDInjectors::fgkDefaultTbMin[kInjLines] = {10,50,100};
const UShort_t AliITSOnlineSDDInjectors::fgkDefaultTbMax[kInjLines] = {20,70,120};

//______________________________________________________________________
AliITSOnlineSDDInjectors::AliITSOnlineSDDInjectors():AliITSOnlineSDD(),fHisto(),fTbZero(0.),fParam(),fPolOrder(0),fMinDriftSpeed(0.),fMaxDriftSpeed(0.),fMaxDriftSpeedErr(0.),fLowThreshold(0.),fHighThreshold(0.),fFirstPadForFit(0),fLastPadForFit(0),fPadStatusCutForFit(0),fTimeStep(0.)
{
  // default constructor
  SetPositions();
  SetDefaults();
  SetTimeStep(fgkDefaultTimeStep);
}
//______________________________________________________________________
AliITSOnlineSDDInjectors::AliITSOnlineSDDInjectors(Int_t nddl, Int_t ncarlos, Int_t sid):AliITSOnlineSDD(nddl,ncarlos,sid),fHisto(),fTbZero(0.),fParam(),fPolOrder(0),fMinDriftSpeed(0.),fMaxDriftSpeed(0.),fMaxDriftSpeedErr(0.),fLowThreshold(0.),fHighThreshold(0.),fFirstPadForFit(0),fLastPadForFit(0),fPadStatusCutForFit(0),fTimeStep(0.)
{ 
// standard constructor
  SetPositions();
  SetDefaults();
  SetTimeStep(fgkDefaultTimeStep);
}
//______________________________________________________________________
AliITSOnlineSDDInjectors::~AliITSOnlineSDDInjectors(){
  // Destructor
  // fHisto should not be deleted here because it points to an histo created 
  // by the external code which calls the method AnalyzeEvent
  // if(fHisto) delete fHisto;  
  if(fParam) delete [] fParam;
}
//______________________________________________________________________
void AliITSOnlineSDDInjectors::SetDefaults(){
  for(Int_t i=0;i<kInjLines;i++) 
    SetInjLineRange(i,fgkDefaultTbMin[i],fgkDefaultTbMax[i]);
  SetThresholds(fgkDefaultLThreshold,fgkDefaultHThreshold);
  SetPolOrder(fgkDefaultPolOrder);
  SetMinDriftSpeed(fgkDefaultMinSpeed);
  SetMaxDriftSpeed(fgkDefaultMaxSpeed);
  SetMaxDriftSpeedErr(fgkDefaultMaxErr);
  SetFitLimits(1,kInjPads-2); // exclude first and last pad
  SetPadStatusCutForFit();
}
//______________________________________________________________________
void AliITSOnlineSDDInjectors::SetPositions(){
  // 
  Float_t xLinFromCenterUm[kInjLines]={31860.,17460.,660.};
  Float_t xAnodeFromCenterUm=35085;
  for(Int_t i=0;i<kInjLines;i++){
    fPosition[i]=xAnodeFromCenterUm-xLinFromCenterUm[i];
    fPosition[i]/=10000.; // from microns to cm
  }
}
//______________________________________________________________________
void AliITSOnlineSDDInjectors::Reset(){
  //
  for(Int_t i=0;i<kInjPads;i++){ 
    fDriftSpeed[i]=0.;
    fDriftSpeedErr[i]=0.;
  }
  for(Int_t i=0;i<kInjPads;i++){
    for(Int_t j=0;j<kInjLines;j++){
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
  for(Int_t j=0;j<kInjPads;j++) CalcDriftSpeed(j);
  FitDriftSpeedVsAnode();
}
//______________________________________________________________________
TGraphErrors* AliITSOnlineSDDInjectors::GetTimeVsDistGraph(Int_t jpad) const{
  // 
  const Int_t kPts=kInjLines+1;
  Float_t x[kPts],y[kPts],ex[kPts],ey[kPts];
  x[0]=0.;
  ex[0]=0.;
  y[0]=fTbZero;
  ey[0]=0.;
  for(Int_t i=0;i<kInjLines;i++){
    x[i+1]=fPosition[i];
    ex[i+1]=0.;
    y[i+1]=fCentroid[jpad][i];
    ey[i+1]=fRMSCentroid[jpad][i];
  }
  TGraphErrors *g=new TGraphErrors(4,x,y,ex,ey);
  return g;
}

//______________________________________________________________________
TGraphErrors* AliITSOnlineSDDInjectors::GetDriftSpeedGraph() const{
  // 
  Int_t ipt=0;
  TGraphErrors *g=new TGraphErrors(0);
  for(Int_t i=0;i<kInjPads;i++){
    if(fDriftSpeed[i]>0){ 
      g->SetPoint(ipt,GetAnodeNumber(i),fDriftSpeed[i]);
      g->SetPointError(ipt,0,fDriftSpeedErr[i]);
      ipt++;
    }
  }
  return g;
}
//______________________________________________________________________
TGraphErrors* AliITSOnlineSDDInjectors::GetSelectedDriftSpeedGraph(Int_t minAcceptStatus) const{
  // TGraphErrors with only pads with status of injector >= minAcceptStatus
  Int_t ipt=0;
  TGraphErrors *g=new TGraphErrors(0);
  for(Int_t i=0;i<kInjPads;i++){
    Int_t padStatus = GetInjPadStatus(i);
    if(fDriftSpeed[i]>0 && padStatus >= minAcceptStatus ){
      g->SetPoint(ipt,GetAnodeNumber(i),fDriftSpeed[i]);
      g->SetPointError(ipt,0,fDriftSpeedErr[i]);
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
      if(cont>fLowThreshold){
	tzero+=cont*float(itb);
	intCont+=cont;
      }
    }
  }
  if(intCont>0) fTbZero=tzero/intCont;
}
//______________________________________________________________________
void AliITSOnlineSDDInjectors::FitDriftSpeedVsAnode(){
  // fits the anode dependence of drift speed with a polynomial function
  const Int_t kNn=fPolOrder+1;
  Float_t **mat = new Float_t*[kNn];
  for(Int_t i=0; i < kNn; i++) mat[i] = new Float_t[kNn];
  Float_t *vect = new Float_t[kNn];

  for(Int_t k1=0;k1<kNn;k1++){
    vect[k1]=0;
    for(Int_t k2=0;k2<kNn;k2++){
      mat[k1][k2]=0;
    }
  }
  Int_t npts = 0;
  for(Int_t k1=0;k1<kNn;k1++){
    for(Int_t jpad=fFirstPadForFit; jpad<=fLastPadForFit; jpad++){
      Float_t x=(Float_t)GetAnodeNumber(jpad);
      if(fDriftSpeed[jpad]>0 && GetInjPadStatus(jpad)>fPadStatusCutForFit){
	  vect[k1]+=fDriftSpeed[jpad]*TMath::Power(x,k1)/TMath::Power(fDriftSpeedErr[jpad],2);	
	  if(k1==0) npts++;
	  for(Int_t k2=0;k2<kNn;k2++){
	    mat[k1][k2]+=TMath::Power(x,k1+k2)/TMath::Power(fDriftSpeedErr[jpad],2);
	  }
      }
    }
  }
  if(npts<fPolOrder+1){ 
    if(fParam) delete [] fParam;
    fParam=new Float_t[kNn];
    for(Int_t i=0; i<kNn;i++)fParam[i]=0;
  }else{
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
  }

  for(Int_t i=0; i < kNn; i++) delete [] mat[i];
  delete [] mat;
  delete [] vect;
}
//______________________________________________________________________
void AliITSOnlineSDDInjectors::CalcDriftSpeed(Int_t jpad){
  // 
  Float_t sumY=0,sumX=0,sumXX=0,sumYY=0.,sumXY=0,sumWEI=0.;
  Int_t npt=0;
  Float_t y[kInjLines],ey[kInjLines];
  Float_t tzero=0,erry=0;
  for(Int_t i=0;i<kInjLines;i++){ 
    y[i]=fCentroid[jpad][i];
    ey[i]=fRMSCentroid[jpad][i];
  }
  for(Int_t i=0;i<kInjLines;i++){
    if(fGoodInj[jpad][i] && ey[i]!=0){
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
    if(slope!=0 && fTimeStep>0.){
      vel=1./slope*10000./fTimeStep;// micron/ns
      evel=eslope/slope/slope*10000./fTimeStep;// micron/ns
    }
  }
  if(npt==1){
    Float_t slope=(sumY-tzero)/sumX;
    Float_t eslope=erry/sumX;
    if(slope!=0 && fTimeStep>0.){
      vel=1./slope*10000./fTimeStep;// micron/ns    
      evel=eslope/slope/slope*10000./fTimeStep;// micron/ns
    }
  }
  if(vel>fMaxDriftSpeed||vel<fMinDriftSpeed || evel>fMaxDriftSpeedErr){ 
    vel=0.;
    evel=0.;
  }
  fDriftSpeed[jpad]=vel;
  fDriftSpeedErr[jpad]=evel;
}
//______________________________________________________________________
Int_t AliITSOnlineSDDInjectors::GetAnodeNumber(Int_t iInjPad) const{
  // Injectors location along anodes:
  // Side left  (UP)   - channel 0: injectors on anodes 0,7,15,...,247,255 
  // Side right (DOWN) - channel 1: injectors on anodes 0,8,16,...,248,255 
  Int_t ian=-1;
  if(iInjPad>=kInjPads) return ian;
  if(fSide==1){  // right side
    ian=iInjPad*8;
    if(iInjPad==32) ian--;
  }else{         // left side
    ian=iInjPad*8-1;
    if(iInjPad==0) ian=0;
  }
  return ian;
}
//______________________________________________________________________
Int_t AliITSOnlineSDDInjectors::GetInjPadNumberFromAnode(Int_t nAnode) const{
  //
  Int_t iInjPad=-1;
  if(fSide==1){  // right side
    if(nAnode%8==0) iInjPad=nAnode/8;
    if(nAnode==255) iInjPad=32;
  }else{         // left side
    if(nAnode%8==7) iInjPad=1+nAnode/8;
    if(nAnode==0) iInjPad=0;
  }
  if(nAnode>=256) iInjPad=-1;
  return iInjPad;
}
//______________________________________________________________________
Int_t AliITSOnlineSDDInjectors::GetInjPadStatus(Int_t jpad) const{
  // returns an integer value with status of injector lines for given pad/anode
  // status=7  -->  111  all injector are good
  // status=6  -->  110  1st line (close to anodes) is bad, other two are good
  // ....
  // status=1  -->  001  only 1st line (close to anodes) good
  // status=0  -->  000  all lines are bad
  Int_t istatus=0;
  if(jpad>=0 && jpad<kInjPads){
    for(Int_t jlin=0;jlin<kInjLines;jlin++) istatus+=fGoodInj[jpad][jlin]<<jlin;
  }
  return istatus;
}
//______________________________________________________________________
void AliITSOnlineSDDInjectors::FindGoodInjectors(){
  // 
  for(Int_t jpad=0;jpad<kInjPads;jpad++){
    Int_t ian=GetAnodeNumber(jpad);
    for(Int_t jlin=0;jlin<kInjLines;jlin++){
      for(Int_t jjj=fTbMin[jlin];jjj<fTbMax[jlin];jjj++){
	Float_t c1=fHisto->GetBinContent(jjj,ian+1);
	Float_t c2=fHisto->GetBinContent(jjj+1,ian+1);
	//	Float_t c3=fHisto->GetBinContent(jjj+2,ian+1);
	if(c1>fLowThreshold && c2>fLowThreshold){ 
	  if(c1>fHighThreshold || c2>fHighThreshold){
	    fGoodInj[jpad][jlin]=1;
	    break;
	  }
	}
      }
    }
  }
}
//______________________________________________________________________
void AliITSOnlineSDDInjectors::FindCentroids(){
  // 
  for(Int_t jpad=0;jpad<kInjPads;jpad++){
    Int_t ian=GetAnodeNumber(jpad);
    for(Int_t jlin=0;jlin<kInjLines;jlin++){
      if(!fGoodInj[jpad][jlin]) continue;
      Float_t maxcont=0;
      Int_t ilmax=-1;
      for(Int_t jjj=fTbMin[jlin];jjj<fTbMax[jlin];jjj++){
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
	if(cont<fLowThreshold) break;
	if(cont<fgkSaturation){
	  fCentroid[jpad][jlin]+=cont*(Float_t)jjj;
	  fRMSCentroid[jpad][jlin]+=cont*TMath::Power((Float_t)jjj,2);
	  intCont+=cont;
	}
	jjj--;
      }
      jjj=ilmax+1;
      while(1){
	Float_t cont=fHisto->GetBinContent(jjj,ian+1);
	if(cont<fLowThreshold) break;
	if(cont<fgkSaturation){
	  fCentroid[jpad][jlin]+=cont*float(jjj);
	  fRMSCentroid[jpad][jlin]+=cont*TMath::Power((Float_t)jjj,2);
	  intCont+=cont;
	}
	jjj++;
      }
      if(intCont>0){ 
	fCentroid[jpad][jlin]/=intCont;
	fRMSCentroid[jpad][jlin]=TMath::Sqrt(fRMSCentroid[jpad][jlin]/intCont-fCentroid[jpad][jlin]*fCentroid[jpad][jlin]);
      }
      else{ 
	fCentroid[jpad][jlin]=0.;
	fRMSCentroid[jpad][jlin]=0.;
	fGoodInj[jpad][jlin]=0;
      }
      if(fRMSCentroid[jpad][jlin]==0) fGoodInj[jpad][jlin]=0;
    }
  }
}
//______________________________________________________________________
void AliITSOnlineSDDInjectors::PrintInjectorStatus(){
  //
  for(Int_t jpad=0;jpad<kInjPads;jpad++){
    printf("Line%d-Anode%d: %d %d %d\n",jpad,GetAnodeNumber(jpad),fGoodInj[jpad][0],fGoodInj[jpad][1],fGoodInj[jpad][2]);
  }
}
//______________________________________________________________________
void AliITSOnlineSDDInjectors::PrintCentroids(){
  //
  for(Int_t jpad=0;jpad<kInjPads;jpad++){
    printf("Line%d-Anode%d: %f+-%f %f+-%f %f+-%f\n",jpad,GetAnodeNumber(jpad),fCentroid[jpad][0],fRMSCentroid[jpad][0],fCentroid[jpad][1],fRMSCentroid[jpad][1],fCentroid[jpad][2],fRMSCentroid[jpad][2]);
  }
}
//______________________________________________________________________
void AliITSOnlineSDDInjectors::WriteToASCII(Int_t evNumb, UInt_t timeStamp, Int_t optAppend){
  //
  Char_t outfilnam[100];
  sprintf(outfilnam,"SDDinj_ddl%02dc%02d_sid%d.data",fDDL,fCarlos,fSide);  
  FILE* outf;
  if(optAppend==0){ 
    outf=fopen(outfilnam,"w");
    fprintf(outf,"%d\n",fPolOrder);
  }
  else outf=fopen(outfilnam,"a");
  fprintf(outf,"%d   %d   ",evNumb,timeStamp);
  for(Int_t ic=0;ic<fPolOrder+1;ic++){
    fprintf(outf,"%G ",fParam[ic]);
  }
  fprintf(outf,"\n");
  fclose(outf);  
}
