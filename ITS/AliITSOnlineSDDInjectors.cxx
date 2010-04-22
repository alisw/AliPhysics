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
#include "AliITSOnlineSDDInjectors.h"
#include "AliLog.h"
#include <TH2F.h>
#include <TF1.h>
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
const Int_t   AliITSOnlineSDDInjectors::fgkDefaultPolDegree = 3;
const Float_t AliITSOnlineSDDInjectors::fgkDefaultTimeStep = 50.;
const UShort_t AliITSOnlineSDDInjectors::fgkDefaultTbMin[kInjLines] = {10,50,100};
const UShort_t AliITSOnlineSDDInjectors::fgkDefaultTbMax[kInjLines] = {20,70,120};

//______________________________________________________________________
AliITSOnlineSDDInjectors::AliITSOnlineSDDInjectors():AliITSOnlineSDD(),fHisto(),fTbZero(0.),fRMSTbZero(0.),fNEvents(0),fParam(),fPolDegree(0),fActualPolDegree(0),fMinDriftSpeed(0.),fMaxDriftSpeed(0.),fMaxDriftSpeedErr(0.),fLowThreshold(0.),fHighThreshold(0.),fFirstPadForFit(0),fLastPadForFit(0),fPadStatusCutForFit(0),fTimeStep(0.),fUseTimeZeroSignal(kFALSE)
{
  // default constructor
  SetPositions();
  SetDefaults();
  SetTimeStep(fgkDefaultTimeStep);
}
//______________________________________________________________________
AliITSOnlineSDDInjectors::AliITSOnlineSDDInjectors(Int_t nddl, Int_t ncarlos, Int_t sid):AliITSOnlineSDD(nddl,ncarlos,sid),fHisto(),fTbZero(0.),fRMSTbZero(0.),fNEvents(0),fParam(),fPolDegree(0),fActualPolDegree(0),fMinDriftSpeed(0.),fMaxDriftSpeed(0.),fMaxDriftSpeedErr(0.),fLowThreshold(0.),fHighThreshold(0.),fFirstPadForFit(0),fLastPadForFit(0),fPadStatusCutForFit(0),fTimeStep(0.),fUseTimeZeroSignal(kFALSE)
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
  for(Int_t i=0;i<kInjLines;i++) {
    SetInjLineRange(i,fgkDefaultTbMin[i],fgkDefaultTbMax[i]);
    SetUseLine(i,kTRUE);
  }
  SetThresholds(fgkDefaultLThreshold,fgkDefaultHThreshold);
  SetPolDegree(fgkDefaultPolDegree);
  SetMinDriftSpeed(fgkDefaultMinSpeed);
  SetMaxDriftSpeed(fgkDefaultMaxSpeed);
  SetMaxDriftSpeedErr(fgkDefaultMaxErr);
  SetFitLimits(1,kInjPads-2); // exclude first and last pad
  SetPadStatusCutForFit();
}
//______________________________________________________________________
void AliITSOnlineSDDInjectors::SetPositions(){
  // 
  Double_t xLinFromCenterUm[kInjLines]={31860.,17460.,660.};
  Double_t xAnodeFromCenterUm=35085;
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
  AddEvent(his);
  FitDriftSpeedVsAnode();
}
//______________________________________________________________________
void AliITSOnlineSDDInjectors::AddEvent(TH2F* his){
  // Add the drift speed from current event to the average value
  if(fNEvents==0){
    for(Int_t i=0;i<kInjPads;i++){ 
      fSumDriftSpeed[i]=0.;
      fSumSqDriftSpeed[i]=0.;
      fSumPadStatus[i]=0;
      fSumPadStatusCut[i]=0;
      fNEventsInPad[i]=0;
    }
  }
  Reset();
  fHisto=his;
  FindGoodInjectors();
  FindCentroids();
  CalcTimeBinZero();
  for(Int_t j=0;j<kInjPads;j++){ 
    CalcDriftSpeed(j);
    Int_t padStatus=GetInjPadStatus(j);
    fSumPadStatus[j]+=padStatus;
    if(padStatus>fPadStatusCutForFit){
      fSumDriftSpeed[j]+=fDriftSpeed[j];
      fSumSqDriftSpeed[j]+=fDriftSpeed[j]*fDriftSpeed[j];
      fSumPadStatusCut[j]+=padStatus;
      fNEventsInPad[j]++;
    }
  }
  ++fNEvents;
}
//______________________________________________________________________
Double_t AliITSOnlineSDDInjectors::GetRMSDriftSpeed(Int_t ipad) const {
  // 
  if(fNEventsInPad[ipad]<=1) return 0.;
  Double_t mean=fSumDriftSpeed[ipad]/(Double_t)fNEventsInPad[ipad];
  Double_t diff=fSumSqDriftSpeed[ipad]/(Double_t)fNEventsInPad[ipad]-mean*mean;
  if(diff<0.) diff=0.;
  return TMath::Sqrt(diff);
}

//______________________________________________________________________
void AliITSOnlineSDDInjectors::FitMeanDriftSpeedVsAnode(){
  // Calculates
  if(fNEvents==0) return;
  for(Int_t i=0;i<kInjPads;i++){ 
    fDriftSpeed[i]=GetMeanDriftSpeed(i);
    Int_t padStatusCut=(Int_t)(GetMeanPadStatusCut(i)+0.5);
    for(Int_t ilin=0; ilin<kInjLines ; ilin++) fGoodInj[i][ilin]=(padStatusCut&1<<ilin)>>ilin;
    if(fNEventsInPad[i]>1){
      Double_t rms=GetRMSDriftSpeed(i);
      if(rms>0.) fDriftSpeedErr[i]=rms/TMath::Sqrt(fNEventsInPad[i]);
    }else{
      for(Int_t ilin=0; ilin<kInjLines ; ilin++) fGoodInj[i][ilin]=0;
    }
  }
  FitDriftSpeedVsAnode();
  for(Int_t i=0;i<kInjPads;i++){ 
    Int_t padStatus=(Int_t)(GetMeanPadStatusCut(i)+0.5);
    for(Int_t ilin=0; ilin<kInjLines ; ilin++) fGoodInj[i][ilin]=(padStatus&1<<ilin)>>ilin;
  }
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
  // Get time zero from trigger signal
  Double_t tzero=0.,intCont=0.,rmsPeak=0.;
  Bool_t isTbUsed[256];
  Int_t nTbUsed=0;
  for(Int_t i=0;i<256;i++) isTbUsed[i]=0;
  for(Int_t ian=0;ian<fgkNAnodes;ian++){
    for(Int_t itb=1;itb<fTbMin[0];itb++){
      Double_t cont=fHisto->GetBinContent(itb,ian+1);
      Double_t contm1=fHisto->GetBinContent(itb+1,ian+1);
      Double_t contp1=fHisto->GetBinContent(itb-1,ian+1);
      if(cont>fLowThreshold){
	if(contm1>fHighThreshold || cont>fHighThreshold || contp1>fHighThreshold){
	  tzero+=cont*float(itb);
	  rmsPeak+=cont*float(itb)*float(itb);
	  intCont+=cont;
	  if(!isTbUsed[itb]){
	    isTbUsed[itb]=1;
	    ++nTbUsed;
	  }
	}
      }
    }
  }
  if(intCont>0){ 
    fTbZero=tzero/intCont;
    fRMSTbZero=TMath::Sqrt(rmsPeak/intCont-fTbZero*fTbZero);
  }
  if(nTbUsed==1) fRMSTbZero=0.5; 
}

//______________________________________________________________________
void AliITSOnlineSDDInjectors::FitDriftSpeedVsAnode(){
  // fits the anode dependence of drift speed

  Float_t rangeForMax[2]={78.,178.};
  PolyFit(fPolDegree);
  fActualPolDegree=fPolDegree;
  if(fPolDegree==3){
    Double_t deltasq=fParam[2]*fParam[2]-3*fParam[1]*fParam[3];
    Double_t zero1=-999.;
    Double_t zero2=-999.;
    if(deltasq>=0. && TMath::Abs(fParam[3])>0.){
      Double_t delta=TMath::Sqrt(deltasq);
      zero1=(-fParam[2]+delta)/3./fParam[3];
      zero2=(-fParam[2]-delta)/3./fParam[3];
    }
    Bool_t twoZeroes=kFALSE;
    Bool_t oneZero=kFALSE;
    if(zero1>0. && zero1<256. && zero2>0. && zero2<256.) twoZeroes=kTRUE;
    if(zero1>rangeForMax[0] && zero1<rangeForMax[1]) oneZero=kTRUE;
    if(zero2>rangeForMax[0] && zero2<rangeForMax[1]) oneZero=kTRUE;
    if(!oneZero || twoZeroes){
      PolyFit(2);
      Double_t xmax=-999.;
      if(fParam[2]<0.) xmax=-fParam[1]/2./fParam[2];
      if(xmax>rangeForMax[0] && xmax<rangeForMax[1]){
	fActualPolDegree=2;
      }else{
	Double_t averSpeed=0.;
	Double_t sumWei=0.;
	Int_t nUsedPts=0;
	for(Int_t jpad=fFirstPadForFit; jpad<=fLastPadForFit; jpad++){
	  if(fDriftSpeed[jpad]>0 && GetInjPadStatus(jpad)>fPadStatusCutForFit){
	    Double_t wei=1./fDriftSpeedErr[jpad]/fDriftSpeedErr[jpad];
	    averSpeed+=wei*fDriftSpeed[jpad];
	    sumWei+=wei;
	    nUsedPts++;
	  }
	}
	if(sumWei>0.) averSpeed/=sumWei;
	if(nUsedPts<fPolDegree+1) averSpeed=0;
	fParam[0]=averSpeed;
	for(Int_t i=1; i < fPolDegree+1; i++) fParam[i]=0.;
	fActualPolDegree=0;
      }
    }
  }
}
//______________________________________________________________________
void AliITSOnlineSDDInjectors::PolyFit(Int_t degree){
  // fits the anode dependence of drift speed with a polynomial function
  const Int_t kNn=degree+1;
  const Int_t kDimens=fPolDegree+1;

  Double_t **mat = new Double_t*[kNn];
  for(Int_t i=0; i < kNn; i++) mat[i] = new Double_t[kNn];
  Double_t *vect = new Double_t[kNn];

  for(Int_t k1=0;k1<kNn;k1++){
    vect[k1]=0;
    for(Int_t k2=0;k2<kNn;k2++){
      mat[k1][k2]=0;
    }
  }
  Int_t npts = 0;
  for(Int_t k1=0;k1<kNn;k1++){
    for(Int_t jpad=fFirstPadForFit; jpad<=fLastPadForFit; jpad++){
      Double_t x=(Double_t)GetAnodeNumber(jpad);
      if(fDriftSpeed[jpad]>0 && GetInjPadStatus(jpad)>fPadStatusCutForFit){
	vect[k1]+=fDriftSpeed[jpad]*TMath::Power(x,k1)/TMath::Power(fDriftSpeedErr[jpad],2);	
	if(k1==0) npts++;
	for(Int_t k2=0;k2<kNn;k2++){
	  mat[k1][k2]+=TMath::Power(x,k1+k2)/TMath::Power(fDriftSpeedErr[jpad],2);
	}
      }
    }
  }
  if(npts<fPolDegree+1){ 
    if(fParam) delete [] fParam;
    fParam=new Double_t[kDimens];
    for(Int_t i=0; i<kDimens;i++)fParam[i]=0;
  }else{
    Int_t *iPivot = new Int_t[kNn];
    Int_t *indxR = new Int_t[kNn];
    Int_t *indxC = new Int_t[kNn];
    for(Int_t i=0;i<kNn;i++) iPivot[i]=0;
    Int_t iCol=-1,iRow=-1;
    for(Int_t i=0;i<kNn;i++){
      Double_t big=0.;
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
      Double_t aux;
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
      Double_t pivinv=1./mat[iCol][iCol];
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
    fParam=new Double_t[kDimens];
    for(Int_t i=0; i<kNn;i++)fParam[i]=vect[i];
    if(degree<fPolDegree) for(Int_t i=kNn; i<kDimens;i++)fParam[i]=0.;
  }

  for(Int_t i=0; i < kNn; i++) delete [] mat[i];
  delete [] mat;
  delete [] vect;
}
//______________________________________________________________________
void AliITSOnlineSDDInjectors::CalcDriftSpeed(Int_t jpad){
  // 
  Double_t sumY=0,sumX=0,sumXX=0,sumYY=0.,sumXY=0,sumWEI=0.;
  Int_t npt=0;
  Double_t y[kInjLines],ey[kInjLines];
  Double_t tzero=0,erry=0;
  for(Int_t i=0;i<kInjLines;i++){ 
    y[i]=fCentroid[jpad][i];
    ey[i]=fRMSCentroid[jpad][i];
  }
  for(Int_t i=0;i<kInjLines;i++){
    if(!fUseLine[i]) continue;
    if(fGoodInj[jpad][i] && ey[i]!=0){
      sumY+=y[i]/ey[i]/ey[i];
      sumX+=fPosition[i]/ey[i]/ey[i];
      sumXX+=fPosition[i]*fPosition[i]/ey[i]/ey[i];
      sumYY+=y[i]*y[i]/ey[i]/ey[i];
      sumXY+=fPosition[i]*y[i]/ey[i]/ey[i];
      sumWEI+=1./ey[i]/ey[i];
      tzero=fTbZero/ey[i]/ey[i];
      erry=ey[i]/ey[i]/ey[i];
      npt++;
    }
  }
  Double_t slope=0.,eslope=0.;
  if(npt==1){
    slope=(sumY-tzero)/sumX;
    eslope=erry/sumX;
  }
  if(npt>1){ 
    if(fUseTimeZeroSignal){
      sumY+=fTbZero/fRMSTbZero/fRMSTbZero;
      sumX+=0.;
      sumXX+=0.;
      sumYY+=fTbZero*fTbZero/fRMSTbZero/fRMSTbZero;
      sumXY+=0.;
      sumWEI+=1./fRMSTbZero/fRMSTbZero;
    }
    slope=(sumWEI*sumXY-sumY*sumX)/(sumWEI*sumXX-sumX*sumX);
    eslope=TMath::Sqrt(sumWEI/(sumWEI*sumXX-sumX*sumX));
  }

  Double_t vel=0,evel=0;
  if(slope!=0. && fTimeStep>0.){
    vel=1./slope*10000./fTimeStep;// micron/ns
    evel=eslope/slope/slope*10000./fTimeStep;// micron/ns
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
      Double_t maxcont=0;
      Int_t ilmax=-1;
      for(Int_t jjj=fTbMin[jlin];jjj<fTbMax[jlin];jjj++){
	Double_t cont=fHisto->GetBinContent(jjj,ian+1);
	if(cont>maxcont){
	  maxcont=cont;
	  ilmax=jjj;
	}
      }
      Double_t intCont=0;
      Int_t jjj=ilmax;
      while(1){
	Double_t cont=fHisto->GetBinContent(jjj,ian+1);
	if(cont<fLowThreshold) break;
	if(cont<fgkSaturation){
	  fCentroid[jpad][jlin]+=cont*(Double_t)jjj;
	  fRMSCentroid[jpad][jlin]+=cont*(Double_t)jjj*(Double_t)jjj;
	  intCont+=cont;
	}
	jjj--;
      }
      jjj=ilmax+1;
      while(1){
	Double_t cont=fHisto->GetBinContent(jjj,ian+1);
	if(cont<fLowThreshold) break;
	if(cont<fgkSaturation){
	  fCentroid[jpad][jlin]+=cont*float(jjj);
	  fRMSCentroid[jpad][jlin]+=cont*(Double_t)jjj*(Double_t)jjj;
	  intCont+=cont;
	}
	jjj++;
      }
      if(intCont>0){ 
	fCentroid[jpad][jlin]/=intCont;
	fRMSCentroid[jpad][jlin]=TMath::Sqrt(fRMSCentroid[jpad][jlin]/intCont-fCentroid[jpad][jlin]*fCentroid[jpad][jlin])/TMath::Sqrt(intCont);
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
    fprintf(outf,"%d\n",fActualPolDegree);
  }
  else outf=fopen(outfilnam,"a");
  fprintf(outf,"%d   %d   ",evNumb,timeStamp);
  for(Int_t ic=0;ic<fPolDegree+1;ic++){
    fprintf(outf,"%G ",fParam[ic]);
  }
  fprintf(outf,"\n");
  fclose(outf);  
}
//______________________________________________________________________
TH1F* AliITSOnlineSDDInjectors::GetMeanDriftSpeedVsPadHisto() const{
  Char_t hisnam[20];
  sprintf(hisnam,"hdrsp%02dc%02ds%d",fDDL,fCarlos,fSide);
  TH1F* h=new TH1F(hisnam,"",kInjPads,-0.5,kInjPads-0.5);
  if(fNEvents>0){
    for(Int_t i=0;i<kInjPads;i++){ 
      h->SetBinContent(i+1,GetMeanDriftSpeed(i));    
      Double_t rms=GetRMSDriftSpeed(i);
      Double_t err=0.;
      if(rms>0.) err=rms/TMath::Sqrt(fNEventsInPad[i]);
      h->SetBinError(i+1,err);
    }
  }
  return h;
}
//______________________________________________________________________
Bool_t AliITSOnlineSDDInjectors::WriteToROOT(TFile *fil) const {
  //
  if(fil==0){ 
    AliWarning("Invalid pointer to ROOT file");
    return kFALSE;    
  }  
  Char_t hisnam[20];
  fil->cd();
  sprintf(hisnam,"hdrsp%02dc%02ds%d",fDDL,fCarlos,fSide);
  TH1F hdsp(hisnam,"",kInjPads,-0.5,kInjPads-0.5);
  if(fNEvents==0){
    AliWarning("Zero analyzed events");
    return kFALSE;    
  }  
    
  for(Int_t i=0;i<kInjPads;i++){ 
    hdsp.SetBinContent(i+1,GetMeanDriftSpeed(i));    
    Double_t rms=GetRMSDriftSpeed(i);
    Double_t err=0.;
    if(rms>0.) err=rms/TMath::Sqrt(fNEventsInPad[i]);
    hdsp.SetBinError(i+1,err);
  }
  hdsp.Write();
  return kTRUE;    
}
//______________________________________________________________________
void AliITSOnlineSDDInjectors::WriteInjectorStatusToASCII(){
  // dump status of injectors encoded into UInt_t
  // 5 bits (value 0-31) to store number of pads with given status
  Char_t outfilnam[100];
  sprintf(outfilnam,"SDDinj_ddl%02dc%02d_sid%d.data",fDDL,fCarlos,fSide);  
  FILE* outf=fopen(outfilnam,"a");
  Int_t n[8]={0,0,0,0,0,0,0,0};
  for(Int_t jpad=fFirstPadForFit; jpad<=fLastPadForFit; jpad++){
    Int_t statusPad=GetInjPadStatus(jpad);
    ++n[statusPad];
  }
  UInt_t statusInj=0;
  statusInj+=(n[7]&0x1F)<<25; // bits 25-29: n. of pads with status 7
  statusInj+=(n[6]&0x1F)<<20; // bits 20-24: n. of pads with status 6
  statusInj+=(n[5]&0x1F)<<15; // bits 15-19: n. of pads with status 5
  statusInj+=(n[4]&0x1F)<<10; // bits 10-14: n. of pads with status 4
  statusInj+=(n[3]&0x1F)<<5;  // bits  5- 9: n. of pads with status 3
  statusInj+=(n[2]&0x1F);     // bits  0- 4: n. of pads with status 2

  fprintf(outf,"-99 %u\n",statusInj); // -99 used in preprocessor to find line
                                      // with injector status info
  fclose(outf);  
  
}
