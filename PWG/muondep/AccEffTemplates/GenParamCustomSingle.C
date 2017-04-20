/*
 *  MuonGenerator.C
 *  aliroot_dev
 *
 *  Created by philippe pillot on 05/03/13.
 *  Copyright 2013 SUBATECH. All rights reserved.
 *
 */


#if !defined(__CINT__) || defined(__MAKECINT__)
#include <Riostream.h>
#include "TRandom.h"
#include "AliGenerator.h"
#include "AliGenParam.h"
#endif


static Int_t IpMuon( TRandom *ran);
static Double_t PtMuon( const Double_t *px, const Double_t */*dummy*/ );
static Double_t YMuon( const Double_t *py, const Double_t */*dummy*/ );
static Double_t V2Muon( const Double_t *pv, const Double_t */*dummy*/ );


//-------------------------------------------------------------------------
AliGenerator* GenParamCustomSingle()
{
  AliGenParam *singleMu = new AliGenParam(2,-1,PtMuon,YMuon,V2Muon,IpMuon);
  singleMu->SetMomentumRange(0., 1.e6);
  singleMu->SetPtRange(VAR_GENPARAMCUSTOMSINGLE_PTMIN, 999.);
  singleMu->SetYRange(-4.3, -2.3);
  singleMu->SetPhiRange(0., 360.);
  singleMu->SetForceDecay(kNoDecay);
  singleMu->SetTrackingFlag(1);
  return singleMu;
}

//-------------------------------------------------------------------------
Int_t IpMuon(TRandom *ran)
{
  // muon composition

  if (ran->Rndm() < 0.5 )
  {
    return 13;
  }
  else
  {
    return -13;
  }
}

//-------------------------------------------------------------------------
Double_t PtMuon( const Double_t *px, const Double_t */*dummy*/ )
{
  // muon pT

  Double_t x=*px;
  Double_t p0,p1,p2,p3,p4,p5;

  //Default param. = tuned MSL LHC15n
  p0 = VAR_GENPARAMCUSTOMSINGLE_PT_P0;//135.137
  p1 = VAR_GENPARAMCUSTOMSINGLE_PT_P1;//0.555323
  p2 = VAR_GENPARAMCUSTOMSINGLE_PT_P2;//0.578374
  p3 = VAR_GENPARAMCUSTOMSINGLE_PT_P3;//10.1345
  p4 = VAR_GENPARAMCUSTOMSINGLE_PT_P4;//0.000232233
  p5 = VAR_GENPARAMCUSTOMSINGLE_PT_P5;//-0.924726

  return p0 * (1. / TMath::Power(p1 + TMath::Power(x,p2), p3) + p4 * TMath::Exp(p5*x));
}

//-------------------------------------------------------------------------
Double_t YMuon( const Double_t *py, const Double_t */*dummy*/ )
{
  // muon y

  Double_t y = *py;
  Double_t p0,p1,p2,p3,p4;

  //Default param. = tuned MSL LHC15n
  p0 = VAR_GENPARAMCUSTOMSINGLE_Y_P0;//1.95551
  p1 = VAR_GENPARAMCUSTOMSINGLE_Y_P1;//0.
  p2 = VAR_GENPARAMCUSTOMSINGLE_Y_P2;//-0.104761
  p3 = VAR_GENPARAMCUSTOMSINGLE_Y_P3;//0.
  p4 = VAR_GENPARAMCUSTOMSINGLE_Y_P4;//0.00311324

  return p0 * (1. + p1*y + p2*y*y + p3*y*y*y + p4*y*y*y*y);
}

//-------------------------------------------------------------------------
Double_t V2Muon( const Double_t */*dummy*/, const Double_t */*dummy*/ )
{
  //muon v2
  return 0.;
}

