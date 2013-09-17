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
  AliGenParam *singleMu = new AliGenParam(1,-1,PtMuon,YMuon,V2Muon,IpMuon);
  singleMu->SetMomentumRange(0,1e6);
  singleMu->SetPtRange(VAR_GENPARAMCUSTOMSINGLE_PTMIN,999.);
  singleMu->SetYRange(-4.2, -2.3);
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
  Float_t p0,p1,p2,p3;
  p0 = VAR_GENPARAMCUSTOMSINGLE_PT_P0; //4.05962;
  p1 = VAR_GENPARAMCUSTOMSINGLE_PT_P1; //1;
  p2 = VAR_GENPARAMCUSTOMSINGLE_PT_P2; //2.46187;
  p3 = VAR_GENPARAMCUSTOMSINGLE_PT_P3; //2.08644;
  return p0 / TMath::Power( p1 + TMath::Power(x,p2), p3 );
}

//-------------------------------------------------------------------------
Double_t YMuon( const Double_t *py, const Double_t */*dummy*/ )
{
  // muon y
  
  Double_t y = *py;
  //pol4 only valid in y= -4;-2.5
  Float_t p0,p1,p2,p3;
  p0 = VAR_GENPARAMCUSTOMSINGLE_Y_P0; //0.729545;
  p1 = VAR_GENPARAMCUSTOMSINGLE_Y_P1; //0.53837;
  p2 = VAR_GENPARAMCUSTOMSINGLE_Y_P2; //0.141776;
  p3 = VAR_GENPARAMCUSTOMSINGLE_Y_P3; //0.0130173;
  return p0 * (1. + p1*y + p2*y*y + p3*y*y*y);
}

//-------------------------------------------------------------------------
Double_t V2Muon( const Double_t */*dummy*/, const Double_t */*dummy*/ )
{
  //muon v2
  return 0.;
}

