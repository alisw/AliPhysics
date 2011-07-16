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

#include "AliGRPRecoParam.h"

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Class with GRP reconstruction parameters                                  //
// Origin: andrea.dainese@lnl.infn.it                                        //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////



ClassImp(AliGRPRecoParam)

//_____________________________________________________________________________
AliGRPRecoParam::AliGRPRecoParam() : AliDetectorRecoParam(),
fMostProbablePt(0.350),
fVertexerTracksConstraintITS(kTRUE),
fVertexerTracksConstraintTPC(kTRUE),
fVertexerTracksNCuts(12),
fVertexerTracksITSdcacut(0.1),
fVertexerTracksITSdcacutIter0(0.1),
fVertexerTracksITSmaxd0z0(0.5),
fVertexerTracksITSminCls(5),
fVertexerTracksITSmintrks(1),
fVertexerTracksITSnsigma(3.),
fVertexerTracksITSnindetfitter(100.),
fVertexerTracksITSmaxtgl(1000.), 
fVertexerTracksITSfidR(3.),
fVertexerTracksITSfidZ(30.),
fVertexerTracksITSalgo(1.),
fVertexerTracksITSalgoIter0(4.),
fVertexerTracksTPCdcacut(0.1),
fVertexerTracksTPCdcacutIter0(1.0),
fVertexerTracksTPCmaxd0z0(5.),
fVertexerTracksTPCminCls(10),
fVertexerTracksTPCmintrks(1),
fVertexerTracksTPCnsigma(3.),
fVertexerTracksTPCnindetfitter(0.1),
fVertexerTracksTPCmaxtgl(1.5), 
fVertexerTracksTPCfidR(3.),
fVertexerTracksTPCfidZ(30.),
fVertexerTracksTPCalgo(1.),
fVertexerTracksTPCalgoIter0(4.),
fVertexerV0NCuts(7),
fVertexerV0Chi2max(33.),
fVertexerV0DNmin(0.05),
fVertexerV0DPmin(0.05),
fVertexerV0DCAmax(1.5),
fVertexerV0CPAmin(0.9),
fVertexerV0Rmin(0.2),
fVertexerV0Rmax(200.),
fVertexerCascadeNCuts(8),
fVertexerCascadeChi2max(33.),
fVertexerCascadeDV0min(0.01),
fVertexerCascadeMassWin(0.008),
fVertexerCascadeDBachMin(0.01),
fVertexerCascadeDCAmax(2.0),
fVertexerCascadeCPAmin(0.98),
fVertexerCascadeRmin(0.2),
fVertexerCascadeRmax(100.)
{
  //
  // constructor
  //
  SetName("GRP");
  SetTitle("GRP");
}

//_____________________________________________________________________________
AliGRPRecoParam::~AliGRPRecoParam() 
{
  //
  // destructor
  //  
}

AliGRPRecoParam::AliGRPRecoParam(const AliGRPRecoParam& par) :
  AliDetectorRecoParam(par),
  fMostProbablePt(par.fMostProbablePt),
  fVertexerTracksConstraintITS(par.fVertexerTracksConstraintITS),
  fVertexerTracksConstraintTPC(par.fVertexerTracksConstraintTPC),
  fVertexerTracksNCuts(par.fVertexerTracksNCuts),
  fVertexerTracksITSdcacut(par.fVertexerTracksITSdcacut),
  fVertexerTracksITSdcacutIter0(par.fVertexerTracksITSdcacutIter0),
  fVertexerTracksITSmaxd0z0(par.fVertexerTracksITSmaxd0z0),
  fVertexerTracksITSminCls(par.fVertexerTracksITSminCls),
  fVertexerTracksITSmintrks(par.fVertexerTracksITSmintrks),
  fVertexerTracksITSnsigma(par.fVertexerTracksITSnsigma),
  fVertexerTracksITSnindetfitter(par.fVertexerTracksITSnindetfitter),
  fVertexerTracksITSmaxtgl(par.fVertexerTracksITSmaxtgl), 
  fVertexerTracksITSfidR(par.fVertexerTracksITSfidR),
  fVertexerTracksITSfidZ(par.fVertexerTracksITSfidZ),
  fVertexerTracksITSalgo(par.fVertexerTracksITSalgo),
  fVertexerTracksITSalgoIter0(par.fVertexerTracksITSalgoIter0),
  fVertexerTracksTPCdcacut(par.fVertexerTracksTPCdcacut),
  fVertexerTracksTPCdcacutIter0(par.fVertexerTracksTPCdcacutIter0),
  fVertexerTracksTPCmaxd0z0(par.fVertexerTracksTPCmaxd0z0),
  fVertexerTracksTPCminCls(par.fVertexerTracksTPCminCls),
  fVertexerTracksTPCmintrks(par.fVertexerTracksTPCmintrks),
  fVertexerTracksTPCnsigma(par.fVertexerTracksTPCnsigma),
  fVertexerTracksTPCnindetfitter(par.fVertexerTracksTPCnindetfitter),
  fVertexerTracksTPCmaxtgl(par.fVertexerTracksTPCmaxtgl), 
  fVertexerTracksTPCfidR(par.fVertexerTracksTPCfidR),
  fVertexerTracksTPCfidZ(par.fVertexerTracksTPCfidZ),
  fVertexerTracksTPCalgo(par.fVertexerTracksTPCalgo),
  fVertexerTracksTPCalgoIter0(par.fVertexerTracksTPCalgoIter0),
  fVertexerV0NCuts(par.fVertexerV0NCuts),
  fVertexerV0Chi2max(par.fVertexerV0Chi2max),
  fVertexerV0DNmin(par.fVertexerV0DNmin),
  fVertexerV0DPmin(par.fVertexerV0DPmin),
  fVertexerV0DCAmax(par.fVertexerV0DCAmax),
  fVertexerV0CPAmin(par.fVertexerV0CPAmin),
  fVertexerV0Rmin(par.fVertexerV0Rmin),
  fVertexerV0Rmax(par.fVertexerV0Rmax),
  fVertexerCascadeNCuts(par.fVertexerCascadeNCuts),
  fVertexerCascadeChi2max(par.fVertexerCascadeChi2max),
  fVertexerCascadeDV0min(par.fVertexerCascadeDV0min),
  fVertexerCascadeMassWin(par.fVertexerCascadeMassWin),
  fVertexerCascadeDBachMin(par.fVertexerCascadeDBachMin),
  fVertexerCascadeDCAmax(par.fVertexerCascadeDCAmax),
  fVertexerCascadeCPAmin(par.fVertexerCascadeCPAmin),
  fVertexerCascadeRmin(par.fVertexerCascadeRmin),
  fVertexerCascadeRmax(par.fVertexerCascadeRmax)
{
  // copy constructor
}

//_____________________________________________________________________________
AliGRPRecoParam& AliGRPRecoParam::operator = (const AliGRPRecoParam& par)
{
  // assignment operator

  if(&par == this) return *this;

  this->~AliGRPRecoParam();
  new(this) AliGRPRecoParam(par);
  return *this;
}

//_____________________________________________________________________________
AliGRPRecoParam *AliGRPRecoParam::GetHighFluxParam() 
{
  //
  // make default reconstruction  parameters for high flux env.
  //
  AliGRPRecoParam *param = new AliGRPRecoParam();

  // to speed up the vertexing in PbPb
  param->fVertexerTracksITSalgoIter0 = 1.;
  param->fVertexerTracksTPCalgoIter0 = 1.;

  // tighter selections for V0s
  param->fVertexerV0Chi2max = 33.;
  param->fVertexerV0DNmin   = 0.1;
  param->fVertexerV0DPmin   = 0.1;
  param->fVertexerV0DCAmax  = 1.0;
  param->fVertexerV0CPAmin  = 0.998;
  param->fVertexerV0Rmin    = 0.9;
  param->fVertexerV0Rmax    = 100.;

  // tighter selections for Cascades
  param->fVertexerCascadeChi2max  = 33.; 
  param->fVertexerCascadeDV0min   = 0.05;  
  param->fVertexerCascadeMassWin  = 0.008; 
  param->fVertexerCascadeDBachMin = 0.030;
  param->fVertexerCascadeDCAmax   = 0.3;  
  param->fVertexerCascadeCPAmin   = 0.999;  
  param->fVertexerCascadeRmin     = 0.9;    
  param->fVertexerCascadeRmax     = 100.;    

  return param;
}
//_____________________________________________________________________________
AliGRPRecoParam *AliGRPRecoParam::GetLowFluxParam() 
{
  //
  // make default reconstruction  parameters for low  flux env.
  //
  AliGRPRecoParam *param = new AliGRPRecoParam();

  return param;
}
//_____________________________________________________________________________
AliGRPRecoParam *AliGRPRecoParam::GetCosmicTestParam() 
{
  //
  // make default reconstruction  parameters for cosmics env.
  //
  AliGRPRecoParam *param = new AliGRPRecoParam();

  param->SetVertexerTracksConstraintITS(kFALSE);
  param->SetVertexerTracksConstraintTPC(kFALSE);
  param->SetMostProbablePt(3.0);

  return param;
}
//_____________________________________________________________________________
void AliGRPRecoParam::GetVertexerTracksCuts(Int_t mode,Double_t *cuts) const {
  //
  // get cuts for ITS (0) or TPC (1) mode
  //
  if(mode==1) {
    cuts[0] = fVertexerTracksTPCdcacut;
    cuts[1] = fVertexerTracksTPCdcacutIter0;
    cuts[2] = fVertexerTracksTPCmaxd0z0;
    cuts[3] = fVertexerTracksTPCminCls;
    cuts[4] = fVertexerTracksTPCmintrks;
    cuts[5] = fVertexerTracksTPCnsigma;
    cuts[6] = fVertexerTracksTPCnindetfitter;
    cuts[7] = fVertexerTracksTPCmaxtgl; 
    cuts[8] = fVertexerTracksTPCfidR;
    cuts[9] = fVertexerTracksTPCfidZ;
    cuts[10]= fVertexerTracksTPCalgo;
    cuts[11]= fVertexerTracksTPCalgoIter0;
  } else {
    cuts[0] = fVertexerTracksITSdcacut;
    cuts[1] = fVertexerTracksITSdcacutIter0;
    cuts[2] = fVertexerTracksITSmaxd0z0;
    cuts[3] = fVertexerTracksITSminCls;
    cuts[4] = fVertexerTracksITSmintrks;
    cuts[5] = fVertexerTracksITSnsigma;
    cuts[6] = fVertexerTracksITSnindetfitter;
    cuts[7] = fVertexerTracksITSmaxtgl; 
    cuts[8] = fVertexerTracksITSfidR;
    cuts[9] = fVertexerTracksITSfidZ;
    cuts[10]= fVertexerTracksITSalgo;
    cuts[11]= fVertexerTracksITSalgoIter0;
  }

  return;
}
//_____________________________________________________________________________
void AliGRPRecoParam::SetVertexerTracksCuts(Int_t mode,Int_t ncuts,Double_t cuts[12]) {
  //
  // set cuts for ITS (0) or TPC (1) mode
  //
  if(ncuts!=fVertexerTracksNCuts) {
    printf("AliGRPRecoParam: Number of AliVertexerTracks cuts is %d\n",fVertexerTracksNCuts);
    return;
  }

  if(mode==1) {
    fVertexerTracksTPCdcacut = cuts[0];
    fVertexerTracksTPCdcacutIter0 = cuts[1];
    fVertexerTracksTPCmaxd0z0 = cuts[2];
    fVertexerTracksTPCminCls = cuts[3];
    fVertexerTracksTPCmintrks = cuts[4];
    fVertexerTracksTPCnsigma = cuts[5];
    fVertexerTracksTPCnindetfitter = cuts[6];
    fVertexerTracksTPCmaxtgl = cuts[7]; 
    fVertexerTracksTPCfidR = cuts[8];
    fVertexerTracksTPCfidZ = cuts[9];
    fVertexerTracksTPCalgo = cuts[10];
    fVertexerTracksTPCalgoIter0 = cuts[11];
  } else {
    fVertexerTracksITSdcacut = cuts[0];
    fVertexerTracksITSdcacutIter0 = cuts[1];
    fVertexerTracksITSmaxd0z0 = cuts[2];
    fVertexerTracksITSminCls = cuts[3];
    fVertexerTracksITSmintrks = cuts[4];
    fVertexerTracksITSnsigma = cuts[5];
    fVertexerTracksITSnindetfitter = cuts[6];
    fVertexerTracksITSmaxtgl = cuts[7]; 
    fVertexerTracksITSfidR = cuts[8];
    fVertexerTracksITSfidZ = cuts[9];
    fVertexerTracksITSalgo = cuts[10];
    fVertexerTracksITSalgoIter0 = cuts[11];
  }

  return;
}
//_____________________________________________________________________________
void AliGRPRecoParam::GetVertexerV0Cuts(Double_t *cuts) const {
  //
  // get cuts for AliV0vertexer
  //
  cuts[0] = fVertexerV0Chi2max;
  cuts[1] = fVertexerV0DNmin;
  cuts[2] = fVertexerV0DPmin;
  cuts[3] = fVertexerV0DCAmax;
  cuts[4] = fVertexerV0CPAmin;
  cuts[5] = fVertexerV0Rmin;
  cuts[6] = fVertexerV0Rmax;
  return;
}
//_____________________________________________________________________________
void AliGRPRecoParam::SetVertexerV0Cuts(Int_t ncuts,Double_t cuts[7]) {
  //
  // set cuts for AliV0vertexer
  //
  if(ncuts!=fVertexerV0NCuts) {
    printf("AliGRPRecoParam: Number of AliV0vertexer cuts is %d\n",fVertexerV0NCuts);
    return;
  }
  fVertexerV0Chi2max = cuts[0];
  fVertexerV0DNmin   = cuts[1];
  fVertexerV0DPmin   = cuts[2];
  fVertexerV0DCAmax  = cuts[3];
  fVertexerV0CPAmin  = cuts[4];
  fVertexerV0Rmin    = cuts[5];
  fVertexerV0Rmax    = cuts[6];
  return;
}
//_____________________________________________________________________________
void AliGRPRecoParam::GetVertexerCascadeCuts(Double_t *cuts) const {
  //
  // get cuts for AliCascadevertexer
  //
  cuts[0] = fVertexerCascadeChi2max;
  cuts[1] = fVertexerCascadeDV0min;
  cuts[2] = fVertexerCascadeMassWin;
  cuts[3] = fVertexerCascadeDBachMin;
  cuts[4] = fVertexerCascadeDCAmax;
  cuts[5] = fVertexerCascadeCPAmin;
  cuts[6] = fVertexerCascadeRmin;
  cuts[7] = fVertexerCascadeRmax;
  return;
}
//_____________________________________________________________________________
void AliGRPRecoParam::SetVertexerCascadeCuts(Int_t ncuts,Double_t cuts[8]) {
  //
  // set cuts for AliCascadeVertexer
  //
  if(ncuts!=fVertexerCascadeNCuts) {
    printf("AliGRPRecoParam: Number of AliCascadeVertexer cuts is %d\n",fVertexerCascadeNCuts);
    return;
  }
  fVertexerCascadeChi2max  = cuts[0];
  fVertexerCascadeDV0min   = cuts[1];
  fVertexerCascadeMassWin  = cuts[2];
  fVertexerCascadeDBachMin = cuts[3];
  fVertexerCascadeDCAmax   = cuts[4];
  fVertexerCascadeCPAmin   = cuts[5];
  fVertexerCascadeRmin     = cuts[6];
  fVertexerCascadeRmax     = cuts[7];
  return;
}
