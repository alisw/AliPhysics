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
fVertexerTracksNCuts(10),
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
fVertexerTracksTPCdcacut(0.1),
fVertexerTracksTPCdcacutIter0(1.0),
fVertexerTracksTPCmaxd0z0(5.),
fVertexerTracksTPCminCls(10),
fVertexerTracksTPCmintrks(1),
fVertexerTracksTPCnsigma(3.),
fVertexerTracksTPCnindetfitter(0.1),
fVertexerTracksTPCmaxtgl(1.5), 
fVertexerTracksTPCfidR(3.),
fVertexerTracksTPCfidZ(30.)
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
  fVertexerTracksTPCdcacut(par.fVertexerTracksTPCdcacut),
  fVertexerTracksTPCdcacutIter0(par.fVertexerTracksTPCdcacutIter0),
  fVertexerTracksTPCmaxd0z0(par.fVertexerTracksTPCmaxd0z0),
  fVertexerTracksTPCminCls(par.fVertexerTracksTPCminCls),
  fVertexerTracksTPCmintrks(par.fVertexerTracksTPCmintrks),
  fVertexerTracksTPCnsigma(par.fVertexerTracksTPCnsigma),
  fVertexerTracksTPCnindetfitter(par.fVertexerTracksTPCnindetfitter),
  fVertexerTracksTPCmaxtgl(par.fVertexerTracksTPCmaxtgl), 
  fVertexerTracksTPCfidR(par.fVertexerTracksTPCfidR),
  fVertexerTracksTPCfidZ(par.fVertexerTracksTPCfidZ)
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
  // make default reconstruction  parameters for hig  flux env.
  //
  AliGRPRecoParam *param = new AliGRPRecoParam();

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
  }

  return;
}
//_____________________________________________________________________________
void AliGRPRecoParam::SetVertexerTracksCuts(Int_t mode,Int_t ncuts,Double_t cuts[10]) {
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
  }

  return;
}
