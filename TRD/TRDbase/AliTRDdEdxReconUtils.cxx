/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
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

//
// TRD dEdx recon utils
// xx
// xx
// xx
// xx
//
//  Xianguo Lu 
//  lu@physi.uni-heidelberg.de
//  Xianguo.Lu@cern.ch
//  
//

#include "TF1.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "THnSparse.h"
#include "TMath.h"
#include "TMatrixD.h"
#include "TMinuit.h"
#include "TObjArray.h"
#include "TRandom3.h"
#include "TStopwatch.h"
#include "TVectorD.h"

#include "TTreeStream.h"

#include "AliCDBId.h"
#include "AliCDBMetaData.h"
#include "AliCDBStorage.h"
#include "AliESDEvent.h"
#include "AliESDfriendTrack.h"
#include "AliESDtrack.h"
#include "AliTRDcalibDB.h"
#include "AliTRDCalROC.h"
#include "AliTRDtrackV1.h"

#include "AliTRDdEdxBaseUtils.h"
#include "AliTRDdEdxReconUtils.h"

#define EPSILON 1e-8

Int_t AliTRDdEdxReconUtils::ApplyCalib(const Int_t nc0, TVectorD *arrayQ, TVectorD *arrayX, const TObjArray *cobj)
{
  //
  //apply calibration on arrayQ
  //
  if(!cobj){ printf("AliTRDdEdxReconUtils::ApplyCalib error gain array null!!\n"); exit(1);}

  TVectorD tmpq(arrayQ->GetNrows());
  TVectorD tmpx(arrayX->GetNrows());
  Int_t ncls = 0;

  const TVectorD * gain = (TVectorD*) cobj->At(0); 
  for(Int_t ii=0; ii<nc0; ii++){
    const Double_t dq = (*arrayQ)[ii];
    const Int_t xx = (Int_t)(*arrayX)[ii];
    const Double_t gg = (*gain)[xx];

    if(gg<EPSILON){
      continue;
    }

    tmpq[ncls] = dq*gg;
    tmpx[ncls] = xx;
    ncls++;
  }

  (*arrayQ)=tmpq;
  (*arrayX)=tmpx;

  return ncls;
}

Double_t AliTRDdEdxReconUtils::ToyCook(const Bool_t kinvq, Int_t &ncluster, TVectorD *arrayQ, TVectorD *arrayX, const TObjArray *cobj)
{
  //
  //template for cookdedx
  //
  if(cobj){
    if(arrayQ && arrayX){
      ncluster = ApplyCalib(ncluster, arrayQ, arrayX, cobj);
    }
    else{
      printf("AliTRDdEdxReconUtils::ToyCook arrayQ arrayX null, applycalib can not be applied!\n"); exit(1);
    }
  }

  Double_t lowFrac =-999, highFrac = -999;
  if(kinvq){
    lowFrac = AliTRDdEdxBaseUtils::Q1Frac(); highFrac = 0.99;
  }
  else{
    lowFrac = 0.01; highFrac = AliTRDdEdxBaseUtils::Q0Frac();
  }

  Double_t meanQ = AliTRDdEdxBaseUtils::TruncatedMean(ncluster, arrayQ->GetMatrixArray(), lowFrac, highFrac);
  if(kinvq){
    if(meanQ>EPSILON){
      meanQ = 1/meanQ;
    }
  }

  return meanQ;
}

Double_t AliTRDdEdxReconUtils::CombineddEdx(const Bool_t kinvq, Int_t &concls, TVectorD *coarrayQ, TVectorD *coarrayX, const Int_t tpcncls, const TVectorD *tpcarrayQ, const TVectorD *tpcarrayX, const Int_t trdncls, const TVectorD *trdarrayQ, const TVectorD *trdarrayX)
{
  //
  //combine tpc and trd dedx
  //

  for(Int_t iq=0; iq<tpcncls; iq++){
    (*coarrayQ)[iq]=(*tpcarrayQ)[iq];
    if(tpcarrayX && trdarrayX && coarrayX){
      (*coarrayX)[iq]=(*tpcarrayX)[iq];
    }
  }
  for(Int_t iq=0; iq<trdncls; iq++){
    (*coarrayQ)[tpcncls+iq]=(*trdarrayQ)[iq];
    if(tpcarrayX && trdarrayX && coarrayX){
      (*coarrayX)[tpcncls+iq]=159+(*trdarrayX)[iq];
    }
  }

  concls=trdncls+tpcncls;

  const Double_t coQ = ToyCook(kinvq, concls, coarrayQ, coarrayX);

  return coQ;
}

Double_t AliTRDdEdxReconUtils::GetPadGain(const Int_t det, const Int_t icol, const Int_t irow)
{
  //
  //get pad calibration
  //
  AliTRDcalibDB* calibration = AliTRDcalibDB::Instance();
  if (!calibration) {
    printf("AliTRDdEdxReconUtils::GetPadCalib No AliTRDcalibDB instance available\n"); exit(1);
  }
  AliTRDCalROC * calGainFactorROC = calibration->GetGainFactorROC(det);
  if(!calGainFactorROC){
    printf("AliTRDdEdxReconUtils::GetPadCalib no calGainFactorROC!\n"); exit(1);
  }

  Double_t padgain = -999;
  if( icol >= 0 && 
      icol < calGainFactorROC->GetNcols() && 
      irow >=0 && 
      irow < calGainFactorROC->GetNrows()){
    padgain = calGainFactorROC->GetValue(icol, irow);
    if(padgain<EPSILON){
      printf("AliTRDdEdxReconUtils::GetPadGain padgain 0! %f %f -- %d %d %d -- %d %d\n", padgain, EPSILON, det, icol, irow, calGainFactorROC->GetNcols(), calGainFactorROC->GetNrows()); exit(1);
    }
  }
  else{
    //printf("\nAliTRDdEdxReconUtils::GetPadGain warning!! indices out of range %d %d %d -- %d %d\n\n", det, icol, irow, calGainFactorROC->GetNcols(), calGainFactorROC->GetNrows() );  
  }

  return padgain;
}

Double_t AliTRDdEdxReconUtils::GetRNDClusterQ(AliTRDcluster *cl, const Double_t baseline)
{
  //
  //get cluter q from GetRawQ, apply baseline and Kr pad-calibration
  //

  const Int_t det     = cl->GetDetector();
  const Int_t pad3col = cl->GetPadCol();
  const Int_t padrow  = cl->GetPadRow();

  Double_t rndqsum = 0;
  for(Int_t ii=0; ii<7; ii++){
    if(cl->GetSignals()[ii] < EPSILON){//bad pad marked by electronics
      continue;
    }

    const Int_t icol = pad3col+(ii-3);
    const Double_t padgain = GetPadGain(det, icol, padrow);
    if(padgain<0){//indices out of range, pad3col near boundary case
      continue;
    }

    const Double_t rndsignal = (cl->GetSignals()[ii] - baseline )/(AliTRDdEdxBaseUtils::IsPadGainOn()? padgain : 1);

    //sum it anyway even if signal below baseline, as long as the total is positive
    rndqsum += rndsignal;
  }

  return rndqsum;
}

Double_t AliTRDdEdxReconUtils::GetClusterQ(const Bool_t kinvq, const AliTRDseedV1 * seed, const Int_t itb)
{
  //
  //get cluster charge
  //
  Double_t dq = 0;
  AliTRDcluster *cl = 0x0;
      
  const Double_t baseline = 10;

  //GetRNDClusterQ(cl)>0 ensures that the total sum of q is above baseline*NsignalPhysical. 
  cl = seed->GetClusters(itb);                    if(cl && GetRNDClusterQ(cl, baseline)>0 ) dq+= GetRNDClusterQ(cl, baseline);
  cl = seed->GetClusters(itb+AliTRDseedV1::kNtb); if(cl && GetRNDClusterQ(cl, baseline)>0 ) dq+= GetRNDClusterQ(cl, baseline);

  dq /= AliTRDdEdxBaseUtils::Getdldx(seed);
  
  dq /= AliTRDdEdxBaseUtils::QScale();
      
  if(kinvq){
    if(dq>EPSILON){
      dq = 1/dq;
    }
  }

  return dq;
}

Int_t AliTRDdEdxReconUtils::GetArrayClusterQ(const Bool_t kinvq, TVectorD *arrayQ, TVectorD *arrayX, const AliTRDtrackV1 *trdtrack, Int_t timeBin0, Int_t timeBin1, Int_t tstep)
{
  //
  //return nclustter
  //(if kinvq, return 1/q array), size of array must be larger than 31*6
  //
  if(!arrayQ || arrayQ->GetNrows()< (AliTRDseedV1::kNtb*AliTRDtrackV1::kNplane)){
    printf("AliTRDdEdxReconUtils::GetArrayClusterQ error arrayQ null or size too small! %d\n", arrayQ? arrayQ->GetNrows() : -999); exit(1);
  }
  if(!arrayX || arrayX->GetNrows()< (AliTRDseedV1::kNtb*AliTRDtrackV1::kNplane)){
    printf("AliTRDdEdxReconUtils::GetArrayClusterQ error arrayX null or size too small! %d\n", arrayX? arrayX->GetNrows() : -999); exit(1);
  }

  const Int_t mintb = 0;
  const Int_t maxtb = AliTRDseedV1::kNtb-1;
  if(timeBin0<mintb) timeBin0=mintb;
  if(timeBin1>maxtb) timeBin1=maxtb;
  if(tstep<=0) tstep=1;

  //============
  Int_t tbN=0;
  Double_t tbQ[200];
  Int_t tbBin[200];
    
  for(Int_t ichamber=0; ichamber < AliTRDtrackV1::kNplane; ichamber++){
    const AliTRDseedV1 * seed = trdtrack->GetTracklet(ichamber);
    if(!seed)
      continue;
    
    const Int_t det = seed->GetDetector();

    for(Int_t itb=timeBin0; itb<=timeBin1; itb+=tstep){
      const Double_t dq = GetClusterQ(kinvq, seed, itb);
      if(dq<EPSILON)
        continue;

      const Int_t gtb = det * AliTRDseedV1::kNtb + itb;

      tbQ[tbN]=dq;
      tbBin[tbN]=gtb;
      tbN++;
    }
  }

  Int_t ncls = 0;
  for(Int_t iq=0; iq<tbN; iq++){
    if(tbQ[iq]<EPSILON)
      continue;

    (*arrayQ)[ncls] = tbQ[iq];
    (*arrayX)[ncls] = tbBin[iq];

    ncls++;
  }

  static Int_t kprint = 100;
  if(kprint<0){
    printf("\nAliTRDdEdxReconUtils::GetArrayClusterQ raw cluster-Q\n");
    for(Int_t iq=0; iq<ncls; iq++){
      const Int_t ichamber =  AliTRDdEdxBaseUtils::ToLayer((*arrayX)[iq]);
      const AliTRDseedV1 * seed = trdtrack->GetTracklet(ichamber);
      if(!seed){
        printf("error seed null!!\n"); exit(1);
      }
      const Double_t rawq =  (*arrayQ)[iq] * 45. * AliTRDdEdxBaseUtils::Getdldx(seed);
      printf("esdid=%d; chamber=%d; timebin=%d; rawq= %.3f; myq[%d]= %e;\n", trdtrack->GetESDid(), ichamber, AliTRDdEdxBaseUtils::ToTimeBin((*arrayX)[iq]), rawq, iq, (*arrayQ)[iq]);
    }
    printf("\n");
  }
  kprint++;

  return ncls;
}

Int_t AliTRDdEdxReconUtils::UpdateArrayX(const Int_t ncls, TVectorD* arrayX)
{
  //
  //arrayX det*Ntb+itb -> itb
  //

  TVectorD countChamber(6);
  for(Int_t ii = 0; ii<ncls; ii++){
    const Int_t xx = (Int_t)(*arrayX)[ii];
    const Int_t idet = AliTRDdEdxBaseUtils::ToDetector(xx);
    
    const Double_t ich = AliTRDgeometry::GetLayer(idet);
    const Double_t itb = AliTRDdEdxBaseUtils::ToTimeBin(xx);
    (*arrayX)[ii] = ich*AliTRDseedV1::kNtb+itb;

    countChamber[ich] = 1;
  }

  const Double_t nch = countChamber.Sum();
  return (Int_t) nch;
}


