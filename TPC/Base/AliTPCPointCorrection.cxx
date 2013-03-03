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

/*  
  
  Unlinearities fitting:

  Model for Outer field cage:
  Unlinearities at the edge aproximated using two exponential decays.
 
  dz = dz0(r,z) +dr(r,z)*tan(theta) 
  dy =          +dr(r,z)*tan(phi)

   
  
    
*/

#include "AliTPCcalibDB.h"
#include "TLinearFitter.h"
#include "Riostream.h"
#include "TList.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TF1.h"
#include "TVectorD.h"
#include "AliLog.h"
#include "AliTPCROC.h"
#include "AliTPCClusterParam.h"
#include "AliTPCPointCorrection.h"

AliTPCPointCorrection* AliTPCPointCorrection::fgInstance = 0;

ClassImp(AliTPCPointCorrection)

AliTPCPointCorrection::AliTPCPointCorrection():
  TNamed(),
  fParamsOutR(38),
  fParamsOutZ(38),
  fParamOutRVersion(0),
  fErrorsOutR(38),
  fErrorsOutZ(38),
  fParamOutZVersion(0),
  //
  fXIO(0),
  fXmiddle(0),
  fXquadrant(0),
  fArraySectorIntParam(36), // array of sector alignment parameters
  fArraySectorIntCovar(36), // array of sector alignment covariances 
  //
  // Kalman filter for global alignment
  //
  fSectorParam(0),     // Kalman parameter   
  fSectorCovar(0)     // Kalman covariance  

{
  //
  // Default constructor
  //
  AliTPCROC * roc = AliTPCROC::Instance();
  fXquadrant = roc->GetPadRowRadii(36,53);
  fXmiddle   = ( roc->GetPadRowRadii(0,0)+roc->GetPadRowRadii(36,roc->GetNRows(36)-1))*0.5;
  fXIO       = ( roc->GetPadRowRadii(0,roc->GetNRows(0)-1)+roc->GetPadRowRadii(36,0))*0.5;
}

AliTPCPointCorrection::AliTPCPointCorrection(const Text_t *name, const Text_t *title):
  TNamed(name,title),
  fParamsOutR(38),
  fParamsOutZ(38),
  fParamOutRVersion(0),
  fErrorsOutR(38),
  fErrorsOutZ(38),
  fParamOutZVersion(0),
  //
  // 
  //
  fXIO(0),
  fXmiddle(0),
  fXquadrant(0),
  fArraySectorIntParam(36), // array of sector alignment parameters
  fArraySectorIntCovar(36), // array of sector alignment covariances 
  //
  // Kalman filter for global alignment
  //
  fSectorParam(0),     // Kalman parameter   for A side
  fSectorCovar(0)     // Kalman covariance  for A side 
  
{
  //
  // Non default constructor
  //
  AliTPCROC * roc = AliTPCROC::Instance();
  fXquadrant = roc->GetPadRowRadii(36,53);
  fXmiddle   = ( roc->GetPadRowRadii(0,0)+roc->GetPadRowRadii(36,roc->GetNRows(36)-1))*0.5;
  fXIO       = ( roc->GetPadRowRadii(0,roc->GetNRows(0)-1)+roc->GetPadRowRadii(36,0))*0.5;

}

AliTPCPointCorrection::~AliTPCPointCorrection(){
  //
  //
  //
}


AliTPCPointCorrection* AliTPCPointCorrection::Instance()
{
  //
  // Singleton implementation
  // Returns an instance of this class, it is created if neccessary
  //
  if (fgInstance == 0){
    fgInstance = new AliTPCPointCorrection();
  }
  return fgInstance;
}



Double_t AliTPCPointCorrection::GetDrOut(Bool_t isGlobal, Bool_t type, Double_t cx, Double_t cy, Double_t cz, Int_t sector){
  //
  //  return radial correction
  //
  if (fParamOutRVersion==0) return CorrectionOutR0(isGlobal, type,cx,cy,cz,sector);
  return 0;
}

Double_t      AliTPCPointCorrection::SGetDrOut(Bool_t isGlobal, Bool_t type, Double_t cx, Double_t cy, Double_t cz, Int_t sector){
  //
  // return radial correction - static function
  // 
  return fgInstance->GetDrOut(isGlobal, type,cx,cy,cz,sector);
}




Double_t AliTPCPointCorrection::GetDzOut(Bool_t isGlobal, Bool_t type, Double_t cx, Double_t cy, Double_t cz, Int_t sector){
  //
  //
  //
  if (fParamOutZVersion==0) return CorrectionOutZ0(isGlobal, type,cx,cy,cz,sector);
  return 0;
}


Double_t      AliTPCPointCorrection::SGetDzOut(Bool_t isGlobal, Bool_t type, Double_t cx, Double_t cy, Double_t cz, Int_t sector){
  //
  //
  //
  return fgInstance->GetDzOut(isGlobal, type,cx,cy,cz,sector);
}




Double_t AliTPCPointCorrection::CorrectionOutR0(Bool_t isGlobal, Bool_t type,  Double_t cx, Double_t cy, Double_t cz, Int_t sector){
  //
  // return dR correction - for correction version 0 
  // Parameters:
  // isGlobal   - is the x in global frame
  // type       - kTRUE - use sectors - kFALSE - only Side param
  // cx, cy,cz  - cluster position
  // sector     - sector number
  if (isGlobal){
    // recalculate sector if in global frame
    Double_t alpha    = TMath::ATan2(cy,cx);
    if (alpha<0)  alpha+=TMath::Pi()*2;
    sector = Int_t(18*alpha/(TMath::Pi()*2));
  }

  if (type==kFALSE) sector=36+(sector%36>=18);  // side level parameters
  // dout
  Double_t radius = TMath::Sqrt(cx*cx+cy*cy);  
  AliTPCROC * roc = AliTPCROC::Instance();
  Double_t xout = roc->GetPadRowRadiiUp(roc->GetNRows(37)-1);
  Double_t dout = xout-radius;  
  if (dout<0) return 0;
  //drift
  Double_t dr   = 0.5 - TMath::Abs(cz/250.);
  //
  //
  TVectorD * vec = GetParamOutR(sector);
  if (!vec) return 0;
  Double_t eout10 = TMath::Exp(-dout/10.);
  Double_t eout5 = TMath::Exp(-dout/5.);		    
  Double_t eoutp  = (eout10+eout5)*0.5;
  Double_t eoutm  = (eout10-eout5)*0.5;
  //
  Double_t result=0;
  result+=(*vec)[1]*eoutp;
  result+=(*vec)[2]*eoutm;
  result+=(*vec)[3]*eoutp*dr;
  result+=(*vec)[4]*eoutm*dr;
  result+=(*vec)[5]*eoutp*dr*dr;
  result+=(*vec)[6]*eoutm*dr*dr;
  return result;
}

Double_t AliTPCPointCorrection::RPhiCOGCorrection(Int_t sector, Int_t padrow, Float_t pad, Float_t cy, Float_t y, Float_t z, Float_t ky,Float_t qMax, Float_t threshold){
  //
  // Calculates COG corection in RPHI direction
  // cluster and track position  y is supposed to be corrected before for other effects
  // (e.g ExB and alignemnt)
  // Rphi distortion dependeds on the distance to the edge-pad, distance to the wire edge and
  // relative distance to the center of the pad. Therefore the y position is trnsfromed to the 
  // pad coordiante frame (correction offset (ExB alignemnt) substracted). 
  //   
  // Input parameters:
  //
  // sector - sector number - 0-71  - cl.GetDetector()
  // padrow - padrow number -0-63 -IROC 0-95 OROC - cl->GetRow()
  // pad    - mean pad number  - cl->GetPad()
  // cy     - cluster y        - cl->GetY()
  //  y     - track -or cluster y
  //  qMax  - cluster max charge - cl-.GetMax()
  //  threshold - clusterer threshold
  //
  AliTPCClusterParam * clparam = AliTPCcalibDB::Instance()->GetClusterParam();
  if (!clparam) {
    AliFatal("TPC OCDB not initialized"); 
    return 0;
  }
  Int_t padtype=0;
  if (sector>=36) padtype = (padrow>64)?2:1;
  Double_t padwidth=(padtype==0)? 0.4:0.6;

  Double_t sigma = clparam->GetRMS0(0,padtype,250-TMath::Abs(z),ky)/padwidth;
  //
  Int_t   npads   =  AliTPCROC::Instance()->GetNPads(sector,padrow);
  Float_t yshift  =  TMath::Abs(cy)-TMath::Abs((-npads*0.5+pad)*padwidth);   // the clusters can be corrected before
                                            // shift to undo correction
  
  y -= yshift*((y>0)?1.:-1.);                             // y in the sector coordinate
  Double_t y0     = npads*0.5*padwidth;
  Double_t dy     = (TMath::Abs(y0)-TMath::Abs(y))/padwidth-0.5;  // distance to  the center of pad0   
  //
  Double_t padcenter = TMath::Nint(dy);
  Double_t sumw=0;
  Double_t sumwi=0;
  for (Double_t ip=padcenter-2.5; ip<=padcenter+2.5;ip++){
    Double_t weightGaus = qMax*TMath::Exp(-(dy-ip)*(dy-ip)/(2*sigma*sigma));
    Double_t ypad       = (ip-npads*0.5)*padwidth;
    Double_t weightGain = GetEdgeQ0(sector,padrow,ypad);
    //
    Double_t weight = TMath::Nint(weightGaus*weightGain);
    if (ip<0 &&weight> threshold) weight = threshold;  // this is done in cl finder
    if (weight<0) continue;
    sumw+=weight;
    sumwi+=weight*(ip);    
  }
  Double_t result =0;
  if (sumw>0) result = sumwi/sumw;
  result = (result-dy)*padwidth;
  return result;
}




Double_t AliTPCPointCorrection::CorrectionOutZ0(Bool_t isGlobal, Bool_t type,  Double_t cx, Double_t cy, Double_t cz, Int_t sector){
  //
  // return dR correction - for correction version 0 
  // Parameters:
  // isGlobal   - is the x in global frame
  // type       - kTRUE - use sectors - kFALSE - only Side param
  // cx, cy,cz  - cluster position
  // sector     - sector number
  if (isGlobal){
    // recalculate sector if in global frame
    Double_t alpha    = TMath::ATan2(cy,cx);
    if (alpha<0)  alpha+=TMath::Pi()*2;
    sector = Int_t(18*alpha/(TMath::Pi()*2));
  }

  if (type==kFALSE) sector=36+(sector%36>=18);  // side level parameters
  // dout
  Double_t radius = TMath::Sqrt(cx*cx+cy*cy);  
  AliTPCROC * roc = AliTPCROC::Instance();
  Double_t xout = roc->GetPadRowRadiiUp(roc->GetNRows(37)-1);
  Double_t dout = xout-radius;
  if (dout<0) return 0;
  //drift
  Double_t dr   = 0.5 - TMath::Abs(cz/250.);
  //
  //
  TVectorD * vec = GetParamOutR(sector);
  if (!vec) return 0;
  Double_t eout10 = TMath::Exp(-dout/10.);
  Double_t eout5 = TMath::Exp(-dout/5.);		    
  Double_t eoutp  = (eout10+eout5)*0.5;
  Double_t eoutm  = (eout10-eout5)*0.5;
  //
  Double_t result=0;
  result+=(*vec)[1]*eoutp;
  result+=(*vec)[2]*eoutm;
  result+=(*vec)[3]*eoutp*dr;
  result+=(*vec)[4]*eoutm*dr;
  result+=(*vec)[5]*eoutp*dr*dr;
  result+=(*vec)[6]*eoutm*dr*dr;
  return result;

}

Double_t  AliTPCPointCorrection::GetEdgeQ0(Int_t sector, Int_t padrow, Float_t y){
  //
  /* TF1 fexp("fexp","1-exp(-[0]*(x-[1]))",0,20)
          | param [0] | param [1]
     IROC | 4.71413   | 1.39558
     OROC1| 2.11437   | 1.52643
     OROC2| 2.15082   | 1.53537 
  */

  Double_t params[2]={100,0};
  if (sector<36){
    params[0]=4.71413; params[1]=1.39558;
  }
  if (sector>=36){
    if (padrow<64) {  params[0]=2.114; params[1]=1.526;}
    if (padrow>=64){  params[0]=2.15; params[1]=1.535;}
  }
  Double_t result= 1;
  Double_t xrow  = AliTPCROC::Instance()->GetPadRowRadii(sector,padrow);
  Double_t ymax  = TMath::Tan(TMath::Pi()/18.)*xrow;
  Double_t dedge = ymax-TMath::Abs(y);
  if (dedge-params[1]<0)             return 0;
  if (dedge>10.*params[1]) return 1;
  result = 1.-TMath::Exp(-params[0]*(dedge-params[1]));
  return result;
}

Double_t AliTPCPointCorrection::SRPhiCOGCorrection(Int_t sector, Int_t padrow, Float_t pad, Float_t cy, Float_t y, Float_t z, Float_t ky,Float_t qMax, Float_t threshold){
  return fgInstance->RPhiCOGCorrection(sector, padrow, pad, cy, y,z, ky, qMax, threshold);
}

Double_t AliTPCPointCorrection::SGetEdgeQ0(Int_t sector, Int_t padrow, Float_t y){
  //
  //
  return fgInstance->GetEdgeQ0(sector, padrow, y);
}

void     AliTPCPointCorrection::AddCorrectionSector(TObjArray & sideAPar, TObjArray &sideCPar, TObjArray & sideACov, TObjArray &sideCCov, Bool_t reset){
  //
  //
  //
  for (Int_t isec=0;isec<36;isec++){
    TMatrixD * mat1 = (TMatrixD*)(sideAPar.At(isec));
    TMatrixD * mat1Covar = (TMatrixD*)(sideACov.At(isec));
    if (!mat1) continue;
    TMatrixD * mat0 = (TMatrixD*)(fArraySectorIntParam.At(isec));
    TMatrixD * mat0Covar = (TMatrixD*)(fArraySectorIntCovar.At(isec));
    if (!mat0) {
      fArraySectorIntParam.AddAt(mat1->Clone(),isec); 
      fArraySectorIntCovar.AddAt(mat1Covar->Clone(),isec); 
      continue;
    }
    (*mat0Covar)=(*mat1Covar);      
    if (reset)   (*mat0)=(*mat1);
    if (!reset) (*mat0)+=(*mat1);
  }

  for (Int_t isec=0;isec<36;isec++){
    TMatrixD * mat1 = (TMatrixD*)(sideCPar.At(isec));
    TMatrixD * mat1Covar = (TMatrixD*)(sideCCov.At(isec));
    if (!mat1) continue;
    TMatrixD * mat0 = (TMatrixD*)(fArraySectorIntParam.At(isec));
    TMatrixD * mat0Covar = (TMatrixD*)(fArraySectorIntCovar.At(isec));
    if (!mat0) {
      fArraySectorIntParam.AddAt(mat1->Clone(),isec); 
      fArraySectorIntCovar.AddAt(mat1Covar->Clone(),isec); 
      continue;
    }
    (*mat0Covar)=(*mat1Covar);      
    if (reset)   (*mat0)=(*mat1);
    if (!reset) (*mat0)+=(*mat1);
  }
}


Double_t AliTPCPointCorrection::GetCorrectionSector(Int_t coord, Int_t sector, Double_t lx, Double_t ly, Double_t /*lz*/, Int_t quadrant){
  //
  // Get position correction for given sector
  //

  TMatrixD * param = (TMatrixD*)fArraySectorIntParam.At(sector%36);
  if (!param) return 0;
  if (quadrant<0){   //recalc quadrant if not specified
    if (lx<fXIO) quadrant=0;
    if(lx>fXIO){
      if (lx<fXquadrant) {
	if (ly<0) quadrant=1;
	if (ly>0) quadrant=2;
      }
      if (lx>=fXquadrant) {
	if (ly<0) quadrant=3;
	if (ly>0) quadrant=4;
      }
    }    
  }
  Double_t a10 = (*param)(quadrant*6+0,0);
  Double_t a20 = (*param)(quadrant*6+1,0);
  Double_t a21 = (*param)(quadrant*6+2,0);
  Double_t dx  = (*param)(quadrant*6+3,0);
  Double_t dy  = (*param)(quadrant*6+4,0);
  Double_t dz  = (*param)(quadrant*6+5,0);
  Double_t deltaX = lx-fXmiddle;
  if (coord==0) return dx;
  if (coord==1) return dy+deltaX*a10;
  if (coord==2) return dz+deltaX*a20+ly*a21;
  if (coord==3) return a10;
  if (coord==4) return a20;
  if (coord==5) return a21;
  //
  return 0;
}

Double_t AliTPCPointCorrection::SGetCorrectionSector(Int_t coord, Int_t sector, Double_t lx, Double_t ly, Double_t lz, Int_t quadrant){
  //
  //
  //
  if (!Instance()) return 0;
  return Instance()->GetCorrectionSector(coord,sector,lx,ly,lz, quadrant);
}



Double_t AliTPCPointCorrection::GetCorrection(Int_t coord, Int_t sector, Double_t lx, Double_t ly, Double_t /*lz*/){
  //
  // Get position correction for given sector
  //
  if (!fSectorParam) return 0;
  
  Double_t a10 = (*fSectorParam)(sector*6+0,0);
  Double_t a20 = (*fSectorParam)(sector*6+1,0);
  Double_t a21 = (*fSectorParam)(sector*6+2,0);
  Double_t dx  = (*fSectorParam)(sector*6+3,0);
  Double_t dy  = (*fSectorParam)(sector*6+4,0);
  Double_t dz  = (*fSectorParam)(sector*6+5,0);
  Double_t deltaX = lx-fXIO;
  //
  if (coord==0) return dx;
  if (coord==1) return dy+deltaX*a10;
  if (coord==2) return dz+deltaX*a20+ly*a21;
  if (coord==3) return a10;
  if (coord==4) return a20;
  if (coord==5) return a21;
  return 0;
}

Double_t AliTPCPointCorrection::SGetCorrection(Int_t coord, Int_t sector, Double_t lx, Double_t ly, Double_t lz){
  //
  //
  //
  if (!Instance()) return 0;
  return Instance()->GetCorrection(coord,sector,lx,ly,lz);
}







