/**************************************************************************
 * Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
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

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  ITS cluster error and shape parameterization                             //
//                                                                           //
//  andrea.dainese@lnl.infn.it                                               //
///////////////////////////////////////////////////////////////////////////////
//#include "TFile.h"
//#include "TTree.h"
//#include <TVectorF.h>
//#include <TLinearFitter.h>
//#include <TH1F.h>
//#include <TProfile2D.h>
#include <TVector3.h>
#include "TMath.h"
#include "AliTracker.h"
#include "AliGeomManager.h"
#include "AliITSRecPoint.h"
#include "AliITSClusterParam.h"
#include "AliITSReconstructor.h"
#include "AliExternalTrackParam.h"
#include "AliCheb3DCalc.h"

ClassImp(AliITSClusterParam)


AliITSClusterParam* AliITSClusterParam::fgInstance = 0;


/*
  Example usage fitting parameterization (NOT YET...):
  TFile fres("resol.root");    //tree with resolution and shape 
  TTree * treeRes =(TTree*)fres.Get("Resol");
  
  AliITSClusterParam param;
  param.SetInstance(&param);
  param.FitResol(treeRes);
  param.FitRMS(treeRes);
  TFile fparam("ITSClusterParam.root","recreate");
  param.Write("Param");
  //
  //
  TFile fparam("ITSClusterParam.root");
  AliITSClusterParam *param2  =  (AliITSClusterParam *) fparam.Get("Param"); 
  param2->SetInstance(param2);
  param2->Test(treeRes);
  

  treeRes->Draw("(Resol-AliITSClusterParam::SGetError0(Dim,Pad,Zm,AngleM))/Resol","Dim==0&&QMean<0")
*/

//-------------------------------------------------------------------------
AliITSClusterParam* AliITSClusterParam::Instance()
{
  //
  // Singleton implementation
  // Returns an instance of this class, it is created if neccessary
  //
  if (fgInstance == 0){
    fgInstance = new AliITSClusterParam();
  }
  return fgInstance;
}
//-------------------------------------------------------------------------
void AliITSClusterParam::GetNTeor(Int_t layer,const AliITSRecPoint* /*cl*/, 
				  Float_t tgl,Float_t tgphitr,
				  Float_t &ny,Float_t &nz)
{
  //
  // Get "mean shape" (original parametrization from AliITStrackerMI)
  //
  tgl = TMath::Abs(tgl);
  tgphitr = TMath::Abs(tgphitr);

  // SPD
  if (layer==0) {
    ny = 1.+tgphitr*3.2;
    nz = 1.+tgl*0.34;
    return;
  }
  if (layer==1) {
    ny = 1.+tgphitr*3.2;
    nz = 1.+tgl*0.28;
    return;
  }
  // SSD
  if (layer==4 || layer==5) {
    ny = 2.02+tgphitr*1.95;
    nz = 2.02+tgphitr*2.35;
    return;
  }
  // SDD
  ny  = 6.6-2.7*tgphitr;
  nz  = 2.8-3.11*tgphitr+0.45*tgl;
  return;
}
//--------------------------------------------------------------------------
Int_t AliITSClusterParam::GetError(Int_t layer,
				   const AliITSRecPoint *cl,
				   Float_t tgl,Float_t tgphitr,Float_t expQ,
				   Float_t &erry,Float_t &errz,Float_t &covyz,
				   Bool_t addMisalErr)
{
  //
  // Calculate cluster position error
  //
  static Double_t bz = (Double_t)AliTracker::GetBz();
  Int_t retval=0;
  covyz=0.;
  switch(AliITSReconstructor::GetRecoParam()->GetClusterErrorsParam()) {
  case 0: 
    retval = GetErrorOrigRecPoint(cl,erry,errz,covyz);
    break;
  case 1: 
    retval = GetErrorParamMI(layer,cl,tgl,tgphitr,expQ,erry,errz);
    break;
  case 2: 
    retval = GetErrorParamAngle(layer,cl,tgl,tgphitr,erry,errz);
    break;
  case 3: 
    retval = GetErrorParamAngleOld(layer,cl,tgl,tgphitr,erry,errz);
    break;
  default: 
    retval = GetErrorParamMI(layer,cl,tgl,tgphitr,expQ,erry,errz);
    break;
  }

  // for SSD use the error provided by the cluster finder 
  // if single-sided clusters are enabled
  if(layer>=4 && AliITSReconstructor::GetRecoParam()->GetUseBadChannelsInClusterFinderSSD()) { 
    //printf("error 1 erry errz covyz %10.7f %10.7f %15.13f\n",erry,errz,covyz);
    retval = GetErrorOrigRecPoint(cl,erry,errz,covyz);
    //printf("type %d erry errz covyz %10.7f %10.7f %15.13f\n",cl->GetType(),erry,errz,covyz);
  }
  
  if(addMisalErr) {
    // add error due to misalignment (to be improved)
    Float_t errmisalY2 = AliITSReconstructor::GetRecoParam()->GetClusterMisalErrorY(layer,bz)
      *AliITSReconstructor::GetRecoParam()->GetClusterMisalErrorY(layer,bz);
    Float_t errmisalZ2 = AliITSReconstructor::GetRecoParam()->GetClusterMisalErrorZ(layer,bz)
      *AliITSReconstructor::GetRecoParam()->GetClusterMisalErrorZ(layer,bz);
    erry = TMath::Sqrt(erry*erry+errmisalY2);
    errz = TMath::Sqrt(errz*errz+errmisalZ2);
  }

  return retval;

}
//--------------------------------------------------------------------------
Int_t AliITSClusterParam::GetErrorOrigRecPoint(const AliITSRecPoint*cl,
				   Float_t &erry,Float_t &errz,Float_t &covyz)
{
  //
  // Calculate cluster position error (just take error from AliITSRecPoint)
  //
  erry   = TMath::Sqrt(cl->GetSigmaY2()); 
  errz   = TMath::Sqrt(cl->GetSigmaZ2()); 
  covyz  = cl->GetSigmaYZ();

  return 1;
}
//--------------------------------------------------------------------------
Int_t AliITSClusterParam::GetErrorParamMI(Int_t layer,const AliITSRecPoint*cl,
					  Float_t tgl,Float_t tgphitr,
					  Float_t expQ,
					  Float_t &erry,Float_t &errz)
{
  //
  // Calculate cluster position error (original parametrization from 
  // AliITStrackerMI)
  //
  Float_t nz,ny;
  GetNTeor(layer, cl,tgl,tgphitr,ny,nz);  
  erry   = TMath::Sqrt(cl->GetSigmaY2()); 
  errz   = TMath::Sqrt(cl->GetSigmaZ2()); 
  //
  // PIXELS
  if (layer<2){
    if (TMath::Abs(ny-cl->GetNy())>0.6)  {
      if (ny<cl->GetNy()){
	erry*=0.4+TMath::Abs(ny-cl->GetNy());
	errz*=0.4+TMath::Abs(ny-cl->GetNy());
      }else{
	erry*=0.7+0.5*TMath::Abs(ny-cl->GetNy());
	errz*=0.7+0.5*TMath::Abs(ny-cl->GetNy());
      }
    }
    if (TMath::Abs(nz-cl->GetNz())>1.)  {
      erry*=TMath::Abs(nz-cl->GetNz());
      errz*=TMath::Abs(nz-cl->GetNz());	      
    }
    erry*=0.85;
    errz*=0.85;
    erry= TMath::Min(erry,float(0.0050));
    errz= TMath::Min(errz,float(0.0300));
    return 10;
  }

  //STRIPS
  if (layer>3){ 
    //factor 1.8 appears in new simulation
    //
    Float_t scale=1.8;
    if (cl->GetNy()==100||cl->GetNz()==100){
      erry = 0.004*scale;
      errz = 0.2*scale;
      return 100;
    }
    if (cl->GetNy()+cl->GetNz()>12){
      erry = 0.06*scale;
      errz = 0.57*scale;
      return 100;
    }
    Float_t normq = cl->GetQ()/(TMath::Sqrt(1+tgl*tgl+tgphitr*tgphitr));
    Float_t chargematch = TMath::Max(double(normq/expQ),2.);
    //
    if (cl->GetType()==1 || cl->GetType()==10 ){     							       
      if (chargematch<1.0 || (cl->GetNy()+cl->GetNz()<nz+ny+0.5)){
	errz = 0.043*scale;
	erry = 0.00094*scale;
	return 101;
      }
      if (cl->GetNy()+cl->GetNz()<nz+ny+1.2){
	errz = 0.06*scale;
	erry =0.0013*scale;
	return 102;
      }
      erry = 0.0027*scale;
      errz = TMath::Min(0.028*(chargematch+cl->GetNy()+cl->GetNz()-nz+ny),0.15)*scale;
      return 103;
    }
    if (cl->GetType()==2 || cl->GetType()==11 ){ 
      erry = TMath::Min(0.0010*(1+chargematch+cl->GetNy()+cl->GetNz()-nz+ny),0.05)*scale;
      errz = TMath::Min(0.025*(1+chargematch+cl->GetNy()+cl->GetNz()-nz+ny),0.5)*scale;
      return 104;
    }
    
    if (cl->GetType()>100 ){     							       
      if ((chargematch+cl->GetNy()+cl->GetNz()-nz-ny<1.5)){
	errz = 0.05*scale;
	erry = 0.00096*scale;
	return 105;
      }
      if (cl->GetNy()+cl->GetNz()-nz-ny<1){
	errz = 0.10*scale;
	erry = 0.0025*scale;
	return 106;
      }

      errz = TMath::Min(0.05*(chargematch+cl->GetNy()+cl->GetNz()-nz-ny),0.4)*scale;
      erry = TMath::Min(0.003*(chargematch+cl->GetNy()+cl->GetNz()-nz-ny),0.05)*scale;
      return 107;
    }    
    Float_t diff = cl->GetNy()+cl->GetNz()-ny-nz;
    if (diff<1) diff=1;
    if (diff>4) diff=4;
        
    if (cl->GetType()==5||cl->GetType()==6||cl->GetType()==7||cl->GetType()==8){
      errz = 0.14*diff;
      erry = 0.003*diff;
      return 108;
    }  
    erry = 0.04*diff;
    errz = 0.06*diff;
    return 109;
  }
  //DRIFTS
  Float_t normq = cl->GetQ()/(TMath::Sqrt(1+tgl*tgl+tgphitr*tgphitr));
  Float_t chargematch = normq/expQ;
  chargematch/=2.4; // F. Prino Sept. 2007: SDD charge conversion keV->ADC
  Float_t factorz=1;
  Int_t   cnz = cl->GetNz()%10;
  //charge match
  if (cl->GetType()==1){
    if (chargematch<1.25){
      erry =  0.0028*(1.+6./cl->GetQ());  // gold clusters
    }
    else{
      erry = 0.003*chargematch;
      if (cl->GetNz()==3) erry*=1.5;
    }
    if (chargematch<1.0){
      errz =  0.0011*(1.+6./cl->GetQ());
    }
    else{
      errz = 0.002*(1+2*(chargematch-1.));
    }
    if (cnz>nz+0.6) {
      erry*=(cnz-nz+0.5);
      errz*=1.4*(cnz-nz+0.5);
    }
  }
  if (cl->GetType()>1){
    if (chargematch<1){
      erry =  0.00385*(1.+6./cl->GetQ());  // gold clusters
      errz =  0.0016*(1.+6./cl->GetQ());
    }
    else{
      errz = 0.0014*(1+3*(chargematch-1.));
      erry = 0.003*(1+3*(chargematch-1.));
    } 
    if (cnz>nz+0.6) {
      erry*=(cnz-nz+0.5);
      errz*=1.4*(cnz-nz+0.5);
    }
  }

  if (TMath::Abs(cl->GetY())>2.5){
    factorz*=1+2*(TMath::Abs(cl->GetY())-2.5);
  }
  if (TMath::Abs(cl->GetY())<1){
    factorz*=1.+0.5*TMath::Abs(TMath::Abs(cl->GetY())-1.);
  }
  factorz= TMath::Min(factorz,float(4.));  
  errz*=factorz;

  erry= TMath::Min(erry,float(0.05));
  errz= TMath::Min(errz,float(0.05));  
  return 200;
}
//--------------------------------------------------------------------------
Int_t AliITSClusterParam::GetErrorParamAngle(Int_t layer,
					     const AliITSRecPoint *cl,
					     Float_t tgl,Float_t tgphitr,
					     Float_t &erry,Float_t &errz)
{
  //
  // Calculate cluster position error (parametrization extracted from rp-hit
  // residuals, as a function of angle between track and det module plane.
  // Origin: M.Lunardon, S.Moretto)
  //
  const int   kNcfSPDResX = 21;
  const float kCfSPDResX[kNcfSPDResX] = {+1.1201e+01,+2.0903e+00,-2.2909e-01,-2.6413e-01,+4.2135e-01,-3.7190e-01,
					 +4.2339e-01,+1.8679e-01,-5.1249e-01,+1.8421e-01,+4.8849e-02,-4.3127e-01,
					 -1.1148e-01,+3.1984e-03,-2.5743e-01,-6.6408e-02,+3.0756e-01,+2.6809e-01,
					 -5.0339e-03,-1.4964e-01,-1.1001e-01};
  const float kSPDazMax=56.000000;
  //
  const int   kNcfSPDMeanX = 16;
  const float kCfSPDMeanX[kNcfSPDMeanX] = {-1.2532e+00,-3.8185e-01,-8.9039e-01,+2.6648e+00,+7.0361e-01,+1.2298e+00,
					   +3.2871e-01,+7.8487e-02,-1.6792e-01,-1.3966e-01,-3.1670e-01,-2.1795e-01,
					   -1.9451e-01,-4.9347e-02,-1.9186e-01,-1.9195e-01};
  //
  const int   kNcfSPDResZ = 5;
  const float kCfSPDResZ[kNcfSPDResZ] = {+9.2384e+01,+3.4352e-01,-2.7317e+01,-1.4642e-01,+2.0868e+00};
  const float kSPDpolMin=34.358002, kSPDpolMax=145.000000;
  //
  const Double_t kMaxSigmaSDDx=100.;
  const Double_t kMaxSigmaSDDz=400.;
  const Double_t kMaxSigmaSSDx=100.;
  const Double_t kMaxSigmaSSDz=1000.;
  //  
  const Double_t kParamSDDx[2]={30.93,0.059};
  const Double_t kParamSDDz[2]={33.09,0.011};
  const Double_t kParamSSDx[2]={18.64,-0.0046};
  const Double_t kParamSSDz[2]={784.4,-0.828};
  Double_t sigmax=1000.0,sigmaz=1000.0;
  Double_t biasx = 0.0;

  Int_t volId = (Int_t)cl->GetVolumeId();
  Double_t rotMA[9]; AliGeomManager::GetRotation(volId,rotMA);      // misaligned rotation
  Double_t rotOR[9]; AliGeomManager::GetOrigRotation(volId,rotOR);  // original rotation
  // difference in phi of original and misaligned sensors
  double cross = rotOR[1]*rotMA[4]-rotOR[4]*rotMA[1];
  cross /= TMath::Sqrt( (1.-rotOR[7]*rotOR[7]) * (1.-rotMA[7]*rotMA[7]) );
  Double_t angleAzi = TMath::Abs(TMath::ATan(tgphitr) - TMath::ASin(cross) );
  Double_t anglePol = TMath::Abs(TMath::ATan(tgl));

  if(angleAzi>0.5*TMath::Pi()) angleAzi = TMath::Pi()-angleAzi;
  if(anglePol>0.5*TMath::Pi()) anglePol = TMath::Pi()-anglePol;
  Double_t angleAziDeg = angleAzi*180./TMath::Pi();
  Double_t anglePolDeg = anglePol*180./TMath::Pi();
  
  if(layer==0 || layer==1) { // SPD
    //
    float phiInt    = angleAziDeg/kSPDazMax; // mapped to -1:1
    if (phiInt>1) phiInt = 1; else if (phiInt<-1) phiInt = -1;
    float phiAbsInt = (TMath::Abs(angleAziDeg+angleAziDeg) - kSPDazMax)/kSPDazMax; // mapped to -1:1
    if (phiAbsInt>1) phiAbsInt = 1; else if (phiAbsInt<-1) phiAbsInt = -1;
    anglePolDeg += 90; // the parameterization was provided in polar angle (90 deg - normal to sensor)
    float polInt   = (anglePolDeg+anglePolDeg - (kSPDpolMax+kSPDpolMin))/(kSPDpolMax-kSPDpolMin); // mapped to -1:1
    if (polInt>1) polInt = 1; else if (polInt<-1) polInt = -1;
    //
    sigmax = AliCheb3DCalc::ChebEval1D(phiAbsInt, kCfSPDResX , kNcfSPDResX);
    biasx  = AliCheb3DCalc::ChebEval1D(phiInt   , kCfSPDMeanX, kNcfSPDMeanX);
    sigmaz = AliCheb3DCalc::ChebEval1D(polInt   , kCfSPDResZ , kNcfSPDResZ);
    //
    // for the moment for the SPD only, need to decide where to put it
    biasx *= 1e-4;
    
  } else if(layer==2 || layer==3) { // SDD

    sigmax = angleAziDeg*kParamSDDx[1]+kParamSDDx[0];
    sigmaz = kParamSDDz[0]+kParamSDDz[1]*anglePolDeg;
    if(sigmax > kMaxSigmaSDDx) sigmax = kMaxSigmaSDDx;
    if(sigmaz > kMaxSigmaSDDz) sigmax = kMaxSigmaSDDz;
    
  } else if(layer==4 || layer==5) { // SSD

    sigmax = angleAziDeg*kParamSSDx[1]+kParamSSDx[0];
    sigmaz = kParamSSDz[0]+kParamSSDz[1]*anglePolDeg;
    if(sigmax > kMaxSigmaSSDx) sigmax = kMaxSigmaSSDx;
    if(sigmaz > kMaxSigmaSSDz) sigmax = kMaxSigmaSSDz;
    
  }

  // convert from micron to cm
  erry = 1.e-4*sigmax; 
  errz = 1.e-4*sigmaz;
  

  return 1;
}
//--------------------------------------------------------------------------
Int_t AliITSClusterParam::GetErrorParamAngleOld(Int_t layer,
						const AliITSRecPoint *cl,
						Float_t tgl,Float_t tgphitr,
						Float_t &erry,Float_t &errz)
{
  //
  // Calculate cluster position error (parametrization extracted from rp-hit
  // residuals, as a function of angle between track and det module plane.
  // Origin: M.Lunardon, S.Moretto)
  //

  Double_t maxSigmaSPDx=100.;
  Double_t maxSigmaSPDz=400.;
  Double_t maxSigmaSDDx=100.;
  Double_t maxSigmaSDDz=400.;
  Double_t maxSigmaSSDx=100.;
  Double_t maxSigmaSSDz=1000.;
  
  Double_t paramSPDx[3]={-6.417,0.18,11.14};
  Double_t paramSPDz[2]={118.,-0.155};
  Double_t paramSDDx[2]={30.93,0.059};
  Double_t paramSDDz[2]={33.09,0.011};
  Double_t paramSSDx[2]={18.64,-0.0046};
  Double_t paramSSDz[2]={784.4,-0.828};
  Double_t sigmax=1000.0,sigmaz=1000.0;
  
  Int_t volId = (Int_t)cl->GetVolumeId();
  Double_t rotMA[9]; AliGeomManager::GetRotation(volId,rotMA);      // misaligned rotation
  Double_t rotOR[9]; AliGeomManager::GetOrigRotation(volId,rotOR);  // original rotation
  // difference in phi of original and misaligned sensors
  double cross = rotOR[1]*rotMA[4]-rotOR[4]*rotMA[1];
  cross /= TMath::Sqrt( (1.-rotOR[7]*rotOR[7]) * (1.-rotMA[7]*rotMA[7]) );
  Double_t angleAzi = TMath::Abs(TMath::ATan(tgphitr) - TMath::ASin(cross) );
  Double_t anglePol = TMath::Abs(TMath::ATan(tgl));

  if(angleAzi>0.5*TMath::Pi()) angleAzi = TMath::Pi()-angleAzi;
  if(anglePol>0.5*TMath::Pi()) anglePol = TMath::Pi()-anglePol;
  Double_t angleAziDeg = angleAzi*180./TMath::Pi();
  Double_t anglePolDeg = anglePol*180./TMath::Pi();
  
  if(layer==0 || layer==1) { // SPD

    sigmax = TMath::Exp(angleAziDeg*paramSPDx[1]+paramSPDx[0])+paramSPDx[2];
    sigmaz = paramSPDz[0]+paramSPDz[1]*anglePolDeg;
    if(sigmax > maxSigmaSPDx) sigmax = maxSigmaSPDx;
    if(sigmaz > maxSigmaSPDz) sigmax = maxSigmaSPDz;

  } else if(layer==2 || layer==3) { // SDD

    sigmax = angleAziDeg*paramSDDx[1]+paramSDDx[0];
    sigmaz = paramSDDz[0]+paramSDDz[1]*anglePolDeg;
    if(sigmax > maxSigmaSDDx) sigmax = maxSigmaSDDx;
    if(sigmaz > maxSigmaSDDz) sigmax = maxSigmaSDDz;
    
  } else if(layer==4 || layer==5) { // SSD

    sigmax = angleAziDeg*paramSSDx[1]+paramSSDx[0];
    sigmaz = paramSSDz[0]+paramSSDz[1]*anglePolDeg;
    if(sigmax > maxSigmaSSDx) sigmax = maxSigmaSSDx;
    if(sigmaz > maxSigmaSSDz) sigmax = maxSigmaSSDz;
    
  }

  // convert from micron to cm
  erry = 1.e-4*sigmax; 
  errz = 1.e-4*sigmaz;
  
  return 1;
}
//--------------------------------------------------------------------------
void AliITSClusterParam::Print(Option_t* /*option*/) const {
  //
  // Print param Information
  //

  //
  // Error parameterization
  //
  printf("NOT YET...\n");
  return;
}





