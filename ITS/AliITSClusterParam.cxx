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
#include "TMath.h"
#include "AliITSRecPoint.h"
#include "AliITSClusterParam.h"

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
				  Float_t theta,Float_t phi,
				  Float_t &ny,Float_t &nz)
{
  //
  //get "mean shape"
  //
  if (layer==0){
    ny = 1.+TMath::Abs(phi)*3.2;
    nz = 1.+TMath::Abs(theta)*0.34;
    return;
  }
  if (layer==1){
    ny = 1.+TMath::Abs(phi)*3.2;
    nz = 1.+TMath::Abs(theta)*0.28;
    return;
  }
  
  if (layer>3){
    ny = 2.02+TMath::Abs(phi)*1.95;
    nz = 2.02+TMath::Abs(phi)*2.35;
    return;
  }
  ny  = 6.6-2.7*TMath::Abs(phi);
  nz  = 2.8-3.11*TMath::Abs(phi)+0.45*TMath::Abs(theta);
}
//--------------------------------------------------------------------------
Int_t AliITSClusterParam::GetError(Int_t layer,const AliITSRecPoint*cl,
				   Float_t theta,Float_t phi,Float_t expQ,
				   Float_t &erry,Float_t &errz)
{
  //calculate cluster position error
  //
  Float_t nz,ny;
  GetNTeor(layer, cl,theta,phi,ny,nz);  
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
    Float_t normq = cl->GetQ()/(TMath::Sqrt(1+theta*theta+phi*phi));
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
  Float_t normq = cl->GetQ()/(TMath::Sqrt(1+theta*theta+phi*phi));
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





