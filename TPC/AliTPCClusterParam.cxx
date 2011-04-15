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


///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  TPC cluster error, shape and charge parameterization as function
//  of drift length, and inclination angle                                   //
//
//  Following notation is used in following
//  Int_t dim 0 - y direction
//            1 - z direction
//
//  Int_t type 0 - short pads 
//             1 - medium pads
//             2 - long pads
//  Float_t z    - drift length
// 
//  Float_t angle - tangent of inclination angle at given dimension 
//
//  Implemented parameterization
//
//
//  1. Resolution as function of drift length and inclination angle
//     1.a) GetError0(Int_t dim, Int_t type, Float_t z, Float_t angle)
//          Simple error parameterization as derived from analytical formula
//          only linear term in drift length and angle^2
//          The formula is valid only with precission +-5%
//          Separate parameterization for differnt pad geometry
//     1.b) GetError0Par
//          Parabolic term correction - better precision
//
//     1.c) GetError1 - JUST FOR Study
//          Similar to GetError1
//          The angular and diffusion effect is scaling with pad length
//          common parameterization for different pad length
//
//  2. Error parameterization using charge 
//     2.a) GetErrorQ
//          GetError0+
//          adding 1/Q component to diffusion and angluar part
//     2.b) GetErrorQPar
//          GetError0Par+
//          adding 1/Q component to diffusion and angluar part
//     2.c) GetErrorQParScaled - Just for study
//          One parameterization for all pad shapes
//          Smaller precission as previous one
//
//
//  Example how to retrieve the paramterization:
/*    
      AliCDBManager::Instance()->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
      AliCDBManager::Instance()->SetRun(0) 
      AliTPCClusterParam * param = AliTPCcalibDB::Instance()->GetClusterParam();

      //
      //
      AliTPCClusterParam::SetInstance(param);
      TF1 f1("f1","AliTPCClusterParam::SGetError0Par(1,0,x,0)",0,250);
*/      

// EXAMPLE hot to create parameterization
/*
// Note resol is the resolution tree created by AliTPCcalibTracks
//
AliTPCClusterParam  *param = new AliTPCClusterParam;
param->FitData(Resol);
AliTPCClusterParam::SetInstance(param);
 
*/

//
//                                                                     //
///////////////////////////////////////////////////////////////////////////////
#include "AliTPCClusterParam.h"
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include <TVectorF.h>
#include <TLinearFitter.h>
#include <TH1F.h>
#include <TH3F.h>
#include <TProfile2D.h>
#include <TVectorD.h>
#include <TObjArray.h>
#include "AliTPCcalibDB.h"
#include "AliTPCParam.h"

#include "AliMathBase.h"

ClassImp(AliTPCClusterParam)


AliTPCClusterParam* AliTPCClusterParam::fgInstance = 0;


/*
  Example usage fitting parameterization:
  TFile fres("resol.root");    //tree with resolution and shape 
  TTree * treeRes =(TTree*)fres.Get("Resol");
  
  AliTPCClusterParam param;
  param.SetInstance(&param);
  param.FitResol(treeRes);
  param.FitRMS(treeRes);
  TFile fparam("TPCClusterParam.root","recreate");
  param.Write("Param");
  //
  //
  TFile fparam("TPCClusterParam.root");
  AliTPCClusterParam *param2  =  (AliTPCClusterParam *) fparam.Get("Param"); 
  param2->SetInstance(param2);
  param2->Test(treeRes);
  

  treeRes->Draw("(Resol-AliTPCClusterParam::SGetError0(Dim,Pad,Zm,AngleM))/Resol","Dim==0&&QMean<0")

*/




//_ singleton implementation __________________________________________________
AliTPCClusterParam* AliTPCClusterParam::Instance()
{
  //
  // Singleton implementation
  // Returns an instance of this class, it is created if neccessary
  //
  if (fgInstance == 0){
    fgInstance = new AliTPCClusterParam();
  }
  return fgInstance;
}


AliTPCClusterParam::AliTPCClusterParam():
  TObject(),
  fRatio(0),
  fQNorm(0),
  fQNormCorr(0),
  fQNormHis(0),
  fQpadTnorm(0),           // q pad normalization - Total charge
  fQpadMnorm(0)            // q pad normalization - Max charge
  //
{
  //
  // Default constructor
  //
  fPosQTnorm[0] = 0;   fPosQTnorm[1] = 0;   fPosQTnorm[2] = 0; 
  fPosQMnorm[0] = 0;   fPosQMnorm[1] = 0;   fPosQMnorm[2] = 0; 
  //
  fPosYcor[0]   = 0;   fPosYcor[1]   = 0;   fPosYcor[2]   = 0; 
  fPosZcor[0]   = 0;   fPosZcor[1]   = 0;   fPosZcor[2]   = 0; 
  fErrorRMSSys[0]=0;   fErrorRMSSys[1]=0; 
}

AliTPCClusterParam::AliTPCClusterParam(const AliTPCClusterParam& param):
  TObject(param),
  fRatio(0),
  fQNorm(0),
  fQNormCorr(0),
  fQNormHis(0),
  fQpadTnorm(new TVectorD(*(param.fQpadTnorm))),           // q pad normalization - Total charge
  fQpadMnorm(new TVectorD(*(param.fQpadMnorm)))            // q pad normalization - Max charge

{
  //
  // copy constructor
  //
  memcpy(this, &param,sizeof(AliTPCClusterParam));
  if (param.fQNorm) fQNorm = (TObjArray*) param.fQNorm->Clone();
  if (param.fQNormHis) fQNormHis = (TObjArray*) param.fQNormHis->Clone();
  //
  if (param.fPosQTnorm[0]){
    fPosQTnorm[0] = new TVectorD(*(param.fPosQTnorm[0]));
    fPosQTnorm[1] = new TVectorD(*(param.fPosQTnorm[1]));
    fPosQTnorm[2] = new TVectorD(*(param.fPosQTnorm[2]));
    //
    fPosQMnorm[0] = new TVectorD(*(param.fPosQMnorm[0]));
    fPosQMnorm[1] = new TVectorD(*(param.fPosQMnorm[1]));
    fPosQMnorm[2] = new TVectorD(*(param.fPosQMnorm[2]));
  }
  if (param.fPosYcor[0]){
    fPosYcor[0] = new TVectorD(*(param.fPosYcor[0]));
    fPosYcor[1] = new TVectorD(*(param.fPosYcor[1]));
    fPosYcor[2] = new TVectorD(*(param.fPosYcor[2]));
    //
    fPosZcor[0] = new TVectorD(*(param.fPosZcor[0]));
    fPosZcor[1] = new TVectorD(*(param.fPosZcor[1]));
    fPosZcor[2] = new TVectorD(*(param.fPosZcor[2]));
  }
  
}


AliTPCClusterParam & AliTPCClusterParam::operator=(const AliTPCClusterParam& param){
  //
  // Assignment operator
  //
  if (this != &param) {
    memcpy(this, &param,sizeof(AliTPCClusterParam));
    if (param.fQNorm) fQNorm = (TObjArray*) param.fQNorm->Clone();
    if (param.fQNormHis) fQNormHis = (TObjArray*) param.fQNormHis->Clone();
    if (param.fPosQTnorm[0]){
      fPosQTnorm[0] = new TVectorD(*(param.fPosQTnorm[0]));
      fPosQTnorm[1] = new TVectorD(*(param.fPosQTnorm[1]));
      fPosQTnorm[2] = new TVectorD(*(param.fPosQTnorm[2]));
      //
      fPosQMnorm[0] = new TVectorD(*(param.fPosQMnorm[0]));
      fPosQMnorm[1] = new TVectorD(*(param.fPosQMnorm[1]));
      fPosQMnorm[2] = new TVectorD(*(param.fPosQMnorm[2]));
    }
    if (param.fPosYcor[0]){
      fPosYcor[0] = new TVectorD(*(param.fPosYcor[0]));
      fPosYcor[1] = new TVectorD(*(param.fPosYcor[1]));
      fPosYcor[2] = new TVectorD(*(param.fPosYcor[2]));
      //
      fPosZcor[0] = new TVectorD(*(param.fPosZcor[0]));
      fPosZcor[1] = new TVectorD(*(param.fPosZcor[1]));
      fPosZcor[2] = new TVectorD(*(param.fPosZcor[2]));
    }
  }
  return *this;
}


AliTPCClusterParam::~AliTPCClusterParam(){
  //
  // destructor
  //
  if (fQNorm) fQNorm->Delete();
  if (fQNormCorr) delete fQNormCorr;
  if (fQNormHis) fQNormHis->Delete();
  delete fQNorm;
  delete fQNormHis;
  if (fPosQTnorm[0]){
    delete fPosQTnorm[0];
    delete fPosQTnorm[1];
    delete fPosQTnorm[2];
    //
    delete fPosQMnorm[0];
    delete fPosQMnorm[1];
    delete fPosQMnorm[2];
  }
  if (fPosYcor[0]){
    delete fPosYcor[0];
    delete fPosYcor[1];
    delete fPosYcor[2];
    //
    delete fPosZcor[0];
    delete fPosZcor[1];
    delete fPosZcor[2];
  }
}


void AliTPCClusterParam::FitResol0(TTree * tree, Int_t dim, Int_t type, Float_t *param0, Float_t *error){
  //
  // Fit z - angular dependence of resolution 
  //
  // Int_t dim=0, type=0;
  TString varVal;
  varVal="Resol:AngleM:Zm";
  TString varErr;
  varErr="Sigma:AngleS:Zs";
  TString varCut;
  varCut=Form("Dim==%d&&Pad==%d&&QMean<0",dim,type);
  //
  Int_t entries = tree->Draw(varVal.Data(),varCut);
  Float_t px[10000], py[10000], pz[10000];
  Float_t ex[10000], ey[10000], ez[10000];
  //
  tree->Draw(varErr.Data(),varCut);  
  for (Int_t ipoint=0; ipoint<entries; ipoint++){
    ex[ipoint]= tree->GetV3()[ipoint];
    ey[ipoint]= tree->GetV2()[ipoint];
    ez[ipoint]= tree->GetV1()[ipoint];
  } 
  tree->Draw(varVal.Data(),varCut);
  for (Int_t ipoint=0; ipoint<entries; ipoint++){
    px[ipoint]= tree->GetV3()[ipoint];
    py[ipoint]= tree->GetV2()[ipoint];
    pz[ipoint]= tree->GetV1()[ipoint];
  }
  
  //  
  TLinearFitter fitter(3,"hyp2");
  for (Int_t ipoint=0; ipoint<entries; ipoint++){
    Float_t val = pz[ipoint]*pz[ipoint];
    Float_t err = 2*pz[ipoint]*TMath::Sqrt(ez[ipoint]*ez[ipoint]+fRatio*fRatio*pz[ipoint]*pz[ipoint]);
    Double_t x[2];
    x[0] = px[ipoint];
    x[1] = py[ipoint]*py[ipoint];
    fitter.AddPoint(x,val,err);
  }
  fitter.Eval();
  TVectorD param(3);
  fitter.GetParameters(param);
  param0[0] = param[0];
  param0[1] = param[1];
  param0[2] = param[2];
  Float_t chi2 =  fitter.GetChisquare()/entries;
  param0[3] = chi2;
  error[0] = (fitter.GetParError(0)*TMath::Sqrt(chi2));
  error[1] = (fitter.GetParError(1)*TMath::Sqrt(chi2));
  error[2] = (fitter.GetParError(2)*TMath::Sqrt(chi2));
}


void AliTPCClusterParam::FitResol0Par(TTree * tree, Int_t dim, Int_t type, Float_t *param0, Float_t *error){
  //
  // Fit z - angular dependence of resolution 
  //
  // Int_t dim=0, type=0;
 TString varVal;
  varVal="Resol:AngleM:Zm";
 TString varErr;
  varErr="Sigma:AngleS:Zs";
 TString varCut;
  varCut=Form("Dim==%d&&Pad==%d&&QMean<0",dim,type);
  //
  Int_t entries = tree->Draw(varVal.Data(),varCut);
  Float_t px[10000], py[10000], pz[10000];
  Float_t ex[10000], ey[10000], ez[10000];
  //
  tree->Draw(varErr.Data(),varCut);  
  for (Int_t ipoint=0; ipoint<entries; ipoint++){
    ex[ipoint]= tree->GetV3()[ipoint];
    ey[ipoint]= tree->GetV2()[ipoint];
    ez[ipoint]= tree->GetV1()[ipoint];
  } 
  tree->Draw(varVal.Data(),varCut);
  for (Int_t ipoint=0; ipoint<entries; ipoint++){
    px[ipoint]= tree->GetV3()[ipoint];
    py[ipoint]= tree->GetV2()[ipoint];
    pz[ipoint]= tree->GetV1()[ipoint];
  }
  
  //  
  TLinearFitter fitter(6,"hyp5");
  for (Int_t ipoint=0; ipoint<entries; ipoint++){
    Float_t val = pz[ipoint]*pz[ipoint];
    Float_t err = 2*pz[ipoint]*TMath::Sqrt(ez[ipoint]*ez[ipoint]+fRatio*fRatio*pz[ipoint]*pz[ipoint]);
    Double_t x[6];
    x[0] = px[ipoint];
    x[1] = py[ipoint]*py[ipoint];
    x[2] = x[0]*x[0];
    x[3] = x[1]*x[1];
    x[4] = x[0]*x[1];
    fitter.AddPoint(x,val,err);
  }
  fitter.Eval();
  TVectorD param(6);
  fitter.GetParameters(param);
  param0[0] = param[0];
  param0[1] = param[1];
  param0[2] = param[2];
  param0[3] = param[3];
  param0[4] = param[4];
  param0[5] = param[5];
  Float_t chi2 =  fitter.GetChisquare()/entries;
  param0[6] = chi2;
  error[0] = (fitter.GetParError(0)*TMath::Sqrt(chi2));
  error[1] = (fitter.GetParError(1)*TMath::Sqrt(chi2));
  error[2] = (fitter.GetParError(2)*TMath::Sqrt(chi2));
  error[3] = (fitter.GetParError(3)*TMath::Sqrt(chi2));
  error[4] = (fitter.GetParError(4)*TMath::Sqrt(chi2));
  error[5] = (fitter.GetParError(5)*TMath::Sqrt(chi2));
}





void AliTPCClusterParam::FitResol1(TTree * tree, Int_t dim, Float_t *param0, Float_t *error){
  //
  // Fit z - angular dependence of resolution - pad length scaling 
  //
  // Int_t dim=0, type=0;
 TString varVal;
  varVal="Resol:AngleM*sqrt(Length):Zm/Length";
 TString varErr;
  varErr="Sigma:AngleS:Zs";
 TString varCut;
  varCut=Form("Dim==%d&&QMean<0",dim);
  //
  Int_t entries = tree->Draw(varVal.Data(),varCut);
  Float_t px[10000], py[10000], pz[10000];
  Float_t ex[10000], ey[10000], ez[10000];
  //
  tree->Draw(varErr.Data(),varCut);  
  for (Int_t ipoint=0; ipoint<entries; ipoint++){
    ex[ipoint]= tree->GetV3()[ipoint];
    ey[ipoint]= tree->GetV2()[ipoint];
    ez[ipoint]= tree->GetV1()[ipoint];
  } 
  tree->Draw(varVal.Data(),varCut);
  for (Int_t ipoint=0; ipoint<entries; ipoint++){
    px[ipoint]= tree->GetV3()[ipoint];
    py[ipoint]= tree->GetV2()[ipoint];
    pz[ipoint]= tree->GetV1()[ipoint];
  }
  
  //  
  TLinearFitter fitter(3,"hyp2");
  for (Int_t ipoint=0; ipoint<entries; ipoint++){
    Float_t val = pz[ipoint]*pz[ipoint];
    Float_t err = 2*pz[ipoint]*TMath::Sqrt(ez[ipoint]*ez[ipoint]+fRatio*fRatio*pz[ipoint]*pz[ipoint]);
    Double_t x[2];
    x[0] = px[ipoint];
    x[1] = py[ipoint]*py[ipoint];
    fitter.AddPoint(x,val,err);
  }
  fitter.Eval();
  TVectorD param(3);
  fitter.GetParameters(param);
  param0[0] = param[0];
  param0[1] = param[1];
  param0[2] = param[2];
  Float_t chi2 =  fitter.GetChisquare()/entries;
  param0[3] = chi2;
  error[0] = (fitter.GetParError(0)*TMath::Sqrt(chi2));
  error[1] = (fitter.GetParError(1)*TMath::Sqrt(chi2));
  error[2] = (fitter.GetParError(2)*TMath::Sqrt(chi2));
}

void AliTPCClusterParam::FitResolQ(TTree * tree, Int_t dim, Int_t type, Float_t *param0, Float_t *error){
  //
  // Fit z - angular dependence of resolution - Q scaling 
  //
  // Int_t dim=0, type=0;
 TString varVal;
  varVal="Resol:AngleM/sqrt(QMean):Zm/QMean";
  char varVal0[100];
  snprintf(varVal0,100,"Resol:AngleM:Zm");
  //
 TString varErr;
  varErr="Sigma:AngleS:Zs";
 TString varCut;
  varCut=Form("Dim==%d&&Pad==%d&&QMean>0",dim,type);
  //
  Int_t entries = tree->Draw(varVal.Data(),varCut);
  Float_t px[20000], py[20000], pz[20000], pu[20000], pt[20000];
  Float_t ex[20000], ey[20000], ez[20000];
  //
  tree->Draw(varErr.Data(),varCut);  
  for (Int_t ipoint=0; ipoint<entries; ipoint++){
    ex[ipoint]= tree->GetV3()[ipoint];
    ey[ipoint]= tree->GetV2()[ipoint];
    ez[ipoint]= tree->GetV1()[ipoint];
  } 
  tree->Draw(varVal.Data(),varCut);
  for (Int_t ipoint=0; ipoint<entries; ipoint++){
    px[ipoint]= tree->GetV3()[ipoint];
    py[ipoint]= tree->GetV2()[ipoint];
    pz[ipoint]= tree->GetV1()[ipoint];
  }
  tree->Draw(varVal0,varCut);
  for (Int_t ipoint=0; ipoint<entries; ipoint++){
    pu[ipoint]= tree->GetV3()[ipoint];
    pt[ipoint]= tree->GetV2()[ipoint];
  }
  
  //  
  TLinearFitter fitter(5,"hyp4");
  for (Int_t ipoint=0; ipoint<entries; ipoint++){
    Float_t val = pz[ipoint]*pz[ipoint];
    Float_t err = 2*pz[ipoint]*TMath::Sqrt(ez[ipoint]*ez[ipoint]+fRatio*fRatio*pz[ipoint]*pz[ipoint]);
    Double_t x[4];
    x[0] = pu[ipoint];
    x[1] = pt[ipoint]*pt[ipoint];
    x[2] = px[ipoint];
    x[3] = py[ipoint]*py[ipoint];
    fitter.AddPoint(x,val,err);
  }

  fitter.Eval();
  TVectorD param(5);
  fitter.GetParameters(param);
  param0[0] = param[0];
  param0[1] = param[1];
  param0[2] = param[2];
  param0[3] = param[3];
  param0[4] = param[4];
  Float_t chi2 =  fitter.GetChisquare()/entries;
  param0[5] = chi2;
  error[0] = (fitter.GetParError(0)*TMath::Sqrt(chi2));
  error[1] = (fitter.GetParError(1)*TMath::Sqrt(chi2));
  error[2] = (fitter.GetParError(2)*TMath::Sqrt(chi2));
  error[3] = (fitter.GetParError(3)*TMath::Sqrt(chi2));
  error[4] = (fitter.GetParError(4)*TMath::Sqrt(chi2));
}

void AliTPCClusterParam::FitResolQPar(TTree * tree, Int_t dim, Int_t type, Float_t *param0, Float_t *error){
  //
  // Fit z - angular dependence of resolution - Q scaling  - parabolic correction
  //
  // Int_t dim=0, type=0;
 TString varVal;
  varVal="Resol:AngleM/sqrt(QMean):Zm/QMean";
  char varVal0[100];
  snprintf(varVal0,100,"Resol:AngleM:Zm");
  //
 TString varErr;
  varErr="Sigma:AngleS:Zs";
 TString varCut;
  varCut=Form("Dim==%d&&Pad==%d&&QMean>0",dim,type);
  //
  Int_t entries = tree->Draw(varVal.Data(),varCut);
  Float_t px[20000], py[20000], pz[20000], pu[20000], pt[20000];
  Float_t ex[20000], ey[20000], ez[20000];
  //
  tree->Draw(varErr.Data(),varCut);  
  for (Int_t ipoint=0; ipoint<entries; ipoint++){
    ex[ipoint]= tree->GetV3()[ipoint];
    ey[ipoint]= tree->GetV2()[ipoint];
    ez[ipoint]= tree->GetV1()[ipoint];
  } 
  tree->Draw(varVal.Data(),varCut);
  for (Int_t ipoint=0; ipoint<entries; ipoint++){
    px[ipoint]= tree->GetV3()[ipoint];
    py[ipoint]= tree->GetV2()[ipoint];
    pz[ipoint]= tree->GetV1()[ipoint];
  }
  tree->Draw(varVal0,varCut);
  for (Int_t ipoint=0; ipoint<entries; ipoint++){
    pu[ipoint]= tree->GetV3()[ipoint];
    pt[ipoint]= tree->GetV2()[ipoint];
  }
  
  //  
  TLinearFitter fitter(8,"hyp7");
  for (Int_t ipoint=0; ipoint<entries; ipoint++){
    Float_t val = pz[ipoint]*pz[ipoint];
    Float_t err = 2*pz[ipoint]*TMath::Sqrt(ez[ipoint]*ez[ipoint]+fRatio*fRatio*pz[ipoint]*pz[ipoint]);
    Double_t x[7];
    x[0] = pu[ipoint];
    x[1] = pt[ipoint]*pt[ipoint];
    x[2] = x[0]*x[0];
    x[3] = x[1]*x[1];
    x[4] = x[0]*x[1];
    x[5] = px[ipoint];
    x[6] = py[ipoint]*py[ipoint];
    //
    fitter.AddPoint(x,val,err);
  }

  fitter.Eval();
  TVectorD param(8);
  fitter.GetParameters(param);
  param0[0] = param[0];
  param0[1] = param[1];
  param0[2] = param[2];
  param0[3] = param[3];
  param0[4] = param[4];
  param0[5] = param[5];
  param0[6] = param[6];
  param0[7] = param[7];

  Float_t chi2 =  fitter.GetChisquare()/entries;
  param0[8] = chi2;
  error[0] = (fitter.GetParError(0)*TMath::Sqrt(chi2));
  error[1] = (fitter.GetParError(1)*TMath::Sqrt(chi2));
  error[2] = (fitter.GetParError(2)*TMath::Sqrt(chi2));
  error[3] = (fitter.GetParError(3)*TMath::Sqrt(chi2));
  error[4] = (fitter.GetParError(4)*TMath::Sqrt(chi2));
  error[5] = (fitter.GetParError(5)*TMath::Sqrt(chi2));
  error[6] = (fitter.GetParError(6)*TMath::Sqrt(chi2));
  error[7] = (fitter.GetParError(7)*TMath::Sqrt(chi2));
}



void AliTPCClusterParam::FitRMS0(TTree * tree, Int_t dim, Int_t type, Float_t *param0, Float_t *error){
  //
  // Fit z - angular dependence of resolution 
  //
  // Int_t dim=0, type=0;
 TString varVal;
  varVal="RMSm:AngleM:Zm";
 TString varErr;
  varErr="sqrt((1./(100.*sqrt(12.))^2)+RMSe0^2):AngleS:Zs";
 TString varCut;
  varCut=Form("Dim==%d&&Pad==%d&&QMean<0",dim,type);
  //
  Int_t entries = tree->Draw(varVal.Data(),varCut);
  Float_t px[10000], py[10000], pz[10000];
  Float_t ex[10000], ey[10000], ez[10000];
  //
  tree->Draw(varErr.Data(),varCut);  
  for (Int_t ipoint=0; ipoint<entries; ipoint++){
    ex[ipoint]= tree->GetV3()[ipoint];
    ey[ipoint]= tree->GetV2()[ipoint];
    ez[ipoint]= tree->GetV1()[ipoint];
  } 
  tree->Draw(varVal.Data(),varCut);
  for (Int_t ipoint=0; ipoint<entries; ipoint++){
    px[ipoint]= tree->GetV3()[ipoint];
    py[ipoint]= tree->GetV2()[ipoint];
    pz[ipoint]= tree->GetV1()[ipoint];
  }
  
  //  
  TLinearFitter fitter(3,"hyp2");
  for (Int_t ipoint=0; ipoint<entries; ipoint++){
    Float_t val = pz[ipoint]*pz[ipoint];
    Float_t err = 2*pz[ipoint]*TMath::Sqrt(ez[ipoint]*ez[ipoint]+fRatio*fRatio*pz[ipoint]*pz[ipoint]);
    Double_t x[2];
    x[0] = px[ipoint];
    x[1] = py[ipoint]*py[ipoint];
    fitter.AddPoint(x,val,err);
  }
  fitter.Eval();
  TVectorD param(3);
  fitter.GetParameters(param);
  param0[0] = param[0];
  param0[1] = param[1];
  param0[2] = param[2];
  Float_t chi2 =  fitter.GetChisquare()/entries;
  param0[3] = chi2;
  error[0] = (fitter.GetParError(0)*TMath::Sqrt(chi2));
  error[1] = (fitter.GetParError(1)*TMath::Sqrt(chi2));
  error[2] = (fitter.GetParError(2)*TMath::Sqrt(chi2));
}

void AliTPCClusterParam::FitRMS1(TTree * tree, Int_t dim, Float_t *param0, Float_t *error){
  //
  // Fit z - angular dependence of resolution - pad length scaling 
  //
  // Int_t dim=0, type=0;
 TString varVal;
  varVal="RMSm:AngleM*Length:Zm";
 TString varErr;
  varErr="sqrt((1./(100.*sqrt(12.))^2)+RMSe0^2):AngleS:Pad";
 TString varCut;
  varCut=Form("Dim==%d&&QMean<0",dim);
  //
  Int_t entries = tree->Draw(varVal.Data(),varCut);
  Float_t px[10000], py[10000], pz[10000];
  Float_t type[10000], ey[10000], ez[10000];
  //
  tree->Draw(varErr.Data(),varCut);  
  for (Int_t ipoint=0; ipoint<entries; ipoint++){
    type[ipoint] = tree->GetV3()[ipoint];
    ey[ipoint]   = tree->GetV2()[ipoint];
    ez[ipoint]   = tree->GetV1()[ipoint];
  } 
  tree->Draw(varVal.Data(),varCut);
  for (Int_t ipoint=0; ipoint<entries; ipoint++){
    px[ipoint]= tree->GetV3()[ipoint];
    py[ipoint]= tree->GetV2()[ipoint];
    pz[ipoint]= tree->GetV1()[ipoint];
  }
  
  //  
  TLinearFitter fitter(4,"hyp3");
  for (Int_t ipoint=0; ipoint<entries; ipoint++){
    Float_t val = pz[ipoint]*pz[ipoint];
    Float_t err = 2*pz[ipoint]*TMath::Sqrt(ez[ipoint]*ez[ipoint]+fRatio*fRatio*pz[ipoint]*pz[ipoint]);
    Double_t x[3];
    x[0] = (type[ipoint]<0.5)? 0.:1.;
    x[1] = px[ipoint];
    x[2] = py[ipoint]*py[ipoint];
    fitter.AddPoint(x,val,err);
  }
  fitter.Eval();
  TVectorD param(4);
  fitter.GetParameters(param);
  param0[0] = param[0];
  param0[1] = param[0]+param[1];
  param0[2] = param[2];
  param0[3] = param[3];
  Float_t chi2 =  fitter.GetChisquare()/entries;
  param0[4] = chi2;
  error[0] = (fitter.GetParError(0)*TMath::Sqrt(chi2));
  error[1] = (fitter.GetParError(1)*TMath::Sqrt(chi2));
  error[2] = (fitter.GetParError(2)*TMath::Sqrt(chi2));
  error[3] = (fitter.GetParError(3)*TMath::Sqrt(chi2));
}

void AliTPCClusterParam::FitRMSQ(TTree * tree, Int_t dim, Int_t type, Float_t *param0, Float_t *error){
  //
  // Fit z - angular dependence of resolution - Q scaling 
  //
  // Int_t dim=0, type=0;
 TString varVal;
  varVal="RMSm:AngleM/sqrt(QMean):Zm/QMean";
  char varVal0[100];
  snprintf(varVal0,100,"RMSm:AngleM:Zm");
  //
 TString varErr;
  varErr="sqrt((1./(100.*sqrt(12.))^2)+RMSe0^2):AngleS:Zs";
 TString varCut;
  varCut=Form("Dim==%d&&Pad==%d&&QMean>0",dim,type);
  //
  Int_t entries = tree->Draw(varVal.Data(),varCut);
  Float_t px[20000], py[20000], pz[20000], pu[20000], pt[20000];
  Float_t ex[20000], ey[20000], ez[20000];
  //
  tree->Draw(varErr.Data(),varCut);  
  for (Int_t ipoint=0; ipoint<entries; ipoint++){
    ex[ipoint]= tree->GetV3()[ipoint];
    ey[ipoint]= tree->GetV2()[ipoint];
    ez[ipoint]= tree->GetV1()[ipoint];
  } 
  tree->Draw(varVal.Data(),varCut);
  for (Int_t ipoint=0; ipoint<entries; ipoint++){
    px[ipoint]= tree->GetV3()[ipoint];
    py[ipoint]= tree->GetV2()[ipoint];
    pz[ipoint]= tree->GetV1()[ipoint];
  }
  tree->Draw(varVal0,varCut);
  for (Int_t ipoint=0; ipoint<entries; ipoint++){
    pu[ipoint]= tree->GetV3()[ipoint];
    pt[ipoint]= tree->GetV2()[ipoint];
  }
  
  //  
  TLinearFitter fitter(5,"hyp4");
  for (Int_t ipoint=0; ipoint<entries; ipoint++){
    Float_t val = pz[ipoint]*pz[ipoint];
    Float_t err = 2*pz[ipoint]*TMath::Sqrt(ez[ipoint]*ez[ipoint]+fRatio*fRatio*pz[ipoint]*pz[ipoint]);
    Double_t x[4];
    x[0] = pu[ipoint];
    x[1] = pt[ipoint]*pt[ipoint];
    x[2] = px[ipoint];
    x[3] = py[ipoint]*py[ipoint];
    fitter.AddPoint(x,val,err);
  }

  fitter.Eval();
  TVectorD param(5);
  fitter.GetParameters(param);
  param0[0] = param[0];
  param0[1] = param[1];
  param0[2] = param[2];
  param0[3] = param[3];
  param0[4] = param[4];
  Float_t chi2 =  fitter.GetChisquare()/entries;
  param0[5] = chi2;
  error[0] = (fitter.GetParError(0)*TMath::Sqrt(chi2));
  error[1] = (fitter.GetParError(1)*TMath::Sqrt(chi2));
  error[2] = (fitter.GetParError(2)*TMath::Sqrt(chi2));
  error[3] = (fitter.GetParError(3)*TMath::Sqrt(chi2));
  error[4] = (fitter.GetParError(4)*TMath::Sqrt(chi2));
}


void AliTPCClusterParam::FitRMSSigma(TTree * tree, Int_t dim, Int_t type, Float_t *param0, Float_t */*error*/){
  //
  // Fit z - angular dependence of resolution - Q scaling 
  //
  // Int_t dim=0, type=0;
  TString varVal;
  varVal="RMSs:RMSm";
  //
 TString varCut;
  varCut=Form("Dim==%d&&Pad==%d&&QMean<0",dim,type);
  //
  Int_t entries = tree->Draw(varVal.Data(),varCut);
  Float_t px[20000], py[20000];
  //
  tree->Draw(varVal.Data(),varCut);
  for (Int_t ipoint=0; ipoint<entries; ipoint++){
    px[ipoint]= tree->GetV2()[ipoint];
    py[ipoint]= tree->GetV1()[ipoint];
  }
  TLinearFitter fitter(2,"pol1");
  for (Int_t ipoint=0; ipoint<entries; ipoint++){
    Float_t val = py[ipoint];
    Float_t err = fRatio*px[ipoint];
    Double_t x[4];
    x[0] = px[ipoint];
    if (err>0) fitter.AddPoint(x,val,err);
  }
  fitter.Eval();
  param0[0]= fitter.GetParameter(0);
  param0[1]= fitter.GetParameter(1);
}



Float_t  AliTPCClusterParam::GetError0(Int_t dim, Int_t type, Float_t z, Float_t angle) const {
  //
  //
  //
  Float_t value=0;
  value += fParamS0[dim][type][0];
  value += fParamS0[dim][type][1]*z;
  value += fParamS0[dim][type][2]*angle*angle;
  value  = TMath::Sqrt(TMath::Abs(value)); 
  return value;
}


Float_t  AliTPCClusterParam::GetError0Par(Int_t dim, Int_t type, Float_t z, Float_t angle) const {
  //
  //
  //
  Float_t value=0;
  value += fParamS0Par[dim][type][0];
  value += fParamS0Par[dim][type][1]*z;
  value += fParamS0Par[dim][type][2]*angle*angle;
  value += fParamS0Par[dim][type][3]*z*z;
  value += fParamS0Par[dim][type][4]*angle*angle*angle*angle;
  value += fParamS0Par[dim][type][5]*z*angle*angle;
  value  = TMath::Sqrt(TMath::Abs(value)); 
  return value;
}



Float_t  AliTPCClusterParam::GetError1(Int_t dim, Int_t type, Float_t z, Float_t angle) const {
  //
  //
  //
  Float_t value=0;
  Float_t length=0.75;
  if (type==1) length=1;
  if (type==2) length=1.5;
  value += fParamS1[dim][0];
  value += fParamS1[dim][1]*z/length;
  value += fParamS1[dim][2]*angle*angle*length;
  value  = TMath::Sqrt(TMath::Abs(value)); 
  return value;
}

Float_t  AliTPCClusterParam::GetErrorQ(Int_t dim, Int_t type, Float_t z, Float_t angle, Float_t Qmean) const {
  //
  //
  //
  Float_t value=0;
  value += fParamSQ[dim][type][0];
  value += fParamSQ[dim][type][1]*z;
  value += fParamSQ[dim][type][2]*angle*angle;
  value += fParamSQ[dim][type][3]*z/Qmean;
  value += fParamSQ[dim][type][4]*angle*angle/Qmean;
  value  = TMath::Sqrt(TMath::Abs(value)); 
  return value;


}

Float_t  AliTPCClusterParam::GetErrorQPar(Int_t dim, Int_t type, Float_t z, Float_t angle, Float_t Qmean) const {
  //
  //
  //
  Float_t value=0;
  value += fParamSQPar[dim][type][0];
  value += fParamSQPar[dim][type][1]*z;
  value += fParamSQPar[dim][type][2]*angle*angle;
  value += fParamSQPar[dim][type][3]*z*z;
  value += fParamSQPar[dim][type][4]*angle*angle*angle*angle;
  value += fParamSQPar[dim][type][5]*z*angle*angle;
  value += fParamSQPar[dim][type][6]*z/Qmean;
  value += fParamSQPar[dim][type][7]*angle*angle/Qmean;
  value  = TMath::Sqrt(TMath::Abs(value)); 
  return value;


}

Float_t  AliTPCClusterParam::GetErrorQParScaled(Int_t dim, Int_t type, Float_t z, Float_t angle, Float_t Qmean) const {
  //
  //
  //
  Float_t value=0;
  value += fParamSQPar[dim][type][0];
  value += fParamSQPar[dim][type][1]*z;
  value += fParamSQPar[dim][type][2]*angle*angle;
  value += fParamSQPar[dim][type][3]*z*z;
  value += fParamSQPar[dim][type][4]*angle*angle*angle*angle;
  value += fParamSQPar[dim][type][5]*z*angle*angle;
  value += fParamSQPar[dim][type][6]*z/Qmean;
  value += fParamSQPar[dim][type][7]*angle*angle/Qmean;
  Float_t valueMean = GetError0Par(dim,type,z,angle);
  value -= 0.35*0.35*valueMean*valueMean; 
  value  = TMath::Sqrt(TMath::Abs(value)); 
  return value;


}

Float_t  AliTPCClusterParam::GetRMS0(Int_t dim, Int_t type, Float_t z, Float_t angle) const {
  //
  // calculate mean RMS of cluster - z,angle - parameters for each pad and dimension separatelly
  //
  Float_t value=0;
  value += fParamRMS0[dim][type][0];
  value += fParamRMS0[dim][type][1]*z;
  value += fParamRMS0[dim][type][2]*angle*angle;
  value  = TMath::Sqrt(TMath::Abs(value)); 
  return value;
}

Float_t  AliTPCClusterParam::GetRMS1(Int_t dim, Int_t type, Float_t z, Float_t angle) const {
  //
  // calculate mean RMS of cluster - z,angle - pad length scalling
  //
  Float_t value=0;
  Float_t length=0.75;
  if (type==1) length=1;
  if (type==2) length=1.5;
  if (type==0){
    value += fParamRMS1[dim][0];
  }else{
    value += fParamRMS1[dim][1];
  }
  value += fParamRMS1[dim][2]*z;
  value += fParamRMS1[dim][3]*angle*angle*length*length;
  value  = TMath::Sqrt(TMath::Abs(value)); 
  return value;
}

Float_t  AliTPCClusterParam::GetRMSQ(Int_t dim, Int_t type, Float_t z, Float_t angle, Float_t Qmean) const {
  //
  // calculate mean RMS of cluster - z,angle, Q dependence
  //
  Float_t value=0;
  value += fParamRMSQ[dim][type][0];
  value += fParamRMSQ[dim][type][1]*z;
  value += fParamRMSQ[dim][type][2]*angle*angle;
  value += fParamRMSQ[dim][type][3]*z/Qmean;
  value += fParamRMSQ[dim][type][4]*angle*angle/Qmean;
  value  = TMath::Sqrt(TMath::Abs(value)); 
  return value;
}

Float_t  AliTPCClusterParam::GetRMSSigma(Int_t dim, Int_t type, Float_t z, Float_t angle, Float_t Qmean) const {
  //
  // calculates RMS of signal shape fluctuation
  //
  Float_t mean = GetRMSQ(dim,type,z,angle,Qmean);
  Float_t value  = fRMSSigmaFit[dim][type][0];
  value+=  fRMSSigmaFit[dim][type][1]*mean;
  return value;
}

Float_t  AliTPCClusterParam::GetShapeFactor(Int_t dim, Int_t type, Float_t z, Float_t angle, Float_t Qmean, Float_t rmsL, Float_t rmsM) const {
  //
  // calculates vallue - sigma distortion contribution
  //
  Double_t value =0;
  //
  Float_t rmsMeanQ  = GetRMSQ(dim,type,z,angle,Qmean);
  if (rmsL<rmsMeanQ) return value;
  //
  Float_t rmsSigma  = GetRMSSigma(dim,type,z,angle,Qmean);
  //
  if ((rmsM-rmsMeanQ)>2.0*(rmsSigma+fErrorRMSSys[dim])){
    //1.5 sigma cut on mean
    value+= rmsL*rmsL+2*rmsM*rmsM-3*rmsMeanQ*rmsMeanQ;
  }else{
    if ((rmsL-rmsMeanQ)>3.*(rmsSigma+fErrorRMSSys[dim])){
      //3 sigma cut on local
      value+= rmsL*rmsL-rmsMeanQ*rmsMeanQ;
    }
  }
  return TMath::Sqrt(TMath::Abs(value));
}



void AliTPCClusterParam::FitData(TTree * tree){
  //
  // make fits for error param and shape param
  //
  FitResol(tree);
  FitRMS(tree);

}

void AliTPCClusterParam::FitResol(TTree * tree){
  //
  SetInstance(this);
  for (Int_t idir=0;idir<2; idir++){    
    for (Int_t itype=0; itype<3; itype++){
      Float_t param0[10];
      Float_t error0[10];
      // model error param
      FitResol0(tree, idir, itype,param0,error0);
      printf("\nResol\t%d\t%d\tchi2=%f\n",idir,itype,param0[3]);
      printf("%f\t%f\t%f\n", param0[0],param0[1],param0[2]);
      printf("%f\t%f\t%f\n", error0[0],error0[1],error0[2]);
      for (Int_t ipar=0;ipar<4; ipar++){
	fParamS0[idir][itype][ipar] = param0[ipar];	
	fErrorS0[idir][itype][ipar] = param0[ipar];	
      } 
      // error param with parabolic correction
      FitResol0Par(tree, idir, itype,param0,error0);
      printf("\nResolPar\t%d\t%d\tchi2=%f\n",idir,itype,param0[6]);
      printf("%f\t%f\t%f\t%f\t%f\t%f\n", param0[0],param0[1],param0[2],param0[3],param0[4],param0[5]);
      printf("%f\t%f\t%f\t%f\t%f\t%f\n", error0[0],error0[1],error0[2],error0[3],error0[4],error0[5]);
      for (Int_t ipar=0;ipar<7; ipar++){
	fParamS0Par[idir][itype][ipar] = param0[ipar];	
	fErrorS0Par[idir][itype][ipar] = param0[ipar];	
      }
      //
      FitResolQ(tree, idir, itype,param0,error0);
      printf("\nResolQ\t%d\t%d\tchi2=%f\n",idir,itype,param0[5]);
      printf("%f\t%f\t%f\t%f\t%f\n", param0[0],param0[1],param0[2],param0[3],param0[4]);
      printf("%f\t%f\t%f\t%f\t%f\n", error0[0],error0[1],error0[2],error0[3],error0[4]);
      for (Int_t ipar=0;ipar<6; ipar++){
	fParamSQ[idir][itype][ipar] = param0[ipar];	
	fErrorSQ[idir][itype][ipar] = param0[ipar];	
      }
      //
      FitResolQPar(tree, idir, itype,param0,error0);
      printf("\nResolQ\t%d\t%d\tchi2=%f\n",idir,itype,param0[8]);
      printf("%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", param0[0],param0[1],param0[2],param0[3],param0[4],param0[5],param0[6],param0[7]);
      printf("%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", error0[0],error0[1],error0[2],error0[3],error0[4],error0[5],error0[6],error0[7]);
      for (Int_t ipar=0;ipar<9; ipar++){
	fParamSQPar[idir][itype][ipar] = param0[ipar];	
	fErrorSQPar[idir][itype][ipar] = param0[ipar];	
      }
    }
  }
  //
  printf("Resol z-scaled\n");
  for (Int_t idir=0;idir<2; idir++){    
    Float_t param0[4];
    Float_t error0[4];
    FitResol1(tree, idir,param0,error0);
    printf("\nResol\t%d\tchi2=%f\n",idir,param0[3]);
    printf("%f\t%f\t%f\n", param0[0],param0[1],param0[2]);
    printf("%f\t%f\t%f\n", error0[0],error0[1],error0[2]);
    for (Int_t ipar=0;ipar<4; ipar++){
      fParamS1[idir][ipar] = param0[ipar];	
      fErrorS1[idir][ipar] = param0[ipar];	
    }
  }

  for (Int_t idir=0;idir<2; idir++){
    printf("\nDirection %d\n",idir);
    printf("%d\t%f\t%f\t%f\n", -1,fParamS1[idir][0],fParamS1[idir][1],fParamS1[idir][2]);
    for (Int_t itype=0; itype<3; itype++){
      Float_t length=0.75;
      if (itype==1) length=1;
      if (itype==2) length=1.5;
      printf("%d\t%f\t%f\t%f\n", itype,fParamS0[idir][itype][0],fParamS0[idir][itype][1]*TMath::Sqrt(length),fParamS0[idir][itype][2]/TMath::Sqrt(length));
    }
  }  
}



void AliTPCClusterParam::FitRMS(TTree * tree){
  //
  SetInstance(this);
  for (Int_t idir=0;idir<2; idir++){    
    for (Int_t itype=0; itype<3; itype++){
      Float_t param0[6];
      Float_t error0[6];
      FitRMS0(tree, idir, itype,param0,error0);
      printf("\nRMS\t%d\t%d\tchi2=%f\n",idir,itype,param0[3]);
      printf("%f\t%f\t%f\n", param0[0],param0[1],param0[2]);
      printf("%f\t%f\t%f\n", error0[0],error0[1],error0[2]);
      for (Int_t ipar=0;ipar<4; ipar++){
	fParamRMS0[idir][itype][ipar] = param0[ipar];	
	fErrorRMS0[idir][itype][ipar] = param0[ipar];	
      }
      FitRMSQ(tree, idir, itype,param0,error0);
      printf("\nRMSQ\t%d\t%d\tchi2=%f\n",idir,itype,param0[5]);
      printf("%f\t%f\t%f\t%f\t%f\n", param0[0],param0[1],param0[2],param0[3],param0[4]);
      printf("%f\t%f\t%f\t%f\t%f\n", error0[0],error0[1],error0[2],error0[3],error0[4]);
      for (Int_t ipar=0;ipar<6; ipar++){
	fParamRMSQ[idir][itype][ipar] = param0[ipar];	
	fErrorRMSQ[idir][itype][ipar] = param0[ipar];	
      }
    }
  }
  //
  printf("RMS z-scaled\n");
  for (Int_t idir=0;idir<2; idir++){    
    Float_t param0[5];
    Float_t error0[5];
    FitRMS1(tree, idir,param0,error0);
    printf("\nRMS\t%d\tchi2=%f\n",idir,param0[4]);
    printf("%f\t%f\t%f\t%f\n", param0[0],param0[1],param0[2], param0[3]);
    printf("%f\t%f\t%f\t%f\n", error0[0],error0[1],error0[2], error0[3]);
    for (Int_t ipar=0;ipar<5; ipar++){
      fParamRMS1[idir][ipar] = param0[ipar];	
      fErrorRMS1[idir][ipar] = param0[ipar];	
    }
  }

  for (Int_t idir=0;idir<2; idir++){
    printf("\nDirection %d\n",idir);
    printf("%d\t%f\t%f\t%f\t%f\n", -1,fParamRMS1[idir][0],fParamRMS1[idir][1],fParamRMS1[idir][2], fParamRMS1[idir][3]);
    for (Int_t itype=0; itype<3; itype++){
      Float_t length=0.75;
      if (itype==1) length=1;
      if (itype==2) length=1.5;
      if (itype==0) printf("%d\t%f\t\t\t%f\t%f\n", itype,fParamRMS0[idir][itype][0],fParamRMS0[idir][itype][1],fParamRMS0[idir][itype][2]/length);
      if (itype>0) printf("%d\t\t\t%f\t%f\t%f\n", itype,fParamRMS0[idir][itype][0],fParamRMS0[idir][itype][1],fParamRMS0[idir][itype][2]/length);
    }
  }  
  //
  // Fit RMS sigma
  //
  printf("RMS fluctuation  parameterization \n");
  for (Int_t idir=0;idir<2; idir++){    
    for (Int_t itype=0; itype<3; itype++){ 
      Float_t param0[5];
      Float_t error0[5];
      FitRMSSigma(tree, idir,itype,param0,error0); 
      printf("\t%d\t%d\t%f\t%f\n", idir, itype, param0[0],param0[1]);
      for (Int_t ipar=0;ipar<2; ipar++){
	fRMSSigmaFit[idir][itype][ipar] = param0[ipar];	
      }
    }
  } 
  //
  // store systematic error end RMS fluctuation parameterization
  //
  TH1F hratio("hratio","hratio",100,-0.1,0.1);
  tree->Draw("(RMSm-AliTPCClusterParam::SGetRMSQ(Dim,Pad,Zm,AngleM,QMean))/RMSm>>hratio","Dim==0&&QMean>0");
  fErrorRMSSys[0] = hratio.GetRMS();
  tree->Draw("(RMSm-AliTPCClusterParam::SGetRMSQ(Dim,Pad,Zm,AngleM,QMean))/RMSm>>hratio","Dim==1&&QMean>0");
  fErrorRMSSys[1] = hratio.GetRMS();
  TH1F hratioR("hratioR","hratioR",100,0,0.2);
  tree->Draw("RMSs/RMSm>>hratioR","Dim==0&&QMean>0");
  fRMSSigmaRatio[0][0]=hratioR.GetMean();
  fRMSSigmaRatio[0][1]=hratioR.GetRMS();
  tree->Draw("RMSs/RMSm>>hratioR","Dim==1&&QMean>0");
  fRMSSigmaRatio[1][0]=hratioR.GetMean();
  fRMSSigmaRatio[1][1]=hratioR.GetRMS();
}

void AliTPCClusterParam::Test(TTree * tree, const char *output){
  //
  // Draw standard quality histograms to output file
  //
  SetInstance(this);
  TFile f(output,"recreate");
  f.cd();
  //
  // 1D histograms - resolution
  //
  for (Int_t idim=0; idim<2; idim++){
    for (Int_t ipad=0; ipad<3; ipad++){
      char hname1[300];
      char hcut1[300];
      char hexp1[300];
      //
      snprintf(hname1,300,"Delta0 Dir %d Pad %d",idim,ipad);
      snprintf(hcut1,300,"Dim==%d&&QMean<0&&Pad==%d",idim,ipad);
      snprintf(hexp1,300,"(Resol-AliTPCClusterParam::SGetError0(Dim,Pad,Zm,AngleM))/Resol>>%s",hname1);
      TH1F  his1DRel0(hname1, hname1, 100,-0.2, 0.2);
      snprintf(hname1,300,"Dim==%d&&QMean<0&&Pad=%d",idim,ipad);
      tree->Draw(hexp1,hcut1,"");
      his1DRel0.Write();
      //
      snprintf(hname1,300,"Delta0Par Dir %d Pad %d",idim,ipad);
      snprintf(hcut1,300,"Dim==%d&&QMean<0&&Pad==%d",idim,ipad);
      snprintf(hexp1,300,"(Resol-AliTPCClusterParam::SGetError0Par(Dim,Pad,Zm,AngleM))/Resol>>%s",hname1);
      TH1F  his1DRel0Par(hname1, hname1, 100,-0.2, 0.2);
      snprintf(hname1,300,"Dim==%d&&QMean<0&&Pad=%d",idim,ipad);
      tree->Draw(hexp1,hcut1,"");
      his1DRel0Par.Write();
      //
    }
  }
  //
  // 2D histograms - resolution
  //
  for (Int_t idim=0; idim<2; idim++){
    for (Int_t ipad=0; ipad<3; ipad++){
      char hname1[300];
      char hcut1[300];
      char hexp1[300];
      //
      snprintf(hname1,300,"2DDelta0 Dir %d Pad %d",idim,ipad);
      snprintf(hcut1,300,"Dim==%d&&QMean<0&&Pad==%d",idim,ipad);
      snprintf(hexp1,300,"(Resol-AliTPCClusterParam::SGetError0(Dim,Pad,Zm,AngleM))/Resol:AngleM:Zm>>%s",hname1);
      TProfile2D  profDRel0(hname1, hname1, 6,0,250,6,0,1);
      snprintf(hname1,300,"Dim==%d&&QMean<0&&Pad=%d",idim,ipad);
      tree->Draw(hexp1,hcut1,"");
      profDRel0.Write();
      //
      snprintf(hname1,300,"2DDelta0Par Dir %d Pad %d",idim,ipad);
      snprintf(hcut1,300,"Dim==%d&&QMean<0&&Pad==%d",idim,ipad);
      snprintf(hexp1,300,"(Resol-AliTPCClusterParam::SGetError0Par(Dim,Pad,Zm,AngleM))/Resol:AngleM:Zm>>%s",hname1);
      TProfile2D profDRel0Par(hname1, hname1,6,0,250,6,0,1);
      snprintf(hname1,300,"Dim==%d&&QMean<0&&Pad=%d",idim,ipad);
      tree->Draw(hexp1,hcut1,"");
      profDRel0Par.Write();
      //
    }
  }
}



void AliTPCClusterParam::Print(Option_t* /*option*/) const{
  //
  // Print param Information
  //

  //
  // Error parameterization
  //
  printf("\nResolution Scaled factors\n");
  printf("Dir\tPad\tP0\t\tP1\t\tP2\t\tchi2\n");
  printf("Y\tall\t%f\t%f\t%f\t%f\n", TMath::Sqrt(TMath::Abs(fParamS1[0][0])),TMath::Sqrt(TMath::Abs(fParamS1[0][1])),
	 TMath::Sqrt(TMath::Abs(fParamS1[0][2])),TMath::Sqrt(TMath::Abs(fParamS1[0][3])));
  for (Int_t ipad=0; ipad<3; ipad++){
    Float_t length=0.75;
    if (ipad==1) length=1;
    if (ipad==2) length=1.5;    
    printf("\t%d\t%f\t%f\t%f\t%f\n", ipad, 
	   TMath::Sqrt(TMath::Abs(fParamS0[0][ipad][0])),
	   TMath::Sqrt(TMath::Abs(fParamS0[0][ipad][1]*length)),
	   TMath::Sqrt(TMath::Abs(fParamS0[0][ipad][2]/length)),
	   TMath::Sqrt(TMath::Abs(fParamS0[0][ipad][3])));
  }
  for (Int_t ipad=0; ipad<3; ipad++){
    Float_t length=0.75;
    if (ipad==1) length=1;
    if (ipad==2) length=1.5;
    printf("\t%dPar\t%f\t%f\t%f\t%f\n", ipad, 
	   TMath::Sqrt(TMath::Abs(fParamS0Par[0][ipad][0])),
	   TMath::Sqrt(TMath::Abs(fParamS0Par[0][ipad][1]*length)),
	   TMath::Sqrt(TMath::Abs(fParamS0Par[0][ipad][2]/length)),
	   TMath::Sqrt(TMath::Abs(fParamS0Par[0][ipad][6])));
  }
  printf("Z\tall\t%f\t%f\t%f\t%f\n", TMath::Sqrt(TMath::Abs(fParamS1[1][0])),TMath::Sqrt(fParamS1[1][1]),
	 TMath::Sqrt(fParamS1[1][2]), TMath::Sqrt(fParamS1[1][3]));
  
  for (Int_t ipad=0; ipad<3; ipad++){
    Float_t length=0.75;
    if (ipad==1) length=1;
    if (ipad==2) length=1.5;    
    printf("\t%d\t%f\t%f\t%f\t%f\n", ipad, 
	   TMath::Sqrt(TMath::Abs(fParamS0[1][ipad][0])),
	   TMath::Sqrt(TMath::Abs(fParamS0[1][ipad][1]*length)),
	   TMath::Sqrt(TMath::Abs(fParamS0[1][ipad][2]/length)),
	   TMath::Sqrt(TMath::Abs(fParamS0[1][ipad][3])));
  }
  for (Int_t ipad=0; ipad<3; ipad++){
    Float_t length=0.75;
    if (ipad==1) length=1;
    if (ipad==2) length=1.5;        
    printf("\t%dPar\t%f\t%f\t%f\t%f\n", ipad, 
	   TMath::Sqrt(TMath::Abs(TMath::Abs(fParamS0Par[1][ipad][0]))),
	   TMath::Sqrt(TMath::Abs(fParamS0Par[1][ipad][1]*length)),
	   TMath::Sqrt(TMath::Abs(fParamS0Par[1][ipad][2]/length)),
	   TMath::Sqrt(TMath::Abs(fParamS0Par[1][ipad][6])));
  }
  
  //
  // RMS scaling
  //
  printf("\n");
  printf("\nRMS Scaled factors\n");
  printf("Dir\tPad\tP00\t\tP01\t\tP1\t\tP2\t\tchi2\n");
  printf("Y\tall\t%f\t%f\t%f\t%f\t%f\n", 
	 TMath::Sqrt(TMath::Abs(fParamRMS1[0][0])),
	 TMath::Sqrt(TMath::Abs(fParamRMS1[0][1])),
	 TMath::Sqrt(TMath::Abs(fParamRMS1[0][2])),
	 TMath::Sqrt(TMath::Abs(fParamRMS1[0][3])),
	 TMath::Sqrt(TMath::Abs(fParamRMS1[0][4])));
  for (Int_t ipad=0; ipad<3; ipad++){
    Float_t length=0.75;
    if (ipad==1) length=1;
    if (ipad==2) length=1.5;    
    if (ipad==0){
      printf("\t%d\t%f\t%f\t%f\t%f\t%f\n", ipad, 
	     TMath::Sqrt(TMath::Abs(fParamRMS0[0][ipad][0])),
	     0.,
	     TMath::Sqrt(TMath::Abs(fParamRMS0[0][ipad][1])),
	     TMath::Sqrt(TMath::Abs(fParamRMS0[0][ipad][2]/(length*length))),
	     TMath::Sqrt(TMath::Abs(fParamRMS0[0][ipad][3])));
    }else{
      printf("\t%d\t%f\t%f\t%f\t%f\t%f\n", ipad, 
	     0.,
	     TMath::Sqrt(TMath::Abs(fParamRMS0[0][ipad][0])),
	     TMath::Sqrt(TMath::Abs(fParamRMS0[0][ipad][1])),
	     TMath::Sqrt(TMath::Abs(fParamRMS0[0][ipad][2]/(length*length))),
	     TMath::Sqrt(TMath::Abs(fParamRMS0[0][ipad][3])));	
    }
  }
  printf("\n");
  printf("Z\tall\t%f\t%f\t%f\t%f\t%f\n", 
	 TMath::Sqrt(TMath::Abs(fParamRMS1[1][0])),
	 TMath::Sqrt(TMath::Abs(fParamRMS1[1][1])),
	 TMath::Sqrt(TMath::Abs(fParamRMS1[1][2])),
	 TMath::Sqrt(TMath::Abs(fParamRMS1[1][3])),
	 TMath::Sqrt(TMath::Abs(fParamRMS1[1][4])));
  for (Int_t ipad=0; ipad<3; ipad++){
    Float_t length=0.75;
    if (ipad==1) length=1;
    if (ipad==2) length=1.5;    
    if (ipad==0){
      printf("\t%d\t%f\t%f\t%f\t%f\t%f\n", ipad, 
	     TMath::Sqrt(TMath::Abs(fParamRMS0[1][ipad][0])),
	     0.,
	     TMath::Sqrt(TMath::Abs(fParamRMS0[1][ipad][1])),
	     TMath::Sqrt(TMath::Abs(fParamRMS0[1][ipad][2]/(length*length))),
	     TMath::Sqrt(TMath::Abs(fParamRMS0[1][ipad][3])));
    }else{
      printf("\t%d\t%f\t%f\t%f\t%f\t%f\n", ipad, 
	     0.,
	     TMath::Sqrt(TMath::Abs(fParamRMS0[1][ipad][0])),
	     TMath::Sqrt(TMath::Abs(fParamRMS0[1][ipad][1])),
	     TMath::Sqrt(TMath::Abs(fParamRMS0[1][ipad][2]/(length*length))),
	     TMath::Sqrt(TMath::Abs(fParamRMS0[1][ipad][3])));	
    }
  }
}





Float_t AliTPCClusterParam::Qnorm(Int_t ipad, Int_t itype, Float_t dr, Float_t ty, Float_t tz){
  // get Q normalization
  // type - 0 Qtot 1 Qmax
  // ipad - 0 (0.75 cm) ,1 (1 cm), 2 (1.5 cm)
  //
  //expession formula - TString *strq0 = toolkit.FitPlane(chain,"dedxQ.fElements[2]","dr++ty++tz++dr*ty++dr*tz++++dr*dr++ty*tz++ty^2++tz^2","IPad==0",chi2,npoints,param,covar,0,100000);

  if (fQNorm==0) return 0;
  TVectorD * norm = (TVectorD*)fQNorm->At(3*itype+ipad);
  if (!norm) return 0;
  TVectorD &no  = *norm;
  Float_t   res = 
    no[0]+
    no[1]*dr+
    no[2]*ty+
    no[3]*tz+
    no[4]*dr*ty+
    no[5]*dr*tz+
    no[6]*ty*tz+
    no[7]*dr*dr+
    no[8]*ty*ty+
    no[9]*tz*tz;
  res/=no[0];
  return res;
}



Float_t AliTPCClusterParam::QnormHis(Int_t ipad, Int_t itype, Float_t dr, Float_t p2, Float_t p3){
  // get Q normalization
  // type - 0 Qtot 1 Qmax
  // ipad - 0 (0.75 cm) ,1 (1 cm), 2 (1.5 cm)
  //

  if (fQNormHis==0) return 0;
  TH3F * norm = (TH3F*)fQNormHis->At(4*itype+ipad);
  if (!norm) return 1;
  p2=TMath::Abs(p2);
  dr=TMath::Min(dr,Float_t(norm->GetXaxis()->GetXmax()-norm->GetXaxis()->GetBinWidth(0)));
  dr=TMath::Max(dr,Float_t(norm->GetXaxis()->GetXmin()+norm->GetXaxis()->GetBinWidth(0)));
  //
  p2=TMath::Min(p2,Float_t(norm->GetYaxis()->GetXmax()-norm->GetYaxis()->GetBinWidth(0)));
  p2=TMath::Max(p2,Float_t(norm->GetYaxis()->GetXmin()+norm->GetYaxis()->GetBinWidth(0)));
  //
  p3=TMath::Min(p3,Float_t(norm->GetZaxis()->GetXmax()-norm->GetZaxis()->GetBinWidth(0)));
  p3=TMath::Max(p3,Float_t(norm->GetZaxis()->GetXmin()+norm->GetZaxis()->GetBinWidth(0)));
  //
  Double_t res = norm->GetBinContent(norm->FindBin(dr,p2,p3));
  if (res==0) res = norm->GetBinContent(norm->FindBin(0.5,0.5,0.5));  // This is just hack - to be fixed entries without 

  return res;
}



void AliTPCClusterParam::SetQnorm(Int_t ipad, Int_t itype, const TVectorD *const norm){
  //
  // set normalization
  //
  // type - 0 Qtot 1 Qmax
  // ipad - 0 (0.75 cm) ,1 (1 cm), 2 (1.5 cm)
  //

  if (fQNorm==0) fQNorm = new TObjArray(6);
  fQNorm->AddAt(new TVectorD(*norm), itype*3+ipad);
}

void AliTPCClusterParam::ResetQnormCorr(){
  //
  //
  //
  if (!fQNormCorr) fQNormCorr= new TMatrixD(12,6);
  for (Int_t irow=0;irow<12; irow++)
    for (Int_t icol=0;icol<6; icol++){
      (*fQNormCorr)(irow,icol)=1.;             // default - no correction
      if (irow>5) (*fQNormCorr)(irow,icol)=0.; // default - no correction
    } 
}

void AliTPCClusterParam::SetQnormCorr(Int_t ipad, Int_t itype, Int_t corrType, Float_t val){
  //
  // ipad        - pad type
  // itype       - 0- qtot 1-qmax
  // corrType    - 0 - s0y corr     - eff. PRF corr
  //             - 1 - s0z corr     - eff. TRF corr
  //             - 2 - d0y          - eff. diffusion correction y
  //             - 3 - d0z          - eff. diffusion correction
  //             - 4 - eff length   - eff.length - wire pitch + x diffsion
  //             - 5 - pad type normalization
  if (!fQNormCorr) {
    ResetQnormCorr();
  }
  //
  // eff shap parameterization matrix
  //
  // rows
  // itype*3+ipad  - itype=0 qtot itype=1 qmax ipad=0
  // 
  if (itype<2) (*fQNormCorr)(itype*3+ipad, corrType) *= val;  // multiplicative correction
  if (itype>=2) (*fQNormCorr)(itype*3+ipad, corrType)+= val;  // additive       correction  
}

Double_t  AliTPCClusterParam::GetQnormCorr(Int_t ipad, Int_t itype, Int_t corrType) const{
  //
  // see AliTPCClusterParam::SetQnormCorr
  //
  if (!fQNormCorr) return 0;
  return  (*fQNormCorr)(itype*3+ipad, corrType);
}


Float_t AliTPCClusterParam::QnormPos(Int_t ipad,Bool_t isMax, Float_t pad, Float_t time, Float_t z, Float_t sy2, Float_t sz2, Float_t qm, Float_t qt){
  //
  // Make Q normalization as function of following parameters
  // Fit parameters to be used in corresponding correction function extracted in the AliTPCclaibTracksGain - Taylor expansion
  // 1 - dp   - relative pad position 
  // 2 - dt   - relative time position
  // 3 - di   - drift length (norm to 1);
  // 4 - dq0  - Tot/Max charge
  // 5 - dq1  - Max/Tot charge
  // 6 - sy   - sigma y - shape
  // 7 - sz   - sigma z - shape
  //  
  
  //The results can be visualized using the debug streamer information of the AliTPCcalibTracksGain - 
  // Following variable used - correspondance to the our variable conventions  
  //chain0->SetAlias("dp","((Cl.fPad-int(Cl.fPad)-0.5)/0.5)");
  Double_t dp = ((pad-int(pad)-0.5)*2.);
  //chain0->SetAlias("dt","((Cl.fTimeBin-int(Cl.fTimeBin)-0.5)/0.5)");
  Double_t dt = ((time-int(time)-0.5)*2.);
  //chain0->SetAlias("di","(sqrt(1.-abs(Cl.fZ)/250.))");
  Double_t di = TMath::Sqrt(1-TMath::Abs(z)/250.);
  //chain0->SetAlias("dq0","(0.2*(Cl.fQ+2)/(Cl.fMax+2))");
  Double_t dq0 = 0.2*(qt+2.)/(qm+2.);
  //chain0->SetAlias("dq1","(5*(Cl.fMax+2)/(Cl.fQ+2))");
  Double_t dq1 = 5.*(qm+2.)/(qt+2.);
  //chain0->SetAlias("sy","(0.32/sqrt(0.01^2+Cl.fSigmaY2))");
  Double_t sy  = 0.32/TMath::Sqrt(0.01*0.01+sy2);
  //chain0->SetAlias("sz","(0.32/sqrt(0.01^2+Cl.fSigmaZ2))");
  Double_t sz  = 0.32/TMath::Sqrt(0.01*0.01+sz2);
  //
  //
  //
  TVectorD * pvec = 0;
  if (isMax){
    pvec = fPosQMnorm[ipad];
  }else{
    pvec = fPosQTnorm[ipad];    
  }
  TVectorD &param = *pvec;
  //
  // Eval part  - in correspondance with fit part from debug streamer
  // 
  Double_t result=param[0];
  Int_t index =1;
  //
  result+=dp*param[index++];                               //1
  result+=dt*param[index++];                               //2
  result+=dp*dp*param[index++];                             //3
  result+=dt*dt*param[index++];                             //4
  result+=dt*dt*dt*param[index++];                             //5
  result+=dp*dt*param[index++];                            //6
  result+=dp*dt*dt*param[index++];                          //7
  result+=(dq0)*param[index++];                            //8
  result+=(dq1)*param[index++];                            //9
  //
  //
  result+=dp*dp*(di)*param[index++];                        //10
  result+=dt*dt*(di)*param[index++];                        //11
  result+=dp*dp*sy*param[index++];                          //12
  result+=dt*sz*param[index++];                          //13
  result+=dt*dt*sz*param[index++];                          //14
  result+=dt*dt*dt*sz*param[index++];                          //15
  //
  result+=dp*dp*1*sy*sz*param[index++];                     //16
  result+=dt*sy*sz*param[index++];                       //17
  result+=dt*dt*sy*sz*param[index++];                       //18
  result+=dt*dt*dt*sy*sz*param[index++];                       //19
  //
  result+=dp*dp*(dq0)*param[index++];                       //20
  result+=dt*1*(dq0)*param[index++];                       //21
  result+=dt*dt*(dq0)*param[index++];                       //22
  result+=dt*dt*dt*(dq0)*param[index++];                       //23
  //
  result+=dp*dp*(dq1)*param[index++];                       //24
  result+=dt*(dq1)*param[index++];                       //25
  result+=dt*dt*(dq1)*param[index++];                       //26
  result+=dt*dt*dt*(dq1)*param[index++];                       //27

  if (result<0.75) result=0.75;
  if (result>1.25) result=1.25;

  return result;
  
}





Float_t AliTPCClusterParam::PosCorrection(Int_t type, Int_t ipad,  Float_t pad, Float_t time, Float_t z, Float_t /*sy2*/, Float_t /*sz2*/, Float_t /*qm*/){

  //
  // Make postion correction
  // type - 0 - y correction
  //        1 - z correction
  // ipad - 0, 1, 2 - short, medium long pads 
  // pad  - float pad number          
  // time - float time bin number
  //    z - z of the cluster
  
  //
  //chainres->SetAlias("dp","(-1+(Cl.fZ>0)*2)*((Cl.fPad-int(Cl.fPad))-0.5)");
  //chainres->SetAlias("dt","(-1+(Cl.fZ>0)*2)*((Cl.fTimeBin-0.66-int(Cl.fTimeBin-0.66))-0.5)");
  //chainres->SetAlias("sp","(sin(dp*pi)-dp*pi)");
  //chainres->SetAlias("st","(sin(dt)-dt)");
  //
  //chainres->SetAlias("di","sqrt(1.-abs(Cl.fZ/250.))");

  //
  // Derived variables
  //
  Double_t dp = (-1+(z>0)*2)*((pad-int(pad))-0.5);
  Double_t dt = (-1+(z>0)*2)*((time-0.66-int(time-0.66))-0.5);
  Double_t sp = (TMath::Sin(dp*TMath::Pi())-dp*TMath::Pi());
  Double_t st = (TMath::Sin(dt)-dt);
  //
  Double_t di = TMath::Sqrt(TMath::Abs(1.-TMath::Abs(z/250.)));
  //
  //
  //
  TVectorD * pvec = 0;
  if (type==0){
    pvec = fPosYcor[ipad];
  }else{
    pvec = fPosZcor[ipad];    
  }
  TVectorD &param = *pvec;
  //
  Double_t result=0;
  Int_t index =1;

  if (type==0){
    // y corr
    result+=(dp)*param[index++];             //1
    result+=(dp)*di*param[index++];          //2
    //
    result+=(sp)*param[index++];             //3
    result+=(sp)*di*param[index++];          //4
  }
  if (type==1){
    result+=(dt)*param[index++];             //1
    result+=(dt)*di*param[index++];          //2
    //
    result+=(st)*param[index++];             //3
    result+=(st)*di*param[index++];          //4
  }
  if (TMath::Abs(result)>0.05) return 0;
  return result;
}



Double_t  AliTPCClusterParam::GaussConvolution(Double_t x0, Double_t x1, Double_t k0, Double_t k1, Double_t s0, Double_t s1){
  //
  // 2 D gaus convoluted with angular effect
  // See in mathematica: 
  //Simplify[Integrate[Exp[-(x0-k0*xd)*(x0-k0*xd)/(2*s0*s0)-(x1-k1*xd)*(x1-k1*xd)/(2*s1*s1)]/(s0*s1),{xd,-1/2,1/2}]]
  // 
  //TF1 f1("f1","AliTPCClusterParam::GaussConvolution(x,0,1,0,0.1,0.1)",-2,2)
  //TF2 f2("f2","AliTPCClusterParam::GaussConvolution(x,y,1,1,0.1,0.1)",-2,2,-2,2)
  //
  const Double_t kEpsilon = 0.0001;
  const Double_t twoPi = TMath::TwoPi();
  const Double_t hnorm = 0.5/TMath::Sqrt(twoPi);
  const Double_t sqtwo = TMath::Sqrt(2.);

  if ((TMath::Abs(k0)+TMath::Abs(k1))<kEpsilon*(s0+s1)){
    // small angular effect
    Double_t val = TMath::Gaus(x0,0,s0)*TMath::Gaus(x1,0,s1)/(s0*s1*twoPi);
    return val;
  }
  Double_t sigma2 = k1*k1*s0*s0+k0*k0*s1*s1;
  Double_t sigma = TMath::Sqrt(sigma2);
  Double_t exp0 = TMath::Exp(-(k1*x0-k0*x1)*(k1*x0-k0*x1)/(2.*sigma2));
  //
  Double_t sigmaErf =  1./(2.*s0*s1*sqtwo*sigma);
  Double_t k0s1s1 = 2.*k0*s1*s1;
  Double_t k1s0s0 = 2.*k1*s0*s0;
  Double_t erf0 = AliMathBase::ErfFast((sigma2-k0s1s1*x0-k1s0s0*x1)*sigmaErf);
  Double_t erf1 = AliMathBase::ErfFast((sigma2+k0s1s1*x0+k1s0s0*x1)*sigmaErf);
  Double_t norm = hnorm/sigma;
  Double_t val  = norm*exp0*(erf0+erf1);
  return val;
}


Double_t  AliTPCClusterParam::GaussConvolutionTail(Double_t x0, Double_t x1, Double_t k0, Double_t k1, Double_t s0, Double_t s1, Double_t tau){
  //
  // 2 D gaus convoluted with angular effect and exponential tail in z-direction
  // tail integrated numerically 
  // Integral normalized to one
  // Mean at 0
  // 
  // TF1 f1t("f1t","AliTPCClusterParam::GaussConvolutionTail(0,x,0,0,0.5,0.5,0.9)",-5,5)
  Double_t sum =1, mean=0;
  // the COG of exponent
  for (Float_t iexp=0;iexp<5;iexp+=0.2){
    mean+=iexp*TMath::Exp(-iexp/tau);
    sum +=TMath::Exp(-iexp/tau);
  }
  mean/=sum;
  //
  sum = 1;
  Double_t val = GaussConvolution(x0,x1+mean, k0, k1 , s0,s1);
  for (Float_t iexp=0;iexp<5;iexp+=0.2){
    val+=GaussConvolution(x0,x1+mean-iexp, k0, k1 , s0,s1)*TMath::Exp(-iexp/tau);
    sum+=TMath::Exp(-iexp/tau);
  }
  return val/sum;
}

Double_t  AliTPCClusterParam::GaussConvolutionGamma4(Double_t x0, Double_t x1, Double_t k0, Double_t k1, Double_t s0, Double_t s1, Double_t tau){
  //
  // 2 D gaus convoluted with angular effect and exponential tail in z-direction
  // tail integrated numerically 
  // Integral normalized to one
  // Mean at 0
  // 
  // TF1 f1g4("f1g4","AliTPCClusterParam::GaussConvolutionGamma4(0,x,0,0,0.5,0.2,1.6)",-5,5)
  // TF2 f2g4("f2g4","AliTPCClusterParam::GaussConvolutionGamma4(y,x,0,0,0.5,0.2,1.6)",-5,5,-5,5)
  Double_t sum =0, mean=0;
  // the COG of G4
  for (Float_t iexp=0;iexp<5;iexp+=0.2){
    Double_t g4 = TMath::Exp(-4.*iexp/tau)*TMath::Power(iexp/tau,4.);
    mean+=iexp*g4;
    sum +=g4;
  }
  mean/=sum;
  //
  sum = 0;
  Double_t val = 0;
  for (Float_t iexp=0;iexp<5;iexp+=0.2){ 
    Double_t g4 = TMath::Exp(-4.*iexp/tau)*TMath::Power(iexp/tau,4.);
    val+=GaussConvolution(x0,x1+mean-iexp, k0, k1 , s0,s1)*g4;
    sum+=g4;
  }
  return val/sum;
}

Double_t  AliTPCClusterParam::QmaxCorrection(Int_t sector, Int_t row, Float_t cpad, Float_t ctime, Float_t ky, Float_t kz, Float_t rmsy0, Float_t rmsz0, Float_t effPad, Float_t effDiff){
  //
  //
  // cpad      - pad (y) coordinate
  // ctime     - time(z) coordinate
  // ky        - dy/dx
  // kz        - dz/dx
  // rmsy0     - RF width in pad units
  // rmsz0     - RF width in time bin  units
  // effLength - contibution of PRF and diffusion
  // effDiff   - overwrite diffusion

  // Response function aproximated by convolution of gaussian with angular effect (integral=1)
  //  
  // Gaus width sy and sz is determined by RF width and diffusion 
  // Integral of Q is equal 1
  // Q max is calculated at position cpad, ctime
  // Example function:         
  //  TF1 f1("f1", "AliTPCClusterParam::QmaxCorrection(0,0.5,x,0,0,0.5,0.6)",0,1000) 
  //
  AliTPCParam * param   = AliTPCcalibDB::Instance()->GetParameters(); 
  Double_t padLength= param->GetPadPitchLength(sector,row);
  Double_t padWidth = param->GetPadPitchWidth(sector);
  Double_t zwidth   = param->GetZWidth();
  Double_t effLength= padLength+(param->GetWWPitch(0)+TMath::Sqrt(ctime*zwidth)*param->GetDiffT())*effPad;

  // diffusion in pad, time bin  units
  Double_t diffT=TMath::Sqrt(ctime*zwidth)*param->GetDiffT()/padWidth;
  Double_t diffL=TMath::Sqrt(ctime*zwidth)*param->GetDiffL()/zwidth;
  diffT*=effDiff;  //
  diffL*=effDiff;  //
  //
  // transform angular effect to pad units
  //
  Double_t pky   = ky*effLength/padWidth;
  Double_t pkz   = kz*effLength/zwidth;
  // position in pad unit
  Double_t py = (cpad+0.5)-TMath::Nint(cpad+0.5);
  Double_t pz = (ctime+0.5)-TMath::Nint(ctime+0.5);
  //
  //
  Double_t sy = TMath::Sqrt(rmsy0*rmsy0+diffT*diffT);
  Double_t sz = TMath::Sqrt(rmsz0*rmsz0+diffL*diffL); 
  //return GaussConvolutionGamma4(py,pz, pky,pkz,sy,sz,tau);
  Double_t length = padLength*TMath::Sqrt(1+ky*ky+kz*kz);
  return GaussConvolution(py,pz, pky,pkz,sy,sz)*length;
}

Double_t  AliTPCClusterParam::QtotCorrection(Int_t sector, Int_t row, Float_t cpad, Float_t ctime, Float_t ky, Float_t kz, Float_t rmsy0, Float_t rmsz0,  Float_t qtot, Float_t thr, Float_t effPad, Float_t effDiff){
  //
  //
  // cpad      - pad (y) coordinate
  // ctime     - time(z) coordinate
  // ky        - dy/dx
  // kz        - dz/dx
  // rmsy0     - RF width in pad units
  // rmsz0     - RF width in time bin  units
  // qtot      - the sum of signal in cluster - without thr correction
  // thr       - threshold
  // effLength - contibution of PRF and diffusion
  // effDiff   - overwrite diffusion

  // Response function aproximated by convolution of gaussian with angular effect (integral=1)
  //  
  // Gaus width sy and sz is determined by RF width and diffusion 
  // Integral of Q is equal 1
  // Q max is calculated at position cpad, ctime
  //          
  //  
  //
  AliTPCParam * param   = AliTPCcalibDB::Instance()->GetParameters(); 
  Double_t padLength= param->GetPadPitchLength(sector,row);
  Double_t padWidth = param->GetPadPitchWidth(sector);
  Double_t zwidth   = param->GetZWidth();
  Double_t effLength= padLength+(param->GetWWPitch(0)+TMath::Sqrt(ctime*zwidth)*param->GetDiffT())*effPad;
  //
  // diffusion in pad units
  Double_t diffT=TMath::Sqrt(ctime*zwidth)*param->GetDiffT()/padWidth;
  Double_t diffL=TMath::Sqrt(ctime*zwidth)*param->GetDiffL()/zwidth;
  diffT*=effDiff;  //
  diffL*=effDiff;  //
  //
  // transform angular effect to pad units 
  Double_t pky   = ky*effLength/padWidth;
  Double_t pkz   = kz*effLength/zwidth;
  // position in pad unit
  //  
  Double_t py = (cpad+0.5)-TMath::Nint(cpad+0.5);
  Double_t pz = (ctime+0.5)-TMath::Nint(ctime+0.5);
  //
  Double_t sy = TMath::Sqrt(rmsy0*rmsy0+diffT*diffT);
  Double_t sz = TMath::Sqrt(rmsz0*rmsz0+diffL*diffL); 
  //
  //
  //
  Double_t sumAll=0,sumThr=0;
  //
  Double_t corr =1;
  Double_t qnorm=qtot;
  for (Float_t iy=-3;iy<=3;iy+=1.)
    for (Float_t iz=-4;iz<=4;iz+=1.){
      //      Double_t val = GaussConvolutionGamma4(py-iy,pz-iz, pky,pkz, sy,sz,tau);      
      Double_t val = GaussConvolution(py-iy,pz-iz, pky,pkz, sy,sz);      
      Double_t qlocal =qnorm*val;
      if (TMath::Abs(iy)<1.5&&TMath::Abs(iz)<1.5){
      	sumThr+=qlocal;   // Virtual charge used in cluster finder
      }
      else{
	if (qlocal>thr && TMath::Abs(iz)<2.5&&TMath::Abs(iy)<2.5) sumThr+=qlocal;
      }
      sumAll+=qlocal;
    }
  if (sumAll>0&&sumThr>0) {
    corr=(sumThr)/sumAll;
  }
  //
  Double_t length = padLength*TMath::Sqrt(1+ky*ky+kz*kz);
  return corr*length;
}







