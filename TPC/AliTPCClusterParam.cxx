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
//  TPC cluster error and shape parameterization                             //
//                                                                           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////
#include "AliTPCClusterParam.h"
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include <TVectorF.h>
#include <TLinearFitter.h>
#include <TH1F.h>
#include <TProfile2D.h>

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




void AliTPCClusterParam::FitResol0(TTree * tree, Int_t dim, Int_t type, Float_t *param0, Float_t *error){
  //
  // Fit z - angular dependence of resolution 
  //
  // Int_t dim=0, type=0;
  char varVal[100];
  sprintf(varVal,"Resol:AngleM:Zm");
  char varErr[100];
  sprintf(varErr,"Sigma:AngleS:Zs");
  char varCut[100];
  sprintf(varCut,"Dim==%d&&Pad==%d&&QMean<0",dim,type);
  //
  Int_t entries = tree->Draw(varVal,varCut);
  Float_t px[10000], py[10000], pz[10000];
  Float_t ex[10000], ey[10000], ez[10000];
  //
  tree->Draw(varErr,varCut);  
  for (Int_t ipoint=0; ipoint<entries; ipoint++){
    ex[ipoint]= tree->GetV3()[ipoint];
    ey[ipoint]= tree->GetV2()[ipoint];
    ez[ipoint]= tree->GetV1()[ipoint];
  } 
  tree->Draw(varVal,varCut);
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
  char varVal[100];
  sprintf(varVal,"Resol:AngleM:Zm");
  char varErr[100];
  sprintf(varErr,"Sigma:AngleS:Zs");
  char varCut[100];
  sprintf(varCut,"Dim==%d&&Pad==%d&&QMean<0",dim,type);
  //
  Int_t entries = tree->Draw(varVal,varCut);
  Float_t px[10000], py[10000], pz[10000];
  Float_t ex[10000], ey[10000], ez[10000];
  //
  tree->Draw(varErr,varCut);  
  for (Int_t ipoint=0; ipoint<entries; ipoint++){
    ex[ipoint]= tree->GetV3()[ipoint];
    ey[ipoint]= tree->GetV2()[ipoint];
    ez[ipoint]= tree->GetV1()[ipoint];
  } 
  tree->Draw(varVal,varCut);
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
  char varVal[100];
  sprintf(varVal,"Resol:AngleM*sqrt(Length):Zm/Length");
  char varErr[100];
  sprintf(varErr,"Sigma:AngleS:Zs");
  char varCut[100];
  sprintf(varCut,"Dim==%d&&QMean<0",dim);
  //
  Int_t entries = tree->Draw(varVal,varCut);
  Float_t px[10000], py[10000], pz[10000];
  Float_t ex[10000], ey[10000], ez[10000];
  //
  tree->Draw(varErr,varCut);  
  for (Int_t ipoint=0; ipoint<entries; ipoint++){
    ex[ipoint]= tree->GetV3()[ipoint];
    ey[ipoint]= tree->GetV2()[ipoint];
    ez[ipoint]= tree->GetV1()[ipoint];
  } 
  tree->Draw(varVal,varCut);
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
  char varVal[100];
  sprintf(varVal,"Resol:AngleM/sqrt(QMean):Zm/QMean");
  char varVal0[100];
  sprintf(varVal0,"Resol:AngleM:Zm");
  //
  char varErr[100];
  sprintf(varErr,"Sigma:AngleS:Zs");
  char varCut[100];
  sprintf(varCut,"Dim==%d&&Pad==%d&&QMean>0",dim,type);
  //
  Int_t entries = tree->Draw(varVal,varCut);
  Float_t px[20000], py[20000], pz[20000], pu[20000], pt[20000];
  Float_t ex[20000], ey[20000], ez[20000];
  //
  tree->Draw(varErr,varCut);  
  for (Int_t ipoint=0; ipoint<entries; ipoint++){
    ex[ipoint]= tree->GetV3()[ipoint];
    ey[ipoint]= tree->GetV2()[ipoint];
    ez[ipoint]= tree->GetV1()[ipoint];
  } 
  tree->Draw(varVal,varCut);
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
  char varVal[100];
  sprintf(varVal,"Resol:AngleM/sqrt(QMean):Zm/QMean");
  char varVal0[100];
  sprintf(varVal0,"Resol:AngleM:Zm");
  //
  char varErr[100];
  sprintf(varErr,"Sigma:AngleS:Zs");
  char varCut[100];
  sprintf(varCut,"Dim==%d&&Pad==%d&&QMean>0",dim,type);
  //
  Int_t entries = tree->Draw(varVal,varCut);
  Float_t px[20000], py[20000], pz[20000], pu[20000], pt[20000];
  Float_t ex[20000], ey[20000], ez[20000];
  //
  tree->Draw(varErr,varCut);  
  for (Int_t ipoint=0; ipoint<entries; ipoint++){
    ex[ipoint]= tree->GetV3()[ipoint];
    ey[ipoint]= tree->GetV2()[ipoint];
    ez[ipoint]= tree->GetV1()[ipoint];
  } 
  tree->Draw(varVal,varCut);
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
  char varVal[100];
  sprintf(varVal,"RMSm:AngleM:Zm");
  char varErr[100];
  sprintf(varErr,"sqrt((1./(100.*sqrt(12.))^2)+RMSe0^2):AngleS:Zs");
  char varCut[100];
  sprintf(varCut,"Dim==%d&&Pad==%d&&QMean<0",dim,type);
  //
  Int_t entries = tree->Draw(varVal,varCut);
  Float_t px[10000], py[10000], pz[10000];
  Float_t ex[10000], ey[10000], ez[10000];
  //
  tree->Draw(varErr,varCut);  
  for (Int_t ipoint=0; ipoint<entries; ipoint++){
    ex[ipoint]= tree->GetV3()[ipoint];
    ey[ipoint]= tree->GetV2()[ipoint];
    ez[ipoint]= tree->GetV1()[ipoint];
  } 
  tree->Draw(varVal,varCut);
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
  char varVal[100];
  sprintf(varVal,"RMSm:AngleM*Length:Zm");
  char varErr[100];
  sprintf(varErr,"sqrt((1./(100.*sqrt(12.))^2)+RMSe0^2):AngleS:Pad");
  char varCut[100];
  sprintf(varCut,"Dim==%d&&QMean<0",dim);
  //
  Int_t entries = tree->Draw(varVal,varCut);
  Float_t px[10000], py[10000], pz[10000];
  Float_t type[10000], ey[10000], ez[10000];
  //
  tree->Draw(varErr,varCut);  
  for (Int_t ipoint=0; ipoint<entries; ipoint++){
    type[ipoint] = tree->GetV3()[ipoint];
    ey[ipoint]   = tree->GetV2()[ipoint];
    ez[ipoint]   = tree->GetV1()[ipoint];
  } 
  tree->Draw(varVal,varCut);
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
  char varVal[100];
  sprintf(varVal,"RMSm:AngleM/sqrt(QMean):Zm/QMean");
  char varVal0[100];
  sprintf(varVal0,"RMSm:AngleM:Zm");
  //
  char varErr[100];
  sprintf(varErr,"sqrt((1./(100.*sqrt(12.))^2)+RMSe0^2):AngleS:Zs");
  char varCut[100];
  sprintf(varCut,"Dim==%d&&Pad==%d&&QMean>0",dim,type);
  //
  Int_t entries = tree->Draw(varVal,varCut);
  Float_t px[20000], py[20000], pz[20000], pu[20000], pt[20000];
  Float_t ex[20000], ey[20000], ez[20000];
  //
  tree->Draw(varErr,varCut);  
  for (Int_t ipoint=0; ipoint<entries; ipoint++){
    ex[ipoint]= tree->GetV3()[ipoint];
    ey[ipoint]= tree->GetV2()[ipoint];
    ez[ipoint]= tree->GetV1()[ipoint];
  } 
  tree->Draw(varVal,varCut);
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
  char varVal[100];
  sprintf(varVal,"RMSs:RMSm");
  //
  char varCut[100];
  sprintf(varCut,"Dim==%d&&Pad==%d&&QMean<0",dim,type);
  //
  Int_t entries = tree->Draw(varVal,varCut);
  Float_t px[20000], py[20000];
  //
  tree->Draw(varVal,varCut);
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
    fitter.AddPoint(x,val,err);
  }
  fitter.Eval();
  param0[0]= fitter.GetParameter(0);
  param0[1]= fitter.GetParameter(1);
}



Float_t  AliTPCClusterParam::GetError0(Int_t dim, Int_t type, Float_t z, Float_t angle){
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


Float_t  AliTPCClusterParam::GetError0Par(Int_t dim, Int_t type, Float_t z, Float_t angle){
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



Float_t  AliTPCClusterParam::GetError1(Int_t dim, Int_t type, Float_t z, Float_t angle){
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

Float_t  AliTPCClusterParam::GetErrorQ(Int_t dim, Int_t type, Float_t z, Float_t angle, Float_t Qmean){
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

Float_t  AliTPCClusterParam::GetErrorQPar(Int_t dim, Int_t type, Float_t z, Float_t angle, Float_t Qmean){
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

Float_t  AliTPCClusterParam::GetErrorQParScaled(Int_t dim, Int_t type, Float_t z, Float_t angle, Float_t Qmean){
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

Float_t  AliTPCClusterParam::GetRMS0(Int_t dim, Int_t type, Float_t z, Float_t angle){
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

Float_t  AliTPCClusterParam::GetRMS1(Int_t dim, Int_t type, Float_t z, Float_t angle){
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

Float_t  AliTPCClusterParam::GetRMSQ(Int_t dim, Int_t type, Float_t z, Float_t angle, Float_t Qmean){
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

Float_t  AliTPCClusterParam::GetRMSSigma(Int_t dim, Int_t type, Float_t z, Float_t angle, Float_t Qmean){
  //
  // calculates RMS of signal shape fluctuation
  //
  Float_t mean = GetRMSQ(dim,type,z,angle,Qmean);
  Float_t value  = fRMSSigmaFit[dim][type][0];
  value+=  fRMSSigmaFit[dim][type][1]*mean;
  return value;
}

Float_t  AliTPCClusterParam::GetShapeFactor(Int_t dim, Int_t type, Float_t z, Float_t angle, Float_t Qmean, Float_t rmsL, Float_t rmsM){
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
  return TMath::Sqrt(value);
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
      sprintf(hname1,"Delta0 Dir %d Pad %d",idim,ipad);
      sprintf(hcut1,"Dim==%d&&QMean<0&&Pad==%d",idim,ipad);
      sprintf(hexp1,"(Resol-AliTPCClusterParam::SGetError0(Dim,Pad,Zm,AngleM))/Resol>>%s",hname1);
      TH1F  his1DRel0(hname1, hname1, 100,-0.2, 0.2);
      sprintf(hname1,"Dim==%d&&QMean<0&&Pad=%d",idim,ipad);
      tree->Draw(hexp1,hcut1,"");
      his1DRel0.Write();
      //
      sprintf(hname1,"Delta0Par Dir %d Pad %d",idim,ipad);
      sprintf(hcut1,"Dim==%d&&QMean<0&&Pad==%d",idim,ipad);
      sprintf(hexp1,"(Resol-AliTPCClusterParam::SGetError0Par(Dim,Pad,Zm,AngleM))/Resol>>%s",hname1);
      TH1F  his1DRel0Par(hname1, hname1, 100,-0.2, 0.2);
      sprintf(hname1,"Dim==%d&&QMean<0&&Pad=%d",idim,ipad);
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
      sprintf(hname1,"2DDelta0 Dir %d Pad %d",idim,ipad);
      sprintf(hcut1,"Dim==%d&&QMean<0&&Pad==%d",idim,ipad);
      sprintf(hexp1,"(Resol-AliTPCClusterParam::SGetError0(Dim,Pad,Zm,AngleM))/Resol:AngleM:Zm>>%s",hname1);
      TProfile2D  profDRel0(hname1, hname1, 6,0,250,6,0,1);
      sprintf(hname1,"Dim==%d&&QMean<0&&Pad=%d",idim,ipad);
      tree->Draw(hexp1,hcut1,"");
      profDRel0.Write();
      //
      sprintf(hname1,"2DDelta0Par Dir %d Pad %d",idim,ipad);
      sprintf(hcut1,"Dim==%d&&QMean<0&&Pad==%d",idim,ipad);
      sprintf(hexp1,"(Resol-AliTPCClusterParam::SGetError0Par(Dim,Pad,Zm,AngleM))/Resol:AngleM:Zm>>%s",hname1);
      TProfile2D profDRel0Par(hname1, hname1,6,0,250,6,0,1);
      sprintf(hname1,"Dim==%d&&QMean<0&&Pad=%d",idim,ipad);
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
  printf("Y\tall\t%f\t%f\t%f\t%f\n", TMath::Sqrt(TMath::Abs(fParamS1[0][0])),TMath::Sqrt(fParamS1[0][1]),
	 TMath::Sqrt(fParamS1[0][2]),TMath::Sqrt(fParamS1[0][3]));
  for (Int_t ipad=0; ipad<3; ipad++){
    Float_t length=0.75;
    if (ipad==1) length=1;
    if (ipad==2) length=1.5;    
    printf("\t%d\t%f\t%f\t%f\t%f\n", ipad, 
	   TMath::Sqrt(TMath::Abs(fParamS0[0][ipad][0])),
	   TMath::Sqrt(fParamS0[0][ipad][1]*length),
	   TMath::Sqrt(fParamS0[0][ipad][2]/length),
	   TMath::Sqrt(fParamS0[0][ipad][3]));
  }
  for (Int_t ipad=0; ipad<3; ipad++){
    Float_t length=0.75;
    if (ipad==1) length=1;
    if (ipad==2) length=1.5;
    printf("\t%dPar\t%f\t%f\t%f\t%f\n", ipad, 
	   TMath::Sqrt(TMath::Abs(fParamS0Par[0][ipad][0])),
	   TMath::Sqrt(fParamS0Par[0][ipad][1]*length),
	   TMath::Sqrt(fParamS0Par[0][ipad][2]/length),
	   TMath::Sqrt(fParamS0Par[0][ipad][6]));
  }
  printf("Z\tall\t%f\t%f\t%f\t%f\n", TMath::Sqrt(TMath::Abs(fParamS1[1][0])),TMath::Sqrt(fParamS1[1][1]),
	 TMath::Sqrt(fParamS1[1][2]), TMath::Sqrt(fParamS1[1][3]));
  
  for (Int_t ipad=0; ipad<3; ipad++){
    Float_t length=0.75;
    if (ipad==1) length=1;
    if (ipad==2) length=1.5;    
    printf("\t%d\t%f\t%f\t%f\t%f\n", ipad, 
	   TMath::Sqrt(TMath::Abs(fParamS0[1][ipad][0])),
	   TMath::Sqrt(fParamS0[1][ipad][1]*length),
	   TMath::Sqrt(fParamS0[1][ipad][2]/length),
	   TMath::Sqrt(fParamS0[1][ipad][3]));
  }
  for (Int_t ipad=0; ipad<3; ipad++){
    Float_t length=0.75;
    if (ipad==1) length=1;
    if (ipad==2) length=1.5;        
    printf("\t%dPar\t%f\t%f\t%f\t%f\n", ipad, 
	   TMath::Sqrt(TMath::Abs(fParamS0Par[1][ipad][0])),
	   TMath::Sqrt(fParamS0Par[1][ipad][1]*length),
	   TMath::Sqrt(fParamS0Par[1][ipad][2]/length),
	   TMath::Sqrt(fParamS0Par[1][ipad][6]));
  }
  
  //
  // RMS scaling
  //
  printf("\n");
  printf("\nRMS Scaled factors\n");
  printf("Dir\tPad\tP00\t\tP01\t\tP1\t\tP2\t\tchi2\n");
  printf("Y\tall\t%f\t%f\t%f\t%f\t%f\n", TMath::Sqrt(TMath::Abs(fParamRMS1[0][0])),TMath::Sqrt(fParamRMS1[0][1]),
	 TMath::Sqrt(fParamRMS1[0][2]),TMath::Sqrt(fParamRMS1[0][3]),TMath::Sqrt(fParamRMS1[0][4]));
  for (Int_t ipad=0; ipad<3; ipad++){
    Float_t length=0.75;
    if (ipad==1) length=1;
    if (ipad==2) length=1.5;    
    if (ipad==0){
      printf("\t%d\t%f\t%f\t%f\t%f\t%f\n", ipad, 
	     TMath::Sqrt(TMath::Abs(fParamRMS0[0][ipad][0])),
	     0.,
	     TMath::Sqrt(fParamRMS0[0][ipad][1]),
	     TMath::Sqrt(fParamRMS0[0][ipad][2]/(length*length)),
	     TMath::Sqrt(fParamRMS0[0][ipad][3]));
    }else{
      printf("\t%d\t%f\t%f\t%f\t%f\t%f\n", ipad, 
	     0.,
	     TMath::Sqrt(TMath::Abs(fParamRMS0[0][ipad][0])),
	     TMath::Sqrt(fParamRMS0[0][ipad][1]),
	     TMath::Sqrt(fParamRMS0[0][ipad][2]/(length*length)),
	     TMath::Sqrt(fParamRMS0[0][ipad][3]));	
    }
  }
  printf("\n");
  printf("Z\tall\t%f\t%f\t%f\t%f\t%f\n", TMath::Sqrt(TMath::Abs(fParamRMS1[1][0])),TMath::Sqrt(fParamRMS1[1][1]),
	 TMath::Sqrt(fParamRMS1[1][2]),TMath::Sqrt(fParamRMS1[1][3]),TMath::Sqrt(fParamRMS1[1][4]));
  for (Int_t ipad=0; ipad<3; ipad++){
    Float_t length=0.75;
    if (ipad==1) length=1;
    if (ipad==2) length=1.5;    
    if (ipad==0){
      printf("\t%d\t%f\t%f\t%f\t%f\t%f\n", ipad, 
	     TMath::Sqrt(TMath::Abs(fParamRMS0[1][ipad][0])),
	     0.,
	     TMath::Sqrt(fParamRMS0[1][ipad][1]),
	     TMath::Sqrt(fParamRMS0[1][ipad][2]/(length*length)),
	     TMath::Sqrt(fParamRMS0[1][ipad][3]));
    }else{
      printf("\t%d\t%f\t%f\t%f\t%f\t%f\n", ipad, 
	     0.,
	     TMath::Sqrt(TMath::Abs(fParamRMS0[1][ipad][0])),
	     TMath::Sqrt(fParamRMS0[1][ipad][1]),
	     TMath::Sqrt(fParamRMS0[1][ipad][2]/(length*length)),
	     TMath::Sqrt(fParamRMS0[1][ipad][3]));	
    }
  }
}





