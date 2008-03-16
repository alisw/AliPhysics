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

//************************************************************************
//
//   ESD track and V0 resolution parameterization
//   
//             
//    Origin: Marian Ivanov marian.ivanov@cern.ch
//-------------------------------------------------------------------------




#include "TVectorD.h"
#include "../STAT/TStatToolkit.h"
#include "TMath.h"
#include "TCut.h"
#include "TTree.h"

#include "AliESDresolParams.h"


ClassImp(AliESDresolParams)


 AliESDresolParams*    AliESDresolParams::fgInstance = 0; //! Instance of this class (singleton implementation)


AliESDresolParams::AliESDresolParams() :
  TObject(),  
  fResolDCAyy(0),            // resolution Y parameterization
  fResolDCAzz(0),            // resolution Z parameterization   
  fResolDCAphi(0),           // resolution phi parameterization - pt-theta
  fResolDCAth(0),            // resolution theta parameterization -pt-theta
  fResolDCA1pt(0),           // resolution 1/pt parameterization - pt-theta
  //
  fResolCyy(0),              // DCA resolution Y parameterization - r-pt
  fResolCzz(0),              // DCA resolution Z parameterization - r-pt
  fResolCphi(0),             // DCA resolution phi parameterization - r-pt
  fResolCth(0),              // DCA resolution theta parameterization - r-pt
  fResolC1pt(0)             // DCA resolution 1/pt parameterization - r-pt
{
  //
  // Default constructor
  //
}

Double_t AliESDresolParams::GetResolPrim(Int_t param, Float_t onept, Float_t tanth){
  //
  // Resolution at primary vertex
  // simple Resolution parameterization 
  // polynom of second order in 2D 
  // 
  if (!fResolDCAyy) return 0;
  TVectorD * pvec     = fResolDCAyy;
  if (param==1)  pvec = fResolDCAzz;
  if (param==2)  pvec = fResolDCAphi;
  if (param==3)  pvec = fResolDCAth;
  if (param==4)  pvec = fResolDCA1pt;
  TVectorD &vec = *pvec;
  //
  Float_t val = vec[0];
  val+= vec[1]*TMath::Abs(onept);
  val+= vec[2]*TMath::Abs(onept*onept);
  val+= vec[3]*TMath::Abs(tanth);
  val+= vec[4]*TMath::Abs(tanth*tanth);
  val+= vec[5]*TMath::Abs(onept*tanth);
  val*=val;
  return val;
}

Double_t AliESDresolParams::GetResolR(Int_t param, Float_t onept, Float_t radius){
  //
  // simple DCA resolution parameterization
  // polynom of second order in 2D 
  // 
  if (!fResolCyy) return 0;
  TVectorD * pvec     = fResolCyy;
  if (param==1)  pvec = fResolCzz;
  if (param==2)  pvec = fResolCphi;
  if (param==3)  pvec = fResolCth;
  if (param==4)  pvec = fResolC1pt;
  TVectorD &vec = *pvec;
  //
  Float_t val = vec[0];
  val+= vec[1]*TMath::Abs(onept);
  val+= vec[2]*TMath::Abs(radius);
  val+= vec[3]*TMath::Abs(onept*onept);
  val+= vec[4]*TMath::Abs(radius*radius);
  val+= vec[5]*TMath::Abs(radius*onept);
  val+= vec[6]*TMath::Abs(radius*radius*onept);
  val*=val;
  return val;
}



void AliESDresolParams::MakeParamPrimFast(TTree * tree,Float_t fraction, Int_t entries){
  //
  // DCA resolution linear parameterization
  // Only valid in ITS acceptance
  //
  Double_t chi2=0;
  Int_t    npoints=0;
  TVectorD fitParam;
  TMatrixD covMatrix;
  TCut  cutDCA = "sqrt(Tracks[].fC[0])<1&&Tracks[].fCchi2<100&&abs(Tracks[].fP[4])<8&&abs(Tracks[].fP[3])<1";  
  //
  // y param
  //
  TString * dcayParam= 
    TStatToolkit::FitPlane(tree,"sqrt(sqrt(Tracks[].fC[0]))","(abs(Tracks[].fP[4]))++(abs(Tracks[].fP[4]))^2++(abs(Tracks[].fP[3]))++(abs(Tracks[].fP[3]))^2++(abs(Tracks[].fP[4]))*(abs(Tracks[].fP[3]))", cutDCA, chi2,npoints,fitParam,covMatrix,fraction, 0, entries);  
  fResolDCAyy = new TVectorD(fitParam);
  printf("Y resol\t%s\n",dcayParam->Data());
  //
  // z param
  //
  TString * dcazParam= 
    TStatToolkit::FitPlane(tree, "sqrt(sqrt(Tracks[].fC[2]))", "(abs(Tracks[].fP[4]))++(abs(Tracks[].fP[4]))^2++(abs(Tracks[].fP[3]))++(abs(Tracks[].fP[3]))^2++(abs(Tracks[].fP[4]))*(abs(Tracks[].fP[3]))",  cutDCA, chi2,npoints,fitParam,covMatrix,fraction,0,entries);  
  printf("Z resol\t%s\n",dcazParam->Data());  
  fResolDCAzz = new TVectorD(fitParam);
  //
  // Phi param
  //
  TString * dcaphiParam= 
    TStatToolkit::FitPlane(tree, "sqrt(sqrt(Tracks[].fC[5]))", "(abs(Tracks[].fP[4]))++(abs(Tracks[].fP[4]))^2++(abs(Tracks[].fP[3]))++(abs(Tracks[].fP[3]))^2++(abs(Tracks[].fP[4]))*(abs(Tracks[].fP[3]))",  cutDCA, chi2,npoints,fitParam,covMatrix,fraction,0,entries);  
  printf("Phi resol\t%s\n",dcaphiParam->Data());  
  fResolDCAphi = new TVectorD(fitParam);
  //
  // theta param
  //
  TString * dcathParam= 
    TStatToolkit::FitPlane(tree, "sqrt(sqrt(Tracks[].fC[9]))", "(abs(Tracks[].fP[4]))++(abs(Tracks[].fP[4]))^2++(abs(Tracks[].fP[3]))++(abs(Tracks[].fP[3]))^2++(abs(Tracks[].fP[4]))*(abs(Tracks[].fP[3]))",  cutDCA, chi2,npoints,fitParam,covMatrix,fraction,0,entries);  
  printf("Theta resol\t%s\n",dcathParam->Data());  
  fResolDCAth = new TVectorD(fitParam);
  //
  // 1pt param
  //
  TString * dca1ptParam= 
    TStatToolkit::FitPlane(tree, "sqrt(sqrt(Tracks[].fC[14]))", "(abs(Tracks[].fP[4]))++(abs(Tracks[].fP[4]))^2++(abs(Tracks[].fP[3]))++(abs(Tracks[].fP[3]))^2++(abs(Tracks[].fP[4]))*(abs(Tracks[].fP[3]))",  cutDCA, chi2,npoints,fitParam,covMatrix,fraction,0,entries);  
  printf("1pt resol\t%s\n",dca1ptParam->Data());  
  fResolDCA1pt = new TVectorD(fitParam);

}



void AliESDresolParams::MakeParamRFast(TTree * tree, Float_t fraction, Int_t entries){
  //
  // 
  //
  Double_t chi2=0;
  Int_t    npoints=0;
  TVectorD fitParam;
  TMatrixD covMatrix;
  //  fraction=0.8;
  //entries=100000000;
  //
  //
  TCut  cutV0 = "abs((V0s[].fParamP.fC[0]))<3&&abs(V0s[].fParamP.fP[3])<1&&abs(V0s[].fParamP.fP[4])<8";//
  //
  //
  //
  TString * v0sigmaY= TStatToolkit::FitPlane(tree,"sqrt(sqrt((V0s[].fParamP.fC[0])))","abs(V0s[].fParamP.fP[4])++V0s[].fParamP.fX++abs(V0s[].fParamP.fP[4])^2++V0s[].fParamP.fX^2++abs(V0s[].fParamP.fP[4])*V0s[].fParamP.fX++abs(V0s[].fParamP.fP[4])*V0s[].fParamP.fX^2", cutV0, chi2,npoints,fitParam,covMatrix,fraction,0,entries);
  tree->SetAlias("v0sigmaY",v0sigmaY->Data());
  fResolCyy = new TVectorD(fitParam);
  //
  //
  //
  TString * v0sigmaZ= TStatToolkit::FitPlane(tree,"sqrt(sqrt((V0s[].fParamP.fC[2])))","abs(V0s[].fParamP.fP[4])++V0s[].fParamP.fX++abs(V0s[].fParamP.fP[4])^2++V0s[].fParamP.fX^2++abs(V0s[].fParamP.fP[4])*V0s[].fParamP.fX++abs(V0s[].fParamP.fP[4])*V0s[].fParamP.fX^2", cutV0, chi2,npoints,fitParam,covMatrix,fraction,0,entries);
  tree->SetAlias("v0sigmaZ",v0sigmaZ->Data());
  fResolCzz = new TVectorD(fitParam);
  //
  //
  //
  TString * v0sigmaPhi= TStatToolkit::FitPlane(tree,"sqrt(sqrt((V0s[].fParamP.fC[5])))","abs(V0s[].fParamP.fP[4])++V0s[].fParamP.fX++abs(V0s[].fParamP.fP[4])^2++V0s[].fParamP.fX^2++abs(V0s[].fParamP.fP[4])*V0s[].fParamP.fX++abs(V0s[].fParamP.fP[4])*V0s[].fParamP.fX^2", cutV0, chi2,npoints,fitParam,covMatrix,fraction,0,entries);
  tree->SetAlias("v0sigmaPhi",v0sigmaPhi->Data());
  fResolCphi = new TVectorD(fitParam);
  //
  //
  //
  TString * v0sigmaTh= TStatToolkit::FitPlane(tree,"sqrt(sqrt((V0s[].fParamP.fC[9])))","abs(V0s[].fParamP.fP[4])++V0s[].fParamP.fX++abs(V0s[].fParamP.fP[4])^2++V0s[].fParamP.fX^2++abs(V0s[].fParamP.fP[4])*V0s[].fParamP.fX++abs(V0s[].fParamP.fP[4])*V0s[].fParamP.fX^2", cutV0, chi2,npoints,fitParam,covMatrix,fraction,0,entries);
  tree->SetAlias("v0sigmaTh",v0sigmaTh->Data());
  fResolCth = new TVectorD(fitParam);
  //
  //
  //
  TString * v0sigma1pt= TStatToolkit::FitPlane(tree,"sqrt(sqrt((V0s[].fParamP.fC[14])))","abs(V0s[].fParamP.fP[4])++V0s[].fParamP.fX++abs(V0s[].fParamP.fP[4])^2++V0s[].fParamP.fX^2++abs(V0s[].fParamP.fP[4])*V0s[].fParamP.fX++abs(V0s[].fParamP.fP[4])*V0s[].fParamP.fX^2", cutV0, chi2,npoints,fitParam,covMatrix,fraction,0,entries);
  tree->SetAlias("v0sigma1pt",v0sigma1pt->Data());
  fResolC1pt = new TVectorD(fitParam);


}
