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
//   ESD track and V0 fast resolution parameterization
//   Fast algorithm to make parameterization
//   Track covariance used
//   
//             
//    Origin: Marian Ivanov marian.ivanov@cern.ch
//-------------------------------------------------------------------------

/*
  EXAMPLE USAGE:
  //
  // Make esd chain
  //
  .x ~/rootlogon.C
  gSystem->Load("libSTAT.so");
  gSystem->AddIncludePath("-I$ALICE_ROOT/TPC/macros");
  gROOT->LoadMacro("$ALICE_ROOT/TPC/macros/AliXRDPROOFtoolkit.cxx+");  
  AliXRDPROOFtoolkit tool;
  TChain * tree = tool.MakeChain("esd.txt","esdTree",0,1000);
  tree->Lookup();
  //
  // Load macros
  //
  .L $ALICE_ROOT/PWG1/AliESDresolParams.cxx+
  .L $ALICE_ROOT/PWG1/AliESDresolMakerFast.cxx+
  AliESDresolParams params; 
  TCut cutDCA="Tracks[].fCchi2<100&&abs(Tracks[].fP[4])<8&&abs(Tracks[].fP[3])<1&&sqrt(Tracks[].fC[0])/(0.2+abs(Tracks[].fP[4]))<0.02&&abs(Tracks[].fX)<3&&Tracks[].fITSncls>4&&Tracks.fTPCncls>40"
  //
  TObjArray * array = AliESDresolMakerFast::MakeParamPrimFast(tree,cutDCA,0.95,2000);


*/



#include "TVectorD.h"
#include "../STAT/TStatToolkit.h"
#include "TMath.h"
#include "TCut.h"
#include "TTree.h"

#include "AliESDresolParams.h"
#include "AliESDresolMakerFast.h"


ClassImp(AliESDresolMakerFast)




AliESDresolMakerFast::AliESDresolMakerFast() :
  TObject()
{
  //
  // Default constructor
  //
}


TObjArray * AliESDresolMakerFast::MakeParamPrimFast(TTree * tree, TCut & cutDCA, Float_t fraction, Int_t entries){
  //
  // DCA resolution linear parameterization
  // Only valid in ITS acceptance
  // Arguments:
  // tree     - esdTree or chain
  // cutDCA   - track selection criteria
  // fraction - robust fitter fraction
  // entries  - total number of entries (tracks) used

  /* Example
    
     TCut cutDCA="Tracks[].fCchi2<100&&abs(Tracks[].fP[4])<8&&abs(Tracks[].fP[3])<1&&sqrt(Tracks[].fC[0])/(0.2+abs(Tracks[].fP[4]))<0.02&&abs(Tracks[].fX)<3&&Tracks[].fITSncls>4&&Tracks.fTPCncls>40"
     fraction =0.95;
     entries  = 1000;
  */

  TObjArray *array = new TObjArray();
  Double_t chi2=0;
  Int_t    npoints=0;
  TVectorD fitParam;
  TMatrixD covMatrix;
  //
  // y param
  //
  TString * dcayParam=  TStatToolkit::FitPlane(tree,"sqrt(Tracks[].fC[0])/(0.2+abs(Tracks[].fP[4]))","(abs(Tracks[].fP[4]))++(abs(Tracks[].fP[4]))^2++(abs(Tracks[].fP[3]))++(abs(Tracks[].fP[3]))^2++(abs(Tracks[].fP[4]))*(abs(Tracks[].fP[3]))", cutDCA, chi2,npoints,fitParam,covMatrix,fraction, 0, entries);
  tree->SetAlias("dcayParam",dcayParam->Data());
  array->AddAt(new TVectorD(fitParam),0);
  printf("Y resol\t%s\n",dcayParam->Data());
  //
  // z param
  //
  TString * dcazParam= TStatToolkit::FitPlane(tree, "sqrt(Tracks[].fC[2])/(0.2+abs(Tracks[].fP[4]))", "(abs(Tracks[].fP[4]))++(abs(Tracks[].fP[4]))^2++(abs(Tracks[].fP[3]))++(abs(Tracks[].fP[3]))^2++(abs(Tracks[].fP[4]))*(abs(Tracks[].fP[3]))",  cutDCA, chi2,npoints,fitParam,covMatrix,fraction,0,entries);  
  printf("Z resol\t%s\n",dcazParam->Data());  
  tree->SetAlias("dcazParam",dcazParam->Data());
  array->AddAt(new TVectorD(fitParam),1);
  //
  // Phi param
  //
  TString * dcaphiParam= TStatToolkit::FitPlane(tree, "sqrt(Tracks[].fC[5])/(0.1+abs(Tracks[].fP[4]))", "(abs(Tracks[].fP[4]))++(abs(Tracks[].fP[4]))^2++(abs(Tracks[].fP[3]))++(abs(Tracks[].fP[3]))^2++(abs(Tracks[].fP[4]))*(abs(Tracks[].fP[3]))",  cutDCA, chi2,npoints,fitParam,covMatrix,fraction,0,entries);  
  printf("Phi resol\t%s\n",dcaphiParam->Data());   
  tree->SetAlias("dcaphiParam",dcaphiParam->Data());
  array->AddAt(new TVectorD(fitParam),2);
  //
  // theta param
  //
  TString * dcathParam=  TStatToolkit::FitPlane(tree, "sqrt(Tracks[].fC[9])/(0.1+abs(Tracks[].fP[4]))", "(abs(Tracks[].fP[4]))++(abs(Tracks[].fP[4]))^2++(abs(Tracks[].fP[3]))++(abs(Tracks[].fP[3]))^2++(abs(Tracks[].fP[4]))*(abs(Tracks[].fP[3]))",  cutDCA, chi2,npoints,fitParam,covMatrix,fraction,0,entries);  
  printf("Theta resol\t%s\n",dcathParam->Data());   
  tree->SetAlias("dcathParam",dcathParam->Data());
  array->AddAt(new TVectorD(fitParam),3);
  //
  // 1pt param
  //
  TString * dca1ptParam=  TStatToolkit::FitPlane(tree, "sqrt(Tracks[].fC[14])/(1+abs(Tracks[].fP[4]))^2", "(abs(Tracks[].fP[4]))++(abs(Tracks[].fP[4]))^2++(abs(Tracks[].fP[3]))++(abs(Tracks[].fP[3]))^2++(abs(Tracks[].fP[4]))*(abs(Tracks[].fP[3]))",  cutDCA, chi2,npoints,fitParam,covMatrix,fraction,0,entries);  
  printf("1pt resol\t%s\n",dca1ptParam->Data());  
  tree->SetAlias("dca1ptParam",dca1ptParam->Data());
  array->AddAt(new TVectorD(fitParam),4);
  return array;
}



TObjArray * AliESDresolMakerFast::MakeParamRFast(TTree * tree,  TCut &cutV0, Float_t fraction, Int_t entries){
  //
  // 
  //
  Double_t chi2=0;
  Int_t    npoints=0;
  TVectorD fitParam;
  TMatrixD covMatrix;
  //
  //
  /*
    TCut  cutV0 = "abs((V0s[].fParamP.fC[0]))<3&&abs(V0s[].fParamP.fP[3])<1&&abs(V0s[].fParamP.fP[4])<8";//
  */
  TObjArray *array = new TObjArray;
  //
  //
  //
  TString * v0sigmaY= TStatToolkit::FitPlane(tree,"sqrt(sqrt((V0s[].fParamP.fC[0])))","abs(V0s[].fParamP.fP[4])++V0s[].fParamP.fX++abs(V0s[].fParamP.fP[4])^2++V0s[].fParamP.fX^2++abs(V0s[].fParamP.fP[4])*V0s[].fParamP.fX++abs(V0s[].fParamP.fP[4])*V0s[].fParamP.fX^2", cutV0, chi2,npoints,fitParam,covMatrix,fraction,0,entries);
  tree->SetAlias("v0sigmaY",v0sigmaY->Data());
  array->AddLast(new TVectorD(fitParam));
  //
  //
  //
  TString * v0sigmaZ= TStatToolkit::FitPlane(tree,"sqrt(sqrt((V0s[].fParamP.fC[2])))","abs(V0s[].fParamP.fP[4])++V0s[].fParamP.fX++abs(V0s[].fParamP.fP[4])^2++V0s[].fParamP.fX^2++abs(V0s[].fParamP.fP[4])*V0s[].fParamP.fX++abs(V0s[].fParamP.fP[4])*V0s[].fParamP.fX^2", cutV0, chi2,npoints,fitParam,covMatrix,fraction,0,entries);
  tree->SetAlias("v0sigmaZ",v0sigmaZ->Data());
  array->AddLast(new TVectorD(fitParam));
  //
  //
  //
  TString * v0sigmaPhi= TStatToolkit::FitPlane(tree,"sqrt(sqrt((V0s[].fParamP.fC[5])))","abs(V0s[].fParamP.fP[4])++V0s[].fParamP.fX++abs(V0s[].fParamP.fP[4])^2++V0s[].fParamP.fX^2++abs(V0s[].fParamP.fP[4])*V0s[].fParamP.fX++abs(V0s[].fParamP.fP[4])*V0s[].fParamP.fX^2", cutV0, chi2,npoints,fitParam,covMatrix,fraction,0,entries);
  tree->SetAlias("v0sigmaPhi",v0sigmaPhi->Data());
  array->AddLast(new TVectorD(fitParam));
  //
  //
  //
  TString * v0sigmaTh= TStatToolkit::FitPlane(tree,"sqrt(sqrt((V0s[].fParamP.fC[9])))","abs(V0s[].fParamP.fP[4])++V0s[].fParamP.fX++abs(V0s[].fParamP.fP[4])^2++V0s[].fParamP.fX^2++abs(V0s[].fParamP.fP[4])*V0s[].fParamP.fX++abs(V0s[].fParamP.fP[4])*V0s[].fParamP.fX^2", cutV0, chi2,npoints,fitParam,covMatrix,fraction,0,entries);
  tree->SetAlias("v0sigmaTh",v0sigmaTh->Data());
  array->AddLast(new TVectorD(fitParam));
  //
  //
  //
  TString * v0sigma1pt= TStatToolkit::FitPlane(tree,"sqrt(sqrt((V0s[].fParamP.fC[14])))","abs(V0s[].fParamP.fP[4])++V0s[].fParamP.fX++abs(V0s[].fParamP.fP[4])^2++V0s[].fParamP.fX^2++abs(V0s[].fParamP.fP[4])*V0s[].fParamP.fX++abs(V0s[].fParamP.fP[4])*V0s[].fParamP.fX^2", cutV0, chi2,npoints,fitParam,covMatrix,fraction,0,entries);
  tree->SetAlias("v0sigma1pt",v0sigma1pt->Data());
  array->AddLast(new TVectorD(fitParam));
  return array;

}
