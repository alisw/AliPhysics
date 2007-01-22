#ifndef ALITPCCLUSTERPARAM_H
#define ALITPCCLUSTERPARAM_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliTPCClusterParam.h,v */

////////////////////////////////////////////////////
//                                                //
//  TPC cluster error and shape parameterization  //
//                                                //
////////////////////////////////////////////////////


#include <TObject.h>
class TTree;

//_____________________________________________________________________________
class AliTPCClusterParam : public TObject {
 public:
  static AliTPCClusterParam* Instance();
  AliTPCClusterParam(){fRatio=0.01;}
  virtual           ~AliTPCClusterParam(){;}
  virtual void	Print(Option_t* option = "") const;
  void SetInstance(AliTPCClusterParam*param){fgInstance = param;}
  void FitData(TTree * tree);
  Float_t GetError0(Int_t dim, Int_t type, Float_t z, Float_t angle);
  Float_t GetError0Par(Int_t dim, Int_t type, Float_t z, Float_t angle);
  Float_t GetError1(Int_t dim, Int_t type, Float_t z, Float_t angle);
  Float_t GetErrorQ(Int_t dim, Int_t type, Float_t z, Float_t angle, Float_t Qmean);
  Float_t GetErrorQPar(Int_t dim, Int_t type, Float_t z, Float_t angle, Float_t Qmean);
  Float_t GetErrorQParScaled(Int_t dim, Int_t type, Float_t z, Float_t angle, Float_t Qmean);
  Float_t GetRMS0(Int_t dim, Int_t type, Float_t z, Float_t angle);
  Float_t GetRMS1(Int_t dim, Int_t type, Float_t z, Float_t angle);
  Float_t GetRMSQ(Int_t dim, Int_t type, Float_t z, Float_t angle, Float_t Qmean);
  Float_t GetRMSSigma(Int_t dim, Int_t type, Float_t z, Float_t angle, Float_t Qmean);
  Float_t GetShapeFactor(Int_t dim, Int_t type, Float_t z, Float_t angle, Float_t Qmean, Float_t rmsL, Float_t rmsM);

  void Test(TTree * tree, const char *output="TestClusterParam.root");
  //
  // static methods equivalents  - use instance of param object - useful for tree draw and TF2 visualization 
  static Float_t SGetError0(Int_t dim, Int_t type, Float_t z, Float_t angle){
    return fgInstance->GetError0(dim,type,z,angle);
  }
  static Float_t SGetError0Par(Int_t dim, Int_t type, Float_t z, Float_t angle){
    return fgInstance->GetError0Par(dim,type,z,angle);
  }
  static Float_t SGetError1(Int_t dim, Int_t type, Float_t z, Float_t angle){
    return fgInstance->GetError1(dim,type,z,angle);
  }
  static Float_t SGetErrorQ(Int_t dim, Int_t type, Float_t z, Float_t angle, Float_t Qmean){
    return fgInstance->GetErrorQ(dim,type,z,angle,Qmean);
  }
  static Float_t SGetErrorQPar(Int_t dim, Int_t type, Float_t z, Float_t angle, Float_t Qmean){
    return fgInstance->GetErrorQPar(dim,type,z,angle,Qmean);
  }
  static Float_t SGetErrorQParScaled(Int_t dim, Int_t type, Float_t z, Float_t angle, Float_t Qmean){
    return fgInstance->GetErrorQParScaled(dim,type,z,angle,Qmean);
  }

  static Float_t SGetRMS0(Int_t dim, Int_t type, Float_t z, Float_t angle){
    return fgInstance->GetRMS0(dim,type,z,angle);
  }
  static Float_t SGetRMS1(Int_t dim, Int_t type, Float_t z, Float_t angle){
    return fgInstance->GetRMS1(dim,type,z,angle);
  }
  static Float_t SGetRMSQ(Int_t dim, Int_t type, Float_t z, Float_t angle, Float_t Qmean){
    return fgInstance->GetRMSQ(dim,type,z,angle,Qmean);
  }
  static Float_t SGetRMSSigma(Int_t dim, Int_t type, Float_t z, Float_t angle, Float_t Qmean){
    return fgInstance->GetRMSSigma(dim,type,z,angle,Qmean);
  }
  static Float_t SGetShapeFactor(Int_t dim, Int_t type, Float_t z, Float_t angle, Float_t Qmean, Float_t rmsL, Float_t rmsM){
    return fgInstance->GetShapeFactor(dim,type,z,angle,Qmean, rmsL, rmsM);
  }
  //
  //
  void FitResol(TTree * tree);
  void FitRMS(TTree * tree);
 protected: 
  void FitResol0(TTree * tree, Int_t dim, Int_t type, Float_t *param0, Float_t *error);
  void FitResol0Par(TTree * tree, Int_t dim, Int_t type, Float_t *param0, Float_t *error);
  void FitResol1(TTree * tree, Int_t dim, Float_t *param0, Float_t *error);
  void FitResolQ(TTree * tree, Int_t dim, Int_t type, Float_t *param0, Float_t *error);
  void FitResolQPar(TTree * tree, Int_t dim, Int_t type, Float_t *param0, Float_t *error);
  void FitRMS0(TTree * tree, Int_t dim, Int_t type, Float_t *param0, Float_t *error);
  void FitRMS1(TTree * tree, Int_t dim, Float_t *param0, Float_t *error);
  void FitRMSQ(TTree * tree, Int_t dim, Int_t type, Float_t *param0, Float_t *error);  
  void FitRMSSigma(TTree * tree, Int_t dim, Int_t type, Float_t *param0, Float_t *error);  
  //
  Float_t fRatio;               //ratio of values constibution to error
  Float_t fParamS0[2][3][4];    //error parameterization coeficients
  Float_t fErrorS0[2][3][4];    //error parameterization coeficients
  Float_t fParamS0Par[2][3][7];    //error parameterization coeficients
  Float_t fErrorS0Par[2][3][7];    //error parameterization coeficients  
  Float_t fParamSQ[2][3][6];    //error parameterization coeficients
  Float_t fErrorSQ[2][3][6];    //error parameterization coeficients
  Float_t fParamSQPar[2][3][9];    //error parameterization coeficients
  Float_t fErrorSQPar[2][3][9];    //error parameterization coeficients
  Float_t fParamS1[2][4];       //error parameterization coeficients
  Float_t fErrorS1[2][4];       //error parameterization coeficients
  //
  Float_t fParamRMS0[2][3][4];   //shape parameterization coeficients
  Float_t fErrorRMS0[2][3][4];   //shape parameterization coeficients
  Float_t fParamRMSQ[2][3][6];   //shape parameterization coeficients
  Float_t fErrorRMSQ[2][3][6];   //shape parameterization coeficients
  Float_t fParamRMS1[2][5];      //shape parameterization coeficients
  Float_t fErrorRMS1[2][5];      //shape parameterization coeficients
  Float_t fErrorRMSSys[2];        // systematic relative error of the parametererization
  Float_t fRMSSigmaRatio[2][2];   // mean value of the varation of RMS to RMS
  Float_t fRMSSigmaFit[2][3][2];   // mean value of the varation of RMS to RMS
  //
 protected:
  static AliTPCClusterParam*   fgInstance; //! Instance of this class (singleton implementation)
  ClassDef(AliTPCClusterParam,1)    //  TPC ROC class
};

#endif
