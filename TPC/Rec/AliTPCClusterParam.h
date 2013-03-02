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
#include <TVectorDfwd.h>
#include <TMatrixDfwd.h>

class TTree;
class TObjArray;
class TH1;
class THnBase;
//_____________________________________________________________________________
class AliTPCClusterParam : public TObject {
 public:
  static AliTPCClusterParam* Instance();
  AliTPCClusterParam();
  AliTPCClusterParam(const AliTPCClusterParam& param);
  AliTPCClusterParam & operator=(const AliTPCClusterParam& param);
  virtual           ~AliTPCClusterParam();
  virtual void	Print(Option_t* option = "") const;
  void SetInstance(AliTPCClusterParam *const param){fgInstance = param;}
  //
  // Seting functions
  //
  void FitData(TTree * tree);
  void FitResol(TTree * tree);
  void FitRMS(TTree * tree);
  void SetQnorm(Int_t ipad, Int_t itype,  const TVectorD *const norm); 
  void SetQnormCorr(Int_t ipad, Int_t itype, Int_t corrType, Float_t val); 
  Double_t  GetQnormCorr(Int_t ipad, Int_t itype, Int_t corrType) const;
  TMatrixD *GetQnormCorrMatrix(){return fQNormCorr;};
  void ResetQnormCorr(); 
  void SetWaveCorrectionMap( THnBase *WaveCorrectionMap);
  void SetResolutionYMap( THnBase *ResolutionYMap);
  //
  // Charge parameterization
  //
  Float_t Qnorm(Int_t ipad, Int_t itype, Float_t dr, Float_t ty, Float_t tz);   
  Float_t QnormHis(Int_t ipad, Int_t itype, Float_t dr, Float_t ty, Float_t tz);


  Float_t QnormPos(Int_t ipad, Bool_t isMax,  Float_t pad, Float_t time, Float_t z, Float_t sy2, Float_t sz2, Float_t qm, Float_t qt);
  static Float_t SQnormPos(Int_t ipad, Bool_t isMax,  Float_t pad, Float_t time, Float_t z, Float_t sy2, Float_t sz2, Float_t qm, Float_t qt){ return fgInstance->QnormPos(ipad,isMax,pad,time,z,sy2,sz2,qm,qt);;}
 
  Float_t PosCorrection(Int_t type, Int_t ipad,  Float_t pad, Float_t time, Float_t z, Float_t sy2, Float_t sz2, Float_t qm);
  static Float_t  SPosCorrection(Int_t type, Int_t ipad,  Float_t pad, Float_t time, Float_t z, Float_t sy2, Float_t sz2, Float_t qm){ return fgInstance->PosCorrection(type,ipad,pad,time,z,sy2,sz2,qm);}
  //
  // Error parameterization
  //
  Float_t GetError0(Int_t dim, Int_t type, Float_t z, Float_t angle) const;
  Float_t GetError0Par(Int_t dim, Int_t type, Float_t z, Float_t angle) const;
  Float_t GetError1(Int_t dim, Int_t type, Float_t z, Float_t angle) const;
  Float_t GetErrorQ(Int_t dim, Int_t type, Float_t z, Float_t angle, Float_t Qmean) const;
  Float_t GetErrorQPar(Int_t dim, Int_t type, Float_t z, Float_t angle, Float_t Qmean) const;
  Float_t GetErrorQParScaled(Int_t dim, Int_t type, Float_t z, Float_t angle, Float_t Qmean) const;
  //
  // Shape parameterization
  //
  Float_t GetRMS0(Int_t dim, Int_t type, Float_t z, Float_t angle) const;
  Float_t GetRMS1(Int_t dim, Int_t type, Float_t z, Float_t angle) const;
  Float_t GetRMSQ(Int_t dim, Int_t type, Float_t z, Float_t angle, Float_t Qmean) const;
  Float_t GetRMSSigma(Int_t dim, Int_t type, Float_t z, Float_t angle, Float_t Qmean) const;
  Float_t GetShapeFactor(Int_t dim, Int_t type, Float_t z, Float_t angle, Float_t Qmean, Float_t rmsL, Float_t rmsM) const;
  // 
  // Correction and resolution maps
  //
  const THnBase *GetWaveCorrectionMap() const { return fWaveCorrectionMap; }
  const THnBase *GetResolutionYMap() const { return  fResolutionYMap; }
  Float_t GetWaveCorrection(Int_t Type, Float_t Z, Int_t QMax, Float_t Pad, Float_t angleY ) const;
  static Float_t SGetWaveCorrection(Int_t Type, Float_t Z, Int_t QMax, Float_t Pad, Float_t angleY ){return (fgInstance) ? fgInstance->GetWaveCorrection(Type,Z,QMax,Pad,angleY):0;}
  //
  //
  //
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
  static Float_t SQnorm(Int_t ipad, Int_t itype,Float_t dr, Float_t ty, Float_t tz) {return fgInstance->Qnorm(ipad, itype, dr,ty,tz);}
  static Float_t SQnormHis(Int_t ipad, Int_t itype,Float_t dr, Float_t ty, Float_t tz) {return fgInstance->QnormHis(ipad, itype, dr,ty,tz);}

  //
  // Analytical position angular correction
  //
  static Double_t  GaussConvolution(Double_t x0, Double_t x1, Double_t k0, Double_t k1, Double_t s0, Double_t s1);
  static Double_t  GaussConvolutionTail(Double_t x0, Double_t x1, Double_t k0, Double_t k1, Double_t s0, Double_t s1, Double_t tau);
  static Double_t  GaussConvolutionGamma4(Double_t x0, Double_t x1, Double_t k0, Double_t k1, Double_t s0, Double_t s1, Double_t tau);
  static Double_t QmaxCorrection(Int_t sector, Int_t row, Float_t cpad, Float_t ctime, Float_t ky, Float_t kz, Float_t rmsy0, Float_t rmsz0,  Float_t effLength=0, Float_t effDiff=1);
  static Double_t QtotCorrection(Int_t sector, Int_t row, Float_t cpad, Float_t ctime, Float_t ky, Float_t kz, Float_t rmsy0, Float_t rmsz0, Float_t qtot, Float_t thr,  Float_t effLength=0, Float_t effDiff=1);

  //
  //
  //
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
  TVectorD*& PosYcor(Int_t ind) {return fPosYcor[ind];}
  TVectorD*& PosZcor(Int_t ind) {return fPosZcor[ind];}
  Float_t ParamS0Par(Int_t i, Int_t j, Int_t k) const {return fParamS0Par[i][j][k];}
  TVectorD* QpadTnorm() const {return fQpadTnorm;}
  TVectorD* QpadMnorm() const {return fQpadMnorm;}

protected:
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
  // charge normalization parametrization
  //
  TObjArray *fQNorm;              // q norm paramters
  TMatrixD  *fQNormCorr;          // q norm correction for analytica  correction
  TObjArray *fQNormHis;           // q norm correction for analytical correction 
  //
  TVectorD  *fPosQTnorm[3];       // q position normalization
  TVectorD  *fPosQMnorm[3];       // q position normalization
  TVectorD  *fQpadTnorm;          // q pad normalization - Total charge
  TVectorD  *fQpadMnorm;          // q pad normalization - Max charge
  //
  // Position corrections
  // 
  TVectorD  *fPosYcor[3];       //  position correction parameterization 
  TVectorD  *fPosZcor[3];       //  position correction parameterization
  //
  // Wave Correction Map 
  //
  THnBase* fWaveCorrectionMap; //dY with respect to the distance to the center of the pad
  Bool_t   fWaveCorrectionMirroredPad;   // flag is the cog axis mirrored at 0.5
  Bool_t   fWaveCorrectionMirroredZ;     // flag is the Z axis mirrored at 0
  Bool_t   fWaveCorrectionMirroredAngle; // flag is the Angle axis mirrored at 0
  //
  // Resolution Map
  //
  THnBase* fResolutionYMap; // Map of resolution in Y
  //
  static AliTPCClusterParam*   fgInstance; //! Instance of this class (singleton implementation)
  ClassDef(AliTPCClusterParam,7)    //  TPC Cluster parameter class
};

#endif
