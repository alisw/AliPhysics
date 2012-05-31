#ifndef ALITPCTRANSFORMATION_H
#define ALITPCTRANSFORMATION_H

//-------------------------------------------------------
//                       TPC transformations
//   
//
//   Origin: marian.ivanov@cern.ch
//           Code is not used anymore for the TPC corrections
//           Obsolete - will be removed soon 
// 
//-------------------------------------------------------


#include "TNamed.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TObjArray.h"
class TTreeSRedirector;
class AliTrackPointArray;
class AliTrackPoint;
class TFormula;
class TBits;
class TString;

class AliTPCTransformation: public TNamed{
public:
  typedef Double_t (*GenFuncG)(const Double_t*,const Double_t*);
  AliTPCTransformation();
  AliTPCTransformation(const char *name,TBits *mask, const char *fx, const char *fy, const char  *fz, Int_t coord);
  AliTPCTransformation(const AliTPCTransformation&trafo);
  ~AliTPCTransformation();
  //
  virtual Double_t GetDeltaXYZ(Int_t coord, Int_t volID, Double_t param, Double_t x, Double_t y, Double_t z);
  void SetParams(Double_t param, Double_t sigma, Double_t sigma2Time, const TVectorD *const fixedParams);
  Bool_t Init();
  void   SetActive(Bool_t flag){ fIsActive = flag;}
  Bool_t IsActive() const {return fIsActive;}
  //
  Double_t GetParam() const {return fParam;}
  void SetParam(Double_t param) {fParam=param;}
  Double_t GetSigma() const {return fSigma;}
  Double_t GetSigmaMax() const {return fSigmaMax;}
  Double_t GetSigma2Time() const {return fSigma2Time;}
  //
  static TBits * BitsSide(Bool_t aside);
  static TBits * BitsAll();
  static void RegisterFormula(const char * name, GenFuncG formula);
  static AliTPCTransformation::GenFuncG FindFormula(const char * name);
  static Double_t Eval(const char * name, const Double_t*x,const Double_t*par);

 private:
  //
  TString  * fNameX;         // x formula
  TString  * fNameY;         // y formula
  TString  * fNameZ;         // z formula  
  //  
  TBits    * fBitMask;       // bitmaps - transformation only for specified volID
  Int_t      fCoordSystem;   // coord system of  output deltas 
  Double_t   fParam;         // free parameter of transformation
  Double_t   fSigma;         // error of the parameter
  Double_t   fSigmaMax;      // maximal sigma (Not allowed to increase in propagate time by bigger factor)
  Double_t   fSigma2Time;    // change of the error in time (per hour) - (For kalman filter) 
  TVectorD  *fFixedParam;    // fixed parameters of tranformation
  Bool_t     fIsActive;      // switch - is transformation active
  //
  // predefined formulas
  //
  static  Int_t          BuildBasicFormulas(); //build list of basic formulas
  static  Double_t       TPCscalingRPol(Double_t *xyz, const Double_t * const param);
  static  Double_t       TPCscalingZDrift(Double_t *xyz, const Double_t * const param);
  static  Double_t       TPCscalingZDriftGy(Double_t *xyz, const Double_t * const param);
  static  Double_t       TPCscalingZDriftT0(Double_t *xyz, const Double_t * const param);
  static  Double_t       TPCscalingPhiLocal(Double_t *xyz, const Double_t * const param);
  static  Double_t       TPClocalRPhiEdge(Double_t *xyz, const Double_t *const param);
  //
  // TPC Field cage + ROC misalingment induced distortion
  //
  static  Double_t       TPCscalingRIFC(Double_t *xyz, const Double_t * const param); // inner field cage r distorion
  static  Double_t       TPCscalingROFC(Double_t *xyz, const Double_t * const param); // outer field cage r distorion
  //
  // TPC field cage + ROC misalignemnt induced distortion
  //
  static  Double_t       TPCdeltaFCROC(Double_t *xyz, const Double_t *const param); 
  static  Double_t       TPCdeltaFCCE(Double_t *xyz, const Double_t *const param); 

  //
  // TPC local misalignment
  //
  static  Double_t       TPClocaldLxdGX(Double_t *xyz, const Double_t *const param);
  static  Double_t       TPClocaldLxdGY(Double_t *xyz, const Double_t *const param);
  static  Double_t       TPClocaldLydGX(Double_t *xyz, const Double_t *const param);
  static  Double_t       TPClocaldLydGY(Double_t *xyz, const Double_t *const param);
  static  Double_t       TPClocaldRzdGX(Double_t *xyz, const Double_t *const param);
  static  Double_t       TPClocaldRzdGY(Double_t *xyz, const Double_t *const param);

  //
  // TPC  quadrant misalignment
  //
  //  static  Double_t       TPCQuadrantDr(Double_t *xyz, Double_t * param){return 0;}
  //static  Double_t       TPCQuadrantDrphi(Double_t *xyz, Double_t * param){return 0;}
  //
  // Z shift -
  //
  static  Double_t       TPCDeltaZ(Double_t *xyz, const Double_t *const param);
  static  Double_t       TPCDeltaZMediumLong(Double_t *xyz, Double_t * param);
  static  Double_t       TPCTiltingZ(Double_t *xyz, const Double_t *const param);
  //
  Bool_t    fInit;          // initialization flag
  GenFuncG  fFormulaX;      //! x formula
  GenFuncG  fFormulaY;      //! y formula
  GenFuncG  fFormulaZ;      //! z formula
  static  GenFuncG    fgFormulas[10000];   //! array of pointers to formula
  static  TObjArray*  fgFormulasName;      //! array of formalas name

  AliTPCTransformation &operator=(const AliTPCTransformation&);    // not implemented

  ClassDef(AliTPCTransformation,2);
};

#endif

