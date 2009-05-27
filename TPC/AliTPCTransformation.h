#ifndef ALITPCTRANSFORMATION_H
#define ALITPCTRANSFORMATION_H

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
  AliTPCTransformation(const char *name,TBits *mask, const char *fx, const char *fy, const char  *fz, Int_t coord, Double_t param, Double_t sigma, TVectorD *fixedParams);
  Bool_t Init();
  //
  virtual Double_t GetDeltaXYZ(Int_t coord, Int_t volID, Double_t param, Double_t x, Double_t y, Double_t z);
  static TBits * BitsSide(Bool_t aside);
  static TBits * BitsAll();

  static void RegisterFormula(const char * name, GenFuncG formula);
  static AliTPCTransformation::GenFuncG FindFormula(const char * name);
  static Double_t Eval(const char * name, const Double_t*x,const Double_t*par);
public:
  //
  TString  * fNameX;         // x formula
  TString  * fNameY;         // y formula
  TString  * fNameZ;         // z formula  
  //  
  TBits    * fBitMask;       // bitmaps - transformation only for specified volID
  Int_t      fCoordSystem;   // coord system of  output deltas 
  Double_t   fParam;         // free parameter of transformation
  Double_t   fSigma;         // error of the parameter
  TVectorD  *fFixedParam;   // fixed parameters of tranformation
  
  //
  // predefined formulas
  //
  static  Int_t          BuildBasicFormulas(); //build list of basic formulas
  static  Double_t       TPCscalingRPol(Double_t *xyz, Double_t * param);
  static  Double_t       TPCscalingZDr(Double_t *xyz, Double_t * param);
  static  Double_t       TPCscalingPhiLocal(Double_t *xyz, Double_t * param);
  //
  Bool_t    fInit;          // initialization flag
  GenFuncG  fFormulaX;      //! x formula
  GenFuncG  fFormulaY;      //! y formula
  GenFuncG  fFormulaZ;      //! z formula
  static  GenFuncG    fgFormulas[10000];   //! array of pointers to formula
  static  TObjArray*  fgFormulasName;      //! array of formalas name
private:

  ClassDef(AliTPCTransformation,1);
};

#endif

