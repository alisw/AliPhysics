#ifndef ALITPCEXB_H
#define ALITPCEXB_H

class AliMagF;
#include "TObject.h"
#include "TVectorDfwd.h"

class AliTPCExB:public TObject {
public:
  AliTPCExB();
  AliTPCExB& operator=(const AliTPCExB &exb);
  AliTPCExB(const AliTPCExB& exb);
  virtual ~AliTPCExB() {};
  virtual void Correct(const Double_t *position,Double_t *corrected)=0;
  virtual void CorrectInverse(const Double_t *position,Double_t *corrected) {
    Correct(position,corrected);
    for (Int_t i=0;i<3;++i)
      corrected[i]=position[i]-(corrected[i]-position[i]);
  }
  //
  // Test and visualization
  //
  void TestExB(const char* fileName);
  static Double_t GetDr(Double_t r, Double_t phi, Double_t z, Double_t bz=5);
  static Double_t GetDrphi(Double_t r, Double_t phi, Double_t z, Double_t bz=5);
  static Double_t GetDphi(Double_t r, Double_t phi, Double_t z, Double_t bz=5);
  static Double_t GetDz(Double_t r, Double_t phi, Double_t z, Double_t bz=5);
  static AliTPCExB*  Instance(){return fgInstance;}
  static void SetInstance(AliTPCExB *const param){fgInstance = param;}
  //
  // Mag field scans
  //
  static  void RegisterField(Int_t index, AliMagF * magf);
  static  Double_t GetBx(Double_t r, Double_t phi, Double_t z,Int_t index=0);
  static  Double_t GetBy(Double_t r, Double_t phi, Double_t z,Int_t index=0);
  static  Double_t GetBz(Double_t r, Double_t phi, Double_t z,Int_t index=0);
  static  Double_t GetBr(Double_t r, Double_t phi, Double_t z,Int_t index=0);
  static  Double_t GetBrfi(Double_t r, Double_t phi, Double_t z,Int_t index=0);
  //
  static  Double_t GetBxI(Double_t r, Double_t phi, Double_t z,Int_t index=0);
  static  Double_t GetByI(Double_t r, Double_t phi, Double_t z,Int_t index=0);
  static  Double_t GetBzI(Double_t r, Double_t phi, Double_t z,Int_t index=0);
  static  Double_t GetBrI(Double_t r, Double_t phi, Double_t z,Int_t index=0);
  static  Double_t GetBrfiI(Double_t r, Double_t phi, Double_t z,Int_t index=0);
  //
  //
  Double_t Eval(Int_t type, Double_t r, Double_t phi, Double_t z);
  Double_t SEval(Int_t type, Double_t r, Double_t phi, Double_t z){return Instance()->Eval(type,r,phi,z);}
  static Double_t EvalMat(const TVectorD &vec, Double_t r, Double_t phi, Double_t z);     // evalute parameterization
 
 private:
  TVectorD *          fMatBrBz;       //param matrix Br/Bz
  TVectorD *          fMatBrfiBz;     //param matrix Br/Bz
  TVectorD *          fMatBrBzI0;     //param matrix Br/Bz integral  z>0 
  TVectorD *          fMatBrBzI1;     //param matrix Br/Bz integral  z<0 
  TVectorD *          fMatBrfiBzI0;   //param matrix Br/Bz integral  z>0 
  TVectorD *          fMatBrfiBzI1;   //param matrix Br/Bz integral  z<0
  
  static AliTPCExB*   fgInstance;  //! Instance of this class (singleton implementation)
  static TObjArray    fgArray;     //! array of magnetic fields
  //
  ClassDef(AliTPCExB,2)
};

#endif
