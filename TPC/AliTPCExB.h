#ifndef ALITPC_EXB
#define ALITPC_EXB

#include "TObject.h"

class AliTPCExB:public TObject {
public:
  virtual ~AliTPCExB() {};
  virtual void Correct(const Double_t *position,Double_t *corrected)=0;
  virtual void CorrectInverse(const Double_t *position,Double_t *corrected) {
    Correct(position,corrected);
    for (Int_t i=0;i<3;++i)
      corrected[i]=position[i]-(corrected[i]-position[i]);
  }
  //
  // Test and visulaization
  //
  void TestExB(const char* fileName);
  static Double_t GetDr(Double_t r, Double_t phi, Double_t z);
  static Double_t GetDrphi(Double_t r, Double_t phi, Double_t z);
  static Double_t GetDphi(Double_t r, Double_t phi, Double_t z);
  static Double_t GetDz(Double_t r, Double_t phi, Double_t z);
  static AliTPCExB*  Instance(){return fgInstance;}
  static void SetInstance(AliTPCExB*param){fgInstance = param;}
 protected:
  static AliTPCExB*   fgInstance; //! Instance of this class (singleton implementation)
  ClassDef(AliTPCExB,0)
};

#endif
