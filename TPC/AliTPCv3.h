#ifndef TPCv3_H
#define TPCv3_H
////////////////////////////////////////////////
//  Version 3 for TPC                         //
////////////////////////////////////////////////
 
#include "AliTPC.h"

class AliTPCv3 : public AliTPC {

public:
  AliTPCv3() {}
  AliTPCv3(const char *name, const char *title);
  virtual      ~AliTPCv3() {}
  virtual void  CreateGeometry();
  virtual void  CreateMaterials();
  virtual void  Init();
  virtual Int_t IsVersion() const {return 3;}
  virtual void  StepManager();
  virtual void  DrawDetector();

protected:

  Int_t fIdSens1; // sensitive volume (entire drift gas)   

private:

  Float_t BetheBloch(Float_t bg);
  
  ClassDef(AliTPCv3,1)  // Time Projection Chamber version 3
};

#endif
