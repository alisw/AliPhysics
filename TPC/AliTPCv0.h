#ifndef TPCv0_H
#define TPCv0_H
////////////////////////////////////////////////
//  Version 0 for TPC                         //
////////////////////////////////////////////////
 
#include "AliTPC.h"

class AliTPCv0 : public AliTPC {

public:
  AliTPCv0() {}
  AliTPCv0(const char *name, const char *title);
  virtual      ~AliTPCv0() {}
  virtual void  CreateGeometry();
  virtual void  CreateMaterials();
  virtual void  Init();
  virtual Int_t IsVersion() const {return 0;}
  virtual void  StepManager();
  virtual void  DrawModule();

protected:
  
  ClassDef(AliTPCv0,1)  // Time Projection Chamber version 0 - coarse
};

#endif

