#ifndef SHIL_H
#define SHIL_H
////////////////////////////////////////////////
//  Manager class for detector: SHIL          //
////////////////////////////////////////////////
 
#include "AliDetector.h"
 
 
class AliSHIL : public AliDetector {
  
public:
  AliSHIL();
  AliSHIL(const char *name, const char *title);
  virtual      ~AliSHIL() {}
  virtual void  BuildGeometry();
  virtual void  CreateGeometry();
  virtual void  CreateMaterials();
  virtual void  Init();
  virtual Int_t IsVersion() const {return 0;}
  virtual void  StepManager();
  virtual void  DrawDetector();
 
  ClassDef(AliSHIL,1)  // Muon Shield Class
};

#endif
