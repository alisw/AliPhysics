#ifndef ABSO_H
#define ABSO_H
////////////////////////////////////////////////
//  Manager class for detector: ABSO          //
////////////////////////////////////////////////
 
#include "AliDetector.h"
 
 
class AliABSO : public AliDetector {
 
public:
  AliABSO();
  AliABSO(const char *name, const char *title);
  virtual      ~AliABSO() {}
  virtual void  BuildGeometry();
  virtual void  CreateGeometry();
  virtual void  CreateMaterials();
  virtual void  Init();
  virtual Int_t IsVersion() const {return 0;}
  virtual void  StepManager();
  virtual void  DrawDetector();
  
  
  ClassDef(AliABSO,1)  // Muon Absorber Class
};

#endif
