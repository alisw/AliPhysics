#ifndef SHIL_H
#define SHIL_H
////////////////////////////////////////////////
//  Manager class for Module: SHIL          //
////////////////////////////////////////////////
 
#include "AliModule.h"
 
 
class AliSHIL : public AliModule {
  
public:
  AliSHIL();
  AliSHIL(const char *name, const char *title);
  virtual      ~AliSHIL() {}
  virtual void  CreateGeometry();
  virtual void  CreateMaterials();
  virtual void  Init();
  virtual Int_t IsVersion() const {return 0;}
  virtual void  DrawModule();
 
  ClassDef(AliSHIL,1)  // Muon Shield Class
};

#endif
