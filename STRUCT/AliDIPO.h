#ifndef DIPO_H
#define DIPO_H
////////////////////////////////////////////////
//  Manager class for detector: DIPO          //
////////////////////////////////////////////////
 
#include "AliDetector.h"
 
 
class AliDIPO : public AliDetector {
 
public:
  AliDIPO();
  AliDIPO(const char *name, const char *title);
  virtual      ~AliDIPO() {}
  virtual void  BuildGeometry();
  virtual void  CreateGeometry()=0;
  virtual void  CreateMaterials()=0;
  virtual void  Init();
  virtual Int_t IsVersion() const {return 0;}
  virtual void  DrawDetector()=0;
  virtual void  StepManager()=0;
  
  ClassDef(AliDIPO,1)  //Class for the dipole magnet
};

#endif
