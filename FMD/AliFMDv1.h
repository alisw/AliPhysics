#ifndef FMDV1_H
#define FMDV1_H
////////////////////////////////////////////////
//  Manager and hits classes for set:FMD     //
////////////////////////////////////////////////
 
#include "AliFMD.h"
 
class AliFMDv1 : public AliFMD {
  
public:
  AliFMDv1();
  AliFMDv1(const char *name, const char *title);
  virtual       ~AliFMDv1() {}
  virtual void   CreateGeometry();
  virtual void   CreateMaterials();
  virtual void   DrawDetector();
  virtual Int_t  IsVersion() const {return 1;}
  
  ClassDef(AliFMDv1,1)  //Class for FMD version 1
};
 
#endif
