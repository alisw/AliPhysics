#ifndef STARTV0_H
#define STARTV0_H
////////////////////////////////////////////////
//  Manager and hits classes for set:START     //
////////////////////////////////////////////////
 
#include "AliSTART.h"
 
class AliSTARTv0 : public AliSTART {
  
public:
  AliSTARTv0() {};
  AliSTARTv0(const char *name, const char *title);
  virtual       ~AliSTARTv0() {}
  virtual void   CreateGeometry();
  virtual void   CreateMaterials();
  virtual void   DrawDetector();
  virtual void   Init();
  virtual Int_t  IsVersion() const {return 0;}
  virtual void   StepManager();
  
protected:
   Int_t fIdSens1; // Sensetive volume  in START
 
  ClassDef(AliSTARTv0,1)  //Class for START version 0
};

#endif


