#ifndef PIPE_H
#define PIPE_H
////////////////////////////////////////////////
//  Manager class for detector: PIPE          //
////////////////////////////////////////////////
 
#include "AliDetector.h"
 
 
class AliPIPE : public AliDetector {
 
public:
  AliPIPE();
  AliPIPE(const char *name, const char *title);
  virtual      ~AliPIPE() {}
  virtual void  BuildGeometry();
  virtual void  CreateGeometry(){}
  virtual void  CreateMaterials(){}
  virtual Int_t IsVersion() const =0;
  virtual void  DrawDetector(){}
  virtual void  StepManager();
  
  ClassDef(AliPIPE,1)  //Beam Pipe base Class
};

#endif
