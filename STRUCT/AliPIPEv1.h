#ifndef PIPEv1_H
#define PIPEv1_H
////////////////////////////////////////////////
//  Manager class for detector: PIPE          //
////////////////////////////////////////////////
 
#include "AliPIPE.h"
 
 
class AliPIPEv1 : public AliPIPE {
 
public:
  AliPIPEv1();
  AliPIPEv1(const char *name, const char *title);
  virtual      ~AliPIPEv1() {}
  virtual void  CreateGeometry();
  virtual void  CreateMaterials();
  virtual Int_t IsVersion() const {return 1;}
  virtual void  DrawDetector();
  
  ClassDef(AliPIPEv1,1)  //Class for PIPE version 1
};

#endif
