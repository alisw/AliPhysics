#ifndef PIPEv3_H
#define PIPEv3_H
////////////////////////////////////////////////
//  Manager class for detector: PIPE          //
////////////////////////////////////////////////
 
#include "AliPIPE.h"
 
 
class AliPIPEv3 : public AliPIPE {
 
public:
  AliPIPEv3();
  AliPIPEv3(const char *name, const char *title);
  virtual      ~AliPIPEv3() {}
  virtual void  CreateGeometry();
  virtual void  CreateMaterials();
  virtual Int_t IsVersion() const {return 3;}
  virtual void  DrawModule();
  virtual void  Undulation(char *, Float_t, Float_t, Float_t, Float_t,
                           char (*)[5]);
  ClassDef(AliPIPEv3,1)  //Class for PIPE version 3
};

#endif
