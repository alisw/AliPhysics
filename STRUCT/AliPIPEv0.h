#ifndef PIPEv0_H
#define PIPEv0_H
/////////////////////////////////////////////////////////
//  Manager and class for detector: PIPE  version 0    //
/////////////////////////////////////////////////////////
 
#include "AliPIPE.h"

class AliPIPEv0 : public AliPIPE {
  
public:
  AliPIPEv0();
  AliPIPEv0(const char *name, const char *title);
  virtual       ~AliPIPEv0() {}
  virtual void   CreateGeometry();
  virtual void   CreateMaterials();
  virtual Int_t  IsVersion() const {return 0;}
  virtual void  DrawDetector();
   
   ClassDef(AliPIPEv0,1)  //Class for PIPE version 0
};
 
#endif
