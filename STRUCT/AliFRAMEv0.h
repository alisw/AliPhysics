#ifndef FRAMEv0_H
#define FRAMEv0_H
/////////////////////////////////////////////////////////
//  Manager and class for detector: FRAME  version 0    //
/////////////////////////////////////////////////////////
 
#include "AliFRAME.h"

class AliFRAMEv0 : public AliFRAME {
  
public:
  AliFRAMEv0();
  AliFRAMEv0(const char *name, const char *title);
  virtual       ~AliFRAMEv0() {}
  virtual void   CreateGeometry();
  virtual void   CreateMaterials();
  virtual Int_t  IsVersion() const {return 0;}
  virtual void  DrawDetector();
   
   ClassDef(AliFRAMEv0,1)  //Class for FRAME version 0
};
 
#endif
