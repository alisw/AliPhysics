#ifndef FRAMEv1_H
#define FRAMEv1_H
////////////////////////////////////////////////
//  Manager class for detector: FRAME          //
////////////////////////////////////////////////
 
#include "AliFRAME.h"
 
 
class AliFRAMEv1 : public AliFRAME {
 
public:
  AliFRAMEv1();
  AliFRAMEv1(const char *name, const char *title);
  virtual      ~AliFRAMEv1() {}
  virtual void  CreateGeometry();
  virtual void  CreateMaterials();
  virtual Int_t IsVersion() const {return 1;}
  virtual void  DrawDetector();
  
  ClassDef(AliFRAMEv1,1)  //Class for FRAME version 1
};

#endif
