#ifndef FRAME_H
#define FRAME_H
////////////////////////////////////////////////
//  Manager class for detector: FRAME         //
////////////////////////////////////////////////
 
#include "AliDetector.h"


class AliFRAME : public AliDetector {
  
public:
  AliFRAME();
  AliFRAME(const char *name, const char *title);
  virtual      ~AliFRAME() {}
  virtual void  BuildGeometry();
  virtual void  CreateGeometry(){}
  virtual void  CreateMaterials(){}
  virtual Int_t IsVersion() const =0;
  virtual void  DrawDetector(){}
  virtual void  StepManager();
  virtual void  Init(){}
 
   ClassDef(AliFRAME,1)  //Class for Space Frame
};

#endif



