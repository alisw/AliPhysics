#ifndef FRAME_H
#define FRAME_H
////////////////////////////////////////////////
//  Manager class for detector: FRAME         //
////////////////////////////////////////////////
 
#include "AliModule.h"


class AliFRAME : public AliModule {
  
public:
  AliFRAME();
  AliFRAME(const char *name, const char *title);
  virtual      ~AliFRAME() {}
  virtual Int_t IsVersion() const =0;
 
   ClassDef(AliFRAME,1)  //Class for Space Frame
};

#endif



