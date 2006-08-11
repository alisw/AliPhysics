#ifndef ALITRDTRIGGERL1_H
#define ALITRDTRIGGERL1_H

///////////////////////////////////////////////////////
//                                                   //
//  TRD trigger for L1                               //
//                                                   //
///////////////////////////////////////////////////////

#include "AliTriggerDetector.h"

class AliTRDTriggerL1 : public AliTriggerDetector
{

 public:

  AliTRDTriggerL1();
  virtual        ~AliTRDTriggerL1() {} 
  virtual void    CreateInputs();
  virtual void    Trigger();

  ClassDef(AliTRDTriggerL1,1)  // TRD Trigger Detector class

};

#endif
