#ifndef HALL_H
#define HALL_H
////////////////////////////////////////////////
//  Manager class for detector: HALL          //
////////////////////////////////////////////////
 
#include "AliDetector.h"
 
 
class AliHALL : public AliDetector {
 
public:
   AliHALL();
   AliHALL(const char *name, const char *title);
   virtual      ~AliHALL() {}
   virtual void  BuildGeometry();
   virtual void  CreateGeometry();
   virtual void  CreateMaterials();
   virtual void  Init();
   virtual Int_t IsVersion() const {return 0;}
   virtual void  DrawDetector();
   virtual void  StepManager();
 
   ClassDef(AliHALL,1)  //Class for ALICE experimental hall
};

#endif
