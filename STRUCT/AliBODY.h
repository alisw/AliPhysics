#ifndef BODY_H
#define BODY_H
////////////////////////////////////////////////
//  Manager class for detector: BODY          //
//   This is the envelop for Alice            //
////////////////////////////////////////////////
 
#include "AliDetector.h"
 
 
class AliBODY : public AliDetector {
 
public:
   AliBODY();
   AliBODY(const char *name, const char *title);
   virtual      ~AliBODY() {}
   virtual void  BuildGeometry();
   virtual void  CreateGeometry();
   virtual void  CreateMaterials();
   virtual Int_t IsVersion() const {return 0;}
   virtual void  DrawDetector();
   virtual void  StepManager();
 
   ClassDef(AliBODY,1)  //Class manager for the ALICE body
};

#endif
