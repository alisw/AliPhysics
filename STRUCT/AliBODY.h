#ifndef BODY_H
#define BODY_H
////////////////////////////////////////////////
//  Manager class for detector: BODY          //
//   This is the envelop for Alice            //
////////////////////////////////////////////////
 
#include "AliModule.h"
 
 
class AliBODY : public AliModule {
 
public:
   AliBODY();
   AliBODY(const char *name, const char *title);
   virtual      ~AliBODY() {}
   virtual void  CreateGeometry();
   virtual void  CreateMaterials();
   virtual Int_t IsVersion() const {return 0;}
   virtual void  DrawModule();

   ClassDef(AliBODY,1)  //Class manager for the ALICE body
};

#endif
