#ifndef HALL_H
#define HALL_H
////////////////////////////////////////////////
//  Manager class for detector: HALL          //
////////////////////////////////////////////////
 
#include "AliModule.h"
 
 
class AliHALL : public AliModule {
 
public:
   AliHALL();
   AliHALL(const char *name, const char *title);
   virtual      ~AliHALL() {}
   virtual void  CreateGeometry();
   virtual void  CreateMaterials();
   virtual void  Init();
   virtual Int_t IsVersion() const {return 0;}
   virtual void  DrawModule();
 
   ClassDef(AliHALL,1)  //Class for ALICE experimental hall
};

#endif
