#ifndef PHOSv3_H
#define PHOSv3_H
////////////////////////////////////////////////////////
//  Manager and hits classes for set:PHOS version 1   //
////////////////////////////////////////////////////////

// --- galice header files ---
#include "AliPHOS.h"
 
class AliPHOSv3 : public AliPHOS {

 public:
                        AliPHOSv3();
                        AliPHOSv3(const char *name, const char *title);
  virtual              ~AliPHOSv3(){}
  virtual void          CreateGeometry();
  virtual Int_t         IsVersion() const {return 3;}
  virtual void          StepManager();

 ClassDef(AliPHOSv3,1)  //Hits manager for set:PHOS version 3
};
 
#endif

