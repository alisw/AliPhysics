#ifndef PHOSv0_H
#define PHOSv0_H
/////////////////////////////////////////////////////////
//  Manager and hits classes for set:PHOS version 0    //
/////////////////////////////////////////////////////////

// --- galice header files ---
#include "AliPHOS.h"
 
class AliPHOSv0 : public AliPHOS
{

 protected:

  Int_t fIdSens; //Sensitive volume for phos

 public:
                        AliPHOSv0();
                        AliPHOSv0(const char *name, const char *title);
  virtual              ~AliPHOSv0(){}
  virtual void          CreateGeometry();
  virtual void          CreateMaterials();
  virtual void          Init();
  virtual Int_t         IsVersion() const {return 0;}
  virtual void          StepManager();

 ClassDef(AliPHOSv0,1)  //Hits manager for set:PHOS version 0
};
 
#endif

