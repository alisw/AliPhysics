#ifndef MAG_H
#define MAG_H
////////////////////////////////////////////////
//  Manager class for detector: MAG           //
////////////////////////////////////////////////
 
#include "AliDetector.h"
 
 
class AliMAG : public AliDetector {
 
public:
   AliMAG();
   AliMAG(const char *name, const char *title);
   virtual      ~AliMAG() {}
   virtual void  BuildGeometry();
   virtual void  CreateGeometry();
   virtual void  CreateMaterials();
   virtual void  Init();
   virtual Int_t IsVersion() const {return 0;}
   virtual void  DrawDetector();
   virtual void  StepManager();
 
   ClassDef(AliMAG,1)  //Class manager for detector:MAG
};

#endif
