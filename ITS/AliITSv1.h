#ifndef ITSv1_H
#define ITSv1_H
/////////////////////////////////////////////////////////
//  Manager and hits classes for set: ITS version 1    //
/////////////////////////////////////////////////////////
 
#include "AliITS.h"
 
class AliITSv1 : public AliITS {
 
public:
  AliITSv1();
  AliITSv1(const char *name, const char *title);
  virtual       ~AliITSv1() {}
  virtual void   CreateGeometry();
  virtual void   CreateMaterials();
  virtual void   Init(); 
  virtual Int_t  IsVersion() const {return 1;}
  virtual void   DrawDetector();
  virtual void   StepManager();
  
   ClassDef(AliITSv1,1)  //Hits manager for set:ITS version 1
};
 
#endif
