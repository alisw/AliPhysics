#ifndef ITSv5_H
#define ITSv5_H
/////////////////////////////////////////////////////////
//  Manager and hits classes for set: ITS version 4    //
/////////////////////////////////////////////////////////
 
#include "AliITS.h"
#include "TString.h"
 
class AliITSv5 : public AliITS {

private:

public:
  AliITSv5();
  AliITSv5(const char *name, const char *title);
  virtual       ~AliITSv5() {}
  virtual void   CreateGeometry();
  virtual void   CreateMaterials();
  virtual void   Init();   
  virtual Int_t  IsVersion() const {return 5;}
  virtual void   StepManager();
  
  ClassDef(AliITSv5,1)  //Hits manager for set:ITS version 4
};
 
#endif
