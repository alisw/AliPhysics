#ifndef DIPO_H
#define DIPO_H
////////////////////////////////////////////////
//  Manager class for Module: DIPO          //
////////////////////////////////////////////////
 
#include "AliModule.h"
 
 
class AliDIPO : public AliModule {
 
public:
  AliDIPO();
  AliDIPO(const char *name, const char *title);
  virtual      ~AliDIPO() {}
  virtual void  Init();
  
  ClassDef(AliDIPO,1)  //Class for the dipole magnet
};

#endif
