#ifndef PIPE_H
#define PIPE_H
////////////////////////////////////////////////
//  Manager class for detector: PIPE          //
////////////////////////////////////////////////
 
#include "AliModule.h"
 
 
class AliPIPE : public AliModule {
 
public:
  AliPIPE();
  AliPIPE(const char *name, const char *title);
  virtual      ~AliPIPE() {}
  
  ClassDef(AliPIPE,1)  //Beam Pipe base Class
};

#endif
