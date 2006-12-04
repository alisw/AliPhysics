// $Header$

#ifndef REVE_RMacro_H
#define REVE_RMacro_H

#include <Reve/Reve.h>

#include <TMacro.h>

namespace Reve {

class RMacro : public TMacro
{
protected:

public:
  RMacro();
  RMacro(const RMacro&);
  RMacro(const char* name);
  virtual ~RMacro() {}

  virtual void Exec(const char* params = "0");

  void ResetRoot();

  ClassDef(RMacro, 1);
}; // endclass RMacro

}

#endif
