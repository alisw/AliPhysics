// $Header$

#ifndef AliEveCLASS_H
#define AliEveCLASS_H

#include <Reve/Reve.h>

#include <TObject.h>

namespace Alieve {

class CLASS
{
private:
  CLASS(const CLASS&);            // Not implemented
  CLASS& operator=(const CLASS&); // Not implemented

protected:

public:
  CLASS();
  virtual ~CLASS() {}

  ClassDef(CLASS, 1);
}; // endclass CLASS

}

#endif
