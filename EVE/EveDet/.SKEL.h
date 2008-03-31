// $Header$

#ifndef AliEveCLASS_H
#define AliEveCLASS_H

#include <Reve/Reve.h>

#include <TObject.h>

namespace Alieve {

class CLASS
{
public:
  CLASS();
  virtual ~CLASS() {}

protected:

private:
  CLASS(const CLASS&);            // Not implemented
  CLASS& operator=(const CLASS&); // Not implemented

  ClassDef(CLASS, 0);  // Short desc/purpose of CLASS
};

}

#endif
