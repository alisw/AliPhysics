// $Header$

#ifndef ALIEVE_AliEVEHOMERSource_H
#define ALIEVE_AliEVEHOMERSource_H

#include <TEveElement.h>

#include <TNamed.h>

class AliHLTHOMERSourceDesc;

class AliEVEHOMERSource : public TEveElement,
  public TNamed
{
private:
  AliEVEHOMERSource(const AliEVEHOMERSource&);            // Not implemented
  AliEVEHOMERSource& operator=(const AliEVEHOMERSource&); // Not implemented

protected:
  AliHLTHOMERSourceDesc *fSource;

public:
  AliEVEHOMERSource(const Text_t* n="HOMER Source", const Text_t* t="");
  AliEVEHOMERSource(AliHLTHOMERSourceDesc* src, const Text_t* n="HOMER Source", const Text_t* t="");
  virtual ~AliEVEHOMERSource() {}

  AliHLTHOMERSourceDesc* GetSource() const   { return fSource; }
  void SetSource(AliHLTHOMERSourceDesc* src) { fSource = src; }

  ClassDef(AliEVEHOMERSource, 1);
}; // endclass AliEVEHOMERSource

#endif
