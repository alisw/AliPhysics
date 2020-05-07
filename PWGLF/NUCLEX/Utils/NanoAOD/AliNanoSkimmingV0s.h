#ifndef ALINANOSKIMMINGV0S_H
#define ALINANOSKIMMINGV0S_H

#include <Rtypes.h>
#include <AliAnalysisCuts.h>


class TObject;
class TList;

class AliESDtrackCuts;

class AliNanoSkimmingV0s : public AliAnalysisCuts
{
public:
  AliNanoSkimmingV0s();

  virtual bool IsSelected(TObject *obj);
  virtual bool IsSelected(TList *);

private:
  ClassDef(AliNanoSkimmingV0s, 2)
};

#endif