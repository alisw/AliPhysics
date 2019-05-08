#ifndef ALINANOSKIMMINGPID_H
#define ALINANOSKIMMINGPID_H

#include <Rtypes.h>

#include <AliAnalysisCuts.h>
#include "AliNanoFilterPID.h"

class TObject;
class TList;

class AliESDtrackCuts;

class AliNanoSkimmingPID : public AliAnalysisCuts {
public:
  AliNanoSkimmingPID();

  virtual bool IsSelected(TObject *obj);
  virtual bool IsSelected(TList *);

  AliNanoFilterPID fTrackFilter;

private:
  ClassDef(AliNanoSkimmingPID, 1)
};

#endif