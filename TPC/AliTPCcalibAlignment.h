#ifndef ALI_TPC_CALIB_ALIGNMENT_H
#define ALI_TPC_CALIB_ALIGNMENT_H

#include "TNamed.h"
// ugly, but...
class AliTPCseed;
class AliExternalTrackParam;
class AliESDtrack;
class AliRieman;
class AliTPCtracker;
class TTreeSRedirector;

class AliTPCcalibAlignment:public TNamed {
public:
  AliTPCcalibAlignment();
  virtual ~AliTPCcalibAlignment();

  virtual void Process(AliTPCseed *track);

private:
  AliTPCcalibAlignment(const AliTPCcalibAlignment&);
  AliTPCcalibAlignment& operator=(const AliTPCcalibAlignment&);

  TTreeSRedirector *fDebugStream;

  ClassDef(AliTPCcalibAlignment,1);
};

#endif
