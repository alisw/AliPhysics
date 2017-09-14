#ifndef _ALINANOAODCUSTOMSETTER_H_
#define _ALINANOAODCUSTOMSETTER_H_

// AliNanoAODCustomSetter

// Virtual class which implements the basic interface for setting
// custom variables in special tracks and headers

// Author: Michele Floris, michele.floris@cern.ch

#include "TNamed.h"

class AliAODEvent;
class AliAODTrack;
class AliNanoAODHeader;
class AliNanoAODTrack;


class AliNanoAODCustomSetter : public TNamed
{
public:
  AliNanoAODCustomSetter(const char * name = "AliNanoAODCustomSetter") : TNamed(name,name) {;}
  virtual ~AliNanoAODCustomSetter() {;}
  virtual void SetNanoAODHeader(const AliAODEvent * event   , AliNanoAODHeader * head , TString varListHeader  ) =0;
  virtual void SetNanoAODTrack (const AliAODTrack * aodTrack, AliNanoAODTrack * spTrack) =0;

  ClassDef(AliNanoAODCustomSetter, 1)
};



#endif /* _ALINANOAODCUSTOMSETTER_H_ */
