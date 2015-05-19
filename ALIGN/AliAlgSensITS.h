#ifndef ALIALGSENSITS_H
#define ALIALGSENSITS_H

#include "AliAlgSens.h"


class TObjArray;


/*--------------------------------------------------------
  ITS sensor
  -------------------------------------------------------*/

// Author: ruben.shahoyan@cern.ch


class AliAlgSensITS : public AliAlgSens
{
 public:
  AliAlgSensITS(const char* name=0, Int_t vid=0, Int_t iid=0);
  virtual ~AliAlgSensITS();
  //
  virtual void   SetTrackingFrame();
  //
 protected:
  //
  ClassDef(AliAlgSensITS,1)
};


#endif
