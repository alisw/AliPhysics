#ifndef ALIITSUSEED_H
#define ALIITSUSEED_H

#include "AliExternalTrackParam.h"


class AliITSUSeed: public AliExternalTrackParam
{
 public:
  AliITSUSeed();
  AliITSUSeed(const AliITSUSeed& src);
  AliITSUSeed &operator=(const AliITSUSeed &src);
  virtual ~AliITSUSeed();
  //
 protected:
  Double_t         fMass;
  
  ClassDef(AliITSUSeed,1)
};


#endif
