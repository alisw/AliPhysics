#ifndef ALIITSUCLUSTERPIX_H
#define ALIITSUCLUSTERPIX_H

#include "AliITSMFTClusterPix.h"

class AliITSUClusterPix : public AliITSMFTClusterPix
{
public:
  AliITSUClusterPix() : AliITSMFTClusterPix() {}
  AliITSUClusterPix(const AliITSUClusterPix& c) : AliITSMFTClusterPix(c) {}
  virtual ~AliITSUClusterPix() {}

protected:
  AliITSUClusterPix &operator=(const AliITSUClusterPix& cluster);

  ClassDef(AliITSUClusterPix,4)
};

#endif
