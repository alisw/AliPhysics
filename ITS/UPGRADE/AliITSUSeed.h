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
  void            SetClusterID(UInt_t id)  {fClID = id;}
  void            SetParent(TObject* par)  {fParent = par;}
  //
  UInt_t          GetClusterID()     const {return fClID;}
  TObject*        GetParent()        const {return fParent;}
  //
 protected:
  UInt_t                fClID;              // packed cluster info (see AliITSUAux::PackCluster)
  TObject*              fParent;            // parent track (in higher tree hierarchy)
  
  ClassDef(AliITSUSeed,1)
};


#endif
