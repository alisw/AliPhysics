#ifndef ALIVVFRIENDTRACK_H
#define ALIVVFRIENDTRACK_H

//_____________________________________________________________________________
#include "AliVVMisc.h"

class AliVVTPCseed;
class AliVVTRDseed;
class AliTrackPointArray;
class AliExternalTrackParam;

class AliVVfriendTrack {
public:
  AliVVfriendTrack(){}
  AliVVfriendTrack(const AliVVfriendTrack &){}
  AliVVfriendTrack& operator=(const AliVVfriendTrack& ){}
  // constructor for reinitialisation of vtable
  AliVVfriendTrack( AliVVConstructorReinitialisationFlag ){}

  virtual ~AliVVfriendTrack(){}

  //used in calibration
  virtual AliVVTPCseed* GetTPCseed() const {return NULL;}
  virtual AliVVTRDseed* GetTRDseed() const {return NULL;}
  virtual const AliTrackPointArray *GetTrackPointArray() const {return NULL;}
  virtual const AliExternalTrackParam * GetITSOut() const {return NULL;} 
  virtual const AliExternalTrackParam * GetTPCOut() const {return  NULL;} 
  virtual const AliExternalTrackParam * GetTRDIn()  const {return NULL;} 

private: 
};

#endif


