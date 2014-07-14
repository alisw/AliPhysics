#ifndef ALIVVFRIENDTRACK_H
#define ALIVVFRIENDTRACK_H

#include "Rtypes.h"
class AliVVTPCseed;
class AliVVtrackPointArray;
class AliVVtrack;

//_____________________________________________________________________________
class AliVVfriendTrack {
public:
  AliVVfriendTrack() {}
  virtual ~AliVVfriendTrack() {}

  //used in calibration
  virtual AliVVtrack* GetTPCseed() const {return NULL;}
  virtual const AliVVtrackPointArray *GetTrackPointArray() const {return NULL;}
  virtual const AliVVtrack * GetITSOut() const {return NULL;} 
  virtual const AliVVtrack * GetTPCOut() const {return  NULL;} 
  virtual const AliVVtrack * GetTRDIn()  const {return NULL;} 

private: 
  AliVVfriendTrack(const AliVVfriendTrack &);
  AliVVfriendTrack& operator=(const AliVVfriendTrack& esd);  

  ClassDef(AliVVfriendTrack,1);
};

#endif


