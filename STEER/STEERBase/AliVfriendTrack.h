#ifndef ALIVFRIENDTRACK_H
#define ALIVFRIENDTRACK_H

//_____________________________________________________________________________
#include "AliVMisc.h"

#include "Rtypes.h"
class AliTPCseed;
class AliVVtrackPointArray;

//_____________________________________________________________________________
class AliVfriendTrack {
public:

  AliVfriendTrack(){}
  // constructor for reinitialisation of vtable
  AliVfriendTrack( AliVConstructorReinitialisationFlag ){}

  virtual ~AliVfriendTrack(){}

  //used in calibration
  
  virtual Int_t GetTPCseed( AliTPCseed &) const = 0;
  //virtual const AliVtrackPointArray *GetTrackPointArray() const {return NULL;}
  //virtual const AliVtrack * GetITSOut() const {return NULL;} 
  //virtual const AliVtrack * GetTPCOut() const {return  NULL;} 
  //virtual const AliVtrack * GetTRDIn()  const {return NULL;} 
  // virtual TObject* GetCalibObject(Int_t /*index*/) const = 0;

private: 
  AliVfriendTrack(const AliVfriendTrack &);
  AliVfriendTrack& operator=(const AliVfriendTrack& esd);  

  //ClassDef(AliVfriendTrack,0);
};

#endif


