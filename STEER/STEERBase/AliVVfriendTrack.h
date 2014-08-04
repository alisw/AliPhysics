#ifndef ALIVVFRIENDTRACK_H
#define ALIVVFRIENDTRACK_H

//_____________________________________________________________________________
#include "AliVVMisc.h"

#include "Rtypes.h"
class AliTPCseed;
class AliVVtrackPointArray;

//_____________________________________________________________________________
class AliVVfriendTrack {
public:

  AliVVfriendTrack(){}
  // constructor for reinitialisation of vtable
  AliVVfriendTrack( AliVVConstructorReinitialisationFlag ){}

  virtual ~AliVVfriendTrack(){}

  //used in calibration
  
  virtual Int_t GetTPCseed( AliTPCseed &) const = 0;
  //virtual const AliVVtrackPointArray *GetTrackPointArray() const {return NULL;}
  //virtual const AliVVtrack * GetITSOut() const {return NULL;} 
  //virtual const AliVVtrack * GetTPCOut() const {return  NULL;} 
  //virtual const AliVVtrack * GetTRDIn()  const {return NULL;} 
  // virtual TObject* GetCalibObject(Int_t /*index*/) const = 0;

private: 
  AliVVfriendTrack(const AliVVfriendTrack &);
  AliVVfriendTrack& operator=(const AliVVfriendTrack& esd);  

  ClassDef(AliVVfriendTrack,1);
};

#endif


