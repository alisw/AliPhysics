#ifndef AliD0_Trigger
#define AliD0_Trigger

#include "AliL3RootTypes.h"

class AliITStrackV2;

class AliD0Trigger {
 
 private:
  AliITStrackV2 * posTrack; //!
  AliITStrackV2 * negTrack; //!
  
  double momenta[6];
  double bestV0[3],primaryVertex[3];
  double cutV0low, cutV0high, cutInvMass, cutPointAngle;
  double Bfield;

 public:
  AliD0Trigger();
  AliD0Trigger(double c[4],double Bfield,double pv[3]);
  AliD0Trigger(AliITStrackV2 * posT, AliITStrackV2 * negT);
  virtual ~AliD0Trigger();

  void SetTracks(AliITStrackV2 * posT, AliITStrackV2 * negT);
  bool FindInvMass();
  bool FindV0();
  void FindMomentaAtVertex();
  bool PointingAngle();
  
  ClassDef(AliD0Trigger,1) 

};

#endif
