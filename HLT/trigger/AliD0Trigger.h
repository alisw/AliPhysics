#ifndef AliD0_Trigger
#define AliD0_Trigger

#include "AliHLTRootTypes.h"
#include <math.h>
#include <TObject.h>
#include <TObjArray.h>

class AliITStrackV2;

class AliD0Trigger : public TObjArray{
 
 private:
  AliITStrackV2 * posTrack; //!
  AliITStrackV2 * negTrack; //!
  
  double momenta[6];
  double bestV0[3],primaryVertex[3];
  double cutV0low, cutV0high, cutInvMass, cutPointAngle, cutd0d0,cutCosThetaStar,cutpTchild;
  double Bfield;

 public:
  AliD0Trigger();
  AliD0Trigger(double c[7],double Bfield,double pv[3]);
  AliD0Trigger(AliITStrackV2 * posT, AliITStrackV2 * negT);
  AliD0Trigger(double c[7],double Bfield,double pv[3],AliITStrackV2 * posT, AliITStrackV2 * negT);
  virtual ~AliD0Trigger();

  void SetTracks(AliITStrackV2 * posT, AliITStrackV2 * negT);
  void SetV0(double v[3]);
  bool FindInvMass();
  bool FindV0();
  bool FindV0offline(double v[3]);
  void FindMomentaAtVertex();
  void FindMomentaOffline();
  bool PointingAngle();
  void SetMomenta(double m[6]);
  bool d0d0();
  bool CosThetaStar();
  double P(){return sqrt(Pt()*Pt()+Pz()*Pz());} 
  double Pt(){return sqrt(Px()*Px()+Py()*Py());}
  double Px(){return (momenta[0]+momenta[3]);}
  double Py(){return (momenta[1]+momenta[4]);}
  double Pz(){return (momenta[2]+momenta[5]);}
  double Energy();
  //double Eta(){return atanh(cos(atan(Pt()/Pz())));}
  double Eta(){return 0.5*(log((P()+Pz())/(P()-Pz())));}
  bool pTchild();

  ClassDef(AliD0Trigger,1) 

};

#endif
