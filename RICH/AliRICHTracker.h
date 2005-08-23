#ifndef AliRICHTracker_h
#define AliRICHTracker_h

#include <AliTracker.h>
#include <AliLog.h>
#include "TNtupleD.h"

class AliCluster;
class AliESD;
class AliRICH;
class TTree;

class AliRICHTracker : public AliTracker
{
public:
           AliRICHTracker() :AliTracker()      {AliDebug(1,"Start.");}
  virtual ~AliRICHTracker()                    {AliDebug(1,"Stop.");}
  Int_t Clusters2Tracks(AliESD *)              {AliDebug(1,"Start.");return 0;} //pure virtual from AliTracker 
  Int_t RefitInward(AliESD *)                  {AliDebug(1,"Start.");return 0;} //pure virtual from AliTracker 
  void UnloadClusters()                        {AliDebug(1,"Start.");}          //pure virtual from AliTracker 
  AliCluster *GetCluster(Int_t )const          {AliDebug(1,"Start.");return 0;} //pure virtual from AliTracker 
  Int_t PropagateBack(AliESD *);                                                //pure virtual from AliTracker 
  void RecWithESD(AliESD *,AliRICH *,Int_t iTrackN);                            //recon with ESD
  void RecWithStack(TNtupleD *hn);                                              //recon from Stack in case ESD empty
  void CalcProb(Double_t thetaCer,Double_t pmod,Double_t *sigmaPID, Double_t *richPID);             // calculate pid for RICH
  Int_t LoadClusters(TTree *);                                                  //pure virtual from AliTracker 

protected:

  Double_t fField; // magnetic field stored

  ClassDef(AliRICHTracker,0)
};//class AliRICHTracker

#endif//AliRICHTracker_h
