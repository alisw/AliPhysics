#ifndef AliRICHTracker_h
#define AliRICHTracker_h

#include <AliTracker.h>
#include <AliLog.h>

class AliCluster;
class AliESD;
class TTree;

class AliRICHTracker : public AliTracker
{
public:
           AliRICHTracker() :AliTracker()      {AliDebug(1,"Start.");}
  virtual ~AliRICHTracker()                    {AliDebug(1,"Start.");}
  Int_t Clusters2Tracks(AliESD *)              {AliDebug(1,"Start.");return 0;} //pure virtual from AliTracker 
  Int_t RefitInward(AliESD *)                  {AliDebug(1,"Start.");return 0;} //pure virtual from AliTracker 
  void UnloadClusters()                        {AliDebug(1,"Start.");}          //pure virtual from AliTracker 
  AliCluster *GetCluster(Int_t )const          {AliDebug(1,"Start.");return 0;} //pure virtual from AliTracker 
  Int_t PropagateBack(AliESD *);                                                //pure virtual from AliTracker 
  Int_t LoadClusters(TTree *);                                                  //pure virtual from AliTracker 

protected:

  ClassDef(AliRICHTracker,0)
};//class AliRICHTracker

#endif//AliRICHTracker_h
