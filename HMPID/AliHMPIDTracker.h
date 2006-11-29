#ifndef AliHMPIDTracker_h
#define AliHMPIDTracker_h

#include <AliTracker.h> //base class

class TNtupleD;         //RecWithStack()   
class AliESD;           //Clusters2Tracks(), RefitInward(), PropagateBack(), RecWithESD()

class AliHMPIDTracker : public AliTracker
{
public:
           AliHMPIDTracker(); 
  virtual ~AliHMPIDTracker()                    {}
//framework part  
  AliCluster *GetCluster     (Int_t                      )const  {return 0;} //pure virtual from AliTracker 
  Bool_t      GetTrackPoint  (Int_t idx,AliTrackPoint &pt)const;             //             from AliTracker  
  Int_t       Clusters2Tracks(AliESD *                   )       {return 0;} //pure virtual from AliTracker 
  Int_t       LoadClusters   (TTree *pCluTr              );                  //pure virtual from AliTracker   
  Int_t       PropagateBack  (AliESD *                   );                  //pure virtual from AliTracker invoked from AliReconstruction::RunTracking()
  Int_t       RefitInward    (AliESD *                   )       {return 0;} //pure virtual from AliTracker 
  void        UnloadClusters (                           )       {         } //pure virtual from AliTracker 
//private part  
  enum ETrackingFlags {kMipDistCut=-9,kMipQdcCut=-5};
protected:
  ClassDef(AliHMPIDTracker,0)
};//class AliHMPIDTracker

typedef AliHMPIDTracker AliRICHTracker; // for backward compatibility

#endif//AliHMPIDTracker_h
