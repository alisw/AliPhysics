#ifndef AliHMPIDTracker_h
#define AliHMPIDTracker_h

#include <AliTracker.h> //base class
#include "AliHMPID.h"   //Recon()
#include <AliRun.h>     //Recon()

class AliESD;           

class AliHMPIDTracker : public AliTracker
{
public:
           AliHMPIDTracker():AliTracker()                               {} 
  virtual ~AliHMPIDTracker()                                            {}
//framework part  
         AliCluster *GetCluster     (Int_t                      )const  {return 0;} //pure virtual from AliTracker 
         Bool_t      GetTrackPoint  (Int_t idx,AliTrackPoint &pt)const;             //             from AliTracker  
         Int_t       Clusters2Tracks(AliESD *                   )       {return 0;} //pure virtual from AliTracker 
         Int_t       LoadClusters   (TTree *pCluTr              );                  //pure virtual from AliTracker   
  inline Int_t       PropagateBack  (AliESD *                   );                  //pure virtual from AliTracker   
         Int_t       RefitInward    (AliESD *                   )       {return 0;} //pure virtual from AliTracker 
         void        UnloadClusters (                           )       {         } //pure virtual from AliTracker 
//private part  
  enum ETrackingFlags {kMipDistCut=-9,kMipQdcCut=-5};
  static Int_t       Recon(AliESD *pEsd,TObjArray *pCluAll);                        //do actual job 
protected:
  ClassDef(AliHMPIDTracker,0)
};//class AliHMPIDTracker
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Int_t AliHMPIDTracker::PropagateBack(AliESD *pEsd)
{
// This method defined as pure virtual in AliTracker. It is invoked from AliReconstruction::RunTracking() after invocation of AliTracker::LoadClusters()
// Agruments: pEsd - pointer to ESD
//   Returns: error code    
  AliHMPID *pHmpid=((AliHMPID*)gAlice->GetDetector("HMPID"));  
  return Recon(pEsd,pHmpid->CluLst());  
}



typedef AliHMPIDTracker AliRICHTracker; // for backward compatibility

#endif//AliHMPIDTracker_h
