#ifndef AliHMPIDTracker_h
#define AliHMPIDTracker_h

#include <AliTracker.h> //base class
#include "AliHMPID.h"   //Recon()
#include <AliRun.h>     //Recon()
#include <TF1.h>        //field 
class AliESD;      //Recon()     
class AliESDtrack; //IntTrkCha()
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
         Int_t       PropagateBack  (AliESD *                   );                  //pure virtual from AliTracker   
         Int_t       RefitInward    (AliESD *                   )       {return 0;} //pure virtual from AliTracker 
         void        UnloadClusters (                           )       {         } //pure virtual from AliTracker 
//private part  
  static Int_t       IntTrkCha(AliESDtrack *pTrk,Float_t &x,Float_t &y);                    //find track-chamber intersection, retuns chamber ID
  static Int_t       Recon    (AliESD *pEsd,TObjArray *pCluAll,TObjArray *pNmean=0);        //do actual job, returns status code  
protected:
  ClassDef(AliHMPIDTracker,0)
};//class AliHMPIDTracker
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


typedef AliHMPIDTracker AliRICHTracker; // for backward compatibility

#endif//AliHMPIDTracker_h
