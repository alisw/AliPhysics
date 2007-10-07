#ifndef AliHMPIDTracker_h
#define AliHMPIDTracker_h

#include <AliTracker.h> //base class
#include "AliHMPID.h"   //Recon()
#include <AliRun.h>     //Recon()
#include <TF1.h>        //field
//.
// HMPID base class fo tracking
//.

class AliESDEvent;      //Recon()     
class AliESDtrack; //IntTrkCha()
class AliHMPIDTracker : public AliTracker
{
public:
           AliHMPIDTracker();
  virtual ~AliHMPIDTracker()                                            {delete fClu;}
//framework part  
         AliCluster *GetCluster     (Int_t                      )const  {return 0;} //pure virtual from AliTracker 
         Bool_t      GetTrackPoint  (Int_t idx,AliTrackPoint &pt)const;             //             from AliTracker  
         Int_t       Clusters2Tracks(AliESDEvent *                   )       {return 0;} //pure virtual from AliTracker 
         Int_t       LoadClusters   (TTree *pCluTr              );                  //pure virtual from AliTracker   
         Int_t       PropagateBack  (AliESDEvent *pEsd               );                  //pure virtual from AliTracker   
         Int_t       RefitInward    (AliESDEvent *                   )       {return 0;} //pure virtual from AliTracker 
         void        UnloadClusters (                           )       {         } //pure virtual from AliTracker 
//private part  
  static Int_t       IntTrkCha     (AliESDtrack *pTrk,Float_t &xPc,Float_t &yPc        );              //find track-PC intersection, retuns chamber ID
  static Int_t       Recon         (AliESDEvent *pEsd,TObjArray *pCluAll,TObjArray *pNmean=0,TObjArray *pQthre=0);//do actual job, returns status code  
  static Int_t       ReconHiddenTrk(Int_t iCh,AliESDtrack *pTrk,TClonesArray *pClus,TObjArray *pNmean, TObjArray *pQthre);//do actual job with Hidden Track Algorithm    
protected:
  TObjArray            *fClu;                     //! each chamber holds it's one list of clusters 
ClassDef(AliHMPIDTracker,0)
};//class AliHMPIDTracker
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


#endif//AliHMPIDTracker_h
