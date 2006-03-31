#ifndef AliRICHTracker_h
#define AliRICHTracker_h

#include <AliTracker.h> //base class

class TNtupleD;         //RecWithStack()   
class AliESD;           //Clusters2Tracks(), RefitInward(), PropagateBack(), RecWithESD()

class AliRICHTracker : public AliTracker
{
public:
           AliRICHTracker(); 
  virtual ~AliRICHTracker()                    {}
//framework part  
  AliCluster *GetCluster     (Int_t                      )const  {return 0;} //pure virtual from AliTracker 
  Bool_t      GetTrackPoint  (Int_t idx,AliTrackPoint &pt)const;             //             from AliTracker  
  Int_t       Clusters2Tracks(AliESD *                   )       {return 0;} //pure virtual from AliTracker 
  Int_t       LoadClusters   (TTree *pCluTr              );                  //pure virtual from AliTracker   
  Int_t       PropagateBack  (AliESD *                   );                  //pure virtual from AliTracker invoked from AliReconstruction::RunTracking()
  Int_t       RefitInward    (AliESD *                   )       {return 0;} //pure virtual from AliTracker 
  void        UnloadClusters (                           )       {         } //pure virtual from AliTracker 
//private part  
         void RecWithStack(TNtupleD *hn                                                     );   //recon from Stack in case ESD empty
  static void CalcProb    (Double_t thetaCer,Double_t pmod,Double_t *pidsigma, Double_t *pid);   //calculate pid for RICH
  static void EsdPrint    (                                                                 );   //print ESD status 
  static void MatrixPrint (Double_t probCut=0.7                                             );   //print prob matrix with cut on probability    
  Double_t fErrPar[5];                                                                       //Temporary stored for debug purpose
  enum ETrackingFlags {kMipDistCut=-990,kMipQdcCut=-999};
protected:
  ClassDef(AliRICHTracker,0)
};//class AliRICHTracker

#endif//AliRICHTracker_h
