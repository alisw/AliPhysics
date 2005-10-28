#ifndef AliRICHReconstructor_h
#define AliRICHReconstructor_h
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <AliReconstructor.h>       //base class
#include "AliRICHTracker.h"         //CreateTracker()
#include "AliRICHClusterFinder.h"   //Reconstruct()

class AliRawReader;
class TTree;

class AliRICHReconstructor: public AliReconstructor 
{
public:
           AliRICHReconstructor(): AliReconstructor()              {}//default ctor
  virtual ~AliRICHReconstructor()                                  {}//dtor  
//framework part  
  AliTracker*  CreateTracker         (AliRunLoader*                      )const{return new AliRICHTracker;}                 //interface from AliReconstructor
  void         Reconstruct           (AliRunLoader* pAL                  )const{AliRICHClusterFinder cf(pAL);    cf.Exec();}//from AliReconstruction for digits->clusters
  void         Reconstruct           (AliRunLoader* pAL,AliRawReader *pRR)const;                                            //from AliReconstruction for raw->clusters
  using AliReconstructor::Reconstruct;                                                                            //to get rid of virtual hidden warning 
//private part  
         void          Dig2Clu    (TClonesArray*pDigList,TClonesArray *pCluList                                    )const;//form clusters out of provided digits list
         void          FormCluster(AliRICHCluster *pClu,AliRICHDigit *pDig,TClonesArray *pDigList,TMatrixF *pDigMap)const;//form cluster recursive algorithm
  inline AliRICHDigit *UseDig     (Int_t padX,Int_t padY,TClonesArray *pDigList,TMatrixF *pDigMap                  )const;//use this pad's digit to form a cluster
  static void          CheckPR    (                                                                                );     //utility-> run staff for stack
  static void          RichAna    (Int_t iNevMax=99999,Bool_t isPatRec=kFALSE                                      );     //utility-> create ntuples for analysis
  
protected:
  ClassDef(AliRICHReconstructor, 0)   //class for the RICH reconstruction
};
//__________________________________________________________________________________________________
AliRICHDigit* AliRICHReconstructor::UseDig(Int_t padX,Int_t padY,TClonesArray *pDigList,TMatrixF *pDigMap)const
{
//Digit map contains a matrix if digit numbers.
//Main operation in forming initial cluster is done here. Requested digit pointer is returned and this digit marked as taken.
//Arguments: padX,padY - pad number
//           pDigList  - list of digits for one sector
//           pDigMap   - map of those digits
//  Returns: pointer to digit if not yet used or 0 if used
  Int_t iDig=(Int_t)(*pDigMap)(padX,padY);(*pDigMap)(padX,padY)=-1;//take digit number from the map and reset this map cell to -1
  if(iDig!=-1)
    return (AliRICHDigit*)pDigList->At(iDig);    //digit pointer
  else      
    return 0;
}

#endif
