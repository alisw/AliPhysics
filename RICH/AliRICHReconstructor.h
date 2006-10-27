#ifndef AliRICHReconstructor_h
#define AliRICHReconstructor_h
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <AliReconstructor.h>       //base class
#include "AliRICHTracker.h"         //CreateTracker()
#include <TMatrixF.h>               //UseDig()
#include <TClonesArray.h>           //UseDig()
class AliRawReader;                 //Reconstruct() with raw data   
class AliRICHDigit;                 //Dig2Clu(), UseDig()
class AliRICHCluster;               //Dig2Clu()

class AliRICHReconstructor: public AliReconstructor 
{
public:
           AliRICHReconstructor(): AliReconstructor()              {}//default ctor
  virtual ~AliRICHReconstructor()                                  {}//dtor  
//framework part  
  AliTracker*  CreateTracker         (AliRunLoader*                      )const{return new AliRICHTracker;}            //from AliReconstructor for clusters->PID
  void         Reconstruct           (AliRunLoader* pAL                  )const;                                       //from AliReconstruction for digits->clusters
  void         Reconstruct           (AliRunLoader* pAL,AliRawReader *pRR)const;                                       //from AliReconstruction for raws->clusters
  virtual void FillESD               (AliRunLoader* pAL,AliESD *pESD)const;                                    //calculate pid for RICH
  virtual void FillESD(AliRunLoader*, AliRawReader*, AliESD*) const { };
  virtual void FillESD(AliRawReader*, TTree*, AliESD*) const { };
  virtual void FillESD(TTree*, TTree*, AliESD*) const { };

  
   using AliReconstructor::Reconstruct;                                                                                 //to get rid of virtual hidden warning 

  //private part  
  static        void          Dig2Clu (TClonesArray*pDigLst,TClonesArray *pCluLst,Bool_t isTryUnfold=kTRUE            );//digits list -> clusters list
  static        void          CluQA   (AliRunLoader* pAL                                                              );//QA for clusters
  static        void          FormClu (AliRICHCluster *pClu,AliRICHDigit *pDig,TClonesArray *pDigLst,TMatrixF *pDigMap);//cluster formation recursive algorithm
  static inline AliRICHDigit* UseDig  (Int_t padX,Int_t padY,TClonesArray *pDigList,TMatrixF *pDigMap                 );//use this pad's digit to form a cluster

  protected:
  ClassDef(AliRICHReconstructor, 0)   //class for the RICH reconstruction
};
//__________________________________________________________________________________________________
AliRICHDigit* AliRICHReconstructor::UseDig(Int_t padX,Int_t padY,TClonesArray *pDigLst,TMatrixF *pDigMap)
{
//Digit map contains a matrix if digit numbers.
//Main operation in forming initial cluster is done here. Requested digit pointer is returned and this digit marked as taken.
//Arguments: padX,padY - pad number
//           pDigLst   - list of digits for one sector
//           pDigMap   - map of those digits
//  Returns: pointer to digit if not yet used or 0 if used
  Int_t iDig=(Int_t)(*pDigMap)(padX,padY);(*pDigMap)(padX,padY)=-1;//take digit number from the map and reset this map cell to -1
  if(iDig!=-1)    return (AliRICHDigit*)pDigLst->At(iDig);         //digit pointer
  else            return 0;
}

#endif
