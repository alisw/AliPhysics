#ifndef AliHMPIDReconstructor_h
#define AliHMPIDReconstructor_h
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
//.
// HMPID base class to reconstruct an event
//.
#include <AliReconstructor.h>        //base class
#include "AliHMPIDTracker.h"         //CreateTracker()
#include "AliHMPIDDigit.h"           //Dig2Clu(), UseDig()
#include "AliHMPIDRecoParamV1.h"
#include <TMatrixF.h>                //UseDig()
#include <TClonesArray.h>            //UseDig()
#include <TObjArray.h>               //SigConv()

class AliRawReader;                  //Reconstruct() with raw data   
class AliHMPIDCluster;               //Dig2Clu()

class AliHMPIDReconstructor: public AliReconstructor 
{
public:
           AliHMPIDReconstructor();              
  virtual ~AliHMPIDReconstructor()                                  {delete fDig;delete fClu;delete [] fUserCut;}//dtor  
//framework part  
  AliTracker*  CreateTracker         () const {return new AliHMPIDTracker;}            //from AliReconstructor for clusters->PID
  void         ConvertDigits         (AliRawReader *pRR, TTree *pDigTree) const;                                        //from AliReconstruction for raw->digit
  Bool_t       HasDigitConversion()   const {return kTRUE;}                                                             //HMPID digits converted with ConvertDigits 
  void         Reconstruct           (TTree* digitsTree, TTree* clustersTree) const;                                    //from AliReconstruction for digit->cluster
  void         FillESD               (TTree* /*digitsTree*/, TTree* /*clustersTree*/, AliESDEvent *pESD)const;                                        //calculate pid for HMPID
  static Int_t StreamLevel()               { return fgStreamLevel;}
  static void  SetStreamLevel(Int_t level) { fgStreamLevel = level;}

    
  using AliReconstructor::FillESD;                                                                                      //
  using AliReconstructor::Reconstruct;                                                                                  // 

  //private part  
  static        void           Dig2Clu (TObjArray *pDigLst,TObjArray *pCluLst,Int_t *pUserCut,Bool_t isUnfold=kTRUE     );//digits->clusters
  static        void           FormClu (AliHMPIDCluster *pClu,AliHMPIDDigit *pDig,TClonesArray *pDigLst,TMatrixF *pPadMap);//cluster formation recursive algorithm
  static inline AliHMPIDDigit* UseDig  (Int_t padX,Int_t padY,                    TClonesArray *pDigLst,TMatrixF *pDigMap);//use this pad's digit to form a cluster
  inline Bool_t                IsDigSurvive(AliHMPIDDigit *pDig                                                     )const;//check for sigma cut
  static const AliHMPIDRecoParamV1* GetRecoParam() { return dynamic_cast<const AliHMPIDRecoParamV1*>(AliReconstructor::GetRecoParam(5)); }  //5 is the HMPID detector code
 
  protected:
  Int_t     *fUserCut;                 // n sigmas for pedestals decided by the User for each chamber(if in OCDB)
  TObjArray *fDaqSig;                  // container for the pad pedestal sigmas
  TObjArray *fDig;                     // tmp list of digits
  TObjArray *fClu;                     // tmp list of clusters
//
  private:
  AliHMPIDReconstructor(const AliHMPIDReconstructor&);              //Not implemented
  AliHMPIDReconstructor &operator=(const AliHMPIDReconstructor&);   //Not implemented
  static Int_t               fgStreamLevel; // flag for streaming   - for HMPID reconstruction  
//  
  ClassDef(AliHMPIDReconstructor, 3)        // class for the HMPID reconstruction
};

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
AliHMPIDDigit* AliHMPIDReconstructor::UseDig(Int_t padX,Int_t padY,TClonesArray *pDigLst,TMatrixF *pPadMap)
{
//Digit map contains a matrix if digit numbers.
//Main operation in forming initial cluster is done here. Requested digit pointer is returned and this digit marked as taken.
//Arguments: padX,padY - pad number
//           pDigLst   - list of digits for one sector
//           pDigMap   - map of those digits
//  Returns: pointer to digit if not yet used or 0 if used
  Int_t iDig=(Int_t)(*pPadMap)(padX,padY);(*pPadMap)(padX,padY)=-1;//take digit number from the map and reset this map cell to -1
  if(iDig!=-1)    return (AliHMPIDDigit*)pDigLst->At(iDig);        //digit pointer
  else            return 0;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Bool_t AliHMPIDReconstructor::IsDigSurvive(AliHMPIDDigit *pDig)const
{
//Check if the current digit survive to a riapllied sigma cut
//Arguments: pDig pointer to the current digit
//  Returns: kTRUE if charge > mean+n*sigma
  Int_t iCh = pDig->Ch();
  Int_t iDaqSigCut =(Int_t)fDaqSig->At(iCh)->GetUniqueID(); 
  if(fUserCut[iCh]<=iDaqSigCut) return kTRUE;
  TMatrixF *pM = (TMatrixF*)fDaqSig->At(pDig->Ch());
  Float_t sig = (*pM)(pDig->PadChX(),pDig->PadChY());
  if(pDig->Q()>fUserCut[iCh]*sig) return kTRUE;
  else return kFALSE;
}

#endif
