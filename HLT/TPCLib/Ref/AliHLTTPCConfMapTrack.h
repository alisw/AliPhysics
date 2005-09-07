// @(#) $Id$

#ifndef ALIHLTTPC_ConfMapTrack
#define ALIHLTTPC_ConfMapTrack

#include <string.h>

#include "AliHLTTPCTrack.h"

#include "AliHLTTPCRootTypes.h"
#include "AliHLTTPCConfMapPoint.h"

class AliHLTTPCVertex;

class AliHLTTPCConfMapTrack :public AliHLTTPCTrack {
  
 private:

  
 public:
  
  AliHLTTPCConfMapTrack();
  virtual ~AliHLTTPCConfMapTrack();
  void Fill(AliHLTTPCVertex *vertex,Double_t max_Dca);
  //void UpdateToFirstPoint();
  void Reset();
  void UpdateParam(AliHLTTPCConfMapPoint *hit);
  void DeleteCandidate();

  void StartLoop() {currentHit = firstHit;}  //!
  void GetNextHit() {currentHit = ((AliHLTTPCConfMapPoint*)currentHit)->nextTrackHit;} //!
  Int_t LoopDone() {return currentHit != 0;}  //!
  AliHLTTPCConfMapPoint *currentHit;  //!
  AliHLTTPCConfMapPoint *lastHit;  //!
  AliHLTTPCConfMapPoint *firstHit;  //!


  // getter
  Double_t const  *GetChiSq()   const { return fChiSq;}
  Int_t GetMCLabel();

  // setter   
  void  SetChiSq1(Double_t f) {fChiSq[0]=f;} 
  void  SetChiSq2(Double_t f) {fChiSq[1]=f;}
  void   SetProperties(Bool_t fUsage);

  Double_t fChiSq[2];
  
  //fit parameters. Bad naming convention, i know...
  Double_t    s11Xy  ; 
  Double_t    s12Xy  ;
  Double_t    s22Xy  ;
  Double_t    g1Xy   ;
  Double_t    g2Xy   ;       
  Double_t    s11Sz  ;
  Double_t    s12Sz  ;
  Double_t    s22Sz  ;
  Double_t    g1Sz   ;
  Double_t    g2Sz   ; 
  
  Double_t    ddXy, a1Xy, a2Xy ;    /*fit par in xy */
  Double_t    ddSz, a1Sz, a2Sz ;    /*fit par in sz */
  
  ClassDef(AliHLTTPCConfMapTrack,1) //Conformal mapping track class
};
    
#endif
    
