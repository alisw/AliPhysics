#ifndef ALIL3_ConfMapTrack
#define ALIL3_ConfMapTrack

#include <string.h>

#include "AliL3Track.h"

#include "AliL3RootTypes.h"
#include "AliL3ConfMapPoint.h"

class AliL3Vertex;

class AliL3ConfMapTrack :public AliL3Track {
  
 private:

  
 public:
  
  AliL3ConfMapTrack();
  virtual ~AliL3ConfMapTrack();
  void Fill(AliL3Vertex *vertex,Double_t max_Dca);
  void UpdateToFirstPoint();
  void Reset();
  void UpdateParam(AliL3ConfMapPoint *hit);
  void DeleteCandidate();

  void StartLoop() {currentHit = firstHit;}  //!
  void GetNextHit() {currentHit = ((AliL3ConfMapPoint*)currentHit)->nextTrackHit;} //!
  Int_t LoopDone() {return currentHit != 0;}  //!
  AliL3ConfMapPoint *currentHit;  //!
  AliL3ConfMapPoint *lastHit;  //!
  AliL3ConfMapPoint *firstHit;  //!


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
  
  ClassDef(AliL3ConfMapTrack,1) //Conformal mapping track class
};
    
#endif
    
