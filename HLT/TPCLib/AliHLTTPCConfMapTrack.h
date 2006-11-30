// @(#) $Id$
// Original: AliHLTConfMapTrack.h,v 1.4 2004/06/16 11:25:06 loizides 

#ifndef ALIHLTTPC_ConfMapTrack
#define ALIHLTTPC_ConfMapTrack

#include <string.h>

#include "AliHLTTPCTrack.h"

#include "AliHLTTPCRootTypes.h"
#include "AliHLTTPCConfMapPoint.h"

class AliHLTTPCVertex;

class AliHLTTPCConfMapTrack :public AliHLTTPCTrack {

 public:
  
  AliHLTTPCConfMapTrack();
  virtual ~AliHLTTPCConfMapTrack();
  void Fill(AliHLTTPCVertex *vertex,Double_t max_Dca);
  void Reset();
  void UpdateParam(AliHLTTPCConfMapPoint *hit);
  void DeleteCandidate();

  void StartLoop() {fCurrentHit = fFirstHit;}
  void GetNextHit() {fCurrentHit = ((AliHLTTPCConfMapPoint*)fCurrentHit)->GetNextTrackHit();}
  Int_t LoopDone() const {return fCurrentHit != 0;}

  // setter   
  void  SetChiSq1(Double_t f) {fChiSq[0]=f;} 
  void  SetChiSq2(Double_t f) {fChiSq[1]=f;}
  void  SetProperties(Bool_t fUsage);

  // getter
  Double_t const  *GetChiSq()   const { return fChiSq;}
  Double_t GetChiSq1() const { return fChiSq[0]; }
  Double_t GetChiSq2() const { return fChiSq[1]; }

  /*
  Double_t GetS11Xy() const {return fs11Xy;}
  Double_t GetS12Xy() const {return fs12Xy;}
  Double_t GetS22Xy() const {return fs22Xy;}
  Double_t GetG1Xy()  const {return fg1Xy;}
  Double_t GetG2Xy()  const {return fg2Xy;}
  Double_t GetS11Sz() const {return fs11Sz;}
  Double_t GetS12Sz() const {return fs12Sz;}
  Double_t GetS22z()  const {return fs22Sz;}
  Double_t GetG1Sz()  const {return fg1Sz;}
  Double_t GetG2Sz()  const { return fg2Sz;}
  */

  Double_t GetDDXy() const {return fddXy;} 
  Double_t GetA1Xy() const {return fa1Xy;} 
  Double_t GetA2Xy() const {return fa2Xy;}   
  Double_t GetDDSz() const {return fddSz;} 
  Double_t GetA1Sz() const {return fa1Sz;} 
  Double_t GetA2Sz() const {return fa2Sz;}  

  AliHLTTPCConfMapPoint* GetFirstHit()   const {return fFirstHit;}
  AliHLTTPCConfMapPoint* GetLastHit()    const {return fLastHit;}
  AliHLTTPCConfMapPoint* GetCurrentHit() const {return fCurrentHit;}
  Int_t GetMCLabel();

 protected:

  AliHLTTPCConfMapPoint *fCurrentHit;  //!
  AliHLTTPCConfMapPoint *fLastHit;  //!
  AliHLTTPCConfMapPoint *fFirstHit;  //!


  Double_t fChiSq[2]; //chi squared
  
  //fit parameters. Bad naming convention, i know...
  Double_t    fs11Xy; //helper
  Double_t    fs12Xy; //helper
  Double_t    fs22Xy; //helper
  Double_t    fg1Xy;  //helper
  Double_t    fg2Xy;  //helper       
  Double_t    fs11Sz; //helper
  Double_t    fs12Sz; //helper
  Double_t    fs22Sz; //helper
  Double_t    fg1Sz;  //helper
  Double_t    fg2Sz;  //helper 
  
  Double_t    fddXy, fa1Xy, fa2Xy ;    /*fit par in xy */
  Double_t    fddSz, fa1Sz, fa2Sz ;    /*fit par in sz */
  
  ClassDef(AliHLTTPCConfMapTrack,1) //Conformal mapping track class
};
    
#endif
    
