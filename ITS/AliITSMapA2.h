#ifndef ALIITSMAPA2_H
#define ALIITSMAPA2_H


#include "AliITSMapA1.h"

class AliITSMapA2 :
public AliITSMapA1 
{

public:
  AliITSMapA2(AliITSsegmentation *seg);
  AliITSMapA2(AliITSsegmentation *seg,Int_t scalesizeX,Int_t scalesizeZ);
  AliITSMapA2(AliITSsegmentation *seg, TObjArray *hist,Double_t thresh);
  virtual ~AliITSMapA2();
  AliITSMapA2(const AliITSMapA2 &source); // copy constructor
   // assignment operator
  AliITSMapA2& operator=(const AliITSMapA2 &source);
   // fill pad signals into map 
  virtual  void  FillMap();
   // clear map
  virtual  void  ClearMap();    
   // set hit
  virtual  void  SetHit(Int_t iz, Int_t ix, Int_t signal){}
    // Flag a hit as used
  virtual  void  FlagHit(Int_t iz, Int_t ix);    
  virtual  void  DeleteHit(Int_t iz, Int_t ix);
    // Get index in the map
  virtual Int_t  GetHitIndex(Int_t iz, Int_t ix);
    // Get object (1D histogram)
  virtual TObject *GetHit(Int_t iz, Int_t dummy);
    // Test hit status
  virtual FlagType TestHit(Int_t iz, Int_t ix);
    // Get signal
  virtual Double_t  GetSignal(Int_t iz, Int_t ix);
   // set hit
  void  SetHit(Int_t iz, Int_t ix, Double_t signal);
    // Get signal
  Double_t  GetSignal(Int_t index);

private:
  void  FillMapFromHist();
  void  FillHist();
  void  ResetHist();

  Double_t *fHitMap;         //! [fMaxIndex]
  Double_t fMapThreshold;    // threshold for signal
  Int_t    fScaleSizeX;      // scale factor on x
  Int_t    fScaleSizeZ;      // scale factor on z

  ClassDef(AliITSMapA2,1) // Implements Signal Map
};


#endif	
