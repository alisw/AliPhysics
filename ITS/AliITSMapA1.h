#ifndef ALIITSMAPA1_H
#define ALIITSMAPA1_H


#include "AliITSMap.h"
class AliITSsegmentation;
class TObjArray;


class AliITSMapA1 :
  public AliITSMap 
{
  
public:
  AliITSMapA1() {
    // constructor
  }
  AliITSMapA1(AliITSsegmentation *seg);
  AliITSMapA1(AliITSsegmentation *seg, TObjArray *dig);
  AliITSMapA1(const AliITSMapA1 &source);
    // Assignment operator
  AliITSMapA1& operator=(const AliITSMapA1 &source);
  
  virtual ~AliITSMapA1();
    // Fill hits from list of digits into hit map
  virtual  void  FillMap();
    // Clear the hit map
  virtual  void  ClearMap();    
    // Set a single hit
  virtual  void  SetHit(Int_t iz, Int_t ix, Int_t idigit);
    // Delete a single hit
  virtual  void  DeleteHit(Int_t iz, Int_t ix);
    // Get index of hit in the list of digits
  virtual Int_t  GetHitIndex(Int_t iz, Int_t ix);
    // Get pointer to digit
  virtual TObject* GetHit(Int_t iz, Int_t ix);
    // Flag a hit as used
  virtual  void  FlagHit(Int_t iz, Int_t ix);    
    // Test hit status
  virtual FlagType TestHit(Int_t iz, Int_t ix);
    // Get signal from map
  virtual Double_t  GetSignal(Int_t iz, Int_t ix); 
    // Get max index inmap
  Int_t   MaxIndex() {return fMaxIndex;}
   // Set the array of objects
  void SetArray(TObjArray *obj);
  
protected:
    // Check index
  Int_t   CheckedIndex(Int_t iz, Int_t ix);

  AliITSsegmentation *fSegmentation;   // segmentation class
  Int_t fNpx;                          // fNpx
  Int_t fNpz;                          // fNpz
  TObjArray  *fObjects;                // object
  Int_t fNobjects;                     // number of objects
  Int_t fMaxIndex;                     // max index in map
  
private:
  Int_t *fHitMap;                      //! [fMaxIndex]

  ClassDef(AliITSMapA1,1)              // Implements Hit/Digit Map 
};

#endif	

