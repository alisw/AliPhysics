#ifndef ALIITSMAP_H
#define ALIITSMAP_H


#include "AliITS.h"

typedef enum {kEmpty, kUsed, kUnused} Flag_t;

//___________________________________________________________________________

class AliITSMap :
  public TObject {

public:
  virtual ~AliITSMap() {}
  virtual  void  FillMap()                                       =0;
  virtual  void  ClearMap()                                      =0;
  virtual  void  SetHit(Int_t iz, Int_t ix, Int_t idigit)        =0;
  virtual  void  DeleteHit(Int_t iz, Int_t ix)                   =0;
  virtual  void   FlagHit(Int_t iz, Int_t ix)                    =0;    
  virtual Int_t  GetHitIndex(Int_t iz, Int_t ix)                 =0;
  virtual TObject * GetHit(Int_t iz, Int_t ix)                   =0;
  virtual Flag_t TestHit(Int_t iz, Int_t ix)                     =0;
  virtual Double_t  GetSignal(Int_t iz, Int_t ix)                =0;
  
  ClassDef(AliITSMap,1) //virtual base class for ITS Hit/Digit Map

    };


class AliITSMapA1 :
  public AliITSMap 
{
  
public:
  AliITSMapA1() {
    // constructor
  }
  AliITSMapA1(AliITSsegmentation *seg);
  AliITSMapA1(AliITSsegmentation *seg, TObjArray *dig);
  AliITSMapA1(const AliITSMapA1 &source); // copy constructor
  AliITSMapA1& operator=(const AliITSMapA1 &source); // assignment operator
  
  virtual ~AliITSMapA1();
  virtual  void  FillMap();
  virtual  void  ClearMap();    
  virtual Double_t  GetSignal(Int_t iz, Int_t ix) {
    // get signal
    return 0.;
  }
  virtual  void  SetHit(Int_t iz, Int_t ix, Int_t idigit);
  virtual void  DeleteHit(Int_t iz, Int_t ix);
  virtual Int_t  GetHitIndex(Int_t iz, Int_t ix);
  virtual TObject*  GetHit(Int_t iz, Int_t ix);
  virtual  void  FlagHit(Int_t iz, Int_t ix);    
  virtual Flag_t TestHit(Int_t iz, Int_t ix);
  Int_t  CheckedIndex(Int_t iz, Int_t ix);
  Int_t   MaxIndex() {
    // max index
    return fMaxIndex;
  }
  void SetArray(TObjArray *obj);
  
protected:
  AliITSsegmentation *fSegmentation;   // segmentation class
  Int_t fNpx;                          // fNpx
  Int_t fNpz;                          // fNpz
  TObjArray  *fObjects;                // object
  Int_t fNobjects;                     // nu of object
  Int_t fMaxIndex;                     // max index
  
private:
  Int_t *fHitMap;                      // hit map

  ClassDef(AliITSMapA1,1) // Implements Hit/Digit Map for SDD - read tree
    };


class AliITSMapA2 :
public AliITSMapA1 
{

public:
  AliITSMapA2(AliITSsegmentation *seg);
  AliITSMapA2(AliITSsegmentation *seg, TObjArray *hist,Double_t thresh);
  virtual ~AliITSMapA2();
  AliITSMapA2(const AliITSMapA2 &source); // copy constructor
  AliITSMapA2& operator=(const AliITSMapA2 &source); // assignment operator
  virtual  void  FillMap();
  virtual  void  ClearMap();    
  virtual  void  SetHit(Int_t iz, Int_t ix, Int_t signal){
    // set hit
  }
  virtual  void  FlagHit(Int_t iz, Int_t ix);    
  virtual  void  DeleteHit(Int_t iz, Int_t ix);
  virtual Int_t  GetHitIndex(Int_t iz, Int_t ix);
  virtual TObject * GetHit(Int_t iz, Int_t dummy);
  virtual Flag_t TestHit(Int_t iz, Int_t ix);
  virtual Double_t  GetSignal(Int_t iz, Int_t ix);
  void  SetHit(Int_t iz, Int_t ix, Double_t signal);
  Double_t  GetSignal(Int_t index);

private:
  Double_t *fHitMap;        // fHitMap
  Double_t fMapThreshold;   // fMapThreshold

  void  FillMapFromHist();
  void  FillHist();
  void  ResetHist();
  
  ClassDef(AliITSMapA2,1) // Implements Signal Map for SDD -fill or read hist
    };


#endif	




