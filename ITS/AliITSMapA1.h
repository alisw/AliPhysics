#ifndef ALIITSMAPA1_H
#define ALIITSMAPA1_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */
/* $Id$ */
////////////////////////////////////////////////////////////////////////
//  Map Class for ITS. Implementation A1. In this implementation, the //
// 2 dimensional (iz,ix) map is filled with integers values. For each //
// cell a corresponding TObject, a hit, can also be stored.           //
////////////////////////////////////////////////////////////////////////
#include "AliITSMap.h"
#include "TArrayI.h"

class AliITSsegmentation;
class TObjArray;

class AliITSMapA1 : public AliITSMap{

 public:
    AliITSMapA1(); // default creator
    // Standard reator using only a segmentation class
    AliITSMapA1(AliITSsegmentation *seg);
    // Standard reator using only a segmentation class and pointer to digits
    AliITSMapA1(AliITSsegmentation *seg, TObjArray *dig);
    // Standard reator using only a segmentation class and pointer to digits
    // and a threshold value
    AliITSMapA1(AliITSsegmentation *seg, TObjArray *dig, Int_t thr);
    AliITSMapA1(AliITSsegmentation *seg, TObjArray *dig, TArrayI thr);
    // Distructor
    virtual ~AliITSMapA1();
    // Fill hits from list of digits into hit map
    virtual  void  FillMap();
    virtual  void  FillMap2();
    // Clear the hit map
    virtual  void  ClearMap();    
    // Set a single hit
    virtual  void  SetHit(Int_t iz, Int_t ix, Int_t idigit);
    // Set threshold for the signal
    virtual  void  SetThreshold(Int_t thresh) {fMapThreshold=thresh;}
    virtual  void  SetThresholdArr(TArrayI th) {fMapThresholdArr=th;} 
    // Delete a single hit
    virtual  void  DeleteHit(Int_t iz, Int_t ix);
    // Get index of hit in the list of digits
    virtual Int_t  GetHitIndex(Int_t iz, Int_t ix) const ;
    // Get pointer to digit
    virtual TObject* GetHit(Int_t iz, Int_t ix) const;
    // Flag a hit as used
    virtual  void  FlagHit(Int_t iz, Int_t ix);
    // Test hit status
    virtual FlagType TestHit(Int_t iz, Int_t ix);
    // Get signal from map
    virtual Double_t  GetSignal(Int_t iz, Int_t ix) const;
    // Get max index inmap
    Int_t   MaxIndex() const  {return fMaxIndex;}
    // Set the array of objects
    void SetArray(TObjArray *obj);

 protected:
    // Copy Constructor
    AliITSMapA1(const AliITSMapA1 &source);
    // Assignment operator
    AliITSMapA1& operator=(const AliITSMapA1 &source);
    // Check index
    Int_t   CheckedIndex(Int_t iz, Int_t ix) const;

    // Data members
    AliITSsegmentation *fSegmentation;   // segmentation class
    Int_t fNpx;                          // fNpx
    Int_t fNpz;                          // fNpz
    TObjArray  *fObjects;                // object
    Int_t fNobjects;                     // number of objects
    Int_t fMaxIndex;                     // max index in map
    TArrayI fMapThresholdArr;            // array with thresholds
 private:
    Int_t *fHitMap;                      //! [fMaxIndex]
    Int_t fMapThreshold;                 // signal threshold (ADC)

    ClassDef(AliITSMapA1,2)  // Implements Hit/Digit Map 
};

#endif	

