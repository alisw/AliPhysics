#ifndef ALIITSMAPA2_H
#define ALIITSMAPA2_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */
/* $Id$ */

////////////////////////////////////////////////////////////////////////
//  Map Class for ITS. Implementation A2. In this implementation, the //
// 2 dimensional (iz,ix) map is filled with Double precision floating //
// point values. For each cell a corresponding TObject, a hit, can    //
// also be stored. This class is derived from AliITSMapA1 so that is  //
// retains all of the functionality of that map class as well.        //
////////////////////////////////////////////////////////////////////////

#include "AliITSMapA1.h"

class AliITSMapA2 : public AliITSMapA1{

 public:
    AliITSMapA2(); // default creator
    // Standard reator using only a segmentation class
    AliITSMapA2(AliITSsegmentation *seg);
    // Standard reator using only a segmentation class and X and Z scale sizes
    AliITSMapA2(AliITSsegmentation *seg,Int_t scalesizeX,Int_t scalesizeZ);
    // Standard reator using only a segmentation class pointer to hits, and
    // a threshold value
    AliITSMapA2(AliITSsegmentation *seg, TObjArray *hist,Double_t thresh);
    virtual ~AliITSMapA2(); // destructor
    // fill pad signals into map 
    virtual void FillMap();
    // clear map
    virtual void ClearMap();    
    // set hit. Over written with a null function. See Double version below.
    virtual void SetHit(Int_t,Int_t,Int_t){}
    // set signal at a certain position in array
    void  SetHit(Int_t iz, Int_t ix, Double_t signal){
	fHitMapD[CheckedIndex(iz, ix)]=signal;}
    // set signal at a certain position in array
    void  SetHit(Int_t index, Double_t signal){fHitMapD[index]=signal;}
    // Flag a hit as used
    // Set threshold for the signal
    virtual void SetThreshold(Int_t thresh) {fMapThresholdD=(Double_t)thresh;}
    // flags hit in map
    virtual  void  FlagHit(Int_t iz, Int_t ix);
    //set the entry value to zero
    virtual  void  DeleteHit(Int_t iz, Int_t ix){
	fHitMapD[CheckedIndex(iz, ix)]=0;}
    //return the index of an entry in array
    virtual Int_t  GetHitIndex(Int_t iz, Int_t ix) const {
	return CheckedIndex(iz, ix);};
    // Get object (1D histogram)
    virtual TObject *GetHit(Int_t iz, Int_t /* dummy */) const;
    // Test hit status
    virtual FlagType TestHit(Int_t iz, Int_t ix);
    // Get signal using two dim. index
    virtual Double_t GetSignal(Int_t iz, Int_t ix) const
	{return GetSignal(GetHitIndex(iz,ix));}
    // Get signal
    Double_t GetSignal(Int_t index) const ;
    // Add new value to Map at cell
    virtual void AddSignal(Int_t iz, Int_t ix, Double_t sig);

 private:
    AliITSMapA2(const AliITSMapA2 &source); // copy constructor
    // assignment operator
    AliITSMapA2& operator=(const AliITSMapA2 &source);
    void  FillMapFromHist(); // fills the map from a historgram
    void  FillHist(); // fills a histogram from the map
    void  ResetHist(); // resets the histogram

    Double_t *fHitMapD;         //! [fMaxIndex]
    Double_t fMapThresholdD;    // threshold for signal
    Int_t    fScaleSizeX;       // scale factor on x
    Int_t    fScaleSizeZ;       // scale factor on z

    ClassDef(AliITSMapA2,1) // Implements Signal Map
};

#endif	
