#ifndef ALIITSPLIST_H
#define ALIITSPLIST_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */

#include "TClonesArray.h"
#include "AliITSMap.h"
#include "AliITSpListItem.h"

class AliITSpList: public AliITSMap {

 public:
    // Default Constructor
    AliITSpList();
    // Standard Constructor
    AliITSpList(Int_t imax,Int_t jmax);
    // Class destrutor
    virtual ~AliITSpList();
    // Copy constructor
    AliITSpList(const AliITSpList &source);
    // = Operator
    virtual AliITSpList& operator=(const AliITSpList &source);
    // Returns the max mape index values
    void GetMaxMapIndex(Int_t &ni,Int_t &nj) const {ni=fNi;nj=fNj;return;}
    // returns the max index value.
    Int_t GetMaxIndex() const {return fNi*fNj;}
    // returns the largest non-zero entry kept in the array fa.
    Int_t GetEntries() const {return fEntries;}
    // returns the max number of track/hit entries per cell.
    Int_t GetNEntries() const {return AliITSpListItem::GetMaxKept();}
    // for a given TClonesArray index it returns the corresponding map index
    void  GetMapIndex(Int_t index,Int_t &i,Int_t &j) const {
	i = index/fNj;j = index - fNj*i;
	if(i<0||i>=fNi || j<0||j>=fNj){i=-1;j=-1; return;}
    }
    // Returns the signal+noise for a give map coordinate
    Double_t GetSignal(Int_t index) {
	if(GetpListItem(index)==0) return 0.0;
	return GetpListItem(index)->GetSumSignal();
    }
    // Returns the signal+noise for a give map coordinate
    virtual Double_t GetSignal(Int_t i,Int_t j)  {
	if(GetpListItem(i,j)==0) return 0.0;
	return GetpListItem(i,j)->GetSumSignal();
    }
    // Returns the signal only for a give map coordinate
    Double_t GetSignalOnly(Int_t i,Int_t j) {
	if(GetpListItem(i,j)==0) return 0.0;
	return GetpListItem(i,j)->GetSignal();
    }
    // Returns the noise for a give map coordinate
    Double_t GetNoise(Int_t i,Int_t j) {
	if(GetpListItem(i,j)==0) return 0.0;
	return GetpListItem(i,j)->GetNoise();
    }
    // returns the track number which generated the signal at a given map
    // coordinate. If there is no signal or only noise, then -2 is returned.
    // k is the track rank number.
    Double_t GetTSignal(Int_t i,Int_t j,Int_t k) {
	if(GetpListItem(i,j)==0) return 0.0;
	return GetpListItem(i,j)->GetSignal(k);
    }
    // returns the track number which generated the signal at a given map
    // coordinate. If there is no signal or only noise, then -2 is returned.
    // k is the track rank number.
    Int_t GetTrack(Int_t i,Int_t j,Int_t k) {
	if(GetpListItem(i,j)==0) return -2;
	return GetpListItem(i,j)->GetTrack(k);
    }
    // returns the hit number which generated the signal at a given map
    // coordinate. If there is no signal or only noise, then -2 is returned.
    // k is the hit rank number.
    Int_t GetHit(Int_t i,Int_t j,Int_t k) {
	if(GetpListItem(i,j)==0) return -2;
	return GetpListItem(i,j)->GetHit(k);
    }
    // returns the number of Signal values
    Int_t GetNSignals(Int_t i,Int_t j) {
	if(GetpListItem(i,j)==0) return 0;
	return GetpListItem(i,j)->GetNsignals();
    }
    // Adds the contents of pl to the list with track number off set given by
    // fileIndex.
    virtual void AddItemTo(Int_t fileIndex, AliITSpListItem *pl);
    // Adds a Signal value to the map. Creating and expanding arrays as needed.
    void AddSignal(Int_t i,Int_t j,Int_t trk,Int_t ht,Int_t mod,Double_t sig);
    // Adds a Noise value to the map. Creating and expanding arrays as needed.
    void AddNoise(Int_t i,Int_t j,Int_t mod,Double_t noise);
    // Delete all AliITSpListItems and zero the TClonesArray
    virtual void ClearMap();
    // Delete a particular AliITSpListItem in the TClonesArray.
    virtual void DeleteHit(Int_t i,Int_t j);
    // returns hit index in TClonesArray
    virtual Int_t GetHitIndex(Int_t i,Int_t j) const {return GetIndex(i,j);}
    // returns "hit" AliITSpListItem as a TObject.
    TObject * GetHit(Int_t i,Int_t j){return (TObject*)GetpListItem(i,j);}
    // tests hit status.
    virtual FlagType TestHit(Int_t i,Int_t j){if(GetpListItem(i,j)==0) return kEmpty;
    else if(GetSignal(i,j)<=0) return kUnused; else return kUsed;}
    // Returns the pointer to the TClonesArray of pList Items
    TClonesArray * GetpListItems(){return fa;}
    // returns the pList Item stored in the TClonesArray
    AliITSpListItem* GetpListItem(Int_t index){
	if(fa!=0)return (AliITSpListItem*) (fa->At(index));
	else return 0;}
    // returns the pList Item stored in the TObject array
    AliITSpListItem* GetpListItem(Int_t i,Int_t j) const {
	if(fa!=0)return (AliITSpListItem*) (fa->At(GetIndex(i,j)));
	else return 0;}

    // Fill pList from digits. Not functional yet
    virtual void FillMap(){NotImplemented("FillMap");}
    // Sets threshold for significance. Not of relavance in this case.
    virtual void SetThreshold(Int_t /* i */){NotImplemented("SetThreshold");}
    // Sets a single hit. Not of relavance in this case.
    virtual void SetHit(Int_t /* i */,Int_t /* j */,Int_t /* k */){NotImplemented("SetHit");}
    // Flags a hit. Not of relavence in this case.
    virtual void FlagHit(Int_t /* i */,Int_t /* j */){NotImplemented("FlagHit");}
    virtual void GetCell(Int_t index,Int_t &i,Int_t &j) const;

 private:

// private methods
    Int_t GetIndex(Int_t i,Int_t j) const;
    void NotImplemented(const char *method) const {if(gDebug>0)
         Warning(method,"This method is not implemented for this class");}
// data members
    Int_t     fNi,fNj;   // The max index in i,j.
    TClonesArray *fa;       // array of pList items
    Int_t     fEntries; // keepts track of the number of non-zero entries.

    ClassDef(AliITSpList,4) // list of signals and track numbers
};	
#endif
