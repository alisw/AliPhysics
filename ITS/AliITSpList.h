#ifndef ALIITSPLISTSSDITEM_H
#define ALIITSPLISTSSDITEM_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */

/* $Id$ */

#include <TObject.h>


class AliITSpListItem: public TObject {
 public:
    // Default Constructor
    AliITSpListItem();
    // Standard Signal Constructor
    AliITSpListItem(Int_t track,Int_t hit,Int_t module,Int_t index,
		   Double_t signal);
    // Standard Noise Constructor
    AliITSpListItem(Int_t module,Int_t index,Double_t signal);
    // Class destrutor
    virtual ~AliITSpListItem();
    // Copy Oporator
    AliITSpListItem(AliITSpListItem &source);
    // = Opoerator
    virtual AliITSpListItem& operator=(const AliITSpListItem &source);
    // Returns the signal value in the list of signals
    virtual Double_t GetSignal(Int_t i){
	                    return ( (i>=0&&i<fkSize) ? fSignal[i] : 0.0);}
    virtual Double_t GetSignal(){
	                    return fTsignal;}
    virtual Double_t GetSignalAfterElect(){
	                    return fSignalAfterElect;}
    // Returns the Sum/Total signal
    virtual Double_t GetSumSignal() const {return fTsignal+fNoise;}
    // Returns the  noise
    virtual Double_t GetNoise() const {return fNoise;}
    // Returns the number of stored singals.
    virtual Int_t GetNsignals() const {return fkSize;}
    // Addes track number and signal to this existing list.
    virtual void AddSignal(Int_t track,Int_t hit,Int_t module,
			   Int_t index,Double_t signal);
    // Adds signal after electronics to this existing list.
    virtual void AddSignalAfterElect(Int_t module,Int_t index,Double_t signal);
    // Addes noise to this existing list.
    virtual void AddNoise(Int_t module,Int_t index,Double_t noise);
    // Returns track number.
    virtual Int_t GetTrack(Int_t i){
	                    return ((i>=0&&i<fkSize) ? fTrack[i] : 0);}
    // Returns hit number.
    virtual Int_t GetHit(Int_t i){
	                    return ((i>=0&&i<fkSize) ? fHits[i] : 0);}
    // Returns module number.
    virtual Int_t GetModule(){
	                    return fmodule;}
    // Returns index number.
    virtual Int_t GetIndex(){
	                    return findex;}
    // Adds the contents of pl to this 
    virtual void Add(AliITSpListItem *pl);
    // Adds the contents of pl to this with track number off set given by
    // fileIndex.
    virtual void AddTo(Int_t fileIndex,AliITSpListItem *pl);
    // Shift an index number to occupy the upper four bits.
    virtual Int_t ShiftIndex(Int_t in,Int_t trk);
    // Standard ascii class print function
    void Print(ostream *os);
    // Standard ascii class read function
    void Read(istream *is);
    // Returns max size of array for for Tracks, Hits, and signals.
    static Int_t GetMaxKept() {return fkSize;};

 private:
    static const Int_t fkSize = 10; // Array sizes
    Int_t    fmodule;         // module number
    Int_t    findex;          // Strip/row,col number linearlized.
    Int_t    fTrack[fkSize];  //[fkSize] track Number
    Int_t    fHits[fkSize];   //[fkSize] hit number
    Double_t fSignal[fkSize]; //[fkSize] Signals
    Double_t fTsignal;        // Total signal (no noise)
    Double_t fNoise;          // Total noise, coupling, ...
    Double_t fSignalAfterElect; // Signal after electronics

    ClassDef(AliITSpListItem,3) // Item list of signals and track numbers
};	
// Input and output functions for standard C++ input/output.
ostream & operator<<(ostream &os,AliITSpListItem &source);
istream & operator>>(istream &is,AliITSpListItem &source);


#endif

#ifndef ALIITSPLISTSSD_H
#define ALIITSPLISTSSD_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */

#include "AliITSMap.h"

class TObjArray;
class AliITSpListItem;

class AliITSpList: public AliITSMap {
 private:
    // returns the TObjArray index for a give set of map indecies.
    Int_t GetIndex(Int_t i,Int_t j){
	if(i<0||i>=fNi || j<0||j>=fNj){
	    Warning("GetIndex","Index out of range 0<i=%d<%d and 0<0j=%d<%d",
		    i,fNi,j,fNj);
	    return -1;
	} // end if
	return fNj*i+j;
    }

 public:
    // Default Constructor
    AliITSpList();
    // Standard Constructor
    AliITSpList(Int_t imax,Int_t jmax);
    // Class destrutor
    virtual ~AliITSpList();
    // Copy Oporator
    AliITSpList(AliITSpList &source);
    // = Opoerator
    virtual AliITSpList& operator=(const AliITSpList &source);
    // Returns the max mape index values
    void GetMaxMapIndex(Int_t &ni,Int_t &nj){ni=fNi;nj=fNj;return;}
    // returns the max index value.
    Int_t GetMaxIndex(){return fNi*fNj;}
    // returns the largest non-zero entry kept in the array fa.
    Int_t GetEntries(){return fEnteries;}
    // returns the max number of track/hit entries per cell.
    Int_t GetNEnteries(){return AliITSpListItem::GetMaxKept();}
    // for a give TObjArray index it returns the corresponding map index
    void  GetMapIndex(Int_t index,Int_t &i,Int_t &j){
	i = index/fNj;j = index - fNj*i;
	if(i<0||i>=fNi || j<0||j>=fNj){i=-1;j=-1; return;}
    }
    // Returns the signal+noise for a give map coordinate
    Double_t GetSignal(Int_t index){
	if(GetpListItem(index)==0) return 0.0;
	return GetpListItem(index)->GetSumSignal();
    }
    // Returns the signal+noise for a give map coordinate
    Double_t GetSignal(Int_t i,Int_t j){
	if(GetpListItem(i,j)==0) return 0.0;
	return GetpListItem(i,j)->GetSumSignal();
    }
    // Returns the signal only for a give map coordinate
    Double_t GetSignalOnly(Int_t i,Int_t j){
	if(GetpListItem(i,j)==0) return 0.0;
	return GetpListItem(i,j)->GetSignal();
    }
    // Returns the noise for a give map coordinate
    Double_t GetNoise(Int_t i,Int_t j){
	if(GetpListItem(i,j)==0) return 0.0;
	return GetpListItem(i,j)->GetNoise();
    }
    // returns the track number which generated the signal at a given map
    // coordinate. If there is no signal or only noise, then -2 is returned.
    // k is the track rank number.
    Double_t GetTSignal(Int_t i,Int_t j,Int_t k){
	if(GetpListItem(i,j)==0) return 0.0;
	return GetpListItem(i,j)->GetSignal(k);
    }
    // returns the track number which generated the signal at a given map
    // coordinate. If there is no signal or only noise, then -2 is returned.
    // k is the track rank number.
    Int_t GetTrack(Int_t i,Int_t j,Int_t k){
	if(GetpListItem(i,j)==0) return -2;
	return GetpListItem(i,j)->GetTrack(k);
    }
    // returns the hit number which generated the signal at a given map
    // coordinate. If there is no signal or only noise, then -2 is returned.
    // k is the hit rank number.
    Int_t GetHit(Int_t i,Int_t j,Int_t k){
	if(GetpListItem(i,j)==0) return -2;
	return GetpListItem(i,j)->GetHit(k);
    }
    // returns the number of Signal values
    Int_t GetNSignals(Int_t i,Int_t j){
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
    // Delete all AliITSpListItems and zero the TObjArray
    void ClearMap();
    // Delete a particular AliITSpListItem in the TObjArray.
    void DeleteHit(Int_t i,Int_t j);
    // returns hit index in TObjArray
    Int_t GetHitIndex(Int_t i,Int_t j){return GetIndex(i,j);}
    // returns "hit" AliITSpListItem as a TObject.
    TObject * GetHit(Int_t i,Int_t j){return (TObject*)GetpListItem(i,j);}
    // tests hit status.
    FlagType TestHit(Int_t i,Int_t j){if(GetpListItem(i,j)==0) return kEmpty;
    else if(GetSignal(i,j)<=0) return kUnused; else return kUsed;}
    // Returns the pointer to the TObjArray of pList Items
    TObjArray * GetpListItems(){return fa;}
    // returns the pList Item stored in the TObject array
    AliITSpListItem* GetpListItem(Int_t index){
	if(fa!=0)return (AliITSpListItem*) (fa->At(index));
	else return 0;}
    // returns the pList Item stored in the TObject array
    AliITSpListItem* GetpListItem(Int_t i,Int_t j){
	if(fa!=0)return (AliITSpListItem*) (fa->At(GetIndex(i,j)));
	else return 0;}

    // Fill pList from digits. Not functional yet
    void FillMap(){;}
    // Sets threshold for significance. Not of relavance in this case.
    void SetThreshold(Int_t i){i=0;}
    // Sets a single hit. Not of relavance in this case.
    void SetHit(Int_t i,Int_t j,Int_t k){i=j=k=0;}
    // Flags a hit. Not of relavence in this case.
    void FlagHit(Int_t i,Int_t j){i=j=0;}
    // returns the i,j index numbers from the liniarized index computed
    // from GetIndex above.
    void GetCell(Int_t index,Int_t &i,Int_t &j){
	if(index<0 || index>=fNi*fNj){
	    Warning("GetCell","Index out of range 0<=index=%d<%d",
		    index,fNi*fNj);
	    i=-1;j=-1;
	    return;
	} // end if
	i = index/fNj;
	j = index - fNj*i;
	return;
    }

 private:
    Int_t     fNi,fNj;   // The max index in i,j.
    TObjArray *fa;       // array of pList items
    Int_t     fEnteries; // keepts track of the number of non-zero entries.

    ClassDef(AliITSpList,3) // list of signals and track numbers
};	
// Input and output functions for standard C++ input/output.
ostream & operator<<(ostream &os,AliITSpList &source);
istream & operator>>(istream &is,AliITSpList &source);

#endif
