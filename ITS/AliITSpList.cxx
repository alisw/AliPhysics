/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/
/* $Id$ */

#include <stdio.h>
#include <stdlib.h>
#include <iostream.h>
#include <iomanip.h>
#include <TObjArray.h>
#include <TRandom.h>
#include <TMath.h>

#include "AliITSpList.h"

//______________________________________________________________________

ClassImp(AliITSpList);
//______________________________________________________________________
AliITSpList::AliITSpList(){
    // Default constructor

    fNi = 0;
    fNj = 0;
    fa  = 0;
}
//______________________________________________________________________
AliITSpList::AliITSpList(Int_t imax,Int_t jmax){
    // Standard constructor

    fNi = imax;
    fNj = jmax;
    fa  = new TObjArray(fNi*fNj); // elements are zeroed by 
                                  // TObjArray creator
}
//______________________________________________________________________
AliITSpList::~AliITSpList(){
    // Default destructor

    for(Int_t i=0;i<GetMaxIndex();i++) if(fa->At(i)!=0){
	delete fa->At(i);
	fa->AddAt(0,i); // zero content
    } // end for i && if
    fNi = 0;
    fNj = 0;
    delete fa;
    fa  = 0;
}
//______________________________________________________________________
void AliITSpList::ClearMap(){
    // Delete all AliITSpListItems and zero TObjArray.

    for(Int_t i=0;i<GetMaxIndex();i++) if(fa->At(i)!=0){
	delete fa->At(i);
	fa->AddAt(0,i); // zero content
    } // end for i && if
}
//______________________________________________________________________
void AliITSpList::DeleteHit(Int_t i,Int_t j){
    // Delete a particular AliITSpListItems and zero TObjArray.
    Int_t k = GetIndex(i,j);

    if(fa->At(k)!=0){
	delete fa->At(k);
	fa->AddAt(0,k); // zero content
    } // end for i && if
}
//______________________________________________________________________
AliITSpList& AliITSpList::operator=(const AliITSpList &source){
    // = operator

    if(this == &source) return *this;

    if(this->fa!=0){ // if this->fa exists delete it first.
	for(Int_t i=0;i<GetMaxIndex();i++) if(fa->At(i)!=0){
	    delete fa->At(i);
	    fa->AddAt(0,i); // zero content
	} // end for i && if
	delete this->fa;
    } // end if this->fa!=0
    this->fNi = source.fNi;
    this->fNj = source.fNj;
    this->fa = new TObjArray(*(source.fa));

    return *this;
}
//______________________________________________________________________
AliITSpList::AliITSpList(AliITSpList &source){
    // Copy operator

    *this = source;
}
//______________________________________________________________________
void AliITSpList::AddSignal(Int_t i,Int_t j,Int_t trk,Int_t ht,Int_t mod,
		       Double_t signal){
    // Adds a Signal value to the TObjArray at i,j. Creates the AliITSpListItem
    // if needed.

    if(GetpListItem(i,j)==0){ // most create AliITSpListItem
	fa->AddAt(new AliITSpListItem(trk,ht,mod,GetIndex(i,j),signal),
		  GetIndex(i,j));
    }else{ // AliITSpListItem exists, just add signal to it.
	GetpListItem(i,j)->AddSignal(trk,ht,mod,GetIndex(i,j),signal);
    } // end if
}
//______________________________________________________________________
void AliITSpList::AddNoise(Int_t i,Int_t j,Int_t mod,Double_t noise){
    // Adds a noise value to the TObjArray at i,j. Creates the AliITSpListItem
    // if needed.

    if(GetpListItem(i,j)==0){ // most create AliITSpListItem
	fa->AddAt(new AliITSpListItem(mod,GetIndex(i,j),noise),
		  GetIndex(i,j));
    }else{ // AliITSpListItem exists, just add signal to it.
	GetpListItem(i,j)->AddNoise(mod,GetIndex(i,j),noise);
    } // end if
}
//______________________________________________________________________

ClassImp(AliITSpListItem)
//______________________________________________________________________
AliITSpListItem::AliITSpListItem(){
    // Default constructor

    fmodule = -1;
    findex  = -1;
    for(Int_t i=0;i<this->fkSize;i++){
	this->fTrack[i]  = -2;
	this->fHits[i]   = -2;
	this->fSignal[i] = 0.0;
    } // end if i
    fTsignal = 0.0;
    fNoise   = 0.0;
}
//______________________________________________________________________
AliITSpListItem::AliITSpListItem(Int_t module,Int_t index,Double_t noise){
    // Standard noise constructor

    this->fmodule    = module;
    this->findex     = index;
    for(Int_t i=0;i<this->fkSize;i++){
	this->fTrack[i]  = -2;
	this->fSignal[i] = 0.0;
	this->fHits[i]   = 0;
    } // end if i
    this->fTsignal = 0.0;
    this->fNoise   = noise;
}
//______________________________________________________________________
AliITSpListItem::AliITSpListItem(Int_t track,Int_t hit,Int_t module,
			       Int_t index,Double_t signal){
    // Standard signal constructor

    this->fmodule    = module;
    this->findex     = index;
    this->fTrack[0]  = track;
    this->fHits[0]   = hit;
    this->fSignal[0] = signal;
    for(Int_t i=1;i<this->fkSize;i++){
	this->fTrack[i]  = -2;
	this->fSignal[i] = 0.0;
	this->fHits[i]   = 0;
    } // end if i
    this->fTsignal = signal;
    this->fNoise   = 0.0;
}
//______________________________________________________________________
AliITSpListItem::~AliITSpListItem(){
    // Denstructor

    this->fmodule = 0;
    this->findex  =0;
    for(Int_t i=0;i<=this->GetNsignals();i++){
	this->fTrack[i]  = 0;
	this->fSignal[i] = 0.0;
	this->fHits[i]   = 0;
    } // end if i
    this->fTsignal = 0.0;
    this->fNoise   =0.0;
}
//______________________________________________________________________
AliITSpListItem& AliITSpListItem::operator=(const AliITSpListItem &source){
    // = operator

    if(this == &source) return *this;

    this->fmodule = source.fmodule;
    this->findex  = source.findex;
    for(Int_t i=0;i<this->fkSize;i++){
	this->fTrack[i]  = source.fTrack[i];
	this->fSignal[i] = source.fSignal[i];
	this->fHits[i]   = source.fHits[i];
    } // end if i
    this->fTsignal = source.fTsignal;
    this->fNoise   = source.fNoise;

    return *this;
}
//______________________________________________________________________
AliITSpListItem::AliITSpListItem(AliITSpListItem &source){
    // Copy operator

    *this = source;
}
//______________________________________________________________________
void AliITSpListItem::AddSignal(Int_t track,Int_t hit,Int_t module,
			       Int_t index,Double_t signal){
    // Adds this track number and sinal to the pList and orders them
    Int_t    i,j,trk,hts;
    Double_t sig;
    Bool_t   flg=kFALSE;

    if(findex!=index || fmodule!=module) 
	Warning("AddSignal","index=%d != findex=%d or module=%d != fmodule=%d",
	    index,findex,module,fmodule);
    fTsignal += signal; // Keep track of sum signal.
    if(signal<=fSignal[fkSize-1]) return; // smaller than smallest
    for(i=0;i<fkSize;i++)if(track==fTrack[i] && hit ==fHits[i]){
	fSignal[i] += signal;
	flg = kTRUE;
    } // end for i & if.
    if(flg){ // resort arrays.  fkSize is small use Insertin sort.
	for(i=1;i<fkSize;i++){
	    trk = fTrack[i];
	    hts = fHits[i];
	    sig = fSignal[i];
	    j = i-1;
	    while(j>=0 && fSignal[i]>signal){
		fTrack[j+1]  = fTrack[j];
		fHits[j+1]   = fHits[j];
		fSignal[j+1] = fSignal[j];
		j--;
	    } // end while
	    fTrack[j+1]  = trk;
	    fHits[j+1]   = hts;
	    fSignal[j+1] = sig;
	} // end if i
    }else{ // new entry add it in order.
	if(!(signal <= fSignal[fkSize-1])) for(i=fkSize-2;i>=0;i--){
	    if(signal > fSignal[i]){
		fSignal[i+1] = fSignal[i];
		fTrack[i+1]  = fTrack[i];
		fHits[i+1]   = fHits[i];
	    }else{
		fSignal[i] = signal;
		fTrack[i]  = track;
		fHits[i]   = hit;
	    } //  end if
	} // end if; end for i
    } // end if flg
}
//______________________________________________________________________
void AliITSpListItem::AddNoise(Int_t module,Int_t index,Double_t noise){
    // Addes noise to this existing list.

    if(findex!=index || fmodule!=module) 
	Warning("AddSignal","index=%d != findex=%d or module=%d != fmodule=%d",
	    index,findex,module,fmodule);
    fNoise += noise; // Keep track of sum signal.
}
//______________________________________________________________________
void AliITSpListItem::AddTo(Int_t fileIndex,AliITSpListItem *pl){
    // Adds the contents of pl to this with track number off set given by
    // fileIndex.
    Int_t i,trk;
    Double_t sig=0.0;

    for(i=0;i<pl->GetNsignals()&&i<this->GetNsignals();i++){
	trk = pl->GetTrack(i);
	trk = pl->ShiftIndex(fileIndex,trk);
	this->AddSignal(trk,pl->GetHit(i),pl->GetModule(),pl->GetIndex(),pl->GetSignal(i));
	sig += pl->GetSignal(i);
    } // end for i
    this->fNoise   += pl->fNoise;
    return;
}
//______________________________________________________________________
Int_t AliITSpListItem::ShiftIndex(Int_t in,Int_t trk){
    // Shift an index number to occupy the upper four bits.
    Int_t si = sizeof(Int_t) * 8;
    UInt_t uin,utrk; // use UInt_t to avoid interger overflow-> goes negitive.

    uin = in;
    utrk = trk;
    for(Int_t i=0;i<si-4;i++) uin *= 2;
    uin += utrk;
    in = uin;
    return in;
}
//______________________________________________________________________
void AliITSpListItem::Print(ostream *os){
    //Standard output format for this class
    Int_t i;

    *os << fmodule <<","<<findex<<",";
    *os << fkSize <<",";
    for(i=0;i<fkSize;i++) *os << fTrack[i] <<",";
    for(i=0;i<fkSize;i++) *os << fHits[i] <<",";
    for(i=0;i<fkSize;i++) *os << fSignal[i] <<",";
    *os << fTsignal <<","<< fNoise;
}
//______________________________________________________________________
void AliITSpListItem::Read(istream *is){
    // Standard output streaming function.
    Int_t i,iss;

    *is >> fmodule >> findex;
    *is >> iss; // read in fkSize
    for(i=0;i<fkSize&&i<iss;i++) *is >> fTrack[i];
    for(i=0;i<fkSize&&i<iss;i++) *is >> fHits[i];
    for(i=0;i<fkSize&&i<iss;i++) *is >> fSignal[i];
    *is >> fTsignal >> fNoise;
}
//______________________________________________________________________
ostream &operator<<(ostream &os,AliITSpListItem &source){
    // Standard output streaming function.

    source.Print(&os);
    return os;
}
//______________________________________________________________________
istream &operator>>(istream &os,AliITSpListItem &source){
    // Standard output streaming function.

    source.Read(&os);
    return os;
}
