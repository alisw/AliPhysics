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

////////////////////////////////////////////////
//  Digits classes for all ITS detectors      //
////////////////////////////////////////////////
#include <TObjArray.h>
#include <TArrayI.h>
#include <TArrayF.h>
#include <TMath.h>
#include "AliITSdigit.h"

//______________________________________________________________________
ClassImp(AliITSdigit)
AliITSdigit::AliITSdigit(const Int_t *digits) {
  // Creates a real data digit object

  fCoord1       = digits[0];
  fCoord2       = digits[1];
  fSignal       = digits[2];
}
//______________________________________________________________________
void AliITSdigit::Print(ostream *os){
    //Standard output format for this class

    *os << fCoord1 <<","<< fCoord2 <<","<< fSignal;
}
//______________________________________________________________________
void AliITSdigit::Read(istream *os){
    //Standard input for this class

    *os >> fCoord1 >> fCoord2 >> fSignal;
}
//______________________________________________________________________
ostream &operator<<(ostream &os,AliITSdigit &source){
    // Standard output streaming function.

    source.Print(&os);
    return os;
}
//______________________________________________________________________
istream &operator>>(istream &os,AliITSdigit &source){
    // Standard output streaming function.

    source.Read(&os);
    return os;
}
//______________________________________________________________________
ClassImp(AliITSdigitSPD)
AliITSdigitSPD::AliITSdigitSPD():AliITSdigit(){
    // default constructor, zero coordinates and set array
    // elements to clearly unphysical values. A value of 0 may
    // be a valide track of hit number.
    Int_t i;

    for(i=0;i<fkSspd;i++) fTracks[i]  = -3;
    for(i=0;i<fkSspd;i++) fHits[i]    = -1;
}
//______________________________________________________________________
AliITSdigitSPD::AliITSdigitSPD(const Int_t *digits){
    // Creates a SPD digit object
    Int_t i;

    for(i=0;i<fkSspd;i++) fTracks[i]  = -3;
    for(i=0;i<fkSspd;i++) fHits[i]    = -1;
    fCoord1       = digits[0];
    fCoord2       = digits[1];
    fSignal       = 1;
    fSignalSPD    = digits[2];
}
//______________________________________________________________________
AliITSdigitSPD::AliITSdigitSPD(const Int_t *digits,const Int_t *tracks,
			       const Int_t *hits){
    // Creates a simulated SPD digit object

    for(Int_t i=0; i<fkSspd; i++) {
	fTracks[i] = tracks[i];
	fHits[i]   = hits[i];
    } // end for i
    fCoord1       = digits[0];
    fCoord2       = digits[1];
    fSignal       = 1;
    fSignalSPD    = digits[2];
}
//______________________________________________________________________
Int_t AliITSdigitSPD::GetListOfTracks(TArrayI &t){
    // Fills the TArrayI t with the tracks found in fTracks removing
    // duplicated tracks, but otherwise in the same order. It will return
    // the number of tracks and fill the remaining elements to the array
    // t with -1.
    // Inputs:
    //   TArrayI  &t Reference to a TArrayI to contain the list of
    //               nonduplicated track numbers.
    // Output:
    //   TArrayI  &t The input array filled with the nonduplicated track
    //               numbers.
    // Return:
    //   Int_t The number of none -1 entries in the TArrayI t.
    Int_t nt = t.GetSize();
    Int_t nth = this->GetNTracks();
    Int_t n = 0,i,j;
    Bool_t inlist = kFALSE;

    t.Reset(-1); // -1 array.
    for(i=0;i<nth;i++) {
	if(this->GetTrack(i) == -1) continue;
	inlist = kFALSE;
	for(j=0;j<n;j++)if(this->GetTrack(i) == t.At(j)) inlist = kTRUE;
	if(!inlist){ // add to end of list
	    t.AddAt(this->GetTrack(i),n);
	    if(n<nt) n++;
	} // end if
    } // end for i
    return n;
}
//______________________________________________________________________
void AliITSdigitSPD::Print(ostream *os){
    //Standard output format for this class
    Int_t i;

    AliITSdigit::Print(os);
    for(i=0;i<fkSspd;i++) *os <<","<< fTracks[i];
    for(i=0;i<fkSspd;i++) *os <<","<< fHits[i];
    *os << "," << fSignalSPD;
}
//______________________________________________________________________
void AliITSdigitSPD::Read(istream *os){
    //Standard input for this class
    Int_t i;

    AliITSdigit::Read(os);
    for(i=0;i<fkSspd;i++) *os >> fTracks[i];
    for(i=0;i<fkSspd;i++) *os >> fHits[i];
    *os >> fSignalSPD;
}
//______________________________________________________________________
ostream &operator<<(ostream &os,AliITSdigitSPD &source){
    // Standard output streaming function.

    source.Print(&os);
    return os;
}
//______________________________________________________________________
istream &operator>>(istream &os,AliITSdigitSPD &source){
    // Standard output streaming function.

    source.Read(&os);
    return os;
}
//______________________________________________________________________
ClassImp(AliITSdigitSDD)
AliITSdigitSDD::AliITSdigitSDD():AliITSdigit(){
    // default constructor, zero coordinates and set array
    // elements to clearly unphysical values. A value of 0 may
    // be a valide track of hit number.
    Int_t i;

    for(i=0;i<fkSsdd;i++) fTracks[i] = -3;
    for(i=0;i<fkSsdd;i++) fHits[i]   = -1;
    fPhysics = 0;
    for(i=0;i<fkSsdd;i++) fTcharges[i] = 0;
}
//________________________________________________________________________
AliITSdigitSDD::AliITSdigitSDD(Float_t phys,const Int_t *digits):
    AliITSdigit(digits){
    // Creates a simulated SDD digit object to be updated

    fPhysics = phys;
}
//_____________________________________________________________________________
AliITSdigitSDD::AliITSdigitSDD(Float_t phys,const Int_t *digits,
			       const Int_t *tracks,const Int_t *hits,
			       const Float_t *charges):
    AliITSdigit(digits){
    // Creates a simulated SDD digit object

    fPhysics = phys;
    for(Int_t i=0; i<fkSsdd; i++) {
	fTcharges[i] = charges[i];
	fTracks[i]   = tracks[i];
	fHits[i]     = hits[i];
    } // end for i
}
//______________________________________________________________________
Int_t AliITSdigitSDD::GetListOfTracks(TArrayI &t,TArrayF &c){
    // Fills the TArrayI t with the tracks found in fTracks removing
    // duplicated tracks, summing up their charge, and ordering the tracks
    // by the charge contributed to this digit. It will return
    // the number of tracks and fill the remaining elements to the array
    // t with -1.
    // Inputs:
    //   TArrayI  &t Reference to a TArrayI to contain the list of
    //               nonduplicated track numbers.
    //   TArrayF  &c Reference to a TArrayF to contain the summed charge
    //               contributed by each track.
    // Output:
    //   TArrayI  &t The input array filled with the nonduplicated track
    //               numbers.
    //   TArrayF  &c The input array filled with the summed charge 
    //               contributed by the corresponding track in the array t.
    // Return:
    //   Int_t The number of none -1 entries in the TArrayI t.
    Int_t nt = t.GetSize();
    nt = TMath::Min(nt,c.GetSize());
    Int_t nth = this->GetNTracks();
    Int_t n = 0,i,j;
    Bool_t inlist = kFALSE;

    t.Reset(-1); // -1 array.
    c.Reset(0.0); // zero array.
    for(i=0;i<nth;i++) {
	if(this->GetTrack(i) == -1) continue;
	inlist = kFALSE;
	for(j=0;j<n;j++)if(this->GetTrack(i) == t.At(j)){
	    inlist = kTRUE;
	    c.AddAt(this->GetCharge(i)+c.At(j),j);
	} // end for j/end if
	if(!inlist){ // add to end of list
	    t.AddAt(this->GetTrack(i),n);
	    c.AddAt(this->GetCharge(i),n);
	    if(n<nt) n++;
	} // end if
    } // end for i

    // Now lets sort the TArrays according to the charge. This algorithm
    // is based on the method from Chapter 8 section 1 Straight Insertion
    // sort. Wiliam H. Press, Saul A. Teukolsky, William T. Vetterling
    // and Brian P. Flannery, "Numerical Recipeis in C, The Art of Scientific
    // Computing", second Edition page 330 (1997).
    Int_t   tr;
    Float_t ch;
    for(i=0;i<n;i++){
	tr = t.At(i);
	ch = c.At(i);
	j = i-1;
	while(j>-1 && c.At(j)>ch){
	    t.AddAt(t.At(j+1),j);
	    c.AddAt(c.At(j+1),j);
	    j--;
	} // end while
	t.AddAt(tr,j+1);
	c.AddAt(ch,j+1);
    } // end for i
    //
    return n;
}
//______________________________________________________________________
void AliITSdigitSDD::Print(ostream *os){
    //Standard output format for this class
    Int_t i;

    AliITSdigit::Print(os);
    *os <<","<< fPhysics;
    for(i=0; i<fkSsdd; i++) *os <<","<< fTcharges[i];
    for(i=0; i<fkSsdd; i++) *os <<","<< fTracks[i];
    for(i=0; i<fkSsdd; i++) *os <<","<< fHits[i];
}
//______________________________________________________________________
void AliITSdigitSDD::Read(istream *os){
    //Standard input for this class
    Int_t i;

    AliITSdigit::Read(os);
    *os >>fPhysics;
    for(i=0; i<fkSsdd; i++) *os >> fTcharges[i];
    for(i=0; i<fkSsdd; i++) *os >> fTracks[i];
    for(i=0; i<fkSsdd; i++) *os >> fHits[i];
}
//______________________________________________________________________
ostream &operator<<(ostream &os,AliITSdigitSDD &source){
    // Standard output streaming function.

    source.Print(&os);
    return os;
}
//______________________________________________________________________
istream &operator>>(istream &os,AliITSdigitSDD &source){
    // Standard output streaming function.

    source.Read(&os);
    return os;
}
//______________________________________________________________________
ClassImp(AliITSTransientDigit)
AliITSTransientDigit::AliITSTransientDigit(Float_t phys,const Int_t *digits): 
    AliITSdigitSDD(phys,digits) {
    // Creates a digit object in a list of digits to be updated

    fTrackList   = new TObjArray;  
}
//__________________________________________________________________________
AliITSTransientDigit::AliITSTransientDigit(const AliITSTransientDigit &source):
 AliITSdigitSDD(source){
    // Copy Constructor 
    if(&source == this) return;
    this->fTrackList = source.fTrackList;
    return;
}
//_________________________________________________________________________
AliITSTransientDigit& AliITSTransientDigit::operator=(
    const AliITSTransientDigit &source) {
    // Assignment operator
    if(&source == this) return *this;
    this->fTrackList = source.fTrackList;
    return *this;
}
//______________________________________________________________________
void AliITSTransientDigit::Print(ostream *os){
    //Standard output format for this class

    AliITSdigitSDD::Print(os);
}
//______________________________________________________________________
void AliITSTransientDigit::Read(istream *os){
    //Standard input for this class

    AliITSdigitSDD::Read(os);
}
//______________________________________________________________________
ostream &operator<<(ostream &os,AliITSTransientDigit &source){
    // Standard output streaming function.

    source.Print(&os);
    return os;
}
//______________________________________________________________________
istream &operator>>(istream &os,AliITSTransientDigit &source){
    // Standard output streaming function.

    source.Read(&os);
    return os;
}
//______________________________________________________________________
ClassImp(AliITSdigitSSD)
AliITSdigitSSD::AliITSdigitSSD():AliITSdigit(){
    // default constructor
    Int_t i;

    for(i=0; i<fkSssd; i++) fTracks[i] = -3;
    for(i=0; i<fkSssd; i++) fHits[i] = -1;
}
//__________________________________________________________________________
AliITSdigitSSD::AliITSdigitSSD(const Int_t *digits):AliITSdigit(digits){
    // Creates a real SSD digit object
}
//_____________________________________________________________________________
AliITSdigitSSD::AliITSdigitSSD(const Int_t *digits,const Int_t *tracks,
			       const Int_t *hits):AliITSdigit(digits){
    // Creates a simulated SSD digit object

    for(Int_t i=0; i<fkSssd; i++) {
	fTracks[i] = tracks[i];
	fHits[i]   = hits[i];
    } // end for i
}
//______________________________________________________________________
Int_t AliITSdigitSSD::GetListOfTracks(TArrayI &t){
    // Fills the TArrayI t with the tracks found in fTracks removing
    // duplicated tracks, but otherwise in the same order. It will return
    // the number of tracks and fill the remaining elements to the array
    // t with -1.
    // Inputs:
    //   TArrayI  &t Reference to a TArrayI to contain the list of
    //               nonduplicated track numbers.
    // Output:
    //   TArrayI  &t The input array filled with the nonduplicated track
    //               numbers.
    // Return:
    //   Int_t The number of none -1 entries in the TArrayI t.
    Int_t nt = t.GetSize();
    Int_t nth = this->GetNTracks();
    Int_t n = 0,i,j;
    Bool_t inlist = kFALSE;

    t.Reset(-1); // -1 array.
    for(i=0;i<nth;i++) {
	if(this->GetTrack(i) == -1) continue;
	inlist = kFALSE;
	for(j=0;j<n;j++)if(this->GetTrack(i) == t.At(j)) inlist = kTRUE;
	if(!inlist){ // add to end of list
	    t.AddAt(this->GetTrack(i),n);
	    if(n<nt) n++;
	} // end if
    } // end for i
    return n;
}
//______________________________________________________________________
void AliITSdigitSSD::Print(ostream *os){
    //Standard output format for this class
    Int_t i;

    AliITSdigit::Print(os);
    for(i=0; i<fkSssd; i++) *os <<","<< fTracks[i];
    for(i=0; i<fkSssd; i++) *os <<","<< fHits[i];
}
//______________________________________________________________________
void AliITSdigitSSD::Read(istream *os){
    //Standard input for this class
    Int_t i;

    AliITSdigit::Read(os);
    for(i=0; i<fkSssd; i++) *os >> fTracks[i];
    for(i=0; i<fkSssd; i++) *os >> fHits[i];
}
//______________________________________________________________________
ostream &operator<<(ostream &os,AliITSdigitSSD &source){
    // Standard output streaming function.

    source.Print(&os);
    return os;
}
//______________________________________________________________________
istream &operator>>(istream &os,AliITSdigitSSD &source){
    // Standard output streaming function.

    source.Read(&os);
    return os;
}
