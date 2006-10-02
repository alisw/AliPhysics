/**************************************************************************
 * Copyright(c) 2004-2006, ALICE Experiment at CERN, All rights reserved. *
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

#include <AliITSdigitSPD.h>
#include <TArrayI.h>

///////////////////////////////////////////////////////////////////
//                                                               //
// Class defining the digit object
// for SPD
// Inherits from AliITSdigit
//                                                               //
///////////////////////////////////////////////////////////////////

ClassImp(AliITSdigitSPD)

//______________________________________________________________________
AliITSdigitSPD::AliITSdigitSPD():AliITSdigit(),
fSignalSPD(0){
    // default constructor, zero coordinates and set array
    // elements to clearly unphysical values. A value of 0 may
    // be a valide track of hit number.
    Int_t i;

    for(i=0;i<fgkSspd;i++) fTracks[i]  = -3;
    for(i=0;i<fgkSspd;i++) fHits[i]    = -1;
}
//______________________________________________________________________
AliITSdigitSPD::AliITSdigitSPD(const Int_t *digits):
fSignalSPD(digits[2]){
    // Creates a SPD digit object
    Int_t i;

    for(i=0;i<fgkSspd;i++) fTracks[i]  = -3;
    for(i=0;i<fgkSspd;i++) fHits[i]    = -1;
    fCoord1       = digits[0];
    fCoord2       = digits[1];
    fSignal       = 1;
 }
//______________________________________________________________________
AliITSdigitSPD::AliITSdigitSPD(const Int_t *digits,const Int_t *tracks,
			       const Int_t *hits):
fSignalSPD(digits[2]){
    // Creates a simulated SPD digit object

    for(Int_t i=0; i<fgkSspd; i++) {
	fTracks[i] = tracks[i];
	fHits[i]   = hits[i];
    } // end for i
    fCoord1       = digits[0];
    fCoord2       = digits[1];
    fSignal       = 1;
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
    for(i=0;i<fgkSspd;i++) *os <<","<< fTracks[i];
    for(i=0;i<fgkSspd;i++) *os <<","<< fHits[i];
    *os << "," << fSignalSPD;
}
//______________________________________________________________________
void AliITSdigitSPD::Read(istream *os){
    //Standard input for this class
    Int_t i;

    AliITSdigit::Read(os);
    for(i=0;i<fgkSspd;i++) *os >> fTracks[i];
    for(i=0;i<fgkSspd;i++) *os >> fHits[i];
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
