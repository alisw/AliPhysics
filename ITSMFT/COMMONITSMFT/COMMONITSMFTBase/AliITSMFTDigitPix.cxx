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
#include <Riostream.h>
#include <TArrayI.h>

#include <AliITSMFTDigitPix.h>

///////////////////////////////////////////////////////////////////
//                                                               //
// Class defining the digit object
// for pixel
//                                                               //
///////////////////////////////////////////////////////////////////

ClassImp(AliITSMFTDigitPix)

//______________________________________________________________________
AliITSMFTDigitPix::AliITSMFTDigitPix()
:  TObject()
  ,fCoord1(0)
  ,fCoord2(0)
  ,fSignal(0)
  ,fSignalPix(0)
  ,fROCycle(0)
{
    // default constructor, zero coordinates and set array
    // elements to clearly unphysical values. A value of 0 may
    // be a valide track of hit number.
    Int_t i;

    for(i=0;i<fgkSize;i++) fTracks[i]  = -3;
    for(i=0;i<fgkSize;i++) fHits[i]    = -1;
}

//______________________________________________________________________
AliITSMFTDigitPix::AliITSMFTDigitPix(const Int_t *digits)
:  TObject()
  ,fCoord1(0)
  ,fCoord2(0)
  ,fSignal(0)
  ,fSignalPix(digits[2])
  ,fROCycle(digits[3])
{
    // Creates a pixel digit object
    Int_t i;

    for(i=0;i<fgkSize;i++) fTracks[i]  = -3;
    for(i=0;i<fgkSize;i++) fHits[i]    = -1;
    fCoord1       = digits[0];
    fCoord2       = digits[1];
    fSignal       = 1;
}

//______________________________________________________________________
AliITSMFTDigitPix::AliITSMFTDigitPix(const Int_t *digits,const Int_t *tracks,const Int_t *hits)
: TObject()
 ,fCoord1(0)
 ,fCoord2(0)
 ,fSignal(0)
 ,fSignalPix(digits[2])
 ,fROCycle(digits[3])
{
    // Creates a simulated pixel digit object

    for(Int_t i=0; i<fgkSize; i++) {
	fTracks[i] = tracks[i];
	fHits[i]   = hits[i];
    } // end for i
    fCoord1       = digits[0];
    fCoord2       = digits[1];
    fSignal       = 1;
}
//______________________________________________________________________

AliITSMFTDigitPix::AliITSMFTDigitPix(const AliITSMFTDigitPix &h)
: TObject(h)
, fCoord1(h.fCoord1)
, fCoord2(h.fCoord2)
, fSignal(h.fSignal)
, fSignalPix(h.fSignalPix)
, fROCycle(h.fROCycle)

{
    for ( int i = 0 ; i < fgkSize; i++)
    {
        fTracks[i] = h.fTracks[i];
        fHits[i] = h.fHits[i];
    }
  
    
}
//______________________________________________________________________
AliITSMFTDigitPix& AliITSMFTDigitPix::operator=(const AliITSMFTDigitPix& source)
{
    // Assignment operator.
    if (this==&source) return *this;
    this->~AliITSMFTDigitPix();
    new(this) AliITSMFTDigitPix(source);
    return *this;
    //
}
//______________________________________________________________________
AliITSMFTDigitPix::~AliITSMFTDigitPix()
{
    // dtor

}
//______________________________________________________________________
Int_t AliITSMFTDigitPix::GetListOfTracks(TArrayI &t)
{
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
void AliITSMFTDigitPix::Print(std::ostream *os) const
{
    //Standard output format for this class
    Int_t i;

    *os << fCoord1 <<","<< fCoord2 <<","<< fSignal;
    for(i=0;i<fgkSize;i++) *os <<","<< fTracks[i];
    for(i=0;i<fgkSize;i++) *os <<","<< fHits[i];
    *os << "," << fSignalPix;
    *os << "," << fROCycle;
}

//______________________________________________________________________
void AliITSMFTDigitPix::Read(std::istream *os)
{
    //Standard input for this class
    Int_t i;

    *os >> fCoord1 >> fCoord2 >> fSignal;
    for(i=0;i<fgkSize;i++) *os >> fTracks[i];
    for(i=0;i<fgkSize;i++) *os >> fHits[i];
    *os >> fSignalPix;
    *os >> fROCycle;
}

//______________________________________________________________________
std::ostream &operator<<(std::ostream &os,AliITSMFTDigitPix &source)
{
    // Standard output streaming function.

    source.Print(&os);
    return os;
}

//______________________________________________________________________
std::istream &operator>>(std::istream &os,AliITSMFTDigitPix &source)
{
    // Standard output streaming function.

    source.Read(&os);
    return os;
}
