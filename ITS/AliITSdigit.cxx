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

    fTracks[0] = fTracks[1] = fTracks[2] = -3;
    fHits[0] = fHits[1] = fHits[2] = -1;
}
//______________________________________________________________________
AliITSdigitSPD::AliITSdigitSPD(const Int_t *digits){
    // Creates a SPD digit object

    fTracks[0] = fTracks[1] = fTracks[2] = -3;
    fHits[0] = fHits[1] = fHits[2] = -1;
    fCoord1       = digits[0];
    fCoord2       = digits[1];
    fSignal       = 1;
    fSignalSPD    = digits[2];
}
//______________________________________________________________________
AliITSdigitSPD::AliITSdigitSPD(const Int_t *digits,const Int_t *tracks,
			       const Int_t *hits){
    // Creates a simulated SPD digit object

    for(Int_t i=0; i<3; i++) {
	fTracks[i] = tracks[i];
	fHits[i]   = hits[i];
    } // end for i
    fCoord1       = digits[0];
    fCoord2       = digits[1];
    fSignal       = 1;
    fSignalSPD    = digits[2];
}
//______________________________________________________________________
void AliITSdigitSPD::Print(ostream *os){
    //Standard output format for this class

    AliITSdigit::Print(os);
    *os <<","<< fTracks[0] <<","<< fTracks[1] <<","<< fTracks[2];
    *os <<","<< fHits[0] <<","<< fHits[1] <<","<< fHits[2];
    *os << "," << fSignalSPD;
}
//______________________________________________________________________
void AliITSdigitSPD::Read(istream *os){
    //Standard input for this class

    AliITSdigit::Read(os);
    *os >> fTracks[0] >> fTracks[1] >> fTracks[2];
    *os >> fHits[0] >> fHits[1] >> fHits[2] >> fSignalSPD;
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

    fTracks[0] = fTracks[1] = fTracks[2] = -3;
    fHits[0] = fHits[1] = fHits[2] = -1;
    fPhysics = 0;
    fTcharges[0] = fTcharges[1] = fTcharges[2] = 0;
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
    for(Int_t i=0; i<3; i++) {
	fTcharges[i] = charges[i];
	fTracks[i]   = tracks[i];
	fHits[i]     = hits[i];
    } // end for i
}
//______________________________________________________________________
void AliITSdigitSDD::Print(ostream *os){
    //Standard output format for this class

    AliITSdigit::Print(os);
    *os <<","<< fPhysics;
    *os <<","<< fTcharges[0] <<","<< fTcharges[1] <<","<< fTcharges[2];
    *os <<","<< fTracks[0] <<","<< fTracks[1] <<","<< fTracks[2];
    *os <<","<< fHits[0] <<","<< fHits[1] <<","<< fHits[2];
}
//______________________________________________________________________
void AliITSdigitSDD::Read(istream *os){
    //Standard input for this class

    AliITSdigit::Read(os);
    *os >>fPhysics;
    *os >> fTcharges[0] >> fTcharges[1] >> fTcharges[2];
    *os >> fTracks[0] >> fTracks[1] >> fTracks[2];
    *os >> fHits[0] >> fHits[1] >> fHits[2];
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
AliITSTransientDigit::AliITSTransientDigit(const AliITSTransientDigit &source){
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

    fTracks[0] = fTracks[1] = fTracks[2] = -3;
    fHits[0] = fHits[1] = fHits[2] = -1;
}
//__________________________________________________________________________
AliITSdigitSSD::AliITSdigitSSD(const Int_t *digits):AliITSdigit(digits){
    // Creates a real SSD digit object
}
//_____________________________________________________________________________
AliITSdigitSSD::AliITSdigitSSD(const Int_t *digits,const Int_t *tracks,
			       const Int_t *hits):AliITSdigit(digits){
    // Creates a simulated SSD digit object

    for(Int_t i=0; i<3; i++) {
	fTracks[i] = tracks[i];
	fHits[i]   = hits[i];
    } // end for i
}
//______________________________________________________________________
void AliITSdigitSSD::Print(ostream *os){
    //Standard output format for this class

    AliITSdigit::Print(os);
    *os <<","<< fTracks[0] <<","<< fTracks[1] <<","<< fTracks[2];
    *os <<","<< fHits[0] <<","<< fHits[1] <<","<< fHits[2];
}
//______________________________________________________________________
void AliITSdigitSSD::Read(istream *os){
    //Standard input for this class

    AliITSdigit::Read(os);
    *os >> fTracks[0] >> fTracks[1] >> fTracks[2];
    *os >> fHits[0] >> fHits[1] >> fHits[2];
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
