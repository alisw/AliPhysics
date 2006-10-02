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

#include <TObjArray.h>
#include "AliITSTransientDigit.h"

///////////////////////////////////////////////////////////////////
//                                                               //
// Class used internally by AliITSsimulationSDD
// for SDD digitisation
// It is not currently used any longer
// The methods in ALiITSsimulationSDD using it are currently commented out
//                                                               //
///////////////////////////////////////////////////////////////////

ClassImp(AliITSTransientDigit)

//______________________________________________________________________
AliITSTransientDigit::AliITSTransientDigit(Float_t phys,const Int_t *digits): 
    AliITSdigitSDD(phys,digits),
fTrackList(0) {
    // Creates a digit object in a list of digits to be updated

    fTrackList   = new TObjArray;  
}
//__________________________________________________________________________
AliITSTransientDigit::AliITSTransientDigit(const AliITSTransientDigit &source):
 AliITSdigitSDD(source),
fTrackList(source.fTrackList){
    // Copy Constructor 
}
//_________________________________________________________________________
AliITSTransientDigit& AliITSTransientDigit::operator=(
    const AliITSTransientDigit &source) {
    // Assignment operator
  this->~AliITSTransientDigit();
  new(this) AliITSTransientDigit(source);
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
