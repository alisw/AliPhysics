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

/* $Id$ */

////////////////////////////////////////////////
//  Digits classes for all ITS detectors      //
//                                            //
//                                            //
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
void AliITSdigit::Print(ostream *os) {
    //Standard output format for this class

    *os << fCoord1 <<","<< fCoord2 <<","<< fSignal;
}
//______________________________________________________________________
void AliITSdigit::Read(istream *os) {
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


